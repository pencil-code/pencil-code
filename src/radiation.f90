! $Id: radiation.f90,v 1.3 2002-07-19 12:41:35 dobler Exp $

!  This modules deals with all aspects of radiation; if no
!  radiation are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  radiationally relevant subroutines listed in here.



module Radiation

  use Cparam

  implicit none

  real :: c_gam=100
  real :: opas=1e-8
  real :: mbar=1.67e-24  !mbar*m_unit in to not get to big numbers
  real :: k_B=3.99e-25      !k_B*m_unit in to not get to big numbers
  real :: a_SB=5.3e8
  real :: kappa_es=4e4

  ! init parameteres
  character (len=labellen) :: initrad='zero'

  ! input parameters
  namelist /radiation_init_pars/ &
       initrad,c_gam,opas,kappa_es,mbar,k_B,a_SB
  ! run parameters
  namelist /radiation_run_pars/ &
       c_gam,opas,kappa_es,mbar,k_B,a_SB

  ! other variables (needs to be consistent with reset list below)
  integer :: i_frms=0,i_fmax=0,i_Erms=0,i_Emax=0

  contains

!***********************************************************************
    subroutine register_rad()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!
!  15-jul-02/nils: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_rad called twice')
      first = .false.
!
      lradiation = .true.
      ie = nvar+1
      iff = nvar+2
      ifx = iff
      ify = iff+1
      ifz = iff+2
      nvar = nvar+4             ! added 4 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_rad:  nvar = ', nvar
        print*, 'ie,iff,ifx,ify,ifz = ', ie,iff,ifx,ify,ifz
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: radiation.f90,v 1.3 2002-07-19 12:41:35 dobler Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_rad: nvar > mvar')
      endif
!
    endsubroutine register_rad
!***********************************************************************
    subroutine init_rad(f,xx,yy,zz)
!
!  initialise radiation; called from start.f90
!  We have an init parameter (initrad) to stear radiation i.c. independently.
!
!   15-jul-2002/nils: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      select case(initrad)

      case('zero', '0') 
         f(:,:,:,ifx:ifz) = 0.
         f(:,:,:,ie     ) = 1.
      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'No such such value for initrad: ', trim(initrad)
        call stop_it("")

      endselect
!
    endsubroutine init_rad
!********************************************************************
    subroutine de_dt(f,df,rho1,divu,uu,uij,TT1,gamma)
!
!
!  13-Dec-01/nils: coded
!  15-Jul-02/nils: adapted from pencil_mpi
!
      use Sub
      use Cdata
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: gradE,uu,tmp,n_vec
      real, dimension (nx,3,3) :: uij, P_tens,f_mat,n_mat
      real, dimension (nx) :: E_rad,divu,rho1,source,Edivu,ugradE,divF
      real, dimension (nx) :: graduP,cooling,lgamma,RF,DFF,absgradE
      real, dimension (nx) :: f_sc,kappa_abs,kappa,E_gas,TT1,f2
      real :: gamma1,gamma
      integer :: i,j
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE dee_dt'
      if (headtt) then
        call identify_bcs('e ',ie )
        call identify_bcs('fx',ifx)
        call identify_bcs('fy',ify)
        call identify_bcs('fz',ifz)
      endif
!
      gamma1=gamma-1
      E_rad=f(l1:l2,m,n,iE)
      E_gas=1.5*k_B/(rho1*mbar*TT1)
      if (opas == 0) then
         kappa_abs=0
      else
         kappa_abs=opas*(1./rho1)**(9./2)*E_gas**(-7./2)
      endif
      Edivu=E_rad*divu
      call grad(f,iE,gradE)
      call dot_mn(uu,gradE,ugradE)
      call div(f,iff,divF)
      source=a_SB*TT1**(-4)
      kappa=kappa_abs+kappa_es
      lgamma=rho1/kappa
      cooling=(source-E_rad)*c_gam*kappa_abs/rho1
      call dot2_mn(gradE,absgradE)
      absgradE=sqrt(absgradE)
      RF=lgamma*absgradE/E_rad
      DFF=(9+RF**2)**(-0.5)
      f_sc=DFF+DFF**2*RF**2
      call multvs_mn(gradE,1./absgradE,n_vec)
      call multvv_mat_mn(n_vec,n_vec,n_mat)
!
!  calculate P_tens
!
      f_mat=0
      do i=1,3
         f_mat(:,i,i)=0.5*(1-f_sc)
      end do
      do i=1,3
         do j=1,3
            f_mat(:,i,j)=f_mat(:,i,j)+0.5*(3*f_sc-1)*n_mat(:,i,j)
         end do
      end do
      do i=1,3
         do j=1,3
            P_tens(:,i,j)=E_rad*f_mat(:,i,j)
         enddo
      enddo
      if (maxval(E_rad)==minval(E_rad)) then !Uniform rad. density
         do i=1,3
            do j=1,3
               P_tens(:,i,j)=0
               if (j==i) P_tens(:,i,i)=E_rad/3.
            enddo
         enddo
      endif
!
!  calculate graduP
!
      call multmm_sc_mn(P_tens,uij,graduP)
!
!  calculate the flux
!
      call multvs_mn(gradE,-DFF*lgamma*c_gam,tmp)
      f(l1:l2,m,n,ifx:ifz)=tmp
!
!  calculate dE/dt
!
      df(l1:l2,m,n,iE)=df(l1:l2,m,n,iE)-ugradE-divF-graduP+cooling
!
!  add (kappa F)/c to momentum equation
!
      df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+kappa*f(l1:l2,m,n,iFx)/c_gam
      df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+kappa*f(l1:l2,m,n,iFy)/c_gam
      df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+kappa*f(l1:l2,m,n,iFz)/c_gam
!
!  add cooling to entropy equation
!
!      print*,maxval(df(l1:l2,m,n,ient)),maxval(cooling*rho1*TT1/(cpour))
      df(l1:l2,m,n,ient)=df(l1:l2,m,n,ient)-(source-E_rad)*kappa_abs*c_gam*TT1
!
!  Calculate diagnostic values
!
      if (ldiagnos) then
        f2=f(l1:l2,m,n,ifx)**2+f(l1:l2,m,n,ify)**2+f(l1:l2,m,n,ifz)**2
        if (headtt.or.ldebug) print*,'Calculate maxima and rms values...'
        if (i_frms/=0) call sum_mn_name(f2,i_frms,lsqrt=.true.)
        if (i_fmax/=0) call max_mn_name(f2,i_fmax,lsqrt=.true.)
        if (i_erms/=0) call sum_mn_name(E_rad,i_erms)
        if (i_emax/=0) call max_mn_name(E_rad,i_emax)
      endif
!
!  Calculate UUmax for use in determinationof time step length
!
      if (UUmax>c_gam) call stop_it('Speed of light too small')
      UUmax=max(UUmax,c_gam)
!
    end subroutine de_dt
!*******************************************************************
    subroutine rprint_radiation(lreset)
!
!  reads and registers print parameters relevant for radiative part
!
!  16-jul-02/nils: adapted from rprint_hydro
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_frms=0; i_fmax=0; i_Erms=0; i_Emax=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'frms',i_frms)
        call parse_name(iname,cname(iname),cform(iname),'fmax',i_fmax)
        call parse_name(iname,cname(iname),cform(iname),'Erms',i_Erms)
        call parse_name(iname,cname(iname),cform(iname),'Emax',i_Emax)
      enddo
!
!  write column where which radiative variable is stored
!
      write(3,*) 'i_frms=',i_frms
      write(3,*) 'i_fmax=',i_fmax
      write(3,*) 'i_Erms=',i_Erms
      write(3,*) 'i_Emax=',i_Emax
      write(3,*) 'nname=',nname
      write(3,*) 'ie=',ie
      write(3,*) 'ifx=',ifx
      write(3,*) 'ify=',ify
      write(3,*) 'ifz=',ifz
!
    endsubroutine rprint_radiation
!***********************************************************************
end module Radiation
