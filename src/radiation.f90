! $Id: radiation.f90,v 1.1 2002-07-16 08:09:33 nilshau Exp $

!  This modules deals with all aspects of radiation; if no
!  radiation are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  radiationally relevant subroutines listed in here.



module Radiation

  use Cparam

  real :: c=100

  implicit none

  ! input parameters
  namelist /radiation_init_pars/ &
       c
  ! run parameters
  namelist /radiation_run_pars/ &
       c
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
      ifx = iFF
      ify = iFF+1
      ifz = iFF+2
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
           "$Id: radiation.f90,v 1.1 2002-07-16 08:09:33 nilshau Exp $")
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
         f(:,:,:,ie     ) = 0.
      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'No such such value for initaa: ', trim(initaa)
        call stop_it("")

      endselect
!
    endsubroutine init_rad
!********************************************************************
    subroutine de_dt(f,df,rho1,divu,uu,uij,TT1)
!
!
!  13-Dec-01/nils: coded
!  15-Jul-02/nils: adapted from pencil_mpi
!
      use Sub
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: gradE,uu,tmp,n_vec
      real, dimension (nx,3,3) :: uij, P_tens,f_mat,n_mat
      real, dimension (nx) :: E_rad,divu,rho1,source,Edivu,ugradE,divF
      real, dimension (nx) :: graduP,cooling,lgamma,RF,DFF,absgradE
      real, dimension (nx) :: f_sc,kappa_abs,kappa,shearE,E_gas
      real :: gamma1
      integer :: i,j
!
      gamma1=gamma-1
      E_rad=f(l1:l2,m,n,iE)
      E_gas=1.5*k_B/(rho1*mbar*TT1)
      kappa_abs=1e-8*(1./rho1)**(9./2)*E_gas**(-7./2)
      kappa_abs=kappa_abs*kappamul
      Edivu=E_rad*divu
      call grad(f,iE,gradE)
      call dot_mn(uu,gradE,ugradE)
      call div(f,iF,divF)
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
!
!  calculate graduP
!
      call multmm_sc_mn(P_tens,uij,graduP)
!
!  calculate the flux
!
      call multvs_mn(gradE,-DFF*lgamma*c_gam,tmp)
      f(l1:l2,m,n,iF:iF+2)=tmp
!
!  calculate dE/dt
!
      df(l1:l2,m,n,iE)=df(l1:l2,m,n,iE)-ugradE-shearE-divF-graduP+cooling
!
!  add (kappa F)/c to momentum equation
!
      df(l1:l2,m,n,iu+0)=df(l1:l2,m,n,iu+0)+kappa*f(l1:l2,m,n,iF+0)/c_gam
      df(l1:l2,m,n,iu+1)=df(l1:l2,m,n,iu+1)+kappa*f(l1:l2,m,n,iF+1)/c_gam
      df(l1:l2,m,n,iu+2)=df(l1:l2,m,n,iu+2)+kappa*f(l1:l2,m,n,iF+2)/c_gam
!
!  add cooling to entropy equation
!
!      print*,maxval(df(l1:l2,m,n,iss)),maxval(cooling*rho1/(TT*cpour))
      df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)-(source-E_rad)*kappa_abs*c_gam/TT
!
    end subroutine de_dt
!*******************************************************************

end module Radiation
