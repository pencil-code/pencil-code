! $Id: cosmicray_nolog.f90,v 1.13 2004-04-30 09:30:49 ajohan Exp $

!  This modules solves the cosmic ray energy density equation.
!  It follows the description of Hanasz & Lesch (2002,2003) as used in their
!  ZEUS 3D implementation.
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************

module CosmicRay

  use Cparam
  use Cdata

  implicit none

  character (len=labellen) :: initecr='zero', initecr2='zero'

  ! input parameters
  real :: gammacr=4./3.,gammacr1
  real :: amplecr=.1,widthecr=.5,ecr_min=0.,ecr_const=0.
  real :: x_pos_cr=.0,y_pos_cr=.0,z_pos_cr=.0
  real :: x_pos_cr2=.0,y_pos_cr2=.0,z_pos_cr2=.0
  real :: amplecr2=0.,kx_ecr=1.,ky_ecr=1.,kz_ecr=1.,radius_ecr=0.,epsilon_ecr=0.

  logical :: lnegl = .false.
  logical :: lvariable_tensor_diff = .false.
  
  namelist /cosmicray_init_pars/ &
       initecr,initecr2,amplecr,amplecr2,kx_ecr,ky_ecr,kz_ecr, &
       radius_ecr,epsilon_ecr,widthecr,ecr_const, &
       gammacr, lnegl, lvariable_tensor_diff,x_pos_cr,y_pos_cr,z_pos_cr, &
       x_pos_cr2, y_pos_cr2, z_pos_cr2

  ! run parameters
  real :: cosmicray_diff=0., Kperp=0., Kpara=0., ampl_Qcr=0.
  real :: limiter_cr=1.
  logical :: simplified_cosmicray_tensor=.false.
  logical :: luse_diff_constants = .false.

  namelist /cosmicray_run_pars/ &
       cosmicray_diff,Kperp,Kpara, &
       gammacr,simplified_cosmicray_tensor,lnegl,lvariable_tensor_diff, &
       luse_diff_constants,ampl_Qcr, &
       limiter_cr

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ecrm=0,i_ecrmax=0
  integer :: i_ecrdivum=0
  integer :: i_kmax=0

  contains

!***********************************************************************
    subroutine register_cosmicray()
!
!  Initialise variables which should know that we solve for active 
!  scalar: iecr - the cosmic ray energy density; increase nvar accordingly
!
!  09-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_cosmicray called twice')
      first = .false.
!
      lcosmicray = .true.
      iecr = nvar+1            ! index to access icr
      nvar = nvar+1            ! added 1 variable
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_cosmicray: nvar = ', nvar
        print*, 'register_cosmicray: iecr = ', iecr
      endif
!
!  Put variable name in array
!
      varname(iecr) = 'ecr'
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: cosmicray_nolog.f90,v 1.13 2004-04-30 09:30:49 ajohan Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_cosmicray: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',ecr $'
          if (nvar == mvar) write(4,*) ',ecr'
        else
          write(4,*) ',ecr $'
        endif
        write(15,*) 'ecr = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_cosmicray
!***********************************************************************
    subroutine initialize_cosmicray(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  09-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!  initialize gammacr1
!
      gammacr1=gammacr-1.
      if(lroot) print*,'gammacr1=',gammacr1
!
      if(ip==0) print*,'f=',f
    endsubroutine initialize_cosmicray
!***********************************************************************
    subroutine init_ecr(f,xx,yy,zz)
!
!  initialise cosmic ray energy density field; called from start.f90
!
!   09-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Density
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz,prof
!
      select case(initecr)
        case('zero'); f(:,:,:,iecr)=0.
        case('const_ecr'); f(:,:,:,iecr)=ecr_const
        case('blob'); call blob(amplecr,f,iecr,radius_ecr,x_pos_cr,y_pos_cr,0.)
        case('gaussian-x'); call gaussian(amplecr,f,iecr,kx=kx_ecr)
        case('gaussian-y'); call gaussian(amplecr,f,iecr,ky=ky_ecr)
        case('gaussian-z'); call gaussian(amplecr,f,iecr,kz=kz_ecr)
        case('parabola-x'); call parabola(amplecr,f,iecr,kx=kx_ecr)
        case('parabola-y'); call parabola(amplecr,f,iecr,ky=ky_ecr)
        case('parabola-z'); call parabola(amplecr,f,iecr,kz=kz_ecr)
        case('gaussian-noise'); call gaunoise(amplecr,f,iecr,iecr)
        case('wave-x'); call wave(amplecr,f,iecr,kx=kx_ecr)
        case('wave-y'); call wave(amplecr,f,iecr,ky=ky_ecr)
        case('wave-z'); call wave(amplecr,f,iecr,kz=kz_ecr)
        case('propto-ux'); call wave_uu(amplecr,f,iecr,kx=kx_ecr)
        case('propto-uy'); call wave_uu(amplecr,f,iecr,ky=ky_ecr)
        case('propto-uz'); call wave_uu(amplecr,f,iecr,kz=kz_ecr)
        case('tang-discont-z')
           print*,'init_ecr: widthecr=',widthecr
        prof=.5*(1.+tanh(zz/widthecr))
        f(:,:,:,iecr)=-1.+2.*prof
        case('hor-tube'); call htube2(amplecr,f,iecr,iecr,xx,yy,zz,radius_ecr,epsilon_ecr)
        case default; call stop_it('init_ecr: bad initecr='//trim(initecr))
      endselect
!
!  superimpose something else
!
      select case(initecr2)
        case('wave-x'); call wave(amplecr2,f,iecr,ky=5.)
        case('const_ecr'); f(:,:,:,iecr)=f(:,:,:,iecr)+ecr_const
        case('blob2'); call blob(amplecr,f,iecr,radius_ecr,x_pos_cr2,y_pos_cr2,0.)
      endselect
!
      if(ip==0) print*,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_ecr
!***********************************************************************
    subroutine decr_dt(f,df,uu,rho1,divu,bij,bb)
!
!  cosmic ray evolution
!  calculate decr/dt + div(u.ecr - flux) = -pcr*divu = -(gammacr-1)*ecr*divu
!  solve as decr/dt + u.grad(ecr) = -gammacr*ecr*divu + div(flux)
!  add du = ... - (1/rho)*grad(pcr) to momentum equation
!
!   09-oct-03/tony: coded
!
      use Sub
!
      real, intent(in), dimension (mx,my,mz,mvar+maux) :: f
      real, intent(inout), dimension (mx,my,mz,mvar) :: df
      real, intent(in), dimension (nx,3,3) :: bij
      real, intent(in), dimension (nx,3) :: uu,bb
      real, intent(in), dimension (nx) :: divu,rho1
!
      real, dimension (nx,3) :: gecr
      real, dimension (nx) :: ecr,del2ecr,ugecr,vKperp,vKpara
      integer :: j
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE decr_dt'
      if (headtt) call identify_bcs('ecr',iecr)
!
!  calculate advection term
!
      call grad(f,iecr,gecr)
      call dot_mn(uu,gecr,ugecr)
!
!  Evolution equation of cosmic ray energy density
!
      ecr=f(l1:l2,m,n,iecr)
      df(l1:l2,m,n,iecr)=df(l1:l2,m,n,iecr)-ugecr-gammacr*ecr*divu
!
!  effect on the momentum equation, (1/rho)*grad(pcr)
!  cosmic ray pressure is: pcr=(gammacr-1)*ecr
!
      if(.not.lnegl)then
        do j=0,2
          df(l1:l2,m,n,iux+j)=df(l1:l2,m,n,iux+j)-gammacr1*rho1*gecr(:,1+j)
        enddo
      endif
!
!  source term added at every time step; constant for now.
!
      if (ampl_Qcr/=0.) df(l1:l2,m,n,iecr)=df(l1:l2,m,n,iecr)+ampl_Qcr
!
!  tensor diffusion, or, alternatively scalar diffusion or no diffusion
!
      if (Kperp/=0. .or. Kpara/=0. .or. lvariable_tensor_diff) then
        if(headtt) print*,'decr_dt: Kperp,Kpara=',Kperp,Kpara
        call tensor_diffusion(f,df,gecr,bij,bb,vKperp,vKpara)
      elseif (cosmicray_diff/=0.) then
        if(headtt) print*,'decr_dt: cosmicray_diff=',cosmicray_diff
        call del2(f,iecr,del2ecr)
        df(l1:l2,m,n,iecr)=df(l1:l2,m,n,iecr)+cosmicray_diff*del2ecr
      else
        if(headtt) print*,'decr_dt: no diffusion'
      endif
!
!  For the timestep calculation, need maximum diffusion
!
      if (lfirst.and.ldt) then
        if(lvariable_tensor_diff)then
           call max_for_dt(cosmicray_diff,maxval(vKperp),maxval(vKpara),maxdiffus)
        else
           call max_for_dt(cosmicray_diff,Kperp,Kpara,maxdiffus)
        endif
      endif
!
!  diagnostics
!
!  output for double and triple correlators (assume z-gradient of cc)
!  <u_k u_j d_j c> = <u_k c uu.gradecr>
!
      if (ldiagnos) then
        ecr=f(l1:l2,m,n,iecr)
        if (i_ecrdivum/=0) call sum_mn_name(ecr*divu,i_ecrdivum)
        if (i_ecrm/=0) call sum_mn_name(ecr,i_ecrm)
        if (i_ecrmax/=0) call max_mn_name(ecr,i_ecrmax)
        if (i_kmax/=0) call max_mn_name(vKperp,i_kmax)
      endif
!
    endsubroutine decr_dt
!***********************************************************************
    subroutine rprint_cosmicray(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   6-jul-02/axel: coded
!
      use Sub
!
      integer :: iname,inamez
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_ecrm=0; i_ecrdivum=0; i_ecrmax=0; i_kmax=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ecrm',i_ecrm)
        call parse_name(iname,cname(iname),cform(iname),'ecrdivum',i_ecrdivum)
        call parse_name(iname,cname(iname),cform(iname),'ecrmax',i_ecrmax)
        call parse_name(iname,cname(iname),cform(iname),'kmax',i_kmax)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
!        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ecrmz',i_ecrmz)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_ecrm=',i_ecrm
        write(3,*) 'i_ecrdivum=',i_ecrdivum
        write(3,*) 'i_ecrmax=',i_ecrmax
        write(3,*) 'i_kmax=',i_kmax
        write(3,*) 'iecr=',iecr
      endif
!
    endsubroutine rprint_cosmicray
!***********************************************************************
    subroutine tensor_diffusion(f,df,gecr,bij,bb,vKperp,vKpara)
!
!  calculates tensor diffusion with variable tensor (or constant tensor)
!  
!  vKperp*del2ecr + d_i(vKperp)d_i(gecr) + (vKpara-vKperp) d_i ( n_i n_j d_j ecr)
!      + n_i n_j d_i(ecr)d_j(vKpara-vKperp)
!   
!  = vKperp*del2ecr + gKpara.gecr + (vKpara-vKperp) (H.G + ni*nj*Gij) 
!      + ni*nj*Gi*(vKpara_j - vKperp_j),
!  where H_i = (nj bij - 2 ni nj nk bk,j)/|b| and vKperp, vKpara are variable
!  diffusion coefficients
!
!  10-oct-03/axel: adapted from pscalar
!  30-nov-03/snod: adapted from tensor_diff without variable diffusion
!  20-mar-04/axel: implemented limiter for CR advection speed such that |H|<1
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: ecr_ij,bij
      real, dimension (nx,3) :: gecr,bb,bunit,hhh,gvKperp,gvKpara
      real, dimension (nx) :: tmp,b2,b1,del2ecr,tmpj,vKperp,vKpara,tmpi
      real, dimension (nx) :: hhh2,quenchfactor
!
!  use global Kperp, Kpara ?
! 
!      real :: Kperp,Kpara
!
      integer :: i,j,k
!
!
      if (Kpara==(0.0).and.Kperp==(0.0).and.luse_diff_constants) then
          print *,"cosmicray: no diffusion"
          stop
      endif

!
!  calculate unit vector of bb
!
      call dot2_mn(bb,b2)
      b1=1./amax1(tiny(b2),sqrt(b2))
      call multsv_mn(b1,bb,bunit)
!
!  calculate first H_i (unless we use simplified_cosmicray_tensor)
!  where H_i = (nj bij - 2 ni nj nk bk,j)/|b| and vKperp, vKpara are variable
!
      if(simplified_cosmicray_tensor) then
        tmp=0.
      else
        do i=1,3
          hhh(:,i)=0.
          do j=1,3
            tmpj(:)=0.
            do k=1,3
              tmpj(:)=tmpj(:)-2.*bunit(:,k)*bij(:,k,j)
            enddo
            hhh(:,i)=hhh(:,i)+bunit(:,j)*(bij(:,i,j)+bunit(:,i)*tmpj(:))
          enddo
        enddo
        call multsv_mn(b1,hhh,hhh)
!
!  limit the length of H such that dxmin*H < 1, so we also multiply
!  by 1/sqrt(1.+dxmin^2*H^2).
!  and dot H with ecr gradient
!
        call dot2_mn(hhh,hhh2)
        quenchfactor=1./sqrt(1.+(limiter_cr*dxmin)**2*hhh2)
        call multsv_mn(quenchfactor,hhh,hhh)
        call dot_mn(hhh,gecr,tmp)
      endif
!
!  calculate Hessian matrix of ecr, dot with bi*bj, and add into tmp
!
      call g2ij(f,iecr,ecr_ij)
!
      del2ecr=0.
      do j=1,3 
        del2ecr=del2ecr+ecr_ij(:,j,j)
        do i=1,3
          tmp(:)=tmp(:)+bunit(:,i)*bunit(:,j)*ecr_ij(:,i,j)
        enddo
      enddo
!
!  if variable tensor, add extra terms and add result into decr/dt 
!
      if(lvariable_tensor_diff)then
!
!  set vKpara, vKperp
!
!  if(luse_diff  _coef)
!  
        vKpara(:)=Kpara
        vKperp(:)=Kperp
!
!  set gvKpara, gvKperp
!
        gvKperp(:,:)=0.0
        gvKpara(:,:)=0.0
!
!  put d_i ecr d_i vKperp into tmpj
!
        call dot_mn(gvKperp,gecr,tmpj)
!
!  add further terms into tmpj     
!
        do i=1,3
          tmpi(:)=bunit(:,i)*(gvKpara(:,i)-gvKperp(:,i))
          do j=1,3
            tmpj(:)=tmpj(:)+bunit(:,j)*gecr(:,j)*tmpi
          enddo
        enddo           
!
!  
!
        df(l1:l2,m,n,iecr)=df(l1:l2,m,n,iecr) & 
        + vKperp*del2ecr + (vKpara-vKperp)*tmp + tmpj 
      else
!
!  for constant tensor (or otherwise), just add result into 
!  the decr/dt equation without tmpj
!
        df(l1:l2,m,n,iecr)=df(l1:l2,m,n,iecr) &
        + Kperp*del2ecr + (Kpara-Kperp)*tmp
 
      endif     
!
    endsubroutine tensor_diffusion
!***********************************************************************

endmodule CosmicRay

