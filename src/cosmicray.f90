! $Id: cosmicray.f90,v 1.7 2003-10-12 22:13:17 mee Exp $

!  This modules solves the cosmic ray energy density advection difussion equation
!  it follows the description of Hanasz & Lesch (2002,2003) as used in their
!  ZEUS 3D implementation
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxilliary variables added by this module
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

  real, parameter :: gammacr=4./3., gammacr1=gammacr-1.

  ! input parameters
  real :: amplecr=.1, widthecr=.5, ecr_min=0., ecr_const=0.
  real :: amplecr2=0.,kx_ecr=1.,ky_ecr=1.,kz_ecr=1.,radius_ecr=0.,epsilon_ecr=0.

  namelist /cosmicray_init_pars/ &
       initecr,initecr2,amplecr,amplecr2,kx_ecr,ky_ecr,kz_ecr, &
       radius_ecr,epsilon_ecr,widthecr,ecr_const

  ! run parameters
  real :: cosmicray_diff=0., Kperp=0., Kpara=0.

  namelist /cosmicray_run_pars/ &
       cosmicray_diff,Kperp,Kpara

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ecrm=0,i_ecrmax=0

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
      nvar = nvar+1             ! added 1 variable
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_cosmicray: nvar = ', nvar
        print*, 'register_cosmicray: iecr = ', iecr
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: cosmicray.f90,v 1.7 2003-10-12 22:13:17 mee Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_cosmicray: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (maux == 0) then
         if (nvar < mvar) write(4,*) ',ecr $'
         if (nvar == mvar) write(4,*) ',ecr'
      else
         write(4,*) ',ecr $'
      endif
      write(5,*) 'ecr = fltarr(mx,my,mz)*one'
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
!  set to zero and then call the same initial condition
!  that was used in start.csh
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
        case('blob'); call blob(amplecr,f,iecr,radius_ecr,0.,0.,0.)
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
      endselect
!
      if(ip==0) print*,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_ecr
!***********************************************************************
    subroutine decr_dt(f,df,uu,divu,rho1,bij,bb)
!
!  cosmic ray evolution
!  calculate decr/dt = -uu.gecr -gammacr*ecr*divu
!  + cosmicray_diff_perp*[del2ecr + ...]
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
      real, dimension (nx) :: ecr,del2ecr,ugecr
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
!  effect on the momentum equation
!
      do j=0,2
        df(l1:l2,m,n,iux+j)=df(l1:l2,m,n,iux+j)-gammacr1*rho1*gecr(:,1+j)
      enddo
!
!  tensor diffusion, or, alternatively scalar diffusion or no diffusion
!
      if (Kperp/=0. .or. Kpara/=0.) then
        if(headtt) print*,'decr_dt: Kperp,Kpara=',Kperp,Kpara
        call tensor_diff(f,df,Kperp,Kpara,gecr,bij,bb)
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
        maxdiffus=amax1(maxdiffus,cosmicray_diff,Kperp,Kpara)
      endif
!
!  diagnostics
!
!  output for double and triple correlators (assume z-gradient of cc)
!  <u_k u_j d_j c> = <u_k c uu.gradecr>
!
      if (ldiagnos) then
        ecr=f(l1:l2,m,n,iecr)
        if (i_ecrm/=0) call sum_mn_name(ecr,i_ecrm)
        if (i_ecrmax/=0) call max_mn_name(ecr,i_ecrmax)
      endif
!
    endsubroutine decr_dt
!***********************************************************************
    subroutine rprint_cosmicray(lreset)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   6-jul-02/axel: coded
!
      use Sub
!
      integer :: iname,inamez
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_ecrm=0; i_ecrmax=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ecrm',i_ecrm)
        call parse_name(iname,cname(iname),cform(iname),'ecrmax',i_ecrmax)
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
      write(3,*) 'i_ecrm=',i_ecrm
      write(3,*) 'i_ecrmax=',i_ecrmax
      write(3,*) 'iecr=',iecr
!
    endsubroutine rprint_cosmicray
!***********************************************************************
    subroutine tensor_diff(f,df,Kperp,Kpara,gecr,bij,bb)
!
!  calculates tensor diffusion
!  Kperp*del2ecr + (Kpara-Kperp) d_i ( n_i n_j d_j ecr) 
!  = Kperp*del2ecr + (Kpara-Kperp) (H.G + ni*nj*Gij),
!  where H_i = (nj bij - 2 ni nj nk bk,j)/|b|
!
!  10-oct-03/axel: adapted from pscalar
!  10-oct-03/snod: fixed term in momentum equation
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: ecr_ij,bij
      real, dimension (nx,3) :: gecr,bb,bunit,hhh,tmpj
      real, dimension (nx) :: tmp,b2,b1,del2ecr
      real :: Kperp,Kpara
      integer :: i,j,k
!
!  calculate unit vector of bb
!
      call dot2_mn(bb,b2)
      b1=1./amax1(tiny(b2),sqrt(b2))
      call multsv_mn(b1,bb,bunit)
!
!  calculate first H_i
!
      do i=1,3
        hhh(:,i)=0.
        do j=1,3
          tmpj(:,j)=0.
          do k=1,3
            tmpj(:,j)=tmpj(:,j)-2.*bunit(:,k)*bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+bunit(:,j)*(bij(:,i,j)+bunit(:,i)*tmpj(:,j))
        enddo
      enddo
      call multsv_mn(b1,hhh,hhh)
!
!  dot H with ecr gradient
!
      call dot_mn(hhh,gecr,tmp)
!
!  calculate Hessian matrix of ecr, dot with bi*bj, and add into tmp
!
      call g2ij(f,iecr,ecr_ij)
!
      del2ecr=0.
      do j=1,3 
        del2ecr=del2ecr+ecr_ij(:,j,j)
        do i=1,3
          tmp=tmp+bunit(:,i)*bunit(:,j)*ecr_ij(:,i,j)
        enddo
      enddo
!
!  and add result to the decr/dt equation
!
      df(l1:l2,m,n,iecr)=df(l1:l2,m,n,iecr)+Kperp*del2ecr+(Kpara-Kperp)*tmp
!
    endsubroutine tensor_diff
!***********************************************************************

endmodule CosmicRay

