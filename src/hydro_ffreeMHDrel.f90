! $Id: hydro_ffreeMHDrel.f90,v 1.7 2003-10-16 12:50:25 mee Exp $

!  This module solve the momentum equation for relativistic force-free MHD
!  dS/dt = curlB x B +  curlE x E + divE E
!  where E = (BxS)/B^2

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxilliary variables added by this module
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Hydro

  use Cparam
  use Density
  use Viscosity  

  implicit none

  ! init parameters
  real :: ampluu=0., widthuu=.1, urand=0., kx_uu=1., ky_uu=1., kz_uu=1.
  real :: uu_left=0.,uu_right=0.,uu_lower=1.,uu_upper=1.
  real :: uy_left=0.,uy_right=0.
  real :: initpower=1.
  character (len=labellen) :: inituu='zero'


  namelist /hydro_init_pars/ &
       ampluu,inituu,widthuu,urand, &
       uu_left,uu_right,uu_lower,uu_upper,kx_uu,ky_uu,kz_uu, &
       uy_left,uy_right, &
       Omega,initpower

  ! run parameters
!ajwm - sij declaration moved to cdata.f90
  real :: theta=0.
  real :: tdamp=0.,dampu=0.,dampuext=0.,rdamp=1.2,wdamp=0.2
  real :: tau_damp_ruxm=0.,tau_damp_ruym=0.
! geodynamo
!       original line replaced and split in two
  namelist /hydro_run_pars/ &
       nu,ivisc, &            !ajwm - kept for backward comp. should 
       Omega,theta, &         ! remove and use viscosity_run_pars only
       tdamp,dampu,dampuext,dampuint,rdampext,rdampint,wdamp, &
       tau_damp_ruxm,tau_damp_ruym
! end geodynamo

!  namelist /hydro_run_pars/ &
!       nu,ivisc, &            !ajwm - kept for backward comp. should 
!       Omega,theta, &         ! remove and use viscosity_run_pars only
!       tdamp,dampu,dampuext,rdamp,wdamp, &
!       tau_damp_ruxm,tau_damp_ruym

  ! other variables (needs to be consistent with reset list below)
  integer :: i_u2m=0,i_um2=0,i_oum=0,i_o2m=0
  integer :: i_urms=0,i_umax=0,i_orms=0,i_omax=0
  integer :: i_ux2m=0, i_uy2m=0, i_uz2m=0
  integer :: i_ruxm=0,i_ruym=0,i_ruzm=0
  integer :: i_uxmz=0,i_uymz=0,i_uzmz=0,i_umx=0,i_umy=0,i_umz=0
  integer :: i_uxmxy=0,i_uymxy=0,i_uzmxy=0
  integer :: i_Marms=0,i_Mamax=0
  integer :: i_divu2m=0,i_epsK=0

  contains

!***********************************************************************
    subroutine register_hydro()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_hydro called twice')
      first = .false.
!
      lhydro = .true.
!
      iuu = nvar+1             ! indices to access uu
      iux = iuu
      iuy = iuu+1
      iuz = iuu+2
      nvar = nvar+3             ! added 3 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_hydro: nvar = ', nvar
        print*, 'register_hydro: iux,iuy,iuz = ', iux,iuy,iuz
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: hydro_ffreeMHDrel.f90,v 1.7 2003-10-16 12:50:25 mee Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_hydro: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (maux == 0) then
         if (nvar < mvar) write(4,*) ',uu $'
         if (nvar == mvar) write(4,*) ',uu'
      else
         write(4,*) ',uu $'
      endif
      write(5,*) 'uu = fltarr(mx,my,mz,3)*one'
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine initialize_hydro()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!
!  do nothing
!
    endsubroutine initialize_hydro
!***********************************************************************
    subroutine init_uu(f,xx,yy,zz)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Hydro module, if there was one.
!
!  07-nov-01/wolf: coded
!  24-nov-02/tony: renamed for consistance (i.e. init_[variable name])
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
      use Global
      use Gravity
      use Initcond
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: r,p,tmp,xx,yy,zz,prof
      real :: kabs,crit
      integer :: i
!
!  inituu corresponds to different initializations of uu (called from start).
!
      select case(inituu)

      case('zero', '0'); if(lroot) print*,'init_uu: zero velocity'
      case('gaussian-noise'); call gaunoise(ampluu,f,iux,iuz)
      case('gaussian-noise-x'); call gaunoise(ampluu,f,iux,iux)
      case('xjump'); call jump(f,iux,uu_left,uu_right,widthuu,'x')
                     call jump(f,iuy,uy_left,uy_right,widthuu,'x')
      case('Beltrami-x'); call beltrami(ampluu,f,iuu,KX=1.)
      case('Beltrami-y'); call beltrami(ampluu,f,iuu,KY=1.)
      case('Beltrami-z'); call beltrami(ampluu,f,iuu,KZ=1.)
      case('trilinear-x'); call trilinear(ampluu,f,iux,xx,yy,zz)
      case('trilinear-y'); call trilinear(ampluu,f,iuy,xx,yy,zz)
      case('trilinear-z'); call trilinear(ampluu,f,iuz,xx,yy,zz)
      case('cos-cos-sin-uz'); call cos_cos_sin(ampluu,f,iuz,xx,yy,zz)
      case('tor_pert'); call tor_pert(ampluu,f,iux,xx,yy,zz)
      case('diffrot'); call diffrot(ampluu,f,iuy,xx,yy,zz)

      case('sound-wave', '11')
        !
        !  sound wave (should be consistent with density module)
        !
        if (lroot) print*,'init_uu: x-wave in uu; ampluu=',ampluu
        f(:,:,:,iux)=ampluu*sin(kx_uu*xx)

      case('sound-wave2')
        !
        !  sound wave (should be consistent with density module)
        !
        crit=cs20-grav_const/kx_uu**2
        if (lroot) print*,'init_uu: x-wave in uu; crit,ampluu=',crit,ampluu
        if (crit>0.) then
          f(:,:,:,iux)=+ampluu*cos(kx_uu*xx)*sqrt(abs(crit))
        else
          f(:,:,:,iux)=-ampluu*sin(kx_uu*xx)*sqrt(abs(crit))
        endif

      case('shock-tube', '13')
        !
        !  shock tube test (should be consistent with density module)
        !
        if (lroot) print*,'init_uu: polytopic standing shock'
        prof=.5*(1.+tanh(xx/widthuu))
        f(:,:,:,iux)=uu_left+(uu_right-uu_left)*prof

      case('bullets')
        !
        !  blob-like velocity perturbations (bullets)
        !
        if (lroot) print*,'init_uu: velocity blobs'
        !f(:,:,:,iux)=f(:,:,:,iux)+ampluu*exp(-(xx**2+yy**2+(zz-1.)**2)/widthuu)
        f(:,:,:,iuz)=f(:,:,:,iuz)-ampluu*exp(-(xx**2+yy**2+zz**2)/widthuu)

      case('Alfven-circ-x')
        !
        !  circularly polarised Alfven wave in x direction
        !
        if (lroot) print*,'init_uu: circular Alfven wave -> x'
        f(:,:,:,iuy) = ampluu*sin(kx_uu*xx)
        f(:,:,:,iuz) = ampluu*cos(kx_uu*xx)

      case('const-ux')
        !
        !  constant x-velocity
        !
        if (lroot) print*,'init_uu: constant x-velocity'
        f(:,:,:,iux) = ampluu

      case('const-uy')
        !
        !  constant y-velocity
        !
        if (lroot) print*,'init_uu: constant y-velocity'
        f(:,:,:,iuy) = ampluu

      case('tang-discont-z')
        !
        !  tangential discontinuity: velocity is directed along x,
        !  ux=uu_lower for z<0 and ux=uu_upper for z>0. This can
        !  be set up together with 'rho-jump' in density.
        !
        if (lroot) print*,'init_uu: tangential discontinuity of uux at z=0'
        if (lroot) print*,'init_uu: uu_lower=',uu_lower,' uu_upper=',uu_upper
        if (lroot) print*,'init_uu: widthuu=',widthuu
        prof=.5*(1.+tanh(zz/widthuu))
        f(:,:,:,iux)=uu_lower+(uu_upper-uu_lower)*prof

!  Add some random noise to see the development of instability
!WD: Can't we incorporate this into the urand stuff?
        print*, 'init_uu: ampluu=',ampluu
        call random_number_wrapper(r)
        call random_number_wrapper(p)
!        tmp=sqrt(-2*alog(r))*sin(2*pi*p)*exp(-zz**2*10.)
        tmp=exp(-zz**2*10.)*cos(2.*xx+sin(4.*xx))
        f(:,:,:,iuz)=f(:,:,:,iuz)+ampluu*tmp
  
      case('Fourier-trunc')
        !
        !  truncated simple Fourier series as nontrivial initial profile
        !  for convection. The corresponding stream function is
        !    exp(-(z-z1)^2/(2w^2))*(cos(kk)+2*sin(kk)+3*cos(3kk)),
        !    with kk=k_x*x+k_y*y
        !  Not a big success (convection starts much slower than with
        !  random or 'up-down' ..
        !
        if (lroot) print*,'init_uu: truncated Fourier'
        prof = ampluu*exp(-0.5*(zz-z1)**2/widthuu**2) ! vertical Gaussian
        tmp = kx_uu*xx + ky_uu*yy               ! horizontal phase
        kabs = sqrt(kx_uu**2+ky_uu**2)
        f(:,:,:,iuz) = prof * kabs*(-sin(tmp) + 4*cos(2*tmp) - 9*sin(3*tmp))
        tmp = (zz-z1)/widthuu**2*prof*(cos(tmp) + 2*sin(2*tmp) + 3*cos(3*tmp))
        f(:,:,:,iux) = tmp*kx_uu/kabs
        f(:,:,:,iuy) = tmp*ky_uu/kabs
  
      case('up-down')
        !
        !  flow upwards in one spot, downwards in another; not soneloidal
        ! 
        if (lroot) print*,'init_uu: up-down'
        prof = ampluu*exp(-0.5*(zz-z1)**2/widthuu**2) ! vertical profile
        tmp = sqrt((xx-(x0+0.3*Lx))**2+(yy-(y0+0.3*Ly))**2) ! dist. from spot 1
        f(:,:,:,iuz) = prof*exp(-0.5*(tmp**2)/widthuu**2)
        tmp = sqrt((xx-(x0+0.5*Lx))**2+(yy-(y0+0.8*Ly))**2) ! dist. from spot 1
        f(:,:,:,iuz) = f(:,:,:,iuz) - 0.7*prof*exp(-0.5*(tmp**2)/widthuu**2)

      case('powern') 
        ! initial spectrum k^power
        call powern(ampluu,initpower,f,iux,iuz)
  
      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'init_uu: No such such value for inituu = ', trim(inituu)
        call stop_it("")

      endselect

!
!  This allows an extra random velocity perturbation on
!  top of the initialization so far.
!
      if (urand /= 0) then
        if (lroot) print*, 'init_uu: Adding random uu fluctuations'
        if (urand > 0) then
          do i=iux,iuz
            call random_number_wrapper(tmp)
            f(:,:,:,i) = f(:,:,:,i) + urand*(tmp-0.5)
          enddo
        else
          if (lroot) print*, 'init_uu:   ... multiplicative fluctuations'
          do i=iux,iuz
            call random_number_wrapper(tmp)
            f(:,:,:,i) = f(:,:,:,i) * urand*(tmp-0.5)
          enddo
        endif
      endif
!
!     if (ip==0) print*,yy,zz !(keep compiler from complaining)
    endsubroutine init_uu
!***********************************************************************
    subroutine duu_dt(f,df,uu,glnrho,divS,rho1,u2,uij,shock,gshock)
!
!  dS/dt = curlB x B +  curlE x E + divE E
!  where E = (BxS)/B^2
!
!  21-jul-03/axel: coded
!
      use Cdata
      use Sub
      use IO
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: Bij,uij
      real, dimension (nx,3) :: uu,SS,BB,CC,EE,divS,curlS,curlB,del2A,curlE
      real, dimension (nx,3) :: SgB,BgS,BdivS,CxB,curlBxB,curlExE,divEE
      real, dimension (nx,3) :: glnrho,oo,gshock
      real, dimension (nx) :: u2,B2,B21,divE,ou,o2,sij2,rho1,shock
      real, dimension (nx) :: ux,uy,uz,ux2,uy2,uz2
      real :: c2=1,B2min=1e-12
!
      intent(in) :: f,rho1
      intent(out) :: df,uu,glnrho,u2,uij,shock,gshock
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'duu_dt: SOLVE (ffreeMHDrel)'
      if (headtt) then
        call identify_bcs('ux',iux)
        call identify_bcs('uy',iuy)
        call identify_bcs('uz',iuz)
      endif
!
!  abbreviations
!
      uu=f(l1:l2,m,n,iux:iuz)
      call dot2_mn(uu,u2)
!
!  the actual calculation happens all in magnetic_ffreeMHDrel.f90
!
!  maximum squared avection speed
!
      if (headtt.or.ldebug) print*,'duu_dt: maxadvec2,u2=',maxval(maxadvec2),maxval(u2)
      if (lfirst.and.ldt) maxadvec2=amax1(maxadvec2,u2+c2)
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'duu_dt: Calculate maxima and rms values...'
        if (i_urms/=0) call sum_mn_name(u2,i_urms,lsqrt=.true.)
        if (i_umax/=0) call max_mn_name(u2,i_umax,lsqrt=.true.)
        if (i_u2m/=0) call sum_mn_name(u2,i_u2m)
        if (i_um2/=0) call max_mn_name(u2,i_um2)
        if (i_divu2m/=0) call sum_mn_name(divS**2,i_divu2m)
        if (i_ux2m/=0) then
           ux2 = uu(:,1)*uu(:,1)
           call sum_mn_name(ux2,i_ux2m)
        endif
        if (i_uy2m/=0) then
           uy2 = uu(:,2)*uu(:,2)
           call sum_mn_name(uy2,i_uy2m)
        endif
        if (i_uz2m/=0) then
           uz2 = uu(:,3)*uu(:,3)
           call sum_mn_name(uz2,i_uz2m)
        endif
!
!  mean heating term
!
        if (i_epsK/=0) then
          call multm2_mn(sij,sij2)
          call sum_mn_name(sij2,i_epsK)
        endif
!
!  this doesn't need to be as frequent (check later)
!
        if (i_uxmz/=0.or.i_uxmxy/=0) ux=uu(:,1)
        if (i_uymz/=0.or.i_uymxy/=0) uy=uu(:,2)
        if (i_uzmz/=0.or.i_uzmxy/=0) uz=uu(:,3)
        if (i_uxmz/=0) call xysum_mn_name_z(ux,i_uxmz)
        if (i_uymz/=0) call xysum_mn_name_z(uy,i_uymz)
        if (i_uzmz/=0) call xysum_mn_name_z(uz,i_uzmz)
        if (i_uxmxy/=0) call zsum_mn_name_xy(ux,i_uxmxy)
        if (i_uymxy/=0) call zsum_mn_name_xy(uy,i_uymxy)
        if (i_uzmxy/=0) call zsum_mn_name_xy(uz,i_uzmxy)
        !
        !  mean momenta
        !
        if (i_ruxm/=0) then; ux=uu(:,1); call sum_mn_name(ux,i_ruxm); endif
        if (i_ruym/=0) then; uy=uu(:,2); call sum_mn_name(uy,i_ruym); endif
        if (i_ruzm/=0) then; uz=uu(:,3); call sum_mn_name(uz,i_ruzm); endif
        !
        !  things related to vorticity
        !
        if (i_oum/=0 .or. i_o2m/=0 .or. i_omax/=0 .or. i_orms/=0) then
          oo(:,1)=Sij(:,3,2)-Sij(:,2,3)
          oo(:,2)=Sij(:,1,3)-Sij(:,3,1)
          oo(:,3)=Sij(:,2,1)-Sij(:,1,2)
          !
          if (i_oum/=0) then
            call dot_mn(oo,uu,ou)
            call sum_mn_name(ou,i_oum)
          endif
          !
          if (i_orms/=0.or.i_omax/=0.or.i_o2m/=0) then
            call dot2_mn(oo,o2)
            if(i_orms/=0) call sum_mn_name(o2,i_orms,lsqrt=.true.)
            if(i_omax/=0) call max_mn_name(o2,i_omax,lsqrt=.true.)
            if(i_o2m/=0)  call sum_mn_name(o2,i_o2m)
          endif
        endif
      endif
!
!  make sure compiler doesn't complain, so need to set them
!
      uij=0.
      glnrho=0.
      shock=0.
      gshock=0.
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine rprint_hydro(lreset)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Sub
!
      integer :: iname,inamez,ixy
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_u2m=0; i_um2=0; i_oum=0; i_o2m=0
        i_urms=0; i_umax=0; i_orms=0; i_omax=0
        i_ruxm=0; i_ruym=0; i_ruzm=0
        i_ux2m=0; i_uy2m=0; i_uz2m=0
        i_umx=0; i_umy=0; i_umz=0
        i_Marms=0; i_Mamax=0
        i_divu2m=0; i_epsK=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_hydro: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'u2m',i_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',i_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',i_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',i_oum)
        call parse_name(iname,cname(iname),cform(iname),'urms',i_urms)
        call parse_name(iname,cname(iname),cform(iname),'umax',i_umax)
        call parse_name(iname,cname(iname),cform(iname),'ux2m',i_ux2m)
        call parse_name(iname,cname(iname),cform(iname),'uy2m',i_uy2m)
        call parse_name(iname,cname(iname),cform(iname),'uz2m',i_uz2m)
        call parse_name(iname,cname(iname),cform(iname),'orms',i_orms)
        call parse_name(iname,cname(iname),cform(iname),'omax',i_omax)
        call parse_name(iname,cname(iname),cform(iname),'ruxm',i_ruxm)
        call parse_name(iname,cname(iname),cform(iname),'ruym',i_ruym)
        call parse_name(iname,cname(iname),cform(iname),'ruzm',i_ruzm)
        call parse_name(iname,cname(iname),cform(iname),'umx',i_umx)
        call parse_name(iname,cname(iname),cform(iname),'umy',i_umy)
        call parse_name(iname,cname(iname),cform(iname),'umz',i_umz)
        call parse_name(iname,cname(iname),cform(iname),'Marms',i_Marms)
        call parse_name(iname,cname(iname),cform(iname),'Mamax',i_Mamax)
        call parse_name(iname,cname(iname),cform(iname),'divu2m',i_divu2m)
        call parse_name(iname,cname(iname),cform(iname),'epsK',i_epsK)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxmz',i_uxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uymz',i_uymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzmz',i_uzmz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxmxy',i_uxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uymxy',i_uymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uzmxy',i_uzmxy)
      enddo
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_u2m=',i_u2m
      write(3,*) 'i_um2=',i_um2
      write(3,*) 'i_o2m=',i_o2m
      write(3,*) 'i_oum=',i_oum
      write(3,*) 'i_urms=',i_urms
      write(3,*) 'i_umax=',i_umax
      write(3,*) 'i_ux2m=',i_ux2m
      write(3,*) 'i_uy2m=',i_uy2m
      write(3,*) 'i_uz2m=',i_uz2m
      write(3,*) 'i_orms=',i_orms
      write(3,*) 'i_omax=',i_omax
      write(3,*) 'i_ruxm=',i_ruxm
      write(3,*) 'i_ruym=',i_ruym
      write(3,*) 'i_ruzm=',i_ruzm
      write(3,*) 'i_umx=',i_umx
      write(3,*) 'i_umy=',i_umy
      write(3,*) 'i_umz=',i_umz
      write(3,*) 'i_Marms=',i_Marms
      write(3,*) 'i_Mamax=',i_Mamax
      write(3,*) 'i_divu2m=',i_divu2m
      write(3,*) 'i_epsK=',i_epsK
      write(3,*) 'nname=',nname
      write(3,*) 'iuu=',iuu
      write(3,*) 'iux=',iux
      write(3,*) 'iuy=',iuy
      write(3,*) 'iuz=',iuz
      write(3,*) 'i_uxmz=',i_uxmz
      write(3,*) 'i_uymz=',i_uymz
      write(3,*) 'i_uzmz=',i_uzmz
      write(3,*) 'i_uxmxy=',i_uxmxy
      write(3,*) 'i_uymxy=',i_uymxy
      write(3,*) 'i_uzmxy=',i_uzmxy
!
    endsubroutine rprint_hydro
!***********************************************************************
    subroutine calc_mflow
!
!  calculate mean flow field from xy- or z-averages
!
!   8-nov-02/axel: adapted from calc_mfield
!   9-nov-02/axel: allowed mean flow to be compressible
!
      use Cdata
      use Sub
!
      logical,save :: first=.true.
      real, dimension(nx) :: uxmx,uymx,uzmx
      real, dimension(ny,nprocy) :: uxmy,uymy,uzmy
      real :: umx,umy,umz
      integer :: l,j
!
!  Magnetic energy in vertically averaged field
!  The uymxy and uzmxy must have been calculated,
!  so they are present on the root processor.
!
        if (i_umx/=0) then
          if(i_uymxy==0.or.i_uzmxy==0) then
            if(first) print*,"calc_mflow:                WARNING"
            if(first) print*, &
                    "calc_mflow: NOTE: to get umx, uymxy and uzmxy must also be set in zaver"
            if(first) print*, &
                    "calc_mflow:       We proceed, but you'll get umx=0"
            umx=0.
          else
            do l=1,nx
              uxmx(l)=sum(fnamexy(l,:,:,i_uxmxy))/(ny*nprocy)
              uymx(l)=sum(fnamexy(l,:,:,i_uymxy))/(ny*nprocy)
              uzmx(l)=sum(fnamexy(l,:,:,i_uzmxy))/(ny*nprocy)
            enddo
            umx=sqrt(sum(uxmx**2+uymx**2+uzmx**2)/nx)
          endif
          call save_name(umx,i_umx)
        endif
!
!  similarly for umy
!
        if (i_umy/=0) then
          if(i_uxmxy==0.or.i_uzmxy==0) then
            if(first) print*,"calc_mflow:                WARNING"
            if(first) print*, &
                    "calc_mflow: NOTE: to get umy, uxmxy and uzmxy must also be set in zaver"
            if(first) print*, &
                    "calc_mflow:       We proceed, but you'll get umy=0"
            umy=0.
          else
            do j=1,nprocy
            do m=1,ny
              uxmy(m,j)=sum(fnamexy(:,m,j,i_uxmxy))/nx
              uymy(m,j)=sum(fnamexy(:,m,j,i_uymxy))/nx
              uzmy(m,j)=sum(fnamexy(:,m,j,i_uzmxy))/nx
            enddo
            enddo
            umy=sqrt(sum(uxmy**2+uymy**2+uzmy**2)/(ny*nprocy))
          endif
          call save_name(umy,i_umy)
        endif
!
!  Magnetic energy in horizontally averaged field
!  The uxmz and uymz must have been calculated,
!  so they are present on the root processor.
!
        if (i_umz/=0) then
          if(i_uxmz==0.or.i_uymz==0.or.i_uzmz==0) then
            if(first) print*,"calc_mflow:               WARNING"
            if(first) print*, &
                    "calc_mflow: NOTE: to get umz, uxmz, uymz and uzmz must also be set in xyaver"
            if(first) print*, &
                    "calc_mflow:       This may be because we renamed zaver.in into xyaver.in"
            if(first) print*, &
                    "calc_mflow:       We proceed, but you'll get umz=0"
            umz=0.
          else
            umz=sqrt(sum(fnamez(:,:,i_uxmz)**2 &
                        +fnamez(:,:,i_uymz)**2 &
                        +fnamez(:,:,i_uzmz)**2)/(nz*nprocz))
          endif
          call save_name(umz,i_umz)
        endif
!
      first = .false.
    endsubroutine calc_mflow
!***********************************************************************

endmodule Hydro
