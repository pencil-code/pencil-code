! $Id: hydro.f90,v 1.166 2004-05-23 10:33:33 ajohan Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
!***************************************************************

!  This module takes care of everything related to velocity

module Hydro

!  Note that Omega is already defined in cdata.

  use Cparam
!ajwm  use Cdata, only: nu,ivisc
  use Density
  use Viscosity 

  implicit none

  ! init parameters
  real :: ampluu=0., widthuu=.1, urand=0., kx_uu=1., ky_uu=1., kz_uu=1.
  real :: uu_left=0.,uu_right=0.,uu_lower=1.,uu_upper=1.
  real :: uy_left=0.,uy_right=0.
  real :: initpower=1.,cutoff=0.
  real :: nu_turb=0.,nu_turb0=0.,tau_nuturb=0.,nu_turb1=0.
  character (len=labellen) :: inituu='zero'
  real, dimension(3) :: gradH0=(/0.,0.,0./), uu_const=(/0.,0.,0./)
  complex, dimension(3) :: coefuu=(/0.,0.,0./)

  namelist /hydro_init_pars/ &
       ampluu,inituu,widthuu,urand, &
       uu_left,uu_right,uu_lower,uu_upper,kx_uu,ky_uu,kz_uu,coefuu, &
       uy_left,uy_right,uu_const, Omega,initpower,cutoff, &
       nu_turb0, tau_nuturb, nu_turb1

  ! run parameters
  real :: theta=0.
  real :: tdamp=0.,dampu=0.,wdamp=0.2
  real :: dampuint=0.0,dampuext=0.0,rdampint=0.0,rdampext=impossible
  real :: tau_damp_ruxm=0.,tau_damp_ruym=0.,tau_diffrot1=0.
  real :: ampl_diffrot=0.,Omega_int=0.
  real :: othresh=0.,othresh_per_orms=0.,orms=0.,othresh_scl=1.
  integer :: novec,novecmax=nx*ny*nz/4
  logical :: ldamp_fade=.false.,lOmega_int=.false.,lupw_uu=.false.
  logical :: lcalc_turbulence_pars
!
! geodynamo
  namelist /hydro_run_pars/ &
       nu,ivisc, &            !ajwm - kept for backward comp. should 
       Omega,theta, &         ! remove and use viscosity_run_pars only
       tdamp,dampu,dampuext,dampuint,rdampext,rdampint,wdamp, &
       tau_damp_ruxm,tau_damp_ruym,tau_diffrot1,ampl_diffrot,gradH0, &
       lOmega_int,Omega_int, ldamp_fade, lupw_uu, othresh,othresh_per_orms, &
       nu_turb0, tau_nuturb, nu_turb1, lcalc_turbulence_pars
! end geodynamo

  ! other variables (needs to be consistent with reset list below)
  integer :: i_u2m=0,i_um2=0,i_oum=0,i_o2m=0
  integer :: i_uxpt=0,i_uypt=0,i_uzpt=0
  integer :: i_dtu=0,i_dtv=0,i_urms=0,i_umax=0,i_uzrms=0,i_uzmax=0
  integer :: i_orms=0,i_omax=0
  integer :: i_ux2m=0, i_uy2m=0, i_uz2m=0
  integer :: i_ox2m=0, i_oy2m=0, i_oz2m=0
  integer :: i_uxuym=0, i_uxuzm=0, i_uyuzm=0, i_oxoym=0, i_oxozm=0, i_oyozm=0
  integer :: i_ruxm=0,i_ruym=0,i_ruzm=0,i_rumax=0
  integer :: i_uxmz=0,i_uymz=0,i_uzmz=0,i_umx=0,i_umy=0,i_umz=0
  integer :: i_uxmxy=0,i_uymxy=0,i_uzmxy=0
  integer :: i_Marms=0,i_Mamax=0
  integer :: i_divum=0,i_divu2m=0,i_epsK=0
  integer :: i_u2u13m
  integer :: i_urmphi=0,i_upmphi=0,i_uzmphi=0,i_u2mphi=0,i_oumphi=0

! Turbulence parameters
  real :: Hp,cs_ave,alphaSS,ul0,tl0,eps_diss,teta,ueta,tl01,teta1

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
!  Put variable names in array
!
      varname(iux) = 'ux'
      varname(iuy) = 'uy'
      varname(iuz) = 'uz'
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: hydro.f90,v 1.166 2004-05-23 10:33:33 ajohan Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_hydro: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',uu $'
          if (nvar == mvar) write(4,*) ',uu'
        else
          write(4,*) ',uu $'
        endif
        write(15,*) 'uu = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine initialize_hydro()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!  13-oct-03/dave: check parameters and warn (if nec.) about velocity damping
!  
!  r_int and r_ext override rdampint and rdampext if both are set
! 
   if (dampuint /= 0.0) then 
     if (r_int > epsi) then
       rdampint = r_int
     elseif (rdampint <= epsi) then
       write(*,*) 'initialize_hydro: inner radius not yet set, dampuint= ',dampuint
     endif 
   endif    

   if (dampuext /= 0.0) then      
     if (r_ext < impossible) then         
       rdampext = r_ext
     elseif (rdampext == impossible) then
       write(*,*) 'initialize_hydro: outer radius not yet set, dampuext= ',dampuext
     endif
   endif

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

      case('nothing'); if(lroot) print*,'init_uu: nothing'
      case('zero', '0'); 
                     if(lroot) print*,'init_uu: zero velocity'
                     ! Ensure really is zero, as may have used lread_oldsnap
                     f(:,:,:,iux:iuz)=0. 
      case('const_uu'); do i=1,3; f(:,:,:,iuu+i-1) = uu_const(i); enddo
      case('mode'); call modev(ampluu,coefuu,f,iuu,kx_uu,ky_uu,kz_uu,xx,yy,zz)
      case('gaussian-noise'); call gaunoise(ampluu,f,iux,iuz)
      case('gaussian-noise-x'); call gaunoise(ampluu,f,iux)
      case('gaussian-noise-y'); call gaunoise(ampluu,f,iuy)
      case('gaussian-noise-z'); call gaunoise(ampluu,f,iuz)
      case('gaussian-noise-xy'); call gaunoise(ampluu,f,iux,iuy)
      case('gaussian-noise-rprof')
        tmp=sqrt(xx**2+yy**2+zz**2)
        call gaunoise_rprof(ampluu,tmp,prof,f,iux,iuz)
      case('xjump'); call jump(f,iux,uu_left,uu_right,widthuu,'x')
                     call jump(f,iuy,uy_left,uy_right,widthuu,'x')
      case('Beltrami-x'); call beltrami(ampluu,f,iuu,kx=kx_uu)
      case('Beltrami-y'); call beltrami(ampluu,f,iuu,ky=ky_uu)
      case('Beltrami-z'); call beltrami(ampluu,f,iuu,kz=kz_uu)
      case('trilinear-x'); call trilinear(ampluu,f,iux,xx,yy,zz)
      case('trilinear-y'); call trilinear(ampluu,f,iuy,xx,yy,zz)
      case('trilinear-z'); call trilinear(ampluu,f,iuz,xx,yy,zz)
      case('cos-cos-sin-uz'); call cos_cos_sin(ampluu,f,iuz,xx,yy,zz)
      case('tor_pert'); call tor_pert(ampluu,f,iux,xx,yy,zz)
      case('diffrot'); call diffrot(ampluu,f,iuy,xx,yy,zz)
      case('olddiffrot'); call olddiffrot(ampluu,f,iuy,xx,yy,zz)
      case('sinwave-x'); call sinwave(ampluu,f,iux,kx=kx_uu)
      case('sinwave-y'); call sinwave(ampluu,f,iuy,ky=ky_uu)
      case('sinwave-z'); call sinwave(ampluu,f,iuz,kz=kz_uu)
      case('coswave-x'); call coswave(ampluu,f,iux,kx=kx_uu)
      case('coswave-y'); call coswave(ampluu,f,iuy,ky=ky_uu)
      case('coswave-z'); call coswave(ampluu,f,iuz,kz=kz_uu)
      case('soundwave-x'); call soundwave(ampluu,f,iux,kx=kx_uu)
      case('soundwave-y'); call soundwave(ampluu,f,iuy,ky=ky_uu)
      case('soundwave-z'); call soundwave(ampluu,f,iuz,kz=kz_uu)
      case('sound-wave', '11')
        !
        !  sound wave (should be consistent with density module)
        !
        if (lroot) print*,'init_uu: x-wave in uu; ampluu=',ampluu
        f(:,:,:,iux)=uu_const(1)+ampluu*sin(kx_uu*xx)

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
        call powern(ampluu,initpower,cutoff,f,iux,iuz)
  
      case('power_randomphase') 
        ! initial spectrum k^power
        call power_randomphase(ampluu,initpower,cutoff,f,iux,iuz)
  
      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'init_uu: No such such value for inituu: ', trim(inituu)
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
          if (lroot) print*, 'init_uu:  ... multiplicative fluctuations'
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
    subroutine duu_dt(f,df,uu,glnrho,divu,rho1,u2,uij,shock,gshock)
!
!  velocity evolution
!  calculate du/dt = - u.gradu - 2Omega x u + grav + Fvisc
!  pressure gradient force added in density and entropy modules.
!
!   7-jun-02/axel: incoporated from subroutine pde
!  10-jun-02/axel+mattias: added Coriolis force
!  23-jun-02/axel: glnrho and fvisc are now calculated in here
!  17-jun-03/ulf:  ux2, uy2 and uz2 added as diagnostic quantities
!
      use Cdata
      use Sub
      use IO
      use Slices
      use Special, only: special_calc_hydro
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: uij
      real, dimension (nx,3) :: uu,ugu,oo,glnrho,gshock,gui
      real, dimension (nx) :: u2,divu,o2,ou,rho1,rho,ux,uy,uz,sij2,shock,ugui
      real, dimension (nx) :: u2u13
      real :: c2,s2
      integer :: i,j
!
      intent(in) :: f,rho1
      intent(out) :: df,uu,glnrho,divu,u2,shock,gshock
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'duu_dt: SOLVE'
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
!  calculate velocity gradient matrix
!
      if (lroot .and. ip < 5) &
           print*,'duu_dt: call dot2_mn(uu,u2); m,n,iux,iuz,u2=',m,n,iux,iuz,u2
      call gij(f,iuu,uij)
      divu=uij(:,1,1)+uij(:,2,2)+uij(:,3,3)
!
!  write divu-slices for output in wvid in run.f90
!  Note: ix is the index with respect to array with ghost zones.
!
        if(lvid.and.lfirst) then
           divu_yz(m-m1+1,n-n1+1)=divu(ix-l1+1)
           if (m.eq.iy)  divu_xz(:,n-n1+1)=divu
           if (n.eq.iz)  divu_xy(:,m-m1+1)=divu
           if (n.eq.iz2) divu_xy2(:,m-m1+1)=divu
        endif
!
!  calculate rate of strain tensor
!
      if (lneed_sij) then
        do j=1,3
          do i=1,3
            sij(:,i,j)=.5*(uij(:,i,j)+uij(:,j,i))
          enddo
          sij(:,j,j)=sij(:,j,j)-.333333*divu
        enddo
      endif
!
!  advection term
!
      if (.not. lupw_uu) then
        if (ldebug) print*,'duu_dt: call multmv_mn(uij,uu,ugu)'
        call multmv_mn(uij,uu,ugu)
      else ! upwinding of velocity -- experimental and inefficent
        if (headtt) print*,'duu_dt: upwinding advection term; use at own risk!'
        !
        call grad(f,iux,gui)    ! gu=grad ux
        call u_dot_gradf(f,iux,gui,uu,ugui,UPWIND=lupw_uu)
        ugu(:,1) = ugui
        !
        call grad(f,iuy,gui)    ! gu=grad ux
        call u_dot_gradf(f,iuy,gui,uu,ugui,UPWIND=lupw_uu)
        ugu(:,2) = ugui
        !
        call grad(f,iuz,gui)    ! gu=grad ux
        call u_dot_gradf(f,iuz,gui,uu,ugui,UPWIND=lupw_uu)
        ugu(:,3) = ugui
      endif
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-ugu
!
!  Coriolis force, -2*Omega x u
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude
!
      if (Omega/=0.) then
        if (theta==0) then
          if (headtt) print*,'duu_dt: add Coriolis force; Omega=',Omega
          c2=2*Omega
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+c2*uu(:,2)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-c2*uu(:,1)
        else
          if (headtt) print*,'duu_dt: Coriolis force; Omega,theta=',Omega,theta
          c2=2*Omega*cos(theta*pi/180.)
          s2=2*Omega*sin(theta*pi/180.)
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+c2*uu(:,2)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-c2*uu(:,1)+s2*uu(:,3)
          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)           -s2*uu(:,2)
        endif
      endif
!
!  calculate grad(lnrho) here: needed for ivisc='nu-const' and continuity
!
      if(lneed_glnrho) call grad(f,ilnrho,glnrho)
!
! calculate viscous force
!
      if (lviscosity) call calc_viscous_force(f,df,glnrho,divu,rho1,shock,gshock)
!
!  maximum squared avection speed
!
      if (headtt.or.ldebug) print*,'duu_dt: maxadvec2,u2=',maxval(maxadvec2),maxval(u2)
      if (lfirst.and.ldt) call max_for_dt(u2,maxadvec2)
!
!  damp motions in some regions for some time spans if desired
!
! geodynamo
! addition of dampuint evaluation
      if (tdamp/=0.or.dampuext/=0.or.dampuint/=0) call udamping(f,df)
! end geodynamo
!
!  adding differential rotation via a frictional term
!  (should later be moved to a separate routine)
!  15-aug-03/christer: Added amplitude (ampl_diffrot) below
!
      if (tau_diffrot1/=0) then
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy) &
            -tau_diffrot1*(f(l1:l2,m,n,iuy) &
                           -ampl_diffrot*cos(x(l1:l2))*cos(z(n)))
      endif
!
!  add the possibility of removing a mean flow in the y-direction
!
      if (tau_damp_ruxm/=0.) call damp_ruxm(f,df)
      if (tau_damp_ruym/=0.) call damp_ruym(f,df)
!
!  add pressure gradient (e.g., from gas pressure in discs)
!
      do j=1,3
        if (gradH0(j)/=0.) then
          df(l1:l2,m,n,iuu+j-1)=df(l1:l2,m,n,iuu+j-1)-gradH0(j)
        endif
      enddo
!
!  interface for your personal subroutines calls
!
      if (lspecial) call special_calc_hydro(f,df,uu,glnrho,divu,rho1,u2,uij)
!
!  write oo-slices for output in wvid in run.f90
!  This must be done outside the diagnostics loop (accessed at different times).
!  Note: ix is the index with respect to array with ghost zones.
!
      if(lvid.and.lfirst) then
        oo(:,1)=uij(:,3,2)-uij(:,2,3)
        oo(:,2)=uij(:,1,3)-uij(:,3,1)
        oo(:,3)=uij(:,2,1)-uij(:,1,2)
        call dot2_mn(oo,o2)
        do j=1,3
          oo_yz(m-m1+1,n-n1+1,j)=oo(ix-l1+1,j)
          if (m.eq.iy)  oo_xz(:,n-n1+1,j)=oo(:,j)
          if (n.eq.iz)  oo_xy(:,m-m1+1,j)=oo(:,j)
          if (n.eq.iz2) oo_xy2(:,m-m1+1,j)=oo(:,j)
        enddo
        o2_yz(m-m1+1,n-n1+1)=o2(ix-l1+1)
        if (m.eq.iy)  o2_xz(:,n-n1+1)=o2
        if (n.eq.iz)  o2_xy(:,m-m1+1)=o2
        if (n.eq.iz2) o2_xy2(:,m-m1+1)=o2
        if(othresh_per_orms/=0) call calc_othresh
        call vecout(41,trim(directory)//'/ovec',oo,othresh,novec)
      endif
!
!  Calculate maxima and rms values for diagnostic purposes
!  (The corresponding things for magnetic fields etc happen inside magnetic etc)
!  The length of the timestep is not known here (--> moved to prints.f90)
!
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'duu_dt: Calculate maxima and rms values...'
        if (i_dtu/=0)    call max_mn_name(sqrt(u2)/dxmin/cdt,i_dtu,l_dt=.true.)
        if (i_urms/=0)   call sum_mn_name(u2,i_urms,lsqrt=.true.)
        if (i_umax/=0)   call max_mn_name(u2,i_umax,lsqrt=.true.)
        if (i_uzrms/=0)  call sum_mn_name(uu(:,3)**2,i_uzrms,lsqrt=.true.)
        if (i_uzmax/=0)  call max_mn_name(uu(:,3)**2,i_uzmax,lsqrt=.true.)
        if (i_rumax/=0)  call max_mn_name(u2/rho1**2,i_rumax,lsqrt=.true.)
        if (i_u2m/=0)    call sum_mn_name(u2,i_u2m)
        if (i_um2/=0)    call max_mn_name(u2,i_um2)
        if (i_divum/=0)  call sum_mn_name(divu,i_divum)
        if (i_divu2m/=0) call sum_mn_name(divu**2,i_divu2m)
        if (i_ux2m/=0)   call sum_mn_name(uu(:,1)**2,i_ux2m)
        if (i_uy2m/=0)   call sum_mn_name(uu(:,2)**2,i_uy2m)
        if (i_uz2m/=0)   call sum_mn_name(uu(:,3)**2,i_uz2m)
        if (i_uxuym/=0)  call sum_mn_name(uu(:,1)*uu(:,2),i_uxuym)
        if (i_uxuzm/=0)  call sum_mn_name(uu(:,1)*uu(:,3),i_uxuzm)
        if (i_uyuzm/=0)  call sum_mn_name(uu(:,2)*uu(:,3),i_uyuzm)
        !
        !  kinetic field components at one point (=pt)
        !
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (i_uxpt/=0) call save_name(uu(lpoint-nghost,1),i_uxpt)
          if (i_uypt/=0) call save_name(uu(lpoint-nghost,2),i_uypt)
          if (i_uzpt/=0) call save_name(uu(lpoint-nghost,3),i_uzpt)
        endif
!
!  mean heating term
!
        if (i_epsK/=0) then
          if (.not. lvisc_hyper) then
            rho=exp(f(l1:l2,m,n,ilnrho))
            call multm2_mn(sij,sij2)
            call sum_mn_name(2*nu*rho*sij2,i_epsK)
          else
            ! In this case the calculation is done in visc_hyper.f90
            itype_name(i_epsK)=ilabel_sum
          endif
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
        if (i_ruxm+i_ruym+i_ruzm/=0) rho=exp(f(l1:l2,m,n,ilnrho))
        if (i_ruxm/=0) then; ux=uu(:,1); call sum_mn_name(rho*ux,i_ruxm); endif
        if (i_ruym/=0) then; uy=uu(:,2); call sum_mn_name(rho*uy,i_ruym); endif
        if (i_ruzm/=0) then; uz=uu(:,3); call sum_mn_name(rho*uz,i_ruzm); endif
        !
        !  things related to vorticity
        !
        if (i_oum/=0 .or. i_o2m/=0 .or. i_omax/=0 .or. i_orms/=0 .or. i_ox2m/=0 .or. i_oy2m/=0 .or. i_oz2m/=0 .or. lvid) then
          oo(:,1)=uij(:,3,2)-uij(:,2,3)
          oo(:,2)=uij(:,1,3)-uij(:,3,1)
          oo(:,3)=uij(:,2,1)-uij(:,1,2)
          !
          if (i_oum/=0.or.i_oumphi/=0) then
            call dot_mn(oo,uu,ou)
            if (i_oum/=0) call sum_mn_name(ou,i_oum)
          endif
          !
          if (i_orms/=0.or.i_omax/=0.or.i_o2m/=0) then
            call dot2_mn(oo,o2)
            if(i_orms/=0) call sum_mn_name(o2,i_orms,lsqrt=.true.)
            if(i_omax/=0) call max_mn_name(o2,i_omax,lsqrt=.true.)
            if(i_o2m/=0)  call sum_mn_name(o2,i_o2m)
          endif
          !
          if (i_ox2m/=0) call sum_mn_name(oo(:,1)**2,i_ox2m)
          if (i_oy2m/=0) call sum_mn_name(oo(:,2)**2,i_oy2m)
          if (i_oz2m/=0) call sum_mn_name(oo(:,3)**2,i_oz2m)
          if (i_oxoym/=0) call sum_mn_name(oo(:,1)*oo(:,2),i_oxoym)
          if (i_oxozm/=0) call sum_mn_name(oo(:,1)*oo(:,3),i_oxozm)
          if (i_oyozm/=0) call sum_mn_name(oo(:,2)*oo(:,3),i_oyozm)
        endif
        !
        !  < u2 u1,3 >
        !
        if (i_u2u13m/=0) then
          u2u13=uu(:,2)*uij(:,1,3)
          call sum_mn_name(u2u13,i_u2u13m)
        endif
        !
      endif
!
!  phi-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        ux=uu(:,1)
        uy=uu(:,2)
        uz=uu(:,3)
        if (i_urmphi/=0) call phisum_mn_name_rz(ux*pomx+uy*pomy,i_urmphi)
        if (i_upmphi/=0) call phisum_mn_name_rz(ux*phix+uy*phiy,i_upmphi)
        if (i_uzmphi/=0) call phisum_mn_name_rz(uz,i_uzmphi)
        if (i_u2mphi/=0) call phisum_mn_name_rz(u2,i_u2mphi)
        if (i_oumphi/=0) then
          oo(:,1)=uij(:,3,2)-uij(:,2,3)
          oo(:,2)=uij(:,1,3)-uij(:,3,1)
          oo(:,3)=uij(:,2,1)-uij(:,1,2)
          call dot_mn(oo,uu,ou)
          call phisum_mn_name_rz(ou,i_oumphi)
        endif
      endif
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine calc_othresh()
!
!  calculate othresh from orms, give warnings if there are problems
!
!  24-nov-03/axel: adapted from calc_bthresh
!
      use Cdata
!
!  give warning if orms is not set in prints.in
!
      if(i_orms==0) then
        if(lroot.and.lfirstpoint) then
          print*,'calc_othresh: need to set orms in print.in to get othresh'
        endif
      endif
!
!  if nvec exceeds novecmax (=1/4) of points per processor, then begin to
!  increase scaling factor on othresh. These settings will stay in place
!  until the next restart
!
      if(novec>novecmax.and.lfirstpoint) then
        print*,'calc_othresh: processor ',iproc,': othresh_scl,novec,novecmax=', &
                                                   othresh_scl,novec,novecmax
        othresh_scl=othresh_scl*1.2
      endif
!
!  calculate othresh as a certain fraction of orms
!
      othresh=othresh_scl*othresh_per_orms*orms
!
    endsubroutine calc_othresh
!***********************************************************************
    subroutine damp_ruxm(f,df)
!
!  Damps mean x momentum, ruxm, to zero.
!  This can be useful in situations where a mean flow is generated.
!  This tends to be the case when there is linear shear but no rotation
!  and the turbulence is forced. A constant drift velocity in the
!  x-direction is most dangerous, because it leads to a linear increase
!  of <uy> due to the shear term. If there is rotation, the epicyclic
!  motion brings one always back to no mean flow on the average.
!
!  20-aug-02/axel: adapted from damp_ruym
!   7-sep-02/axel: corrected mpireduce_sum (was mpireduce_max)
!   1-oct-02/axel: turned into correction of momentum rather than velocity
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(nx) :: rho,ux
      real, dimension(1) :: fsum_tmp,fsum
      real :: tau_damp_ruxm1
      real, save :: ruxm=0.,rux_sum=0.
!
!  at the beginning of each timestep we calculate ruxm
!  using the sum of rho*ux over all meshpoints, rux_sum,
!  that was calculated at the end of the previous step.
!  This result is only known on the root processor and
!  needs to be broadcasted.
!
      if(lfirstpoint) then
        fsum_tmp(1)=rux_sum
        call mpireduce_sum(fsum_tmp,fsum,1)
        if(lroot) ruxm=fsum(1)/(nw*ncpus)
        call mpibcast_real(ruxm,1)
      endif
!
      ux=f(l1:l2,m,n,iux)
      rho=exp(f(l1:l2,m,n,ilnrho))
      call sum_mn(rho*ux,rux_sum)
      tau_damp_ruxm1=1./tau_damp_ruxm
      df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_damp_ruxm1*ruxm/rho
!
    endsubroutine damp_ruxm
!***********************************************************************
    subroutine damp_ruym(f,df)
!
!  Damps mean y momentum, ruym, to zero.
!  This can be useful in situations where a mean flow is generated.
!
!  18-aug-02/axel: coded
!   1-oct-02/axel: turned into correction of momentum rather than velocity
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(nx) :: rho,uy
      real, dimension(1) :: fsum_tmp,fsum
      real :: tau_damp_ruym1
      real, save :: ruym=0.,ruy_sum=0.
!
!  at the beginning of each timestep we calculate ruym
!  using the sum of rho*uy over all meshpoints, ruy_sum,
!  that was calculated at the end of the previous step.
!  This result is only known on the root processor and
!  needs to be broadcasted.
!
      if(lfirstpoint) then
        fsum_tmp(1)=ruy_sum
        call mpireduce_sum(fsum_tmp,fsum,1)
        if(lroot) ruym=fsum(1)/(nw*ncpus)
        call mpibcast_real(ruym,1)
      endif
!
      uy=f(l1:l2,m,n,iuy)
      rho=exp(f(l1:l2,m,n,ilnrho))
      call sum_mn(rho*uy,ruy_sum)
      tau_damp_ruym1=1./tau_damp_ruym
      df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_damp_ruym1*ruym/rho
!
    endsubroutine damp_ruym
!***********************************************************************
    subroutine udamping(f,df)
!
!  damping terms (artificial, but sometimes useful):
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(nx) :: pdamp
      real :: zbot,ztop,t_infl,t_span,tau,pfade
      integer :: i
!  
!  warn about the damping term
!
        if (headtt .and. (dampu /= 0.) .and. (t < tdamp)) then
          if (ldamp_fade) then
            print*, 'udamping: Damping velocities constantly until time ', tdamp
          else
            print*, 'udamping: Damping velocities smoothly until time ', tdamp
          end if
        endif
!
!  define bottom and top height
!
      zbot=xyz0(3)
      ztop=xyz0(3)+Lxyz(3)
!
!  1. damp motion during time interval 0<t<tdamp.
!  Damping coefficient is dampu (if >0) or |dampu|/dt (if dampu <0).
!  With ldamp_fade=T, damping coefficient is smoothly fading out
!
        if ((dampu .ne. 0.) .and. (t < tdamp)) then
          if (ldamp_fade) then  ! smoothly fade
            !
            ! smoothly fade out damping according to the following
            ! function of time:
            !
            !    ^
            !    |
            !  1 +**************
            !    |              ****
            !    |                  **
            !    |                    *
            !    |                     **
            !    |                       ****
            !  0 +-------------+-------------**********---> t
            !    |             |             |
            !    0          Tdamp/2        Tdamp
            !
            ! i.e. for 0<t<Tdamp/2, full damping is applied, while for
            ! Tdamp/2<t<Tdamp, damping goes smoothly (with continuous
            ! derivatives) to zero.
            !
            t_infl = 0.75*tdamp ! position of inflection point
            t_span = 0.5*tdamp   ! width of transition (1->0) region
            tau = (t-t_infl)/t_span ! normalized t, tr. region is [-0.5,0.5]
            if (tau <= -0.5) then
              pfade = 1.
            elseif (tau <= 0.5) then
              pfade = 0.5*(1-tau*(3-4*tau**2))
            else
              call stop_it("UDAMPING: Never got here.")
            endif
          else                ! don't fade, switch
            pfade = 1.
          endif
          !
          ! damp absolutely or relative to time step
          !
          if (dampu > 0) then   ! absolutely
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) &
                                    - pfade*dampu*f(l1:l2,m,n,iux:iuz)
          else                  ! relative to dt
            if (dt > 0) then    ! dt known and good
              df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) &
                                      + pfade*dampu/dt*f(l1:l2,m,n,iux:iuz)
            else
              call stop_it("UDAMP: dt <=0 -- what does this mean?")
            endif
          endif
        endif
!
!  2. damp motions for r_mn > rdampext or r_ext AND r_mn < rdampint or r_int
!
        if (lgravr) then
! geodynamo
! original block
!          pdamp = step(r_mn,rdamp,wdamp) ! damping profile
!          do i=iux,iuz
!            df(l1:l2,m,n,i) = df(l1:l2,m,n,i) - dampuext*pdamp*f(l1:l2,m,n,i)
!          enddo
!
          pdamp = step(r_mn,rdampext,wdamp) ! outer damping profile
          do i=iux,iuz
            df(l1:l2,m,n,i) = df(l1:l2,m,n,i) - dampuext*pdamp*f(l1:l2,m,n,i)
          enddo

          if (dampuint > 0.0) then
            pdamp = 1 - step(r_mn,rdampint,wdamp) ! inner damping profile
            do i=iux,iuz
              df(l1:l2,m,n,i) = df(l1:l2,m,n,i) - dampuint*pdamp*f(l1:l2,m,n,i)
            enddo
          endif
! end geodynamo 
        endif
!
!  coupling the above internal and external rotation rates to lgravr is not
!  a good idea. So, because of that spherical Couette flow has to be coded
!  separately.
!  ==> reconsider name <==
!
        if (lOmega_int) then
          pdamp = step(r_mn,rdampext,wdamp) ! outer damping profile
          do i=iux,iuz
            df(l1:l2,m,n,i)=df(l1:l2,m,n,i)-dampuext*pdamp*f(l1:l2,m,n,i)
          enddo
!
!  internal angular velocity, uref=(-y,x,0)*Omega_int
!
          if (dampuint > 0.0) then
            pdamp = 1 - step(r_mn,rdampint,wdamp) ! inner damping profile
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux) &
              -dampuint*pdamp*(f(l1:l2,m,n,iux)+y(m)*Omega_int)
            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy) &
              -dampuint*pdamp*(f(l1:l2,m,n,iuy)-x(l1:l2)*Omega_int)
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz) &
              -dampuint*pdamp*(f(l1:l2,m,n,iuz))
          endif
        endif
!
    endsubroutine udamping
!***********************************************************************
    subroutine rprint_hydro(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Sub
!
      integer :: iname,inamez,ixy,irz
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
        i_u2m=0; i_um2=0; i_oum=0; i_o2m=0
        i_uxpt=0; i_uypt=0; i_uzpt=0
        i_dtu=0; i_dtv=0; i_urms=0; i_umax=0; i_uzrms=0; i_uzmax=0
        i_orms=0; i_omax=0
        i_ruxm=0; i_ruym=0; i_ruzm=0; i_rumax=0
        i_ux2m=0; i_uy2m=0; i_uz2m=0; i_uxuym=0; i_uxuzm=0; i_uyuzm=0
        i_ox2m=0; i_oy2m=0; i_oz2m=0; i_oxoym=0; i_oxozm=0; i_oyozm=0
        i_umx=0; i_umy=0; i_umz=0
        i_Marms=0; i_Mamax=0
        i_divum=0; i_divu2m=0; i_epsK=0
        i_u2u13m=0
        i_urmphi=0; i_upmphi=0; i_uzmphi=0; i_u2mphi=0; i_oumphi=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'u2m',i_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',i_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',i_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',i_oum)
        call parse_name(iname,cname(iname),cform(iname),'dtu',i_dtu)
        call parse_name(iname,cname(iname),cform(iname),'dtv',i_dtv)
        call parse_name(iname,cname(iname),cform(iname),'urms',i_urms)
        call parse_name(iname,cname(iname),cform(iname),'umax',i_umax)
        call parse_name(iname,cname(iname),cform(iname),'uzrms',i_uzrms)
        call parse_name(iname,cname(iname),cform(iname),'uzmax',i_uzmax)
        call parse_name(iname,cname(iname),cform(iname),'ux2m',i_ux2m)
        call parse_name(iname,cname(iname),cform(iname),'uy2m',i_uy2m)
        call parse_name(iname,cname(iname),cform(iname),'uz2m',i_uz2m)
        call parse_name(iname,cname(iname),cform(iname),'uxuym',i_uxuym)
        call parse_name(iname,cname(iname),cform(iname),'uxuzm',i_uxuzm)
        call parse_name(iname,cname(iname),cform(iname),'uyuzm',i_uyuzm)
        call parse_name(iname,cname(iname),cform(iname),'ox2m',i_ox2m)
        call parse_name(iname,cname(iname),cform(iname),'oy2m',i_oy2m)
        call parse_name(iname,cname(iname),cform(iname),'oz2m',i_oz2m)
        call parse_name(iname,cname(iname),cform(iname),'oxoym',i_oxoym)
        call parse_name(iname,cname(iname),cform(iname),'oxozm',i_oxozm)
        call parse_name(iname,cname(iname),cform(iname),'oyozm',i_oyozm)
        call parse_name(iname,cname(iname),cform(iname),'orms',i_orms)
        call parse_name(iname,cname(iname),cform(iname),'omax',i_omax)
        call parse_name(iname,cname(iname),cform(iname),'ruxm',i_ruxm)
        call parse_name(iname,cname(iname),cform(iname),'ruym',i_ruym)
        call parse_name(iname,cname(iname),cform(iname),'ruzm',i_ruzm)
        call parse_name(iname,cname(iname),cform(iname),'rumax',i_rumax)
        call parse_name(iname,cname(iname),cform(iname),'umx',i_umx)
        call parse_name(iname,cname(iname),cform(iname),'umy',i_umy)
        call parse_name(iname,cname(iname),cform(iname),'umz',i_umz)
        call parse_name(iname,cname(iname),cform(iname),'Marms',i_Marms)
        call parse_name(iname,cname(iname),cform(iname),'Mamax',i_Mamax)
        call parse_name(iname,cname(iname),cform(iname),'divum',i_divum)
        call parse_name(iname,cname(iname),cform(iname),'divu2m',i_divu2m)
        call parse_name(iname,cname(iname),cform(iname),'epsK',i_epsK)
        call parse_name(iname,cname(iname),cform(iname),'u2u13m',i_u2u13m)
        call parse_name(iname,cname(iname),cform(iname),'uxpt',i_uxpt)
        call parse_name(iname,cname(iname),cform(iname),'uypt',i_uypt)
        call parse_name(iname,cname(iname),cform(iname),'uzpt',i_uzpt)
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
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'urmphi',i_urmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'upmphi',i_upmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uzmphi',i_uzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'u2mphi',i_u2mphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'oumphi',i_oumphi)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_u2m=',i_u2m
        write(3,*) 'i_um2=',i_um2
        write(3,*) 'i_o2m=',i_o2m
        write(3,*) 'i_oum=',i_oum
        write(3,*) 'i_dtu=',i_dtu
        write(3,*) 'i_dtv=',i_dtv
        write(3,*) 'i_urms=',i_urms
        write(3,*) 'i_umax=',i_umax
        write(3,*) 'i_uzrms=',i_uzrms
        write(3,*) 'i_uzmax=',i_uzmax
        write(3,*) 'i_ux2m=',i_ux2m
        write(3,*) 'i_uy2m=',i_uy2m
        write(3,*) 'i_uz2m=',i_uz2m
        write(3,*) 'i_uxuym=',i_uxuym
        write(3,*) 'i_uxuzm=',i_uxuzm
        write(3,*) 'i_uyuzm=',i_uyuzm
        write(3,*) 'i_ox2m=',i_ox2m
        write(3,*) 'i_oy2m=',i_oy2m
        write(3,*) 'i_oz2m=',i_oz2m
        write(3,*) 'i_oxoym=',i_oxoym
        write(3,*) 'i_oxozm=',i_oxozm
        write(3,*) 'i_oyozm=',i_oyozm
        write(3,*) 'i_orms=',i_orms
        write(3,*) 'i_omax=',i_omax
        write(3,*) 'i_ruxm=',i_ruxm
        write(3,*) 'i_ruym=',i_ruym
        write(3,*) 'i_ruzm=',i_ruzm
        write(3,*) 'i_rumax=',i_rumax
        write(3,*) 'i_umx=',i_umx
        write(3,*) 'i_umy=',i_umy
        write(3,*) 'i_umz=',i_umz
        write(3,*) 'i_Marms=',i_Marms
        write(3,*) 'i_Mamax=',i_Mamax
        write(3,*) 'i_divum=',i_divum
        write(3,*) 'i_divu2m=',i_divu2m
        write(3,*) 'i_epsK=',i_epsK
        write(3,*) 'i_u2u13m=',i_u2u13m
        write(3,*) 'nname=',nname
        write(3,*) 'iuu=',iuu
        write(3,*) 'iux=',iux
        write(3,*) 'iuy=',iuy
        write(3,*) 'iuz=',iuz
        write(3,*) 'i_uxpt=',i_uxpt
        write(3,*) 'i_uypt=',i_uypt
        write(3,*) 'i_uzpt=',i_uzpt
        write(3,*) 'i_uxmz=',i_uxmz
        write(3,*) 'i_uymz=',i_uymz
        write(3,*) 'i_uzmz=',i_uzmz
        write(3,*) 'i_uxmxy=',i_uxmxy
        write(3,*) 'i_uymxy=',i_uymxy
        write(3,*) 'i_uzmxy=',i_uzmxy
        write(3,*) 'i_urmphi=',i_urmphi
        write(3,*) 'i_upmphi=',i_upmphi
        write(3,*) 'i_uzmphi=',i_uzmphi
        write(3,*) 'i_u2mphi=',i_u2mphi
        write(3,*) 'i_oumphi=',i_oumphi
      endif
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
            if(first) print*, 'calc_mflow:                    WARNING'
            if(first) print*, &
                    "calc_mflow: NOTE: to get umx, uymxy and uzmxy must also be set in zaver"
            if(first) print*, &
                    "calc_mflow:      We proceed, but you'll get umx=0"
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
            if(first) print*, 'calc_mflow:                    WARNING'
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
            if(first) print*,"calc_mflow:                    WARNING"
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
    subroutine calc_turbulence_pars(f)
!
!  Calculate turbulence parameters for a disc. 
!  Currently only works in parallel when nprocy=1
!
!  18-may-04/anders: programmed
!
      use Cdata
      use Cparam
      use Mpicomm
      use Ionization, only: pressure_gradient,eoscalc,ilnrho_ss
      use Viscosity, only: nu_mol

      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(nx) :: cs2,cp1tilde
      real, dimension(1) :: cs_sum_allprocs_arr,Hp_arr
      real :: cs_sum_thisproc,pp0,pp1,pp2
      integer :: iprocHp
!
!  Calculate turbulent viscosity
!
      if (tau_nuturb == 0.) then
        nu_turb = nu_turb0
      else
        nu_turb = nu_turb0*exp(-t/tau_nuturb)
        if (nu_turb < nu_turb1) nu_turb = nu_turb1
      endif
!
!  Calculate average sound speed in disc
!
      do m=m1,m2
        do n=n1,n2
          call pressure_gradient(f,cs2,cp1tilde)
          cs_sum_thisproc = cs_sum_thisproc + sum(sqrt(cs2))
        enddo
      enddo
!
!  Get sum of cs_sum_thisproc_arr on all procs
!
      call mpireduce_sum((/ cs_sum_thisproc /),cs_sum_allprocs_arr,1)
!
!  Calculate average cs
!        
      if (lroot) cs_ave = cs_sum_allprocs_arr(1)/nwgrid
!
!  Send to all procs
!          
      call mpibcast_real(cs_ave,1)
!
!  Need mid-plane pressure for pressure scale height calculation
!
      if (iproc == nprocz/2) then
        if (nprocz == 2*(nprocz/2)) then  ! Even no. of procs in z
          call eoscalc(ilnrho_ss,f(lpoint,mpoint,n1,ilnrho),&
              f(lpoint,mpoint,n1,iss),pp=pp0)
        else                              ! Odd no. of procs in z
          call eoscalc(ilnrho_ss,f(lpoint,mpoint,npoint,ilnrho),&
              f(lpoint,mpoint,npoint,iss),pp=pp0)
        endif
      endif
      call mpibcast_real(pp0,1,nprocz/2)
!
!  Find pressure scale height and calculate turbulence properties
!
      Hp = 0.
      pp2 = 0.
      do n=n1,n2
        pp1 = pp2
        call eoscalc(ilnrho_ss,f(lpoint,mpoint,n,ilnrho), &
            f(lpoint,mpoint,n,iss),pp=pp2)
        if (pp1 > 0.367879*pp0 .and. pp2 <= 0.367879*pp0) then
!
!  Interpolate linearly between z1 and z2 (P_1+(P_2-P_1)/dz*Delta z = 1/e*P_0)
!          
          Hp = z(n-1) + dz/(pp2-pp1)*(0.367879*pp0-pp1)
          exit
        endif
      enddo
!
!  Broadcast scale height to all processors (Hp is 0 except where Hp is found)
!
      call mpireduce_sum((/ Hp /),Hp_arr,1)
      if (lroot) Hp=Hp_arr(1)
      call mpibcast_real(Hp,1)
!  Shakury-Sunyaev alpha      
      alphaSS = nu_turb/(Hp**2*Omega)
!  Speed of largest scale      
      ul0  = alphaSS*cs_ave
!  Eddy turn over time of largest scale      
      tl0  = Hp/ul0
!  Energy dissipation rate for Kolmogorov spectrum      
      eps_diss = nu_turb*(qshear*Omega)**2
!  Speed of smallest (viscous) scale      
      ueta = (nu_mol*eps_diss)**0.25
!  Eddy turn over time of smallest (viscous) scale      
      teta = (nu_mol/eps_diss)**0.5
!  Auxiliary      
      tl01 = 1/tl0
      teta1 = 1/teta

    endsubroutine calc_turbulence_pars
!***********************************************************************

endmodule Hydro
