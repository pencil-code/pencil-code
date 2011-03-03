! $Id$
!
!  This module takes care of simple types of gravity, i.e. where
!    gx=gx(x) or gy=gy(y) or gz=gz(z)
!  Here the gravity master pencils gravx_xpencil, gravy_ypencil and
!  gravz_zpencil only need to be calculated once, and then these can
!  simply be added to the equations of motion again and again.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lgrav = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED gg(3); epot
!
!***************************************************************
module Gravity
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'gravity.h'
!
  interface potential
    module procedure potential_global
    module procedure potential_penc
    module procedure potential_point
  endinterface
!
  interface acceleration
    module procedure acceleration_penc
    module procedure acceleration_penc_1D
    module procedure acceleration_point
  endinterface
!
  double precision, parameter :: g_B_cgs=6.172d20, g_D_cgs=3.086d21
  double precision :: g_B, g_D
  real, dimension(mx) :: gravx_xpencil=0.0,potx_xpencil=0.0
  real, dimension(my) :: gravy_ypencil=0.0,poty_ypencil=0.0
  real, dimension(mz) :: gravz_zpencil=0.0,potz_zpencil=0.0
  real, dimension(mx) :: xdep=0.0
  real, dimension(mz) :: zdep=0.0
  real, parameter :: g_A_cgs=4.4e-9, g_C_cgs=1.7e-9
  real :: gravx=0.0, gravy=0.0, gravz=0.0
  real :: kx_gg=1.0, ky_gg=1.0, kz_gg=1.0, gravz_const=1.0, reduced_top=1.0
  real :: xgrav=impossible, ygrav=impossible, zgrav=impossible
  real :: xinfty=0.0, yinfty=0.0, zinfty=impossible
  real :: dgravx=0.0, pot_ratio=1.0, qgshear=1.5
  real :: z1=0.0, z2=1.0, zref=impossible, sphere_rad=0.0, g_ref=0.0
  real :: nu_epicycle=1.0, nu_epicycle2=1.0
  real :: nux_epicycle=0.0, nux_epicycle2=0.0
  real :: r0_pot=0.0
  real :: g_A, g_C
  real :: cs0hs=0.0, H0hs=0.0
  integer :: n_pot=10
  character (len=labellen) :: gravx_profile='zero',gravy_profile='zero'
  character (len=labellen) :: gravz_profile='const'
!
!  Parameters used by other modules (only defined for other gravities)
!
  logical :: lnumerical_equilibrium=.false.
  logical :: lxyzdependence=.false.
  logical :: lcalc_zinfty=.false.
  logical :: lboussinesq=.false.
  real :: g0=0.0
  real :: lnrho_bot=0.0, lnrho_top=0.0, ss_bot=0.0, ss_top=0.0
  real :: kappa_x1=0.0, kappa_x2=0.0, kappa_z1=0.0, kappa_z2=0.0
  real :: grav_tilt=0.0, grav_amp=0.0
!
  namelist /grav_init_pars/ &
      gravx_profile, gravy_profile, gravz_profile, gravx, gravy, gravz, &
      xgrav, ygrav, zgrav, kx_gg, ky_gg, kz_gg, dgravx, pot_ratio, z1, z2, &
      nux_epicycle, nu_epicycle, zref, g_ref, sphere_rad, lnrho_bot, lnrho_top, &
      ss_bot, ss_top, lgravx_gas, lgravx_dust, lgravy_gas, lgravy_dust, &
      lgravz_gas, lgravz_dust, xinfty, yinfty, zinfty, lxyzdependence, &
      lcalc_zinfty, kappa_x1, kappa_x2, kappa_z1, kappa_z2, reduced_top, &
      lboussinesq, n_pot, cs0hs, H0hs, grav_tilt, grav_amp
!
  namelist /grav_run_pars/ &
      gravx_profile, gravy_profile, gravz_profile, gravx, gravy, gravz, &
      xgrav, ygrav, zgrav, kx_gg, ky_gg, kz_gg, dgravx, pot_ratio, &
      nux_epicycle, nu_epicycle, zref, g_ref, sphere_rad, &
      lgravx_gas, lgravx_dust, lgravy_gas, lgravy_dust, &
      lgravz_gas, lgravz_dust, xinfty, yinfty, zinfty, lxyzdependence, &
      lcalc_zinfty, kappa_x1, kappa_x2, kappa_z1, kappa_z2, reduced_top, &
      lboussinesq, n_pot, grav_tilt, grav_amp
!
  integer :: idiag_epot=0
  integer :: idiag_epotmx=0
  integer :: idiag_epotmy=0
  integer :: idiag_epotmz=0
  integer :: idiag_epotuzmz=0
!
  contains
!***********************************************************************
    subroutine register_gravity()
!
!  Initialise gravity variables (currently none).
!
!  12-nov-04/anders: coded
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          '$Id$')
!
    endsubroutine register_gravity
!***********************************************************************
    subroutine initialize_gravity(f,lstarting)
!
!  Calculate master pencils for gravity. These are put into gravity pencils
!  in the subroutine calc_pencils_grav.
!
!  12-nov-04/anders: coded, copied init conds from grav_x, grav_y and grav_y.
!
      use SharedVariables, only: put_shared_variable, get_shared_variable
      use Sub, only: notanumber, cubic_step
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, pointer :: cs20,mpoly,gamma
      logical :: lstarting
!
      real, dimension (mz) :: prof
      real :: ztop
      integer :: ierr
!
!  Sanity check.
!
      if (.not. lstarting) then
        if (gravx_profile == 'zero' .and. &
            gravy_profile == 'zero' .and. &
            gravz_profile == 'zero') then
          call fatal_error('initialize_gravity', &
              'You do not need gravity_simple for zero gravity...')
        endif
      endif
!
!  Possibility of specifying zref (if zinfty is not given).
!  In that case we need to know cs20 and mpoly from the EOS module.
!  We do this using shared variables.
!
      if (lcalc_zinfty) then
        if (zinfty==impossible) then
          if (zref==impossible) then
            call fatal_error('initialize_gravity','zref=impossible')
          else
            call get_shared_variable('cs20',cs20,ierr)
            if (ierr/=0) call fatal_error('initialize_gravity','getting cs20')
!
            call get_shared_variable('mpoly',mpoly,ierr)
            if (ierr/=0) call fatal_error('initialize_gravity','getting mpoly')
!
            call get_shared_variable('gamma',gamma,ierr)
            if (ierr/=0) call fatal_error('initialize_gravity','getting gamma')
!
            zinfty=zref+cs20*(mpoly+1)/(-gamma*gravz)
            if (lroot) print*,'initialize_gravity: computed zinfty=',zinfty
          endif
        endif
      else
        zinfty=0.
      endif
!
!  Different x-gravity profiles.
!
      if (gravx/=0) lgravx=.true.
!
      select case (gravx_profile)
!
      case ('zero')
        if (lroot) print*,'initialize_gravity: no x-gravity'
        gravx_xpencil=0.
        lgravx_gas=.false.
        lgravx_dust=.false.
!
      case ('const')
        if (lroot) print*,'initialize_gravity: constant x-grav=', gravx
        gravx_xpencil=gravx
        potx_xpencil=-gravx*(x-xinfty)
        call put_shared_variable('gravx', gravx, ierr)
!
      case ('const_tilt')
        gravx=grav_amp*sin(grav_tilt*pi/180.)
        if (lroot) print*,'initialize_gravity: constant gravx=', gravx
        gravx_xpencil=gravx
        potx_xpencil=-gravx*(x-xinfty)
        call put_shared_variable('gravx', gravx, ierr)
!
!  Linear gravity potential with additional z dependence.
!  Calculate zdep here, but don't multiply it onto gravx_xpencil
!  or potx_xpencil, so they will not be correct yet!
!
      case ('linear_zdep')      !  Linear gravity profile (for accretion discs)
        nux_epicycle2=nux_epicycle**2
        if (lroot) print*,'initialize_gravity: linear x-grav with z dep, nu=', &
          nux_epicycle,kappa_z1,kappa_z2
        zdep=(1.+kappa_z1*z+.5*(kappa_z1*z)**2)
        gravx_xpencil = -nux_epicycle2*x
        potx_xpencil=0.5*nux_epicycle2*(x**2-xinfty**2)
!
!  tanh profile
!  For isothermal EOS, we have 0=-cs2*dlnrho+gravx.
!  pot_ratio gives the resulting ratio in the density.
!
      case ('tanh-pot')
        if (dgravx==0.) call fatal_error("initialize_gravity","dgravx=0 not OK")
        if (lroot) print*,'initialize_gravity: tanh x-grav, gravx=',gravx
        if (lroot) print*,'initialize_gravity: xgrav,dgravx=',xgrav,dgravx
        gravx=-alog(pot_ratio)/dgravx
        gravx_xpencil=gravx*.5/cosh((x-xgrav)/dgravx)**2
        potx_xpencil=-gravx*.5*(1.+tanh((x-xgrav)/dgravx))*dgravx
!
      case ('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal x-grav, gravx=',gravx
        gravx_xpencil = -gravx*sin(kx_gg*x)
!
      case ('Baker74')
        !abag added, need to make adaptable to any width/position
        if (lroot) print*,'initialize_gravity: Baker_74=',gravx
        gravx_xpencil =(tanh((x+pi/3.)/0.1)+tanh(-(x-pi/3.)/0.1))/2.*&
            gravx*sin(2*(x-pi/2.))
        potx_xpencil=(tanh((x+pi/3.)/0.1)+tanh(-(x-pi/3.)/0.1))/2.*&
            gravx*(.5*cos(2*(x-pi/2.))-0.5)
!
      case ('kepler')
        if (lroot) print*,'initialize_gravity: kepler x-grav, gravx=',gravx
        gravx_xpencil=-gravx/x**2
        potx_xpencil=-gravx/x
        g0=gravx
        call put_shared_variable('gravx', gravx, ierr)
!
      case ('loop')
        if (lroot) print*,'initialize_gravity: 1D loop, gravx=',gravx
        gravx_xpencil=cos((x-xyz0(1)) / (xyz0(1)+Lxyz(1)) * pi)
        gravx_xpencil=gravx * gravx_xpencil
        potx_xpencil =sin((x-xyz0(1)) / (xyz0(1)+Lxyz(1)) * pi)
        potx_xpencil =gravx*(xyz0(1)+Lxyz(1))/pi * potx_xpencil
!
      case default
        if (lroot) print*, &
            'initialize_gravity: unknown gravx_profile ', gravx_profile
        call fatal_error('initialize_gravity','chosen gravx_profile not valid')
!
      endselect
!
!  Different y-gravity profiles.
!
      select case (gravy_profile)
!
      case ('zero')
        if (lroot) print*,'initialize_gravity: no y-gravity'
        gravy_ypencil=0.
        lgravy_gas=.false.
        lgravy_dust=.false.
!
      case ('const')
        if (lroot) print*,'initialize_gravity: constant y-grav=', gravy
        gravy_ypencil=gravy
        poty_ypencil=-gravy*(y-yinfty)
!
      case ('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal y-grav, gravy=', gravy
        gravy_ypencil = -gravy*sin(ky_gg*y)
!
      case ('kepler')
        if (lroot) print*,'initialize_gravity: kepler gravy, y-grav=', gravy
        gravy_ypencil = -gravy/y**2
!
      case default
        if (lroot) print*, &
            'initialize_gravity: unknown gravy_profile ', gravy_profile
        call fatal_error('initialize_gravity','chosen gravy_profile not valid')
!
      endselect
!
!  Different z-gravity profiles.
!  Set lgravz=T only when gravz_profile is not zero.
!
      if (gravz_profile/='zero') lgravz=.true.
!
      select case (gravz_profile)
!
      case ('zero')
        if (lroot) print*,'initialize_gravity: no z-gravity'
        gravz_zpencil=0.
        lgravz_gas=.false.
        lgravz_dust=.false.
!
      case ('const')
        if (lroot) print*,'initialize_gravity: constant gravz=', gravz
        gravz_zpencil=gravz
        potz_zpencil=-gravz*(z-zinfty)
!
      case ('const_tilt')
        gravz=-grav_amp*cos(grav_tilt*pi/180.)
        if (lroot) print*,'initialize_gravity: constant gravz=', gravz
        gravz_zpencil=gravz
        potz_zpencil=-gravz*(z-zinfty)
!
      case ('tanh')
        if (lroot) print*,'initialize_gravity: tanh gravz=', gravz
        gravz_zpencil=-gravz*tanh(z/zref)
        potz_zpencil=gravz*zref*alog(cosh(z/zref))
!
      case ('boussinesq')
        if (lroot) print*,'initialize_gravity: boussinesq gravz=', gravz
        gravz_zpencil=gravz
        potz_zpencil=-gravz*(z-zinfty)
        lboussinesq=.true.
!
      case ('const_zero')  !  Const. gravity acc. (but zero for z>zgrav)
        if (headtt) print*,'initialize_gravity: const_zero gravz=', gravz
        if (zgrav==impossible .and. lroot) &
            print*,'initialize_gravity: zgrav is not set!'
        do n=n1,n2
          if (z(n)<=zgrav) gravz_zpencil(n) = gravz
        enddo
!
      case ('linear')      !  Linear gravity profile (for accretion discs)
        nu_epicycle2=nu_epicycle**2
        if (lroot) print*,'initialize_gravity: linear z-grav, nu=', nu_epicycle
        gravz_zpencil=-nu_epicycle2*z
        potz_zpencil=0.5*nu_epicycle2*(z**2-zinfty**2)
!
      ! solid sphere (homogenous) gravity profile
      ! 'zref' determines the location of the sphere border relative to the lower box boundary.
      ! 'g_ref' is the gravity acceleration at zref (implies the mass of the sphere).
      ! 'sphere_rad' is the radius of the sphere (independant from box coordinates).
      ! Together, sphere_rad and zref describe, how much of the sphere lies inside the box.
      case ('solid_sphere')
        if (lroot) print *, 'initialize_gravity: solid shpere zref=', zref, ', g_ref=', g_ref, ', sphere_rad=', sphere_rad
        where (z > zref)
          gravz_zpencil = g_ref * sphere_rad**2 / (z-zref+sphere_rad)**2
        else where (z < zref-2*sphere_rad)
          gravz_zpencil = -g_ref * sphere_rad**2 / (z-zref+sphere_rad)**2
        else where
          gravz_zpencil = g_ref * (z-zref+sphere_rad) / sphere_rad
        end where
        potz_zpencil=-gravz_zpencil*(z-zinfty)
!
      case ('spherical')
        nu_epicycle2=nu_epicycle**2
        if (lroot) print*,'initialize_gravity: spherical z-grav, nu, z1=', nu_epicycle, z1
        gravz_zpencil=-nu_epicycle2*z/(1.0+(z/z1)**2)
        potz_zpencil=0.5*nu_epicycle2*z1**2*alog(1.0+(z/z1)**2)
!
!  Linear gravity potential with additional x dependence.
!  Calculate xdep here, but don't multiply it onto gravz_zpencil
!  or potz_zpencil, so they will not be correct yet!
!
      case ('linear_xdep')
        nu_epicycle2=nu_epicycle**2
        if (lroot) print*,'initialize_gravity: linear z-grav with x dep, nu=', &
          nu_epicycle,kappa_x1,kappa_x2
        xdep=(1.+kappa_x1*x+.5*(kappa_x1*x)**2)
        gravz_zpencil = -nu_epicycle2*z
        potz_zpencil=0.5*nu_epicycle2*(z**2-zinfty**2)
!
      case ('linear_smoothed')
        nu_epicycle2=nu_epicycle**2
        if (lroot) print*,'initialize_gravity: linear z-grav, '// &
                         'smoothed to zero at top/bottom, nu=', nu_epicycle
        prof = 1. + (z/zref)**(2*n_pot)
        gravz_zpencil = -nu_epicycle2*z/prof**(1./n_pot+1.)
        potz_zpencil = 0.5*nu_epicycle2*z**2/prof**(1./n_pot)
!
      case ('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal z-grav, gravz=', gravz
        gravz_zpencil = -gravz*sin(kz_gg*z)
        potz_zpencil = -gravz/kz_gg*cos(kz_gg*z)
!
      case ('kepler')
        if (lroot) print*,'initialize_gravity: kepler z-grav, gravz=', gravz
        gravz_zpencil = -gravz/z**2
!
      case ('Ferriere')
!
!  Set up physical units.
!
      if (unit_system=='cgs') then
        g_A = g_A_cgs/unit_velocity*unit_time
        g_B = g_B_cgs/unit_length
        g_C = g_C_cgs/unit_velocity*unit_time
        g_D = g_D_cgs/unit_length
      else if (unit_system=='SI') then
        call fatal_error('initialize_gravity','SI unit conversions not inplemented')
      endif
!
!  Gravity profile from K. Ferriere, ApJ 497, 759, 1998, eq (34)
!  at solar radius.  (for interstellar runs)
!
!  nb: 331.5 is conversion factor: 10^-9 cm/s^2 -> kpc/Gyr^2)  (/= 321.1 ?!?)
!AB: These numbers should be inserted in the appropriate units.
!AB: As it is now, it can never make much sense.
        gravz_zpencil = -(g_A*z/sqrt(z**2+g_B**2) + g_C*z/g_D)
!
      case ('Galactic-hs')
        if (lroot) print*,'Galactic hydrostatic equilibrium gravity profile'
        if (lroot.and.(cs0hs==0.or.H0hs==0)) &
            call fatal_error('initialize-gravity', &
            'Set cs0hs and H0hs in grav_init_pars!')
        gravz_zpencil = -z*(cs0hs/H0hs)**2/sqrt(1 + (z/H0hs)**2)
!
      case ('reduced_top')
        if (lroot) print*,'initialize_gravity: reduced, gravz=',gravz
        if (zgrav==impossible.and.lroot) print*,'zgrav is not set!'
        ztop = xyz0(3)+Lxyz(3)
        prof = cubic_step(z,(zgrav+ztop)/2,(ztop-zgrav)/2)
        gravz_zpencil = (1 - prof*(1-reduced_top))*gravz
!
      case default
        if (lroot) print*, &
            'initialize_gravity: unknown gravz_profile ', gravz_profile
        call fatal_error('initialize_gravity','chosen gravz_profile not valid')
!
      endselect
!
!  Sanity check.
!
      if (notanumber(gravx_xpencil)) then
        call fatal_error('initialize_gravity','found NaN or +/-Inf in gravx_xpencil')
      endif
      if (notanumber(gravy_ypencil)) then
        call fatal_error('initialize_gravity','found NaN or +/-Inf in gravy_ypencil')
      endif
      if (notanumber(gravz_zpencil)) then
        call fatal_error('initialize_gravity','found NaN or +/-Inf in gravz_zpencil')
      endif
!
      call put_shared_variable('nu_epicycle',nu_epicycle)
      call put_shared_variable('gravz_zpencil',gravz_zpencil)
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine init_gg(f)
!
!  Initialise gravity; called from start.f90.
!
!  12-nov-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Don't do anything
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_gg
!***********************************************************************
    subroutine pencil_criteria_gravity()
!
!  All pencils that the Gravity module depends on are specified here.
!
!  20-11-04/anders: coded
!
      lpenc_requested(i_gg)=.true.
      if (lanelastic) lpenc_requested(i_rhop)=.true.
!
      if (idiag_epot/=0 .or. idiag_epotmx/=0 .or. idiag_epotmy/=0 .or. &
          idiag_epotmz/=0) lpenc_diagnos(i_epot)=.true.
      if (idiag_epotuzmz/=0) then
        lpenc_diagnos(i_epot)=.true.
        lpenc_diagnos(i_uu)=.true.
      endif
!
    endsubroutine pencil_criteria_gravity
!***********************************************************************
    subroutine pencil_interdep_gravity(lpencil_in)
!
!  Interdependency among pencils from the Gravity module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_epot)) lpencil_in(i_rho)=.true.
!
    endsubroutine pencil_interdep_gravity
!***********************************************************************
    subroutine calc_pencils_gravity(f,p)
!
!  Calculate Gravity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  12-nov-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      if (lpencil(i_gg)) then
        p%gg(:,1) = gravx_xpencil(l1:l2)
        p%gg(:,2) = gravy_ypencil(m)
        p%gg(:,3) = gravz_zpencil(n)
      endif
!
      if (lpencil(i_epot)) p%epot=p%rho* &
          (potx_xpencil(l1:l2)+poty_ypencil(m)+potz_zpencil(n))
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_gravity
!***********************************************************************
    subroutine duu_dt_grav(f,df,p)
!
!  Add gravitational acceleration to gas and dust.
!
!  The special option lboussinesq=T is applicable when |z|/H  << 1.
!  However, in the present formulation the resulting equations,
!  du/dt = -lnrho, and dlnrho/dt=-du/dz, lead to an instability
!  with the growth rate lambda = (1+i)*sqrt(k/2).
!
!  12-nov-04/anders: coded
!   5-dec-06/petri: added Boussinesq approximation
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: k
!
      intent(in) :: f,p
      intent(out) :: df
!
!  Add gravity acceleration on gas.
!
      if (lhydro) then
        if (lboussinesq) then
          if (lentropy) then
            if (lgravx_gas) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+p%ss*p%gg(:,1)
            if (lgravy_gas) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+p%ss*p%gg(:,2)
            if (lgravz_gas) df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+p%ss*p%gg(:,3)
          else
            if (headtt) print*,'duu_dt_grav: lboussinesq w/o lentropy not ok!'
          endif
        else if (lanelastic) then
! Now works for the linear anelastic formulation only
                if (lgravx_gas) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+ p%gg(:,1)*&
                                p%rhop
                if (lgravy_gas) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+ p%gg(:,2)*&
                                p%rhop
                if (lgravz_gas) df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+p%gg(:,3)*&
                                 (-p%ss)
!                                p%rhop
        else
          if (lxyzdependence) then
            if (lgravx_gas) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+p%gg(:,1)*zdep(n)
            if (lgravy_gas) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+p%gg(:,2)
            if (lgravz_gas) df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+p%gg(:,3)*xdep(l1:l2)
          else
            if (lgravx_gas) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+p%gg(:,1)
            if (lgravy_gas) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+p%gg(:,2)
            if (lgravz_gas) df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+p%gg(:,3)
          endif
        endif
      endif
!
!  Add gravity acceleration on dust.
!
      if (ldustvelocity) then
        do k=1,ndustspec
          if (lgravx_dust) &
              df(l1:l2,m,n,iudx(k)) = df(l1:l2,m,n,iudx(k)) + p%gg(:,1)
          if (lgravy_dust) &
              df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) + p%gg(:,2)
          if (lgravz_dust) &
              df(l1:l2,m,n,iudz(k)) = df(l1:l2,m,n,iudz(k)) + p%gg(:,3)
        enddo
      endif
!
!  Gravity diagnostics.
!
      if (ldiagnos) then
        if (idiag_epot/=0) call sum_mn_name(p%epot,idiag_epot)
      endif
!
!  Gravity 1-D diagnostics.
!
      if (l1davgfirst) then
        if (idiag_epotmx/=0) call yzsum_mn_name_x(p%epot,idiag_epotmx)
        if (idiag_epotmy/=0) call xzsum_mn_name_y(p%epot,idiag_epotmy)
        if (idiag_epotmz/=0) call xysum_mn_name_z(p%epot,idiag_epotmz)
        if (idiag_epotuzmz/=0) call xysum_mn_name_z(p%epot*p%uu(:,3), &
            idiag_epotuzmz)
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential_global(pot,pot0)
!
!  Calculates gravity potential globally
!
!  13-nov-04/anders: coded
!
      real, dimension (mx,my,mz) :: pot
      real, optional :: pot0
!
      call fatal_error('potential_global','this subroutine has been '// &
          'deprecated for gravity_simple')
!
      call keep_compiler_quiet(pot)
      if (present(pot0)) call keep_compiler_quiet(pot0)
!
    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn,pot,pot0,grav,rmn)
!
!  Calculates gravity potential on a pencil.
!
!  13-nov-04/anders: coded
!
      real, dimension (nx) :: pot
      real, optional :: ymn,zmn,pot0
      real, optional, dimension (nx) :: xmn,rmn
      real, optional, dimension (nx,3) :: grav
!
      intent(in) :: xmn,ymn,zmn,rmn
      intent(out) :: pot
!
!  Calculate potential from master pencils defined in initialize_gravity.
!
      if (lxyzdependence) then
        pot=potx_xpencil(l1:l2)*zdep(n) &
           +poty_ypencil(m) &
           +potz_zpencil(n)*xdep(l1:l2)
      else
        pot=potx_xpencil(l1:l2)+poty_ypencil(m)+potz_zpencil(n)
      endif
!
      if (present(xmn)) call keep_compiler_quiet(xmn)
      if (present(ymn)) call keep_compiler_quiet(ymn)
      if (present(zmn)) call keep_compiler_quiet(zmn)
      if (present(pot0)) call keep_compiler_quiet(pot0)
      if (present(grav)) call keep_compiler_quiet(grav)
      if (present(rmn)) call keep_compiler_quiet(rmn)
!
    endsubroutine potential_penc
!***********************************************************************
    subroutine potential_point(x,y,z,r, pot,pot0, grav)
!
!  Calculates gravity potential in one point.
!
!  13-nov-04/anders: coded
!  24-oct-06/bing: added constant gravity profiles
!
      real :: pot
      real, optional :: x,y,z,r
      real, optional :: pot0,grav
      real :: potx_xpoint,poty_ypoint,potz_zpoint
      real :: prof,xdep,zdep
!
      potx_xpoint=0.0
      poty_ypoint=0.0
      potz_zpoint=0.0
!
      if (present(x)) then
        select case (gravx_profile)
        case ('zero')
          if (lroot) print*,'potential_point: no x-gravity'
        case ('const')
          potx_xpoint=-gravx*(x-xinfty)
        case ('linear_zdep')
          zdep=(1.+kappa_z1*z+.5*(kappa_z1*z)**2)
          potx_xpoint=0.5*(x**2-xinfty**2)*nux_epicycle**2*zdep
        case ('kepler')
          potx_xpoint=-gravx/x
        case default
          call fatal_error('potential_point', &
               'gravx_profile='//gravx_profile//' not implemented')
        endselect
      endif
!
      if (present(y)) then
        select case (gravy_profile)
        case ('zero')
          if (lroot) print*,'potential_point: no y-gravity'
        case ('const')
          poty_ypoint=-gravy*(y-yinfty)
        case default
          call fatal_error('potential_point', &
               'gravy_profile='//gravy_profile//' not implemented')
        endselect
      endif
!
      if (present(z)) then
        select case (gravz_profile)
        case ('zero')
          if (lroot) print*,'potential_point: no z-gravity'
        case ('const')
          potz_zpoint=-gravz*(z-zinfty)
        case ('linear')
          potz_zpoint=0.5*(z**2-zinfty**2)*nu_epicycle**2
        case ('spherical')
          potz_zpoint=0.5*nu_epicycle**2*z1**2*alog(1.0+(z/z1)**2)
        case ('linear_xdep')
          xdep=(1.+kappa_x1*x+.5*(kappa_x1*x)**2)
          potz_zpoint=0.5*(z**2-zinfty**2)*nu_epicycle**2*xdep
        case ('linear_smoothed')
          prof = 1. + (z/zref)**(2*n_pot)
          potz_zpoint = 0.5*(nu_epicycle*z)**2/prof**(1./n_pot)
        case ('tanh')
          potz_zpoint=gravz*zref*alog(cosh(z/zref))
        case default
          call fatal_error('potential_point', &
               'gravz_profile='//gravz_profile//' not implemented')
        endselect
      endif
!
      pot = potx_xpoint + poty_ypoint + potz_zpoint
!
      if (present(r)) call keep_compiler_quiet(r)
      if (present(grav)) call keep_compiler_quiet(grav)
      if (present(pot0)) call keep_compiler_quiet(pot0)
!
    endsubroutine potential_point
!***********************************************************************
    subroutine acceleration_penc(gg)
!
!  Calculates gravitational acceleration on a pencil.
!
!  21-apr-07/tobi: adapted from potential_penc
!
      real, dimension (:,:), intent (out) :: gg
!
!  Calculate acceleration from master pencils defined in initialize_gravity.
!
      if (size(gg,2)/=3) then
        call fatal_error('acceleration_penc','Expecting a 3-vector pencil')
      endif
!
!  Note: the following would not yet work if lxyzdependence is set to true.
!
      select case (size(gg,1))
      case (nx)
        gg(:,1) = gravx_xpencil(l1:l2)
        !ABlater: gg(:,1) = gravx_xpencil(l1:l2)*zdep(n)
      case (mx)
        gg(:,1) = gravx_xpencil
        !ABlater: gg(:,1) = gravx_xpencil*zdep
      case default
        call fatal_error('acceleration_penc','Wrong pencil size.')
      endselect
!
      gg(:,2) = gravy_ypencil(m)
      gg(:,3) = gravz_zpencil(n)
!
      !ABlater: gg(:,3) = gravz_zpencil(n)*xdep(:)
!
    endsubroutine acceleration_penc
!***********************************************************************
    subroutine acceleration_penc_1D(gr)
!
!  Calculates gravitational acceleration on a pencil.
!
!  21-apr-07/tobi: adapted from potential_penc
!
      real, dimension (nx), intent (out) :: gr
!
!  Calculate acceleration from master pencils defined in initialize_gravity.
!
      call fatal_error('acceleration_penc_1D','Not implemented')
!
      call keep_compiler_quiet(gr)
!
    endsubroutine acceleration_penc_1D
!***********************************************************************
    subroutine acceleration_point(x,y,z,r,g_r)
!
!  Gravity in one point
!
!  18-nov-08/wlad: coded.
!
      real :: g_r
      real, optional :: x,y,z,r
!
      intent(in)  :: x,y,z,r
      intent(out) :: g_r
!
      call fatal_error('gravity_simple','acceleration_point not implemented')
!
      g_r=0.0
!
      call keep_compiler_quiet(present(x))
      call keep_compiler_quiet(present(y))
      call keep_compiler_quiet(present(z))
      call keep_compiler_quiet(present(r))
!
    endsubroutine acceleration_point
!***********************************************************************
    subroutine read_gravity_init_pars(unit,iostat)
!
!  Read gravity init parameters.
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=grav_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=grav_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_gravity_init_pars
!***********************************************************************
    subroutine write_gravity_init_pars(unit)
!
!  Write gravity init parameters.
!
      integer, intent(in) :: unit
!
      write(unit,NML=grav_init_pars)
!
    endsubroutine write_gravity_init_pars
!***********************************************************************
    subroutine read_gravity_run_pars(unit,iostat)
!
!  Read gravity run parameters.
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=grav_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=grav_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_gravity_run_pars
!***********************************************************************
    subroutine write_gravity_run_pars(unit)
!
!  Write gravity run parameters.
!
      integer, intent(in) :: unit
!
      write(unit,NML=grav_run_pars)
!
    endsubroutine write_gravity_run_pars
!***********************************************************************
    subroutine rprint_gravity(lreset,lwrite)
!
!  Reads and registers print parameters relevant for gravity advance.
!
!  12-jun-04/axel: adapted from grav_z
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_epot=0
        idiag_epotmx=0; idiag_epotmy=0; idiag_epotmz=0
        idiag_epotuzmz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'epot',idiag_epot)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'epotmx', &
            idiag_epotmx)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'epotmy', &
            idiag_epotmy)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'epotmz', &
            idiag_epotmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'epotuzmz', &
            idiag_epotuzmz)
      enddo
!
!  Write column, idiag_XYZ, where our variable XYZ is stored.
!  IDL needs this even if everything is zero.
!
      if (lwr) then
        write(3,*) 'igg=',igg
        write(3,*) 'igx=',igx
        write(3,*) 'igy=',igy
        write(3,*) 'igz=',igz
      endif
!
    endsubroutine rprint_gravity
!***********************************************************************
    subroutine compute_gravity_star(f,wheat,luminosity,star_cte)
!
!  5-jan-10/boris: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: wheat, luminosity, star_cte
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(wheat,luminosity,star_cte)
!
    endsubroutine compute_gravity_star
!***********************************************************************
    subroutine get_xgravity(xgrav)
!      
!  Used from the initial conditions
!
!  04-oct-10/bing: coded
!
      real, dimension(mx) :: xgrav
!
      xgrav = gravx_xpencil
!
    endsubroutine get_xgravity
!***********************************************************************
endmodule Gravity
