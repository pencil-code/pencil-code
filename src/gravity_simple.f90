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
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
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
  double precision :: g_B, g_D, g_B_factor=1.0, g_D_factor=1.0
  real :: gravitational_const=0., mass_cent_body=0.
  real, dimension(mx) :: gravx_xpencil=0.0, potx_xpencil=0.0
  real, dimension(my) :: gravy_ypencil=0.0, poty_ypencil=0.0
  real, dimension(mz) :: gravz_zpencil=0.0, potz_zpencil=0.0
  real, dimension(mx) :: xdep=0.0
  real, dimension(mz) :: zdep=0.0
  real, parameter :: g_A_cgs=4.4e-9, g_C_cgs=1.7e-9
  real :: gravx=0.0, gravy=0.0, gravz=0.0
  real :: kx_gg=1.0, ky_gg=1.0, kz_gg=1.0, gravz_const=1.0, reduced_top=1.0
  real :: xgrav=impossible, ygrav=impossible, zgrav=impossible
  real :: xinfty=0.0, yinfty=0.0, zinfty=impossible, zclip=impossible
  real :: dgravx=0.0, pot_ratio=1.0
  real :: z1=0.0, z2=1.0, xref=0.0, zref=impossible, sphere_rad=0.0, g_ref=0.0
  real :: nu_epicycle=1.0, nu_epicycle2=1.0
  real :: nux_epicycle=0.0, nux_epicycle2=0.0
  real :: g_A, g_C, g_A_factor=1.0, g_C_factor=1.0
  real :: cs0hs=0.0, H0hs=0.0
  real :: potx_const=0.0, poty_const=0.0, potz_const=0.0
  integer :: n_pot=10
  integer :: n_adjust_sphersym=0
  character (len=labellen) :: gravx_profile='zero',gravy_profile='zero', &
                              gravz_profile='zero'
!
!  Parameters used by other modules (only defined for other gravities)
!
  logical :: lnumerical_equilibrium=.false.
  logical :: lxyzdependence=.false.
  logical :: lcalc_zinfty=.false.
  logical :: lboussinesq_grav=.false.
  logical :: ladjust_sphersym=.false.
 
  real :: g0=0.0
  real :: lnrho_bot=0.0, lnrho_top=0.0, ss_bot=0.0, ss_top=0.0
  real :: kappa_x1=0.0, kappa_x2=0.0, kappa_z1=0.0, kappa_z2=0.0
  real :: grav_tilt=0.0, grav_amp=0.0
!
  namelist /grav_init_pars/ &
      gravx_profile, gravy_profile, gravz_profile, gravx, gravy, gravz, &
      xgrav, ygrav, zgrav, kx_gg, ky_gg, kz_gg, dgravx, pot_ratio, z1, z2, &
      nux_epicycle, nu_epicycle, xref, zref, g_ref, sphere_rad, lnrho_bot, lnrho_top, &
      ss_bot, ss_top, lgravx_gas, lgravx_dust, lgravy_gas, lgravy_dust, &
      lgravz_gas, lgravz_dust, xinfty, yinfty, zinfty, lxyzdependence, &
      lcalc_zinfty, kappa_x1, kappa_x2, kappa_z1, kappa_z2, reduced_top, &
      lboussinesq_grav, n_pot, cs0hs, H0hs, grav_tilt, grav_amp, &
      potx_const,poty_const,potz_const, zclip, n_adjust_sphersym, gravitational_const, &
      mass_cent_body, g_A_factor, g_C_factor, g_B_factor, g_D_factor
!
  namelist /grav_run_pars/ &
      gravx_profile, gravy_profile, gravz_profile, gravx, gravy, gravz, &
      xgrav, ygrav, zgrav, kx_gg, ky_gg, kz_gg, dgravx, pot_ratio, &
      nux_epicycle, nu_epicycle, xref, zref, g_ref, sphere_rad, &
      lgravx_gas, lgravx_dust, lgravy_gas, lgravy_dust, &
      lgravz_gas, lgravz_dust, xinfty, yinfty, zinfty, lxyzdependence, &
      lcalc_zinfty, kappa_x1, kappa_x2, kappa_z1, kappa_z2, reduced_top, &
      lboussinesq_grav, n_pot, grav_tilt, grav_amp, &
      potx_const,poty_const,potz_const, zclip, n_adjust_sphersym, gravitational_const, &
      mass_cent_body, g_A_factor, g_C_factor, g_B_factor, g_D_factor
!
!  Diagnostic variables for print.in
! (needs to be consistent with reset list below)
!
  integer :: idiag_epot=0          ! DIAG_DOC: $\left<\varrho \Phi_{\rm grav}
                                   ! DIAG_DOC: \right>$ \quad(mean potential
                                   ! DIAG_DOC: energy)
  integer :: idiag_epottot=0       ! DIAG_DOC: $\int_V\varrho \Phi_{\rm grav}
                                   ! DIAG_DOC: dV$ \quad(total potential
                                   ! DIAG_DOC: energy)
  integer :: idiag_ugm=0           ! DIAG_DOC: $\left<\uv \cdot \gv\right>$
!
! xy averaged diagnostics given in xyaver.in written every it1d timestep
!
  integer :: idiag_epotmz=0        ! XYAVG_DOC: $\left<\varrho \Phi_{\rm grav}
                                   ! XYAVG_DOC: \right>_{xy}$
  integer :: idiag_epotuzmz=0      ! XYAVG_DOC: $\left<\varrho \Phi_{\rm grav}
                                   ! XYAVG_DOC: u_z \right>_{xy}$
                                   ! XYAVG_DOC: \quad(potential energy flux)
!
! xz averaged diagnostics given in xzaver.in
!
  integer :: idiag_epotmy=0        ! XZAVG_DOC: $\left<\varrho \Phi_{\rm grav}
                                   ! XZAVG_DOC: \right>_{xz}$
!
! yz averaged diagnostics given in yzaver.in
!
  integer :: idiag_epotmx=0        ! YZAVG_DOC: $\left<\varrho \Phi_{\rm grav}
                                   ! YZAVG_DOC: \right>_{yz}$
  integer :: idiag_epotuxmx=0      ! YZAVG_DOC: $\left<\varrho \Phi_{\rm grav}
                                   ! YZAVG_DOC: u_x \right>_{yz}$
                                   ! YZAVG_DOC: \quad(potential energy flux)
!
! z averaged diagnostics given in zaver.in
!
  integer :: idiag_epotmxy=0       ! ZAVG_DOC: $\left<\varrho \Phi_{\rm grav}
                                   ! ZAVG_DOC: \right>_{z}$
  integer :: idiag_epotuxmxy=0     ! ZAVG_DOC: $\left<\varrho \Phi_{\rm grav}
                                   ! ZAVG_DOC: u_x \right>_{z}$
                                   ! ZAVG_DOC: \quad(potential energy flux)
!
! work variables
!
  real, dimension(mx) :: gravx_xpencil_0
  real ::G4pi
!  
  contains
!***********************************************************************
    subroutine register_gravity
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
    subroutine initialize_gravity(f)
!
!  Calculate master pencils for gravity. These are put into gravity pencils
!  in the subroutine calc_pencils_grav.
!
!  12-nov-04/anders: coded, copied init conds from grav_x, grav_y and grav_y.
!   9-jun-15/MR: added parameter n_adjust_sphersym: if > 0 after each n_adjust_sphersym
!                timesteps the spherically symmetric part of gravity is adjusted
!                according to the actual density distribution.
!                Only in effect for spherical co-ordinates.
!  12-jun-15/MR: added (alternative) parameters gravitational_const and mass_cent_body for 
!                gravity adjustment.
!                For Kepler profile, gravitational_const is calculated from gravx and mass_cent_body.
!
      use SharedVariables, only: put_shared_variable, get_shared_variable
      use General, only: notanumber
      use Sub, only: cubic_step
      use Mpicomm, only: mpibcast_real, MPI_COMM_WORLD
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, pointer :: cs20,mpoly,gamma
!
      real, dimension (mz) :: prof
      real :: ztop
!
      character(len=*), parameter :: gravity_z_dat = 'prof_g.dat'
      real, dimension(:), allocatable :: grav_init_z
      integer :: alloc_err
!
!  Sanity check.
!
      if (lrun) then
        if (gravx_profile == 'zero' .and. &
            gravy_profile == 'zero' .and. &
            gravz_profile == 'zero') then
          call warning('initialize_gravity', &
              'You do not need gravity_simple for zero gravity...')
        endif
      endif
!
!  Possibility of specifying zref (if zinfty is not given).
!  In that case we need to know cs20 and mpoly from the EOS module.
!  We do this using shared variables.
!  Normally we give zinfty, but if this is not the case, we set it to zero.
!  zinfty has to be higher than the top of a polytropic atmosphere.
!
      if (lcalc_zinfty) then
        if (zinfty==impossible) then
          if (zref==impossible) then
            call fatal_error('initialize_gravity','zref=impossible')
          else
            call get_shared_variable('cs20',cs20,caller='initialize_gravity')
            call get_shared_variable('gamma',gamma)
            if (ldensity.and..not.lstratz) then
              call get_shared_variable('mpoly',mpoly)
            else
              if (lroot) call warning('initialize_eos','mpoly not obtained from density,'// &
                                      'set impossible')
              allocate(mpoly); mpoly=impossible
            endif
!
            zinfty=zref+cs20*(mpoly+1)/(-gamma*gravz)
            if (lroot) print*,'initialize_gravity: computed zinfty=',zinfty
          endif
        endif
      else
        if (zinfty==impossible) zinfty=0.
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
        call put_shared_variable('gravx', gravx, caller='initialize_gravity')
        call put_shared_variable('gravx_xpencil', gravx_xpencil)
!
      case ('const_tilt')
        gravx=grav_amp*sin(grav_tilt*pi/180.)
        if (lroot) print*,'initialize_gravity: constant gravx=', gravx
        gravx_xpencil=gravx
        potx_xpencil=-gravx*(x-xinfty)
        call put_shared_variable('gravx', gravx, caller='initialize_gravity')
!
! Linear gravity profile in x for planetary core dynamos
!
      case ('linear_x')
        if (lroot) print*,'initialize_gravity: linear x-grav, gravx=', gravx
        gravx_xpencil = -gravx*x
        potx_xpencil  = 0.5*gravx*(x**2-xinfty**2)
        call put_shared_variable('gravx', gravx, caller='initialize_gravity')
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
!AB: these potx_const values seem to be missing in corresponding entries for potx_xpoint.
!
      case ('tanh-pot')
        if (dgravx==0.) call fatal_error("initialize_gravity","dgravx=0 not OK")
        if (lroot) print*,'initialize_gravity: tanh x-grav, gravx=',gravx
        if (lroot) print*,'initialize_gravity: xgrav,dgravx=',xgrav,dgravx
        gravx=-alog(pot_ratio)/dgravx
        gravx_xpencil=gravx*.5/cosh((x-xgrav)/dgravx)**2
        potx_xpencil=-gravx*.5*(1.+tanh((x-xgrav)/dgravx))*dgravx + potx_const
!
      case ('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal x-grav, gravx=',gravx
        gravx_xpencil = -gravx*sin(kx_gg*x)
        potx_xpencil  = -gravx/kx_gg*cos(kx_gg*x) + potx_const
!
      case ('Baker74')
        !abag added, need to make adaptable to any width/position
        if (lroot) print*,'initialize_gravity: Baker_74=',gravx
        gravx_xpencil =(tanh((x+pi/3.)/0.1)+tanh(-(x-pi/3.)/0.1))/2.*&
            gravx*sin(2*(x-pi/2.))
        potx_xpencil=(tanh((x+pi/3.)/0.1)+tanh(-(x-pi/3.)/0.1))/2.*&
            gravx*(.5*cos(2*(x-pi/2.))-0.5) + potx_const
!
      case ('kepler')
        if (lroot) print*,'initialize_gravity: kepler x-grav, gravx=',gravx
        gravx_xpencil=-gravx/x**2
        potx_xpencil=-gravx/x + potx_const
        g0=gravx
        call put_shared_variable('gravx', gravx, caller='initialize_gravity')
        call put_shared_variable('gravx_xpencil', gravx_xpencil)
!
      case ('kepler_2d')
        if (lroot) print*,'initialize_gravity: kepler_2d x-grav, gravx=',gravx
        gravx_xpencil=-gravx/x
        potx_xpencil=-gravx*alog(x) + potx_const
        g0=gravx
        call put_shared_variable('gravx', gravx, caller='initialize_gravity')
        call put_shared_variable('gravx_xpencil', gravx_xpencil)
!
!  Convection zone model, normalized to the bottom of the domain
!
      case ('CZbot1')
        if (lroot) print*,'initialize_gravity: CZbot1 x-grav, gravx=',gravx
        gravx_xpencil=-gravx/x**2
        potx_xpencil=-gravx*(1./x-1./xyz0(1)) + potx_const
        g0=gravx
        call put_shared_variable('gravx', gravx, caller='initialize_gravity')
        call put_shared_variable('gravx_xpencil', gravx_xpencil)
!
!  Convection zone model, normalized to the middle of the domain
!
      case ('CZmid1')
        if (lroot) print*,'initialize_gravity: CZmid1 x-grav, gravx=',gravx
        gravx_xpencil=-gravx/x**2
        potx_xpencil=-gravx*(1./x-2./(xyz0(1)+xyz1(1))) + potx_const
        g0=gravx
        call put_shared_variable('gravx', gravx, caller='initialize_gravity')
        call put_shared_variable('gravx_xpencil', gravx_xpencil)
!
      ! solid sphere (homogenous) gravity profile
      ! 'xref' determines the location of the sphere border relative to the lower box boundary.
      ! 'g_ref' is the gravity acceleration at xref (implies the mass of the sphere).
      ! 'sphere_rad' is the radius of the sphere (independant from box coordinates).
      ! Together, sphere_rad and xref describe, how much of the sphere lies inside the box.
      case ('solid_sphere')
        if (lroot) print *, 'initialize_gravity: solid sphere xref=', xref, ', g_ref=', g_ref, ', sphere_rad=', sphere_rad
        where (x > xref)
          gravx_xpencil = g_ref * sphere_rad**2 / (x-xref+sphere_rad)**2
        else where (x < xref-2*sphere_rad)
          gravx_xpencil = -g_ref * sphere_rad**2 / (x-xref+sphere_rad)**2
        else where
          gravx_xpencil = g_ref * (x-xref+sphere_rad) / sphere_rad
        end where
        potx_xpencil=-gravx_xpencil*(x-xinfty)
!
      case ('loop')
        if (lroot) print*,'initialize_gravity: 1D loop, gravx=',gravx
        gravx_xpencil=cos(x / (2*xyz0(1)+Lxyz(1)) * pi)
        gravx_xpencil=gravx * gravx_xpencil
        potx_xpencil =sin(x / (2*xyz0(1)+Lxyz(1)) * pi)
        potx_xpencil =gravx*(2*xyz0(1)+Lxyz(1))/pi * potx_xpencil + potx_const
!
      case ('half-loop')
        if (lroot) print*,'initialize_gravity: 1D half-loop, gravx=',gravx
        gravx_xpencil=cos(x / (2*xyz0(1)+2*Lxyz(1)) * pi)
        gravx_xpencil=gravx * gravx_xpencil
        potx_xpencil =sin(x / (2*xyz0(1)+2*Lxyz(1)) * pi)
        potx_xpencil =gravx*(2*xyz0(1)+2*Lxyz(1))/pi * potx_xpencil + potx_const
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
        gravy_ypencil = gravy
        poty_ypencil  = -gravy*(y-yinfty)
!
      case ('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal y-grav, gravy=', gravy
        gravy_ypencil = -gravy*sin(ky_gg*y)
        poty_ypencil  = -gravy/ky_gg*cos(ky_gg*y) + poty_const
!
      case ('kepler')
        if (lroot) print*,'initialize_gravity: kepler gravy, y-grav=', gravy
        gravy_ypencil = -gravy/y**2
        poty_ypencil  = -gravy/y + poty_const
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
        if (zclip==impossible) then
          potz_zpencil=-gravz*(z-zinfty)
        else
          potz_zpencil=-gravz*min(z-zinfty,zclip-zinfty)
        endif
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
        lboussinesq_grav=.true.
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
        if (lroot) print *, 'initialize_gravity: solid sphere zref=', zref, ', g_ref=', g_ref, ', sphere_rad=', sphere_rad
        where (z > zref)
          gravz_zpencil = g_ref * sphere_rad**2 / (z-zref+sphere_rad)**2
        else where (z < zref-2*sphere_rad)
          gravz_zpencil = -g_ref * sphere_rad**2 / (z-zref+sphere_rad)**2
        else where
          gravz_zpencil = g_ref * (z-zref+sphere_rad) / sphere_rad
        end where
        potz_zpencil = -gravz_zpencil * (z - zinfty)
!        
      ! gravity profile provided as binary data
      ! 'zref' determines the location of the sphere border relative to the lower box boundary.
      ! 'g_ref' is the gravity acceleration at zref (implies the mass of the sphere).
      case ('profile')
        allocate(grav_init_z(mz), stat = alloc_err)
        if (alloc_err > 0) call fatal_error ('initialize_gravity', &
          'Could not allocate memory for gravity and gravity potential profiles', .true.)
        call read_grav_profile(gravity_z_dat, grav_init_z, real(unit_length/(unit_time*unit_time)), .false.)
        gravz_zpencil = grav_init_z
        potz_zpencil = -gravz_zpencil * (z - zinfty)
        gravz = grav_init_z(n1) !test code to use automatic calculation of hcond0 in boundcond.f90
!	
        deallocate(grav_init_z)
!
      case ('spherical')
        nu_epicycle2=nu_epicycle**2
        if (lroot) print*,'initialize_gravity: spherical z-grav, nu, z1=', nu_epicycle, z1
        gravz_zpencil=-nu_epicycle2*z/(1.0+(z/z1)**2)
        potz_zpencil=0.5*nu_epicycle2*z1**2*alog(1.0+(z/z1)**2) + potz_const
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
        potz_zpencil = 0.5*nu_epicycle2*z**2/prof**(1./n_pot) + potz_const
!
      case ('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal z-grav, gravz=', gravz
        gravz_zpencil = -gravz*sin(kz_gg*z)
        potz_zpencil = -gravz/kz_gg*cos(kz_gg*z) + potz_const
!
      case ('loop')
        if (lroot) print*,'initialize_gravity: loop, gravz=',gravz
        gravz_zpencil=cos(z / (2*xyz0(3)+Lxyz(3)) * pi)
        gravz_zpencil=-1.*gravz * gravz_zpencil
        potz_zpencil =sin(z / (2*xyz0(3)+Lxyz(3)) * pi)
        potz_zpencil =-1*gravz*(2*xyz0(3)+Lxyz(3))/pi * potz_zpencil + potz_const
!
      case ('kepler')
        if (lroot) print*,'initialize_gravity: kepler z-grav, gravz=', gravz
        gravz_zpencil = -gravz/z**2
        potz_zpencil  = -gravz/z + potz_const
!
      case ('Ferriere')
!
!  Set up physical units.
!
      if (unit_system=='cgs') then
        g_A = g_A_factor*g_A_cgs/unit_velocity*unit_time
        g_B = g_B_factor*g_B_cgs/unit_length
        g_C = g_C_factor*g_C_cgs/unit_velocity*unit_time
        g_D = g_D_factor*g_D_cgs/unit_length
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
      if (notanumber(gravx_xpencil)) &
        call fatal_error('initialize_gravity','found NaN or +/-Inf in gravx_xpencil')
      if (notanumber(gravy_ypencil)) &
        call fatal_error('initialize_gravity','found NaN or +/-Inf in gravy_ypencil')
      if (notanumber(gravz_zpencil)) &
        call fatal_error('initialize_gravity','found NaN or +/-Inf in gravz_zpencil')
!
      call put_shared_variable('nu_epicycle',nu_epicycle,caller='initialize_gravity')
      call put_shared_variable('gravz_zpencil',gravz_zpencil)
      if (lreference_state) call put_shared_variable('gravx', gravx)
!
      if (n_adjust_sphersym>0 .and. lspherical_coords) then

        if (gravx_profile=='kepler' .and. mass_cent_body/=0.) then
          if (gravitational_const/=0.) then
            if (lroot) call warning('initialize_gravity', 'central body mass is ignored')
          else
            if (ipx==0) gravitational_const = -gravx_xpencil(l1)*x(l1)**2/mass_cent_body
            if (nprocx>1) call mpibcast_real(gravitational_const,comm=MPI_COMM_WORLD)
          endif
        else
          if (gravitational_const<=0.) &
            call fatal_error('initialize_gravity', 'positive gravitational constant needed')
        endif

        G4pi=gravitational_const*4.*pi
        ladjust_sphersym=.true.
        gravx_xpencil_0 = gravx_xpencil

      endif
  
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine read_grav_profile(filename,profile,data_unit,llog)
!
!  Read vertical profile data.
!  Values are expected in SI units.
!
!  15-sept-2010/Bourdin.KIS: coded
!  17-feb-2021/joh-tsch: copied from solar_corona.f90
!
      use File_io, only: parallel_file_exists, file_size
      use Mpicomm, only: mpibcast_int, mpibcast_real
!
      character(len=*), intent(in) :: filename
      real, dimension(mz), intent(out) :: profile
      real, intent(in) :: data_unit
      logical, intent(in) :: llog
!
      real, dimension(:), allocatable :: data, data_z
      integer :: n_data
      integer :: iter
!
      integer, parameter :: unit=12
      integer :: len_double
!      
      integer :: alloc_err
!
!
      inquire (IOLENGTH=len_double) 1.0d0
!
      if (.not. parallel_file_exists (filename)) &
          call fatal_error ('read_grav_profile', "can't find "//filename)
      ! file access is only done on the MPI root rank
      if (lroot) then
        ! determine the number of data points in the profile
        n_data = (file_size (filename) - 2*2*4) / (len_double * 2)
      endif
      call mpibcast_int (n_data)
!
      ! allocate memory
      allocate (data(n_data), data_z(n_data), stat=alloc_err)
      if (alloc_err > 0) call fatal_error ('read_grav_profile', &
          'Could not allocate memory for data and its z coordinate', .true.)
!
      if (lroot) then
        ! read profile
        open (unit, file=filename, form='unformatted', recl=len_double*n_data)
        read (unit) data
        read (unit) data_z
        close (unit)
!
        if (llog) then
          ! convert data from logarithmic SI to logarithmic Pencil units
          data = data - alog (data_unit)
        else
          ! convert data from SI to Pencil units
          data = data / data_unit
        endif
!
        ! convert z coordinates from SI to Pencil units
        data_z = data_z / unit_length
      endif
!
      ! broadcast profile
      call mpibcast_real (data, n_data)
      call mpibcast_real (data_z, n_data)
!
      ! interpolate logarthmic data to Pencil grid profile
      call interpolate_grav_profile (data, data_z, n_data, profile)
!
      deallocate (data, data_z)
!
    endsubroutine read_grav_profile
!***********************************************************************
    subroutine interpolate_grav_profile(data,data_z,n_data,profile)
!
!  Interpolate profile data to Pencil grid.
!
!  15-sept-2010/Bourdin.KIS: coded
!  17-feb-2021/joh-tsch: copied from solar_corona.f90
!
      use General, only: itoa
!
      integer :: n_data
      real, dimension(n_data), intent(in) :: data, data_z
      real, dimension(mz), intent(out) :: profile
!
      integer :: i, j, num_over, num_below
!
!
      ! linear interpolation of data
      num_below = 0
      num_over = 0
      do j = 1, mz
        if (z(j) < data_z(1) ) then
          ! extrapolate linarily below bottom
          num_below = num_below + 1
          profile(j) = data(1) + (data(2)-data(1))/(data_z(2)-data_z(1)) * (z(j)-data_z(1))
        elseif (z(j) >= data_z(n_data)) then
          ! extrapolate linarily over top
          num_over = num_over + 1
          profile(j) = data(n_data) + (data(n_data)-data(n_data-1))/(data_z(n_data)-data_z(n_data-1)) * (z(j)-data_z(n_data))
        else
          do i = 1, n_data-1
            if ((z(j) >= data_z(i)) .and. (z(j) < data_z(i+1))) then
              ! y = m*(x-x1) + y1
              profile(j) = (data(i+1)-data(i)) / (data_z(i+1)-data_z(i)) * (z(j)-data_z(i)) + data(i)
              exit
            endif
          enddo
        endif
      enddo
!
      if (lfirst_proc_xy .and. (num_below > 0)) then
        call warning ("interpolate_grav_profile", &
            "extrapolated "//trim (itoa (num_below))//" grid points below bottom")
      endif
      if (lfirst_proc_xy .and. (num_over > 0)) then
        call warning ("interpolate_grav_profile", &
            "extrapolated "//trim (itoa (num_over))//" grid points over top")
      endif
!
    endsubroutine interpolate_grav_profile
!***********************************************************************
    subroutine gravity_sphersym(f)
!
!  Adjusts spherically symmetric part of gravitational acceleration
!  according to density.
!
!  9-jun-15/MR: coded
!
      use Sub, only: meanyz
      use Mpicomm, only: mpirecv_real, mpisend_real

      real, dimension(mx,my,mz,mfarray) :: f

      integer :: l, la
      real :: integ
      real, dimension(nx) :: rhomean

      call meanyz(f,ilnrho,rhomean,lexp=.not.ldensity_nolog)

      if (ipx==0) then
        integ=0.; la=l1+1          ! -> rectangle rule applied in integration below
      else
        call mpirecv_real(integ,xlneigh,iproc)
        la=l1
      endif

      do l=la,l2
        integ = integ+x(l)**2*rhomean(l-l1+1)*(x(l)-x(l-1))
        gravx_xpencil(l) = gravx_xpencil_0(l) - G4pi*integ/x(l)**2
      enddo

      if (ipx<nprocx-1) call mpisend_real(integ,xuneigh,xuneigh)

!if (ipy==0 .and. ipz==0) print'(a,38(f6.3,",",1x))', 'gravx_xpencil=', gravx_xpencil
    endsubroutine gravity_sphersym
!***********************************************************************
    subroutine set_consistent_gravity(ginput,gtype,gprofile,lsuccess)
!
!  This subroutine checks, if the gravity paramters as type, profile and values
!  are set consistently with initial condition for example.
!
!  ginput     =     value for the gravity, GM    : 4, 10, 200
!  gtype      =     type of gravity              : 'gravx','gravy','gravz'
!  gprofile   =     profile of the gravity       : 'kepler','const'
!  lsuccess   =     switch, if it was successful : .true., .false.
!
!  13-jun-12/dhruba+joern: coded
!
      real :: ginput
      character (len=labellen) :: gtype,gprofile
      character (len=labellen) :: gprof
      logical :: lsuccess
      logical :: lconsistent=.true.
!
! check for consistency
!
      gprof=trim(gprofile)
      select case(trim(gtype))
        case('gravx')
          if (gprof/=gravx_profile) then
            lconsistent=.false.
            gravx_profile=gprof
          endif
          if (gravx/=ginput) then
            lconsistent=.false.
            gravx=ginput
          endif
        case('gravy')
          if (gprof/=gravy_profile) then
            lconsistent=.false.
            gravy_profile=gprof
          endif
          if (gravy/=ginput) then
            lconsistent=.false.
            gravy=ginput
          endif
        case('gravz')
          if (gprof/=gravz_profile) then
            lconsistent=.false.
            gravz_profile=gprof
          endif
          if (gravz/=ginput) then
            lconsistent=.false.
            gravz=ginput
          endif
      case default
        call fatal_error('set_consistent_gravity','gtype does not match any, aborting')
      endselect
      lsuccess=.true.
!
! gravity parameters set consistently.
!
    endsubroutine set_consistent_gravity
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
    subroutine pencil_criteria_gravity
!
!  All pencils that the Gravity module depends on are specified here.
!
!  20-nov-04/anders: coded
!  20-jan-15/MR: pencil request for rho1 added when reference_state is used
!
      lpenc_requested(i_gg)=.true.
      if (lanelastic) lpenc_requested(i_rho_anel)=.true.
      if (lreference_state) lpenc_requested(i_rho1)=.true.
!
      if (idiag_epot/=0 .or. idiag_epot/=0 .or. idiag_epotmx/=0 .or. &
          idiag_epotmy/=0 .or. idiag_epotmz/=0 .or. &
          idiag_epottot/=0) lpenc_diagnos(i_epot)=.true.
      if (idiag_epotuxmx/=0 .or. idiag_epotuzmz/=0) then
        lpenc_diagnos(i_epot)=.true.
        lpenc_diagnos(i_uu)=.true.
      endif
      if (idiag_ugm/=0) then
        lpenc_diagnos(i_uu)=.true.
        lpenc_diagnos(i_gg)=.true.
      endif
      if (idiag_epotmxy/=0) then
        lpenc_diagnos2d(i_epot)=.true.
      endif
      if (idiag_epotuxmxy/=0) then
        lpenc_diagnos2d(i_epot)=.true.
        lpenc_diagnos2d(i_uu)=.true.
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
!  The special option lboussinesq_grav=T is applicable when |z|/H  << 1.
!  However, in the present formulation the resulting equations,
!  du/dt = -lnrho, and dlnrho/dt=-du/dz, lead to an instability
!  with the growth rate lambda = (1+i)*sqrt(k/2).
!
!  12-nov-04/anders: coded
!   5-dec-06/petri: added Boussinesq approximation
!  20-jan-15/MR: changes for use of reference state
!
      use Diagnostics
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      integer :: k
      real, dimension(nx,3) :: gg
      real, dimension(nx) :: refac
      real, dimension(:,:), pointer :: reference_state
!
!  Add gravity acceleration on gas.
!
      if (lhydro) then
        if (lboussinesq_grav) then
          if (lentropy) then
            if (lgravx_gas) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+p%ss*p%gg(:,1)
            if (lgravy_gas) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+p%ss*p%gg(:,2)
            if (lgravz_gas) df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+p%ss*p%gg(:,3)
          else
            if (headtt) print*,'duu_dt_grav: lboussinesq_grav w/o lentropy not ok!'
          endif
        else if (lanelastic) then
! Now works for the linear anelastic formulation only
                if (lgravx_gas) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+p%gg(:,1)*&
                                p%rho_anel
                if (lgravy_gas) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+p%gg(:,2)*&
                                p%rho_anel
                if (lgravz_gas) df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+p%gg(:,3)*&
                                 (-p%ss)
!                                p%rho_anel
        else
          if (lgravx_gas) gg(:,1)=p%gg(:,1)
          if (lgravy_gas) gg(:,2)=p%gg(:,2)
          if (lgravz_gas) gg(:,3)=p%gg(:,3)
!
! When reference state is used, gravity needs a correction factor rho'/rho=1-rho_0/rho.
! Note that the pencil case contains always the total quantities.
!
          if (lreference_state) then
            call get_shared_variable('reference_state',reference_state,caller='duu_dt_grav')  ! shouldn't be within the mn loop
            refac=1.-reference_state(:,iref_rho)*p%rho1
            if (lgravx_gas) gg(:,1)=gg(:,1)*refac
            if (lgravy_gas) gg(:,2)=gg(:,2)*refac
            if (lgravz_gas) gg(:,3)=gg(:,3)*refac
          endif
!
          if (lxyzdependence) then
            if (lgravx_gas) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+gg(:,1)*zdep(n)
            if (lgravy_gas) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+gg(:,2)
            if (lgravz_gas) df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+gg(:,3)*xdep(l1:l2)
          else
            if (lgravx_gas) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+gg(:,1)
            if (lgravy_gas) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+gg(:,2)
            if (lgravz_gas) df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+gg(:,3)
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
        if (idiag_ugm/=0) &
          call sum_mn_name(p%uu(:,1)*p%gg(:,1) + &
          p%uu(:,2)*p%gg(:,2) + &
          p%uu(:,3)*p%gg(:,3),idiag_ugm)
        if (idiag_epottot/=0) call integrate_mn_name(p%epot,idiag_epottot)
      endif
!
!  Gravity 1-D diagnostics.
!
      if (l1davgfirst) then
        call yzsum_mn_name_x(p%epot,idiag_epotmx)
        call yzsum_mn_name_x(p%epot*p%uu(:,1),idiag_epotuxmx)
        call xzsum_mn_name_y(p%epot,idiag_epotmy)
        call xysum_mn_name_z(p%epot,idiag_epotmz)
        call xysum_mn_name_z(p%epot*p%uu(:,3),idiag_epotuzmz)
      endif
!
!  Gravity 2-D diagnostics.
!
      if (l2davgfirst) then
        call zsum_mn_name_xy(p%epot,idiag_epotmxy)
        call zsum_mn_name_xy(p%epot*p%uu(:,1),idiag_epotuxmxy)
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine gravity_after_boundary(f)
!
!  For actions outside mn-loop. 
!  At the moment only adjustment of spherically symmetric gravity.
!
!  9-jun-15/MR: coded
!
      real, dimension(mx,my,mz,mfarray) :: f

      if (lfirst.and.ladjust_sphersym) then
        if ( mod(it, n_adjust_sphersym) == 0 ) call gravity_sphersym(f)
      endif

    endsubroutine gravity_after_boundary
!***********************************************************************
    subroutine potential_global(pot,pot0)
!
!  Calculates gravity potential globally
!
!  13-nov-04/anders: coded
!
      real, dimension (:,:,:) :: pot
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
      real, dimension (:) :: pot
      real, optional :: ymn,zmn,pot0
      real, optional, dimension (size(pot)) :: xmn,rmn
      real, optional, dimension (size(pot),3) :: grav
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
!  Selection different potentials at one point.
!  These entries should match those for potx_xpencil.
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
        case ('CZbot1')
          potx_xpoint=-gravx*(1./x-1./xyz0(1))
        case ('CZmid1')
          potx_xpoint=-gravx*(1./x-2./(xyz0(1)+xyz1(1)))
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
          if (zclip==impossible) then
            potz_zpoint=-gravz*(z-zinfty)
          else
            potz_zpoint=-gravz*max(z-zinfty,zclip-zinfty)
          endif
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
    subroutine read_gravity_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=grav_init_pars, IOSTAT=iostat)
!
    endsubroutine read_gravity_init_pars
!***********************************************************************
    subroutine write_gravity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=grav_init_pars)
!
    endsubroutine write_gravity_init_pars
!***********************************************************************
    subroutine read_gravity_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=grav_run_pars, IOSTAT=iostat)
!
    endsubroutine read_gravity_run_pars
!***********************************************************************
    subroutine write_gravity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=grav_run_pars)
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
      use FArrayManager, only: farray_index_append
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez, inamexy
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_epot=0; idiag_epottot=0; idiag_ugm=0
        idiag_epotmx=0; idiag_epotuxmx=0; idiag_epotmy=0; idiag_epotmz=0;
        idiag_epotmxy=0; idiag_epotuzmz=0; idiag_epotuxmxy=0
      endif
!
!  iname runs through all possible names that may be listed in print.in.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'epot',idiag_epot)
        call parse_name(iname,cname(iname),cform(iname),'epottot',idiag_epottot)
        call parse_name(iname,cname(iname),cform(iname),'ugm',idiag_ugm)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'epotmx', &
            idiag_epotmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'epotuxmx', &
            idiag_epotuxmx)
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
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy), &
            'epotmxy', idiag_epotmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy), &
            'epotuxmxy', idiag_epotuxmxy)
      enddo
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
    logical function is_constant_zgrav()
!
!  15-apr-15/MR: coded
!
      is_constant_zgrav = gravx_profile=='zero'.and.gravy_profile=='zero'.and.gravz_profile=='const'

    endfunction is_constant_zgrav
!***********************************************************************
endmodule Gravity
