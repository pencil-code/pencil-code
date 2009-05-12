! $Id$
!
!
!  This module takes care of simple types of gravity, i.e. where
!    gx=gx(x) or gy=gy(y) or gz=gz(z)
!  Here the gravity master pencils gravx_xpencil, gravy_ypencil and
!  gravz_zpencil only need to be calculated once, and then these can
!  simply be added to the equations of motion again and again.
!
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
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
  real, dimension(mx) :: gravx_xpencil=0.,potx_xpencil=0.
  real, dimension(my) :: gravy_ypencil=0.,poty_ypencil=0.
  real, dimension(mz) :: gravz_zpencil=0.,potz_zpencil=0.
  real :: gravx=0.,gravy=0.,gravz=0.
  real :: kx_gg=1.,ky_gg=1.,kz_gg=1.,grav_const=1.,reduced_top=1.
  real :: xgrav=impossible,ygrav=impossible,zgrav=impossible
  real :: xinfty=0.,yinfty=0.,zinfty=0.
  real :: dgravx=0.,pot_ratio=1.
  real :: z1=0.,z2=1.,zref=0.,qgshear=1.5
  real :: nu_epicycle=1.,nu_epicycle2=1.
  real :: r0_pot=0.    ! peak radius for smoothed potential
  real :: g_A, g_C
  real, parameter :: g_A_cgs=4.4e-9, g_C_cgs=1.7e-9
  double precision :: g_B, g_D
  double precision, parameter :: g_B_cgs=6.172D20, g_D_cgs=3.086D21
  real :: cs0hs, H0hs
  integer :: n_pot=10  ! exponent for smoothed potential
  character (len=labellen) :: gravx_profile='zero',gravy_profile='zero'
  character (len=labellen) :: gravz_profile='zero'
!
!  Parameters used by other modules (only defined for other gravities)
!
  logical :: lnumerical_equilibrium=.false.
  logical :: lboussinesq=.false.
  real :: g0=0.
  real :: lnrho_bot=0.,lnrho_top=0.,ss_bot=0.,ss_top=0.
  character (len=labellen) :: grav_profile='const'
!
  namelist /grav_init_pars/ &
       gravx_profile,gravy_profile,gravz_profile,gravx,gravy,gravz, &
       xgrav,ygrav,zgrav,kx_gg,ky_gg,kz_gg,dgravx,nu_epicycle,pot_ratio, &
       z1,z2,zref,lnrho_bot,lnrho_top,ss_bot,ss_top, &
       lgravx_gas,lgravx_dust,lgravy_gas,lgravy_dust,lgravz_gas,lgravz_dust, &
       xinfty,yinfty,zinfty, &
       reduced_top,lboussinesq,grav_profile,n_pot, &
       cs0hs,H0hs
!
  namelist /grav_run_pars/ &
       gravx_profile,gravy_profile,gravz_profile,gravx,gravy,gravz, &
       xgrav,ygrav,zgrav,kx_gg,ky_gg,kz_gg,dgravx,nu_epicycle,pot_ratio, &
       lgravx_gas,lgravx_dust,lgravy_gas,lgravy_dust,lgravz_gas,lgravz_dust, &
       xinfty,yinfty,zinfty, &
       zref,reduced_top,lboussinesq,grav_profile,n_pot
!
  integer :: idiag_epot=0
!
  contains
!***********************************************************************
    subroutine register_gravity()
!
!  Initialise gravity variables (currently none)
!
!  12-nov-04/anders: coded
!
      use Mpicomm
      use Sub
!
!  Identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id$")
!
!  Set lgrav and lgravz (the latter for backwards compatibility)
!  Set lgravz only when gravz_profile is set.
!
      lgrav=.true.
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
      use SharedVariables
      use Sub, only: notanumber, cubic_step
!
      real, dimension(mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      real, dimension (mz) :: prof
      real :: ztop
      integer :: ierr
!
!  Sanity check
!
      if (.not. lstarting) then
        if (gravx_profile == 'zero' .and. &
            gravy_profile == 'zero' .and. &
            gravz_profile == 'zero') then
          call fatal_error('initialize_gravity', &
              "You don't need gravity_simple for zero gravity...")
        endif
      endif
!
!  Different x-gravity profiles
!
      if (gravx/=0) lgravx=.true.
      select case (gravx_profile)

      case('zero')
        if (lroot) print*,'initialize_gravity: no x-gravity'
        gravx_xpencil=0.
        lgravx_gas=.false.
        lgravx_dust=.false.

      case('const')
        if (lroot) print*,'initialize_gravity: constant x-grav=',gravx
        gravx_xpencil=gravx
        potx_xpencil=-gravx*(x-xinfty)
!
!  tanh profile
!  for isothermal EOS, we have 0=-cs2*dlnrho+gravx
!  pot_ratio gives the resulting ratio in the density.
!
      case('tanh-pot')
        if (dgravx==0.) call fatal_error("initialize_gravity","dgravx=0 not OK")
        if (lroot) print*,'initialize_gravity: tanh x-grav, gravx=',gravx
        if (lroot) print*,'initialize_gravity: xgrav,dgravx=',xgrav,dgravx
        gravx=-alog(pot_ratio)/dgravx
        gravx_xpencil=gravx*.5/cosh((x-xgrav)/dgravx)**2
        potx_xpencil=-gravx*.5*(1.+tanh((x-xgrav)/dgravx))*dgravx

      case('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal x-grav, gravx=',gravx
        gravx_xpencil = -gravx*sin(kx_gg*x)

      case('kepler')
        if (lroot) print*,'initialize_gravity: kepler x-grav, gravx=',gravx
        gravx_xpencil=-gravx/x**2
        potx_xpencil=-gravx/x
        g0=gravx

      case default
        if (lroot) print*, &
            'initialize_gravity: unknown gravx_profile ', gravx_profile
        call fatal_error('initialize_gravity','chosen gravx_profile not valid')

      endselect
!
!  Different y-gravity profiles
!
      select case (gravy_profile)

      case('zero')
        if (lroot) print*,'initialize_gravity: no y-gravity'
        gravy_ypencil=0.
        lgravy_gas=.false.
        lgravy_dust=.false.

      case('const')
        if (lroot) print*,'initialize_gravity: constant y-grav=', gravy
        gravy_ypencil=gravy
        poty_ypencil=-gravy*(y-yinfty)

      case('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal y-grav, gravy=', gravy
        gravy_ypencil = -gravy*sin(ky_gg*y)

      case('kepler')
        if (lroot) print*,'initialize_gravity: kepler gravy, y-grav=', gravy
        gravy_ypencil = -gravy/y**2

      case default
        if (lroot) print*, &
            'initialize_gravity: unknown gravy_profile ', gravy_profile
        call fatal_error('initialize_gravity','chosen gravy_profile not valid')

      endselect
!
!  Different z-gravity profiles
!  Set lgravz=T only when gravz_profile is not zero
!
      if (gravz_profile/='zero') lgravz=.true.

      select case (gravz_profile)

      case('zero')
        if (lroot) print*,'initialize_gravity: no z-gravity'
        gravz_zpencil=0.
        lgravz_gas=.false.
        lgravz_dust=.false.

      case('const')
        if (lroot) print*,'initialize_gravity: constant gravz=', gravz
        gravz_zpencil=gravz
        potz_zpencil=-gravz*(z-zinfty)

      case('tanh')
        if (lroot) print*,'initialize_gravity: tanh gravz=', gravz
        gravz_zpencil=-gravz*tanh(z/zref)
        potz_zpencil=gravz*zref*alog(cosh(z/zref))

      case('boussinesq')
        if (lroot) print*,'initialize_gravity: boussinesq gravz=', gravz
        gravz_zpencil=gravz
        potz_zpencil=-gravz*(z-zinfty)
        lboussinesq=.true.

      case('const_zero')  !  Const. gravity acc. (but zero for z>zgrav)
        if (headtt) print*,'initialize_gravity: const_zero gravz=', gravz
        if (zgrav==impossible .and. lroot) &
            print*,'initialize_gravity: zgrav is not set!'
        do n=n1,n2
          if (z(n)<=zgrav) gravz_zpencil(n) = gravz
        enddo

      case('linear')      !  Linear gravity profile (for accretion discs)
        nu_epicycle2=nu_epicycle**2
        if (lroot) print*,'initialize_gravity: linear z-grav, nu=', nu_epicycle
        gravz_zpencil = -nu_epicycle2*z
        potz_zpencil=0.5*nu_epicycle2*(z**2-zinfty**2)

      case('linear_smoothed')
        nu_epicycle2=nu_epicycle**2
        if (lroot) print*,'initialize_gravity: linear z-grav, '// &
                         'smoothed to zero at top/bottom, nu=', nu_epicycle
        prof = 1. + (z/zref)**(2*n_pot)
        gravz_zpencil = -nu_epicycle2*z/prof**(1./n_pot+1.)
        potz_zpencil = 0.5*nu_epicycle2*z**2/prof**(1./n_pot)

      case('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal z-grav, gravz=', gravz
        gravz_zpencil = -gravz*sin(kz_gg*z)
        potz_zpencil = -gravz/kz_gg*cos(kz_gg*z)

      case('kepler')
        if (lroot) print*,'initialize_gravity: kepler z-grav, gravz=', gravz
        gravz_zpencil = -gravz/z**2

      case('Ferriere')
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
!  gravity profile from K. Ferriere, ApJ 497, 759, 1998, eq (34)
!   at solar radius.  (for interstellar runs)
!
!  nb: 331.5 is conversion factor: 10^-9 cm/s^2 -> kpc/Gyr^2)  (/= 321.1 ?!?)
!AB: These numbers should be inserted in the appropriate units.
!AB: As it is now, it can never make much sense.
        gravz_zpencil = -(g_A*z/sqrt(z**2+g_B**2) + g_C*z/g_D)

      case('Galactic-hs')
        if (lroot) print*,'Galactic hydrostatic equilibrium gravity profile'
        if (lroot.and.(cs0hs==0.or.H0hs==0)) &
            call fatal_error('initialize-gravity', &
            'Set cs0hs and H0hs in grav_init_pars!')
        gravz_zpencil = -z*(cs0hs/H0hs)**2/sqrt(1 + (z/H0hs)**2)

      case('reduced_top')
        if (lroot) print*,'initialize_gravity: reduced, gravz=',gravz
        if (zgrav==impossible.and.lroot) print*,'zgrav is not set!'
        ztop = xyz0(3)+Lxyz(3)
        prof = cubic_step(z,(zgrav+ztop)/2,(ztop-zgrav)/2)
        gravz_zpencil = (1 - prof*(1-reduced_top))*gravz

      case default
        if (lroot) print*, &
            'initialize_gravity: unknown gravz_profile ', gravz_profile
        call fatal_error('initialize_gravity','chosen gravz_profile not valid')

      endselect
!
!  Sanity check
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
      call put_shared_variable('nu_epicycle',nu_epicycle,ierr)
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine init_gg(f)
!
!  Initialise gravity; called from start.f90
!
!  12-nov-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Don't do anything
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
!
      if (idiag_epot/=0) lpenc_diagnos(i_epot)=.true.
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
      if (lpencil(i_epot)) p%epot = p%rho*(potx_xpencil(l1:l2) + poty_ypencil(m) + potz_zpencil(n))
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_gravity
!***********************************************************************
    subroutine duu_dt_grav(f,df,p)
!
!  Add gravitational acceleration to gas and dust
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
!  Add gravity acceleration on gas
!
      if (lhydro) then
        if (lboussinesq) then
          if (lentropy) then
            if (lgravx_gas) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+p%ss*p%gg(:,1)
            if (lgravy_gas) df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+p%ss*p%gg(:,2)
            if (lgravz_gas) df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+p%ss*p%gg(:,3)
          else
            if (headtt) print*,'duu_dt_grav: lbounssinesq w/o lentropy not ok!'
          endif
        else
          if (lgravx_gas) df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + p%gg(:,1)
          if (lgravy_gas) df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + p%gg(:,2)
          if (lgravz_gas) df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + p%gg(:,3)
        endif
      endif
!
!  Add gravity acceleration on dust
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
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
      call keep_compiler_quiet(pot0)
!
    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn,pot,pot0,grav,rmn)
!
!  Calculates gravity potential on a pencil
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
!  Calculate potential from master pencils defined in initialize_gravity
!
      pot = potx_xpencil(l1:l2) + poty_ypencil(m) + potz_zpencil(n)
!
      call keep_compiler_quiet(xmn)
      call keep_compiler_quiet(ymn)
      call keep_compiler_quiet(zmn)
      call keep_compiler_quiet(pot0)
      call keep_compiler_quiet(grav)
      call keep_compiler_quiet(rmn)
!
    endsubroutine potential_penc
!***********************************************************************
    subroutine potential_point(x,y,z,r, pot,pot0, grav)
!
!  Calculates gravity potential in one point
!
!  13-nov-04/anders: coded
!  24-oct-06bing: added constant gravity profiles
!
      real :: pot
      real, optional :: x,y,z,r
      real, optional :: pot0,grav
      real :: potx_xpoint,poty_ypoint,potz_zpoint
      real :: prof
      integer :: i
!
      potx_xpoint=0.
      poty_ypoint=0.
      potz_zpoint=0.
!
      if (present(x)) then
        select case (gravx_profile)
        case('zero')
          if (lroot) print*,'potential_point: no x-gravity'
        case('const')
          potx_xpoint=-gravx*(x-xinfty)
        case('kepler')
          potx_xpoint=-gravx/x
        case default
          call fatal_error('potential_point', &
               'gravx_profile='//gravx_profile//' not implemented')
        endselect
      endif
!
      if (present(y)) then
        select case (gravy_profile)
        case('zero')
          if (lroot) print*,'potential_point: no y-gravity'
        case('const')
          poty_ypoint=-gravy*(y-yinfty)
        case default
          call fatal_error('potential_point', &
               'gravy_profile='//gravy_profile//' not implemented')
        endselect
      endif
!
      if (present(z)) then
        select case (gravz_profile)
        case('zero')
          if (lroot) print*,'potential_point: no z-gravity'
        case('const')
          potz_zpoint=-gravz*(z-zinfty)
        case('linear')
          potz_zpoint=0.5*(z**2-zinfty**2)*nu_epicycle**2
        case('linear_smoothed')
          prof = 1. + (z/zref)**(2*n_pot)
          potz_zpoint = 0.5*(nu_epicycle*z)**2/prof**(1./n_pot)
        case('tanh')
          potz_zpoint=gravz*zref*alog(cosh(z/zref))
        case default
          call fatal_error('potential_point', &
               'gravz_profile='//gravz_profile//' not implemented')
        endselect
      endif
!
      pot = potx_xpoint + poty_ypoint + potz_zpoint
!
    endsubroutine potential_point
!***********************************************************************
    subroutine acceleration_penc(gg)
!
!  Calculates gravitational acceleration on a pencil
!
!  21-apr-07/tobi: adapted from potential_penc
!
      real, dimension (:,:), intent (out) :: gg
!
!  Calculate acceleration from master pencils defined in initialize_gravity
!
      if (size(gg,2) /= 3) then
        call fatal_error("acceleration_penc","Expecting a 3-vector pencil.")
      endif

      select case (size(gg,1))
      case (nx)
        gg(:,1) = gravx_xpencil(l1:l2)
      case (mx)
        gg(:,1) = gravx_xpencil
      case default
        call fatal_error("acceleration_penc","Wrong pencil size.")
      endselect

      gg(:,2) = gravy_ypencil(m)
      gg(:,3) = gravz_zpencil(n)
!
    endsubroutine acceleration_penc
!***********************************************************************
    subroutine acceleration_penc_1D(gr)
!
!  Calculates gravitational acceleration on a pencil
!
!  21-apr-07/tobi: adapted from potential_penc
!
      real, dimension (nx), intent (out) :: gr
!
!  Calculate acceleration from master pencils defined in initialize_gravity
!
      call fatal_error("acceleration_penc_1D","Not implemented")
!
      call keep_compiler_quiet(gr)
!
    endsubroutine acceleration_penc_1D
!***********************************************************************
    subroutine acceleration_point(x,y,z,r,g_r)
!
!  Gravity in one point
!
!  18-nov-08/wlad: coded
!
      use Mpicomm, only: stop_it
!
      real :: g_r
      real, optional :: x,y,z,r
!
      intent(in)  :: x,y,z,r
      intent(out) :: g_r
!
      call stop_it("gravity_simple: acceleration_point not implemented")
!
      g_r = 0.
!
      call keep_compiler_quiet(x)
      call keep_compiler_quiet(y)
      call keep_compiler_quiet(z)
      call keep_compiler_quiet(r)
      call keep_compiler_quiet(g_r)
!
    endsubroutine acceleration_point
!***********************************************************************
    subroutine read_gravity_init_pars(unit,iostat)
!
!  Read gravity init parameters
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
    endsubroutine read_gravity_init_pars
!***********************************************************************
    subroutine write_gravity_init_pars(unit)
!
!  Write gravity init parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=grav_init_pars)
!
    endsubroutine write_gravity_init_pars
!***********************************************************************
    subroutine read_gravity_run_pars(unit,iostat)
!
!  Read gravity run parameters
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
    endsubroutine read_gravity_run_pars
!***********************************************************************
    subroutine write_gravity_run_pars(unit)
!
!  Write gravity run parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=grav_run_pars)
!
    endsubroutine write_gravity_run_pars
!***********************************************************************
    subroutine rprint_gravity(lreset,lwrite)
!
!  reads and registers print parameters relevant for gravity advance
!  dummy routine
!
!  12-jun-04/axel: adapted from grav_z
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_epot=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'epot',idiag_epot)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!  idl needs this even if everything is zero
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
endmodule Gravity
