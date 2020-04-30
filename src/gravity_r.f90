! $Id: gravity_r.f90,v 1.1 2018/08/24 15:48:10 wlyra Exp $
!
!  Radial gravity
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lgrav = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
! MGLOBAL CONTRIBUTION 3
!
! PENCILS PROVIDED gg(3)
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
  ! coefficients for potential
  real, dimension (5,ninit) :: cpot=0. !=(/ 0., 0., 0., 0., 0. /)
  real, dimension(ninit) :: g01=0.,rpot=0.
  real :: lnrho_bot,lnrho_top,ss_bot,ss_top
  real :: gravz_const=1.,reduced_top=1.
  real :: g0=0.
  real, target :: g1=0.,rp1=1.0,rp1_smooth=impossible
  real, target :: gsum=0.
  real :: rp1_smooth1,frac_smooth=1.0
  real :: r0_pot=0.,r1_pot1=0.    ! peak radius for smoothed potential
  real :: n_pot=10,n_pot1=10   ! exponent for smoothed potential
  real :: qgshear=1.5  ! (global) shear parameter
                       !     1.5 for Keplerian disks, 1.0 for galaxies
!
  character (len=labellen), dimension(ninit) :: ipotential='zero'
!
  ! variables for compatibility with grav_z (used by Entropy and Density):
  real :: z1,z2,zref,zgrav,gravz=0.,zinfty
  real :: nu_epicycle=1.0
  real, target :: t_ramp_mass=impossible
  real :: t1_ramp_mass
  character (len=labellen) :: gravz_profile='zero'
  character (len=labellen) :: ipotential_secondary='plummer'
!
  logical :: lnumerical_equilibrium=.false.
  logical :: lgravity_gas=.true.
  logical :: lgravity_neutrals=.true.
  logical :: lgravity_dust=.true.
  logical :: lindirect_terms=.true.
  logical, target :: lramp_mass=.false.
  integer :: iglobal_gg=0
  logical, target :: lsecondary_wait=.false.
  logical :: lcoriolis_force_gravity=.true.
  logical :: lcentrifugal_force_gravity=.true.
  real, target :: t_start_secondary = -impossible
!
  namelist /grav_init_pars/ &
      ipotential,g0,r0_pot,r1_pot1,n_pot,n_pot1,lnumerical_equilibrium, &
      qgshear,lgravity_gas,g01,rpot,gravz_profile,gravz,nu_epicycle, &
      lgravity_neutrals,g1,rp1_smooth,lindirect_terms,lramp_mass,t_ramp_mass,&
      ipotential_secondary,lsecondary_wait,t_start_secondary, &
      lcoriolis_force_gravity,lcentrifugal_force_gravity,frac_smooth
!
  namelist /grav_run_pars/ &
      ipotential,g0,r0_pot,n_pot,lnumerical_equilibrium, &
      qgshear,lgravity_gas,g01,rpot,gravz_profile,gravz,nu_epicycle, &
      lgravity_neutrals,g1,rp1_smooth,lindirect_terms,lramp_mass,t_ramp_mass, &
      ipotential_secondary,lsecondary_wait,t_start_secondary, &
      lcoriolis_force_gravity,lcentrifugal_force_gravity,frac_smooth
!
  integer :: idiag_torque=0
!
  contains
!***********************************************************************
    subroutine register_gravity()
!
!  initialise gravity flags
!
!  10-jan-02/wolf: coded
!
!  Identify version number.
!
      if (lroot) call svn_id("$Id: gravity_r.f90,v 1.1 2018/08/24 15:48:10 wlyra Exp $")
!
      lgravr=.true.
      lgravr_gas =.true.
      lgravr_dust=.true.
      lgravr_neutrals=.true.
!
    endsubroutine register_gravity
!***********************************************************************
    subroutine initialize_gravity(f)
!
!  Set up cpot according to the value of ipotential, and initialize the
!  global variable gg (gravity field).
!  Needed by both start.f90 and run.f90
!
!  16-jul-02/wolf: coded
!  22-nov-02/tony: renamed
!  15-mar-07/wlad: made it coordinate-independent
!
      use Sub, only: poly, step, get_radial_distance
      use Mpicomm
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: gg_mn=0.0
      real, dimension (nx)   :: g_r,rr_mn,rr_sph,rr_cyl,pot
      logical :: lpade=.true. ! set to false for 1/r potential
      integer :: j,ierr
!
!  Pre-calculate the angular frequency of the rotating frame in the case a
!  secondary stationary body is added to the simulation.
!
      if (lcorotational_frame) then 
        gsum=g0+g1
        Omega_corot = sqrt(gsum/rcorot**3)
        rp1=rcorot
      else
        gsum=g0
        if (g1/=0) call fatal_error("initialize_gravity",&
             "companion gravity coded only for corotational frame")
      endif
!
!  Smoothing radius: default to Hill radius
!
      if (rp1_smooth == impossible) &
           rp1_smooth   = frac_smooth * rp1 * (g1/3.)**(1./3)
      if (rp1_smooth /= 0) then
        rp1_smooth1  = 1./rp1_smooth
      else
        rp1_smooth1 = 0.
      endif
!
!  Share variables related to corotational frame to modules that may need it
!
      call put_shared_variable('gsum',gsum,ierr)
      if (ierr/=0) call fatal_error('initialize_gravity', &
          'there was a problem when putting gsum')
!
      call put_shared_variable('g0',g0,ierr)
      if (ierr/=0) call fatal_error('initialize_gravity', &
          'there was a problem when putting g0')
!
      call put_shared_variable('g1',g1,ierr)
      if (ierr/=0) call fatal_error('initialize_gravity', &
          'there was a problem when putting g1')
!
      call put_shared_variable('rp1',rp1,ierr)
      if (ierr/=0) call fatal_error('initialize_gravity', &
          'there was a problem when putting rp1')
!
      call put_shared_variable('rp1_smooth',rp1_smooth,ierr)
      if (ierr/=0) call fatal_error('initialize_gravity', &
          'there was a problem when putting rp1_smooth')
!
      call put_shared_variable('lramp_mass',lramp_mass,ierr)
      if (ierr/=0) call fatal_error('initialize_gravity', &
          'there was a problem when putting lramp_mass')
!
      call put_shared_variable('t_ramp_mass',t_ramp_mass,ierr)
      if (ierr/=0) call fatal_error('initialize_gravity', &
          'there was a problem when putting t_ramp_mass')
!
      call put_shared_variable('lsecondary_wait',lsecondary_wait,ierr)
      if (ierr/=0) call fatal_error('initialize_gravity', &
          'there was a problem when putting lsecondary_wait')
!
      call put_shared_variable('t_start_secondary',t_start_secondary,ierr)
      if (ierr/=0) call fatal_error('initialize_gravity', &
          'there was a problem when putting t_start_secondary')
!
!  Shortcut for optimization
!
      if (t_ramp_mass/=impossible) t1_ramp_mass=1./t_ramp_mass
!
!  for lpade=.true. set coefficients for potential (coefficients a0, a2, a3,
!  b2, b3) for the rational approximation
!
!              a_0   +   a_2 r^2 + a_3 r^3
!    Phi(r) = ---------------------------------------
!               1    +   b_2 r^2 + b_3 r^3 + a_3 r^4
!
      if (lnumerical_equilibrium) then
!
        if (lroot) then
          print*,'initialize_gravity: numerical exact equilibrium'
          print*,' - gravity will be calculated in density module'
        endif
!
      else
!
!  Initialize gg, so we can later retrieve gravity via get_global.
!  Ensure the reserved array slots are initialized to zero, so we can add
!  ninit different gravity fields.
!
        if (igg==0) call farray_register_global('global_gg',iglobal_gg,vector=3)
        f(l1:l2,m1:m2,n1:n2,iglobal_gg:iglobal_gg+2) = 0.
!
        do j=1,1
!
          lpade=.true.
!
          select case (ipotential(j))
!
          case ('zero')           ! zero potential
            if (lroot) print*, 'initialize_gravity: zero gravity potential'
            cpot(:,j) = 0.
!
          case ('solar')          ! solar case
            if (lroot) print*, 'initialize_gravity: solar gravity potential'
            cpot(:,j) = (/ 5.088, -4.344, 61.36, 10.91, -13.93 /)
!
          case ('M5-dwarf')       ! M5 dwarf
            if (lroot) print*, 'initialize_gravity: M5 dwarf gravity potential'
            cpot(:,j) = (/ 2.3401, 0.44219, 2.5952, 1.5986, 0.20851 /)
!
          case ('M2-sgiant')       ! M super giant
            if (lroot) print*, 'initialize_gravity: M super giant gravity potential'
            cpot(:,j) = (/ 1.100, 0.660, 2.800, 1.400, 0.100 /)
!
          case ('A7-star')       ! Ap star
            if (lroot) print*, 'initialize_gravity: A star gravity potential'
            cpot(:,j) = (/ 4.080, -3.444, 15.2000, 11.2000, -12.1000 /)
!
          case ('A0-star')       ! A0 star
            if (lroot) print*, 'initialize_gravity: A0 star gravity potential'
            !cpot(:,j) = (/ 4.7446,  -1.9456,  0.6884,  4.8007, 1.79877 /)
            cpot(:,j) = (/ 4.3641,  -1.5612,  0.4841, 4.0678, 1.2548 /)
!
          case ('simple')         ! simple potential for tests
            if (lroot) print*, 'initialize_gravity: very simple gravity potential'
            cpot(:,j) =  (/ 1., 0., 0., 1., 0. /)
!
          case ('simple-2')       ! another simple potential for tests
            if (lroot) print*, 'initialize_gravity: simple gravity potential'
            cpot(:,j) =  (/ 1., 1., 0., 1., 1. /)
!
          case ('smoothed-newton')
            if (lroot) print*,'initialize_gravity: smoothed 1/r potential'
            lpade=.false.
!
          case ('no-smooth','newton','newtonian')
            if (lroot) print*,'initialize_gravity: non-smoothed newtonian gravity'
            lpade=.false.
!
          case ('varying-q')
            if (lroot) print*,'initialize_gravity: shear with Omega proto r^-q, q=',qgshear
            lpade=.false.
!
          case ('varying-q-smooth')
            if (lroot) &
                 print*,'initialize_gravity: shear with smoothed Omega proto r^-q, q=',qgshear
            lpade=.false.
!            
          case ('dark-matter-halo')
            if (lroot) &
                 print*,'initialize_gravity: arc-tangent potential generated by a dark matter halo'
            lpade=.false.
!
          case ('light-matter')
            if (lroot) &
                 print*,'initialize_gravity: potential generated by an exponential baryonic disk'
            lpade=.false.
!
            ! geodynamo
          case ('geo-kws-approx')     ! approx. 1/r potential between r=.5 and r=1
            if (lroot) print*, 'initialize_gravity: approximate 1/r potential'
            cpot(:,j) = (/ 0., 2.2679, 0., 0., 1.1697 /)
!
          case ('geo-benchmark')      ! for geodynamo benchmark runs
            if (lroot) print*, 'initialize_gravity: gravity linear in radius'
            cpot(:,j) = (/ 0., -.5, 0., 0., 0. /)
!
          case ('geo-kws')
            if (lroot) print*, 'initialize_gravity: '//&
                 'smoothed 1/r potential in spherical shell'
            if (r0_pot < epsi) print*, 'WARNING: grav_r: r0_pot is too small.'//&
                 'Can be set in grav_r namelists.'
            lpade=.false.
            ! end geodynamo
!
          case default
!
!  Catch unknown values
!
            if (lroot) print*, 'initialize_gravity: '//&
                 'No such value for ipotential: ', trim(ipotential(j))
            call stop_it("")
!
          endselect
!
          do n=n1,n2
            do m=m1,m2
!
!  rr_mn differs depending on the coordinate system
!
              call get_radial_distance(rr_sph,rr_cyl)
!
!  choose between sphere-in-box and cylindrical gravity
!
              if (lcylindrical_gravity) then
                rr_mn=rr_cyl
              else
! sphere in a box
                rr_mn=rr_sph
              endif
!
              if (lpade) then
!
                g_r = - rr_mn * poly( (/ 2*(cpot(1,j)*cpot(4,j)-cpot(2,j)), &
                     3*(cpot(1,j)*cpot(5,j)-cpot(3,j)), &
                     4*cpot(1,j)*cpot(3,j), &
                     cpot(5,j)*cpot(2,j)-cpot(3,j)*cpot(4,j), &
                     2*cpot(2,j)*cpot(3,j), &
                     cpot(3,j)**2  /), rr_mn) &
                     / poly( (/ 1., 0., cpot(4,j), cpot(5,j), &
                     cpot(3,j) /), rr_mn)**2
              else
                if (ipotential(j) == 'no-smooth'.or.&
                    ipotential(j) == 'newton'.or.&
                    ipotential(j) == 'newtonian') then
                  g_r=-g0/rr_mn**2
                elseif (ipotential(j) == 'varying-q') then
                  g_r=-g0/rr_mn**(2*qgshear-1)
                elseif (ipotential(j) == 'varying-q-smooth') then
                  g_r=-g0*rr_mn/(rr_mn**2+r0_pot**2)**qgshear
                elseif (ipotential(j) == 'dark-matter-halo') then
                  g_r=-g01(j)*(1-rpot(j)/rr_mn*atan2(rr_mn,rpot(j)))/rr_mn
                elseif (ipotential(j) == 'light-matter') then
                  !approximation of the bessel functions potential
                  !of Freeman 1970. Reference: Persic et al 1996
                  g_r=-4.134e-4*g01(j)*(.5*rr_mn/rpot(j))**1.22/&
                       ((.5*rr_mn/rpot(j))**2+1.502)**1.43/rr_mn
                else
!
!  smoothed 1/r potential in a spherical shell
!  r0_pot is the smoothing radius, and n_pot the smoothing exponent
!
                  !g_r=-g0*rr_mn**(n_pot-1) &
                  !     *(rr_mn**n_pot+r0_pot**n_pot)**(-1./n_pot-1.)
                  pot = -g0*(1. + (r1_pot1*rr_mn)**n_pot1)**(1.0/n_pot1) &
                           /(rr_mn**n_pot + r0_pot**n_pot)**(1.0/n_pot)
                  g_r = pot*(rr_mn**(n_pot-1)/(rr_mn**n_pot+r0_pot**n_pot) &
                           - r1_pot1*(r1_pot1*rr_mn)**(n_pot1-1) &
                                    /(1 + (r1_pot1*rr_mn)**n_pot1))
                endif
              endif
!
              call get_gravity_field(g_r,gg_mn,rr_mn)
              f(l1:l2,m,n,iglobal_gg:iglobal_gg+2) = &
                  f(l1:l2,m,n,iglobal_gg:iglobal_gg+2) + gg_mn
!
            enddo
          enddo
!
        enddo
      endif
!
!  Add a simple vertical gravity as well.
!
      select case (gravz_profile)
!
      case ('zero')
!
      case ('const')
        if (lroot) print*,'initialize_gravity: constant gravz=', gravz
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iglobal_gg+2)=f(l1:l2,m,n,iglobal_gg+2)+gravz
        enddo; enddo
!
      case ('linear')
        if (lroot) print*,'initialize_gravity: linear z-grav, nu=', nu_epicycle
        do n=n1,n2; do m=m1,m2
          f(l1:l2,m,n,iglobal_gg+2)=f(l1:l2,m,n,iglobal_gg+2)-nu_epicycle**2*z(n)
        enddo; enddo
!
      case default
        if (lroot) print*, &
            'initialize_gravity: unknown gravz_profile ', gravz_profile
        call fatal_error('initialize_gravity','chosen gravz_profile not valid')
!
      endselect
!
    endsubroutine initialize_gravity
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
    subroutine init_gg(f)
!
!  initialise gravity; called from start.f90
!  10-jan-02/wolf: coded
!  24-nov-02/tony: renamed from init_grav for consistancy (i.e. init_[variable name])
!
      real, dimension (mx,my,mz,mfarray) :: f
!
! Not doing anything (this might change if we decide to save gg to a file)
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
      if (lcorotational_frame) lpenc_requested(i_uu)=.true.
      if (idiag_torque/=0)     lpenc_diagnos(i_rho) =.true.
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
      if ((lhydro.and.lgravity_gas).or.&
          (lneutralvelocity.and.lgravity_neutrals).or.&
          (ldustvelocity.and.lgravity_dust)) &
          lpencil_in(i_gg)=.true.
!
      call keep_compiler_quiet(lpencil_in)
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
      use FArrayManager, only: farray_use_global
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3):: ggp
      type (pencil_case) :: p
      integer, pointer :: iglobal_gg
!
      intent(in) :: f
!
      call farray_use_global('global_gg',iglobal_gg)
!
      if (lpencil(i_gg)) then
        call farray_use_global('global_gg',iglobal_gg)
        p%gg = f(l1:l2,m,n,iglobal_gg:iglobal_gg+2)
!
!  If there is a secondary body whose mass is changing in time (ramped-up), then 
!  this gravity is not added to the global array. It is re-calculated instead. 
!
        if (g1/=0 .and. (lramp_mass.or.lsecondary_wait)) then
          call secondary_body_gravity(ggp)
          p%gg = p%gg + ggp
        endif
      endif
!
    endsubroutine calc_pencils_gravity
!***********************************************************************
    subroutine duu_dt_grav(f,df,p)
!
!  add duu/dt according to gravity
!
!  10-jan-02/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      integer :: k
!
! if statement for testing purposes
!
      if (lgravity_gas) then
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%gg
      endif
!
      if (lneutralvelocity.and.lgravity_neutrals) then
        df(l1:l2,m,n,iunx:iunz) = df(l1:l2,m,n,iunx:iunz) + p%gg
      endif
!
      if (ldustvelocity.and.lgravity_dust) then
        do k=1,ndustspec
          df(l1:l2,m,n,iudx(k):iudz(k)) = df(l1:l2,m,n,iudx(k):iudz(k)) + p%gg
        enddo
      endif
!
!  Indirect term for binary systems with origin at the primary.
!
      if (lcorotational_frame) call indirect_plus_inertial_terms(df,p)
!
      if (ldiagnos) then
         if (idiag_torque/=0) call calc_torque(p)
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine gravity_after_boundary(f)
!
!  For actions outside mn-loop. 
!
!  9-jun-15/MR: coded
!
      real, dimension(mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)

    endsubroutine gravity_after_boundary
!***********************************************************************
    subroutine indirect_plus_inertial_terms(df,p)
!
!  Add the indirect terms from the motion of the primary in a non-inertial frame,
!  plus the Coriolis and centrifugal terms from the motion of the secondary.
!
!  The indirect potential is Phi = GMp/rp**3 r dot rp
!
!  Mp is planet mass and rp planet position. For rp=(1,0,0), the potential is
!
!    Phi = gp * x = gp * r*cos(phi)
!
!  with gp = GMp/rp**2. So the resulting acceleration is, in cylindrical: 
!
!    grad = -       dPhi/drad = - gp * cos(phi)
!    gphi = - 1/r * dPhi/dphi =   gp * sin(phi)
!    gzed = -       dPhi/dzed =   0.
!
!  In spherical coordinates, x = r*sin(tht)*cos(phi), so 
!
!    grad = -                  dPhi/drad = - gp * sin(tht)*cos(phi)
!    gtht = - 1/ r           * dPhi/dtht = - gp * cos(tht)*cos(phi)
!    gphi = - 1/[r*sin(tht)] * dPhi/dphi =   gp *          sin(phi)  
!
!  20-jul-14/wlad: coded
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: rrcyl_mn
      real :: c2,s2,g2
      type (pencil_case) :: p
!
      if (lramp_mass.and.(t<=t_ramp_mass)) then 
!
!  Ramp up g1 for the first 5 orbits, to prevent to big an initial impulse
!
        g2 = g1/rp1**2 * t*t1_ramp_mass
      else
        g2 = g1/rp1**2
      endif
!
!  Do not allow secondary before t_start_secondary
!
      if (lsecondary_wait .and. t <= t_start_secondary) then
        g2 = 0. 
      endif
!
      if (lcylindrical_coords) then
!
!  Indirect terms
!
        if (lindirect_terms) then
          df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) - g2*cos(y(m))
          df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + g2*sin(y(m))
        endif
!
!  Coriolis force
!
        if (lcoriolis_force_gravity) then
          c2 = 2*Omega_corot
          df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + c2*p%uu(:,2)
          df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) - c2*p%uu(:,1)
        endif
!
!  Centrifugal force
!
        if (lcentrifugal_force_gravity) &
             df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + x(l1:l2)*Omega_corot**2
!
      else if (lspherical_coords) then
!
!  Indirect terms
!
        if (lindirect_terms) then
          df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) - g2*sinth(m)*cosph(n)
          df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) - g2*costh(m)*cosph(n)
          df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + g2*         sinph(n)
        endif
!
!  Coriolis force
!
        if (lcoriolis_force_gravity) then
          c2 = 2*Omega_corot*costh(m)
          s2 = 2*Omega_corot*sinth(m)
!
          df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + s2*p%uu(:,3)
          df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + c2*p%uu(:,3)
          df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) - c2*p%uu(:,2) - s2*p%uu(:,1)
        endif
!
!  Centrifugal force
!
        if (lcentrifugal_force_gravity) then
          rrcyl_mn=x(l1:l2)*sinth(m)
          df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + rrcyl_mn*sinth(m)*Omega_corot**2
          df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + rrcyl_mn*costh(m)*Omega_corot**2
        endif
!
      endif
!
    endsubroutine indirect_plus_inertial_terms
!***********************************************************************
    subroutine potential_global(pot,pot0)
!
!  gravity potential; version called by init_hydro, which operates on
!  full global coordinate arrays
!
!  16-jul-02/wolf: coded
!
      use Sub,     only: poly
      use Mpicomm, only: stop_it
!
      real, dimension (:,:,:) :: pot
      real, optional :: pot0           ! potential at r=0
      real, dimension (size(pot,1),size(pot,2),size(pot,3)) :: rr
      integer :: j
!
      intent(out) :: pot,pot0
!
      if (lcylindrical_gravity) &
           call stop_it("gravity_r: potential global not implemented "//&
           "for cylindrical gravity approximation")
!
!  remove this if you are sure rr is already calculated elsewhere
!
      if     (coord_system=='cartesian') then
        do n=n1,n2; do m=m1,m2
          rr(l1:l2,m,n)=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
        enddo; enddo
      elseif (coord_system=='cylindric') then
        do n=n1,n2; do m=m1,m2
          rr(l1:l2,m,n)=sqrt(x(l1:l2)**2+z(n)**2)
        enddo; enddo
      elseif (coord_system=='spherical') then
        do n=n1,n2; do m=m1,m2
          rr(l1:l2,m,n)=x(l1:l2)
        enddo; enddo
      endif
!
      pot=0.
      if (present(pot0)) pot0=0.
      do j=1,ninit
        select case (ipotential(j))
!
        case ('geo-kws','smoothed-newton')
          !pot = pot -g0*(rr**n_pot+r0_pot**n_pot)**(-1.0/n_pot)
          pot = pot - g0*(1. + (r1_pot1*rr)**n_pot1)**(1.0/n_pot1) &
                        /(rr**n_pot + r0_pot**n_pot)**(1.0/n_pot)
          if (present(pot0)) pot0=pot0-g0/r0_pot
!
        case ('no-smooth','newton','newtonian')
          pot = pot -g0/rr
!
        case default
          pot = pot - poly((/cpot(1,j), 0., cpot(2,j), cpot(3,j)/), rr) &
               / poly((/1., 0., cpot(4,j), cpot(5,j), cpot(3,j)/), rr)
          if (present(pot0)) pot0 = pot0-cpot(1,j)
!
        endselect
      enddo
!
    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn,pot,pot0,grav,rmn)
!
!  Gravity potential along one pencil
!
!  21-jan-02/wolf: coded
!  29-aug-07/wlad: allow arrays of arbitrary dimension
!
      use Sub,     only: poly
      use Mpicomm, only: stop_it
!
      real, dimension (:) :: pot
      real, dimension (size(pot)) :: rad
      real, optional :: ymn,zmn,pot0
      real, optional, dimension (size(pot)) :: xmn,rmn
      real, optional, dimension (size(pot),3) :: grav
      integer :: j
!
      intent(in)  :: xmn,ymn,zmn,rmn
      intent(out) :: pot,pot0,grav
!
      if (present(rmn)) then
        rad = rmn
      else
        if (present(xmn) .and. present(ymn) .and. present(zmn)) then
!
          if (.not.lcartesian_coords) &
              call stop_it("gravity_r: potential_penc with xmn,ymn,zmn is "//&
              "not yet implemented for non-cartesian coordinates. Fix the call "//&
              "to use radial distance instead")
!
          rad = sqrt(xmn**2+ymn**2+zmn**2)
        else
          call stop_it("POTENTIAL_PENC: Need to specify either x,y,z or r.")
        endif
      endif
!
      pot=0.
      if (present(pot0)) pot0=0.
      do j=1,ninit
        select case (ipotential(j))
!
        case ('geo-kws','smoothed-newton')
          !pot=pot-g0*(rmn**n_pot+r0_pot**n_pot)**(-1.0/n_pot)
          pot = pot - g0*(1. + (r1_pot1*rmn)**n_pot1)**(1.0/n_pot1) &
                        /(rmn**n_pot + r0_pot**n_pot)**(1.0/n_pot)
          if (present(pot0)) pot0=pot0-g0/r0_pot
!
        case ('no-smooth','newton','newtonian')
          pot=pot-g0/rmn
!
        case default
          pot = pot - poly((/cpot(1,j), 0., cpot(2,j), cpot(3,j)/), rad) &
               / poly((/1., 0., cpot(4,j), cpot(5,j), cpot(3,j)/), rad)
          if (present(pot0)) pot0=pot0-cpot(1,j)
!
        endselect
      enddo
!
      if (present(grav)) then
        call not_implemented("potential_penc", "optional argument grav")
        call keep_compiler_quiet(grav)
      endif
!
    endsubroutine potential_penc
!***********************************************************************
    subroutine potential_point(x,y,z,r, pot,pot0, grav)
!
!  Gravity potential in one point
!
!  20-dec-03/wolf: coded
!
      use Sub,     only: poly
      use Mpicomm, only: stop_it
!
      real :: pot,rad
      real, optional :: x,y,z,r
      real, optional :: pot0,grav
      integer :: j
!
      intent(in)  :: x,y,z,r
      intent(out) :: pot,pot0,grav
!
      if (present(r)) then
        rad = r
      else
        if (present(x) .and. present(y) .and. present(z)) then
!
          if (.not.lcartesian_coords) &
               call stop_it("gravity_r: potential_point with x,y,z is "//&
               "not yet implemented for non-cartesiand coordinates. Fix the call "//&
               "to  use radial distance instead")
!
          rad = sqrt(x**2+y**2+z**2)
        else
          call stop_it("Need to specify either x,y,z or r in potential_point()")
        endif
      endif
!
      pot=0.
      if (present(pot0)) pot0=0.
      do j=1,ninit
        select case (ipotential(j))
!
        case ('geo-kws','smoothed-newton')
          !pot=pot-g0*(rad**n_pot+r0_pot**n_pot)**(-1.0/n_pot)
          pot = pot - g0*(1. + (r1_pot1*rad)**n_pot1)**(1.0/n_pot1) &
                        /(rad**n_pot + r0_pot**n_pot)**(1.0/n_pot)
          if (present(pot0)) pot0=pot0-g0/r0_pot
!
        case ('no-smooth','newton','newtonian')
          pot=pot-g0/r
!
        case default
          pot = pot- poly((/cpot(1,j), 0., cpot(2,j), cpot(3,j)/), rad) &
               / poly((/1., 0., cpot(4,j), cpot(5,j), cpot(3,j)/), rad)
          if (present(pot0)) pot0=pot0-cpot(1,j)
!
        endselect
      enddo
!
      if (present(grav)) then
        call not_implemented("potential_point", "optional argument grav")
        call keep_compiler_quiet(grav)
      endif
!
    endsubroutine potential_point
!***********************************************************************
    subroutine acceleration_penc(gg)
!
!  Calculates gravitational acceleration on a pencil
!
!  21-apr-07/tobi: adapted from potential_penc
!
      use Messages, only: fatal_error
!
      real, dimension (:,:), intent (out) :: gg
!
!  Calculate acceleration from master pencils defined in initialize_gravity
!
      call fatal_error("acceleration_penc","Not implemented")
      call keep_compiler_quiet(gg)
!
    endsubroutine acceleration_penc
!***********************************************************************
    subroutine acceleration_penc_1D(g_r)
!
!  Gravitational acceleration along one pencil
!
!  Analogous to potential, but for the radial acceleration.
!   useful for coding initial condition with centrifugal balance
!
!  21-aug-07/wlad: coded
!
      use Mpicomm,only: stop_it
      use Sub,    only: get_radial_distance
!
      real, dimension (:) :: g_r
      real, dimension(size(g_r)) :: rr_mn,rr_sph,rr_cyl,pot
      integer :: j
!
      intent(out) :: g_r
!
      call get_radial_distance(rr_sph,rr_cyl)
      if (lcylindrical_gravity) then
        rr_mn=rr_cyl
      else
        rr_mn=rr_sph
      endif
!
      g_r=0.
      do j=1,ninit
        select case (ipotential(j))
        case ('no-smooth','newton','newtonian')
          g_r=g_r -g0/rr_mn**2
!
        case ('smoothed-newton')
          !g_r=g_r -g0*rr_mn**(n_pot-1) &
          !     *(rr_mn**n_pot+r0_pot**n_pot)**(-1./n_pot-1.)
          pot = -g0*(1. + (r1_pot1*rr_mn)**n_pot1)**(1.0/n_pot1) &
                   /(rr_mn**n_pot + r0_pot**n_pot)**(1.0/n_pot)
          g_r = pot*(rr_mn**(n_pot-1)/(rr_mn**n_pot + r0_pot**n_pot) &
                   - r1_pot1*(r1_pot1*rr_mn)**(n_pot1-1) &
                            /(1 + (r1_pot1*rr_mn)**n_pot1))
!
        case ('varying-q')
          g_r=g_r -g0/rr_mn**(2*qgshear-1)
!
        case ('varying-q-smooth')
          g_r=g_r -g0*rr_mn/(rr_mn**2+r0_pot**2)**qgshear
!
        case ('dark-matter-halo')
          g_r=g_r -g01(j)*(1-rpot(j)/rr_mn*atan2(rr_mn,rpot(j)))/rr_mn
!
        case ('light-matter')
          g_r=g_r-4.134e-4*g01(j)*(.5*rr_mn/rpot(j))**1.22/&
               ((.5*rr_mn/rpot(j))**2+1.502)**1.43/rr_mn
!
        case ('zero')
          g_r=g_r
!
        case default
          if (lroot) print*, 'acceleration: '//&
               'No such value for ipotential: ', trim(ipotential(j))
          call stop_it("")
!
        endselect
      enddo
!
    endsubroutine acceleration_penc_1D
!***********************************************************************
    subroutine acceleration_point(x,y,z,r,g_r)
!
!  Gravitational acceleration in one point
!
!  Analogous to potential, but for the radial acceleration.
!   useful for coding initial condition with centrifugal balance
!
!  18-nov-08/wlad: coded
!
      use Mpicomm,only: stop_it
!
      real :: g_r,rad,pot
      real, optional :: x,y,z,r
      integer :: j
!
      intent(in)  :: x,y,z,r
      intent(out) :: g_r
!
      if (present(r)) then
        rad = r
      else
        if (present(x) .and. present(y) .and. present(z)) then
!
          if (.not.lcartesian_coords) &
              call stop_it("gravity_r: acceleration_point with x,y,z is "//&
              "not yet implemented for non-cartesiand coordinates. Fix  "//&
              "the call to  use radial distance instead")
!
          rad = sqrt(x**2+y**2+z**2)
        else
          call stop_it("Need to specify either x,y,z or r in acceleration_point()")
        endif
      endif
!
      g_r=0.
      do j=1,ninit
        select case (ipotential(j))
        case ('no-smooth','newton','newtonian')
          g_r=g_r -g0/rad**2
!
        case ('smoothed-newton')
          !g_r=g_r -g0*rad**(n_pot-1) &
          !     *(rad**n_pot+r0_pot**n_pot)**(-1./n_pot-1.)
          pot = -g0*(1. + (r1_pot1*rad)**n_pot1)**(1.0/n_pot1) &
                   /(rad**n_pot + r0_pot**n_pot)**(1.0/n_pot)
          g_r = pot*(rad**(n_pot-1)/(rad**n_pot + r0_pot**n_pot) &
                   - r1_pot1*(r1_pot1*rad)**(n_pot1-1) &
                            /(1 + (r1_pot1*rad)**n_pot1))
!
        case ('varying-q')
          g_r=g_r -g0/rad**(2*qgshear-1)
!
        case ('varying-q-smooth')
          g_r=g_r -g0*rad/(rad**2+r0_pot**2)**qgshear
!
        case ('dark-matter-halo')
          g_r=g_r -g01(j)*(1-rpot(j)/rad*atan2(rad,rpot(j)))/rad
!
        case ('light-matter')
          g_r=g_r-4.134e-4*g01(j)*(.5*rad/rpot(j))**1.22/&
               ((.5*rad/rpot(j))**2+1.502)**1.43/rad
!
        case ('zero')
          g_r=g_r
!
        case default
          if (lroot) print*, 'acceleration: '//&
               'No such value for ipotential: ', trim(ipotential(j))
          call stop_it("")
!
        endselect
      enddo
!
    endsubroutine acceleration_point
!***********************************************************************
    subroutine get_gravity_field(gr,gg_mn,rr_mn)
!
!  Calculate gravity field for different coordinate systems
!
!  15-mar-07/wlad: coded
!
      real, dimension(nx),intent(in) :: gr,rr_mn
      real, dimension(nx,3),intent(out) :: gg_mn
      real, dimension(nx,3) :: ggp
!
      if (coord_system=='cartesian') then
        gg_mn(:,1) = x(l1:l2)/rr_mn*gr
        gg_mn(:,2) = y(  m  )/rr_mn*gr
        gg_mn(:,3) = z(  n  )/rr_mn*gr
        if (lcylindrical_gravity) gg_mn(:,3)=0.
      elseif (coord_system=='cylindric') then
        gg_mn(:,1) = x(l1:l2)/rr_mn*gr
        gg_mn(:,2) = 0.
        gg_mn(:,3) = z(  n  )/rr_mn*gr
        if (lcylindrical_gravity) gg_mn(:,3)=0.
      elseif (coord_system=='spherical') then
        gg_mn(:,2)=0.
        gg_mn(:,3)=0.
        if (lcylindrical_gravity) then
          gg_mn(:,1) = gr*sinth(m)
          gg_mn(:,2) = gr*costh(m)
        else
          gg_mn(:,1) = gr
        endif
      endif
!
!  Add the gravity of a stationary secondary body (i.e., following reference frame)
!
      if (g1/=0 .and. (.not.(lramp_mass.or.lsecondary_wait))) then 
!
!  In this case, the mass is constant in time. Add to the global array.
!
        call secondary_body_gravity(ggp)
        gg_mn = gg_mn + ggp 
      endif
!
    endsubroutine get_gravity_field
!***********************************************************************
    subroutine secondary_body_gravity(ggp)
!
!  Add the gravity of a body fixed at position (rp1,0,0).
!
!  20-jul-14/wlad: coded
!  03-may-16/wlad: added boley smoothing      
!
      real, dimension(nx,3), intent(out) :: ggp
      real, dimension(nx) :: rr2_pm,gp
      real :: rhill,rhill1
      integer :: i
!
      if (lcylindrical_coords) then
        if (lcylindrical_gravity) then
          rr2_pm = x(l1:l2)**2 + rp1**2 -2*x(l1:l2)*rp1*cos(y(m))
        else
          rr2_pm = x(l1:l2)**2 + rp1**2 -2*x(l1:l2)*rp1*cos(y(m)) + z(n)**2
        endif
      elseif (lspherical_coords) then
        rr2_pm = x(l1:l2)**2 + rp1**2 - 2*x(l1:l2)*rp1*sinth(m)*cosph(n)
      else
        call fatal_error("secondary_body_gravity","not coded for Cartesian")
      endif
!
!  Select the potential smoothing for the secondary. 
!      
      select case (ipotential_secondary)
!
      case ('plummer')
        gp = -g1*(rr2_pm+rp1_smooth**2)**(-1.5)
!           
     case ('boley')
!
!  Correct potential outside a sphere of radius rsmooth. Default to Hill sphere.
!
        do i=1,nx
          if (rr2_pm(i) .gt. rp1_smooth**2) then
            gp(i) = -g1*rr2_pm(i)**(-1.5)
          else
            gp(i) =  g1*(3*sqrt(rr2_pm(i))*rp1_smooth1 - 4)*rp1_smooth1**3
          endif
        enddo
!           
      case default
!
!  Catch unknown values
!           
        if (lroot) print*, 'secondary_body_gravity: '//&
             'No such value for ipotential: ', trim(ipotential_secondary)
        call fatal_error("","")
!
      endselect
!
      if (lramp_mass     .and.(t<=t_ramp_mass)      ) gp = gp * t*t1_ramp_mass
      if (lsecondary_wait.and.(t<=t_start_secondary)) gp = 0.
!       
!  Set the acceleration
!
      if (lcylindrical_coords) then 
        ggp(:,1) =  gp * (x(l1:l2)-rp1*cos(y(m)))
        ggp(:,2) =  gp *           rp1*sin(y(m))
        ggp(:,3) =  gp *  z(  n  )
        if (lcylindrical_gravity) ggp(:,3)=0.
      else if (lspherical_coords) then
        ggp(:,1) =  gp * (x(l1:l2) - rp1*sinth(m)*cosph(n))
        ggp(:,2) = -gp *             rp1*costh(m)*cosph(n)
        ggp(:,3) =  gp *             rp1*         sinph(n)
      endif
!
    endsubroutine secondary_body_gravity
!***********************************************************************
    subroutine calc_torque(p)
!
!  Output torque diagnostic for nbody particle ks
!
!  24-jan-20/wlad : coded
!
      use Diagnostics
!
      type (pencil_case) :: p
      real, dimension(nx) :: torque
      real, dimension(nx) :: rr2,rpre
!
      if (lcartesian_coords) then
        call fatal_error("calc_torque","not coded for Cartesian")
      elseif (lcylindrical_coords) then
        rpre = x(l1:l2)*rp1*sin(y(m))
        rr2  = x(l1:l2)**2 + rp1**2 -2*x(l1:l2)*rp1*cos(y(m)) + z(n)**2
      elseif (lspherical_coords) then
        rpre = x(l1:l2)*rp1*sinth(m)*sinph(n)
        rr2  = x(l1:l2)**2 + rp1**2 - 2*x(l1:l2)*rp1*sinth(m)*cosph(n)
      else
        call fatal_error("calc_torque",&
             "the world is flat and we should never gotten here")
      endif
!
!  Define separate torques for gas and dust/particles
!
      torque = g1*p%rho*rpre*(rr2 + rp1_smooth**2)**(-1.5)
!
      call integrate_mn_name(torque,idiag_torque)
!
    endsubroutine calc_torque
!***********************************************************************
    subroutine rprint_gravity(lreset,lwrite)
!
!  reads and registers print parameters relevant for gravity advance
!  dummy routine
!
!  26-apr-03/axel: coded
!
      use FArrayManager, only: farray_index_append
      use Diagnostics
!
      logical :: lreset,lwr
      logical, optional :: lwrite
      integer :: iname
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lreset) then
        idiag_torque=0
      endif

      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'torque',idiag_torque)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!  idl needs this even if everything is zero
!
      if (lwr) then
        call farray_index_append('igg',igg)
        call farray_index_append('igx',igx)
        call farray_index_append('igy',igy)
        call farray_index_append('igz',igz)
      endif
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_gravity
!***********************************************************************
    subroutine compute_gravity_star(f, wheat, luminosity, star_cte)
!
!  Compute the global gravity field for a star-in-a-box with a central heating
!  05-jan-2010/dintrans: coded
!
    use Sub, only: get_radial_distance, erfunc
!
    real, dimension (mx,my,mz,mfarray) :: f
    real :: wheat, luminosity, star_cte
    real, dimension (nx) :: g_r, rr_mn, rr_sph, rr_cyl, u_mn, lumi_mn
    real, dimension (nx,3) :: gg_mn=0.0
!
    do n=n1,n2
    do m=m1,m2
      call get_radial_distance(rr_sph, rr_cyl)
      rr_mn=rr_sph
      u_mn=rr_mn/sqrt(2.)/wheat
      if (nzgrid==1) then
        lumi_mn=luminosity*(1.-exp(-u_mn**2))
        g_r=-lumi_mn/(2.*pi*rr_mn)*star_cte
      else
        lumi_mn=luminosity*(erfunc(u_mn)-2.*u_mn/sqrt(pi)*exp(-u_mn**2))
        g_r=-lumi_mn/(4.*pi*rr_mn**2)*star_cte
      endif
      call get_gravity_field(g_r, gg_mn, rr_mn)
      f(l1:l2,m,n,iglobal_gg:iglobal_gg+2) = gg_mn
    enddo
    enddo
!
    endsubroutine compute_gravity_star
!***********************************************************************
    subroutine get_xgravity(xgrav)
!
!  Dummy routine for getting the gravity profile into the
!  intial conditions
!
!  04-oct-10/bing: coded
!
      real, dimension(nx) :: xgrav
!
      call keep_compiler_quiet(xgrav)
!
    endsubroutine get_xgravity
!***********************************************************************
    subroutine set_consistent_gravity(ginput,gtype,gprofile,lsuccess)
!
! Dummy routine
!
      real :: ginput
      character (len=labellen) :: gtype,gprofile
      logical :: lsuccess
!
      call keep_compiler_quiet(ginput)
      call keep_compiler_quiet(gtype,gprofile)
!
! This routine should never be called in the way it is written now.
!
      lsuccess=.false.
!
! gravity parameters set consistently.
!
    endsubroutine set_consistent_gravity
!***********************************************************************
    logical function is_constant_zgrav()
!
!  14-apr-15/MR: coded
!
      is_constant_zgrav=.false.

    endfunction is_constant_zgrav
!***********************************************************************
endmodule Gravity
