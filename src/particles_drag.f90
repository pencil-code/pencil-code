! $Id$
!
!  This module integrates drag forces between particles and gas.
!
!  Reference:
!    Yang, C.-C., & Johansen, A. 2016, ApJS, in press (arXiv:1603.08523)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_drag = .true.
!
!***************************************************************
module Particles_drag
!
  use Cdata
  use Cparam
  use Messages
  use Particles_cdata
  use Particles_map
!
  implicit none
!
  include 'particles_drag.h'
!
!  Initialization parameters
!
  logical :: ldrag_equilibrium_global = .false.  ! Set global/local drag equilibrium when lset_drag_equilibrium = .true.
  logical :: ldrag_on_gas = .false.              ! Turn on/off drag on gas
  logical :: ldrag_on_par = .false.              ! Turn on/off drag on particles
  logical :: lset_drag_equilibrium = .false.     ! Initialize the velocities of the particles and the gas under equilibrium
  real :: gx_gas = 0.0        ! Background acceleration of gas in x direction
  real :: gz_par_coeff = 0.0  ! Coefficient for background acceleration of a particle in z direction
  real :: taus = 0.0          ! Dimensionless stopping time: Omega * tdrag
  real :: tdrag = 0.0         ! Drag timescale
!
  namelist /particles_drag_init_pars/ &
      ldrag_equilibrium_global, ldrag_on_gas, ldrag_on_par, lset_drag_equilibrium, gx_gas, gz_par_coeff, taus, tdrag
!
!  Runtime parameters
!
  logical :: ldrag_pm_back_reaction = .true.  ! Couple the gas by particle-meshing back reaction from particles
!
  namelist /particles_drag_run_pars/ ldrag_on_gas, ldrag_on_par, ldrag_pm_back_reaction, gx_gas, gz_par_coeff, taus, tdrag
!
!  Module variables
!
  real :: dv_gas = 0.0
  real :: epicycle_freq = 0.0
  real :: epicycle_ratio = 0.0
  real :: gzpc = 0.0
  real :: oneplustaus2inv = 0.0
  real :: taus2 = 0.0
  real :: twominusqtaus = 0.0
  real :: twotaus = 0.0
  real :: twoomega1 = 0.0
!
  contains
!***********************************************************************
    subroutine register_particles_drag()
!
!  Register this module.
!
!  14-dec-14/ccyang: coded.
!
      use SharedVariables, only: put_shared_variable
!
      if (lroot) call svn_id("$Id$")
!
      call put_shared_variable("ldrag_on_gas", ldrag_on_gas, caller="register_particles_drag")
      call put_shared_variable("taus", taus, caller="register_particles_drag")
      call put_shared_variable("gx_gas", gx_gas, caller="register_particles_drag")
!
    endsubroutine register_particles_drag
!***********************************************************************
    subroutine initialize_particles_drag()
!
!  Perform any post-parameter-read initialization, i.e., calculate
!  derived parameters.
!
!  11-dec-15/ccyang: coded.
!
!  Sanity checks.
!
      if (coord_system /= 'cartesian' .and. Omega /= 0.0) &
          call fatal_error('initialize_particles_drag', 'non-inertial curvilinear system is not implemented. ')
      if (ldrag_on_gas .and. iuu == 0) call fatal_error('initialize_particles_drag', 'drag on gas requires uu. ')
!
!  Check the drag timescale.
!
      drag: if (ldrag_on_par) then
        settdrag: if (tdrag == 0.0) then
          if (taus == 0.0) call fatal_error('initialize_particles_drag', 'tdrag = 0')
          if (Omega == 0.0) call fatal_error('initialize_particles_drag', 'Omega = 0')
          tdrag = taus / Omega
        elseif (taus == 0.0) then settdrag
          taus = Omega * tdrag
        endif settdrag
      endif drag
      taus2 = 2.0 * (2.0 - qshear) * taus**2
      twotaus = 2.0 * taus
      twominusqtaus = (2.0 - qshear) * taus
      oneplustaus2inv = 1.0 / (1.0 + taus2)
      if (lroot) print *, 'initialize_particles_drag: tdrag = ', tdrag
!
!  Find the frequency and aspect ratio of the epicycles.
!
      epicycle_freq = sqrt(2.0 * (2.0 - qshear)) * Omega
      epicycle_ratio = sqrt(2.0 / (2.0 - qshear))
!
!  Get the velocity reduction of the gas.
!
      dv: if (gx_gas /= 0.0) then
        if (Omega == 0.0) call fatal_error('initialize_particles_drag', 'Omega = 0')
        if (.not. ldrag_on_gas .and. .not. ldrag_on_par) call fatal_error('initialize_particles_drag', 'gx_gas has no effect. ')
        twoomega1 = 0.5 / Omega
        dv_gas = twoomega1 * gx_gas
      endif dv
!
!  Initialization for particle_zaccel().
!
      gzpc = gz_par_coeff**2
!
!  Calculate and overwrite mp_swarm and rhop_swarm.
!
      mprhop: if (eps_dtog > 0.0) then
        mp_swarm = find_mp_swarm(eps_dtog)
        rhop_swarm = mp_swarm / product((/dx,dy,dz/), lactive_dimension)
        info: if (lroot) then
          print *, 'initialize_particles_drag: reset mp_swarm = ', mp_swarm
          print *, 'initialize_particles_drag: reset rhop_swarm = ', rhop_swarm
        endif info
      endif mprhop
!
    endsubroutine initialize_particles_drag
!***********************************************************************
    subroutine init_particles_drag(f, fp)
!
!  Set some initial conditions for the gas and the particles, if any.
!
!  08-may-16/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
!
      type(pic), dimension(nx,ny,nz) :: cell
      type(particle), dimension(:), pointer :: ghost
      integer, dimension(:), pointer :: sendlist
      integer, dimension(0:nproc_comm) :: ngp_send, ngp_recv
      integer :: npsend
!
      real :: ux, uy, vx, vy
!
!  Initialize the gas and the particle velocities with drag equilibrium, if requested.
!
      init: if (lset_drag_equilibrium) then
        global: if (ldrag_equilibrium_global) then
          if (lroot) print *, 'init_particles_drag: Set global drag equilibrium. '
          call get_nsh_solution(dv_gas, eps_dtog, ux, uy, vx, vy)
          f(l1:l2,m1:m2,n1:n2,iux) = f(l1:l2,m1:m2,n1:n2,iux) + ux
          f(l1:l2,m1:m2,n1:n2,iuy) = f(l1:l2,m1:m2,n1:n2,iuy) + uy
          fp(1:npar_loc,ivpx) = fp(1:npar_loc,ivpx) + vx
          fp(1:npar_loc,ivpy) = fp(1:npar_loc,ivpy) + vy
        else global
          if (lroot) print *, 'init_particles_drag: Set local drag equilibrium. '
          if (lparticles_blocks) call fatal_error('init_particles_drag', 'particles_blocks version is not implemented yet. ')
          nullify(ghost)
          call distribute_particles(f, fp, npsend, sendlist, ghost, cell, ngp_send, ngp_recv)
          call set_drag_equilibrium(cell)
          call collect_particles(f, fp, npsend, sendlist, ghost, cell, ngp_send, ngp_recv, &
                                 ldrag_on_par, ldrag_on_gas .or. Omega /= 0.0, .false.)
        endif global
      endif init
!
    endsubroutine init_particles_drag
!***********************************************************************
    subroutine read_particles_drag_init_pars(iostat)
!
!  Read initialization parameters from namelist particles_drag_init_pars.
!
!  27-aug-15/ccyang: coded.
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_drag_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_drag_init_pars
!***********************************************************************
    subroutine write_particles_drag_init_pars(unit)
!
!  Write initialization parameters to namelist particles_drag_init_pars.
!
!  14-dec-15/ccyang: coded.
!
      integer, intent(in) :: unit
!
      integer :: stat
!
      write(unit, NML=particles_drag_init_pars, IOSTAT=stat)
      if (stat /= 0) call fatal_error('write_particles_drag_init_pars', 'cannot write particles_drag_init_pars. ', force=.true.)
!
    endsubroutine write_particles_drag_init_pars
!***********************************************************************
    subroutine read_particles_drag_run_pars(iostat)
!
!  Read runtime parameters from namelist particles_drag_run_pars.
!
!  27-aug-15/ccyang: coded.
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_drag_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_drag_run_pars
!***********************************************************************
    subroutine write_particles_drag_run_pars(unit)
!
!  Write runtime parameters to namelist particles_drag_run_pars.
!
!  14-dec-15/ccyang: coded.
!
      integer, intent(in) :: unit
!
      integer :: stat
!
      write(unit, NML=particles_drag_run_pars, IOSTAT=stat)
      if (stat /= 0) call fatal_error('write_particles_drag_run_pars', 'cannot write particles_drag_run_pars. ', force=.true.)
!
    endsubroutine write_particles_drag_run_pars
!***********************************************************************
    subroutine integrate_drag(f, fp, dt)
!
!  Wrapper for the integration of the drag force between particles and
!  gas.
!
!  08-may-16/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      real, intent(in) :: dt
!
      type(pic), dimension(nx,ny,nz) :: cell
      type(particle), dimension(:), pointer :: ghost
      integer, dimension(:), pointer :: sendlist
      integer, dimension(0:nproc_comm) :: ngp_send, ngp_recv
      integer :: npsend
!
      if (lparticles_blocks) call fatal_error('integrate_drag', 'particles_blocks version is not implemented yet. ')
!
!  Distribute particles to the cells.
!
      nullify(ghost)
      call distribute_particles(f, fp, npsend, sendlist, ghost, cell, ngp_send, ngp_recv)
!
!  Locally integrate the drag force as well as other source terms.
!
      if (ldrag_on_par .and. ldrag_on_gas) then
        call drag_on_both(cell, dt)
      elseif (ldrag_on_par) then
        call drag_on_particles(cell, dt)
      elseif (ldrag_on_gas) then
        call fatal_error('integrate_drag', 'drag on gas only is not implemented. ')
      endif
!
!  Assign the properties of the particles onto the grid.
!
      call map_particles(f, cell)
!
!  Collect the change in velocities.
!
      call collect_particles(f, fp, npsend, sendlist, ghost, cell, ngp_send, ngp_recv, &
                             ldrag_on_par, ldrag_on_gas .or. Omega /= 0.0, ldrag_on_gas .and. ldrag_pm_back_reaction)
!
    endsubroutine integrate_drag
!***********************************************************************
!***********************************************************************
! LOCAL SUBROUTINES START BELOW HERE.
!***********************************************************************
!***********************************************************************
    elemental subroutine drag_on_particles(cell, dt)
!
!  Finds the change in velocity for the particles under the drag.
!
!  11-dec-15/ccyang: coded.
!
      type(pic), intent(inout) :: cell
      real, intent(in) :: dt
!
!  Find the velocity change.
!
      rotation: if (Omega /= 0.0) then
        if (cell%np > 0) call drag_particle_omega(cell%p%v(1), cell%p%v(2), cell%u(1), cell%u(2), dt, cell%p%dv(1), cell%p%dv(2))
        call gas_epicycle(cell%u(1), cell%u(2), dt, cell%du(1), cell%du(2))
      else rotation
        call drag_particle(cell%p%v(1), cell%u(1), dt, cell%p%dv(1))
        call drag_particle(cell%p%v(2), cell%u(2), dt, cell%p%dv(2))
      endif rotation
!
      if (gz_par_coeff /= 0.0) then
        call drag_particle(cell%p%v(3), cell%u(3), dt, cell%p%dv(3), particle_zaccel(cell%p%x(3)))
      else
        call drag_particle(cell%p%v(3), cell%u(3), dt, cell%p%dv(3))
      endif
!
    endsubroutine drag_on_particles
!***********************************************************************
    elemental subroutine drag_on_both(cell, dt)
!
!  Finds the change in velocity for the gas and the particles under
!  mutual drags as well as other source terms.
!
!  13-sep-15/ccyang: coded.
!
      type(pic), intent(inout) :: cell
      real, intent(in) :: dt
!
      rotation: if (Omega /= 0.0) then
        call drag_mutual_omega(cell%p%eps, cell%p%v(1), cell%p%v(2), cell%u(1), cell%u(2), dt, &
                               cell%p%dv(1), cell%p%dv(2), cell%du(1), cell%du(2))
      else rotation
        if (gx_gas /= 0.0) then
          call drag_mutual(cell%p%eps, cell%p%v(1), cell%u(1), dt, cell%p%dv(1), cell%du(1), agas=gx_gas)
        else
          call drag_mutual(cell%p%eps, cell%p%v(1), cell%u(1), dt, cell%p%dv(1), cell%du(1))
        endif
        call drag_mutual(cell%p%eps, cell%p%v(2), cell%u(2), dt, cell%p%dv(2), cell%du(2))
      endif rotation
!
      if (gz_par_coeff /= 0.0) then
        call drag_mutual(cell%p%eps, cell%p%v(3), cell%u(3), dt, cell%p%dv(3), cell%du(3), apar=particle_zaccel(cell%p%x(3)))
      else
        call drag_mutual(cell%p%eps, cell%p%v(3), cell%u(3), dt, cell%p%dv(3), cell%du(3))
      endif
!
    endsubroutine drag_on_both
!***********************************************************************
    elemental subroutine drag_particle(v, u, dt, dv, ap)
!
!  Calculates the velocity change of each particle under the gas drag.
!
!  18-jun-15/ccyang: coded.
!
!  Input Arguments:
!      v        Initial particle velocity.
!      u        Gas velocity.
!      dt       Time step.
!  Output Arguments:
!      dv       Change in particle velocity in time dt.
!  Optional Argument:
!      ap       Background acceleration on the particle.
!
      use Sub, only: one_minus_exp
!
      real, intent(in) :: v, u, dt
      real, intent(out) :: dv
      real, intent(in), optional :: ap
!
      real :: t
!
      t = dt / tdrag
      getdv: if (present(ap)) then
        if (t * t > epsilon(1.0)) then
          dv = (u - v + ap * tdrag) * (1.0 - exp(-t))
        else
          dv = ((u - v) / tdrag + ap) * (1.0 - 0.5 * t) * dt
        endif
      else getdv
        dv = (u - v) * one_minus_exp(t)
      endif getdv
!
    endsubroutine drag_particle
!***********************************************************************
    pure subroutine drag_particle_omega(vx, vy, ux, uy, dt, dvx, dvy)
!
!  Calculates the velocity change of each particle under the gas drag
!  and rotation.
!
!  05-may-18/ccyang: coded.
!
!  Input Arguments:
!      vx       x-component of the initial velocities of the particles.
!      vy       y-component of the initial velocities of the particles.
!      ux       x-component of the gas velocity.
!      uy       y-component of the gas velocity.
!      dt       Time step.
!  Output Arguments:
!      dvx      Changes in x-component of particle velocities in time dt.
!      dvy      Changes in y-component of particle velocities in time dt.
!
      real, dimension(:), intent(in) :: vx, vy
      real, intent(in) :: ux, uy, dt
      real, dimension(:), intent(out) :: dvx, dvy
!
      real :: cosot, sinot1, sinot2
      real :: x, y, vxe, vye, vx0, vy0
!
!  Find the epicycle.
!
      x = epicycle_freq * dt
      cosot = cos(x)
      y = sin(x)
      sinot1 = epicycle_ratio * y
      sinot2 = y / epicycle_ratio
!
!  Find the decay factor.
!
      y = exp(-dt / tdrag)
      x = 1.0 - y
!
!  Find the terminal velocity.
!
      vye = -oneplustaus2inv * dv_gas
      vxe = twotaus * vye
!
      vx0 = vxe + (x * ux - y * vxe) * cosot + (x * (uy + dv_gas) - y * vye) * sinot1
      vy0 = vye + (x * (uy + dv_gas) - y * vye) * cosot - (x * ux - y * vxe) * sinot2
!
!  Find the velocity changes of the particles.
!
      x = 1.0 - y * cosot
      sinot1 = y * sinot1
      sinot2 = y * sinot2
!
      dvx = vx0 - (x * vx - vy * sinot1)
      dvy = vy0 - (x * vy + vx * sinot2)
!
    endsubroutine drag_particle_omega
!***********************************************************************
    pure subroutine drag_mutual(eps, v, u, dt, dv, du, apar, agas)
!
!  Calculates the velocity change of the gas and the particles under
!  mutual drag interaction.
!
!  13-sep-15/ccyang: coded.
!
!  Input Arguments:
!      eps      Particle-to-gas density ratio of each particle.
!      v        Initial particle velocity.
!      u        Initial gas velocity.
!      dt       Time step.
!  Output Arguments:
!      dv       Change in particle velocity in time dt.
!      du       Change in gas velocity in time dt, if
!          ldrag_pm_back_reaction = .false.; (1 + sum(eps)) times the
!          change in center-of-mass velocity, otherwise.
!  Optional Arguments:
!      apar     External acceleration of each of the particles.
!      agas     External acceleartion of the gas.
!
      use Sub, only: one_minus_exp
!
      real, dimension(:), intent(in) :: eps
      real, dimension(size(eps)), intent(in) :: v
      real, intent(in) :: u, dt
      real, dimension(size(eps)), intent(out) :: dv
      real, intent(out) :: du
      real, dimension(size(eps)), intent(in), optional :: apar
      real, intent(in), optional :: agas
!
      real :: epstot, norm, ucm, du0, acm, dag
      real :: t, x, y, z
!
!  Find weights.
!
      epstot = sum(eps)
      norm = 1.0 / (1.0 + epstot)
!
!  Find exponential decays.
!
      t = dt / tdrag
      x = one_minus_exp(t)
      if (.not. ldrag_pm_back_reaction) y = one_minus_exp((1.0 + epstot) * t)
      z = exp(-t) * one_minus_exp(epstot * t) / epstot
!
!  Find center-of-mass velocity.
!
      ucm = norm * (u + sum(eps * v))
      du0 = ucm - u
!
!  Find the velocity changes due to mutual drag.
!
      dv = (ucm - v) * x - du0 * z
      if (ldrag_pm_back_reaction) then
        if (.not. present(apar) .and. .not. present(agas)) du = 0.0
      else
        du = du0 * y
      endif
!
!  Add the velocity changes due to external accelerations.
!
      accel: if (present(apar) .or. present(agas)) then
!       Find center-of-mass acceleration.
        if (present(apar)) then
          acm = sum(eps * apar)
        else
          acm = 0.0
        endif
        if (present(agas)) acm = acm + agas
        if (ldrag_pm_back_reaction) du = acm * dt
        acm = norm * tdrag * acm
!       Find the velocity changes.
        if (present(apar)) dv = dv + (apar * tdrag - acm) * x
        if (present(agas)) then
          dag = norm * (agas * tdrag - acm)
        else
          dag = -norm * acm
        endif
        du0 = acm * t
        dv = dv + (du0 + dag * (x - z))
        if (.not. ldrag_pm_back_reaction) du = du + (du0 + dag * y)
      endif accel
!
    endsubroutine drag_mutual
!***********************************************************************
    pure subroutine drag_mutual_omega(eps, vx, vy, ux, uy, dt, dvx, dvy, dux, duy, agx)
!
!  Calculates the change in the horizontal components of velocity for
!  the gas and the particles under mutual drag interaction, shear,
!  Coriolis force, and background gas pressure gradient.
!
!  12-jun-15/ccyang: coded
!
!  Input Arguments:
!      eps      Particle-to-gas density ratio of each particle.
!      vx       x-component of the initial velocity of each particle.
!      vy       y-component of the initial velocity of each particle.
!      ux       x-component of the initial gas velocity.
!      uy       y-component of the initial gas velocity.
!      dt       Time step.
!  Output Arguments:
!      dvx      Change in x-component of particle velocity in time dt.
!      dvy      Change in y-component of particle velocity in time dt.
!      dux      Change in x-component of gas velocity in time dt, if
!          ldrag_pm_back_reaction = .false.; (1 + sum(eps)) times the
!          change in center-of-mass velocity, otherwise.
!      duy      Change in y-component of gas velocity in time dt, if
!          ldrag_pm_back_reaction = .false.; (1 + sum(eps)) times the
!          change in center-of-mass velocity, otherwise.
!  Optional Argument:
!      agx      x-component of the acceleartion of the gas.
!
      use Sub, only: one_minus_exp
!
      real, dimension(:), intent(in) :: eps
      real, dimension(size(eps)), intent(in) :: vx, vy
      real, intent(in) :: ux, uy, dt
      real, dimension(size(eps)), intent(out) :: dvx, dvy
      real, intent(out) :: dux, duy
      real, intent(in), optional :: agx
!
      real, dimension(size(eps)) :: vx0, vy0
      real :: epstot, dvg
      real :: ux0, uy0, uxe, uye, vxe, vye, vxcm, vycm
      real :: cosot, sinot1, sinot2
      real :: ot, t, ts, a0, a1, a2, a3, a4
!
!  Get the total solid-to-gas ratio.
!
      epstot = sum(eps)
!
!  Get the velocity reduction of the gas.
!
      if (present(agx)) then
        dvg = twoomega1 * (gx_gas + agx)
      else
        dvg = dv_gas
      endif
!
!  Get the deviation from the NSH solution.
!
      call get_nsh_solution(dvg, epstot, uxe, uye, vxe, vye)
      ux0 = ux - uxe
      uy0 = uy - uye
      vx0 = vx - vxe
      vy0 = vy - vye
!
!  Find the center-of-mass velocity of the particles.
!
      vcm: if (epstot > 0.0) then
        vxcm = sum(eps * vx0) / epstot
        vycm = sum(eps * vy0) / epstot
      else vcm
        vxcm = 0.0
        vycm = 0.0
      endif vcm
!
!  Get the epicycle motions.
!
      ot = epicycle_freq * dt
      cosot = cos(ot)
      a0 = sin(ot)
      sinot1 = a0 * epicycle_ratio
      sinot2 = a0 / epicycle_ratio
!
      uxe = ux0 * cosot + uy0 * sinot1
      uye = uy0 * cosot - ux0 * sinot2
!
      vxe = vxcm * cosot + vycm * sinot1
      vye = vycm * cosot - vxcm * sinot2
!
!  Find the velocity change for the particles.
!
      t = dt / tdrag
      a0 = exp(-t)
      a3 = 1.0 + epstot
      ts = a3 * t
      a4 = exp(-ts)
      if (abs(ts**4) > epsilon(1.0)) then
        a1 = (epstot + a4) / a3 - a0
      else
        a1 = 0.5 * epstot * t**2 * (1.0 - (t + ts) / 3.0)
      endif
      a2 = one_minus_exp(ts) / a3
!
      dvx = a1 * vxe + a2 * uxe + (a0 * (vx0 * cosot + vy0 * sinot1) - vx0)
      dvy = a1 * vye + a2 * uye + (a0 * (vy0 * cosot - vx0 * sinot2) - vy0)
!
!  Find the velocity change for the gas.
!
      gas: if (ldrag_pm_back_reaction) then
        uxe = ux0 + epstot * vxcm
        uye = uy0 + epstot * vycm
        a0 = ot * ot
        if (abs(a0) > epsilon(1.0)) then
          a0 = 1.0 - cosot
        else
          a0 = 0.5 * a0 * (1.0 - a0 / 12.0)
        endif
        dux = -a0 * uxe + sinot1 * uye
        duy = -a0 * uye - sinot2 * uxe
      else gas
        if (abs(ts * ts) > epsilon(1.0)) then
          a1 = (1.0 + epstot * a4) / a3
        else
          a1 = 1.0 - epstot * t
        endif
        a2 = epstot * a2
        dux = a1 * uxe + a2 * vxe - ux0
        duy = a1 * uye + a2 * vye - uy0
      endif gas
!
    endsubroutine drag_mutual_omega
!***********************************************************************
    pure subroutine gas_epicycle(ux, uy, dt, dux, duy)
!
!  Find the change in gas velocity due to rotation.
!
!  14-dec-15/ccyang: coded
!
!  Input Arguments
!      ux      x-component of the initial velocity of the gas.
!      uy      y-component of the initial velocity of the gas.
!      dt      Time step.
!  Output Arguments
!      dux     Change in x-component of the gas velocity in time dt.
!      duy     Change in y-component of the gas velocity in time dt.
!
      real, intent(in) :: ux, uy, dt
      real, intent(out) :: dux, duy
!
      real :: x, y, z
!
!  Find the epicycle.
!
      z = epicycle_freq * dt
      epicycle: if (z**4 > epsilon(1.0)) then
        x = 1.0 - cos(z)
        y = sin(z)
      else epicycle
        y = z**2
        x = 0.5 * y * (1.0 - y / 12.0)
        y = z * (1.0 - y / 6.0)
      endif epicycle
!
!  Find the velocity change of the gas.
!
      z = uy + dv_gas
      dux = -x * ux + y * epicycle_ratio * z
      duy = -x * z - y / epicycle_ratio * ux
!
    endsubroutine gas_epicycle
!***********************************************************************
    elemental subroutine get_nsh_solution(dv_gas, eps, ux, uy, vx, vy)
!
!  Finds the Nakagawa-Sekiya-Hayashi (1986) solution.
!
!  17-may-15/ccyang: coded.
!
!  Input Arguments
!      dv_gas  Velocity reduction of the gas due to radial pressure
!              gradient.
!      eps     Solid-to-gas mass ratio.
!
!  Output Arguments
!      ux      x-component of the velocity of the gas.
!      uy      y-component of the velocity of the gas.
!      vx      x-component of the velocity of the solids.
!      vy      y-component of the velocity of the solids.
!
      real, intent(in) :: dv_gas, eps
      real, intent(out) :: ux, uy, vx, vy
!
      real :: a, b
!
      a = 1.0 + eps
      b = dv_gas / (a**2 + taus2)
      vx = -twotaus * b
      vy = -a * b
      ux = -eps * vx
      uy = -(a + taus2) * b
!
    endsubroutine get_nsh_solution
!***********************************************************************
    elemental subroutine set_drag_equilibrium(cell)
!
!  Set the state of drag equilibrium to du and dv.
!
!  11-dec-15/ccyang: coded.
!
      type(pic), intent(inout) :: cell
!
      real, dimension(3) :: u, v
      integer :: i
      real :: eps
!
!  Find the equilibrium according to the drag type.
!
      u = 0.0
      v = 0.0
      eq: if (ldrag_on_par) then
        if (ldrag_on_gas) then
          eps = sum(cell%p%eps)
        else
          eps = 0.0
        endif
        call get_nsh_solution(dv_gas, eps, u(1), u(2), v(1), v(2))
        v = v + cell%u
      endif eq
!
!  Assign the solution to the velocity change.
!
      cell%du = u
      forall(i = 1:3) cell%p%dv(i) = v(i)
!
    endsubroutine set_drag_equilibrium
!***********************************************************************
!***********************************************************************
! LOCAL FUNCTIONS START BELOW HERE.
!***********************************************************************
!***********************************************************************
    function find_mp_swarm(eps_dtog) result(mp_swarm)
!
!  Finds the mass mp_swarm of each super-particle given the global dust-
!  to-gas density ratio eps_dtog.
!
!  13-sep-15/ccyang: coded.
!
      use EquationOfState, only: gamma, rho0, cs0
!
      real :: mp_swarm
      real, intent(in) :: eps_dtog
!
      real :: mass
!
!  Find the total gas mass.
!
      gasmass: if (gz_par_coeff > 0.0) then
!       For stratified gas
        mass = sqrt(twopi / gamma) * (cs0 / Omega) * product(lxyz(1:2), lactive_dimension(1:2))
      else gasmass
!       For uniform gas
        mass = product(lxyz, lactive_dimension)
      endif gasmass
      mass = rho0 * mass
!
!  Uniformly distribute the mass to the particles.
!
      mp_swarm = eps_dtog * mass / real(npar)
!
    endfunction find_mp_swarm
!***********************************************************************
    elemental real function particle_zaccel(zp) result(gz)
!
!  Returns the z-acceleration of a particle.
!
!  22-jun-15/ccyang: coded.
!
!  Currently, only supports -gz_par_coeff^2 * zp.
!
!  Input Argument:
!      zp  Vertical position of the particle.
!
      real, intent(in) :: zp
!
      gz = -gzpc * zp
!
    endfunction particle_zaccel
!***********************************************************************
endmodule Particles_drag
