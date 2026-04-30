! $Id$
!
!  Secondary breakup of Lagrangian liquid parcels.
!
!  Three independent models are available:
!
!    KH (Kelvin-Helmholtz) and RT (Rayleigh-Taylor) - Aguerre & Nigro,
!         Computers and Fluids 190 (2019) 30-48, sections 3.2-3.3 and
!         Appendices B3-B4.  These are designed for very high Weber
!         numbers typical of fuel-injector sprays and can be (and by
!         default now are) gated on a minimum gas Weber number to keep
!         them inside their regime of validity.
!
!    PE (Pilch-Erdman) - Pilch, M. & Erdman, C.A., Use of breakup time
!         data and velocity history data to predict the maximum size of
!         stable fragments for acceleration-induced breakup of a liquid
!         drop, Int. J. Multiphase Flow 13(6) (1987) 741-757.  This is
!         a regime-based (vibrational / bag / bag-stamen / sheet
!         stripping / catastrophic) acceleration-driven model that
!         covers the full Weber number range, including the bag and
!         bag-stamen regime that KH+RT do not capture.
!         The breakup-time scaling and Schiller-Naumann drag closure
!         used here also follow Apte, Gorokhovski & Moin, IJMF 29
!         (2003) 1503-1522.
!
!  When the PE branch is enabled it has priority over KH and RT inside
!  particles_breakup_pencils: at most one of the three models acts on
!  any given parcel in any given step.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! CPARAM logical, parameter :: lparticles_breakup = .true.
!
! MPVAR CONTRIBUTION 0
! MPAUX CONTRIBUTION 3
!
!***************************************************************
module Particles_breakup
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map, only: interpolate_linear
  use Particles_number, only: set_particle_number
  use Particles_radius, only: set_particle_radius
  use Particles_sub, only: append_npaux, sum_par_name, max_par_name
!
  implicit none
!
  include 'particles_breakup.h'
!
  logical :: lparticles_breakup_kh = .false.
  logical :: lparticles_breakup_rt = .false.
  logical :: lparticles_breakup_pe = .false.
!  Verbose flag: when .true., a one-line summary is written to stdout at
!  every individual breakup event, naming the mechanism (PE bag /
!  bag-stamen / sheet stripping / catastrophic; RT; KH stripping or KH
!  full conversion) along with the parcel id, time, parent and child
!  diameters and the local gas Weber number.  Off by default to keep
!  log files quiet for production runs.
  logical :: lprint_breakup_mode = .false.
  real :: tstart_particles_breakup = 0.0
  real :: kh_B0 = 0.61
  real :: kh_B1 = 40.0
  real :: kh_B2 = 0.05
  real :: rt_CRT = 0.1
  real :: rt_Ctau = 1.0
  real :: sigma_liq = 7.2e-2
  real :: kh_min_child_npswarm = 0.0
  real :: nu_liq = 1e-2    !PAR_DOC: kinematic viscosity of the liquid (in cgs units)
!
!  Weber-number floors for the KH and RT branches (Pilch-Erdman regime
!  thresholds).  Default to 0.0, which reproduces the original Aguerre-Nigro
!  behaviour.  Setting kh_we_min~100 restricts KH to the sheet-stripping
!  regime; rt_we_min~350 restricts RT to the catastrophic regime.
  real :: kh_we_min = 0.0
  real :: rt_we_min = 0.0
!
!  Pilch-Erdman parameters.
  real :: pe_we_crit = 12.0          ! base critical Weber number, Eq. (5)
  real :: pe_child_factor = 0.5      ! child-parcel size as a fraction of d_max
                                     !  (~0.5 = mass-median, 1.0 = max stable)
!  Minimum parent diameter at which PE may still fire.  Default is one
!  micrometre (1e-4 cm in cgs) - well below the meteorological raindrop
!  band but still above the integrator-stiff regime that develops once a_p
!  drops to subnormal.  Setting this to 0 reproduces the unguarded
!  cascade.
  real :: pe_min_diameter = 1.0e-12
!  Lower bound on the post-breakup diameter; serves as a final numerical
!  guard so that round-off in compute_pe_scales (e.g. an over-large urmag
!  due to a stiff drag step) cannot drive a_p into subnormal range.
  real :: pe_floor_diameter = 1.0e-12
!  Maximum factor by which a single breakup event is allowed to shrink
!  the parent diameter.  Caps catastrophic-regime events that would
!  otherwise feed the runaway described above; the cascade is then
!  spread over multiple breakup events (each separated by tau_PE) which
!  is what Pilch-Erdman themselves describe as a multistage process.
  real :: pe_max_shrink_ratio = 0.0001
!
  integer :: rt_method = 2
  integer :: idiag_imskhm = 0
  integer :: idiag_tbpem = 0
  integer :: idiag_wegm = 0       ! mean per-parcel gas Weber number
  integer :: idiag_wegmin = 0     ! minimum per-parcel gas Weber number
  integer :: idiag_wegmax = 0     ! maximum per-parcel gas Weber number
!
  namelist /particles_breakup_init_pars/ &
      lparticles_breakup_kh, lparticles_breakup_rt, lparticles_breakup_pe, &
      tstart_particles_breakup, &
      kh_B0, kh_B1, kh_B2, rt_CRT, rt_Ctau, rt_method, sigma_liq, &
      kh_min_child_npswarm, nu_liq, kh_we_min, rt_we_min, &
      pe_we_crit, pe_child_factor, pe_min_diameter, &
      pe_floor_diameter, pe_max_shrink_ratio, &
      lprint_breakup_mode
!
  namelist /particles_breakup_run_pars/ &
      lparticles_breakup_kh, lparticles_breakup_rt, lparticles_breakup_pe, &
      tstart_particles_breakup, &
      kh_B0, kh_B1, kh_B2, rt_CRT, rt_Ctau, rt_method, sigma_liq, &
      kh_min_child_npswarm, nu_liq, kh_we_min, rt_we_min, &
      pe_we_crit, pe_child_factor, pe_min_diameter, &
      pe_floor_diameter, pe_max_shrink_ratio, &
      lprint_breakup_mode
!
contains
!***********************************************************************
    subroutine register_particles_breakup()
!
!  Register particle-local breakup state.
!
      call append_npaux('imskh', imskh)
      call append_npaux('itbrt', itbrt)
      call append_npaux('itbpe', itbpe)
!
    endsubroutine register_particles_breakup
!***********************************************************************
    subroutine initialize_particles_breakup(f)
!
!  Check the breakup configuration.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      if (lparticles_breakup_kh .or. lparticles_breakup_rt .or. &
          lparticles_breakup_pe) then
        if (.not. lparticles_radius) call fatal_error('initialize_particles_breakup', &
            'particle breakup requires Particles_radius')
        if (.not. lparticles_number) call fatal_error('initialize_particles_breakup', &
            'particle breakup requires Particles_number')
        if (lparticles_blocks) call not_implemented('initialize_particles_breakup', &
            'particle breakup for particles_blocks')
      endif
      if (lparticles_breakup_kh) then
        if (kh_B0 <= 0.0) call fatal_error('initialize_particles_breakup', &
            'kh_B0 must be positive')
        if (kh_B1 <= 0.0) call fatal_error('initialize_particles_breakup', &
            'kh_B1 must be positive')
        if (kh_B2 <= 0.0) call fatal_error('initialize_particles_breakup', &
            'kh_B2 must be positive')
      endif
      if (lparticles_breakup_rt) then
        if (rt_CRT <= 0.0) call fatal_error('initialize_particles_breakup', &
            'rt_CRT must be positive')
        if (rt_Ctau <= 0.0) call fatal_error('initialize_particles_breakup', &
            'rt_Ctau must be positive')
        if (rt_method /= 1 .and. rt_method /= 2) call fatal_error('initialize_particles_breakup', &
            'rt_method must be 1 or 2')
      endif
      if (lparticles_breakup_pe) then
        if (pe_we_crit <= 0.0) call fatal_error('initialize_particles_breakup', &
            'pe_we_crit must be positive')
        if (pe_child_factor <= 0.0 .or. pe_child_factor > 1.0) &
            call fatal_error('initialize_particles_breakup', &
            'pe_child_factor must lie in (0,1]')
        if (pe_floor_diameter <= 0.0) call fatal_error( &
            'initialize_particles_breakup', &
            'pe_floor_diameter must be positive (it is the lower bound on a_p)')
        if (pe_max_shrink_ratio <= 0.0 .or. pe_max_shrink_ratio > 1.0) &
            call fatal_error('initialize_particles_breakup', &
            'pe_max_shrink_ratio must lie in (0,1]')
      endif
      if (kh_we_min < 0.0 .or. rt_we_min < 0.0) call fatal_error( &
          'initialize_particles_breakup', 'KH/RT Weber floors must be >= 0')
      if (lparticles_breakup_kh .or. lparticles_breakup_rt .or. &
          lparticles_breakup_pe) then
        if (sigma_liq <= 0.0) call fatal_error('initialize_particles_breakup', &
            'sigma_liq must be positive')
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_breakup
!***********************************************************************
    subroutine read_particles_breakup_init_pars(iostat)
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_breakup_init_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_breakup_init_pars
!***********************************************************************
    subroutine write_particles_breakup_init_pars(unit)
      integer, intent(in) :: unit
      integer :: stat
!
      write(unit, NML=particles_breakup_init_pars, IOSTAT=stat)
      if (stat /= 0) call fatal_error('write_particles_breakup_init_pars', &
          'cannot write particles_breakup_init_pars', force=.true.)
!
    endsubroutine write_particles_breakup_init_pars
!***********************************************************************
    subroutine read_particles_breakup_run_pars(iostat)
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=particles_breakup_run_pars, IOSTAT=iostat)
!
    endsubroutine read_particles_breakup_run_pars
!***********************************************************************
    subroutine write_particles_breakup_run_pars(unit)
      integer, intent(in) :: unit
      integer :: stat
!
      write(unit, NML=particles_breakup_run_pars, IOSTAT=stat)
      if (stat /= 0) call fatal_error('write_particles_breakup_run_pars', &
          'cannot write particles_breakup_run_pars', force=.true.)
!
    endsubroutine write_particles_breakup_run_pars
!***********************************************************************
    subroutine rprint_particles_breakup(lreset,lwrite)
      use Diagnostics

      logical, intent(in) :: lreset
      logical, optional, intent(in) :: lwrite
      integer :: iname
!
      if (lreset) then
        idiag_imskhm = 0
        idiag_tbpem = 0
        idiag_wegm = 0
        idiag_wegmin = 0
        idiag_wegmax = 0
      endif
!
      if (lroot .and. ip < 14) print*,'rprint_particles_breakup: run through parse list'
      do iname = 1,nname
        call parse_name(iname,cname(iname),cform(iname),'imskhm',idiag_imskhm)
        call parse_name(iname,cname(iname),cform(iname),'tbpem',idiag_tbpem)
        call parse_name(iname,cname(iname),cform(iname),'wegm',idiag_wegm)
        call parse_name(iname,cname(iname),cform(iname),'wegmin',idiag_wegmin)
        call parse_name(iname,cname(iname),cform(iname),'wegmax',idiag_wegmax)
      enddo
!
      call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_breakup
!***********************************************************************
    subroutine dbreakup_dt(f,df,fp,dfp,ineargrid)
!
!  Diagnostics for breakup-owned particle state.
!  This follows the standard particle-module pattern: diagnostics are gathered
!  during the regular particle RHS phase, while the actual discrete breakup
!  update is still applied later in particles_discrete_collisions.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mpar_loc,mpvar), intent(inout) :: dfp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
      real, dimension(mpar_loc) :: we_g_arr
      real, dimension(3) :: uup
      real :: rho_gas, lnrho_gas, urmag, dp_loc
      integer :: k
      logical :: lneed_we
!
      if (ldiagnos) then
        if (idiag_imskhm /= 0) call sum_par_name(fp(1:npar_loc,imskh),idiag_imskhm)
        if (idiag_tbpem /= 0) call sum_par_name(fp(1:npar_loc,itbpe),idiag_tbpem)
!
!  Per-parcel gas Weber number (We_g = rho_gas * u_r^2 * d_p / sigma).
!  Computed on demand, since it is not a stored aux variable: each
!  diagnostic call costs one gas-state interpolation per parcel.
        lneed_we = idiag_wegm /= 0 .or. idiag_wegmin /= 0 .or. idiag_wegmax /= 0
        if (lneed_we) then
          we_g_arr = 0.0
          do k = 1, npar_loc
            if (fp(k,iap) <= 0.0) cycle
            if (fp(k,inpswarm) <= 0.0) cycle
            call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup, &
                ineargrid(k,:),0,ipar(k))
            call get_local_gas_density(f,fp(k,ixp:izp),ineargrid(k,:), &
                ipar(k),rho_gas,lnrho_gas)
            urmag = sqrt(sum((uup - fp(k,ivpx:ivpz))**2))
            dp_loc = 2.0*fp(k,iap)
            we_g_arr(k) = rho_gas*urmag**2*dp_loc/max(sigma_liq,tini)
          enddo
          if (idiag_wegm /= 0) call sum_par_name(we_g_arr(1:npar_loc),idiag_wegm)
          if (idiag_wegmax /= 0) call max_par_name(we_g_arr(1:npar_loc),idiag_wegmax)
          if (idiag_wegmin /= 0) call max_par_name(-we_g_arr(1:npar_loc), &
              idiag_wegmin,lneg=.true.)
        endif
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dbreakup_dt
!***********************************************************************
    subroutine particles_breakup_pencils(f,fp,ineargrid)
!
!  Apply secondary-breakup updates after particle advection for the parcels
!  that exist at the start of this step.  Newly created child parcels are
!  appended to fp and take part from the next step on.
!
!  Branch priority (highest first):
!    1. Pilch-Erdman (PE) - regime-aware, covers vibrational / bag /
!       bag-stamen / sheet-stripping / catastrophic regimes.
!    2. Rayleigh-Taylor (RT) - Aguerre & Nigro (Appendix B4) gated by
!       rt_we_min.
!    3. Kelvin-Helmholtz (KH) - Aguerre & Nigro (Appendix B3) gated by
!       kh_we_min, only relevant for sheet-stripping.
!
!  At most one branch acts on any parcel in any step.  When PE is enabled
!  it preempts RT and KH for that parcel; otherwise the legacy KH+RT path
!  is followed exactly as in the original Aguerre-Nigro implementation.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
!
      real, dimension(3) :: uup
      real :: rho_gas, lnrho_gas, nu_gas, urmag
      real :: dp_old, dp_new, dp_child, diameter_ratio
      real :: we_g, we_p, rep, oh, tay, lambda_kh, omega_kh, tau_kh, alpha
      real :: lambda_rt, omega_rt, tau_rt, acc_rt, cd_drag
      real :: tau_pe, dp_max_stable, we_g_pe, oh_pe
      real :: mp_old, mp_child, ms_new, ms_inc, np_child
      integer :: k, npar_loc_old, pe_regime
!
      if (.not. (lparticles_breakup_kh .or. lparticles_breakup_rt .or. &
                 lparticles_breakup_pe)) return
      if (.not. llast) return
      if (t < tstart_particles_breakup) return
!
!  Freeze the loop upper bound so that daughters created during this call do not
!  themselves fragment before the next particle step.
      npar_loc_old = npar_loc
!
      do k = 1, npar_loc_old
!
!  Ignore deleted / empty parcels.
        if (fp(k,iap) <= 0.0) cycle
        if (fp(k,inpswarm) <= 0.0) cycle
!
!  Interpolate local carrier-gas state to the particle position.
        call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:),0,ipar(k))
        call get_local_gas_density(f,fp(k,ixp:izp),ineargrid(k,:),ipar(k),rho_gas,lnrho_gas)
        call get_local_gas_kinematic_viscosity(rho_gas,nu_gas)
!
        urmag = sqrt(sum((uup - fp(k,ivpx:ivpz))**2))
        if (urmag <= tini) cycle
!
!  Work in terms of droplet diameter because the breakup correlations in
!  the underlying papers are written for d_p, not for the stored parcel
!  radius a_p.
        dp_old = 2.0*fp(k,iap)
        if (dp_old <= 0.0) cycle
!
!  Compute the local gas Weber number once; it is the gate for all three
!  branches and the natural specialisation criterion.
        we_g = rho_gas * urmag**2 * dp_old / max(sigma_liq,tini)
!
!  Branch priority (highest first):
!    1. RT (catastrophic specialist) - fires when lparticles_breakup_rt
!       and We_g >= rt_we_min.
!    2. KH (sheet-stripping specialist) - fires when lparticles_breakup_kh
!       and We_g >= kh_we_min and RT did not claim the parcel.
!    3. PE (regime-aware fallback) - fires when lparticles_breakup_pe
!       and neither RT nor KH claimed the parcel.
!  At most one branch acts on any parcel in any step.  When a branch owns
!  the parcel, the OTHER branches' accumulators are reset so that stale
!  values cannot affect a future hand-back.
!
!  -------- Branch 1 : RT (Aguerre & Nigro, App. B4) -----------------
        if (lparticles_breakup_rt .and. we_g >= rt_we_min) then
          fp(k,imskh) = 0.0
          fp(k,itbpe) = 0.0
          call compute_rt_scales(dp_old,rho_gas,nu_gas,urmag,cd_drag,acc_rt,lambda_rt,omega_rt,tau_rt)
!  RT is supposed to be active in this branch; degenerate scales are a
!  red flag, not a "skip silently" condition.
          if (tau_rt <= 0.0 .or. lambda_rt <= 0.0) call fatal_error( &
              'particles_breakup_pencils', &
              'RT scales degenerate (tau_rt or lambda_rt <= 0) above rt_we_min')
          if (dp_old > lambda_rt) then
            fp(k,itbrt) = fp(k,itbrt) + dt
            if (fp(k,itbrt) > tau_rt) then
              if (lprint_breakup_mode) then
                write(*,'(a,i0,a,1pe12.4,a,1pe10.3,a,1pe10.3,a,1pe10.3,a,1pe10.3,a,1pe10.3)') &
                    'particles_breakup: ipar=', ipar(k), ' t=', t, &
                    ' RT catastrophic  d_p: ', dp_old, &
                    ' -> ', lambda_rt, ' cm  urmag=', urmag, &
                    ' cm/s  rho_g=', rho_gas, '  We_g=', we_g
              endif
              call apply_rt_breakup(fp,k,dp_old,lambda_rt)
            endif
          else
            fp(k,itbrt) = 0.0
          endif
          cycle
        endif
!
!  -------- Branch 2 : KH (Aguerre & Nigro, App. B3) -----------------
        if (lparticles_breakup_kh .and. we_g >= kh_we_min) then
          fp(k,itbrt) = 0.0
          fp(k,itbpe) = 0.0
          call compute_kh_scales(dp_old,rho_gas,nu_gas,urmag,we_g,we_p,rep,oh,tay,lambda_kh,omega_kh,tau_kh)
          if (tau_kh <= 0.0 .or. lambda_kh <= 0.0) cycle
!
!  Section 3.2 / Appendix B3: the stable daughter size is d_KH = 2 B0 lambda_KH.
          dp_child = 2.0*kh_B0*lambda_kh
          if (dp_child <= 0.0) cycle
          if (dp_old <= dp_child) cycle
!
!  Backward-Euler relaxation of the mother diameter toward d_KH over tau_KH.
          alpha = dt / tau_kh
          dp_new = (alpha*dp_child + dp_old) / (1.0 + alpha)
          dp_new = min(dp_old, max(dp_child, dp_new))
!
!  The removed mass is not emitted immediately. It is accumulated until it
!  either exceeds the breakup tolerance B2 m_p or strips the full mother mass.
          mp_old = mass_from_diameter(dp_old)
          ms_inc = mass_from_diameter(dp_old) - mass_from_diameter(dp_new)
          ms_new = fp(k,imskh) + ms_inc
!
            if (ms_new > mp_old) then
!
!  Full conversion branch: the parcel is interpreted as having completely
!  transformed into KH-stable droplets of size d_child. We keep one parcel
!  and adjust its internal droplet multiplicity accordingly.
!
            if (lprint_breakup_mode) then
              write(*,'(a,i0,a,1pe12.4,a,1pe10.3,a,1pe10.3,a,1pe10.3,a,1pe10.3,a,1pe10.3)') &
                  'particles_breakup: ipar=', ipar(k), ' t=', t, &
                  ' KH full-conversion  d_p: ', dp_old, &
                  ' -> ', dp_child, ' cm  urmag=', urmag, &
                  ' cm/s  rho_g=', rho_gas, '  We_g=', we_g
            endif
            diameter_ratio = dp_old / dp_child
            fp(k,inpswarm) = fp(k,inpswarm) * diameter_ratio**3
            fp(k,iap) = 0.5*dp_child
            fp(k,imskh) = 0.0
            fp(k,itbrt) = 0.0
            call update_particle_mass_fields(fp,k)
          elseif (ms_new >= kh_B2*mp_old) then
!
!  Child-emission branch: convert the accumulated stripped mass into one new
!  super-droplet parcel, while the mother retains the remaining liquid mass.
!
            mp_child = mass_from_diameter(dp_child)
            np_child = fp(k,inpswarm) * ms_new / mp_child
!
            if (lprint_breakup_mode) then
              write(*,'(a,i0,a,1pe12.4,a,1pe10.3,a,1pe10.3,a,1pe10.3,a,1pe10.3,a,1pe10.3)') &
                  'particles_breakup: ipar=', ipar(k), ' t=', t, &
                  ' KH stripping        d_p: ', dp_old, &
                  ' -> ', dp_child, ' cm (child)  urmag=', urmag, &
                  ' cm/s  rho_g=', rho_gas, '  We_g=', we_g
            endif
            fp(k,iap) = 0.5 * diameter_from_mass(max(mp_old - ms_new, 0.0))
            fp(k,imskh) = 0.0
            fp(k,itbrt) = 0.0
            call update_particle_mass_fields(fp,k)
!
            if (np_child > kh_min_child_npswarm) then
              call insert_child_particle(fp,ineargrid,k,dp_child,np_child)
            endif
          else
!
!  Not enough stripped mass yet: keep carrying it on the mother parcel.
!
            fp(k,imskh) = ms_new
          endif
          cycle    ! KH owned this parcel for the step
        endif
!
!  -------- Branch 3 : PE (Pilch-Erdman) -----------------------------
!  Catches everything below rt_we_min and below kh_we_min - i.e. the bag,
!  bag-stamen and (when KH is off) sheet-stripping regimes.
!
!  Active guards: pe_min_diameter, pe_floor_diameter, pe_max_shrink_ratio
!  (these only constrain the post-breakup size; they never silently skip
!  a breakup event).  If tau_pe < dt the run aborts with a fatal_error
!  rather than silently dropping the breakup, on the assumption that
!  such a regime indicates either an inappropriate time step or a deeper
!  numerical pathology that the user should see.
        if (lparticles_breakup_pe) then
          fp(k,itbrt) = 0.0
          fp(k,imskh) = 0.0
          if (dp_old > max(pe_min_diameter, pe_floor_diameter)) then
            call compute_pe_scales(dp_old,rho_gas,urmag,oh_pe,we_g_pe, &
                                   tau_pe,dp_max_stable,pe_regime)
            if (pe_regime > 0) then
!  Refuse to advance the breakup accumulator with a step longer than the
!  characteristic breakup time; either the time step is too long or the
!  parcel has fallen into a stiff regime that requires investigation.
              if (tau_pe < dt) call fatal_error( &
                  'particles_breakup_pencils', &
                  'PE total breakup time tau_pe < dt; reduce dt or check setup')
              fp(k,itbpe) = fp(k,itbpe) + dt
              if (fp(k,itbpe) > tau_pe) then
                dp_new = pe_child_factor * dp_max_stable
                dp_new = min(dp_new, dp_old)
                dp_new = max(dp_new, pe_max_shrink_ratio * dp_old)
                dp_new = max(dp_new, pe_floor_diameter)
                if (lprint_breakup_mode) then
                  write(*,'(a,i0,a,1pe12.4,a,a,a,1pe10.3,a,1pe10.3,a,1pe10.3,a,1pe10.3,a,1pe10.3)') &
                      'particles_breakup: ipar=', ipar(k), ' t=', t, ' PE ', &
                      trim(pe_regime_name(pe_regime)), &
                      '  d_p: ', dp_old, ' -> ', dp_new, &
                      ' cm  urmag=', urmag, ' cm/s  rho_g=', rho_gas, &
                      '  We_g=', we_g_pe
                endif
                call apply_pe_breakup(fp,k,dp_old,dp_new)
              endif
            else
              fp(k,itbpe) = 0.0
            endif
          else
            fp(k,itbpe) = 0.0
          endif
          cycle
        endif
!
!  No model active for this parcel - nothing to do.  Reset all
!  accumulators so they cannot resurrect when a branch becomes active
!  again later.
        fp(k,itbrt) = 0.0
        fp(k,imskh) = 0.0
        fp(k,itbpe) = 0.0
      enddo
!
    endsubroutine particles_breakup_pencils
!***********************************************************************
    subroutine get_local_gas_density(f,xxp,inear,ipid,rho_gas,lnrho_gas)
      use EquationOfState, only: rho0
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(3), intent(in) :: xxp
      integer, dimension(3), intent(in) :: inear
      integer, intent(in) :: ipid
      real, intent(out) :: rho_gas, lnrho_gas
!
!  Support both density formulations used in the code base. If neither rho nor
!  lnrho is present, fall back to the EOS reference density.
      if (irho /= 0) then
        call interpolate_linear(f,irho,xxp,rho_gas,inear,0,ipid)
        rho_gas = max(rho_gas, tini)
        lnrho_gas = log(rho_gas)
      elseif (ilnrho /= 0) then
        call interpolate_linear(f,ilnrho,xxp,lnrho_gas,inear,0,ipid)
        rho_gas = max(exp(lnrho_gas), tini)
      else
        rho_gas = rho0
        lnrho_gas = log(rho_gas)
      endif
!
    endsubroutine get_local_gas_density
!***********************************************************************
    subroutine get_local_gas_kinematic_viscosity(rho_gas,nu_gas)
      use Viscosity, only: getnu
!
      real, intent(in) :: rho_gas
      real, intent(out) :: nu_gas
!
      character(len=labellen) :: ivis
      real :: nu_input
!
      ivis = ''
      nu_input = 0.0
      call getnu(nu_input=nu_input,ivis=ivis)
!
!  The sample setup uses a constant gas viscosity. For density-scaled constant
!  viscosities, reconstruct the local kinematic viscosity from rho_gas. More
!  elaborate viscosity models fall back to the base coefficient here.
      select case (trim(ivis))
      case ('nu-const')
        nu_gas = nu_input
      case ('rho-nu-const')
        nu_gas = nu_input / max(rho_gas,tini)
      case ('sqrtrho-nu-const')
        nu_gas = nu_input / sqrt(max(rho_gas,tini))
      case default
        nu_gas = nu_input
      end select
!
      nu_gas = max(nu_gas,tini)
!
    endsubroutine get_local_gas_kinematic_viscosity
!***********************************************************************
    subroutine compute_kh_scales(dp,rho_gas,nu_gas,urmag,we_g,we_p,rep,oh,tay,lambda_kh,omega_kh,tau_kh)
      real, intent(in) :: dp, rho_gas, nu_gas, urmag
      real, intent(out) :: we_g, we_p, rep, oh, tay, lambda_kh, omega_kh, tau_kh
      real :: denom, rep_liq
!
!  Dimensionless groups and KH correlations used in the breakup model.
!  These are the quantities needed to reconstruct the stable KH scale d_KH
!  and its response time tau_KH.
      we_g = rho_gas * urmag**2 * dp / sigma_liq
      we_p = rhopmat * urmag**2 * dp / sigma_liq
      rep = urmag * dp / max(nu_gas,tini)
      rep_liq = urmag * dp / max(nu_liq,tini)
!
      oh = sqrt(max(we_p,0.0)) / max(rep_liq,tini)
      tay = oh * sqrt(max(we_g,0.0))
!      
      lambda_kh = 4.51 * dp * (1.0 + 0.45*sqrt(max(oh,0.0))) * (1.0 + 0.4*tay**0.7) / &
           (1.0 + 0.865*we_g**1.67)**0.6
!
      denom = (1.0 + oh) * (1.0 + 1.4*tay**0.6)
      omega_kh = (0.34 + 0.385*we_g**1.5) / max(denom,tini) * &
          sqrt(8.0*sigma_liq / max(rhopmat*dp**3,tini))
!
      tau_kh = 3.788 * kh_B1 * dp / max(2.0*omega_kh*lambda_kh,tini)
!
    endsubroutine compute_kh_scales
!***********************************************************************
    subroutine compute_drag_coefficient(rep,cd_drag)
      real, intent(in) :: rep
      real, intent(out) :: cd_drag
!
      if (rep <= tini) then
        cd_drag = 0.0
      elseif (rep <= 1000.0) then
        cd_drag = 24.0/rep * (1.0 + 0.15*rep**0.687)
      else
        cd_drag = 0.44
      endif
!
    endsubroutine compute_drag_coefficient
!***********************************************************************
    subroutine compute_rt_scales(dp,rho_gas,nu_gas,urmag,cd_drag,acc_rt,lambda_rt,omega_rt,tau_rt)
      real, intent(in) :: dp, rho_gas, nu_gas, urmag
      real, intent(out) :: cd_drag, acc_rt, lambda_rt, omega_rt, tau_rt
!
      real :: rep, density_jump
!
      rep = urmag * dp / max(nu_gas,tini)
      call compute_drag_coefficient(rep,cd_drag)
!
      acc_rt = 0.75 * cd_drag * rho_gas * urmag**2 / max(rhopmat*dp,tini)
      density_jump = max(rhopmat - rho_gas, 0.0)
!
      if (acc_rt <= 0.0 .or. density_jump <= 0.0) then
        lambda_rt = 0.0
        omega_rt = 0.0
        tau_rt = 0.0
        return
      endif
!
      lambda_rt = 2.0*pi*rt_CRT * sqrt(3.0*sigma_liq / max(acc_rt*density_jump,tini))
      omega_rt = sqrt(2.0*(acc_rt*density_jump)**1.5 / max(3.0*sqrt(3.0)*sigma_liq,tini))
      tau_rt = rt_Ctau / max(omega_rt,tini)
!
    endsubroutine compute_rt_scales
!***********************************************************************
    subroutine apply_rt_breakup(fp,k,dp_old,lambda_rt)
!
!  RT breakup follows Appendix B4 of the paper. The parcel remains single and
!  its droplet multiplicity is increased to conserve represented liquid mass.
!
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, intent(in) :: k
      real, intent(in) :: dp_old, lambda_rt
!
      real :: dp_new, diameter_ratio
!
      if (lambda_rt <= 0.0) return
!
      select case (rt_method)
      case (1)
        dp_new = lambda_rt
      case default
        dp_new = dp_old**(2.0/3.0) * lambda_rt**(1.0/3.0)
      end select
!
      dp_new = min(dp_old, max(lambda_rt, dp_new))
      diameter_ratio = dp_old / max(dp_new,tini)
!
      fp(k,inpswarm) = fp(k,inpswarm) * diameter_ratio**3
      fp(k,iap) = 0.5*dp_new
      fp(k,imskh) = 0.0
      fp(k,itbrt) = 0.0
      fp(k,itbpe) = 0.0
      call update_particle_mass_fields(fp,k)
!
    endsubroutine apply_rt_breakup
!***********************************************************************
    subroutine compute_pe_scales(dp,rho_gas,urmag,oh,we_g,tau_pe,dp_max_stable,regime)
!
!  Pilch-Erdman regime classification, total-breakup-time and
!  maximum-stable-fragment-size correlations (Pilch & Erdman, IJMF 13(6),
!  1987; equations referenced are from that paper).  The breakup-time
!  scaling and the Schiller-Naumann drag closure used here are also the
!  ones recommended by Apte, Gorokhovski & Moin, IJMF 29 (2003).
!
!  All quantities are returned in the simulation's physical units.
!
!  Output regimes:
!     0  vibrational / no breakup (We_g < We_c)
!     1  bag breakup                (We_c <= We_g <= 50)
!     2  bag-stamen / multimode    (50  < We_g <= 100)
!     3  sheet stripping            (100 < We_g <= 350)
!     4  catastrophic breakup       (We_g > 350)
!
      real, intent(in) :: dp, rho_gas, urmag
      real, intent(out) :: oh, we_g, tau_pe, dp_max_stable
      integer, intent(out) :: regime
!
      real :: rep_liq, eps, td, vd_dim, ur_at_end, we_c, we_excess
      real, parameter :: cd_pe = 0.5     ! Eq. (23): rigid-sphere C_D
      real, parameter :: bb_pe = 0.0758  ! Eq. (23): empirical B (incompressible)
!
      we_g = rho_gas * urmag**2 * dp / max(sigma_liq,tini)
!  Ohnesorge number (Eq. (2)):  Oh = mu_l / sqrt(rho_l sigma d_p)
!  expressed in cgs-friendly form via the liquid kinematic viscosity.
      rep_liq = urmag * dp / max(nu_liq,tini)
      oh = sqrt(max(rhopmat * urmag**2 * dp / max(sigma_liq,tini), 0.0)) / max(rep_liq,tini)
!
!  Eq. (5): viscosity-corrected critical Weber number.
      we_c = pe_we_crit * (1.0 + 1.077 * oh**1.6)
!
      regime = 0
      tau_pe = 0.0
      dp_max_stable = dp
      if (we_g <= we_c) return
!
!  Density ratio used in the dimensionless time T = t * V/D * sqrt(eps),
!  Eq. (3) with eps = rho_g/rho_p, Eq. (4).
      eps = rho_gas / max(rhopmat,tini)
      we_excess = max(we_g - 12.0, tini)
!
!  Total-breakup-time correlations Eqs. (8)-(12) and regime classification.
      if (we_g <= 18.0) then
        td = 6.0 * we_excess**(-0.25)
        regime = 1
      else if (we_g <= 45.0) then
        td = 2.45 * we_excess**(0.25)
        regime = 1
      else if (we_g <= 100.0) then
        td = 14.1 * we_excess**(-0.25)
        regime = 2
      else if (we_g <= 351.0) then
        td = 14.1 * we_excess**(-0.25)
        regime = 3
      else if (we_g <= 2670.0) then
        td = 0.766 * we_excess**(0.25)
        regime = 4
      else
        td = 5.5
        regime = 4
      endif
!
!  Convert dimensionless time to physical time.
!     T = t * V/D * sqrt(eps)  =>  t = T * D/V / sqrt(eps)
      tau_pe = td * dp / max(urmag,tini) / sqrt(max(eps,tini))
!
!  Improved-velocity correlation, Eqs. (20) and (23) (incompressible flow):
!     V_d / V * sqrt(rho_p/rho_g) = 0.75 C_d T + 3 B T^2
      vd_dim = 0.75*cd_pe*td + 3.0*bb_pe*td**2
!  Drop velocity in the gas frame at the end of the breakup window:
!     V_d = V * sqrt(rho_g/rho_p) * vd_dim,
!  hence the relative velocity at end of breakup is
!     V - V_d = V * (1 - sqrt(eps) * vd_dim).
      ur_at_end = urmag * max(1.0 - sqrt(eps)*vd_dim, 0.0)
!
!  Eq. (33): improved estimate of the maximum stable fragment diameter.
      if (ur_at_end <= tini) then
        dp_max_stable = dp
      else
        dp_max_stable = we_c * sigma_liq / max(rho_gas*ur_at_end**2, tini)
      endif
      dp_max_stable = min(dp_max_stable, dp)
!
    endsubroutine compute_pe_scales
!***********************************************************************
    subroutine apply_pe_breakup(fp,k,dp_old,dp_new)
!
!  Realise a Pilch-Erdman breakup event on parcel k.  We retain a single
!  parcel and inflate its droplet multiplicity to conserve liquid mass,
!  in the same spirit as the RT branch in apply_rt_breakup.  This avoids
!  generating an unbounded number of child parcels for runs where many
!  parcels reach the catastrophic regime simultaneously.
!
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, intent(in) :: k
      real, intent(in) :: dp_old, dp_new
!
      real :: diameter_ratio, d_new
!
!  The d_new value is also clamped at pe_floor_diameter (and tini) to
!  guarantee that the new radius cannot drop into the subnormal regime
!  even if the caller passes a pathologically small target.  This is the
!  last-line numerical guard for the cascade.
      d_new = max(dp_new, pe_floor_diameter, tini)
      d_new = min(d_new, dp_old)
      diameter_ratio = dp_old / d_new
!
      fp(k,inpswarm) = fp(k,inpswarm) * diameter_ratio**3
      fp(k,iap) = 0.5*d_new
      fp(k,imskh) = 0.0
      fp(k,itbrt) = 0.0
      fp(k,itbpe) = 0.0
      call update_particle_mass_fields(fp,k)
!
    endsubroutine apply_pe_breakup
!***********************************************************************
    subroutine insert_child_particle(fp,ineargrid,kmother,dp_child,np_child)
!
!  Append a child parcel to fp that inherits the mother state except for
!  droplet size, child swarm count, and breakup accumulator.
!  The child starts at the same position and with the same velocity as the
!  mother; only the parcel-internal droplet properties are changed.
!
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
      integer, intent(in) :: kmother
      real, intent(in) :: dp_child, np_child
!
      integer :: knew
!
      if (npar_loc >= mpar_loc) then
        call fatal_error('insert_child_particle', &
            'mpar_loc exhausted while creating KH breakup child parcels')
      endif
!
      knew = npar_loc + 1
      fp(knew,:) = fp(kmother,:)
      ineargrid(knew,:) = ineargrid(kmother,:)
!
!  Assign a fresh particle id so the new parcel can be tracked independently.
!
      ipar(knew) = max(maxval(ipar(1:npar_loc)), npar_inserted_tot) + 1
      npar_inserted_tot = ipar(knew)
!
      fp(knew,iap) = 0.5*dp_child
      fp(knew,inpswarm) = np_child
      fp(knew,imskh) = 0.0
      fp(knew,itbrt) = 0.0
      fp(knew,itbpe) = 0.0
!
      if (iapinit /= 0) fp(knew,iapinit) = fp(knew,iap)
      call update_particle_mass_fields(fp,knew,reset_initial=.true.)
!
      npar_loc = knew
!
    endsubroutine insert_child_particle
!***********************************************************************
    subroutine update_particle_mass_fields(fp,k,reset_initial)
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, intent(in) :: k
      logical, optional, intent(in) :: reset_initial
!
      logical :: lreset_initial
      real :: mp
!
      lreset_initial = .false.
      if (present(reset_initial)) lreset_initial = reset_initial
!
!  Only touch explicit mass fields when the particles_mass module is active.
!  Radius and number are always updated directly by this breakup module.
      if (.not. lparticles_mass) return
!
      mp = four_pi_over_three * rhopmat * fp(k,iap)**3
      if (imp /= 0) fp(k,imp) = mp
      if (irhosurf /= 0) fp(k,irhosurf) = rhopmat
      if (lreset_initial) then
        if (impinit /= 0) fp(k,impinit) = mp
        if (iapinit /= 0) fp(k,iapinit) = fp(k,iap)
      endif
!
    endsubroutine update_particle_mass_fields
!***********************************************************************
    real function mass_from_diameter(dp)
      real, intent(in) :: dp
!  Material mass represented by one physical droplet of diameter dp.
      mass_from_diameter = (pi / 6.0) * rhopmat * dp**3
    endfunction mass_from_diameter
!***********************************************************************
    real function diameter_from_mass(mp)
      real, intent(in) :: mp
!  Inverse of mass_from_diameter, with a safe zero-mass branch.
      if (mp <= 0.0) then
        diameter_from_mass = 0.0
      else
        diameter_from_mass = (6.0*mp/(pi*rhopmat))**(1.0/3.0)
      endif
    endfunction diameter_from_mass
!***********************************************************************
    function pe_regime_name(regime) result(name)
!  Short human-readable name for the Pilch-Erdman regime index used by
!  the lprint_breakup_mode diagnostic.
      integer, intent(in) :: regime
      character(len=18) :: name
!
      select case (regime)
      case (1)
        name = 'bag              '
      case (2)
        name = 'bag-stamen       '
      case (3)
        name = 'sheet stripping  '
      case (4)
        name = 'catastrophic     '
      case default
        name = 'vibrational/none '
      end select
!
    endfunction pe_regime_name
!***********************************************************************
endmodule Particles_breakup
