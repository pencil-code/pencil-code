! $Id$
!
!  Kelvin-Helmholtz / Rayleigh-Taylor breakup for Lagrangian particles
!  following Aguerre & Nigro, Computers and Fluids 190 (2019) 30-48,
!  sections 3.2-3.3 and Appendices B3-B4.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! CPARAM logical, parameter :: lparticles_breakup = .true.
!
! MPVAR CONTRIBUTION 0
! MPAUX CONTRIBUTION 2
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
  use Particles_sub, only: append_npaux, sum_par_name
!
  implicit none
!
  include 'particles_breakup.h'
!
  logical :: lparticles_breakup_kh = .false.
  logical :: lparticles_breakup_rt = .false.
  real :: tstart_particles_breakup = 0.0
  real :: kh_B0 = 0.61
  real :: kh_B1 = 40.0
  real :: kh_B2 = 0.05
  real :: rt_CRT = 0.1
  real :: rt_Ctau = 1.0
  real :: kh_sigma = 7.2e-2
  real :: kh_min_child_npswarm = 0.0
  integer :: rt_method = 2
  integer :: idiag_imskhm = 0
!
  namelist /particles_breakup_init_pars/ &
      lparticles_breakup_kh, lparticles_breakup_rt, tstart_particles_breakup, &
      kh_B0, kh_B1, kh_B2, rt_CRT, rt_Ctau, rt_method, kh_sigma, &
      kh_min_child_npswarm
!
  namelist /particles_breakup_run_pars/ &
      lparticles_breakup_kh, lparticles_breakup_rt, tstart_particles_breakup, &
      kh_B0, kh_B1, kh_B2, rt_CRT, rt_Ctau, rt_method, kh_sigma, &
      kh_min_child_npswarm
!
contains
!***********************************************************************
    subroutine register_particles_breakup()
!
!  Register particle-local breakup state.
!
      call append_npaux('imskh', imskh)
      call append_npaux('itbrt', itbrt)
!
    endsubroutine register_particles_breakup
!***********************************************************************
    subroutine initialize_particles_breakup(f)
!
!  Check the breakup configuration.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      if (lparticles_breakup_kh .or. lparticles_breakup_rt) then
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
      if (lparticles_breakup_kh .or. lparticles_breakup_rt) then
        if (kh_sigma <= 0.0) call fatal_error('initialize_particles_breakup', &
            'kh_sigma must be positive')
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
      if (lreset) idiag_imskhm = 0
!
      if (lroot .and. ip < 14) print*,'rprint_particles_breakup: run through parse list'
      do iname = 1,nname
        call parse_name(iname,cname(iname),cform(iname),'imskhm',idiag_imskhm)
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
      if (ldiagnos) then
        if (idiag_imskhm /= 0) call sum_par_name(fp(1:npar_loc,imskh),idiag_imskhm)
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
!  Apply KH/RT breakup after particle advection for the original particles only.
!  Newly created child parcels are appended to fp and take part from the next step on.
!  The implementation mirrors the paper logic:
!    1. compute local gas-particle slip quantities,
!    2. give RT priority whenever its instability criterion is active,
!    3. otherwise apply the KH relaxation / stripping update,
!    4. newly created KH child parcels take part from the next step on.
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
      real :: mp_old, mp_child, ms_new, ms_inc, np_child
      integer :: k, npar_loc_old
!
      if (.not. (lparticles_breakup_kh .or. lparticles_breakup_rt)) return
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
!  Work in terms of droplet diameter because the breakup correlations in the
!  paper are written for d_p, not for the stored parcel radius a_p.
        dp_old = 2.0*fp(k,iap)
        if (dp_old <= 0.0) cycle
!
!  Section 3.3 / Appendix B4: once RT is active, it takes priority over KH.
        if (lparticles_breakup_rt) then
          call compute_rt_scales(dp_old,rho_gas,nu_gas,urmag,cd_drag,acc_rt,lambda_rt,omega_rt,tau_rt)
          if (tau_rt > 0.0 .and. lambda_rt > 0.0 .and. dp_old > lambda_rt) then
            fp(k,itbrt) = fp(k,itbrt) + dt
            if (fp(k,itbrt) > tau_rt) call apply_rt_breakup(fp,k,dp_old,lambda_rt)
            cycle
          else
            fp(k,itbrt) = 0.0
          endif
        endif
!
        if (.not. lparticles_breakup_kh) cycle
!
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
      real :: denom
!
!  Dimensionless groups and KH correlations used in the breakup model.
!  These are the quantities needed to reconstruct the stable KH scale d_KH
!  and its response time tau_KH.
      we_g = rho_gas * urmag**2 * dp / kh_sigma
      we_p = rhopmat * urmag**2 * dp / kh_sigma
      rep = urmag * dp / max(nu_gas,tini)
!
      oh = sqrt(max(we_p,0.0)) / max(rep,tini)
      tay = oh * sqrt(max(we_g,0.0))
!      
      lambda_kh = 4.51 * dp * (1.0 + 0.45*sqrt(max(oh,0.0))) * (1.0 + 0.4*tay**0.7) / &
           (1.0 + 0.865*we_g**1.67)**0.6
!
      denom = (1.0 + oh) * (1.0 + 1.4*tay**0.6)
      omega_kh = (0.34 + 0.385*we_g**1.5) / max(denom,tini) * &
          sqrt(8.0*kh_sigma / max(rhopmat*dp**3,tini))
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
      lambda_rt = 2.0*pi*rt_CRT * sqrt(3.0*kh_sigma / max(acc_rt*density_jump,tini))
      omega_rt = sqrt(2.0*(acc_rt*density_jump)**1.5 / max(3.0*sqrt(3.0)*kh_sigma,tini))
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
      call update_particle_mass_fields(fp,k)
!
    endsubroutine apply_rt_breakup
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
endmodule Particles_breakup
