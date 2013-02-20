! $Id$
!
!  This module takes care of everything related to inertial particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 6
! MAUX CONTRIBUTION 2
! CPARAM logical, parameter :: lparticles=.true.
!
! PENCILS PROVIDED np; rhop
! PENCILS PROVIDED epsp; grhop(3)
!
!***************************************************************
module Particles
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_radius
!
  implicit none
!
  include 'particles.h'
!
  complex, dimension (7) :: coeff=(0.0,0.0)
  real, target, dimension (npar_species) :: tausp_species=0.0
  real, target, dimension (npar_species) :: tausp1_species=0.0
  real, dimension (3) :: temp_grad0=(/0.0,0.0,0.0/)
  real :: xp0=0.0, yp0=0.0, zp0=0.0, vpx0=0.0, vpy0=0.0, vpz0=0.0
  real :: xp1=0.0, yp1=0.0, zp1=0.0, vpx1=0.0, vpy1=0.0, vpz1=0.0
  real :: xp2=0.0, yp2=0.0, zp2=0.0, vpx2=0.0, vpy2=0.0, vpz2=0.0
  real :: xp3=0.0, yp3=0.0, zp3=0.0, vpx3=0.0, vpy3=0.0, vpz3=0.0
  real :: Lx0=0.0, Ly0=0.0, Lz0=0.0
  real :: delta_vp0=1.0, tausp=0.0, tausp1=0.0, eps_dtog=0.0
  real :: nu_epicycle=0.0, nu_epicycle2=0.0
  real :: beta_dPdr_dust=0.0, beta_dPdr_dust_scaled=0.0
  real :: tausg_min=0.0, tausg1_max=0.0, epsp_friction_increase=100.0
  real :: cdtp=0.2, cdtpgrav=0.1
  real :: gravx=0.0, gravz=0.0, gravr=1.0, kx_gg=1.0, kz_gg=1.0
  real :: gravsmooth=0.0, gravsmooth2=0.0, Ri0=0.25, eps1=0.5
  real :: kx_xxp=0.0, ky_xxp=0.0, kz_xxp=0.0, amplxxp=0.0
  real :: kx_vvp=0.0, ky_vvp=0.0, kz_vvp=0.0, amplvvp=0.0
  real :: kx_vpx=0.0, kx_vpy=0.0, kx_vpz=0.0
  real :: ky_vpx=0.0, ky_vpy=0.0, ky_vpz=0.0
  real :: kz_vpx=0.0, kz_vpy=0.0, kz_vpz=0.0
  real :: phase_vpx=0.0, phase_vpy=0.0, phase_vpz=0.0
  real :: tstart_dragforce_par=0.0
  real :: tstart_grav_par=0.0, tstart_grav_x_par=0.0
  real :: tstart_grav_z_par=0.0, tstart_grav_r_par=0.0
  real :: tstart_liftforce_par=0.0
  real :: tstart_brownian_par=0.0
  real :: tstart_collisional_cooling=0.0
  real :: tstart_sink_par=0.0
  real :: tau_coll_min=0.0, tau_coll1_max=0.0
  real :: coeff_restitution=0.5, mean_free_path_gas=0.0
  real :: pdlaw=0.0, rad_sphere=0.0
  real :: taucool=0.0, taucool1=0.0, brownian_T0=0.0, thermophoretic_T0=0.0
  real :: xsinkpoint=0.0, ysinkpoint=0.0, zsinkpoint=0.0, rsinkpoint=0.0
  real :: particles_insert_rate=0.
  real :: avg_n_insert, remaining_particles=0.0
  real :: max_particle_insert_time=huge1
  real :: Deltauy_gas_friction=0.0
  real :: cond_ratio=0.0
  real :: fake_particle_radius=0.0
  integer :: l_hole=0, m_hole=0, n_hole=0
  integer :: iffg=0, ifgx=0, ifgy=0, ifgz=0
  logical :: ldragforce_dust_par=.false., ldragforce_gas_par=.false.
  logical :: ldragforce_stiff=.false., ldragforce_radialonly=.false.
  logical :: ldragforce_heat=.false., lcollisional_heat=.false.
  logical :: lpar_spec=.false., lcompensate_friction_increase=.false.
  logical :: lcollisional_cooling_taucool=.false.
  logical :: lcollisional_cooling_rms=.false.
  logical :: lcollisional_cooling_twobody=.false.
  logical :: lcollisional_dragforce_cooling=.false.
  logical :: ltau_coll_min_courant=.true.
  logical :: ldragforce_equi_global_eps=.false., ldragforce_equi_noback=.false.
  logical :: ldraglaw_epstein=.true., ldraglaw_epstein_stokes_linear=.false.
  logical :: ldraglaw_simple=.false.
  logical :: ldraglaw_steadystate=.false., ldraglaw_variable=.false.
  logical :: ldraglaw_epstein_transonic=.false.
  logical :: ldraglaw_eps_stk_transonic=.false.
  logical :: ldraglaw_variable_density=.false.
  logical :: lcoldstart_amplitude_correction=.false.
  logical :: luse_tau_ap=.true.
  logical :: lbrownian_forces=.false.
  logical :: lthermophoretic_forces=.false.
  logical :: lenforce_policy=.false., lnostore_uu=.true.
  logical :: ldt_grav_par=.true., ldt_adv_par=.true.
  logical :: lsinkpoint=.false., lglobalrandom=.false.
  logical :: lcoriolis_force_par=.true., lcentrifugal_force_par=.false.
  logical :: lcalc_uup=.false.
  logical :: lparticle_gravity=.true.
  logical :: lcylindrical_gravity_par=.false.
!
  character (len=labellen) :: interp_pol_uu ='ngp'
  character (len=labellen) :: interp_pol_oo ='ngp'
  character (len=labellen) :: interp_pol_TT ='ngp'
  character (len=labellen) :: interp_pol_rho='ngp'
  character (len=labellen) :: interp_pol_gradTT='ngp'
!
  character (len=labellen), dimension (ninit) :: initxxp='nothing'
  character (len=labellen), dimension (ninit) :: initvvp='nothing'
  character (len=labellen) :: gravx_profile='', gravz_profile=''
  character (len=labellen) :: gravr_profile=''
  character (len=labellen) :: thermophoretic_eq= 'nothing'
!
  integer :: init_repeat=0       !repeat particle initialization for distance statistics
!
!  Interactions with special/shell
!
  integer :: k_shell=-1            !k associated with minshell (special/shell.f90)
  logical :: l_shell=.false.       !using special/shell.f90 for gas velocities
!
  real, dimension(3) :: uup_shared=0
  real :: turnover_shared=0
  logical :: vel_call=.false., turnover_call=.false.
!
  namelist /particles_init_pars/ &
      initxxp, initvvp, xp0, yp0, zp0, vpx0, vpy0, vpz0, delta_vp0, &
      ldragforce_gas_par, ldragforce_dust_par, bcpx, bcpy, bcpz, tausp, &
      beta_dPdr_dust, np_swarm, mp_swarm, mpmat, rhop_swarm, eps_dtog, &
      nu_epicycle, rp_int, rp_ext, gravx_profile, gravz_profile, &
      gravr_profile, gravx, gravz, gravr, gravsmooth, kx_gg, kz_gg, Ri0, &
      eps1, lmigration_redo, ldragforce_equi_global_eps, coeff, kx_vvp, &
      ky_vvp, kz_vvp, amplvvp, kx_xxp, ky_xxp, kz_xxp, amplxxp, kx_vpx, &
      kx_vpy, kx_vpz, ky_vpx, ky_vpy, ky_vpz, kz_vpx, kz_vpy, kz_vpz, &
      phase_vpx, phase_vpy, phase_vpz, lcoldstart_amplitude_correction, &
      lparticlemesh_cic, lparticlemesh_tsc, linterpolate_spline, &
      tstart_dragforce_par, tstart_grav_par, lparticle_gravity,&
      tstart_grav_x_par, tstart_grav_z_par,tstart_grav_r_par, taucool, &
      lcollisional_cooling_taucool, lcollisional_cooling_rms, &
      lcollisional_cooling_twobody, tausp_species, tau_coll_min, &
      ltau_coll_min_courant, coeff_restitution, tstart_collisional_cooling, &
      tausg_min, l_hole, m_hole, n_hole, &
      epsp_friction_increase,lcollisional_dragforce_cooling, ldragforce_heat, &
      lcollisional_heat, lcompensate_friction_increase, &
      lmigration_real_check, ldraglaw_epstein,ldraglaw_simple,ldraglaw_epstein_stokes_linear, &
      mean_free_path_gas, ldraglaw_epstein_transonic, lcheck_exact_frontier, &
      ldraglaw_eps_stk_transonic, pdlaw, rad_sphere, ldragforce_stiff, &
      ldraglaw_steadystate, tstart_liftforce_par, &
      tstart_brownian_par, tstart_sink_par, &
      lbrownian_forces,lthermophoretic_forces,lenforce_policy, &
      interp_pol_uu,interp_pol_oo,interp_pol_TT,interp_pol_rho, brownian_T0, &
      thermophoretic_T0, lnostore_uu, ldt_grav_par, ldragforce_radialonly, &
      lsinkpoint, xsinkpoint, ysinkpoint, zsinkpoint, rsinkpoint, &
      lcoriolis_force_par, lcentrifugal_force_par, ldt_adv_par, Lx0, Ly0, &
      Lz0, lglobalrandom, linsert_particles_continuously, &
      lrandom_particle_pencils, lnocalc_np, lnocalc_rhop, np_const, &
      rhop_const, ldragforce_equi_noback, rhopmat, Deltauy_gas_friction, xp1, &
      yp1, zp1, vpx1, vpy1, vpz1, xp2, yp2, zp2, vpx2, vpy2, vpz2, &
      xp3, yp3, zp3, vpx3, vpy3, vpz3, lsinkparticle_1, rsinkparticle_1, &
      lcalc_uup, temp_grad0, thermophoretic_eq, cond_ratio, interp_pol_gradTT
!
  namelist /particles_run_pars/ &
      bcpx, bcpy, bcpz, tausp, dsnap_par_minor, beta_dPdr_dust, &
      ldragforce_gas_par, ldragforce_dust_par, np_swarm, mp_swarm, &
      rhop_swarm, eps_dtog, cdtp, cdtpgrav, lpar_spec, linterp_reality_check, &
      nu_epicycle, gravx_profile, gravz_profile, gravr_profile, gravx, gravz, &
      gravr, gravsmooth, kx_gg, kz_gg, lmigration_redo, tstart_dragforce_par, &
      tstart_grav_par, tstart_grav_x_par, &
      tstart_grav_z_par, tstart_grav_r_par, lparticle_gravity, &
      lparticlemesh_cic, lparticlemesh_tsc, taucool, &
      lcollisional_cooling_taucool, lcollisional_cooling_rms, &
      lcollisional_cooling_twobody, lcollisional_dragforce_cooling, &
      tau_coll_min, ltau_coll_min_courant, coeff_restitution, &
      tstart_collisional_cooling, tausg_min, epsp_friction_increase, &
      ldragforce_heat, lcollisional_heat, lcompensate_friction_increase, &
      lmigration_real_check,ldraglaw_variable, luse_tau_ap, ldraglaw_epstein, ldraglaw_simple,&
      ldraglaw_epstein_stokes_linear, mean_free_path_gas, &
      ldraglaw_epstein_transonic, lcheck_exact_frontier, &
      ldraglaw_eps_stk_transonic, ldragforce_stiff, &
      ldraglaw_variable_density, ldraglaw_steadystate, tstart_liftforce_par, &
      tstart_brownian_par, tstart_sink_par, &
      lbrownian_forces, lenforce_policy, &
      interp_pol_uu,interp_pol_oo, interp_pol_TT, interp_pol_rho, &
      brownian_T0,thermophoretic_T0, lnostore_uu, ldt_grav_par, &
      ldragforce_radialonly, lsinkpoint, xsinkpoint, ysinkpoint, zsinkpoint, &
      rsinkpoint, lcoriolis_force_par, lcentrifugal_force_par, ldt_adv_par, &
      linsert_particles_continuously, particles_insert_rate, &
      max_particle_insert_time, lrandom_particle_pencils, lnocalc_np, &
      lnocalc_rhop, np_const, rhop_const, Deltauy_gas_friction, &
      loutput_psize_dist, log_ap_min_dist, log_ap_max_dist, nbin_ap_dist, &
      lsinkparticle_1, rsinkparticle_1, lthermophoretic_forces, temp_grad0, &
      thermophoretic_eq, cond_ratio, interp_pol_gradTT, lcommunicate_rhop, & 
      lcommunicate_np, lcylindrical_gravity_par
!
  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0
  integer :: idiag_xp2m=0, idiag_yp2m=0, idiag_zp2m=0
  integer :: idiag_rpm=0, idiag_rp2m=0
  integer :: idiag_vpxm=0, idiag_vpym=0, idiag_vpzm=0
  integer :: idiag_vpx2m=0, idiag_vpy2m=0, idiag_vpz2m=0, idiag_ekinp=0
  integer :: idiag_vpxmax=0, idiag_vpymax=0, idiag_vpzmax=0, idiag_vpmax=0
  integer :: idiag_vpxvpym=0, idiag_vpxvpzm=0, idiag_vpyvpzm=0
  integer :: idiag_rhopvpxm=0, idiag_rhopvpym=0, idiag_rhopvpzm=0
  integer :: idiag_rhopvpxt=0, idiag_rhopvpyt=0, idiag_rhopvpzt=0
  integer :: idiag_rhopvpysm=0
  integer :: idiag_lpxm=0, idiag_lpym=0, idiag_lpzm=0
  integer :: idiag_lpx2m=0, idiag_lpy2m=0, idiag_lpz2m=0
  integer :: idiag_npm=0, idiag_np2m=0, idiag_npmax=0, idiag_npmin=0
  integer :: idiag_dtdragp=0
  integer :: idiag_nparmin=0, idiag_nparmax=0, idiag_npargone=0
  integer :: idiag_rhopm=0, idiag_rhoprms=0, idiag_rhop2m=0, idiag_rhopmax=0
  integer :: idiag_rhopmin=0, idiag_decollp=0, idiag_rhopmphi=0
  integer :: idiag_epspmin=0, idiag_epspmax=0
  integer :: idiag_npmx=0, idiag_npmy=0, idiag_npmz=0
  integer :: idiag_rhopmx=0, idiag_rhopmy=0, idiag_rhopmz=0
  integer :: idiag_epspmx=0, idiag_epspmy=0, idiag_epspmz=0
  integer :: idiag_mpt=0, idiag_dedragp=0, idiag_rhopmxy=0, idiag_rhopmr=0
  integer :: idiag_dvpx2m=0, idiag_dvpy2m=0, idiag_dvpz2m=0
  integer :: idiag_dvpm=0, idiag_dvpmax=0, idiag_epotpm=0
  integer :: idiag_rhopmxz=0, idiag_nparpmax=0, idiag_npmxy=0
  integer :: idiag_eccpxm=0, idiag_eccpym=0, idiag_eccpzm=0
  integer :: idiag_eccpx2m=0, idiag_eccpy2m=0, idiag_eccpz2m=0
  integer :: idiag_vprms=0, idiag_vpyfull2m=0, idiag_deshearbcsm=0
!
  contains
!***********************************************************************
    subroutine register_particles()
!
!  Set up indices for access to the fp and dfp arrays
!
!  29-dec-04/anders: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Indices for particle position.
!
      ixp=npvar+1
      pvarname(npvar+1)='ixp'
      iyp=npvar+2
      pvarname(npvar+2)='iyp'
      izp=npvar+3
      pvarname(npvar+3)='izp'
!
!  Indices for particle velocity.
!
      ivpx=npvar+4
      pvarname(npvar+4)='ivpx'
      ivpy=npvar+5
      pvarname(npvar+5)='ivpy'
      ivpz=npvar+6
      pvarname(npvar+6)='ivpz'
!
!  Increase npvar accordingly.
!
      npvar=npvar+6
!
!  Set indices for particle assignment.
!
      if (.not. lnocalc_np) call farray_register_auxiliary('np',inp)
      if (.not. lnocalc_rhop) call farray_register_auxiliary('rhop',irhop, &
          communicated=lparticles_sink)
      if (lcalc_uup .or. ldragforce_stiff) then
        call farray_register_auxiliary('uup',iuup,vector=3)
        iupx=iuup; iupy=iuup+1; iupz=iuup+2
      endif
!
!  Special variable for stiff drag force equations.
!
      if (ldragforce_stiff) then
        call farray_register_auxiliary('ffg',iffg,communicated=.true.,vector=3)
        ifgx=iffg; ifgy=iffg+1; ifgz=iffg+2
      endif
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_particles','npvar > mpvar')
      endif
!
    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: rho0, cs0
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      real :: rhom
      integer :: ierr, jspec
!
!  This module is incompatible with particle block domain decomposition.
!
      if (lparticles_blocks) then
        if (lroot) then
          print*, 'initialize_particles: must use PARTICLES =  PARTICLES_DUST_BLOCKS'
          print*, '                      with particle block domain decomposition'
        endif
        call fatal_error('initialize_particles','')
      endif
!
!  The inverse stopping time is needed for drag force and collisional cooling.
!
      if (tausp/=0.0) tausp1=1/tausp
!
!  Inverse cooling time.
!
      if (taucool/=0.0) taucool1=1/taucool
!
!  Inverse material density.
!
      if (rhopmat/=0.0) rhopmat1=1/rhopmat
!
!  Multiple dust species. Friction time is given in the array tausp_species.
!
      if (npar_species>1) then
        if (lroot) then
          print*, 'initialize_particles: '// &
              'Number of particle species = ', npar_species
          print*, 'initialize_particles: tausp_species = ', tausp_species
        endif
!
!  Must have set tausp_species for drag force.
!
        if (ldragforce_dust_par .or. ldragforce_gas_par) then
          if (any(tausp_species==0.0)) then
            if (lroot) print*, &
                'initialize_particles: drag force must have tausp_species/=0.0!'
                call fatal_error('initialize_particles','')
          endif
!
!  Inverse friction time is needed for drag force.
!
          do jspec=1,npar_species
            if (tausp_species(jspec)/=0.0) &
                tausp1_species(jspec)=1/tausp_species(jspec)
          enddo
        endif
      else
!
!  Single dust species => If tausp_species is set, it is probably an error.
!
        if (any(tausp_species/=0.0) .and. tausp/=tausp_species(1)) then
          call fatal_error('initialize_particles', 'When there is only '// &
              '1 particle species, use tausp instead of tausp_species')
        endif
        tausp_species(1)=tausp
        if (tausp_species(1)/=0.0) tausp1_species(1)=1/tausp_species(1)
      endif
!
!  Share friction time (but only if Epstein drag regime!).
!
      if (ldraglaw_epstein .or. ldraglaw_simple) then
        call put_shared_variable( 'tausp_species', tausp_species)
        call put_shared_variable('tausp1_species',tausp1_species)
      endif
!
!  Global gas pressure gradient seen from the perspective of the dust.
!
      if (beta_dPdr_dust/=0.0) then
        beta_dPdr_dust_scaled=beta_dPdr_dust*Omega/cs0
        if (lroot) print*, 'initialize_particles: Global pressure '// &
            'gradient with beta_dPdr_dust=', beta_dPdr_dust
      endif
!
!  Calculate mass density per particle (for back-reaction drag force on gas)
!  based on the dust-to-gas ratio eps:
!
!    rhop_swarm*N_cell = eps*rhom
!
!  where rhop_swarm is the mass density per superparticle, N_cell is the number
!  of particles per grid cell and rhom is the mean gas density in the box.
!
      if (eps_dtog>0.0) then
! For stratification, take into account gas present outside the simulation box.
        if ((lgravz .and. lgravz_gas) .or. gravz_profile=='linear' ) then
          ! rhom = (total mass) / (box volume) = Sigma / Lz
          ! Sigma = sqrt(2pi) * rho0 * H
          !   rho0 = mid-plane density, H = (sound speed) / (epicycle frequency)
          rhom = sqrt(2.0 * pi) / Lz
          if (nu_epicycle > 0.0) rhom = rhom * (rho0 * cs0 / nu_epicycle)
        else
          rhom = rho0
        endif
        if (rhop_swarm==0.0) &
            rhop_swarm = eps_dtog*rhom/(real(npar)/nwgrid)
        if (mp_swarm==0.0) &
            mp_swarm   = eps_dtog*rhom*box_volume/(real(npar))
        if (lroot) print*, 'initialize_particles: '// &
            'dust-to-gas ratio eps_dtog=', eps_dtog
      endif
!
      if (lroot) then
        print*, 'initialize_particles: '// &
            'mass per constituent particle mpmat=', mpmat
        print*, 'initialize_particles: '// &
            'mass per superparticle mp_swarm =', mp_swarm
        print*, 'initialize_particles: '// &
            'number density per superparticle np_swarm=', np_swarm
        print*, 'initialize_particles: '// &
            'mass density per superparticle rhop_swarm=', rhop_swarm
      endif
!
!  Calculate nu_epicycle**2 for gravity.
!
      if (gravz_profile=='' .and. nu_epicycle/=0.0) gravz_profile='linear'
      nu_epicycle2=nu_epicycle**2
!
!  Calculate gravsmooth**2 for gravity.
!
      if (gravsmooth/=0.0) gravsmooth2=gravsmooth**2
!
!  Share Keplerian gravity.
!
      call put_shared_variable('gravr',gravr,ierr)
      if (ierr/=0) call fatal_error('initialize_particles', &
          'there was a problem when sharing gravr')
!
!  Inverse of minimum gas friction time (time-step control).
!
      if (tausg_min/=0.0) tausg1_max=1.0/tausg_min
!
!  Set minimum collisional time-scale so that time-step is not affected.
!
      if ((.not. lstarting) .and. ltau_coll_min_courant) then
        if (cs0==impossible) then
          tau_coll_min=impossible
        else
          tau_coll_min=2*dx/cs0
        endif
        if (lroot) print*, 'initialize particles: set minimum collisional '// &
            'time-scale equal to two times the Courant time-step.'
      endif
!
!  Inverse of minimum collisional time-scale.
!
      if ((.not. lstarting) .and. tau_coll_min>0.0) tau_coll1_max=1/tau_coll_min
!
!  Gas density is needed for back-reaction friction force.
!
      if (ldragforce_gas_par .and. .not. ldensity) then
        if (lroot) then
          print*, 'initialize_particles: friction force on gas only works '
          print*, '                      together with gas density module!'
        endif
        call fatal_error('initialize_particles','')
      endif
!
!  Need to map particles on the grid for dragforce on gas.
!
      if (ldragforce_gas_par) then
!
!  When drag force is smoothed, df is also set in the first ghost zone. This
!  region needs to be folded back into the df array after pde is finished,
!
        if (lparticlemesh_cic .or. lparticlemesh_tsc) lfold_df=.true.
      endif
!
      if (lcollisional_cooling_twobody) then
        allocate(kneighbour(mpar_loc))
        lshepherd_neighbour=.true.
      endif
!
      if (ldraglaw_epstein_stokes_linear) ldraglaw_epstein=.false.
      if (ldraglaw_epstein_transonic    .or.&
          ldraglaw_eps_stk_transonic    .or.&
          ldraglaw_steadystate          .or.&
          ldraglaw_simple) then
        ldraglaw_epstein=.false.
      endif
      if (ldraglaw_epstein_transonic         .and.&
          ldraglaw_eps_stk_transonic) then
        print*,'both epstein and epstein-stokes transonic '//&
               'drag laws are switched on. You cannot have '//&
               'both. Stop and choose only one.'
        call fatal_error('initialize_particles','')
      endif
!
!  Stiff drag force approximation.
!
      if (ldragforce_stiff) then
        if (ldragforce_dust_par .or. ldragforce_gas_par) then
          if (lroot) print*, 'initialize_particles: stiff drag force '// &
              'approximation is incompatible with normal drag'
          call fatal_error('initialize_particles','')
        endif
        f(l1:l2,m1:m2,n1:n2,ifgx:ifgz)=0.0
      endif
!
!  Initialize storage of energy gain released by shearing boundaries.
!
      if (idiag_deshearbcsm/=0) energy_gain_shear_bcs=0.0
!
!  Drag force on gas right now assumed rhop_swarm is the same for all particles.
!
      if (ldragforce_gas_par.and.(lparticles_radius.or.lparticles_number)) then
        if (lroot) print*, 'initialize_particles: drag force on gas is '// &
            'not yet implemented for variable particle radius or number'
        call fatal_error('initialize_particles','')
      endif
!
!  Fatal error if sink particle radius is zero or negative.
!
      if (lsinkparticle_1 .and. rsinkparticle_1<=0.0) then
        if (lroot) print*, 'initialize_particles: sink particle radius is '// &
            'zero or negative: ', rsinkparticle_1
        call fatal_error('initialize_particles','')
      endif
!
!  Particle self gravity for x,z,r direction
!
      if ( tstart_grav_par > 0.0) then
        tstart_grav_x_par = tstart_grav_par
        tstart_grav_z_par = tstart_grav_par
        tstart_grav_r_par = tstart_grav_par
      endif
!
!  Set up interpolation logicals. These logicals can be OR'ed with some logical
!  in the other particle modules' initialization subroutines to enable
!  interpolation based on some condition local to that module.
!  (The particles_spin module will for instance enable interpolation of the
!  vorticity oo)
!
      if (lnostore_uu) then
        interp%luu=.false.
      else
        interp%luu=ldragforce_dust_par.or.ldraglaw_steadystate.or. &
          lparticles_spin
      endif
      interp%loo=.false.
      interp%lTT=(lbrownian_forces.and.(brownian_T0==0.0))&
          .or.(lthermophoretic_forces.and.(thermophoretic_T0==0.0))
      interp%lgradTT=lthermophoretic_forces .and. (temp_grad0(1)==0.0) &
          .and. (temp_grad0(2)==0.0) .and. (temp_grad0(3)==0.0)
      interp%lrho=lbrownian_forces.or.ldraglaw_steadystate &
          .or. lthermophoretic_forces
!
!  Determine interpolation policies:
!   Make sure that interpolation of uu is chosen in a backwards compatible
!   manner. NGP is chosen by default.
!
      if (.not.lenforce_policy) then
        if (lparticlemesh_cic) then
          interp_pol_uu='cic'
        else if (lparticlemesh_tsc) then
          interp_pol_uu='tsc'
        endif
      endif
!
!  Overwrite with new policy variables:
!
      select case (interp_pol_uu)
      case ('tsc')
        interp%pol_uu=tsc
      case ('cic')
        interp%pol_uu=cic
      case ('ngp')
        interp%pol_uu=ngp
      case default
        call fatal_error('initialize_particles','No such such value for '// &
          'interp_pol_uu: '//trim(interp_pol_uu))
      endselect
!
      select case (interp_pol_oo)
      case ('tsc')
        interp%pol_oo=tsc
      case ('cic')
        interp%pol_oo=cic
      case ('ngp')
        interp%pol_oo=ngp
      case default
        call fatal_error('initialize_particles','No such such value for '// &
          'interp_pol_oo: '//trim(interp_pol_oo))
      endselect
!
      select case (interp_pol_TT)
      case ('tsc')
        interp%pol_TT=tsc
      case ('cic')
        interp%pol_TT=cic
      case ('ngp')
        interp%pol_TT=ngp
      case default
        call fatal_error('initialize_particles','No such such value for '// &
          'interp_pol_TT: '//trim(interp_pol_TT))
      endselect
!
      select case (interp_pol_gradTT)
      case ('tsc')
        call fatal_error('initialize_particles','Not implemented gradTT'// &
          'interp_pol_gradTT: '//trim(interp_pol_gradTT))
      case ('cic')
        call fatal_error('initialize_particles','Not implemented gradTT'// &
          'interp_pol_gradTT: '//trim(interp_pol_gradTT))
      case ('ngp')
        interp%pol_gradTT=ngp
      case default
        call fatal_error('initialize_particles','No such such value for '// &
          'interp_pol_gradTT: '//trim(interp_pol_gradTT))
      endselect
!
      select case (interp_pol_rho)
      case ('tsc')
        interp%pol_rho=tsc
      case ('cic')
        interp%pol_rho=cic
      case ('ngp')
        interp%pol_rho=ngp
      case default
        call fatal_error('initialize_particles','No such such value for '// &
          'interp_pol_rho: '//trim(interp_pol_rho))
      endselect
!
      if (l_shell) then
        if ( k_shell < 0) call fatal_error('initialize_particles','Set k_shell')
        call put_shared_variable('uup_shared',uup_shared,ierr)
        call put_shared_variable('vel_call',vel_call,ierr)
        call put_shared_variable('turnover_call',turnover_call,ierr)
        call put_shared_variable('turnover_shared',turnover_shared,ierr)
      endif
!
!  Write constants to disk.
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position='append')
          write (1,*) 'np_swarm=', np_swarm
          write (1,*) 'mpmat=', mpmat
          write (1,*) 'mp_swarm=', mp_swarm
          write (1,*) 'rhop_swarm=', rhop_swarm
        close (1)
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f,fp,ineargrid)
!
!  Initial positions and velocities of dust particles.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: beta_glnrho_global, cs20
      use General, only: random_number_wrapper
      use Mpicomm, only: mpireduce_sum, mpibcast_real
      use InitialCondition, only: initial_condition_xxp,&
                                  initial_condition_vvp
      use Particles_diagnos_dv, only: repeated_init
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (3) :: uup, Lxyz_par, xyz0_par, xyz1_par
      real :: vpx_sum, vpy_sum, vpz_sum
      real :: r, p, q, px, py, pz, eps, cs, k2_xxp, rp2
      real :: dim1, npar_loc_x, npar_loc_y, npar_loc_z, dx_par, dy_par, dz_par
      real :: rad,rad_scl,phi,tmp,OO,xx0,yy0,r2
      integer :: l, j, k, ix0, iy0, iz0
      logical :: lequidistant=.false.
!
      intent (out) :: f, fp, ineargrid
!
!  Use either a local random position or a global random position for certain
!  initial conditions. The default is a local random position, but the equal
!  number of particles per processors means that this is not completely random.
!
      if (lglobalrandom) then
        Lxyz_par=Lxyz
        xyz0_par=xyz0
        xyz1_par=xyz1
      else
        Lxyz_par=Lxyz_loc
        xyz0_par=xyz0_loc
        xyz1_par=xyz1_loc
      endif
!
!  Initial particle position.
!
      do j=1,ninit
!
        select case (initxxp(j))
!
        case ('nothing')
          if (lroot .and. j==1) print*, 'init_particles: nothing'
!
        case ('origin')
          if (lroot) print*, 'init_particles: All particles at origin'
          fp(1:npar_loc,ixp:izp)=0.0
!
        case ('zero-z')
          if (lroot) print*, 'init_particles: Zero z coordinate'
          fp(1:npar_loc,izp)=0.0
!
        case ('constant')
          if (lroot) &
              print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
          fp(1:npar_loc,ixp)=xp0
          fp(1:npar_loc,iyp)=yp0
          fp(1:npar_loc,izp)=zp0
!
        case ('constant-1')
          if (lroot) &
              print*, 'init_particles: Particle 1 at x,y,z=', xp1, yp1, zp1
          do k=1,npar_loc
            if (ipar(k)==1) then
              fp(k,ixp)=xp1
              fp(k,iyp)=yp1
              fp(k,izp)=zp1
            endif
          enddo
!
        case ('constant-2')
          if (lroot) &
              print*, 'init_particles: Particle 2 at x,y,z=', xp2, yp2, zp2
          do k=1,npar_loc
            if (ipar(k)==2) then
              fp(k,ixp)=xp2
              fp(k,iyp)=yp2
              fp(k,izp)=zp2
            endif
          enddo
!
        case ('constant-3')
          if (lroot) &
              print*, 'init_particles: Particle 2 at x,y,z=', xp3, yp3, zp3
          do k=1,npar_loc
            if (ipar(k)==3) then
              fp(k,ixp)=xp3
              fp(k,iyp)=yp3
              fp(k,izp)=zp3
            endif
          enddo
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
            if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
            if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
          enddo
          if (nxgrid/=1) &
              fp(1:npar_loc,ixp)=xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          if (nzgrid/=1) &
              fp(1:npar_loc,izp)=xyz0_par(3)+fp(1:npar_loc,izp)*Lxyz_par(3)
!
        case ('random-circle')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            call random_number_wrapper(r)
            if (zp0>yp0) then
              fp(k,ixp)=xp0*cos((zp0-yp0)*r+yp0)
              fp(k,iyp)=xp0*sin((zp0-yp0)*r+yp0)
            else
              fp(k,ixp)=xp0*cos(2*pi*r)
              fp(k,iyp)=xp0*sin(2*pi*r)
            endif
          enddo
!
        case ('random-sphere')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'in a sphere around (0,0,0) with radius=',rad_sphere
          if (rad_sphere==0) then
            call fatal_error('init_particles','random-sphere '// &
                  'radius needs to be larger than zero')
          endif
          if (-rad_sphere<xyz0(1) .or. rad_sphere>xyz1(1) .or. &
              -rad_sphere<xyz0(2) .or. rad_sphere>xyz1(2) .or. &
              -rad_sphere<xyz0(3) .or. rad_sphere>xyz1(3)) then
            call fatal_error('init_particles','random-sphere '// &
                 'sphere needs to fit in the box')
          endif
          if (lcartesian_coords) then
            do k=1,npar_loc
              rp2=2.*rad_sphere**2
              do while (rp2>rad_sphere**2)
                call random_number_wrapper(fp(k,ixp))
                call random_number_wrapper(fp(k,iyp))
                call random_number_wrapper(fp(k,izp))
                fp(k,ixp)=(fp(k,ixp)-0.5)*2.*rad_sphere
                fp(k,iyp)=(fp(k,iyp)-0.5)*2.*rad_sphere
                fp(k,izp)=(fp(k,izp)-0.5)*2.*rad_sphere
                rp2=fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2
              enddo
            enddo
          else
            call fatal_error('init_particles','random-sphere '// &
                 'only implemented for cartesian coordinates')
          endif
!
        case ('random-line-x')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
          enddo
          if (nxgrid/=1) &
              fp(1:npar_loc,ixp)=xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          fp(1:npar_loc,iyp)=yp0
          fp(1:npar_loc,izp)=zp0
!
        case ('random-line-y')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=1,npar_loc
            if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
          enddo
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          fp(1:npar_loc,ixp)=xp0
          fp(1:npar_loc,izp)=zp0
!
        case ('random-hole')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'with inner hole'
          do k=1,npar_loc
            rp2=-1.0
            do while (rp2<rp_int**2)
              if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
              if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
              if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
              if (nxgrid/=1) fp(k,ixp)=xyz0(1)+fp(k,ixp)*Lxyz(1)
              if (nygrid/=1) fp(k,iyp)=xyz0(2)+fp(k,iyp)*Lxyz(2)
              if (nzgrid/=1) fp(k,izp)=xyz0(3)+fp(k,izp)*Lxyz(3)
              rp2=fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2
            enddo
          enddo
!
        case ('random-box')
          if (lroot) print*, 'init_particles: Random particle positions '// &
               'within a box'
          do k=1,npar_loc
            if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
            if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
            if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
            if (lcylindrical_coords) then
              xx0=xp0+fp(k,ixp)*Lx0
              yy0=yp0+fp(k,iyp)*Ly0
              r2=xx0**2+yy0**2
              if (nxgrid/=1) fp(k,ixp)=sqrt(r2)
              if (nygrid/=1) fp(k,iyp)=atan(yy0/xx0)+pi*(xx0/abs(xx0)-1)*0.5
              if (nzgrid/=1) fp(k,izp)=zp0+fp(k,izp)*Lz0
            else
              if (nxgrid/=1) fp(k,ixp)=xp0+fp(k,ixp)*Lx0
              if (nygrid/=1) fp(k,iyp)=yp0+fp(k,iyp)*Ly0
              if (nzgrid/=1) fp(k,izp)=zp0+fp(k,izp)*Lz0
            endif
          enddo
!
       case ('random-cylindrical','random-cyl')
!
          if (lroot) print*, 'init_particles: Random particle '//&
               'cylindrical positions with power-law pdlaw=',pdlaw
!
          do k=1,npar_loc
!
! Start the particles obeying a power law pdlaw
!
            tmp=2-pdlaw
            call random_number_wrapper(rad_scl)
            rad_scl = rp_int**tmp + rad_scl*(rp_ext**tmp-rp_int**tmp)
            rad = rad_scl**(1./tmp)
!
! Random in azimuth
!
            call random_number_wrapper(phi)
!
            if (lcartesian_coords) then
              phi = 2*pi*phi
              if (nxgrid/=1) fp(k,ixp)=rad*cos(phi)
              if (nygrid/=1) fp(k,iyp)=rad*sin(phi)
            elseif (lcylindrical_coords) then
              phi = xyz0_par(2)+phi*Lxyz_par(2)
              if (nxgrid/=1) fp(k,ixp)=rad
              if (nygrid/=1) fp(k,iyp)=phi
            elseif (lspherical_coords) then
              call fatal_error('init_particles','random-cylindrical '// &
                  'not implemented for spherical coordinates')
            endif
!
            if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
            if (nzgrid/=1) &
                fp(k,izp)=xyz0_par(3)+fp(k,izp)*Lxyz_par(3)
!
          enddo
!
        case ('np-constant')
          if (lroot) print*, 'init_particles: Constant number density'
          k=1
k_loop:   do while (.not. (k>npar_loc))
            do l=l1,l2; do m=m1,m2; do n=n1,n2
              if (nxgrid/=1) call random_number_wrapper(px)
              if (nygrid/=1) call random_number_wrapper(py)
              if (nzgrid/=1) call random_number_wrapper(pz)
              fp(k,ixp)=x(l)+(px-0.5)*dx
              fp(k,iyp)=y(m)+(py-0.5)*dy
              fp(k,izp)=z(n)+(pz-0.5)*dz
              k=k+1
              if (k>npar_loc) exit k_loop
            enddo; enddo; enddo
          enddo k_loop
!
        case ('equidistant')
          if (lroot) print*, 'init_particles: Particles placed equidistantly'
          dim1=1.0/dimensionality
!
!  Number of particles per direction. Found by solving the equation system
!
!    npar_loc_x/npar_loc_y = Lx_loc/Ly_loc
!    npar_loc_x/npar_loc_z = Lx_loc/Lz_loc
!    npar_loc_y/npar_loc_z = Ly_loc/Lz_loc
!    npar_loc_x*npar_loc_y*npar_loc_z = npar_loc
!
!  Found it to be easier to separate in all possible dimensionalities.
!  For a missing direction i, set npar_loc_i=1 in the above equations and
!  ignore any equation that has Li_loc in it.
!
!  Initiate to avoid compiler warnings. Will be overwritten.
          npar_loc_x=1;npar_loc_y=1;npar_loc_z=1
!
          if (dimensionality==3) then
!  3-D
            npar_loc_x=(npar_loc*Lxyz_loc(1)**2/(Lxyz_loc(2)*Lxyz_loc(3)))**dim1
            npar_loc_y=(npar_loc*Lxyz_loc(2)**2/(Lxyz_loc(1)*Lxyz_loc(3)))**dim1
            npar_loc_z=(npar_loc*Lxyz_loc(3)**2/(Lxyz_loc(1)*Lxyz_loc(2)))**dim1
          elseif (dimensionality==2) then
!  2-D
            if (nxgrid==1) then
              npar_loc_x=1
              npar_loc_y=(npar_loc*Lxyz_loc(2)/Lxyz_loc(3))**dim1
              npar_loc_z=(npar_loc*Lxyz_loc(3)/Lxyz_loc(2))**dim1
            elseif (nygrid==1) then
              npar_loc_x=(npar_loc*Lxyz_loc(1)/Lxyz_loc(3))**dim1
              npar_loc_y=1
              npar_loc_z=(npar_loc*Lxyz_loc(3)/Lxyz_loc(2))**dim1
            elseif (nzgrid==1) then
              npar_loc_x=(npar_loc*Lxyz_loc(1)/Lxyz_loc(2))**dim1
              npar_loc_y=(npar_loc*Lxyz_loc(2)/Lxyz_loc(1))**dim1
              npar_loc_z=1
            endif
          elseif (dimensionality==1) then
!  1-D
            if (nxgrid/=1) then
              npar_loc_x=npar_loc
              npar_loc_y=1
              npar_loc_z=1
            elseif (nygrid/=1) then
              npar_loc_x=1
              npar_loc_y=npar_loc
              npar_loc_z=1
            elseif (nzgrid/=1) then
              npar_loc_x=1
              npar_loc_y=1
              npar_loc_z=npar_loc
            endif
          endif
!  Distance between particles.
          dx_par=Lxyz_loc(1)/npar_loc_x
          dy_par=Lxyz_loc(2)/npar_loc_y
          dz_par=Lxyz_loc(3)/npar_loc_z
!  Place first particle.
          fp(1,ixp) = x(l1) ; fp(1,iyp) = y(m1) ; fp(1,izp) = z(n1)
          if (nxgrid/=1) fp(1,ixp) = xyz0_loc(1)+dx_par/2
          if (nygrid/=1) fp(1,iyp) = xyz0_loc(2)+dy_par/2
          if (nzgrid/=1) fp(1,izp) = xyz0_loc(3)+dz_par/2
!  Place all other particles iteratively.
          if (dimensionality==3) then
!  3-D
            do k=2,npar_loc
              fp(k,ixp)=fp(k-1,ixp)+dx_par
              fp(k,iyp)=fp(k-1,iyp)
              fp(k,izp)=fp(k-1,izp)
              if (fp(k,ixp)>xyz1_loc(1)) then
                fp(k,ixp)=fp(1,ixp)
                fp(k,iyp)=fp(k,iyp)+dy_par
              endif
              if (fp(k,iyp)>xyz1_loc(2)) then
                fp(k,iyp)=fp(1,iyp)
                fp(k,izp)=fp(k,izp)+dz_par
              endif
            enddo
          elseif (dimensionality==2) then
!  2-D
            if (nxgrid==1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)
                fp(k,iyp)=fp(k-1,iyp)+dy_par
                fp(k,izp)=fp(k-1,izp)
                if (fp(k,iyp)>xyz1_loc(2)) then
                  fp(k,iyp)=fp(1,iyp)
                  fp(k,izp)=fp(k,izp)+dz_par
                endif
              enddo
            elseif (nygrid==1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)+dx_par
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)
                if (fp(k,ixp)>xyz1_loc(1)) then
                  fp(k,ixp)=fp(1,ixp)
                  fp(k,izp)=fp(k,izp)+dz_par
                endif
              enddo
            elseif (nzgrid==1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)+dx_par
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)
                if (fp(k,ixp)>xyz1_loc(1)) then
                  fp(k,ixp)=fp(1,ixp)
                  fp(k,iyp)=fp(k,iyp)+dy_par
                endif
              enddo
            endif
          elseif (dimensionality==1) then
!  1-D
            if (nxgrid/=1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)+dx_par
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)
              enddo
            elseif (nygrid/=1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)
                fp(k,iyp)=fp(k-1,iyp)+dy_par
                fp(k,izp)=fp(k-1,izp)
              enddo
            elseif (nzgrid/=1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)+dz_par
              enddo
            endif
          else
!  0-D
            fp(2:npar_loc,ixp)=fp(1,ixp)
            fp(2:npar_loc,iyp)=fp(1,iyp)
            fp(2:npar_loc,izp)=fp(1,izp)
          endif
          lequidistant=.true.
!
!  Shift particle locations slightly so that a mode appears.
!
        case ('shift')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) then
            if (lroot) print*, 'init_particles: must place particles equidistantly before shifting!'
            call fatal_error('init_particles','')
          endif
          k2_xxp=kx_xxp**2+ky_xxp**2+kz_xxp**2
          if (k2_xxp==0.0) then
            if (lroot) print*, &
                'init_particles: kx_xxp=ky_xxp=kz_xxp=0.0 is not allowed!'
            call fatal_error('init_particles','')
          endif
          do k=1,npar_loc
            fp(k,ixp) = fp(k,ixp) - kx_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,iyp) = fp(k,iyp) - ky_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,izp) = fp(k,izp) - kz_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
          enddo
!
        case ('gaussian-z')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=1,npar_loc
            do while (.true.)
              if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
              if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
              call random_number_wrapper(r)
              call random_number_wrapper(p)
              if (nprocz==2) then
                if (lfirst_proc_z) fp(k,izp)=-abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
                if (llast_proc_z) fp(k,izp)=abs(zp0*sqrt(-2*alog(r))*cos(2*pi*  p))
              else
                fp(k,izp)= zp0*sqrt(-2*alog(r))*cos(2*pi*p)
              endif
              if ((fp(k,izp)>=xyz0(3)).and.(fp(k,izp)<=xyz1(3))) exit
            enddo
          enddo
          if (nxgrid/=1) &
              fp(1:npar_loc,ixp)=xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
!
        case ('gaussian-x')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=1,npar_loc
            do while (.true.)
              if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
              if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
              call random_number_wrapper(r)
              call random_number_wrapper(p)
              fp(k,ixp)= xp0*sqrt(-2*alog(r))*cos(2*pi*p)
              if ((fp(k,ixp)>=xyz0(1)).and.(fp(k,ixp)<=xyz1(1))) exit
            enddo
          enddo
          if (nygrid/=1) &
              fp(1:npar_loc,iyp)=xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          if (nzgrid/=1) &
              fp(1:npar_loc,izp)=xyz0_par(3)+fp(1:npar_loc,izp)*Lxyz_par(3)
!
        case ('gaussian-z-pure')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=1,npar_loc
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            if (nprocz==2) then
              if (lfirst_proc_z) fp(k,izp)=-abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
              if (llast_proc_z) fp(k,izp)=abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
            else
              fp(k,izp)= zp0*sqrt(-2*alog(r))*cos(2*pi*p)
            endif
          enddo
!
        case ('gaussian-r')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=1,npar_loc
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            call random_number_wrapper(q)
            fp(k,ixp)= xp0*sqrt(-2*alog(r))*cos(2*pi*p)*cos(2*pi*q)
            fp(k,iyp)= yp0*sqrt(-2*alog(r))*cos(2*pi*p)*sin(2*pi*q)
          enddo
!
        case ('hole')
!
          call map_nearest_grid(fp,ineargrid)
          call map_xxp_grid(f,fp,ineargrid)
          call sort_particles_imn(fp,ineargrid,ipar)
          do k=k1_imn(imn_array(m_hole+m1-1,n_hole+n1-1)), &
               k2_imn(imn_array(m_hole+m1-1,n_hole+n1-1))
            if (ineargrid(k,1)==l_hole+l1-1) then
              print*, k
              if (nxgrid/=0) fp(k,ixp)=fp(k,ixp)-dx
            endif
          enddo
!
        case ('streaming')
          call streaming(fp,f)
!
        case ('streaming_coldstart')
          call streaming_coldstart(fp,f)
!
        case ('constant-Ri')
          call constant_richardson(fp,f)
!
        case default
          if (lroot) &
              print*, 'init_particles: No such such value for initxxp: ', &
              trim(initxxp(j))
          call fatal_error('init_particles','')
!
        endselect
!
      enddo ! do j=1,ninit
!
!  Interface for user's own initial condition for position
!
      if (linitial_condition) call initial_condition_xxp(f,fp)
!
!  Particles are not allowed to be present in non-existing dimensions.
!  This would give huge problems with interpolation later.
!
      if (nxgrid==1) fp(1:npar_loc,ixp)=x(nghost+1)
      if (nygrid==1) fp(1:npar_loc,iyp)=y(nghost+1)
      if (nzgrid==1) fp(1:npar_loc,izp)=z(nghost+1)
!
      if (init_repeat/=0) call repeated_init(fp,init_repeat)
!
!  Redistribute particles among processors (now that positions are determined).
!
      call boundconds_particles(fp,ipar)
!
!  Map particle position on the grid.
!
      call map_nearest_grid(fp,ineargrid)
      call map_xxp_grid(f,fp,ineargrid)
!
!  Initial particle velocity.
!
      do j=1,ninit
!
        select case (initvvp(j))
!
        case ('nothing')
          if (lroot.and.j==1) print*, 'init_particles: No particle velocity set'
!
        case ('zero')
          if (lroot) print*, 'init_particles: Zero particle velocity'
          fp(1:npar_loc,ivpx:ivpz)=0.0
!
        case ('zero-shear')
          if (lroot) print*, 'init_particles: Zero particle velocity'
          fp(1:npar_loc,ivpy)=-Sshear*fp(1:npar_loc,ixp)
          fp(1:npar_loc,ivpx)=0.0
          fp(1:npar_loc,ivpz)=0.0
!
        case ('constant')
          if (lroot) print*, 'init_particles: Constant particle velocity'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          if (lcylindrical_coords) then
            fp(1:npar_loc,ivpx)&
                =vpx0*cos(fp(k,iyp))+vpy0*sin(fp(k,iyp))
            fp(1:npar_loc,ivpy)&
                =vpy0*cos(fp(k,iyp))-vpx0*sin(fp(k,iyp))
            fp(1:npar_loc,ivpz)=vpz0
          else
            fp(1:npar_loc,ivpx)=vpx0
            fp(1:npar_loc,ivpy)=vpy0
            fp(1:npar_loc,ivpz)=vpz0
          endif
!
        case ('constant-1')
          if (lroot) &
              print*, 'init_particles: Particle 1 velocity vx,vy,vz=', &
              vpx1, vpy1, vpz1
          do k=1,npar_loc
            if (ipar(k)==1) then
              fp(k,ivpx)=vpx1
              fp(k,ivpy)=vpy1
              fp(k,ivpz)=vpz1
            endif
          enddo
!
        case ('constant-2')
          if (lroot) &
              print*, 'init_particles: Particle 2 velocity vx,vy,vz=', &
              vpx2, vpy2, vpz2
          do k=1,npar_loc
            if (ipar(k)==2) then
              fp(k,ivpx)=vpx2
              fp(k,ivpy)=vpy2
              fp(k,ivpz)=vpz2
            endif
          enddo
!
        case ('constant-3')
          if (lroot) &
              print*, 'init_particles: Particle 3 velocity vx,vy,vz=', &
              vpx3, vpy3, vpz3
          do k=1,npar_loc
            if (ipar(k)==3) then
              fp(k,ivpx)=vpx3
              fp(k,ivpy)=vpy3
              fp(k,ivpz)=vpz3
            endif
          enddo
!
        case ('sinwave-phase')
          if (lroot) print*, 'init_particles: sinwave-phase'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k=1,npar_loc
            fp(k,ivpx)=fp(k,ivpx)+vpx0*sin(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy)=fp(k,ivpy)+vpy0*sin(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz)=fp(k,ivpz)+vpz0*sin(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
!
        case ('coswave-phase')
          if (lroot) print*, 'init_particles: coswave-phase'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k=1,npar_loc
            fp(k,ivpx)=fp(k,ivpx)+vpx0*cos(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy)=fp(k,ivpy)+vpy0*cos(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz)=fp(k,ivpz)+vpz0*cos(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle velocities; '// &
              'delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpx) = fp(k,ivpx) + delta_vp0*(2*r-1)
            call random_number_wrapper(r)
            fp(k,ivpy) = fp(k,ivpy) + delta_vp0*(2*r-1)
            call random_number_wrapper(r)
            fp(k,ivpz) = fp(k,ivpz) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-x')
          if (lroot) print*, 'init_particles: Random particle x-velocity; '// &
              'delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpx) = fp(k,ivpx) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-y')
          if (lroot) print*, 'init_particles: Random particle y-velocity; '// &
              'delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpy) = fp(k,ivpy) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-z')
          if (lroot) print*, 'init_particles: Random particle z-velocity; '// &
              'delta_vp0=', delta_vp0
          do k=1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpz) = fp(k,ivpz) + delta_vp0*(2*r-1)
          enddo
!
        case ('average-to-zero')
          call mpireduce_sum(sum(fp(1:npar_loc,ivpx)),vpx_sum)
          call mpireduce_sum(sum(fp(1:npar_loc,ivpy)),vpy_sum)
          call mpireduce_sum(sum(fp(1:npar_loc,ivpz)),vpz_sum)
          call mpibcast_real(vpx_sum)
          call mpibcast_real(vpy_sum)
          call mpibcast_real(vpz_sum)
          fp(1:npar_loc,ivpx)=fp(1:npar_loc,ivpx)-vpx_sum/npar
          fp(1:npar_loc,ivpy)=fp(1:npar_loc,ivpy)-vpy_sum/npar
          fp(1:npar_loc,ivpz)=fp(1:npar_loc,ivpz)-vpz_sum/npar
!
        case ('follow-gas')
          if (lroot) &
              print*, 'init_particles: Particle velocity equal to gas velocity'
          do k=1,npar_loc
            call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup, &
                ineargrid(k,:),0,0)
            fp(k,ivpx:ivpz) = uup
          enddo
!
        case ('jeans-wave-dustpar-x')
        ! assumes rhs_poisson_const=1 !
          do k=1,npar_loc
            fp(k,ivpx) = fp(k,ivpx) - amplxxp* &
                (sqrt(1+4*1.0*1.0*tausp**2)-1)/ &
                (2*kx_xxp*1.0*tausp)*sin(kx_xxp*(fp(k,ixp)))
          enddo
!
        case ('dragforce_equilibrium','dragforce-equilibrium')
!
!  Equilibrium between drag forces on dust and gas and other forces
!  (from Nakagawa, Sekiya, & Hayashi 1986).
!
          if (lroot) then
            print*, 'init_particles: drag equilibrium'
            print*, 'init_particles: beta_glnrho_global=', beta_glnrho_global
          endif
          cs=sqrt(cs20)
!  Calculate average dust-to-gas ratio in box.
          if (ldensity_nolog) then
            eps = sum(f(l1:l2,m1:m2,n1:n2,irhop))/ &
                sum(f(l1:l2,m1:m2,n1:n2,irho))
          else
            eps = sum(f(l1:l2,m1:m2,n1:n2,irhop))/ &
                sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))
          endif
          if (ldragforce_equi_noback) eps=0.0
!
          if (lroot) print*, 'init_particles: average dust-to-gas ratio=', eps
!  Set gas velocity field.
          if (lhydro) then
            do l=l1,l2; do m=m1,m2; do n=n1,n2
!  Take either global or local dust-to-gas ratio.
              if (.not. ldragforce_equi_global_eps) then
                if (ldensity_nolog) then
                  eps=f(l,m,n,irhop)/f(l,m,n,irho)
                else
                  eps=f(l,m,n,irhop)/exp(f(l,m,n,ilnrho))
                endif
              endif
!
              f(l,m,n,iux) = f(l,m,n,iux) - &
                  beta_glnrho_global(1)*eps*Omega*tausp/ &
                  ((1.0+eps)**2+(Omega*tausp)**2)*cs
              f(l,m,n,iuy) = f(l,m,n,iuy) + &
                  beta_glnrho_global(1)*(1+eps+(Omega*tausp)**2)/ &
                  (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
!
            enddo; enddo; enddo
          endif
!  Set particle velocity field.
          do k=1,npar_loc
!  Take either global or local dust-to-gas ratio.
            if (ldragforce_equi_noback) then
              eps=0.0
            else
              if (.not. ldragforce_equi_global_eps) then
                ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
                if (ldensity_nolog) then
                  eps=f(ix0,iy0,iz0,irhop)/f(ix0,iy0,iz0,irho)
                else
                  eps=f(ix0,iy0,iz0,irhop)/exp(f(ix0,iy0,iz0,ilnrho))
                endif
              endif
            endif
!
            fp(k,ivpx) = fp(k,ivpx) + &
                beta_glnrho_global(1)*Omega*tausp/ &
                ((1.0+eps)**2+(Omega*tausp)**2)*cs
            fp(k,ivpy) = fp(k,ivpy) + &
                beta_glnrho_global(1)*(1+eps)/ &
                (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
!
          enddo
!
        case ('dragforce_equi_nohydro')
!
          do k=1,npar_loc
            fp(k,ivpx) = fp(k,ivpx) - 2*Deltauy_gas_friction* &
                1/(1.0/(Omega*tausp)+Omega*tausp)
            fp(k,ivpy) = fp(k,ivpy) - Deltauy_gas_friction* &
                1/(1.0+(Omega*tausp)**2)
          enddo
!
        case ('dragforce_equi_dust')
!
!  Equilibrium between drag force and Coriolis force on the dust.
!
          if (lroot) then
            print*, 'init_particles: drag equilibrium dust'
            print*, 'init_particles: beta_dPdr_dust=', beta_dPdr_dust
          endif
!  Set particle velocity field.
          cs=sqrt(cs20)
          do k=1,npar_loc
            fp(k,ivpx) = fp(k,ivpx) + &
                1/(Omega*tausp+1/(Omega*tausp))*beta_dPdr_dust*cs
            fp(k,ivpy) = fp(k,ivpy) - &
                1/(1.0+1/(Omega*tausp)**2)*beta_dPdr_dust/2*cs
          enddo
!
       case ('Keplerian','keplerian')
!
!  Keplerian velocity based on gravr.
!
          if (lroot) then
            print*, 'init_particles: Keplerian velocity'
            if (lspherical_coords) call fatal_error('init_particles', &
                 'Keplerian particle initial condition: '// &
                 'not implemented for spherical coordinates')
            if (lshear) call fatal_error("init_particles",&
                 "Keplerian initial condition is for global disks, not shearing boxes")
          endif
          do k=1,npar_loc
            if (lcartesian_coords) then
              rad=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2)
              OO=sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) = -OO*fp(k,iyp)
              fp(k,ivpy) =  OO*fp(k,ixp)
              fp(k,ivpz) =  0.0
            elseif (lcylindrical_coords) then
              rad=fp(k,ixp)
              OO=sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) =  0.0
              fp(k,ivpy) =  OO*rad
              fp(k,ivpz) =  0.0
            endif
          enddo
!
!  Explosion.
!
       case ('explosion')
         do k=1,npar_loc
           rad=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2)
           fp(k,ivpx) = delta_vp0*fp(k,ixp)/rp_ext
           fp(k,ivpy) = delta_vp0*fp(k,iyp)/rp_ext
           fp(k,ivpz) = delta_vp0*fp(k,izp)/rp_ext
         enddo
!
!
        case default
          if (lroot) &
              print*, 'init_particles: No such such value for initvvp: ', &
              trim(initvvp(j))
          call fatal_error('','')
!
        endselect
!
      enddo ! do j=1,ninit
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_vvp(f,fp)
!
!  Map particle velocity on the grid.
!
      call map_vvp_grid(f,fp,ineargrid)
!
!  Sort particles (must happen at the end of the subroutine so that random
!  positions and velocities are not displaced relative to when there is no
!  sorting).
!
      call sort_particles_imn(fp,ineargrid,ipar)
!
    endsubroutine init_particles
!***********************************************************************
    subroutine insert_lost_particles(f,fp,ineargrid)
!
!  14-oct-12/dhruba: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (inout) :: fp,ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine insert_lost_particles
!***********************************************************************
    subroutine insert_particles(f,fp,ineargrid)
!
! Insert particles continuously (when linsert_particles_continuously == T),
! i.e. in each timestep. If number of particles to be inserted are less
! than unity, accumulate number over several timesteps until the integer value
! is larger than one. Keep the remainder and accumulate this to the next insert.
!
! Works only for particles_dust - add neccessary variable
! declarations in particles_tracers to make it work here.
!
      use General, only: random_number_wrapper
      use Particles_diagnos_state, only: insert_particles_diagnos_state
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      logical :: linsertmore=.true.
      real :: xx0, yy0,r2
!
      integer :: j, k, n_insert, npar_loc_old, iii
!
      intent (inout) :: fp,ineargrid
!
! Stop call to this routine when maximum number of particles is reached!
! Since root inserts all new particles, make sure
! npar_total + n_insert < mpar
! so that a processor can not exceed its maximum number of particles.
!
      if (lroot) then
        avg_n_insert=particles_insert_rate*dt
        n_insert=int(avg_n_insert + remaining_particles)
! Remaining particles saved for subsequent timestep:
        remaining_particles=avg_n_insert + remaining_particles - n_insert
        if ((n_insert+npar_total <= mpar_loc) &
            .and. (t<max_particle_insert_time)) then
          linsertmore=.true.
        else
          linsertmore=.false.
        endif
!
        if (linsertmore) then
! Actual (integer) number of particles to be inserted at this timestep:
          do iii=npar_loc+1,npar_loc+n_insert
            ipar(iii)=npar_total+iii-npar_loc
          enddo
          npar_total=npar_total+n_insert
          npar_loc_old=npar_loc
          npar_loc=npar_loc + n_insert
!
! Insert particles in chosen position (as in init_particles).
!
          do j=1,ninit
            select case (initxxp(j))
            case ('random-box')
!
              do k=npar_loc_old+1,npar_loc
                if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
                if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
                if (nzgrid/=1) call random_number_wrapper(fp(k,izp))
                if (lcylindrical_coords) then
                  xx0=xp0+fp(k,ixp)*Lx0
                  yy0=yp0+fp(k,iyp)*Ly0
                  r2=xx0**2+yy0**2
                  if (nxgrid/=1) fp(k,ixp)=sqrt(r2)
                  if (nygrid/=1) fp(k,iyp)=atan(yy0/xx0)+pi*(xx0/abs(xx0)-1)*0.5
                  if (nzgrid/=1) fp(k,izp)=zp0+fp(k,izp)*Lz0
                else
                  if (nxgrid/=1) fp(k,ixp)=xp0+fp(k,ixp)*Lx0
                  if (nygrid/=1) fp(k,iyp)=yp0+fp(k,iyp)*Ly0
                  if (nzgrid/=1) fp(k,izp)=zp0+fp(k,izp)*Lz0
                endif
              enddo
!
            case ('nothing')
              if (lroot .and. j==1) print*, 'init_particles: nothing'
!
            case default
              print*, 'insert_particles: No such such value for initxxp: ', &
                  trim(initxxp(j))
              call fatal_error('init_particles','')
!
            endselect
          enddo
!
!  Initial particle velocity.
!
          do j=1,ninit
            select case (initvvp(j))
            case ('nothing')
              if (j==1) print*, 'init_particles: No particle velocity set'
!
            case ('constant')
              if (lcylindrical_coords) then
                fp(npar_loc_old+1:npar_loc,ivpx)&
                    =vpx0*cos(fp(npar_loc_old+1:npar_loc,iyp))&
                    +vpy0*sin(fp(npar_loc_old+1:npar_loc,iyp))
                fp(npar_loc_old+1:npar_loc,ivpy)&
                    =vpy0*cos(fp(npar_loc_old+1:npar_loc,iyp))&
                    -vpx0*sin(fp(npar_loc_old+1:npar_loc,iyp))
                fp(npar_loc_old+1:npar_loc,ivpz)=vpz0
              else
                fp(npar_loc_old+1:npar_loc,ivpx)=vpx0
                fp(npar_loc_old+1:npar_loc,ivpy)=vpy0
                fp(npar_loc_old+1:npar_loc,ivpz)=vpz0
              endif
!
            case default
              print*, 'insert_particles: No such such value for initvvp: ', &
                  trim(initvvp(j))
              call fatal_error('','')
              !
            endselect
!
          enddo ! do j=1,ninit
!
!  Initialize particle radius
!
          call set_particle_radius(f,fp,npar_loc_old+1,npar_loc)
!
!  Particles are not allowed to be present in non-existing dimensions.
!  This would give huge problems with interpolation later.
!
          if (nxgrid==1) fp(npar_loc_old+1:npar_loc,ixp)=x(nghost+1)
          if (nygrid==1) fp(npar_loc_old+1:npar_loc,iyp)=y(nghost+1)
          if (nzgrid==1) fp(npar_loc_old+1:npar_loc,izp)=z(nghost+1)
!
          if (lparticles_diagnos_state) &
              call insert_particles_diagnos_state(fp, npar_loc_old)
!
        endif
      endif ! if (lroot) then
!
!  Redistribute particles only when t < max_particle_insert_time.
!  Could have included some other tests here aswell......
!
      if (t<max_particle_insert_time) then
!
!  Redistribute particles among processors.
!
        call boundconds_particles(fp,ipar,linsert=.true.)
!
!  Map particle position on the grid.
!
        call map_nearest_grid(fp,ineargrid)
        call map_xxp_grid(f,fp,ineargrid)
!
!  Map particle velocity on the grid.
!
        call map_vvp_grid(f,fp,ineargrid)
!
!  Sort particles (must happen at the end of the subroutine so that random
!  positions and velocities are not displaced relative to when there is no
!  sorting).
!
        call sort_particles_imn(fp,ineargrid,ipar)
      endif
!
    endsubroutine insert_particles
!***********************************************************************
    subroutine streaming_coldstart(fp,f)
!
!  Mode that is unstable to the streaming instability of Youdin & Goodman (2005)
!
!  14-apr-06/anders: coded
!
      use EquationOfState, only: beta_glnrho_global
      use General, only: random_number_wrapper
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: eta_glnrho, v_Kepler, ampluug, dxp, dzp
      integer :: i, i1, i2, j, k, npar_loc_x, npar_loc_z
!
!  The number of particles per grid cell must be a quadratic number.
!
      if ( sqrt(npar/real(nwgrid))/=int(sqrt(npar/real(nwgrid))) .or. &
           sqrt(npar_loc/real(nw))/=int(sqrt(npar_loc/real(nw))) ) then
        if (lroot) then
          print*, 'streaming_coldstart: the number of particles per grid must'
          print*, '                     be a quadratic number!'
        endif
        print*, '                     iproc, npar/nw, npar_loc/nwgrid=', &
            iproc, npar/real(nwgrid), npar_loc/real(nw)
        call fatal_error('streaming_coldstart','')
      endif
!
!  Define a few disc parameters.
!
      eta_glnrho = -0.5*abs(beta_glnrho_global(1))*beta_glnrho_global(1)
      v_Kepler   =  1.0/abs(beta_glnrho_global(1))
      if (lroot) print*, 'streaming: eta, vK=', eta_glnrho, v_Kepler
!
!  Place particles equidistantly.
!
      npar_loc_x=sqrt(npar_loc/(Lxyz_loc(3)/Lxyz_loc(1)))
      npar_loc_z=npar_loc/npar_loc_x
      dxp=Lxyz_loc(1)/npar_loc_x
      dzp=Lxyz_loc(3)/npar_loc_z
      do i=1,npar_loc_x
        i1=(i-1)*npar_loc_z+1; i2=i*npar_loc_z
        fp(i1:i2,ixp)=mod(i*dxp,Lxyz_loc(1))+dxp/2
        do j=i1,i2
          fp(j,izp)=xyz0_loc(3)+dzp/2+(j-i1)*dzp
        enddo
      enddo
!
!  Shift particle locations slightly so that wanted mode appears.
!
      do k=1,npar_loc
        fp(k,ixp) = fp(k,ixp) - &
            amplxxp/(2*(kx_xxp**2+kz_xxp**2))* &
            (kx_xxp*sin(kx_xxp*fp(k,ixp)+kz_xxp*fp(k,izp))+ &
             kx_xxp*sin(kx_xxp*fp(k,ixp)-kz_xxp*fp(k,izp)))
        fp(k,izp) = fp(k,izp) - &
            amplxxp/(2*(kx_xxp**2+kz_xxp**2))* &
            (kz_xxp*sin(kx_xxp*fp(k,ixp)+kz_xxp*fp(k,izp))- &
             kz_xxp*sin(kx_xxp*fp(k,ixp)-kz_xxp*fp(k,izp)))
        fp(k,ixp) = fp(k,ixp) + &
            kx_xxp/(2*(kx_xxp**2+kz_xxp**2))*amplxxp**2* &
            sin(2*(kx_xxp*fp(k,ixp)+kz_xxp*fp(k,izp)))
        fp(k,izp) = fp(k,izp) + &
            kz_xxp/(2*(kx_xxp**2+kz_xxp**2))*amplxxp**2* &
            sin(2*(kx_xxp*fp(k,ixp)+kz_xxp*fp(k,izp)))
      enddo
!  Set particle velocity.
      do k=1,npar_loc
        fp(k,ivpx) = fp(k,ivpx) + eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(1))*cos(kx_xxp*fp(k,ixp)) - &
             aimag(coeff(1))*sin(kx_xxp*fp(k,ixp)))*cos(kz_xxp*fp(k,izp))
        fp(k,ivpy) = fp(k,ivpy) + eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(2))*cos(kx_xxp*fp(k,ixp)) - &
             aimag(coeff(2))*sin(kx_xxp*fp(k,ixp)))*cos(kz_xxp*fp(k,izp))
        fp(k,ivpz) = fp(k,ivpz) + eta_glnrho*v_Kepler*(-amplxxp)* &
            (aimag(coeff(3))*cos(kx_xxp*fp(k,ixp)) + &
              real(coeff(3))*sin(kx_xxp*fp(k,ixp)))*sin(kz_xxp*fp(k,izp))
      enddo
!
!  Change the gas velocity amplitude so that the numerical error on the drag
!  force is corrected (the error is due to the interpolation of the gas
!  velocity field to the positions of the particles). A better way to correct
!  this is to go to a quadratic interpolation scheme.
!
      ampluug=amplxxp
      if (lcoldstart_amplitude_correction) &
          ampluug=amplxxp/(1-dx**2/8*(kx_xxp**2+kz_xxp**2))
!
!  Set fluid fields.
!
      do m=m1,m2; do n=n1,n2
        f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + &
            amplxxp* &
            ( real(coeff(7))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(7))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
        f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + &
            eta_glnrho*v_Kepler*ampluug* &
            ( real(coeff(4))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(4))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
        f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
            eta_glnrho*v_Kepler*ampluug* &
            ( real(coeff(5))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(5))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
        f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + &
            eta_glnrho*v_Kepler*(-ampluug)* &
            (aimag(coeff(6))*cos(kx_xxp*x(l1:l2)) + &
              real(coeff(6))*sin(kx_xxp*x(l1:l2)))*sin(kz_xxp*z(n))
      enddo; enddo
!
    endsubroutine streaming_coldstart
!***********************************************************************
    subroutine streaming(fp,f)
!
!  Mode that is unstable to the streaming instability of Youdin & Goodman (2005)
!
!  30-jan-06/anders: coded
!
      use EquationOfState, only: beta_glnrho_global
      use General, only: random_number_wrapper
      use Particles_mpicomm
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: eta_glnrho, v_Kepler, kx, kz
      real :: r, p, xprob, zprob, dzprob, fprob, dfprob
      integer :: j, k
      logical :: lmigration_redo_org
!
!  Define a few disc parameters.
!
      eta_glnrho = -0.5*abs(beta_glnrho_global(1))*beta_glnrho_global(1)
      v_Kepler   =  1.0/abs(beta_glnrho_global(1))
      if (lroot) print*, 'streaming: eta, vK=', eta_glnrho, v_Kepler
!
!  Place particles according to probability function.
!
!  Invert
!    r = x
!    p = int_0^z f(x,z') dz' = z + A/kz*cos(kx*x)*sin(kz*z)
!  where r and p are random numbers between 0 and 1.
      kx=kx_xxp*Lxyz(1); kz=kz_xxp*Lxyz(3)
      do k=1,npar_loc
!
        call random_number_wrapper(r)
        call random_number_wrapper(p)
!
        fprob = 1.0
        zprob = 0.0
!
        j=0
!  Use Newton-Raphson iteration to invert function.
        do while ( abs(fprob)>0.0001 )
!
          xprob = r
          fprob = zprob + amplxxp/kz*cos(kx*xprob)*sin(kz*zprob) - p
          dfprob= 1.0 + amplxxp*cos(kx*xprob)*cos(kz*zprob)
          dzprob= -fprob/dfprob
          zprob = zprob+0.2*dzprob
!
          j=j+1
!
        enddo
!
        if ( mod(k,npar_loc/100)==0) then
          print '(i7,i3,4f11.7)', k, j, r, p, xprob, zprob
        endif
!
        fp(k,ixp)=xprob*Lxyz(1)+xyz0(1)
        fp(k,izp)=zprob*Lxyz(3)+xyz0(3)
!  Set particle velocity.
        fp(k,ivpx) = fp(k,ivpx) + eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(1))*cos(kx_xxp*fp(k,ixp)) - &
             aimag(coeff(1))*sin(kx_xxp*fp(k,ixp)))*cos(kz_xxp*fp(k,izp))
        fp(k,ivpy) = fp(k,ivpy) + eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(2))*cos(kx_xxp*fp(k,ixp)) - &
             aimag(coeff(2))*sin(kx_xxp*fp(k,ixp)))*cos(kz_xxp*fp(k,izp))
        fp(k,ivpz) = fp(k,ivpz) + eta_glnrho*v_Kepler*(-amplxxp)* &
            (aimag(coeff(3))*cos(kx_xxp*fp(k,ixp)) + &
              real(coeff(3))*sin(kx_xxp*fp(k,ixp)))*sin(kz_xxp*fp(k,izp))
!
      enddo
!
!  Particles were placed randomly in the entire simulation space, so they need
!  to be send to the correct processors now.
!
      if (lmpicomm) then
        lmigration_redo_org=lmigration_redo
        lmigration_redo=.true.
        call migrate_particles(fp,ipar)
        lmigration_redo=lmigration_redo_org
      endif
!
!  Set fluid fields.
!
      do m=m1,m2; do n=n1,n2
        f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + &
            (eta_glnrho*v_Kepler)**2*amplxxp* &
            ( real(coeff(7))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(7))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
        f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + &
            eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(4))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(4))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
        f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
            eta_glnrho*v_Kepler*amplxxp* &
            ( real(coeff(5))*cos(kx_xxp*x(l1:l2)) - &
             aimag(coeff(5))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
        f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + &
            eta_glnrho*v_Kepler*(-amplxxp)* &
            (aimag(coeff(6))*cos(kx_xxp*x(l1:l2)) + &
              real(coeff(6))*sin(kx_xxp*x(l1:l2)))*sin(kz_xxp*z(n))
      enddo; enddo
!
    endsubroutine streaming
!***********************************************************************
    subroutine constant_richardson(fp,f)
!
!  Setup dust density with a constant Richardson number (Sekiya, 1998).
!    eps=1/sqrt(z^2/Hd^2+1/(1+eps1)^2)-1
!
!  14-sep-05/anders: coded
!
      use EquationOfState, only: beta_glnrho_scaled, gamma, cs20
      use General, only: random_number_wrapper
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer, parameter :: nz_inc=10
      real, dimension (nz_inc*nz) :: z_dense, eps
      real :: r, Hg, Hd, frac, rho1, Sigmad, Sigmad_num, Xi, fXi, dfdXi
      real :: dz_dense, eps_point, z00_dense, rho, lnrho
      integer :: nz_dense=nz_inc*nz, npar_bin
      integer :: i, i0, k
!
!  Calculate dust "scale height".
!
      rho1=1.0
      Hg=1.0
      Sigmad=eps_dtog*rho1*Hg*sqrt(2*pi)
      Hd = sqrt(Ri0)*abs(beta_glnrho_scaled(1))/2*1.0
!
!  Need to find eps1 that results in given dust column density.
!
      Xi = sqrt(eps1*(2+eps1))/(1+eps1)
      fXi=-2*Xi + alog((1+Xi)/(1-Xi))-Sigmad/(Hd*rho1)
      i=0
!
!  Newton-Raphson on equation Sigmad/(Hd*rho1)=-2*Xi + alog((1+Xi)/(1-Xi)).
!  Here Xi = sqrt(eps1*(2+eps1))/(1+eps1).
!
      do while (abs(fXi)>=0.00001)
!
        dfdXi=2*Xi**2/(1-Xi**2)
        Xi=Xi-0.1*fXi/dfdXi
!
        fXi=-2*Xi + alog((1+Xi)/(1-Xi))-Sigmad/(Hd*rho1)
!
        i=i+1
        if (i>=1000) stop
!
      enddo
!
!  Calculate eps1 from Xi.
!
      eps1=-1+1/sqrt(-(Xi**2)+1)
      if (lroot) print*, 'constant_richardson: Hd, eps1=', Hd, eps1
!
!  Make z denser for higher resolution in density.
!
      dz_dense=Lxyz_loc(3)/nz_dense
      z00_dense=xyz0_loc(3)+0.5*dz_dense
      do n=1,nz_dense
        z_dense(n)=z00_dense+(n-1)*dz_dense
      enddo
!
!  Dust-to-gas ratio as a function of z (with cutoff).
!
      eps=1/sqrt(z_dense**2/Hd**2+1/(1+eps1)**2)-1
      where (eps<=0.0) eps=0.0
!
!  Calculate the dust column density numerically.
!
      Sigmad_num=sum(rho1*eps*dz_dense)
      if (lroot) print*, 'constant_richardson: Sigmad, Sigmad (numerical) = ', &
          Sigmad, Sigmad_num
!
!  Place particles according to probability function.
!
      i0=0
      do n=1,nz_dense
        frac=eps(n)/Sigmad_num*dz_dense
        npar_bin=int(frac*npar_loc)
        if (npar_bin>=2.and.mod(n,2)==0) npar_bin=npar_bin+1
        do i=i0+1,i0+npar_bin
          if (i<=npar_loc) then
            call random_number_wrapper(r)
            fp(i,izp)=z_dense(n)+(2*r-1.0)*dz_dense/2
          endif
        enddo
        i0=i0+npar_bin
      enddo
      if (lroot) print '(A,i7,A)', 'constant_richardson: placed ', &
          i0, ' particles according to Ri=const'
!
!  Particles left out by round off are just placed randomly.
!
      if (i0+1<=npar_loc) then
        do k=i0+1,npar_loc
          call random_number_wrapper(fp(k,izp))
          fp(k,izp)=xyz0(3)+fp(k,izp)*Lxyz(3)
        enddo
        if (lroot) print '(A,i7,A)', 'constant_richardson: placed ', &
            npar_loc-i0, ' particles randomly.'
      endif
!
!  Random positions in x and y.
!
      do k=1,npar_loc
        if (nxgrid/=1) call random_number_wrapper(fp(k,ixp))
        if (nygrid/=1) call random_number_wrapper(fp(k,iyp))
      enddo
      if (nxgrid/=1) &
          fp(1:npar_loc,ixp)=xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
      if (nygrid/=1) &
          fp(1:npar_loc,iyp)=xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
!
!  Set gas velocity according to dust-to-gas ratio and global pressure gradient.
!
      do imn=1,ny*nz
!
        n=nn(imn); m=mm(imn)
!
        if (abs(z(n))<=Hd*sqrt(1-1/(1+eps1)**2)) then
          lnrho = -sqrt(z(n)**2/Hd**2+1/(1+eps1)**2)* &
              gamma*Omega**2*Hd**2/cs20 + gamma*Omega**2*Hd**2/(cs20*(1+eps1))
        else
          lnrho = -0.5*gamma*Omega**2/cs20*z(n)**2 + &
              gamma*Omega**2*Hd**2/cs20*(1/(1+eps1)-1/(2*(1+eps1)**2) - 0.5)
        endif
!
!  Isothermal stratification.
!
        if (lentropy) f(l1:l2,m,n,iss) = (1/gamma-1.0)*lnrho
!
        rho=exp(lnrho)
!
        if (ldensity_nolog) then
          f(l1:l2,m,n,irho)  =rho
        else
          f(l1:l2,m,n,ilnrho)=lnrho
        endif
!
        eps_point=1/sqrt(z(n)**2/Hd**2+1/(1+eps1)**2)-1
        if (eps_point<=0.0) eps_point=0.0
!
        f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) - &
            cs20*beta_glnrho_scaled(1)*eps_point*tausp/ &
            (1.0+2*eps_point+eps_point**2+(Omega*tausp)**2)
        f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
            cs20*beta_glnrho_scaled(1)*(1+eps_point+(Omega*tausp)**2)/ &
            (2*Omega*(1.0+2*eps_point+eps_point**2+(Omega*tausp)**2))
        f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + 0.0
      enddo
!
!  Set particle velocity.
!
      do k=1,npar_loc
!
        eps_point=1/sqrt(fp(k,izp)**2/Hd**2+1/(1+eps1)**2)-1
        if (eps_point<=0.0) eps_point=0.0
!
        fp(k,ivpx) = fp(k,ivpx) + &
            cs20*beta_glnrho_scaled(1)*tausp/ &
            (1.0+2*eps_point+eps_point**2+(Omega*tausp)**2)
        fp(k,ivpy) = fp(k,ivpy) + &
            cs20*beta_glnrho_scaled(1)*(1+eps_point)/ &
            (2*Omega*(1.0+2*eps_point+eps_point**2+(Omega*tausp)**2))
        fp(k,ivpz) = fp(k,ivpz) - tausp*Omega**2*fp(k,izp)
!
      enddo
!
    endsubroutine constant_richardson
!***********************************************************************
    subroutine particles_dragforce_stiff(f,fp,ineargrid)
!
!  Force stiff drag force equations towards their equilibrium.
!
!  10-june-11/anders: coded
!
      use Boundcond
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx) :: eps
      real, dimension (3) :: vvp
      integer :: imn, i, k, ix0, iy0, iz0
!
      if (ldragforce_stiff .and. .not. lpencil_check_at_work) then
        do imn=1,ny*nz
          n=nn(imn)
          m=mm(imn)
          eps=f(l1:l2,m,n,irhop)/f(l1:l2,m,n,irho)
          do i=0,2
            f(l1:l2,m,n,iux+i)=(f(l1:l2,m,n,iux+i)+eps*f(l1:l2,m,n,iupx+i) + &
                eps/(1.0+eps)*tausp*f(l1:l2,m,n,ifgx+i))/(1.0+eps)
          enddo
          f(l1:l2,m,n,ifgx:ifgz)=f(l1:l2,m,n,iux:iuz)-f(l1:l2,m,n,ifgx:ifgz)
        enddo
        call boundconds_x(f,ifgx,ifgz)
        call initiate_isendrcv_bdry(f,ifgx,ifgz)
        call finalize_isendrcv_bdry(f,ifgx,ifgz)
        call boundconds_y(f,ifgx,ifgz)
        call boundconds_z(f,ifgx,ifgz)
        do k=1,npar_loc
          ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
          if (lparticlemesh_cic) then
            call interpolate_linear(f,ifgx,ifgz, &
                fp(k,ixp:izp),vvp,ineargrid(k,:),0,ipar(k))
          elseif (lparticlemesh_tsc) then
            if (linterpolate_spline) then
              call interpolate_quadratic_spline(f,ifgx,ifgz, &
                  fp(k,ixp:izp),vvp,ineargrid(k,:),0,ipar(k))
            else
              call interpolate_quadratic(f,ifgx,ifgz, &
                  fp(k,ixp:izp),vvp,ineargrid(k,:),0,ipar(k))
            endif
          else
            vvp=f(ix0,iy0,iz0,ifgx:ifgz)
          endif
          fp(k,ivpx:ivpz)=vvp
        enddo
      endif
!
    endsubroutine particles_dragforce_stiff
!***********************************************************************
    subroutine pencil_criteria_particles()
!
!  All pencils that the Particles module depends on are specified here.
!
!  20-04-06/anders: coded
!
      if (ldragforce_gas_par) then
        lpenc_requested(i_epsp)=.true.
        lpenc_requested(i_np)=.true.
      endif
      if (ldragforce_heat .or. lcollisional_heat) then
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
      if (lcollisional_cooling_taucool) then
        lpenc_requested(i_np)=.true.
      endif
      if (lcollisional_cooling_rms) then
        lpenc_requested(i_epsp)=.true.
      endif
      if (lcollisional_cooling_rms .or. lcollisional_dragforce_cooling) then
        lpenc_requested(i_np)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
      if (ldraglaw_epstein_transonic  .or.&
          ldraglaw_eps_stk_transonic) then
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_cs2)=.true.
      endif
      if (ldragforce_stiff) then
        lpenc_requested(i_fpres)=.true.
        lpenc_requested(i_jxbr)=.true.
        lpenc_requested(i_fvisc)=.true.
      endif
      if (lthermophoretic_forces) then
        lpenc_requested(i_gTT)=.true.
      endif
!
      if (idiag_npm/=0 .or. idiag_np2m/=0 .or. idiag_npmax/=0 .or. &
          idiag_npmin/=0 .or. idiag_npmx/=0 .or. idiag_npmy/=0 .or. &
          idiag_npmz/=0 .or. idiag_nparpmax/=0) lpenc_diagnos(i_np)=.true.
      if (idiag_rhopm/=0 .or. idiag_rhoprms/=0 .or. idiag_rhop2m/=0 .or. &
          idiag_rhopmax/=0 .or. idiag_rhopmin/=0 .or. idiag_rhopmphi/=0 .or. &
          idiag_rhopmx/=0 .or. idiag_rhopmy/=0 .or. idiag_rhopmz/=0) &
          lpenc_diagnos(i_rhop)=.true.
      if (idiag_dedragp/=0 .or. idiag_decollp/=0) then
        lpenc_diagnos(i_TT1)=.true.
        lpenc_diagnos(i_rho1)=.true.
      endif
      if (idiag_epspmx/=0 .or. idiag_epspmy/=0 .or. idiag_epspmz/=0 .or. &
          idiag_epspmin/=0 .or. idiag_epspmax/=0) &
          lpenc_diagnos(i_epsp)=.true.
      if (idiag_rhopmxy/=0 .or. idiag_rhopmxz/=0 .or. idiag_rhopmphi/=0) &
          lpenc_diagnos2d(i_rhop)=.true.
      if (idiag_npmxy/=0 ) lpenc_diagnos2d(i_np)=.true.
!
    endsubroutine pencil_criteria_particles
!***********************************************************************
    subroutine pencil_interdep_particles(lpencil_in)
!
!  Interdependency among pencils provided by the Particles module
!  is specified here.
!
!  16-feb-06/anders: dummy
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_rhop) .and. irhop==0) then
        lpencil_in(i_np)=.true.
      endif
!
      if (lpencil_in(i_epsp)) then
        lpencil_in(i_rhop)=.true.
        lpencil_in(i_rho1)=.true.
      endif
!
    endsubroutine pencil_interdep_particles
!***********************************************************************
    subroutine calc_pencils_particles(f,p)
!
      use Sub, only: grad
!
!  Calculate Particles pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  16-feb-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lpencil(i_np)) then
        if (inp/=0) then
          p%np=f(l1:l2,m,n,inp)
        else
          p%np=0.0
        endif
      endif
!
      if (lpencil(i_rhop)) then
        if (irhop/=0) then
          p%rhop=f(l1:l2,m,n,irhop)
        else
          p%rhop=rhop_swarm*f(l1:l2,m,n,inp)
        endif
      endif
!
      if (lpencil(i_grhop)) then
        if (irhop/=0) then
          if ((nprocx/=1).and.(.not.lcommunicate_rhop)) & 
               call fatal_error("calc_pencils_particles",&
               "Switch on lcommunicate_rhop=T in particles_run_pars")  
          call grad(f,irhop,p%grhop)
        else
          if ((nprocx/=1).and.(.not.lcommunicate_np)) & 
               call fatal_error("calc_pencils_particles",&
               "Switch on lcommunicate_np=T in particles_run_pars")  
          call grad(f,inp,p%grhop)
          p%grhop=rhop_swarm*p%grhop
        endif
      endif
!
      if (lpencil(i_epsp)) p%epsp=p%rhop*p%rho1
!
    endsubroutine calc_pencils_particles
!***********************************************************************
    subroutine dxxp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle position.
!
!  02-jan-05/anders: coded
!
      use General, only: random_number_wrapper, random_seed_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
!
!  Identify module and boundary conditions.
!
      if (lheader) print*,'dxxp_dt: Calculate dxxp_dt'
      if (lheader) then
        print*, 'dxxp_dt: Particles boundary condition bcpx=', bcpx
        print*, 'dxxp_dt: Particles boundary condition bcpy=', bcpy
        print*, 'dxxp_dt: Particles boundary condition bcpz=', bcpz
      endif
!
      if (lheader) print*, 'dxxp_dt: Set rate of change of particle '// &
          'position equal to particle velocity.'
!
!  The rate of change of a particle's position is the particle's velocity.
!
      if (lcartesian_coords) then
!
        if (nxgrid/=1) &
            dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
        if (nygrid/=1) &
            dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + fp(1:npar_loc,ivpy)
        if (nzgrid/=1) &
            dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
      elseif (lcylindrical_coords) then
!
        if (nxgrid/=1) &
            dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
        if (nygrid/=1) &
            dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + &
            fp(1:npar_loc,ivpy)/max(fp(1:npar_loc,ixp),tini)
        if (nzgrid/=1) &
            dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
      elseif (lspherical_coords) then
!
        if (nxgrid/=1) &
            dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
        if (nygrid/=1) &
            dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + &
            fp(1:npar_loc,ivpy)/max(fp(1:npar_loc,ixp),tini)
        if (nzgrid/=1) &
            dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + &
            fp(1:npar_loc,ivpz)/(max(fp(1:npar_loc,ixp),tini)*&
            sin(fp(1:npar_loc,iyp)))
!
      endif
!
!  With shear there is an extra term due to the background shear flow.
!
      if (lshear.and.nygrid/=1) dfp(1:npar_loc,iyp) = &
          dfp(1:npar_loc,iyp) - qshear*Omega*fp(1:npar_loc,ixp)
!
      if (lfirstcall) lfirstcall=.false.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine dvvp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle velocity.
!
!  29-dec-04/anders: coded
!
      use Diagnostics
      use EquationOfState, only: cs20
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real :: Omega2
      integer :: npar_found
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
      if (lheader) then
        print*,'dvvp_dt: Calculate dvvp_dt'
      endif
!
!  Add Coriolis force from rotating coordinate frame.
!
      if (Omega/=0.) then
        if (lcoriolis_force_par) then
          if (lheader) print*,'dvvp_dt: Add Coriolis force; Omega=', Omega
          Omega2=2*Omega
          if (.not.lspherical_coords) then
            dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + &
                Omega2*fp(1:npar_loc,ivpy)
            dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) - &
                Omega2*fp(1:npar_loc,ivpx)
          else
            print*,'dvvp_dt: Coriolis force on the particles is '
            print*,'not yet implemented for spherical coordinates.'
            call fatal_error('dvvp_dt','')
          endif
        endif
!
!  Add centrifugal force.
!
        if (lcentrifugal_force_par) then
          if (lheader) print*,'dvvp_dt: Add Centrifugal force; Omega=', Omega
          if (lcartesian_coords) then
!
            dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + &
                Omega**2*fp(1:npar_loc,ixp)
!
            dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) + &
                Omega**2*fp(1:npar_loc,iyp)
!
          elseif (lcylindrical_coords) then
            dfp(1:npar_loc,ivpx) = &
                dfp(1:npar_loc,ivpx) + Omega**2*fp(1:npar_loc,ixp)
          else
            print*,'dvvp_dt: Centrifugal force on the particles is '
            print*,'not implemented for spherical coordinates.'
            call fatal_error('dvvp_dt','')
          endif
        endif
!
!  With shear there is an extra term due to the background shear flow.
!
        if (lshear) dfp(1:npar_loc,ivpy) = &
            dfp(1:npar_loc,ivpy) + qshear*Omega*fp(1:npar_loc,ivpx)
      endif
!
!  Add constant background pressure gradient beta=alpha*H0/r0, where alpha
!  comes from a global pressure gradient P = P0*(r/r0)^alpha.
!  (the term must be added to the dust equation of motion when measuring
!  velocities relative to the shear flow modified by the global pressure grad.)
!
      if (beta_dPdr_dust/=0.0 .and. t>=tstart_dragforce_par) then
        dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + cs20*beta_dPdr_dust_scaled
      endif
!
!  Gravity on the particles.
!  Gravity on particles is implemented only if lparticle_gravity is true which is the default.
!
      if (lparticle_gravity)  call particle_gravity(f,df,fp,dfp,ineargrid)
!
!  Diagnostic output
!
      if (ldiagnos) then
        if (idiag_nparmin/=0) call max_name(-npar_loc,idiag_nparmin,lneg=.true.)
        if (idiag_nparmax/=0) call max_name(+npar_loc,idiag_nparmax)
        if (idiag_nparpmax/=0) call max_name(maxval(npar_imn),idiag_nparpmax)
        if (idiag_xpm/=0)  call sum_par_name(fp(1:npar_loc,ixp),idiag_xpm)
        if (idiag_ypm/=0)  call sum_par_name(fp(1:npar_loc,iyp),idiag_ypm)
        if (idiag_zpm/=0)  call sum_par_name(fp(1:npar_loc,izp),idiag_zpm)
        if (idiag_xp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2,idiag_xp2m)
        if (idiag_yp2m/=0) call sum_par_name(fp(1:npar_loc,iyp)**2,idiag_yp2m)
        if (idiag_zp2m/=0) call sum_par_name(fp(1:npar_loc,izp)**2,idiag_zp2m)
        if (idiag_rpm/=0)  call sum_par_name(sqrt(fp(1:npar_loc,ixp)**2+ &
            fp(1:npar_loc,iyp)**2+fp(1:npar_loc,izp)**2),idiag_rpm)
        if (idiag_rp2m/=0) call sum_par_name(fp(1:npar_loc,ixp)**2+ &
            fp(1:npar_loc,iyp)**2+fp(1:npar_loc,izp)**2,idiag_rp2m)
        if (idiag_vpxm/=0) call sum_par_name(fp(1:npar_loc,ivpx),idiag_vpxm)
        if (idiag_vpym/=0) call sum_par_name(fp(1:npar_loc,ivpy),idiag_vpym)
        if (idiag_vpzm/=0) call sum_par_name(fp(1:npar_loc,ivpz),idiag_vpzm)
        if (idiag_vpxvpym/=0) call sum_par_name( &
            fp(1:npar_loc,ivpx)*fp(1:npar_loc,ivpy),idiag_vpxvpym)
        if (idiag_vpxvpzm/=0) call sum_par_name( &
            fp(1:npar_loc,ivpx)*fp(1:npar_loc,ivpz),idiag_vpxvpzm)
        if (idiag_vpyvpzm/=0) call sum_par_name( &
            fp(1:npar_loc,ivpy)*fp(1:npar_loc,ivpz),idiag_vpyvpzm)
        if (idiag_lpxm/=0) call sum_par_name( &
            fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpz)- &
            fp(1:npar_loc,izp)*fp(1:npar_loc,ivpy),idiag_lpxm)
        if (idiag_lpym/=0) call sum_par_name( &
            fp(1:npar_loc,izp)*fp(1:npar_loc,ivpx)- &
            fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpz),idiag_lpym)
        if (idiag_lpzm/=0) call sum_par_name( &
            fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpy)- &
            fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpx),idiag_lpzm)
        if (idiag_lpx2m/=0) call sum_par_name( &
            (fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpz)- &
            fp(1:npar_loc,izp)*fp(1:npar_loc,ivpy))**2,idiag_lpx2m)
        if (idiag_lpy2m/=0) call sum_par_name( &
            (fp(1:npar_loc,izp)*fp(1:npar_loc,ivpx)- &
            fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpz))**2,idiag_lpy2m)
        if (idiag_lpz2m/=0) call sum_par_name( &
            (fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpy)- &
            fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpx))**2,idiag_lpz2m)
        if (idiag_vpx2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpx)**2,idiag_vpx2m)
        if (idiag_vpy2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpy)**2,idiag_vpy2m)
        if (idiag_vpz2m/=0) &
            call sum_par_name(fp(1:npar_loc,ivpz)**2,idiag_vpz2m)
        if (idiag_vprms/=0) &
            call sum_par_name((fp(1:npar_loc,ivpx)**2 &
                              +fp(1:npar_loc,ivpy)**2 &
                              +fp(1:npar_loc,ivpz)**2),idiag_vprms,lsqrt=.true.)
        if (idiag_vpyfull2m/=0) &
            call sum_par_name((fp(1:npar_loc,ivpy)- &
            qshear*Omega*fp(1:npar_loc,ixp))**2,idiag_vpyfull2m)
        if (idiag_ekinp/=0) then
          if (lparticles_mass) then
            call sum_par_name(0.5*fp(1:npar_loc,irhopswarm)* &
                sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
          else
            if (lcartesian_coords.and.(all(lequidist))) then
              call sum_par_name(0.5*rhop_swarm*npar_per_cell* &
                   sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
            else
              call sum_par_name(0.5*mp_swarm* &
                   sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
            endif
          endif
        endif
        if (idiag_epotpm/=0) call sum_par_name( &
            -gravr/sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_epotpm)
        if (idiag_vpmax/=0) call max_par_name( &
            sqrt(sum(fp(1:npar_loc,ivpx:ivpz)**2,2)),idiag_vpmax)
        if (idiag_vpxmax/=0) call max_par_name(fp(1:npar_loc,ivpx),idiag_vpxmax)
        if (idiag_vpymax/=0) call max_par_name(fp(1:npar_loc,ivpy),idiag_vpymax)
        if (idiag_vpzmax/=0) call max_par_name(fp(1:npar_loc,ivpz),idiag_vpzmax)
        if (idiag_eccpxm/=0) call sum_par_name( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,ixp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpx)/gravr-fp(1:npar_loc,ixp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_eccpxm)
        if (idiag_eccpym/=0) call sum_par_name( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,iyp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpy)/gravr-fp(1:npar_loc,iyp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_eccpym)
        if (idiag_eccpzm/=0) call sum_par_name( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,izp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpz)/gravr-fp(1:npar_loc,izp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_eccpzm)
        if (idiag_eccpx2m/=0) call sum_par_name(( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,ixp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpx)/gravr-fp(1:npar_loc,ixp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)))**2,idiag_eccpx2m)
        if (idiag_eccpy2m/=0) call sum_par_name(( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,iyp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpy)/gravr-fp(1:npar_loc,iyp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)))**2,idiag_eccpy2m)
        if (idiag_eccpz2m/=0) call sum_par_name(( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,izp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpz)/gravr-fp(1:npar_loc,izp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)))**2,idiag_eccpz2m)
        if (idiag_rhopvpxm/=0) then
          if (lparticles_mass) then
            call sum_par_name(fp(1:npar_loc,irhopswarm)*fp(1:npar_loc,ivpx), &
                idiag_rhopvpxm)
          elseif (lparticles_radius.and.lparticles_number) then
            call sum_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpx),idiag_rhopvpxm)
          endif
        endif
        if (idiag_rhopvpym/=0) then
          if (lparticles_mass) then
            call sum_par_name(fp(1:npar_loc,irhopswarm)*fp(1:npar_loc,ivpy), &
                idiag_rhopvpym)
          elseif (lparticles_radius.and.lparticles_number) then
            call sum_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpy),idiag_rhopvpym)
          endif
        endif
        if (idiag_rhopvpzm/=0) then
          if (lparticles_mass) then
            call sum_par_name(fp(1:npar_loc,irhopswarm)*fp(1:npar_loc,ivpz), &
                idiag_rhopvpzm)
          elseif (lparticles_radius.and.lparticles_number) then
            call sum_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpz),idiag_rhopvpzm)
          endif
        endif
        if (idiag_rhopvpxt/=0) then
          if (lparticles_mass) then
            call integrate_par_name(fp(1:npar_loc,irhopswarm)* &
                fp(1:npar_loc,ivpx),idiag_rhopvpxt)
          elseif (lparticles_radius.and.lparticles_number) then
            call integrate_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpx),idiag_rhopvpxt)
          endif
        endif
        if (idiag_rhopvpyt/=0) then
          if (lparticles_mass) then
            call integrate_par_name(fp(1:npar_loc,irhopswarm)* &
                fp(1:npar_loc,ivpy),idiag_rhopvpyt)
          elseif (lparticles_radius.and.lparticles_number) then
            call integrate_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpy),idiag_rhopvpyt)
          endif
        endif
        if (idiag_rhopvpzt/=0) then
          if (lparticles_mass) then
            call integrate_par_name(fp(1:npar_loc,irhopswarm)* &
                fp(1:npar_loc,ivpz),idiag_rhopvpzt)
          elseif (lparticles_radius.and.lparticles_number) then
            call integrate_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpz),idiag_rhopvpzt)
          endif
        endif
        if (idiag_rhopvpysm/=0) then
          if (lparticles_mass) then
            call sum_par_name(fp(1:npar_loc,irhopswarm)* &
                Sshear*fp(1:npar_loc,ixp),idiag_rhopvpysm)
          elseif (lparticles_radius.and.lparticles_number) then
            call sum_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                Sshear*fp(1:npar_loc,ixp),idiag_rhopvpysm)
          endif
        endif
        if (idiag_mpt/=0) then
          if (lparticles_mass) then
            call integrate_par_name((/fp(1:npar_loc,irhopswarm)/),idiag_mpt)
          elseif (lparticles_radius.and.lparticles_number) then
            call integrate_par_name((/four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)/),idiag_mpt)
          endif
        endif
        if (idiag_npargone/=0) then
          call count_particles(ipar,npar_found)
          call save_name(float(npar-npar_found),idiag_npargone)
        endif
        if (idiag_deshearbcsm/=0) then
          call sum_name(energy_gain_shear_bcs/npar,idiag_deshearbcsm)
        endif
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dvvp_dt
!***********************************************************************
    subroutine particle_gravity(f,df,fp,dfp,ineargrid)
!
!  Contribution of gravity to dvvp_dt
!
!  11-oct-12/dhruba: copied from dvvp_dt
!
      use Diagnostics
      use EquationOfState, only: cs20
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension(3) :: ggp
      real :: rr=0, vv=0, OO2
      integer :: k
      logical :: lheader, lfirstcall=.true.
!
      intent (in) :: f, fp, ineargrid
      intent (inout) :: df, dfp
!
      call keep_compiler_quiet(f,df)
!
!  Print out header information in first time step.
!
      lheader=lfirstcall .and. lroot
      if (lheader) then
        print*,'dvvp_dt: Calculating gravity'
      endif
!
!  Gravity in the x-direction.
!
      if (t>=tstart_grav_x_par) then
!
        select case (gravx_profile)
!
          case ('')
            if (lheader) print*, 'dvvp_dt: No gravity in x-direction.'
!
          case ('zero')
            if (lheader) print*, 'dvvp_dt: No gravity in x-direction.'
!
          case ('const','plain')
            if (lheader) print*, 'dvvp_dt: Constant gravity field in x-direction'
            dfp(1:npar_loc,ivpx)=dfp(1:npar_loc,ivpx) + gravx
!
          case ('linear')
            if (lheader) print*, 'dvvp_dt: Linear gravity field in x-direction.'
            dfp(1:npar_loc,ivpx)=dfp(1:npar_loc,ivpx) - &
                nu_epicycle2*fp(1:npar_loc,ixp)
!
          case ('sinusoidal')
            if (lheader) &
                print*, 'dvvp_dt: Sinusoidal gravity field in x-direction.'
            dfp(1:npar_loc,ivpx)=dfp(1:npar_loc,ivpx) + &
                gravx*sin(kx_gg*fp(1:npar_loc,ixp))
!
          case default
            call fatal_error('dvvp_dt','chosen gravx_profile is not valid!')
!
        endselect
!
      endif
!
!  Gravity in the z-direction.
!
      if (t>=tstart_grav_z_par) then
!
        select case (gravz_profile)
!
          case ('')
            if (lheader) print*, 'dvvp_dt: No gravity in z-direction.'
!
          case ('zero')
            if (lheader) print*, 'dvvp_dt: No gravity in z-direction.'
!
          case ('const','plain')
            if (lheader) print*, 'dvvp_dt: Constant gravity field in z-direction.'
            dfp(1:npar_loc,ivpz)=dfp(1:npar_loc,ivpz) + gravz
!
          case ('linear')
            if (lheader) print*, 'dvvp_dt: Linear gravity field in z-direction.'
            dfp(1:npar_loc,ivpz)=dfp(1:npar_loc,ivpz) - &
                nu_epicycle2*fp(1:npar_loc,izp)
!
          case ('sinusoidal')
            if (lheader) &
                print*, 'dvvp_dt: Sinusoidal gravity field in z-direction.'
            dfp(1:npar_loc,ivpz)=dfp(1:npar_loc,ivpz) + &
                gravz*sin(kz_gg*fp(1:npar_loc,izp))
!
          case default
            call fatal_error('dvvp_dt','chosen gravz_profile is not valid!')
!
        endselect
!
      endif
!
!  Radial gravity.
!
      if (t>=tstart_grav_r_par) then
!
        select case (gravr_profile)
!
        case ('')
          if (lheader) print*, 'dvvp_dt: No radial gravity'
!
        case ('zero')
          if (lheader) print*, 'dvvp_dt: No radial gravity'
!
        case ('newtonian-central','newtonian')
          if (lparticles_nbody) &
              call fatal_error('dvvp_dt','You are using massive particles. '//&
              'The N-body code should take care of the stellar-like '// &
              'gravity on the dust. Switch off the '// &
              'gravr_profile=''newtonian'' on particles_init')
          if (lheader) &
               print*, 'dvvp_dt: Newtonian gravity from a fixed central object'
          do k=1,npar_loc
            if (lcartesian_coords) then
              if (lcylindrical_gravity_par) then
                rr=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+gravsmooth2)
              else
                rr=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2+gravsmooth2)
              endif
              OO2=rr**(-3.)*gravr
              ggp(1) = -fp(k,ixp)*OO2
              ggp(2) = -fp(k,iyp)*OO2
              if (lcylindrical_gravity_par) then
                ggp(3) = 0.
              else
                ggp(3) = -fp(k,izp)*OO2
              endif
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
            elseif (lcylindrical_coords) then
              if (lcylindrical_gravity_par) then 
                rr=sqrt(fp(k,ixp)**2+gravsmooth2)
              else
                rr=sqrt(fp(k,ixp)**2+fp(k,izp)**2+gravsmooth2)
              endif
              OO2=rr**(-3.)*gravr
              ggp(1) = -fp(k,ixp)*OO2
              ggp(2) = 0.0
              if (lcylindrical_gravity_par) then
                ggp(3) = 0.
              else
                ggp(3) = -fp(k,izp)*OO2
              endif
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
            elseif (lspherical_coords) then
              rr=sqrt(fp(k,ixp)**2+gravsmooth2)
              OO2=rr**(-3)*gravr
              ggp(1) = -fp(k,ixp)*OO2
              ggp(2) = 0.0; ggp(3) = 0.0
              if (lcylindrical_gravity_par) call fatal_error("dvvp_dt",&
                   "No cylindrical gravity in spherical coordinates.")
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
            endif
!  Limit time-step if particles close to gravity source.
            if (ldt_grav_par.and.(lfirst.and.ldt)) then
              if (lcartesian_coords) then
                vv=sqrt(fp(k,ivpx)**2+fp(k,ivpy)**2+fp(k,ivpz)**2)
              elseif (lcylindrical_coords) then
                vv=sqrt(fp(k,ivpx)**2+fp(k,ivpz)**2)
              elseif (lspherical_coords) then
                vv=abs(fp(k,ivpx))
              endif
              dt1_max(ineargrid(k,1)-nghost)= &
                  max(dt1_max(ineargrid(k,1)-nghost),vv/rr/cdtpgrav)
            endif
          enddo
!
        case default
          call fatal_error('dvvp_dt','chosen gravr_profile is not valid!')
!
        endselect
!
      endif
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine particle_gravity
!***********************************************************************
    subroutine dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle position (called from main pencil loop).
!
!  25-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      integer :: k, ix0, iy0, iz0
      real :: dt1_advpx, dt1_advpy, dt1_advpz
      logical :: lnbody
!
!  Contribution of dust particles to time step.
!
      if (lfirst.and.ldt.and.ldt_adv_par) then
        if (npar_imn(imn)/=0) then
          do k=k1_imn(imn),k2_imn(imn)
            lnbody=(lparticles_nbody.and.any(ipar(k)==ipar_nbody))
            if (.not.lnbody) then
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              dt1_advpx=abs(fp(k,ivpx))*dx_1(ix0)
              if (lshear) then
                dt1_advpy=(-qshear*Omega*fp(k,ixp)+abs(fp(k,ivpy)))*dy_1(iy0)
              else
                dt1_advpy=abs(fp(k,ivpy))*dy_1(iy0)
              endif
              dt1_advpz=abs(fp(k,ivpz))*dz_1(iz0)
              if (l_shell) then
                dt1_advpx=abs(fp(k,ivpx))/k_shell
                dt1_advpy=abs(fp(k,ivpy))/k_shell
                dt1_advpz=abs(fp(k,ivpz))/k_shell
              endif
              dt1_max(ix0-nghost)=max(dt1_max(ix0-nghost), &
                   sqrt(dt1_advpx**2+dt1_advpy**2+dt1_advpz**2)/cdtp)
            endif
          enddo
        endif
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
!
    endsubroutine dxxp_dt_pencil
!***********************************************************************
    subroutine dvvp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of dust particle velocity (called from main pencil loop).
!
!  25-apr-06/anders: coded
!
      use Diagnostics
      use EquationOfState, only: cs20
      use Particles_spin, only: calc_liftforce
      use Particles_diagnos_dv, only: collisions
      use Particles_diagnos_state, only: persistence_check
      use Particles_dragforce
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx) :: dt1_drag, dt1_drag_gas, dt1_drag_dust
      real, dimension (nx) :: drag_heat
      real, dimension (3) :: dragforce, liftforce, bforce,thermforce, uup
      real, dimension(:), allocatable :: rep,stocunn
      real :: rho_point, rho1_point, tausp1_par, up2
      real :: weight, weight_x, weight_y, weight_z
      real :: rhop_swarm_par
      integer :: k, l, ix0, iy0, iz0
      integer :: ixx, iyy, izz, ixx0, iyy0, izz0, ixx1, iyy1, izz1
      logical :: lnbody, lsink
!
      intent (inout) :: f, df, dfp, fp, ineargrid
!
!  Identify module.
!
      if (headtt) then
        if (lroot) print*,'dvvp_dt_pencil: calculate dvvp_dt'
      endif
!
!  Precalculate certain quantities, if necessary.
!
      if (npar_imn(imn)/=0) then
!
!  Precalculate particle Reynolds numbers.
!
        if (ldraglaw_steadystate.or.lparticles_spin) then
          allocate(rep(k1_imn(imn):k2_imn(imn)))
!
          if (.not.allocated(rep)) then
            call fatal_error('dvvp_dt_pencil','unable to allocate sufficient'//&
              ' memory for rep')
          endif
          if (.not. interp%luu) then
            call fatal_error('dvvp_dt_pencil','you must set lnostore_uu=F'//&
              ' when rep is to be calculated')
          endif
!
          call calc_pencil_rep(fp,rep)
        endif
!
!  Precalculate Stokes-Cunningham factor (only if not ldraglaw_simple)
!
        if (.not.ldraglaw_simple) then
          if (ldraglaw_steadystate.or.lbrownian_forces) then
            allocate(stocunn(k1_imn(imn):k2_imn(imn)))
            if (.not.allocated(stocunn)) then
              call fatal_error('dvvp_dt_pencil','unable to allocate sufficient'//&
                'memory for stocunn')
            endif
!
            call calc_stokes_cunningham(fp,stocunn)
          endif
        endif
      endif
!
!  Drag force on particles and on gas.
!
      if (ldragforce_dust_par .and. t>=tstart_dragforce_par) then
        if (headtt) print*,'dvvp_dt: Add drag force; tausp=', tausp
!
        if (ldragforce_heat.or.(ldiagnos.and.idiag_dedragp/=0)) drag_heat=0.0
!
        if (npar_imn(imn)/=0) then
!
          if (lfirst.and.ldt) then
            dt1_drag_dust=0.0
            if (ldragforce_gas_par) dt1_drag_gas=0.0
          endif
!
!  Loop over all particles in current pencil.
!
          do k=k1_imn(imn),k2_imn(imn)
!
!  Exclude massive nbody particles and sink particles from the drag.
!
            lnbody=any(ipar(k)==ipar_nbody)
            lsink =.false.
            if (lparticles_sink) then
              if (fp(k,iaps)>0.0) lsink=.true.
            endif
            if ((.not.lnbody).and.(.not.lsink)) then
              ix0=ineargrid(k,1)
              iy0=ineargrid(k,2)
              iz0=ineargrid(k,3)
!
!  The interpolated gas velocity is either precalculated, and stored in
!  interp_uu, or it must be calculated here.
!
              if (.not.interp%luu) then
                if (lhydro) then
                  if (lparticlemesh_cic) then
                    call interpolate_linear(f,iux,iuz, &
                      fp(k,ixp:izp),uup,ineargrid(k,:),0,ipar(k))
                  elseif (lparticlemesh_tsc) then
                    if (linterpolate_spline) then
                      call interpolate_quadratic_spline(f,iux,iuz, &
                        fp(k,ixp:izp),uup,ineargrid(k,:),0,ipar(k))
                    else
                      call interpolate_quadratic(f,iux,iuz, &
                        fp(k,ixp:izp),uup,ineargrid(k,:),0,ipar(k))
                    endif
                  else
                    uup=f(ix0,iy0,iz0,iux:iuz)
                  endif
                elseif (l_shell) then
                  call calc_gas_velocity_shell_call(k,uup,fp)
                else
                  uup=0.0
                endif
              else
                uup=interp_uu(k,:)
              endif
!
!  Track particle state in terms of local gas velocity
!
              if (lparticles_diagnos_state .and. lfirst) &
                  call persistence_check(fp, k, uup)
!
!  Get the friction time. For the case of |uup| ~> cs, the Epstein drag law
!  is dependent on the relative mach number, hence the need to feed uup as
!  an optional argument to get_frictiontime.
!
              if (ldraglaw_epstein_transonic .or. &
                  ldraglaw_eps_stk_transonic) then
                call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par,uup)
              elseif (ldraglaw_steadystate) then
                call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par,rep=rep(k),&
                  stocunn=stocunn(k))
              else
                call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par)
              endif
!
!  Calculate and add drag force.
!
              dragforce = -tausp1_par*(fp(k,ivpx:ivpz)-uup)
!
!  Consider only (spherical) radial component of drag force (for testing).
!
              if (ldragforce_radialonly) then
                dragforce=fp(k,ixp:izp)*(dragforce(1)*fp(k,ixp)+ &
                    dragforce(2)*fp(k,iyp)+dragforce(3)*fp(k,izp))/ &
                    (fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2)
              endif
!
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + dragforce
!
!  Back-reaction friction force from particles on gas. Three methods are
!  implemented for assigning a particle to the mesh (see Hockney & Eastwood):
!
!    0. NGP (Nearest Grid Point)
!       The entire effect of the particle goes to the nearest grid point.
!    1. CIC (Cloud In Cell)
!       The particle has a region of influence with the size of a grid cell.
!       This is equivalent to a first order (spline) interpolation scheme.
!    2. TSC (Triangular Shaped Cloud)
!       The particle is spread over a length of two grid cells, but with
!       a density that falls linearly outwards.
!       This is equivalent to a second order spline interpolation scheme.
!
              if (ldragforce_gas_par) then
!
!  Cloud In Cell (CIC) scheme.
!
                if (lparticlemesh_cic) then
                  ixx0=ix0; iyy0=iy0; izz0=iz0
                  ixx1=ix0; iyy1=iy0; izz1=iz0
!
!  Particle influences the 8 surrounding grid points. The reference point is
!  the grid point at the lower left corner.
!
                  if ( (x(ix0)>fp(k,ixp)) .and. nxgrid/=1) ixx0=ixx0-1
                  if ( (y(iy0)>fp(k,iyp)) .and. nygrid/=1) iyy0=iyy0-1
                  if ( (z(iz0)>fp(k,izp)) .and. nzgrid/=1) izz0=izz0-1
                  if (nxgrid/=1) ixx1=ixx0+1
                  if (nygrid/=1) iyy1=iyy0+1
                  if (nzgrid/=1) izz1=izz0+1
                  do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
                    weight=1.0
                    if (nxgrid/=1) &
                         weight=weight*( 1.0-abs(fp(k,ixp)-x(ixx))*dx_1(ixx) )
                    if (nygrid/=1) &
                         weight=weight*( 1.0-abs(fp(k,iyp)-y(iyy))*dy_1(iyy) )
                    if (nzgrid/=1) &
                         weight=weight*( 1.0-abs(fp(k,izp)-z(izz))*dz_1(izz) )
!  Save the calculation of rho1 when inside pencil.
                    if ( (iyy/=m).or.(izz/=n).or.(ixx<l1).or.(ixx>l2) ) then
                      rho_point=f(ixx,iyy,izz,ilnrho)
                      if (.not. ldensity_nolog) rho_point=exp(rho_point)
                      rho1_point=1/rho_point
                    else
                      rho1_point=p%rho1(ixx-nghost)
                    endif
!  Add friction force to grid point.
                    call get_rhopswarm(mp_swarm,fp,k,ixx,iyy,izz,rhop_swarm_par)
                    df(ixx,iyy,izz,iux:iuz)=df(ixx,iyy,izz,iux:iuz) - &
                         rhop_swarm_par*rho1_point*dragforce*weight
                  enddo; enddo; enddo
!
!  Triangular Shaped Cloud (TSC) scheme.
!
                elseif (lparticlemesh_tsc) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
                  if (nxgrid/=1) then
                    ixx0=ix0-1; ixx1=ix0+1
                  else
                    ixx0=ix0  ; ixx1=ix0
                  endif
                  if (nygrid/=1) then
                    iyy0=iy0-1; iyy1=iy0+1
                  else
                    iyy0=iy0  ; iyy1=iy0
                  endif
                  if (nzgrid/=1) then
                    izz0=iz0-1; izz1=iz0+1
                  else
                    izz0=iz0  ; izz1=iz0
                  endif
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
                  do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
                    if ( ((ixx-ix0)==-1) .or. ((ixx-ix0)==+1) ) then
                      weight_x=1.125-1.5* abs(fp(k,ixp)-x(ixx))*dx_1(ixx) + &
                                     0.5*(abs(fp(k,ixp)-x(ixx))*dx_1(ixx))**2
                    else
                      if (nxgrid/=1) &
                           weight_x=0.75-((fp(k,ixp)-x(ixx))*dx_1(ixx))**2
                    endif
                    if ( ((iyy-iy0)==-1) .or. ((iyy-iy0)==+1) ) then
                      weight_y=1.125-1.5* abs(fp(k,iyp)-y(iyy))*dy_1(iyy) + &
                                     0.5*(abs(fp(k,iyp)-y(iyy))*dy_1(iyy))**2
                    else
                      if (nygrid/=1) &
                           weight_y=0.75-((fp(k,iyp)-y(iyy))*dy_1(iyy))**2
                    endif
                    if ( ((izz-iz0)==-1) .or. ((izz-iz0)==+1) ) then
                      weight_z=1.125-1.5* abs(fp(k,izp)-z(izz))*dz_1(izz) + &
                                     0.5*(abs(fp(k,izp)-z(izz))*dz_1(izz))**2
                    else
                      if (nzgrid/=1) &
                           weight_z=0.75-((fp(k,izp)-z(izz))*dz_1(izz))**2
                    endif
!
                    weight=1.0
!
                    if (nxgrid/=1) weight=weight*weight_x
                    if (nygrid/=1) weight=weight*weight_y
                    if (nzgrid/=1) weight=weight*weight_z
!  Save the calculation of rho1 when inside pencil.
                    if ( (iyy/=m).or.(izz/=n).or.(ixx<l1).or.(ixx>l2) ) then
                      rho_point=f(ixx,iyy,izz,ilnrho)
                      if (.not. ldensity_nolog) rho_point=exp(rho_point)
                      rho1_point=1/rho_point
                    else
                      rho1_point=p%rho1(ixx-nghost)
                    endif
!  Add friction force to grid point.
                    call get_rhopswarm(mp_swarm,fp,k,ixx,iyy,izz,rhop_swarm_par)
                    df(ixx,iyy,izz,iux:iuz)=df(ixx,iyy,izz,iux:iuz) - &
                         rhop_swarm_par*rho1_point*dragforce*weight
                  enddo; enddo; enddo
                else
!
!  Nearest Grid Point (NGP) scheme.
!
                  l=ineargrid(k,1)
                  call get_rhopswarm(mp_swarm,fp,k,l,m,n,rhop_swarm_par)
                  df(l,m,n,iux:iuz) = df(l,m,n,iux:iuz) - &
                       rhop_swarm_par*p%rho1(l-nghost)*dragforce
                endif
              endif
!
!  Heating of gas due to drag force.
!
              if (ldragforce_heat .or. (ldiagnos .and. idiag_dedragp/=0)) then
                if (ldragforce_gas_par) then
                  up2=sum((fp(k,ivpx:ivpz)-uup)**2)
                else
                  up2=sum(fp(k,ivpx:ivpz)*(fp(k,ivpx:ivpz)-uup))
                endif
!
! WL: Check if this is right. drag_heat is of dimension nx, not a scalar
!
                call get_rhopswarm(mp_swarm,fp,k,ix0,iy0,iz0,rhop_swarm_par)
                drag_heat(ix0-nghost)=drag_heat(ix0-nghost) + &
                     rhop_swarm_par*tausp1_par*up2
              endif
!
!  The minimum friction time of particles in a grid cell sets the local friction
!  time-step when there is only drag force on the dust,
!    dt1_drag = max(1/tausp)
!
!  With drag force on the gas as well, the maximum time-step is set as
!    dt1_drag = Sum_k[eps_k/tau_k]
!
              if (lfirst.and.ldt) then
                dt1_drag_dust(ix0-nghost)= &
                     max(dt1_drag_dust(ix0-nghost),tausp1_par)
                if (ldragforce_gas_par) then
                  if (p%np(ix0-nghost)/=0.0) &
                       dt1_drag_gas(ix0-nghost)=dt1_drag_gas(ix0-nghost)+ &
                       rhop_swarm_par*tausp1_par
                endif
              endif
            endif
          enddo
!
!  Add drag force heating in pencils.
!
          if (lentropy .and. ldragforce_heat) &
              df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%rho1*p%TT1*drag_heat
!
!  Contribution of friction force to time-step. Dust and gas inverse friction
!  time-steps are added up to give a valid expression even when the two are
!  of similar magnitude.
!
          if (lfirst.and.ldt) then
            if (ldragforce_gas_par) then
              dt1_drag=dt1_drag_dust+dt1_drag_gas
            else
              dt1_drag=dt1_drag_dust
            endif
            dt1_drag=dt1_drag/cdtp
            dt1_max=max(dt1_max,dt1_drag)
          endif
        else
!
!  No particles in this pencil.
!
          if (lfirst.and.ldt) dt1_drag=0.0
        endif
      endif
!
!  Add friction force from gas that moves systematically slower. Can be used
!  to mimic e.g. a sub-Keplerian gas flow without using the Hydro module.
!
      if (Deltauy_gas_friction/=0.0 .and. t>=tstart_dragforce_par) then
        if (npar_imn(imn)/=0) then
          do k=k1_imn(imn),k2_imn(imn)
            call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par)
            dfp(k,ivpy) = dfp(k,ivpy) - Deltauy_gas_friction*tausp1_par
          enddo
        endif
      endif
!
!  Collisional cooling is in a separate subroutine.
!
      if ( (lcollisional_cooling_taucool .or. lcollisional_cooling_rms .or. &
          lcollisional_cooling_twobody .or. lcollisional_dragforce_cooling) &
          .and. t>=tstart_collisional_cooling) &
          call collisional_cooling(f,df,fp,dfp,p,ineargrid)
!
!  Compensate for increased friction time by appying extra friction force to
!  particles.
!
      if (lcompensate_friction_increase) &
          call compensate_friction_increase(f,fp,dfp,p,ineargrid)
!
!  Add lift forces.
!
      if (lparticles_spin .and. t>=tstart_liftforce_par) then
        if (npar_imn(imn)/=0) then
          do k=k1_imn(imn),k2_imn(imn)
            call calc_liftforce(fp(k,:), k, rep(k), liftforce)
            dfp(k,ivpx:ivpz)=dfp(k,ivpx:ivpz)+liftforce
          enddo
        endif
      endif
!
!  Add Brownian forces.
!
      if (lbrownian_forces .and. t>=tstart_brownian_par) then
        if (npar_imn(imn)/=0) then
          do k=k1_imn(imn),k2_imn(imn)
            call calc_brownian_force(fp,k,ineargrid(k,:),stocunn(k),bforce)
            dfp(k,ivpx:ivpz)=dfp(k,ivpx:ivpz)+bforce
          enddo
        endif
      endif
!
!  Add Thernophoretic forces.
!
      if (lthermophoretic_forces) then
        if (npar_imn(imn)/=0) then
          do k=k1_imn(imn),k2_imn(imn)
            call calc_thermophoretic_force(fp,k,ineargrid(k,:),thermforce)
            dfp(k,ivpx:ivpz)=dfp(k,ivpx:ivpz)+thermforce
          enddo
        endif
      endif
!
!  For stiff drag force equations we need to store the forces that are
!  unique to the gas.
!
      if (ldragforce_stiff .and. .not. lpencil_check_at_work) then
        f(l1:l2,m,n,ifgx:ifgz)=p%fpres+p%jxbr+p%fvisc
      endif
!
!  Diagnostic output.
!
      if (ldiagnos) then
        if (idiag_npm/=0)      call sum_mn_name(p%np,idiag_npm)
        if (idiag_np2m/=0)     call sum_mn_name(p%np**2,idiag_np2m)
        if (idiag_npmax/=0)    call max_mn_name(p%np,idiag_npmax)
        if (idiag_npmin/=0)    call max_mn_name(-p%np,idiag_npmin,lneg=.true.)
        if (idiag_rhopm/=0)    call sum_mn_name(p%rhop,idiag_rhopm)
        if (idiag_rhop2m/=0 )  call sum_mn_name(p%rhop**2,idiag_rhop2m)
        if (idiag_rhoprms/=0)  call sum_mn_name(p%rhop**2,idiag_rhoprms,lsqrt=.true.)
        if (idiag_rhopmax/=0)  call max_mn_name(p%rhop,idiag_rhopmax)
        if (idiag_rhopmin/=0)  call max_mn_name(-p%rhop,idiag_rhopmin,lneg=.true.)
        if (idiag_epspmax/=0)  call max_mn_name(p%epsp,idiag_epspmax)
        if (idiag_epspmin/=0)  call max_mn_name(-p%epsp,idiag_epspmin,lneg=.true.)
        if (idiag_dedragp/=0)  call sum_mn_name(drag_heat,idiag_dedragp)
        if (idiag_dvpx2m/=0 .or. idiag_dvpx2m/=0 .or. idiag_dvpx2m/=0 .or. &
            idiag_dvpm  /=0 .or. idiag_dvpmax/=0) &
            call calculate_rms_speed(fp,ineargrid,p)
        if (idiag_dtdragp/=0.and.(lfirst.and.ldt))  &
            call max_mn_name(dt1_drag,idiag_dtdragp,l_dt=.true.)
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
        if (idiag_npmx/=0)    call yzsum_mn_name_x(p%np,idiag_npmx)
        if (idiag_npmy/=0)    call xzsum_mn_name_y(p%np,idiag_npmy)
        if (idiag_npmz/=0)    call xysum_mn_name_z(p%np,idiag_npmz)
        if (idiag_rhopmx/=0)  call yzsum_mn_name_x(p%rhop,idiag_rhopmx)
        if (idiag_rhopmy/=0)  call xzsum_mn_name_y(p%rhop,idiag_rhopmy)
        if (idiag_rhopmz/=0)  call xysum_mn_name_z(p%rhop,idiag_rhopmz)
        if (idiag_epspmx/=0)  call yzsum_mn_name_x(p%epsp,idiag_epspmx)
        if (idiag_epspmy/=0)  call xzsum_mn_name_y(p%epsp,idiag_epspmy)
        if (idiag_epspmz/=0)  call xysum_mn_name_z(p%epsp,idiag_epspmz)
        if (idiag_rhopmr/=0)  call phizsum_mn_name_r(p%rhop,idiag_rhopmr)
      endif
!
      if (l2davgfirst) then
        if (idiag_npmxy/=0)    call zsum_mn_name_xy(p%np,idiag_npmxy)
        if (idiag_rhopmphi/=0) call phisum_mn_name_rz(p%rhop,idiag_rhopmphi)
        if (idiag_rhopmxy/=0)  call zsum_mn_name_xy(p%rhop,idiag_rhopmxy)
        if (idiag_rhopmxz/=0)  call ysum_mn_name_xz(p%rhop,idiag_rhopmxz)
      endif
!
!  particle-particle separation and relative velocity diagnostics
!
      if (lparticles_diagnos_dv .and. lfirstpoint .and. lfirst) then
        if (t > t_nextcol) call collisions(fp)
      endif
!
!  Clean up (free allocated memory).
!
      if (allocated(rep)) deallocate(rep)
      if (allocated(stocunn)) deallocate(stocunn)
!
    endsubroutine dvvp_dt_pencil
!***********************************************************************
    subroutine dxxp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle position in blocks.
!
!  29-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dxxp_dt_blocks
!***********************************************************************
    subroutine dvvp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle velocity in blocks.
!
!  29-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dvvp_dt_blocks
!***********************************************************************
    subroutine remove_particles_sink_simple(f,fp,dfp,ineargrid)
!
!  Subroutine for taking particles out of the simulation due to their proximity
!  to a sink particle or sink point.
!
!  25-sep-08/anders: coded
!
      use Mpicomm
      use Solid_Cells
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(3) :: momp_swarm_removed, momp_swarm_removed_send
      real :: rp, rp_box, rhop_swarm_removed, rhop_swarm_removed_send
      real :: xsinkpar, ysinkpar, zsinkpar
      integer :: k, ksink, iproc_sink, iproc_sink_send
      integer :: ix, ix1, ix2, iy, iy1, iy2, iz, iz1, iz2
      integer, parameter :: itag1=100, itag2=101
      real :: particle_radius
!
!  Sinkparticle activated at time tstart_sink_par
      if (t <= tstart_sink_par) return
!
      if (lsinkpoint) then
        k=1
        do while (k<=npar_loc)
          rp=sqrt((fp(k,ixp)-xsinkpoint)**2+(fp(k,iyp)-ysinkpoint)**2+ &
             (fp(k,izp)-zsinkpoint)**2)
          if (rp<rsinkpoint) then
            call remove_particle(fp,ipar,k,dfp,ineargrid)
          else
            k=k+1
          endif
        enddo
      endif
!
      if (lsinkparticle_1) then
!
!  Search for position of particle with index 1.
!
        xsinkpoint=0.0
        ysinkpoint=0.0
        zsinkpoint=0.0
        iproc_sink=-1
        do k=1,npar_loc
          if (ipar(k)==1) then
            ksink=k
            iproc_sink=iproc
            xsinkpar=fp(k,ixp)
            ysinkpar=fp(k,iyp)
            zsinkpar=fp(k,izp)
          endif
        enddo
!
!  Communicate position of sink particle to all processors.
!
        iproc_sink_send=iproc_sink
        call mpireduce_max_int(iproc_sink_send,iproc_sink)
        call mpibcast_int(iproc_sink)
        call mpibcast_real(xsinkpar,proc=iproc_sink)
        call mpibcast_real(ysinkpar,proc=iproc_sink)
        call mpibcast_real(zsinkpar,proc=iproc_sink)
!
!  Remove particles that are too close to sink particle.
!
        if (lparticles_mass) then
          rhop_swarm_removed=0.0
          momp_swarm_removed=0.0
        endif
!
        k=1
        do while (k<=npar_loc)
          rp=-1.0
!
!  Take into account periodic boundary conditions.
!
          if (bcpx=='p' .and. nxgrid/=1) then
            ix1=-1; ix2=+1
          else
            ix1=0; ix2=0
          endif
          if (bcpy=='p' .and. nygrid/=1) then
            iy1=-1; iy2=+1
          else
            iy1=0; iy2=0
          endif
          if (bcpz=='p' .and. nzgrid/=1) then
            iz1=-1; iz2=+1
          else
            iz1=0; iz2=0
          endif
          do iz=iz1,iz2; do iy=iy1,iy2; do ix=ix1,ix2
            if (lshear) then
              rp_box=sqrt((fp(k,ixp)+Lx*ix-xsinkpar)**2+ &
                          (fp(k,iyp)+Ly*iy-deltay*ix-ysinkpar)**2+ &
                          (fp(k,izp)+Lz*iz-zsinkpar)**2)
            else
              rp_box=sqrt((fp(k,ixp)+Lx*ix-xsinkpar)**2+ &
                          (fp(k,iyp)+Ly*iy-ysinkpar)**2+ &
                          (fp(k,izp)+Lz*iz-zsinkpar)**2)
            endif
            if (rp<0.0) then
              rp=rp_box
            else
              rp=min(rp,rp_box)
            endif
          enddo; enddo; enddo
          if (ipar(k)/=1 .and. rp<rsinkparticle_1) then
            if (lparticles_mass) then
              rhop_swarm_removed = rhop_swarm_removed + fp(k,irhopswarm)
              momp_swarm_removed = momp_swarm_removed + &
                  fp(k,ivpx:ivpz)*fp(k,irhopswarm)
!              if (lshear) then
!                momp_swarm_removed(2) = momp_swarm_removed(2) + &
!                    Sshear*fp(k,ixp)*fp(k,irhopswarm)
!              endif
            endif
            call remove_particle(fp,ipar,k,dfp,ineargrid)
          else
            k=k+1
          endif
        enddo
!
!  Add mass and momentum to sink particle.
!
        if (lparticles_mass) then
          rhop_swarm_removed_send=rhop_swarm_removed
          momp_swarm_removed_send=momp_swarm_removed
          call mpireduce_sum(rhop_swarm_removed_send,rhop_swarm_removed)
          call mpireduce_sum(momp_swarm_removed_send,momp_swarm_removed,3)
          if (lroot) then
            call mpisend_real(rhop_swarm_removed,1,iproc_sink,itag1)
            call mpisend_real(momp_swarm_removed,3,iproc_sink,itag2)
          endif
          if (iproc==iproc_sink) then
            call mpirecv_real(rhop_swarm_removed,1,0,itag1)
            call mpirecv_real(momp_swarm_removed,3,0,itag2)
!
!  Need to find sink particle again since particle removal may have shifted
!  its index.
!
            do k=1,npar_loc
              if (ipar(k)==1) then
                ksink=k
                iproc_sink=iproc
                xsinkpar=fp(k,ixp)
                ysinkpar=fp(k,iyp)
                zsinkpar=fp(k,izp)
              endif
            enddo
!            if (lshear) then
!              fp(ksink,ivpx) = (fp(ksink,ivpx)* &
!                  fp(ksink,irhopswarm) + momp_swarm_removed(1))/ &
!                  (fp(ksink,irhopswarm) + rhop_swarm_removed)
!              fp(ksink,ivpy) = ((fp(ksink,ivpy)+Sshear*fp(ksink,ixp))* &
!                  fp(ksink,irhopswarm) + momp_swarm_removed(2))/ &
!                  (fp(ksink,irhopswarm) + rhop_swarm_removed) - &
!                  Sshear*fp(ksink,ixp)
!              fp(ksink,ivpz) = (fp(ksink,ivpz)* &
!                  fp(ksink,irhopswarm) + momp_swarm_removed(3))/ &
!                  (fp(ksink,irhopswarm) + rhop_swarm_removed)
!            else
              fp(ksink,ivpx:ivpz) = (fp(ksink,ivpx:ivpz)* &
                  fp(ksink,irhopswarm) + momp_swarm_removed)/ &
                  (fp(ksink,irhopswarm) + rhop_swarm_removed)
!            endif
            fp(ksink,irhopswarm) = fp(ksink,irhopswarm) + rhop_swarm_removed
          endif
        endif
      endif
!
!  Remove particles if they are within a solid geometry
!
      if (lsolid_cells) then
        k=1
        do while (k<=npar_loc)
          if (lparticles_radius) then
            particle_radius=fp(k,iap)
          else
            particle_radius=fake_particle_radius
          endif
          if (in_solid_cell(fp(k,ixp:izp),particle_radius)) then
            call remove_particle(fp,ipar,k,dfp,ineargrid)
          else
            k=k+1
          endif
        enddo
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine remove_particles_sink_simple
!***********************************************************************
    subroutine create_particles_sink_simple(f,fp,dfp,ineargrid)
!
!  Subroutine for creating new sink particles or sink points.
!
!  Just a dummy routine for now.
!
!  25-sep-08/anders: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mpvar) :: fp, dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine create_particles_sink_simple
!***********************************************************************
    subroutine get_frictiontime(f,fp,p,ineargrid,k,tausp1_par,uup,&
      nochange_opt,rep,stocunn)
!
!  Calculate the friction time.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      type (pencil_case) :: p
      real :: tausp1_par, tmp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: k
      logical, optional :: nochange_opt
!
      real, optional, dimension(3) :: uup
      real, optional :: rep, stocunn
!
      real :: tausg1_point,OO
      integer :: ix0, iy0, iz0, inx0, jspec
      logical :: nochange=.true.
!
      intent(in) :: rep,uup
!
      if (present(nochange_opt)) then
        if (nochange_opt) then
          nochange=.true.
        else
          nochange=.false.
        endif
      else
        nochange=.false.
      endif
!
      ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
      inx0=ix0-nghost
!
!  Epstein drag law.
!
      if (ldraglaw_epstein) then
        if (iap/=0) then
          if (fp(k,iap)/=0.0) tausp1_par = 1/(fp(k,iap)*rhopmat)
        else
!  Check if we are using multiple or single particle species.
          if (npar_species>1) then
            jspec=npar_species*(ipar(k)-1)/npar+1
            tmp=tausp1_species(jspec)
          else
            tmp=tausp1
          endif
!
!  Scale friction time with local density.
!
          if (ldraglaw_variable_density) then
            if (ldensity_nolog) then
              tausp1_par=tmp*f(ix0,iy0,iz0,irho)
            else
              tausp1_par=tmp*exp(f(ix0,iy0,iz0,ilnrho))
            endif
!
!  Discriminate between constant tau and special case for
!  1/tau=omega when omega is not constant (as, for instance,
!  global Keplerian disks, for which omega=rad**(-3/2)
!
          elseif (ldraglaw_variable) then
            if (lcartesian_coords) then
              OO=(fp(k,ixp)**2 + fp(k,iyp)**2)**(-0.75)
            elseif (lcylindrical_coords) then
              OO=fp(k,ixp)**(-1.5)
            elseif (lspherical_coords) then
              call fatal_error("get_frictiontime",&
                   "variable draglaw not implemented for"//&
                   "spherical coordinates")
              OO=0.
            else
              call fatal_error("get_frictiontime", &
                  "no valid coord system")
              OO=0.
            endif
            tausp1_par=tmp*OO
          else
!
!  Constant friction time.
!
            tausp1_par=tmp
          endif
        endif
      else if (ldraglaw_epstein_stokes_linear) then
!
!  When the particle radius is larger than 9/4 times the mean free path
!  of the gas molecules one must use the Stokes drag law rather than the
!  Epstein law.
!
!  We need here to know the mean free path of the gas molecules:
!    lambda = mu_mol/(rhog*sigma_mol)
!
!  The quantities are:
!    mu_mol    = mean molecular weight          [=3.9e-24 g for H_2]
!    rhog      = gas density
!    sigma_mol = cross section of gas molecules [=2e-15 cm^2 for H_2]
!
!  Actually need to know the mean free path in units of the gas scale
!  height H [if H=1]. Inserting the mid-plane expression
!    rhog=Sigmag/[sqrt(2*pi)*H]
!  gives
!    lambda/H = sqrt(2*pi)*mu_mol/(Sigmag*sigma_mol)
!            ~= 4.5e-9/Sigmag
!  when Sigmag is given in g/cm^2.
!
        if (iap==0) then
          if (lroot) print*, 'get_frictiontime: need particle radius as dynamical variable for Stokes law'
          call fatal_error('get_frictiontime','')
        endif
        if (fp(k,iap)<2.25*mean_free_path_gas) then
          tausp1_par = 1/(fp(k,iap)*rhopmat)
        else
          tausp1_par = 1/(fp(k,iap)*rhopmat)*2.25*mean_free_path_gas/fp(k,iap)
        endif
!
      else if (ldraglaw_epstein_transonic) then
!
! Draw laws for intermediate mach number. This is for pure Epstein drag...
!
        call calc_draglaw_parameters(fp,k,uup,p,inx0,tausp1_par)
!
      else if (ldraglaw_eps_stk_transonic) then
!
! ...and this is for a linear combination of Esptein and Stokes drag at
! intermediate mach number. Pure Stokes drag is not implemented.
!
        call calc_draglaw_parameters(fp,k,uup,p,inx0,tausp1_par,lstokes=.true.)
!
      elseif (ldraglaw_steadystate) then
        if (.not.present(rep)) then
          call fatal_error('get_frictiontime','need particle reynolds '// &
                  'number, rep, to calculate the steady state drag '// &
                  'relaxation time!')
        elseif (.not.present(stocunn)) then
          call fatal_error('get_frictiontime','need particle stokes '// &
                  'cunningham factor, stocunn, to calculate the steady '// &
                  ' state drag relaxation time!')
        else
          call calc_draglaw_steadystate(fp,k,rep,stocunn,tausp1_par)
        endif
!
! Simple drag law, drag force = 1/\taup (v-u) where \taup is an input parameter.
!
      else if (ldraglaw_simple) then
!        write(*,*)'DM','simple drag'
!  Check if we are using multiple or single particle species.
        if (npar_species>1) then
          jspec=npar_species*(ipar(k)-1)/npar+1
          tmp=tausp1_species(jspec)
        else
          tmp=tausp1
        endif
        tausp1_par=tmp
      endif
!
!  Change friction time artificially.
!
      if (.not. nochange) then
!
!  Increase friction time to avoid very small time-steps where the
!  dust-to-gas ratio is high.
!
        if (tausg_min/=0.0) then
          tausg1_point=tausp1_par*p%epsp(ix0-nghost)
          if (tausg1_point>tausg1_max) &
              tausp1_par=tausg1_max/p%epsp(ix0-nghost)
        endif
!
!  Increase friction time linearly with dust density where the dust-to-gas
!  ratio is higher than a chosen value. Supposed to mimick decreased cooling
!  when the gas follows the dust.
!
        if (epsp_friction_increase/=0.0) then
          if (p%epsp(ix0-nghost)>epsp_friction_increase) &
              tausp1_par=tausp1_par/(p%epsp(ix0-nghost)/epsp_friction_increase)
        endif
!
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine get_frictiontime
!***********************************************************************
    subroutine calc_draglaw_parameters(fp,k,uup,p,inx0,tausp1_par,lstokes)
!
      use EquationOfState, only: rho0,cs0
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension(3) :: uup,duu
      type (pencil_case) :: p
      real :: tausp1_par,tmp,tmp1
      integer :: k, inx0, jspec
      real :: kd,fd,mach,mach2,fac,OO
      real :: knudsen,reynolds,lambda
      real :: inv_particle_radius,kn_crit
      logical, optional :: lstokes
      logical, save :: lfirstcall
!
!  Epstein drag away from the limit of subsonic particle motion. The drag
!  force is given by (Schaaf 1963)
!
!       Feps=-pi*a**2 * rhog * |Delta(u)| * Delta(u) &                    (1)
!        *[(1+1/m**2+1/(4*m**4))*erf(m) + (1/m+1/(2*m**3)*exp(-m**2)/sqrt(pi)]
!
!  where Delta(u) is the relative dust-to-gas velocity (vector)
!  and m=|Delta(u)|/cs (scalar) is the relative mach number of the flow
!
!  As erf is too cumbersome a function to implement numerically, an interpolation
!  between the limits of
!
!     subsonic:    Feps=-sqrt(128*pi)/3*a**2*rhog*cs*Delta(u)             (2)
!     supersonic:  Feps=-pi*a**2*rhog*|Delta(u)|*Delta(u)                 (3)
!
!  is used, leading to an expression that can be used for arbitrary velocities
!  as derived by Kwok (1975).
!
!     transonic:  Feps=-sqrt(128*pi)/3*a**2*rhog*cs*fd*Delta(u)           (4)
!
!  where fd=sqrt(1 + 9*pi/128*m**2)                                       (5)
!
!  The force Feps is divided by the mass of the particle mp=4/3*pi*a**3*rhopmat
!  to yield the acceleration feps=Feps/mp
!
!         feps = -sqrt(8/pi)*rhog*cs*fd*Delta(u)/[a*rhopmat]              (6)
!
!  Epstein drag ceases to work when the particle diameter becomes comparable
!  to the mean free path (lambda) of the gas molecules. In this case, the force
!  is given by Stokes friction in the viscous case (low dust Reynolds number)
!
!      Fsto=-6*pi*a*mu_kin*Delta(u)                                       (7)
!
!  where mu_kin is the kinematic viscosity of the gas
!
!      mu_kin=1/3*rhog*vth*lambda                                         (8)
!
!  and vth=sqrt(8/pi)*cs is the mean thermal velocity of the gas. For high dust
!  Reynolds numbers the viscosity if uninmportant and the drag force of the tur-
!  bulent flow past the particle is given by Newtonian friction
!
!     Fnew=-1.3*pi*a**2*rhog*|Delta(u)|*Delta(u)
!
!  The two cases are once again connected by an interpolating factor
!
!     F'sto=-6*pi*a*kd*mu_kin*Delta(u)
!
!  where kd is a factor that contains the Reynolds number of the flow over the
!  particle (defined in the code, some lines below).
!
!  The following interpolation then works for flows of arbitrary Knudsen, Mach and Reynolds
!  numbers
!
!     Fdrag = [Kn'/(Kn'+1)]**2 * Feps +  [1/(Kn'+1)]**2 * F'sto
!
!  Where Kn'=3*Kn is the critical Knudsen number where the viscous (Stokes) drag and the subsonic
!  Epstein drag are equal.
!
!  (The discussion above was taken from Paardekooper 2006, Woite & Helling 2003 and Kwok 1975)
!
!  In the 2D case, the density rhog is to be replaced by
!
!     rhog=Sigmag/[sqrt(2*pi)H]
!         =Sigmag*Omega/[sqrt(2*pi)*cs]
!
!  which removes the dependence of (6) on cs. We are left with
!
!         feps = -2/pi*sigmag*Omega*fd*Delta(u)/[a*rhopmat]
!
!  the constant terms are tausp1. The same follows for Stokes drag
!
!  Friction time for different species
!
      if (npar_species==1) then
        tmp=tausp
        tmp1=tausp1
      else
        jspec=npar_species*(ipar(k)-1)/npar+1
        tmp=tausp_species(jspec)
        tmp1=tausp1_species(jspec)
      endif
!
!  Relative velocity
!
      duu=fp(k,ivpx:ivpz)-uup
!
      if (nzgrid==1) then
!  then omega is needed
        if (ldraglaw_variable) then
          !these omegas assume GM=1
          if (lcartesian_coords) then
            OO=(fp(k,ixp)**2 + fp(k,iyp)**2)**(-0.75)
          elseif (lcylindrical_coords) then
            OO=fp(k,ixp)**(-1.5)
          elseif (lspherical_coords) then
            call fatal_error('get_frictiontime', &
                'variable draglaw not implemented for '//&
                'spherical coordinates')
          endif
        else
          OO=nu_epicycle
        endif
      endif
!
!  Possibility to include the transition from Epstein to Stokes drag
!
      if (present(lstokes)) then
!
        if (lfirstcall) &
            print*, 'get_frictiontime: Epstein-Stokes transonic drag law'
!
!  The mach number and the correction fd to flows of arbitrary mach number
!
        mach=sqrt((duu(1)**2+duu(2)**2+duu(3)**2)/p%cs2(inx0))
        fd=sqrt(1+(9.0*pi/128)*mach**2)
!
!  For Stokes drag, the mean free path is needed
!
!   lambda = 1/rhog*(mu/sigma_coll)_H2
!
!  were mu_H2 is the mean molecular weight of the hydrogen molecule (3.9e-24 g),
!  and sigma_coll its cross section (2e-15 cm^2).
!  Assume that (mu/sigma_coll) is the input parameter mean_free_path_gas
!
        if (mean_free_path_gas==0) then
          print*,'You want to use Stokes drag but you forgot to set '//&
               'mean_free_path_gas in the .in files. Stop and check.'
          call fatal_error("calc_draglaw_parameters","")
        endif
!
        if (nzgrid==1) then
          !the sqrt(2pi) factor is inside the mean_free_path_gas constant
          lambda=mean_free_path_gas * sqrt(p%cs2(inx0))*rho0/(p%rho(inx0)*OO*cs0)
        else
          lambda=mean_free_path_gas * rho0/p%rho(inx0)
        endif
!
!  The Knudsen number is the ratio of the mean free path to the particle
!  radius, 2s. To keep consistency with the formulation evolving for radius,
!  tausp1 is C/(s*rhopmat) where C is 2/pi for 2d runs and sqrt(8/pi) for 3D
!  runs (because of the sqrt(2*pi) factor coming from the substitution
!  Sigma=rho/(sqrt(2*pi)*H). 's' is the particle radius
        if (iap/=0) then
          inv_particle_radius=1/fp(k,iap)
        else
          if (luse_tau_ap) then
            ! use tausp as the radius (in meter) to make life easier
            inv_particle_radius=tmp1
          else
            if (nzgrid==1) then
              inv_particle_radius=0.5*pi*tmp1     !rhopmat=1, particle_radius in meters
            else
              inv_particle_radius=sqrt(pi/8)*tmp1 !rhopmat=1, particle_radius in meters
            endif
          endif
        endif
!
        knudsen=0.5*lambda*inv_particle_radius
!
!  The Stokes drag depends non-linearly on
!
!    Re = 2*s*rho_g*|delta(v)|/mu_kin
!
        reynolds=3*sqrt(pi/8)*mach/knudsen
!
!  the Reynolds number of the flow over the particle. It can parameterized by
!
        if (reynolds<=500) then
          kd=1.0+0.15*reynolds**0.687
        elseif ((reynolds>500).and.(reynolds<=1500)) then
          kd=3.96e-6*reynolds**2.4
        elseif (reynolds>1500) then
          kd=0.11*reynolds
        else
          call fatal_error("calc_draglaw_parameters", &
              "something went pretty wrong")
          kd=0.
        endif
!
!  And we finally have the Stokes correction to intermediate Knudsen numbers
!  kn_crit is the critical knudsen number where viscous (low reynolds)
!  Stokes and subsonic Epstein friction are equal (Woitke & Helling, 2003)
!
        kn_crit=3*knudsen
        fac=kn_crit/(kn_crit+1)**2 * (kn_crit*fd + kd)
!
      else
!
!  Only use Epstein drag
!
        if (lfirstcall) &
            print*,'get_frictiontime: Epstein transonic drag law'
!
        mach2=(duu(1)**2+duu(2)**2+duu(3)**2)/p%cs2(inx0)
        fd=sqrt(1+(9.0*pi/128)*mach2)
        fac=fd
!
      endif
!
! Calculate tausp1_par for 2d and 3d cases with and without particle_radius
! as a dynamical variable
!
      if (iap/=0) then
        if (fp(k,iap)/=0.0) then
          if (nzgrid==1) then
            tausp1_par=     2*pi_1*OO* &
                p%rho(inx0)*fac/(fp(k,iap)*rhopmat)
          else
            tausp1_par=sqrt(8*pi_1*p%cs2(inx0))*p%rho(inx0)* &
                fac/(fp(k,iap)*rhopmat)
          endif
        endif
      else
          !normalize to make tausp1 not dependent on cs0 or rho0
          !bad because it comes at the expense of evil divisions
        if (nzgrid==1) then
          if (luse_tau_ap) then
            tausp1_par=tmp1*2*pi_1*OO*p%rho(inx0)*fac/(rho0*rhopmat)
          else
            tausp1_par=tmp1*OO*p%rho(inx0)*fac/ rho0
          endif
        else
          if (luse_tau_ap) then
            tausp1_par=tmp1*sqrt(8*pi_1*p%cs2(inx0))*p%rho(inx0)*fac/(rho0*cs0)
          else
            tausp1_par=tmp1*sqrt(p%cs2(inx0))*p%rho(inx0)*fac/(rho0*cs0)
          endif
        endif
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine calc_draglaw_parameters
!***********************************************************************
    subroutine collisional_cooling(f,df,fp,dfp,p,ineargrid)
!
!  Reduce relative speed between particles due to inelastic collisions.
!
!  23-sep-06/anders: coded
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx,npar_species,npar_species) :: tau_coll1_species
      real, dimension (nx,3,npar_species) :: vvpm_species
      real, dimension (nx,npar_species) :: np_species, vpm_species
      real, dimension (nx,npar_species) :: tau_coll1_tot
      real, dimension (nx,3) :: vvpm
      real, dimension (nx) :: vpm, tau_coll1, tausp1m, vcoll, coll_heat
      real, dimension (nx) :: rhop_swarm_mn
      real, dimension (3) :: deltavp_vec, vbar_jk
      real :: deltavp, tau_cool1_par, dt1_cool
      real :: tausp1_par, tausp1_parj, tausp1_park, tausp_parj, tausp_park
      real :: tausp_parj3, tausp_park3, rhop_swarm_par
      integer :: j, k, l, ix0
      integer :: ispecies, jspecies
!
      if (ldiagnos .or. lentropy .and. lcollisional_heat) coll_heat=0.0
!
!  Add collisional cooling of the rms speed.
!
      if (lcollisional_cooling_taucool) then
        if (npar_imn(imn)/=0) then
          vvpm=0.0
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            vvpm(ix0-nghost,:) = vvpm(ix0-nghost,:) + fp(k,ivpx:ivpz)
          enddo
          do l=1,nx
            if (p%np(l)>1.0) vvpm(l,:)=vvpm(l,:)/p%np(l)
          enddo
!
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - &
                taucool1*(fp(k,ivpx:ivpz)-vvpm(ix0-nghost,:))
          enddo
!
          if (lfirst.and.ldt) dt1_max=max(dt1_max,taucool1/cdtp)
        endif
      endif
!
!  Add collisional cooling of the rms speed, with cooling time-scale based
!  on friction time and local rms speed of particles.
!
      if (lcollisional_cooling_rms) then
        if (npar_imn(imn)/=0) then
!  When multiple friction times are present, the average is used for the
!  number density in each superparticle.
          if (npar_species>1) then
            tausp1m=0.0
          else
            call get_frictiontime(f,fp,p,ineargrid,1,tausp1_par, &
                nochange_opt=.true.)
          endif
!  Need vpm=<|vvp-<vvp>|> to calculate the collisional time-scale.
          vvpm=0.0; vpm=0.0
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            vvpm(ix0-nghost,:) = vvpm(ix0-nghost,:) + fp(k,ivpx:ivpz)
            if (npar_species>1) then
              call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par, &
                  nochange_opt=.true.)
              tausp1m(ix0-nghost) = tausp1m(ix0-nghost) + tausp1_par
            endif
          enddo
          do l=1,nx
            if (p%np(l)>1.0) then
              vvpm(l,:)=vvpm(l,:)/p%np(l)
              if (npar_species>1) tausp1m=tausp1m/p%np
            endif
          enddo
!  vpm
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            vpm(ix0-nghost) = vpm(ix0-nghost) + &
                sqrt( (fp(k,ivpx)-vvpm(ix0-nghost,1))**2 + &
                      (fp(k,ivpy)-vvpm(ix0-nghost,2))**2 + &
                      (fp(k,ivpz)-vvpm(ix0-nghost,3))**2 )
          enddo
          do l=1,nx
            if (p%np(l)>1.0) then
              vpm(l)=vpm(l)/p%np(l)
            endif
          enddo
!  The collisional time-scale is 1/tau_coll=nd*vrms*sigma_coll.
!  Inserting Epstein friction time gives 1/tau_coll=3*rhod/rho*vprms/tauf.
          if (npar_species>1) then
            tau_coll1=(1.0-coeff_restitution)*p%epsp*vpm*tausp1m
          else
            tau_coll1=(1.0-coeff_restitution)*p%epsp*vpm*tausp1_par
          endif
!  Limit inverse time-step of collisional cooling if requested.
          if (tau_coll_min>0.0) then
            where (tau_coll1>tau_coll1_max) tau_coll1=tau_coll1_max
          endif
!
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - &
                tau_coll1(ix0-nghost)*(fp(k,ivpx:ivpz)-vvpm(ix0-nghost,:))
            if (lcollisional_heat .or. ldiagnos) then
              call get_rhopswarm(mp_swarm,fp,k,ineargrid(k,:),rhop_swarm_par)
              coll_heat(ix0-nghost) = coll_heat(ix0-nghost) + &
                   rhop_swarm_par*tau_coll1(ix0-nghost)*&
                   sum(fp(k,ivpx:ivpz)*(fp(k,ivpx:ivpz)-vvpm(ix0-nghost,:)))
            endif
          enddo
!
          if (lfirst.and.ldt) dt1_max=max(dt1_max,tau_coll1/cdtp)
        endif
      endif
!
!  More advanced collisional cooling model. Collisions are considered for
!  every possible two-body process in a grid cell.
!
      if (lcollisional_cooling_twobody) then
        do l=1,nx
! Collisions between particle k and all other particles in the grid cell.
          k=kshepherd(l)
          if (k>0) then
!  Limit inverse time-step of collisional cooling if requested.
            do while (k/=0)
              dt1_cool=0.0
              call get_frictiontime(f,fp,p,ineargrid,k,tausp1_park, &
                  nochange_opt=.true.)
              tausp_park=1/tausp1_park
              tausp_park3=tausp_park**3
              j=k
              do while (kneighbour(j)/=0)
!  Collide with the neighbours of k and their neighbours.
                j=kneighbour(j)
                call get_frictiontime(f,fp,p,ineargrid,j,tausp1_parj, &
                    nochange_opt=.true.)
                tausp_parj=1/tausp1_parj
                tausp_parj3=tausp_parj**3
!  Collision velocity.
                deltavp_vec=fp(k,ivpx:ivpz)-fp(j,ivpx:ivpz)
                deltavp=sqrt( deltavp_vec(1)**2 + deltavp_vec(2)**2 + &
                              deltavp_vec(3)**2 )
                vbar_jk= &
                    (tausp_parj3*fp(k,ivpx:ivpz)+tausp_park3*fp(j,ivpx:ivpz))/ &
                    (tausp_parj3+tausp_park3)
!  Cooling time-scale.
                call get_rhopswarm(mp_swarm,fp,k,ineargrid(k,:),rhop_swarm_par)
                tau_cool1_par= &
                     (1.0-coeff_restitution)* &
                     rhop_swarm_par*deltavp*(tausp_parj+tausp_park)**2/ &
                     (tausp_parj3+tausp_park3)
                dt1_cool=dt1_cool+tau_cool1_par
!                if (tau_coll_min>0.0) then
!                  if (tau_cool1_par>tau_coll1_max) tau_cool1_par=tau_coll1_max
!                endif
                dfp(j,ivpx:ivpz) = dfp(j,ivpx:ivpz) - &
                    tau_cool1_par*(fp(j,ivpx:ivpz)-vbar_jk)
                dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - &
                    tau_cool1_par*(fp(k,ivpx:ivpz)-vbar_jk)
              enddo
              if (lfirst.and.ldt) dt1_max=max(dt1_max(l),dt1_cool/cdtp)
!  Go through all possible k.
              k=kneighbour(k)
            enddo
          endif
        enddo
!
      endif
!
!  Treat collisions as a drag force that damps the rms speed at the same
!  time-scale.
!
      if (lcollisional_dragforce_cooling) then
        if (npar_imn(imn)/=0) then
          vvpm_species=0.0; vpm_species=0.0; np_species=0.0
!  Calculate mean velocity and number of particles for each species.
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            ispecies=npar_species*(ipar(k)-1)/npar+1
            vvpm_species(ix0-nghost,:,ispecies) = &
                vvpm_species(ix0-nghost,:,ispecies) + fp(k,ivpx:ivpz)
            np_species(ix0-nghost,ispecies)  = &
                np_species(ix0-nghost,ispecies) + 1.0
          enddo
          do l=1,nx
            do ispecies=1,npar_species
              if (np_species(l,ispecies)>1.0) then
                vvpm_species(l,:,ispecies)=vvpm_species(l,:,ispecies)/np_species(l,ispecies)
              endif
            enddo
          enddo
!  Calculate rms speed for each species.
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            ispecies=npar_species*(ipar(k)-1)/npar+1
            vpm_species(ix0-nghost,ispecies) = &
                vpm_species(ix0-nghost,ispecies) + sqrt( &
                (fp(k,ivpx)-vvpm_species(ix0-nghost,1,ispecies))**2 + &
                (fp(k,ivpy)-vvpm_species(ix0-nghost,2,ispecies))**2 + &
                (fp(k,ivpz)-vvpm_species(ix0-nghost,3,ispecies))**2 )
          enddo
          do l=1,nx
            do ispecies=1,npar_species
              if (np_species(l,ispecies)>1.0) then
                vpm_species(l,ispecies)=vpm_species(l,ispecies)/np_species(l,ispecies)
              endif
            enddo
          enddo
!
!  Collisional drag force time-scale between particles i and j with R_i < R_j.
!
!    tau_ji = tau_j^3/(tau_i+tau_j)^2/(deltav_ij/cs*rhoi/rhog)
!    tau_ij = tau_ji*rho_j/rho_i
!
          do ispecies=1,npar_species; do jspecies=ispecies,npar_species
            vcoll= &
                sqrt(vpm_species(:,ispecies)**2+vpm_species(:,ispecies)**2 + &
                  (vvpm_species(:,1,ispecies)-vvpm_species(:,1,jspecies))**2 + &
                  (vvpm_species(:,2,ispecies)-vvpm_species(:,2,jspecies))**2 + &
                  (vvpm_species(:,3,ispecies)-vvpm_species(:,3,jspecies))**2)
            call get_rhopswarm(mp_swarm,fp,k,m,n,rhop_swarm_mn)
            tau_coll1_species(:,jspecies,ispecies) = &
                 vcoll*np_species(:,ispecies)*rhop_swarm_mn*p%rho1 / ( &
                 tausp_species(jspecies)**3/ &
                 (tausp_species(ispecies)+tausp_species(jspecies))**2 )
            where (np_species(:,ispecies)/=0.0) &
              tau_coll1_species(:,ispecies,jspecies)= &
                   tau_coll1_species(:,jspecies,ispecies)*np_species(:,jspecies)/np_species(:,ispecies)
          enddo; enddo
!
          tau_coll1_tot=0.0
          do ispecies=1,npar_species; do jspecies=1,npar_species
            tau_coll1_tot(:,ispecies)=tau_coll1_tot(:,ispecies)+tau_coll1_species(:,ispecies,jspecies)
          enddo; enddo
!  Limit inverse time-step of collisional cooling if requested.
          if (tau_coll_min>0.0) then
            do ispecies=1,npar_species; do l=1,nx
              if (tau_coll1_tot(l,ispecies) > tau_coll1_max) then
                tau_coll1_species(l,ispecies,:)=tau_coll1_species(l,ispecies,:)* &
                    tau_coll1_max/tau_coll1_tot(l,ispecies)
              endif
            enddo; enddo
            tau_coll1_tot=0.0
            do ispecies=1,npar_species; do jspecies=1,npar_species
              tau_coll1_tot(:,ispecies)=tau_coll1_tot(:,ispecies)+tau_coll1_species(:,ispecies,jspecies)
            enddo; enddo
          endif
          if (lfirst.and.ldt) then
            do ispecies=1,npar_species
              dt1_max=max(dt1_max,tau_coll1_tot(:,ispecies)/cdtp)
            enddo
          endif
!  Add to equation of motion.
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            ispecies=npar_species*(ipar(k)-1)/npar+1
            do jspecies=1,npar_species
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - &
                  tau_coll1_species(ix0-nghost,ispecies,jspecies)* &
                  (fp(k,ivpx:ivpz)-vvpm_species(ix0-nghost,:,jspecies))
              if (lcollisional_heat .or. ldiagnos) then
                call get_rhopswarm(mp_swarm,fp,k,ineargrid(k,:),rhop_swarm_par)
                coll_heat(ix0-nghost) = coll_heat(ix0-nghost) + &
                     rhop_swarm_par* &
                     tau_coll1_species(ix0-nghost,ispecies,jspecies)* &
                     sum(fp(k,ivpx:ivpz)*(fp(k,ivpx:ivpz) - &
                     vvpm_species(ix0-nghost,:,jspecies)))
              endif
            enddo
          enddo
        endif
      endif
!
!  Heating of the gas due to dissipative collisions.
!
      if (lentropy .and. lcollisional_heat) &
          df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%rho1*p%TT1*coll_heat
!
!  Diagnostics.
      if (ldiagnos) then
        if (idiag_decollp/=0) &
          call sum_mn_name(coll_heat,idiag_decollp)
      endif
!
    endsubroutine collisional_cooling
!***********************************************************************
    subroutine compensate_friction_increase(f,fp,dfp,p,ineargrid)
!
!  Compensate for increased friction time in regions of high solids-to-gas
!  ratio by applying missing friction force to particles only.
!
!  26-feb-07/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, dimension (nx,3) :: vvpm
      real :: tausp1_par, tausp1_par_mod, tausp1_par_org
      integer :: k, l, ix0
!
      if (npar_imn(imn)/=0) then
!  Calculate mean particle velocity.
        vvpm=0.0
        do k=k1_imn(imn),k2_imn(imn)
          ix0=ineargrid(k,1)
          vvpm(ix0-nghost,:) = vvpm(ix0-nghost,:) + fp(k,ivpx:ivpz)
        enddo
        do l=1,nx
          if (p%np(l)>1.0) vvpm(l,:)=vvpm(l,:)/p%np(l)
        enddo
!
!  Compare actual and modified friction time and apply the difference as
!  friction force relative to mean particle velocity.
!
        do k=k1_imn(imn),k2_imn(imn)
          call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par_mod, &
              nochange_opt=.false.)
          call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par_org, &
              nochange_opt=.true.)
          tausp1_par=tausp1_par_org-tausp1_par_mod
          ix0=ineargrid(k,1)
          dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - &
              tausp1_par*(fp(k,ivpx:ivpz)-vvpm(ix0-nghost,:))
        enddo
      endif
!
    endsubroutine compensate_friction_increase
!***********************************************************************
    subroutine calc_gas_velocity_shell_call(k1, uup, fp)
!
      use Special, only: special_calc_particles
!
      real,dimension(mpar_loc,mpvar) :: fp
      real, dimension(3) :: uup
      integer :: k1
!
      vel_call=.true.
      uup_shared=fp(k1,ixp:izp)
      call special_calc_particles(fp)
      uup=uup_shared
!
    endsubroutine calc_gas_velocity_shell_call
!***********************************************************************
    subroutine calculate_rms_speed(fp,ineargrid,p)
!
      use Diagnostics
!
!  Calculate the rms speed dvpm=sqrt(<(vvp-<vvp>)^2>) of the
!  particle for diagnostic purposes
!
!  08-04-08/wlad: coded
!
      real,dimension(mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      real,dimension(nx,3) :: vvpm,dvp2m
      integer :: inx0,k,l
      type (pencil_case) :: p
      logical :: lnbody
!
!  Initialize the variables
!
      vvpm=0.0; dvp2m=0.0
!
!  Calculate the average velocity at each cell
!  if there are particles in the pencil only
!
      if (npar_imn(imn)/=0) then
!
        do k=k1_imn(imn),k2_imn(imn)
          lnbody=any(ipar(k)==ipar_nbody)
          if (.not.lnbody) then
            inx0=ineargrid(k,1)-nghost
            vvpm(inx0,:) = vvpm(inx0,:) + fp(k,ivpx:ivpz)
          endif
        enddo
        do l=1,nx
          if (p%np(l)>1.0) vvpm(l,:)=vvpm(l,:)/p%np(l)
        enddo
!
!  Get the residual in quadrature, dvp2m. Need vvpm calculated above.
!
        do k=k1_imn(imn),k2_imn(imn)
          lnbody=any(ipar(k)==ipar_nbody)
          if (.not.lnbody) then
            inx0=ineargrid(k,1)-nghost
            dvp2m(inx0,1)=dvp2m(inx0,1)+(fp(k,ivpx)-vvpm(inx0,1))**2
            dvp2m(inx0,2)=dvp2m(inx0,2)+(fp(k,ivpy)-vvpm(inx0,2))**2
            dvp2m(inx0,3)=dvp2m(inx0,3)+(fp(k,ivpz)-vvpm(inx0,3))**2
          endif
        enddo
        do l=1,nx
          if (p%np(l)>1.0) dvp2m(l,:)=dvp2m(l,:)/p%np(l)
        enddo
!
      endif
!
!  Output the diagnostics
!
      if (idiag_dvpx2m/=0) call sum_mn_name(dvp2m(:,1),idiag_dvpx2m)
      if (idiag_dvpy2m/=0) call sum_mn_name(dvp2m(:,2),idiag_dvpy2m)
      if (idiag_dvpz2m/=0) call sum_mn_name(dvp2m(:,3),idiag_dvpz2m)
      if (idiag_dvpm/=0)   call sum_mn_name(dvp2m(:,1)+dvp2m(:,2)+dvp2m(:,3),&
                                            idiag_dvpm,lsqrt=.true.)
      if (idiag_dvpmax/=0) call max_mn_name(dvp2m(:,1)+dvp2m(:,2)+dvp2m(:,3),&
                                            idiag_dvpmax,lsqrt=.true.)
!
    endsubroutine calculate_rms_speed
!***********************************************************************
    subroutine calc_pencil_rep(fp,rep)
!
!  Calculate particle Reynolds numbers
!
!  16-jul-08/kapelrud: coded
!
      use Viscosity, only: getnu
!
      real,dimension(mpar_loc,mpvar) :: fp
      real,dimension(k1_imn(imn):k2_imn(imn)) :: rep,nu
      character (len=labellen) :: ivis=''
      intent(in) :: fp
      intent(inout) :: rep
!
      real :: nu_
      integer :: k
!
      call getnu(nu_input=nu_,IVIS=ivis)
      if (ivis=='nu-const') then
        nu=nu_
      elseif (ivis=='rho-nu-const') then
        nu=nu_/interp_rho(k1_imn(imn):k2_imn(imn))
      elseif (ivis=='sqrtrho-nu-const') then
        nu=nu_/sqrt(interp_rho(k1_imn(imn):k2_imn(imn)))
      elseif (ivis=='nu-therm') then
        nu=nu_*sqrt(interp_TT(k1_imn(imn):k2_imn(imn)))
      elseif (ivis=='mu-therm') then
        nu=nu_*sqrt(interp_TT(k1_imn(imn):k2_imn(imn)))&
            /interp_rho(k1_imn(imn):k2_imn(imn))
      else
        call fatal_error('calc_pencil_rep','No such ivis!')
      endif
!
      if (.not.lparticles_radius) then
        print*,'calc_pencil_rep: particle_radius module needs to be '// &
          'enabled to calculate the particles Reynolds numbers.'
        call fatal_error('calc_pencil_rep','')
      elseif (maxval(nu)==0.0) then
        print*,'calc_pencil_rep: nu (kinematic visc.) must be non-zero!'
        call fatal_error('calc_pencil_rep','')
      endif
!
      do k=k1_imn(imn),k2_imn(imn)
        rep(k)=2.0*fp(k,iap)*sqrt(sum((interp_uu(k,:)-fp(k,ivpx:ivpz))**2))/nu(k)
      enddo
!
    endsubroutine calc_pencil_rep
!***********************************************************************
    subroutine calc_stokes_cunningham(fp,stocunn)
!
!  Calculate the Stokes-Cunningham factor
!
!  12-aug-08/kapelrud: coded
!
      use Particles_radius
!
      real,dimension(mpar_loc,mpvar) :: fp
      real,dimension(k1_imn(imn):k2_imn(imn)) :: stocunn
!
      real :: dia
      integer :: k
!
      do k=k1_imn(imn),k2_imn(imn)
!
!  Particle diameter
!
        dia=2.0*fp(k,iap)
!
        stocunn(k)=1.+2.*mean_free_path_gas/dia* &
          (1.257+0.4*exp(-0.55*dia/mean_free_path_gas))
!
      enddo
!
    endsubroutine calc_stokes_cunningham
!***********************************************************************
    subroutine calc_draglaw_steadystate(fp,k,rep,stocunn,tausp1_par)
!
!   Calculate relaxation time for particles under steady state drag.
!
!   15-jul-08/kapelrud: coded
!
      use Viscosity, only: getnu
      use Particles_radius
!
      real, dimension(mpar_loc,mpvar) :: fp
      integer :: k
      real :: rep, stocunn, tausp1_par
      character (len=labellen) :: ivis=''
!
      intent(in) :: fp,k,rep,stocunn
      intent(out) :: tausp1_par
!
      real :: cdrag,dia,nu,nu_
!
!  Find the kinematic viscosity
!
      call getnu(nu_input=nu_,ivis=ivis)
      if (ivis=='nu-const') then
        nu=nu_
      elseif (ivis=='rho-nu-const') then
        nu=nu_/interp_rho(k)
      elseif (ivis=='sqrtrho-nu-const') then
        nu=nu_/sqrt(interp_rho(k))
      elseif (ivis=='nu-therm') then
        nu=nu_*sqrt(interp_TT(k))
      elseif (ivis=='mu-therm') then
        nu=nu_*sqrt(interp_TT(k))&
            /interp_rho(k)
      else
        call fatal_error('calc_pencil_rep','No such ivis!')
      endif
!
!  Particle diameter
!
      if (.not.lparticles_radius) then
        print*,'calc_draglaw_steadystate: need particles_radius module to '// &
            'calculate the relaxation time!'
        call fatal_error('calc_draglaw_steadystate','')
      endif
!
      dia=2.0*fp(k,iap)
!
!  Calculate drag coefficent pre-factor:
!
      if (rep<1) then
        cdrag=1.0
      elseif (rep>1000) then
        cdrag=0.44*rep/24.0
      else
        cdrag=(1.+0.15*rep**0.687)
      endif
!
!  Relaxation time:
!
      tausp1_par=18.0*cdrag*nu/((rhopmat/interp_rho(k))*stocunn*dia**2)
!
    endsubroutine calc_draglaw_steadystate
!***********************************************************************
    subroutine calc_brownian_force(fp,k,ineark,stocunn,force)
!
!  Calculate the Brownian force contribution due to the random thermal motions
!  of the gas molecules.
!
!  28-jul-08/kapelrud: coded
!
      use General, only: normal_deviate
      use Viscosity, only: getnu
!
      real, dimension(mpar_loc,mpvar) :: fp
      real, dimension(3), intent(out) :: force
      character (len=labellen) :: ivis=''
      integer, dimension(3) :: ineark
      integer :: k
      real :: stocunn, rhop_swarm_par
!
      intent(in) :: fp,k
!
      real :: Szero,dia,TT,nu,nu_
!
!  Find kinematic viscosity
!
      call getnu(nu_input=nu_,ivis=ivis)
      if (ivis=='nu-const') then
        nu=nu_
      elseif (ivis=='rho-nu-const') then
        nu=nu_/interp_rho(k)
      elseif (ivis=='sqrtrho-nu-const') then
        nu=nu_/sqrt(interp_rho(k))
      elseif (ivis=='nu-therm') then
        nu=nu_*sqrt(interp_TT(k))
      elseif (ivis=='mu-therm') then
        nu=nu_*sqrt(interp_TT(k))&
            /interp_rho(k)
      else
        call fatal_error('calc_pencil_rep','No such ivis!')
      endif
!
!  Particle diameter:
!
      dia=2.0*fp(k,iap)
!
!  Get zero mean, unit variance Gaussian random numbers:
!
      call normal_deviate(force(1))
      call normal_deviate(force(2))
      call normal_deviate(force(3))
!
      if (interp%lTT) then
        TT=interp_TT(k)
      else
        TT=brownian_T0
      endif
!
      call get_rhopswarm(mp_swarm,fp,k,ineark,rhop_swarm_par)
!
      Szero=216*nu*k_B*TT*pi_1/ &
           (dia**5*stocunn*rhop_swarm_par**2/interp_rho(k))
!
      if (dt==0.0) then
        force=0.0
      else
        force=force*sqrt(Szero/dt)
      endif
!
    endsubroutine calc_brownian_force
!***********************************************************************
 subroutine calc_thermophoretic_force(fp,k,ineark,force)
!
!  Calculate the Thermophoretic force due to a temperature gradient in the gas.
!
!   Henrik Lutro, testing
!
      use Viscosity, only: getnu
!
      real, dimension(mpar_loc,mpvar) :: fp
      real, dimension(3), intent(out) :: force
      real, dimension(3) :: temp_grad
      real TT,mu,nu_,Kn,phi_therm,mass_p
      real Ktc,Ce,Cm,Cint
      character (len=labellen) :: ivis=''
      integer, dimension(3) :: ineark
      integer :: k
      intent(in) :: fp,k
!
      call keep_compiler_quiet(ineark)
!
      Ktc=1.10
      Cm=1.13
      Ce=2.17
      if (interp%lTT) then
        TT=interp_TT(k)
      else
        TT=thermophoretic_T0
      endif
      !Find dynamic viscosity
      call getnu(nu_input=nu_,ivis=ivis)
      if (ivis=='nu-const') then
        mu=nu_*interp_rho(k)
      elseif (ivis=='rho-nu-const') then
        mu=nu_
      elseif (ivis=='sqrtrho-nu-const') then
        mu=nu_*sqrt(interp_rho(k))
      elseif (ivis=='nu-therm') then
        mu=nu_*interp_rho(k)*sqrt(TT)
      elseif (ivis=='mu-therm') then
        mu=nu_*sqrt(TT)
      else
        call fatal_error('calc_pencil_rep','No such ivis!')
      endif
      Cint=0.5
      if (interp%lgradTT) then
        temp_grad=interp_gradTT(k,:)
      else
        temp_grad=temp_grad0
     endif
!
      Kn=mean_free_path_gas/fp(k,iap)
!
      mass_p=(4.0*pi/3.0)*rhopmat*fp(k,iap)**3
      if (thermophoretic_eq=='near_continuum') then
        phi_therm=-9*pi/cond_ratio
      elseif (thermophoretic_eq=='transition') then
        phi_therm=-12.0*pi*(Ktc*(1.0+cond_ratio*Ce*Kn)+3.0*Cm*Kn*(1.0-cond_ratio+cond_ratio*Ce*Kn))&
            /((1.0+3.0*Kn*exp(-Cint/Kn))*(1.0+3.0*Cm*Kn)*(2.0+cond_ratio+2.0*cond_ratio*Ce*Kn))
      elseif (thermophoretic_eq=='free_molecule') then
        phi_therm=0.0
      else
        call fatal_error('calc_pencil_rep','No thermoporetic range chosen')
      endif
!
      force=(fp(k,iap)*temp_grad*mu**2*phi_therm)/(TT*interp_rho(k)*mass_p)
!
    endsubroutine calc_thermophoretic_force
!***********************************************************************
    subroutine read_particles_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_init_pars
!***********************************************************************
    subroutine write_particles_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_init_pars)
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_run_pars
!***********************************************************************
    subroutine write_particles_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_run_pars)
!
    endsubroutine write_particles_run_pars
!***********************************************************************
    subroutine powersnap_particles(f)
!
!  Calculate power spectra of dust particle variables.
!
!  01-jan-06/anders: coded
!
      use Power_spectrum, only: power_1d
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lpar_spec) call power_1d(f,'p',0,irhop)
!
    endsubroutine powersnap_particles
!***********************************************************************
    subroutine rprint_particles(lreset,lwrite)
!
!  Read and register print parameters relevant for particles.
!
!  29-dec-04/anders: coded
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname,inamez,inamey,inamex,inamexy,inamexz,inamer,inamerz
      logical :: lwr
!
!  Write information to index.pro.
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lwr) then
        write(3,*) 'ixp=', ixp
        write(3,*) 'iyp=', iyp
        write(3,*) 'izp=', izp
        write(3,*) 'ivpx=', ivpx
        write(3,*) 'ivpy=', ivpy
        write(3,*) 'ivpz=', ivpz
        write(3,*) 'inp=', inp
        write(3,*) 'irhop=', irhop
        write(3,*) 'iupx=', iupx
        write(3,*) 'iupy=', iupy
        write(3,*) 'iupz=', iupz
        write(3,*) 'ifgx=', ifgx
        write(3,*) 'ifgy=', ifgy
        write(3,*) 'ifgz=', ifgz
      endif
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_xpm=0; idiag_ypm=0; idiag_zpm=0
        idiag_xp2m=0; idiag_yp2m=0; idiag_zp2m=0; idiag_rpm=0; idiag_rp2m=0
        idiag_vpxm=0; idiag_vpym=0; idiag_vpzm=0
        idiag_vpxvpym=0; idiag_vpxvpzm=0; idiag_vpyvpzm=0
        idiag_vpx2m=0; idiag_vpy2m=0; idiag_vpz2m=0; idiag_ekinp=0
        idiag_vpxmax=0; idiag_vpymax=0; idiag_vpzmax=0; idiag_vpmax=0
        idiag_rhopvpxm=0; idiag_rhopvpym=0; idiag_rhopvpzm=0; idiag_rhopvpysm=0
        idiag_rhopvpxt=0; idiag_rhopvpyt=0; idiag_rhopvpzt=0
        idiag_lpxm=0; idiag_lpym=0; idiag_lpzm=0
        idiag_lpx2m=0; idiag_lpy2m=0; idiag_lpz2m=0
        idiag_npm=0; idiag_np2m=0; idiag_npmax=0; idiag_npmin=0
        idiag_dtdragp=0; idiag_dedragp=0
        idiag_rhopm=0; idiag_rhoprms=0; idiag_rhop2m=0; idiag_rhopmax=0
        idiag_rhopmin=0; idiag_decollp=0; idiag_rhopmphi=0
        idiag_epspmin=0; idiag_epspmax=0
        idiag_nparmin=0; idiag_nparmax=0; idiag_nmigmax=0; idiag_mpt=0
        idiag_npmx=0; idiag_npmy=0; idiag_npmz=0; idiag_epotpm=0
        idiag_rhopmx=0; idiag_rhopmy=0; idiag_rhopmz=0
        idiag_epspmx=0; idiag_epspmy=0; idiag_epspmz=0
        idiag_rhopmxy=0; idiag_rhopmxz=0; idiag_rhopmr=0
        idiag_dvpx2m=0; idiag_dvpy2m=0; idiag_dvpz2m=0
        idiag_dvpmax=0; idiag_dvpm=0; idiag_nparpmax=0
        idiag_eccpxm=0; idiag_eccpym=0; idiag_eccpzm=0
        idiag_eccpx2m=0; idiag_eccpy2m=0; idiag_eccpz2m=0
        idiag_npargone=0; idiag_vpyfull2m=0; idiag_deshearbcsm=0
        idiag_npmxy=0; idiag_vprms=0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip<14) print*,'rprint_particles: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'nparmin',idiag_nparmin)
        call parse_name(iname,cname(iname),cform(iname),'nparmax',idiag_nparmax)
        call parse_name(iname,cname(iname),cform(iname),'nparpmax',idiag_nparpmax)
        call parse_name(iname,cname(iname),cform(iname),'xpm',idiag_xpm)
        call parse_name(iname,cname(iname),cform(iname),'ypm',idiag_ypm)
        call parse_name(iname,cname(iname),cform(iname),'zpm',idiag_zpm)
        call parse_name(iname,cname(iname),cform(iname),'xp2m',idiag_xp2m)
        call parse_name(iname,cname(iname),cform(iname),'yp2m',idiag_yp2m)
        call parse_name(iname,cname(iname),cform(iname),'zp2m',idiag_zp2m)
        call parse_name(iname,cname(iname),cform(iname),'rpm',idiag_rpm)
        call parse_name(iname,cname(iname),cform(iname),'rp2m',idiag_rp2m)
        call parse_name(iname,cname(iname),cform(iname),'vpxm',idiag_vpxm)
        call parse_name(iname,cname(iname),cform(iname),'vpym',idiag_vpym)
        call parse_name(iname,cname(iname),cform(iname),'vpzm',idiag_vpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpxvpym',idiag_vpxvpym)
        call parse_name(iname,cname(iname),cform(iname),'vpxvpzm',idiag_vpxvpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpyvpzm',idiag_vpyvpzm)
        call parse_name(iname,cname(iname),cform(iname),'vpx2m',idiag_vpx2m)
        call parse_name(iname,cname(iname),cform(iname),'vpy2m',idiag_vpy2m)
        call parse_name(iname,cname(iname),cform(iname),'vpz2m',idiag_vpz2m)
        call parse_name(iname,cname(iname),cform(iname),'ekinp',idiag_ekinp)
        call parse_name(iname,cname(iname),cform(iname),'vpxmax',idiag_vpxmax)
        call parse_name(iname,cname(iname),cform(iname),'vpymax',idiag_vpymax)
        call parse_name(iname,cname(iname),cform(iname),'vpzmax',idiag_vpzmax)
        call parse_name(iname,cname(iname),cform(iname),'vpmax',idiag_vpmax)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpxm', &
            idiag_rhopvpxm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpym', &
            idiag_rhopvpym)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpzm', &
            idiag_rhopvpzm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpysm', &
            idiag_rhopvpysm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpxt', &
            idiag_rhopvpxt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpyt', &
            idiag_rhopvpyt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpzt', &
            idiag_rhopvpzt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpysm', &
            idiag_rhopvpysm)
        call parse_name(iname,cname(iname),cform(iname),'lpxm',idiag_lpxm)
        call parse_name(iname,cname(iname),cform(iname),'lpym',idiag_lpym)
        call parse_name(iname,cname(iname),cform(iname),'lpzm',idiag_lpzm)
        call parse_name(iname,cname(iname),cform(iname),'lpx2m',idiag_lpx2m)
        call parse_name(iname,cname(iname),cform(iname),'lpy2m',idiag_lpy2m)
        call parse_name(iname,cname(iname),cform(iname),'lpz2m',idiag_lpz2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpxm',idiag_eccpxm)
        call parse_name(iname,cname(iname),cform(iname),'eccpym',idiag_eccpym)
        call parse_name(iname,cname(iname),cform(iname),'eccpzm',idiag_eccpzm)
        call parse_name(iname,cname(iname),cform(iname),'eccpx2m',idiag_eccpx2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpy2m',idiag_eccpy2m)
        call parse_name(iname,cname(iname),cform(iname),'eccpz2m',idiag_eccpz2m)
        call parse_name(iname,cname(iname),cform(iname),'dtdragp',idiag_dtdragp)
        call parse_name(iname,cname(iname),cform(iname),'npm',idiag_npm)
        call parse_name(iname,cname(iname),cform(iname),'np2m',idiag_np2m)
        call parse_name(iname,cname(iname),cform(iname),'npmax',idiag_npmax)
        call parse_name(iname,cname(iname),cform(iname),'npmin',idiag_npmin)
        call parse_name(iname,cname(iname),cform(iname),'rhopm',idiag_rhopm)
        call parse_name(iname,cname(iname),cform(iname),'rhoprms',idiag_rhoprms)
        call parse_name(iname,cname(iname),cform(iname),'rhop2m',idiag_rhop2m)
        call parse_name(iname,cname(iname),cform(iname),'rhopmin',idiag_rhopmin)
        call parse_name(iname,cname(iname),cform(iname),'rhopmax',idiag_rhopmax)
        call parse_name(iname,cname(iname),cform(iname),'epspmin',idiag_epspmin)
        call parse_name(iname,cname(iname),cform(iname),'epspmax',idiag_epspmax)
        call parse_name(iname,cname(iname),cform(iname),'rhopmphi',idiag_rhopmphi)
        call parse_name(iname,cname(iname),cform(iname),'nmigmax',idiag_nmigmax)
        call parse_name(iname,cname(iname),cform(iname),'mpt',idiag_mpt)
        call parse_name(iname,cname(iname),cform(iname),'dvpx2m',idiag_dvpx2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpy2m',idiag_dvpy2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpz2m',idiag_dvpz2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpm',idiag_dvpm)
        call parse_name(iname,cname(iname),cform(iname),'dvpmax',idiag_dvpmax)
        call parse_name(iname,cname(iname),cform(iname), &
            'dedragp',idiag_dedragp)
        call parse_name(iname,cname(iname),cform(iname), &
            'decollp',idiag_decollp)
        call parse_name(iname,cname(iname),cform(iname), &
            'epotpm',idiag_epotpm)
        call parse_name(iname,cname(iname),cform(iname), &
            'npargone',idiag_npargone)
        call parse_name(iname,cname(iname),cform(iname), &
            'vpyfull2m',idiag_vpyfull2m)
        call parse_name(iname,cname(iname),cform(iname),'vprms',idiag_vprms)
        call parse_name(iname,cname(iname),cform(iname), &
            'deshearbcsm',idiag_deshearbcsm)
      enddo
!
!  Check for those quantities for which we want x-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'npmx',idiag_npmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhopmx',idiag_rhopmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'epspmx',idiag_epspmx)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'npmy',idiag_npmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhopmy',idiag_npmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'epspmy',idiag_epspmy)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'npmz',idiag_npmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhopmz',idiag_rhopmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'epspmz',idiag_epspmz)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'npmxy',idiag_npmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhopmxy',idiag_rhopmxy)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'rhopmxz',idiag_rhopmxz)
      enddo
!
!  Check for those quantities for which we want phiz-averages.
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'rhopmr',idiag_rhopmr)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do inamerz=1,nnamerz
        call parse_name(inamerz,cnamerz(inamerz),cformrz(inamerz),'rhopmphi',idiag_rhopmphi)
      enddo
!
    endsubroutine rprint_particles
!***********************************************************************
endmodule Particles
