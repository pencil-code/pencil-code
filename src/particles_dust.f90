! $Id: particles_dust.f90,v 1.1 2018/08/24 15:48:10 wlyra Exp $
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
! CPARAM character (len=20), parameter :: particles_module="dust"
!
! PENCILS PROVIDED np; rhop; vol; peh
! PENCILS PROVIDED uup(3)
! PENCILS PROVIDED np_rad(ndustrad); npvz(ndustrad); npvz2(ndustrad);
! PENCILS PROVIDED npuz(ndustrad); sherwood
! PENCILS PROVIDED epsp; grhop(3)
! PENCILS PROVIDED tauascalar
! PENCILS PROVIDED condensationRate
! PENCILS PROVIDED waterMixingRatio, part_heatcap
!
!***************************************************************
module Particles
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet, indgen
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_radius
  use Particles_potential
!
  implicit none
!
  include 'particles.h'
  include 'particles_common.h'
!
  complex, dimension(7) :: coeff=(0.0,0.0)
  real, target, dimension(npar_species) :: tausp_species=0.0
  real, target, dimension(npar_species) :: tausp1_species=0.0
  real, target, dimension(npar_species) :: rpbeta_species=0.0
  real, dimension(3) :: temp_grad0=(/0.0,0.0,0.0/)
  real, dimension(3) :: pos_sphere=(/0.0,0.0,0.0/)
  real, dimension(3) :: pos_ellipsoid=(/0.0,0.0,0.0/)
  real, target :: tausp = 0.0
  real, pointer :: true_density_cond_spec
  real :: xp0=0.0, yp0=0.0, zp0=0.0, vpx0=0.0, vpy0=0.0, vpz0=0.0
  real :: xp1=0.0, yp1=0.0, zp1=0.0, vpx1=0.0, vpy1=0.0, vpz1=0.0
  real :: xp2=0.0, yp2=0.0, zp2=0.0, vpx2=0.0, vpy2=0.0, vpz2=0.0
  real :: xp3=0.0, yp3=0.0, zp3=0.0, vpx3=0.0, vpy3=0.0, vpz3=0.0
  real :: Lx0=0.0, Ly0=0.0, Lz0=0.0
  real :: delta_vp0=1.0, tausp1=0.0, rpbeta=0.0
  real :: nu_epicycle=0.0, nu_epicycle2=0.0
  real :: beta_dPdr_dust=0.0, beta_dPdr_dust_scaled=0.0
  real :: tausg_min=0.0, tausg1_max=0.0, epsp_friction_increase=100.0
  real :: cdtp=0.2, cdtpgrav=0.1, cdtp_drag=0.2
  real :: gravx=0.0, gravz=0.0, gravr=1.0, kx_gg=1.0, kz_gg=1.0
  real :: gravsmooth=0.0, gravsmooth2=0.0, Ri0=0.25, eps1=0.5
  real :: kx_xxp=0.0, ky_xxp=0.0, kz_xxp=0.0, amplxxp=0.0
  real :: kx_vvp=0.0, ky_vvp=0.0, kz_vvp=0.0, amplvvp=0.0
  real :: kx_vpx=0.0, kx_vpy=0.0, kx_vpz=0.0
  real :: ky_vpx=0.0, ky_vpy=0.0, ky_vpz=0.0
  real :: kz_vpx=0.0, kz_vpy=0.0, kz_vpz=0.0
  real :: phase_vpx=0.0, phase_vpy=0.0, phase_vpz=0.0
  real :: gaussian_x_frac=0.0
  real :: tstart_dragforce_par=0.0
  real :: tstart_grav_par=0.0, tstart_grav_x_par=0.0
  real :: tstart_grav_z_par=0.0, tstart_grav_r_par=0.0
  real :: tstart_liftforce_par=0.0
  real :: tstart_brownian_par=0.0
  real :: tstart_collisional_cooling=0.0
  real :: tstart_sink_par=0.0
  real :: tau_coll_min=0.0, tau_coll1_max=0.0
  real :: coeff_restitution=0.5, mean_free_path_gas=0.0
  real :: rad_sphere=0.0
  real :: a_ellipsoid=0.0, b_ellipsoid=0.0, c_ellipsoid=0.0
  real :: a_ell2=0.0, b_ell2=0.0, c_ell2=0.0
  real :: taucool=0.0, taucool1=0.0, brownian_T0=0.0, thermophoretic_T0=0.0
  real :: xsinkpoint=0.0, ysinkpoint=0.0, zsinkpoint=0.0, rsinkpoint=0.0
  real :: particles_insert_rate=0.
  real :: avg_n_insert, remaining_particles=0.0
  real :: max_particle_insert_time=huge1
  real :: Deltauy_gas_friction=0.0
  real :: cond_ratio=0.0
  real :: pscalar_sink_rate=0.0
  real :: frac_init_particles=1.0, particles_insert_ramp_time=0.0
  real :: tstart_insert_particles=0.0
  real :: birthring_r=1.0, birthring_width=0.1
  real :: tstart_rpbeta=0.0, birthring_lifetime=huge1
  real :: rdiffconst_dragf=0.07, rdiffconst_pass=0.07
  real :: r0gaussz=1.0, qgaussz=0.0
  real :: vapor_mixing_ratio_qvs=0., rhoa=1.0, redfrac=0.9
  real, pointer :: g1, rp1, rp1_smooth, t_ramp_mass, t_start_secondary
  integer :: l_hole=0, m_hole=0, n_hole=0
  integer :: iffg=0, ifgx=0, ifgy=0, ifgz=0, ibrtime=0
  integer :: istep_dragf=3, istep_pass=3
  integer :: it_insert_nuclei=1
  logical, target :: ldragforce_gas_par=.false.
  logical :: ldragforce_dust_par=.false.
  logical :: ldragforce_stiff=.false., ldragforce_radialonly=.false.
  logical :: ldragforce_heat=.false., lcollisional_heat=.false.
  logical :: lpar_spec=.false., lcompensate_friction_increase=.false.
  logical :: lcollisional_cooling_taucool=.false.
  logical :: lcollisional_cooling_rms=.false.
  logical :: lcollisional_cooling_twobody=.false.
  logical :: lcollisional_dragforce_cooling=.false.
  logical :: ltau_coll_min_courant=.true.
  logical :: ldragforce_equi_global_eps=.false.
  logical :: ldraglaw_epstein=.true., ldraglaw_epstein_stokes_linear=.false.
  logical :: ldraglaw_simple=.false.
  logical :: ldraglaw_steadystate=.false., ldraglaw_variable=.false.
  logical :: ldraglaw_purestokes=.false.
  logical :: ldraglaw_stokesschiller=.false.
  logical :: ldraglaw_epstein_transonic=.false.
  logical :: ldraglaw_eps_stk_transonic=.false.
  logical :: ldraglaw_variable_density=.false.
  logical :: lcoldstart_amplitude_correction=.false.
  logical :: luse_tau_ap=.true.
  logical :: lbrownian_forces=.false.
  logical :: lbrownian_forces_Li_Ahmadi=.false.
  logical :: lthermophoretic_forces=.false.
  logical :: lenforce_policy=.false., lnostore_uu=.true.
  logical :: ldt_grav_par=.true., ldt_adv_par=.true.
  logical :: lsinkpoint=.false., lglobalrandom=.false.
  logical :: lcoriolis_force_par=.true., lcentrifugal_force_par=.false.
  logical :: lshear_accel_par = .true.
  logical :: lcalc_uup=.false.
  logical :: lparticle_gravity=.true.
  logical :: lcylindrical_gravity_par=.false.
  logical :: lpscalar_sink=.false.
  logical :: lsherwood_const=.false. ! RUN_DOC: Constant quiescent (2) Sherwood number
  logical :: lbubble=.false.
  logical :: linsert_as_many_as_possible=.false.
  logical :: lwithhold_init_particles=.false.
  logical :: lgaussian_birthring=.false.
  logical :: lvector_gravity=.false.
  logical :: lpeh_radius=.false.
  logical :: lbirthring_depletion=.false.
  logical :: lstocunn1 = .false.
  logical :: lprecalc_cell_volumes=.false.
  logical :: ldiffuse_passive = .false., ldiff_pass=.false.
  logical :: ldiffuse_dragf= .false., ldiff_dragf=.false.
  logical :: lsimple_volume=.false.
  logical :: lnpmin_exclude_zero = .false.
  logical :: ltauascalar = .false., lfollow_gas=.false.
  logical :: lset_df_insert_nucleii=.false.
  logical, pointer :: lramp_mass, lsecondary_wait
!
  character(len=labellen) :: interp_pol_uu ='ngp'
  character(len=labellen) :: interp_pol_oo ='ngp'
  character(len=labellen) :: interp_pol_TT ='ngp'
  character(len=labellen) :: interp_pol_rho='ngp'
  character(len=labellen) :: interp_pol_pp ='ngp'
  character(len=labellen) :: interp_pol_species='ngp'
  character(len=labellen) :: interp_pol_gradTT='ngp'
  character(len=labellen) :: interp_pol_nu='ngp'
!
  character(len=labellen), dimension(ninit) :: initxxp='nothing'
  character(len=labellen), dimension(ninit) :: initvvp='nothing'
  character(len=labellen) :: gravx_profile='', gravz_profile=''
  character(len=labellen) :: gravr_profile=''
  character(len=labellen) :: thermophoretic_eq= 'nothing'
!
  integer :: init_repeat=0       !repeat particle initialization for distance statistics
!
  integer :: icondensationRate=0
  integer :: iwaterMixingRatio=0
  logical :: lcondensation_rate=.false., lnp_ap_as_aux=.false.
!
!  Interactions with special/shell
!
  integer :: k_shell=-1            !k associated with minshell (special/shell.f90)
  logical :: l_shell=.false.       !using special/shell.f90 for gas velocities
!
  real, dimension(3) :: uup_shared=0
  real :: turnover_shared=0, nu_draglaw=0.
  logical :: vel_call=.false., turnover_call=.false.
  logical :: lreassign_strat_rhom=.true., lnu_draglaw=.false.
!
  logical :: lcdtp_shear = .true.
!
  logical :: lcompensate_sedimentation=.false.
  real :: compensate_sedimentation=1.
!
!  real :: A3=0., A2=0., G_condensation=0.
  real :: A3=0., A2=0.
  real, target :: G_condensation=0.0
  real :: nucleation_threshold
  logical :: ascalar_ngp=.false., ascalar_cic=.false.
!
  namelist /particles_init_pars/ &
      initxxp, initvvp, xp0, yp0, zp0, vpx0, vpy0, vpz0, delta_vp0, &
      ldragforce_gas_par, ldragforce_dust_par, bcpx, bcpy, bcpz, tausp, &
      beta_dPdr_dust, np_swarm, mp_swarm, mpmat, rhop_swarm, eps_dtog, &
      nu_epicycle, rp_int, rp_ext, rp_ext_width, gravx_profile, gravz_profile, &
      gravr_profile, gravx, gravz, gravr, gravsmooth, kx_gg, kz_gg, Ri0, &
      eps1, lmigration_redo, ldragforce_equi_global_eps, coeff, kx_vvp, &
      ky_vvp, kz_vvp, amplvvp, kx_xxp, ky_xxp, kz_xxp, amplxxp, kx_vpx, &
      kx_vpy, kx_vpz, ky_vpx, ky_vpy, ky_vpz, kz_vpx, kz_vpy, kz_vpz, &
      phase_vpx, phase_vpy, phase_vpz, lcoldstart_amplitude_correction, &
      particle_mesh, lparticlemesh_cic, lparticlemesh_tsc, &
      gaussian_x_frac, linterpolate_spline, &
      tstart_dragforce_par, tstart_grav_par, lparticle_gravity, &
      tstart_grav_x_par, tstart_grav_z_par,tstart_grav_r_par, taucool, &
      lcollisional_cooling_taucool, lcollisional_cooling_rms, &
      lcollisional_cooling_twobody, tausp_species, tau_coll_min, &
      ltau_coll_min_courant, coeff_restitution, tstart_collisional_cooling, &
      tausg_min, l_hole, m_hole, n_hole, &
      epsp_friction_increase,lcollisional_dragforce_cooling, ldragforce_heat, &
      lcollisional_heat, lcompensate_friction_increase, &
      lmigration_real_check, ldraglaw_epstein,ldraglaw_simple, &
      ldraglaw_epstein_stokes_linear, &
      mean_free_path_gas, ldraglaw_epstein_transonic, lcheck_exact_frontier, &
      ldraglaw_eps_stk_transonic, dustdensity_powerlaw, rad_sphere, &
      pos_sphere, ldragforce_stiff, &
      a_ellipsoid, b_ellipsoid, c_ellipsoid, pos_ellipsoid, &
      ldraglaw_steadystate, tstart_liftforce_par, &
      ldraglaw_purestokes,rpbeta_species, rpbeta, gab_width, &
      tstart_brownian_par, tstart_sink_par, &
      lbrownian_forces,lbrownian_forces_Li_Ahmadi,&
      lthermophoretic_forces,lenforce_policy, &
      interp_pol_uu,interp_pol_oo,interp_pol_TT,interp_pol_rho, &
      interp_pol_pp,interp_pol_species,brownian_T0, &
      thermophoretic_T0, lnostore_uu, ldt_grav_par, ldragforce_radialonly, &
      lsinkpoint, xsinkpoint, ysinkpoint, zsinkpoint, rsinkpoint, &
      lcoriolis_force_par, lcentrifugal_force_par, ldt_adv_par, Lx0, Ly0, &
      Lz0, lglobalrandom, lswap_radius_and_number, linsert_particles_continuously, &
      lrandom_particle_pencils, lnocalc_np, lnocalc_rhop, &
      np_const, rhop_const, particle_radius, lignore_rhop_swarm, &
      rhopmat, Deltauy_gas_friction, xp1, &
      yp1, zp1, vpx1, vpy1, vpz1, xp2, yp2, zp2, vpx2, vpy2, vpz2, &
      xp3, yp3, zp3, vpx3, vpy3, vpz3, lsinkparticle_1, rsinkparticle_1, &
      lcalc_uup, temp_grad0, thermophoretic_eq, cond_ratio, interp_pol_gradTT, &
      lreassign_strat_rhom, lparticlemesh_pqs_assignment, &
      lwithhold_init_particles, frac_init_particles, lvector_gravity, &
      birthring_r, birthring_width, lgaussian_birthring, &
      ldraglaw_stokesschiller, lbirthring_depletion, &
      remove_particle_at_time, remove_particle_criteria, &
      remove_particle_criteria_size, remove_particle_criteria_edtog, &
      lnocollapse_xdir_onecell, lnocollapse_ydir_onecell, &
      lnocollapse_zdir_onecell, qgaussz, r0gaussz, lnp_ap_as_aux,&
      lpartnucleation
!
  namelist /particles_run_pars/ &
      bcpx, bcpy, bcpz, tausp, dsnap_par_minor, beta_dPdr_dust, &
      ldragforce_gas_par, ldragforce_dust_par, np_swarm, mp_swarm, &
      rhop_swarm, eps_dtog, cdtp, cdtpgrav, lpar_spec, linterp_reality_check, &
      nu_epicycle, gravx_profile, gravz_profile, gravr_profile, gravx, gravz, &
      gravr, gravsmooth, kx_gg, kz_gg, lmigration_redo, tstart_dragforce_par, &
      tstart_grav_par, tstart_grav_x_par, cdtp_drag, &
      tstart_grav_z_par, tstart_grav_r_par, lparticle_gravity, &
      particle_mesh, lparticlemesh_cic, lparticlemesh_tsc, taucool, &
      lcollisional_cooling_taucool, lcollisional_cooling_rms, &
      lcollisional_cooling_twobody, lcollisional_dragforce_cooling, &
      tau_coll_min, ltau_coll_min_courant, coeff_restitution, &
      tstart_collisional_cooling, tausg_min, epsp_friction_increase, &
      ldragforce_heat, lcollisional_heat, lcompensate_friction_increase, &
      lmigration_real_check,ldraglaw_variable, luse_tau_ap, &
      ldraglaw_epstein, ldraglaw_simple, rdiffconst_pass, &
      ldraglaw_epstein_stokes_linear, mean_free_path_gas, &
      ldraglaw_epstein_transonic, lcheck_exact_frontier, &
      ldraglaw_eps_stk_transonic, ldragforce_stiff, &
      ldraglaw_variable_density, ldraglaw_steadystate, tstart_liftforce_par, &
      ldraglaw_purestokes, ldiffuse_dragf, ldiffuse_passive, rdiffconst_dragf, &
      tstart_brownian_par, tstart_sink_par, ldiff_dragf, ldiff_pass, &
      lbrownian_forces, lbrownian_forces_Li_Ahmadi, lenforce_policy, &
      interp_pol_uu,interp_pol_oo, interp_pol_TT, interp_pol_rho, &
      interp_pol_pp,interp_pol_species, &
      brownian_T0,thermophoretic_T0, lnostore_uu, ldt_grav_par, &
      ldragforce_radialonly, lsinkpoint, xsinkpoint, ysinkpoint, zsinkpoint, &
      rsinkpoint, lshear_accel_par, lcoriolis_force_par, &
      lcentrifugal_force_par, ldt_adv_par, &
      linsert_particles_continuously, particles_insert_rate, &
      max_particle_insert_time, lrandom_particle_pencils, lnocalc_np, &
      lnocalc_rhop, lstocunn1, istep_dragf, istep_pass, &
      np_const, rhop_const, particle_radius, lprecalc_cell_volumes, &
      Deltauy_gas_friction, loutput_psize_dist, log_ap_min_dist, &
      log_ap_max_dist, nbin_ap_dist, lsinkparticle_1, rsinkparticle_1, &
      lthermophoretic_forces, temp_grad0, &
      thermophoretic_eq, cond_ratio, interp_pol_gradTT, lcommunicate_rhop, &
      lcommunicate_np, lcylindrical_gravity_par, lignore_rhop_swarm, &
      l_shell, k_shell, lparticlemesh_pqs_assignment, pscalar_sink_rate, &
      lpscalar_sink, lsherwood_const, lnu_draglaw, nu_draglaw,lbubble, &
      rpbeta_species, rpbeta, gab_width, initxxp, initvvp, &
      particles_insert_ramp_time, tstart_insert_particles, birthring_r, &
      birthring_width, lsimple_volume, &
      lgaussian_birthring, tstart_rpbeta, linsert_as_many_as_possible, &
      lvector_gravity, lcompensate_sedimentation,compensate_sedimentation, &
      lpeh_radius, A3, A2, ldraglaw_stokesschiller, lbirthring_depletion, birthring_lifetime, &
      remove_particle_at_time, remove_particle_criteria, remove_particle_criteria_size, &
      remove_particle_criteria_edtog, &
      ascalar_ngp, ascalar_cic, rp_int, rp_ext, rp_ext_width, lnpmin_exclude_zero, &
      lcondensation_rate, vapor_mixing_ratio_qvs, lfollow_gas, &
      ltauascalar, rhoa, G_condensation, lpartnucleation, nucleation_threshold, &
      redfrac, lset_df_insert_nucleii, it_insert_nuclei
!
  integer :: idiag_xpm=0, idiag_ypm=0, idiag_zpm=0      ! DIAG_DOC: $x_{part}$
  integer :: idiag_xpmin=0, idiag_ypmin=0, idiag_zpmin=0      ! DIAG_DOC: $x_{part}$
  integer :: idiag_xpmax=0, idiag_ypmax=0, idiag_zpmax=0      ! DIAG_DOC: $x_{part}$
  integer :: idiag_xp2m=0, idiag_yp2m=0, idiag_zp2m=0   ! DIAG_DOC: $x^2_{part}$
  integer :: idiag_vrelpabsm=0                          ! DIAG_DOC: $\rm{Absolute value of mean relative velocity}$
  integer :: idiag_rpm=0, idiag_rp2m=0
  integer :: idiag_vpxm=0, idiag_vpym=0, idiag_vpzm=0   ! DIAG_DOC: $u_{part}$
  integer :: idiag_vpx2m=0, idiag_vpy2m=0, idiag_vpz2m=0 ! DIAG_DOC: $u^2_{part}$
  integer :: idiag_ekinp=0     ! DIAG_DOC: $E_{kin,part}$
!  integer :: idiag_vtherm500=0
  integer :: idiag_vpxmax=0, idiag_vpymax=0, idiag_vpzmax=0, idiag_vpmax=0 ! DIAG_DOC: $MAX(u_{part})$
  integer :: idiag_vpxmin=0, idiag_vpymin=0, idiag_vpzmin=0    ! DIAG_DOC: $MIN(u_{part})$
  integer :: idiag_urel=0
  integer :: idiag_vpxvpym=0, idiag_vpxvpzm=0, idiag_vpyvpzm=0
  integer :: idiag_rhopvpxm=0, idiag_rhopvpym=0, idiag_rhopvpzm=0
  integer :: idiag_rhopvpxt=0, idiag_rhopvpyt=0, idiag_rhopvpzt=0
  integer :: idiag_rhopvpysm=0, idiag_rhopart=0
  integer :: idiag_lpxm=0, idiag_lpym=0, idiag_lpzm=0
  integer :: idiag_lpx2m=0, idiag_lpy2m=0, idiag_lpz2m=0
  integer :: idiag_npm=0, idiag_np2m=0, idiag_npmax=0, idiag_npmin=0 ! DIAG_DOC: $\rm{mean particle number density}$
  integer :: idiag_dtdragp=0
  integer :: idiag_nparmin=0, idiag_nparmax=0, idiag_npargone=0, idiag_npar_found=0
  integer :: idiag_nparsum=0
  integer :: idiag_rhopm=0, idiag_rhoprms=0, idiag_rhop2m=0, idiag_rhopmax=0
  integer :: idiag_rhopmin=0, idiag_decollp=0, idiag_rhopmphi=0
  integer :: idiag_epspmin=0, idiag_epspmax=0, idiag_epspm=0
  integer :: idiag_npmx=0, idiag_npmy=0, idiag_npmz=0
  integer :: idiag_rhopmx=0, idiag_rhopmy=0, idiag_rhopmz=0
  integer :: idiag_rhop2mx=0, idiag_rhop2my=0, idiag_rhop2mz=0
  integer :: idiag_epspmx=0, idiag_epspmy=0, idiag_epspmz=0
  integer :: idiag_mpt=0, idiag_dedragp=0, idiag_rhopmxy=0, idiag_rhopmr=0
  integer :: idiag_sigmap=0
  integer :: idiag_rpvpxmx = 0, idiag_rpvpymx = 0, idiag_rpvpzmx = 0
  integer :: idiag_rpvpx2mx = 0, idiag_rpvpy2mx = 0, idiag_rpvpz2mx = 0
  integer :: idiag_dvpx2m=0, idiag_dvpy2m=0, idiag_dvpz2m=0
  integer :: idiag_dvpm=0, idiag_dvpmax=0, idiag_epotpm=0
  integer :: idiag_rhopmxz=0, idiag_nparpmax=0, idiag_npmxy=0
  integer :: idiag_eccpxm=0, idiag_eccpym=0, idiag_eccpzm=0
  integer :: idiag_eccpx2m=0, idiag_eccpy2m=0, idiag_eccpz2m=0
  integer :: idiag_vprms=0, idiag_vpyfull2m=0, idiag_deshearbcsm=0
  integer :: idiag_Shm=0, idiag_latentheatm=0
  integer :: idiag_ffcondposm=0, idiag_ffcondnegm=0, idiag_ffcondm=0
  integer, dimension(ndustrad) :: idiag_npvzmz=0, idiag_npvz2mz=0, idiag_nptz=0
  integer, dimension(ndustrad) :: idiag_npuzmz=0
!
  real, dimension(:), pointer :: beta_glnrho_global, beta_glnrho_scaled
  real, dimension(nx) :: dt1_drag=0., urel_sum, drag_heat
!
  contains
!***********************************************************************
    subroutine register_particles
!
!  Set up indices for access to the fp and dfp arrays
!
!  29-dec-04/anders: coded
!
      use FArrayManager, only: farray_register_auxiliary, farray_index_append
      use Particles_caustics, only: register_particles_caustics
      use Particles_tetrad, only: register_particles_tetrad
      use Particles_potential, only: register_particles_potential
      use SharedVariables, only: put_shared_variable
!
      integer ind_tmp
!
      if (lroot) call svn_id( &
          "$Id: particles_dust.f90,v 1.1 2018/08/24 15:48:10 wlyra Exp $")
!
!  Indices for particle position.
!
      call append_npvar('ixp',ixp)
      call append_npvar('iyp',iyp)
      call append_npvar('izp',izp)
!
!  Indices for particle velocity.
!
      call append_npvar('ivpx',ivpx)
      call append_npvar('ivpy',ivpy)
      call append_npvar('ivpz',ivpz)
!
!  Set indices for particle assignment.
!
      if (.not. lnocalc_np) call farray_register_auxiliary('np',inp,communicated = lcommunicate_np)
      call put_shared_variable('inp', inp, caller='register_particles')
      if (lpeh_radius) then
        lnocalc_rhop = .false.
        call farray_register_auxiliary('peh',ipeh,communicated=.false.)
      endif
!
!  Store the number of particles per grid cell if requested.
!
      if (np_ap_spec.or.lnp_ap_as_aux) then
        call farray_register_auxiliary('np_ap',ind_tmp,array=ndustrad)
        iapn(1) = ind_tmp
        iapn = iapn(1) + indgen(ndustrad) - 1
        call farray_index_append('n_np_ap',ndustrad)
      end if
      call put_shared_variable('iapn', iapn)
!
!  Rhop
!
      if (.not. lnocalc_rhop) call farray_register_auxiliary('rhop',irhop, &
          communicated = lparticles_sink .or. lcommunicate_rhop)
      call put_shared_variable('irhop', irhop)
      if (lcalc_uup .or. ldragforce_stiff) then
        call farray_register_auxiliary('uup',iuup,communicated=.true.,vector=3)
        iupx = iuup
        iupy = iuup+1
        iupz = iuup+2
      endif
!
!  Register nucleation radius
!
      if (lpartnucleation) then
        call farray_register_auxiliary('nucl_rmin',inucl,communicated=.false.)
        call farray_register_auxiliary('nucl_rate',inucrate,communicated=.false.)
        call farray_register_auxiliary('supersat',isupsat,communicated=.false.)
        call append_npvar('born',iborn)
     endif
!
!  Special variable for stiff drag force equations.
!
      if (ldragforce_stiff) then
        call farray_register_auxiliary('ffg',iffg,communicated=.true.,vector=3)
        ifgx = iffg
        ifgy = iffg+1
        ifgz = iffg+2
      endif
!
!  Special variable for diffusion of the passive scalar consumption
!
      if (ldiffuse_passive .and. lpscalar_sink) then
        call farray_register_auxiliary('dlncc',idlncc,communicated=.true.)
      elseif (ldiffuse_passive .and. .not. lpscalar_sink) then
        call fatal_error('particles_dust','ldiffuse_passive needs lpscalar_sink')
      endif
!
!  Special variable for diffusion of particle-gas dragforce
!
      if (ldiffuse_dragf .and. ldragforce_gas_par) then
        call farray_register_auxiliary('dfg',idfg,communicated=.true.,vector=3)
        idfx = idfg
        idfy = idfg+1
        idfz = idfg+2
      elseif (ldiffuse_dragf .and. .not. ldragforce_gas_par) then
        call fatal_error('particles_dust','ldiffuse_dragf needs ldragforce_gas_par')
      endif
!
!  Relaxation time of supersaturation
      if (lascalar) then
        if (ltauascalar) call farray_register_auxiliary('tauascalar', itauascalar)
        call farray_register_auxiliary('condensationRate', icondensationRate)
        call farray_register_auxiliary('waterMixingRatio', iwaterMixingRatio)
      endif
!
!  Share friction time (but only if Epstein drag regime!).
!
      call put_shared_variable("ldragforce_gas_par", ldragforce_gas_par)
!
      if (ldraglaw_epstein .or. ldraglaw_simple) then
        call put_shared_variable("tausp", tausp)
        call put_shared_variable("tausp_species", tausp_species)
        call put_shared_variable("tausp1_species", tausp1_species)
      endif
!
!  Share Keplerian gravity.
!
      call put_shared_variable('gravr',gravr,caller = 'initialize_particles')

      if (l_shell) then
        call put_shared_variable('uup_shared',uup_shared)
        call put_shared_variable('vel_call',vel_call)
        call put_shared_variable('turnover_call',turnover_call)
        call put_shared_variable('turnover_shared',turnover_shared)
      endif
!
!  Shared variables
!
      if (lchemistry) then
        call put_shared_variable('lpartnucleation',lpartnucleation)
      endif
!
      if (lascalar) call put_shared_variable('G_condensation',G_condensation)
!
!  Kill particles that spend enough time in birth ring
!
      if (lbirthring_depletion) call append_npvar('ibrtime',ibrtime)
!
!  If we are using caustics, here we should call the correspoding register equation:
!
      if (lparticles_caustics) call register_particles_caustics()
!
!  If we are using tetrad, here we should call the correspoding register equation:
!
      if (lparticles_tetrad) call register_particles_tetrad()
!
!  If we are using particles_potential, here we should call the correspoding register equation:
!
      if (lparticles_potential) call register_particles_potential()
!
    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles(f,fp)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-dec-04/anders: coded
!   5-mar-15/MR: reference state included in calculation of mean density
!
      use EquationOfState, only: rho0, cs0
      use SharedVariables, only: get_shared_variable
      use Particles_caustics, only: initialize_particles_caustics
      use Particles_tetrad, only: initialize_particles_tetrad
      use Particles_potential, only: initialize_particles_potential
      use General, only: random_gen, rtoa
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
!
      real :: rhom
      integer :: jspec
      logical, pointer :: lshearadvection_as_shift
      real, pointer :: reference_state_mass
      integer :: ik, icnt
      real :: xxp, yyp, zzp
!
!  This module is incompatible with particle block domain decomposition.
!
      if (lparticles_blocks) &
        call fatal_error('initialize_particles','must use PARTICLES = PARTICLES_DUST_BLOCKS'// &
                         'with particle block domain decomposition')
!
!  Due to an un-resolved problems with persistent and seeds, the default random 
!  number generator can not be used if linsert_particles_continuously=T
!  together with lpersist=T and nproc>1.
!  I have not understood exactly what causes the problem, but it is
!  related to a change in seed on the processor(s) where the particles
!  are being inserted continously. Such a change in seed, triggers the
!  code to save the new seed if lpersist=T (I think) on this/these processors.
!  This results in var.dat files that do not have the same size for all
!  processors, and hence an error occures when we try to read these files.
!
      if (lpersist .and. linsert_particles_continuously .and. (ncpus .gt. 1)) then
        if (random_gen .eq. "min_std") &
          call fatal_error('initialize_particles','cannot use random_gen="min_std" when '// &
                           'linsert_particles_continuously=T and lpersist=T')
      endif
       if (lpersist .and. lpartnucleation .and. (ncpus .gt. 1)) then
        if (random_gen .eq. "min_std") &
          call fatal_error('initialize_particles','cannot use random_gen="min_std" when '// &
                           'lpartnucleation=T and lpersist=T')
      endif
!
! Mark particles in a pencil
! We make a pencil along the x-direction
!
!      ntagged_local = 0
!      do ik=1,npar_loc
!        xxp=fp(ik,ixp);yyp=fp(ik,iyp);zzp=fp(ik,izp)
!        if ((yyp .lt. ypmax) .and.   (yyp .lt. ypmin) ) then
!          if ((zzp .lt. zpmax) .and.   (zzp .lt. zpmin) ) then
!           ntagged_local(iproc)=ntagged_local(iproc)+1
!          endif
!        endif
!      enddo
!
! Sum ntagged_local from each processor to get ntagged_global
!                        call mpiallreduce_sum_int_arr(ntagged_local,ntagged_global,ncpus)
!                        itagged_start = sum(ntagged_global(1:iproc-1))+1
! AKS
! Then at each processor, sum over ntagged_global to obtain npar_tagged
!                        npar_tagged = sum(ntagged_global)
! Then allocate the taggedp_list which as dimension npar_tagged
!                        allocate(taggedp_list_local(npar_tagged)); taggedp_list_local = 0
!                        allocate(taggedp_list_global(npar_tagged)); taggedp_list_global = 0
! Then fill up the taggedp_list
!                        icnt = 0
!      do ik=1,npar_loc
!        xxp=fp(ik,ixp);yyp=fp(ik,iyp);zzp=fp(ik,izp)
!        if ((yyp .lt. ypmax) .and. (yyp .lt. ypmin) ) then
!          if ((zzp .lt. zpmax) .and. (zzp .lt. zpmin) ) then
!           taggedp_list_local(itagged_start+icnt) = ipar(ik)
!                                         icnt = icnt+1
!          endif
!        endif
!      enddo
!
! Each processor writes on different locations at this array.
! Then communicate this array over all processors such that they all have it.
!                        call mpiallreduce_sum_int_arr(taggedp_list_local,taggedp_list_global,npar_tagged)
!
!
!  Hand over Coriolis force and shear acceleration to Particles_drag.
!
      if (lparticles_drag) then
        if (lcoriolis_force_par) then
          lcoriolis_force_par = .false.
          if (lroot) print *, 'initialize_particles: turned off and hand over Coriolis force to Particles_drag. '
        endif
        if (lshear .and. lshear_accel_par) then
          lshear_accel_par = .false.
          if (lroot) print *, 'initialize_particles: turned off and hand over shear acceleration to Particles_drag. '
        endif
      endif
!
!  Check if shear advection is on and decide if it needs to be included in the timestep condition.
!
      if (lshear) then
        call get_shared_variable('lshearadvection_as_shift', lshearadvection_as_shift, caller='initialize_particles')
        lcdtp_shear = .not. lshearadvection_as_shift
        nullify(lshearadvection_as_shift)
      endif
!
!  Report the particle radius, if set.
!
      if (particle_radius /= 0.0) then
        if (lparticles_radius) &
          call fatal_error('initialize_particles','particle_radius /= 0 has no effect when Particles_radius is on')
        if (lroot) print *, 'initialize_particles: radius of each constituent particle = ', particle_radius
      endif
!
!  The inverse stopping time is needed for drag force and collisional cooling.
!
      if (tausp /= 0.0) tausp1 = 1/tausp
!
!  Inverse cooling time.
!
      if (taucool /= 0.0) taucool1 = 1/taucool
!
!  Inverse material density.
!
      if (rhopmat /= 0.0) rhopmat1 = 1/rhopmat
!
!  Multiple dust species. Friction time is given in the array tausp_species.
!
      if (npar_species > 1) then
        if (lroot) then
          print*, 'initialize_particles: '//'Number of particle species = ', npar_species
          print*, 'initialize_particles: tausp_species = ', tausp_species
        endif
!
!  Must have set tausp_species for drag force.
!
        if (ldragforce_dust_par .or. ldragforce_gas_par) then
          if (any(tausp_species == 0.0)) &
            call fatal_error('initialize_particles','drag force must have tausp_species/=0.0')
!
!  Inverse friction time is needed for drag force.
!
          where (tausp_species /= 0.0) tausp1_species=1/tausp_species
        endif
      else
!
!  Single dust species => If tausp_species is set, it is probably an error.
!
        if (any(tausp_species /= 0.0) .and. tausp /= tausp_species(1)) then
          call fatal_error('initialize_particles', 'When there is only '// &
                           'one particle species, use tausp instead of tausp_species')
        endif
        tausp_species(1) = tausp
        if (tausp_species(1) /= 0.0) tausp1_species(1) = 1/tausp_species(1)
      endif
!
!  Global gas pressure gradient seen from the perspective of the dust.
!
      if (beta_dPdr_dust /= 0.0) then
        beta_dPdr_dust_scaled = beta_dPdr_dust*Omega/cs0
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
      if (eps_dtog > 0.0) then
! For stratification, take into account gas present outside the simulation box.
        if (lreassign_strat_rhom .and.((lgravz.and.lgravz_gas) .or. gravz_profile == 'linear')) then
          ! rhom = (total mass) / (box volume) = Sigma / Lz
          ! Sigma = sqrt(2pi) * rho0 * H
          !   rho0 = mid-plane density, H = (sound speed) / (epicycle frequency)
          rhom = sqrt(2.0 * pi) / Lz
          if (nu_epicycle > 0.0) rhom = rhom * (rho0 * cs0 / nu_epicycle)
        else
          rhom = rho0
          if (lreference_state) then
            call get_shared_variable('reference_state_mass',reference_state_mass, &
                 caller = 'initialize_particles')
            rhom = rhom+reference_state_mass/box_volume
          endif
        endif
        if (rhop_swarm == 0.0) rhop_swarm = eps_dtog*rhom/(real(npar)/nwgrid)
        if (mp_swarm == 0.0) mp_swarm   = eps_dtog*rhom*box_volume/(real(npar))
        if (lroot) print*, 'initialize_particles: dust-to-gas ratio eps_dtog=', eps_dtog
      endif
!
      if (lroot) then
        print*, 'initialize_particles: mass per constituent particle mpmat=', mpmat
        print*, 'initialize_particles: mass per superparticle mp_swarm =', mp_swarm
        print*, 'initialize_particles: number density per superparticle np_swarm=', np_swarm
        print*, 'initialize_particles: mass density per superparticle rhop_swarm=', rhop_swarm
      endif
!
!  Calculate nu_epicycle**2 for gravity.
!
      if (gravz_profile == '' .and. nu_epicycle /= 0.0) gravz_profile = 'linear'
      nu_epicycle2 = nu_epicycle**2
!
!  Calculate gravsmooth**2 for gravity.
!
      if (gravsmooth /= 0.0) gravsmooth2 = gravsmooth**2
!
!  Inverse of minimum gas friction time (time-step control).
!
      if (tausg_min /= 0.0) tausg1_max = 1.0/tausg_min
!
!  Set minimum collisional time-scale so that time-step is not affected.
!
      if (lrun .and. ltau_coll_min_courant) then
        if (cs0 == impossible) then
          tau_coll_min = impossible
        else
          tau_coll_min = 2*dx/cs0
        endif
        if (lroot) print*, 'initialize particles: set minimum collisional '// &
                           'time-scale equal to two times the Courant time-step.'
      endif
!
!  Inverse of minimum collisional time-scale.
!
      if (lrun .and. tau_coll_min > 0.0) tau_coll1_max = 1/tau_coll_min
!
!  Gas density is needed for back-reaction friction force.
!
      if (ldragforce_gas_par .and. .not. ldensity) &
        call fatal_error('initialize_particles','friction force on gas only works together with gas density')
!
!  Need to map particles on the grid for dragforce on gas.
!
      if (ldragforce_gas_par) then
!
!  When drag force is smoothed, df is also set in the first ghost zone. This
!  region needs to be folded back into the df array after pde is finished,
!
        if (lparticlemesh_cic .or. lparticlemesh_tsc) lfold_df = .true.
      endif
!
      if (lcollisional_cooling_twobody) then
        allocate(kneighbour(mpar_loc))
        lshepherd_neighbour = .true.
      endif
!
      if (ldraglaw_epstein_stokes_linear) ldraglaw_epstein = .false.
      if (ldraglaw_epstein_transonic    .or. &
          ldraglaw_eps_stk_transonic    .or. &
          ldraglaw_steadystate          .or. &
          ldraglaw_purestokes          .or. &
          ldraglaw_stokesschiller      .or. &
          ldraglaw_simple) then
        ldraglaw_epstein = .false.
      endif
      if (ldraglaw_epstein_transonic .and. ldraglaw_eps_stk_transonic) then
        print*,'both epstein and epstein-stokes transonic '// &
            'drag laws are switched on. You cannot have '// &
            'both. Stop and choose only one.'
        call fatal_error('initialize_particles','')
      endif
      if (ldraglaw_steadystate .and. mean_free_path_gas==0. .and. .not.lstocunn1) &
        call fatal_error('initialize_particles','mean_free_path_gas==0. for ldraglaw_steadystate=T, lstocunn1=F')
!
!  Stiff drag force approximation.
!
      if (ldragforce_stiff) then
        if (ldragforce_dust_par .or. ldragforce_gas_par) &
          call fatal_error('initialize_particles','stiff drag force approximation incompatible with normal drag')
        f(l1:l2,m1:m2,n1:n2,ifgx:ifgz) = 0.0
      endif
!
!  Initialize storage of energy gain released by shearing boundaries.
!
      if (idiag_deshearbcsm /= 0) energy_gain_shear_bcs = 0.0
!
!  Drag force on gas right now assumed rhop_swarm is the same for all particles.
!
! NILS: Commented out the check below (after discussion with Anders)
! NILS: as it does not seem to be required any longer.
!
!      if (ldragforce_gas_par.and.(lparticles_radius.or.lparticles_number) &
!          .and..not.lparticles_density) then
!        if (lroot) print*, 'initialize_particles: drag force on gas is '// &
!            'not yet implemented for variable particle radius or number'
!        call fatal_error('initialize_particles','')
!      endif
!
!  Fatal error if sink particle radius is zero or negative.
!
      if (lsinkparticle_1 .and. rsinkparticle_1 <= 0.0) &
        call fatal_error('initialize_particles','sink particle radius is <=0:'//trim(rtoa(rsinkparticle_1)))
!
!  Particle self gravity for x,z,r direction
!
      if ( tstart_grav_par > 0.0) then
        tstart_grav_x_par = tstart_grav_par
        tstart_grav_z_par = tstart_grav_par
        tstart_grav_r_par = tstart_grav_par
      endif
!
!  Die if lcompensate_sedimentation is used and the vertical direction is present.
!
      if (lcompensate_sedimentation .and. &
          ( (nzgrid /= 1 .and. (lcartesian_coords .or. lcylindrical_coords)).or. &
            (nygrid /= 1 .and. lspherical_coords) &
          )) call fatal_error("initialize_particles", &
          "compensate_sedimentation can only be used when vertical dimension is not present")
!
!  Set up interpolation logicals. These logicals can be OR'ed with some logical
!  in the other particle modules' initialization subroutines to enable
!  interpolation based on some condition local to that module.
!  (The particles_spin module will for instance enable interpolation of the
!  vorticity oo)
!
      if (lnostore_uu) then
        if (ldraglaw_steadystate .or. lparticles_spin .or. lsolid_ogrid) &
            call fatal_error('initialize_particles', 'lnostore_uu = F required')
        interp%luu = .false.
      else
        interp%luu = ldragforce_dust_par .or. ldraglaw_steadystate .or. lparticles_spin .or. lsolid_ogrid
      endif
      interp%loo = .false.
      interp%lTT = (lbrownian_forces .and.(brownian_T0 == 0.0)) &
          .or.(lthermophoretic_forces .and.(thermophoretic_T0 == 0.0)) &
          .or. lparticles_temperature
      interp%lgradTT = lthermophoretic_forces .and. (temp_grad0(1) == 0.0) &
          .and. (temp_grad0(2) == 0.0) .and. (temp_grad0(3) == 0.0)
      interp%lrho = lbrownian_forces .or. ldraglaw_steadystate &
          .or. lthermophoretic_forces .or. ldraglaw_purestokes &
          .or. ldraglaw_stokesschiller .or. lsolid_ogrid
      interp%lcs = lbrownian_forces .or. ldraglaw_epstein
      interp%lnu = lchemistry
      interp%lpp = lparticles_chemistry
      interp%lspecies = lparticles_surfspec
!
!  Determine interpolation policies:
!  Make sure that interpolation of uu is chosen in a backwards compatible
!  manner. NGP is chosen by default.
!
      if (.not. lenforce_policy) then
        if (lparticlemesh_cic) then
          interp_pol_uu = 'cic'
        elseif (lparticlemesh_tsc) then
          interp_pol_uu = 'tsc'
        endif
      endif
!
!  Overwrite with new policy variables:
!
      select case (interp_pol_uu)
      case ('tsc')
        interp%pol_uu = tsc
      case ('cic')
        interp%pol_uu = cic
      case ('ngp')
        interp%pol_uu = ngp
      case default
        call fatal_error('initialize_particles','no such such interp_pol_uu: '//trim(interp_pol_uu))
      endselect
!
      select case (interp_pol_oo)
      case ('tsc')
        interp%pol_oo = tsc
      case ('cic')
        interp%pol_oo = cic
      case ('ngp')
        interp%pol_oo = ngp
      case default
        call fatal_error('initialize_particles','no such such interp_pol_oo: '//trim(interp_pol_oo))
      endselect
!
      select case (interp_pol_TT)
      case ('tsc')
        interp%pol_TT = tsc
      case ('cic')
        interp%pol_TT = cic
      case ('ngp')
        interp%pol_TT = ngp
      case default
        call fatal_error('initialize_particles','no such interp_pol_TT: '//trim(interp_pol_TT))
      endselect
!
      select case (interp_pol_gradTT)
      case ('tsc')
        call not_implemented('initialize_particles','interp_pol_gradTT: '//trim(interp_pol_gradTT))
      case ('cic')
        call not_implemented('initialize_particles','interp_pol_gradTT: '//trim(interp_pol_gradTT))
      case ('ngp')
        interp%pol_gradTT = ngp
      case default
        call fatal_error('initialize_particles','no such interp_pol_gradTT: '//trim(interp_pol_gradTT))
      endselect
!
      select case (interp_pol_rho)
      case ('tsc')
        interp%pol_rho = tsc
      case ('cic')
        interp%pol_rho = cic
      case ('ngp')
        interp%pol_rho = ngp
      case default
        call fatal_error('initialize_particles','no such interp_pol_rho: '//trim(interp_pol_rho))
      endselect
!
      select case (interp_pol_pp)
      case ('tsc')
        call not_implemented('initialize_particles','interp_pol_pp: '//trim(interp_pol_pp))
      case ('cic')
        call fatal_error('initialize_particles','interp_pol_pp: '//trim(interp_pol_pp))
      case ('ngp')
        interp%pol_pp = ngp
      case default
        call fatal_error('initialize_particles','no such such interp_pol_pp: '//trim(interp_pol_pp))
      endselect
!
      select case (interp_pol_nu)
      case ('tsc')
        call not_implemented('initialize_particles','interp_pol_nu: '//trim(interp_pol_nu))
      case ('cic')
        call not_implemented('initialize_particles','interp_pol_nu: '//trim(interp_pol_nu))
      case ('ngp')
        interp%pol_nu = ngp
      case default
        call fatal_error('initialize_particles','no such interp_pol_nu: '//trim(interp_pol_nu))
      endselect
!
      if (l_shell) then
        if ( k_shell < 0) call fatal_error('initialize_particles','set k_shell >= 0')
      endif
!
!  Write constants to disk.
!
     if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position='append')
        write (1,*) 'np_swarm=', np_swarm
        write (1,*) 'mpmat=', 0.0
        write (1,*) 'mp_swarm=', mp_swarm
        write (1,*) 'rhop_swarm=', rhop_swarm
        close (1)
     endif
!
!  if we are using caustics:
!
     if (lparticles_caustics) call initialize_particles_caustics(f)
!
!  if we are using tetrad:
!
      if (lparticles_tetrad) call initialize_particles_tetrad(f)
!
!  if we are using particles_potential:
!
      if (lparticles_potential) call initialize_particles_potential(fp)
!
!  Variables needed for corotational frame on particles.
!
      if (lgravr .and. lcorotational_frame) then
        call get_shared_variable('g1',g1)
        call get_shared_variable('rp1',rp1)
        call get_shared_variable('rp1_smooth',rp1_smooth)
        call get_shared_variable('t_ramp_mass',t_ramp_mass)
        call get_shared_variable('t_start_secondary',t_start_secondary)
        call get_shared_variable('lramp_mass',lramp_mass)
        call get_shared_variable('lsecondary_wait',lsecondary_wait)
      endif
!
      if (lparticles_grad.and.iguij==0) &
        call fatal_error('initialize_particles','particles_grad demands iguij /= 0')
!
      if (lthermophoretic_forces.and.cond_ratio==0.0) & 
        call fatal_error('initialize_particles','Cond ratio=0 with thermophoretic force')

      if (any(gravr_profile==(/'newtonian-central','newtonian        '/).and.lpointmasses.and. &
          lparticle_gravity)) &
        call fatal_error('initialize_particles','you are using massive particles. '// &
             achar(10)//'The N-body code should take care of the stellar-like gravity on the dust.'// &
             achar(10)//'Switch off the gravr_profile="newtonian*" in particles_init')
!
      if (ldiffuse_passive.and.ilncc == 0) &
        call fatal_error('initialize_particles','ldiffuse_passive needs pscalar_nolog=F')
!
      if (ldragforce_dust_par .and. lparticlemesh_cic .and. &
          lpscalar_sink .and. lpscalar .and. ilncc == 0) &
        call fatal_error('initialize_particles','lpscalar_sink not allowed for pscalar_nolog')
!
      if (lchemistry) then
        if (lpartnucleation .or. lcondensing_species) then
          call get_shared_variable('true_density_cond_spec',&
               true_density_cond_spec)
          if (rhopmat /= 1.0 .and. rhopmat /= true_density_cond_spec) &
               call fatal_error("initialize_particles",&
               "Can not set rhopmat if true_density_cond_spec is given.")
          rhopmat=true_density_cond_spec
        endif
      endif
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f,fp,ineargrid)
!
!  Initial positions and velocities of dust particles.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: cs20, rho0
      use General, only: random_number_wrapper, normal_deviate
      use Mpicomm, only: mpireduce_sum, mpibcast_real
      use InitialCondition, only: initial_condition_xxp, initial_condition_vvp
      use Particles_diagnos_dv, only: repeated_init
      use Particles_caustics, only: init_particles_caustics
      use Particles_tetrad, only: init_particles_tetrad
      use SharedVariables, only: get_shared_variable
!
      real, dimension(mx,my,mz,mfarray), intent(out) :: f
      real, dimension(mpar_loc,mparray), intent(out) :: fp
      integer, dimension(mpar_loc,3), intent(out) :: ineargrid
      real, dimension(mpar_loc) :: rr_tmp, az_tmp
!
      real, dimension(3) :: uup, Lxyz_par, xyz0_par, xyz1_par
      real :: vpx_sum, vpy_sum, vpz_sum
      real :: r, p, q, px, py, pz, eps, cs, k2_xxp, rp2
      real :: dim1, dx_par, dy_par, dz_par
      real :: rad, rad_scl, phi, tht, tmp, OO, xx0, yy0, r2
      real :: rpar_int, rpar_ext, tausp_par
      integer :: npar_loc_x, npar_loc_y, npar_loc_z
      integer :: l, j, k, ix0, iy0, iz0
      logical :: lequidistant=.false.
!     
      call get_shared_variable('beta_glnrho_global',beta_glnrho_global,caller='init_particles')
      call get_shared_variable('beta_glnrho_scaled',beta_glnrho_scaled)
!
!  Optionally withhold some number of particles, to be inserted in
!  insert_particles. The particle indices to be removed are not randomized,
!  so any randomization needs to be taken care of in the above initxxp cases.
!
      if (lwithhold_init_particles .and. frac_init_particles < 1.0-tini) then
        npar_loc = nint(frac_init_particles*real(npar_loc))
      endif
!
!  Use either a local random position or a global random position for certain
!  initial conditions. The default is a local random position, but the equal
!  number of particles per processors means that this is not completely random.
!
      if (lglobalrandom) then
        Lxyz_par = Lxyz
        xyz0_par = xyz0
        xyz1_par = xyz1
      else
        Lxyz_par = Lxyz_loc
        xyz0_par = xyz0_loc
        xyz1_par = xyz1_loc
      endif
!
!  Initial particle position.
!
      do j = 1,ninit
!
        select case (initxxp(j))
!
        case ('nothing')
          if (lroot .and. j == 1) print*, 'init_particles: nothing'
!
        case ('origin')
          if (lroot) print*, 'init_particles: All particles at origin'
          fp(1:npar_loc,ixp:izp) = 0.0
!
        case ('zero-z')
          if (lroot) print*, 'init_particles: Zero z coordinate'
          fp(1:npar_loc,izp) = 0.0
!
        case ('constant')
          if (lroot) &
              print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
          fp(1:npar_loc,ixp) = xp0
          fp(1:npar_loc,iyp) = yp0
          fp(1:npar_loc,izp) = zp0
!
        case ('constant-1')
          if (lroot) &
              print*, 'init_particles: Particle 1 at x,y,z=', xp1, yp1, zp1
          do k = 1,npar_loc
            if (ipar(k) == 1) then
              fp(k,ixp) = xp1
              fp(k,iyp) = yp1
              fp(k,izp) = zp1
            endif
          enddo
!
        case ('constant-2')
          if (lroot) &
              print*, 'init_particles: Particle 2 at x,y,z=', xp2, yp2, zp2
          do k = 1,npar_loc
            if (ipar(k) == 2) then
              fp(k,ixp) = xp2
              fp(k,iyp) = yp2
              fp(k,izp) = zp2
            endif
          enddo
!
        case ('constant-3')
          if (lroot) &
              print*, 'init_particles: Particle 2 at x,y,z=', xp3, yp3, zp3
          do k = 1,npar_loc
            if (ipar(k) == 3) then
              fp(k,ixp) = xp3
              fp(k,iyp) = yp3
              fp(k,izp) = zp3
            endif
          enddo
!
        case ('random-constz')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k = 1,npar_loc
            if (nxgrid /= 1) then
              call random_number_wrapper(r)
              fp(k,ixp) = r
            endif
            if (nygrid /= 1) then
              call random_number_wrapper(r)
              fp(k,iyp) = r
            endif
          enddo
          if (nxgrid /= 1) &
              fp(1:npar_loc,ixp) = xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          if (nygrid /= 1) &
              fp(1:npar_loc,iyp) = xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          if (nzgrid /= 1) fp(1:npar_loc,izp) = zp0
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k = 1,npar_loc
            if (nxgrid /= 1 .or. lnocollapse_xdir_onecell) then
              call random_number_wrapper(r)
              fp(k,ixp) = r
            endif
            if (nygrid /= 1 .or. lnocollapse_ydir_onecell) then
              call random_number_wrapper(r)
              fp(k,iyp) = r
            endif
            if (nzgrid /= 1 .or. lnocollapse_zdir_onecell) then
              call random_number_wrapper(r)
              fp(k,izp) = r
            endif
          enddo
          if (nxgrid /= 1 .or. lnocollapse_xdir_onecell) &
              fp(1:npar_loc,ixp) = xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          if (nygrid /= 1 .or. lnocollapse_ydir_onecell) &
              fp(1:npar_loc,iyp) = xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          if (nzgrid /= 1 .or. lnocollapse_zdir_onecell) &
              fp(1:npar_loc,izp) = xyz0_par(3)+fp(1:npar_loc,izp)*Lxyz_par(3)
!
!  Special thing where Lxyz_par(2) has been replaced by Lxyz_loc(2)
!
        case ('random_proc')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k = 1,npar_loc
            if (nxgrid /= 1 .or. lnocollapse_xdir_onecell) then
              call random_number_wrapper(r)
              fp(k,ixp) = r
            endif
            if (nygrid /= 1 .or. lnocollapse_ydir_onecell) then
              call random_number_wrapper(r)
              fp(k,iyp) = r
            endif
            if (nzgrid /= 1 .or. lnocollapse_zdir_onecell) then
              call random_number_wrapper(r)
              fp(k,izp) = r
            endif
          enddo
          if (nxgrid /= 1 .or. lnocollapse_xdir_onecell) &
              fp(1:npar_loc,ixp) = xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          if (nygrid /= 1 .or. lnocollapse_ydir_onecell) &
              fp(1:npar_loc,iyp) = xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
          if (nzgrid /= 1 .or. lnocollapse_zdir_onecell) &
              fp(1:npar_loc,izp) = xyz0_par(3)+fp(1:npar_loc,izp)*Lxyz_par(3)
!
        case ('random-shift')
          if (lroot) print*, 'init_particles: Randomly distribute particle positions, then shift them'
          k2_xxp = kx_xxp**2 + ky_xxp**2 + kz_xxp**2
          do k = 1, npar_loc
            if (nxgrid /= 1) then
              call random_number_wrapper(r)
              fp(k, ixp) = xyz0_par(1) + r * Lxyz_par(1) - &
                  kx_xxp/k2_xxp * amplxxp * sin(kx_xxp*fp(k, ixp) + ky_xxp*fp(k, iyp) + kz_xxp*fp(k, izp))
            endif
            if (nygrid /= 1) then
              call random_number_wrapper(r)
              fp(k, iyp) = xyz0_par(2) + r * Lxyz_par(2) - &
                  ky_xxp/k2_xxp * amplxxp * sin(kx_xxp*fp(k, ixp) + ky_xxp*fp(k, iyp) + kz_xxp*fp(k, izp))
            endif
            if (nzgrid /= 1) then
              call random_number_wrapper(r)
              fp(k, izp) = xyz0_par(3) + r * Lxyz_par(3) - &
                  kz_xxp/k2_xxp * amplxxp * sin(kx_xxp*fp(k, ixp) + ky_xxp*fp(k, iyp) + kz_xxp*fp(k, izp))
            endif
          enddo
!
        case ('random-circle')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k = 1,npar_loc
            call random_number_wrapper(r)
            if (zp0 > yp0) then
              fp(k,ixp) = xp0*cos((zp0-yp0)*r+yp0)
              fp(k,iyp) = xp0*sin((zp0-yp0)*r+yp0)
            else
              fp(k,ixp) = xp0*cos(2*pi*r)
              fp(k,iyp) = xp0*sin(2*pi*r)
            endif
          enddo
!
        case ('random-sphere')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'in a sphere around (0,0,0) with radius=',rad_sphere
          if (rad_sphere == 0) &
            call fatal_error('init_particles','random-sphere radius needs to be larger than zero')

          if (-rad_sphere+pos_sphere(1) < xyz0(1) .or. &
              +rad_sphere+pos_sphere(1) > xyz1(1) .or. &
              -rad_sphere+pos_sphere(2) < xyz0(2) .or. &
              +rad_sphere+pos_sphere(2) > xyz1(2) .or. &
              -rad_sphere+pos_sphere(3) < xyz0(3) .or. &
              +rad_sphere+pos_sphere(3) > xyz1(3)) then
            call fatal_error('init_particles','random-sphere needs to fit in the box')
          endif
          if (lcartesian_coords) then
            do k = 1,npar_loc
              rp2 = 2.*rad_sphere**2
              do while (rp2 > rad_sphere**2)
                call random_number_wrapper(r)
                fp(k,ixp) = (r-0.5)*2.*rad_sphere
                call random_number_wrapper(r)
                fp(k,iyp) = (r-0.5)*2.*rad_sphere
                call random_number_wrapper(r)
                fp(k,izp) = (r-0.5)*2.*rad_sphere
                rp2 = fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2
              enddo
              fp(k,ixp) = fp(k,ixp)+pos_sphere(1)
              fp(k,iyp) = fp(k,iyp)+pos_sphere(2)
              fp(k,izp) = fp(k,izp)+pos_sphere(3)
            enddo
          else
            call not_implemented('init_particles','random-sphere for non-Cartesian coordinates')
          endif
!
        case ('random-ellipsoid')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'in an ellipsoid around ', pos_ellipsoid, ' with ' // &
              'semi-principal axes a,b,c =',a_ellipsoid,b_ellipsoid,c_ellipsoid
          if ((a_ellipsoid == 0) .or. (b_ellipsoid == 0) .or. (c_ellipsoid == 0)) then
            call fatal_error('init_particles','random-ellipsoid '// &
                             'all semi-principal axes need to be larger than zero')
          endif
          if (-a_ellipsoid+pos_ellipsoid(1) < xyz0(1) .or. &
              +a_ellipsoid+pos_ellipsoid(1) > xyz1(1) .or. &
              -b_ellipsoid+pos_ellipsoid(2) < xyz0(2) .or. &
              +b_ellipsoid+pos_ellipsoid(2) > xyz1(2) .or. &
              -c_ellipsoid+pos_ellipsoid(3) < xyz0(3) .or. &
              +c_ellipsoid+pos_ellipsoid(3) > xyz1(3)) then
            call fatal_error('init_particles','random-ellipsoid ellipsoid needs to fit in the box')
          endif
          if (lcartesian_coords) then
            a_ell2 = a_ellipsoid**2
            b_ell2 = b_ellipsoid**2
            c_ell2 = c_ellipsoid**2
            do k = 1,npar_loc
              rp2 = 2.
              do while (rp2 > 1.)
                call random_number_wrapper(r)
                fp(k,ixp) = (r-0.5)*2.*a_ellipsoid
                call random_number_wrapper(r)
                fp(k,iyp) = (r-0.5)*2.*b_ellipsoid
                call random_number_wrapper(r)
                fp(k,izp) = (r-0.5)*2.*c_ellipsoid
                rp2 = fp(k,ixp)**2/a_ell2+fp(k,iyp)**2/b_ell2+fp(k,izp)**2/c_ell2
              enddo
              fp(k,ixp) = fp(k,ixp)+pos_ellipsoid(1)
              fp(k,iyp) = fp(k,iyp)+pos_ellipsoid(2)
              fp(k,izp) = fp(k,izp)+pos_ellipsoid(3)
            enddo
          else
            call not_implemented('init_particles','random-ellipsoid for non-Cartesian coordinates')
          endif
!
        case ('random-line-x')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k = 1,npar_loc
            if (nxgrid /= 1) then
              call random_number_wrapper(r)
              fp(k,ixp) = r
            endif
          enddo
          if (nxgrid /= 1) &
              fp(1:npar_loc,ixp) = xyz0_par(1)+fp(1:npar_loc,ixp)*Lxyz_par(1)
          fp(1:npar_loc,iyp) = yp0
          fp(1:npar_loc,izp) = zp0
!
        case ('random-line-y')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k = 1,npar_loc
            if (nygrid /= 1) then
              call random_number_wrapper(r)
              fp(k,iyp) = r
            endif
          enddo
          if (nygrid /= 1) &
              fp(1:npar_loc,iyp) = xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          fp(1:npar_loc,ixp) = xp0
          fp(1:npar_loc,izp) = zp0
!
        case ('random-hole')
          if (lroot) print*, 'init_particles: Random particle positions with inner hole'
          do k = 1,npar_loc
            rp2 = -1.0
            do while (rp2 < rp_int**2)
              rp2 = 0.0
              if (nxgrid /= 1) then
                call random_number_wrapper(r)
                fp(k,ixp) = xyz0(1)+r*Lxyz(1)
                rp2 = rp2+fp(k,ixp)**2
              endif
              if (nygrid /= 1) then
                call random_number_wrapper(r)
                fp(k,iyp) = xyz0(2)+r*Lxyz(2)
                rp2 = rp2+fp(k,iyp)**2
              endif
              if (nzgrid /= 1) then
                call random_number_wrapper(r)
                fp(k,izp) = xyz0(3)+r*Lxyz(3)
                rp2 = rp2+fp(k,izp)**2
              endif
            enddo
          enddo
!
        case ('random-box')
          if (lroot) print*, 'init_particles: Random particle positions within a box'
          do k = 1,npar_loc
            if (nxgrid /= 1) then
              call random_number_wrapper(r)
              fp(k,ixp) = r
            endif
            if (nygrid /= 1) then
              call random_number_wrapper(r)
              fp(k,iyp) = r
            endif
            if (nzgrid /= 1) then
              call random_number_wrapper(r)
              fp(k,izp) = r
            endif
            if (lcylindrical_coords) then
              xx0 = xp0+fp(k,ixp)*Lx0
              yy0 = yp0+fp(k,iyp)*Ly0
              r2 = xx0**2+yy0**2
              if (nxgrid /= 1) fp(k,ixp) = sqrt(r2)
              if (nygrid /= 1) fp(k,iyp) = atan(yy0/xx0)+pi*(xx0/abs(xx0)-1)*0.5
              if (nzgrid /= 1) fp(k,izp) = zp0+fp(k,izp)*Lz0
            else
              if (nxgrid /= 1) fp(k,ixp) = xp0+fp(k,ixp)*Lx0
              if (nygrid /= 1) fp(k,iyp) = yp0+fp(k,iyp)*Ly0
              if (nzgrid /= 1) fp(k,izp) = zp0+fp(k,izp)*Lz0
            endif
          enddo
!
        case ('random-cylindrical','random-cyl')
!
          if (lroot) print*, 'init_particles: Random particle '// &
              'cylindrical positions with power-law = ', dustdensity_powerlaw
!
          do k = 1,npar_loc
!
! Start the particles obeying a power law
!
            if (lcylindrical_coords .or. lcartesian_coords) then
              tmp = 2-dustdensity_powerlaw
            elseif (lspherical_coords) then
              tmp = 3-dustdensity_powerlaw
            else
              call fatal_error("init_particles","unsupported coordinate system")
            endif
!
            if (lcartesian_coords) then
              if (nprocx==1) then
                rpar_int=rp_int
                rpar_ext=rp_ext
              else
                call not_implemented("init_particles", &
                                     "random-cyl for nprocx/=1 in Cartesian. Parallelize in y or z")
              endif
            else
              rpar_int = xyz0_loc(1)
              rpar_ext = xyz1_loc(1)
            endif
            call random_number_wrapper(rad_scl)
            rad_scl = rpar_int**tmp + rad_scl*(rpar_ext**tmp-rpar_int**tmp)
            rad = rad_scl**(1./tmp)
!
! Random in azimuth
!
            if (lcartesian_coords) then
              call random_number_wrapper(phi)
              phi = 2*pi*phi
              if (nxgrid /= 1) fp(k,ixp) = rad*cos(phi)
              if (nygrid /= 1) fp(k,iyp) = rad*sin(phi)
              if (nzgrid /= 1) then
                call random_number_wrapper(r)
                fp(k,izp) = xyz0_par(3)+r*Lxyz_par(3)
              endif
            elseif (lcylindrical_coords) then
              if (nxgrid /= 1) fp(k,ixp) = rad
              if (nygrid /= 1) then
                call random_number_wrapper(phi)
                fp(k,iyp) = xyz0_par(2)+phi*Lxyz_par(2)
              endif
              if (nzgrid /= 1) then
                call random_number_wrapper(r)
                fp(k,izp) = xyz0_par(3)+r*Lxyz_par(3)
              endif
            elseif (lspherical_coords) then
              if (nxgrid /= 1) fp(k,ixp) = rad
              if (nygrid /= 1) then
                call random_number_wrapper(tht)
                fp(k,iyp) = xyz0_par(2)+tht*Lxyz_par(2)
              endif
              if (nzgrid /= 1) then
                call random_number_wrapper(phi)
                fp(k,izp) = xyz0_par(3)+phi*Lxyz_par(3)
              endif
            endif
!
          enddo
!
        case ('np-constant')
          if (lroot) print*, 'init_particles: Constant number density'
          k = 1
          do while (.not. (k > npar_loc))
            do l = l1,l2
              do m = m1,m2
                do n = n1,n2
                  if (nxgrid /= 1) then
                    call random_number_wrapper(px)
                    fp(k,ixp) = x(l)+(px-0.5)*dx
                  endif
                  if (nygrid /= 1) then
                    call random_number_wrapper(py)
                    fp(k,iyp) = y(m)+(py-0.5)*dy
                  endif
                  if (nzgrid /= 1) then
                    call random_number_wrapper(pz)
                    fp(k,izp) = z(n)+(pz-0.5)*dz
                  endif
                  k = k+1
                  if (k > npar_loc) exit
                enddo
              enddo
            enddo
          enddo
!
        case ('equidistant')
          if (lroot) print*, 'init_particles: Particles placed equidistantly'
          dim1 = 1.0/max(dimensionality,1)
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
          npar_loc_x = 1
          npar_loc_y = 1
          npar_loc_z = 1
!
          if (dimensionality == 3) then
!  3-D
            npar_loc_x = nint((npar_loc*Lxyz_loc(1)**2/(Lxyz_loc(2)*Lxyz_loc(3)))**dim1)
            npar_loc_y = nint((npar_loc*Lxyz_loc(2)**2/(Lxyz_loc(1)*Lxyz_loc(3)))**dim1)
            npar_loc_z = nint((npar_loc*Lxyz_loc(3)**2/(Lxyz_loc(1)*Lxyz_loc(2)))**dim1)
          elseif (dimensionality == 2) then
!  2-D
            if (nxgrid == 1) then
              npar_loc_x = 1
              npar_loc_y = nint((npar_loc*Lxyz_loc(2)/Lxyz_loc(3))**dim1)
              npar_loc_z = nint((npar_loc*Lxyz_loc(3)/Lxyz_loc(2))**dim1)
            elseif (nygrid == 1) then
              npar_loc_x = nint((npar_loc*Lxyz_loc(1)/Lxyz_loc(3))**dim1)
              npar_loc_y = 1
              npar_loc_z = nint((npar_loc*Lxyz_loc(3)/Lxyz_loc(1))**dim1)
            elseif (nzgrid == 1) then
              npar_loc_x = nint((npar_loc*Lxyz_loc(1)/Lxyz_loc(2))**dim1)
              npar_loc_y = nint((npar_loc*Lxyz_loc(2)/Lxyz_loc(1))**dim1)
              npar_loc_z = 1
            endif
          elseif (dimensionality == 1) then
!  1-D
            if (nxgrid /= 1) then
              npar_loc_x = npar_loc
              npar_loc_y = 1
              npar_loc_z = 1
            elseif (nygrid /= 1) then
              npar_loc_x = 1
              npar_loc_y = npar_loc
              npar_loc_z = 1
            elseif (nzgrid /= 1) then
              npar_loc_x = 1
              npar_loc_y = 1
              npar_loc_z = npar_loc
            endif
          endif
!  Distance between particles.
          dx_par = Lxyz_loc(1)/npar_loc_x
          dy_par = Lxyz_loc(2)/npar_loc_y
          dz_par = Lxyz_loc(3)/npar_loc_z
!  Place first particle.
          fp(1,ixp) = x(l1)
          fp(1,iyp) = y(m1)
          fp(1,izp) = z(n1)
          if (nxgrid /= 1) fp(1,ixp) = xyz0_loc(1)+dx_par/2
          if (nygrid /= 1) fp(1,iyp) = xyz0_loc(2)+dy_par/2
          if (nzgrid /= 1) fp(1,izp) = xyz0_loc(3)+dz_par/2
!  Place all other particles iteratively.
          if (dimensionality == 3) then
!  3-D
            do k = 2,npar_loc
              fp(k,ixp) = fp(1,ixp) + mod(k-1, npar_loc_x) * dx_par
              fp(k,iyp) = fp(1,iyp) + mod((k-1)/npar_loc_x, npar_loc_y) * dy_par
              fp(k,izp) = fp(1,izp) + (k-1) / (npar_loc_x*npar_loc_y) * dz_par
            enddo
!
          elseif (dimensionality == 2) then
!  2-D
            if (nxgrid == 1) then
              do k = 2,npar_loc
                fp(k,ixp) = fp(1,ixp)
                fp(k,iyp) = fp(1,iyp) + mod(k-1, npar_loc_y) * dy_par
                fp(k,izp) = fp(1,izp) + (k-1) / npar_loc_y * dz_par
              enddo
            elseif (nygrid == 1) then
              do k = 2,npar_loc
                fp(k,ixp) = fp(1,ixp) + mod(k-1, npar_loc_x) * dx_par
                fp(k,iyp) = fp(1,iyp)
                fp(k,izp) = fp(1,izp) + (k-1) / npar_loc_x * dz_par
              enddo
            elseif (nzgrid == 1) then
              do k = 2,npar_loc
                fp(k,ixp) = fp(1,ixp) + mod(k-1, npar_loc_x) * dx_par
                fp(k,iyp) = fp(1,iyp) + (k-1) / npar_loc_x * dy_par
                fp(k,izp) = fp(1,izp)
              enddo
            endif
          elseif (dimensionality == 1) then
!  1-D
            if (nxgrid /= 1) then
              do k = 2,npar_loc
                fp(k,ixp) = fp(1,ixp)+(k-1)*dx_par
                fp(k,iyp) = fp(1,iyp)
                fp(k,izp) = fp(1,izp)
              enddo
            elseif (nygrid /= 1) then
              do k = 2,npar_loc
                fp(k,ixp) = fp(1,ixp)
                fp(k,iyp) = fp(1,iyp)+(k-1)*dy_par
                fp(k,izp) = fp(1,izp)
              enddo
            elseif (nzgrid /= 1) then
              do k = 2,npar_loc
                fp(k,ixp) = fp(1,ixp)
                fp(k,iyp) = fp(1,iyp)
                fp(k,izp) = fp(1,izp)+(k-1)*dz_par
              enddo
            endif
          else
!  0-D
            fp(2:npar_loc,ixp) = fp(1,ixp)
            fp(2:npar_loc,iyp) = fp(1,iyp)
            fp(2:npar_loc,izp) = fp(1,izp)
          endif
          lequidistant = .true.
!
!  Shift particle locations slightly so that a mode appears.
!
        case ('shift')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) &
            call fatal_error('init_particles','must place particles equidistantly before shifting')
          k2_xxp = kx_xxp**2+ky_xxp**2+kz_xxp**2
          if (k2_xxp == 0.0) call fatal_error('init_particles','kx_xxp=ky_xxp=kz_xxp=0.0 not allowed')
          do k = 1,npar_loc
            fp(k,ixp) = fp(k,ixp) - kx_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,iyp) = fp(k,iyp) - ky_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,izp) = fp(k,izp) - kz_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
          enddo
!
!  Shift to egg crate mode 2d, cos(x)cos(z)
!
        case ('cosxcosz')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) &
            call fatal_error('init_particles','must place particles equidistantly before shifting')
          k2_xxp = kx_xxp**2+kz_xxp**2
          if (k2_xxp == 0.0) call fatal_error('init_particles','kx_xxp=ky_xxp=kz_xxp=0.0 not allowed')
          do k = 1,npar_loc
            fp(k,ixp) = fp(k,ixp) - kx_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp))*cos(kz_xxp*fp(k,izp))
            fp(k,izp) = fp(k,izp) - kz_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp))*cos(kz_xxp*fp(k,izp))
          enddo
!
!  Shift to egg crate mode 2d, sin(x)sin(z)
!
        case ('sinxsinz')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) &
            call fatal_error('init_particles','must place particles equidistantly before shifting')
          k2_xxp = kx_xxp**2+kz_xxp**2
          if (k2_xxp == 0.0) call fatal_error('init_particles','kx_xxp=ky_xxp=kz_xxp=0.0 not allowed')
          do k = 1,npar_loc
            fp(k,ixp) = fp(k,ixp) + kx_xxp/k2_xxp*amplxxp* &
                cos(kx_xxp*fp(k,ixp))*sin(kz_xxp*fp(k,izp))
            fp(k,izp) = fp(k,izp) + kz_xxp/k2_xxp*amplxxp* &
                cos(kx_xxp*fp(k,ixp))*sin(kz_xxp*fp(k,izp))
          enddo
!
        case ('gaussian-z')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k = 1,npar_loc
            do while (.true.)
              if (nxgrid /= 1) then
                call random_number_wrapper(r)
                fp(k,ixp) = r
                fp(k,ixp) = xyz0_par(1)+fp(k,ixp)*Lxyz_par(1)
              endif
              if (nygrid /= 1) then
                call random_number_wrapper(r)
                fp(k,iyp) = r
                fp(k,iyp) = xyz0_par(2)+fp(k,iyp)*Lxyz_par(2)
              endif
              call random_number_wrapper(r)
              call random_number_wrapper(p)
              if (nprocz == 2) then
                if (lfirst_proc_z) fp(k,izp) = -abs(zp0*(fp(k,ixp)/r0gaussz)**qgaussz*sqrt(-2*alog(r))*cos(2*pi*p))
                if (llast_proc_z) fp(k,izp) = abs(zp0*(fp(k,ixp)/r0gaussz)**qgaussz*sqrt(-2*alog(r))*cos(2*pi*p))
              else
                fp(k,izp) = zp0*(fp(k,ixp)/r0gaussz)**qgaussz*sqrt(-2*alog(r))*cos(2*pi*p)
              endif
              if ((fp(k,izp) >= xyz0(3)) .and. (fp(k,izp) <= xyz1(3))) exit
            enddo
          enddo
!
        case ('gaussian-x')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k = 1,npar_loc
            do while (.true.)
              if (nygrid /= 1) then
                call random_number_wrapper(r)
                fp(k,iyp) = r
              endif
              if (nzgrid /= 1) then
                call random_number_wrapper(r)
                fp(k,izp) = r
              endif
              call random_number_wrapper(r)
              call random_number_wrapper(p)
              fp(k,ixp) = xp0*sqrt(-2*alog(r))*cos(2*pi*p)
              if ((fp(k,ixp) >= xyz0(1)).and.(fp(k,ixp) <= xyz1(1))) exit
            enddo
          enddo
          if (nygrid /= 1) fp(1:npar_loc,iyp) = xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          if (nzgrid /= 1) fp(1:npar_loc,izp) = xyz0_par(3)+fp(1:npar_loc,izp)*Lxyz_par(3)
!
        case ('gaussian-x+random')
          if (lroot) print*, 'init_particles: Gaussian particle positions+random'
          do k = 1,npar_loc
            if (nygrid /= 1) then
              call random_number_wrapper(r)
              fp(k,iyp) = r
            endif
            if (nzgrid /= 1) then
              call random_number_wrapper(r)
              fp(k,izp) = r
            endif
          enddo
          do k = 1,nint(npar_loc*gaussian_x_frac)
            do while (.true.)
              call random_number_wrapper(r)
              call random_number_wrapper(p)
              fp(k,ixp) = xp0*sqrt(-2*alog(r))*cos(2*pi*p)
              if ((fp(k,ixp) >= xyz0(1)).and.(fp(k,ixp) <= xyz1(1))) exit
            enddo
          enddo
          do k = nint(npar_loc*gaussian_x_frac)+1, npar_loc
            call random_number_wrapper(r)
            fp(k,ixp) = xyz0_par(1)+r*Lxyz_par(1)
          enddo
          if (nygrid /= 1) fp(1:npar_loc,iyp) = xyz0_par(2)+fp(1:npar_loc,iyp)*Lxyz_par(2)
          if (nzgrid /= 1) fp(1:npar_loc,izp) = xyz0_par(3)+fp(1:npar_loc,izp)*Lxyz_par(3)
!
        case ('gaussian-z-pure')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k = 1,npar_loc
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            if (nprocz == 2) then
              if (lfirst_proc_z) fp(k,izp) = -abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
              if (llast_proc_z) fp(k,izp) = abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
            else
              fp(k,izp) = zp0*sqrt(-2*alog(r))*cos(2*pi*p)   ! generates a random gaussion number*zp0
            endif
          enddo
!
        case ('gaussian-r')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k = 1,npar_loc
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            call random_number_wrapper(q)
            fp(k,ixp) = xp0*sqrt(-2*alog(r))*cos(2*pi*p)*cos(2*pi*q)
            fp(k,iyp) = yp0*sqrt(-2*alog(r))*cos(2*pi*p)*sin(2*pi*q)
          enddo
!
        case ('hole')
!
          call map_nearest_grid(fp,ineargrid)
          call map_xxp_grid(f,fp,ineargrid)
          call sort_particles_imn(fp,ineargrid,ipar)
          do k = k1_imn(imn_array(m_hole+m1-1,n_hole+n1-1)), &
                k2_imn(imn_array(m_hole+m1-1,n_hole+n1-1))
            if (ineargrid(k,1) == l_hole+l1-1) then
              print*, k
              if (nxgrid /= 0) fp(k,ixp) = fp(k,ixp)-dx
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
        case ('birthring')
          if (birthring_width > tini) then
            if (lgaussian_birthring) then
              do k = 1,npar_loc
                call normal_deviate(rr_tmp(k))
              enddo
            else
              call random_number_wrapper(rr_tmp(1:npar_loc))
            endif
            rr_tmp(1:npar_loc) = birthring_r+rr_tmp(1:npar_loc)*birthring_width
          else
            rr_tmp(1:npar_loc) = birthring_r
          endif
          call random_number_wrapper(az_tmp(1:npar_loc))
          az_tmp(1:npar_loc) = -pi+az_tmp(1:npar_loc)*2.0*pi
          if (lcartesian_coords) then
            fp(1:npar_loc,ixp) = rr_tmp(1:npar_loc)*cos(az_tmp(1:npar_loc))
            fp(1:npar_loc,iyp) = rr_tmp(1:npar_loc)*sin(az_tmp(1:npar_loc))
            fp(1:npar_loc,izp) = 0.0
          else
            fp(1:npar_loc,ixp) = rr_tmp(1:npar_loc)
            if (lcylindrical_coords) then
              fp(1:npar_loc,iyp) = az_tmp(1:npar_loc)
              fp(1:npar_loc,izp) = 0.0
            elseif (lspherical_coords) then
              fp(1:npar_loc,iyp) = pi/2.0
              fp(1:npar_loc,izp) = az_tmp(1:npar_loc)
            endif
          endif
          if (nzgrid /= 0) call warning('init_particles',"birthring only implemented for 2D(xy)")
!
        case default
          call fatal_error('init_particles','no such initxxp: "'//trim(initxxp(j))//'"')
!
        endselect
!
      enddo ! do j=1,ninit
!
!  Interface for user's own initial condition for position
!
      if (linitial_condition) call initial_condition_xxp(f,fp)
!
!  Particles are by default not allowed to be present in non-existing
!  dimensions. This could give huge problems with interpolation later.
!
      if (nxgrid == 1 .and. .not. lnocollapse_xdir_onecell) &
          fp(1:npar_loc,ixp) = x(nghost+1)
      if (nygrid == 1 .and. .not. lnocollapse_ydir_onecell) &
          fp(1:npar_loc,iyp) = y(nghost+1)
      if (nzgrid == 1 .and. .not. lnocollapse_zdir_onecell) &
          fp(1:npar_loc,izp) = z(nghost+1)
!
      if (init_repeat /= 0) call repeated_init(fp,init_repeat)
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
      do j = 1,ninit
!
        select case (initvvp(j))
!
        case ('nothing')
          if (lroot .and. j == 1) print*, 'init_particles: No particle velocity set'
!
        case ('zero')
          if (lroot) print*, 'init_particles: Zero particle velocity'
          fp(1:npar_loc,ivpx:ivpz) = 0.0
!
        case ('zero-shear')
          if (lroot) print*, 'init_particles: Zero particle velocity'
          fp(1:npar_loc,ivpy) = -Sshear*fp(1:npar_loc,ixp)
          fp(1:npar_loc,ivpx) = 0.0
          fp(1:npar_loc,ivpz) = 0.0
!
        case ('constant')
          if (lroot) print*, 'init_particles: Constant particle velocity'
          if (lroot) print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          if (lcylindrical_coords) then
            fp(1:npar_loc,ivpx) = vpx0*cos(fp(k,iyp))+vpy0*sin(fp(k,iyp))
            fp(1:npar_loc,ivpy) = vpy0*cos(fp(k,iyp))-vpx0*sin(fp(k,iyp))
            fp(1:npar_loc,ivpz) = vpz0
          else
            fp(1:npar_loc,ivpx) = vpx0
            fp(1:npar_loc,ivpy) = vpy0
            fp(1:npar_loc,ivpz) = vpz0
          endif
!
        case ('constant-1')
          if (lroot) print*, 'init_particles: Particle 1 velocity; vx,vy,vz=', vpx1, vpy1, vpz1
          do k = 1,npar_loc
            if (ipar(k) == 1) then
              fp(k,ivpx) = vpx1
              fp(k,ivpy) = vpy1
              fp(k,ivpz) = vpz1
            endif
          enddo
!
        case ('constant-2')
          if (lroot) print*, 'init_particles: Particle 2 velocity; vx,vy,vz=', vpx2, vpy2, vpz2
          do k = 1,npar_loc
            if (ipar(k) == 2) then
              fp(k,ivpx) = vpx2
              fp(k,ivpy) = vpy2
              fp(k,ivpz) = vpz2
            endif
          enddo
!
        case ('constant-3')
          if (lroot) print*, 'init_particles: Particle 3 velocity; vx,vy,vz=', vpx3, vpy3, vpz3
          do k = 1,npar_loc
            if (ipar(k) == 3) then
              fp(k,ivpx) = vpx3
              fp(k,ivpy) = vpy3
              fp(k,ivpz) = vpz3
            endif
          enddo
!
        case ('sinwave-phase')
          if (lroot) print*, 'init_particles: sinwave-phase; vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k = 1,npar_loc
            fp(k,ivpx) = fp(k,ivpx)+vpx0*sin(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy) = fp(k,ivpy)+vpy0*sin(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz) = fp(k,ivpz)+vpz0*sin(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
!
        case ('coswave-phase')
          if (lroot) print*, 'init_particles: coswave-phase; vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k = 1,npar_loc
            fp(k,ivpx) = fp(k,ivpx)+vpx0*cos(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy) = fp(k,ivpy)+vpy0*cos(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz) = fp(k,ivpz)+vpz0*cos(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle velocities; delta_vp0=', delta_vp0
          do k = 1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpx) = fp(k,ivpx) + delta_vp0*(2*r-1)
            call random_number_wrapper(r)
            fp(k,ivpy) = fp(k,ivpy) + delta_vp0*(2*r-1)
            call random_number_wrapper(r)
            fp(k,ivpz) = fp(k,ivpz) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-x')
          if (lroot) print*, 'init_particles: Random particle x-velocity; delta_vp0=', delta_vp0
          do k = 1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpx) = fp(k,ivpx) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-y')
          if (lroot) print*, 'init_particles: Random particle y-velocity; delta_vp0=', delta_vp0
          do k = 1,npar_loc
            call random_number_wrapper(r)
            fp(k,ivpy) = fp(k,ivpy) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-z')
          if (lroot) print*, 'init_particles: Random particle z-velocity; delta_vp0=', delta_vp0
          do k = 1,npar_loc
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
          fp(1:npar_loc,ivpx) = fp(1:npar_loc,ivpx)-vpx_sum/npar
          fp(1:npar_loc,ivpy) = fp(1:npar_loc,ivpy)-vpy_sum/npar
          fp(1:npar_loc,ivpz) = fp(1:npar_loc,ivpz)-vpz_sum/npar
!
        case ('follow-gas')
          if (lroot) &
              print*, 'init_particles: Particle velocity equal to gas velocity'
          do k = 1,npar_loc
            call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:),0,0)
            fp(k,ivpx:ivpz) = uup
          enddo
!
        case ('jeans-wave-dustpar-x')
! assumes rhs_poisson_const=1 !
          do k = 1,npar_loc
            fp(k,ivpx) = fp(k,ivpx) - amplxxp*(sqrt(1+4*1.0*1.0*tausp**2)-1)/ &
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
          cs = sqrt(cs20)
!
          if (ldragforce_gas_par) then
            eps = eps_dtog ! will be recomputed if ldragforce_equi_global_eps = .true.
          else
            eps = 0.0
          endif
!
          if (lroot) print*, 'init_particles: average dust-to-gas ratio=', eps
!
          uudrag: if (lhydro .and. (.not. lread_oldsnap .or. lread_oldsnap_nohydro)) then
!  Set gas velocity field.
            do l = l1,l2
              do m = m1,m2
                do n = n1,n2
!  Take either global or local dust-to-gas ratio.
                  if (ldragforce_gas_par .and. .not. ldragforce_equi_global_eps) eps = f(l,m,n,irhop) / get_gas_density(f,l,m,n)
!
                  f(l,m,n,iux) = f(l,m,n,iux) - beta_glnrho_global(1)*eps*Omega*tausp/ &
                                                ((1.0+eps)**2+(Omega*tausp)**2)*cs
                  f(l,m,n,iuy) = f(l,m,n,iuy) + beta_glnrho_global(1)*(1+eps+(Omega*tausp)**2)/ &
                                                (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
                enddo
              enddo
            enddo
          endif uudrag
!  Set particle velocity field.
          do k = 1,npar_loc
!  Take either global or local dust-to-gas ratio.
            if (.not. ldragforce_gas_par) then
              eps = 0.0
            else if (.not. ldragforce_equi_global_eps) then
              ix0 = ineargrid(k,1)
              iy0 = ineargrid(k,2)
              iz0 = ineargrid(k,3)
              eps = f(ix0,iy0,iz0,irhop) / get_gas_density(f,ix0,iy0,iz0)
            endif
!
            fp(k,ivpx) = fp(k,ivpx) + beta_glnrho_global(1)*Omega*tausp/ &
                         ((1.0+eps)**2+(Omega*tausp)**2)*cs
            fp(k,ivpy) = fp(k,ivpy) + beta_glnrho_global(1)*(1+eps)/ &
                         (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
          enddo
!
        case ('dragforce_equi_nohydro')
!
          do k = 1,npar_loc
            if (lparticles_radius) then
              tausp_par = rhopmat*fp(k,iap)/(sqrt(cs20)*rho0)
            else
              tausp_par = tausp
            endif
            cs = sqrt(cs20)
            if (Deltauy_gas_friction /= 0.0) then
              fp(k,ivpx) = fp(k,ivpx) - 2*Deltauy_gas_friction/(1.0/(Omega*tausp_par)+Omega*tausp_par)
              fp(k,ivpy) = fp(k,ivpy) - Deltauy_gas_friction/(1.0+(Omega*tausp_par)**2)
            else
              fp(k,ivpx) = fp(k,ivpx) + beta_glnrho_global(1)*cs/(1.0/(Omega*tausp_par)+Omega*tausp_par)
              fp(k,ivpy) = fp(k,ivpy) + beta_glnrho_global(1)/2*cs/(1.0+(Omega*tausp_par)**2)
            endif
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
          cs = sqrt(cs20)
          do k = 1,npar_loc
            fp(k,ivpx) = fp(k,ivpx) + 1/(Omega*tausp+1/(Omega*tausp))*beta_dPdr_dust*cs
            fp(k,ivpy) = fp(k,ivpy) - 1/(1.0+1/(Omega*tausp)**2)*beta_dPdr_dust/2*cs
          enddo
!
        case ('Keplerian','keplerian')
!
!  Keplerian velocity based on gravr.
!
          if (lroot) then
            print*, 'init_particles: Keplerian velocity'
            if (lshear) call not_implemented("init_particles", &
                "Keplerian initial condition for shearing boxes")
          endif
!
          do k = 1,npar_loc
            if (lcartesian_coords) then
              rad = sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2)
              OO = sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) = -OO*fp(k,iyp)
              fp(k,ivpy) =  OO*fp(k,ixp)
              fp(k,ivpz) =  0.0
            elseif (lcylindrical_coords) then
              rad = fp(k,ixp)
              OO = sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) =  0.0
              fp(k,ivpy) =  (OO-Omega_corot)*rad
              fp(k,ivpz) =  0.0
            elseif (lspherical_coords) then
              rad = fp(k,ixp)*sin(fp(k,iyp))
              OO = sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) =  0.0
              fp(k,ivpy) =  0.0
              fp(k,ivpz) =  (OO-Omega_corot)*rad
            endif
          enddo
!
!  Explosion.
!
        case ('explosion')
          do k = 1,npar_loc
            rad = sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2)
            fp(k,ivpx) = delta_vp0*fp(k,ixp)/rp_ext
            fp(k,ivpy) = delta_vp0*fp(k,iyp)/rp_ext
            fp(k,ivpz) = delta_vp0*fp(k,izp)/rp_ext
          enddo
!
        case default
          call fatal_error('init_particles','no such initvvp: "'//trim(initvvp(j))//'"')
        endselect
!
      enddo ! do j=1,ninit
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_vvp(f, fp, ineargrid)
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
!  If we are solving for caustics, then we call their initial condition here:
!
      if (lparticles_caustics) call init_particles_caustics(f,fp,ineargrid)
!
!  If we are solving for tetrad, then we call their initial condition here:
!
      if (lparticles_tetrad) call init_particles_tetrad(f,fp,ineargrid)
!
!  Set the initial auxiliary array for the passive scalar to zero
!
      if (ldiffuse_passive) f(:,:,:,idlncc) = 0.0
      if (ldiffuse_dragf) f(:,:,:,idfx:idfz) = 0.0
!
    endsubroutine init_particles
!***********************************************************************
    subroutine insert_lost_particles(f,fp,ineargrid)
!
!  14-oct-12/dhruba: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
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
      use General, only: random_number_wrapper, normal_deviate
      use Particles_diagnos_state, only: insert_particles_diagnos_state
      use Mpicomm, only: mpireduce_sum_int
      use Particles_number, only: set_particle_number
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
      real, dimension(mpar_loc) :: rr_tmp, az_tmp, OO_tmp
      real, dimension(3) :: uup
!
      logical, save :: linsertmore=.true.
      real :: xx0, yy0, r2, r, tmp
      integer :: j, k, n_insert, npar_loc_old, iii
      integer (kind=8) :: particles_insert_rate_tmp
      real, pointer :: gravr
!
! Insertion of particles is stopped when maximum number of particles is reached,
! unless linsert_as_many_as_possible is set.
! Maximum numer of particles allowed in system is defined by max_particles,
! initialized to npar. Note that this may cause errors at a processor further
! downstream, if particles accumulate and mpar_loc is too small.
! Since root inserts all new particles, make sure
! npar_loc + n_insert < mpar_loc
! so that a processor can not exceed its maximum number of particles.
!
      if (t >= tstart_insert_particles+particles_insert_ramp_time) then
        particles_insert_rate_tmp = particles_insert_rate
      else
        particles_insert_rate_tmp = int(real(particles_insert_rate) &
            *(t-tstart_insert_particles)/particles_insert_ramp_time)
      endif
      call mpireduce_sum_int(npar_loc,npar_total)
!
      if (lroot) then
        avg_n_insert = particles_insert_rate_tmp*dt
        n_insert = int(avg_n_insert + remaining_particles)
!        
! Remaining particles saved for subsequent timestep:
!
        remaining_particles = avg_n_insert + remaining_particles - n_insert
!
! Insert particles if maximum count is not reached
!
        if ((n_insert+npar_loc <= mpar_loc) .and. (npar_inserted_tot+n_insert <= max_particles) &
            .and. (t < max_particle_insert_time) .and. (t > tstart_insert_particles)) then
          linsertmore = .true.
        else
          linsertmore = .false.
        endif
!
! Continue inserting even if  max_particles is reached, if linsert_as_many_as_possible
! is set.
!
        if (linsert_as_many_as_possible) then
          n_insert = min((mpar_loc-npar_loc),n_insert)
          if ((n_insert+npar_loc <= mpar_loc)  &
              .and. (t < max_particle_insert_time) .and. (t > tstart_insert_particles)) then
            linsertmore = .true.
          endif
        endif
!
        if (linsertmore) then
!
! Actual (integer) number of particles to be inserted at this timestep:
!
          do iii = npar_loc+1,npar_loc+n_insert
            ipar(iii) = npar_inserted_tot+iii-npar_loc
          enddo
          npar_loc_old = npar_loc
          npar_loc = npar_loc + n_insert
!
! Update total number of inserted particles, npar_inserted_tot.
! Not the same as npar_total, which is the number of particles in the system,
! without counting removed particles
!
          npar_inserted_tot = n_insert + npar_inserted_tot          
!
! Insert particles in chosen position (as in init_particles).
!
          do j = 1,ninit
            select case (initxxp(j))
            case ('random-box')
!
              do k = npar_loc_old+1,npar_loc
                if (nxgrid /= 1) then
                  call random_number_wrapper(r)
                  fp(k,ixp) = r
                endif
                if (nygrid /= 1) then
                  call random_number_wrapper(r)
                  fp(k,iyp) = r
                endif
                if (nzgrid /= 1) then
                  call random_number_wrapper(r)
                  fp(k,izp) = r
                endif
                if (lcylindrical_coords) then
                  xx0 = xp0+fp(k,ixp)*Lx0
                  yy0 = yp0+fp(k,iyp)*Ly0
                  r2 = xx0**2+yy0**2
                  if (nxgrid /= 1) fp(k,ixp) = sqrt(r2)
                  if (nygrid /= 1) fp(k,iyp) = atan(yy0/xx0)+pi*(xx0/abs(xx0)-1)*0.5
                  if (nzgrid /= 1) fp(k,izp) = zp0+fp(k,izp)*Lz0
                else
                  if (nxgrid /= 1) fp(k,ixp) = xp0+fp(k,ixp)*Lx0
                  if (nygrid /= 1) fp(k,iyp) = yp0+fp(k,iyp)*Ly0
                  if (nzgrid /= 1) fp(k,izp) = zp0+fp(k,izp)*Lz0
                endif
              enddo
!
! Maybe random-cylindrical case should be combined with normal initxxp case
!
            case ('random-cylindrical','random-cyl')
              if (lcylindrical_coords .or. lcartesian_coords) then
                tmp = 2-dustdensity_powerlaw
              elseif (lspherical_coords) then
                tmp = 3-dustdensity_powerlaw
              else
                call fatal_error("init_particles","unsupported coordinate system")
              endif
!
              call random_number_wrapper(rr_tmp(npar_loc_old+1:npar_loc))
              rr_tmp(npar_loc_old+1:npar_loc) = rp_int**tmp + &
                  rr_tmp(npar_loc_old+1:npar_loc)*(rp_ext**tmp-rp_int**tmp)
              rr_tmp(npar_loc_old+1:npar_loc) = rr_tmp(npar_loc_old+1:npar_loc)**(1./tmp)

              if ((lcartesian_coords) .or. (lcylindrical_coords.and.nygrid /= 1) .or. (lspherical_coords.and.nzgrid /= 1)) &
                call random_number_wrapper(az_tmp(npar_loc_old+1:npar_loc))
              if ((lcartesian_coords) .or. (lcylindrical_coords.and.nzgrid /= 1) .or. (lspherical_coords.and.nygrid /= 1)) &
                call random_number_wrapper(fp(npar_loc_old+1:npar_loc,izp))
!
              if (lcartesian_coords) then
                fp(npar_loc_old+1:npar_loc,iyp) = -pi+az_tmp(npar_loc_old+1:npar_loc)*2.0*pi
                if (nxgrid /= 1) fp(npar_loc_old+1:npar_loc,ixp) = rr_tmp(npar_loc_old+1:npar_loc) &
                    *cos(az_tmp(npar_loc_old+1:npar_loc))
                if (nygrid /= 1) fp(npar_loc_old+1:npar_loc,iyp) = rr_tmp(npar_loc_old+1:npar_loc) &
                    *sin(az_tmp(npar_loc_old+1:npar_loc))
                if (nzgrid /= 1) fp(npar_loc_old+1:npar_loc,izp) = xyz0(3)+fp(npar_loc_old+1:npar_loc,izp)*Lxyz(3)
              elseif (lcylindrical_coords) then
                if (nxgrid /= 1) fp(npar_loc_old+1:npar_loc,ixp) = rr_tmp(npar_loc_old+1:npar_loc)
                if (nygrid /= 1) fp(npar_loc_old+1:npar_loc,iyp) = xyz0(2)+az_tmp(npar_loc_old+1:npar_loc)*Lxyz(2)
                if (nzgrid /= 1) fp(npar_loc_old+1:npar_loc,izp) = xyz0(3)+fp(npar_loc_old+1:npar_loc,izp)*Lxyz(3)
              elseif (lspherical_coords) then
                if (nxgrid /= 1) fp(npar_loc_old+1:npar_loc,ixp) = rr_tmp(npar_loc_old+1:npar_loc)
                if (nygrid /= 1) fp(npar_loc_old+1:npar_loc,iyp) = xyz0(2)+az_tmp(npar_loc_old+1:npar_loc)*Lxyz(2)
                if (nzgrid /= 1) fp(npar_loc_old+1:npar_loc,izp) = xyz0(3)+fp(npar_loc_old+1:npar_loc,izp)*Lxyz(3)
              endif
!
            case ('birthring')
              if (birthring_width > tini) then
                if (lgaussian_birthring) then
                  do k = npar_loc_old+1,npar_loc
                    call normal_deviate(rr_tmp(k))
                  enddo
                else
                  call random_number_wrapper(rr_tmp(npar_loc_old+1:npar_loc))
                endif
                rr_tmp(npar_loc_old+1:npar_loc) = birthring_r+rr_tmp(npar_loc_old+1:npar_loc)*birthring_width
              else
                rr_tmp(npar_loc_old+1:npar_loc) = birthring_r
              endif
              call random_number_wrapper(az_tmp(npar_loc_old+1:npar_loc))
              az_tmp(npar_loc_old+1:npar_loc) = -pi+az_tmp(npar_loc_old+1:npar_loc)*2.0*pi
              if (lcartesian_coords) then
                fp(npar_loc_old+1:npar_loc,ixp) = rr_tmp(npar_loc_old+1:npar_loc)*cos(az_tmp(npar_loc_old+1:npar_loc))
                fp(npar_loc_old+1:npar_loc,iyp) = rr_tmp(npar_loc_old+1:npar_loc)*sin(az_tmp(npar_loc_old+1:npar_loc))
                fp(npar_loc_old+1:npar_loc,izp) = 0.0
              else
                fp(npar_loc_old+1:npar_loc,ixp) = rr_tmp(npar_loc_old+1:npar_loc)
                if (lcylindrical_coords) then
                  fp(npar_loc_old+1:npar_loc,iyp) = az_tmp(npar_loc_old+1:npar_loc)
                  fp(npar_loc_old+1:npar_loc,izp) = 0.0
                elseif (lspherical_coords) then
                  fp(npar_loc_old+1:npar_loc,iyp) = pi/2.0
                  fp(npar_loc_old+1:npar_loc,izp) = az_tmp(npar_loc_old+1:npar_loc)
                endif
              endif
!
            case ('nothing')
              if (j == 1) print*, 'insert_particles: nothing'
!
            case default
              call fatal_error_local('insert_particles','no such initxxp: "'//trim(initxxp(j))//'"')
!
            endselect
          enddo
!
!  Initial particle velocity.
!
          do j = 1,ninit
            select case (initvvp(j))
            case ('nothing')
              if (j == 1) print*, 'insert_particles: No particle velocity set'
!
            case ('Keplerian','keplerian')
              if (lcartesian_coords) then
                rr_tmp(npar_loc_old+1:npar_loc) = sqrt(fp(npar_loc_old+1:npar_loc,ixp)**2+ &
                    fp(npar_loc_old+1:npar_loc,iyp)**2+fp(npar_loc_old+1:npar_loc,izp)**2)
                OO_tmp(npar_loc_old+1:npar_loc) = sqrt(gravr)*rr_tmp(npar_loc_old+1:npar_loc)**(-1.5)
                fp(npar_loc_old+1:npar_loc,ivpx) = -OO_tmp(npar_loc_old+1:npar_loc)*fp(npar_loc_old+1:npar_loc,iyp)
                fp(npar_loc_old+1:npar_loc,ivpy) =  OO_tmp(npar_loc_old+1:npar_loc)*fp(npar_loc_old+1:npar_loc,ixp)
                fp(npar_loc_old+1:npar_loc,ivpz) =  0.0
              elseif (lcylindrical_coords) then
                rr_tmp(npar_loc_old+1:npar_loc) = fp(npar_loc_old+1:npar_loc,ixp)
                OO_tmp(npar_loc_old+1:npar_loc) = sqrt(gravr)*rr_tmp(npar_loc_old+1:npar_loc)**(-1.5)
                fp(npar_loc_old+1:npar_loc,ivpx) =  0.0
                fp(npar_loc_old+1:npar_loc,ivpy) =  OO_tmp(npar_loc_old+1:npar_loc)*rr_tmp(npar_loc_old+1:npar_loc)
                fp(npar_loc_old+1:npar_loc,ivpz) =  0.0
              elseif (lspherical_coords) then
                rr_tmp(npar_loc_old+1:npar_loc) = fp(npar_loc_old+1:npar_loc,ixp)*sin(fp(npar_loc_old+1:npar_loc,iyp))
                OO_tmp(npar_loc_old+1:npar_loc) = sqrt(gravr)*rr_tmp(npar_loc_old+1:npar_loc)**(-1.5)
                fp(npar_loc_old+1:npar_loc,ivpx) =  0.0
                fp(npar_loc_old+1:npar_loc,ivpy) =  0.0
                fp(npar_loc_old+1:npar_loc,ivpz) =  OO_tmp(npar_loc_old+1:npar_loc)*rr_tmp(npar_loc_old+1:npar_loc)
              endif
!
            case ('constant')
              if (lcylindrical_coords) then
                fp(npar_loc_old+1:npar_loc,ivpx) &
                    = vpx0*cos(fp(npar_loc_old+1:npar_loc,iyp)) &
                     +vpy0*sin(fp(npar_loc_old+1:npar_loc,iyp))
                fp(npar_loc_old+1:npar_loc,ivpy) &
                    = vpy0*cos(fp(npar_loc_old+1:npar_loc,iyp)) &
                     -vpx0*sin(fp(npar_loc_old+1:npar_loc,iyp))
                fp(npar_loc_old+1:npar_loc,ivpz) = vpz0
              else
                fp(npar_loc_old+1:npar_loc,ivpx) = vpx0
                fp(npar_loc_old+1:npar_loc,ivpy) = vpy0
                fp(npar_loc_old+1:npar_loc,ivpz) = vpz0
              endif
!
            case ('zero')
!              if (lroot) print*, 'insert_particles: Zero particle velocity'
              fp(1:npar_loc,ivpx:ivpz) = 0.0
!
            case ('follow-gas')
              if (lroot) &
                  print*, 'insert_particles: Particle velocity equal to gas velocity'
              do k = 1,npar_loc
                call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:),0,0)
                fp(k,ivpx:ivpz) = uup
              enddo
!
            case default
              call fatal_error_local('insert_particles','no such initvvp: "'//trim(initvvp(j))//'"')
            endselect
!
          enddo ! do j=1,ninit
!
!  Initialize particle radius
!
          if (lparticles_radius) call set_particle_radius(f,fp,npar_loc_old+1,npar_loc)
          if (lparticles_number) call set_particle_number(f,fp,npar_loc_old+1,npar_loc)
          if (lbirthring_depletion) fp(npar_loc_old+1:npar_loc,ibrtime) = 0.0
!
!  Particles are not allowed to be present in non-existing dimensions.
!  This would give huge problems with interpolation later.
!
          if (nxgrid == 1) fp(npar_loc_old+1:npar_loc,ixp) = x(nghost+1)
          if (nygrid == 1) fp(npar_loc_old+1:npar_loc,iyp) = y(nghost+1)
          if (nzgrid == 1) fp(npar_loc_old+1:npar_loc,izp) = z(nghost+1)
!
          if (lparticles_diagnos_state) call insert_particles_diagnos_state(fp, npar_loc_old)
!
        endif
      endif ! if (lroot) then
!
!  Redistribute particles only when t < max_particle_insert_time
!  and t>tstart_insert_particles.
!  Could have included some other tests here aswell......
!
      if (t < max_particle_insert_time .and. t > tstart_insert_particles) then
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
      if (lbirthring_depletion) then
        if (lcartesian_coords) then
          rr_tmp(1:npar_loc) = sqrt(fp(1:npar_loc,ixp)**2.0+fp(1:npar_loc,iyp)**2.0)
        else
          rr_tmp(1:npar_loc) = fp(1:npar_loc,ixp)
        endif
        do k = 1,npar_loc
          if (lgaussian_birthring) then
            fp(k,ibrtime) = fp(k,ibrtime) + &
                dt*exp(-((rr_tmp(k)-birthring_r)**2.0)/(2.0*birthring_width**2.0))
          else
            if ((rr_tmp(k)  >=  birthring_r-birthring_width) .and. &
                (rr_tmp(k)  <=  birthring_r+birthring_width)) &
                fp(k,ibrtime) = fp(k,ibrtime)+dt
          endif
          if (fp(k,ibrtime)  >=  birthring_lifetime) call remove_particle(fp,ipar,k)
        enddo
      endif
!
    endsubroutine insert_particles
!***********************************************************************
    subroutine insert_nucleii(f,fp,ineargrid,df)
!
! Insert particles nucleii continuously (when lnucleation == T),
! A particle is inserted whenever the mass fraction of a scalar
! that corresponds to the nucleated mass is above a certain threshold
! in a grid cell.
!
      use General, only: random_number_wrapper, normal_deviate
      use Particles_diagnos_state, only: insert_particles_diagnos_state
      use Mpicomm, only: mpireduce_sum_int, mpibarrier, mpisend_int, mpirecv_int
      use Particles_number, only: set_particle_number
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
      real, dimension(3) :: uup
!
      logical, save :: linsertmore=.true.
      real :: xx0, yy0, r2, r, mass_nucleii, part_mass, TTp
      integer :: j, k, n_insert, npar_loc_old, iii
      integer :: ii,jj,kk
      integer :: jproc,tag_id,tag0=283
!
! In this subroutine we loop over all processors one-by-one. This takes
! time. It may therefore be beneficial not to do this at every timestep.
! The interval by which this is done is given by it_insert_nuclei.
!
      if (mod(it, it_insert_nuclei) == 0) then
!
! Insertion of particles is stopped when maximum number of particles is reached,
! unless linsert_as_many_as_possible is set.
! Maximum numer of particles allowed in system is defined by max_particles,
! initialized to npar. Note that this may cause errors at a processor further
! downstream, if particles accumulate and mpar_loc is too small.
! Since root inserts all new particles, make sure
! npar_loc + n_insert < mpar_loc
! so that a processor can not exceed its maximum number of particles.
!
         call mpireduce_sum_int(npar_loc,npar_total)
         npar_loc_old=npar_loc
!
!  Loop over all processors jproc, and only proceed of jproc==iproc.
!  Receive npar_inserted_tot from previous one, unless we are on zero.
!
         do jproc=0,ncpus-1
            if (iproc==jproc) then
               if (iproc/=0) then
                  tag_id=tag0+jproc
                  call mpirecv_int(npar_inserted_tot,mod(jproc-1,ncpus),tag_id)
               endif
               !
               ! Check if we want to insert particles
               !
               if (t < max_particle_insert_time .and. t > tstart_insert_particles) then
                  !
                  ! Loop over all grid cells to identify those where nucleii should
                  ! be inserted
                  !
                  do ii=l1,l2
                     do jj=m1,m2
                        do kk=n1,n2
                           !
                           ! Insert nucleii if scalar concentration is above threshold value
                           !
                           if (ldensity_nolog) then
                              mass_nucleii=f(ii,jj,kk,icc)*f(ii,jj,kk,irho)
                           else
                              mass_nucleii=f(ii,jj,kk,icc)*exp(f(ii,jj,kk,ilnrho))
                           endif
                           if (mass_nucleii .gt. nucleation_threshold) then
                              if (1+npar_loc <= mpar_loc) then
                                 linsertmore = .true.
                              else
                                 linsertmore = .false.
                                 call fatal_error("insert_nucleii","mpar_loc is too small!")
                              endif
                              if (linsertmore) then
                                 !
                                 ! Insert nucleii:
                                 !
                                 iii = npar_loc+1
                                 ipar(iii) = npar_inserted_tot+1
                                 npar_loc = npar_loc + 1
                                 k=npar_loc
!
! Update total number of inserted particles, npar_inserted_tot.
! Not the same as npar_total, which is the number of particles in the system,
! without counting removed particles
!
                                 npar_inserted_tot = 1 + npar_inserted_tot
                                 !
                                 ! Put the particle in the center of the local grid cell
                                 !
                                 fp(k,ixp) = x(ii)
                                 fp(k,iyp) = y(jj)
                                 fp(k,izp) = z(kk)
                                 !
                                 ! Save the time when the particle was born
                                 !
                                 fp(k,iborn) = t
                                 !
                                 ! Give the particle the same velocity as the local fluid cell
                                 !
                                 ineargrid(k,1)=ii
                                 ineargrid(k,2)=jj
                                 ineargrid(k,3)=kk
                                 call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:),0,0)
                                 fp(k,ivpx:ivpz) = uup
                                 ! 
                                 !  Initialize particle radius
                                 !
                                 if (lparticles_radius) then
                                    fp(k,iap)=f(ii,jj,kk,icc+1)/f(ii,jj,kk,icc)
                                    if (lparticles_number) then
                                       part_mass=4.*pi*fp(k,iap)**3/3.*true_density_cond_spec
                                       fp(k,inpswarm)=mass_nucleii*redfrac/part_mass
                                    endif
                                 endif
                                 !
                                 ! Initialize temperature
                                 !
                                 if (lparticles_temperature) then
                                    if (ltemperature_nolog) then
                                       call interpolate_linear(f,iTT,fp(k,ixp:izp),TTp,ineargrid(k,:),0,0)
                                    else
                                       call interpolate_linear(f,ilnTT,fp(k,ixp:izp),TTp,ineargrid(k,:),0,0)
                                       TTp=exp(TTp)
                                    endif
                                    fp(k,iTp) = TTp
                                 endif
!
!  Particles are not allowed to be present in non-existing dimensions.
!  This would give huge problems with interpolation later.
!
                                 if (nxgrid == 1) fp(k,ixp) = x(nghost+1)
                                 if (nygrid == 1) fp(k,iyp) = y(nghost+1)
                                 if (nzgrid == 1) fp(k,izp) = z(nghost+1)
                                 !
                                 if (lparticles_diagnos_state) call insert_particles_diagnos_state(fp, npar_loc_old)
                                 !
                                 ! Set the scalar to zero since the nucleii have now been moved to the
                                 ! particle phase
                                 !
                                 if (lset_df_insert_nucleii) then
                                    df(ii,jj,kk,icc) = df(ii,jj,kk,icc) - redfrac*f(ii,jj,kk,icc)/dt
                                    df(ii,jj,kk,icc+1) = df(ii,jj,kk,icc+1) - redfrac*f(ii,jj,kk,icc+1)/dt
                                 else
                                    f(ii,jj,kk,icc)   = (1.-redfrac)*f(ii,jj,kk,icc)
                                    f(ii,jj,kk,icc+1) = (1.-redfrac)*f(ii,jj,kk,icc+1)
                                 endif
                              endif
                           endif
                        enddo
                     enddo
                  enddo
               endif  !(over time)
               !
               !  send to next processor, or to zero if on the last one.
               !
               tag_id=tag0+jproc+1
               call mpisend_int(npar_inserted_tot,mod(jproc+1,ncpus),tag_id)
            endif  !(iproc==jproc)
            !
            !  apply barrier, because this is sequential.
            !
            if (jproc .lt. ncpus-1) call mpibarrier
         enddo
         !
         !  root receives from last processor to be ready for the next time step
         !
         tag_id=tag0+ncpus
         if (iproc==0) call mpirecv_int(npar_inserted_tot,ncpus-1,tag_id)
         !
         !  Redistribute particles only when t < max_particle_insert_time
         !  and t>tstart_insert_particles.
         !  Could have included some other tests here aswell......
         !
         if (t < max_particle_insert_time .and. t > tstart_insert_particles) then
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
      endif
!
    endsubroutine insert_nucleii
!***********************************************************************
    subroutine streaming_coldstart(fp,f)
!
!  Mode that is unstable to the streaming instability of Youdin & Goodman (2005)
!
!  14-apr-06/anders: coded
!
      use EquationOfState, only: cs0
!
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mx,my,mz,mfarray) :: f
!
      real :: eta_vK, ampluug, dxp, dzp
      integer :: i, i1, i2, j, k, npar_loc_x, npar_loc_z
!
!  The number of particles per grid cell must be a quadratic number.
!
      if ( sqrt(npar/real(nwgrid)) /= int(sqrt(npar/real(nwgrid))) .or. &
          sqrt(npar_loc/real(nw)) /= int(sqrt(npar_loc/real(nw))) ) then
        print*, '   iproc, npar/nw, npar_loc/nwgrid=', iproc, npar/real(nwgrid), npar_loc/real(nw)
        call fatal_error('streaming_coldstart','number of particles per grid must be a square number')
      endif
!
!  Define a few disc parameters.
!
      eta_vK = -0.5 * beta_glnrho_global(1) * cs0
      if (lroot) print*, 'streaming: eta * v_K = ', eta_vK
!
!  Place particles equidistantly.
!
      npar_loc_x = sqrt(npar_loc/(Lxyz_loc(3)/Lxyz_loc(1)))
      npar_loc_z = npar_loc/npar_loc_x
      dxp = Lxyz_loc(1)/npar_loc_x
      dzp = Lxyz_loc(3)/npar_loc_z
      do i = 1,npar_loc_x
        i1 = (i-1)*npar_loc_z+1
        i2 = i*npar_loc_z
        fp(i1:i2,ixp) = xyz0_loc(1) + (real(i) - 0.5) * dxp
        do j = i1,i2
          fp(j,izp) = xyz0_loc(3)+dzp/2+(j-i1)*dzp
        enddo
      enddo
!
!  Shift particle locations slightly so that wanted mode appears.
!
      do k = 1,npar_loc
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
      do k = 1,npar_loc
        fp(k,ivpx) = fp(k,ivpx) + eta_vK*amplxxp* &
            ( real(coeff(1))*cos(kx_xxp*fp(k,ixp)) - &
             aimag(coeff(1))*sin(kx_xxp*fp(k,ixp)))*cos(kz_xxp*fp(k,izp))
        fp(k,ivpy) = fp(k,ivpy) + eta_vK*amplxxp* &
            ( real(coeff(2))*cos(kx_xxp*fp(k,ixp)) - &
             aimag(coeff(2))*sin(kx_xxp*fp(k,ixp)))*cos(kz_xxp*fp(k,izp))
        fp(k,ivpz) = fp(k,ivpz) + eta_vK*(-amplxxp)* &
            (aimag(coeff(3))*cos(kx_xxp*fp(k,ixp)) + &
              real(coeff(3))*sin(kx_xxp*fp(k,ixp)))*sin(kz_xxp*fp(k,izp))
      enddo
!
!  Change the gas velocity amplitude so that the numerical error on the drag
!  force is corrected (the error is due to the interpolation of the gas
!  velocity field to the positions of the particles). A better way to correct
!  this is to go to a quadratic interpolation scheme.
!
      ampluug = amplxxp
      if (lcoldstart_amplitude_correction) &
          ampluug = amplxxp/(1-dx**2/8*(kx_xxp**2+kz_xxp**2))
!
!  Set fluid fields.
!
      do m = m1,m2
        do n = n1,n2
          f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + &
              amplxxp* &
              ( real(coeff(7))*cos(kx_xxp*x(l1:l2)) - &
               aimag(coeff(7))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
          f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + &
              eta_vK*ampluug* &
              ( real(coeff(4))*cos(kx_xxp*x(l1:l2)) - &
               aimag(coeff(4))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
          f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
              eta_vK*ampluug* &
              ( real(coeff(5))*cos(kx_xxp*x(l1:l2)) - &
               aimag(coeff(5))*sin(kx_xxp*x(l1:l2)))*cos(kz_xxp*z(n))
!
          f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + &
              eta_vK*(-ampluug)* &
              (aimag(coeff(6))*cos(kx_xxp*x(l1:l2)) + &
                real(coeff(6))*sin(kx_xxp*x(l1:l2)))*sin(kz_xxp*z(n))
        enddo
      enddo
!
    endsubroutine streaming_coldstart
!***********************************************************************
    subroutine streaming(fp,f)
!
!  Mode that is unstable to the streaming instability of Youdin & Goodman (2005)
!
!  30-jan-06/anders: coded
!
      use General, only: random_number_wrapper
!
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mx,my,mz,mfarray) :: f
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
!
      kx = kx_xxp*Lxyz(1)
      kz = kz_xxp*Lxyz(3)
      do k = 1,npar_loc
!
        call random_number_wrapper(r)
        call random_number_wrapper(p)
!
        fprob = 1.0
        zprob = 0.0
!
        j = 0
!
!  Use Newton-Raphson iteration to invert function.
!
        do while ( abs(fprob) > 0.0001 )
!
          xprob = r
          fprob = zprob + amplxxp/kz*cos(kx*xprob)*sin(kz*zprob) - p
          dfprob = 1.0 + amplxxp*cos(kx*xprob)*cos(kz*zprob)
          dzprob = -fprob/dfprob
          zprob = zprob+0.2*dzprob
!
          j = j+1
!
        enddo
!
        if ( mod(k,npar_loc/100) == 0) then
          print '(i7,i3,4f11.7)', k, j, r, p, xprob, zprob
        endif
!
        fp(k,ixp) = xprob*Lxyz(1)+xyz0(1)
        fp(k,izp) = zprob*Lxyz(3)+xyz0(3)
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
        lmigration_redo_org = lmigration_redo
        lmigration_redo = .true.
        call migrate_particles(fp,ipar)
        lmigration_redo = lmigration_redo_org
      endif
!
!  Set fluid fields.
!
      do m = m1,m2
        do n = n1,n2
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
        enddo
      enddo
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
      use EquationOfState, only: cs20, get_gamma_etc
      use General, only: random_number_wrapper
!
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mx,my,mz,mfarray) :: f
!
      integer, parameter :: nz_inc=10
      real, dimension(nz_inc*nz) :: z_dense, eps
      real :: r, Hg, Hd, frac, rho1, Sigmad, Sigmad_num, Xi, fXi, dfdXi
      real :: dz_dense, eps_point, z00_dense, rho, lnrho
      integer :: nz_dense=nz_inc*nz, npar_bin
      integer :: i, i0, k
      real :: gamma

      call get_gamma_etc(gamma)
!
!  Calculate dust "scale height".
!
      rho1 = 1.0
      Hg = 1.0
      Sigmad = eps_dtog*rho1*Hg*sqrt(2*pi)
      Hd = sqrt(Ri0)*abs(beta_glnrho_scaled(1))/2*1.0
!
!  Need to find eps1 that results in given dust column density.
!
      Xi = sqrt(eps1*(2+eps1))/(1+eps1)
      fXi = -2*Xi + alog((1+Xi)/(1-Xi))-Sigmad/(Hd*rho1)
      i = 0
!
!  Newton-Raphson on equation Sigmad/(Hd*rho1)=-2*Xi + alog((1+Xi)/(1-Xi)).
!  Here Xi = sqrt(eps1*(2+eps1))/(1+eps1).
!
      do while (abs(fXi) >= 0.00001)
!
        dfdXi = 2*Xi**2/(1-Xi**2)
        Xi = Xi-0.1*fXi/dfdXi
!
        fXi = -2*Xi + alog((1+Xi)/(1-Xi))-Sigmad/(Hd*rho1)
!
        i = i+1
        if (i >= 1000) call fatal_error('constant_richardson','1000 Newton-Raphson iterations reached')
!
      enddo
!
!  Calculate eps1 from Xi.
!
      eps1 = -1+1/sqrt(-(Xi**2)+1)
      if (lroot) print*, 'constant_richardson: Hd, eps1=', Hd, eps1
!
!  Make z denser for higher resolution in density.
!
      dz_dense = Lxyz_loc(3)/nz_dense
      z00_dense = xyz0_loc(3)+0.5*dz_dense
      do n = 1,nz_dense
        z_dense(n) = z00_dense+(n-1)*dz_dense
      enddo
!
!  Dust-to-gas ratio as a function of z (with cutoff).
!
      eps = 1/sqrt(z_dense**2/Hd**2+1/(1+eps1)**2)-1
      where (eps <= 0.0) eps = 0.0
!
!  Calculate the dust column density numerically.
!
      Sigmad_num = sum(rho1*eps*dz_dense)
      if (lroot) print*, 'constant_richardson: Sigmad, Sigmad (numerical) = ', Sigmad, Sigmad_num
!
!  Place particles according to probability function.
!
      i0 = 0
      do n = 1,nz_dense
        frac = eps(n)/Sigmad_num*dz_dense
        npar_bin = int(frac*npar_loc)
        if (npar_bin >= 2 .and. mod(n,2) == 0) npar_bin = npar_bin+1
        do i = i0+1,i0+npar_bin
          if (i <= npar_loc) then
            call random_number_wrapper(r)
            fp(i,izp) = z_dense(n)+(2*r-1.0)*dz_dense/2
          endif
        enddo
        i0 = i0+npar_bin
      enddo
      if (lroot) print '(A,i7,A)', 'constant_richardson: placed ', &
          i0, ' particles according to Ri=const'
!
!  Particles left out by round off are just placed randomly.
!
      if (i0+1 <= npar_loc) then
        do k = i0+1,npar_loc
          call random_number_wrapper(r)
          fp(k,izp) = xyz0(3)+r*Lxyz(3)
        enddo
        if (lroot) print '(A,i7,A)', 'constant_richardson: placed ', npar_loc-i0, ' particles randomly.'
      endif
!
!  Random positions in x and y.
!
      do k = 1,npar_loc
        if (nxgrid /= 1) then
          call random_number_wrapper(r)
          fp(k,ixp) = r
        endif
        if (nygrid /= 1) then
          call random_number_wrapper(r)
          fp(k,iyp) = r
        endif
      enddo
      if (nxgrid /= 1) &
          fp(1:npar_loc,ixp) = xyz0_loc(1)+fp(1:npar_loc,ixp)*Lxyz_loc(1)
      if (nygrid /= 1) &
          fp(1:npar_loc,iyp) = xyz0_loc(2)+fp(1:npar_loc,iyp)*Lxyz_loc(2)
!
!  Set gas velocity according to dust-to-gas ratio and global pressure gradient.
!
      do imn = 1,ny*nz
!
        n = nn(imn)
        m = mm(imn)
!
        if (abs(z(n)) <= Hd*sqrt(1-1/(1+eps1)**2)) then
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
        rho = exp(lnrho)
!
        if (ldensity_nolog) then
          f(l1:l2,m,n,irho)  = rho
        else
          f(l1:l2,m,n,ilnrho) = lnrho
        endif
!
        eps_point = 1/sqrt(z(n)**2/Hd**2+1/(1+eps1)**2)-1
        if (eps_point <= 0.0) eps_point = 0.0
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
      do k = 1,npar_loc
!
        eps_point = 1/sqrt(fp(k,izp)**2/Hd**2+1/(1+eps1)**2)-1
        if (eps_point <= 0.0) eps_point = 0.0
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
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(nx) :: eps
      real, dimension(3) :: vvp, force
      integer :: imn, i, k, ix0, iy0, iz0
!
      if (ldragforce_stiff .and. .not. lpencil_check_at_work) then
        do imn = 1,ny*nz
          n = nn(imn)
          m = mm(imn)
          eps = f(l1:l2,m,n,irhop)/f(l1:l2,m,n,irho)
          do i = 0,2
            f(l1:l2,m,n,iux+i) = (f(l1:l2,m,n,iux+i)+eps*f(l1:l2,m,n,iupx+i) + &
                eps/(1.0+eps)*tausp*f(l1:l2,m,n,ifgx+i))/(1.0+eps)
          enddo
          f(l1:l2,m,n,ifgx:ifgz) = f(l1:l2,m,n,iux:iuz)-f(l1:l2,m,n,ifgx:ifgz)
        enddo
        call boundconds_x(f,ifgx,ifgz)
        call initiate_isendrcv_bdry(f,ifgx,ifgz)
        call finalize_isendrcv_bdry(f,ifgx,ifgz)
        call boundconds_y(f,ifgx,ifgz)
        call boundconds_z(f,ifgx,ifgz)
        do k = 1,npar_loc
          ix0 = ineargrid(k,1)
          iy0 = ineargrid(k,2)
          iz0 = ineargrid(k,3)
          if (lparticlemesh_cic) then
            call interpolate_linear(f,ifgx,ifgz, &
                fp(k,ixp:izp),vvp,ineargrid(k,:),0,ipar(k))
          elseif (lparticlemesh_tsc) then
            if (linterpolate_spline) then
              call interpolate_quadratic_spline(f,ifgx,ifgz,fp(k,ixp:izp),vvp,ineargrid(k,:),0,ipar(k))
            else
              call interpolate_quadratic(f,ifgx,ifgz,fp(k,ixp:izp),vvp,ineargrid(k,:),0,ipar(k))
            endif
          else
            vvp = f(ix0,iy0,iz0,ifgx:ifgz)
          endif
          fp(k,ivpx:ivpz) = vvp
        enddo
      !elseif (lfollow_gas) then
      !  do k = 1,npar_loc
      !    ix0 = ineargrid(k,1)
      !    iy0 = ineargrid(k,2)
      !    iz0 = ineargrid(k,3)
      !    if (lparticlemesh_cic) then
      !      call interpolate_linear(f,iux,iuz, &
      !           fp(k,ixp:izp),vvp,ineargrid(k,:),0,ipar(k))
      !    elseif (lparticlemesh_tsc) then
      !      if (linterpolate_spline) then
      !        call interpolate_quadratic_spline(f,iux,iuz,fp(k,ixp:izp),vvp,ineargrid(k,:),0,ipar(k))
      !      else
      !        call interpolate_quadratic(f,iux,iuz,fp(k,ixp:izp),vvp,ineargrid(k,:),0,ipar(k))
      !      endif
      !    else
      !      vvp = f(ix0,iy0,iz0,iux:iuz)
      !    endif
      !    fp(k,ivpx:ivpz) = vvp
      !    !
      !    ! Add Brownian motion on top of the fluid velocity
      !    !
      !    if (lbrownian_forces) then
      !      !if (lbrownian_forces_Li_Ahmadi) then
      !      !  call calc_brownian_force_Li_Ahmadi(f,fp,k,ineargrid(k,:),stocunn(k),force)
      !      !else
      !      call calc_brownian_force(f,fp,k,ineargrid(k,:),force)
      !      !endif
      !      fp(k,ivpx:ivpz) = fp(k,ivpx:ivpz) + force*dt
      !    endif
      !  enddo
      endif
!
    endsubroutine particles_dragforce_stiff
!***********************************************************************
    subroutine pencil_criteria_particles
!
!  All pencils that the Particles module depends on are specified here.
!
!  20-04-06/anders: coded
!
      if (ldragforce_gas_par) then
        lpenc_requested(i_epsp) = .true.
        lpenc_requested(i_np) = .true.
        lpenc_requested(i_rho1) = .true.
      endif
      if (ldraglaw_epstein .and. &
          (lcollisional_cooling_rms .or. lcollisional_cooling_twobody .or. lcompensate_friction_increase)) then
        lpenc_requested(i_cs2)=.true.
        lpenc_requested(i_rho)=.true.
      endif
      if (ldragforce_heat .or. lcollisional_heat) then
        lpenc_requested(i_TT1) = .true.
        lpenc_requested(i_rho1) = .true.
      endif
      if (lcollisional_cooling_taucool) then
        lpenc_requested(i_np) = .true.
      endif
      if (lcollisional_cooling_rms) then
        lpenc_requested(i_epsp) = .true.
      endif
      if (lcollisional_cooling_rms .or. lcollisional_dragforce_cooling) then
        lpenc_requested(i_np) = .true.
        lpenc_requested(i_rho1) = .true.
      endif
      if (ldraglaw_epstein_transonic .or. &
          ldraglaw_eps_stk_transonic) then
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_cs2)=.true.
      endif
      if (ldragforce_stiff) then
        lpenc_requested(i_fpres) = .true.
        lpenc_requested(i_jxbr) = .true.
        lpenc_requested(i_fvisc) = .true.
      endif
      if (lthermophoretic_forces) then
        lpenc_requested(i_gTT) = .true.
      endif
!
      if (lascalar) then
        if (ltauascalar) lpenc_requested(i_tauascalar) = .true.
        lpenc_requested(i_condensationRate) = .true.
        lpenc_requested(i_waterMixingRatio) = .true.
        lpenc_requested(i_ssat) = .true.
      endif
!
      if (idiag_npm /= 0 .or. idiag_np2m /= 0 .or. idiag_npmax /= 0 .or. &
          idiag_npmin /= 0 .or. idiag_npmx /= 0 .or. idiag_npmy /= 0 .or. &
          idiag_npmz /= 0 .or. idiag_nparpmax /= 0) lpenc_diagnos(i_np) = .true.
      if (idiag_rhopm /= 0 .or. idiag_rhoprms /= 0 .or. idiag_rhop2m /= 0 .or. &
          idiag_rhopmax /= 0 .or. idiag_rhopmin /= 0 .or. idiag_rhopmphi /= 0 .or. &
          idiag_rhopmx /= 0 .or. idiag_rhopmy /= 0 .or. idiag_rhopmz /= 0) &
          lpenc_diagnos(i_rhop) = .true.
      if (idiag_rhop2mx /= 0 .or. idiag_rhop2my /= 0 .or. idiag_rhop2mz /= 0) lpenc_diagnos(i_rhop) = .true.
      if (idiag_dedragp /= 0 .or. idiag_decollp /= 0) then
        lpenc_diagnos(i_TT1) = .true.
        lpenc_diagnos(i_rho1) = .true.
      endif
      if (idiag_epspmx /= 0 .or. idiag_epspmy /= 0 .or. idiag_epspmz /= 0 .or. &
          idiag_epspmin /= 0 .or. idiag_epspmax /= 0 .or. idiag_epspm /= 0) &
          lpenc_diagnos(i_epsp) = .true.
      if (idiag_rhopmxy /= 0 .or. idiag_rhopmxz /= 0 .or. idiag_rhopmphi /= 0) &
          lpenc_diagnos2d(i_rhop) = .true.
      if (idiag_npmxy /= 0 ) lpenc_diagnos2d(i_np) = .true.
      if (idiag_sigmap /= 0) lpenc_diagnos2d(i_rhop) = .true.
!
      if (idiag_rpvpxmx /= 0 .or. idiag_rpvpymx /= 0 .or. idiag_rpvpzmx /= 0 .or. &
          idiag_rpvpx2mx /= 0 .or. idiag_rpvpy2mx /= 0 .or. idiag_rpvpz2mx /= 0) then
        lpenc_diagnos(i_rhop) = .true.
        lpenc_diagnos(i_uup) = .true.
      endif
      !
      if (ltemp_equip_part_gas) lpenc_requested(i_part_heatcap) = .true.
      if (maxval(idiag_npvzmz) > 0) lpenc_requested(i_npvz) = .true.   !MR: not pencil_diagnos?
      if (maxval(idiag_npvz2mz) > 0) lpenc_requested(i_npvz2) = .true.
      if (maxval(idiag_npuzmz) > 0) lpenc_requested(i_npuz) = .true.
      if (maxval(idiag_nptz) > 0)   lpenc_requested(i_np_rad) = .true.
      if (idiag_Shm /= 0)  lpenc_requested(i_sherwood) = .true.
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
      if (lpencil_in(i_rhop) .and. irhop == 0) then
        lpencil_in(i_np) = .true.
      endif
!
      if (lpencil_in(i_epsp)) then
        lpencil_in(i_rhop) = .true.
        lpencil_in(i_rho1) = .true.
      endif
!
      if (lpencil_in(i_uup) .and. iuup == 0) &
        call fatal_error("pencil_interdep_particles","p%uup requested but not calculated")
!
      if (lascalar .and. ltauascalar) lpencil_in(i_tauascalar) = .true.
      if (lascalar) lpencil_in(i_condensationRate) = .true.
      if (lascalar) lpencil_in(i_waterMixingRatio) = .true.
      if (lascalar) lpencil_in(i_ssat) = .true.
!
      if (lpencil_in(i_grhop)) then
        if (irhop /= 0) then
          if (nprocx /= 1.and.(.not.lcommunicate_rhop)) &
            call fatal_error("pencil_interdep_particles", &
                             "Switch on lcommunicate_rhop=T in particles_run_pars")
        else
          if (nprocx /= 1.and.(.not.lcommunicate_np)) &
            call fatal_error("pencil_interdep_particles", &
                             "Switch on lcommunicate_np=T in particles_run_pars")
        endif
      endif
!
    endsubroutine pencil_interdep_particles
!***********************************************************************
    subroutine calc_pencils_particles(f,p,fp,ineargrid)
!
      use Sub, only: grad
!
!  Calculate Particles pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  16-feb-06/anders: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      integer :: k, inx0
      type (pencil_case) :: p
!
      if (lpencil(i_np)) then
        if (inp /= 0) then
          p%np = f(l1:l2,m,n,inp)
        else
          p%np = 0.0
        endif
      endif
!
      if (lpencil(i_rhop)) then
        if (irhop /= 0) then
          p%rhop = f(l1:l2,m,n,irhop)
        else
          p%rhop = rhop_swarm*f(l1:l2,m,n,inp)
        endif
      endif
!
      if (lpencil(i_grhop)) then
        if (irhop /= 0) then
          call grad(f,irhop,p%grhop)
        else
          call grad(f,inp,p%grhop)
          p%grhop = rhop_swarm*p%grhop
        endif
      endif
!
      if (lpencil(i_epsp)) p%epsp = p%rhop*p%rho1
!
      if (lpencil(i_uup) .and. iuup > 0) p%uup = f(l1:l2,m,n,iupx:iupz)
!
      if (ipeh > 0) p%peh = f(l1:l2,m,n,ipeh)
!
      if (itauascalar > 0) p%tauascalar = f(l1:l2,m,n,itauascalar)
      if (icondensationRate > 0) p%condensationRate = f(l1:l2,m,n,icondensationRate)
      if (iwaterMixingRatio > 0) p%waterMixingRatio = f(l1:l2,m,n,iwaterMixingRatio)
      if (lpencil(i_part_heatcap) .and. iap/=0) then
        p%part_heatcap=0
        do k = 1,npar_loc
          inx0 = ineargrid(k,1)-nghost
          if (irhopswarm /= 0) then
            p%part_heatcap(inx0) = p%part_heatcap(inx0)+4*pi*fp(k,iap)**3*rhopmat*fp(k,irhopswarm)/3.
          else
            p%part_heatcap(inx0) = p%part_heatcap(inx0)+4*pi*fp(k,iap)**3*rhopmat/3.
          endif
        enddo
      endif
!
    endsubroutine calc_pencils_particles
!***********************************************************************
    subroutine dxxp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of dust particle position.
!
!  02-jan-05/anders: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mpar_loc,mpvar), intent(inout) :: dfp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
      logical :: lheader, lfirstcall=.true.
!
!  Print out header information in first time step.
!
      lheader = lfirstcall .and. lroot
!
!  Identify module and boundary conditions.
!
      if (lheader) then
        print*,'dxxp_dt: Calculate dxxp_dt'
        print*, 'dxxp_dt: Particles boundary condition bcpx=', bcpx
        print*, 'dxxp_dt: Particles boundary condition bcpy=', bcpy
        print*, 'dxxp_dt: Particles boundary condition bcpz=', bcpz
        print*, 'dxxp_dt: Set rate of change of particle '// &
                'position equal to particle velocity.'
      endif
!
!  The rate of change of a particle's position is the particle's velocity.
!  If pointmasses are used do the evolution in that module, for better
!  conservation of the Jacobi constant.
!
      if (.not. lpointmasses) then
        if (lcartesian_coords) then
!
          if (nxgrid /= 1 .or. lnocollapse_xdir_onecell) &
              dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
          if (nygrid /= 1 .or. lnocollapse_ydir_onecell) &
              dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + fp(1:npar_loc,ivpy)
          if (nzgrid /= 1 .or. lnocollapse_zdir_onecell) &
              dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
        elseif (lcylindrical_coords) then
!
          if (nxgrid /= 1 .or. lnocollapse_xdir_onecell) &
              dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
          if (nygrid /= 1 .or. lnocollapse_ydir_onecell) &
              dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + &
              fp(1:npar_loc,ivpy)/max(fp(1:npar_loc,ixp),tini)
          if (nzgrid /= 1 .or. lnocollapse_zdir_onecell) &
              dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + fp(1:npar_loc,ivpz)
!
        elseif (lspherical_coords) then
!
          if (nxgrid /= 1 .or. lnocollapse_xdir_onecell) &
              dfp(1:npar_loc,ixp) = dfp(1:npar_loc,ixp) + fp(1:npar_loc,ivpx)
          if (nygrid /= 1 .or. lnocollapse_ydir_onecell) &
              dfp(1:npar_loc,iyp) = dfp(1:npar_loc,iyp) + &
              fp(1:npar_loc,ivpy)/max(fp(1:npar_loc,ixp),tini)
          if (nzgrid /= 1 .or. lnocollapse_zdir_onecell) &
              dfp(1:npar_loc,izp) = dfp(1:npar_loc,izp) + &
              fp(1:npar_loc,ivpz)/(max(fp(1:npar_loc,ixp),tini)*sin(fp(1:npar_loc,iyp)))
        endif
      endif
!
!  With shear there is an extra term due to the background shear flow.
!
      if (lshear .and. nygrid /= 1) dfp(1:npar_loc,iyp) = &
          dfp(1:npar_loc,iyp) - qshear*Omega*fp(1:npar_loc,ixp)
!
      if (lfirstcall) lfirstcall = .false.
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
      use Particles_caustics, only: dcaustics_dt
      use Particles_tetrad, only: dtetrad_dt
      use Particles_potential, only: dvvp_dt_potential
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      real, dimension(mpar_loc,mpvar), intent(inout) :: dfp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
      real :: Omega2
      integer :: npar_found
      logical :: lheader, lfirstcall=.true.
!
!  Print out header information in first time step.
!
      lheader = lfirstcall .and. lroot
      if (lheader) print*,'dvvp_dt: Calculate dvvp_dt'
!
! Calculate the acceleration of particles only if lfollow_gas is not true
!
      if (.not. lfollow_gas) then
!
!  Add Coriolis force from rotating coordinate frame.
!
        if (Omega /= 0.) then
          if (lcoriolis_force_par) then
            if (lheader) print*,'dvvp_dt: Add Coriolis force; Omega=', Omega
            Omega2 = 2*Omega
            if (.not. lspherical_coords) then
              dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + Omega2*fp(1:npar_loc,ivpy)
              dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) - Omega2*fp(1:npar_loc,ivpx)
            else
              call not_implemented('dvvp_dt','Coriolis force on particles for spherical coordinates')
            endif
          endif
!
!  Add centrifugal force.
!
          if (lcentrifugal_force_par) then
            if (lheader) print*,'dvvp_dt: Add Centrifugal force; Omega=', Omega
            if (lcartesian_coords) then
              !
              dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + Omega**2*fp(1:npar_loc,ixp)
              !
              dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) + Omega**2*fp(1:npar_loc,iyp)
              !
            elseif (lcylindrical_coords) then
              dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + Omega**2*fp(1:npar_loc,ixp)
            else
              call not_implemented('dvvp_dt','centrifugal force on particles for spherical coordinates')
            endif
          endif
!
!  With shear there is an extra term due to the background shear flow.
!
          if (lshear .and. lshear_accel_par) &
               dfp(1:npar_loc,ivpy) = dfp(1:npar_loc,ivpy) + qshear * Omega * fp(1:npar_loc,ivpx)
        endif
!
!  Add constant background pressure gradient beta=alpha*H0/r0, where alpha
!  comes from a global pressure gradient P = P0*(r/r0)^alpha.
!  (the term must be added to the dust equation of motion when measuring
!  velocities relative to the shear flow modified by the globalfpressure grad.)
!
        if (beta_dPdr_dust /= 0.0 .and. t >= tstart_dragforce_par) &
             dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + cs20*beta_dPdr_dust_scaled
!
!  Gravity on the particles.
!  Gravity on particles is implemented only if lparticle_gravity is true which is the default.
!
        if (lparticle_gravity) call particle_gravity(f,df,fp,dfp,ineargrid)
!
        if (lcorotational_frame) call indirect_inertial_particles(f,df,fp,dfp,ineargrid)
!
!  The auxiliary has to be set to zero afterwards
!
        call diffuse_backreaction(f,df)
    endif
!
!  Diagnostic output
!
    if (ldiagnos) then
        call sum_name(npar_loc,idiag_nparsum)
        if (idiag_nparmin /= 0) call max_name(-npar_loc,idiag_nparmin,lneg=.true.)
        call max_name(+npar_loc,idiag_nparmax)
        if (idiag_nparpmax /= 0) call max_name(maxval(npar_imn),idiag_nparpmax)
        if (ixp/=0) then
          call sum_par_name(fp(1:npar_loc,ixp),idiag_xpm)
          call max_par_name(fp(1:npar_loc,ixp),idiag_xpmax)
        endif
        if (iyp/=0) then
          call sum_par_name(fp(1:npar_loc,iyp),idiag_ypm)
          call max_par_name(fp(1:npar_loc,iyp),idiag_ypmax)
        endif
        if (izp/=0) then
          call sum_par_name(fp(1:npar_loc,izp),idiag_zpm)
          call max_par_name(fp(1:npar_loc,izp),idiag_zpmax)
        endif
        if (idiag_xpmin /= 0)  call max_par_name(-fp(1:npar_loc,ixp),idiag_xpmin,lneg=.true.)
        if (idiag_ypmin /= 0)  call max_par_name(-fp(1:npar_loc,iyp),idiag_ypmin,lneg=.true.)
        if (idiag_zpmin /= 0)  call max_par_name(-fp(1:npar_loc,izp),idiag_zpmin,lneg=.true.)
        if (idiag_xp2m /= 0) call sum_par_name(fp(1:npar_loc,ixp)**2,idiag_xp2m)
        if (idiag_yp2m /= 0) call sum_par_name(fp(1:npar_loc,iyp)**2,idiag_yp2m)
        if (idiag_zp2m /= 0) call sum_par_name(fp(1:npar_loc,izp)**2,idiag_zp2m)
        if (idiag_rpm /= 0)  call sum_par_name(sqrt(fp(1:npar_loc,ixp)**2+ &
            fp(1:npar_loc,iyp)**2+fp(1:npar_loc,izp)**2),idiag_rpm)
        if (idiag_rp2m /= 0) call sum_par_name(fp(1:npar_loc,ixp)**2+ &
            fp(1:npar_loc,iyp)**2+fp(1:npar_loc,izp)**2,idiag_rp2m)
        if (ivpx/=0) call sum_par_name(fp(1:npar_loc,ivpx),idiag_vpxm)
        if (ivpy/=0) call sum_par_name(fp(1:npar_loc,ivpy),idiag_vpym)
        if (ivpz/=0) call sum_par_name(fp(1:npar_loc,ivpz),idiag_vpzm)
        if (idiag_vpxvpym /= 0) call sum_par_name( &
            fp(1:npar_loc,ivpx)*fp(1:npar_loc,ivpy),idiag_vpxvpym)
        if (idiag_vpxvpzm /= 0) call sum_par_name( &
            fp(1:npar_loc,ivpx)*fp(1:npar_loc,ivpz),idiag_vpxvpzm)
        if (idiag_vpyvpzm /= 0) call sum_par_name( &
            fp(1:npar_loc,ivpy)*fp(1:npar_loc,ivpz),idiag_vpyvpzm)
        if (idiag_lpxm /= 0) call sum_par_name( &
            fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpz)- &
            fp(1:npar_loc,izp)*fp(1:npar_loc,ivpy),idiag_lpxm)
        if (idiag_lpym /= 0) call sum_par_name( &
            fp(1:npar_loc,izp)*fp(1:npar_loc,ivpx)- &
            fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpz),idiag_lpym)
        if (idiag_lpzm /= 0) call sum_par_name( &
            fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpy)- &
            fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpx),idiag_lpzm)
        if (idiag_lpx2m /= 0) call sum_par_name( &
            (fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpz)- &
            fp(1:npar_loc,izp)*fp(1:npar_loc,ivpy))**2,idiag_lpx2m)
        if (idiag_lpy2m /= 0) call sum_par_name( &
            (fp(1:npar_loc,izp)*fp(1:npar_loc,ivpx)- &
            fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpz))**2,idiag_lpy2m)
        if (idiag_lpz2m /= 0) call sum_par_name( &
            (fp(1:npar_loc,ixp)*fp(1:npar_loc,ivpy)- &
            fp(1:npar_loc,iyp)*fp(1:npar_loc,ivpx))**2,idiag_lpz2m)
        if (idiag_vpx2m /= 0) &
            call sum_par_name(fp(1:npar_loc,ivpx)**2,idiag_vpx2m)
        if (idiag_vpy2m /= 0) &
            call sum_par_name(fp(1:npar_loc,ivpy)**2,idiag_vpy2m)
        if (idiag_vpz2m /= 0) &
            call sum_par_name(fp(1:npar_loc,ivpz)**2,idiag_vpz2m)
        if (idiag_vprms /= 0) &
            call sum_par_name((fp(1:npar_loc,ivpx)**2 &
            +fp(1:npar_loc,ivpy)**2 &
            +fp(1:npar_loc,ivpz)**2),idiag_vprms,lsqrt = .true.)
        if (idiag_vpyfull2m /= 0) &
            call sum_par_name((fp(1:npar_loc,ivpy)- &
            qshear*Omega*fp(1:npar_loc,ixp))**2,idiag_vpyfull2m)
        if (idiag_ekinp /= 0) then
          if (lparticles_density) then
            call sum_par_name(0.5*fp(1:npar_loc,irhopswarm)* &
                sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
          else
            if (mp_swarm == 0.) then
              if (lparticles_radius) then
                call sum_par_name(0.5*(4./3.)*pi*(fp(1:npar_loc,iap)**3)*rhopmat* &
                    ( fp(1:npar_loc,ivpx)**2 &
                    +fp(1:npar_loc,ivpy)**2 &
                    +fp(1:npar_loc,ivpz)**2),idiag_ekinp)
              endif
            else
              if (lcartesian_coords .and.(all(lequidist))) then
                call sum_par_name(0.5*rhop_swarm*npar_per_cell* &
                    sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
              else
                call sum_par_name(0.5*mp_swarm*sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2),idiag_ekinp)
              endif
            endif
          endif
        endif
        if (idiag_epotpm /= 0) call sum_par_name( &
            -gravr/sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_epotpm)
        ! Better formulate with lsqrt
        if (idiag_vpmax /= 0) call max_par_name( &
            sqrt(sum(fp(1:npar_loc,ivpx:ivpz)**2,2)),idiag_vpmax)
        if (idiag_vrelpabsm /= 0) call calc_relative_velocity(f,fp,ineargrid)
        if (ivpx/=0) call max_par_name(fp(1:npar_loc,ivpx),idiag_vpxmax)
        if (ivpy/=0) call max_par_name(fp(1:npar_loc,ivpy),idiag_vpymax)
        if (ivpz/=0) call max_par_name(fp(1:npar_loc,ivpz),idiag_vpzmax)
        if (idiag_vpxmin /= 0) call max_par_name(-fp(1:npar_loc,ivpx),idiag_vpxmin,lneg=.true.)
        if (idiag_vpymin /= 0) call max_par_name(-fp(1:npar_loc,ivpy),idiag_vpymin,lneg=.true.)
        if (idiag_vpzmin /= 0) call max_par_name(-fp(1:npar_loc,ivpz),idiag_vpzmin,lneg=.true.)
        if (idiag_eccpxm /= 0) call sum_par_name( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,ixp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpx)/gravr-fp(1:npar_loc,ixp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_eccpxm)
        if (idiag_eccpym /= 0) call sum_par_name( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,iyp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpy)/gravr-fp(1:npar_loc,iyp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_eccpym)
        if (idiag_eccpzm /= 0) call sum_par_name( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,izp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpz)/gravr-fp(1:npar_loc,izp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)),idiag_eccpzm)
        if (idiag_eccpx2m /= 0) call sum_par_name(( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,ixp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpx)/gravr-fp(1:npar_loc,ixp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)))**2,idiag_eccpx2m)
        if (idiag_eccpy2m /= 0) call sum_par_name(( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,iyp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpy)/gravr-fp(1:npar_loc,iyp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)))**2,idiag_eccpy2m)
        if (idiag_eccpz2m /= 0) call sum_par_name(( &
            sum(fp(1:npar_loc,ivpx:ivpz)**2,dim=2)*fp(1:npar_loc,izp)/gravr- &
            sum(fp(1:npar_loc,ixp:izp)*fp(1:npar_loc,ivpx:ivpz),dim=2)* &
            fp(1:npar_loc,ivpz)/gravr-fp(1:npar_loc,izp)/ &
            sqrt(sum(fp(1:npar_loc,ixp:izp)**2,dim=2)))**2,idiag_eccpz2m)
        if (idiag_rhopvpxm /= 0) then
          if (lparticles_density) then
            call sum_par_name(fp(1:npar_loc,irhopswarm)*fp(1:npar_loc,ivpx), &
                idiag_rhopvpxm)
          elseif (lparticles_radius .and. lparticles_number) then
            call sum_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpx),idiag_rhopvpxm)
          endif
        endif

        if (idiag_rhopart /= 0) then
          if (lparticles_radius .and. lparticles_number) then
            call sum_par_name(four_pi_rhopmat_over_three* &
                 fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)*npar_loc,idiag_rhopart,len=npar_loc)
          endif
        endif
        
        if (idiag_rhopvpym /= 0) then
          if (lparticles_density) then
            call sum_par_name(fp(1:npar_loc,irhopswarm)*fp(1:npar_loc,ivpy), &
                idiag_rhopvpym)
          elseif (lparticles_radius .and. lparticles_number) then
            call sum_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpy),idiag_rhopvpym)
          endif
        endif
        if (idiag_rhopvpzm /= 0) then
          if (lparticles_density) then
            call sum_par_name(fp(1:npar_loc,irhopswarm)*fp(1:npar_loc,ivpz), &
                idiag_rhopvpzm)
          elseif (lparticles_radius .and. lparticles_number) then
            call sum_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpz),idiag_rhopvpzm)
          endif
        endif
        if (idiag_rhopvpxt /= 0) then
          if (lparticles_density) then
            call integrate_par_name(fp(1:npar_loc,irhopswarm)* &
                fp(1:npar_loc,ivpx),idiag_rhopvpxt)
          elseif (lparticles_radius .and. lparticles_number) then
            call integrate_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpx),idiag_rhopvpxt)
          endif
        endif
        if (idiag_rhopvpyt /= 0) then
          if (lparticles_density) then
            call integrate_par_name(fp(1:npar_loc,irhopswarm)* &
                fp(1:npar_loc,ivpy),idiag_rhopvpyt)
          elseif (lparticles_radius .and. lparticles_number) then
            call integrate_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpy),idiag_rhopvpyt)
          endif
        endif
        if (idiag_rhopvpzt /= 0) then
          if (lparticles_density) then
            call integrate_par_name(fp(1:npar_loc,irhopswarm)* &
                fp(1:npar_loc,ivpz),idiag_rhopvpzt)
          elseif (lparticles_radius .and. lparticles_number) then
            call integrate_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                fp(1:npar_loc,ivpz),idiag_rhopvpzt)
          endif
        endif
        if (idiag_rhopvpysm /= 0) then
          if (lparticles_density) then
            call sum_par_name(fp(1:npar_loc,irhopswarm)*Sshear*fp(1:npar_loc,ixp),idiag_rhopvpysm)
          elseif (lparticles_radius .and. lparticles_number) then
            call sum_par_name(four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)* &
                Sshear*fp(1:npar_loc,ixp),idiag_rhopvpysm)
          endif
        endif
        if (idiag_mpt /= 0) then
          if (lparticles_density) then
            call integrate_par_name((/fp(1:npar_loc,irhopswarm)/),idiag_mpt)
          elseif (lparticles_radius .and. lparticles_number) then
            call integrate_par_name((/four_pi_rhopmat_over_three* &
                fp(1:npar_loc,iap)**3*fp(1:npar_loc,inpswarm)/),idiag_mpt)
          endif
        endif
        if (idiag_npargone /= 0) then
          call count_particles(ipar,npar_found)
          call save_name(float(npar-npar_found),idiag_npargone)
        endif
        if (idiag_npar_found /= 0) then
          call count_particles(ipar,npar_found)
          call save_name(float(npar_found),idiag_npar_found)
        endif
        if (idiag_deshearbcsm /= 0) &
          call sum_name(energy_gain_shear_bcs/npar,idiag_deshearbcsm) !MR: What to sum here?

      endif
      if (lfirstcall) lfirstcall = .false.
!
!  If we are using caustics :
!
      if (lparticles_caustics) call dcaustics_dt(f,df,fp,dfp,ineargrid)
!
!  If we are using tetrad :
!
      if (lparticles_tetrad) call dtetrad_dt(f,df,fp,dfp,ineargrid)
!
! and potential
!
      if (lparticles_potential) call dvvp_dt_potential(f,df,fp,dfp,ineargrid)
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
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mpar_loc,mpvar), intent(inout) :: dfp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      real, dimension(mpar_loc) :: rpbeta_tmp_arr, OO2_arr, rr_arr, vv_arr
      real, dimension(mpar_loc,3) :: gpp_arr
      integer, dimension(mpar_loc) :: jspec_arr
!
      real, dimension(3) :: ggp
      real :: rr=0, vv=0, OO2, rpbeta_tmp=0
      integer :: k, jspec
      logical :: lheader, lfirstcall=.true.
!
      call keep_compiler_quiet(f,df)
!
!  Print out header information in first time step.
!
      lheader = lfirstcall .and. lroot
      if (lheader) print*,'particle_gravity: Calculating gravity'
!
!  Gravity in the x-direction.
!
      if (t >= tstart_grav_x_par) then
!
        select case (gravx_profile)
!
        case ('','zero')
          if (lheader) print*, 'particle_gravity: No gravity in x-direction.'
!
        case ('const','plain')
          if (lheader) print*, 'particle_gravity: Constant gravity field in x-direction'
          dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + gravx
!
        case ('linear')
          if (lheader) print*, 'particle_gravity: Linear gravity field in x-direction.'
          dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) - nu_epicycle2*fp(1:npar_loc,ixp)
!
        case ('sinusoidal')
          if (lheader) &
              print*, 'particle_gravity: Sinusoidal gravity field in x-direction.'
          dfp(1:npar_loc,ivpx) = dfp(1:npar_loc,ivpx) + gravx*sin(kx_gg*fp(1:npar_loc,ixp))
!
        case default
          call fatal_error('particle_gravity','no such gravx_profile: '//trim(gravx_profile))
        endselect
!
      endif
!
!  Gravity in the z-direction.
!
      if (t >= tstart_grav_z_par) then
!
        select case (gravz_profile)
!
        case ('','zero')
          if (lheader) print*, 'particle_gravity: No gravity in z-direction.'
!
        case ('const','plain')
          if (lheader) print*, 'particle_gravity: Constant gravity field in z-direction.'
          dfp(1:npar_loc,ivpz) = dfp(1:npar_loc,ivpz) + gravz
!
        case ('linear')
          if (lheader) print*, 'particle_gravity: Linear gravity field in z-direction.'
          dfp(1:npar_loc,ivpz) = dfp(1:npar_loc,ivpz) - nu_epicycle2*fp(1:npar_loc,izp)
!
        case ('sinusoidal')
          if (lheader) print*, 'particle_gravity: Sinusoidal gravity field in z-direction.'
          dfp(1:npar_loc,ivpz) = dfp(1:npar_loc,ivpz) + gravz*sin(kz_gg*fp(1:npar_loc,izp))
!
        case default
          call fatal_error('particle_gravity','no such gravz_profile: '//trim(gravz_profile))
        endselect
!
      endif
!
!  Radial gravity.
!
      if (t >= tstart_grav_r_par) then
!
        select case (gravr_profile)
!
        case ('','zero')
          if (lheader) print*, 'particle_gravity: No radial gravity'
!
        case ('newtonian-central','newtonian')
          if (lheader) print*, 'particle_gravity: Newtonian gravity from a fixed central object'
          if (lvector_gravity) then
            if (t >= tstart_rpbeta) then
              if (lparticles_radius .and. lparticles_radius_rpbeta) then
                rpbeta_tmp_arr(1:npar_loc) = fp(1:npar_loc,irpbeta)
              elseif (npar_species > 1) then
                jspec_arr(1:npar_loc) = assign_species(ipar(1:npar_loc))
              else
                rpbeta_tmp_arr(1:npar_loc) = rpbeta
              endif
            else
              rpbeta_tmp_arr(1:npar_loc) = 0.0
            endif
            if (lcartesian_coords) then
              if (lcylindrical_gravity_par) then
                rr_arr(1:npar_loc) = sqrt(fp(1:npar_loc,ixp)**2+fp(1:npar_loc,iyp)**2+gravsmooth2)
              else
                rr_arr(1:npar_loc) = sqrt(fp(1:npar_loc,ixp)**2+fp(1:npar_loc,iyp)**2+fp(1:npar_loc,izp)**2+gravsmooth2)
              endif
              OO2_arr(1:npar_loc) = rr_arr(1:npar_loc)**(-3.)*gravr*(1.0-rpbeta_tmp_arr(1:npar_loc))
              gpp_arr(1:npar_loc,1) = -fp(1:npar_loc,ixp)*OO2_arr(1:npar_loc)
              gpp_arr(1:npar_loc,2) = -fp(1:npar_loc,iyp)*OO2_arr(1:npar_loc)
              if (lcylindrical_gravity_par) then
                gpp_arr(1:npar_loc,3) = 0.
              else
                gpp_arr(1:npar_loc,3) = -fp(1:npar_loc,izp)*OO2_arr(1:npar_loc)
              endif
              dfp(1:npar_loc,ivpx:ivpz) = dfp(1:npar_loc,ivpx:ivpz) + gpp_arr(1:npar_loc,:)
            elseif (lcylindrical_coords) then
              if (lcylindrical_gravity_par) then
                rr_arr(1:npar_loc) = sqrt(fp(1:npar_loc,ixp)**2+gravsmooth2)
              else
                rr_arr(1:npar_loc) = sqrt(fp(1:npar_loc,ixp)**2+fp(1:npar_loc,izp)**2+gravsmooth2)
              endif
              OO2_arr(1:npar_loc) = rr_arr(1:npar_loc)**(-3.)*gravr*(1.0-rpbeta_tmp_arr(1:npar_loc))
              gpp_arr(1:npar_loc,1) = -fp(1:npar_loc,ixp)*OO2_arr(1:npar_loc)
              gpp_arr(1:npar_loc,2) = 0.0
              if (lcylindrical_gravity_par) then
                gpp_arr(1:npar_loc,3) = 0.
              else
                gpp_arr(1:npar_loc,3) = -fp(1:npar_loc,izp)*OO2_arr(1:npar_loc)
              endif
              dfp(1:npar_loc,ivpx:ivpz) = dfp(1:npar_loc,ivpx:ivpz) + gpp_arr(1:npar_loc,:)
            elseif (lspherical_coords) then
              rr_arr(1:npar_loc) = sqrt(fp(1:npar_loc,ixp)**2+gravsmooth2)
              OO2_arr(1:npar_loc) = rr_arr(1:npar_loc)**(-3)*gravr*(1.0-rpbeta_tmp_arr(1:npar_loc))
              gpp_arr(1:npar_loc,1) = -fp(1:npar_loc,ixp)*OO2_arr(1:npar_loc)
              gpp_arr(1:npar_loc,2) = 0.0
              gpp_arr(1:npar_loc,3) = 0.0
              if (lcylindrical_gravity_par) call fatal_error("particle_gravity", &
                  "no cylindrical gravity in spherical coordinates")
              dfp(1:npar_loc,ivpx:ivpz) = dfp(1:npar_loc,ivpx:ivpz) + gpp_arr(1:npar_loc,:)
            endif
!  Limit time-step if particles close to gravity source.
            if (ldt_grav_par .and.(lupdate_courant_dt)) then
              if (lcartesian_coords) then
                vv_arr(1:npar_loc) = sqrt(fp(1:npar_loc,ivpx)**2+fp(1:npar_loc,ivpy)**2+fp(1:npar_loc,ivpz)**2)
              elseif (lcylindrical_coords) then
                vv_arr(1:npar_loc) = sqrt(fp(1:npar_loc,ivpx)**2+fp(1:npar_loc,ivpz)**2)
              elseif (lspherical_coords) then
                vv_arr(1:npar_loc) = abs(fp(1:npar_loc,ivpx))
              endif
              do k = 1,npar_loc
                dt1_max(ineargrid(k,1)-nghost) = max(dt1_max(ineargrid(k,1)-nghost),vv_arr(k)/rr_arr(k)/cdtpgrav)
              enddo
            endif
          else
            do k = 1,npar_loc
              if (t >= tstart_rpbeta) then
                if (npar_species > 1) then
                  jspec = assign_species(ipar(k))
                  rpbeta_tmp = rpbeta_species(jspec)
                else
                  rpbeta_tmp = rpbeta
                endif
              else
                rpbeta_tmp = 0.0
              endif
              if (lcartesian_coords) then
                if (lcylindrical_gravity_par) then
                  rr = sqrt(fp(k,ixp)**2+fp(k,iyp)**2+gravsmooth2)
                else
                  rr = sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2+gravsmooth2)
                endif
                OO2 = rr**(-3.)*gravr*(1.0-rpbeta_tmp)
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
                  rr = sqrt(fp(k,ixp)**2+gravsmooth2)
                else
                  rr = sqrt(fp(k,ixp)**2+fp(k,izp)**2+gravsmooth2)
                endif
                OO2 = rr**(-3.)*gravr*(1.0-rpbeta_tmp)
                ggp(1) = -fp(k,ixp)*OO2
                ggp(2) = 0.0
                if (lcylindrical_gravity_par) then
                  ggp(3) = 0.
                else
                  ggp(3) = -fp(k,izp)*OO2
                endif
                dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
              elseif (lspherical_coords) then
                rr = sqrt(fp(k,ixp)**2+gravsmooth2)
                OO2 = rr**(-3)*gravr*(1.0-rpbeta_tmp)
                ggp(1) = -fp(k,ixp)*OO2
                ggp(2) = 0.0
                ggp(3) = 0.0
                if (lcylindrical_gravity_par) call fatal_error("particle_gravity", &
                    "no cylindrical gravity in spherical coordinates.")
                dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
              endif
!
!  Limit time-step if particles close to gravity source.
!
              if (ldt_grav_par .and.(lupdate_courant_dt)) then
                if (lcartesian_coords) then
                  vv = sqrt(fp(k,ivpx)**2+fp(k,ivpy)**2+fp(k,ivpz)**2)
                elseif (lcylindrical_coords) then
                  vv = sqrt(fp(k,ivpx)**2+fp(k,ivpz)**2)
                elseif (lspherical_coords) then
                  vv = abs(fp(k,ivpx))
                endif
                dt1_max(ineargrid(k,1)-nghost) = &
                    max(dt1_max(ineargrid(k,1)-nghost),vv/rr/cdtpgrav)
              endif
            enddo
          endif
!
        case default
          call fatal_error('particle_gravity','no such gravr_profile: '//trim(gravr_profile))
        endselect
!
      endif
      if (lfirstcall) lfirstcall = .false.
!
    endsubroutine particle_gravity
!***********************************************************************
    subroutine indirect_inertial_particles(f,df,fp,dfp,ineargrid)
!
!  Clone of gas routine in gravity_r, but for particles.
!
!  24-aug-18/wlad: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mpar_loc,mpvar), intent(inout) :: dfp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      real :: c2, s2, g2, rrcyl, rr2, gp
      real :: sint, cost, sinp, cosp
      logical :: lcentrifugal_force_gravity=.true.,lcoriolis_force_gravity=.true.,lindirect_terms=.true.
      real, dimension(3) :: ggp
      integer :: k
!
      if (lramp_mass .and.(t <= t_ramp_mass)) then
!
!  Ramp up g1 for the first 5 orbits, to prevent to big an initial impulse
!
        g2 = g1/rp1**2 * t/t_ramp_mass
      else
        g2 = g1/rp1**2
      endif
!
!  Do not allow secondary before t_start_secondary
!
      if (lsecondary_wait .and. t <= t_start_secondary) g2 = 0.
!
      do k = 1,npar_loc
!
!  Trigonometric shortcuts
!
        if (lcylindrical_coords) then
          cosp = cos(fp(k,iyp))
          sinp = sin(fp(k,iyp))
        elseif (lspherical_coords) then
          sint = sin(fp(k,iyp))
          cost = cos(fp(k,iyp))
          sinp = sin(fp(k,izp))
          cosp = cos(fp(k,izp))
        endif
!
!  Secondary gravity
!
        if (lcylindrical_coords) then
          if (lcylindrical_gravity) then
            rr2 = fp(k,ixp)**2 + rp1**2 -2*fp(k,ixp)*rp1*cosp
          else
            rr2 = fp(k,ixp)**2 + rp1**2 -2*fp(k,ixp)*rp1*cosp + fp(k,izp)**2
          endif
        elseif (lspherical_coords) then
          rr2 = fp(k,ixp)**2 + rp1**2 - 2*fp(k,ixp)*rp1*sint*cosp
        else
          call not_implemented("secondary_body_gravity","for Cartesian")
        endif
        gp = -g1*(rr2+rp1_smooth**2)**(-1.5)
!
        if (lcylindrical_coords) then
          ggp(1) =  gp * (fp(k,ixp)-rp1*cosp)
          ggp(2) =  gp *            rp1*sinp
          ggp(3) =  gp *  fp(k,izp)
          if (lcylindrical_gravity) ggp(3) = 0.
        elseif (lspherical_coords) then
          ggp(1) =  gp * (fp(k,ixp) - rp1*sint*cosp)
          ggp(2) = -gp *              rp1*cost*cosp
          ggp(3) =  gp *              rp1*     sinp
        endif
        dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + ggp
!
! Corrections coordinates
!
        if (lcylindrical_coords) then
!
!  Indirect terms
!
          if (lindirect_terms) then
            dfp(k,ivpx) = dfp(k,ivpx) - g2*cosp
            dfp(k,ivpy) = dfp(k,ivpy) + g2*sinp
          endif
!
!  Coriolis force
!
          if (lcoriolis_force_gravity) then
            c2 = 2*Omega_corot
            dfp(k,ivpx) = dfp(k,ivpx) + c2*fp(k,ivpy)
            dfp(k,ivpy) = dfp(k,ivpy) - c2*fp(k,ivpx)
          endif
!
!  Centrifugal force
!
          if (lcentrifugal_force_gravity) dfp(k,ivpx) = dfp(k,ivpx) + fp(k,ixp)*Omega_corot**2
!
        elseif (lspherical_coords) then
!
!  Indirect terms
!
          if (lindirect_terms) then
            dfp(k,ivpx) = dfp(k,ivpx) - g2*sint*cosp
            dfp(k,ivpy) = dfp(k,ivpy) - g2*cost*cosp
            dfp(k,ivpz) = dfp(k,ivpz) + g2*     sinp
          endif
!
!  Coriolis force
!
          if (lcoriolis_force_gravity) then
            c2 = 2*Omega_corot*cost
            s2 = 2*Omega_corot*sint
!
            dfp(k,ivpx) = dfp(k,ivpx) + s2*fp(k,ivpz)
            dfp(k,ivpy) = dfp(k,ivpy) + c2*fp(k,ivpz)
            dfp(k,ivpz) = dfp(k,ivpz) - c2*fp(k,ivpy) - s2*fp(k,ivpx)
          endif
!
!  Centrifugal force
!
          if (lcentrifugal_force_gravity) then
            rrcyl = fp(k,ixp)*sint
            dfp(k,ivpx) = dfp(k,ivpx) + rrcyl*sint*Omega_corot**2
            dfp(k,ivpy) = dfp(k,ivpy) + rrcyl*cost*Omega_corot**2
          endif
!
        endif
      enddo
!
    endsubroutine indirect_inertial_particles
!***********************************************************************
    subroutine dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle position (called from main pencil loop).
!
!  25-apr-06/anders: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3) :: ineargrid
!
      integer :: k, ix0, iy0, iz0
      real :: dt1_advpx, dt1_advpy, dt1_advpz
!
!  Contribution of dust particles to time step.
!
      if (lupdate_courant_dt .and. ldt_adv_par) then
        if (npar_imn(imn) /= 0) then
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            iy0 = ineargrid(k,2)
            iz0 = ineargrid(k,3)
            if (nxgrid /= 1) then
              dt1_advpx = abs(fp(k,ivpx))*dx_1(ix0)
            else
              dt1_advpx = 0.0
            endif
            if (nygrid /= 1) then
              if (lshear .and. lcdtp_shear) then
                dt1_advpy = (-qshear*Omega*fp(k,ixp)+abs(fp(k,ivpy)))*dy_1(iy0)
              else
                dt1_advpy = abs(fp(k,ivpy))*dy_1(iy0)
              endif
            else
              dt1_advpy = 0.0
            endif
            if (nzgrid /= 1) then
              dt1_advpz = abs(fp(k,ivpz))*dz_1(iz0)
            else
              dt1_advpz = 0.0
            endif
            if (l_shell) then
              dt1_advpx = abs(fp(k,ivpx))/k_shell
              dt1_advpy = abs(fp(k,ivpy))/k_shell
              dt1_advpz = abs(fp(k,ivpz))/k_shell
            endif
            dt1_max(ix0-nghost) = max(dt1_max(ix0-nghost), &
                sqrt(dt1_advpx**2+dt1_advpy**2+dt1_advpz**2)/cdtp)
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
      use Particles_spin, only: calc_liftforce
      use Particles_diagnos_dv, only: collisions
      use Particles_diagnos_state, only: persistence_check
      use SharedVariables, only: get_shared_variable
      use Viscosity, only: getnu
      use Particles_caustics, only: dcaustics_dt_pencil
      use Particles_tetrad, only: dtetrad_dt_pencil
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case) :: p
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      real, dimension(mpar_loc,mpvar), intent(inout) :: dfp
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
!
      real, dimension(nx) :: dt1_drag_gas, dt1_drag_dust
      real, dimension(3) :: dragforce, liftforce, bforce, thermforce, uup, xxp
      real, dimension(3) :: adv_der_uup = 0.0
      real, dimension(:), allocatable :: rep, stocunn
      real :: added_mass_beta = 0.0
      real :: rho1_point, tausp1_par, up2, urel
      real :: weight, weight_x, weight_y, weight_z
      real :: dxp, dyp, dzp, volume_cell
      integer :: k, l, ix0, iy0, iz0, irad
      integer :: ixx, iyy, izz, ixx0, iyy0, izz0, ixx1, iyy1, izz1
      integer, dimension(3) :: inear
      logical :: lsink, lvapour
      real, pointer :: pscalar_diff, ap0(:)
      real :: Sherwood, mass_trans_coeff, lambda_tilde
      real :: dthetadt, mp_vcell
!
      real, dimension(k1_imn(imn):k2_imn(imn)) :: nu
      character(len=labellen) :: ivis=''
      real :: nu_
      real :: inversetau
!
!  Identify module.
!
      if (headtt) then
        print*,'dvvp_dt_pencil: calculate dvvp_dt'
        print*, 'dvvp_dt_pencil: ldraglaw_purestokes=', ldraglaw_purestokes
        if (lhydro .and. ldragforce_gas_par) then
          print*,'dvvp_dt_pencil: adding feedback from dust to gas'
          if (eps_dtog == 0.) then
            print*,'dvvp_dt_pencil: eps_dtog=',eps_dtog
            print*, 'dvvp_dt_pencil: mp_swarm=', mp_swarm
            if (.not. lparticles_number) then
              print*,'dvvp_dt_pencil: lparticles_number=',lparticles_number
              print*, 'dvvp_dt_pencil: mp_vcell=4.*pi*fp(k,iap)**3*rhopmat/(3.*volume_cell)'
              print*, 'dvvp_dt_pencil: df(ixx,iyy,izz,iux:iuz)= mp_vcell*rho1_point*dragforce*weight'
            endif
          endif
        endif
      endif
!
!  When the field based handling of passive scalar consumption is enabled, set
!  the current pencil of the lncc auxiliary to zero
!
! supersat
!      if (lascalar) f(:,m,n,itauascalar) = 0.0
!
      if (lascalar .and. ltauascalar) f(:,m,n,itauascalar) = 0.0
      if (lascalar) f(:,m,n,icondensationRate) = 0.0
      if (lascalar) f(:,m,n,iwaterMixingRatio) = 0.0
      if (ldiffuse_passive) f(:,m,n,idlncc) = 0.0
      if (ldiffuse_passive .and. ilncc == 0) &
          call fatal_error('particles_dust','ldiffuse_passive needs pscalar_nolog=F')
!
!  Initialize the pencils that are calculated within this subroutine
!
      if (lpenc_requested(i_npvz))     p%npvz = 0.
      if (lpenc_requested(i_npvz2))    p%npvz2 = 0.
      if (lpenc_requested(i_npuz))     p%npuz = 0.
      if (lpenc_requested(i_np_rad))   p%np_rad = 0.
      if (lpenc_requested(i_sherwood)) p%sherwood = 0.
      if (lpenc_requested(i_tauascalar)) p%tauascalar = 0.
      if (lpenc_requested(i_condensationRate)) p%condensationRate = 0.
      if (lpenc_requested(i_waterMixingRatio)) p%waterMixingRatio = 0.
!
      if (idiag_urel /= 0) urel_sum = 0.
!
!  If lfollow_gas=T, the evolution equation for the particle velocity is not
!  calculated. Instead we set the particle velocity equal to the velocity of
!  the gas. Brownian velocities might be added on top of that. This is useful
!  for particles with very short response times.
!
      if (lfollow_gas) then
        !
        ! Since we set the velocity directly in the fp-array here, we
        ! should do this only at the first sub-timestep.
        !
        if (lfirst) then
          if (npar_imn(imn) /= 0) then
            if (lbrownian_forces_Li_Ahmadi) then
              allocate(stocunn(k1_imn(imn):k2_imn(imn)))
              call calc_stokes_cunningham(fp,stocunn)
            endif
            do k = k1_imn(imn),k2_imn(imn)
              fp(k,ivpx:ivpz) = interp_uu(k,:)
              !
              ! Add Brownian motion on top of the fluid velocity
              !
              if (lbrownian_forces) then
                if (lbrownian_forces_Li_Ahmadi) then
                  call calc_brownian_force_Li_Ahmadi(fp,k,ineargrid(k,:),stocunn(k),bforce)
                else
                  call calc_brownian_force(f,fp,k,ineargrid(k,:),bforce,tausp1_par)
                endif
                fp(k,ivpx:ivpz) = fp(k,ivpx:ivpz) + bforce/tausp1_par
              endif
            enddo
            if (lbrownian_forces_Li_Ahmadi) then
              if (allocated(stocunn)) deallocate(stocunn)
            endif
          endif
        endif
      else
!
!  Precalculate certain quantities, if necessary.
!
        if (npar_imn(imn) /= 0) then          
!
!  Precalculate particle Reynolds numbers.
!
        if (ldraglaw_steadystate .or. lparticles_spin .or. ldraglaw_stokesschiller) then
          allocate(rep(k1_imn(imn):k2_imn(imn)))
          if (.not. allocated(rep)) call fatal_error('dvvp_dt_pencil','unable to allocate rep',.true.)
          call calc_pencil_rep(fp, rep)
        endif
!
!  Precalculate Stokes-Cunningham factor (only if not ldraglaw_simple  or ldraglaw_purestokes)
!
        if (lbrownian_forces .or. .not. (ldraglaw_simple .or. ldraglaw_purestokes .or. ldraglaw_stokesschiller)) then
          if (ldraglaw_steadystate .or. lbrownian_forces_Li_Ahmadi) then
            allocate(stocunn(k1_imn(imn):k2_imn(imn)))
            if (.not. allocated(stocunn)) &
              call fatal_error('dvvp_dt_pencil','unable to allocate stocunn',.true.)
!
            call calc_stokes_cunningham(fp,stocunn)
          endif
        endif
      endif
!
!  Drag force on particles and on gas.
!
      if (ldragforce_heat .or. (ldiagnos .and. idiag_dedragp /= 0)) drag_heat = 0.0
!
      if (ldragforce_dust_par .and. t >= tstart_dragforce_par) then
        if (headtt) print*,'dvvp_dt: Add drag force; tausp=', tausp
!
        if (npar_imn(imn) /= 0) then
!
          if (lupdate_courant_dt) then
            dt1_drag_dust = 0.0
            if (ldragforce_gas_par) dt1_drag_gas = 0.0
          endif
!
!  Get viscosity used to calculate the pscalar Schmidt number
!
          if (lpscalar .and. lpscalar_sink) then
            if (.not. lsherwood_const) then
              call getnu(nu_input=nu_,IVIS=ivis)
              if (ivis == 'nu-const') then
                nu = nu_
              elseif (ivis == 'nu-mixture') then
                nu = interp_nu
              elseif (ivis == 'rho-nu-const') then
                nu = nu_/interp_rho(k1_imn(imn):k2_imn(imn))
              elseif (ivis == 'sqrtrho-nu-const') then
                nu = nu_/sqrt(interp_rho(k1_imn(imn):k2_imn(imn)))
              elseif (ivis == 'nu-therm') then
                nu = nu_*sqrt(interp_TT(k1_imn(imn):k2_imn(imn)))
              elseif (ivis == 'mu-therm') then
                nu = nu_*sqrt(interp_TT(k1_imn(imn):k2_imn(imn)))/interp_rho(k1_imn(imn):k2_imn(imn))
              else
                call fatal_error('dvvp_dt_pencil','no such ivis: '//trim(ivis))
              endif
            endif
!
!  Get the passive scalar diffusion rate
!
            call get_shared_variable('pscalar_diff',pscalar_diff,caller='dvvp_dt_pencil')
!
          endif
!
!  Loop over all particles in current pencil.
!         
          do k = k1_imn(imn),k2_imn(imn)
!
!  Vapour particles acquire same speed as the gas.
!
            lvapour = .false.
            if (lparticles_radius) then
              if (fp(k,iap) == 0.0) then
                lvapour = .true.
                inear = ineargrid(k,:)
                xxp = fp(k,ixp:izp)
                if (lparticlemesh_cic) then
                  call interpolate_linear(f,iux,iuz,xxp,uup,inear,0,ipar(k))
                elseif (lparticlemesh_tsc) then
                  if (linterpolate_spline) then
                    call interpolate_quadratic_spline(f,iux,iuz,xxp,uup,inear,0,ipar(k))
                  else
                    call interpolate_quadratic(f,iux,iuz,xxp,uup,inear,0,ipar(k))
                  endif
                else
                  ix0 = ineargrid(k,1)
                  iy0 = ineargrid(k,2)
                  iz0 = ineargrid(k,3)
                  uup = f(ix0,iy0,iz0,iux:iuz)
                endif
                fp(k,ivpx:ivpz) = uup
                dfp(k,ivpx:ivpz) = 0.0
              endif
            endif
!
!  Exclude sink and vapour particles from the drag.
!
            lsink = .false.
            if (lparticles_sink) then
              if (fp(k,iaps) > 0.0) lsink = .true.
            endif
!
            if (.not. lsink .and. .not. lvapour) then
              ix0 = ineargrid(k,1)
              iy0 = ineargrid(k,2)
              iz0 = ineargrid(k,3)
!
!  Calculate required pencils
!  NILS: Could this be moved to calc_pencils_particles
!
              if (lpenc_requested(i_npvz) .or. lpenc_requested(i_npvz2) .or. &
                  lpenc_requested(i_np_rad) .or. lpenc_requested(i_npuz)) then
                call get_shared_variable('ap0',ap0,caller='dvvp_dt_pencil')
                do irad = 1,npart_radii
                  if ((fp(k,iap) > ap0(irad)*0.99) .and. (fp(k,iap) < ap0(irad)*1.01)) then
                    p%npvz(ix0-nghost,irad) = p%npvz(ix0-nghost,irad)+fp(k,ivpz)
                    p%npvz2(ix0-nghost,irad) = p%npvz2(ix0-nghost,irad)+fp(k,ivpz)**2
                    p%np_rad(ix0-nghost,irad) = p%np_rad(ix0-nghost,irad)+1.
                    if (lpenc_requested(i_npuz)) then
!
! DM : I understand that the following is not efficient, it interplolates twice, but seems
! DM: to be the quickest solution now.
!
                      inear = ineargrid(k,:)
                      xxp = fp(k,ixp:izp)
                      if (lparticlemesh_cic) then
                        call interpolate_linear(f,iux,iuz,xxp,uup,inear,0,ipar(k))
                      else
                        if (linterpolate_spline) then
                          call interpolate_quadratic_spline(f,iux,iuz,xxp,uup,inear,0,ipar(k))
                        else
                          call interpolate_quadratic(f,iux,iuz,xxp,uup,inear,0,ipar(k))
                        endif
                      endif
                      p%npuz(ix0-nghost,irad) = p%npuz(ix0-nghost,irad)+uup(3)
                    endif
                  endif
                enddo
              endif
!
!  The interpolated gas velocity is either precalculated, and stored in
!  interp_uu, or it must be calculated here.
!
              if (.not. interp%luu) then
                if (lhydro) then
                  inear = ineargrid(k,:)
                  xxp = fp(k,ixp:izp)
                  if (lparticlemesh_cic) then
                    call interpolate_linear(f,iux,iuz,xxp,uup,inear,0,ipar(k))
                  elseif (lparticlemesh_tsc) then
                    if (linterpolate_spline) then
                      call interpolate_quadratic_spline(f,iux,iuz,xxp,uup,inear,0,ipar(k))
                    else
                      call interpolate_quadratic(f,iux,iuz,xxp,uup,inear,0,ipar(k))
                    endif
                  else
                    uup = f(ix0,iy0,iz0,iux:iuz)
                  endif
                elseif (l_shell) then
                  call calc_gas_velocity_shell_call(k,uup,fp)
                else
                  uup = 0.0
                endif
              else
                uup = interp_uu(k,:)
              endif
!
!  Track particle state in terms of local gas velocity
!
              if (lparticles_diagnos_state .and. lfirst) call persistence_check(fp, k, uup)
!
!  Get the friction time. For the case of |uup| ~> cs, the Epstein drag law
!  is dependent on the relative mach number, hence the need to feed uup as
!  an optional argument to get_frictiontime.
!
              if (ldraglaw_epstein_transonic .or. ldraglaw_eps_stk_transonic) then
                call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par,uup)
              elseif (ldraglaw_steadystate) then
                call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par,rep=rep(k),stocunn=stocunn(k))
              elseif (ldraglaw_stokesschiller) then
                call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par,rep=rep(k))
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
                dragforce = fp(k,ixp:izp)*(dragforce(1)*fp(k,ixp)+ &
                    dragforce(2)*fp(k,iyp)+dragforce(3)*fp(k,izp))/ &
                    (fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2)
              endif
!
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + dragforce
!
! If we are using the module particles_caustics, then we call for them here:
!
              if (lparticles_caustics) call dcaustics_dt_pencil(f,df,fp,dfp,p,ineargrid,k,tausp1_par)
!
! If we are using the module particles_tetrad, then we call for them here:
!
              if (lparticles_tetrad) call dtetrad_dt_pencil(f,df,fp,dfp,p,ineargrid,k,tausp1_par)
!
!  Account for added mass term beta
!  JONAS: The advective derivative of velocity is interpolated for each
!  particle
!
              if (lbubble) then
                if (lhydro) then
                  inear = ineargrid(k,:)
                  xxp = fp(k,ixp:izp)
                  if (lparticlemesh_cic) then
                    call interpolate_linear(f,i_adv_derx,i_adv_derz,xxp,adv_der_uup,inear,0,ipar(k))
                  elseif (lparticlemesh_tsc) then
                    if (linterpolate_spline) then
                      call interpolate_quadratic_spline(f,i_adv_derx,i_adv_derz,xxp,adv_der_uup,inear,0,ipar(k))
                    else
                      call interpolate_quadratic(f,i_adv_derx,i_adv_derz,xxp,adv_der_uup,inear,0,ipar(k))
                    endif
                  else
                    adv_der_uup = f(ix0,iy0,iz0,i_adv_derx:i_adv_derz)
                  endif
                else
                  adv_der_uup = 0.0
                endif
!
!  Calculate the beta for the current particle
!
                call calc_added_mass_beta(fp,k,added_mass_beta)
!
!  Add the contribution of the added mass/virtual mass term to the velocity evolution
!
                dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) + added_mass_beta * adv_der_uup
              endif
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
              if (ldragforce_gas_par .or. (lpscalar_sink .and. lpscalar)) then
!
!  Check if the particles consume passive scalar, and calculate the
!  consumption rate
!
                if (lpscalar_sink .and. lpscalar) then
!
!  JONAS: michaelides 2006,p122
!  Nu/Sh = 0.922+Pe**0.33+0.1*Pe**0.33*Re**0.33
!  Present: rep(k), needed: Pe(k)
!  From Multiphase flows with Droplets and particles, p.62:
!  Sh = 2+0.69*Re_rel**0.5 * Sc**0.33
!  The long number is 0.7**(1/3)
!  Sc = nu/pscalar_diff implement with getnu
!
                  if (lsherwood_const) then
                    Sherwood = 2.
                  else
                    Sherwood = 2.0 + 0.69*sqrt(rep(k))*(nu(k)/pscalar_diff)**(1./3.)
                  endif
!
                  if (lpenc_requested(i_sherwood)) then
                    p%sherwood(ix0-nghost) = p%sherwood(ix0-nghost)+Sherwood
                  endif
                  if (idiag_urel /= 0) then
                    urel = sqrt(sum((interp_uu(k,:) - fp(k,ivpx:ivpz))**2))
                    urel_sum(ix0-nghost) = urel_sum(ix0-nghost)+urel
                  endif
!
                  mass_trans_coeff = Sherwood*pscalar_diff/ (2*fp(k,iap))
                  lambda_tilde = pscalar_sink_rate*mass_trans_coeff/(pscalar_sink_rate+mass_trans_coeff)
                  dthetadt = lambda_tilde*4.*pi*fp(k,iap)**2
                endif
!
!  Cloud In Cell (CIC) scheme.
!
                if (lparticlemesh_cic) then
!
                  ixx0 = ix0
                  iyy0 = iy0
                  izz0 = iz0
                  ixx1 = ix0
                  iyy1 = iy0
                  izz1 = iz0
!
!  Particle influences param_iothe 8 surrounding grid points. The reference point is
!  the grid point at the lower left corner.
!
                  if ( (x(ix0) > fp(k,ixp)) .and. nxgrid /= 1) ixx0 = ixx0-1
                  if ( (y(iy0) > fp(k,iyp)) .and. nygrid /= 1) iyy0 = iyy0-1
                  if ( (z(iz0) > fp(k,izp)) .and. nzgrid /= 1) izz0 = izz0-1
                  if (nxgrid /= 1) ixx1 = ixx0+1
                  if (nygrid /= 1) iyy1 = iyy0+1
                  if (nzgrid /= 1) izz1 = izz0+1
!
                  do izz = izz0,izz1
                    do iyy = iyy0,iyy1
                      do ixx = ixx0,ixx1
                        weight = 1.0
                        if (nxgrid /= 1) weight = weight*( 1.0-abs(fp(k,ixp)-x(ixx))*dx_1(ixx) )
                        if (nygrid /= 1) weight = weight*( 1.0-abs(fp(k,iyp)-y(iyy))*dy_1(iyy) )
                        if (nzgrid /= 1) weight = weight*( 1.0-abs(fp(k,izp)-z(izz))*dz_1(izz) )
!  Save the calculation of rho1 when inside pencil.
                        if ( (iyy /= m).or.(izz /= n).or.(ixx < l1).or.(ixx > l2) ) then
                          rho1_point = 1.0 / get_gas_density(f,ixx,iyy,izz)
                        else
                          rho1_point = p%rho1(ixx-nghost)
                        endif
!
!  Add friction force to grid point.
!NILS: The grid volume should be put into a pencil when required
!                    if ((lpscalar_sink .and. lpscalar) .or. &
!                        (ldragforce_gas_par .and. ldraglaw_steadystate)) &
!DM : the above if condition is always true.
!DM : the following call is being made twice.
!                        call find_grid_volume(ixx,iyy,izz,volume_cell)
! alexrichert: above call to find_grid_volume is superfluous, not sure why
! conditions are different from call below. Perhaps lparticles_radius or
! iap>0 would be a better condition than eps_dtog/ldraglaw_steadystate?
!
                        if (lhydro .and. ldragforce_gas_par) then
                          if ((eps_dtog == 0.) .or. ldraglaw_steadystate) then
                            call find_grid_volume(ixx,iyy,izz,volume_cell)
                            mp_vcell = 4.*pi*fp(k,iap)**3*rhopmat/(3.*volume_cell)
                            if (lparticles_number) then
                              mp_vcell = mp_vcell*fp(k,inpswarm)
                            elseif (np_swarm  >  0) then
                              mp_vcell = mp_vcell*np_swarm
                            endif
                          else
                            call get_rhopswarm(mp_swarm,fp,k,ixx,iyy,izz,mp_vcell)
                          endif
                          if (.not. ldiffuse_dragf) then
                            df(ixx,iyy,izz,iux:iuz) = df(ixx,iyy,izz,iux:iuz) - &
                                mp_vcell*rho1_point*dragforce*weight
                          endif
                        endif
!
                        if (lpscalar_sink .and. lpscalar) then
                          if (.not. ldiffuse_passive) then
                            if (ilncc == 0) then
                              call fatal_error('dvvp_dt_pencil','lpscalar_sink not allowed for pscalar_nolog')
                            else
                              call find_grid_volume(ixx,iyy,izz,volume_cell)
                              df(ixx,iyy,izz,ilncc) = df(ixx,iyy,izz,ilncc) - weight*dthetadt/volume_cell
                            endif
                          endif
                        endif
                      enddo
                    enddo
                  enddo
!
!  Triangular Shaped Cloud (TSC) scheme.
!
                elseif (lparticlemesh_tsc) then
                  if (.not. lparticlemesh_pqs_assignment) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
                    if (nxgrid /= 1) then
                      ixx0 = ix0-1
                      ixx1 = ix0+1
                    else
                      ixx0 = ix0
                      ixx1 = ix0
                    endif
                    if (nygrid /= 1) then
                      iyy0 = iy0-1
                      iyy1 = iy0+1
                    else
                      iyy0 = iy0
                      iyy1 = iy0
                    endif
                    if (nzgrid /= 1) then
                      izz0 = iz0-1
                      izz1 = iz0+1
                    else
                      izz0 = iz0
                      izz1 = iz0
                    endif
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
                    do izz = izz0,izz1
                      do iyy = iyy0,iyy1
                        do ixx = ixx0,ixx1
                          if ( ((ixx-ix0) == -1) .or. ((ixx-ix0) == +1) ) then
                            weight_x = 1.125-1.5* abs(fp(k,ixp)-x(ixx))*dx_1(ixx) + &
                                0.5*(abs(fp(k,ixp)-x(ixx))*dx_1(ixx))**2
                          else
                            if (nxgrid /= 1) &
                                weight_x = 0.75-((fp(k,ixp)-x(ixx))*dx_1(ixx))**2
                          endif
                          if ( ((iyy-iy0) == -1) .or. ((iyy-iy0) == +1) ) then
                            weight_y = 1.125-1.5* abs(fp(k,iyp)-y(iyy))*dy_1(iyy) + &
                                0.5*(abs(fp(k,iyp)-y(iyy))*dy_1(iyy))**2
                          else
                            if (nygrid /= 1) &
                                weight_y = 0.75-((fp(k,iyp)-y(iyy))*dy_1(iyy))**2
                          endif
                          if ( ((izz-iz0) == -1) .or. ((izz-iz0) == +1) ) then
                            weight_z = 1.125-1.5* abs(fp(k,izp)-z(izz))*dz_1(izz) + &
                                0.5*(abs(fp(k,izp)-z(izz))*dz_1(izz))**2
                          else
                            if (nzgrid /= 1) &
                                weight_z = 0.75-((fp(k,izp)-z(izz))*dz_1(izz))**2
                          endif
!
                          weight = 1.0
!
                          if (nxgrid /= 1) weight = weight*weight_x
                          if (nygrid /= 1) weight = weight*weight_y
                          if (nzgrid /= 1) weight = weight*weight_z
!  Save the calculation of rho1 when inside pencil.
                          if ( (iyy /= m).or.(izz /= n).or.(ixx < l1).or.(ixx > l2) ) then
                            rho1_point = 1.0 / get_gas_density(f,ixx,iyy,izz)
                          else
                            rho1_point = p%rho1(ixx-nghost)
                          endif
!NILS: The grid volume should be put into a pencil when required
                          if ((lpscalar_sink .and. lpscalar) .or. &
                              (ldragforce_gas_par .and. ldraglaw_steadystate)) &
                              call find_grid_volume(ixx,iyy,izz,volume_cell)
!  Add friction force to grid point.
                          if (lhydro .and. ldragforce_gas_par) then
!  Calculate the particle mass divided by the cell volume
                            if ((eps_dtog == 0.) .or. ldraglaw_steadystate) then
                              call find_grid_volume(ixx,iyy,izz,volume_cell)
                              mp_vcell = 4.*pi*fp(k,iap)**3*rhopmat/(3.*volume_cell)
                              if (lparticles_number) then
                                mp_vcell = mp_vcell*fp(k,inpswarm)
                              elseif (np_swarm  >  0) then
                                mp_vcell = mp_vcell*np_swarm
                              endif
                            else
                              call get_rhopswarm(mp_swarm,fp,k,ixx,iyy,izz,mp_vcell)
                            endif
                            if (.not. lcompensate_sedimentation) then
                              if (.not. ldiffuse_dragf) then
                                df(ixx,iyy,izz,iux:iuz) = df(ixx,iyy,izz,iux:iuz) - &
                                    mp_vcell*rho1_point*dragforce*weight
                              endif
                            else
                              df(ixx,iyy,izz,iux:iuz) = df(ixx,iyy,izz,iux:iuz) - &
                                  mp_vcell*rho1_point*dragforce*weight*compensate_sedimentation
                            endif
                          endif
                          if (lpscalar_sink .and. lpscalar) then
                            if (.not. ldiffuse_passive) then
                              if (ilncc == 0) then
                                call fatal_error('dvvp_dt_pencil', &
                                                 'lpscalar_sink not allowed for pscalar_nolog!')
                              else
                                call find_grid_volume(ixx,iyy,izz,volume_cell)
                                df(ixx,iyy,izz,ilncc) = df(ixx,iyy,izz,ilncc)-weight*dthetadt/volume_cell
                              endif
                            endif
                          endif
                        enddo
                      enddo
                    enddo
                  else
!
!  QPS assignment, one order higher than TSC. Experimental.
!
                    if (nxgrid /= 1) then
                      ixx0 = ix0-2
                      ixx1 = ix0+2
                    else
                      ixx0 = ix0
                      ixx1 = ix0
                    endif
                    if (nygrid /= 1) then
                      iyy0 = iy0-2
                      iyy1 = iy0+2
                    else
                      iyy0 = iy0
                      iyy1 = iy0
                    endif
                    if (nzgrid /= 1) then
                      izz0 = iz0-2
                      izz1 = iz0+2
                    else
                      izz0 = iz0
                      izz1 = iz0
                    endif
                    do izz = izz0,izz1
                      do iyy = iyy0,iyy1
                        do ixx = ixx0,ixx1
                          dxp = (x(ixx)-fp(k,ixp))*dx_1(ixx)
                          dyp = (y(iyy)-fp(k,iyp))*dy_1(iyy)
                          dzp = (z(izz)-fp(k,izp))*dz_1(izz)
!
                          if (abs(dxp) <= 1.0) then
                            weight_x = (2/3.0)+dxp**2*(abs(dxp)/2-1.0)
                          elseif (abs(dxp) <= 2.0) then
                            weight_x = (1/6.0)*(2.0-abs(dxp))**3
                          else
                            weight_x = 0.0
                          endif
!
                          if (abs(dyp) <= 1.0) then
                            weight_y = (2/3.0)+dyp**2*(abs(dyp)/2-1.0)
                          elseif (abs(dyp) <= 2.0) then
                            weight_y = (1/6.0)*(2.0-abs(dyp))**3
                          else
                            weight_y = 0.0
                          endif
!
                          if (abs(dzp) <= 1.0) then
                            weight_z = (2/3.0)+dzp**2*(abs(dzp)/2-1.0)
                          elseif (abs(dzp) <= 2.0) then
                            weight_z = (1/6.0)*(2.0-abs(dzp))**3
                          else
                            weight_z = 0.0
                          endif
!
                          weight = 1.0
                          if (nxgrid /= 1) weight = weight*weight_x
                          if (nygrid /= 1) weight = weight*weight_y
                          if (nzgrid /= 1) weight = weight*weight_z
!
                          if ( (iyy /= m).or.(izz /= n).or.(ixx < l1).or.(ixx > l2) ) then
                            rho1_point = 1.0 / get_gas_density(f,ixx,iyy,izz)
                          else
                            rho1_point = p%rho1(ixx-nghost)
                          endif
!
!NILS: The grid volume should be put into a pencil when required
                          if (ldragforce_gas_par .and. ldraglaw_steadystate) then
                            call find_grid_volume(ixx,iyy,izz,volume_cell)
                          endif
                          if (lhydro .and. ldragforce_gas_par) then
!  Calculate the particle mass divided by the cell volume
                            if ((eps_dtog == 0.) .or. ldraglaw_steadystate) then
                              call find_grid_volume(ixx,iyy,izz,volume_cell)
                              mp_vcell = 4.*pi*fp(k,iap)**3*rhopmat/(3.*volume_cell)
                              if (lparticles_number) then
                                mp_vcell = mp_vcell*fp(k,inpswarm)
                              elseif (np_swarm  >  0) then
                                mp_vcell = mp_vcell*np_swarm
                              endif
                            else
                              call get_rhopswarm(mp_swarm,fp,k,ixx,iyy,izz, mp_vcell)
                            endif
                            if (.not. lcompensate_sedimentation) then
                              df(ixx,iyy,izz,iux:iuz) = df(ixx,iyy,izz,iux:iuz) - &
                                  mp_vcell*rho1_point*dragforce*weight
                            else
                              df(ixx,iyy,izz,iux:iuz) = df(ixx,iyy,izz,iux:iuz) - &
                                  mp_vcell*rho1_point*dragforce*weight*compensate_sedimentation
                            endif
                          endif
                          if (lpscalar_sink .and. lpscalar) then
                            if (.not. ldiffuse_passive) then
                              if (ilncc == 0) then
                                call fatal_error('dvvp_dt_pencil', &
                                                 'lpscalar_sink not allowed for pscalar_nolog')
                              else
                                call find_grid_volume(ixx,iyy,izz,volume_cell)
                                df(ixx,iyy,izz,ilncc)=df(ixx,iyy,izz,ilncc)-weight*dthetadt/volume_cell
                              endif
                            endif
                          endif
                        enddo
                      enddo
                    enddo
                  endif
!
!  GAussian Box (GAB) scheme.
!
                elseif (lparticlemesh_gab) then
!
                  call find_grid_volume(ix0,iy0,iz0,volume_cell)
!
!  In case the diffusion of particle-fluid interaction is activated, the values are only taken from
!  the nearest grid point
!
                  if (ldiffuse_dragf .and. ldiffuse_passive) then
                    if ( (iy0 /= m).or.(iz0 /= n).or.(ix0 < l1).or.(ix0 > l2) ) then
                      rho1_point = 1.0 / get_gas_density(f,ix0,iy0,iz0)
                    else
                      rho1_point = p%rho1(ix0-nghost)
                    endif
                    if (lhydro .and. ldragforce_gas_par) then
                      if ((eps_dtog == 0.) .or. ldraglaw_steadystate) then
                        mp_vcell = 4.*pi*fp(k,iap)**3*rhopmat/(3.*volume_cell)
                      endif
                    else
                      call fatal_error('particles_dust', &
                        'diffusion of dragforce needs ldragforce_gas_par and ldraglaw_steadystate')
                    endif
                  else
!
                    ixx0 = ix0-3
                    ixx1 = ix0+3
                    iyy0 = iy0-3
                    iyy1 = iy0+3
                    izz0 = iz0-3
                    izz1 = iz0+3
!
                    do izz = izz0,izz1
                      do iyy = iyy0,iyy1
                        do ixx = ixx0,ixx1
                          weight = 1.
                          if (nxgrid /= 1) weight = weight*gab_weights(abs(ixx-ix0)+1)
                          if (nygrid /= 1) weight = weight*gab_weights(abs(iyy-iy0)+1)
                          if (nzgrid /= 1) weight = weight*gab_weights(abs(izz-iz0)+1)
!
                          if ( (iyy /= m).or.(izz /= n).or.(ixx < l1).or.(ixx > l2) ) then
                            rho1_point = 1.0 / get_gas_density(f,ixx,iyy,izz)
                          else
                            rho1_point = p%rho1(ixx-nghost)
                          endif
!  Add friction force to grid point.
! alexrichert: above call to find_grid_volume is superfluous, not sure why
! conditions are different from call below. Perhaps lparticles_radius or
! iap>0 would be a better condition than eps_dtog/ldraglaw_steadystate?
                          if (lhydro .and. ldragforce_gas_par) then
                            if ((eps_dtog == 0.) .or. ldraglaw_steadystate) then
                              mp_vcell = 4.*pi*fp(k,iap)**3*rhopmat/(3.*volume_cell)
                              if (lparticles_number) then
                                mp_vcell = mp_vcell*fp(k,inpswarm)
                              elseif (np_swarm  >  0) then
                                mp_vcell = mp_vcell*np_swarm
                              endif
                            else
                              call get_rhopswarm(mp_swarm,fp,k,ixx,iyy,izz,mp_vcell)
                            endif
                            if (.not. ldiffuse_dragf) then
                              df(ixx,iyy,izz,iux:iuz) = df(ixx,iyy,izz,iux:iuz) - &
                                                        mp_vcell*rho1_point*dragforce*weight
                            endif
                          endif
                          if (lpscalar_sink .and. lpscalar) then
                            if (.not. ldiffuse_passive) then
                              if (ilncc == 0) then
                                call fatal_error('dvvp_dt_pencil','lpscalar_sink not allowed for pscalar_nolog')
                              else
                                df(ixx,iyy,izz,ilncc) = df(ixx,iyy,izz,ilncc) - weight*dthetadt/volume_cell
                              endif
                            endif
                          endif
                        enddo
                      enddo
                    enddo
                  endif
                else
!
!  Nearest Grid Point (NGP) scheme.
!
                  l = ineargrid(k,1)
!NILS: The grid volume should be put into a pencil when required
                  if ((lpscalar_sink .and. lpscalar) .or. &
                      (ldragforce_gas_par .and. ldraglaw_steadystate)) &
                      call find_grid_volume(ix0,iy0,iz0,volume_cell)
                  if (lhydro .and. ldragforce_gas_par) then
!  Calculate the particle mass divided by the cell volume
                    if ((eps_dtog == 0.) .or. ldraglaw_steadystate) then
                      call find_grid_volume(ix0,iy0,iz0,volume_cell)
                      mp_vcell = 4.*pi*fp(k,iap)**3*rhopmat/(3.*volume_cell)
                      if (lparticles_number) then
                        mp_vcell = mp_vcell*fp(k,inpswarm)
                      elseif (np_swarm  >  0) then
                        mp_vcell = mp_vcell*np_swarm
                      endif
                    else
                      call get_rhopswarm(mp_swarm,fp,k,l,m,n,mp_vcell)
                    endif
                    if (.not. lcompensate_sedimentation) then
                      df(l,m,n,iux:iuz) = df(l,m,n,iux:iuz) - mp_vcell*p%rho1(l-nghost)*dragforce
                    else
                      df(l,m,n,iux:iuz) = df(l,m,n,iux:iuz) - &
                                          mp_vcell*p%rho1(l-nghost)*dragforce*compensate_sedimentation
                    endif
                  endif
                  if (lpscalar_sink .and. lpscalar) then
                    if (.not. ldiffuse_passive) then
                      if (ilncc == 0) then
                        call fatal_error('dvvp_dt_pencil','lpscalar_sink not allowed for pscalar_nolog')
                      else
                        call find_grid_volume(l,m,n,volume_cell)
                        df(l,m,n,ilncc) = df(l,m,n,ilncc) - dthetadt/volume_cell
                      endif
                    endif
                  endif
                endif
!
!  Adding the passive scalar contribution to the auxiliary grid
!
                if (ldiffuse_passive .and. lpscalar_sink .and. .not. lpencil_check_at_work) then
                  call find_grid_volume(ix0,iy0,iz0,volume_cell)
                  f(ix0,iy0,iz0,idlncc) = f(ix0,iy0,iz0,idlncc) - dthetadt/volume_cell
                endif
!
!  Adding the particle dragforce contribution to the auxiliary grid
!
                if (ldiffuse_dragf .and. ldragforce_gas_par .and. .not. lpencil_check_at_work) then
                  call find_grid_volume(ix0,iy0,iz0,volume_cell)
                  f(ix0,iy0,iz0,idfx:idfz) = f(ix0,iy0,iz0,idfx:idfz) - mp_vcell*rho1_point*dragforce
                endif
!
              endif
!
!  Calculate particle mass density in grid cell
!
              if ((eps_dtog == 0.) .or. ldraglaw_steadystate) then
                call find_grid_volume(ix0,iy0,iz0,volume_cell)
                mp_vcell = 4.*pi*fp(k,iap)**3*rhopmat/(3.*volume_cell)
                if (lparticles_number) then
                  mp_vcell = mp_vcell*fp(k,inpswarm)
                elseif (np_swarm  >  0) then
                  mp_vcell = mp_vcell*np_swarm
                endif
              else
                call get_rhopswarm(mp_swarm,fp,k,ix0,iy0,iz0,mp_vcell)
              endif
!
!  Heating of gas due to drag force.
!
              if (ldragforce_heat .or. (ldiagnos .and. idiag_dedragp /= 0)) then
                if (ldragforce_gas_par) then
                  up2 = sum((fp(k,ivpx:ivpz)-uup)**2)
                else
                  up2 = sum(fp(k,ivpx:ivpz)*(fp(k,ivpx:ivpz)-uup))
                endif
!
! WL: Check if this is right. drag_heat is of dimension nx, not a scalar
!
                drag_heat(ix0-nghost) = drag_heat(ix0-nghost) + mp_vcell*tausp1_par*up2
              endif
!
!  The minimum friction time of particles in a grid cell sets the local friction
!  time-step when there is only drag force on the dust,
!    dt1_drag = max(1/tausp)
!
!  With drag force on the gas as well, the maximum time-step is set as
!    dt1_drag = Sum_k[eps_k/tau_k]
!
              if (lupdate_courant_dt) then
                dt1_drag_dust(ix0-nghost) = max(dt1_drag_dust(ix0-nghost), tausp1_par)
                if (ldragforce_gas_par) &
                    dt1_drag_gas(ix0-nghost) = dt1_drag_gas(ix0-nghost) + mp_vcell * p%rho1(ix0-nghost) * tausp1_par
              endif
!
!  Particle growth by condensation in a active scalar field,
!  calculate relaxation time.
!  14-June-16/Xiang-Yu: coded
!
              if (lascalar) then
                inversetau = 4.*pi*rhopmat/rhoa*fp(k,iap)*fp(k,inpswarm)
                if (ascalar_ngp) then
                  l = ineargrid(k,1)
                  if (ltauascalar) f(l,m,n,itauascalar) = f(l,m,n,itauascalar) + A3*A2*inversetau
                  if (lcondensation_rate) then
                    f(l,m,n,icondensationRate) = f(l,m,n,icondensationRate)+f(l,m,n,issat)*inversetau*G_condensation
                    f(l,m,n,iwaterMixingRatio) = f(l,m,n,iwaterMixingRatio)+(4.*pi*rhopmat/3./rhoa)*(fp(k,iap))**3*fp(k,inpswarm)
                  endif
                elseif (ascalar_cic) then
                  ixx0 = ix0
                  iyy0 = iy0
                  izz0 = iz0
                  ixx1 = ix0
                  iyy1 = iy0
                  izz1 = iz0
                  if ( (x(ix0) > fp(k,ixp)) .and. nxgrid /= 1) ixx0 = ixx0-1
                  if ( (y(iy0) > fp(k,iyp)) .and. nygrid /= 1) iyy0 = iyy0-1
                  if ( (z(iz0) > fp(k,izp)) .and. nzgrid /= 1) izz0 = izz0-1
                  if (nxgrid /= 1) ixx1 = ixx0+1
                  if (nygrid /= 1) iyy1 = iyy0+1
                  if (nzgrid /= 1) izz1 = izz0+1
                  do izz = izz0,izz1
                    do iyy = iyy0,iyy1
                      do ixx = ixx0,ixx1
                        weight = 1.0
                        if (nxgrid /= 1) weight = weight*( 1.0-abs(fp(k,ixp)-x(ixx))*dx_1(ixx) )
                        if (nygrid /= 1) weight = weight*( 1.0-abs(fp(k,iyp)-y(iyy))*dy_1(iyy) )
                        if (nzgrid /= 1) weight = weight*( 1.0-abs(fp(k,izp)-z(izz))*dz_1(izz) )
                        f(ixx,iyy,izz,itauascalar) = f(ixx,iyy,izz,itauascalar) - weight*inversetau
                      enddo
                    enddo
                  enddo
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
          if (lupdate_courant_dt) then
            if (ldragforce_gas_par) then
              dt1_drag = dt1_drag_dust+dt1_drag_gas
            else
              dt1_drag = dt1_drag_dust
            endif
            dt1_drag = dt1_drag/cdtp_drag
            dt1_max = max(dt1_max,dt1_drag)
          endif
        else
!
!  No particles in this pencil.
!
          if (lupdate_courant_dt) dt1_drag = 0.0
        endif
      endif
!
!  Add friction force from gas that moves systematically slower. Can be used
!  to mimic e.g. a sub-Keplerian gas flow without using the Hydro module.
!
      if (Deltauy_gas_friction /= 0.0 .and. t >= tstart_dragforce_par) then
        if (npar_imn(imn) /= 0) then
          do k = k1_imn(imn),k2_imn(imn)
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
          .and. t >= tstart_collisional_cooling) &
          call collisional_cooling(f,df,fp,dfp,p,ineargrid)
!
!  Compensate for increased friction time by appying extra friction force to
!  particles.
!
      if (lcompensate_friction_increase) call compensate_friction_increase(f,fp,dfp,p,ineargrid)
!
!  Add lift forces.
!
      if (lparticles_spin .and. t >= tstart_liftforce_par) then
        if (npar_imn(imn) /= 0) then
          do k = k1_imn(imn),k2_imn(imn)
            call calc_liftforce(fp(k,:), k, rep(k), liftforce)
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz)+liftforce
          enddo
        endif
      endif
!
!  Add Brownian forces.
!
      if (lbrownian_forces .and. t >= tstart_brownian_par) then
        if (npar_imn(imn) /= 0) then
          do k = k1_imn(imn),k2_imn(imn)
            if (lbrownian_forces_Li_Ahmadi) then
              call calc_brownian_force_Li_Ahmadi(fp,k,ineargrid(k,:),stocunn(k),bforce)
            else
              call calc_brownian_force(f,fp,k,ineargrid(k,:),bforce,tausp1)
            endif
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz)+bforce
          enddo
        endif
      endif
!
!  Add Thermophoretic forces.
!
      if (lthermophoretic_forces) then
        if (npar_imn(imn) /= 0) then
          do k = k1_imn(imn),k2_imn(imn)
            call calc_thermophoretic_force(fp,k,ineargrid(k,:),thermforce)
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz)+thermforce
         enddo
        endif
      endif
!
!  For stiff drag force equations we need to store the forces that are
!  unique to the gas.
!
      if (ldragforce_stiff .and. .not. lpencil_check_at_work) &
          f(l1:l2,m,n,ifgx:ifgz) = p%fpres+p%jxbr+p%fvisc
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
endif
!
      call calc_diagnostics_particles(fp,p,ineargrid)
      
    endsubroutine dvvp_dt_pencil
!***********************************************************************
    subroutine calc_diagnostics_particles(fp,p,ineargrid)

      use Diagnostics

      real, dimension (mpar_loc,mparray) :: fp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid

      integer :: k
!
!  Diagnostic output.
!
      if (ldiagnos) then
        call sum_mn_name(p%np,idiag_npm)
        if (idiag_np2m /= 0) call sum_mn_name(p%np**2,idiag_np2m)
        call max_mn_name(p%np,idiag_npmax)
!
        if (idiag_npmin /= 0) then
          if (lnpmin_exclude_zero) then
            call max_mn_name(-merge(p%np,1.0,p%np /= 0),idiag_npmin,lneg=.true.)
          else
            call max_mn_name(-p%np,idiag_npmin,lneg=.true.)
          endif
        endif
!
        call sum_mn_name(p%rhop,idiag_rhopm)
        if (idiag_rhop2m /= 0 )  call sum_mn_name(p%rhop**2,idiag_rhop2m)
        if (idiag_rhoprms /= 0)  call sum_mn_name(p%rhop**2,idiag_rhoprms,lsqrt=.true.)
        call max_mn_name(p%rhop,idiag_rhopmax)
        if (idiag_rhopmin /= 0)  call max_mn_name(-p%rhop,idiag_rhopmin,lneg=.true.)
        call max_mn_name(p%epsp,idiag_epspmax)
        if (idiag_epspmin /= 0)  call max_mn_name(-p%epsp,idiag_epspmin,lneg=.true.)
        call sum_mn_name(p%epsp,idiag_epspm)
        call sum_mn_name(drag_heat,idiag_dedragp)
        if (idiag_Shm /= 0)      call sum_mn_name(p%sherwood/npar*nwgrid,idiag_Shm)
        if (idiag_urel /= 0)     call sum_mn_name(urel_sum/npar*nwgrid,idiag_urel)
        if (idiag_dvpx2m /= 0 .or. idiag_dvpx2m /= 0 .or. idiag_dvpx2m /= 0 .or. &
            idiag_dvpm  /= 0 .or. idiag_dvpmax /= 0) call calculate_rms_speed(fp,ineargrid,p)
        if (lupdate_courant_dt) call max_mn_name(dt1_drag,idiag_dtdragp,l_dt=.true.)
        if (idiag_latentheatm/= 0) &
             call sum_mn_name(p%latent_heat,idiag_latentheatm)
        if (idiag_ffcondm/= 0) &
             call sum_mn_name(p%ff_cond,idiag_ffcondm)
        if (idiag_ffcondposm/= 0) &
             call sum_mn_name(max(0.,p%ff_cond),idiag_ffcondposm)
        if (idiag_ffcondnegm/= 0) &
             call sum_mn_name(min(0.,p%ff_cond),idiag_ffcondnegm)
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
        call yzsum_mn_name_x(p%np,idiag_npmx)
        call xzsum_mn_name_y(p%np,idiag_npmy)
        call xysum_mn_name_z(p%np,idiag_npmz)
        call yzsum_mn_name_x(p%rhop,idiag_rhopmx)
        call xzsum_mn_name_y(p%rhop,idiag_rhopmy)
        call xysum_mn_name_z(p%rhop,idiag_rhopmz)
        if (idiag_rhop2mx /= 0) call yzsum_mn_name_x(p%rhop**2,idiag_rhop2mx)
        if (idiag_rhop2my /= 0) call xzsum_mn_name_y(p%rhop**2,idiag_rhop2my)
        if (idiag_rhop2mz /= 0) call xysum_mn_name_z(p%rhop**2,idiag_rhop2mz)
        call yzsum_mn_name_x(p%epsp,idiag_epspmx)
        call xzsum_mn_name_y(p%epsp,idiag_epspmy)
        call xysum_mn_name_z(p%epsp,idiag_epspmz)
        call phizsum_mn_name_r(p%rhop,idiag_rhopmr)
!
        if (idiag_rpvpxmx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,1), idiag_rpvpxmx)
        if (idiag_rpvpymx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,2), idiag_rpvpymx)
        if (idiag_rpvpzmx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,3), idiag_rpvpzmx)
        if (idiag_rpvpx2mx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,1)**2, idiag_rpvpx2mx)
        if (idiag_rpvpy2mx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,2)**2, idiag_rpvpy2mx)
        if (idiag_rpvpz2mx /= 0) call yzsum_mn_name_x(p%rhop * p%uup(:,3)**2, idiag_rpvpz2mx)
!
        do k = 1,ndustrad
          call xysum_mn_name_z(p%npvz(:,k),idiag_npvzmz(k))
          call xysum_mn_name_z(p%npvz2(:,k),idiag_npvz2mz(k))
          call xysum_mn_name_z(p%npuz(:,k),idiag_npuzmz(k))
          call xysum_mn_name_z(p%np_rad(:,k),idiag_nptz(k))
        enddo
      endif
!
      if (l2davgfirst) then
        call zsum_mn_name_xy(p%np,idiag_npmxy)
        call phisum_mn_name_rz(p%rhop,idiag_rhopmphi)
        call zsum_mn_name_xy(p%rhop,idiag_rhopmxy)
        call ysum_mn_name_xz(p%rhop,idiag_rhopmxz)
        call zsum_mn_name_xy(p%rhop, idiag_sigmap, lint=.true.)
      endif

    endsubroutine calc_diagnostics_particles
!***********************************************************************
    subroutine dxxp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle position in blocks.
!
!  29-nov-09/anders: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
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
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
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
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(3) :: momp_swarm_removed, momp_swarm_removed_send
      real :: rp, rp_box, rhop_swarm_removed, rhop_swarm_removed_send
      real :: xsinkpar, ysinkpar, zsinkpar
      integer :: k, ksink, iproc_sink, iproc_sink_send
      integer :: ix, ix1, ix2, iy, iy1, iy2, iz, iz1, iz2
      integer, parameter :: itag1=100, itag2=101
      real :: apar
!
!  Sinkparticle activated at time tstart_sink_par
      if (t <= tstart_sink_par) return
!
      if (lsinkpoint) then
        k = 1
        do while (k <= npar_loc)
          rp = sqrt((fp(k,ixp)-xsinkpoint)**2+(fp(k,iyp)-ysinkpoint)**2+(fp(k,izp)-zsinkpoint)**2)
          if (rp < rsinkpoint) then
            call remove_particle(fp,ipar,k,dfp,ineargrid)
          else
            k = k+1
          endif
        enddo
      endif
!
      if (lsinkparticle_1) then
!
!  Search for position of particle with index 1.
!
        xsinkpoint = 0.0
        ysinkpoint = 0.0
        zsinkpoint = 0.0
        iproc_sink = -1
        do k = 1,npar_loc
          if (ipar(k) == 1) then
            ksink = k
            iproc_sink = iproc
            xsinkpar = fp(k,ixp)
            ysinkpar = fp(k,iyp)
            zsinkpar = fp(k,izp)
          endif
        enddo
!
!  Communicate position of sink particle to all processors.
!
        iproc_sink_send = iproc_sink
        call mpireduce_max_int(iproc_sink_send,iproc_sink)
        call mpibcast_int(iproc_sink)
        call mpibcast_real(xsinkpar,proc=iproc_sink)
        call mpibcast_real(ysinkpar,proc=iproc_sink)
        call mpibcast_real(zsinkpar,proc=iproc_sink)
!
!  Remove particles that are too close to sink particle.
!
        if (lparticles_density) then
          rhop_swarm_removed = 0.0
          momp_swarm_removed = 0.0
        endif
!
        k = 1
        do while (k <= npar_loc)
          rp = -1.0
!
!  Take into account periodic boundary conditions.
!
          if (bcpx == 'p' .and. nxgrid /= 1) then
            ix1 = -1
            ix2 = +1
          else
            ix1 = 0
            ix2 = 0
          endif
          if (bcpy == 'p' .and. nygrid /= 1) then
            iy1 = -1
            iy2 = +1
          else
            iy1 = 0
            iy2 = 0
          endif
          if (bcpz == 'p' .and. nzgrid /= 1) then
            iz1 = -1
            iz2 = +1
          else
            iz1 = 0
            iz2 = 0
          endif
          do iz = iz1,iz2
            do iy = iy1,iy2
              do ix = ix1,ix2
                if (lshear) then
                  rp_box = sqrt((fp(k,ixp)+Lx*ix-xsinkpar)**2 + &
                      (fp(k,iyp)+Ly*iy-deltay*ix-ysinkpar)**2 + (fp(k,izp)+Lz*iz-zsinkpar)**2)
                else
                  rp_box = sqrt((fp(k,ixp)+Lx*ix-xsinkpar)**2 + &
                           (fp(k,iyp)+Ly*iy-ysinkpar)**2+(fp(k,izp)+Lz*iz-zsinkpar)**2)
                endif
                if (rp < 0.0) then
                  rp = rp_box
                else
                  rp = min(rp,rp_box)
                endif
              enddo
            enddo
          enddo
          if (ipar(k) /= 1 .and. rp < rsinkparticle_1) then
            if (lparticles_density) then
              rhop_swarm_removed = rhop_swarm_removed + fp(k,irhopswarm)
              momp_swarm_removed = momp_swarm_removed + fp(k,ivpx:ivpz)*fp(k,irhopswarm)
!              if (lshear) then
!                momp_swarm_removed(2) = momp_swarm_removed(2) + &
!                    Sshear*fp(k,ixp)*fp(k,irhopswarm)
!              endif
            endif
            call remove_particle(fp,ipar,k,dfp,ineargrid)
          else
            k = k+1
          endif
        enddo
!
!  Add mass and momentum to sink particle.
!
        if (lparticles_density) then
          rhop_swarm_removed_send = rhop_swarm_removed
          momp_swarm_removed_send = momp_swarm_removed
          call mpireduce_sum(rhop_swarm_removed_send,rhop_swarm_removed)
          call mpireduce_sum(momp_swarm_removed_send,momp_swarm_removed,3)
          if (lroot) then
            call mpisend_real(rhop_swarm_removed,iproc_sink,itag1)
            call mpisend_real(momp_swarm_removed,3,iproc_sink,itag2)
          endif
          if (iproc == iproc_sink) then
            call mpirecv_real(rhop_swarm_removed,0,itag1)
            call mpirecv_real(momp_swarm_removed,3,0,itag2)
!
!  Need to find sink particle again since particle removal may have shifted
!  its index.
!
            do k = 1,npar_loc
              if (ipar(k) == 1) then
                ksink = k
                iproc_sink = iproc
                xsinkpar = fp(k,ixp)
                ysinkpar = fp(k,iyp)
                zsinkpar = fp(k,izp)
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
            fp(ksink,ivpx:ivpz) = (fp(ksink,ivpx:ivpz)*fp(ksink,irhopswarm) + momp_swarm_removed)/ &
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
        k = 1
        do while (k <= npar_loc)
          if (lparticles_radius) then
            apar = fp(k,iap)
          else
            apar = particle_radius
          endif
          if (in_solid_cell(fp(k,ixp:izp),apar)) then
            call remove_particle(fp,ipar,k,dfp,ineargrid)
          else
            k = k+1
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
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine create_particles_sink_simple
!***********************************************************************
    subroutine get_frictiontime(f,fp,p,ineargrid,k,tausp1_par,uup,nochange_opt,rep,stocunn)
!
!  Calculate the friction time.
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mparray) :: fp
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3) :: ineargrid
      integer :: k
      real :: tausp1_par
      real, optional, dimension(3), intent(in) :: uup
      logical, optional :: nochange_opt
      real, optional, intent(in) :: rep
      real, optional :: stocunn
!
      real :: tausg1_point, OO, tmp
      integer :: ix0, iy0, iz0, inx0, jspec
      logical :: nochange=.true.
!
      if (present(nochange_opt)) then
        if (nochange_opt) then
          nochange = .true.
        else
          nochange = .false.
        endif
      else
        nochange = .false.
      endif
!
      ix0 = ineargrid(k,1)
      iy0 = ineargrid(k,2)
      iz0 = ineargrid(k,3)
      inx0 = ix0-nghost
!
!  Epstein drag law.
!
      if (ldraglaw_epstein) then
        if (iap /= 0) then
          if (fp(k,iap) /= 0.0) tausp1_par = (interp_cs(k)*interp_rho(k))/(fp(k,iap)*rhopmat)
!
! DM : 10 Nov  2016
! For the usual Epstein drag we need also the thermal velocity and the
! density at this point. None of them fluctuate by large amounts, so it is better
! to get them from the same pencil without any interpolation. This may be even
! a better prescription for flows with shocks if the particle sits very close to
! a shock.
!
        else
!  Check if we are using multiple or single particle species.
          if (npar_species > 1) then
            jspec = assign_species(ipar(k))
            tmp = tausp1_species(jspec)
          else
            tmp = tausp1
          endif
!
!  Scale friction time with local density.
!
          if (ldraglaw_variable_density) then
            tausp1_par = tmp * get_gas_density(f,ix0,iy0,iz0)
!
!  Discriminate between constant tau and special case for
!  1/tau=omega when omega is not constant (as, for instance,
!  global Keplerian disks, for which omega=rad**(-3/2)
!
          elseif (ldraglaw_variable) then
            if (lcartesian_coords) then
              OO = (fp(k,ixp)**2 + fp(k,iyp)**2)**(-0.75)
            elseif (lcylindrical_coords) then
              OO = fp(k,ixp)**(-1.5)
            elseif (lspherical_coords) then
              OO = (fp(k,ixp)*sin(fp(k,iyp)))**(-1.5)
            else
              call fatal_error("get_frictiontime", "no valid coordinate system")
              OO = 0.
            endif
            tausp1_par = tmp*OO
          else
!
!  Constant friction time.
!
            tausp1_par = tmp
          endif
        endif
      elseif (ldraglaw_epstein_stokes_linear) then
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
        if (iap == 0) &
          call fatal_error('get_frictiontime','need particle radius as dynamical variable for Stokes law')

        if (fp(k,iap) < 2.25*mean_free_path_gas) then
          tausp1_par = 1/(fp(k,iap)*rhopmat)
        else
          tausp1_par = 1/(fp(k,iap)*rhopmat)*2.25*mean_free_path_gas/fp(k,iap)
        endif
!
      elseif (ldraglaw_epstein_transonic) then
!
! Draw laws for intermediate mach number. This is for pure Epstein drag...
!
        call calc_draglaw_parameters(fp,k,uup,p,inx0,tausp1_par)
!
      elseif (ldraglaw_eps_stk_transonic) then
!
! ...and this is for a linear combination of Esptein and Stokes drag at
! intermediate mach number. Pure Stokes drag is not implemented. (implemented now, see below)
!
        call calc_draglaw_parameters(fp,k,uup,p,inx0,tausp1_par,lstokes=.true.)
!
!  draglaw_purestokes implemented below, it is simple stokes drag with no
!  dependence on partile's reynolds number.
!
      elseif (ldraglaw_purestokes) then
        call calc_draglaw_purestokes(fp,k,tausp1_par)
!
!  stokes drag with schiller_nauman, but without epstein
!
      elseif (ldraglaw_stokesschiller) then
        if (.not. present(rep)) then
          call fatal_error('get_frictiontime','need particle reynolds '// &
              'number, rep, to calculate the stokes_schiller drag relaxation time')
        else
          call calc_draglaw_steadystate(fp,k,rep,1.0,tausp1_par)
        endif
      elseif (ldraglaw_steadystate) then
        if (.not. present(rep)) then
          call fatal_error('get_frictiontime','need particle reynolds '// &
              'number, rep, to calculate the steady state drag relaxation time')
        elseif (.not. present(stocunn)) then
          call fatal_error('get_frictiontime','need particle stokes '// &
              'cunningham factor, stocunn, to calculate the steady state drag relaxation time')
        else
          call calc_draglaw_steadystate(fp,k,rep,stocunn,tausp1_par)
        endif
!
! Simple drag law, drag force = 1/\taup (v-u) where \taup is an input parameter.
!
      elseif (ldraglaw_simple) then
!        write(*,*)'DM','simple drag'
!  Check if we are using multiple or single particle species.
        if (npar_species > 1) then
          jspec = assign_species(ipar(k))
          tmp = tausp1_species(jspec)
        else
          tmp = tausp1
        endif
        tausp1_par = tmp
      endif
!
!  Change friction time artificially.
!
      if (.not. nochange) then
!
!  Increase friction time to avoid very small time-steps where the
!  dust-to-gas ratio is high.
!
        if (tausg_min /= 0.0) then
          tausg1_point = tausp1_par*p%epsp(ix0-nghost)
          if (tausg1_point > tausg1_max) tausp1_par = tausg1_max/p%epsp(ix0-nghost)
        endif
!
!  Increase friction time linearly with dust density where the dust-to-gas
!  ratio is higher than a chosen value. Supposed to mimick decreased cooling
!  when the gas follows the dust.
!
        if (epsp_friction_increase /= 0.0) then
          if (p%epsp(ix0-nghost) > epsp_friction_increase) &
              tausp1_par = tausp1_par/(p%epsp(ix0-nghost)/epsp_friction_increase)
        endif
!
      endif

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine get_frictiontime
!***********************************************************************
    subroutine calc_draglaw_parameters(fp,k,uup,p,inx0,tausp1_par,lstokes)
!
      use EquationOfState, only: rho0,cs0
!
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(3) :: uup, duu
      type (pencil_case) :: p
      real :: tausp1_par, tmp, tmp1
      integer :: k, inx0, jspec
      real :: kd, fd, mach, mach2, fac, OO
      real :: knudsen, reynolds, lambda
      real :: inv_particle_radius, kn_crit
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
      if (npar_species == 1) then
        tmp = tausp
        tmp1 = tausp1
      else
        jspec = assign_species(ipar(k))
        tmp = tausp_species(jspec)
        tmp1 = tausp1_species(jspec)
      endif
!
!  Relative velocity
!
      duu = fp(k,ivpx:ivpz)-uup
!
      if (nzgrid == 1) then
!  then omega is needed
        if (ldraglaw_variable) then
          !these omegas assume GM=1
          if (lcartesian_coords) then
            OO = (fp(k,ixp)**2 + fp(k,iyp)**2)**(-0.75)
          elseif (lcylindrical_coords) then
            OO = fp(k,ixp)**(-1.5)
          elseif (lspherical_coords) then
            OO = (fp(k,ixp)*sin(fp(k,iyp)))**(-1.5)
          else
            call fatal_error("calc_draglaw_parameters", "no valid coord system")
            OO = 0.
          endif
        else
          OO = nu_epicycle
        endif
      endif
!
!  Possibility to include the transition from Epstein to Stokes drag
!
      if (present(lstokes)) then
!
        if (headtt) print*, 'calc_draglaw_parameters: Epstein-Stokes transonic drag law'
!
!  The mach number and the correction fd to flows of arbitrary mach number
!
        mach = sqrt((duu(1)**2+duu(2)**2+duu(3)**2)/p%cs2(inx0))
        fd = sqrt(1+(9.0*pi/128)*mach**2)
!
!  For Stokes drag, the mean free path is needed
!
!   lambda = 1/rhog*(mu/sigma_coll)_H2
!
!  were mu_H2 is the mean molecular weight of the hydrogen molecule (3.9e-24 g),
!  and sigma_coll its cross section (2e-15 cm^2).
!  Assume that (mu/sigma_coll) is the input parameter mean_free_path_gas
!
        if (mean_free_path_gas == 0) &
          call fatal_error("calc_draglaw_parameters","You use Stokes drag"// &
              "but you forgot to set 'mean_free_path_gas' in the *.in files")
!
        if (nzgrid == 1) then
          !the sqrt(2pi) factor is inside the mean_free_path_gas constant
          lambda = mean_free_path_gas * sqrt(p%cs2(inx0))*rho0/(p%rho(inx0)*OO*cs0)
        else
          lambda = mean_free_path_gas * rho0/p%rho(inx0)
        endif
!
!  The Knudsen number is the ratio of the mean free path to the particle
!  radius, 2s. To keep consistency with the formulation evolving for radius,
!  tausp1 is C/(s*rhopmat) where C is 2/pi for 2d runs and sqrt(8/pi) for 3D
!  runs (because of the sqrt(2*pi) factor coming from the substitution
!  Sigma=rho/(sqrt(2*pi)*H). 's' is the particle radius
        if (iap /= 0) then
          inv_particle_radius = 1/fp(k,iap)
        else
          if (luse_tau_ap) then
            ! use tausp as the radius (in meter) to make life easier
            inv_particle_radius = tmp1
          else
            if (nzgrid == 1) then
              inv_particle_radius = 0.5*pi*tmp1     !rhopmat=1, particle_radius in meters
            else
              inv_particle_radius = sqrt(pi/8)*tmp1 !rhopmat=1, particle_radius in meters
            endif
          endif
        endif
!
        knudsen = 0.5*lambda*inv_particle_radius
!
!  The Stokes drag depends non-linearly on
!
!    Re = 2*s*rho_g*|delta(v)|/mu_kin
!
        reynolds = 3*sqrt(pi/8)*mach/knudsen
!
!  the Reynolds number of the flow over the particle. It can parameterized by
!
        if (reynolds <= 500) then
          kd = 1.0+0.15*reynolds**0.687
        elseif ((reynolds > 500).and.(reynolds <= 1500)) then
          kd = 3.96e-6*reynolds**2.4
        elseif (reynolds > 1500) then
          kd = 0.11*reynolds
        else
          call fatal_error_local("calc_draglaw_parameters", "'reynolds' seems to be NaN")
          kd = 0.
        endif
!
!  And we finally have the Stokes correction to intermediate Knudsen numbers
!  kn_crit is the critical knudsen number where viscous (low reynolds)
!  Stokes and subsonic Epstein friction are equal (Woitke & Helling, 2003)
!
        kn_crit = 3*knudsen
        fac = kn_crit/(kn_crit+1)**2 * (kn_crit*fd + kd)
!
      else
!
!  Only use Epstein drag
!
        if (lfirstcall) print*,'calc_draglaw_parameters: Epstein transonic drag law'
!
        mach2 = (duu(1)**2+duu(2)**2+duu(3)**2)/p%cs2(inx0)
        fd = sqrt(1+(9.0*pi/128)*mach2)
        fac = fd
!
      endif
!
! Calculate tausp1_par for 2d and 3d cases with and without particle_radius
! as a dynamical variable
! In 1-D (nzgrid=1, for example), we still want to check whether we want
! to excape to do some accretion disc experiment where OO /= 0.
! Otherwise, continue as in other cases.
!
      if (iap /= 0) then
        if (fp(k,iap) /= 0.0) then
          if (nzgrid == 1 .and. OO/=0.) then
            tausp1_par = 2*pi_1*OO* p%rho(inx0)*fac/(fp(k,iap)*rhopmat)
          else
            tausp1_par = sqrt(8*pi_1*p%cs2(inx0))*p%rho(inx0)*fac/(fp(k,iap)*rhopmat)
          endif
        endif
      else
! normalize to make tausp1 not dependent on cs0 or rho0
! bad because it comes at the expense of evil divisions
        if (nzgrid == 1 .and. OO/=0.) then
          if (luse_tau_ap) then
            tausp1_par = tmp1*2*pi_1*OO*p%rho(inx0)*fac/(rho0*rhopmat)
          else
            tausp1_par = tmp1*OO*p%rho(inx0)*fac/ rho0
          endif
        else
          if (luse_tau_ap) then
            tausp1_par = tmp1*sqrt(8*pi_1*p%cs2(inx0))*p%rho(inx0)*fac/(rho0*cs0)
          else
            tausp1_par = tmp1*sqrt(p%cs2(inx0))*p%rho(inx0)*fac/(rho0*cs0)
          endif
        endif
      endif
!
      lfirstcall = .false.
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
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(nx,npar_species,npar_species) :: tau_coll1_species
      real, dimension(nx,3,npar_species) :: vvpm_species
      real, dimension(nx,npar_species) :: np_species, vpm_species
      real, dimension(nx,npar_species) :: tau_coll1_tot
      real, dimension(nx,3) :: vvpm
      real, dimension(nx) :: vpm, tau_coll1, tausp1m, vcoll, coll_heat
      real, dimension(nx) :: rhop_swarm_mn
      real, dimension(3) :: deltavp_vec, vbar_jk
      real :: deltavp, tau_cool1_par, dt1_cool
      real :: tausp1_par, tausp1_parj, tausp1_park, tausp_parj, tausp_park
      real :: tausp_parj3, tausp_park3, rhop_swarm_par
      integer :: j, k, l, ix0
      integer :: ispecies, jspecies
!
      if (ldiagnos .or. lentropy .and. lcollisional_heat) coll_heat = 0.0
!
!  Add collisional cooling of the rms speed.
!
      if (lcollisional_cooling_taucool) then
        if (npar_imn(imn) /= 0) then
          vvpm = 0.0
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            vvpm(ix0-nghost,:) = vvpm(ix0-nghost,:) + fp(k,ivpx:ivpz)
          enddo
          do l = 1,nx
            if (p%np(l) > 1.0) vvpm(l,:) = vvpm(l,:)/p%np(l)
          enddo
!
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - taucool1*(fp(k,ivpx:ivpz)-vvpm(ix0-nghost,:))
          enddo
!
          if (lupdate_courant_dt) dt1_max = max(dt1_max,taucool1/cdtp)
        endif
      endif
!
!  Add collisional cooling of the rms speed, with cooling time-scale based
!  on friction time and local rms speed of particles.
!
      if (lcollisional_cooling_rms) then
        if (npar_imn(imn) /= 0) then
!  When multiple friction times are present, the average is used for the
!  number density in each superparticle.
          if (npar_species > 1) then
            tausp1m = 0.0
          else
            call get_frictiontime(f,fp,p,ineargrid,1,tausp1_par, nochange_opt = .true.)
          endif
!  Need vpm=<|vvp-<vvp>|> to calculate the collisional time-scale.
          vvpm = 0.0
          vpm = 0.0
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            vvpm(ix0-nghost,:) = vvpm(ix0-nghost,:) + fp(k,ivpx:ivpz)
            if (npar_species > 1) then
              call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par,nochange_opt = .true.)
              tausp1m(ix0-nghost) = tausp1m(ix0-nghost) + tausp1_par
            endif
          enddo
          do l = 1,nx
            if (p%np(l) > 1.0) then
              vvpm(l,:) = vvpm(l,:)/p%np(l)
              if (npar_species > 1) tausp1m = tausp1m/p%np
            endif
          enddo
!  vpm
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            vpm(ix0-nghost) = vpm(ix0-nghost) + sqrt( (fp(k,ivpx)-vvpm(ix0-nghost,1))**2 + &
                                                      (fp(k,ivpy)-vvpm(ix0-nghost,2))**2 + &
                                                      (fp(k,ivpz)-vvpm(ix0-nghost,3))**2 )
          enddo
          do l = 1,nx
            if (p%np(l) > 1.0) then
              vpm(l) = vpm(l)/p%np(l)
            endif
          enddo
!  The collisional time-scale is 1/tau_coll=nd*vrms*sigma_coll.
!  Inserting Epstein friction time gives 1/tau_coll=3*rhod/rho*vprms/tauf.
          if (npar_species > 1) then
            tau_coll1 = (1.0-coeff_restitution)*p%epsp*vpm*tausp1m
          else
            tau_coll1 = (1.0-coeff_restitution)*p%epsp*vpm*tausp1_par
          endif
!  Limit inverse time-step of collisional cooling if requested.
          if (tau_coll_min > 0.0) then
            where (tau_coll1 > tau_coll1_max) tau_coll1 = tau_coll1_max
          endif
!
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - &
                tau_coll1(ix0-nghost)*(fp(k,ivpx:ivpz)-vvpm(ix0-nghost,:))
            if (lcollisional_heat .or. ldiagnos) then
              call get_rhopswarm(mp_swarm,fp,k,ineargrid(k,:),rhop_swarm_par)
              coll_heat(ix0-nghost) = coll_heat(ix0-nghost) + rhop_swarm_par*tau_coll1(ix0-nghost)* &
                                      sum(fp(k,ivpx:ivpz)*(fp(k,ivpx:ivpz)-vvpm(ix0-nghost,:)))
            endif
          enddo
!
          if (lupdate_courant_dt) dt1_max = max(dt1_max,tau_coll1/cdtp)
        endif
      endif
!
!  More advanced collisional cooling model. Collisions are considered for
!  every possible two-body process in a grid cell.
!
      if (lcollisional_cooling_twobody) then
        do l = 1,nx
! Collisions between particle k and all other particles in the grid cell.
          k = kshepherd(l)
          if (k > 0) then
!  Limit inverse time-step of collisional cooling if requested.
            do while (k /= 0)
              dt1_cool = 0.0
              call get_frictiontime(f,fp,p,ineargrid,k,tausp1_park, nochange_opt = .true.)
              tausp_park = 1/tausp1_park
              tausp_park3 = tausp_park**3
              j = k
              do while (kneighbour(j) /= 0)
!  Collide with the neighbours of k and their neighbours.
                j = kneighbour(j)
                call get_frictiontime(f,fp,p,ineargrid,j,tausp1_parj, nochange_opt = .true.)
                tausp_parj = 1/tausp1_parj
                tausp_parj3 = tausp_parj**3
!  Collision velocity.
                deltavp_vec = fp(k,ivpx:ivpz)-fp(j,ivpx:ivpz)
                deltavp = sqrt( deltavp_vec(1)**2 + deltavp_vec(2)**2 + deltavp_vec(3)**2 )
                vbar_jk = (tausp_parj3*fp(k,ivpx:ivpz)+tausp_park3*fp(j,ivpx:ivpz))/ &
                    (tausp_parj3+tausp_park3)
!  Cooling time-scale.
                call get_rhopswarm(mp_swarm,fp,k,ineargrid(k,:),rhop_swarm_par)
                tau_cool1_par = (1.0-coeff_restitution)* &
                    rhop_swarm_par*deltavp*(tausp_parj+tausp_park)**2/(tausp_parj3+tausp_park3)
                dt1_cool = dt1_cool+tau_cool1_par
!                if (tau_coll_min>0.0) then
!                  if (tau_cool1_par>tau_coll1_max) tau_cool1_par=tau_coll1_max
!                endif
                dfp(j,ivpx:ivpz) = dfp(j,ivpx:ivpz) - tau_cool1_par*(fp(j,ivpx:ivpz)-vbar_jk)
                dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - tau_cool1_par*(fp(k,ivpx:ivpz)-vbar_jk)
              enddo
              if (lupdate_courant_dt) dt1_max = max(dt1_max(l),dt1_cool/cdtp)
!  Go through all possible k.
              k = kneighbour(k)
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
        if (npar_imn(imn) /= 0) then
          vvpm_species = 0.0
          vpm_species = 0.0
          np_species = 0.0
!
!  Calculate mean velocity and number of particles for each species.
!
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            ispecies = assign_species(ipar(k))
            vvpm_species(ix0-nghost,:,ispecies) = vvpm_species(ix0-nghost,:,ispecies) + fp(k,ivpx:ivpz)
            np_species(ix0-nghost,ispecies) = np_species(ix0-nghost,ispecies) + 1.0
          enddo
          do l = 1,nx
            do ispecies = 1,npar_species
              if (np_species(l,ispecies) > 1.0) then
                vvpm_species(l,:,ispecies) = vvpm_species(l,:,ispecies)/np_species(l,ispecies)
              endif
            enddo
          enddo
!  Calculate rms speed for each species.
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            ispecies = assign_species(ipar(k))
            vpm_species(ix0-nghost,ispecies) = vpm_species(ix0-nghost,ispecies) + sqrt( &
                (fp(k,ivpx)-vvpm_species(ix0-nghost,1,ispecies))**2 + &
                (fp(k,ivpy)-vvpm_species(ix0-nghost,2,ispecies))**2 + &
                (fp(k,ivpz)-vvpm_species(ix0-nghost,3,ispecies))**2 )
          enddo
          do l = 1,nx
            do ispecies = 1,npar_species
              if (np_species(l,ispecies) > 1.0) then
                vpm_species(l,ispecies) = vpm_species(l,ispecies)/np_species(l,ispecies)
              endif
            enddo
          enddo
!
!  Collisional drag force time-scale between particles i and j with R_i < R_j.
!
!    tau_ji = tau_j^3/(tau_i+tau_j)^2/(deltav_ij/cs*rhoi/rhog)
!    tau_ij = tau_ji*rho_j/rho_i
!
          do ispecies = 1,npar_species
            do jspecies = ispecies,npar_species
              vcoll = &
                  sqrt(vpm_species(:,ispecies)**2+vpm_species(:,ispecies)**2 + &
                  (vvpm_species(:,1,ispecies)-vvpm_species(:,1,jspecies))**2 + &
                  (vvpm_species(:,2,ispecies)-vvpm_species(:,2,jspecies))**2 + &
                  (vvpm_species(:,3,ispecies)-vvpm_species(:,3,jspecies))**2)
              call get_rhopswarm(mp_swarm,fp,k,m,n,rhop_swarm_mn)
              tau_coll1_species(:,jspecies,ispecies) = &
                  vcoll*np_species(:,ispecies)*rhop_swarm_mn*p%rho1 / (tausp_species(jspecies)**3/ &
                  (tausp_species(ispecies)+tausp_species(jspecies))**2 )
              where (np_species(:,ispecies) /= 0.0) &
                  tau_coll1_species(:,ispecies,jspecies) = &
                  tau_coll1_species(:,jspecies,ispecies)*np_species(:,jspecies)/np_species(:,ispecies)
            enddo
          enddo
!
          tau_coll1_tot = 0.0
          do ispecies = 1,npar_species
            do jspecies = 1,npar_species
              tau_coll1_tot(:,ispecies) = tau_coll1_tot(:,ispecies)+tau_coll1_species(:,ispecies,jspecies)
            enddo
          enddo
!  Limit inverse time-step of collisional cooling if requested.
          if (tau_coll_min > 0.0) then
            do ispecies = 1,npar_species
              do l = 1,nx
                if (tau_coll1_tot(l,ispecies) > tau_coll1_max) then
                  tau_coll1_species(l,ispecies,:) = tau_coll1_species(l,ispecies,:)* &
                                                    tau_coll1_max/tau_coll1_tot(l,ispecies)
                endif
              enddo
            enddo
            tau_coll1_tot = 0.0
            do ispecies = 1,npar_species
              do jspecies = 1,npar_species
                tau_coll1_tot(:,ispecies) = tau_coll1_tot(:,ispecies)+tau_coll1_species(:,ispecies,jspecies)
              enddo
            enddo
          endif
          if (lupdate_courant_dt) then
            do ispecies = 1,npar_species
              dt1_max = max(dt1_max,tau_coll1_tot(:,ispecies)/cdtp)
            enddo
          endif
!  Add to equation of motion.
          do k = k1_imn(imn),k2_imn(imn)
            ix0 = ineargrid(k,1)
            ispecies = assign_species(ipar(k))
            do jspecies = 1,npar_species
              dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - &
                  tau_coll1_species(ix0-nghost,ispecies,jspecies)* &
                  (fp(k,ivpx:ivpz)-vvpm_species(ix0-nghost,:,jspecies))
              if (lcollisional_heat .or. ldiagnos) then
                call get_rhopswarm(mp_swarm,fp,k,ineargrid(k,:),rhop_swarm_par)
                coll_heat(ix0-nghost) = coll_heat(ix0-nghost) + rhop_swarm_par* &
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
      if (lentropy .and. lcollisional_heat) df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%rho1*p%TT1*coll_heat
!
!  Diagnostics.
!
      if (ldiagnos) call sum_mn_name(coll_heat,idiag_decollp)
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
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3) :: ineargrid
!
      real, dimension(nx,3) :: vvpm
      real :: tausp1_par, tausp1_par_mod, tausp1_par_org
      integer :: k, l, ix0
!
      if (npar_imn(imn) /= 0) then
!
!  Calculate mean particle velocity.
!
        vvpm = 0.0
        do k = k1_imn(imn),k2_imn(imn)
          ix0 = ineargrid(k,1)
          vvpm(ix0-nghost,:) = vvpm(ix0-nghost,:) + fp(k,ivpx:ivpz)
        enddo
        do l = 1,nx
          if (p%np(l) > 1.0) vvpm(l,:) = vvpm(l,:)/p%np(l)
        enddo
!
!  Compare actual and modified friction time and apply the difference as
!  friction force relative to mean particle velocity.
!
        do k = k1_imn(imn),k2_imn(imn)
          call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par_mod,nochange_opt = .false.)
          call get_frictiontime(f,fp,p,ineargrid,k,tausp1_par_org,nochange_opt = .true.)
          tausp1_par = tausp1_par_org-tausp1_par_mod
          ix0 = ineargrid(k,1)
          dfp(k,ivpx:ivpz) = dfp(k,ivpx:ivpz) - tausp1_par*(fp(k,ivpx:ivpz)-vvpm(ix0-nghost,:))
        enddo
      endif
!
    endsubroutine compensate_friction_increase
!***********************************************************************
    subroutine calc_gas_velocity_shell_call(k1, uup, fp)
!
      use Special, only: special_calc_particles
!
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(3) :: uup
      integer :: k1
!
      vel_call = .true.
      uup_shared = fp(k1,ixp:izp)
!      call special_calc_particles(fp,ineargrid)
      uup = uup_shared
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
      real, dimension(mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
      real, dimension(nx,3) :: vvpm, dvp2m
      integer :: inx0, k,l
      type (pencil_case) :: p
!
!  Initialize the variables
!
      vvpm = 0.0
      dvp2m = 0.0
!
!  Calculate the average velocity at each cell
!  if there are particles in the pencil only
!
      if (npar_imn(imn) /= 0) then
!
        do k = k1_imn(imn),k2_imn(imn)
          inx0 = ineargrid(k,1)-nghost
          vvpm(inx0,:) = vvpm(inx0,:) + fp(k,ivpx:ivpz)
        enddo
        do l = 1,nx
          if (p%np(l) > 1.0) vvpm(l,:) = vvpm(l,:)/p%np(l)
        enddo
!
!  Get the residual in quadrature, dvp2m. Need vvpm calculated above.
!
        do k = k1_imn(imn),k2_imn(imn)
          inx0 = ineargrid(k,1)-nghost
          dvp2m(inx0,1) = dvp2m(inx0,1)+(fp(k,ivpx)-vvpm(inx0,1))**2
          dvp2m(inx0,2) = dvp2m(inx0,2)+(fp(k,ivpy)-vvpm(inx0,2))**2
          dvp2m(inx0,3) = dvp2m(inx0,3)+(fp(k,ivpz)-vvpm(inx0,3))**2
        enddo
        do l = 1,nx
          if (p%np(l) > 1.0) dvp2m(l,:) = dvp2m(l,:)/p%np(l)
        enddo
!
      endif
!
!  Output the diagnostics
!
      call sum_mn_name(dvp2m(:,1),idiag_dvpx2m)
      call sum_mn_name(dvp2m(:,2),idiag_dvpy2m)
      call sum_mn_name(dvp2m(:,3),idiag_dvpz2m)
      if (idiag_dvpm /= 0)   call sum_mn_name(dvp2m(:,1)+dvp2m(:,2)+dvp2m(:,3), &
          idiag_dvpm,lsqrt = .true.)
      if (idiag_dvpmax /= 0) call max_mn_name(dvp2m(:,1)+dvp2m(:,2)+dvp2m(:,3), &
          idiag_dvpmax,lsqrt = .true.)
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
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(k1_imn(imn):k2_imn(imn)), intent(inout) :: rep
!
      real, dimension(k1_imn(imn):k2_imn(imn)) :: nu
      character(len=labellen) :: ivis=''
      real :: nu_
      integer :: k
!
!  Find the kinematic viscosity.
!  Check whether we want to override the usual viscosity for the drag law.
!
      if (lnu_draglaw) then
        nu = nu_draglaw
      else
        call getnu(nu_input=nu_,IVIS=ivis)
        if (ivis == 'nu-const') then
          nu = nu_
        elseif (ivis == 'nu-mixture') then
          nu = interp_nu
        elseif (ivis == 'rho-nu-const') then
          nu = nu_/interp_rho(k1_imn(imn):k2_imn(imn))
        elseif (ivis == 'sqrtrho-nu-const') then
          nu = nu_/sqrt(interp_rho(k1_imn(imn):k2_imn(imn)))
        elseif (ivis == 'nu-therm') then
          nu = nu_*sqrt(interp_TT(k1_imn(imn):k2_imn(imn)))
        elseif (ivis == 'mu-therm') then
          nu = nu_*sqrt(interp_TT(k1_imn(imn):k2_imn(imn)))/interp_rho(k1_imn(imn):k2_imn(imn))
        else
          call fatal_error('calc_pencil_rep','no such ivis: '//trim(ivis))
        endif
      endif
!
      if (maxval(nu) == 0.0) call fatal_error('calc_pencil_rep', 'nu (kinematic visc.) must >0')
!
      do k = k1_imn(imn),k2_imn(imn)
        rep(k) = 2.0 * sqrt(sum((interp_uu(k,:) - fp(k,ivpx:ivpz))**2)) / nu(k)
      enddo
!
      if (lparticles_radius) then
        rep = rep * fp(k1_imn(imn):k2_imn(imn),iap)
      elseif (particle_radius > 0.0) then
        rep = rep * particle_radius
      else
        call fatal_error('calc_pencil_rep', &
             'unable to calculate the particle Reynolds number without a particle radius')
      endif
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
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(k1_imn(imn):k2_imn(imn)) :: stocunn
!
      real :: dia
      integer :: k
!
      do k = k1_imn(imn),k2_imn(imn)
!
!  Particle diameter
!
        dia = 2.0*fp(k,iap)
!
        if (lstocunn1) then
          stocunn(k) = 1.
        else
          stocunn(k) = 1.+2.*mean_free_path_gas/dia*(1.257+0.4*exp(-0.55*dia/mean_free_path_gas))
        endif
!
      enddo
!
    endsubroutine calc_stokes_cunningham
!**********************************************************************
    subroutine calc_added_mass_beta(fp,k,added_mass_beta)
!
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, intent(in) :: k
      real, intent(out) :: added_mass_beta
!
!  beta for added mass according to beta=3rho_fluid/(2rho_part+rho_fluid)
!  problem: we would have to calculate beta every time for every particle
!
      added_mass_beta = 3*interp_rho(k)/(2*rhopmat+interp_rho(k))
!
    endsubroutine calc_added_mass_beta
!***********************************************************************
    subroutine calc_draglaw_purestokes(fp,k,tausp1_par)
!
!   Calculate relaxation time for particles under Pure Stokes drag
!
!   6-Aug-15/nils+dhruba: coded
!
      use Viscosity, only: getnu
      use Particles_radius
!
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, intent(in) :: k
      real, intent(out) :: tausp1_par
!
      character(len=labellen) :: ivis=''
      real :: dia, nu, nu_
!
!  Find the kinematic viscosity
!
      call getnu(nu_input=nu_,ivis=ivis)
      if (ivis == 'nu-const') then
        nu = nu_
      elseif (ivis == 'nu-mixture') then
        nu = interp_nu(k)
      elseif (ivis == 'rho-nu-const') then
        nu = nu_/interp_rho(k)
      elseif (ivis == 'sqrtrho-nu-const') then
        nu = nu_/sqrt(interp_rho(k))
      elseif (ivis == 'nu-therm') then
        nu = nu_*sqrt(interp_TT(k))
      elseif (ivis == 'mu-therm') then
        nu = nu_*sqrt(interp_TT(k)) /interp_rho(k)
      else
        call fatal_error('calc_draglaw_purestokes','no such ivis: '//trim(ivis))
      endif
!
!  Particle diameter
!
      if (.not. lparticles_radius) call fatal_error('calc_draglaw_purestokes', &
               'need particles_radius module to calculate the relaxation time')
!
      dia = 2.0*fp(k,iap)
!
!  Relaxation time:
!
      tausp1_par = 18.0*nu/((rhopmat/interp_rho(k))*dia**2)
!
    endsubroutine calc_draglaw_purestokes
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
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, intent(in) :: k
      real, intent(in) :: rep, stocunn
      real, intent(out) :: tausp1_par
!
      character(len=labellen) :: ivis=''
      real :: cdrag, dia, nu, nu_
!
!  Find the kinematic viscosity.
!  Check whether we want to override the usual viscosity for the drag law.
!
      if (lnu_draglaw) then
        nu = nu_draglaw
      else
!
!  Use usual viscosity for the drag law.
!
        call getnu(nu_input=nu_,ivis=ivis)
        if (ivis == 'nu-const') then
          nu = nu_
        elseif (ivis == 'nu-mixture') then
          nu = interp_nu(k)
        elseif (ivis == 'rho-nu-const') then
          nu = nu_/interp_rho(k)
        elseif (ivis == 'sqrtrho-nu-const') then
          nu = nu_/sqrt(interp_rho(k))
        elseif (ivis == 'nu-therm') then
          nu = nu_*sqrt(interp_TT(k))
        elseif (ivis == 'mu-therm') then
          nu = nu_*sqrt(interp_TT(k)) /interp_rho(k)
        else
          call fatal_error('calc_draglaw_steadystate','no such ivis: '//trim(ivis))
        endif
      endif
!
!  Particle diameter
!
      if (.not. lparticles_radius) call fatal_error('calc_draglaw_steadystate', &
            'need particles_radius module to calculate the relaxation time')
!
      dia = 2.0*fp(k,iap)
!
!  Calculate drag coefficent pre-factor:
!
      if (rep < 1) then
        cdrag = 1.0
      elseif (rep > 1000) then
        cdrag = 0.44*rep/24.0
      else
        cdrag = (1.+0.15*rep**0.687)
      endif
!
!  Relaxation time:
!
      tausp1_par = 18.0*cdrag*nu/((rhopmat/interp_rho(k))*stocunn*dia**2)
!
    endsubroutine calc_draglaw_steadystate
!***********************************************************************
    subroutine calc_brownian_force_Li_Ahmadi(fp,k,ineark,stocunn,force)
!
!  Calculate the Brownian force contribution due to the random thermal motions
!  of the gas molecules.
!      
!  See A. Li and G. Ahmadi. "Dispersion and Deposition of Spherical
!      Particles from Point Sources in a Turbulent Channel Flow".
!      Aerosol Science and Technology. 16. 209–226. 1992.
!
!  28-jul-08/kapelrud: coded
!  13-nov-24/axel+nils: renamed
!
      use General, only: normal_deviate
      use Viscosity, only: getnu
!
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, intent(in) :: k
      integer, dimension(3) :: ineark
      real :: stocunn
      real, dimension(3), intent(out) :: force
!
      character(len=labellen) :: ivis=''
      real :: Szero, dia, TT, rhop_swarm_par, nu, nu_
!
!  Find kinematic viscosity
!
      call getnu(nu_input=nu_,IVIS=ivis)
      if (ivis == 'nu-const') then
        nu = nu_
      elseif (ivis == 'nu-mixture') then
        nu = interp_nu(k)
      elseif (ivis == 'rho-nu-const') then
        nu = nu_/interp_rho(k)
      elseif (ivis == 'sqrtrho-nu-const') then
        nu = nu_/sqrt(interp_rho(k))
      elseif (ivis == 'nu-therm') then
        nu = nu_*sqrt(interp_TT(k))
      elseif (ivis == 'mu-therm') then
        nu = nu_*sqrt(interp_TT(k)) /interp_rho(k)
      else
        call fatal_error('calc_brownian_force_Li_Ahmadi','no such ivis: '//trim(ivis))
      endif
!
!  Particle diameter:
!
      dia = 2.0*fp(k,iap)
!
!  Get zero mean, unit variance Gaussian random numbers:
!
      call normal_deviate(force(1))
      call normal_deviate(force(2))
      call normal_deviate(force(3))
!
      if (interp%lTT) then
        TT = interp_TT(k)
      else
        TT = brownian_T0
      endif
!
      call get_rhopswarm(mp_swarm,fp,k,ineark,rhop_swarm_par)
!
!      Szero = 216*nu*k_B*TT*pi_1/ &
!          (dia**5*stocunn*rhop_swarm_par**2/interp_rho(k))
      Szero = 216*nu*k_B*TT*pi_1/(dia**5*stocunn*rhopmat**2/interp_rho(k))
!
!  https://en.wikipedia.org/wiki/Brownian_motion
!  .3 * urms * lmfp = D=kB*T/(6pi*eta*a)
!  .3 * urms^2 * tau = D=kB*T/(6pi*eta*a)
!
!  du/dt = sqrt(dt)
!  .3 * urms^2 * tau = D=kB*T/(6pi*eta*a)
!
      if (dt == 0.0) then
        force = 0.0
      else
        force = force*sqrt(Szero/dt)
      endif
!
    endsubroutine calc_brownian_force_Li_Ahmadi
    !***********************************************************************
    subroutine calc_brownian_force(f,fp,k,ineark,force,tausp1)
!
!  Calculate the Brownian force contribution due to the random thermal motions
!  of the gas molecules.
!      
!  13-nov-24/(axel + nils): coded
!
      use General, only: gaunoise_number
      use Viscosity, only: getnu
      !
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, intent(in) :: k
      integer, dimension(3) :: ineark
      real, dimension(2) :: gn
      real, dimension(3), intent(out) :: force
!
      character(len=labellen) :: ivis=''
      real :: TT, tausp1, Szero, mp, eta
      !
      !  Get zero mean, unit variance Gaussian random numbers:
      !
      call gaunoise_number(gn)
      force(1)=gn(1)
      force(2)=gn(2)
      call gaunoise_number(gn)
      force(3)=gn(1)
      !
      if (interp%lTT) then
        TT = interp_TT(k)
      else
        TT = brownian_T0
      endif
      !
      if (ldraglaw_epstein) then
        tausp1 = interp_cs(k)*interp_rho(k)/(fp(k,iap)*rhopmat)
      elseif (ldraglaw_purestokes .or. ldraglaw_steadystate) then
        ! NILS: This should probably be used for all non-epstein drag-laws
        tausp1 = 18*interp_nu(k)*interp_rho(k)/(rhopmat*4*fp(k,iap)**2)
      else
        call fatal_error("calc_brownian_force","no such draglaw")
      endif
      mp=four_pi_over_three*fp(k,iap)**3*rhopmat
! Calling it Szero to be consistent with the other brownian routine
      Szero=2*tausp1*k_B*TT/mp
!
      if (dt == 0.0) then
        force = 0.0
      else
        !
        ! NILS: The Browian forces are correctly implemented only when
        ! lfollow_gas=T. For lfollow_gas=F, everything is correct for
        ! itorder=1 (and almost correct also for itorder=2), but for itorder=3
        ! the motions are somewhat underpredicted. I have tried to compensate
        ! for this by compensating the time step in the line below with
        ! beta_ts, but it is not enough. I do not know how to resolve this
        ! issue.
        !
        if (lfollow_gas) then
          force = force*sqrt(Szero/dt)
        else
          force = force*sqrt(Szero/(beta_ts(itorder)*dt))
        endif
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
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, intent(in) :: k
      integer, dimension(3) :: ineark
      real, dimension(3), intent(out) :: force
      integer :: i
      real :: Inf=1e14
!
      real, dimension(3) :: temp_grad
      real TT,dyn_visc,nu_,Kn,phi_therm,mass_p
      real Ktc,Ce,Cm,Cint
      character(len=labellen) :: ivis=''
!
      call keep_compiler_quiet(ineark)
!
      Ktc = 1.10
      Cm = 1.13
      Ce = 2.17
      if (interp%lTT) then
        TT = interp_TT(k)
      else
        TT = thermophoretic_T0
      endif
      !Find dynamic viscosity
      call getnu(nu_input=nu_,ivis=ivis)
      if (ivis == 'nu-const') then
        dyn_visc = nu_*interp_rho(k)
      elseif (ivis == 'nu-mixture') then
        dyn_visc = interp_nu(k)*interp_rho(k)
      elseif (ivis == 'rho-nu-const') then
        dyn_visc = nu_
      elseif (ivis == 'sqrtrho-nu-const') then
        dyn_visc = nu_*sqrt(interp_rho(k))
      elseif (ivis == 'nu-therm') then
        dyn_visc = nu_*interp_rho(k)*sqrt(TT)
      elseif (ivis == 'mu-therm') then
        dyn_visc = nu_*sqrt(TT)
      else
        call fatal_error('calc_thermophoretic_force','no such ivis: '//trim(ivis))
      endif
      Cint = 0.5
      if (interp%lgradTT) then
        temp_grad=interp_gradTT(k,:)
      else
        temp_grad = temp_grad0
      endif
!
      Kn = mean_free_path_gas/fp(k,iap)
!
      mass_p = (4.0*pi/3.0)*rhopmat*fp(k,iap)**3
      if (thermophoretic_eq == 'near_continuum') then
        phi_therm = -9*pi/cond_ratio
      elseif (thermophoretic_eq == 'transition') then
        phi_therm = -12.0*pi*(Ktc*(1.0+cond_ratio*Ce*Kn)+3.0*Cm*Kn*(1.0-cond_ratio+cond_ratio*Ce*Kn)) &
                    /((1.0+3.0*Kn*exp(-Cint/Kn))*(1.0+3.0*Cm*Kn)*(2.0+cond_ratio+2.0*cond_ratio*Ce*Kn))
      elseif (thermophoretic_eq == 'free_molecule') then
        phi_therm = 0.0
      else
        call fatal_error('calc_thermophoretic_force','no such thermophoretic equation: '// &
                         trim(thermophoretic_eq))
      endif
!
      force = (fp(k,iap)*temp_grad*dyn_visc**2*phi_therm)/(TT*interp_rho(k)*mass_p)
      
      do i = 1,3
        if (force(i) < -Inf .or. force(i) > Inf) then
          print*, 'Force xyz:', force, 'Temp_grad xyz:',temp_grad
          print*, 'Thermophoretic term ',phi_therm,'Dyn_Visc ',dyn_visc,'TT',TT
          print*, 'Rho_gas', interp_rho(k), 'M_p',mass_p
          call fatal_error('calc_thermophoretic_force','infs in force')
        endif
      enddo
!
    endsubroutine calc_thermophoretic_force
!***********************************************************************
    subroutine read_particles_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read (parallel_unit, NML=particles_init_pars, IOSTAT=iostat)
!
! if we are using particles_potential
!
      if (lparticles_potential) call read_particles_pot_init_pars(iostat)
!
    endsubroutine read_particles_init_pars
!***********************************************************************
    subroutine write_particles_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit, NML=particles_init_pars)
!
! if we are using particles_potential
!
      if (lparticles_potential) call write_particles_pot_init_pars(unit)
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read (parallel_unit, NML=particles_run_pars, IOSTAT=iostat)
!
! if we are using particles_potential
!
      if (lparticles_potential) call read_particles_pot_run_pars(iostat)
!
!  If we have bubbles, the advective derivative has to be saved in
!  an auxiliary variable
!  COMMENT: This would be better to do in a step between registering and
!  initializing. Such a hook does not exist at the moment.
!
      if (lbubble) ladv_der_as_aux = .true.
!
    endsubroutine read_particles_run_pars
!***********************************************************************
    subroutine write_particles_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit, NML=particles_run_pars)
!
! if we are using particles_potential
!
      if (lparticles_potential) call write_particles_pot_run_pars(unit)
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
      real, dimension(mx,my,mz,mfarray) :: f
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
      use General, only: itoa
      use Particles_caustics, only: rprint_particles_caustics
      use Particles_tetrad, only: rprint_particles_tetrad
      use Particles_potential, only: rprint_particles_potential
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamez, inamey, inamex, inamexy, inamexz, inamer, inamerz
      integer :: k
      character(len=intlen) :: srad
!
!  Reset everything in case of reset.
!
      if (lreset) then
        idiag_xpm = 0
        idiag_ypm = 0
        idiag_zpm = 0
        idiag_xpmax = 0
        idiag_ypmax = 0
        idiag_zpmax = 0
        idiag_xpmin = 0
        idiag_ypmin = 0
        idiag_zpmin = 0
        idiag_vrelpabsm = 0
        idiag_xp2m = 0
        idiag_yp2m = 0
        idiag_zp2m = 0
        idiag_rpm = 0
        idiag_rp2m = 0
        idiag_vpxm = 0
        idiag_vpym = 0
        idiag_vpzm = 0
        idiag_urel = 0
        idiag_vpxvpym = 0
        idiag_vpxvpzm = 0
        idiag_vpyvpzm = 0
        idiag_vpx2m = 0
        idiag_vpy2m = 0
        idiag_vpz2m = 0
        idiag_ekinp = 0
        idiag_vpxmax = 0
        idiag_vpymax = 0
        idiag_vpzmax = 0
        idiag_vpxmin = 0
        idiag_vpymin = 0
        idiag_vpzmin = 0
        idiag_rhopart = 0
        idiag_rhopvpxm = 0
        idiag_rhopvpym = 0
        idiag_rhopvpzm = 0
        idiag_rhopvpysm = 0
        idiag_rhopvpxt = 0
        idiag_rhopvpyt = 0
        idiag_rhopvpzt = 0
        idiag_lpxm = 0
        idiag_lpym = 0
        idiag_lpzm = 0
        idiag_lpx2m = 0
        idiag_lpy2m = 0
        idiag_lpz2m = 0
        idiag_npm = 0
        idiag_np2m = 0
        idiag_npmax = 0
        idiag_npmin = 0
        idiag_dtdragp = 0
        idiag_dedragp = 0
        idiag_rhopm = 0
        idiag_rhoprms = 0
        idiag_rhop2m = 0
        idiag_rhopmax = 0
        idiag_rhopmin = 0
        idiag_decollp = 0
        idiag_rhopmphi = 0
        idiag_epspmin = 0
        idiag_epspmax = 0
        idiag_epspm = 0
!
        idiag_nparmin = 0
        idiag_nparmax = 0
        idiag_nparsum = 0
        idiag_nmigmax = 0
        idiag_nmigmmax = 0
        idiag_mpt = 0
        idiag_npmx = 0
        idiag_npmy = 0
        idiag_npmz = 0
        idiag_epotpm = 0
        idiag_rhopmx = 0
        idiag_rhopmy = 0
        idiag_rhopmz = 0
        idiag_rhop2mx = 0
        idiag_rhop2my = 0
        idiag_rhop2mz = 0
        idiag_epspmx = 0
        idiag_epspmy = 0
        idiag_epspmz = 0
        idiag_rhopmxy = 0
        idiag_rhopmxz = 0
        idiag_rhopmr = 0
        idiag_sigmap = 0
        idiag_rpvpxmx = 0
        idiag_rpvpymx = 0
        idiag_rpvpzmx = 0
        idiag_rpvpx2mx = 0
        idiag_rpvpy2mx = 0
        idiag_rpvpz2mx = 0
        idiag_dvpx2m = 0
        idiag_dvpy2m = 0
        idiag_dvpz2m = 0
        idiag_dvpmax = 0
        idiag_dvpm = 0
        idiag_nparpmax = 0
        idiag_eccpxm = 0
        idiag_eccpym = 0
        idiag_eccpzm = 0
        idiag_eccpx2m = 0
        idiag_eccpy2m = 0
        idiag_eccpz2m = 0
        idiag_npargone = 0
        idiag_vpyfull2m = 0
        idiag_deshearbcsm = 0
        idiag_npar_found = 0
        idiag_npmxy = 0
        idiag_vprms = 0
        idiag_npvzmz = 0
        idiag_npvz2mz = 0
        idiag_nptz = 0
        idiag_Shm = 0
        idiag_npuzmz = 0

        idiag_latentheatm = 0
        idiag_ffcondposm = 0
        idiag_ffcondm = 0
        idiag_ffcondnegm = 0
      endif
!
!  Run through all possible names that may be listed in print.in.
!
      if (lroot .and. ip < 14) print*,'rprint_particles: run through parse list'
      do iname = 1,nname
        call parse_name(iname,cname(iname),cform(iname),'nparsum',idiag_nparsum)
        call parse_name(iname,cname(iname),cform(iname),'nparmin',idiag_nparmin)
        call parse_name(iname,cname(iname),cform(iname),'nparmax',idiag_nparmax)
        call parse_name(iname,cname(iname),cform(iname),'nparpmax',idiag_nparpmax)
        call parse_name(iname,cname(iname),cform(iname),'vrelpabsm',idiag_vrelpabsm)
        call parse_name(iname,cname(iname),cform(iname),'xpm',idiag_xpm)
        call parse_name(iname,cname(iname),cform(iname),'ypm',idiag_ypm)
        call parse_name(iname,cname(iname),cform(iname),'zpm',idiag_zpm)
        call parse_name(iname,cname(iname),cform(iname),'xpmin',idiag_xpmin)
        call parse_name(iname,cname(iname),cform(iname),'ypmin',idiag_ypmin)
        call parse_name(iname,cname(iname),cform(iname),'zpmin',idiag_zpmin)
        call parse_name(iname,cname(iname),cform(iname),'xpmax',idiag_xpmax)
        call parse_name(iname,cname(iname),cform(iname),'ypmax',idiag_ypmax)
        call parse_name(iname,cname(iname),cform(iname),'zpmax',idiag_zpmax)
        call parse_name(iname,cname(iname),cform(iname),'xp2m',idiag_xp2m)
        call parse_name(iname,cname(iname),cform(iname),'yp2m',idiag_yp2m)
        call parse_name(iname,cname(iname),cform(iname),'zp2m',idiag_zp2m)
        call parse_name(iname,cname(iname),cform(iname),'rpm',idiag_rpm)
        call parse_name(iname,cname(iname),cform(iname),'rp2m',idiag_rp2m)
        call parse_name(iname,cname(iname),cform(iname),'urel',idiag_urel)
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
        call parse_name(iname,cname(iname),cform(iname),'vpxmin',idiag_vpxmin)
        call parse_name(iname,cname(iname),cform(iname),'vpymin',idiag_vpymin)
        call parse_name(iname,cname(iname),cform(iname),'vpzmin',idiag_vpzmin)
        call parse_name(iname,cname(iname),cform(iname),'vpmax',idiag_vpmax)
        call parse_name(iname,cname(iname),cform(iname),'rhopart',idiag_rhopart)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpxm',idiag_rhopvpxm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpym',idiag_rhopvpym)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpzm',idiag_rhopvpzm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpysm',idiag_rhopvpysm)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpxt',idiag_rhopvpxt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpyt',idiag_rhopvpyt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpzt',idiag_rhopvpzt)
        call parse_name(iname,cname(iname),cform(iname),'rhopvpysm',idiag_rhopvpysm)
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
        call parse_name(iname,cname(iname),cform(iname),'epspm',idiag_epspm)
        call parse_name(iname,cname(iname),cform(iname),'epspmin',idiag_epspmin)
        call parse_name(iname,cname(iname),cform(iname),'epspmax',idiag_epspmax)
        call parse_name(iname,cname(iname),cform(iname),'rhopmphi',idiag_rhopmphi)
        call parse_name(iname,cname(iname),cform(iname),'nmigmax',idiag_nmigmax)
        call parse_name(iname,cname(iname),cform(iname),'nmigmmax',idiag_nmigmmax)
        call parse_name(iname,cname(iname),cform(iname),'mpt',idiag_mpt)
        call parse_name(iname,cname(iname),cform(iname),'dvpx2m',idiag_dvpx2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpy2m',idiag_dvpy2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpz2m',idiag_dvpz2m)
        call parse_name(iname,cname(iname),cform(iname),'dvpm',idiag_dvpm)
        call parse_name(iname,cname(iname),cform(iname),'dvpmax',idiag_dvpmax)
        call parse_name(iname,cname(iname),cform(iname),'dedragp',idiag_dedragp)
        call parse_name(iname,cname(iname),cform(iname),'decollp',idiag_decollp)
        call parse_name(iname,cname(iname),cform(iname),'epotpm',idiag_epotpm)
        call parse_name(iname,cname(iname),cform(iname),'npargone',idiag_npargone)
        call parse_name(iname,cname(iname),cform(iname),'npar_found',idiag_npar_found)
        call parse_name(iname,cname(iname),cform(iname),'vpyfull2m',idiag_vpyfull2m)
        call parse_name(iname,cname(iname),cform(iname),'vprms',idiag_vprms)
        call parse_name(iname,cname(iname),cform(iname),'Shm',idiag_Shm)
        call parse_name(iname,cname(iname),cform(iname),'deshearbcsm',idiag_deshearbcsm)
        call parse_name(iname,cname(iname),cform(iname),'latentheatm',idiag_latentheatm)
        call parse_name(iname,cname(iname),cform(iname),'ffcondposm',idiag_ffcondposm)
        call parse_name(iname,cname(iname),cform(iname),'ffcondm',idiag_ffcondm)
        call parse_name(iname,cname(iname),cform(iname),'ffcondnegm',idiag_ffcondnegm)
      enddo
!
!  Check for those quantities for which we want x-averages.
!
      do inamex = 1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'npmx',idiag_npmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhopmx',idiag_rhopmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhop2mx',idiag_rhop2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'epspmx',idiag_epspmx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpxmx', idiag_rpvpxmx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpymx', idiag_rpvpymx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpzmx', idiag_rpvpzmx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpx2mx', idiag_rpvpx2mx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpy2mx', idiag_rpvpy2mx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rpvpz2mx', idiag_rpvpz2mx)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamey = 1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'npmy',idiag_npmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhopmy',idiag_rhopmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhop2my',idiag_rhop2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'epspmy',idiag_epspmy)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamez = 1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'npmz',idiag_npmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhopmz',idiag_rhopmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhop2mz',idiag_rhop2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'epspmz',idiag_epspmz)
        do k = 1,ndustrad
          srad = itoa(k)
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'npvzmz'//trim(srad),idiag_npvzmz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'npvz2mz'//trim(srad),idiag_npvz2mz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'npuzmz'//trim(srad),idiag_npuzmz(k))
          call parse_name(inamez,cnamez(inamez),cformz(inamez),'nptz'//trim(srad),idiag_nptz(k))
        enddo
!
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamexy = 1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'npmxy',idiag_npmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhopmxy',idiag_rhopmxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'sigmap',idiag_sigmap)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamexz = 1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'rhopmxz',idiag_rhopmxz)
      enddo
!
!  Check for those quantities for which we want phiz-averages.
!
      do inamer = 1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'rhopmr',idiag_rhopmr)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do inamerz = 1,nnamerz
        call parse_name(inamerz,cnamerz(inamerz),cformrz(inamerz),'rhopmphi',idiag_rhopmphi)
      enddo
!
!    If we are using caustics
!
      if (lparticles_caustics) call rprint_particles_caustics(lreset,lwrite)
!
!    If we are using tetrad
!
      if (lparticles_tetrad) call rprint_particles_tetrad(lreset,lwrite)
!
! ... and potential
!
      if (lparticles_potential) call rprint_particles_potential(lreset,lwrite)
!
    endsubroutine rprint_particles
!***********************************************************************
    subroutine particles_final_clean_up
!
!  cleanup
!
      use Particles_potential, only: particles_potential_clean_up
!
      if (lparticles_potential) call particles_potential_clean_up()
!
    endsubroutine particles_final_clean_up
!***********************************************************************
    subroutine periodic_boundcond_on_aux(f)
!
! Impose periodic boundary condition on gradu as auxiliary variable
!
      use Boundcond, only: set_periodic_boundcond_on_aux
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      if (lparticles_grad) then
        if (iguij /=  0) then
          call set_periodic_boundcond_on_aux(f,igradu11)
          call set_periodic_boundcond_on_aux(f,igradu12)
          call set_periodic_boundcond_on_aux(f,igradu13)
          call set_periodic_boundcond_on_aux(f,igradu21)
          call set_periodic_boundcond_on_aux(f,igradu22)
          call set_periodic_boundcond_on_aux(f,igradu23)
          call set_periodic_boundcond_on_aux(f,igradu31)
          call set_periodic_boundcond_on_aux(f,igradu32)
          call set_periodic_boundcond_on_aux(f,igradu33)
        endif
      endif
!
    endsubroutine periodic_boundcond_on_aux
!***********************************************************************
    subroutine calc_relative_velocity(f,fp,ineargrid)
!
      use Diagnostics
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid

      real, dimension(3) :: uup, rel_vel_sing
      real, dimension(npar_loc) :: rel_vel
      integer :: k
!
!  Calculate particle relative velocity
!
      rel_vel = 0.0
!
      do k = 1,npar_loc
        call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup,ineargrid(k,:),0,0)
        rel_vel_sing = (fp(k,ivpx:ivpz)-uup)**2
        rel_vel(k) = sqrt(sum(rel_vel_sing))
      enddo
!
      call sum_par_name(rel_vel,idiag_vrelpabsm)
!
    endsubroutine calc_relative_velocity
!***********************************************************************
    subroutine diffuse_backreaction(f,df)
!
      use Boundcond
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
!
      integer :: i
!
!  passive scalar consumption
!
      if (.not. lpencil_check_at_work) then
        if (ldiffuse_passive .and. lpscalar_sink) then
          do i = 1,istep_pass
            call boundconds_x(f,idlncc,idlncc)
            call initiate_isendrcv_bdry(f,idlncc,idlncc)
            call finalize_isendrcv_bdry(f,idlncc,idlncc)
            call boundconds_y(f,idlncc,idlncc)
            call boundconds_z(f,idlncc,idlncc)
            call diffuse_interaction(f(:,:,:,idlncc),ldiff_pass,.False.,rdiffconst_pass)
          enddo
          df(l1:l2,m1:m2,n1:n2,ilncc) = df(l1:l2,m1:m2,n1:n2,ilncc) + f(l1:l2,m1:m2,n1:n2,idlncc)
          f(:,:,:,idlncc) = 0.0
        endif
!
!  particle fluid dragforce
!
        if (ldiffuse_dragf .and. ldragforce_gas_par) then
          do i = 1,istep_dragf
            call boundconds_x(f,idfx,idfz)
            call initiate_isendrcv_bdry(f,idfx,idfz)
            call finalize_isendrcv_bdry(f,idfx,idfz)
            call boundconds_y(f,idfx,idfz)
            call boundconds_z(f,idfx,idfz)
            call diffuse_interaction(f(:,:,:,idfx),ldiff_dragf,.False.,rdiffconst_dragf)
            call diffuse_interaction(f(:,:,:,idfy),ldiff_dragf,.False.,rdiffconst_dragf)
            call diffuse_interaction(f(:,:,:,idfz),ldiff_dragf,.False.,rdiffconst_dragf)
          enddo
          df(l1:l2,m1:m2,n1:n2,iux:iuz)=df(l1:l2,m1:m2,n1:n2,iux:iuz) + f(l1:l2,m1:m2,n1:n2,idfx:idfz)
          f(:,:,:,idfx:idfz) = 0.0
        endif
      endif
!
    endsubroutine diffuse_backreaction
!***********************************************************************
endmodule Particles
