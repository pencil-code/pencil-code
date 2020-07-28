! $Id$
!
!  This modules deals with all aspects of magnetic fields; if no
!  magnetic fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  magnetically relevant subroutines listed in here.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagnetic = .true.
! CPARAM logical, parameter :: lbfield = .false.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED aa(3); a2; aij(3,3); bb(3); bbb(3); ab; ua; exa(3); aps
! PENCILS PROVIDED b2; bf2; bij(3,3); del2a(3); graddiva(3); jj(3); e3xa(3)
! PENCILS PROVIDED bijtilde(3,3),bij_cov_corr(3,3)
! PENCILS PROVIDED j2; jb; va2; jxb(3); jxbr(3); jxbr2; ub; uxb(3); uxb2
! PENCILS PROVIDED uxj(3); chibp; beta; beta1; uga(3); uuadvec_gaa(3); djuidjbi; jo
! PENCILS PROVIDED StokesI; StokesQ; StokesU; StokesQ1; StokesU1
! PENCILS PROVIDED ujxb; oxu(3); oxuxb(3); jxbxb(3); jxbrxb(3)
! PENCILS PROVIDED glnrhoxb(3); del4a(3); del6a(3); oxj(3); diva
! PENCILS PROVIDED jij(3,3); sj; ss12; d6ab
! PENCILS PROVIDED etava; etaj; etaj2; etajrho
! PENCILS PROVIDED cosjb; jparallel; jperp
! PENCILS PROVIDED cosub; bunit(3)
! PENCILS PROVIDED hjj(3); hj2; hjb; coshjb
! PENCILS PROVIDED hjparallel; hjperp; nu_ni1
! PENCILS PROVIDED gamma_A2; clight2; gva(3); vmagfric(3)
! PENCILS PROVIDED bb_sph(3)
!***************************************************************
module Magnetic
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet, loptest
  use Magnetic_meanfield
  use Messages, only: fatal_error,inevitably_fatal_error,warning,svn_id,timing
  use EquationOfState, only: gamma1
  use SharedVariables, only: get_shared_variable  
  use Mpicomm, only: stop_it
!
  implicit none
!
  include 'record_types.h'
  include 'magnetic.h'
!
  interface input_persistent_magnetic
     module procedure input_persist_magnetic_id
     module procedure input_persist_magnetic
  endinterface
!
! Slice precalculation buffers
!
  real, target, dimension (:,:,:), allocatable :: bb_xy, jj_xy, poynting_xy ,bb_sph_xy
  real, target, dimension (:,:,:), allocatable :: bb_xy2,jj_xy2,poynting_xy2,bb_sph_xy2
  real, target, dimension (:,:,:), allocatable :: bb_xy3,jj_xy3,poynting_xy3,bb_sph_xy3
  real, target, dimension (:,:,:), allocatable :: bb_xy4,jj_xy4,poynting_xy4,bb_sph_xy4
  real, target, dimension (:,:,:), allocatable :: bb_xz, jj_xz, poynting_xz ,bb_sph_xz
  real, target, dimension (:,:,:), allocatable :: bb_yz, jj_yz, poynting_yz ,bb_sph_yz
  real, target, dimension (:,:,:), allocatable :: bb_xz2,jj_xz2,poynting_xz2,bb_sph_xz2
!
  real, target, dimension (:,:), allocatable :: b2_xy, jb_xy, j2_xy,  ab_xy
  real, target, dimension (:,:), allocatable :: b2_xy2,jb_xy2,j2_xy2, ab_xy2
  real, target, dimension (:,:), allocatable :: b2_xy3,jb_xy3,j2_xy3, ab_xy3
  real, target, dimension (:,:), allocatable :: b2_xy4,jb_xy4,j2_xy4, ab_xy4
  real, target, dimension (:,:), allocatable :: b2_yz, jb_yz, j2_yz,  ab_yz
  real, target, dimension (:,:), allocatable :: b2_xz, jb_xz, j2_xz,  ab_xz
  real, target, dimension (:,:), allocatable :: b2_xz2,jb_xz2,j2_xz2, ab_xz2
!
  real, target, dimension (:,:), allocatable :: beta1_xy
  real, target, dimension (:,:), allocatable :: beta1_xy2
  real, target, dimension (:,:), allocatable :: beta1_xy3
  real, target, dimension (:,:), allocatable :: beta1_xy4
  real, target, dimension (:,:), allocatable :: beta1_yz
  real, target, dimension (:,:), allocatable :: beta1_xz
  real, target, dimension (:,:), allocatable :: beta1_xz2
!
  real, target, dimension (:,:), allocatable :: aps_xy, aps_xz, aps_yz, aps_xz2
!
!  xy-averaged field
!
  real, dimension (mz,3) :: aamz
  real, dimension (nz,3) :: bbmz,jjmz
  real, dimension (:,:), allocatable :: aamxy
!
! Parameters
!
  integer, parameter :: nresi_max=4
!
  real, dimension (ninit) :: amplaa=0.0, kx_aa=1.0, ky_aa=1.0, kz_aa=1.0
  real, dimension (ninit) :: ampl_ax=0.0, ampl_ay=0.0, ampl_az=0.0
  real, dimension (ninit) :: kx_ax=0.0, kx_ay=0.0, kx_az=0.0
  real, dimension (ninit) :: ky_ax=0.0, ky_ay=0.0, ky_az=0.0
  real, dimension (ninit) :: kz_ax=0.0, kz_ay=0.0, kz_az=0.0
  real, dimension (ninit) :: phase_ax=0.0, phase_ay=0.0, phase_az=0.0
  real, dimension (ninit) :: amplaaJ=0.0, amplaaB=0.0, RFPrad=1.0
  real, dimension (ninit) :: phasex_aa=0.0, phasey_aa=0.0, phasez_aa=0.0
  real, dimension (ninit) :: phase_aa=0.0
  integer, dimension (ninit) :: ll_sh=0, mm_sh=0
  integer :: nzav=0,indzav=0,izav_start=1
  character (len=fnlen) :: source_zav=''
  character (len=labellen), dimension(ninit) :: initaa='nothing'
  character (len=labellen), dimension(3) :: borderaa='nothing'
  character (len=labellen), dimension(nresi_max) :: iresistivity=''
  character (len=labellen) :: ihall_term='const'
!
! Input parameters
!
  complex, dimension(3) :: coefaa=(/0.0,0.0,0.0/), coefbb=(/0.0,0.0,0.0/)
  real, dimension(3) :: B_ext = 0.0, B0_ext = 0.0, ABCaa=1., widthaa=0.5
  real, dimension(3) :: B1_ext, B_ext_inv
  real, dimension(3) :: J_ext=(/0.0,0.0,0.0/)
  real, dimension(3) :: eta_aniso_hyper3=0.0
  real, dimension(2) :: magnetic_xaver_range=(/-max_real,max_real/)
  real, dimension(2) :: magnetic_zaver_range=(/-max_real,max_real/)
  real, dimension(nx) :: xmask_mag
  real, dimension(nz) :: zmask_mag
  real :: sheet_position=1.,sheet_thickness=0.1,sheet_hyp=1.
  real :: t_bext = 0.0, t0_bext = 0.0
  real :: radius=0.1, epsilonaa=0.01, x0aa=0.0, y0aa=0.0, z0aa=0.0
  real :: by_left=0.0, by_right=0.0, bz_left=0.0, bz_right=0.0
  real :: relhel_aa=1.
  real :: bthresh=0.0, bthresh_per_brms=0.0, brms=0.0, bthresh_scl=1.0
  real :: eta1_aniso_ratio=impossible, eta1_aniso=impossible
  real :: eta1_aniso_r=0.0, eta1_aniso_d=0.0
  real :: eta_shock=0.0, eta_shock2=0.0, alp_aniso=0.0, eta_aniso_BB=0.0
  real :: quench_aniso=impossible
  real :: eta_va=0., eta_j=0., eta_j2=0., eta_jrho=0., eta_min=0., &
          etaj20=0., va_min=0., vArms=1.
  real :: rhomin_jxb=0.0, va2max_jxb=0.0, va2max_boris=0.0,cmin=0.0
  real :: omega_Bz_ext=0.0
  real :: mu_r=-0.5 !(still needed for backwards compatibility)
  real :: mu_ext_pot=-0.5,inclaa=0.0
  real :: mu012=0.5 !(=1/2mu0)
  real :: rescale_aa=0.0
  real :: ampl_B0=0.0, D_smag=0.17, B_ext2, B_ext21, B_ext11
  real :: nu_ni=0.0, nu_ni1, hall_term=0.0, battery_term=0.0
  real :: hall_tdep_t0=0.0, hall_tdep_exponent=0.0
  real :: Hhall=0., hall_zdep_exponent=4.0
  real :: initpower_aa=0.0, initpower2_aa=-11./3., cutoff_aa=0.0, ncutoff_aa=1.
  real :: kpeak_aa=10., kgaussian_aa=0., brms_target=1.0, rescaling_fraction=1.0
  real :: phase_beltrami=0.0, ampl_beltrami=0.0
  real :: bmz=0, bmz_beltrami_phase=0.0
  real :: taareset=0.0, daareset=0.0
  real :: center1_x=0.0, center1_y=0.0, center1_z=0.0
  real :: fluxtube_border_width=impossible
  real :: eta_jump=0.0, eta_jump2=0.0, damp=0., two_step_factor=1.
  real :: radRFP=1.
  real :: rnoise_int=impossible,rnoise_ext=impossible
  real :: znoise_int=impossible,znoise_ext=impossible
  real :: mix_factor=0.
  real :: RFPradB=1., RFPradJ=1.
  real :: th_spot=PI/4
  real :: non_ffree_factor=1.
  real :: etaB=0.
  real :: tau_relprof=0.0, tau_relprof1, amp_relprof=1.0 , k_relprof=1.0
  real, pointer :: cp
  real :: dipole_moment=0.0
  real :: eta_power_x=0., eta_power_z=0.
  real :: z1_aa=0., z2_aa=0.
  real :: Pm_smag1=1., k1hel=0., k2hel=max_real
  integer, target :: va2power_jxb = 5
  integer :: nbvec, nbvecmax=nx*ny*nz/4, iua=0
  integer :: N_modes_aa=1, naareset
  logical, pointer :: lrelativistic_eos
  logical :: lpress_equil=.false., lpress_equil_via_ss=.false.
  logical :: lpress_equil_alt=.false., lset_AxAy_zero=.false.
  logical :: llorentzforce=.true., llorentz_rhoref=.false., linduction=.true.
  logical :: ldiamagnetism=.false., lcovariant_magnetic=.false.
  logical :: ladd_global_field=.false. 
  logical :: lresi_eta_const=.false.
  logical :: lresi_eta_tdep=.false., lresi_eta_ztdep=.false.
  logical :: lresi_eta_tdep_t0_norm=.false.
  logical :: lresi_sqrtrhoeta_const=.false.
  logical :: lresi_eta_aniso=.false., lquench_eta_aniso=.false.
  logical :: lresi_etaSS=.false.
  logical :: lresi_hyper2=.false.
  logical :: lresi_hyper3=.false.
  logical :: lresi_hyper3_polar=.false.
  logical :: lresi_hyper3_mesh=.false.
  logical :: lresi_hyper3_csmesh=.false.
  logical :: lresi_hyper3_strict=.false.
  logical :: lresi_zdep, lresi_ydep, lresi_xdep, lresi_rdep, lresi_xydep
  logical, dimension(7) :: lresi_dep=.false.
  logical :: lresi_dust=.false.
  logical :: lresi_hyper3_aniso=.false.
  logical :: lresi_eta_shock=.false.
  logical :: lresi_eta_shock2=.false.
  logical :: lresi_eta_shock_profz=.false.
  logical :: lresi_eta_shock_profr=.false.
  logical :: lresi_eta_shock_perp=.false.
  logical :: lresi_etava=.false.
  logical :: lresi_etaj=.false.
  logical :: lresi_etaj2=.false.
  logical :: lresi_etajrho=.false.
  logical :: lresi_shell=.false.
  logical :: lresi_smagorinsky=.false.
  logical :: lresi_smagorinsky_nusmag=.false.
  logical :: lresi_smagorinsky_cross=.false.
  logical :: lresi_anomalous=.false.
  logical :: lresi_spitzer=.false.
  logical :: lresi_cspeed=.false.
  logical :: lresi_vAspeed=.false.,lalfven_as_aux=.false.
  logical :: lresi_magfield=.false.
  logical :: lresi_eta_proptouz=.false.
  logical, target, dimension (3) :: lfrozen_bb_bot=(/.false.,.false.,.false./)
  logical, target, dimension (3) :: lfrozen_bb_top=(/.false.,.false.,.false./)
  logical :: lohmic_heat=.true., lneutralion_heat=.true.
  logical :: reinitialize_aa=.false.
  logical :: lB_ext_pot=.false., lJ_ext=.false.
  logical :: lforce_free_test=.false.
  logical :: lforcing_cont_aa_local=.false.
  logical :: lEE_as_aux=.false.
  logical :: lbb_as_aux=.false., ljj_as_aux=.false., ljxb_as_aux=.false.
  logical :: lbbt_as_aux=.false., ljjt_as_aux=.false., lua_as_aux=.false.
  logical :: letasmag_as_aux=.false.,ljj_as_comaux=.false.
  logical :: lbb_as_comaux=.false., lB_ext_in_comaux=.true.
  logical :: lbb_sph_as_aux=.false.
  logical :: lbext_curvilinear=.true., lcheck_positive_va2=.false.
  logical :: lreset_aa=.false., lsmooth_jj=.false.
  logical :: lbx_ext_global=.false.,lby_ext_global=.false.,&
             lbz_ext_global=.false.
  logical :: lax_ext_global=.false.,lay_ext_global=.false.,&
             laz_ext_global=.false.
  logical :: lambipolar_diffusion=.false.
  logical :: lskip_projection_aa=.false., lno_second_ampl_aa=.true.
  logical :: lscale_tobox=.true.
!
  namelist /magnetic_init_pars/ &
      B_ext, B0_ext, t_bext, t0_bext, J_ext, lohmic_heat, radius, epsilonaa, &
      ABCaa, x0aa, y0aa, z0aa, widthaa, &
      RFPradB, RFPradJ, by_left, by_right, bz_left, bz_right, relhel_aa, &
      initaa, amplaa, kx_aa, ky_aa, kz_aa, amplaaJ, amplaaB, RFPrad, radRFP, &
      coefaa, coefbb, phase_aa, phasex_aa, phasey_aa, phasez_aa, inclaa, &
      lpress_equil, lpress_equil_via_ss, lset_AxAy_zero, &
      mu_r, mu_ext_pot, lB_ext_pot, &
      alp_aniso, ljj_as_comaux, lsmooth_jj, &
      lforce_free_test, ampl_B0, N_modes_aa, &
      initpower_aa, initpower2_aa, cutoff_aa, ncutoff_aa, kpeak_aa, &
      lscale_tobox, kgaussian_aa, z1_aa, z2_aa, &
      lcheck_positive_va2, lskip_projection_aa, lno_second_ampl_aa, &
      lbb_as_aux, lbb_as_comaux, lB_ext_in_comaux, lEE_as_aux,& 
      ljxb_as_aux, ljj_as_aux, lbext_curvilinear, lbbt_as_aux, ljjt_as_aux, &
      lua_as_aux, lneutralion_heat, center1_x, center1_y, center1_z, &
      fluxtube_border_width, va2max_jxb, va2max_boris, cmin,va2power_jxb, eta_jump, &
      lpress_equil_alt, rnoise_int, rnoise_ext, mix_factor, damp, &
      two_step_factor, th_spot, non_ffree_factor, etaB, ampl_ax, ampl_ay, &
      ampl_az, kx_ax, kx_ay, kx_az, ky_ax, ky_ay, ky_az, kz_ax, kz_ay, kz_az, &
      phase_ax, phase_ay, phase_az, magnetic_xaver_range, amp_relprof, k_relprof, &
      tau_relprof, znoise_int, znoise_ext, magnetic_zaver_range, &
      lbx_ext_global,lby_ext_global,lbz_ext_global, dipole_moment, &
      lax_ext_global,lay_ext_global,laz_ext_global, eta_jump2, &
      sheet_position,sheet_thickness,sheet_hyp,ll_sh,mm_sh, &
      source_zav,nzav,indzav,izav_start, k1hel, k2hel, lbb_sph_as_aux
!
! Run parameters
!
  real :: eta=0.0, eta1=0.0, eta_hyper2=0.0, eta_hyper3=0.0
  real :: eta_tdep_exponent=0.0, eta_tdep_t0=0.0, eta_tdep_toffset=0.0
  real :: eta_hyper3_mesh=5.0, eta_spitzer=0., eta_anom=0.0,& 
          eta_anom_thresh=0.0
  real :: eta_int=0.0, eta_ext=0.0, wresistivity=0.01, eta_xy_max=1.0
  real :: height_eta=0.0, eta_out=0.0, eta_cspeed=0.5
  real :: tau_aa_exterior=0.0, tauAD=0.0, alev=1.0
  real :: sigma_ratio=1.0, eta_z0=1.0, eta_z1=1.0
  real :: eta_xwidth=0.0, eta_ywidth=0.0, eta_zwidth=0.0, eta_zwidth2=0.0
  real :: eta_rwidth=0.0
  real :: eta_width_shock=0.0, eta_zshock=1.0, eta_jump_shock=1.0
  real :: eta_xwidth0=0.0, eta_xwidth1=0.0, eta_rwidth0=0.0, eta_rwidth1=0.0
  real :: eta_xshock=1.0
  real :: eta_x0=1.0, eta_x1=1.0, eta_y0=1.0, eta_y1=1.0
  real :: eta_r0=1.0, eta_r1=1.0
  real :: alphaSSm=0.0, J_ext_quench=0.0, B2_diamag=0.0
  real :: k1_ff=1.0, ampl_ff=1.0, swirl=1.0
  real :: k1x_ff=1.0, k1y_ff=1.0, k1z_ff=1.0
  real :: inertial_length=0.0, linertial_2
  real :: forcing_continuous_aa_phasefact=1.0
  real :: forcing_continuous_aa_amplfact=1.0, ampl_fcont_aa=1.0
  real :: LLambda_aa=0.0, vcrit_anom=1.0
  real :: numag=0.0, B0_magfric=1.0, ekman_friction_aa=0.0
  real :: gamma_epspb=2.4, exp_epspb, ncr_quench=0.
  real :: ampl_eta_uz=0.0
  real :: no_ohmic_heat_z0=1.0, no_ohmic_heat_zwidth=0.0
  real, target :: betamin_jxb = 0.0
  real, dimension(mx,my) :: eta_xy
  real, dimension(mx,my,3) :: geta_xy
  real, dimension(nx,ny,nz,3) :: A_relprof
  real, dimension(mz) :: coskz,sinkz,eta_z,geta_z
  real, dimension(mx) :: eta_x,geta_x
  real, dimension(my) :: eta_y,geta_y
  real, dimension(nx) :: eta_r
  real, dimension(nx,3) :: geta_r
  real, dimension(nx) :: va2max_beta=1.
  real, dimension(nz) :: clight2_zdep=1.
  logical :: lfreeze_aint=.false., lfreeze_aext=.false.
  logical :: lweyl_gauge=.false., ladvective_gauge=.false.
  logical :: lupw_aa=.false., ladvective_gauge2=.false.
  logical :: lcalc_aameanz=.false.,lcalc_aamean
  equivalence (lcalc_aameanz,lcalc_aamean)     ! for compatibility
  logical :: lforcing_cont_aa=.false.
  integer :: iforcing_cont_aa=0
  logical :: lelectron_inertia=.false.
  logical :: lkinematic=.false.
  logical :: lignore_Bext_in_b2=.false., luse_Bext_in_b2=.true.
  logical :: lmean_friction=.false.,llocal_friction=.false.
  logical :: lambipolar_strong_coupling=.false.
  logical :: lhalox=.false., lno_ohmic_heat_bound_z=.false.
  logical :: lrun_initaa=.false.,lmagneto_friction=.false.
  logical :: limplicit_resistivity=.false.
  logical :: lncr_correlated=.false., lncr_anticorrelated=.false.
  logical :: lpropagate_borderaa=.true.
  logical :: lremove_meanaz=.false., lremove_meanax=.false., &
             lremove_meanaxy=.false.,lremove_meanaxz=.false.
  logical :: ladd_efield=.false.
  logical :: lsld_bb=.false.
  logical :: lA_relprof_global=.false.
  logical :: lmagnetic_slope_limited=.false.
  logical :: lboris_correction=.false.
  logical :: lnoinduction=.false.
  logical :: lkeplerian_gauge=.false.
  logical :: lremove_volume_average=.false.
  logical :: lrhs_max=.false.
  real :: h_sld_magn=2.0,nlf_sld_magn=1.0,fac_sld_magn=1.0
  real :: ampl_efield=0.
  real :: w_sldchar_mag=1.
  real :: rhoref=impossible, rhoref1
  character (len=labellen) :: A_relaxprofile='0,coskz,0'
  character (len=labellen) :: zdep_profile='fs'
  character (len=labellen) :: ydep_profile='two-step'
  character (len=labellen) :: xdep_profile='two-step'
  character (len=labellen) :: rdep_profile='two-step'
  character (len=labellen) :: eta_xy_profile='schnack89'
  character (len=labellen) :: iforcing_continuous_aa='fixed_swirl'
  character (len=labellen) :: ambipolar_diffusion='constant'
  character (len=labellen) :: div_sld_magn='2nd'
!
  namelist /magnetic_run_pars/ &
      eta, eta1, eta_hyper2, eta_hyper3, eta_anom, eta_anom_thresh, &
      B_ext, B0_ext, t_bext, t0_bext, J_ext, &
      J_ext_quench, omega_Bz_ext, nu_ni, hall_term, Hhall, battery_term, &
      ihall_term, hall_tdep_t0, hall_tdep_exponent, hall_zdep_exponent, &
      eta_hyper3_mesh, eta_tdep_exponent, eta_tdep_t0, &
      eta_tdep_toffset, lresi_eta_tdep_t0_norm, &
      tau_aa_exterior, tauAD, kx_aa, ky_aa, kz_aa, lcalc_aamean,lohmic_heat, &
      lforcing_cont_aa, lforcing_cont_aa_local, iforcing_continuous_aa, &
      forcing_continuous_aa_phasefact, forcing_continuous_aa_amplfact, k1_ff, &
      ampl_ff, swirl, radius, epsilonaa, k1x_ff, k1y_ff, k1z_ff, &
      center1_x, center1_y, center1_z, lcheck_positive_va2, &
      lmean_friction, llocal_friction, LLambda_aa, bthresh, bthresh_per_brms, &
      iresistivity, lweyl_gauge, ladvective_gauge, ladvective_gauge2, lupw_aa, &
      alphaSSm,eta_int, eta_ext, eta_shock, eta_va,eta_j, eta_j2, eta_jrho, &
      eta_min, wresistivity, eta_xy_max, rhomin_jxb, va2max_jxb, va2max_boris, &
      va_min, cmin,va2power_jxb, llorentzforce, linduction, ldiamagnetism, &
      B2_diamag, reinitialize_aa, rescale_aa, initaa, amplaa, lcovariant_magnetic, &
      lB_ext_pot, D_smag, brms_target, rescaling_fraction, lfreeze_aint, &
      lfreeze_aext, sigma_ratio, zdep_profile, ydep_profile, xdep_profile, &
      rdep_profile, height_eta, eta_out, &
      eta_xwidth, eta_ywidth, eta_zwidth, eta_rwidth, &
      eta_zwidth2, eta_xwidth0, eta_xwidth1, eta_rwidth0, eta_rwidth1, &
      eta_z0, eta_z1, eta_y0, eta_y1, eta_x0, eta_x1, eta_r0, eta_r1, &
      eta1_aniso_ratio, eta1_aniso, eta1_aniso_r, eta1_aniso_d, alp_aniso, quench_aniso, &
      eta_aniso_BB, &
      eta_spitzer, borderaa, ljj_as_comaux, lsmooth_jj, &
      eta_aniso_hyper3, lelectron_inertia, inertial_length, &
      lbext_curvilinear, lbb_as_aux, lbb_as_comaux, lB_ext_in_comaux, ljj_as_aux, &
      lkinematic, lbbt_as_aux, ljjt_as_aux, lua_as_aux, ljxb_as_aux, &
      lneutralion_heat, lreset_aa, daareset, eta_shock2, &
      lignore_Bext_in_b2, luse_Bext_in_b2, ampl_fcont_aa, &
      lhalox, vcrit_anom, eta_jump, eta_jump2, lrun_initaa, two_step_factor, &
      magnetic_xaver_range, A_relaxprofile, tau_relprof, amp_relprof, &
      k_relprof,lmagneto_friction,numag, magnetic_zaver_range,&
      lncr_correlated, lncr_anticorrelated, ncr_quench, B0_magfric, ekman_friction_aa, &
      lbx_ext_global,lby_ext_global,lbz_ext_global, &
      lax_ext_global,lay_ext_global,laz_ext_global, &
      limplicit_resistivity,ambipolar_diffusion, betamin_jxb, gamma_epspb, &
      lpropagate_borderaa, lremove_meanaz, lremove_meanax, lremove_meanaxy, lremove_meanaxz, &
      eta_jump_shock, eta_zshock, &
      eta_width_shock, eta_xshock, ladd_global_field, eta_power_x, eta_power_z, & 
      ladd_efield,ampl_efield, h_sld_magn,w_sldchar_mag, lsld_bb, eta_cspeed, &
      lboris_correction,lkeplerian_gauge,lremove_volume_average, &
      rhoref, lambipolar_strong_coupling,letasmag_as_aux,Pm_smag1, &
      ampl_eta_uz, lalfven_as_aux, lno_ohmic_heat_bound_z, &
      no_ohmic_heat_z0, no_ohmic_heat_zwidth, alev, lrhs_max, &
      lnoinduction, lA_relprof_global, nlf_sld_magn, fac_sld_magn, div_sld_magn, &
      lbb_sph_as_aux
!
! Diagnostic variables (need to be consistent with reset list below)
!
  integer :: idiag_eta_tdep=0   ! DIAG_DOC: $t$-dependent $\eta$
  integer :: idiag_ab_int=0     ! DIAG_DOC: $\int\Av\cdot\Bv\;dV$
  integer :: idiag_jb_int=0     ! DIAG_DOC: $\int\jv\cdot\Bv\;dV$
  integer :: idiag_b2tm=0       ! DIAG_DOC: $\left<\bv(t)\cdot\int_0^t\bv(t')
                                ! DIAG_DOC:   dt'\right>$
  integer :: idiag_bjtm=0       ! DIAG_DOC: $\left<\bv(t)\cdot\int_0^t\jv(t')
                                ! DIAG_DOC:   dt'\right>$
  integer :: idiag_jbtm=0       ! DIAG_DOC: $\left<\jv(t)\cdot\int_0^t\bv(t')
                                ! DIAG_DOC:   dt'\right>$
  integer :: idiag_b2ruzm=0     ! DIAG_DOC: $\left<\Bv^2\rho u_z\right>$
  integer :: idiag_b2uzm=0      ! DIAG_DOC: $\left<\Bv^2u_z\right>$
  integer :: idiag_ubbzm=0      ! DIAG_DOC: $\left<(\uv\cdot\Bv)B_z\right>$
  integer :: idiag_b1m=0        ! DIAG_DOC: $\left<|\Bv|\right>$
  integer :: idiag_b2m=0        ! DIAG_DOC: $\left<\Bv^2\right>$
  integer :: idiag_EEM=0        ! DIAG_DOC: $\left<\Bv^2\right>/2$
  integer :: idiag_b4m=0        ! DIAG_DOC: $\left<\Bv^4\right>$
  integer :: idiag_bm2=0        ! DIAG_DOC: $\max(\Bv^2)$
  integer :: idiag_j2m=0        ! DIAG_DOC: $\left<\jv^2\right>$
  integer :: idiag_jm2=0        ! DIAG_DOC: $\max(\jv^2)$
  integer :: idiag_abm=0        ! DIAG_DOC: $\left<\Av\cdot\Bv\right>$
  integer :: idiag_abumx=0      ! DIAG_DOC: $\left<u_x\Av\cdot\Bv\right>$
  integer :: idiag_abumy=0      ! DIAG_DOC: $\left<u_y\Av\cdot\Bv\right>$
  integer :: idiag_abumz=0      ! DIAG_DOC: $\left<u_z\Av\cdot\Bv\right>$
  integer :: idiag_abmh=0       ! DIAG_DOC: $\left<\Av\cdot\Bv\right>$ (temp)
  integer :: idiag_abmn=0       ! DIAG_DOC: $\left<\Av\cdot\Bv\right>$ (north)
  integer :: idiag_abms=0       ! DIAG_DOC: $\left<\Av\cdot\Bv\right>$ (south)
  integer :: idiag_abrms=0      ! DIAG_DOC: $\left<(\Av\cdot\Bv)^2\right>^{1/2}$
  integer :: idiag_jbrms=0      ! DIAG_DOC: $\left<(\jv\cdot\Bv)^2\right>^{1/2}$
  integer :: idiag_ajm=0        ! DIAG_DOC: $\left<\jv\cdot\Av\right>$
  integer :: idiag_jbm=0        ! DIAG_DOC: $\left<\jv\cdot\Bv\right>$
  integer :: idiag_hjbm=0       ! DIAG_DOC:
  integer :: idiag_jbmh=0       ! DIAG_DOC: $\left<\Jv\cdot\Bv\right>$ (temp)
  integer :: idiag_jbmn=0       ! DIAG_DOC: $\left<\Jv\cdot\Bv\right>$ (north)
  integer :: idiag_jbms=0       ! DIAG_DOC: $\left<\Jv\cdot\Bv\right>$ (south)
  integer :: idiag_ubm=0        ! DIAG_DOC: $\left<\uv\cdot\Bv\right>$
  integer :: idiag_dubrms=0     ! DIAG_DOC: $\left<(\uv-\Bv)^2\right>^{1/2}$
  integer :: idiag_dobrms=0     ! DIAG_DOC: $\left<(\boldsymbol{\omega}-\Bv)^2
                                ! DIAG_DOC: \right>^{1/2}$
  integer :: idiag_uxbxm=0      ! DIAG_DOC: $\left<u_xB_x\right>$
  integer :: idiag_uybxm=0      ! DIAG_DOC: $\left<u_yB_x\right>$
  integer :: idiag_uzbxm=0      ! DIAG_DOC: $\left<u_zB_x\right>$
  integer :: idiag_uxbym=0      ! DIAG_DOC: $\left<u_xB_y\right>$
  integer :: idiag_uybym=0      ! DIAG_DOC: $\left<u_yB_y\right>$
  integer :: idiag_uzbym=0      ! DIAG_DOC: $\left<u_zB_y\right>$
  integer :: idiag_uxbzm=0      ! DIAG_DOC: $\left<u_xB_z\right>$
  integer :: idiag_uybzm=0      ! DIAG_DOC: $\left<u_yB_z\right>$
  integer :: idiag_uzbzm=0      ! DIAG_DOC: $\left<u_zB_z\right>$
  integer :: idiag_uxjxm=0      ! DIAG_DOC: $\left<u_xJ_x\right>$
  integer :: idiag_uxjym=0      ! DIAG_DOC: $\left<u_xJ_y\right>$
  integer :: idiag_uxjzm=0      ! DIAG_DOC: $\left<u_xJ_z\right>$
  integer :: idiag_uyjxm=0      ! DIAG_DOC: $\left<u_yJ_x\right>$
  integer :: idiag_uyjym=0      ! DIAG_DOC: $\left<u_yJ_y\right>$
  integer :: idiag_uyjzm=0      ! DIAG_DOC: $\left<u_yJ_z\right>$
  integer :: idiag_uzjxm=0      ! DIAG_DOC: $\left<u_zJ_x\right>$
  integer :: idiag_uzjym=0      ! DIAG_DOC: $\left<u_zJ_y\right>$
  integer :: idiag_uzjzm=0      ! DIAG_DOC: $\left<u_zJ_z\right>$
  integer :: idiag_cosubm=0     ! DIAG_DOC: $\left<\Uv\cdot\Bv/(|\Uv|\,|\Bv|)
                                ! DIAG_DOC: \right>$
  integer :: idiag_jxbxm=0      ! DIAG_DOC: $\left<j_xB_x\right>$
  integer :: idiag_jybxm=0      ! DIAG_DOC: $\left<j_yB_x\right>$
  integer :: idiag_jzbxm=0      ! DIAG_DOC: $\left<j_zB_x\right>$
  integer :: idiag_jxbym=0      ! DIAG_DOC: $\left<j_xB_y\right>$
  integer :: idiag_jybym=0      ! DIAG_DOC: $\left<j_yB_y\right>$
  integer :: idiag_jzbym=0      ! DIAG_DOC: $\left<j_zB_y\right>$
  integer :: idiag_jxbzm=0      ! DIAG_DOC: $\left<j_xB_z\right>$
  integer :: idiag_jybzm=0      ! DIAG_DOC: $\left<j_yB_z\right>$
  integer :: idiag_jzbzm=0      ! DIAG_DOC: $\left<j_zB_z\right>$

  integer :: idiag_uam=0        ! DIAG_DOC: $\left<\uv\cdot\Av\right>$
  integer :: idiag_ujm=0        ! DIAG_DOC: $\left<\uv\cdot\Jv\right>$
  integer :: idiag_fbm=0        ! DIAG_DOC: $\left<\fv\cdot\Bv\right>$
  integer :: idiag_fxbxm=0      ! DIAG_DOC: $\left<f_x B_x\right>$
  integer :: idiag_epsM=0       ! DIAG_DOC: $\left<\eta\mu_0\jv^2\right>$
  integer :: idiag_epsAD=0      ! DIAG_DOC: $\left<\rho^{-1} t_{\rm AD}
                                ! DIAG_DOC: (\vec{J}\times\vec{B})^2\right>$
                                ! DIAG_DOC: (heating by ion-neutrals friction)
  integer :: idiag_bxpt=0       ! DIAG_DOC: $B_x(x_1,y_1,z_1,t)$
  integer :: idiag_bypt=0       ! DIAG_DOC: $B_y(x_1,y_1,z_1,t)$
  integer :: idiag_bzpt=0       ! DIAG_DOC: $B_z(x_1,y_1,z_1,t)$
  integer :: idiag_jxpt=0       ! DIAG_DOC: $J_x(x_1,y_1,z_1,t)$
  integer :: idiag_jypt=0       ! DIAG_DOC: $J_y(x_1,y_1,z_1,t)$
  integer :: idiag_jzpt=0       ! DIAG_DOC: $J_z(x_1,y_1,z_1,t)$
  integer :: idiag_Expt=0       ! DIAG_DOC: ${\cal E}_x(x_1,y_1,z_1,t)$
  integer :: idiag_Eypt=0       ! DIAG_DOC: ${\cal E}_y(x_1,y_1,z_1,t)$
  integer :: idiag_Ezpt=0       ! DIAG_DOC: ${\cal E}_z(x_1,y_1,z_1,t)$
  integer :: idiag_axpt=0       ! DIAG_DOC: $A_x(x_1,y_1,z_1,t)$
  integer :: idiag_aypt=0       ! DIAG_DOC: $A_y(x_1,y_1,z_1,t)$
  integer :: idiag_azpt=0       ! DIAG_DOC: $A_z(x_1,y_1,z_1,t)$
  integer :: idiag_bxp2=0       ! DIAG_DOC: $B_x(x_2,y_2,z_2,t)$
  integer :: idiag_byp2=0       ! DIAG_DOC: $B_y(x_2,y_2,z_2,t)$
  integer :: idiag_bzp2=0       ! DIAG_DOC: $B_z(x_2,y_2,z_2,t)$
  integer :: idiag_jxp2=0       ! DIAG_DOC: $J_x(x_2,y_2,z_2,t)$
  integer :: idiag_jyp2=0       ! DIAG_DOC: $J_y(x_2,y_2,z_2,t)$
  integer :: idiag_jzp2=0       ! DIAG_DOC: $J_z(x_2,y_2,z_2,t)$
  integer :: idiag_Exp2=0       ! DIAG_DOC: ${\cal E}_x(x_2,y_2,z_2,t)$
  integer :: idiag_Eyp2=0       ! DIAG_DOC: ${\cal E}_y(x_2,y_2,z_2,t)$
  integer :: idiag_Ezp2=0       ! DIAG_DOC: ${\cal E}_z(x_2,y_2,z_2,t)$
  integer :: idiag_axp2=0       ! DIAG_DOC: $A_x(x_2,y_2,z_2,t)$
  integer :: idiag_ayp2=0       ! DIAG_DOC: $A_y(x_2,y_2,z_2,t)$
  integer :: idiag_azp2=0       ! DIAG_DOC: $A_z(x_2,y_2,z_2,t)$
  integer :: idiag_epsM_LES=0   ! DIAG_DOC:
  integer :: idiag_aybym2=0     ! DIAG_DOC:
  integer :: idiag_exaym2=0     ! DIAG_DOC:
  integer :: idiag_exabot=0     ! DIAG_DOC: $\int\Ev\times\Av\,dS|_{\rm bot}$
  integer :: idiag_exatop=0     ! DIAG_DOC: $\int\Ev\times\Av\,dS|_{\rm top}$
  integer :: idiag_exjm2=0      ! DIAG_DOC:
  integer :: idiag_emag=0       ! DIAG_DOC: $\int_V{1\over2\mu_0}\Bv^2\, dV$
  integer :: idiag_brms=0       ! DIAG_DOC: $\left<\Bv^2\right>^{1/2}$
  integer :: idiag_bfrms=0      ! DIAG_DOC: $\left<{\Bv'}^2\right>^{1/2}$
  integer :: idiag_bf2m=0       ! DIAG_DOC: $\left<{\Bv'}^2\right>$
  integer :: idiag_bf4m=0       ! DIAG_DOC: $\left<{\Bv'}^4\right>$
  integer :: idiag_bmax=0       ! DIAG_DOC: $\max(|\Bv|)$
  integer :: idiag_bxmin=0      ! DIAG_DOC: $\min(|B_x|)$
  integer :: idiag_bymin=0      ! DIAG_DOC: $\min(|B_y|)$
  integer :: idiag_bzmin=0      ! DIAG_DOC: $\min(|B_z|)$
  integer :: idiag_bxmax=0      ! DIAG_DOC: $\max(|B_x|)$
  integer :: idiag_bymax=0      ! DIAG_DOC: $\max(|B_y|)$
  integer :: idiag_bzmax=0      ! DIAG_DOC: $\max(|B_z|)$
  integer :: idiag_bbxmax=0     ! DIAG_DOC: $\max(|B_x|) excluding Bv_{ext}$
  integer :: idiag_bbymax=0     ! DIAG_DOC: $\max(|B_y|) excluding Bv_{ext}$
  integer :: idiag_bbzmax=0     ! DIAG_DOC: $\max(|B_z|) excluding Bv_{ext}$
  integer :: idiag_jxmax=0      ! DIAG_DOC: $\max(|jv_x|)$
  integer :: idiag_jymax=0      ! DIAG_DOC: $\max(|jv_y|)$
  integer :: idiag_jzmax=0      ! DIAG_DOC: $\max(|jv_z|)$
  integer :: idiag_jrms=0       ! DIAG_DOC: $\left<\jv^2\right>^{1/2}$
  integer :: idiag_hjrms=0      ! DIAG_DOC: $\left<\jv^2\right>^{1/2}$
  integer :: idiag_jmax=0       ! DIAG_DOC: $\max(|\jv|)$
  integer :: idiag_vArms=0      ! DIAG_DOC: $\left<\Bv^2/\varrho\right>^{1/2}$
  integer :: idiag_vAmax=0      ! DIAG_DOC: $\max(\Bv^2/\varrho)^{1/2}$
  integer :: idiag_dtb=0        ! DIAG_DOC: $\delta t / [c_{\delta t}\,\delta x
                                ! DIAG_DOC:   /v_{\rm A,max}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   Alfv{\'e}n time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_dteta=0      ! DIAG_DOC: $\delta t/[c_{\delta t,{\rm v}}\,
                                ! DIAG_DOC:   \delta x^2/\eta_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   resistive time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_dteta3=0     ! DIAG_DOC: $\delta t/[c_{\delta t,{\rm v3}}\,
                                ! DIAG_DOC:   \delta x^6/\eta^{\rm hyper}_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   hyper resistive time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_dtHr=0       ! DIAG_DOC:
  integer :: idiag_dtFr=0       ! DIAG_DOC:
  integer :: idiag_dtBr=0       ! DIAG_DOC:
  integer :: idiag_axm=0        ! DIAG_DOC:
  integer :: idiag_aym=0        ! DIAG_DOC:
  integer :: idiag_azm=0        ! DIAG_DOC:
  integer :: idiag_a2m=0        ! DIAG_DOC: $\left<\Av^2\right>$
  integer :: idiag_arms=0       ! DIAG_DOC: $\left<\Av^2\right>^{1/2}$
  integer :: idiag_amax=0       ! DIAG_DOC: $\max(|\Av|)$
  integer :: idiag_divarms=0    ! DIAG_DOC: $\langle(\nabla\cdot\Av)^2\rangle^{1/2}$
  integer :: idiag_beta1m=0     ! DIAG_DOC: $\left<\Bv^2/(2\mu_0 p)\right>$
                                ! DIAG_DOC:   \quad(mean inverse plasma beta)
  integer :: idiag_beta1max=0   ! DIAG_DOC: $\max[\Bv^2/(2\mu_0 p)]$
                                ! DIAG_DOC:   \quad(maximum inverse plasma beta)
  integer :: idiag_betam = 0    ! DIAG_DOC: $\langle\beta\rangle$
  integer :: idiag_betamax = 0  ! DIAG_DOC: $\max\beta$
  integer :: idiag_betamin = 0  ! DIAG_DOC: $\min\beta$
  integer :: idiag_bxm=0        ! DIAG_DOC: $\left<B_x\right>$
  integer :: idiag_bym=0        ! DIAG_DOC: $\left<B_y\right>$
  integer :: idiag_bzm=0        ! DIAG_DOC: $\left<B_z\right>$
  integer :: idiag_bxbym=0      ! DIAG_DOC: $\left<B_x B_y\right>$
  integer :: idiag_bxbzm=0      ! DIAG_DOC:
  integer :: idiag_bybzm=0      ! DIAG_DOC:
  integer :: idiag_djuidjbim=0  ! DIAG_DOC:
  integer :: idiag_bmx=0        ! DIAG_DOC: $\left<\left<\Bv\right>_{yz}^2
                                ! DIAG_DOC:   \right>^{1/2}$
                                ! DIAG_DOC:   \quad(energy of $yz$-averaged
                                ! DIAG_DOC:   mean field)
  integer :: idiag_bmy=0        ! DIAG_DOC: $\left<\left<\Bv\right>_{xz}^2
                                ! DIAG_DOC:   \right>^{1/2}$
                                ! DIAG_DOC:   \quad(energy of $xz$-averaged
                                ! DIAG_DOC:   mean field)
  integer :: idiag_bmz=0        ! DIAG_DOC: $\left<\left<\Bv\right>_{xy}^2
                                ! DIAG_DOC:   \right>^{1/2}$
                                ! DIAG_DOC:   \quad(energy of $xy$-averaged
                                ! DIAG_DOC:   mean field)
  integer :: idiag_bmzS2=0      ! DIAG_DOC: $\left<\left<\Bv_S\right>_{xy}^2\right>$
  integer :: idiag_bmzA2=0      ! DIAG_DOC: $\left<\left<\Bv_A\right>_{xy}^2\right>$
  integer :: idiag_jmx=0        ! DIAG_DOC: $\left<\left<\Jv\right>_{yz}^2
                                ! DIAG_DOC:   \right>^{1/2}$
                                ! DIAG_DOC:   \quad(energy of $yz$-averaged
                                ! DIAG_DOC:   mean current density)
  integer :: idiag_jmy=0        ! DIAG_DOC: $\left<\left<\Jv\right>_{xz}^2
                                ! DIAG_DOC:   \right>^{1/2}$
                                ! DIAG_DOC:   \quad(energy of $xz$-averaged
                                ! DIAG_DOC:   mean current density)
  integer :: idiag_jmz=0        ! DIAG_DOC: $\left<\left<\Jv\right>_{xy}^2
                                ! DIAG_DOC:   \right>^{1/2}$
                                ! DIAG_DOC:   \quad(energy of $xy$-averaged
                                ! DIAG_DOC:   mean current density)
  integer :: idiag_bmzph=0      ! DIAG_DOC: Phase of a Beltrami field
  integer :: idiag_bmzphe=0     ! DIAG_DOC: Error of phase of a Beltrami field
  integer :: idiag_bsinphz=0    ! DIAG_DOC: sine of phase of a Beltrami field
  integer :: idiag_bcosphz=0    ! DIAG_DOC: cosine of phase of a Beltrami field
  integer :: idiag_emxamz3=0    ! DIAG_DOC: $\left<\left<\Ev\right>_{xy}\times\left<\Av\right>_{xy}
                                ! DIAG_DOC:   \right>$ \quad($xy$-averaged
                                ! DIAG_DOC:   mean field helicity flux)
  integer :: idiag_embmz=0      ! DIAG_DOC: $\left<\left<\Ev\right>_{xy}\cdot\left<\Bv\right>_{xy}
                                ! DIAG_DOC:   \right>$ \quad($xy$-averaged
                                ! DIAG_DOC:   mean field helicity production )
  integer :: idiag_ambmz=0      ! DIAG_DOC: $\left<\left<\Av\right>_{xy}\cdot\left<\Bv\right>_{xy}\right>$
                                ! DIAG_DOC:   \quad (magnetic helicity of $xy$-averaged mean field)
  integer :: idiag_ambmzh=0     ! DIAG_DOC: $\left<\left<\Av\right>_{xy}\cdot\left<\Bv\right>_{xy}\right>$
                                ! DIAG_DOC:   \quad (magnetic helicity of $xy$-averaged mean field, temp)
  integer :: idiag_ambmzn=0     ! DIAG_DOC: $\left<\left<\Av\right>_{xy}\cdot\left<\Bv\right>_{xy}\right>$
                                ! DIAG_DOC:   \quad (magnetic helicity of $xy$-averaged mean field, north)
  integer :: idiag_ambmzs=0     ! DIAG_DOC: $\left<\left<\Av\right>_{xy}\cdot\left<\Bv\right>_{xy}\right>$
                                ! DIAG_DOC:   \quad (magnetic helicity of $xy$-averaged mean field, south)
  integer :: idiag_jmbmz=0      ! DIAG_DOC: $\left<\left<\Jv\right>_{xy}\cdot\left<\Bv\right>_{xy}
                                ! DIAG_DOC:   \right>$ \quad(current helicity
                                ! DIAG_DOC:   of $xy$-averaged mean field)
  integer :: idiag_Rmmz=0       ! DIAG_DOC: $\left<\frac{|\uv\times\Bv|}{|\eta\Jv|}
                                ! DIAG_DOC: \right>_{xy}$
  integer :: idiag_kx_aa=0      ! DIAG_DOC: $k_x$
  integer :: idiag_kmz=0        ! DIAG_DOC: $\left<\left<\Jv\right>_{xy}\cdot\left<\Bv\right>_{xy}\right>/
                                ! DIAG_DOC:  \left<\left<\Bv\right>_{xy}^2\right>$
  integer :: idiag_bx2m=0       ! DIAG_DOC: $\left< B_x^2 \right>$
  integer :: idiag_by2m=0       ! DIAG_DOC: $\left< B_y^2 \right>$
  integer :: idiag_bz2m=0       ! DIAG_DOC: $\left< B_z^2 \right>$
  integer :: idiag_bx4m=0       ! DIAG_DOC: $\left< B_x^4 \right>$
  integer :: idiag_by4m=0       ! DIAG_DOC: $\left< B_y^4 \right>$
  integer :: idiag_bz4m=0       ! DIAG_DOC: $\left< B_z^4 \right>$
  integer :: idiag_jx2m=0       ! DIAG_DOC: $\left< J_x^2 \right>$
  integer :: idiag_jy2m=0       ! DIAG_DOC: $\left< J_y^2 \right>$
  integer :: idiag_jz2m=0       ! DIAG_DOC: $\left< J_z^2 \right>$
  integer :: idiag_jx4m=0       ! DIAG_DOC: $\left< J_x^4 \right>$
  integer :: idiag_jy4m=0       ! DIAG_DOC: $\left< J_y^4 \right>$
  integer :: idiag_jz4m=0       ! DIAG_DOC: $\left< J_z^4 \right>$
  integer :: idiag_uxbm=0       ! DIAG_DOC: $\left<\uv\times\Bv\right>\cdot\Bv_0/B_0^2$
  integer :: idiag_jxbm=0       ! DIAG_DOC: $\left<\jv\times\Bv\right>\cdot\Bv_0/B_0^2$
  integer :: idiag_vmagfricmax=0 ! DIAG_DOC: $\max(1/\nu_{\rm mag}|\jv\times\Bv/\Bv^2|)$
  integer :: idiag_vmagfricrms=0 ! DIAG_DOC: $\left<1/\nu_{\rm mag}|\jv\times\Bv/\Bv^2|^2\right>^{1/2}$
  integer :: idiag_oxuxbm=0     ! DIAG_DOC:
  integer :: idiag_jxbxbm=0     ! DIAG_DOC:
  integer :: idiag_gpxbm=0      ! DIAG_DOC:
  integer :: idiag_uxDxuxbm=0   ! DIAG_DOC:
  integer :: idiag_b3b21m=0     ! DIAG_DOC: $\left<B_3 B_{2,1} \right>$
  integer :: idiag_b3b12m=0     ! DIAG_DOC: $\left<B_3 B_{1,2} \right>$
  integer :: idiag_b1b32m=0     ! DIAG_DOC: $\left<B_1 B_{3,2} \right>$
  integer :: idiag_b1b23m=0     ! DIAG_DOC: $\left<B_1 B_{2,3} \right>$
  integer :: idiag_b2b13m=0     ! DIAG_DOC: $\left<B_2 B_{1,3} \right>$
  integer :: idiag_b2b31m=0     ! DIAG_DOC: $\left<B_2 B_{3,1} \right>$
  integer :: idiag_udotxbm=0    ! DIAG_DOC:
  integer :: idiag_uxbdotm=0    ! DIAG_DOC:
  integer :: idiag_uxbmx=0      ! DIAG_DOC: $\left<(\uv\times\Bv)_x\right>$
  integer :: idiag_uxbmy=0      ! DIAG_DOC: $\left<(\uv\times\Bv)_y\right>$
  integer :: idiag_uxbmz=0      ! DIAG_DOC: $\left<(\uv\times\Bv)_z\right>$
  integer :: idiag_jxbmx=0      ! DIAG_DOC: $\left<(\jv\times\Bv)_x\right>$
  integer :: idiag_jxbmy=0      ! DIAG_DOC: $\left<(\jv\times\Bv)_y\right>$
  integer :: idiag_jxbmz=0      ! DIAG_DOC: $\left<(\jv\times\Bv)_z\right>$
  integer :: idiag_uxbcmx=0     ! DIAG_DOC:
  integer :: idiag_uxbcmy=0     ! DIAG_DOC:
  integer :: idiag_uxbsmx=0     ! DIAG_DOC:
  integer :: idiag_uxbsmy=0     ! DIAG_DOC:
  integer :: idiag_examx=0      ! DIAG_DOC: $\left<\Ev\times\Av\right>|_x$
  integer :: idiag_examy=0      ! DIAG_DOC: $\left<\Ev\times\Av\right>|_y$
  integer :: idiag_examz=0      ! DIAG_DOC: $\left<\Ev\times\Av\right>|_z$
  integer :: idiag_exjmx=0      ! DIAG_DOC: $\left<\Ev\times\Jv\right>|_x$
  integer :: idiag_exjmy=0      ! DIAG_DOC: $\left<\Ev\times\Jv\right>|_y$
  integer :: idiag_exjmz=0      ! DIAG_DOC: $\left<\Ev\times\Jv\right>|_z$
  integer :: idiag_dexbmx=0     ! DIAG_DOC: $\left<\nabla\times\Ev\times\Bv\right>|_x$
  integer :: idiag_dexbmy=0     ! DIAG_DOC: $\left<\nabla\times\Ev\times\Bv\right>|_y$
  integer :: idiag_dexbmz=0     ! DIAG_DOC: $\left<\nabla\times\Ev\times\Bv\right>|_z$
  integer :: idiag_phibmx=0     ! DIAG_DOC: $\left<\phi\Bv\right>|_x$
  integer :: idiag_phibmy=0     ! DIAG_DOC: $\left<\phi\Bv\right>|_y$
  integer :: idiag_phibmz=0     ! DIAG_DOC: $\left<\phi\Bv\right>|_z$
  integer :: idiag_uxjm=0       ! DIAG_DOC:
  integer :: idiag_b2divum=0    ! DIAG_DOC: $\left<\Bv^2\nabla\cdot\uv\right>$
  integer :: idiag_jdel2am=0    ! DIAG_DOC: $\left<\Jv\cdot\nabla^2\Av)\right>$
  integer :: idiag_ujxbm=0      ! DIAG_DOC: $\left<\uv\cdot(\Jv\times\Bv)\right>$
  integer :: idiag_jxbrxm=0     ! DIAG_DOC:
  integer :: idiag_jxbrym=0     ! DIAG_DOC:
  integer :: idiag_jxbrzm=0     ! DIAG_DOC:
  integer :: idiag_jxbrmax=0    ! DIAG_DOC: $\max(|\Jv\times\Bv/\rho|)$
  integer :: idiag_jxbr2m=0     ! DIAG_DOC: $\left<(\Jv\times\Bv/\rho)^2\right>$
  integer :: idiag_uxBrms=0     ! DIAG_DOC:
  integer :: idiag_Bresrms=0    ! DIAG_DOC:
  integer :: idiag_Rmrms=0      ! DIAG_DOC:
  integer :: idiag_jfm=0        ! DIAG_DOC:
  integer :: idiag_brbpmr=0     ! DIAG_DOC:
  integer :: idiag_vA2m=0       ! DIAG_DOC:
  integer :: idiag_b2mr=0       ! DIAG_DOC:
  integer :: idiag_brmr=0       ! DIAG_DOC:
  integer :: idiag_bpmr=0       ! DIAG_DOC:
  integer :: idiag_bzmr=0       ! DIAG_DOC:
  integer :: idiag_armr=0       ! DIAG_DOC:
  integer :: idiag_apmr=0       ! DIAG_DOC:
  integer :: idiag_azmr=0       ! DIAG_DOC:
  integer :: idiag_mflux_x=0    ! DIAG_DOC:
  integer :: idiag_mflux_y=0    ! DIAG_DOC:
  integer :: idiag_mflux_z=0    ! DIAG_DOC:
  integer :: idiag_bmxy_rms=0   ! DIAG_DOC: $\sqrt{[\left<b_x\right>_z(x,y)]^2 +
                                ! DIAG_DOC: [\left<b_y\right>_z(x,y)]^2 +
                                ! DIAG_DOC: [\left<b_z\right>_z(x,y)]^2} $
  integer :: idiag_etasmagm=0   ! DIAG_DOC: Mean of Smagorinsky resistivity
  integer :: idiag_etasmagmin=0 ! DIAG_DOC: Min of Smagorinsky resistivity
  integer :: idiag_etasmagmax=0 ! DIAG_DOC: Max of Smagorinsky resistivity
  integer :: idiag_etavamax=0   ! DIAG_DOC: Max of artificial resistivity
                                ! DIAG_DOC: $\eta\sim v_A$
  integer :: idiag_etajmax=0    ! DIAG_DOC: Max of artificial resistivity
                                ! DIAG_DOC: $\eta\sim J / \sqrt{\rho}$
  integer :: idiag_etaj2max=0   ! DIAG_DOC: Max of artificial resistivity
                                ! DIAG_DOC: $\eta\sim J^2 / \rho$
  integer :: idiag_etajrhomax=0 ! DIAG_DOC: Max of artificial resistivity
                                ! DIAG_DOC: $\eta\sim J / \rho$
  integer :: idiag_etaaniso=0   ! DIAG_DOC: $\eta_1$
  integer :: idiag_cosjbm=0     ! DIAG_DOC: $\left<\Jv\cdot\Bv/(|\Jv|\,|\Bv|)\right>$
  integer :: idiag_coshjbm=0    ! DIAG_DOC:
  integer :: idiag_jparallelm=0 ! DIAG_DOC: Mean value of the component
                                ! DIAG_DOC: of J parallel to B
  integer :: idiag_jperpm=0     ! DIAG_DOC: Mean value of the component
                                ! DIAG_DOC: of J perpendicular to B
  integer :: idiag_hjparallelm=0 ! DIAG_DOC: Mean value of the component
                                ! DIAG_DOC: of $J_{\rm hyper}$ parallel to B
  integer :: idiag_hjperpm=0    ! DIAG_DOC: Mean value of the component
                                ! DIAG_DOC: of $J_{\rm hyper}$ perpendicular to B
  integer :: idiag_brmsn=0,idiag_brmss=0,idiag_brmsh=0
  integer :: idiag_brmsx=0      ! DIAG_DOC: $\left<\Bv^2\right>^{1/2}$ for
                                ! DIAG_DOC: the magnetic_xaver_range
  integer :: idiag_brmsz=0      ! DIAG_DOC: $\left<\Bv^2\right>^{1/2}$ for
                                ! DIAG_DOC: the magnetic_zaver_range
  integer :: idiag_Exmxy=0      ! DIAG_DOC: $\left<{\cal E}_x\right>_{z}$
  integer :: idiag_Eymxy=0      ! DIAG_DOC: $\left<{\cal E}_y\right>_{z}$
  integer :: idiag_Ezmxy=0      ! DIAG_DOC: $\left<{\cal E}_z\right>_{z}$
!
! phi averaged diagnostics given in phiaver.in
!
  integer :: idiag_jxbrmphi=0   ! PHIAVG_DOC:
  integer :: idiag_jxbpmphi=0   ! PHIAVG_DOC:
  integer :: idiag_jxbzmphi=0   ! PHIAVG_DOC:
  integer :: idiag_jbmphi=0     ! PHIAVG_DOC: $\left<\Jv\cdot\Bv\right>_\varphi$
  integer :: idiag_armphi=0     ! PHIAVG_DOC:
  integer :: idiag_apmphi=0     ! PHIAVG_DOC:
  integer :: idiag_azmphi=0     ! PHIAVG_DOC:
  integer :: idiag_brmphi=0     ! PHIAVG_DOC: $\left<B_\varpi\right>_\varphi$
                                ! PHIAVG_DOC: [cyl.\ polar coords
                                ! PHIAVG_DOC:  $(\varpi,\varphi,z)$]
  integer :: idiag_bpmphi=0     ! PHIAVG_DOC: $\left<B_\varphi\right>_\varphi$
  integer :: idiag_bzmphi=0     ! PHIAVG_DOC: $\left<B_z\right>_\varphi$
  ! For the manual: bbmphi      ! PHIAVG_DOC: shorthand for \var{brmphi},
                                ! PHIAVG_DOC: \var{bpmphi} and \var{bzmphi}
                                ! PHIAVG_DOC: together
  ! For the manual: bbsphmphi   ! PHIAVG_DOC: shorthand for \var{brsphmphi},
                                ! PHIAVG_DOC: \var{bthmphi} and \var{bpmphi}
                                ! PHIAVG_DOC: together
  integer :: idiag_b2mphi=0     ! PHIAVG_DOC: $\left<\Bv^2\right>_\varphi$
  integer :: idiag_brsphmphi=0  ! PHIAVG_DOC: $\left<B_r\right>_\varphi$
  integer :: idiag_bthmphi=0    ! PHIAVG_DOC: $\left<B_\vartheta\right>_\varphi$
  integer :: idiag_uxbrmphi=0   ! PHIAVG_DOC:
  integer :: idiag_uxbpmphi=0   ! PHIAVG_DOC:
  integer :: idiag_uxbzmphi=0   ! PHIAVG_DOC:
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_axmz=0       ! XYAVG_DOC: $\left<{\cal A}_x\right>_{xy}$
  integer :: idiag_aymz=0       ! XYAVG_DOC: $\left<{\cal A}_y\right>_{xy}$
  integer :: idiag_azmz=0       ! XYAVG_DOC: $\left<{\cal A}_z\right>_{xy}$
  integer :: idiag_abuxmz=0     ! XYAVG_DOC: $\left<(\Av \cdot \Bv) u_x \right>_{xy}$
  integer :: idiag_abuymz=0     ! XYAVG_DOC: $\left<(\Av \cdot \Bv) u_y \right>_{xy}$
  integer :: idiag_abuzmz=0     ! XYAVG_DOC: $\left<(\Av \cdot \Bv) u_z \right>_{xy}$
  integer :: idiag_uabxmz=0     ! XYAVG_DOC: $\left<(\uv \cdot \Av) B_x \right>_{xy}$
  integer :: idiag_uabymz=0     ! XYAVG_DOC: $\left<(\uv \cdot \Av) B_y \right>_{xy}$
  integer :: idiag_uabzmz=0     ! XYAVG_DOC: $\left<(\uv \cdot \Av) B_z \right>_{xy}$
  integer :: idiag_bbxmz=0      ! XYAVG_DOC: $\left<{\cal B}'_x\right>_{xy}$
  integer :: idiag_bbymz=0      ! XYAVG_DOC: $\left<{\cal B}'_y\right>_{xy}$
  integer :: idiag_bbzmz=0      ! XYAVG_DOC: $\left<{\cal B}'_z\right>_{xy}$
  integer :: idiag_bxmz=0       ! XYAVG_DOC: $\left<{\cal B}_x\right>_{xy}$
  integer :: idiag_bymz=0       ! XYAVG_DOC: $\left<{\cal B}_y\right>_{xy}$
  integer :: idiag_bzmz=0       ! XYAVG_DOC: $\left<{\cal B}_z\right>_{xy}$
  integer :: idiag_jxmz=0       ! XYAVG_DOC: $\left<{\cal J}_x\right>_{xy}$
  integer :: idiag_jymz=0       ! XYAVG_DOC: $\left<{\cal J}_y\right>_{xy}$
  integer :: idiag_jzmz=0       ! XYAVG_DOC: $\left<{\cal J}_z\right>_{xy}$
  integer :: idiag_Exmz=0       ! XYAVG_DOC: $\left<{\cal E}_x\right>_{xy}$
  integer :: idiag_Eymz=0       ! XYAVG_DOC: $\left<{\cal E}_y\right>_{xy}$
  integer :: idiag_Ezmz=0       ! XYAVG_DOC: $\left<{\cal E}_z\right>_{xy}$
  integer :: idiag_bx2mz=0      ! XYAVG_DOC: $\left< B_x^2 \right>_{xy}$
  integer :: idiag_by2mz=0      ! XYAVG_DOC: $\left< B_y^2 \right>_{xy}$
  integer :: idiag_bz2mz=0      ! XYAVG_DOC: $\left< B_z^2 \right>_{xy}$
  integer :: idiag_bx2rmz=0     ! XYAVG_DOC: $\left< B_x^2/\varrho \right>_{xy}$
  integer :: idiag_by2rmz=0     ! XYAVG_DOC: $\left< B_y^2/\varrho \right>_{xy}$
  integer :: idiag_bz2rmz=0     ! XYAVG_DOC: $\left< B_z^2/\varrho \right>_{xy}$
  integer :: idiag_beta1mz=0    ! XYAVG_DOC: $\left< (B^2 / 2\mu_0) / p \right>_{xy}$
  integer :: idiag_betamz = 0   ! XYAVG_DOC: $\langle\beta\rangle_{xy}$
  integer :: idiag_beta2mz = 0  ! XYAVG_DOC: $\langle\beta^2\rangle_{xy}$
  integer :: idiag_jbmz=0       ! XYAVG_DOC: $\left<\Jv\cdot\Bv\right>|_{xy}$
  integer :: idiag_d6abmz=0     ! XYAVG_DOC: $\left<\nabla^6 \Av\cdot\Bv\right>|_{xy}$
  integer :: idiag_d6amz1=0     ! XYAVG_DOC: $\left<\nabla^6 \Av \right>_{xy}|_x$
  integer :: idiag_d6amz2=0     ! XYAVG_DOC: $\left<\nabla^6 \Av \right>_{xy}|_y$
  integer :: idiag_d6amz3=0     ! XYAVG_DOC: $\left<\nabla^6 \Av \right>_{xy}|_z$
  integer :: idiag_abmz=0       ! XYAVG_DOC: $\left<\Av\cdot\Bv\right>|_{xy}$
  integer :: idiag_ubmz=0       ! XYAVG_DOC: $\left<\uv\cdot\Bv\right>|_{xy}$
  integer :: idiag_uamz=0       ! XYAVG_DOC: $\left<\uv\cdot\Av\right>|_{xy}$
  integer :: idiag_bzdivamz=0   ! XYAVG_DOC: $\left<B_z\nabla\cdot\Av\right>|_{xy}$
  integer :: idiag_divamz=0     ! XYAVG_DOC: $\left<\nabla\cdot\Av\right>|_{xy}$
  integer :: idiag_uxbxmz=0     ! XYAVG_DOC: $\left<u_x b_x\right>|_{xy}$
  integer :: idiag_uybxmz=0     ! XYAVG_DOC: $\left<u_y b_x\right>|_{xy}$
  integer :: idiag_uzbxmz=0     ! XYAVG_DOC: $\left<u_z b_x\right>|_{xy}$
  integer :: idiag_uxbymz=0     ! XYAVG_DOC: $\left<u_x b_y\right>|_{xy}$
  integer :: idiag_uybymz=0     ! XYAVG_DOC: $\left<u_y b_y\right>|_{xy}$
  integer :: idiag_uzbymz=0     ! XYAVG_DOC: $\left<u_z b_y\right>|_{xy}$
  integer :: idiag_uxbzmz=0     ! XYAVG_DOC: $\left<u_x b_z\right>|_{xy}$
  integer :: idiag_uybzmz=0     ! XYAVG_DOC: $\left<u_y b_z\right>|_{xy}$
  integer :: idiag_uzbzmz=0     ! XYAVG_DOC: $\left<u_z b_z\right>|_{xy}$
  integer :: idiag_examz1=0     ! XYAVG_DOC: $\left<\Ev\times\Av\right>_{xy}|_x$
  integer :: idiag_examz2=0     ! XYAVG_DOC: $\left<\Ev\times\Av\right>_{xy}|_y$
  integer :: idiag_examz3=0     ! XYAVG_DOC: $\left<\Ev\times\Av\right>_{xy}|_z$
  integer :: idiag_e3xamz1=0    ! XYAVG_DOC: $\left<\Ev_{hyper3}\times\Av\right>_{xy}|_x$
  integer :: idiag_e3xamz2=0    ! XYAVG_DOC: $\left<\Ev_{hyper3}\times\Av\right>_{xy}|_y$
  integer :: idiag_e3xamz3=0    ! XYAVG_DOC: $\left<\Ev_{hyper3}\times\Av\right>_{xy}|_z$
  integer :: idiag_etatotalmz=0 ! XYAVG_DOC: $\left<\eta\right>_{xy}$
  integer :: idiag_bxbymz=0     ! XYAVG_DOC: $\left< B_x B_y \right>_{xy}$
  integer :: idiag_bxbzmz=0     ! XYAVG_DOC: $\left< B_x B_z \right>_{xy}$
  integer :: idiag_bybzmz=0     ! XYAVG_DOC: $\left< B_y B_z \right>_{xy}$
  integer :: idiag_jxbrxmz=0    ! XYAVG_DOC:
  integer :: idiag_jxbrymz=0    ! XYAVG_DOC:
  integer :: idiag_jxbrzmz=0    ! XYAVG_DOC:
  integer :: idiag_b2mz=0       ! XYAVG_DOC: $\left<\Bv^2\right>_{xy}$
  integer :: idiag_bf2mz=0      ! XYAVG_DOC: $\left<\Bv'^2\right>_{xy}$
  integer :: idiag_j2mz=0       ! XYAVG_DOC: $\left<\jv^2\right>_{xy}$
  integer :: idiag_poynzmz=0    ! XYAVG_DOC: Averaged poynting flux in z direction
  integer :: idiag_epsMmz=0     ! XYAVG_DOC: $\left<\eta\mu_0\jv^2\right>_{xy}$
  integer :: idiag_vmagfricmz=0 ! XYAVG_DOC: $\left<1/\nu_{\rm mag}|\jv\times\Bv/\Bv^2|\right>_{xy}$
!
! xz averaged diagnostics given in xzaver.in
!
  integer :: idiag_bxmy=0       ! XZAVG_DOC: $\left< B_x \right>_{xz}$
  integer :: idiag_bymy=0       ! XZAVG_DOC: $\left< B_y \right>_{xz}$
  integer :: idiag_bzmy=0       ! XZAVG_DOC: $\left< B_z \right>_{xz}$
  integer :: idiag_bx2my=0      ! XZAVG_DOC: $\left< B_x^2 \right>_{xz}$
  integer :: idiag_by2my=0      ! XZAVG_DOC: $\left< B_y^2 \right>_{xz}$
  integer :: idiag_bz2my=0      ! XZAVG_DOC: $\left< B_z^2 \right>_{xz}$
  integer :: idiag_bxbymy=0     ! XZAVG_DOC: $\left< B_x B_y \right>_{xz}$
  integer :: idiag_bxbzmy=0     ! XZAVG_DOC: $\left< B_x B_z \right>_{xz}$
  integer :: idiag_bybzmy=0     ! XZAVG_DOC: $\left< B_y B_z \right>_{xz}$
  integer :: idiag_jxbrxmy=0    ! XZAVG_DOC:
  integer :: idiag_jxbrymy=0    ! XZAVG_DOC:
  integer :: idiag_jxbrzmy=0    ! XZAVG_DOC:
!
! yz averaged diagnostics given in yzaver.in
!
  integer :: idiag_b2mx = 0     ! YZAVG_DOC: $\langle B^2\rangle_{yz}$
  integer :: idiag_bxmx=0       ! YZAVG_DOC: $\left< B_x \right>_{yz}$
  integer :: idiag_bymx=0       ! YZAVG_DOC: $\left< B_y \right>_{yz}$
  integer :: idiag_bzmx=0       ! YZAVG_DOC: $\left< B_z \right>_{yz}$
  integer :: idiag_bx2mx=0      ! YZAVG_DOC: $\left< B_x^2 \right>_{yz}$
  integer :: idiag_by2mx=0      ! YZAVG_DOC: $\left< B_y^2 \right>_{yz}$
  integer :: idiag_bz2mx=0      ! YZAVG_DOC: $\left< B_z^2 \right>_{yz}$
  integer :: idiag_bxbymx=0     ! YZAVG_DOC: $\left<B_x B_y\right>_{yz}$
  integer :: idiag_bxbzmx = 0   ! YZAVG_DOC: $\langle B_x B_z\rangle_{yz}$
  integer :: idiag_bybzmx = 0   ! YZAVG_DOC: $\langle B_y B_z\rangle_{yz}$
  integer :: idiag_betamx = 0   ! YZAVG_DOC: $\langle\beta\rangle_{yz}$
  integer :: idiag_beta2mx = 0  ! YZAVG_DOC: $\langle\beta^2\rangle_{yz}$
  integer :: idiag_etatotalmx=0 ! YZAVG_DOC: $\left<\eta\right>_{yz}$
  integer :: idiag_jxbrxmx=0    ! YZAVG_DOC:
  integer :: idiag_jxbrymx=0    ! YZAVG_DOC:
  integer :: idiag_jxbrzmx=0    ! YZAVG_DOC:
!
! y averaged diagnostics given in yaver.in
!
  integer :: idiag_b2mxz=0      ! YAVG_DOC: $\left< \Bv^2 \right>_{y}$
  integer :: idiag_axmxz=0      ! YAVG_DOC: $\left< A_x \right>_{y}$
  integer :: idiag_aymxz=0      ! YAVG_DOC: $\left< A_y \right>_{y}$
  integer :: idiag_azmxz=0      ! YAVG_DOC: $\left< A_z \right>_{y}$
  integer :: idiag_bx1mxz=0     ! YAVG_DOC: $\left<|B_x|\right>_{y}$
  integer :: idiag_by1mxz=0     ! YAVG_DOC: $\left<|B_y|\right>_{y}$
  integer :: idiag_bz1mxz=0     ! YAVG_DOC: $\left<|B_z|\right>_{y}$
  integer :: idiag_bxmxz=0      ! YAVG_DOC: $\left< B_x \right>_{y}$
  integer :: idiag_bymxz=0      ! YAVG_DOC: $\left< B_y \right>_{y}$
  integer :: idiag_bzmxz=0      ! YAVG_DOC: $\left< B_z \right>_{y}$
  integer :: idiag_jxmxz=0      ! YAVG_DOC: $\left< J_x \right>_{y}$
  integer :: idiag_jymxz=0      ! YAVG_DOC: $\left< J_y \right>_{y}$
  integer :: idiag_jzmxz=0      ! YAVG_DOC: $\left< J_z \right>_{y}$
  integer :: idiag_bx2mxz=0     ! YAVG_DOC: $\left< B_x^2 \right>_{y}$
  integer :: idiag_by2mxz=0     ! YAVG_DOC: $\left< B_y^2 \right>_{y}$
  integer :: idiag_bz2mxz=0     ! YAVG_DOC: $\left< B_z^2 \right>_{y}$
  integer :: idiag_bxbymxz=0    ! YAVG_DOC: $\left< B_x B_y \right>_{y}$
  integer :: idiag_bxbzmxz=0    ! YAVG_DOC: $\left< B_x B_z \right>_{y}$
  integer :: idiag_bybzmxz=0    ! YAVG_DOC: $\left< B_y B_z \right>_{y}$
  integer :: idiag_uybxmxz=0    ! YAVG_DOC: $\left< U_y B_x \right>_{y}$
  integer :: idiag_uybzmxz=0    ! YAVG_DOC: $\left< U_y B_z \right>_{y}$
  integer :: idiag_Exmxz=0      ! YAVG_DOC: $\left<{\cal E}_x\right>_{y}$
  integer :: idiag_Eymxz=0      ! YAVG_DOC: $\left<{\cal E}_y\right>_{y}$
  integer :: idiag_Ezmxz=0      ! YAVG_DOC: $\left<{\cal E}_z\right>_{y}$
  integer :: idiag_vAmxz=0      ! YAVG_DOC: $\left<v_A^2\right>_{y}$
!
! z averaged diagnostics given in zaver.in
!
  integer :: idiag_bxmxy=0      ! ZAVG_DOC: $\left< B_x \right>_{z}$
  integer :: idiag_bymxy=0      ! ZAVG_DOC: $\left< B_y \right>_{z}$
  integer :: idiag_bzmxy=0      ! ZAVG_DOC: $\left< B_z \right>_{z}$
  integer :: idiag_jxmxy=0      ! ZAVG_DOC: $\left< J_x \right>_{z}$
  integer :: idiag_jymxy=0      ! ZAVG_DOC: $\left< J_y \right>_{z}$
  integer :: idiag_jzmxy=0      ! ZAVG_DOC: $\left< J_z \right>_{z}$
  integer :: idiag_axmxy=0      ! ZAVG_DOC: $\left< A_x \right>_{z}$
  integer :: idiag_aymxy=0      ! ZAVG_DOC: $\left< A_y \right>_{z}$
  integer :: idiag_azmxy=0      ! ZAVG_DOC: $\left< A_z \right>_{z}$
  integer :: idiag_bx2mxy=0     ! ZAVG_DOC: $\left< B_x^2 \right>_{z}$
  integer :: idiag_by2mxy=0     ! ZAVG_DOC: $\left< B_y^2 \right>_{z}$
  integer :: idiag_bz2mxy=0     ! ZAVG_DOC: $\left< B_z^2 \right>_{z}$
  integer :: idiag_bxbymxy=0    ! ZAVG_DOC: $\left< B_x B_y \right>_{z}$
  integer :: idiag_bxbzmxy=0    ! ZAVG_DOC: $\left< B_x B_z \right>_{z}$
  integer :: idiag_bybzmxy=0    ! ZAVG_DOC: $\left< B_y B_z \right>_{z}$
  integer :: idiag_poynxmxy=0   ! ZAVG_DOC: $\left< \Ev\times\Bv \right>_{x}$
  integer :: idiag_poynymxy=0   ! ZAVG_DOC: $\left< \Ev\times\Bv \right>_{y}$
  integer :: idiag_poynzmxy=0   ! ZAVG_DOC: $\left< \Ev\times\Bv \right>_{z}$
  integer :: idiag_jbmxy=0      ! ZAVG_DOC: $\left< \Jv\cdot\Bv \right>_{z}$
  integer :: idiag_abmxy=0      ! ZAVG_DOC: $\left< \Av\cdot\Bv \right>_{z}$
  integer :: idiag_ubmxy=0      ! ZAVG_DOC: $\left< \Uv\cdot\Bv \right>_{z}$
  integer :: idiag_examxy1=0    ! ZAVG_DOC: $\left< \Ev\times\Av \right>_{z}|_x$
  integer :: idiag_examxy2=0    ! ZAVG_DOC: $\left< \Ev\times\Av \right>_{z}|_y$
  integer :: idiag_examxy3=0    ! ZAVG_DOC: $\left< \Ev\times\Av \right>_{z}|_z$
  integer :: idiag_StokesImxy=0 ! ZAVG_DOC: $\left< \epsilon_{B\perp} \right>_{z}|_z$
  integer :: idiag_StokesQmxy=0 ! ZAVG_DOC: $-\left<\epsilon_{B\perp} \cos2\chi \right>_{z}|_z$
  integer :: idiag_StokesUmxy=0 ! ZAVG_DOC: $-\left<\epsilon_{B\perp} \sin2\chi \right>_{z}|_z$
  integer :: idiag_StokesQ1mxy=0! ZAVG_DOC: $+\left<F\epsilon_{B\perp} \sin2\chi \right>_{z}|_z$
  integer :: idiag_StokesU1mxy=0! ZAVG_DOC: $-\left<F\epsilon_{B\perp} \cos2\chi \right>_{z}|_z$
  integer :: idiag_beta1mxy=0   ! ZAVG_DOC: $\left< \Bv^2/(2\mu_0 p) \right>_{z}|_z$
  integer :: idiag_dbxdxmxy=0
  integer :: idiag_dbxdymxy=0
  integer :: idiag_dbxdzmxy=0
  integer :: idiag_dbydxmxy=0
  integer :: idiag_dbydymxy=0
  integer :: idiag_dbydzmxy=0
  integer :: idiag_dbzdxmxy=0
  integer :: idiag_dbzdymxy=0
  integer :: idiag_dbzdzmxy=0
!
!  Video data.
!
  integer :: ivid_aps=0, ivid_bb=0, ivid_jj=0, ivid_b2=0, ivid_j2=0, ivid_ab=0, &
             ivid_jb=0, ivid_beta1=0, ivid_poynting=0, ivid_bb_sph=0
!
! Module Variables
!
  real, dimension(nx) :: etatotal=0.,eta_smag=0.,Fmax=0.,dAmax=0.,ss0=0., &
                         diffus_eta=0.,diffus_eta2=0.,diffus_eta3=0.,advec_va2=0.
  real, dimension(nx,3) :: fres,uxbb
  real, dimension(nzgrid) :: eta_zgrid=0.0
  real, dimension(mz) :: feta_ztdep=0.0
  real :: eta_shock_jump1=1.0, eta_tdep=0.0, Arms=0.0
  real, dimension(-nghost:nghost,-nghost:nghost,-nghost:nghost) :: kern_jjsmooth
!
  contains
!***********************************************************************
    subroutine register_magnetic
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!
!  01-may-02/wolf: coded
!  15-oct-15/MR: changes for slope-limited diffusion
!  03-apr-20/joern: restructured and fixed slope-limited diffusion
!
      use FArrayManager, only: farray_register_pde,farray_register_auxiliary
      use SharedVariables, only: get_shared_variable
!
      call farray_register_pde('aa',iaa,vector=3)
      iax = iaa; iay = iaa+1; iaz = iaa+2
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',aa $'
          if (nvar == mvar) write(4,*) ',aa'
        else
          write(4,*) ',aa $'
        endif
        write(15,*) 'aa = fltarr(mx,my,mz,3)*one'
      endif
!
! register EE as auxilliary array if asked for.
!
      if (lEE_as_aux) then
        call farray_register_auxiliary('EE',iEE,vector=3)
        iEEx=iEE; iEEy=iEE+1; iEEz=iEE+2
!
!  Writing files for use with IDL
!
        if (lroot) write(4,*) ',ee $'
        if (lroot) write(15,*) 'ee = fltarr(mx,my,mz,3)*one'
      endif
!
      if (any(iresistivity=='eta-slope-limited')) then
        lslope_limit_diff = .true.
        if (dimensionality<3)lisotropic_advection=.true.
        lbb_as_comaux=lsld_bb
        if (isld_char == 0) then
          call farray_register_auxiliary('sld_char',isld_char,communicated=.true.)
          if (lroot) write(15,*) 'sld_char= fltarr(mx,my,mz)*one'
          aux_var(aux_count)=',sld_char'
          if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=trim(aux_var(aux_count))//' $'
          aux_count=aux_count+1
        endif
      endif

!
!  Register nusmag as auxilliary variable
!
      if (letasmag_as_aux.and.any(iresistivity=='smagorinsky')) then
        call farray_register_auxiliary('etasmag',ietasmag,communicated=.true.)
        if (lroot) write(15,*) 'etasmag = fltarr(mx,my,mz)*one'
        aux_var(aux_count)=',etasmag'
        if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=trim(aux_var(aux_count))//' $'
        aux_count=aux_count+1
      endif
!
!  Check if we are solving the relativistic eos equations.
!  In that case we'd need to get lrelativistic_eos from density.
!
      if (ldensity) &
        call get_shared_variable('lrelativistic_eos', &
            lrelativistic_eos, caller='register_magnetic')
!
!  register the mean-field module
!
      if (lmagn_mf) call register_magn_mf
!
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic(f)
!
!  Perform any post-parameter-read initialization
!
!  24-nov-02/tony: dummy routine - nothing to do at present
!  20-may-03/axel: reinitialize_aa added
!  13-jan-11/MR: use subroutine 'register_report_aux' instead of repeated code
!  26-feb-13/axel: reinitialize_aa added
!  21-jan-15/MR: avoided double put_shared_variable for B_ext
!   7-jun-16/MR: modifications in z average removal for Yin-Yang, yet inoperational
!  24-jun-17/MR: moved calculation of clight2_zdep from calc_pencils to initialize
!  28-feb-18/piyali: moved back the calculation of clight2_zdep to calc_pencils to use va2 pencil
!
      use Sub, only: register_report_aux, write_zprof, step, get_smooth_kernel
      use Magnetic_meanfield, only: initialize_magn_mf
      use BorderProfiles, only: request_border_driving
      use FArrayManager
      use SharedVariables, only: put_shared_variable
      use EquationOfState, only: cs0
      use Initcond
      use Forcing, only: n_forcing_cont
      use Yinyang_mpi, only: initialize_zaver_yy
      !use Slices_methods, only: alloc_slice_buffers
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,myl,nycap
      real :: J_ext2, eta_zdep_exponent
!
!  Set ljj_as_comaux=T and get kernels
!   if lsmooth_jj is used
!
      if(lsmooth_jj) then
        ljj_as_comaux=lsmooth_jj
        call get_smooth_kernel(kern_jjsmooth,LGAUSS=.true.)
      endif
!
!  Set initial value for alfven speed in farray if used
!
      if (lalfven_as_aux) f(:,:,:,ialfven)=0.0
!
!  Share lbb_as_comaux with gravitational wave module.
!
      call put_shared_variable('lbb_as_comaux', lbb_as_comaux, caller='initialize_magnetic')
!
!  Share the external magnetic field with module Shear.
!
      if (lmagn_mf.or.lshock .or. leos .or. lspecial) &
        call put_shared_variable('B_ext', B_ext)
!
!  Share the external magnetic field with mean field module.
!
      if (lmagn_mf) &
        call put_shared_variable('B_ext2', B_ext2)
!
!  Share several parameters for Alfven limiter with module Shock.
!
      if (lshock) then
        call put_shared_variable('va2power_jxb', va2power_jxb)
        call put_shared_variable('betamin_jxb', betamin_jxb)
      endif
!
      call put_shared_variable('rhoref', rhoref)
      call put_shared_variable('eta', eta)
!
!  Shear of B_ext,x is not implemented.
!
      if (lshear .and. B_ext(1) /= 0.0) &
        call fatal_error('initialize_magnetic', 'B_ext,x /= 0 with shear is not implemented.')
!
!  Compute mask for x-averaging where x is in magnetic_xaver_range.
!  Normalize such that the average over the full domain
!  gives still unity.
!
      if (l1 == l2) then
        xmask_mag = 1.
      else
        where (      x(l1:l2) >= magnetic_xaver_range(1) &
               .and. x(l1:l2) <= magnetic_xaver_range(2))
          xmask_mag = 1.
        elsewhere
          xmask_mag = 0.
        endwhere
        magnetic_xaver_range(1) = max(magnetic_xaver_range(1), xyz0(1))
        magnetic_xaver_range(2) = min(magnetic_xaver_range(2), xyz1(1))
        if (lspherical_coords) then
          xmask_mag = xmask_mag * (xyz1(1)**3 - xyz0(1)**3) &
              / (magnetic_xaver_range(2)**3 - magnetic_xaver_range(1)**3)
        elseif (lcylindrical_coords) then
          xmask_mag = xmask_mag * (xyz1(1)**2 - xyz0(1)**2)&
              / (magnetic_xaver_range(2)**2 - magnetic_xaver_range(1)**2)
        else
          xmask_mag = xmask_mag*Lxyz(1) &
              / (magnetic_xaver_range(2) - magnetic_xaver_range(1))
        endif
      endif
!
!  Compute mask for z-averaging where z is in magnetic_zaver_range.
!  Normalize such that the average over the full domain
!  gives still unity.
!
      if (n1 == n2) then
        zmask_mag = 1.
      else
        where (z(n1:n2) >= magnetic_zaver_range(1) .and. z(n1:n2) <= magnetic_zaver_range(2))
          zmask_mag = 1.
        elsewhere
          zmask_mag = 0.
        endwhere
        magnetic_zaver_range(1) = max(magnetic_zaver_range(1), xyz0(3))
        magnetic_zaver_range(2) = min(magnetic_zaver_range(2), xyz1(3))
        zmask_mag = zmask_mag*Lxyz(3)/(magnetic_zaver_range(2) - magnetic_zaver_range(1))
      endif
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'xmask_mag=',xmask_mag
        print*,'zmask_mag=',zmask_mag
      endif
!
!  Precalculate 1/mu (moved here from register.f90)
!
      mu01=1./mu0
      mu012=.5*mu01
!
!  Precalculate eta if 1/eta (==eta1) is given instead
!
      if (eta1/=0.0) then
        eta=1./eta1
      endif
!
!  Precalculate 1/inertial_length^2
!
      if (inertial_length/=0.0) then
        linertial_2 = inertial_length**(-2)
      else
        linertial_2 = 0.0
        ! make sure not to use this value by checking that
        ! (inertial_length /= 0.)...
      endif
!
!  Precalculate 1/nu_ni
!
      if (nu_ni/=0.0) then
        lambipolar_diffusion=.true.
        nu_ni1=1./nu_ni
      else
        nu_ni1=0.0
      endif
!
!  calculate B_ext21 = 1/B_ext**2 and the unit vector B1_ext = B_ext/|B_ext|
!  Also calculate B_ext_inv = B_ext/|B_ext|^2
!
      B_ext2=B_ext(1)**2+B_ext(2)**2+B_ext(3)**2
      if (B_ext2/=0.0) then
        B_ext21=1/B_ext2
      else
        B_ext21=1.0
      endif
      B_ext11=sqrt(B_ext21)
      B1_ext=B_ext*B_ext11
      B_ext_inv=B_ext*B_ext21
!
!  Compute exp_epspb=(gamma_epspb+1.)/4.
!  Note that the extra 1/2 factor is because we work with B^2.
!
      exp_epspb=(gamma_epspb+1.)/4.
!
!  Calculate lJ_ext (true if any of the 3 components in true).
!
      J_ext2=J_ext(1)**2+J_ext(2)**2+J_ext(3)**2
      lJ_ext=(J_ext2/=0)
!
!  Reinitialize magnetic field using a small selection of perturbations
!  that were mostly also available as initial conditions.
!
      if (reinitialize_aa) then
        do j=1,ninit
          select case (initaa(j))
          case ('rescale'); f(:,:,:,iax:iaz)=rescale_aa*f(:,:,:,iax:iaz)
          case ('gaussian-noise'); call gaunoise(amplaa(j),f,iax,iaz)
          case ('hor-tube'); call htube(amplaa(j),f,iax,iaz,radius,epsilonaa, &
              center1_x,center1_z)
          case ('cosxcosy'); call cosx_cosy_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),kz_aa(j))
          case ('coswave-Ay-kx'); call coswave(amplaa(j),f,iay,kx=kx_aa(j))
          case ('sinwave-Ax-kz'); call sinwave(amplaa(j),f,iax,kz=kz_aa(j))
          case default
          endselect
        enddo
      endif
!
!  set lrescaling_magnetic=T if linit_aa=T
!
      if (lreset_aa) then
        lrescaling_magnetic=.true.
      endif
!
      if (lfreeze_aint) lfreeze_varint(iax:iaz) = .true.
      if (lfreeze_aext) lfreeze_varext(iax:iaz) = .true.
!
!     Store spatially dependent external field in a global array
!
      if (lbx_ext_global) &
        call farray_register_global("global_bx_ext",iglobal_bx_ext)
      if (lby_ext_global) &
        call farray_register_global("global_by_ext",iglobal_by_ext)
      if (lbz_ext_global) &
        call farray_register_global("global_bz_ext",iglobal_bz_ext)
!
!     Store spatially dependent external potential field in a global array
!
      if (lax_ext_global) &
        call farray_register_global("global_ax_ext",iglobal_ax_ext)
      if (lay_ext_global) &
        call farray_register_global("global_ay_ext",iglobal_ay_ext)
      if (laz_ext_global) &
        call farray_register_global("global_az_ext",iglobal_az_ext)
!
!  Initialize resistivity.
!
      if (iresistivity(1)=='') iresistivity(1)='eta-const'  ! default
      lresi_eta_const=.false.
      lresi_eta_tdep=.false.
      lresi_sqrtrhoeta_const=.false.
      lresi_eta_aniso=.false.
      lresi_hyper2=.false.
      lresi_hyper3=.false.
      lresi_hyper3_polar=.false.
      lresi_hyper3_mesh=.false.
      lresi_hyper3_csmesh=.false.
      lresi_hyper3_strict=.false.
      lresi_hyper3_aniso=.false.
      lresi_eta_shock=.false.
      lresi_eta_shock2=.false.
      lresi_eta_shock_profz=.false.
      lresi_eta_shock_profr=.false.
      lresi_eta_shock_perp=.false.
      lresi_etava=.false.
      lresi_etaj=.false.
      lresi_etaj2=.false.
      lresi_etajrho=.false.
      lresi_smagorinsky=.false.
      lresi_smagorinsky_nusmag=.false.
      lresi_smagorinsky_cross=.false.
      lresi_anomalous=.false.
      lresi_spitzer=.false.
      lresi_cspeed=.false.
      lresi_vAspeed=.false.
      lresi_eta_proptouz=.false.
      lmagnetic_slope_limited=.false.
!
      do i=1,nresi_max
        select case (iresistivity(i))
        case ('eta-const')
          if (lroot) print*, 'resistivity: constant eta'
          lresi_eta_const=.true.
        case ('eta-tdep')
          if (lroot) print*, 'resistivity: time-dependent eta'
          lresi_eta_tdep=.true.
        case ('eta-ztdep')
          if (lroot) print*, 'resistivity: time-dependent eta'
          lresi_eta_tdep=.true.
          lresi_eta_ztdep=.true.
        case ('sqrtrhoeta-const')
          if (lroot) print*, 'resistivity: constant sqrt(rho)*eta'
          lresi_sqrtrhoeta_const=.true.
        case ('eta-aniso')
          if (lroot) print*, 'resistivity: eta-aniso'
          lresi_eta_aniso=.true.
          if (eta1_aniso==impossible.and.eta1_aniso_ratio==impossible) then
            call fatal_error('initialize_magnetic','eta1_aniso and eta1_aniso_ratio undefined')
          elseif (eta1_aniso==impossible.and.eta1_aniso_ratio/=impossible) then
            eta1_aniso=eta1_aniso_ratio*eta
          endif
        case ('etaSS')
          if (lroot) print*, 'resistivity: etaSS (Shakura-Sunyaev)'
          lresi_etaSS=.true.
        case ('hyper2')
          if (lroot) print*, 'resistivity: hyper2'
          lresi_hyper2=.true.
        case ('hyper3')
          if (lroot) print*, 'resistivity: hyper3'
          lresi_hyper3=.true.
        case ('hyper3_cyl','hyper3-cyl','hyper3_sph','hyper3-sph')
          if (lroot) print*, 'resistivity: hyper3 curvilinear'
          lresi_hyper3_polar=.true.
        case ('hyper3-mesh','hyper3_mesh')
          if (lroot) print*, 'resistivity: hyper3 resolution-invariant'
          lresi_hyper3_mesh=.true.
        case ('hyper3-csmesh')
          if (lroot) print*, 'resistivity: hyper3 cspeed resol-invariant'
          lresi_hyper3_csmesh=.true.
        case ('hyper3_strict')
          if (lroot) print*, 'resistivity: strict hyper3 with positive definite heating rate'
          lresi_hyper3_strict=.true.
        case ('xydep')
          if (lroot) print*, 'resistivity: xy-dependent'
          lresi_xydep=.true.
          call eta_xy_dep(eta_xy,geta_xy,eta_xy_profile)
        case ('xdep','eta-xdep')
          if (lroot) print*, 'resistivity: x-dependent'
          lresi_xdep=.true.
          call eta_xdep(eta_x,geta_x,xdep_profile)
        case ('ydep','eta-ydep')
          if (lroot) print*, 'resistivity: y-dependent'
          lresi_ydep=.true.
          call eta_ydep(ydep_profile, my, y, eta_y, geta_y)
        case ('zdep','eta-zdep')
          if (lroot) print*, 'resistivity: z-dependent'
          lresi_zdep=.true.
          call eta_zdep(zdep_profile, mz, z, eta_z, geta_z)
          if (limplicit_resistivity) call eta_zdep(zdep_profile, nzgrid, zgrid, eta_zgrid)
        case ('rdep','eta-rdep')
          if (lroot) print*, 'resistivity: r-dependent'
          lresi_rdep=.true.
        case ('dust')
          if (lroot) print*, 'resistivity: depending on dust density'
          lresi_dust=.true.
        case ('hyper3-aniso')
          if (lroot) print*, 'resistivity: hyper3_aniso'
          lresi_hyper3_aniso=.true.
        case ('shell')
          if (lroot) print*, 'resistivity: shell'
          lresi_shell=.true.
        case ('shock','eta-shock')
          if (lroot) print*, 'resistivity: shock'
          lresi_eta_shock=.true.
          if (.not. lshock) &
              call fatal_error('initialize_magnetic', &
              'shock resistivity, but module setting SHOCK=noshock')
        case ('eta-shock2')
          if (lroot) print*, 'resistivity: shock'
          lresi_eta_shock2=.true.
          if (.not. lshock) &
              call fatal_error('initialize_magnetic', &
              'shock resistivity, but module setting SHOCK=noshock')
        case ('eta-shock-profz')
          if (lroot) print*, 'resistivity: shock with a vertical profile'
          lresi_eta_shock_profz=.true.
          if (.not. lshock) &
              call fatal_error('initialize_magnetic', &
              'shock resistivity, but module setting SHOCK=noshock')
        case ('eta-shock-profr')
          if (lroot) print*, 'resistivity: shock with a radial profile'
          lresi_eta_shock_profr=.true.
          if (.not. lshock) &
              call fatal_error('initialize_magnetic', &
              'shock resistivity, but module setting SHOCK=noshock')
        case ('shock-perp')
          if (lroot) print*, 'resistivity: shock perpendicular to B'
          lresi_eta_shock_perp=.true.
          if (.not. lshock) &
              call fatal_error('initialize_magnetic', &
              'shock-perp resistivity, but module setting SHOCK=noshock')
          if (.not. ldivu_perp) &
              call fatal_error('initialize_magnetic', &
              'shock-perp resistivity, but not ldivu_perp=.true.')
        case ('eta_va')
          if (lroot) print*, 'resistivity: eta_va'
          lresi_etava=.true.
        case ('eta_j')
          if (lroot) print*, 'resistivity: eta_j'
          lresi_etaj=.true.
        case ('eta_j2')
          if (lroot) print*, 'resistivity: eta_j2'
          lresi_etaj2=.true.
          etaj20 = eta_j2 * mu0**2 * dxmax**3 / cs0
        case ('eta_jrho')
          if (lroot) print*, 'resistivity: eta_jrho'
          lresi_etajrho=.true.
        case ('smagorinsky')
          if (lroot) print*, 'resistivity: smagorinsky'
          lresi_smagorinsky=.true.
        case ('smagorinsky-nusmag','smagorinsky_nusmag')
          if (lroot) print*, 'resistivity: smagorinsky_nusmag'
          lresi_smagorinsky_nusmag=.true.
        case ('smagorinsky-cross')
          if (lroot) print*, 'resistivity: smagorinsky_cross'
          lresi_smagorinsky_cross=.true.
        case ('anomalous')
          if (lroot) print*, 'resistivity: anomalous'
          lresi_anomalous=.true.
        case ('spitzer','eta-spitzer')
          if (lroot) print*, 'resistivity: temperature dependent (Spitzer 1969)'
          lresi_spitzer=.true.
        case ('eta-cspeed')
          if (lroot) print*, 'resistivity: sound speed dependent e.g. SN driven ISM'
          lresi_cspeed=.true.
        case ('eta-vAspeed')
          if (lroot) print*, 'resistivity: Alfven speed dependent e.g. SN driven ISM'
          if (.not. lalfven_as_aux) &
              call fatal_error('initialize_magnetic', &
              'Alfven speed dependent resistivity, but not lalfven_as_aux=.true.')
          lresi_vAspeed=.true.
        case ('magfield')
          if (lroot) print*, 'resistivity: magnetic field dependent'
          lresi_magfield=.true.
        case ('eta-proptouz')
          if (lroot) print*, 'resistivity: eta proportional to uz'
          lresi_eta_proptouz=.true.
        case ('eta-slope-limited')
          if (lroot) then
            if (lsld_bb) then
              print*,'resistivity: slope-limited diffusion on bb'
            else
              print*,'resistivity: slope-limited diffusion on aa'
            endif
            print*, 'resistivity: using ',trim(div_sld_magn),' order'
          endif
          lmagnetic_slope_limited=.true.
        case ('none')
          ! do nothing
        case ('')
          ! do nothing
        case default
          if (lroot) print*, 'No such value for iresistivity(',i,'): ', &
              trim(iresistivity(i))
          call fatal_error('initialize_magnetic','')
        endselect
      enddo
!
!  If lresi_eta_ztdep, compute z-dependent fraction here:
!  It is called feta_ztdep, because it works only if lresi_eta_tdep.
!  eta_zdep_exponent = (2/3)*hall_zdep_exponent; see Gourgouliatos+20.
!
      if (lresi_eta_ztdep) then
        if (Hhall==0.) call fatal_error('initialize_magnetic','Hhall=0 not allowed.')
        eta_zdep_exponent=(2./3.)*hall_zdep_exponent
        feta_ztdep=1./(1.-(z-xyz1(3))/Hhall)**eta_zdep_exponent
      endif
!
      if (lyinyang) then
        if (lresi_eta_shock_profz.or.lresi_xydep.or.lresi_ydep.or.lresi_zdep) &
          call fatal_error('initialize_magnetic','y or z dependent profiles not implemented on Yin-Yang grid.')
      endif
      if (lresi_eta_shock_profz .or. lresi_eta_shock_profr) then
        eta_shock_jump1 = eta_shock*(eta_jump_shock-1.)
        if (lresi_eta_shock_profz) &
          call write_zprof('resi_shock', eta_shock + eta_shock_jump1*step(z(n1:n2), eta_zshock, -eta_width_shock))
      endif
!
!  for communication with testfield_general
!
      lresi_dep(1:4) = (/lresi_xdep,lresi_ydep,lresi_zdep,lresi_xydep/)
!
!  If we're timestepping, die or warn if the the resistivity coefficient that
!  corresponds to the chosen resistivity type is not set.
!
      if (lrun) then
        if (lroot) then
          if ((lresi_eta_const.or.lresi_eta_tdep).and.(eta==0.0)) &
              call warning('initialize_magnetic', &
              'Resistivity coefficient eta is zero!')
          if (lresi_sqrtrhoeta_const.and.(eta==0.0)) &
              call warning('initialize_magnetic', &
              'Resistivity coefficient eta is zero!')
          if (lresi_hyper2.and.eta_hyper2==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_hyper2 is zero!')
          if (lresi_hyper3.and.eta_hyper3==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_hyper3 is zero!')
          if (lresi_hyper3_polar.and.eta_hyper3==0.0) &
               call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_hyper3 is zero!')
          if (lresi_hyper3_mesh.and.eta_hyper3_mesh==0.0) &
               call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_hyper3_mesh is zero!')
          if (lresi_hyper3_csmesh.and.eta_hyper3_mesh==0.0) &
               call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_hyper3_mesh is zero!')
          if (lresi_hyper3_strict.and.eta_hyper3==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_hyper3 is zero!')
          if ( (lresi_hyper3_aniso) .and.  &
               ((eta_aniso_hyper3(1)==0.0 .and. nxgrid/=1 ).or. &
                (eta_aniso_hyper3(2)==0.0 .and. nygrid/=1 ).or. &
                (eta_aniso_hyper3(3)==0.0 .and. nzgrid/=1 )) ) &
              call fatal_error('initialize_magnetic', &
              'A resistivity coefficient of eta_aniso_hyper3 is zero!')
          if (lresi_eta_shock.and.eta_shock==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_shock is zero!')
          if (lresi_eta_shock2.and.eta_shock2==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_shock is zero!')
          if (lresi_eta_shock_profz.and.eta_shock==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_shock is zero!')
           if (lresi_eta_shock_profr.and.eta_shock==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_shock is zero!')
          if (lresi_eta_shock_perp.and.eta_shock==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_shock is zero!')
          if ((lresi_etava.or.lresi_vAspeed).and.eta_va==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_va is zero!')
          if (lresi_vAspeed.and.idiag_vArms==0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity requires vArms be included in print.in!')
          if (lresi_etaj .and. eta_j==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_j is zero!')
          if (lresi_etaj2 .and. eta_j2==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_j2 is zero!')
          if (lresi_etajrho .and. eta_jrho==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_jrho is zero!')
          if (lresi_anomalous.and.eta_anom==0.0) &
              call fatal_error('initialize_magnetic', &
              'Resistivity coefficient eta_anom is zero!')
!          if (lentropy .and. lohmic_heat .and. .not. lresi_eta_const) &
!              call fatal_error('initialize_magnetic', &
!            'Resistivity heating only works with regular resistivity!')
          if (lresi_hyper2.and.lresi_hyper3) &
              call warning('initialize_magnetic', &
              '4th & 6th order hyperdiffusion are both set. ' // &
              'Timestep is currently only sensitive to fourth order.')
        endif
!
!  if meanfield theory is invoked, we need to tell the other routines
!  eta is also needed with the chiral fluids procedure.
!
        if (lmagn_mf .or. lspecial) call put_shared_variable('eta',eta)

      endif
!
!  Quenching of \eta by rms of magnetic vector potential?
!
      lquench_eta_aniso=(quench_aniso/=impossible)
      if (.not.lquench_eta_aniso) idiag_etaaniso=0
!
!  precalculating fixed (on timescales larger than tau) vectorpotential
!
      if (tau_relprof/=0.0) then

        if (lyinyang) then
          if (A_relaxprofile/='') &
            call fatal_error('initialize_magnetic', &
            'z dependent relaxation profiles for A not implemented on Yin-Yang grid.')
        endif

        tau_relprof1=1./tau_relprof
        select case (A_relaxprofile)
        case('0,coskz,0')
          A_relprof(:,:,:,1)=0.
          do i=1,nx
            do j=1,ny
              A_relprof(i,j,:,2)=amp_relprof*cos(k_relprof*z(n1:n2))
            enddo
          enddo
          A_relprof(:,:,:,3)=0.
        case('sinkz,coskz,0')
          do i=1,nx
            do j=1,ny
              A_relprof(i,j,:,1)=amp_relprof*sin(k_relprof*z(n1:n2))
              A_relprof(i,j,:,2)=amp_relprof*cos(k_relprof*z(n1:n2))
            enddo
          enddo
          A_relprof(:,:,:,3)=0.
        case('aa_from_global')
        if (lroot) print*, &
             'initialize_mag: Set A_relaxprofile to: ', A_relaxprofile
          if (iglobal_ax_ext/=0 .or. iglobal_ay_ext/=0 .or. iglobal_az_ext/=0) &
          A_relprof(:,:,:,1:3)=amp_relprof*f(l1:l2,m1:m2,n1:n2,iglobal_ax_ext:iglobal_az_ext)
        endselect
      endif
!
!  Write profile (uncomment for debugging)
!
!     if (lfirst_proc_xy) then
!       do n=1,mz
!         write (100+iproc,*) z(n), eta_z(n), geta_z(n)
!       enddo
!     endif
!
!  Share lweyl_gauge
!
      if (lspecial.or.lmagn_mf) &
        call put_shared_variable('lweyl_gauge',lweyl_gauge,caller='initialize_magnetic')
!
!  share eta profile with test-field procedure
!
      if (ltestfield) then
        if (lresi_xdep) then
          call put_shared_variable('eta_x',eta_x,caller='initialize_magnetic')
          call put_shared_variable('geta_x',geta_x)
        endif
        if (lresi_ydep) then
          call put_shared_variable('eta_y',eta_y,caller='initialize_magnetic')
          call put_shared_variable('geta_y',geta_y)
        endif
        if (lresi_zdep) then
          call put_shared_variable('eta_z',eta_z,caller='initialize_magnetic')
          call put_shared_variable('geta_z',geta_z)
        endif
        if (lresi_xydep) then
          call put_shared_variable('eta_xy',eta_xy,caller='initialize_magnetic')
          call put_shared_variable('geta_xy',geta_xy)
        endif
      endif
!
!  Border profile backward compatibility. For a vector, if only the first
!  borderaa is set, then the other components get the same value.
!
      if (lpropagate_borderaa     .and. &
           borderaa(1)/='nothing' .and. &
           borderaa(2)=='nothing' .and. &
           borderaa(3)=='nothing') then
        borderaa(2)=borderaa(1)
        borderaa(3)=borderaa(1)
      endif
!
!  Tell the BorderProfiles module if we intend to use border driving, so
!  that the module can request the right pencils.
!
      do j=1,3
        if (borderaa(j)/='nothing') call request_border_driving(borderaa(j))
      enddo
!
!  Register an extra aux slot for bb if requested (so bb and jj are written
!  to snapshots and can be easily analyzed later). For this to work you
!  must reserve enough auxiliary workspace by setting, for example,
!     ! MAUX CONTRIBUTION 6
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!
!  After a reload, we need to rewrite index.pro, but the auxiliary
!  arrays are already allocated and must not be allocated again.
!
      if (lbb_as_aux .or. lbb_as_comaux) &
        call register_report_aux('bb', ibb, ibx, iby, ibz, communicated=lbb_as_comaux)
      if (ljj_as_aux .or. ljj_as_comaux) &
        call register_report_aux('jj', ijj, ijx, ijy, ijz, communicated=ljj_as_comaux)
!
      if (lbbt_as_aux) then
        call register_report_aux('bbt',ibbt,ibxt,ibyt,ibzt)
        ltime_integrals=.true.
      endif
!
      if (ljjt_as_aux) then
        call register_report_aux('jjt',ijjt,ijxt,ijyt,ijzt)
        ltime_integrals=.true.
      endif
!
      if (lua_as_aux ) call register_report_aux('ua',iua)
      if (ljxb_as_aux) call register_report_aux('jxb',ijxb,ijxbx,ijxby,ijxbz)
!
      if (lbb_sph_as_aux) &
        call register_report_aux('bb_sph', ibb_sph, ibb_sphr, ibb_spht, ibb_sphp)
!
!  Register va as auxilliary array if asked for also requires
!  ! MAUX CONTRIBUTION 1
!  ! COMMUNICATED AUXILIARIES 1
!  in cparam.local
!
      if (lalfven_as_aux) call register_report_aux('alfven',ialfven,communicated=.true.)
!
!  Initialize individual modules, but need to do this only if
!  lmagn_mf is true.
!
      if (lmagn_mf) call initialize_magn_mf(f)
!
      call put_shared_variable('lfrozen_bb_bot',lfrozen_bb_bot,caller='initialize_magnetic')
      call put_shared_variable('lfrozen_bb_top',lfrozen_bb_top)
!
!  Calculate cosz and sinz for calculating the phase of a Beltrami field
!  The choice to use k1_ff may not be optimal, but keep it for now.
!
      if (idiag_bsinphz/=0 .or. idiag_bcosphz/=0 &
          .or. idiag_uxbcmx/=0 .or. idiag_uxbcmy/=0 &
          .or. idiag_uxbsmx/=0 .or. idiag_uxbsmy/=0 ) then
        sinkz=sin(k1_ff*z)
        coskz=cos(k1_ff*z)
      endif
!
!  When adding a magnetic field to a snapshot of a nomagnetic simulation,
!  the code allows only the initialization of the field to zero. This
!  hack allows a init_aa (from start.in) to be read and added to the
!  field upon executing run.csh
!
      if (lread_oldsnap_nomag.and.lrun_initaa) then
        if (lroot) then
          print*,'Adding a magnetic field to a previously '//&
              'non-magnetic simulation. The field is given by initaa=',initaa
        endif
        call init_aa(f)
      endif
!
!  Break if Galilean-invariant advection (fargo) is used without
!  the advective gauge (only in run-time)
!
      if (lrun) then
        if (lfargo_advection.and..not.ladvective_gauge) &
             call fatal_error('initialize_magnetic',&
             'For fargo advection you need the advective gauge. '//&
             'You may want to switch ladvective_gauge=T in magnetic_run_pars')
      endif
!
!  Write constants to disk. In future we may want to deal with this
!  using an include file or another module.
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,'(a,1pd26.16)') 'mu0=',mu0
        close (1)
      endif
!
      if (.not.lforcing_cont) lforcing_cont_aa=.false.
      if (lforcing_cont_aa) then
        iforcing_cont_aa=min(n_forcing_cont,2)
        if (iforcing_cont_aa==0) &
          call fatal_error('initialize_magnetic','no valid continuous forcing available')
      endif
!
      if (lremove_meanaxy) then
        myl=my
        if (lyinyang) then
          call fatal_error('initialize_magnetic','Removal of z average not implmented for Yin-Yang')
          call initialize_zaver_yy(myl,nycap)
        endif
        allocate(aamxy(mx,myl))
      endif
!
      llorentz_rhoref = llorentzforce .and. rhoref/=impossible .and. rhoref>0.
      if (llorentz_rhoref) rhoref1=1./rhoref

      if (ivid_aps/=0) then
        if (.not.lcylindrical_coords.and..not.lspherical_coords.and. &
            .not.(dimensionality==2.and..not.lactive_dimension(3)))  then
          call warning('initialize_magnetic','aa_phi x (axis distance) on slices only implemented for axisymmetric setups')
          ivid_aps=0
        endif
        if (lwrite_slice_xy .and..not.allocated(aps_xy ) ) allocate(aps_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(aps_xz ) ) allocate(aps_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(aps_yz ) ) allocate(aps_yz (ny,nz))
        if (lwrite_slice_xz2.and..not.allocated(aps_xz2) ) allocate(aps_xz2(nx,nz))
      endif

      if (ivid_bb/=0) then
        !call alloc_slice_buffers(bb_xy,bb_xz,bb_yz,bb_xy2,bb_xy3,bb_xy4,bb_xz2)
        if (lwrite_slice_xy .and..not.allocated(bb_xy) ) allocate(bb_xy (nx,ny,3))
        if (lwrite_slice_xz .and..not.allocated(bb_xz) ) allocate(bb_xz (nx,nz,3))
        if (lwrite_slice_yz .and..not.allocated(bb_yz) ) allocate(bb_yz (ny,nz,3))
        if (lwrite_slice_xy2.and..not.allocated(bb_xy2)) allocate(bb_xy2(nx,ny,3))
        if (lwrite_slice_xy3.and..not.allocated(bb_xy3)) allocate(bb_xy3(nx,ny,3))
        if (lwrite_slice_xy4.and..not.allocated(bb_xy4)) allocate(bb_xy4(nx,ny,3))
        if (lwrite_slice_xz2.and..not.allocated(bb_xz2)) allocate(bb_xz2(nx,nz,3))
      endif
      if (ivid_jj/=0) then
        !call alloc_slice_buffers(jj_xy,jj_xz,jj_yz,jj_xy2,jj_xy3,jj_xy4,jj_xz2)
        if (lwrite_slice_xy .and..not.allocated(jj_xy) ) allocate(jj_xy (nx,ny,3))
        if (lwrite_slice_xz .and..not.allocated(jj_xz) ) allocate(jj_xz (nx,nz,3))
        if (lwrite_slice_yz .and..not.allocated(jj_yz) ) allocate(jj_yz (ny,nz,3))
        if (lwrite_slice_xy2.and..not.allocated(jj_xy2)) allocate(jj_xy2(nx,ny,3))
        if (lwrite_slice_xy3.and..not.allocated(jj_xy3)) allocate(jj_xy3(nx,ny,3))
        if (lwrite_slice_xy4.and..not.allocated(jj_xy4)) allocate(jj_xy4(nx,ny,3))
        if (lwrite_slice_xz2.and..not.allocated(jj_xz2)) allocate(jj_xz2(nx,nz,3))
      endif
      if (ivid_b2/=0) then
        !call alloc_slice_buffers(b2_xy,b2_xz,b2_yz,b2_xy2,b2_xy3,b2_xy4,b2_xz2)
        if (lwrite_slice_xy .and..not.allocated(b2_xy) ) allocate(b2_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(b2_xz) ) allocate(b2_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(b2_yz) ) allocate(b2_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(b2_xy2)) allocate(b2_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(b2_xy3)) allocate(b2_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(b2_xy4)) allocate(b2_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(b2_xz2)) allocate(b2_xz2(nx,nz))
      endif
      if (ivid_j2/=0) then
        !call alloc_slice_buffers(j2_xy,j2_xz,j2_yz,j2_xy2,j2_xy3,j2_xy4,j2_xz2)
        if (lwrite_slice_xy .and..not.allocated(j2_xy) ) allocate(j2_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(j2_xz) ) allocate(j2_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(j2_yz) ) allocate(j2_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(j2_xy2)) allocate(j2_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(j2_xy3)) allocate(j2_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(j2_xy4)) allocate(j2_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(j2_xz2)) allocate(j2_xz2(nx,nz))
      endif
      if (ivid_bb_sph/=0) then
        !call alloc_slice_buffers(j2_xy,j2_xz,j2_yz,j2_xy2,j2_xy3,j2_xy4,j2_xz2)
        if (lwrite_slice_xy .and..not.allocated(bb_sph_xy) ) allocate(bb_sph_xy (nx,ny,3))
        if (lwrite_slice_xz .and..not.allocated(bb_sph_xz) ) allocate(bb_sph_xz (nx,nz,3))
        if (lwrite_slice_yz .and..not.allocated(bb_sph_yz) ) allocate(bb_sph_yz (ny,nz,3))
        if (lwrite_slice_xy2.and..not.allocated(bb_sph_xy2)) allocate(bb_sph_xy2(nx,ny,3))
        if (lwrite_slice_xy3.and..not.allocated(bb_sph_xy3)) allocate(bb_sph_xy3(nx,ny,3))
        if (lwrite_slice_xy4.and..not.allocated(bb_sph_xy4)) allocate(bb_sph_xy4(nx,ny,3))
        if (lwrite_slice_xz2.and..not.allocated(bb_sph_xz2)) allocate(bb_sph_xz2(nx,nz,3))
      endif
      if (ivid_ab/=0) then
        !call alloc_slice_buffers(ab_xy,ab_xz,ab_yz,ab_xy2,ab_xy3,ab_xy4,ab_xz2)
        if (lwrite_slice_xy .and..not.allocated(ab_xy) ) allocate(ab_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(ab_xz) ) allocate(ab_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(ab_yz) ) allocate(ab_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(ab_xy2)) allocate(ab_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(ab_xy3)) allocate(ab_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(ab_xy4)) allocate(ab_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(ab_xz2)) allocate(ab_xz2(nx,nz))
      endif
      if (ivid_jb/=0) then
        !call alloc_slice_buffers(jb_xy,jb_xz,jb_yz,jb_xy2,jb_xy3,jb_xy4,jb_xz2)
        if (lwrite_slice_xy .and..not.allocated(jb_xy) ) allocate(jb_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(jb_xz) ) allocate(jb_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(jb_yz) ) allocate(jb_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(jb_xy2)) allocate(jb_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(jb_xy3)) allocate(jb_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(jb_xy4)) allocate(jb_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(jb_xz2)) allocate(jb_xz2(nx,nz))
      endif
      if (ivid_beta1/=0) then
        !call alloc_slice_buffers(beta1_xy,beta1_xz,beta1_yz,beta1_xy2, &
        !                         beta1_xy3,beta1_xy4,beta1_xz2)
        if (lwrite_slice_xy .and..not.allocated(beta1_xy) ) allocate(beta1_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(beta1_xz) ) allocate(beta1_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(beta1_yz) ) allocate(beta1_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(beta1_xy2)) allocate(beta1_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(beta1_xy3)) allocate(beta1_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(beta1_xy4)) allocate(beta1_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(beta1_xz2)) allocate(beta1_xz2(nx,nz))
      endif
      if (ivid_poynting/=0) then
        !call alloc_slice_buffers(poynting_xy,poynting_xz,poynting_yz,poynting_xy2,&
        !                         poynting_xy3,poynting_xy4,poynting_xz2)
        if (lwrite_slice_xy .and..not.allocated(poynting_xy) ) allocate(poynting_xy (nx,ny,3))
        if (lwrite_slice_xz .and..not.allocated(poynting_xz) ) allocate(poynting_xz (nx,nz,3))
        if (lwrite_slice_yz .and..not.allocated(poynting_yz) ) allocate(poynting_yz (ny,nz,3))
        if (lwrite_slice_xy2.and..not.allocated(poynting_xy2)) allocate(poynting_xy2(nx,ny,3))
        if (lwrite_slice_xy3.and..not.allocated(poynting_xy3)) allocate(poynting_xy3(nx,ny,3))
        if (lwrite_slice_xy4.and..not.allocated(poynting_xy4)) allocate(poynting_xy4(nx,ny,3))
        if (lwrite_slice_xz2.and..not.allocated(poynting_xz2)) allocate(poynting_xz2(nx,nz,3))
      endif
!
    endsubroutine initialize_magnetic
!***********************************************************************
    subroutine init_aa(f)
!
!  initialise magnetic field; called from start.f90
!  We have an init parameter (initaa) to stear magnetic i.c. independently.
!
!   7-nov-2001/wolf: coded
!
      use EquationOfState
      use FArrayManager
      use IO, only: input_snap, input_snap_finalize
      use Gravity, only: gravz, z1, z2
      use Initcond
      use Boundcond
      use InitialCondition, only: initial_condition_aa
      use Mpicomm
      use SharedVariables
      use Sub
      use General, only: yin2yang_coors, transform_thph_yy
      use File_io, only: read_zaver
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (mz) :: tmp
      real, dimension (nx,3) :: bb
      real, dimension (nx) :: b2,fact,cs2,lnrho_old,ssold,cs2old,x1,x2
      real, dimension (nx) :: beq2_pencil, prof, tmpx
      real, dimension (nx,ny) :: ax, ay
      real, dimension(3) :: b_ext
      real, dimension (:,:,:,:), allocatable :: ap
      real, dimension (:,:), allocatable :: yz

      real :: beq2,RFPradB12,RFPradJ12
      real :: s,c,sph_har_der
      real :: cosalp, sinalp
      integer :: j,iyz,llp1
      logical :: lvectorpotential=.true.
!
      do j=1,ninit
!
        select case (initaa(j))
        case ('nothing'); if (lroot .and. j==1) print*,'init_aa: nothing'
        case ('zero', '0'); f(:,:,:,iax:iaz) = 0.0
        case ('rescale'); f(:,:,:,iax:iaz)=amplaa(j)*f(:,:,:,iax:iaz)
        case ('tanhxy'); call tanh_hyperbola(amplaa(j),f,iaa,sheet_position,sheet_thickness,sheet_hyp)
        case ('exponential'); call exponential(amplaa(j),f,iaa,kz_aa(j))
        case ('bsiny'); call acosy(amplaa(j),f,iaa,ky_aa(j))
        case ('mode'); call modev(amplaa(j),coefaa,f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('modeb'); call modeb(amplaa(j),coefbb,f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('sph_constb'); call sph_constb(amplaa(j),f,iaa)
        case ('const_lou'); call const_lou(amplaa(j),f,iaa)
        case ('power_randomphase')
          call power_randomphase(amplaa(j),initpower_aa,cutoff_aa,f,iax,iaz,lscale_tobox)
        case ('power_randomphase_hel')
          call power_randomphase_hel(amplaa(j),initpower_aa,initpower2_aa, &
            cutoff_aa,ncutoff_aa,kpeak_aa,f,iax,iaz,relhel_aa,kgaussian_aa, &
            lskip_projection_aa,lno_second_ampl_aa,lvectorpotential, &
            lscale_tobox, k1hel=k1hel, k2hel=k2hel)
        case ('random-isotropic-KS')
          call random_isotropic_KS(initpower_aa,f,iax,N_modes_aa)
        case ('random_isotropic_shell')
          call random_isotropic_shell(f,iax,amplaa(j),z1_aa,z2_aa)
        case ('gaussian-noise'); call gaunoise(amplaa(j),f,iax,iaz)
        case ('gaussian-noise-z'); call gaunoise(amplaa(j),f,iaz,iaz)
        case ('gaussian-noise-rprof')
          call gaunoise_rprof(amplaa(j),f,iax,iaz,rnoise_int,rnoise_ext)
        case ('gaussian-noise-zprof')
          tmp=amplaa(1)*0.5*(tanh((z-z1)/0.05)-tanh((z-z2)/0.05))
          call gaunoise(tmp,f,iax,iaz)
        case ('gaussian-noise-zprof2')
          tmp=amplaa(1)*0.5*(tanh((z-znoise_int)/0.05)-tanh((z-znoise_ext)/0.05))
          call gaunoise(tmp,f,iax,iaz)
!
!  ABC field (includes Beltrami fields when only one coefficient /= 0)
!
        case ('ABC_field')
          call ABC_field(f,iaa,kx_aa(j),ky_aa(j),kz_aa(j),ABCaa,x0aa,y0aa,z0aa,widthaa)
!
!  Beltrami fields, put k=-k to make sure B=curl(A) has the right phase
!
        case ('Beltrami-general')
               call beltrami_general(amplaa(j),f,iaa,kx_aa(j),ky_aa(j),kz_aa(j),phase_aa(j))
        case ('Beltrami-x')
               call beltrami(amplaa(j),f,iaa,KX=kx_aa(j),phase=phasex_aa(j),sigma=relhel_aa)
        case ('Beltrami-xy-samehel')
               call beltrami(amplaa(j),f,iaa,KX=kx_aa(j),phase=phasex_aa(j),sigma=relhel_aa)
               call beltrami(amplaa(j),f,iaa,KY=kx_aa(j),phase=phasex_aa(j),sigma=relhel_aa)
        case ('Beltrami-xy-diffhel')
               call beltrami(-amplaa(j),f,iaa,KX=kx_aa(j),phase=phasex_aa(j),sigma=relhel_aa)
               call beltrami(amplaa(j),f,iaa,KY=kx_aa(j),phase=phasex_aa(j),sigma=relhel_aa)
        case ('Beltrami-xy-mixed')
               call beltrami(-amplaa(j)*mix_factor,f,iaa,KX=kx_aa(j),phase=phasex_aa(j),sigma=relhel_aa)
               call beltrami(amplaa(j),f,iaa,KY=kx_aa(j),phase=phasex_aa(j),sigma=relhel_aa)
        case ('Beltrami-yy')
               call beltrami(amplaa(j),f,iaa,KX=kx_aa(j),phase=phasex_aa(j),sigma=relhel_aa)
               call beltrami(amplaa(j),f,iaa,KX=2*kx_aa(j),phase=phasex_aa(j),sigma=relhel_aa)
        case ('Beltrami-y')
               call beltrami(amplaa(j),f,iaa,KY=ky_aa(j),phase=phasey_aa(j),sigma=relhel_aa)
        case ('Beltrami-z')
               call beltrami(amplaa(j),f,iaa,KZ=kz_aa(j),phase=phasez_aa(j), &
                   sigma=relhel_aa,z0=z0aa,width=widthaa(1))
        case ('bihelical-z')
               call bihelical(amplaa(j),f,iaa,KZ=kz_aa(j),phase=phasez_aa(j))
        case ('bihelical-z-sym')
               call bihelical(amplaa(j),f,iaa,KZ=kz_aa(j),phase=phasez_aa(j),sym=.true.)
        case ('Beltrami-z-old'); call beltrami_old(amplaa(j),f,iaa,KZ=-kz_aa(j),phase=phasez_aa(j))
        case ('Beltrami-z-complex'); call beltrami_complex(amplaa(j),f,iaa,KZ=-kz_aa(j),phase=phasez_aa(j))
        case ('Beltrami-x-equ'); call beltrami(amplaa(j),f,iaa,KZ=-kz_aa(j),phase=phasez_aa(j),KX2=2*pi/Lxyz(1))
        case ('Beltrami-y-equ'); call beltrami(amplaa(j),f,iaa,KZ=-kz_aa(j),phase=phasez_aa(j),KY2=2*pi/Lxyz(2))
        case ('Beltrami-z-equ'); call beltrami(amplaa(j),f,iaa,KZ=-kz_aa(j),phase=phasez_aa(j),KZ2=2*pi/Lxyz(3))
!
        case ('Bessel-x'); call bessel_x(amplaa(j),f,iaa,kx_aa(j))
        case ('Bessel_Az-x'); call bessel_az_x(amplaa(j),f,iaa,kx_aa(j))
        case ('propto-ux'); call wave_uu(amplaa(j),f,iaa,kx=kx_aa(j))
        case ('propto-uy'); call wave_uu(amplaa(j),f,iaa,ky=ky_aa(j))
        case ('propto-uz'); call wave_uu(amplaa(j),f,iaa,kz=kz_aa(j))
        case ('diffrot'); call diffrot(amplaa(j),f,iay)
        case ('ver-tube'); call vtube(amplaa(j),f,iax,iaz,radius)
        case ('ver-tube-peri'); call vtube_peri(amplaa(j),f,iax,iaz,radius)
        case ('hor-tanh'); call htanh(amplaa(j),f,iaz,epsilonaa)
        case ('hor-tube'); call htube(amplaa(j),f,iax,iaz,radius,epsilonaa, &
                                     center1_x,center1_z)
        case ('hor-tube-x'); call htube_x(amplaa(j),f,iax,iaz,radius,epsilonaa, &
                                     center1_y,center1_z)
        case ('hor-tube_erf'); call htube_erf(amplaa(j),f,iax,iaz,radius,epsilonaa, &
                                     center1_x,center1_z,fluxtube_border_width)
        case ('hor-fluxlayer'); call hfluxlayer(amplaa(j),f,iaa,z0aa,widthaa(1))
        case ('hor-fluxlayer-y'); call hfluxlayer_y(amplaa(j),f,iaa,z0aa,widthaa(1))
        case ('hor-fluxlayer-y-theta'); call hfluxlayer_y_theta(amplaa(j),f,iaa)
        case ('ver-fluxlayer'); call vfluxlayer(amplaa(j),f,iaa,x0aa,widthaa(1))
        case ('mag-support'); call magsupport(amplaa(j),f,gravz,cs0,rho0)
        case ('arcade-x'); call arcade_x(amplaa(j),f,iaa,kx_aa(j),kz_aa(j))
        case ('halfcos-Bx'); call halfcos_x(amplaa(j),f,iaa)
        case ('halfcos-Bz'); call halfcos_z(amplaa(j),f,iaa)
        case ('uniform-Bx'); call uniform_x(amplaa(j),f,iaa)
        case ('uniform-By'); call uniform_y(amplaa(j),f,iaa)
        case ('uniform-Bz'); call uniform_z(amplaa(j),f,iaa)
        case ('uniform-Bphi'); call uniform_phi(amplaa(j),f,iaa)
        case ('phi_comp_over_r'); call phi_comp_over_r(amplaa(j),f,iaa)
        case ('phi_comp_over_r_noise')
          call phi_comp_over_r(amplaa(j),f,iaa)
          call gaunoise(1.0e-5*amplaa(j),f,iax,iaz)
        case ('Bz(x)', '3'); call vfield(amplaa(j),f,iaa)
        case ('vfield2'); call vfield2(amplaa(j),f,iaa)
        case ('bipolar'); call bipolar(amplaa(j),f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('bipolar_restzero'); call bipolar_restzero(amplaa(j),f,iaa,kx_aa(j),ky_aa(j))
        case ('vecpatternxy'); call vecpatternxy(amplaa(j),f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('xjump'); call bjump(f,iaa,by_left,by_right,bz_left,bz_right,widthaa(1),'x')
        case ('x-point_xy'); call xpoint(amplaa(j),f,iaz,center1_x,center1_y)
        case ('x-point_xy2'); call xpoint2(amplaa(j),f,iaz,center1_x,center1_y)
        case ('sinxsinz'); call sinxsinz(amplaa(j),f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('bhyperz'); call bhyperz(amplaa(j),f,iaa,kz_aa(j),non_ffree_factor)
        case ('sinxsinz_Hz'); call sinxsinz(amplaa(j),f,iaa,kx_aa(j),ky_aa(j),kz_aa(j),KKz=kz_aa(j))
        case ('sin2xsin2y'); call sin2x_sin2y_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case ('cosxcosy'); call cosx_cosy_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
!PJK
        case ('Bzcosxcosy'); call cosx_cosy_cosz(amplaa(j),f,iay,kx_aa(j),ky_aa(j),0.)
        case ('sinxsiny'); call sinx_siny_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case ('xsiny'); call x_siny_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case ('x1siny'); call x1_siny_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.,phasey_aa(j))
        case ('sinxcosz'); call sinx_siny_cosz(amplaa(j),f,iay,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('sinycosz'); call cosx_siny_cosz(amplaa(j),f,iax,kx_aa(j),ky_aa(j),0.)
        case ('cosysinz'); call cosy_sinz(amplaa(j),f,iax,ky_aa(j),kz_aa(j))
        case ('x3sinycosy'); call x3_siny_cosz(amplaa(j),f,iay,xyz0(1),xyz1(1),ky_aa(j),kz_aa(j))
        case ('x3cosycosz'); call x3_cosy_cosz(amplaa(j),f,iax,ky_aa(j),kz_aa(j))
        case ('Ax=cosysinz'); call cosy_sinz(amplaa(j),f,iax,ky_aa(j),kz_aa(j))
        case ('magnetogram'); call mag_init(f)
        case ('Bz_Az_file'); call mag_Az_init(f)
        case ('Axyz_file'); call file_init(f)
        case ('Bz-floor'); call mdi_init(f,.true.,z0aa)
        case ('magnetogram_nonperiodic'); call mdi_init(f,.false.,z0aa)
        case ('cosxcoscosy'); call cosx_coscosy_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case ('Bphi_cosy')
          do n=n1,n2; do m=m1,m2
             f(l1:l2,m,n,iax)=amplaa(j)*(cos(2*pi*ky_aa(j)*(y(m)-y0)/Lxyz(2)))/Lxyz(2)
          enddo; enddo
        case ('Az_cosh2x1')
          do n=n1,n2; do m=m1,m2
             f(l1:l2,m,n,iax)=0. 
             f(l1:l2,m,n,iay)=0.
             f(l1:l2,m,n,iaz)=amplaa(j)/(cosh(x(l1:l2))*cosh(x(l1:l2))) 
          enddo; enddo
        case ('By_sinx')
          do n=n1,n2; do m=m1,m2
             f(l1:l2,m,n,iax)=0. 
             f(l1:l2,m,n,iay)=0.
             f(l1:l2,m,n,iaz)=amplaa(j)*( cos(x(l1:l2))  )
          enddo; enddo
        case ('By_step')
          do n=n1,n2; do m=m1,m2
             f(l1:l2,m,n,iax)=0. 
             f(l1:l2,m,n,iay)=0.
             f(l1:l2,m,n,iaz)=2*amplaa(j)*step(x(l1:l2),xyz0(1)+Lxyz(1)/2.,widthaa(1)) - amplaa(j)
          enddo; enddo
        case ('crazy', '5'); call crazy(amplaa(j),f,iaa)
        case ('strange'); call strange(amplaa(j),f,iaa)
        case ('read_arr_file'); call read_outside_vec_array(f, "aa.arr", iaa)
        case ('read_bin_file'); call read_outside_vec_array(f, "ap.dat", iaa,.true.,amplaa(j))
        case ('sinwave-phase')
          call sinwave_phase(f,iax,ampl_ax(j),kx_ax(j),ky_ax(j),kz_ax(j),phase_ax(j))
          call sinwave_phase(f,iay,ampl_ay(j),kx_ay(j),ky_ay(j),kz_ay(j),phase_ay(j))
          call sinwave_phase(f,iaz,ampl_az(j),kx_az(j),ky_az(j),kz_az(j),phase_az(j))
        case ('coswave-phase')
          call coswave_phase(f,iax,ampl_ax(j),kx_ax(j),ky_ax(j),kz_ax(j),phase_ax(j))
          call coswave_phase(f,iay,ampl_ay(j),kx_ay(j),ky_ay(j),kz_ay(j),phase_ay(j))
          call coswave_phase(f,iaz,ampl_az(j),kx_az(j),ky_az(j),kz_az(j),phase_az(j))
        case ('sinwave-x'); call sinwave(amplaa(j),f,iaa,kx=kx_aa(j))
        case ('coswave-Ax-kx'); call coswave(amplaa(j),f,iax,kx=kx_aa(j))
        case ('coswave-Ax-ky'); call coswave(amplaa(j),f,iax,ky=ky_aa(j))
        case ('coswave-Ax-kz'); call coswave(amplaa(j),f,iax,kz=kz_aa(j))
        case ('coswave-Ay-kx'); call coswave(amplaa(j),f,iay,kx=kx_aa(j))
        case ('coswave-Ay-ky'); call coswave(amplaa(j),f,iay,ky=ky_aa(j))
        case ('coswave-Ay-kz'); call coswave(amplaa(j),f,iay,kz=kz_aa(j))
        case ('coswave-Az-kx'); call coswave(amplaa(j),f,iaz,kx=kx_aa(j))
        case ('coswave-Az-ky'); call coswave(amplaa(j),f,iaz,ky=ky_aa(j))
        case ('coswave-Az-kz'); call coswave(amplaa(j),f,iaz,kz=kz_aa(j))
        case ('sinwave-Ay-kz'); call sinwave_phase(f,iay,amplaa(j),kx_aa(j),ky_aa(j),kz_aa(j),phasez_aa(j))
        case ('dipole'); call dipole(f,iax,amplaa(j))
        case ('dipole_tor'); call dipole_tor(f,iax,amplaa(j))    !,ladd=.true.)
        case ('linear-zx')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iay)=-0.5*amplaa(j)*z(n)**2/Lxyz(3)
          enddo; enddo
        case ('spot_xy')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iaz)=amplaa(j)*exp(-((x(l1:l2)-xyz1(1))**2+(y(m))**2)/radius**2)
          enddo; enddo
        case ('spot_xz')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iay)=amplaa(j)*exp(-(x(l1:l2)**2+(z(n)-xyz1(3))**2)/radius**2)
          enddo; enddo
        case ('spot_spherical')
          do n=n1,n2; do m=m1,m2
             f(l1:l2,m,n,iaz)=amplaa(j)*exp(-((x(l1:l2)-xyz1(1))**2)/0.04**2)* &
                  exp(-(y(m)-th_spot*pi/180.)**2/radius**2)
          enddo; enddo
        case ('Az=x2')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iaz)=.25*pi_1*amplaa(j)*(x(l1:l2)/Lxyz(1))**2
          enddo; enddo
        case ('Az=x4')
!
!  Initial  Az=(r/2)^2 [1-(r/2)^2]  corresponds to Bphi=r/2 and Jz=1.
!
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iaz)=.25*amplaa(j)*x(l1:l2)**2*(1.-.25*x(l1:l2)**2)
          enddo; enddo
        case ('JzBz_cyl_ct')
!
!  Initial  Az=(r/2)^2 [1-(r/2)^2]  corresponds to Bphi=r/2 and Jz=1.
!
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iaz)=-amplaaJ(j)*(x(l1:l2)**2)/4
            f(l1:l2,m,n,iay)=amplaaB(j)*(x(l1:l2))/2
          enddo; enddo
        case ('JzBz_cyl_ct_sq')
!
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iaz)=-amplaaJ(j)*(x(l1:l2)**2)/4
            f(l1:l2,m,n,iay)=amplaaB(j)*x(l1:l2)/2*(RFPrad(j)**2-x(l1:l2)**2/2)
          enddo; enddo
!
!  (work in progress, koen 14/03/11)
!
        case ('JzBz_cyl_4_4')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iaz)=-amplaaJ(j)/4*x(l1:l2)**2*(RFPrad(j)**4-x(l1:l2)**4/9)
            f(l1:l2,m,n,iay)=amplaaB(j)/2*x(l1:l2)*(RFPrad(j)**4-x(l1:l2)**4/3)
          enddo; enddo
!
!  Uniform B-field in cylindrical coordinates
!
        case ('Bz=const(cyl)')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iay)=.5*x(l1:l2)*amplaaB(j)
          enddo; enddo
!
!  Logarithmic B-spiral in cylindrical coordinates
!
        case ('Bspiral(cyl)')
          cosalp=cos(alp_aniso*dtor)
          sinalp=sin(alp_aniso*dtor)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iax)=+amplaa(j)*sinalp*x(l1:l2)*cos(kz_aa(j)*z(n))
            f(l1:l2,m,n,iay)=-amplaa(j)*cosalp*x(l1:l2)*cos(kz_aa(j)*z(n))
            f(l1:l2,m,n,iaz)=0.
          enddo; enddo
!
!  generalized
!
        case ('JzBz_cyl_RFPradJB')
          RFPradB12=1./RFPradB
          RFPradJ12=1./RFPradJ
          x1=x(l1:l2)
          x2=x(l1:l2)**2
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iay)=+.50*amplaaB(j)*x1*(1.-.50*x2*RFPradB12)
            f(l1:l2,m,n,iaz)=-.25*amplaaJ(j)*x2*(1.-.25*x2*RFPradJ12)
          enddo; enddo
!
        case ('Alfven-x'); call alfven_x(amplaa(j),f,iuu,iaa,ilnrho,kx_aa(j))
        case ('Alfven-y'); call alfven_y(amplaa(j),f,iuu,iaa,ky_aa(j),mu0)
        case ('Alfven-z'); call alfven_z(amplaa(j),f,iuu,iaa,kz_aa(j),mu0)
        case ('Alfven-xy'); call alfven_xy(amplaa(j),f,iuu,iaa,kx_aa(j),ky_aa(j))
        case ('Alfven-xz'); call alfven_xz(amplaa(j),f,iuu,iaa,kx_aa(j),kz_aa(j))
        case ('Alfvenz-rot'); call alfvenz_rot(amplaa(j),f,iuu,iaa,kz_aa(j),Omega)
        case ('Alfvenz-bell'); call alfvenz_bell(amplaa(j),f,iuu,iaa,kz_aa(j),B_ext(3),J_ext(3))
        case ('Alfvenz-rot-shear'); call alfvenz_rot_shear(amplaa(j),f,iuu,iaa,kz_aa(j),Omega)
        case ('piecewise-dipole'); call piecew_dipole_aa (amplaa(j),inclaa,f,iaa)
        case ('Ferriere-uniform-Bx'); call ferriere_uniform_x(amplaa(j),f,iaa)
        case ('Ferriere-uniform-By'); call ferriere_uniform_y(amplaa(j),f,iaa)
        case ('robertsflow'); call robertsflow(amplaa(j),f,iaa,relhel_aa)
        case ('tony-nohel')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iay)=amplaa(j)/kz_aa(j)*cos(kz_aa(j)*2.*pi/Lz*z(n))
          enddo;enddo
        case ('tony-nohel-yz')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iay)=amplaa(j)/kx_aa(j)*sin(kx_aa(j)*2.*pi/Lx*x(l1:l2))
          enddo;enddo
        case ('tony-hel-xy')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iax)=amplaa(j)/kz_aa(j)*sin(kz_aa(j)*2.*pi/Lz*z(n))
            f(l1:l2,m,n,iay)=amplaa(j)/kz_aa(j)*cos(kz_aa(j)*2.*pi/Lz*z(n))
          enddo;enddo
        case ('tony-hel-yz')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iay)=amplaa(j)/kx_aa(j)*sin(kx_aa(j)*2.*pi/Lx*x(l1:l2))
            f(l1:l2,m,n,iaz)=amplaa(j)/kx_aa(j)*cos(kx_aa(j)*2.*pi/Lx*x(l1:l2))
          enddo;enddo
        case ('force-free-jet')
          lB_ext_pot=.true.
          call force_free_jet(mu_ext_pot)
        case ('Alfven-circ-x')
!
!  Circularly polarised Alfven wave in x direction.
!
          if (lroot) print*,'init_aa: circular Alfven wave -> x'
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iay) = amplaa(j)/kx_aa(j)*sin(kx_aa(j)*x(l1:l2))
            f(l1:l2,m,n,iaz) = amplaa(j)/kx_aa(j)*cos(kx_aa(j)*x(l1:l2))
          enddo;enddo
        case ('geo-benchmark-case1','geo-benchmark-case2'); call geo_benchmark_B(f)
!
        case ('torus-test'); call torus_test(amplaa(j),f)
!
! test case horizontal dipole for spherical shell polar boundary conditions
!
        case ('horizontal_dipole')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iax) =   0.
            f(l1:l2,m,n,iay) = - amplaa(j)*x(l1:l2)*sin(z(n))
            f(l1:l2,m,n,iaz) = - amplaa(j)*x(l1:l2)*cos(y(m))*cos(z(n))
          enddo; enddo
!
! test case vertical dipole for spherical shell polar boundary conditions
!
        case ('vertical_dipole')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iax) = 0.
            f(l1:l2,m,n,iay) = 0.
            f(l1:l2,m,n,iaz) = amplaa(j)*x(l1:l2)*sin(y(m))
          enddo; enddo
!
        case ('relprof')
          f(l1:l2,m1:m2,n1:n2,iax:iay)=A_relprof
!
        case ('inclined-dipole')
!
!  Inclined dipole initial condition. In principle can use precession as well (though for that it should be moved to
!  runtime). Works only for spherical coordinates, and needs global external storing of fields.
!
          if (.not.(lbx_ext_global.and.lby_ext_global.and.lbz_ext_global)) &
               call fatal_error("init_aa",&
               "inclined-dipole: switch lb[xyz]_ext_global=T in magnetic_start_pars")
          if (.not.lspherical_coords) &
               call fatal_error("init_aa",&
               "inclined-dipole: so far only implemented for spherical coordinates")
!
          c=cos(inclaa*pi/180); s=sin(inclaa*pi/180)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iglobal_bx_ext) = dipole_moment * 2*(c*costh(m) + s*sinth(m)*cos(z(n)))/x(l1:l2)**3
            f(l1:l2,m,n,iglobal_by_ext) = dipole_moment *   (c*sinth(m) - s*costh(m)*cos(z(n)))/x(l1:l2)**3
            f(l1:l2,m,n,iglobal_bz_ext) = dipole_moment *   (s*sin(z(n)))                      /x(l1:l2)**3
          enddo;enddo
!
        case ('Axy_from_file') !(prelim version to set Ax,Ay on one proc)
          open(1,file='Axy.dat',form='unformatted')
          read(1) ax,ay
          close(1)
          f(l1:l2,m1:m2,n1,iax)=ax*amplaa(1)
          f(l1:l2,m1:m2,n1,iay)=ay*amplaa(1)
!
        case ('A_from_file') !(prelim version to set Ax,Ay on one proc)
          allocate(ap(nx,ny,nz,3))
          open(1,file='aa.dat',form='unformatted')
          read(1) ap
          close(1)
          f(l1:l2,m1:m2,n1:n2,iax)=ap(:,:,:,1)*amplaa(1)
          f(l1:l2,m1:m2,n1:n2,iay)=ap(:,:,:,2)*amplaa(1)
          f(l1:l2,m1:m2,n1:n2,iaz)=ap(:,:,:,3)*amplaa(1)
          if (allocated(ap)) deallocate(ap)
!
        case ('B_ext_from_file')
          allocate(ap(mx,my,mz,6))
          call input_snap('ap.dat',ap,6,0)
          call input_snap_finalize
!
          if (iglobal_ax_ext/=0) f(:,:,:,iglobal_ax_ext) = ap(:,:,:,1)
          if (iglobal_ay_ext/=0) f(:,:,:,iglobal_ay_ext) = ap(:,:,:,2)
          if (iglobal_az_ext/=0) f(:,:,:,iglobal_az_ext) = ap(:,:,:,3)
!
          if (iglobal_bx_ext/=0) f(:,:,:,iglobal_bx_ext) = ap(:,:,:,4)
          if (iglobal_by_ext/=0) f(:,:,:,iglobal_by_ext) = ap(:,:,:,5)
          if (iglobal_bz_ext/=0) f(:,:,:,iglobal_bz_ext) = ap(:,:,:,6)
          call initiate_isendrcv_bdry(f)
          call finalize_isendrcv_bdry(f)
          if (allocated(ap)) deallocate(ap)
!
        case('spher-harm-poloidal')
          if (.not.lspherical_coords) call fatal_error("init_uu", &
              "spher-harm-poloidal only meaningful for spherical coordinates"//trim(initaa(j)))
          tmpx=(x(l1:l2)-xyz0(1))*(x(l1:l2)-xyz1(1))/x(l1:l2) + (xyz1(1) - 0.5*xyz0(1))         ! S/r
          prof=3.*(x(l1:l2)-xyz0(1))                                                            ! S' + S/r
          llp1=(ll_sh(j)+1)*ll_sh(j)
          if (lyang) then
            allocate(yz(2,ny*nz))
            call yin2yang_coors(costh(m1:m2),sinth(m1:m2),cosph(n1:n2),sinph(n1:n2),yz)
            iyz=1
            do m=m1,m2
              do n=n1,n2
!if (iproc_world==55) print*, 'm,n,yz=', m,n,yz(:,iyz)
                bb(:,1)=amplaa(j)*llp1*tmpx*ylm_other(yz(1,iyz),yz(2,iyz),ll_sh(j),mm_sh(j),sph_har_der)
                bb(:,2)=amplaa(j)*prof*sph_har_der
                if (mm_sh(j)/=0) then
                  bb(:,3) = -amplaa(j)*prof*mm_sh(j)/sin(yz(1,iyz))*sin(yz(2,iyz))/cos(yz(2,iyz))
                else
                  bb(:,3) = 0.
                endif
                call transform_thph_yy( bb, (/1,1,1/), f(l1:l2,m,n,iax:iaz), yz(1,iyz), yz(2,iyz) )
                iyz=iyz+1
              enddo
            enddo
          else
            do n=n1,n2
              do m=m1,m2
                f(l1:l2,m,n,iax) = amplaa(j)*llp1*tmpx*ylm(ll_sh(j),mm_sh(j),sph_har_der)
                f(l1:l2,m,n,iay) = amplaa(j)*prof*sph_har_der
                if (mm_sh(j)/=0) &      ! tb improved!
                  f(l1:l2,m,n,iaz) = -amplaa(j)*prof*mm_sh(j)/sinth(m)*sinph(n)/cosph(n)
              enddo
            enddo
          endif

        case('from-zaverage')
          call read_zaver(f,iax,iaz,source_zav,nzav,indzav)
        case default
!
!  Catch unknown values.
!
          call fatal_error('init_aa', &
              'init_aa value "' // trim(initaa(j)) // '" not recognised')
!
        endselect
!
!  End loop over initial conditions.
!
      enddo
!
!  Initialize individual modules, but need to do this only if
!  lmagn_mf is true.
!
      if (lmagn_mf) call init_aa_mf(f)
!
!  Interface for user's own initial condition.
!
      if (linitial_condition) call initial_condition_aa(f)
!
!  In 2-D with nz=1, setting Ax=Ay=0 makes sense, but shouldn't
!  be compulsory, so allow for this possibility in 2-D.
!
      if (lset_AxAy_zero) then
        if (nz==1) then
          f(:,:,:,iax:iay)=0.0
        else
          call fatal_error("init_aa","lset_AxAy_zero=T only allowed with nz=1")
        endif
      endif
!
!  Allow for pressure equilibrium (for isothermal tube)
!  assume that ghost zones have already been set.
!  corrected expression below for gamma /= 1 case.
!  The beq2 expression for 2*mu0*p is not general yet.
!
      if (lpress_equil.or.lpress_equil_via_ss) then
        if (lroot) print*,'init_aa: adjust lnrho to have pressure equilib; cs0=',cs0
        call boundconds_x(f)
        call initiate_isendrcv_bdry(f)
        call finalize_isendrcv_bdry(f)
        call boundconds_y(f)
        call boundconds_z(f)
!
        call get_shared_variable('cp', cp, caller='init_aa')
!
        do n=n1,n2
        do m=m1,m2
          call curl(f,iaa,bb)
          call dot2_mn(bb,b2)
          if (gamma==1.0) then
            f(l1:l2,m,n,ilnrho)=log(exp(f(l1:l2,m,n,ilnrho))-b2/(2.*cs0**2))
          else
            beq2=2.*rho0*cs0**2
            fact=max(1.0e-6,1.0-b2/beq2)
            if (lentropy.and.lpress_equil_via_ss) then
              if (lpress_equil_alt) then
                ssold=f(l1:l2,m,n,iss)
                cs2old=cs20*exp(gamma_m1*(f(l1:l2,m,n,ilnrho)-lnrho0) &
                  + gamma/cp*ssold)   ! generalised for cp /= 1
                cs2=cs2old-gamma*b2/(beq2*exp(f(l1:l2,m,n,ilnrho)-lnrho0))
                f(l1:l2,m,n,iss)=ssold+cp*gamma1*(log(cs2/cs20)-log(cs2old/cs20))  ! generalised for cp /= 1
              else
                f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+fact/gamma
              endif
            else
              if (lpress_equil_alt) then
                lnrho_old=f(l1:l2,m,n,ilnrho)
                cs2=cs20*exp(gamma_m1*(lnrho_old-lnrho0) &
                  +gamma/cp*f(l1:l2,m,n,iss))   ! generalised for cp /= 1
                f(l1:l2,m,n,ilnrho)=log(exp(lnrho_old)-b2*gamma/(2.*cs2))
                f(l1:l2,m,n,iss)=cp*gamma1*(log(cs2/cs20)- &  ! generalised for cp /= 1
                  gamma_m1*(f(l1:l2,m,n,ilnrho)-lnrho0))      ! lnrho0 added for generality
              else
                !f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)+fact/gamma_m1
                beq2_pencil=2.*rho0*cs0**2*exp(gamma*(f(l1:l2,m,n,ilnrho)-lnrho0))
                fact=max(1.0e-6,1.0-b2/beq2_pencil)
                f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)+alog(fact)/gamma
              endif
            endif
          endif
        enddo
        enddo
!
      endif
!
    endsubroutine init_aa
!***********************************************************************
    subroutine pencil_criteria_magnetic
!
!   All pencils that the Magnetic module depends on are specified here.
!
!  19-nov-04/anders: coded
!
      use Mpicomm, only: stop_it
!
      lpenc_requested(i_bb)=.true.
      if (.not.ladvective_gauge) lpenc_requested(i_uxb)=.true.
!
!  need uga always for advective gauge.
!
      if (ladvective_gauge) then
        if (lfargo_advection) then
          lpenc_requested(i_uuadvec_gaa)=.true.
        else
          lpenc_requested(i_uga)=.true.
        endif
      endif
!
      if (tauAD/=0.0) then
        lpenc_requested(i_jxb)=.true.
        lpenc_requested(i_jxbxb)=.true.
        lpenc_requested(i_b2)=.true.
      endif
!
!  anisotropic B-dependent diffusivity
!
      if (eta_aniso_BB/=0.0) then
        lpenc_requested(i_jb)=.true.
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_b2)=.true.
      endif
!
!  magneto friction model
!
      if (lmagneto_friction) then
        lpenc_requested(i_jxbxb)=.true.
        lpenc_requested(i_b2)=.true.
        lpenc_requested(i_jxb)=.true.
      endif
!
      if (lwrite_slices) then
        if (ivid_aps/=0) then
          lpenc_video(i_aps)=.true.
          lpenc_video(i_rcyl_mn)=.true.
        endif
        if (ivid_bb/=0) lpenc_video(i_bb)=.true.
        if (ivid_jj/=0) lpenc_video(i_jj)=.true.
        if (ivid_b2/=0) lpenc_video(i_b2)=.true.
        if (ivid_j2/=0) lpenc_video(i_j2)=.true.
        if (ivid_bb_sph/=0) lpenc_video(i_bb_sph)=.true.
        if (ivid_jb/=0) lpenc_video(i_jb)=.true.
        if (ivid_ab/=0) lpenc_video(i_ab)=.true.
        if (ivid_beta1/=0) lpenc_video(i_beta1)=.true.
        if (ivid_poynting/=0) then
          lpenc_video(i_uxb)=.true.
          lpenc_video(i_jxb)=.true.
        endif
      endif
!
      if ((lresi_eta_const.or.lresi_eta_tdep) .and. .not. lweyl_gauge &
          .and. .not. limplicit_resistivity) lpenc_requested(i_del2a) = .true.
!
      if (lresi_zdep) then
        if (.not. limplicit_resistivity) lpenc_requested(i_del2a) = .true.
        lpenc_requested(i_diva) = .true.
      endif
!
!  jj pencil always needed when in Weyl gauge
!
      if ((hall_term/=0.0.and.ldt).or.height_eta/=0.0.or.ip<=4.or. &
          lweyl_gauge.or.lspherical_coords.or.lJ_ext.or.ljj_as_aux.or. &
          lresi_eta_aniso) &
          lpenc_requested(i_jj)=.true.
      if (battery_term/=0.0) then
        lpenc_requested(i_fpres)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
      if ((.not.lweyl_gauge).and.(lresi_shell.or. &
          lresi_eta_shock.or.lresi_smagorinsky.or.lresi_smagorinsky_nusmag.or. &
          lresi_eta_shock2.or.lresi_xdep.or.lresi_ydep.or.lresi_xydep.or. &
          lresi_rdep.or.lresi_eta_shock_profz.or.lresi_eta_shock_profr.or. &
          lresi_smagorinsky_cross.or.lresi_spitzer.or. &
          lresi_eta_proptouz)) &
          lpenc_requested(i_del2a)=.true.
      if (lresi_rdep) then
         lpenc_requested(i_r_mn)=.true.
         lpenc_requested(i_r_mn1)=.true.
      endif
      if ((.not.lweyl_gauge).and.(lresi_eta_proptouz)) &
         lpenc_requested(i_diva)=.true.
         lpenc_requested(i_uij)=.true.
      if (lresi_eta_proptouz) then
         lpenc_requested(i_uu)=.true.
      endif
      if (lresi_sqrtrhoeta_const) then
        lpenc_requested(i_jj)=.true.
        lpenc_requested(i_rho1)=.true.
        if (.not.lweyl_gauge) then
          lpenc_requested(i_del2a)=.true.
          lpenc_requested(i_glnrho)=.true.
        endif
      endif
      if (lresi_eta_shock.or.lresi_eta_shock_profz.or. &
          lresi_eta_shock2.or.lresi_eta_shock_profr) then
        lpenc_requested(i_shock)=.true.
        if (.not.lweyl_gauge) then
          lpenc_requested(i_gshock)=.true.
          lpenc_requested(i_diva)=.true.
        endif
      endif
      if (lresi_etava) then
        lpenc_requested(i_etava)=.true.
        lpenc_requested(i_jj)=.true.
      endif
      if (lresi_etaj) then
        lpenc_requested(i_etaj)=.true.
        lpenc_requested(i_jj)=.true.
      endif
      if (lresi_etaj2) then
        lpenc_requested(i_etaj2)=.true.
        lpenc_requested(i_jj)=.true.
      endif
      if (lresi_etajrho) then
        lpenc_requested(i_etajrho)=.true.
        lpenc_requested(i_jj)=.true.
      endif
      if (lresi_shell) then
        lpenc_requested(i_r_mn)=.true.
        lpenc_requested(i_evr)=.true.
      endif
      if (lresi_eta_shock_profr) then
        lpenc_requested(i_r_mn)=.true.
      endif
      if (lresi_eta_shock_perp) then
        lpenc_requested(i_shock_perp)=.true.
        if (.not.lweyl_gauge) then
          lpenc_requested(i_gshock_perp)=.true.
          lpenc_requested(i_diva)=.true.
        endif
      endif
      if (lresi_spitzer.or.lresi_cspeed) then
        lpenc_requested(i_lnTT)=.true.
        if (lweyl_gauge) then
          lpenc_requested(i_jj)=.true.
        else
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2a)=.true.
          lpenc_requested(i_diva)=.true.
        endif
      endif
      if (lresi_vAspeed) then
        lpenc_requested(i_etava)=.true.
        lpenc_requested(i_va2)=.true.
        if (lweyl_gauge) then
          lpenc_requested(i_jj)=.true.
        else
          lpenc_requested(i_gva)=.true.
          lpenc_requested(i_del2a)=.true.
          lpenc_requested(i_diva)=.true.
        endif
      endif
!
!  Pencils requested for diamagnetism
!
      if (ldiamagnetism) then
        lpenc_requested(i_b2)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_jj)=.true.
      endif
!
      if (lupw_aa) then
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_aij)=.true.
      endif
      if (tau_relprof/=0 .or. ekman_friction_aa/=0) then
        lpenc_requested(i_aa)=.true.
      endif
      if (ladvective_gauge) then
        lpenc_requested(i_aa)=.true.
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_aij)=.true.
        lpenc_requested(i_uij)=.true.
      endif
      if (lresi_shell.or.lresi_xdep.or.lresi_ydep.or.lresi_xydep.or. &
          lresi_rdep.or.lresi_smagorinsky.or.lresi_smagorinsky_nusmag) then
           lpenc_requested(i_diva)=.true.
      endif
      if (lresi_smagorinsky_nusmag) then
         lpenc_requested(i_nu_smag)=.true.
      endif
      if (lresi_smagorinsky_cross) lpenc_requested(i_jo)=.true.
      if (lresi_hyper2) lpenc_requested(i_del4a)=.true.
      if (lresi_hyper3) lpenc_requested(i_del6a)=.true.
      if (lresi_hyper3_csmesh) lpenc_requested(i_cs2)=.true.
!
!  Note that for the cylindrical case, according to lpencil_check,
!  graddiva is not needed. We still need it for the lspherical_coords
!  case, although we should check this.
!  del2a now computed directly in all spherical so not reuired
!      if (lspherical_coords) lpenc_requested(i_graddiva)=.true.
      if (lentropy .or. lresi_smagorinsky .or. ltemperature) then
        lpenc_requested(i_j2)=.true.
      endif
      if (lresi_dust) lpenc_requested(i_rhop)=.true.
      if (lresi_anomalous) then
        lpenc_requested(i_jj)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
!
!  when b2 is needed for quenching factors
!
      if (J_ext_quench/=0) lpenc_requested(i_b2)=.true.
!
!  Pencils needed by thermodynamics.
!
      if (lentropy .or. ltemperature .or. lhydro) lpenc_requested(i_rho1)=.true.
      if (lentropy .or. ltemperature) lpenc_requested(i_TT1)=.true.
      if (lrhs_max.and.lentropy) lpenc_requested(i_cv1)=.true.
      if (ltemperature) lpenc_requested(i_cv1)=.true.
      if (lenergy .and. .not. lkinematic .and. lohmic_heat) then
        lpenc_requested(i_j2)=.true.
        if (lentropy .and. pretend_lnTT) lpenc_requested(i_cv1)=.true.
      endif
!
!  Ambipolar diffusion.
!
      if (lambipolar_diffusion) then
        lpenc_requested(i_nu_ni1)=.true.
        lpenc_requested(i_va2)=.true.
        lpenc_requested(i_jxbrxb)=.true.
        lpenc_requested(i_jxbr2)=.true.
        lpenc_requested(i_jxbr)=.true.
        if (ambipolar_diffusion=="ionization-equilibrium") &
          lpenc_requested(i_rho1)=.true.
        if (ambipolar_diffusion=="ionization-yH") &
          lpenc_requested(i_yH)=.true.
          lpenc_requested(i_rho1)=.true.
      endif
!
      if (hall_term/=0.0) lpenc_requested(i_jxb)=.true.
!
!  Take care of Lorentz force for potential flows.
!  In that case, use only the magnetic pressure.
!
      if (lhydro .and. llorentzforce) then
        if (iphiuu==0) then
          lpenc_requested(i_jxbr)=.true.
        else
          lpenc_requested(i_b2)=.true.
        endif
      endif
!
      if (lresi_smagorinsky_cross) lpenc_requested(i_oo)=.true.
!
      if (dipole_moment/=0.0) lpenc_requested(i_r_mn1)=.true.
!
! pencils for semirelativistic Boris correction to velocity eq.
!
      if (lboris_correction) then
        lpenc_requested(i_jxbr)=.true.
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_u2)=.true.
        lpenc_requested(i_cs2)=.true.
        lpenc_requested(i_gamma_A2)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_rho1gpp)=.true.
        lpenc_requested(i_ugu)=.true.
      endif
!
      if (lmagneto_friction.and.(.not.lhydro)) lpenc_requested(i_vmagfric)=.true.
!
!  ua pencil if lua_as_aux
!
      if (lua_as_aux) lpenc_diagnos(i_ua)=.true.
!
!  Request unit vectors for transformation of magnetic field from
!  Cartesian to spherical coordinates.
!
      if (lbb_sph_as_aux.and.lsphere_in_a_box) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_evr)=.true.
        lpenc_requested(i_evth)=.true.
        lpenc_requested(i_phix)=.true.
        lpenc_requested(i_phiy)=.true.
      endif
!
!  diagnostics pencils
!
      if (idiag_jxmax/=0 .or. idiag_jymax/=0 .or. idiag_jzmax/=0) lpenc_diagnos(i_jj)=.true.
      if (idiag_jxbxm/=0 .or. idiag_jybxm/=0 .or. idiag_jzbxm/=0 .or. idiag_jxbym/=0 .or. & 
          idiag_jxbzm/=0 .or. idiag_jybzm/=0 .or. idiag_jzbzm/=0) lpenc_diagnos(i_jj)=.true.
      if (idiag_jxbrxm/=0 .or. idiag_jxbrym/=0 .or. idiag_jxbrzm/=0) &
          lpenc_diagnos(i_jxbr)=.true.
      if (idiag_jxbrmax/=0) lpenc_diagnos(i_jxbr2)=.true.
      if (idiag_poynzmz/=0.or.idiag_jxbm/=0) lpenc_diagnos(i_jxb)=.true.
      if (idiag_jxbr2m/=0) lpenc_diagnos(i_jxbr2)=.true.
      if (idiag_jxbrxmx/=0 .or. idiag_jxbrymx/=0 .or. idiag_jxbrzmx/=0 .or. &
          idiag_jxbrxmy/=0 .or. idiag_jxbrymy/=0 .or. idiag_jxbrzmy/=0 .or. &
          idiag_jxbrxmz/=0 .or. idiag_jxbrymz/=0 .or. idiag_jxbrzmz/=0) &
          lpenc_diagnos(i_jxbr)=.true.
!
      if (idiag_b2ruzm/=0) &
          lpenc_diagnos(i_rho)=.true.
!
      if (idiag_hjbm/=0) lpenc_diagnos(i_hjb)=.true.
!
      if (     idiag_brmphi/=0  .or. idiag_uxbrmphi/=0 .or. idiag_jxbrmphi/=0 &
          .or. idiag_armphi/=0  .or. idiag_brmr/=0     .or. idiag_armr/=0 ) then
        lpenc_diagnos(i_pomx)=.true.
        lpenc_diagnos(i_pomy)=.true.
      endif
!
      if (     idiag_bpmphi/=0  .or. idiag_uxbpmphi/=0 .or. idiag_jxbpmphi/=0 &
          .or. idiag_bpmr/=0    .or. idiag_brbpmr/=0   .or. idiag_apmphi/=0  &
          .or. idiag_apmr/=0 ) then
        lpenc_diagnos(i_phix)=.true.
        lpenc_diagnos(i_phiy)=.true.
      endif
!
      if (idiag_armr/=0 .or. idiag_apmr/=0 .or. idiag_azmr/=0) &
          lpenc_diagnos(i_aa)=.true.
!
      if (idiag_armphi/=0 .or. idiag_apmphi/=0 .or. idiag_azmphi/=0 .or. &
          idiag_axmxz/=0 .or. idiag_aymxz/=0 .or. idiag_azmxz/=0 .or. &
          idiag_axmxy/=0 .or. idiag_aymxy/=0 .or. idiag_azmxy/=0) &
           lpenc_diagnos2d(i_aa)=.true.
!
      if (idiag_vAmxz/=0) lpenc_diagnos2d(i_va2)=.true.
!
      if (idiag_aybym2/=0 .or. idiag_exaym2/=0 .or. &
          idiag_examx/=0 .or. idiag_examy/=0 .or. idiag_examz/=0 .or. &
          idiag_examz1/=0 .or. idiag_examz2/=0 .or. idiag_examz3/=0 .or. &
          idiag_e3xamz1/=0 .or. idiag_e3xamz2/=0 .or. idiag_e3xamz3/=0 &
         ) lpenc_diagnos(i_aa)=.true.
!
      if (idiag_examxy1/=0 .or. idiag_examxy2/=0 .or. idiag_examxy3/=0 &
         ) lpenc_diagnos2d(i_aa)=.true.
!
      if (idiag_poynxmxy/=0 .or. idiag_poynymxy/=0 .or. idiag_poynzmxy/=0 &
         ) lpenc_diagnos2d(i_jxb)=.true.
      if (idiag_poynxmxy/=0 .or. idiag_poynymxy/=0 .or. idiag_poynzmxy/=0 .or. &
          idiag_Exmxy/=0 .or. idiag_Eymxy/=0 .or. idiag_Ezmxy/=0 &
         ) lpenc_diagnos2d(i_uxb)=.true.
!
      if (idiag_StokesImxy/=0) lpenc_diagnos2d(i_StokesI)=.true.
      if (idiag_StokesQmxy/=0) lpenc_diagnos2d(i_StokesQ)=.true.
      if (idiag_StokesUmxy/=0) lpenc_diagnos2d(i_StokesU)=.true.
      if (idiag_StokesQ1mxy/=0) lpenc_diagnos2d(i_StokesQ1)=.true.
      if (idiag_StokesU1mxy/=0) lpenc_diagnos2d(i_StokesU1)=.true.
      if (idiag_beta1mxy/=0) lpenc_diagnos2d(i_beta1)=.true.
!
      if (idiag_a2m/=0 .or. idiag_arms/=0 .or. idiag_amax/=0 &
      ) lpenc_diagnos(i_a2)=.true.
      if (idiag_divarms /= 0) lpenc_diagnos(i_diva) = .true.
      if (idiag_ab_int/=0 .or. idiag_abm/=0 .or. idiag_abmh/=0 &
          .or. idiag_abmz/=0 .or. idiag_abrms/=0 &
          .or. idiag_abumx/=0 .or. idiag_abumy/=0 .or. idiag_abumz/=0 &
          .or. idiag_abuxmz/=0 .or. idiag_abuymz/=0 .or. idiag_abuzmz/=0 &
         ) lpenc_diagnos(i_ab)=.true.
      if (idiag_abmxy/=0) lpenc_diagnos2d(i_ab)=.true.
      if (idiag_ubmxy/=0) lpenc_diagnos2d(i_ub)=.true.
!
      if (idiag_uam/=0 .or. idiag_uamz/=0) lpenc_diagnos(i_ua)=.true.
      if (idiag_djuidjbim/=0 &
          .or. idiag_dexbmx/=0 .or. idiag_dexbmy/=0 .or. idiag_dexbmz/=0 &
          .or. idiag_b3b12m/=0 .or. idiag_b3b21m/=0 &
          .or. idiag_b1b32m/=0 .or. idiag_b1b23m/=0 &
          .or. idiag_b2b13m/=0 .or. idiag_b2b31m/=0 ) &
          lpenc_diagnos(i_bij)=.true.
!
      if (idiag_divamz/=0 .or. idiag_bzdivamz/=0) lpenc_diagnos(i_diva)=.true.
      if (idiag_bzdivamz/=0) lpenc_diagnos(i_bb)=.true.
!
      if (idiag_j2m/=0 .or. idiag_jm2/=0 .or. idiag_jrms/=0 .or. &
          idiag_jmax/=0 .or. idiag_epsM/=0 .or. idiag_epsM_LES/=0 .or. &
          idiag_ajm/=0 .or. idiag_j2mz/=0 .or. idiag_epsMmz/=0) &
          lpenc_diagnos(i_j2)=.true.
!
      if (idiag_hjrms/=0 ) lpenc_diagnos(i_hj2)= .true.
      if (idiag_epsAD/=0 .and. &
          .not. (lambipolar_strong_coupling.and.tauAD/=0.0)) &
          lpenc_diagnos(i_jxbr2)=.true.
      if (idiag_jb_int/=0 .or. idiag_jbm/=0.or. idiag_jbmn/=0  .or. idiag_jbms/=0 &
          .or. idiag_jbmz/=0 .or. idiag_jbrms/=0 &
         ) lpenc_diagnos(i_jb)=.true.
      if (idiag_d6abmz/=0) lpenc_diagnos(i_d6ab)=.true.
      if (idiag_d6amz1/=0 .or. idiag_d6amz2 /=0 .or. idiag_d6amz3/=0) lpenc_diagnos(i_del6a)=.true.
      if (idiag_hjbm/=0 ) lpenc_diagnos(i_hjb)=.true.
      if (idiag_jbmphi/=0 .or. idiag_jbmxy/=0) lpenc_diagnos2d(i_jb)=.true.
      if (idiag_vArms/=0 .or. idiag_vAmax/=0 .or. idiag_vA2m/=0) lpenc_diagnos(i_va2)=.true.
      if (idiag_cosubm/=0) lpenc_diagnos(i_cosub)=.true.
      if (idiag_ubm/=0 .or. idiag_ubmz/=0 &
          .or. idiag_ubbzm/=0) lpenc_diagnos(i_ub)=.true.
!
      if (idiag_djuidjbim/=0 .or. idiag_uxDxuxbm/=0) lpenc_diagnos(i_uij)=.true.
      if (idiag_uxjm/=0) lpenc_diagnos(i_uxj)=.true.
      if (idiag_uxBrms/=0 .or. idiag_Rmrms/=0 .or. idiag_Rmmz/=0) &
          lpenc_diagnos(i_uxb2)=.true.
      if (idiag_beta1m/=0 .or. idiag_beta1max/=0 .or. idiag_beta1mz/=0) &
          lpenc_diagnos(i_beta1)=.true.
      if (idiag_betam /= 0 .or. idiag_betamax /= 0 .or. idiag_betamin /= 0 .or. &
          idiag_betamz /= 0 .or. idiag_beta2mz /= 0 .or. &
          idiag_betamx /= 0 .or. idiag_beta2mx /= 0) lpenc_diagnos(i_beta) = .true.
      if (idiag_bxmz/=0 .or. idiag_bymz/=0) lpenc_diagnos(i_bb)=.true.
      if (idiag_djuidjbim/=0) lpenc_diagnos(i_djuidjbi)=.true.
      if (idiag_b2divum/=0) lpenc_diagnos(i_divu)=.true.
      if (idiag_b2divum/=0) lpenc_diagnos(i_b2)=.true.
      if (idiag_ujxbm/=0) lpenc_diagnos(i_ujxb)=.true.
      if (idiag_gpxbm/=0) lpenc_diagnos(i_glnrhoxb)=.true.
      if (idiag_jxbxbm/=0) lpenc_diagnos(i_jxbxb)=.true.
      if (idiag_oxuxbm/=0) lpenc_diagnos(i_oxuxb)=.true.
      if (idiag_exaym2/=0 .or. idiag_exjm2/=0 &
          .or. idiag_jxmz/=0 .or. idiag_jymz/=0 &
          .or. idiag_jxpt/=0 .or. idiag_jypt/=0 .or. idiag_jzpt/=0 &
          .or. idiag_jxp2/=0 .or. idiag_jyp2/=0 .or. idiag_jzp2/=0 &
          .or. idiag_jmx/=0 .or. idiag_jmy/=0 .or. idiag_jmz/=0 &
          .or. idiag_ambmz/=0 .or. idiag_jmbmz/=0 .or. idiag_kmz/=0 &
          .or. idiag_examx/=0 .or. idiag_examy/=0 .or. idiag_examz/=0 &
          .or. idiag_examz1/=0 .or. idiag_examz2/=0 .or. idiag_examz3/=0 &
          .or. idiag_exjmx/=0 .or. idiag_exjmy/=0 .or. idiag_exjmz/=0 &
         ) lpenc_diagnos(i_jj)=.true.
!
       if (idiag_examxy1/=0 .or. idiag_examxy2/=0 .or. idiag_examxy3/=0 &
         ) lpenc_diagnos2d(i_jj)=.true.
!
      if (idiag_examz1/=0 .or. idiag_examz2/=0 .or. idiag_examz3/=0 .or. &
          idiag_exatop/=0 .or. idiag_exabot/=0 &
         ) lpenc_diagnos(i_exa)=.true.
      if (idiag_e3xamz1/=0 .or. idiag_e3xamz2/=0 .or. idiag_e3xamz3/=0 &
         ) lpenc_diagnos(i_e3xa)=.true.
!
      if (idiag_examxy1/=0 .or. idiag_examxy2/=0 .or. idiag_examxy3/=0 &
         ) lpenc_diagnos2d(i_exa)=.true.
!
      if (idiag_phibmx/=0 .or. idiag_phibmy/=0 .or. idiag_phibmz/=0 &
         ) lpenc_diagnos(i_diva)=.true.
      if (idiag_b2uzm/=0 .or. idiag_b2ruzm/=0 .or. &
          idiag_b1m/=0 .or. idiag_b2m/=0 .or. idiag_b4m/=0 .or. &
          idiag_bm2/=0 .or. idiag_EEM/=0 .or. &
          idiag_brmsh/=0 .or. idiag_brmsn/=0 .or. idiag_brmss/=0 .or. &
          idiag_brmsx/=0 .or. idiag_brmsz/=0 .or. &
          idiag_brms/=0 .or. idiag_bmax/=0 .or. &
          idiag_emag/=0 .or. idiag_b2mx /= 0 .or. idiag_b2mz/=0) &
          lpenc_diagnos(i_b2)=.true.
      if (idiag_bfrms/=0 .or.idiag_bf2m/=0 .or.  idiag_bf4m/=0 .or. &
          idiag_bf2mz/=0) lpenc_diagnos(i_bf2)=.true.
      if (idiag_etavamax/=0) lpenc_diagnos(i_etava)=.true.
      if (idiag_etajmax/=0) lpenc_diagnos(i_etaj)=.true.
      if (idiag_etaj2max/=0) lpenc_diagnos(i_etaj2)=.true.
      if (idiag_etajrhomax/=0) lpenc_diagnos(i_etajrho)=.true.
      if (idiag_cosjbm/=0) lpenc_diagnos(i_cosjb)=.true.
      if (idiag_coshjbm/=0) lpenc_diagnos(i_coshjb)=.true.
      if (idiag_cosubm/=0) lpenc_diagnos(i_cosub)=.true.
      if ((idiag_jparallelm/=0).or.(idiag_jperpm/=0)) then
        lpenc_diagnos(i_jparallel)=.true.
        lpenc_diagnos(i_jperp)=.true.
      endif
      if ((idiag_hjparallelm/=0).or.(idiag_hjperpm/=0)) then
        lpenc_diagnos(i_hjparallel)=.true.
        lpenc_diagnos(i_hjperp)=.true.
      endif
      if (idiag_b2mphi/=0 .or. idiag_b2mxz/=0) lpenc_diagnos2d(i_b2)=.true.
      if (idiag_brsphmphi/=0) lpenc_diagnos2d(i_evr)=.true.
      if (idiag_bthmphi/=0) lpenc_diagnos2d(i_evth)=.true.
      if (idiag_dbxdxmxy/=0.or.idiag_dbxdymxy/=0.or.idiag_dbxdzmxy/=0 .or. &
          idiag_dbydxmxy/=0.or.idiag_dbydymxy/=0.or.idiag_dbydzmxy/=0 .or. &
          idiag_dbzdxmxy/=0.or.idiag_dbzdymxy/=0.or.idiag_dbzdzmxy/=0)     &
        lpenc_diagnos2d(i_bijtilde)=.true.
      if (lisotropic_advection) lpenc_requested(i_va2)=.true.
      if (idiag_abumx/=0 .or. idiag_abumy/=0 .or. idiag_abumz/=0 &
          .or. idiag_abuxmz/=0 .or. idiag_abuymz/=0 .or. idiag_abuzmz/=0) &
          lpenc_diagnos(i_uu)=.true.
!
!  Check whether right variables are set for half-box calculations.
!
      if (idiag_brmsn/=0 .or. idiag_abmn/=0 .or. idiag_ambmzn/=0 &
          .or. idiag_jbmn/= 0 ) then
        if ((.not.lequatory).and.(.not.lequatorz)) then
          call stop_it("You have to set either of lequatory or lequatorz to true to calculate averages over half the box")
        else
          if (lequatory) write(*,*) 'pencil-criteria_magnetic: box divided along y dirn'
          if (lequatorz) write(*,*) 'pencil-criteria_magnetic: box divided along z dirn'
        endif
      endif
!
!  check for pencil_criteria_magn_mf
!
      if (lmagn_mf) call pencil_criteria_magn_mf
!
    endsubroutine pencil_criteria_magnetic
!***********************************************************************
    subroutine pencil_interdep_magnetic(lpencil_in)
!
!  Interdependency among pencils from the Magnetic module is specified here.
!
!  19-nov-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (.not.lcartesian_coords.and.(lpencil_in(i_bij).or.lpencil_in(i_del2a))) &
        lpencil_in(i_aij)=.true.

      if (lpencil_in(i_hjparallel).or.lpencil_in(i_hjperp)) then
        lpencil_in(i_coshjb)=.true.
        lpencil_in(i_hj2)=.true.
        lpencil_in(i_b2)=.true.
      endif
!
      if (lpencil_in(i_jparallel).or.lpencil_in(i_jperp)) then
        lpencil_in(i_cosjb)=.true.
        lpencil_in(i_jxb)=.true.
        lpencil_in(i_j2)=.true.
        lpencil_in(i_b2)=.true.
      endif
!
      if (lpencil_in(i_cosjb)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_j2)=.true.
        lpencil_in(i_jb)=.true.
      endif
      if (lpencil_in(i_coshjb)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_hj2)=.true.
        lpencil_in(i_hjb)=.true.
      endif
!
      if (lpencil_in(i_hjb)) then
        lpencil_in(i_bb)=.true.
        lpencil_in(i_hjj)=.true.
      endif
!
      if (lpencil_in(i_hj2)) lpencil_in(i_hjj)=.true.
!
      if (lpencil_in(i_hjj)) lpencil_in(i_del4a)=.true.
!
      if (lpencil_in(i_a2)) lpencil_in(i_aa)=.true.
!
      if (lpencil_in(i_ab)) then
        lpencil_in(i_aa)=.true.
        lpencil_in(i_bb)=.true.
      endif
!
      if (lpencil_in(i_ua)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_aa)=.true.
      endif
!
      if (lpencil_in(i_etava)) lpencil_in(i_va2)=.true.
      if (lpencil_in(i_jxbr) .and. va2max_jxb>0) lpencil_in(i_va2)=.true.
!
      if (lpencil_in(i_jxbr) .and. betamin_jxb>0) then
        lpencil_in(i_va2)=.true.
        lpencil_in(i_cs2)=.true.
      endif
!
      if (lpencil_in(i_va2)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_rho1)=.true.
      endif
!
      if (lpencil_in(i_etaj) .or. lpencil_in(i_etaj2) .or. lpencil_in(i_etajrho)) then
        lpencil_in(i_j2)=.true.
        lpencil_in(i_rho1)=.true.
      endif
!
      if (lpencil_in(i_j2)) lpencil_in(i_jj)=.true.
!
      if (lpencil_in(i_uxj)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_jj)=.true.
      endif
!
      if (lpencil_in(i_jb)) then
        lpencil_in(i_bb)=.true.
        lpencil_in(i_jj)=.true.
      endif
!
      if (lpencil_in(i_hjb)) then
        lpencil_in(i_bb)=.true.
        lpencil_in(i_hjj)=.true.
      endif
!
      if (lpencil_in(i_jxbr2)) lpencil_in(i_jxbr)=.true.
!
      if (lpencil_in(i_jxbr)) then
        lpencil_in(i_jxb)=.true.
        lpencil_in(i_rho1)=.true.
      endif
!
      if (lpencil_in(i_jxb)) then
        lpencil_in(i_jj)=.true.
        lpencil_in(i_bb)=.true.
      endif
!
      if (lpencil_in(i_uxb2)) lpencil_in(i_uxb)=.true.
!
      if (lpencil_in(i_uxb)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_bb)=.true.
      endif
!
      if (lpencil_in(i_cosub)) then
        lpencil_in(i_ub)=.true.
        lpencil_in(i_u2)=.true.
        lpencil_in(i_b2)=.true.
      endif
!
      if (lpencil_in(i_ub)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_bb)=.true.
      endif
!
      if (lpencil_in(i_chibp)) lpencil_in(i_bb)=.true.
      if (lpencil_in(i_StokesI)) lpencil_in(i_bb)=.true.
!
      if (lpencil_in(i_StokesQ).or.lpencil_in(i_StokesU).or.&
          lpencil_in(i_StokesQ1).or.lpencil_in(i_StokesU1)) then
        lpencil_in(i_chibp)=.true.
        lpencil_in(i_StokesI)=.true.
      endif
!
      if (lpencil_in(i_beta1) .or. lpencil_in(i_beta)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_pp)=.true.
      endif
!
      if (lpencil_in(i_bunit)) then
        lpencil_in(i_bb)=.true.
        lpencil_in(i_b2)=.true.
      endif
!
      if (lpencil_in(i_b2)) lpencil_in(i_bb)=.true.
      if ((lpencil_in(i_jj)) .and. .not. (ljj_as_comaux)) lpencil_in(i_bij)=.true.
!
      if (lpencil_in(i_djuidjbi)) then
        lpencil_in(i_uij)=.true.
        lpencil_in(i_bij)=.true.
      endif
!
      if (lpencil_in(i_jo)) then
        lpencil_in(i_oo)=.true.
        lpencil_in(i_jj)=.true.
      endif
!
      if (lpencil_in(i_ujxb)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_jxb)=.true.
      endif
!
      if (lpencil_in(i_oxu)) then
        lpencil_in(i_oo)=.true.
        lpencil_in(i_uu)=.true.
      endif
!
      if (lpencil_in(i_oxuxb)) then
        lpencil_in(i_oxu)=.true.
        lpencil_in(i_bb)=.true.
      endif
!
      if (lpencil_in(i_jxbxb)) then
        lpencil_in(i_jxb)=.true.
        lpencil_in(i_bb)=.true.
      endif
!
      if (lpencil_in(i_jxbrxb)) then
        lpencil_in(i_jxbr)=.true.
        lpencil_in(i_bb)=.true.
      endif
!
      if (lpencil_in(i_glnrhoxb)) then
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_bb)=.true.
      endif
!
      if (lpencil_in(i_oxj)) then
        lpencil_in(i_oo)=.true.
        lpencil_in(i_jj)=.true.
      endif
!
      if (lpencil_in(i_jij)) lpencil_in(i_bij)=.true.
!
      if (lpencil_in(i_sj)) then
        lpencil_in(i_sij)=.true.
        lpencil_in(i_jij)=.true.
      endif
      if (lpencil_in(i_ss12)) lpencil_in(i_sij)=.true.
      if (lpencil_in(i_del2a)) then
        if (lspherical_coords) then
!WL: for the cylindrical case, lpencil_check says these pencils are not needed
!AB: don't understand; this comment is in lspherical_coords.
!AB: I suggest you just fix it when you get to this point again.
          lpencil_in(i_jj)=.true.
          lpencil_in(i_graddivA)=.true.
        endif
      endif
      if (lpencil_in(i_uga)) then
        lpencil_in(i_aij)=.true.
        lpencil_in(i_uu)=.true.
      endif
!
      if (lpencil_in(i_uuadvec_gaa)) then
        lpencil_in(i_uu_advec)=.true.
        lpencil_in(i_aij)=.true.
        lpencil_in(i_uu)=.true.
        lpencil_in(i_aa)=.true.
      endif
!
      if (lpencil_in(i_d6ab) .or. lpencil_in(i_e3xa)) lpencil_in(i_del6a)=.true.
!
      if (lpencil_in(i_ss12)) lpencil_in(i_sj)=.true.
!
      if (lpencil_in(i_gva)) lpencil_in(i_va2)=.true.
!
! Interdependencies for Boris Correction
!
      if (lpencil_in(i_gamma_A2)) then
        lpencil_in(i_va2)=.true.
        lpencil_in(i_clight2)=.true.
      endif
!
! For magnetic friction
!
      if (lpencil_in(i_vmagfric)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_jxb)=.true.
      endif
!
! The dependence of diva or bb on aa relies on if aij is requested above.
!
      if ((lpencil_in(i_diva) .or. (lpencil_in(i_bb) .and. .not. lbb_as_comaux)) .and. &
          (.not. lcartesian_coords .and. lpencil_in(i_aij))) &
          lpencil_in(i_aa) = .true.
!
!  check for pencil_interdep_magn_mf
!
      if (lmagn_mf) call pencil_interdep_magn_mf(lpencil_in)
!
    endsubroutine pencil_interdep_magnetic
!***********************************************************************
    subroutine magnetic_before_boundary(f)
!
!  Conduct pre-processing required before boundary conditions and pencil
!  calculations.
!
!  30-may-14/ccyang: coded
!
      use Boundcond, only: update_ghosts, zero_ghosts
      use Sub, only: gij, gij_etc, curl_mn, dot2_mn
      use EquationOfState, only: rho0
!
      real, dimension(mx,my,mz,mfarray), intent(inout):: f
!
      real, dimension(nx,3,3) :: aij, bij
      real, dimension(nx,3) :: bb, jj
      real, dimension(nx) :: rho1, b2, tmp
      real, dimension(3) :: b_ext
      integer :: j
!
!  Find bb and jj if as communicated auxiliary.
!
      getbb: if (lbb_as_comaux .or. ljj_as_comaux .or. &
                 lalfven_as_aux.or. (lslope_limit_diff .and. llast)) then
        call zero_ghosts(f, iax, iaz)
        call update_ghosts(f, iax, iaz)
        mn_loop: do imn = 1, ny * nz
          m = mm(imn)
          n = nn(imn)
          call gij(f, iaa, aij, 1)
          call curl_mn(aij, bb, A=f(:,m,n,iax:iaz))
!
!  calculate jj if requested
!
          if(ljj_as_comaux) then
            if (lcartesian_coords) then
              call gij_etc(f,iaa,BIJ=bij)
              call curl_mn(bij,jj)
            else
              call gij_etc(f,iaa,AA=f(:,m,n,iax:iaz),AIJ=aij,BIJ=bij, &
                           LCOVARIANT_DERIVATIVE=lcovariant_magnetic)
              call curl_mn(bij,jj, A=bb,LCOVARIANT_DERIVATIVE=lcovariant_magnetic)
            endif
            f(l1:l2,m,n,ijx:ijz) = jj
          endif
!
!  Add imposed field, if any
!
          bext: if (lB_ext_in_comaux) then
            call get_bext(b_ext)
            forall(j = 1:3, b_ext(j) /= 0.0) bb(:,j) = bb(:,j) + b_ext(j)
            if (headtt .and. imn == 1) print *, 'magnetic_before_boundary: B_ext = ', b_ext
          endif bext
          if (lbb_as_comaux) f(l1:l2,m,n,ibx:ibz) = bb
!
!  Find alfven speed as communicated auxiliary
!
          if (lalfven_as_aux.or. (lslope_limit_diff .and. llast)) then
            if (ldensity) then
              if (ldensity_nolog) then
                rho1=1./f(l1:l2,m,n,irho)
              else
                rho1=1./exp(f(l1:l2,m,n,ilnrho))
              endif
            else
               rho1=1./rho0
            endif
            call dot2_mn(bb,b2)
            tmp= b2*mu01*rho1
            if (lcheck_positive_va2 .and. minval(tmp)<0.0) then
              print*, 'magnetic_before_boundary: Alfven speed is imaginary!'
              print*, 'magnetic_before_boundary: it, itsub, iproc=', it, itsub, iproc_world
              print*, 'magnetic_before_boundary: m, y(m), n, z(n)=', m, y(m), n, z(n)
              tmp=abs(tmp)
            endif
            if (lalfven_as_aux) f(l1:l2,m,n,ialfven)= tmp
            if (lslope_limit_diff .and. llast) then
              f(l1:l2,m,n,isld_char)=f(l1:l2,m,n,isld_char)+w_sldchar_mag*sqrt(tmp)
           endif
          endif
        enddo mn_loop
      endif getbb
!
    endsubroutine magnetic_before_boundary
!***********************************************************************
    subroutine update_char_vel_magnetic(f)
!
!   Add the Alfven speed to the characteritic velocity
!   for slope limited diffusion.
!
!   25-sep-15/MR+joern: coded
!
!      use General, only: staggered_mean_scal
      use General, only: staggered_max_scal
      use Density, only: calc_pencils_density

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      type(pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc=.false.
      real, dimension(nx) :: tmp
!
      if (lslope_limit_diff) then
!
!  Calculate Alfven speed and store temporarily in first slot of diffusive fluxes.
!
        lpenc_loc((/i_va2,i_bb,i_b2,i_rho1/))=.true.
        do n=n1,n2; do m=m1,m2

          call calc_pencils_density(f,p,lpenc_loc)
          call calc_pencils_magnetic(f,p,lpenc_loc)

          if (lboris_correction .and. va2max_boris>0) then
            tmp=(1+(p%va2/va2max_boris)**2.)**(-1.0/2.0)
            f(l1:l2,m,n,iFF_diff) = sqrt(p%va2*tmp)
          else
            f(l1:l2,m,n,iFF_diff) = sqrt(p%va2)
          endif

        enddo; enddo
!
        !call staggered_mean_scal(f,iFF_diff,iFF_char_c,w_sldchar_mag)
        call staggered_max_scal(f,iFF_diff,iFF_char_c,w_sldchar_mag)
      endif
!
    endsubroutine update_char_vel_magnetic
!***********************************************************************
    subroutine calc_pencils_magnetic_std(f,p)
!
!  Standard version (_std): global variable lpencil contains information about needed pencils.
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      type (pencil_case),                 intent(out)  :: p
!
      call calc_pencils_magnetic_pencpar(f,p,lpencil)
!
    endsubroutine calc_pencils_magnetic_std
!***********************************************************************
    subroutine calc_pencils_magnetic_pencpar(f,p,lpenc_loc)
!
!  Calculate Magnetic pencils.
!  Most basic pencils should come first, as others may depend on them.
! 
!  Version with formal parameter lpencil_loc instead of global lpencil for cases
!  in which not necessarily all generally needed pencil are to be calculated.
!
!  19-nov-04/anders: coded
!  18-jun-13/axel: b2 now includes B_ext by default (luse_Bext_in_b2=T is kept)
!  20-jun-16/fred: added derivative tensor option and streamlined gij_etc
!
      use Sub
      use EquationOfState, only: rho0
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      type (pencil_case),                 intent(out)  :: p
      logical, dimension(:),              intent(in)   :: lpenc_loc
!
      real, dimension (nx,3) :: tmp ! currently unused: bb_ext_pot
      real, dimension (nx) :: rho1_jxb, quench, StokesI_ncr, tmp1
      real, dimension(3) :: B_ext
      real :: c,s
      integer :: i,j,ix
! aa
      if (lpenc_loc(i_aa)) p%aa=f(l1:l2,m,n,iax:iaz)
! a2
      if (lpenc_loc(i_a2)) call dot2_mn(p%aa,p%a2)
! aij
      if (lpenc_loc(i_aij)) call gij(f,iaa,p%aij,1)
! diva
      if (lpenc_loc(i_diva)) then
!     ccyang: Note that the following two methods do not give exactly
!             the same results.
        if (lpenc_loc(i_aij) .and. .not. lpencil_check_at_work) then
          call div_mn(p%aij,p%diva,p%aa)
        else
          call div_other(f(:,:,:,iax:iaz),p%diva)
        endif
      endif
! aps
      if (lpenc_loc(i_aps)) p%aps=f(l1:l2,m,n,iaz)*p%rcyl_mn
! bb
      if (lpenc_loc(i_bb)) then
        if (lbb_as_comaux) then
          p%bb = f(l1:l2,m,n,ibx:ibz)
!     ccyang: Note that the following two methods do not give exactly
!             the same results.
        elseif (lpenc_loc(i_aij) .and. .not. lpencil_check_at_work) then
          call curl_mn(p%aij,p%bb,A=p%aa)
        else
          call curl_other(f(:,:,:,iax:iaz),p%bb)
        endif
!
!  Save field before adding imposed field (for diagnostics).
!  ##ccyang: Note that p%bb already contains B_ext if lbb_as_comaux 
!      .and. lB_ext_in_comaux = .true., which needs to be fixed.
!
        p%bbb = p%bb
!
!  Add a uniform background field, optionally precessing. 
!
        addBext: if (.not. (lbb_as_comaux .and. lB_ext_in_comaux) .and. &
                 (.not. ladd_global_field)) then
          call get_bext(B_ext)
          if (any(B_ext/=0.)) then
            forall(j = 1:3, B_ext(j) /= 0.0) p%bb(:,j) = p%bb(:,j) + B_ext(j)
            if (headtt) print *, 'calc_pencils_magnetic_pencpar: B_ext = ', B_ext
            if (headtt) print *, 'calc_pencils_magnetic_pencpar: logic = ', &
                        (lbb_as_comaux .and. lB_ext_in_comaux .and. ladd_global_field)
          endif
        endif addBext
!
!  Add a precessing dipole not in the Bext field
!
        if (dipole_moment .ne. 0) then
          c=cos(inclaa*pi/180); s=sin(inclaa*pi/180)
          p%bb(:,1) = p%bb(:,1) + dipole_moment * 2*(c*costh(m) + s*sinth(m)*cos(z(n)-omega_Bz_ext*t))*p%r_mn1**3
          p%bb(:,2) = p%bb(:,2) + dipole_moment *   (c*sinth(m) - s*costh(m)*cos(z(n)-omega_Bz_ext*t))*p%r_mn1**3
          p%bb(:,3) = p%bb(:,3) + dipole_moment *   (             s*         sin(z(n)-omega_Bz_ext*t))*p%r_mn1**3
          if (ladd_global_field) call fatal_error("calc_pencils_magnetic_pencpar",&
               "Switch ladd_global_field=F in magnetic_run_pars in run.in")
        endif
!
!  Add the external potential field.
!
!        if (lB_ext_pot) then
!          call get_global(bb_ext_pot,m,n,'B_ext_pot')
!          p%bb=p%bb+bb_ext_pot
!        endif
!
!  Add external B-field.
!
        if (ladd_global_field) then
          call get_bext(B_ext)
! Only need the second component to scale the global field and first and third to add a
! const oblique component at arbitrary deg inclination.
          if (iglobal_bx_ext/=0) p%bb(:,1)=p%bb(:,1)+B_ext(2)*f(l1:l2,m,n,iglobal_bx_ext)+B_ext(1)
          if (iglobal_by_ext/=0) p%bb(:,2)=p%bb(:,2)+B_ext(2)*f(l1:l2,m,n,iglobal_by_ext)
          if (iglobal_bz_ext/=0) p%bb(:,3)=p%bb(:,3)+B_ext(2)*f(l1:l2,m,n,iglobal_bz_ext)+B_ext(3)
        endif
      endif
!
!  b2 now (since 18 June 2013) includes B_ext by default.
!  This can be changed by setting lignore_Bext_in_b2=T
!
      if (lignore_Bext_in_b2 .or. (.not.luse_Bext_in_b2) ) then
        if (lpenc_loc(i_b2)) call dot2_mn(p%bbb,p%b2)
      else
        if (lpenc_loc(i_b2)) call dot2_mn(p%bb,p%b2)
      endif
      if (lpenc_loc(i_bf2)) call dot2_mn(p%bbb,p%bf2)
!
! rho=(rho0/10+B^2) !!!MR: below it is /100!
!
      if (lmagneto_friction.and.lpenc_loc(i_rho1)) then
        p%rho=rho0*1.0e-2+p%b2
        p%rho1=1./p%rho
      endif
! bunit
      if (lpenc_loc(i_bunit)) then
        quench = 1.0/max(tini,sqrt(p%b2))
        if (lignore_Bext_in_b2 .or. (.not.luse_Bext_in_b2) ) then
          do j=1,3
            p%bunit(:,j) = p%bbb(:,j)*quench
          enddo
        else
          do j=1,3
            p%bunit(:,j) = p%bb(:,j)*quench
          enddo
        endif
      endif
! ab
      if (lpenc_loc(i_ab)) call dot_mn(p%aa,p%bbb,p%ab)
      if (lpenc_loc(i_ua)) call dot_mn(p%uu,p%aa,p%ua)
! uxb
      if (lpenc_loc(i_uxb)) then
        call cross_mn(p%uu,p%bb,p%uxb)
        call cross_mn(p%uu,p%bbb,uxbb)
!  add external e-field.
        do j=1,3
          if (iglobal_eext(j)/=0) p%uxb(:,j)=p%uxb(:,j)+f(l1:l2,m,n,iglobal_eext(j))
        enddo
      endif
! uga
      if (lpenc_loc(i_uga)) call u_dot_grad(f,iaa,p%aij,p%uu,p%uga,UPWIND=lupw_aa)
!
! uga for fargo
!
      if (lpenc_loc(i_uuadvec_gaa)) then
        do j=1,3
          ! This is calling scalar h_dot_grad, that does not add
          ! the inertial terms. They will be added here.
          tmp = p%aij(:,j,:)
          call h_dot_grad(p%uu_advec,tmp,p%uuadvec_gaa(:,j))
        enddo
        if (lcylindrical_coords) then
          p%uuadvec_gaa(:,1) = p%uuadvec_gaa(:,1) - rcyl_mn1*p%uu(:,2)*p%aa(:,2)
          p%uuadvec_gaa(:,2) = p%uuadvec_gaa(:,2) + rcyl_mn1*p%uu(:,2)*p%aa(:,1)
        elseif (lspherical_coords) then
          call fatal_error("uuadvec_gaa","not implemented yet for spherical coordinates")
        endif
      endif
!
!  bij, del2a, graddiva
!  For non-cartesian coordinates jj is always required for del2a=graddiva-jj
!  fred: del2a now calculated directly if required and gradient tensor available
!  reduced calls to exclude unnecessary calculation of unwanted variables
!      if (lpenc_loc(i_bij) .or. lpenc_loc(i_del2a) .or. lpenc_loc(i_graddiva) .or. &
!          lpenc_loc(i_jj) ) then
!        if (lcartesian_coords) then
!          call gij_etc(f,iaa,p%aa,p%aij,p%bij,p%del2a,p%graddiva)
!          if (.not. lpenc_loc(i_bij)) p%bij=0.0      ! Avoid warnings from pencil
!          if (.not. lpenc_loc(i_del2a)) p%del2a=0.0  ! consistency check...
!          if (.not. lpenc_loc(i_graddiva)) p%graddiva=0.0
!!          if (lpenc_loc(i_jj)) call curl_mn(p%bij,p%jj,p%bb)
!!DM curl in cartesian does not need p%bb, then it is better not
!! to give it.
!          if (lpenc_loc(i_jj)) call curl_mn(p%bij,p%jj)
!        else
!          call gij_etc(f,iaa,AA=p%aa,AIJ=p%aij,BIJ=p%bij,&
!                               GRADDIV=p%graddiva,&
!                               LCOVARIANT_DERIVATIVE=lcovariant_magnetic)
!          if (.not. lpenc_loc(i_bij)) p%bij=0.0      ! Avoid warnings from pencil
!          call curl_mn(p%bij,p%jj,A=p%bb,LCOVARIANT_DERIVATIVE=lcovariant_magnetic)            ! consistency check...
!          if (lpenc_loc(i_del2a)) p%del2a=p%graddiva-p%jj
!!           if (lpenc_loc(i_del2a)) call del2v(f,iaa,p%del2a,p%aij,p%aa)
!        endif
!      endif
      if (lpenc_loc(i_bij).and.lpenc_loc(i_del2a)) then
        if (lcartesian_coords) then
          call gij_etc(f,iaa,BIJ=p%bij,DEL2=p%del2a)
          if (lpenc_loc(i_jj) .and. .not. ljj_as_comaux) call curl_mn(p%bij,p%jj)
        else
          call gij_etc(f,iaa,AA=p%aa,AIJ=p%aij,BIJ=p%bij,DEL2=p%del2a,&
!                               GRADDIV=p%graddiva,&
                               LCOVARIANT_DERIVATIVE=lcovariant_magnetic)
          if (lpenc_loc(i_jj) .and. .not. ljj_as_comaux) &
              call curl_mn(p%bij,p%jj,A=p%bb,LCOVARIANT_DERIVATIVE=lcovariant_magnetic)
        endif
      elseif (lpenc_loc(i_bij).and..not.lpenc_loc(i_del2a)) then
        if (lcartesian_coords) then
          call gij_etc(f,iaa,BIJ=p%bij)
          if (lpenc_loc(i_jj).and. .not. ljj_as_comaux) call curl_mn(p%bij,p%jj)
        else
          call gij_etc(f,iaa,AA=p%aa,AIJ=p%aij,BIJ=p%bij,&
                       LCOVARIANT_DERIVATIVE=lcovariant_magnetic)
          if (lpenc_loc(i_jj).and. .not. ljj_as_comaux) &
              call curl_mn(p%bij,p%jj,A=p%bb,LCOVARIANT_DERIVATIVE=lcovariant_magnetic)
        endif
      elseif (lpenc_loc(i_del2a).and..not.lpenc_loc(i_bij)) then
        if (lcartesian_coords) then
          call gij_etc(f,iaa,DEL2=p%del2a)
        else
          call gij_etc(f,iaa,AA=p%aa,AIJ=p%aij,DEL2=p%del2a,&
                       LCOVARIANT_DERIVATIVE=lcovariant_magnetic)
        endif
      endif
      if (lpenc_loc(i_bijtilde)) then
        if (lcovariant_magnetic) then
          call bij_tilde(f,p%bb,p%bijtilde,p%bij_cov_corr)
        else
          call bij_tilde(f,p%bb,p%bijtilde)
        endif
      endif
!
!     Possibility that jj is already calculated as comaux
!
      if (ljj_as_comaux .and. lpenc_loc(i_jj)) then
!
!  Apply a gaussian smoothing using 3 points before and after
!  make sense to use with SLD on magnetic field or vector potential
!  uses the kern_jjsmooth as kernels
!
        if (lsmooth_jj) then
          do j=1,3
            call smooth_mn(f,ijx+(j-1),kern_jjsmooth,p%jj(:,j))
          enddo
        else
          p%jj=f(l1:l2,m,n,ijx:ijz)
        endif
      endif
!
!  possibility of diamagnetism
!
      if (ldiamagnetism) call diamagnetism(p)
!
! jj
      if (lpenc_loc(i_jj)) then
! consistency check...
        p%jj=mu01*p%jj
!
!  Add external j-field.
!
        do j=1,3
          if (iglobal_jext(j)/=0) p%jj(:,j)=p%jj(:,j)+f(l1:l2,m,n,iglobal_jext(j))
        enddo
!
!  Add an external J-field (for the Bell instability).
!
        if (lJ_ext) then
          if (J_ext_quench/=0) then
            quench=1./(1.+J_ext_quench*p%b2)
            do j=1,3
              p%jj(:,j)=p%jj(:,j)-J_ext(j)*quench
            enddo
          else
            do j=1,3
              p%jj(:,j)=p%jj(:,j)-J_ext(j)
            enddo
          endif
        endif
      endif
! exa
      if (lpenc_loc(i_exa)) then
        call cross_mn(-p%uxb+eta*p%jj,p%aa,p%exa)
      endif
! j2
      if (lpenc_loc(i_j2)) call dot2_mn(p%jj,p%j2)
! jb
      if (lpenc_loc(i_jb)) call dot_mn(p%jj,p%bbb,p%jb)
! va2
      if (lpenc_loc(i_va2)) then
        p%va2=p%b2*mu01*p%rho1
        if (lcheck_positive_va2 .and. minval(p%va2)<0.0) then
          print*, 'calc_pencils_magnetic: Alfven speed is imaginary!'
          print*, 'calc_pencils_magnetic: it, itsub, iproc=', it, itsub, iproc_world
          print*, 'calc_pencils_magnetic: m, y(m), n, z(n)=', m, y(m), n, z(n)
          p%va2=abs(p%va2)
        endif
      endif
! eta_va
      if (lpenc_loc(i_etava)) then
        if (lresi_vAspeed) then
          p%etava = eta_va * sqrt(p%va2)/vArms
          if (va_min > 0.) where (p%etava < va_min) p%etava = va_min
        else
          p%etava = mu0 * eta_va * dxmax * sqrt(p%va2)
          if (eta_min > 0.) where (p%etava < eta_min) p%etava = 0.
        endif
      endif
! gradient of va
      if (lpenc_loc(i_gva).and.lalfven_as_aux) then
        call grad(f,ialfven,p%gva)
        if (lresi_vAspeed) then
          do i=1,3
            if (va_min > 0.) where (p%etava < va_min) p%gva(:,i) = 0.
          enddo
        endif
      endif
! eta_j
      if (lpenc_loc(i_etaj)) then
        p%etaj = mu0 * eta_j * dxmax**2 * sqrt(mu0 * p%j2 * p%rho1)
        if (eta_min > 0.) where (p%etaj < eta_min) p%etaj = 0.
      endif
! eta_j2
      if (lpenc_loc(i_etaj2)) then
        p%etaj2 = etaj20 * p%j2 * p%rho1
        if (eta_min > 0.) where (p%etaj2 < eta_min) p%etaj2 = 0.
      endif
! eta_jrho
      if (lpenc_loc(i_etajrho)) then
        p%etajrho = mu0 * eta_jrho * dxmax * sqrt(p%j2) * p%rho1
        if (eta_min > 0.) where (p%etajrho < eta_min) p%etajrho = 0.
      endif
! jxb
      if (lpenc_loc(i_jxb)) call cross_mn(p%jj,p%bb,p%jxb)
!
! cosjb
      if (lpenc_loc(i_cosjb)) then
        do ix=1,nx
          if ((abs(p%j2(ix))<=tini).or.(abs(p%b2(ix))<=tini))then
            p%cosjb(ix)=0.
          else
            p%cosjb(ix)=p%jb(ix)/sqrt(p%j2(ix)*p%b2(ix))
          endif
        enddo
        if (lpencil_check_at_work) then
        ! map penc0 value back to interval [-1,1]
          p%cosjb = modulo(p%cosjb + 1.0, 2.0) - 1
        endif
      endif
! jparallel and jperp
      if (lpenc_loc(i_jparallel).or.lpenc_loc(i_jperp)) then
        p%jparallel=sqrt(p%j2)*p%cosjb
        p%jperp=sqrt(p%j2)*sqrt(abs(1-p%cosjb**2))
      endif
! jxbr
      if (lpenc_loc(i_jxbr)) then
        rho1_jxb=p%rho1
!
!  Set rhomin_jxb>0 in order to limit the jxb term at very low densities.
!  Set va2max_jxb>0 in order to limit the jxb term at very high Alfven speeds.
!  Set va2power_jxb to an integer value in order to specify the power of the
!  limiting term,
!
        if (rhomin_jxb>0) rho1_jxb=min(rho1_jxb,1/rhomin_jxb)
        if (va2max_jxb>0 .and. (.not. betamin_jxb>0)) then
          rho1_jxb = rho1_jxb &
                   * (1+(p%va2/va2max_jxb)**va2power_jxb)**(-1.0/va2power_jxb)
        endif
        if (betamin_jxb>0) then
          va2max_beta = p%cs2/betamin_jxb*2.0*gamma1
          if (va2max_jxb > 0) va2max_beta=min(va2max_beta,va2max_jxb)
          rho1_jxb = rho1_jxb &
                   * (1+(p%va2/va2max_beta)**va2power_jxb)**(-1.0/va2power_jxb)
        endif
        call multsv_mn(rho1_jxb,p%jxb,p%jxbr)
      endif
! jxbr2
      if (lpenc_loc(i_jxbr2)) call dot2_mn(p%jxbr,p%jxbr2)
! ub
      if (lpenc_loc(i_ub)) call dot_mn(p%uu,p%bb,p%ub)
! cosub
      if (lpenc_loc(i_cosub)) then
        do ix=1,nx
          if ((abs(p%u2(ix))<=tini).or.(abs(p%b2(ix))<=tini)) then
            p%cosub(ix)=0.
          else
            p%cosub(ix)=p%ub(ix)/sqrt(p%u2(ix)*p%b2(ix))
          endif
        enddo
        if (lpencil_check) then
        ! map penc0 value back to interval [-1,1]
          p%cosub = modulo(p%cosub + 1.0, 2.0) - 1
        endif
      endif
! uxb2
      if (lpenc_loc(i_uxb2)) call dot2_mn(p%uxb,p%uxb2)
! uxj
      if (lpenc_loc(i_uxj)) call cross_mn(p%uu,p%jj,p%uxj)
! chibp
      if (lpenc_loc(i_chibp)) p%chibp=atan2(p%bb(:,2),p%bb(:,1))+.5*pi
! StokesI
      if (lpenc_loc(i_StokesI)) p%StokesI=(p%bb(:,1)**2+p%bb(:,2)**2)**exp_epspb
!
! StokesQ, StokesU, StokesQ1, and StokesU1
!
      if (lncr_correlated) then
        StokesI_ncr=p%StokesI*p%b2
        if (lpenc_loc(i_StokesQ)) p%StokesQ=-StokesI_ncr*cos(2.*p%chibp)
        if (lpenc_loc(i_StokesU)) p%StokesU=-StokesI_ncr*sin(2.*p%chibp)
        if (lpenc_loc(i_StokesQ1)) p%StokesQ1=+StokesI_ncr*sin(2.*p%chibp)*p%bb(:,3)
        if (lpenc_loc(i_StokesU1)) p%StokesU1=-StokesI_ncr*cos(2.*p%chibp)*p%bb(:,3)
      elseif (lncr_anticorrelated) then
        StokesI_ncr=p%StokesI/(1.+ncr_quench*p%b2)
        if (lpenc_loc(i_StokesQ)) p%StokesQ=-StokesI_ncr*cos(2.*p%chibp)
        if (lpenc_loc(i_StokesU)) p%StokesU=-StokesI_ncr*sin(2.*p%chibp)
        if (lpenc_loc(i_StokesQ1)) p%StokesQ1=+StokesI_ncr*sin(2.*p%chibp)*p%bb(:,3)
        if (lpenc_loc(i_StokesU1)) p%StokesU1=-StokesI_ncr*cos(2.*p%chibp)*p%bb(:,3)
      else
        if (lpenc_loc(i_StokesQ)) p%StokesQ=-p%StokesI*cos(2.*p%chibp)
        if (lpenc_loc(i_StokesU)) p%StokesU=-p%StokesI*sin(2.*p%chibp)
        if (lpenc_loc(i_StokesQ1)) p%StokesQ1=+p%StokesI*sin(2.*p%chibp)*p%bb(:,3)
        if (lpenc_loc(i_StokesU1)) p%StokesU1=-p%StokesI*cos(2.*p%chibp)*p%bb(:,3)
      endif
!
! beta1
      if (lpenc_loc(i_beta1)) p%beta1=0.5*p%b2*mu01/p%pp
! beta
      if (lpenc_loc(i_beta)) p%beta = 2.0 * mu0 * p%pp / max(p%b2, epsilon(1.0))
! djuidjbi
      if (lpenc_loc(i_djuidjbi)) call multmm_sc(p%uij,p%bij,p%djuidjbi)
! jo
      if (lpenc_loc(i_jo)) call dot(p%jj,p%oo,p%jo)
! ujxb
      if (lpenc_loc(i_ujxb)) call dot_mn(p%uu,p%jxb,p%ujxb)
! oxu
      if (lpenc_loc(i_oxu)) call cross_mn(p%oo,p%uu,p%oxu)
! oxuxb
      if (lpenc_loc(i_oxuxb)) call cross_mn(p%oxu,p%bb,p%oxuxb)
! jxbxb
      if (lpenc_loc(i_jxbxb)) call cross_mn(p%jxb,p%bb,p%jxbxb)
! jxbrxb
      if (lpenc_loc(i_jxbrxb)) call cross_mn(p%jxbr,p%bb,p%jxbrxb)
! glnrhoxb
      if (lpenc_loc(i_glnrhoxb)) call cross_mn(p%glnrho,p%bb,p%glnrhoxb)
! del4a
      if (lpenc_loc(i_del4a)) call del4v(f,iaa,p%del4a)
! hjj
      if (lpenc_loc(i_hjj)) p%hjj = p%del4a
! hj2
      if (lpenc_loc(i_hj2)) call dot2_mn(p%hjj,p%hj2)
! hjb
      if (lpenc_loc(i_hjb)) call dot_mn(p%hjj,p%bb,p%hjb)
! coshjb
      if (lpenc_loc(i_coshjb)) then
        do ix=1,nx
          if ((abs(p%hj2(ix))<=tini).or.(abs(p%b2(ix))<=tini))then
            p%coshjb(ix)=0.
          else
            !p%coshjb(ix)=p%hjb(ix)/sqrt(p%hj2(ix)*p%b2(ix))
            p%coshjb(ix)=p%hjb(ix)/(sqrt(p%hj2(ix))*sqrt(p%b2(ix)))
          endif
        enddo
        if (lpencil_check_at_work) then
! map penc0 value back to interval [-1,1]
          p%coshjb = modulo(p%coshjb + 1.0, 2.0) - 1
        endif
      endif
! hjparallel and hjperp
      if (lpenc_loc(i_hjparallel).or.lpenc_loc(i_hjperp)) then
        p%hjparallel=sqrt(p%hj2)*p%coshjb
        p%hjperp=sqrt(p%hj2)*sqrt(abs(1-p%coshjb**2))
      endif
! del6a
      if (lpenc_loc(i_del6a)) call del6v(f,iaa,p%del6a)
! e3xa
      if (lpenc_loc(i_e3xa)) then
        call cross_mn(-p%uxb+eta_hyper3*p%del6a,p%aa,p%e3xa)
      endif
! oxj
      if (lpenc_loc(i_oxj)) call cross_mn(p%oo,p%jj,p%oxJ)
! jij
      if (lpenc_loc(i_jij)) then
        do j=1,3
          do i=1,3
            p%jij(:,i,j)=.5*(p%bij(:,i,j)+p%bij(:,j,i))
          enddo
        enddo
      endif
! d6ab
      if (lpenc_loc(i_d6ab)) call dot_mn(p%del6a,p%bb,p%d6ab)
! sj
      if (lpenc_loc(i_sj)) call multmm_sc(p%sij,p%jij,p%sj)
! ss12
      if (lpenc_loc(i_ss12)) p%ss12=sqrt(abs(p%sj))
! vmagfric
      if (lpenc_loc(i_vmagfric).and.numag/=0.0) then
        tmp1=mu01/(numag*(B0_magfric/unit_magnetic**2+p%b2))
        do i=1,3
          p%vmagfric(:,i)=abs(p%jxb(:,i))*tmp1
        enddo
      endif
!
!  Store bb, jj or jxb in auxiliary variable if requested.
!  Just neccessary immediately before writing snapshots, but how would we
!  know we are?
!
      if (lbb_as_aux .and. .not. lbb_as_comaux) f(l1:l2,m,n,ibx:ibz) = p%bb
      if (ljj_as_aux .and. .not. ljj_as_comaux) f(l1:l2,m,n,ijx:ijz) = p%jj
      if (ljxb_as_aux) f(l1:l2,m,n,ijxbx:ijxbz)=p%jxb
      if (lpenc_loc(i_bb_sph).and.lbb_sph_as_aux) p%bb_sph(:,1:3)=f(l1:l2,m,n,ibb_sphr:ibb_sphp)
!
!  Calculate magnetic mean-field pencils.
!  This should always be done after calculating magnetic pencils.
!
      if (lmagn_mf) call calc_pencils_magn_mf(f,p)
!
!  Ambipolar diffusion pencil
!
      if (lpenc_loc(i_nu_ni1)) call set_ambipolar_diffusion(p)
!
! reduced speed of light pencil
!
     if (lpenc_loc(i_clight2)) then
       if (va2max_boris > 0) then
         p%clight2=spread(va2max_boris,1,nx)
       else
         if (lcartesian_coords) then
           clight2_zdep(n-n1+1) = max(dble(cmin)**2,c_light**2/(1.+max(z(n),0.0)**8)+max(25.0*maxval(p%u2),maxval(p%cs2)))
           p%clight2=clight2_zdep(n-n1+1)
         else if (lspherical_coords) then
           p%clight2=spread(max(cmin**2,25*maxval(p%u2),maxval(p%cs2)),1,nx)
         endif
       endif
       p%gamma_A2=p%clight2/(p%clight2+p%va2+tini)
     endif
!
    endsubroutine calc_pencils_magnetic_pencpar
!***********************************************************************
    subroutine set_ambipolar_diffusion(p)
!
!  Subroutine to choose between various models of ion density
!  to enter in the ambipolar diffusion calculation. The "scale-density"
!  model assumes equilibrium between ionization and recombination:
!
!    ksi * rho = beta * rho_i * rho_e
!              ~ beta * rho_i^2
!
!  wkere ksi is ionization rate, beta the recombination rate, and
!  rho, rho_i, and rho_e are the neutral, ion, and electron density,
!  respectively, with rho >> rho_i ~ rho_e. Solving for rho_i,
!
!    rho_i   =    sqrt(ksi*rho/beta)
!          propto sqrt(rho)
!
!  This is implemented as ionization-equilibrium below. The proportionality
!  coefficients are incorporated into nu_ni1.
!
!  08-mar-13/wlad: coded
!
      type (pencil_case) :: p
!
      select case (ambipolar_diffusion)
!
      case('constant'); p%nu_ni1=nu_ni1
      case('ionization-equilibrium'); p%nu_ni1=nu_ni1*sqrt(p%rho1)
      case('ionization-yH'); p%nu_ni1=nu_ni1*sqrt(p%rho1)*(1.-p%yH)/p%yH
      case default
         write(unit=errormsg,fmt=*) &
              'set_ambipolar_diffusion: No such value for ambipolar_diffusion: ', &
              trim(ambipolar_diffusion)
         call fatal_error('set_ambipolar_diffusion',errormsg)
!
      endselect
!
    endsubroutine set_ambipolar_diffusion
!***********************************************************************
    subroutine diamagnetism(p)
!
!  Compute diamagnetism
!
      use Sub
!
      real, dimension (nx) :: chi_diamag
      real, dimension (nx,3) :: gchi_diamag, Bk_Bki, jj_diamag, tmp
      type (pencil_case) :: p
!
      intent(inout)  :: p
!
!  cmpute chi, and gradchi
!
      chi_diamag=B2_diamag/p%b2
!
!  Add (1/2)*grad[qp*B^2]. This initializes p%jxb_mf.
!
      call multmv_transp(p%bij,p%bb,Bk_Bki) !=1/2 grad B^2
      call multsv(-.5*chi_diamag/p%b2,Bk_Bki,gchi_diamag)
      call cross(gchi_diamag,p%bb,jj_diamag)
      call multsv_add(jj_diamag,chi_diamag,p%jj,tmp)
      jj_diamag=tmp
!
!  update current density
!
      p%jj=p%jj+jj_diamag
!
    endsubroutine diamagnetism
!***********************************************************************
    subroutine daa_dt(f,df,p)
!
!  Magnetic field evolution.
!
!  Calculate dA/dt=uxB+3/2 Omega_0 A_y x_dir -eta mu_0 J.
!  Add jxb/rho to momentum equation.
!  Add eta mu_0 j2/rho to entropy equation.
!
!  22-nov-01/nils: coded
!   1-may-02/wolf: adapted for pencil_modular
!  17-jun-03/ulf:  added bx^2, by^2 and bz^2 as separate diagnostics
!   8-aug-03/axel: introduced B_ext21=1./B_ext**2, and set =1 to avoid div. by 0
!  12-aug-03/christer: added alpha effect (alpha in the equation above)
!  26-may-04/axel: ambipolar diffusion added
!  18-jun-04/axel: Hall term added
!   9-apr-12/MR: upwinding for ladvective_gauge=F generalized
!  31-mar-13/axel: Stokes parameter integration from synchrotron emission
!  25-aug-13/MR: simplified calls of save_name_sound
!  15-oct-15/MR: changes for slope-limited diffusion
!  14-apr-16/MR: changes for Yin-Yang: only yz slices at the moment!
!   4-aug-17/axel: implemented terms for ultrarelativistic EoS
!  03-apr-20/joern: restructured and fixed slope-limited diffusion
!
      use Debug_IO, only: output_pencil
      use Deriv, only: der6
      use Mpicomm, only: stop_it
      use Special, only: special_calc_magnetic
      use Sub
      use General, only: transform_thph_yy, notanumber

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)   :: p
      intent(inout):: f,df
!
      real, dimension (nx,3) :: ujiaj,gua,ajiuj
      real, dimension (nx,3) :: aa_xyaver
      real, dimension (nx,3) :: geta,uxb_upw,tmp2
      real, dimension (nx,3) :: dAdt, gradeta_shock
      real, dimension (nx,3,3) :: d_sld_flux
      real, dimension (nx) :: ftot, dAtot
      real, dimension (nx) :: peta_shock
      real, dimension (nx) :: sign_jo,rho1_jxb,tmp1
      real, dimension (nx) :: eta_mn,etaSS
      real, dimension (nx) :: vdrift
      real, dimension (nx) :: del2aa_ini,tanhx2,advec_hall,advec_hypermesh_aa
      real, dimension(nx) :: eta_BB, prof
      real, dimension(3) :: B_ext
      real :: tmp, eta_out1, maxetaBB=0., cosalp, sinalp, hall_term_
      real, parameter :: OmegaSS=1.0
      integer :: i,j,k,ju,ix
      integer, parameter :: nxy=nxgrid*nygrid
!
!  Identify module and boundary conditions.
!
      call timing('daa_dt','entered',mnloop=.true.)
      if (headtt.or.ldebug) print*,'daa_dt: SOLVE'
      if (headtt) then
        call identify_bcs('Ax',iax)
        call identify_bcs('Ay',iay)
        call identify_bcs('Az',iaz)
        if(lbb_as_comaux) then
          call identify_bcs('Bx',ibx)
          call identify_bcs('By',iby)
          call identify_bcs('Bz',ibz)
        endif
        if(ljj_as_comaux) then
          call identify_bcs('Jx',ijx)
          call identify_bcs('Jy',ijy)
          call identify_bcs('Jz',ijz)
        endif
        if(lslope_limit_diff) then
          call identify_bcs('sld_char',isld_char)
        endif
      endif
!
! set dAdt to zero at the beginning of each execution of this routine
!
      dAdt=0.
      Fmax=tini
      dAmax=tini
!
!  Replace B_ext locally to accommodate its time dependence.
!
      call get_bext(B_ext)
!
!  Add jxb/rho to momentum equation.
!
      if (lhydro) then
        if (.not.lkinematic) then
          if (llorentzforce) then
            if (lboris_correction) then
!
!  Following Eq. 34 of Gombosi et al. 2002 for Boris correction. Can work with
!  only const gravity at present
!
                df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+p%gamma_A2*p%jxbr(:,1)+&
                                  (p%ugu(:,1)+p%rho1gpp(:,1)-p%gg(:,1))*(1-p%gamma_A2)-&
                                  mu01*(p%gamma_A2**2*p%rho1/p%clight2)* &
                                  (p%bb(:,1)**2*(p%ugu(:,1)+p%rho1gpp(:,1)-p%gg(:,1))+&
                                  p%bb(:,1)*p%bb(:,2)*(p%ugu(:,2)+p%rho1gpp(:,2)-p%gg(:,2))+&
                                  p%bb(:,1)*p%bb(:,3)*(p%ugu(:,3)+p%rho1gpp(:,3)-p%gg(:,3)))
                df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+p%gamma_A2*p%jxbr(:,2)+&
                                  (p%ugu(:,2)+p%rho1gpp(:,2)-p%gg(:,2))*(1-p%gamma_A2)-&
                                  mu01*(p%gamma_A2**2*p%rho1/p%clight2)* &
                                  (p%bb(:,2)**2*(p%ugu(:,2)+p%rho1gpp(:,2)-p%gg(:,2))+&
                                  p%bb(:,2)*p%bb(:,1)*(p%ugu(:,1)+p%rho1gpp(:,1)-p%gg(:,1))+&
                                  p%bb(:,2)*p%bb(:,3)*(p%ugu(:,3)+p%rho1gpp(:,3)-p%gg(:,3)))
                df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+p%gamma_A2*p%jxbr(:,3)+&
                                  (p%ugu(:,3)+p%rho1gpp(:,3)-p%gg(:,3))*(1-p%gamma_A2)-&
                                  mu01*(p%gamma_A2**2*p%rho1/p%clight2)* &
                                  (p%bb(:,3)**2*(p%ugu(:,3)+p%rho1gpp(:,3)-p%gg(:,3))+&
                                  p%bb(:,3)*p%bb(:,1)*(p%ugu(:,1)+p%rho1gpp(:,1)-p%gg(:,1))+&
                                  p%bb(:,3)*p%bb(:,2)*(p%ugu(:,2)+p%rho1gpp(:,2)-p%gg(:,2)))
            elseif (llorentz_rhoref) then
              df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+p%jxb*rhoref1
            else
              if (ldensity.and.lrelativistic_eos) then
                call dot(p%uu,p%jxbr,tmp1)
                df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+tmp1
                df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+.75*p%jxbr
              else
                if (iphiuu==0) then
                  df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+p%jxbr
                else
                  df(l1:l2,m,n,iphiuu)=df(l1:l2,m,n,iphiuu)-.5*(p%b2-B_ext2)
                endif
              endif
            endif
          endif
        endif
      endif
!
!  Restivivity term
!
!  Because of gauge invariance, we can add the gradient of an arbitrary scalar
!  field Phi to the induction equation without changing the magnetic field,
!    dA/dt = u x B - eta j + grad(Phi).
!
!  If lweyl_gauge=T, we choose Phi = const. and solve
!    dA/dt = u x B - eta mu0 j.
!
!  Else, if lweyl_gauge=F, we make the gauge choice Phi = eta div(A)
!  and thus solve
!    dA/dt = u x B + eta laplace(A) + div(A) grad(eta).
!
!  Note: lweyl_gauge=T is so far only implemented for some resistivity types.
!
      if (headtt) print*, 'daa_dt: iresistivity=', iresistivity
!
      fres=0.
      etatotal=0.
      diffus_eta=0.; diffus_eta2=0.; diffus_eta3=0.
!
!  Uniform resistivity
!
      eta_const: if (lresi_eta_const) then
        exp_const: if (.not. limplicit_resistivity) then
          if (lweyl_gauge) then
            fres = fres - eta * mu0 * p%jj
          else
            fres = fres + eta * p%del2a
          endif
!
! whatever the gauge is, add an external space-varying electric field
!
          if (ladd_efield) then
             tanhx2 = tanh( x(l1:l2) )*tanh( x(l1:l2) )
             del2aa_ini = ampl_efield*(-2 + 8*tanhx2 - 6*tanhx2*tanhx2 ) 
             fres(:,3) = fres(:,3) - eta*mu0*del2aa_ini
          endif
          if (lfirst .and. ldt) diffus_eta = diffus_eta + eta
        end if exp_const
        etatotal = etatotal + eta
      endif eta_const
!
!  Time-dependent resistivity
!  If both z and t dependent, then use eta_tdep for del2 (in non-Weyl),
!  and -(eta_zdep-1.)*eta_tdep*mu0*p%jj, where eta_zdep < 1 is assumed.
!
      if (lresi_eta_tdep) then
        if (lresi_eta_ztdep) then
          if (lweyl_gauge) then
            fres = fres                 -eta_tdep* feta_ztdep(n)    *mu0*p%jj
          else
            fres = fres+eta_tdep*p%del2a-eta_tdep*(feta_ztdep(n)-1.)*mu0*p%jj
          endif
        else
          if (lweyl_gauge) then
            fres = fres - eta_tdep * mu0 * p%jj
          else
            fres = fres + eta_tdep * p%del2a
          endif
        endif
        if (lfirst .and. ldt) diffus_eta = diffus_eta + eta_tdep
        etatotal = etatotal + eta_tdep
      endif
!
!  z-dependent resistivity
!
      eta_zdep: if (lresi_zdep) then
        exp_zdep: if (.not. limplicit_resistivity) then
          if (lweyl_gauge) then
            fres = fres - eta_z(n) * mu0 * p%jj
          else
            forall(j = 1:3) fres(:,j) = fres(:,j) + eta_z(n) * p%del2a(:,j)
            fres(:,3) = fres(:,3) + geta_z(n) * p%diva
          endif
          if (lfirst .and. ldt) diffus_eta = diffus_eta + eta_z(n)
        else exp_zdep   !MR: What about Weyl gauge here?
          ! Assuming geta_z(:,1) = geta_z(:,2) = 0
          fres(:,3) = fres(:,3) + geta_z(n) * p%diva
          if (lfirst .and. ldt) maxadvec = maxadvec + abs(geta_z(n)) * dz_1(n)
        endif exp_zdep
        etatotal = etatotal + eta_z(n)
      endif eta_zdep
!
      if (lresi_sqrtrhoeta_const) then
        if (lweyl_gauge) then
          do j=1,3
            fres(:,j)=fres(:,j)-eta*sqrt(p%rho1)*mu0*p%jj(:,j)
          enddo
        else
          do j=1,3
            fres(:,j)=fres(:,j)+eta*sqrt(p%rho1)*&
                      (p%del2a(:,j)-0.5*p%diva*p%glnrho(:,j))
          enddo
        endif
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta*sqrt(p%rho1)
        etatotal=etatotal+eta*sqrt(p%rho1)
      endif
!
!  Anisotropic tensor, eta_ij = eta*delta_ij + eta1*qi*qj; see
!  Ruderman & Ruzmaikin (1984) and Plunian & Alboussiere (2020).
!
      if (lresi_eta_aniso) then
        cosalp=cos(alp_aniso*dtor)
        sinalp=sin(alp_aniso*dtor)
        if (eta1_aniso_r==0.) then
          prof=eta1_aniso
        else
          prof=eta1_aniso*(1.-step(x(l1:l2),eta1_aniso_r,eta1_aniso_d))
        endif
        if (lquench_eta_aniso) prof=prof/(1.+quench_aniso*Arms)
        fres(:,1)=fres(:,1)-prof*cosalp*(cosalp*p%jj(:,1)+sinalp*p%jj(:,2))
        fres(:,2)=fres(:,2)-prof*sinalp*(cosalp*p%jj(:,1)+sinalp*p%jj(:,2))
        if (lfirst.and.ldt) diffus_eta=diffus_eta+abs(eta1_aniso)
      endif
!
!  Shakura-Sunyaev type resistivity (mainly just as a demo to show
!  how resistivity can be made depend on temperature.
!  Since etaSS is nonuniform, we use this contribution only for -etaSS*JJ
!  and keep the constant piece with +eta*del2A. (The divA term is eliminated
!  by a suitable gauge transformation.) A sample run is checked in under
!  pencil-runs/1d-tests/bdecay
!
      if (lresi_etaSS) then
        etaSS=alphaSSm*p%cs2/OmegaSS
        do j=1,3
          fres(:,j)=fres(:,j)+eta*p%del2a(:,j)-etaSS*p%jj(:,j)
        enddo
        if (lfirst.and.ldt) diffus_eta=diffus_eta+etaSS+eta
        etatotal=etatotal+etaSS+eta
      endif
!
      if (lresi_xydep) then
        do j=1,3
          fres(:,j)=fres(:,j)+eta_xy(l1:l2,m)*p%del2a(:,j)+geta_xy(l1:l2,m,j)*p%diva
        enddo
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_xy_max
        etatotal=etatotal+eta_xy(l1,m)
      endif
!
      if (lresi_xdep) then
        do j=1,3
          fres(:,j)=fres(:,j)+eta_x(l1:l2)*p%del2a(:,j)
        enddo
        fres(:,1)=fres(:,1)+geta_x(l1:l2)*p%diva
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_x(l1:l2)
        etatotal=etatotal+eta_x(l1:l2)
      endif
!
      if (lresi_rdep) then
        call eta_rdep(eta_r,geta_r,rdep_profile)
        do j=1,3
          fres(:,j)=fres(:,j)+eta_r*p%del2a(:,j)+geta_r(:,j)*p%diva
        enddo
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_r
        etatotal=etatotal+eta_r
      endif
!
      if (lresi_ydep) then
        do j=1,3
          fres(:,j)=fres(:,j)+eta_y(m)*p%del2a(:,j)
        enddo
        fres(:,2)=fres(:,2)+geta_y(m)*p%diva
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_y(m)
        etatotal=etatotal+eta_y(m)
      endif
!
      if (lresi_hyper2) then
        fres=fres+eta_hyper2*p%del4a
        if (lfirst.and.ldt) diffus_eta2=diffus_eta2+eta_hyper2
      endif
!
      if (lresi_hyper3) then
        fres=fres+eta_hyper3*p%del6a
        if (lfirst.and.ldt) diffus_eta3=diffus_eta3+eta_hyper3
      endif
!
      if (lresi_hyper3_polar) then
        do j=1,3
          ju=j+iaa-1
          do i=1,3
            call der6(f,ju,tmp1,i,IGNOREDX=.true.)
            fres(:,j)=fres(:,j)+eta_hyper3*pi4_1*tmp1*dline_1(:,i)**2
          enddo
        enddo
        if (lfirst.and.ldt) &
             diffus_eta3=diffus_eta3+eta_hyper3*pi4_1*dxmin_pencil**4
      endif
!
      if (lresi_hyper3_mesh) then
        do j=1,3
          ju=j+iaa-1
          do i=1,3
            call der6(f,ju,tmp1,i,IGNOREDX=.true.)
            if (ldynamical_diffusion) then
              fres(:,j) = fres(:,j) + eta_hyper3_mesh * tmp1 * dline_1(:,i)
            else
              fres(:,j)=fres(:,j)+eta_hyper3_mesh*pi5_1/60.*tmp1*dline_1(:,i)
            endif
          enddo
        enddo
        if (lfirst.and.ldt) then
          if (ldynamical_diffusion) then
            diffus_eta3 = diffus_eta3 + eta_hyper3_mesh
            advec_hypermesh_aa = 0.0
          else
            advec_hypermesh_aa=eta_hyper3_mesh*pi5_1*sqrt(dxyz_2)
          endif
          advec2_hypermesh=advec2_hypermesh+advec_hypermesh_aa**2
        endif
      endif
!
      if (lresi_hyper3_csmesh) then
        do j=1,3
          ju=j+iaa-1
          do i=1,3
            call der6(f,ju,tmp1,i,IGNOREDX=.true.)
            if (ldynamical_diffusion) then
              fres(:,j)=fres(:,j)+eta_hyper3_mesh*sqrt(p%cs2) &
                       *tmp1*dline_1(:,i)
            else
              fres(:,j)=fres(:,j)+eta_hyper3_mesh*sqrt(p%cs2) &
                       *pi5_1/60.*tmp1*dline_1(:,i)
            endif
          enddo
        enddo
        if (lfirst.and.ldt) then
          if (ldynamical_diffusion) then
            diffus_eta3=diffus_eta3+eta_hyper3_mesh*sqrt(p%cs2)
            advec_hypermesh_aa=0.0
          else
            advec_hypermesh_aa=eta_hyper3_mesh*pi5_1*sqrt(dxyz_2*p%cs2)
          endif
          advec2_hypermesh=advec2_hypermesh+advec_hypermesh_aa**2
        endif
      endif
!
      if (lresi_hyper3_strict) then
        fres=fres+eta_hyper3*f(l1:l2,m,n,ihypres:ihypres+2)
        if (lfirst.and.ldt) diffus_eta3=diffus_eta3+eta_hyper3
      endif
!
      if (lresi_hyper3_aniso) then
         call del6fjv(f,eta_aniso_hyper3,iaa,tmp2)
         fres=fres+tmp2
!  Must divide by dxyz_6 here, because it is multiplied on later.
         if (lfirst.and.ldt) diffus_eta3=diffus_eta3 + &
             (eta_aniso_hyper3(1)*dline_1(:,1)**6 + &
              eta_aniso_hyper3(2)*dline_1(:,2)**6 + &
              eta_aniso_hyper3(3)*dline_1(:,3)**6)/dxyz_6
      endif
!
      if (lresi_shell) then
        call eta_shell(p,eta_mn,geta)
        do j=1,3
          fres(:,j)=fres(:,j)+eta_mn*p%del2a(:,j)+geta(:,j)*p%diva
        enddo
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_mn
        etatotal=etatotal+eta_mn
      endif
!
      if (lresi_eta_shock) then
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i)=fres(:,i)-eta_shock*p%shock*mu0*p%jj(:,i)
          enddo
        else
          do i=1,3
            fres(:,i)=fres(:,i)+ &
                eta_shock*(p%shock*p%del2a(:,i)+p%diva*p%gshock(:,i))
          enddo
        endif
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_shock*p%shock
        etatotal=etatotal+eta_shock*p%shock
      endif
!
      if (lresi_eta_shock2) then
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i)=fres(:,i)-eta_shock2*p%shock**2*mu0*p%jj(:,i)
          enddo
        else
          do i=1,3
            fres(:,i)=fres(:,i)+ &
                eta_shock2*(p%shock**2*p%del2a(:,i)+2*p%shock*p%diva*p%gshock(:,i))
          enddo
        endif
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_shock2*p%shock**2
        etatotal=etatotal+eta_shock2*p%shock**2
      endif
!
! diffusivity: eta-shock with vertical profile
!
      if (lresi_eta_shock_profz) then
        peta_shock = eta_shock + eta_shock_jump1*step(p%z_mn,eta_zshock,-eta_width_shock)
!
! MR: the following only correct in Cartesian geometry!
!
        gradeta_shock(:,1) = 0.
        gradeta_shock(:,2) = 0.
        gradeta_shock(:,3) = eta_shock_jump1*der_step(p%z_mn,eta_zshock,-eta_width_shock)
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i)=fres(:,i)-peta_shock*p%shock*mu0*p%jj(:,i)
          enddo
        else
          do i=1,3
            fres(:,i)=fres(:,i)+ &
                peta_shock*(p%shock*p%del2a(:,i)+p%diva*p%gshock(:,i))+ &
                p%diva*p%shock*gradeta_shock(:,i)
          enddo
        endif
        if (lfirst.and.ldt) diffus_eta=diffus_eta+peta_shock*p%shock
        etatotal=etatotal+peta_shock*p%shock
      endif
!
! diffusivity: eta-shock with radial profile
!
      if (lresi_eta_shock_profr) then
        if (lspherical_coords.or.lsphere_in_a_box) then
          tmp1=p%r_mn
        else
          tmp1=p%rcyl_mn
        endif
        peta_shock = eta_shock + eta_shock_jump1*step(tmp1,eta_xshock,eta_width_shock)
!
        gradeta_shock(:,1) = eta_shock_jump1*der_step(tmp1,eta_xshock,eta_width_shock)
        gradeta_shock(:,2) = 0.
        gradeta_shock(:,3) = 0.
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i)=fres(:,i)-peta_shock*p%shock*mu0*p%jj(:,i)
          enddo
        else
          do i=1,3
            fres(:,i)=fres(:,i)+ &
                peta_shock*(p%shock*p%del2a(:,i)+p%diva*p%gshock(:,i))+ &
                p%diva*p%shock*gradeta_shock(:,i)
          enddo
        endif
        if (lfirst.and.ldt) diffus_eta=diffus_eta+peta_shock*p%shock
        etatotal=etatotal+peta_shock*p%shock
      endif
!
      if (lresi_eta_shock_perp) then
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i)=fres(:,i)-eta_shock*p%shock_perp*mu0*p%jj(:,i)
          enddo
        else
          do i=1,3
            fres(:,i)=fres(:,i)+ &
                eta_shock*(p%shock_perp*p%del2a(:,i)+p%diva*p%gshock_perp(:,i))
          enddo
        endif
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_shock*p%shock_perp
        etatotal=etatotal+eta_shock*p%shock_perp
      endif
!
      if (lresi_etava) then
        if (lweyl_gauge) then
          forall (i = 1:3) fres(:,i) = fres(:,i) - p%etava * p%jj(:,i)
        else
          call fatal_error('daa_dt','eta_va not implemented for resistive gauge')
        endif
        if (lfirst.and.ldt) diffus_eta = diffus_eta + p%etava
        etatotal = etatotal + p%etava
      endif
!
!  Generalised alfven speed dependent resistivity
!
      if (lresi_vAspeed) then
        etatotal = etatotal + p%etava
        if (lweyl_gauge) then
          forall (i = 1:3) fres(:,i) = fres(:,i) - p%etava * p%jj(:,i)
        else
          do i=1,3
            fres(:,i) = fres(:,i) + mu0 * p%etava * p%del2a(:,i) + &
                        eta_va/vArms * p%diva * p%gva(:,i)
          enddo
        endif
        if (lfirst.and.ldt) diffus_eta = diffus_eta + p%etava
      endif
!
      if (lresi_etaj) then
        forall (i = 1:3) fres(:,i) = fres(:,i) - p%etaj * p%jj(:,i)
        if (lfirst.and.ldt) diffus_eta = diffus_eta + p%etaj
        etatotal = etatotal + p%etaj
      endif
!
      if (lresi_etaj2) then
        forall (i = 1:3) fres(:,i) = fres(:,i) - p%etaj2 * p%jj(:,i)
        if (lfirst.and.ldt) diffus_eta = diffus_eta + p%etaj2
        etatotal = etatotal + p%etaj2
      endif
!
      if (lresi_etajrho) then
        forall (i = 1:3) fres(:,i) = fres(:,i) - p%etajrho * p%jj(:,i)
        if (lfirst.and.ldt) diffus_eta = diffus_eta + p%etajrho
        etatotal = etatotal + p%etajrho
      endif
!
      if (lresi_smagorinsky) then
        if (.not.lweyl_gauge) then
          if (letasmag_as_aux) then
             eta_smag=(D_smag*dxmax)**2.*sqrt(p%j2)
             call multsv(eta_smag+eta,p%del2a,fres)
             call grad(f,ietasmag,geta)
!
             do j=1,3 
               fres(:,j)=fres(:,j)+geta(:,j)*p%diva
             enddo
!
          else  
!
!  Term grad(eta_smag) divA not implemented with pencils!
!
            eta_smag=(D_smag*dxmax)**2.*sqrt(p%j2)
            call multsv(eta_smag+eta,p%del2a,fres)
!
          endif
        else
          eta_smag=(D_smag*dxmax)**2.*sqrt(p%j2)
!
          do j=1,3
            fres(:,j)=fres(:,j)-eta_smag*mu0*p%jj(:,j)
          enddo
!
        endif
!
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_smag+eta
        etatotal=etatotal+eta_smag+eta
      endif
!
      if (lresi_smagorinsky_nusmag) then
         eta_smag=Pm_smag1*p%nu_smag
         call multsv(eta_smag+eta,p%del2a,fres)
!
         call grad(f,inusmag,geta)
         do j=1,3 
           fres(:,j)=fres(:,j)+Pm_smag1*geta(:,j)*p%diva
         enddo
!
         if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_smag+eta
         etatotal=etatotal+eta_smag+eta
      endif
!
      if (lresi_smagorinsky_cross) then
        sign_jo=1.
        do i=1,nx
          if (p%jo(i) < 0) sign_jo(i)=-1.
        enddo
        eta_smag=(D_smag*dxmax)**2.*sign_jo*sqrt(p%jo*sign_jo)
        call multsv(eta_smag+eta,p%del2a,fres)
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_smag+eta
        etatotal=etatotal+eta_smag+eta
      endif
!
!  Anomalous resistivity. Sets in when the ion-electron drift speed is
!  larger than some critical value.
!
      if (lresi_anomalous) then
        vdrift=sqrt(sum(p%jj**2,2))*p%rho1
        if (lweyl_gauge) then
          do i=1,3
            if (eta_anom_thresh/=0) then
              where (eta_anom*vdrift > eta_anom_thresh*vcrit_anom) 
                fres(:,i)=fres(:,i)-eta_anom_thresh*mu0*p%jj(:,i)
              elsewhere
                fres(:,i)=fres(:,i)-eta_anom*vdrift/vcrit_anom*mu0*p%jj(:,i)
              endwhere
            else
              where (vdrift>vcrit_anom) &
                fres(:,i)=fres(:,i)-eta_anom*vdrift/vcrit_anom*mu0*p%jj(:,i)
            endif
          enddo
        else
          if (lroot) print*, 'daa_dt: must have Weyl gauge for '// &
              'anomalous resistivity'
          call fatal_error('daa_dt','')
        endif
        if (lfirst.and.ldt) then
          if (eta_anom_thresh/=0) then
            where (eta_anom*vdrift > eta_anom_thresh*vcrit_anom) 
              diffus_eta=diffus_eta+eta_anom_thresh
            elsewhere
              diffus_eta=diffus_eta+eta_anom*vdrift/vcrit_anom
            endwhere
          else
            where (vdrift>vcrit_anom) &
              diffus_eta=diffus_eta+eta_anom*vdrift/vcrit_anom
            endif
        endif
        if (eta_anom_thresh/=0) then
            where (eta_anom*vdrift > eta_anom_thresh*vcrit_anom) 
              etatotal=etatotal+eta_anom_thresh
            elsewhere
              etatotal=etatotal+eta_anom*vdrift/vcrit_anom
            endwhere
        else
          where (vdrift>vcrit_anom) etatotal=etatotal+eta_anom*vdrift/vcrit_anom
        endif
      endif
!
! Temperature dependent resistivity for the solar corona (Spitzer 1969)
!
      if (lresi_spitzer) then
        etatotal = etatotal + eta_spitzer*exp(-1.5*p%lnTT)
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i)=fres(:,i)-eta_spitzer*exp(-1.5*p%lnTT)*mu0*p%jj(:,i)
          enddo
        else
          do i=1,3
            fres(:,i)=fres(:,i)+eta_spitzer*exp(-1.5*p%lnTT)* &
                (p%del2a(:,i)-1.5*p%diva*p%glnTT(:,i))
          enddo
        endif
        if (lfirst.and.ldt) then
          diffus_eta=diffus_eta+eta_spitzer*exp(-1.5*p%lnTT)
        endif
      endif
!
! Resistivity proportional to sound speed for stability of SN Turbulent ISM
! fred: 23.9.17 replaced 0.5 with eta_cspeed so exponent can be generalised
!
      if (lresi_cspeed) then
        etatotal = etatotal + eta*exp(eta_cspeed*p%lnTT)
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i)=fres(:,i)-eta*exp(eta_cspeed*p%lnTT)*mu0*p%jj(:,i)
          enddo
        else
          do i=1,3
            fres(:,i)=fres(:,i)+eta*exp(eta_cspeed*p%lnTT)* &
                (p%del2a(:,i)+0.5*p%diva*p%glnTT(:,i))
          enddo
        endif
        if (lfirst.and.ldt) then
          diffus_eta=diffus_eta+eta*exp(eta_cspeed*p%lnTT)
        endif
      endif
!
! Resistivity proportional to vertical velocity
!
      if (lresi_eta_proptouz) then
        etatotal = etatotal + eta*ampl_eta_uz*p%uu(:,3)
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i)=fres(:,i)-eta*ampl_eta_uz*p%uu(:,3)*mu0*p%jj(:,i)
          enddo
        else
          do i=1,3
            fres(:,i)=fres(:,i)+eta*ampl_eta_uz* &
                 (p%uu(:,3)*p%del2a(:,i)+p%uij(:,3,i)*p%diva)
          enddo
        endif
        if (lfirst.and.ldt) then
          diffus_eta=diffus_eta+eta*ampl_eta_uz*p%uu(:,3)
        endif
      endif
!
! Magnetic field dependent resistivity
!
      if (lresi_magfield) then
        eta_BB(:)=eta/(1.+etaB*p%bb(:,2)**2)
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i)=fres(:,i)-mu0*eta_BB(:)*p%jj(:,i)
          enddo
        else
          call fatal_error('daa_dt','lweyl_gauge=T for lresi_magfield')
        endif
        if (lfirst.and.ldt) then
          call max_mn(eta_BB,maxetaBB)
          diffus_eta=diffus_eta+maxetaBB
        endif
      endif
!
!  anisotropic B-dependent diffusivity
!
      if (eta_aniso_BB/=0.0) then
        where (p%b2==0.)
          tmp1=0.
        elsewhere
          tmp1=eta_aniso_BB/p%b2
        endwhere
        do j=1,3
          ju=j-1+iaa
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-tmp1*p%jb*p%bb(:,j)
        enddo
      endif
!
!  Ambipolar diffusion in the strong coupling approximation.
!
      if (lambipolar_diffusion) then
        do j=1,3
          ju=j-1+iaa
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)+p%nu_ni1*p%jxbrxb(:,j)
        enddo
        if (lentropy .and. lneutralion_heat) then
          if (pretend_lnTT) then
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%cv1*p%TT1*p%nu_ni1*p%jxbr2
          else
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + p%TT1*p%nu_ni1*p%jxbr2
          endif
        elseif (ltemperature .and. lneutralion_heat) then
            df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + p%cv1*p%TT1*p%nu_ni1*p%jxbr2
        endif
        if (lfirst.and.ldt) diffus_eta=diffus_eta+p%nu_ni1*p%va2
      endif
!
!  Consider here the action of a mean friction term, -LLambda*Abar.
!  Works only on one processor. Note that there is also the analogous case
!  to to Ekman friction; see below.
!
      if (lmean_friction) then
        if (nprocxy==1) then
          do j=1,3
            aa_xyaver(:,j)=sum(f(l1:l2,m1:m2,n,j+iax-1))/nxy
          enddo
          dAdt = dAdt-LLambda_aa*aa_xyaver
        else
          call stop_it("magnetic: lmean_friction works only for nprocxy=1")
        endif
      elseif (llocal_friction) then
        dAdt = dAdt-LLambda_aa*p%aa
      endif
!
!   Slope limited diffusion for magnetic field
!
        if (lmagnetic_slope_limited.and.llast) then
!       if (lmagnetic_slope_limited) then
          if (lsld_bb) then
!
!   Using diffusive flux of B on A
!   Idea: DA_i/dt  = ... - e_ikl Dsld_k B_l
!     where Dsld is the SLD operator
!   old way:  DA_i/dt = ... partial_j Dsld_j A_l
!
!
            do j=1,3
              call calc_slope_diff_flux(f,ibx+(j-1),p,h_sld_magn,nlf_sld_magn,tmp1,div_sld_magn, &
                                        FLUX1=d_sld_flux(:,1,j),FLUX2=d_sld_flux(:,2,j),FLUX3=d_sld_flux(:,3,j))
            enddo
!
            tmp2(:,1)= (-d_sld_flux(:,2,3) + d_sld_flux(:,3,2))*fac_sld_magn
            tmp2(:,2)= (-d_sld_flux(:,3,1) + d_sld_flux(:,1,3))*fac_sld_magn
            tmp2(:,3)= (-d_sld_flux(:,1,2) + d_sld_flux(:,2,1))*fac_sld_magn
!
            fres=fres + tmp2
          else
!
            if (lcylindrical_coords .or. lspherical_coords) then
              do j=1,3
                call calc_slope_diff_flux(f,iax+(j-1),p,h_sld_magn,nlf_sld_magn,tmp2(:,j),div_sld_magn, &
                                          FLUX1=d_sld_flux(:,1,j),FLUX2=d_sld_flux(:,2,j),FLUX3=d_sld_flux(:,3,j))
              enddo
!
              if (lcylindrical_coords) then
                fres(:,1)=fres(:,1)+tmp2(:,1)-(d_sld_flux(:,2,2))/x(l1:l2)
                fres(:,2)=fres(:,2)+tmp2(:,2)+(d_sld_flux(:,2,1))/x(l1:l2)
                fres(:,3)=fres(:,3)+tmp2(:,3)
              elseif(lspherical_coords) then
                fres(:,1)=fres(:,1)+tmp2(:,1)-(d_sld_flux(:,2,2)+d_sld_flux(:,3,3))/x(l1:l2)
                fres(:,2)=fres(:,2)+tmp2(:,2)+(d_sld_flux(:,2,1)-d_sld_flux(:,3,3)*cotth(m))/x(l1:l2)
                fres(:,3)=fres(:,3)+tmp2(:,3)+(d_sld_flux(:,3,1)+d_sld_flux(:,3,2)*cotth(m))/x(l1:l2)
              endif
            else
              do j=1,3
                call calc_slope_diff_flux(f,iax+(j-1),p,h_sld_magn,nlf_sld_magn,tmp2(:,j),div_sld_magn)
              enddo
                fres=fres+tmp2
            endif
          endif
!
!     Heating is jj*divF_sld
!     or Heating is just jj*(-e_ijk Dsld_k B_l) (for lsld_bb=T)
!
          if (lohmic_heat) then
            call dot(tmp2,p%jj,tmp1)
            if (lentropy) then
              if (pretend_lnTT) then
                df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + &
                  p%cv1*max(0.0,tmp1)*p%rho1*p%TT1
              else
                df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + &
                        max(0.0,tmp1)*p%rho1*p%TT1
              endif
            else if (ltemperature) then
              if (ltemperature_nolog) then
                df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT) + &
                        p%cv1*max(0.0,tmp1)*p%rho1
              else
                df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + &
                        p%cv1*max(0.0,tmp1)*p%rho1*p%TT1
              endif
            else if (lthermal_energy) then
             df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + max(0.0,tmp1)
            endif
          endif
        endif
!
!  Special contributions to this module are called here.
!
      if (lspecial) call special_calc_magnetic(f,df,p)
! 
! possibility to reduce ohmic heating near the boundary
! currently implemented only for a profile in z above
! a value no_ohmic_heat_z0
! with the width no_ohmic_heat_zwidth
! for reduction, width has to tbe negative.
!
      if (lno_ohmic_heat_bound_z.and.lohmic_heat) &
         etatotal=etatotal*cubic_step(z(n),no_ohmic_heat_z0,no_ohmic_heat_zwidth)
!
!  Add Ohmic heat to entropy or temperature equation.
!
      if (.not.lkinematic.and.lohmic_heat) then
        if (lentropy) then
          if (pretend_lnTT) then
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + &
                p%cv1*etatotal*mu0*p%j2*p%rho1*p%TT1
          else
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + &
                      etatotal*mu0*p%j2*p%rho1*p%TT1
          endif
        else if (ltemperature) then
          if (ltemperature_nolog) then
            df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT) + &
                 p%cv1*etatotal*mu0*p%j2*p%rho1
          else
            df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + &
                 p%cv1*etatotal*mu0*p%j2*p%rho1*p%TT1
          endif
        else if (lthermal_energy) then
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + etatotal*mu0*p%j2
        endif
      endif
!
!  Switch off diffusion in boundary slice if requested by boundconds.
!
!  Only need to do this on bottommost (topmost) processors
!  and in bottommost (topmost) pencils.
!
      do j=1,3
        if (lfrozen_bb_bot(j).and.lfirst_proc_z.and.n==n1) fres(:,j)=0.
        if (lfrozen_bb_top(j).and.llast_proc_z.and.n==n2) fres(:,j)=0.
      enddo
!
!  Induction equation.
!
      if (.not.lupw_aa) then
        if (linduction) then
          if (ladvective_gauge) then
!
!  Take care of possibility of imposed field.
!
            if (any(B_ext/=0.)) then
              if (lfargo_advection) call fatal_error("daadt","fargo advection with external field not tested")
              call cross(p%uu,B_ext,ujiaj)
            else
              if (lfargo_advection) then
                ajiuj=0.
              else
                ujiaj=0.
              endif
            endif
!
!  Calculate ujiaj (=aj uj;i) or ajiuj for fargo advection
!
            do j=1,3
              do k=1,3
                if (lfargo_advection) then
                  ajiuj(:,j)=ajiuj(:,j)+p%uu(:,k)*p%aij(:,k,j)
                else
                  ujiaj(:,j)=ujiaj(:,j)+p%aa(:,k)*p%uij(:,k,j)
                endif
              enddo
            enddo
!
!  Curvature terms on ujiaj
!
            if (lcylindrical_coords) then
              if (lfargo_advection) then
                ajiuj(:,2) = ajiuj(:,2) + (p%aa(:,1)*p%uu(:,2) - p%aa(:,2)*p%uu(:,1))*rcyl_mn1
              else
                ujiaj(:,2) = ujiaj(:,2) + (p%uu(:,1)*p%aa(:,2) - p%uu(:,2)*p%aa(:,1))*rcyl_mn1
              endif
!
            else if (lspherical_coords) then
              if (lfargo_advection) call fatal_error("daadt",&
                   "curvature terms on ajiuj not added for spherical coordinates yet.")
              ujiaj(:,2) = ujiaj(:,2) + (p%uu(:,1)*p%aa(:,2) - p%uu(:,2)*p%aa(:,1))*r1_mn
              ujiaj(:,3) = ujiaj(:,3) + (p%uu(:,1)*p%aa(:,3)          - &
                                         p%uu(:,3)*p%aa(:,1)          + &
                                         p%uu(:,2)*p%aa(:,3)*cotth(m) - &
                                         p%uu(:,3)*p%aa(:,3)*cotth(m))*r1_mn
            endif
!
            if (.not.lfargo_advection) then
              dAdt = dAdt-p%uga-ujiaj+fres
            else
              ! the gauge above, with -ujiaj is unstable due to the buildup of the irrotational term 
              ! Candelaresi et al. 2011. The gauge below does not have the irrotational term. On the 
              ! other hand it cancels out the full advection term if fargo isn't used.
              dAdt = dAdt-p%uuadvec_gaa+ajiuj+fres
            endif
!            df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-p%uga-ujiaj+fres
!
!  ladvective_gauge2
!
          elseif (ladvective_gauge2) then
            if (lua_as_aux) then
              call grad(f,iua,gua)
              dAdt = dAdt + p%uxb+fres-gua
            else
              call fatal_error('daa_dt','must put lua_as_aux=T')
            endif
!
!  ladvective_gauge=F, so just the normal uxb term plus resistive term.
!
          else
            !print*,'this, right?'
            dAdt = dAdt+ p%uxb+fres
          endif
!
!NS: added lnoinduction switch to suppress uxb term when needed
!
           if (lnoinduction) then
             !print*,'no induction'
              dAdt = dAdt - p%uxb
           endif
        endif
      else
!
!  Use upwinding for the advection term.
!
!  We only do upwinding for advection-like terms, u_i f_k,j,
!  for which i=j. This means that for instance in the evolution
!  equation for A_x,
!
!  d(A_x)/dt + u_y A_x,y + u_z A_x,z = u_y A_y,x + u_z A_z,x
!
!  we only do upwinding for the advection-type terms on the
!  left hand side.
!
        if (lupw_aa.and.headtt) then
          print *,'daa_dt: use upwinding in advection term'
        endif
!
!  Add Lorentz force that results from the external field.
!  Note: For now, this only works for uniform external fields.
!
        if (any(B_ext/=0.)) then
          call cross(p%uu,B_ext,uxb_upw)
        else
          uxb_upw=0.
        endif
!
!  Add u_k A_k,j and `upwinded' advection term.
!  Note: this works currently only in cartesian geometry!
!
        if (ladvective_gauge) then
          ujiaj=0.
          do j=1,3
            do k=1,3
              ujiaj(:,j)=ujiaj(:,j)+p%aa(:,k)*p%uij(:,k,j)
            enddo
            uxb_upw(:,j)=uxb_upw(:,j)-p%uga(:,j)-ujiaj(:,j)
          enddo
!
!  ladvective_gauge=F, so just the normal uxb term plus resistive term.
!
        else
          do j=1,3
!
            do k=1,3
              if (k/=j) then
                uxb_upw(:,j)=uxb_upw(:,j)+p%uu(:,k)*(p%aij(:,k,j)-p%aij(:,j,k))
              endif
            enddo
!
            call doupwind(f,iaa+j-1,p%uu,uxb_upw(1,j),mask=j)
!
          enddo
!
        endif
!
!  Full right hand side of the induction equation.
!
        if (linduction) dAdt= dAdt + uxb_upw + fres
      endif
!
!  Add Hall term.
!
      if (hall_term/=0.0) then
        select case (ihall_term)
          case ('const'); hall_term_=hall_term
          case ('t-dep'); hall_term_=hall_term*max(real(t),hall_tdep_t0)**hall_tdep_exponent
          case ('z-dep'); hall_term_=hall_term/(1.-(z(n)-xyz1(3))/Hhall)**hall_zdep_exponent
        endselect
        if (headtt) print*,'daa_dt: hall_term=',hall_term_
        dAdt=dAdt-hall_term_*p%jxb
        if (lfirst.and.ldt) then
          advec_hall=abs(p%uu(:,1)-hall_term_*p%jj(:,1))*dx_1(l1:l2)+ &
                     abs(p%uu(:,2)-hall_term_*p%jj(:,2))*dy_1(  m  )+ &
                     abs(p%uu(:,3)-hall_term_*p%jj(:,3))*dz_1(  n  )
          if (notanumber(advec_hall)) print*, 'advec_hall =',advec_hall
          advec2=advec2+advec_hall**2
          if (headtt.or.ldebug) print*,'daa_dt: max(advec_hall) =',&
                                        maxval(advec_hall)
        endif
      endif
!
!  Add Battery term.
!  corrected by Mikhail Modestov
!
      if (battery_term/=0.0) then
        if (headtt) print*,'daa_dt: battery_term=',battery_term
!---    call multsv_mn(p%rho2,p%fpres,baroclinic)
!AB: corrected by Patrick Adams
        dAdt = dAdt-battery_term*p%fpres
        if (headtt.or.ldebug) print*,'daa_dt: max(battery_term) =',&
!MR; corrected for the time being to fix the auto-test
!            battery_term*maxval(baroclinic)
            battery_term*maxval(p%fpres)
      endif
!
!  Add ambipolar diffusion in strong coupling approximation
!
      if (lambipolar_strong_coupling.and.tauAD/=0.0) then
        if (lfirst.and.ldt) diffus_eta=diffus_eta+tauAD*mu01*p%b2
        dAdt=dAdt+tauAD*p%jxbxb
      endif
!
!  Add jxb/(b^2\nu) magneto-frictional velocity to uxb term
!  Note that this is similar to lambipolar_strong_coupling, but here
!  there is a division by b^2.
!AB: Piyali, I think the mu01 should be removed
!
      if (lmagneto_friction.and.(.not.lhydro).and.numag/=0.0) then
         tmp1=mu01/(numag*(B0_magfric/unit_magnetic**2+p%b2))
         do i=1,3
           dAdt(:,i) = dAdt(:,i) + p%jxbxb(:,i)*tmp1
         enddo
         if (.not. linduction) dAdt = dAdt + fres
      endif
!
!  Possibility of adding extra diffusivity in some halo of given geometry.
!  eta_out is now the diffusivity in the outer halo.
!
      if (height_eta/=0.0) then
        if (headtt) print*,'daa_dt: height_eta,eta_out,lhalox=',height_eta,eta_out,lhalox
        if (lhalox) then
          do ix=1,nx
            tmp=(x(ix)/height_eta)**2
            eta_out1=eta_out*(1.0-exp(-tmp**5/max(1.0-tmp,1.0e-5)))-eta
          enddo
        else
          !eta_out1=eta_out*0.5*(1.-erfunc((z(n)-height_eta)/eta_zwidth))-eta
!AB: 2018-12-18 changed to produce change *above* height_eta.
          eta_out1=(eta_out-eta)*.5*(1.+erfunc((z(n)-height_eta)/eta_zwidth))
        endif
        dAdt = dAdt-(eta_out1*mu0)*p%jj
      endif
!
!  Ekman Friction, used only in two dimensional runs.
!
      if (ekman_friction_aa/=0) &
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-ekman_friction_aa*p%aa
!
!  Add possibility of forcing that is not delta-correlated in time.
!
      if (lforcing_cont_aa) dAdt=dAdt+ ampl_fcont_aa*p%fcont(:,:,iforcing_cont_aa)
!
!  Add possibility of local forcing that is also not delta-correlated in time.
!
      if (lforcing_cont_aa_local) call forcing_continuous(df,p)
!
!  Possibility of relaxation of A in exterior region.
!
      if (tau_aa_exterior/=0.0) call calc_tau_aa_exterior(f,df)
!
!  Relaxing A towards a given profile A_0 on a timescale tau_relprof,
!  note that tau_relprof*u_rms*kf>>1  for this relaxation to affect only the mean fields.
!
      if (tau_relprof/=0.0) then
!        dAdt= dAdt-(p%aa-A_relprof(:,m,n,:))*tau_relprof1
! Piyali: The above is not right as dimension of A_relprof(nx,ny,nz,3), 
! so m,n indices should be the following:
!
        if (lA_relprof_global) then
!
!  use the directly the global external vector potential
!
          dAdt= dAdt-(p%aa-f(l1:l2,m,n,iglobal_ax_ext:iglobal_az_ext))*tau_relprof1
        else
          dAdt= dAdt-(p%aa-A_relprof(:,m-m1+1,n-n1+1,:))*tau_relprof1
        endif
      endif
!
!  Add ``va^2/dx^2'' contribution to timestep.
!  Consider advective timestep only when lhydro=T.
!
      if (lfirst.and.ldt) then
        advec_va2=0.
        if (lhydro) then
          rho1_jxb=p%rho1
          if (rhomin_jxb>0) rho1_jxb=min(rho1_jxb,1/rhomin_jxb)
          if (va2max_jxb>0 .and. (.not. betamin_jxb>0)) then
            rho1_jxb = rho1_jxb &
                     * (1+(p%va2/va2max_jxb)**va2power_jxb)**(-1.0/va2power_jxb)
          endif
          if (betamin_jxb>0) then
            va2max_beta = p%cs2/betamin_jxb*2.0*gamma1
            if (va2max_jxb > 0) va2max_beta=min(va2max_beta,va2max_jxb)
            rho1_jxb = rho1_jxb &
                     * (1+(p%va2/va2max_beta)**va2power_jxb)**(-1.0/va2power_jxb)
          endif
          if (lboris_correction .and. va2max_boris>0) then
            rho1_jxb = rho1_jxb &
                     * (1+(p%va2/va2max_boris)**2.)**(-1.0/2.0)
          endif
          if (lboris_correction .and. cmin>0) then
            rho1_jxb = rho1_jxb &
                     * (1+(p%va2/p%clight2)**2.)**(-1.0/2.0)
          endif
          if (lspherical_coords) then
            advec_va2=((p%bb(:,1)*dx_1(l1:l2))**2+ &
                       (p%bb(:,2)*dy_1(  m  )*r1_mn)**2+ &
                       (p%bb(:,3)*dz_1(  n  )*r1_mn*sin1th(m))**2)*mu01*rho1_jxb
          elseif (lcylindrical_coords) then
            advec_va2=((p%bb(:,1)*dx_1(l1:l2))**2+ &
                       (p%bb(:,2)*dy_1(  m  )*rcyl_mn1)**2+ &
                       (p%bb(:,3)*dz_1(  n  ))**2)*mu01*rho1_jxb
          else
            advec_va2=((p%bb(:,1)*dx_1(l1:l2))**2+ &
                       (p%bb(:,2)*dy_1(  m  ))**2+ &
                       (p%bb(:,3)*dz_1(  n  ))**2)*mu01*rho1_jxb
          endif
        endif
!
!WL: don't know if this is correct, but it's the only way I can make
!    some 1D and 2D samples work when the non-existent direction has the
!    largest velocity (like a 2D rz slice of a Keplerian disk that rotates
!    on the phi direction)
!    Please check
!
        if (lisotropic_advection) then
          if ((nxgrid==1).or.(nygrid==1).or.(nzgrid==1)) &
!            advec_va2=sqrt(p%va2*dxyz_2)
            advec_va2=p%va2*dxyz_2
!JW: no sqrt needed, advec_va2 is quadratic
        endif
!
!mcnallcp: If hall_term is on, the fastest alfven-type mode is the Whistler wave at the grid scale.
! Since the Alfven waves split into the fast whistler mode, the old advec_va2 is not right anymore.
!  This is the generalization for Hall-MHD.
!  This is not used in EMHD simulations.
!
        if (lhydro.and.hall_term/=0.0) then
          if (lcartesian_coords) then
            advec_va2 = ( &
                  (p%bb(:,1)*dx_1(l1:l2)*( &
                    hall_term*pi*dx_1(l1:l2)*mu01 &
                    +sqrt(mu01*p%rho1 + (hall_term*pi*dx_1(l1:l2)*mu01)**2 ) ) )**2 &
                 +(p%bb(:,2)*dy_1(  m  )*( &
                    hall_term*pi*dy_1(  m  )*mu01 &
                    +sqrt(mu01*p%rho1 + (hall_term*pi*dy_1(  m  )*mu01)**2 ) ) )**2 &
                 +(p%bb(:,3)*dz_1(  n  )*( &
                    hall_term*pi*dz_1(  n  )*mu01 &
                    +sqrt(mu01*p%rho1 + (hall_term*pi*dz_1(  n  )*mu01)**2 ) ) )**2 &
                        )
          else
            call fatal_error('daa_dt', 'Timestep advec_va2 with hall_term is not implemented in these coordinates.')
          endif
        endif
        if (notanumber(advec_va2)) print*, 'advec_va2  =',advec_va2
        advec2=advec2+advec_va2
        if (lmagneto_friction) then
          call dot2(p%vmagfric,tmp1)
          advec2=advec2 + tmp1
        endif
      endif
!
!  Apply border profiles.
!
      if (lborder_profiles) call set_border_magnetic(f,df,p)
!
!  Option to constrain timestep for large forces and heat sources to include
!  Lorentz force and Ohmic heating terms 
!  can set in entropy lthdiff_Hmax=F & in hydro lcdt_tauf=F as included here
!
      if (lfirst.and.ldt.and.lrhs_max) then
        do j=1,3
          dAtot=abs(dAdt(:,j))
          dt1_max=max(dt1_max,dAtot/(cdts*alev))
          dAmax=max(dAmax,dAtot/alev)
          ftot=abs(df(l1:l2,m,n,iux+j-1))
          dt1_max=max(dt1_max,ftot/(cdts*alev))
          Fmax=max(Fmax,ftot/alev)
        enddo
        ss0 = abs(df(l1:l2,m,n,iss))
        dt1_max=max(dt1_max,ss0*p%cv1/cdts)
      endif
!
! Electric field E = -dA/dt, store the Electric field in f-array if asked for.
!
      if (lEE_as_aux ) f(l1:l2,m,n,iEEx :iEEz  )= -dAdt
!
!  Magnetic field in spherical coordinates from a Cartesian simulation
!  for sphere-in-a-box setups
!
      if (lbb_sph_as_aux.and.lsphere_in_a_box) then
        f(l1:l2,m,n,ibb_sphr) = p%bb(:,1)*p%evr(:,1)+p%bb(:,2)*p%evr(:,2)+p%bb(:,3)*p%evr(:,3)
        f(l1:l2,m,n,ibb_spht) = p%bb(:,1)*p%evth(:,1)+p%bb(:,2)*p%evth(:,2)+p%bb(:,3)*p%evth(:,3)
        f(l1:l2,m,n,ibb_sphp) = p%bb(:,1)*p%phix+p%bb(:,2)*p%phiy
     endif
!
! Now add all the contribution to dAdt so far into df.
! This is done here, such that contribution from mean-field models are not added to
! the electric field. This may need review later.
!
      df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+dAdt
!
!
!  Call right-hand side for mean-field stuff (do this just before ldiagnos)
!
      if (lmagn_mf) call daa_dt_meanfield(f,df,p)
!
!  Multiply resistivity by Nyquist scale, for resistive time-step.
!
      if (lfirst.and.ldt) then
!
        diffus_eta =diffus_eta *dxyz_2
        diffus_eta2=diffus_eta2*dxyz_4
!
        if (ldynamical_diffusion .and. lresi_hyper3_mesh) then
          diffus_eta3 = diffus_eta3 * (abs(dline_1(:,1)) + abs(dline_1(:,2)) + abs(dline_1(:,3)))
        else
          diffus_eta3 = diffus_eta3*dxyz_6
        endif
        if (ietat/=0) diffus_eta=diffus_eta+maxval(f(l1:l2,m,n,ietat))*dxyz_2
!
        if (headtt.or.ldebug) then
          print*, 'daa_dt: max(diffus_eta)  =', maxval(diffus_eta)
          print*, 'daa_dt: max(diffus_eta2) =', maxval(diffus_eta2)
          print*, 'daa_dt: max(diffus_eta3) =', maxval(diffus_eta3)
        endif

        maxdiffus=max(maxdiffus,diffus_eta)
        maxdiffus2=max(maxdiffus2,diffus_eta2)
        maxdiffus3=max(maxdiffus3,diffus_eta3)
!
      endif

      call calc_diagnostics_magnetic(f,p)
!
!  Debug output.
!
      if (headtt .and. lfirst .and. ip<=4) then
        call output_pencil('aa.dat',p%aa,3)
        call output_pencil('bb.dat',p%bb,3)
        call output_pencil('jj.dat',p%jj,3)
        call output_pencil('del2A.dat',p%del2a,3)
        call output_pencil('JxBr.dat',p%jxbr,3)
        call output_pencil('JxB.dat',p%jxb,3)
        call output_pencil('df.dat',df(l1:l2,m,n,:),mvar)
      endif
      call timing('daa_dt','finished',mnloop=.true.)
!
    endsubroutine daa_dt
!******************************************************************************
    subroutine calc_diagnostics_magnetic(f,p)

      use Diagnostics, only: save_name_sound
      use Slices_methods, only: store_slices
      use Sub, only: cross,vecout

      real, dimension(:,:,:,:) :: f
      type(pencil_case) :: p

      integer :: isound,lspoint,mspoint,nspoint,j
      real, dimension (nx,3) :: uxbxb,poynting
!
! Magnetic field components at the list of points written out in sound.dat
! lwrite_sound is false if either no sound output is required, or if none of
! the desired sound output location occur in the subdomain in this processor.
!
      if (lwrite_sound.and.lout_sound) then
!
! go through the list of lpoints and mpoints
!
        do isound=1,ncoords_sound
!
          lspoint=sound_coords_list(isound,1)
          mspoint=sound_coords_list(isound,2)
          nspoint=sound_coords_list(isound,3)
!
          if ((m==mspoint).and.(n==nspoint)) then
            call save_name_sound(f(lspoint,mspoint,nspoint,iax),idiag_axpt,isound)
            call save_name_sound(f(lspoint,mspoint,nspoint,iay),idiag_aypt,isound)
            call save_name_sound(f(lspoint,mspoint,nspoint,iaz),idiag_azpt,isound)
            call save_name_sound(p%bb(lspoint-nghost,1),idiag_bxpt,isound)
            call save_name_sound(p%bb(lspoint-nghost,2),idiag_bypt,isound)
            call save_name_sound(p%bb(lspoint-nghost,3),idiag_bzpt,isound)
            call save_name_sound(p%jj(lspoint-nghost,1),idiag_jxpt,isound)
            call save_name_sound(p%jj(lspoint-nghost,2),idiag_jypt,isound)
            call save_name_sound(p%jj(lspoint-nghost,3),idiag_jzpt,isound)
            call save_name_sound(uxbb(lspoint-nghost,1),idiag_Expt,isound)
            call save_name_sound(uxbb(lspoint-nghost,2),idiag_Eypt,isound)
            call save_name_sound(uxbb(lspoint-nghost,3),idiag_Ezpt,isound)
          endif
        enddo
      endif

      call calc_2d_diagnostics_magnetic(p)
      call calc_1d_diagnostics_magnetic(p)
      if (ldiagnos) call calc_0d_diagnostics_magnetic(f,p)

!
!  Write B-slices for output in wvid in run.f90.
!  Note: ix_loc is the index with respect to array with ghost zones.
!
      if (lvideo.and.lfirst) then
!
        if (ivid_aps/=0) call store_slices(p%aps,aps_xy,aps_xz,aps_yz,xz2=aps_xz2)
        if (ivid_bb/=0) call store_slices(p%bb,bb_xy,bb_xz,bb_yz,bb_xy2,bb_xy3,bb_xy4,bb_xz2)
        if (ivid_jj/=0) call store_slices(p%jj,jj_xy,jj_xz,jj_yz,jj_xy2,jj_xy3,jj_xy4,jj_xz2)
        if (ivid_b2/=0) call store_slices(p%b2,b2_xy,b2_xz,b2_yz,b2_xy2,b2_xy3,b2_xy4,b2_xz2)
        if (ivid_j2/=0) call store_slices(p%j2,j2_xy,j2_xz,j2_yz,j2_xy2,j2_xy3,j2_xy4,j2_xz2)
        if (ivid_bb_sph/=0) call store_slices(p%bb_sph,bb_sph_xy,bb_sph_xz,bb_sph_yz,bb_sph_xy2,bb_sph_xy3,bb_sph_xy4,bb_sph_xz2)
        if (ivid_jb/=0) call store_slices(p%jb,jb_xy,jb_xz,jb_yz,jb_xy2,jb_xy3,jb_xy4,jb_xz2)
        if (ivid_ab/=0) call store_slices(p%ab,ab_xy,ab_xz,ab_yz,ab_xy2,ab_xy3,ab_xy4,ab_xz2)
        if (ivid_beta1/=0) call store_slices(p%beta1,beta1_xy,beta1_xz,beta1_yz,beta1_xy2, &
                                             beta1_xy3,beta1_xy4,beta1_xz2)
        if (.not.lgpu.and.ivid_poynting/=0) then
          call cross(p%uxb,p%bb,uxbxb)
          do j=1,3
            poynting(:,j) = etatotal*p%jxb(:,j) - mu01*uxbxb(:,j)  !!! etatotal invalid outside daa_dt!
          enddo
          call store_slices(poynting,poynting_xy,poynting_xz,poynting_yz, &
                            poynting_xy2,poynting_xy3,poynting_xy4,poynting_xz2)
        endif
!
        if (bthresh_per_brms/=0) then
          call calc_bthresh
          call vecout(41,trim(directory)//'/bvec',p%bb,bthresh,nbvec)
        endif
!
      endif

    endsubroutine calc_diagnostics_magnetic
!******************************************************************************
    subroutine calc_0d_diagnostics_magnetic(f,p)
!
!  Calculate diagnostic quantities.
!
      use Diagnostics
      use Sub, only: dot, dot2, multvs,cross_mn, multmv_transp,dot2_mn,cross_mn, dot2_mn, dot_mn_sv

      real, dimension(:,:,:,:) :: f
      type(pencil_case) :: p

      real, dimension (nx,3) :: exj,dexb,phib,jxbb,uxDxuxb
      real, dimension (nx) :: uxj_dotB0,b3b21,b3b12,b1b32,b1b23,b2b13,b2b31
      real, dimension (nx) :: jxb_dotB0,uxb_dotB0
      real, dimension (nx) :: oxuxb_dotB0,jxbxb_dotB0,uxDxuxb_dotB0
      real, dimension (nx) :: aj,tmp1,fres2
      real, dimension (nx) :: B1dot_glnrhoxb,fb,fxbx
      real, dimension (nx) :: b2t,bjt,jbt
      real, dimension (nx) :: uj,phi,dub,dob,jdel2a,epsAD

      call sum_mn_name(p%beta1,idiag_beta1m)
      call max_mn_name(p%beta1,idiag_beta1max)
      call sum_mn_name(p%beta, idiag_betam)
      call max_mn_name(p%beta, idiag_betamax)

      if (idiag_betamin /= 0) call max_mn_name(-p%beta, idiag_betamin, lneg=.true.)

      if (.not.lgpu) then
!
!  These diagnostics rely upon mn-dependent quantities which are not in the pencil case.
!
        if (lrhs_max) then
          if (idiag_dtHr/=0) call max_mn_name(ss0*p%cv1/cdts,idiag_dtHr,l_dt=.true.)
          if (idiag_dtFr/=0) call max_mn_name(Fmax/cdts,idiag_dtFr,l_dt=.true.)
          if (idiag_dtBr/=0) call max_mn_name(dAmax/cdts,idiag_dtBr,l_dt=.true.)
        endif
!
!  Integrate velocity in time, to calculate correlation time later.
!
        if (idiag_b2tm/=0) then
          call dot(p%bb,f(l1:l2,m,n,ibxt:ibzt),b2t)
          call sum_mn_name(b2t,idiag_b2tm)
        endif
!
!  Integrate velocity in time, to calculate correlation time later.
!
        if (idiag_jbtm/=0) then
          call dot(p%jj,f(l1:l2,m,n,ibxt:ibzt),jbt)
          call sum_mn_name(jbt,idiag_jbtm)
        endif
!
!  Integrate velocity in time, to calculate correlation time later.
!
        if (idiag_bjtm/=0) then
          call dot(p%bb,f(l1:l2,m,n,ijxt:ijzt),bjt)
          call sum_mn_name(bjt,idiag_bjtm)
        endif

        if (idiag_Bresrms/=0 .or. idiag_Rmrms/=0) then
          call dot2_mn(fres,fres2)
          call sum_mn_name(fres2,idiag_Bresrms,lsqrt=.true.)
          if (idiag_Rmrms/=0) &
              call sum_mn_name(p%uxb2/fres2,idiag_Rmrms,lsqrt=.true.)
        endif
      endif
!
!  Contributions to vertical Poynting vector. Consider them here
!  separately, because the contribution b^2*uz may be a good indicator
!  of magnetic buoyancy.
!
      if (idiag_b2ruzm/=0) call sum_mn_name(p%b2*p%rho*p%uu(:,3),idiag_b2ruzm)
      if (idiag_b2uzm/=0) call sum_mn_name(p%b2*p%uu(:,3),idiag_b2uzm)
      if (idiag_ubbzm/=0) call sum_mn_name(p%ub*p%bb(:,3),idiag_ubbzm)
!
!  Mean squared and maximum squared magnetic field.
!
      if (idiag_b1m/=0) call sum_mn_name(sqrt(p%b2),idiag_b1m)
      call sum_mn_name(p%b2,idiag_b2m)
      if (idiag_EEM/=0) call sum_mn_name(.5*p%b2,idiag_EEM)
      if (idiag_b4m/=0) call sum_mn_name(p%b2**2,idiag_b4m)
      call max_mn_name(p%b2,idiag_bm2)
      call sum_mn_name(p%b2,idiag_brms,lsqrt=.true.)
      call sum_mn_name(p%bf2,idiag_bfrms,lsqrt=.true.)
      call sum_mn_name(p%bf2,idiag_bf2m)
      if (idiag_bf4m/=0) call sum_mn_name(p%bf2**2,idiag_bf4m)
      if (idiag_emag/=0) call integrate_mn_name(mu012*p%b2,idiag_emag)
      if (idiag_brmsh/=0) then
        if (lequatory) call sum_mn_name_halfy(p%b2,idiag_brmsh)
        if (lequatorz) call sum_mn_name_halfz(p%b2,idiag_brmsh)
        fname(idiag_brmsn)=fname_half(idiag_brmsh,1)
        fname(idiag_brmss)=fname_half(idiag_brmsh,2)
        itype_name(idiag_brmsn)=ilabel_sum_sqrt
        itype_name(idiag_brmss)=ilabel_sum_sqrt
      endif
      if (idiag_brmsx/=0) call sum_mn_name(p%b2*xmask_mag,idiag_brmsx,lsqrt=.true.)
      if (idiag_brmsz/=0) call sum_mn_name(p%b2*zmask_mag(n-n1+1),idiag_brmsz,lsqrt=.true.)
      call max_mn_name(p%b2,idiag_bmax,lsqrt=.true.)
      if (idiag_bxmin/=0) call max_mn_name(-p%bb(:,1),idiag_bxmin,lneg=.true.)
      if (idiag_bymin/=0) call max_mn_name(-p%bb(:,2),idiag_bymin,lneg=.true.)
      if (idiag_bzmin/=0) call max_mn_name(-p%bb(:,3),idiag_bzmin,lneg=.true.)
      call max_mn_name(p%bb(:,1),idiag_bxmax)
      call max_mn_name(p%bb(:,2),idiag_bymax)
      call max_mn_name(p%bb(:,3),idiag_bzmax)
      call max_mn_name(abs(p%bbb(:,1)),idiag_bbxmax)
      call max_mn_name(abs(p%bbb(:,2)),idiag_bbymax)
      call max_mn_name(abs(p%bbb(:,3)),idiag_bbzmax)
      call max_mn_name(abs(p%jj(:,1)),idiag_jxmax)
      call max_mn_name(abs(p%jj(:,2)),idiag_jymax)
      call max_mn_name(abs(p%jj(:,3)),idiag_jzmax)
      if (idiag_aybym2/=0) &
          call sum_mn_name(2.*p%aa(:,2)*p%bb(:,2),idiag_aybym2)
      call sum_mn_name(p%ab,idiag_abm)
      if (idiag_abumx/=0) call sum_mn_name(p%uu(:,1)*p%ab,idiag_abumx)
      if (idiag_abumy/=0) call sum_mn_name(p%uu(:,2)*p%ab,idiag_abumy)
      if (idiag_abumz/=0) call sum_mn_name(p%uu(:,3)*p%ab,idiag_abumz)
      if (idiag_abrms/=0) call sum_mn_name(p%ab**2,idiag_abrms,lsqrt=.true.)
      if (idiag_jbrms/=0) call sum_mn_name(p%jb**2,idiag_jbrms,lsqrt=.true.)
!
!  Hemispheric magnetic helicity of total field.
!  North means 1 and south means 2.
!
      if (idiag_abmh/=0) then
        if (lequatory) call sum_mn_name_halfy(p%ab,idiag_abmh)
        if (lequatorz) call sum_mn_name_halfz(p%ab,idiag_abmh)
        fname(idiag_abmn)=fname_half(idiag_abmh,1)
        fname(idiag_abms)=fname_half(idiag_abmh,2)
        itype_name(idiag_abmn)=ilabel_sum
        itype_name(idiag_abms)=ilabel_sum
      endif
!
!  Hemispheric current helicity of total field.
!  North means 1 and south means 2.
!
      if (idiag_jbmh/=0) then
        if (lequatory) call sum_mn_name_halfy(p%jb,idiag_jbmh)
        if (lequatorz) call sum_mn_name_halfz(p%jb,idiag_jbmh)
        fname(idiag_jbmn)=fname_half(idiag_jbmh,1)
        fname(idiag_jbms)=fname_half(idiag_jbmh,2)
        itype_name(idiag_jbmn)=ilabel_sum
        itype_name(idiag_jbms)=ilabel_sum
      endif
!
!  Mean dot product of forcing and magnetic field, <f.b>.
!
      if (idiag_fbm/=0) then
        call dot(p%fcont(:,:,iforcing_cont_aa),p%bb,fb)
        call sum_mn_name(fb,idiag_fbm)
      endif
!
      if (idiag_fxbxm/=0) then
        fxbx=p%fcont(:,1,iforcing_cont_aa)*p%bb(:,1)
        call sum_mn_name(fxbx,idiag_fxbxm)
      endif
!
!  Cross helicity (linkage between vortex tubes and flux tubes).
!
      call sum_mn_name(p%ub,idiag_ubm)
      if (idiag_uxbxm/=0) call sum_mn_name(p%uu(:,1)*p%bb(:,1),idiag_uxbxm)
      if (idiag_uybxm/=0) call sum_mn_name(p%uu(:,2)*p%bb(:,1),idiag_uybxm)
      if (idiag_uzbxm/=0) call sum_mn_name(p%uu(:,3)*p%bb(:,1),idiag_uzbxm)
      if (idiag_uxbym/=0) call sum_mn_name(p%uu(:,1)*p%bb(:,2),idiag_uxbym)
      if (idiag_uybym/=0) call sum_mn_name(p%uu(:,2)*p%bb(:,2),idiag_uybym)
      if (idiag_uzbym/=0) call sum_mn_name(p%uu(:,3)*p%bb(:,2),idiag_uzbym)
      if (idiag_uxbzm/=0) call sum_mn_name(p%uu(:,1)*p%bb(:,3),idiag_uxbzm)
      if (idiag_uybzm/=0) call sum_mn_name(p%uu(:,2)*p%bb(:,3),idiag_uybzm)
      if (idiag_uzbzm/=0) call sum_mn_name(p%uu(:,3)*p%bb(:,3),idiag_uzbzm)
      call sum_mn_name(p%cosub,idiag_cosubm)
!
!  Current helicity tensor (components)
!
      !if (idiag_jbm/=0) call sum_mn_name(p%ub,idiag_ubm)
      if (idiag_jxbxm/=0) call sum_mn_name(p%jj(:,1)*p%bb(:,1),idiag_jxbxm)
      if (idiag_jybxm/=0) call sum_mn_name(p%jj(:,2)*p%bb(:,1),idiag_jybxm)
      if (idiag_jzbxm/=0) call sum_mn_name(p%jj(:,3)*p%bb(:,1),idiag_jzbxm)
      if (idiag_jxbym/=0) call sum_mn_name(p%jj(:,1)*p%bb(:,2),idiag_jxbym)
      if (idiag_jybym/=0) call sum_mn_name(p%jj(:,2)*p%bb(:,2),idiag_jybym)
      if (idiag_jzbym/=0) call sum_mn_name(p%jj(:,3)*p%bb(:,2),idiag_jzbym)
      if (idiag_jxbzm/=0) call sum_mn_name(p%jj(:,1)*p%bb(:,3),idiag_jxbzm)
      if (idiag_jybzm/=0) call sum_mn_name(p%jj(:,2)*p%bb(:,3),idiag_jybzm)
      if (idiag_jzbzm/=0) call sum_mn_name(p%jj(:,3)*p%bb(:,3),idiag_jzbzm)
!
!  Velocity-current density tensor (components)
!
      if (idiag_uxjxm/=0) call sum_mn_name(p%uu(:,1)*p%jj(:,1),idiag_uxjxm)
      if (idiag_uxjym/=0) call sum_mn_name(p%uu(:,1)*p%jj(:,2),idiag_uxjym)
      if (idiag_uxjzm/=0) call sum_mn_name(p%uu(:,1)*p%jj(:,3),idiag_uxjzm)
      if (idiag_uyjxm/=0) call sum_mn_name(p%uu(:,2)*p%jj(:,1),idiag_uyjxm)
      if (idiag_uyjym/=0) call sum_mn_name(p%uu(:,2)*p%jj(:,2),idiag_uyjym)
      if (idiag_uyjzm/=0) call sum_mn_name(p%uu(:,2)*p%jj(:,3),idiag_uyjzm)
      if (idiag_uzjxm/=0) call sum_mn_name(p%uu(:,3)*p%jj(:,1),idiag_uzjxm)
      if (idiag_uzjym/=0) call sum_mn_name(p%uu(:,3)*p%jj(:,2),idiag_uzjym)
      if (idiag_uzjzm/=0) call sum_mn_name(p%uu(:,3)*p%jj(:,3),idiag_uzjzm)
!
!  compute rms value of difference between u and b    !!!MR: units?
!
      if (idiag_dubrms/=0) then
        call dot2(p%uu-p%bb,dub)
        call sum_mn_name(dub,idiag_dubrms,lsqrt=.true.)
      endif
!
!  compute rms value of difference between vorticity and b   !!!MR: units?
!
      if (idiag_dobrms/=0) then
        call dot2(p%oo-p%bb,dob)
        call sum_mn_name(dob,idiag_dobrms,lsqrt=.true.)
      endif
!
!  Field-velocity cross helicity (linkage between velocity and magnetic tubes).
!
      call sum_mn_name(p%ua,idiag_uam)
!
!  Current-vortex cross helicity (linkage between vortex and current tubes).
!
      if (idiag_ujm/=0) then
        call dot(p%uu,p%jj,uj)
        call sum_mn_name(uj,idiag_ujm)
      endif
!
!  Mean field <B_i>, and mean components of the correlation matrix <B_i B_j>.
!  Note that this quantity does not include any imposed field!
!
      call sum_mn_name(p%bbb(:,1),idiag_bxm)
      call sum_mn_name(p%bbb(:,2),idiag_bym)
      call sum_mn_name(p%bbb(:,3),idiag_bzm)
      if (idiag_bx2m/=0) call sum_mn_name(p%bbb(:,1)**2,idiag_bx2m)
      if (idiag_by2m/=0) call sum_mn_name(p%bbb(:,2)**2,idiag_by2m)
      if (idiag_bz2m/=0) call sum_mn_name(p%bbb(:,3)**2,idiag_bz2m)
      if (idiag_bx4m/=0) call sum_mn_name(p%bbb(:,1)**4,idiag_bx4m)
      if (idiag_by4m/=0) call sum_mn_name(p%bbb(:,2)**4,idiag_by4m)
      if (idiag_bz4m/=0) call sum_mn_name(p%bbb(:,3)**4,idiag_bz4m)
      if (idiag_jx2m/=0) call sum_mn_name(p%jj(:,1)**2,idiag_jx2m)
      if (idiag_jy2m/=0) call sum_mn_name(p%jj(:,2)**2,idiag_jy2m)
      if (idiag_jz2m/=0) call sum_mn_name(p%jj(:,3)**2,idiag_jz2m)
      if (idiag_jx4m/=0) call sum_mn_name(p%jj(:,1)**4,idiag_jx4m)
      if (idiag_jy4m/=0) call sum_mn_name(p%jj(:,2)**4,idiag_jy4m)
      if (idiag_jz4m/=0) call sum_mn_name(p%jj(:,3)**4,idiag_jz4m)
      if (idiag_bxbym/=0) call sum_mn_name(p%bbb(:,1)*p%bbb(:,2),idiag_bxbym)
      if (idiag_bxbzm/=0) call sum_mn_name(p%bbb(:,1)*p%bbb(:,3),idiag_bxbzm)
      if (idiag_bybzm/=0) call sum_mn_name(p%bbb(:,2)*p%bbb(:,3),idiag_bybzm)
      call sum_mn_name(p%djuidjbi,idiag_djuidjbim)
!
!  Calculate B*sin(phi) = -<Bx*sinkz> + <By*coskz>.
!
      if (idiag_bsinphz/=0) call sum_mn_name(-p%bbb(:,1)*sinkz(n)+p%bbb(:,2)*coskz(n),idiag_bsinphz)
!
!  Calculate B*cos(phi) = <Bx*coskz> + <By*sinkz>.
!
      if (idiag_bcosphz/=0) call sum_mn_name(+p%bbb(:,1)*coskz(n)+p%bbb(:,2)*sinkz(n),idiag_bcosphz)
!
!  v_A = |B|/sqrt(rho); in units where mu_0=1
!
      call sum_mn_name(p%va2,idiag_vA2m)
      call sum_mn_name(p%va2,idiag_vArms,lsqrt=.true.)
      call max_mn_name(p%va2,idiag_vAmax,lsqrt=.true.)
      if (.not.lgpu) then
        if (idiag_dtb/=0) &
          call max_mn_name(sqrt(advec_va2)/cdt,idiag_dtb,l_dt=.true.)
      endif
!
!  Lorentz force.
!
      call sum_mn_name(p%jxbr(:,1),idiag_jxbrxm)
      call sum_mn_name(p%jxbr(:,2),idiag_jxbrym)
      call sum_mn_name(p%jxbr(:,3),idiag_jxbrzm)
      call sum_mn_name(p%jxbr2,idiag_jxbr2m)
      call max_mn_name(p%jxbr2,idiag_jxbrmax,lsqrt=.true.)
!
!  <J.A> for calculating k_effective, for example.
!
      if (idiag_ajm/=0) then
        call dot(p%aa,p%jj,aj)
        call sum_mn_name(aj,idiag_ajm)
      endif
!
!  Helicity integrals.
!
      call integrate_mn_name(p%ab,idiag_ab_int)
      call integrate_mn_name(p%jb,idiag_jb_int)
!
! <J.B>
!
      call sum_mn_name(p%jb,idiag_jbm)
      call sum_mn_name(p%hjb,idiag_hjbm)
      call sum_mn_name(p%j2,idiag_j2m)
      call max_mn_name(p%j2,idiag_jm2)
      call sum_mn_name(p%j2,idiag_jrms,lsqrt=.true.)
      call sum_mn_name(p%hj2,idiag_hjrms,lsqrt=.true.)
      call max_mn_name(p%j2,idiag_jmax,lsqrt=.true.)
      if (.not.lgpu) then
        if (idiag_epsM_LES/=0) call sum_mn_name(eta_smag*p%j2,idiag_epsM_LES)
        if (idiag_dteta/=0)  call max_mn_name(diffus_eta/cdtv,idiag_dteta,l_dt=.true.)
        if (idiag_dteta3/=0)  call max_mn_name(diffus_eta3/cdtv3,idiag_dteta3,l_dt=.true.)
      endif
      call sum_mn_name(p%cosjb,idiag_cosjbm)
      call sum_mn_name(p%coshjb,idiag_coshjbm)
      call sum_mn_name(p%jparallel,idiag_jparallelm)
      call sum_mn_name(p%jperp,idiag_jperpm)
      call sum_mn_name(p%hjparallel,idiag_hjparallelm)
      call sum_mn_name(p%hjperp,idiag_hjperpm)
!
!  Resistivity.
!
      if (.not.lgpu) then
        call sum_mn_name(eta_smag,idiag_etasmagm)
        if (idiag_etasmagmin/=0) call max_mn_name(-eta_smag,idiag_etasmagmin,lneg=.true.)
        call max_mn_name(eta_smag,idiag_etasmagmax)
        call save_name(eta1_aniso/(1.+quench_aniso*Arms),idiag_etaaniso)
      endif
      call max_mn_name(p%etava,idiag_etavamax)
      call max_mn_name(p%etaj,idiag_etajmax)
      call max_mn_name(p%etaj2,idiag_etaj2max)
      call max_mn_name(p%etajrho,idiag_etajrhomax)
!
!  Not correct for hyperresistivity:
!
      if (.not.lgpu) then
        if (idiag_epsM/=0) call sum_mn_name(etatotal*mu0*p%j2,idiag_epsM)
      endif
!
!  Heating by ion-neutrals friction.
!
      if (idiag_epsAD/=0) then
        if (lambipolar_strong_coupling.and.tauAD/=0.0) then
          call dot(p%jj,p%jxbxb,epsAD)
          call sum_mn_name(-tauAD*epsAD,idiag_epsAD)
        else
          call sum_mn_name(p%nu_ni1*p%rho*p%jxbr2,idiag_epsAD)
        endif
      endif
!
!  <A>'s, <A^2> and A^2|max
!
      call sum_mn_name(p%aa(:,1),idiag_axm)
      call sum_mn_name(p%aa(:,2),idiag_aym)
      call sum_mn_name(p%aa(:,3),idiag_azm)
      call sum_mn_name(p%a2,idiag_a2m)
      call sum_mn_name(p%a2,idiag_arms,lsqrt=.true.)
      call max_mn_name(p%a2,idiag_amax,lsqrt=.true.)
!
!  Divergence of A
!
      if (idiag_divarms /= 0) call sum_mn_name(p%diva**2, idiag_divarms, lsqrt=.true.)
!
!  Calculate surface integral <2ExA>*dS.
!
      if (idiag_exaym2/=0) call helflux(p%aa,p%uxb,p%jj)
!
!  Calculate surface integral <2ExJ>*dS.
!
      if (idiag_exjm2/=0) call curflux(p%uxb,p%jj)
!--     if (idiag_exjm2/=0) call curflux_dS(p%uxb,p%jj)
!
!  Calculate emf for alpha effect (for imposed field).
!  Note that uxbm means <EMF.B0>/B0^2, so it gives already alpha=EMF/B0.
!
      if (idiag_uxbm/=0 .or. idiag_uxbmx/=0 .or. idiag_uxbmy/=0 &
          .or. idiag_uxbcmx/=0 .or. idiag_uxbcmy/=0 &
          .or. idiag_uxbsmx/=0 .or. idiag_uxbsmy/=0 &
          .or. idiag_uxbmz/=0) then
        if (idiag_uxbm/=0) then
          call dot(B_ext_inv,p%uxb,uxb_dotB0)
          call sum_mn_name(uxb_dotB0,idiag_uxbm)
        endif
        call sum_mn_name(uxbb(:,1),idiag_uxbmx)
        call sum_mn_name(uxbb(:,2),idiag_uxbmy)
        call sum_mn_name(uxbb(:,3),idiag_uxbmz)
        if (idiag_uxbcmx/=0) call sum_mn_name(uxbb(:,1)*coskz(n),idiag_uxbcmx)
        if (idiag_uxbcmy/=0) call sum_mn_name(uxbb(:,2)*coskz(n),idiag_uxbcmy)
        if (idiag_uxbsmx/=0) call sum_mn_name(uxbb(:,1)*sinkz(n),idiag_uxbsmx)
        if (idiag_uxbsmy/=0) call sum_mn_name(uxbb(:,2)*sinkz(n),idiag_uxbsmy)
      endif
!
!  Calculate part I of magnetic helicity flux (ExA contribution).
!
      if (idiag_examx/=0 .or. idiag_examy/=0 .or. idiag_examz/=0 .or. &
          idiag_exatop/=0 .or. idiag_exabot/=0) then
        call sum_mn_name(p%exa(:,1),idiag_examx)
        call sum_mn_name(p%exa(:,2),idiag_examy)
        call sum_mn_name(p%exa(:,3),idiag_examz)
!
        if (idiag_exabot/=0) then
          if (n==n1.and.lfirst_proc_z) call integrate_mn_name(p%exa(:,3),idiag_exabot)
        endif
!
        if (idiag_exatop/=0) then
          if (n==n2.and.llast_proc_z) call integrate_mn_name(p%exa(:,3),idiag_exatop)
        endif
!
      endif
!
!  Calculate part II of magnetic helicity flux (phi*B contribution).
!
      if (idiag_phibmx/=0 .or. idiag_phibmy/=0 .or. idiag_phibmz/=0) then
        if (lweyl_gauge) then
          phi=0.
        elseif (ladvective_gauge) then
          phi=p%ua
        else
          phi=eta*p%diva
        endif
        call multvs(p%bb,phi,phib)
        call sum_mn_name(phib(:,1),idiag_phibmx)
        call sum_mn_name(phib(:,2),idiag_phibmy)
        call sum_mn_name(phib(:,3),idiag_phibmz)
      endif
!
!  Calculate part I of current helicity flux (for imposed field).
!
      if (idiag_exjmx/=0 .or. idiag_exjmy/=0 .or. idiag_exjmz/=0) then
        call cross_mn(-p%uxb+eta*p%jj,p%jj,exj)
        call sum_mn_name(exj(:,1),idiag_exjmx)
        call sum_mn_name(exj(:,2),idiag_exjmy)
        call sum_mn_name(exj(:,3),idiag_exjmz)
      endif
!
!  Calculate part II of current helicity flux (for imposed field).
!  < curlE x B >|_i  =  < B_{j,i} E_j >
!  Use the full B (with B_ext)
!
      if (idiag_dexbmx/=0 .or. idiag_dexbmy/=0 .or. idiag_dexbmz/=0) then
        call multmv_transp(p%bij,-p%uxb+eta*p%jj,dexb)
        call sum_mn_name(dexb(:,1),idiag_dexbmx)
        call sum_mn_name(dexb(:,2),idiag_dexbmy)
        call sum_mn_name(dexb(:,3),idiag_dexbmz)
      endif
!
!  Calculate <uxj>.B0/B0^2.
!
      if (idiag_uxjm/=0) then
        call dot(B_ext_inv,p%uxj,uxj_dotB0)
        call sum_mn_name(uxj_dotB0,idiag_uxjm)
      endif
!
!  Calculate <u x B>_rms, <resistive terms>_rms, <ratio ~ Rm>_rms.
!
      call sum_mn_name(p%uxb2,idiag_uxBrms,lsqrt=.true.)
!
!  Calculate <b^2*divu>, which is part of <u.(jxb)>.
!  Note that <u.(jxb)>=1/2*<b^2*divu>+<u.bgradb>.
!
      if (idiag_b2divum/=0) call sum_mn_name(p%b2*p%divu,idiag_b2divum)
!
!  Calculate <J.del2a>.
!
      if (idiag_jdel2am/=0) then
        call dot(p%jj,p%del2a,jdel2a)
        call sum_mn_name(jdel2a,idiag_jdel2am)
      endif
!
!  Calculate <u.(jxb)>.
!
      call sum_mn_name(p%ujxb,idiag_ujxbm)
!
!  Calculate <jxb>.B_0/B_0^2.
!
      if (idiag_jxbm/=0) then
        call dot(B_ext_inv,p%jxb,jxb_dotB0)
        call sum_mn_name(jxb_dotB0,idiag_jxbm)
      endif
      if (idiag_jxbmx/=0.or.idiag_jxbmy/=0.or.idiag_jxbmz/=0) then
        call cross_mn(p%jj,p%bbb,jxbb)
        call sum_mn_name(jxbb(:,1),idiag_jxbmx)
        call sum_mn_name(jxbb(:,2),idiag_jxbmy)
        call sum_mn_name(jxbb(:,3),idiag_jxbmz)
      endif
      call max_mn_name(p%vmagfric,idiag_vmagfricmax)
      if (idiag_vmagfricrms/=0) then
        call dot2_mn(p%vmagfric,tmp1)
        call sum_mn_name(tmp1,idiag_vmagfricrms,lsqrt=.true.)
      endif
!
!  Magnetic triple correlation term (for imposed field).
!
      if (idiag_jxbxbm/=0) then
        call dot(B_ext_inv,p%jxbxb,jxbxb_dotB0)
        call sum_mn_name(jxbxb_dotB0,idiag_jxbxbm)
      endif
!
!  Triple correlation from Reynolds tensor (for imposed field).
!
      if (idiag_oxuxbm/=0) then
        call dot(B_ext_inv,p%oxuxb,oxuxb_dotB0)
        call sum_mn_name(oxuxb_dotB0,idiag_oxuxbm)
      endif
!
!  Triple correlation from pressure gradient (for imposed field).
!  (assume cs2=1, and that no entropy evolution is included)
!  This is ok for all applications currently under consideration.
!
      if (idiag_gpxbm/=0) then
        call dot_mn_sv(B1_ext,p%glnrhoxb,B1dot_glnrhoxb)
        call sum_mn_name(B1dot_glnrhoxb,idiag_gpxbm)
      endif
!
!  < u x curl(uxB) > = < E_i u_{j,j} - E_j u_{j,i} >
!   ( < E_1 u2,2 + E1 u3,3 - E2 u2,1 - E3 u3,1 >
!     < E_2 u1,1 + E2 u3,3 - E1 u2,1 - E3 u3,2 >
!     < E_3 u1,1 + E3 u2,2 - E1 u3,1 - E2 u2,3 > )
!
      if (idiag_uxDxuxbm/=0) then
        uxDxuxb(:,1)=p%uxb(:,1)*(p%uij(:,2,2)+p%uij(:,3,3))-p%uxb(:,2)*p%uij(:,2,1)-p%uxb(:,3)*p%uij(:,3,1)
        uxDxuxb(:,2)=p%uxb(:,2)*(p%uij(:,1,1)+p%uij(:,3,3))-p%uxb(:,1)*p%uij(:,1,2)-p%uxb(:,3)*p%uij(:,3,2)
        uxDxuxb(:,3)=p%uxb(:,3)*(p%uij(:,1,1)+p%uij(:,2,2))-p%uxb(:,1)*p%uij(:,1,3)-p%uxb(:,2)*p%uij(:,2,3)
        call dot(B_ext_inv,uxDxuxb,uxDxuxb_dotB0)
        call sum_mn_name(uxDxuxb_dotB0,idiag_uxDxuxbm)
      endif
!
!  alpM11=<b3*b2,1>
!
      if (idiag_b3b21m/=0) then
        b3b21=p%bb(:,3)*p%bij(:,2,1)
        call sum_mn_name(b3b21,idiag_b3b21m)
      endif
!
!  alpM11=<b3*b1,2>
!
      if (idiag_b3b12m/=0) then
        b3b12=p%bb(:,3)*p%bij(:,1,2)
        call sum_mn_name(b3b12,idiag_b3b12m)
      endif
!
!  alpM22=<b1*b3,2>
!
      if (idiag_b1b32m/=0) then
        b1b32=p%bb(:,1)*p%bij(:,3,2)
        call sum_mn_name(b1b32,idiag_b1b32m)
      endif
!
!  alpM22=<b1*b2,3>
!
      if (idiag_b1b23m/=0) then
        b1b23=p%bb(:,1)*p%bij(:,2,3)
        call sum_mn_name(b1b23,idiag_b1b23m)
      endif
!
!  alpM33=<b2*b1,3>
!
      if (idiag_b2b13m/=0) then
        b2b13=p%bb(:,2)*p%bij(:,1,3)
        call sum_mn_name(b2b13,idiag_b2b13m)
      endif
!
!  alpM33=<b2*b3,1>
!
      if (idiag_b2b31m/=0) then
        b2b31=p%bb(:,2)*p%bij(:,3,1)
        call sum_mn_name(b2b31,idiag_b2b31m)
      endif
!
!  current density components at one point (=pt).
!
      if (lroot.and.m==mpoint.and.n==npoint) then  
        !MR: i.e., only pointwise data from root proc domain can be obtained!
        if (idiag_bxpt/=0) call save_name(p%bb(lpoint-nghost,1),idiag_bxpt)
        if (idiag_bypt/=0) call save_name(p%bb(lpoint-nghost,2),idiag_bypt)
        if (idiag_bzpt/=0) call save_name(p%bb(lpoint-nghost,3),idiag_bzpt)
        if (idiag_jxpt/=0) call save_name(p%jj(lpoint-nghost,1),idiag_jxpt)
        if (idiag_jypt/=0) call save_name(p%jj(lpoint-nghost,2),idiag_jypt)
        if (idiag_jzpt/=0) call save_name(p%jj(lpoint-nghost,3),idiag_jzpt)
      endif
!
!  current density components at point 2 (=p2).
!
      if (lroot.and.m==mpoint2.and.n==npoint2) then
        if (idiag_bxp2/=0) call save_name(p%bb(lpoint2-nghost,1),idiag_bxp2)
        if (idiag_byp2/=0) call save_name(p%bb(lpoint2-nghost,2),idiag_byp2)
        if (idiag_bzp2/=0) call save_name(p%bb(lpoint2-nghost,3),idiag_bzp2)
        if (idiag_jxp2/=0) call save_name(p%jj(lpoint2-nghost,1),idiag_jxp2)
        if (idiag_jyp2/=0) call save_name(p%jj(lpoint2-nghost,2),idiag_jyp2)
        if (idiag_jzp2/=0) call save_name(p%jj(lpoint2-nghost,3),idiag_jzp2)
      endif
!
    endsubroutine calc_0d_diagnostics_magnetic
!******************************************************************************
    subroutine calc_1d_diagnostics_magnetic(p)
!
!  2-D averages.
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      use Diagnostics
      use Sub, only: dot2_mn

      type(pencil_case) :: p

      real, dimension(nx) :: fres2, tmp1, Rmmz 
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1.
!
      if (l1davgfirst .or. (ldiagnos .and. ldiagnos_need_zaverages)) then

        call yzsum_mn_name_x(p%b2, idiag_b2mx)
        call yzsum_mn_name_x(p%bb(:,1),idiag_bxmx)
        call yzsum_mn_name_x(p%bb(:,2),idiag_bymx)
        call yzsum_mn_name_x(p%bb(:,3),idiag_bzmx)
        if (idiag_bx2mx/=0) call yzsum_mn_name_x(p%bb(:,1)**2,idiag_bx2mx)
        if (idiag_by2mx/=0) call yzsum_mn_name_x(p%bb(:,2)**2,idiag_by2mx)
        if (idiag_bz2mx/=0) call yzsum_mn_name_x(p%bb(:,3)**2,idiag_bz2mx)
        if (idiag_bxbymx/=0) call yzsum_mn_name_x(p%bbb(:,1)*p%bbb(:,2),idiag_bxbymx)
        if (idiag_bxbzmx/=0) call yzsum_mn_name_x(p%bbb(:,1)*p%bbb(:,3),idiag_bxbzmx)
        if (idiag_bybzmx/=0) call yzsum_mn_name_x(p%bbb(:,2)*p%bbb(:,3),idiag_bybzmx)
        call yzsum_mn_name_x(p%beta, idiag_betamx)
        call xzsum_mn_name_y(p%bb(:,1),idiag_bxmy)
        call xzsum_mn_name_y(p%bb(:,2),idiag_bymy)
        call xzsum_mn_name_y(p%bb(:,3),idiag_bzmy)
        if (idiag_beta2mx/=0) call yzsum_mn_name_x(p%beta**2, idiag_beta2mx)
        if (idiag_bx2my/=0) call xzsum_mn_name_y(p%bb(:,1)**2,idiag_bx2my)
        if (idiag_by2my/=0) call xzsum_mn_name_y(p%bb(:,2)**2,idiag_by2my)
        if (idiag_bz2my/=0) call xzsum_mn_name_y(p%bb(:,3)**2,idiag_bz2my)
        call xysum_mn_name_z(p%aa(:,1),idiag_axmz)
        call xysum_mn_name_z(p%aa(:,2),idiag_aymz)
        call xysum_mn_name_z(p%aa(:,3),idiag_azmz)
        if (idiag_abuxmz/=0) call xysum_mn_name_z(p%ab*p%uu(:,1),idiag_abuxmz)
        if (idiag_abuymz/=0) call xysum_mn_name_z(p%ab*p%uu(:,2),idiag_abuymz)
        if (idiag_abuzmz/=0) call xysum_mn_name_z(p%ab*p%uu(:,3),idiag_abuzmz)
        if (idiag_uabxmz/=0) call xysum_mn_name_z(p%ua*p%bb(:,1),idiag_uabxmz)
        if (idiag_uabymz/=0) call xysum_mn_name_z(p%ua*p%bb(:,2),idiag_uabymz)
        if (idiag_uabzmz/=0) call xysum_mn_name_z(p%ua*p%bb(:,3),idiag_uabzmz)
        call xysum_mn_name_z(p%bbb(:,1),idiag_bbxmz)
        call xysum_mn_name_z(p%bbb(:,2),idiag_bbymz)
        call xysum_mn_name_z(p%bbb(:,3),idiag_bbzmz)
        call xysum_mn_name_z(p%bb(:,1),idiag_bxmz)
        call xysum_mn_name_z(p%bb(:,2),idiag_bymz)
        call xysum_mn_name_z(p%bb(:,3),idiag_bzmz)
        call xysum_mn_name_z(p%jj(:,1),idiag_jxmz)
        call xysum_mn_name_z(p%jj(:,2),idiag_jymz)
        call xysum_mn_name_z(p%jj(:,3),idiag_jzmz)
        call xysum_mn_name_z(p%uxb(:,1),idiag_Exmz)
        call xysum_mn_name_z(p%uxb(:,2),idiag_Eymz)
        call xysum_mn_name_z(p%uxb(:,3),idiag_Ezmz)
        if (idiag_bx2mz/=0) call xysum_mn_name_z(p%bb(:,1)**2,idiag_bx2mz)
        if (idiag_by2mz/=0) call xysum_mn_name_z(p%bb(:,2)**2,idiag_by2mz)
        if (idiag_bz2mz/=0) call xysum_mn_name_z(p%bb(:,3)**2,idiag_bz2mz)
        if (idiag_bx2rmz/=0) call xysum_mn_name_z(p%bb(:,1)**2*p%rho1,idiag_bx2rmz)
        if (idiag_by2rmz/=0) call xysum_mn_name_z(p%bb(:,2)**2*p%rho1,idiag_by2rmz)
        if (idiag_bz2rmz/=0) call xysum_mn_name_z(p%bb(:,3)**2*p%rho1,idiag_bz2rmz)
        if (idiag_beta2mz/=0) call xysum_mn_name_z(p%beta**2, idiag_beta2mz)
        call xysum_mn_name_z(p%beta1,idiag_beta1mz)
        call xysum_mn_name_z(p%beta, idiag_betamz)
        call xysum_mn_name_z(p%jb,idiag_jbmz)
        call xysum_mn_name_z(p%d6ab,idiag_d6abmz)
        call xysum_mn_name_z(p%del6a(:,1),idiag_d6amz1)
        call xysum_mn_name_z(p%del6a(:,2),idiag_d6amz2)
        call xysum_mn_name_z(p%del6a(:,3),idiag_d6amz3)
        call xysum_mn_name_z(p%ab,idiag_abmz)
        call xysum_mn_name_z(p%ub,idiag_ubmz)
        call xysum_mn_name_z(p%ua,idiag_uamz)
        call xysum_mn_name_z(p%diva,idiag_divamz)
        call xysum_mn_name_z(p%bb(:,3)*p%diva,idiag_bzdivamz)
        if (idiag_uxbxmz/=0) call xysum_mn_name_z(p%uu(:,1)*p%bb(:,1),idiag_uxbxmz)
        if (idiag_uybxmz/=0) call xysum_mn_name_z(p%uu(:,2)*p%bb(:,1),idiag_uybxmz)
        if (idiag_uzbxmz/=0) call xysum_mn_name_z(p%uu(:,3)*p%bb(:,1),idiag_uzbxmz)
        if (idiag_uxbymz/=0) call xysum_mn_name_z(p%uu(:,1)*p%bb(:,2),idiag_uxbymz)
        if (idiag_uybymz/=0) call xysum_mn_name_z(p%uu(:,2)*p%bb(:,2),idiag_uybymz)
        if (idiag_uzbymz/=0) call xysum_mn_name_z(p%uu(:,3)*p%bb(:,2),idiag_uzbymz)
        if (idiag_uxbzmz/=0) call xysum_mn_name_z(p%uu(:,1)*p%bb(:,3),idiag_uxbzmz)
        if (idiag_uybzmz/=0) call xysum_mn_name_z(p%uu(:,2)*p%bb(:,3),idiag_uybzmz)
        if (idiag_uzbzmz/=0) call xysum_mn_name_z(p%uu(:,3)*p%bb(:,3),idiag_uzbzmz)
        if (.not.lgpu) then
          if (idiag_epsMmz/=0) call xysum_mn_name_z(etatotal*mu0*p%j2,idiag_epsMmz)
          call yzsum_mn_name_x(etatotal,idiag_etatotalmx)
          call xysum_mn_name_z(etatotal,idiag_etatotalmz)
        endif
        if (idiag_vmagfricmz/=0) then
          call dot2_mn(p%vmagfric,tmp1)
          call xysum_mn_name_z(tmp1,idiag_vmagfricmz)
        endif
!
!  Calculate magnetic helicity flux (ExA contribution).
!
        call xysum_mn_name_z(p%exa(:,1),idiag_examz1)
        call xysum_mn_name_z(p%exa(:,2),idiag_examz2)
        call xysum_mn_name_z(p%exa(:,3),idiag_examz3)
!
!  Calculate magnetic helicity flux for n=3 hyperdiffusion (E^{(3)}xA contribution).
!
        call xysum_mn_name_z(p%e3xa(:,1),idiag_e3xamz1)
        call xysum_mn_name_z(p%e3xa(:,2),idiag_e3xamz2)
        call xysum_mn_name_z(p%e3xa(:,3),idiag_e3xamz3)
!
!  Maxwell stress components.
!
        if (idiag_bxbymy/=0) call xzsum_mn_name_y(p%bbb(:,1)*p%bbb(:,2),idiag_bxbymy)
        if (idiag_bxbzmy/=0) call xzsum_mn_name_y(p%bbb(:,1)*p%bbb(:,3),idiag_bxbzmy)
        if (idiag_bybzmy/=0) call xzsum_mn_name_y(p%bbb(:,2)*p%bbb(:,3),idiag_bybzmy)
        if (idiag_bxbymz/=0) call xysum_mn_name_z(p%bbb(:,1)*p%bbb(:,2),idiag_bxbymz)
        if (idiag_bxbzmz/=0) call xysum_mn_name_z(p%bbb(:,1)*p%bbb(:,3),idiag_bxbzmz)
        if (idiag_bybzmz/=0) call xysum_mn_name_z(p%bbb(:,2)*p%bbb(:,3),idiag_bybzmz)
        call yzsum_mn_name_x(p%jxbr(:,1),idiag_jxbrxmx)
        call yzsum_mn_name_x(p%jxbr(:,2),idiag_jxbrymx)
        call yzsum_mn_name_x(p%jxbr(:,3),idiag_jxbrzmx)
        call xzsum_mn_name_y(p%jxbr(:,1),idiag_jxbrxmy)
        call xzsum_mn_name_y(p%jxbr(:,2),idiag_jxbrymy)
        call xzsum_mn_name_y(p%jxbr(:,3),idiag_jxbrzmy)
        call xysum_mn_name_z(p%jxbr(:,1),idiag_jxbrxmz)
        call xysum_mn_name_z(p%jxbr(:,2),idiag_jxbrymz)
        call xysum_mn_name_z(p%jxbr(:,3),idiag_jxbrzmz)
        call xysum_mn_name_z(p%b2,idiag_b2mz)
        call xysum_mn_name_z(p%bf2,idiag_bf2mz)
        call xysum_mn_name_z(p%j2,idiag_j2mz)
        if (.not.lgpu) then
          if (idiag_poynzmz/=0) call xysum_mn_name_z(etatotal*p%jxb(:,3)-mu01* &
            (p%uxb(:,1)*p%bb(:,2)-p%uxb(:,2)*p%bb(:,1)),idiag_poynzmz)
        endif
        call phizsum_mn_name_r(p%b2,idiag_b2mr)
        if (idiag_brmr/=0)   &
             call phizsum_mn_name_r(p%bb(:,1)*p%pomx+p%bb(:,2)*p%pomy,idiag_brmr)
        if (idiag_bpmr/=0)   &
             call phizsum_mn_name_r(p%bb(:,1)*p%phix+p%bb(:,2)*p%phiy,idiag_bpmr)
        call phizsum_mn_name_r(p%bb(:,3),idiag_bzmr)
        if (idiag_armr/=0)   &
             call phizsum_mn_name_r(p%aa(:,1)*p%pomx+p%aa(:,2)*p%pomy,idiag_armr)
        if (idiag_apmr/=0)   &
             call phizsum_mn_name_r(p%aa(:,1)*p%phix+p%aa(:,2)*p%phiy,idiag_apmr)
        call phizsum_mn_name_r(p%aa(:,3),idiag_azmr)
        call yzintegrate_mn_name_x(p%bb(:,1),idiag_mflux_x)
        call xzintegrate_mn_name_y(p%bb(:,2),idiag_mflux_y)
        call xyintegrate_mn_name_z(p%bb(:,3),idiag_mflux_z)

        if (.not.lgpu) then
!
!  This diagnostic relies upon mn-dependent quantities which are not in the pencil case.
!
          if (idiag_Rmmz/=0) then
            call dot2_mn(fres,fres2)
            Rmmz=sqrt(p%uxb2/fres2)
            where (fres2 < tini) Rmmz = 0.
            call xysum_mn_name_z(Rmmz,idiag_Rmmz)
          endif
        endif
      endif

    endsubroutine calc_1d_diagnostics_magnetic
!***********************************************************************
    subroutine calc_2d_diagnostics_magnetic(p)
!
!  2-D averages.
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      use Diagnostics

      type(pencil_case) :: p

      real, dimension(nx,3) :: tmp2

      if (l2davgfirst) then
        if (idiag_brmphi/=0) call phisum_mn_name_rz(p%bb(:,1)*p%pomx+p%bb(:,2)*p%pomy,&
                                                    idiag_brmphi)
        if (idiag_brsphmphi/=0) call phisum_mn_name_rz(p%bb(:,1)*p%evr(:,1)+&
            p%bb(:,2)*p%evr(:,2)+p%bb(:,3)*p%evr(:,3),idiag_brsphmphi)
        if (idiag_bthmphi/=0) call phisum_mn_name_rz(p%bb(:,1)*p%evth(:,1)+&
            p%bb(:,2)*p%evth(:,2)+p%bb(:,3)*p%evth(:,3),idiag_bthmphi)
        if (idiag_bpmphi/=0) call phisum_mn_name_rz(p%bb(:,1)*p%phix+p%bb(:,2)*p%phiy,&
                                                    idiag_bpmphi)
        call phisum_mn_name_rz(p%bb(:,3),idiag_bzmphi)
        call phisum_mn_name_rz(p%b2,idiag_b2mphi)
        call phisum_mn_name_rz(p%jb,idiag_jbmphi)
        if (any((/idiag_uxbrmphi,idiag_uxbpmphi,idiag_uxbzmphi/) /= 0)) then
          if (idiag_uxbrmphi/=0) call phisum_mn_name_rz(p%uxb(:,1)*p%pomx+p%uxb(:,2)*p%pomy,idiag_uxbrmphi)
          if (idiag_uxbpmphi/=0) call phisum_mn_name_rz(p%uxb(:,1)*p%phix+p%uxb(:,2)*p%phiy,idiag_uxbpmphi)
          call phisum_mn_name_rz(p%uxb(:,3),idiag_uxbzmphi)
        endif
        if (any((/idiag_jxbrmphi,idiag_jxbpmphi,idiag_jxbzmphi/) /= 0)) then
          if (idiag_jxbrmphi/=0) call phisum_mn_name_rz(p%jxb(:,1)*p%pomx+p%jxb(:,2)*p%pomy,idiag_jxbrmphi)
          if (idiag_jxbpmphi/=0) call phisum_mn_name_rz(p%jxb(:,1)*p%phix+p%jxb(:,2)*p%phiy,idiag_jxbpmphi)
          call phisum_mn_name_rz(p%jxb(:,3),idiag_jxbzmphi)
        endif
        if (any((/idiag_armphi,idiag_apmphi,idiag_azmphi/) /= 0)) then
          if (idiag_armphi/=0) call phisum_mn_name_rz(p%aa(:,1)*p%pomx+p%aa(:,2)*p%pomy,idiag_armphi)
          if (idiag_apmphi/=0) call phisum_mn_name_rz(p%aa(:,1)*p%phix+p%aa(:,2)*p%phiy,idiag_apmphi)
          call phisum_mn_name_rz(p%aa(:,3),idiag_azmphi)
        endif
        call zsum_mn_name_xy(p%bb(:,1),idiag_bxmxy)
        call zsum_mn_name_xy(p%bb,idiag_bymxy,(/0,1,0/))
        call zsum_mn_name_xy(p%bb,idiag_bzmxy,(/0,0,1/))
        call zsum_mn_name_xy(p%jj(:,1),idiag_jxmxy)
        call zsum_mn_name_xy(p%jj,idiag_jymxy,(/0,1,0/))
        call zsum_mn_name_xy(p%jj,idiag_jzmxy,(/0,0,1/))
        call zsum_mn_name_xy(p%aa(:,1),idiag_axmxy)
        call zsum_mn_name_xy(p%aa,idiag_aymxy,(/0,1,0/))
        call zsum_mn_name_xy(p%aa,idiag_azmxy,(/0,0,1/))
        if (lcovariant_magnetic) then
          if (idiag_dbxdxmxy/=0) call zsum_mn_name_xy(p%bijtilde(:,1,1)+p%bij_cov_corr(:,1,1),idiag_dbxdxmxy)
          if (idiag_dbxdymxy/=0) call zsum_mn_name_xy(p%bijtilde(:,1,2)+p%bij_cov_corr(:,1,2),idiag_dbxdymxy)
          if (idiag_dbxdzmxy/=0) call zsum_mn_name_xy(p%bijtilde(:,1,3)+p%bij_cov_corr(:,1,3),idiag_dbxdzmxy)
          if (idiag_dbydxmxy/=0) call zsum_mn_name_xy(p%bijtilde(:,2,1)+p%bij_cov_corr(:,2,1),idiag_dbydxmxy)
          if (idiag_dbydymxy/=0) call zsum_mn_name_xy(p%bijtilde(:,2,2)+p%bij_cov_corr(:,2,2),idiag_dbydymxy)
          if (idiag_dbydzmxy/=0) call zsum_mn_name_xy(p%bijtilde(:,2,3)+p%bij_cov_corr(:,2,3),idiag_dbydzmxy)
          if (idiag_dbzdxmxy/=0) call zsum_mn_name_xy(p%bijtilde(:,3,1)+p%bij_cov_corr(:,3,1),idiag_dbzdxmxy)
          if (idiag_dbzdymxy/=0) call zsum_mn_name_xy(p%bijtilde(:,3,2)+p%bij_cov_corr(:,3,2),idiag_dbzdymxy)
          if (idiag_dbzdzmxy/=0) call zsum_mn_name_xy(p%bijtilde(:,3,3)+p%bij_cov_corr(:,3,3),idiag_dbzdzmxy)
          !if (idiag_dbxdxmxy/=0)&
          !  call zsum_mn_name_xy(p%bijtilde(:,1,1)+p%bij_cov_corr(:,1,1)-p%bij(:,1,1),idiag_dbxdxmxy)
          !if (idiag_dbxdymxy/=0)&
          !  call zsum_mn_name_xy(p%bijtilde(:,1,2)+p%bij_cov_corr(:,1,2)-p%bij(:,1,2),idiag_dbxdymxy)
          !if (idiag_dbxdzmxy/=0)&
          !  call zsum_mn_name_xy(p%bijtilde(:,1,3)+p%bij_cov_corr(:,1,3)-p%bij(:,1,3),idiag_dbxdzmxy)
          !if (idiag_dbydxmxy/=0)&
          !  call zsum_mn_name_xy(p%bijtilde(:,2,1)+p%bij_cov_corr(:,2,1)-p%bij(:,2,1),idiag_dbydxmxy)
          !if (idiag_dbydymxy/=0)&
          !  call zsum_mn_name_xy(p%bijtilde(:,2,2)+p%bij_cov_corr(:,2,2)-p%bij(:,2,2),idiag_dbydymxy)
          !if (idiag_dbydzmxy/=0)&
          !  call zsum_mn_name_xy(p%bijtilde(:,2,3)+p%bij_cov_corr(:,2,3)-p%bij(:,2,3),idiag_dbydzmxy)
          !if (idiag_dbzdxmxy/=0)&
          !  call zsum_mn_name_xy(p%bijtilde(:,3,1)+p%bij_cov_corr(:,3,1)-p%bij(:,3,1),idiag_dbzdxmxy)
          !if (idiag_dbzdymxy/=0)& 
          !  call zsum_mn_name_xy(p%bijtilde(:,3,2)+p%bij_cov_corr(:,3,2)-p%bij(:,3,2),idiag_dbzdymxy)
          !if (idiag_dbzdzmxy/=0)& 
          !  call zsum_mn_name_xy(p%bijtilde(:,3,3)+p%bij_cov_corr(:,3,3)-p%bij(:,3,3),idiag_dbzdzmxy)

        else
          call zsum_mn_name_xy(p%bijtilde(:,1,1),idiag_dbxdxmxy)
          call zsum_mn_name_xy(p%bijtilde(:,1,2),idiag_dbxdymxy)
          call zsum_mn_name_xy(p%bijtilde(:,2,1),idiag_dbydxmxy)
          call zsum_mn_name_xy(p%bijtilde(:,2,2),idiag_dbydymxy)
          call zsum_mn_name_xy(p%bijtilde(:,3,1),idiag_dbzdxmxy)
          call zsum_mn_name_xy(p%bijtilde(:,3,2),idiag_dbzdymxy)
        endif
!
        call ysum_mn_name_xz(p%b2,idiag_b2mxz)
        call ysum_mn_name_xz(p%aa(:,1),idiag_axmxz)
        call ysum_mn_name_xz(p%aa(:,2),idiag_aymxz)
        call ysum_mn_name_xz(p%aa(:,3),idiag_azmxz)
        if (idiag_bx1mxz/=0) call ysum_mn_name_xz(abs(p%bb(:,1)),idiag_bx1mxz)
        if (idiag_by1mxz/=0) call ysum_mn_name_xz(abs(p%bb(:,2)),idiag_by1mxz)
        if (idiag_bz1mxz/=0) call ysum_mn_name_xz(abs(p%bb(:,3)),idiag_bz1mxz)
        call ysum_mn_name_xz(p%bb(:,1),idiag_bxmxz)
        call ysum_mn_name_xz(p%bb(:,2),idiag_bymxz)
        call ysum_mn_name_xz(p%bb(:,3),idiag_bzmxz)
        call ysum_mn_name_xz(p%jj(:,1),idiag_jxmxz)
        call ysum_mn_name_xz(p%jj(:,2),idiag_jymxz)
        call ysum_mn_name_xz(p%jj(:,3),idiag_jzmxz)
        if (idiag_bx2mxz/=0) call ysum_mn_name_xz(p%bb(:,1)**2,idiag_bx2mxz)
        if (idiag_by2mxz/=0) call ysum_mn_name_xz(p%bb(:,2)**2,idiag_by2mxz)
        if (idiag_bz2mxz/=0) call ysum_mn_name_xz(p%bb(:,3)**2,idiag_bz2mxz)
!
        if (idiag_bx2mxy/=0) call zsum_mn_name_xy(p%bb(:,1)**2,idiag_bx2mxy)
        call zsum_mn_name_xy(p%bb,idiag_by2mxy,(/0,2,0/))
        call zsum_mn_name_xy(p%bb,idiag_bz2mxy,(/0,0,2/))
        call zsum_mn_name_xy(p%jb,idiag_jbmxy)
        call zsum_mn_name_xy(p%ab,idiag_abmxy)
        call zsum_mn_name_xy(p%ub,idiag_ubmxy)
        call zsum_mn_name_xy(p%exa(:,1),idiag_examxy1)
        call zsum_mn_name_xy(p%exa,idiag_examxy2,(/0,1,0/))
        call zsum_mn_name_xy(p%exa,idiag_examxy3,(/0,0,1/))
        call zsum_mn_name_xy(p%uxb(:,1),idiag_Exmxy)
        call zsum_mn_name_xy(p%uxb,idiag_Eymxy,(/0,1,0/))
        call zsum_mn_name_xy(p%uxb,idiag_Ezmxy,(/0,0,1/))
        if (.not.lgpu) then
          if (idiag_poynxmxy/=0) &
              call zsum_mn_name_xy(etatotal*p%jxb(:,1)-mu01* &
              (p%uxb(:,2)*p%bb(:,3)-p%uxb(:,3)*p%bb(:,2)),idiag_poynxmxy)
          if (idiag_poynymxy/=0.or.idiag_poynzmxy/=0) then
            tmp2(:,1)=0.
            tmp2(:,2)=etatotal*p%jxb(:,2)-mu01*(p%uxb(:,3)*p%bb(:,1)-p%uxb(:,1)*p%bb(:,3))
            tmp2(:,3)=etatotal*p%jxb(:,3)-mu01*(p%uxb(:,1)*p%bb(:,2)-p%uxb(:,2)*p%bb(:,1))
          endif
          call zsum_mn_name_xy(tmp2,idiag_poynymxy,(/0,1,0/))
          call zsum_mn_name_xy(tmp2,idiag_poynzmxy,(/0,0,1/))
        endif
        call zsum_mn_name_xy(p%beta1,idiag_beta1mxy)
!
! Stokes parameters correct for Yin-Yang?
!
        call zsum_mn_name_xy(p%StokesI,idiag_StokesImxy)
        call zsum_mn_name_xy(p%StokesQ,idiag_StokesQmxy)
        call zsum_mn_name_xy(p%StokesU,idiag_StokesUmxy)
        call zsum_mn_name_xy(p%StokesQ1,idiag_StokesQ1mxy)
        call zsum_mn_name_xy(p%StokesU1,idiag_StokesU1mxy)
        call zsum_mn_name_xy(p%bb,idiag_bxbymxy,(/1,1,0/))
        call zsum_mn_name_xy(p%bb,idiag_bxbzmxy,(/1,0,1/))
        call zsum_mn_name_xy(p%bb,idiag_bybzmxy,(/0,1,1/))
!
        if (idiag_bxbymxz/=0) &
            call ysum_mn_name_xz(p%bb(:,1)*p%bb(:,2),idiag_bxbymxz)
        if (idiag_bxbzmxz/=0) &
            call ysum_mn_name_xz(p%bb(:,1)*p%bb(:,3),idiag_bxbzmxz)
        if (idiag_bybzmxz/=0) &
            call ysum_mn_name_xz(p%bb(:,2)*p%bb(:,3),idiag_bybzmxz)
        if (idiag_uybxmxz/=0) &
            call ysum_mn_name_xz(p%uu(:,2)*p%bb(:,1),idiag_uybxmxz)
        if (idiag_uybzmxz/=0) &
            call ysum_mn_name_xz(p%uu(:,2)*p%bb(:,3),idiag_uybzmxz)
        call ysum_mn_name_xz(p%uxb(:,1),idiag_Exmxz)
        call ysum_mn_name_xz(p%uxb(:,2),idiag_Eymxz)
        call ysum_mn_name_xz(p%uxb(:,3),idiag_Ezmxz)
        call ysum_mn_name_xz(p%va2,idiag_vAmxz)
      else  !MR: Why else?
!
!  idiag_bxmxy and idiag_bymxy also need to be calculated when
!  ldiagnos and idiag_bmx and/or idiag_bmy, so
!
!  We may need to calculate bxmxy without calculating bmx. The following
!  if condition was messing up calculation of bmxy_rms
!
        if (ldiagnos) then
          call zsum_mn_name_xy(p%bb(:,1),idiag_bxmxy)
          call zsum_mn_name_xy(p%bb,idiag_bymxy,(/0,1,0/))
          call zsum_mn_name_xy(p%bb,idiag_bzmxy,(/0,0,1/))
          call zsum_mn_name_xy(p%jj(:,1),idiag_jxmxy)
          call zsum_mn_name_xy(p%jj,idiag_jymxy,(/0,1,0/))
          call zsum_mn_name_xy(p%jj,idiag_jzmxy,(/0,0,1/))
        endif
      endif
    
    endsubroutine calc_2d_diagnostics_magnetic
!***********************************************************************
    subroutine time_integrals_magnetic(f,p)
!
!  Calculate time_integrals within each pencil (as long as each
!  pencil case p still contains the current data). This routine
!  is now being called at the end of equ.
!
!  28-jun-07/axel+mreinhard: coded
!  24-jun-08/axel: moved call to this routine to the individual pde routines
!   1-jul-08/axel: moved this part to magnetic
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f
      intent(in) :: p
!
      if (ibbt/=0) f(l1:l2,m,n,ibxt:ibzt)=f(l1:l2,m,n,ibxt:ibzt)+dt*p%bb
      if (ijjt/=0) f(l1:l2,m,n,ijxt:ijzt)=f(l1:l2,m,n,ijxt:ijzt)+dt*p%jj
!
    endsubroutine time_integrals_magnetic
!***********************************************************************
    subroutine df_diagnos_magnetic(df,p)
!
!  calculate diagnostics that involves df
!  Here we calculate <du/dt x b> and <u x db/dt>.
!  The latter is calculated as <divu dai/dt> -  <uji daj/dt>
!  This is used in dynamo theory for checking the minimal tau approximation.
!
!  10-oct-06/axel: coded
!
      use Diagnostics, only: sum_mn_name
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: uudot,aadot,udotxb,B1_gradu
      real, dimension (nx) :: B1dot_udotxb,B1dot_uxbdot,B1dot_aadot,uxbdot2
!
      intent(in)  :: df, p
!
!  this routine is only called when ldiagnos=T
!  start with <du/dt x b>
!
      if (idiag_udotxbm/=0) then
        uudot=df(l1:l2,m,n,iux:iuz)
        call cross_mn(uudot,p%bb,udotxb)
        call dot_mn_sv(B1_ext,udotxb,B1dot_udotxb)
        call sum_mn_name(B1dot_udotxb,idiag_udotxbm)
      endif
!
!  next, do <divu dai/dt> -  <uji daj/dt>
!
      if (idiag_uxbdotm/=0) then
        aadot=df(l1:l2,m,n,iax:iaz)
        call dot_mn_sv(B1_ext,aadot,B1dot_aadot)
        call dot_mn_sm(B1_ext,p%uij,B1_gradu)
        call dot_mn(B1_gradu,aadot,uxbdot2)
        B1dot_uxbdot=p%divu*B1dot_aadot-uxbdot2
        call sum_mn_name(B1dot_uxbdot,idiag_uxbdotm)
      endif
!
    endsubroutine df_diagnos_magnetic
!***********************************************************************
    subroutine magnetic_after_boundary(f)
!
!  Calculate <A>, which is needed for test-field methods.
!
!   2-jan-10/axel: adapted from hydro_after_boundary
!  10-jan-13/MR: added possibility to remove evolving mean field
!  15-oct-15/MR: changes for slope-limited diffusion
!   7-jun-16/MR: modifications in z average removal for Yin-Yang, yet incomplete
!  03-apr-20/joern: restructured and fixed slope-limited diffusion, now in daa_dt
!
      use Deriv, only: der_z,der2_z
      use Sub, only: finalize_aver, div, calc_all_diff_fluxes, dot2_mn
      use Yinyang_mpi, only: zsum_yy
      use Boundcond, only: update_ghosts
      use Diagnostics, only: save_name
      use Density, only: calc_pencils_density
      use Mpicomm, only: mpiallreduce_sum

      real, dimension(mx,my,mz,mfarray), intent(inout) :: f

      real :: fact
      integer :: l,n,j,ml,nl
      real, dimension(nz,3) :: gaamz,d2aamz
      real, dimension(:,:,:), allocatable :: buffer
      real, dimension(nx) :: tmp
      real, dimension(mx,3) :: aamx
      real, dimension(mx,mz) :: aamxz
      type(pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc=.false.
!
!  Compute mean field (xy verage) for each component. Include the ghost zones,
!  because they have just been set.
!
      if (lrmv) then
        if (lremove_meanax) then

          fact=1./nyzgrid
          do j=1,3
            do l=1,mx
              aamx(l,j)=fact*sum(f(l,m1:m2,n1:n2,iax+j-1))
            enddo
          enddo
          call finalize_aver(nprocyz,23,aamx)
!
          do j=1,3
            do l=1,mx
              f(l,:,:,iax+j-1) = f(l,:,:,iax+j-1)-aamx(l,j)
            enddo
          enddo

        endif
      endif
!
      if (lcalc_aameanz .or. lrmv.and.lremove_meanaz) then
!
        fact=1./nxygrid
        do j=1,3
          do n=1,mz
            aamz(n,j)=fact*sum(f(l1:l2,m1:m2,n,iax+j-1))
          enddo
        enddo
        call finalize_aver(nprocxy,12,aamz)
!
        if (lrmv.and.lremove_meanaz) then
          do j=1,3
            do n=1,mz
              f(:,:,n,iax+j-1) = f(:,:,n,iax+j-1)-aamz(n,j)
            enddo
          enddo
        endif

        if (lcalc_aameanz) then
!
!  Compute first and second derivatives.
!
          do j=1,3
            call der_z(aamz(:,j),gaamz(:,j))
            call der2_z(aamz(:,j),d2aamz(:,j))
          enddo
!
!  Compute mean magnetic field and current density.
!
          bbmz(:,1)=-gaamz(:,2)
          bbmz(:,2)=+gaamz(:,1)
          bbmz(:,3)=0.
          jjmz(:,1)=-d2aamz(:,1)
          jjmz(:,2)=-d2aamz(:,2)
          jjmz(:,3)=0.
        endif
      endif
!
!  Calculate Arms for quenching.
!
      if (lquench_eta_aniso) then
        Arms=sum(f(l1:l2,m1:m2,n1:n2,iaa:iaa+2)**2)  ! requires equidistant grid
        call mpiallreduce_sum(Arms,fact)
        Arms=sqrt(fact/nwgrid)
      endif
!
!  Remove mean field (y average).
!
      if (lrmv) then

        if (lremove_meanaxz) then
!
          fact=1./nygrid
          do j=1,3

            aamxz=fact*sum(f(:,m1:m2,:,iaa+j-1),2)  ! requires equidistant grid
            call finalize_aver(nprocy,2,aamxz)
!
            do m=1,my
              f(:,m,:,iaa+j-1) = f(:,m,:,iaa+j-1)-aamxz
            enddo

          enddo
        endif
!
!  Remove mean field (z average).
!
        if (lremove_meanaxy) then
!
          fact=1./nzgrid_eff
          if (lyang) allocate(buffer(1,mx,my))

          do j=1,3

            if (lyang) then
!
!  On Yang grid:
!
              do nl=n1,n2
                do ml=1,my
                  call zsum_yy(buffer,1,ml,nl,f(:,ml,nl,iaa+j-1))
                enddo
              enddo
              aamxy=fact*buffer(1,:,:)
            else
!
! Normal summing-up in Yin procs.
! 
              aamxy=fact*sum(f(:,:,n1:n2,iaa+j-1),3)  ! requires equidistant grid
            endif

            call finalize_aver(nprocz,3,aamxy)
!
            do n=1,mz
              f(:,:,n,iaa+j-1) = f(:,:,n,iaa+j-1)-aamxy
            enddo
          enddo

        endif
      endif
!
!  put u.a into auxiliary array
!
      if (lua_as_aux) then
        do n=1,mz
          do m=1,my
            f(:,m,n,iua)=sum(f(:,m,n,iux:iuz)*f(:,m,n,iax:iaz),2)
          enddo
        enddo
      endif
!
!  Compute eta_smag and put into auxilliary variable
!
      if (letasmag_as_aux) then
        if (ljj_as_aux) then
          do n=n1,n2; do m=m1,m2
!
            call dot2_mn(f(l1:l2,m,n,ijx:ijz),tmp)
            f(l1:l2,m,n,ietasmag)=(D_smag*dxmax)**2.*sqrt(tmp)
!
          enddo; enddo
!
          call update_ghosts(f,ietasmag)
!
        endif
      endif
!
!  The following allows us to let eta change with time, t-eta_tdep_toffset.
!  The eta_tdep_toffset is used in cosmology where time starts at t=1.
!  lresi_eta_tdep_t0_norm is not the default because of backward compatbility.
!  The default is problematic because then eta_tdep /= eta for t < eta_tdep_t0.
!
      if (lresi_eta_tdep) then
        if (lresi_eta_tdep_t0_norm) then
          eta_tdep=eta*max(real(t-eta_tdep_toffset)/eta_tdep_t0,1.)**eta_tdep_exponent
        else
          eta_tdep=eta*max(real(t-eta_tdep_toffset),eta_tdep_t0)**eta_tdep_exponent
        endif
        if (lroot.and.ldiagnos) call save_name(eta_tdep,idiag_eta_tdep)
      endif
!
!  Output kx_aa for calculating k_effective.
!
      if (lroot.and.ldiagnos) call save_name(kx_aa(1),idiag_kx_aa)
!
!     if (lmagn_mf) call magnetic_after_boundary
!
    endsubroutine magnetic_after_boundary
!***********************************************************************
    subroutine set_border_magnetic(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the aa variable.
!
!  28-jul-06/wlad: coded
!
      use BorderProfiles, only: border_driving,set_border_initcond
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: f_target
      integer :: ju,j
!
!  select for different target profiles
!
      do j=1,3
!
        select case (borderaa(j))
!
        case ('zero','0')
          f_target(:,j)=0.
!
        case ('initial-condition')
          ju=j+iaa-1
          call set_border_initcond(f,ju,f_target(:,j))
!
        case ('nothing')
          if (lroot.and.ip<=5) &
               print*,"set_border_magnetic: borderaa='nothing'"
!
        case default
          write(unit=errormsg,fmt=*) &
               'set_border_magnetic: No such value for borderaa: ', &
               trim(borderaa(j))
          call fatal_error('set_border_magnetic',errormsg)
!
        endselect
!
!  apply border profile
!
        if (borderaa(j) /= 'nothing') then
          ju=j+iaa-1
          call border_driving(f,df,p,f_target(:,j),ju)
        endif
!
      enddo
!
    endsubroutine set_border_magnetic
!***********************************************************************
    subroutine eta_shell(p,eta_mn,geta)
!
!  24-nov-03/dave: coded
!  23-jun-09/axel: generalized to lcylinder_in_a_box
!
      use Sub, only: step, der_step
      use Mpicomm, only: stop_it
!
      type (pencil_case) :: p
      real, dimension (nx) :: eta_mn
      real, dimension (nx) :: prof,eta_r
      real, dimension (nx,3) :: geta
      real :: d_int,d_ext
!
      eta_r=0.
!
      if (eta_int > 0.) then
        d_int = eta_int - eta
      else
        d_int = 0.
      endif
      if (eta_ext > 0.) then
        d_ext = eta_ext - eta
      else
        d_ext = 0.
      endif
!
!  calculate steps in resistivity
!  make this dependent on the geometry used
!
!  (i) lcylinder_in_a_box
!
      if (lcylinder_in_a_box.or.lcylindrical_coords) then
        prof=step(p%rcyl_mn,r_int,wresistivity)
        eta_mn=d_int*(1-prof)
        prof=step(p%rcyl_mn,r_ext,wresistivity)
        eta_mn=eta+eta_mn+d_ext*prof
!
!     calculate radial derivative of steps and gradient of eta
!
        prof=der_step(p%rcyl_mn,r_int,wresistivity)
        eta_r=-d_int*prof
        prof=der_step(p%rcyl_mn,r_ext,wresistivity)
        eta_r=eta_r+d_ext*prof
        geta=p%evr*spread(eta_r,2,3)
!
!  (ii) lsphere_in_a_box
!
      elseif (lsphere_in_a_box.or.lspherical_coords) then
        prof=step(p%r_mn,r_int,wresistivity)
        eta_mn=d_int*(1-prof)
        prof=step(p%r_mn,r_ext,wresistivity)
        eta_mn=eta+eta_mn+d_ext*prof
!
!     calculate radial derivative of steps and gradient of eta
!
        prof=der_step(p%r_mn,r_int,wresistivity)
        eta_r=-d_int*prof
        prof=der_step(p%r_mn,r_ext,wresistivity)
        eta_r=eta_r+d_ext*prof
        geta=p%evr*spread(eta_r,2,3)
!
!  (iii) other cases are not implemented yet
!
      else
        call stop_it("eta_shell works only for spheres or cylinders")
      endif
!
    endsubroutine eta_shell
!***********************************************************************
    subroutine calc_bthresh
!
!  calculate bthresh from brms, give warnings if there are problems
!
!   6-aug-03/axel: coded
!
!  give warning if brms is not set in prints.in
!
      if (idiag_brms==0) then
        if (lroot.and.lfirstpoint) then
          print*,'calc_bthresh: need to set brms in print.in to get bthresh'
        endif
      endif
!
!  if nvec exceeds nbvecmax (=1/4) of points per processor, then begin to
!  increase scaling factor on bthresh. These settings will stay in place
!  until the next restart
!
      if (nbvec>nbvecmax.and.lfirstpoint) then
        print*,'calc_bthresh: processor ',iproc_world,': bthresh_scl,nbvec,nbvecmax=', &
                                          bthresh_scl,nbvec,nbvecmax
        bthresh_scl=bthresh_scl*1.2
      endif
!
!  calculate bthresh as a certain fraction of brms
!  MR: can hardly be correct as at this moment the updated brms is not yet
!      available.
!
      bthresh=bthresh_scl*bthresh_per_brms*brms
!
    endsubroutine calc_bthresh
!***********************************************************************
    subroutine rescaling_magnetic(f)
!
!  Rescale magnetic field by factor rescale_aa,
!
!  22-feb-05/axel: coded
!  10-feb-09/petri: adapted from testfield
!
      use Sub, only: update_snaptime, read_snaptime
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=fnlen) :: file
      logical :: lmagnetic_out
      logical, save :: lfirst_call=.true.
!
      intent(inout) :: f
!
!  Reinitialize aa periodically if requested
!
      if (lreset_aa) then
        file=trim(datadir)//'/treset_aa.dat'
        if (lfirst_call) then
          call read_snaptime(trim(file),taareset,naareset,daareset,t)
          if (taareset==0 .or. taareset < t-daareset) then
            taareset=t+daareset
          endif
          lfirst_call=.false.
        endif
!
!  Rescale when the time has come
!  (Note that lmagnetic_out and ch are not used here)
!
        if (t >= taareset) then
          f(:,:,:,iax:iaz)=rescale_aa*f(:,:,:,iax:iaz)
          call update_snaptime(file,taareset,naareset,daareset,t,lmagnetic_out)
        endif
      endif
!
    endsubroutine rescaling_magnetic
!***********************************************************************
    subroutine calc_tau_aa_exterior(f,df)
!
!  magnetic field relaxation to zero on time scale tau_aa_exterior within
!  exterior region. For the time being this means z > zgrav.
!
!  29-jul-02/axel: coded
!
      use Gravity, only: zgrav
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: scl
      integer :: j
!
      intent(in) :: f
      intent(inout) :: df
!
      if (headtt) print*,'calc_tau_aa_exterior: tau=',tau_aa_exterior
      if (z(n)>zgrav) then
        scl=1./tau_aa_exterior
        do j=iax,iaz
          df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-scl*f(l1:l2,m,n,j)
        enddo
      endif
!
    endsubroutine calc_tau_aa_exterior
!***********************************************************************
    subroutine helflux(aa,uxb,jj)
!
!  magnetic helicity flux (preliminary)
!
!  14-aug-03/axel: coded
!
      use Diagnostics
!
      real, dimension (nx,3), intent(in) :: aa,uxb,jj
      real, dimension (nx,3) :: ee
      real, dimension (nx) :: FHx,FHz
      real :: FH
!
      ee=eta*jj-uxb
!
!  calculate magnetic helicity flux in the X and Z directions
!
      FHx=-2*ee(:,3)*aa(:,2)*dsurfyz
      FHz=+2*ee(:,1)*aa(:,2)*dsurfxy
!
!  sum up contribution per pencil
!  and then stuff result into surf_mn_name for summing up all processors.
!
      FH=FHx(nx)-FHx(1)
      if (lfirst_proc_z.and.n==n1) FH=FH-sum(FHz)
      if (llast_proc_z .and.n==n2) FH=FH+sum(FHz)
      call surf_mn_name(FH,idiag_exaym2)
!
    endsubroutine helflux
!***********************************************************************
    subroutine curflux_dS(uxb,jj)
!
!  current helicity flux (preliminary)
!
!  27-nov-03/axel: adapted from helflux
!
      use Diagnostics
!
      real, dimension (nx,3), intent(in) :: uxb,jj
      real, dimension (nx,3) :: ee
      real, dimension (nx) :: FCx,FCz
      real :: FC
!
      ee=eta*jj-uxb
!
!  calculate current helicity flux in the X and Z directions
!  Could speed up by only calculating here boundary points!
!
      FCx=2*(ee(:,2)*jj(:,3)-ee(:,3)*jj(:,2))*dsurfyz
      FCz=2*(ee(:,1)*jj(:,2)-ee(:,2)*jj(:,1))*dsurfxy
!
!  sum up contribution per pencil
!  and then stuff result into surf_mn_name for summing up all processors.
!
      FC=FCx(nx)-FCx(1)
      if (lfirst_proc_z.and.n==n1) FC=FC-sum(FCz)
      if (llast_proc_z .and.n==n2) FC=FC+sum(FCz)
      call surf_mn_name(FC,idiag_exjm2)
!
    endsubroutine curflux_dS
!***********************************************************************
    subroutine curflux(uxb,jj)
!
!  current helicity flux (preliminary)
!
!  27-nov-03/axel: adapted from helflux
!
      use Diagnostics
!
      real, dimension (nx,3), intent(in) :: uxb,jj
      real, dimension (nx,3) :: ee
      real, dimension (nx) :: FCz
!
      ee=eta*jj-uxb
!
!  calculate current helicity flux in the Z direction
!  exj = e1*j2 - e2*j1
!
      FCz=2*(ee(:,1)*jj(:,2)-ee(:,2)*jj(:,1))
      call sum_mn_name(FCz,idiag_exjm2)
!
    endsubroutine curflux
!***********************************************************************
    subroutine read_magnetic_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=magnetic_init_pars, IOSTAT=iostat)
!
!  read namelist for mean-field theory (if invoked)
!
      if (lmagn_mf) call read_magn_mf_init_pars(iostat)
!
    endsubroutine read_magnetic_init_pars
!***********************************************************************
    subroutine write_magnetic_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=magnetic_init_pars)
!
!  write namelist for mean-field theory (if invoked)
!
      if (lmagn_mf) call write_magn_mf_init_pars(unit)
!
    endsubroutine write_magnetic_init_pars
!***********************************************************************
    subroutine read_magnetic_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=magnetic_run_pars, IOSTAT=iostat)
!
!  read namelist for mean-field theory (if invoked)
!
      if (lmagn_mf) call read_magn_mf_run_pars(iostat)
!
    endsubroutine read_magnetic_run_pars
!***********************************************************************
    subroutine write_magnetic_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=magnetic_run_pars)
!
!  write namelist for mean-field theory (if invoked)
!
      if (lmagn_mf) call write_magn_mf_run_pars(unit)
!
    endsubroutine write_magnetic_run_pars
!***********************************************************************
    subroutine forcing_continuous(df,p)
!
!  add a continuous forcing term (here currently only for localized rotors)
!
!  21-jan-07/axel: adapted from hydro
!  24-feb-09/axel: calls to this routine are now replaced by adding p$fcont
!
      use Diagnostics
      use Mpicomm, only: mpibcast_real
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: forcing_rhs
      real, dimension (nx) :: jf,phi
      real, dimension (mx), save :: phix,sinx,cosx
      real, dimension (my), save :: phiy,siny,cosy
      real, dimension (mz), save :: phiz,sinz,cosz
      real, save :: phase_beltrami_before
      real, save :: R2,R12
      type (pencil_case) :: p
      logical, save :: lfirst_call=.true.
      integer :: j,jfff,ifff
      real :: fact
!
!  at the first step, the sin and cos functions are calculated for all
!  x,y,z points and are then saved and used for all subsequent steps
!  and pencils
!
      if (lfirst_call) then
        if (ip<=6) print*,'forcing_continuous: lfirst_call=',lfirst_call
        if (iforcing_continuous_aa=='fixed_swirl') then
          if (lroot) print*,'forcing_continuous: fixed_swirl; swirl=',swirl
          R2=radius**2
          R12=1./R2
          phix=exp(-R12*x**2)
          phiy=exp(-R12*y**2)
          phiz=exp(-R12*z**2)
        elseif (iforcing_continuous_aa=='cosxcosz') then
          cosx=cos(k1x_ff*x)
          cosz=cos(k1z_ff*z)
        elseif (iforcing_continuous_aa=='Azsinx') then
          sinx=cos(k1z_ff*x)
        elseif (iforcing_continuous_aa=='Aycosz') then
          cosz=cos(k1z_ff*z)
        elseif (iforcing_continuous_aa=='RobertsFlow') then
          if (lroot) print*,'forcing_continuous: RobertsFlow'
          sinx=sin(k1_ff*x); cosx=cos(k1_ff*x)
          siny=sin(k1_ff*y); cosy=cos(k1_ff*y)
        elseif (iforcing_continuous_aa=='Beltrami-z') then
          if (lroot) print*,'forcing_continuous: Beltrami-z'
          ampl_beltrami=ampl_ff
          sinz=sin(k1_ff*z+phase_beltrami)
          cosz=cos(k1_ff*z+phase_beltrami)
        endif
        lfirst_call=.false.
      endif
      if (ip<=6) print*,'forcing_continuous: dt, lfirst_call=',dt,lfirst_call
!
!  at the first meshpoint, and if phase_beltrami is finite,
!  recalculate sinz and cosz for the phase correction of an
!  imposed Beltrami field.
!
      if (lfirstpoint) then
        if (iforcing_continuous_aa=='Beltrami-z') then
          if (phase_beltrami/=phase_beltrami_before) then
            phase_beltrami_before=phase_beltrami
            sinz=sin(k1_ff*z+phase_beltrami)
            cosz=cos(k1_ff*z+phase_beltrami)
          endif
        endif
      endif
!
!  calculate forcing
!
      if (iforcing_continuous_aa=='fixed_swirl') then
        fact=ampl_ff
        phi=2.*R12*fact*phix(l1:l2)*phiy(m)*phiz(n)
        forcing_rhs(:,1)=(-swirl*y(m    )+2.*x(l1:l2)*z(n))*phi
        forcing_rhs(:,2)=(+swirl*x(l1:l2)+2.*y(m    )*z(n))*phi
        forcing_rhs(:,3)=(R2-x(l1:l2)**2-y(m)**2)*2.*R12*phi
      elseif (iforcing_continuous_aa=='cosxcosz') then
        fact=ampl_ff
        forcing_rhs(:,1)=0.
        forcing_rhs(:,2)=fact*cosx(l1:l2)*cosz(n)
        forcing_rhs(:,3)=0.
      elseif (iforcing_continuous_aa=='Azsinx') then
        fact=ampl_ff
        forcing_rhs(:,1)=0.
        forcing_rhs(:,2)=0.
        forcing_rhs(:,3)=fact*sinx(l1:l2)
      elseif (iforcing_continuous_aa=='Aycosz') then
        fact=ampl_ff
        forcing_rhs(:,1)=0.
        forcing_rhs(:,2)=fact*cosz(n)
        forcing_rhs(:,3)=0.
      elseif (iforcing_continuous_aa=='RobertsFlow') then
        fact=ampl_ff
        forcing_rhs(:,1)=-fact*cosx(l1:l2)*siny(m)
        forcing_rhs(:,2)=+fact*sinx(l1:l2)*cosy(m)
        forcing_rhs(:,3)=+fact*cosx(l1:l2)*cosy(m)*sqrt(2.)
      elseif (iforcing_continuous_aa=='Beltrami-z') then
        fact=-eta*k1_ff*ampl_beltrami
        forcing_rhs(:,1)=fact*cosz(n)
        forcing_rhs(:,2)=fact*sinz(n)
        forcing_rhs(:,3)=0.
      endif
!
!  apply forcing in uncurled induction equation
!
      ifff=iax
      do j=1,3
        jfff=j+ifff-1
        df(l1:l2,m,n,jfff)=df(l1:l2,m,n,jfff)+forcing_rhs(:,j)
      enddo
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_jfm/=0) then
          call dot_mn(p%jj,forcing_rhs,jf)
          call sum_mn_name(jf,idiag_jfm)
        endif
      endif
!
    endsubroutine forcing_continuous
!***********************************************************************
    subroutine get_slices_magnetic(f,slices)
!
!  Write slices for animation of Magnetic variables.
!
!  26-jul-06/tony: coded
!  14-apr-16/MR: changes for Yin-Yang
!
      use General, only: transform_thph_yy_other
      use Slices_methods, only: assign_slices_vec, assign_slices_scal

      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Magnetic vector potential (code variable)
!
        case ('aa')
          call assign_slices_vec(slices,f,iaa)
!
!  Phi component of magnetic vector potential times axis distance (derived variable)
!
        case ('aps')
          call assign_slices_scal(slices,aps_xy,aps_xz,aps_yz,xz2=aps_xz2)
!
!  Magnetic field (derived variable)
!
        case ('bb')
          call assign_slices_vec(slices,bb_xy,bb_xz,bb_yz,bb_xy2,bb_xy3,bb_xy4,bb_xz2)
!
!  Current density (derived variable)
!
        case ('jj')
          call assign_slices_vec(slices,jj_xy,jj_xz,jj_yz,jj_xy2,jj_xy3,jj_xy4,jj_xz2)
!
!  Magnetic field squared (derived variable)
!
        case ('b2')
          call assign_slices_scal(slices,b2_xy,b2_xz,b2_yz,b2_xy2,b2_xy3,b2_xy4,b2_xz2)
!
!  Current squared (derived variable)
!
        case ('j2')
          call assign_slices_scal(slices,j2_xy,j2_xz,j2_yz,j2_xy2,j2_xy3,j2_xy4,j2_xz2)
!
!  Magnetic field in spherical coordinates (derived variable)
!
        case ('bb_sph')
          call assign_slices_vec(slices,bb_sph_xy,bb_sph_xz,bb_sph_yz,bb_sph_xy2,bb_sph_xy3,bb_sph_xy4,bb_sph_xz2)
!
!  Current density times magnetic field (derived variable)
!
        case ('jb')
          call assign_slices_scal(slices,jb_xy,jb_xz,jb_yz,jb_xy2,jb_xy3,jb_xy4,jb_xz2)
!
!  Plasma beta
!
       case ('beta1')
          call assign_slices_scal(slices,beta1_xy,beta1_xz,beta1_yz,beta1_xy2,beta1_xy3,&
                                  beta1_xy4,beta1_xz2)
!
       case ('beta')
         if (lroot) then
           print*,"The 'beta' slice was renamed to 'beta1'. Please "
           print*,"update the label in your video.in file."
         endif
         call fatal_error("get_slices_magnetic","")
!
! Poynting vector
!
        case ('poynting')
          call assign_slices_vec(slices,poynting_xy,poynting_xz,poynting_yz,poynting_xy2,&
                                 poynting_xy3,poynting_xy4,poynting_xz2)
!
!  Magnetic helicity density
!
        case ('ab')
          call assign_slices_scal(slices,ab_xy,ab_xz,ab_yz,ab_xy2,ab_xy3,ab_xy4,ab_xz2)
!
      endselect
!
    endsubroutine get_slices_magnetic
!***********************************************************************
    subroutine calc_mfield
!
!  calculate mean magnetic field from xy- or z-averages
!
!  19-jun-02/axel: moved from print to here
!   9-nov-02/axel: corrected bxmy(m,j); it used bzmy instead!
!
      use Mpicomm
      use Sub
!
!  For vector output (of bb vectors) we need brms
!  on all processors. It suffices to have this for times when lout=.true.,
!  but we need to broadcast the result to all procs.
!
!  calculate brms (this requires that brms is set in print.in)
!  broadcast result to other processors
!
      if (idiag_brms/=0) then
        if (iproc==0) brms=fname(idiag_brms)
        call mpibcast_real(brms)
      endif
!
!  The following calculation involving spatial averages
!
      if (idiag_bmx/=0) call calc_bmx
      if (idiag_bmy/=0) call calc_bmy
      if (idiag_bmz/=0) call calc_bmz
      if (idiag_bmzS2/=0) call calc_bmzS2
      if (idiag_bmzA2/=0) call calc_bmzA2
      if (idiag_jmx/=0) call calc_jmx
      if (idiag_jmy/=0) call calc_jmy
      if (idiag_jmz/=0) call calc_jmz
      if (idiag_emxamz3/=0) call calc_emxamz3
      if (idiag_embmz/=0) call calc_embmz
      if (idiag_ambmz/=0) call calc_ambmz
      if (idiag_ambmzh/=0) call calc_ambmzh
      if (idiag_jmbmz/=0.or.idiag_kmz/=0) call calc_jmbmz
      if (idiag_bmxy_rms/=0) call calc_bmxy_rms
      if (idiag_bmzph/=0) call calc_bmz_beltrami_phase
!
!  Set the phase of the Beltrami forcing equal to the actual phase
!  of the magnetic field (times forcing_continuous_aa_phasefact).
!
      phase_beltrami=forcing_continuous_aa_phasefact*bmz_beltrami_phase
!
!  set amplitude to ampl_ff minus a correction term that is
!  proportional to the actual field minus the target field strength,
!  scaled by some forcing_continuous_aa_amplfact, and broadcast, ie
!  A = Atarget - factor*(Aactual-Atarget).
!
      ampl_beltrami=ampl_ff-forcing_continuous_aa_amplfact*(bmz-ampl_ff)
      call mpibcast_real(ampl_beltrami,comm=MPI_COMM_WORLD)
!
    endsubroutine calc_mfield
!***********************************************************************
    subroutine calc_bmx
!
!  Magnetic energy in the yz-averaged field. The bymxy and bzmxy must have
!  been calculated, so they are present on the z-root processors.
!
!   6-apr-08/axel: moved from calc_mfield to here
!  26-aug-09/anders: used mpireduce_sum to remove need for nyproc arrays
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real, dimension(nx) :: bymx, bzmx, bmx2
      real, dimension(nx,ny) :: fsumxy
      real :: bmx
!
!  This only works if bymxy and bzmxy are in zaver, so print warning if this is
!  not ok.
!
      if (idiag_bymxy==0.or.idiag_bzmxy==0) then
        if (first) then
          print*, 'calc_mfield: WARNING'
          print*, 'NOTE: to get bmx, set bymxy and bzmxy in zaver'
          print*, 'We proceed, but you will get bmx=0'
        endif
        bmx2=0.0
      else
        if (lfirst_proc_z) then
          call mpireduce_sum(fnamexy(idiag_bymxy,:,:),fsumxy,(/nx,ny/),idir=2)
          bymx=sum(fsumxy,dim=2)/nygrid
          call mpireduce_sum(fnamexy(idiag_bzmxy,:,:),fsumxy,(/nx,ny/),idir=2)
          bzmx=sum(fsumxy,dim=2)/nygrid
        endif
        if (lfirst_proc_yz) then
          call mpireduce_sum(bymx**2+bzmx**2,bmx2,nx,idir=1)
        endif
      endif
!
!  Save the name in the idiag_bmx slot and set first to false.
!  Compute final result only on the root processor.
!  This is not current working with x-parallelization!
!
      if (lroot) then
        bmx=sqrt(sum(bmx2)/nxgrid)
        call save_name(bmx,idiag_bmx)
      endif
      first=.false.
!
    endsubroutine calc_bmx
!***********************************************************************
    subroutine calc_bmy
!
!  Magnetic energy in the xz-averaged field. The bxmxy and bzmxy must have
!  been calculated, so they are present on the z-root processors.
!
!   6-apr-08/axel: moved from calc_mfield to here
!  26-aug-09/axel: adapted change of Anders to used mpireduce_sum
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real, dimension(ny) :: bxmy, bzmy, bmy2
      real, dimension(nx,ny) :: fsumxy
      real :: bmy
!
!  This only works if bxmxy and bzmxy are in zaver, so print warning if this is
!  not ok.
!
      if (idiag_bxmxy==0.or.idiag_bzmxy==0) then
        if (first) then
          print*, 'calc_mfield: WARNING'
          print*, 'NOTE: to get bmy, set bxmxy and bzmxy in zaver'
          print*, 'We proceed, but you will get bmy=0'
        endif
        bmy2=0.0
      else
        if (lfirst_proc_z) then
          call mpireduce_sum(fnamexy(idiag_bxmxy,:,:),fsumxy,(/nx,ny/),idir=1)
          bxmy=sum(fsumxy,dim=1)/nxgrid
          call mpireduce_sum(fnamexy(idiag_bzmxy,:,:),fsumxy,(/nx,ny/),idir=1)
          bzmy=sum(fsumxy,dim=1)/nxgrid
        endif
        if (lfirst_proc_xz) then
          call mpireduce_sum(bxmy**2+bzmy**2,bmy2,ny,idir=2)
        endif
      endif
!
!  Save the name in the idiag_bmy slot and set first to false.
!  Compute final result only on the root processor.
!
      if (lroot) then
        bmy=sqrt(sum(bmy2)/nygrid)
        call save_name(bmy,idiag_bmy)
      endif
      first=.false.
!
    endsubroutine calc_bmy
!***********************************************************************
    subroutine calc_bmzS2
!
!  Magnetic energy in anisymmetric part of the horizontally averaged field.
!  The bxmz and bymz must have been calculated, and present on root processor.
!
!   8-mar-10/axel: adapted from bmz
!
      use Diagnostics, only: save_name
!
      logical,save :: first=.true.
      integer :: n_reverse,ipz_reverse
      real :: bxmS,bymS
!
!  This only works if bxmz and bzmz are in xyaver, so print warning if this is
!  not ok.
!
      if (idiag_bxmz==0.or.idiag_bymz==0) then
        if (first) then
          print*, 'calc_mfield: WARNING'
          print*, 'NOTE: to get bmzS2, set bxmz and bymz in xyaver'
          print*, 'We proceed, but you will get bmzS2=0'
        endif
        bxmS=0.
        bymS=0.
      else
        bxmS=0.
        bymS=0.
        do n=1,nz
          n_reverse=nz-n+1
          ipz_reverse=ipz !(for now)
          bxmS=bxmS+fnamez(n,ipz+1,idiag_bxmz)+fnamez(n_reverse,ipz_reverse+1,idiag_bxmz)
          bymS=bymS+fnamez(n,ipz+1,idiag_bymz)+fnamez(n_reverse,ipz_reverse+1,idiag_bymz)
        enddo
      endif
!
!  Save the name in the idiag_bmzS slot and set first to false.
!
      call save_name((bxmS**2+bymS**2)/nz**2,idiag_bmzS2)
      first=.false.
!
    endsubroutine calc_bmzS2
!***********************************************************************
    subroutine calc_bmzA2
!
!  Magnetic energy in anisymmetric part of the horizontally averaged field.
!  The bxmz and bymz must have been calculated, and present on root processor.
!
!   8-mar-10/axel: adapted from bmz
!
      use Diagnostics, only: save_name
!
      logical,save :: first=.true.
      integer :: n_reverse,ipz_reverse
      real :: bxmA,bymA
!
!  This only works if bxmz and bzmz are in xyaver, so print warning if this is
!  not ok.
!
      if (idiag_bxmz==0.or.idiag_bymz==0) then
        if (first) then
          print*, 'calc_mfield: WARNING'
          print*, 'NOTE: to get bmzA2, set bxmz and bymz in xyaver'
          print*, 'We proceed, but you will get bmzA2=0'
        endif
        bxmA=0.
        bymA=0.
      else
        bxmA=0.
        bymA=0.
        do n=1,nz
          n_reverse=nz-n+1
          ipz_reverse=ipz !(for now)
          bxmA=bxmA+fnamez(n,ipz+1,idiag_bxmz)-fnamez(n_reverse,ipz_reverse+1,idiag_bxmz)
          bymA=bymA+fnamez(n,ipz+1,idiag_bymz)-fnamez(n_reverse,ipz_reverse+1,idiag_bymz)
        enddo
      endif
!
!  Save the name in the idiag_bmzA slot and set first to false.
!
      call save_name((bxmA**2+bymA**2)/nz**2,idiag_bmzA2)
      first=.false.
!
    endsubroutine calc_bmzA2
!***********************************************************************
    subroutine calc_bmz
!
!  Magnetic energy in horizontally averaged field. The bxmz and bymz must have
!  been calculated, so they are present on the root processor.
!
!   6-apr-08/axel: moved from calc_mfield to here
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
!
!  This only works if bxmz and bzmz are in xyaver, so print warning if this is
!  not ok.
!
      if (idiag_bxmz==0.or.idiag_bymz==0) then
        if (first) then
          print*, 'calc_mfield: WARNING'
          print*, 'NOTE: to get bmz, set bxmz and bymz in xyaver'
          print*, 'We proceed, but you will get bmz=0'
        endif
        bmz=0.0
      else
        bmz=sqrt(sum(fnamez(:,:,idiag_bxmz)**2 &
                    +fnamez(:,:,idiag_bymz)**2)/(nz*nprocz))
      endif
!
!  Save the name in the idiag_bmz slot and set first to false.
!
      call save_name(bmz,idiag_bmz)
      first=.false.
!
    endsubroutine calc_bmz
!***********************************************************************
    subroutine calc_jmx
!
!  Magnetic energy in the yz-averaged field. The jymxy and jzmxy must have
!  been calculated, so they are present on the z-root processors.
!
!   6-apr-08/axel: moved from calc_mfield to here
!  28-feb-10/axel: final jmx can only be computed on root processor
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real, dimension(nx) :: jymx,jzmx,jmx2
      real, dimension(nx,ny) :: fsumxy
      real :: jmx
!
!  This only works if jymxy and jzmxy are in zaver, so print warning if this is
!  not ok.
!
      if (idiag_jymxy==0.or.idiag_jzmxy==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get jmx, set jymxy and jzmxy in zaver"
          print*,"We proceed, jut you'll get jmx=0"
        endif
        jmx2=0.
      else
        if (lfirst_proc_z) then
          call mpireduce_sum(fnamexy(idiag_jymxy,:,:),fsumxy,(/nx,ny/),idir=2)
          jymx=sum(fsumxy,dim=2)/nygrid
          call mpireduce_sum(fnamexy(idiag_jzmxy,:,:),fsumxy,(/nx,ny/),idir=2)
          jzmx=sum(fsumxy,dim=2)/nygrid
        endif
        if (lfirst_proc_yz) then
          call mpireduce_sum(jymx**2+jzmx**2,jmx2,nx,idir=1)
        endif
      endif
!
!  Save the name in the idiag_jmx slot and set first to false.
!  Compute final result only on the root processor.
!
      if (lroot) then
        jmx=sqrt(sum(jmx2)/nxgrid)
        call save_name(jmx,idiag_jmx)
      endif
      first=.false.
!
    endsubroutine calc_jmx
!***********************************************************************
    subroutine calc_jmy
!
!  Magnetic energy in the xz-averaged field. The jxmxy and jzmxy must have
!  been calculated, so they are present on the z-root processors.
!
!   6-apr-08/axel: moved from calc_mfield to here
!  28-feb-10/axel: final jmy can only be computed on root processor
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real, dimension(ny) :: jxmy,jzmy,jmy2
      real, dimension(nx,ny) :: fsumxy
      real :: jmy
!
!  This only works if jxmxy and jzmxy are in zaver, so print warning if this is
!  not ok.
!
      if (idiag_jxmxy==0.or.idiag_jzmxy==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get jmy, set jxmxy and jzmxy in zaver"
          print*,"We proceed, but you'll get jmy=0"
        endif
        jmy2=0.
      else
        if (lfirst_proc_z) then
          call mpireduce_sum(fnamexy(idiag_jxmxy,:,:),fsumxy,(/nx,ny/),idir=1)
          jxmy=sum(fsumxy,dim=1)/nxgrid
          call mpireduce_sum(fnamexy(idiag_jzmxy,:,:),fsumxy,(/nx,ny/),idir=1)
          jzmy=sum(fsumxy,dim=1)/nxgrid
        endif
        if (lfirst_proc_xz) then
          call mpireduce_sum(jxmy**2+jzmy**2,jmy2,ny,idir=2)
        endif
      endif
!
!  Save the name in the idiag_jmy slot and set first to false.
!  Compute final result only on the root processor.
!
      if (lroot) then
        jmy=sqrt(sum(jmy2)/nygrid)
        call save_name(jmy,idiag_jmy)
      endif
      first=.false.
!
    endsubroutine calc_jmy
!***********************************************************************
    subroutine calc_jmz
!
!  Magnetic energy in horizontally averaged field. The jxmz and jymz must have
!  been calculated, so they are present on the root processor.
!
!   6-apr-08/axel: moved from calc_mfield to here
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real :: jmz
!
!  This only works if jxmz and jzmz are in xyaver, so print warning if this is
!  not ok.
!
      if (idiag_jxmz==0.or.idiag_jymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get jmz, set jxmz and jymz in xyaver"
          print*,"We proceed, but you'll get jmz=0"
        endif
        jmz=0.
      else
        jmz=sqrt(sum(fnamez(:,:,idiag_jxmz)**2 &
                    +fnamez(:,:,idiag_jymz)**2)/(nz*nprocz))
      endif
!
!  Save the name in the idiag_jmz slot and set first to false.
!
      if (lroot) call save_name(jmz,idiag_jmz)
      first=.false.
!
    endsubroutine calc_jmz
!***********************************************************************
    subroutine calc_embmz
!
!  Magnetic helicity production of mean field. The bxmz and bymz as well as Exmz
!  and Eymz must have been calculated, so they are present on the root
!  processor.
!
!   6-apr-08/axel: moved from calc_mfield to here
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real :: embmz
!
!  This only works if bxmz and bzmz are in xyaver, so print warning if this is
!  not ok.
!
      if (idiag_Exmz==0.or.idiag_Eymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get embmz, set bxmz, bymz, Exmz, and Eymz in xyaver"
          print*,"We proceed, but you'll get embmz=0"
        endif
        embmz=0.
      else
        embmz=sum(fnamez(:,:,idiag_bxmz)*fnamez(:,:,idiag_Exmz) &
                 +fnamez(:,:,idiag_bymz)*fnamez(:,:,idiag_Eymz))/(nz*nprocz)
      endif
!
!  Save the name in the idiag_embmz slot and set first to false.
!
      if (lroot) call save_name(embmz,idiag_embmz)
      first=.false.
!
    endsubroutine calc_embmz
!***********************************************************************
    subroutine calc_emxamz3
!
!  Volume average of magnetic helicity flux of the mean field. The axmz and
!  aymz as well as Exmz and Eymz must have been calculated, so they are present
!  on the root processor.
!
!   6-apr-08/axel: moved from calc_mfield to here
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real :: emxamz3
!
!  This only works if bxmz and bzmz are in xyaver, so print warning if this is
!  not ok.
!
      if (idiag_Exmz==0.or.idiag_Eymz==0.or. &
          idiag_axmz==0.or.idiag_aymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get emxamz3, set axmz, aymz, Exmz, and Eymz in xyaver"
          print*,"We proceed, but you'll get emxamz3=0"
        endif
        emxamz3=0.
      else
        emxamz3=sum(fnamez(:,:,idiag_Exmz)*fnamez(:,:,idiag_aymz) &
                   -fnamez(:,:,idiag_Eymz)*fnamez(:,:,idiag_axmz))/(nz*nprocz)
      endif
!
!  Save the name in the idiag_emxamz3 slot and set first to false.
!
      if (lroot) call save_name(emxamz3,idiag_emxamz3)
      first=.false.
!
    endsubroutine calc_emxamz3
!***********************************************************************
    subroutine calc_ambmz
!
!  Magnetic helicity of the xy-averaged mean field. The bxmz and bymz as well as
!  axmz and aymz must have been calculated, so they are present on the root
!  processor.
!
!  16-may-09/axel: adapted from calc_jmbmz
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real :: ambmz
!
!  This only works if bxmz and bzmz are in xyaver, so print warning if this is
!  not ok.
!
      if (idiag_axmz==0.or.idiag_aymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get ambmz, set bxmz, bymz, axmz, and aymz in xyaver"
          print*,"We proceed, but you'll get ambmz=0"
        endif
        ambmz=0.
      else
        ambmz=sum(fnamez(:,:,idiag_bxmz)*fnamez(:,:,idiag_axmz) &
                 +fnamez(:,:,idiag_bymz)*fnamez(:,:,idiag_aymz))/(nz*nprocz)
      endif
!
!  Save the name in the idiag_ambmz slot and set first to false.
!
      if (lroot) call save_name(ambmz,idiag_ambmz)
      first=.false.
!
    endsubroutine calc_ambmz
!***********************************************************************
    subroutine calc_ambmzh
!
!  Hemispheric magnetic helicity of the xy-averaged mean field. The bxmz and
!  bymz as well as axmz and aymz must have been calculated, so they are present
!  on the root processor.
!
!  16-may-09/axel: adapted from calc_ambmz
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real, dimension(2) :: ambmzh
      real :: ambmz_tmp,fact
      integer :: iprocz
!
!  initialize ambmzh to zero each time, because its two elements
!  are integration counters. If idiag_axmz etc are not set, the
!  routine escapes, but even then it needs to be initialized.
!
      ambmzh=0.
!
!  This only works if bxmz and bzmz are in xyaver,
!  so print warning if this is not ok.
!  Loop over all processors, but don't use (overwrite) ipz for that
!
      if (idiag_axmz==0.or.idiag_aymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get ambmzh, set bxmz, bymz, axmz, and aymz in xyaver"
          print*,"We proceed, but you'll get ambmzh=0"
        endif
      else
        fact=1./(nz*nprocz)
        do n=1,nz
        do iprocz=1,nprocz
          ambmz_tmp=fact*( &
             fnamez(n,iprocz,idiag_bxmz)*fnamez(n,iprocz,idiag_axmz)&
            +fnamez(n,iprocz,idiag_bymz)*fnamez(n,iprocz,idiag_aymz))
          if (z_allprocs(n,iprocz)>=zequator) then
            ambmzh(2)=ambmzh(2)+ambmz_tmp
          else
            ambmzh(1)=ambmzh(1)+ambmz_tmp
          endif
        enddo
        enddo
      endif
!
!  North means 1 and south means 2.
!  save the name in the idiag_ambmz slot
!  and set first to false
!
      if (lroot) call save_name_halfz(ambmzh,idiag_ambmzh)
      fname(idiag_ambmzn)=fname_half(idiag_ambmzh,1)
      fname(idiag_ambmzs)=fname_half(idiag_ambmzh,2)
      itype_name(idiag_ambmzn)=ilabel_save
      itype_name(idiag_ambmzs)=ilabel_save
      first=.false.
!
    endsubroutine calc_ambmzh
!***********************************************************************
    subroutine calc_jmbmz
!
!  Current helicity of the xy-averaged mean field
!  The bxmz and bymz as well as jxmz and jymz must have been calculated,
!  so they are present on the root processor.
!
!  21-apr-08/axel: adapted from calc_embmz
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real :: jmbmz
!
!  This only works if bxmz and bzmz are in xyaver, so print warning if this is
!  not ok.
!
      if (idiag_jxmz==0.or.idiag_jymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get jmbmz, set bxmz, bymz, jxmz, and jymz in xyaver"
          print*,"We proceed, but you'll get jmbmz=0"
        endif
        jmbmz=0.
      else
        jmbmz=sum(fnamez(:,:,idiag_bxmz)*fnamez(:,:,idiag_jxmz) &
                 +fnamez(:,:,idiag_bymz)*fnamez(:,:,idiag_jymz))/(nz*nprocz)
      endif
!
!  Save the name in the idiag_jmbmz slot and set first to false.
!
      if (lroot) then
        if (idiag_jmbmz/=0) call save_name(jmbmz,idiag_jmbmz)
        if (idiag_kmz/=0) call save_name(jmbmz/bmz**2,idiag_kmz)
      endif
      first=.false.
!
    endsubroutine calc_jmbmz
!***********************************************************************
    subroutine calc_bmxy_rms
!
!  Magnetic energy in z averaged field. The bxmxy, bymxy and bzmxy must have
!  been calculated, so they are present on the root processor.
!
!   6-apr-08/axel: moved from calc_mfield to here
!
      use Diagnostics
      use Mpicomm
!
      logical, save :: first = .true.
      real, dimension(2) :: b2mxy_local, b2mxy
      real :: bmxy_rms, nVol2d_local, btemp
      integer :: l
!
      if (lcylindrical_coords) &
          call stop_it("bmxy_rms not yet implemented for cylindrical")
!
      if (.not. lfirst_proc_z) return
!
!  The quantities bxmxy,bymxy,bzmxy are required in "zaver.in"
!
      bmxy_rms = 0.0
      if ((idiag_bxmxy == 0) .or. (idiag_bymxy == 0) .or. (idiag_bzmxy == 0)) then
        if (first) then
          print*,"calc_bmxy_rms: WARNING"
          print*,"NOTE: to get bmxy_rms, set bxmxy, bymxy and bzmxy in 'zaver.in'"
          print*,"We proceed, but you'll get bmxy_rms=0"
        endif
      else
        b2mxy_local = 0.0
        do l=1, nx
          do m=1, ny
            btemp = fnamexy(idiag_bxmxy,l,m)**2 + fnamexy(idiag_bymxy,l,m)**2 + fnamexy(idiag_bzmxy,l,m)**2
            if (lspherical_coords) then
               btemp = btemp*r2_weight(l)*sinth_weight(m)
               nvol2d_local = r2_weight(l)*sinth_weight(m)
            endif
            b2mxy_local(1) = b2mxy_local(1) + btemp
            b2mxy_local(2) = b2mxy_local(2) + nVol2d_local
          enddo
        enddo
        call mpireduce_sum(b2mxy_local(:),b2mxy,2,idir=2)
        if (lfirst_proc_x) then
          if (lcartesian_coords) bmxy_rms = sqrt(b2mxy(1) / (nxgrid*nygrid))
          if (lspherical_coords) bmxy_rms = sqrt(b2mxy(1) / b2mxy(2))
        endif
      endif
!
!  Save the name in the idiag_bmxy_rms slot and set first to false.
!
      call save_name(bmxy_rms,idiag_bmxy_rms)
      first = .false.
!
    endsubroutine calc_bmxy_rms
!***********************************************************************
    subroutine calc_bmz_beltrami_phase
!
!  The following is useful if the xy-averaged field is a Beltrami field
!  Determine its phase as in B ~ [ cos(kz+phi), sin(kz+phi), 0 ].
!  This means that for positive phi the wave is shifted to the left.
!
!  bxmz, bymz must have been calculated,
!  so they are present on the root processor.
!
!   2-apr-08/MR: introduced phase calculation for Beltrami mean fields
!   6-apr-08/axel: moved from calc_mfield to here
!
      use Diagnostics
      use Mpicomm
!
      logical,save :: first=.true.
      real :: bmz_belphase1,bmz_belphase2
      real, dimension (nz,nprocz), save :: sinz,cosz
      real ::  c, s
      integer :: jprocz
!
      if (first) then
        sinz=sin(k1_ff*z_allprocs); cosz=cos(k1_ff*z_allprocs)
      endif
!
!  print warning if bxmz and bymz are not calculated
!
      if (idiag_bxmz==0.or.idiag_bymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"to get bmz_beltrami_phase, set bxmz, bymz in zaver"
          print*,"We proceed, but you'll get Beltrami phase bmzpb=0"
        endif
        bmz_beltrami_phase=0.
!
!  add up c = <B_x> cos(kz) and s = <B_x> sin(kz)
!  and determine phase of Beltrami field from <B_x>
!
      else
        c=0.; s=0.
        do jprocz=1,nprocz
          c=c+dot_product(fnamez(:,jprocz,idiag_bxmz),cosz(:,jprocz))
          s=s+dot_product(fnamez(:,jprocz,idiag_bxmz),sinz(:,jprocz))
        enddo
        bmz_belphase1=atan2(-s,c)
!
!  add up c = <B_y> cos(kz) and s = <B_y> sin(kz)
!  and determine phase of Beltrami field from <B_y>
!
        c=0.; s=0.
        do jprocz=1,nprocz
          c=c+dot_product(fnamez(:,jprocz,idiag_bymz),cosz(:,jprocz))
          s=s+dot_product(fnamez(:,jprocz,idiag_bymz),sinz(:,jprocz))
        enddo
        bmz_belphase2=atan2(c,s)
!
!  Difference of both determinations (to estimate error)
!  and take the mean of both calculations (called bmz_beltrami_phase
!  and bmzph in the print.in file, for brevity)
!
        bmz_beltrami_phase=.5*(bmz_belphase1+bmz_belphase2)
      endif
!
!  Save the name in the idiag_bmzph slot; as estimate of the error,
!  calculate also the difference between the two.
!  Finally, set first to false
!
      call save_name(bmz_beltrami_phase,idiag_bmzph)
      if (idiag_bmzphe/=0) &
          call save_name(abs(bmz_belphase1-bmz_belphase2),idiag_bmzphe)
      first=.false.
!
    endsubroutine calc_bmz_beltrami_phase
!***********************************************************************
    subroutine alfven_x(ampl,f,iuu,iaa,ilnrho,kx)
!
!  Alfven wave propagating in the x-direction
!
!  uy = +sink(x-vA*t)
!  Az = -cosk(x-vA*t)*sqrt(rho*mu0)/k
!
!  Alfven and slow magnetosonic are the same here and both incompressible, and
!  a fast magnetosonic (compressible) wave is also excited, but decoupled.
!
!  satisfies the four equations
!  dlnrho/dt = -ux'
!  dux/dt = -cs2*(lnrho)'
!  duy/dt = B0*By'  ==>  duy/dt = -B0*Az''
!  dBy/dt = B0*uy'  ==>  dAz/dt = -B0*ux
!
!   8-nov-03/axel: coded
!  29-apr-03/axel: added sqrt(rho*mu0)/k factor
!   7-aug-17/axel: added sqrt(.75) for lrelativistic_eos=T
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rho,ampl_Az
      real :: ampl,kx,ampl_lr,ampl_ux,ampl_uy
      integer :: iuu,iaa,ilnrho
!
!  Amplitude factors
!
      ampl_lr=+0.
      ampl_ux=+0.
      if (ldensity.and.lrelativistic_eos) then
        ampl_uy=+ampl*sqrt(.75)
      else
        ampl_uy=+ampl
      endif
!
!  ux and Ay.
!  Don't overwrite the density, just add to the log of it.
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,ilnrho)=ampl_lr*(sin(kx*x(l1:l2))+f(l1:l2,m,n,ilnrho))
        f(l1:l2,m,n,iuu+0 )=ampl_ux*sin(kx*x(l1:l2))
        f(l1:l2,m,n,iuu+1 )=ampl_uy*sin(kx*x(l1:l2))
        rho=exp(f(l1:l2,m,n,ilnrho))
        ampl_Az=-ampl*sqrt(rho*mu0)/kx
        f(l1:l2,m,n,iaa+2 )=ampl_Az*cos(kx*x(l1:l2))
      enddo; enddo
!
    endsubroutine alfven_x
!***********************************************************************
    subroutine alfven_y(ampl,f,iuu,iaa,ky,mu0)
!
!  Alfven wave propagating in the y-direction; can be used in 2-d runs.
!  ux = cos(ky-ot), for B0y=1 and rho=1.
!  Az = sin(ky-ot), ie Bx=-cos(ky-ot)
!
!  [wd nov-2006: There should be a 1/ky in the aa term here and in
!  alfven_x, I think]
!
!  satisfies the equations
!  dux/dt = Bx'  ==>  dux/dt = -Az''
!  dBx/dt = ux'  ==>  dAz/dt = -ux.
!
!  06-dec-06/wolf: adapted from alfven_z
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,ky,mu0
      integer :: iuu,iaa
!
!  ux and Az
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+0) = +ampl*cos(ky*y(m))
        f(l1:l2,m,n,iaa+2) = -ampl*sin(ky*y(m))*sqrt(mu0)/ky
      enddo; enddo
!
    endsubroutine alfven_y
!***********************************************************************
    subroutine alfven_z(ampl,f,iuu,iaa,kz,mu0)
!
!  Alfven wave propagating in the z-direction
!  ux = cos(kz-ot), for B0z=1 and rho=1.
!  Ay = sin(kz-ot), ie Bx=-cos(kz-ot)
!
!  satisfies the equations
!  dux/dt = Bx'  ==>  dux/dt = -Ay''
!  dBx/dt = ux'  ==>  dAy/dt = -ux.
!
!  18-aug-02/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kz,mu0
      integer :: iuu,iaa
!
!  ux and Ay
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+0)=+ampl*cos(kz*z(n))
        f(l1:l2,m,n,iaa+1)=+ampl*sin(kz*z(n))*sqrt(mu0)
      enddo; enddo
!
    endsubroutine alfven_z
!***********************************************************************
    subroutine alfven_xy(ampl,f,iuu,iaa,kx,ky)
!
!  Alfven wave propagating in the xy-direction; can be used in 2-d runs.
!  uz = cos(kx*x+ky*y-ot), for B0=(1,1,0) and rho=1.
!  Ax = sin(kx*x+ky*y-ot),
!  Ay = sin(kx*x+ky*y-ot),
!
!  16-jun-07/axel: adapted from alfven_y
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,ky,om
      real, parameter :: mu0=1.
      integer :: iuu,iaa
!
!  set ux, Ax, and Ay
!
      om=B_ext(1)*kx+B_ext(2)*ky
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+2)=+ampl*cos(kx*x(l1:l2)+ky*y(m))
        f(l1:l2,m,n,iaa+0)=+ampl*sin(kx*x(l1:l2)+ky*y(m))*sqrt(mu0)/om*B_ext(2)
        f(l1:l2,m,n,iaa+1)=-ampl*sin(kx*x(l1:l2)+ky*y(m))*sqrt(mu0)/om*B_ext(1)
      enddo; enddo
!
    endsubroutine alfven_xy
!***********************************************************************
    subroutine alfven_xz(ampl,f,iuu,iaa,kx,kz)
!
!  Alfven wave propagating in the xz-direction; can be used in 2-d runs.
!  uz = cos(kx*x+kz*z-ot), for B0=(1,1,0) and rho=1.
!  Ax = sin(kx*x+kz*z-ot),
!  Az = sin(kx*x+kz*z-ot),
!
!  16-jun-07/axel: adapted from alfven_xy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,kz,om
      real, parameter :: mu0=1.
      integer :: iuu,iaa
!
!  set ux, Ax, and Az
!
      om=B_ext(1)*kx+B_ext(3)*kz
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+2)=+ampl*cos(kx*x(l1:l2)+kz*z(n))
        f(l1:l2,m,n,iaa+0)=+ampl*sin(kx*x(l1:l2)+kz*z(n))*sqrt(mu0)/om*B_ext(2)
        f(l1:l2,m,n,iaa+2)=-ampl*sin(kx*x(l1:l2)+kz*z(n))*sqrt(mu0)/om*B_ext(1)
      enddo; enddo
!
    endsubroutine alfven_xz
!***********************************************************************
    subroutine alfvenz_rot(ampl,f,iuu,iaa,kz,O)
!
!  Alfven wave propagating in the z-direction (with Coriolis force)
!  ux = cos(kz-ot), for B0z=1 and rho=1.
!  Ay = sin(kz-ot), ie Bx=-cos(kz-ot)
!
!  satisfies the equations
!  dux/dt - 2Omega*uy = -Ay''
!  duy/dt + 2Omega*ux = +Ax''
!  dAx/dt = +uy
!  dAy/dt = -ux
!
!  18-aug-02/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kz,O,fac
      integer :: iuu,iaa
!
!  ux, uy, Ax and Ay
!
      if (lroot) print*,'alfvenz_rot: Alfven wave with rotation; O,kz=',O,kz
      fac=-O+sqrt(O**2+kz**2)
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+0)=-ampl*sin(kz*z(n))*fac/kz
        f(l1:l2,m,n,iuu+1)=-ampl*cos(kz*z(n))*fac/kz
        f(l1:l2,m,n,iaa+0)=+ampl*sin(kz*z(n))/kz
        f(l1:l2,m,n,iaa+1)=+ampl*cos(kz*z(n))/kz
      enddo; enddo
!
    endsubroutine alfvenz_rot
!***********************************************************************
    subroutine alfvenz_bell(ampl,f,iuu,iaa,kz,B0,J0)
!
!  Bell instability for the pressureless equations
!  du/dt = j x B0 - J0 x b
!  da/dt = u x B0
!
!  13-jan-12/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kz,B0,J0,lam,oA
      integer :: iuu,iaa
!
!  ux, uy, Ax and Ay
!
      oA=kz*B0
      lam=sqrt(oA*(J0-oA))
      if (lroot) print*,'alfvenz_bell: Bell inst., lam,kz,oA,B0,J0=',lam,kz,oA,B0,J0
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+0)=-ampl*sin(kz*z(n))*lam/kz
        f(l1:l2,m,n,iuu+1)=+ampl*cos(kz*z(n))*lam/kz
        f(l1:l2,m,n,iaa+0)=+ampl*cos(kz*z(n))*B0/kz
        f(l1:l2,m,n,iaa+1)=+ampl*sin(kz*z(n))*B0/kz
      enddo; enddo
!
    endsubroutine alfvenz_bell
!***********************************************************************
    subroutine alfvenz_rot_shear(ampl,f,iuu,iaa,kz,OO)
!
!  Alfven wave propagating in the z-direction (with Coriolis force and shear)
!
!  satisfies the equations
!  dux/dt - 2*Omega*uy = -Ay''
!  duy/dt + (2-q)*Omega*ux = +Ax''
!  dAx/dt = q*Omega*Ay + uy
!  dAy/dt = -ux
!
!  Assume B0=rho0=mu0=1
!
!  28-june-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kz,OO
      complex :: fac
      integer :: iuu,iaa
!
!  ux, uy, Ax and Ay
!
      if (lroot) print*,'alfvenz_rot_shear: '// &
          'Alfven wave with rotation and shear; OO,kz=',OO,kz
      fac=cmplx(OO-sqrt(16*kz**2+OO**2),0.)
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+0)=f(l1:l2,m,n,iuu+0)+ampl*fac/(4*kz)*sin(kz*z(n))
        f(l1:l2,m,n,iuu+1)=f(l1:l2,m,n,iuu+1)+ampl*real(exp(cmplx(0,z(n)*kz))* &
            fac*sqrt(2*kz**2+OO*fac)/(sqrt(2.)*kz*(-6*OO-fac)))
        f(l1:l2,m,n,iaa+0)=f(l1:l2,m,n,iaa+0)+ampl*sin(kz*z(n))/kz
        f(l1:l2,m,n,iaa+1)=f(l1:l2,m,n,iaa+1)-ampl*2*sqrt(2.)*aimag(exp(cmplx(0,z(n)*kz))* &
            sqrt(2*kz**2+OO*fac)/(-6*OO-fac)/(cmplx(0,kz)))
      enddo; enddo
!
    endsubroutine alfvenz_rot_shear
!***********************************************************************
    subroutine torus_test(ampl,f)
!
!  Initial field concentrated along torus well inside the computational
!  domain.
!  Implements the same field for cartesian and spherical cordinates.
!  The field is of mixed parity (bb_pol symmetric, bb_tor antisymmetric)
!  and the relative contributions of toroidal and poloidal field are
!  determined by
!    ampl(1) -- bb_pol (through aa_tor)
!    ampl(3) -- bb_tor (through aa_pol)
!  Uses x_max as reference radius.
!
!   05-may-2008/wolf: coded
!
      real :: ampl
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx) :: xxi2,ee
      real, dimension (nx) :: costh,sinth,cosphi,sinphi,ss,rr,aar,aap
      real :: radius,width,r_cent
!
      intent(in)    :: ampl
      intent(inout) :: f
!
      radius = xyz1(1)
      width  = 0.1*radius
      r_cent = 0.6*radius
!
      if (lspherical_coords) then
        do n=n1,n2; do m=m1,m2
          xxi2 = (x(l1:l2)*sin(y(m))-r_cent)**2 + x(l1:l2)**2*cos(y(m))**2
          ee = ampl * exp(-0.5 * xxi2 / width**2)
          f(l1:l2,m,n,iax) = f(l1:l2,m,n,iax) + ee * x(l1:l2)*cos(y(m))
          f(l1:l2,m,n,iaz) = f(l1:l2,m,n,iaz) + ee
        enddo; enddo
      else
        do n=n1,n2; do m=m1,m2
          xxi2 = (sqrt(x(l1:l2)**2+y(m)**2) - r_cent)**2 + z(n)**2
          ee = ampl * exp(-0.5 * xxi2 / width**2)
          aar = z(n) * ee
          aap = ee
          ss = sqrt(x(l1:l2)**2+y(m)**2)
          rr = sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          ss = max(ss, tini)
          rr = max(rr, tini)
          costh = z(n)/rr
          sinth = ss/rr
          cosphi   = x(l1:l2)/ss
          sinphi   = y(m)/ss
          f(l1:l2,m,n,iax) = f(l1:l2,m,n,iax) + aar*sinth*cosphi - aap*sinphi
          f(l1:l2,m,n,iay) = f(l1:l2,m,n,iay) + aar*sinth*sinphi + aap*cosphi
          f(l1:l2,m,n,iaz) = f(l1:l2,m,n,iaz) + aar*costh
        enddo; enddo
      endif
!
    endsubroutine torus_test
!***********************************************************************
    subroutine force_free_jet(mu)
!
!  Force free magnetic field configuration for jet simulations
!  with a fixed accretion disc at the bottom boundary.
!
!  The input parameter mu specifies the radial dependency of
!  the magnetic field in the disc.
!
!  Solves the laplace equation in cylindrical coordinates for the
!  phi-component of the vector potential. A_r and A_z are taken to
!  be zero.
!
!    nabla**2 A_phi - A_phi / r**2 = 0
!
!  For the desired boundary condition in the accretion disc
!
!    B_r=B0*r**(mu-1)  (z == 0)
!
!  the solution is
!
!    A_phi = Hypergeometric2F1( (1-mu)/2, (2+mu)/2, 2, xi**2 )
!            *xi*(r**2+z**2)**(mu/2)
!
!  where xi = sqrt(r**2/(r**2+z**2))
!
!  30-may-04/tobi: coded
!
      use Sub, only: hypergeometric2F1,gamma_function
      use Deriv, only: der
      use Debug_IO, only: output
!
      real, intent(in) :: mu
      real :: xi2,A_phi
      real :: r2
      real :: B1r_,B1z_,B1
      real, parameter :: tol=10*epsilon(1.0)
      integer :: l
      real, dimension(mx,my,mz) :: Ax_ext,Ay_ext
      !real, dimension(nx,3) :: bb_ext_pot
      !real, dimension(nx) :: bb_x,bb_y,bb_z
!
!  calculate un-normalized |B| at r=r_ref and z=0 for later normalization
!
      if (lroot.and.ip<=5) print*,'FORCE_FREE_JET: calculating normalization'
!
      B1r_=sin(pi*mu/2)*gamma_function(   abs(mu) /2) / &
                        gamma_function((1+abs(mu))/2)
!
      B1z_=cos(pi*mu/2)*gamma_function((1+abs(mu))/2) / &
                        gamma_function((2+abs(mu))/2)
!
      B1=sqrt(4/pi)*r_ref**(mu-1)*sqrt(B1r_**2+B1z_**2)
!
!  calculate external vector potential
!
      if (lroot) print*,'FORCE_FREE_JET: calculating external vector potential'
!
      if (lforce_free_test) then
!
        if (lroot) print*,'FORCE_FREE_JET: using analytic solution for mu=-1'
        do l=1,mx; do m=1,my; do n=1,mz
          Ax_ext=-2*y(m)*(1-z(n)/sqrt(x(l)**2+y(m)**2+z(n)**2))/(x(l)**2+y(m)**2)/B1
          Ay_ext= 2*x(l)*(1-z(n)/sqrt(x(l)**2+y(m)**2+z(n)**2))/(x(l)**2+y(m)**2)/B1
        enddo; enddo; enddo
!
      else
!
        do l=1,mx; do m=1,my; do n=1,mz
          r2=x(l)**2+y(m)**2
          xi2=r2/(r2+z(n)**2)
          A_phi=hypergeometric2F1((1-mu)/2,(2+mu)/2,2.0,xi2,tol) &
               *sqrt(xi2)*sqrt(r2+z(n)**2)**mu/B1
!
          Ax_ext(l,m,n)=-y(m)*A_phi/sqrt(r2)
          Ay_ext(l,m,n)= x(l)*A_phi/sqrt(r2)
        enddo; enddo; enddo
!
      endif
!
!  calculate external magnetic field
!
      if (lroot.and.ip<=5) &
        print*,'FORCE_FREE_JET: calculating the external magnetic field'
!
      do n=n1,n2
      do m=m1,m2
!        call der(Ay_ext,bb_x,3)
!        bb_ext_pot(:,1)=-bb_x
!        call der(Ax_ext,bb_y,3)
!        bb_ext_pot(:,2)= bb_y
!        call der(Ay_ext,bb_z,1)
!        bb_ext_pot(:,3)= bb_z
!        call der(Ax_ext,bb_z,2)
!        bb_ext_pot(:,3)=bb_ext_pot(:,3)-bb_z
!        call set_global(bb_ext_pot,m,n,'B_ext_pot',nx)
      enddo
      enddo
!
      if (ip<=5) then
        call output(trim(directory)//'/Ax_ext.dat',Ax_ext,1)
        call output(trim(directory)//'/Ay_ext.dat',Ay_ext,1)
      endif
!
    endsubroutine force_free_jet
!***********************************************************************
    subroutine piecew_dipole_aa(ampl,inclaa,f,ivar)
!
!  A field that is vertical uniform for r<R_int, inclined dipolar for
!  r>R_ext, and potential in the shell R_int<r<R_ext.
!  This mimics a neutron star just after the Meissner effect forced the
!  internal field to become vertical (aligned with rotation axis).
!
!  AMPL represents mu/4 pi, where  mu = 1/2 Int rr jj dV  is the
!  magnetic moment of the external dipole field.
!  INCLAA is the inclination of the dipolar field.
!
!  Pencilized in order to minimize memory consumption with all the
!  auxiliary variables used.
!
!  23-jul-05/wolf:coded
!
      real, intent(inout), dimension (mx,my,mz,mfarray) :: f
      real, intent(in) :: ampl,inclaa
      real, dimension (nx) :: r_1_mn,r_2_mn,sigma0,sigma1, r_mn
      real :: fact
      real, dimension(2) :: beta(0:1)
      real, dimension(2,3) :: a(0:1,1:3),b(0:1,1:3)
      integer :: ivar
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
        r_1_mn = 1./max(r_mn,tini)
        r_2_mn = 1./max(r_mn**2,tini)
!
        fact = ampl
        ! beta = [beta_1^0, beta_1^1] combines coefficients for m=0, m=1
        beta =  fact * (/ cos(inclaa), -sin(inclaa)/sqrt(2.) /)
        ! a and b for m=0, m=1 (index 1) and interior, shell, exterior index 2)
        a(0,:) = (/ 1./r_ext**3, 1./r_ext**3                  , 0. /) * beta(0)
        a(1,:) = (/ 0.         , 1./(r_ext**3-r_int**3)       , 0. /) * beta(1)
        !
        b(0,:) = (/ 0.         , 0.                           , 1. /) * beta(0)
        b(1,:) = (/ 0.         , -r_int**3/(r_ext**3-r_int**3), 1. /) * beta(1)
        !
        ! The following could be coded much clearer using elsewhere, but
        ! that is not in F90 (and pgf90 doesn't support F95)
        ! r_in < r < r_ext
        sigma0 = a(0,2)*r_mn + b(0,2)*r_2_mn
        sigma1 = a(1,2)*r_mn + b(1,2)*r_2_mn
        where(r_mn>r_ext) ! r > r_ext
          sigma0 = a(0,3)*r_mn + b(0,3)*r_2_mn
          sigma1 = a(1,3)*r_mn + b(1,3)*r_2_mn
        endwhere
        where(r_mn<r_int) ! r < r_int
          sigma0 = a(0,1)*r_mn + b(0,1)*r_2_mn
          sigma1 = a(1,1)*r_mn + b(1,1)*r_2_mn
        endwhere
        sigma1 = sigma1*sqrt(2.)
        f(l1:l2,m,n,ivar+0) = -sigma0*y(m)*r_1_mn
        f(l1:l2,m,n,ivar+1) =  sigma0*x(l1:l2)*r_1_mn + sigma1*z(n)*r_1_mn
        f(l1:l2,m,n,ivar+2) =                     - sigma1*y(m)*r_1_mn
      enddo
!
    endsubroutine piecew_dipole_aa
!***********************************************************************
    subroutine geo_benchmark_B(f)
!
!  30-june-04/grs: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx) :: theta_mn,ar,atheta,aphi,r_mn,phi_mn
      real :: C_int,C_ext,A_int,A_ext
      integer :: j
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
        theta_mn=acos(spread(z(n),1,nx)/r_mn)
        phi_mn=atan2(spread(y(m),1,nx),x(l1:l2))
!
! calculate ax,ay,az (via ar,atheta,aphi) inside shell (& leave zero outside shell)
!
        do j=1,ninit
           select case (initaa(j))
           case ('geo-benchmark-case1')
              if (lroot .and. imn==1) print*, 'geo_benchmark_B: geo-benchmark-case1'
              C_int=-( -1./63.*r_int**4 + 11./84.*r_int**3*r_ext            &
                     + 317./1050.*r_int**2*r_ext**2                         &
                     - 1./5.*r_int**2*r_ext**2*log(r_int) )
              C_ext=-( -1./63.*r_ext**9 + 11./84.*r_ext**8*r_int            &
                     + 317./1050.*r_ext**7*r_int**2                         &
                     - 1./5.*r_ext**7*r_int**2*log(r_ext) )
              A_int=5./2.*(r_ext-r_int)
              A_ext=5./8.*(r_ext**4-r_int**4)
!
              where (r_mn < r_int)
                ar=C_int*ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)*r_mn
                atheta=3.*C_int*ampl_B0*80.*sin(2.*theta_mn)*r_mn
                aphi=ampl_B0*A_int*r_mn*sin(theta_mn)
              endwhere
!
              where (r_mn <= r_ext .and. r_mn >= r_int)
                ar=ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)*                 &
                   (   1./36.*r_mn**5 - 1./12.*(r_int+r_ext)*r_mn**4        &
                     + 1./14.*(r_int**2+4.*r_int*r_ext+r_ext**2)*r_mn**3    &
                     - 1./3.*(r_int**2*r_ext+r_int*r_ext**2)*r_mn**2        &
                     - 1./25.*r_int**2*r_ext**2*r_mn                        &
                     + 1./5.*r_int**2*r_ext**2*r_mn*log(r_mn) )
                atheta=-ampl_B0*80.*sin(2.*theta_mn)*                       &
                   (   7./36.*r_mn**5 - 1./2.*(r_int+r_ext)*r_mn**4         &
                     + 5./14.*(r_int**2+4.*r_int*r_ext+r_ext**2)*r_mn**3    &
                     - 4./3.*(r_int**2*r_ext+r_int*r_ext**2)*r_mn**2        &
                     + 2./25.*r_int**2*r_ext**2*r_mn                        &
                     + 3./5.*r_int**2*r_ext**2*r_mn*log(r_mn) )
                aphi=ampl_B0*5./8.*sin(theta_mn)*                           &
                   ( 4.*r_ext*r_mn - 3.*r_mn**2 - r_int**4/r_mn**2 )
              endwhere
!
              where (r_mn > r_ext)
                ar=C_ext*ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)/r_mn**4
                atheta=-2.*C_ext*ampl_B0*80.*sin(2.*theta_mn)/r_mn**4
                aphi=ampl_B0*A_ext/r_mn**2*sin(theta_mn)
              endwhere
!
          ! debug checks -- look at a pencil near the centre...
              if (ip<=4 .and. imn==(ny+1)*nz/2) then
                 print*,'r_int,r_ext',r_int,r_ext
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(r_mn), imn, iproc:', &
                      iproc, imn, minval(r_mn), maxval(r_mn)
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(theta_mn), imn, iproc:', &
                      iproc, imn, minval(theta_mn), maxval(theta_mn)
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(phi_mn), imn, iproc:', &
                      iproc, imn, minval(phi_mn), maxval(phi_mn)
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(ar), imn, iproc:', &
                      iproc, imn, minval(ar), maxval(ar)
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(atheta), imn, iproc:', &
                      iproc, imn, minval(atheta), maxval(atheta)
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(aphi), imn, iproc:', &
                      iproc, imn, minval(aphi), maxval(aphi)
              endif
!
           case ('geo-benchmark-case2')
              if (lroot .and. imn==1) print*, 'geo_benchmark_B: geo-benchmark-case2 not yet coded.'
!
           case ('nothing')
              if (lroot .and. imn==1) print*, 'geo_benchmark_B: nothing more to add (in ninit loop).'
              exit
!
           case default
              if (lroot .and. imn==1) print*,'geo_benchmark_B: case not defined!'
              call stop_it("")
           endselect
        enddo
        f(l1:l2,m,n,iax)=sin(theta_mn)*cos(phi_mn)*ar + cos(theta_mn)*cos(phi_mn)*atheta - sin(phi_mn)*aphi
        f(l1:l2,m,n,iay)=sin(theta_mn)*sin(phi_mn)*ar + cos(theta_mn)*sin(phi_mn)*atheta + cos(phi_mn)*aphi
        f(l1:l2,m,n,iaz)=cos(theta_mn)*ar - sin(theta_mn)*atheta
     enddo
!
!
     if (ip<=14) then
        print*,'geo_benchmark_B: minmax(ax) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iax)),maxval(f(l1:l2,m1:m2,n1:n2,iax))
        print*,'geo_benchmark_B: minmax(ay) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iay)),maxval(f(l1:l2,m1:m2,n1:n2,iay))
        print*,'geo_benchmark_B: minmax(az) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iaz)),maxval(f(l1:l2,m1:m2,n1:n2,iaz))
     endif
!
    endsubroutine geo_benchmark_B
!***********************************************************************
    subroutine eta_xy_dep(eta_xy,geta_xy,eta_xy_profile)
!
!   2-jul-2009/koen: creates an xy-dependent resistivity (for RFP studies)
!   (under reconstruction)
!
      real, dimension(mx,my) :: eta_xy,r2,gradr_eta_xy
      real, dimension(mx,my,3)  :: geta_xy
      character (len=labellen) :: eta_xy_profile
      real :: rmax2,a,w
      integer :: i,j
!
      intent(out) :: eta_xy,geta_xy
!
      select case (eta_xy_profile)
      case ('schnack89')
      do i=1,mx
        do j=1,my
          r2(i,j)=x(i)**2+y(j)**2
        enddo
      enddo
!
!  define eta_xy: radial resistivity profile from Y.L. Ho, S.C. Prager &
!              D.D. Schnack, Phys rev letters vol 62 nr 13 1989
!  and define gradr_eta_xy: 1/r *d_r(eta_xy))
!
!  to prevent numerically impossible diffusivities the value outside rmax2
!  the diffusivity is set to stay below an input value eta_xy_max (> 100*eta),
!  keeping the transition continuous in the first derivative
!
      rmax2=1.
      a=(eta_xy_max-100.*eta)/eta_xy_max
      w=(eta_xy_max-100.*eta)/(5400.*eta)
!
      do i=1,mx
      do j=1,my
!  inside
        if (r2(i,j) < rmax2) then
          eta_xy(i,j) = eta*(1.+9.*(r2(i,j)/rmax2)**15)**2
          gradr_eta_xy= 540.*eta*(1.+9.*(r2(i,j)/rmax2)**15)*(r2(i,j)/rmax2)**14/rmax2**0.5
!  outside
        else
          eta_xy(i,j) = eta_xy_max*(1.-a*exp(-(r2(i,j)**0.5-rmax2**0.5)/w))
          gradr_eta_xy(i,j)= eta_xy_max*(a*exp(-(r2(i,j)**0.5-rmax2**0.5)/w))/w/rmax2**0.5
        endif
!  gradient
          geta_xy(i,j,1) = x(i)*gradr_eta_xy(i,j)
          geta_xy(i,j,2) = y(j)*gradr_eta_xy(i,j)
          geta_xy(i,j,3) = 0.
      enddo
      enddo
!
      endselect
!
    endsubroutine eta_xy_dep
!***********************************************************************
    subroutine eta_zdep(zdep_profile, nz, z, eta_z, geta_z)
!
!  creates a z-dependent resistivity for protoplanetary disk studies
!
!  12-jul-2005/joishi: coded
!
!  Input Arguments
!    zdep_profile: name/code of the z-dependent profile
!    nz: number of elements in array z
!    z: z coordinates
!  Output Arguments
!    eta_z: resistivity at the corresponding z
!  Optional Output Arguments
!    geta_z: gradient of resistivity at the corresponding z
!
      use General, only: erfcc
      use Sub, only: step, der_step, cubic_step, cubic_der_step,erfunc
      use EquationOfState, only: cs0
!
      character(len=labellen), intent(in) :: zdep_profile
      integer, intent(in) :: nz
      real, dimension(nz), intent(in) :: z
      real, dimension(nz), intent(out) :: eta_z
      real, dimension(nz), intent(out), optional :: geta_z
!
      real, dimension(nz) :: zoh, z2
      real :: h
!
      select case (zdep_profile)
        case ('fs')
          if (cs0 > 0. .and. Omega > 0.) then
            h = sqrt(2.) * cs0 / Omega
          else
            h = 1.
          endif
          zoh = z / h
          z2 = zoh**2
!  resistivity profile from Fleming & Stone (ApJ 585:908-920)
!          eta_z = eta*exp(-z2/2.+sigma_ratio*erfcc(abs(z))/4.)
          eta_z = eta * exp(-0.5 * z2 + 0.25 * sigma_ratio * erfcc(abs(zoh)))
!
! its gradient:
          if (present(geta_z)) then
            geta_z = -eta_z * (zoh + (0.5 * sigma_ratio / sqrt(pi)) * sign(1.,z) * exp(-z2)) / h
          endif
!
        case ('tanh')
!  default to spread gradient over ~5 grid cells.
          if (eta_zwidth == 0.) eta_zwidth = 5.*dz
          eta_z = eta*0.5*(tanh((z + eta_z0)/eta_zwidth) &
                - tanh((z - eta_z0)/eta_zwidth))
!
! its gradient:
          if (present(geta_z)) then
             geta_z = -eta/(2.*eta_zwidth) * &
               ((tanh((z + eta_z0)/eta_zwidth))**2. &
               -(tanh((z - eta_z0)/eta_zwidth))**2.)
          endif
!
!  Single tanh step function
!
        case ('step')
!
!  default to spread gradient over ~5 grid cells.
!
          if (eta_zwidth == 0.) eta_zwidth = 5.*dz
          eta_z = eta + eta*(eta_jump-1.)*step(z,eta_z0,-eta_zwidth)
!
! its gradient:
          if (present(geta_z)) then
            geta_z = eta*(eta_jump-1.)*der_step(z,eta_z0,-eta_zwidth)
          endif
!
!
!  Cubic-step profile
!
        case ('cubic_step')
!
          if (eta_zwidth == 0.) eta_zwidth = 5.*dz
          eta_z = eta + eta*(eta_jump-1.)*cubic_step(z,eta_z0,-eta_zwidth)
!
! its gradient:
          if (present(geta_z)) then
            geta_z = eta*(eta_jump-1.)*cubic_der_step(z,eta_z0,-eta_zwidth)
          endif
!
! zlayer: eta is peaked within eta_z0 and eta_z1 and decreases outside. 
! the profile is made exactly the same as the similarly named profile in 
! the forcing module.
!
        case ('zlayer')
!
!  Default to spread gradient over ~5 grid cells,
!
          if (eta_zwidth == 0.) eta_zwidth = 5.*dz
          eta_z=eta*eta_jump + (eta-eta*eta_jump)*(.5*(1.+erfunc((z-eta_z0)/eta_zwidth)) &
                  -.5*(1.+erfunc((z-eta_z1)/eta_zwidth)) )
!
! and its gradient
!
          if (present(geta_z)) then
            geta_z=(eta-eta*eta_jump)*(exp(((z-eta_z0)**2)/(eta_zwidth**2)) &
                                - exp(((z-eta_z1)**2)/(eta_zwidth**2)))/(eta_zwidth*sqrt(pi))
          endif
!
!
!  Two-step function
!
        case ('two-step','two_step')
!
!  Default to spread gradient over ~5 grid cells,
!
          if (eta_zwidth == 0.) eta_zwidth = 5.*dz
          eta_z = eta*eta_jump-eta*(eta_jump-two_step_factor)* &
            (step(z,eta_z0,eta_zwidth)-step(z,eta_z1,eta_zwidth))
!
!  ... and its gradient. Note that the sign of the second term enters
!  with the opposite sign, because we have used negative eta_zwidth.
!
          if (present(geta_z)) then
            geta_z = eta*(eta_jump-two_step_factor)*( &
              der_step(z,eta_z0,-eta_zwidth)+der_step(z,eta_z1,eta_zwidth))
          endif
!
!
!  Cubic-two-step profile
!
        case ('cubic_two_step')
!
          if (eta_zwidth == 0.) eta_zwidth = 5.*dz
          if (eta_zwidth2 == 0.) eta_zwidth2 = 5.*dz
          eta_z = (eta + eta*(eta_jump-1.)*cubic_step(z,eta_z0,-eta_zwidth))* &
                  (1+(eta_jump2-1.)*cubic_step(z,eta_z1,-eta_zwidth2))
!
! its gradient:
          if (present(geta_z)) then
            geta_z = (eta*(eta_jump-1.)*cubic_der_step(z,eta_z0,-eta_zwidth))* &
                  (1+(eta_jump2-1.)*cubic_step(z,eta_z1,-eta_zwidth2)) &
                + (eta + eta*(eta_jump-1.)*cubic_step(z,eta_z0,-eta_zwidth))* &
                  ((eta_jump2-1.)*cubic_der_step(z,eta_z1,-eta_zwidth2))
          endif
!
!  Powerlaw-z
!
        case ('powerlaw-z','powerlaw_z')
!
!  eta proportional to chosen power of z...
!
          eta_z = eta*(1.+(z-z0)/eta_z0)**eta_power_z
!
!  ... and its gradient.
!
          if (present(geta_z)) geta_z = eta_power_z*eta_z/(eta_z0+z)
      endselect
!
    endsubroutine eta_zdep
!***********************************************************************
    subroutine eta_ydep(ydep_profile, ny, y, eta_y, geta_y)
!
!  creates a z-dependent resistivity for protoplanetary disk studies
!
!  19-aug-2013/wlad: adapted from eta_zdep
!
      use General, only: erfcc
      use Sub, only: step, der_step
!
      character(len=labellen), intent(in) :: ydep_profile
      integer, intent(in) :: ny
      real, dimension(ny), intent(in) :: y
      real, dimension(ny), intent(out) :: eta_y
      real, dimension(ny), intent(out), optional :: geta_y
!
      select case (ydep_profile)
!
!  Two-step function
!
        case ('two-step','two_step')
!
!  Default to spread gradient over ~5 grid cells,
!
          if (eta_ywidth == 0.) eta_ywidth = 5.*dy
          eta_y = eta*eta_jump-eta*(eta_jump-two_step_factor)* &
            (step(y,eta_y0,eta_ywidth)-step(y,eta_y1,eta_ywidth))
!
!  ... and its gradient. Note that the sign of the second term enters
!  with the opposite sign, because we have used negative eta_ywidth.
!
          if (present(geta_y)) then
            geta_y = eta*(eta_jump-two_step_factor)*( &
              der_step(y,eta_y0,-eta_ywidth)+der_step(y,eta_y1,eta_ywidth))
          endif
!
      endselect
!
    endsubroutine eta_ydep
!***********************************************************************
    subroutine eta_xdep(eta_x,geta_x,xdep_profile)
!
!  creates a x-dependent resistivity for protoplanetary disk studies
!
!  09-mar-2011/wlad: adapted from eta_zdep
!  10-mar-2011/axel: corrected gradient term: should point in x
!
      use General, only: erfcc
      use Sub, only: step, der_step
!
      real, dimension(mx) :: eta_x,x2
      real, dimension(mx) :: geta_x
      character (len=labellen) :: xdep_profile
      integer :: l
!
      intent(out) :: eta_x,geta_x
!
      select case (xdep_profile)
        case ('fs')
          x2 = x**2.
!
!  resistivity profile from Fleming & Stone (ApJ 585:908-920)
!
          eta_x = eta*exp(-x2/2.+sigma_ratio*erfcc(abs(x))/4.)
!
! its gradient:
!
          geta_x = eta_x*(-x-sign(1.,x)*sigma_ratio*exp(-x2)/(2.*sqrt(pi)))
!
!  tanh profile
!
        case ('tanh')
!
!  default to spread gradient over ~5 grid cells.
!
           if (eta_xwidth == 0.) eta_xwidth = 5.*dx
           eta_x = eta*0.5*(tanh((x + eta_x0)/eta_xwidth) &
             - tanh((x - eta_x0)/eta_xwidth))
!
! its gradient:
!
           geta_x = -eta/(2.*eta_xwidth) * ((tanh((x + eta_x0)/eta_xwidth))**2. &
             - (tanh((x - eta_x0)/eta_xwidth))**2.)
!
!  linear profile
!
        case ('linear')
!
!  default to spread gradient over ~5 grid cells.
!
           if (eta_xwidth == 0.) eta_xwidth = Lxyz(1)
           eta_x = eta*(1.+(x-xyz1(1))/eta_xwidth)
!
! its gradient:
!
           geta_x = eta/eta_xwidth
!
        case ('quadratic')

           eta_x = eta*((2.*x-(xyz0(1)+xyz1(1)))/Lxyz(1))**2+eta_min

           geta_x = (4./Lxyz(1))*eta*((2.*x-(xyz0(1)+xyz1(1)))/Lxyz(1))
!
!  Single step function
!  Note that eta_x increases with increasing x when eta_xwidth is negative (!)
!
        case ('step')
!
!  default to spread gradient over ~5 grid cells.
!
           if (eta_xwidth == 0.) eta_xwidth = 5.*dx
           eta_x = eta + eta*(eta_jump-1.)*step(x,eta_x0,-eta_xwidth)
!
!  its gradient:
!  Note that geta_x points then only in the x direction.
!
           geta_x = eta*(eta_jump-1.)*der_step(x,eta_x0,-eta_xwidth)
!
        case ('RFP_1D')
!
           eta_x = eta*(1+9*(x/radRFP)**30)**2
!
! its gradient:
!
           geta_x = 2*eta*(1+9*(x/radRFP)**30)*270*(x/radRFP)**29/radRFP
!
!
!  Two-step function
!
        case ('two_step','two-step')
!
!  Allow for the each step to have its width. If they are
!  not specified, then eta_xwidth takes precedence.
!
           if (eta_xwidth .ne. 0.) then
             eta_xwidth0 =eta_xwidth
             eta_xwidth1 =eta_xwidth
           endif
!
!  Default to spread gradient over ~5 grid cells,
!
           if (eta_xwidth0 == 0.) eta_xwidth0 = 5.*dx
           if (eta_xwidth1 == 0.) eta_xwidth1 = 5.*dx
!
           eta_x = eta*eta_jump-eta*(eta_jump-two_step_factor)* &
             (step(x,eta_x0,eta_xwidth0)-step(x,eta_x1,eta_xwidth1))
!
!  ... and its gradient. Note that the sign of the second term enters
!  with the opposite sign, because we have used negative eta_xwidth.
!  Note that geta_x points then only in the x direction.
!
           geta_x = eta*(eta_jump-two_step_factor)*( &
             der_step(x,eta_x0,-eta_xwidth0)+der_step(x,eta_x1,eta_xwidth1))
!
!  Powerlaw-x
!
        case ('powerlaw-x','powerlaw_x')
!
!  eta proportional to chosen power of x...
!
          eta_x = eta*(1.+(x-x0)/eta_x0)**eta_power_x
!
!  ... and its gradient.
!
          geta_x = eta_power_x*eta_x/(eta_x0+x)
!
!  Powerlaw-x2: same profile as for viscosity (will produce eta=0 if 
!  x crosses zero)
!
        case ('powerlaw-x2','powerlaw_x2')
!
!  eta proportional to chosen power of x...
!
          eta_x = eta*(x/eta_x0)**eta_power_x
!
!  ... and its gradient.
!
          geta_x = eta_power_x*eta_x/eta_x0
!
      endselect
!
!  debug output (currently only on root processor)
!
      if (lroot.and.ldebug) then
        print*
        print*,'x, eta_x, geta_x'
        do l=l1,l2
          write(*,'(1p,3e11.3)') x(l),eta_x(l),geta_x(l)
        enddo
      endif
!
    endsubroutine eta_xdep
!***********************************************************************
    subroutine eta_rdep(eta_r,geta_r,rdep_profile)
!
!  creates a r-dependent resistivity
!
!  11-mar-2019/pete: adapted from eta_xdep
!
      use Sub, only: step, der_step
!
      real, dimension(nx) :: eta_r,tmp1,tmp2
      real, dimension(nx,3) :: geta_r
      type (pencil_case) :: p
      character (len=labellen), intent(in) :: rdep_profile
      integer :: l
!
      intent(out) :: eta_r,geta_r
!
      select case (rdep_profile)
!
!  Single step function
!  Note that eta_r increases with increasing r when eta_rwidth is negative (!)
!
        case ('step')
!
!  default to spread gradient over ~5 grid cells.
!
           if (eta_rwidth == 0.) eta_rwidth = 5.*dx
           tmp1=p%r_mn
           eta_r = eta + eta*(eta_jump-1.)*step(tmp1,eta_r0,-eta_rwidth)
!
!  its gradient:
!
           tmp2 = eta*(eta_jump-1.)*der_step(tmp1,eta_r0,-eta_rwidth)
           geta_r(:,1)=tmp2*x(l1:l2)*p%r_mn1
           geta_r(:,2)=tmp2*y(  m  )*p%r_mn1
           geta_r(:,3)=tmp2*z(  n  )*p%r_mn1
!
!  Two-step function
!
        case ('two_step','two-step')
!
!  Allow for the each step to have its width. If they are
!  not specified, then eta_xwidth takes precedence.
!
           if (eta_rwidth .ne. 0.) then
             eta_rwidth0 =eta_rwidth
             eta_rwidth1 =eta_rwidth
           endif
!
!  Default to spread gradient over ~5 grid cells,
!
           if (eta_rwidth0 == 0.) eta_rwidth0 = 5.*dx
           if (eta_rwidth1 == 0.) eta_rwidth1 = 5.*dx
!
           eta_r = eta*eta_jump-eta*(eta_jump-two_step_factor)* &
             (step(p%r_mn,eta_r0,eta_rwidth0)-step(p%r_mn,eta_r1,eta_rwidth1))
!
!  ... and its gradient. Note that the sign of the second term enters
!  with the opposite sign, because we have used negative eta_rwidth.
!
           tmp1 = eta*(eta_jump-two_step_factor)*( &
             der_step(p%r_mn,eta_r0,-eta_rwidth0)+der_step(p%r_mn,eta_r1,eta_rwidth1))
           geta_r(:,1)=tmp1*x(l1:l2)*p%r_mn1(l1:l2)
           geta_r(:,2)=tmp1*y(  m  )*p%r_mn1(l1:l2)
           geta_r(:,3)=tmp1*z(  n  )*p%r_mn1(l1:l2)
      endselect
!
!  debug output (currently only on root processor)
!
      if (lroot.and.ldebug) then
        print*
        print*,'p%r_mn, eta_r, geta_r'
        do l=l1,l2
          write(*,'(1p,3e11.3)') p%r_mn(l),eta_r(l),geta_r(l,1),geta_r(l,2),geta_r(l,3)
        enddo
      endif
!
    endsubroutine eta_rdep
!***********************************************************************
    subroutine bb_unitvec_shock(f,bb_hat)
!
!  Compute unit vector along the magnetic field.
!  Accurate to 2nd order.
!  Tries to avoid division by zero.
!  Taken from http://nuclear.llnl.gov/CNP/apt/apt/aptvunb.html.
!  If anybody knows a more accurate way of doing this, please modify.
!
!  16-aug-06/tobi: coded
!  19-jun-18/fred: added high order option with lbb_as_aux
!
      use Sub, only: dot2
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,3), intent (out) :: bb_hat
!
      !Tobi: Not sure about this value
      real, parameter :: tol=1e-11
!
      real, dimension (mx,3) :: bb,bb2
      real, dimension (mx) :: bb_len,aerr2
      real :: fac
      integer :: j
!
!  Compute magnetic field from vector potential.
!
      bb=0.
!
      if (.not. lbb_as_aux) then
!
        if (nxgrid/=1) then
          fac = 1/(2*dx)
          bb(l1-2:l2+2,3) = bb(l1-2:l2+2,3) + fac*( f(l1-1:l2+3,m  ,n  ,iay)   &
                                                  - f(l1-3:l2+1,m  ,n  ,iay) )
          bb(l1-2:l2+2,2) = bb(l1-2:l2+2,2) - fac*( f(l1-1:l2+3,m  ,n  ,iaz)   &
                                                  - f(l1-3:l2+1,m  ,n  ,iaz) )
        endif
!
        if (nygrid/=1) then
          fac = 1/(2*dy)
          bb(l1-2:l2+2,1) = bb(l1-2:l2+2,1) + fac*( f(l1-2:l2+2,m+1,n  ,iaz)   &
                                                  - f(l1-2:l2+2,m-1,n  ,iaz) )
          bb(l1-2:l2+2,3) = bb(l1-2:l2+2,3) - fac*( f(l1-2:l2+2,m+1,n  ,iax)   &
                                                  - f(l1-2:l2+2,m-1,n  ,iax) )
        endif
!
        if (nzgrid/=1) then
          fac = 1/(2*dz)
          bb(l1-2:l2+2,2) = bb(l1-2:l2+2,2) + fac*( f(l1-2:l2+2,m  ,n+1,iax)   &
                                                  - f(l1-2:l2+2,m  ,n-1,iax) )
          bb(l1-2:l2+2,1) = bb(l1-2:l2+2,1) - fac*( f(l1-2:l2+2,m  ,n+1,iay)   &
                                                  - f(l1-2:l2+2,m  ,n-1,iay) )
        endif
!
!  Add external magnetic field.
!
        do j=1,3; bb(:,j) = bb(:,j) + B_ext(j); enddo
!
!  Truncate small components to zero.
!
        bb2 = bb**2
!
        aerr2 = tol**2 * max(sum(bb2,2),1.)
!
        do j=1,3
          where (bb2(:,j) < aerr2)
            bb_hat(:,j) = 0.
          elsewhere
            bb_hat(:,j) = bb(:,j)
          endwhere
        enddo
!
!  Get unit vector.
!
        bb_len = sqrt(sum(bb_hat**2,2))
!
        do j=1,3; bb_hat(:,j) = bb_hat(:,j)/(bb_len+tini); enddo
!
! else if lbb_as_aux
!
      else
        bb(l1:l2,:) = f(l1:l2,m,n,ibx:ibz)
        call dot2(bb(l1:l2,:),bb_len(l1:l2),PRECISE_SQRT=.true.)
        aerr2 = tol * max(bb_len,1.)
!
!  Truncate small components to zero.
!
        do j=1,3
          where (bb_len(l1:l2) <= aerr2(l1:l2))
            bb_hat(l1:l2,j) = 0.
          elsewhere
            bb_hat(l1:l2,j) = bb(l1:l2,j)
          endwhere
        enddo
!
!  Get unit vector.
!
        do j=1,3; bb_hat(l1:l2,j) = bb_hat(l1:l2,j)/(bb_len(l1:l2)+tini); enddo
!
      endif
!
    endsubroutine bb_unitvec_shock
!***********************************************************************
    subroutine input_persist_magnetic_id(id,done)
!
!  Read in the stored phase and amplitude for the correction of the Beltrami
!  wave forcing.
!
!   5-apr-08/axel: adapted from input_persist_forcing
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use IO, only: read_persist
!
      integer :: id
      logical :: done
!
      select case (id)
        case (id_record_MAGNETIC_PHASE)
          if (read_persist ('MAGNETIC_PHASE', phase_beltrami)) return
          if (lroot) print *, 'input_persist_magnetic: phase_beltrami = ', phase_beltrami
          done = .true.
        case (id_record_MAGNETIC_AMPL)
          if (read_persist ('MAGNETIC_AMPL', ampl_beltrami)) return
          if (lroot) print *, 'input_persist_magnetic: ampl_beltrami = ', ampl_beltrami
          done = .true.
      endselect
!
    endsubroutine input_persist_magnetic_id
!***********************************************************************
    subroutine input_persist_magnetic()
!
!  Read in the stored phase and amplitude for the correction of the Beltrami
!  wave forcing.
!
!  12-Oct-2019/PABourdin: coded
!
      use IO, only: read_persist
!
      logical :: error
!
      error = read_persist ('MAGNETIC_PHASE', phase_beltrami)
      if (lroot .and. .not. error) print *, 'input_persist_magnetic: phase_beltrami = ', phase_beltrami
!
      error = read_persist ('MAGNETIC_AMPL', ampl_beltrami)
      if (lroot .and. .not. error) print *, 'input_persist_magnetic: ampl_beltrami = ', ampl_beltrami
!
    endsubroutine input_persist_magnetic
!***********************************************************************
    logical function output_persistent_magnetic()
!
!  Write the stored phase and amplitude for the
!  correction of the Beltrami wave forcing
!
!    5-apr-08/axel: adapted from output_persistent_forcing
!   13-Dec-2011/Bourdin.KIS: reworked
!
      use IO, only: write_persist
!
      if (lroot .and. (ip < 14) .and. lforcing_cont_aa_local) then
        print *, 'output_persistent_magnetic: ', phase_beltrami, ampl_beltrami
      endif
!
!  write details
!
      output_persistent_magnetic = .true.
!
      if (lforcing_cont_aa_local) then
        if (write_persist ('MAGNETIC_PHASE', id_record_MAGNETIC_PHASE, phase_beltrami)) return
        if (write_persist ('MAGNETIC_AMPL', id_record_MAGNETIC_AMPL, ampl_beltrami)) return
      endif
!
      output_persistent_magnetic = .false.
!
    endfunction output_persistent_magnetic
!***********************************************************************
    subroutine rprint_magnetic(lreset,lwrite)
!
!  Reads and registers print parameters relevant for magnetic fields.
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics
!
      integer :: iname,inamex,inamey,inamez,ixy,ixz,irz,inamer,iname_half,iname_sound,inamev
      logical :: lreset
      logical, intent(in), optional :: lwrite
!
      integer :: idum
!
!  Reset everything in case of RELOAD.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_eta_tdep=0
        idiag_ab_int=0; idiag_jb_int=0; idiag_b2tm=0; idiag_bjtm=0; idiag_jbtm=0
        idiag_b2uzm=0; idiag_b2ruzm=0; idiag_ubbzm=0
        idiag_b1m=0; idiag_b2m=0; idiag_EEM=0; idiag_b4m=0
        idiag_bm2=0; idiag_j2m=0; idiag_jm2=0
        idiag_abm=0; idiag_abrms=0; idiag_jbrms=0; idiag_abmh=0
        idiag_abumx=0; idiag_abumy=0; idiag_abumz=0
        idiag_abmn=0; idiag_abms=0; idiag_jbmh=0; idiag_jbmn=0; idiag_jbms=0
        idiag_ajm=0; idiag_cosubm=0; idiag_jbm=0; idiag_hjbm=0
        idiag_uam=0; idiag_ubm=0; idiag_dubrms=0; idiag_dobrms=0; idiag_ujm=0
        idiag_uxbxm=0; idiag_uybxm=0; idiag_uzbxm=0
        idiag_uxbym=0; idiag_uybym=0; idiag_uzbym=0
        idiag_uxbzm=0; idiag_uybzm=0; idiag_uzbzm=0
        idiag_jxbxm=0; idiag_jybxm=0; idiag_jzbxm=0
        idiag_jxbym=0; idiag_jybym=0; idiag_jzbym=0
        idiag_jxbzm=0; idiag_jybzm=0; idiag_jzbzm=0
        idiag_uxjxm=0; idiag_uyjxm=0; idiag_uzjxm=0
        idiag_uxjym=0; idiag_uyjym=0; idiag_uzjym=0
        idiag_uxjzm=0; idiag_uyjzm=0; idiag_uzjzm=0
        idiag_fbm=0; idiag_fxbxm=0; idiag_epsM=0; idiag_epsM_LES=0
        idiag_epsAD=0; idiag_epsMmz=0
        idiag_bxpt=0; idiag_bypt=0; idiag_bzpt=0
        idiag_jxpt=0; idiag_jypt=0; idiag_jzpt=0
        idiag_Expt=0; idiag_Eypt=0; idiag_Ezpt=0
        idiag_axpt=0; idiag_aypt=0; idiag_azpt=0
        idiag_bxp2=0; idiag_byp2=0; idiag_bzp2=0
        idiag_jxp2=0; idiag_jyp2=0; idiag_jzp2=0
        idiag_Exp2=0; idiag_Eyp2=0; idiag_Ezp2=0
        idiag_axp2=0; idiag_ayp2=0; idiag_azp2=0
        idiag_aybym2=0; idiag_exaym2=0; idiag_exjm2=0
        idiag_brms=0; idiag_bfrms=0; idiag_bf2m=0; idiag_bf4m=0
        idiag_bmax=0; idiag_jrms=0; idiag_jmax=0
        idiag_vArms=0; idiag_emag=0; idiag_bxmin=0; idiag_bymin=0; idiag_bzmin=0
        idiag_bxmax=0; idiag_bymax=0; idiag_bzmax=0; idiag_vAmax=0; idiag_dtb=0
        idiag_dtFr=0;idiag_dtHr=0;idiag_dtBr=0;
        idiag_bbxmax=0; idiag_bbymax=0; idiag_bbzmax=0
        idiag_jxmax=0; idiag_jymax=0; idiag_jzmax=0
        idiag_a2m=0; idiag_arms=0; idiag_amax=0; idiag_beta1m=0; idiag_beta1mz=0
        idiag_divarms = 0
        idiag_beta1max=0; idiag_bxm=0; idiag_bym=0; idiag_bzm=0; idiag_axm=0
        idiag_betam = 0; idiag_betamax = 0; idiag_betamin = 0
        idiag_betamz = 0; idiag_beta2mz = 0
        idiag_betamx = 0; idiag_beta2mx = 0
        idiag_aym=0; idiag_azm=0
        idiag_bx2m=0; idiag_by2m=0; idiag_bz2m=0
        idiag_bx4m=0; idiag_by4m=0; idiag_bz4m=0
        idiag_jx2m=0; idiag_jy2m=0; idiag_jz2m=0
        idiag_jx4m=0; idiag_jy4m=0; idiag_jz4m=0
        idiag_bxbymx = 0; idiag_bxbzmx = 0; idiag_bybzmx = 0
        idiag_bxbymy=0; idiag_bxbzmy=0; idiag_bybzmy=0; idiag_bxbymz=0
        idiag_bxbzmz=0; idiag_bybzmz=0
        idiag_b2mx=0; idiag_b2mz=0; idiag_bf2mz=0; idiag_j2mz=0
        idiag_jbmz=0; idiag_abmz=0; idiag_ubmz=0; idiag_uamz=0
        idiag_bzdivamz=0; idiag_divamz=0; idiag_d6abmz=0
        idiag_uxbxmz=0; idiag_uybxmz=0; idiag_uzbxmz=0
        idiag_uxbymz=0; idiag_uybymz=0; idiag_uzbymz=0
        idiag_uxbzmz=0; idiag_uybzmz=0; idiag_uzbzmz=0
        idiag_d6amz3=0; idiag_d6amz2=0; idiag_d6amz1=0; idiag_poynzmz=0
        idiag_bxbym=0; idiag_bxbzm=0; idiag_bybzm=0; idiag_djuidjbim=0
        idiag_axmz=0; idiag_aymz=0; idiag_azmz=0; idiag_bxmz=0; idiag_bymz=0
        idiag_jxmz=0; idiag_jymz=0; idiag_jzmz=0
        idiag_abuxmz=0; idiag_abuymz=0; idiag_abuzmz=0
        idiag_uabxmz=0; idiag_uabymz=0; idiag_uabzmz=0
        idiag_bzmz=0; idiag_bmx=0; idiag_bmy=0; idiag_bmz=0; idiag_embmz=0
        idiag_bmzS2=0; idiag_bmzA2=0
        idiag_emxamz3=0; idiag_jmx=0; idiag_jmy=0; idiag_jmz=0; idiag_ambmz=0
        idiag_jmbmz=0; idiag_kmz=0; idiag_kx_aa=0
        idiag_ambmzh=0;idiag_ambmzn=0;idiag_ambmzs=0; idiag_bmzph=0
        idiag_bmzphe=0; idiag_bcosphz=0; idiag_bsinphz=0
        idiag_bx2mz=0; idiag_by2mz=0; idiag_bz2mz=0
        idiag_bx2rmz=0; idiag_by2rmz=0; idiag_bz2rmz=0
        idiag_bxmxy=0; idiag_bymxy=0; idiag_bzmxy=0
        idiag_dbxdxmxy=0; idiag_dbxdymxy=0; idiag_dbxdzmxy=0
        idiag_dbydxmxy=0; idiag_dbydymxy=0; idiag_dbydzmxy=0
        idiag_dbzdxmxy=0; idiag_dbzdymxy=0; idiag_dbzdzmxy=0
        idiag_jxmxy=0; idiag_jymxy=0; idiag_jzmxy=0
        idiag_poynxmxy=0; idiag_poynymxy=0; idiag_poynzmxy=0
        idiag_bx2mxy=0; idiag_by2mxy=0; idiag_bz2mxy=0; idiag_bxbymxy=0
        idiag_examxy1=0; idiag_examxy3=0; idiag_examxy2=0
        idiag_bxbzmxy=0; idiag_bybzmxy=0; idiag_bxbymxz=0; idiag_bxbzmxz=0
        idiag_Exmz=0 ; idiag_Eymz=0; idiag_Ezmz=0
        idiag_Exmxy=0 ; idiag_Eymxy=0; idiag_Ezmxy=0; idiag_beta1mxy=0
        idiag_StokesImxy=0; idiag_StokesQmxy=0; idiag_StokesUmxy=0
        idiag_StokesQ1mxy=0; idiag_StokesU1mxy=0
        idiag_bybzmxz=0; idiag_uybxmxz=0; idiag_uybzmxz=0
        idiag_bx1mxz=0; idiag_by1mxz=0; idiag_bz1mxz=0
        idiag_bxmxz=0; idiag_bymxz=0; idiag_bzmxz=0; idiag_jbmxy=0
        idiag_jxmxz=0; idiag_jymxz=0; idiag_jzmxz=0
        idiag_abmxy=0; idiag_b2mxz=0; idiag_ubmxy=0;
        idiag_axmxz=0; idiag_aymxz=0; idiag_azmxz=0; idiag_Exmxz=0
        idiag_axmxy=0; idiag_aymxy=0; idiag_azmxy=0
        idiag_Eymxz=0; idiag_Ezmxz=0; idiag_jxbm=0; idiag_uxbm=0; idiag_oxuxbm=0
        idiag_jxbxbm=0; idiag_gpxbm=0; idiag_uxDxuxbm=0; idiag_vAmxz=0
        idiag_uxbmx=0; idiag_uxbmy=0; idiag_uxbmz=0
        idiag_jxbmx=0; idiag_jxbmy=0; idiag_jxbmz=0
        idiag_uxbcmx=0; idiag_uxbsmx=0
        idiag_uxbcmy=0; idiag_uxbsmy=0; idiag_examz1=0; idiag_examz2=0
        idiag_examz3=0; idiag_examx=0; idiag_examy=0; idiag_examz=0
        idiag_e3xamz1=0; idiag_e3xamz2=0; idiag_e3xamz3=0
        idiag_exjmx=0; idiag_exjmy=0; idiag_exjmz=0; idiag_dexbmx=0
        idiag_dexbmy=0; idiag_dexbmz=0; idiag_phibmx=0; idiag_phibmy=0
        idiag_phibmz=0; idiag_uxjm=0; idiag_jdel2am=0
        idiag_ujxbm=0; idiag_b2divum=0
        idiag_b3b21m=0; idiag_b3b12m=0; idiag_b1b32m=0; idiag_b1b23m=0
        idiag_b2b13m=0 ; idiag_b2b31m=0
        idiag_udotxbm=0; idiag_uxbdotm=0; idiag_brmphi=0; idiag_bpmphi=0
        idiag_bzmphi=0; idiag_b2mphi=0; idiag_jbmphi=0; idiag_uxbrmphi=0
        idiag_uxbpmphi=0; idiag_uxbzmphi=0; idiag_jxbrmphi=0; idiag_jxbpmphi=0
        idiag_jxbzmphi=0; idiag_jxbrxm=0; idiag_jxbrym=0; idiag_jxbrzm=0
        idiag_jxbr2m=0; idiag_jxbrxmx=0; idiag_jxbrymx=0; idiag_jxbrzmx=0; idiag_jxbrmax=0;
        idiag_jxbrxmy=0; idiag_jxbrymy=0; idiag_jxbrzmy=0; idiag_jxbrxmz=0
        idiag_jxbrymz=0; idiag_jxbrzmz=0; idiag_armphi=0; idiag_apmphi=0
        idiag_azmphi=0; idiag_dteta=0; idiag_uxBrms=0; idiag_Bresrms=0
        idiag_Rmrms=0; idiag_jfm=0; idiag_brbpmr=0; idiag_va2m=0; idiag_b2mr=0
        idiag_brmr=0; idiag_bpmr=0; idiag_bzmr=0; idiag_armr=0; idiag_apmr=0
        idiag_azmr=0; idiag_bxmx=0; idiag_bymx=0; idiag_bzmx=0; idiag_bxmy=0
        idiag_bymy=0; idiag_bzmy=0; idiag_bx2my=0; idiag_by2my=0; idiag_bz2my=0
        idiag_mflux_x=0; idiag_mflux_y=0; idiag_mflux_z=0; idiag_bmxy_rms=0
        idiag_brsphmphi=0; idiag_bthmphi=0; idiag_brmsh=0; idiag_brmsn=0
        idiag_brmss=0; idiag_etatotalmx=0; idiag_etatotalmz=0; idiag_Rmmz=0
        idiag_brmsx=0; idiag_brmsz=0
        idiag_etavamax=0; idiag_etajmax=0; idiag_etaj2max=0; idiag_etajrhomax=0
        idiag_hjrms=0;idiag_hjbm=0;idiag_coshjbm=0
        idiag_etasmagm=0; idiag_etaaniso=0
        idiag_cosjbm=0;idiag_jparallelm=0;idiag_jperpm=0
        idiag_hjparallelm=0;idiag_hjperpm=0
        idiag_vmagfricmax=0; idiag_vmagfricrms=0; idiag_vmagfricmz=0
        ivid_aps=0; ivid_bb=0; ivid_jj=0; ivid_b2=0; ivid_j2=0; ivid_ab=0
        ivid_jb=0; ivid_beta1=0; ivid_poynting=0; ivid_bb_sph=0; idiag_dteta3=0
      endif
!
!  Check for those quantities that we want to evaluate online.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'eta_tdep',idiag_eta_tdep)
        call parse_name(iname,cname(iname),cform(iname),'ab_int',idiag_ab_int)
        call parse_name(iname,cname(iname),cform(iname),'jb_int',idiag_jb_int)
        call parse_name(iname,cname(iname),cform(iname),'dteta',idiag_dteta)
        call parse_name(iname,cname(iname),cform(iname),'dteta3',idiag_dteta3)
        call parse_name(iname,cname(iname),cform(iname),'aybym2',idiag_aybym2)
        call parse_name(iname,cname(iname),cform(iname),'exaym2',idiag_exaym2)
        call parse_name(iname,cname(iname),cform(iname),'exabot',idiag_exabot)
        call parse_name(iname,cname(iname),cform(iname),'exatop',idiag_exatop)
        call parse_name(iname,cname(iname),cform(iname),'exjm2',idiag_exjm2)
        call parse_name(iname,cname(iname),cform(iname),'b2tm',idiag_b2tm)
        call parse_name(iname,cname(iname),cform(iname),'bjtm',idiag_bjtm)
        call parse_name(iname,cname(iname),cform(iname),'jbtm',idiag_jbtm)
        call parse_name(iname,cname(iname),cform(iname),'abm',idiag_abm)
        call parse_name(iname,cname(iname),cform(iname),'abmn',idiag_abmn)
        call parse_name(iname,cname(iname),cform(iname),'abms',idiag_abms)
        call parse_name(iname,cname(iname),cform(iname),'abrms',idiag_abrms)
        call parse_name(iname,cname(iname),cform(iname),'jbrms',idiag_jbrms)
        call parse_name(iname,cname(iname),cform(iname),'abumx',idiag_abumx)
        call parse_name(iname,cname(iname),cform(iname),'abumy',idiag_abumy)
        call parse_name(iname,cname(iname),cform(iname),'abumz',idiag_abumz)
        call parse_name(iname,cname(iname),cform(iname),'ajm',idiag_ajm)
        call parse_name(iname,cname(iname),cform(iname),'jbm',idiag_jbm)
        call parse_name(iname,cname(iname),cform(iname),'hjbm',idiag_hjbm)
        call parse_name(iname,cname(iname),cform(iname),'jbmn',idiag_jbmn)
        call parse_name(iname,cname(iname),cform(iname),'jbms',idiag_jbms)
        call parse_name(iname,cname(iname),cform(iname),'ubm',idiag_ubm)
        call parse_name(iname,cname(iname),cform(iname),'dubrms',idiag_dubrms)
        call parse_name(iname,cname(iname),cform(iname),'dobrms',idiag_dobrms)
        call parse_name(iname,cname(iname),cform(iname),'uxbxm',idiag_uxbxm)
        call parse_name(iname,cname(iname),cform(iname),'uybxm',idiag_uybxm)
        call parse_name(iname,cname(iname),cform(iname),'uzbxm',idiag_uzbxm)
        call parse_name(iname,cname(iname),cform(iname),'uxbym',idiag_uxbym)
        call parse_name(iname,cname(iname),cform(iname),'uybym',idiag_uybym)
        call parse_name(iname,cname(iname),cform(iname),'uzbym',idiag_uzbym)
        call parse_name(iname,cname(iname),cform(iname),'uxbzm',idiag_uxbzm)
        call parse_name(iname,cname(iname),cform(iname),'uybzm',idiag_uybzm)
        call parse_name(iname,cname(iname),cform(iname),'uzbzm',idiag_uzbzm)
        call parse_name(iname,cname(iname),cform(iname),'uxjxm',idiag_uxjxm)
        call parse_name(iname,cname(iname),cform(iname),'uyjxm',idiag_uyjxm)
        call parse_name(iname,cname(iname),cform(iname),'uzjxm',idiag_uzjxm)
        call parse_name(iname,cname(iname),cform(iname),'uxjym',idiag_uxjym)
        call parse_name(iname,cname(iname),cform(iname),'uyjym',idiag_uyjym)
        call parse_name(iname,cname(iname),cform(iname),'uzjym',idiag_uzjym)
        call parse_name(iname,cname(iname),cform(iname),'uxjzm',idiag_uxjzm)
        call parse_name(iname,cname(iname),cform(iname),'uyjzm',idiag_uyjzm)
        call parse_name(iname,cname(iname),cform(iname),'uzjzm',idiag_uzjzm)
        call parse_name(iname,cname(iname),cform(iname),'cosubm',idiag_cosubm)
        call parse_name(iname,cname(iname),cform(iname),'jxbxm',idiag_jxbxm)
        call parse_name(iname,cname(iname),cform(iname),'jybxm',idiag_jybxm)
        call parse_name(iname,cname(iname),cform(iname),'jzbxm',idiag_jzbxm)
        call parse_name(iname,cname(iname),cform(iname),'jxbym',idiag_jxbym)
        call parse_name(iname,cname(iname),cform(iname),'jybym',idiag_jybym)
        call parse_name(iname,cname(iname),cform(iname),'jzbym',idiag_jzbym)
        call parse_name(iname,cname(iname),cform(iname),'jxbzm',idiag_jxbzm)
        call parse_name(iname,cname(iname),cform(iname),'jybzm',idiag_jybzm)
        call parse_name(iname,cname(iname),cform(iname),'jzbzm',idiag_jzbzm)
        call parse_name(iname,cname(iname),cform(iname),'uam',idiag_uam)
        call parse_name(iname,cname(iname),cform(iname),'ujm',idiag_ujm)
        call parse_name(iname,cname(iname),cform(iname),'fbm',idiag_fbm)
        call parse_name(iname,cname(iname),cform(iname),'fxbxm',idiag_fxbxm)
        call parse_name(iname,cname(iname),cform(iname),'b2ruzm',idiag_b2ruzm)
        call parse_name(iname,cname(iname),cform(iname),'b2uzm',idiag_b2uzm)
        call parse_name(iname,cname(iname),cform(iname),'ubbzm',idiag_ubbzm)
        call parse_name(iname,cname(iname),cform(iname),'b1m',idiag_b1m)
        call parse_name(iname,cname(iname),cform(iname),'b2m',idiag_b2m)
        call parse_name(iname,cname(iname),cform(iname),'EEM',idiag_EEM)
        call parse_name(iname,cname(iname),cform(iname),'b4m',idiag_b4m)
        call parse_name(iname,cname(iname),cform(iname),'bm2',idiag_bm2)
        call parse_name(iname,cname(iname),cform(iname),'j2m',idiag_j2m)
        call parse_name(iname,cname(iname),cform(iname),'jm2',idiag_jm2)
        call parse_name(iname,cname(iname),cform(iname),'epsM',idiag_epsM)
        call parse_name(iname,cname(iname),cform(iname),&
            'epsM_LES',idiag_epsM_LES)
        call parse_name(iname,cname(iname),cform(iname),'epsAD',idiag_epsAD)
        call parse_name(iname,cname(iname),cform(iname),'emag',idiag_emag)
        call parse_name(iname,cname(iname),cform(iname),'brms',idiag_brms)
        call parse_name(iname,cname(iname),cform(iname),'bfrms',idiag_bfrms)
        call parse_name(iname,cname(iname),cform(iname),'bf2m',idiag_bf2m)
        call parse_name(iname,cname(iname),cform(iname),'bf4m',idiag_bf4m)
        call parse_name(iname,cname(iname),cform(iname),'brmsn',idiag_brmsn)
        call parse_name(iname,cname(iname),cform(iname),'brmss',idiag_brmss)
        call parse_name(iname,cname(iname),cform(iname),'brmsx',idiag_brmsx)
        call parse_name(iname,cname(iname),cform(iname),'brmsz',idiag_brmsz)
        call parse_name(iname,cname(iname),cform(iname),'bmax',idiag_bmax)
        call parse_name(iname,cname(iname),cform(iname),'bxmin',idiag_bxmin)
        call parse_name(iname,cname(iname),cform(iname),'bymin',idiag_bymin)
        call parse_name(iname,cname(iname),cform(iname),'bzmin',idiag_bzmin)
        call parse_name(iname,cname(iname),cform(iname),'bxmax',idiag_bxmax)
        call parse_name(iname,cname(iname),cform(iname),'bymax',idiag_bymax)
        call parse_name(iname,cname(iname),cform(iname),'bzmax',idiag_bzmax)
        call parse_name(iname,cname(iname),cform(iname),'bbxmax',idiag_bbxmax)
        call parse_name(iname,cname(iname),cform(iname),'bbymax',idiag_bbymax)
        call parse_name(iname,cname(iname),cform(iname),'bbzmax',idiag_bbzmax)
        call parse_name(iname,cname(iname),cform(iname),'jxmax',idiag_jxmax)
        call parse_name(iname,cname(iname),cform(iname),'jymax',idiag_jymax)
        call parse_name(iname,cname(iname),cform(iname),'jzmax',idiag_jzmax)
        call parse_name(iname,cname(iname),cform(iname),'jrms',idiag_jrms)
        call parse_name(iname,cname(iname),cform(iname),'hjrms',idiag_hjrms)
        call parse_name(iname,cname(iname),cform(iname),'jmax',idiag_jmax)
        call parse_name(iname,cname(iname),cform(iname),'axm',idiag_axm)
        call parse_name(iname,cname(iname),cform(iname),'aym',idiag_aym)
        call parse_name(iname,cname(iname),cform(iname),'azm',idiag_azm)
        call parse_name(iname,cname(iname),cform(iname),'a2m',idiag_a2m)
        call parse_name(iname,cname(iname),cform(iname),'arms',idiag_arms)
        call parse_name(iname,cname(iname),cform(iname),'amax',idiag_amax)
        call parse_name(iname,cname(iname),cform(iname),'divarms',idiag_divarms)
        call parse_name(iname,cname(iname),cform(iname),'vArms',idiag_vArms)
        call parse_name(iname,cname(iname),cform(iname),'vAmax',idiag_vAmax)
        call parse_name(iname,cname(iname),cform(iname),'vA2m',idiag_vA2m)
        call parse_name(iname,cname(iname),cform(iname),'beta1m',idiag_beta1m)
        call parse_name(iname,cname(iname),cform(iname),'beta1max',idiag_beta1max)
        call parse_name(iname,cname(iname),cform(iname),'betam',idiag_betam)
        call parse_name(iname,cname(iname),cform(iname),'betamax',idiag_betamax)
        call parse_name(iname,cname(iname),cform(iname),'betamin',idiag_betamin)
        call parse_name(iname,cname(iname),cform(iname),'dtb',idiag_dtb)
        call parse_name(iname,cname(iname),cform(iname),'dtHr',idiag_dtHr)
        call parse_name(iname,cname(iname),cform(iname),'dtFr',idiag_dtFr)
        call parse_name(iname,cname(iname),cform(iname),'dtBr',idiag_dtBr)
        call parse_name(iname,cname(iname),cform(iname),'bxm',idiag_bxm)
        call parse_name(iname,cname(iname),cform(iname),'bym',idiag_bym)
        call parse_name(iname,cname(iname),cform(iname),'bzm',idiag_bzm)
        call parse_name(iname,cname(iname),cform(iname),'bx2m',idiag_bx2m)
        call parse_name(iname,cname(iname),cform(iname),'by2m',idiag_by2m)
        call parse_name(iname,cname(iname),cform(iname),'bz2m',idiag_bz2m)
        call parse_name(iname,cname(iname),cform(iname),'bx4m',idiag_bx4m)
        call parse_name(iname,cname(iname),cform(iname),'by4m',idiag_by4m)
        call parse_name(iname,cname(iname),cform(iname),'bz4m',idiag_bz4m)
        call parse_name(iname,cname(iname),cform(iname),'jx2m',idiag_jx2m)
        call parse_name(iname,cname(iname),cform(iname),'jy2m',idiag_jy2m)
        call parse_name(iname,cname(iname),cform(iname),'jz2m',idiag_jz2m)
        call parse_name(iname,cname(iname),cform(iname),'jx4m',idiag_jx4m)
        call parse_name(iname,cname(iname),cform(iname),'jy4m',idiag_jy4m)
        call parse_name(iname,cname(iname),cform(iname),'jz4m',idiag_jz4m)
        call parse_name(iname,cname(iname),cform(iname),'bxbym',idiag_bxbym)
        call parse_name(iname,cname(iname),cform(iname),'bxbzm',idiag_bxbzm)
        call parse_name(iname,cname(iname),cform(iname),'bybzm',idiag_bybzm)
        call parse_name(iname,cname(iname),cform(iname),'djuidjbim',idiag_djuidjbim)
        call parse_name(iname,cname(iname),cform(iname),'jxbrxm',idiag_jxbrxm)
        call parse_name(iname,cname(iname),cform(iname),'jxbrym',idiag_jxbrym)
        call parse_name(iname,cname(iname),cform(iname),'jxbrzm',idiag_jxbrzm)
        call parse_name(iname,cname(iname),cform(iname),'jxbrmax',idiag_jxbrmax)
        call parse_name(iname,cname(iname),cform(iname),'jxbr2m',idiag_jxbr2m)
        call parse_name(iname,cname(iname),cform(iname),'jxbm',idiag_jxbm)
        call parse_name(iname,cname(iname),cform(iname),'uxbm',idiag_uxbm)
        call parse_name(iname,cname(iname),cform(iname),'uxbmx',idiag_uxbmx)
        call parse_name(iname,cname(iname),cform(iname),'uxbmy',idiag_uxbmy)
        call parse_name(iname,cname(iname),cform(iname),'uxbmz',idiag_uxbmz)
        call parse_name(iname,cname(iname),cform(iname),'jxbmx',idiag_jxbmx)
        call parse_name(iname,cname(iname),cform(iname),'jxbmy',idiag_jxbmy)
        call parse_name(iname,cname(iname),cform(iname),'jxbmz',idiag_jxbmz)
        call parse_name(iname,cname(iname),cform(iname),'vmagfricmax',&
                       idiag_vmagfricmax)
        call parse_name(iname,cname(iname),cform(iname),'vmagfricrms',&
                       idiag_vmagfricrms)
        call parse_name(iname,cname(iname),cform(iname),'uxbcmx',idiag_uxbcmx)
        call parse_name(iname,cname(iname),cform(iname),'uxbcmy',idiag_uxbcmy)
        call parse_name(iname,cname(iname),cform(iname),'uxbsmx',idiag_uxbsmx)
        call parse_name(iname,cname(iname),cform(iname),'uxbsmy',idiag_uxbsmy)
        call parse_name(iname,cname(iname),cform(iname),'examx',idiag_examx)
        call parse_name(iname,cname(iname),cform(iname),'examy',idiag_examy)
        call parse_name(iname,cname(iname),cform(iname),'examz',idiag_examz)
        call parse_name(iname,cname(iname),cform(iname),'exjmx',idiag_exjmx)
        call parse_name(iname,cname(iname),cform(iname),'exjmy',idiag_exjmy)
        call parse_name(iname,cname(iname),cform(iname),'exjmz',idiag_exjmz)
        call parse_name(iname,cname(iname),cform(iname),'dexbmx',idiag_dexbmx)
        call parse_name(iname,cname(iname),cform(iname),'dexbmy',idiag_dexbmy)
        call parse_name(iname,cname(iname),cform(iname),'dexbmz',idiag_dexbmz)
        call parse_name(iname,cname(iname),cform(iname),'phibmx',idiag_phibmx)
        call parse_name(iname,cname(iname),cform(iname),'phibmy',idiag_phibmy)
        call parse_name(iname,cname(iname),cform(iname),'phibmz',idiag_phibmz)
        call parse_name(iname,cname(iname),cform(iname),'uxjm',idiag_uxjm)
        call parse_name(iname,cname(iname),cform(iname),'jdel2am',idiag_jdel2am)
        call parse_name(iname,cname(iname),cform(iname),'ujxbm',idiag_ujxbm)
        call parse_name(iname,cname(iname),cform(iname),'b2divum',idiag_b2divum)
        call parse_name(iname,cname(iname),cform(iname),'jxbxbm',idiag_jxbxbm)
        call parse_name(iname,cname(iname),cform(iname),'oxuxbm',idiag_oxuxbm)
        call parse_name(iname,cname(iname),cform(iname),'gpxbm',idiag_gpxbm)
        call parse_name(iname,cname(iname),cform(iname),&
            'uxDxuxbm',idiag_uxDxuxbm)
        call parse_name(iname,cname(iname),cform(iname),'b3b21m',idiag_b3b21m)
        call parse_name(iname,cname(iname),cform(iname),'b3b12m',idiag_b3b12m)
        call parse_name(iname,cname(iname),cform(iname),'b1b32m',idiag_b1b32m)
        call parse_name(iname,cname(iname),cform(iname),'b1b23m',idiag_b1b23m)
        call parse_name(iname,cname(iname),cform(iname),'b2b13m',idiag_b2b13m)
        call parse_name(iname,cname(iname),cform(iname),'b2b31m',idiag_b2b31m)
        call parse_name(iname,cname(iname),cform(iname),'udotxbm',idiag_udotxbm)
        call parse_name(iname,cname(iname),cform(iname),'uxbdotm',idiag_uxbdotm)
        call parse_name(iname,cname(iname),cform(iname),'bmx',idiag_bmx)
        call parse_name(iname,cname(iname),cform(iname),'bmy',idiag_bmy)
        call parse_name(iname,cname(iname),cform(iname),'bmz',idiag_bmz)
        call parse_name(iname,cname(iname),cform(iname),'bmzS2',idiag_bmzS2)
        call parse_name(iname,cname(iname),cform(iname),'bmzA2',idiag_bmzA2)
        call parse_name(iname,cname(iname),cform(iname),'jmx',idiag_jmx)
        call parse_name(iname,cname(iname),cform(iname),'jmy',idiag_jmy)
        call parse_name(iname,cname(iname),cform(iname),'jmz',idiag_jmz)
        call parse_name(iname,cname(iname),cform(iname),'bmzph',idiag_bmzph)
        call parse_name(iname,cname(iname),cform(iname),'bmzphe',idiag_bmzphe)
        call parse_name(iname,cname(iname),cform(iname),'bcosphz',idiag_bcosphz)
        call parse_name(iname,cname(iname),cform(iname),'bsinphz',idiag_bsinphz)
        call parse_name(iname,cname(iname),cform(iname),'emxamz3',idiag_emxamz3)
        call parse_name(iname,cname(iname),cform(iname),'embmz',idiag_embmz)
        call parse_name(iname,cname(iname),cform(iname),'ambmz',idiag_ambmz)
        call parse_name(iname,cname(iname),cform(iname),'ambmzn',idiag_ambmzn)
        call parse_name(iname,cname(iname),cform(iname),'ambmzs',idiag_ambmzs)
        call parse_name(iname,cname(iname),cform(iname),'jmbmz',idiag_jmbmz)
        call parse_name(iname,cname(iname),cform(iname),'kmz',idiag_kmz)
        call parse_name(iname,cname(iname),cform(iname),'kx_aa',idiag_kx_aa)
        call parse_name(iname,cname(iname),cform(iname),'bxpt',idiag_bxpt)
        call parse_name(iname,cname(iname),cform(iname),'bypt',idiag_bypt)
        call parse_name(iname,cname(iname),cform(iname),'bzpt',idiag_bzpt)
        call parse_name(iname,cname(iname),cform(iname),'jxpt',idiag_jxpt)
        call parse_name(iname,cname(iname),cform(iname),'jypt',idiag_jypt)
        call parse_name(iname,cname(iname),cform(iname),'jzpt',idiag_jzpt)
        call parse_name(iname,cname(iname),cform(iname),'Expt',idiag_Expt)
        call parse_name(iname,cname(iname),cform(iname),'Eypt',idiag_Eypt)
        call parse_name(iname,cname(iname),cform(iname),'Ezpt',idiag_Ezpt)
        call parse_name(iname,cname(iname),cform(iname),'axpt',idiag_axpt)
        call parse_name(iname,cname(iname),cform(iname),'aypt',idiag_aypt)
        call parse_name(iname,cname(iname),cform(iname),'azpt',idiag_azpt)
        call parse_name(iname,cname(iname),cform(iname),'bxp2',idiag_bxp2)
        call parse_name(iname,cname(iname),cform(iname),'byp2',idiag_byp2)
        call parse_name(iname,cname(iname),cform(iname),'bzp2',idiag_bzp2)
        call parse_name(iname,cname(iname),cform(iname),'jxp2',idiag_jxp2)
        call parse_name(iname,cname(iname),cform(iname),'jyp2',idiag_jyp2)
        call parse_name(iname,cname(iname),cform(iname),'jzp2',idiag_jzp2)
        call parse_name(iname,cname(iname),cform(iname),'Exp2',idiag_Exp2)
        call parse_name(iname,cname(iname),cform(iname),'Eyp2',idiag_Eyp2)
        call parse_name(iname,cname(iname),cform(iname),'Ezp2',idiag_Ezp2)
        call parse_name(iname,cname(iname),cform(iname),'axp2',idiag_axp2)
        call parse_name(iname,cname(iname),cform(iname),'ayp2',idiag_ayp2)
        call parse_name(iname,cname(iname),cform(iname),'azp2',idiag_azp2)
        call parse_name(iname,cname(iname),cform(iname),'uxBrms',idiag_uxBrms)
        call parse_name(iname,cname(iname),cform(iname),'Bresrms',idiag_Bresrms)
        call parse_name(iname,cname(iname),cform(iname),'Rmrms',idiag_Rmrms)
        call parse_name(iname,cname(iname),cform(iname),'jfm',idiag_jfm)
        call parse_name(iname,cname(iname),cform(iname),'bmxy_rms',idiag_bmxy_rms)
        call parse_name(iname,cname(iname),cform(iname),'etasmagm',idiag_etasmagm)
        call parse_name(iname,cname(iname),cform(iname),'etasmagmin',idiag_etasmagmin)
        call parse_name(iname,cname(iname),cform(iname),'etasmagmax',idiag_etasmagmax)
        call parse_name(iname,cname(iname),cform(iname),'etavamax',idiag_etavamax)
        call parse_name(iname,cname(iname),cform(iname),'etajmax',idiag_etajmax)
        call parse_name(iname,cname(iname),cform(iname),'etaj2max',idiag_etaj2max)
        call parse_name(iname,cname(iname),cform(iname),'etajrhomax',idiag_etajrhomax)
        call parse_name(iname,cname(iname),cform(iname),'etaaniso',idiag_etaaniso)
        call parse_name(iname,cname(iname),cform(iname),'cosjbm',idiag_cosjbm)
        call parse_name(iname,cname(iname),cform(iname),'jparallelm',idiag_jparallelm)
        call parse_name(iname,cname(iname),cform(iname),'jperpm',idiag_jperpm)
        call parse_name(iname,cname(iname),cform(iname),'coshjbm',idiag_coshjbm)
        call parse_name(iname,cname(iname),cform(iname),'hjparallelm',idiag_hjparallelm)
        call parse_name(iname,cname(iname),cform(iname),'hjperpm',idiag_hjperpm)
      enddo
!
      if (idiag_exabot/=0) call set_type(idiag_exabot,lint=.true.)
      if (idiag_exatop/=0) call set_type(idiag_exatop,lint=.true.)

      if (idiag_b2tm/=0) then
        if (ibbt==0) call stop_it("Cannot calculate b2tm if ibbt==0")
        idiag_b2tm=0
      endif
      if (idiag_jbtm/=0) then
        if (ibbt==0) call stop_it("Cannot calculate jbtm if ibbt==0")
        idiag_jbtm=0
      endif
      if (idiag_bjtm/=0) then
        if (ijjt==0) call stop_it("Cannot calculate bjtm if ijjt==0")
        idiag_bjtm=0
      endif
!
!  Quantities which are averaged over half (north-south) the box.
!
      iname_half=name_half_max
!
!  Magnetic helicity (north and south) of total field.
!
      if ((idiag_abmn/=0).or.(idiag_abms/=0)) then
        iname_half=iname_half+1
        idiag_abmh=iname_half
      endif
!
!  Current helicity (north and south) of total field.
!
      if ((idiag_jbmn/=0).or.(idiag_jbms/=0)) then
        iname_half=iname_half+1
        idiag_jbmh=iname_half
      endif
!
!  Magnetic energy (north and south) of total field.
!
      if ((idiag_brmsn/=0).or.(idiag_brmss/=0)) then
        iname_half=iname_half+1
        idiag_brmsh=iname_half
      endif
!
!  Magnetic helicity (north and south) of mean field.
!
      if ((idiag_ambmzn/=0).or.(idiag_ambmzs/=0)) then
        iname_half=iname_half+1
        idiag_ambmzh=iname_half
      endif
!
!  Update name_half_max.
!
      name_half_max=iname_half
!
!  Currently need to force zaverage calculation at every lout step for
!  bmx and bmy and bmxy_rms.
!
      if ((idiag_bmx+idiag_bmy+idiag_bmxy_rms)>0) ldiagnos_need_zaverages=.true.
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'b2mx',idiag_b2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bxmx',idiag_bxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bymx',idiag_bymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bzmx',idiag_bzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bx2mx',idiag_bx2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'by2mx',idiag_by2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bz2mx',idiag_bz2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bxbymx',idiag_bxbymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bxbzmx',idiag_bxbzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bybzmx',idiag_bybzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'betamx',idiag_betamx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'beta2mx',idiag_beta2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'jxbrxmx',idiag_jxbrxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'jxbrymx',idiag_jxbrymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'jxbrzmx',idiag_jxbrzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'mflux_x',idiag_mflux_x)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'etatotalmx',idiag_etatotalmx)
     enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'bxmy',idiag_bxmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'bymy',idiag_bymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'bzmy',idiag_bzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'bx2my',idiag_bx2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'by2my',idiag_by2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'bz2my',idiag_bz2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'bxbymy',idiag_bxbymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'bxbzmy',idiag_bxbzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'bybzmy',idiag_bybzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'jxbrxmy',idiag_jxbrxmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'jxbrymy',idiag_jxbrymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'jxbrzmy',idiag_jxbrzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'mflux_y',idiag_mflux_y)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'axmz',idiag_axmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'aymz',idiag_aymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'azmz',idiag_azmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'abuxmz',idiag_abuxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'abuymz',idiag_abuymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'abuzmz',idiag_abuzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uabxmz',idiag_uabxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uabymz',idiag_uabymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uabzmz',idiag_uabzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bbxmz',idiag_bbxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bbymz',idiag_bbymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bbzmz',idiag_bbzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bxmz',idiag_bxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bymz',idiag_bymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bzmz',idiag_bzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'jxmz',idiag_jxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'jymz',idiag_jymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'jzmz',idiag_jzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Exmz',idiag_Exmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Eymz',idiag_Eymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Ezmz',idiag_Ezmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bx2mz',idiag_bx2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'by2mz',idiag_by2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bz2mz',idiag_bz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bx2rmz',idiag_bx2rmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'by2rmz',idiag_by2rmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bz2rmz',idiag_bz2rmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'beta1mz',idiag_beta1mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'betamz',idiag_betamz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'beta2mz',idiag_beta2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'bxbymz',idiag_bxbymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'bxbzmz',idiag_bxbzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'bybzmz',idiag_bybzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'b2mz',idiag_b2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bf2mz',idiag_bf2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'j2mz',idiag_j2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'poynzmz',idiag_poynzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'jxbrxmz',idiag_jxbrxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'jxbrymz',idiag_jxbrymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'jxbrzmz',idiag_jxbrzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'mflux_z',idiag_mflux_z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'jbmz',idiag_jbmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'d6abmz',idiag_d6abmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'d6amz1',idiag_d6amz1)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'d6amz2',idiag_d6amz2)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'d6amz3',idiag_d6amz3)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'abmz',idiag_abmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ubmz',idiag_ubmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uamz',idiag_uamz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bzdivamz',idiag_bzdivamz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'divamz',idiag_divamz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxbxmz',idiag_uxbxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uybxmz',idiag_uybxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzbxmz',idiag_uzbxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxbymz',idiag_uxbymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uybymz',idiag_uybymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzbymz',idiag_uzbymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxbzmz',idiag_uxbzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uybzmz',idiag_uybzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzbzmz',idiag_uzbzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'examz1',idiag_examz1)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'examz2',idiag_examz2)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'examz3',idiag_examz3)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'e3xamz1',idiag_e3xamz1)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'e3xamz2',idiag_e3xamz2)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'e3xamz3',idiag_e3xamz3)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'etatotalmz',idiag_etatotalmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'epsMmz',idiag_epsMmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Rmmz',idiag_Rmmz)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do ixz=1,nnamexz
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'b2mxz',idiag_b2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'axmxz',idiag_axmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'aymxz',idiag_aymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'azmxz',idiag_azmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bx1mxz',idiag_bx1mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'by1mxz',idiag_by1mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bz1mxz',idiag_bz1mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bxmxz',idiag_bxmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bymxz',idiag_bymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bzmxz',idiag_bzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'jxmxz',idiag_jxmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'jymxz',idiag_jymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'jzmxz',idiag_jzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bx2mxz',idiag_bx2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'by2mxz',idiag_by2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bz2mxz',idiag_bz2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bxbymxz',idiag_bxbymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bxbzmxz',idiag_bxbzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bybzmxz',idiag_bybzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uybxmxz',idiag_uybxmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uybzmxz',idiag_uybzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'Exmxz',idiag_Exmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'Eymxz',idiag_Eymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'Ezmxz',idiag_Ezmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'vAmxz',idiag_vAmxz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'axmxy',idiag_axmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'aymxy',idiag_aymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'azmxy',idiag_azmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxmxy',idiag_bxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bymxy',idiag_bymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bzmxy',idiag_bzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'jxmxy',idiag_jxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'jymxy',idiag_jymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'jzmxy',idiag_jzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bx2mxy',idiag_bx2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'by2mxy',idiag_by2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bz2mxy',idiag_bz2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxbymxy',idiag_bxbymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxbzmxy',idiag_bxbzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bybzmxy',idiag_bybzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'jbmxy',idiag_jbmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'abmxy',idiag_abmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ubmxy',idiag_ubmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'examxy1',idiag_examxy1)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'examxy2',idiag_examxy2)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'examxy3',idiag_examxy3)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Exmxy',idiag_Exmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Eymxy',idiag_Eymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Ezmxy',idiag_Ezmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'StokesImxy',idiag_StokesImxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'StokesQmxy',idiag_StokesQmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'StokesUmxy',idiag_StokesUmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'StokesQ1mxy',idiag_StokesQ1mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'StokesU1mxy',idiag_StokesU1mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'beta1mxy',idiag_beta1mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'poynxmxy',idiag_poynxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'poynymxy',idiag_poynymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'poynzmxy',idiag_poynzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'dbxdxmxy',idiag_dbxdxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'dbxdymxy',idiag_dbxdymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'dbxdzmxy',idiag_dbxdzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'dbydxmxy',idiag_dbydxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'dbydymxy',idiag_dbydymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'dbydzmxy',idiag_dbydzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'dbzdxmxy',idiag_dbzdxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'dbzdymxy',idiag_dbzdymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'dbzdzmxy',idiag_dbzdzmxy)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'brmphi'  ,idiag_brmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'brsphmphi',idiag_brsphmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'bthmphi',idiag_bthmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'bpmphi'  ,idiag_bpmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'bzmphi'  ,idiag_bzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'b2mphi'  ,idiag_b2mphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'jbmphi'  ,idiag_jbmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uxbrmphi',idiag_uxbrmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uxbpmphi',idiag_uxbpmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uxbzmphi',idiag_uxbzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'jxbrmphi',idiag_jxbrmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'jxbpmphi',idiag_jxbpmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'jxbzmphi',idiag_jxbzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'armphi'  ,idiag_armphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'apmphi'  ,idiag_apmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'azmphi'  ,idiag_azmphi)
!
      enddo
!
!  Check for those quantities for which we want phiz-averages.
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'brmr',  idiag_brmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'bpmr',  idiag_bpmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'bzmr',  idiag_bzmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'armr',  idiag_armr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'apmr',  idiag_apmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'azmr',  idiag_azmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'b2mr',  idiag_b2mr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'brbpmr',idiag_brbpmr)
      enddo
!
!  Check for those quantities for which we want to have in the sound.dat file.
!
      do iname_sound=1,nname_sound
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'bxpt',idiag_bxpt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'bypt',idiag_bypt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'bzpt',idiag_bzpt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'axpt',idiag_axpt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'aypt',idiag_aypt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'azpt',idiag_azpt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'jxpt',idiag_jxpt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'jypt',idiag_jypt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'jzpt',idiag_jzpt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'Expt',idiag_Expt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'Eypt',idiag_Eypt)
        call parse_name(iname_sound,cname_sound(iname_sound),cform_sound(iname_sound),'Ezpt',idiag_Ezpt)
      enddo
!
!  check for those quantities for which we want video slices
!
      idum=0
      do inamev=1,nnamev
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'aa',idum)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'aps',ivid_aps)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'bb',ivid_bb)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'jj',ivid_jj)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'b2',ivid_b2)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'j2',ivid_j2)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'bb_sph',ivid_bb_sph)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'ab',ivid_ab)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'jb',ivid_jb)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'beta1',ivid_beta1)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'poynting',ivid_poynting)
      enddo
!
!  call corresponding mean-field routine
!
      if (lmagn_mf) call rprint_magn_mf(lreset,lwrite)
!
    endsubroutine rprint_magnetic
!***********************************************************************
    subroutine dynamical_resistivity(uc)
!
!  Dynamically set resistivity coefficient given fixed mesh Reynolds number.
!
!  27-jul-11/ccyang: coded
!
!  Input Argument
!      uc
!          Characteristic velocity of the system.
!
      real, intent(in) :: uc
!
!  Hyper-resistivity coefficient
!
      if (eta_hyper3 /= 0.0) eta_hyper3 = pi5_1 * uc * dxmax**5 / re_mesh
      if (eta_hyper3_mesh /= 0.0) eta_hyper3_mesh = pi5_1 * uc / re_mesh / sqrt(real(dimensionality))
!
    endsubroutine dynamical_resistivity
!***********************************************************************
    subroutine split_update_magnetic(f)
!
!  Update the magnetic potential by integrating the operator split
!  magnetic terms.
!
!  23-aug-13/ccyang: coded
!
      use ImplicitDiffusion, only: integrate_diffusion
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
!  Implicitly solve the resistive term.
!
      if (limplicit_resistivity) call integrate_diffusion(get_resistivity_implicit, f, iax, iaz)
!
    endsubroutine split_update_magnetic
!***********************************************************************
   subroutine magnetic_after_timestep(f,df,dtsub)
!
      use Mpicomm , only: mpibcast_real
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
      if (lfargo_advection) then
        if (lkeplerian_gauge) call keplerian_gauge(f)
        if (lrmv.and.lremove_volume_average) call remove_volume_average(f)
      endif
!
      if (lresi_vAspeed) then
        if (lroot) then
          if (fname(idiag_vArms)/=0.0) vArms=fname(idiag_vArms)
        endif
        call mpibcast_real(vArms)
      endif
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine magnetic_after_timestep
!***********************************************************************
    subroutine keplerian_gauge(f)
!
      use Mpicomm , only: mpiallreduce_sum
      use Deriv, only: der
      use Boundcond, only: update_ghosts
!
!  Substract mean emf from the radial component of the induction
!  equation. Activated only when large Bz fields and are present
!  keplerian advection. Due to this u_phi x Bz term, the radial
!  component of the magnetic potential
!  develops a divergence that grows linearly in time. Since it is
!  purely divergent, it is okay analytically. But numerically it leads to
!  problems if this divergent grows bigger than the curl, which it does
!  eventually.
!
!  This is a cylindrical version of the rtime_phiavg special file.
!
!  13-sep-07/wlad: adapted from remove_mean_momenta
!  28-mar-17/MR: reinstated update_ghosts.
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mx,mz) :: fsum_tmp,glambda_rz
      real, dimension (mx,my,mz) :: lambda
      real, dimension (nx) :: glambda_z
      real :: fac
      integer :: i
!
      !if (.not.lupdate_bounds_before_special) then
      !  print*,'The boundaries have not been updated prior '
      !  print*,'to calling this subroutine. This may lead '
      !  print*,'to troubles since it needs derivatives '
      !  print*,'and integrals, thus properly set ghost zones. '
      !  print*,'Use lupdate_bounds_before_special=T in '
      !  print*,'the run_pars of run.in.'
      !  call fatal_error("apply_keplerian_gauge","")
      !endif
!
      fac = 1.0/nygrid
!
! Set ghost zones of iax.
!
      call update_ghosts(f,iax)
!
! Average over phi - the result is a (mr=mx,mz) array
!
      fsum_tmp = 0.
      do m=m1,m2; do n=1,mz
        fsum_tmp(:,n) = fsum_tmp(:,n) + fac*f(:,m,n,iax)
      enddo; enddo
!
! The sum has to be done processor-wise
! Sum over processors of same ipz, and different ipy
!
      call mpiallreduce_sum(fsum_tmp,glambda_rz,(/mx,mz/),idir=2)
!
! Gauge-transform radial A
!
      do m=m1,m2
        f(l1:l2,m,n1:n2,iax) = f(l1:l2,m,n1:n2,iax) - glambda_rz(l1:l2,n1:n2)
      enddo
!
! Integrate in R to get lambda, using N=6 composite Simpson's rule.
! Ghost zones in r needed for glambda_r.
!
      do i=l1,l2 ; do n=1,mz
        lambda(i,:,n) = dx/6.*(   glambda_rz(i-3,n)+glambda_rz(i+3,n)+&
                               4*(glambda_rz(i-2,n)+glambda_rz(i  ,n)+glambda_rz(i+2,n))+&
                               2*(glambda_rz(i-1,n)+glambda_rz(i+1,n)))
      enddo; enddo
!
!  Gauge-transform vertical A. Ghost zones in z needed for lambda.
!
      do m=m1,m2; do n=n1,n2
        call der(lambda,glambda_z,3)
        f(l1:l2,m,n,iaz) = f(l1:l2,m,n,iaz) - glambda_z
      enddo; enddo
!
    endsubroutine keplerian_gauge
!********************************************************************
    subroutine remove_volume_average(f)
!
      use Mpicomm , only: mpiallreduce_sum
!
!  Substract mean emf from the radial component of the induction
!  equation. Activated only when large Bz fields and are present
!  keplerian advection. Due to this u_phi x Bz term, the radial
!  component of the magnetic potential
!  develops a divergence that grows linearly in time. Since it is
!  purely divergent, it is okay analytically. But numerically it leads to
!  problems if this divergent grows bigger than the curl, which it does
!  eventually.
!
!  This is a cylindrical version of the rtime_phiavg special file.
!
!  13-sep-07/wlad: adapted from remove_mean_momenta
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real :: fsum_tmp,mean_ax
      integer :: i
!
! Average over phi - the result is a (nr,nz) array
!
      fsum_tmp = 0.
      do m=m1,m2; do n=n1,n2 ; do i=l1,l2
        fsum_tmp = fsum_tmp + f(i,m,n,iax)
      enddo; enddo; enddo
      fsum_tmp = fsum_tmp/nwgrid
!
! The sum has to be done processor-wise
! Sum over processors of same ipz, and different ipy
!
      call mpiallreduce_sum(fsum_tmp,mean_ax)
!
! Gauge-transform radial A
!
      f(l1:l2,m1:m2,n1:n2,iax) = f(l1:l2,m1:m2,n1:n2,iax) - mean_ax
!
    endsubroutine remove_volume_average
!********************************************************************
    subroutine get_resistivity_implicit(ndc, diffus_coeff, iz)
!
!  Gets the diffusion coefficient along a given pencil for the implicit algorithm.
!
!  23-aug-13/ccyang: coded.
!
      integer, intent(in) :: ndc
      real, dimension(ndc), intent(out) :: diffus_coeff
      integer, intent(in), optional :: iz
!
!  Uniform resistivity
!
      if (lresi_eta_const) then
        diffus_coeff = eta
      else
        diffus_coeff = 0.0
      endif
!
!  z-dependent resistivity
!
      zdep: if (lresi_zdep) then
        if (present(iz)) then
          diffus_coeff = diffus_coeff + eta_zgrid(iz)
        else
          diffus_coeff = diffus_coeff + eta_zgrid
        endif
      endif zdep
!
    endsubroutine get_resistivity_implicit
!***********************************************************************
    subroutine expand_shands_magnetic
!
!  Expands shorthand labels of magnetic diagnostics.
!
!  16-may-12/MR: coded
!
      use Diagnostics, only : name_is_present, expand_cname
!
      if (nnamerz>0) then
!
        call expand_cname(cnamerz,nnamerz,name_is_present(cnamerz,'bbmphi'),&
                          'brmphi','bpmphi','bzmphi')
        if (name_is_present(cnamerz,'bpmphi')>0) then
          call expand_cname(cnamerz,nnamerz,name_is_present(cnamerz,'bbsphmphi'),&
                            'brsphmphi','bthmphi')
        else
          call expand_cname(cnamerz,nnamerz,name_is_present(cnamerz,'bbsphmphi'),&
                            'brsphmphi','bthmphi','bpmphi')
        endif
!
        call expand_cname(cnamerz,nnamerz,name_is_present(cnamerz,'uxbmphi'),&
                          'uxbrmphi','uxbpmphi','uxbzmphi')
        call expand_cname(cnamerz,nnamerz,name_is_present(cnamerz,'jxbmphi'),&
                          'jxbrmphi','jxbpmphi','jxbzmphi')
!
      endif
!
    endsubroutine expand_shands_magnetic
!***********************************************************************
    subroutine get_bext(B_ext_out)
!
!  Get the external magnetic field at current time step.
!
!  lbext_curvilinear = .true. is default.  The B_ext the user defines in 
!  magnetic_init_pars respects the coordinate system of preference which 
!  means that B_ext=(0.0,1.0,0.0) is an azimuthal field in cylindrical 
!  coordinates and a polar one in spherical.
!
!  01-nov-14/ccyang: coded
!
      real, dimension(3), intent(out) :: B_ext_out
!
      real :: c, s
!
      addBext: if (any(B_ext /= 0.0)) then
!
        precess: if (omega_Bz_ext /= 0.0) then
!
!  Allow external field to precess about z-axis with frequency omega_Bz_ext.
!
          coord1: if (lcartesian_coords .or. lbext_curvilinear) then
            c = cos(omega_Bz_ext * t)
            s = sin(omega_Bz_ext * t)
            B_ext_out(1) = B_ext(1) * c - B_ext(2) * s
            B_ext_out(2) = B_ext(1) * s + B_ext(2) * c
            B_ext_out(3) = B_ext(3)
          else coord1
            call fatal_error('get_bext', 'precession of the external field not implemented for curvilinear coordinates')
          endif coord1
        else precess
!
!  Or add uniform background field.
!
          coord2: if (lcartesian_coords .or. lbext_curvilinear) then
            B_ext_out = B_ext
          elseif (lcylindrical_coords) then coord2
            B_ext_out(1) =  B_ext(1) * cos(y(m)) + B_ext(2) * sin(y(m))
            B_ext_out(2) = -B_ext(1) * sin(y(m)) + B_ext(2) * cos(y(m))
            B_ext_out(3) =  B_ext(3)
          elseif (lspherical_coords) then coord2
            B_ext_out(1) =  B_ext(1) * sinth(m) * cos(z(n)) + B_ext(2) * sinth(m) * sin(z(n)) + B_ext(3) * costh(m)
            B_ext_out(2) =  B_ext(1) * costh(m) * cos(z(n)) + B_ext(2) * costh(m) * sin(z(n)) - B_ext(3) * sinth(m)
            B_ext_out(3) = -B_ext(1)            * sin(z(n)) + B_ext(2)            * cos(z(n))
          endif coord2
        endif precess
      else addBext
!
!  Or no background field.
!
        B_ext_out = 0.0
!
      endif addBext
!
!  Make the field gently increasing.
!
      gentle: if (t_bext > 0.0 .and. t < t_bext) then
        if (t <= t0_bext) then
          B_ext_out = B0_ext
        else
          B_ext_out = B0_ext + 0.5 * (1.0 - cos(pi * (t - t0_bext) / (t_bext - t0_bext))) * (B_ext_out - B0_ext)
        endif
      endif gentle
!
    endsubroutine get_bext
!***********************************************************************
    subroutine pushpars2c(p_par)

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr_c(eta,p_par(1))

    endsubroutine pushpars2c
!***********************************************************************
    subroutine pushdiags2c(p_diag)

    use Diagnostics, only: set_type

    integer, parameter :: n_diags=8
    integer(KIND=ikind8), dimension(n_diags) :: p_diag

    call copy_addr_c(idiag_brms,p_diag(1))
    call set_type(idiag_brms,lsqrt=.true.)
    call copy_addr_c(idiag_bmax,p_diag(2))
    call set_type(idiag_bmax,lmax=.true.)
    call copy_addr_c(idiag_bxmin,p_diag(3))
    call set_type(idiag_bxmin,lmin=.true.)
    call copy_addr_c(idiag_bymin,p_diag(4))
    call set_type(idiag_bymin,lmin=.true.)
    call copy_addr_c(idiag_bzmin,p_diag(5))
    call set_type(idiag_bzmin,lmin=.true.)
    call copy_addr_c(idiag_bxmax,p_diag(6))
    call set_type(idiag_bxmax,lmax=.true.)
    call copy_addr_c(idiag_bymax,p_diag(7))
    call set_type(idiag_bymax,lmax=.true.)
    call copy_addr_c(idiag_bzmax,p_diag(8))
    call set_type(idiag_bzmax,lmax=.true.)

    endsubroutine pushdiags2c
!***********************************************************************
endmodule Magnetic
