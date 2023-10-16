! $Id$
!
!  This module takes care of the continuity equation.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .true.
! CPARAM logical, parameter :: lanelastic = .false.
! CPARAM logical, parameter :: lboussinesq = .false.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED lnrho; rho; rho1; glnrho(3); grho(3); uglnrho; ugrho
! PENCILS PROVIDED glnrho2; del2lnrho; del2rho; del6lnrho; del6rho
! PENCILS PROVIDED hlnrho(3,3); sglnrho(3); uij5glnrho(3); transprho
! PENCILS PROVIDED ekin, uuadvec_glnrho; uuadvec_grho
! PENCILS PROVIDED rhos1; glnrhos(3)
! PENCILS PROVIDED totenergy_rel; divss
!
!***************************************************************
module Density
!
  use Cdata
  use General, only: keep_compiler_quiet, itoa
  use Messages
  use EquationOfState, only: cs0, cs20, cs2bot, cs2top, rho0, lnrho0
  use DensityMethods
!
  implicit none
!
  include 'density.h'
!
  real, dimension (ninit) :: ampllnrho=0.0, widthlnrho=0.1
  real, dimension (ninit) :: rho_left=1.0, rho_right=1.0
  real, dimension (ninit) :: amplrho=0.0, phase_lnrho=0.0, radius_lnrho=0.5
  real, dimension (ninit) :: xblob=0.0, yblob=0.0, zblob=0.0
  real, dimension (ninit) :: kx_lnrho=1.0, ky_lnrho=1.0, kz_lnrho=1.0
  real, dimension (ninit) :: kxx_lnrho=0.0, kyy_lnrho=0.0, kzz_lnrho=0.0
  real, dimension (mz,1) :: lnrhomz
  real, dimension (nz) :: lnrho_init_z_nz=0.0
  real, dimension (mz) :: lnrho_init_z=0.0
  real, dimension (mz) :: dlnrhodz_init_z=0.0, del2lnrho_glnrho2_init_z=0.0
  real, dimension (3) :: diffrho_hyper3_aniso=0.0
  real, dimension (nx) :: profx_ffree=1.0, dprofx_ffree=0.0
  real, dimension (my) :: profy_ffree=1.0, dprofy_ffree=0.0
  real, dimension (mz) :: profz_ffree=1.0, dprofz_ffree=0.0
  real, dimension(mz) :: profz_eos=1.0,dprofz_eos=0.0
  real, target :: mpoly=impossible
  real, pointer :: mpoly0, mpoly1, mpoly2, eps_hless
  real, dimension(nx) :: xmask_den
  real, dimension(nx) :: fprofile_x=1.
  real, dimension(nz) :: fprofile_z=1.
  real, dimension(nz) :: zmask_den
  real, dimension(nx) :: reduce_cs2_profx = 1.0
  real, dimension(mz) :: reduce_cs2_profz = 1.0
  real :: width_eos_prof=0.2
  character(LEN=labellen) :: ireference_state='nothing', ieos_profile='nothing'
  real :: reference_state_mass=0.
!
!  Schur flow quantities
!
  real, dimension(1) :: Schur_dlnrho_RHS_xyzaver
  real, dimension (nz) :: Schur_dlnrho_RHS_xyaver_z
  real, dimension (nx,ny) :: Schur_dlnrho_RHS_zaver_xy
!
! reference state, components:  1       2          3              4            5      6     7         8            9
!                              rho, d rho/d z, d^2 rho/d z^2, d^6 rho/d z^6, d p/d z, s, d s/d z, d^2 s/d z^2, d^6 s/d z^6
  real, dimension(nx,nref_vars) :: reference_state=0.
  real, dimension(2) :: density_xaver_range=(/-max_real,max_real/)
  real, dimension(2) :: density_zaver_range=(/-max_real,max_real/)
  real :: lnrho_const=0.0, rho_const=1.0, Hrho=1., ggamma=impossible
  real :: cdiffrho=0.0, diffrho=0.0, diff_cspeed=0.5
  real :: diffrho_hyper3=0.0, diffrho_hyper3_mesh=5.0, diffrho_shock=0.0
  real :: eps_planet=0.5, q_ell=5.0, hh0=0.0
  real :: mass_source_omega=0.
  real :: co1_ss=0.0, co2_ss=0.0, Sigma1=150.0
  real :: lnrho_int=0.0, lnrho_ext=0.0, damplnrho_int=0.0, damplnrho_ext=0.0
  real :: wdamp=0.0, density_floor=-1.0, density_floor_exp=0.
  real :: mass_source_Mdot=0.0, mass_source_sigma=0.0, mass_source_offset=0.0
  real :: tstart_mass_source=0.0, tstop_mass_source=-1.0
  real :: lnrho_z_shift=0.0
  real :: powerlr=3.0, zoverh=1.5, hoverr=0.05
  real :: rzero_ffree=0.,wffree=0.
  real :: rho_top=1.,rho_bottom=1.
  real :: rmax_mass_source=0., fnorm=1.
  real :: r0_rho=impossible
  real :: invgrav_ampl = 0.0
  real :: rnoise_int=impossible,rnoise_ext=impossible
  real :: mass_source_tau1=1.
  real :: mass_cloud=0.0, T_cloud=0.0, T_cloud_out_rel=1.0, xi_coeff=1.0
  real :: temp_coeff=1.0, dens_coeff=1.0, temp_trans = 0.0
  real :: temp_coeff_out = 1.0
  real :: rss_coef1=1.0, rss_coef2=1.0
  real :: total_mass=-1.
  real :: rescale_rho=1.0
  real :: xjump_mid=0.0,yjump_mid=0.0,zjump_mid=0.0
  real, target :: reduce_cs2 = 1.0
  complex :: coeflnrho=0.0
  integer, parameter :: ndiff_max=4
  integer :: iglobal_gg=0
  logical :: lrelativistic_eos=.false., ladvection_density=.true.
  logical, pointer :: lconservative, lhiggsless
  logical :: lisothermal_fixed_Hrho=.false.
  logical :: lmass_source=.false., lmass_source_random=.false., lcontinuity_gas=.true.
  logical :: lupw_lnrho=.false.,lupw_rho=.false.
  logical :: ldiff_normal=.false.,ldiff_hyper3=.false.,ldiff_shock=.false.
  logical :: ldiff_cspeed=.false.
  logical :: ldiff_hyper3lnrho=.false.,ldiff_hyper3_aniso=.false.
  logical :: ldiff_hyper3_polar=.false.,lanti_shockdiffusion=.false.
  logical :: ldiff_hyper3_mesh=.false.,ldiff_hyper3_strict=.false.
  logical :: ldiff_hyper3lnrho_strict=.false.
  logical :: lfreeze_lnrhoint=.false.,lfreeze_lnrhoext=.false.
  logical :: lfreeze_lnrhosqu=.false.
  logical :: lrho_as_aux=.false., ldiffusion_nolog=.false.
  logical :: lrho_flucz_as_aux=.false.
  logical :: lmassdiff_fix = .true.,lmassdiff_fixmom = .false.,lmassdiff_fixkin = .false.
  logical :: lcheck_negative_density=.false.
  logical :: lcalc_lnrhomean=.false.
  logical :: ldensity_profile_masscons=.false.
  logical :: lffree=.false.
  logical :: lSchur_3D3D1D=.false.
  logical, target :: lreduced_sound_speed=.false.
  logical, target :: lscale_to_cs2top=.false.
  logical :: lconserve_total_mass=.false.
  real :: density_ceiling=-1.
  logical :: lreinitialize_lnrho=.false., lreinitialize_rho=.false.
  logical :: lsubtract_init_stratification=.false., lwrite_stratification=.false.
  character (len=labellen), dimension(ninit) :: initlnrho='nothing'
  character (len=labellen) :: strati_type='lnrho_ss'
  character (len=labellen), dimension(ndiff_max) :: idiff=''
  character (len=labellen) :: borderlnrho='nothing'
  character (len=labellen) :: mass_source_profile='nothing'
  character (len=intlen) :: iinit_str
  character (len=labellen) :: div_sld_dens='2nd'
  character (len=labellen) :: ffree_profile='none'
  character (len=fnlen) :: datafile='dens_temp.dat'
  character (len=labellen) :: cloud_mode='isothermal'
  character (len=labellen) :: density_floor_profile='uniform'
  logical :: ldensity_slope_limited=.false.
  real :: h_sld_dens=2.0, nlf_sld_dens=1.0
  real, dimension(3) :: beta_glnrho_global = 0., beta_glnrho_scaled=0.
!
  namelist /density_init_pars/ &
      ampllnrho, initlnrho, widthlnrho, rho_left, rho_right, lnrho_const, &
      Hrho, rho_const, cs2bot, cs2top, radius_lnrho, eps_planet, xblob, &
      yblob, zblob, b_ell, q_ell, hh0, rbound, lwrite_stratification, &
      mpoly, ggamma, &
      strati_type, beta_glnrho_global, kx_lnrho, ky_lnrho, kz_lnrho, &
      amplrho, phase_lnrho, coeflnrho, kxx_lnrho, kyy_lnrho,  kzz_lnrho, &
      co1_ss, co2_ss, Sigma1, idiff, ldensity_nolog, wdamp, lcontinuity_gas, &
      lisothermal_fixed_Hrho, density_floor, lanti_shockdiffusion, &
      density_floor_profile, density_floor_exp, &
      lmassdiff_fix, lmassdiff_fixmom, lmassdiff_fixkin,  &
      lrho_as_aux, ldiffusion_nolog, lnrho_z_shift, powerlr, zoverh, hoverr, &
      lffree, ffree_profile, rzero_ffree, wffree, rho_top, rho_bottom, &
      r0_rho, invgrav_ampl, rnoise_int, rnoise_ext, datafile, mass_cloud, &
      T_cloud, cloud_mode, T_cloud_out_rel, xi_coeff, density_xaver_range, &
      dens_coeff, temp_coeff, temp_trans, temp_coeff_out, reduce_cs2, &
      lreduced_sound_speed, lrelativistic_eos, &
      lscale_to_cs2top, density_zaver_range, &
      ieos_profile, width_eos_prof, &
      lconserve_total_mass, total_mass, ireference_state, lrho_flucz_as_aux,&
      ldensity_linearstart, xjump_mid, yjump_mid, zjump_mid
!
  namelist /density_run_pars/ &
      cdiffrho, diffrho, diffrho_hyper3, diffrho_hyper3_mesh, diffrho_shock, &
      cs2bot, cs2top, lupw_lnrho, lupw_rho, idiff, &
      lmass_source, lmass_source_random, diff_cspeed, &
      mass_source_profile, mass_source_Mdot, mass_source_sigma, &
      mass_source_offset, rmax_mass_source, lnrho_int, lnrho_ext, &
      damplnrho_int, damplnrho_ext, wdamp, lfreeze_lnrhoint, lfreeze_lnrhoext, &
      lnrho_const, lcontinuity_gas, borderlnrho, diffrho_hyper3_aniso, &
      lfreeze_lnrhosqu, density_floor, lanti_shockdiffusion, lrho_as_aux, &
      density_floor_profile, density_floor_exp, &
      ldiffusion_nolog, lcheck_negative_density, &
      lmassdiff_fix, lmassdiff_fixmom, lmassdiff_fixkin,&
      lcalc_lnrhomean, ldensity_profile_masscons, lffree, ffree_profile, &
      rzero_ffree, wffree, tstart_mass_source, tstop_mass_source, &
      density_xaver_range, mass_source_tau1, reduce_cs2, &
      lreduced_sound_speed, lrelativistic_eos, ladvection_density, &
      xblob, yblob, zblob, mass_source_omega, lscale_to_cs2top, &
      density_zaver_range, rss_coef1, rss_coef2, &
      ieos_profile, width_eos_prof, beta_glnrho_global, &
      lconserve_total_mass, total_mass, density_ceiling, &
      lreinitialize_lnrho, lreinitialize_rho, initlnrho, rescale_rho, &
      lsubtract_init_stratification, ireference_state, &
      h_sld_dens, lrho_flucz_as_aux, nlf_sld_dens, div_sld_dens, &
      lSchur_3D3D1D
!
!  Diagnostic variables (need to be consistent with reset list below).
!  Note: drho2m is based on rho0, while rhof2m is based on <rho>(z).
!
  integer :: idiag_rhom=0       ! DIAG_DOC: $\left<\varrho\right>$
                                ! DIAG_DOC:   \quad(mean density)
  integer :: idiag_rhomxmask=0  ! DIAG_DOC: $\left<\varrho\right>$ for
                                ! DIAG_DOC: the density_xaver_range
  integer :: idiag_rhomzmask=0  ! DIAG_DOC: $\left<\varrho\right>$ for
                                ! DIAG_DOC: the density_zaver_range
  integer :: idiag_rho2m=0      ! DIAG_DOC: $\left<\varrho^2\right>$
  integer :: idiag_rho4m=0      ! DIAG_DOC: $\left<\varrho^4\right>$
  integer :: idiag_rho6m=0      ! DIAG_DOC: $\left<\varrho^6\right>$
  integer :: idiag_rho12m=0     ! DIAG_DOC: $\left<\varrho^{12}\right>$
  integer :: idiag_rhof2m=0     ! DIAG_DOC: $\left<\varrho'^2\right>$
  integer :: idiag_lnrho2m=0    ! DIAG_DOC:
  integer :: idiag_drho2m=0     ! DIAG_DOC: $<(\varrho-\varrho_0)^2>$
  integer :: idiag_drhom=0      ! DIAG_DOC: $<\varrho-\varrho_0>$
  integer :: idiag_rhomin=0     ! DIAG_DOC: $\min(\rho)$
  integer :: idiag_rhomax=0     ! DIAG_DOC: $\max(\rho)$
  integer :: idiag_lnrhomin=0   ! DIAG_DOC: $\min(\log\rho)$
  integer :: idiag_lnrhomax=0   ! DIAG_DOC: $\max(\log\rho)$
  integer :: idiag_rhorms=0     ! DIAG_DOC: $\sqrt{<\varrho^2>}$
  integer :: idiag_lnrhorms=0   ! DIAG_DOC: $\sqrt{<(\ln\varrho)^2>}$
  integer :: idiag_ugrhom=0     ! DIAG_DOC: $\left<\uv\cdot\nabla\varrho\right>$
  integer :: idiag_uglnrhom=0   ! DIAG_DOC: $\left<\uv\cdot\nabla\ln\varrho\right>$
  integer :: idiag_lnrhomphi=0  ! PHIAVG_DOC: $\left<\ln\varrho\right>_\varphi$
  integer :: idiag_rhomphi=0    ! PHIAVG_DOC: $\left<\varrho\right>_\varphi$
  integer :: idiag_dtd=0        ! DIAG_DOC:
  integer :: idiag_dtd3=0       ! DIAG_DOC:
  integer :: idiag_rhomr=0      ! DIAG_DOC:
  integer :: idiag_totmass=0    ! DIAG_DOC: $\int\varrho\,dV$
  integer :: idiag_mass=0       ! DIAG_DOC: $\int\varrho\,dV$
  integer :: idiag_sphmass=0    ! DIAG_DOC: $\int\varrho\,dV$ inside $r < r_{\rm diag}$.
  integer :: idiag_inertiaxx=0  ! DIAG_DOC: $xx$ component of the inertia tensor (spherical coordinates)
  integer :: idiag_inertiayy=0  ! DIAG_DOC: $yy$ component of the inertia tensor (spherical coordinates)
  integer :: idiag_inertiazz=0  ! DIAG_DOC: $zz$ component of the inertia tensor (spherical coordinates)
  integer :: idiag_inertiaxx_car=0  ! DIAG_DOC: $xx$ component of the inertia tensor (Cartesian coordinates)
  integer :: idiag_inertiayy_car=0  ! DIAG_DOC: $xx$ component of the inertia tensor (Cartesian coordinates)
  integer :: idiag_inertiazz_car=0  ! DIAG_DOC: $xx$ component of the inertia tensor (Cartesian coordinates)
  integer :: idiag_vol=0        ! DIAG_DOC: $\int\,dV$ (volume)
  integer :: idiag_grhomax=0    ! DIAG_DOC: $\max (|\nabla \varrho|)$
!
! xy averaged diagnostics given in xyaver.in
  integer :: idiag_rhomz=0      ! XYAVG_DOC: $\left<\varrho\right>_{xy}$
  integer :: idiag_rhoupmz=0    ! XYAVG_DOC: $\left<\varrho_\uparrow\right>_{xy}$
  integer :: idiag_rhodownmz=0  ! XYAVG_DOC: $\left<\varrho_\downarrow\right>_{xy}$
  integer :: idiag_rho2mz=0     ! XYAVG_DOC: $\left<\varrho^2\right>_{xy}$
  integer :: idiag_rho2upmz=0   ! XYAVG_DOC: $\left<\varrho_\uparrow^2\right>_{xy}$
  integer :: idiag_rho2downmz=0 ! XYAVG_DOC: $\left<\varrho_\downarrow^2\right>_{xy}$
  integer :: idiag_rhof2mz=0    ! XYAVG_DOC: $\left<\varrho'^2\right>_{xy}$
  integer :: idiag_rhof2upmz=0  ! XYAVG_DOC: $\left<\varrho'^2_\uparrow\right>_{xy}$
  integer :: idiag_rhof2downmz=0 ! XYAVG_DOC: $\left<\varrho'^2_\downarrow\right>_{xy}$
  integer :: idiag_gzlnrhomz=0  ! XYAVG_DOC: $\left<\nabla_z\ln\varrho\right>_{xy}$
  integer :: idiag_uglnrhomz=0  ! XYAVG_DOC: $\left<\uv\cdot\nabla\ln\varrho\right>_{xy}$
  integer :: idiag_ugrhomz=0    ! XYAVG_DOC: $\left<\uv\cdot\nabla\varrho\right>_{xy}$
  integer :: idiag_uygzlnrhomz=0! XYAVG_DOC: $\left<u_y\nabla_z\ln\varrho\right>_{xy}$
  integer :: idiag_uzgylnrhomz=0! XYAVG_DOC: $\left<u_z\nabla_y\ln\varrho\right>_{xy}$
!
! xz averaged diagnostics given in xzaver.in
  integer :: idiag_rhomy=0      ! XZAVG_DOC: $\left<\varrho\right>_{xz}$
!
! yz averaged diagnostics given in yzaver.in
  integer :: idiag_rhomx=0      ! YZAVG_DOC: $\left<\varrho\right>_{yz}$
  integer :: idiag_rho2mx=0     ! XYAVG_DOC: $\left<\varrho^2\right>_{yz}$
!
! y averaged diagnostics given in yaver.in
  integer :: idiag_rhomxz=0     ! YAVG_DOC: $\left<\varrho\right>_{y}$
!
! z averaged diagnostics given in zaver.in
  integer :: idiag_rhomxy=0     ! ZAVG_DOC: $\left<\varrho\right>_{z}$
  integer :: idiag_rho2mxy=0    ! ZAVG_DOC: $\left<\varrho^2\right>_{z}$
  integer :: idiag_sigma=0      ! ZAVG_DOC; $\Sigma\equiv\int\varrho\,\mathrm{d}z$
!
  interface calc_pencils_density
    module procedure calc_pencils_density_pnc
    module procedure calc_pencils_density_std
  endinterface calc_pencils_density
!
  interface calc_pencils_linear_density
    module procedure calc_pencils_linear_density_pnc
    module procedure calc_pencils_linear_density_std
  endinterface calc_pencils_linear_density
!
!  module auxiliaries
!
  logical :: lupdate_mass_source
  real, dimension(nx) :: diffus_diffrho
  real, dimension(nx) :: diffus_diffrho3
  real :: density_floor_log, density_ceiling_log
  real :: gamma, gamma1, gamma_m1, cp1
!
  contains
!***********************************************************************
    subroutine register_density
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrho; increase nvar accordingly.
!
!   4-jun-02/axel: adapted from hydro
!   21-oct-15/MR: changes for slope-limited diffusion
!  03-apr-20/joern: restructured and fixed slope-limited diffusion
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
      if (ldensity_nolog) then
        call farray_register_pde('rho',irho)
        ilnrho=irho
      else
        call farray_register_pde('lnrho',ilnrho)
      endif

      if (any(idiff=='density-slope-limited')) then
        if (dimensionality<3) lisotropic_advection=.true.
        lslope_limit_diff = .true.
        if (isld_char == 0) then
          call farray_register_auxiliary('sld_char',isld_char,communicated=.true.)
          if (lroot) write(15,*) 'sld_char = fltarr(mx,my,mz)*one'
          aux_var(aux_count)=',sld_char'
          if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=trim(aux_var(aux_count))//' $'
          aux_count=aux_count+1
        endif
      endif
!
!  Fluctuating density = \rho - \mean_xy(\rho)
!
      if (lrho_flucz_as_aux) then
        if (irho_flucz==0) then
          call farray_register_auxiliary('rho_flucz',irho_flucz)
        else
          if (lroot) print*, 'register_energy: irho_run_aver = ', irho_flucz
          call farray_index_append('irho_flucz',irho_flucz)
        endif
        if (lroot) write(15,*) 'rho_flucz = fltarr(mx,my,mz)*one'
        aux_var(aux_count)=',rho_flucz'
        if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=trim(aux_var(aux_count))//' $'
        aux_count=aux_count+1
      endif
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
!
! mpoly needs to be put here as in initialize_density it were
! too late for other modules (especially the nomodule noentropy).
!
      call put_shared_variable('mpoly',mpoly,caller='register_density')
!
!  Communicate lrelativistic_eos to magnetic too.
!
      call put_shared_variable('lrelativistic_eos',lrelativistic_eos)
!
!  Communicate lffree to entropy too.
!
      call put_shared_variable('lffree',lffree)
!
!  Put the force-free profiles to entropy and hydro too.
!
      if (lffree) then
        call put_shared_variable('profx_ffree',profx_ffree)
        call put_shared_variable('profy_ffree',profy_ffree)
        call put_shared_variable('profz_ffree',profz_ffree)
      endif
!
!  Communicate whether reduced sound speed is used to entropy.
!
      call put_shared_variable('lreduced_sound_speed',lreduced_sound_speed)
!
      if (lreduced_sound_speed) then
        call put_shared_variable('reduce_cs2',reduce_cs2)
        call put_shared_variable('lscale_to_cs2top',lscale_to_cs2top)
      endif

      call put_shared_variable('beta_glnrho_global',beta_glnrho_global)
      call put_shared_variable('beta_glnrho_scaled',beta_glnrho_scaled)

      if (lreference_state) then
        call put_shared_variable('reference_state',reference_state)
        call put_shared_variable('reference_state_mass',reference_state_mass)
      endif


    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  For compatibility with other applications, we keep the possibility
!  of giving diffrho units of dxmin*cs0, but cs0 is not well defined general
!
!  24-nov-02/tony: coded
!  31-aug-03/axel: normally, diffrho should be given in absolute units
!  14-apr-14/axel: lreinitialize_lnrho and lreinitialize_rho added
!  21-jan-15/MR: changes for use of reference state.
!  10-feb-15/MR: added getting reference state from initial condition.
!                upwind switch according to log/nolog (with warning).
!  15-nov-16/fred: option to apply z-profile to reinitialize_*
!  25-may-18/fred: definitive test of mass diffusion correction implemented
!
      use EquationOfState, only: select_eos_variable, get_gamma_etc
      use BorderProfiles, only: request_border_driving
      use Deriv, only: der,der2
      use FArrayManager
      use Gravity, only: lnumerical_equilibrium
      use Sub, only: stepdown,der_stepdown, erfunc,step
      use SharedVariables, only: put_shared_variable, get_shared_variable
      use InitialCondition, only: initial_condition_all
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      real, dimension (nzgrid) :: tmpz
!
      integer :: i,j,m,n, stat
      logical :: lnothing, exist
      real :: rho_bot,sref
      real, dimension(:), pointer :: gravx_xpencil
      real :: cp
!
!  Prevent this module when background stratification is on.
!
      if (lstratz) call fatal_error('initialize_density','lstratz=.true. -> use density_stratified instead')
!
!  Check the switches.
!
      if (ldensity_nolog.and.lupw_lnrho) then
        lupw_rho=.true.
        call warning('initialize_density','enabled upwinding for linear density')
      elseif (.not.ldensity_nolog.and.lupw_rho) then
        lupw_lnrho=.true.
        call warning('initialize_density','enabled upwinding for log density')
      endif
!
      if (.not.ldensity_nolog.and.lweno_transport) then
        lweno_transport=.false.
        !call warning('initialize_density','disabled WENO transport for logarithmic density')
        call fatal_error('initialize_density','cannot do WENO transport for logarithmic density!')
      endif

      lreference_state = ireference_state/='nothing'
      lfullvar_in_slices = lfullvar_in_slices.and.lreference_state
      lsubstract_reference_state = lsubstract_reference_state.and.lreference_state
!
!  If density variable is actually deviation from reference state, log(density) cannot be used.
!
      if (lreference_state) then
        if (.not.ldensity_nolog) &
          call fatal_error('initialize_density','use of reference state requires use of linear density')

        lcheck_negative_density=.false.

        if (.not.lentropy) then
!
!  Use of reference state at the moment implemented only for entropy.
!
          call fatal_error('initialize_density','use of reference state requires use of entropy')
          if (ltemperature) then
            if (.not.ltemperature_nolog) &
              call fatal_error('initialize_density','use of reference state requires use of linear temperature')
          endif
        endif
      endif
!
!  Initialize cs2cool to cs20
!  (currently disabled, because it causes problems with mdarf auto-test)
!     cs2cool=cs20
!
!  Made to work by adding diffrho + cdiffrho to the rprint reset list.
!
      if (diffrho==0.0) diffrho=cdiffrho*dxmin*cs0
!
!  Define stratification Gamma (upper case, as in the manual).
!  Set mpoly to default value, unless ggamma is set.
!
      if (ggamma==impossible) then
        if (mpoly==impossible) mpoly=1.5
        ggamma=1.+1./mpoly
      else
        if (ggamma/=1.) mpoly=1./(ggamma-1.)
      endif
!
!  Turn off continuity equation if we are not using hydrodynamics.
!
      if (.not.lhydro .and. .not.lhydro_kinematic) then
        lcontinuity_gas=.false.
        if (lroot) print*, 'initialize_density: no hydro, turned off continuity equation'
      endif
!
!  Turn off continuity equation term for 0-D runs.  (MR: What about Schur?)
!
      if (nwgrid==1) then
        lcontinuity_gas=.false.
        if (lroot) print*, 'initialize_density: 0-D run, turned off continuity equation'
      endif
!
!  Rescale density by a factor rescale_rho.
!
      if (lreinitialize_lnrho.or.lreinitialize_rho) then
        do j=1,ninit
          if (ldensity_nolog) then
            select case (initlnrho(j))
            case ('rescale'); f(:,:,:,irho)=f(:,:,:,irho)*rescale_rho
            case ('zprofile')
              inquire(file='zprof.txt',exist=exist)
              if (exist) then
                open(31,file='zprof.txt')
              else
                inquire(file=trim(directory)//'/zprof.ascii',exist=exist)
                if (exist) then
                  open(31,file=trim(directory)//'/zprof.ascii')
                else
                  call fatal_error('reinitialize_rho','error - no zprof.txt input file')
                endif
              endif
              do n=1,nzgrid
                read(31,*,iostat=stat) tmpz(n)
                if (stat<0) exit
              enddo
              do n=n1,n2
                f(:,:,n,irho)=f(:,:,n,irho)*rescale_rho*tmpz(n-nghost+nz*ipz)
              enddo
            endselect
          else
            select case (initlnrho(j))
            case ('rescale'); f(:,:,:,ilnrho)=f(:,:,:,ilnrho)+log(rescale_rho)
            case ('zprofile')
              inquire(file='zprof.txt',exist=exist)
              if (exist) then
                open(31,file='zprof.txt')
              else
                inquire(file=trim(directory)//'/zprof.ascii',exist=exist)
                if (exist) then
                  open(31,file=trim(directory)//'/zprof.ascii')
                else
                  call fatal_error('reinitialize_rho','error - no zprof.txt file')
                endif
              endif
              do n=1,nzgrid
                read(31,*,iostat=stat) tmpz(n)
                if (stat<0) exit
              enddo
              do n=n1,n2
                f(:,:,n,ilnrho)=f(:,:,n,ilnrho)+log(rescale_rho*tmpz(n-nghost+nz*ipz))
              enddo
            case ('particle_profile')
              if (iapn(1)/=0) f(:,:,:,ilnrho)=alog(f(:,:,:,iapn(1)))
            endselect
          endif
        enddo
      endif
!
!  Compute mask for x-averaging where x is in density_xaver_range.
!  Normalize such that the average over the full domain
!  gives still unity.
!
      if (l1 == l2) then
        xmask_den = 1.
      else
        where (x(l1:l2) >= density_xaver_range(1) .and. x(l1:l2) <= density_xaver_range(2))
          xmask_den = 1.
        elsewhere
          xmask_den = 0.
        endwhere
        density_xaver_range(1) = max(density_xaver_range(1), xyz0(1))
        density_xaver_range(2) = min(density_xaver_range(2), xyz1(1))
        if (lspherical_coords) then
          xmask_den = xmask_den * (xyz1(1)**3 - xyz0(1)**3) &
              / (density_xaver_range(2)**3 - density_xaver_range(1)**3)
        elseif (lcylindrical_coords) then
          xmask_den = xmask_den * (xyz1(1)**2 - xyz0(1)**2) &
              / (density_xaver_range(2)**2 - density_xaver_range(1)**2)
        else
          xmask_den = xmask_den * Lxyz(1) / (density_xaver_range(2) - density_xaver_range(1))
        endif
      endif
!
!  Compute mask for z-averaging where z is in density_zaver_range.
!  Normalize such that the average over the full domain
!  gives still unity.
!
      if (n1 == n2) then
        zmask_den = 1.
      else
        where (z(n1:n2) >= density_zaver_range(1) .and. z(n1:n2) <= density_zaver_range(2))
          zmask_den = 1.
        elsewhere
          zmask_den = 0.
        endwhere
        density_zaver_range(1) = max(density_zaver_range(1), xyz0(3))
        density_zaver_range(2) = min(density_zaver_range(2), xyz1(3))
        zmask_den = zmask_den * Lxyz(3) / (density_zaver_range(2) - density_zaver_range(1))
      endif
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'xmask_den=',xmask_den
        print*,'zmask_den=',zmask_den
      endif
!
!  Initialize mass diffusion.
!
      ldiff_normal=.false.
      ldiff_cspeed=.false.
      ldiff_shock=.false.
      ldiff_hyper3=.false.
      ldiff_hyper3lnrho=.false.
      ldiff_hyper3_aniso=.false.
      ldiff_hyper3_polar=.false.
      ldiff_hyper3_strict=.false.
      ldiff_hyper3lnrho_strict=.false.
      ldiff_hyper3_mesh=.false.
      ldensity_slope_limited=.false.
!
!  initialize lnothing. It is needed to prevent multiple output.
!
      lnothing=.false.
!
!  Different choices of mass diffusion (if any).
!
      do i=1,ndiff_max
        select case (idiff(i))
        case ('normal')
          if (lroot) print*,'diffusion: div(D*grad(rho))'
          ldiff_normal=.true.
        case ('cspeed')
          if (lroot) print*,'diffusion: div(D*grad(rho))'
          ldiff_cspeed=.true.
        case ('hyper3')
          if (lroot) print*,'diffusion: (d^6/dx^6+d^6/dy^6+d^6/dz^6)rho'
          ldiff_hyper3=.true.
        case ('hyper3lnrho','hyper3-lnrho')
          if (lroot) print*,'diffusion: (d^6/dx^6+d^6/dy^6+d^6/dz^6)lnrho'
          ldiff_hyper3lnrho=.true.
       case ('hyper3_aniso','hyper3-aniso')
          if (lroot.and.ldensity_nolog) &
               print*,'diffusion: (Dx*d^6/dx^6 + Dy*d^6/dy^6 + Dz*d^6/dz^6)rho'
          if (lroot.and..not.ldensity_nolog) &
               print*,'diffusion: (Dx*d^6/dx^6 + Dy*d^6/dy^6 + Dz*d^6/dz^6)lnrho'
          ldiff_hyper3_aniso=.true.
        case ('hyper3_cyl','hyper3-cyl','hyper3_sph','hyper3-sph')
          if (lroot) print*,'diffusion: Dhyper/pi^4 *(Delta(rho))^6/Deltaq^2'
          ldiff_hyper3_polar=.true.
        case ('hyper3-strict','hyper3_strict')
          if (lroot) print*,'diffusion: Dhyper*del2(del2(del2(rho)))'
          ldiff_hyper3_strict=.true.
        case ('hyper3-lnrho-strict','hyper3_lnrho_strict')
          if (lroot) print*,'diffusion: Dhyper*del2(del2(del2(lnrho)))'
          ldiff_hyper3lnrho_strict=.true.
        case ('hyper3_mesh','hyper3-mesh')
          if (lroot) print*,'diffusion: mesh hyperdiffusion'
          ldiff_hyper3_mesh=.true.
        case ('shock','diff-shock','diffrho-shock')
          if (lroot) print*,'diffusion: shock diffusion'
          ldiff_shock=.true.
        case ('density-slope-limited')
          if (lroot) print*,'mass diffusion: slope limited diffusion'
          if (lroot) print*,'mass diffusion: using ',trim(div_sld_dens),' order'
            ldensity_slope_limited=.true.
            lmassdiff_fix=.true.
        case ('','none')
          if (lroot .and. (.not. lnothing)) print*,'diffusion: nothing (i.e. no mass diffusion)'
        case default
          call fatal_error('initialize_density','No such idiff('//trim(itoa(i))//'): '//trim(idiff(i)))
        endselect
        lnothing=.true.
      enddo
!
!  If we're timestepping, die or warn if the the diffusion coefficient that
!  corresponds to the chosen diffusion type is not set.
!
      if (lrun) then
        if ((ldiff_normal.or.ldiff_cspeed).and.diffrho==0.0) &
            call warning('initialize_density','Diffusion coefficient diffrho is zero!', 0)
        if (ldiff_cspeed.and..not.(lentropy.or.ltemperature)) &
            call warning('initialize_density','Diffusion with cspeed can only be used with lenergy!', 0)
        if ( (ldiff_hyper3 .or. ldiff_hyper3lnrho .or. ldiff_hyper3_strict .or. &
              ldiff_hyper3lnrho_strict) .and. diffrho_hyper3==0.0) &
            call fatal_error('initialize_density','Diffusion coefficient diffrho_hyper3 is zero!')
        if ( (ldiff_hyper3_aniso) .and.  &
             ((diffrho_hyper3_aniso(1)==0. .and. nxgrid/=1 ).or. &
              (diffrho_hyper3_aniso(2)==0. .and. nygrid/=1 ).or. &
              (diffrho_hyper3_aniso(3)==0. .and. nzgrid/=1 )) ) &
            call fatal_error('initialize_density','A diffusion coefficient of diffrho_hyper3 is zero!')
        if (ldiff_shock .and. diffrho_shock==0.0) &
            call fatal_error('initialize_density','diffusion coefficient diffrho_shock is zero!')
        if (ldiff_normal.or.ldiff_cspeed.or.ldiff_shock) then
          if (.not.lmassdiff_fix) call warning('initialize_density', &
            'For diffusion energy/momentum correction should use lmassdiff_fix=T', 0)
        endif
        if (lmassdiff_fixkin.or.lmassdiff_fixmom) then
          lmassdiff_fix=.true.
          call information('initialize_density','lmassdiff_fix=T now the default', 0)
        endif
!
!  Dynamical hyper-diffusivity operates only for mesh formulation of hyper-diffusion
!
        if (ldynamical_diffusion.and..not.ldiff_hyper3_mesh) then
          call fatal_error("initialize_density", &
               "Dynamical diffusion requires mesh hyper-diffusion, switch idiff='hyper3-mesh'")
        endif
!
      endif
!
      if (lfreeze_lnrhoint) lfreeze_varint(ilnrho)    = .true.
      if (lfreeze_lnrhoext) lfreeze_varext(ilnrho)    = .true.
      if (lfreeze_lnrhosqu) lfreeze_varsquare(ilnrho) = .true.
!
!  Tell the equation of state that we're here and what f variable we use.
!
      if (ldensity_nolog) then
        call select_eos_variable('rho',irho)
      else
        call select_eos_variable('lnrho',ilnrho)
      endif
!
!  Do not allow inconsistency between rho0 (from eos) and rho_const
!  or lnrho0 and lnrho_const. But this can only be of concern if
!  const_lnrho or const_rho are set in initlnrho.
!
      if (rho0/=rho_const .and. rho0/=impossible) then
        if (any(initlnrho=='const_lnrho') .or. any(initlnrho=='const_rho')) then
          call warning("initialize_density","are rho_const or lnrho_const ok?")
          if (lroot) then
            print*,"inconsistency between the density constants from eos  "
            print*,"(rho0 or lnrho0) and the ones from the density module "
            print*,"(rho_const or lnrho_const). It may damage your        "
            print*,"simulation if you are using them in different places. "
          endif
        else
          if (lroot) print*,'Note: rho_const or lnrho_const are not used'
        endif
      endif
!
      if (lnumerical_equilibrium) then
         if (lroot) print*,'initializing global gravity in density'
         call farray_register_global('gg',iglobal_gg,vector=3)
      endif
!
!  To counteract the shock diffusion of the mean stratification, the initial stratification
!  needs to be read from file.
!
      if (lanti_shockdiffusion) then
        if ( .not. lwrite_stratification) &
          call fatal_error('initialize_density','must have lwrite_stratification for anti shock diffusion')
      else
        call information('initialize_density','no need to read initial stratification for '// &
                                              'lanti_shockdiffusion=F')
      endif
!
!  Possible to read initial stratification from file.
!
      if (lrun .and. lwrite_stratification) then
        if (ldensity_nolog) then
          call fatal_error('initialize_density', &
            'currently only possible to read *logarithmic* stratification from file')
        else
          if (lroot) print*, 'initialize_density: reading original stratification from stratification.dat'
          open(19,file=trim(directory_dist)//'/stratification.dat')
          read(19,*) lnrho_init_z_nz
          close(19)
        endif
!
!  Need to precalculate some terms for anti shock diffusion.
!
        if (lanti_shockdiffusion) then
          call der(3,lnrho_init_z,dlnrhodz_init_z)
          call der2(3,lnrho_init_z,del2lnrho_glnrho2_init_z(n1:n2))
          del2lnrho_glnrho2_init_z=del2lnrho_glnrho2_init_z+dlnrhodz_init_z**2
        endif
      endif
!
!  Possible to store non log rho as auxiliary variable.
!
      if (lrho_as_aux) then
        if (ldensity_nolog) then
          call warning('initialize_density', &
            'makes no sense to have lrho_as_aux=T if evolving linear rho - ignored')
          lrho_as_aux=.false.
        else
          call farray_register_auxiliary('rho',irho,communicated=.true.)
        endif
      endif
!
!  For diffusion term with non-logarithmic density we need to save rho
!  as an auxiliary variable.
!
      if (ldiffusion_nolog .and. .not. lrho_as_aux) then
        if (lroot) then
          print*, 'initialize_density: must have lrho_as_aux=T for non-logarithmic diffusion'
          print*, '  (consider setting lrho_as_aux=T and'
          print*, '   !  MAUX CONTRIBUTION 1'
          print*, '   !  COMMUNICATED AUXILIARIES 1'
          print*, '   in cparam.local)'
        endif
        call fatal_error('initialize_density','')
      endif
!
!  Tell the BorderProfiles module if we intend to use border driving, so
!  that the module can request the right pencils.
!
      select case (borderlnrho)
!
      case ('zero','0','constant','initial-condition')
        call request_border_driving(borderlnrho)
      case ('nothing')
        if (headtt.and.ip<=5) print*,"initialize_density: borderlnrho='nothing'"
      case default
        call fatal_error('set_border_density',"No such borderlnrho: '"//trim(borderlnrho)//"'")
      end select
!
!  Check if we are solving partially force-free equations.
!
      if (lffree) then
        select case(ffree_profile)
        case('radial_stepdown')
          profx_ffree=1.+stepdown(x(l1:l2),rzero_ffree,wffree)
          dprofx_ffree=der_stepdown(x(l1:l2),rzero_ffree,wffree)
        case('y_stepdown')
          profy_ffree=1.+stepdown(y,rzero_ffree,wffree)
          dprofy_ffree=der_stepdown(y,rzero_ffree,wffree)
        case('z_stepdown')
          profz_ffree=1.+stepdown(z,rzero_ffree,wffree)
          dprofz_ffree=der_stepdown(z,rzero_ffree,wffree)
        case('surface_x')
          profx_ffree=0.5*(1.0-erfunc((x(l1:l2)-rzero_ffree)/wffree))
          dprofx_ffree=-exp(-((x(l1:l2)-rzero_ffree)/wffree)**2)/(sqrtpi*wffree)
        case('surface_y')
          profy_ffree=0.5*(1.0-erfunc((y-rzero_ffree)/wffree))
          dprofy_ffree=-exp(-((y-rzero_ffree)/wffree)**2)/(sqrtpi*wffree)
        case('surface_z')
          profz_ffree=0.5*(1.0-erfunc((z-rzero_ffree)/wffree))
          dprofz_ffree=-exp(-((z-rzero_ffree)/wffree)**2)/(sqrtpi*wffree)
        case('none')
          profx_ffree=1.
          profy_ffree=1.
          profz_ffree=1.
          dprofx_ffree=0.
          dprofy_ffree=0.
          dprofz_ffree=0.
        case default
          call fatal_error('initialize_density', 'lffree=T but no profile selected')
        endselect
      endif
!
      if (lreduced_sound_speed) then
!
        if (lscale_to_cs2top) then
          if (.not.leos_idealgas.and.lroot) &
            call warning('initialize density','lscale_to_cs2top not implemented for this eos')
          if (lgravx) reduce_cs2_profx=1./(rss_coef1*((x0+Lxyz(1))/x(l1:l2)-rss_coef2))
          if (lgravz) reduce_cs2_profz(n1:n2)=cs2top/(rss_coef1-rss_coef2*(z(n1:n2)-z0))
        else
          reduce_cs2_profx=1.
          reduce_cs2_profz=1.
        endif
!
        reduce_cs2_profx=reduce_cs2*reduce_cs2_profx
!
      endif

      if (lmass_source) then
!
!  precalculate mass-source profiles.
!
        if (mass_source_profile(1:4)=='bump') fnorm=(2.*pi*mass_source_sigma**2)**1.5
        select case (mass_source_profile)
          case ('nothing')
            call not_implemented('initialize_density','mass source with no profile')
          case('bump2')
            fprofile_z=(mass_source_Mdot/fnorm)*exp(-.5*(z(n1:n2)/mass_source_sigma)**2)
          case('bumpx')
            fprofile_x=(mass_source_Mdot/fnorm)*exp(-.5*((x(l1:l2)-mass_source_offset)/mass_source_sigma)**2)
          case ('sph-step-down')
            if (.not.lspherical_coords) call fatal_error('initialize_density', &
                'you have chosen mass-source profiles "sph-step-down",'//achar(10)// &
                ' but coordinate system is not spherical!')
            fprofile_x=mass_source_Mdot*stepdown(r1_mn,rmax_mass_source,wdamp)
          case ('exponential','const','cylindric','bump')
!
! default to catch unknown values
!
          case default
            call fatal_error('initialize_density','no such mass_source_profile')
        endselect
      endif
!
!  Calculate profile functions (used as prefactors to turn off pressure
!  gradient term).
!
      if (ieos_profile=='nothing') then
        profz_eos=1.0
        dprofz_eos=0.0
      elseif (ieos_profile=='surface_z') then
        profz_eos=0.5*(1.0-erfunc(z/width_eos_prof))
        dprofz_eos=-exp(-(z/width_eos_prof)**2)/(sqrtpi*width_eos_prof)
      endif
!
!  For global density gradient beta=H/r*dlnrho/dlnr, calculate actual
!  gradient dlnrho/dr = beta/H.
!
      if (any(beta_glnrho_global/=0.0)) then
        beta_glnrho_scaled=beta_glnrho_global*Omega/cs0
        if (lroot) print*, 'initialize_density: Global density gradient with beta_glnrho_global=', &
                           beta_glnrho_global
      endif
!
      call get_gamma_etc(gamma,cp)
      gamma1=1./gamma; gamma_m1=gamma-1.; cp1=1./cp
!
      if (lreference_state) then
!
        select case(ireference_state)
        case ('adiabatic_simple_x')
!
!  Simple adiabatic reference state s=sref=const., p ~ rho^gamma for gravity ~ 1/x^2 --> rho/rho_bottom = (x/x_bottom)^(1/(1-gamma))
!
          rho_bot=1.
          sref=1.
          reference_state(:,iref_rho  ) =  rho_bot*(x(l1:l2)/xyz0(1))**(-1/gamma_m1)
          reference_state(:,iref_grho ) = -rho_bot/(gamma_m1*xyz0(1)) * (x(l1:l2)/xyz0(1))**(-gamma/gamma_m1)
          reference_state(:,iref_s    ) =  sref
          reference_state(:,iref_d2rho) =  rho_bot*gamma/(gamma_m1*xyz0(1))**2 * (x(l1:l2)/xyz0(1))**(-gamma/gamma_m1-1)
          reference_state(:,iref_d6rho) =  0.    ! yet missing
!
        case ('from_initcond')
!
!  Get from initial condition.
!
          call initial_condition_all(f,profiles=reference_state)
          reference_state(:,iref_d6rho) = 0.    ! yet missing
!
        case ('file')
!
!  Read from a file.
!
          call read_reference_state
!
        case default
          call not_implemented('initialize_density','type of reference state "'//trim(ireference_state)//'"')
        end select
!
!  Fetch gravity pencil. Requires that the gravity module is initialized before the density module
!  which is presently the case.
!
        call get_shared_variable('gravx_xpencil',gravx_xpencil)
        reference_state(:,iref_gp) = gravx_xpencil(l1:l2)*reference_state(:,iref_rho)
!
!  Compute the mass of the reference state for later use
!
        do n=n1,n2
          tmp=0.
          do m=m1,m2
            tmp=tmp+sum(reference_state(:,iref_rho)*dVol_x(l1:l2))*dVol_y(m)
          enddo
          reference_state_mass=reference_state_mass+tmp*dVol_z(n)
        enddo
!
        if (ncpus>1) then
          call mpiallreduce_sum(reference_state_mass,tmp)
          reference_state_mass=tmp
        endif
!
      endif
!
      call initialize_density_methods
!
!  Find the total mass for later use if lconserve_total_mass = .true.
!
      if (lrun .and. total_mass < 0.0) then
        total_mass = mean_density(f) * box_volume
        if (lreference_state) total_mass = total_mass + reference_state_mass
      endif
!
!  Check if we are solving for relativistic bulk motions, not just EoS.
!
      if (lhydro.and..not.lhydro_potential.and.iphiuu==0) then
        call get_shared_variable('lconservative', lconservative)
      else
        allocate(lconservative)
        lconservative=.false.
      endif
!
      if (lhydro.and..not.lhydro_potential) then
        call get_shared_variable('lhiggsless', lhiggsless)
      else
        allocate(lhiggsless)
        lhiggsless=.false.
      endif
!
       if (lhydro.and.lhiggsless) call get_shared_variable('eps_hless',eps_hless)

      if (lcontinuity_gas.and..not.lweno_transport.and.ldensity_nolog.and.lconservative.and..not.lhydro) &
        call fatal_error_local('initialize_density','divss not available without hydro')
!
      if ((idiag_rhof2mz/=0 .or. idiag_rhof2upmz/=0 .or. idiag_rhof2downmz/=0) &
          .and..not. lrho_flucz_as_aux) &
         call fatal_error('initialize_density','Need to set lrho_flucz_as_aux=T'// &
                          ' in density_run_pars for some density fluctuation diagnostics')
!
      if (density_floor>0.) density_floor_log=alog(density_floor)
      if (density_ceiling>0.) density_ceiling_log=alog(density_ceiling)
!
    endsubroutine initialize_density
!***********************************************************************
    subroutine init_lnrho(f)
!
!  Initialise logarithmic or non-logarithmic density.
!
!   7-nov-01/wolf: coded
!  28-jun-02/axel: added isothermal
!  15-oct-03/dave: added spherical shell (kws)
!
      use EquationOfState, only: eoscalc, ilnrho_TT
      use General, only: itoa,rtoa,complex_phase,notanumber
      use Gravity, only: zref,z1,z2,gravz,nu_epicycle,potential
      use Initcond
      use Mpicomm
      use Sub, only: blob
      use InitialCondition, only: initial_condition_lnrho
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: pot,prof
      real, dimension (ninit) :: lnrho_left,lnrho_right
      real :: lnrhoint,cs2int,pot0
      real :: pot_ext,lnrho_ext,cs2_ext,tmp1,k_j2
      real :: zbot,ztop,haut
      real, dimension (nx) :: r_mn,lnrho,TT,ss
      real, pointer :: gravx, rhs_poisson_const,fac_cs,cs2cool
      integer, pointer :: isothmid, isothtop
      complex :: omega_jeans
      integer :: j,ix,iy
      logical :: lnothing
!
      intent(inout) :: f
!
!  Sanity check.
!
      if (lread_oldsnap .and. ldensity_nolog .and. .not. all(initlnrho == 'nothing')) &
        call fatal_error('init_lnrho', 'cannot add initial conditions to the old snapshot')
        !MR: Why? should just be possible for linear density.
!
!  Define bottom and top height.
!
      zbot=xyz0(3)
      ztop=xyz0(3)+Lxyz(3)
!
!  Set default values for sound speed at top and bottom if nothing was
!  specified in density_init_pars.
!  These may be updated in one of the following initialization routines.
!
      if (cs2top==impossible) cs2top=cs20
      if (cs2bot==impossible) cs2bot=cs20
!
!  Different initializations of lnrho (called from start).
!
      lnrho0      = log(rho0)
      lnrho_left  = log(rho_left)
      lnrho_right = log(rho_right)

      if (lentropy) then
        call get_shared_variable('mpoly0', mpoly0, caller='init_lnrho')
        call get_shared_variable('mpoly1', mpoly1)
        call get_shared_variable('mpoly2', mpoly2)
      elseif ( any(initlnrho=='piecew-poly') .or. any(initlnrho=='piecew-disc') .or. any(initlnrho=='polytropic') .or. &
               any(initlnrho(:)(1:1)=='4') .or. any(initlnrho=='5') ) then
        call warning('init_lnrho','mpoly[0-2] not provided by entropy, take default 1.5')
        allocate(mpoly0,mpoly1,mpoly2)
        mpoly0=1.5; mpoly1=1.5; mpoly2=1.5
      endif
!
!  Initialize lnothing and cycle ninit (=4) times through the list of
!  initial conditions with the various options.
!
      lnothing=.true.
      do j=1,ninit
!
        if (initlnrho(j)=='nothing') cycle
        lnothing=.false.
!
        iinit_str=itoa(j)
!
        select case (initlnrho(j))
!
        case ('zero', '0'); f(:,:,:,ilnrho)=0.
        case ('const_lnrho'); f(:,:,:,ilnrho)=lnrho_const
        case ('const_rho'); f(:,:,:,ilnrho)=log(rho_const)
        case ('constant'); f(:,:,:,ilnrho)=log(rho_left(j))
        case ('linear_lnrho'); f(:,:,:,ilnrho)=lnrho_const-spread(spread(z,1,mx),2,my)/Hrho
        case ('exp_zbot'); f(:,:,:,ilnrho)=alog(rho_left(j))-spread(spread(z-zbot,1,mx),2,my)/Hrho
        case ('exp_rbot'); f(:,:,:,ilnrho)=lnrho_const-spread(spread(x-xyz0(1),2,my),3,mz)/Hrho
        case ('invsqr')
          do ix=1,mx
            if (x(ix)<=r0_rho) then
              f(ix,:,:,ilnrho)=0.0
            else
              f(ix,:,:,ilnrho)=2*log(r0_rho)-2*log(x(ix))
            endif
          enddo
        case ('invgravx')
          do ix=1,mx
            if (x(ix)<=r0_rho) then
              f(ix,:,:,ilnrho)=0.0
            else
              f(ix,:,:,ilnrho)=lnrho_const+invgrav_ampl/(x(ix))-invgrav_ampl/(r0_rho)
            endif
          enddo
        case ('x-point_xy')
          f(:,:,:,ilnrho)=f(:,:,:,ilnrho)-.5*ampllnrho(j)/cs20*( spread(spread(x**2,2,my),3,mz) &
                                                                +spread(spread(y**2,1,mx),3,mz) )
        case ('mode')
          call modes(ampllnrho(j),coeflnrho,f,ilnrho,kx_lnrho(j), ky_lnrho(j),kz_lnrho(j))
        case ('blob')
          call blob(ampllnrho(j),f,ilnrho,radius_lnrho(j),xblob(j),yblob(j),zblob(j))
        case ('blob_normalized')
          call blob(ampllnrho(j)/(sqrt2pi*radius_lnrho(j))**3,f,ilnrho, &
              radius_lnrho(j)*sqrt2,xblob(j),yblob(j),zblob(j))
        case ('blob_hs')
          if (lroot) print*, 'init_lnrho: put a blob in hydrostatic equilibrium: '// &
                             'radius_lnrho, ampllnrho, position=',radius_lnrho(j), &
                             ampllnrho(j), xblob(j), yblob(j), zblob(j)
          call blob(ampllnrho(j),f,ilnrho,radius_lnrho(j),xblob(j),yblob(j),zblob(j))
          call blob(-ampllnrho(j),f,iss,radius_lnrho(j),xblob(j),yblob(j),zblob(j))
        case ('pre-stellar'); call pre_stellar_cloud(f, datafile, mass_cloud, &
                     cloud_mode, T_cloud_out_rel, dens_coeff, temp_coeff, temp_trans, temp_coeff_out)
        case ('read_arr_file'); call read_outside_scal_array(f, "lnrho.arr", ilnrho)
        case ('isothermal'); call isothermal_density(f)
        case ('stratification'); call stratification(f,strati_type)
        case ('stratification-x'); call stratification_x(f,strati_type)
        case ('stratification-xz'); call stratification_xz(f,strati_type)
        case ('polytropic_simple'); call polytropic_simple(f)
        case ('stratification_tsallis'); call stratification_tsallis(f)
        case ('hydrostatic_TT'); call temp_hydrostatic(f,rho_const,gamma)
        case ('hydrostatic-z', '1')
          if (lroot) print*, 'init_lnrho: use polytropic_simple instead!'
        case ('xjump')
          call jump(f,ilnrho,lnrho_left(j),lnrho_right(j),widthlnrho(j),xjump_mid,yjump_mid,zjump_mid,'x')
        case ('yjump')
          call jump(f,ilnrho,lnrho_left(j),lnrho_right(j),widthlnrho(j),xjump_mid,yjump_mid,zjump_mid,'y')
        case ('zjump')
          call jump(f,ilnrho,lnrho_left(j),lnrho_right(j),widthlnrho(j),xjump_mid,yjump_mid,zjump_mid,'z')
        case ('xyjump')
          call jump(f,ilnrho,lnrho_left(j),lnrho_right(j),widthlnrho(j),xjump_mid,yjump_mid,zjump_mid,'xy')
        case ('x-y-jump')
           call jump(f,ilnrho,lnrho_left(j),lnrho_right(j),widthlnrho(j),xjump_mid,yjump_mid,zjump_mid,'x-y')
        case ('soundwave-x')
          call soundwave(ampllnrho(j),f,ilnrho,kx=kx_lnrho(j),width=widthlnrho(j))
        case ('soundwave-y')
          call soundwave(ampllnrho(j),f,ilnrho,ky=ky_lnrho(j))
        case ('soundwave-z')
          call soundwave(ampllnrho(j),f,ilnrho,kz=kz_lnrho(j))
        case ('sinwave-phase')
          call sinwave_phase(f,ilnrho,ampllnrho(j),kx_lnrho(j), &
                             ky_lnrho(j),kz_lnrho(j),phase_lnrho(j))
        case ('sinwave-phase-nolog')
          do m=m1,m2; do n=n1,n2
            f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + alog(1+amplrho(j)*sin(kx_lnrho(j)*x(l1:l2)+ &
                                  ky_lnrho(j)*y(m)+kz_lnrho(j)*z(n)+phase_lnrho(j)))
          enddo; enddo
        case ('coswave-phase')
          call coswave_phase(f,ilnrho,ampllnrho(j),kx_lnrho(j), &
                             ky_lnrho(j),kz_lnrho(j),phase_lnrho(j))
        case ('sinwave-x')
          call sinwave(ampllnrho(j),f,ilnrho,kx=kx_lnrho(j))
        case ('sinwave-y')
          call sinwave(ampllnrho(j),f,ilnrho,ky=ky_lnrho(j))
        case ('sinwave-z')
          call sinwave(ampllnrho(j),f,ilnrho,kz=kz_lnrho(j))
        case ('coswave-x')
          call coswave(ampllnrho(j),f,ilnrho,kx=kx_lnrho(j))
        case ('coswave-y')
          call coswave(ampllnrho(j),f,ilnrho,ky=ky_lnrho(j))
        case ('coswave-z')
          call coswave(ampllnrho(j),f,ilnrho,kz=kz_lnrho(j))
        case ('triquad')
          call triquad(ampllnrho(j),f,ilnrho,kx_lnrho(j), &
              ky_lnrho(j),kz_lnrho(j), kxx_lnrho(j), kyy_lnrho(j),kzz_lnrho(j))
        case ('isotdisk')
          call isotdisk(powerlr,f,ilnrho,zoverh, hoverr)
          f(1:mx,1:my,1:mz,iss)=-(gamma-1)/gamma*f(1:mx,1:my,1:mz,ilnrho)
!          call isotdisk(powerlr,f,iss,zoverh,hoverr, -(gamma-1)/gamma)
        case ('sinx_siny_sinz')
          call sinx_siny_sinz(ampllnrho(j),f,ilnrho,kx_lnrho(j),ky_lnrho(j),kz_lnrho(j))
        case ('corona'); call corona_init(f,gamma)
        case ('gaussian3d')
          call gaussian3d(ampllnrho(j),f,ilnrho,radius_lnrho(j))
        case ('gaussian-z')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) - z(n)**2/(2*radius_lnrho(j)**2)
          enddo; enddo
        case ('gauss-z-offset')
          do n=n1,n2
             f(:,:,n,ilnrho) = f(:,:,n,ilnrho) + alog(exp(f(:,:,n,ilnrho))+ &
                ampllnrho(j)*(exp(-(z(n)+lnrho_z_shift)**2/(2*radius_lnrho(j)**2))))
          enddo
        case ('gaussian-noise')
          if (lnrho_left(j) /= 0.) f(:,:,:,ilnrho)=lnrho_left(j)
          call gaunoise(ampllnrho(j),f,ilnrho,ilnrho)
        case ('gaussian-noise-rprof')
          call gaunoise_rprof(ampllnrho(j),f,ilnrho,rnoise_int,rnoise_ext)
!
!  1/cosh^2 profile
!
        case ('cosh21-z')
          f(:,:,:,ilnrho)=spread(spread(max(lnrho_const,alog(1./cosh(kz_lnrho(j)*z/sqrt2)**2)),1,mx),2,my)
!
!  use code to plot EoS
!
        case ('lnrho_vs_lnT')
          if (ilnTT==0) call fatal_error("init_lnrho","ilnTT==0")
          f(:,:,:,ilnrho)=spread(spread(y,1,mx),3,mz)
          f(:,:,:,ilnTT)=spread(spread(x,2,my),3,mz)
!
!  use code to plot EoS
!
        case ('lnrho_vs_ss')
          if (iss==0) call fatal_error("init_lnrho","iss==0")
          f(:,:,:,ilnrho)=spread(spread(y,1,mx),3,mz)
          f(:,:,:,iss)=spread(spread(x,2,my),3,mz)
!
!  Noise, but just x-dependent.
!
        case ('gaussian-noise-x')
          call gaunoise(ampllnrho(j),f,ilnrho,ilnrho)
          f(:,:,:,ilnrho)=spread(spread(f(:,4,4,ilnrho),2,my),3,mz) !(watch 1-D)
!
!  Density jump (for shocks).
!
        case ('rho-jump-z', '2')
          if (lroot) print*, 'init_lnrho: density jump; rho_left,right=', rho_left(j), rho_right(j)
          if (lroot) print*, 'init_lnrho: density jump; widthlnrho=', widthlnrho(j)
          do n=n1,n2; do m=m1,m2
            prof=0.5*(1.0+tanh(z(n)/widthlnrho(j)))
            f(l1:l2,m,n,ilnrho)=log(rho_left(j))+log(rho_left(j)/rho_right(j))*prof
          enddo; enddo
!
!  A*tanh(y/d) profile
!
        case ('tanhy')
          if (lroot) print*,'init_lnrho: tangential discontinuity'
          do m=m1,m2
            prof=ampllnrho(j)*tanh(y(m)/widthlnrho(j))
            do n=n1,n2
              f(l1:l2,m,n,ilnrho)=prof
            enddo
          enddo
!
!  Hydrostatic density stratification for isentropic atmosphere.
!
        case ('hydrostatic-z-2', '3')
          if (lgravz) then
            if (lroot) print*,'init_lnrho: vertical density stratification'
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,ilnrho) = -grads0*z(n) + 1./gamma_m1*log( 1 + gamma_m1*gravz/grads0/cs20 &
                                    *(1-exp(-grads0*z(n))) )
            enddo; enddo
          endif
        case ('hydrostatic-r');   call init_hydrostatic_r (f)
!
!  Hydrostatic radial density stratification for isentropic (or
!  isothermal) sphere.
!
        case ('sph_isoth'); call init_sph_isoth (f)
!
        case ('cylind_isoth')
          call get_shared_variable('gravx', gravx, caller='init_lnrho')
          if (lroot) print*, 'init_lnrho: isothermal cylindrical ring with gravx=', gravx
          haut=-cs20/gamma/gravx
          TT=spread(cs20/gamma_m1,1,nx)
          do n=n1,n2
          do m=m1,m2
            lnrho=lnrho0-(x(l1:l2)-r_ext)/haut
            f(l1:l2,m,n,ilnrho)=lnrho
            call eoscalc(ilnrho_TT,lnrho,TT,ss=ss)
            f(l1:l2,m,n,iss)=ss
          enddo
          enddo
        case ('isentropic-star')
!
!  Isentropic/isothermal hydrostatic sphere"
!    ss  = 0       for r<R,
!    cs2 = const   for r>R
!
!  Only makes sense if both initlnrho=initss='isentropic-star'
!
          if (lgravr) then
            if (lentropy) then
              call get_shared_variable('cs2cool', cs2cool, caller='init_lnrho')
            else
              call warning('init_lnrho','cs2cool not provided by entropy, set zero')
              allocate(cs2cool); cs2cool=0.
            endif
            if (lroot) print*, 'init_lnrho: isentropic star with isothermal atmosphere'
            do n=n1,n2; do m=m1,m2
              r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
              call potential(POT=pot,POT0=pot0,RMN=r_mn) ! gravity potential
!
!  rho0, cs0, pot0 are the values in the centre
!
              if (gamma /= 1) then
!  Note:
!  (a) `where' is expensive, but this is only done at
!      initialization.
!  (b) Comparing pot with pot_ext instead of r with r_ext will
!      only work if grav_r<=0 everywhere -- but that seems
!      reasonable.
                call potential(R=r_ext,POT=pot_ext) ! get pot_ext=pot(r_ext)
!  Do consistency check before taking the log() of a potentially
!  negative number
                tmp1 = 1 - gamma_m1*(pot_ext-pot0)/cs20
                if (tmp1 <= 0.) &
                  call fatal_error('init_lnrho', 'Imaginary density values: '//achar(10)// &
                                   'Trying to calculate log('//rtoa(tmp1)// &
                                   ') for r_ext -- need to increase cs20?')
                lnrho_ext = lnrho0 + log(tmp1) / gamma_m1
                cs2_ext   = cs20*tmp1
!  Adjust for given cs2cool (if given) or set cs2cool (otherwise)
                if (cs2cool/=0) then
                  lnrho_ext = lnrho_ext - log(cs2cool/cs2_ext)
                else
                  cs2cool   = cs2_ext
                endif
!
!  Add temperature and entropy jump (such that pressure
!  remains continuous) if cs2cool was specified in start.in:
!
                where (pot <= pot_ext) ! isentropic for r<r_ext
                  f(l1:l2,m,n,ilnrho) = lnrho0 + log(1 - gamma_m1*(pot-pot0)/cs20) / gamma_m1
                elsewhere           ! isothermal for r>r_ext
                  f(l1:l2,m,n,ilnrho) = lnrho_ext - gamma*(pot-pot_ext)/cs2cool
                endwhere
              else                  ! gamma=1 --> simply isothermal (I guess [wd])
                f(l1:l2,m,n,ilnrho) = lnrho0 - (pot-pot0)/cs20
              endif
            enddo; enddo
          endif
!
        case ('piecew-poly', '4')
!
!  Piecewise polytropic for stellar convection models.
!
          if (lroot) print*, 'init_lnrho: piecewise polytropic vertical stratification (lnrho)'
!  Top region.
          cs2int = cs0**2
          lnrhoint = lnrho0
!
          call get_shared_variable('isothmid', isothmid, caller='init_lnrho')
          call get_shared_variable('fac_cs', fac_cs)
          if (lentropy) then
            call get_shared_variable('isothtop', isothtop)
          else
            call warning('init_lnrho','isothtop not provided by entropy, set to zero')
            allocate(isothtop); isothtop=0
          endif
!
          f(:,:,:,ilnrho) = lnrho0 ! just in case
          call polytropic_lnrho_z(f,mpoly2,zref,z2,ztop+Lz,isothtop,cs2int,lnrhoint,fac_cs)
!  Unstable layer.
          call polytropic_lnrho_z(f,mpoly0,z2,z1,z2,isothmid,cs2int,lnrhoint)
!  Stable layer.
          call polytropic_lnrho_z(f,mpoly1,z1,z0,z1,0,cs2int,lnrhoint)
!
!  Calculate cs2bot and cs2top for run.x (boundary conditions).
!
          cs2bot = cs2int + gamma/gamma_m1*gravz/(mpoly2+1)*(zbot-z0  )
          if (isothtop /= 0) then
            cs2top = cs20
          else
            cs2top = cs20 + gamma/gamma_m1*gravz/(mpoly0+1)*(ztop-zref)
          endif
        case ('piecew-disc', '41')
!
!  Piecewise polytropic for accretion discs.
!
          if (lroot) print*, 'init_lnrho: piecewise polytropic disc stratification (lnrho)'
!  Bottom region.
          cs2int = cs0**2
          lnrhoint = lnrho0
          f(:,:,:,ilnrho) = lnrho0 ! just in case
          call polytropic_lnrho_disc(f,mpoly1,zref,z1,z1,0,cs2int,lnrhoint)
!  Unstable layer.
          call polytropic_lnrho_disc(f,mpoly0,z1,z2,z2,0,cs2int,lnrhoint)
!  Stable layer (top).
          call polytropic_lnrho_disc(f,mpoly2,z2,ztop,ztop,isothtop,cs2int,lnrhoint)
!
!  Calculate cs2bot and cs2top for run.x (boundary conditions).
!
!  cs2bot = cs2int + gamma/gamma_m1*gravz*nu_epicycle**2/(mpoly2+1)* &
!         (zbot**2-z0**2)/2.
          cs2bot = cs20
          if (isothtop /= 0) then
            cs2top = cs20
          else
            cs2top = cs20 + gamma/gamma_m1*gravz*nu_epicycle**2/(mpoly0+1)*(ztop**2-zref**2)/2.
          endif
        case ('polytropic', '5')
!
!  Polytropic stratification.
!  cs0, rho0 and ss0=0 refer to height z=zref
!
          if (lroot) print*, 'init_lnrho: polytropic vertical stratification (lnrho)'
!
          cs2int = cs20
          lnrhoint = lnrho0
          f(:,:,:,ilnrho) = lnrho0 ! just in case
!  Only one layer.
          call polytropic_lnrho_z(f,mpoly0,zref,z0,z0+2*Lz,0,cs2int,lnrhoint)
!
!  Calculate cs2bot and cs2top for run.x (boundary conditions).
!
          cs2bot = cs20 + gamma*gravz/(mpoly0+1)*(zbot-zref)
          cs2top = cs20 + gamma*gravz/(mpoly0+1)*(ztop-zref)
        case ('sound-wave', '11')
!
!  Sound wave (should be consistent with hydro module).
!
          if (lroot) print*,'init_lnrho: x-wave in lnrho; ampllnrho=', ampllnrho(j)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=lnrho_const+ampllnrho(j)*sin(kx_lnrho(j)*x(l1:l2))
          enddo; enddo
        case ('sound-wave-exp')
!
!  Sound wave (should be consistent with hydro module).
!
          if (lroot) print*,'init_lnrho: x-wave in rho; ampllnrho=', ampllnrho(j)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=log(rho_const+amplrho(j)*sin(kx_lnrho(j)*x(l1:l2)))
          enddo; enddo
        case ('sound-wave2')
!
!  Sound wave (should be consistent with hydro module).
!
          if (lroot) print*,'init_lnrho: x-wave in lnrho; ampllnrho=', ampllnrho(j)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=lnrho_const+ampllnrho(j)*cos(kx_lnrho(j)*x(l1:l2))
          enddo; enddo
        case ('shock-tube', '13')
!
!  Shock tube test (should be consistent with hydro module).
!  Remember that for ldensity_nolog=T, we say exp(...) at the end of this routine!
!
          call information('init_lnrho','polytropic standing shock')
          do n=n1,n2; do m=m1,m2
            prof=0.5*(1.+tanh(x(l1:l2)/widthlnrho(j)))
            f(l1:l2,m,n,ilnrho)=log(rho_left(j)) + (log(rho_right(j))-log(rho_left(j)))*prof
          enddo; enddo
        case ('nolog-shock-tube')
!
!  Shock tube test (should be consistent with hydro module).
!  Remember that for ldensity_nolog=T, we say exp(...) at the end of this routine!
!
          call information('init_lnrho','standing shock')
          do n=n1,n2; do m=m1,m2
            prof=0.5*(1.+tanh(x(l1:l2)/widthlnrho(j)))
            f(l1:l2,m,n,ilnrho)=log(rho_left(j)+(rho_right(j)-rho_left(j))*prof)
          enddo; enddo
        case ('sin-xy')
!
!  sin profile in x and y.
!
          call information('init_lnrho','lnrho=sin(x)*sin(y)')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=log(rho0) + ampllnrho(j)*sin(kx_lnrho(j)*x(l1:l2))*sin(ky_lnrho(j)*y(m))
          enddo; enddo
        case ('sin-xy-rho')
!
!  sin profile in x and y, but in rho, not ln(rho).
!
          call information('init_lnrho','rho=sin(x)*sin(y)')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=log(rho0*(1+ampllnrho(j)*sin(kx_lnrho(j)*x(l1:l2))*sin(ky_lnrho(j)*y(m))))
          enddo; enddo
        case ('Schur_rho')
!
!  3d3d1d Schur flow init density, but in rho, not ln(rho).
!
          call information('init_lnrho','rho=(1+0.5*Sin(x+y))*(1+0.5*Sin(z))')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho)=log(rho0*(1+ampllnrho(j)*sin((kx_lnrho(j)*x(l1:l2))+ky_lnrho(j)*y(m)) &
                                *(1+ampllnrho(j)*sin(kz_lnrho(j)*z(n)))))
          enddo; enddo
        case ('linear')
!
!  Linear profile in kk.xxx.
!
          call information('init_lnrho','linear profile')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho) = log(rho0) + ampllnrho(j)*(kx_lnrho(j)*x(l1:l2)+ &
                ky_lnrho(j)*y(m)+kz_lnrho(j)*z(n))/sqrt(kx_lnrho(j)**2+ky_lnrho(j)**2+kz_lnrho(j)**2)
          enddo; enddo
        case ('planet')
!
!  Planet solution of Goodman, Narayan & Goldreich (1987).
!  (Simple 3-D)
!
          call planet(rbound,f,eps_planet,radius_lnrho(j),gamma,cs20,rho0,widthlnrho(j),hh0)
        case ('planet_hc')
!
!  Planet solution of Goodman, Narayan & Goldreich (1987).
!  (3-D with hot corona)
!
          call planet_hc(amplrho(j),f,eps_planet,radius_lnrho(j), gamma,cs20,rho0,widthlnrho(j))
        case ('Ferriere')
          call information('init_lnrho','Ferriere set in entropy')
        case ('Galactic-hs')
          call information('init_lnrho','Galactic hydrostatic equilibrium setup done in entropy')
        case ('thermal-hs')
          call information('init_lnrho','thermal hydrostatic equilibrium setup done in interstellar')
        case ('geo-kws')
!
!  Radial hydrostatic profile in shell region only.
!
          call information('init_lnrho','kws hydrostatic in spherical shell region')
          call shell_lnrho(f)
        case ('geo-kws-constant-T','geo-benchmark')
!
!  Radial hydrostatic profile throughout box, which is consistent
!  with constant temperature in exterior regions, and gives continuous
!  density at shell boundaries.
!
          call information('init_lnrho','kws hydrostatic in spherical shell and exterior')
          call shell_lnrho(f)
        case ('step_xz')
          call fatal_error('init_lnrho','neutron_star initial condition '// &
                                        'is now in special/neutron_star.f90')
        case ('jeans-wave-x')
!
!  Soundwave + self gravity.
!
          call get_shared_variable('rhs_poisson_const', rhs_poisson_const,caller='init_lnrho')
          omega_jeans = sqrt(cmplx(cs20*kx_lnrho(j)**2 - rhs_poisson_const*rho0,0.))/(rho0*kx_lnrho(j))
          if (lroot) print*,'Re(omega_jeans), Im(omega_jeans), Abs(omega_jeans)',&
                            real(omega_jeans),aimag(omega_jeans),abs(omega_jeans)
!
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho) = lnrho_const + ampllnrho(j)*sin(kx_lnrho(j)*x(l1:l2)+phase_lnrho(j))
            if (abs(omega_jeans)/=0.0) &
              f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + abs(omega_jeans*ampllnrho(j)) * &
                 sin(kx_lnrho(j)*x(l1:l2)+phase_lnrho(j) + complex_phase(omega_jeans*ampllnrho(j)))
          enddo; enddo
        case ('jeans-wave-oblique')
!
!  Soundwave + self gravity.
!
          call get_shared_variable('rhs_poisson_const', rhs_poisson_const,caller='init_lnrho')
          k_j2 = kx_lnrho(j)**2 + ky_lnrho(j)**2 + kz_lnrho(j)**2
          omega_jeans = sqrt(cmplx(cs20*k_j2 - rhs_poisson_const*rho0,0.))/(rho0*sqrt(k_j2))
          if (lroot) print*,'Re(omega_jeans), Im(omega_jeans), Abs(omega_jeans)', &
              real(omega_jeans),aimag(omega_jeans),abs(omega_jeans)
!
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho) = lnrho_const + ampllnrho(j)*sin(kx_lnrho(j)*x(l1:l2) + &
                                  ky_lnrho(j)*y(m) + kz_lnrho(j)*z(n))
            if (kx_lnrho(j)/=0) &
                f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + abs(omega_jeans*ampllnrho(j)) * &
                                   sin(kx_lnrho(j)*x(l1:l2)+complex_phase(omega_jeans*ampllnrho(j)))
            if (ky_lnrho(j)/=0) &
                f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + abs(omega_jeans*ampllnrho(j)) * &
                                   sin(ky_lnrho(j)*y(m)+complex_phase(omega_jeans*ampllnrho(j)))
            if (kz_lnrho(j)/=0) &
                f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + abs(omega_jeans*ampllnrho(j)) * &
                                   sin(kz_lnrho(j)*z(n)+complex_phase(omega_jeans*ampllnrho(j)))
          enddo; enddo
!
        case ('toomre-wave-x')
!
!  Soundwave + self gravity + (differential) rotation.
!
          call get_shared_variable('rhs_poisson_const', rhs_poisson_const,caller='init_lnrho')
          omega_jeans = sqrt(cmplx(cs20*kx_lnrho(j)**2 + &
              Omega**2 - rhs_poisson_const*rho0,0.))/(rho0*kx_lnrho(j))
!
          if (lroot) print*,'Re(omega_jeans), Im(omega_jeans), Abs(omega_jeans)', &
              real(omega_jeans),aimag(omega_jeans),abs(omega_jeans)
!
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnrho) = lnrho_const + &
              ampllnrho(j)*sin(kx_lnrho(j)*x(l1:l2))
            f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + &
                abs(omega_jeans*ampllnrho(j)) * &
                sin(kx_lnrho(j)*x(l1:l2)+complex_phase(omega_jeans*ampllnrho(j)))
            f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + &
                 abs(ampllnrho(j)* &
                 cmplx(0,-0.5*Omega/(kx_lnrho(j)*rho0))) * &
                 sin(kx_lnrho(j)*x(l1:l2)+complex_phase(ampllnrho(j)* &
                 cmplx(0,-0.5*Omega/(kx_lnrho(j)*rho0))))
          enddo; enddo
!
        case ('compressive-shwave')
!  Should be consistent with density
          f(:,:,:,ilnrho) = log(rho_const + f(:,:,:,ilnrho))
!
        case ('linz')
! linearly increasing with z
          do ix=l1,l2;do iy=m1,m2
            f(ix,iy,1:mz,ilnrho) = log( rho_bottom+((rho_top-rho_bottom)/(Lxyz(3)))*(z(1:mz)-xyz0(3)) )
          enddo;enddo
!
        case default
!
!  Catch unknown values
!
          call fatal_error('init_lnrho','No such initlnrho('//trim(iinit_str)//'): '//trim(initlnrho(j)))
!
        endselect
!
        if (lroot) print*,'init_lnrho: initlnrho('//trim(iinit_str)//') = '//trim(initlnrho(j))
!
      enddo  ! End loop over initial conditions.
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_lnrho(f)
!
      if (lnothing.and.lroot) print*,'init_lnrho: nothing'
!
!  check that cs2bot,cs2top are ok
!  for runs with ionization or fixed ionization, don't print them
!
      if (leos_ionization .or. leos_fixed_ionization) then
        cs2top=impossible
        cs2bot=impossible
      else
        if (lroot) print*,'init_lnrho: cs2bot,cs2top=',cs2bot,cs2top
      endif
!
!  If unlogarithmic density considered, take exp of lnrho resulting from
!  initlnrho.
!
      if (ldensity_nolog .and. (.not.lread_oldsnap) .and. (.not.ldensity_linearstart)) &
           f(:,:,:,irho)=exp(f(:,:,:,ilnrho))   !???
!AB: Wlad, why the question marks on Dec 7, 2020?
!
!  Impose density floor or ceiling if requested.
!
      call impose_density_floor(f)
      call impose_density_ceiling(f)
!
!  sanity check
!
      if (notanumber(f(l1:l2,m1:m2,n1:n2,ilnrho))) call error('init_lnrho','imaginary density values')
!
    endsubroutine init_lnrho
!***********************************************************************
    subroutine density_before_boundary(f)
!
!  Actions to take before boundary conditions are applied.
!
!   31-aug-09/MR: adapted from hydro_after_boundary
!    3-oct-12/MR: global averaging corrected
!    3-mar-14/MR: correction of total mass added
!   17-aug-14/MR: finalize_aver used
!   10-feb-15/MR: extended calc of <grad(log(rho))> to linear density;
!                 mass conservation slightly simplified
!   21-oct-15/MR: changes for slope-limited diffusion
!   03-apr-20/joern: restructured and fixed slope-limited diffusion:
!                    removed from here.
!
      use Sub, only: finalize_aver
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      real :: fact,cur_mass
      real, dimension (nx) :: tmp
!
      if ( (.not.ldensity_nolog) .and. (irho/=0) ) &
          f(l1:l2,m1:m2,n1:n2,irho)=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
!
!  Calculate mean (= xy-average) of lnrho.
!
      if (lcalc_lnrhomean) then
!
        fact=1./nxygrid
        lnrhomz=0.
        do n=n1,n2
          if (ldensity_nolog) then
            lnrhomz(n,1)=lnrhomz(n,1)+sum(alog(f(l1:l2,m1:m2,n,irho)))
          else
            lnrhomz(n,1)=lnrhomz(n,1)+sum(f(l1:l2,m1:m2,n,ilnrho))
          endif
        enddo
        lnrhomz=fact*lnrhomz
        call finalize_aver(nprocxy,12,lnrhomz)

        if (lrho_flucz_as_aux) then
          do n=n1,n2
            f(l1:l2,m1:m2,n,irho_flucz)=exp(f(l1:l2,m1:m2,n,ilnrho))-exp(lnrhomz(n,1))
          enddo
        endif
      endif
!
      if (lrmv) then
!
!  Force mass conservation if requested
!
        if (lconserve_total_mass .and. total_mass > 0.0) then
!
          cur_mass=box_volume*mean_density(f)
!
          if (lreference_state) then
            fact=total_mass/(cur_mass+reference_state_mass)
            tmp=(fact - 1.)*reference_state(:,iref_rho)
          else
            fact=total_mass/cur_mass
          endif
!
          if (ldensity_nolog) then
            if (lreference_state) then
              do n=n1,n2
                do m=m1,m2
                  f(l1:l2,m,n,irho) = f(l1:l2,m,n,irho)*fact + tmp
                enddo
              enddo
            else
              f(l1:l2,m1:m2,n1:n2,irho) = f(l1:l2,m1:m2,n1:n2,irho)*fact
            endif
          else
            f(l1:l2,m1:m2,n1:n2,ilnrho) = f(l1:l2,m1:m2,n1:n2,ilnrho)+alog(fact)
          endif
!
!  Conserve the momentum.
!
          if (lhydro) f(l1:l2,m1:m2,n1:n2,iux:iuz) = f(l1:l2,m1:m2,n1:n2,iux:iuz)/fact
!
        endif
      endif

      lupdate_mass_source = lmass_source .and. t>=tstart_mass_source .and. &
                            (tstop_mass_source==-1.0 .or. t<=tstop_mass_source)
!
   endsubroutine density_before_boundary
!***********************************************************************
    subroutine update_char_vel_density(f)
!
!  Updates characteristic velocity for slope-limited diffusion.
!
!  25-sep-15/MR+joern: coded
!   9-oct-15/MR: added updating of characteristic velocity by density
!                not yet activated (weight=0), tb reconsidered
!
      use General, only: staggered_mean_scal
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      real, parameter :: weight=.0
!
      if (lslope_limit_diff) call staggered_mean_scal(f,ilnrho,iFF_char_c,weight)
!
    endsubroutine update_char_vel_density
!***********************************************************************
    subroutine polytropic_lnrho_z(f,mpoly,zint,zbot,zblend,isoth,cs2int,lnrhoint,fac_cs)
!
!  Implement a polytropic profile in ss above zbot. If this routine is
!  called several times (for a piecewise polytropic atmosphere), on needs
!  to call it from top to bottom.
!
!  zint    -- height of (previous) interface, where cs2int and lnrhoint
!             are set
!  zbot    -- z at bottom of layer
!  zblend  -- smoothly blend (with width whcond) previous ss (for z>zblend)
!             with new profile (for z<zblend)
!  isoth   -- flag for isothermal stratification;
!  lnrhoin -- value of lnrho at the interface, i.e. at the zint on entry,
!             at the zbot on exit
!  cs2int  -- same for cs2
!
      use Sub, only: step
      use Gravity, only: gravz
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: stp
      real :: tmp,mpoly,zint,zbot,zblend,beta1,cs2int,lnrhoint
      integer :: isoth
!
      intent(in)    :: mpoly,zint,zbot,zblend,isoth
      real, intent(in), optional    :: fac_cs
      intent(inout) :: cs2int,lnrhoint,f
!
      stp = step(z,zblend,widthlnrho(1))
      do n=n1,n2; do m=m1,m2
! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
        if (isoth/=0.0) then ! isothermal layer
          beta1 = 0.0
          tmp = lnrhoint+gamma*gravz/cs2int*(z(n)-zint)
        else
          beta1 = gamma*gravz/(mpoly+1)
          tmp = 1.0 + beta1*(z(n)-zint)/cs2int
! Abort if args of log() are negative
          if ( (tmp<=0.0) .and. (z(n)<=zblend) ) &
            call fatal_error_local('polytropic_lnrho_z','Imaginary density values -- z_inf is too low')
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          tmp = lnrhoint + mpoly*log(tmp)
        endif
!
! smoothly blend the old value (above zblend) and the new one (below
! zblend) for the two regions:
!
        f(l1:l2,m,n,ilnrho) = stp(n)*f(l1:l2,m,n,ilnrho) + (1-stp(n))*tmp
!
      enddo; enddo
!
      if (isoth/=0.0) then
        lnrhoint = lnrhoint + gamma*gravz/cs2int*(zbot-zint)
      else
        lnrhoint = lnrhoint + mpoly*log(1 + beta1*(zbot-zint)/cs2int)
      endif
      if (isoth.ne.0 .and. present(fac_cs)) then
        cs2int = fac_cs**2*cs2int ! cs2 at layer interface (bottom)
      else
        cs2int = cs2int + beta1*(zbot-zint) ! cs2 at layer interface (bottom)
      endif
!
    endsubroutine polytropic_lnrho_z
!***********************************************************************
    subroutine polytropic_lnrho_disc(f,mpoly,zint,zbot,zblend,isoth,cs2int,lnrhoint)
!
!  Implement a polytropic profile in a disc. If this routine is
!  called several times (for a piecewise polytropic atmosphere), on needs
!  to call it from bottom (middle of disc) upwards.
!
!  zint    -- height of (previous) interface, where cs2int and lnrhoint
!             are set
!  zbot    -- z at top of layer (name analogous with polytropic_lnrho_z)
!  zblend  -- smoothly blend (with width whcond) previous ss (for z>zblend)
!             with new profile (for z<zblend)
!  isoth   -- flag for isothermal stratification;
!  lnrhoint -- value of lnrho at the interface, i.e. at the zint on entry,
!             at the zbot on exit
!  cs2int  -- same for cs2
!
!  24-jun-03/ulf:  coded
!
      use Sub, only: step
      use Gravity, only: gravz, nu_epicycle
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: stp
      real :: tmp,mpoly,zint,zbot,zblend,beta1,cs2int,lnrhoint,nu_epicycle2
      integer :: isoth
!
      do n=n1,n2; do m=m1,m2
! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
        nu_epicycle2 = nu_epicycle**2
        if (isoth/=0.0) then ! isothermal layer
          beta1 = 0.0
          tmp = gamma*gravz*nu_epicycle2/cs2int*(z(n)**2-zint**2)/2.
        else
          beta1 = gamma*gravz*nu_epicycle2/(mpoly+1)
          tmp = 1.0 + beta1*(z(n)**2-zint**2)/cs2int/2.
! Abort if args of log() are negative
          if ( (tmp<=0.0) .and. (z(n)<=zblend) ) &
            call fatal_error_local('polytropic_lnrho_disc','Imaginary density values -- z_inf is too low')
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          tmp = lnrhoint + mpoly*log(tmp)
        endif
!
! smoothly blend the old value (above zblend) and the new one (below
! zblend) for the two regions:
!
        stp = step(z,zblend,widthlnrho(1))
        f(l1:l2,m,n,ilnrho) = stp(n)*f(l1:l2,m,n,ilnrho) + (1-stp(n))*tmp
!
      enddo; enddo
!
      if (isoth/=0.0) then
        lnrhoint = lnrhoint + gamma*gravz*nu_epicycle2/cs2int*(zbot**2-zint**2)/2.
      else
        lnrhoint = lnrhoint + mpoly*log(1 + beta1*(zbot**2-zint**2)/cs2int/2.)
      endif
      cs2int = cs2int + beta1*(zbot**2-zint**2)/2.
!
    endsubroutine polytropic_lnrho_disc
!***********************************************************************
    subroutine shell_lnrho(f)
!
!  Initialize density based on specified radial profile in
!  a spherical shell
!
!  22-oct-03/dave -- coded
!  21-aug-08/dhruba -- added spherical coordinates
!
      use Gravity, only: g0,potential
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: pot, r_mn
      real :: beta1,lnrho_int,lnrho_ext,pot_int,pot_ext
!
      beta1=g0/(mpoly+1)*gamma/gamma_m1  ! gamma_m1/gamma=R_{*} (for cp=1)
!
      if (lspherical_coords) then
!     densities at shell boundaries
        lnrho_int=lnrho0+mpoly*log(1+beta1*(x(l2)/x(l1)-1.))
        lnrho_ext=lnrho0
!
! always inside the fluid shell
        do imn=1,nyz
          n=nn(imn)
          m=mm(imn)
          f(l1:l2-1,m,n,ilnrho)=lnrho0+mpoly*log(1+beta1*(x(l2)/x(l1:l2-1)-1.))
          f(l2,m,n,ilnrho)=lnrho_ext
        enddo
!
      elseif (lcylindrical_coords) then
        call fatal_error('shell_lnrho','not consistent with cylindrical coords')
      else
!     densities at shell boundaries
        lnrho_int=lnrho0+mpoly*log(1+beta1*(r_ext/r_int-1.))
        lnrho_ext=lnrho0
!
        do imn=1,nyz
          n=nn(imn)
          m=mm(imn)
!
          r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
!
        ! in the fluid shell
          where (r_mn < r_ext .AND. r_mn > r_int) f(l1:l2,m,n,ilnrho)=lnrho0+mpoly*log(1+beta1*(r_ext/r_mn-1.))
        ! outside the fluid shell
            if (initlnrho(1)=='geo-kws') then
              where (r_mn >= r_ext) f(l1:l2,m,n,ilnrho)=lnrho_ext
              where (r_mn <= r_int) f(l1:l2,m,n,ilnrho)=lnrho_int
            elseif (initlnrho(1)=='geo-kws-constant-T'.or.initlnrho(1)=='geo-benchmark') then
              call potential(R=r_int,POT=pot_int)
              call potential(R=r_ext,POT=pot_ext)
              call potential(RMN=r_mn,POT=pot)
! gamma/gamma_m1=1/R_{*} (for cp=1)
              where (r_mn >= r_ext) f(l1:l2,m,n,ilnrho)=lnrho_ext+(pot_ext-pot)*exp(-lnrho_ext/mpoly)*gamma/gamma_m1
              where (r_mn <= r_int) f(l1:l2,m,n,ilnrho)=lnrho_int+(pot_int-pot)*exp(-lnrho_int/mpoly)*gamma/gamma_m1
            endif
        enddo
      endif
!
    endsubroutine shell_lnrho
!***********************************************************************
    subroutine numerical_equilibrium(f)
!
!  sets gravity gg in order to achieve an numerical exact equilbrium
!  at t=0. This is only valid for the polytropic case, i.e.
!
!    (1/rho) grad(P) = cs20 (rho/rho0)^(gamma-2) grad(rho)
!
      use Sub, only: grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: lnrho,cs2
      real, dimension (nx,3) :: glnrho
      real, dimension (nx,3) :: gg_mn
      integer :: j
!
      if (ldensity_nolog) call not_implemented('numerical_equilibrium','for linear density')
!
      do m=m1,m2
        do n=n1,n2
!
          lnrho=f(l1:l2,m,n,ilnrho)
          cs2=cs20*exp(gamma_m1*(lnrho-lnrho0))
          call grad(f,ilnrho,glnrho)
          do j=1,3
            gg_mn(:,j)=cs2*glnrho(:,j)
          enddo
          f(l1:l2,m,n,iglobal_gg:iglobal_gg+2)=gg_mn
!
        enddo
      enddo
!
    endsubroutine numerical_equilibrium
!***********************************************************************
    subroutine pencil_criteria_density
!
!  All pencils that the Density module depends on are specified here.
!
!  19-11-04/anders: coded
!
      if (ldensity_nolog) lpenc_requested(i_rho)=.true.
      if (lcontinuity_gas) then
        if (lweno_transport) then
          lpenc_requested(i_transprho)=.true.
        else
          if (.not.lconservative) lpenc_requested(i_divu)=.true.
          if (ldensity_nolog) then
            if (lconservative) then
              if (lhydro) lpenc_requested(i_divss)=.true.
            else
              lpenc_requested(i_ugrho)=.true.
!            lpenc_requested(i_uglnrho)=.false.
            endif
          else
            if (lconservative) call fatal_error('pencil_criteria_density', 'must use ldensity_nolog=T')
            lpenc_requested(i_uglnrho)=.true.
!            lpenc_requested(i_ugrho)=.false.
          endif
        endif
      endif
      if (lSchur_3D3D1D) then
        lpenc_requested(i_divu)=.true.
        lpenc_requested(i_uglnrho)=.true.
      endif
      if (ldiff_shock) then
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        if (ldensity_nolog .or. ldiffusion_nolog) then
          lpenc_requested(i_grho)=.true.
          lpenc_requested(i_del2rho)=.true.
          if (ldiffusion_nolog) lpenc_requested(i_rho1)=.true.
        else
          lpenc_requested(i_glnrho)=.true.
          lpenc_requested(i_glnrho2)=.true.
          lpenc_requested(i_del2lnrho)=.true.
        endif
      endif
      if (ldiff_normal.or.ldiff_cspeed) then
        if (ldensity_nolog .or. ldiffusion_nolog) then
          lpenc_requested(i_del2rho)=.true.
          if (ldiffusion_nolog) lpenc_requested(i_rho1)=.true.
        else
          lpenc_requested(i_glnrho2)=.true.
          lpenc_requested(i_del2lnrho)=.true.
        endif
      endif
      if (ldiff_cspeed) lpenc_requested(i_TT)=.true.
      if ( ldiff_hyper3.or.ldiff_hyper3_strict) lpenc_requested(i_del6rho)=.true.
      if ((ldiff_hyper3.or.ldiff_hyper3_strict).and..not.ldensity_nolog) lpenc_requested(i_rho)=.true.
      if (ldiff_hyper3_polar.and..not.ldensity_nolog) lpenc_requested(i_rho1)=.true.
      if (ldiff_hyper3lnrho .or. ldiff_hyper3lnrho_strict) lpenc_requested(i_del6lnrho)=.true.
!
      if (lmass_source) then
        if ((mass_source_profile=='bump').or.(mass_source_profile=='sph-step-down')) &
            lpenc_requested(i_r_mn)=.true.
        if (mass_source_profile=='cylindric') lpenc_requested(i_rcyl_mn)=.true.
      endif
!
      if (lmassdiff_fix) then
        if ((lhydro.or.lentropy.or.ltemperature).and.ldensity_nolog) lpenc_requested(i_rho1)=.true.
        if (lhydro) lpenc_requested(i_uu)=.true.
        if (lentropy) lpenc_requested(i_cv)=.true.
        if (ltemperature.and.ltemperature_nolog) lpenc_requested(i_TT)=.true.
        if (lthermal_energy) lpenc_requested(i_u2) = .true.
      endif
!
      if (lfargo_advection) then
        lpenc_requested(i_uu_advec)=.true.
        if (ldensity_nolog) then
          lpenc_requested(i_grho)=.true.
          lpenc_requested(i_uuadvec_grho)=.true.
        else
          lpenc_requested(i_glnrho)=.true.
          lpenc_requested(i_uuadvec_glnrho)=.true.
        endif
      endif
!
      if (lreference_state) lpenc_requested(i_rho1) = .true.
      lpenc_diagnos2d(i_lnrho)=.true.
      lpenc_diagnos2d(i_rho)=.true.
!
!  Diagnostic pencils.
!
      if (idiag_rhom/=0 .or. idiag_rho2m/=0 .or.idiag_rho4m/=0 .or.idiag_rho6m/=0 .or. &
           idiag_rho12m/=0 .or. idiag_rhof2m/=0 .or. idiag_rhomy/=0 .or. &
           idiag_rhomx/=0 .or. idiag_rho2mx/=0 .or. idiag_rhomz/=0 .or. idiag_rho2mz/=0 .or. &
           idiag_rhomin/=0 .or.  idiag_rhomax/=0 .or. idiag_rhomxz/=0 .or. &
           idiag_totmass/=0 .or. idiag_mass/=0 .or. idiag_drho2m/=0 .or. idiag_rhorms/=0 .or. &
           idiag_inertiaxx/=0 .or. idiag_inertiayy/=0 .or. idiag_inertiazz/=0 .or. &
           idiag_inertiaxx_car/=0 .or. idiag_inertiayy_car/=0 .or. idiag_inertiazz_car/=0 .or. &
           idiag_drhom/=0 .or. idiag_rhomxmask/=0 .or. idiag_sigma/=0 .or. idiag_rhomzmask/=0 .or. &
           idiag_rhoupmz/=0 .or. idiag_rhodownmz/=0 .or. idiag_rho2upmz/=0 .or. &
           idiag_rho2downmz/=0 .or. idiag_rhof2mz/=0 .or. idiag_rhof2upmz/=0 .or. &
           idiag_rhof2downmz/=0 .or. idiag_sphmass/=0) &
           lpenc_diagnos(i_rho)=.true.
      if (idiag_rhomxy/=0 .or. idiag_rho2mxy/=0) lpenc_diagnos2d(i_rho)=.true.
!AB: idiag_rhof2mz, idiag_rhodownmz, etc, shouldn't be here, right?
!PJK: The up/down averages should actually appear above as well and they
!PJK: need the velocity too. idiag_rhof2mz does not, however.
      if (idiag_rhoupmz/=0 .or. idiag_rhodownmz/=0 .or. idiag_rho2upmz/=0 .or. &
           idiag_rho2downmz/=0 .or. idiag_rhof2upmz/=0 .or. idiag_rhof2downmz/=0) &
           lpenc_diagnos(i_uu)=.true.
      if (idiag_lnrho2m/=0.or.idiag_lnrhomin/=0.or.idiag_lnrhomax/=0.or.idiag_lnrhorms/=0) lpenc_diagnos(i_lnrho)=.true.
      if (idiag_ugrhom/=0) lpenc_diagnos(i_ugrho)=.true.
      if (idiag_uglnrhom/=0) lpenc_diagnos(i_uglnrho)=.true.
      if (idiag_uglnrhomz/=0) lpenc_diagnos(i_uglnrho)=.true.
      if (idiag_ugrhomz/=0) lpenc_diagnos(i_ugrho)=.true.
      if (idiag_uygzlnrhomz/=0 .or. idiag_uzgylnrhomz/=0) lpenc_diagnos(i_glnrho)=.true.
      if (idiag_grhomax/=0) lpenc_diagnos(i_grho)=.true.
      if (idiag_inertiaxx_car/=0 .or. idiag_inertiayy_car/=0 .or. &
          idiag_inertiazz_car/=0 .or. idiag_sphmass/=0) lpenc_diagnos(i_r_mn)=.true.
!
      if (lconservative) then
        lpenc_requested(i_lorentz)=.true.
        if (lhiggsless)  lpenc_requested(i_hless)=.true.
      endif
!
    endsubroutine pencil_criteria_density
!***********************************************************************
    subroutine pencil_interdep_density(lpencil_in)
!
!  Interdependency among pencils from the Density module is specified here.
!
!  19-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (ldensity_nolog) then
        if (lpencil_in(i_rho1)) lpencil_in(i_rho)=.true.
      else
        if (lpencil_in(i_rho)) lpencil_in(i_rho1)=.true.
      endif
      if (lpencil_in(i_grho)) then
        if (.not.ldensity_nolog) lpencil_in(i_rho)=.true.
      endif
      if (lpencil_in(i_glnrho)) then
        if (ldensity_nolog) lpencil_in(i_rho1)=.true.
      endif
      if (lpencil_in(i_uglnrho)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_ugrho)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_grho)=.true.
      endif
      if (lpencil_in(i_glnrho2)) lpencil_in(i_glnrho)=.true.
      if (lpencil_in(i_sglnrho)) then
        lpencil_in(i_sij)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_uij5glnrho)) then
        lpencil_in(i_uij5)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_del2lnrho)) then
        if (ldensity_nolog) then
          lpencil_in(i_rho1)=.true.
          lpencil_in(i_del2rho)=.true.
          lpencil_in(i_glnrho2)=.true.
        endif
      endif
      if (lpencil_in(i_ekin)) then
        lpencil_in(i_rho)=.true.
        lpencil_in(i_u2)=.true.
      endif
!
    endsubroutine pencil_interdep_density
!***********************************************************************
    subroutine calc_pencils_density_pnc(f,p,lpenc_loc)
!
!  Calculate Density pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  19-11-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(:), intent(IN) :: lpenc_loc
      intent(in) :: f
      intent(inout) :: p
!
!  Differentiate between log density and linear density.
!
      if (ldensity_nolog) then
        call calc_pencils_linear_density_pnc(f,p,lpenc_loc)
      else
        call calc_pencils_log_density_pnc(f,p,lpenc_loc)
      endif
! ekin
      if (lpenc_loc(i_ekin)) then
        if (lconservative) then
          p%ekin=fourthird*p%rho*p%lorentz*p%u2
        else
          p%ekin=0.5*p%rho*p%u2
        endif
      endif
!
!  Dummy pencils.
!
      if (lpenc_loc(i_rhos1)) call not_implemented('calc_pencils_density_pnc','rhos1')
      if (lpenc_loc(i_glnrhos)) call not_implemented('calc_pencils_density_pnc','glnrhos')
!
    endsubroutine calc_pencils_density_pnc
!***********************************************************************
    subroutine calc_pencils_density_std(f,p)
!
! Envelope adjusting calc_pencils_density_pnc to the standard use with
! lpenc_loc=lpencil
!
! 21-sep-13/MR    : coded
!
      real, dimension (mx,my,mz,mfarray),intent(IN)   :: f
      type (pencil_case),                intent(INOUT):: p
!
      call calc_pencils_density_pnc(f,p,lpencil)
!
      endsubroutine calc_pencils_density_std
!***********************************************************************
    subroutine calc_pencils_linear_density_std(f,p)
!
! Envelope adjusting calc_pencils_density_pnc to the standard use with
! lpenc_loc=lpencil
!
! 21-sep-13/MR    : coded
!
      real, dimension (mx,my,mz,mfarray),intent(IN)   :: f
      type (pencil_case),                intent(INOUT):: p
!
      call calc_pencils_linear_density_pnc(f,p,lpencil)
!
      endsubroutine calc_pencils_linear_density_std
!***********************************************************************
    subroutine calc_pencils_linear_density_pnc(f,p,lpenc_loc)
!
!  Calculate Density pencils for linear density.
!  Most basic pencils should come first, as others may depend on them.
!
!  13-05-10/dhruba: stolen parts of earlier calc_pencils_density.
!  21-jan-15/MR: changes for use of reference state.
!  22-jan-15/MR: removed unneeded get_shared_variable('reference_state,...
!  19-jan-15/MR: adapted weno transport for use with reference state.
!                suppressed weno for log density
!
      use WENO_transport
      use Sub, only: div,grad,dot,dot2,u_dot_grad,del2,del6,multmv,g2ij,dot_mn,h_dot_grad, &
                     del6_strict,calc_del6_for_upwind

      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(:), intent(IN) :: lpenc_loc
      intent(in) :: f
      intent(inout) :: p
      real, dimension(nx) :: tmp
      integer :: i
!
! rho
!
      p%rho=f(l1:l2,m,n,irho)
      if (lreference_state) p%rho=p%rho+reference_state(:,iref_rho)
! rho1
      if (lpenc_loc(i_rho1)) p%rho1=1.0/p%rho
! lnrho
      if (lpenc_loc(i_lnrho)) p%lnrho=log(p%rho)
! glnrho and grho
      if (lpenc_loc(i_glnrho).or.lpenc_loc(i_grho)) then
!
        call grad(f,irho,p%grho)
        if (lreference_state) p%grho(:,1)=p%grho(:,1)+reference_state(:,iref_grho)
!
        if (lpenc_loc(i_glnrho)) then
          do i=1,3
            p%glnrho(:,i)=p%rho1*p%grho(:,i)
          enddo
        endif
      endif
! uglnrho
      if (lpenc_loc(i_uglnrho)) call fatal_error('calc_pencils_density', &
          'uglnrho not available for linear mass density')   ! Why not implementing it?
! ugrho
      if (lpenc_loc(i_ugrho)) call u_dot_grad(f,ilnrho,p%grho,p%uu,p%ugrho,UPWIND=lupw_rho)
! glnrho2
      if (lpenc_loc(i_glnrho2)) call dot2(p%glnrho,p%glnrho2)
! del2rho
      if (lpenc_loc(i_del2rho)) then
        call del2(f,irho,p%del2rho)
        if (lreference_state) p%del2rho=p%del2rho+reference_state(:,iref_d2rho)
      endif
! del2lnrho
      if (lpenc_loc(i_del2lnrho)) p%del2lnrho=p%rho1*p%del2rho-p%glnrho2
! del6rho
      if (lpenc_loc(i_del6rho)) then
        if (ldiff_hyper3) then
          call del6(f,irho,p%del6rho)
          if (lreference_state) p%del6rho=p%del6rho+reference_state(:,iref_d6rho)
        elseif (ldiff_hyper3_strict) then
          call del6_strict(f,irho,p%del6rho)
        endif
      endif
! del6lnrho
      if (lpenc_loc(i_del6lnrho)) then
        if (ldiff_hyper3lnrho) then
          call not_implemented('calc_pencils_density','del6lnrho for linear mass density')
        elseif (ldiff_hyper3lnrho_strict) then
          call not_implemented('calc_pencils_density','del6lnrho_strict for linear mass density')
        endif
      endif
! hlnrho
      if (lpenc_loc(i_hlnrho)) call not_implemented('calc_pencils_density','hlnrho linear mass density')
! sglnrho
      if (lpenc_loc(i_sglnrho)) call multmv(p%sij,p%glnrho,p%sglnrho)
! uij5glnrho
      if (lpenc_loc(i_uij5glnrho)) call multmv(p%uij5,p%glnrho,p%uij5glnrho)
! transprho
      if (lpenc_loc(i_transprho)) then
        if (lreference_state) then
          call weno_transp(f,m,n,irho,-1,iux,iuy,iuz,p%transprho,dx_1,dy_1,dz_1, &
                           reference_state(:,iref_rho))
        else
          call weno_transp(f,m,n,irho,-1,iux,iuy,iuz,p%transprho,dx_1,dy_1,dz_1)
        endif
      endif
!
      if (lpenc_loc(i_uuadvec_grho)) then
        call h_dot_grad(p%uu_advec,p%grho,p%uuadvec_grho)
        if (lupw_rho) then
          call calc_del6_for_upwind(f,irho,p%uu_advec,tmp)
          p%uuadvec_grho = p%uuadvec_grho - tmp
        endif
      endif
!
! divS, needed for relativistic calculations
!
      if (lpenc_loc(i_divss)) call div(f,iux,p%divss)
!
      if (lconservative) then
        if (lhiggsless) then
          where(real(t) < p%hless) p%rho=p%rho-eps_hless
          p%rho=p%rho/(fourthird*p%lorentz*(1.-.25/p%lorentz))
        else
          p%rho=p%rho/(fourthird*p%lorentz*(1.-.25/p%lorentz))
        endif
      endif
!
    endsubroutine calc_pencils_linear_density_pnc
!***********************************************************************
    subroutine calc_pencils_log_density_pnc(f,p,lpenc_loc)
!
!  Calculate Density pencils for logarithmic density.
!  Most basic pencils should come first, as others may depend on them.
!
!  13-05-10/dhruba: stolen parts of earlier calc_pencils_density
!
      use WENO_transport
      use General, only: notanumber
      use Sub, only: grad,dot,dot2,u_dot_grad,del2,del6,del6_strict,multmv,g2ij,dot_mn,h_dot_grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(:) :: lpenc_loc
!
      intent(in) :: f
      intent(inout) :: p
!
      integer :: i
!
! lnrho
      p%lnrho=f(l1:l2,m,n,ilnrho)
! rho1
      if (lpenc_loc(i_rho1)) p%rho1=exp(-f(l1:l2,m,n,ilnrho))
! rho
      if (lpenc_loc(i_rho)) p%rho=1.0/p%rho1
! glnrho and grho
      if (lpenc_loc(i_glnrho).or.lpenc_loc(i_grho)) then
        call grad(f,ilnrho,p%glnrho)
        if (notanumber(p%glnrho)) then
          print*, 'density:iproc,it,m,n=', iproc_world,it,m,n
          !print*, 'it,m,n=', it,m,n
          !write(iproc+10,*) "density:f(:,m,n,ilnrho)=",f(:,m,n,ilnrho)
          call fatal_error_local('calc_pencils_density','NaNs in p%glnrho')
        endif
        if (lpenc_loc(i_grho)) then
          do i=1,3
            p%grho(:,i)=p%rho*p%glnrho(:,i)
          enddo
        endif
      endif
! uglnrho
      if (lpenc_loc(i_uglnrho)) then
        if (lupw_lnrho) then
          call u_dot_grad(f,ilnrho,p%glnrho,p%uu,p%uglnrho,UPWIND=lupw_lnrho)
        else
          call dot(p%uu,p%glnrho,p%uglnrho)
        endif
      endif
! ugrho
      if (lpenc_loc(i_ugrho)) call not_implemented('calc_pencils_density', &
          'ugrho for logarithmic mass density')
! glnrho2
      if (lpenc_loc(i_glnrho2)) call dot2(p%glnrho,p%glnrho2)
! del2rho
      if (lpenc_loc(i_del2rho)) call not_implemented('calc_pencils_density', &
          'del2rho for logarithmic mass density')
! del2lnrho
      if (lpenc_loc(i_del2lnrho)) call del2(f,ilnrho,p%del2lnrho)
! del6rho
      if (lpenc_loc(i_del6rho)) call not_implemented('calc_pencils_density', &
          'del6rho for logarithmic mass density')
! del6lnrho
      if (lpenc_loc(i_del6lnrho)) then
        if (ldiff_hyper3lnrho) then
          call del6(f,ilnrho,p%del6lnrho)
        elseif (ldiff_hyper3lnrho_strict) then
          call del6_strict(f,ilnrho,p%del6lnrho)
        endif
      endif

! hlnrho
      if (lpenc_loc(i_hlnrho))  call g2ij(f,ilnrho,p%hlnrho)
! sglnrho
      if (lpenc_loc(i_sglnrho)) call multmv(p%sij,p%glnrho,p%sglnrho)
! uij5glnrho
      if (lpenc_loc(i_uij5glnrho)) call multmv(p%uij5,p%glnrho,p%uij5glnrho)
!
      if (lpenc_loc(i_uuadvec_glnrho)) call h_dot_grad(p%uu_advec,p%glnrho,p%uuadvec_glnrho)
!
    endsubroutine calc_pencils_log_density_pnc
!***********************************************************************
    subroutine density_after_boundary(f)
!
!  Actions to take after boundary conditions are set.
!
!   2-apr-08/anders: coded
!
      !use Sub, only: div, grad, dot_mn, finalize_aver
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      if (lSchur_3D3D1D) then
!
!  For Schur flows, initialize the averages of the RHS of the density equation
!
        Schur_dlnrho_RHS_xyzaver=0.
        Schur_dlnrho_RHS_xyaver_z=0.
        Schur_dlnrho_RHS_zaver_xy=0.
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine density_after_boundary
!***********************************************************************
    subroutine dlnrho_dt(f,df,p)
!
!  Continuity equation.
!  Calculate dlnrho/dt = - u.gradlnrho - divu
!
!   7-jun-02/axel: incorporated from subroutine pde
!  21-oct-15/MR: changes for slope-limited diffusion
!   4-aug-17/axel: implemented terms for ultrarelativistic EoS
!  25-may-18/fred: updated mass diffusion correction and set default
!                  not default for hyperdiff, but correction applies
!                  to all fdiff if lmassdiff_fix=T
!  03-apr-20/joern: restructured and fixed slope-limited diffusion
!
      use Deriv, only: der6
      use Special, only: special_calc_density
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  :: p
      intent(inout) :: df,f
!
      real, dimension (nx) :: fdiff
      real, dimension (nx) :: tmp
      real, dimension (nx,3) :: tmpv
      real, dimension (nx) :: density_rhs,advec_hypermesh_rho
      integer :: j
      logical :: ldt_up
!
!  Identify module and boundary conditions.
!
      call timing('dlnrho_dt','entered',mnloop=.true.)
      if (headtt.or.ldebug) print*,'dlnrho_dt: SOLVE'
      if (headtt) call identify_bcs('lnrho',ilnrho)
!
      if (lSchur_3D3D1D) then
!
!  Accumulatively calculate the RHS of Schur flow equations, but only finalize after the mn loop.
!
        density_rhs=p%uglnrho+p%divu
!
        Schur_dlnrho_RHS_xyaver_z(n-nghost) = Schur_dlnrho_RHS_xyaver_z(n-nghost)+sum(density_rhs)/nxygrid
        Schur_dlnrho_RHS_zaver_xy(:,m-nghost) = Schur_dlnrho_RHS_zaver_xy(:,m-nghost)+density_rhs/nzgrid
        Schur_dlnrho_RHS_xyzaver = Schur_dlnrho_RHS_xyzaver+sum(density_rhs)/nwgrid
      else
!
!  Continuity equation.
!
      if (lcontinuity_gas) then
        if (.not. lweno_transport .and. .not. lffree .and. .not. lreduced_sound_speed .and. &
            ieos_profile=='nothing' .and. .not. lfargo_advection) then
!
!  Evolution of rho; set and initiate density_rhs
!
          if (ldensity_nolog) then
            if (lconservative) then
              density_rhs=-p%divss
            else
              density_rhs=-p%rho*p%divu
              if (ladvection_density) density_rhs = density_rhs - p%ugrho
              if (lrelativistic_eos) density_rhs=fourthird*density_rhs
            endif
!
!  Evolution of lnrho
!
          else
            density_rhs= - p%divu
            if (ladvection_density) density_rhs = density_rhs - p%uglnrho
!
!  The following few lines only enter without lconservative,
!  and also only without ldensity_nolog, but with lrelativistic_eos.
!
            if (lrelativistic_eos.and..not.lconservative) then
              if (lhydro) then
                call multvs(p%uu,density_rhs,tmpv)
                df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-onethird*tmpv
              endif
              density_rhs=fourthird*density_rhs
            endif
          endif

        else
          density_rhs=0.
        endif
!
!  WENO transport.
!
        if (lweno_transport) density_rhs= density_rhs - p%transprho
!
!  Choice of vertical profile in front of density evolution.
!  Default is off. This is useful to simulate outer halo regions.
!  There is an additional option of doing this by obeying mass
!  conservation, which is not currently the default.
!
        if (ieos_profile=='surface_z') then
          if (ldensity_nolog) then
            density_rhs= density_rhs - profz_eos(n)*(p%ugrho + p%rho*p%divu)
            if (ldensity_profile_masscons) density_rhs = density_rhs-dprofz_eos(n)*p%rho*p%uu(:,3)
          else
            density_rhs= density_rhs - profz_eos(n)*(p%uglnrho + p%divu)
            if (ldensity_profile_masscons) density_rhs = density_rhs -dprofz_eos(n)*p%uu(:,3)
          endif
        endif
!
!  If we are solving the force-free equation in parts of our domain.
!
        if (lffree) then
          if (ldensity_nolog) then
            density_rhs= density_rhs - profx_ffree*profy_ffree(m)*profz_ffree(n)*(p%ugrho + p%rho*p%divu)
            if (ldensity_profile_masscons) density_rhs=density_rhs - p%rho*( dprofx_ffree   *p%uu(:,1) &
                                                                            +dprofy_ffree(m)*p%uu(:,2) &
                                                                            +dprofz_ffree(n)*p%uu(:,3))
          else
            density_rhs= density_rhs - profx_ffree*(profy_ffree(m)*profz_ffree(n))*(p%uglnrho + p%divu)
            if (ldensity_profile_masscons) density_rhs=density_rhs-dprofx_ffree   *p%uu(:,1) &
                                                                  -dprofy_ffree(m)*p%uu(:,2) &
                                                                  -dprofz_ffree(n)*p%uu(:,3)
          endif
        endif
!
!  Use reduced sound speed according to Rempel (2005), ApJ, 622, 1320;
!  see also Hotta et al. (2012), A&A, 539, A30.
!
        if (lreduced_sound_speed) then
          if (ldensity_nolog) then
            density_rhs = density_rhs - reduce_cs2_profx*reduce_cs2_profz(n)*(p%ugrho + p%rho*p%divu)
          else
            density_rhs = density_rhs - reduce_cs2_profx*reduce_cs2_profz(n)*(p%uglnrho + p%divu)
          endif
        endif
!
        if (lfargo_advection) then
          if (ldensity_nolog) then
            density_rhs = density_rhs - p%uuadvec_grho   - p%rho*p%divu
          else
            density_rhs = density_rhs - p%uuadvec_glnrho - p%divu
          endif
        endif
!
!  Add the continuity equation terms to the RHS of the density df.
!
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + density_rhs
!
      endif
!
!  Mass sources and sinks.
!
      if (lupdate_mass_source) call mass_source(f,df,p)
!
!  Mass diffusion.
!
      diffus_diffrho=0.; diffus_diffrho3=0.
      fdiff=0.0
!
      ldt_up = lfirst.and.ldt
!
      if (ldiff_normal) then  ! Normal diffusion operator
        if (ldensity_nolog) then
          fdiff = fdiff + diffrho*p%del2rho
        else
          if (ldiffusion_nolog) then
            fdiff = fdiff + diffrho*p%rho1*p%del2rho
          else
            fdiff = fdiff + diffrho*(p%del2lnrho+p%glnrho2)
          endif
        endif
        if (ldt_up) diffus_diffrho=diffus_diffrho+diffrho
        if (headtt) print*,'dlnrho_dt: diffrho=', diffrho
      endif
!
      if (ldiff_cspeed) then  ! Normal diffusion operator
        if (ldensity_nolog) then
          fdiff = fdiff + diffrho*p%TT**diff_cspeed*p%del2rho
        else
          if (ldiffusion_nolog) then
            fdiff = fdiff + diffrho*p%TT**diff_cspeed*p%rho1*p%del2rho
          else
            fdiff = fdiff + diffrho*p%TT**diff_cspeed*(p%del2lnrho+p%glnrho2)
          endif
        endif
        if (lfirst.and.ldt) diffus_diffrho=diffus_diffrho+diffrho
        if (headtt) print*,'dlnrho_dt: diffrho=', diffrho
      endif
!
!  Shock diffusion
!
      if (ldiff_shock) then
        if (ldensity_nolog) then
          call dot_mn(p%gshock,p%grho,tmp)
          fdiff = fdiff + diffrho_shock * (p%shock * p%del2rho + tmp)
        else
          if (ldiffusion_nolog) then
            call dot_mn(p%gshock,p%grho,tmp)
            fdiff = fdiff + p%rho1 * diffrho_shock * (p%shock * p%del2rho + tmp)
          else
            call dot_mn(p%gshock,p%glnrho,tmp)
            fdiff = fdiff + diffrho_shock * (p%shock * (p%del2lnrho + p%glnrho2) + tmp)
!
!  Counteract the shock diffusion of the mean stratification. Must set
!  lwrite_stratification=T in start.in for this to work.
!
            if (lanti_shockdiffusion) &
              fdiff = fdiff - diffrho_shock * (p%shock*(del2lnrho_glnrho2_init_z(n) + &
                      2*(p%glnrho(:,3)-dlnrhodz_init_z(n))*dlnrhodz_init_z(n)) + &
                      p%gshock(:,3)*dlnrhodz_init_z(n) )
          endif
        endif
        if (ldt_up) diffus_diffrho=diffus_diffrho+diffrho_shock*p%shock
        if (headtt) print*,'dlnrho_dt: diffrho_shock=', diffrho_shock
      endif
!
!   Slope limited diffusion for density
!
      if (ldensity_slope_limited.and.llast) then
        if (ldensity_nolog) then
          call calc_slope_diff_flux(f,irho,p,h_sld_dens,nlf_sld_dens,tmp,div_sld_dens)
          fdiff=fdiff+tmp
        else
          call calc_slope_diff_flux(f,ilnrho,p,h_sld_dens,nlf_sld_dens,tmp,div_sld_dens)
          fdiff=fdiff+tmp*p%rho1
        endif
      endif
!
!  Interface for your personal subroutines calls
!
      if (lspecial) call special_calc_density(f,df,p)
!
!  Improve energy and momentum conservation by compensating for mass diffusion
!
      if (lmassdiff_fix.and..not.lconservative) then
        if (ldensity_nolog) then
          tmp = fdiff*p%rho1
        else
          tmp = fdiff
        endif

        if (lhydro) then
          forall(j = iux:iuz) df(l1:l2,m,n,j) = df(l1:l2,m,n,j) - p%uu(:,j-iuu+1) * tmp
        endif

        if (lentropy.and.(.not.pretend_lnTT)) then
          df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%cv*tmp
        elseif (lentropy.and.pretend_lnTT) then
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - tmp
        elseif (ltemperature.and.(.not. ltemperature_nolog)) then
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - tmp
        elseif (ltemperature.and.ltemperature_nolog) then
          df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - tmp*p%TT
        elseif (lthermal_energy) then
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + 0.5 * fdiff * p%u2
        endif
      endif
!
!  Hyper diffusion.
!
      if (ldiff_hyper3.or.ldiff_hyper3_strict) then
        if (ldensity_nolog) then
          fdiff = fdiff + diffrho_hyper3*p%del6rho
        else
          if (ldiffusion_nolog) fdiff = fdiff + diffrho_hyper3*p%rho1*p%del6rho
        endif
        if (ldt_up) diffus_diffrho3=diffus_diffrho3+diffrho_hyper3
        if (headtt) print*,'dlnrho_dt: diffrho_hyper3=', diffrho_hyper3
      endif
!
      if (ldiff_hyper3_polar) then
        do j=1,3
          !for ldensity_nolog it is doing del6lnrho, as it should for a simple hyperdiffusion
          if (ldensity_nolog) then
            call der6(f,irho,tmp,j,IGNOREDX=.true.)
          else
            call der6(f,ilnrho,tmp,j,IGNOREDX=.true.)
          endif
          fdiff = fdiff + diffrho_hyper3*pi4_1*tmp*dline_1(:,j)**2
!
!AB: wouldn't the following line make more sense?
!AB: But then the meaning of diffrho_hyper3 changes, so we
!AB: should not just add it to the previous diffus_diffrho3.
          !fdiff = fdiff + diffrho_hyper3*pi5_1*tmp*dline_1(:,j)**5
        enddo
        if (ldt_up) diffus_diffrho3=diffus_diffrho3+diffrho_hyper3*pi4_1*dxmin_pencil**4
        if (headtt) print*,'dlnrho_dt: diffrho_hyper3=', diffrho_hyper3
      endif
!
!  Mesh-hyperdiffusion. The parameter diffrho_hyper3_mesh has currently the
!  unit of velocity and should later be changed to maxadvec (to be checked as
!  diagnostics). A good value for diffrho_hyper3_mesh is 5 (which is currently
!  the default). This method should also work in all coordinate systems.
!
!WL: Why should diffrho_hyper3_mesh be 5? An explanation should go in the manual.
!
!
      if (ldiff_hyper3_mesh) then
        do j=1,3
          ! for rho this is d6rho, for lnrho this is d6lnrho
           call der6(f,ilnrho,tmp,j,IGNOREDX=.true.)
          if (ldynamical_diffusion) then
            fdiff = fdiff + diffrho_hyper3_mesh * tmp * dline_1(:,j)
          else
            fdiff = fdiff + diffrho_hyper3_mesh*pi5_1/60.*tmp*dline_1(:,j)
!AB: at the meeting in Toulouse we discussed that it should rather
!AB: involve the actual time step dt, but its value is not known during
!AB: the first execution of a 3rd order substep. I keep this comment,
!AB: so we don't forget about this issue. A plauble thing to try is to
!AB: check if dt==0, then leave fdiff unchanged (hoping that this is
!AB: the default value when it hasn't been set yet).
            !fdiff = fdiff + diffrho_hyper3_mesh*pi5_1/60.*tmp/dt
          endif
        enddo
        if (ldt_up) then
          if (ldynamical_diffusion) then
            diffus_diffrho3 = diffus_diffrho3 + diffrho_hyper3_mesh
            advec_hypermesh_rho=0.
          else
            advec_hypermesh_rho=diffrho_hyper3_mesh*pi5_1*sqrt(dxyz_2)
          endif
          advec2_hypermesh=advec2_hypermesh+advec_hypermesh_rho**2
        endif
        if (headtt) print*,'dlnrho_dt: diffrho_hyper3_mesh=', diffrho_hyper3_mesh
      endif
!
      if (ldiff_hyper3_aniso) then
        call del6fj(f,diffrho_hyper3_aniso,ilnrho,tmp)
        fdiff = fdiff + tmp
        if (lsubtract_init_stratification) then
          call del6fj(f,diffrho_hyper3_aniso,iglobal_lnrho0,tmp)
          fdiff = fdiff - tmp
        endif
!  Must divide by dxyz_6 here, because it is multiplied on later.
        if (ldt_up) diffus_diffrho3=diffus_diffrho3 + &
                 (diffrho_hyper3_aniso(1)*dline_1(:,1)**6 + &
                  diffrho_hyper3_aniso(2)*dline_1(:,2)**6 + &
                  diffrho_hyper3_aniso(3)*dline_1(:,3)**6)/dxyz_6
        if (headtt) print*,'dlnrho_dt: diffrho_hyper3=(Dx,Dy,Dz)=',diffrho_hyper3_aniso
      endif
!
      if (ldiff_hyper3lnrho .or. ldiff_hyper3lnrho_strict) then
        if (.not. ldensity_nolog) fdiff = fdiff + diffrho_hyper3*p%del6lnrho
        if (ldt_up) diffus_diffrho3=diffus_diffrho3+diffrho_hyper3
        if (headtt .and. ldiff_hyper3lnrho ) print*,'dlnrho_dt: diffrho_hyper3=', diffrho_hyper3
        if (headtt .and. ldiff_hyper3lnrho_strict ) print*, &
           'dlnrho_dt: diffrho_hyper3_strict=', diffrho_hyper3
      endif
!
!  Add diffusion term to continuity equation
!
      if (ldensity_nolog) then
        df(l1:l2,m,n,irho)   = df(l1:l2,m,n,irho)   + fdiff
      else
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + fdiff
      endif
!
!  Multiply diffusion coefficient by Nyquist scale.
!
      if (ldt_up) then
        diffus_diffrho = diffus_diffrho*dxyz_2
        if (ldynamical_diffusion .and. ldiff_hyper3_mesh) then
          diffus_diffrho3 = diffus_diffrho3 * sum(abs(dline_1),2)
        else
          diffus_diffrho3 = diffus_diffrho3*dxyz_6
        endif
        if (headtt.or.ldebug) then
          print*,'dlnrho_dt: max(diffus_diffrho ) =', maxval(diffus_diffrho)
          print*,'dlnrho_dt: max(diffus_diffrho3) =', maxval(diffus_diffrho3)
        endif
        maxdiffus=max(maxdiffus,diffus_diffrho)
        maxdiffus3=max(maxdiffus3,diffus_diffrho3)
      endif
!
!  Apply border profile
!
      if (lborder_profiles) call set_border_density(f,df,p)
!
      endif   !  if (lSchur_3D3D1D) then ... else
!
      call timing('dlnrho_dt','before l2davgfirst',mnloop=.true.)
      call calc_diagnostics_density(f,p)
      call timing('dlnrho_dt','finished',mnloop=.true.)
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine calc_diagnostics_density(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p
!
      call calc_2d_diagnostics_density(p)
      call calc_1d_diagnostics_density(f,p)
      call calc_0d_diagnostics_density(f,p)
!
    endsubroutine calc_diagnostics_density
!***********************************************************************
    subroutine calc_2d_diagnostics_density(p)
!
!  2d-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      use Diagnostics

      type(pencil_case) :: p

      if (l2davgfirst) then
        call phisum_mn_name_rz(p%lnrho,idiag_lnrhomphi)
        call phisum_mn_name_rz(p%rho,idiag_rhomphi)
        call ysum_mn_name_xz(p%rho,idiag_rhomxz)
        call zsum_mn_name_xy(p%rho,idiag_rhomxy)
        if (idiag_rho2mxy/=0) call zsum_mn_name_xy(p%rho**2,idiag_rho2mxy)
        call zsum_mn_name_xy(p%rho,idiag_sigma,lint=.true.)
      endif

    endsubroutine calc_2d_diagnostics_density
!***********************************************************************
    subroutine calc_1d_diagnostics_density(f,p)
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p
!
      real, dimension (nx) :: uzmask
!
      if (l1davgfirst) then
!
        call phizsum_mn_name_r(p%rho,idiag_rhomr)
        call xysum_mn_name_z(p%rho,idiag_rhomz)
        if (idiag_rho2mz/=0) call xysum_mn_name_z(p%rho**2,idiag_rho2mz)
        call xysum_mn_name_z(p%glnrho(:,3),idiag_gzlnrhomz)
        call xysum_mn_name_z(p%uglnrho,idiag_uglnrhomz)
        call xysum_mn_name_z(p%ugrho,idiag_ugrhomz)
        if (idiag_uygzlnrhomz/=0) call xysum_mn_name_z(p%uu(:,2)*p%glnrho(:,3),idiag_uygzlnrhomz)
        if (idiag_uzgylnrhomz/=0) call xysum_mn_name_z(p%uu(:,3)*p%glnrho(:,2),idiag_uzgylnrhomz)
        call yzsum_mn_name_x(p%rho,idiag_rhomx)
        if (idiag_rho2mx/=0) call yzsum_mn_name_x(p%rho**2,idiag_rho2mx)
        call xzsum_mn_name_y(p%rho,idiag_rhomy)
        if (lrho_flucz_as_aux) then
          if (idiag_rhof2mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,irho_flucz)**2,idiag_rhof2mz)
        endif
        if (idiag_rhoupmz/=0 .or. idiag_rho2upmz/=0 .or. idiag_rhof2upmz/=0) then
          where (p%uu(:,3) > 0.)
            uzmask = p%uu(:,3)/abs(p%uu(:,3))
          elsewhere
            uzmask=0.
          endwhere
          call xysum_mn_name_z(uzmask*p%rho,idiag_rhoupmz)
          call xysum_mn_name_z(uzmask*p%rho**2,idiag_rho2upmz)
          if (lrho_flucz_as_aux) then
            if (idiag_rhof2upmz/=0) call xysum_mn_name_z(uzmask*f(l1:l2,m,n,irho_flucz)**2,idiag_rhof2upmz)
          endif
        endif
        if (idiag_rhodownmz/=0 .or. idiag_rho2downmz/=0 .or. idiag_rhof2downmz/=0) then
          where (p%uu(:,3) < 0.)
            uzmask = -p%uu(:,3)/abs(p%uu(:,3))
          elsewhere
            uzmask=0.
          endwhere
          call xysum_mn_name_z(uzmask*p%rho,idiag_rhodownmz)
          call xysum_mn_name_z(uzmask*p%rho**2,idiag_rho2downmz)
          if (lrho_flucz_as_aux) then
            if (idiag_rhof2downmz/=0) call xysum_mn_name_z(uzmask*f(l1:l2,m,n,irho_flucz)**2,idiag_rhof2downmz)
          endif
        endif
      endif

    endsubroutine calc_1d_diagnostics_density
!***********************************************************************
    subroutine calc_0d_diagnostics_density(f,p)
!
!  Calculate density diagnostics
!
!  21-dec-18/felipe: added inertiaxx - inertiazz for Applegate mechanism
!
      use Diagnostics

      use Sub,only: dot2
      use General
!
      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p
!
      real, dimension(nx), parameter :: unitpencil=1.
      real, dimension(nx) :: tmp, rmask
!
!  The inertiaxx - inertiazz terms are needed for computing the star's
!  quadrupole moment, relevant for the Applegate mechanism; see
!  Navarrete et al. (MNRAS 2019).
!
      if (ldiagnos) then
        call sum_mn_name(p%rho,idiag_rhom)
        if (idiag_rhomxmask/=0)call sum_mn_name(p%rho*xmask_den,idiag_rhomxmask)
        if (idiag_rhomzmask/=0)call sum_mn_name(p%rho*zmask_den(n-n1+1),idiag_rhomzmask)
        call sum_mn_name(p%rho,idiag_totmass,lint=.true.)
        call integrate_mn_name(p%rho,idiag_mass)
        if (idiag_inertiaxx/=0)call integrate_mn_name(p%rho*x(l1:l2)**2*sin(y(m))**2*cos(z(n))**2,idiag_inertiaxx)
        if (idiag_inertiayy/=0)call integrate_mn_name(p%rho*x(l1:l2)**2*sin(y(m))**2*sin(z(n))**2,idiag_inertiayy)
        if (idiag_inertiazz/=0)call integrate_mn_name(p%rho*x(l1:l2)**2*cos(y(m))**2,idiag_inertiazz)
        if (idiag_inertiaxx_car/=0 .or. idiag_inertiayy_car/=0 .or. &
            idiag_inertiazz_car/=0 .or. idiag_sphmass/=0) then
          where (p%r_mn <= radius_diag)
            rmask = 1.
          elsewhere
            rmask = 0.
          endwhere
          if (idiag_inertiaxx_car/=0) call integrate_mn_name(rmask*p%rho*x(l1:l2)**2,idiag_inertiaxx_car)
          if (idiag_inertiayy_car/=0) call integrate_mn_name(rmask*p%rho*y(m)**2,idiag_inertiayy_car)
          if (idiag_inertiazz_car/=0) call integrate_mn_name(rmask*p%rho*z(n)**2,idiag_inertiazz_car)
          if (idiag_sphmass/=0) call integrate_mn_name(rmask*p%rho,idiag_sphmass)
        endif
        call integrate_mn_name(unitpencil,idiag_vol)
        if (idiag_rhomin/=0)   call max_mn_name(-p%rho,idiag_rhomin,lneg=.true.)
        call max_mn_name(p%rho,idiag_rhomax)
        if (idiag_lnrhomin/=0) call max_mn_name(-p%lnrho,idiag_lnrhomin,lneg=.true.)
        call max_mn_name(p%lnrho,idiag_lnrhomax)
        if (idiag_rho2m/=0)    call sum_mn_name(p%rho**2,idiag_rho2m)
        if (idiag_rho4m/=0)    call sum_mn_name(p%rho**4,idiag_rho4m)
        if (idiag_rho6m/=0)    call sum_mn_name(p%rho**6,idiag_rho6m)
        if (idiag_rho12m/=0)    call sum_mn_name(p%rho**12,idiag_rho12m)
        if (idiag_rhof2m/=0.and.lrho_flucz_as_aux) call sum_mn_name(f(l1:l2,m,n,irho_flucz)**2,idiag_rhof2m)
        if (idiag_rhorms/=0)   call sum_mn_name(p%rho**2,idiag_rhorms,lsqrt=.true.)
        if (idiag_lnrhorms/=0) call sum_mn_name(p%lnrho**2,idiag_lnrhorms,lsqrt=.true.)
        if (idiag_lnrho2m/=0)  call sum_mn_name(p%lnrho**2,idiag_lnrho2m)
        if (idiag_drho2m/=0)   call sum_mn_name((p%rho-rho0)**2,idiag_drho2m)
        if (idiag_drhom/=0)    call sum_mn_name(p%rho-rho0,idiag_drhom)
        call sum_mn_name(p%ugrho,idiag_ugrhom)
        call sum_mn_name(p%uglnrho,idiag_uglnrhom)
        if (.not.lgpu) then
          if (idiag_dtd/=0) call max_mn_name(diffus_diffrho/cdtv,idiag_dtd,l_dt=.true.)
          if (idiag_dtd3/=0) call max_mn_name(diffus_diffrho3/cdtv3,idiag_dtd3,l_dt=.true.)
        endif
        if (idiag_grhomax/=0) then
          call dot2(p%grho,tmp); tmp=sqrt(tmp)
          call max_mn_name(tmp,idiag_grhomax)
        endif
      endif

    endsubroutine calc_0d_diagnostics_density
!***********************************************************************
    subroutine split_update_density(f)
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine split_update_density
!***********************************************************************
    subroutine density_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dtsub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine density_after_timestep
!***********************************************************************
    subroutine set_border_density(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the lnrho variable.
!
!  28-jul-06/wlad: coded
!
      use BorderProfiles, only: border_driving,set_border_initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: f_target
      type (pencil_case)  :: p
!
      select case (borderlnrho)
!
      case ('zero','0')
        if (ldensity_nolog) then
          f_target=0.
        else
          f_target=1.
        endif
!
      case ('constant')
        if (ldensity_nolog) then
          f_target=rho_const
        else
          f_target=lnrho_const
        endif
!
      case ('initial-condition')
        call set_border_initcond(f,ilnrho,f_target)
!
      case ('nothing')
        return
!
      endselect
!
      call border_driving(f,df,p,f_target,ilnrho)
!
    endsubroutine set_border_density
!***********************************************************************
!  Here comes a collection of different density stratification routines
!***********************************************************************
    subroutine isothermal_density(f)
!
!  Isothermal stratification (for lnrho and ss).
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!  Sound speed (and hence Temperature), and density (at infinity) are
!  initialised to their respective reference values:
!           sound speed: cs^2_0            from start.in
!           density: rho0 = exp(lnrho0)
!
!   8-jul-02/axel: incorporated/adapted from init_lnrho
!  11-jul-02/axel: fixed sign; should be tmp=gamma*pot/cs20
!  02-apr-03/tony: made entropy explicit rather than using tmp/-gamma
!  11-jun-03/tony: moved entropy initialisation to separate routine
!                  to allow isothermal condition for arbitrary density
!  26-sep-10/axel: added lisothermal_fixed_Hrho to allow for isothermal density
!                  stratification to remain unchanged when gamma is changed.
!
      use Gravity
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx) :: pot,tmp
!
!  Stratification depends on the gravity potential.
!
      if (lroot) print*, 'isothermal_density: isothermal stratification'
      if (gamma/=1.0) then
        if ((.not. lentropy) .and. (.not. ltemperature)) call fatal_error('isothermal_density', &
            'for gamma/=1.0, you need entropy or temperature!');
      endif
!
      do n=n1,n2
        do m=m1,m2
          call potential(x(l1:l2),y(m),z(n),pot=pot)
          if (lisothermal_fixed_Hrho) then
            tmp=-pot/cs20
          else
            tmp=-gamma*pot/cs20
          endif
          f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + lnrho0 + tmp
          if (lentropy.and..not.pretend_lnTT) then
!
!  note: the initial condition is always produced for lnrho
!
            f(l1:l2,m,n,iss) = f(l1:l2,m,n,iss) - (1/cp1)*gamma_m1*(f(l1:l2,m,n,ilnrho)-lnrho0)/gamma
          elseif (lentropy.and.pretend_lnTT) then
            f(l1:l2,m,n,ilnTT)=log(cs20*cp1/gamma_m1)
          elseif (ltemperature.and..not.ltemperature_nolog) then
            f(l1:l2,m,n,ilnTT)=log(cs20*cp1/gamma_m1)
          elseif (ltemperature.and.ltemperature_nolog) then
            f(l1:l2,m,n,iTT)=cs20*cp1/gamma_m1
          endif
        enddo
      enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      cs2bot=cs20
      cs2top=cs20
!
    endsubroutine isothermal_density
!***********************************************************************
    subroutine stratification_tsallis(f)
!
!  Tsallis stratification (for lnrho and ss).
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!  17-apr-13/axel+illa: adapted from isothermal
!
      use Gravity
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx) :: pot,tmp
      real :: pot1,tmp1
!
!  Stratification depends on the gravity potential;
!
      if (lroot) print*, 'stratification_tsallis: Tsallis stratification'
!
      do n=n1,n2
        do m=m1,m2
          call potential(x(l1:l2),y(m),z(n),pot=pot)
          if (gamma==1.) then
            tmp=lnrho0-pot/cs20
          else
            tmp=lnrho0+alog(1.+(gamma-1.)*(-pot/cs20))/(gamma-1.)
          endif
          f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)+tmp
!
!  Add corresponding profile to the other thermodynamic variable.
!
          if (lentropy.and..not.pretend_lnTT) then
            f(l1:l2,m,n,iss) = f(l1:l2,m,n,iss) + &
                (1/cp1)*(f(l1:l2,m,n,ilnrho)-lnrho0)*(ggamma/gamma-1.)
          elseif (lentropy.and.pretend_lnTT) then
            call not_implemented("stratification_tsallis","for pretend_lnTT=T")
          elseif (ltemperature.and..not.ltemperature_nolog) then
            call not_implemented("stratification_tsallis","for ltemperature=T")
          elseif (ltemperature.and.ltemperature_nolog) then
            call not_implemented("stratification_tsallis","for ltemperature_nolog=T")
          endif
        enddo
      enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      if (gamma==1.) then
        cs2bot=cs20
        cs2top=cs20
      else
        call potential(z=xyz0(3),pot=pot1)
        tmp1=lnrho0+alog(1.+(gamma-1.)*(-pot1/cs20))/(gamma-1.)
        cs2bot=cs20*exp(gamma_m1*tmp1)
!
        call potential(z=xyz1(3),pot=pot1)
        tmp1=lnrho0+alog(1.+(gamma-1.)*(-pot1/cs20))/(gamma-1.)
        cs2top=cs20*exp(gamma_m1*tmp1)
      endif
!
    endsubroutine stratification_tsallis
!***********************************************************************
    subroutine polytropic_simple(f)
!
!  Polytropic stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!
!  To maintain continuity with respect to the isothermal case,
!  one may want to specify cs20 (=1), and so zinfty is calculated from that.
!  On the other hand, for polytropic atmospheres it may be more
!  physical to specify zinfty (=1), ie the top of the polytropic atmosphere.
!  This is done if zinfty is different from 0.
!
!   8-jul-02/axel: incorporated/adapted from init_lnrho
!
      use Gravity, only: gravz_profile,gravz,zinfty,zref,zgrav,potential,nu_epicycle
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: pot,dlncs2,r_mn
      real :: ztop,zbot,zref2,pot_ext,lnrho_ref,ptop,pbot
!
!  identifier
!
      if (lroot) print*,'polytropic_simple: mpoly=',mpoly
!
!  The following is specific only to cases with gravity in the z direction
!  zref is calculated such that rho=rho0 and cs2=cs20 at z=zref.
!  Note: gravz is normally negative!
!
      if (lgravz) then
        if (gravz_profile=='const') then
          if (lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zref=zinfty-(mpoly+1.)*cs20/(-gamma*gravz)
        elseif (gravz_profile=='const_zero') then
          if (lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zref=zinfty-(mpoly+1.)*cs20/(-gamma*gravz)
        elseif (gravz_profile=='linear') then
          if (lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zref2=zinfty**2-(mpoly+1.)*cs20/(0.5*gamma*nu_epicycle**2)
          if (zref2<0) then
            if (lroot) print*,'polytropic_simple: zref**2<0 is not ok'
            zref2=0. !(and see what happens)
          endif
          zref=sqrt(zref2)
        else
          if (lroot) print*,'polytropic_simple: zref not prepared!'
        endif
        if (lroot) print*,'polytropic_simple: zref=',zref
!
!  check whether zinfty lies outside the domain (otherwise density
!  would vanish within the domain). At the moment we are not properly
!  testing the lower boundary on the case of a disc (commented out).
!
        ztop=xyz0(3)+Lxyz(3)
        zbot=xyz0(3)
        !-- if (zinfty<min(ztop,zgrav) .or. (-zinfty)>min(zbot,zgrav)) then
        if (zinfty<min(ztop,zgrav)) then
          if (lroot) print*,'polytropic_simple: domain too big; zinfty=',zinfty
          !call stop_it( &
          !         'polytropic_simply: rho and cs2 will vanish within domain')
        endif
      endif
!
!  polytropic sphere with isothermal exterior
!  calculate potential at the stellar surface, pot_ext
!
      if (lgravr) then
        call potential(R=r_ext,POT=pot_ext)
        cs2top=-gamma/(mpoly+1.)*pot_ext
        lnrho_ref=mpoly*log(cs2top)-(mpoly+1.)
        if (lroot) print*,'polytropic_simple: pot_ext=',pot_ext
        do n=n1,n2; do m=m1,m2
          r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          call potential(x(l1:l2),y(m),z(n),pot=pot)
!
!  density
!  these formulae assume lnrho0=0 and cs0=1
!
          where (r_mn > r_ext)
            f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)+lnrho_ref-gamma*pot/cs2top
          elsewhere
            dlncs2=log(-gamma*pot/((mpoly+1.)*cs20))
            f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)+lnrho0+mpoly*dlncs2
          endwhere
!
!  entropy
!
          if (lentropy) then
            where (r_mn > r_ext)
              f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss) &
                -(1.-1./gamma)*f(l1:l2,m,n,ilnrho)+log(cs2top)/gamma
            elsewhere
              dlncs2=log(-gamma*pot/((mpoly+1.)*cs20))
              f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+mpoly*(ggamma/gamma-1.)*dlncs2
            endwhere
          endif
        enddo; enddo
      else
!
!  cartesian case with gravity in the z direction
!
        do n=n1,n2; do m=m1,m2
          call potential(x(l1:l2),y(m),z(n),pot=pot)
          dlncs2=log(-gamma*pot/((mpoly+1.)*cs20))
          f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)+lnrho0+mpoly*dlncs2
          if (lentropy) f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+mpoly*(ggamma/gamma-1.)*dlncs2
!         if (ltemperature) f(l1:l2,m,n,ilnTT)=dlncs2-log(gamma_m1)
          if (ltemperature) f(l1:l2,m,n,ilnTT)=log(-gamma*pot/(mpoly+1.)/gamma_m1)
        enddo; enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  In spherical geometry, ztop is z at the outer edge of the box,
!  so this calculation still makes sense.
!
        call potential(xyz0(1),xyz0(2),ztop,pot=ptop)
        cs2top=-gamma*ptop/(mpoly+1.)
!
!  In spherical geometry ztop should never be used.
!  Even in slab geometry ztop is not normally used.
!
        call potential(xyz0(1),xyz0(2),zbot,pot=pbot)
        cs2bot=-gamma*pbot/(mpoly+1.)
      endif
!
    endsubroutine polytropic_simple
!***********************************************************************
    subroutine mass_source(f,df,p)
!
!  Add (isothermal) mass sources and sinks.
!
!  28-apr-2005/axel: coded
!  14-may-2019/axel: changed xblob -> xblob(1) for now
!
      use General, only: random_number_wrapper
      use Sub, only: step
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (2) :: fran
      real :: tmp
      real, dimension (nx) :: dlnrhodt, pdamp, fprofile, radius2
!
      if (ldebug) print*,'mass_source: cs20,cs0=',cs20,cs0
!
!  Choose between different possibilities.
!
      select case (mass_source_profile)
        case ('nothing')
          call not_implemented('mass_source','mass source with no profile')
        case ('exponential')
          dlnrhodt=mass_source_Mdot
        case('bump')
          dlnrhodt=(mass_source_Mdot/fnorm)*exp(-.5*(p%r_mn/mass_source_sigma)**2)
        case('bump2')
          dlnrhodt=fprofile_z(n-n1+1)
        case('bumpr')
          radius2=(x(l1:l2)-xblob(1))**2+(y(m)-yblob(1))**2+(z(n)-zblob(1))**2
          fprofile=(mass_source_Mdot/fnorm)*exp(-.5*radius2/mass_source_sigma**2)
          if (lmass_source_random) then
            call random_number_wrapper(fran)
            tmp=sqrt(-2*log(fran(1)))*sin(2*pi*fran(2))
            dlnrhodt=fprofile*cos(mass_source_omega*t)*tmp
          else
            dlnrhodt=fprofile*cos(mass_source_omega*t)
          endif
        case('bumpx','sph-step-down')
          dlnrhodt=fprofile_x
        case('const' )
          dlnrhodt=-mass_source_tau1*(f(l1:l2,m,n,ilnrho)-lnrho0)
        case('cylindric')
!
!  Cylindrical profile for inner cylinder.
!
          pdamp=1.-step(p%rcyl_mn,r_int,wdamp) ! inner damping profile
          dlnrhodt=-damplnrho_int*pdamp*(f(l1:l2,m,n,ilnrho)-lnrho_int)
!
!  Cylindrical profile for outer cylinder.
!
          pdamp=step(p%rcyl_mn,r_ext,wdamp) ! outer damping profile
          dlnrhodt=dlnrhodt-damplnrho_ext*pdamp*(f(l1:l2,m,n,ilnrho)-lnrho_ext)
!
! default to catch unknown values
!
        case default
          call fatal_error('mass_source','no such mass_source_profile: '//trim(mass_source_profile))
        endselect
!
!  Add mass source.
!
      if (ldensity_nolog) then
        df(l1:l2,m,n,irho)=df(l1:l2,m,n,irho)+p%rho*dlnrhodt
      else
        df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+dlnrhodt
      endif
!
!  Change entropy to keep temperature constant.
!
      if (lentropy) df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+(gamma1-1.0)*dlnrhodt
!
    endsubroutine mass_source
!***********************************************************************
    subroutine impose_density_floor(f)
!
!  Impose a minimum density by setting all lower densities to the minimum
!  value (density_floor). Useful for debugging purposes.
!
!  13-aug-07/anders: implemented.
!  10-feb-15/MR: adaptations for reference state
!  17-may-23/hongzhe: implemented non-uniform density floor (for lnrho only)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real :: density_floor_local
      integer :: i, j, k
!
!  Impose the density floor.
!
      if (density_floor>0.) then
!
        if (ldensity_nolog) then
          if (.not.lreference_state) &     ! else?
            where (f(:,:,:,irho)<density_floor) f(:,:,:,irho)=density_floor
        else
          select case (density_floor_profile)
          case ('uniform')
            where (f(:,:,:,ilnrho)<density_floor_log) f(:,:,:,ilnrho)=density_floor_log
          case('sphr_powerlaw')
            do i=1,mx
              density_floor_local = density_floor_log + density_floor_exp*log(x(i))
              where (f(i,:,:,ilnrho)<density_floor_local) f(i,:,:,ilnrho) = density_floor_local
            enddo
          case default
            call fatal_error_local('impose_density_floor','no such density_floor_profile implemented')
          endselect
        endif
      else
!
!  Trap any negative density if no density floor is set.
!
        if (lcheck_negative_density .and. ldensity_nolog) then
          if (any(f(l1:l2,m1:m2,n1:n2,irho) <= 0.)) then
            do k = n1, n2
              do j = m1, m2
                do i = l1, l2
                  if (f(i,j,k,irho) <= 0.) print 10, f(i,j,k,irho), x(i), y(j), z(k)
                  10 format (1x,'rho = ',es13.6,' at x = ',es13.6,', y = ',es13.6,', z = ',es13.6)
                enddo
              enddo
            enddo
            call fatal_error_local('impose_density_floor', 'negative density detected')
          endif
        endif
!
      endif
!
    endsubroutine impose_density_floor
!***********************************************************************
    subroutine read_density_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=density_init_pars, IOSTAT=iostat)
!
    endsubroutine read_density_init_pars
!***********************************************************************
    subroutine write_density_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=density_init_pars)
!
    endsubroutine write_density_init_pars
!***********************************************************************
    subroutine read_density_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=density_run_pars, IOSTAT=iostat)
!
    endsubroutine read_density_run_pars
!***********************************************************************
    subroutine write_density_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=density_run_pars)
!
    endsubroutine write_density_run_pars
!***********************************************************************
    subroutine rprint_density(lreset,lwrite)
!
!  Reads and registers print parameters relevant for continuity equation.
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
      use General, only: loptest
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez, inamexy, inamexz, irz, inamer
!
!  Reset everything in case of reset.
!  (This needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rhom=0; idiag_rho2m=0; idiag_rho4m=0; idiag_rho6m=0; idiag_rho12m=0;
        idiag_rhof2m=0; idiag_lnrho2m=0
        idiag_drho2m=0; idiag_drhom=0; idiag_rhorms=0; idiag_lnrhorms=0
        idiag_ugrhom=0; idiag_ugrhomz=0; idiag_uglnrhom=0
        idiag_rhomin=0; idiag_rhomax=0; idiag_dtd=0; idiag_dtd3=0
        idiag_lnrhomin=0; idiag_lnrhomax=0;
        idiag_lnrhomphi=0; idiag_rhomphi=0
        idiag_rhomz=0; idiag_rho2mz=0; idiag_rhomy=0; idiag_rhomx=0; idiag_rho2mx=0
        idiag_rhoupmz=0; idiag_rhodownmz=0; idiag_rho2upmz=0; idiag_rho2downmz=0
        idiag_rhof2mz=0; idiag_rhof2upmz=0; idiag_rhof2downmz=0
        idiag_gzlnrhomz=0; idiag_uglnrhomz=0; idiag_uygzlnrhomz=0; idiag_uzgylnrhomz=0
        idiag_rhomxy=0; idiag_rhomr=0; idiag_totmass=0; idiag_mass=0; idiag_vol=0
        idiag_rhomxz=0; idiag_grhomax=0; idiag_inertiaxx=0; idiag_inertiayy=0
        idiag_inertiazz=0; idiag_inertiaxx_car=0; idiag_inertiayy_car=0
        idiag_inertiazz_car=0; idiag_rhomxmask=0; idiag_rhomzmask=0
        idiag_sigma=0; idiag_rho2mxy=0; idiag_sphmass=0
      endif
!
!  iname runs through all possible names that may be listed in print.in.
!
      if (lroot.and.ip<14) print*,'rprint_density: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rhom',idiag_rhom)
        call parse_name(iname,cname(iname),cform(iname),'rhomxmask',idiag_rhomxmask)
        call parse_name(iname,cname(iname),cform(iname),'rhomzmask',idiag_rhomzmask)
        call parse_name(iname,cname(iname),cform(iname),'rho2m',idiag_rho2m)
        call parse_name(iname,cname(iname),cform(iname),'rho4m',idiag_rho4m)
        call parse_name(iname,cname(iname),cform(iname),'rho6m',idiag_rho6m)
        call parse_name(iname,cname(iname),cform(iname),'rho12m',idiag_rho12m)
        call parse_name(iname,cname(iname),cform(iname),'rhof2m',idiag_rhof2m)
        call parse_name(iname,cname(iname),cform(iname),'rhorms',idiag_rhorms)
        call parse_name(iname,cname(iname),cform(iname),'lnrhorms',idiag_lnrhorms)
        call parse_name(iname,cname(iname),cform(iname),'drho2m',idiag_drho2m)
        call parse_name(iname,cname(iname),cform(iname),'drhom',idiag_drhom)
        call parse_name(iname,cname(iname),cform(iname),'rhomin',idiag_rhomin)
        call parse_name(iname,cname(iname),cform(iname),'rhomax',idiag_rhomax)
        call parse_name(iname,cname(iname),cform(iname),'lnrhomin',idiag_lnrhomin)
        call parse_name(iname,cname(iname),cform(iname),'lnrhomax',idiag_lnrhomax)
        call parse_name(iname,cname(iname),cform(iname),'lnrho2m',idiag_lnrho2m)
        call parse_name(iname,cname(iname),cform(iname),'ugrhom',idiag_ugrhom)
        call parse_name(iname,cname(iname),cform(iname),'uglnrhom',idiag_uglnrhom)
        call parse_name(iname,cname(iname),cform(iname),'dtd',idiag_dtd)
        call parse_name(iname,cname(iname),cform(iname),'dtd3',idiag_dtd3)
        call parse_name(iname,cname(iname),cform(iname),'totmass',idiag_totmass)
        call parse_name(iname,cname(iname),cform(iname),'mass',idiag_mass)
        call parse_name(iname,cname(iname),cform(iname),'sphmass',idiag_sphmass)
        call parse_name(iname,cname(iname),cform(iname),'inertiaxx',idiag_inertiaxx)
        call parse_name(iname,cname(iname),cform(iname),'inertiayy',idiag_inertiayy)
        call parse_name(iname,cname(iname),cform(iname),'inertiazz',idiag_inertiazz)
        call parse_name(iname,cname(iname),cform(iname),'inertiaxx_car',idiag_inertiaxx_car)
        call parse_name(iname,cname(iname),cform(iname),'inertiayy_car',idiag_inertiayy_car)
        call parse_name(iname,cname(iname),cform(iname),'inertiazz_car',idiag_inertiazz_car)
        call parse_name(iname,cname(iname),cform(iname),'vol',idiag_vol)
        call parse_name(iname,cname(iname),cform(iname),'grhomax',idiag_grhomax)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhomz',idiag_rhomz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhoupmz',idiag_rhoupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhodownmz',idiag_rhodownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rho2upmz',idiag_rho2upmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rho2downmz',idiag_rho2downmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhof2mz',idiag_rhof2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhof2upmz',idiag_rhof2upmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhof2downmz',idiag_rhof2downmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rho2mz',idiag_rho2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gzlnrhomz',idiag_gzlnrhomz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uglnrhomz',idiag_uglnrhomz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ugrhomz',idiag_ugrhomz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uygzlnrhomz',idiag_uygzlnrhomz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzgylnrhomz',idiag_uzgylnrhomz)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhomy', idiag_rhomy)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rhomx', idiag_rhomx)
        call parse_name(inamex, cnamex(inamex), cformx(inamex), 'rho2mx', idiag_rho2mx)
      enddo
!
!  Check for those quantities for which we want phiz-averages.
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'rhomr',idiag_rhomr)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'rhomxz',idiag_rhomxz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhomxy',idiag_rhomxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rho2mxy',idiag_rho2mxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'sigma',idiag_sigma)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'lnrhomphi', idiag_lnrhomphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rhomphi',idiag_rhomphi)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='rho'.or.cnamev=='lnrho') cformv='DEFINED'
      endif
!
!  Write column where which density variable is stored.
!
      if (loptest(lwrite)) then
        if (ldensity_nolog) then
          call farray_index_append('ilnrho',0)
        else
          call farray_index_append('irho',0)
        endif
      endif
!
    endsubroutine rprint_density
!***********************************************************************
    subroutine init_hydrostatic_r (f)
!
      use Gravity, only: potential,lnumerical_equilibrium
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: pot
      real :: pot0
      real, dimension (nx) :: r_mn
!
      intent(inout) :: f
!
      if (lgravr) then
         if (lroot) print*, 'init_hydrostatic_r: radial density stratification (assumes s=const)'
!
         do n=n1,n2; do m=m1,m2
            r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
            call potential(RMN=r_mn,POT=pot)
            call potential(R=r_ref,POT=pot0)
!  MEMOPT/AJ: commented out, since we do not use global potential anymore.
!            call output(trim(directory)//'/pot.dat',pot,1)
!
!  rho0, cs0, pot0 are the values at r=r_ref
!
            if (gamma/=1.0) then  ! isentropic
               f(l1:l2,m,n,ilnrho) = lnrho0 + log(1 - gamma_m1*(pot-pot0)/cs20) / gamma_m1
            else                  ! isothermal
               f(l1:l2,m,n,ilnrho) = lnrho0 - (pot-pot0)/cs20
            endif
         enddo; enddo
!
!  The following sets gravity gg in order to achieve numerical
!  exact equilibrium at t=0.
!
         if (lnumerical_equilibrium) call numerical_equilibrium(f)
      endif
!
    endsubroutine init_hydrostatic_r
!***********************************************************************
    subroutine init_sph_isoth (f)
!
!  Initialize isothermal sphere
!
!  14-may-10/dhruba: coded
!
      use EquationOfState, only: eoscalc,ilnrho_TT

      real, dimension (mx,my,mz,mfarray) :: f
      real :: haut
      real, dimension (nx) :: lnrho,TT,ss
!
      intent(inout) :: f
!
      if (lgravr) then
        if (lroot) print*, 'init_sph_isoth: isothermal sphere'
        haut=cs20/gamma
        TT=spread(cs20/gamma_m1,1,nx)
        do n=n1,n2
          do m=m1,m2
            r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
            f(l1:l2,m,n,ilnrho)=lnrho0-r_mn/haut
            lnrho=f(l1:l2,m,n,ilnrho)
            call eoscalc(ilnrho_TT,lnrho,TT,ss=ss)
            f(l1:l2,m,n,iss)=ss
          enddo
        enddo
      endif
!
    endsubroutine init_sph_isoth
!***********************************************************************
    subroutine get_slices_density(f,slices)
!
!  Write slices for animation of Density variables.
!
!  26-jul-06/tony: coded
!  10-feb-15/MR: adaptations for reference state
!
      use Slices_methods

      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
      character(LEN=labellen) :: name
!
!  Loop over slices
!
      name=trim(slices%name)
      if (name=='rho' .or. name=='lnrho') then
!
        call assign_slices_scal(slices,f,ilnrho)

        if (lfullvar_in_slices) &      ! implies ldensity_nolog=T
          call addto_slices(slices,reference_state(:,iref_rho))
!
        select case (trim(slices%name))
!
!  Density.
!
          case ('rho')
            if (.not.ldensity_nolog) call process_slices(slices,exp2d)
!
!  Logarithmic density.
!
          case ('lnrho')
            if (ldensity_nolog) then
!
              if (lreference_state.and..not.lfullvar_in_slices) call process_slices(slices,abs2d)
              call process_slices(slices,log2d)
!
            endif
        endselect
!
      endif
!
    endsubroutine get_slices_density
!***********************************************************************
    subroutine get_slices_pressure(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_pressure
!***********************************************************************
    subroutine get_init_average_density(f,init_average_density)
!
!  10-dec-09/piyali: added to pass initial average density
!
    real, dimension (mx,my,mz,mfarray):: f
    real:: init_average_density
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(init_average_density)
!
    endsubroutine get_init_average_density
!***********************************************************************
    subroutine density_after_mn(f,df,mass_per_proc)
!
!  20-apr-23/hongzhe: adapted from anelastic_after_mn
!
      use Sub, only: finalize_aver
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (1) :: mass_per_proc
!
      if (lSchur_3D3D1D) then
!
        call finalize_aver(nprocxy,12,Schur_dlnrho_RHS_xyaver_z)
        call finalize_aver(nprocz,3,Schur_dlnrho_RHS_zaver_xy)
        call finalize_aver(ncpus,123,Schur_dlnrho_RHS_xyzaver)
!
        do n=n1,n2
          df(l1:l2,m1:m2,n,ilnrho) = df(l1:l2,m1:m2,n,ilnrho)+Schur_dlnrho_RHS_xyzaver(1)  &
                                                             -Schur_dlnrho_RHS_xyaver_z(n-nghost) &
                                                             -Schur_dlnrho_RHS_zaver_xy
        enddo
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(mass_per_proc)
!
    endsubroutine density_after_mn
!***********************************************************************
    subroutine dynamical_diffusion(uc)
!
!  Dynamically set mass diffusion coefficient given fixed mesh Reynolds number.
!
!  27-jul-11/ccyang: coded
!
!  Input Argument
!      uc
!          Characteristic velocity of the system.
!
      real, intent(in) :: uc
!
!  Hyper-diffusion coefficient
!
      if (diffrho_hyper3 /= 0.0) diffrho_hyper3 = pi5_1 * uc * dxmax**5/re_mesh
      if (diffrho_hyper3_mesh /= 0.0) diffrho_hyper3_mesh = pi5_1 * uc/re_mesh/sqrt(real(dimensionality))
!
    endsubroutine dynamical_diffusion
!***********************************************************************
    subroutine boussinesq(f)
!
!  23-mar-2012/dintrans: coded
!  dummy routine for the Boussinesq approximation
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine boussinesq
!***********************************************************************
    subroutine update_reference_state
!
    endsubroutine update_reference_state
!***********************************************************************
    subroutine read_reference_state
!
!  Read reference state from a file
!
!  23-jan-2015/pete: aped from read_hcond
!
      use Mpicomm, only: mpibcast_real_arr
!
      integer, parameter :: ntotal=nx*nprocx
      real, dimension(ntotal) :: tmp1, tmp2, tmp3, tmp4
      real :: var1, var2, var3, var4
      logical :: exist
      integer :: stat, nn
!
!  Read reference_state and write into an array.
!  If file is not found in run directory, search under trim(directory).
!
      if (lroot ) then
        inquire(file='reference_state.dat',exist=exist)
        if (exist) then
          open(36,file='reference_state.dat')
        else
          inquire(file=trim(directory)//'/reference_state.ascii',exist=exist)
          if (exist) then
            open(36,file=trim(directory)//'/reference_state.ascii')
          else
            call fatal_error('read_reference_state','no input file')
          endif
        endif
!
!  Read profiles.
!
        do nn=1,ntotal
          read(36,*,iostat=stat) var1,var2,var3,var4
          if (stat<0) exit
          if (ip<5) print*,'rho, grho, del2rho, ss: ',var1,var2,var3,var4
          tmp1(nn)=var1
          tmp2(nn)=var2
          tmp3(nn)=var3
          tmp4(nn)=var4
        enddo
        close(36)
!
      endif
!
      call mpibcast_real_arr(tmp1, ntotal)
      call mpibcast_real_arr(tmp2, ntotal)
      call mpibcast_real_arr(tmp3, ntotal)
      call mpibcast_real_arr(tmp4, ntotal)
!
!  Assuming no ghost zones in cooling_profile.dat
!
      do nn=l1,l2
        reference_state(nn-nghost,iref_rho)=tmp1(ipx*nx+nn-nghost)
        reference_state(nn-nghost,iref_grho)=tmp2(ipx*nx+nn-nghost)
        reference_state(nn-nghost,iref_d2rho)=tmp3(ipx*nx+nn-nghost)
        reference_state(nn-nghost,iref_s)=tmp4(ipx*nx+nn-nghost)
      enddo
!
    endsubroutine read_reference_state
!***********************************************************************
    function mean_density(f)
!
!  Calculate mean density from f-array. With lreference_state=T it is the mean
!  density deviation.
!
!  3-mar-14/MR: coded
!
      use Mpicomm, only: mpiallreduce_sum
!
      real :: mean_density
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      integer :: n,m
      real :: tmp
!
      mean_density=0.
!
      do n=n1,n2
        tmp=0.
        do m=m1,m2
          if (ldensity_nolog) then
            tmp=tmp+sum(f(l1:l2,m,n,irho)*dVol_x(l1:l2))*dVol_y(m)
          else
            tmp=tmp+sum(exp(f(l1:l2,m,n,ilnrho))*dVol_x(l1:l2))*dVol_y(m)
          endif
        enddo
        mean_density=mean_density+tmp*dVol_z(n)
      enddo
!
      mean_density = mean_density/box_volume
!
      if (ncpus>1) then
        call mpiallreduce_sum(mean_density,tmp)
        mean_density=tmp
      endif
!
    endfunction mean_density
!***********************************************************************
    subroutine impose_density_ceiling(f)
!
!  Impose a maximum (log) density by setting all higher (log) densities to the maximum
!  value (density_ceiling).
!
!  3-mar-2017/MR: implemented.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  Impose the density ceiling
!
      if (density_ceiling>0.) then
!
        if (ldensity_nolog) then
          if (.not.lreference_state) f(:,:,:,irho) = min(f(:,:,:,irho),density_ceiling)
        else
         f(:,:,:,ilnrho) = min(f(:,:,:,ilnrho),density_ceiling_log)
        endif
      endif
!
    endsubroutine impose_density_ceiling
!***********************************************************************s
    subroutine write_z_stratification(f)

      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      if (lwrite_stratification) then
        open(19,file=trim(directory_dist)//'/stratification.dat')
        write(19,*) f(l1,m1,:,ilnrho)
        close(19)
      endif

    endsubroutine write_z_stratification
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=0
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    endsubroutine pushpars2c
!***********************************************************************
endmodule Density
