! $Id$
!
! MODULE_DOC: This module takes care of most of the things related to velocity.
! MODULE_DOC: Pressure, for example, is added in the energy (entropy) module.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lhydro = .true.
! CPARAM logical, parameter :: lhydro_kinematic = .false.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED divu; oo(3); o2; ou; oxu2; oxu(3); u2; uij(3,3); uu(3); curlo(3)
! PENCILS PROVIDED sij(3,3); sij2; uij5(3,3); ugu(3); ugu2; oij(3,3)
! PENCILS PROVIDED d2uidxj(3,3), uijk(3,3,3); ogu(3)
! PENCILS PROVIDED u3u21; u1u32; u2u13; del2u(3); del4u(3); del6u(3)
! PENCILS PROVIDED u2u31; u3u12; u1u23
! PENCILS PROVIDED graddivu(3); del6u_bulk(3); grad5divu(3)
! PENCILS PROVIDED rhougu(3); der6u(3); transpurho(3)
! PENCILS PROVIDED divu0; u0ij(3,3); uu0(3)
! PENCILS PROVIDED uu_advec(3); uuadvec_guu(3)
! PENCILS PROVIDED del6u_strict(3); del4graddivu(3); uu_sph(3)
! PENCILS PROVIDED der6u_res(3,3)
! PENCILS PROVIDED lorentz_gamma2; lorentz_gamma; ss_rel2; ss_rel(3)
! PENCILS PROVIDED ss_rel_ij(3,3); ss_rel_factor; divss_rel
!***************************************************************
!
module Hydro
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Viscosity, only: calc_viscous_force
!
  implicit none
!
  include 'hydro.h'
!
!  Slice precalculation buffers.
!
  real, target, dimension (:,:,:), allocatable :: oo_xy,oo_xy2,oo_xy3,oo_xy4
  real, target, dimension (:,:,:), allocatable :: oo_xz,oo_yz,oo_xz2
  real, target, dimension (:,:,:), allocatable :: uu_sph_xy,uu_sph_xy2,uu_sph_xy3,uu_sph_xy4
  real, target, dimension (:,:,:), allocatable :: uu_sph_xz,uu_sph_yz,uu_sph_xz2
  real, target, dimension (:,:), allocatable :: divu_xy,u2_xy,o2_xy,mach_xy
  real, target, dimension (:,:), allocatable :: divu_xy2,u2_xy2,o2_xy2,mach_xy2
  real, target, dimension (:,:), allocatable :: divu_xy3,divu_xy4,u2_xy3,u2_xy4,mach_xy4
  real, target, dimension (:,:), allocatable :: o2_xy3,o2_xy4,mach_xy3
  real, target, dimension (:,:), allocatable :: divu_xz,u2_xz,o2_xz,mach_xz
  real, target, dimension (:,:), allocatable :: divu_xz2,u2_xz2,o2_xz2,mach_xz2
  real, target, dimension (:,:), allocatable :: divu_yz,u2_yz,o2_yz,mach_yz
  real, target, dimension (:,:), allocatable :: ou_xy,ou_xy2,ou_xy3,ou_xy4
  real, target, dimension (:,:), allocatable :: ou_xz,ou_yz,ou_xz2
  real, target, dimension (:,:,:,:,:), allocatable :: divu_r,u2_r,o2_r,mach_r,ou_r
  real, target, dimension (:,:,:,:,:,:), allocatable :: oo_r,uu_sph_r

  real, dimension (mz,3) :: uumz=0.0, ruumz=0.0
  real, dimension (mx,3) :: uumx=0.0
  real, dimension (my,3) :: uumy=0.0
  real, dimension (:,:,:), allocatable :: uumxy, ruumxy
  real, dimension (mx,mz,3) :: uumxz=0.0
!
!  phi-averaged arrays for orbital advection
!
  real, dimension (mx,mz) :: uu_average_cyl=0.
  real, dimension (mx,my) :: uu_average_sph=0.
!
!  Cosine and sine function for setting test fields and analysis.
!
  real, dimension (mz) :: c2z,s2z,cz,sz
!
!  Profiles for setting differential rotation
!  (analogously to hydro_kinematic.f90).
!
  real, dimension(nx) :: profx_diffrot1=1., profx_diffrot2=1., profx_diffrot3=1.
  real, dimension(my) :: profy_diffrot1=1., profy_diffrot2=1., profy_diffrot3=1.
  real, dimension(mz) :: profz_diffrot1=1.
!
!  Precession matrices.
!
  real, dimension (3,3) :: mat_cori=0.,mat_cent=0.
!
!  Init parameters.
!
  real :: widthuu=.1, radiusuu=1., urand=0., kx_uu=1., ky_uu=1., kz_uu=1.
  real :: relhel_uu=1.,urandi=0.
  real :: uu_left=0.,uu_right=0.,uu_lower=1.,uu_upper=1.
  real :: uy_left=0.,uy_right=0.
  real :: initpower=1.,initpower2=-5./3.,cutoff=0.,ncutoff=1., kpeak=10.
  real :: xhalf, kgaussian_uu=0., nfact_uu=4.
  real, dimension (ninit) :: ampl_ux=0.0, ampl_uy=0.0, ampl_uz=0.0
  real, dimension (ninit) :: kx_ux=0.0, kx_uy=0.0, kx_uz=0.0
  real, dimension (ninit) :: ky_ux=0.0, ky_uy=0.0, ky_uz=0.0
  real, dimension (ninit) :: kz_ux=0.0, kz_uy=0.0, kz_uz=0.0
  real, dimension (ninit) :: phase_ux=0.0, phase_uy=0.0, phase_uz=0.0
  real :: omega_precession=0., alpha_precession=0.
  real, dimension (ninit) :: ampluu=0.0, uu_xz_angle=0.0
  character (len=labellen), dimension(ninit) :: inituu='nothing'
  character (len=labellen), dimension(3) :: borderuu='nothing'
  real, dimension (3) :: uu_const=(/0.,0.,0./), mean_momentum=(/0.,0.,0./)
  complex, dimension (3) :: coefuu=(/0.,0.,0./)
  real, dimension(nx) :: xmask_hyd
  real, dimension(nz) :: zmask_hyd
  real, dimension(nx) :: prof_om
  real, dimension(2) :: hydro_xaver_range=(/-max_real,max_real/)
  real, dimension(2) :: hydro_zaver_range=(/-max_real,max_real/)
  real :: u_out_kep=0.0, velocity_ceiling=.0, w_sldchar_hyd=1.0
  real :: mu_omega=0., gap=0., r_omega=0., w_omega=0.
  real :: z1_uu=0., z2_uu=0.
  real :: ABC_A=1., ABC_B=1., ABC_C=1.
  integer :: nb_rings=0
  integer :: neddy=0
!
! variables for expansion into spherical harmonics
!
  integer,parameter :: lSH_max=2
  integer, parameter :: Nmodes_SH=(lSH_max+1)*(lSH_max+1)
  integer :: index_rSH=1   !ceiling(nx/2.)
  real, dimension(nx) :: profile_SH=0.
!
  real, dimension (5) :: om_rings=0.
  integer :: N_modes_uu=0
  logical :: llinearized_hydro=.false.
  logical :: ladvection_velocity=.true.
  logical :: lprecession=.false.
  logical :: lshear_rateofstrain=.false.
  logical :: loo_as_aux = .false.
  logical :: luut_as_aux=.false., luust_as_aux=.false.
  logical :: loot_as_aux=.false., loost_as_aux=.false.
  logical :: luuk_as_aux=.false., look_as_aux=.false.
  logical :: luu_fluc_as_aux=.false.
  logical :: luu_sph_as_aux=.false.
  logical :: lscale_tobox=.true.
  logical :: lfactors_uu=.false.
  logical, target :: lpressuregradient_gas=.true.
  logical :: lcoriolis_force=.true.
  logical :: lshear_in_coriolis=.false.
  logical :: lcentrifugal_force=.false.
  logical, pointer :: lffree
  logical :: lreflecteddy=.false.,louinit=.false.
  logical :: lskip_projection=.false.
  logical :: lrelativistic=.false.
  real, pointer :: profx_ffree(:),profy_ffree(:),profz_ffree(:)
  real :: incl_alpha = 0.0, rot_rr = 0.0
  real :: xsphere = 0.0, ysphere = 0.0, zsphere = 0.0
  real :: amp_meri_circ = 0.0
  real :: max_uu = 0., delta_u=1
! The following is useful to debug the forcing - Dhruba
  real :: outest
  real :: ampl_Omega=0.0
  real :: omega_ini=0.0
  logical :: loutest,ldiffrot_test=.false.
  real :: r_cyl = 1.0, skin_depth = 1e-1
  real :: rnoise_int=impossible,rnoise_ext=impossible
  real :: PrRa  !preliminary
  real :: amp_factor=0.,kx_uu_perturb=0.
  integer, dimension(ninit) :: ll_sh=0, mm_sh=0, n_xprof=-1
!
  namelist /hydro_init_pars/ &
      ampluu, ampl_ux, ampl_uy, ampl_uz, phase_ux, phase_uy, phase_uz, &
      inituu, widthuu, radiusuu, urand, urandi, lpressuregradient_gas, &
      uu_xz_angle, relhel_uu, coefuu, r_omega, w_omega,&
      uu_left, uu_right, uu_lower, uu_upper, kx_uu, ky_uu, kz_uu, &
      kx_ux, ky_ux, kz_ux, kx_uy, ky_uy, kz_uy, kx_uz, ky_uz, kz_uz, &
      uy_left, uy_right, uu_const, Omega, u_out_kep, &
      initpower, initpower2, cutoff, ncutoff, kpeak, kgaussian_uu, &
      lrelativistic, lskip_projection, z1_uu, z2_uu, &
      N_modes_uu, lcoriolis_force, lcentrifugal_force, ladvection_velocity, &
      lprecession, omega_precession, alpha_precession, velocity_ceiling, &
      loo_as_aux, luut_as_aux, luust_as_aux, loot_as_aux, loost_as_aux, &
      luuk_as_aux, look_as_aux, &
      mu_omega, nb_rings, om_rings, gap, &
      lscale_tobox, ampl_Omega, omega_ini, r_cyl, skin_depth, incl_alpha, &
      rot_rr, xsphere, ysphere, zsphere, neddy, amp_meri_circ, &
      rnoise_int, rnoise_ext, lreflecteddy, louinit, hydro_xaver_range, max_uu,&
      amp_factor,kx_uu_perturb,llinearized_hydro, hydro_zaver_range, index_rSH, &
      ll_sh, mm_sh, delta_u, n_xprof, luu_fluc_as_aux, luu_sph_as_aux, nfact_uu, &
      lfactors_uu
!
!  Run parameters.
!
  real :: tdamp=0.,tfade_start=-1.,dampu=0.,wdamp=0.
  real :: dampuint=0.0,dampuext=0.0,rdampint=impossible,rdampext=impossible
  real :: ydampint=impossible,ydampext=impossible
  real :: ruxm=0.,ruym=0.,ruzm=0.
  real :: tau_damp_ruxm1=0.,tau_damp_ruym1=0.,tau_damp_ruzm1=0.
  real :: tau_damp_ruxm=0.,tau_damp_ruym=0.,tau_damp_ruzm=0.,tau_diffrot1=0.
  real :: ampl1_diffrot=0.,ampl2_diffrot=0., ampl_wind=0.
  real :: Omega_int=0.,xexp_diffrot=1.,kx_diffrot=1.,kz_diffrot=0., phase_diffrot=0.
  real :: othresh=0.,othresh_per_orms=0.,orms=0.,othresh_scl=1.
  real :: omega_out=0., omega_in=0., omega_fourier=0.
  real :: width_ff_uu=1.,x1_ff_uu=0.,x2_ff_uu=0.
  real :: ekman_friction=0.0, uzjet=0.0
  real :: ampl_forc=0., k_forc=impossible, w_forc=0., x_forc=0., dx_forc=0.1
  real :: ampl_fcont_uu=1., k_diffrot=1., amp_centforce=1.
  real :: uphi_rbot=1., uphi_rtop=1., uphi_step_width=0.
  integer :: novec,novecmax=nx*ny*nz/4
  logical :: ldamp_fade=.false.,lOmega_int=.false.,lupw_uu=.false.
  logical :: lfreeze_uint=.false.,lfreeze_uext=.false.
  logical :: lremove_mean_angmom=.false.
  logical :: lremove_mean_momenta=.false.
  logical :: lremove_mean_flow=.false.
  logical :: lremove_uumeanxy=.false.
  logical :: lremove_uumeanx=.false., lremove_uumeany=.false., lremove_uumeanz=.false.
  logical :: lremove_uumeanz_horizontal=.false.
  logical :: lreinitialize_uu=.false.
  logical :: lalways_use_gij_etc=.false.
  logical :: lcalc_uumeanz=.false.,lcalc_uumeanxy=.false.,lcalc_uumean
  logical :: lcalc_uumeanx=.false.,lcalc_uumeany=.false.,lcalc_uumeanxz=.false.
  logical :: lcalc_ruumeanz=.false.,lcalc_ruumeanxy=.false.
  logical :: lforcing_cont_uu=.false.
  logical :: lcoriolis_xdep=.false.
  logical :: lno_meridional_flow=.false.
  logical :: lrotation_xaxis=.false.
  logical :: lpropagate_borderuu=.true.
  logical :: lgradu_as_aux=.false.
  logical :: lOmega_cyl_xy=.false.
  logical :: lno_radial_advection=.false.
  logical :: lfargoadvection_as_shift=.true.
  logical :: lhelmholtz_decomp=.false.
  logical :: limpose_only_horizontal_uumz=.false.
  logical :: ltime_integrals_always=.true.
  logical :: lvart_in_shear_frame=.false.
  real :: dtcor=0.
  character (len=labellen) :: uuprof='nothing'
!
!  Parameters for interior boundary conditions.
!
  character (len=labellen) :: interior_bc_hydro_profile='nothing'
  logical :: lhydro_bc_interior=.false.
  real :: z1_interior_bc_hydro=0.,kz_analysis=1.
  real :: Shearx=0., rescale_uu=0.
  real :: Ra=0.0, Pr=0.0 ! Boussinesq approximation
  real :: Om_inner=0.
!
!  Option to constrain time for large df.
!
  real :: cdt_tauf=1.0, ulev=1.0
  logical :: lcdt_tauf=.false.
  logical, target :: lcalc_uuavg=.false.
!
  namelist /hydro_run_pars/ &
      Omega, theta, tdamp, dampu, dampuext, dampuint, rdampext, rdampint, &
      wdamp, tau_damp_ruxm, tau_damp_ruym, tau_damp_ruzm, tau_diffrot1, &
      inituu, ampluu, kz_uu, ampl1_diffrot, ampl2_diffrot, uuprof, &
      xexp_diffrot, kx_diffrot, kz_diffrot, kz_analysis, phase_diffrot, ampl_wind, &
      lreinitialize_uu, lremove_mean_momenta, lremove_mean_flow, lremove_uumeanx,&
      lremove_uumeany,lremove_uumeanxy,lremove_uumeanz,lremove_uumeanz_horizontal, &
      ldamp_fade, tfade_start, lOmega_int, Omega_int, lupw_uu, othresh, &
      othresh_per_orms, borderuu, lfreeze_uint, lpressuregradient_gas, &
      lfreeze_uext, lcoriolis_force, lcentrifugal_force, ladvection_velocity, &
      omega_out, omega_in, lprecession, omega_precession, omega_fourier, &
      alpha_precession, lshear_rateofstrain, r_omega, w_omega, &
      lrelativistic, lalways_use_gij_etc, amp_centforce, &
      lcalc_uumean,lcalc_uumeanx,lcalc_uumeanxy,lcalc_uumeanxz,lcalc_uumeanz, &
      lcalc_ruumeanz, lcalc_ruumeanxy, &
      lforcing_cont_uu, width_ff_uu, x1_ff_uu, x2_ff_uu, &
      loo_as_aux, luut_as_aux, luust_as_aux, loot_as_aux, loost_as_aux, &
      loutest, ldiffrot_test, &
      interior_bc_hydro_profile, lhydro_bc_interior, z1_interior_bc_hydro, &
      velocity_ceiling, ekman_friction, ampl_Omega, lcoriolis_xdep, &
      ampl_forc, k_forc, w_forc, x_forc, dx_forc, ampl_fcont_uu, &
      lno_meridional_flow, lrotation_xaxis, k_diffrot,Shearx, rescale_uu, &
      hydro_xaver_range, Ra, Pr, llinearized_hydro, lremove_mean_angmom, &
      lpropagate_borderuu, hydro_zaver_range, index_rSH, &
      uzjet, ydampint, ydampext, mean_momentum, lshear_in_coriolis, &
      lcdt_tauf, cdt_tauf, ulev, &
      w_sldchar_hyd, uphi_rbot, uphi_rtop, uphi_step_width, lOmega_cyl_xy, &
      lno_radial_advection, lfargoadvection_as_shift, lhelmholtz_decomp, &
      limpose_only_horizontal_uumz, luu_fluc_as_aux, Om_inner, luu_sph_as_aux, &
      ltime_integrals_always, dtcor, lvart_in_shear_frame
!
!  Diagnostic variables (need to be consistent with reset list below).
!
  integer :: idiag_u2tm=0       ! DIAG_DOC: $\left<\uv(t)\cdot\int_0^t\uv(t')
                                ! DIAG_DOC:   dt'\right>$
  integer :: idiag_uotm=0       ! DIAG_DOC: $\left<\uv(t)\cdot\int_0^t\omv(t')
                                ! DIAG_DOC:   dt'\right>$
  integer :: idiag_outm=0       ! DIAG_DOC: $\left<\omv(t)\cdot\int_0^t\uv(t')
                                ! DIAG_DOC:   dt'\right>$
  integer :: idiag_fkinzm=0     ! DIAG_DOC: $\left<{1\over2} \varrho\uv^2 u_z\right>$
  integer :: idiag_u2m=0        ! DIAG_DOC: $\left<\uv^2\right>$
  integer :: idiag_u2sphm=0     ! DIAG_DOC: $\int_{r=0}^{r=1} \uv^2 dV$,
                                ! DIAG_DOC:   where $r=\sqrt{x^2+y^2+z^2}$
  integer :: idiag_um2=0        ! DIAG_DOC:
  integer :: idiag_uxpt=0       ! DIAG_DOC: $u_x(x_1,y_1,z_1,t)$
  integer :: idiag_uypt=0       ! DIAG_DOC: $u_y(x_1,y_1,z_1,t)$
  integer :: idiag_uzpt=0       ! DIAG_DOC: $u_z(x_1,y_1,z_1,t)$
  integer :: idiag_uxp2=0       ! DIAG_DOC: $u_x(x_2,y_2,z_2,t)$
  integer :: idiag_uyp2=0       ! DIAG_DOC: $u_y(x_2,y_2,z_2,t)$
  integer :: idiag_uzp2=0       ! DIAG_DOC: $u_z(x_2,y_2,z_2,t)$
  integer :: idiag_uxuypt=0     ! DIAG_DOC: $(u_x u_y)(x_1,y_1,z_1,t)$
  integer :: idiag_uyuzpt=0     ! DIAG_DOC: $(u_y u_z)(x_1,y_1,z_1,t)$
  integer :: idiag_uzuxpt=0     ! DIAG_DOC: $(u_z u_x)(x_1,y_1,z_1,t)$
  integer :: idiag_urms=0       ! DIAG_DOC: $\left<\uv^2\right>^{1/2}$
  integer :: idiag_urmsx=0      ! DIAG_DOC: $\left<\uv^2\right>^{1/2}$ for
                                ! DIAG_DOC: the hydro_xaver_range
  integer :: idiag_urmsz=0      ! DIAG_DOC: $\left<\uv^2\right>^{1/2}$ for
                                ! DIAG_DOC: the hydro_zaver_range
  integer :: idiag_durms=0      ! DIAG_DOC: $\left<\delta\uv^2\right>^{1/2}$
  integer :: idiag_umax=0       ! DIAG_DOC: $\max(|\uv|)$
  integer :: idiag_umin=0       ! DIAG_DOC: $\min(|\uv|)$
  integer :: idiag_uxrms=0      ! DIAG_DOC: $\left<u_x^2\right>^{1/2}$
  integer :: idiag_uyrms=0      ! DIAG_DOC: $\left<u_y^2\right>^{1/2}$
  integer :: idiag_uzrms=0      ! DIAG_DOC: $\left<u_z^2\right>^{1/2}$
  integer :: idiag_uzrmaxs=0    ! DIAG_DOC:
  integer :: idiag_uxmin=0      ! DIAG_DOC: $\min(|u_x|)$
  integer :: idiag_uymin=0      ! DIAG_DOC: $\min(|u_y|)$
  integer :: idiag_uzmin=0      ! DIAG_DOC: $\min(|u_z|)$
  integer :: idiag_uxmax=0      ! DIAG_DOC: $\max(|u_x|)$
  integer :: idiag_uymax=0      ! DIAG_DOC: $\max(|u_y|)$
  integer :: idiag_uzmax=0      ! DIAG_DOC: $\max(|u_z|)$
  integer :: idiag_uxm=0        ! DIAG_DOC: $\left<u_x\right>$
  integer :: idiag_uym=0        ! DIAG_DOC: $\left<u_y\right>$
  integer :: idiag_uzm=0        ! DIAG_DOC: $\left<u_z\right>$
  integer :: idiag_uduum=0      ! DIAG_DOC: $\left<\boldsymbol{u}\cdot\boldsymbol{u}\right>$
  integer :: idiag_ux2m=0       ! DIAG_DOC: $\left<u_x^2\right>$
  integer :: idiag_uy2m=0       ! DIAG_DOC: $\left<u_y^2\right>$
  integer :: idiag_uz2m=0       ! DIAG_DOC: $\left<u_z^2\right>$
  integer :: idiag_ux3m=0       ! DIAG_DOC: $\left<u_x^3\right>$
  integer :: idiag_uy3m=0       ! DIAG_DOC: $\left<u_y^3\right>$
  integer :: idiag_uz3m=0       ! DIAG_DOC: $\left<u_z^3\right>$
  integer :: idiag_ux4m=0       ! DIAG_DOC: $\left<u_x^4\right>$
  integer :: idiag_uy4m=0       ! DIAG_DOC: $\left<u_y^4\right>$
  integer :: idiag_uz4m=0       ! DIAG_DOC: $\left<u_z^4\right>$
  integer :: idiag_uxuy2m=0     ! DIAG_DOC: $\left<u_x^2u_y^2\right>$
  integer :: idiag_uyuz2m=0     ! DIAG_DOC: $\left<u_y^2u_z^2\right>$
  integer :: idiag_uzux2m=0     ! DIAG_DOC: $\left<u_z^2u_x^2\right>$
  integer :: idiag_ux2ccm=0     ! DIAG_DOC: $\left<u_x^2\cos^2kz\right>$
  integer :: idiag_ux2ssm=0     ! DIAG_DOC: $\left<u_x^2\sin^2kz\right>$
  integer :: idiag_uy2ccm=0     ! DIAG_DOC: $\left<u_y^2\cos^2kz\right>$
  integer :: idiag_uy2ssm=0     ! DIAG_DOC: $\left<u_y^2\sin^2kz\right>$
  integer :: idiag_uxuycsm=0    ! DIAG_DOC: $\left<u_xu_y\cos kz\sin kz\right>$
  integer :: idiag_uxuym=0      ! DIAG_DOC: $\left<u_x u_y\right>$
  integer :: idiag_uxuzm=0      ! DIAG_DOC: $\left<u_x u_z\right>$
  integer :: idiag_uyuzm=0      ! DIAG_DOC: $\left<u_y u_z\right>$
  integer :: idiag_umx=0        ! DIAG_DOC: $\left< u_x \right>$
  integer :: idiag_umy=0        ! DIAG_DOC: $\left< u_y \right>$
  integer :: idiag_umz=0        ! DIAG_DOC: $\left< u_z \right>$
  integer :: idiag_omumz=0      ! DIAG_DOC: $\left<\left<\Wv\right>_{xy}
                                ! DIAG_DOC:   \cdot\left<\Uv\right>_{xy}
                                ! DIAG_DOC:   \right>$ \quad($xy$-averaged
                                ! DIAG_DOC:   mean cross helicity production)
  integer :: idiag_umamz=0      ! DIAG_DOC: $\left<\left<\uv\right>_{xy}\cdot\left<\Av\right>_{xy}\right>$
  integer :: idiag_umbmz=0      ! DIAG_DOC: $\left<\left<\Uv\right>_{xy}
                                ! DIAG_DOC:   \cdot\left<\Bv\right>_{xy}
                                ! DIAG_DOC:   \right>$ \quad($xy$-averaged
                                ! DIAG_DOC:   mean cross helicity production)
  integer :: idiag_umxbmz=0     ! DIAG_DOC: $\left<\left<\Uv\right>_{xy}
                                ! DIAG_DOC:   \times\left<\Bv\right>_{xy}
                                ! DIAG_DOC:   \right>_z$ \quad($xy$-averaged
                                ! DIAG_DOC:   mean emf)
  integer :: idiag_rux2m=0      ! DIAG_DOC: $\left<\rho u_x^2\right>$
  integer :: idiag_ruy2m=0      ! DIAG_DOC: $\left<\rho u_y^2\right>$
  integer :: idiag_ruz2m=0      ! DIAG_DOC: $\left<\rho u_z^2\right>$
  integer :: idiag_divum=0      ! DIAG_DOC: $\left<{\rm div}\uv)\right>$
  integer :: idiag_rdivum=0     ! DIAG_DOC: $\left<\varrho{\rm div}\uv)\right>$
  integer :: idiag_divu2m=0     ! DIAG_DOC: $\left<({\rm div}\uv)^2\right>$
  integer :: idiag_gdivu2m=0    ! DIAG_DOC: $\left<({\rm grad\,div}\uv)^2\right>$
  integer :: idiag_u3u21m=0     ! DIAG_DOC: $\left<u_3 u_{2,1}\right>$
  integer :: idiag_u1u32m=0     ! DIAG_DOC: $\left<u_1 u_{3,2}\right>$
  integer :: idiag_u2u13m=0     ! DIAG_DOC: $\left<u_2 u_{1,3}\right>$
  integer :: idiag_u2u31m=0     ! DIAG_DOC: $\left<u_2 u_{3,1}\right>$
  integer :: idiag_u3u12m=0     ! DIAG_DOC: $\left<u_3 u_{1,2}\right>$
  integer :: idiag_u1u23m=0     ! DIAG_DOC: $\left<u_1 u_{2,3}\right>$
  integer :: idiag_urmphi=0     ! PHIAVG_DOC: $\left<u_\varpi\right>_\varphi$
                                ! PHIAVG_DOC: [cyl.\ polar coords
                                ! PHIAVG_DOC:  $(\varpi,\varphi,z)$]
  integer :: idiag_upmphi=0     ! PHIAVG_DOC: $\left<u_\varphi\right>_\varphi$
  integer :: idiag_uzmphi=0     ! PHIAVG_DOC: $\left<u_z\right>_\varphi$
  integer :: idiag_rurmphi=0    ! PHIAVG_DOC: $\left<\rho u_\varpi\right>_\varphi$
  integer :: idiag_rupmphi=0    ! PHIAVG_DOC: $\left<\rho u_\varphi\right>_\varphi$
  integer :: idiag_ruzmphi=0    ! PHIAVG_DOC: $\left<\rho u_z\right>_\varphi$
  integer :: idiag_ur2mphi=0    ! PHIAVG_DOC: $\left<u_\varpi^2\right>_\varphi$
  integer :: idiag_up2mphi=0    ! PHIAVG_DOC: $\left<u_\varphi^2\right>_\varphi$
  integer :: idiag_uz2mphi=0    ! PHIAVG_DOC: $\left<u_z^2\right>_\varphi$
  integer :: idiag_urupmphi=0   ! PHIAVG_DOC: $\left<u_\varpi u_\varphi\right>_\varphi$
  integer :: idiag_uruzmphi=0   ! PHIAVG_DOC: $\left<u_\varpi u_z \right>_\varphi$
  integer :: idiag_upuzmphi=0   ! PHIAVG_DOC: $\left<u_\varphi u_z \right>_\varphi$
  integer :: idiag_rurupmphi=0  ! PHIAVG_DOC: $\left<\rho u_\varpi u_\varphi\right>_\varphi$
  integer :: idiag_ruruzmphi=0  ! PHIAVG_DOC: $\left<\rho u_\varpi u_z \right>_\varphi$
  integer :: idiag_rupuzmphi=0  ! PHIAVG_DOC: $\left<\rho u_\varphi u_z \right>_\varphi$
  integer :: idiag_ursphmphi=0  ! PHIAVG_DOC: $\left<u_r\right>_\varphi$
  integer :: idiag_uthmphi=0    ! PHIAVG_DOC: $\left<u_\vartheta\right>_\varphi$
  integer :: idiag_rursphmphi=0 ! PHIAVG_DOC: $\left<\rho u_r\right>_\varphi$
  integer :: idiag_ruthmphi=0   ! PHIAVG_DOC: $\left<\rho u_\vartheta\right>_\varphi$
  ! For the manual: uumphi      ! PHIAVG_DOC: shorthand for \var{urmphi},
                                ! PHIAVG_DOC: \var{upmphi} and \var{uzmphi}
                                ! PHIAVG_DOC: together
  ! For the manual: uusphmphi   ! PHIAVG_DOC: shorthand for \var{ursphmphi},
                                ! PHIAVG_DOC: \var{uthmphi} and \var{upmphi}
                                ! PHIAVG_DOC: together
  integer :: idiag_u2mphi=0     ! PHIAVG_DOC: $\left<\uv^2\right>_\varphi$
  integer :: idiag_fkinrsphmphi=0 ! PHIAVG_DOC: $\left<{1\over2}\varrho\uv^2
                                ! PHIAVG_DOC: u_r\right>_{\varphi}$
  integer :: idiag_u2mr=0       ! DIAG_DOC:
  integer :: idiag_urmr=0       ! DIAG_DOC:
  integer :: idiag_upmr=0       ! DIAG_DOC:
  integer :: idiag_uzmr=0       ! DIAG_DOC:
  integer :: idiag_uxfampm=0    ! DIAG_DOC:
  integer :: idiag_uyfampm=0    ! DIAG_DOC:
  integer :: idiag_uzfampm=0    ! DIAG_DOC:
  integer :: idiag_uxfampim=0   ! DIAG_DOC:
  integer :: idiag_uyfampim=0   ! DIAG_DOC:
  integer :: idiag_uzfampim=0   ! DIAG_DOC:
  integer :: idiag_ruxm=0       ! DIAG_DOC: $\left<\varrho u_x\right>$
                                ! DIAG_DOC:   \quad(mean $x$-momentum density)
  integer :: idiag_ruym=0       ! DIAG_DOC: $\left<\varrho u_y\right>$
                                ! DIAG_DOC:   \quad(mean $y$-momentum density)
  integer :: idiag_ruzm=0       ! DIAG_DOC: $\left<\varrho u_z\right>$
                                ! DIAG_DOC:   \quad(mean $z$-momentum density)
  integer :: idiag_ruxtot=0     ! DIAG_DOC: $\left<\rho |u|\right>$
                                ! DIAG_DOC:   \quad(mean absolute $x$-momentum density)
  integer :: idiag_rumax=0      ! DIAG_DOC: $\max(\varrho |\uv|)$
                                ! DIAG_DOC:   \quad(maximum modulus of momentum)
  integer :: idiag_ruxuym=0     ! DIAG_DOC: $\left<\varrho u_x u_y\right>$
                                ! DIAG_DOC:   \quad(mean Reynolds stress)
  integer :: idiag_ruxuzm=0     ! DIAG_DOC: $\left<\varrho u_x u_z\right>$
                                ! DIAG_DOC:   \quad(mean Reynolds stress)
  integer :: idiag_ruyuzm=0     ! DIAG_DOC: $\left<\varrho u_y u_z\right>$
                                ! DIAG_DOC:   \quad(mean Reynolds stress)
  integer :: idiag_divrhourms=0 ! DIAG_DOC: $\left|\nabla\cdot(\varrho\uv)\right|_{\rm rms}$
  integer :: idiag_divrhoumax=0 ! DIAG_DOC: $\left|\nabla\cdot(\varrho\uv)\right|_{\rm max}$
  integer :: idiag_rlxm=0       ! DIAG_DOC: $\left< \rho y u_z - z u_y \right>$
  integer :: idiag_rlym=0       ! DIAG_DOC: $\left< \rho z u_x - x u_z \right>$
  integer :: idiag_rlzm=0       ! DIAG_DOC: $\left< \rho x u_y - y u_x \right>$
  integer :: idiag_rlx2m=0      ! DIAG_DOC: $\left<(\rho y u_z-z u_y)^2\right>$
  integer :: idiag_rly2m=0      ! DIAG_DOC: $\left<(\rho z u_x-x u_z)^2\right>$
  integer :: idiag_rlz2m=0      ! DIAG_DOC: $\left<(\rho x u_y-y u_x)^2\right>$
  integer :: idiag_tot_ang_mom=0! DIAG_DOC: Total angular momentum in spherical
                                ! DIAG_DOC: coordinates about the axis.
  integer :: idiag_dtu=0        ! DIAG_DOC: $\delta t/[c_{\delta t}\,\delta x
                                ! DIAG_DOC:  /\max|\mathbf{u}|]$
                                ! DIAG_DOC:  \quad(time step relative to
                                ! DIAG_DOC:   advective time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_oum=0        ! DIAG_DOC: $\left<\boldsymbol{\omega}\cdot\uv\right>$
  integer :: idiag_oxum=0       ! DIAG_DOC: $\left<\boldsymbol{\omega}\times\uv\right>$
  integer :: idiag_ourms=0      ! DIAG_DOC: $\left<(\boldsymbol{\omega}\cdot\uv)^2\right>^{1/2}$
  integer :: idiag_oxurms=0     ! DIAG_DOC: $\left<(\boldsymbol{\omega}\times\uv)^2\right>^{1/2}$
  integer :: idiag_ou_int=0     ! DIAG_DOC: $\int_V\boldsymbol{\omega}\cdot\uv\,dV$
  integer :: idiag_fum=0        ! DIAG_DOC: $\left<\fv\cdot\uv\right>$
  integer :: idiag_odel2um=0    ! DIAG_DOC: $\left<\boldsymbol{\omega}\nabla^2\uv\right>$
  integer :: idiag_o2m=0        ! DIAG_DOC: $\left<\boldsymbol{\omega}^2\right>
                                ! DIAG_DOC:   \equiv \left<(\curl\uv)^2\right>$
  integer :: idiag_o2u2m=0      ! DIAG_DOC: $\left<\uv^2\boldsymbol{\omega}^2\right>$
  integer :: idiag_o2sphm=0     ! DIAG_DOC: $\int_{r=0}^{r=1} \boldsymbol{\omega}^2 dV$,
                                ! DIAG_DOC:   where $r=\sqrt{x^2+y^2+z^2}$
  integer :: idiag_orms=0       ! DIAG_DOC: $\left<\boldsymbol{\omega}^2\right>^{1/2}$
  integer :: idiag_omax=0       ! DIAG_DOC: $\max(|\boldsymbol{\omega}|)$
  integer :: idiag_ox2m=0       ! DIAG_DOC: $\left<\omega_x^2\right>$
  integer :: idiag_oy2m=0       ! DIAG_DOC: $\left<\omega_y^2\right>$
  integer :: idiag_oz2m=0       ! DIAG_DOC: $\left<\omega_z^2\right>$
  integer :: idiag_ox3m=0       ! DIAG_DOC: $\left<\omega_x^3\right>$
  integer :: idiag_oy3m=0       ! DIAG_DOC: $\left<\omega_y^3\right>$
  integer :: idiag_oz3m=0       ! DIAG_DOC: $\left<\omega_z^3\right>$
  integer :: idiag_ox4m=0       ! DIAG_DOC: $\left<\omega_x^4\right>$
  integer :: idiag_oy4m=0       ! DIAG_DOC: $\left<\omega_y^4\right>$
  integer :: idiag_oz4m=0       ! DIAG_DOC: $\left<\omega_z^4\right>$
  integer :: idiag_oxm=0        ! DIAG_DOC:
  integer :: idiag_oym=0        ! DIAG_DOC:
  integer :: idiag_ozm=0        ! DIAG_DOC:
  integer :: idiag_oxuzxm=0     ! DIAG_DOC: $\left<\omega_x u_{z,x} \right>$
  integer :: idiag_oyuzym=0     ! DIAG_DOC: $\left<\omega_y u_{z,y} \right>$
  integer :: idiag_oxoym=0      ! DIAG_DOC: $\left<\omega_x\omega_y\right>$
  integer :: idiag_oxozm=0      ! DIAG_DOC: $\left<\omega_x\omega_z\right>$
  integer :: idiag_oyozm=0      ! DIAG_DOC: $\left<\omega_y\omega_z\right>$
  integer :: idiag_qfm=0        ! DIAG_DOC: $\left<\qv\cdot\fv\right>$
  integer :: idiag_q2m=0        ! DIAG_DOC: $\left<\qv^2\right>$
  integer :: idiag_qrms=0       ! DIAG_DOC: $\left<\qv^2\right>^{1/2}$
  integer :: idiag_qmax=0       ! DIAG_DOC: $\max(|\qv|)$
  integer :: idiag_qom=0        ! DIAG_DOC: $\left<\qv\cdot\omv\right>$
  integer :: idiag_quxom=0      ! DIAG_DOC: $\left<\qv\cdot(\uv\times\omv)\right>$
  integer :: idiag_pvzm=0       ! DIAG_DOC: $\left<\omega_z + 2\Omega/\varrho\right>$
                                ! DIAG_DOC: \quad(z component of potential vorticity)
  integer :: idiag_oumphi=0     ! DIAG_DOC: $\left<\omv\cdot\uv\right>_\varphi$
  integer :: idiag_ozmphi=0     ! DIAG_DOC:
  integer :: idiag_ormr=0       ! DIAG_DOC:
  integer :: idiag_opmr=0       ! DIAG_DOC:
  integer :: idiag_ozmr=0       ! DIAG_DOC:
  integer :: idiag_uguxm=0      ! DIAG_DOC:
  integer :: idiag_uguym=0      ! DIAG_DOC:
  integer :: idiag_uguzm=0      ! DIAG_DOC:
  integer :: idiag_ugurmsx=0    ! DIAG_DOC: $\left<\left(\uv\nabla\uv\right)^2\right>^{1/2}$
                                ! DIAG_DOC: for the hydro_xaver_range
  integer :: idiag_ugu2m=0      ! DIAG_DOC: $\left<\uv\nabla\uv\right>^2$
  integer :: idiag_dudx=0        ! DIAG_DOC: $\left<\frac{\delta \uv}{\delta x}\right>$
  integer :: idiag_Marms=0      ! DIAG_DOC: $\left<\uv^2/\cs^2\right>$
                                ! DIAG_DOC:   \quad(rms Mach number)
  integer :: idiag_Mamax=0      ! DIAG_DOC: $\max |\uv|/\cs$
                                ! DIAG_DOC:   \quad(maximum Mach number)
  integer :: idiag_fintm=0      ! DIAG_DOC:
  integer :: idiag_fextm=0      ! DIAG_DOC:
  integer :: idiag_duxdzma=0    ! DIAG_DOC:
  integer :: idiag_duydzma=0    ! DIAG_DOC:
  integer :: idiag_EEK=0        ! DIAG_DOC: $\left<\varrho\uv^2\right>/2$
  integer :: idiag_ekin=0       ! DIAG_DOC: $\left<{1\over2}\varrho\uv^2\right>$
  integer :: idiag_ekintot=0    ! DIAG_DOC: $\int_V{1\over2}\varrho\uv^2\, dV$
  integer :: idiag_totangmom=0  ! DIAG_DOC:
  integer :: idiag_uxglnrym=0   ! DIAG_DOC: $\left<u_x\partial_y\ln\varrho\right>$
  integer :: idiag_uyglnrxm=0   ! DIAG_DOC: $\left<u_y\partial_x\ln\varrho\right>$
  integer :: idiag_uzdivum=0    ! DIAG_DOC: $\left<u_z\nabla\cdot\uv\right>$
  integer :: idiag_uxuydivum=0  ! DIAG_DOC: $\left<u_x u_y\nabla\cdot\uv\right>$
  integer :: idiag_divuHrms=0   ! DIAG_DOC: $(\nabla_{\rm H}\cdot\uv_{\rm H})^{\rm rms}$
  integer :: idiag_uxxrms=0     ! DIAG_DOC: $u_{x,x}^{\rm rms}$
  integer :: idiag_uyyrms=0     ! DIAG_DOC: $u_{y,y}^{\rm rms}$
  integer :: idiag_uxzrms=0     ! DIAG_DOC: $u_{x,z}^{\rm rms}$
  integer :: idiag_uyzrms=0     ! DIAG_DOC: $u_{y,z}^{\rm rms}$
  integer :: idiag_uzyrms=0     ! DIAG_DOC: $u_{z,y}^{\rm rms}$
  integer :: idiag_urmsn=0,idiag_urmss=0,idiag_urmsh=0
  integer :: idiag_ormsn=0,idiag_ormss=0,idiag_ormsh=0
  integer :: idiag_oumn=0,idiag_oums=0,idiag_oumh=0
  integer :: idiag_taufmin=0
  integer :: idiag_dtF=0        ! DIAG_DOC: $\delta t/[c_{\delta t}\,\delta x
                                ! DIAG_DOC:  /\max|\mathbf{F}|]$
                                ! DIAG_DOC:  \quad(time step relative to
                                ! DIAG_DOC:   max force time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
!
  integer, dimension(Nmodes_SH) :: idiag_urlm=0 ! DIAG_DOC: $ \int u_r(\theta,\phi)Y^m_{\ell}(\theta,\phi)\sin(\theta)d\theta d\phi$
  integer :: idiag_udpxxm=0, &  ! DIAG_DOC: components of symmetric tensor
             idiag_udpyym=0, &  ! DIAG_DOC: $\left< u_i \partial_j p + u_j \partial_i p \right>$
             idiag_udpzzm=0, &
             idiag_udpxym=0, &
             idiag_udpyzm=0, &
             idiag_udpxzm=0
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_u2mz=0       ! XYAVG_DOC: $\left< \uv^2 \right>_{xy}$
  integer :: idiag_o2mz=0       ! XYAVG_DOC: $\left< \Wv^2 \right>_{xy}$
  integer :: idiag_divu2mz=0    ! XYAVG_DOC: $\left<(\nabla\cdot\uv)^2\right>_{xy}$
  integer :: idiag_curlru2mz=0  ! XYAVG_DOC: $\left<(\nabla\times\varrho\Uv)^2 \right>_{xy}$
  integer :: idiag_divru2mz=0   ! XYAVG_DOC: $\left<(\nabla\cdot\varrho\uv)^2\right>_{xy}$
  integer :: idiag_fmasszmz=0   ! XYAVG_DOC: $\left< \varrho u_z \right>_{xy}$
  integer :: idiag_fkinzmz=0    ! XYAVG_DOC: $\left<{1\over2}\varrho\uv^2 u_z\right>_{xy}$
  integer :: idiag_fkinzupmz=0  ! XYAVG_DOC: $\left<{1\over2}\varrho\uv^2 u_{z\uparrow}\right>_{xy}$
  integer :: idiag_fkinzdownmz=0 ! XYAVG_DOC: $\left<{1\over2}\varrho\uv^2 u_{z\downarrow}\right>_{xy}$
  integer :: idiag_uxmz=0       ! XYAVG_DOC: $\left< u_x \right>_{xy}$
                                ! XYAVG_DOC:   \quad(horiz. averaged $x$
                                ! XYAVG_DOC:   velocity)
  integer :: idiag_uymz=0       ! XYAVG_DOC: $\left< u_y \right>_{xy}$
  integer :: idiag_uzmz=0       ! XYAVG_DOC: $\left< u_z \right>_{xy}$
  integer :: idiag_ffdownmz=0   ! XYAVG_DOC: Filling factor of downflows
  integer :: idiag_uzupmz=0     ! XYAVG_DOC: $\left< u_{z\uparrow} \right>_{xy}$
  integer :: idiag_uzdownmz=0   ! XYAVG_DOC: $\left< u_{z\downarrow} \right>_{xy}$
  integer :: idiag_ruzupmz=0    ! XYAVG_DOC: $\left< \varrho u_{z\uparrow} \right>_{xy}$
  integer :: idiag_ruzdownmz=0  ! XYAVG_DOC: $\left< \varrho u_{z\downarrow} \right>_{xy}$
  integer :: idiag_divumz=0     ! XYAVG_DOC: $\left< {\rm div}\uv \right>_{xy}$
  integer :: idiag_uzdivumz=0   ! XYAVG_DOC: $\left< u_z{\rm div}\uv \right>_{xy}$
  integer :: idiag_oxmz=0       ! XYAVG_DOC: $\left< \omega_x \right>_{xy}$
  integer :: idiag_oymz=0       ! XYAVG_DOC: $\left< \omega_y \right>_{xy}$
  integer :: idiag_ozmz=0       ! XYAVG_DOC: $\left< \omega_z \right>_{xy}$
  integer :: idiag_ux2mz=0      ! XYAVG_DOC: $\left<u_x^2\right>_{xy}$
  integer :: idiag_uy2mz=0      ! XYAVG_DOC: $\left<u_y^2\right>_{xy}$
  integer :: idiag_uz2mz=0      ! XYAVG_DOC: $\left<u_z^2\right>_{xy}$
  integer :: idiag_uz2upmz=0    ! XYAVG_DOC: $\left<u_{z\uparrow}^2\right>_{xy}$
  integer :: idiag_uz2downmz=0  ! XYAVG_DOC: $\left<u_{z\downarrow}^2\right>_{xy}$
  integer :: idiag_ox2mz=0      ! XYAVG_DOC: $\left< \omega_x^2 \right>_{xy}$
  integer :: idiag_oy2mz=0      ! XYAVG_DOC: $\left< \omega_y^2 \right>_{xy}$
  integer :: idiag_oz2mz=0      ! XYAVG_DOC: $\left< \omega_z^2 \right>_{xy}$
  integer :: idiag_ruxmz=0      ! XYAVG_DOC: $\left<\varrho u_x \right>_{xy}$
  integer :: idiag_ruymz=0      ! XYAVG_DOC: $\left<\varrho u_y \right>_{xy}$
  integer :: idiag_ruzmz=0      ! XYAVG_DOC: $\left<\varrho u_z \right>_{xy}$
  integer :: idiag_rux2mz=0     ! XYAVG_DOC: $\left<\varrho u_x^2\right>_{xy}$
  integer :: idiag_ruy2mz=0     ! XYAVG_DOC: $\left<\varrho u_y^2\right>_{xy}$
  integer :: idiag_ruz2mz=0     ! XYAVG_DOC: $\left<\varrho u_z^2\right>_{xy}$
  integer :: idiag_uxuymz=0     ! XYAVG_DOC: $\left<u_x u_y\right>_{xy}$
  integer :: idiag_uxuzmz=0     ! XYAVG_DOC: $\left<u_x u_z\right>_{xy}$
  integer :: idiag_uyuzmz=0     ! XYAVG_DOC: $\left<u_y u_z\right>_{xy}$
  integer :: idiag_Rxymz=0      ! XYAVG_DOC: $\left<u_x' u_y'\right>_{xy}$
  integer :: idiag_Rxyupmz=0    ! XYAVG_DOC: $\left<u_{x\downarrow}' u_{y\downarrow}'\right>_{xy}$
  integer :: idiag_Rxydownmz=0  ! XYAVG_DOC: $\left<u_{x\uparrow}' u_{y\uparrow}'\right>_{xy}$
  integer :: idiag_Rxzmz=0      ! XYAVG_DOC: $\left<u_x' u_z'\right>_{xy}$
  integer :: idiag_Rxzupmz=0    ! XYAVG_DOC: $\left<u_{x\downarrow}' u_{z\downarrow}'\right>_{xy}$
  integer :: idiag_Rxzdownmz=0  ! XYAVG_DOC: $\left<u_{x\uparrow}' u_{z\uparrow}'\right>_{xy}$
  integer :: idiag_Ryzmz=0      ! XYAVG_DOC: $\left<u_y' u_z'\right>_{xy}$
  integer :: idiag_Ryzupmz=0    ! XYAVG_DOC: $\left<u_{y\downarrow}' u_{z\downarrow}'\right>_{xy}$
  integer :: idiag_Ryzdownmz=0  ! XYAVG_DOC: $\left<u_{y\uparrow}' u_{z\uparrow}'\right>_{xy}$
  integer :: idiag_ruxuymz=0    ! XYAVG_DOC: $\langle\rho u_x u_y\rangle_{xy}$
  integer :: idiag_ruxuzmz=0    ! XYAVG_DOC: $\langle\rho u_x u_z\rangle_{xy}$
  integer :: idiag_ruyuzmz=0    ! XYAVG_DOC: $\langle\rho u_y u_z\rangle_{xy}$
  integer :: idiag_ruxuy2mz=0   ! XYAVG_DOC: $\langle\rho u_x u_y\rangle_{xy}$
  integer :: idiag_ruxuz2mz=0   ! XYAVG_DOC: $\langle\rho u_x u_z\rangle_{xy}$
  integer :: idiag_ruyuz2mz=0   ! XYAVG_DOC: $\langle\rho u_y u_z\rangle_{xy}$
  integer :: idiag_oxuxxmz=0    ! XYAVG_DOC: $\left<\omega_x u_{x,x}\right>_{xy}$
  integer :: idiag_oyuxymz=0    ! XYAVG_DOC: $\left<\omega_y u_{x,y}\right>_{xy}$
  integer :: idiag_oxuyxmz=0    ! XYAVG_DOC: $\left<\omega_x u_{y,x}\right>_{xy}$
  integer :: idiag_oyuyymz=0    ! XYAVG_DOC: $\left<\omega_y u_{y,y}\right>_{xy}$
  integer :: idiag_oxuzxmz=0    ! XYAVG_DOC: $\left<\omega_x u_{z,x}\right>_{xy}$
  integer :: idiag_oyuzymz=0    ! XYAVG_DOC: $\left<\omega_y u_{z,y}\right>_{xy}$
  integer :: idiag_uyxuzxmz=0   ! XYAVG_DOC: $\left<u_{y,x} u_{z,x}\right>_{xy}$
  integer :: idiag_uyyuzymz=0   ! XYAVG_DOC: $\left<u_{y,y} u_{z,y}\right>_{xy}$
  integer :: idiag_uyzuzzmz=0   ! XYAVG_DOC: $\left<u_{y,z} u_{z,z}\right>_{xy}$
  integer :: idiag_ekinmz=0     ! XYAVG_DOC: $\left<{1\over2}\varrho\uv^2\right>_{xy}$
  integer :: idiag_oumz=0       ! XYAVG_DOC: $\left<\boldsymbol{\omega}
                                ! XYAVG_DOC: \cdot\uv\right>_{xy}$
  integer :: idiag_uguxmz=0     ! XYAVG_DOC:
  integer :: idiag_uguymz=0     ! XYAVG_DOC:
  integer :: idiag_uguzmz=0     ! XYAVG_DOC:
  integer :: idiag_Remz=0       ! XYAVG_DOC: $\langle\frac{|\uv\cdot\uv|}{\left|
                                ! XYAVG_DOC: \frac{\partial}{\partial x_j}
                                ! XYAVG_DOC: (\nu{\sf S}_{ij})\right|}\rangle_{xy}$
  integer :: idiag_oguxmz=0     ! XYAVG_DOC: $\left<(\boldsymbol{\omega}
                                ! XYAVG_DOC: \cdot\nabla \uv)_x\right>_{xy}$
  integer :: idiag_oguymz=0     ! XYAVG_DOC: $\left<(\boldsymbol{\omega}
                                ! XYAVG_DOC: \cdot\nabla \uv)_y\right>_{xy}$
  integer :: idiag_oguzmz=0     ! XYAVG_DOC: $\left<(\boldsymbol{\omega}
                                ! XYAVG_DOC: \cdot\nabla \uv)_z\right>_{xy}$
  integer :: idiag_ogux2mz=0    ! XYAVG_DOC: $\left<(\boldsymbol{\omega}
                                ! XYAVG_DOC: \cdot\nabla \uv)_x^2\right>_{xy}$
  integer :: idiag_oguy2mz=0    ! XYAVG_DOC: $\left<(\boldsymbol{\omega}
                                ! XYAVG_DOC: \cdot\nabla \uv)_y^2\right>_{xy}$
  integer :: idiag_oguz2mz=0    ! XYAVG_DOC: $\left<(\boldsymbol{\omega}
                                ! XYAVG_DOC:\cdot\nabla \uv)_z^2\right>_{xy}$
  integer :: idiag_oxdivumz=0   ! XYAVG_DOC: $\left<\omega_x\nabla\cdot\uv\right>_{xy}$
  integer :: idiag_oydivumz=0   ! XYAVG_DOC: $\left<\omega_y\nabla\cdot\uv\right>_{xy}$
  integer :: idiag_ozdivumz=0   ! XYAVG_DOC: $\left<\omega_z\nabla\cdot\uv\right>_{xy}$
  integer :: idiag_oxdivu2mz=0  ! XYAVG_DOC: $\left<(\omega_x\nabla\cdot\uv)^2\right>_{xy}$
  integer :: idiag_oydivu2mz=0  ! XYAVG_DOC: $\left<(\omega_y\nabla\cdot\uv)^2\right>_{xy}$
  integer :: idiag_ozdivu2mz=0  ! XYAVG_DOC: $\left<(\omega_z\nabla\cdot\uv)^2\right>_{xy}$
  integer :: idiag_u3u21mz=0    ! XYAVG_DOC:
  integer :: idiag_u1u32mz=0    ! XYAVG_DOC:
  integer :: idiag_u2u13mz=0    ! XYAVG_DOC:
  integer :: idiag_u2u31mz=0    ! XYAVG_DOC:
  integer :: idiag_u3u12mz=0    ! XYAVG_DOC:
  integer :: idiag_u1u23mz=0    ! XYAVG_DOC:
  integer :: idiag_acczmz=0     ! XYAVG_DOC: $\left<Du_z/Dt\right>_{xy}$
  integer :: idiag_acczupmz=0   ! XYAVG_DOC: $\left<Du_z/Dt\right>_{xy+}$
  integer :: idiag_acczdownmz=0 ! XYAVG_DOC: $\left<Du_z/Dt\right>_{xy-}$
  integer :: idiag_accpowzmz=0  ! XYAVG_DOC: $\left<(u_z Du_z/Dt)^2\right>_{xy}$
  integer :: idiag_accpowzupmz=0  ! XYAVG_DOC: $\left<(u_z Du_z/Dt)^2\right>_{xy+}$
  integer :: idiag_accpowzdownmz=0! XYAVG_DOC: $\left<(u_z Du_z/Dt)^2\right>_{xy-}$
  integer :: idiag_totalforcezmz=0     ! XYAVG_DOC: $\left<\varrho Du_z/Dt\right>_{xy}$
  integer :: idiag_totalforcezupmz=0   ! XYAVG_DOC: $\left<\varrho Du_z/Dt\right>_{xy+}$
  integer :: idiag_totalforcezdownmz=0 ! XYAVG_DOC: $\left<\varrho Du_z/Dt\right>_{xy-}$
!
! xz averaged diagnostics given in xzaver.in
!
  integer :: idiag_uxmy=0       ! XZAVG_DOC: $\left< u_x \right>_{xz}$
  integer :: idiag_uymy=0       ! XZAVG_DOC: $\left< u_y \right>_{xz}$
  integer :: idiag_uzmy=0       ! XZAVG_DOC: $\left< u_z \right>_{xz}$
  integer :: idiag_ux2my=0      ! XZAVG_DOC:
  integer :: idiag_uy2my=0      ! XZAVG_DOC:
  integer :: idiag_uz2my=0      ! XZAVG_DOC:
  integer :: idiag_uxuymy=0     ! XZAVG_DOC:
  integer :: idiag_uxuzmy=0     ! XZAVG_DOC:
  integer :: idiag_uyuzmy=0     ! XZAVG_DOC:
  integer :: idiag_oumy=0       ! XZAVG_DOC: $\left<\boldsymbol{\omega}
                                ! XZAVG_DOC: \cdot\uv\right>_{xz}$
  integer :: idiag_uguxmy=0     ! XZAVG_DOC:
  integer :: idiag_uguymy=0     ! XZAVG_DOC:
  integer :: idiag_uguzmy=0     ! XZAVG_DOC:
!
! yz averaged diagnostics given in yzaver.in
!
  integer :: idiag_uxmx=0       ! YZAVG_DOC: $\left< u_x \right>_{yz}$
  integer :: idiag_uymx=0       ! YZAVG_DOC: $\left< u_y \right>_{yz}$
  integer :: idiag_uzmx=0       ! YZAVG_DOC: $\left< u_z \right>_{yz}$
  integer :: idiag_ruxmx=0      ! YZAVG_DOC: $\left<\varrho u_x \right>_{yz}$
  integer :: idiag_ruymx=0      ! YZAVG_DOC: $\left<\varrho u_y \right>_{yz}$
  integer :: idiag_ruzmx=0      ! YZAVG_DOC: $\left<\varrho u_z \right>_{yz}$
  integer :: idiag_rux2mx = 0   ! YZAVG_DOC: $\langle\rho u_x^2\rangle_{yz}$
  integer :: idiag_ruy2mx = 0   ! YZAVG_DOC: $\langle\rho u_y^2\rangle_{yz}$
  integer :: idiag_ruz2mx = 0   ! YZAVG_DOC: $\langle\rho u_z^2\rangle_{yz}$
  integer :: idiag_ruxuymx = 0  ! YZAVG_DOC: $\langle\rho u_x u_y\rangle_{yz}$
  integer :: idiag_ruxuzmx = 0  ! YZAVG_DOC: $\langle\rho u_x u_z\rangle_{yz}$
  integer :: idiag_ruyuzmx = 0  ! YZAVG_DOC: $\langle\rho u_y u_z\rangle_{yz}$
  integer :: idiag_ux2mx=0      ! YZAVG_DOC: $\left<u_x^2\right>_{yz}$
  integer :: idiag_uy2mx=0      ! YZAVG_DOC: $\left<u_y^2\right>_{yz}$
  integer :: idiag_uz2mx=0      ! YZAVG_DOC: $\left<u_z^2\right>_{yz}$
  integer :: idiag_ox2mx=0      ! YZAVG_DOC: $\left<\omega_x^2\right>_{yz}$
  integer :: idiag_oy2mx=0      ! YZAVG_DOC: $\left<\omega_y^2\right>_{yz}$
  integer :: idiag_oz2mx=0      ! YZAVG_DOC: $\left<\omega_z^2\right>_{yz}$
  integer :: idiag_uxuymx=0     ! YZAVG_DOC: $\langle u_x u_y\rangle_{yz}$
  integer :: idiag_uxuzmx=0     ! YZAVG_DOC: $\langle u_x u_z\rangle_{yz}$
  integer :: idiag_uyuzmx=0     ! YZAVG_DOC: $\langle u_y u_z\rangle_{yz}$
  integer :: idiag_oumx=0       ! YZAVG_DOC: $\left<\boldsymbol{\omega}
                                ! YZAVG_DOC: \cdot\uv\right>_{yz}$
  integer :: idiag_uguxmx=0     ! YZAVG_DOC:
  integer :: idiag_uguymx=0     ! YZAVG_DOC:
  integer :: idiag_uguzmx=0     ! YZAVG_DOC:
  integer :: idiag_ekinmx=0     ! YZAVG_DOC: $\langle{1\over2}\rho u^2\rangle_{xy}$
  integer :: idiag_fkinxmx=0    ! XYAVG_DOC: $\left<{1\over2}\varrho\uv^2 u_x\right>_{yz}$
!
! y averaged diagnostics given in yaver.in
!
  integer :: idiag_uxmxz=0      ! YAVG_DOC: $\left< u_x \right>_{y}$
  integer :: idiag_uymxz=0      ! YAVG_DOC: $\left< u_y \right>_{y}$
  integer :: idiag_uzmxz=0      ! YAVG_DOC: $\left< u_z \right>_{y}$
  integer :: idiag_ux2mxz=0     ! YAVG_DOC: $\left< u_x^2 \right>_{y}$
  integer :: idiag_uy2mxz=0     ! YAVG_DOC: $\left< u_y^2 \right>_{y}$
  integer :: idiag_uz2mxz=0     ! YAVG_DOC: $\left< u_z^2 \right>_{y}$
  integer :: idiag_uxuymxz=0    ! YAVG_DOC: $\left< u_x u_y \right>_{y}$
  integer :: idiag_uxuzmxz=0    ! YAVG_DOC: $\left< u_x u_z \right>_{y}$
  integer :: idiag_uyuzmxz=0    ! YAVG_DOC: $\left< u_y u_z \right>_{y}$
  integer :: idiag_oumxz=0      ! YAVG_DOC: $\left<\boldsymbol{\omega}
                                ! YAVG_DOC: \cdot\uv\right>_{y}$
!
! z averaged diagnostics given in zaver.in
!
  integer :: idiag_uxmxy=0      ! ZAVG_DOC: $\left< u_x \right>_{z}$
  integer :: idiag_uymxy=0      ! ZAVG_DOC: $\left< u_y \right>_{z}$
  integer :: idiag_uzmxy=0      ! ZAVG_DOC: $\left< u_z \right>_{z}$
  integer :: idiag_uxupmxy=0    ! ZAVG_DOC: $\left< u_{x\uparrow} \right>_{z}$
  integer :: idiag_uxdownmxy=0  ! ZAVG_DOC: $\left< u_{x\downarrow} \right>_{z}$
  integer :: idiag_ruxupmxy=0   ! ZAVG_DOC: $\left<\rho u_{x\uparrow} \right>_{z}$
  integer :: idiag_ruxdownmxy=0 ! ZAVG_DOC: $\left<\rho u_{x\downarrow} \right>_{z}$
  integer :: idiag_ux2upmxy=0   ! ZAVG_DOC: $\left< u^2_{x\uparrow} \right>_{z}$
  integer :: idiag_ux2downmxy=0 ! ZAVG_DOC: $\left< u^2_{x\downarrow} \right>_{z}$
  integer :: idiag_ffdownmxy=0  ! ZAVG_DOC: Filling factor of downflows
  integer :: idiag_uxuymxy=0    ! ZAVG_DOC: $\left< u_x u_y \right>_{z}$
  integer :: idiag_uxuzmxy=0    ! ZAVG_DOC: $\left< u_x u_z \right>_{z}$
  integer :: idiag_uyuzmxy=0    ! ZAVG_DOC: $\left< u_y u_z \right>_{z}$
  integer :: idiag_Rxymxy=0     ! ZAVG_DOC: $\left<u_x' u_y'\right>_{z}$
  integer :: idiag_Rxyupmxy=0   ! ZAVG_DOC: $\left<(u_x' u_y')_\uparrow\right>_{z}$
  integer :: idiag_Rxydownmxy=0 ! ZAVG_DOC: $\left<(u_x' u_y')_\downarrow\right>_{z}$
  integer :: idiag_Rxzmxy=0     ! ZAVG_DOC: $\left<u_x' u_z'\right>_{z}$
  integer :: idiag_Rxzupmxy=0   ! ZAVG_DOC: $\left<(u_x' u_z')_\uparrow\right>_{z}$
  integer :: idiag_Rxzdownmxy=0 ! ZAVG_DOC: $\left<(u_x' u_z')_\downarrow\right>_{z}$
  integer :: idiag_Ryzmxy=0     ! ZAVG_DOC: $\left<u_y' u_z'\right>_{z}$
  integer :: idiag_Ryzupmxy=0   ! ZAVG_DOC: $\left<(u_y' u_z')_\uparrow\right>_{z}$
  integer :: idiag_Ryzdownmxy=0 ! ZAVG_DOC: $\left<(u_y' u_z')_\downarrow\right>_{z}$
  integer :: idiag_oxmxy=0      ! ZAVG_DOC: $\left< \omega_x \right>_{z}$
  integer :: idiag_oymxy=0      ! ZAVG_DOC: $\left< \omega_y \right>_{z}$
  integer :: idiag_ozmxy=0      ! ZAVG_DOC: $\left< \omega_z \right>_{z}$
  integer :: idiag_oumxy=0      ! ZAVG_DOC: $\left<\boldsymbol{\omega}
                                ! ZAVG_DOC: \cdot\uv\right>_{z}$
  integer :: idiag_pvzmxy=0     ! ZAVG_DOC: $\left< (\omega_z+2\Omega)/\varrho
                                ! ZAVG_DOC: \right>_{z}$ \quad(z component of
                                ! ZAVG_DOC: potential vorticity)
  integer :: idiag_uguxmxy=0    ! ZAVG_DOC: $\left< (\boldsymbol{u}\cdot\boldsymbol{\nabla} \boldsymbol{u})_x \right>_{z}$
  integer :: idiag_uguymxy=0    ! ZAVG_DOC: $\left< (\boldsymbol{u}\cdot\boldsymbol{\nabla} \boldsymbol{u})_y \right>_{z}$
  integer :: idiag_uguzmxy=0    ! ZAVG_DOC: $\left< (\boldsymbol{u}\cdot\boldsymbol{\nabla} \boldsymbol{u})_z \right>_{z}$
  integer :: idiag_ruxmxy=0     ! ZAVG_DOC: $\left< \rho u_x \right>_{z}$
  integer :: idiag_ruymxy=0     ! ZAVG_DOC: $\left< \rho u_y \right>_{z}$
  integer :: idiag_ruzmxy=0     ! ZAVG_DOC: $\left< \rho u_z \right>_{z}$
  integer :: idiag_ux2mxy=0     ! ZAVG_DOC: $\left< u_x^2 \right>_{z}$
  integer :: idiag_uy2mxy=0     ! ZAVG_DOC: $\left< u_y^2 \right>_{z}$
  integer :: idiag_uz2mxy=0     ! ZAVG_DOC: $\left< u_z^2 \right>_{z}$
  integer :: idiag_rux2mxy=0    ! ZAVG_DOC: $\left< \rho u_x^2 \right>_{z}$
  integer :: idiag_ruy2mxy=0    ! ZAVG_DOC: $\left< \rho u_y^2 \right>_{z}$
  integer :: idiag_ruz2mxy=0    ! ZAVG_DOC: $\left< \rho u_z^2 \right>_{z}$
  integer :: idiag_ruxuymxy=0   ! ZAVG_DOC: $\left< \rho u_x u_y \right>_{z}$
  integer :: idiag_ruxuzmxy=0   ! ZAVG_DOC: $\left< \rho u_x u_z \right>_{z}$
  integer :: idiag_ruyuzmxy=0   ! ZAVG_DOC: $\left< \rho u_y u_z \right>_{z}$
  integer :: idiag_fkinxmxy=0   ! ZAVG_DOC: $\left<{1\over2}\varrho\uv^2
                                ! ZAVG_DOC: u_x\right>_{z}$
  integer :: idiag_fkinymxy=0   ! ZAVG_DOC: $\left<{1\over2}\varrho\uv^2
                                ! ZAVG_DOC: u_y\right>_{z}$
  integer :: idiag_fkinxupmxy=0 ! ZAVG_DOC: $\left<{1\over2}\varrho\uv^2
                                ! ZAVG_DOC: u_{x\uparrow}\right>_{z}$
  integer :: idiag_fkinxdownmxy=0 ! ZAVG_DOC: $\left<{1\over2}\varrho\uv^2
                                ! ZAVG_DOC: u_{x\downarrow}\right>_{z}$
  integer :: idiag_nshift=0
!
!  Video data.
!
  integer :: ivid_oo=0, ivid_o2=0, ivid_divu=0, ivid_u2=0, ivid_Ma2=0, ivid_uu_sph=0
  integer :: ivid_ou=0
!
!  Auxiliary variables
!
  real, dimension(:,:), pointer :: reference_state
  real, dimension(3) :: Omegav=0.
  real, dimension(nx) :: Fmax,advec_uu=0.
  real :: t_vart=0.
!$omp THREADPRIVATE(advec_uu)
!
  real, dimension (nx) :: prof_amp1, prof_amp2
  real, dimension (mz) :: prof_amp3
  real, dimension (my) :: prof_amp4
  real, dimension (nz,3) :: uumz_prof
  real, dimension (nx,ny) :: omega_prof

  contains
!***********************************************************************
    subroutine register_hydro
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
! 18-may-12/MR: put Pr*Ra as a shared variable for use in
!               temperature_idealgas
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
      use Sub, only: register_report_aux
!
!  indices to access uu
!
      call farray_register_pde('uu',iuu,vector=3)
      iux = iuu; iuy = iuu+1; iuz = iuu+2
      if (llinearized_hydro) then
        call farray_register_auxiliary('uu0',iuu0,vector=3)
        iu0x = iuu0; iu0y = iuu0+1; iu0z = iuu0+2
      endif
!
!  Register an extra aux slot for uut if requested. This is needed
!  for calculating the correlation time from <u.intudt>. For this to work
!  you must reserve enough auxiliary workspace by setting, for example,
!     ! MAUX CONTRIBUTION 3
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!  29-oct-20/hongzhe: added uust
!
      if (luut_as_aux) then
        call register_report_aux('uut', iuut, iuxt, iuyt, iuzt)
        ltime_integrals=.true.
      endif
      if (luust_as_aux) then
        call register_report_aux('uust', iuust, iuxst, iuyst, iuzst)
        ltime_integrals=.true.
      endif
      if (loot_as_aux) then
        call register_report_aux('oot', ioot, ioxt, ioyt, iozt)
        ltime_integrals=.true.
      endif
      if (loost_as_aux) then
        call register_report_aux('oost', ioost, ioxst, ioyst, iozst)
        ltime_integrals=.true.
      endif
!
!  omega as aux
!
      !if (loo_as_aux) call register_report_aux('oo', ioo, iox, ioy, ioz, communicated=.true.)
      if (loo_as_aux) call register_report_aux('oo', ioo, iox, ioy, ioz)
!!
!!  Fourier transformed uu as aux
!!
!!     if (luuk_as_aux) call register_report_aux('uuk', iuuk)
!!     if (look_as_aux) call register_report_aux('ook', iook)
!
!  To compute the added mass term for particle drag,
!  the advective derivative is needed.
!  The advective derivative is set as an auxiliary
!
      if (ladv_der_as_aux) then
        call farray_register_auxiliary('adv_der_uu',i_adv_der,vector=3)
        i_adv_derx = i_adv_der;  i_adv_dery = i_adv_der+1; i_adv_derz = i_adv_der+2
      endif
!
!  Velocity fluctuation for computation of turbulent Reynolds stress etc.
!  directly in the code.
!
      if (luu_fluc_as_aux) then
        if (iuu_fluc==0) then
          call farray_register_auxiliary('uu_fluc',iuu_fluc,vector=3)
          iuu_flucx = iuu_fluc; iuu_flucy = iuu_fluc+1; iuu_flucz = iuu_fluc+2;
        else
          if (lroot) print*, 'initialize_hydro: iuu_fluc = ', iuu_fluc
          call farray_index_append('iuu_fluc',iuu_fluc)
        endif
      endif
!
!  Velocity in spherical coordinates from Cartesian simulation, useful
!  for sphere-in-a-box models.
!
      if (luu_sph_as_aux.and.lsphere_in_a_box) then
        if (iuu_sph==0) then
          call farray_register_auxiliary('uu_sph',iuu_sph,vector=3)
          iuu_sphr = iuu_sph; iuu_spht = iuu_sph+1; iuu_sphp = iuu_sph+2;
        else
          if (lroot) print*, 'initialize_hydro: iuu_sph = ', iuu_sph
          call farray_index_append('iuu_sph',iuu_sph)
        endif
      endif
!
!  For Helmholtz decomposition of uu the potential of the curl-free part is registered  as an auxiliary.
!
      if (lhelmholtz_decomp) then
        if (dsnap_down==0.) call fatal_error('register_hydro','Helmholtz decomposition requires dsnap_down>0')
        call farray_register_auxiliary('phiuu',iphiuu)
      endif
!
!  Share lpressuregradient_gas so the entropy module knows whether to apply
!  pressure gradient or not. But hydro wants pressure gradient only when
!  the density is computed, i.e. not with lboussinesq nor lanelastic.
!
      if  (.not.ldensity.or.lanelastic) lpressuregradient_gas=.false.
      call put_shared_variable('lpressuregradient_gas', &
          lpressuregradient_gas,caller='register_hydro')
!
!  Special settings for lboussinesq.
!
      if  (lboussinesq) then
        PrRa=Pr*Ra
        call put_shared_variable('PrRa',PrRa)
        call put_shared_variable('Pr',Pr)
      endif
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Writing files for use with IDL.
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',uu $'
          if (nvar == mvar) write(4,*) ',uu'
        else
          write(4,*) ',uu $'
        endif
        write(15,*) 'uu = fltarr(mx,my,mz,3)*one'
      endif
!
!  shared variable of lrelativistic for density
!
      call put_shared_variable('lrelativistic',lrelativistic, caller='register_hydro')
!
! If we are to solve for gradient of dust particle velocity, we must store gradient
! of gas velocity as auxiliary
!
      if (lparticles_grad) lgradu_as_aux=.true.
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine initialize_hydro(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!  13-oct-03/dave: check parameters and warn (if nec.) about velocity damping
!  26-mar-10/axel: lreinitialize_uu added
!  23-dec-15/MR: Cartesian vector Omegav intro'd; ltime_integrals set; rotation
!                of \vec{Omega} on Yang grid added.
!   7-jun.16/MR: modifications for calculation of z average on Yin-Yang grid, not yet operational
!
      use BorderProfiles, only: request_border_driving
      use Initcond
      use SharedVariables, only: put_shared_variable,get_shared_variable
      use Sub, only: step, erfunc, register_report_aux
      use Slices_methods, only: alloc_slice_buffers
      use Yinyang_mpi, only: initialize_zaver_yy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: c, s
      integer :: j,myl ! currently unused: nycap
      real :: slope,uinn,uext,zbot
!
! set the right point in profile to unity.
!
      profile_SH(index_rSH)=dx_1(index_rSH)
!
!  Block use of uninitalised p%fcont
!
      if (.not.lforcing_cont) lforcing_cont_uu=.false.
!
!  Calculate cosz*sinz, cos^2, and sinz^2, to take moments with
!  of ux2, uxuy, etc.
!
      c=cos(kz_analysis*z)
      s=sin(kz_analysis*z)
      cz=c
      sz=s
      c2z=c**2
      s2z=s**2
!
!  Rescale velocity by a factor reinitialize_uu.
!
      if (lreinitialize_uu) then
        do j=1,ninit
          select case (inituu(j))
          case ('Beltrami-z'); call beltrami(ampluu(j),f,iuu,kz=kz_uu)
          case ('rescale'); f(:,:,:,iux:iuz)=rescale_uu*f(:,:,:,iux:iuz)
          case ('gaussian-noise'); call gaunoise(ampluu(j),f,iux,iuz)
          case ('gaussian-noise-z'); call gaunoise(ampluu(j),f,iuz)
          case ('no-uy'); f(:,:,:,iuy)=0.
          case ('flip-ux'); f(:,:,:,iux)=-f(:,:,:,iux)
          case ('flip-uy'); f(:,:,:,iuy)=-f(:,:,:,iuy)
          case ('mult-uz-lower-xbdry'); if (ipx==0) f(1:l1,:,:,iuz)=rescale_uu*f(1:l1,:,:,iuz)
          case ('Om_inner'); f(:,:,:,iuz)=Om_inner*spread(spread(xyz0(1)**2/x,2,my)*spread(sin(y),1,mx),3,mz)
          endselect
        enddo
      endif
!
      ! Default value of 'tfade_start' is tdamp/2 for faded damping
      if (.not. ldamp_fade .and. (tfade_start >= 0.0) .and. (tdamp > 0.0)) ldamp_fade = .true.
      if (ldamp_fade .and. (tfade_start == -1.0)) tfade_start = 0.5 * tdamp
      if (ldamp_fade .and. (tfade_start >= tdamp) .and. (tdamp > 0.0)) &
          call fatal_error ('initialize_hydro', 'Please set tfade_start < tdamp')
      call put_shared_variable ('dampu', dampu, caller='initialize_hydro')
      call put_shared_variable ('tdamp', tdamp)
      call put_shared_variable ('ldamp_fade', ldamp_fade)
      call put_shared_variable ('tfade_start', tfade_start)
!
!  r_int and r_ext override rdampint and rdampext if both are set
!
      if (dampuint /= 0.) then
        if (r_int > epsi) then
          rdampint = r_int
        elseif (rdampint <= epsi) then
          write(*,*) 'initialize_hydro: inner radius not yet set, dampuint= ',dampuint
        endif
      endif

      if (Omega/=0.) then
!
!  defining an r-dependent profile for Omega. The Coriolis force will be suppressed
!  in r < r_omega with the width w_omega, for having the supression for r> r_omega,
!  choose a negativ w_omega.
!
        prof_om = 1.0
        if (r_omega /= 0.0) &
          prof_om = step(x(l1:l2),r_omega,w_omega)
!
!  Cartesian components of \vec{Omega} (not yet used).
!
        if (lyang) then
          Omegav(1) = -Omega*sin(theta*dtor)*cos(phi*dtor)
          Omegav(2) = -Omega*cos(theta*dtor)
          Omegav(3) = -Omega*sin(theta*dtor)*sin(phi*dtor)
        else
          Omegav(1) = Omega*sin(theta*dtor)*cos(phi*dtor)
          Omegav(2) = Omega*sin(theta*dtor)*sin(phi*dtor)
          Omegav(3) = Omega*cos(theta*dtor)
        endif
      endif
!
!  damping parameters for damping velocities outside an embedded sphere
!  04-feb-2008/dintrans: corrected because otherwise rdampext=r_ext all the time
!
      if (dampuext /= 0.0) then
!       if (r_ext < impossible) then
!         rdampext = r_ext
!       elseif (rdampext == impossible) then
!         write(*,*) 'initialize_hydro: outer radius not yet set, dampuext= ',dampuext
!       endif
        if (rdampext == impossible) then
          if (r_ext < impossible) then
            if (lroot) write(*,*) 'initialize_hydro: set outer radius rdampext=r_ext'
            rdampext = r_ext
          else
            if (lroot) write(*,*) 'initialize_hydro: cannot set outer radius rdampext=r_ext'
          endif
        else
          if (lroot) write(*,*) 'initialize_hydro: outer radius rdampext=',rdampext
        endif
      endif
!
!  Calculate inverse damping times for damping momentum in the
!  x and y directions.
!
      if (tau_damp_ruxm /= 0.) tau_damp_ruxm1=1./tau_damp_ruxm
      if (tau_damp_ruym /= 0.) tau_damp_ruym1=1./tau_damp_ruym
      if (tau_damp_ruzm /= 0.) tau_damp_ruzm1=1./tau_damp_ruzm
!
!  Set freezing arrays.
!
      if (lfreeze_uint) lfreeze_varint(iux:iuz) = .true.
      if (lfreeze_uext) lfreeze_varext(iux:iuz) = .true.
!
!  Turn off advection for 0-D runs.
!
      if (nwgrid==1) then
        ladvection_velocity=.false.
        if (lroot) print*, &
             'initialize_hydro: 0-D run, turned off advection of velocity'
      endif
!
!  Border profile backward compatibility. For a vector, if only the first
!  borderuu is set, then the other components get the same value.
!
      if (lpropagate_borderuu     .and. &
           borderuu(1)/='nothing' .and. &
           borderuu(2)=='nothing' .and. &
           borderuu(3)=='nothing') then
        borderuu(2)=borderuu(1)
        borderuu(3)=borderuu(1)
      endif
!
!  Tell the BorderProfiles module if we intend to use border driving, so
!  that the module can request the right pencils.
!
      do j=1,3
        if (borderuu(j)/='nothing') call request_border_driving(borderuu(j))
      enddo
!
!  Hand over Coriolis force to Particles_drag.
!
      drag: if (lparticles_drag .and. lcoriolis_force) then
        lcoriolis_force = .false.
        if (lroot) print *, 'initialize_hydro: turned off and hand over Coriolis force to Particles_drag. '
      endif drag
!
      lshear_in_coriolis=lshear_in_coriolis.and.lcoriolis_force.and.lshear
!
!  Compute mask for x-averaging where x is in hydro_xaver_range.
!  Normalize such that the average over the full domain
!  gives still unity.
!
      if (l1 == l2) then
        xmask_hyd = 1.
      else
        where (x(l1:l2) >= hydro_xaver_range(1) .and. x(l1:l2) <= hydro_xaver_range(2))
          xmask_hyd = 1.
        elsewhere
          xmask_hyd = 0.
        endwhere
        hydro_xaver_range(1) = max(hydro_xaver_range(1), xyz0(1))
        hydro_xaver_range(2) = min(hydro_xaver_range(2), xyz1(1))
        if (lspherical_coords) then
          xmask_hyd = xmask_hyd * (xyz1(1)**3 - xyz0(1)**3) &
              / (hydro_xaver_range(2)**3 - hydro_xaver_range(1)**3)
        elseif (lcylindrical_coords) then
          xmask_hyd = xmask_hyd * (xyz1(1)**2 - xyz0(1)**2) &
              / (hydro_xaver_range(2)**2 - hydro_xaver_range(1)**2)
        else
          xmask_hyd = xmask_hyd*Lxyz(1) &
              / (hydro_xaver_range(2) - hydro_xaver_range(1))
        endif
      endif
!
!  Compute mask for z-averaging where z is in hydro_zaver_range.
!  Normalize such that the average over the full domain
!  gives still unity.
!
      if (n1 == n2) then
        zmask_hyd = 1.
      else
        where (z(n1:n2) >= hydro_zaver_range(1) .and. z(n1:n2) <= hydro_zaver_range(2))
          zmask_hyd = 1.
        elsewhere
          zmask_hyd = 0.
        endwhere
        hydro_zaver_range(1) = max(hydro_zaver_range(1), xyz0(3))
        hydro_zaver_range(2) = min(hydro_zaver_range(2), xyz1(3))
        zmask_hyd = zmask_hyd*Lxyz(3)/(hydro_zaver_range(2) - hydro_zaver_range(1))
      endif
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'xmask_hyd=',xmask_hyd
        print*,'zmask_hyd=',zmask_hyd
      endif
!
      if (force_lower_bound == 'vel_time' .or. force_upper_bound == 'vel_time') then
        call put_shared_variable('ampl_forc', ampl_forc)
        call put_shared_variable('k_forc', k_forc)
        call put_shared_variable('w_forc', w_forc)
        call put_shared_variable('x_forc', x_forc)
        call put_shared_variable('dx_forc', dx_forc)
      endif

      call put_shared_variable('lshear_rateofstrain',lshear_rateofstrain)
!
!  Check if we are solving the force-free equations in parts of domain.
!  This is currently possible with density (including anelastic, but not
!  boussinesq).
!
      if (ldensity) then
        call get_shared_variable('lffree',lffree)
        if (lffree) then
          call get_shared_variable('profx_ffree',profx_ffree)
          call get_shared_variable('profy_ffree',profy_ffree)
          call get_shared_variable('profz_ffree',profz_ffree)
        endif
      endif
!
      if (lviscosity) call put_shared_variable ('lcalc_uuavg',lcalc_uuavg)
!
!  Get the reference state if requested
!
      if (lreference_state) &
        call get_shared_variable('reference_state',reference_state)
!
      lcalc_uumeanz = lcalc_uumeanz .or. lcalc_uumean .or. ltestfield_xz .or. &      ! lcalc_uumean for compatibility
                      lremove_uumeanz .or. lremove_uumeanz_horizontal
      lcalc_uumeanx = lcalc_uumeanx.or.lremove_uumeanx
      lcalc_uumeany = lcalc_uumeany.or.lremove_uumeany
      lcalc_uumeanxy = lcalc_uumeanxy .or. lremove_uumeanxy.or.ltestfield_xy
!
      if (lremove_uumeanz.or.lremove_uumeanx.or.lremove_uumeany) then
        if (lremove_mean_flow) call warning('initialize_hydro', &
            'lremove_mean_flow=T may interfere with lremove_uumean[xyz]=T')
        if (lremove_mean_momenta) call warning('initialize_hydro', &
            'lremove_mean_momenta=T may interfere with lremove_uumean[xyz]=T')
        if (lremove_uumeanxy) call warning('initialize_hydro', &
            'lremove_uumeanxy=T may interfere with lremove_uumean[xyz]=T')
      endif
!
      if (Omega/=0. .and. lyinyang) then
        if (phi==0.) then
!
!  Rotate \vec{Omega} on Yang grid.
!
          if (lyang) then
            Omega = -Omega
            phi   = 90.-theta
            theta = 90.
          endif
        else
          call fatal_error('initialize_hydro', 'phi /= 0. not allowed for Yin-Yang grid')
        endif
      endif
!
      if (lcalc_uumeanxy .or. lcalc_ruumeanxy) then
        myl=my
        if (lyinyang) then
          call fatal_error('initialize_hydro','Calculation of z average not implmented for Yin-Yang')
          !call initialize_zaver_yy(myl,nycap)
        endif
        allocate(uumxy(mx,myl,3))
        allocate(ruumxy(mx,myl,3))
        uumxy=0.0
        ruumxy=0.0
      endif
!
!  Preparations for adding/removing mean flows.
!  Set profiles for forcing differential rotation.
!
      select case (uuprof)

      case ('BS04')
        if (wdamp/=0.) then
          prof_amp1=1.-step(x(l1:l2),rdampint,wdamp)
        else
          prof_amp1=1.
        endif
        prof_amp1=ampl1_diffrot*prof_amp1*cos(kx_diffrot*x(l1:l2))**xexp_diffrot
        prof_amp3=cos(z)

      case ('BS04c','BS04c1','HP09')

        if (wdamp/=0.) then
          prof_amp3=ampl1_diffrot*0.5*(1.+tanh((z-rdampint)/(wdamp)))
        else
          prof_amp3=ampl1_diffrot
        endif

        if (uuprof=='BS04c') then
          prof_amp1=sin(0.5*pi*((x(l1:l2))-x0)/Lx)**xexp_diffrot
        elseif (uuprof=='BS04c1') then
          prof_amp1=sin(pi*((x(l1:l2))-x0)/Lx)**xexp_diffrot
        elseif(uuprof=='HP09') then
          prof_amp1=cos(kx_diffrot*x(l1:l2))
!or       prof_amp1=cos(2.*pi*kx_diffrot*(x(l1:l2)-x0)/Lx)
        endif

      case ('BS04m')
        if (wdamp/=0.) then
          prof_amp1=1.-step(x(l1:l2),rdampint,wdamp)
        else
          prof_amp1=1.
        endif
        prof_amp1=ampl1_diffrot*prof_amp1*sin((pi/(2.*x(l2)))*x(l1:l2))
        prof_amp4=cos(pi/(2.*y(m2))*y)

      case ('solar_DC99')
        prof_amp1=(1.-ampl1_diffrot*step(x(l1:l2),rdampext,wdamp))*step(x(l1:l2),rdampint,wdamp)*x(l1:l2)
        prof_amp4=ampl2_diffrot*(1.064-0.145*costh**2-0.155*costh**4-1.)*sinth

      case ('vertical_shear')
        zbot=xyz0(3)
        prof_amp3=ampl1_diffrot*cos(kz_diffrot*(z-zbot)-phase_diffrot)

      case ('vertical_compression','vertical_shear_x')
        zbot=xyz0(3)
        prof_amp3=ampl1_diffrot*cos(kz_diffrot*(z-zbot))

      case ('remove_vertical_shear')
        if (.not.lcalc_uumean) &
          call fatal_error("initialize_hydro","you need to set lcalc_uumean=T for uuprof='remove_vertical_shear'")

      case ('vertical_shear_x_sinz')
        zbot=xyz0(3)
        where (z <= 0.) 
          prof_amp3=ampl1_diffrot*sin(.5*pi/abs(zbot)*z)
        elsewhere
          prof_amp3=0.
        endwhere

      case ('vertical_shear_z')
        prof_amp3=ampl1_diffrot*tanh((z-rdampint)/width_ff_uu)

      case ('vertical_shear_z2')
        if (.not.lcalc_uumeanxz) &
          call fatal_error("initialize_hydro","you need to set lcalc_uumeanxz=T for uuprof='vertical_shear_z2'")
        prof_amp3=ampl1_diffrot*tanh((z-rdampint)/width_ff_uu)

      case ('vertical_shear_linear')
        if (.not.lcalc_uumeanxz) &
          call fatal_error("initialize_hydro","you need to set lcalc_uumeanxz=T for uuprof='vertical_shear_linear'")
        prof_amp3=ampl1_diffrot*z

      case ('tachocline')
        if (wdamp/=0.) then
          prof_amp1=1.-step(x(l1:l2),rdampint,wdamp)
        else
          prof_amp1=1.
        endif
      case ('solar_simple')
        if (lspherical_coords) then
          prof_amp1=ampl1_diffrot*step(x(l1:l2),x1_ff_uu,width_ff_uu)
          prof_amp4=1.5-7.5*costh*costh
        elseif (lcartesian_coords) then
          prof_amp1=ampl1_diffrot*cos(x(l1:l2))
          prof_amp4=cos(y)*cos(y)
        !prof_amp2=1.-step(x(l1:l2),x2_ff_uu,width_ff_uu)
        else
          call fatal_error("initialize_hydro", &
                          "uuprof='solar_simple' not implemented for other than spherical or Cartesian coordinates") 
        endif
      case ('radial_uniform_shear')
        uinn = omega_in*x(l1)
        uext = omega_out*x(l2)
        slope = (uext - uinn)/(x(l2)-x(l1))
        prof_amp1=slope*x(l1:l2)+(uinn*x(l2)- uext*x(l1))/(x(l2)-x(l1))

      case ('breeze')
        prof_amp3=ampl_wind*z/(2.*pi)

      case ('slow_wind')
        prof_amp3=ampl_wind*(1.+tanh((z-rdampext)/wdamp))

      case ('radial_shear')
        if (.not.lcalc_uumeanxy) &
          call fatal_error("initialize_hydro","you need to set lcalc_uumeanxy=T for uuprof='radial_shear'")
        prof_amp1=ampl1_diffrot*cos(2*pi*k_diffrot*(x(l1:l2)-x0)/Lx)

      case ('radial_shear_damp')
        if (.not.lcalc_uumeanxy) &
          call fatal_error("initialize_hydro","you need to set lcalc_uumeanxy=T for uuprof='radial_shear_damp'")
        prof_amp1=ampl1_diffrot*tanh((x(l1:l2)-rdampint)/wdamp)

      case ('damp_corona')
        if (lspherical_coords) then
          if (.not.lcalc_uumeanxy) &
            call fatal_error("initialize_hydro","you need to set lcalc_uumeanxy=T for uuprof='damp_corona'")
          prof_amp1=0.5*(tanh((x(l1:l2)-rdampext)/wdamp)+1.)
        elseif (lcartesian_coords) then
          prof_amp3=0.5*(tanh((z-zbot)/wdamp)+1.)
        endif

      case ('damp_horiz_vel')
        prof_amp3=0.5*(tanh((z-rdampext)/wdamp)+1.)

      case ('latitudinal_shear')
        if (.not.lcalc_uumeanxy) &
          call fatal_error("initialize_hydro","you need to set lcalc_uumeanxy=T for uuprof='latitudinal_shear'")
        prof_amp4=ampl1_diffrot*cos(2.*pi*k_diffrot*(y-y0)/Ly)

      case ('damp_jets')
        if (.not.lcalc_uumeanxy) &
          call fatal_error("initialize_hydro","you need to set lcalc_uumeanxy=T for uuprof='damp_jets'")
        prof_amp4=1.-0.5*(1.+tanh((y-(y0+ydampint))/wdamp)-(1.+tanh((y-(y0+Lxyz(2)-ydampext))/wdamp)))

      case ('spoke-like-NSSL')
        if (.not.lspherical_coords) &
          call warning("initialize_hydro", &
                       "uuprof='spoke-like-NSSL' only meningful for spherical coordinates")
        if (.not.lcalc_uumeanxy) &
          call fatal_error("initialize_hydro","you need lcalc_uumeanxy=T for uuprof='spoke-like-NSSL'")

        prof_amp1=ampl1_diffrot*x(l1:l2)
        profx_diffrot1=+0.5*(1.+erfunc(((x(l1:l2)-uphi_rbot)/uphi_step_width)))
        profx_diffrot2=+0.5*(1.-erfunc(((x(l1:l2)-uphi_rtop)/uphi_step_width)))
        profx_diffrot3=+0.5*(1.+erfunc(((x(l1:l2)-uphi_rtop)/uphi_step_width)))
        profx_diffrot2=(x(l1:l2)-uphi_rbot)*profx_diffrot1*profx_diffrot2 !(redefined)
        profy_diffrot1=-1.5*(5.*costh**2-1.)
        profy_diffrot2=-1.0*(4.*costh**2-3.)
        profy_diffrot3=-1.0
        profz_diffrot1=+1.
!
      case ('uumz_profile')
        if (.not.lcalc_uumeanz) then
          call fatal_error("initialize_hydro","you need to set lcalc_uumean=T for uuprof='uumz_profile'")
        else
          if (.not.lgravz) &
            call fatal_error("initialize_hydro","gravitation in z-direction (lgravz=T) needed for uuprof='uumz_profile'")
          call read_uumz_profile(uumz_prof)
        endif

      case ('omega_profile')
        if (.not.lspherical_coords) &
          call warning("initialize_hydro", &
                       "uuprof='omega_profile' only meaningful for spherical coordinates")
        if (.not.lcalc_uumeanxy) &
          call fatal_error("initialize_hydro","you need to set lcalc_uumeanxy=T for uuprof='omega_profile'")
        call read_omega_profile(omega_prof)

      case ('nothing')

      case default
         if (lroot) print*,'initialize_hydro: No profile of mean flow "'//trim(uuprof)//'"'

      endselect
!
      if (ivid_oo/=0) call alloc_slice_buffers(oo_xy,oo_xz,oo_yz,oo_xy2,oo_xy3,oo_xy4,oo_xz2,oo_r)
      if (ivid_o2/=0) call alloc_slice_buffers(o2_xy,o2_xz,o2_yz,o2_xy2,o2_xy3,o2_xy4,o2_xz2,o2_r)
      if (ivid_ou/=0) call alloc_slice_buffers(ou_xy,ou_xz,ou_yz,ou_xy2,ou_xy3,ou_xy4,ou_xz2,ou_r)
      if (ivid_u2/=0) call alloc_slice_buffers(u2_xy,u2_xz,u2_yz,u2_xy2,u2_xy3,u2_xy4,u2_xz2,u2_r)
      if (ivid_divu/=0) call alloc_slice_buffers(divu_xy,divu_xz,divu_yz,divu_xy2,divu_xy3,divu_xy4,divu_xz2,divu_r)
      if (ivid_Ma2/=0) call alloc_slice_buffers(mach_xy,mach_xz,mach_yz,mach_xy2,mach_xy3,mach_xy4,mach_xz2,mach_r)
      if (ivid_uu_sph/=0) &
        call alloc_slice_buffers(uu_sph_xy,uu_sph_xz,uu_sph_yz,uu_sph_xy2,uu_sph_xy3,uu_sph_xy4,uu_sph_xz2,uu_sph_r)
!
!  give warning if orms is not set in prints.in
!
      if (othresh_per_orms/=0..and.idiag_orms==0) then
        if (lroot) print*,'calc_othresh: need to set orms in print.in to get othresh.'
      endif
!
      call keep_compiler_quiet(f)
!
      endsubroutine initialize_hydro
!***********************************************************************
      subroutine calc_means_hydro(f)
!
! calculates various means
!
!  14-oct-13/MR: outsourced from hydro_after_boundary
!  13-feb-15/MR: changes for use of reference_state
!
      use Mpicomm, only: mpiallreduce_sum
      use Sub, only: finalize_aver
      use Deriv, only: der_z
      use DensityMethods, only: getrho
      use Yinyang_mpi, only: zsum_yy
!
      real, dimension (mx,my,mz,mfarray), intent(IN) :: f
!
      real, dimension (nx) :: rho,rux,ruy,ruz
      integer, parameter :: nreduce=3
      real, dimension (nreduce) :: fsum_tmp,fsum
      real, dimension (:,:,:), allocatable :: buffer
      real :: fact
      integer :: j,nnz,l,m,n
!
!  calculate averages of rho*ux and rho*uy
!
      if (ldensity) then
        if (tau_damp_ruxm/=0. .or. tau_damp_ruym/=0. .or. tau_damp_ruzm/=0.) then
          ruxm=0.
          ruym=0.
          ruzm=0.
          fact=1./nwgrid
          do n=n1,n2
          do m=m1,m2
            call getrho(f(:,m,n,ilnrho),rho)
            rux=rho*f(l1:l2,m,n,iux)
            ruy=rho*f(l1:l2,m,n,iuy)
            ruz=rho*f(l1:l2,m,n,iuz)
            ruxm=ruxm+fact*sum(rux)
            ruym=ruym+fact*sum(ruy)
            ruzm=ruzm+fact*sum(ruz)
          enddo
          enddo
!
!  communicate to the other processors
!
          fsum_tmp(1)=ruxm
          fsum_tmp(2)=ruym
          fsum_tmp(3)=ruzm
          call mpiallreduce_sum(fsum_tmp,fsum,nreduce)
          ruxm=fsum(1)
          ruym=fsum(2)
          ruzm=fsum(3)
        endif
      endif
!
!  do xy-averaged mean flow for each component
!
      if (lcalc_uumeanz) then
        fact=1./nxygrid
        do nnz=1,mz
          do j=1,3
            uumz(nnz,j)=fact*sum(f(l1:l2,m1:m2,nnz,iux+j-1))
          enddo
        enddo
        call finalize_aver(nprocxy,12,uumz)
      endif
!
!  do xy-averaged mean momentum for each component
!
      if (lcalc_ruumeanz) then
        fact=1./nxygrid
        do nnz=1,mz
          do j=1,3
            ruumz(nnz,j)=fact*sum(exp(f(l1:l2,m1:m2,nnz,ilnrho))*f(l1:l2,m1:m2,nnz,iux+j-1))
          enddo
        enddo
        call finalize_aver(nprocxy,12,ruumz)
!
      endif
!
!  do yz-averaged mean flow for each component
!
      if (lcalc_uumeanx) then
        fact=1./nyzgrid
        do l=1,mx
          do j=1,3
            uumx(l,j)=fact*sum(f(l,m1:m2,n1:n2,iux+j-1))
          enddo
        enddo
        call finalize_aver(nprocyz,23,uumx)
      endif
!
!  do xz-averaged mean flow for each component
!
      if (lcalc_uumeany) then
        fact=1./nxzgrid
        do m=1,my
          do j=1,3
            uumy(m,j)=fact*sum(f(l1:l2,m,n1:n2,iux+j-1))
          enddo
        enddo
        call finalize_aver(nprocxz,13,uumy)
      endif
!
!  Do mean 2D field in (x,y)-plane for each component
!
      if (lcalc_uumeanxy) then
!
        fact=1./nzgrid_eff
        if (lyang) allocate(buffer(1,mx,my))
        do j=1,3
          if (lyang) then
!
!  On Yang grid:
!
            do n=n1,n2
              do m=1,my
                call zsum_yy(buffer,1,m,n,f(:,m,n,iuu+j-1))
              enddo
            enddo
            uumxy(:,:,j)=fact*buffer(1,:,:)
          else
!
! Normal summing-up in Yin procs.
!
            uumxy(:,:,j)=fact*sum(f(:,:,n1:n2,iuu+j-1),3)  ! requires equidistant grid
          endif
        enddo
        call finalize_aver(nprocz,3,uumxy)
!
      endif
!
!  Do mean 2D momentum in (x,y)-plane for each component
!
      if (lcalc_ruumeanxy) then
!
        fact=1./nzgrid_eff
        if (lyang) allocate(buffer(1,mx,my))
        do j=1,3
          if (lyang) then
!
!  On Yang grid:
!
            do n=n1,n2
              do m=1,my
                call zsum_yy(buffer,1,m,n,exp(f(:,m,n,ilnrho))*f(:,m,n,iuu+j-1))
              enddo
            enddo
            ruumxy(:,:,j)=fact*buffer(1,:,:)
          else
!
! Normal summing-up in Yin procs.
!
            ruumxy(:,:,j)=fact*sum(exp(f(:,:,n1:n2,ilnrho))*f(:,:,n1:n2,iuu+j-1),3)  ! requires equidistant grid
          endif
        enddo
        call finalize_aver(nprocz,3,ruumxy)
!
      endif
!
!  Do mean 2D field in (x,z)-plane for each component
!
      if (lcalc_uumeanxz) then
!
        uumxz = sum(f(:,m1:m2,:,iux:iuz),2)/nygrid
        call finalize_aver(nprocy,2,uumxz)
!
      endif
!
!  do communication for array of size nz*nprocz*3*njtest
!
!     if (nprocy>1) then
!       uum2=reshape(uumz1,shape=(/nz*nprocz*3/))
!       call mpireduce_sum(uumz2,uum3,(/nz,nprocz,3/))
!       call mpibcast_real(uumz3,nz*nprocz*3)
!       uum1=reshape(uum3,shape=(/nz,nprocz,3/))
!       do n=n1,n2
!         do j=1,3
!           uumz(n,j)=uumz1(n-n1+1,ipz+1,j)
!         enddo
!       enddo
!     endif
!
      endsubroutine calc_means_hydro
!***********************************************************************
    subroutine init_uu(f)
!
!  initialise velocity field ; called from start.f90
!
!  07-nov-01/wolf: coded
!  24-nov-02/tony: renamed for consistence (i.e. init_[variable name])
!  13-feb-15/MR: changes for use of reference_state
!
      use Boundcond, only: update_ghosts
      use Density, only: beta_glnrho_scaled
      use DensityMethods, only: getrho, putlnrho
      use EquationOfState, only: cs20
      use General
      use Gravity, only: gravz_const,z1
      use Initcond
      use InitialCondition, only: initial_condition_uu
      use Sub
      use Mpicomm, only: lyang
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (3) :: tmpvec
!
      real, dimension (nx,3) :: tmp_nx3
      real, dimension (mx) :: tmpmx
      real, dimension (nx) :: r,p1,tmp,prof,xc0,yc0,ur,lnrhor
      real, dimension (:,:), allocatable :: yz
      real :: kabs,crit,eta_sigma,tmp0,phi0
      real :: a2, rr2, wall_smoothing
      real :: dis, xold,yold,uprof, factx, factz, sph, sph_har_der, der
      integer :: j,i,l,ixy,ix,iy,iz,iz0,iyz,iter,niter=100
      logical :: lvectorpotential=.false.
!
!  inituu corresponds to different initializations of uu (called from start).
!
      do j=1,ninit
!
        select case (inituu(j))
!
        case ('nothing'); if (lroot .and. j==1) print*,'init_uu: nothing'
        case ('zero', '0')
          if (lroot) print*,'init_uu: zero velocity'
          ! Ensure really is zero, as may have used lread_oldsnap
          f(:,:,:,iux:iuz)=0.
        case ('const_uu','const-uu'); do i=1,3; f(:,:,:,iuu+i-1) = uu_const(i); enddo
        case('smooth_step_ux')
          xhalf= 0.5*(xyz1(1)+xyz0(1))
          do iy=m1,m2;do iz=n1,n2
            f(:,iy,iz,iuy)= -ampluu(j)+2.*ampluu(j)*step(x,xhalf,widthuu)
          enddo;enddo
        case('parabola_x')
          do ix=l1,l2;do iy=m1,m2;do iz=n1,n2
            f(ix,iy,iz,iuu)=0
            f(ix,iy,iz,iuu+1)=max_uu*(1-(x(ix)/x(l1))**2)
            f(ix,iy,iz,iuu+2)=0
          enddo;enddo;enddo
        case ('mode'); call modev(ampluu(j),coefuu,f,iuu,kx_uu,ky_uu,kz_uu)
        case ('ortho')
          do ix=l1,l2;do iy=m1,m2;do iz=n1,n2
            f(ix,iy,iz,iuu)=-2.*ampluu(j)*sin(y(iy))
            f(ix,iy,iz,iuu+1)=ampluu(j)*sin(x(ix))
            f(ix,iy,iz,iuu+2)=0.001*ampluu(j)*sin(y(iy))
          enddo;enddo;enddo
        case ('Kolmogorov-x')
          do iy=m1,m2; do iz=n1,n2
            f(:,iy,iz,iuy)=ampluu(j)*cos(kx_uu*x)
          enddo; enddo
        case ('Kolmogorov-x-perturb')
          print*,'ampluu,kx_uu,amp_factor,kx_uu_perturb',ampluu,kx_uu,amp_factor,kx_uu_perturb
          do iy=m1,m2; do iz=n1,n2
            f(:,iy,iz,iuy)=ampluu(j)*(cos(kx_uu*x)+amp_factor*cos(kx_uu_perturb*x))
          enddo; enddo
        case ('uxsinx-deltaz')
          iz0=nint((zsphere-xyz0(3))/dz)
          print*, 'Ux=sin(kx*x)\delta(z)'
          print*, 'ampluu,kx_uu,zsphere,iz0',ampluu,kx_uu,zsphere,iz0
          print*, 'dz=',dz
          f(:,:,:,iux:iuz) = 0.
          do ix=l1,l2
             if ((iz0 .le. n2) .and. (iz0 .gt. n1) )  &
             f(ix,:,iz0,iux) = ampluu(j)*sin(kx_uu*x(ix))
          enddo
        case ('random_isotropic_shell')
          call random_isotropic_shell(f,iux,ampluu(j),z1_uu,z2_uu)
        case ('gaussian-noise'); call gaunoise(ampluu(j),f,iux,iuz)
        case ('gaussian-noise-x'); call gaunoise(ampluu(j),f,iux)
        case ('gaussian-noise-y'); call gaunoise(ampluu(j),f,iuy)
        case ('gaussian-noise-z'); call gaunoise(ampluu(j),f,iuz)
        case ('gaussian-noise-xy'); call gaunoise(ampluu(j),f,iux,iuy)
        case ('gaussian-noise-rprof')
          call gaunoise_rprof(ampluu(j),f,iux,iuz,rnoise_int,rnoise_ext)
        case ('xjump')
          call jump(f,iux,uu_left,uu_right,widthuu,'x')
          call jump(f,iuy,uy_left,uy_right,widthuu,'x')
        case ('Beltrami-x'); call beltrami(ampluu(j),f,iuu,kx=kx_uu,sigma=relhel_uu)
        case ('Beltrami-y'); call beltrami(ampluu(j),f,iuu,ky=ky_uu,sigma=relhel_uu)
        case ('Beltrami-z'); call beltrami(ampluu(j),f,iuu,kz=kz_uu,sigma=relhel_uu)
        case ('Straining'); call straining(ampluu(j),f,iuu,kx_uu,ky_uu,kz_uu,dimensionality)
        case ('rolls'); call rolls(ampluu(j),f,iuu,kx_uu,kz_uu)
        case ('trilinear-x'); call trilinear(f,iux,ampl_ux(j),ampl_uy(j),ampl_uz(j))
        case ('trilinear-y'); call trilinear(f,iuy,ampl_ux(j),ampl_uy(j),ampl_uz(j))
        case ('trilinear-z'); call trilinear(f,iuz,ampl_ux(j),ampl_uy(j),ampl_uz(j))
        case ('cos-cos-sin-uz'); call cos_cos_sin(ampluu(j),f,iuz)
        case ('tor_pert'); call tor_pert(ampluu(j),f,iux)
        case ('rotblob'); call rotblob(ampluu(j),incl_alpha,f,iux,rot_rr,xsphere,ysphere,zsphere)
        case ('rotblob_yz'); call rotblob_yz(ampluu(j),f,iux,rot_rr,xsphere,ysphere,zsphere)
        case ('read_arr_file'); call read_outside_vec_array(f, "uu.arr", iuu)
        case ('diffrot'); call diffrot(ampluu(j),f,iuy)
        case ('olddiffrot'); call olddiffrot(ampluu(j),f,iuy)
        case ('sinwave-phase')
          call sinwave_phase(f,iux,ampl_ux(j),kx_ux(j),ky_ux(j),kz_ux(j),phase_ux(j))
          call sinwave_phase(f,iuy,ampl_uy(j),kx_uy(j),ky_uy(j),kz_uy(j),phase_uy(j))
          call sinwave_phase(f,iuz,ampl_uz(j),kx_uz(j),ky_uz(j),kz_uz(j),phase_uz(j))
        case ('coswave-phase')
          call coswave_phase(f,iux,ampl_ux(j),kx_ux(j),ky_ux(j),kz_ux(j),phase_ux(j))
          call coswave_phase(f,iuy,ampl_uy(j),kx_uy(j),ky_uy(j),kz_uy(j),phase_uy(j))
          call coswave_phase(f,iuz,ampl_uz(j),kx_uz(j),ky_uz(j),kz_uz(j),phase_uz(j))
        case ('sinwave-x'); call sinwave(ampluu(j),f,iux,kx=kx_uu)
        case ('sinwave-y'); call sinwave(ampluu(j),f,iuy,ky=ky_uu)
        case ('sinwave-z'); call sinwave(ampluu(j),f,iuz,kz=kz_uu)
        case ('sinwave-ux-kx'); call sinwave(ampluu(j),f,iux,kx=kx_uu)
        case ('sinwave-ux-ky'); call sinwave(ampluu(j),f,iux,ky=ky_uu)
        case ('sinwave-ux-kz'); call sinwave(ampluu(j),f,iux,kz=kz_uu)
        case ('sinwave-uy-kx'); call sinwave(ampluu(j),f,iuy,kx=kx_uu)
        case ('sinwave-uy-ky'); call sinwave(ampluu(j),f,iuy,ky=ky_uu)
        case ('sinwave-uy-kz'); call sinwave(ampluu(j),f,iuy,kz=kz_uu)
        case ('sinwave-uz-kx'); call sinwave(ampluu(j),f,iuz,kx=kx_uu)
        case ('sinwave-uz-ky'); call sinwave(ampluu(j),f,iuz,ky=ky_uu)
        case ('sinwave-uz-kz'); call sinwave(ampluu(j),f,iuz,kz=kz_uu)
        case ('sinwave-y-z')
          if (lroot) print*, 'init_uu: sinwave-y-z, ampluu=', ampluu(j)
          call sinwave(ampluu(j),f,iuy,kz=kz_uu)
        case ('sinwave-z-y')
          if (lroot) print*, 'init_uu: sinwave-z-y, ampluu=', ampluu(j)
          call sinwave(ampluu(j),f,iuz,ky=ky_uu)
        case ('sinwave-z-x')
          if (lroot) print*, 'init_uu: sinwave-z-x, ampluu=', ampluu(j)
          call sinwave(ampluu(j),f,iuz,kx=kx_uu)
        case ('damped_sinwave-z-x')
          if (lroot) print*, 'init_uu: damped_sinwave-z-x, ampluu=', ampluu(j)
          do m=m1,m2; do n=n1,n2
            f(:,m,n,iuz)=f(:,m,n,iuz)+ampluu(j)*sin(kx_uu*x)*exp(-10*z(n)**2)
          enddo; enddo
        !case ('hatwave-x'); call hatwave(ampluu(j),f,iux,widthuu,kx=kx_uu,power=initpower)
        case ('45deg-sinwave-x-y')
          if (lroot) print*, 'init_uu: 45deg_sinwave-x-y, ampluu=', ampluu(j)
          do m=m1,m2; do n=n1,n2
            f(:,m,n,iux)=f(:,m,n,iux)+ampluu(j)*sin(kx_uu*x+ky_uu*y(m))
            f(:,m,n,iuy)=f(:,m,n,iuy)+ampluu(j)*sin(kx_uu*x+ky_uu*y(m))
          enddo; enddo
        case ('45deg-sinwave-y-z')
          if (lroot) print*, 'init_uu: 45deg_sinwave-y-z, ampluu=', ampluu(j)
          do m=m1,m2; do n=n1,n2
            f(:,m,n,iuy)=f(:,m,n,iuy)+ampluu(j)*sin(ky_uu*y(m)+kz_uu*z(n))
            f(:,m,n,iuz)=f(:,m,n,iuz)+ampluu(j)*sin(ky_uu*y(m)+kz_uu*z(n))
          enddo; enddo
        case ('coswave-x'); call coswave(ampluu(j),f,iux,kx=kx_uu,ky=ky_uu,kz=kz_uu)
        case ('coswave-y'); call coswave(ampluu(j),f,iuy,kx=kx_uu,ky=ky_uu,kz=kz_uu)
        case ('coswave-z'); call coswave(ampluu(j),f,iuz,kz=kz_uu)
        case ('coswave-x-z'); call coswave(ampluu(j),f,iux,kz=kz_uu)
        case ('coswave-z-x'); call coswave(ampluu(j),f,iuz,kx=kx_uu)
        case ('x1cosycosz'); call x1_cosy_cosz(ampluu(j),f,iuy,ky=ky_uu,kz=kz_uu)
        case ('couette'); call couette(ampluu(j),mu_omega,f,iuy)
        case ('couette_rings'); call couette_rings(ampluu(j),mu_omega,nb_rings,om_rings,gap,f,iuy)
        case ('soundwave-x'); call soundwave(ampluu(j),f,iux,kx=kx_uu,width=widthuu)
        case ('soundwave-y'); call soundwave(ampluu(j),f,iuy,ky=ky_uu)
        case ('soundwave-z'); call soundwave(ampluu(j),f,iuz,kz=kz_uu)
        case ('robertsflow'); call robertsflow(ampluu(j),f,iuu,relhel_uu)
        case ('hawley-et-al'); call hawley_etal99a(ampluu(j),f,iuy,Lxyz)
        case ('meri_circ'); call meri_circ(f)
        case ('sound-wave', '11')
!
!  sound wave (should be consistent with density module)
!
          if (lroot) print*,'init_uu: x-wave in uu; ampluu(j)=',ampluu(j)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iux)=uu_const(1)+ampluu(j)*sin(kx_uu*x(l1:l2))
          enddo; enddo
!
        case ('sound-wave2')
!
!  sound wave (should be consistent with density module)
!
          crit=cs20-gravz_const/kx_uu**2
          if (lroot) print*,'init_uu: x-wave in uu; crit,ampluu(j)=',crit,ampluu(j)
          do n=n1,n2; do m=m1,m2
            if (crit>0.) then
              f(l1:l2,m,n,iux)=+ampluu(j)*cos(kx_uu*x(l1:l2))*sqrt(abs(crit))
            else
              f(l1:l2,m,n,iux)=-ampluu(j)*sin(kx_uu*x(l1:l2))*sqrt(abs(crit))
            endif
          enddo; enddo
!
        case ('ABC')
          if (headtt) print*,'ABC flow'
! uu
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iux)=ampluu(j)*(ABC_A*sin(kz_uu*z(n))    +ABC_C*cos(ky_uu*y(m))    )
            f(l1:l2,m,n,iuy)=ampluu(j)*(ABC_B*sin(kx_uu*x(l1:l2))+ABC_A*cos(kz_uu*z(n))    )
            f(l1:l2,m,n,iuz)=ampluu(j)*(ABC_C*sin(ky_uu*y(m))    +ABC_B*cos(kx_uu*x(l1:l2)))
          enddo; enddo
!
        case ('double_sine')
          if (headtt) print*,'double sine flow'
! uu
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iux)=ampluu(j)*(ABC_A*sin(kz_uu*z(n))    +ABC_C*sin(ky_uu*y(m))    )
            f(l1:l2,m,n,iuy)=ampluu(j)*(ABC_B*sin(kx_uu*x(l1:l2))+ABC_A*sin(kz_uu*z(n))    )
            f(l1:l2,m,n,iuz)=ampluu(j)*(ABC_C*sin(ky_uu*y(m))    +ABC_B*sin(kx_uu*x(l1:l2)))
          enddo; enddo
!
        case ('shock-tube', '13')
!
!  shock tube test (should be consistent with density module)
!
          if (lroot) print*,'init_uu: polytopic standing shock'
          do n=n1,n2; do m=m1,m2
            prof=.5*(1.+tanh(x(l1:l2)/widthuu))
            f(l1:l2,m,n,iux)=uu_left+(uu_right-uu_left)*prof
          enddo; enddo
!
        case ('tanhx')
!
!  Burgers shock
!
          if (lroot) print*,'init_uu: Burgers shock'
          prof=-ampluu(j)*tanh(.5*x(l1:l2)/widthuu)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iux)=prof
          enddo; enddo
!
        case ('tanhy')
!
!  A*tanh(y/d) profile in the x-direction
!
          if (lroot) print*,'init_uu: tangential discontinuity'
          do m=m1,m2
            prof=ampluu(j)*tanh(y(m)/widthuu)
            do n=n1,n2
              f(l1:l2,m,n,iux)=prof
            enddo
          enddo
!
        case ('double_shear_layer')
!
!  Double shear layer
!
          if (lroot) print*,'init_uu: Double shear layer'
          do m=m1,m2
            der=2./delta_u
            prof=ampluu(j)*(&
                tanh(der*(y(m)+widthuu))-&
                tanh(der*(y(m)-widthuu)))/2.
            do n=n1,n2
              f(l1:l2,m,n,iux)=prof
            enddo
          enddo
!
        case ('cos_exp_perturbed')

          if (lroot) print*,'init_uu: tanhy_perturbed'
          do l=l1,l2;  do m=m1,m2
            uprof=amp_factor*cos(kx_uu*x(l))*exp(-abs(y(m))/widthuu)
            do n=n1,n2
              f(l,m,n,iuy)=uprof
              f(l,m,n,iux)=0.0
            enddo
          enddo; enddo
!
        case ('sin_exp_perturbed')

          if (lroot) print*,'init_uu: tanhy_perturbed'
          do l=l1,l2;  do m=m1,m2
            uprof=amp_factor*sin(kx_uu*x(l))*exp(-abs(y(m))/widthuu)
            do n=n1,n2
              f(l,m,n,iuy)=uprof
              f(l,m,n,iux)=0.0
            enddo
          enddo; enddo
!
        case ('tanhy_perturbed')
!
!  A*( tanh(y/d)+amp_factor*sin(kx_uu*x)*exp(-abs(y/d)) ) profile in the x-direction
!
          if (lroot) print*,'init_uu: tanhy_perturbed'
          do l=l1,l2;  do m=m1,m2
            uprof=ampluu(j)*tanh(y(m)/widthuu)
            do n=n1,n2
              f(l,m,n,iux)=uprof
            enddo
            uprof=amp_factor*sin(kx_uu*x(l))*exp(-abs(y(m))/widthuu)
!            uprof=0.0001*sin(94.2477796*x(l))*exp(-0.1*abs(y(m))/widthuu)
            do n=n1,n2
              f(l,m,n,iuy)=uprof
            enddo
          enddo; enddo
!
        case ('shock-sphere')
!
!  shock tube test (should be consistent with density module)
!
          if (lroot) print*,'init_uu: spherical shock, widthuu=',widthuu,' radiusuu=',radiusuu
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iux)=0.5*x(l1:l2)/radiusuu*ampluu(j)*(1.-tanh((sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)-radiusuu)/widthuu))
            f(l1:l2,m,n,iuy)=0.5*y(m)/radiusuu*ampluu(j)*(1.-tanh((sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)-radiusuu)/widthuu))
            f(l1:l2,m,n,iuz)=0.5*z(n)/radiusuu*ampluu(j)*(1.-tanh((sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)-radiusuu)/widthuu))
          enddo; enddo
!
!
        case ('bullets')
!
!  blob-like velocity perturbations (bullets)
!
          if (lroot) print*,'init_uu: velocity blobs'
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)+ampluu(j)*exp(-(x(l1:l2)**2+y(m)**2+z(n)**2)/widthuu)
          enddo; enddo
!
!  blob-like velocity perturbations in x-direction (bullets)
!
        case ('bullets_x')
          if (lroot) print*,'init_uu: velocity blobs in x-direction'
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iux)=uu_const(1)+f(l1:l2,m,n,iux) &
              +ampluu(j)*exp(-((x(l1:l2)-xsphere)**2+ &
                                   (y(m)-ysphere)**2+ &
                                   (z(n)-zsphere)**2)/widthuu**2)
          enddo; enddo
!
!  blob-like velocity perturbations in x-direction (bullets)
!
        case ('bullets_xz_plane')
          factx=ampluu(j)*sin(uu_xz_angle(j)*pi/180.)
          factz=ampluu(j)*cos(uu_xz_angle(j)*pi/180.)
          if (lroot) print*,'init_uu: velocity blobs in xz-plane'
          do n=n1,n2; do m=m1,m2
            tmp=exp(-((x(l1:l2)-xsphere)**2+ &
              (y(m)-ysphere)**2+(z(n)-zsphere)**2)/widthuu**2)
            f(l1:l2,m,n,iux)=f(l1:l2,m,n,iux)+factx*tmp
            f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)+factz*tmp
          enddo; enddo
!
!  X-point, xy plane
!
        case ('x-point_xy')
          if (lroot) print*,'init_uu: X-point'
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iux)=f(l1:l2,m,n,iux)+ampluu(j)*x(l1:l2)
            f(l1:l2,m,n,iuy)=f(l1:l2,m,n,iuy)-ampluu(j)*y(m)
          enddo; enddo
!
!  X-point, xz plane
!
        case ('x-point_xz')
          if (lroot) print*,'init_uu: X-point'
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iux)=f(l1:l2,m,n,iux)+ampluu(j)*x(l1:l2)
            f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)-ampluu(j)*z(n)
          enddo; enddo
!
        case ('Alfven-circ-x')
!
!  circularly polarised Alfven wave in x direction
!
          do n=n1,n2; do m=m1,m2
          if (lroot) print*,'init_uu: circular Alfven wave -> x'
            f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + ampluu(j)*sin(kx_uu*x(l1:l2))
            f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + ampluu(j)*cos(kx_uu*x(l1:l2))
          enddo; enddo
!
        case ('coszsiny-uz')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)- &
                ampluu(j)*cos(pi*z(n)/Lxyz(3))*sin(2*pi*y(m)/Lxyz(2))
          enddo; enddo
!
        case ('siny/x2')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iuy)=f(l1:l2,m,n,iuy)+ &
                ampluu(j)*sin(y(m))/x(l1:l2)**2
          enddo; enddo
!
        case ('linear-shear')
!
!  Linear shear
!
          if (lroot) print*,'init_uu: linear-shear, ampluu=', ampluu(j)
          do l=l1,l2; do m=m1,m2
            f(l,m,n1:n2,iuy) = ampluu(j)*z(n1:n2)
          enddo; enddo
!
        case ('linear-x_shear-uy')
          if (lroot) print*,'init_uu: linear-x_shear-uy, ampluu=', ampluu(j)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iuy) = ampluu(j)*x(l1:l2)
          enddo; enddo
!
        case ('tanh-x-z')
          if (lroot) print*, &
              'init_uu: tanh-x-z, widthuu, ampluu=', widthuu, ampluu(j)
          do l=l1,l2; do m=m1,m2
            f(l,m,n1:n2,iux) = ampluu(j)*tanh(z(n1:n2)/widthuu)
          enddo; enddo
!
        case ('tanh-y-z')
          if (lroot) print*, &
              'init_uu: tanh-y-z, widthuu, ampluu=', widthuu, ampluu(j)
          do l=l1,l2; do m=m1,m2
            f(l,m,n1:n2,iuy) = ampluu(j)*tanh(z(n1:n2)/widthuu)
          enddo; enddo
!
        case ('gauss-x-z')
          if (lroot) print*, &
              'init_uu: gauss-x-z, widthuu, ampluu=', widthuu, ampluu(j)
          do l=l1,l2; do m=m1,m2
            f(l,m,n1:n2,iux) = ampluu(j)*exp(-z(n1:n2)**2/widthuu**2)
          enddo; enddo
!
        case ('const-ux')
!
!  constant x-velocity
!
          if (lroot) print*,'init_uu: constant x-velocity'
          f(:,:,:,iux) = ampluu(j)
!
        case ('const-uy')
!
!  constant y-velocity
!
          if (lroot) print*,'init_uu: constant y-velocity'
          f(:,:,:,iuy) = ampluu(j)
!
        case ('const-uz')
!
!  constant z-velocity
!
          if (lroot) print*,'init_uu: constant z-velocity'
          f(:,:,:,iuz) = ampluu(j)
!
        case ('omega-z')
!
!  constant x-velocity
!
          if (lroot) print*,'init_uu: constant angular velocity omega_ini=',omega_ini
          f(:,:,:,iux) = 0
          f(:,:,:,iuy) = 0
          do n=n1,n2
            do m=m1,m2
              f(l1:l2,m,n,iuz) = omega_ini*x(l1:l2)*sinth(m)
            enddo
          enddo
!
        case ('tang-discont-z')
!
!  tangential discontinuity: velocity is directed along x,
!  ux=uu_lower for z<0 and ux=uu_upper for z>0. This can
!  be set up together with 'rho-jump' in density.
!
          if (lroot) print*,'init_uu: tangential discontinuity of uux at z=0'
          if (lroot) print*,'init_uu: uu_lower=',uu_lower,' uu_upper=',uu_upper
          if (lroot) print*,'init_uu: widthuu=',widthuu
          do n=n1,n2; do m=m1,m2
            prof=.5*(1.+tanh(z(n)/widthuu))
            f(l1:l2,m,n,iux)=uu_lower+(uu_upper-uu_lower)*prof
!
!  Add some random noise to see the development of instability
!WD: Can't we incorporate this into the urand stuff?
            print*, 'init_uu: ampluu(j)=',ampluu(j)
            call random_number_wrapper(r)
            call random_number_wrapper(p1)
!          tmp=sqrt(-2*log(r))*sin(2*pi*p1)*exp(-z(n)**2*10.)
            tmp=exp(-z(n)**2*10.)*cos(2.*x(l1:l2)+sin(4.*x(l1:l2)))
            f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)+ampluu(j)*tmp
          enddo; enddo
!
        case ('Fourier-trunc')
!
!  truncated simple Fourier series as nontrivial initial profile
!  for convection. The corresponding stream function is
!    exp(-(z-z1)^2/(2w^2))*(cos(kk)+2*sin(kk)+3*cos(3kk)),
!    with kk=k_x*x+k_y*y
!  Not a big success (convection starts much slower than with
!  random or 'up-down') ..
!
          if (lroot) print*,'init_uu: truncated Fourier'
          do n=n1,n2; do m=m1,m2
            prof = ampluu(j)*exp(-0.5*(z(n)-z1)**2/widthuu**2)!vertical Gaussian
            tmp = kx_uu*x(l1:l2) + ky_uu*y(m)                 ! horizontal phase
            kabs = sqrt(kx_uu**2+ky_uu**2)
            f(l1:l2,m,n,iuz) = prof * kabs*(-sin(tmp) + 4*cos(2*tmp) - 9*sin(3*tmp))
            tmp = (z(n)-z1)/widthuu**2*prof*(cos(tmp) + 2*sin(2*tmp) + 3*cos(3*tmp))
            f(l1:l2,m,n,iux) = tmp*kx_uu/kabs
            f(l1:l2,m,n,iuy) = tmp*ky_uu/kabs
          enddo; enddo
!
        case ('up-down')
!
!  flow upwards in one spot, downwards in another; not soneloidal
!
          if (lroot) print*,'init_uu: up-down'
          do n=n1,n2; do m=m1,m2

            prof = ampluu(j)*exp(-0.5*(z(n)-z1)**2/widthuu**2) ! vertical profile

            tmp = sqrt((x(l1:l2)-(x0+0.3*Lx))**2+(y(m)-(y0+0.3*Ly))**2)! dist. from spot 1
            f(l1:l2,m,n,iuz) = prof*exp(-0.5*(tmp**2)/widthuu**2)

            tmp = sqrt((x(l1:l2)-(x0+0.5*Lx))**2+(y(m)-(y0+0.8*Ly))**2)! dist. from spot 1
            f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) - 0.7*prof*exp(-0.5*(tmp**2)/widthuu**2)

          enddo; enddo
!
        case ('powern')
! initial spectrum k^power
          call powern(ampluu(j),initpower,cutoff,f,iux,iuz)
!
! initial spectrum k^power
!
        case ('power_randomphase')
          call power_randomphase(ampluu(j),initpower,kgaussian_uu,kpeak,cutoff,f,iux,iuz,lscale_tobox)
!
! initial spectrum k^power
!
        case ('power_randomphase_hel')
          call power_randomphase_hel(ampluu(j),initpower,initpower2, &
            cutoff,ncutoff,kpeak,f,iux,iuz,relhel_uu,kgaussian_uu, &
            lskip_projection, lvectorpotential,lscale_tobox, &
            nfact0=nfact_uu, lfactors0=lfactors_uu)
!
        case ('random-isotropic-KS')
          call random_isotropic_KS(initpower,f,iux,N_modes_uu)
!
        case ('vortex_2d')
! Vortex solution of Goodman, Narayan, & Goldreich (1987)
          call vortex_2d(f,b_ell,widthuu,rbound)
!
        case ('sub-Keplerian')
          if (lroot) print*, 'init_hydro: set sub-Keplerian gas velocity'
          f(:,:,:,iux) = f(:,:,:,iux) - 1/(2*Omega)*cs20*beta_glnrho_scaled(2)
          f(:,:,:,iuy) = f(:,:,:,iuy) + 1/(2*Omega)*cs20*beta_glnrho_scaled(1)
          ! superimpose here for the case of pressure bump special module chaning f as well
!
        case ('rigid')
          do n=n1,n2; do m=m1,m2; do l=l1,l2
            if (x(l)**2+y(m)**2+z(n)**2<=radiusuu**2) then
              f(l,m,n,iux)=-ampluu(j)*y(m)
              f(l,m,n,iuy)=+ampluu(j)*x(l)
            endif
          enddo; enddo; enddo
!
!  compressive (non-vortical) shear wave of Johnson & Gammie (2005a)
!
        case ('compressive-shwave')
          if (ldensity.or.lanelastic) then
            call coswave_phase(f,iux,ampl_ux(j),kx_uu,ky_uu,kz_uu,phase_ux(j))
            call coswave_phase(f,iuy,ampl_uy(j),kx_uu,ky_uu,kz_uu,phase_uy(j))
            eta_sigma = (2. - qshear)*Omega
            do n=n1,n2; do m=m1,m2
              tmp = -kx_uu*ampl_uy(j)*eta_sigma* &
                    (cos(kx_uu*x(l1:l2)+ky_uu*y(m)+kz_uu*z(n)) + &
                     sin(kx_uu*x(l1:l2)+ky_uu*y(m)+kz_uu*z(n)))
              call putlnrho(f(:,m,n,ilnrho),tmp)
            enddo; enddo
          endif
        case ( 'random-2D-eddies')
          if (lroot) &
            print*, "random-2D-eddies: ampluu,kx_uu,ky_uu = ", ampluu(j),kx_uu,ky_uu
          f(:,:,:,iuz)=0.
          call random_number_wrapper(xc0)
! Introduce both counter clockwise and clockwise eddies
          do ixy=1,neddy
            if (xc0(ixy)<=0.5) then
              tmp(ixy)=-1.0
            else
              tmp(ixy)=1.0
            endif
          enddo
!
          call random_number_wrapper(xc0)
          xc0=(1.-2*xc0)*Lxyz(1)/2
          call random_number_wrapper(yc0)
          yc0=(1.-2*yc0)*Lxyz(2)/2
! need to initialize xold, yold
! bing: suggestion use pos of last eddy
          xold = xc0(neddy)
          yold = yc0(neddy)
!
          do n=n1,n2; do m=m1,m2
! Check for nearest neighbour eddies and change their sign
            do ixy=1,neddy
              dis=sqrt((xold-xc0(ixy))**2+(yold-yc0(ixy))**2)
              if (dis<5*sqrt(1./kx_uu**2+1./ky_uu**2)) then
                tmp(ixy)=-tmp(ixy-1)
                if (lroot) &
                  write(*,*) 'init_uu: random-2D-eddies have come very close!'
              endif
              f(l1:l2,m,n,iuz)=f(l1:l2,m,n,iuz)+tmp(ixy)*ampluu(j) &
              *exp(-kx_uu*(x(l1:l2)-xc0(ixy))**2-ky_uu*(y(m)-yc0(ixy))**2) &
              *exp(-kz_uu*z(n)**2)
!
              xold=xc0(ixy)
              yold=yc0(ixy)
            enddo
          enddo; enddo
          call update_ghosts(f)
! 2D curl
          do n=n1,n2;do m=m1,m2
            call grad(f,iuz,tmp_nx3)
            f(l1:l2,m,n,iux) = -tmp_nx3(:,2)
            f(l1:l2,m,n,iuy) =  tmp_nx3(:,1)
          enddo;enddo
          f(:,:,:,iuz)=0.
          do m=m1,m2
            call random_number_wrapper(f(l1:l2,m,n1,iuz))
            do n=n1,n2
              if (louinit) then
                f(l1:l2,m,n,iuz)=100*ampluu(j)*(2*f(l1:l2,m,n1,iuz)-1)
              else
                f(l1:l2,m,n,iuz)=0.0d0
              endif
            enddo
          enddo
!  Transformation-reflection x -> -x and ux -> -ux
          if (lreflecteddy) then
           do iz=1,mz; do iy=1,my;do ix=1, mx/2
              tmpvec = f(mx-ix+1,iy,iz,iux:iuz)
              f(mx-ix+1,iy,iz,iux)= -f(ix,iy,iz,iux)
              f(ix,iy,iz,iux)=-tmpvec(1)
              f(mx-ix+1,iy,iz,iuy)= f(ix,iy,iz,iuy)
              f(ix,iy,iz,iuy)=tmpvec(2)
              f(mx-ix+1,iy,iz,iuz)= f(ix,iy,iz,iuz)
              f(ix,iy,iz,iuz)=tmpvec(3)
            enddo; enddo; enddo
          endif
          close(15)
!
        case ( 'anelastic-nlin')
          print*, "anelastic-2dxz: ampl_uy,kx_uu,kz_uu = ", ampl_uy(j),kx_uu,ky_uu,kz_uu
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iuy)=ampl_uy(j)*exp(-kx_uu*x(l1:l2)**2-kz_uu*z(n)**2)
          enddo; enddo
          call update_ghosts(f)
! 2D curl
          do n=n1,n2;do m=m1,m2
            call grad(f,iuy,tmp_nx3)
            call getrho(f(:,m,n,ilnrho),tmp)
            f(l1:l2,m,n,iux) = -tmp_nx3(:,3)/tmp
            f(l1:l2,m,n,iuz) =  tmp_nx3(:,1)/tmp
          enddo;enddo
          f(:,:,:,iuy)=0.
!
        case ( 'anelastic-lin')
          print*, "anelastic-2dxz: ampl_ux,kx_uu,kz_uu = ", ampl_ux(j),kx_uu,kz_uu
          f(:,:,:,iuy)=0.
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iux)=ampl_ux(j)*sin(x(l1:l2))*cos(z(n))
            f(l1:l2,m,n,iuz)=-ampl_ux(j)*cos(x(l1:l2))*sin(z(n))
          enddo; enddo
          f(:,:,:,iuz)=0.   !!!
!
        case ('incompressive-shwave')
! incompressible shear wave of Johnson & Gammine (2005a)
          print*, "incomp-shwave: ampl_ux/ky_uu = ", ampl_ux(j)/ky_uu
! Get the streamfunction, save it in the iuz slot
          call sinwave_phase(f,iuz,ampl_ux(j)/ky_uu,&
              kx_uu,ky_uu,kz_uu,phase_ux(j))
! Set the boundaries before taking the curl
          call update_ghosts(f)
! 2D curl
          do n=n1,n2;do m=m1,m2
            call grad(f,iuz,tmp_nx3)
            f(l1:l2,m,n,iux) = -tmp_nx3(:,2)
            f(l1:l2,m,n,iuy) =  tmp_nx3(:,1)
          enddo;enddo
          f(:,:,:,iuz)=0.
!
        case ('cylcoords-stream-x')
! Constant velocity in x-direction in cylinder coordinates
          do l=l1,l2; do m=m1,m2
            wall_smoothing=1-exp(-((x(l)-r_cyl)/skin_depth)**2)
            f(l,m,:,iux)=cos(y(m))*ampluu(j)*wall_smoothing
            f(l,m,:,iuy)=-sin(y(m))*ampluu(j)*wall_smoothing
          enddo; enddo
          f(:,:,:,iuz)=0.
          f(1:l1,:,:,iux:iuz)=0.
!
        case ('cylinderstream_cyl')
!   Stream functions for flow around a cylinder as initial condition.
!   Cylindrical coordinates. Flow in x-direction.
          a2 = r_cyl**2
          do l=l1,l2
            if (x(l) < r_cyl) then
              f(l,:,:,iux:iuz) = 0
            else
              rr2 = x(l)**2
              wall_smoothing=1-exp(-((x(l)-r_cyl)/skin_depth)**2)
              do m=m1,m2
                f(l,m,:,iux) = ampluu(j)*cos(y(m))&
                               *(1. - a2/rr2)*wall_smoothing
                f(l,m,:,iuy) = -ampluu(j)*sin(y(m))&
                               *(1. + a2/rr2)*wall_smoothing
              enddo
            endif
          enddo
!
        case('Gressel-hs')
          call information('init_uu', &
              'Gressel hydrostatic equilibrium setup done in interstellar')
!
        case('spher-harm-poloidal')
!
!  Poloidal flow, defined by spherical harmonic ll_sh(j), mm_sh(j) 
!  and radial profile of the generating scalar tmp, depending on n_xprof(j) (see below)
!
          if (.not.lspherical_coords) call fatal_error("init_uu", &
              "spher-harm-poloidal only meaningful for spherical coordinates")

          if (n_xprof(j)==-1) then
            tmp=(x(l1:l2)-xyz0(1))*(x(l1:l2)-xyz1(1))                 ! parabola, vanishing at top and bottom
            prof=tmp/x(l1:l2) + 2.*x(l1:l2)-(xyz0(1)+xyz1(1))         ! f/r + d f/dr
          else
            tmp=sin((2.*pi/(Lxyz(1))*n_xprof(j))*(x(l1:l2)-xyz0(1)))  ! sine with wavenumber n_xprof(j)
            prof=tmp/x(l1:l2) + (2.*pi/(Lxyz(1))*n_xprof(j))*cos((2.*pi/(Lxyz(1))*n_xprof(j))*(x(l1:l2)-xyz0(1)))
          endif

          if (lyang) then
            allocate(yz(2,ny*nz))
            call yin2yang_coors(costh(m1:m2),sinth(m1:m2),cosph(n1:n2),sinph(n1:n2),yz)
            iyz=1
            do m=m1,m2
              do n=n1,n2
                sph=ampluu(j)*ylm_other(yz(1,iyz),yz(2,iyz),ll_sh(j),mm_sh(j),sph_har_der)
                tmp_nx3(:,1)=2.*tmp*sph
                tmp_nx3(:,2)=ampluu(j)*prof*sph_har_der
                if (mm_sh(j)/=0) then
                  tmp_nx3(:,3) = -sph*prof*mm_sh(j)/sin(yz(1,iyz))*sin(mm_sh(j)*yz(2,iyz))/cos(mm_sh(j)*yz(2,iyz))
                else
                  tmp_nx3(:,3) = 0.
                endif
                call transform_thph_yy( tmp_nx3, (/1,1,1/), f(l1:l2,m,n,iux:iuz), yz(1,iyz), yz(2,iyz) )
                iyz=iyz+1
              enddo
            enddo
          else
            do n=n1,n2
              do m=m1,m2
                sph=ampluu(j)*ylm(ll_sh(j),mm_sh(j),sph_har_der)
                f(l1:l2,m,n,iux) = ll_sh(j)*(ll_sh(j)+1)*tmp*sph
                f(l1:l2,m,n,iuy) = ampluu(j)*prof*sph_har_der
                if (mm_sh(j)/=0) &      ! tb improved!
                  f(l1:l2,m,n,iuz) = -prof*sph*mm_sh(j)/sinth(m)* sin(mm_sh(j)*z(n))/cos(mm_sh(j)*z(n))
              enddo
            enddo
          endif
!
        case('parker_wind')
!
!  Parker wind
!
          if (.not.lspherical_coords) call fatal_error("init_uu", &
              "parker_wind only meaningful for spherical coordinates")
!
!  initial iteration step
!
          where (x(l1:l2) < .5) 
            ur=2.*x(l1:l2)
          elsewhere
            ur=2.
          endwhere
          phi0=-1.5+alog(4.)
!
!  iterate
!
          do iter=1,niter
            where (x(l1:l2) < .5) 
              ur=exp(.5*ur**2-2.*alog(x(l1:l2))-1./x(l1:l2)-phi0)
            elsewhere
              ur=sqrt(2.*(alog(ur*x(l1:l2)**2)+1./x(l1:l2)+phi0))
            endwhere
          enddo
          lnrhor=-alog(4.*ur*x(l1:l2)**2)
!
          do n=n1,n2
            do m=m1,m2
              f(l1:l2,m,n,iux)=ur
              f(l1:l2,m,n,ilnrho)=lnrhor
            enddo
          enddo
!
        case default
          !
          !  Catch unknown values
          !
          call fatal_error("init_uu", &
              "No such value for inituu: "//trim(inituu(j)))
!
        endselect
!
!  End loop over initial conditions
!
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_uu(f)
!
!  This allows an extra random velocity perturbation on
!  top of the initialization so far.
!
      if (urand /= 0) then
        if (lroot) print*, 'init_uu: Adding random uu fluctuations.'
        if (urand > 0) then
          do i=iux,iuz
            do n=1,mz; do m=1,my
              call random_number_wrapper(tmpmx)
              f(:,m,n,i) = f(:,m,n,i) + urand*(tmpmx-0.5)
            enddo; enddo
          enddo
        else
          if (lroot) print*, 'init_uu: Multiplicative fluctuations.'
          do i=iux,iuz
            do n=1,mz; do m=1,my
              call random_number_wrapper(tmpmx)
              f(:,m,n,i) = f(:,m,n,i) * urand*(tmpmx-0.5)
            enddo; enddo
          enddo
        endif
      endif
!
! mgellert, add random fluctuation only inside domain, not on boundary
!           (to be able to use the 'freeze' option for BCs)
      if (urandi /= 0) then
        if (urandi > 0) then
          if (lroot) print*, 'init_uu: Adding random uu fluctuations (not on boundary), urandi=',urandi
          do i=iux,iuz
            do n=n1+1,n2-1; do m=m1,m2; do l=l1+1,l2-1
              call random_number_wrapper(tmp0)
              f(l,m,n,i) = f(l,m,n,i) + urandi*(tmp0-0.5)
            enddo; enddo; enddo
          enddo
        else
          if (lroot) print*, 'init_uu: Multiplicative fluctuations (not on boundary), urandi=',urandi
          do i=iux,iuz
            do n=n1+1,n2-1; do m=m1,m2; do l=l1+1,l2-1
              call random_number_wrapper(tmp0)
              f(l,m,n,i) = f(l,m,n,i) * urandi*(tmp0-0.5)
            enddo; enddo; enddo
          enddo
        endif
      endif
!
    endsubroutine init_uu
!***********************************************************************
    subroutine pencil_criteria_hydro
!
!  All pencils that the Hydro module depends on are specified here.
!
!  20-nov-04/anders: coded
!
      if (lparticles_lyapunov .or. lparticles_caustics .or. lparticles_tetrad) & 
        lpenc_requested(i_uij) = .true.
      if (ladvection_velocity) then
        if (lweno_transport) then
          lpenc_requested(i_uu)=.true.
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_transprho)=.true.
          lpenc_requested(i_transpurho)=.true.
        else
          lpenc_requested(i_ugu)=.true.
        endif
      endif
      if (lprecession) lpenc_requested(i_rr)=.true.
      if (ldt.or.(ekman_friction/=0)) lpenc_requested(i_uu)=.true.
      if (Omega/=0.0) lpenc_requested(i_uu)=.true.
!
!  Damping terms for lcylinder_in_a_box
!
      if (tdamp/=0.or.dampuext/=0.or.dampuint/=0) then
        if (lOmega_int.or.rdampext/=impossible.or.rdampint/=impossible) &
            lpenc_requested(i_r_mn)=.true.
        if (lcylinder_in_a_box) lpenc_requested(i_rcyl_mn)=.true.
      endif
!
!  1/rho needed for correcting the damping term
!
      if (tau_damp_ruxm/=0..or.tau_damp_ruym/=0..or.tau_damp_ruzm/=0.) &
        lpenc_requested(i_rho1)=.true.
!
      if (lisotropic_advection) lpenc_requested(i_u2)=.true.
!
!  The following only play a role if lsphere_in_a_box
!  and there is a buoyancy term. (If only hydro, there is no buoyancy.)
!
      if (lboussinesq.and.ltemperature.and.lsphere_in_a_box) then
        lpenc_requested(i_evr)=.true.
        lpenc_requested(i_r_mn)=.true.
      endif
!
      if (lfargo_advection) then
        lpenc_requested(i_uu_advec)=.true.
        lpenc_requested(i_uuadvec_guu)=.true.
      endif
!
!  Request unit vectors for transformation of velocity from Cartesian
!  to spherical coordinates.
!
      if (luu_sph_as_aux.and.lsphere_in_a_box) then
        lpenc_requested(i_evr)=.true.
        lpenc_requested(i_evth)=.true.
        lpenc_requested(i_phix)=.true.
        lpenc_requested(i_phiy)=.true.
      endif
!
!  request rho1 if relativistic
!
      if (lrelativistic) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_ss_rel)=.true.
        lpenc_requested(i_ss_rel2)=.true.
        lpenc_requested(i_ss_rel_ij)=.true.
   !    lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_grho)=.true.
      endif
!
!  video pencils
!
      if (lwrite_slices) then
        if (ivid_oo  /=0) lpenc_video(i_oo)=.true.
        if (ivid_ou  /=0) lpenc_video(i_ou)=.true.
        if (ivid_o2  /=0) lpenc_video(i_o2)=.true.
        if (ivid_divu/=0) lpenc_video(i_divu)=.true.
        if (ivid_u2  /=0) lpenc_video(i_u2)=.true.
        if (ivid_Ma2 /=0) lpenc_video(i_Ma2)=.true.
        if (ivid_uu_sph/=0) lpenc_video(i_uu_sph)=.true.
      endif
!
!  diagnostic pencils
!
      lpenc_diagnos(i_uu)=.true.
      if (idiag_oumphi/=0 .or. idiag_oumxy/=0 .or. &
          idiag_oumxz/=0) lpenc_diagnos2d(i_ou)=.true.
      if (idiag_ozmphi/=0) lpenc_diagnos2d(i_oo)=.true.
      if (idiag_u2mphi/=0) lpenc_diagnos2d(i_u2)=.true.
      if (idiag_ox2m/=0 .or. idiag_oy2m/=0 .or. idiag_oz2m/=0 .or. &
          idiag_ox3m/=0 .or. idiag_oy3m/=0 .or. idiag_oz3m/=0 .or. &
          idiag_ox4m/=0 .or. idiag_oy4m/=0 .or. idiag_oz4m/=0 .or. &
          idiag_oxm /=0 .or. idiag_oym /=0 .or. idiag_ozm /=0 .or. &
          idiag_oxoym/=0 .or. idiag_oxozm/=0 .or. idiag_oyozm/=0 .or. &
          idiag_oxuzxm/=0 .or. idiag_oyuzym/=0 .or. &
          idiag_pvzm /=0 .or. idiag_quxom/=0) &
          lpenc_diagnos(i_oo)=.true.
      if (idiag_orms/=0 .or. idiag_omax/=0 .or. idiag_o2m/=0 .or. idiag_o2u2m/=0 .or. idiag_o2sphm/=0 .or. &
          idiag_ormsh/=0 .or. idiag_o2mz/=0 ) lpenc_diagnos(i_o2)=.true.
      if (idiag_q2m/=0 .or. idiag_qrms/=0 .or. idiag_qmax/=0 .or. &
          idiag_qfm/=0 .or. idiag_qom/=0 .or. idiag_quxom/=0) &
          lpenc_diagnos(i_curlo)=.true.
      if (idiag_divu2m/=0 .or. idiag_divu2mz/=0 .or. idiag_divru2mz/=0 .or. &
          idiag_divum/=0 .or. idiag_rdivum/=0) lpenc_diagnos(i_divu)=.true.
      if (idiag_rdivum/=0 .or. idiag_fkinzm/=0 .or. idiag_ruzupmz/=0 .or. idiag_ruzdownmz/=0) &
          lpenc_diagnos(i_rho)=.true.
      if (idiag_ruxupmxy/=0 .or. idiag_ruxdownmxy/=0) &
          lpenc_diagnos2d(i_rho)=.true.
      if (idiag_gdivu2m/=0) lpenc_diagnos(i_graddivu)=.true.
      if (idiag_oum/=0 .or. idiag_ourms/=0 .or. idiag_oumx/=0 .or. idiag_oumy/=0 .or. &
           idiag_oumz/=0 .or. idiag_oumh/=0 .or. idiag_ou_int/=0 ) lpenc_diagnos(i_ou)=.true.
      if (idiag_oxum/=0) lpenc_diagnos(i_oxu)=.true.
      if (idiag_oxurms/=0) lpenc_diagnos(i_oxu2)=.true.
      if (idiag_divrhourms/=0.or.idiag_divrhoumax/=0) then
        lpenc_diagnos(i_ugrho)=.true.
        lpenc_diagnos(i_divu)=.true.
        lpenc_diagnos(i_rho)=.true.
      endif
      if (idiag_Marms/=0 .or. idiag_Mamax/=0) lpenc_diagnos(i_Ma2)=.true.
      if (idiag_u3u21m/=0 .or. idiag_u3u21mz/=0) lpenc_diagnos(i_u3u21)=.true.
      if (idiag_u1u32m/=0 .or. idiag_u1u32mz/=0) lpenc_diagnos(i_u1u32)=.true.
      if (idiag_u2u13m/=0 .or. idiag_u2u13mz/=0) lpenc_diagnos(i_u2u13)=.true.
      if (idiag_u2u31m/=0 .or. idiag_u2u31mz/=0) lpenc_diagnos(i_u2u31)=.true.
      if (idiag_u3u12m/=0 .or. idiag_u3u12mz/=0) lpenc_diagnos(i_u3u12)=.true.
      if (idiag_u1u23m/=0 .or. idiag_u1u23mz/=0) lpenc_diagnos(i_u1u23)=.true.
      if (idiag_acczmz/=0 .or. idiag_acczupmz/=0 .or. idiag_acczdownmz/=0 .or. &
          idiag_accpowzmz/=0 .or. idiag_accpowzupmz/=0 .or. idiag_accpowzdownmz/=0 .or. &
          idiag_totalforcezmz/=0 .or. idiag_totalforcezupmz/=0 .or. &
          idiag_totalforcezdownmz/=0) then
          lpenc_diagnos(i_fpres)=.true.; lpenc_diagnos(i_fvisc)=.true.
          if (lgrav) lpenc_diagnos(i_gg)=.true.
      endif
      if (idiag_urms/=0 .or. idiag_durms/=0 .or. &
          idiag_umax/=0 .or. idiag_rumax/=0 .or. &
          idiag_fkinzm/=0 .or. idiag_u2m/=0 .or. idiag_um2/=0 .or. idiag_u2mz/=0 .or. &
          idiag_urmsh/=0 .or. idiag_urmsx/=0 .or. idiag_urmsz/=0 .or. idiag_u2sphm/=0) &
          lpenc_diagnos(i_u2)=.true.
      if (idiag_u2sphm/=0 .or. idiag_o2sphm/=0) lpenc_diagnos(i_r_mn)=.true.
      if (idiag_duxdzma/=0 .or. idiag_duydzma/=0 .or. lgradu_as_aux) lpenc_diagnos(i_uij)=.true.
      if (idiag_fmasszmz/=0 .or. idiag_ruxuym/=0 .or. &
          idiag_ruxm/=0 .or. idiag_ruym/=0 .or. idiag_ruzm/=0 .or. &
          idiag_ruxuzm/=0 .or. idiag_ruyuzm/=0 .or. idiag_pvzm/=0 .or. &
          idiag_ruxtot/=0) lpenc_diagnos(i_rho)=.true.
      if (idiag_ruxmz/=0 .or. idiag_ruymz/=0 .or. idiag_ruzmz/=0 .or. &
          idiag_rux2mz/=0 .or. idiag_ruy2mz/=0 .or. idiag_ruz2mz/=0 .or. &
          idiag_ruxmx/=0 .or. idiag_ruymx/=0 .or. idiag_ruzmx/=0 .or. &
          idiag_ruxuymz/=0 .or. idiag_ruxuzmz/=0 .or. idiag_ruyuzmz/=0 .or. &
          idiag_ruxuy2mz/=0 .or. idiag_ruxuz2mz/=0 .or. idiag_ruyuz2mz/=0 .or. idiag_uduum/=0) &
          lpenc_diagnos(i_rho)=.true.
      if (idiag_rux2mx /= 0 .or. idiag_ruy2mx /= 0 .or. idiag_ruz2mx /= 0 .or. &
          idiag_ruxuymx /= 0 .or. idiag_ruxuzmx /= 0 .or. idiag_ruyuzmx /= 0 .or. &
          idiag_rurmphi/=0 .or. idiag_rupmphi/=0 .or. idiag_ruzmphi/=0 .or. &
          idiag_rurupmphi/=0 .or. idiag_ruruzmphi/=0 .or. idiag_rursphmphi/=0 .or. &
          idiag_ruthmphi/=0) lpenc_diagnos(i_rho) = .true.
      if (idiag_ormr/=0 .or. idiag_opmr/=0 .or. idiag_ozmr/=0) &
          lpenc_diagnos(i_oo)=.true.
      if (idiag_oxmxy/=0 .or. idiag_oymxy/=0 .or. idiag_ozmxy/=0 .or. &
          idiag_oxmz/=0 .or. idiag_oymz/=0 .or. idiag_ozmz/=0 .or. &
          idiag_ox2mz/=0 .or. idiag_oy2mz/=0 .or. idiag_oz2mz/=0 .or. &
          idiag_pvzmxy/=0) &
          lpenc_diagnos2d(i_oo)=.true.
      if (idiag_pvzmxy/=0) lpenc_diagnos2d(i_rho)=.true.
      if (idiag_totangmom/=0 ) lpenc_diagnos(i_rcyl_mn)=.true.
      if (idiag_urmr/=0 .or.  idiag_ormr/=0 .or. idiag_urmphi/=0 .or. idiag_rurmphi/=0 .or. &
          idiag_ur2mphi/=0 .or. idiag_urupmphi/=0 .or. idiag_uruzmphi/=0 .or. &
          idiag_rurupmphi/=0 .or. idiag_ruruzmphi/=0) then
        lpenc_diagnos(i_pomx)=.true.
        lpenc_diagnos(i_pomy)=.true.
      endif
      if (idiag_ursphmphi/=0 .or. idiag_rursphmphi/=0 .or. idiag_fkinrsphmphi/=0) &
        lpenc_diagnos2d(i_evr)=.true.
      if (idiag_uthmphi/=0 .or. idiag_ruthmphi/=0) lpenc_diagnos2d(i_evth)=.true.
      if (idiag_upmr/=0 .or. idiag_opmr/=0 .or. idiag_upmphi/=0 .or. idiag_rupmphi/=0 .or. &
          idiag_up2mphi/=0 .or. idiag_urupmphi/=0 .or. idiag_upuzmphi/=0 .or. &
          idiag_rurupmphi/=0 .or. idiag_rupuzmphi/=0) then
        lpenc_diagnos(i_phix)=.true.
        lpenc_diagnos(i_phiy)=.true.
      endif
      if (idiag_EEK/=0 .or. idiag_ekin/=0 .or. idiag_ekintot/=0 .or. idiag_fkinzmz/=0 .or. &
           idiag_fkinzupmz/=0 .or. idiag_fkinzdownmz/=0 .or. &
           idiag_ekinmx /= 0 .or. idiag_ekinmz/=0 .or. idiag_fkinxmx/=0) then
        lpenc_diagnos(i_ekin)=.true.
      endif
      if (idiag_fkinxmxy/=0 .or. idiag_fkinymxy/=0 .or. &
          idiag_fkinxupmxy/=0 .or. idiag_fkinxdownmxy/=0 .or. idiag_fkinrsphmphi/=0) then
        lpenc_diagnos2d(i_uu)=.true.
        lpenc_diagnos2d(i_ekin)=.true.
      endif
      if (idiag_uxmxy/=0 .or. idiag_uymxy/=0 .or. idiag_uzmxy/=0 .or. &
          idiag_ux2mxy/=0 .or. idiag_uy2mxy/=0 .or. idiag_uz2mxy/=0 .or. &
          idiag_ruxmxy/=0 .or. idiag_ruymxy/=0 .or. idiag_ruzmxy/=0 .or.  &
          idiag_uxuymxy/=0 .or. idiag_uxuzmxy/=0 .or. idiag_uyuzmxy/=0 .or. &
          idiag_ruxuymxy/=0 .or. idiag_ruxuzmxy/=0 .or. &
          idiag_uxupmxy/=0 .or. idiag_uxdownmxy/=0 .or. &
          idiag_ux2upmxy/=0 .or. idiag_ux2downmxy/=0 .or. &
          idiag_ruxupmxy/=0 .or. idiag_ruxdownmxy/=0 .or. &
          idiag_ruyuzmxy/=0 .or. idiag_ffdownmxy/=0) then
        lpenc_diagnos2d(i_uu)=.true.
      endif
      if (idiag_ruxmxy/=0 .or. idiag_ruymxy/=0 .or. idiag_ruzmxy/=0 .or. &
          idiag_ruxuymxy/=0 .or. idiag_ruxuzmxy/=0 .or. idiag_ruyuzmxy/=0) &
        lpenc_diagnos2d(i_rho)=.true.
      if (idiag_uguxm/=0 .or. idiag_uguym/=0 .or. idiag_uguzm/=0) &
          lpenc_diagnos(i_ugu)=.true.
      if (idiag_uguxmxy/=0 .or. idiag_uguymxy/=0 .or. idiag_uguzmxy/=0) &
          lpenc_diagnos2d(i_ugu)=.true.
!
      if (idiag_ugu2m/=0 .or. idiag_ugurmsx/=0) lpenc_diagnos(i_ugu2)=.true.
      if (idiag_uguxmx/=0 .or. idiag_uguymx/=0 .or. idiag_uguzmx/=0 .or. &
          idiag_uguxmy/=0 .or. idiag_uguymy/=0 .or. idiag_uguzmy/=0 .or. &
          idiag_uguxmz/=0 .or. idiag_uguymz/=0 .or. idiag_uguzmz/=0) &
          lpenc_diagnos(i_ugu)=.true.
!
      if (idiag_Remz/=0) then
        lpenc_diagnos(i_ugu2)=.true.
        lpenc_diagnos(i_diffus_total)=.true.
      endif
!
      if (idiag_oguxmz/=0 .or. idiag_oguymz/=0 .or. idiag_oguzmz/=0 .or. &
          idiag_ogux2mz/=0 .or. idiag_oguy2mz/=0 .or. idiag_oguz2mz/=0) &
          lpenc_diagnos(i_ogu)=.true.
!
      if (idiag_oxdivumz/=0  .or. idiag_oydivumz/=0  .or. idiag_ozdivumz/=0 .or. &
          idiag_oxdivu2mz/=0 .or. idiag_oydivu2mz/=0 .or. idiag_ozdivu2mz/=0) then
        lpenc_diagnos(i_oo)=.true.
        lpenc_diagnos(i_divu)=.true.
      endif
!
      if (idiag_udpxxm/=0 .or. idiag_udpyym/=0 .or. idiag_udpzzm/=0 .or. &
          idiag_udpxym/=0 .or. idiag_udpyzm/=0 .or. idiag_udpxzm/=0) then
        lpenc_diagnos(i_uu)=.true.
        lpenc_diagnos(i_fpres)=.true.
      endif
!
! check whether right variables are set for half-box calculations.
!
      if (idiag_urmsn/=0 .or. idiag_ormsn/=0 .or. idiag_oumn/=0) then
        if ((.not.lequatory).and.(.not.lequatorz)) then
          call fatal_error("pencil_criteria_hydro","You have to set either of"// &
              "lequatory or lequatorz to true to calculate averages over half the box")
        else
          if (lequatory) write(*,*) 'pencil-criteria_hydro: box divided along y dirn'
          if (lequatorz) write(*,*) 'pencil-criteria_hydro: box divided along z dirn'
        endif
      else
      endif
!
    endsubroutine pencil_criteria_hydro
!***********************************************************************
    subroutine pencil_interdep_hydro(lpencil_in)
!
!  Interdependency among pencils from the Hydro module is specified here.
!
!  20-nov-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      if (lpencil_in(i_u2)) lpencil_in(i_uu)=.true.
      if (lpencil_in(i_divu)) lpencil_in(i_uij)=.true.
      if (lpencil_in(i_sij)) then
        lpencil_in(i_uij)=.true.
        lpencil_in(i_divu)=.true.
      endif
      !if (.not.lcartesian_coords.or.lalways_use_gij_etc.or.lpencil(i_oij)) then
      if (.not.lcartesian_coords.or.lalways_use_gij_etc) then
        if (lpencil_in(i_del2u)) then
          lpencil_in(i_graddivu)=.true.
          lpencil_in(i_curlo)=.true.
        endif
      endif
      if (lpencil_in(i_oo) .and. ioo == 0) lpencil_in(i_uij) = .true.
      if (lpencil_in(i_o2)) lpencil_in(i_oo)=.true.
      if (lpencil_in(i_oxu2)) lpencil_in(i_oxu)=.true.
      if (lpencil_in(i_ou) .or. lpencil_in(i_oxu)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_oo)=.true.
      endif
      if (lpencil_in(i_ugu)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_uij)=.true.
      endif
      if (lpencil_in(i_ogu)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_oo)=.true.
        lpencil_in(i_uij)=.true.
      endif
      if (lpencil_in(i_u3u21) .or. &
          lpencil_in(i_u1u32) .or. &
          lpencil_in(i_u2u13) .or. &
          lpencil_in(i_u2u31) .or. &
          lpencil_in(i_u3u12) .or. &
          lpencil_in(i_u1u23)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_uij)=.true.
      endif
!
!  other general interdependencies
!
      if (lpencil_in(i_sij2)) lpencil_in(i_sij)=.true.
!     if (lpencil_in(i_del2u)) lpencil_in(i_curlo)=.true.
      if (lpencil_in(i_curlo)) lpencil_in(i_oij)=.true.
      if (lpencil_in(i_oij).and.(.not.lcartesian_coords)) &
          lpencil_in(i_oo)=.true.
!
      if (lpencil_in(i_uu_advec)) lpencil_in(i_uu)=.true.
      if (lpencil_in(i_uuadvec_guu)) then
        lpencil_in(i_uu_advec)=.true.
        lpencil_in(i_uij)=.true.
      endif
!
    endsubroutine pencil_interdep_hydro
!***********************************************************************
    subroutine calc_pencils_hydro_pencpar(f,p,lpenc_loc)
!
! Calls linearized or nonlinear hydro routine depending on whether llinearized_hydro is
! selected or not (the default is nonlinear).
!
! 24-jun-13/dhruba: coded
! 20-sep-13/MR    : added optional list of indices in lpencil, penc_inds,
!                   for which pencils are calculated, default: all
! 21-sep-13/MR    : returned to pencil mask as parameter lpenc_loc
!
      real, dimension (mx,my,mz,mfarray),intent(INOUT):: f
      type (pencil_case),                intent(INOUT):: p
      logical, dimension(:),             intent(IN) :: lpenc_loc
!
      if (llinearized_hydro) then
        call calc_pencils_hydro_linearized(f,p,lpenc_loc)
      else
        call calc_pencils_hydro_nonlinear(f,p,lpenc_loc)
      endif
!
      endsubroutine calc_pencils_hydro_pencpar
!***********************************************************************
    subroutine calc_pencils_hydro_std(f,p)
!
! Envelope adjusting calc_pencils_hydro_pencpar to the standard use with
! lpenc_loc=lpencil
!
! 21-sep-13/MR    : coded
!
      real, dimension (mx,my,mz,mfarray),intent(INOUT):: f
      type (pencil_case),                intent(INOUT):: p
!
      call calc_pencils_hydro_pencpar(f,p,lpencil)
!
      endsubroutine calc_pencils_hydro_std
!***********************************************************************
    subroutine calc_pencils_hydro_nonlinear(f,p,lpenc_loc)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  08-nov-04/tony: coded
!  26-mar-07/axel: started using the gij_etc routine
!  24-jun-13/dhruba: the earlier calc_pencils_hydro routine is copied here.
!  20-sep-13/MR    : added optional list of indices in lpencil, penc_inds,
!                    for which pencils are calculated, default: all
!  21-sep-13/MR    : instead pencil mask as parameter lpenc_loc
!  19-feb-15/MR    : adapted weno transport to use of reference state
!
      use Deriv
      use Sub
      use WENO_transport
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      real, dimension (nx) :: tmp, c_sld_im12, c_sld_ip12, sqrt_ss_term, ratio_mom2en
      real, dimension (nx,3) :: tmp3
      integer :: i, j, ju, jj, kk, jk
!
      intent(in)   :: lpenc_loc
      intent(inout):: f,p
!
! uu, is either directly p%uu=f(l1:l2,m,n,iux:iuz), or, if lrelativistic
! p%ss_rel=f(l1:l2,m,n,iux:iuz), and we compute p%lorentz_gamma and
! p%lorentz_gamma2, so uu = ss_rel_factor * ss_rel
! from Eq.(11) of Brandenburg, Enqvist, & Olesen (1996).
!
      if (lpenc_loc(i_uu)) then
        if (lrelativistic) then
          p%ss_rel=f(l1:l2,m,n,iux:iuz)
          call dot2_mn(p%ss_rel,p%ss_rel2)
          !sqrt_ss_term=sqrt(1.+2.25*p%ss_rel2*p%rho1**2)
          !p%lorentz_gamma2=.5*(1.+sqrt_ss_term)
          ratio_mom2en=p%ss_rel2/p%totenergy_rel
          p%lorentz_gamma2=.5/(1.-ratio_mom2en)*(1.-.5*ratio_mom2en+sqrt(1.-.75*ratio_mom2en))
          p%ss_rel_factor=.75*p%rho1/p%lorentz_gamma2
          call multsv_mn(p%ss_rel_factor,p%ss_rel,p%uu)
        else
          p%uu=f(l1:l2,m,n,iux:iuz)
        endif
      endif
! u2
      if (lpenc_loc(i_u2)) call dot2_mn(p%uu,p%u2)
!
! uij = u_i,j = ss_rel_factor * ss_rel_i,j + ss_rel_factor_{,j} * ss_rel_i
! 
      if (lpenc_loc(i_uij)) then
        if (lrelativistic) then
          if (.not.lpenc_loc(i_uu)) &
            call fatal_error('calc_pencils_hydro_nonlinear','need lpenc_loc(i_uu)')
          call gij(f,iuu,p%ss_rel_ij,1)
          call multsv_mn(p%ss_rel_factor,p%ss_rel_ij,p%uij)
!
!  add ss_rel_factor_{,j} = ss_rel_factor*(-lnrho_{,j}-lngam_{,j}) * ss_rel_i
!
          call multmv_transp(p%ss_rel_ij,p%ss_rel,tmp3)
          call multsv_add(2.*tmp3,-p%ss_rel2,p%glnrho,tmp3)
          call multsv_mn(.5/sqrt_ss_term*(9./16.)*p%rho1**2/p%lorentz_gamma2,tmp3,tmp3)
          call multvv_smat_add(p%ss_rel_factor,-p%glnrho-tmp3,p%ss_rel,p%uij)
        else
          call gij(f,iuu,p%uij,1)
        endif
!
!  if gradu is to be stored as auxiliary then we store it now
!
        if (lgradu_as_aux .or. lparticles_lyapunov .or. lparticles_caustics .or. lparticles_tetrad) then
          jk=0
          do jj=1,3; do kk=1,3
            f(l1:l2,m,n,iguij+jk) = p%uij(:,jj,kk)
            jk=jk+1
          enddo;enddo
        endif
      endif
!      if (.not.lpenc_loc_check_at_work) then
!        write(*,*) 'uurad,rad',p%uij(1:6,1,1)
!      endif
! divS
      if (lpenc_loc(i_divss_rel)) then
        call div_mn(p%ss_rel_ij,p%divss_rel,p%ss_rel)
      endif
! divu
      if (lpenc_loc(i_divu)) then
        call div_mn(p%uij,p%divu,p%uu)
      endif
! sij
      if (lpenc_loc(i_sij)) call traceless_strain(p%uij,p%divu,p%sij,p%uu,lshear_rateofstrain)
! sij2
      if (lpenc_loc(i_sij2)) call multm2_sym_mn(p%sij,p%sij2)
! uij5
      if (lpenc_loc(i_uij5)) call gij(f,iuu,p%uij5,5)
! oo (=curlu)
      if (lpenc_loc(i_oo)) then
        if (ioo /= 0) then
          p%oo = f(l1:l2,m,n,iox:ioz)
        else
          call curl_mn(p%uij,p%oo,p%uu)
        endif
      endif
! o2 and oxu2
      if (lpenc_loc(i_o2)) call dot2_mn(p%oo,p%o2)
! ou and oxu
      if (lpenc_loc(i_ou)) call dot_mn(p%oo,p%uu,p%ou)
      if (lpenc_loc(i_oxu)) call cross(p%oo,p%uu,p%oxu)
      if (lpenc_loc(i_oxu2)) call dot2_mn(p%oxu,p%oxu2)
! Useful to debug forcing - Dhruba
      if (loutest.and.lpenc_loc(i_ou))then
!      write(*,*) lpenc_loc(i_ou)
        outest = minval(p%ou)
        if (outest<(-1.0d-8))then
          write(*,*) m,n,outest,maxval(p%ou),lpenc_loc(i_ou)
          write(*,*)'WARNING : hydro:ou has different sign than relhel'
        endif
      endif
!
! ugu. In the relativistic case, we compute u.gradS, so we need to input p%ss_rel_ij,
! but p%uu (as opposed to p%ss_rel) is neededm because it enteres as prefactor for upwinding.
!
      if (lpenc_loc(i_ugu)) then
        if (headtt.and.lupw_uu) print *,'calc_pencils_hydro: upwinding advection term'
        if (lrelativistic) then
          call u_dot_grad(f,iuu,p%ss_rel_ij,p%uu,p%ugu,UPWIND=lupw_uu)
        else
          call u_dot_grad(f,iuu,p%uij,p%uu,p%ugu,UPWIND=lupw_uu)
        endif
!
!      if (.not.lpenc_loc_check_at_work) then
!        write(*,*) 'ugu',p%ugu(1:6,1)
!      endif
!        if (.not.lpenc_loc_check_at_work) then
!          write(*,*) 'DM',x(l1:l2)
!          write(*,*) 'DM',p%uu(:,1)
!          write(*,*) 'DM',p%uij(:,1,1)
!        endif
!
!  If lffree switched is used, we need to turn off the u.gradu term
!  to ensure momentum conservation.
!
        if (ldensity) then
          if (lffree) then
            tmp=profx_ffree*profy_ffree(m)*profz_ffree(n)
            do j=1,3
              p%ugu(:,j)=p%ugu(:,j)*tmp
            enddo
          endif
        endif
      endif
! ugu2
      if (lpenc_loc(i_ugu2)) call dot2_mn(p%ugu,p%ugu2)
! ogu ... ogu2
      if (lpenc_loc(i_ogu)) call u_dot_grad(f,iuu,p%uij,p%oo,p%ogu,UPWIND=lupw_uu)
! u3u21, u1u32, u2u13, u2u31, u3u12, u1u23
      if (lpenc_loc(i_u3u21)) p%u3u21=p%uu(:,3)*p%uij(:,2,1)
      if (lpenc_loc(i_u1u32)) p%u1u32=p%uu(:,1)*p%uij(:,3,2)
      if (lpenc_loc(i_u2u13)) p%u2u13=p%uu(:,2)*p%uij(:,1,3)
      if (lpenc_loc(i_u2u31)) p%u2u31=p%uu(:,2)*p%uij(:,3,1)
      if (lpenc_loc(i_u3u12)) p%u3u12=p%uu(:,3)*p%uij(:,1,2)
      if (lpenc_loc(i_u1u23)) p%u1u23=p%uu(:,1)*p%uij(:,2,3)
! del4u, del6u, del4graddivu, and del6u_strict
      if (lpenc_loc(i_del4u)) call del4v(f,iuu,p%del4u)
      if (lpenc_loc(i_del6u)) call del6v(f,iuu,p%del6u)
      if (lpenc_loc(i_del6u_strict)) call del6v(f,iuu,p%del6u_strict,LSTRICT=.true.)
      if (lpenc_loc(i_del4graddivu)) call del4graddiv(f,iuu,p%del4graddivu)
! del6u_bulk
      if (lpenc_loc(i_del6u_bulk)) then
        call der6(f,iux,p%del6u_bulk(:,1),1)
        call der6(f,iuy,p%del6u_bulk(:,2),2)
        call der6(f,iuz,p%del6u_bulk(:,3),3)
      endif
! der6u_res
      if (lpenc_loc(i_der6u_res)) then
        if (lcartesian_coords) call fatal_error("calc_pencils_hydro_nonlinear",&
                  "der6u_res: This pencil is not implemented for Cartesian coordinates")
        do j=1,3
          ju=j+iuu-1
          do i=1,3
            if (lcylindrical_coords.and.ju==iuy.and.i==1) then  
              call der6(i,f(:,m,n,iuy)-uu_average_cyl(:,n),p%der6u_res(:,i,j),IGNOREDX=.true.)                
            elseif (lspherical_coords.and.ju==iuz.and.i==1) then  
              call der6(i,f(:,m,n,iuz)-uu_average_sph(:,m),p%der6u_res(:,i,j),IGNOREDX=.true.)                
            else   
              call der6(f,ju,p%der6u_res(:,i,j),i,IGNOREDX=.true.)
            endif
          enddo
        enddo
      endif
!
! del2u, graddivu
!
      if (.not.lcartesian_coords.or.lalways_use_gij_etc) then
        if (lpenc_loc(i_graddivu)) then
          if (headtt.or.ldebug) print*,'calc_pencils_hydro: call gij_etc'
          call gij_etc(f,iuu,p%uu,p%uij,p%oij,GRADDIV=p%graddivu)
        endif
        if (lpenc_loc(i_del2u)) then
          call curl_mn(p%oij,p%curlo,p%oo)
          p%del2u=p%graddivu-p%curlo
        endif
      else
!
!  all 3 together
!
        if (lpenc_loc(i_del2u).and.lpenc_loc(i_graddivu).and.lpenc_loc(i_curlo)) then
          call del2v_etc(f,iuu,DEL2=p%del2u,GRADDIV=p%graddivu,CURLCURL=p%curlo)
!
!  all 3 possible pairs
!
        elseif (lpenc_loc(i_del2u).and.lpenc_loc(i_graddivu)) then
          call del2v_etc(f,iuu,DEL2=p%del2u,GRADDIV=p%graddivu)
        elseif (lpenc_loc(i_del2u).and.lpenc_loc(i_curlo)) then
          call del2v_etc(f,iuu,DEL2=p%del2u,CURLCURL=p%curlo)
        elseif (lpenc_loc(i_graddivu).and.lpenc_loc(i_curlo)) then
          call del2v_etc(f,iuu,GRADDIV=p%graddivu,CURLCURL=p%curlo)
!
!  all 3 individually
!
        elseif (lpenc_loc(i_del2u)) then
          call del2v_etc(f,iuu,DEL2=p%del2u)
        elseif (lpenc_loc(i_graddivu)) then
          call del2v_etc(f,iuu,GRADDIV=p%graddivu)
        elseif (lpenc_loc(i_curlo)) then
          call del2v_etc(f,iuu,CURLCURL=p%curlo)
        endif
      endif
!
!del2uj, d^u/dx^2 etc
!
      if (lpenc_loc(i_d2uidxj)) call d2fi_dxj(f,iuu,p%d2uidxj)
!
! deluidxjk
!
      if (lpenc_loc(i_uijk)) call del2fi_dxjk(f,iuu,p%uijk)
!
! grad5divu
!
      if (lpenc_loc(i_grad5divu)) then
        do i=1,3
          p%grad5divu(:,i) = 0.0
          do j=1,3
            ju=iuu+j-1
            call der5i1j(f,ju,tmp,i,j)
            p%grad5divu(:,i) = p%grad5divu(:,i) + tmp
          enddo
        enddo
      endif
! transpurho
      if (lpenc_loc(i_transpurho)) then
        if (lreference_state) then
          call weno_transp(f,m,n,iux,irho,iux,iuy,iuz,p%transpurho(:,1),dx_1,dy_1,dz_1, &
                           ref1=reference_state(:,iref_rho))
          call weno_transp(f,m,n,iuy,irho,iux,iuy,iuz,p%transpurho(:,2),dx_1,dy_1,dz_1, &
                           ref1=reference_state(:,iref_rho))
          call weno_transp(f,m,n,iuz,irho,iux,iuy,iuz,p%transpurho(:,3),dx_1,dy_1,dz_1, &
                           ref1=reference_state(:,iref_rho))
        else
          call weno_transp(f,m,n,iux,irho,iux,iuy,iuz,p%transpurho(:,1),dx_1,dy_1,dz_1)
          call weno_transp(f,m,n,iuy,irho,iux,iuy,iuz,p%transpurho(:,2),dx_1,dy_1,dz_1)
          call weno_transp(f,m,n,iuz,irho,iux,iuy,iuz,p%transpurho(:,3),dx_1,dy_1,dz_1)
        endif
      endif
!
      if (lpenc_loc(i_uu_advec)) then
!
! Advect by the original radial and vertical, but residual azimuthal (= relative) velocity
!
        p%uu_advec(:,1)=p%uu(:,1)
        if (lcylindrical_coords) then
          p%uu_advec(:,2)=p%uu(:,2)-uu_average_cyl(l1:l2,n)
          p%uu_advec(:,3)=p%uu(:,3)
        elseif (lspherical_coords) then
          p%uu_advec(:,2)=p%uu(:,2)
          p%uu_advec(:,3)=p%uu(:,3)-uu_average_sph(l1:l2,m)
        endif
!
        do j=1,3
          ! This is calling scalar h_dot_grad, that does not add
          ! the inertial terms. They will be added here.
          tmp3 = p%uij(:,j,:)
          call h_dot_grad(p%uu_advec,tmp3,p%uuadvec_guu(:,j))
        enddo
        if (lcylindrical_coords) then
          p%uuadvec_guu(:,1)=p%uuadvec_guu(:,1)-rcyl_mn1*p%uu(:,2)*p%uu(:,2)
          p%uuadvec_guu(:,2)=p%uuadvec_guu(:,2)+rcyl_mn1*p%uu(:,2)*p%uu(:,1)
        elseif (lspherical_coords) then
          p%uuadvec_guu(:,1)=p%uuadvec_guu(:,1)-r1_mn*(p%uu(:,2)*p%uu(:,2)+p%uu(:,3)*p%uu(:,3))
          p%uuadvec_guu(:,2)=p%uuadvec_guu(:,2)+r1_mn*(p%uu(:,2)*p%uu(:,1)-p%uu(:,3)*p%uu(:,3)*cotth(m))
          p%uuadvec_guu(:,3)=p%uuadvec_guu(:,3)+r1_mn*(p%uu(:,3)*p%uu(:,1)+p%uu(:,3)*p%uu(:,2)*cotth(m))
        endif
      endif
!
      if (lpenc_loc(i_uu_sph).and.luu_sph_as_aux) p%uu_sph=f(l1:l2,m,n,iuu_sphr:iuu_sphp)
!
    endsubroutine calc_pencils_hydro_nonlinear
!***********************************************************************
    subroutine calc_pencils_hydro_linearized(f,p,lpenc_loc)
!
!  Calculate linearized hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-jun-13/dhruba: aped from calc_pencils_hydro from previous version
!  20-sep-13/MR    : added optional list of indices in lpencil, penc_inds,
!                    for which pencils are calculated, default: all
!  21-sep-13/MR    : instead pencil mask as parameter lpenc_loc
!
      use Deriv
      use Sub
      use WENO_transport
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      real, dimension (nx) :: tmp, tmp2
      real, dimension (nx,3) :: ugu0,u0gu
      integer :: i, j, ju
!
      intent(in) :: f, lpenc_loc
      intent(out) :: p
! uu
      if (lpenc_loc(i_uu)) p%uu=f(l1:l2,m,n,iux:iuz)
      if (lpenc_loc(i_uu0)) p%uu0=f(l1:l2,m,n,iu0x:iu0z)
! u2, should not be calculated
      if (lpenc_loc(i_u2)) call fatal_error('calc_pencils_hydro_linearized:','does not calculate u2 pencil')
! uij
      if (lpenc_loc(i_uij)) call gij(f,iuu,p%uij,1)
      if (lpenc_loc(i_u0ij)) call gij(f,iuu0,p%u0ij,1)
! divu
      if (lpenc_loc(i_divu)) call div_mn(p%uij,p%divu,p%uu)
      if (lpenc_loc(i_divu0)) call div_mn(p%u0ij,p%divu0,p%uu0)
! sij
      if (lpenc_loc(i_sij)) call traceless_strain(p%uij,p%divu,p%sij,p%uu,lshear_rateofstrain)
! sij2
      if (lpenc_loc(i_sij2)) call fatal_error('calc_pencils_hydro_linearized:','does not calculate sij2 pencil')
! uij5
      if (lpenc_loc(i_uij5)) call fatal_error('calc_pencils_hydro_linearized:','does not calculate uij5 pencil')
! oo (=curlu)
      if (lpenc_loc(i_oo)) then
        call curl_mn(p%uij,p%oo,p%uu)
      endif
! o2
      if (lpenc_loc(i_o2) .or. lpenc_loc(i_oxu2)) &
        call fatal_error('calc_pencils_hydro_linearized:','does not calculate o2 or oxu2 pencils')
! ou
      if (lpenc_loc(i_ou) .or. lpenc_loc(i_oxu)) &
        call fatal_error('calc_pencils_hydro_linearized:','does not calculate ou or oxu pencils')
! ugu
      if (lpenc_loc(i_ugu)) then
        if (headtt.and.lupw_uu) print *,'calc_pencils_hydro: upwinding advection term'
!        call u_dot_grad(f,iuu,p%uij,p%uu,p%ugu,UPWIND=lupw_uu)
        if (lrelativistic) then
          call fatal_error('calc_pencils_hydro_linearized:','does not calculate ugu pencil')
        else
          call u_dot_grad(f,iuu0,p%u0ij,p%uu,ugu0,UPWIND=lupw_uu)
          call u_dot_grad(f,iuu,p%uij,p%uu0,u0gu,UPWIND=lupw_uu)
        endif
        p%ugu = ugu0+u0gu
      endif
! ugu2
      if (lpenc_loc(i_ugu2)) call fatal_error('calc_pencils_hydro_linearized:','does not calculate ugu2 pencil')
! u3u21, u1u32, u2u13, u2u31, u3u12, u1u23
      if (lpenc_loc(i_u3u21).or.lpenc_loc(i_u1u32).or.lpenc_loc(i_u2u13)  &
        .or.lpenc_loc(i_u2u31).or.lpenc_loc(i_u3u12).or.lpenc_loc(i_u1u23) ) then
        call fatal_error('calc_pencils_hydro_linearized:','ujukl pencils not calculated')
      endif
! del4u, del6u, and del4graddivu
      if (lpenc_loc(i_del4u)) call del4v(f,iuu,p%del4u)
      if (lpenc_loc(i_del6u)) call del6v(f,iuu,p%del6u)
      if (lpenc_loc(i_del6u_strict)) call del6v(f,iuu,p%del6u_strict,LSTRICT=.true.)
      if (lpenc_loc(i_del4graddivu)) call del4graddiv(f,iuu,p%del4graddivu)
! del6u_bulk
      if (lpenc_loc(i_del6u_bulk)) then
        call der6(f,iux,tmp,1)
        p%del6u_bulk(:,1)=tmp
        call der6(f,iuy,tmp,2)
        p%del6u_bulk(:,2)=tmp
        call der6(f,iuz,tmp,3)
        p%del6u_bulk(:,3)=tmp
      endif
!
! del2u, graddivu
!
      if (.not.lcartesian_coords.or.lalways_use_gij_etc) then
        if (lpenc_loc(i_graddivu)) then
          if (headtt.or.ldebug) print*,'calc_pencils_hydro: call gij_etc'
          call gij_etc(f,iuu,p%uu,p%uij,p%oij,GRADDIV=p%graddivu)
        endif
        if (lpenc_loc(i_del2u)) then
          call curl_mn(p%oij,p%curlo,p%oo)
          p%del2u=p%graddivu-p%curlo
        endif
      else
!
!  all 3 together
!
        if (lpenc_loc(i_del2u).and.lpenc_loc(i_graddivu).and.lpenc_loc(i_curlo)) then
          call del2v_etc(f,iuu,DEL2=p%del2u,GRADDIV=p%graddivu,CURLCURL=p%curlo)
!
!  all 3 possible pairs
!
        elseif (lpenc_loc(i_del2u).and.lpenc_loc(i_graddivu)) then
          call del2v_etc(f,iuu,DEL2=p%del2u,GRADDIV=p%graddivu)
        elseif (lpenc_loc(i_del2u).and.lpenc_loc(i_curlo)) then
          call del2v_etc(f,iuu,DEL2=p%del2u,CURLCURL=p%curlo)
        elseif (lpenc_loc(i_graddivu).and.lpenc_loc(i_curlo)) then
          call del2v_etc(f,iuu,GRADDIV=p%graddivu,CURLCURL=p%curlo)
!
!  all 3 individually
!
        elseif (lpenc_loc(i_del2u)) then
          call del2v_etc(f,iuu,DEL2=p%del2u)
        elseif (lpenc_loc(i_graddivu)) then
          call del2v_etc(f,iuu,GRADDIV=p%graddivu)
        elseif (lpenc_loc(i_curlo)) then
          call del2v_etc(f,iuu,CURLCURL=p%curlo)
        endif
      endif
!
!del2uj, d^u/dx^2 etc
!
      if (lpenc_loc(i_d2uidxj)) &
        call d2fi_dxj(f,iuu,p%d2uidxj)
!
! deluidxjk
!
      if (lpenc_loc(i_uijk)) then
        call del2fi_dxjk(f,iuu,p%uijk)
      endif
!
! grad5divu
!
      if (lpenc_loc(i_grad5divu)) then
        do i=1,3
          tmp=0.0
          do j=1,3
            ju=iuu+j-1
            call der5i1j(f,ju,tmp2,i,j)
            tmp=tmp+tmp2
          enddo
          p%grad5divu(:,i)=tmp
        enddo
      endif
!
! transpurho
!
      if (lpenc_loc(i_transpurho)) &
        call fatal_error('calc_pencils_hydro_linearized:','no linearized weno transport ')
!
    endsubroutine calc_pencils_hydro_linearized
!***********************************************************************
    subroutine hydro_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!  15-dec-10/MR: adapted from density for homogeneity
!  19-oct-15/ccyang: add calculation of the vorticity field.
!
      use Sub, only: curl
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(nx,3) :: pv
!
      real, dimension (mx,mz) :: fsum_tmp_cyl
      real, dimension (mx,my) :: fsum_tmp_sph
      real, dimension (mx) :: uphi
      real :: nygrid1,nzgrid1
!
!  Remove mean momenta or mean flows if desired.
!  Useful to avoid unphysical winds, for example in shearing box simulations.
!
      if (lrmv) then
        if (lremove_mean_momenta) then
          call remove_mean_momenta(f,iux,ilnrho)
        else
          if (lremove_mean_flow) call remove_mean_flow(f,iux)
          !if (lremove_mean_flow) call remove_mean_value(f,iux,iuz)  !(could use this one)
          if (lremove_mean_angmom) call remove_mean_angmom(f,iuz)
        endif
      endif
!
!  Calculate the vorticity field if required.
!
      if (ioo /= 0) then
        do n = n1, n2
          do m = m1, m2
            call curl(f, iux, pv)
            f(l1:l2,m,n,iox:ioz) = pv
          enddo 
        enddo
      endif
!
!    Slope limited diffusion: update characteristic speed
!    Not staggered yet, happpens later
!
     if (lslope_limit_diff .and. llast) then
!     if (lslope_limit_diff) then
       do m=1,my
       do n=1,mz
           f(:,m,n,isld_char)=w_sldchar_hyd* &
!            (f(:,m,n,iux)**2.+f(:,m,n,iuy)**2.+f(:,m,n,iuz)**2.)
            sqrt(f(:,m,n,iux)**2.+f(:,m,n,iuy)**2.+f(:,m,n,iuz)**2.)
       enddo
       enddo
     endif
!
!  For FARGO (orbital advection) algorithm.
!  Calculate the average velocity at the first sub-timestep.
!
     if ((lfargo_advection.and.lfirst).or.lcalc_uuavg) then
!
!  Pre-calculate the average large scale (= phi-averaged) speed of the flow
!
        if (lcylindrical_coords) then
          fsum_tmp_cyl=0.
          nygrid1=1./nygrid
          do n=1,mz
            do m=m1,m2
              uphi=f(:,m,n,iuy)
              fsum_tmp_cyl(:,n)=fsum_tmp_cyl(:,n)+uphi*nygrid1
            enddo
          enddo
!
! The sum has to be done processor-wise
! Sum over processors of same ipz, and different ipy
! --only relevant for 3D, but is here for generality
!
          call mpiallreduce_sum(fsum_tmp_cyl,uu_average_cyl,&
               (/mx,mz/),idir=2) !idir=2 is equal to old LSUMY=.true.
!
        elseif (lspherical_coords) then
          nzgrid1=1./nzgrid
          fsum_tmp_sph=0.
          do n=n1,n2
            do m=1,my
              uphi=f(:,m,n,iuz)
              fsum_tmp_sph(:,m)=fsum_tmp_sph(:,m)+uphi*nzgrid1
            enddo
          enddo
          call mpiallreduce_sum(fsum_tmp_sph,uu_average_sph,&
               (/mx,my/),idir=3) !idir=2 is equal to old LSUMY=.true.
!
        endif
     endif
!
    endsubroutine hydro_before_boundary
!***********************************************************************
    subroutine update_char_vel_hydro(f)
!
!   25-sep-15/MR+joern: coded
!   23-dec-15/joern: changed to staggered_max_vec
!   26-jan-16/MR: removed if clause around call as it is already in the caller
!
!   calculation of characteristic velocity
!   for slope limited diffusion
!
      use General, only: staggered_mean_vec,staggered_max_vec

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!      call staggered_mean_vec(f,iux,iFF_char_c,w_sldchar_hyd)
      call staggered_max_vec(f,iux,iFF_char_c,w_sldchar_hyd)
!
    endsubroutine update_char_vel_hydro
!***********************************************************************
    subroutine duu_dt(f,df,p)
!
!  velocity evolution
!  calculate du/dt = - u.gradu - 2Omega x u + grav + Fvisc
!  pressure gradient force added in density and entropy modules.
!
!   7-jun-02/axel: incoporated from subroutine pde
!  10-jun-02/axel+mattias: added Coriolis force
!  23-jun-02/axel: glnrho and fvisc are now calculated in here
!  17-jun-03/ulf: ux2, uy2 and uz2 added as diagnostic quantities
!  27-jun-07/dhruba: differential rotation as subroutine call
!  12-apr-16/MR: changes for Yin-Yang: only yz slices at the moment!
!  22-apr-16/MR: Changed calls to zsum_mn_name_xy for compatibility with Yin-Yang grid.
!
      use Diagnostics
      use Special, only: special_calc_hydro
      use Sub, only: dot, dot2, identify_bcs, cross, multsv, multsv_mn_add
      use General, only: transform_thph_yy, notanumber
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: f,df

      real, dimension (nx,3) :: tmpv
      real, dimension (nx) :: ftot
      integer :: j, ju
!
      Fmax=tini
!
!  Identify module and boundary conditions.
!
      call timing('duu_dt','entered',mnloop=.true.)
      if (headtt.or.ldebug) print*,'duu_dt: SOLVE'
      if (headtt) then
        call identify_bcs('ux',iux)
        call identify_bcs('uy',iuy)
        call identify_bcs('uz',iuz)
        if (lslope_limit_diff) then
          call identify_bcs('sld_char',isld_char)
        endif
      endif
!
!  Advection term.
!
      if (ladvection_velocity) then
        if (.not. lweno_transport .and. &
            .not. lno_meridional_flow .and. .not. lfargo_advection) &
!
!  Subtract u.gradu. In the relativistic case, this is automatically
!  u.gradS, and then we also need to subtract (divu)*S.
!
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-p%ugu
      if (ldensity.and.lrelativistic) then
        call multsv(p%divu,p%ss_rel,tmpv)
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-tmpv
      endif
!
!  WENO transport.
!
        if (lweno_transport) then
          do j=1,3
            df(l1:l2,m,n,iux-1+j)=df(l1:l2,m,n,iux-1+j) &
                - (p%transpurho(:,j) - p%uu(:,j)*p%transprho)*p%rho1
          enddo
        endif
!
!  No meridional flow : turn off the meridional flow (in spherical)
!  useful for debugging.
!
!  18-Mar-2010/AJ: this should probably go in a special module.
!  12-Mar-2017/WL: Agree, looks very specific.
!
        if (lno_meridional_flow) then
          f(l1:l2,m,n,iux:iuy)=0.0
          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-p%ugu(:,3)
        endif
!
        if (lfargo_advection) &
          df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-p%uuadvec_guu

      endif
!
!  Coriolis force, -2*Omega x u (unless lprecession=T)
!  Omega=(-sin_theta, 0, cos_theta), where theta corresponds to
!  colatitude. theta=0 places the box to the north pole and theta=90
!  the equator. Cartesian coordinates (x,y,z) now correspond to
!  (theta,phi,r) i.e. (south,east,up), in spherical polar coordinates
!
      if (Omega/=0.) then
!
        if (lcylindrical_coords) then
          call coriolis_cylindrical(df,p)
        elseif (lspherical_coords) then
          call coriolis_spherical(df,p)
        elseif (lprecession) then
          call precession(df,p)
        elseif (lrotation_xaxis) then
          call coriolis_cartesian_xaxis(df,p%uu,iux)
        else
          call coriolis_cartesian(df,p%uu,iux)
        endif
!
      endif
!
!  Coriolis force with in Cartesian domain with Omega=Omega(x)
!
      if (lcoriolis_xdep) call coriolis_xdep(df,p)
!
!  Interface for your personal subroutines calls
!
      if (lspecial) call special_calc_hydro(f,df,p)
!
! calculate viscous force
!
      if (lviscosity) call calc_viscous_force(df,p)
!
!  ``uu/dx'' for timestep
!
      if (lfirst.and.ldt.and.ladvection_velocity) then
        if (lmaximal_cdt) then
          advec_uu=max(abs(p%uu(:,1))*dline_1(:,1),&
                       abs(p%uu(:,2))*dline_1(:,2),&
                       abs(p%uu(:,3))*dline_1(:,3))
        else
          if (lfargo_advection) then
            advec_uu=sum(abs(p%uu_advec)*dline_1,2)
          else
            advec_uu=sum(abs(p%uu)*dline_1,2)
          endif
        endif
!
!  Empirically, it turns out that we need to take the full 3-D velocity
!  into account for computing the time step. It is not clear why.
!  Wlad finds that, otherwise, the code blows up for 2-D disk problem.
!  Some 1D and 2D samples work when the non-existent direction has the
!  largest velocity (like a 2D rz slice of a Keplerian disk that rotates
!  on the phi direction).
!
        if (lisotropic_advection) then
          if (dimensionality<3) advec_uu=sqrt(p%u2*dxyz_2)
        endif

        if (notanumber(advec_uu)) print*, 'advec_uu   =',advec_uu
        maxadvec=maxadvec+advec_uu
        if (headtt.or.ldebug) print*,'duu_dt: max(advec_uu) =',maxval(advec_uu)
      endif
!
!  Ekman Friction, used only in two dimensional runs.
!
      if (ekman_friction/=0) &
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-ekman_friction*p%uu
!
!  Boussinesq approximation: -g_z*alpha*(T-T_0) added.
!  Use Rayleigh number only with ltemperature.
!  Note: the buoyancy term is currently scaled with Ra*Pr, but Pr is also
!  regulated though mu and K, so Pr should eventually be eliminated.
!
      if (lboussinesq.and.ltemperature) then
        if (lsphere_in_a_box) then
          do j=1,3
            ju=j+iuu-1
            df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)+ &
            p%r_mn*Ra*Pr*f(l1:l2,m,n,iTT)*p%evr(:,j)
          enddo
        else
          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+Ra*Pr*f(l1:l2,m,n,iTT) !  gravity in the opposite z direction
        endif
      endif
!
!  Add possibility of forcing that is not delta-correlated in time.
!
      if (lforcing_cont_uu) &
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+ &
                              ampl_fcont_uu*p%fcont(:,:,1)
!
!  Damp motions in some regions for some time spans if desired.
!  For geodynamo: addition of dampuint evaluation.
!
      if (tdamp/=0.or.dampuext/=0.or.dampuint/=0) call udamping(f,df,p)
!
!  adding differential rotation via a frictional term
!
      if (tau_diffrot1/=0.) &
        call impose_profile_diffrot(f,df,uuprof,ldiffrot_test)
!
!  impose velocity ceiling
!
      if (velocity_ceiling/=0.) call impose_velocity_ceiling(f)
!
!  Save the advective derivative as an auxiliary.
!  MR: not the full story if there is a Lorentz force
!
      if (ladv_der_as_aux) then
        f(l1:l2,m,n,i_adv_derx:i_adv_derz) = p%fpres + p%fvisc
        if (lgrav) f(l1:l2,m,n,i_adv_derx:i_adv_derz) = f(l1:l2,m,n,i_adv_derx:i_adv_derz) + p%gg
      endif
!
!  Velocity in spherical coordinates from a Cartesian simulation
!  for sphere-in-a-box setups
!
      if (luu_sph_as_aux.and.lsphere_in_a_box) then
        f(l1:l2,m,n,iuu_sphr) = p%uu(:,1)*p%evr(:,1)+p%uu(:,2)*p%evr(:,2)+p%uu(:,3)*p%evr(:,3)
        f(l1:l2,m,n,iuu_spht) = p%uu(:,1)*p%evth(:,1)+p%uu(:,2)*p%evth(:,2)+p%uu(:,3)*p%evth(:,3)
        f(l1:l2,m,n,iuu_sphp) = p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy
      endif
!
!  Possibility to damp mean x momentum, ruxm, to zero.
!  This can be useful in situations where a mean flow is generated.
!  This tends to be the case when there is linear shear but no rotation
!  and the turbulence is forced. A constant drift velocity in the
!  x-direction is most dangerous, because it leads to a linear increase
!  of <uy> due to the shear term. If there is rotation, the epicyclic
!  motion brings one always back to no mean flow on the average.
!
      if (tau_damp_ruxm/=0.) &
        df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-ruxm*p%rho1*tau_damp_ruxm1
      if (tau_damp_ruym/=0.) &
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-ruym*p%rho1*tau_damp_ruym1
      if (tau_damp_ruzm/=0.) &
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-ruzm*p%rho1*tau_damp_ruzm1
!
!  Apply border profiles
!
      if (lborder_profiles) call set_border_hydro(f,df,p)
!
!  Fred: Option to constrain timestep for large forces
!
      if (lfirst.and.ldt.and.lcdt_tauf) then
        do j=1,3
          ftot=abs(df(l1:l2,m,n,iux+j-1))
          dt1_max=max(dt1_max,ftot/(cdt_tauf*ulev))
          Fmax=max(Fmax,ftot/ulev)
        enddo
      endif

      call calc_diagnostics_hydro(f,p)
!
      call timing('duu_dt','finished',mnloop=.true.)
!
    endsubroutine duu_dt
!*******************************************************************************
    subroutine calc_0d_diagnostics_hydro(f,p)
!
!   6-sep-19/MR: taken out from duu_dt
!
      use Diagnostics
      use Sub, only: dot, dot2, cross

      real, dimension(:,:,:,:) :: f
      type(pencil_case), intent(in) :: p
!
      real, dimension (nx,3) :: uxo
      real, dimension (nx) :: space_part_re,space_part_im,u2t,uot,out,fu
      real, dimension (nx) :: odel2um,uref,curlo2,qo,quxo,graddivu2
      real, dimension (nx,Nmodes_SH) :: urlm
      real, dimension (nx) :: rmask
      real :: kx
      integer :: k
      logical, save :: lcorr_zero_dt=.false.
!
!  Calculate maxima and rms values for diagnostic purposes
!
      call timing('calc_0d_diagnostics_hydro','before ldiagnos',mnloop=.true.)

      if (ldiagnos) then
        call sum_mn_name(p%u2,idiag_urms,lsqrt=.true.)
        if (idiag_durms/=0) then
          uref=ampluu(1)*cos(kx_uu*x(l1:l2))
          call sum_mn_name(p%u2-2.*p%uu(:,2)*uref+uref**2,idiag_durms)
        endif
        if (.not.lgpu) then
          if (ladvection_velocity.and.idiag_dtu/=0) call max_mn_name(advec_uu/cdt,idiag_dtu,l_dt=.true.)
          if (idiag_dtF/=0) call max_mn_name(Fmax/cdt_tauf,idiag_dtF,l_dt=.true.)
          if ((idiag_taufmin/=0).and.lcdt_tauf) &
            call max_mn_name(Fmax,idiag_taufmin,lreciprocal=.true.)
        endif
!
! urlm
!
        if (lspherical_coords.and.any(idiag_urlm/=0)) then
          call amp_lm(p%uu(:,1),urlm,profile_SH)    ! MR: tb restricted to needed urlm
          do k=1,Nmodes_SH
            call sum_mn_name(urlm(:,k),idiag_urlm(k))
          enddo
        endif
        if (idiag_urmsh/=0) then
          if (lequatory) call sum_mn_name_halfy(p%u2,idiag_urmsh)
          if (lequatorz) call sum_mn_name_halfz(p%u2,idiag_urmsh)
          fname(idiag_urmsn)=fname_half(idiag_urmsh,1)
          fname(idiag_urmss)=fname_half(idiag_urmsh,2)
          itype_name(idiag_urmsn)=ilabel_sum_sqrt
          itype_name(idiag_urmss)=ilabel_sum_sqrt
        endif
        if (idiag_urmsx/=0) call sum_mn_name(p%u2*xmask_hyd,idiag_urmsx,lsqrt=.true.)
        if (idiag_urmsz/=0) call sum_mn_name(p%u2*zmask_hyd(n-n1+1),idiag_urmsz,lsqrt=.true.)
        call max_mn_name(p%u2,idiag_umax,lsqrt=.true.)
        if (idiag_umin /=0) call max_mn_name(-sqrt(p%u2),idiag_umin,lneg=.true.)
        if (idiag_uxrms/=0) &
            call sum_mn_name(p%uu(:,1)**2,idiag_uxrms,lsqrt=.true.)
        if (idiag_uyrms/=0) &
            call sum_mn_name(p%uu(:,2)**2,idiag_uyrms,lsqrt=.true.)
        if (idiag_uzrms/=0) &
            call sum_mn_name(p%uu(:,3)**2,idiag_uzrms,lsqrt=.true.)
        if (idiag_uzrmaxs/=0) &
            call max_mn_name(p%uu(:,3)**2,idiag_uzrmaxs,lsqrt=.true.)
        if (idiag_uxmin/=0) call max_mn_name(-p%uu(:,1),idiag_uxmin,lneg=.true.)
        if (idiag_uymin/=0) call max_mn_name(-p%uu(:,2),idiag_uymin,lneg=.true.)
        if (idiag_uzmin/=0) call max_mn_name(-p%uu(:,3),idiag_uzmin,lneg=.true.)
        call max_mn_name(p%uu(:,1),idiag_uxmax)
        call max_mn_name(p%uu(:,2),idiag_uymax)
        call max_mn_name(p%uu(:,3),idiag_uzmax)
        if (idiag_rumax/=0) call max_mn_name(p%u2*p%rho**2,idiag_rumax,lsqrt=.true.)
        call sum_mn_name(p%ugu(:,1),idiag_uguxm)
        call sum_mn_name(p%ugu(:,2),idiag_uguym)
        call sum_mn_name(p%ugu(:,3),idiag_uguzm)
        call sum_mn_name(p%uij(:,1,1),idiag_dudx)
        call sum_mn_name(p%ugu2,idiag_ugu2m)
        if (idiag_ugurmsx/=0) call sum_mn_name(p%ugu2*xmask_hyd,idiag_ugurmsx,lsqrt=.true.)
        if (idiag_fkinzm/=0) call sum_mn_name(.5*p%rho*p%u2*p%uu(:,3),idiag_fkinzm)
        call max_mn_name(p%u2,idiag_um2)
        call sum_mn_name(p%u2,idiag_u2m)
        if (idiag_u2sphm/=0 .or. idiag_o2sphm/=0) then
          where (p%r_mn <= 1.)
            rmask = 1.
          elsewhere
            rmask = 0.
          endwhere
          call integrate_mn_name(rmask*p%u2,idiag_u2sphm)
          call integrate_mn_name(rmask*p%o2,idiag_o2sphm)
        endif
        call sum_mn_name(p%divu,idiag_divum)
        if (idiag_rdivum/=0)  call sum_mn_name(p%rho*p%divu,idiag_rdivum)
        if (idiag_divu2m/=0)  call sum_mn_name(p%divu**2,idiag_divu2m)
        if (idiag_gdivu2m/=0) then
          call dot2(p%graddivu,graddivu2)
          call sum_mn_name(graddivu2,idiag_gdivu2m)
        endif
        if (idiag_divrhourms/=0) call sum_mn_name((p%rho*p%divu+p%ugrho)**2,idiag_divrhourms,lsqrt=.true.)
        if (idiag_divrhoumax/=0) call max_mn_name(p%rho*p%divu+p%ugrho,idiag_divrhoumax)
        call sum_mn_name(p%uu(:,1),idiag_uxm)
        call sum_mn_name(p%uu(:,2),idiag_uym)
        call sum_mn_name(p%uu(:,3),idiag_uzm)
        if (idiag_ux2m/=0)    call sum_mn_name(p%uu(:,1)**2,idiag_ux2m)
        if (idiag_uy2m/=0)    call sum_mn_name(p%uu(:,2)**2,idiag_uy2m)
        if (idiag_uz2m/=0)    call sum_mn_name(p%uu(:,3)**2,idiag_uz2m)
        if (idiag_ux3m/=0)    call sum_mn_name(p%uu(:,1)**3,idiag_ux3m)
        if (idiag_uy3m/=0)    call sum_mn_name(p%uu(:,2)**3,idiag_uy3m)
        if (idiag_uz3m/=0)    call sum_mn_name(p%uu(:,3)**3,idiag_uz3m)
        if (idiag_ux4m/=0)    call sum_mn_name(p%uu(:,1)**4,idiag_ux4m)
        if (idiag_uy4m/=0)    call sum_mn_name(p%uu(:,2)**4,idiag_uy4m)
        if (idiag_uz4m/=0)    call sum_mn_name(p%uu(:,3)**4,idiag_uz4m)
        if (idiag_uxuy2m/=0)  call sum_mn_name(p%uu(:,1)**2*p%uu(:,2)**2,idiag_uxuy2m)
        if (idiag_uyuz2m/=0)  call sum_mn_name(p%uu(:,2)**2*p%uu(:,3)**2,idiag_uyuz2m)
        if (idiag_uzux2m/=0)  call sum_mn_name(p%uu(:,3)**2*p%uu(:,1)**2,idiag_uzux2m)
        if (idiag_ux2ccm/=0)  call sum_mn_name(c2z(n)*p%uu(:,1)**2,idiag_ux2ccm)
        if (idiag_ux2ssm/=0)  call sum_mn_name(s2z(n)*p%uu(:,1)**2,idiag_ux2ssm)
        if (idiag_uy2ccm/=0)  call sum_mn_name(c2z(n)*p%uu(:,2)**2,idiag_uy2ccm)
        if (idiag_uy2ssm/=0)  call sum_mn_name(s2z(n)*p%uu(:,2)**2,idiag_uy2ssm)
        if (idiag_uxuycsm/=0) &
            call sum_mn_name(cz(n)*sz(n)*p%uu(:,1)*p%uu(:,2),idiag_uxuycsm)
        if (idiag_uxuym/=0)   call sum_mn_name(p%uu(:,1)*p%uu(:,2),idiag_uxuym)
        if (idiag_uxuzm/=0)   call sum_mn_name(p%uu(:,1)*p%uu(:,3),idiag_uxuzm)
        if (idiag_uyuzm/=0)   call sum_mn_name(p%uu(:,2)*p%uu(:,3),idiag_uyuzm)
        if (idiag_ruxuym/=0) &
            call sum_mn_name(p%rho*p%uu(:,1)*p%uu(:,2),idiag_ruxuym)
        if (idiag_ruxuzm/=0) &
            call sum_mn_name(p%rho*p%uu(:,1)*p%uu(:,3),idiag_ruxuzm)
        if (idiag_ruyuzm/=0) &
            call sum_mn_name(p%rho*p%uu(:,2)*p%uu(:,3),idiag_ruyuzm)
        if (idiag_divuHrms/=0) call sum_mn_name((p%uij(:,1,1)+p%uij(:,2,2))**2,idiag_divuHrms,lsqrt=.true.)
        if (idiag_uxxrms/=0) call sum_mn_name(p%uij(:,1,1)**2,idiag_uxxrms,lsqrt=.true.)
        if (idiag_uyyrms/=0) call sum_mn_name(p%uij(:,2,2)**2,idiag_uyyrms,lsqrt=.true.)
        if (idiag_uxzrms/=0) call sum_mn_name(p%uij(:,1,3)**2,idiag_uxzrms,lsqrt=.true.)
        if (idiag_uyzrms/=0) call sum_mn_name(p%uij(:,2,3)**2,idiag_uyzrms,lsqrt=.true.)
        if (idiag_uzyrms/=0) call sum_mn_name(p%uij(:,3,2)**2,idiag_uzyrms,lsqrt=.true.)
        if (idiag_duxdzma/=0) call sum_mn_name(abs(p%uij(:,1,3)),idiag_duxdzma)
        if (idiag_duydzma/=0) call sum_mn_name(abs(p%uij(:,2,3)),idiag_duydzma)
!
        if (idiag_rux2m/=0) call sum_mn_name(p%rho*p%uu(:,1)**2,idiag_rux2m)
        if (idiag_ruy2m/=0) call sum_mn_name(p%rho*p%uu(:,2)**2,idiag_ruy2m)
        if (idiag_ruz2m/=0) call sum_mn_name(p%rho*p%uu(:,3)**2,idiag_ruz2m)
        if (idiag_ekin/=0) call sum_mn_name(p%ekin,idiag_ekin)
        if (idiag_EEK/=0) call sum_mn_name(p%ekin,idiag_EEK)
        if (idiag_ekintot/=0) call integrate_mn_name(p%ekin,idiag_ekintot)
        if (idiag_totangmom/=0) &
            call sum_lim_mn_name(p%rho*(p%uu(:,2)*x(l1:l2)-p%uu(:,1)*y(m)),&
            idiag_totangmom,p)
        if (idiag_uxglnrym/=0) &
            call sum_mn_name(p%uu(:,1)*p%glnrho(:,2),idiag_uxglnrym)
        if (idiag_uyglnrxm/=0) &
            call sum_mn_name(p%uu(:,2)*p%glnrho(:,1),idiag_uyglnrxm)
        if (idiag_uzdivum/=0) call sum_mn_name(p%uu(:,3)*p%divu,idiag_uzdivum)
        if (idiag_uxuydivum/=0) &
            call sum_mn_name(p%uu(:,1)*p%uu(:,2)*p%divu,idiag_uxuydivum)
!
        if (idiag_udpxxm/=0) call sum_mn_name(2.*(p%uu(:,1)*p%fpres(:,1)),idiag_udpxxm)
        if (idiag_udpyym/=0) call sum_mn_name(2.*(p%uu(:,2)*p%fpres(:,2)),idiag_udpyym)
        if (idiag_udpzzm/=0) call sum_mn_name(2.*(p%uu(:,3)*p%fpres(:,3)),idiag_udpzzm)
        if (idiag_udpxym/=0) call sum_mn_name(p%uu(:,1)*p%fpres(:,2)+p%uu(:,2)*p%fpres(:,1),idiag_udpxym)
        if (idiag_udpyzm/=0) call sum_mn_name(p%uu(:,2)*p%fpres(:,3)+p%uu(:,3)*p%fpres(:,2),idiag_udpyzm)
        if (idiag_udpxzm/=0) call sum_mn_name(p%uu(:,1)*p%fpres(:,3)+p%uu(:,3)*p%fpres(:,1),idiag_udpxzm)
!
!  Velocity components at one point (=pt).
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_uxpt/=0) call save_name(p%uu(lpoint-nghost,1),idiag_uxpt)
          if (idiag_uypt/=0) call save_name(p%uu(lpoint-nghost,2),idiag_uypt)
          if (idiag_uzpt/=0) call save_name(p%uu(lpoint-nghost,3),idiag_uzpt)
          if (idiag_uxuypt/=0) call save_name(p%uu(lpoint-nghost,1)*p%uu(lpoint-nghost,2),idiag_uxuypt)
          if (idiag_uyuzpt/=0) call save_name(p%uu(lpoint-nghost,2)*p%uu(lpoint-nghost,3),idiag_uyuzpt)
          if (idiag_uzuxpt/=0) call save_name(p%uu(lpoint-nghost,3)*p%uu(lpoint-nghost,1),idiag_uzuxpt)
        endif
!
!  Velocity components at point 2 (=p2).
!
        if (lroot.and.m==mpoint2.and.n==npoint2) then
          if (idiag_uxp2/=0) call save_name(p%uu(lpoint2-nghost,1),idiag_uxp2)
          if (idiag_uyp2/=0) call save_name(p%uu(lpoint2-nghost,2),idiag_uyp2)
          if (idiag_uzp2/=0) call save_name(p%uu(lpoint2-nghost,3),idiag_uzp2)
        endif
!
!  Mean momenta.
!
        if (idiag_ruxm/=0) call sum_mn_name(p%rho*p%uu(:,1),idiag_ruxm)
        if (idiag_ruym/=0) call sum_mn_name(p%rho*p%uu(:,2),idiag_ruym)
        if (idiag_ruzm/=0) call sum_mn_name(p%rho*p%uu(:,3),idiag_ruzm)
        if (idiag_ruxtot/=0) call sum_mn_name(p%rho*abs(p%uu(:,1)),idiag_ruxtot)
!
!  Mean angular momenta.
!
        if (idiag_rlxm/=0) call sum_mn_name( &
            p%rho*(y(m)*p%uu(:,3)-z(n)*p%uu(:,2)),idiag_rlxm)
        if (idiag_rlym/=0) call sum_mn_name( &
            p%rho*(z(n)*p%uu(:,1)-x(l1:l2)*p%uu(:,3)),idiag_rlym)
        if (idiag_rlzm/=0) call sum_mn_name( &
            p%rho*(x(l1:l2)*p%uu(:,2)-y(m)*p%uu(:,1)),idiag_rlzm)
        if (idiag_rlx2m/=0) call sum_mn_name( &
            (p%rho*(y(m)*p%uu(:,3)-z(n)*p%uu(:,2)))**2,idiag_rlx2m)
        if (idiag_rly2m/=0) call sum_mn_name( &
            (p%rho*(z(n)*p%uu(:,1)-x(l1:l2)*p%uu(:,3)))**2,idiag_rly2m)
        if (idiag_rlz2m/=0) call sum_mn_name( &
            (p%rho*(x(l1:l2)*p%uu(:,2)-y(m)*p%uu(:,1)))**2,idiag_rlz2m)
!
!  Total angular momentum in spherical coordinates
!
        if (idiag_tot_ang_mom/=0) call integrate_mn_name( &
            p%rho*x(l1:l2)*sin(y(m))*p%uu(:,3),idiag_tot_ang_mom)
!
!  Mean dot product of forcing and velocity field, <f.u>.
!
        if (idiag_fum/=0) then
          call dot(p%fcont(:,:,1),p%uu,fu)
          call sum_mn_name(ampl_fcont_uu*fu,idiag_fum)
        endif
!
!  Mean dot product of forcing and velocity field, <f.u>.
!
  !     if (idiag_rufm/=0) then
  !       call dot(p%fcont(:,:,1),p%uu,fu)
  !       call sum_mn_name(ampl_fcont_uu*fu,idiag_rufm)
  !     endif
!
!  Things related to kinetic helicity.
!
        if (idiag_ou_int/=0) call integrate_mn_name(p%ou,idiag_ou_int)
        if (idiag_oum/=0) call sum_mn_name(p%ou,idiag_oum)
        if (idiag_oxum/=0) call sum_mn_name(p%oxu(:,1),idiag_oxum)
        if (idiag_ourms/=0) call sum_mn_name(p%ou**2,idiag_ourms,lsqrt=.true.)
        if (idiag_oxurms/=0) call sum_mn_name(p%oxu2,idiag_oxurms,lsqrt=.true.)
!
!  Things related to vorticity.
!
        if (idiag_oumh/=0) then
          if (lequatory) call sum_mn_name_halfy(p%ou,idiag_oumh)
          if (lequatorz) call sum_mn_name_halfz(p%ou,idiag_oumh)
          fname(idiag_oumn)=fname_half(idiag_oumh,1)
          fname(idiag_oums)=fname_half(idiag_oumh,2)
          itype_name(idiag_oumn)=ilabel_sum
          itype_name(idiag_oums)=ilabel_sum
        endif
        call sum_mn_name(p%o2,idiag_orms,lsqrt=.true.)
        if (idiag_ormsh/=0) then
          if (lequatory) call sum_mn_name_halfy(p%o2,idiag_ormsh)
          if (lequatorz) call sum_mn_name_halfz(p%o2,idiag_ormsh)
          fname(idiag_ormsn)=fname_half(idiag_ormsh,1)
          fname(idiag_ormss)=fname_half(idiag_ormsh,2)
          itype_name(idiag_ormsn)=ilabel_sum_sqrt
          itype_name(idiag_ormss)=ilabel_sum_sqrt
        endif
!
!  <o.del2u>
!
        if (idiag_odel2um/=0) then
          call dot(p%oo,p%del2u,odel2um)
          call sum_mn_name(odel2um,idiag_odel2um)
        endif
!
!  various vorticity diagnostics
!
        if (idiag_omax/=0) call max_mn_name(p%o2,idiag_omax,lsqrt=.true.)
        if (idiag_o2m/=0)  call sum_mn_name(p%o2,idiag_o2m)
        if (idiag_o2u2m/=0)call sum_mn_name(p%o2*p%u2,idiag_o2u2m)
        if (idiag_ox2m/=0) call sum_mn_name(p%oo(:,1)**2,idiag_ox2m)
        if (idiag_oy2m/=0) call sum_mn_name(p%oo(:,2)**2,idiag_oy2m)
        if (idiag_oz2m/=0) call sum_mn_name(p%oo(:,3)**2,idiag_oz2m)
        if (idiag_ox3m/=0) call sum_mn_name(p%oo(:,1)**3,idiag_ox3m)
        if (idiag_oy3m/=0) call sum_mn_name(p%oo(:,2)**3,idiag_oy3m)
        if (idiag_oz3m/=0) call sum_mn_name(p%oo(:,3)**3,idiag_oz3m)
        if (idiag_ox4m/=0) call sum_mn_name(p%oo(:,1)**4,idiag_ox4m)
        if (idiag_oy4m/=0) call sum_mn_name(p%oo(:,2)**4,idiag_oy4m)
        if (idiag_oz4m/=0) call sum_mn_name(p%oo(:,3)**4,idiag_oz4m)
        call sum_mn_name(p%oo(:,1),idiag_oxm)
        call sum_mn_name(p%oo(:,2),idiag_oym)
        call sum_mn_name(p%oo(:,3),idiag_ozm)
        if (idiag_oxuzxm/=0) call sum_mn_name(p%oo(:,1)*p%uij(:,3,1),idiag_oxuzxm)
        if (idiag_oyuzym/=0) call sum_mn_name(p%oo(:,2)*p%uij(:,3,2),idiag_oyuzym)
        if (idiag_oxoym/=0) call sum_mn_name(p%oo(:,1)*p%oo(:,2),idiag_oxoym)
        if (idiag_oxozm/=0) call sum_mn_name(p%oo(:,1)*p%oo(:,3),idiag_oxozm)
        if (idiag_oyozm/=0) call sum_mn_name(p%oo(:,2)*p%oo(:,3),idiag_oyozm)
        if (idiag_pvzm/=0) call sum_mn_name((p%oo(:,3) + 2.*Omega)/p%rho,idiag_pvzm)
!
!  diagnostics involving curlo [ =curl(omega) ]
!
        if (idiag_q2m/=0 .or. idiag_qrms/=0 .or. idiag_qmax/=0 ) then
          call dot2(p%curlo,curlo2)
          call sum_mn_name(curlo2,idiag_q2m)
          call sum_mn_name(curlo2,idiag_qrms,lsqrt=.true.)
          call max_mn_name(curlo2,idiag_qmax,lsqrt=.true.)
        endif
!
!  <q.o>
!
        if (idiag_qom/=0) then
          call dot(p%curlo,p%oo,qo)
          call sum_mn_name(qo,idiag_qom)
        endif
!
!  <q.(uxo)>
!
        if (idiag_quxom/=0) then
          call cross(p%uu,p%oo,uxo)
          call dot(p%curlo,uxo,quxo)
          call sum_mn_name(quxo,idiag_quxom)
        endif
!
!  Mach number, rms and max
!
        call sum_mn_name(p%Ma2,idiag_Marms,lsqrt=.true.)
        call max_mn_name(p%Ma2,idiag_Mamax,lsqrt=.true.)
!
!  Diagonal components of alpha using FOSA:
!    alp11=<u3*u2,1>-<u2*u3,1>
!    alp22=<u1*u3,2>-<u3*u1,2>
!    alp33=<u2*u1,3>-<u1*u2,3>
!  For fully periodic domains it is sufficient to compute, e.g., only:
!    alp11=<u3*u2,1>,  alp22=<u1*u3,2>,  alp33=<u2*u1,3>
!
        call sum_mn_name(p%u3u21,idiag_u3u21m)
        call sum_mn_name(p%u1u32,idiag_u1u32m)
        call sum_mn_name(p%u2u13,idiag_u2u13m)
        call sum_mn_name(p%u2u31,idiag_u2u31m)
        call sum_mn_name(p%u3u12,idiag_u3u12m)
        call sum_mn_name(p%u1u23,idiag_u1u23m)
!
! fourier amplitude f(t) for non-axisymmetric waves:
!         u_x = f(t)*exp[i(kx*x+ky*y+kz*z)]
!
        if (idiag_uxfampm/=0 .or. idiag_uyfampm/=0 .or. idiag_uzfampm/=0 .or.&
            idiag_uxfampim/=0 .or. idiag_uxfampim/=0 .or. idiag_uzfampim/=0) then
          kx = kx_uu + qshear*Omega*ky_uu*t
          space_part_re = cos(kx*x(l1:l2)+ky_uu*y(m)+kz_uu*z(n))
          space_part_im = -sin(kx*x(l1:l2)+ky_uu*y(m)+kz_uu*z(n))
          if (idiag_uxfampm/=0) &
            call sum_mn_name(p%uu(:,1)*space_part_re,idiag_uxfampm)
          if (idiag_uyfampm/=0) &
            call sum_mn_name(p%uu(:,2)*space_part_re,idiag_uyfampm)
          if (idiag_uzfampm/=0) &
            call sum_mn_name(p%uu(:,3)*space_part_re,idiag_uzfampm)
          if (idiag_uxfampim/=0) &
            call sum_mn_name(p%uu(:,1)*space_part_im,idiag_uxfampim)
          if (idiag_uyfampim/=0) &
            call sum_mn_name(p%uu(:,2)*space_part_im,idiag_uyfampim)
          if (idiag_uzfampim/=0) &
            call sum_mn_name(p%uu(:,3)*space_part_im,idiag_uzfampim)
!
        endif
!
!  integrate velocity in time, to calculate correlation time later
!
        if (idiag_u2tm/=0) then
          call dot(p%uu,f(l1:l2,m,n,iuxt:iuzt),u2t)
          call sum_mn_name(u2t,idiag_u2tm)
        endif
!
!  integrate velocity in time, to calculate correlation time later
!
        if (idiag_outm/=0) then
          call dot(p%oo,f(l1:l2,m,n,iuxt:iuzt),out)
          call sum_mn_name(out,idiag_outm)
        endif
!
!  integrate velocity in time, to calculate correlation time later
!
        if (idiag_uotm/=0) then
          call dot(p%uu,f(l1:l2,m,n,ioxt:iozt),uot)
          call sum_mn_name(uot,idiag_uotm)
        endif
!
        if (lfargo_advection.and.idiag_nshift/=0) then
          if (lcylindrical_coords) then
            !phidot=uu_average_cyl(:,n)*rcyl_mn1
            !nshift=phidot*dt*dy_1(m)
            !call max_mn_name(nshift,idiag_nshift)
            if (dt==0.) then
              lcorr_zero_dt=.true.
              call max_mn_name(uu_average_cyl(l1:l2,n)*rcyl_mn1*dy_1(m),idiag_nshift)
            else
              call max_mn_name(uu_average_cyl(l1:l2,n)*rcyl_mn1*dt*dy_1(m),idiag_nshift)
            endif
          elseif (lspherical_coords) then
            !mnghost=m-nghost
            !phidot=uu_average_sph(:,n)*rcyl_mn1  ! rcyl = r*sinth(m)
            !nshift=phidot*dt*dz_1(n)
            !call max_mn_name(nshift,idiag_nshift)
            if (dt==0.) then
              lcorr_zero_dt=.true.
              call max_mn_name(uu_average_sph(l1:l2,m)*rcyl_mn1*dz_1(n),idiag_nshift)
            else
              call max_mn_name(uu_average_sph(l1:l2,m)*rcyl_mn1*dt*dz_1(n),idiag_nshift)
            endif
          endif
        endif

      elseif (lcorr_zero_dt) then
!
!  Here all quantities should be updated the calculation of which requires dt which
!  is zero at the very first diagnostics output time. (Doesn't work for itorder=1.)
! 
        lcorr_zero_dt=.false.
        if (lroot) fname(idiag_nshift)=fname(idiag_nshift)*dt
      endif  ! if (ldiagnos)

    endsubroutine calc_0d_diagnostics_hydro
!*******************************************************************************
    subroutine calc_1d_diagnostics_hydro(f,p)
!
!   6-sep-19/MR: taken out from duu_dt
!
      use Diagnostics
      use Sub, only: dot, dot2, cross, multsv_mn_add

      real, dimension(:,:,:,:) :: f
      type(pencil_case) :: p

      real, dimension (nx,3) :: curlru
      real, dimension (nx) :: uus, curlru2, Remz, uzmask
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1.
!
      if (l1davgfirst) then
        if (idiag_fmasszmz/=0) call xysum_mn_name_z(p%rho*p%uu(:,3), idiag_fmasszmz)
        if (idiag_fkinzmz /=0) call xysum_mn_name_z(p%ekin*p%uu(:,3),idiag_fkinzmz )
        if (idiag_fkinxmx /=0) call yzsum_mn_name_x(p%ekin*p%uu(:,1),idiag_fkinxmx )
        call xysum_mn_name_z(p%uu(:,1),idiag_uxmz)
        call xysum_mn_name_z(p%uu(:,2),idiag_uymz)
        call xysum_mn_name_z(p%uu(:,3),idiag_uzmz)
        call xysum_mn_name_z(p%divu,idiag_divumz)
        if (idiag_uzdivumz/=0) call xysum_mn_name_z(p%uu(:,3)*p%divu,idiag_uzdivumz)
        call xysum_mn_name_z(p%oo(:,1),idiag_oxmz)
        call xysum_mn_name_z(p%oo(:,2),idiag_oymz)
        call xysum_mn_name_z(p%oo(:,3),idiag_ozmz)
        if (idiag_ox2mz/=0) call xysum_mn_name_z(p%oo(:,1)**2,idiag_ox2mz)
        if (idiag_oy2mz/=0) call xysum_mn_name_z(p%oo(:,2)**2,idiag_oy2mz)
        if (idiag_oz2mz/=0) call xysum_mn_name_z(p%oo(:,3)**2,idiag_oz2mz)
        call xzsum_mn_name_y(p%uu(:,1),idiag_uxmy)
        call xzsum_mn_name_y(p%uu(:,2),idiag_uymy)
        call xzsum_mn_name_y(p%uu(:,3),idiag_uzmy)
        call yzsum_mn_name_x(p%uu(:,1),idiag_uxmx)
        call yzsum_mn_name_x(p%uu(:,2),idiag_uymx)
        call yzsum_mn_name_x(p%uu(:,3),idiag_uzmx)
        if (idiag_ruxmx/=0) call yzsum_mn_name_x(p%rho*p%uu(:,1),idiag_ruxmx)
        if (idiag_ruymx/=0) call yzsum_mn_name_x(p%rho*p%uu(:,2),idiag_ruymx)
        if (idiag_ruzmx/=0) call yzsum_mn_name_x(p%rho*p%uu(:,3),idiag_ruzmx)
        if (idiag_rux2mx /= 0) call yzsum_mn_name_x(p%rho*p%uu(:,1)**2, idiag_rux2mx)
        if (idiag_ruy2mx /= 0) call yzsum_mn_name_x(p%rho*p%uu(:,2)**2, idiag_ruy2mx)
        if (idiag_ruz2mx /= 0) call yzsum_mn_name_x(p%rho*p%uu(:,3)**2, idiag_ruz2mx)
        if (idiag_ruxuymx/=0) call yzsum_mn_name_x(p%rho*p%uu(:,1)*p%uu(:,2),idiag_ruxuymx)
        if (idiag_ruxuzmx/=0) call yzsum_mn_name_x(p%rho*p%uu(:,1)*p%uu(:,3),idiag_ruxuzmx)
        if (idiag_ruyuzmx/=0) call yzsum_mn_name_x(p%rho*p%uu(:,2)*p%uu(:,3),idiag_ruyuzmx)
        if (idiag_ux2mz/=0) call xysum_mn_name_z(p%uu(:,1)**2,idiag_ux2mz)
        if (idiag_uy2mz/=0) call xysum_mn_name_z(p%uu(:,2)**2,idiag_uy2mz)
        if (idiag_uz2mz/=0) call xysum_mn_name_z(p%uu(:,3)**2,idiag_uz2mz)
        if (idiag_uzupmz/=0 .or. idiag_ruzupmz/=0 .or. idiag_uz2upmz/=0 .or. &
            idiag_fkinzupmz/=0 .or. idiag_Rxyupmz/=0 .or. idiag_Rxzupmz/=0 .or. &
            idiag_Ryzupmz/=0) then
          where (p%uu(:,3) > 0.)
            uus = p%uu(:,3)
            uzmask = p%uu(:,3)/abs(p%uu(:,3))
          elsewhere
            uus=0.
            uzmask = 0.
          endwhere
          call xysum_mn_name_z(uus,idiag_uzupmz)
          if (idiag_ruzupmz/=0) call xysum_mn_name_z(p%rho*uus,idiag_ruzupmz)
          if (idiag_uz2upmz/=0) call xysum_mn_name_z(uus**2,idiag_uz2upmz)
          if (idiag_fkinzupmz/=0) call xysum_mn_name_z(p%ekin*uus,idiag_fkinzupmz)
          if (idiag_Rxyupmz/=0) call &
              xysum_mn_name_z(uzmask*f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucy),idiag_Rxyupmz)
          if (idiag_Rxzupmz/=0) call &
              xysum_mn_name_z(uzmask*f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucz),idiag_Rxzupmz)
          if (idiag_Ryzupmz/=0) call &
              xysum_mn_name_z(uzmask*f(l1:l2,m,n,iuu_flucy)*f(l1:l2,m,n,iuu_flucz),idiag_Ryzupmz)
        endif
        if (idiag_ffdownmz/=0 .or. idiag_uzupmz/=0 .or. idiag_ruzupmz/=0 .or. &
            idiag_uz2upmz/=0 .or. idiag_fkinzupmz/=0 .or. idiag_Rxydownmz/=0 .or. &
            idiag_Rxzdownmz/=0 .or. idiag_Ryzdownmz/=0) then
          where (p%uu(:,3) < 0.)
            uus = p%uu(:,3)
            uzmask = -p%uu(:,3)/abs(p%uu(:,3))
          elsewhere
            uus = 0.
            uzmask = 0.
          endwhere
          if (idiag_ffdownmz/=0) call xysum_mn_name_z(-uus/abs(p%uu(:,3)),idiag_ffdownmz)
          call xysum_mn_name_z(uus,idiag_uzdownmz)
          if (idiag_ruzdownmz/=0) call xysum_mn_name_z(p%rho*uus,idiag_ruzdownmz)
          if (idiag_uz2downmz/=0) call xysum_mn_name_z(uus**2,idiag_uz2downmz)
          if (idiag_fkinzdownmz/=0) call xysum_mn_name_z(p%ekin*uus,idiag_fkinzdownmz)
          if (idiag_Rxydownmz/=0) call &
              xysum_mn_name_z(uzmask*f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucy),idiag_Rxydownmz)
          if (idiag_Rxzdownmz/=0) call &
              xysum_mn_name_z(uzmask*f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucz),idiag_Rxzdownmz)
          if (idiag_Ryzdownmz/=0) call &
              xysum_mn_name_z(uzmask*f(l1:l2,m,n,iuu_flucy)*f(l1:l2,m,n,iuu_flucz),idiag_Ryzdownmz)
        endif
!
!  mean squared velocity and vorticity
!
        call xysum_mn_name_z(p%u2,idiag_u2mz)
        call xysum_mn_name_z(p%o2,idiag_o2mz)
        if (idiag_divu2mz/=0) call xysum_mn_name_z(p%divu**2,idiag_divu2mz)
!
!  mean squared mass flux divergence
!
        if (idiag_divru2mz/=0) call xysum_mn_name_z((p%rho*p%divu+p%ugrho)**2,idiag_divru2mz)
!
!  mean squared curl of mass flux
!
        if (idiag_curlru2mz/=0) then
          call cross(p%grho,p%uu,curlru)
          call multsv_mn_add(p%rho,p%oo,curlru)
          call dot2(curlru,curlru2)
          call xysum_mn_name_z(curlru2,idiag_curlru2mz)
        endif

        if (idiag_ruxmz/=0) call xysum_mn_name_z(p%rho*p%uu(:,1),idiag_ruxmz)
        if (idiag_ruymz/=0) call xysum_mn_name_z(p%rho*p%uu(:,2),idiag_ruymz)
        if (idiag_ruzmz/=0) call xysum_mn_name_z(p%rho*p%uu(:,3),idiag_ruzmz)
        if (idiag_rux2mz/=0) call xysum_mn_name_z(p%rho*p%uu(:,1)**2,idiag_rux2mz)
        if (idiag_ruy2mz/=0) call xysum_mn_name_z(p%rho*p%uu(:,2)**2,idiag_ruy2mz)
        if (idiag_ruz2mz/=0) call xysum_mn_name_z(p%rho*p%uu(:,3)**2,idiag_ruz2mz)
        if (idiag_ux2my/=0) call xzsum_mn_name_y(p%uu(:,1)**2,idiag_ux2my)
        if (idiag_uy2my/=0) call xzsum_mn_name_y(p%uu(:,2)**2,idiag_uy2my)
        if (idiag_uz2my/=0) call xzsum_mn_name_y(p%uu(:,3)**2,idiag_uz2my)
        if (idiag_ux2mx/=0) call yzsum_mn_name_x(p%uu(:,1)**2,idiag_ux2mx)
        if (idiag_uy2mx/=0) call yzsum_mn_name_x(p%uu(:,2)**2,idiag_uy2mx)
        if (idiag_uz2mx/=0) call yzsum_mn_name_x(p%uu(:,3)**2,idiag_uz2mx)
        if (idiag_ox2mx/=0) call yzsum_mn_name_x(p%oo(:,1)**2,idiag_ox2mx)
        if (idiag_oy2mx/=0) call yzsum_mn_name_x(p%oo(:,2)**2,idiag_oy2mx)
        if (idiag_oz2mx/=0) call yzsum_mn_name_x(p%oo(:,3)**2,idiag_oz2mx)
        if (idiag_uxuymz/=0) call xysum_mn_name_z(p%uu(:,1)*p%uu(:,2),idiag_uxuymz)
        if (idiag_uxuzmz/=0) call xysum_mn_name_z(p%uu(:,1)*p%uu(:,3),idiag_uxuzmz)
        if (idiag_uyuzmz/=0) call xysum_mn_name_z(p%uu(:,2)*p%uu(:,3),idiag_uyuzmz)
        if (idiag_Rxymz/=0) call xysum_mn_name_z(f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucy),idiag_Rxymz)
        if (idiag_Rxzmz/=0) call xysum_mn_name_z(f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucz),idiag_Rxzmz)
        if (idiag_Ryzmz/=0) call xysum_mn_name_z(f(l1:l2,m,n,iuu_flucy)*f(l1:l2,m,n,iuu_flucz),idiag_Ryzmz)
        if (idiag_ruxuymz/=0) call xysum_mn_name_z(p%rho*p%uu(:,1)*p%uu(:,2),idiag_ruxuymz)
        if (idiag_ruxuzmz/=0) call xysum_mn_name_z(p%rho*p%uu(:,1)*p%uu(:,3),idiag_ruxuzmz)
        if (idiag_ruyuzmz/=0) call xysum_mn_name_z(p%rho*p%uu(:,2)*p%uu(:,3),idiag_ruyuzmz)
        if (idiag_ruxuy2mz/=0) call xysum_mn_name_z((p%rho*p%uu(:,1)*p%uu(:,2))**2,idiag_ruxuy2mz)
        if (idiag_ruxuz2mz/=0) call xysum_mn_name_z((p%rho*p%uu(:,1)*p%uu(:,3))**2,idiag_ruxuz2mz)
        if (idiag_ruyuz2mz/=0) call xysum_mn_name_z((p%rho*p%uu(:,2)*p%uu(:,3))**2,idiag_ruyuz2mz)
        if (idiag_oxuxxmz/=0) call xysum_mn_name_z(p%oo(:,1)*p%uij(:,1,1),idiag_oxuxxmz)
        if (idiag_oyuxymz/=0) call xysum_mn_name_z(p%oo(:,2)*p%uij(:,1,2),idiag_oyuxymz)
        if (idiag_oxuyxmz/=0) call xysum_mn_name_z(p%oo(:,1)*p%uij(:,2,1),idiag_oxuyxmz)
        if (idiag_oyuyymz/=0) call xysum_mn_name_z(p%oo(:,2)*p%uij(:,2,2),idiag_oyuyymz)
        if (idiag_oxuzxmz/=0) call xysum_mn_name_z(p%oo(:,1)*p%uij(:,3,1),idiag_oxuzxmz)
        if (idiag_oyuzymz/=0) call xysum_mn_name_z(p%oo(:,2)*p%uij(:,3,2),idiag_oyuzymz)
        if (idiag_uyxuzxmz/=0) call xysum_mn_name_z(p%uij(:,2,1)*p%uij(:,3,1),idiag_uyxuzxmz)
        if (idiag_uyyuzymz/=0) call xysum_mn_name_z(p%uij(:,2,2)*p%uij(:,3,2),idiag_uyyuzymz)
        if (idiag_uyzuzzmz/=0) call xysum_mn_name_z(p%uij(:,2,3)*p%uij(:,3,3),idiag_uyzuzzmz)
        if (idiag_uxuymy/=0) call xzsum_mn_name_y(p%uu(:,1)*p%uu(:,2),idiag_uxuymy)
        if (idiag_uxuzmy/=0) call xzsum_mn_name_y(p%uu(:,1)*p%uu(:,3),idiag_uxuzmy)
        if (idiag_uyuzmy/=0) call xzsum_mn_name_y(p%uu(:,2)*p%uu(:,3),idiag_uyuzmy)
        if (idiag_uxuymx/=0) call yzsum_mn_name_x(p%uu(:,1)*p%uu(:,2),idiag_uxuymx)
        if (idiag_uxuzmx/=0) call yzsum_mn_name_x(p%uu(:,1)*p%uu(:,3),idiag_uxuzmx)
        if (idiag_uyuzmx/=0) call yzsum_mn_name_x(p%uu(:,2)*p%uu(:,3),idiag_uyuzmx)
        call yzsum_mn_name_x(p%ekin,idiag_ekinmx)
        call xysum_mn_name_z(p%ekin,idiag_ekinmz)
        call yzsum_mn_name_x(p%ou,idiag_oumx)
        call xzsum_mn_name_y(p%ou,idiag_oumy)
        call xysum_mn_name_z(p%ou,idiag_oumz)
        call yzsum_mn_name_x(p%ugu(:,1),idiag_uguxmx)
        call yzsum_mn_name_x(p%ugu(:,2),idiag_uguymx)
        call yzsum_mn_name_x(p%ugu(:,3),idiag_uguzmx)
        call xzsum_mn_name_y(p%ugu(:,1),idiag_uguxmy)
        call xzsum_mn_name_y(p%ugu(:,2),idiag_uguymy)
        call xzsum_mn_name_y(p%ugu(:,3),idiag_uguzmy)
        call xysum_mn_name_z(p%ugu(:,1),idiag_uguxmz)
        call xysum_mn_name_z(p%ugu(:,2),idiag_uguymz)
        call xysum_mn_name_z(p%ugu(:,3),idiag_uguzmz)
        call xysum_mn_name_z(p%ogu(:,1),idiag_oguxmz)
        call xysum_mn_name_z(p%ogu(:,2),idiag_oguymz)
        call xysum_mn_name_z(p%ogu(:,3),idiag_oguzmz)
        if (idiag_ogux2mz/=0) call xysum_mn_name_z(p%ogu(:,1)**2,idiag_ogux2mz)
        if (idiag_oguy2mz/=0) call xysum_mn_name_z(p%ogu(:,2)**2,idiag_oguy2mz)
        if (idiag_oguz2mz/=0) call xysum_mn_name_z(p%ogu(:,3)**2,idiag_oguz2mz)
        if (idiag_oxdivumz/=0) call xysum_mn_name_z(p%oo(:,1)*p%divu,idiag_oxdivumz)
        if (idiag_oydivumz/=0) call xysum_mn_name_z(p%oo(:,2)*p%divu,idiag_oydivumz)
        if (idiag_ozdivumz/=0) call xysum_mn_name_z(p%oo(:,3)*p%divu,idiag_ozdivumz)
        if (idiag_oxdivu2mz/=0) call xysum_mn_name_z((p%oo(:,1)*p%divu)**2,idiag_oxdivu2mz)
        if (idiag_oydivu2mz/=0) call xysum_mn_name_z((p%oo(:,2)*p%divu)**2,idiag_oydivu2mz)
        if (idiag_ozdivu2mz/=0) call xysum_mn_name_z((p%oo(:,3)*p%divu)**2,idiag_ozdivu2mz)
        call xysum_mn_name_z(p%u3u21,idiag_u3u21mz)
        call xysum_mn_name_z(p%u1u32,idiag_u1u32mz)
        call xysum_mn_name_z(p%u2u13,idiag_u2u13mz)
        call xysum_mn_name_z(p%u2u31,idiag_u2u31mz)
        call xysum_mn_name_z(p%u3u12,idiag_u3u12mz)
        call xysum_mn_name_z(p%u1u23,idiag_u1u23mz)
!
        if (idiag_totalforcezmz/=0) then
          uus = p%rho*(p%fpres(:,3) + p%fvisc(:,3))
          if (lgrav) uus = uus + p%rho*p%gg(:,3)
          call xysum_mn_name_z(uus,idiag_totalforcezmz)
        endif
        if (idiag_totalforcezupmz/=0) then
          where (p%uu(:,3) > 0.) 
            uus = p%rho*(p%fpres(:,3) + p%fvisc(:,3) + p%gg(:,3))
          elsewhere
            uus=0.
          endwhere
          call xysum_mn_name_z(uus,idiag_totalforcezupmz)
        endif
        if (idiag_totalforcezdownmz/=0) then
          where (p%uu(:,3) < 0.)
            uus = p%rho*(p%fpres(:,3) + p%fvisc(:,3) + p%gg(:,3))
          elsewhere
            uus=0.
          endwhere
          call xysum_mn_name_z(uus,idiag_totalforcezdownmz)
        endif
!
        if (idiag_acczmz/=0) then
          uus = p%fpres(:,3) + p%fvisc(:,3)
          if (lgrav) uus = uus + p%gg(:,3)
          call xysum_mn_name_z(uus,idiag_acczmz)   ! yet incorrect for Yin-Yang
        endif
        if (idiag_acczupmz/=0) then
          where (p%uu(:,3) > 0.) 
            uus = p%fpres(:,3) + p%fvisc(:,3) + p%gg(:,3)
          elsewhere
            uus=0.
          endwhere
          call xysum_mn_name_z(uus,idiag_acczupmz)   ! yet incorrect for Yin-Yang
        endif
        if (idiag_acczdownmz/=0) then
          where (p%uu(:,3) < 0.)
            uus = p%fpres(:,3) + p%fvisc(:,3) + p%gg(:,3)
          elsewhere
            uus=0.
          endwhere
          call xysum_mn_name_z(uus,idiag_acczdownmz)   ! yet incorrect for Yin-Yang
        endif
!
        if (idiag_accpowzmz/=0) then
          uus = p%fpres(:,3) + p%fvisc(:,3)
          if (lgrav) uus = uus + p%gg(:,3)
          uus=p%uu(:,3)*uus
          call xysum_mn_name_z(uus,idiag_accpowzmz)   ! yet incorrect for Yin-Yang
        endif
        if (idiag_accpowzupmz/=0) then
          where (p%uu(:,3) > 0.) 
            uus = p%fpres(:,3) + p%fvisc(:,3) + p%gg(:,3)
            uus = p%uu(:,3)*uus
          elsewhere
            uus=0.
          endwhere
          call xysum_mn_name_z(uus,idiag_accpowzupmz)   ! yet incorrect for Yin-Yang
        endif
        if (idiag_accpowzdownmz/=0) then
          where (p%uu(:,3) < 0.)
            uus = p%fpres(:,3) + p%fvisc(:,3) + p%gg(:,3)
            uus = p%uu(:,3)*uus
          elsewhere
            uus=0.
          endwhere
          call xysum_mn_name_z(uus,idiag_accpowzdownmz)   ! yet incorrect for Yin-Yang
        endif
!
!  phi-z averages
!
        if (idiag_Remz/=0) then
          Remz = sqrt(p%ugu2/p%diffus_total**2)
          where (p%diffus_total < tini) Remz = 0.
          call xysum_mn_name_z(Remz,idiag_Remz)
        endif
        if (idiag_u2mr/=0)   call phizsum_mn_name_r(p%u2,idiag_u2mr)
        if (idiag_urmr/=0) &
            call phizsum_mn_name_r(p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy,idiag_urmr)
        if (idiag_upmr/=0) &
            call phizsum_mn_name_r(p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy,idiag_upmr)
        if (idiag_uzmr/=0) &
             call phizsum_mn_name_r(p%uu(:,3),idiag_uzmr)
        if (idiag_ormr/=0) &
            call phizsum_mn_name_r(p%oo(:,1)*p%pomx+p%oo(:,2)*p%pomy,idiag_ormr)
        if (idiag_opmr/=0) &
            call phizsum_mn_name_r(p%oo(:,1)*p%phix+p%oo(:,2)*p%phiy,idiag_opmr)
        if (idiag_ozmr/=0) &
             call phizsum_mn_name_r(p%oo(:,3),idiag_ozmr)
      endif

    endsubroutine calc_1d_diagnostics_hydro
!******************************************************************************
    subroutine calc_2d_diagnostics_hydro(f,p)
!
!   6-sep-19/MR: taken out from duu_dt
!
      use Diagnostics
      use Sub, only: dot, dot2, cross

      real, dimension(:,:,:,:) :: f
      type(pencil_case), intent(in) :: p
!
      real, dimension (nx) :: uus, uxmask
!
!  2-D averages.
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        if (idiag_urmphi/=0) &
            call phisum_mn_name_rz(p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy,idiag_urmphi)
        if (idiag_ur2mphi/=0) &
            call phisum_mn_name_rz((p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy)**2,idiag_ur2mphi)
        if (idiag_ursphmphi/=0) &
            call phisum_mn_name_rz(p%uu(:,1)*p%evr(:,1)+ &
              p%uu(:,2)*p%evr(:,2)+p%uu(:,3)*p%evr(:,3),idiag_ursphmphi)
        if (idiag_uthmphi/=0) &
            call phisum_mn_name_rz(p%uu(:,1)*p%evth(:,1)+ &
              p%uu(:,2)*p%evth(:,2)+p%uu(:,3)*p%evth(:,3),idiag_uthmphi)
        if (idiag_rursphmphi/=0) &
            call phisum_mn_name_rz(p%rho*(p%uu(:,1)*p%evr(:,1)+ &
              p%uu(:,2)*p%evr(:,2)+p%uu(:,3)*p%evr(:,3)),idiag_rursphmphi)
        if (idiag_ruthmphi/=0) &
            call phisum_mn_name_rz(p%rho*(p%uu(:,1)*p%evth(:,1)+ &
              p%uu(:,2)*p%evth(:,2)+p%uu(:,3)*p%evth(:,3)),idiag_ruthmphi)
        if (idiag_upmphi/=0) &
            call phisum_mn_name_rz(p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy,idiag_upmphi)
        if (idiag_up2mphi/=0) &
            call phisum_mn_name_rz((p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy)**2,idiag_up2mphi)
        call phisum_mn_name_rz(p%uu(:,3),idiag_uzmphi)
        call phisum_mn_name_rz(p%uu(:,3)**2,idiag_uz2mphi)
        call phisum_mn_name_rz(p%u2,idiag_u2mphi)
        if (idiag_urupmphi/=0) &
            call phisum_mn_name_rz((p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy)*(p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy),idiag_urupmphi)
        if (idiag_uruzmphi/=0) &
            call phisum_mn_name_rz((p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy)*p%uu(:,3),idiag_uruzmphi)
        if (idiag_upuzmphi/=0) &
            call phisum_mn_name_rz((p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy)*p%uu(:,3),idiag_upuzmphi)
        if (idiag_rurupmphi/=0) &
            call phisum_mn_name_rz(p%rho*(p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy)*(p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy),idiag_rurupmphi)
        if (idiag_ruruzmphi/=0) &
            call phisum_mn_name_rz(p%rho*(p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy)*p%uu(:,3),idiag_ruruzmphi)
        if (idiag_rupuzmphi/=0) &
            call phisum_mn_name_rz(p%rho*(p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy)*p%uu(:,3),idiag_rupuzmphi)
        if (idiag_rurmphi/=0) &
            call phisum_mn_name_rz(p%rho*(p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy),idiag_rurmphi)
        if (idiag_rupmphi/=0) &
            call phisum_mn_name_rz(p%rho*(p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy),idiag_rupmphi)
        if (idiag_rupmphi/=0) &
            call phisum_mn_name_rz(p%rho*p%uu(:,3),idiag_ruzmphi)
        call phisum_mn_name_rz(p%oo(:,3),idiag_ozmphi)
        call phisum_mn_name_rz(p%ou,idiag_oumphi)
        if (idiag_fkinrsphmphi/=0) &
            call phisum_mn_name_rz(p%ekin*(p%uu(:,1)*p%evr(:,1)+ &
              p%uu(:,2)*p%evr(:,2)+p%uu(:,3)*p%evr(:,3)),idiag_fkinrsphmphi)
!
        call ysum_mn_name_xz(p%uu(:,1),idiag_uxmxz)
        call ysum_mn_name_xz(p%uu(:,2),idiag_uymxz)
        call ysum_mn_name_xz(p%uu(:,3),idiag_uzmxz)
        if (idiag_ux2mxz/=0) &
            call ysum_mn_name_xz(p%uu(:,1)**2,idiag_ux2mxz)
        if (idiag_uy2mxz/=0) &
            call ysum_mn_name_xz(p%uu(:,2)**2,idiag_uy2mxz)
        if (idiag_uz2mxz/=0) &
            call ysum_mn_name_xz(p%uu(:,3)**2,idiag_uz2mxz)
        if (idiag_uxuymxz/=0) &
            call ysum_mn_name_xz(p%uu(:,1)*p%uu(:,2),idiag_uxuymxz)
        if (idiag_uxuzmxz/=0) &
            call ysum_mn_name_xz(p%uu(:,1)*p%uu(:,3),idiag_uxuzmxz)
        if (idiag_uyuzmxz/=0) &
            call ysum_mn_name_xz(p%uu(:,2)*p%uu(:,3),idiag_uyuzmxz)
        call ysum_mn_name_xz(p%ou,idiag_oumxz)
!
        call zsum_mn_name_xy(p%uu(:,1),idiag_uxmxy)
!
!  Changed calls for compatibility with Yin-Yang grid:
!  all non-scalars in which y or z components of a vector are used must
!  be treated as below,
!
        call zsum_mn_name_xy(p%uu,idiag_uymxy,(/0,1,0/))
        call zsum_mn_name_xy(p%uu,idiag_uzmxy,(/0,0,1/))
        call zsum_mn_name_xy(p%uu,idiag_uxuymxy,(/1,1,0/))
        call zsum_mn_name_xy(p%uu,idiag_uxuzmxy,(/1,0,1/))
        call zsum_mn_name_xy(p%uu,idiag_uyuzmxy,(/0,1,1/))
        if (idiag_Rxymxy/=0) call zsum_mn_name_xy(f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucy),idiag_Rxymxy)
        if (idiag_Rxzmxy/=0) call zsum_mn_name_xy(f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucz),idiag_Rxzmxy)
        if (idiag_Ryzmxy/=0) call zsum_mn_name_xy(f(l1:l2,m,n,iuu_flucy)*f(l1:l2,m,n,iuu_flucz),idiag_Ryzmxy)
        call zsum_mn_name_xy(p%oo(:,1),idiag_oxmxy)
        call zsum_mn_name_xy(p%oo,idiag_oymxy,(/0,1,0/))
        call zsum_mn_name_xy(p%oo,idiag_ozmxy,(/0,0,1/))
        call zsum_mn_name_xy(p%ou,idiag_oumxy)
        if (idiag_pvzmxy/=0) &
            call zsum_mn_name_xy((p%oo(:,3)+2.*Omega)/p%rho,idiag_pvzmxy)    ! yet incorrect for Yin-Yang
        if (idiag_ruxmxy/=0) call zsum_mn_name_xy(p%rho*p%uu(:,1),idiag_ruxmxy)
        call zsum_mn_name_xy(p%uu,idiag_ruymxy,(/0,1,0/),p%rho)
        call zsum_mn_name_xy(p%uu,idiag_ruzmxy,(/0,0,1/),p%rho)
        if (idiag_ux2mxy/=0) &
            call zsum_mn_name_xy(p%uu(:,1)**2,idiag_ux2mxy)
        call zsum_mn_name_xy(p%uu,idiag_uy2mxy,(/0,2,0/))
        call zsum_mn_name_xy(p%uu,idiag_uz2mxy,(/0,0,2/))
        if (idiag_rux2mxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,1)**2,idiag_rux2mxy)
        call zsum_mn_name_xy(p%uu,idiag_ruy2mxy,(/0,2,0/),p%rho)
        call zsum_mn_name_xy(p%uu,idiag_ruz2mxy,(/0,0,2/),p%rho)
        call zsum_mn_name_xy(p%uu,idiag_ruxuymxy,(/1,1,0/),p%rho)
        call zsum_mn_name_xy(p%uu,idiag_ruxuzmxy,(/1,0,1/),p%rho)
        call zsum_mn_name_xy(p%uu,idiag_ruyuzmxy,(/0,1,1/),p%rho)
        if (idiag_fkinxmxy/=0) &
            call zsum_mn_name_xy(p%ekin*p%uu(:,1),idiag_fkinxmxy)
        call zsum_mn_name_xy(p%uu,idiag_fkinymxy,(/0,1,0/),p%ekin)
        call zsum_mn_name_xy(p%ugu(:,1),idiag_uguxmxy)
        call zsum_mn_name_xy(p%ugu,idiag_uguymxy,(/0,1,0/))
        call zsum_mn_name_xy(p%ugu,idiag_uguzmxy,(/0,0,1/))
        if (idiag_uxupmxy/=0 .or. idiag_ruxupmxy/=0 .or. idiag_ux2upmxy/=0 .or. &
            idiag_fkinxupmxy/=0 .or. idiag_Rxyupmxy/=0 .or. idiag_Rxzupmxy/=0 .or. &
            idiag_Ryzupmxy/=0) then
          where (p%uu(:,1) > 0.)
            uus = p%uu(:,1)
            uxmask = p%uu(:,1)/abs(p%uu(:,1))
          elsewhere
            uus=0.
            uxmask = 0.
          endwhere
          call zsum_mn_name_xy(uus,idiag_uxupmxy)
          if (idiag_ruxupmxy/=0) call zsum_mn_name_xy(p%rho*uus,idiag_ruxupmxy)
          if (idiag_ux2upmxy/=0) call zsum_mn_name_xy(uus**2,idiag_ux2upmxy)
          if (idiag_fkinxupmxy/=0) call zsum_mn_name_xy(p%ekin*uus,idiag_fkinxupmxy)
          if (idiag_Rxyupmxy/=0) call &
              zsum_mn_name_xy(uxmask*f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucy),idiag_Rxyupmxy)
          if (idiag_Rxzupmxy/=0) call &
              zsum_mn_name_xy(uxmask*f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucz),idiag_Rxzupmxy)
          if (idiag_Ryzupmxy/=0) call &
              zsum_mn_name_xy(uxmask*f(l1:l2,m,n,iuu_flucy)*f(l1:l2,m,n,iuu_flucz),idiag_Ryzupmxy)
        endif
        if (idiag_ffdownmxy/=0  .or. idiag_uxdownmxy/=0 .or. idiag_ruxdownmxy/=0 .or. &
            idiag_ux2downmxy/=0 .or. idiag_fkinxdownmxy/=0 .or. idiag_Rxydownmxy/=0 .or. &
            idiag_Rxzdownmxy/=0 .or. idiag_Ryzdownmxy/=0) then
          where (p%uu(:,1) < 0.)
            uus = p%uu(:,1)
            uxmask = -p%uu(:,1)/abs(p%uu(:,1))
          elsewhere
            uus = 0.
            uxmask = 0.
          endwhere
          if (idiag_ffdownmxy/=0) call zsum_mn_name_xy(-uus/abs(p%uu(:,1)),idiag_ffdownmxy)
          Call zsum_mn_name_xy(uus,idiag_uxdownmxy)
          if (idiag_ruxdownmxy/=0) call zsum_mn_name_xy(p%rho*uus,idiag_ruxdownmxy)
          if (idiag_ux2downmxy/=0) call zsum_mn_name_xy(uus**2,idiag_ux2downmxy)
          if (idiag_fkinxdownmxy/=0) call zsum_mn_name_xy(p%ekin*uus,idiag_fkinxdownmxy)
          if (idiag_Rxydownmxy/=0) call &
              zsum_mn_name_xy(uxmask*f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucy),idiag_Rxydownmxy)
          if (idiag_Rxzdownmxy/=0) call &
              zsum_mn_name_xy(uxmask*f(l1:l2,m,n,iuu_flucx)*f(l1:l2,m,n,iuu_flucz),idiag_Rxzdownmxy)
          if (idiag_Ryzdownmxy/=0) call &
              zsum_mn_name_xy(uxmask*f(l1:l2,m,n,iuu_flucy)*f(l1:l2,m,n,iuu_flucz),idiag_Ryzdownmxy)
        endif
      else
!
!  idiag_uxmxy and idiag_uymxy also need to be calculated when
!  ldiagnos and idiag_umx and/or idiag_umy, so
!
!  We may need to calculate uxmxy without calculating umx. The following
!  if condition was messing up calculation of umxy_rms
!
        if (ldiagnos) then
          call zsum_mn_name_xy(p%uu(:,1),idiag_uxmxy)
          call zsum_mn_name_xy(p%uu,idiag_uymxy,(/0,1,0/))
          call zsum_mn_name_xy(p%uu,idiag_uzmxy,(/0,0,1/))
        endif
      endif

    endsubroutine calc_2d_diagnostics_hydro
!**************************************************************************************
    subroutine calc_diagnostics_hydro(f,p)

      use Slices_methods, only: store_slices
      use Sub, only: vecout

      real, dimension(:,:,:,:) :: f
      type(pencil_case), intent(in) :: p

      call calc_2d_diagnostics_hydro(f,p)
      call calc_1d_diagnostics_hydro(f,p)
      call calc_0d_diagnostics_hydro(f,p)
!
!  store slices for output in wvid in run.f90
!  This must be done outside the diagnostics loop (accessed at different times).
!
      if (lvideo.and.lfirst) then
        if (ivid_divu/=0) call store_slices(p%divu,divu_xy,divu_xz,divu_yz,divu_xy2,divu_xy3,divu_xy4,divu_xz2,divu_r)
        if (ivid_oo  /=0) call store_slices(p%oo,oo_xy,oo_xz,oo_yz,oo_xy2,oo_xy3,oo_xy4,oo_xz2,oo_r)
        if (ivid_u2  /=0) call store_slices(p%u2,u2_xy,u2_xz,u2_yz,u2_xy2,u2_xy3,u2_xy4,u2_xz2,u2_r)
        if (ivid_o2  /=0) call store_slices(p%o2,o2_xy,o2_xz,o2_yz,o2_xy2,o2_xy3,o2_xy4,o2_xz2,o2_r)
        if (ivid_ou  /=0) call store_slices(p%ou,ou_xy,ou_xz,ou_yz,ou_xy2,ou_xy3,ou_xy4,ou_xz2,ou_r)
        if (ivid_Ma2 /=0) call store_slices(p%Ma2,mach_xy,mach_xz,mach_yz,mach_xy2,mach_xy3,mach_xy4,mach_xz2,mach_r)
        if (ivid_uu_sph/=0) call store_slices(p%uu_sph,uu_sph_xy,uu_sph_xz,uu_sph_yz,uu_sph_xy2, &
                                                       uu_sph_xy3,uu_sph_xy4,uu_sph_xz2,uu_sph_r)
        if (othresh_per_orms/=0) then
          call calc_othresh
          call vecout(41,trim(directory)//'/ovec',p%oo,othresh,novec)
        endif
      endif

    endsubroutine calc_diagnostics_hydro
!******************************************************************************
    subroutine df_diagnos_hydro(df,p)

      use Diagnostics, only: sum_mn_name

      type(pencil_case), intent(in) :: p
      real, dimension(:,:,:,:) :: df

      real, dimension (nx) :: uduu
      integer :: i

      if (idiag_uduum/=0) then
        uduu=0.
        do i = 1,3
          uduu=uduu+p%uu(:,i)*df(l1:l2,m,n,iux-1+i)
        enddo
        call sum_mn_name(p%rho*uduu,idiag_uduum)
      endif

    endsubroutine df_diagnos_hydro
!******************************************************************************
    subroutine time_integrals_hydro(f,p)
!
!  Calculate time_integrals within each pencil (as long as each
!  pencil case p still contains the current data). This routine
!  is now being called at the end of equ.
!
!  28-jun-07/axel+mreinhard: coded
!  24-jun-08/axel: moved call to this routine to the individual pde routines
!   1-jul-08/axel: moved this part to hydro
!  29-oct-20/hongzhe: coding for computing frequency-fft'ed u(x,y,z,omega_fourier)
!  19-may-21/axel: possibility of ltime_integrals_always=F to compute <u(t,x).u(t0,x)>
!  11-dec-21/hongzhe: uut and oot are now always in lab frame, but their update time
!                     is write to file for future shear-frame transformation.
!
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f
      intent(in) :: p
!
      real :: fact_cos,fact_sin,t_cor=0.
      logical :: lreset_vart=.false.
!
      fact_cos=cos(omega_fourier*t)
      fact_sin=sin(omega_fourier*t)
!
!  reset uut etc. if lreset_vart=.true.
!
      if (ltime_integrals_always .or. dtcor<=0.) then
        if (it==1) lreset_vart=.true.
      else
        if (t>t_vart) then
          lreset_vart=.true.
          t_cor=t
          call put_shared_variable('t_cor',t_cor,caller='time_integrals_hydro')
!
!  If uut and oot are updated, write t to file and advance t_var.
!  These only happen at the last point in the mn loop.
!
          if (llastpoint) then
            if (lroot) then
              open(1,file=trim(datadir)//'/tvart.dat',status='unknown',position='append')
              write(1,'(4f14.7)') t_cor
              close (1)
            endif
            t_vart=t_vart+dtcor
          endif
        endif
      endif
!
!  assign values to uut etc.
!
      if (ltime_integrals_always) then
        if (lreset_vart) then
          if (iuut/=0)  f(l1:l2,m,n,iuxt:iuzt)  =0.
          if (iuust/=0) f(l1:l2,m,n,iuxst:iuzst)=0.
          if (ioot/=0)  f(l1:l2,m,n,ioxt:iozt)  =0.
          if (ioost/=0) f(l1:l2,m,n,ioxst:iozst)=0.
        else
          if (iuut/=0)  f(l1:l2,m,n,iuxt:iuzt)  =f(l1:l2,m,n,iuxt:iuzt)  +dt*p%uu*fact_cos
          if (iuust/=0) f(l1:l2,m,n,iuxst:iuzst)=f(l1:l2,m,n,iuxst:iuzst)+dt*p%uu*fact_sin
          if (ioot/=0)  f(l1:l2,m,n,ioxt:iozt)  =f(l1:l2,m,n,ioxt:iozt)  +dt*p%oo*fact_cos
          if (ioost/=0) f(l1:l2,m,n,ioxst:iozst)=f(l1:l2,m,n,ioxst:iozst)+dt*p%oo*fact_sin
        endif
      else
        if (iuust/=0) f(l1:l2,m,n,iuxst:iuzst)=0.
        if (ioost/=0) f(l1:l2,m,n,ioxst:iozst)=0.
        if (lreset_vart) then
          if (lvart_in_shear_frame) &
              call fatal_error('time_integrals_hydro',&
                  'lvart_in_shear_frame is no longer supported in hydro. uut etc. are always &
                   &in lab frame. lshear_frame_correlation=T now transforms both uu and uut.')
          if (iuxt/=0)              f(l1:l2,m,n,iuxt:iuzt)  =f(l1:l2,m,n,iux:iuz)
          if (ioxt/=0 .and. iox/=0) f(l1:l2,m,n,ioxt:iozt)  =f(l1:l2,m,n,iox:ioz)
        endif
      endif
      lreset_vart=.false.
!
    endsubroutine time_integrals_hydro
!***********************************************************************
    subroutine hydro_after_boundary(f)
!
!  Calculate <rho*ux> and <rho*uy> when tau_damp_ruxm, tau_damp_ruym,
!  or tau_damp_ruzm are different from zero. Was used to remove net
!  momenta in any of the three directions. A better method is now
!  to set lremove_mean_momenta=T in the call to remove_mean_momenta.
!  Calculates <U>, when lcalc_uumean=.true.
!
!   9-nov-06/axel: adapted from calc_ltestfield_pars
!  31-jul-08/axel: Poincare force with O=(sinalp*cosot,sinalp*sinot,cosalp)
!  12-sep-13/MR  : use finalize_aver
!
      use Sub, only: finalize_aver

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (3,3) :: mat_cent1=0.,mat_cent2=0.,mat_cent3=0.
      real, dimension (3) :: OO, dOO, uum0
      real :: c,s,sinalp,cosalp,OO2,alpha_precession_rad
      integer :: i,j,l
!
      intent(inout) :: f
!
!  possibility of setting interior boundary conditions
!
      if (lhydro_bc_interior) call interior_bc_hydro(f)
!
      call calc_means_hydro(f)
!
!  calculate precession matrices, assume that alpha_precession is given
!  in degrees.
!
      if (lprecession) then
        c=cos(omega_precession*t)
        s=sin(omega_precession*t)
        alpha_precession_rad=alpha_precession
        cosalp=cos(alpha_precession_rad)
        sinalp=sin(alpha_precession_rad)
!
!  Components of Omega vector
!
        OO(1)=Omega*sinalp*c
        OO(2)=Omega*sinalp*s
        OO(3)=Omega*cosalp
!
!  Components of time derivative of Omega vector
!
        dOO(1)=-Omega*sinalp*s
        dOO(2)=+Omega*sinalp*c
        dOO(3)=0.
!
!  Coriolis matrix
!
        mat_cori(1,2)=+2.*OO(3); mat_cori(2,1)=-2.*OO(3)
        mat_cori(2,3)=+2.*OO(1); mat_cori(3,2)=-2.*OO(1)
        mat_cori(3,1)=+2.*OO(2); mat_cori(1,3)=-2.*OO(2)
!
!  1st centrifugal matrix
!
        do i=1,3
        do j=1,3
          mat_cent1(i,j)=-OO(i)*OO(j)
        enddo
        enddo
!
!  2nd centrifugal matrix
!
        OO2=OO(1)+OO(2)+OO(3)
        do j=1,3
          mat_cent2(j,j)=OO2
        enddo
!
!  3rd centrifugal matrix
!
        mat_cent3(1,2)=+dOO(3); mat_cent3(2,1)=-dOO(3)
        mat_cent3(2,3)=+dOO(1); mat_cent3(3,2)=-dOO(1)
        mat_cent3(3,1)=+dOO(2); mat_cent3(1,3)=-dOO(2)
!
!  adding the three centrifugal matrixes together
!
        mat_cent=mat_cent1+mat_cent2+mat_cent3
      endif
!
!  calculate precession matrices
!
!     if (lprecession) then
!       c=cos(omega_precession*t)
!       s=sin(omega_precession*t)
!       mat_cori1(2,3)=+1.
!       mat_cori1(3,2)=-1.
!
!
!       mat_cori2(1,2)=+c
!       mat_cori2(1,3)=-s
!       mat_cori2(2,1)=-c
!       mat_cori2(3,1)=+s
!       mat_cent1(2,2)=+1.
!       mat_cent1(3,3)=+1.
!       mat_cent2(1,1)=+1.
!       mat_cent2(2,2)=+c**2
!       mat_cent2(3,3)=+s**2
!       mat_cent2(2,3)=-2.*s*c
!       mat_cent2(3,2)=-2.*s*c
!       mat_cent3(1,2)=-s+c
!       mat_cent3(2,1)=-s-c
!       mat_cent3(1,3)=-s-c
!       mat_cent3(3,1)=+s-c
!       mat_cori=2.*(omega_precession*mat_cori1+Omega*mat_cori2)
!       mat_cent=omega_precession**2*mat_cent1+Omega**2*mat_cent2 &
!         +2.*omega_precession*Omega*mat_cent3
!     endif
!
!
!  Remove mean flow (z average).
!
      if (lrmv) then
        if (lremove_uumeanxy) then
!
          do j=1,3
            do n=1,mz
              f(:,:,n,iuu+j-1) = f(:,:,n,iuu+j-1)-uumxy(:,:,j)
            enddo
          enddo
        endif
!
!  Remove mean flow (xy average).
!
        if (lremove_uumeanz) then
          do j=1,3
            do n=1,mz
              f(:,:,n,iuu+j-1) = f(:,:,n,iuu+j-1)-uumz(n,j)
            enddo
          enddo
        elseif (lremove_uumeanz_horizontal) then
!
!  Remove only xy-averaged horizontal flows
!
          do j=1,2
            do n=1,mz
              f(:,:,n,iuu+j-1) = f(:,:,n,iuu+j-1)-uumz(n,j)
            enddo
          enddo
        endif
!
!  Remove mean flow (xz average).
!
        if (lremove_uumeany) then
!
!  uum0 is formed to avoid double substraction of <uu>_xyz.
!
          if (lremove_uumeanz.or.lremove_uumeanz_horizontal) then
            uum0=sum(uumy(m1:m2,:),1)/nygrid
            call finalize_aver(nprocy,2,uum0)
            if (lremove_uumeanz_horizontal) uum0(3)=0.
          endif

          do j=1,3
            do m=1,my
              f(:,m,:,iuu+j-1) = f(:,m,:,iuu+j-1)-uumy(m,j)
            enddo 
            if (lremove_uumeanz.or.lremove_uumeanz_horizontal) &
              f(:,:,:,iuu+j-1)=f(:,:,:,iuu+j-1)+uum0(j)  ! compensation as uum0 is already substracted once
          enddo

        endif
!
!  Remove mean flow (yz average).
!
        if (lremove_uumeanx) then
!
!  uum0 is formed to avoid double substraction of <uu>_xyz.
!
          if (lremove_uumeanz.or.lremove_uumeanz_horizontal.or.lremove_uumeany) then
            uum0=sum(uumx(l1:l2,:),1)/nxgrid
            call finalize_aver(nprocx,1,uum0)
            if (lremove_uumeanz_horizontal) uum0(3)=0.
          endif

          do j=1,3
            do l=1,mx
              f(l,:,:,iuu+j-1) = f(l,:,:,iuu+j-1)-uumx(l,j)
            enddo
            if (lremove_uumeanz.or.lremove_uumeanz_horizontal.or.lremove_uumeany) &
              f(:,:,:,iuu+j-1)=f(:,:,:,iuu+j-1)+uum0(j)  ! compensation as uum0 is already substracted once
          enddo
        endif
      endif
!
!  Compute fluctuating velocity and put in an auxilliary array
!
      if (luu_fluc_as_aux) then
        if (lcalc_uumeanz) then
          do j=1,3
            do n=1,mz
              f(l1:l2,m1:m2,n,iuu_fluc+j-1) = f(l1:l2,m1:m2,n,iuu+j-1) - uumz(n,j)
            enddo
          enddo
        endif
!
        if (lcalc_uumeanxy) then
          do j=1,3
            do n=1,mz
              f(l1:l2,m1:m2,n,iuu_fluc+j-1) = f(l1:l2,m1:m2,n,iuu+j-1) - uumxy(l1:l2-nghost,m1:m2-nghost,j)
            enddo
          enddo
        endif
      endif

    endsubroutine hydro_after_boundary
!***********************************************************************
    subroutine set_border_hydro(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the uu variable.
!
!  28-jul-06/wlad: coded
!
      use BorderProfiles,  only: border_driving,set_border_initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: f_target
      integer :: j,ju
!
      do j=1,3

        select case (borderuu(j))
!
        case ('zero','0')
          f_target(:,j)=0.
!
        case ('constant')
          f_target(:,j) = uu_const(j)
!
        case ('initial-condition')
          ju=j+iuu-1
          call set_border_initcond(f,ju,f_target(:,j))
!
        case ('nothing')
          if (lroot.and.ip<=5) &
               print*,"set_border_hydro: borderuu='nothing'"
!
        case default
          write(unit=errormsg,fmt=*) &
               'set_border_hydro: No such value for borderuu: ', &
               trim(borderuu(j))
          call fatal_error('set_border_hydro',errormsg)
!
        endselect

        if (borderuu(j)/='nothing') then
          ju=j+iuu-1
          call border_driving(f,df,p,f_target(:,j),ju)
        endif
!
      enddo
!
    endsubroutine set_border_hydro
!***********************************************************************
    subroutine calc_othresh
!
!  calculate othresh from orms, give warnings if there are problems
!
!  24-nov-03/axel: adapted from calc_bthresh
!
!  if nvec exceeds novecmax (=1/4) of points per processor, then begin to
!  increase scaling factor on othresh. These settings will stay in place
!  until the next restart
!
      if (novec>novecmax.and.lfirstpoint) then
        print*,'calc_othresh: processor ',iproc_world,': othresh_scl,novec,novecmax=', &
                                          othresh_scl,novec,novecmax
        othresh_scl=othresh_scl*1.2
      endif
!
!  calculate othresh as a certain fraction of orms
!
      othresh=othresh_scl*othresh_per_orms*orms  !!!MR: orms not yet valid when used inside duu_dt
!
    endsubroutine calc_othresh
!***********************************************************************
    subroutine precession(df,p)
!
!  precession terms
!
!  19-jan-07/axel: added terms derived by Gailitis
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: i,j,k
!
!  info about precession term
!
      if (headtt) then
        print*, 'precession: omega_precession=', omega_precession
      endif
!
!  matrix multiply
!
      do j=1,3
      do i=1,3
        k=iuu-1+i
        df(l1:l2,m,n,k)=df(l1:l2,m,n,k) &
          +mat_cent(i,j)*p%rr(:,j) &
          +mat_cori(i,j)*p%uu(:,j)
      enddo
      enddo
!
    endsubroutine precession
!***********************************************************************
   subroutine coriolis_cartesian(df,uu,velind)
!
!  Coriolis terms for Cartesian geometry.
!
!  30-oct-09/MR: outsourced, parameter velind added
!  15-feb-15/MR: calculation of Coriolis force of shear flow added
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension (nx,3),          intent(in)  :: uu
      integer,                         intent(in)  :: velind
!
!  velind is start index for velocity variable to which Coriolis force
!  corresponds
!  x,y,z components are referred to by velind, velind+1, velind+2
!
      real :: c2, s2
!
      if (Omega==0.) return
!
      if (theta==0) then
!
        if (lcoriolis_force) then
!
          if (headtt) print*,'coriolis_cartesian: add Coriolis force; Omega=',Omega
!
          c2=2*Omega
          df(l1:l2,m,n,velind  )=df(l1:l2,m,n,velind  )+c2*uu(:,2)
          df(l1:l2,m,n,velind+1)=df(l1:l2,m,n,velind+1)-c2*uu(:,1)
!
        endif
!
!  Add centrifugal force (doing this with periodic boundary
!  conditions in x and y would not be compatible, so it is
!  therefore usually ignored in those cases!)
!
        if (lcentrifugal_force) then
!
          if (headtt) print*,'coriolis_cartesian: add Centrifugal force; Omega=',Omega
          if (headtt) print*,'coriolis_cartesian: Centrifugal force amplitude; amp_centforce=',amp_centforce
!
          df(l1:l2,m,n,velind  )=df(l1:l2,m,n,velind  )+x(l1:l2)*amp_centforce*Omega**2
          df(l1:l2,m,n,velind+1)=df(l1:l2,m,n,velind+1)+y(  m  )*amp_centforce*Omega**2
!
        endif
!
      else

        if (phi/=0.) &
          call fatal_error("coriolis_cartesian:","not coded if Omega has a y component")
!
!  Add Coriolis force with an angle (defined such that theta=60,
!  for example, would correspond to 30 degrees latitude).
!  Omega=(-sin(theta), 0, cos(theta)).
!
        if (lcoriolis_force) then
!
          if (headtt) &
              print*,'coriolis_cartesian: Coriolis force; Omega, theta=', Omega, theta
!
!  Note the minus sign in front of the sin_theta term!
!
          c2= 2*Omega*cos(theta*pi/180.)
          s2=-2*Omega*sin(theta*pi/180.)
!
          df(l1:l2,m,n,velind  )=df(l1:l2,m,n,velind  )+c2*uu(:,2)
          df(l1:l2,m,n,velind+1)=df(l1:l2,m,n,velind+1)-c2*uu(:,1)+s2*uu(:,3)
          df(l1:l2,m,n,velind+2)=df(l1:l2,m,n,velind+2)           -s2*uu(:,2)
!
!  Add -2 \Omega \times U^shear, if requested.
!
          if (lshear_in_coriolis) then
            df(l1:l2,m,n,velind  )=df(l1:l2,m,n,velind  )+c2*Sshear*x(l1:l2)
            df(l1:l2,m,n,velind+2)=df(l1:l2,m,n,velind+2)-s2*Sshear*x(l1:l2)
          endif
        endif
!
      endif
!
   endsubroutine coriolis_cartesian
!***********************************************************************
   subroutine coriolis_cartesian_xaxis(df,uu,velind)
!
!  Coriolis force for a box where the rotation axis is in the x-direction
!  In this case the axis of the box represent to x=r, y=theta and  z=phi,
!  so that the z-averages in spherical coordinates corresponds to averages
!  in the phi direction in cartesian coordinates too.
!
!  Adapted from coriolis_cartesian.
!
!  09-aug-10/GG:
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension (nx,3),          intent(in)  :: uu
      integer,                         intent(in)  :: velind
!
      real :: c2, s2
!
      if (Omega==0.) return
!
      if (lcoriolis_force) then
!
        if (headtt) &
          print*,'coriolis_cartesian_xaxis: Coriolis force; Omega, theta=', Omega, theta
!
        c2= 2*Omega*cos(theta*pi/180.)
        s2=-2*Omega*sin(theta*pi/180.)
!
        df(l1:l2,m,n,velind  )=df(l1:l2,m,n,velind  )-s2*uu(:,3)
        df(l1:l2,m,n,velind+1)=df(l1:l2,m,n,velind+1)+c2*uu(:,3)
        df(l1:l2,m,n,velind+2)=df(l1:l2,m,n,velind+2)-c2*uu(:,2)+s2*uu(:,1)
!
      endif
!
   endsubroutine coriolis_cartesian_xaxis
!***********************************************************************
   subroutine coriolis_spherical(df,p)
!
!  coriolis_spherical terms using spherical polars
!
!  21-feb-07/axel+dhruba: coded
!  22-dec-15/MR: extended for situation with Omega along y axis (relevant for Yin-Yang grid).
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: c2,s2,Om2,cp2,cs2,ss2
!
      intent(in) :: p
      intent(inout) :: df
!
!  info about coriolis_spherical term
!
      if (headtt) print*, 'coriolis_spherical: Omega=', Omega
!
! Not yet coded for angular velocity at an angle with z or y axis.
!
      if (.not.(theta==.0.or.theta==90..and.phi==90.)) then
        if (lroot) print*, 'coriolis_spherical: Omega,theta,phi=', Omega,theta, phi
        call fatal_error("coriolis_spherical:","not coded if Omega isn't aligned with z or y axis")
      endif

      if (lcoriolis_force) then
        c2= 2*Omega*costh(m)
        s2=-2*Omega*sinth(m)
      endif

      if (theta==.0) then
!
!  Omega along z axis: In (r,theta,phi) coords, we have \vec{Omega}=Omega*(costh, -sinth, 0). Thus,
!
!                                ( costh)   (u1)           (+sinth*u3)
!  -2*\vec{Omega} x U = -2*Omega*(-sinth) X (u2) = 2*Omega*(+costh*u3)
!                                (   0  )   (u3)           (-costh*u2-sinth*u1)
!
!  With c2=2*Omega*costh and s2=-2*Omega*sinth we have then
!
!                       (-s2*u3)
!  -2*\vec{Omega} x U = (+c2*u3)
!                       (-c2*u2+s2*u1)
!
        if (lcoriolis_force) then
          if (r_omega /= 0.) then
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)- s2*p%uu(:,3)              *prof_om
            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+ c2*p%uu(:,3)              *prof_om
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-(c2*p%uu(:,2)-s2*p%uu(:,1))*prof_om
          else
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-s2*p%uu(:,3)
            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+c2*p%uu(:,3)
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-c2*p%uu(:,2)+s2*p%uu(:,1)
          endif
        endif
      else
!
!  Omega along y axis: In (r,theta,phi) coords, we have \vec{Omega}=Omega*(sinth*sinph, costh*sinph, cosph). Thus,
!
!                                (sinth*sinph)   (u1)           (-costh*sinph*u3+cosph*      u2)   (-cs2*u3 + cp2*u2)
!  -2*\vec{Omega} x U = -2*Omega*(costh*sinph) X (u2) = 2*Omega*(-cosph*      u1+sinth*sinph*u3) = (-cp2*u1 - ss2*u3)
!                                (      cosph)   (u3)           (-sinth*sinph*u2+costh*sinph*u1)   ( ss2*u2 + cs2*u1)

        cp2=2*Omega*cosph(n); cs2=c2*sinph(n); ss2=s2*sinph(n)

        if (lcoriolis_force) then
          if (r_omega /= 0.) then
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+(-cs2*p%uu(:,3) + cp2*p%uu(:,2))*prof_om
            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+(-cp2*p%uu(:,1) - ss2*p%uu(:,3))*prof_om
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+(+ss2*p%uu(:,2) + cs2*p%uu(:,1))*prof_om
          else
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux) - cs2*p%uu(:,3) + cp2*p%uu(:,2)
            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy) - cp2*p%uu(:,1) - ss2*p%uu(:,3)
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz) + ss2*p%uu(:,2) + cs2*p%uu(:,1)
          endif
        endif
      endif
!
!  Centrifugal force
!  The term added is F_{centrifugal} = - \Omega X \Omega X r
!
      if (lcentrifugal_force) then

        Om2=amp_centforce*Omega**2
        if (theta==.0) then
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-Om2*x(l1:l2)*sinth(m)
        else
          call fatal_error("coriolis_spherical:","Centrifugal force not coded if Omega isn't aligned with z axis.")
        endif

      endif
!
    endsubroutine coriolis_spherical
!***********************************************************************
    subroutine coriolis_spherical_del2p(f,p)
!
!  coriolis_spherical terms using spherical polars
!
!  21-feb-07/axel+dhruba: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
!  info about coriolis_spherical term
!
      if (headtt) print*, 'coriolis_spherical: Omega=', Omega
!
! Not yet coded for angular velocity at an angle with the z axis.
!
      if (theta/=0) then
        print*, 'coriolis_spherical: Omega=,theta=', Omega,theta
        call fatal_error("coriolis_spherical:","not coded if Omega is at an angle with the z axis. ")
      endif
!
!  In (r,theta,phi) coords, we have Omega=(costh, -sinth, 0). Thus,
!
!                    ( costh)   (u1)      (+sinth*u3)
!  -2*Omega x U = -2*(-sinth) X (u2) = 2*(+costh*u3)
!                    (   0  )   (u3)      (-costh*u2-sinth*u1)
!
!  With c2=2*Omega*costh and s2=-2*Omega*sinth we have then
!
!                (-s2*u3)
!  -2*Omega x U = (+c2*u3)
!                (-c2*u2+s2*u1)
!
!
!  Centrifugal force
!
      if (lcentrifugal_force) &
          call fatal_error("coriolis_spherical_del2p","Centrifugal force not "//&
          "implemented in spherical coordinates")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine coriolis_spherical_del2p
!***********************************************************************
    subroutine coriolis_cylindrical(df,p)
!
!  Coriolis terms using cylindrical coords
!  The formulation is the same as in cartesian, but it is better to
!  keep it here because precession is not implemented for
!  cylindrical coordinates.
!
!  19-sep-07/steveb: coded
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: c2, s2
!
!  info about coriolis_cylindrical term
!
      if (headtt) &
          print*, 'coriolis_cylindrical: Omega=', Omega
!
! Not yet coded for angular velocity at an angle with the z axis.
!
      if (theta/=0) then
         print*, 'coriolis_cylindrical: Omega=,theta=', Omega,theta
         call fatal_error("coriolis_cylindrical:","not coded if the angular velocity is at an angle to the z axis. ")
      endif
!
!  -2 Omega x u
!
      if (lcoriolis_force) then
        if (lOmega_cyl_xy) then
          c2= 2*Omega*cos(y(m))
          s2=-2*Omega*sin(y(m))
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-c2*p%uu(:,3)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-s2*p%uu(:,3)
          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+c2*p%uu(:,1)+s2*p%uu(:,2)
        else
          c2=2*Omega
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+c2*p%uu(:,2)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-c2*p%uu(:,1)
        endif
      endif
!
!  Centrifugal force
!
      if (lcentrifugal_force) &
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+x(l1:l2)*Omega**2
!
!  Note, there is no z-component
!
    endsubroutine coriolis_cylindrical
!***********************************************************************
    subroutine coriolis_cylindrical_del2p(f,p)
!
!  Coriolis terms using cylindrical coords
!  The formulation is the same as in cartesian, but it is better to
!  keep it here because precession is not implemented for
!  cylindrical coordinates.
!
!  19-sep-07/steveb: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
!  info about coriolis_cylindrical term
!
      if (headtt) &
          print*, 'coriolis_cylindrical: Omega=', Omega
!
! Not yet coded for angular velocity at an angle with the z axis.
!
      if (theta/=0) then
         print*, 'coriolis_cylindrical: Omega=,theta=', Omega,theta
         call fatal_error("coriolis_cylindrical:","not coded if the angular velocity is at an angle to the z axis. ")
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine coriolis_cylindrical_del2p
!***********************************************************************
    subroutine coriolis_xdep(df,p)
!
!  Coriolis terms in Cartesian coordinates with Omega depending
!  on x, i.e. Omega=Omega0*(-sin(k_x*x),0,cos(l_x*x)) with k_x=2pi/Lx.
!
!  28-may-09/PJK: coded
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: c1, c2
!
!  info about coriolis_cylindrical term
!
      if (headtt) &
          print*, 'coriolis_xdep: ampl_Omega=', ampl_Omega
!
!  -2 Omega x u
!
      c1=-2*ampl_Omega*sin(pi*((x(l1:l2))-x0)/Lx)
      c2= 2*ampl_Omega*cos(pi*((x(l1:l2))-x0)/Lx)
      df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)             +c2*p%uu(:,2)
      df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+c1*p%uu(:,3)-c2*p%uu(:,1)
      df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-c1*p%uu(:,2)
!
!  Centrifugal force not coded yet
!
    endsubroutine coriolis_xdep
!***********************************************************************
    subroutine udamping(f,df,p)
!
!  damping terms (artificial, but sometimes useful):
!
!  20-nov-04/axel: added cylindrical Couette flow
!
      use Diagnostics, only: sum_mn_name
      use Sub, only: step
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: pdamp,fint_work,fext_work
      real, dimension (nx,3) :: fint,fext
      real, save :: tau, fade_fact, last_t = -1.0
      integer :: i,j
!
!  warn about the damping term
!
        if (headtt .and. (dampu /= 0.) .and. (t < tdamp)) then
          if (ldamp_fade) then
            print*, 'udamping: Damping velocities until time ', tdamp
            print*, 'udamping: with a smooth fade starting at ', tfade_start
          else
            print*, 'udamping: Damping velocities constantly until time ', tdamp
          endif
        endif
!
!  1. damp motion during time interval 0<t<tdamp.
!  Damping coefficient is dampu (if >0) or |dampu|/dt (if dampu <0).
!  With ldamp_fade=T, damping coefficient is smoothly fading out
!
!
        if ((dampu /= 0.) .and. (t < tdamp)) then
!
          if (.not. ldamp_fade) then
            ! no fading => full damping:
            fade_fact = 1.
          elseif (t <= tfade_start) then
            ! before transition => full damping:
            fade_fact = 1.
          else
            ! inside transition => smooth fading:
            if (last_t /= t) then
              last_t = t
!
!  smoothly fade out damping according to the following
!  function of time:
!
!    ^
!    |
!  1 +**************
!    |              ****
!    |                  **
!    |                    *
!    |                     **
!    |                       ****
!  0 +-------------+-------------**********---> t
!    |             |             |
!    0            Tfade_start   Tdamp
!
!  For 0 < t < Tfade_start, full damping is applied.
!  In the interval Tfade_start < t < Tdamp, damping goes smoothly to zero
!  with continuous derivatives. (The default value for Tfade_start is Tdamp/2.)
!
              ! tau is a normalized t, the transition interval is [-0.5, 0.5]:
              tau = (t-tfade_start) / (tdamp-tfade_start) - 0.5
              if (tau <= -0.5) then
                fade_fact = 1.
              elseif (tau <= 0.5) then
                fade_fact = 0.5 - tau * (1.5 - 2.0*tau**2)
              else
                call fatal_error("udamping","tau is invalid (tau > 0.5).")
              endif
            endif
          endif
!
          if (dampu > 0.0) then
            ! absolute damping per time unit
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) &
                                    - fade_fact*dampu*f(l1:l2,m,n,iux:iuz)
          else
            ! dampu < 0: damping per time-step (dt is multiplied in timestep)
            df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) &
                                    + fade_fact*dampu/dt*f(l1:l2,m,n,iux:iuz)
          endif
        endif
!
!  2. damp motions for p%r_mn > rdampext or r_ext AND p%r_mn < rdampint or r_int
!
        if (dampuext > 0.0 .and. rdampext /= impossible) then
          ! outer damping profile
          pdamp = step(p%r_mn,rdampext,wdamp)
          do i=iux,iuz
            df(l1:l2,m,n,i) = df(l1:l2,m,n,i) - dampuext*pdamp*f(l1:l2,m,n,i)
          enddo
        endif
!
        if (dampuint > 0.0 .and. rdampint /= impossible) then
          ! inner damping profile
          pdamp = 1 - step(p%r_mn,rdampint,wdamp)
          do i=iux,iuz
            df(l1:l2,m,n,i) = df(l1:l2,m,n,i) - dampuint*pdamp*f(l1:l2,m,n,i)
          enddo
        endif
!
!  coupling the above internal and external rotation rates to lgravr is not
!  a good idea. So, because of that, spherical Couette flow has to be coded
!  separately.
!  ==> reconsider name <==
!  Allow now also for cylindical Couette flow (if lcylinder_in_a_box=T)
!
        if (lOmega_int) then
!
!  relax outer angular velocity to zero, and
!  calculate work done to sustain zero rotation on outer cylinder/sphere
!
!
          if (lcylinder_in_a_box) then
            pdamp = step(p%rcyl_mn,rdampext,wdamp) ! outer damping profile
          else
            pdamp = step(p%r_mn,rdampext,wdamp) ! outer damping profile
          endif
!
          do i=1,3
            j=iux-1+i
            fext(:,i)=-dampuext*pdamp*f(l1:l2,m,n,j)
            df(l1:l2,m,n,j)=df(l1:l2,m,n,j)+fext(:,i)
          enddo
          if (ldiagnos.and.idiag_fextm/=0) then
            fext_work=f(l1:l2,m,n,iux)*fext(:,1)&
                     +f(l1:l2,m,n,iuy)*fext(:,2)&
                     +f(l1:l2,m,n,iuz)*fext(:,3)
            call sum_mn_name(fext_work,idiag_fextm)
          endif
!
!  internal angular velocity, uref=(-y,x,0)*Omega_int, and
!  calculate work done to sustain uniform rotation on inner cylinder/sphere
!
          if (dampuint > 0.0) then
            if (lcylinder_in_a_box) then
              pdamp = 1 - step(p%rcyl_mn,rdampint,wdamp) ! inner damping profile
            else
              pdamp = 1 - step(p%r_mn,rdampint,wdamp) ! inner damping profile
            endif
            fint(:,1)=-dampuint*pdamp*(f(l1:l2,m,n,iux)+y(m)*Omega_int)
            fint(:,2)=-dampuint*pdamp*(f(l1:l2,m,n,iuy)-x(l1:l2)*Omega_int)
            fint(:,3)=-dampuint*pdamp*(f(l1:l2,m,n,iuz))
            df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+fint(:,1)
            df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+fint(:,2)
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+fint(:,3)
            if (ldiagnos.and.idiag_fintm/=0) then
              fint_work=f(l1:l2,m,n,iux)*fint(:,1)&
                       +f(l1:l2,m,n,iuy)*fint(:,2)&
                       +f(l1:l2,m,n,iuz)*fint(:,3)
              call sum_mn_name(fint_work,idiag_fintm)
            endif
          endif
        endif
!
    endsubroutine udamping
!***********************************************************************
    subroutine read_hydro_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=hydro_init_pars, IOSTAT=iostat)
!
    endsubroutine read_hydro_init_pars
!***********************************************************************
    subroutine write_hydro_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=hydro_init_pars)
!
    endsubroutine write_hydro_init_pars
!***********************************************************************
    subroutine read_hydro_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=hydro_run_pars, IOSTAT=iostat)
!
    endsubroutine read_hydro_run_pars
!***********************************************************************
    subroutine write_hydro_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=hydro_run_pars)
!
    endsubroutine write_hydro_run_pars
!***********************************************************************
    subroutine input_persistent_hydro(id,done)
!
!  Dummy. Reads the hydro persistent variables only in 'hydro_kinematic'.
!
      integer, optional :: id
      logical, optional :: done
!
      if (present (id)) call keep_compiler_quiet(id)
      if (present (done)) call keep_compiler_quiet(done)
!
    endsubroutine input_persistent_hydro
!***********************************************************************
    logical function output_persistent_hydro()
!
      output_persistent_hydro = .false.
!
    endfunction output_persistent_hydro
!***********************************************************************
    subroutine rprint_hydro(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
      use General, only: itoa
      use Mpicomm, only: stop_it
!
      integer :: k
      character (len=intlen) :: smode
      integer :: iname,inamez,inamey,inamex,ixy,ixz,irz,inamer,iname_half,inamev,idum
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_u2tm=0
        idiag_uotm=0
        idiag_outm=0
        idiag_fkinzm=0
        idiag_u2m=0
        idiag_u2sphm=0
        idiag_um2=0
        idiag_uxpt=0
        idiag_uypt=0
        idiag_uzpt=0
        idiag_uxp2=0
        idiag_uyp2=0
        idiag_uzp2=0
        idiag_uxuypt=0
        idiag_uyuzpt=0
        idiag_uzuxpt=0
        idiag_urms=0
        idiag_durms=0
        idiag_urmsx=0
        idiag_urmsz=0
        idiag_umax=0
        idiag_umin=0
        idiag_uxrms=0
        idiag_uyrms=0
        idiag_uzrms=0
        idiag_uzrmaxs=0
        idiag_uxmin=0
        idiag_uymin=0
        idiag_uzmin=0
        idiag_uxmax=0
        idiag_uymax=0
        idiag_uzmax=0
        idiag_uxm=0
        idiag_uym=0
        idiag_uzm=0
        idiag_ux2m=0
        idiag_uy2m=0
        idiag_uz2m=0
        idiag_ux3m=0
        idiag_uy3m=0
        idiag_uz3m=0
        idiag_ux4m=0
        idiag_uy4m=0
        idiag_uz4m=0
        idiag_uxuy2m=0
        idiag_uyuz2m=0
        idiag_uzux2m=0
        idiag_ux2ccm=0
        idiag_ux2ssm=0
        idiag_uy2ccm=0
        idiag_uy2ssm=0
        idiag_uxuycsm=0
        idiag_rux2m=0
        idiag_ruy2m=0
        idiag_ruz2m=0
        idiag_uduum=0
        idiag_uxmx=0
        idiag_uymx=0
        idiag_uzmx=0
        idiag_ux2mx=0
        idiag_uy2mx=0
        idiag_uz2mx=0
        idiag_ox2mx=0
        idiag_oy2mx=0
        idiag_oz2mx=0
        idiag_ux2my=0
        idiag_uy2my=0
        idiag_uz2my=0
        idiag_ux2mz=0
        idiag_uy2mz=0
        idiag_uz2mz=0
        idiag_uz2upmz=0
        idiag_uz2downmz=0
        idiag_ruxmx=0
        idiag_ruymx=0
        idiag_ruzmx=0
        idiag_rux2mx = 0
        idiag_ruy2mx = 0
        idiag_ruz2mx = 0
        idiag_ruxuymx = 0
        idiag_ruxuzmx = 0
        idiag_ruyuzmx = 0
        idiag_u2mz=0
        idiag_o2mz=0
        idiag_ruxmz=0
        idiag_ruymz=0
        idiag_ruzmz=0
        idiag_rux2mz=0
        idiag_ruy2mz=0
        idiag_ruz2mz=0
        idiag_uxmz=0
        idiag_uymz=0
        idiag_uzmz=0
        idiag_uzupmz=0
        idiag_uzdownmz=0
        idiag_ffdownmz=0
        idiag_ruzupmz=0
        idiag_ruzdownmz=0
        idiag_divumz=0
        idiag_divu2mz=0
        idiag_uzdivumz=0
        idiag_divrhourms=0
        idiag_divrhoumax=0
        idiag_oxmz=0
        idiag_oymz=0
        idiag_ozmz=0
        idiag_ox2mz=0
        idiag_oy2mz=0
        idiag_oz2mz=0
        idiag_uxuym=0
        idiag_uxuzm=0
        idiag_uyuzm=0
        idiag_uxuymx=0
        idiag_uxuzmx=0
        idiag_uyuzmx=0
        idiag_uxuymz=0
        idiag_uxuzmz=0
        idiag_uyuzmz=0
        idiag_uxuymz=0
        idiag_Rxymz=0
        idiag_Rxyupmz=0
        idiag_Rxydownmz=0
        idiag_Rxzmz=0
        idiag_Rxzupmz=0
        idiag_Rxzdownmz=0
        idiag_Ryzmz=0
        idiag_Ryzupmz=0
        idiag_Ryzdownmz=0
        idiag_oxuxxmz=0
        idiag_oyuxymz=0
        idiag_oxuyxmz=0
        idiag_oyuyymz=0
        idiag_oxuzxmz=0
        idiag_oyuzymz=0
        idiag_uyxuzxmz=0
        idiag_uyyuzymz=0
        idiag_uyzuzzmz=0
        idiag_umx=0
        idiag_umy=0
        idiag_umz=0
        idiag_omumz=0
        idiag_umamz=0
        idiag_umbmz=0
        idiag_umxbmz=0
        idiag_divum=0
        idiag_rdivum=0
        idiag_divu2m=0
        idiag_gdivu2m=0
        idiag_u3u21m=0
        idiag_u1u32m=0
        idiag_u2u13m=0
        idiag_u2u31m=0
        idiag_u3u12m=0
        idiag_u1u23m=0
        idiag_u3u21mz=0
        idiag_u1u32mz=0
        idiag_u2u13mz=0
        idiag_u2u31mz=0
        idiag_u3u12mz=0
        idiag_u1u23mz=0
        idiag_acczmz=0
        idiag_acczupmz=0
        idiag_acczdownmz=0
        idiag_accpowzmz=0
        idiag_accpowzupmz=0
        idiag_accpowzdownmz=0
        idiag_totalforcezmz=0
        idiag_totalforcezupmz=0
        idiag_totalforcezdownmz=0
        idiag_urmphi=0
        idiag_ursphmphi=0
        idiag_uthmphi=0
        idiag_rursphmphi=0
        idiag_ruthmphi=0
        idiag_upmphi=0
        idiag_uzmphi=0
        idiag_rurmphi=0
        idiag_rupmphi=0
        idiag_ruzmphi=0
        idiag_u2mphi=0
        idiag_ur2mphi=0
        idiag_up2mphi=0
        idiag_uz2mphi=0
        idiag_urupmphi=0
        idiag_uruzmphi=0
        idiag_upuzmphi=0
        idiag_rurupmphi=0
        idiag_ruruzmphi=0
        idiag_rupuzmphi=0
        idiag_fkinrsphmphi=0
        idiag_uxmy=0
        idiag_uymy=0
        idiag_uzmy=0
        idiag_uxuymy=0
        idiag_uxuzmy=0
        idiag_uyuzmy=0
        idiag_u2mr=0
        idiag_urmr=0
        idiag_upmr=0
        idiag_uzmr=0
        idiag_uxfampm=0
        idiag_uyfampm=0
        idiag_uzfampm=0
        idiag_uxmxz=0
        idiag_uymxz=0
        idiag_uzmxz=0
        idiag_ux2mxz=0
        idiag_uy2mxz=0
        idiag_uz2mxz=0
        idiag_uxuymxz=0
        idiag_uxuzmxz=0
        idiag_uyuzmxz=0
        idiag_uxmxy=0
        idiag_uymxy=0
        idiag_uzmxy=0
        idiag_uxupmxy=0
        idiag_uxdownmxy=0
        idiag_ruxupmxy=0
        idiag_ruxdownmxy=0
        idiag_ux2upmxy=0
        idiag_ux2downmxy=0
        idiag_ffdownmxy=0
        idiag_uxuymxy=0
        idiag_uxuzmxy=0
        idiag_uyuzmxy=0
        idiag_Rxymxy=0
        idiag_Rxyupmxy=0
        idiag_Rxydownmxy=0
        idiag_Rxzmxy=0
        idiag_Rxzupmxy=0
        idiag_Rxzdownmxy=0
        idiag_Ryzmxy=0
        idiag_Ryzupmxy=0
        idiag_Ryzdownmxy=0
        idiag_oxmxy=0
        idiag_oymxy=0
        idiag_ozmxy=0
        idiag_pvzmxy=0
        idiag_ruxmxy=0
        idiag_ruymxy=0
        idiag_ruzmxy=0
        idiag_ux2mxy=0
        idiag_uy2mxy=0
        idiag_uz2mxy=0
        idiag_rux2mxy=0
        idiag_ruy2mxy=0
        idiag_ruz2mxy=0
        idiag_ruxuymxy=0
        idiag_ruxuzmxy=0
        idiag_ruyuzmxy=0
        idiag_ruxm=0
        idiag_ruym=0
        idiag_ruzm=0
        idiag_ruxtot=0
        idiag_rlxm=0
        idiag_rlym=0
        idiag_rlzm=0
        idiag_rlx2m=0
        idiag_rly2m=0
        idiag_rlz2m=0
        idiag_tot_ang_mom=0
        idiag_rumax=0
        idiag_dtu=0
        idiag_oum=0
        idiag_oxum=0
        idiag_ourms=0
        idiag_oxurms=0
        idiag_ou_int=0
        idiag_fum=0
        idiag_odel2um=0
        idiag_o2m=0
        idiag_o2u2m=0
        idiag_o2sphm=0
        idiag_orms=0
        idiag_omax=0
        idiag_ox2m=0
        idiag_oy2m=0
        idiag_oz2m=0
        idiag_ox3m=0
        idiag_oy3m=0
        idiag_oz3m=0
        idiag_ox4m=0
        idiag_oy4m=0
        idiag_oz4m=0
        idiag_oxm=0
        idiag_oym=0
        idiag_ozm=0
        idiag_oxuzxm=0
        idiag_oyuzym=0
        idiag_oxoym=0
        idiag_oxozm=0
        idiag_oyozm=0
        idiag_qfm=0
        idiag_q2m=0
        idiag_qrms=0
        idiag_qmax=0
        idiag_qom=0
        idiag_quxom=0
        idiag_pvzm=0
        idiag_oumx=0
        idiag_oumy=0
        idiag_oumz=0
        idiag_oumxy=0
        idiag_oumxz=0
        idiag_oumphi=0
        idiag_ozmphi=0
        idiag_ormr=0
        idiag_opmr=0
        idiag_ozmr=0
        idiag_Marms=0
        idiag_Mamax=0
        idiag_fintm=0
        idiag_fextm=0
        idiag_divuHrms=0
        idiag_uxxrms=0
        idiag_uyyrms=0
        idiag_uxzrms=0
        idiag_uyzrms=0
        idiag_uzyrms=0
        idiag_duxdzma=0
        idiag_duydzma=0
        idiag_EEK=0
        idiag_ekin=0
        idiag_totangmom=0
        idiag_ekintot=0
        idiag_ekinmx=0
        idiag_ekinmz=0
        idiag_fmasszmz=0
        idiag_fkinzmz=0
        idiag_fkinzdownmz=0
        idiag_fkinzupmz=0
        idiag_fkinxmx=0
        idiag_fkinxmxy=0
        idiag_fkinymxy=0
        idiag_fkinxupmxy=0
        idiag_fkinxdownmxy=0
        idiag_uguxmxy=0
        idiag_uguymxy=0
        idiag_uguzmxy=0
        idiag_ruxuym=0
        idiag_ruxuzm=0
        idiag_ruyuzm=0
        idiag_ruxuymz=0
        idiag_ruxuzmz=0
        idiag_ruyuzmz=0
        idiag_ruxuy2mz=0
        idiag_ruxuz2mz=0
        idiag_ruyuz2mz=0
        idiag_uguxm=0
        idiag_uguym=0
        idiag_uguzm=0
        idiag_ugu2m=0
        idiag_dudx=0
        idiag_ugurmsx=0
        idiag_uguxmx=0
        idiag_uguymx=0
        idiag_uguzmx=0
        idiag_uguxmy=0
        idiag_uguymy=0
        idiag_uguzmy=0
        idiag_uguxmz=0
        idiag_uguymz=0
        idiag_uguzmz=0
        idiag_Remz=0
        idiag_oguxmz=0
        idiag_oguymz=0
        idiag_oguzmz=0
        idiag_ogux2mz=0
        idiag_oguy2mz=0
        idiag_oguz2mz=0
        idiag_oxdivumz=0
        idiag_oydivumz=0
        idiag_ozdivumz=0
        idiag_oxdivu2mz=0
        idiag_oydivu2mz=0
        idiag_ozdivu2mz=0
        idiag_uxglnrym=0
        idiag_uyglnrxm=0
        idiag_uzdivum=0
        idiag_uxuydivum=0
        idiag_urmsh=0;idiag_urmsn=0;idiag_urmss=0
        idiag_ormsh=0;idiag_ormsn=0;idiag_ormss=0
        idiag_oumh=0;idiag_oumn=0;idiag_oums=0
        idiag_udpxxm=0;idiag_udpyym=0;idiag_udpzzm=0
        idiag_udpxym=0;idiag_udpyzm=0;idiag_udpxzm=0
        idiag_taufmin=0
        idiag_dtF=0
        idiag_nshift=0
        ivid_oo=0; ivid_o2=0; ivid_ou=0; ivid_divu=0; ivid_u2=0; ivid_Ma2=0; ivid_uu_sph=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_hydro: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'EEK',idiag_EEK)
        call parse_name(iname,cname(iname),cform(iname),'ekin',idiag_ekin)
        call parse_name(iname,cname(iname),cform(iname),'ekintot',idiag_ekintot)
        call parse_name(iname,cname(iname),cform(iname),'u2tm',idiag_u2tm)
        call parse_name(iname,cname(iname),cform(iname),'uotm',idiag_uotm)
        call parse_name(iname,cname(iname),cform(iname),'outm',idiag_outm)
        call parse_name(iname,cname(iname),cform(iname),'fkinzm',idiag_fkinzm)
        call parse_name(iname,cname(iname),cform(iname),'u2m',idiag_u2m)
        call parse_name(iname,cname(iname),cform(iname),'u2sphm',idiag_u2sphm)
        call parse_name(iname,cname(iname),cform(iname),'um2',idiag_um2)
        call parse_name(iname,cname(iname),cform(iname),'odel2um',idiag_odel2um)
        call parse_name(iname,cname(iname),cform(iname),'o2m',idiag_o2m)
        call parse_name(iname,cname(iname),cform(iname),'o2u2m',idiag_o2u2m)
        call parse_name(iname,cname(iname),cform(iname),'o2sphm',idiag_o2sphm)
        call parse_name(iname,cname(iname),cform(iname),'oum',idiag_oum)
        call parse_name(iname,cname(iname),cform(iname),'oxum',idiag_oxum)
        call parse_name(iname,cname(iname),cform(iname),'ourms',idiag_ourms)
        call parse_name(iname,cname(iname),cform(iname),'oxurms',idiag_oxurms)
        call parse_name(iname,cname(iname),cform(iname),'ou_int',idiag_ou_int)
        call parse_name(iname,cname(iname),cform(iname),'fum',idiag_fum)
        call parse_name(iname,cname(iname),cform(iname),'oumn',idiag_oumn)
        call parse_name(iname,cname(iname),cform(iname),'oums',idiag_oums)
        call parse_name(iname,cname(iname),cform(iname),'dtu',idiag_dtu)
        call parse_name(iname,cname(iname),cform(iname),'urms',idiag_urms)
        call parse_name(iname,cname(iname),cform(iname),'durms',idiag_durms)
        call parse_name(iname,cname(iname),cform(iname),'urmsx',idiag_urmsx)
        call parse_name(iname,cname(iname),cform(iname),'urmsz',idiag_urmsz)
        call parse_name(iname,cname(iname),cform(iname),'urmsn',idiag_urmsn)
        call parse_name(iname,cname(iname),cform(iname),'urmss',idiag_urmss)
        call parse_name(iname,cname(iname),cform(iname),'umax',idiag_umax)
        call parse_name(iname,cname(iname),cform(iname),'umin',idiag_umin)
        call parse_name(iname,cname(iname),cform(iname),'uxmin',idiag_uxmin)
        call parse_name(iname,cname(iname),cform(iname),'uymin',idiag_uymin)
        call parse_name(iname,cname(iname),cform(iname),'uzmin',idiag_uzmin)
        call parse_name(iname,cname(iname),cform(iname),'uxmax',idiag_uxmax)
        call parse_name(iname,cname(iname),cform(iname),'uymax',idiag_uymax)
        call parse_name(iname,cname(iname),cform(iname),'uzmax',idiag_uzmax)
        call parse_name(iname,cname(iname),cform(iname),'uxrms',idiag_uxrms)
        call parse_name(iname,cname(iname),cform(iname),'uyrms',idiag_uyrms)
        call parse_name(iname,cname(iname),cform(iname),'uzrms',idiag_uzrms)
        call parse_name(iname,cname(iname),cform(iname),'uzrmaxs',idiag_uzrmaxs)
        call parse_name(iname,cname(iname),cform(iname),'uxm',idiag_uxm)
        call parse_name(iname,cname(iname),cform(iname),'uym',idiag_uym)
        call parse_name(iname,cname(iname),cform(iname),'uzm',idiag_uzm)
        call parse_name(iname,cname(iname),cform(iname),'ux2m',idiag_ux2m)
        call parse_name(iname,cname(iname),cform(iname),'uy2m',idiag_uy2m)
        call parse_name(iname,cname(iname),cform(iname),'uz2m',idiag_uz2m)
        call parse_name(iname,cname(iname),cform(iname),'ux3m',idiag_ux3m)
        call parse_name(iname,cname(iname),cform(iname),'uy3m',idiag_uy3m)
        call parse_name(iname,cname(iname),cform(iname),'uz3m',idiag_uz3m)
        call parse_name(iname,cname(iname),cform(iname),'ux4m',idiag_ux4m)
        call parse_name(iname,cname(iname),cform(iname),'uy4m',idiag_uy4m)
        call parse_name(iname,cname(iname),cform(iname),'uz4m',idiag_uz4m)
        call parse_name(iname,cname(iname),cform(iname),'uxuy2m',idiag_uxuy2m)
        call parse_name(iname,cname(iname),cform(iname),'uyuz2m',idiag_uyuz2m)
        call parse_name(iname,cname(iname),cform(iname),'uzux2m',idiag_uzux2m)
        call parse_name(iname,cname(iname),cform(iname),'ux2ccm',idiag_ux2ccm)
        call parse_name(iname,cname(iname),cform(iname),'ux2ssm',idiag_ux2ssm)
        call parse_name(iname,cname(iname),cform(iname),'uy2ccm',idiag_uy2ccm)
        call parse_name(iname,cname(iname),cform(iname),'uy2ssm',idiag_uy2ssm)
        call parse_name(iname,cname(iname),cform(iname),'uxuycsm',idiag_uxuycsm)
        call parse_name(iname,cname(iname),cform(iname),'rux2m',idiag_rux2m)
        call parse_name(iname,cname(iname),cform(iname),'ruy2m',idiag_ruy2m)
        call parse_name(iname,cname(iname),cform(iname),'ruz2m',idiag_ruz2m)
        call parse_name(iname,cname(iname),cform(iname),'uxuym',idiag_uxuym)
        call parse_name(iname,cname(iname),cform(iname),'uxuzm',idiag_uxuzm)
        call parse_name(iname,cname(iname),cform(iname),'uyuzm',idiag_uyuzm)
        call parse_name(iname,cname(iname),cform(iname),'ruxuym',idiag_ruxuym)
        call parse_name(iname,cname(iname),cform(iname),'ruxuzm',idiag_ruxuzm)
        call parse_name(iname,cname(iname),cform(iname),'ruyuzm',idiag_ruyuzm)
        call parse_name(iname,cname(iname),cform(iname),'ox2m',idiag_ox2m)
        call parse_name(iname,cname(iname),cform(iname),'oy2m',idiag_oy2m)
        call parse_name(iname,cname(iname),cform(iname),'oz2m',idiag_oz2m)
        call parse_name(iname,cname(iname),cform(iname),'ox3m',idiag_ox3m)
        call parse_name(iname,cname(iname),cform(iname),'oy3m',idiag_oy3m)
        call parse_name(iname,cname(iname),cform(iname),'oz3m',idiag_oz3m)
        call parse_name(iname,cname(iname),cform(iname),'ox4m',idiag_ox4m)
        call parse_name(iname,cname(iname),cform(iname),'oy4m',idiag_oy4m)
        call parse_name(iname,cname(iname),cform(iname),'oz4m',idiag_oz4m)
        call parse_name(iname,cname(iname),cform(iname),'oxm',idiag_oxm)
        call parse_name(iname,cname(iname),cform(iname),'oym',idiag_oym)
        call parse_name(iname,cname(iname),cform(iname),'ozm',idiag_ozm)
        call parse_name(iname,cname(iname),cform(iname),'oxuzxm',idiag_oxuzxm)
        call parse_name(iname,cname(iname),cform(iname),'oyuzym',idiag_oyuzym)
        call parse_name(iname,cname(iname),cform(iname),'oxoym',idiag_oxoym)
        call parse_name(iname,cname(iname),cform(iname),'oxozm',idiag_oxozm)
        call parse_name(iname,cname(iname),cform(iname),'oyozm',idiag_oyozm)
        call parse_name(iname,cname(iname),cform(iname),'pvzm',idiag_pvzm)
        call parse_name(iname,cname(iname),cform(iname),'orms',idiag_orms)
        call parse_name(iname,cname(iname),cform(iname),'qfm',idiag_qfm)
        call parse_name(iname,cname(iname),cform(iname),'q2m',idiag_q2m)
        call parse_name(iname,cname(iname),cform(iname),'qrms',idiag_qrms)
        call parse_name(iname,cname(iname),cform(iname),'qmax',idiag_qmax)
        call parse_name(iname,cname(iname),cform(iname),'qom',idiag_qom)
        call parse_name(iname,cname(iname),cform(iname),'quxom',idiag_quxom)
        call parse_name(iname,cname(iname),cform(iname),'ormsn',idiag_ormsn)
        call parse_name(iname,cname(iname),cform(iname),'ormss',idiag_ormss)
        call parse_name(iname,cname(iname),cform(iname),'omax',idiag_omax)
        call parse_name(iname,cname(iname),cform(iname),'ruxm',idiag_ruxm)
        call parse_name(iname,cname(iname),cform(iname),'ruym',idiag_ruym)
        call parse_name(iname,cname(iname),cform(iname),'ruzm',idiag_ruzm)
        call parse_name(iname,cname(iname),cform(iname),'ruxtot',idiag_ruxtot)
        call parse_name(iname,cname(iname),cform(iname),'rlxm',idiag_rlxm)
        call parse_name(iname,cname(iname),cform(iname),'rlym',idiag_rlym)
        call parse_name(iname,cname(iname),cform(iname),'rlzm',idiag_rlzm)
        call parse_name(iname,cname(iname),cform(iname),'rlx2m',idiag_rlx2m)
        call parse_name(iname,cname(iname),cform(iname),'rly2m',idiag_rly2m)
        call parse_name(iname,cname(iname),cform(iname),'rlz2m',idiag_rlz2m)
        call parse_name(iname,cname(iname),cform(iname),'tot_ang_mom',idiag_tot_ang_mom)
        call parse_name(iname,cname(iname),cform(iname),'rumax',idiag_rumax)
        call parse_name(iname,cname(iname),cform(iname),'umx',idiag_umx)
        call parse_name(iname,cname(iname),cform(iname),'umy',idiag_umy)
        call parse_name(iname,cname(iname),cform(iname),'umz',idiag_umz)
        call parse_name(iname,cname(iname),cform(iname),'omumz',idiag_omumz)
        call parse_name(iname,cname(iname),cform(iname),'umamz',idiag_umamz)
        call parse_name(iname,cname(iname),cform(iname),'umbmz',idiag_umbmz)
        call parse_name(iname,cname(iname),cform(iname),'umxbmz',idiag_umxbmz)
        call parse_name(iname,cname(iname),cform(iname),'Marms',idiag_Marms)
        call parse_name(iname,cname(iname),cform(iname),'Mamax',idiag_Mamax)
        call parse_name(iname,cname(iname),cform(iname),'divum',idiag_divum)
        call parse_name(iname,cname(iname),cform(iname),'rdivum',idiag_rdivum)
        call parse_name(iname,cname(iname),cform(iname),'divu2m',idiag_divu2m)
        call parse_name(iname,cname(iname),cform(iname),'gdivu2m',idiag_gdivu2m)
        call parse_name(iname,cname(iname),cform(iname),'drurms',idiag_divrhourms)
        call parse_name(iname,cname(iname),cform(iname),'drumax',idiag_divrhoumax)
        call parse_name(iname,cname(iname),cform(iname),'u3u21m',idiag_u3u21m)
        call parse_name(iname,cname(iname),cform(iname),'u1u32m',idiag_u1u32m)
        call parse_name(iname,cname(iname),cform(iname),'u2u13m',idiag_u2u13m)
        call parse_name(iname,cname(iname),cform(iname),'u2u31m',idiag_u2u31m)
        call parse_name(iname,cname(iname),cform(iname),'u3u12m',idiag_u3u12m)
        call parse_name(iname,cname(iname),cform(iname),'u1u23m',idiag_u1u23m)
        call parse_name(iname,cname(iname),cform(iname),'uxpt',idiag_uxpt)
        call parse_name(iname,cname(iname),cform(iname),'uypt',idiag_uypt)
        call parse_name(iname,cname(iname),cform(iname),'uzpt',idiag_uzpt)
        call parse_name(iname,cname(iname),cform(iname),'uxp2',idiag_uxp2)
        call parse_name(iname,cname(iname),cform(iname),'uyp2',idiag_uyp2)
        call parse_name(iname,cname(iname),cform(iname),'uzp2',idiag_uzp2)
        call parse_name(iname,cname(iname),cform(iname),'uxuypt',idiag_uxuypt)
        call parse_name(iname,cname(iname),cform(iname),'uyuzpt',idiag_uyuzpt)
        call parse_name(iname,cname(iname),cform(iname),'uzuxpt',idiag_uzuxpt)
        call parse_name(iname,cname(iname),cform(iname),'fintm',idiag_fintm)
        call parse_name(iname,cname(iname),cform(iname),'fextm',idiag_fextm)
        call parse_name(iname,cname(iname),cform(iname),'divuHrms',idiag_divuHrms)
        call parse_name(iname,cname(iname),cform(iname),'uxxrms',idiag_uxxrms)
        call parse_name(iname,cname(iname),cform(iname),'uyyrms',idiag_uyyrms)
        call parse_name(iname,cname(iname),cform(iname),'uxzrms',idiag_uxzrms)
        call parse_name(iname,cname(iname),cform(iname),'uyzrms',idiag_uyzrms)
        call parse_name(iname,cname(iname),cform(iname),'uzyrms',idiag_uzyrms)
        call parse_name(iname,cname(iname),cform(iname),'duxdzma',idiag_duxdzma)
        call parse_name(iname,cname(iname),cform(iname),'duydzma',idiag_duydzma)
        call parse_name(iname,cname(iname),cform(iname),'totangmom',idiag_totangmom)
        call parse_name(iname,cname(iname),cform(iname),'uxfampm',idiag_uxfampm)
        call parse_name(iname,cname(iname),cform(iname),'uyfampm',idiag_uyfampm)
        call parse_name(iname,cname(iname),cform(iname),'uzfampm',idiag_uzfampm)
        call parse_name(iname,cname(iname),cform(iname),'uxfampim',idiag_uxfampim)
        call parse_name(iname,cname(iname),cform(iname),'uyfampim',idiag_uyfampim)
        call parse_name(iname,cname(iname),cform(iname),'uzfampim',idiag_uzfampim)
        call parse_name(iname,cname(iname),cform(iname),'uguxm',idiag_uguxm)
        call parse_name(iname,cname(iname),cform(iname),'uguym',idiag_uguym)
        call parse_name(iname,cname(iname),cform(iname),'uguzm',idiag_uguzm)
        call parse_name(iname,cname(iname),cform(iname),'dudx',idiag_dudx)
        call parse_name(iname,cname(iname),cform(iname),'ugu2m',idiag_ugu2m)
        call parse_name(iname,cname(iname),cform(iname),'ugurmsx',idiag_ugurmsx)
        call parse_name(iname,cname(iname),cform(iname),'uxglnrym',idiag_uxglnrym)
        call parse_name(iname,cname(iname),cform(iname),'uyglnrxm',idiag_uyglnrxm)
        call parse_name(iname,cname(iname),cform(iname),'uzdivum',idiag_uzdivum)
        call parse_name(iname,cname(iname),cform(iname),'uxuydivum',idiag_uxuydivum)
        call parse_name(iname,cname(iname),cform(iname),'udpxxm',idiag_udpxxm)
        call parse_name(iname,cname(iname),cform(iname),'udpyym',idiag_udpyym)
        call parse_name(iname,cname(iname),cform(iname),'udpzzm',idiag_udpzzm)
        call parse_name(iname,cname(iname),cform(iname),'udpxym',idiag_udpxym)
        call parse_name(iname,cname(iname),cform(iname),'udpyzm',idiag_udpyzm)
        call parse_name(iname,cname(iname),cform(iname),'udpxzm',idiag_udpxzm)
        call parse_name(iname,cname(iname),cform(iname),'taufmin',idiag_taufmin)
        call parse_name(iname,cname(iname),cform(iname),'dtF',idiag_dtF)
        call parse_name(iname,cname(iname),cform(iname),'nshift',idiag_nshift)
        call parse_name(iname,cname(iname),cform(iname),'uduum',idiag_uduum)
      enddo
!
      if (idiag_u2tm/=0) then
        if (iuut==0) call stop_it("Cannot calculate u2tm if iuut==0")
      endif
      if (idiag_outm/=0) then
        if (iuut==0) call stop_it("Cannot calculate outm if iuut==0")
      endif
      if (idiag_uotm/=0) then
        if (ioot==0) call stop_it("Cannot calculate uotm if ioot==0")
      endif
!
!  Loop over spherical harmonic modes.
!
      if (lroot.and.ip<14) print*,'rprint_hydro: run through parse list'

      do k=1,Nmodes_SH
        smode=itoa(k)
!
!  iname runs through all possible names that may be listed in print.in
!
        do iname=1,nname
          call parse_name(iname,cname(iname),cform(iname), &
              'urlm'//trim(smode),idiag_urlm(k))
        enddo
      enddo
!
! Quantities which are averaged over half (north-south) the box
!
      iname_half=name_half_max
      if ((idiag_urmsn/=0).or.(idiag_urmss/=0))then
        iname_half=iname_half+1
        idiag_urmsh=iname_half
      endif
      if ((idiag_ormsn/=0).or.(idiag_ormss/=0))then
        iname_half=iname_half+1
        idiag_ormsh=iname_half
      endif
      if ((idiag_oumn/=0).or.(idiag_oums/=0))then
        iname_half=iname_half+1
        idiag_oumh=iname_half
      endif
      name_half_max=iname_half
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uxmx',idiag_uxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uymx',idiag_uymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uzmx',idiag_uzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ruxmx',idiag_ruxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ruymx',idiag_ruymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ruzmx',idiag_ruzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rux2mx',idiag_rux2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ruy2mx',idiag_ruy2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ruz2mx',idiag_ruz2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ruxuymx',idiag_ruxuymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ruxuzmx',idiag_ruxuzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ruyuzmx',idiag_ruyuzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'ux2mx',idiag_ux2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uy2mx',idiag_uy2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uz2mx',idiag_uz2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'ox2mx',idiag_ox2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'oy2mx',idiag_oy2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'oz2mx',idiag_oz2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uxuymx',idiag_uxuymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uxuzmx',idiag_uxuzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uyuzmx',idiag_uyuzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'oumx',idiag_oumx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uguxmx',idiag_uguxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uguymx',idiag_uguymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'uguzmx',idiag_uguzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ekinmx',idiag_ekinmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'fkinxmx',idiag_fkinxmx)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'uxmy',idiag_uxmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'uymy',idiag_uymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'uzmy',idiag_uzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'ux2my',idiag_ux2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uy2my',idiag_uy2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uz2my',idiag_uz2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uxuymy',idiag_uxuymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uxuzmy',idiag_uxuzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uyuzmy',idiag_uyuzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'oumy',idiag_oumy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uguxmy',idiag_uguxmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uguymy',idiag_uguymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
            'uguzmy',idiag_uguzmy)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxmz',idiag_uxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uymz',idiag_uymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzmz',idiag_uzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzupmz',idiag_uzupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzdownmz',idiag_uzdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ffdownmz',idiag_ffdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ruzupmz',idiag_ruzupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ruzdownmz',idiag_ruzdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'divumz',idiag_divumz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzdivumz',idiag_uzdivumz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oxmz',idiag_oxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oymz',idiag_oymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ozmz',idiag_ozmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ux2mz',idiag_ux2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uy2mz',idiag_uy2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uz2mz',idiag_uz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uz2upmz',idiag_uz2upmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uz2downmz',idiag_uz2downmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ox2mz',idiag_ox2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oy2mz',idiag_oy2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oz2mz',idiag_oz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ruxmz',idiag_ruxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ruymz',idiag_ruymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ruzmz',idiag_ruzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rux2mz',idiag_rux2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ruy2mz',idiag_ruy2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ruz2mz',idiag_ruz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uxuymz',idiag_uxuymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uxuzmz',idiag_uxuzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uyuzmz',idiag_uyuzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'Rxymz',idiag_Rxymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'Rxyupmz',idiag_Rxyupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'Rxydownmz',idiag_Rxydownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'Rxzmz',idiag_Rxzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'Rxzupmz',idiag_Rxzupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'Rxzdownmz',idiag_Rxzdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'Ryzmz',idiag_Ryzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'Ryzupmz',idiag_Ryzupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'Ryzdownmz',idiag_Ryzdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'ruxuymz',idiag_ruxuymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'ruxuzmz',idiag_ruxuzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'ruyuzmz',idiag_ruyuzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'ruxuy2mz',idiag_ruxuy2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'ruxuz2mz',idiag_ruxuz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), 'ruyuz2mz',idiag_ruyuz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oxuxxmz',idiag_oxuxxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oyuxymz',idiag_oyuxymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oxuyxmz',idiag_oxuyxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oyuyymz',idiag_oyuyymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oxuzxmz',idiag_oxuzxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oyuzymz',idiag_oyuzymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uyxuzxmz',idiag_uyxuzxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uyyuzymz',idiag_uyyuzymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uyzuzzmz',idiag_uyzuzzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'fmasszmz',idiag_fmasszmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'fkinzmz',idiag_fkinzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'fkinzdownmz',idiag_fkinzdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'fkinzupmz',idiag_fkinzupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'ekinmz',idiag_ekinmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'u2mz',idiag_u2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'o2mz',idiag_o2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'curlru2mz',idiag_curlru2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'divru2mz',idiag_divru2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'divu2mz',idiag_divu2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oumz',idiag_oumz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uguxmz',idiag_uguxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uguymz',idiag_uguymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'uguzmz',idiag_uguzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Remz',idiag_Remz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oguxmz', idiag_oguxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oguymz', idiag_oguymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oguzmz', idiag_oguzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ogux2mz',idiag_ogux2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oguy2mz',idiag_oguy2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oguz2mz',idiag_oguz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oxdivumz', idiag_oxdivumz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oydivumz', idiag_oydivumz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ozdivumz', idiag_ozdivumz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oxdivu2mz',idiag_oxdivu2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oydivu2mz',idiag_oydivu2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ozdivu2mz',idiag_ozdivu2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'u3u21mz',idiag_u3u21mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'u1u32mz',idiag_u1u32mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'u2u13mz',idiag_u2u13mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'u2u31mz',idiag_u2u31mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'u3u12mz',idiag_u3u12mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'u1u23mz',idiag_u1u23mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'acczmz',idiag_acczmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'acczupmz',idiag_acczupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'acczdownmz',idiag_acczdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'accpowzmz',idiag_accpowzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'accpowzupmz',idiag_accpowzupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'accpowzdownmz',idiag_accpowzdownmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'totalforcezmz',idiag_totalforcezmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'totalforcezupmz',idiag_totalforcezupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'totalforcezdownmz',idiag_totalforcezdownmz)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do ixz=1,nnamexz
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uxmxz',idiag_uxmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uymxz',idiag_uymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uzmxz',idiag_uzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'ux2mxz',idiag_ux2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uy2mxz',idiag_uy2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uz2mxz',idiag_uz2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uxuymxz',idiag_uxuymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uxuzmxz',idiag_uxuzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'uyuzmxz',idiag_uyuzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'oumxz',idiag_oumxz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxmxy',idiag_uxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uymxy',idiag_uymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uzmxy',idiag_uzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ffdownmxy',idiag_ffdownmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxupmxy',idiag_uxupmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxdownmxy',idiag_uxdownmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruxupmxy',idiag_ruxupmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruxdownmxy',idiag_ruxdownmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ux2upmxy',idiag_uxupmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ux2downmxy',idiag_uxdownmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxuymxy',idiag_uxuymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxuzmxy',idiag_uxuzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uyuzmxy',idiag_uyuzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Rxymxy',idiag_Rxymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Rxyupmxy',idiag_Rxyupmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Rxydownmxy',idiag_Rxydownmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Rxzmxy',idiag_Rxzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Rxzupmxy',idiag_Rxzupmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Rxzdownmxy',idiag_Rxzdownmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Ryzmxy',idiag_Ryzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Ryzupmxy',idiag_Ryzupmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'Ryzdownmxy',idiag_Ryzdownmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'oxmxy',idiag_oxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'oymxy',idiag_oymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ozmxy',idiag_ozmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'oumxy',idiag_oumxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'pvzmxy',idiag_pvzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruxmxy',idiag_ruxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruymxy',idiag_ruymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruzmxy',idiag_ruzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ux2mxy',idiag_ux2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uy2mxy',idiag_uy2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uz2mxy',idiag_uz2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'rux2mxy',idiag_rux2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruy2mxy',idiag_ruy2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruz2mxy',idiag_ruz2mxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruxuymxy',idiag_ruxuymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruxuzmxy',idiag_ruxuzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ruyuzmxy',idiag_ruyuzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'fkinxmxy',idiag_fkinxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'fkinxupmxy',idiag_fkinxupmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'fkinxdownmxy',idiag_fkinxdownmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'fkinymxy',idiag_fkinymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uguxmxy',idiag_uguxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uguymxy',idiag_uguymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uguzmxy',idiag_uguzmxy)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'urmphi',idiag_urmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ursphmphi',idiag_ursphmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uthmphi',idiag_uthmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rursphmphi',idiag_rursphmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ruthmphi',idiag_ruthmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'upmphi',idiag_upmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uzmphi',idiag_uzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rurmphi',idiag_rurmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rupmphi',idiag_rupmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ruzmphi',idiag_ruzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ur2mphi',idiag_ur2mphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'up2mphi',idiag_up2mphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uz2mphi',idiag_uz2mphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'u2mphi',idiag_u2mphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'urupmphi',idiag_urupmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uruzmphi',idiag_uruzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'upuzmphi',idiag_upuzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rurupmphi',idiag_rurupmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ruruzmphi',idiag_ruruzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rupuzmphi',idiag_rupuzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'oumphi',idiag_oumphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ozmphi',idiag_ozmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'fkinrsphmphi',idiag_fkinrsphmphi)
      enddo
!
!  check for those quantities for which we want phiz-averages
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'urmr',  idiag_urmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'upmr',  idiag_upmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'uzmr',  idiag_uzmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'ormr',  idiag_ormr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'opmr',  idiag_opmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'ozmr',  idiag_ozmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'u2mr',  idiag_u2mr)
      enddo
!
!  check for those quantities for which we want video slices
!
      idum=0
      do inamev=1,nnamev
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'uu',  idum)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'oo',  ivid_oo)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'o2',  ivid_o2)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'ou',  ivid_ou)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'divu',ivid_divu)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'u2',  ivid_u2)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'Ma2', ivid_Ma2)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'uu_sph',ivid_uu_sph)
      enddo
!
!  write column where which hydro variable is stored
!
      if (lwr) then
        if (lhelmholtz_decomp) call farray_index_append('iphiuu',iphiuu)
      endif
!
    endsubroutine rprint_hydro
!***********************************************************************
    subroutine get_slices_hydro(f,slices)
!
!  Write slices for animation of Hydro variables.
!
!  26-jul-06/tony: coded
!  12-apr-16/MR: modifications for Yin-Yang grid
!
      use General, only: transform_thph_yy_other
      use Slices_methods, only: assign_slices_scal, assign_slices_vec

      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Velocity field.
!
        case ('uu'); call assign_slices_vec(slices,f,iuu)
!
!  Divergence of velocity.
!
        case ('divu')
          call assign_slices_scal(slices,divu_xy,divu_xz,divu_yz,divu_xy2,divu_xy3,divu_xy4,divu_xz2,divu_r)
!
!  Velocity squared.
!
        case ('u2')
          call assign_slices_scal(slices,u2_xy,u2_xz,u2_yz,u2_xy2,u2_xy3,u2_xy4,u2_xz2,u2_r)
!
!  Vorticity.
!
        case ('oo')
          call assign_slices_vec(slices,oo_xy,oo_xz,oo_yz,oo_xy2,oo_xy3,oo_xy4,oo_xz2,oo_r)
!
!  Vorticity squared.
!
        case ('o2')
          call assign_slices_scal(slices,o2_xy,o2_xz,o2_yz,o2_xy2,o2_xy3,o2_xy4,o2_xz2,o2_r)
!
!  kinetic helicity.
!
        case ('ou')
          call assign_slices_scal(slices,ou_xy,ou_xz,ou_yz,ou_xy2,ou_xy3,ou_xy4,ou_xz2,ou_r)
!
!  Mach number squared.
!
        case ('Ma2')
          call assign_slices_scal(slices,mach_xy,mach_xz,mach_yz,mach_xy2,mach_xy3,mach_xy4,mach_xz2,mach_r)
!
!  Velocity in spherical coordinates
!
        case ('uu_sph')
          call assign_slices_vec(slices,uu_sph_xy,uu_sph_xz,uu_sph_yz,uu_sph_xy2, &
                                        uu_sph_xy3,uu_sph_xy4,uu_sph_xz2,uu_sph_r)
!
      endselect
!
    endsubroutine get_slices_hydro
!***********************************************************************
    function decomp_prepare() result (ldecomp)
!
!  Prepare for Helmholtz decomposition.
!
!  20-oct-97/axel: coded
!
      use Sub, only: read_snaptime, update_snaptime
!
      logical :: ldecomp

      character (len=fnlen) :: file
      integer :: ndummy
      real :: tdummy
!
!  Perform the decomposition in dsnap_down time intervals.
!
      file = trim(datadir)//'/tsnap_down.dat'
!
!  This routine sets ldecomp=T whenever its time to perform the decomposition.
!
      call update_snaptime(file,tdummy,ndummy,dsnap_down,t,ldecomp,nowrite=.true.)
!
    endfunction decomp_prepare
!***********************************************************************
    subroutine hydro_after_timestep(f,df,dt_sub)
!
!  Hook for modification of the f and df arrays
!  according to the hydro module, after the
!  timestep is performed.
!
!  12-mar-17/wlyra: coded.
!  28-mar-17/MR: reinstated update_ghosts.
!
      use Boundcond, only: update_ghosts
      use Sub, only: div
      use Poisson, only: inverse_laplacian, inverse_laplacian_fft_z    !, inverse_laplacian_z_2nd_neumann
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dt_sub
!
      logical :: lwrite_debug=.false.
      integer :: iorder_z=2
!
      fargo: if (lfargo_advection) then
!
!  Call update ghosts as derivatives are needed.
!
        call update_ghosts(f)
!
        fourier: if (lfargoadvection_as_shift) then
          call fourier_shift_fargo(f,df,dt_sub)
        else
          if (lroot) then
            print*,'Fargo advection without Fourier shift'
            print*,'is not functional.'
            call fatal_error("hydro_after_timestep","")
          endif
        endif fourier
!
!  To disable radial advection, intended for tests in cylindrical coordinates
!
        if (lno_radial_advection) then
          f(:,:,:,iux) = 0.
          df(:,:,:,iux) = 0.
        endif
!
      endif fargo

      if (lhelmholtz_decomp) then 

        call fatal_error("hydro_after_timestep","Helmholtz decomposition not yet operational")

        if (decomp_prepare()) then
!
!  Find the divergence of uu and put it into the f-slot which is later used for the flow potential.
!
          do n=n1,n2; do m=m1,m2
            call div(f,iuu,f(l1:l2,m,n,iphiuu))
          enddo; enddo
!
!print*, 'minmax(div)=', minval(f(l1:l2,m1:m2,n1:n2,iphiuu)),maxval(f(l1:l2,m1:m2,n1:n2,iphiuu))
!          if (lwrite_debug) write(31) f(l1:l2,m1:m2,n1:n2,iphiuu)
!
          if (lperi(3)) then
            call inverse_laplacian(f(l1:l2,m1:m2,n1:n2,iphiuu))
          else
            if (iorder_z==2) then
              !call inverse_laplacian_z_2nd_neumann(f)
            else
              ! call inverse_laplacian_fft_z(f(l1:l2,m1:m2,n1:n2,iphiuu))
              ! call inverse_laplacian_z(f(l1:l2,m1:m2,n1:n2,iphiuu),(/'n1s','n1s'/),f(:,:,n1,iuz),f(:,:,n2,iuz))
            endif
          endif
          if (lwrite_debug) write(32) f(l1:l2,4,n1:n2,iphiuu)

        endif
      endif
!
    endsubroutine hydro_after_timestep
!***********************************************************************
    subroutine fourier_shift_fargo(f,df,dt_)
!
!  Add the fargo shift to the f and df-array, in
!  fourier space.
!
!  12-mar-17/wlyra: moved here from special
!
      use Sub
      use Fourier, only: fft_y_parallel,fft_z_parallel
      use Cdata
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,ny) :: acyl_re,acyl_im
      real, dimension (nz) :: asph_re,asph_im
      real, dimension (nx) :: phidot
      integer :: ivar,ig,i
      real :: dt_
!
!  Pencil uses linear velocity. Fargo will shift based on
!  angular velocity. Get phidot from uphi.
!
      ifcoordinates: if (lcylindrical_coords) then
        zloopcyl: do n=n1,n2
          phidot=uu_average_cyl(l1:l2,n)*rcyl_mn1
!
          varloopcyl: do ivar=1,mvar
!
            acyl_re=f(l1:l2,m1:m2,n,ivar)
            acyl_im=0.
!
!  Forward transform. No need for computing the imaginary part.
!  The transform is just a shift in y, so no need to compute
!  the x-transform either.
!
            call fft_y_parallel(acyl_re,acyl_im,SHIFT_Y=phidot*dt_,lneed_im=.false.)
!
!  Inverse transform of the shifted array back into real space.
!  No need again for either imaginary part of x-transform.
!
            call fft_y_parallel(acyl_re,acyl_im,linv=.true.)
            f(l1:l2,m1:m2,n,ivar)=acyl_re
!
!  Also shift df, unless it is the last subtimestep.
!
            if (.not.llast) then
              acyl_re=df(l1:l2,m1:m2,n,ivar)
              acyl_im=0.
              call fft_y_parallel(acyl_re,acyl_im,SHIFT_Y=phidot*dt_,lneed_im=.false.)
              call fft_y_parallel(acyl_re,acyl_im,linv=.true.)
              df(l1:l2,m1:m2,n,ivar)=acyl_re
            endif
!
          enddo varloopcyl
        enddo zloopcyl
      elseif (lspherical_coords) then
        yloopsph: do m=m1,m2
          xloopsph: do i=l1,l2
            ig=i-l1+1
            phidot=uu_average_sph(l1:l2,m)*rcyl_mn1
!
            varloopsph: do ivar=1,mvar
!
              asph_re=f(i,m,n1:n2,ivar)
              asph_im=0.
!
!  Forward transform. No need for computing the imaginary part.
!  The transform is just a shift in z, so no need to compute
!  the x-transform either.
!
              call fft_z_parallel(asph_re,asph_im,SHIFT_Z=phidot(ig)*dt_,lneed_im=.false.)
!
!  Inverse transform of the shifted array back into real space.
!  No need again for either imaginary part of x-transform.
!
              call fft_z_parallel(asph_re,asph_im,linv=.true.)
              f(i,m,n1:n2,ivar)=asph_re
!
!  Also shift df, unless it is the last subtimestep.
!
              if (.not.llast) then
                asph_re=df(i,m,n1:n2,ivar)
                asph_im=0.
                call fft_z_parallel(asph_re,asph_im,SHIFT_Z=phidot(ig)*dt_,lneed_im=.false.)
                call fft_z_parallel(asph_re,asph_im,linv=.true.)
                df(i,m,n1:n2,ivar)=asph_re
              endif
!
            enddo varloopsph
          enddo xloopsph
        enddo yloopsph
      endif ifcoordinates
!
    endsubroutine fourier_shift_fargo
!***********************************************************************
    subroutine calc_mflow
!
!  calculate mean flow field from xy- or z-averages
!
!   8-nov-02/axel: adapted from calc_mfield
!   9-nov-02/axel: allowed mean flow to be compressible
!  24-aug-15/MR: corrected declaration of umx2
!
      use Diagnostics, only: save_name
      use Mpicomm, only: mpibcast_real, mpireduce_sum
!
      logical,save :: first=.true.
      real, dimension (nx,ny) :: fsumxy
      real, dimension (nx) :: uxmx,uymx,uzmx,umx2
      real, dimension (ny) :: uxmy,uymy,uzmy,umy2
      real :: umx,umy,umz
!
!  For vector output (of oo vectors) we need orms
!  on all processors. It suffices to have this for times when lout=.true.,
!  but we need to broadcast the result to all procs.
!
!  calculate orms (this requires that orms is set in print.in)
!  broadcast result to other processors
!
      if (idiag_orms/=0) then
        if (iproc==0) orms=fname(idiag_orms)
        call mpibcast_real(orms)
      endif
!
!  Magnetic energy in vertically averaged field. The uymxy and uzmxy must
!  have been calculated, so they are present on the z-root processors.
!
      if (idiag_umx/=0) then
        if (idiag_uymxy==0.or.idiag_uzmxy==0) then
          if (first) print*, 'calc_mflow:                    WARNING'
          if (first) print*, &
                  "calc_mflow: NOTE: to get umx, uymxy and uzmxy must also be set in zaver"
          if (first) print*, &
                  "calc_mflow:      We proceed, but you'll get umx=0"
          umx=0.
        else
          if (lfirst_proc_z) then
            call mpireduce_sum(fnamexy(idiag_uxmxy,:,:),fsumxy,(/nx,ny/),idir=2)
            uxmx=sum(fsumxy,dim=2)/nygrid
            call mpireduce_sum(fnamexy(idiag_uymxy,:,:),fsumxy,(/nx,ny/),idir=2)
            uymx=sum(fsumxy,dim=2)/nygrid
            call mpireduce_sum(fnamexy(idiag_uzmxy,:,:),fsumxy,(/nx,ny/),idir=2)
            uzmx=sum(fsumxy,dim=2)/nygrid
          endif
          if (lfirst_proc_yz) &
            call mpireduce_sum(uxmx**2+uymx**2+uzmx**2,umx2,nx,idir=1)
          umx=sqrt(sum(umx2)/nxgrid)
        endif
        call save_name(umx,idiag_umx)
      endif
!
!  Similarly for umy.
!
      if (idiag_umy/=0) then
        if (idiag_uxmxy==0.or.idiag_uzmxy==0) then
          if (first) print*, 'calc_mflow:                    WARNING'
          if (first) print*, &
                  "calc_mflow: NOTE: to get umy, uxmxy and uzmxy must also be set in zaver"
          if (first) print*, &
                  "calc_mflow:       We proceed, but you'll get umy=0"
          umy=0.
        else
          if (lfirst_proc_z) then
            call mpireduce_sum(fnamexy(idiag_uxmxy,:,:),fsumxy,(/nx,ny/),idir=1)
            uxmy=sum(fsumxy,dim=1)/nxgrid
            call mpireduce_sum(fnamexy(idiag_uymxy,:,:),fsumxy,(/nx,ny/),idir=1)
            uymy=sum(fsumxy,dim=1)/nxgrid
            call mpireduce_sum(fnamexy(idiag_uzmxy,:,:),fsumxy,(/nx,ny/),idir=1)
            uzmy=sum(fsumxy,dim=1)/nxgrid
          endif
          if (lfirst_proc_xz) &
            call mpireduce_sum(uxmy**2+uymy**2+uzmy**2,umy2,ny,idir=2)
          umy=sqrt(sum(umy2)/nygrid)
        endif
        call save_name(umy,idiag_umy)
      endif
!
!  Kinetic energy in horizontally averaged flow. The uxmz and uymz must
!  have been calculated, so they are present on the root processor.
!
      if (idiag_umz/=0) then
        if (idiag_uxmz==0.or.idiag_uymz==0.or.idiag_uzmz==0) then
          if (first) print*,"calc_mflow:                    WARNING"
          if (first) print*, &
                  "calc_mflow: NOTE: to get umz, uxmz, uymz and uzmz must also be set in xyaver"
          if (first) print*, &
                  "calc_mflow:       This may be because we renamed zaver.in into xyaver.in"
          if (first) print*, &
                  "calc_mflow:       We proceed, but you'll get umz=0"
          umz=0.
        else
          umz=sqrt(sum(fnamez(:,:,idiag_uxmz)**2 &
                      +fnamez(:,:,idiag_uymz)**2 &
                      +fnamez(:,:,idiag_uzmz)**2)/(nz*nprocz))
        endif
        call save_name(umz,idiag_umz)
      endif
!
!  calculation of <U>.<B> in separate subroutine.
!  Should do the same for <Ux^2>, <Uy^2>, <Uz^2> later (as in magnetic)
!
      if (idiag_omumz/=0) call calc_omumz
      if (idiag_umamz/=0) call calc_umamz
      if (idiag_umbmz/=0) call calc_umbmz
      if (idiag_umxbmz/=0) call calc_umxbmz
!
      first = .false.
!
    endsubroutine calc_mflow
!***********************************************************************
    subroutine calc_omumz
!
!  Calculate kinetic helicity of mean field
!  The oxmz and oymz as well as uxmz and uymz must have been calculated,
!  so they are present on the root processor.
!
!  14-feb-09/axel: adapted from calc_umbmz
!
      use Diagnostics, only: save_name
!
      logical,save :: first=.true.
      real :: omumz
!
!  This only works if uxmz, uymz, bxmz, bymz, are in xyaver,
!  so print warning if this is not ok.
!
      if (idiag_oxmz==0.or.idiag_oymz==0.or.idiag_uxmz==0.or.idiag_uymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get omumz, set uxmz, uymz, oxmz, and oymz in xyaver"
          print*,"We proceed, but you'll get omumz=0"
        endif
        omumz=0.
      else
        omumz=sum(fnamez(:,:,idiag_oxmz)*fnamez(:,:,idiag_uxmz) &
                 +fnamez(:,:,idiag_oymz)*fnamez(:,:,idiag_uymz))/(nz*nprocz)
      endif
!
!  save the name in the idiag_omumz slot
!  and set first to false
!
      call save_name(omumz,idiag_omumz)
      first=.false.
!
    endsubroutine calc_omumz
!***********************************************************************
    subroutine calc_umamz
!
!  Cross helicity production of mean field
!  The uxmz and uymz as well as axmz and aymz must have been calculated,
!  so they are present on the root processor.
!
!   5-mar-10/axel: adapted from calc_umbmz
!
      use Diagnostics, only: save_name
      use Magnetic, only: idiag_axmz, idiag_aymz
!
      logical,save :: first=.true.
      real :: umamz
!
!  This only works if uxmz, uymz, axmz, aymz, are in xyaver,
!  so print warning if this is not ok.
!
      if (idiag_uxmz==0.or.idiag_uymz==0.or.idiag_axmz==0.or.idiag_aymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get umamz, set uxmz, uymz, axmz, and aymz in xyaver"
          print*,"We proceed, but you'll get umamz=0"
        endif
        umamz=0.
      else
        umamz=sum(fnamez(:,:,idiag_uxmz)*fnamez(:,:,idiag_axmz) &
                 +fnamez(:,:,idiag_uymz)*fnamez(:,:,idiag_aymz))/(nz*nprocz)
      endif
!
!  save the name in the idiag_umamz slot
!  and set first to false
!
      call save_name(umamz,idiag_umamz)
      first=.false.
!
    endsubroutine calc_umamz
!***********************************************************************
    subroutine calc_umbmz
!
!  Cross helicity production of mean field
!  The uxmz and uymz as well as bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
!  26-jan-09/axel: adapted from calc_ebmz
!
      use Diagnostics, only: save_name
      use Magnetic, only: idiag_bxmz,idiag_bymz
!
      logical,save :: first=.true.
      real :: umbmz
!
!  This only works if uxmz, uymz, bxmz, bymz, are in xyaver,
!  so print warning if this is not ok.
!
      if (idiag_uxmz==0.or.idiag_uymz==0.or.idiag_bxmz==0.or.idiag_bymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get umbmz, set uxmz, uymz, bxmz, and bymz in xyaver"
          print*,"We proceed, but you'll get umbmz=0"
        endif
        umbmz=0.
      else
        umbmz=sum(fnamez(:,:,idiag_uxmz)*fnamez(:,:,idiag_bxmz) &
                 +fnamez(:,:,idiag_uymz)*fnamez(:,:,idiag_bymz))/(nz*nprocz)
      endif
!
!  save the name in the idiag_umbmz slot
!  and set first to false
!
      call save_name(umbmz,idiag_umbmz)
      first=.false.
!
    endsubroutine calc_umbmz
!***********************************************************************
    subroutine calc_umxbmz
!
!  EMF of xy-averaged mean velocity and magnetic fields
!  The uxmz and uymz as well as bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
!  17-mar-09/axel: adapted from calc_umbmz
!
      use Diagnostics, only: save_name
      use Magnetic, only: idiag_bxmz,idiag_bymz
!
      logical,save :: first=.true.
      real :: umxbmz
!
!  This only works if uxmz, uymz, bxmz, bymz, are in xyaver,
!  so print warning if this is not ok.
!
      if (idiag_uxmz==0.or.idiag_uymz==0.or.idiag_bxmz==0.or.idiag_bymz==0) then
        if (first) then
          print*,"calc_mfield: WARNING"
          print*,"NOTE: to get umxbmz, set uxmz, uymz, bxmz, and bymz in xyaver"
          print*,"We proceed, but you'll get umxbmz=0"
        endif
        umxbmz=0.
      else
        umxbmz=sum(fnamez(:,:,idiag_uxmz)*fnamez(:,:,idiag_bymz) &
                  -fnamez(:,:,idiag_uymz)*fnamez(:,:,idiag_bxmz))/(nz*nprocz)
      endif
!
!  save the name in the idiag_umxbmz slot
!  and set first to false
!
      call save_name(umxbmz,idiag_umxbmz)
      first=.false.
!
    endsubroutine calc_umxbmz
!***********************************************************************
    subroutine remove_mean_momenta(f,indux,indrho)
!
!  Substract mean x-momentum over density from the x-velocity field.
!  Useful to avoid unphysical winds in shearing box simulations.
!  Note: this is possibly not useful when there is rotation, because
!  then epicyclic motions don't usually grow catastrophically.
!
!  15-nov-06/tobi: coded
!  15-dec-10/MR  : added parameters indux, indrho to make routine applicable
!                  to other velocities/densities
!  13-feb-15/MR  : changes for use of reference_state (used for main run density only)
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray), intent(inout)        :: f
      integer,                            intent(in)           :: indux
      integer,                            intent(in), optional :: indrho
!
      real, dimension (nx) :: rho,rho1,mm
      real :: fac
      real, dimension (indux:indux+2) :: rum, rum_tmp
      integer :: m,n,j,indrhol
      logical :: lref
!
!  check if ldensity=T. Otherwise switch to remove_mean_flow.
!
      if (ldensity) then
!
!  initialize mean momentum, rum, to zero
!
        lref=.false.
        if (present(indrho)) then
          indrhol = indrho
        else if (indux==iux) then
!
          if (ldensity_nolog) then
            indrhol = irho
          else
            indrhol = ilnrho
          endif
          lref=lreference_state
!
        else
!
! for testflow
!
          indrhol = indux+3
        endif
!
        rum = 0.0
        fac = 1.0/nwgrid
!
!  Go through all pencils.
!
        do n = n1,n2
        do m = m1,m2
!
!  Compute density from the f-array.
!
          if (ldensity_nolog) then
            rho = f(l1:l2,m,n,indrhol)
            if (lref) rho=rho+reference_state(:,iref_rho)
          else
            rho = exp(f(l1:l2,m,n,indrhol))
          endif
!
!  Compute mean momentum in each of the 3 directions.
!
          do j=indux,indux+2
            mm = rho*f(l1:l2,m,n,j)
            rum(j) = rum(j) + fac*sum(mm)
          enddo
        enddo
        enddo
!
!  Compute total sum for all processors.
!  Allow here for the possibility to add mean_momentum.
!  It needs to be subtracted here, because "rum" is removed.
!
        call mpiallreduce_sum(rum,rum_tmp,3)
        rum = rum_tmp - mean_momentum
!
!  Compute inverse density, rho1.
!
        do n = n1,n2
        do m = m1,m2
          if (ldensity_nolog) then
            if (lref) then
              rho1 = 1./(f(l1:l2,m,n,indrhol)+reference_state(:,iref_rho))
            else
              rho1 = 1./f(l1:l2,m,n,indrhol)
            endif
          else
            rho1 = exp(-f(l1:l2,m,n,indrhol))
          endif
!
!  Subtract out the mean momentum separately for each direction.
!
          do j=indux,indux+2
            f(l1:l2,m,n,j) = f(l1:l2,m,n,j) - rho1*rum(j)
          enddo
        enddo
        enddo
        if (lroot.and.ip<6) print*,'remove_mean_momenta: rum=',rum
      else
        call remove_mean_flow(f,indux)       ! as this is equivalent to remove
                                             ! mean momenta for constant density
      endif
    endsubroutine remove_mean_momenta
!***********************************************************************
    subroutine remove_mean_flow(f,indux)
!
!  Substract mean x-flow from the x-velocity field.
!  Useful to avoid unphysical winds in shearing box simulations.
!  Note: this is possibly not useful when there is rotation, because
!  then epicyclic motions don't usually grow catastrophically.
!
!  22-may-07/axel: adapted from remove_mean_momenta
!  15-dec-10/MR  : added parameters indux to make routine applicable
!                  to other velocities
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer,                            intent (in)    :: indux
!
      real, dimension (indux:indux+2) :: um, um_tmp
      integer :: m,n,j
      real    :: fac
!
!  initialize um and compute normalization factor fac
!
        um = 0.0
        fac = 1.0/nwgrid
!
!  Go through all pencils.
!
        do n = n1,n2
        do m = m1,m2
!
!  Compute mean flow in each of the 3 directions.
!
          do j=indux,indux+2
            um(j) = um(j) + fac*sum(f(l1:l2,m,n,j))
          enddo
        enddo
        enddo
!
!  Compute total sum for all processors
!
        call mpiallreduce_sum(um,um_tmp,3)
        um = um_tmp
!
!  Go through all pencils and subtract out the mean flow
!  separately for each direction.
!
        do n = n1,n2
        do m = m1,m2
          do j=indux,indux+2
            f(l1:l2,m,n,j) = f(l1:l2,m,n,j) - um(j)
          enddo
        enddo
        enddo
!
        if (lroot.and.ip<6) print*,'remove_mean_flow: um=',um
!
    endsubroutine remove_mean_flow
!***********************************************************************
    subroutine remove_mean_angmom(f,induz)
!
!  Substract <L_z>/<rho*sin(theta)> from z-flow. Useful to avoid
!  unphysical accumulation of angular momentum in spherical
!  coordinates.
!
!  29-aug-13/pete: adapted from remove_mean_flow
!  30-jan-15/pete: take reference state into account
!  13-feb-15/MR  : some optimizations
!
      use Mpicomm, only: mpiallreduce_sum
      use DensityMethods, only: getrho
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer,                            intent (in)    :: induz
!
      real, dimension (nx) :: tmp, rho, wx
      real :: um, angmom, angmom_tmp, rhosint, rhosint_tmp, fac, wmn
      integer :: m,n
!
!  Initialize um and compute normalization factor fac
!
      angmom = 0.0
      rhosint = 0.0
      fac = 1./(Lxyz(1)*(cos(y0)-cos(y0+Lxyz(2)))*Lxyz(3))
      wx=x(l1:l2)*dVol_x(l1:l2)
!
!  Go through all pencils.
!
      do n = n1,n2
        do m = m1,m2
!
!  Compute volume integrals of angular momentum and rho*sin(theta)
!
          call getrho(f(:,m,n,ilnrho),rho)
!
          tmp=rho*wx
          wmn=sinth(m)*dVol_y(m)*dVol_z(n)
          angmom=angmom+sum(tmp*f(l1:l2,m,n,induz))*wmn
          rhosint=rhosint+sum(tmp)*wmn
!
        enddo
      enddo
!
      angmom=fac*angmom; rhosint=fac*rhosint
!
!  Compute total sum for all processors
!
      call mpiallreduce_sum(angmom,angmom_tmp)
      call mpiallreduce_sum(rhosint,rhosint_tmp)
      um=angmom_tmp/rhosint_tmp
!
!  Go through all pencils and subtract out the excess u_phi
!
      f(l1:l2,m1:m2,n1:n2,induz) = f(l1:l2,m1:m2,n1:n2,induz) - um
!
      if (lroot.and.ip<6) print*,'remove_mean_angmom: um=',um
!
    endsubroutine remove_mean_angmom
!***********************************************************************
    subroutine interior_bc_hydro(f)
!
!  Set interior boundary condition within the domain
!
!  11-jun-08/axel: coded
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer :: l1bc,l2bc
!
      select case (interior_bc_hydro_profile)
!
!  single propeller blade
!
      case ('blade')
        z1_interior_bc_hydro=2.
        l1bc=(l1+l2)/2
        l2bc=l1bc+1
        do n=n1,n2
          if (z(n)<z1_interior_bc_hydro) then
            do m=m1,m2
              f(l1bc:l2bc,m,n,iux:iuz)=0.
            enddo
          endif
        enddo
!
!  no profile
!
      case ('nothing')
!
!  no profile matches
!
      case default
        if (lroot) print*,'interior_bc_hydro: No such profile ',interior_bc_hydro_profile
      endselect
!
    endsubroutine interior_bc_hydro
!***********************************************************************
    subroutine impose_profile_diffrot(f,df,prof_diffrot,ldiffrot_test)
!
!  forcing of differential rotation with a -(1/tau)*(u-uref) method
!
!  27-june-2007 dhruba: coded
!
      use Sub, only: step
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=labellen) :: prof_diffrot
      logical :: ldiffrot_test
!
      real, dimension(nx) :: local_Omega

      select case (prof_diffrot)
!
!  diffrot profile from Brandenburg & Sandin (2004, A&A)
!
      case ('BS04')
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(f(l1:l2,m,n,iuy)-prof_amp1*prof_amp3(n))
!
!  diffrot profile from Brandenburg & Sandin (2004, A&A), modified
!  for convection.
!
      case ('BS04c')
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(f(l1:l2,m,n,iuy)-prof_amp1*prof_amp3(n))
!
!  same as above but with an equator
!
      case ('BS04c1')
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(f(l1:l2,m,n,iuy)-prof_amp1*prof_amp3(n))
!
!  modified diffrot profile from Brandenburg & Sandin (2004, A&A)
!
      case ('BS04m')
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(f(l1:l2,m,n,iuz)-prof_amp1*prof_amp4(m))
!
!  Shear profile from Hughes & Proctor (2009, PRL). Force either the
!  full velocity field or only the y-average.
!
      case ('HP09')
      if (.not.lcalc_uumeanxz) then
        !df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_diffrot1* f(l1:l2,m,n,iux)
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(f(l1:l2,m,n,iuy)-prof_amp1*prof_amp3(n))
      else
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(uumxz(l1:l2,n,2)-prof_amp1*prof_amp3(n))
      endif
!
!  Shear profile linear in x
!
      case ('Sx')
      if (lcalc_uumeanx) then
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(uumx(l1:l2,iuy)-Shearx*x(l1:l2))
      else
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(f(l1:l2,m,n,iuy)-Shearx*x(l1:l2))
      endif
!
!  Solar rotation profile from Dikpati & Charbonneau (1999, ApJ)
!
      case ('solar_DC99')
      if (lcalc_uumeanxy) then 
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(uumxy(l1:l2,m,3)-prof_amp1*prof_amp4(m))
      else
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(f(l1:l2,m,n,iuz)-prof_amp1*prof_amp4(m))
      endif
!
!  vertical shear profile
!
      case ('vertical_shear')
      if (.not.lcalc_uumean) then
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(f(l1:l2,m,n,iuy)-prof_amp3(n))
      else
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(uumz(n,2)-prof_amp3(n))
      endif
!
!  vertical compression profile
!
      case ('vertical_compression')
      if (.not.lcalc_uumean) then
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(f(l1:l2,m,n,iuz)-prof_amp3(n))
      else
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(uumz(n,3)-prof_amp3(n))
      endif
!
!  Remove vertical shear profile
!
      case ('remove_vertical_shear')
        f(l1:l2,m,n,iux)=f(l1:l2,m,n,iux)-uumz(n,1)
        f(l1:l2,m,n,iuy)=f(l1:l2,m,n,iuy)-uumz(n,2)
!
!  vertical shear profile
!
      case ('vertical_shear_x')
        df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_diffrot1*(f(l1:l2,m,n,iux)-prof_amp3(n))
!
!  vertical shear profile with sinz profile
!
      case ('vertical_shear_x_sinz')
        df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_diffrot1*(f(l1:l2,m,n,iux)-prof_amp3(n))
!
!  Vertical shear profile U_y(z) centred around given z, forcing the
!  horizontally averaged flow.
!
      case ('vertical_shear_z')
      if (.not.lcalc_uumean) then
         df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(f(l1:l2,m,n,iuy)-prof_amp3(n))
      else
         df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(uumz(n,2)-prof_amp3(n))
      endif
!
!  Vertical shear profile U_y(z) centred around given z, forcing the
!  y-averaged averaged flow.
!
      case ('vertical_shear_z2')
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(uumxz(l1:l2,n,2)-prof_amp3(n))
!
!  Linear vertical shear profile U_y(z), forcing the y-averaged flow.
!
      case ('vertical_shear_linear')
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(uumxz(l1:l2,n,2)-prof_amp3(n))
!
!  set u_phi=0 below given radius, i.e. enforce a tachocline in
!  spherical convection setup
!
      case ('tachocline')
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*prof_amp1*f(l1:l2,m,n,iuz)
!
!  write differential rotation in terms of Gegenbauer polynomials
!  Omega = Omega0 + Omega2*P31(costh)/sinth + Omega4*P51(costh)/sinth + ...
!  Note that P31(theta)/sin(theta) = (3/2) * [1 - 5*cos(theta)^2 ]
!
      case ('solar_simple')
      if (lspherical_coords.or.lcartesian_coords) then
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(f(l1:l2,m,n,iuz) &
                                                          -prof_amp1*prof_amp4(m))
!            -prof_amp1*cos(20.*x(llx))*cos(20.*y(m)) )
      endif
      if (ldiffrot_test) then
        f(l1:l2,m,n,iux:iuy) = 0.
        if (lspherical_coords.or.lcartesian_coords) f(l1:l2,m,n,iuz) = prof_amp1*prof_amp4(m)
      endif
!
!  radial_uniform_shear
!  uphi = slope*x + uoffset
!
      case ('radial_uniform_shear')
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(f(l1:l2,m,n,iuz)-prof_amp1)
!
! vertical shear in uz with respect to uz (slow wind)
!
      case ('breeze')
      if (.not.lcalc_uumean) then
         df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(f(l1:l2,m,n,iuz)-prof_amp3(n))
      else
         df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(uumz(n,3)-prof_amp3(n))
      endif
!
!  uz with respect to z (slow wind) like 'breeze' case but ampl_wind a function of z
!
      case ('slow_wind')
      if (.not.lcalc_uumean) then
         df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(f(l1:l2,m,n,iuz)-prof_amp3(n))
      else
         df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(uumz(n,3)-prof_amp3(n))
      endif
!
!  Radial shear profile
!
      case ('radial_shear')
         df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(uumxy(l1:l2,m,3)-prof_amp1)
!
!  Radial shear profile with damping of other mean velocities
!
      case ('radial_shear_damp')
         df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_diffrot1*(uumxy(l1:l2,m,1)-0.)
         df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(uumxy(l1:l2,m,2)-0.)
         df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(uumxy(l1:l2,m,3)-prof_amp1)
!
!  Damp mean velocities in the corona above rdampext
!
      case ('damp_corona')
      if (lspherical_coords) then
        df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_diffrot1*prof_amp1*(uumxy(l1:l2,m,1)-0.)
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*prof_amp1*(uumxy(l1:l2,m,2)-0.)
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*prof_amp1*(uumxy(l1:l2,m,3)-0.)
      else if (lcartesian_coords) then
        if (.not.lcalc_uumeanz) then
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_diffrot1*prof_amp3(n)*(f(l1:l2,m,n,iux)-0.)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*prof_amp3(n)*(f(l1:l2,m,n,iuy)-0.)
          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*prof_amp3(n)*(f(l1:l2,m,n,iuz)-0.)
        else
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_diffrot1*prof_amp3(n)*(uumz(n,1)-0.)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*prof_amp3(n)*(uumz(n,2)-0.)
          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*prof_amp3(n)*(uumz(n,3)-0.)
        endif
      endif
!
!  Damp horizontal motion
!
      case ('damp_horiz_vel')
        if (lcartesian_coords) then
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_diffrot1*prof_amp3(n)*f(l1:l2,m,n,iux)
          df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*prof_amp3(n)*f(l1:l2,m,n,iuy)
        endif
!
!  Latitudinal shear profile
!
      case ('latitudinal_shear')
         df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(uumxy(l1:l2,m,3)-prof_amp4(m))
!
!  Damp/enforce high latitude jets in spherical coordinates
!
      case ('damp_jets')
         df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*prof_amp4(m)*(uumxy(l1:l2,m,3)-uzjet)
!
!  spoke-like-NSSL (analogous to that in hydro_kinematic)
!
      case ('spoke-like-NSSL')
        local_Omega=profx_diffrot1*profy_diffrot1(m)+profx_diffrot2*profy_diffrot2(m)+profx_diffrot3*profy_diffrot3(m)
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(uumxy(l1:l2,m,3)-prof_amp1*local_Omega*sinth(m))
!
!  Relax toward (1D-z) profile read from file
!
      case ('uumz_profile')
        df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-tau_diffrot1*(uumz(n,1)-uumz_prof(n-nghost,1))
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-tau_diffrot1*(uumz(n,2)-uumz_prof(n-nghost,2))
        if (.not.limpose_only_horizontal_uumz) &
           df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(uumz(n,3)-uumz_prof(n-nghost,3))
!
!  Relax toward (2D-xy) profile read from file
!
      case ('omega_profile')
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-tau_diffrot1*(uumxy(l1:l2,m,3)-omega_prof(:,m-nghost)*x(l1:l2)*sinth(m))
!
      case default

      endselect
!
    endsubroutine impose_profile_diffrot
!***********************************************************************
    subroutine read_uumz_profile(uumz_prof)
!
!  Read vertical profiles of horizontally averaged flows from file.
!
!  23-oct-18/pjk: added
!  20-apr-21/MR: use mpiscatter
!
      use Mpicomm, only: mpibcast_real_arr, mpiscatter, &
                         MPI_COMM_ZBEAM, MPI_COMM_XYPLANE
!
      real, dimension(nz,3), intent(out) :: uumz_prof
!
      real, dimension(nzgrid) :: tmp1z, tmp2z, tmp3z
      real :: var1,var2,var3
      logical :: exist
      integer :: stat, n
!
      uumz_prof = 0.

      if (lroot) then
!
!  Read mean velocity (uumz) and write into an array.
!  If file is not found in run directory, search under trim(directory).
!
        inquire(file='uumz.dat',exist=exist)
        if (exist) then
          open(31,file='uumz.dat')
        else
          inquire(file=trim(directory)//'/uumz.ascii',exist=exist)
          if (exist) then
            open(31,file=trim(directory)//'/uumz.ascii')
          else
            call fatal_error('read_uumz_profile','no input file uumz.dat or '//trim(directory)//' uumz.ascii')
          endif
        endif
!
!  Read profiles. Assuming no ghost zones in uumz.dat.
!
        do n=1,nzgrid
          read(31,*,iostat=stat) var1,var2,var3
          if (stat<0) exit      !MR: then data missing!
          if (ip<5) print*,'uxmz, uymz, uzmz: ',var1,var2,var3
          tmp1z(n)=var1
          tmp2z(n)=var2
          tmp3z(n)=var3
        enddo
!
        close(31)

      endif
!
!  Distribute tmp*z along z-beams, then replicate across xy-planes.
!
      if (ipx==0.and.ipy==0) then
        call mpiscatter(tmp1z, uumz_prof(:,1), root, comm=MPI_COMM_ZBEAM)
        call mpiscatter(tmp2z, uumz_prof(:,2), root, comm=MPI_COMM_ZBEAM)
        call mpiscatter(tmp3z, uumz_prof(:,3), root, comm=MPI_COMM_ZBEAM)
      endif
      call mpibcast_real_arr(uumz_prof, 3*nz, proc=0, comm=MPI_COMM_XYPLANE)
!
!print*, 'iproc,uumz_prof=', iproc,uumz_prof(:,1)
!
    endsubroutine read_uumz_profile
!***********************************************************************
    subroutine read_omega_profile(omega_prof)
!
!  Read profile of angular velocity as functio nof x and y from file
!  omega.dat or omega.ascii .
!
!  20-apr-21/MR: coded
!
      use Mpicomm, only: mpibcast_real_arr, mpiscatter, &
                         MPI_COMM_ZBEAM, MPI_COMM_XYPLANE
!
      real, dimension(nx,ny), intent(out) :: omega_prof

      real, dimension(:,:), allocatable :: omega_prof_glob
      logical :: exist
      integer :: stat
!
      omega_prof = 0.
      if (lroot) then
!
!  Read angular velocity omega and write into an array.
!  If file is not found in run directory, search under trim(directory).
!
        inquire(file='omega.dat',exist=exist)
        if (exist) then
          open(31,file='omega.dat')
        else
          inquire(file=trim(directory)//'/omega.ascii',exist=exist)
          if (exist) then
            open(31,file=trim(directory)//'/omega.ascii')
          else
            call fatal_error('read_omega_profile','no input file omega.dat or' &
                             //trim(directory)//' omega.ascii')
          endif
        endif
!
!  Read profile. Assuming no ghost zones in omega.dat.
!
        allocate(omega_prof_glob(nxgrid,nygrid))
        read(31,*,iostat=stat) omega_prof_glob
        close(31)

      endif
!
!  Distribute omega across xy-planes, then replicate along z-beams.
!
      if (ipz==0) &
        call mpiscatter(omega_prof_glob, omega_prof, root, comm=MPI_COMM_XYPLANE)

      call mpibcast_real_arr(omega_prof, nx*ny, proc=0, comm=MPI_COMM_ZBEAM)
!write(iproc+40,'(8(f8.1,1x))') omega_prof

    endsubroutine read_omega_profile
!***********************************************************************
    subroutine impose_velocity_ceiling(f)
!
!  Impose a maximum velocity by setting all higher velocities to the maximum
!  value (velocity_ceiling). Useful for debugging purposes.
!
!  13-aug-2007/anders: implemented.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      if (velocity_ceiling>0.0) then
        where (f(l1:l2,m,n,iux:iuz)> velocity_ceiling) &
            f(l1:l2,m,n,iux:iuz)= velocity_ceiling
        where (f(l1:l2,m,n,iux:iuz)<-velocity_ceiling) &
            f(l1:l2,m,n,iux:iuz)=-velocity_ceiling
      endif
!
    endsubroutine impose_velocity_ceiling
!***********************************************************************
    subroutine meri_circ(f)
!
! Meridional circulation as initial condition.
!
!  26-apr-2010/dhruba: coded.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: m,n
      real :: rone,theta,theta1
!
      do n=n1,n2
        do m=m1,m2
          rone=xyz0(1)
          theta=y(m)
          theta1=xyz0(2)
!
          f(l1:l2,m,n,iux)=amp_meri_circ*(r1_mn**2)*(sin1th(m))*(&
              2*sin(theta-theta1)*cos(theta-theta1)*cos(theta)&
              -sin(theta)*sin(theta-theta1)**2)*&
              (x(l1:l2)-1.)*(x(l1:l2)-rone)**2
          f(l1:l2,m,n,iuy)=-amp_meri_circ*r1_mn*sin1th(m)*(&
              cos(theta)*sin(theta-theta1)**2)*&
              (x(l1:l2)-rone)*(3*x(l1:l2)-rone-2.)
          f(l1:l2,m,n,iuz)=0.
        enddo
      enddo
!
    endsubroutine  meri_circ
!***********************************************************************
    subroutine hydro_clean_up
!
!  This is a dummy routine.
!
!  8-sep-2009/dhruba: coded
!
      if (ldebug) call warning('hydro_clean_up','Nothing to do for hydro.f90')
!
    endsubroutine hydro_clean_up
!***********************************************************************
    subroutine expand_shands_hydro
!
!  Expands shorthand labels of hydro diagnostics.
!
!  16-may-12/MR: coded
!
      use Diagnostics, only : name_is_present, expand_cname
!
      if (nnamerz>0) then
!
        call expand_cname(cnamerz,nnamerz,name_is_present(cnamerz,'uumphi'),&
                          'urmphi','upmphi','uzmphi')
!
        if (name_is_present(cnamerz,'upmphi')>0) then
          call expand_cname(cnamerz,nnamerz,name_is_present(cnamerz,'uusphmphi'),&
                            'ursphmphi','uthmphi')
        else
          call expand_cname(cnamerz,nnamerz,name_is_present(cnamerz,'uusphmphi'),&
                            'ursphmphi','uthmphi','upmphi')
        endif
      endif
!
    endsubroutine expand_shands_hydro
!***********************************************************************
    subroutine amp_lm(psi,psilm,rselect)

      use Sub, only: ylm

      real,dimension(nx),intent(in):: psi,rselect
      real,dimension(nx,Nmodes_SH),intent(out) :: psilm
      real :: sph_har
      integer :: ell,emm,imode
      real,dimension(nx) :: one_by_rsqr

      one_by_rsqr=1./(x(l1:l2)*x(l1:l2))
      if ((m.eq.1).and.(n.eq.1)) write(*,*) rselect

      do ell=0,lSH_max
        do emm=-ell,ell
          imode=(ell+1)*(ell+1)-ell+emm
          sph_har= ylm(ell,emm)
          psilm(:,imode) = psi*rselect*sph_har*one_by_rsqr
        enddo
      enddo

    endsubroutine amp_lm
!***********************************************************************
    subroutine calc_gradu(f)
!
    use Sub, only : gij
    real, dimension (mx,my,mz,mfarray) :: f
    integer :: imn,jk,jj,kk,nyz
    real, dimension(nx,3,3) :: gradu
!
! Calculates gradu and stores it as an auxiliary. This is expected to be called
! only once either during initialization or post-processing. 
!
    nyz=ny*nz
    do imn=1,nyz
       n=nn(imn)
       m=mm(imn)
       lfirstpoint=(imn==1)      ! true for very first iteration of m-n loop
       llastpoint=(imn==nyz) 
       call gij(f,iuu,gradu,1)
       jk=0
       do jj=1,3; do kk=1,3
         f(l1:l2,m,n,iguij+jk) = gradu(:,jj,kk)
         jk=jk+1
       enddo;enddo
    enddo
!
    endsubroutine calc_gradu
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(lpressuregradient_gas,p_par(1))  ! int

    endsubroutine pushpars2c
!***********************************************************************
    subroutine pushdiags2c(p_diag)

    use Syscalls, only: copy_addr
    use Diagnostics, only: set_type

    integer, parameter :: n_diags=12
    integer(KIND=ikind8), dimension(n_diags) :: p_diag

    call copy_addr(idiag_urms,p_diag(1))
    call set_type(idiag_urms,lsqrt=.true.)
    call copy_addr(idiag_uxrms,p_diag(2))
    call set_type(idiag_uxrms,lsqrt=.true.)
    call copy_addr(idiag_uyrms,p_diag(3))
    call set_type(idiag_uyrms,lsqrt=.true.)
    call copy_addr(idiag_uzrms,p_diag(4))
    call set_type(idiag_uzrms,lsqrt=.true.)
    call copy_addr(idiag_umax,p_diag(5))
    call set_type(idiag_umax,lmax=.true.)
    call copy_addr(idiag_umin,p_diag(6))
    call set_type(idiag_umin,lmin=.true.)
    call copy_addr(idiag_uxmin,p_diag(7))
    call set_type(idiag_uxmin,lmin=.true.)
    call copy_addr(idiag_uymin,p_diag(8))
    call set_type(idiag_uymin,lmin=.true.)
    call copy_addr(idiag_uzmin,p_diag(9))
    call set_type(idiag_uzmin,lmin=.true.)
    call copy_addr(idiag_uxmax,p_diag(10))
    call set_type(idiag_uxmax,lmax=.true.)
    call copy_addr(idiag_uymax,p_diag(11))
    call set_type(idiag_uymax,lmax=.true.)
    call copy_addr(idiag_uzmax,p_diag(12))
    call set_type(idiag_uzmax,lmax=.true.)

    endsubroutine pushdiags2c
!***********************************************************************
endmodule Hydro
