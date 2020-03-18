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
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 3
!
! PENCILS PROVIDED divu; oo(3); o2; ou; u2; uij(3,3); uu(3); curlo(3)
! PENCILS PROVIDED sij(3,3); sij2; uij5(3,3); ugu(3); ugu2; oij(3,3)
! PENCILS PROVIDED d2uidxj(3,3), uijk(3,3,3); ogu(3)
! PENCILS PROVIDED u3u21; u1u32; u2u13; del2u(3); del4u(3); del6u(3)
! PENCILS PROVIDED u2u31; u3u12; u1u23
! PENCILS PROVIDED graddivu(3); del6u_bulk(3); grad5divu(3)
! PENCILS PROVIDED rhougu(3); der6u(3); transpurho(3)
! PENCILS PROVIDED divu0; u0ij(3,3); uu0(3)
! PENCILS PROVIDED uu_advec(3); uuadvec_guu(3)
! PENCILS PROVIDED del6u_strict(3); del4graddivu(3)
!
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
  real, target, dimension (:,:,:), allocatable :: oo_xz, oo_yz, oo_xz2
  real, target, dimension (:,:), allocatable :: divu_xy,u2_xy,o2_xy,mach_xy
  real, target, dimension (:,:), allocatable :: divu_xy2,u2_xy2,o2_xy2,mach_xy2
  real, target, dimension (:,:), allocatable :: divu_xy3,divu_xy4,u2_xy3,u2_xy4,mach_xy4
  real, target, dimension (:,:), allocatable :: o2_xy3,o2_xy4,mach_xy3
  real, target, dimension (:,:), allocatable :: divu_xz,u2_xz,o2_xz,mach_xz
  real, target, dimension (:,:), allocatable :: divu_xz2,u2_xz2,o2_xz2,mach_xz2
  real, target, dimension (:,:), allocatable :: divu_yz,u2_yz,o2_yz,mach_yz

  real, dimension (mz,3) :: uumz
  real, dimension (nz,3) :: guumz=0.0
  real, dimension (mx,3) :: uumx=0.0
  real, dimension (:,:,:), allocatable :: uumxy
  real, dimension (mx,mz,3) :: uumxz=0.0
!
!  phi-averaged arrays for orbital advection
!
  real, dimension (nx,nz) :: uu_average_cyl=0.
  real, dimension (nx,ny) :: uu_average_sph=0.
!
!  Cosine and sine function for setting test fields and analysis.
!
  real, dimension (mz) :: c2z,s2z,cz,sz
  real, dimension (nx) :: Remz
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
  real :: xhalf, kgaussian_uu=0.
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
  real :: u_out_kep=0.0, w_sldchar_hyd=1.0
  real :: mu_omega=0., gap=0., r_omega=0., w_omega=0.
  real :: z1_uu=0., z2_uu=0.
  real :: ABC_A=1., ABC_B=1., ABC_C=1.
!
  integer :: igradu
  logical :: llinearized_hydro=.false.
  logical :: ladvection_velocity=.true.
  logical :: lprecession=.false.
  logical :: lshear_rateofstrain=.false.
  logical :: loo_as_aux = .false.
  logical :: luut_as_aux=.false.,loot_as_aux=.false.
  logical :: lscale_tobox=.true.
  logical, target :: lpressuregradient_gas=.true.
  logical :: lcoriolis_force=.true.
  logical :: lshear_in_coriolis=.false.
  logical :: lcentrifugal_force=.false.
  logical, pointer :: lffree
  logical :: lreflecteddy=.false.,louinit=.false.
  logical :: lskip_projection=.false., lno_second_ampl=.true.
  real, pointer :: profx_ffree(:),profy_ffree(:),profz_ffree(:)
  real :: incl_alpha = 0.0, rot_rr = 0.0
  real :: xsphere = 0.0, ysphere = 0.0, zsphere = 0.0
  real :: amp_meri_circ = 0.0
  real :: max_uu = 0.
! The following is useful to debug the forcing - Dhruba
  real :: outest
  real :: omega_ini=0.0
  logical :: loutest
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
      uy_left, uy_right, uu_const, u_out_kep, &
      initpower, initpower2, cutoff, ncutoff, kpeak, kgaussian_uu, &
      lskip_projection, lno_second_ampl, z1_uu, z2_uu, &
      lcoriolis_force, lcentrifugal_force, ladvection_velocity, &
      lprecession, omega_precession, alpha_precession, &
      loo_as_aux, luut_as_aux, loot_as_aux, mu_omega, gap, &
      lscale_tobox, omega_ini, r_cyl, skin_depth, incl_alpha, &
      rot_rr, xsphere, ysphere, zsphere, amp_meri_circ, &
      rnoise_int, rnoise_ext, lreflecteddy, louinit, hydro_xaver_range, max_uu,&
      amp_factor,kx_uu_perturb,llinearized_hydro, hydro_zaver_range, &
      ll_sh, mm_sh, delta_u, n_xprof
!
!  Run parameters.
!
  real :: ruxm=0.,ruym=0.,ruzm=0.
  real :: ampl1_diffrot=0.,ampl2_diffrot=0., ampl_wind=0.
  real :: nu=0.,xexp_diffrot=1.,kx_diffrot=1.,kz_diffrot=0., phase_diffrot=0.
  real :: othresh=0.,othresh_per_orms=0.,orms=0.,othresh_scl=1.
  real :: utop=0.,ubot=0.,omega_out=0.,omega_in=0.
  real :: width_ff_uu=1.,x1_ff_uu=0.,x2_ff_uu=0.
  real :: ekman_friction=0.0, uzjet=0.0
  real :: ampl_forc=0., k_forc=impossible, w_forc=0., x_forc=0., dx_forc=0.1
  real :: ampl_fcont_uu=1., k_diffrot=1., amp_centforce=1.
  real :: uphi_rbot=1., uphi_rtop=1., uphi_step_width=0., delta_u=1.
  integer :: novec,novecmax=nx*ny*nz/4
  logical :: lupw_uu=.false.
  logical :: lremove_mean_angmom=.false.
  logical :: lremove_mean_momenta=.false.
  logical :: lremove_mean_flow=.false.
  logical :: lremove_uumeanxy=.false.
  logical :: lremove_uumeanz=.false.
  logical :: lremove_uumeanz_horizontal=.false.
  logical :: lreinitialize_uu=.false.
  logical :: lalways_use_gij_etc=.false.
  logical :: lcalc_uumeanz=.false.,lcalc_uumeanxy=.false.,lcalc_uumean
  logical :: lcalc_uumeanx=.false.,lcalc_uumeanxz=.false.
  logical :: lforcing_cont_uu=.false.
  logical :: lcoriolis_xdep=.false.
  logical :: lno_meridional_flow=.false.
  logical :: lrotation_xaxis=.false.
  logical :: lpropagate_borderuu=.true.
  logical :: lgradu_as_aux=.false.
  logical :: lno_radial_advection=.false.
  logical :: lfargoadvection_as_shift=.true.
  logical :: lhelmholtz_decomp=.false.
  character (len=labellen) :: uuprof='nothing'
!
!  Parameters for interior boundary conditions.
!
  character (len=labellen) :: interior_bc_hydro_profile='nothing'
  logical :: lhydro_bc_interior=.false.
  real :: z1_interior_bc_hydro=0.,kz_analysis=1.
  real :: Shearx=0., rescale_uu=0.
  real :: Ra=0.0, Pr=0.0 ! Boussinesq approximation
!
  namelist /hydro_run_pars/ &
      nu, inituu, ampluu, kz_uu, ampl1_diffrot, ampl2_diffrot, uuprof, &
      xexp_diffrot, kx_diffrot, kz_diffrot, kz_analysis, phase_diffrot, ampl_wind, &
      lreinitialize_uu, lremove_mean_momenta, lremove_mean_flow, &
      lremove_uumeanxy,lremove_uumeanz,lremove_uumeanz_horizontal, &
      lupw_uu, othresh, &
      othresh_per_orms, borderuu, lpressuregradient_gas, &
      ladvection_velocity, &
      utop, ubot, omega_out, omega_in, lprecession, omega_precession, &
      alpha_precession, lshear_rateofstrain, r_omega, w_omega, &
      lalways_use_gij_etc, amp_centforce, &
      lcalc_uumean,lcalc_uumeanx,lcalc_uumeanxy,lcalc_uumeanxz,lcalc_uumeanz, &
      lforcing_cont_uu, width_ff_uu, x1_ff_uu, x2_ff_uu, &
      loo_as_aux, luut_as_aux, loot_as_aux, loutest, &
      interior_bc_hydro_profile, lhydro_bc_interior, z1_interior_bc_hydro, &
      ekman_friction, &
      ampl_forc, k_forc, w_forc, x_forc, dx_forc, ampl_fcont_uu, &
      lno_meridional_flow, lrotation_xaxis, k_diffrot,Shearx, rescale_uu, &
      hydro_xaver_range, Ra, Pr, llinearized_hydro, lremove_mean_angmom, &
      lpropagate_borderuu, hydro_zaver_range, &
      uzjet, mean_momentum, lshear_in_coriolis, &
      w_sldchar_hyd, uphi_rbot, uphi_rtop, uphi_step_width, &
      lno_radial_advection, lfargoadvection_as_shift, lhelmholtz_decomp
!
!  Diagnostic variables (need to be consistent with reset list below).
!
  integer :: idiag_u2tm=0       ! DIAG_DOC: $\left<\uv(t)\cdot\int_0^t\uv(t')
                                ! DIAG_DOC:   dt'\right>$
  integer :: idiag_fkinzm=0     ! DIAG_DOC: $\left<{1\over2} \varrho\uv^2 u_z\right>$
  integer :: idiag_u2m=0        ! DIAG_DOC: $\left<\uv^2\right>$
  integer :: idiag_um2=0        ! DIAG_DOC:
  integer :: idiag_uxpt=0       ! DIAG_DOC: $u_x(x_1,y_1,z_1,t)$
  integer :: idiag_uypt=0       ! DIAG_DOC: $u_y(x_1,y_1,z_1,t)$
  integer :: idiag_uzpt=0       ! DIAG_DOC: $u_z(x_1,y_1,z_1,t)$
  integer :: idiag_uxp2=0       ! DIAG_DOC: $u_x(x_2,y_2,z_2,t)$
  integer :: idiag_uyp2=0       ! DIAG_DOC: $u_y(x_2,y_2,z_2,t)$
  integer :: idiag_uzp2=0       ! DIAG_DOC: $u_z(x_2,y_2,z_2,t)$
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
  integer :: idiag_ux2m=0       ! DIAG_DOC: $\left<u_x^2\right>$
  integer :: idiag_uy2m=0       ! DIAG_DOC: $\left<u_y^2\right>$
  integer :: idiag_uz2m=0       ! DIAG_DOC: $\left<u_z^2\right>$
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
  integer :: idiag_ursphmphi=0  ! PHIAVG_DOC: $\left<u_r\right>_\varphi$
  integer :: idiag_uthmphi=0    ! PHIAVG_DOC: $\left<u_\vartheta\right>_\varphi$
  ! For the manual: uumphi      ! PHIAVG_DOC: shorthand for \var{urmphi},
                                ! PHIAVG_DOC: \var{upmphi} and \var{uzmphi}
                                ! PHIAVG_DOC: together
  ! For the manual: uusphmphi   ! PHIAVG_DOC: shorthand for \var{ursphmphi},
                                ! PHIAVG_DOC: \var{uthmphi} and \var{upmphi}
                                ! PHIAVG_DOC: together
  integer :: idiag_u2mphi=0     ! PHIAVG_DOC: $\left<\uv^2\right>_\varphi$
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
  integer :: idiag_oum=0        ! DIAG_DOC: $\left<\boldsymbol{\omega}
                                ! DIAG_DOC:   \cdot\uv\right>$
  integer :: idiag_ou_int=0     ! DIAG_DOC: $\int_V\boldsymbol{\omega}\cdot\uv\,dV$
  integer :: idiag_fum=0        ! DIAG_DOC: $\left<\fv\cdot\uv\right>$
  integer :: idiag_odel2um=0    ! DIAG_DOC: $\left<\boldsymbol{\omega}\nabla^2\uv\right>$
  integer :: idiag_o2m=0        ! DIAG_DOC: $\left<\boldsymbol{\omega}^2\right>
                                ! DIAG_DOC:   \equiv \left<(\curl\uv)^2\right>$
  integer :: idiag_orms=0       ! DIAG_DOC: $\left<\boldsymbol{\omega}^2\right>^{1/2}$
  integer :: idiag_omax=0       ! DIAG_DOC: $\max(|\boldsymbol{\omega}|)$
  integer :: idiag_ox2m=0       ! DIAG_DOC: $\left<\omega_x^2\right>$
  integer :: idiag_oy2m=0       ! DIAG_DOC: $\left<\omega_y^2\right>$
  integer :: idiag_oz2m=0       ! DIAG_DOC: $\left<\omega_z^2\right>$
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
  integer :: idiag_oumphi=0     ! DIAG_DOC: $\left<\omv\cdot\uv\right>_\varphi$
  integer :: idiag_ozmphi=0     ! DIAG_DOC:
  integer :: idiag_ormr=0       ! DIAG_DOC:
  integer :: idiag_opmr=0       ! DIAG_DOC:
  integer :: idiag_ozmr=0       ! DIAG_DOC:
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
  integer :: idiag_uxmz=0       ! XYAVG_DOC: $\left< u_x \right>_{xy}$
                                ! XYAVG_DOC:   \quad(horiz. averaged $x$
                                ! XYAVG_DOC:   velocity)
  integer :: idiag_uymz=0       ! XYAVG_DOC: $\left< u_y \right>_{xy}$
  integer :: idiag_uzmz=0       ! XYAVG_DOC: $\left< u_z \right>_{xy}$
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
  integer :: idiag_oxdivu2mz=0  ! XYAVG_DOC: $\left<(\omega_x nabla\cdot\uv)^2\right>_{xy}$
  integer :: idiag_oydivu2mz=0  ! XYAVG_DOC: $\left<(\omega_y\nabla\cdot\uv)^2\right>_{xy}$
  integer :: idiag_ozdivu2mz=0  ! XYAVG_DOC: $\left<(\omega_z\nabla\cdot\uv)^2\right>_{xy}$
  integer :: idiag_u3u21mz=0    ! XYAVG_DOC:
  integer :: idiag_u1u32mz=0    ! XYAVG_DOC:
  integer :: idiag_u2u13mz=0    ! XYAVG_DOC:
  integer :: idiag_u2u31mz=0    ! XYAVG_DOC:
  integer :: idiag_u3u12mz=0    ! XYAVG_DOC:
  integer :: idiag_u1u23mz=0    ! XYAVG_DOC:
  integer :: idiag_accpowzmz=0  ! XYAVG_DOC: $\left<(u_z Du_z/Dt)^2\right>_{xy}$
  integer :: idiag_accpowzupmz=0  ! XYAVG_DOC: $\left<(u_z Du_z/Dt)^2\right>_{xy+}$
  integer :: idiag_accpowzdownmz=0! XYAVG_DOC: $\left<(u_z Du_z/Dt)^2\right>_{xy-}$
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
  integer :: idiag_uxuymxy=0    ! ZAVG_DOC: $\left< u_x u_y \right>_{z}$
  integer :: idiag_uxuzmxy=0    ! ZAVG_DOC: $\left< u_x u_z \right>_{z}$
  integer :: idiag_uyuzmxy=0    ! ZAVG_DOC: $\left< u_y u_z \right>_{z}$
  integer :: idiag_oxmxy=0      ! ZAVG_DOC: $\left< \omega_x \right>_{z}$
  integer :: idiag_oymxy=0      ! ZAVG_DOC: $\left< \omega_y \right>_{z}$
  integer :: idiag_ozmxy=0      ! ZAVG_DOC: $\left< \omega_z \right>_{z}$
  integer :: idiag_oumxy=0      ! ZAVG_DOC: $\left<\boldsymbol{\omega}
                                ! ZAVG_DOC: \cdot\uv\right>_{z}$
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
  integer :: idiag_nshift=0
!
!  Video data.
!
  integer :: ivid_oo=0, ivid_o2=0, ivid_divu=0, ivid_u2=0, ivid_Ma2=0
!
!  Auxiliary variables
!
  real, dimension(:,:), pointer :: reference_state
!
  contains
!***********************************************************************
    subroutine register_hydro
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
! 26-dec-18/axel: adapted from hydro
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
!  indices to access uu
!
      call farray_register_pde('varphi',iphiuu)
      call farray_register_auxiliary('uu',iuu,vector=3)
      iux = iuu; iuy = iuu+1; iuz = iuu+2
!
!!  For Helmholtz decomposition of uu the potential of the curl-free part is registered  as an auxiliary.
!!
!      if (lhelmholtz_decomp) then
!        if (dsnap_down==0.) call fatal_error('register_hydro','Helmholtz decomposition requires dsnap_down>0')
!        call farray_register_auxiliary('phiuu',iphiuu)
!      endif
!
!  Share lpressuregradient_gas so the entropy module knows whether to apply
!  pressure gradient or not. But hydro wants pressure gradient only when
!  the density is computed, i.e. not with lboussinesq nor lanelastic.
!
      if  (.not.ldensity.or.lanelastic) lpressuregradient_gas=.false.
      call put_shared_variable('lpressuregradient_gas', &
          lpressuregradient_gas,caller='register_hydro')
!
!!  Special settings for lboussinesq.
!!
!      if  (lboussinesq) then
!        PrRa=Pr*Ra
!        call put_shared_variable('PrRa',PrRa)
!        call put_shared_variable('Pr',Pr)
!      endif
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Writing files for use with IDL.
!
      if (lroot) then
!       if (maux == 0) then
!         if (nvar < mvar) write(4,*) ',uu $'
!         if (nvar == mvar) write(4,*) ',uu'
!       else
          write(4,*) ',phiuu $'
          write(4,*) ',uu $'
!       endif
        write(15,*) 'uu = fltarr(mx,my,mz,3)*one'
        write(15,*) 'phiuu = fltarr(mx,my,mz,1)*one'
      endif
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
!  26-dec-18/axel: adapted from hydro
!
      use BorderProfiles, only: request_border_driving
      use Initcond
      use SharedVariables, only: put_shared_variable,get_shared_variable
      use Sub, only: step, erfunc, register_report_aux
      !use Slices_methods, only: alloc_slice_buffers
      use Yinyang_mpi, only: initialize_zaver_yy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: c, s
      integer :: j,myl ! currently unused: nycap
!
!     ! MAUX CONTRIBUTION 3
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!
      if (luut_as_aux) then
        call register_report_aux('uut', iuut, iuxt, iuyt, iuzt)
        ltime_integrals=.true.
      endif
      if (loot_as_aux) then
        call register_report_aux('oot', ioot, ioxt, ioyt, iozt)
        ltime_integrals=.true.
      endif
      if (loo_as_aux) call register_report_aux('oo', ioo, iox, ioy, ioz, communicated=.true.)
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
!  Get the reference state if requested
!
      if (lreference_state) &
        call get_shared_variable('reference_state',reference_state)
!
      lcalc_uumeanz = lcalc_uumeanz .or. lcalc_uumean .or. ltestfield_xz      ! lcalc_uumean for compatibility
!
      lcalc_uumeanxy=lremove_uumeanxy .or. lcalc_uumeanxy .or. ltestfield_xy
!
      if (lcalc_uumeanxy) then
        myl=my
        if (lyinyang) then
          call fatal_error('initialize_hydro','Calculation of z average not implmented for Yin-Yang')
          !call initialize_zaver_yy(myl,nycap)
        endif
        allocate(uumxy(mx,myl,3))
        uumxy=0.0
      endif
!
      if (ivid_oo/=0) then
        !call alloc_slice_buffers(oo_xy,oo_xz,oo_yz,oo_xy2,oo_xy3,oo_xy4,oo_xz2)
        if (lwrite_slice_xy .and..not.allocated(oo_xy) ) allocate(oo_xy (nx,ny,3))
        if (lwrite_slice_xz .and..not.allocated(oo_xz) ) allocate(oo_xz (nx,nz,3))
        if (lwrite_slice_yz .and..not.allocated(oo_yz) ) allocate(oo_yz (ny,nz,3))
        if (lwrite_slice_xy2.and..not.allocated(oo_xy2)) allocate(oo_xy2(nx,ny,3))
        if (lwrite_slice_xy3.and..not.allocated(oo_xy3)) allocate(oo_xy3(nx,ny,3))
        if (lwrite_slice_xy4.and..not.allocated(oo_xy4)) allocate(oo_xy4(nx,ny,3))
        if (lwrite_slice_xz2.and..not.allocated(oo_xz2)) allocate(oo_xz2(nx,nz,3))
      endif
      if (ivid_o2/=0) then
        !call alloc_slice_buffers(o2_xy,o2_xz,o2_yz,o2_xy2,o2_xy3,o2_xy4,o2_xz2)
        if (lwrite_slice_xy .and..not.allocated(o2_xy) ) allocate(o2_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(o2_xz) ) allocate(o2_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(o2_yz) ) allocate(o2_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(o2_xy2)) allocate(o2_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(o2_xy3)) allocate(o2_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(o2_xy4)) allocate(o2_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(o2_xz2)) allocate(o2_xz2(nx,nz))
      endif
      if (ivid_u2/=0) then
        !call alloc_slice_buffers(u2_xy,u2_xz,u2_yz,u2_xy2,u2_xy3,u2_xy4,u2_xz2)
        if (lwrite_slice_xy .and..not.allocated(u2_xy) ) allocate(u2_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(u2_xz) ) allocate(u2_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(u2_yz) ) allocate(u2_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(u2_xy2)) allocate(u2_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(u2_xy3)) allocate(u2_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(u2_xy4)) allocate(u2_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(u2_xz2)) allocate(u2_xz2(nx,nz))
      endif
      if (ivid_divu/=0) then
        !call alloc_slice_buffers(divu_xy,divu_xz,divu_yz,divu_xy2,divu_xy3,divu_xy4,divu_xz2)
        if (lwrite_slice_xy .and..not.allocated(divu_xy) ) allocate(divu_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(divu_xz) ) allocate(divu_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(divu_yz) ) allocate(divu_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(divu_xy2)) allocate(divu_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(divu_xy3)) allocate(divu_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(divu_xy4)) allocate(divu_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(divu_xz2)) allocate(divu_xz2(nx,nz))
      endif
      if (ivid_Ma2 /=0) then
        !call alloc_slice_buffers(mach_xy,mach_xz,mach_yz,mach_xy2,mach_xy3,mach_xy4,mach_xz2)
        if (lwrite_slice_xy .and..not.allocated(mach_xy) ) allocate(mach_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(mach_xz) ) allocate(mach_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(mach_yz) ) allocate(mach_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(mach_xy2)) allocate(mach_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(mach_xy3)) allocate(mach_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(mach_xy4)) allocate(mach_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(mach_xz2)) allocate(mach_xz2(nx,nz))
      endif
!
!  give warning if orms is not set in prints.in
!
      if (othresh_per_orms/=0..and.idiag_orms==0) then
        if (lroot) &
          print*,'calc_othresh: need to set orms in print.in to get othresh.'
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
!  do xy-averaged mean field for each component
!
      if (lcalc_uumeanz) then
        call fatal_error('hydro_potential','duu_dt')
      endif
!
!  do yz-averaged mean field for each component
!
      if (lcalc_uumeanx) then
        call fatal_error('hydro_potential','duu_dt')
      endif
!
!  Do mean 2D field in (x,y)-plane for each component
!
      if (lcalc_uumeanxy) then
        call fatal_error('hydro_potential','duu_dt')
      endif
!
!  Do mean 2D field in (x,z)-plane for each component
!
      if (lcalc_uumeanxz) then
        call fatal_error('hydro_potential','duu_dt')
      endif
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
      use Boundcond, only:update_ghosts
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
      real, dimension (nx) :: r,p1,tmp,prof,xc0,yc0
      real, dimension (:,:), allocatable :: yz
      real :: kabs,crit,eta_sigma,tmp0
      real :: a2, rr2, wall_smoothing
      real :: dis, xold,yold,uprof, factx, factz, sph, sph_har_der, der
      integer :: j,i,l,ixy,ix,iy,iz,iz0,iyz
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
          f(:,:,:,iphiuu)=0.
        case ('gaussian-noise'); call gaunoise(ampluu(j),f,iphiuu)
        case ('sinwave-phase')
          call sinwave_phase(f,iphiuu,ampl_ux(j),kx_ux(j),ky_ux(j),kz_ux(j),phase_ux(j))
        case ('coswave-phase')
          call coswave_phase(f,iphiuu,ampl_ux(j),kx_ux(j),ky_ux(j),kz_ux(j),phase_ux(j))
        case ('sound-wave', '11')
!
!  Burgers shock
!
          if (lroot) print*,'init_uu: Burgers shock'
          prof=-ampluu(j)*tanh(.5*x(l1:l2)/widthuu)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iphiuu)=prof
          enddo; enddo
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
      if (lparticles_lyapunov .or. lparticles_caustics) & 
        lpenc_requested(i_uij) = .true.
!AXEL if (ladvection_velocity) then
!       if (lweno_transport) then
!         lpenc_requested(i_uu)=.true.
!         lpenc_requested(i_rho1)=.true.
!         lpenc_requested(i_transprho)=.true.
!         lpenc_requested(i_transpurho)=.true.
!       endif
!     endif
!     if (lprecession) lpenc_requested(i_rr)=.true.
!     if (ldt.or.(ekman_friction/=0)) lpenc_requested(i_uu)=.true.
!
      if (lisotropic_advection) lpenc_requested(i_u2)=.true.
!
!  The following only play a role if lsphere_in_a_box
!  and there is a buoyancy term. (If only hydro, there is no buoyancy.)
!
!     if (lboussinesq.and.ltemperature.and.lsphere_in_a_box) then
!       lpenc_requested(i_evr)=.true.
!       lpenc_requested(i_r_mn)=.true.
!     endif
!
!  video pencils
!
      if (lwrite_slices) then
        if (ivid_oo  /=0) lpenc_video(i_oo)=.true.
        if (ivid_o2  /=0) lpenc_video(i_o2)=.true.
        if (ivid_divu/=0) lpenc_video(i_divu)=.true.
        if (ivid_u2  /=0) lpenc_video(i_u2)=.true.
        if (ivid_Ma2 /=0) lpenc_video(i_Ma2)=.true.
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
          idiag_oxm /=0 .or. idiag_oym /=0 .or. idiag_ozm /=0 .or. &
          idiag_oxoym/=0 .or. idiag_oxozm/=0 .or. idiag_oyozm/=0 .or. &
          idiag_oxuzxm/=0 .or. idiag_oyuzym/=0 .or. &
          idiag_quxom/=0) &
          lpenc_diagnos(i_oo)=.true.
      if (idiag_orms/=0 .or. idiag_omax/=0 .or. idiag_o2m/=0 .or. &
          idiag_ormsh/=0 .or. idiag_o2mz/=0 ) lpenc_diagnos(i_o2)=.true.
      if (idiag_q2m/=0 .or. idiag_qrms/=0 .or. idiag_qmax/=0 .or. &
          idiag_qfm/=0 .or. idiag_qom/=0 .or. idiag_quxom/=0) &
          lpenc_diagnos(i_curlo)=.true.
      if (idiag_divu2m/=0 .or. idiag_divu2mz/=0 .or. idiag_divru2mz/=0 .or. &
          idiag_divum/=0 .or. idiag_rdivum/=0) lpenc_diagnos(i_divu)=.true.
      if (idiag_rdivum/=0 .or. idiag_fkinzm/=0 .or. idiag_ruzupmz/=0 .or. idiag_ruzdownmz/=0) &
          lpenc_diagnos(i_rho)=.true.
      if (idiag_gdivu2m/=0) lpenc_diagnos(i_graddivu)=.true.
      if (idiag_oum/=0 .or. idiag_oumx/=0.or.idiag_oumy/=0.or.idiag_oumz/=0 .or. &
           idiag_oumh/=0 .or. idiag_ou_int/=0 ) lpenc_diagnos(i_ou)=.true.
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
      if (idiag_accpowzmz/=0 .or. idiag_accpowzupmz/=0 .or. idiag_accpowzdownmz/=0) then
        lpenc_diagnos(i_fpres)=.true.; lpenc_diagnos(i_fvisc)=.true.
        if (lgrav) lpenc_diagnos(i_gg)=.true.
      endif
      if (idiag_urms/=0 .or. idiag_durms/=0 .or. &
          idiag_umax/=0 .or. idiag_rumax/=0 .or. &
          idiag_fkinzm/=0 .or. idiag_u2m/=0 .or. idiag_um2/=0 .or. idiag_u2mz/=0 .or. &
          idiag_urmsh/=0 .or. idiag_urmsx/=0 .or. idiag_urmsz/=0) lpenc_diagnos(i_u2)=.true.
      if (idiag_duxdzma/=0 .or. idiag_duydzma/=0 .or. lgradu_as_aux) lpenc_diagnos(i_uij)=.true.
      if (idiag_fmasszmz/=0 .or. idiag_ruxuym/=0 .or. &
          idiag_ruxm/=0 .or. idiag_ruym/=0 .or. idiag_ruzm/=0 .or. &
          idiag_ruxuzm/=0 .or. idiag_ruyuzm/=0 .or. &
          idiag_ruxtot/=0) lpenc_diagnos(i_rho)=.true.
      if (idiag_ruxmz/=0 .or. idiag_ruymz/=0 .or. idiag_ruzmz/=0 .or. &
          idiag_rux2mz/=0 .or. idiag_ruy2mz/=0 .or. idiag_ruz2mz/=0 .or. &
          idiag_ruxmx/=0 .or. idiag_ruymx/=0 .or. idiag_ruzmx/=0 .or. &
          idiag_ruxuymz/=0 .or. idiag_ruxuzmz/=0 .or. idiag_ruyuzmz/=0 .or. &
          idiag_ruxuy2mz/=0 .or. idiag_ruxuz2mz/=0 .or. idiag_ruyuz2mz/=0) &
          lpenc_diagnos(i_rho)=.true.
      if (idiag_rux2mx /= 0 .or. idiag_ruy2mx /= 0 .or. idiag_ruz2mx /= 0 .or. &
          idiag_ruxuymx /= 0 .or. idiag_ruxuzmx /= 0 .or. idiag_ruyuzmx /= 0) lpenc_diagnos(i_rho) = .true.
      if (idiag_ormr/=0 .or. idiag_opmr/=0 .or. idiag_ozmr/=0) &
          lpenc_diagnos(i_oo)=.true.
      if (idiag_oxmxy/=0 .or. idiag_oymxy/=0 .or. idiag_ozmxy/=0 .or. &
          idiag_oxmz/=0 .or. idiag_oymz/=0 .or. idiag_ozmz/=0 .or. &
          idiag_ox2mz/=0 .or. idiag_oy2mz/=0 .or. idiag_oz2mz/=0 ) &
          lpenc_diagnos2d(i_oo)=.true.
      if (idiag_totangmom/=0 ) lpenc_diagnos(i_rcyl_mn)=.true.
      if (idiag_urmr/=0 .or.  idiag_ormr/=0 .or. idiag_urmphi/=0) then
        lpenc_diagnos(i_pomx)=.true.
        lpenc_diagnos(i_pomy)=.true.
      endif
      if (idiag_ursphmphi/=0) lpenc_diagnos2d(i_evr)=.true.
      if (idiag_uthmphi/=0) lpenc_diagnos2d(i_evth)=.true.
      if (idiag_upmr/=0 .or. idiag_opmr/=0 .or. idiag_upmphi/=0) then
        lpenc_diagnos(i_phix)=.true.
        lpenc_diagnos(i_phiy)=.true.
      endif
      if (idiag_EEK/=0 .or. idiag_ekin/=0 .or. idiag_ekintot/=0 .or. idiag_fkinzmz/=0 .or. &
          idiag_ekinmx /= 0 .or. idiag_ekinmz/=0 .or. idiag_fkinxmx/=0) then
        lpenc_diagnos(i_ekin)=.true.
      endif
      if (idiag_fkinxmxy/=0 .or. idiag_fkinymxy/=0) then
        lpenc_diagnos2d(i_uu)=.true.
        lpenc_diagnos2d(i_ekin)=.true.
      endif
      if (idiag_uxmxy/=0 .or. idiag_uymxy/=0 .or. idiag_uzmxy/=0 .or. &
          idiag_ux2mxy/=0 .or. idiag_uy2mxy/=0 .or. idiag_uz2mxy/=0 .or. &
          idiag_ruxmxy/=0 .or. idiag_ruymxy/=0 .or. idiag_ruzmxy/=0 .or.  &
          idiag_uxuymxy/=0 .or. idiag_uxuzmxy/=0 .or. idiag_uyuzmxy/=0 .or. &
          idiag_ruxuymxy/=0 .or. idiag_ruxuzmxy/=0 .or. &
          idiag_ruyuzmxy/=0) then
        lpenc_diagnos2d(i_uu)=.true.
      endif
      if (idiag_ruxmxy/=0 .or. idiag_ruymxy/=0 .or. idiag_ruzmxy/=0 .or. &
          idiag_ruxuymxy/=0 .or. idiag_ruxuzmxy/=0 .or. idiag_ruyuzmxy/=0) &
        lpenc_diagnos2d(i_rho)=.true.
      if (idiag_Remz/=0) then
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
      if (lpencil_in(i_ou)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_oo)=.true.
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
      real, dimension (mx,my,mz,mfarray),intent(IN) :: f
      type (pencil_case),                intent(OUT):: p
      logical, dimension(:),             intent(IN) :: lpenc_loc
!
      integer :: kk

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
      real, dimension (mx,my,mz,mfarray),intent(IN) :: f
      type (pencil_case),                intent(OUT):: p
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
      real, dimension (nx) :: tmp
      real, dimension (nx,3) :: tmp3
      integer :: i, j, ju, ij, jj, kk, jk
!
      intent(in) :: lpenc_loc
      intent(out):: p

! uu
      if (lpenc_loc(i_uu)) then
        call grad(f,iphiuu,p%uu)
        if (iux/=0) f(l1:l2,m,n,iux:iuz)=p%uu
      endif
! u2
      if (lpenc_loc(i_u2)) call dot2_mn(p%uu,p%u2)
! uij
      if (lpenc_loc(i_uij)) then
        call g2ij(f,iphiuu,p%uij)
        if (lparticles_lyapunov .or. lparticles_caustics) then
          jk=0
          do jj=1,3; do kk=1,3
            jk=jk+1
            f(l1:l2,m,n,iguij+jk-1) = p%uij(:,jj,kk)
          enddo;enddo
        endif
      endif
!      if (.not.lpenc_loc_check_at_work) then
!        write(*,*) 'uurad,rad',p%uij(1:6,1,1)
!      endif
!
!  if gradu is to be stored as auxiliary, then we store it now
!
      if (lgradu_as_aux) then
        ij=igradu-1
        do i=1,3
          do j=1,3
            ij=ij+1
            f(l1:l2,m,n,ij) = p%uij(:,i,j)
          enddo
        enddo
      endif
! divu
      if (lpenc_loc(i_divu)) call div_mn(p%uij,p%divu,p%uu)
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
! o2
      if (lpenc_loc(i_o2)) call dot2_mn(p%oo,p%o2)
! ou
      if (lpenc_loc(i_ou)) call dot_mn(p%oo,p%uu,p%ou)
! Useful to debug forcing - Dhruba
      if (loutest.and.lpenc_loc(i_ou))then
!      write(*,*) lpenc_loc(i_ou)
        outest = minval(p%ou)
        if (outest<(-1.0d-8))then
          write(*,*) m,n,outest,maxval(p%ou),lpenc_loc(i_ou)
          write(*,*)'WARNING : hydro:ou has different sign than relhel'
        else
        endif
      else
      endif
! ugu
      if (lpenc_loc(i_ugu)) then
        if (headtt.and.lupw_uu) print *,'calc_pencils_hydro: upwinding advection term'
        call fatal_error('calc_pencils_hydro_nonlinear:','does not calculate ugu pencil')
      endif
! ugu2
      if (lpenc_loc(i_ugu2)) call fatal_error('calc_pencils_hydro_nonlinear:','does not calculate ugu2 pencil')
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
          p%uu_advec(:,2)=p%uu(:,2)-uu_average_cyl(:,n-nghost)
          p%uu_advec(:,3)=p%uu(:,3)
        elseif (lspherical_coords) then
          p%uu_advec(:,2)=p%uu(:,2)
          p%uu_advec(:,3)=p%uu(:,3)-uu_average_sph(:,m-nghost)
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
      if (lpenc_loc(i_o2)) call fatal_error('calc_pencils_hydro_linearized:','does not calculate o2 pencil')
! ou
      if (lpenc_loc(i_ou)) call fatal_error('calc_pencils_hydro_linearized:','does not calculate ou pencil')
! ugu
      if (lpenc_loc(i_ugu)) then
        if (headtt.and.lupw_uu) print *,'calc_pencils_hydro: upwinding advection term'
!        call u_dot_grad(f,iuu,p%uij,p%uu,p%ugu,UPWIND=lupw_uu)
        call u_dot_grad(f,iuu0,p%u0ij,p%uu,ugu0,UPWIND=lupw_uu)
        call u_dot_grad(f,iuu,p%uij,p%uu0,u0gu,UPWIND=lupw_uu)
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
    subroutine calc_diagnostics_hydro(f,p)

      real, dimension(:,:,:,:) :: f
      type(pencil_case), intent(in) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)

    endsubroutine calc_diagnostics_hydro
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
      real, dimension (nx,nz) :: fsum_tmp_cyl
      real, dimension (nx,ny) :: fsum_tmp_sph
      real, dimension (nx) :: uphi
      real :: nygrid1,nzgrid1
      integer :: nnghost,mnghost
!
!  Remove mean momenta or mean flows if desired.
!  Useful to avoid unphysical winds, for example in shearing box simulations.
!
      if (lremove_mean_momenta) then
        call remove_mean_momenta(f,iux)
      else
        if (lremove_mean_flow) call remove_mean_flow(f,iux)
        !if (lremove_mean_flow) call remove_mean_value(f,iux,iuz)  !(could use this one)
        if (lremove_mean_angmom) call remove_mean_angmom(f,iuz)
      endif
!
!  Calculate the vorticity field if required.
!
      getoo: if (ioo /= 0) then
        nloop: do n = n1, n2
          mloop: do m = m1, m2
            call curl(f, iux, pv)
            f(l1:l2,m,n,iox:ioz) = pv
          enddo mloop
        enddo nloop
      endif getoo
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
!  calculate dvarphi/dt = f - h + nu*divu
!
!  26-dec-18/axel: adapted from hydro
!
      use Diagnostics
      use Special, only: special_calc_hydro
      use Sub, only: vecout, dot, dot2, identify_bcs, cross, multsv_mn_add
      use General, only: transform_thph_yy, notanumber
      use Slices_methods, only: store_slices
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: f,df
!
      real, dimension (nx,3) :: curlru,uxo
      real, dimension (nx) :: space_part_re,space_part_im,u2t,uot,out,fu
      real, dimension (nx) :: odel2um,curlru2,uref,curlo2,qo,quxo,graddivu2
      real, dimension (nx) :: uus,ftot,Fmax,advec_uu
      real :: kx
      integer :: j, ju, k
      integer, dimension(nz), save :: nuzup=0, nuzdown=0, nruzup=0, nruzdown=0
!
      Fmax=tini
!
!  Identify module and boundary conditions.
!
      call timing('duu_dt','entered',mnloop=.true.)
      if (headtt.or.ldebug) print*,'duu_dt: SOLVE'
      if (headtt) then
        call identify_bcs('phiuu',iphiuu)
      endif
!
!  Advection term.
!
      df(l1:l2,m,n,iphiuu)=df(l1:l2,m,n,iphiuu)-.5*p%u2+nu*p%divu
!
!  Interface for your personal subroutines calls
!
      if (lspecial) call special_calc_hydro(f,df,p)
!
!  ``uu/dx'' for timestep
!
      if (lfirst.and.ldt.and.ladvection_velocity) then
        if (lmaximal_cdt) then
          advec_uu=max(abs(p%uu(:,1))*dline_1(:,1),&
                       abs(p%uu(:,2))*dline_1(:,2),&
                       abs(p%uu(:,3))*dline_1(:,3))
        endif
      else
        advec_uu=0.0
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
         if (lfirst.and.ldt) then
            if ((nxgrid==1).or.(nygrid==1).or.(nzgrid==1)) &
                 advec_uu=sqrt(p%u2*dxyz_2)
         endif
      endif
      if (lfirst.and.ldt) then
        if (notanumber(advec_uu)) print*, 'advec_uu   =',advec_uu
        maxadvec=maxadvec+advec_uu
        if (headtt.or.ldebug) print*,'duu_dt: max(advec_uu) =',maxval(advec_uu)
      endif
!
!  Ekman Friction, used only in two dimensional runs.
!
      if (ekman_friction/=0) &
        df(l1:l2,m,n,iphiuu)=df(l1:l2,m,n,iphiuu)-ekman_friction*f(l1:l2,m,n,iphiuu)
!
!  Add possibility of forcing that is not delta-correlated in time.
!
      if (lforcing_cont_uu) &
        df(l1:l2,m,n,iphiuu)=df(l1:l2,m,n,iphiuu)+ &
            ampl_fcont_uu*p%fcont(:,1,1)
!
!  store slices for output in wvid in run.f90
!  This must be done outside the diagnostics loop (accessed at different times).
!
      if (lvideo.and.lfirst) then
        if (ivid_divu/=0) call store_slices(p%divu,divu_xy,divu_xz,divu_yz,divu_xy2,divu_xy3,divu_xy4,divu_xz2)
        if (ivid_oo  /=0) call store_slices(p%oo,oo_xy,oo_xz,oo_yz,oo_xy2,oo_xy3,oo_xy4,oo_xz2)
        call vecout(41,trim(directory)//'/ovec',p%oo,othresh,novec)

        if (ivid_u2  /=0) call store_slices(p%u2,u2_xy,u2_xz,u2_yz,u2_xy2,u2_xy3,u2_xy4,u2_xz2)
        if (ivid_o2  /=0) call store_slices(p%o2,o2_xy,o2_xz,o2_yz,o2_xy2,o2_xy3,o2_xy4,o2_xz2)
        if (othresh_per_orms/=0) call calc_othresh
        if (ivid_Ma2 /=0) call store_slices(p%Ma2,mach_xy,mach_xz,mach_yz,mach_xy2,mach_xy3,mach_xy4,mach_xz2)
      endif
!
!  Calculate maxima and rms values for diagnostic purposes
!
      call timing('duu_dt','just before ldiagnos',mnloop=.true.)
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'duu_dt: Calculate maxima and rms values...'
        if (idiag_dtu/=0) call max_mn_name(advec_uu/cdt,idiag_dtu,l_dt=.true.)
        if (idiag_urms/=0) call sum_mn_name(p%u2,idiag_urms,lsqrt=.true.)
        if (idiag_durms/=0) then
          uref=ampluu(1)*cos(kx_uu*x(l1:l2))
          call sum_mn_name(p%u2-2.*p%uu(:,2)*uref+uref**2,idiag_durms)
        endif
        endif
!
        if (idiag_urmsx/=0) call sum_mn_name(p%u2*xmask_hyd,idiag_urmsx,lsqrt=.true.)
        if (idiag_urmsz/=0) call sum_mn_name(p%u2*zmask_hyd(n-n1+1),idiag_urmsz,lsqrt=.true.)
        if (idiag_umax /=0) call max_mn_name(p%u2,idiag_umax,lsqrt=.true.)
        if (idiag_umin /=0) call max_mn_name(-sqrt(p%u2),idiag_umin,lneg=.true.)
        if (idiag_uxrms/=0) &
            call sum_mn_name(p%uu(:,1)**2,idiag_uxrms,lsqrt=.true.)
        if (idiag_uzrms/=0) &
            call sum_mn_name(p%uu(:,3)**2,idiag_uzrms,lsqrt=.true.)
        if (idiag_uzrmaxs/=0) &
            call max_mn_name(p%uu(:,3)**2,idiag_uzrmaxs,lsqrt=.true.)
        if (idiag_uxmin/=0) call max_mn_name(-p%uu(:,1),idiag_uxmin,lneg=.true.)
        if (idiag_uymin/=0) call max_mn_name(-p%uu(:,2),idiag_uymin,lneg=.true.)
        if (idiag_uzmin/=0) call max_mn_name(-p%uu(:,3),idiag_uzmin,lneg=.true.)
        if (idiag_uxmax/=0) call max_mn_name(p%uu(:,1),idiag_uxmax)
        if (idiag_uymax/=0) call max_mn_name(p%uu(:,2),idiag_uymax)
        if (idiag_uzmax/=0) call max_mn_name(p%uu(:,3),idiag_uzmax)
        if (idiag_rumax/=0) call max_mn_name(p%u2*p%rho**2,idiag_rumax,lsqrt=.true.)
        if (idiag_dudx/=0) then
          call sum_mn_name(p%uij(:,1,1),idiag_dudx)
        endif
        if (idiag_fkinzm/=0) call sum_mn_name(.5*p%rho*p%u2*p%uu(:,3),idiag_fkinzm)
        if (idiag_u2m/=0)     call sum_mn_name(p%u2,idiag_u2m)
        if (idiag_um2/=0)     call max_mn_name(p%u2,idiag_um2)
        if (idiag_divum/=0)   call sum_mn_name(p%divu,idiag_divum)
        if (idiag_rdivum/=0)  call sum_mn_name(p%rho*p%divu,idiag_rdivum)
        if (idiag_divu2m/=0)  call sum_mn_name(p%divu**2,idiag_divu2m)
        if (idiag_gdivu2m/=0) then
          call dot2(p%graddivu,graddivu2)
          call sum_mn_name(graddivu2,idiag_gdivu2m)
        endif
        if (idiag_divrhourms/=0) call sum_mn_name((p%rho*p%divu+p%ugrho)**2,idiag_divrhourms,lsqrt=.true.)
        if (idiag_divrhoumax/=0) call max_mn_name(p%rho*p%divu+p%ugrho,idiag_divrhoumax)
        if (idiag_uxm/=0)     call sum_mn_name(p%uu(:,1),idiag_uxm)
        if (idiag_uym/=0)     call sum_mn_name(p%uu(:,2),idiag_uym)
        if (idiag_uzm/=0)     call sum_mn_name(p%uu(:,3),idiag_uzm)
        if (idiag_ux2m/=0)    call sum_mn_name(p%uu(:,1)**2,idiag_ux2m)
        if (idiag_uy2m/=0)    call sum_mn_name(p%uu(:,2)**2,idiag_uy2m)
        if (idiag_uz2m/=0)    call sum_mn_name(p%uu(:,3)**2,idiag_uz2m)
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
        if (idiag_rux2m/=0)   call sum_mn_name(p%rho*p%uu(:,1)**2,idiag_rux2m)
        if (idiag_ruy2m/=0)   call sum_mn_name(p%rho*p%uu(:,2)**2,idiag_ruy2m)
        if (idiag_ruz2m/=0)   call sum_mn_name(p%rho*p%uu(:,3)**2,idiag_ruz2m)
        if (idiag_EEK/=0)     call sum_mn_name(p%ekin,idiag_EEK)
        if (idiag_ekin/=0)    call sum_mn_name(p%ekin,idiag_ekin)
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
!  mean squared velocity and vorticity
!
        call xysum_mn_name_z(p%u2,idiag_u2mz)
        call xysum_mn_name_z(p%o2,idiag_o2mz)
        call xysum_mn_name_z(p%divu**2,idiag_divu2mz)
!
!  mean squared mass flux divergence
!
        call xysum_mn_name_z((p%rho*p%divu+p%ugrho)**2,idiag_divru2mz)
!
!  mean squared curl of mass flux
!
        if (idiag_curlru2mz/=0) then
          call cross(p%grho,p%uu,curlru)
          call multsv_mn_add(p%rho,p%oo,curlru)
          call dot2(curlru,curlru2)
          call xysum_mn_name_z(curlru2,idiag_curlru2mz)
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
!  Things related to vorticity.
!
        if (idiag_ou_int/=0)  call integrate_mn_name(p%ou,idiag_ou_int)
        if (idiag_oum/=0) call sum_mn_name(p%ou,idiag_oum)
        if (idiag_oumh/=0) then
          if (lequatory) call sum_mn_name_halfy(p%ou,idiag_oumh)
          if (lequatorz) call sum_mn_name_halfz(p%ou,idiag_oumh)
          fname(idiag_oumn)=fname_half(idiag_oumh,1)
          fname(idiag_oums)=fname_half(idiag_oumh,2)
          itype_name(idiag_oumn)=ilabel_sum
          itype_name(idiag_oums)=ilabel_sum
        endif
        if (idiag_orms/=0) call sum_mn_name(p%o2,idiag_orms,lsqrt=.true.)
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
        if (idiag_ox2m/=0) call sum_mn_name(p%oo(:,1)**2,idiag_ox2m)
        if (idiag_oy2m/=0) call sum_mn_name(p%oo(:,2)**2,idiag_oy2m)
        if (idiag_oz2m/=0) call sum_mn_name(p%oo(:,3)**2,idiag_oz2m)
        if (idiag_oxm /=0) call sum_mn_name(p%oo(:,1)   ,idiag_oxm)
        if (idiag_oym /=0) call sum_mn_name(p%oo(:,2)   ,idiag_oym)
        if (idiag_ozm /=0) call sum_mn_name(p%oo(:,3)   ,idiag_ozm)
        if (idiag_oxuzxm/=0) call sum_mn_name(p%oo(:,1)*p%uij(:,3,1),idiag_oxuzxm)
        if (idiag_oyuzym/=0) call sum_mn_name(p%oo(:,2)*p%uij(:,3,2),idiag_oyuzym)
        if (idiag_oxoym/=0) call sum_mn_name(p%oo(:,1)*p%oo(:,2),idiag_oxoym)
        if (idiag_oxozm/=0) call sum_mn_name(p%oo(:,1)*p%oo(:,3),idiag_oxozm)
        if (idiag_oyozm/=0) call sum_mn_name(p%oo(:,2)*p%oo(:,3),idiag_oyozm)
!
!  diagnostics involving curlo [ =curl(omega) ]
!
        if (idiag_q2m/=0 .or. idiag_qrms/=0 .or. idiag_qmax/=0 ) then
          call dot2(p%curlo,curlo2)
          if (idiag_q2m/=0)  call sum_mn_name(curlo2,idiag_q2m)
          if (idiag_qrms/=0) call sum_mn_name(curlo2,idiag_qrms,lsqrt=.true.)
          if (idiag_qmax/=0) call max_mn_name(curlo2,idiag_qmax,lsqrt=.true.)
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
        if (idiag_Marms/=0) call sum_mn_name(p%Ma2,idiag_Marms,lsqrt=.true.)
        if (idiag_Mamax/=0) call max_mn_name(p%Ma2,idiag_Mamax,lsqrt=.true.)
!
!  Diagonal components of alpha using FOSA:
!    alp11=<u3*u2,1>-<u2*u3,1>
!    alp22=<u1*u3,2>-<u3*u1,2>
!    alp33=<u2*u1,3>-<u1*u2,3>
!  For fully periodic domains it is sufficient to compute, e.g., only:
!    alp11=<u3*u2,1>,  alp22=<u1*u3,2>,  alp33=<u2*u1,3>
!
        if (idiag_u3u21m/=0) call sum_mn_name(p%u3u21,idiag_u3u21m)
        if (idiag_u1u32m/=0) call sum_mn_name(p%u1u32,idiag_u1u32m)
        if (idiag_u2u13m/=0) call sum_mn_name(p%u2u13,idiag_u2u13m)
        if (idiag_u2u31m/=0) call sum_mn_name(p%u2u31,idiag_u2u31m)
        if (idiag_u3u12m/=0) call sum_mn_name(p%u3u12,idiag_u3u12m)
        if (idiag_u1u23m/=0) call sum_mn_name(p%u1u23,idiag_u1u23m)
!
!  integrate velocity in time, to calculate correlation time later
!
        if (idiag_u2tm/=0) then
          if (iuut==0) call fatal_error("duu_dt","Cannot calculate u2tm if iuut==0")
          call dot(p%uu,f(l1:l2,m,n,iuxt:iuzt),u2t)
          call sum_mn_name(u2t,idiag_u2tm)
        endif
!
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
        call sign_masked_xyaver(p%uu(:,3),idiag_uzupmz,nuzup)
        call sign_masked_xyaver(-p%uu(:,3),idiag_uzdownmz,nuzdown)
        if (idiag_ruzupmz>0) call sign_masked_xyaver(p%rho*p%uu(:,3),idiag_ruzupmz,nruzup)
        if (idiag_ruzdownmz>0) call sign_masked_xyaver(-p%rho*p%uu(:,3),idiag_ruzdownmz,nruzdown)
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
        if (idiag_rux2mx /= 0) call yzsum_mn_name_x(p%rho * p%uu(:,1)**2, idiag_rux2mx)
        if (idiag_ruy2mx /= 0) call yzsum_mn_name_x(p%rho * p%uu(:,2)**2, idiag_ruy2mx)
        if (idiag_ruz2mx /= 0) call yzsum_mn_name_x(p%rho * p%uu(:,3)**2, idiag_ruz2mx)
        if (idiag_ruxuymx/=0) call yzsum_mn_name_x(p%rho*p%uu(:,1)*p%uu(:,2),idiag_ruxuymx)
        if (idiag_ruxuzmx/=0) call yzsum_mn_name_x(p%rho*p%uu(:,1)*p%uu(:,3),idiag_ruxuzmx)
        if (idiag_ruyuzmx/=0) call yzsum_mn_name_x(p%rho*p%uu(:,2)*p%uu(:,3),idiag_ruyuzmx)
        if (idiag_ux2mz/=0) call xysum_mn_name_z(p%uu(:,1)**2,idiag_ux2mz)
        if (idiag_uy2mz/=0) call xysum_mn_name_z(p%uu(:,2)**2,idiag_uy2mz)
        if (idiag_uz2mz/=0) call xysum_mn_name_z(p%uu(:,3)**2,idiag_uz2mz)
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
!
!  2-D averages.
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        if (idiag_urmphi/=0) &
            call phisum_mn_name_rz(p%uu(:,1)*p%pomx+p%uu(:,2)*p%pomy,idiag_urmphi)
        if (idiag_ursphmphi/=0) &
            call phisum_mn_name_rz(p%uu(:,1)*p%evr(:,1)+ &
              p%uu(:,2)*p%evr(:,2)+p%uu(:,3)*p%evr(:,3),idiag_ursphmphi)
        if (idiag_uthmphi/=0) &
            call phisum_mn_name_rz(p%uu(:,1)*p%evth(:,1)+ &
              p%uu(:,2)*p%evth(:,2)+p%uu(:,3)*p%evth(:,3),idiag_uthmphi)
        if (idiag_upmphi/=0) &
            call phisum_mn_name_rz(p%uu(:,1)*p%phix+p%uu(:,2)*p%phiy,idiag_upmphi)
        if (idiag_uzmphi/=0) &
            call phisum_mn_name_rz(p%uu(:,3),idiag_uzmphi)
        if (idiag_u2mphi/=0) &
            call phisum_mn_name_rz(p%u2,idiag_u2mphi)
        if (idiag_ozmphi/=0) &
            call phisum_mn_name_rz(p%oo(:,3),idiag_ozmphi)
        if (idiag_oumphi/=0) call phisum_mn_name_rz(p%ou,idiag_oumphi)
!
        if (idiag_uxmxz/=0) &
            call ysum_mn_name_xz(p%uu(:,1),idiag_uxmxz)
        if (idiag_uymxz/=0) &
            call ysum_mn_name_xz(p%uu(:,2),idiag_uymxz)
        if (idiag_uzmxz/=0) &
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
        if (idiag_oumxz/=0) &
            call ysum_mn_name_xz(p%ou,idiag_oumxz)
!
        if (idiag_uxmxy/=0) call zsum_mn_name_xy(p%uu(:,1),idiag_uxmxy)
!
!  Changed calls for compatibility with Yin-Yang grid:
!  all non-scalars in which y or z components of a vector are used must
!  be treated as below,
!
        if (idiag_uymxy/=0) call zsum_mn_name_xy(p%uu,idiag_uymxy,(/0,1,0/))
        if (idiag_uzmxy/=0) call zsum_mn_name_xy(p%uu,idiag_uzmxy,(/0,0,1/))
        if (idiag_uxuymxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_uxuymxy,(/1,1,0/))
        if (idiag_uxuzmxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_uxuzmxy,(/1,0,1/))
        if (idiag_uyuzmxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_uyuzmxy,(/0,1,1/))
        if (idiag_oxmxy/=0) call zsum_mn_name_xy(p%oo(:,1),idiag_oxmxy)
        if (idiag_oymxy/=0) call zsum_mn_name_xy(p%oo,idiag_oymxy,(/0,1,0/))
        if (idiag_ozmxy/=0) call zsum_mn_name_xy(p%oo,idiag_ozmxy,(/0,0,1/))
        if (idiag_oumxy/=0) call zsum_mn_name_xy(p%ou,idiag_oumxy)
        if (idiag_ruxmxy/=0) call zsum_mn_name_xy(p%rho*p%uu(:,1),idiag_ruxmxy)
        if (idiag_ruymxy/=0) call zsum_mn_name_xy(p%uu,idiag_ruymxy,(/0,1,0/),p%rho)
        if (idiag_ruzmxy/=0) call zsum_mn_name_xy(p%uu,idiag_ruzmxy,(/0,0,1/),p%rho)
        if (idiag_ux2mxy/=0) &
            call zsum_mn_name_xy(p%uu(:,1)**2,idiag_ux2mxy)
        if (idiag_uy2mxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_uy2mxy,(/0,2,0/))
        if (idiag_uz2mxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_uz2mxy,(/0,0,2/))
        if (idiag_rux2mxy/=0) &
            call zsum_mn_name_xy(p%rho*p%uu(:,1)**2,idiag_rux2mxy)
        if (idiag_ruy2mxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_ruy2mxy,(/0,2,0/),p%rho)
        if (idiag_ruz2mxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_ruz2mxy,(/0,0,2/),p%rho)
        if (idiag_ruxuymxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_ruxuymxy,(/1,1,0/),p%rho)
        if (idiag_ruxuzmxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_ruxuzmxy,(/1,0,1/),p%rho)
        if (idiag_ruyuzmxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_ruyuzmxy,(/0,1,1/),p%rho)
        if (idiag_fkinxmxy/=0) &
            call zsum_mn_name_xy(p%ekin*p%uu(:,1),idiag_fkinxmxy)
        if (idiag_fkinymxy/=0) &
            call zsum_mn_name_xy(p%uu,idiag_fkinymxy,(/0,1,0/),p%ekin)
      else
!
!  idiag_uxmxy and idiag_uymxy also need to be calculated when
!  ldiagnos and idiag_umx and/or idiag_umy, so
!
!  We may need to calculate uxmxy without calculating umx. The following
!  if condition was messing up calculation of umxy_rms
!
        if (ldiagnos) then
          if (idiag_uxmxy/=0) call zsum_mn_name_xy(p%uu(:,1),idiag_uxmxy)
          if (idiag_uymxy/=0) call zsum_mn_name_xy(p%uu,idiag_uymxy,(/0,1,0/))
          if (idiag_uzmxy/=0) call zsum_mn_name_xy(p%uu,idiag_uzmxy,(/0,0,1/))
        endif
      endif
      call timing('duu_dt','finished',mnloop=.true.)
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine time_integrals_hydro(f,p)
!
!  Calculate time_integrals within each pencil (as long as each
!  pencil case p still contains the current data). This routine
!  is now being called at the end of equ.
!
!  28-jun-07/axel+mreinhard: coded
!  24-jun-08/axel: moved call to this routine to the individual pde routines
!   1-jul-08/axel: moved this part to hydro
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f
      intent(in) :: p
!
      if (iuut/=0) f(l1:l2,m,n,iuxt:iuzt)=f(l1:l2,m,n,iuxt:iuzt)+dt*p%uu
      if (ioot/=0) f(l1:l2,m,n,ioxt:iozt)=f(l1:l2,m,n,ioxt:iozt)+dt*p%oo
!
    endsubroutine time_integrals_hydro
!***********************************************************************
   subroutine coriolis_cartesian(df,uu,velind)
!
!  coriolis terms for cartesian geometry
!
!  30-oct-09/MR: outsourced, parameter velind added
!  checked to be an equivalent change by auto-test conv-slab-noequi, mdwarf
!
      real, dimension (mx,my,mz,mvar), intent(out) :: df
      real, dimension (nx,3),          intent(in)  :: uu
      integer,                         intent(in)  :: velind
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(uu)
      call keep_compiler_quiet(velind)
!
   endsubroutine coriolis_cartesian
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (3,3) :: mat_cent1=0.,mat_cent2=0.,mat_cent3=0.
      real, dimension (3) :: OO, dOO
      real :: c,s,sinalp,cosalp,OO2,alpha_precession_rad
      integer :: i,j
!
      intent(inout) :: f
!
!  possibility of setting interior boundary conditions
!
      if (lhydro_bc_interior) call interior_bc_hydro(f)
!
!    Slope limited diffusion: update characteristic speed
!    Not staggered yet
!
     if (lslope_limit_diff .and. llast) then
       do m=1,my
       do n=1,mz
           f(:,m,n,isld_char)=w_sldchar_hyd* &
            sqrt(f(:,m,n,iux)**2.+f(:,m,n,iuy)**2.+f(:,m,n,iuz)**2.)
       enddo
       enddo
     endif
!
    endsubroutine hydro_after_boundary
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
      othresh=othresh_scl*othresh_per_orms*orms
!
    endsubroutine calc_othresh
!***********************************************************************
    subroutine input_persistent_hydro(id,done)
!
      integer, intent(in), optional :: id
      logical, intent(inout), optional :: done
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
        idiag_fkinzm=0
        idiag_u2m=0
        idiag_um2=0
        idiag_uxpt=0
        idiag_uypt=0
        idiag_uzpt=0
        idiag_uxp2=0
        idiag_uyp2=0
        idiag_uzp2=0
        idiag_urms=0
        idiag_durms=0
        idiag_urmsx=0
        idiag_urmsz=0
        idiag_umax=0
        idiag_umin=0
        idiag_uxrms=0
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
        idiag_ux2ccm=0
        idiag_ux2ssm=0
        idiag_uy2ccm=0
        idiag_uy2ssm=0
        idiag_uxuycsm=0
        idiag_rux2m=0
        idiag_ruy2m=0
        idiag_ruz2m=0
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
        idiag_ruxmx=0
        idiag_ruymx=0
        idiag_ruzmx=0
        idiag_rux2mx = 0
        idiag_ruy2mx = 0
        idiag_ruz2mx = 0
        idiag_ruxuymx = 0
        idiag_ruxuzmx = 0
        idiag_ruyuzmx = 0
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
        idiag_accpowzmz=0
        idiag_accpowzupmz=0
        idiag_accpowzdownmz=0
        idiag_urmphi=0
        idiag_ursphmphi=0
        idiag_uthmphi=0
        idiag_upmphi=0
        idiag_uzmphi=0
        idiag_u2mphi=0
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
        idiag_uxuymxy=0
        idiag_uxuzmxy=0
        idiag_uyuzmxy=0
        idiag_oxmxy=0
        idiag_oymxy=0
        idiag_ozmxy=0
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
        idiag_ou_int=0
        idiag_fum=0
        idiag_odel2um=0
        idiag_o2m=0
        idiag_orms=0
        idiag_omax=0
        idiag_ox2m=0
        idiag_oy2m=0
        idiag_oz2m=0
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
        idiag_fkinxmx=0
        idiag_fkinxmxy=0
        idiag_fkinymxy=0
        idiag_ruxuym=0
        idiag_ruxuzm=0
        idiag_ruyuzm=0
        idiag_ruxuymz=0
        idiag_ruxuzmz=0
        idiag_ruyuzmz=0
        idiag_ruxuy2mz=0
        idiag_ruxuz2mz=0
        idiag_ruyuz2mz=0
        idiag_dudx=0
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
        idiag_nshift=0
        ivid_oo=0; ivid_o2=0; ivid_divu=0; ivid_u2=0; ivid_Ma2=0
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
        call parse_name(iname,cname(iname),cform(iname),'fkinzm',idiag_fkinzm)
        call parse_name(iname,cname(iname),cform(iname),'u2m',idiag_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',idiag_um2)
        call parse_name(iname,cname(iname),cform(iname),'odel2um',idiag_odel2um)
        call parse_name(iname,cname(iname),cform(iname),'o2m',idiag_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',idiag_oum)
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
        call parse_name(iname,cname(iname),cform(iname),'uzrms',idiag_uzrms)
        call parse_name(iname,cname(iname),cform(iname),'uzrmaxs',idiag_uzrmaxs)
        call parse_name(iname,cname(iname),cform(iname),'uxm',idiag_uxm)
        call parse_name(iname,cname(iname),cform(iname),'uym',idiag_uym)
        call parse_name(iname,cname(iname),cform(iname),'uzm',idiag_uzm)
        call parse_name(iname,cname(iname),cform(iname),'ux2m',idiag_ux2m)
        call parse_name(iname,cname(iname),cform(iname),'uy2m',idiag_uy2m)
        call parse_name(iname,cname(iname),cform(iname),'uz2m',idiag_uz2m)
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
        call parse_name(iname,cname(iname),cform(iname),'oxm',idiag_oxm)
        call parse_name(iname,cname(iname),cform(iname),'oym',idiag_oym)
        call parse_name(iname,cname(iname),cform(iname),'ozm',idiag_ozm)
        call parse_name(iname,cname(iname),cform(iname),'oxuzxm',idiag_oxuzxm)
        call parse_name(iname,cname(iname),cform(iname),'oyuzym',idiag_oyuzym)
        call parse_name(iname,cname(iname),cform(iname),'oxoym',idiag_oxoym)
        call parse_name(iname,cname(iname),cform(iname),'oxozm',idiag_oxozm)
        call parse_name(iname,cname(iname),cform(iname),'oyozm',idiag_oyozm)
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
        call parse_name(iname,cname(iname),cform(iname),'dudx',idiag_dudx)
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
        call parse_name(iname,cname(iname),cform(iname),'nshift',idiag_nshift)
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
            'ekinmz',idiag_ekinmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'u2mz',idiag_u2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'o2mz',idiag_o2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'curlru2mz',idiag_curlru2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'divru2mz',idiag_divru2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'divu2mz',idiag_divu2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'oumz',idiag_oumz)
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
             'accpowzmz',idiag_accpowzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'accpowzupmz',idiag_accpowzupmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
             'accpowzdownmz',idiag_accpowzdownmz)
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
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxuymxy',idiag_uxuymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uxuzmxy',idiag_uxuzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'uyuzmxy',idiag_uyuzmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'oxmxy',idiag_oxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'oymxy',idiag_oymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'ozmxy',idiag_ozmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'oumxy',idiag_oumxy)
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
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'fkinymxy',idiag_fkinymxy)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'urmphi',idiag_urmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ursphmphi',idiag_ursphmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uthmphi',idiag_uthmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'upmphi',idiag_upmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uzmphi',idiag_uzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'u2mphi',idiag_u2mphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'oumphi',idiag_oumphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ozmphi',idiag_ozmphi)
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
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'divu',ivid_divu)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'u2',  ivid_u2)
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'Ma2', ivid_Ma2)
      enddo
!
!  write column where which hydro variable is stored
!
!     if (lwr) then
!       if (lhelmholtz_decomp) call farray_index_append('iphiuu',iphiuu)
!     endif
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
          call assign_slices_scal(slices,divu_xy,divu_xz,divu_yz,divu_xy2,divu_xy3,divu_xy4,divu_xz2)
!
!  Velocity squared.
!
        case ('u2')
          call assign_slices_scal(slices,u2_xy,u2_xz,u2_yz,u2_xy2,u2_xy3,u2_xy4,u2_xz2)
!
!  Vorticity.
!
        case ('oo')
          call assign_slices_vec(slices,oo_xy,oo_xz,oo_yz,oo_xy2,oo_xy3,oo_xy4,oo_xz2)
!
!  Vorticity squared.
!
        case ('o2')
          call assign_slices_scal(slices,o2_xy,o2_xz,o2_yz,o2_xy2,o2_xy3,o2_xy4,o2_xz2)
!
!  Mach number squared.
!
        case ('mach')
          call assign_slices_scal(slices,mach_xy,mach_xz,mach_yz,mach_xy2,mach_xy3,mach_xy4,mach_xz2)
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
      integer :: ivar,ng,mg,ig,i
      real :: dt_
!
!  Pencil uses linear velocity. Fargo will shift based on
!  angular velocity. Get phidot from uphi.
!
      ifcoordinates: if (lcylindrical_coords) then
        zloopcyl: do n=n1,n2
          ng=n-n1+1
          phidot=uu_average_cyl(:,ng)*rcyl_mn1
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
          mg=m-m1+1
          xloopsph: do i=l1,l2
            ig=i-l1+1
            phidot=uu_average_sph(ig,mg)*rcyl_mn1
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
        !call remove_mean_flow(f,iux)         ! as this is equivalent to remove
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
    subroutine read_uumz_profile(uumz_prof)
!
!  Read radial profiles of horizontally averaged flows from file.
!
!  23-oct-18/pjk: added
!
      real, dimension(nz,3), intent(out) :: uumz_prof
      integer, parameter :: nztotal=nz*nprocz
      real, dimension(nz*nprocz) :: tmp1z,tmp2z,tmp3z
      real :: var1,var2,var3
      logical :: exist
      integer :: stat
!
!  Read hcond and glhc and write into an array.
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
          call fatal_error('read_uumz_profile','*** error *** - no input file')
        endif
      endif
!
!  Read profiles.
!
!  Gravity in the z-direction
!
      if (lgravz) then
        do n=1,nztotal
          read(31,*,iostat=stat) var1,var2,var3
          if (stat<0) exit
          if (ip<5) print*,'uxmz, uymz, uzmz: ',var1,var2,var3
          tmp1z(n)=var1
          tmp2z(n)=var2
          tmp3z(n)=var3
        enddo
!
!  Assuming no ghost zones in uumz.dat.
!
        do n=n1,n2
          uumz_prof(n-nghost,1)=tmp1z(ipz*nz+n-nghost)
          uumz_prof(n-nghost,2)=tmp2z(ipz*nz+n-nghost)
          uumz_prof(n-nghost,3)=tmp3z(ipz*nz+n-nghost)
        enddo
!
        close(31)
!
      endif
!
    endsubroutine read_uumz_profile
!***********************************************************************
    subroutine impose_velocity_ceiling(f)
!
!  13-aug-2007/anders: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_velocity_ceiling
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
    subroutine kinematic_random_phase
!
!  This is a dummy routine.
!
!  16-feb-2010/bing: coded
!
      call fatal_error('kinematic_random_phase', &
          'Use HYDRO=hydro_kinematic in Makefile.local instead')
!
    endsubroutine kinematic_random_phase
!***********************************************************************
    subroutine kinematic_random_ampl
!
!  This is a dummy routine.
!
!  26-jun-2019/axel: coded
!
      call fatal_error('kinematic_random_ampl', &
          'Use HYDRO=hydro_kinematic in Makefile.local instead')
!
    endsubroutine kinematic_random_ampl
!***********************************************************************
    subroutine kinematic_random_wavenumber
!
!  This is a dummy routine.
!
!  26-jun-2019/axel: coded
!
      call fatal_error('kinematic_random_wavenumber', &
          'Use HYDRO=hydro_kinematic in Makefile.local instead')
!
    endsubroutine kinematic_random_wavenumber
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
    subroutine calc_gradu(f)
!
    use Sub, only : gij
    real, dimension (mx,my,mz,mfarray) :: f
    integer :: imn,jk,jj,kk,nyz
    real, dimension(nx,3,3) :: gradu
!
! Calculated gradu and stores it as an auxiliary. This is expected to be called
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

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr_c(lpressuregradient_gas,p_par(1))  ! int

    endsubroutine pushpars2c
!***********************************************************************
    subroutine pushdiags2c(p_diag)

    use Diagnostics, only: set_type

    integer, parameter :: n_diags=12
    integer(KIND=ikind8), dimension(n_diags) :: p_diag

    call copy_addr_c(idiag_urms,p_diag(1))
    call set_type(idiag_urms,lsqrt=.true.)
    call copy_addr_c(idiag_uxrms,p_diag(2))
    call set_type(idiag_uxrms,lsqrt=.true.)
    call copy_addr_c(idiag_uyrms,p_diag(3))
    call set_type(idiag_uyrms,lsqrt=.true.)
    call copy_addr_c(idiag_uzrms,p_diag(4))
    call set_type(idiag_uzrms,lsqrt=.true.)
    call copy_addr_c(idiag_umax,p_diag(5))
    call set_type(idiag_umax,lmax=.true.)
    call copy_addr_c(idiag_umin,p_diag(6))
    call set_type(idiag_umin,lmin=.true.)
    call copy_addr_c(idiag_uxmin,p_diag(7))
    call set_type(idiag_uxmin,lmin=.true.)
    call copy_addr_c(idiag_uymin,p_diag(8))
    call set_type(idiag_uymin,lmin=.true.)
    call copy_addr_c(idiag_uzmin,p_diag(9))
    call set_type(idiag_uzmin,lmin=.true.)
    call copy_addr_c(idiag_uxmax,p_diag(10))
    call set_type(idiag_uxmax,lmax=.true.)
    call copy_addr_c(idiag_uymax,p_diag(11))
    call set_type(idiag_uymax,lmax=.true.)
    call copy_addr_c(idiag_uzmax,p_diag(12))
    call set_type(idiag_uzmax,lmax=.true.)

    endsubroutine pushdiags2c
!***********************************************************************
endmodule Hydro
