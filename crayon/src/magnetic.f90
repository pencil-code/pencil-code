! $Id: magnetic.f90 13608 2010-04-09 11:56:42Z mppiyali $
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
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED aa(3); a2; aij(3,3); bb(3); bbb(3); ab; ua; uxb(3); exa(3)
! PENCILS PROVIDED b2; bij(3,3); del2a(3); graddiva(3);jj(3)
! PENCILS PROVIDED j2; jb; va2; jxb(3); jxbr(3); jxbr2; ub; uxb(3); uxb2
! PENCILS PROVIDED uxj(3); beta; djuidjbi; jo
! PENCILS PROVIDED ujxb; oxu(3); oxuxb(3); jxbxb(3); jxbrxb(3)
! PENCILS PROVIDED glnrhoxb(3); oxj(3); diva
! PENCILS PROVIDED jij(3,3); sj; ss12
! PENCILS PROVIDED etava, etaj, etaj2, etajrho
! PENCILS PROVIDED cosjb,jparallel;jperp
! PENCILS PROVIDED cosub
!
!***************************************************************
module Magnetic
!
  use Cdata
  use Cparam
  use Messages, only: fatal_error,inevitably_fatal_error,warning,svn_id,timing
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'record_types.h'
  include 'magnetic.h'
!
! Slice precalculation buffers
!
  real, target, dimension (nx,ny,3) :: bb_xy, jj_xy
  real, target, dimension (nx,ny,3) :: bb_xy2, jj_xy2
  real, target, dimension (nx,nz,3) :: bb_xz, jj_xz
  real, target, dimension (ny,nz,3) :: bb_yz, jj_yz
!
  real, target, dimension (nx,ny) :: b2_xy, jb_xy
  real, target, dimension (nx,ny) :: b2_xy2, jb_xy2
  real, target, dimension (ny,nz) :: b2_yz, jb_yz
  real, target, dimension (nx,nz) :: b2_xz, jb_xz
!
  real, target, dimension (nx,ny) :: beta_xy
  real, target, dimension (nx,ny) :: beta_xy2
  real, target, dimension (ny,nz) :: beta_yz
  real, target, dimension (nx,nz) :: beta_xz
!
!  xy-averaged field
!
  real, dimension (mz,3) :: aamz
  real, dimension (nz,3) :: bbmz,jjmz
!
! Parameters
!
  integer, parameter :: nresi_max=4
!
  real, dimension (ninit) :: amplaa=0.0, kx_aa=1.0, ky_aa=1.0, kz_aa=1.0
  real, dimension (ninit) :: phasex_aa=0.0, phasey_aa=0.0, phasez_aa=0.0
  character (len=labellen), dimension(ninit) :: initaa='nothing'
  character (len=labellen), dimension(nresi_max) :: iresistivity=''
  character (len=labellen) :: fring_profile='tanh'
!
! Input parameters
!
  complex, dimension(3) :: coefaa=(/0.0,0.0,0.0/), coefbb=(/0.0,0.0,0.0/)
  real, dimension(3) :: B_ext=(/0.0,0.0,0.0/), B1_ext, B_ext_inv, B_ext_tmp
  real, dimension(3) :: axisr1=(/0,0,1/), dispr1=(/0.0,0.5,0.0/)
  real, dimension(3) :: axisr2=(/1,0,0/), dispr2=(/0.0,-0.5,0.0/)
  real, dimension(3) :: axisr3=(/1,0,0/), dispr3=(/0.0,-0.5,0.0/)
  real, dimension(nx,3) :: uxbb
  real, target :: zmode=1.0 !(temporary)
  real :: fring1=0.0, Iring1=0.0, Rring1=1.0, wr1=0.3
  real :: fring2=0.0, Iring2=0.0, Rring2=1.0, wr2=0.3
  real :: fring3=0.0, Iring3=0.0, Rring3=1.0, wr3=0.3
  real :: radius=0.1, epsilonaa=0.01, widthaa=0.5, x0aa=0.0, z0aa=0.0
  real :: by_left=0.0, by_right=0.0, bz_left=0.0, bz_right=0.0
  real :: bthresh=0.0, bthresh_per_brms=0.0, brms=0.0, bthresh_scl=1.0
  real :: eta_shock=0.0
  real :: eta_va=0., eta_j=0., eta_j2=0., eta_jrho=0., eta_min=0., etaj20=0.
  real :: rhomin_jxb=0.0, va2max_jxb=0.0
  real :: omega_Bz_ext=0.0
  real :: mu_r=-0.5 !(still needed for backwards compatibility)
  real :: mu_ext_pot=-0.5,inclaa=0.0
  real :: mu012=0.5 !(=1/2mu0)
  real :: rescale_aa=0.0
  real :: ampl_B0=0.0, D_smag=0.17, B_ext21, B_ext11
  real :: Omega_ampl=0.0
  real :: rmode=1.0, rm_int=0.0, rm_ext=0.0
  real :: nu_ni=0.0, nu_ni1
  real :: displacement_gun=0.0
  real :: initpower_aa=0.0, cutoff_aa=0.0, brms_target=1.0
  real :: rescaling_fraction=1.0
  real :: phase_beltrami=0.0, ampl_beltrami=0.0
  real :: bmz=0, bmz_beltrami_phase=0.0
  real :: taareset=0.0, daareset=0.0
  real :: center1_x=0.0, center1_y=0.0, center1_z=0.0
  real :: fluxtube_border_width=impossible
  real :: eta_jump=0.0
  integer :: nbvec,nbvecmax=nx*ny*nz/4, va2power_jxb=5
  integer :: N_modes_aa=1, naareset
  integer :: nrings=2
  logical :: lpress_equil=.false., lpress_equil_via_ss=.false.
  logical :: lpress_equil_alt=.false.
  logical :: llorentzforce=.true., linduction=.true.
  logical :: lresi_eta_const=.false.
  logical :: lresi_etaSS=.false.
  logical :: lresi_eta_shock=.false.
  logical :: lresi_etava=.false.
  logical :: lresi_etaj=.false.
  logical :: lresi_etaj2=.false.
  logical :: lresi_etajrho=.false.
  logical :: lresi_anomalous=.false.
  logical, target, dimension (3) :: lfrozen_bb_bot=(/.false.,.false.,.false./)
  logical, target, dimension (3) :: lfrozen_bb_top=(/.false.,.false.,.false./)
  logical :: lohmic_heat=.true., lneutralion_heat=.true.
  logical :: reinitialize_aa=.false.
  logical :: lB_ext_pot=.false.
  logical :: lforce_free_test=.false.
  logical :: lgauss=.false.
  logical :: lbb_as_aux=.false., ljj_as_aux=.false.
  logical :: lbbt_as_aux=.false., ljjt_as_aux=.false.
  logical :: lcheck_positive_va2=.false.
  logical :: lreset_aa=.false.
!
  namelist /magnetic_init_pars/ &
      B_ext, lohmic_heat, fring1, Iring1, Rring1, wr1, axisr1, dispr1, fring2, &
      Iring2, Rring2, wr2, axisr2, dispr2, fring3, Iring3, Rring3, wr3,  &
      axisr3, dispr3, fring_profile, nrings, radius, epsilonaa, x0aa, z0aa, &
      widthaa, by_left, by_right, bz_left, bz_right, initaa, amplaa, kx_aa, &
      ky_aa, kz_aa, coefaa, coefbb, phasex_aa, phasey_aa, phasez_aa, inclaa, &
      lpress_equil, lpress_equil_via_ss, mu_r, mu_ext_pot, lB_ext_pot, &
      lforce_free_test, ampl_B0, initpower_aa, cutoff_aa, N_modes_aa, rmode, &
      zmode, rm_int, rm_ext, lgauss, lcheck_positive_va2, lbb_as_aux, &
      ljj_as_aux, lbbt_as_aux, ljjt_as_aux, &
      lneutralion_heat, center1_x, center1_y, center1_z, &
      fluxtube_border_width, va2max_jxb, va2power_jxb, eta_jump,&
      lpress_equil_alt
!
! Run parameters
!
  real :: eta=0.0, eta1=0.0, eta_hyper2=0.0, eta_hyper3=0.0, eta_anom=0.0
  real :: eta_int=0.0, eta_ext=0.0, wresistivity=0.01, eta_xy_max=1.0
  real :: height_eta=0.0, eta_out=0.0
  real :: z_surface=0.0
  real :: tau_aa_exterior=0.0
  real :: sigma_ratio=1.0, eta_width=0.0, eta_z0=1.0, eta_z1=1.0
  real :: alphaSSm=0.0
  real :: alpha_rmax=0.0, alpha_width=0.0
  real :: Omega_rmax=0.0, Omega_rwidth=0.0
  real :: k1_ff=1.0, ampl_ff=1.0, swirl=1.0
  real :: k1x_ff=1.0, k1y_ff=1.0, k1z_ff=1.0
  real :: inertial_length=0.0, linertial_2
  real :: vcrit_anom=1.0
  real, dimension(mx,my) :: eta_xy
  real, dimension(mx,my,3) :: geta_xy
  real, dimension(mz) :: coskz,sinkz,eta_z
  real, dimension(mz,3) :: geta_z
  logical :: lcalc_aamean=.false.
  logical :: luse_Bext_in_b2=.false.
  logical :: lmean_friction=.false.
  logical :: llarge_scale_velocity=.false.
  logical :: lhalox=.false.
  logical :: lrun_initaa=.false.
  character (len=labellen) :: zdep_profile='fs'
  character (len=labellen) :: eta_xy_profile='schnack89'
!
  namelist /magnetic_run_pars/ &
      eta, eta1, eta_anom, B_ext, omega_Bz_ext, nu_ni, &
      lohmic_heat,eta_out,z_surface, &
      tau_aa_exterior, kx_aa, ky_aa, kz_aa, &
      lcalc_aamean, k1_ff, &
      ampl_ff, swirl, radius, k1x_ff, k1y_ff, k1z_ff, lcheck_positive_va2, &
      lmean_friction, bthresh, bthresh_per_brms, iresistivity, &
      alphaSSm, alpha_rmax, alpha_width, eta_int, &
      eta_ext, eta_shock, eta_va,eta_j, eta_j2, eta_jrho, eta_min, &
      wresistivity, eta_xy_max, rhomin_jxb, va2max_jxb, &
      va2power_jxb, llorentzforce, linduction, reinitialize_aa, rescale_aa, &
      lB_ext_pot, displacement_gun, D_smag, brms_target, &
      rescaling_fraction, sigma_ratio, zdep_profile, eta_width, &
      eta_z0, eta_z1, inertial_length, &
      lbb_as_aux, ljj_as_aux, &
      lbbt_as_aux, ljjt_as_aux, lneutralion_heat, lreset_aa, daareset, &
      luse_Bext_in_b2, llarge_scale_velocity, &
      lhalox, vcrit_anom, eta_jump,&
      Omega_rmax,Omega_rwidth,lrun_initaa
!
! Diagnostic variables (need to be consistent with reset list below)
!
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
  integer :: idiag_ajm=0        ! DIAG_DOC: $\left<\jv\cdot\Av\right>$
  integer :: idiag_jbm=0        ! DIAG_DOC: $\left<\jv\cdot\Bv\right>$
  integer :: idiag_jbmh=0       ! DIAG_DOC: $\left<\Av\cdot\Bv\right>$ (temp)
  integer :: idiag_jbmn=0       ! DIAG_DOC: $\left<\Av\cdot\Bv\right>$ (north)
  integer :: idiag_jbms=0       ! DIAG_DOC: $\left<\Av\cdot\Bv\right>$ (south)
  integer :: idiag_ubm=0        ! DIAG_DOC: $\left<\uv\cdot\Bv\right>$
  integer :: idiag_uxbxm=0      ! DIAG_DOC: $\left<u_xB_x\right>$
  integer :: idiag_uybym=0      ! DIAG_DOC: $\left<u_yB_y\right>$
  integer :: idiag_uzbzm=0      ! DIAG_DOC: $\left<u_zB_z\right>$
  integer :: idiag_cosubm=0     ! DIAG_DOC: $\left<\Uv\cdot\Bv/(|\Uv|\,|\Bv|)\right>$
  integer :: idiag_uam=0        ! DIAG_DOC: $\left<\uv\cdot\Av\right>$
  integer :: idiag_ujm=0        ! DIAG_DOC: $\left<\uv\cdot\Jv\right>$
  integer :: idiag_fbm=0        ! DIAG_DOC: $\left<\fv\cdot\Bv\right>$
  integer :: idiag_fxbxm=0      ! DIAG_DOC: $\left<f_x B_x\right>$
  integer :: idiag_epsM=0       ! DIAG_DOC: $\left<2\eta\mu_0\jv^2\right>$
  integer :: idiag_epsAD=0      ! DIAG_DOC: $\left<\rho^{-1} t_{\rm AD} (\vec{J}\times\vec{B})^2\right>$ (heating by ion-neutrals friction)
  integer :: idiag_bxpt=0       ! DIAG_DOC: $B_x(x_0,y_0,z_0,t)$
  integer :: idiag_bypt=0       ! DIAG_DOC: $B_y(x_0,y_0,z_0,t)$
  integer :: idiag_bzpt=0       ! DIAG_DOC: $B_z(x_0,y_0,z_0,t)$
  integer :: idiag_Expt=0       ! DIAG_DOC: ${\cal E}_x(x_0,y_0,z_0,t)$
  integer :: idiag_Eypt=0       ! DIAG_DOC: ${\cal E}_y(x_0,y_0,z_0,t)$
  integer :: idiag_Ezpt=0       ! DIAG_DOC: ${\cal E}_z(x_0,y_0,z_0,t)$
  integer :: idiag_epsM_LES=0   ! DIAG_DOC:
  integer :: idiag_aybym2=0     ! DIAG_DOC:
  integer :: idiag_exaym2=0     ! DIAG_DOC:
  integer :: idiag_exabot=0     ! DIAG_DOC: $\int\Ev\times\Av\,dS|_{\rm bot}$
  integer :: idiag_exatop=0     ! DIAG_DOC: $\int\Ev\times\Av\,dS|_{\rm top}$
  integer :: idiag_exjm2=0      ! DIAG_DOC:
  integer :: idiag_emag=0       ! DIAG_DOC: $\int_V{1\over2\mu_0}\Bv^2\, dV$
  integer :: idiag_brms=0       ! DIAG_DOC: $\left<\Bv^2\right>^{1/2}$
  integer :: idiag_bmax=0       ! DIAG_DOC: $\max(|\Bv|)$
  integer :: idiag_bxmin=0      ! DIAG_DOC: $\min(|B_x|)$
  integer :: idiag_bymin=0      ! DIAG_DOC: $\min(|B_y|)$
  integer :: idiag_bzmin=0      ! DIAG_DOC: $\min(|B_z|)$
  integer :: idiag_bxmax=0      ! DIAG_DOC: $\max(|B_x|)$
  integer :: idiag_bymax=0      ! DIAG_DOC: $\max(|B_y|)$
  integer :: idiag_bzmax=0      ! DIAG_DOC: $\max(|B_z|)$
  integer :: idiag_jrms=0       ! DIAG_DOC: $\left<\jv^2\right>^{1/2}$
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
  integer :: idiag_axm=0        ! DIAG_DOC:
  integer :: idiag_aym=0        ! DIAG_DOC:
  integer :: idiag_azm=0        ! DIAG_DOC:
  integer :: idiag_arms=0       ! DIAG_DOC:
  integer :: idiag_amax=0       ! DIAG_DOC:
  integer :: idiag_Qsm=0        ! DIAG_DOC: $\left<Q_p(\overline{B})\right>$
  integer :: idiag_Qpm=0        ! DIAG_DOC: $\left<Q_p(\overline{B})\right>$
  integer :: idiag_qem=0        ! DIAG_DOC: $\left<q_e(\overline{B})\right>$
  integer :: idiag_beta1m=0     ! DIAG_DOC: $\left<\Bv^2/(2\mu_0 p)\right>$
                                ! DIAG_DOC:   \quad(mean inverse plasma beta)
  integer :: idiag_beta1max=0   ! DIAG_DOC: $\max[\Bv^2/(2\mu_0 p)]$
                                ! DIAG_DOC:   \quad(maximum inverse plasma beta)
  integer :: idiag_bxm=0        ! DIAG_DOC: $\left<\left<B\right>_{yz}^2\right>_x^{1/2}$
  integer :: idiag_bym=0        ! DIAG_DOC: $\left<\left<B\right>_{zx}^2\right>_y^{1/2}$
  integer :: idiag_bzm=0        ! DIAG_DOC: $\left<\left<B\right>_{xy}^2\right>_z^{1/2}$
  integer :: idiag_bx2m=0       ! DIAG_DOC: $\left<B_x^2\right>$
  integer :: idiag_by2m=0       ! DIAG_DOC: $\left<B_y^2\right>$
  integer :: idiag_bz2m=0       ! DIAG_DOC: $\left<B_z^2\right>$
  integer :: idiag_bx2mx=0      ! DIAG_DOC: $\left<B_x^2\right>_{yz}$
  integer :: idiag_by2mx=0      ! DIAG_DOC: $\left<B_y^2\right>_{yz}$
  integer :: idiag_bz2mx=0      ! DIAG_DOC: $\left<B_z^2\right>_{yz}$
  integer :: idiag_bxbym=0      ! DIAG_DOC: $\left<B_x B_y\right>$
  integer :: idiag_bxbymx=0     ! DIAG_DOC: $\left<B_x B_y\right>_{yz}$
  integer :: idiag_bxbzm=0      ! DIAG_DOC:
  integer :: idiag_bybzm=0      ! DIAG_DOC:
  integer :: idiag_djuidjbim=0  ! DIAG_DOC:
  integer :: idiag_bxbymy=0     ! DIAG_DOC:
  integer :: idiag_bxbzmy=0     ! DIAG_DOC:
  integer :: idiag_bybzmy=0     ! DIAG_DOC:
  integer :: idiag_bxbymz=0     ! DIAG_DOC:
  integer :: idiag_bxbzmz=0     ! DIAG_DOC:
  integer :: idiag_bybzmz=0     ! DIAG_DOC:
  integer :: idiag_b2mz=0       ! DIAG_DOC:
  integer :: idiag_axmz=0       ! DIAG_DOC: $\left<{\cal A}_x\right>_{xy}$
  integer :: idiag_aymz=0       ! DIAG_DOC: $\left<{\cal A}_y\right>_{xy}$
  integer :: idiag_azmz=0       ! DIAG_DOC: $\left<{\cal A}_z\right>_{xy}$
  integer :: idiag_abuxmz=0     ! DIAG_DOC: $\left<(\Av \cdot \Bv) u_x \right>_{xy}$
  integer :: idiag_abuymz=0     ! DIAG_DOC: $\left<(\Av \cdot \Bv) u_y \right>_{xy}$
  integer :: idiag_abuzmz=0     ! DIAG_DOC: $\left<(\Av \cdot \Bv) u_z \right>_{xy}$
  integer :: idiag_uabxmz=0     ! DIAG_DOC: $\left<(\uv \cdot \Av) B_x \right>_{xy}$
  integer :: idiag_uabymz=0     ! DIAG_DOC: $\left<(\uv \cdot \Av) B_y \right>_{xy}$
  integer :: idiag_uabzmz=0     ! DIAG_DOC: $\left<(\uv \cdot \Av) B_z \right>_{xy}$
  integer :: idiag_bxmz=0       ! DIAG_DOC: $\left<{\cal B}_x\right>_{xy}$
  integer :: idiag_bymz=0       ! DIAG_DOC: $\left<{\cal B}_y\right>_{xy}$
  integer :: idiag_bzmz=0       ! DIAG_DOC: $\left<{\cal B}_z\right>_{xy}$
  integer :: idiag_jxmz=0       ! DIAG_DOC: $\left<{\cal J}_x\right>_{xy}$
  integer :: idiag_jymz=0       ! DIAG_DOC: $\left<{\cal J}_y\right>_{xy}$
  integer :: idiag_jzmz=0       ! DIAG_DOC: $\left<{\cal J}_z\right>_{xy}$
  integer :: idiag_Exmz=0       ! DIAG_DOC: $\left<{\cal E}_x\right>_{xy}$
  integer :: idiag_Eymz=0       ! DIAG_DOC: $\left<{\cal E}_y\right>_{xy}$
  integer :: idiag_Ezmz=0       ! DIAG_DOC: $\left<{\cal E}_z\right>_{xy}$
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
  integer :: idiag_kx_aa=0      ! DIAG_DOC: $k_x$
  integer :: idiag_kmz=0        ! DIAG_DOC: $\left<\left<\Jv\cdot\Bv\right>_{xy}\right>/
                                ! DIAG_DOC:  \left<\left<\Jv\cdot\Bv\right>_{xy}\right>$
  integer :: idiag_bx2my=0      ! DIAG_DOC: $\left< B_x^2 \right>_{xz}$
  integer :: idiag_by2my=0      ! DIAG_DOC: $\left< B_y^2 \right>_{xz}$
  integer :: idiag_bz2my=0      ! DIAG_DOC: $\left< B_z^2 \right>_{xz}$
  integer :: idiag_bx2mz=0      ! DIAG_DOC: $\left< B_x^2 \right>_{xy}$
  integer :: idiag_by2mz=0      ! DIAG_DOC: $\left< B_y^2 \right>_{xy}$
  integer :: idiag_bz2mz=0      ! DIAG_DOC: $\left< B_z^2 \right>_{xy}$
  integer :: idiag_bx2rmz=0     ! DIAG_DOC: $\left< B_x^2/\varrho \right>_{xy}$
  integer :: idiag_by2rmz=0     ! DIAG_DOC: $\left< B_y^2/\varrho \right>_{xy}$
  integer :: idiag_bz2rmz=0     ! DIAG_DOC: $\left< B_z^2/\varrho \right>_{xy}$
  integer :: idiag_jbmz=0       ! DIAG_DOC: $\left<\Jv\cdot\Bv\right>|_{xy}$
  integer :: idiag_abmz=0       ! DIAG_DOC: $\left<\Av\cdot\Bv\right>|_{xy}$
  integer :: idiag_uamz=0       ! DIAG_DOC: $\left<\uv\cdot\Av\right>|_{xy}$
  integer :: idiag_examz1=0     ! DIAG_DOC: $\left<\Ev\times\Av\right>_{xy}|_x$
  integer :: idiag_examz2=0     ! DIAG_DOC: $\left<\Ev\times\Av\right>_{xy}|_y$
  integer :: idiag_examz3=0     ! DIAG_DOC: $\left<\Ev\times\Av\right>_{xy}|_z$
  integer :: idiag_bxmxy=0      ! DIAG_DOC: $\left< B_x \right>_{xy}$
  integer :: idiag_bymxy=0      ! DIAG_DOC: $\left< B_y \right>_{xy}$
  integer :: idiag_bzmxy=0      ! DIAG_DOC: $\left< B_z \right>_{xy}$
  integer :: idiag_jxmxy=0      ! DIAG_DOC: $\left< J_x \right>_{xy}$
  integer :: idiag_jymxy=0      ! DIAG_DOC: $\left< J_y \right>_{xy}$
  integer :: idiag_jzmxy=0      ! DIAG_DOC: $\left< J_z \right>_{xy}$
  integer :: idiag_axmxz=0      ! DIAG_DOC: $\left< A_x \right>_{xz}$
  integer :: idiag_aymxz=0      ! DIAG_DOC: $\left< A_y \right>_{xz}$
  integer :: idiag_azmxz=0      ! DIAG_DOC: $\left< A_z \right>_{xz}$
  integer :: idiag_bxmxz=0      ! DIAG_DOC: $\left< B_x \right>_{xz}$
  integer :: idiag_bymxz=0      ! DIAG_DOC: $\left< B_y \right>_{xz}$
  integer :: idiag_bzmxz=0      ! DIAG_DOC: $\left< B_z \right>_{xz}$
  integer :: idiag_bx2mxz=0     ! DIAG_DOC: $\left< B_x^2 \right>_{xz}$
  integer :: idiag_by2mxz=0     ! DIAG_DOC: $\left< B_y^2 \right>_{xz}$
  integer :: idiag_bz2mxz=0     ! DIAG_DOC: $\left< B_z^2 \right>_{xz}$
  integer :: idiag_bx2mxy=0     ! DIAG_DOC: $\left< B_x^2 \right>_{z}$
  integer :: idiag_by2mxy=0     ! DIAG_DOC: $\left< B_y^2 \right>_{z}$
  integer :: idiag_bz2mxy=0     ! DIAG_DOC: $\left< B_z^2 \right>_{z}$
  integer :: idiag_bxbymxy=0    ! DIAG_DOC: $\left< B_x B_y \right>_{z}$
  integer :: idiag_bxbzmxy=0    ! DIAG_DOC: $\left< B_x B_z \right>_{z}$
  integer :: idiag_bybzmxy=0    ! DIAG_DOC: $\left< B_y B_z \right>_{z}$
  integer :: idiag_bxbymxz=0    ! DIAG_DOC: $\left< B_x B_y \right>_{y}$
  integer :: idiag_bxbzmxz=0    ! DIAG_DOC: $\left< B_x B_z \right>_{y}$
  integer :: idiag_bybzmxz=0    ! DIAG_DOC: $\left< B_y B_z \right>_{y}$
  integer :: idiag_uxbm=0       ! DIAG_DOC: $\left<\uv\times\Bv\right>\cdot\Bv_0/B_0^2$
  integer :: idiag_jxbm=0       ! DIAG_DOC: $\left<\jv\times\Bv\right>\cdot\Bv_0/B_0^2$
  integer :: idiag_oxuxbm=0     ! DIAG_DOC:
  integer :: idiag_jxbxbm=0     ! DIAG_DOC:
  integer :: idiag_gpxbm=0      ! DIAG_DOC:
  integer :: idiag_uxDxuxbm=0   ! DIAG_DOC:
  integer :: idiag_b3b21m=0     ! DIAG_DOC: $\left<B_3 B_{2,1} \right>$
  integer :: idiag_b1b32m=0     ! DIAG_DOC: $\left<B_1 B_{3,2} \right>$
  integer :: idiag_b2b13m=0     ! DIAG_DOC: $\left<B_2 B_{1,3} \right>$
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
  integer :: idiag_ujxbm=0      ! DIAG_DOC: $\left<\uv\cdot(\Jv\times\Bv\right>$
  integer :: idiag_jxbrxm=0     ! DIAG_DOC:
  integer :: idiag_jxbrym=0     ! DIAG_DOC:
  integer :: idiag_jxbrzm=0     ! DIAG_DOC:
  integer :: idiag_jxbr2m=0     ! DIAG_DOC: $\left<(\Jv\times\Bv/\rho)^2\right>$
  integer :: idiag_jxbrxmx=0    ! DIAG_DOC:
  integer :: idiag_jxbrymx=0    ! DIAG_DOC:
  integer :: idiag_jxbrzmx=0    ! DIAG_DOC:
  integer :: idiag_jxbrxmy=0    ! DIAG_DOC:
  integer :: idiag_jxbrymy=0    ! DIAG_DOC:
  integer :: idiag_jxbrzmy=0    ! DIAG_DOC:
  integer :: idiag_jxbrxmz=0    ! DIAG_DOC:
  integer :: idiag_jxbrymz=0    ! DIAG_DOC:
  integer :: idiag_jxbrzmz=0    ! DIAG_DOC:
  integer :: idiag_uxBrms=0     ! DIAG_DOC:
  integer :: idiag_Bresrms=0    ! DIAG_DOC:
  integer :: idiag_Rmrms=0      ! DIAG_DOC:
  integer :: idiag_jfm=0        ! DIAG_DOC:
  integer :: idiag_vA2m=0       ! DIAG_DOC:
  integer :: idiag_bxmx=0       ! DIAG_DOC:
  integer :: idiag_bymx=0       ! DIAG_DOC:
  integer :: idiag_bzmx=0       ! DIAG_DOC:
  integer :: idiag_bxmy=0       ! DIAG_DOC:
  integer :: idiag_bymy=0       ! DIAG_DOC:
  integer :: idiag_bzmy=0       ! DIAG_DOC:
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
  integer :: idiag_cosjbm=0     ! DIAG_DOC: $\left<\Jv\cdot\Bv/(|\Jv|\,|\Bv|)\right>$
  integer :: idiag_jparallelm=0 ! DIAG_DOC: Mean value of the component
                                ! DIAG_DOC: of J parallel to B
  integer :: idiag_jperpm=0     ! DIAG_DOC: Mean value of the component
                                ! DIAG_DOC: of J perpendicular to B
  integer :: idiag_brmsn=0,idiag_brmss=0,idiag_brmsh=0
  integer :: idiag_Exmxz=0      ! DIAG_DOC: $\left<{\cal E}_x\right>_{y}$
  integer :: idiag_Eymxz=0      ! DIAG_DOC: $\left<{\cal E}_y\right>_{y}$
  integer :: idiag_Ezmxz=0      ! DIAG_DOC: $\left<{\cal E}_z\right>_{y}$
  integer :: idiag_etatotalmx=0 ! DIAG_DOC: $\left<\eta\right>_{yz}$
  integer :: idiag_etatotalmz=0 ! DIAG_DOC: $\left<\eta\right>_{xy}$
!
  contains
!***********************************************************************
    subroutine register_magnetic()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!
!  1-may-02/wolf: coded
!
      use FArrayManager, only: farray_register_pde
!
      call farray_register_pde('aa',iaa,vector=3)
      iax = iaa; iay = iaa+1; iaz = iaa+2
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id: magnetic.f90 13608 2010-04-09 11:56:42Z mppiyali $")
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
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic(f,lstarting)
!
!  Perform any post-parameter-read initialization
!
!  24-nov-02/tony: dummy routine - nothing to do at present
!  20-may-03/axel: reinitialize_aa added
!
      use FArrayManager
      use SharedVariables
      use EquationOfState, only: cs0
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      integer :: i
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
        nu_ni1=1./nu_ni
      else
        nu_ni1=0.0
      endif
!
!  calculate B_ext21 = 1/B_ext**2 and the unit vector B1_ext = B_ext/|B_ext|
!  Also calculate B_ext_inv = B_ext/|B_ext|^2
!
      B_ext21=B_ext(1)**2+B_ext(2)**2+B_ext(3)**2
      if (B_ext21/=0.0) then
        B_ext21=1/B_ext21
      else
        B_ext21=1.0
      endif
      B_ext11=sqrt(B_ext21)
      B1_ext=B_ext*B_ext11
      B_ext_inv=B_ext*B_ext21
!
!  Initialize resistivity.
!
      if (iresistivity(1)=='') iresistivity(1)='eta-const'  ! default
      lresi_eta_const=.false.
      lresi_eta_shock=.false.
!
      do i=1,nresi_max
        select case (iresistivity(i))
        case ('eta-const')
          if (lroot) print*, 'resistivity: constant eta'
          lresi_eta_const=.true.
        case ('shock','eta-shock')
          if (lroot) print*, 'resistivity: shock'
          lresi_eta_shock=.true.
          if (.not. lshock) &
              call fatal_error('initialize_magnetic', &
              'shock resistivity, but module setting SHOCK=noshock')
        case ('none')
          ! do nothing
        case ('')
          ! do nothing
        case default
          if (lroot) print*, 'No such such value for iresistivity(',i,'): ', &
              trim(iresistivity(i))
          call fatal_error('initialize_magnetic','')
        endselect
      enddo
!
!  If we're timestepping, die or warn if the the resistivity coefficient that
!  corresponds to the chosen resistivity type is not set.
!
      if (lrun) then
        if (lresi_eta_const.and.(eta==0.0)) &
            call warning('initialize_magnetic', &
            'Resistivity coefficient eta is zero!')
        if (lresi_eta_shock.and.eta_shock==0.0) &
            call fatal_error('initialize_magnetic', &
            'Resistivity coefficient eta_shock is zero!')
      endif
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
      if (lbb_as_aux) then
        if (ibb==0) then
          call farray_register_auxiliary('bb',ibb,vector=3)
          ibx=ibb
          iby=ibb+1
          ibz=ibb+2
        endif
        if (ibb/=0.and.lroot) then
          print*, 'initialize_magnetic: ibb = ', ibb
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'ibb=',ibb
          write(3,*) 'ibx=',ibx
          write(3,*) 'iby=',iby
          write(3,*) 'ibz=',ibz
          close(3)
        endif
      endif
!
!  do the same for jj (current density)
!
      if (ljj_as_aux) then
        if (ijj==0) then
          call farray_register_auxiliary('jj',ijj,vector=3)
          ijx=ijj
          ijy=ijj+1
          ijz=ijj+2
        endif
        if (ijj/=0.and.lroot) then
          print*, 'initialize_magnetic: ijj = ', ijj
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'ijj=',ijj
          write(3,*) 'ijx=',ijx
          write(3,*) 'ijy=',ijy
          write(3,*) 'ijz=',ijz
          close(3)
        endif
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_magnetic
!***********************************************************************
    subroutine init_aa(f)
!
!  initialise magnetic field; called from start.f90
!  AB: maybe we should here call different routines (such as rings)
!  AB: and others, instead of accummulating all this in a huge routine.
!  We have an init parameter (initaa) to stear magnetic i.c. independently.
!
!   7-nov-2001/wolf: coded
!
      use EquationOfState
      use FArrayManager
      use Gravity, only: gravz, z1, z2
      use Initcond
      use Boundcond
      use InitialCondition, only: initial_condition_aa
      use Mpicomm
      use SharedVariables
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (mz) :: tmp
      integer :: j
!
      do j=1,ninit
!
        select case (initaa(j))
        case ('nothing'); if (lroot .and. j==1) print*,'init_aa: nothing'
        case ('zero', '0'); f(:,:,:,iax:iaz) = 0.0
        case ('rescale'); f(:,:,:,iax:iaz)=amplaa(j)*f(:,:,:,iax:iaz)
        case ('bsiny'); call acosy(amplaa(j),f,iaa,ky_aa(j))
        case ('mode'); call modev(amplaa(j),coefaa,f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('modeb'); call modeb(amplaa(j),coefbb,f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('sph_constb'); call sph_constb(amplaa(j),f,iaa)
        case ('random-isotropic-KS')
          call random_isotropic_KS(initpower_aa,f,iax,N_modes_aa)
        case ('gaussian-noise'); call gaunoise(amplaa(j),f,iax,iaz)
        case ('gaussian-noise-zprof')
          tmp=amplaa(1)*0.5*(tanh((z-z1)/0.05)-tanh((z-z2)/0.05))
          call gaunoise(tmp,f,iax,iaz)
!
!  Beltrami fields, put k=-k to make sure B=curl(A) has the right phase
!
        case ('Beltrami-x')
               call beltrami(amplaa(j),f,iaa,KX=-kx_aa(j),phase=phasex_aa(j))
        case ('Beltrami-xy')
               call beltrami(amplaa(j),f,iaa,KX=kx_aa(j),phase=phasex_aa(j))
               call beltrami(amplaa(j),f,iaa,KY=-kx_aa(j),phase=phasex_aa(j))
        case ('Beltrami-y'); call beltrami(amplaa(j),f,iaa,KY=-ky_aa(j),phase=phasey_aa(j))
        case ('Beltrami-z'); call beltrami(amplaa(j),f,iaa,KZ=-kz_aa(j),phase=phasez_aa(j))
!
        case ('propto-ux'); call wave_uu(amplaa(j),f,iaa,kx=kx_aa(j))
        case ('propto-uy'); call wave_uu(amplaa(j),f,iaa,ky=ky_aa(j))
        case ('propto-uz'); call wave_uu(amplaa(j),f,iaa,kz=kz_aa(j))
        case ('diffrot'); call diffrot(amplaa(j),f,iay)
        case ('hor-tube'); call htube(amplaa(j),f,iax,iaz,radius,epsilonaa, &
                                     center1_x,center1_z)
        case ('hor-tube-x'); call htube_x(amplaa(j),f,iax,iaz,radius,epsilonaa, &
                                     center1_y,center1_z)
        case ('hor-tube_erf'); call htube_erf(amplaa(j),f,iax,iaz,radius,epsilonaa, &
                                     center1_x,center1_z,fluxtube_border_width)
        case ('hor-fluxlayer'); call hfluxlayer(amplaa(j),f,iaa,z0aa,widthaa)
        case ('hor-fluxlayer-y'); call hfluxlayer_y(amplaa(j),f,iaa,z0aa,widthaa)
        case ('ver-fluxlayer'); call vfluxlayer(amplaa(j),f,iaa,x0aa,widthaa)
        case ('mag-support'); call magsupport(amplaa(j),f,gravz,cs0,rho0)
        case ('arcade-x'); call arcade_x(amplaa(j),f,iaa,kx_aa(j),kz_aa(j))
        case ('halfcos-Bx'); call halfcos_x(amplaa(j),f,iaa)
        case ('uniform-Bx'); call uniform_x(amplaa(j),f,iaa)
        case ('uniform-By'); call uniform_y(amplaa(j),f,iaa)
        case ('uniform-Bz'); call uniform_z(amplaa(j),f,iaa)
        case ('uniform-Bphi'); call uniform_phi(amplaa(j),f,iaa)
        case ('Bz(x)', '3'); call vfield(amplaa(j),f,iaa)
        case ('vfield2'); call vfield2(amplaa(j),f,iaa)
        case ('bipolar'); call bipolar(amplaa(j),f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('bipolar_restzero'); call bipolar_restzero(amplaa(j),f,iaa,kx_aa(j),ky_aa(j))
        case ('vecpatternxy'); call vecpatternxy(amplaa(j),f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('xjump'); call bjump(f,iaa,by_left,by_right,bz_left,bz_right,widthaa,'x')
        case ('sinxsinz'); call sinxsinz(amplaa(j),f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('sinxsinz_Hz'); call sinxsinz(amplaa(j),f,iaa,kx_aa(j),ky_aa(j),kz_aa(j),KKz=kz_aa(j))
        case ('sin2xsin2y'); call sin2x_sin2y_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case ('cosxcosy'); call cosx_cosy_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case ('sinxsiny'); call sinx_siny_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case ('xsiny'); call x_siny_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case ('x1siny'); call x1_siny_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.,phasey_aa(j))
        case ('sinxcosz'); call sinx_siny_cosz(amplaa(j),f,iay,kx_aa(j),ky_aa(j),kz_aa(j))
        case ('sinycosz'); call cosx_siny_cosz(amplaa(j),f,iax,kx_aa(j),ky_aa(j),0.)
        case ('cosysinz'); call cosy_sinz(amplaa(j),f,iax,ky_aa(j),kz_aa(j))
        case ('x3cosycosz'); call x3_cosy_cosz(amplaa(j),f,iax,ky_aa(j),kz_aa(j))
        case ('Ax=cosysinz'); call cosy_sinz(amplaa(j),f,iax,ky_aa(j),kz_aa(j))
        case ('cosxcoscosy'); call cosx_coscosy_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case ('crazy', '5'); call crazy(amplaa(j),f,iaa)
!DM
        case ('strange'); call strange(amplaa(j),f,iaa)
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
        case ('linear-zx')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iay)=-0.5*amplaa(j)*z(n)**2/Lxyz(3)
          enddo; enddo
        case ('Ferriere-uniform-Bx'); call ferriere_uniform_x(amplaa(j),f,iaa)
        case ('Ferriere-uniform-By'); call ferriere_uniform_y(amplaa(j),f,iaa)
!
        case default
!
!  Catch unknown values
!
          call fatal_error('init_aa', &
              'init_aa value "' // trim(initaa(j)) // '" not recognised')

        endselect
!
!  End loop over initial conditions
!
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_aa(f)
!
    endsubroutine init_aa
!***********************************************************************
    subroutine pencil_criteria_magnetic()
!
!   All pencils that the Magnetic module depends on are specified here.
!
!  19-nov-04/anders: coded
!
      use Mpicomm, only: stop_it
!
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_uxb)=.true.
!
      if (dvid/=0.0) lpenc_video(i_b2)=.true.
      if (dvid/=0.0) lpenc_video(i_jb)=.true.
!
!  jj pencil always needed when in Weyl gauge
!
      if (ietat/=0) lpenc_requested(i_del2a)=.true.
      if (lresi_eta_const.or.lresi_eta_shock) &
          lpenc_requested(i_del2a)=.true.
      if (lresi_eta_shock) then
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_diva)=.true.
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
      if (lentropy.and.lohmic_heat) lpenc_requested(i_j2)=.true.
      if (lresi_anomalous) then
        lpenc_requested(i_jj)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
      if (lentropy.or.ldt) lpenc_requested(i_rho1)=.true.
      if (lentropy) lpenc_requested(i_TT1)=.true.
!
!  ambipolar diffusion
!
      if (nu_ni/=0.0) then
        lpenc_requested(i_va2)=.true.
        lpenc_requested(i_jxbrxb)=.true.
        lpenc_requested(i_jxbr2)=.true.
      endif
      if ((lhydro .and. llorentzforce) .or. nu_ni/=0.0) &
          lpenc_requested(i_jxbr)=.true.
      if (nu_ni/=0.0) lpenc_requested(i_va2)=.true.
!
      if (idiag_jxbrxm/=0 .or. idiag_jxbrym/=0 .or. idiag_jxbrzm/=0) &
          lpenc_diagnos(i_jxbr)=.true.
      if (idiag_jxbr2m/=0) lpenc_diagnos(i_jxbr2)=.true.
      if (idiag_jxbrxmx/=0 .or. idiag_jxbrymx/=0 .or. idiag_jxbrzmx/=0 .or. &
          idiag_jxbrxmy/=0 .or. idiag_jxbrymy/=0 .or. idiag_jxbrzmy/=0 .or. &
          idiag_jxbrxmz/=0 .or. idiag_jxbrymz/=0 .or. idiag_jxbrzmz/=0) &
          lpenc_diagnos(i_jxbr)=.true.
!
      if (idiag_b2ruzm/=0) &
          lpenc_diagnos(i_rho)=.true.
!
      if (idiag_axmxz/=0 .or. idiag_aymxz/=0 .or. idiag_azmxz/=0) &
           lpenc_diagnos2d(i_aa)=.true.
!
      if (idiag_aybym2/=0 .or. idiag_exaym2/=0 .or. &
          idiag_examx/=0 .or. idiag_examy/=0 .or. idiag_examz/=0 .or. &
          idiag_examz1/=0 .or. idiag_examz2/=0 .or. idiag_examz3/=0 &
         ) lpenc_diagnos(i_aa)=.true.
      if (idiag_arms/=0 .or. idiag_amax/=0) lpenc_diagnos(i_a2)=.true.
      if (idiag_ab_int/=0 .or. idiag_abm/=0 .or. idiag_abmh/=0 &
          .or. idiag_abmz/=0 .or. idiag_abrms/=0 &
          .or. idiag_abumx/=0 .or. idiag_abumy/=0 .or. idiag_abumz/=0 &
          .or. idiag_abuxmz/=0 .or. idiag_abuymz/=0 .or. idiag_abuzmz/=0 &
         ) lpenc_diagnos(i_ab)=.true.
      if (idiag_uam/=0 .or. idiag_uamz/=0) lpenc_diagnos(i_ua)=.true.
      if (idiag_djuidjbim/=0 .or. idiag_b3b21m/=0 &
          .or. idiag_dexbmx/=0 .or. idiag_dexbmy/=0 .or. idiag_dexbmz/=0 &
          .or. idiag_b1b32m/=0 .or.  idiag_b2b13m/=0) &
          lpenc_diagnos(i_bij)=.true.
      if (idiag_j2m/=0 .or. idiag_jm2/=0 .or. idiag_jrms/=0 .or. &
          idiag_jmax/=0 .or. idiag_epsM/=0 .or. idiag_epsM_LES/=0 .or. &
          idiag_ajm/=0 ) &
          lpenc_diagnos(i_j2)=.true.
      if (idiag_epsAD/=0) lpenc_diagnos(i_jxbr2)=.true.
      if (idiag_jb_int/=0 .or. idiag_jbm/=0 .or. idiag_jbmz/=0) lpenc_diagnos(i_jb)=.true.
      if (idiag_vArms/=0 .or. idiag_vAmax/=0 .or. idiag_vA2m/=0) lpenc_diagnos(i_va2)=.true.
      if (idiag_cosubm/=0) lpenc_diagnos(i_cosub)=.true.
      if (idiag_ubm/=0 .or. idiag_ubbzm/=0) lpenc_diagnos(i_ub)=.true.
      if (idiag_djuidjbim/=0 .or. idiag_uxDxuxbm/=0) lpenc_diagnos(i_uij)=.true.
      if (idiag_uxjm/=0) lpenc_diagnos(i_uxj)=.true.
      if (idiag_vArms/=0 .or. idiag_vAmax/=0) lpenc_diagnos(i_va2)=.true.
      if (idiag_uxBrms/=0 .or. idiag_Rmrms/=0) lpenc_diagnos(i_uxb2)=.true.
      if (idiag_beta1m/=0 .or. idiag_beta1max/=0) lpenc_diagnos(i_beta)=.true.
      if (idiag_djuidjbim/=0) lpenc_diagnos(i_djuidjbi)=.true.
      if (idiag_b2divum/=0) lpenc_diagnos(i_divu)=.true.
      if (idiag_b2divum/=0) lpenc_diagnos(i_b2)=.true.
      if (idiag_ujxbm/=0) lpenc_diagnos(i_ujxb)=.true.
      if (idiag_gpxbm/=0) lpenc_diagnos(i_glnrhoxb)=.true.
      if (idiag_jxbxbm/=0) lpenc_diagnos(i_jxbxb)=.true.
      if (idiag_oxuxbm/=0) lpenc_diagnos(i_oxuxb)=.true.
      if (idiag_exaym2/=0 .or. idiag_exjm2/=0 &
          .or. idiag_jmx/=0 .or. idiag_jmy/=0 .or. idiag_jmz/=0 &
          .or. idiag_ambmz/=0 .or. idiag_jmbmz/=0 .or. idiag_kmz/=0 &
          .or. idiag_examx/=0 .or. idiag_examy/=0 .or. idiag_examz/=0 &
          .or. idiag_examz1/=0 .or. idiag_examz2/=0 .or. idiag_examz3/=0 &
          .or. idiag_exjmx/=0 .or. idiag_exjmy/=0 .or. idiag_exjmz/=0 &
         ) lpenc_diagnos(i_jj)=.true.
      if (idiag_examz1/=0 .or. idiag_examz2/=0 .or. idiag_examz3/=0 &
         ) lpenc_diagnos(i_exa)=.true.
      if (idiag_phibmx/=0 .or. idiag_phibmy/=0 .or. idiag_phibmz/=0 &
         ) lpenc_diagnos(i_diva)=.true.
      if (idiag_b2uzm/=0 .or. idiag_b2ruzm/=0 .or. &
          idiag_b1m/=0 .or. idiag_b2m/=0 .or. idiag_bm2/=0 .or. &
          idiag_brmsh/=0 .or. idiag_brmsn/=0 .or. idiag_brmss/=0 .or. &
          idiag_brms/=0 .or. idiag_bmax/=0 .or. &
          idiag_emag/=0 ) &
          lpenc_diagnos(i_b2)=.true.
      if (idiag_etavamax/=0) lpenc_diagnos(i_etava)=.true.
      if (idiag_etajmax/=0) lpenc_diagnos(i_etaj)=.true.
      if (idiag_etaj2max/=0) lpenc_diagnos(i_etaj2)=.true.
      if (idiag_etajrhomax/=0) lpenc_diagnos(i_etajrho)=.true.
      if (idiag_cosjbm/=0) lpenc_diagnos(i_cosjb)=.true.
      if (idiag_cosubm/=0) then
        lpenc_diagnos(i_cosub)=.true.
      endif
      if ((idiag_jparallelm/=0).or.(idiag_jperpm/=0)) then
        lpenc_diagnos(i_cosjb)=.true.
        lpenc_diagnos(i_jparallel)=.true.
        lpenc_diagnos(i_jperp)=.true.
      endif
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
      else
      endif
!
    endsubroutine pencil_criteria_magnetic
!***********************************************************************
    subroutine pencil_interdep_magnetic(lpencil_in)
!
!  Interdependency among pencils from the Magnetic module is specified here.
!
!  19-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_jparallel).or.lpencil_in(i_jperp)) then
        lpencil_in(i_cosjb)=.true.
        lpencil_in(i_jxb)=.true.
      endif
      if (lpencil_in(i_cosjb)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_j2)=.true.
        lpencil_in(i_jb)=.true.
      endif
      if (lpencil_in(i_a2)) lpencil_in(i_aa)=.true.
      if (lpencil_in(i_ab)) then
        lpencil_in(i_aa)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_ua)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_aa)=.true.
      endif
      if (lpencil_in(i_va2)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_rho1)=.true.
      endif
      if (lpencil_in(i_etava)) lpencil_in(i_va2)=.true.
      if (lpencil_in(i_etaj) .or. lpencil_in(i_etaj2) .or. lpencil_in(i_etajrho)) then
        lpencil_in(i_j2)=.true.
        lpencil_in(i_rho1)=.true.
      endif
      if (lpencil_in(i_j2)) lpencil_in(i_jj)=.true.
      if (lpencil_in(i_uxj)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_jj)=.true.
      endif
      if (lpencil_in(i_jb)) then
        lpencil_in(i_bb)=.true.
        lpencil_in(i_jj)=.true.
      endif
      if (lpencil_in(i_jxbr) .and. va2max_jxb>0) lpencil_in(i_va2)=.true.
      if (lpencil_in(i_jxbr)) then
        lpencil_in(i_jxb)=.true.
        lpencil_in(i_rho1)=.true.
      endif
      if (lpencil_in(i_jxb)) then
        lpencil_in(i_jj)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_uxb2)) lpencil_in(i_uxb)=.true.
      if (lpencil_in(i_uxb)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_cosub)) then
        lpencil_in(i_ub)=.true.
        lpencil_in(i_u2)=.true.
        lpencil_in(i_b2)=.true.
      endif
      if (lpencil_in(i_ub)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_beta)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_pp)=.true.
      endif
      if (lpencil_in(i_b2)) lpencil_in(i_bb)=.true.
      if (lpencil_in(i_jj)) lpencil_in(i_bij)=.true.
      if (lpencil_in(i_bb)) then
        lpencil_in(i_aij)=.true.
      endif
      if (lpencil_in(i_djuidjbi)) then
        lpencil_in(i_uij)=.true.
        lpencil_in(i_bij)=.true.
      endif
      if (lpencil_in(i_jo)) then
        lpencil_in(i_oo)=.true.
        lpencil_in(i_jj)=.true.
      endif
      if (lpencil_in(i_ujxb)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_jxb)=.true.
      endif
      if (lpencil_in(i_oxu)) then
        lpencil_in(i_oo)=.true.
        lpencil_in(i_uu)=.true.
      endif
      if (lpencil_in(i_oxuxb)) then
        lpencil_in(i_oxu)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_jxbxb)) then
        lpencil_in(i_jxb)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_jxbrxb)) then
        lpencil_in(i_jxbr)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_glnrhoxb)) then
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_oxj)) then
        lpencil_in(i_oo)=.true.
        lpencil_in(i_jj)=.true.
      endif
      if (lpencil_in(i_jij)) lpencil_in(i_bij)=.true.
      if (lpencil_in(i_sj)) then
        lpencil_in(i_sij)=.true.
        lpencil_in(i_jij)=.true.
      endif
      if (lpencil_in(i_ss12)) lpencil_in(i_sij)=.true.
      if (lpencil_in(i_del2A)) then
        lpencil_in(i_jj)=.true.
        lpencil_in(i_graddivA)=.true.
      endif
      if (lpencil_in(i_ss12)) lpencil_in(i_sj)=.true.
!
    endsubroutine pencil_interdep_magnetic
!***********************************************************************
    subroutine calc_pencils_magnetic(f,p)
!
!  Calculate Magnetic pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  19-nov-04/anders: coded
!
      use Sub
      use Deriv
      use Diagnostics, only: sum_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
!      real, dimension (nx,3) :: bb_ext_pot
      real, dimension (nx) :: rho1_jxb
      real, dimension (nx) :: jcrossb2
      real :: B2_ext,c,s
      integer :: i,j,ix
!
      intent(inout) :: f,p
! aa
      if (lpencil(i_aa)) p%aa=f(l1:l2,m,n,iax:iaz)
! a2
      if (lpencil(i_a2)) call dot2_mn(p%aa,p%a2)
! aij
      if (lpencil(i_aij)) call gij(f,iaa,p%aij,1)
! diva
      if (lpencil(i_diva)) call div_mn(p%aij,p%diva)
! bb
      if (lpencil(i_bb)) then
        call curl_mn(p%aij,p%bb)
!
!  Save field before adding imposed field (for diagnostics).
!
        p%bbb=p%bb
        B2_ext=B_ext(1)**2+B_ext(2)**2+B_ext(3)**2
!
!  Allow external field to precess about z-axis with frequency omega_Bz_ext.
!
        if (B2_ext/=0.0) then
!
!  luse_curvilinear_bext is default. The B_ext the user defines in
!  magnetic_init_pars respects the coordinate system of preference
!  which means that B_ext=(0.0,1.0,0.0) is an azimuthal field in cylindrical
!  coordinates and a polar one in spherical.
!
          if (omega_Bz_ext==0.0) then
            B_ext_tmp=B_ext
          elseif (omega_Bz_ext/=0.0) then
            c=cos(omega_Bz_ext*t)
            s=sin(omega_Bz_ext*t)
            B_ext_tmp(1)=B_ext(1)*c-B_ext(2)*s
            B_ext_tmp(2)=B_ext(1)*s+B_ext(2)*c
            B_ext_tmp(3)=B_ext(3)
          endif
!
!  Add the external field.
!
          if (B_ext_tmp(1)/=0.0) p%bb(:,1)=p%bb(:,1)+B_ext_tmp(1)
          if (B_ext_tmp(2)/=0.0) p%bb(:,2)=p%bb(:,2)+B_ext_tmp(2)
          if (B_ext_tmp(3)/=0.0) p%bb(:,3)=p%bb(:,3)+B_ext_tmp(3)
          if (headtt) print*,'calc_pencils_magnetic: B_ext=',B_ext
          if (headtt) print*,'calc_pencils_magnetic: B_ext_tmp=',B_ext_tmp
        endif
!
!  Add the external potential field.
!
        if (lB_ext_pot) then
!          call get_global(bb_ext_pot,m,n,'B_ext_pot')
!          p%bb=p%bb+bb_ext_pot
        endif
!
!  Add external B-field.
!
        if (iglobal_bx_ext/=0) p%bb(:,1)=p%bb(:,1)+f(l1:l2,m,n,iglobal_bx_ext)
        if (iglobal_by_ext/=0) p%bb(:,2)=p%bb(:,2)+f(l1:l2,m,n,iglobal_by_ext)
        if (iglobal_bz_ext/=0) p%bb(:,3)=p%bb(:,3)+f(l1:l2,m,n,iglobal_bz_ext)
      endif
!
! b2 (default is that B_ext is not included), but this can be changed
! by setting luse_Bext_in_b2=.true.
!
      if (luse_Bext_in_b2) then
        if (lpencil(i_b2)) call dot2_mn(p%bb,p%b2)
      else
        if (lpencil(i_b2)) call dot2_mn(p%bbb,p%b2)
      endif
! ab
      if (lpencil(i_ab)) call dot_mn(p%aa,p%bbb,p%ab)
      if (lpencil(i_ua)) call dot_mn(p%uu,p%aa,p%ua)
! uxb
      if (lpencil(i_uxb)) then
        call cross_mn(p%uu,p%bb,p%uxb)
        call cross_mn(p%uu,p%bbb,uxbb)
!  add external e-field.
        if (iglobal_ex_ext/=0) p%uxb(:,1)=p%uxb(:,1)+f(l1:l2,m,n,iglobal_ex_ext)
        if (iglobal_ey_ext/=0) p%uxb(:,2)=p%uxb(:,2)+f(l1:l2,m,n,iglobal_ey_ext)
        if (iglobal_ez_ext/=0) p%uxb(:,3)=p%uxb(:,3)+f(l1:l2,m,n,iglobal_ez_ext)
      endif
! exa
      if (lpencil(i_exa)) call cross_mn(-p%uxb+eta*p%jj,p%aa,p%exa)
!
!  bij, del2a, graddiva
!  For non-cartesian coordinates jj is always required for del2a=graddiva-jj
!
      if (lpencil(i_bij) .or. lpencil(i_del2a) .or. lpencil(i_graddiva) .or. &
          lpencil(i_jj) ) then
        call gij_etc(f,iaa,p%bij,p%del2a,p%graddiva)
        if (.not. lpencil(i_bij)) p%bij=0.0      ! Avoid warnings from pencil
        if (.not. lpencil(i_del2A)) p%del2A=0.0  ! consistency check...
        if (.not. lpencil(i_graddiva)) p%graddiva=0.0
        if (lpencil(i_jj)) call curl_mn(p%bij,p%jj)
      endif
! jj
      if (lpencil(i_jj)) then
        p%jj=mu01*p%jj
!
!  Add external j-field.
!
        if (iglobal_jx_ext/=0) p%jj(:,1)=p%jj(:,1)+f(l1:l2,m,n,iglobal_jx_ext)
        if (iglobal_jy_ext/=0) p%jj(:,2)=p%jj(:,2)+f(l1:l2,m,n,iglobal_jy_ext)
        if (iglobal_jz_ext/=0) p%jj(:,3)=p%jj(:,3)+f(l1:l2,m,n,iglobal_jz_ext)
      endif
! j2
      if (lpencil(i_j2)) call dot2_mn(p%jj,p%j2)
! jb
      if (lpencil(i_jb)) call dot_mn(p%jj,p%bbb,p%jb)
! va2
      if (lpencil(i_va2)) then
        p%va2=p%b2*mu01*p%rho1
        if (lcheck_positive_va2 .and. minval(p%va2)<0.0) then
          print*, 'calc_pencils_magnetic: Alfven speed is imaginary!'
          print*, 'calc_pencils_magnetic: it, itsub, iproc=', it, itsub, iproc
          print*, 'calc_pencils_magnetic: m, y(m), n, z(n)=', m, y(m), n, z(n)
          p%va2=abs(p%va2)
        endif
      endif
! eta_va
      if (lpencil(i_etava)) then
        p%etava = mu0 * eta_va * dxmax * sqrt(p%va2)
        if (eta_min > 0.) where (p%etava < eta_min) p%etava = 0.
      endif
! eta_j
      if (lpencil(i_etaj)) then
        p%etaj = mu0 * eta_j * dxmax**2 * sqrt(mu0 * p%j2 * p%rho1)
        if (eta_min > 0.) where (p%etaj < eta_min) p%etaj = 0.
      endif
! eta_j2
      if (lpencil(i_etaj2)) then
        p%etaj2 = etaj20 * p%j2 * p%rho1
        if (eta_min > 0.) where (p%etaj2 < eta_min) p%etaj2 = 0.
      endif
! eta_jrho
      if (lpencil(i_etajrho)) then
        p%etajrho = mu0 * eta_jrho * dxmax * sqrt(p%j2) * p%rho1
        if (eta_min > 0.) where (p%etajrho < eta_min) p%etajrho = 0.
      endif
! jxb
      if (lpencil(i_jxb)) call cross_mn(p%jj,p%bb,p%jxb)
! jxbr
      if (lpencil(i_jxbr)) rho1_jxb=p%rho1
! cosjb
      if (lpencil(i_cosjb)) then
        do ix=1,nx
          if ((abs(p%j2(ix)).le.tini).or.(abs(p%b2(ix)).le.tini))then
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
      if (lpencil(i_jparallel).or.lpencil(i_jperp)) then
        p%jparallel=sqrt(p%j2)*p%cosjb
        call dot2_mn(p%jxb,jcrossb2)
        do ix=1,nx
          if ((abs(p%j2(ix)).le.tini).or.(abs(p%b2(ix)).le.tini))then
            p%jperp=0
          else
            p%jperp=sqrt(jcrossb2(ix))/sqrt(p%b2(ix))
          endif
        enddo
      endif
! jxbr
      if (lpencil(i_jxbr)) then
        rho1_jxb=p%rho1
!
!  Set rhomin_jxb>0 in order to limit the jxb term at very low densities.
!  Set va2max_jxb>0 in order to limit the jxb term at very high Alfven speeds.
!  Set va2power_jxb to an integer value in order to specify the power of the
!  limiting term,
!
        if (rhomin_jxb>0) rho1_jxb=min(rho1_jxb,1/rhomin_jxb)
        if (va2max_jxb>0) then
          rho1_jxb = rho1_jxb &
                   * (1+(p%va2/va2max_jxb)**va2power_jxb)**(-1.0/va2power_jxb)
        endif
        call multsv_mn(rho1_jxb,p%jxb,p%jxbr)
      endif
! jxbr2
      if (lpencil(i_jxbr2)) call dot2_mn(p%jxbr,p%jxbr2)
! ub
      if (lpencil(i_ub)) call dot_mn(p%uu,p%bb,p%ub)
! cosub
      if (lpencil(i_cosub)) then
        do ix=1,nx
          if ((abs(p%u2(ix)).le.tini).or.(abs(p%b2(ix)).le.tini)) then
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
      if (lpencil(i_uxb2)) call dot2_mn(p%uxb,p%uxb2)
! uxj
      if (lpencil(i_uxj)) call cross_mn(p%uu,p%jj,p%uxj)
! beta
      if (lpencil(i_beta)) p%beta=0.5*p%b2/p%pp
! djuidjbi
      if (lpencil(i_djuidjbi)) call multmm_sc(p%uij,p%bij,p%djuidjbi)
! jo
      if (lpencil(i_jo)) call dot(p%jj,p%oo,p%jo)
! ujxb
      if (lpencil(i_ujxb)) call dot_mn(p%uu,p%jxb,p%ujxb)
! oxu
      if (lpencil(i_oxu)) call cross_mn(p%oo,p%uu,p%oxu)
! oxuxb
      if (lpencil(i_oxuxb)) call cross_mn(p%oxu,p%bb,p%oxuxb)
! jxbxb
      if (lpencil(i_jxbxb)) call cross_mn(p%jxb,p%bb,p%jxbxb)
! jxbrxb
      if (lpencil(i_jxbrxb)) call cross_mn(p%jxbr,p%bb,p%jxbrxb)
! glnrhoxb
      if (lpencil(i_glnrhoxb)) call cross_mn(p%glnrho,p%bb,p%glnrhoxb)
! oxj
      if (lpencil(i_oxj)) call cross_mn(p%oo,p%jj,p%oxJ)
! jij
      if (lpencil(i_jij)) then
        do j=1,3
          do i=1,3
            p%jij(:,i,j)=.5*(p%bij(:,i,j)+p%bij(:,j,i))
          enddo
        enddo
      endif
! sj
      if (lpencil(i_sj)) call multmm_sc(p%sij,p%jij,p%sj)
! ss12
      if (lpencil(i_ss12)) p%ss12=sqrt(abs(p%sj))
!
!  Store bb in auxiliary variable if requested.
!  Just neccessary immediately before writing snapshots, but how would we
!  know we are?
!
     if (lbb_as_aux) f(l1:l2,m,n,ibx:ibz)=p%bb
     if (ljj_as_aux) f(l1:l2,m,n,ijx:ijz)=p%jj
!
    endsubroutine calc_pencils_magnetic
!***********************************************************************
    subroutine daa_dt(df,p)
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
!  26-may-04/axel: ambipolar diffusion added
!
      use Diagnostics
      use EquationOfState, only: eoscalc,gamma_m1,cs0
      use Io, only: output_pencil
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: uxDxuxb,fres
      real, dimension (nx,3) :: exj,dexb,phib,jxbb
      real, dimension (nx) :: exabot,exatop
      real, dimension (nx) :: jxb_dotB0,uxb_dotB0
      real, dimension (nx) :: oxuxb_dotB0,jxbxb_dotB0,uxDxuxb_dotB0
      real, dimension (nx) :: uj,aj,phi
      real, dimension (nx) :: uxj_dotB0,b3b21,b1b32,b2b13
      real, dimension (nx) :: rho1_jxb
      real, dimension (nx) :: B1dot_glnrhoxb
      real, dimension (nx) :: eta_smag,etatotal
      real, dimension (nx) :: fres2
      real, parameter :: OmegaSS=1.0
      integer :: i
      integer, parameter :: nxy=nxgrid*nygrid
!
      intent(inout)  :: p
      intent(inout)  :: df
!
!  Identify module and boundary conditions.
!
      call timing('daa_dt','entered',mnloop=.true.)
      if (headtt.or.ldebug) print*,'daa_dt: SOLVE'
      if (headtt) then
        call identify_bcs('Ax',iax)
        call identify_bcs('Ay',iay)
        call identify_bcs('Az',iaz)
      endif
!
!  Add jxb/rho to momentum equation.
!
      if (lhydro.and.llorentzforce) &
           df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+p%jxbr
!
!  Restivivity term
!
!  Because of gauge invariance, we can add the gradient of an arbitrary scalar
!  field Phi to the induction equation without changing the magnetic field,
!    dA/dt = u x B - eta j + grad(Phi).
!
!  If lweyl_gauge=T, we choose Phi = const. and solve
!  Else, if lweyl_gauge=F, we make the gauge choice Phi = eta div(A)
!  and thus solve
!    dA/dt = u x B + eta laplace(A) + div(A) grad(eta).
!
!  Note: lweyl_gauge=T is so far only implemented for some resistivity types.
!
      if (headtt) print*, 'daa_dt: iresistivity=', iresistivity
!
      fres=0.0
      etatotal=0.0
!
      if (lresi_eta_const) then
        fres=fres+eta*p%del2a
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta
        etatotal=etatotal+eta
      endif
!
      if (lresi_eta_shock) then
        do i=1,3
          fres(:,i)=fres(:,i)+ &
               eta_shock*(p%shock*p%del2a(:,i)+p%diva*p%gshock(:,i))
        enddo
        if (lfirst.and.ldt) diffus_eta=diffus_eta+eta_shock*p%shock
        etatotal=etatotal+eta_shock*p%shock
      endif
!
!  Multiply resistivity by Nyquist scale, for resistive time-step.
!
      if (lfirst.and.ldt) then
        diffus_eta =diffus_eta*dxyz_2
        if (headtt.or.ldebug) then
          print*, 'daa_dt: max(diffus_eta)  =', maxval(diffus_eta)
        endif
      endif
!
!  Add Ohmic heat to entropy or temperature equation.
!
      if (lohmic_heat.and.lentropy) &
           df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + &
           etatotal*mu0*p%j2*p%rho1*p%TT1
!
!  Induction equation.
!
      if (linduction) &
           df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+p%uxb+fres
!
!  Add ``va^2/dx^2'' contribution to timestep.
!  Consider advective timestep only when lhydro=T.
!
      if (lfirst.and.ldt) then
        if (lhydro) then
          rho1_jxb=p%rho1
          if (rhomin_jxb>0) rho1_jxb=min(rho1_jxb,1/rhomin_jxb)
          if (va2max_jxb>0) then
            rho1_jxb = rho1_jxb &
                     * (1+(p%va2/va2max_jxb)**va2power_jxb)**(-1.0/va2power_jxb)
          endif
          advec_va2=((p%bb(:,1)*dx_1(l1:l2))**2+ &
               (p%bb(:,2)*dy_1(  m  ))**2+ &
               (p%bb(:,3)*dz_1(  n  ))**2)*mu01*rho1_jxb
        endif
      endif
!
!  Calculate diagnostic quantities.
!
      if (ldiagnos) then
        if (idiag_beta1m/=0) call sum_mn_name(p%beta,idiag_beta1m)
        if (idiag_beta1max/=0) call max_mn_name(p%beta,idiag_beta1max)
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
        if (idiag_b2m/=0) call sum_mn_name(p%b2,idiag_b2m)
        if (idiag_bm2/=0) call max_mn_name(p%b2,idiag_bm2)
        if (idiag_brms/=0) call sum_mn_name(p%b2,idiag_brms,lsqrt=.true.)
        if (idiag_emag/=0) call integrate_mn_name(mu012*p%b2,idiag_emag)
        if (idiag_brmsh/=0) then
          if (lequatory) call sum_mn_name_halfy(p%b2,idiag_brmsh)
          if (lequatorz) call sum_mn_name_halfz(p%b2,idiag_brmsh)
          fname(idiag_brmsn)=fname_half(idiag_brmsh,1)
          fname(idiag_brmss)=fname_half(idiag_brmsh,2)
          itype_name(idiag_brmsn)=ilabel_sum_sqrt
          itype_name(idiag_brmss)=ilabel_sum_sqrt
        endif
        if (idiag_bmax/=0) call max_mn_name(p%b2,idiag_bmax,lsqrt=.true.)
        if (idiag_bxmin/=0) call max_mn_name(-p%bb(:,1),idiag_bxmin,lneg=.true.)
        if (idiag_bymin/=0) call max_mn_name(-p%bb(:,2),idiag_bymin,lneg=.true.)
        if (idiag_bzmin/=0) call max_mn_name(-p%bb(:,3),idiag_bzmin,lneg=.true.)
        if (idiag_bxmax/=0) call max_mn_name(p%bb(:,1),idiag_bxmax)
        if (idiag_bymax/=0) call max_mn_name(p%bb(:,2),idiag_bymax)
        if (idiag_bzmax/=0) call max_mn_name(p%bb(:,3),idiag_bzmax)
        if (idiag_aybym2/=0) &
            call sum_mn_name(2*p%aa(:,2)*p%bb(:,2),idiag_aybym2)
        if (idiag_abm/=0) call sum_mn_name(p%ab,idiag_abm)
        if (idiag_abumx/=0) call sum_mn_name(p%uu(:,1)*p%ab,idiag_abumx)
        if (idiag_abumy/=0) call sum_mn_name(p%uu(:,2)*p%ab,idiag_abumy)
        if (idiag_abumz/=0) call sum_mn_name(p%uu(:,3)*p%ab,idiag_abumz)
        if (idiag_abrms/=0) call sum_mn_name(p%ab**2,idiag_abrms,lsqrt=.true.)
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
!  Cross helicity (linkage between vortex tubes and flux tubes).
!
        if (idiag_ubm/=0) call sum_mn_name(p%ub,idiag_ubm)
        if (idiag_uxbxm/=0) call sum_mn_name(p%uu(:,1)*p%bb(:,1),idiag_uxbxm)
        if (idiag_uybym/=0) call sum_mn_name(p%uu(:,2)*p%bb(:,2),idiag_uybym)
        if (idiag_uzbzm/=0) call sum_mn_name(p%uu(:,3)*p%bb(:,3),idiag_uzbzm)
        if (idiag_cosubm/=0) call sum_mn_name(p%cosub,idiag_cosubm)
!
!  Field-velocity cross helicity (linkage between velocity and magnetic tubes).
!
        if (idiag_uam/=0) call sum_mn_name(p%ua,idiag_uam)
!
!  Current-vortex cross helicity (linkage between vortex and current tubes).
!
        if (idiag_ujm/=0) then
          call dot (p%uu,p%jj,uj)
          call sum_mn_name(uj,idiag_ujm)
        endif
!
!  rms values of mean field.
!
        if (idiag_bxm/=0) call sum_mn_name(p%bbb(:,1),idiag_bxm)
        if (idiag_bym/=0) call sum_mn_name(p%bbb(:,2),idiag_bym)
        if (idiag_bzm/=0) call sum_mn_name(p%bbb(:,3),idiag_bzm)
        if (idiag_bx2m/=0) call sum_mn_name(p%bbb(:,1)**2,idiag_bx2m)
        if (idiag_by2m/=0) call sum_mn_name(p%bbb(:,2)**2,idiag_by2m)
        if (idiag_bz2m/=0) call sum_mn_name(p%bbb(:,3)**2,idiag_bz2m)
        if (idiag_bxbym/=0) call sum_mn_name(p%bbb(:,1)*p%bbb(:,2),idiag_bxbym)
        if (idiag_bxbzm/=0) call sum_mn_name(p%bbb(:,1)*p%bbb(:,3),idiag_bxbzm)
        if (idiag_bybzm/=0) call sum_mn_name(p%bbb(:,2)*p%bbb(:,3),idiag_bybzm)
        if (idiag_djuidjbim/=0) call sum_mn_name(p%djuidjbi,idiag_djuidjbim)
!
!  Calculate B*sin(phi) = -<Bx*sinkz> + <By*coskz>.
!
        if (idiag_bsinphz/=0) call sum_mn_name(-p%bbb(:,1)*sinkz(n)+p%bbb(:,2)*coskz(n),idiag_bsinphz)
!
!  Calculate B*cos(phi) = <Bx*coskz> + <By*sinkz>.
!
        if (idiag_bcosphz/=0) call sum_mn_name(+p%bbb(:,1)*coskz(n)+p%bbb(:,2)*sinkz(n),idiag_bcosphz)
!
!  Magnetic field components at one point (=pt).
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_bxpt/=0) call save_name(p%bb(lpoint-nghost,1),idiag_bxpt)
          if (idiag_bypt/=0) call save_name(p%bb(lpoint-nghost,2),idiag_bypt)
          if (idiag_bzpt/=0) call save_name(p%bb(lpoint-nghost,3),idiag_bzpt)
          if (idiag_Expt/=0) call save_name(uxbb(lpoint-nghost,1),idiag_Expt)
          if (idiag_Eypt/=0) call save_name(uxbb(lpoint-nghost,2),idiag_Eypt)
          if (idiag_Ezpt/=0) call save_name(uxbb(lpoint-nghost,3),idiag_Ezpt)
        endif
!
!  v_A = |B|/sqrt(rho); in units where mu_0=1
!
        if (idiag_vA2m/=0)  call sum_mn_name(p%va2,idiag_vA2m)
        if (idiag_vArms/=0) call sum_mn_name(p%va2,idiag_vArms,lsqrt=.true.)
        if (idiag_vAmax/=0) call max_mn_name(p%va2,idiag_vAmax,lsqrt=.true.)
        if (idiag_dtb/=0) &
            call max_mn_name(sqrt(advec_va2)/cdt,idiag_dtb,l_dt=.true.)
!
!  Lorentz force.
!
        if (idiag_jxbrxm/=0) call sum_mn_name(p%jxbr(:,1),idiag_jxbrxm)
        if (idiag_jxbrym/=0) call sum_mn_name(p%jxbr(:,2),idiag_jxbrym)
        if (idiag_jxbrzm/=0) call sum_mn_name(p%jxbr(:,3),idiag_jxbrzm)
        if (idiag_jxbr2m/=0) call sum_mn_name(p%jxbr2,idiag_jxbr2m)
!
!  <J.A> for calculating k_effective, for example.
!
        if (idiag_ajm/=0) then
          call dot (p%aa,p%jj,aj)
          call sum_mn_name(aj,idiag_ajm)
        endif
!
!  Output kx_aa for calculating k_effective.
!
        if (idiag_kx_aa/=0) call save_name(kx_aa(1),idiag_kx_aa)
!
!  Helicity integrals.
!
        if (idiag_ab_int/=0) call integrate_mn_name(p%ab,idiag_ab_int)
        if (idiag_jb_int/=0) call integrate_mn_name(p%jb,idiag_jb_int)
!
! <J.B>
!
        if (idiag_jbm/=0) call sum_mn_name(p%jb,idiag_jbm)
        if (idiag_j2m/=0) call sum_mn_name(p%j2,idiag_j2m)
        if (idiag_jm2/=0) call max_mn_name(p%j2,idiag_jm2)
        if (idiag_jrms/=0) call sum_mn_name(p%j2,idiag_jrms,lsqrt=.true.)
        if (idiag_jmax/=0) call max_mn_name(p%j2,idiag_jmax,lsqrt=.true.)
        if (idiag_epsM_LES/=0) call sum_mn_name(eta_smag*p%j2,idiag_epsM_LES)
        if (idiag_dteta/=0)  call max_mn_name(diffus_eta/cdtv,idiag_dteta,l_dt=.true.)
        if (idiag_cosjbm/=0) call sum_mn_name(p%cosjb,idiag_cosjbm)
        if (idiag_jparallelm/=0) call sum_mn_name(p%jparallel,idiag_jparallelm)
        if (idiag_jperpm/=0) call sum_mn_name(p%jperp,idiag_jperpm)
!
!  Resistivity.
!
        if (idiag_etasmagm/=0)   call sum_mn_name(eta_smag,idiag_etasmagm)
        if (idiag_etasmagmin/=0) call max_mn_name(-eta_smag,idiag_etasmagmin,lneg=.true.)
        if (idiag_etasmagmax/=0) call max_mn_name(eta_smag,idiag_etasmagmax)
        if (idiag_etavamax/=0) call max_mn_name(p%etava,idiag_etavamax)
        if (idiag_etajmax/=0) call max_mn_name(p%etaj,idiag_etajmax)
        if (idiag_etaj2max/=0) call max_mn_name(p%etaj2,idiag_etaj2max)
        if (idiag_etajrhomax/=0) call max_mn_name(p%etajrho,idiag_etajrhomax)
!
!  Not correct for hyperresistivity:
!
        if (idiag_epsM/=0) call sum_mn_name(eta*p%j2,idiag_epsM)
!
!  Heating by ion-neutrals friction.
!
        if (idiag_epsAD/=0) call sum_mn_name(nu_ni1*p%rho*p%jxbr2,idiag_epsAD)
!
!  <A>'s, <A^2> and A^2|max
!
        if (idiag_axm/=0) call sum_mn_name(p%aa(:,1),idiag_axm)
        if (idiag_aym/=0) call sum_mn_name(p%aa(:,2),idiag_aym)
        if (idiag_azm/=0) call sum_mn_name(p%aa(:,3),idiag_azm)
        if (idiag_arms/=0) call sum_mn_name(p%a2,idiag_arms,lsqrt=.true.)
        if (idiag_amax/=0) call max_mn_name(p%a2,idiag_amax,lsqrt=.true.)
!
        if (idiag_uxbm /=0 .or. & 
            idiag_uxbmx/=0 .or. idiag_uxbmy/=0 .or. idiag_uxbmz/=0) then
          call dot(B_ext_inv,p%uxb,uxb_dotB0)
          if (idiag_uxbm/=0) call sum_mn_name(uxb_dotB0,idiag_uxbm)
          if (idiag_uxbmx/=0) call sum_mn_name(uxbb(:,1),idiag_uxbmx)
          if (idiag_uxbmy/=0) call sum_mn_name(uxbb(:,2),idiag_uxbmy)
          if (idiag_uxbmz/=0) call sum_mn_name(uxbb(:,3),idiag_uxbmz)
        endif
!
!  Calculate part I of magnetic helicity flux (ExA contribution).
!
        if (idiag_examx/=0 .or. idiag_examy/=0 .or. idiag_examz/=0 .or. &
            idiag_exatop/=0 .or. idiag_exabot/=0) then
          if (idiag_examx/=0) call sum_mn_name(p%exa(:,1),idiag_examx)
          if (idiag_examy/=0) call sum_mn_name(p%exa(:,2),idiag_examy)
          if (idiag_examz/=0) call sum_mn_name(p%exa(:,3),idiag_examz)
!
          if (idiag_exabot/=0) then
            if (z(n)==xyz0(3)) then
              exabot=p%exa(:,3)
            else
              exabot=0.
            endif
            call integrate_mn_name(exabot,idiag_exabot)
          endif
!
          if (idiag_exatop/=0) then
            if (z(n)==xyz1(3)) then
              exatop=p%exa(:,3)
            else
              exatop=0.
            endif
            call integrate_mn_name(exatop,idiag_exatop)
          endif
!
        endif
!
!  Calculate part II of magnetic helicity flux (phi*B contribution).
!
        if (idiag_phibmx/=0 .or. idiag_phibmy/=0 .or. idiag_phibmz/=0) then
          phi=eta*p%diva
          call multvs(p%bb,phi,phib)
          if (idiag_phibmx/=0) call sum_mn_name(phib(:,1),idiag_phibmx)
          if (idiag_phibmy/=0) call sum_mn_name(phib(:,2),idiag_phibmy)
          if (idiag_phibmz/=0) call sum_mn_name(phib(:,3),idiag_phibmz)
        endif
!
!  Calculate part I of current helicity flux (for imposed field).
!
        if (idiag_exjmx/=0 .or. idiag_exjmy/=0 .or. idiag_exjmz/=0) then
          call cross_mn(-p%uxb+eta*p%jj,p%jj,exj)
          if (idiag_exjmx/=0) call sum_mn_name(exj(:,1),idiag_exjmx)
          if (idiag_exjmy/=0) call sum_mn_name(exj(:,2),idiag_exjmy)
          if (idiag_exjmz/=0) call sum_mn_name(exj(:,3),idiag_exjmz)
        endif
!
!  Calculate part II of current helicity flux (for imposed field).
!  < curlE x B >|_i  =  < B_{j,i} E_j >
!  Use the full B (with B_ext)
!
        if (idiag_dexbmx/=0 .or. idiag_dexbmy/=0 .or. idiag_dexbmz/=0) then
          call multmv_transp(p%bij,-p%uxb+eta*p%jj,dexb)
          if (idiag_dexbmx/=0) call sum_mn_name(dexb(:,1),idiag_dexbmx)
          if (idiag_dexbmy/=0) call sum_mn_name(dexb(:,2),idiag_dexbmy)
          if (idiag_dexbmz/=0) call sum_mn_name(dexb(:,3),idiag_dexbmz)
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
        if (idiag_uxBrms/=0) call sum_mn_name(p%uxb2,idiag_uxBrms,lsqrt=.true.)
        if (idiag_Bresrms/=0 .or. idiag_Rmrms/=0) then
          call dot2_mn(fres,fres2)
          if (idiag_Bresrms/=0) &
              call sum_mn_name(fres2,idiag_Bresrms,lsqrt=.true.)
          if (idiag_Rmrms/=0) &
              call sum_mn_name(p%uxb2/fres2,idiag_Rmrms,lsqrt=.true.)
        endif
!
!  Calculate <b^2*divu>, which is part of <u.(jxb)>.
!  Note that <u.(jxb)>=1/2*<b^2*divu>+<u.bgradb>.
!
        if (idiag_b2divum/=0) call sum_mn_name(p%b2*p%divu,idiag_b2divum)
!
!  Calculate <u.(jxb)>.
!
        if (idiag_ujxbm/=0) call sum_mn_name(p%ujxb,idiag_ujxbm)
!
!  Calculate <jxb>.B_0/B_0^2.
!
        call dot(B_ext_inv,p%jxb,jxb_dotB0)
        if (idiag_jxbm/=0) call sum_mn_name(jxb_dotB0,idiag_jxbm)
        if (idiag_jxbmx/=0.or.idiag_jxbmy/=0.or.idiag_jxbmz/=0) then
          call cross_mn(p%jj,p%bbb,jxbb)
          if (idiag_jxbmx/=0) call sum_mn_name(jxbb(:,1),idiag_jxbmx)
          if (idiag_jxbmy/=0) call sum_mn_name(jxbb(:,2),idiag_jxbmy)
          if (idiag_jxbmz/=0) call sum_mn_name(jxbb(:,3),idiag_jxbmz)
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
!  alpM22=<b1*b3,2>
!
        if (idiag_b1b32m/=0) then
          b1b32=p%bb(:,1)*p%bij(:,3,2)
          call sum_mn_name(b1b32,idiag_b1b32m)
        endif
!
!  alpM33=<b2*b1,3>
!
        if (idiag_b2b13m/=0) then
          b2b13=p%bb(:,2)*p%bij(:,1,3)
          call sum_mn_name(b2b13,idiag_b2b13m)
        endif
!
      endif ! endif (ldiagnos)
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1.
!
      if (l1davgfirst .or. (ldiagnos .and. ldiagnos_need_zaverages)) then
        if (idiag_bxmx/=0)   call yzsum_mn_name_x(p%bb(:,1),idiag_bxmx)
        if (idiag_bymx/=0)   call yzsum_mn_name_x(p%bb(:,2),idiag_bymx)
        if (idiag_bzmx/=0)   call yzsum_mn_name_x(p%bb(:,3),idiag_bzmx)
        if (idiag_bx2mx/=0)  call yzsum_mn_name_x(p%bb(:,1)**2,idiag_bx2mx)
        if (idiag_by2mx/=0)  call yzsum_mn_name_x(p%bb(:,2)**2,idiag_by2mx)
        if (idiag_bz2mx/=0)  call yzsum_mn_name_x(p%bb(:,3)**2,idiag_bz2mx)
        if (idiag_bxbymx/=0) &
            call yzsum_mn_name_x(p%bbb(:,1)*p%bbb(:,2),idiag_bxbymx)
        if (idiag_bxmy/=0)   call xzsum_mn_name_y(p%bb(:,1),idiag_bxmy)
        if (idiag_bymy/=0)   call xzsum_mn_name_y(p%bb(:,2),idiag_bymy)
        if (idiag_bzmy/=0)   call xzsum_mn_name_y(p%bb(:,3),idiag_bzmy)
        if (idiag_bx2my/=0)  call xzsum_mn_name_y(p%bb(:,1)**2,idiag_bx2my)
        if (idiag_by2my/=0)  call xzsum_mn_name_y(p%bb(:,2)**2,idiag_by2my)
        if (idiag_bz2my/=0)  call xzsum_mn_name_y(p%bb(:,3)**2,idiag_bz2my)
        if (idiag_axmz/=0)   call xysum_mn_name_z(p%aa(:,1),idiag_axmz)
        if (idiag_aymz/=0)   call xysum_mn_name_z(p%aa(:,2),idiag_aymz)
        if (idiag_azmz/=0)   call xysum_mn_name_z(p%aa(:,3),idiag_azmz)
        if (idiag_abuxmz/=0) call xysum_mn_name_z(p%ab*p%uu(:,1),idiag_abuxmz)
        if (idiag_abuymz/=0) call xysum_mn_name_z(p%ab*p%uu(:,2),idiag_abuymz)
        if (idiag_abuzmz/=0) call xysum_mn_name_z(p%ab*p%uu(:,3),idiag_abuzmz)
        if (idiag_uabxmz/=0) call xysum_mn_name_z(p%ua*p%bb(:,1),idiag_uabxmz)
        if (idiag_uabymz/=0) call xysum_mn_name_z(p%ua*p%bb(:,2),idiag_uabymz)
        if (idiag_uabzmz/=0) call xysum_mn_name_z(p%ua*p%bb(:,3),idiag_uabzmz)
        if (idiag_bxmz/=0)   call xysum_mn_name_z(p%bb(:,1),idiag_bxmz)
        if (idiag_bymz/=0)   call xysum_mn_name_z(p%bb(:,2),idiag_bymz)
        if (idiag_bzmz/=0)   call xysum_mn_name_z(p%bb(:,3),idiag_bzmz)
        if (idiag_jxmz/=0)   call xysum_mn_name_z(p%jj(:,1),idiag_jxmz)
        if (idiag_jymz/=0)   call xysum_mn_name_z(p%jj(:,2),idiag_jymz)
        if (idiag_jzmz/=0)   call xysum_mn_name_z(p%jj(:,3),idiag_jzmz)
        if (idiag_Exmz/=0)   call xysum_mn_name_z(p%uxb(:,1),idiag_Exmz)
        if (idiag_Eymz/=0)   call xysum_mn_name_z(p%uxb(:,2),idiag_Eymz)
        if (idiag_Ezmz/=0)   call xysum_mn_name_z(p%uxb(:,3),idiag_Ezmz)
        if (idiag_bx2mz/=0)  call xysum_mn_name_z(p%bb(:,1)**2,idiag_bx2mz)
        if (idiag_by2mz/=0)  call xysum_mn_name_z(p%bb(:,2)**2,idiag_by2mz)
        if (idiag_bz2mz/=0)  call xysum_mn_name_z(p%bb(:,3)**2,idiag_bz2mz)
        if (idiag_bx2rmz/=0) call xysum_mn_name_z(p%bb(:,1)**2*p%rho1,idiag_bx2rmz)
        if (idiag_by2rmz/=0) call xysum_mn_name_z(p%bb(:,2)**2*p%rho1,idiag_by2rmz)
        if (idiag_bz2rmz/=0) call xysum_mn_name_z(p%bb(:,3)**2*p%rho1,idiag_bz2rmz)
        if (idiag_jbmz/=0)   call xysum_mn_name_z(p%jb,idiag_jbmz)
        if (idiag_abmz/=0)   call xysum_mn_name_z(p%ab,idiag_abmz)
        if (idiag_uamz/=0)   call xysum_mn_name_z(p%ua,idiag_uamz)
        if (idiag_etatotalmx/=0) call yzsum_mn_name_x(etatotal,idiag_etatotalmx)
        if (idiag_etatotalmz/=0) call xysum_mn_name_z(etatotal,idiag_etatotalmz)
!
!  Calculate magnetic helicity flux (ExA contribution).
!
        if (idiag_examz1/=0) call xysum_mn_name_z(p%exa(:,1),idiag_examz1)
        if (idiag_examz2/=0) call xysum_mn_name_z(p%exa(:,2),idiag_examz2)
        if (idiag_examz3/=0) call xysum_mn_name_z(p%exa(:,3),idiag_examz3)
!
!  Maxwell stress components.
!
        if (idiag_bxbymy/=0) &
            call xzsum_mn_name_y(p%bbb(:,1)*p%bbb(:,2),idiag_bxbymy)
        if (idiag_bxbzmy/=0) &
            call xzsum_mn_name_y(p%bbb(:,1)*p%bbb(:,3),idiag_bxbzmy)
        if (idiag_bybzmy/=0) &
            call xzsum_mn_name_y(p%bbb(:,2)*p%bbb(:,3),idiag_bybzmy)
        if (idiag_bxbymz/=0) &
            call xysum_mn_name_z(p%bbb(:,1)*p%bbb(:,2),idiag_bxbymz)
        if (idiag_bxbzmz/=0) &
            call xysum_mn_name_z(p%bbb(:,1)*p%bbb(:,3),idiag_bxbzmz)
        if (idiag_bybzmz/=0) &
            call xysum_mn_name_z(p%bbb(:,2)*p%bbb(:,3),idiag_bybzmz)
        if (idiag_jxbrxmx/=0) call yzsum_mn_name_x(p%jxbr(:,1),idiag_jxbrxmx)
        if (idiag_jxbrymx/=0) call yzsum_mn_name_x(p%jxbr(:,2),idiag_jxbrymx)
        if (idiag_jxbrzmx/=0) call yzsum_mn_name_x(p%jxbr(:,3),idiag_jxbrzmx)
        if (idiag_jxbrxmy/=0) call xzsum_mn_name_y(p%jxbr(:,1),idiag_jxbrxmy)
        if (idiag_jxbrymy/=0) call xzsum_mn_name_y(p%jxbr(:,2),idiag_jxbrymy)
        if (idiag_jxbrzmy/=0) call xzsum_mn_name_y(p%jxbr(:,3),idiag_jxbrzmy)
        if (idiag_jxbrxmz/=0) call xysum_mn_name_z(p%jxbr(:,1),idiag_jxbrxmz)
        if (idiag_jxbrymz/=0) call xysum_mn_name_z(p%jxbr(:,2),idiag_jxbrymz)
        if (idiag_jxbrzmz/=0) call xysum_mn_name_z(p%jxbr(:,3),idiag_jxbrzmz)
        if (idiag_b2mz/=0) call xysum_mn_name_z(p%b2,idiag_b2mz)
        if (idiag_mflux_x/=0) &
          call yzintegrate_mn_name_x(p%bb(:,1),idiag_mflux_x)
        if (idiag_mflux_y/=0) &
          call xzintegrate_mn_name_y(p%bb(:,2),idiag_mflux_y)
        if (idiag_mflux_z/=0) &
          call xyintegrate_mn_name_z(p%bb(:,3),idiag_mflux_z)
      endif
!
!  2-D averages.
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        if (idiag_bxmxy/=0)  call zsum_mn_name_xy(p%bb(:,1),idiag_bxmxy)
        if (idiag_bymxy/=0)  call zsum_mn_name_xy(p%bb(:,2),idiag_bymxy)
        if (idiag_bzmxy/=0)  call zsum_mn_name_xy(p%bb(:,3),idiag_bzmxy)
        if (idiag_jxmxy/=0)  call zsum_mn_name_xy(p%jj(:,1),idiag_jxmxy)
        if (idiag_jymxy/=0)  call zsum_mn_name_xy(p%jj(:,2),idiag_jymxy)
        if (idiag_jzmxy/=0)  call zsum_mn_name_xy(p%jj(:,3),idiag_jzmxy)
        if (idiag_axmxz/=0)  call ysum_mn_name_xz(p%aa(:,1),idiag_axmxz)
        if (idiag_aymxz/=0)  call ysum_mn_name_xz(p%aa(:,2),idiag_aymxz)
        if (idiag_azmxz/=0)  call ysum_mn_name_xz(p%aa(:,3),idiag_azmxz)
        if (idiag_bxmxz/=0)  call ysum_mn_name_xz(p%bb(:,1),idiag_bxmxz)
        if (idiag_bymxz/=0)  call ysum_mn_name_xz(p%bb(:,2),idiag_bymxz)
        if (idiag_bzmxz/=0)  call ysum_mn_name_xz(p%bb(:,3),idiag_bzmxz)
        if (idiag_bx2mxz/=0) call ysum_mn_name_xz(p%bb(:,1)**2,idiag_bx2mxz)
        if (idiag_by2mxz/=0) call ysum_mn_name_xz(p%bb(:,2)**2,idiag_by2mxz)
        if (idiag_bz2mxz/=0) call ysum_mn_name_xz(p%bb(:,3)**2,idiag_bz2mxz)
        if (idiag_bx2mxy/=0) call zsum_mn_name_xy(p%bb(:,1)**2,idiag_bx2mxy)
        if (idiag_by2mxy/=0) call zsum_mn_name_xy(p%bb(:,2)**2,idiag_by2mxy)
        if (idiag_bz2mxy/=0) call zsum_mn_name_xy(p%bb(:,3)**2,idiag_bz2mxy)
        if (idiag_bxbymxy/=0) &
            call zsum_mn_name_xy(p%bb(:,1)*p%bb(:,2),idiag_bxbymxy)
        if (idiag_bxbzmxy/=0) &
            call zsum_mn_name_xy(p%bb(:,1)*p%bb(:,3),idiag_bxbzmxy)
        if (idiag_bybzmxy/=0) &
            call zsum_mn_name_xy(p%bb(:,2)*p%bb(:,3),idiag_bybzmxy)
        if (idiag_bxbymxz/=0) &
            call ysum_mn_name_xz(p%bb(:,1)*p%bb(:,2),idiag_bxbymxz)
        if (idiag_bxbzmxz/=0) &
            call ysum_mn_name_xz(p%bb(:,1)*p%bb(:,3),idiag_bxbzmxz)
        if (idiag_bybzmxz/=0) &
            call ysum_mn_name_xz(p%bb(:,2)*p%bb(:,3),idiag_bybzmxz)
        if (idiag_Exmxz/=0) call ysum_mn_name_xz(p%uxb(:,1),idiag_Exmxz)
        if (idiag_Eymxz/=0) call ysum_mn_name_xz(p%uxb(:,2),idiag_Eymxz)
        if (idiag_Ezmxz/=0) call ysum_mn_name_xz(p%uxb(:,3),idiag_Ezmxz)
      else
!
!  idiag_bxmxy and idiag_bymxy also need to be calculated when
!  ldiagnos and idiag_bmx and/or idiag_bmy, so
!
!  We may need to calculate bxmxy without calculating bmx. The following
!  if condition was messing up calculation of bmxy_rms
!
        if (ldiagnos) then
          if (idiag_bxmxy/=0) call zsum_mn_name_xy(p%bb(:,1),idiag_bxmxy)
          if (idiag_bymxy/=0) call zsum_mn_name_xy(p%bb(:,2),idiag_bymxy)
          if (idiag_bzmxy/=0) call zsum_mn_name_xy(p%bb(:,3),idiag_bzmxy)
          if (idiag_jxmxy/=0) call zsum_mn_name_xy(p%jj(:,1),idiag_jxmxy)
          if (idiag_jymxy/=0) call zsum_mn_name_xy(p%jj(:,2),idiag_jymxy)
          if (idiag_jzmxy/=0) call zsum_mn_name_xy(p%jj(:,3),idiag_jzmxy)
        endif
      endif
!
!  Debug output.
!
      if (headtt .and. lfirst .and. ip<=4) then
        call output_pencil(trim(directory)//'/aa.dat',p%aa,3)
        call output_pencil(trim(directory)//'/bb.dat',p%bb,3)
        call output_pencil(trim(directory)//'/jj.dat',p%jj,3)
        call output_pencil(trim(directory)//'/del2A.dat',p%del2a,3)
        call output_pencil(trim(directory)//'/JxBr.dat',p%jxbr,3)
        call output_pencil(trim(directory)//'/JxB.dat',p%jxb,3)
        call output_pencil(trim(directory)//'/df.dat',df(l1:l2,m,n,:),mvar)
      endif
!
      call timing('daa_dt','finished',mnloop=.true.)
!
    endsubroutine daa_dt
!***********************************************************************
    subroutine read_magnetic_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=magnetic_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=magnetic_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_magnetic_init_pars
!***********************************************************************
    subroutine write_magnetic_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=magnetic_init_pars)
!
    endsubroutine write_magnetic_init_pars
!***********************************************************************
    subroutine read_magnetic_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=magnetic_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=magnetic_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_magnetic_run_pars
!***********************************************************************
    subroutine write_magnetic_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=magnetic_run_pars)
!
    endsubroutine write_magnetic_run_pars
!***********************************************************************
    subroutine get_slices_magnetic(f,slices)
!
!  Write slices for animation of Magnetic variables.
!
!  26-jul-06/tony: coded
!
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
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =f(slices%ix,m1:m2    ,n1:n2     ,iax-1+slices%index)
            slices%xz =f(l1:l2    ,slices%iy,n1:n2     ,iax-1+slices%index)
            slices%xy =f(l1:l2    ,m1:m2    ,slices%iz ,iax-1+slices%index)
            slices%xy2=f(l1:l2    ,m1:m2    ,slices%iz2,iax-1+slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
!
!  Magnetic field (derived variable)
!
        case ('bb')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =>bb_yz(:,:,slices%index)
            slices%xz =>bb_xz(:,:,slices%index)
            slices%xy =>bb_xy(:,:,slices%index)
            slices%xy2=>bb_xy2(:,:,slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
!
!  Current density (derived variable)
!
        case ('jj')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =>jj_yz(:,:,slices%index)
            slices%xz =>jj_xz(:,:,slices%index)
            slices%xy =>jj_xy(:,:,slices%index)
            slices%xy2=>jj_xy2(:,:,slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
!
!  Magnetic field squared (derived variable)
!
        case ('b2')
          slices%yz =>b2_yz
          slices%xz =>b2_xz
          slices%xy =>b2_xy
          slices%xy2=>b2_xy2
          slices%ready=.true.
!
!  Current density (derived variable)
!
        case ('jb')
          slices%yz =>jb_yz
          slices%xz =>jb_xz
          slices%xy =>jb_xy
          slices%xy2=>jb_xy2
          slices%ready=.true.
!
!  Plasma beta
!
       case ('beta')
          slices%yz =>jb_yz
          slices%xz =>jb_xz
          slices%xy =>jb_xy
          slices%xy2=>jb_xy2
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_magnetic
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
      integer :: iname,inamex,inamey,inamez,ixy,ixz,iname_half
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of RELOAD.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_ab_int=0; idiag_jb_int=0; idiag_b2tm=0; idiag_bjtm=0; idiag_jbtm=0
        idiag_b2uzm=0; idiag_b2ruzm=0; idiag_ubbzm=0; idiag_b1m=0; idiag_b2m=0
        idiag_bm2=0; idiag_j2m=0; idiag_jm2=0
        idiag_abm=0; idiag_abrms=0; idiag_abmh=0
        idiag_abumx=0; idiag_abumy=0; idiag_abumz=0
        idiag_abmn=0; idiag_abms=0; idiag_jbmh=0; idiag_jbmn=0; idiag_jbms=0
        idiag_ajm=0; idiag_cosubm=0; idiag_jbm=0
        idiag_uam=0; idiag_ubm=0; idiag_ujm=0
        idiag_uxbxm=0; idiag_uybym=0; idiag_uzbzm=0
        idiag_fbm=0; idiag_fxbxm=0; idiag_epsM=0; idiag_epsM_LES=0
        idiag_epsAD=0; idiag_bxpt=0; idiag_bypt=0; idiag_bzpt=0; idiag_Expt=0
        idiag_Eypt=0; idiag_Ezpt=0; idiag_aybym2=0; idiag_exaym2=0
        idiag_exjm2=0; idiag_brms=0; idiag_bmax=0; idiag_jrms=0; idiag_jmax=0
        idiag_vArms=0; idiag_emag=0; idiag_bxmin=0; idiag_bymin=0; idiag_bzmin=0
        idiag_bxmax=0; idiag_bymax=0; idiag_bzmax=0; idiag_vAmax=0; idiag_dtb=0
        idiag_arms=0; idiag_amax=0; idiag_Qsm=0; idiag_Qpm=0; idiag_qem=0; idiag_beta1m=0
        idiag_beta1max=0; idiag_bxm=0; idiag_bym=0; idiag_bzm=0; idiag_axm=0
        idiag_aym=0; idiag_azm=0; idiag_bx2m=0; idiag_by2m=0; idiag_bz2m=0
        idiag_bxbymy=0; idiag_bxbzmy=0; idiag_bybzmy=0; idiag_bxbymz=0
        idiag_bxbzmz=0; idiag_bybzmz=0; idiag_b2mz=0
        idiag_jbmz=0; idiag_abmz=0; idiag_uamz=0
        idiag_bxbym=0; idiag_bxbzm=0; idiag_bybzm=0; idiag_djuidjbim=0
        idiag_axmz=0; idiag_aymz=0; idiag_azmz=0; idiag_bxmz=0; idiag_bymz=0
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
        idiag_jxmxy=0; idiag_jymxy=0; idiag_jzmxy=0
        idiag_bx2mxy=0; idiag_by2mxy=0; idiag_bz2mxy=0; idiag_bxbymxy=0
        idiag_bxbzmxy=0; idiag_bybzmxy=0; idiag_bxbymxz=0; idiag_bxbzmxz=0
        idiag_bybzmxz=0; idiag_bxmxz=0; idiag_bymxz=0; idiag_bzmxz=0
        idiag_axmxz=0; idiag_aymxz=0; idiag_azmxz=0; idiag_Exmxz=0
        idiag_Eymxz=0; idiag_Ezmxz=0; idiag_jxbm=0; idiag_uxbm=0; idiag_oxuxbm=0
        idiag_jxbxbm=0; idiag_gpxbm=0; idiag_uxDxuxbm=0
        idiag_uxbmx=0; idiag_uxbmy=0; idiag_uxbmz=0
        idiag_jxbmx=0; idiag_jxbmy=0; idiag_jxbmz=0
        idiag_uxbcmx=0; idiag_uxbsmx=0
        idiag_uxbcmy=0; idiag_uxbsmy=0; idiag_examz1=0; idiag_examz2=0
        idiag_examz3=0; idiag_examx=0; idiag_examy=0; idiag_examz=0
        idiag_exjmx=0; idiag_exjmy=0; idiag_exjmz=0; idiag_dexbmx=0
        idiag_dexbmy=0; idiag_dexbmz=0; idiag_phibmx=0; idiag_phibmy=0
        idiag_phibmz=0; idiag_uxjm=0; idiag_ujxbm=0; idiag_b2divum=0
        idiag_b3b21m=0; idiag_b1b32m=0; idiag_b2b13m=0
        idiag_udotxbm=0; idiag_uxbdotm=0
        idiag_jxbrxm=0; idiag_jxbrym=0; idiag_jxbrzm=0
        idiag_jxbr2m=0; idiag_jxbrxmx=0; idiag_jxbrymx=0; idiag_jxbrzmx=0
        idiag_jxbrxmy=0; idiag_jxbrymy=0; idiag_jxbrzmy=0; idiag_jxbrxmz=0
        idiag_jxbrymz=0; idiag_jxbrzmz=0
        idiag_dteta=0; idiag_uxBrms=0; idiag_Bresrms=0
        idiag_Rmrms=0; idiag_jfm=0; idiag_va2m=0
        idiag_bxmx=0; idiag_bymx=0; idiag_bzmx=0; idiag_bxmy=0
        idiag_bymy=0; idiag_bzmy=0; idiag_bx2my=0; idiag_by2my=0; idiag_bz2my=0
        idiag_mflux_x=0; idiag_mflux_y=0; idiag_mflux_z=0; idiag_bmxy_rms=0
        idiag_brmsh=0; idiag_brmsn=0
        idiag_brmss=0; idiag_etatotalmx=0; idiag_etatotalmz=0
        idiag_etavamax=0; idiag_etajmax=0; idiag_etaj2max=0; idiag_etajrhomax=0
      endif
!
!  Check for those quantities that we want to evaluate online.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ab_int',idiag_ab_int)
        call parse_name(iname,cname(iname),cform(iname),'jb_int',idiag_jb_int)
        call parse_name(iname,cname(iname),cform(iname),'dteta',idiag_dteta)
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
        call parse_name(iname,cname(iname),cform(iname),'abumx',idiag_abumx)
        call parse_name(iname,cname(iname),cform(iname),'abumy',idiag_abumy)
        call parse_name(iname,cname(iname),cform(iname),'abumz',idiag_abumz)
        call parse_name(iname,cname(iname),cform(iname),'ajm',idiag_ajm)
        call parse_name(iname,cname(iname),cform(iname),'jbm',idiag_jbm)
        call parse_name(iname,cname(iname),cform(iname),'jbmn',idiag_jbmn)
        call parse_name(iname,cname(iname),cform(iname),'jbms',idiag_jbms)
        call parse_name(iname,cname(iname),cform(iname),'ubm',idiag_ubm)
        call parse_name(iname,cname(iname),cform(iname),'uxbxm',idiag_uxbxm)
        call parse_name(iname,cname(iname),cform(iname),'uybym',idiag_uybym)
        call parse_name(iname,cname(iname),cform(iname),'uzbzm',idiag_uzbzm)
        call parse_name(iname,cname(iname),cform(iname),'cosubm',idiag_cosubm)
        call parse_name(iname,cname(iname),cform(iname),'uam',idiag_uam)
        call parse_name(iname,cname(iname),cform(iname),'ujm',idiag_ujm)
        call parse_name(iname,cname(iname),cform(iname),'fbm',idiag_fbm)
        call parse_name(iname,cname(iname),cform(iname),'fxbxm',idiag_fxbxm)
        call parse_name(iname,cname(iname),cform(iname),'b2ruzm',idiag_b2ruzm)
        call parse_name(iname,cname(iname),cform(iname),'b2uzm',idiag_b2uzm)
        call parse_name(iname,cname(iname),cform(iname),'ubbzm',idiag_ubbzm)
        call parse_name(iname,cname(iname),cform(iname),'b1m',idiag_b1m)
        call parse_name(iname,cname(iname),cform(iname),'b2m',idiag_b2m)
        call parse_name(iname,cname(iname),cform(iname),'bm2',idiag_bm2)
        call parse_name(iname,cname(iname),cform(iname),'j2m',idiag_j2m)
        call parse_name(iname,cname(iname),cform(iname),'jm2',idiag_jm2)
        call parse_name(iname,cname(iname),cform(iname),'epsM',idiag_epsM)
        call parse_name(iname,cname(iname),cform(iname),&
            'epsM_LES',idiag_epsM_LES)
        call parse_name(iname,cname(iname),cform(iname),'epsAD',idiag_epsAD)
        call parse_name(iname,cname(iname),cform(iname),'emag',idiag_emag)
        call parse_name(iname,cname(iname),cform(iname),'brms',idiag_brms)
        call parse_name(iname,cname(iname),cform(iname),'brmsn',idiag_brmsn)
        call parse_name(iname,cname(iname),cform(iname),'brmss',idiag_brmss)
        call parse_name(iname,cname(iname),cform(iname),'bmax',idiag_bmax)
        call parse_name(iname,cname(iname),cform(iname),'bxmin',idiag_bxmin)
        call parse_name(iname,cname(iname),cform(iname),'bymin',idiag_bymin)
        call parse_name(iname,cname(iname),cform(iname),'bzmin',idiag_bzmin)
        call parse_name(iname,cname(iname),cform(iname),'bxmax',idiag_bxmax)
        call parse_name(iname,cname(iname),cform(iname),'bymax',idiag_bymax)
        call parse_name(iname,cname(iname),cform(iname),'bzmax',idiag_bzmax)
        call parse_name(iname,cname(iname),cform(iname),'jrms',idiag_jrms)
        call parse_name(iname,cname(iname),cform(iname),'jmax',idiag_jmax)
        call parse_name(iname,cname(iname),cform(iname),'axm',idiag_axm)
        call parse_name(iname,cname(iname),cform(iname),'aym',idiag_aym)
        call parse_name(iname,cname(iname),cform(iname),'azm',idiag_azm)
        call parse_name(iname,cname(iname),cform(iname),'arms',idiag_arms)
        call parse_name(iname,cname(iname),cform(iname),'amax',idiag_amax)
        call parse_name(iname,cname(iname),cform(iname),'vArms',idiag_vArms)
        call parse_name(iname,cname(iname),cform(iname),'vAmax',idiag_vAmax)
        call parse_name(iname,cname(iname),cform(iname),'vA2m',idiag_vA2m)
        call parse_name(iname,cname(iname),cform(iname),'Qsm',idiag_Qsm)
        call parse_name(iname,cname(iname),cform(iname),'Qpm',idiag_Qpm)
        call parse_name(iname,cname(iname),cform(iname),'qem',idiag_qem)
        call parse_name(iname,cname(iname),cform(iname),&
            'beta1m',idiag_beta1m)
        call parse_name(iname,cname(iname),cform(iname),&
            'beta1max',idiag_beta1max)
        call parse_name(iname,cname(iname),cform(iname),'dtb',idiag_dtb)
        call parse_name(iname,cname(iname),cform(iname),'bxm',idiag_bxm)
        call parse_name(iname,cname(iname),cform(iname),'bym',idiag_bym)
        call parse_name(iname,cname(iname),cform(iname),'bzm',idiag_bzm)
        call parse_name(iname,cname(iname),cform(iname),'bx2m',idiag_bx2m)
        call parse_name(iname,cname(iname),cform(iname),'by2m',idiag_by2m)
        call parse_name(iname,cname(iname),cform(iname),'bz2m',idiag_bz2m)
        call parse_name(iname,cname(iname),cform(iname),'bxbym',idiag_bxbym)
        call parse_name(iname,cname(iname),cform(iname),'bxbzm',idiag_bxbzm)
        call parse_name(iname,cname(iname),cform(iname),'bybzm',idiag_bybzm)
        call parse_name(iname,cname(iname),cform(iname),&
            'djuidjbim',idiag_djuidjbim)
        call parse_name(iname,cname(iname),cform(iname),'jxbrxm',idiag_jxbrxm)
        call parse_name(iname,cname(iname),cform(iname),'jxbrym',idiag_jxbrym)
        call parse_name(iname,cname(iname),cform(iname),'jxbrzm',idiag_jxbrzm)
        call parse_name(iname,cname(iname),cform(iname),'jxbr2m',idiag_jxbr2m)
        call parse_name(iname,cname(iname),cform(iname),'jxbm',idiag_jxbm)
        call parse_name(iname,cname(iname),cform(iname),'uxbm',idiag_uxbm)
        call parse_name(iname,cname(iname),cform(iname),'uxbmx',idiag_uxbmx)
        call parse_name(iname,cname(iname),cform(iname),'uxbmy',idiag_uxbmy)
        call parse_name(iname,cname(iname),cform(iname),'uxbmz',idiag_uxbmz)
        call parse_name(iname,cname(iname),cform(iname),'jxbmx',idiag_jxbmx)
        call parse_name(iname,cname(iname),cform(iname),'jxbmy',idiag_jxbmy)
        call parse_name(iname,cname(iname),cform(iname),'jxbmz',idiag_jxbmz)
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
        call parse_name(iname,cname(iname),cform(iname),'ujxbm',idiag_ujxbm)
        call parse_name(iname,cname(iname),cform(iname),'b2divum',idiag_b2divum)
        call parse_name(iname,cname(iname),cform(iname),'jxbxbm',idiag_jxbxbm)
        call parse_name(iname,cname(iname),cform(iname),'oxuxbm',idiag_oxuxbm)
        call parse_name(iname,cname(iname),cform(iname),'gpxbm',idiag_gpxbm)
        call parse_name(iname,cname(iname),cform(iname),&
            'uxDxuxbm',idiag_uxDxuxbm)
        call parse_name(iname,cname(iname),cform(iname),'b3b21m',idiag_b3b21m)
        call parse_name(iname,cname(iname),cform(iname),'b1b32m',idiag_b1b32m)
        call parse_name(iname,cname(iname),cform(iname),'b2b13m',idiag_b2b13m)
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
        call parse_name(iname,cname(iname),cform(iname),'Expt',idiag_Expt)
        call parse_name(iname,cname(iname),cform(iname),'Eypt',idiag_Eypt)
        call parse_name(iname,cname(iname),cform(iname),'Ezpt',idiag_Ezpt)
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
        call parse_name(iname,cname(iname),cform(iname),'cosjbm',idiag_cosjbm)
        call parse_name(iname,cname(iname),cform(iname),'jparallelm',idiag_jparallelm)
        call parse_name(iname,cname(iname),cform(iname),'jperpm',idiag_jperpm)
      enddo
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
!  Check for those quantities for which we want xy-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bxmx',idiag_bxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bymx',idiag_bymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bzmx',idiag_bzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'bx2mx',idiag_bx2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'by2mx',idiag_by2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'bz2mx',idiag_bz2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
            'bxbymx',idiag_bxbymx)
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
!  Check for those quantities for which we want yz-averages.
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
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'bxbymz',idiag_bxbymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'bxbzmz',idiag_bxbzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'bybzmz',idiag_bybzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'b2mz',idiag_b2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'jxbrxmz',idiag_jxbrxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'jxbrymz',idiag_jxbrymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'jxbrzmz',idiag_jxbrzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'mflux_z',idiag_mflux_z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'jbmz',idiag_jbmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'abmz',idiag_abmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uamz',idiag_uamz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'examz1',idiag_examz1)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'examz2',idiag_examz2)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'examz3',idiag_examz3)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
            'etatotalmz',idiag_etatotalmz)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do ixz=1,nnamexz
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'axmxz',idiag_axmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'aymxz',idiag_aymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'azmxz',idiag_azmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bxmxz',idiag_bxmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bymxz',idiag_bymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bzmxz',idiag_bzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bx2mxz',idiag_bx2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'by2mxz',idiag_by2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bz2mxz',idiag_bz2mxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bxbymxz',idiag_bxbymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bxbzmxz',idiag_bxbzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bybzmxz',idiag_bybzmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'Exmxz',idiag_Exmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'Eymxz',idiag_Eymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'Ezmxz',idiag_Ezmxz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do ixy=1,nnamexy
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
      enddo
!
!  Write column, idiag_XYZ, where our variable XYZ is stored.
!
      if (lwr) then
        write(3,*) 'i_ab_int=',idiag_ab_int
        write(3,*) 'i_jb_int=',idiag_jb_int
        write(3,*) 'i_dteta=',idiag_dteta
        write(3,*) 'i_aybym2=',idiag_aybym2
        write(3,*) 'i_exaym2=',idiag_exaym2
        write(3,*) 'i_exabot=',idiag_exabot
        write(3,*) 'i_exatop=',idiag_exatop
        write(3,*) 'i_exjm2=',idiag_exjm2
        write(3,*) 'i_b2tm=',idiag_b2tm
        write(3,*) 'i_bjtm=',idiag_bjtm
        write(3,*) 'i_jbtm=',idiag_jbtm
        write(3,*) 'i_abm=',idiag_abm
        write(3,*) 'i_abumx=',idiag_abumx
        write(3,*) 'i_abumy=',idiag_abumy
        write(3,*) 'i_abumz=',idiag_abumz
        write(3,*) 'i_abmh=',idiag_abmh
        write(3,*) 'i_abmn=',idiag_abmn
        write(3,*) 'i_abms=',idiag_abms
        write(3,*) 'i_abrms=',idiag_abrms
        write(3,*) 'i_ajm=',idiag_ajm
        write(3,*) 'i_brmsn=',idiag_brmsn
        write(3,*) 'i_brmss=',idiag_brmss
        write(3,*) 'i_brmsh=',idiag_brmsh
        write(3,*) 'i_jbm=',idiag_jbm
        write(3,*) 'i_jbmh=',idiag_abmh
        write(3,*) 'i_jbmn=',idiag_abmn
        write(3,*) 'i_jbms=',idiag_abms
        write(3,*) 'i_ubm=',idiag_ubm
        write(3,*) 'i_uxbxm=',idiag_uxbxm
        write(3,*) 'i_uybym=',idiag_uybym
        write(3,*) 'i_uzbzm=',idiag_uzbzm
        write(3,*) 'i_cosubm=',idiag_cosubm
        write(3,*) 'i_uam=',idiag_uam
        write(3,*) 'i_ujm=',idiag_ujm
        write(3,*) 'i_fbm=',idiag_fbm
        write(3,*) 'i_fxbxm=',idiag_fxbxm
        write(3,*) 'i_b2ruzm=',idiag_b2ruzm
        write(3,*) 'i_b2uzm=',idiag_b2uzm
        write(3,*) 'i_ubbzm=',idiag_ubbzm
        write(3,*) 'i_b1m=',idiag_b1m
        write(3,*) 'i_b2m=',idiag_b2m
        write(3,*) 'i_bm2=',idiag_bm2
        write(3,*) 'i_j2m=',idiag_j2m
        write(3,*) 'i_jm2=',idiag_jm2
        write(3,*) 'i_epsM=',idiag_epsM
        write(3,*) 'i_epsM_LES=',idiag_epsM_LES
        write(3,*) 'i_emag=',idiag_emag
        write(3,*) 'i_brms=',idiag_brms
        write(3,*) 'i_bmax=',idiag_bmax
        write(3,*) 'i_jrms=',idiag_jrms
        write(3,*) 'i_jmax=',idiag_jmax
        write(3,*) 'i_axm=',idiag_axm
        write(3,*) 'i_aym=',idiag_aym
        write(3,*) 'i_azm=',idiag_azm
        write(3,*) 'i_arms=',idiag_arms
        write(3,*) 'i_amax=',idiag_amax
        write(3,*) 'i_vArms=',idiag_vArms
        write(3,*) 'i_vAmax=',idiag_vAmax
        write(3,*) 'i_vA2m=',idiag_vA2m
        write(3,*) 'i_Qsm=',idiag_Qsm
        write(3,*) 'i_Qpm=',idiag_Qpm
        write(3,*) 'i_qem=',idiag_qem
        write(3,*) 'i_beta1m=',idiag_beta1m
        write(3,*) 'i_beta1max=',idiag_beta1max
        write(3,*) 'i_dtb=',idiag_dtb
        write(3,*) 'i_bxm=',idiag_bxm
        write(3,*) 'i_bym=',idiag_bym
        write(3,*) 'i_bzm=',idiag_bzm
        write(3,*) 'i_bx2m=',idiag_bx2m
        write(3,*) 'i_by2m=',idiag_by2m
        write(3,*) 'i_bz2m=',idiag_bz2m
        write(3,*) 'i_bxbym=',idiag_bxbym
        write(3,*) 'i_bxbzm=',idiag_bxbzm
        write(3,*) 'i_bybzm=',idiag_bybzm
        write(3,*) 'i_djuidjbim=',idiag_djuidjbim
        write(3,*) 'i_uxbm=',idiag_uxbm
        write(3,*) 'i_jxbm=',idiag_jxbm
        write(3,*) 'i_uxbmx=',idiag_uxbmx
        write(3,*) 'i_uxbmy=',idiag_uxbmy
        write(3,*) 'i_uxbmz=',idiag_uxbmz
        write(3,*) 'i_jxbmx=',idiag_jxbmx
        write(3,*) 'i_jxbmy=',idiag_jxbmy
        write(3,*) 'i_jxbmz=',idiag_jxbmz
        write(3,*) 'i_examx=',idiag_examx
        write(3,*) 'i_examy=',idiag_examy
        write(3,*) 'i_examz=',idiag_examz
        write(3,*) 'i_exjmx=',idiag_exjmx
        write(3,*) 'i_exjmy=',idiag_exjmy
        write(3,*) 'i_exjmz=',idiag_exjmz
        write(3,*) 'i_dexbmx=',idiag_dexbmx
        write(3,*) 'i_dexbmy=',idiag_dexbmy
        write(3,*) 'i_dexbmz=',idiag_dexbmz
        write(3,*) 'i_phibmx=',idiag_phibmx
        write(3,*) 'i_phibmy=',idiag_phibmy
        write(3,*) 'i_phibmz=',idiag_phibmz
        write(3,*) 'i_uxjm=',idiag_uxjm
        write(3,*) 'i_ujxbm=',idiag_ujxbm
        write(3,*) 'i_b2divum=',idiag_b2divum
        write(3,*) 'i_oxuxbm=',idiag_oxuxbm
        write(3,*) 'i_jxbxbm=',idiag_jxbxbm
        write(3,*) 'i_gpxbm=',idiag_gpxbm
        write(3,*) 'i_uxDxuxbm=',idiag_uxDxuxbm
        write(3,*) 'i_b3b21m=',idiag_b3b21m
        write(3,*) 'i_b1b32m=',idiag_b1b32m
        write(3,*) 'i_b2b13m=',idiag_b2b13m
        write(3,*) 'i_udotxbm=',idiag_udotxbm
        write(3,*) 'i_uxbdotm=',idiag_uxbdotm
        write(3,*) 'i_axmz=',idiag_axmz
        write(3,*) 'i_aymz=',idiag_aymz
        write(3,*) 'i_azmz=',idiag_azmz
        write(3,*) 'i_abuxmz=',idiag_abuxmz
        write(3,*) 'i_abuymz=',idiag_abuymz
        write(3,*) 'i_abuzmz=',idiag_abuzmz
        write(3,*) 'i_uabxmz=',idiag_uabxmz
        write(3,*) 'i_uabymz=',idiag_uabymz
        write(3,*) 'i_uabzmz=',idiag_uabzmz
        write(3,*) 'i_bxmz=',idiag_bxmz
        write(3,*) 'i_bymz=',idiag_bymz
        write(3,*) 'i_bzmz=',idiag_bzmz
        write(3,*) 'i_jxmz=',idiag_jxmz
        write(3,*) 'i_jymz=',idiag_jymz
        write(3,*) 'i_jzmz=',idiag_jzmz
        write(3,*) 'i_Exmz=',idiag_Exmz
        write(3,*) 'i_Eymz=',idiag_Eymz
        write(3,*) 'i_Ezmz=',idiag_Ezmz
        write(3,*) 'i_bx2mz=',idiag_bx2mz
        write(3,*) 'i_by2mz=',idiag_by2mz
        write(3,*) 'i_bz2mz=',idiag_bz2mz
        write(3,*) 'i_bx2rmz=',idiag_bx2rmz
        write(3,*) 'i_by2rmz=',idiag_by2rmz
        write(3,*) 'i_bz2rmz=',idiag_bz2rmz
        write(3,*) 'i_bxbymz=',idiag_bxbymz
        write(3,*) 'i_b2mz=',idiag_b2mz
        write(3,*) 'i_bmx=',idiag_bmx
        write(3,*) 'i_bmy=',idiag_bmy
        write(3,*) 'i_bmz=',idiag_bmz
        write(3,*) 'i_bmzS2=',idiag_bmzS2
        write(3,*) 'i_bmzA2=',idiag_bmzA2
        write(3,*) 'i_jmx=',idiag_jmx
        write(3,*) 'i_jmy=',idiag_jmy
        write(3,*) 'i_jmz=',idiag_jmz
        write(3,*) 'i_bmzph=',idiag_bmzph
        write(3,*) 'i_bmzphe=',idiag_bmzphe
        write(3,*) 'i_bcosphz=',idiag_bcosphz
        write(3,*) 'i_bsinphz=',idiag_bsinphz
        write(3,*) 'i_emxamz3=',idiag_emxamz3
        write(3,*) 'i_embmz=',idiag_embmz
        write(3,*) 'i_ambmz=',idiag_ambmz
        write(3,*) 'i_ambmzh=',idiag_ambmzh
        write(3,*) 'i_ambmzn=',idiag_ambmzn
        write(3,*) 'i_ambmzs=',idiag_ambmzs
        write(3,*) 'i_jmbmz=',idiag_jmbmz
        write(3,*) 'i_kmz=',idiag_kmz
        write(3,*) 'i_kx_aa=',idiag_kx_aa
        write(3,*) 'i_bxpt=',idiag_bxpt
        write(3,*) 'i_bypt=',idiag_bypt
        write(3,*) 'i_bzpt=',idiag_bzpt
        write(3,*) 'i_Expt=',idiag_Expt
        write(3,*) 'i_Eypt=',idiag_Eypt
        write(3,*) 'i_Ezpt=',idiag_Ezpt
        write(3,*) 'i_bxmxy=',idiag_bxmxy
        write(3,*) 'i_bymxy=',idiag_bymxy
        write(3,*) 'i_bzmxy=',idiag_bzmxy
        write(3,*) 'i_jxmxy=',idiag_jxmxy
        write(3,*) 'i_jymxy=',idiag_jymxy
        write(3,*) 'i_jzmxy=',idiag_jzmxy
        write(3,*) 'i_axmxz=',idiag_axmxz
        write(3,*) 'i_aymxz=',idiag_aymxz
        write(3,*) 'i_azmxz=',idiag_azmxz
        write(3,*) 'i_bxmxz=',idiag_bxmxz
        write(3,*) 'i_bymxz=',idiag_bymxz
        write(3,*) 'i_bzmxz=',idiag_bzmxz
        write(3,*) 'i_bx2mxz=',idiag_bx2mxz
        write(3,*) 'i_by2mxz=',idiag_by2mxz
        write(3,*) 'i_bz2mxz=',idiag_bz2mxz
        write(3,*) 'i_Exmxz=',idiag_Exmxz
        write(3,*) 'i_Eymxz=',idiag_Eymxz
        write(3,*) 'i_Ezmxz=',idiag_Ezmxz
        write(3,*) 'i_jbmz=',idiag_jbmz
        write(3,*) 'i_abmz=',idiag_abmz
        write(3,*) 'i_uamz=',idiag_uamz
        write(3,*) 'i_examz1=',idiag_examz1
        write(3,*) 'i_examz2=',idiag_examz2
        write(3,*) 'i_examz3=',idiag_examz3
        write(3,*) 'i_bx2mxy=',idiag_bx2mxy
        write(3,*) 'i_by2mxy=',idiag_by2mxy
        write(3,*) 'i_bz2mxy=',idiag_bz2mxy
        write(3,*) 'i_bxbymxy=',idiag_bxbymxy
        write(3,*) 'i_bxbzmxy=',idiag_bxbzmxy
        write(3,*) 'i_bybzmxy=',idiag_bybzmxy
        write(3,*) 'i_bxbymxz=',idiag_bxbymxz
        write(3,*) 'i_bxbzmxz=',idiag_bxbzmxz
        write(3,*) 'i_bybzmxz=',idiag_bybzmxz
        write(3,*) 'i_uxBrms=',idiag_uxBrms
        write(3,*) 'i_Bresrms=',idiag_Bresrms
        write(3,*) 'i_Rmrms=',idiag_Rmrms
        write(3,*) 'i_jfm=',idiag_jfm
        write(3,*) 'i_bxmx=',idiag_bxmx
        write(3,*) 'i_bxmy=',idiag_bxmy
        write(3,*) 'i_bymy=',idiag_bymy
        write(3,*) 'i_bzmy=',idiag_bzmy
        write(3,*) 'i_bx2my=',idiag_bx2my
        write(3,*) 'i_by2my=',idiag_by2my
        write(3,*) 'i_bz2my=',idiag_bz2my
        write(3,*) 'i_mflux_x=',idiag_mflux_x
        write(3,*) 'i_mflux_y=',idiag_mflux_y
        write(3,*) 'i_mflux_z=',idiag_mflux_z
        write(3,*) 'nname=',nname
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamexz=',nnamexz
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'iaa=',iaa
        write(3,*) 'iax=',iax
        write(3,*) 'iay=',iay
        write(3,*) 'iaz=',iaz
        write(3,*) 'idiag_bmxy_rms=',idiag_bmxy_rms
        write(3,*) 'idiag_cosjbm=',idiag_cosjbm
        write(3,*) 'i_etavamax=', idiag_etavamax
        write(3,*) 'i_etajmax=', idiag_etajmax
        write(3,*) 'i_etaj2max=', idiag_etaj2max
        write(3,*) 'i_etajrhomax=', idiag_etajrhomax
        write(3,*) 'idiag_jparallel=',idiag_jparallelm
        write(3,*) 'idiag_jperp=',idiag_jperpm
      endif
!
    endsubroutine rprint_magnetic
!***********************************************************************
endmodule Magnetic
