! $Id$
!
!  Module containing global parameters (constants).
!
module Cparam
!
  implicit none
!
  integer, parameter :: ikind8=selected_int_kind(14)  ! 8-byte integer kind
  integer, parameter :: ikind4=selected_int_kind(9)   ! 4-byte integer kind
  integer, parameter :: ikind1=selected_int_kind(2)   ! 1-byte integer kind
  integer, parameter :: rkind8=selected_real_kind(12) ! 8-byte real kind
  integer, parameter :: rkind4=selected_real_kind(6)  ! 4-byte real kind
  integer, parameter :: rkind16 = selected_real_kind(33, 4931) ! 16-byte real kind
  !integer, parameter :: rkind16 = rkind8
!
  include 'cparam.local'
!
!
  integer, parameter :: nx=nxgrid/nprocx,ny=nygrid/nprocy,nz=nzgrid/nprocz,nyz=ny*nz
  integer, parameter :: nxygrid=nxgrid*nygrid,nxzgrid=nxgrid*nzgrid,nyzgrid=nygrid*nzgrid
  integer, parameter :: nprocxy=nprocx*nprocy
  integer, parameter :: nprocyz=nprocy*nprocz
  integer, parameter :: nprocxz=nprocx*nprocz
  integer, parameter :: n_forcing_cont_max=2
  integer, parameter :: ndustspec0=8
  character, dimension(3), parameter :: coornames=(/'x','y','z'/)
  character(LEN=2), dimension(12), parameter :: compnames=(/'x ','y ','z ','xx','xy','xz','yx','yy','yz','zx','zy','zz'/)
  integer, dimension(6),parameter :: compinds_6=(/1,2,3,5,6,9/)
  logical, dimension(3), parameter :: lactive_dimension = (/ nxgrid > 1, nygrid > 1, nzgrid > 1 /)
  integer, parameter :: dimensionality=min(nxgrid-1,1)+min(nygrid-1,1)+min(nzgrid-1,1)
  integer, dimension(3), parameter :: grid_dims=(/nx,ny,nz/)
!
  include 'cparam.inc'
  logical, parameter :: lenergy=lentropy.or.ltemperature.or.lthermal_energy
!
  integer, parameter :: penc_name_len=16
!
  include 'cparam_pencils.inc'
!
!  Derived and fixed parameters.
!
! BEGIN CHANGE FOR DYNAMICAL ALLOCATION
  integer, parameter :: mfarray=mvar+maux+mglobal+mscratch
  integer, parameter :: mcom=mvar+maux_com
  integer, parameter :: mparray=mpvar+mpaux
  integer, parameter :: mpcom=mpvar+mpaux
  integer, parameter :: mqarray=mqvar+mqaux
! END CHANGE FOR DYNAMICAL ALLOCATION
!
  integer(KIND=ikind8), parameter :: nw=nx*ny*nz
!
!!!  integer, parameter :: mx=nx+2*nghost,l1=1+nghost,l2=mx-nghost
!!!  integer, parameter :: my=ny+2*nghost,m1=1+nghost,m2=my-nghost
!!!  integer, parameter :: mz=nz+2*nghost,n1=1+nghost,n2=mz-nghost
  integer, parameter :: mx=nx+2*nghost,l1=1+nghost
  integer, parameter :: my=ny+2*nghost,m1=1+nghost
  integer, parameter :: mz=nz+2*nghost,n1=1+nghost
  integer, parameter :: mxgrid=nxgrid+2*nghost
  integer, parameter :: mygrid=nygrid+2*nghost
  integer, parameter :: mzgrid=nzgrid+2*nghost
  integer, parameter :: mw=mx*my*mz
  integer(KIND=ikind8), parameter :: nwgrid=int(nxgrid,kind=ikind8)* &
                                            int(nygrid,kind=ikind8)* &
                                            int(nzgrid,kind=ikind8)
!
!!!  integer, parameter :: l1i=l1+nghost-1,l2i=l2-nghost+1
!!!  integer, parameter :: m1i=m1+nghost-1,m2i=m2-nghost+1
!!!  integer, parameter :: n1i=n1+nghost-1,n2i=n2-nghost+1
  integer, parameter :: l1i=l1+nghost-1
  integer, parameter :: m1i=m1+nghost-1
  integer, parameter :: n1i=n1+nghost-1


  integer, parameter :: nrcyl=nxgrid/2
  integer, parameter :: nrcylrun=max(nx/20,1)
!
!  Number of bins for Pulsar Timing Array
!
  integer, parameter :: nbin_angular=19*2
!
!  Array dimension for reduce operation (maxima and sums).
!  Use here symbol mreduce, use nreduce in call.
!
  integer, parameter :: mreduce=6
!
!  Number of slots in initlnrho etc.
!
  integer, parameter :: ninit=5
!
!  Name:          Maximum string length of a:
!  --------------------------------------------
!  fnlen          file name
!  intlen         integer (64 bit plus sign)
!  bclen          string for boundary condition
!  labellen       label (eg. initss, initaa)
!  linelen        line to be read in
!  datelen        date-and-time string
!  max_col_width  diagnostic column
!  nscbc_len      ?
!
  integer, parameter :: fnlen=135,intlen=21,bclen=3,labellen=40,linelen=256
  integer, parameter :: datelen=30,max_col_width=30,nscbc_len=24,fmtlen=30
!
!  Significant length of random number generator state.
!  Different compilers have different lengths:
!    NAG: 1, Compaq: 2, Intel: 47, SGI: 64, NEC: 256
!
  integer, parameter :: mseed=256
!
!  Predefine maximum possible numbers.
!
  integer(KIND=ikind4), parameter :: int_sgl=0
  integer, parameter :: max_int=huge(int_sgl)
  real, parameter :: huge_real=huge(0.)
  real(KIND=rkind8), parameter :: zero_double=0., huge_double=huge(zero_double)
  real, parameter :: max_real=huge_real/10.    ! division necessary as INTEL compiler considers
                                               ! huge(0.) illegal when reading it from a namelist
!
!  Tiny and huge numbers.
!
  real, parameter :: one_real=1.0
  real, parameter :: epsi=5*epsilon(one_real),tini=5*tiny(one_real)
  real, parameter :: huge1=0.2*huge_real
!
!  A marker value that is highly unlikely ("impossible") to ever occur
!  during a meaningful run: use a very large number.
!  We use numbers ~ 2 orders of magnitude below the maximum possible
!  values, as they may still get multiplied by some moderate numbers.
!
!  This value is a marker for some variable being uninitialized, and it is
!  tempting to replace the mechanism used here by NaN.
!  This may or may not work (need to reliably create NaN [including
!  REAL_PRECISION=double], some compilers seem to trap assignment of NaN
!  values, etc.
!  Also, there is no NaN concept for integers.
!
  real, parameter :: impossible=3.9085e37
  integer, parameter :: impossible_int=-max_int/100
!
! MPI
!
  integer, parameter :: root=0
!
!  Diagnostic variable types.
!
!  Values > 0 get maxed across all processors before any
!  transformation using mpi_reduce_max;
!  values < 0 get summed over all processors before any transformation
!  using mpi_reduce_sum;
!  value 0 causes the value simply to be used from the root processor.
!
  integer, parameter :: ilabel_max=-1,ilabel_sum=1,ilabel_save=0
  integer, parameter :: ilabel_max_sqrt=-2,ilabel_sum_sqrt=2
  integer, parameter :: ilabel_sum_log10=10, ilabel_sum_masked=11
  integer, parameter :: ilabel_max_dt=-3,ilabel_max_neg=-4
  integer, parameter :: ilabel_max_reciprocal=-5
  integer, parameter :: ilabel_integrate=3,ilabel_integrate_sqrt=30, ilabel_integrate_log10=40
  integer, parameter :: ilabel_surf=4
  integer, parameter :: ilabel_sum_par=5,ilabel_sum_sqrt_par=6, ilabel_sum_log10_par=20, ilabel_sum_plain=21
  integer, parameter :: ilabel_sum_weighted=7,ilabel_sum_weighted_sqrt=8
  integer, parameter :: ilabel_sum_lim=9,ilabel_complex=100
!
  real, parameter :: lntwo=0.69314718055995d0
!
!  first zeros of Bessel functions of order 0 and 1
!  k2bessel0 is the second zero of Bessel function of order 0
!
  real, parameter :: k1bessel0=2.4048255577, k1bessel1=3.8317060
  real, parameter :: k2bessel0=5.5200781
!
!  pi and its derivatives.
!
  real, parameter :: pi=3.14159265358979323846264338327950d0
  real, parameter :: pi_1=1./pi,pi4_1=(1.0)/(pi*pi*pi*pi),pi5_1=1.0/(pi*pi*pi*pi*pi)
  real, parameter :: sqrtpi=1.77245385090551602729816748334115d0
  real, parameter :: sqrt2=1.41421356237309504880168872420970d0
  real, parameter :: sqrt2pi=sqrt2*sqrtpi
  real, parameter :: four_pi_over_three=4.0/3.0*pi
  real, parameter :: onethird=1./3., twothird=2./3., fourthird=4./3., onesixth=1./6.
  real, parameter :: one_over_sqrt3=0.577350269189625764509148780501958d0
  real, parameter :: twopi = 6.2831853071795864769252867665590d0
  real, parameter :: dtor = pi/180.d0
!
!  Physical constants, taken from
!  http://physics.nist.gov/cuu/Constants/index.html.
!
  real(KIND=rkind8), parameter :: hbar_cgs=1.054571596d-27  ! [erg*s]
  real(KIND=rkind8), parameter :: k_B_cgs=1.3806505d-16     ! [erg/K]
  real(KIND=rkind8), parameter :: m_u_cgs=1.66053886d-24    ! [g]
  real(KIND=rkind8), parameter :: mu0_cgs=4*pi              ! [cgs]
  ! Better express R_cgs as a derived quantity (i.e. don't define here...)
  ! (Not done yet since it breaks the interstellar test)
  !real(KIND=rkind8), parameter :: R_cgs=k_B_cgs/m_u_cgs    ! [erg/g/K]
  real(KIND=rkind8), parameter :: R_cgs=8.3144D7            ! [erg/g/K]
  ! It would be better to specify the following masses in units of m_u:
  real(KIND=rkind8), parameter :: m_p_cgs=1.67262158d-24    ! [g]
  real(KIND=rkind8), parameter :: m_e_cgs=9.10938188d-28    ! [g]
  real(KIND=rkind8), parameter :: m_H_cgs=m_e_cgs+m_p_cgs   ! [g]
  real(KIND=rkind8), parameter :: eV_cgs=1.602176462d-12    ! [erg]
  real(KIND=rkind8), parameter :: sigmaSB_cgs=5.670400d-5   ! [erg/cm^2/s/K^4]
! unclear source (probably just guessing?)
  real(KIND=rkind8), parameter :: sigmaH_cgs=4.d-17         ! [cm^2]
  real(KIND=rkind8), parameter :: kappa_es_cgs=3.4d-1       ! [cm^2/g]
  real(KIND=rkind8), parameter :: c_light_cgs=2.99792458d10 ! [cm/s]
  real(KIND=rkind8), parameter :: G_Newton_cgs=6.6742d-8    ! [cm3/g/s2]
  real(KIND=rkind8), parameter :: density_scale_cgs=1.2435d21 ![cm] 403pc Reynolds 91, etc
  real(KIND=rkind8), parameter :: N_avogadro_cgs=6.022d23 ![1/mol]
!
  logical, parameter :: ALWAYS_FALSE=.false.
!
!  Data structure used to gather slice information from the various modules.
!
  type slice_data
    character (LEN=labellen) :: name
    integer :: index
    logical :: ready
    real, pointer, dimension (:,:) :: xy
    real, pointer, dimension (:,:) :: xz
    real, pointer, dimension (:,:) :: xz2
    real, pointer, dimension (:,:) :: yz
    real, pointer, dimension (:,:) :: xy2
    real, pointer, dimension (:,:) :: xy3
    real, pointer, dimension (:,:) :: xy4
    real, pointer, dimension (:,:) :: r
  endtype slice_data
!
!  Data structure used to allow module specific boundary conditions.
!
  type boundary_condition
    character (len=bclen) :: bcname
    integer :: ivar
    integer :: location
    logical :: done
    real :: value1
    real :: value2
  endtype boundary_condition
!
  integer, parameter :: iBC_X_TOP=1
  integer, parameter :: iBC_X_BOT=-1
  integer, parameter :: iBC_Y_TOP=2
  integer, parameter :: iBC_Y_BOT=-2
  integer, parameter :: iBC_Z_TOP=3
  integer, parameter :: iBC_Z_BOT=-3
  integer, parameter :: BOT=1, TOP=2, BOTH=3
!
!  Indices of rho, d rho/d x, d^2 rho/d x^2, d^6 rho/d x^6, d p/d x, s, d s/d x, &
!             d^2 s/d x^2, d^6 s/d x^6, d lnrho/d z in array reference_state.
!
  integer, parameter :: iref_rho=1, iref_grho=2, iref_d2rho=3, iref_d6rho=4, &
                        iref_gp=5, iref_s=6, iref_gs=7, iref_d2s=8, iref_d6s=9
  integer, parameter :: nref_vars=9
!
!  Symbolic constants for Yin-Yang grid.
!
  integer, parameter :: BILIN=1, BIQUAD=2, BICUB=3, QUADSPLINE=4, BIQUIN=5
!
!  Symbolic constants for Cubed Sphere grid.
!  The order of the patches is the same as in MATINS.
!
  integer, parameter :: XPLUS=1, YPLUS=2, XMINUS=3, YMINUS=4, ZPLUS=5, ZMINUS=6
!
  integer, parameter :: max_threads_possible = 200
  integer, parameter :: PERF_DIAGS=1, PERF_WSNAP=2, PERF_POWERSNAP=3, PERF_WSNAP_DOWN=4
  integer, parameter :: n_helperflags=4
  integer, parameter :: n_xy_specs_max=10,nk_max=10, nz_max=10
  integer, parameter :: mname=100
  integer, parameter :: mname_half=20

  !TP: strings to enums
  integer, parameter :: string_enum_unknown_string_string = 0
  integer, parameter :: string_enum_pde_string = 1
  integer, parameter :: string_enum_before_lanelastic_string = 2
  integer, parameter :: string_enum_calc_pencils_grid_string = 3
  !integer, parameter :: string_enum_pomxZ_pomyZ_phix_and_phiy_for_spherical_polars_string = 4
  integer, parameter :: string_enum_position_vector_for__string = 5
  integer, parameter :: string_enum_nonZcartesian_coordinates_string = 6
  !integer, parameter :: string_enum_radial_unit_vector_for_nonZcartesian_coordinates_string = 7
  integer, parameter :: string_enum_coZlatitudinal_unit_vector_for__string = 8
  integer, parameter :: string_enum_calc_pencils_hydro_linearized_string = 9
  integer, parameter :: string_enum_u2_pencil_not_calculated_string = 10
  integer, parameter :: string_enum_sij2_pencil_not_calculated_string = 11
  integer, parameter :: string_enum_uij5_pencil_not_calculated_string = 12
  integer, parameter :: string_enum_o2_or_oxu2_pencils_not_calculate_string = 13
  integer, parameter :: string_enum_ou_or_oxu_pencils_not_calculated_string = 14
  integer, parameter :: string_enum_calc_pencils_hydroZ_upwinding_advection_term_string = 15
  !integer, parameter :: string_enum_ugu_pencil_not_calculated_in_lconservative_case_string = 16
  integer, parameter :: string_enum_ugu2_pencil_not_calculated_string = 17
  integer, parameter :: string_enum_ujukl_pencils_not_calculated_string = 18
  integer, parameter :: string_enum_calc_pencils_hydroZ_call_gij_etc_string = 19
  integer, parameter :: string_enum_no_linearized_weno_transport_string = 20
  !integer, parameter :: string_enum_warning_Z_hydroZou_has_different_sign_than_relhel_string = 21
  integer, parameter :: string_enum_calc_pencils_density_string = 22
  !integer, parameter :: string_enum_uglnrho_not_available_for_linear_mass_density_string = 23
  integer, parameter :: string_enum_del6lnrho_for_linear_mass_density_string = 24
  integer, parameter :: string_enum_del6lnrho_strict_for_linear_mass_density_string = 25
  integer, parameter :: string_enum_hlnrho_linear_mass_density_string = 26
  integer, parameter :: string_enum_densityZiprocZitZmZnZ_string = 27
  integer, parameter :: string_enum_nans_in_ac_transformed_pencil_glnrho_string = 28
  integer, parameter :: string_enum_ugrho_for_logarithmic_mass_density_string = 29
  integer, parameter :: string_enum_del2rho_for_logarithmic_mass_density_string = 30
  integer, parameter :: string_enum_del6rho_for_logarithmic_mass_density_string = 31
  integer, parameter :: string_enum_calc_pencils_density_pnc_string = 32
  integer, parameter :: string_enum_rhos1_string = 33
  integer, parameter :: string_enum_glnrhos_string = 34
  integer, parameter :: string_enum_calc_pencils_eos_string = 35
  integer, parameter :: string_enum_rho1gpp_not_available_string = 36
  !integer, parameter :: string_enum_leos_localisothermal_for_ilnrho_ssZ_try_ilnrho_cs2_string = 37
  integer, parameter :: string_enum_rho1gpp_not_available_2_string = 38
  integer, parameter :: string_enum_del6ss_for_ilnrho_lntt_string = 39
  !integer, parameter :: string_enum_leos_isentropic_for_ilnrho_cs2Z_try_ilnrho_ss_string = 40
  !integer, parameter :: string_enum_temperature_not_needed_for_localisothermal_string = 41
  !integer, parameter :: string_enum_no_gradients_yet_for_localisothermal_string = 42
  !integer, parameter :: string_enum_entropy_not_needed_for_localisothermal_string = 43
  !integer, parameter :: string_enum_entropy_gradient_not_needed_for_localisothermal_string = 44
  integer, parameter :: string_enum_full_equation_of_state_for_ilnrho_cs2_string = 45
  integer, parameter :: string_enum_local_isothermal_case_for_ipp_ss_string = 46
  integer, parameter :: string_enum_isentropic_for_ZppZlnttZ_string = 47
  integer, parameter :: string_enum_local_isothermal_case_for_ipp_cs2_string = 48
  integer, parameter :: string_enum_del6ss_for_ilnrho_cs2_string = 49
  integer, parameter :: string_enum_geth_is_not_available_string = 50
  integer, parameter :: string_enum_del2eth_is_not_available_string = 51
  integer, parameter :: string_enum_eths_is_not_available_string = 52
  integer, parameter :: string_enum_geths_is_not_available_string = 53
  integer, parameter :: string_enum_hlntt_for_ilnrho_eth_or_irho_eth_string = 54
  integer, parameter :: string_enum_unknown_combination_of_eos_vars_string = 55
  integer, parameter :: string_enum_calc_pencils_energyZ_maxZadvec_cs2Z_Z_string = 56
  integer, parameter :: string_enum_carreau_string = 57
  integer, parameter :: string_enum_step_string = 58
  integer, parameter :: string_enum_getnu_non_newtonianZ_string = 59
  integer, parameter :: string_enum_no_such_nnewton_typeZ__string = 60
  !integer, parameter :: string_enum_powerZlaw_viscosity_with_luse_nu_rmn_profZt_for_other_than_spherical_coordinates_string = 61
  integer, parameter :: string_enum_calc_pencils_viscosity_string = 62
  !integer, parameter :: string_enum_shock_heating_not_implemented_for_lvisc_shock_simpleZt_string = 63
  integer, parameter :: string_enum_viscous_heating__string = 64
  integer, parameter :: string_enum_not_implemented_for_lvisc_hyper2_simplified_string = 65
  integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_simplified_string = 66
  integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_polar_string = 67
  integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_mesh_string = 68
  integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_csmesh_string = 69
  !integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_rho_nu_const_string = 70
  !integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_mu_const_strict_string = 71
  !integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_mu_const_strict_otf_string = 72
  !integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_nu_const_strict_string = 73
  integer, parameter :: string_enum_del2fjv_string = 74
  !integer, parameter :: string_enum_not_implemented_for_nonZcartesian_coordinates_string = 75
  !integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_rho_nu_const_aniso_string = 76
  !integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_nu_const_aniso_string = 77
  !integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_rho_nu_const_bulk_string = 78
  integer, parameter :: string_enum_not_implemented_for_lvisc_hyper3_nu_const_string = 79
  integer, parameter :: string_enum_viscous_heating_term__string = 80
  !integer, parameter :: string_enum_is_not_implemented_for_lvisc_smag_cross_simplified_string = 81
  integer, parameter :: string_enum_viscose_string = 82
  integer, parameter :: string_enum_get_bext_string = 83
  !integer, parameter :: string_enum_precession_of_external_field_for_curvilinear_coordinates_string = 84
  integer, parameter :: string_enum_step_scalar_string = 85
  integer, parameter :: string_enum_width_must_not_be_zero_string = 86
  integer, parameter :: string_enum_calc_pencils_magnetic_pencparZ_b_ext_Z__string = 87
  integer, parameter :: string_enum_calc_pencils_magnetic_pencparZ_logic_Z__string = 88
  integer, parameter :: string_enum_calc_pencils_magnetic_pencpar_string = 89
  integer, parameter :: string_enum_standard_string = 90
  integer, parameter :: string_enum_standard2_string = 91
  integer, parameter :: string_enum_logZswitchZon_string = 92
  integer, parameter :: string_enum_linearZsigma_string = 93
  integer, parameter :: string_enum_eta_table_string = 94
  integer, parameter :: string_enum_magnetic_after_boundary_string = 95
  integer, parameter :: string_enum_eta_table_not_yet_completed_string = 96
  integer, parameter :: string_enum_meanZfield_string = 97
  integer, parameter :: string_enum_need_leedot_as_auxZt_in_specialZdisp_current_string = 98
  !integer, parameter :: string_enum_calc_pencils_magneticZ_alfven_speed_is_imaginaryZ_string = 99
  integer, parameter :: string_enum_calc_pencils_magneticZ_itZ_itsubZ_iprocZ_string = 100
  integer, parameter :: string_enum_calc_pencils_magneticZ_mZ_yZmZZ_nZ_zZnZZ_string = 101
  integer, parameter :: string_enum_coulomb_gauge_needs_to_be_invoked_for_lam_string = 102
  integer, parameter :: string_enum_constant_string = 103
  integer, parameter :: string_enum_ionizationZequilibrium_string = 104
  integer, parameter :: string_enum_ionizationZyh_string = 105
  integer, parameter :: string_enum_set_ambipolar_diffusion_string = 106
  integer, parameter :: string_enum_no_such_ambipolar_diffusionZ__string = 107
  integer, parameter :: string_enum_duu_dt_string = 108
  integer, parameter :: string_enum_entered_string = 109
  integer, parameter :: string_enum_duu_dtZ_solve_string = 110
  !integer, parameter :: string_enum_ZaZa10ZZZ__xZ_ZZZa8ZZZZ_yZ_ZZZa8ZZZZ__zZ_ZZZa8ZZZZZ_string = 111
  integer, parameter :: string_enum_bcs_for__string = 112
  integer, parameter :: string_enum_ux_string = 113
  integer, parameter :: string_enum_uy_string = 114
  integer, parameter :: string_enum_uz_string = 115
  integer, parameter :: string_enum_sld_char_string = 116
  !integer, parameter :: string_enum_coriolis_cylindricalZ_omegaZ_string = 117
  !integer, parameter :: string_enum_coriolis_cylindricalZ_omegaZZthetaZ_string = 118
  !integer, parameter :: string_enum_coriolis_sphericalZ_omegaZ_string = 119
  !integer, parameter :: string_enum_coriolis_sphericalZ_omegaZthetaZphiZ_string = 120
  !integer, parameter :: string_enum_precessionZ_omega_precessionZ_string = 121
  !integer, parameter :: string_enum_coriolis_cartesian_xaxisZ_coriolis_forceZ_omegaZ_thetaZ_string = 122
  !integer, parameter :: string_enum_coriolis_cartesianZ_add_coriolis_forceZ_omegaZ_string = 123
  !integer, parameter :: string_enum_coriolis_cartesianZ_add_centrifugal_forceZ_omegaZ_string = 124
  !integer, parameter :: string_enum_coriolis_cartesianZ_centrifugal_force_amplitudeZ_amp_centforceZ_string = 125
  !integer, parameter :: string_enum_coriolis_cartesianZ_coriolis_forceZ_omegaZ_thetaZ_string = 126
  integer, parameter :: string_enum_coriolis_xdepZ_ampl_omegaZ_string = 127
  integer, parameter :: string_enum_duu_dtZ_maxZadvec_uuZ_Z_string = 128
  integer, parameter :: string_enum_nothing_string = 129
  integer, parameter :: string_enum_linear_string = 130
  integer, parameter :: string_enum_inverse_string = 131
  integer, parameter :: string_enum_current_string = 132
  integer, parameter :: string_enum_bs04_string = 133
  integer, parameter :: string_enum_bs04c_string = 134
  integer, parameter :: string_enum_bs04c1_string = 135
  integer, parameter :: string_enum_bs04m_string = 136
  integer, parameter :: string_enum_hp09_string = 137
  integer, parameter :: string_enum_sx_string = 138
  integer, parameter :: string_enum_solar_dc99_string = 139
  integer, parameter :: string_enum_vertical_shear_string = 140
  integer, parameter :: string_enum_vertical_compression_string = 141
  integer, parameter :: string_enum_remove_vertical_shear_string = 142
  integer, parameter :: string_enum_vertical_shear_x_string = 143
  integer, parameter :: string_enum_vertical_shear_x_sinz_string = 144
  integer, parameter :: string_enum_vertical_shear_z_string = 145
  integer, parameter :: string_enum_vertical_shear_z2_string = 146
  integer, parameter :: string_enum_vertical_shear_linear_string = 147
  integer, parameter :: string_enum_tachocline_string = 148
  integer, parameter :: string_enum_solar_simple_string = 149
  integer, parameter :: string_enum_radial_uniform_shear_string = 150
  integer, parameter :: string_enum_breeze_string = 151
  integer, parameter :: string_enum_slow_wind_string = 152
  integer, parameter :: string_enum_radial_shear_string = 153
  integer, parameter :: string_enum_radial_shear_damp_string = 154
  integer, parameter :: string_enum_damp_corona_string = 155
  integer, parameter :: string_enum_damp_horiz_vel_string = 156
  integer, parameter :: string_enum_latitudinal_shear_string = 157
  integer, parameter :: string_enum_damp_jets_string = 158
  integer, parameter :: string_enum_spokeZlikeZnssl_string = 159
  integer, parameter :: string_enum_uumz_profile_string = 160
  integer, parameter :: string_enum_omega_profile_string = 161
  integer, parameter :: string_enum_zero_string = 162
  integer, parameter :: string_enum_0_string = 163
  integer, parameter :: string_enum_initialZcondition_string = 164
  integer, parameter :: string_enum_finished_string = 165
  integer, parameter :: string_enum_dlnrho_dt_string = 166
  integer, parameter :: string_enum_dlnrho_dtZ_solve_string = 167
  integer, parameter :: string_enum_lnrho_string = 168
  integer, parameter :: string_enum_surface_z_string = 169
  integer, parameter :: string_enum_mass_sourceZ_cs20Zcs0Z_string = 170
  integer, parameter :: string_enum_mass_source_string = 171
  integer, parameter :: string_enum_mass_source_with_no_profile_string = 172
  integer, parameter :: string_enum_exponential_string = 173
  integer, parameter :: string_enum_bump_string = 174
  integer, parameter :: string_enum_bump2_string = 175
  integer, parameter :: string_enum_bumpr_string = 176
  integer, parameter :: string_enum_bumpx_string = 177
  integer, parameter :: string_enum_sphZstepZdown_string = 178
  integer, parameter :: string_enum_const_string = 179
  integer, parameter :: string_enum_cylindric_string = 180
  integer, parameter :: string_enum_no_such_mass_source_profileZ__string = 181
  integer, parameter :: string_enum_dlnrho_dtZ_diffrhoZ_string = 182
  integer, parameter :: string_enum_dlnrho_dtZ_diffrho_shockZ_string = 183
  integer, parameter :: string_enum_dlnrho_dtZ_diffrho_hyper3Z_string = 184
  integer, parameter :: string_enum_dlnrho_dtZ_diffrho_hyper3_meshZ_string = 185
  integer, parameter :: string_enum_dlnrho_dtZ_diffrho_hyper3ZZdxZdyZdzZZ_string = 186
  integer, parameter :: string_enum_dlnrho_dtZ_diffrho_hyper3_strictZ_string = 187
  integer, parameter :: string_enum_dlnrho_dtZ_maxZdiffus_diffrho_Z_Z_string = 188
  integer, parameter :: string_enum_dlnrho_dtZ_maxZdiffus_diffrho3Z_Z_string = 189
  integer, parameter :: string_enum_before_calc_diagnostics_string = 190
  !integer, parameter :: string_enum_denergy_dtZ_adding_global_pressure_gradient_force_string = 191
  integer, parameter :: string_enum_daa_dt_string = 192
  integer, parameter :: string_enum_daa_dtZ_solve_string = 193
  integer, parameter :: string_enum_ax_string = 194
  integer, parameter :: string_enum_ay_string = 195
  integer, parameter :: string_enum_az_string = 196
  integer, parameter :: string_enum_bx_string = 197
  integer, parameter :: string_enum_by_string = 198
  integer, parameter :: string_enum_bz_string = 199
  integer, parameter :: string_enum_jx_string = 200
  integer, parameter :: string_enum_jy_string = 201
  integer, parameter :: string_enum_jz_string = 202
  integer, parameter :: string_enum_daa_dtZ_iresistivityZ_string = 203
  integer, parameter :: string_enum_two_step_string = 204
  integer, parameter :: string_enum_twoZstep_string = 205
  integer, parameter :: string_enum_two_step2_string = 206
  integer, parameter :: string_enum_twoZstep2_string = 207
  integer, parameter :: string_enum_ac_transformed_pencil_r_mnZ_eta_rZ_geta_r_string = 208
  integer, parameter :: string_enum_Z1pZ5e11Z3Z_string = 209
  !integer, parameter :: string_enum_must_have_weyl_gauge_for_anomalous_resistivity_string = 210
  integer, parameter :: string_enum_must_put_lua_as_auxZt_for_advective_gauge2_string = 211
  integer, parameter :: string_enum_daa_dtZ_use_upwinding_in_advection_term_string = 212
  integer, parameter :: string_enum_tZdep_string = 213
  integer, parameter :: string_enum_zZdep_string = 214
  integer, parameter :: string_enum_daa_dtZ_hall_termZ_string = 215
  integer, parameter :: string_enum_daa_dtZ_maxZadvec_hallZ_Z_string = 216
  integer, parameter :: string_enum_daa_dtZ_battery_termZ_string = 217
  integer, parameter :: string_enum_daa_dtZ_maxZbattery_termZ_Z_string = 218
  integer, parameter :: string_enum_daa_dtZ_height_etaZeta_outZlhaloxZ_string = 219
  integer, parameter :: string_enum_calc_tau_aa_exteriorZ_tauZ_string = 220
  integer, parameter :: string_enum_fZl1Zl2ZmZnZiexZiezZZZdadt_is_set_string = 221
  integer, parameter :: string_enum_aaZdat_string = 222
  integer, parameter :: string_enum_bbZdat_string = 223
  integer, parameter :: string_enum_jjZdat_string = 224
  integer, parameter :: string_enum_del2aZdat_string = 225
  integer, parameter :: string_enum_jxbrZdat_string = 226
  integer, parameter :: string_enum_jxbZdat_string = 227
  integer, parameter :: string_enum_dfZdat_string = 228
  integer, parameter :: string_enum_particles_pde_pencil_string = 229
  integer, parameter :: string_enum_shepherd_neighbour_pencil_string = 230
  integer, parameter :: string_enum_not_implementedZ__string = 231
  integer, parameter :: string_enum_interpolate_quantities_string = 232
  integer, parameter :: string_enum_interp_is_not_implementedZ__string = 233
  integer, parameter :: string_enum_cleanup_interpolated_quantities_string = 234
  integer, parameter :: string_enum_pdeZ_maxadvec_contains_a_nan_at_iprocZ_string = 235
  integer, parameter :: string_enum_advec_cs2__Z_string = 236
  integer, parameter :: string_enum_set_dt1_max_string = 237
  integer, parameter :: string_enum__string = 238
  integer, parameter :: string_enum_sum_mn_string = 239
  integer, parameter :: string_enum_not_implemented_for_cylindrical_string = 240
  integer, parameter :: string_enum_rhs_cpu_string = 241
  integer, parameter :: string_enum_end_of_mn_loop_string = 242
  integer, parameter :: string_enum_ZtvartZdat_string = 243
  integer, parameter :: string_enum_unknown_string = 244
  integer, parameter :: string_enum_append_string = 245
  integer, parameter :: string_enum_Z4f14Z7Z_string = 246

endmodule Cparam
