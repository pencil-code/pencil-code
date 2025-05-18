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
 ! integer, parameter :: rkind16 = selected_real_kind(33, 4931) ! 16-byte real kind - not accepted by all compilers
  integer, parameter :: rkind16 = rkind8
!
  include 'cparam.local'
!
  integer, parameter :: nx=nxgrid/nprocx,ny=nygrid/nprocy,nz=nzgrid/nprocz,nyz=ny*nz
  integer, parameter :: max_n = max(nx,max(ny,nz))
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
  real, parameter :: sqrt21=1./sqrt2
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
  real(KIND=rkind8), parameter :: alpha_fine=7.2973525643d-3
  real(KIND=rkind8), parameter :: sigma_Thomson_cgs=6.652458732160d-25 ![cm^2]
  real(KIND=rkind8), parameter :: e_cgs=4.8032047d-10  ![statcoulombs]
  real(KIND=rkind8), parameter :: Chypercharge=41./12. !(nondimensional)
  real(KIND=rkind8), parameter :: mass_zboson=7.48e-18 !(in Mpl, not reduced)
  real(KIND=rkind8), parameter :: mass_zboson_GeV=91.2 ![in GeV]
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

integer, parameter :: enum_unknown_string_string = 0
integer, parameter :: enum_pde_string = 1
integer, parameter :: enum_before_lanelastic_string = 2
integer, parameter :: enum_calc_pencils_grid_string = 3
integer, parameter :: enum_position_vector_for__string = 4
integer, parameter :: enum_nonZcartesian_coordinates_string = 5
integer, parameter :: enum_coZlatitudinal_unit_vector_for__string = 6
integer, parameter :: enum_calc_pencils_hydro_linearized_string = 7
integer, parameter :: enum_u2_pencil_not_calculated_string = 8
integer, parameter :: enum_sij2_pencil_not_calculated_string = 9
integer, parameter :: enum_uij5_pencil_not_calculated_string = 10
integer, parameter :: enum_o2_or_oxu2_pencils_not_calculate_string = 11
integer, parameter :: enum_ou_or_oxu_pencils_not_calculated_string = 12
integer, parameter :: enum_ugu2_pencil_not_calculated_string = 13
integer, parameter :: enum_ujukl_pencils_not_calculated_string = 14
integer, parameter :: enum_calc_pencils_hydroZ_call_gij_etc_string = 15
integer, parameter :: enum_no_linearized_weno_transport_string = 16
integer, parameter :: enum_calc_pencils_hydro_nonlinear_string = 17
integer, parameter :: enum_calc_pencils_density_string = 18
integer, parameter :: enum_del6lnrho_for_linear_mass_density_string = 19
integer, parameter :: enum_hlnrho_linear_mass_density_string = 20
integer, parameter :: enum_densityZiprocZitZmZnZ_string = 21
integer, parameter :: enum_nans_in_ac_transformed_pencil_glnrho_string = 22
integer, parameter :: enum_ugrho_for_logarithmic_mass_density_string = 23
integer, parameter :: enum_del2rho_for_logarithmic_mass_density_string = 24
integer, parameter :: enum_del6rho_for_logarithmic_mass_density_string = 25
integer, parameter :: enum_calc_pencils_density_pnc_string = 26
integer, parameter :: enum_rhos1_string = 27
integer, parameter :: enum_glnrhos_string = 28
integer, parameter :: enum_calc_pencils_eos_string = 29
integer, parameter :: enum_rho1gpp_not_available_string = 30
integer, parameter :: enum_rho1gpp_not_available_2_string = 31
integer, parameter :: enum_del6ss_for_ilnrho_lntt_string = 32
integer, parameter :: enum_no_gradients_yet_for_localisothermal_string = 33
integer, parameter :: enum_entropy_not_needed_for_localisothermal_string = 34
integer, parameter :: enum_full_equation_of_state_for_ilnrho_cs2_string = 35
integer, parameter :: enum_local_isothermal_case_for_ipp_ss_string = 36
integer, parameter :: enum_isentropic_for_ZppZlnttZ_string = 37
integer, parameter :: enum_local_isothermal_case_for_ipp_cs2_string = 38
integer, parameter :: enum_del6ss_for_ilnrho_cs2_string = 39
integer, parameter :: enum_geth_is_not_available_string = 40
integer, parameter :: enum_del2eth_is_not_available_string = 41
integer, parameter :: enum_eths_is_not_available_string = 42
integer, parameter :: enum_geths_is_not_available_string = 43
integer, parameter :: enum_hlntt_for_ilnrho_eth_or_irho_eth_string = 44
integer, parameter :: enum_unknown_combination_of_eos_vars_string = 45
integer, parameter :: enum_calc_pencils_energyZ_maxZadvec_cs2Z_Z_string = 46
integer, parameter :: enum_carreau_string = 47
integer, parameter :: enum_step_string = 48
integer, parameter :: enum_getnu_non_newtonianZ_string = 49
integer, parameter :: enum_no_such_nnewton_typeZ__string = 50
integer, parameter :: enum_calc_pencils_viscosity_string = 51
integer, parameter :: enum_viscous_heating__string = 52
integer, parameter :: enum_not_implemented_for_lvisc_hyper3_polar_string = 53
integer, parameter :: enum_not_implemented_for_lvisc_hyper3_mesh_string = 54
integer, parameter :: enum_not_implemented_for_lvisc_hyper3_csmesh_string = 55
integer, parameter :: enum_del2fjv_string = 56
integer, parameter :: enum_viscous_heating_term__string = 57
integer, parameter :: enum_viscose_string = 58
integer, parameter :: enum_init_uu_string = 59
integer, parameter :: enum_get_bext_string = 60
integer, parameter :: enum_step_scalar_string = 61
integer, parameter :: enum_width_must_not_be_zero_string = 62
integer, parameter :: enum_calc_pencils_magnetic_pencparZ_b_ext_Z__string = 63
integer, parameter :: enum_calc_pencils_magnetic_pencparZ_logic_Z__string = 64
integer, parameter :: enum_calc_pencils_magnetic_pencpar_string = 65
integer, parameter :: enum_uuadvec_gaa_for_spherical_coordinates_string = 66
integer, parameter :: enum_constant_string = 67
integer, parameter :: enum_ionizationZequilibrium_string = 68
integer, parameter :: enum_ionizationZyh_string = 69
integer, parameter :: enum_set_ambipolar_diffusion_string = 70
integer, parameter :: enum_no_such_ambipolar_diffusionZ__string = 71
integer, parameter :: enum_duu_dt_string = 72
integer, parameter :: enum_entered_string = 73
integer, parameter :: enum_duu_dtZ_solve_string = 74
integer, parameter :: enum_bcs_for__string = 75
integer, parameter :: enum_ux_string = 76
integer, parameter :: enum_uy_string = 77
integer, parameter :: enum_uz_string = 78
integer, parameter :: enum_sld_char_string = 79
integer, parameter :: enum_coriolis_cylindricalZ_omegaZ_string = 80
integer, parameter :: enum_coriolis_cylindricalZ_omegaZZthetaZ_string = 81
integer, parameter :: enum_coriolis_cylindrical_string = 82
integer, parameter :: enum_coriolis_sphericalZ_omegaZ_string = 83
integer, parameter :: enum_coriolis_sphericalZ_omegaZthetaZphiZ_string = 84
integer, parameter :: enum_coriolis_spherical_string = 85
integer, parameter :: enum_for_omega_not_aligned_with_z_or_y_axis_string = 86
integer, parameter :: enum_precessionZ_omega_precessionZ_string = 87
integer, parameter :: enum_coriolis_cartesian_string = 88
integer, parameter :: enum_if_omega_has_y_component_string = 89
integer, parameter :: enum_coriolis_xdepZ_ampl_omegaZ_string = 90
integer, parameter :: enum_duu_dtZ_maxZadvec_uuZ_Z_string = 91
integer, parameter :: enum_nothing_string = 92
integer, parameter :: enum_linear_string = 93
integer, parameter :: enum_inverse_string = 94
integer, parameter :: enum_current_string = 95
integer, parameter :: enum_lmagnetic_must_be_true_string = 96
integer, parameter :: enum_bs04_string = 97
integer, parameter :: enum_bs04c_string = 98
integer, parameter :: enum_bs04c1_string = 99
integer, parameter :: enum_bs04m_string = 100
integer, parameter :: enum_hp09_string = 101
integer, parameter :: enum_sx_string = 102
integer, parameter :: enum_solar_dc99_string = 103
integer, parameter :: enum_vertical_shear_string = 104
integer, parameter :: enum_vertical_compression_string = 105
integer, parameter :: enum_remove_vertical_shear_string = 106
integer, parameter :: enum_vertical_shear_x_string = 107
integer, parameter :: enum_vertical_shear_x_sinz_string = 108
integer, parameter :: enum_vertical_shear_z_string = 109
integer, parameter :: enum_vertical_shear_z2_string = 110
integer, parameter :: enum_vertical_shear_linear_string = 111
integer, parameter :: enum_tachocline_string = 112
integer, parameter :: enum_solar_simple_string = 113
integer, parameter :: enum_radial_uniform_shear_string = 114
integer, parameter :: enum_breeze_string = 115
integer, parameter :: enum_slow_wind_string = 116
integer, parameter :: enum_radial_shear_string = 117
integer, parameter :: enum_radial_shear_damp_string = 118
integer, parameter :: enum_damp_corona_string = 119
integer, parameter :: enum_damp_horiz_vel_string = 120
integer, parameter :: enum_latitudinal_shear_string = 121
integer, parameter :: enum_damp_jets_string = 122
integer, parameter :: enum_spokeZlikeZnssl_string = 123
integer, parameter :: enum_uumz_profile_string = 124
integer, parameter :: enum_omega_profile_string = 125
integer, parameter :: enum_zero_string = 126
integer, parameter :: enum_0_string = 127
integer, parameter :: enum_initialZcondition_string = 128
integer, parameter :: enum_finished_string = 129
integer, parameter :: enum_dlnrho_dt_string = 130
integer, parameter :: enum_dlnrho_dtZ_solve_string = 131
integer, parameter :: enum_lnrho_string = 132
integer, parameter :: enum_surface_z_string = 133
integer, parameter :: enum_mass_sourceZ_cs20Zcs0Z_string = 134
integer, parameter :: enum_mass_source_string = 135
integer, parameter :: enum_mass_source_with_no_profile_string = 136
integer, parameter :: enum_exponential_string = 137
integer, parameter :: enum_bump_string = 138
integer, parameter :: enum_bump2_string = 139
integer, parameter :: enum_bumpr_string = 140
integer, parameter :: enum_bumpx_string = 141
integer, parameter :: enum_sphZstepZdown_string = 142
integer, parameter :: enum_const_string = 143
integer, parameter :: enum_cylindric_string = 144
integer, parameter :: enum_no_such_mass_source_profileZ__string = 145
integer, parameter :: enum_dlnrho_dtZ_diffrhoZ_string = 146
integer, parameter :: enum_dlnrho_dtZ_diffrho_shockZ_string = 147
integer, parameter :: enum_dlnrho_dtZ_diffrho_hyper3Z_string = 148
integer, parameter :: enum_dlnrho_dtZ_diffrho_hyper3_meshZ_string = 149
integer, parameter :: enum_dlnrho_dtZ_diffrho_hyper3ZZdxZdyZdzZZ_string = 150
integer, parameter :: enum_dlnrho_dtZ_diffrho_hyper3_strictZ_string = 151
integer, parameter :: enum_dlnrho_dtZ_maxZdiffus_diffrho_Z_Z_string = 152
integer, parameter :: enum_dlnrho_dtZ_maxZdiffus_diffrho3Z_Z_string = 153
integer, parameter :: enum_before_calc_diagnostics_string = 154
integer, parameter :: enum_daa_dt_string = 155
integer, parameter :: enum_daa_dtZ_solve_string = 156
integer, parameter :: enum_ax_string = 157
integer, parameter :: enum_ay_string = 158
integer, parameter :: enum_az_string = 159
integer, parameter :: enum_bx_string = 160
integer, parameter :: enum_by_string = 161
integer, parameter :: enum_bz_string = 162
integer, parameter :: enum_jx_string = 163
integer, parameter :: enum_jy_string = 164
integer, parameter :: enum_jz_string = 165
integer, parameter :: enum_daa_dtZ_iresistivityZ_string = 166
integer, parameter :: enum_two_step_string = 167
integer, parameter :: enum_twoZstep_string = 168
integer, parameter :: enum_two_step2_string = 169
integer, parameter :: enum_twoZstep2_string = 170
integer, parameter :: enum_Z1pZ5e11Z3Z_string = 171
integer, parameter :: enum_eta_shell_string = 172
integer, parameter :: enum_daa_dtZ_use_upwinding_in_advection_term_string = 173
integer, parameter :: enum_tZdep_string = 174
integer, parameter :: enum_zZdep_string = 175
integer, parameter :: enum_daa_dtZ_hall_termZ_string = 176
integer, parameter :: enum_daa_dtZ_maxZadvec_hallZ_Z_string = 177
integer, parameter :: enum_daa_dtZ_battery_termZ_string = 178
integer, parameter :: enum_daa_dtZ_maxZbattery_termZ_Z_string = 179
integer, parameter :: enum_daa_dtZ_height_etaZeta_outZlhaloxZ_string = 180
integer, parameter :: enum_calc_tau_aa_exteriorZ_tauZ_string = 181
integer, parameter :: enum_fZl1Zl2ZmZnZiexZiezZZZdadt_is_set_string = 182
integer, parameter :: enum_aaZdat_string = 183
integer, parameter :: enum_bbZdat_string = 184
integer, parameter :: enum_jjZdat_string = 185
integer, parameter :: enum_del2aZdat_string = 186
integer, parameter :: enum_jxbrZdat_string = 187
integer, parameter :: enum_jxbZdat_string = 188
integer, parameter :: enum_dfZdat_string = 189
integer, parameter :: enum_dspecial_dtZ_solve_dspecial_dt_string = 190
integer, parameter :: enum_rhs_cpu_string = 191
integer, parameter :: enum_end_of_mn_loop_string = 192
integer, parameter :: enum_denergy_dtZ_solve_denergy_dt_string = 193
integer, parameter :: enum_ss_string = 194
integer, parameter :: enum_denergy_dtZ_lnttZcs2Zcp1Z_string = 195
integer, parameter :: enum_ac_transformed_pencil_fpres_Z_string = 196
integer, parameter :: enum_denergy_dt_string = 197
integer, parameter :: enum__string = 198
integer, parameter :: enum_calc_heatcondZ_hcond0Z_string = 199
integer, parameter :: enum_calc_heatcondZ_lgravzZ_string = 200
integer, parameter :: enum_calc_heatcondZ_fbotZftopZ_string = 201
integer, parameter :: enum_calc_heatcond_string = 202
integer, parameter :: enum_nans_in_ac_transformed_pencil_glntt_string = 203
integer, parameter :: enum_calc_heatcondZ__string = 204
integer, parameter :: enum_calc_heatcondZ_nans_in_rho1_string = 205
integer, parameter :: enum_calc_heatcondZ_nans_in_del2ss_string = 206
integer, parameter :: enum_calc_heatcondZ_nans_in_hcond_string = 207
integer, parameter :: enum_calc_heatcondZ_nans_in_1Zhcond_string = 208
integer, parameter :: enum_calc_heatcondZ_nans_in_glhc_string = 209
integer, parameter :: enum_calc_heatcondZ_nans_in_chix_string = 210
integer, parameter :: enum_calc_heatcondZ_nans_in_glnthcond_string = 211
integer, parameter :: enum_calc_heatcondZ_nans_in_g2_string = 212
integer, parameter :: enum_chiZdat_string = 213
integer, parameter :: enum_hcondZdat_string = 214
integer, parameter :: enum_glhcZdat_string = 215
integer, parameter :: enum_heatcondZdat_string = 216
integer, parameter :: enum_calc_heatcondZ_added_thdiff_string = 217
integer, parameter :: enum_calc_heatcond_constkZ_hcondZ_string = 218
integer, parameter :: enum_calc_heatcond_constkZ_added_thdiff_string = 219
integer, parameter :: enum_calc_heatcond_sfluctZ_chi_tZ_string = 220
integer, parameter :: enum_calc_heatcond_constchiZ_chiZ_string = 221
integer, parameter :: enum_calc_heatcond_constchiZ_added_thdiff_string = 222
integer, parameter :: enum_calc_heatcond_cspeed_chiZ_chiZ_string = 223
integer, parameter :: enum_calc_heatcond_cspeed_chiZ_added_thdiff_string = 224
integer, parameter :: enum_calc_heatcond_sqrtrhochiZ_chi_rhoZ_string = 225
integer, parameter :: enum_calc_heatcond_sqrtrhochiZ_added_thdiff_string = 226
integer, parameter :: enum_calc_heatcond_shockZ_chi_shockZ_string = 227
integer, parameter :: enum_calc_heatcond_shockZ_added_thdiff_string = 228
integer, parameter :: enum_calc_heatcond_shock_profrZ_added_thdiff_string = 229
integer, parameter :: enum_calc_heatcond_hyper3Z_chi_hyper3Z_string = 230
integer, parameter :: enum_calc_heatcond_hyper3Z_added_thdiff_string = 231
integer, parameter :: enum_spitzerZdat_string = 232
integer, parameter :: enum_viscousZdat_string = 233
integer, parameter :: enum_enter_heatcond_hubeny_string = 234
integer, parameter :: enum_calc_heatcond_kramersZ_nans_in_rho1_string = 235
integer, parameter :: enum_calc_heatcond_kramersZ_nans_in_kZrho_string = 236
integer, parameter :: enum_calc_heatcond_kramersZ_nans_in_del2ss_string = 237
integer, parameter :: enum_calc_heatcond_kramersZ_nans_in_tt_string = 238
integer, parameter :: enum_calc_heatcond_kramersZ_nans_in_glnt_string = 239
integer, parameter :: enum_calc_heatcond_kramersZ_nans_in_g2_string = 240
integer, parameter :: enum_calc_heatcond_kramersZ_nans_in_thdiff_string = 241
integer, parameter :: enum_calc_heatcond_kramersZ_added_thdiff_string = 242
integer, parameter :: enum_calc_heatcond_chitZ_chi_t0Z_string = 243
integer, parameter :: enum_calc_heatcond_chitZ_chi_t1Z_string = 244
integer, parameter :: enum_calc_heatcond_smagorinskyZ_nans_in_rho1_string = 245
integer, parameter :: enum_calc_heatcond_smagorinskyZ_nans_in_chix_string = 246
integer, parameter :: enum_calc_heatcond_smagorinskyZ_nans_in_tt_string = 247
integer, parameter :: enum_calc_heatcond_smagorinskyZ_nans_in_glnt_string = 248
integer, parameter :: enum_calc_heatcond_smagorinskyZ_nans_in_g2_string = 249
integer, parameter :: enum_calc_heatcond_smagorinskyZ_added_thdiff_string = 250
integer, parameter :: enum_newtonZdat_string = 251
integer, parameter :: enum_calc_heat_cool_rtv_string = 252
integer, parameter :: enum_for_pretend_lntt_Z_t_string = 253
integer, parameter :: enum_cgs_string = 254
integer, parameter :: enum_rtvZdat_string = 255
integer, parameter :: enum_calc_heatcond_hyper3_polarZ_chi_hyper3Z_string = 256
integer, parameter :: enum_calc_heatcond_hyper3_meshZ_chi_hyper3Z_string = 257
integer, parameter :: enum_gaussianZz_string = 258
integer, parameter :: enum_linZz_string = 259
integer, parameter :: enum_sinZz_string = 260
integer, parameter :: enum_surface_x_string = 261
integer, parameter :: enum_twoZlayer_string = 262
integer, parameter :: enum_squareZwell_string = 263
integer, parameter :: enum_cubic_step_string = 264
integer, parameter :: enum_cubic_step_topbot_string = 265
integer, parameter :: enum_surface_pp_string = 266
integer, parameter :: enum_plain_string = 267
integer, parameter :: enum_corona_string = 268
integer, parameter :: enum_temp_string = 269
integer, parameter :: enum_get_cool_generalZ_cs20Zcs2coolZ_string = 270
integer, parameter :: enum_temp2_string = 271
integer, parameter :: enum_rho_cs2_string = 272
integer, parameter :: enum_twoZlayerZmean_string = 273
integer, parameter :: enum_get_cool_general_string = 274
integer, parameter :: enum_no_such_cooltypeZ__string = 275
integer, parameter :: enum_cooling_profileZz2ZwcoolZcs2coolZ_string = 276
integer, parameter :: enum_gaussian_string = 277
integer, parameter :: enum_step2_string = 278
integer, parameter :: enum_surfcool_string = 279
integer, parameter :: enum_volheat_surfcool_string = 280
integer, parameter :: enum_cs2Zrho_string = 281
integer, parameter :: enum_get_heat_cool_gravr_string = 282
integer, parameter :: enum_no_such_heattypeZ__string = 283
integer, parameter :: enum_heatZdat_string = 284
integer, parameter :: enum_cs2_string = 285
integer, parameter :: enum_tempZrho_string = 286
integer, parameter :: enum_entropy_string = 287
integer, parameter :: enum_pressure_string = 288
integer, parameter :: enum_shell_string = 289
integer, parameter :: enum_calc_heat_coolZ_deltat_poleqZ_string = 290
integer, parameter :: enum_ac_transformed_pencil_rcyl_mnZ_string = 291
integer, parameter :: enum_ac_transformed_pencil_z_mnZ_string = 292
integer, parameter :: enum_shell2_string = 293
integer, parameter :: enum_shell3_string = 294
integer, parameter :: enum_shell_mean_yz_string = 295
integer, parameter :: enum_shell_mean_yz2_string = 296
integer, parameter :: enum_shell_mean_downflow_string = 297
integer, parameter :: enum_latheat_string = 298
integer, parameter :: enum_shellZlatheat_string = 299
integer, parameter :: enum_shellZlatss_string = 300
integer, parameter :: enum_top_layer_string = 301
integer, parameter :: enum_calc_heat_cool_gravx_cartesian_string = 302
integer, parameter :: enum_eoscalc_pencil_string = 303
integer, parameter :: enum_eoscalc_point_string = 304
integer, parameter :: enum_thermodynamic_variable_combination_string = 305
integer, parameter :: enum_calc_tau_ss_exteriorZ_tauZ_string = 306
integer, parameter :: enum_initialZtemperature_string = 307
integer, parameter :: enum_daa_dtZ_maxZdiffus_etaZ__Z_string = 308
integer, parameter :: enum_daa_dtZ_maxZdiffus_eta2Z_Z_string = 309
integer, parameter :: enum_daa_dtZ_maxZdiffus_eta3Z_Z_string = 310
integer, parameter :: enum_pdeZ_maxadvec_contains_a_nan_at_iprocZ_string = 311
integer, parameter :: enum_advec_cs2__Z_string = 312
integer, parameter :: enum_set_dt1_max_string = 313
integer, parameter :: enum_cst_string = 314
integer, parameter :: enum_wolfire_string = 315
integer, parameter :: enum_wolfire_min_string = 316
integer, parameter :: enum_thermalZhs_string = 317
integer, parameter :: enum_off_string = 318
integer, parameter :: enum_calc_heat_cool_interstellarZ_enter_string = 319
integer, parameter :: enum_calc_pencils_dustdensity_string = 320
integer, parameter :: enum_average_string = 321
integer, parameter :: enum_neighbor_string = 322
integer, parameter :: enum_neighbor_asymmetric_string = 323
integer, parameter :: enum_dustdensityZcoag_kernel_string = 324
integer, parameter :: enum_no_such_self_collisionsZ__string = 325
integer, parameter :: enum_coag_kernel_string = 326
integer, parameter :: enum_this_should_never_happen_string = 327
integer, parameter :: enum_duud_dtZ_solve_duud_dt_string = 328
integer, parameter :: enum_udx_string = 329
integer, parameter :: enum_udy_string = 330
integer, parameter :: enum_udz_string = 331
integer, parameter :: enum_epstein_cst_string = 332
integer, parameter :: enum_epstein_cst_b_string = 333
integer, parameter :: enum_stokes_cst_tausd_string = 334
integer, parameter :: enum_stokes_varmass_string = 335
integer, parameter :: enum_epstein_var_string = 336
integer, parameter :: enum_epstein_gaussian_z_string = 337
integer, parameter :: enum_get_stoppingtime_string = 338
integer, parameter :: enum_no_such_drag_lawZ__string = 339
integer, parameter :: enum_duud_dtZ_add_coriolis_forceZ_omegaZ_string = 340
integer, parameter :: enum_duud_dtZ_coriolis_forceZ_omegaZthetaZ_string = 341
integer, parameter :: enum_duud_dtZ_maxZdiffus_nudZ_Z_string = 342
integer, parameter :: enum_duud_dtZ_calculate_diagnostic_valuesZZZ_string = 343
integer, parameter :: enum_dndmd_dtZ_solve_dnd_dtZ_dmd_dtZ_dmi_dt_string = 344
integer, parameter :: enum_nd_string = 345
integer, parameter :: enum_md_string = 346
integer, parameter :: enum_mi_string = 347
integer, parameter :: enum_simplified_string = 348
integer, parameter :: enum_pscalar_string = 349
integer, parameter :: enum_ice_string = 350
integer, parameter :: enum_aerosol_string = 351
integer, parameter :: enum_condensing_species_test_string = 352
integer, parameter :: enum_condensing_species_string = 353
integer, parameter :: enum_hatZomZtZ_string = 354
integer, parameter :: enum_cosZomZtZ_string = 355
integer, parameter :: enum_get_mfluxcond_string = 356
integer, parameter :: enum_no_such_dust_chemistryZ__string = 357
integer, parameter :: enum_lin_radius_string = 358
integer, parameter :: enum_log_radius_string = 359
integer, parameter :: enum_log_mass_string = 360
integer, parameter :: enum_dndmd_dt_string = 361
integer, parameter :: enum_not_implemented_for_llog_massbins_yet_string = 362
integer, parameter :: enum_register_dustvelocity_string = 363
integer, parameter :: enum_no_valid_dust_binning_string = 364
integer, parameter :: enum_kk_is_too_large_string = 365
integer, parameter :: enum_dnd_dtZ_diffnd_hyper3_meshZ_string = 366
integer, parameter :: enum_ZZdataZreactZout_string = 367
integer, parameter :: enum_ac_transformed_pencil_ttZ0_string = 368
integer, parameter :: enum_chemkin_string = 369
integer, parameter :: enum_get_reaction_rate_string = 370
integer, parameter :: enum_unit_system_must_be_cgsZ_string = 371
integer, parameter :: enum_oZo3_string = 372
integer, parameter :: enum_o1dZo_string = 373
integer, parameter :: enum_ohZcoZho2_string = 374
integer, parameter :: enum_2ho2Zh2o2_string = 375
integer, parameter :: enum_ohZhno3Zno3_string = 376
integer, parameter :: enum_calc_extra_react_string = 377
integer, parameter :: enum_no_such_reaction_nameZ__Z_string = 378
integer, parameter :: enum_Z_string = 379
integer, parameter :: enum_roux_string = 380
integer, parameter :: enum_nreactions_should_always_be_1_string = 381
integer, parameter :: enum_global_phi_must_be_given_as_input_string = 382
integer, parameter :: enum_o2_string = 383
integer, parameter :: enum_c3h8_string = 384
integer, parameter :: enum_o2_is_not_defined_string = 385
integer, parameter :: enum_c3h8_is_not_defined_string = 386
integer, parameter :: enum_i_io2fZ_i_c3h8Z_ichem_o2Z_ichem_c3h8Z_string = 387
integer, parameter :: enum_lo2Z_lc3h8Z_string = 388
integer, parameter :: enum_init_c3h8Zinit_o2Zmo2Zmc3h8Z_string = 389
integer, parameter :: enum_calc_pencils_chemistry_string = 390
integer, parameter :: enum_dchemistry_dtZ_kZdiff_kZ_string = 391
integer, parameter :: enum_clausius_string = 392
integer, parameter :: enum_cond_spec_sat_conc_string = 393
integer, parameter :: enum_no_such_iconc_sat_spec_string = 394
integer, parameter :: enum_kingery_string = 395
integer, parameter :: enum_cond_spec_nucl_rate_string = 396
integer, parameter :: enum_no_such_isurf_energy_string = 397
integer, parameter :: enum_oxtoby_string = 398
integer, parameter :: enum_calc_pencils_energy_string = 399
integer, parameter :: enum_lscale_to_cs2top_not_possible_string = 400
integer, parameter :: enum_lntt_string = 401
integer, parameter :: enum_denergy_dtZ_cs2_Z_string = 402
integer, parameter :: enum_dchemistry_dt_string = 403
integer, parameter :: enum_dchemistry_dtZ_solve_dchemistry_dt_string = 404
integer, parameter :: enum_fixed_swirl_string = 405
integer, parameter :: enum_cosxcosz_string = 406
integer, parameter :: enum_azsinx_string = 407
integer, parameter :: enum_aycosz_string = 408
integer, parameter :: enum_robertsflow_string = 409
integer, parameter :: enum_beltramiZz_string = 410
integer, parameter :: enum_shearingZ_sshearZsshear1Z_string = 411
integer, parameter :: enum_shearingZ_qshearZqshear0Z_string = 412
integer, parameter :: enum_finished_boundconds_z_string = 413
integer, parameter :: enum_accretor_string = 414
integer, parameter :: enum_default_string = 415
integer, parameter :: enum_calc_pencils_gravity_string = 416
integer, parameter :: enum_no_such_grav_type_string = 417
integer, parameter :: enum_denergy_dtZ_it_string = 418
integer, parameter :: enum_t_string = 419
integer, parameter :: enum_calc_heatcondZ_nans_in_thdiff_string = 420
integer, parameter :: enum_calc_heatcondZ_mZnZyZmZZzZnZZ_string = 421
integer, parameter :: enum_nans_in_thdiff_string = 422
integer, parameter :: enum_calc_heatcond_kramersZ_mZnZyZmZZzZnZZ_string = 423
integer, parameter :: enum_calc_heatcond_kramers_string = 424
integer, parameter :: enum_calc_heatcond_smagorinsky_string = 425
integer, parameter :: enum_get_lnq_string = 426
integer, parameter :: enum_tabulated_values_in_lntt_are_invalid_string = 427
integer, parameter :: enum_too_few_tabulated_values_in_lntt_string = 428
integer, parameter :: enum_sum_mn_string = 429
integer, parameter :: enum_not_implemented_for_cylindrical_string = 430
integer, parameter :: enum_ZtvartZdat_string = 431
integer, parameter :: enum_unknown_string = 432
integer, parameter :: enum_append_string = 433
integer, parameter :: enum_Z4f14Z7Z_string = 434
integer, parameter :: enum_standard_string = 435
integer, parameter :: enum_standard2_string = 436
integer, parameter :: enum_logZswitchZon_string = 437
integer, parameter :: enum_linearZsigma_string = 438
integer, parameter :: enum_eta_table_string = 439
integer, parameter :: enum_magnetic_after_boundary_string = 440
integer, parameter :: enum_eta_table_not_yet_completed_string = 441
integer, parameter :: enum_meanZfield_string = 442
integer, parameter :: enum_lrho_chi_string = 443
integer, parameter :: enum_initialize_magnetic_string = 444
integer, parameter :: enum_e2m_all_string = 445
integer, parameter :: enum_b2m_all_string = 446
integer, parameter :: enum_meanZfieldZlocal_string = 447
integer, parameter :: enum_electric_field_must_be_computed_string = 448
integer, parameter :: enum_calc_pencils_magneticZ_advec_va2__Z_string = 449
integer, parameter :: enum_thomson_string = 450
integer, parameter :: enum_fatal_error_wZo_force_string = 451
integer, parameter :: enum_lambdaZconstant_string = 452
integer, parameter :: enum_daa_dtZ_advec_hall_Z_string = 453
integer, parameter :: enum_tau_jj_must_be_finite_and_positive_string = 454
integer, parameter :: enum_forcingZ_add_continuous_forcing_string = 455
integer, parameter :: enum_fyZconst_string = 456
integer, parameter :: enum_fzZconst_string = 457
integer, parameter :: enum_abc_string = 458
integer, parameter :: enum_schur_nonhelical_string = 459
integer, parameter :: enum_schur_helical_string = 460
integer, parameter :: enum_abctdep_string = 461
integer, parameter :: enum_aka_string = 462
integer, parameter :: enum_grav_z_string = 463
integer, parameter :: enum_uniform_vorticity_string = 464
integer, parameter :: enum_kolmogorovflowZx_string = 465
integer, parameter :: enum_kolmogorovflowZz_string = 466
integer, parameter :: enum_nocos_string = 467
integer, parameter :: enum_straining_string = 468
integer, parameter :: enum_strainingexcact_string = 469
integer, parameter :: enum_forcing_cont_string = 470
integer, parameter :: enum_shock_string = 471
integer, parameter :: enum_hyper_string = 472
integer, parameter :: enum_getnu_string = 473
integer, parameter :: enum_some_viscosity_string = 474
integer, parameter :: enum_robertsflowii_string = 475
integer, parameter :: enum_robertsflowmask_string = 476
integer, parameter :: enum_robertsflow2d_string = 477
integer, parameter :: enum_robertsflow_exact_string = 478
integer, parameter :: enum_robertsflowZzdep_string = 479
integer, parameter :: enum_zZdependent_roberts_flowZ_eps_fcontZ_string = 480
integer, parameter :: enum_elevatorZflow_string = 481
integer, parameter :: enum_zZdependent_elevatorZflowZ_eps_fcontZ_string = 482
integer, parameter :: enum_robertsZforZssd_string = 483
integer, parameter :: enum_sinx_string = 484
integer, parameter :: enum_Z0Z0ZcosxZ_string = 485
integer, parameter :: enum_Z0Z0ZcosxcosyZ_string = 486
integer, parameter :: enum_bZZ0Z0ZcosxcosyZ_string = 487
integer, parameter :: enum_Z0ZxZ0Z_string = 488
integer, parameter :: enum_Z0ZsinxsintZ0Z_string = 489
integer, parameter :: enum_Z0ZsinxZ0Z_string = 490
integer, parameter :: enum_Z0ZcosxZcoszZ0Z_string = 491
integer, parameter :: enum_Z0ZsinxZexpZZzZ2ZZ0Z_string = 492
integer, parameter :: enum_Z0Zaycont_zZ0Z_string = 493
integer, parameter :: enum_ZsinzZcoszZ0Z_string = 494
integer, parameter :: enum_tg_string = 495
integer, parameter :: enum_tgZrandomZnonhel_string = 496
integer, parameter :: enum_tgZrandomZhel_string = 497
integer, parameter :: enum_cosxZcosyZcosz_string = 498
integer, parameter :: enum_gp_string = 499
integer, parameter :: enum_gallowayZproctorZ92_string = 500
integer, parameter :: enum_gp_tc13_string = 501
integer, parameter :: enum_gp_tc13_yzx_string = 502
integer, parameter :: enum_mbi_emf_string = 503
integer, parameter :: enum_j0_k1x_string = 504
integer, parameter :: enum_fluxring_cylindrical_string = 505
integer, parameter :: enum_counter_centrifugal_string = 506
integer, parameter :: enum_vortex_string = 507
integer, parameter :: enum_blob_string = 508
integer, parameter :: enum_zblob_string = 509
integer, parameter :: enum_vert_field_blob_string = 510
integer, parameter :: enum_ks_string = 511
integer, parameter :: enum_expZZx2Zy2Z_string = 512
integer, parameter :: enum_xz_string = 513
integer, parameter :: enum_1ZZ4Z3ZZ1ZrZ2Z4ZZrZ2_string = 514
integer, parameter :: enum_tidal_string = 515
integer, parameter :: enum_from_file_string = 516
integer, parameter :: enum_no_such_iforcing_contZ__string = 517
integer, parameter :: enum_axelZ_should_not_be_here_ZetaZ_ZZZ__string = 518
integer, parameter :: enum_superconformal_string = 519
integer, parameter :: enum_axel2Z_should_not_be_here_ZetaZ_ZZZ__string = 520

endmodule Cparam
