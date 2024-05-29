! $Id$
!
!  This modules contains the routines for SNe-driven ISM simulations.
!  Replaces old module relabelled as interstellar_old Jul 15 FAG.
!
!***************** AUTOMATIC CPARAM.INC GENERATION ***************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linterstellar = .true.
!
! MAUX CONTRIBUTION 2
!!! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED heat; cool; heatcool
!
!*****************************************************************************
module Interstellar
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'interstellar.h'
  include 'record_types.h'
  include 'eos_params.h'
!
  interface input_persistent_interstellar
    module procedure input_persist_interstellar_id
    module procedure input_persist_interstellar
  endinterface
!
  type ExplosionSite
    real :: rho, lnrho, yH, lnTT, TT, ss, ee
  endtype
!
  type RemnantFeature
    real :: x, y, z, t    ! Time and location
    real :: MM, EE, CR    ! Mass, energy & CR injected
    real :: rhom          ! Local mean density and energy at explosion time
    real :: radius, dr    ! Injection radius and resolution
    real :: t_sedov, t_SF ! Sedov time shell formation time
  endtype
!
  type RemnantIndex
    integer :: l,m,n           ! Grid position
    integer :: iproc,ipx,ipy,ipz
    integer :: SN_type
    integer :: state
  endtype
!
  type SNRemnant
    type (RemnantFeature) :: feat
    type (RemnantIndex) :: indx
    type (ExplosionSite) :: site
  endtype
!
!  Add cluster type
!
  type ClusterFeature
    real :: x, y, z, t ! Location and end time
    real :: radius     ! OB radius
    real :: rhom0      ! Mean density at OB initiation
    real :: rhom       ! Mean density in OB currently
    real :: tnextOB    ! Sedov time shell formation time
  endtype
!
  type ClusterIndex
    integer :: l,m,n             ! Grid position
    integer :: iproc,ipx,ipy,ipz ! Processor location
    integer :: SN_number, iSN_remain ! Total SN in OB and SN remaining
    integer :: iOB               ! OB identifier
    integer :: state
  endtype
!
  type SNCluster
    type (ClusterFeature) :: feat
    type (ClusterIndex) :: indx
  endtype
!
!  Enumeration of Supernovae states.
!
  integer, parameter :: SNstate_invalid   = 0
  integer, parameter :: SNstate_waiting   = 1
  integer, parameter :: SNstate_exploding = 2
  integer, parameter :: SNstate_finished  = 3
!
!  Enumeration of Explosion Errors
!
  integer, parameter :: iEXPLOSION_OK           = 0
  integer, parameter :: iEXPLOSION_TOO_HOT      = 1
  integer, parameter :: iEXPLOSION_TOO_RARIFIED = 2
  integer, parameter :: iEXPLOSION_TOO_UNEVEN   = 3
!
!  Enumeration of OB cluster states.
!
  integer, parameter :: OBstate_invalid   = 0
  integer, parameter :: OBstate_waiting   = 1
  integer, parameter :: OBstate_active    = 2
  integer, parameter :: OBstate_finished  = 3
!
!  04-sep-09/fred: amended xsi_sedov
!  ref Dyson & Williams Ch7 value = (25/3/pi)**(1/5)=1.215440704 for gamma=5/3
!  2nd ref Ostriker & McKee 1988 Rev.Mod.Phys 60,1 1.15166956^5=2.206
!  Est'd value for similarity variable at shock
!
!  real :: xsi_sedov=1.215440704
!  real :: xsi_sedov=1.15166956, mu
  real :: xsi_sedov=2.026, mu, SFt_norm, SFr_norm, sedov_norm, kfrac_norm
!
!  'Current' SN Explosion site parameters
!
  integer, parameter :: mSNR = 100, mOB = 20
  integer :: nSNR = 0, nOB = 0
  type (SNRemnant), dimension(mSNR) :: SNRs
  type (SNCluster), dimension(mOB) :: OBs
  integer, dimension(mSNR) :: SNR_index
  integer, dimension(mOB) :: OB_index
  integer, parameter :: npreSN = 5
  integer, dimension(4,npreSN) :: preSN
!
!  Squared distance to the SNe site along the current pencil
!  Outward normal vector from SNe site along the current pencil
!
  real, dimension(nx) :: dr2_SN, dr2_OB
  !$omp threadprivate(dr2_SN, dr2_OB)
!
!  Allocate time of next SNI/II and intervals until next
!
  real :: t_next_SNI=0.0, t_next_SNII=0.0
  real :: x_cluster=0.0, y_cluster=0.0, z_cluster=0.0, t_cluster=0.0
  real, dimension(2) :: t_interval=impossible
  real :: t_interval_OB=impossible
  real :: zdisk, maxrho !varying location of centre of mass of the disk
  logical :: lfirst_zdisk
!
  logical :: lfirst_warning=.true.
!  normalisation factors for 1-d, 2-d, and 3-d profiles like exp(-r^6)
!  ( 1d: 2    int_0^infty exp(-(r/a)^6)     dr) / a
!    2d: 2 pi int_0^infty exp(-(r/a)^6) r   dr) / a^2
!    3d: 4 pi int_0^infty exp(-(r/a)^6) r^2 dr) / a^3 )
!  ( cf. 3.128289613 -- from where ?!? )
!  NB: 1d and 2d results just from numerical integration -- calculate
!      exact integrals at some point...
!  3-D  was 3.71213666 but replaced with Maple result....
!
  real, parameter, dimension(3) :: &
      cnorm_gaussian_SN    = (/ 0.8862269254527579, pi,                5.568327996831708  /)
  real, parameter, dimension(3) :: &
      cnorm_gaussian2_SN   = (/ 0.9064024770554771, 2.784163998415854, 3.849760110050832  /)
  real, parameter, dimension(3) :: &
      cnorm_SN             = (/ 0.9277193336300392, 2.805377873352155, 3.712218664554472  /)
  real, parameter, dimension(3) :: &
      cnorm_para_SN        = (/ fourthird,          1.570796326794897, 1.6755161          /)
  real, parameter, dimension(3) :: &
      cnorm_quar_SN        = (/  0.,                2.0943951,         0.                 /)
!  kinetic energy with lmass_SN=F
  real, parameter, dimension(3) :: &
      vnormEj_gaussian_SN  = (/ 0.5116633539732443, 1.047197551196598, 1.0716252226356386 /)
!  kinetic energy addition lmass_SN=T
  real, parameter, dimension(3) :: &
      vnorm_gaussian_SN    = (/ 0.6266570686577501, 1.570796326794897, 1.9687012432153024 /)
  real, parameter, dimension(3) :: &
      vnormEj_gaussian2_SN = (/ 0.6887169476297503, 1.607437833953458, 1.6888564123130090 /)
!  kinetic energy addition lmass_SN=T
  real, parameter, dimension(3) :: &
      vnorm_gaussian2_SN   = (/ 0.7621905937330379, 1.968701243215302, 2.2890810569630537 /)
  real, parameter, dimension(3) :: &
      vnormEj_SN           = (/ 0.7724962826996008, 1.945140377302524, 2.1432504452712773 /)
!  kinetic energy addition lmass_SN=T
  real, parameter, dimension(3) :: &
      vnorm_SN             = (/ 0.8265039651250117, 2.226629893663761, 2.624934990953737  /)
!
!  cp1=1/cp used to convert TT (and ss) into interstellar code units
!  (useful, as many conditions conveniently expressed in terms of TT)
!  code units based on:
!    [length]  = 1kpc  = 3.09 10^21 cm
!    [time]    = 1Gyr  = 3.15 10^16 s             !no, on [u]=1km/s...
!    [rho]     =       = 1.00 10^-24 g/cm^3
!  Lambdaunits converts coolH into interstellar code units.
!
  real :: unit_Lambda, unit_Gamma
!
!  SNe placement limitations (for code stability)
!  Minimum resulting central temperature of a SN explosion.
!  If this is not reached then consider moving mass to achieve this.
!  10-aug-10/fred:
!  As per joung et al apj653 2005 min temp 1E6 to avoid excess radiative
!  energy losses in early stages.
!
  real, parameter :: rho_SN_min_cgs=1E-28,rho_SN_max_cgs=2.364E-23
!  fred: max rho intended to avoid explosion sites that are difficult to
!  resolve, but can lead to persistent high density structures that cannot be
!  destroyed by SN, so may be better to allow unrestricted
  real, parameter :: TT_SN_min_cgs=1., TT_SN_max_cgs=5E6
  real :: rho_SN_min=impossible, rho_SN_max=impossible
  real :: TT_SN_min=impossible, TT_SN_max=impossible
  real :: SN_rho_ratio=1e4, SN_TT_ratio=1.0e1
!
!  Enumeration of Supernovae types.
!
  integer, parameter :: SNI=1, SNII=2
!
!  SNI per (x,y)-area explosion rate
!
  real(KIND=rkind8), parameter :: SNI_area_rate_cgs=1.330982784D-56
  real(KIND=rkind8), parameter :: OB_area_rate_cgs=1.576417151D-57
  real :: SNI_area_rate=impossible, SNII_area_rate=impossible
  real :: OB_area_rate=impossible
  real :: SNI_factor=1.0, SNII_factor=1.0
!
!  SNII rate=5.E-12 mass(H1+HII)/solar_mass
!  van den Bergh/Tammann Annu. Rev Astron. Astrophys. 1991 29:363-407
!  SNI rate=4.7E-14/solar_mass + 0.35 x SNII rate
!  Mannucci et al A&A 433, 807-814 (2005)
!
  real, parameter :: SNII_mass_rate_cgs=1.584434515E-19
  real, parameter :: SNI_mass_rate_cgs=1.489368444E-21
  real :: SNII_mass_rate, SNI_mass_rate, Lxyzmin
  real :: SN_interval_rhom=impossible, SN_interval_rhom_cgs=2.8e-25
  logical :: lSN_mass_rate=.false., lscale_SN_interval=.false.
  real :: iSNdx=4.
!
!  Some useful constants
!
  real, parameter :: kpc_cgs=3.086E+21, pc_cgs=3.086E+18 ! [cm]
  real, parameter :: yr_cgs=3.155692E7, kyr_cgs=3.155692E10, &
          Myr_cgs=3.155692E13 ! [s]
  real, parameter :: solar_mass_cgs=1.989E33 ! [g]
  real :: solar_mass
!
!  Scale heights for SNI/II with Gaussian z distributions
!
  real, parameter :: h_SNI_cgs=1.00295E21, h_SNII_cgs=2.7774E20
  real :: h_SNI=impossible, h_SNII=impossible
  logical :: lh_SNII_adjust=.false.
!
!  Self regulating SNII explosion coefficients
!
  real, parameter :: cloud_rho_cgs=1.67262158E-24, cloud_TT_cgs=4000.
  real, parameter :: cloud_tau_cgs=2.E7 * yr_cgs, minTT_cgs = 0.75E2
  real, parameter :: mass_SN_progenitor_cgs=10.*solar_mass_cgs
  real, parameter :: frac_converted=0.02, frac_heavy=0.10
!  real, parameter :: tosolarMkpc3=1.483E7
  real :: cloud_rho=impossible, cloud_TT=impossible, lncloud_TT
  !!$omp target if(loffload) enter data map(cloud_rho,lncloud_TT)
  real :: cloud_tau=impossible
  real :: mass_SN_progenitor=impossible
!
!  Total SNe energy
!
  real(KIND=rkind8), parameter :: ampl_SN_cgs=1D51
  real :: frac_ecr=0.0, frac_eth=1.0, frac_kin=0.0, kin_max=0.075
  real :: ampl_SN=impossible, campl_SN=0.0, eampl_SN=0.0, kampl_SN=0.0, &
          kperp=0.05, kpara=0.025
!
!  SNe composition
!
  logical :: lSN_eth=.true., lSN_ecr=.false., lSN_mass=.false., &
             lSN_velocity=.false., lSN_fcr=.false., lSN_autofrackin=.true., &
             lSN_momentum=.true., lSN_coolingmass=.false.
!
!  Total mass added by a SNe
!
  real, parameter :: mass_SN_cgs=10.*solar_mass_cgs
  real :: mass_SN=impossible
!
!  Size of SN insertion site (energy and mass) and shell in mass movement
!
  real :: sigma_SN, sigma_SN1
  real, parameter :: width_SN_cgs=6.172E19
  real :: energy_width_ratio=1.
  real :: mass_width_ratio=1.
  real :: velocity_width_ratio=1.
  real :: outer_shell_proportion = 1.2
  real :: inner_shell_proportion = 1.
  real :: width_SN=impossible
  real :: energy_Nsigma=impossible, energy_Nsigma2
!
!  Set SN by predefined list of time and coordinates for direct comparison
!  between models
!
  logical :: lSN_list=.false.
  real, dimension(:,:), allocatable :: SN_list
  integer, dimension(:), allocatable :: SN_type
  integer :: type_list, nlist, SNfirst
!
!  Parameters for 'averaged'-SN heating
!
  real :: r_SNI_yrkpc2=4.E-6, r_SNII_yrkpc2=3.E-5
  real :: r_SNI, r_SNII
  real :: average_SNI_heating=impossible, average_SNII_heating=impossible
!
!  Limit placed of minimum density resulting from cavity creation and
!  parameters for thermal_hse(hydrostatic equilibrium) assuming RBNr
!
  real, parameter :: rho_min_cgs=1.E-34, rho0ts_cgs=3.5E-24, T_init_cgs=1.E3
  real :: rho0ts=impossible, T_init=impossible, rho_min=impossible
!
!  Time of most recent SNII event
!
  real :: last_SN_t=0.
!
!  Time to wait before allowing SN to begin firing
!
  real :: t_settle=0.
!
!  Initial dist'n of explosions
!
  character (len=labellen), dimension(ninit) :: initinterstellar='nothing'
!
!  Number of randomly placed SNe to generate to mix the initial medium
!
  integer :: initial_SNI=0
!
!  Parameters for UV heating of Wolfire et al.
!
  real, parameter :: rhoUV_cgs=0.1
  real, parameter :: GammaUV_cgs=0.0147
  real, parameter :: TUV_cgs=7000., T0UV_cgs=20000., cUV_cgs=5.E-4
  real :: GammaUV=impossible, T0UV=impossible, cUV=impossible
!
!  04-jan-10/fred:
!  Amended cool dim from 7 to 11 to accomodate WSW dimension.
!  Appended null last term to all arrays for RBN and SS cooling
!
  integer, parameter :: len_cool=11
  real(KIND=rkind8), dimension(len_cool) :: lncoolH, coolH_cgs
  real, dimension(len_cool) :: coolT_cgs, coolB, lncoolT
  integer :: ncool
!
!  TT & z-dependent uv-heating profile
!
  real, dimension(mz) :: heat_z
  logical :: lthermal_hse=.false., lheatz_min=.true.
!
  real :: coolingfunction_scalefactor=1.
  real :: heatingfunction_scalefactor=1.
  real, dimension(2) :: heatingfunction_scale
  real :: heatingfunction_fadefactor=1.
!
  real :: heating_rate = 0.015
  real :: heating_rate_code = impossible
!
  real :: heatcool_shock_cutoff_rate = 0.
  real :: heatcool_shock_cutoff_rate1 = 0.0
!
!  Set .true. to smoothly turn off the heating and cooling where the
!  shock_profile is > 0 if  heatcool_shock_cutoff_rate/=0
!  Fred: revealed too high timestep for large df(u) source of code instability
!  rather than cooling timestep or shock instability, so this device no longer
!  recommended.
!  2009 simulations used mu=0.62 for partially ionized medium. Direct comparison
!  using alternate mu with ideal gas can be used with lcooling_revert=.true., now
!  somewhat redundant as more recent simulations more robust.
!
  logical :: lheatcool_shock_cutoff = .false., lcooling_revert=.false.
!
!  SN type flags
!
  logical :: lSNI=.true., lSNII=.true.
!
!  Cooling & heating flags
!
  logical :: laverage_SNI_heating = .false.
  logical :: laverage_SNII_heating = .false.
!
!  Remnant location flags
!
  logical :: luniform_zdist_SNI = .false., lSNII_gaussian=.true.
  logical :: lOB_cluster = .false. ! SN clustering
  real    :: p_OB=0.7 ! Probability that an SN is in a cluster
  real    :: SN_clustering_radius=impossible
  real    :: SN_clustering_time=impossible
  real    :: SN_clustering_radius_cgs=1.25e20 ! cm (~40 pc)
  real    :: SN_clustering_time_cgs=1.6e14 ! cm (~5 Myr)
!
!  Adjust SNR%feat%radius inversely with density
!
  logical :: lSN_scale_rad=.false.
  real :: N_mass=250.0, eps_mass=0.05, rfactor_SN=5.0, Nsol_added=30.
!
!  Requested SNe location (used for test SN)
!
  real :: center_SN_x = impossible
  real :: center_SN_y = impossible
  real :: center_SN_z = impossible
!
!  Cooling time diagnostic
!
  integer :: idiag_taucmin=0  ! DIAG_DOC: $\min(\tau_{\rm cool})$
  integer :: idiag_Hmax_ism=0 ! DIAG_DOC: $\max(\Gamma-\rho\Lambda)$
  integer :: idiag_Lamm=0     ! DIAG_DOC: $\left<\Lambda\right>$
  integer :: idiag_nrhom=0    ! DIAG_DOC: TBC
  integer :: idiag_rhoLm=0    ! DIAG_DOC: $\left<\rho\Lambda\right>$
  integer :: idiag_Gamm=0     ! DIAG_DOC: $\left<\Gamma\right>$
  integer :: idiag_rhoHCmz=0  ! XYAVG_DOC: $\left<\rho\Gamma-\rho^2\Lambda\right>_{xy}$
!
!  Heating function, cooling function and mass movement
!  method selection.
!
  character (len=labellen) :: cooling_select  = 'WSW'
  character (len=labellen) :: heating_select  = 'wolfire'
  character (len=labellen) :: thermal_profile = 'gaussian'
  character (len=labellen) :: velocity_profile= 'gaussian'
  character (len=labellen) :: mass_profile    = 'gaussian'
!
!  Variables required for returning mass to disk given no inflow
!  boundary condition used in addmassflux
!
  real :: boldmass=0.0, old_rhom=1.0
  logical :: ladd_massflux = .false.
!  switches required to override the persistent values when continuing a run
  logical :: l_persist_overwrite_lSNI=.false., l_persist_overwrite_lSNII=.false.
  logical :: l_persist_overwrite_tSNI=.false., l_persist_overwrite_tSNII=.false.
  logical :: l_persist_overwrite_tcluster=.false., l_persist_overwrite_xcluster=.false.
  logical :: l_persist_overwrite_ycluster=.false., l_persist_overwrite_zcluster=.false.
  logical :: lreset_ism_seed=.false.
  integer :: seed_reset=1963
!
!  Gravity constansts - required for pre-2015 vertical heating profile
!
  real :: a_S, a_D
  real, parameter :: a_S_cgs=4.4e-9, a_D_cgs=1.7e-9
  real :: z_S, z_D, H_z
  real, parameter :: z_S_cgs=6.172e20, z_D_cgs=3.086e21, H_z_cgs=9.258E20
!
  logical :: loffload=.false.
!
!  start parameters
!
  namelist /interstellar_init_pars/ &
      initinterstellar, initial_SNI, h_SNI, h_SNII, lSNII, lSNI, &
      lSN_scale_rad, ampl_SN, mass_SN, width_SN, &
      mass_width_ratio, energy_width_ratio, velocity_width_ratio, &
      t_next_SNI, t_next_SNII, center_SN_x, center_SN_y, center_SN_z, &
      lSN_eth, lSN_ecr, lSN_fcr, lSN_mass, lSN_list, &
      frac_ecr, frac_kin, thermal_profile, velocity_profile, mass_profile, &
      luniform_zdist_SNI, inner_shell_proportion, outer_shell_proportion, &
      SNI_factor, SNII_factor, lSN_autofrackin, kin_max, &
      cooling_select, heating_select, heating_rate, GammaUV, rho0ts, &
      T_init, TT_SN_max, rho_SN_min, N_mass, lSNII_gaussian, rho_SN_max, &
      lthermal_hse, lheatz_min, kperp, kpara, average_SNII_heating, &
      average_SNI_heating, SN_rho_ratio, SN_TT_ratio, &
      eps_mass, rfactor_SN, &
      energy_Nsigma, lSN_momentum, lSN_coolingmass, Nsol_added
!
!  run parameters
!
  namelist /interstellar_run_pars/ &
      ampl_SN, mass_SN, t_next_SNI, t_next_SNII, lSN_list, &
      mass_width_ratio, energy_width_ratio, velocity_width_ratio, &
      lSN_eth, lSN_ecr, lSN_fcr, lSN_mass, width_SN, lSNI, lSNII, &
      luniform_zdist_SNI, SNI_area_rate, SNII_area_rate, &
      SNI_factor, SNII_factor, lSN_autofrackin, kin_max, &
      frac_ecr, frac_kin, thermal_profile,velocity_profile, mass_profile, &
      h_SNI, h_SNII, TT_SN_min, lSN_scale_rad, lh_SNII_adjust, &
      mass_SN_progenitor, cloud_tau, cloud_rho, cloud_TT, &
      laverage_SNI_heating, laverage_SNII_heating, coolingfunction_scalefactor, &
      heatingfunction_scalefactor, heatingfunction_fadefactor, t_settle, &
      center_SN_x, center_SN_y, center_SN_z, rho_SN_min, TT_SN_max, &
      cooling_select, heating_select, heating_rate, GammaUV, &
      heatcool_shock_cutoff_rate, ladd_massflux, lcooling_revert, &
      N_mass, T_init, rho0ts, &
      lSNII_gaussian, rho_SN_max, lSN_mass_rate, lthermal_hse, lheatz_min, &
      p_OB, SN_clustering_time, SN_clustering_radius, lOB_cluster, kperp, &
      average_SNII_heating, average_SNI_heating, seed_reset, &
      l_persist_overwrite_lSNI, l_persist_overwrite_lSNII, &
      l_persist_overwrite_tSNI, l_persist_overwrite_tSNII, &
      l_persist_overwrite_tcluster, l_persist_overwrite_xcluster, &
      l_persist_overwrite_ycluster, l_persist_overwrite_zcluster, &
      lreset_ism_seed, SN_rho_ratio, SN_TT_ratio, eps_mass, &
      lscale_SN_interval, SN_interval_rhom, rfactor_SN, iSNdx, &
      energy_Nsigma, lSN_momentum, lSN_coolingmass, Nsol_added, loffload
!
  real :: gamma
!
  contains
!
!***********************************************************************
    subroutine register_interstellar()
!
!  19-nov-02/tony: coded
!
      use FArrayManager
!
!      call farray_register_auxiliary('netheat',inetheat,communicated=.true.)
      call farray_register_auxiliary('netheat',inetheat)
      call farray_register_auxiliary('cooling',icooling)
!
!  identify version number
!
      if (lroot) call svn_id( &
              "$Id$")
!
!  Invalidate all SNRs
!
      nSNR=0
      SNRs%indx%state=SNstate_invalid
!
!  Writing files for use with IDL
!
!      if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=',netheat $'
!      if (naux+naux_com == maux+maux_com) aux_var(aux_count)=',netheat'
!      aux_count=aux_count+1
!      if (lroot) write(15,*) 'netheat = fltarr(mx,my,mz)*one'
!
    endsubroutine register_interstellar
!***********************************************************************
    subroutine initialize_interstellar(f)
!
!  Perform any post-parameter-read initialization eg. set derived
!  parameters
!
!  24-nov-02/tony: coded
!
!  10-feb-23/fred: timing and location of series of SNe shall be read
!                  from formated file ./sn_series.in
!                  containing same format as data/sn_series.dat typically
!                  from previous run, if lSN_list.
!
      use EquationOfState, only: getmu, get_gamma_etc
      use General, only: random_seed_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i, int1_list, stat
      integer, dimension(4) :: int4_list
      real, dimension(13) :: real13_list
      real :: t_list,x_list,y_list,z_list
      logical :: exist
!
      f(:,:,:,icooling)=0.0
      f(:,:,:,inetheat)=0.0
!
      if (lroot) print"(1x,'initialize_interstellar: t_next_SNI =',e15.8)",t_next_SNI
!
      if (lroot.and.luniform_zdist_SNI) &
          print*,'initialize_interstellar: using UNIFORM z-distribution of SNI'

      lSN_ecr=lcosmicray.and.lSN_ecr
      lSN_fcr=lcosmicrayflux.and.lSN_fcr.and.lSN_ecr
!
      call get_gamma_etc(gamma)
      if (leos_idealgas) then
        call getmu(f,mu)
      else
        mu = 0.62
        if (lroot) print*, 'initialize_interstellar: mu not set by eos:',mu
      endif
!
      if (unit_system=='cgs') then
!
!  this Lambda as such enters as n^2*Lambda(T) on the rhs of the
!  energy equation per unit volume (n Gamma(T))
!
        unit_Lambda = unit_velocity**2 / unit_density / unit_time
        unit_Gamma  = unit_velocity**3 / unit_length
!
        T0UV=T0UV_cgs / unit_temperature
        cUV=cUV_cgs * unit_temperature
        if (GammaUV==impossible) GammaUV=GammaUV_cgs / unit_Gamma
!
        if (rho_SN_min==impossible) rho_SN_min=rho_SN_min_cgs / unit_density
        if (rho_SN_max==impossible) rho_SN_max=rho_SN_max_cgs / unit_density
        if (TT_SN_max==impossible) TT_SN_max=TT_SN_max_cgs/unit_temperature
        if (TT_SN_min==impossible) TT_SN_min=TT_SN_min_cgs / unit_temperature
        if (OB_area_rate==impossible) OB_area_rate=OB_area_rate_cgs * unit_length**2 * unit_time
        if (SNI_area_rate==impossible) SNI_area_rate=SNI_area_rate_cgs * unit_length**2 * unit_time
        if (SNII_area_rate==impossible) &
            SNII_area_rate=7.5*SNI_area_rate_cgs * unit_length**2 * unit_time
        if (h_SNI==impossible) h_SNI=h_SNI_cgs / unit_length
        if (h_SNII==impossible) h_SNII=h_SNII_cgs / unit_length
        SNII_mass_rate=SNII_mass_rate_cgs*unit_time
        SNI_mass_rate=SNI_mass_rate_cgs*unit_time
        solar_mass=solar_mass_cgs / unit_mass
        if (lroot.and.ip==1963) print"(1x,'initialize_interstellar: solar_mass =',e11.4)",solar_mass
        r_SNI =r_SNI_yrkpc2  * (unit_time/yr_cgs) * (unit_length/kpc_cgs)**2
        r_SNII=r_SNII_yrkpc2 * (unit_time/yr_cgs) * (unit_length/kpc_cgs)**2
!
!  set SN parameters self-consistently
!
        if (ampl_SN==impossible) ampl_SN=ampl_SN_cgs / unit_energy
!  dimensional norm for Sedov-Taylor relations
        sedov_norm=unit_density/1e-24*ampl_SN_cgs/unit_energy
!  parameters for energy losses prior to SN initialisation
!  ref Kim & Ostriker 2015 Eq 7 dimensional norm shell formation time
!  ref Simpson et al. 2015 Eq 17
        SFt_norm = 26.5*kyr_cgs/unit_time*&
                   (1.4*m_H_cgs/unit_density)**(4./7)*&
                   (unit_energy/ampl_SN_cgs)**(3./14)
!  ref Kim & Ostriker 2015 Eq 8 dimensional norm shell formation radius
!  ref Simpson et al. 2015 Eq 18
        SFr_norm = 18.5*pc_cgs/unit_length*&
                   (unit_energy/ampl_SN_cgs)**(2./7)*&
                   (1.4*m_H_cgs/unit_density)**(3./7)
!  ref Simpson 15 Eq 16 dimensional norm kinetic energy fraction
        kfrac_norm=3.97e-6*mu/1.4/m_H_cgs*unit_density*& !RPDS
                   ampl_SN_cgs/unit_energy*(unit_length/pc_cgs)**5*&
                   (kyr_cgs/unit_time)**2
        if (lroot.and.lSN_autofrackin.and.ip==1963) then
            print "(1x,'initialize_interstellar:   SFt_norm =',e11.4)",SFt_norm
            print "(1x,'initialize_interstellar:   SFr_norm =',e11.4)",SFr_norm
            print "(1x,'initialize_interstellar: sedov_norm =',e11.4)",sedov_norm
            print "(1x,'initialize_interstellar: kfrac_norm =',e11.4)",kfrac_norm
        endif
        if (lSN_coolingmass) then
          if (.not.lSN_eth) then
            if (lroot) print*,'initialize_interstellar: lSN_coolingmass false if not lSN_eth'
            lSN_coolingmass = .false.
          else
            lSN_mass = .false.
          endif
        endif
        if (lSN_autofrackin) lSN_velocity = .true.
        if (lSN_eth) then
          if (lSN_ecr) then
            if (frac_ecr==0) frac_ecr=0.005
            campl_SN=frac_ecr*ampl_SN
          endif
          frac_eth=1.-frac_ecr-frac_kin
          if (frac_eth<0.) call fatal_error('initialize_interstellar','energy fractions must sum to 1')
          kampl_SN=frac_kin*ampl_SN
          eampl_SN=frac_eth*ampl_SN
          if (frac_kin >0.0) lSN_velocity = .true.
        else
          if (lSN_ecr) then
            if (frac_ecr==0) frac_ecr=0.1
            campl_SN=frac_ecr*ampl_SN
          endif
          if (frac_kin >0.0) lSN_velocity = .true.
          kampl_SN=frac_kin*ampl_SN
          if (frac_kin+frac_ecr>1) &
            call fatal_error('initialize_interstellar','energy fractions not to be > 1')
        endif
        if (lroot) print &
            "(1x,'initialize_interstellar: eampl_SN, kampl_SN = ',2e11.4)", eampl_SN, kampl_SN
        if (energy_Nsigma==impossible) then
          if (thermal_profile=="gaussian3") then
            energy_Nsigma=1.25
          elseif (thermal_profile=="gaussian2") then
            energy_Nsigma=1.75
          elseif (thermal_profile=="gaussian") then
            energy_Nsigma=2.25
          else
            energy_Nsigma=1.5
          endif
        endif
        energy_Nsigma2=energy_Nsigma**2
        if (lSN_scale_rad) then
           if (thermal_profile=="gaussian")  &
               Lxyzmin=min(minval(Lxyz)/2, minval(Lxyz)/cnorm_gaussian_SN(dimensionality))
           if (thermal_profile=="gaussian2") &
               Lxyzmin=min(minval(Lxyz)/2, minval(Lxyz)/cnorm_gaussian2_SN(dimensionality))
           if (thermal_profile=="gaussian3") &
               Lxyzmin=min(minval(Lxyz)/2, minval(Lxyz)/cnorm_SN(dimensionality))
           if (lroot) print "(1x,'initialize_interstellar: Lxyzmin = ',e11.4)", Lxyzmin
        endif
        if (rho_min == impossible) rho_min=rho_min_cgs/unit_temperature
        if (T_init == impossible) T_init=T_init_cgs/unit_temperature
        if (rho0ts == impossible) rho0ts=rho0ts_cgs/unit_density
        if (cloud_rho==impossible) cloud_rho=cloud_rho_cgs / unit_density
        if (cloud_TT==impossible) cloud_TT=cloud_TT_cgs / unit_temperature
        lncloud_TT = log(cloud_TT)
        if (cloud_tau==impossible) cloud_tau=cloud_tau_cgs / unit_time
        if (mass_SN==impossible) mass_SN=mass_SN_cgs / unit_mass
        if (mass_SN_progenitor==impossible) &
            mass_SN_progenitor=mass_SN_progenitor_cgs / unit_mass
        if (width_SN==impossible) width_SN= &
            min(width_SN_cgs / real(unit_length),rfactor_SN*dxmin)
        if (SN_clustering_radius==impossible) &
            SN_clustering_radius=SN_clustering_radius_cgs / unit_length
        if (SN_clustering_time==impossible) &
            SN_clustering_time=SN_clustering_time_cgs / unit_time
        t_interval_OB   = 1./(OB_area_rate * Lx * Ly)
      else
        call not_implemented('initialize_interstellar','SI unit conversions')
      endif
!
      call select_cooling(cooling_select,lncoolT,lncoolH,coolB)
!
      if (lroot) print "(1x,'initialize_interstellar: unit_Lambda',e11.4)",unit_Lambda
      if (lroot) print "(1x,'initialize_interstellar: unit_Gamma',e11.4)",unit_Gamma
!
      heating_rate_code=heating_rate*real(unit_length/unit_velocity**3)
!
      if (heating_select == 'thermal-hs') call heat_interstellar(f,heat_z)
!
!  Cooling cutoff in shocks
!
      if (heatcool_shock_cutoff_rate/=0.) then
        lheatcool_shock_cutoff=.true.
        heatcool_shock_cutoff_rate1=1.0/heatcool_shock_cutoff_rate
      else
        lheatcool_shock_cutoff=.false.
      endif
!
!  Slopeyness used for tanh rounding profiles etc.
!
      sigma_SN=dxmax*3
      sigma_SN1=1./sigma_SN
!
      preSN=0
!
      if (SN_interval_rhom==impossible) SN_interval_rhom=SN_interval_rhom_cgs/unit_density
      t_interval(SNI )= 1./( SNI_factor *  SNI_area_rate * Lx * Ly)
      t_interval(SNII)= 1./(SNII_factor * SNII_area_rate * Lx * Ly)
      if (average_SNI_heating == impossible) average_SNI_heating = &
          r_SNI *ampl_SN/(sqrt(2*pi)*h_SNI*SN_interval_rhom)
      if (average_SNII_heating == impossible) average_SNII_heating = &
          r_SNII*ampl_SN/(sqrt(2*pi)*h_SNII*SN_interval_rhom)
      if (lroot.and.ip==1963) print &
          "(1x,'initialize_interstellar: t_interval_SNI, SNI rate =',2e11.4)", &
          t_interval(SNI),SNI_factor*SNI_area_rate
      if (laverage_SNI_heating) then
        if (lSNI.or.lSNII) then
          if (lroot.and.ip==1963) print &
              "(1x,'initialize_interstellar: average_SNI_heating =',e11.4)", &
              average_SNI_heating*sqrt(2*pi)*h_SNI*SN_interval_rhom* &
              t_interval(SNI)/(t_interval(SNI)+t*heatingfunction_fadefactor)* &
              heatingfunction_scalefactor
        else
          if (lroot.and.ip==1963) print &
              "(1x,'initialize_interstellar: average_SNI_heating =',e11.4)", &
              average_SNI_heating*sqrt(2*pi)*h_SNI*SN_interval_rhom* &
              heatingfunction_scalefactor
        endif
      else
        if (lroot.and.ip==1963) print*, 'initialize_interstellar: average_SNI_heating = 0'
      endif
      if (laverage_SNII_heating) then
        if (lSNI.or.lSNII) then
          if (lroot.and.ip==1963) print &
              "(1x,'initialize_interstellar: average_SNII_heating =',e11.4)", &
              average_SNII_heating*sqrt(2*pi)*h_SNII*SN_interval_rhom* &
              t_interval(SNII)/(t_interval(SNII)+t*heatingfunction_fadefactor)* &
              heatingfunction_scalefactor
        else
          if (lroot.and.ip==1963) print &
              "(1x,'initialize_interstellar: average_SNII_heating =',e11.4)", &
              average_SNII_heating*sqrt(2*pi)*h_SNII*SN_interval_rhom* &
              heatingfunction_scalefactor
        endif
      else
        if (lroot.and.ip==1963) print*,'initialize_interstellar: average_SNII_heating =0'
      endif
      if ((laverage_SNI_heating.or.laverage_SNII_heating).and..not.(lSNI.or.lSNII)) &
          heatingfunction_scale = heatingfunction_scalefactor
!
      if (lroot.and.ip==1963) then
        print*,'initialize_interstellar: nseed,seed',nseed,seed(1:nseed)
        print*,'initialize_interstellar: finished'
      endif
      if (ladd_massflux) call warning('initialize_interstellar','ladd_massflux now redundant,' // &
                                    'set lconserve_total_mass, total_mass in density_run_pars instead')
!
!  Fred: 06-Nov-17 added SN_rate column and changed site_mass to site_Nsol
!        Note changes may affect reading/meaning of pre-existing file contents
!
      if (lroot .and. lstart) then
        open(1,file=trim(datadir)//'/sn_series.dat',position='append')
        write(1,'("#",5A)')  &
            '---it-----------t------------itype-iproc---l-----m-----n-------x-',&
            '-----------y------------z-----------rho---------rhom---------',&
            '-TT-----------EE----------Ekin----------Ecr--------t_sedov---',&
            '---radius------site_Nsol----added_Nsol-------maxTT---', &
            '---t_interval----SN_rate--'
        close(1)
      endif
!
      if (lSN_list) then
        inquire(file='sn_series.in',exist=exist)
        if (exist) then
          open(33,file='sn_series.in')
        else
          call fatal_error('initialize_interstellar','error - no sn_series input file')
        endif
!
!  Read profiles.
!
        nlist=0
        read(33,*,iostat=stat)
        do while(1==1)
          read(33,*,iostat=stat) int1_list,t_list,type_list,int4_list,x_list,y_list,z_list,real13_list
          !Ignore int1_list=0, applying under START and t_list<t to avoid repeating previous explosions
          if (t<=t_list.and.int1_list/=0) nlist=nlist+1
          if (stat<0) exit
        enddo
        close(33)
        if (lroot) print*,"initialize_interstellar: nlist =",nlist
        if (nlist<2) call fatal_error('initialize_interstellar',&
            'sn_series.in list needs extending or set lSN_list=F to continue')
        if (allocated(SN_list) ) deallocate(SN_list)
        allocate(SN_list(4,nlist))
        if (allocated(SN_type) ) deallocate(SN_type)
        allocate(SN_type(  nlist))
!
        open(33,file='sn_series.in')
        read(33,*,iostat=stat)
        i=1
        do while(i<=nlist)
          read(33,*,iostat=stat) int1_list,t_list,type_list,int4_list,x_list,y_list,z_list,real13_list
          if (stat<0) exit
          !Ignore int1_list=0, applying under START and t_list<t to avoid repeating previous explosions
          if (t_list>=t.and.int1_list/=0) then
            SN_list(1,i)=t_list
            SN_list(2,i)=x_list
            SN_list(3,i)=y_list
            SN_list(4,i)=z_list
            SN_type(  i)=type_list
            if (lroot.and.i==1) print*,"i,int1_list,t_list,t",i,int1_list,t_list,t
            if (lroot.and.i>=nlist-10) print*,"i,int1_list,t_list,t",i,int1_list,t_list,t
            i=i+1
          endif
        enddo
        close(33)
        SNfirst=1
        t_next_SNI = SN_list(1,1)
        t_next_SNII = SN_list(1,1)
      endif

      if (leos_ionization) &
          call warning('initialize_interstellar','using T/e instead of cv1 '// &
          'for diagnostics. Not yet implemented cv1 for ionization')
 !
 !  Write unit_Lambda to pc_constants file
 !
      if (lroot) then
        print "(1x,'initialize_interstellar: t_next_SNI, t_next_SNII=',2e15.8)", &
              t_next_SNI, t_next_SNII
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,'(a,1pd26.16)') 'unit_Lambda=',unit_Lambda
        write (1,'(a,1pd26.16)') 'unit_Gamma=',unit_Gamma
        close (1)
      endif

    endsubroutine initialize_interstellar
!*****************************************************************************
    subroutine select_cooling(cooling_select,lncoolT,lncoolH,coolB)
!
!  Routine for selecting parameters for temperature dependent cooling
!  Lambda.
!
      character (len=labellen), intent(IN) :: cooling_select
      real, dimension (:), intent(OUT)  :: lncoolT, coolB
      real(KIND=rkind8), dimension (:), intent(OUT)  :: lncoolH
      real :: lnmu2
!
!  Scale rho^2 to gas number density^2 with mu^2
!
      if (lcooling_revert) then
        lnmu2 = 2*log(mu/0.62)
      else
        lnmu2 = 0.
      endif
!
!  Mara: Initialize cooling parameters according to selection
!  Default selection 'RBN' Rosen & Bregman (1993)
!  Alternative selection 'SS' Sanchez-Salcedo et al. (2002)
!  Turn off cooling: cooling_select='off'
!  cooling_select in interstellar_init_pars added
!
      if (cooling_select == 'RBN') then
        if (lroot) print*,'select_cooling: default RBN cooling fct'
        coolT_cgs = (/  100.0,       &
                        2000.0,      &
                        8000.0,      &
                        1.0E5,       &
                        4.0E7,       &
                        1.0E9,       &
                        tiny(0E0),   &
                        tiny(0E0),   &
                        tiny(0E0),   &
                        tiny(0E0),   &
                        tiny(0E0) /)
        coolH_cgs = (/  2.238751968D-32,  &
                        1.0012D-30,  &
                        4.6240D-36,  &
                        1.7800D-18,  &
                        3.2217D-27,  &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0) /) / ( m_p_cgs )**2
        coolB =     (/  2.0,         &
                        1.5,         &
                        2.867,       &
                       -0.65,        &
                        0.5,         &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.)  /)
        ncool = 5
!
!  04-jan-10/fred:
!  Reset above to original RBN parameters. Altered coolB(5) and coolT(5,6) in
!  RBNr to ensure continuity and increase cooling at high temperatures later in
!  diffuse remnant cores.
!  RBNr: new lower terms for smooth cooling below 300K
!  and extended range to 1E13 in SSrr to deal with temperature spiking
!
      else if (cooling_select == 'RBNr') then
        if (lroot) print*,'select_cooling: RBN cooling fct (revised)'
        coolT_cgs = (/  10.0,         &
                        2000.0,       &
                        8000.0,       &
                        1E5,          &
                        1E6,          &
                        1E17,         &
                        tiny(0E0),    &
                        tiny(0E0),    &
                        tiny(0E0),    &
                        tiny(0E0),    &
                        tiny(0E0) /)
        coolH_cgs = (/  2.2380D-32,   &
                        1.0012D-30,   &
                        4.6240D-36,   &
                        1.7783524D-18,&
                        2.238814D-25, &
                        tiny(0D0),    &
                        tiny(0D0),    &
                        tiny(0D0),    &
                        tiny(0D0),    &
                        tiny(0D0),    &
                        tiny(0D0) /)  / ( m_p_cgs )**2
        coolB =     (/  2.0,          &
                        1.5,          &
                        2.867,        &
                       -0.65,         &
                        0.5,          &
                        tiny(0.),     &
                        tiny(0.),     &
                        tiny(0.),     &
                        tiny(0.),     &
                        tiny(0.),     &
                        tiny(0.)  /)
        ncool=5
      else if (cooling_select == 'SS') then
!
!  These are the SS et al (2002) coefficients multiplied by m_proton**2
!
        coolT_cgs = (/  10.0,        &
                        141.0,       &
                        313.0,       &
                        6102.0,      &
                        1.0E5,       &
                        1.0E17,      &
                        tiny(0E0),   &
                        tiny(0E0),   &
                        tiny(0E0),   &
                        tiny(0E0),   &
                        tiny(0E0) /)
        coolH_cgs = (/  3.42D16,     &
                        9.10D18,     &
                        1.11D20,     &
                        2.00D8,      &
                        7.962D29,    &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0) /)
        coolB =     (/  2.12,        &
                        1.0,         &
                        0.56,        &
                        3.67,        &
                        -0.65,       &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.)  /)
        ncool=5
!
      else if (cooling_select == 'SSr') then
        if (lroot) print*,'select_cooling: revised SS cooling fct'
        coolT_cgs = (/  10.0,        &
                        141.0,       &
                        313.0,       &
                        6102.0,      &
                        1E5,         &
                        1E9,         &
                        1E17,        &
                        tiny(0E0),   &
                        tiny(0E0),   &
                        tiny(0E0),   &
                        tiny(0E0) /)
        coolH_cgs = (/  3.70D16,     &
                        9.46D18,     &
                        1.185D20,    &
                        2.00D8,      &
                        7.96D29,     &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0),   &
                        tiny(0D0) /)
        coolB =     (/  2.12,        &
                        1.0,         &
                        0.56,        &
                        3.67,        &
                        -0.65,       &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.),    &
                        tiny(0.)  /)
        ncool=5
!
!  Revised to make continuous
!
      else if (cooling_select == 'SSrr') then
        if (lroot) print*,'select_cooling: revised SS cooling fct'
        coolT_cgs = (/  10.0,                 &
                        141.0,                &
                        313.0,                &
                        6102.0,               &
                        1E5,                  &
                        4E7,                  &
                        1E17,                 &
                        tiny(0E0),            &
                        tiny(0E0),            &
                        tiny(0E0),            &
                        tiny(0.0)  /)
        coolH_cgs = (/  3.703109927416290D16, &
                        9.455658188464892D18, &
                        1.185035244783337D20, &
                        1.9994576479D8,       &
                        7.96D29,              &
                        1.440602814622207D21, &
                        tiny(0D0),            &
                        tiny(0D0),            &
                        tiny(0D0),            &
                        tiny(0D0),            &
                        tiny(0D0)  /)
        coolB =     (/  2.12,                 &
                        1.0,                  &
                        0.56,                 &
                        3.67,                 &
                        -0.65,                &
                        0.5,                  &
                        tiny(0.),             &
                        tiny(0.),             &
                        tiny(0.),             &
                        tiny(0.),             &
                        tiny(0.)   /)
        ncool=6
!
!  26-Jan-10/fred
!  Combines Sanchez-Salcedo (2002) with Slyz et al (2005) above 1E5K
!  as Gressel simulation (2008) with constants revised for continuity
!
      else if (cooling_select == 'WSW') then
        if (lroot) print*,'select_cooling: WSW cooling fct'
        coolT_cgs = (/  90.0,                 &
                        141.0,                &
                        313.0,                &
                        6102.0,               &
                        1E5,                  &
                        2.88E5,               &
                        4.73E5,               &
                        2.11E6,               &
                        3.98E6,               &
                        2.0E7,                &
                        1.0E17      /)
        coolH_cgs = (/  3.703109927416290D16, &
                        9.455658188464892D18, &
                        1.185035244783337D20, &
                        1.102120336D10,       &
                        1.236602671D27,       &
                        2.390722374D42,       &
                        4.003272698D26,       &
                        1.527286104D44,       &
                        1.608087849D22,       &
                        9.228575532D20,       &
                        tiny(0D0) /)
        coolB =     (/  2.12,                 &
                        1.0,                  &
                        0.56,                 &
                        3.21,                 &
                       -0.20,                 &
                       -3.0,                  &
                       -0.22,                 &
                       -3.00,                 &
                        0.33,                 &
                        0.50,                 &
                        tiny(0.) /)
        ncool=10
!
!  As above but with cooling truncated below 2e4 K for Cioffi comparison
!
      else if (cooling_select == 'WSWr') then
        if (lroot) print*,'select_cooling: WSWr cooling fct'
        coolT_cgs = (/  2E2,                  &
                        1.2E4,                &
                        1E5,                  &
                        2.88E5,               &
                        4.73E5,               &
                        2.11E6,               &
                        3.98E6,               &
                        2.0E7,                &
                        1E17,                 &
                        2E17,                 &
                        2E17      /)
        coolH_cgs = (/  7.762822459260039D-19,&
                        1.102120336D10,       &
                        1.236602671D27,       &
                        2.390722374D42,       &
                        4.003272698D26,       &
                        1.527286104D44,       &
                        1.608087849D22,       &
                        9.228575532D20,       &
                        tiny(0D0),            &
                        tiny(0D0),            &
                        tiny(0D0) /)
        coolB =     (/  10.0,                 &
                        3.21,                 &
                       -0.20,                 &
                       -3.0,                  &
                       -0.22,                 &
                       -3.00,                 &
                        0.33,                 &
                        0.5,                  &
                        tiny(0.),             &
                        tiny(0.),             &
                        tiny(0.) /)
        ncool=10
      else if (cooling_select == 'off') then
        if (lroot) print*,'select_cooling: no cooling applied'
        coolT_cgs=tiny(0.0)
        coolH_cgs=tiny(0.0)
        coolB=tiny(0.)
      endif
!
!  BEGIN TEMPORARY
      if (any(coolH_cgs(1:ncool) == 0) .or. any(coolT_cgs(1:ncool+1) == 0)) &
        call fatal_error('select_cooling', &
                         'Calculating lncoolH and lncoolT: One of the cooling coefficient is zero')
!  END TEMPORARY
      lncoolH(1:ncool) = real(log(coolH_cgs(1:ncool)) - log(unit_Lambda) &
                              + log(unit_temperature**coolB(1:ncool)) &
                              - lnmu2 &
                              + log(coolingfunction_scalefactor))
      lncoolT(1:ncool+1) = real(log(coolT_cgs(1:ncool+1) / unit_temperature))
!
    endsubroutine select_cooling
!*****************************************************************************
    subroutine input_persist_interstellar_id(id,done)
!
!  Read in the parameters of the next SNI
!
!  13-Dec-2011/Bourdin.KIS: reworked
!  14-jul-2015/fred: removed obsolete Remnant persistant variable from current
!  read and added new cluster variables. All now consistent with any io
!
      use IO, only: read_persist, lun_input
!
      integer, intent(in) :: id
      logical, intent(inout) :: done
!
      integer :: i
      real :: dummy
      logical :: ldummy
!
!      if (lcollective_IO) call fatal_error ('input_persist_interstellar_id', &
!          "The interstellar persistent variables can't be read collectively!")
!
      select case (id)
        ! for backwards-compatibility (deprecated):
        case (id_record_ISM_T_NEXT_OLD)
          read (lun_input) t_next_SNI, t_next_SNII
          done = .false.
        case (id_record_ISM_POS_NEXT_OLD)
          read (lun_input) x_cluster, y_cluster
          done = .false.
        case (id_record_ISM_TOGGLE_OLD)
          read (lun_input) lSNI, lSNII
          done = .false.
        case (id_record_ISM_BOLD_MASS)
          if (read_persist ('ISM_BOLD_MASS', boldmass)) return
          done = .false.
        case (id_record_ISM_SNRS_OLD)
          ! Forget any existing SNRs.
          SNRs%indx%state = SNstate_invalid
          read (lun_input) nSNR
          do i = 1, nSNR
            read (lun_input) SNRs(i)
            SNR_index(i) = i
          enddo
          done = .false.
        ! currently active tags:
        case (id_record_ISM_T_NEXT_SNI)
          if (l_persist_overwrite_tSNI) then
            done=read_persist ('ISM_T_NEXT_SNI', dummy)
          else
            call warning('input_persist_interstellar','t_next_SNI from run.in ' // &
                         'overwritten. Set l_persist_overwrite_tSNI=T to update')
            done=read_persist ('ISM_T_NEXT_SNI', t_next_SNI)
          endif
        case (id_record_ISM_T_NEXT_SNII)
          if (l_persist_overwrite_tSNII) then
            done=read_persist ('ISM_T_NEXT_SNII', dummy)
          else
            call warning('input_persist_interstellar','t_next_SNII from run.in ' // &
                         'overwritten. Set l_persist_overwrite_tSNII=T to update')
            done=read_persist ('ISM_T_NEXT_SNII', t_next_SNII)
          endif
        case (id_record_ISM_X_CLUSTER)
          if (l_persist_overwrite_xcluster) then
            done=read_persist ('ISM_X_CLUSTER', dummy)
          else
            call warning('input_persist_interstellar','x_cluster from run.in ' // &
                         'overwritten. Set l_persist_overwrite_xcluster=T to update')
            done=read_persist ('ISM_X_CLUSTER', x_cluster)
            if (lroot.and.lOB_cluster.and..not.done) print *, 'input_persist_interstellar: x_cluster = ', x_cluster
          endif
        case (id_record_ISM_Y_CLUSTER)
          if (l_persist_overwrite_ycluster) then
            done=read_persist ('ISM_Y_CLUSTER', dummy)
          else
            call warning('input_persist_interstellar','y_cluster from run.in ' // &
                         'overwritten. Set l_persist_overwrite_ycluster=T to update')
            done=read_persist ('ISM_Y_CLUSTER', y_cluster)
            if (lroot.and.lOB_cluster.and..not.done) print *, 'input_persist_interstellar: y_cluster = ', y_cluster
          endif
        case (id_record_ISM_Z_CLUSTER)
          if (l_persist_overwrite_zcluster) then
             done=read_persist ('ISM_Z_CLUSTER', dummy)
          else
            call warning('input_persist_interstellar','z_cluster from run.in ' // &
                         'overwritten. Set l_persist_overwrite_zcluster=T to update')
            done=read_persist ('ISM_Z_CLUSTER', z_cluster)
            if (lroot.and.lOB_cluster.and..not.done) print *, 'input_persist_interstellar: z_cluster = ', z_cluster
          endif
        case (id_record_ISM_T_CLUSTER)
          if (l_persist_overwrite_tcluster) then
            done=read_persist ('ISM_T_CLUSTER', dummy)
          else
            call warning('input_persist_interstellar','t_cluster from run.in ' // &
                         'overwritten. Set l_persist_overwrite_tcluster=T to update')
            done=read_persist ('ISM_T_CLUSTER', t_cluster)
            if (lroot.and.lOB_cluster.and..not.done) print *, 'input_persist_interstellar: t_cluster = ', t_cluster
          endif
        case (id_record_ISM_TOGGLE_SNI)
          if (l_persist_overwrite_lSNI) then
            done=read_persist ('ISM_TOGGLE_SNI', ldummy)
          else
            call warning('input_persist_interstellar','lSNI from run.in ' // &
                         'overwritten. Set l_persist_overwrite_lSNI=T to update')
            done=read_persist ('ISM_TOGGLE_SNI', lSNI)
            if (lroot.and..not.done) print *, 'input_persist_interstellar: lSNI  = ',  lSNI, ' t_next_SNI  = ', t_next_SNI
          endif
        case (id_record_ISM_TOGGLE_SNII)
          if (l_persist_overwrite_lSNII) then
            done=read_persist ('ISM_TOGGLE_SNII', ldummy)
          else
            call warning('input_persist_interstellar','lSNII from run.in ' // &
                         'overwritten. Set l_persist_overwrite_lSNII=T to update')
            done=read_persist ('ISM_TOGGLE_SNII', lSNII)
            if (lroot.and..not.done) print *, 'input_persist_interstellar: lSNII = ', lSNII, &
                                              ' t_next_SNII = ', t_next_SNII
          endif
      endselect
!
    endsubroutine input_persist_interstellar_id
!*****************************************************************************
    subroutine input_persist_interstellar()
!
!  Read in the stored time of the next SNI.
!
!  12-Oct-2019/PABourdin: coded
!
      use IO, only: read_persist
!
      logical :: error
!
      if (.not. l_persist_overwrite_lSNI) then
        call warning('input_persist_interstellar','lSNI from run.in overwritten. ' // &
                     'Set l_persist_overwrite_lSNI=T to update', 0)
        error = read_persist ('ISM_TOGGLE_SNI', lSNI)
        if (lroot .and. .not. error) print *, 'input_persist_interstellar: lSNI = ', lSNI
      endif
!
      if (.not. l_persist_overwrite_lSNII) then
        call warning('input_persist_interstellar','lSNII from run.in overwritten. ' // &
                     'Set l_persist_overwrite_lSNII=T to update', 0)
        error = read_persist ('ISM_TOGGLE_SNII', lSNII)
        if (lroot .and. .not. error) print *, 'input_persist_interstellar: lSNII = ', lSNII
      endif
!
      if (.not. l_persist_overwrite_tSNI) then
        call warning('input_persist_interstellar','t_next_SNI from run.in overwritten. ' // &
                     'Set l_persist_overwrite_tSNI=T to update', 0)
        error = read_persist ('ISM_T_NEXT_SNI', t_next_SNI)
        if (lroot .and. .not. error) print*, 'input_persist_interstellar: t_next_SNI = ', t_next_SNI
      endif
!
      if (.not. l_persist_overwrite_tSNII) then
        call warning('input_persist_interstellar','t_next_SNII from run.in overwritten. ' // &
                     'Set l_persist_overwrite_tSNII=T to update', 0)
        error = read_persist ('ISM_T_NEXT_SNII', t_next_SNII)
        if (lroot .and. .not. error) print *, 'input_persist_interstellar: t_next_SNII = ', t_next_SNII
      endif
!
      if (.not. l_persist_overwrite_xcluster) then
        call warning('input_persist_interstellar','x_cluster from run.in overwritten. ' // &
                     'Set l_persist_overwrite_xcluster=T to update', 0)
        error = read_persist ('ISM_X_CLUSTER', x_cluster)
        if (lOB_cluster .and. lroot .and. .not. error) print*,'input_persist_interstellar: x_cluster = ',x_cluster
      endif
!
      if (.not. l_persist_overwrite_ycluster) then
        call warning('input_persist_interstellar','y_cluster from run.in overwritten. ' // &
                     'Set l_persist_overwrite_ycluster=T to update', 0)
        error = read_persist ('ISM_Y_CLUSTER', y_cluster)
        if (lOB_cluster .and. lroot .and. .not. error) print*,'input_persist_interstellar: y_cluster = ',y_cluster
      endif
!
      if (.not. l_persist_overwrite_zcluster) then
        call warning('input_persist_interstellar','z_cluster from run.in overwritten. ' // &
                     'Set l_persist_overwrite_zcluster=T to update', 0)
        error = read_persist ('ISM_Z_CLUSTER', z_cluster)
        if (lOB_cluster .and. lroot .and. .not. error) print*,'input_persist_interstellar: z_cluster = ',z_cluster
      endif
!
      if (.not. l_persist_overwrite_tcluster) then
        call warning('input_persist_interstellar','t_cluster from run.in overwritten. '// &
                     'Set l_persist_overwrite_tcluster=T to update', 0)
        error = read_persist ('ISM_T_CLUSTER', t_cluster)
        if (lOB_cluster .and. lroot .and. .not. error) print*,'input_persist_interstellar: t_cluster = ',t_cluster
      endif
!
    endsubroutine input_persist_interstellar
!*****************************************************************************
    logical function output_persistent_interstellar()
!
!  Writes out the time of the next SNI
!
!  13-Dec-2011/Bourdin.KIS: reworked
!  14-jul-2015/fred: removed obsolete Remnant persistant variable from current
!  write and added new cluster variables. All now consistent with any io
!
      use IO, only: write_persist
!
!      if (lcollective_IO) call fatal_error ('output_persistent_interstellar', &
!          "The interstellar persistent variables can't be written collectively!")
!
      output_persistent_interstellar = .true.
!
      if (write_persist ('ISM_T_NEXT_SNI', id_record_ISM_T_NEXT_SNI, t_next_SNI)) return
      if (write_persist ('ISM_T_NEXT_SNII', id_record_ISM_T_NEXT_SNII, t_next_SNII)) return
      if (write_persist ('ISM_X_CLUSTER', id_record_ISM_X_CLUSTER, x_cluster)) return
      if (write_persist ('ISM_Y_CLUSTER', id_record_ISM_Y_CLUSTER, y_cluster)) return
      if (write_persist ('ISM_Z_CLUSTER', id_record_ISM_Z_CLUSTER, z_cluster)) return
      if (write_persist ('ISM_T_CLUSTER', id_record_ISM_T_CLUSTER, t_cluster)) return
      if (write_persist ('ISM_TOGGLE_SNI', id_record_ISM_TOGGLE_SNI, lSNI)) return
      if (write_persist ('ISM_TOGGLE_SNII', id_record_ISM_TOGGLE_SNII, lSNII)) return
!
      output_persistent_interstellar = .false.
!
    endfunction output_persistent_interstellar
!*****************************************************************************
    subroutine rprint_interstellar(lreset,lwrite)
!
!  Reads and registers print parameters relevant to interstellar
!
!  01-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics, only: parse_name
!
      integer :: iname,inamez
      logical :: lreset
      logical, intent(in), optional :: lwrite
!
!  Reset everything in case of reset
!  (This needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_taucmin=0
        idiag_Hmax_ism=0
        idiag_Lamm=0
        idiag_nrhom=0
        idiag_rhoLm=0
        idiag_Gamm=0
        idiag_rhoHCmz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'taucmin',idiag_taucmin)
        call parse_name(iname,cname(iname),cform(iname),'Hmax_ism',idiag_Hmax_ism)
        call parse_name(iname,cname(iname),cform(iname),'Lamm',idiag_Lamm)
        call parse_name(iname,cname(iname),cform(iname),'nrhom',idiag_nrhom)
        call parse_name(iname,cname(iname),cform(iname),'rhoLm',idiag_rhoLm)
        call parse_name(iname,cname(iname),cform(iname),'Gamm',idiag_Gamm)
      enddo
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhoHCmz',idiag_rhoHCmz)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='ism_cool'.or.cnamev=='ism_netheat') cformv='DEFINED'
      endif
!
    endsubroutine rprint_interstellar
!*****************************************************************************
    subroutine get_slices_interstellar(f,slices)
!
!  Write slices for animation of Interstellar variables.
!
!  26-jul-06/tony: coded
!
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Shock profile
!
        case ('ism_cool'); call assign_slices_scal(slices,f,icooling)
        case ('ism_netheat'); call assign_slices_scal(slices,f,inetheat)
!
      endselect
!
    endsubroutine get_slices_interstellar
!***********************************************************************
    subroutine read_interstellar_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=interstellar_init_pars, IOSTAT=iostat)
!
    endsubroutine read_interstellar_init_pars
!***********************************************************************
    subroutine write_interstellar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=interstellar_init_pars)
!
    endsubroutine write_interstellar_init_pars
!***********************************************************************
    subroutine read_interstellar_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=interstellar_run_pars, IOSTAT=iostat)
!
    endsubroutine read_interstellar_run_pars
!***********************************************************************
    subroutine write_interstellar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=interstellar_run_pars)
!
    endsubroutine write_interstellar_run_pars
!***********************************************************************
    subroutine init_interstellar(f)
!
!  Initialise some explosions etc.
!  24-nov-2002/tony: coded
!
      use General, only: itoa, random_seed_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      logical :: lnothing=.true.
      character (len=intlen) :: iinit_str
      integer :: i,j,iSNR
!
      intent(inout) :: f
!
!  Fred: 07Oct17 The random numbers must be synchronized on all processors or
!        else the treatment of SN explosions shall diverge and the MPI will
!        break or worse hang
!
      seed(1)=seed_reset
      call random_seed_wrapper(PUT=seed)
!
      do j=1,ninit
!
      if (initinterstellar(j)/='nothing') then
!
      lnothing=.false.
      iinit_str=itoa(j)
!
!  Select different initial conditions
!
      select case (initinterstellar(j))
!
        case ('single')
          iSNR=get_free_SNR()
          SNRs(iSNR)%site%TT=1E20
          SNRs(iSNR)%site%rho=0.
          SNRs(iSNR)%feat%t=t
          SNRs(iSNR)%indx%SN_type=1
          SNRs(iSNR)%feat%radius=width_SN
          call position_SN_testposition(f,SNRs(iSNR))
          call explode_SN(f,SNRs(iSNR))
          lSNI=.false.
          lSNII=.false.
        case ('sedov')
          iSNR=get_free_SNR()
          SNRs(iSNR)%site%TT=1E20
          SNRs(iSNR)%site%rho=0.
          SNRs(iSNR)%feat%t=t
          SNRs(iSNR)%indx%SN_type=1
          SNRs(iSNR)%feat%radius=width_SN
          center_SN_x=0.
          center_SN_y=0.
          center_SN_z=0.
          call position_SN_testposition(f,SNRs(iSNR))
          call explode_SN(f,SNRs(iSNR))
          lSNI=.false.
          lSNII=.false.
        case ('courant-friedricks')
          iSNR=get_free_SNR()
          SNRs(iSNR)%site%TT=1E20
          SNRs(iSNR)%site%rho=0.
          SNRs(iSNR)%feat%t=t
          SNRs(iSNR)%indx%SN_type=1
          SNRs(iSNR)%feat%radius=width_SN
          center_SN_x=0.
          center_SN_y=0.
          center_SN_z=-0.015
          call position_SN_testposition(f,SNRs(iSNR))
          call explode_SN(f,SNRs(iSNR))
          iSNR=get_free_SNR()
          SNRs(iSNR)%site%TT=1E20
          SNRs(iSNR)%site%rho=0.
          SNRs(iSNR)%feat%t=t
          SNRs(iSNR)%feat%radius=width_SN
          center_SN_x=0.
          center_SN_y=0.
          center_SN_z=0.015
          call position_SN_testposition(f,SNRs(iSNR))
          call explode_SN(f,SNRs(iSNR))
          lSNI=.false.
          lSNII=.false.
        case ('kompaneets')
          iSNR=get_free_SNR()
          SNRs(iSNR)%site%TT=1E20
          SNRs(iSNR)%site%rho=0.
          SNRs(iSNR)%feat%t=0
          SNRs(iSNR)%indx%SN_type=1
          SNRs(iSNR)%feat%radius=width_SN
          center_SN_x=0.
          center_SN_y=0.
          center_SN_z=0.
          call position_SN_testposition(f,SNRs(iSNR))
          call explode_SN(f,SNRs(iSNR))
          lSNI=.false.
          lSNII=.false.
        case ('multiple')
          do i = 1,initial_SNI
            iSNR=get_free_SNR()
            SNRs(iSNR)%site%TT=1E20
            SNRs(iSNR)%site%rho=0.
            SNRs(iSNR)%feat%t=0.
            SNRs(iSNR)%indx%SN_type=1
            SNRs(iSNR)%feat%radius=width_SN
            if (luniform_zdist_SNI) then
              call position_SN_uniformz(f,SNRs(i))
            else
              lfirst_zdisk=.false.
              zdisk=0.
              call position_SN_gaussianz(f,h_SNII,SNRs(i))
            endif
            call explode_SN(f,SNRs(i))
          enddo
!
        case default
!
!  Catch unknown values
!
        call fatal_error('init_interstellar','no such initinterstellar(' &
                         //trim(iinit_str)//'): '//trim(initinterstellar(j)))
      endselect
!
      if (lroot) print*,'init_interstellar: initinterstellar(' &
                        //trim(iinit_str)//') = ',trim(initinterstellar(j))
      endif
!
      enddo
!
      if (lnothing.and.lroot) print*,'init_interstellar: nothing'
!
      call tidy_SNRs
!
    endsubroutine init_interstellar
!*****************************************************************************
    subroutine pencil_criteria_interstellar()
!
!  All pencils that the Interstellar module depends on are specified here.
!
!  26-mar-05/tony: coded
!  11-mar-06/axel: added idiag_nrhom
!
      lpenc_requested(i_heat)=.true.
      lpenc_requested(i_cool)=.true.
      lpenc_requested(i_lnrho)=.true.
      lpenc_requested(i_lnTT)=.true.
      if (ltemperature.and..not.ltemperature_nolog) lpenc_requested(i_ee)=.true.
      if (.not.(ltemperature_nolog.and.ltemperature).or.pretend_lnTT) then
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_TT1)=.true.
      endif
      if (ltemperature_nolog) lpenc_requested(i_cv1)=.true.
      if (lheatcool_shock_cutoff) lpenc_requested(i_gshock)=.true.
!
      if (idiag_taucmin/=0.or.idiag_Hmax_ism/=0) then
        if (ltemperature_nolog.and.ltemperature) then
          lpenc_diagnos(i_TT1)=.true.
        else
          lpenc_diagnos(i_ee)=.true.
          lpenc_diagnos(i_TT)=.true.
        endif
        lpenc_diagnos(i_heatcool)=.true.
      endif
      if (idiag_rhoHCmz/=0) lpenc_diagnos(i_heatcool)=.true.
      if (idiag_Lamm/=0) lpenc_diagnos(i_rho1)=.true.
      if (idiag_nrhom/=0) lpenc_diagnos(i_ee)=.true.
!
    endsubroutine pencil_criteria_interstellar
!***********************************************************************
    subroutine pencil_interdep_interstellar(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_cool)) then
        lpencil_in(i_lnTT)=.true.
        lpencil_in(i_lnrho)=.true.
      endif
      if (lpencil_in(i_heat)) then
        lpencil_in(i_lnTT)=.true.
      endif
      if (lpencil_in(i_heatcool)) then
        lpencil_in(i_cool)=.true.
        lpencil_in(i_heat)=.true.
      endif
!
    endsubroutine pencil_interdep_interstellar
!***********************************************************************
    subroutine calc_pencils_interstellar(f,p)
!
!  Calculate Interstellar pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      real, dimension(mx,my,mz,mfarray), intent(IN)   :: f
      type(pencil_case),                 intent(INOUT):: p
!
      if (lpencil(i_cool)) call calc_cool_func(p%cool,p%lnTT,p%lnrho)
!
      if (lpencil(i_heat)) call calc_heat(p%heat,p%TT)
!
!  For clarity with lentropy we have constructed the rhs in erg/s/g [=T*Ds/Dt]
!  so therefore we now need to multiply by TT1, or otherwise.
!
      if (ltemperature) then
        if (ltemperature_nolog) then
          p%heat=p%heat*p%cv1
          p%cool=p%cool*p%cv1
        else
          p%heat=p%heat/p%ee
          p%cool=p%cool/p%ee
        endif
      elseif (pretend_lnTT) then
        p%heat=p%heat*gamma
        p%cool=p%cool*gamma
      else   !? if (lentropy) then
        p%heat=p%heat*p%TT1
        p%cool=p%cool*p%TT1
        !heatcool=p%TT1*(heat-cool)
      endif
      if (lpencil(i_heatcool)) p%heatcool=p%heat-p%cool
!
    endsubroutine calc_pencils_interstellar
!***********************************************************************
    subroutine interstellar_after_boundary(f)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  01-aug-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine interstellar_after_boundary
!***********************************************************************
    subroutine interstellar_before_boundary(f)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  01-aug-06/tony: coded
!
!$    use OMP_lib
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!$    call OMP_set_num_threads(num_threads)
      if (lfirst.and..not.lpencil_check_at_work) call check_SN(f)
!
      if (lSNI.or.lSNII) then
        if (laverage_SNI_heating) &
          heatingfunction_scale(SNI) = t_interval(SNI) /(t_interval(SNI) +t*heatingfunction_fadefactor) &
                                      *heatingfunction_scalefactor
        if (laverage_SNII_heating) &
          heatingfunction_scale(SNII)= t_interval(SNII)/(t_interval(SNII)+t*heatingfunction_fadefactor) &
                                      *heatingfunction_scalefactor
      endif

    endsubroutine interstellar_before_boundary
!***********************************************************************
    subroutine calc_0d_diag_interstellar(p)
!
!  Calculate entropy related diagnostics.
!
      use Diagnostics
!
      type(pencil_case) :: p
!
      if (ldiagnos) then
        if (idiag_Hmax_ism/=0) then
          if (ltemperature.and.ltemperature_nolog) then
            call max_mn_name(p%heatcool*p%TT1,idiag_Hmax_ism)
          elseif (pretend_lnTT) then
            call max_mn_name(p%heatcool,idiag_Hmax_ism)
          else
            call max_mn_name(p%heatcool*p%TT/p%ee,idiag_Hmax_ism)
          endif
        endif
        if (idiag_taucmin/=0) then
          if (ltemperature.and.ltemperature_nolog) then
            call max_mn_name(-p%heatcool*p%TT1,idiag_taucmin,lreciprocal=.true.)
          elseif (pretend_lnTT) then
            call max_mn_name(-p%heatcool,idiag_taucmin,lreciprocal=.true.)
          else
            call max_mn_name(-p%heatcool*p%TT/p%ee,idiag_taucmin,lreciprocal=.true.)
          endif
        endif
        if (idiag_Lamm/=0) call sum_mn_name(p%rho1*p%cool,idiag_Lamm)
        if (idiag_nrhom/=0) call sum_mn_name(p%cool/p%ee,idiag_nrhom)
        call sum_mn_name(p%cool,idiag_rhoLm)
        call sum_mn_name(p%heat,idiag_Gamm)
      endif
!
    endsubroutine calc_0d_diag_interstellar
!***********************************************************************
    subroutine calc_1d_diag_interstellar(p)
!
      use Diagnostics, only: max_mn_name, sum_mn_name, xysum_mn_name_z
!
      type(pencil_case) :: p
!
!  1-D averages.
!
      if (l1davgfirst) then
        call xysum_mn_name_z(p%heatcool,idiag_rhoHCmz)
      endif
!
    endsubroutine calc_1d_diag_interstellar
!***********************************************************************
    subroutine calc_diagnostics_interstellar(p)

      type(pencil_case) :: p

      call calc_1d_diag_interstellar(p)
      call calc_0d_diag_interstellar(p)

    endsubroutine calc_diagnostics_interstellar
!*****************************************************************************
    subroutine heat_interstellar(f,zheat)
!
!  This routine calculates a vertical profile for uv-heating designed to
!  satisfy an initial condition with heating and cooling approximately balanced
!  for an isothermal hydrostatic equilibrium.
!  Requires: gravz_profile='Ferriere' in gravity_simple.f90
!            heating_select='thermal-hs' in interstellar.f90
!  Using here a similar method to O. Gressel 2008 (PhD) lthermal_hse=T
!  or similar to Joung & Mac Low Apj 653 Dec 2006 without hse
!
!  22-mar-10/fred:
!  adapted from galactic-hs,ferriere-hs
!  13-jul-15/fred: requires initial_condition/hs_equilibrium_ism.f90
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mz), intent(out) :: zheat
!
      real, dimension(mz) :: lambda=0.0, lnTT, zrho
      real :: logrho, TT
      integer :: j
!
!  Identifier
!
      if (lroot.and.ip==1963) print*, &
         'heat_interstellar: calculating z-dependent uv-heating'// &
         'function for initial hydrostatic and thermal equilibrium'
!
!  Set up physical units.
!
      if (unit_system=='cgs') then
        a_S = a_S_cgs/unit_velocity*unit_time
        a_D = a_D_cgs/unit_velocity*unit_time
        z_D = z_D_cgs/unit_length
        z_S = z_S_cgs/unit_length
        H_z = H_z_cgs/unit_length
      else if (unit_system=='SI') then
        call not_implemented('heat_interstellar', 'SI unit conversions')
      endif
!
      do n=1,mz
        if (lthermal_hse) then
          logrho = log(rho0ts)+(a_S*z_S*m_u*mu/k_B/T_init)*(log(T_init)- &
              log(T_init/(a_S*z_S)* &
              (a_S*sqrt(z_S**2+(z(n))**2)+0.5*a_D*(z(n))**2/z_D)))
        else
          logrho = log(rho0ts)-0.015*(- &
              a_S*z_S+ &
              a_S*sqrt(z_S**2+(z(n))**2)+0.5*a_D*(z(n))**2/z_D)
        endif
        if (logrho < -40.0) logrho=-40.0
        zrho(n)=exp(logrho)
        TT=T_init/(a_S*z_S)*(a_S*sqrt(z_S**2+(z(n))**2)+0.5*a_D*(z(n))**2/z_D)
        lnTT(n)=log(TT)
        zheat(n)=GammaUV*exp(-abs(z(n))/H_z)
      enddo
      do j=1,ncool
        where (lncoolT(j)<=lnTT.and.lnTT<lncoolT(j+1)) lambda=lambda+exp(lncoolH(j)+lnTT*coolB(j))
      enddo
      if (lthermal_hse) zheat=lambda*zrho

      if (lheatz_min) then
        where (zheat<1E-5*GammaUV) zheat=1E-5*GammaUV
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine heat_interstellar
!*****************************************************************************
    subroutine calc_heat_cool_interstellar(f,df,p,Hmax)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  This public subroutine is called by entropy.f90. If other energy/temp/entopy
!  modules are selected an equivalent call should be included.
!  We may want to move it to the entropy module for good, because its use
!  is not restricted to interstellar runs (could be used for solar corona).
!  Also, it doesn't pose an extra load on memory usage or compile time.
!  (We should allow that UV heating can be turned off; so rhoUV should
!  be made an input parameter.)
!
!  19-nov-02/graeme: adapted from calc_heat_cool
!  10-aug-03/axel: TT is used as input
!   3-apr-06/axel: add ltemperature switch
!
      use Sub, only: dot2
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx), intent(inout) :: Hmax

      real, dimension (nx) :: heat,cool,heatcool
      real, dimension (nx) :: damp_profile, gsh2
!
!  Identifier
!
      if (headtt) print*,'calc_heat_cool_interstellar: ENTER'
!
!  Not all eos define cv1, and this routine is called after
!  calc_all_pencils, so p%cv1 cannot be modified here
!
      if (ltemperature_nolog.and.any(p%cv1 == impossible)) &
        call fatal_error("calc_heat_cool_interstellar","p%cv1 not set by eos")
!
      heat=p%heat
      cool=p%cool
!
!  Average SN heating (due to SNI and SNII)
!  The amplitudes of both types is assumed the same (=ampl_SN)
!  Added option to gradually remove average heating as SN are introduced when
!  initial condition is in equilibrium prepared in 1D
!  Division by density to balance LHS of entropy equation
!
      if (laverage_SNI_heating) &
          heat=heat+average_SNI_heating *exp(-(2.0*z(n)/h_SNI )**2)*heatingfunction_scale(SNI)
      if (laverage_SNII_heating) &
          heat=heat+average_SNII_heating*exp(-(2.0*z(n)/h_SNII)**2)*heatingfunction_scale(SNII)
!
!  Prevent unresolved heating/cooling in shocks. This is deprecated as
!  use of RKF timestep control or the mor coarse RHS timestep control
!  appended to the Courant timestep control is effective without
!  unphysical artefacts. Fred
!
      if (lheatcool_shock_cutoff) then
        call dot2(p%gshock,gsh2)
!
        damp_profile=exp(-(gsh2*heatcool_shock_cutoff_rate1))
!
        cool=cool*damp_profile
        heat=heat*damp_profile
      endif
!
      heatcool=heat-cool
!
!  Save result in aux variables
!  cool=rho*Lambda/TT, heatcool=(Gamma-rho*Lambda)/TT
!
      f(l1:l2,m,n,icooling) = cool
      f(l1:l2,m,n,inetheat) = heatcool
!
      call calc_diagnostics_interstellar(p)
!
!  Limit timestep by the cooling time (having subtracted any heating)
!
!  FG: todo this applies a primarily negative heatcool to positive Hmax
!      which controls the ldt time step in energy RHS. Require to sum
!      Hmax with abs(heatcool) or to compute cooling timestep separately
!      Hmax may need to be max(abs(Hmax)) for ldt??
!
      if (ldt.and.lfirst .or. ldiagnos) then
        if (ltemperature.or.pretend_lnTT) then
          Hmax=Hmax+heatcool
        else
          Hmax=Hmax+heatcool*p%TT
        endif
      endif
!
!  Apply heating/cooling to temperature/entropy variable
!
      if (ltemperature) then
        df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+heatcool
      else
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+heatcool
      endif
!
    endsubroutine calc_heat_cool_interstellar
!*****************************************************************************
    subroutine calc_cool_func(cool,lnTT,lnrho)
!
!  This routine calculates the temperature dependent radiative cooling.
!  Applies Rosen et al., ApJ, 413, 137, 1993 ('RBN') OR
!  Sanchez-Salcedo et al. ApJ, 577, 768, 2002 ('SS') OR
!  Slyz et al MNRAS, 356 2005 ('WSW') fit of Wolfire with Sarazin & White
!  ApJ, 443:152-168, 1985 and ApJ, 320:32-48, 1987
!
!
!  Cooling is Lambda*rho^2, with (eq 7) & Lambda=coolH(i)*TT**coolB(i),
!  for coolT(i) <= TT < coolT(i+1). Nb: our coefficients coolH(i) differ from
!  those in Rosen et al. by a factor (mu mp)^2, with mu=1.2, since Rosen works
!  in number density, n. (their cooling = Lambda*n^2,  rho=mu mp n)
!  The factor Lambdaunits converts from cgs units to code units.
!
!  [Currently, coolT(1) is not modified, but this may be necessary
!  to avoid creating gas too cold to resolve.]
!
      real, dimension (nx), intent(out) :: cool
      real, dimension (nx), intent(in) :: lnTT, lnrho
      integer :: i
!
      cool=0.0
      do i=1,ncool
        where (lncoolT(i) <= lnTT .and. lnTT < lncoolT(i+1))
          cool=cool+exp(lncoolH(i)+lnrho+lnTT*coolB(i))
        endwhere
      enddo

    endsubroutine calc_cool_func
!*****************************************************************************
    subroutine calc_heat(heat,TT)
!
!  This routine adds UV heating, cf. Wolfire et al., ApJ, 443, 152, 1995
!  with the values above, this gives about 0.012 erg/g/s (T < ~1.E4 K)
!  Nb: need rho0 from density_[init/run]_pars to implement the arm/interarm
!  scaling.
!
!  Control with heating_select in interstellar_init_pars/run_pars.
!  Default heating_rate GammaUV = 0.015.
!
      real, dimension (nx), intent(out) :: heat
      real, dimension (nx), intent(in) :: TT
!
!  Constant heating with a rate heating_rate[erg/g/s].
!
      if (heating_select == 'cst') then
        heat = heating_rate_code
      else if (heating_select == 'wolfire') then
        heat = GammaUV*0.5*(1.0+tanh(cUV*(T0UV-TT)))
      else if (heating_select == 'wolfire_min') then
        heat = GammaUV*0.5*(1.0+tanh(cUV*(T0UV-TT)))
        heat = max(heat,heating_rate_code)
!
!  If using thermal-hs in initial entropy this must also be specified for
!  thermal equilibrium and applies for vertically stratified density supported
!  by vertical gravity profile 'Ferriere'.
!
      else if (heating_select == 'thermal-hs') then
        heat = heat_z(n)*0.5*(1.0+tanh(cUV*(T0UV-TT)))
      else if (heating_select == 'off') then
        heat = 0.
      endif
!
    endsubroutine calc_heat
!*****************************************************************************
    subroutine check_SN(f)
!
!  Checks for SNe, and implements appropriately:
!  relevant subroutines in entropy.f90
!
      use General, only: touch_file
!
      real, dimension(mx,my,mz,mfarray) :: f
!
!  Only allow SNII if no SNI this step (may not be worth keeping).
!  This is depricated, both types can occur in same step
!
      logical :: l_SNI=.false.
!
      intent(inout) :: f
      integer :: i
!
!  Identifier
!
      if (headtt) print*,'check_SN: ENTER'
!
!  If SN are listed in source file then obtain parameters from list
!
      if (lSN_list) then
        if (t>=t_next_SNI.or.t>=t_next_SNII) then
          call tidy_SNRs
          do i=SNfirst,nlist-1
            if (SN_list(1,i)>=t_next_SNI.or.SN_list(1,i)>=t_next_SNII) then
              center_SN_x=SN_list(2,i)
              center_SN_y=SN_list(3,i)
              center_SN_z=SN_list(4,i)
              type_list=SN_type(i)
              SN_list(1,i)=0.
              if (i==nlist-1) then
                !Exit run early if next on the list is the last
                call touch_file('ENDTIME')
                tmax=t
              endif
              t_next_SNI=SN_list(1,i+1)
              t_next_SNII=SN_list(1,i+1)
              if (lroot) print"(1x,'check_SN: t_next_SNI on list =',e16.8)",t_next_SNI
              SNfirst=i+1
              if (lroot) print"(1x,'check_SN: SNfirst at t =',i7,e16.8)",SNfirst,t
              if (i==nlist-2) then
                !Check whether the list has been extended since last reading
                call touch_file('RELOAD')
              endif
              exit
            endif
          enddo
          call check_SNI(f,l_SNI)
          if (t>=tmax) then
            if (lroot) print*, 'check_SN: sn_series.in list needs', &
                               ' extending or set lSN_list=F to continue'
          endif
        endif
      else
        if (t < t_settle) return
        call tidy_SNRs
        if (lSNI) call check_SNI(f,l_SNI)
!
!  Do separately for SNI (simple scheme) and SNII (Boris' scheme).
!
        if (lSNII) then
          if (lSNII_gaussian) then
            call check_SNIIb(f,l_SNI)
          else
            call check_SNII(f,l_SNI)
          endif
        endif
      endif
!
    endsubroutine check_SN
!*****************************************************************************
    subroutine check_SNI(f,l_SNI)
!
!  If time for next SNI, then implement, and calculate time of subsequent SNI.
!
      real, dimension(mx,my,mz,mfarray) :: f
      logical :: l_SNI
      integer :: try_count, iSNR, ierr
!
      intent(inout) :: f,l_SNI
!
!  Identifier
!
      if (headtt.and.ip==1963) print*,'check_SNI: ENTER'
!
      if (lSN_list) then
        iSNR=get_free_SNR()
        SNRs(iSNR)%site%TT=1E20
        SNRs(iSNR)%site%rho=0.0
        SNRs(iSNR)%feat%t=t
        SNRs(iSNR)%indx%SN_type=type_list
        SNRs(iSNR)%feat%radius=width_SN
        call position_SN_testposition(f,SNRs(iSNR))
        ierr=iEXPLOSION_OK
        call explode_SN(f,SNRs(iSNR),ierr,preSN)
!
!  Free up slots in case loop fails repeatedly over many time steps.
!
        call free_SNR(iSNR)
        ierr=iEXPLOSION_OK
        return
      endif
      if (t >= t_next_SNI) then
        iSNR=get_free_SNR()
        SNRs(iSNR)%site%TT=1E20
        SNRs(iSNR)%site%rho=0.0
        SNRs(iSNR)%feat%t=t
        SNRs(iSNR)%indx%SN_type=1
        SNRs(iSNR)%feat%radius=width_SN
        try_count=10
!
        do while (try_count>0)
!
          ierr=iEXPLOSION_OK
          try_count=try_count-1
!
          if (luniform_zdist_SNI) then
            call position_SN_uniformz(f,SNRs(iSNR))
          else
            lfirst_zdisk=.false.
            zdisk=0.
            call position_SN_gaussianz(f,h_SNI,SNRs(iSNR))
          endif
!
          if (.not.lSN_scale_rad) then
            if ((SNRs(iSNR)%site%rho < rho_SN_min) .or. &
                (SNRs(iSNR)%site%TT > TT_SN_max)) cycle
          else
!
!  Avoid sites with dense mass shown to excessively cool remnants, given the
!  high thermal conductivity we are forced to adopt.
!
            if (SNRs(iSNR)%site%rho > rho_SN_max) cycle
          endif
!
          call explode_SN(f,SNRs(iSNR),ierr,preSN)
          if (ierr==iEXPLOSION_OK) then
            l_SNI=.true.
            exit
          elseif (ierr==iEXPLOSION_TOO_HOT) then
            if (lroot.and.ip==1963) print &
                "(1x,'check_SNI: TOO HOT, (x,y,z) =',3f8.4,', rho, rad =',2e11.4)", &
                SNRs(iSNR)%feat%x, SNRs(iSNR)%feat%y, SNRs(iSNR)%feat%z, &
                SNRs(iSNR)%site%rho, SNRs(iSNR)%feat%radius
          elseif (ierr==iEXPLOSION_TOO_UNEVEN) then
            if (lroot.and.ip==1963) print &
                "(1x,'check_SNI: TOO UNEVEN, (x,y,z) =',3f8.4,', rho, rad =',2e11.4)", &
                SNRs(iSNR)%feat%x, SNRs(iSNR)%feat%y, SNRs(iSNR)%feat%z,&
                SNRs(iSNR)%site%rho, SNRs(iSNR)%feat%radius
          elseif (ierr==iEXPLOSION_TOO_RARIFIED) then
            if (lroot.and.ip==1963) print &
                "(1x,'check_SNI: TOO RARIFIED, (x,y,z) =',3f8.4,', rho, rad =',2e11.4)",&
                SNRs(iSNR)%feat%x, SNRs(iSNR)%feat%y, SNRs(iSNR)%feat%z, &
                SNRs(iSNR)%site%rho, SNRs(iSNR)%feat%radius
          endif
        enddo
!
        if (try_count==0) then
          if (lroot.and.ip==1963) print*,"check_SNI: 10 RETRIES OCCURED - skipping SNI insertion"
        endif
!
!  Free up slots in case loop fails repeatedly over many time steps.
!
        call free_SNR(iSNR)
!
!  Reset ierr or else explode_SN may terminate erroneously on subsequent
!  timesteps never to be reset.
!
        ierr=iEXPLOSION_OK
      endif
!
    endsubroutine check_SNI
!*****************************************************************************
    subroutine check_SNIIb(f,l_SNI)
!
!  If time for next SNI, then implement, and calculate time of subsequent SNI.
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: try_count, iSNR, ierr
      logical :: l_SNI
!
      intent(inout) :: f, l_SNI
!
!  Identifier
!
      if (headtt.and.ip==1963) print*,'check_SNIIb: ENTER'
!
      if (t >= t_next_SNII) then
        iSNR=get_free_SNR()
        SNRs(iSNR)%site%TT=1E20
        SNRs(iSNR)%site%rho=0.0
        SNRs(iSNR)%feat%t=t
        SNRs(iSNR)%indx%SN_type=2
        SNRs(iSNR)%feat%radius=width_SN
        try_count=10
        lfirst_zdisk=.true.
        zdisk=0.
!
        do while (try_count>0)
          ierr=iEXPLOSION_OK
          try_count=try_count-1
!
          if (luniform_zdist_SNI) then
            call position_SN_uniformz(f,SNRs(iSNR))
          else
            call position_SN_gaussianz(f,h_SNII,SNRs(iSNR))
            lfirst_zdisk=.false.
          endif
!
          if (.not.lSN_scale_rad) then
            if ((SNRs(iSNR)%site%rho < rho_SN_min) .or. &
                (SNRs(iSNR)%site%TT > TT_SN_max)) then
              cycle
            endif
          else
!
!  Avoid sites with dense mass shown to excessively cool remnants, given the
!  high thermal conductivity we are forced to adopt.
!
            if ((SNRs(iSNR)%site%rho > rho_SN_max) .or. &
                (SNRs(iSNR)%site%TT < TT_SN_min)) then
              cycle
            endif
          endif
!
          call explode_SN(f,SNRs(iSNR),ierr,preSN)
          if (ierr==iEXPLOSION_OK) then
            exit
          elseif (ierr==iEXPLOSION_TOO_HOT) then
            if (lroot.and.ip==1963) print &
                "(1x,'check_SNIIb: TOO HOT, (x,y,z) =',3f8.4,', rho, rad =',2e11.4)", &
                SNRs(iSNR)%feat%x, SNRs(iSNR)%feat%y, SNRs(iSNR)%feat%z, &
                SNRs(iSNR)%site%rho, SNRs(iSNR)%feat%radius
          elseif (ierr==iEXPLOSION_TOO_UNEVEN) then
            if (lroot.and.ip==1963) print &
                "(1x,'check_SNIIb: TOO UNEVEN, (x,y,z) =',3f8.4,', rho, rad =',2e11.4)", &
                SNRs(iSNR)%feat%x, SNRs(iSNR)%feat%y, SNRs(iSNR)%feat%z,&
                SNRs(iSNR)%site%rho, SNRs(iSNR)%feat%radius
          elseif (ierr==iEXPLOSION_TOO_RARIFIED) then
            if (lroot.and.ip==1963) print &
                "(1x,'check_SNIIb: TOO RARIFIED, (x,y,z) =',3f8.4,', rho, rad =',2e11.4)", &
                SNRs(iSNR)%feat%x, SNRs(iSNR)%feat%y, SNRs(iSNR)%feat%z,&
                SNRs(iSNR)%site%rho, SNRs(iSNR)%feat%radius
          endif
        enddo
!
        if (try_count==0) then
          if (lroot.and.ip==1963) print*,"check_SNIIb: 10 RETRIES OCCURED - skipping SNII insertion"
        endif
!
!  Free up slots in case loop fails repeatedly over many time steps.
!
        call free_SNR(iSNR)
!
!  Reset ierr or else explode_SN may terminate erroneously on subsequent
!  timesteps never to be reset.
!
        ierr=iEXPLOSION_OK
      endif
      l_SNI=.true.
!
    endsubroutine check_SNIIb
!***********************************************************************
    subroutine set_next_SNI(scaled_interval)
!
      use General, only:  random_seed_wrapper, random_number_wrapper
!
      real :: franSN, scaled_interval
!
      intent(out) :: scaled_interval
!
!  Pre-determine time for next SNI.
!
      if (lroot) print"(1x,'set_next_SNI: Old t_next_SNI =',e15.8)",t_next_SNI
      if (lreset_ism_seed) then
        seed=seed_reset
        call random_seed_wrapper(PUT=seed)
        lreset_ism_seed=.false.
      endif
      !Fred: plan to time and locate SN via stellar mass, add variable
      !      stellar mass, routine to accrete stellar mass function of
      !      ISM density, explode SN based on stellar mass distribution
      !if (lSN_mass_rate) then
      !  call set_interval(f,t_interval(SNI),l_SNI)
      !endif
      call random_number_wrapper(franSN)
!
!  Time interval follows Poisson process with rate 1/interval_SNI
!
      scaled_interval=-log(franSN)*t_interval(SNI)
!
      t_next_SNI=t+scaled_interval
      if (lroot) print"(1x,'set_next_SNI: Next SNI at time =',e15.8)",t_next_SNI
!
    endsubroutine set_next_SNI
!*****************************************************************************
    subroutine set_next_SNII(f,scaled_interval)
!
      use General, only:  random_seed_wrapper, random_number_wrapper
      use Mpicomm, only: mpiallreduce_sum
      use Grid, only: get_dVol
!
      real, dimension(mx,my,mz,mfarray) :: f
      real :: franSN, rhom, scaled_interval, tmp_interval
!
      intent(in) :: f
      intent(out) :: scaled_interval

      real, dimension(nx) :: dV
      integer :: m,n
!
!  Pre-determine time for next SNII
!  Check_SNII has a selfregulating random rate governed by the parameters of
!  cloud mass and cloud temperature, but this is very hard to regulate to test
!  different regimes, so this acts as a contraint on the rate.
!
      if (lroot) print "(1x,'check_SNII: Old t_next_SNII =',e15.8)",t_next_SNII
      if (lreset_ism_seed) then
        seed=seed_reset
        call random_seed_wrapper(PUT=seed)
        lreset_ism_seed=.false.
      endif
      call random_number_wrapper(franSN)
!
!  fred: the SN rate can be varied to reduce with mass outflows and increase
!        with inflows to regulate the mass in the disk
!        a high index iSNdx can be used to increase the sensitivity, but once
!        net mass loss occurs the index is set to 1, to avoid overheating
!
      rhom=0.
      if (lscale_SN_interval) then
        if (ldensity_nolog) then
          if (lcart_equi) then
            rhom=sum(f(l1:l2,m1:m2,n1:n2,irho))*dVol_glob
          else
            !$omp target if(loffload) map(from:rhom) has_device_addr(f) !globals: irho
            !$omp teams distribute parallel do collapse(2) private(dV) reduction(+:rhom)
            do n=n1,n2; do m=m1,m2
              call get_dVol(m,n,dV)
              rhom=rhom+sum(f(l1:l2,m,n,irho)*dV)
            enddo; enddo
            !$omp end teams distribute parallel do
            !$omp end target
          endif
        else
          if (lcart_equi) then
            rhom=sum(exp(f(l1:l2,m1:m2,n1:n2,ilnrho)))*dVol_glob
          else
            !$omp target if(loffload) map(from:rhom) has_device_addr(f)   ! globals: ilnrho
            !$omp teams distribute parallel do collapse(2) private(dV) reduction(+:rhom)
            do n=n1,n2; do m=m1,m2
              call get_dVol(m,n,dV)
              rhom=rhom+sum(exp(f(l1:l2,m,n,ilnrho))*dV)
            enddo; enddo
            !$omp end teams distribute parallel do
            !$omp end target
          endif
        endif
        call mpiallreduce_sum(rhom/box_volume,rhom)
        tmp_interval=t_interval(SNII)*(SN_interval_rhom/rhom)**iSNdx
        old_rhom=rhom !not used FG: to be removed
      else
        tmp_interval=t_interval(SNII)
      endif
      !Fred: plan to time and locate SN via stellar mass, add variable
      !      stellar mass, routine to accrete stellar mass function of
      !      ISM density, explode SN based on stellar mass distribution
      !if (lSN_mass_rate) then
      !  call set_interval(f,t_interval(SNII),l_SNI)
      !endif
!
!  Time interval follows Poisson process with rate 1/interval_SNII
!
      scaled_interval=-log(franSN)*tmp_interval
!
      t_next_SNII=t+scaled_interval
      if (lroot) print"(1x,'check_SNII: Next SNII at time =',e15.8)",t_next_SNII
!
    endsubroutine set_next_SNII
!*****************************************************************************
    subroutine set_interval(f,t_interval,l_SNI)
!
      use General, only: find_index_range
      use Grid, only: get_dVol
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension(mx,my,mz,mfarray) :: f
      real :: t_interval
      logical :: l_SNI
!
      intent(IN) :: f, l_SNI
      intent(OUT) :: t_interval
!
      real :: surface_massII, disk_massII
      integer :: m,n,nz1,nz2
      real, dimension(nx) :: dV
!
!  Identifier
!
      if (headtt) print*,'set_interval: ENTER'
!
!  Adapt expected time interval until next SN depending on the ISM mass
!  within 2x h_SN of midplane
!
!  SNII rate=5.E-12 mass(H1+HII)/solar_mass
!  van den Bergh/Tammann Annu. Rev Astron. Astrophys. 1991 29:363-407
!  SNI rate=4.7E-14 mass(H1+HII)/solar_mass + 0.35 x SNII rate
!  Mannucci et al A&A 433, 807-814 (2005)
!
      call find_index_range(z,mz,-2.*h_SNII,2.*h_SNII,nz1,nz2)
      disk_massII=0.
      if (nz1<=nz2) then
!
        if (ldensity_nolog) then
          if (lcart_equi) then
            disk_massII=sum(f(l1:l2,m1:m2,nz1:nz2,irho))*dVol_glob
          else
            !$omp target if(loffload) map(from: disk_massII) has_device_addr(f)  ! globals: irho
            !$omp teams distribute parallel do collapse(2) reduction(+:disk_massII)
            do n=nz1,nz2; do m=m1,m2
              call get_dVol(m,n,dV)
              disk_massII=disk_massII+sum(f(l1:l2,m,n,irho)*dV)
            enddo; enddo
            !$omp end teams distribute parallel do
            !$omp end target
          endif
        else
          if (lcart_equi) then
            disk_massII=sum(exp(f(l1:l2,m1:m2,nz1:nz2,ilnrho)))*dVol_glob
          else
            !$omp target if(loffload) map(from: disk_massII) has_device_addr(f) ! globals: ilnrho
            !$omp teams distribute parallel do collapse(2) private(dV) reduction(+:disk_massII)
            do n=nz1,nz2; do m=m1,m2
              call get_dVol(m,n,dV)
              disk_massII=disk_massII+sum(exp(f(l1:l2,m,n,ilnrho))*dV)
            enddo; enddo
            !$omp end teams distribute parallel do
            !$omp end target
          endif
        endif
!
      endif
!
      call mpiallreduce_sum(disk_massII,surface_massII)
!
      if (l_SNI) then
!        t_interval=solar_mass/(SNI_mass_rate+0.35*SNII_mass_rate)/surface_massII/mu
        t_interval=7.5*solar_mass/SNII_mass_rate/surface_massII/mu
        if (lroot.and.ip==1963) print"(1x,'set_interval: expected interval for SNI =',e11.4)",t_interval
      else
        t_interval=solar_mass/surface_massII/SNII_mass_rate/mu
        if (lroot.and.ip==1963) print"(1x,'set_interval: expected interval for SNII =',e11.4)",t_interval
      endif
!
    endsubroutine set_interval
!*****************************************************************************
    subroutine check_SNII(f,l_SNI)
!
!  Check for SNII, via self-regulating scheme.
!
!  03-feb-10/fred: Tested and working correctly.
!
      use General, only: random_seed_wrapper,  random_number_wrapper
      use Mpicomm, only: mpiallreduce_sum, mpibcast_real
      use EquationOfState, only: eoscalc
      use Grid, only: get_dVol
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx) :: rho, lnTT, dV
      real :: cloud_mass, cloud_mass_dim, freq_SNII, prob_SNII
      real :: franSN, fmpi1
      real, dimension(ncpus) :: cloud_mass_byproc
      integer :: icpu, m, n, iSNR, ierr
      logical :: l_SNI
      real :: dtsn
!
      intent(inout) :: f,l_SNI
!
!  Identifier
!
      if (headtt.and.ip==1963) print*,'check_SNII: ENTER'
!
      if (l_SNI) return         ! Only do if no SNI this step.
!
!  13-jul-15/fred:
!  Location by mass was found to lose too much energy due to the numerically
!  necessarily high thermal conductivity coefficients. SNe in OBs mainly explode
!  into diffuse bubbles left by their neighbours anyway. This routine left for
!  reference and possible later applications for clustering and feedback through
!  star formation, which to date has been neglected
!
      iSNR=get_free_SNR()
!
!  Determine and sum all cells comprising dense cooler clouds where type II
!  SNe are prevalent. Only check if t_next_SNII exceeded (see set_next_SNII).
!
      if (t >= t_next_SNII) then
        cloud_mass=0.0
!
!  Calculate the total mass in locations where the temperature is below
!  cloud_TT and the density is above cloud_rho, i.e. cold and dense.
!
        !$omp target if(loffload) map(from: cloud_mass) & !needs: cloud_rho, lncloud_TT 
        !$omp        has_device_addr(f) ! globals: irho, ilnrho, iss, ilnTT, FG: ldensity_nolog
        !$omp teams distribute parallel do collapse(2) private(rho,lnTT,dV) reduction(+:cloud_mass)
        do n=n1,n2
        do m=m1,m2
          call get_dVol(m,n,dV)
          if (ldensity_nolog) then
            rho=f(l1:l2,m,n,irho)
          else
            rho=exp(f(l1:l2,m,n,ilnrho))
          endif
          if (ltemperature.or.(lentropy.and.ilnTT/=0)) then
            lnTT=f(l1:l2,m,n,ilnTT)
          elseif (lentropy) then
            call eoscalc(irho_ss,rho,f(l1:l2,m,n,iss),lnTT=lnTT)
          endif
!
!  Multiply by volume element dV to find total mass.
!
          cloud_mass=cloud_mass+sum(rho*dV,mask=(rho >= cloud_rho .and. lnTT <= lncloud_TT))
        enddo
        enddo
        !$omp end teams distribute parallel do
        !$omp end target
!
!  Sum the total over all processors to find total mass.
!
        call mpiallreduce_sum(cloud_mass,cloud_mass_dim)
!
        if (ip==1963) print*,'check_SNII: cloud_mass,it,iproc=',cloud_mass,it,iproc
!
        if (lroot .and. ip ==1963) print*,'check_SNII: cloud_mass_dim,fsum(1):', cloud_mass_dim
!
!  Additional contraint on the interval between SNII events. The total time
!  elapsed since last SNII is dtsn. Probability of next event increases with
!  time and availabilty of cold dense cloud material. Calculate probability.
!  (SNI distribution is independent of mass distribution.)
!  This probability is close to 1 for solar neighbourhood rate so this check
!  could be discarded, but worth keeping as may be interesting tool.
!
        dtsn=t-last_SN_t
        freq_SNII=frac_heavy*frac_converted*cloud_mass_dim/mass_SN_progenitor/cloud_tau
        prob_SNII=freq_SNII*dtsn*5.
        if (lreset_ism_seed) then
          seed=seed_reset
          call random_seed_wrapper(PUT=seed)
          lreset_ism_seed=.false.
        endif
        call random_number_wrapper(franSN)
!
        if (lroot.and.ip==1963) then
          if (cloud_mass_dim>0.0.and.franSN<=2.0*prob_SNII) then
            print*,'check_SNII: freq,prob,rnd,dtsn:',freq_SNII,prob_SNII,franSN,dtsn
            print*,'check_SNII: frac_heavy,frac_converted,cloud_mass_dim,mass_SN,cloud_tau', &
                                frac_heavy,frac_converted,cloud_mass_dim,mass_SN,cloud_tau
          endif
        endif
!
!  If likelihood of SNII greater than random number locate SNII.
!
        if (franSN <= prob_SNII) then
!
!  The position_SN_bycloudmass needs the cloud_masses for each processor.
!  Communicate and store them here to avoid recalculation.
!
          cloud_mass_byproc=0.
!
!  Use non-root broadcasts for the communication...
!
          do icpu=1,ncpus
            fmpi1=cloud_mass
            call mpibcast_real(fmpi1,icpu-1)
            cloud_mass_byproc(icpu)=fmpi1
          enddo
!
!  Locate the next explosion.
!
          if (lroot.and.ip==1963) print*,'check_SNII: cloud_mass_byproc:',cloud_mass_byproc
          call position_SN_bycloudmass(f,cloud_mass_byproc,SNRs(iSNR),preSN,ierr)
!
!  If location too hot reset ierr and return to program.
!
          if (ierr == iEXPLOSION_TOO_HOT) then
            call free_SNR(iSNR)
            ierr=iEXPLOSION_OK
            return
          endif
!
!  Try to explode SNII and if successful reset time of most recent (last_SN_t)
!  and next (t_next_SNII) explosion.
!
          SNRs(iSNR)%feat%t=t
          SNRs(iSNR)%indx%SN_type=2
          call explode_SN(f,SNRs(iSNR),ierr,preSN)
          if (ierr==iEXPLOSION_OK) then
            if (lSN_mass_rate) call set_interval(f,t_interval(SNII),l_SNI)
            last_SN_t=t
          endif
!
        endif
      endif
!
!  If returned unexploded stop code running out of free slots & reset ierr.
!
      ierr=iEXPLOSION_OK
      call free_SNR(iSNR)
      l_SNI=.true.
!
    endsubroutine check_SNII
!*****************************************************************************
    subroutine position_SN_testposition(f,SNR)
!
!  Determine position for next SN (w/ fixed scale-height).
!
      use General, only: find_proc

      real, intent(in), dimension(mx,my,mz,mfarray) :: f
      type (SNRemnant), intent(inout) :: SNR
!
      real :: z00, x00, y00
      integer :: i
!
      if (headtt) print*,'position_SN_testposition: ENTER'
!
!  Calculate the global (nzgrid) lower z-coordinate.
!
      if (lperi(1)) then; x00=xyz0(1)+.5*dx; else; x00=xyz0(1); endif
      if (lperi(2)) then; y00=xyz0(2)+.5*dy; else; y00=xyz0(2); endif
      if (lperi(3)) then; z00=xyz0(3)+.5*dz; else; z00=xyz0(3); endif
!
!  Pick SN position (SNR%indx%l,SNR%indx%m,SNR%indx%n).
!
      if (lroot) then
        if (center_SN_x==impossible) then
          i=max(int(nxgrid/2)+1,1)
        else
          i=int((center_SN_x-x00)/dx)+1
        endif
        SNR%indx%ipx=(i-1)/nx ! uses integer division
        SNR%indx%l=i-(SNR%indx%ipx*nx)+nghost
!
        if (center_SN_y==impossible) then
          i=max(int(nygrid/2)+1,1)
        else
          i=int((center_SN_y-y00)/dy)+1
        endif
        SNR%indx%ipy=(i-1)/ny ! uses integer division
        SNR%indx%m=i-(SNR%indx%ipy*ny)+nghost
!
        if (center_SN_z==impossible) then
          i=max(int(nzgrid/2)+1,1)
        else
          i=int((center_SN_z-z00)/dz)+1
        endif
        SNR%indx%ipz=(i-1)/nz   ! uses integer division
        SNR%indx%n=i-(SNR%indx%ipz*nz)+nghost
        SNR%indx%iproc=find_proc(SNR%indx%ipx,SNR%indx%ipy,SNR%indx%ipz)
      endif
      call share_SN_parameters(f,SNR)
!
    endsubroutine position_SN_testposition
!*****************************************************************************
    subroutine position_SN_gaussianz(f,h_SN,SNR)
!
!  Determine position for next SN (w/ fixed scale-height).
!  15-jul-14/luiz+fred: added SN horizontal clustering algorithm for SNII
!                       70% probability SNII located within 300pc of previous
!  27-oct-16/fred: z-location revised to use non-equidistant grid needs z from
!                  each processor - use mpi on all procs not just root
!
      use General, only: random_number_wrapper, random_seed_wrapper, find_proc
      use Mpicomm, only: mpiallreduce_max, mpireduce_min, mpireduce_max, &
                         mpireduce_sum, mpibcast_real, mpigather_z_1D
      use Grid, only: get_dVol
!
      real, dimension(nx) :: dV
!
      real, intent(in), dimension(mx,my,mz,mfarray) :: f
      real, intent(in) :: h_SN
      type (SNRemnant), intent(inout) :: SNR
!
!  parameters required to determine the vertical centre of mass of the disk
!
      real, dimension(nzgrid) :: rhotot
      real, dimension(nz) :: rhosum, mpisum
      real :: rhomax, hSN
      real, dimension(2) :: mpizspan
      integer :: lm_range, previous_SNl, previous_SNm, previous_SNn
!
!  parameters for random location of SN - about zdisk
!
      real, dimension(nzgrid) :: cum_prob_SN
      real, dimension(3) :: fran3
      real, dimension(4) :: mpicluster
      logical :: lgauss, lnew_cluster
      integer :: i,j,k,m,n,indrhomax, nzskip=10 !prevent SN from being too close to boundaries
!
      if (headtt.and.ip==1963) print*,'position_SN_gaussianz: ENTER'
!
!  The disk oscillates. To keep the SN random dist centered about centre of mass find
!  zdisk where the peak mean density(z) resides and shift gaussian up/down
!
      Get_zdisk: if (lfirst_zdisk) then
!
!  sum the mass on each processor
!
        rhosum=0.0
        if (.not.lcart_equi) then
          !$omp target if(loffload) map(from: rhosum) has_device_addr(f)   ! globals: irho, ilnrho, ldensity_nolog
          !$omp teams distribute parallel do collapse(2) private(dV) reduction(+:rhosum)
          do n=n1,n2; do m=m1,m2
            call get_dVol(m,n,dV)
            if (ldensity_nolog) then
              rhosum(n-n1+1)=rhosum(n-n1+1)+sum(f(l1:l2,m,n,irho)*dV)
            else
              rhosum(n-n1+1)=rhosum(n-n1+1)+sum(exp(f(l1:l2,m,n,ilnrho))*dV)
            endif
          enddo; enddo
          !$omp end teams distribute parallel do
          !$omp end target
        else
          !$omp target if(loffload) map(from: rhosum) has_device_addr(f)   ! globals: irho, ilnrho, ldensity_nolog, FG: dVol_glob
          !$omp teams distribute parallel do
          do n=n1,n2
            if (ldensity_nolog) then
              rhosum(n-n1+1)=sum(f(l1:l2,m1:m2,n,irho))*dVol_glob
            else
              rhosum(n-n1+1)=sum(exp(f(l1:l2,m1:m2,n,ilnrho)))*dVol_glob
            endif
          enddo
          !$omp end teams distribute parallel do
          !$omp end target
        endif
!
!  broadcast the mass on xy-plane and then collect rhosum on first processor
!  in each plane.
!
!  mpi sum rhosum(iz) on each horizontal slice to mpisum(iz) then
!  collect each mpisum into rhotot along z to complete array over nzgrid
!
        call mpireduce_sum(rhosum,mpisum,nz,idir=12)
        if (lfirst_proc_xy) then
          call mpigather_z_1D(mpisum,rhotot)
          if (lroot) then
            if (mod(it,2)==0.and..not.lstart) then ! for reproducibility: it=1 -> it=0 initially
              indrhomax=maxloc(rhotot,dim=1); rhomax=rhotot(indrhomax)
            else
              indrhomax=maxloc(rhotot(nzgrid:1:-1),dim=1); rhomax=rhotot(indrhomax)
            endif
            mpizspan(1)=zgrid(indrhomax)
            if (lh_SNII_adjust) mpizspan(2) = rhomax*dz_1(indrhomax)/(Lxyz(1)*Lxyz(2))
          endif
        endif
        call mpibcast_real(mpizspan,2)
        zdisk=mpizspan(1)
        if (lh_SNII_adjust) maxrho=mpizspan(2)
      endif Get_zdisk
!
      if (lroot.and.ip==1963) print"(1x,'position_SN_gaussianz: zdisk =',f8.4)",zdisk
      if (lh_SNII_adjust .and. h_SN==h_SNII) then
        hSN = h_SN * cloud_rho/maxrho
        if (lroot.and.ip==1963) then
          print "(1x,'position_SN_gaussianz: hSN vs h_SNII =',2e11.4)",hSN,h_SN
          print "(1x,'position_SN_gaussianz: maxrho =',e11.4)",maxrho
          print "(1x,'position_SN_gaussianz: cloud_rho =',e11.4)",cloud_rho
        endif
      else
        hSN = h_SN
      endif
!
!  Pick SN position (SNR%indx%l,SNR%indx%m,SNR%indx%n).
!
!  Get 3 random numbers on all processors to keep rnd. generators in sync.
!
      if (lreset_ism_seed) then
        seed=seed_reset
        call random_seed_wrapper(PUT=seed)
        lreset_ism_seed=.false.
      endif
      call random_number_wrapper(fran3)
!
!  13-jul-15/fred: NB need to revisit OB clustering x,y not updated or time
!  constrained. May need to include z also
!
      Get_z: if (lroot) then
        if (lOB_cluster .and. h_SN==h_SNII) then
!  If OB clustering for SNII, while within time span of current cluster
          if (t < t_cluster) then ! still using current cluster coords
            if (ip==1963) print "(1x,'position_SN_gaussianz: cluster lifetime until',e11.4)",t_cluster
            previous_SNl = int(( x_cluster - xyz0(1) )/Lx)*nxgrid +1
            previous_SNm = int(( y_cluster - xyz0(2) )/Ly)*nygrid +1
            previous_SNn = int(( z_cluster - xyz0(3) )/Lz)*nzgrid +1
            lm_range = 2*int(SN_clustering_radius*nxgrid/Lx)
            if (fran3(1) < p_OB) then ! checks whether the SN is in a cluster
              if (ip==1963) print "(1x,'position_SN_gaussianz: in cluster x, y, z =',3e11.4)", &
                                             x_cluster,y_cluster,z_cluster
              i=int(fran3(1)*lm_range/p_OB)+previous_SNl+1
              j=int(fran3(2)*lm_range/p_OB)+previous_SNm+1
              k=int(fran3(3)*lm_range/p_OB)+previous_SNn+1
              if (k>nxgrid) k=k-nz
              SNR%indx%ipz=(k-1)/nz  ! uses integer division
              SNR%indx%n=k-(SNR%indx%ipz*nz)+nghost
              lgauss = .true.
              lnew_cluster = .false.
            else ! outside cluster
              i=int(fran3(1)*(nxgrid-lm_range)/(1.0-p_OB))+previous_SNl+1
              j=int(fran3(2)*(nygrid-lm_range)/(1.0-p_OB))+previous_SNm+1
              lgauss = .true.
              lnew_cluster = .false.
            endif
          else
!  If OB clustering for SNII, time to set new cluster location and duration
            i=int(fran3(1)*nxgrid)+1
            j=int(fran3(2)*nygrid)+1
            lgauss = .true.
            lnew_cluster = .true.
            call set_next_OB(t_interval_OB)
          endif
        else ! clustering not used
          i=int(fran3(1)*nxgrid)+1
          j=int(fran3(2)*nygrid)+1
          lgauss = .true.
          lnew_cluster = .false.
        endif
!
!  Cumulative probability function in z calculated each time for moving zdisk.
!
        if (lgauss) then
          cum_prob_SN=0.0
          do k=1,nzgrid
            if (k<nzskip) then
              cum_prob_SN(k)=0.0
            elseif (k>nzgrid-nzskip) then
              cum_prob_SN(k)=cum_prob_SN(k-1)
            else
              cum_prob_SN(k)=cum_prob_SN(k-1)+exp(-0.5*((k * Lz/nzgrid + xyz0(3)-zdisk)/hSN)**2)
            endif
          enddo
          cum_prob_SN = cum_prob_SN / max(cum_prob_SN(nzgrid-nzskip), tini)
!
!  The following should never be needed, but just in case floating point
!  errors ever lead to cum_prob_SNI(nzgrid-nzskip) < rnd < 1.
!
          cum_prob_SN(nzgrid-nzskip+1:nzgrid)=1.0
!
          do k=nzskip+1,nzgrid-nzskip
            if (cum_prob_SN(k-1)<=fran3(3) .and. fran3(3)<cum_prob_SN(k)) then
              SNR%indx%ipz=(k-1)/nz  ! uses integer division
              SNR%indx%n=k-(SNR%indx%ipz*nz)+nghost
              exit
            endif
          enddo
        endif
!
        if (i>nxgrid) i=i-nxgrid
        SNR%indx%ipx=(i-1)/nx  ! uses integer division
        SNR%indx%l=i-(SNR%indx%ipx*nx)+nghost
        if (j>nygrid) j=j-nygrid
        SNR%indx%ipy=(j-1)/ny  ! uses integer division
        SNR%indx%m=j-(SNR%indx%ipy*ny)+nghost
        SNR%indx%iproc=find_proc(SNR%indx%ipx,SNR%indx%ipy,SNR%indx%ipz)
        if (lnew_cluster) then
          x_cluster = (SNR%indx%l-1) * Lx/nxgrid + xyz0(1)
          y_cluster = (SNR%indx%m-1) * Ly/nygrid + xyz0(2)
          z_cluster = zdisk
        endif
      endif Get_z

      mpicluster = (/x_cluster,y_cluster,z_cluster,t_cluster/)
      call mpibcast_real(mpicluster,4)
      x_cluster=mpicluster(1)
      y_cluster=mpicluster(2)
      z_cluster=mpicluster(3)
      t_cluster=mpicluster(4)
!
      call share_SN_parameters(f,SNR)
!
    endsubroutine position_SN_gaussianz
!*****************************************************************************
    subroutine position_SN_uniformz(f,SNR)
!
!  Determine position for next SN (w/ fixed scale-height).
!
      use General, only: random_seed_wrapper, random_number_wrapper, find_proc
!
      real, intent(in), dimension(mx,my,mz,mfarray) :: f
      type (SNRemnant), intent(inout) :: SNR
!
      real, dimension(3) :: fran3
      integer :: i   !prevent SN from being too close to boundaries
!
      if (headtt) print*,'position_SN_uniformz: ENTER'
!
!  Pick SN position (SNR%indx%l,SNR%indx%m,SNR%indx%n).
!
      if (lreset_ism_seed) then
        seed=seed_reset
        call random_seed_wrapper(PUT=seed)
        lreset_ism_seed=.false.
      endif
      call random_number_wrapper(fran3)
!
!  Get 3 random numbers on all processors to keep rnd. generators in sync.
!
      if (lroot) then
        i=int(fran3(1)*nxgrid)+1
        if (nxgrid==1) i=1
        SNR%indx%ipx=(i-1)/nx ! uses integer division
        SNR%indx%l=i-(SNR%indx%ipx*nx)+nghost
!
        i=int(fran3(2)*nygrid)+1
        if (nygrid==1) i=1
        SNR%indx%ipy=(i-1)/ny  ! uses integer division
        SNR%indx%m=i-(SNR%indx%ipy*ny)+nghost
!
        i=int(fran3(3)*nzgrid)+1
        if (nzgrid==1) i=1
        SNR%indx%ipz=(i-1)/nz   ! uses integer division
        SNR%indx%n=i-(SNR%indx%ipz*nz)+nghost
        SNR%indx%iproc=find_proc(SNR%indx%ipx,SNR%indx%ipy,SNR%indx%ipz)
      endif
!
      call share_SN_parameters(f,SNR)
!
    endsubroutine position_SN_uniformz
!*****************************************************************************
    subroutine position_SN_bycloudmass(f,cloud_mass_byproc,SNR,preSN,ierr)
!
!  Determine position for next SNII (using Boris' scheme). It seems impractical
!  to sort all high density points across all processors; instead, we just
!  construct cumulative pdfs that allow us to pick a processor, and then a
!  point on that processor, with probability proportional to rho. As a result,
!  the SN position is *not* independent of ncpus (nor of nprocy and nprocz).
!  It is repeatable given fixed nprocy/z though.
!
      use General, only: random_seed_wrapper,  random_number_wrapper
      use Grid, only: get_dVol
      use EquationOfState, only: eoscalc
      use Mpicomm, only: mpibcast_int, mpibcast_real
      use Grid, only: get_dVol
      !$ use omp_lib
!
      real, intent(in), dimension(mx,my,mz,mfarray) :: f
      real, intent(in) , dimension(ncpus) :: cloud_mass_byproc
      type (SNRemnant), intent(inout) :: SNR
      integer, intent(in), dimension(4,npreSN)::preSN
      integer :: ierr

      real, dimension(0:ncpus) :: cum_prob_byproc
      integer, dimension(4) :: tmpsite
      real :: franSN, cloud_mass,cum_mass,cloud_mass_proc
      real, dimension(nx) :: dV,rho,lnTT
      integer :: icpu,l,m,n, ipsn
      logical :: lfound
!
!  identifier
!
      if (lroot.and.ip==1963) print*,'position_SN_bycloudmass: ENTER'
!
!  Construct cumulative distribution function, using cloud_mass_byproc.
!  NB: icpu=iproc+1 (iproc in [0,ncpus-1], icpu in [1,ncpus] ).
!
      ierr=iEXPLOSION_OK
      cloud_mass=0.
      cum_prob_byproc=0.
      do icpu=1,ncpus
        cloud_mass=cloud_mass+cloud_mass_byproc(icpu)
        cum_prob_byproc(icpu)=cum_prob_byproc(icpu-1)+cloud_mass_byproc(icpu)
      enddo
      cum_prob_byproc=cum_prob_byproc/cum_prob_byproc(ncpus)
      if (lroot.and.ip==1963) then
        print*,'position_SN_bycloudmass: cloud_mass_byproc=',cloud_mass_byproc
        print*,'position_SN_bycloudmass: cum_prob_byproc=',cum_prob_byproc
        print*,'position_SN_bycloudmass: cloud_mass=',cloud_mass
      endif
!
!  Use random number to detemine which processor SN is on.
!  (Use root processor for rand to ensure repeatability.)
!  All processors must be at the same rand seed or MPI will hang so above line incorrect
!
      if (lreset_ism_seed) then
        seed=seed_reset
        call random_seed_wrapper(PUT=seed)
        lreset_ism_seed=.false.
      endif
      call random_number_wrapper(franSN)
      do icpu=1,ncpus
        if (cum_prob_byproc(icpu-1)<=franSN .and. franSN < cum_prob_byproc(icpu)) then
          SNR%indx%iproc=icpu-1
          exit
        endif
      enddo
      if (lroot.and.ip==1963) &
          print*, 'position_SN_bycloudmass: franSN,SNR%indx%iproc=',franSN,SNR%indx%iproc
!
!  Use random number to pick SNII location on the right processor.
!  (No obvious reason to re-use the original random number for this.)
!
      call random_number_wrapper(franSN)
      if (iproc==SNR%indx%iproc) then
        cum_mass=0.
        cloud_mass_proc=cloud_mass_byproc(SNR%indx%iproc+1)
        lfound = .false.
        !$omp target if(loffload) map(from: ierr,cum_mass) map(tofrom: SNR) & !needs: cloud_rho,lncloud_TT,preSN) &
        !$omp        has_device_addr(f)  !globals: ip, irho,ilnrho,ilnTT,iss,ldensity_nolog
        !$omp teams distribute parallel do collapse(2) private(dV,rho,lnTT,l,ipsn) reduction(+:cum_mass)  !!!reduction unclear
mn_loop:do n=n1,n2
        do m=m1,m2
          !$ if (lfound) cycle
          if (ldensity_nolog) then
            rho=f(l1:l2,m,n,irho)
          else
            rho=exp(f(l1:l2,m,n,ilnrho))
          endif
          if (ilnTT/=0) then
            lnTT=f(l1:l2,m,n,ilnTT)
          elseif (lentropy) then
            call eoscalc(irho_ss,rho,f(l1:l2,m,n,iss), lnTT=lnTT)
          endif
          call get_dVol(m,n,dV)
          do l=1,nx
            if (rho(l)>=cloud_rho .and. lnTT(l)<=lncloud_TT) then
              cum_mass=cum_mass+rho(l)*dV(l)
              if (franSN <= cum_mass/cloud_mass_proc) then
! hit -> leave the loop
                !$omp critical
                SNR%indx%l=l+l1-1; SNR%indx%m=m; SNR%indx%n=n
                !$omp end critical
!
!  Check that the same site is not being used repeatedly. If used recently
!  skip and get a new random number next time step.
!
                do ipsn=1,npreSN
                  if ((SNR%indx%l==preSN(1,ipsn)) .and. &
                      (SNR%indx%m==preSN(2,ipsn)) .and. &
                      (SNR%indx%n==preSN(3,ipsn)) .and. &
                      (SNR%indx%iproc==preSN(4,ipsn))) then
                    !$omp atomic
                    ierr=iEXPLOSION_TOO_HOT
                    if (.not.lgpu.and.ip==1963) &
                      print*,'position_by_cloudmass: iEXPLOSION_TOO_HOT, iproc,it,preSN=',iproc,it,preSN(:,ipsn)
                  endif
                enddo
                !!DEC$ if OPENMP
                !$ lfound = .true.
                !!DEC$ else
                exit mn_loop
                !!DEC$ endif
              endif
            endif
          enddo
        enddo
        enddo mn_loop   !!!MR: whatif lfound=.false.? FG: cum_mass in (0,1) as is franSN(1) but needs to be checked before use
        !$omp end teams distribute parallel do
        !$omp end target
        if (ip==1963) then
          tmpsite=(/SNR%indx%l,SNR%indx%m,SNR%indx%n,SNR%indx%iproc/)
          print*, 'position_by_cloudmass: preSN,iproc,it =',preSN,iproc,it
          print*, 'position_SN_bycloudmass: l,m,n,SNRproc =',tmpsite
          if (lroot) print*, 'position_SN_bycloudmass:cum_mass,cum_prob_onproc,franSN=', &
                             cum_mass,cum_mass/cloud_mass_proc,franSN
        endif
      endif
!
      call mpibcast_int(ierr,SNR%indx%iproc)
      if (ierr==iEXPLOSION_TOO_HOT) then
        if (lroot.and.ip==1963) print*,'position_SN_bycloudmass: iEXPLOSION_TOO_HOT,ierr',ierr
        return
      endif
!
      call mpibcast_int(tmpsite,4,SNR%indx%iproc)
      SNR%indx%l=tmpsite(1);SNR%indx%m=tmpsite(2)
      SNR%indx%n=tmpsite(3);SNR%indx%iproc=tmpsite(4)
      if (ip==1963) print*,'position_SN_bycloudmass: MPI tmpsite,iproc,it =',tmpsite,iproc,it
      call share_SN_parameters(f,SNR)
      if (ip==1963) print*,'position_SN_bycloudmass: SN_param,iproc,it =', &
          SNR%indx%l,SNR%indx%m,SNR%indx%n,SNR%indx%iproc,iproc,it
!
!  Reset status for next explosion.
!
      ierr=iEXPLOSION_OK
!
    endsubroutine position_SN_bycloudmass
!*****************************************************************************
    subroutine share_SN_parameters(f,SNR)
!
!  Handle common SN positioning processor communications.
!
!  27-aug-2003/tony: coded
!
      use EquationOfState, only: eoscalc
      use Mpicomm, only: mpibcast_int, mpibcast_real
!
      real, intent(in), dimension(mx,my,mz,mfarray) :: f
      type (SNRemnant), intent(inout) :: SNR
!
      real, dimension(nx) :: lnTT
      real, dimension(7) :: fmpi7
      integer, dimension(4) :: impi4
      real :: sndx, sndy, sndz
!
!  Broadcast position to all processors from root; also broadcast SNR%indx%iproc,
!  needed for later broadcast of SNR%site%rho.
!
      impi4=(/ SNR%indx%iproc, SNR%indx%l, SNR%indx%m, SNR%indx%n /)
      call mpibcast_int(impi4,4)
      SNR%indx%iproc=impi4(1)
      SNR%indx%l=impi4(2)
      SNR%indx%m=impi4(3)
      SNR%indx%n=impi4(4)
!
!  With current SN scheme, we need rho at the SN location.
!
      if (iproc==SNR%indx%iproc) then
        !$omp target if(loffload) map(tofrom: SNR) has_device_addr(f)  ! globals: irho, ilnrho  !?single
        if (ldensity_nolog) then
          SNR%site%lnrho=log(f(SNR%indx%l,SNR%indx%m,SNR%indx%n,irho))
        else
          SNR%site%lnrho=f(SNR%indx%l,SNR%indx%m,SNR%indx%n,ilnrho)
        endif
        !$omp end target
        SNR%site%rho=exp(SNR%site%lnrho);
!
        m=SNR%indx%m
        n=SNR%indx%n
        !$omp target if(loffload) map(tofrom: SNR) has_device_addr(f) ! map(to: m,n)  ! globals: ilnTT   !?single
        if (leos_ionization.or.leos_temperature_ionization) then
          SNR%site%lnTT=f(SNR%indx%l+1,m,n,ilnTT)
        else
          call eoscalc(irho_ss,SNR%site%rho,f(SNR%indx%l,SNR%indx%m,SNR%indx%n,iss),lnTT=SNR%site%lnTT)
        endif
        !$omp end target
        SNR%feat%x=0.; SNR%feat%y=0.; SNR%feat%z=0.
        sndx=0.; sndy=0.; sndz=0.
        if (nxgrid/=1) then
          SNR%feat%x=x(SNR%indx%l)
          sndx=x(SNR%indx%l+1)-x(SNR%indx%l)
        endif
        if (nygrid/=1) then
          SNR%feat%y=y(SNR%indx%m)
          sndy=y(SNR%indx%m+1)-y(SNR%indx%m)
        endif
        if (nzgrid/=1) then
          SNR%feat%z=z(SNR%indx%n)
          sndz=z(SNR%indx%n+1)-z(SNR%indx%n)
        endif
        SNR%feat%dr=max( sndx,sndy,sndz )
        if (center_SN_x/=impossible) SNR%feat%x=center_SN_x
        if (center_SN_y/=impossible) SNR%feat%y=center_SN_y
        if (center_SN_z/=impossible) SNR%feat%z=center_SN_z
!
!  10-Jun-10/fred:
!  Adjust radius according to density of explosion site to concentrate energy
!  when locations are dense.
!
        SNR%feat%radius=width_SN
        if (lSN_scale_rad) then
          SNR%feat%radius=(0.75*solar_mass/SNR%site%rho*pi_1*N_mass)**onethird
          SNR%feat%radius=min(SNR%feat%radius,Lxyzmin)
          SNR%feat%radius=max(SNR%feat%radius,rfactor_SN*SNR%feat%dr) ! minimum grid resolution
        endif
!
!  Better initialise these to something on the other processors
!
      else
        SNR%site%lnrho=0.
        SNR%site%lnTT=0.
        SNR%feat%x=0.
        SNR%feat%y=0.
        SNR%feat%z=0.
        SNR%feat%radius=0.
        SNR%feat%dr=0.
      endif
!
!    Broadcast to all processors.
!
      fmpi7=(/ SNR%feat%x, SNR%feat%y, SNR%feat%z, SNR%site%lnrho, SNR%site%lnTT, SNR%feat%radius, SNR%feat%dr /)
      call mpibcast_real(fmpi7,7,SNR%indx%iproc)
!
      SNR%feat%x=fmpi7(1); SNR%feat%y=fmpi7(2); SNR%feat%z=fmpi7(3);
      SNR%site%lnrho=fmpi7(4); SNR%site%lnTT=fmpi7(5)
      SNR%feat%radius=fmpi7(6); SNR%feat%dr=fmpi7(7)
!
      SNR%site%rho=exp(SNR%site%lnrho);
!
      call eoscalc(ilnrho_lnTT,SNR%site%lnrho,SNR%site%lnTT,yH=SNR%site%yH,ss=SNR%site%ss,ee=SNR%site%ee)
      SNR%site%TT=exp(SNR%site%lnTT)
!
      if (lroot.and.ip==1963) then
        print "(1x,'share_SN_parameters: iproc, l, m, n =',i8,3i6)", &
              SNR%indx%iproc,SNR%indx%l,SNR%indx%m,SNR%indx%n
        print "(1x,'share_SN_parameters: x_SN, y_SN, z_SN =',3f8.4)", SNR%feat%x,SNR%feat%y,SNR%feat%z
        print "(1x,'share_SN_parameters: rho, ss, TT =',3e11.4)", SNR%site%rho,SNR%site%ss,SNR%site%TT
        print "(1x,'share_SN_parameters: SN radius, SN dr =',2f8.5)", SNR%feat%radius,SNR%feat%dr
      endif
!
    endsubroutine share_SN_parameters
!*****************************************************************************
    subroutine explode_SN(f,SNR,ierr,preSN)
!
!  Implement SN (of either type), at pre-calculated position.
!
!  ??-nov-02/grs : coded from GalaxyCode
!  20-may-03/tony: pencil formulation and broken into subroutines
!
      use EquationOfState, only: eoscalc
      use Mpicomm, only: mpiallreduce_max, mpiallreduce_sum
      use General, only: keep_compiler_quiet
      use Grid, only: get_dVol
!
      real, intent(inout), dimension(mx,my,mz,mfarray) :: f
      type (SNRemnant), intent(inout) :: SNR
      integer, intent(inout), optional, dimension(4,npreSN) :: preSN
      integer, intent(inout), optional :: ierr
!
      real :: c_SN,cmass_SN,cvelocity_SN,ecr_SN,cmass_tmp,c_SNmax
      real :: width_energy, width_mass, width_velocity
      real :: rhom, rhomin, ekintot, radius_best
      real :: ekintot_new, ambient_mass
      real :: Nsol_ratio, Nsol_ratio_best, radius_min, radius_max, sol_mass_tot
      real :: uu_sedov, rad_hot, deltarho_hot
      real :: radius2, radius3, SNvol, radius2mass
!
      real, dimension(nx) :: deltarho, deltaEE, deltaCR, rho_new, dV
      real, dimension(nx,3) :: deltafcr, deltauu, outward_normal_SN
      real, dimension(3) :: dmpi2, dmpi2_tmp
      real, dimension(nx) :: yH, lnTT, rho_old, ee_old
      real :: maxlnTT, site_mass, maxTT, mmpi, etmp, ktmp, max_cmass
      real :: t_interval_SN, SNrate, frackin, RPDS, SN_TT_ratio_max
      integer :: i, mpiierr, ind_maxTT
      logical :: lfound
!
      SNR%indx%state=SNstate_exploding
      c_SN=0.;cmass_SN=0.;cvelocity_SN=0.;ecr_SN=0.;cmass_tmp=0.
!
!  Identifier
!
      if (lroot.and.ip==1963) print*,'explode_SN: SN type =',SNR%indx%SN_type
!
!  Calculate explosion site mean density.
!
      call get_properties(f,SNR,rhom,ekintot,rhomin)
      SNR%feat%rhom=rhom
!
!  Initialize SN volume
!
      sol_mass_tot=solar_mass*N_mass
      SNvol=fourthird*pi/sol_mass_tot
!
!  Rescale injection radius to contain only N_mass solar masses or at least
!  min radius. Iterate a few times to improve radius match to N_mass.
!
      if (lSN_scale_rad) then
        radius_min=rfactor_SN*SNR%feat%dr
        radius_max=Lxyzmin
        radius_best=SNR%feat%radius
        Nsol_ratio=SNvol*rhom*SNR%feat%radius**3
        Nsol_ratio_best=impossible
        if ((Nsol_ratio<1).or.(SNR%feat%radius/=radius_min)) then
          do i=1,25
            if (Nsol_ratio<1) then
              radius_min=SNR%feat%radius
            else
              radius_max=SNR%feat%radius
            endif
            SNR%feat%radius=0.25*(3*radius_min+radius_max)
            call get_properties(f,SNR,rhom,ekintot,rhomin)
            Nsol_ratio=SNvol*rhom*SNR%feat%radius**3
            if ((Nsol_ratio>=1).and.(abs(Nsol_ratio-1)<Nsol_ratio_best)) then
              Nsol_ratio_best=Nsol_ratio
              radius_best=SNR%feat%radius
            endif
            if (lroot.and.ip==1963) then
              print "(1x,'explode_SN: i          ',   i6)",i
              print "(1x,'explode_SN: radius_min ', f9.6)",radius_min
              print "(1x,'explode_SN: radius_max ', f9.6)",radius_max
              print "(1x,'explode_SN: Rmax-Rmin  ',e11.4)",radius_max-radius_min
              print "(1x,'explode_SN: radius_best', f9.6)",radius_best
              print "(1x,'explode_SN: Nsol       ',e11.4)",Nsol_ratio*sol_mass_tot/solar_mass
            endif
            if (radius_max-radius_min<SNR%feat%dr*0.01) exit
          enddo
          SNR%feat%radius=radius_best
          call get_properties(f,SNR,rhom,ekintot,rhomin)
        endif
      endif

      if (present(ierr)) then
        call get_properties(f,SNR,rhom,ekintot,rhomin,ierr)
        if (ierr==iEXPLOSION_TOO_UNEVEN.and..not.lSN_list) return
        ambient_mass=SNvol*rhom*SNR%feat%radius**3
        if (ambient_mass/sol_mass_tot<eps_mass) then
          ierr=iEXPLOSION_TOO_RARIFIED
          if (.not.lSN_list) return
        endif
        SNR%feat%rhom=rhom
        call get_properties(f,SNR,rhom,ekintot,rhomin,ierr)
      else
        call get_properties(f,SNR,rhom,ekintot,rhomin)
      endif

      SNR%feat%rhom=rhom
      if (lroot.and.ip==1963) print "(1x,'explode_SN: total old kinetic energy =',e12.5)", ekintot
!
!  Calculate effective Sedov evolution time and shell speed diagnostic.
!
      SNR%feat%t_sedov=sqrt(SNR%feat%radius**(2+dimensionality)* &
                            SNR%feat%rhom/ampl_SN/xsi_sedov*sedov_norm)
      uu_sedov=2./(2.+dimensionality)*SNR%feat%radius/SNR%feat%t_sedov
!
!  Calculate the end of the Sedov-Taylor phase (shell formation) t_SF
!  Ref Kim & Ostriker 2015 ApJ 802:99 Eq. 7
!  ref Simpson et al. 2015 ApJ 809:69 Eq. 17, 18
!
      SNR%feat%t_SF = SFt_norm/SNR%feat%rhom**(4./7)*ampl_SN**(3./14)
      RPDS=SFr_norm*ampl_SN**(2./7)/SNR%feat%rhom**(3./7)
      if (lroot.and.ip==1963) then
        print "(1x,'explode_SN: Shell forming start t_SF',e12.5)",SNR%feat%t_SF
        print "(1x,'explode_SN: Elapsed shell formation',e13.5)",SNR%feat%t_sedov-SNR%feat%t_SF
        print "(1x,'explode_SN: Shell forming radius RPDS vs radius',2f7.4)",RPDS,SNR%feat%radius
      endif
!
!  Calculate the SN kinetic energy fraction for shell formation energy
!  losses correction ref Simpson et al. 2015 ApJ 809:69 Eq. 16.
!
      etmp=eampl_SN; ktmp=kampl_SN
      if (RPDS<SNR%feat%radius.and.lSN_autofrackin) then
        if (SNR%feat%rhom>0.8*m_H_cgs/unit_density .and. SNR%feat%dr>pc_cgs/unit_length) then
          frackin = kfrac_norm*SNR%feat%rhom*RPDS**7/ampl_SN/(SNR%feat%t_SF*SNR%feat%dr)**2
          if (frackin < 1e-3) frackin = 0.
        else
          frackin = 0.
        endif
        frackin=min(kin_max,frackin)
        etmp=(1.-frackin-frac_ecr)*ampl_SN
        ktmp=frackin*ampl_SN
        if (lroot.and.ip==1963) print "(1x,'explode_SN: Reset fractions SNE frackin',f7.4)",frackin
      endif
      if (lroot.and.ip==1963) print &
          "(1x,'explode_SN: SNE fractional energy kampl_SN, eampl_SN',2e9.2)", ktmp, etmp
!
!  Adjust radial scale if different from SNR%feat%radius.
!
      width_energy  =SNR%feat%radius*energy_width_ratio
      width_mass    =SNR%feat%radius*mass_width_ratio
      width_velocity=SNR%feat%radius*velocity_width_ratio
!
!  Energy insertion normalization.
!
      if (thermal_profile=="gaussian3") then
        c_SN=etmp/(cnorm_SN(dimensionality)*width_energy**dimensionality)
        if (frac_ecr>0.) ecr_SN=campl_SN/(cnorm_SN(dimensionality)*width_energy**dimensionality)
        c_SNmax=ampl_SN/(cnorm_SN(dimensionality)*rfactor_SN*SNR%feat%dr**dimensionality)
      elseif (thermal_profile=="gaussian2") then
        c_SN=etmp/(cnorm_gaussian2_SN(dimensionality)*width_energy**dimensionality)
        if (frac_ecr>0.) &
            ecr_SN=campl_SN/(cnorm_gaussian2_SN(dimensionality)*width_energy**dimensionality)
        c_SNmax=ampl_SN/(cnorm_gaussian2_SN(dimensionality)*rfactor_SN*SNR%feat%dr**dimensionality)
      elseif (thermal_profile=="gaussian") then
        c_SN=etmp/(cnorm_gaussian_SN(dimensionality)*width_energy**dimensionality)
        if (frac_ecr>0.) &
            ecr_SN=campl_SN/(cnorm_gaussian_SN(dimensionality)*width_energy**dimensionality)
        c_SNmax=ampl_SN/(cnorm_gaussian_SN(dimensionality)*rfactor_SN*SNR%feat%dr**dimensionality)
      endif
      if (lroot.and.ip==1963) print "(1x,'explode_SN: c_SN =',e11.4)",c_SN
!
!  Mass insertion normalization.
!
      if (lSN_mass) then
        if (mass_profile=="gaussian3") then
          cmass_SN=mass_SN/(cnorm_SN(dimensionality)*width_mass**dimensionality)
        elseif (mass_profile=="gaussian2") then
          cmass_SN=mass_SN/(cnorm_gaussian2_SN(dimensionality)*width_mass**dimensionality)
        elseif (mass_profile=="gaussian") then
          cmass_SN=mass_SN/(cnorm_gaussian_SN(dimensionality)*width_mass**dimensionality)
        endif
      else
        cmass_SN=0.
      endif
      if (lroot.and.ip==1963) print "(1x,'explode_SN: cmass_SN =',e12.5)",cmass_SN
!
!  Velocity insertion normalization.
!  26-aug-10/fred:
!  E_k=int(0.5*rho*vel^2)=approx 2pi*rhom*V0^2*int(r^2*v(r)dr).
!  Total energy=kinetic (kampl_SN) + thermal (eampl_SN).
!
      if (lSN_velocity) then
        if (velocity_profile=="gaussian3") then
          cvelocity_SN=sqrt(2*ktmp/( &
                            SNR%feat%rhom*vnorm_SN(dimensionality)* &
                            width_velocity**dimensionality))
        elseif (velocity_profile=="gaussian2") then
          cvelocity_SN=sqrt(2*ktmp/( &
                            SNR%feat%rhom*vnorm_gaussian2_SN(dimensionality)* &
                            width_velocity**dimensionality))
        elseif (velocity_profile=="gaussian") then
          cvelocity_SN=sqrt(2*ktmp/( &
                            SNR%feat%rhom*vnorm_gaussian_SN(dimensionality)* &
                            width_velocity**dimensionality))
        endif
        if (lroot) print*,'explode_SN: cvelocity_SN is uu_sedov, velocity profile = ', &
                          velocity_profile
      else
        cvelocity_SN=0.
      endif
      if (lroot.and.ip==1963) print"(1x,'explode_SN: cvelocity_SN =',e12.5)",cvelocity_SN
!
!  Validate the explosion.
!
      radius2=energy_Nsigma2*SNR%feat%radius**2
      radius3=2.25*energy_Nsigma2*SNR%feat%radius**2
      radius2mass=SNR%feat%radius**2

      site_mass=0.
      maxlnTT=-10.
      max_cmass=0.
      SNR%feat%EE=0.
      SNR%feat%MM=0.
      SNR%feat%CR=0.
      SNR%indx%state=SNstate_waiting
      SN_TT_ratio_max = SN_TT_ratio*TT_SN_max
      lfound = .false.

      !$omp target if(loffload) map(tofrom: SNR) has_device_addr(f)  ! needs:
      !irho,ilnrho,iss,ilnTT,frac_eth,dr2_SN,rfactor_SN,TT_SN_max, switches FG radius2mass, c_SN, width_energy, width_mass,&
      !irho_lnTT, irho_ss, irho_ee, SN_TT_ratio_max
      !$omp teams distribute parallel do collapse(2) &
      !$omp private(dV,rho_old,rho_new,deltarho,deltauu,deltaEE,deltaCR,ee_old,lnTT,outward_normal_SN,maxTT,rad_hot,deltarho_hot,ind_maxTT,cmass_tmp) &
      !$omp reduction(+:site_mass,SNR%feat%MM,SNR%feat%EE) reduction(max:maxlnTT,max_cmass)
mnloop:do n=n1,n2
      do m=m1,m2
!
        !$ if (lfound) cycle
        call get_dVol(m,n,dV)
!
!  Calculate the distances to the SN origin for all points in the current
!  pencil and store in the dr2_SN global array.
!
        call proximity_SN(SNR)
!
!  Calculate the unperturbed mass and multiply by volume element to derive
!  mass, and sum for the remnant ambient mass, and add ejecta if lSN_mass.
!
        if (ldensity_nolog) then
          rho_old=f(l1:l2,m,n,irho)
        else
          rho_old=exp(f(l1:l2,m,n,ilnrho))
        endif
        rho_new=rho_old
        site_mass=site_mass+sum(rho_old*dV,mask=dr2_SN<=radius2mass)
        if (lSN_mass.and.cmass_SN>0.) then
          call injectmass_SN(deltarho,width_mass,cmass_SN)
          SNR%feat%MM = SNR%feat%MM + sum(deltarho*dV)
          rho_new=rho_old+deltarho
        endif
!
!  Get the unperturbed energy and then add thermal energy if lSN_eth.
!  Check max temperature in dense remnant does not exceed TT_SN_max.
!  Check max temperature within Nsigma of any remant does
!  not exceed TT_SN_max*SN_TT_ratio.
!
        if (lSN_eth) then
          call injectenergy_SN(deltaEE,width_energy,c_SN)
          SNR%feat%EE=SNR%feat%EE+sum(deltaEE*dV)
          if (ltemperature) then
            call eoscalc(irho_lnTT,rho_old,f(l1:l2,m,n,ilnTT),ee=ee_old)
          elseif (lentropy) then
            call eoscalc(irho_ss,rho_old,f(l1:l2,m,n,iss),ee=ee_old)
          endif
          call eoscalc(irho_ee,rho_new,(ee_old*rho_old+deltaEE*frac_eth)/rho_new,lnTT=lnTT)
          ind_maxTT=maxloc( lnTT, mask = (dr2_SN<=radius3), dim=1 )
          if (ind_maxTT>0) then
            maxTT=exp(lnTT(ind_maxTT))
            maxlnTT=max(lnTT(ind_maxTT),maxlnTT)
            if (SNR%feat%radius<=1.01*rfactor_SN*SNR%feat%dr) then
              !dense remnant
              if (maxTT>TT_SN_max) then
                if (present(ierr)) then
                  !$omp atomic
                  ierr=iEXPLOSION_TOO_HOT
                  if (.not.lSN_list) then
                    !$ lfound = .true.
                    exit mnloop
                  endif
                endif
              endif
            else
              !diffuse remnants
              if (maxTT>SN_TT_ratio_max) then
                if (lSN_coolingmass) then
                  rad_hot=dr2_SN(ind_maxTT)
                  deltarho_hot=rho_old(ind_maxTT)*(maxTT/SN_TT_ratio_max - 1.)
                  call getmass_SN(deltarho_hot,rad_hot,width_mass,cmass_tmp)
                  max_cmass=max(max_cmass,cmass_tmp)
                  maxTT=SN_TT_ratio_max
                  maxlnTT=max(log(maxTT),maxlnTT)
                else
                  if (present(ierr)) then
                    !$omp atomic
                    ierr=iEXPLOSION_TOO_HOT
                    if (.not.lSN_list) then
                      !$ lfound = .true.
                      exit mnloop
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif
      enddo
      enddo mnloop
      !$omp end teams distribute parallel do
      !$omp end target
!
!  Broadcast maxlnTT from remnant to all processors so all take the same path
!  after these checks.
!
      if (lSN_eth.or.lSN_coolingmass) then
        mmpi=maxlnTT
        call mpiallreduce_max(mmpi,maxlnTT)
        maxTT=exp(maxlnTT)
      else
        maxTT=cloud_TT
      endif
      if (lSN_eth.and..not.lSN_coolingmass) then
        if (present(ierr)) then
          mpiierr=ierr
          call mpiallreduce_max(mpiierr,ierr)
          if (ierr==iEXPLOSION_TOO_HOT.and..not.lSN_list) return
        endif
      endif
!
!  If cooling mass used to limit max T broadcast cmass_SN and check mass does
!  not exceed Nsol_added solar masses
!
      if (lSN_coolingmass) then
        call mpiallreduce_max(max_cmass,cmass_SN)
        if (lroot.and.ip==1963) print "(1x,'explode_SN: validating cmass_SN =',e12.5)",cmass_SN
        if (cmass_SN>0) then
          SNR%feat%MM=0.
          !$omp target if(loffload) map(tofrom: SNR) has_device_addr(f)  ! gives: m,n FG global width_mass, cmass_SN
          !$omp teams distribute parallel do collapse(2) private(deltarho,dV) reduction(+:SNR%feat%MM)
          do n=n1,n2
          do m=m1,m2
            call get_dVol(m,n,dV)
            call proximity_SN(SNR)
            call injectmass_SN(deltarho,width_mass,cmass_SN)
            SNR%feat%MM=SNR%feat%MM+sum(deltarho*dV)
          enddo
          enddo
          !$omp end teams distribute parallel do
          !$omp end target
          dmpi2_tmp=(/ SNR%feat%MM, SNR%feat%EE, SNR%feat%CR /)
          call mpiallreduce_sum(dmpi2_tmp,dmpi2,3)
          SNR%feat%MM=dmpi2(1)
          SNR%feat%EE=dmpi2(2) !without added kinetic energy
          if (SNR%feat%MM>Nsol_added*solar_mass) then
            if (lroot.and.ip==1963) print '("explode_SN: SNR%feat%MM > ",f5.1," solar mass",f12.4)', &
                                          Nsol_added,SNR%feat%MM/solar_mass
            if (lroot.and.ip==1963) print'("explode_SN: SNR%feat%EE = ",f8.4," SNE")',SNR%feat%EE/ampl_SN
            ierr=iEXPLOSION_TOO_HOT
            if (.not.lSN_list) return
            cmass_SN = cmass_SN*Nsol_added*solar_mass/SNR%feat%MM
          else
            ierr=iEXPLOSION_OK
            if (lroot.and.ip==1963) &
                print "(1x,'explode_SN: SNR%feat%MM = ',f8.4,' solar mass')",SNR%feat%MM/solar_mass
          endif
        endif
        call get_props_check(f,SNR,rhom,ekintot_new,cvelocity_SN,cmass_SN)
        if (ktmp>0) then
          if (ekintot_new-ekintot>ktmp) then
            cvelocity_SN=cvelocity_SN*sqrt(ktmp/(ekintot_new-ekintot))
          else
            if (ekintot_new-ekintot>0) then
              c_SN = min(c_SNmax,c_SN*(ktmp+ekintot-ekintot_new+etmp)/(ktmp+etmp))
            else
              c_SN=c_SNmax
            endif
          endif
        else
          if (ekintot_new-ekintot>0) c_SN = c_SN*etmp/(etmp+ekintot_new-ekintot)
        endif
      endif
!
!  Rescale cvelocity_SN if density distribution yields excessively high
!  kinetic energy to approximate kinetic energy = ktmp.
!
      if (lSN_velocity.and.ktmp>0) then
        call get_props_check(f,SNR,rhom,ekintot_new,cvelocity_SN,cmass_SN)
        if (ekintot_new-ekintot<ktmp) then
          if (lSN_eth) then
            if (ekintot_new-ekintot>0) then
              c_SN = min(c_SNmax,c_SN*(ktmp+ekintot-ekintot_new+etmp)/(ktmp+etmp))
            else
              c_SN = c_SNmax
            endif
          endif
        else
          cvelocity_SN=cvelocity_SN*sqrt(ktmp/(ekintot_new-ekintot))
        endif
      endif
!
!  Remnant parameters pass, so now implement the explosion
!
      if (lroot.and.ip==1963) then
        print "(1x,'explode_SN:         c_SN finally =',e12.5)",c_SN
        print "(1x,'explode_SN:     cmass_SN finally =',e12.5)",cmass_SN
        print "(1x,'explode_SN: cvelocity_SN finally =',e12.5)",cvelocity_SN
      endif

      SNR%feat%EE=0.
      SNR%feat%MM=0.
      SNR%feat%CR=0.

      !$omp target if(loffload) map(tofrom: SNR) has_device_addr(f)  ! needs: irho,ilnrho,iss,ilnTT,iyH,iecr,ifcr, frac_eth, switches;  gives: m,n
      !$omp teams distribute parallel do collapse(2) reduction(+:SNR%feat%MM,SNR%feat%CR,SNR%feat%EE) &
      !$omp private(rho_old,rho_new,deltauu,deltarho,deltaEE,ee_old,yH,lnTT,deltafcr,deltaCR,outward_normal_SN,dV)
!FG: global width_mass, cmass_SN, width_energy, c_SN, width_velocity, cvelocity_SN, ecr_SN (cfcr_SN? to do)
      do n=n1,n2
      do m=m1,m2
!
        call get_dVol(m,n,dV)
!
!  Calculate the distances to the SN origin for all points in the current
!  pencil and store in the dr2_SN global array.
!
        if (lSN_velocity.and.cvelocity_SN>0. .or. lSN_fcr) then
          call proximity_SN(SNR,outward_normal_SN)
        else
          call proximity_SN(SNR)
        endif
!
!  Calculate the unperturbed mass and multiply by volume element to derive
!  mass, and sum for the remnant ambient mass, and add ejecta if lSN_mass.
!  Save changes to f-array.
!
        if (ldensity_nolog) then
          rho_old=f(l1:l2,m,n,irho)
        else
          rho_old=exp(f(l1:l2,m,n,ilnrho))
        endif
        rho_new=rho_old

        if ((lSN_mass.and.cmass_SN>0).or.(lSN_coolingmass.and.cmass_SN>0)) then
          call injectmass_SN(deltarho,width_mass,cmass_SN)
          SNR%feat%MM=SNR%feat%MM+sum(deltarho*dV)
          rho_new=rho_old+deltarho
          if (ldensity_nolog) then
            f(l1:l2,m,n,irho)=rho_new
          else
            f(l1:l2,m,n,ilnrho)=log(rho_new)
          endif
        endif
!
!  Get the unperturbed energy and then add thermal energy if lSN_eth.
!  Save changes to f-array.
!
        if (lSN_eth) then
          call injectenergy_SN(deltaEE,width_energy,c_SN)
          SNR%feat%EE=SNR%feat%EE+sum(deltaEE*dV)
          if (ltemperature) then
            call eoscalc(irho_lnTT,rho_old,f(l1:l2,m,n,ilnTT),ee=ee_old)
            call eoscalc(irho_ee,rho_new,real((ee_old*rho_old+deltaEE*frac_eth)/rho_new),lnTT=f(l1:l2,m,n,ilnTT),yH=yH)
          elseif (lentropy) then
            call eoscalc(irho_ss,rho_old,f(l1:l2,m,n,iss),ee=ee_old)
            call eoscalc(irho_ee,rho_new,real((ee_old*rho_old+deltaEE*frac_eth)/rho_new),ss=f(l1:l2,m,n,iss),yH=yH,lnTT=lnTT)
            if (ilnTT/=0) f(l1:l2,m,n,ilnTT)=lnTT
          endif
          if (iyH/=0) f(l1:l2,m,n,iyH)=yH
        endif
!
!  Apply changes to the velocity field if lSN_velocity.
!  Save changes to f-array.
!
        if (lSN_velocity.and.cvelocity_SN>0.) then
          call injectvelocity_SN(deltauu,width_velocity,cvelocity_SN,SNR,rho_new,outward_normal_SN)
          f(l1:l2,m,n,iux:iuz)=f(l1:l2,m,n,iux:iuz)+deltauu
        endif
!
!  Apply changes to the cosmic ray energy density if lSN_ecr and to
!  cosmicray flux if lSN_ecr.
!  Optionally add cosmicray flux, consistent with addition to ecr, via
!  delta fcr = -K grad (delta ecr) .
!  Still experimental/in testing.
!  Save changes to f-array.
!
        if (lSN_ecr) then
          call injectenergy_SN(deltaCR,width_energy,ecr_SN)
          SNR%feat%CR=SNR%feat%CR+sum(deltaCR*dV)
          f(l1:l2,m,n,iecr) = f(l1:l2,m,n,iecr) + deltaCR
          if (lSN_fcr) then
            call injectfcr_SN(deltafcr,width_energy,ecr_SN,outward_normal_SN)
! cfcr_SN needs to be computed independently as cvelocity_SN and scaling included
            f(l1:l2,m,n,ifcr:ifcr+2)=f(l1:l2,m,n,ifcr:ifcr+2) + deltafcr
          endif
        endif
      enddo
      enddo  !  mn-loop
      !$omp end teams distribute parallel do
      !$omp end target
!
      call get_properties(f,SNR,rhom,ekintot_new,rhomin,ierr)
      if (lroot.and.ip==1963) then
        print"(1x,'explode_SN: total new kinetic energy =',e12.5)", ekintot_new
        print"(1x,'explode_SN: total kinetic energy change =',e12.5)",ekintot_new-ekintot
      endif
!
!  Sum and share diagnostics etc. amongst processors.
!
      dmpi2_tmp=(/ SNR%feat%MM, SNR%feat%EE, SNR%feat%CR /)
      call mpiallreduce_sum(dmpi2_tmp,dmpi2,3)
      SNR%feat%MM=dmpi2(1)
      SNR%feat%EE=dmpi2(2)+ekintot_new-ekintot !include added kinetic energy
      SNR%feat%CR=dmpi2(3)
      call mpiallreduce_sum(site_mass,mmpi)
      site_mass=mmpi
!
! FAG need to consider effect of CR and fcr on total energy for data collection
! and the energy budget applied to the SNR similar to kinetic energy?
!
      if (lroot.and.ip==1963) then
        print "(1x,'explode_SN: SNR%feat%MM = ',f8.4)",SNR%feat%MM/solar_mass
        print "(1x,'explode_SN: SNR%feat%EE = ',f8.4)",SNR%feat%EE/ampl_SN
        print "(1x,'explode_SN: SNR%feat%EK = ',f8.4)",(ekintot_new-ekintot)/ampl_SN
        print "(1x,'explode_SN: SNR%feat%CR = ',f8.4)",SNR%feat%CR/ampl_SN
      endif
      if (.not. lSN_list) then
        if (SNR%indx%SN_type==1) then
          call set_next_SNI(t_interval_SN)
          SNrate = t_interval(SNI)
        else
          call set_next_SNII(f,t_interval_SN)
          SNrate = t_interval(SNII)
        endif
      else
        if (SNR%indx%SN_type==1) then
          t_interval_SN=t_interval(SNI)
          SNrate = t_interval(SNI)
        else
          t_interval_SN=t_interval(SNII)
          SNrate = t_interval(SNII)
        endif
      endif
!
      if (lroot) then
!
        if (lOB_cluster) then
          print "(1x,'explode_SN: t_cluster',e12.5)", t_cluster
          print "(1x,'explode_SN: x_cluster', f8.4)", x_cluster
          print "(1x,'explode_SN: z_cluster', f8.4)", z_cluster
        endif
!
        if (lfirst_warning) then
          call warning('sn_series.dat','new column SN_rate added 19.10.17 '// &
          'continuation of old data may need header and extra column appended')
          call warning('sn_series.dat','new columns rhom, Ekin, Ecr 20.03.19 '// &
          'continuation of old data may need header and extra columns appended')
          call warning('sn_series.dat','new columns added_Nsol 27.05.21 '// &
          'continuation of old data may need header and extra columns appended')
        endif

        open(1,file=trim(datadir)//'/sn_series.dat',position='append')
        print "(1x,'explode_SN:    step, time = ',i8,e12.5)",it,t
        print "(1x,'explode_SN:          dVol = ',   e12.5)",dVol_glob
        print "(1x,'explode_SN:       SN type = ',      i3)",SNR%indx%SN_type
        print "(1x,'explode_SN: proc, l, m, n = ',     4i6)",SNR%indx%iproc,SNR%indx%l,SNR%indx%m,SNR%indx%n
        print "(1x,'explode_SN:       x, y, z = ',   3f8.4)",SNR%feat%x,SNR%feat%y,SNR%feat%z
        print "(1x,'explode_SN:remnant radius = ',    f8.4)",SNR%feat%radius
        print "(1x,'explode_SN:       rho, TT = ',  2e11.4)",SNR%site%rho,SNR%site%TT
        print "(1x,'explode_SN:    maximum TT = ',   e12.5)",maxTT
        print "(1x,'explode_SN:  Mean density = ',   e12.5)",SNR%feat%rhom
        print "(1x,'explode_SN:  Total energy = ',   e12.5)",SNR%feat%EE+SNR%feat%CR
        print "(1x,'explode_SN:    CR energy  = ',   e12.5)",SNR%feat%CR
        print "(1x,'explode_SN:    Added mass = ',   e12.5)",SNR%feat%MM/solar_mass
        print "(1x,'explode_SN:  Ambient Nsol = ',   e12.5)",site_mass/solar_mass
        print "(1x,'explode_SN:    Sedov time = ',   e12.5)", SNR%feat%t_sedov
        print "(1x,'explode_SN:   Shell speed = ',   e12.5)",uu_sedov
        write(1,'(i10,E16.8,5i6,16E13.5)') &
            it, t, SNR%indx%SN_type, SNR%indx%iproc, SNR%indx%l, SNR%indx%m, &
            SNR%indx%n, SNR%feat%x, SNR%feat%y, SNR%feat%z, SNR%site%rho, &
            SNR%feat%rhom, SNR%site%TT, SNR%feat%EE+SNR%feat%CR, &
            ekintot_new-ekintot, SNR%feat%CR, SNR%feat%t_sedov, &
            SNR%feat%radius, site_mass/solar_mass, SNR%feat%MM/solar_mass, &
            maxTT, t_interval_SN, SNrate
        close(1)
      endif
      lfirst_warning=.false.
!
      if (present(preSN)) then
        do i=2,npreSN
          preSN(:,i-1)= preSN(:,i)
        enddo
        preSN(1,npreSN)= SNR%indx%l
        preSN(2,npreSN)= SNR%indx%m
        preSN(3,npreSN)= SNR%indx%n
        preSN(4,npreSN)= SNR%indx%iproc
      endif
      SNR%indx%state=SNstate_finished
!
      if (present(ierr).and.lSN_list) then
        select case (ierr)
          case (iEXPLOSION_TOO_HOT)
            if (lroot.and.ip==1963) print &
                "(1x,'explode_SN: TOO HOT, (x,y,z) =',3f8.4,', rho, rad =',2e11.4)", &
                SNR%feat%x, SNR%feat%y, SNR%feat%z,SNR%site%rho, SNR%feat%radius
          case (iEXPLOSION_TOO_UNEVEN)
            if (lroot.and.ip==1963) print &
                "(1x,'explode_SN: TOO UNEVEN, (x,y,z) =',3f8.4,', rho, rad =',2e11.4)", &
                SNR%feat%x, SNR%feat%y, SNR%feat%z,SNR%site%rho, SNR%feat%radius
          case (iEXPLOSION_TOO_RARIFIED)
            if (lroot.and.ip==1963) print &
                "(1x,'explode_SN: TOO RARIFIED, (x,y,z) =',3f8.4,', rho, rad =',2e11.4)", &
                SNR%feat%x, SNR%feat%y, SNR%feat%z,SNR%site%rho, SNR%feat%radius
        endselect
      endif
!
      if (present(ierr)) ierr=iEXPLOSION_OK
!
    endsubroutine explode_SN
!***********************************************************************
    subroutine get_properties(f,remnant,rhom,ekintot,rhomin,ierr)   !needs dr2_SN
!
!  Calculate integral of mass cavity profile and total kinetic energy.
!
!  22-may-03/tony: coded
!  10-oct-09/axel: return zero density if the volume is zero
!
      use Sub, only: dot2
      use Mpicomm, only: mpiallreduce_sum,mpiallreduce_min,mpiallreduce_max
      use Grid, only: get_dVol
!
      real, intent(in), dimension(mx,my,mz,mfarray) :: f
      type (SNRemnant), intent(in) :: remnant
      real, intent(out) :: rhom, ekintot, rhomin
      integer, optional :: ierr

      real :: rhotmp, rhomax, radius2
      real, dimension(nx) :: rho, u2, dV
      real, dimension(nx,3) :: uu
      logical, dimension(nx) :: lmask
      real, dimension(2) :: tmp,tmp2
!
!  inner rad defined to determine mean density inside rad and smooth if desired
!
      radius2 = remnant%feat%radius**2
      tmp=0.0
      rhomin=1e20
      rhomax=0.0
!
!  Obtain distance to SN and sum all points inside SNR radius and
!  divide by number of points.
!
      !$omp target if(loffload) map(from: rhomin,rhomax,tmp) map(to: remnant) has_device_addr(f)  ! needs: irho,ilnrho,iuu,dr2_SN, switches;  gives: m,n
! FG: radius2
      !$omp teams distribute parallel do collapse(2) private(rho,uu,u2,lmask,dV), reduction(min:rhomin) reduction(max:rhomax) reduction(+:tmp)
      do n=n1,n2
      do m=m1,m2
        call get_dVol(m,n,dV)
        call proximity_SN(remnant)
!
!  get rho from existing ambient density everywhere
!
        if (ldensity_nolog) then
          rho=f(l1:l2,m,n,irho)
        else
          rho=exp(f(l1:l2,m,n,ilnrho))
        endif
!
!  if radius scaled to total mass mask mass outside inner rad
!
        if (lSN_scale_rad) then
          lmask = (dr2_SN <= radius2)
          if (any(lmask)) then
            rhomin=min(rhomin,minval(rho,mask=lmask))
            rhomax=max(rhomax,maxval(rho,mask=lmask))
          endif
        endif
        uu=f(l1:l2,m,n,iuu:iuu+2)
!
!  avoid NaN where uu less than tini
!
        where (abs(f(l1:l2,m,n,iuu:iuu+2))<sqrt(tini)) uu=0.0
        call dot2(uu,u2)
!
!  compute kinetic energy everywhere before applying the mask
!
        tmp(2)=tmp(2)+sum(rho*u2*dV)
        tmp(1)=tmp(1)+sum(rho*dV,mask=(dr2_SN <= radius2))
      enddo
      enddo  !  mn-loop
      !$omp end teams distribute parallel do
      !$omp end target
!
!  Calculate mean density inside the remnant and return error if the volume is
!  zero.
!
      call mpiallreduce_sum(tmp,tmp2,2)
      ekintot=0.5*tmp2(2)
      if ((lSN_velocity).and.(abs(tmp2(1)) < tini)) then
        if (lroot) print*,'enclosed mass, total kinetic energy = ', tmp2(1), ekintot
        call fatal_error("interstellar.get_properties","dividing by zero?")
      else
        rhom=tmp2(1)*0.75*pi_1/radius2**1.5
      endif
!
!  Determine the density rarification ratio in order to avoid excessive spikes
!
      rhotmp=rhomin
      call mpiallreduce_min(rhotmp,rhomin)
      rhotmp=rhomax
      call mpiallreduce_max(rhotmp,rhomax)
      if (present(ierr)) then
        if (rhomax/rhomin > SN_rho_ratio) ierr=iEXPLOSION_TOO_UNEVEN
        if (lroot.and.ip==1963) then
          print "(1x,'get_properties: rhomax =',e11.4)",rhomax
          print "(1x,'get_properties: rhomin =',e11.4)",rhomin
          print "(1x,'get_properties: ierr   =',   i3)",ierr
          print "(1x,'get_properties: radius =', f8.4)",remnant%feat%radius
          print "(1x,'get_properties: ekintot =', f8.4)",ekintot
        endif
      endif
!
    endsubroutine get_properties
!*****************************************************************************
    subroutine get_props_check(f,remnant,rhom,ekintot,cvelocity_SN,cmass_SN)   !needs dr2_SN
!
!  Calculate integral of mass cavity profile and total kinetic energy.
!
!  22-may-03/tony: coded
!  10-oct-09/axel: return zero density if the volume is zero
!
      use Sub, only: dot2
      use Mpicomm, only: mpiallreduce_sum
      use Grid, only: get_dVol
!
      real, intent(in), dimension(mx,my,mz,mfarray) :: f
      type (SNRemnant), intent(in) :: remnant
      real, intent(in) :: cvelocity_SN, cmass_SN
      real, intent(out) :: ekintot, rhom

      real :: radius2, width_mass, width_velocity
      real, dimension(nx) :: rho, u2, deltarho, dV
      real, dimension(nx,3) :: uu,deltauu, outward_normal_SN
      real, dimension(2) :: tmp,tmp2
!
!  inner rad defined to determine mean density inside rad and smooth if desired
!
      width_mass     = remnant%feat%radius*mass_width_ratio
      width_velocity = remnant%feat%radius*velocity_width_ratio
      radius2        = energy_Nsigma2*remnant%feat%radius**2
!
!  Obtain distance to SN and sum all points inside SNR radius and divide by number of points.
!
!NEW CHANGES
      tmp=0.
      !$omp target if(loffload) map(from: tmp) map(to: remnant) has_device_addr(f)  ! needs: irho,ilnrho,iuu,dr2_SN, switches;  gives: m,n
!FG: cvelocity_SN, cmass_SN, width_mass, width_velocity, radius2
      !$omp teams distribute parallel do collapse(2) private(uu,u2,rho,deltauu,deltarho,dV,outward_normal_SN) reduction(+:tmp)
      do n=n1,n2
      do m=m1,m2
        call get_dVol(m,n,dV)
        if (lSN_velocity.and.cvelocity_SN>0) then
          call proximity_SN(remnant,outward_normal_SN)
        else
          call proximity_SN(remnant)
        endif
!
!  get rho from existing ambient density everywhere
!
        if (ldensity_nolog) then
          rho=f(l1:l2,m,n,irho)
        else
          rho=exp(f(l1:l2,m,n,ilnrho))
        endif
        if ((lSN_mass.or.lSN_coolingmass).and.cmass_SN>0) then
          call injectmass_SN(deltarho,width_mass,cmass_SN)
          rho=rho+deltarho
        endif
        uu=f(l1:l2,m,n,iuu:iuu+2)
        if (lSN_velocity.and.cvelocity_SN>0) then
          call injectvelocity_SN(deltauu,width_velocity,cvelocity_SN,remnant,rho,outward_normal_SN)
          uu=uu+deltauu
        endif
!
!  avoid NaN where uu less than tini  MR: needed?
!  FG: experienced underflow errors for velocities close to zero for dot2,
!      particularly for sedov test with momentum injection
!
        where (abs(f(l1:l2,m,n,iuu:iuu+2))<sqrt(tini)) uu=0.
        call dot2(uu,u2)
!
!  compute kinetic energy everywhere before applying the mask
!
        tmp(2)=tmp(2)+sum(rho*u2*dV)
        tmp(1)=tmp(1)+sum(rho*dV,mask=(dr2_SN <= radius2))
      enddo
      enddo  !  mn-loop
      !$omp end teams distribute parallel do
      !$omp end target
!
!  Calculate mean density inside the remnant and return zero if the volume is
!  zero.
!
      call mpiallreduce_sum(tmp,tmp2,2)
      ekintot=0.5*tmp2(2)
      if ((lSN_velocity).and.(abs(tmp2(1)) < tini)) then
        call fatal_error("interstellar.get_props_check","dividing by zero?")
      else
        rhom=tmp2(1)*0.75*pi_1/radius2**1.5
      endif
      if (lroot.and.ip==1963) print"(1x,'get_props_check: rhom =',e11.4,' ekintot =',f8.4)",rhom,ekintot
!
    endsubroutine get_props_check
!*****************************************************************************
    subroutine get_lowest_rho(f,SNR,radius,rho_lowest)   !needs dr2_SN
!
!  Calculate integral of mass cavity profile.
!
!  22-may-03/tony: coded
!
      use Mpicomm, only: mpiallreduce_max
!
      real, intent(in), dimension(mx,my,mz,mfarray) :: f
      type (SNRemnant), intent(inout) :: SNR
      real, intent(in) :: radius
      real, intent(out) :: rho_lowest
      real :: tmp
      real :: radius2
      real, dimension(nx) :: rho
!
!  Find lowest rho value in the surronding cavity.
!
      rho_lowest=1E10
      radius2 = energy_Nsigma2*radius**2
      do n=n1,n2
      do m=m1,m2
        call proximity_SN(SNR)
        if (ldensity_nolog) then
          rho=f(l1:l2,m,n,irho)
        else
          rho=exp(f(l1:l2,m,n,ilnrho))
        endif
        rho_lowest=min(rho_lowest,minval(rho,mask=(dr2_SN <= radius2)))
      enddo
      enddo
!
      tmp=rho_lowest
      call mpiallreduce_max(tmp,rho_lowest)
!
    endsubroutine get_lowest_rho
!*****************************************************************************
    subroutine proximity_OB(OB,outward_normal)    ! make pars: m,n, dr2_SN
!$omp declare target    ! needs: x,y,z, L[xyz, lperi; gives: dr2_OB
!
!  Calculate pencil of distance to OB cluster origin.
!
!  06-aug-19/fred: cloned from proximity_SN
!
      type (SNCluster),      intent(inout) :: OB
      real, dimension(nx,3), intent(out), optional :: outward_normal
!
      real,dimension(nx) :: dx_OB, dr_OB
      real :: dy_OB
      real :: dz_OB
!
!  Obtain distance to OB cluster origin
!
      dx_OB=x(l1:l2)-OB%feat%x
      if (lperi(1)) then
        where (dx_OB >  Lx/2) dx_OB=dx_OB-Lx
        where (dx_OB < -Lx/2) dx_OB=dx_OB+Lx
      endif
!
      dy_OB=y(m)-OB%feat%y
      if (lperi(2)) then
        if (dy_OB >  Ly/2) dy_OB=dy_OB-Ly
        if (dy_OB < -Ly/2) dy_OB=dy_OB+Ly
      endif
!
      dz_OB=z(n)-OB%feat%z
      if (lperi(3)) then
        if (dz_OB >  Lz/2) dz_OB=dz_OB-Lz
        if (dz_OB < -Lz/2) dz_OB=dz_OB+Lz
      endif
!
      dr2_OB=dx_OB**2 + dy_OB**2 + dz_OB**2
!
      if (present(outward_normal)) then
        dr_OB=sqrt(dr2_OB)
        dr_OB=max(dr_OB,tiny(0.0))
!
!  Avoid dr_SN = 0 above to avoid div by zero below.
!
        outward_normal(:,1)=dx_OB/dr_OB
        where (dr2_OB < tini) outward_normal(:,1)=0.0
        outward_normal(:,2)=dy_OB/dr_OB
        where (dr2_OB < tini) outward_normal(:,2)=0.0
        outward_normal(:,3)=dz_OB/dr_OB
        where (dr2_OB < tini) outward_normal(:,3)=0.0
      endif
!
    endsubroutine proximity_OB
!*****************************************************************************
    subroutine proximity_SN(SNR,outward_normal)    ! make pars: m,n, dr2_SN
!$omp declare target    ! needs: x,y,z, L[xyz], deltay
!
!  Calculate pencil of distance to SN explosion site.
!
!  20-may-03/tony: extracted from explode_SN code written by grs
!  22-may-03/tony: pencil formulation
!
      type (SNRemnant),                intent(in) :: SNR
      real, dimension(nx,3), optional, intent(out) :: outward_normal
!
      real,dimension(nx) :: dx_SN, dr_SN
      real :: dy_SN
      real :: dz_SN
      integer :: l
!
!  Obtain distance to SN
!
      dx_SN=x(l1:l2)-SNR%feat%x
      if (lperi(1)) then
        where (dx_SN >  Lx/2) dx_SN=dx_SN-Lx
        where (dx_SN < -Lx/2) dx_SN=dx_SN+Lx
      endif
!
      dy_SN=y(m)-SNR%feat%y
      if (lperi(2)) then
        if (dy_SN >  Ly/2) dy_SN=dy_SN-Ly
        if (dy_SN < -Ly/2) dy_SN=dy_SN+Ly
      endif
!
      dz_SN=z(n)-SNR%feat%z
      if (lperi(3)) then
        if (dz_SN >  Lz/2) dz_SN=dz_SN-Lz
        if (dz_SN < -Lz/2) dz_SN=dz_SN+Lz
      endif
!
      if (lshear) then
        do l=l1,l2
          if (x(l)-SNR%feat%x> Lx/2) then
            if (y(m)-SNR%feat%y >  Ly/2) then
               dr2_SN(l-nghost)=dx_SN(l-nghost)**2 + (y(m)-SNR%feat%y+deltay-Ly)**2 + dz_SN**2
            elseif (y(m)-SNR%feat%y < -Ly/2) then
               dr2_SN(l-nghost)=dx_SN(l-nghost)**2 + (y(m)-SNR%feat%y+deltay+Ly)**2 + dz_SN**2
            else
               dr2_SN(l-nghost)=dx_SN(l-nghost)**2 + (dy_SN+deltay)**2 + dz_SN**2
            endif
          elseif (x(l)-SNR%feat%x< -Lx/2) then
            if (y(m)-SNR%feat%y >  Ly/2) then
               dr2_SN(l-nghost)=dx_SN(l-nghost)**2 + (y(m)-SNR%feat%y-deltay-Ly)**2 + dz_SN**2
            elseif (y(m)-SNR%feat%y < -Ly/2) then
               dr2_SN(l-nghost)=dx_SN(l-nghost)**2 + (y(m)-SNR%feat%y-deltay+Ly)**2 + dz_SN**2
            else
               dr2_SN(l-nghost)=dx_SN(l-nghost)**2 + (dy_SN-deltay)**2 + dz_SN**2
            endif
          else
            dr2_SN(l-nghost)=dx_SN(l-nghost)**2 + dy_SN**2 + dz_SN**2
          endif
        enddo
      else
        dr2_SN=dx_SN**2 + dy_SN**2 + dz_SN**2
      endif
!
      if (present(outward_normal)) then
!
        dr_SN=sqrt(dr2_SN)
        dr_SN=max(dr_SN,tiny(0.))
!
!  Avoid dr_SN = 0 above to avoid div by zero below.
!
        outward_normal(:,1)=dx_SN/dr_SN
        where (dr2_SN < tini) outward_normal(:,1)=0.0
        if (lshear) then
          do l=l1,l2
            if (x(l)-SNR%feat%x> Lx/2) then
              if (y(m)-SNR%feat%y >  Ly/2) then
                outward_normal(l-nghost,2)=(y(m)-SNR%feat%y+deltay-Ly)/ &
                    sqrt(dx_SN(l-nghost)**2 + (y(m)-SNR%feat%y+deltay-Ly)**2 + dz_SN**2)
              elseif (y(m)-SNR%feat%y < -Ly/2) then
                outward_normal(l-nghost,2)=(y(m)-SNR%feat%y+deltay+Ly)/ &
                    sqrt(dx_SN(l-nghost)**2 + (y(m)-SNR%feat%y+deltay+Ly)**2 + dz_SN**2)
              else
                if (dr2_SN(l-nghost)==0) then
                  outward_normal(l-nghost,2)=0.0
                else
                  outward_normal(l-nghost,2)=dy_SN/dr_SN(l-nghost)
                endif
              endif
            elseif (x(l)-SNR%feat%x< -Lx/2) then
              if (y(m)-SNR%feat%y >  Ly/2) then
                outward_normal(l-nghost,2)=(y(m)-SNR%feat%y-deltay-Ly)/ &
                    sqrt(dx_SN(l-nghost)**2 + (y(m)-SNR%feat%y-deltay-Ly)**2 + dz_SN**2)
              elseif (y(m)-SNR%feat%y < -Ly/2) then
                outward_normal(l-nghost,2)=(y(m)-SNR%feat%y-deltay+Ly)/ &
                    sqrt(dx_SN(l-nghost)**2 + (y(m)-SNR%feat%y-deltay+Ly)**2 + dz_SN**2)
              else
                if (dr2_SN(l-nghost)==0) then
                  outward_normal(l-nghost,2)=0.0
                else
                  outward_normal(l-nghost,2)=dy_SN/dr_SN(l-nghost)
                endif
              endif
            else
              dr2_SN(l-nghost)=dx_SN(l-nghost)**2 + dy_SN**2 + dz_SN**2
            endif
          enddo
        else
          outward_normal(:,2)=dy_SN/dr_SN
          where (dr2_SN < tini) outward_normal(:,2)=0.0
        endif
        outward_normal(:,3)=dz_SN/dr_SN
        where (dr2_SN <tini) outward_normal(:,3)=0.0
      endif
!
    endsubroutine proximity_SN
!*****************************************************************************
    subroutine injectenergy_SN(deltaEE,width,c_SN)
!$omp declare target    ! needs: dr2_SN, sigma_SN1
!
      real, intent(in) :: width,c_SN
      real, intent(out), dimension(nx) :: deltaEE
!
      real, dimension(nx) :: profile_SN
!
!  Whether mass is moved or not, inject energy.
!
      if (thermal_profile=="gaussian3") then
        profile_SN=exp(-(dr2_SN/width**2)**3)
      elseif (thermal_profile=="gaussian2") then
        profile_SN=exp(-(dr2_SN/width**2)**2)
      elseif (thermal_profile=="gaussian") then
        profile_SN=exp(-(dr2_SN/width**2))
      elseif (thermal_profile=="quadratic") then
        profile_SN=max(1.-(dr2_SN/width**2),0.)
      elseif (thermal_profile=="quadratictanh") then
        profile_SN=max(1.-(dr2_SN/width**2),0.)*0.5*(1.-tanh((sqrt(dr2_SN)-width)*sigma_SN1))
      elseif (thermal_profile=="quartictanh") then
        profile_SN=max(1.-(dr2_SN/width**2)**2,0.)*0.5*(1.-tanh((sqrt(dr2_SN)-width)*sigma_SN1))
      elseif (thermal_profile=="tanh") then
        profile_SN=(1.-tanh((sqrt(dr2_SN)-width)*sigma_SN1))*0.5
      endif
!
      deltaEE=c_SN*profile_SN ! spatial energy density
!
    endsubroutine injectenergy_SN
!*****************************************************************************
    subroutine injectmass_SN(deltarho,width,cmass_SN)
!$omp declare target     ! needs: dr2_SN, sigma_SN1
!
      real, intent(in) :: width,cmass_SN
      real, intent(out), dimension(nx) :: deltarho
      real, dimension(nx) :: profile_SN
!
!  Inject mass.
!
      if (mass_profile=="gaussian3") then
        profile_SN=exp(-(dr2_SN/width**2)**3)
      elseif (mass_profile=="gaussian2") then
        profile_SN=exp(-(dr2_SN/width**2)**2)
      elseif (mass_profile=="gaussian") then
        profile_SN=exp(-(dr2_SN/width**2))
      elseif (mass_profile=="quadratic") then
        profile_SN=max(1.0-(dr2_SN/width**2),0.0)
      elseif (mass_profile=="tanh") then
!
!  This is normally handled in the mass movement section
!
        profile_SN=(1.-tanh((sqrt(dr2_SN)-width)*sigma_SN1))*0.5
      endif
!
      deltarho=cmass_SN*profile_SN ! spatial mass density
!
    endsubroutine injectmass_SN
!*****************************************************************************
    subroutine getmass_SN(deltarho,rad,width,cmass_SN)
!$omp declare target     ! needs: sigma_SN1; switches
!
      real, intent(in) :: width, rad, deltarho
      real, intent(out) :: cmass_SN
!
      real :: profile_SN
!
!  Determine prefactor to increase density locally by deltarho.
!
      if (mass_profile=="gaussian3") then
        profile_SN=exp(-(rad**2/width**2)**3)
      elseif (mass_profile=="gaussian2") then
        profile_SN=exp(-(rad**2/width**2)**2)
      elseif (mass_profile=="gaussian") then
        profile_SN=exp(-(rad**2/width**2))
      elseif (mass_profile=="quadratic") then
        profile_SN=max(1.0-(rad**2/width**2),0.0)
      elseif (mass_profile=="tanh") then
!
!  This is normally handled in the mass movement section
!
        profile_SN=(1.-tanh((rad-width)*sigma_SN1))*0.5
      endif
!
      cmass_SN=deltarho/profile_SN !mass prefactor
!
    endsubroutine getmass_SN
!***********************************************************************
    subroutine injectvelocity_SN(deltauu,width,cvelocity_SN,SNR,rho,outward_normal)
!$omp declare target    ! needs: dr2_SN; switches
!
      real, intent(in) :: width,cvelocity_SN
      real, intent(in), dimension(nx,3) :: outward_normal
      real, intent(out), dimension(nx,3) :: deltauu
      type (SNRemnant), intent(in) :: SNR
      real, dimension(nx), intent(in) :: rho
!
      real, dimension(nx) :: profile_SN, ekin_scale
!
      integer :: j
!
!  Calculate deltauu.
!
      if (velocity_profile=="gaussian") then
        profile_SN=exp(-(dr2_SN/width**2))
!
      elseif (velocity_profile=="gaussian2") then
        profile_SN=exp(-(dr2_SN/width**2)**2)
!
      elseif (velocity_profile=="gaussian3") then
        profile_SN=exp(-(dr2_SN/width**2)**3)
!
      endif
!
!  lSN_momentum: kinetic energy is applied with a velocity_profile
!                given the turbulent ambient velocity depends on
!                sqrt(rho), else outward normal velocity profile
!                independent of density distribution
!
      if (lSN_momentum) then
        ekin_scale=cvelocity_SN*sqrt(SNR%feat%rhom/rho)
      else
        ekin_scale=cvelocity_SN
      endif
      do j=1,3
        deltauu(:,j)=ekin_scale*profile_SN*outward_normal(:,j) ! spatial outflow
      enddo
!
    endsubroutine injectvelocity_SN
!***********************************************************************
    subroutine injectfcr_SN(deltafcr,width,cfcr_SN,outward_normal)
!$omp declare target   ! needs: dr2_SN,kperp, switches
!
      real, intent(in) :: width,cfcr_SN
      real, intent(in), dimension(nx,3) :: outward_normal
      real, intent(out), dimension(nx,3) :: deltafcr
!
      real, dimension(nx) :: profile_SN
!
      integer :: j
!
!  Calculate cosmicray flux, delta fcr, consistent with addition to ecr,
!  via  delta fcr = -K grad (delta ecr).  (Currently using K=kperp.)
!
      if (thermal_profile=="gaussian") then
        profile_SN=2.*sqrt(dr2_SN)/width**2*exp(-(dr2_SN/width**2))
!
      elseif (thermal_profile=="gaussian2") then
        profile_SN=4.*(sqrt(dr2_SN)**3)/width**4*exp(-(dr2_SN/width**2)**2)
!
      elseif (thermal_profile=="gaussian3") then
        profile_SN=6.*(sqrt(dr2_SN)**5)/width**6*exp(-(dr2_SN/width**2)**3)
!
      endif
!
      do j=1,3
        deltafcr(:,j)=cfcr_SN*profile_SN*kperp * outward_normal(:,j) ! spatial CR flux
      enddo
!
    endsubroutine injectfcr_SN
!*****************************************************************************
    function get_free_SNR()
!
      integer :: get_free_SNR
      integer :: i,iSNR
!
      if (nSNR>=mSNR) call fatal_error("get_free_SNR","run out of SNR slots. Increase mSNR")
!
      iSNR=-1
      do i=1,mSNR
        if (SNRs(i)%indx%state==SNstate_invalid) then
          iSNR=i
          exit
        endif
      enddo
!
      if (iSNR<0) call fatal_error("get_free_SNR", &
          "could not find an empty SNR slot. Slots were not properly freed")
!
      nSNR=nSNR+1
      SNRs(iSNR)%indx%state=SNstate_waiting
      SNR_index(nSNR)=iSNR
      get_free_SNR=iSNR
!
    endfunction get_free_SNR
!*****************************************************************************
    subroutine free_SNR(iSNR)
!
      integer :: i,iSNR
!
      if (SNRs(iSNR)%indx%state==SNstate_invalid) then
        if (lroot) print*,"Tried to free an already invalid SNR"
        return
      endif
!
      nSNR=nSNR-1
      SNRs(iSNR)%indx%state=SNstate_invalid
!
      do i=iSNR,nSNR
        SNR_index(i)=SNR_index(i+1)
      enddo
!
    endsubroutine free_SNR
!*****************************************************************************
    subroutine tidy_SNRs
!
      integer :: i
!
      do i=1,mSNR
        if (SNRs(i)%indx%state==SNstate_finished) call free_SNR(i)
      enddo
!
    endsubroutine tidy_SNRs
!*****************************************************************************
    function get_free_OB()
!
      integer :: get_free_OB
      integer :: i,iOB
!
      if (nOB>=mOB) call fatal_error("get_free_OB","run out of OB slots. Increase mOB")
!
      iOB=-1
      do i=1,mOB
        if (OBs(i)%indx%state==OBstate_invalid) then
          iOB=i
          exit
        endif
      enddo
!
      if (iOB<0) call fatal_error("get_free_OB", &
          "Could not find an empty OB slot. Slots were not properly freed")
!
      nOB=nOB+1
      OBs(iOB)%indx%state=OBstate_waiting
      OB_index(nOB)=iOB
      get_free_OB=iOB
!
    endfunction get_free_OB
!*****************************************************************************
    subroutine free_OB(iOB)
!
      integer :: i,iOB
!
      if (OBs(iOB)%indx%state==OBstate_invalid) then
        if (lroot) print*,"Tried to free an already invalid OB"
        return
      endif
!
      nOB=nOB-1
      OBs(iOB)%indx%state=OBstate_invalid
!
      do i=iOB,nOB
        OB_index(i)=OB_index(i+1)
      enddo
!
    endsubroutine free_OB
!*****************************************************************************
    subroutine tidy_OBs
!
      integer :: i
!
      do i=1,mOB
        if (OBs(i)%indx%state==OBstate_finished) call free_OB(i)
      enddo
!
    endsubroutine tidy_OBs
!***********************************************************************
    subroutine set_next_OB(t_interval_OB)
!
      use General, only:  random_seed_wrapper, random_number_wrapper
!
      real, intent(in) :: t_interval_OB
!
      real :: franSN
!
!  Pre-determine time for next OB.
!
      if (lroot) print"(1x,'set_next_OB: Old t_next_OB =',e11.4)",t_cluster
      call random_number_wrapper(franSN)
!
!  Time interval follows Poisson process with rate 1/t_interval_OB
!
      t_cluster=t-log(franSN)*t_interval_OB
!
      if (lroot) print"(1x,'set_next_OB: Next OB at time =',e11.4)",t_cluster
!
    endsubroutine set_next_OB
!*****************************************************************************
    subroutine pushpars2c(p_par)

      use Syscalls, only: copy_addr, copy_addr_dble_1D

      integer, parameter :: n_pars=15
      integer(KIND=ikind8), dimension(n_pars) :: p_par
!
      call copy_addr(GammaUV,p_par(1))
      call copy_addr(cUV,p_par(2))
      call copy_addr(T0UV,p_par(3))
      call copy_addr(ncool,p_par(4))                  ! int
!
      call copy_addr(laverage_SNI_heating,p_par(5))   ! int
      call copy_addr(laverage_SNII_heating,p_par(6))  ! int
      call copy_addr(heatingfunction_scale,p_par(7))  ! (2)
      call copy_addr(average_SNI_heating,p_par(8))
      call copy_addr(average_SNII_heating,p_par(9))
      call copy_addr(h_SNI,p_par(10))
      call copy_addr(h_SNII,p_par(11))

      call copy_addr(lncoolT,p_par(12))          ! (ncool+1)
      call copy_addr_dble_1D(lncoolH,p_par(13))  ! (ncool)
      call copy_addr(coolB,p_par(14))            ! (ncool))
      call copy_addr(heating_rate_code,p_par(15))
!      call copy_addr(heat_z,p_par(15))   ! (mz)

    endsubroutine pushpars2c
!*******************************************************************
endmodule Interstellar
