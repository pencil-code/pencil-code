! $Id$
!
!  MODULE_DOC: This modules adds chemical species and reactions.
!  MODULE_DOC: The units used in the chem.in files are cm3,mole,sec,kcal and K
!  This was found out by comparing the mechanism found in
!  samples/0D/chemistry_H2_ignition_delay
!  with Flow Reactor Studies and Kinetic Modeling of the ReactionH/O22
!  of  A. MUELLER, T. J. KIM, R. A. YETTER, F. L. DRYER
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this
! CPARAM logical, parameter :: lchemistry = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 1
!
! PENCILS PROVIDED cv; cv1; cp; cp1; glncp(3);  gXXk(3,nchemspec); gamma
! PENCILS PROVIDED nu; gradnu(3); gYYk(3,nchemspec); chem_conc(nchemspec)
! PENCILS PROVIDED DYDt_reac(nchemspec); DYDt_diff(nchemspec)
! PENCILS PROVIDED lambda; glambda(3); lambda1; gdiffk(3,nchemspec)
! PENCILS PROVIDED g2XXk(nchemspec)
! PENCILS PROVIDED Diff_penc_add(nchemspec); H0_RT(nchemspec); hhk_full(nchemspec)
! PENCILS PROVIDED ghhk(3,nchemspec); S0_R(nchemspec); cs2
! PENCILS PROVIDED glnpp(3); del2pp; mukmu1(nchemspec)
! PENCILS PROVIDED ccondens; ppwater
! PENCILS PROVIDED Ywater, nucl_rate, nucl_rmin, conc_sat_spec, ff_nucl,ff_cond
! PENCILS PROVIDED latent_heat
!
!***************************************************************
module Chemistry
!
  use Cdata
  use General, only: keep_compiler_quiet, itoa
  use EquationOfState, only: cp_const,lpres_grad,imass,&
                             getdensity,gettemperature,getpressure

  use Messages, only: svn_id,timing,fatal_error,inevitably_fatal_error,warning,not_implemented
!
  implicit none
!
  include 'chemistry.h'
!
  real :: Rgas, Rgas_unit_sys=1.
  real, dimension(mx,my,mz) :: cp_full, cv_full
  real, dimension(:,:,:), pointer :: mu1_full
  real, dimension(mx,my,mz) :: lambda_full, rho_full, TT_full
  real, dimension(mx,my,mz,nchemspec) :: cv_R_spec_full
!real, dimension (mx,my,mz) ::  e_int_full,cp_R_spec

  real, dimension(mx,my,mz,nchemspec) :: cp_spec_glo

! parameters for simplified cases
  real :: lambda_const=impossible
  real :: visc_const=impossible
  real :: Diff_coef_const=impossible
  real :: Sc_number=0.7!, Pr_number=0.7
!  real :: Cp_const=impossible
  real :: Cv_const=impossible
  logical :: lfix_Sc=.false., lfix_Pr=.false.
  logical :: init_from_file, reinitialize_chemistry=.false.
  logical :: lnucleation, lcorr_vel=.false., lnucl_dynamic=.false.
  logical :: lchem_detailed=.true.
  logical :: lgradP_terms=.true.
  character(len=30) :: reac_rate_method = 'chemkin'
  character(len=2) :: mixture_fraction_element
! --AB--  character(len=30) :: inucl_pre_exp="const"
! parameters for initial conditions
  real :: init_x1=-0.2, init_x2=0.2
  real :: init_y1=-0.2, init_y2=0.2
  real :: init_z1=-0.2, init_z2=0.2
  real :: init_TT1=298., init_TT2=2400., init_ux=0., init_uy=0., init_uz=0.
  real :: init_zz1=0.01, init_zz2=0.2
  real :: flame_pos=0.
  real :: init_rho2=1.
  real :: init_rho=1.
  real :: str_thick=0.02
  real :: init_pressure=10.13e5
  real :: global_phi=impossible
  real :: delta_chem=0., press=0.
!
  logical :: lone_spec=.false., lfilter_strict=.false.
!
!  parameters related to chemical reactions, diffusion and advection
!
  logical :: lreactions=.true.
  logical :: ladvection=.true.
  logical :: ldiffusion=.true.
!  logical :: lnospec_eqns=.true.
!
  logical :: lheatc_chemistry=.true.
  logical :: lspecies_cond_simplified=.false.
  logical :: lDiff_simple=.false.
  logical :: lDiff_lewis=.false.
  logical :: lThCond_simple=.false.
  logical :: lT_const=.false.
  logical :: lDiff_fick=.false.
  logical :: lFlux_simple=.false.
  logical :: ldiff_corr=.false.
  logical, save :: tran_exist=.false.
  logical, save :: lew_exist=.false.
  logical :: lSmag_heat_transport=.false.
  logical :: lSmag_diffusion=.false.
  logical :: lnormalize_chemspec=.false., lnormalize_chemspec_N2=.false.
  logical :: lFlame_index_as_aux=.false., lmixture_fraction_as_aux=.false.
  !
  logical :: lfilter=.false., lupw_chemspec=.false.
  logical :: lkreactions_profile=.false., lkreactions_alpha=.false.
  integer :: nreactions=0, nreactions1=0, nreactions2=0
  integer :: ll1, ll2, mm1, mm2, nn1, nn2
  real, allocatable, dimension(:) :: kreactions_profile_width, kreactions_alpha
  character(len=5) :: flameind_spec1="H2"
  character(len=5) :: flameind_spec2="CO2"
!
  integer :: mreactions, iadv=0
  real, allocatable, dimension(:,:) :: vreactions_p,vreactions_m
  !$omp threadprivate(vreactions_p,vreactions_m)
!
!  The stociometric factors need to be reals for arbitrary reaction orders
!
  real, allocatable, dimension(:,:) :: stoichio, Sijm, Sijp, Sijm_, Sijp_, Sijm_mod, Sijp_mod
  real, allocatable, dimension(:,:) :: kreactions_z
  real, allocatable, dimension(:) :: kreactions_m, kreactions_p
  logical, allocatable, dimension(:) :: back
  logical, allocatable, dimension(:,:,:) :: lnucleii_generated
  character(len=30), allocatable, dimension(:) :: reaction_name
  character (len=labellen) :: self_collisions='nothing',condensing_species='nothing'
  logical :: lT_tanh=.false.
  logical :: ldamp_zone_for_NSCBC=.false.
  logical :: linit_temperature=.false., linit_density=.false.
  logical :: lreac_as_aux=.false.
  integer :: ifuel_flow=0.
!
! 1step_test case
! possible, will be removed later
!
  logical :: l1step_test=.false.
  logical :: lflame_front=.false., ltriple_flame=.false.
  logical :: lFlameMaster=.false.
  integer :: ipr=2
  real :: Tc=440., Tinf=2000., beta=1.09
  real :: z_cloud=0., mix_frac_IH=0.
!
!  hydro-related parameters
!
  real, dimension(nchemspec) :: amplchemk=0., amplchemk2=0.
  real, dimension(nchemspec) :: chem_diff_prefactor=1., initial_massfractions
  real :: amplchem=1., kx_chem=1., ky_chem=1., kz_chem=1., widthchem=1.
  real :: chem_diff=0.
  character(len=labellen), dimension(ninit) :: initchem='nothing'
  character(len=labellen), allocatable, dimension(:) :: kreactions_profile
  character(len=60) :: prerun_directory='nothing'
  character(len=60) :: file_name='nothing'
!
  real, allocatable, dimension(:,:,:,:,:) :: Bin_Diff_coef
  real, allocatable, dimension(:,:,:,:) :: Diff_full, Diff_full_add
  real, dimension(mx,my,mz,nchemspec) :: XX_full
  real, dimension(mx,my,mz,nchemspec) :: species_viscosity
  real, dimension(mx,my,mz,nchemspec) :: RHS_Y_full
  real, dimension(nchemspec) :: nu_spec=0., mobility=1.
!
!  Chemkin related parameters
!
  logical :: lcheminp=.false., lchem_cdtc=.false.
  logical :: lmobility=.false.
  integer :: iTemp1=2, iTemp2=3, iTemp3=4
  integer, parameter :: iaa1_offset=4,iaa2_offset=11
  real, allocatable, dimension(:) :: B_n, alpha_n, E_an
  real, allocatable, dimension(:,:) :: low_coeff, high_coeff, troe_coeff, a_k4
  real, allocatable, dimension(:) :: low_coeff_abs_max, high_coeff_abs_max, troe_coeff_abs_max, a_k4_min
  logical, allocatable, dimension(:) :: Mplus_case
  logical, allocatable, dimension(:) :: photochem_case
  real :: lamb_low, lamb_up, Pr_turb=0.7
  !
  ! Condensing species parameters
  !
  integer :: i_cond_spec,ichem_cond_spec
  integer :: imassH=19, imassO=20, imassC=21, imassN=22, imassS=23, imassTI=24
  integer, pointer :: it_insert_nuclei, Ntau
  real, pointer :: redfrac
  real :: true_density_cond_spec_cgs=2.196, true_density_cond_spec
  real :: gam_surf_energy_cgs=32.
  real :: nucleation_rate_coeff_cgs=1e19
  real :: molar_mass_spec, atomic_m_spec, A_spec
  real :: deltaH_cgs = 7.07e12 ! in erg/mol
  !real :: deltaH_cgs = 3.55e12 ! in erg/mol (based on Tboil and 1025C)
  real :: gam_surf_energy_mul_fac=1.0
  real :: conc_sat_spec_cgs=1e-8 !units of mol/cm^3
  real :: min_nucl_radius_cgs=1e-8 ! units of cm
  logical, pointer :: ldustnucleation, lpartnucleation, lcondensing_species
  character(len=labellen) :: isurf_energy="const"
  character(len=labellen) :: iconc_sat_spec="const"
  character(len=30) :: inucl_pre_exp="const"      
  logical :: lnoevap=.false., lnolatentheat=.true.
  !
  character(len=30), dimension(nchemspec) :: init_premixed_fuel="nothing"
  real, dimension(nchemspec) :: init_fuel_molar_ratio=0.
  real, dimension(nchemspec) :: init_fuel_O2_demand=0.
  real :: init_temp_fuel=0., init_temp_oxidizer=0., init_phi=0.
!
!   Atmospheric physics
!
  logical :: latmchem=.false.
  logical :: lcloud=.false.
  integer, SAVE :: index_O2=0., index_N2=0., index_O2N2=0., index_H2O=0.
!
!   Species constants
!
  real, dimension(nchemspec,24), target :: species_constants
!
!   Lewis coefficients
!
 real, dimension(nchemspec) :: Lewis_coef, Lewis_coef1
!
!   Transport data
!
 real, dimension(nchemspec,7) :: tran_data
!
!   Diagnostics
!
  real, allocatable, dimension(:,:), target :: net_react_m, net_react_p
  !$omp threadprivate(net_react_m,net_react_p)
! For concurrency
  real, dimension(:,:), pointer :: p_net_react_m, p_net_react_p
!
  real, dimension(nchemspec) :: Ythresh=0.
  logical :: lchemistry_diag=.false.
!
!   Hot spot problem
!
  logical :: lhotspot=.false.
!
! mmx is 1 if nx=1, and it is mx if nx > 1
!
  integer :: mmx=2*nghost*min(1,nx-1)+nx
!
  integer, dimension(nchemspec) :: ireaci=0
!
! input parameters
  namelist /chemistry_init_pars/ &
      initchem, amplchem, kx_chem, ky_chem, kz_chem, widthchem, &
      amplchemk,amplchemk2, chem_diff,nu_spec,lDiff_simple, &
      lDiff_lewis,lFlux_simple, &
      lThCond_simple,lambda_const, visc_const,Cp_const,Cv_const,Diff_coef_const, &
      init_x1,init_x2,init_y1,init_y2,init_z1,init_z2,init_TT1,&
      init_TT2,init_rho, &
      init_ux,init_uy,init_uz,l1step_test,Sc_number,init_pressure,lfix_Sc, &
      str_thick,lfix_Pr, lT_tanh,lT_const, lheatc_chemistry, lspecies_cond_simplified, &
      ldamp_zone_for_NSCBC, latmchem, lcloud, prerun_directory, &
      lchemistry_diag,lfilter_strict,linit_temperature, &
      linit_density, init_rho2, &
      file_name, lreac_as_aux, init_zz1, init_zz2, flame_pos, &
      reac_rate_method,global_phi, lSmag_heat_transport, Pr_turb, lSmag_diffusion, z_cloud, &
      lhotspot, lchem_detailed, condensing_species, conc_sat_spec_cgs, &
      true_density_cond_spec_cgs, delta_chem, press, init_premixed_fuel, &
      init_fuel_molar_ratio,init_fuel_O2_demand,init_temp_fuel, init_temp_oxidizer, init_phi, &
      lFlame_index_as_aux, lmixture_fraction_as_aux, mixture_fraction_element, ifuel_flow, &
      flameind_spec1, flameind_spec2, lnucl_dynamic
!
!
! run parameters
  namelist /chemistry_run_pars/ &
      lkreactions_profile, lkreactions_alpha, &
      chem_diff,chem_diff_prefactor, nu_spec, ldiffusion, ladvection, &
      lreactions, lchem_cdtc, lheatc_chemistry, lspecies_cond_simplified, lchemistry_diag, &
      lmobility,mobility, lfilter,lT_tanh,lDiff_simple,lDiff_lewis,lFlux_simple, &
      lThCond_simple,visc_const,cp_const,reinitialize_chemistry,init_from_file, &
      lfilter_strict,init_TT1,init_TT2,init_x1,init_x2, linit_temperature, &
      linit_density, &
      ldiff_corr, lDiff_fick, lreac_as_aux, reac_rate_method,global_phi, &
      Ythresh, lchem_detailed, conc_sat_spec_cgs, inucl_pre_exp, lcorr_vel, &
      lgradP_terms, lnormalize_chemspec, lnormalize_chemspec_N2, &
      gam_surf_energy_cgs, isurf_energy, iconc_sat_spec, nucleation_rate_coeff_cgs, &
      lnoevap, lnolatentheat, gam_surf_energy_mul_fac, deltaH_cgs,&
      min_nucl_radius_cgs, lupw_chemspec, &
      lFlame_index_as_aux, lmixture_fraction_as_aux, mixture_fraction_element, &
      flameind_spec1, flameind_spec2
!
! diagnostic variables (need to be consistent with reset list below)
!
  integer, dimension(nchemspec) :: idiag_Ym=0      ! DIAG_DOC: $\left<Y_x\right>$
  integer, dimension(nchemspec) :: idiag_rhoYm=0   ! DIAG_DOC: $\left<\rho Y_x\right>$
  integer, dimension(nchemspec) :: idiag_TYm=0     ! DIAG_DOC: $\left<Y_{\rm thresh}-Y_x\right>$
  integer, dimension(nchemspec) :: idiag_dYm=0     ! DIAG_DOC: $\delta\left<Y_x\right>/\delta t$
  integer, dimension(nchemspec) :: idiag_dYmax=0   ! DIAG_DOC: $max\delta\left<Y_x\right>/\delta t$
  integer, dimension(nchemspec) :: idiag_Ymax=0    ! DIAG_DOC: $\left<Y_{x,max}\right>$
  integer, dimension(nchemspec) :: idiag_Ymin=0    ! DIAG_DOC: $\left<Y_{x,min}\right>$
  integer, dimension(nchemspec) :: idiag_hm=0      ! DIAG_DOC: $\left<H_{x,max}\right>$
  integer, dimension(nchemspec) :: idiag_cpm=0     ! DIAG_DOC: $\left<c_{p,x}\right>$
  integer, dimension(nchemspec) :: idiag_diffm=0   ! DIAG_DOC: $\left<D_{x}\right>$
  integer, dimension(nchemspec) :: idiag_diffmax=0 ! DIAG_DOC: $\left<D_{x,max}\right>$
  integer, dimension(nchemspec) :: idiag_diffmin=0 ! DIAG_DOC: $\left<D_{x,min}\right>$
  integer, dimension(nchemspec) :: idiag_Ymz=0     ! DIAG_DOC: $\left<Y_x\right>_{xy}(z)$
  integer :: idiag_dtchem=0     ! DIAG_DOC: $dt_{chem}$
  integer :: idiag_nuclrmin=0! DIAG_DOC: $\left< r_{\min} \right>$
  integer :: idiag_nuclrate=0! DIAG_DOC: $\left< J \right>$
  integer :: idiag_conc_satm=0
!
  integer :: idiag_cpfull=0
  integer :: idiag_cvfull=0
  integer :: idiag_e_intm=0
!
  integer :: idiag_lambdam=0,idiag_lambdamax=0,idiag_lambdamin=0
  integer :: idiag_alpham=0,idiag_alphamax=0,idiag_alphamin=0
  integer :: idiag_ffnucl=0, idiag_supersat=0, idiag_latent_heat=0
!
  integer :: idiag_mixfracmax=0
  integer :: idiag_flameindmax=0
  integer :: idiag_mixfracmin=0
  integer :: idiag_flameindmin=0
!
!  Auxiliaries.
!
  integer :: ireac=0
!
  integer :: enum_reac_rate_method = 0
  integer, dimension(:), allocatable :: enum_reaction_name
  integer :: enum_iconc_sat_spec = 0
  integer :: enum_isurf_energy = 0
  integer :: enum_inucl_pre_exp = 0

  integer :: i_O2_glob,ichem_O2, i_C3H8_glob,ichem_C3H8
  logical :: lO2, lC3H8
  real    :: mO2, mC3H8
  contains
!
!***********************************************************************
    subroutine register_chemistry
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  13-aug-07/steveb: coded
!   8-jan-08/axel: added modifications analogously to dustdensity
!   5-mar-08/nils: Read thermodynamical data from chem.inp
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
      integer :: k, ichemspec_tmp, stat
      character(len=fnlen) :: input_file
      logical ::  chemin, cheminp
!
!  Set ichemistry to consecutive numbers nvar+1, nvar+2, ..., nvar+nchemspec.
!
      call farray_register_pde('chemspec',ichemspec_tmp,array=nchemspec)
      do k = 1,nchemspec
        ichemspec(k) = ichemspec_tmp+k-1
      enddo
!
!  Register viscosity
!
      call farray_register_auxiliary('viscosity',iviscosity,communicated=.false.)
!
!  Writing files for use with IDL
!
      if (naux+naux_com <  maux+maux_com) aux_var(aux_count) = ',viscosity $'
      if (naux+naux_com == maux+maux_com) aux_var(aux_count) = ',viscosity'
      aux_count = aux_count+1
      if (lroot) write (4,*) ',visocsity $'
      if (lroot) write (15,*) 'viscosity = fltarr(mx,my,mz)*one'
!
!  Read species to be used from chem.inp (if the file exists).
!
      inquire (FILE='chem.inp', EXIST=lcheminp)
      if (.not. lcheminp) inquire (FILE='chem.in', EXIST=lcheminp)

      inquire (FILE='chem.inp', EXIST=cheminp)
      inquire (FILE='chem.in', EXIST=chemin)
      if (cheminp .and. chemin) call fatal_error('chemistry', &
          'chem.inp and chem.in found, please decide for one')
      if (cheminp) input_file='chem.inp'
      if (chemin) input_file='chem.in'
!
      if (lcheminp) call read_species(input_file)
!
!  Read data on the thermodynamical properties of the different species.
!  All these data are stored in the array species_constants.
!
      if (lcheminp) call read_thermodyn(input_file)
!
!  Write all data on species and their thermodynamics to file.
!
      if (lcheminp) call write_thermodyn
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( "$Id$")
!
      call put_shared_variable('species_constants',species_constants,caller='register_chemistry')
      if (lparticles) then
        call put_shared_variable('true_density_cond_spec',true_density_cond_spec)
        call put_shared_variable('lnucl_dynamic',lnucl_dynamic,caller='register_chemistry')
        if (lnucl_dynamic .and. npscalar < 3) then
          call fatal_error("register_chemistry",&
               "npscalar must be at least 3 if lnucl_dynamic is true.")
        endif
        if (lnucl_dynamic) then
          if (.not. allocated(lnucleii_generated)) then
            allocate(lnucleii_generated(nx,ny,nz),STAT=stat)
            if (stat > 0) call fatal_error("register_chemistry","Couldn't allocate lnucleii_generated!")
          endif
!          call put_shared_variable('lnucleii_generated',lnucl_dynamic,caller='register_chemistry')
        endif
      endif
!
    endsubroutine register_chemistry
!***********************************************************************
    subroutine initialize_chemistry(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  13-aug-07/steveb: coded
!  19-feb-08/axel: reads in chemistry.dat file
!  21-nov-10/julien: added the reaction rates as optional auxiliary variables
!                    in the f array for output.
!
      use FArrayManager
      use SharedVariables, only: get_shared_variable
      use Messages, only: warning
!
      real, dimension(mx,my,mz,mfarray) :: f
      logical :: data_file_exit=.false.
      logical :: exist, exist1, exist2
      character(len=15) :: file1='chemistry_m.dat', file2='chemistry_p.dat'
      integer :: ind_glob, ind_chem, stat
      integer :: i, ichem, imass_spec
      logical :: found_specie=.false.
      character(len=2) :: car2
      real, dimension(nchemspec) :: Y_init_mix_frac=0.
!
!  initialize chemistry
      !
      if (lcheminp) then
        if (unit_temperature /= 1) &
             call fatal_error('initialize_chemistry','unit_temperature must be unity for chemistry')
      endif
!
!  Find indices for oxygen and propane and their masses
!
      call find_species_index('O2',i_O2_glob,ichem_O2,lO2)
      call find_species_index('C3H8',i_C3H8_glob,ichem_C3H8,lC3H8)
      if (lO2)   mO2 = species_constants(ichem_O2,imass)
      if (lC3H8) mC3H8 = species_constants(ichem_C3H8,imass)
!
!  calculate universal gas constant based on Boltzmann constant
!  and the proton mass
!
      if (unit_system == 'cgs') then
        Rgas_unit_sys = k_B_cgs/m_u_cgs
        Rgas = Rgas_unit_sys/unit_energy
      endif
!
      if (nchemspec == 1) then
        lone_spec = .true.
        lreactions = .false.
        ladvection = .false.
        ldiffusion = .false.
      endif
!
!  check for the existence of chemistry input files
!
      inquire (file=file1,exist=exist1)
      inquire (file=file2,exist=exist2)
      inquire (file='chemistry.dat',exist=exist)
!
!  Read in data file in ChemKin format
!
      if (lcheminp) then
        call chemkin_data(f)
        data_file_exit = .true.
        if (latmchem) then
          call find_species_index('O2',ind_glob,ind_chem,found_specie)
          if (found_specie) then
            index_O2 = ind_chem
          else
            call fatal_error('initialize_chemistry', 'no O2 has been found')
          endif
          call find_species_index('N2',ind_glob,ind_chem,found_specie)
          if (found_specie) then
            index_N2 = ind_chem
          else
            call fatal_error('initialize_chemistry', 'no N2 has been found')
          endif
          call find_species_index('O2N2',ind_glob,ind_chem,found_specie)
          if (found_specie) then
            index_O2N2 = ind_chem
          else
            call fatal_error('initialize_chemistry', 'no O2N2 has been found')
          endif
          call find_species_index('H2O',ind_glob,ind_chem,found_specie)
          if (found_specie) then
            index_H2O = ind_chem
          else
            call fatal_error('initialize_chemistry', 'no H2O has been found')
          endif
        endif
        if (lcloud) then
          call find_species_index('H2O',ind_glob,ind_chem,found_specie)
          if (found_specie) then
            index_H2O = ind_chem
          else
            call fatal_error('initialize_chemistry', 'no H2O has been found')
          endif
        endif
        if (ldustdensity .or. lparticles_radius) then
          call get_shared_variable('lcondensing_species',lcondensing_species)

          if (lcondensing_species .and. (lparticles_radius .or. ldustdensity)) then
            call find_species_index(condensing_species,i_cond_spec,ichem_cond_spec,found_specie)
            if (.not. found_specie) then
              print*,"condensing_species=",condensing_species
              call fatal_error('initialize_chemistry','there is no such species')
            else
              if (lroot) then
                print*,"Condensing species:"
                print*,"i_cond_spec=",i_cond_spec
                print*,"ichem_cond_spec=",ichem_cond_spec
                print*,"species=",varname(ichemspec(ichem_cond_spec))
              endif
            endif
          endif
        endif

        if (lfix_Sc .and. any(idiag_diffm/=0)) then
          call warning('initialize_chemistry', &
                       'diagnostic diffm not available for lfix_Sc=T - switched off')
          idiag_diffm=0
        endif
      endif
!
!  Alternatively, read in stoichiometric matrices in explicit format.
!  For historical reasons this is referred to as "astrobiology_data"
!
      if (exist1.and.exist2 .or. exist) then
        call astrobiology_data
        data_file_exit = .true.
      endif
!
! check the existence of a data file
!
      if (.not. data_file_exit) &
        call fatal_error('initialize_chemistry','there is no chemistry data file')
!
      if (dimensionality==0) then
        ll1 = l1
        ll2 = l2
        mm1 = m1
        mm2 = m2
        nn1 = n1
        nn2 = n2
      else
        if (nxgrid == 1) then
          ll1 = l1
          ll2 = l2
        else
          ll1 = 1
          ll2 = mx
        endif
!
        if (nygrid == 1) then
          mm1 = m1
          mm2 = m2
        else
          mm1 = 1
          mm2 = my
        endif
!
        if (nzgrid == 1) then
          nn1 = n1
          nn2 = n2
        else
          nn1 = 1
          nn2 = mz
        endif
      endif
!
!  Reinitialize if required
!
      if (reinitialize_chemistry) then
        if (lroot) print*,'Reinitializing chemistry.'
        call init_chemistry(f)
      endif
!
!  Find initial mass fractions
!
      if (reac_rate_method=='roux') then
        do i = 1, nchemspec
          initial_massfractions(i)=f(l1,m1,n1,ichemspec(i))
        enddo
!
        if (lroot) print*,'initial_massfractions=',initial_massfractions
      endif
!
!  allocate memory for net_reaction diagnostics
!
      if (allocated(net_react_p) .and. .not. lreloading) print*, 'this should not be here'

      if (lchemistry_diag .and. .not. lreloading) then
        allocate(net_react_p(nchemspec,nreactions),STAT=stat)
        if (stat > 0) call fatal_error("initialize_chemistry","Couldn't allocate net_react_p")
        net_react_p = 0.
        allocate(net_react_m(nchemspec,nreactions),STAT=stat)
        if (stat > 0) call fatal_error("initialize_chemistry","Couldn't allocate net_react_m")
        net_react_m = 0.
      endif
!
!  Define the chemical reaction rates as auxiliary variables for output
!
      if (lreac_as_aux) then
        if (ireac == 0) then
          call farray_register_auxiliary('reac',ireac,array=nchemspec)
          do i = 0, nchemspec-1
            ireaci(i+1) = ireac+i
          enddo
        else
          if (lroot) print*, 'initialize_chemistry: ireac = ', ireac
          call farray_index_append('ireac',ireac)
          do i = 1, nchemspec
            write (car2,'(i2)') i
            call farray_index_append('ireac'//trim(adjustl(car2)),ireaci(i))
          enddo
        endif
      endif
!
! Allocate space for flame index and mixture fraction
!
      if (lFlame_index_as_aux) then
        call farray_register_auxiliary('flameind',iflameind,communicated=.false.)
      endif
      if (lmixture_fraction_as_aux) then
        call farray_register_auxiliary('mixfrac',imixfrac,communicated=.false.)
!
! Define initial mass fraction used to find mixture fractions
!
        select case (initchem(1))
        case ("double_shear_layer", "double_shear_layer_x")
           if (ifuel_flow == 1) then
              Y_init_mix_frac=amplchemk
           elseif (ifuel_flow == 2) then
              Y_init_mix_frac=amplchemk2
           else
              call fatal_error("initialize_chemistry","Please specify ifuel_flow!")
           endif
        case default
           call fatal_error("initialize_chemistry","initial mass fractions are not defined for mixt. frac.")
        end select
!
! Initialize variable used for calculation of Bilger mixture fraction based on Hydrogen
!
        mix_frac_IH=0.
        select case (mixture_fraction_element)
        case ("H")
          imass_spec=imassH
        case ("C")
          imass_spec=imassC
        case ("O")
          imass_spec=imassO
        case ("S")
           imass_spec=imassS
        case ("TI")
           imass_spec=imassTI
        case default
          print*,"mixture_fraction_element=",mixture_fraction_element
          call fatal_error("make_mixture_fraction","No such mixture_fraction_element")
        end select
        do ichem=1,nchemspec
          mix_frac_IH=mix_frac_IH+Y_init_mix_frac(ichem)*species_constants(ichem,imass_spec)/species_constants(ichem,imass)
       enddo
      endif
      !
      if (leos) then
        call get_shared_variable('mu1_full',mu1_full,caller='initialize_chemistry')
      else
        call warning('initialize_chemistry','mu1_full not provided by eos')
      endif
!
!  04-18-11/Julien: Modified the computation of simple heat conductivity according
!                   to Smooke & Giovangigli 1991. lambda_const is now equal
!                   to 2.58e-4 instead of 1e4, and represents the ratio \lambda0/cp0.
!                   Formula:
!                   \lambda = lambda_const*cp*(T/T0)**0.7, with T0 now = 298K
!
      if (lThCond_simple .and. lambda_const == impossible) lambda_const = 2.58e-4
!MR: Is it guaranteed that lambda_const is nowhere used with value impossible?
!
!  04-18-11/Julien: Changed the value of Diff_coef_const from 10 to 2.58e-4
!                   according to Smooke & Giovangigli 1991. Now Diff_coef_const
!                   does not represent a constant diffusion coefficient,
!                   but \rho0 D0.
!                   Diff_penc_add are still the diffusion coefficients. Formula:
!                   D = Diff_coef_const/\rho*(T/T0)**n0.7, with T0 now = 298K.
!
      if (lDiff_simple .and. Diff_coef_const == impossible) Diff_coef_const = 2.58e-4  !MR: the same as lambda_const?
!
!  true_density_cond_spec
!
      if (ldustdensity) then
        true_density_cond_spec=true_density_cond_spec_cgs/unit_density
        call get_shared_variable('ldustnucleation',ldustnucleation)
        lnucleation=ldustnucleation
      endif
      if (lparticles) then
        call get_shared_variable('lpartnucleation',lpartnucleation)
        lnucleation=lpartnucleation
        if  (lparticles_radius) then
          true_density_cond_spec=true_density_cond_spec_cgs/unit_density
        endif
        call get_shared_variable('it_insert_nuclei',it_insert_nuclei)
        call get_shared_variable('Ntau',Ntau)
        call get_shared_variable('redfrac',redfrac)
      endif
!
! Define some constants used for condensing species
!
      if (lcondensing_species .and. (lparticles_radius .or. ldustdensity)) then
        molar_mass_spec = species_constants(ichem_cond_spec,imass)
        atomic_m_spec=molar_mass_spec*m_u
        A_spec=sqrt(8.*k_B/(pi*atomic_m_spec))*molar_mass_spec/(4.*true_density_cond_spec)
      endif
!
! Do we want upwinding
!
      if (lupw_chemspec) iadv=1
!
!  write array dimension to chemistry diagnostics file
!
      open (1,file=trim(datadir)//'/net_reactions.dat',position='append')
      write (1,*) nchemspec,nreactions
      close (1)
!
    if(allocated(a_k4))       a_k4_min = minval(a_k4,1)
    if(allocated(low_coeff))  low_coeff_abs_max = maxval(abs(low_coeff),1)
    if(allocated(high_coeff)) high_coeff_abs_max = maxval(abs(high_coeff),1)
    if(allocated(troe_coeff)) troe_coeff_abs_max = maxval(abs(troe_coeff),1)
    if(.not. lmultithread) call chemistry_allocate_rhs_arrays

    endsubroutine initialize_chemistry
!***********************************************************************
    subroutine chemistry_allocate_rhs_arrays

    if(.not. allocated(vreactions_p)) allocate(vreactions_p(nx,mreactions))
    if(.not. allocated(vreactions_m)) allocate(vreactions_m(nx,mreactions))

    endsubroutine chemistry_allocate_rhs_arrays
!***********************************************************************
    subroutine init_chemistry(f)
!
!  initialise chemistry initial condition; called from start.f90
!
!  13-aug-07/steveb: coded
!     jul-10/julien: Added some new initial cases
!
      use Initcond
      use InitialCondition, only: initial_condition_chemistry
!
      real, dimension(mx,my,mz,mfarray) :: f
      real :: PP, prof, der, tmp1, tmp2, mol1, mol2, mean_molar_mass
      real :: rho, TT, Ysum
      integer :: l,j,k
      logical :: lnothing, air_exist
!
      intent(inout) :: f
!
!  different initializations of nd (called from start)
!
      lnothing = .false.
      do j = 1,ninit
        select case (initchem(j))
!
        case ('nothing')
          if (lroot .and. .not. lnothing) print*, 'initchem: nothing '
          lnothing = .true.
        case ('constant')
          do k = 1,nchemspec
            f(:,:,:,ichemspec(k)) = f(:,:,:,ichemspec(k)) + amplchemk(k)
            print*,k,f(:,:,:,ichemspec(k))
          enddo
        case ('positive-noise')
          do k = 1,nchemspec
            call posnoise(amplchemk(k),f,ichemspec(k))
          enddo
        case ('positive-noise-rel')
          do k = 1,nchemspec
            call posnoise_rel(amplchemk(k),amplchemk2(k),f,ichemspec(k))
          enddo
        case ('innerbox')
          do k = 1,nchemspec
            call innerbox(amplchemk(k),amplchemk2(k),f,ichemspec(k),widthchem)
          enddo
        case ('cos2x_cos2y_cos2z')
          do k = 1,nchemspec
            call cos2x_cos2y_cos2z(amplchemk(k),f,ichemspec(k))
          enddo
        case ('coswave-x')
          do k = 1,nchemspec
            call coswave(amplchem,f,ichemspec(k),kx=kx_chem)
          enddo
        case ('gaussian-x')
          do k = 1,nchemspec
            call gaussian(amplchem,f,ichemspec(k),kx=kx_chem)
          enddo
        case ('gaussian-pos')
          do k = 1,nchemspec
            call gaussianpos(amplchemk(k),f,ichemspec(k),widthchem,kx_chem,ky_chem,kz_chem)
            print*,"c(",x(l1),",",y(m1),",",z(n1),")=", f(l1,m1,n1,ichemspec(k))
          enddo
        case ('hatwave-x')
          do k = 1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),widthchem,kx=kx_chem,pos=1.)
          enddo
        case ('hatwave-y')
          do k = 1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),widthchem,ky=ky_chem)
          enddo
        case ('hatwave-z')
          do k = 1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),widthchem,kz=kz_chem)
          enddo
        case ('air')
          if (lroot ) print*, 'init_chem: air '
          inquire (file='air.dat',exist=air_exist)
          if (.not. air_exist) inquire (file='air.in',exist=air_exist)
          if (air_exist) then
            call air_field(f,PP)
          else
            call fatal_error('init_chemistry','no air.in or air.dat file found')
          endif
        case ('flame_front')
          call flame_front(f)
        case ('TTD')
          call TTD(f)
        case ('triple_flame')
          call triple_flame(f)
        case ('flame')
          if (lroot) print*, 'initchem: flame '
          call flame(f)
        case ('flame_blob')
          call flame_blob(f)
        case ('opposite_flames')
          call opposite_flames(f)
        case ('opposite_ignitions')
          call opposite_ignitions(f)
        case ('prerun_1D')
          call prerun_1D(f,prerun_directory)
        case ('prerun_1D_opp')
          call prerun_1D_opp(f,prerun_directory)
        case ('FlameMaster')
          call FlameMaster_ini(f,file_name)
        case ('flame_front_new')
          call flame_front_new(f)
        case ('premixed_equiv_ratio')
          call premixed_equiv_ratio(f)
        case ('double_shear_layer', 'double_shear_layer_x')
!
!  Double shear layer
!
          tmp1=0.
          tmp2=0.
          if (lroot) then
            print*,'init_chemistry: Double shear layer'
            write (*,'(A10, 2A23)')"species","mass frac jet","mass frac co-flow"
            print*,"********************COMPOSITIONS****************************"
          endif
          do k=1,nchemspec
            if (lroot) write (*,'(A10,2F23.6)') trim(varname(ichemspec(k))), &
                 amplchemk(k),amplchemk2(k)
            tmp1=tmp1+amplchemk(k)  /species_constants(k,imass)
            tmp2=tmp2+amplchemk2(k)/species_constants(k,imass)
          enddo
          if (lroot) write (*,'(A10,2F23.6)') "SUM",sum(amplchemk),sum(amplchemk2)
          
          if (lroot) print*," "
          
          if (lroot) then
            write (*,'(A10, 2A23)')"species","mol frac jet","mol frac co-flow"
            print*,"********************COMPOSITIONS****************************"
          endif
          do k=1,nchemspec
            mol1=amplchemk(k)/species_constants(k,imass)/tmp1
            mol2=amplchemk2(k)/species_constants(k,imass)/tmp2
            if (lroot) write (*,'(A10,2F23.6)') trim(varname(ichemspec(k))), &
                 mol1,mol2
          enddo
          if (lroot) write (*,'(A10,2F23.6)') "SUM",sum(amplchemk),sum(amplchemk2)
!
!  Different cases for x and y extent.
!
          if (initchem(j)=='double_shear_layer_x') then
            do l=1,mx
              der=2./delta_chem
              prof=(tanh(der*(x(l)+widthchem))-   &
                    tanh(der*(x(l)-widthchem)))/2.
              do k=1,nchemspec
                do n=1,mz
                do m=1,my
                  f(l,m,n,ichemspec(k))=amplchemk(k)*prof +amplchemk2(k)*(1-prof) 
                enddo
                enddo
              enddo
            enddo
          else
            do m=1,my
              der=2./delta_chem
              prof=(tanh(der*(y(m)+widthchem))-   &
                   tanh(der*(y(m)-widthchem)))/2.
              do k=1,nchemspec
                do n=1,mz
                  f(1:mx,m,n,ichemspec(k))=amplchemk(k)*prof +amplchemk2(k)*(1-prof) 
                enddo
              enddo
            enddo
          endif
          !
          call getmu_array(f,mu1_full,linit=.true.)
          !
          ! Must also set density such that the pressure is correct
          !
          do n=1,mz
            do m=1,my
              do l=1,mx
                !tmp1=0
                !do k=1,nchemspec
                !  tmp1=tmp1+f(l,m,n,ichemspec(k))/species_constants(k,imass)
                !enddo
                mean_molar_mass=1./tmp1
                if (ltemperature_nolog) then
                  TT=f(l,m,n,iTT)
                else
                  TT=exp(f(l,m,n,ilnTT))
                endif
                rho=press/(mu1_full(l,m,n)*Rgas*TT)
                if (ldensity_nolog) then
                  f(l,m,n,irho)=rho
                else
                  f(l,m,n,ilnrho)=alog(rho)
                endif

                Ysum=sum(f(l,m,n,ichemspec(1:nchemspec)))
                if ((Ysum  > 1+1e-3) .or. (Ysum < 1-1e-3)) then
                  print*,"l,m,n,Ysum=",l,m,n,Ysum
                  call fatal_error("init_chemistry","The sum of all mass fractions should be unit!")
                else
                  ! Normalize the sum of all species to unity
                  f(l,m,n,ichemspec(1:nchemspec))=&
                       f(l,m,n,ichemspec(1:nchemspec))/Ysum
                endif
              enddo
            enddo
          enddo
!
        case default
!
!  Catch unknown values
!
          call fatal_error('init_chemistry','no such initchem: '//trim(initchem(j)))
        endselect
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_chemistry(f)
!
!   The following lines are kept temporally
!
!      if (lone_spec) then
!        f(:,:,:,ichemspec(1))=1.
!        if (lroot) print*, 'initchem: this is one specie case'
!      endif
!
    endsubroutine init_chemistry    
!***********************************************************************
    subroutine pencil_criteria_chemistry
!
!  All pencils that this chemistry module depends on are specified here.
!
!  13-aug-07/steveb: coded
      !
      lpenc_requested(i_gXXk) = .true.
      lpenc_requested(i_gYYk) = .true.
!      if (lreactions)
      lpenc_requested(i_ghhk) = .true.
!
      if (lreactions) lpenc_requested(i_DYDt_reac) = .true.
      lpenc_requested(i_DYDt_diff) = .true.
!
      if (lcheminp) then
        lpenc_requested(i_rho) = .true.
        lpenc_requested(i_cv) = .true.
        lpenc_requested(i_cp) = .true.
        lpenc_requested(i_cv1) = .true.
        lpenc_requested(i_cp1) = .true.
!         if (lreactions)
        lpenc_requested(i_H0_RT) = .true.
!         if (lreactions)
        lpenc_requested(i_S0_R) = .true.
        lpenc_requested(i_nu) = .true.
        lpenc_requested(i_gradnu) = .true.
        lpenc_requested(i_cs2) = .true.
!
!         if (lreactions)
        lpenc_requested(i_hhk_full) = .true.
        if (lThCond_simple) lpenc_requested(i_glncp) = .true.
!
        if (lheatc_chemistry) then
          lpenc_requested(i_lambda) = .true.
          lpenc_requested(i_glambda) = .true.
          lpenc_requested(i_lambda1) = .true.
          if (lSmag_heat_transport) lpenc_requested(i_sij2) = .true.
        endif
!
        if (latmchem .or. lcloud) then
          lpenc_requested(i_ppwater) = .true.
          lpenc_requested(i_Ywater) = .true.
        endif
!
        if (ldiffusion .or. lparticles_chemistry) then
          lpenc_requested(i_Diff_penc_add) = .true.
          if (.not. lDiff_fick) then
            lpenc_requested(i_mukmu1) = .true.
            lpenc_requested(i_glnmu) = .true.
          endif
          lpenc_requested(i_gdiffk) = .true.
          lpenc_requested(i_g2XXk) = .true.
          lpenc_requested(i_glnrho) = .true.
          if (lcorr_vel) then
            lpenc_requested(i_grho) = .true.
          endif
        endif
!
      endif
!
!  Needed for condensing species
!
      if (ldustdensity .or. lparticles_radius) then
        lpenc_requested(i_chem_conc) = .true.
        if (lnucleation) then
          lpenc_requested(i_nucl_rate) = .true.
          lpenc_requested(i_nucl_rmin) = .true.
        endif
        if (lnucleation .or. lcondensing_species) then
          lpenc_requested(i_conc_sat_spec) = .true.
        endif
      endif
!
    endsubroutine pencil_criteria_chemistry
!***********************************************************************
    subroutine pencil_interdep_chemistry(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here
!
!  02-03-08/Natalia: coded
!
      logical, dimension(npencils) :: lpencil_in
      !
      if (lpencil_in(i_conc_sat_spec)) then
        lpencil_in(i_TT) = .true.
        lpencil_in(i_TT1) = .true.
      endif
      if (lpencil_in(i_nucl_rate) .or. lpencil_in(i_nucl_rmin))  then
        lpencil_in(i_chem_conc) = .true.
        lpencil_in(i_TT) = .true.
      endif
      if (lpencil_in(i_cv1))    lpencil_in(i_cv) = .true.
      if (lpencil_in(i_cp1))    lpencil_in(i_cp) = .true.
      if (lpencil_in(i_glambda).or. lpencil_in(i_lambda1)) &
          lpencil_in(i_lambda) = .true.
!
      if (lpencil_in(i_H0_RT).or. lpencil_in(i_S0_R))  then
        lpencil_in(i_TT) = .true.
        lpencil_in(i_TT1) = .true.
      endif
      if (lpencil_in(i_S0_R)) then
        lpencil_in(i_lnTT) = .true.
      endif
      if (lpencil_in(i_DYDt_reac))  then
        lpencil_in(i_TT) = .true.
        lpencil_in(i_TT1) = .true.
        lpencil_in(i_lnTT) = .true.
        lpencil_in(i_H0_RT) = .true.
        lpencil_in(i_mu1) = .true.
        lpencil_in(i_rho) = .true.
      endif
      if (lpencil_in(i_hhk_full))  then
        lpencil_in(i_H0_RT) = .true.
        lpencil_in(i_TT) = .true.
      endif
      if (lpencil_in(i_ghhk))  then
        lpencil_in(i_glnTT) = .true.
        lpencil_in(i_TT) = .true.
      endif
      if (lpencil_in(i_ppwater))  then
        lpencil_in(i_TT) = .true.
        lpencil_in(i_rho) = .true.
      endif
      if (lpencil_in(i_chem_conc))  then
        lpencil_in(i_mu1) = .true.
        lpencil_in(i_rho) = .true.
      endif
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_chemistry
!***********************************************************************
    subroutine get_cs2_cheminp(p)
      type(pencil_case) :: p
      !TP: any calls across pencil cannot be faithfully transpiled to the GPU
      if (any(p%cv == 0.0)) then
      else
        p%cs2 = p%cp*p%cv1*p%mu1*p%TT*Rgas
      endif
    endsubroutine get_cs2_cheminp
!***********************************************************************
    subroutine calc_pencils_chemistry(f,p)
!
!  Calculate chemistry pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   13-aug-07/steveb: coded
!   10-jan-11/julien: adapted for the case where chemistry is solved by LSODE
!
      use Sub, only: grad, del2, dot_mn
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(nx,3) :: glncp_tmp
      real, dimension (nx) :: nucleation_rmin, nucleation_rate, conc_sat_spec
!
      intent(inout) :: f
      intent(inout) :: p
      integer :: k,i
      integer, parameter :: ii1=1, ii2=2, ii3=3, ii4=4, ii5=5, ii6=6, ii7=7
      real :: T_low, T_up, T_mid
      real, dimension(nx) :: T_loc, TT_2, TT_3, TT_4
      logical :: ldiffusion2
!
      ldiffusion2 = ldiffusion .and. (.not. lchemonly)
!
      if (lpencil(i_gYYk)) then
        do k = 1,nchemspec
          call grad(f,ichemspec(k),p%gYYk(:,:,k))
        enddo
      endif
!
      if (lpencil(i_gXXk)) then
        do k = 1,nchemspec
          call grad(XX_full(:,:,:,k),p%gXXk(:,:,k))
        enddo
      endif
!
      if (lpencil(i_g2XXk)) then
        do k = 1,nchemspec
          call del2(XX_full(:,:,:,k),p%g2XXk(:,k))
        enddo
      endif
!
      if (lpencil(i_mukmu1)) then
        do k = 1,nchemspec
          p%mukmu1(:,k) = species_constants(k,imass)/unit_mass*p%mu1(:)
        enddo
      endif
!
!  compute chem_conc pencil
!
      if (lpencil(i_chem_conc)) then
        do k = 1,nchemspec
          p%chem_conc(:,k)=XX_full(l1:l2,m,n,k)*p%mu1*p%rho
        enddo
      endif
!
      if (lcheminp) then
!
        if (lpencil(i_glncp) .and. lThCond_simple) then
          call grad(cp_full,glncp_tmp)
          do i = 1,3
            p%glncp(:,i) = glncp_tmp(:,i)/cp_full(l1:l2,m,n)
          enddo
        endif
!
!  Specific heat at constant volume (i.e. density)
!
        if (lpencil(i_cv)) p%cv = cv_full(l1:l2,m,n)
        if (lpencil(i_cv1)) p%cv1 = 1./p%cv
        if (lpencil(i_cp)) p%cp = cp_full(l1:l2,m,n)
        if (lpencil(i_cp1)) p%cp1 = 1./p%cp
        if (lpencil(i_gamma)) p%gamma = p%cp/p%cv
!
        TT_2 = p%TT*p%TT
        TT_3 = TT_2*p%TT
        TT_4 = TT_2*TT_2
!
!  Viscosity of a mixture
!
        if (lpencil(i_nu).and.(.not. lchemonly)) then
          if (visc_const < impossible) then
            p%nu = visc_const
          else
            p%nu = f(l1:l2,m,n,iviscosity)
          endif
          if (lpencil(i_gradnu)) then
            if (visc_const < impossible) then
              p%gradnu = 0.
            else
              call grad(f(:,:,:,iviscosity),p%gradnu)
            endif
          endif
        endif
!
      endif
!
!  Dimensionless Standard-state molar enthalpy H0/RT
!
      if (lpencil(i_H0_RT)) then
        if ((.not. lT_const).and.(ilnTT /= 0)) then
          do k = 1,nchemspec
            T_low = species_constants(k,iTemp1)
            T_mid = species_constants(k,iTemp2)
            T_up = species_constants(k,iTemp3)
!
! Natalia:pencil_check
! if lpencil_check and full compiler settings then
! the problem appears for T_loc= p%TT
! Does anybody know why it is so?
! While this problem is not resolved
! I use T_loc= exp(f(l1:l2,m,n,ilnTT))
!  NILS: This is going to fail if nologtemperature=T. The real error
!  NILS: should be found instead of making a quick fix.
!  NILS: I am not able to reproduce the error natalia reported.
!
!              T_loc= exp(f(l1:l2,m,n,ilnTT))
            T_loc = p%TT
            where (T_loc <= T_mid)
              p%H0_RT(:,k) = species_constants(k,iaa2_offset+ii1) &
                            +species_constants(k,iaa2_offset+ii2)*T_loc/2 &
                            +species_constants(k,iaa2_offset+ii3)*TT_2/3 &
                            +species_constants(k,iaa2_offset+ii4)*TT_3/4 &
                            +species_constants(k,iaa2_offset+ii5)*TT_4/5 &
                            +species_constants(k,iaa2_offset+ii6)/T_loc
            elsewhere
              p%H0_RT(:,k) = species_constants(k,iaa1_offset+ii1) &
                            +species_constants(k,iaa1_offset+ii2)*T_loc/2 &
                            +species_constants(k,iaa1_offset+ii3)*TT_2/3 &
                            +species_constants(k,iaa1_offset+ii4)*TT_3/4 &
                            +species_constants(k,iaa1_offset+ii5)*TT_4/5 &
                            +species_constants(k,iaa1_offset+ii6)/T_loc
            endwhere
          enddo
!
!  Enthalpy flux
!
          if (lpencil(i_hhk_full) ) then
            do k = 1,nchemspec
              if (species_constants(k,imass) > 0.)  then
                p%hhk_full(:,k) = p%H0_RT(:,k)*Rgas*T_loc/species_constants(k,imass)
              endif
            enddo
          endif
!
          if (lpencil(i_ghhk)  .and. (.not. lchemonly)) then
            do k = 1,nchemspec
              if (species_constants(k,imass) > 0.)  then
                do i = 1,3
                  p%ghhk(:,i,k) = (cv_R_spec_full(l1:l2,m,n,k)+1) &
                      /species_constants(k,imass)*Rgas*p%glnTT(:,i)*T_loc(:)
                enddo
              endif
            enddo
          endif
!
        endif
      endif
!
!  Find the entropy by using fifth order temperature fitting function
!
      if (lpencil(i_S0_R) .and. (.not. llsode .or. lchemonly)) then
!AB: Natalia, maybe we should ask earlier for lentropy?
        if ((.not. lT_const).and.(ilnTT /= 0)) then
          do k = 1,nchemspec
            T_low = species_constants(k,iTemp1)
            T_mid = species_constants(k,iTemp2)
            T_up = species_constants(k,iTemp3)
!
! Natalia:pencil_check
! if lpencil_check and full compiler settings then
! the problem appears for T_loc= p%TT
! Does anybody know why it is so?
! While this problem is not resolved
! I use T_loc= exp(f(l1:l2,m,n,ilnTT))
!  NILS: This is going to fail if nologtemperature=T. The real error
!  NILS: should be found instead of making a quick fix.
!  NILS: I am not able to reproduce the error natalia reported.
!
            T_loc = p%TT
!T_loc= exp(f(l1:l2,m,n,ilnTT))
            where (T_loc <= T_mid .and. T_low <= T_loc)
              p%S0_R(:,k) = species_constants(k,iaa2_offset+ii1)*p%lnTT &
                           +species_constants(k,iaa2_offset+ii2)*T_loc &
                           +species_constants(k,iaa2_offset+ii3)*TT_2/2 &
                           +species_constants(k,iaa2_offset+ii4)*TT_3/3 &
                           +species_constants(k,iaa2_offset+ii5)*TT_4/4 &
                           +species_constants(k,iaa2_offset+ii7)
            elsewhere (T_mid <= T_loc .and. T_loc <= T_up)
              p%S0_R(:,k) = species_constants(k,iaa1_offset+ii1)*p%lnTT &
                           +species_constants(k,iaa1_offset+ii2)*T_loc &
                           +species_constants(k,iaa1_offset+ii3)*TT_2/2 &
                           +species_constants(k,iaa1_offset+ii4)*TT_3/3 &
                           +species_constants(k,iaa1_offset+ii5)*TT_4/4 &
                           +species_constants(k,iaa1_offset+ii7)
            endwhere
          enddo
        endif
      endif
!
! Calculate the reaction term and the corresponding pencil
!
      if (lreactions .and. lpencil(i_DYDt_reac)) then
        if (.not. llsode .or. lchemonly) then
          call calc_reaction_term(f,p)
        else
          p%DYDt_reac = 0.
        endif
      else
        p%DYDt_reac = 0.
      endif
!
! Calculate the thermal diffusivity
!
      if (lpencil(i_lambda) .and. lheatc_chemistry) then
        if ((lThCond_simple) .or. (lambda_const < impossible)) then
          if (lThCond_simple) then
!
            p%lambda = lambda_const*p%cp*exp(0.7*log(p%TT(:)/298.))
            if (lpencil(i_glambda))  then
              do i = 1,3
                p%glambda(:,i) = p%lambda(:)*(0.7*p%glnTT(:,i)+p%glncp(:,i))
              enddo
            endif
          else
            p%lambda = lambda_const
            if (lpencil(i_glambda)) p%glambda = 0.
          endif
        elseif (lSmag_heat_transport) then
!
! Natalia
! turbulent heat transport in Smagorinsky case
! probably it should be moved to viscosity module
!
          p%lambda=(0.15*dxmax)**2.*sqrt(2*p%sij2)/Pr_turb*p%cv*p%rho
          if (lpencil(i_glambda)) &
            call not_implemented('calc_pencils_chemistry','glambda pencil for lSmag_heat_transport=T')
        else
          p%lambda = lambda_full(l1:l2,m,n)
          if (lpencil(i_glambda)) call grad(lambda_full,p%glambda)
        endif
        if (lpencil(i_lambda1)) p%lambda1 = 1./max(tini,p%lambda)
      endif
!
! Calculate the diffusion term and the corresponding pencil
!
      if (lcheminp) then
!
! There are 4 cases:
! 1) the case of simplifyed expression for the difusion coef. (Oran paper,)
! 2) the case of constant diffusion coefficients
! 3) the case of constant Lewis numbers with diffusion coef. depending
!    on heat diffusivity
! 4) full complex diffusion coefficients
!
        if (lpencil(i_Diff_penc_add)) then
          if (lDiff_simple) then
            do k = 1,nchemspec
              p%Diff_penc_add(:,k) = Diff_coef_const*p%rho1*exp(0.7*log(p%TT(:)/298.))
              if (lew_exist) p%Diff_penc_add(:,k) = p%Diff_penc_add(:,k)*Lewis_coef1(k)
            enddo
!
!  Constant diffusion coefficients
!
          elseif ((.not. lDiff_simple) .and. (Diff_coef_const < impossible)) then
            if (lSmag_diffusion) then
              if (lcloud) then
                if (z(n)>=z_cloud) then
                  do k = 1,nchemspec
                    p%Diff_penc_add(:,k) = (0.15*dxmax)**2.*sqrt(2*p%sij2)/Sc_number
                  enddo
                else
                  do k = 1,nchemspec
                    p%Diff_penc_add(:,k) = Diff_coef_const*p%rho1
                  enddo
                endif
              else
                do k = 1,nchemspec
                  p%Diff_penc_add(:,k) = (0.15*dxmax)**2.*sqrt(2*p%sij2)/Sc_number
                enddo
              endif
            else
              do k = 1,nchemspec
                p%Diff_penc_add(:,k) = Diff_coef_const*p%rho1
              enddo
            endif
!
!  Diffusion coefficient of a mixture with constant Lewis numbers and
!  given heat conductivity
!
          elseif (lDiff_lewis .and. lew_exist) then
            do k = 1,nchemspec
              p%Diff_penc_add(:,k) = p%lambda*p%rho1*p%cp1*Lewis_coef1(k)
            enddo
          elseif (lDiff_lewis .and. l1step_test) then
            p%Diff_penc_add(:,k) = p%lambda*p%rho1*p%cp1
!
!  Full diffusion coefficient case
!
          else
            do k = 1,nchemspec
              p%Diff_penc_add(:,k) = Diff_full_add(l1:l2,m,n,k)
            enddo
          endif
        endif
        !
        ! Calculate gradient of diffusivity
        !
        if (lpencil(i_gdiffk)) then
          do k=1,nchemspec
            if (lDiff_simple) then
              do i = 1,3
                p%gdiffk(:,i,k) = p%Diff_penc_add(:,k) *(0.7*p%glnTT(:,i)-p%glnrho(:,i))
              enddo
            elseif (lDiff_lewis .and. lew_exist) then
              do i = 1,3
                p%gdiffk(:,i,k) = &
                     (p%glambda(:,i)*p%lambda1(:)-p%glncp(:,i)-p%glnrho(:,i)) &
                     *p%Diff_penc_add(:,k)
              enddo
            elseif (lDiff_lewis .and. l1step_test) then
              do i = 1,3
                p%gdiffk(:,i,k) = &
                     (p%glambda(:,i)*p%lambda1(:)-p%glncp(:,i)-p%glnrho(:,i)) &
                     *p%Diff_penc_add(:,k)
              enddo
            else
              call grad(Diff_full_add(:,:,:,k),p%gdiffk(:,:,k))
            endif
          enddo
        endif
      endif
!
!  More initialization of pencils
!
      if (ldiffusion2 .and. lpencil(i_DYDt_diff)) then
        if (.not. lchemonly) then
          call calc_diffusion_term(f,p)
        else
          p%DYDt_diff = 0.
        endif
      else
        p%DYDt_diff = 0.
      endif
!
      if (latmchem) then
        RHS_Y_full(l1:l2,m,n,:) = p%DYDt_diff
      else
        RHS_Y_full(l1:l2,m,n,:) = p%DYDt_reac+p%DYDt_diff
      endif
!
      if (lpencil(i_cs2) .and. lcheminp) then
        call get_cs2_cheminp(p)
      endif
!
      if (lpencil(i_ppwater) .and. .not. lchemonly) then
        if (index_H2O > 0) p%ppwater = p%rho*Rgas*p%TT/18.*f(l1:l2,m,n,ichemspec(index_H2O))
      endif

      if (lpencil(i_Ywater) .and. .not. lchemonly) then
        if (index_H2O > 0) p%Ywater = f(l1:l2,m,n,ichemspec(index_H2O))
      endif
!
!  Energy per unit mass
!
      if (lpencil(i_ee)) p%ee = p%cv*p%TT
      !
      !  Nucleation rate of condensing species
      !
      if (ldustdensity .or. lparticles) then
        !
        ! Calculate saturation concentration of condensing species
        ! (This must be done before nucleation radius and rate are calculated)
        !
        if (lnucleation .or. lcondensing_species) then
           call cond_spec_sat_conc(p,conc_sat_spec)
           p%conc_sat_spec=conc_sat_spec
        endif
        !
        !  Calculate nucleation rate and corresponding radius
        !
        if (lnucleation) then
           if (lpencil(i_nucl_rate) .or. lpencil(i_nucl_rmin)) then
              call cond_spec_nucl_rate(p,nucleation_rmin,nucleation_rate)
              p%nucl_rate=nucleation_rate
              p%nucl_rmin=nucleation_rmin
            if (.not. lnolatentheat .and. it == 1) then
              ! Initialize some pencil here at the first time-step. This is not
              ! an ideal solution, but will do it like this now to make the
              ! simulations reproducable (avoid *e-314 for the first time step).
              p%latent_heat=0.
              p%ff_cond=0.
            endif
          endif
        endif
      endif
!
    endsubroutine calc_pencils_chemistry
!***********************************************************************
    subroutine flame_front(f)
!
!  05-jun-09/Nils Erland L. Haugen: adapted from similar
!                                   routine in special/chem_stream.f90
!  24-jun-10/Julien Savre: Modifications for lean methane/air combustion
!  This routine set up the initial profiles used in 1D flame speed measurments
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i, j,k
!
      real :: mO2=0., mH2=0., mN2=0., mH2O=0., mCH4=0., mCO2=0., mSIO=0., mSIO2=0.
      real :: log_inlet_density, del, PP
      integer :: i_H2=0, i_O2=0, i_H2O=0, i_N2=0, i_SIO=0, i_SIO2=0
      integer :: ichem_H2=0, ichem_O2=0, ichem_N2=0, ichem_H2O=0, ichem_SIO=0, ichem_SIO2=0
      integer :: i_CH4=0, i_CO2=0, ichem_CH4=0, ichem_CO2=0
      real :: initial_mu1, final_massfrac_O2, final_massfrac_CH4, &
          final_massfrac_H2O, final_massfrac_CO2, &
          final_massfrac_SIO, final_massfrac_SIO2
      real :: init_H2, init_O2, init_N2, init_H2O, init_CO2, init_CH4
      real :: init_SIO, init_SIO2
      logical :: lH2=.false., lO2=.false., lN2=.false., lH2O=.false.
      logical :: lCH4=.false., lCO2=.false.
      logical :: lSIO=.false., lSIO2=.false.
!
      lflame_front = .true.
!
!  Cal air_field to set logicals
!
      call air_field(f,PP)
!
      if (ltemperature_nolog) f(:,:,:,ilnTT) = log(f(:,:,:,ilnTT))
!
! Initialize some indexes
!
      call find_species_index('H2',i_H2,ichem_H2,lH2)
      if (lH2) then
        mH2 = species_constants(ichem_H2,imass)
        init_H2 = initial_massfractions(ichem_H2)
      endif
      call find_species_index('O2',i_O2,ichem_O2,lO2)
      if (lO2) then
        mO2 = species_constants(ichem_O2,imass)
        init_O2 = initial_massfractions(ichem_O2)
      endif
      call find_species_index('N2',i_N2,ichem_N2,lN2)
      if (lN2) then
        mN2 = species_constants(ichem_N2,imass)
        init_N2 = initial_massfractions(ichem_N2)
      else
        init_N2 = 0
      endif
      call find_species_index('SIO',i_SIO,ichem_SIO,lSIO)
      if (lSIO) then
        mSIO = species_constants(ichem_SIO,imass)
        init_SIO = initial_massfractions(ichem_SIO)
      else
        init_SIO = 0
      endif
      call find_species_index('SIO2',i_SIO2,ichem_SIO2,lSIO2)
      if (lSIO2) then
        mSIO2 = species_constants(ichem_SIO2,imass)
        print*,"mSIO2=", mSIO2
        init_SIO2 = initial_massfractions(ichem_SIO2)
      else
        init_SIO2 = 0
      endif
      call find_species_index('H2O',i_H2O,ichem_H2O,lH2O)
      if (lH2O) then
        mH2O = species_constants(ichem_H2O,imass)
        print*,"mH2O=", mH2O
        init_H2O = initial_massfractions(ichem_H2O)
      endif
      call find_species_index('CH4',i_CH4,ichem_CH4,lCH4)
      if (lCH4) then
        mCH4 = species_constants(ichem_CH4,imass)
        init_CH4 = initial_massfractions(ichem_CH4)
      endif
      call find_species_index('CO2',i_CO2,ichem_CO2,lCO2)
      if (lCO2) then
        mCO2 = species_constants(ichem_CO2,imass)
        init_CO2 = initial_massfractions(ichem_CO2)
        final_massfrac_CO2 = init_CO2
      else
        init_CO2 = 0
        final_massfrac_CO2 = init_CO2
      endif
!
! Find approximate value for the mass fraction of O2 after the flame front
! Warning: These formula are only correct for lean fuel/air mixtures. They
!          must be modified under rich conditions to account for the excess
!          of fuel.
!
      final_massfrac_O2 = 0.
      if (lH2 .and. .not. lCH4) then
        final_massfrac_H2O = mH2O/mH2 * init_H2
        final_massfrac_O2 = 1. - final_massfrac_H2O- init_N2
      elseif (lCH4) then
        final_massfrac_CH4 = 0.
        final_massfrac_H2O = 2.*mH2O/mCH4 * init_CH4
        final_massfrac_CO2 = mCO2/mCH4 * init_CH4
        final_massfrac_O2 = 1. - final_massfrac_CO2 - final_massfrac_H2O - init_N2
      elseif (lSIO) then
        final_massfrac_SIO = 0.
        final_massfrac_SIO2 = mSIO2/mSIO * init_SIO
        final_massfrac_O2 = 1. - final_massfrac_SIO2 - init_N2
      endif
!
      if (final_massfrac_O2 < 0.) final_massfrac_O2 = 0.
      if (lroot) then
        print*, '          init                      final'
        if (lH2 .and. .not. lCH4) print*, 'H2 :', init_H2, 0.
        if (lCH4) print*, 'CH4 :',  init_CH4, 0.
        if (lO2) print*, 'O2 :',    init_O2, final_massfrac_O2
        if (lH2O) print*, 'H2O :',  init_H2O, final_massfrac_H2O
        if (lCO2)  print*, 'CO2 :', 0., final_massfrac_CO2
        if (lSIO) print*, 'SiO :',  init_SIO, 0.
        if (lSIO2) print*, 'SiO2 :',  init_SIO2, final_massfrac_SIO2
     endif
!
!  Initialize temperature and species
!
      do k = 1,mx
!
!  Initialize temperature
!
        if (lT_tanh) then
          del = init_x2-init_x1
          f(k,:,:,ilnTT) = f(k,:,:,ilnTT)+log((init_TT2+init_TT1)*0.5  &
              +((init_TT2-init_TT1)*0.5)  &
              *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del)))
        else
          if (x(k) <= init_x1) f(k,:,:,ilnTT) = log(init_TT1)
          if (x(k) >= init_x2) f(k,:,:,ilnTT) = log(init_TT2)
          if (x(k) > init_x1 .and. x(k) < init_x2) &
              f(k,:,:,ilnTT) = log((x(k)-init_x1)/(init_x2-init_x1)*(init_TT2-init_TT1)+init_TT1)
        endif
!
!  Initialize fuel
!
        if (lT_tanh) then
          if (lH2 .and. .not. lCH4) then
            del = (init_x2-init_x1)
            f(k,:,:,i_H2) = (0.+f(l1,:,:,i_H2))*0.5+(0.-f(l1,:,:,i_H2))*0.5  &
                *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
          endif
          if (lCH4) then
            if (k==1) call warning('flame_front','No tanh initial function available for CH4 combustion')
          endif
!
        else
          if (x(k) > init_x1) then
            if (lH2 .and. .not. lCH4) f(k,:,:,i_H2) = init_H2* &
                (exp(f(k,:,:,ilnTT))-init_TT2)/(init_TT1-init_TT2)
            if (lCH4) f(k,:,:,i_CH4) = init_CH4*(exp(f(k,:,:,ilnTT))-init_TT2)/(init_TT1-init_TT2)
            if (lSIO) f(k,:,:,i_SIO) = init_SIO*(exp(f(k,:,:,ilnTT))-init_TT2)/(init_TT1-init_TT2)
          endif
        endif
!
!  Initialize oxygen
!
        if (lT_tanh) then
          del = (init_x2-init_x1)
          f(k,:,:,i_O2) = (f(l2,:,:,i_O2)+f(l1,:,:,i_O2))*0.5  &
              +((f(l2,:,:,i_O2)-f(l1,:,:,i_O2))*0.5)  &
              *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
        else
          if (x(k) > init_x2) f(k,:,:,i_O2) = final_massfrac_O2
          if (x(k) > init_x1 .and. x(k) <= init_x2) &
              f(k,:,:,i_O2) = (x(k)-init_x1)/(init_x2-init_x1)*(final_massfrac_O2-init_O2)+init_O2
        endif
      enddo
!
! Initialize products
!
      if (lT_tanh) then
        do k = 1,mx
          if (lH2 .and. .not. lCH4) then
            del = (init_x2-init_x1)
            f(k,:,:,i_H2O) = (f(l1,:,:,i_H2)/2.*18.+f(l1,:,:,i_H2O))*0.5  &
                +((f(l1,:,:,i_H2)/2.*18.-f(l1,:,:,i_H2O))*0.5)  &
                *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
          endif
          if (lCH4) then
            if (k==1) call warning('flame_front','No tanh initial function available for CH4 combustion')
          endif
        enddo
!
      elseif (.not. l1step_test) then
        if (lH2O) then
           do k = 1,mx
              if (x(k) >= init_x1 .and. x(k) < init_x2) then
                 f(k,:,:,i_H2O) = (x(k)-init_x1)/(init_x2-init_x1) &
                      *final_massfrac_H2O
                 if (lCO2) f(k,:,:,i_CO2) = (x(k)-init_x1)/(init_x2-init_x1) &
                      *final_massfrac_CO2
              elseif (x(k) >= init_x2) then
                 if (lCO2) f(k,:,:,i_CO2) = final_massfrac_CO2
                 if (lH2O) f(k,:,:,i_H2O) = final_massfrac_H2O
              endif
            enddo
          endif
          if (lSIO2) then
           do k = 1,mx
              if (x(k) >= init_x1 .and. x(k) < init_x2) then
                 f(k,:,:,i_SIO2) = (x(k)-init_x1)/(init_x2-init_x1) &
                      *final_massfrac_SIO2
                 if (lCO2) f(k,:,:,i_CO2) = (x(k)-init_x1)/(init_x2-init_x1) &
                      *final_massfrac_CO2
              elseif (x(k) >= init_x2) then
                 if (lCO2) f(k,:,:,i_CO2) = final_massfrac_CO2
                 if (lSIO2) f(k,:,:,i_SIO2) = final_massfrac_SIO2
              endif
           enddo
        end if
      endif
!
   !   if (unit_system == 'cgs') then
   !     Rgas_unit_sys = k_B_cgs/m_u_cgs
   !     Rgas = Rgas_unit_sys/unit_energy
   !  endif
!
!  Find logaritm of density at inlet
!
     initial_mu1 = &
          initial_massfractions(ichem_O2)/(mO2) &
          +initial_massfractions(ichem_N2)/(mN2)
     if (lH2 .and. .not. lCH4) initial_mu1 = initial_mu1+ &
          initial_massfractions(ichem_H2)/(mH2)
     if (lCO2) initial_mu1 = initial_mu1+init_CO2/(mCO2)
     if (lCH4) initial_mu1 = initial_mu1+init_CH4/(mCH4)
     if (lH2O) initial_mu1 = initial_mu1+init_H2O/(mH2O)
     if (lSIO) initial_mu1 = initial_mu1+init_SIO/(mSIO)
     if (lSIO2) initial_mu1 = initial_mu1+init_SIO2/(mSIO2)
     log_inlet_density = log(init_pressure)-log(Rgas)-log(init_TT1)-log(initial_mu1)
!
!  Initialize density
!
      call getmu_array(f,mu1_full)
      f(l1:l2,m1:m2,n1:n2,ilnrho) = log(init_pressure)-log(Rgas)  &
          -f(l1:l2,m1:m2,n1:n2,ilnTT)-log(mu1_full(l1:l2,m1:m2,n1:n2))
      !
!  Initialize velocity
!
!      f(l1:l2,m1:m2,n1:n2,iux)=exp(log_inlet_density - f(l1:l2,m1:m2,n1:n2,ilnrho)) &
!          * (f(l1:l2,m1:m2,n1:n2,iux)+init_ux)
      f(l1:l2,m1:m2,n1:n2,iux) = f(l1:l2,m1:m2,n1:n2,iux)+init_ux
!
!  Check if we want nolog of density or nolog of temperature
!
      if (ldensity_nolog) f(l1:l2,m1:m2,n1:n2,irho) = exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ltemperature_nolog) f(l1:l2,m1:m2,n1:n2,iTT) = exp(f(l1:l2,m1:m2,n1:n2,ilnTT))
!
! Renormalize all species to be sure that the sum of all mass fractions
! are unity
!
      do i = 1,mx
        do j = 1,my
          do k = 1,mz
            f(i,j,k,ichemspec) = f(i,j,k,ichemspec)/sum(f(i,j,k,ichemspec))
          enddo
        enddo
      enddo
!
    endsubroutine flame_front
!***********************************************************************
    subroutine flame_front_new(f)
!
!  10-nov-18/chengeng: adapted from flame_front, but flame front on the left,
!                      and added initial condition for hotspot problem.
!                      This version replaces experimental/test_chemistry.
!
!  This routine set up the initial profiles used in 1D flame speed measurments
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i, j,k
!
      real :: mO2=0., mH2=0., mN2=0., mH2O=0., mCH4=0., mCO2=0.

      real :: log_inlet_density, del, PP
      integer :: i_H2=0, i_O2=0, i_H2O=0, i_N2=0
      integer :: ichem_H2=0, ichem_O2=0, ichem_N2=0, ichem_H2O=0
      integer :: i_CH4=0, i_CO2=0, ichem_CH4=0, ichem_CO2=0
      real :: initial_mu1, final_massfrac_O2, final_massfrac_CH4, &
          final_massfrac_H2O, final_massfrac_CO2,final_massfrac_H2
      real :: init_H2, init_O2, init_N2, init_H2O, init_CO2, init_CH4
      logical :: lH2=.false., lO2=.false., lN2=.false., lH2O=.false.
      logical :: lCH4=.false., lCO2=.false.
      real :: theta
!
      lflame_front = .true.
!
      call air_field(f,PP)
!
      if (ltemperature_nolog) f(:,:,:,ilnTT) = log(f(:,:,:,ilnTT))
!
! Initialize some indexes
!
      call find_species_index('H2',i_H2,ichem_H2,lH2)
      if (lH2) then
        mH2 = species_constants(ichem_H2,imass)
        init_H2 = initial_massfractions(ichem_H2)
      endif
      call find_species_index('O2',i_O2,ichem_O2,lO2)
      if (lO2) then
        mO2 = species_constants(ichem_O2,imass)
        init_O2 = initial_massfractions(ichem_O2)
      endif
      call find_species_index('N2',i_N2,ichem_N2,lN2)
      if (lN2) then
        mN2 = species_constants(ichem_N2,imass)
        init_N2 = initial_massfractions(ichem_N2)
      else
        init_N2 = 0
      endif
      call find_species_index('H2O',i_H2O,ichem_H2O,lH2O)
      if (lH2O) then
        mH2O = species_constants(ichem_H2O,imass)
        init_H2O = initial_massfractions(ichem_H2O)
      endif
!
! Find approximate value for the mass fraction of O2 after the flame front
! Warning: These formula are only correct for lean fuel/air mixtures. They
!          must be modified under rich conditions to account for the excess
!          of fuel.
!
      final_massfrac_O2 = 0.
      final_massfrac_H2 = 0.
      if ( lH2 ) then
        final_massfrac_H2O = mH2O/mH2 * init_H2
        final_massfrac_O2 = 1. - final_massfrac_H2O- init_N2
      endif
!
      if (final_massfrac_O2 < 0.) final_massfrac_O2 = 0.
      if (lroot) then
        print*, '          init                      final'
        if (lH2 ) print*, 'H2 :' , init_H2,  final_massfrac_H2
        if (lO2)  print*, 'O2 :' , init_O2,  final_massfrac_O2
        if (lH2O) print*, 'H2O :', init_H2O, final_massfrac_H2O
      endif
!
!  Initialize temperature and species
!
      do k = 1,mx
!
!  Initialize temperature
!
      if (lhotspot )then
        if( x(k)<init_x2 )then
          theta = ( x(k)-init_x1 )/( init_x2-init_x1 )
          f(k,:,:,ilnTT) = log( init_TT1-theta*(init_TT1-init_TT2) )
        else
          f(k,:,:,ilnTT) = log( init_TT2 )
        end if
        f(k,:,:,i_H2)  = init_H2
        f(k,:,:,i_O2)  = init_O2
        f(k,:,:,i_H2O) = init_H2O
        f(k,:,:,iux) = 0.0
      else
        if ( x(k)<=init_x1 )then
          f(k,:,:,ilnTT) = log( init_TT1 )
          f(k,:,:,i_H2)  = final_massfrac_H2
          f(k,:,:,i_O2)  = final_massfrac_O2
          f(k,:,:,i_H2O) = final_massfrac_H2O
          f(k,:,:,iux) = 0.0
        else if ( x(k)>init_x2 )then
          f(k,:,:,ilnTT) = log( init_TT2 )
          f(k,:,:,i_H2)  = init_H2
          f(k,:,:,i_O2)  = init_O2
          f(k,:,:,i_H2O) = init_H2O
          f(k,:,:,iux) = init_ux*(init_TT1/init_TT2-1.0)
        else
          theta = ( x(k)-init_x1 )/( init_x2-init_x1 )
          f(k,:,:,ilnTT) = log( init_TT1-theta*(init_TT1-init_TT2) )
          f(k,:,:,i_H2)  = init_H2*theta+final_massfrac_H2
          f(k,:,:,i_O2)  = init_O2*theta+final_massfrac_O2
          f(k,:,:,i_H2O) = init_H2O*(1.0-theta)+final_massfrac_H2O
          f(k,:,:,iux) = init_ux*(init_TT1/init_TT2-1.0)*theta
        end if
      end if

      end do
!
     ! if (unit_system == 'cgs') then
     !   Rgas_unit_sys = k_B_cgs/m_u_cgs
     !   Rgas = Rgas_unit_sys/unit_energy
     ! endif
!
!  Find logaritm of density at inlet
!
      initial_mu1 = &
          initial_massfractions(ichem_O2)/(mO2) &
          +initial_massfractions(ichem_H2O)/(mH2O) &
          +initial_massfractions(ichem_N2)/(mN2)
      if (lH2 .and. .not. lCH4) initial_mu1 = initial_mu1+initial_massfractions(ichem_H2)/(mH2)
      if (lCO2) initial_mu1 = initial_mu1+init_CO2/(mCO2)
      if (lCH4) initial_mu1 = initial_mu1+init_CH4/(mCH4)
      log_inlet_density = log(init_pressure)-log(Rgas)-log(init_TT1)-log(initial_mu1)
!
!  Initialize density
!
      call getmu_array(f,mu1_full)
      f(l1:l2,m1:m2,n1:n2,ilnrho) = log(init_pressure)-log(Rgas)  &
          -f(l1:l2,m1:m2,n1:n2,ilnTT)-log(mu1_full(l1:l2,m1:m2,n1:n2))
!
!  Initialize velocity
!
!      f(l1:l2,m1:m2,n1:n2,iux)=exp(log_inlet_density - f(l1:l2,m1:m2,n1:n2,ilnrho)) &
!          * (f(l1:l2,m1:m2,n1:n2,iux)+init_ux)
      !f(l1:l2,m1:m2,n1:n2,iux) = f(l1:l2,m1:m2,n1:n2,iux)+init_ux
!
! Renormalize all species to be sure that the sum of all mass fractions
! are unity
!
      do i = 1,mx
        do j = 1,my
          do k = 1,mz
            f(i,j,k,ichemspec) = f(i,j,k,ichemspec)/sum(f(i,j,k,ichemspec))
          enddo
        enddo
      enddo
!
    endsubroutine flame_front_new
!***********************************************************************
    subroutine TTD(f)
!
!  15-may-03/Nils Erland L. Haugen: adapted from flame_front
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i, j,k
!
      real :: initial_mu1, ksi_TTD, dTdr_c, deltaT, PP
!
      call air_field(f,PP)
!
      ksi_TTD = 1.
      dTdr_c = 2000 ![K/m]
!      dTdr_c=20000 ![K/m]
!
!  Initialize temperature
!
      do k = l1,l2
        if (ltemperature_nolog) then
          f(k,:,:,iTT) = f(k,:,:,iTT)+ksi_TTD*dTdr_c*(xyz1(1)-x(k))/100.
        else
          deltaT = ksi_TTD*dTdr_c*(xyz1(1)-x(k))/100.
          f(k,:,:,ilnTT) = log(exp(f(k,:,:,ilnTT))+deltaT)
!          print*,'deltaT=',deltaT, exp(f(k,m1,n1,ilnTT)), f(k,m1,n1,ilnTT)
        endif
      enddo
!
!  Initialize density
!
      call getmu_array(f,mu1_full)
      f(l1:l2,m1:m2,n1:n2,ilnrho) = log(PP)-log(Rgas)  &
          -f(l1:l2,m1:m2,n1:n2,ilnTT)-log(mu1_full(l1:l2,m1:m2,n1:n2))
!
!  Initialize velocity
!
      f(l1:l2,m1:m2,n1:n2,iux) = f(l1:l2,m1:m2,n1:n2,iux)+init_ux
!
!  Check if we want nolog of density or nolog of temperature
!
      if (ldensity_nolog) f(l1:l2,m1:m2,n1:n2,irho) = exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ltemperature_nolog) f(l1:l2,m1:m2,n1:n2,iTT) = exp(f(l1:l2,m1:m2,n1:n2,ilnTT))
!
! Renormalize all species to be sure that the sum of all mass fractions
! are unity
!
      do i = 1,mx
        do j = 1,my
          do k = 1,mz
            f(i,j,k,ichemspec) = f(i,j,k,ichemspec)/sum(f(i,j,k,ichemspec))
          enddo
        enddo
      enddo
!
    endsubroutine TTD
!***********************************************************************
    subroutine triple_flame(f)
!
! 26-jul-10/Julien Savre: Copy from the flame_front case
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i, j,k
!
      real :: mO2=0., mH2=0., mN2=0., mH2O=0., mCH4=0., mCO2=0.
      real :: del, PP
      integer :: i_H2=0, i_O2=0, i_H2O=0, i_N2=0
      integer :: ichem_H2=0, ichem_O2=0, ichem_N2=0, ichem_H2O=0
      integer :: i_CH4=0, i_CO2=0, ichem_CH4=0, ichem_CO2=0
      real :: final_massfrac_O2, final_massfrac_CH4, final_massfrac_H2O, final_massfrac_CO2
      real :: init_H2, init_O2, init_N2, init_H2O, init_CO2, init_CH4
      real :: beta
      real :: init_y1, init_y2
      real :: init_rho, init_m1
      real, dimension(ny) :: dim
      logical :: lH2=.false., lO2=.false., lN2=.false., lH2O=.false.
      logical :: lCH4=.false., lCO2=.false.
!
      ltriple_flame = .true.
!
      call air_field(f,PP)
!
      init_y1 = xyz0(2) + Lxyz(2)/3.
      init_y2 = xyz0(2) + 2.*Lxyz(2)/3.
!
      if (ltemperature_nolog) f(:,:,:,ilnTT) = log(f(:,:,:,ilnTT))
      if (lroot) print*, 'init_chem: triple_flame '
!
! Initialize some indexes
!
      call find_species_index('H2',i_H2,ichem_H2,lH2)
      if (lH2) then
        mH2 = species_constants(ichem_H2,imass)
        init_H2 = initial_massfractions(ichem_H2)
      endif
      call find_species_index('O2',i_O2,ichem_O2,lO2)
      if (lO2) then
        mO2 = species_constants(ichem_O2,imass)
        init_O2 = initial_massfractions(ichem_O2)
      endif
      call find_species_index('N2',i_N2,ichem_N2,lN2)
      if (lN2) then
        mN2 = species_constants(ichem_N2,imass)
        init_N2 = initial_massfractions(ichem_N2)
      endif
      call find_species_index('H2O',i_H2O,ichem_H2O,lH2O)
      if (lH2O) then
        mH2O = species_constants(ichem_H2O,imass)
        init_H2O = initial_massfractions(ichem_H2O)
      endif
      call find_species_index('CH4',i_CH4,ichem_CH4,lCH4)
      if (lCH4) then
        mCH4 = species_constants(ichem_CH4,imass)
        init_CH4 = initial_massfractions(ichem_CH4)
      endif
      call find_species_index('CO2',i_CO2,ichem_CO2,lCO2)
      if (lCO2) then
        mCO2 = species_constants(ichem_CO2,imass)
        init_CO2 = initial_massfractions(ichem_CO2)
      endif
!
! Find approximate value for the mass fraction of O2 after the flame front
!
      final_massfrac_O2 = 0.
      if (lH2) then
        final_massfrac_H2O = mH2O/(2.*mH2) * init_H2
        final_massfrac_O2 = 1. - final_massfrac_H2O - init_N2
      elseif (lCH4) then
        final_massfrac_CH4 = 0.
        final_massfrac_H2O = 2.*mH2O/mCH4 * init_CH4
        final_massfrac_CO2 = mCO2/mCH4 * init_CH4
        final_massfrac_O2 = 1. - final_massfrac_CO2 - final_massfrac_H2O - init_N2
      endif
!
      if (final_massfrac_O2 < 0.) final_massfrac_O2 = 0.
      if (lroot) then
        print*, '          init                      final'
        if (lH2) print*, 'H2 :', init_H2, 0.
        if (lCH4) print*, 'CH4 :', init_CH4, 0.
        if (lO2) print*, 'O2 :', init_O2, final_massfrac_O2
        if (lH2O) print*, 'H2O :', 0., final_massfrac_H2O
        if (lCO2)  print*, 'CO2 :', 0., final_massfrac_CO2
      endif
!
!  Initialize temperature and species
!
      if (lT_tanh) then
        if (lroot) print*, 'Temperature initialization: tanh function.'
      else
        if (lroot) print*, 'Temperature initialization: linear.'
      endif
!
      do k = 1,mx
!
!  Initialize temperature
!
        if (lT_tanh) then
          del = init_x2-init_x1
          f(k,:,:,ilnTT) = f(k,:,:,ilnTT)+log((init_TT2+init_TT1)*0.5  &
              +((init_TT2-init_TT1)*0.5)  &
              *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del)))
        else
          if (x(k) <= init_x1) f(k,:,:,ilnTT) = log(init_TT1)
          if (x(k) >= init_x2) f(k,:,:,ilnTT) = log(init_TT2)
          if (x(k) > init_x1 .and. x(k) < init_x2) &
            f(k,:,:,ilnTT) = log((x(k)-init_x1)/(init_x2-init_x1) &
                *(init_TT2-init_TT1)+init_TT1)
        endif
!
!  Initialize fuel
!
        if (lT_tanh) then
          if (lH2) then
            del = (init_x2-init_x1)
            f(k,:,:,i_H2) = (0.+f(l1,:,:,i_H2))*0.5  &
                +(0.-f(l1,:,:,i_H2))*0.5  &
                *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
          endif
          if (lCH4) then
            if (k==1) call warning('flame_front','No tanh initial function available for CH4 combustion')
          endif
!
        else
          if (x(k) > init_x1) then
            if (lH2) &
              f(k,:,:,i_H2) = init_H2*(exp(f(k,:,:,ilnTT))-init_TT2)/(init_TT1-init_TT2)
            if (lCH4) &
              f(k,:,:,i_CH4) = init_CH4*(exp(f(k,:,:,ilnTT))-init_TT2)/(init_TT1-init_TT2)
          endif
        endif
!
!  Initialize oxygen
!
        if (lT_tanh) then
          del = (init_x2-init_x1)
          f(k,:,:,i_O2) = (f(l2,:,:,i_O2)+f(l1,:,:,i_O2))*0.5  &
              +((f(l2,:,:,i_O2)-f(l1,:,:,i_O2))*0.5)  &
              *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
        else
!
          if (x(k) > init_x2) f(k,:,:,i_O2) = final_massfrac_O2
          if (x(k) > init_x1 .and. x(k) <= init_x2) &
              f(k,:,:,i_O2) = (x(k)-init_x1)/(init_x2-init_x1)*(final_massfrac_O2-init_O2)+init_O2
        endif
      enddo
!
! Initialize products
!
      if (lT_tanh) then
        do k = 1,mx
          if (lH2) then
            del = (init_x2-init_x1)
            f(k,:,:,i_H2O) = (f(l1,:,:,i_H2)/2.*18.+f(l1,:,:,i_H2O))*0.5  &
                +((f(l1,:,:,i_H2)/2.*18.-f(l1,:,:,i_H2O))*0.5)  &
                *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
          endif
          if (lCH4) then
            if (k==1) call warning('flame_front','no tanh initial function available for CH4 combustion')
          endif
        enddo
      else
        do k = 1,mx
          if (x(k) >= init_x1 .and. x(k) < init_x2) then
            f(k,:,:,i_H2O) = (x(k)-init_x1)/(init_x2-init_x1)*final_massfrac_H2O
            if (lCO2) f(k,:,:,i_CO2) = (x(k)-init_x1)/(init_x2-init_x1)*final_massfrac_CO2
          elseif (x(k) >= init_x2) then
            if (lCO2) f(k,:,:,i_CO2) = final_massfrac_CO2
            if (lH2O) f(k,:,:,i_H2O) = final_massfrac_H2O
          endif
        enddo
      endif
!
    !  if (unit_system == 'cgs') then
    !    Rgas_unit_sys = k_B_cgs/m_u_cgs
    !    Rgas = Rgas_unit_sys/unit_energy
    !  endif
!
!  Set the initial equivalence ratio gradient in the fresh gases
!
      dim(1:ny) = y(m1:m2)
      beta = init_N2/init_O2
      do i = 1, mx
        if (x(i) <= init_x1) then
          do j = 1, ny
            if (dim(j) >= init_y1 .and. dim(j) <= init_y2) then
              if (lH2) then
                f(i,m1+j-1,:,i_H2) = init_zz2 - (init_zz2-init_zz1) * (init_y2 - dim(j)) / &
                    (init_y2 - init_y1)
              elseif (lCH4) then
                f(i,m1+j-1,:,i_CH4) = init_zz2 - (init_zz2-init_zz1) * (init_y2 - dim(j)) / &
                    (init_y2 - init_y1)
              endif
            elseif (dim(j) <= init_y1) then
              if (lH2) then
                f(i,m1+j-1,:,i_H2) = init_zz1
              elseif (lCH4) then
                f(i,m1+j-1,:,i_CH4) = init_zz1
              endif
            elseif (dim(j) >= init_y2) then
              if (lH2) then
                f(i,m1+j-1,:,i_H2) = init_zz2
              elseif (lCH4) then
                f(i,m1+j-1,:,i_CH4) = init_zz2
              endif
            endif
          enddo
!
          if (lH2) then
            f(i,ny:my,:,i_H2) = init_zz2
            if (lO2) f(i,:,:,i_O2) = (1.-f(i,:,:,i_H2)) / (1.+beta)
            if (lN2) f(i,:,:,i_N2) = 1.-f(i,:,:,i_O2)-f(i,:,:,i_H2)
          elseif (lCH4) then
            f(i,ny:my,:,i_CH4) = init_zz2
            if (lO2) f(i,:,:,i_O2) = (1.-f(i,:,:,i_CH4)) / (1.+beta)
            if (lN2) f(i,:,:,i_N2) = 1.-f(i,:,:,i_O2)-f(i,:,:,i_CH4)
            if (lCO2) f(i,:,:,i_CO2) = 0.
            if (lH2O) f(i,:,:,i_H2O) = 0.
          endif
        endif
      enddo
!
! Renormalize all species to be sure that the sum of all mass fractions
! are unity
!
      do i = 1,mx
        do j = 1,my
          do k = 1,mz
            f(i,j,k,ichemspec) = f(i,j,k,ichemspec)/sum(f(i,j,k,ichemspec))
          enddo
        enddo
      enddo
!
!  Initialize density
!
      call getmu_array(f,mu1_full)
      f(l1:l2,m1:m2,n1:n2,ilnrho) = log(init_pressure)-log(Rgas)  &
          -f(l1:l2,m1:m2,n1:n2,ilnTT)-log(mu1_full(l1:l2,m1:m2,n1:n2))
!
!  Initialize velocity
!
      if (lCH4) then
        init_m1 = init_CH4/mCH4 + init_O2/mO2 + init_N2/mN2
        init_rho = init_pressure/(init_TT1 * init_m1 * Rgas)
        f(l1:l2,m1:m2,n1:n2,iux) = f(l1:l2,m1:m2,n1:n2,iux) +   &
            init_ux * init_rho / exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
      else
        f(l1:l2,m1:m2,n1:n2,iux) = f(l1:l2,m1:m2,n1:n2,iux)+init_ux
      endif
!
!  Check if we want nolog of density or nolog of temperature
!
      if (ldensity_nolog) f(l1:l2,m1:m2,n1:n2,irho) = exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ltemperature_nolog) f(l1:l2,m1:m2,n1:n2,iTT) = exp(f(l1:l2,m1:m2,n1:n2,ilnTT))
!
    endsubroutine triple_flame
!***********************************************************************
    subroutine flame(f)
!
! 05-jun-09/Nils Erland L. Haugen: adapted from similar
!                                   routine in special/chem_stream.f90
! This routine set up the initial profiles used in 1D flame speed measurments
! NILS: This routine is essentially the samw as flame_front, but I leave
! NILS: flame_front as it is for now for backwards compatibility.
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i, j,k
!
      real :: mO2, mH2, mN2, mH2O, mCH4, mCO2
      real :: log_inlet_density, del, PP
      integer :: i_H2, i_O2, i_H2O, i_N2, ichem_H2, ichem_O2, ichem_N2, ichem_H2O
      integer :: i_CH4, i_CO2, ichem_CH4, ichem_CO2
      real :: initial_mu1, final_massfrac_O2
      real :: init_H2, init_O2, init_N2, init_H2O, init_CO2, init_CH4
      logical :: lH2, lO2, lN2, lH2O, lCH4, lCO2
!
      lflame_front = .true.
!
      call air_field(f,PP)
!
      if (ltemperature_nolog) f(:,:,:,ilnTT) = log(f(:,:,:,ilnTT))
!
! Initialize some indexes
!
      call find_species_index('H2',i_H2,ichem_H2,lH2)
      if (lH2) then
        mH2 = species_constants(ichem_H2,imass)
        init_H2 = initial_massfractions(ichem_H2)
      endif
      call find_species_index('O2',i_O2,ichem_O2,lO2)
      if (lO2) then
        mO2 = species_constants(ichem_O2,imass)
        init_O2 = initial_massfractions(ichem_O2)
      endif
      call find_species_index('N2',i_N2,ichem_N2,lN2)
      if (lN2) then
        mN2 = species_constants(ichem_N2,imass)
        init_N2 = initial_massfractions(ichem_N2)
      endif
      call find_species_index('H2O',i_H2O,ichem_H2O,lH2O)
      if (lH2O) then
        mH2O = species_constants(ichem_H2O,imass)
        init_H2O = initial_massfractions(ichem_H2O)
      endif
      call find_species_index('CH4',i_CH4,ichem_CH4,lCH4)
      if (lCH4) then
        mCH4 = species_constants(ichem_CH4,imass)
        init_CH4 = initial_massfractions(ichem_CH4)
      endif
      call find_species_index('CO2',i_CO2,ichem_CO2,lCO2)
      if (lCO2) then
        mCO2 = species_constants(ichem_CO2,imass)
        init_CO2 = initial_massfractions(ichem_CO2)
      endif
!
! Find approximate value for the mass fraction of O2 after the flame front
!
      final_massfrac_O2 = (init_O2/mO2 -init_H2/(2.*mH2) -init_CH4*2/(mCH4))*mO2
!
      if (final_massfrac_O2 < 0.) final_massfrac_O2 = 0.
!
!  Initialize temperature and species
!
      do k = 1,mx
!
!  Initialize temperature
!
        if (lT_tanh) then
          del = init_x2-init_x1
          f(k,:,:,ilnTT) = f(k,:,:,ilnTT)+log((init_TT2+init_TT1)*0.5  &
              +((init_TT2-init_TT1)*0.5)  &
              *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del)))
        else
          if (x(k) <= init_x1) f(k,:,:,ilnTT) = log(init_TT1)
          if (x(k) >= init_x2) f(k,:,:,ilnTT) = log(init_TT2)
          if (x(k) > init_x1 .and. x(k) < init_x2) &
            f(k,:,:,ilnTT) = log((x(k)-init_x1)/(init_x2-init_x1)*(init_TT2-init_TT1)+init_TT1)
        endif
!
!  Initialize steam and hydrogen
!
        if (lT_tanh) then
          del = (init_x2-init_x1)
          if (lH2) then
            f(k,:,:,i_H2) = (0.+f(l1,:,:,i_H2))*0.5  &
                +(0.-f(l1,:,:,i_H2))*0.5  &
                *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
            f(k,:,:,i_H2O) = (f(l1,:,:,i_H2)/2.*18.+f(l1,:,:,i_H2O))*0.5  &
                +((f(l1,:,:,i_H2)/2.*18.-f(l1,:,:,i_H2O))*0.5)  &
                *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
          endif
          if (lCH4) then
            f(k,:,:,i_CH4) = init_CH4*0.5-init_CH4*0.5  &
                *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
            f(k,:,:,i_H2O) = init_H2O+(init_CH4-f(k,:,:,i_CH4))*2.*mH2O/mCH4
            f(k,:,:,i_CO2) = init_CO2+(init_CH4-f(k,:,:,i_CH4))*1.*mCO2/mCH4
            f(k,:,:,i_O2) = init_O2 -(init_CH4-f(k,:,:,i_CH4))*2.*mO2 /mCH4
          endif
        else
          if (x(k) > init_x1) then
            if (lH2) &
              f(k,:,:,i_H2) = init_H2*(exp(f(k,:,:,ilnTT))-init_TT2)/(init_TT1-init_TT2)
            if (lCH4) &
              f(k,:,:,i_CH4) = init_CH4*(exp(f(k,:,:,ilnTT))-init_TT2)/(init_TT1-init_TT2)
          endif
        endif
!
!  Initialize oxygen
!
        if (lT_tanh) then
! $          del=(init_x2-init_x1)
! $          f(k,:,:,i_O2)=(f(l2,:,:,i_O2)+f(l1,:,:,i_O2))*0.5  &
! $              +((f(l2,:,:,i_O2)-f(l1,:,:,i_O2))*0.5)  &
! $              *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
        else
          if (x(k) > init_x2) then
            f(k,:,:,i_O2) = final_massfrac_O2
          endif
          if (x(k) > init_x1 .and. x(k) < init_x2) then
            f(k,:,:,i_O2) = (x(k)-init_x2)/(init_x1-init_x2) &
                *(init_O2-final_massfrac_O2)+final_massfrac_O2
          endif
        endif
      enddo
!
! Initialize steam and CO2
!
      if (.not. lT_tanh) then
        do k = 1,mx
          if (x(k) >= init_x1) then
            if (final_massfrac_O2 > 0.) then
              f(k,:,:,i_H2O) = initial_massfractions(ichem_H2)/mH2*mH2O &
                  *(exp(f(k,:,:,ilnTT))-init_TT1)/(init_TT2-init_TT1)+init_H2O
            else
              if (x(k) >= init_x2) then
                if (lCO2) f(k,:,:,i_CO2) = init_CO2+(init_CH4-f(k,:,:,i_CH4))
                f(k,:,:,i_H2O) = 1.-f(k,:,:,i_N2)-f(k,:,:,i_H2)-f(k,:,:,i_O2)
                if (lCH4) f(k,:,:,i_H2O) = f(k,:,:,i_H2O)-f(k,:,:,i_CH4)
                if (lCO2) f(k,:,:,i_H2O) = f(k,:,:,i_H2O)-f(k,:,:,i_CO2)
              else
                f(k,:,:,i_H2O) = (x(k)-init_x1)/(init_x2-init_x1) &
                    *((1.-f(l2,:,:,i_N2)-f(l2,:,:,i_H2)) &
                    -initial_massfractions(ichem_H2O)) &
                    +initial_massfractions(ichem_H2O)
              endif
            endif
          endif
        enddo
      endif
!
    !  if (unit_system == 'cgs') then
    !    Rgas_unit_sys = k_B_cgs/m_u_cgs
    !    Rgas = Rgas_unit_sys/unit_energy
    !  endif
!
!  Find logaritm of density at inlet
!
      initial_mu1 = initial_massfractions(ichem_H2)/(mH2) &
                   +initial_massfractions(ichem_O2)/(mO2) &
                   +initial_massfractions(ichem_H2O)/(mH2O) &
                   +initial_massfractions(ichem_N2)/(mN2)
      if (lCO2) initial_mu1 = initial_mu1+init_CO2/(mCO2)
      if (lCH4) initial_mu1 = initial_mu1+init_CH4/(mCH4)
      log_inlet_density = log(init_pressure)-log(Rgas)-log(init_TT1)-log(initial_mu1)
!
      call getmu_array(f,mu1_full)
!
!  Initialize density
!
      f(l1:l2,m1:m2,n1:n2,ilnrho) = log(init_pressure)-log(Rgas)  &
          -f(l1:l2,m1:m2,n1:n2,ilnTT)-log(mu1_full(l1:l2,m1:m2,n1:n2))
!
!  Initialize velocity
!
      f(l1:l2,m1:m2,n1:n2,iux) = f(l1:l2,m1:m2,n1:n2,iux)+init_ux
!
!  Check if we want nolog of density or nolog of temperature
!
      if (ldensity_nolog) f(l1:l2,m1:m2,n1:n2,irho) = exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ltemperature_nolog) f(l1:l2,m1:m2,n1:n2,iTT) = exp(f(l1:l2,m1:m2,n1:n2,ilnTT))
!
! Renormalize all species too be sure that the sum of all mass fractions
! are unity
!
      do i = 1,mx
        do j = 1,my
          do k = 1,mz
            f(i,j,k,ichemspec) = f(i,j,k,ichemspec)/sum(f(i,j,k,ichemspec))
          enddo
        enddo
      enddo
!
    endsubroutine flame
!***********************************************************************
    subroutine flame_blob(f)
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: j1, j2, j3
!
      real :: mO2, mH2, mN2, mH2O, PP
      integer :: i_H2, i_O2, i_H2O, i_N2, ichem_H2, ichem_O2, ichem_N2, ichem_H2O
      real :: initial_mu1, final_massfrac_O2
      logical :: found_specie
!
      real :: Rad, sz1, sz2
!
      lflame_front = .true.
!
      call air_field(f,PP)
!
! Initialize some indexes
!
      call find_species_index('H2',i_H2,ichem_H2,found_specie)
      call find_species_index('O2',i_O2,ichem_O2,found_specie)
      call find_species_index('N2',i_N2,ichem_N2,found_specie)
      call find_species_index('H2O',i_H2O,ichem_H2O,found_specie)
      mO2 = species_constants(ichem_O2,imass)
      mH2 = species_constants(ichem_H2,imass)
      mH2O = species_constants(ichem_H2O,imass)
      mN2 = species_constants(ichem_N2,imass)
!
! Find approximate value for the mass fraction of O2 after the flame front
!
      final_massfrac_O2 = (initial_massfractions(ichem_O2)/mO2 &
                          -initial_massfractions(ichem_H2)/(2*mH2))*mO2
!
!  Initialize temperature and species in air_field(f)
!
    !  if (unit_system == 'cgs') then
    !    Rgas_unit_sys = k_B_cgs/m_u_cgs
    !    Rgas = Rgas_unit_sys/unit_energy
    !  endif
!
!  Find logaritm of density at inlet
!
      initial_mu1 = initial_massfractions(ichem_H2)/(mH2) &
                   +initial_massfractions(ichem_O2)/(mO2) &
                   +initial_massfractions(ichem_H2O)/(mH2O) &
                   +initial_massfractions(ichem_N2)/(mN2)
!
      call getmu_array(f,mu1_full)
!
      do j3 = 1,mz
        do j2 = 1,my
          do j1 = 1,mx
!
            Rad = 0.
            if (nxgrid > 1) Rad =     x(j1)*x(j1)
            if (nygrid > 1) Rad = Rad+y(j2)*y(j2)
            if (nzgrid > 1) Rad = Rad+z(j3)*z(j3)
!
            Rad = sqrt(Rad)
!
            f(j1,j2,j3,ilnTT) = log((init_TT2-init_TT1)*exp(-(Rad/init_x2)**2)+init_TT1)
            mu1_full(j1,j2,j3) = f(j1,j2,j3,i_H2)/(mH2)+f(j1,j2,j3,i_O2)/(mO2) &
                                +f(j1,j2,j3,i_H2O)/(mH2O)+f(j1,j2,j3,i_N2)/(mN2)
!
            f(j1,j2,j3,ilnrho) = log(init_pressure)-log(Rgas)-f(j1,j2,j3,ilnTT)  &
                                -log(mu1_full(j1,j2,j3))
!
!  Initialize velocity
!
            f(j1,j2,j3,iux) = f(j1,j2,j3,iux)  &
                             +init_ux!*exp(log_inlet_density)/exp(f(j1,j2,j3,ilnrho))
            f(j1,j2,j3,iuy) = f(j1,j2,j3,iuy)+ init_uy
            f(j1,j2,j3,iuz) = f(j1,j2,j3,iuz)+ init_uz
!
            if (nxgrid == 1) f(j1,j2,j3,iux) = 0.
            if (nygrid == 1) f(j1,j2,j3,iuy) = 0.
            if (nzgrid == 1) f(j1,j2,j3,iuz) = 0.
!
            if (nxgrid /= 1) then
              sz1 = (xyz0(1)+Lxyz(1)*0.15)
              sz2 = (xyz0(1)+Lxyz(1)*(1.-0.15))
            endif
!
            if (nygrid /= 1)   then
              sz1 = (xyz0(2)+Lxyz(2)*0.15)
              sz2 = (xyz0(2)+Lxyz(2)*(1.-0.15))
            endif
!
            if (nzgrid /= 1)  then
              sz1 = (xyz0(3)+Lxyz(3)*0.15)
              sz2 = (xyz0(3)+Lxyz(3)*(1.-0.15))
            endif
!
          enddo
        enddo
      enddo
!
!  Check if we want nolog of density
!
      if (ldensity_nolog) f(:,:,:,irho) = exp(f(:,:,:,ilnrho))
!
    endsubroutine flame_blob
!***********************************************************************
    subroutine opposite_flames(f)
!
!  03-jan-10/nilshau: adapted from opposite_ignitions
!
!  Set up two oppositely directed flame fronts in the x-direction.
!  The two fronts have fresh gas between them.
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: j1, j2, j3
!
      real :: mO2, mH2, mN2, mH2O, lower, upper, PP
      integer :: i_H2, i_O2, i_H2O, i_N2, ichem_H2, ichem_O2, ichem_N2, ichem_H2O
      integer :: i_C3H8, ichem_C3H8, i_CO2, ichem_CO2
      real :: final_massfrac_O2, mu1, phi, delta_O2, mC3H8, mCO2
      logical :: found_specie, lH2, lCO2, lC3H8, lH2O
      real :: norm, flat_range
!
      lflame_front = .true.
!
      call air_field(f,PP)
!
! Initialize some indexes
!
      call find_species_index('H2',i_H2,ichem_H2,lH2)
      call find_species_index('O2',i_O2,ichem_O2,found_specie)
      call find_species_index('N2',i_N2,ichem_N2,found_specie)
      call find_species_index('H2O',i_H2O,ichem_H2O,lH2O)
      call find_species_index('C3H8',i_C3H8,ichem_C3H8,lC3H8)
      call find_species_index('CO2',i_CO2,ichem_CO2,lCO2)
      mO2 = species_constants(ichem_O2,imass)
      mH2O = species_constants(ichem_H2O,imass)
      mN2 = species_constants(ichem_N2,imass)
      if (lC3H8) mC3H8 = species_constants(ichem_C3H8,imass)
      if (lH2)   mH2   = species_constants(ichem_H2,imass)
      if (lCO2)  mCO2 = species_constants(ichem_CO2,imass)
!
! Find approximate value for the mass fraction of O2 after the flame front
!
      final_massfrac_O2 = initial_massfractions(ichem_O2)
      if (lH2) then
        delta_O2 = initial_massfractions(ichem_H2)/(2*mH2)*mO2
      else
        delta_O2 = 0.
      endif
      final_massfrac_O2 = final_massfrac_O2-delta_O2
!
      if (lC3H8) then
        delta_O2 = 5*initial_massfractions(ichem_C3H8)/mC3H8*mO2
      else
        delta_O2 = 0.
      endif
      final_massfrac_O2 = final_massfrac_O2-delta_O2
!
      if (ltemperature_nolog) call not_implemented('opposite_flames','for ltemperature_nolog=T')
!
!  Loop over all grid points
!
      do j3 = 1,mz
        do j2 = 1,my
          do j1 = 1,mx
!
!  First define the distance from the lower and upper domain boundary.
!
            lower = x(j1)-xyz0(1)
            upper = xyz1(1)-x(j1)
!
!  Find progress variable phi based on distance from boundaries.
!
            flat_range = init_x2*0.1
            phi = exp(-((lower-flat_range)/init_x2)**2) &
                 +exp(-((upper-flat_range)/init_x2)**2)
            if (phi > 1.0) phi = 1.
            if (flat_range > lower) phi = 1.
            if (flat_range > upper) phi = 1.
!
!  Find temperature, species and density based on progress variable
!
            f(j1,j2,j3,ilnTT) = log(init_TT1+phi*(init_TT2-init_TT1))
            if (lH2) then
              f(j1,j2,j3,i_H2) = (1-phi)*initial_massfractions(ichem_H2)
              f(j1,j2,j3,i_H2O) = phi*initial_massfractions(ichem_H2)*mH2O/mH2
            endif
            if (lC3H8) then
              f(j1,j2,j3,i_C3H8) = (1-phi)*initial_massfractions(ichem_C3H8)
              f(j1,j2,j3,i_H2O) = 4*phi*initial_massfractions(ichem_C3H8)*mH2O/mC3H8
              f(j1,j2,j3,i_CO2) = 3*phi*initial_massfractions(ichem_C3H8)*mCO2/mC3H8
            endif
            f(j1,j2,j3,i_O2) = (1-phi)*(initial_massfractions(ichem_O2) &
                -final_massfrac_O2)+final_massfrac_O2
!
!  Re-normalize mass fractions
!
            norm = f(j1,j2,j3,i_O2)+f(j1,j2,j3,i_N2)
            if (lC3H8) norm = norm + f(j1,j2,j3,i_C3H8)
            if (lCO2)  norm = norm + f(j1,j2,j3,i_CO2)
            if (lH2O)  norm = norm + f(j1,j2,j3,i_H2O)
            if (lH2)   norm = norm + f(j1,j2,j3,i_H2)
            f(j1,j2,j3,minval(ichemspec):maxval(ichemspec)) &
                = f(j1,j2,j3,minval(ichemspec):maxval(ichemspec))/norm
!
!  Find mean molecular weight and density
!
            mu1 = f(j1,j2,j3,i_O2 )/mO2 &
                 +f(j1,j2,j3,i_H2O)/mH2O &
                 +f(j1,j2,j3,i_N2 )/mN2
            if (lH2)   mu1 = mu1+f(j1,j2,j3,i_H2  )/mH2
            if (lCO2)  mu1 = mu1+f(j1,j2,j3,i_CO2 )/mCO2
            if (lC3H8) mu1 = mu1+f(j1,j2,j3,i_C3H8)/mC3H8
            f(j1,j2,j3,ilnrho) = log(init_pressure)-log(Rgas)-f(j1,j2,j3,ilnTT)-log(mu1)
          enddo
        enddo
      enddo
!
!  Check if we want nolog of density
!
      if (ldensity_nolog)     f(:,:,:,irho) = exp(f(:,:,:,ilnrho))
      if (ltemperature_nolog) f(:,:,:,iTT) = exp(f(:,:,:,ilnTT))
!
    endsubroutine opposite_flames
!***********************************************************************
    subroutine opposite_ignitions(f)
!
!  03-jan-10/nilshau: adapted from flame_blob
!
!  Set up two oppositely directed flame fronts in the x-direction.
!  The two fronts have fresh gas between them.
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: j1, j2, j3
!
      real :: lower, upper, PP
      real :: T0, T1, rho0, rho1
!
      lflame_front = .true.
!
      call air_field(f,PP)
!
! Initialize some indexes
!
      do j3 = n1,n2
        do j2 = m1,m2
          do j1 = l1,l2
!
!  First define the distance from the lower and upper domain boundary.
!
            lower = x(j1)-xyz0(1)
            upper = xyz1(1)-x(j1)
!
!  Find temperature and density based on distance from boundaries.
!
            if (ltemperature_nolog) then
              T0 = f(j1,j2,j3,ilnTT)
            else
              T0 = exp(f(j1,j2,j3,ilnTT))
            endif
            if (ldensity_nolog) then
              rho0 = f(j1,j2,j3,ilnrho)
            else
              rho0 = exp(f(j1,j2,j3,ilnrho))
            endif
            T1 = (init_TT2-init_TT1)*exp(-(lower/init_x2)**2)+ &
                 (init_TT2-init_TT1)*exp(-(upper/init_x2)**2)+ &
                  init_TT1
            if (ltemperature_nolog) then
              f(j1,j2,j3,ilnTT) = T1
            else
              f(j1,j2,j3,ilnTT) = log(T1)
            endif
            rho1 = rho0*T0/T1
            if (ldensity_nolog) then
              f(j1,j2,j3,ilnrho) = rho1
            else
              f(j1,j2,j3,ilnrho) = log(rho1)
            endif
          enddo
        enddo
      enddo
!
    endsubroutine opposite_ignitions
!***********************************************************************
    subroutine calc_for_chem_mixture(f)
!
!  Calculate quantities for a mixture
!
!  22-jun-10/julien: Added evaluation of diffusion coefficients using constant
!                     Lewis numers Di = lambda/(rho*Cp*Lei)
!  10-jan-11/julien: Modified for a resolution with LSODE
!  26-jui-11/julien: Replaced fatal_error by inevitably_fatal_error to allow
!                    proper exit when T_loc<T_low or T_loc>T_up
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mmx) ::  tmp_sum, tmp_sum2, nu_dyn, nuk_nuj, Phi
      real, dimension(mmx) :: cp_R_spec, T_loc, T_loc_2, T_loc_3, T_loc_4
!
      integer :: k, j, j2, j3
      real :: T_up, T_mid, T_low
      real :: mk_mj
      real :: EE=0., TT=0., yH=1.
!
      logical, save :: lwrite=.true.
!
      character(len=fnlen) :: output_file="./data/mix_quant.out"
      character(len=15) :: writeformat
      integer :: file_id=123
      integer, parameter :: ii1=1, ii2=2, ii3=3, ii4=4, ii5=5
!
!  Density and temperature
!
      call timing('calc_for_chem_mixture','entered')
!
      call getdensity(f,EE,TT,yH,rho_full)
      call gettemperature(f,TT_full)
!
! Now this routine is only for chemkin data !!!
!
      if (lcheminp) then
!
        if (unit_system == 'cgs') then
       !   Rgas_unit_sys = k_B_cgs/m_u_cgs
       !   Rgas = Rgas_unit_sys/unit_energy
!
          call getmu_array(f,mu1_full)
!
          if (l1step_test) then
            species_constants(:,imass) = 1.
            mu1_full = 1.
          endif
!
!  Mole fraction XX
!
          do k = 1,nchemspec
            if (species_constants(k,imass) > 0.) &
              XX_full(ll1:ll2,mm1:mm2,nn1:nn2,k) = f(ll1:ll2,mm1:mm2,nn1:nn2,ichemspec(k))*unit_mass &
                      /(species_constants(k,imass)*mu1_full(ll1:ll2,mm1:mm2,nn1:nn2))
          enddo
!
!          do m=m1,m2
!            do n=n1,n2
!              call  getpressure(pp_full(l1:l2,m,n),TT_full(l1:l2,m,n),&
!                  rho_full(l1:l2,m,n),mu1_full(l1:l2,m,n))
!            enddo
!          enddo
!
!  Specific heat at constant pressure
!
          if ((Cp_const < impossible) .or. (Cv_const < impossible)) then
!
            if (Cp_const < impossible) then
              cp_full = Cp_const*mu1_full
              cv_full = (Cp_const-Rgas)*mu1_full
            endif
!
            if (Cv_const < impossible) then
              cv_full = Cv_const*mu1_full
            endif
          else
            cp_full = 0.
            cv_full = 0.
            do j3 = nn1,nn2
              do j2 = mm1,mm2
                T_loc = TT_full(ll1:ll2,j2,j3)
                T_loc_2 = T_loc*T_loc
                T_loc_3 = T_loc_2*T_loc
                T_loc_4 = T_loc_3*T_loc
                do k = 1,nchemspec
                  if (species_constants(k,imass) > 0.) then
                    T_low = species_constants(k,iTemp1)-1.
                    T_mid = species_constants(k,iTemp2)
                    T_up = species_constants(k,iTemp3)
!
! $                  if (j1<=l1 .or. j2>=l2) then
! $                    T_low=0.
! $                    T_up=1e10
! $                  endif
!
                    if (j2 <= m1 .or. j2 >= m2) then
                      T_low = 0.
                      T_up = 1e10
                    endif
!
                    if (j3 <= n1 .or. j3 >= n2) then
                      T_low = 0.
                      T_up = 1e10
                    endif
!
                    if (lcloud) T_low = 20.
                    where (T_loc >= T_low .and. T_loc <= T_mid)
                      cp_R_spec = species_constants(k,iaa2_offset+ii1) &
                          +species_constants(k,iaa2_offset+ii2)*T_loc &
                          +species_constants(k,iaa2_offset+ii3)*T_loc_2 &
                          +species_constants(k,iaa2_offset+ii4)*T_loc_3 &
                          +species_constants(k,iaa2_offset+ii5)*T_loc_4
                    elsewhere (T_loc >= T_mid .and. T_loc <= T_up)
                      cp_R_spec = species_constants(k,iaa1_offset+ii1) &
                          +species_constants(k,iaa1_offset+ii2)*T_loc &
                          +species_constants(k,iaa1_offset+ii3)*T_loc_2 &
                          +species_constants(k,iaa1_offset+ii4)*T_loc_3 &
                          +species_constants(k,iaa1_offset+ii5)*T_loc_4
                    endwhere
                    cv_R_spec_full(ll1:ll2,j2,j3,k) = cp_R_spec-1.
!
! Check if the temperature are within bounds
!
                    if (maxval(T_loc) > T_up .or. minval(T_loc) < T_low) then
                      print*,'iproc=',iproc
                      print*,'TT_full(ll1:ll2,j2,j3)=',T_loc
                      print*,'j2,j3=',j2,j3
                      call inevitably_fatal_error('calc_for_chem_mixture', &
                          'TT_full(ll1:ll2,j2,j3) is outside range', .true.)
                    endif
!
! Find cp and cv for the mixture for the full domain
!
                    cp_full(ll1:ll2,j2,j3) = cp_full(ll1:ll2,j2,j3)+f(ll1:ll2,j2,j3,ichemspec(k))  &
                                      *cp_R_spec/species_constants(k,imass)*Rgas
                    cv_full(ll1:ll2,j2,j3) = cv_full(ll1:ll2,j2,j3)+f(ll1:ll2,j2,j3,ichemspec(k))  &
                                      *cv_R_spec_full(ll1:ll2,j2,j3,k)/species_constants(k,imass)*Rgas

                    cp_spec_glo(ll1:ll2,j2,j3,k)=cp_R_spec/species_constants(k,imass)*Rgas

                  endif
                enddo
              enddo
            enddo
          endif
!
!  All the transport properties are calculated only if we are not using LSODE
!  to solve chemistry or during the transport substep
!
          if (.not. lchemonly) then
!
!  Viscosity of a mixture
!
            if (tran_exist) call calc_diff_visc_coef(f)
!
            if (visc_const == impossible) then
              do j3 = nn1,nn2
                do j2 = mm1,mm2
!
                  if  (lone_spec) then
                    f(ll1:ll2,j2,j3,iviscosity) = species_viscosity(ll1:ll2,j2,j3,1)/rho_full(ll1:ll2,j2,j3)
                  else
                    nu_dyn = 0.
                    do k = 1,nchemspec
                      if (species_constants(k,imass) > 0.) then
                        tmp_sum2 = 0.
                        do j = 1,nchemspec
                          mk_mj = species_constants(k,imass)/species_constants(j,imass)
                          nuk_nuj = species_viscosity(ll1:ll2,j2,j3,k)/species_viscosity(ll1:ll2,j2,j3,j)
                          Phi = 1./sqrt(8.)*1./sqrt(1.+mk_mj)*(1.+sqrt(nuk_nuj)*mk_mj**(-0.25))**2
                          tmp_sum2 = tmp_sum2 + XX_full(ll1:ll2,j2,j3,j)*Phi
                        enddo
                        nu_dyn = nu_dyn + XX_full(ll1:ll2,j2,j3,k)*species_viscosity(ll1:ll2,j2,j3,k)/tmp_sum2
                      endif
                    enddo
                    f(ll1:ll2,j2,j3,iviscosity) = nu_dyn/rho_full(ll1:ll2,j2,j3)
                  endif
!
                enddo
              enddo
            endif
!
!  Diffusion coefficient of a mixture from tran.dat file
!
            if ((.not. lDiff_simple).and.(.not. lDiff_lewis)) then
!
              do j3 = nn1,nn2
                do j2 = mm1,mm2
!
                  Diff_full(ll1:ll2,j2,j3,:) = 0.
                  if (.not. lone_spec) then
!
                    if (lfix_Sc) then
                      do k = 1,nchemspec
                        if (species_constants(k,imass) > 0.) &
                          Diff_full(ll1:ll2,j2,j3,k) = species_viscosity(ll1:ll2,j2,j3,k)/rho_full(ll1:ll2,j2,j3)/Sc_number
                      enddo
                    elseif (ldiffusion) then
!
! The mixture diffusion coefficient as described in eq. 5-45 of the Chemkin
! manual. Previously eq. 5-44 was used, but due to problems in the limit
! when the mixture becomes a pure specie we changed to the more robust eq. 5-45.
!
                      do k = 1,nchemspec
                        tmp_sum = 0.
                        tmp_sum2 = 0.
                        do j = 1,nchemspec
                          if (species_constants(k,imass) > 0.) then
                            if (j /= k) then
                              tmp_sum = tmp_sum + XX_full(ll1:ll2,j2,j3,j)/Bin_Diff_coef(ll1:ll2,j2,j3,j,k)
                              tmp_sum2 = tmp_sum2 + XX_full(ll1:ll2,j2,j3,j)*species_constants(j,imass)
                            endif
                          endif
                        enddo
                        Diff_full_add(ll1:ll2,j2,j3,k) = mu1_full(ll1:ll2,j2,j3)*tmp_sum2/tmp_sum
                      enddo
                    endif
                  endif
                enddo
              enddo
            endif
!
!  Thermal diffusivity
!
            if (lheatc_chemistry .and. (.not. lThCond_simple)) call calc_therm_diffus_coef
!
          endif
        else
          call fatal_error('calc_for_chem_mixture','lcheminp=T works only for cgs units')
        endif
      endif
!
!  Write block
!
      if (lwrite) then
        open (file_id,file=output_file)
        write (file_id,*) 'Mixture quantities'
        write (file_id,*) '*******************'
        write (file_id,*) ''
        write (file_id,*) 'Mass, g/mole'
        write (file_id,'(7E12.4)') 1./maxval(mu1_full/unit_mass)
        write (file_id,*) ''
        write (file_id,*) 'Density, g/cm^3'
        write (file_id,'(7E12.4)') rho_full(l1,m1,n1)*unit_mass/unit_length**3, &
            rho_full(l2,m2,n2)*unit_mass/unit_length**3
        write (file_id,*) ''
        write (file_id,*) 'Temperature, K'
! Commented the next line out because
! samples/2d-tests/chemistry_GrayScott apparently has no f(:,:,:,5)
        if (iTT > 0) then
          write (file_id,'(7E12.4)')  &
               f(l1,m1,n1,iTT)*unit_temperature, &
               f(l2,m2,n2,iTT)*unit_temperature
        else if (ilnTT > 0) then
          write (file_id,'(7E12.4)')  &
               exp(f(l1,m1,n1,ilnTT))*unit_temperature, &
               exp(f(l2,m2,n2,ilnTT))*unit_temperature
        endif
        write (file_id,*) ''
        write (file_id,*) 'Cp,  erg/mole/K'
        write (file_id,'(7E12.4)') cp_full(l1,m1,n1)/Rgas* &
            Rgas_unit_sys/mu1_full(l1,m1,n1)/unit_mass,cp_full(l2,m2,n2)/Rgas* &
            Rgas_unit_sys/mu1_full(l2,m2,n2)/unit_mass
        write (file_id,*) ''
        write (file_id,*) 'cp, erg/g/K'
        write (file_id,'(7E12.4)') cp_full(l1,m1,n1)/Rgas*Rgas_unit_sys,cp_full(l2,m2,n2)/Rgas*Rgas_unit_sys
        write (file_id,*) ''
        write (file_id,*) 'gamma,max,min'
        write (file_id,'(7E12.4)') cp_full(l1,m1,n1)/cv_full(l1,m1,n1), &
            cp_full(l2,m2,n2)/cv_full(l2,m2,n2)
        if (.not. lchemonly) then
          write (file_id,*) ''
          write (file_id,*) 'Species viscosity, g/cm/s,'
          do k = 1,nchemspec
             write (file_id,'(7E12.4)') species_viscosity(l1,m1,n1,k),species_viscosity(l2,m2,n2,k)
          enddo
          write (file_id,*) ''
          write (file_id,*) 'Thermal cond, erg/(cm K s),'
          write (file_id,'(7E12.4)') (lambda_full(l1,m1,n1)* &
              unit_energy/unit_time/unit_length/unit_temperature), &
              (lambda_full(l2,m2,n2)* &
              unit_energy/unit_time/unit_length/unit_temperature)
          write (file_id,*) ''
          write (file_id,*) 'Species  Diffusion coefficient, cm^2/s'
          if (.not. lDiff_simple) then
            do k = 1,nchemspec
              write (file_id,'(7E12.4)') &
                  Diff_full_add(l1,m1,n1,k)*unit_length**2/unit_time, &
                  Diff_full_add(l2,m2,n2,k)*unit_length**2/unit_time
            enddo
          endif
          if (lparticles_chemistry) then
            write (file_id,*) ''
            writeformat = '(  E12.4)'
            write (writeformat(2:3),'(I2)') nchemspec
            write (file_id,*) 'Mass fraction, -'
            write (file_id,writeformat) f(l1,m1,n1,ichemspec(1):ichemspec(nchemspec))
            write (file_id,*) ''
            write (file_id,writeformat) species_constants(:,imass)
          endif
        endif
        write (file_id,*) ''
        if (lroot) print*,'calc_for_chem_mixture: writing mix_quant.out file'
        close (file_id)
        lwrite = .false.
      endif
!
      call timing('calc_for_chem_mixture','finished')
!
    endsubroutine calc_for_chem_mixture
!***********************************************************************
    subroutine chemistry_before_boundary(f)
!
!  Calculate quantities for a chemical mixture
!
      real, dimension(mx,my,mz,mfarray) :: f

      if (ldustdensity .or. lnormalize_chemspec) &
           call chemspec_normalization(f)
!
!  Remove unphysical values of the mass fractions. This must be done
!  before the call to update_solid_cells in order to avoid corrections
!  within the solid structure.
!  There is no reason why this call should only be active in the
!  combination with solid cells....
!
      if (lsolid_cells .or. lnormalize_chemspec_N2) &
           call chemspec_normalization_N2(f)

    endsubroutine chemistry_before_boundary
!***********************************************************************
    subroutine chemspec_normalization(f)
!
!   20-sep-10/Natalia: coded
!   renormalization of the species
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz) :: sum_Y
      integer :: k
!
      sum_Y = 0.
      do k = 1,nchemspec
        sum_Y = sum_Y+f(:,:,:,ichemspec(k))
      enddo

      do k = 1,nchemspec
        f(:,:,:,ichemspec(k)) = f(:,:,:,ichemspec(k))/sum_Y
      enddo
!
    endsubroutine chemspec_normalization
!***********************************************************************
    subroutine astrobiology_data
!
!  Proceedure to read in stoichiometric matrices in explicit format for
!  forward and backward reations. For historical reasons this is referred
!  to as "astrobiology_data".
!
!  28-feb-08/axel: coded
!
      character(len=80) :: chemicals=''
      ! Careful, limits the absolut size of the input matrix !!!
      character(len=15) :: file1='chemistry_m.dat', file2='chemistry_p.dat'
!      character (len=fnlen) :: input_file='chem.inp'
      real :: dummy
      logical :: exist, exist1, exist2
      integer :: i, j, stat
      integer :: nchemspectemp
      character :: tmpchar
      logical :: inside
!
!  Find number of reactions by reading how many lines we have in file2
!
      j = 1
      open (19,file=file2)
      read (19,*) chemicals
      do while (.true.)
        read (19,*,end=996) dummy
        j = j+1
      enddo
996   close (19)
      mreactions = j-1
      if (lroot) print*,'Number of reactions=',mreactions
!
!  Find number of compounds by reading how many columns we have in file1
!
      open (19,file=file1)
      read (19,fmt="(a80)") chemicals
      nchemspectemp = 0
      inside = .true.
      do i = 1,len_trim(chemicals)
        tmpchar = chemicals(i:i)
        if (tmpchar == ' ') then
          if (.not. inside) then
            inside = .true.
            nchemspectemp = nchemspectemp+1
          endif
        else
          inside = .false.
        endif
      enddo
      if (inside) nchemspectemp = nchemspectemp-1
      close (19)
      if (lroot) print*,'Number of compounds=',nchemspectemp
      if (nchemspectemp > nchemspec) call fatal_error("astrobiology_data", &
        "Too many chemicals! Change NCHEMSPEC in src/cparam.local")
!
!  Allocate reaction arrays (but not during reloading!)
!
      if (.not. lreloading) then
        allocate(stoichio(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate stoichio")
        allocate(Sijm(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate Sijm")
        allocate(Sijp(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate Sijp")
        allocate(kreactions_z(mz,mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate kreactions_z")
        allocate(kreactions_p(mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate kreactions_p")
        allocate(kreactions_m(mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate kreactions_m")
        allocate(reaction_name(mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate reaction_name")
        allocate(enum_reaction_name(mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate enum_reaction_name")
        enum_reaction_name = 0
        allocate(kreactions_profile(mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate kreactions_profile")
        allocate(kreactions_profile_width(mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate kreactions_profile_width")
        allocate(kreactions_alpha(mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate kreactions_alpha")
        kreactions_alpha = 0.
        allocate(back(mreactions),STAT=stat)
        if (stat > 0) call fatal_error("astrobiology_data","Couldn't allocate back")
      endif
!
!  Initialize data
!
      kreactions_z = 1.
      Sijp = 0
      Sijm = 0
      back = .true.
!
!  read chemistry data
!
      inquire (file=file1,exist=exist1)
      inquire (file=file2,exist=exist2)
!
      if (exist1 .and. exist2) then
!
!  if both chemistry1.dat and chemistry2.dat are present,
!  then read Sijp and Sijm, and calculate their sum
!
!  file1
!
        open (19,file=file1)
        read (19,*) chemicals
        do j = 1,mreactions
          read (19,*,end=994) kreactions_m(j),(Sijm(i,j),i=1,nchemspectemp)
        enddo
994     close (19)
        nreactions1 = j-1
!
!  file2
!
        open (19,file=file2)
        read (19,*) chemicals
        do j = 1,mreactions
          if (lkreactions_profile) then
            if (lkreactions_alpha) then
              read (19,*) kreactions_p(j),(Sijp(i,j),i=1,nchemspectemp),kreactions_profile(j), &
                  kreactions_profile_width(j),kreactions_alpha(j)
            else
              read (19,*) kreactions_p(j),(Sijp(i,j),i=1,nchemspectemp),kreactions_profile(j), &
                  kreactions_profile_width(j)
            endif
          else
            if (lkreactions_alpha) then
              read (19,*) kreactions_p(j),(Sijp(i,j),i=1,nchemspectemp),kreactions_alpha(j)
            else
              read (19,*) kreactions_p(j),(Sijp(i,j),i=1,nchemspectemp)
            endif
          endif
        enddo
        close (19)
        nreactions2 = j-1
!
!  calculate stoichio and nreactions
!
        if (nreactions1 == nreactions2) then
          nreactions = nreactions1
          stoichio = Sijp-Sijm
        else
          call fatal_error("astrobiology_data",'nreactions1/=nreactions2')
        endif
        if (nreactions /= mreactions) call fatal_error("astrobiology_data",'nreactions/=mreactions')
!
      else
!
!  old method: read chemistry data, if present
!
        inquire (file='chemistry.dat',exist=exist)
        if (exist) then
          open (19,file='chemistry.dat')
          read (19,*) chemicals
          do j = 1,mreactions
            read (19,*,end=990) kreactions_p(j),(stoichio(i,j),i=1,nchemspec)
          enddo
990       close (19)
          nreactions = j-1
          Sijm = -min(stoichio,0.0)
          Sijp = +max(stoichio,0.0)
        else
          if (lroot) print*,'no chemistry.dat file to be read.'
          lreactions = .false.
        endif
      endif
!
!  print input data for verification
!
      if (lroot) then
!        print*,'chemicals=',chemicals
        print*,'kreactions_m=',kreactions_m(1:nreactions)
        print*,'kreactions_p=',kreactions_p(1:nreactions)
        print*,'Sijm:'
        do i = 1,nreactions
          print*,Sijm(:,i)
        enddo
        print*,'Sijp:'
        do i = 1,nreactions
          print*,Sijp(:,i)
        enddo
        print*,'stoichio='
        do i = 1,nreactions
          print*,stoichio(:,i)
        enddo
      endif
!
!  possibility of z-dependent kreactions_z profile
!
      if (lkreactions_profile) then
        do j = 1,nreactions
          if (kreactions_profile(j) == 'cosh') then
            do n = 1,mz
              kreactions_z(n,j) = 1./cosh(z(n)/kreactions_profile_width(j))**2
            enddo
          elseif (kreactions_profile(j) == 'gauss') then
            do n = 1,mz
              kreactions_z(n,j) = exp(-((z(n)/kreactions_profile_width(j))**2))
            enddo
          elseif (kreactions_profile(j) == 'square') then
            do n = 1,mz
              if (n < mz/2) then
                kreactions_z(n,j) = kreactions_profile_width(j)
              else
                kreactions_z(n,j) = 0.
              endif
            enddo
          elseif (kreactions_profile(j) == 'saw') then
            do n = 1,mz
              kreactions_z(n,j) = 0.51+(sin(pi*z(n)/kreactions_profile_width(j)) &
                  +sin(2*pi*z(n)/kreactions_profile_width(j))/2 &
                  +sin(3*pi*z(n)/kreactions_profile_width(j))/3 &
                  +sin(4*pi*z(n)/kreactions_profile_width(j))/4)/3
            enddo
          elseif (kreactions_profile(j) == 'sin') then
            do n = 1,mz
              kreactions_z(n,j) = 0.5*(1+cos(pi*z(n)/kreactions_profile_width(j)))
            enddo
          elseif (kreactions_profile(j) == 'sin-bg') then
            do n = 1,mz
              kreactions_z(n,j) = 0.5*(1.1+cos(pi*z(n)/kreactions_profile_width(j)))
            enddo
          elseif (kreactions_profile(j) == 'spike') then
            do n = 1,mz
              if (cos(pi*z(n)/kreactions_profile_width(j)) > 0.99) then
                kreactions_z(n,j) = 1
              else
                kreactions_z(n,j) = 0
              endif
            enddo
          endif
        enddo
      endif
!
    endsubroutine astrobiology_data
!***********************************************************************
    subroutine get_sum_DYDts(p,sum_DYDt,sum_hhk_DYDt_reac)

      type(pencil_case), intent(IN) :: p
      real, dimension(nx), intent(OUT) :: sum_DYDt, sum_hhk_DYDt_reac
      integer :: k

      sum_DYDt = 0.
      sum_hhk_DYDt_reac = 0.
      do k = 1,nchemspec
        if (species_constants(k,imass) > 0.) then
          sum_DYDt = sum_DYDt+Rgas/species_constants(k,imass)*(p%DYDt_reac(:,k)+p%DYDt_diff(:,k))
          if (lreactions) then
            sum_hhk_DYDt_reac = sum_hhk_DYDt_reac-p%hhk_full(:,k)*p%DYDt_reac(:,k)
          endif
        endif
      enddo
     endsubroutine
!***********************************************************************
     subroutine get_1step_test_sum_DYDts(f,p,sum_DYDt)
!
!  26-oct-25/TP: carved from dchemistry_dt
!
      real, dimension(mx,my,mz,mfarray) :: f
      type(pencil_case), intent(IN) :: p
      real, dimension(nx), intent(OUT) :: sum_DYDt

      integer :: i

      sum_DYDt = 0.
      !TP: cannot access 1 on the GPU
      do i = 1,nx
        sum_DYDt(i) = -p%rho(1)*(p%TT(i)-Tinf)*p%TT1(i) &
            *Cp_const/lambda_const*beta*(beta-1.)*f(l1,m,n,iux)*f(l1,m,n,iux)
      enddo
     endsubroutine get_1step_test_sum_DYDts
!***********************************************************************
    subroutine dchemistry_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   13-aug-07/steveb: coded
!    8-jan-08/natalia: included advection/diffusion
!   20-feb-08/axel: included reactions
!   22-jun-10/julien: modified evaluation of enthalpy fluxes with
!                     constant Lewis numbers
!   10-jan-11/julien: modified to solve chemistry with LSODE
!
      use Sub, only: grad,dot_mn, u_dot_grad_alt
      use Special, only: special_calc_chemistry
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx,3) :: dk_D, sum_diff=0., corr_vel, tmp1
      real, dimension(nx) :: ugchemspec, sum_DYDT, sum_dhhk=0., gXkgDk, gXkgmu1
      real, dimension(nx) :: sum_dk_ghk, dk_dhhk, sum_hhk_DYDt_reac, div_corr_vel
      real, dimension(nx) :: div_rho_Vc_Yk, rho_Yk_divVc, vc_grad_rho_yk
      type (pencil_case) :: p
      real, dimension(nx) :: RHS_T_full, diffus_chem
!
!  indices
!
      integer :: j, k,i,ii
      integer, parameter :: i1=1, i2=2, i3=3, i4=4, i5=5, i6=6, i7=7, i8=8, i9=9, i10=10
      integer, parameter :: i11=11, i12=12, i13=13, i14=14, i15=15, i16=16, i17=17, i18=18, i19=19
      integer, parameter :: iz1=1, iz2=2, iz3=3, iz4=4, iz5=5, iz6=6, iz7=7, iz8=8, iz9=9, iz10=10
      integer, parameter :: iz11=11, iz12=12, iz13=13, iz14=14, iz15=15, iz16=16, iz17=17
      integer, parameter :: iz18=18, iz19=19
!
      intent(in) :: p,f
      intent(inout) :: df
!
      logical :: ldiffusion2
      ldiffusion2 = ldiffusion .and. (.not. lchemonly)
!
!  identify module and boundary conditions
!
      call timing('dchemistry_dt','entered',mnloop=.true.)
      if (headtt .or. ldebug) print*,'dchemistry_dt: SOLVE dchemistry_dt'
!
!  Interface for your personal subroutines calls
      !
      if (lspecial) call special_calc_chemistry(f,df,p)
!
! Calculate correction velocity to ensure that the sum of all mass fractions is unit
!
      if (lcorr_vel) then
        corr_vel=0.0
        div_corr_vel=0.0
        do k = 1,nchemspec
          do i=1,3
            corr_vel(:,i)=corr_vel(:,i)+p%gXXk(:,i,k)*p%mukmu1(:,k)*p%Diff_penc_add(:,k)
          enddo
          call dot_mn(p%gXXk(:,:,k),p%gdiffk(:,:,k),gXkgDk)
          call dot_mn(p%gXXk(:,:,k),p%gmu1,gXkgmu1)
          div_corr_vel=div_corr_vel&
               +p%mukmu1(:,k)*gXkgDk&
               +p%mukmu1(:,k)*p%Diff_penc_add(:,k)*p%g2XXk(:,k) &
               +p%Diff_penc_add(:,k)*gXkgmu1*species_constants(k,imass)
        enddo
      endif
!
!  loop over all chemicals
!
      do k = 1,nchemspec
!
!  advection terms
!
        if (lhydro .and. ladvection .and.(.not. lchemonly)) then
!          call grad(f,ichemspec(k),gchemspec)
           !call dot_mn(p%uu,p%gYYk(:,:,k),ugchemspec)

           call u_dot_grad_alt(f,ichemspec(k),p%gYYk(:,:,k),p%uu,ugchemspec,iadv)
           
          if (lmobility) ugchemspec = ugchemspec*mobility(k)
          df(l1:l2,m,n,ichemspec(k)) = df(l1:l2,m,n,ichemspec(k))-ugchemspec
        endif
!
!  diffusion operator
!
!  Temporary we check the existence of chem.imp data,
!  further one should check the existence of a file with
!  binary diffusion coefficients!
!
        if (ldiffusion2) &
             df(l1:l2,m,n,ichemspec(k)) = df(l1:l2,m,n,ichemspec(k))+p%DYDt_diff(:,k)
!
!  Add correction velocity to account for the fact that the species mass fractions
!  are not guaranteed to sum to unity for any of the implemented diffusion operators.
!  If the correction velocity does not ensure unit mass fractions, one must set
!  lnormalize_chemspec or lnormalize_chemspec_N2 to true in run.in.
!
        if (lcorr_vel) then
          rho_Yk_divVc=p%rho*f(l1:l2,m,n,ichemspec(k))*div_corr_vel
          do i=1,3
            tmp1(:,i)=f(l1:l2,m,n,ichemspec(k))*p%grho(:,i)+p%rho*p%gYYk(:,i,k)
          enddo
          call dot_mn(corr_vel,tmp1,Vc_grad_rho_Yk)
          div_rho_Vc_Yk=rho_Yk_divVc+Vc_grad_rho_Yk
          df(l1:l2,m,n,ichemspec(k)) = df(l1:l2,m,n,ichemspec(k))-div_rho_Vc_Yk*p%rho1
        endif
!
!  chemical reactions:
!  multiply with stoichiometric matrix with reaction speed
!  d/dt(x_i) = S_ij v_j
!
        if (lreactions) then
          if (lchemonly) then
!  If chemistry is solved in a separate step, we want df to contain only the
!  chemical contribution, that's why no sum is required
            df(l1:l2,m,n,ichemspec(k)) = p%DYDt_reac(:,k)
          else
            df(l1:l2,m,n,ichemspec(k)) = df(l1:l2,m,n,ichemspec(k))+p%DYDt_reac(:,k)
          endif
        endif
!
!  Add filter for negative concentrations
!
        if (lfilter .and. .not. lfilter_strict) then
          do i = 1,mx
            if ((f(i,m,n,ichemspec(k))+df(i,m,n,ichemspec(k))*dt) < -1e-25 ) then
              df(i,m,n,ichemspec(k)) = -1e-25*dt
               !df(i,m,n,ichemspec(k)) = -0.99*f(i,m,n,ichemspec(k))/dt
            endif
            if ((f(i,m,n,ichemspec(k))+df(i,m,n,ichemspec(k))*dt) > 1. ) df(i,m,n,ichemspec(k)) = 1.*dt
          enddo
        endif
!
!  Add strict filter for negative concentrations
!
        if (lfilter_strict) then
          do i = 1,mx
            if ((f(i,m,n,ichemspec(k))+df(i,m,n,ichemspec(k))*dt) < 0.0 ) then
              if (df(i,m,n,ichemspec(k)) < 0.) df(i,m,n,ichemspec(k)) = 0.
            endif
            if ((f(i,m,n,ichemspec(k))+df(i,m,n,ichemspec(k))*dt) > 1. ) df(i,m,n,ichemspec(k)) = 1.*dt
          enddo
        endif
!
      enddo
!
      if (ldensity .and. lcheminp) then
!
        if (l1step_test) then
          call get_1step_test_sum_DYDts(f,p,sum_DYDt)
        else
          sum_dk_ghk = 0.
          call get_sum_DYDts(p,sum_DYDt,sum_hhk_DYDt_reac)
!
!  Reaction terms
!
          do k = 1,nchemspec
            if (species_constants(k,imass) > 0.) then
!
!  Sum over all species of diffusion terms
!
              if (ldiffusion2) then
                if (lDiff_fick) then
                  !call grad(f,ichemspec(k),gchemspec)
                  do i = 1,3
                    dk_D(:,i) = p%gYYk(:,i,k)*p%Diff_penc_add(:,k)
                  enddo
                elseif (lFlux_simple) then
                  do i = 1,3
                    dk_D(:,i) = p%gXXk(:,i,k)*p%Diff_penc_add(:,k)*p%mukmu1(:,k)
                  enddo
                else
                  do i = 1,3
                    dk_D(:,i) = (p%gXXk(:,i,k) &
                        +(XX_full(l1:l2,m,n,k)-f(l1:l2,m,n,ichemspec(k)))*p%glnpp(:,i)) &
                        *p%Diff_penc_add(:,k)*p%mukmu1(:,k)
                  enddo
                endif
!
                call dot_mn(dk_D,p%ghhk(:,:,k),dk_dhhk)
                sum_dk_ghk = sum_dk_ghk+dk_dhhk
                if (ldiff_corr) sum_diff(:,k) = sum_diff(:,k)+dk_D(:,k)
              endif
!
            endif
          enddo
!
! If the correction velocity is added
!
          if (ldiff_corr .and. ldiffusion2) then
            call fatal_error('dchemistry_dt', &
                'correction velocity is not properly implemented - please fix')
            do k = 1,nchemspec
              call dot_mn(sum_diff,p%ghhk(:,:,k),sum_dhhk)
              sum_dk_ghk(:) = sum_dk_ghk(:)-f(l1:l2,m,n,ichemspec(k))*sum_dhhk(:)
            enddo
          endif
        endif
!
        if (l1step_test) then
          RHS_T_full = sum_DYDt(:)
        else
          if (ltemperature_nolog) then
              RHS_T_full = p%cv1*((sum_DYDt(:)-Rgas*p%mu1*p%divu)*p%TT &
                  +sum_dk_ghk+sum_hhk_DYDt_reac)
   !       call stop_it('ltemperature_nolog case does not work now!')
          else
            if (lchemonly) then
              RHS_T_full = (sum_DYDt(:)+sum_hhk_DYDt_reac*p%TT1(:))*p%cv1
            else
              RHS_T_full = (sum_DYDt(:)-Rgas*p%mu1*p%divu)*p%cv1 &
                  +sum_dk_ghk*p%TT1(:)*p%cv1+sum_hhk_DYDt_reac*p%TT1(:)*p%cv1
            endif
          endif
        endif
!
        if (.not. lT_const) then
          if (lchemonly) then
!  If chemistry is solved in a separate step, we want df to contain only the
!  chemical contribution, that's why no sum is required
            df(l1:l2,m,n,ilnTT) = RHS_T_full
          else
            df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + RHS_T_full
          endif
        endif
!
        if (lheatc_chemistry .and.(.not. lchemonly)) call calc_heatcond_chemistry(f,df,p)
      endif
!
!  Atmosphere case
!
      if (lcloud .and.(.not. lchemonly)) then
!
        df(l1:l2,m,n,ichemspec(index_H2O)) = df(l1:l2,m,n,ichemspec(index_H2O)) - p%ccondens
        do i = 1,mx
          if ((f(i,m,n,ichemspec(index_H2O)) &
              +df(i,m,n,ichemspec(index_H2O))*dt) >= 1. ) df(i,m,n,ichemspec(index_H2O)) = 0.
          if ((f(i,m,n,ichemspec(index_H2O)) &
              +df(i,m,n,ichemspec(index_H2O))*dt) < 0. ) df(i,m,n,ichemspec(index_H2O)) = 0.
        enddo
      endif
!
! this damping zone is needed in a case of NSCBC
!
      if (ldamp_zone_for_NSCBC) call damp_zone_for_NSCBC(f,df)
!
!  For the timestep calculation, need maximum diffusion
!
      if (lupdate_courant_dt .and. (.not. lchemonly)) then
        if (.not. lcheminp) then
          diffus_chem = chem_diff*maxval(chem_diff_prefactor)*dxyz_2
        else
          diffus_chem=0.
          do j = 1,nx
            if (ldiffusion .and. .not. ldiff_simple) then
!
!--------------------------------------
!  This expression should be discussed
!--------------------------------------
!
              diffus_chem(j) = diffus_chem(j)+maxval(Diff_full_add(l1+j-1,m,n,1:nchemspec))*dxyz_2(j)
            else
              diffus_chem(j) = 0.
            endif
          enddo
        endif
        maxdiffus=max(maxdiffus,diffus_chem)
      endif
!
! NB: it should be discussed
!
      if (lupdate_courant_dt) then
        if (lreactions .and.(.not. llsode .or. lchemonly)) then
!
!  calculate maximum of *relative* reaction rate if decaying,
!  or maximum of absolute rate, if growing.
!
          if (lchem_cdtc) then
            reac_chem = 0.
            do k = 1,nchemspec
              reac_chem = max(reac_chem, &
                  abs(p%DYDt_reac(:,k)/max(f(l1:l2,m,n,ichemspec(k)),.001)))
            enddo
!
          elseif (lcheminp) then
            reac_chem = 0.
            !sum_reac_rate=0.
            do k = 1,nchemspec
              reac_chem = reac_chem+abs(p%DYDt_reac(:,k)/max(f(l1:l2,m,n,ichemspec(k)),0.001))
              !sum_reac_rate=sum_reac_rate+p%DYDt_reac(:,k)
            enddo
            if (maxval(reac_chem) > 1e11) reac_chem = 1e11   !MR: not where(...)?
          endif
        endif
      endif

      call timing('dchemistry_dt','before ldiagnos',mnloop=.true.)
      call calc_diagnostics_chemistry(f,p)
      call timing('dchemistry_dt','finished',mnloop=.true.)
      !
    endsubroutine dchemistry_dt
!***********************************************************************
    subroutine calc_diagnostics_chemistry(f,p)
!
!  Calculate diagnostic quantities
!
      use Diagnostics
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(nx) :: ff_condm,sum_DYDt,sum_hhk_DYDt_reac

      integer :: ii,k,j
!
      if (ldiagnos.and.lchemistry_diag) then
        do k = 1,nchemspec
          do j = 1,nreactions
            net_react_p(k,j) = net_react_p(k,j)+stoichio(k,j)*sum(vreactions_p(:,j))
            net_react_m(k,j) = net_react_m(k,j)+stoichio(k,j)*sum(vreactions_m(:,j))
          enddo
        enddo
      endif

      if (lreactions .and. ireac /= 0 .and. ((.not. llsode).or. lchemonly)) then
        !TP: sum_hhk_DYDt_reac is needed only if maux == nchemspec+1
        if (maux == nchemspec+1) call get_sum_DYDts(p,sum_DYDt,sum_hhk_DYDt_reac)
        if (llast) call get_reac_rate(sum_hhk_DYDt_reac,f,p)   ! updates f!!!
      endif

      if (ldustdensity .or. lparticles) then
         if (lnucleation .or. lcondensing_species) then
            f(l1:l2,m,n,isupsat)=p%chem_conc(:,ichem_cond_spec)&
               /max(p%conc_sat_spec,1e-20)
        endif
        if (lnucleation) then
           if (lpencil(i_nucl_rate) .or. lpencil(i_nucl_rmin)) then
            !
            ! Fill auxilliary array with radius and rate of nucleii. For use in
            ! insert_nucleii in particles_dust.f90 and for visualization
            !
            f(l1:l2,m,n,inucl)=p%nucl_rmin
            f(l1:l2,m,n,inucrate)=p%nucl_rate
          endif
        endif
      endif

      if (ldiagnos) then
        if (idiag_dtchem /= 0) call max_mn_name(reac_chem/cdtc,idiag_dtchem,l_dt=.true.)
!
!  WL: instead of hardcoding Y1-Y9, wouldn't it be possible
!      to have them all in the same array? The nbody
!      module, for instance, has idiag_xxq and idiag_vvq, which
!      allows the user to add output the positions and velocities
!      of as many particle they wants.
!  RP: Totally agree... I still have to expand manually these hard-coded
!      Y1-Y9 and chemspec-chemspec9 when needed, but this is just work-around...
!  AB: I also agree!
!
        do ii = 1,nchemspec
          call sum_mn_name(f(l1:l2,m,n,ichemspec(ii)),idiag_Ym(ii))
          if (idiag_rhoYm(ii)/=0) call sum_mn_name(p%rho*f(l1:l2,m,n,ichemspec(ii)),idiag_rhoYm(ii))
          call max_mn_name(f(l1:l2,m,n,ichemspec(ii)),idiag_Ymax(ii))
          if (idiag_Ymin(ii)/= 0) call max_mn_name(-f(l1:l2,m,n,ichemspec(ii)),idiag_Ymin(ii),lneg=.true.)
          if (idiag_TYm(ii)/= 0) &
            call sum_mn_name(max(1.-f(l1:l2,m,n,ichemspec(ii))/Ythresh(ii),0.),idiag_TYm(ii))
          if (idiag_diffm(ii)/= 0)   call sum_mn_name( Diff_full_add(l1:l2,m,n,ii),idiag_diffm(ii))
          if (idiag_diffmax(ii)/= 0) call max_mn_name( Diff_full_add(l1:l2,m,n,ii),idiag_diffmax(ii))
          if (idiag_diffmin(ii)/= 0) call max_mn_name(-Diff_full_add(l1:l2,m,n,ii),idiag_diffmin(ii),lneg=.true.)
        enddo
!
        call sum_mn_name(cp_full(l1:l2,m,n),idiag_cpfull)
        call sum_mn_name(cv_full(l1:l2,m,n),idiag_cvfull)
!
        if (idiag_mixfracmax/=0) call max_mn_name(f(l1:l2,m,n,imixfrac),idiag_mixfracmax)
        if (idiag_mixfracmin/=0) call max_mn_name(-f(l1:l2,m,n,imixfrac),idiag_mixfracmin,lneg=.true.)
        if (idiag_flameindmax/=0) call max_mn_name(f(l1:l2,m,n,iflameind),idiag_flameindmax)
        if (idiag_flameindmin/=0) call max_mn_name(-f(l1:l2,m,n,iflameind),idiag_flameindmin,lneg=.true.)
        call sum_mn_name(lambda_full(l1:l2,m,n),idiag_lambdam)
        call max_mn_name(lambda_full(l1:l2,m,n),idiag_lambdamax)
        if (idiag_lambdamin/=0) call max_mn_name(-lambda_full(l1:l2,m,n),idiag_lambdamin,lneg=.true.)
        if (idiag_alpham/=0) &
            call sum_mn_name(lambda_full(l1:l2,m,n)/(p%rho*cp_full(l1:l2,m,n)),idiag_alpham)
        if (idiag_alphamax/=0) &
            call max_mn_name(lambda_full(l1:l2,m,n)/(p%rho*cp_full(l1:l2,m,n)),idiag_alphamax)
        if (idiag_alphamin/=0) &
            call max_mn_name(-lambda_full(l1:l2,m,n)/(p%rho*cp_full(l1:l2,m,n)),idiag_alphamin,lneg=.true.)
        if (lnucleation) then
          call sum_mn_name(p%nucl_rmin,idiag_nuclrmin)
          call sum_mn_name(p%nucl_rate,idiag_nuclrate)
          call sum_mn_name(p%conc_sat_spec,idiag_conc_satm)
          !if (idiag_ffcondposm/= 0) call sum_mn_name(max(0.,p%ff_cond),idiag_ffcondposm)
          !if (idiag_ffcondnegm/= 0) call sum_mn_name(min(0.,p%ff_cond),idiag_ffcondnegm)
          !call sum_mn_name(p%ff_cond,idiag_ffcondm)
          call sum_mn_name(p%ff_nucl,idiag_ffnucl)
        endif
        if (.not. lnolatentheat) then
          if (idiag_latent_heat/= 0) call sum_mn_name(p%latent_heat,idiag_latent_heat)
        endif
        if (isupsat/=0) then
          if (idiag_supersat/= 0) call sum_mn_name(f(l1:l2,m,n,isupsat),idiag_supersat)
        endif
!
!  Sample for hard coded diffusion diagnostics
!
!        call sum_mn_name(Diff_full(l1:l2,m,n,i1),idiag_diff1m)

        if (lreactions .and. lpencil(i_DYDt_reac) .and. (.not. llsode .or. lchemonly)) then
          do ii=1,nchemspec
            call sum_mn_name(p%DYDt_reac(:,ii),idiag_dYm(ii))
            if (idiag_dYmax(ii) /= 0) call max_mn_name(abs(p%DYDt_reac(:,ii)),idiag_dYmax(ii))
            if (idiag_hm(ii)    /= 0) call sum_mn_name(p%H0_RT(:,ii)*Rgas* &
                                           p%TT(:)/species_constants(ii,imass),idiag_hm(ii))
          enddo
        endif
!
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
        do ii = 1,nchemspec
          call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(ii)),idiag_Ymz(ii))
        enddo
      endif

    endsubroutine calc_diagnostics_chemistry
!***********************************************************************
    subroutine rprint_chemistry(lreset,lwrite)
!
!  reads and registers print parameters relevant to chemistry
!
!  13-aug-07/steveb: coded
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
      use General, only: itoa, get_species_nr
!
      integer :: iname, inamez,ii
      logical :: lreset, lwr
      logical, optional :: lwrite
      character(len=6) :: diagn_Ym, number
      character(len=6) :: diagn_Ymax
      character(len=6) :: diagn_Ymin
      character(len=7) :: diagn_dYmax
      character(len=7) :: diagn_rhoYm
      character(len=6) :: diagn_TYm
      character(len=6) :: diagn_dYm
      character(len=6) :: diagn_hm
      character(len=6) :: diagn_cpm
      character(len=7) :: diagn_diffm
      character(len=8) :: diagn_diffmax
      character(len=8) :: diagn_diffmin
      character(len=fmtlen) :: sname
!
      lwr = .false.
      if (present(lwrite)) lwr = lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtchem = 0
        idiag_Ym = 0
        idiag_rhoYm = 0
        idiag_TYm = 0
        idiag_dYm = 0
        idiag_Ymax = 0
        idiag_Ymin = 0
        idiag_dYmax = 0
        idiag_hm = 0
        idiag_cpm = 0
        idiag_diffm = 0
        idiag_diffmax = 0
        idiag_diffmin = 0
!
        idiag_cpfull = 0
        idiag_cvfull = 0
        idiag_e_intm = 0
        idiag_Ymz = 0
!
        idiag_mixfracmax = 0
        idiag_mixfracmin = 0
        idiag_flameindmax = 0
        idiag_flameindmin = 0
        idiag_lambdam = 0
        idiag_lambdamax = 0
        idiag_lambdamin = 0
        idiag_alpham = 0
        idiag_alphamax = 0
        idiag_alphamin = 0
        idiag_ffnucl = 0
        idiag_supersat = 0
        idiag_latent_heat = 0
!
        idiag_nuclrmin=0
        idiag_nuclrate=0
        idiag_conc_satm=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname = 1,nname
        do ii=1,nchemspec
          write (number,'(I2)') ii
          diagn_Ym = 'Y'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_Ym),idiag_Ym(ii))
          diagn_rhoYm = 'rhoY'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_rhoYm),idiag_rhoYm(ii))
          diagn_Ymax = 'Y'//trim(adjustl(number))//'max'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_Ymax),idiag_Ymax(ii))
          diagn_Ymin = 'Y'//trim(adjustl(number))//'min'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_Ymin),idiag_Ymin(ii))
          diagn_dYm = 'dY'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_dYm),idiag_dYm(ii))
          diagn_TYm = 'TY'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_TYm),idiag_TYm(ii))
          diagn_dYmax = 'dY'//trim(adjustl(number))//'max'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_dYmax),idiag_dYmax(ii))
          diagn_hm = 'h'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_hm),idiag_hm(ii))
          diagn_cpm = 'cp'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_cpm),idiag_cpm(ii))
          diagn_diffm = 'diff'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_diffm),idiag_diffm(ii))
          diagn_diffmax = 'diff'//trim(adjustl(number))//'max'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_diffmax),idiag_diffmax(ii))
          diagn_diffmin = 'diff'//trim(adjustl(number))//'min'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_diffmin),idiag_diffmin(ii))
        enddo
        call parse_name(iname,cname(iname),cform(iname),'dtchem',idiag_dtchem)
        call parse_name(iname,cname(iname),cform(iname),'cpfull',idiag_cpfull)
        call parse_name(iname,cname(iname),cform(iname),'cvfull',idiag_cvfull)
        call parse_name(iname,cname(iname),cform(iname),'nuclrmin',idiag_nuclrmin)
        call parse_name(iname,cname(iname),cform(iname),'nuclrate',idiag_nuclrate)
        call parse_name(iname,cname(iname),cform(iname),'conc_satm',idiag_conc_satm)
!
!   Sample for hard-coded heat capacity diagnostics
!
!        call parse_name(iname,cname(iname),cform(iname),'cp1m',idiag_cp1m)
        call parse_name(iname,cname(iname),cform(iname),'mixfracmax',idiag_mixfracmax)
        call parse_name(iname,cname(iname),cform(iname),'mixfracmin',idiag_mixfracmin)
        call parse_name(iname,cname(iname),cform(iname),'flameindmax',idiag_flameindmax)
        call parse_name(iname,cname(iname),cform(iname),'flameindmin',idiag_flameindmin)
        call parse_name(iname,cname(iname),cform(iname),'e_intm',idiag_e_intm)
        call parse_name(iname,cname(iname),cform(iname),'lambdam',idiag_lambdam)
        call parse_name(iname,cname(iname),cform(iname),'lambdamax',idiag_lambdamax)
        call parse_name(iname,cname(iname),cform(iname),'lambdamin',idiag_lambdamin)
        call parse_name(iname,cname(iname),cform(iname),'alpham',idiag_alpham)
        call parse_name(iname,cname(iname),cform(iname),'alphamax',idiag_alphamax)
        call parse_name(iname,cname(iname),cform(iname),'alphamin',idiag_alphamin)
        call parse_name(iname,cname(iname),cform(iname),'ffnucl',idiag_ffnucl)
        call parse_name(iname,cname(iname),cform(iname),'supersat',idiag_supersat)
        call parse_name(iname,cname(iname),cform(iname),'latent_heat',idiag_latent_heat)
!
!  Sample for hard-coded diffusion diagnostics
!
!        call parse_name(iname,cname(iname),cform(iname),'diff1m',idiag_diff1m)
      enddo
!
!  xy-averages
!
      do inamez = 1,nnamez
        do ii=1,nchemspec
          write (number,'(I2)') ii
          diagn_Ym = 'Y'//trim(adjustl(number))//'mz'
          call parse_name(inamez,cnamez(inamez),cformz(inamez),trim(diagn_Ym),idiag_Ymz(ii))
        enddo
      enddo
!
!  check for those quantities for which we want video slices
!
      do iname=1,nnamev
        sname=trim(cnamev(iname))
        if (sname(1:8)=='chemspec') then
          if (get_species_nr(sname,'chemspec',nchemspec,'rprint_chemistry')>0) cformv(iname)='DEFINED'
        endif
      enddo
!
      if (lwrite_slices) then
        where(cnamev=='nuclrmin') cformv='DEFINED'
        where(cnamev=='nuclrate') cformv='DEFINED'
        where(cnamev=='supersat') cformv='DEFINED'
        where(cnamev=='mixfrac') cformv='DEFINED'
        where(cnamev=='flameind') cformv='DEFINED'
      endif

!
    endsubroutine rprint_chemistry
!***********************************************************************
    subroutine get_slices_chemistry(f,slices)
!
!  Write slices for animation of Chemistry variables.
!
!  13-aug-07/steveb: dummy
!  16-may-09/raphael: added more slices
!
      use Slices_methods, only: assign_slices_scal

      real, dimension(mx,my,mz,mfarray) :: f
      type (slice_data) :: slices

      character(len=fmtlen) :: sname
      integer :: ispec
!
!  Chemical species mass fractions.
!
      sname=trim(slices%name)
      if (sname(1:8)=='chemspec') then
        if (sname(9:)==' ') then    ! 9=len('chemspec')+1
          ispec=1
        else
          read(sname(9:),'(i3)') ispec
        endif
!
        call assign_slices_scal(slices,f,ichemspec(ispec))
      endif
!
!  Nucleation radius
!
      if (sname=='nuclrmin') then
        call assign_slices_scal(slices,f,inucl)
      endif
!
!  Nucleation radte
!
      if (sname=='nuclrate') then
        call assign_slices_scal(slices,f,inucrate)
      endif
!
!  Supersaturation
!
      if (sname=='supersat') then
        call assign_slices_scal(slices,f,isupsat)
      endif
!
!  Mixture fraction
!
      if (sname=='mixfrac') then
        call assign_slices_scal(slices,f,imixfrac)
      endif
!
!  Flame index
!
      if (sname=='flameind') then
        call assign_slices_scal(slices,f,iflameind)
      endif

!
    endsubroutine get_slices_chemistry
!***********************************************************************
    subroutine read_chemistry_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read (parallel_unit, NML=chemistry_init_pars, IOSTAT=iostat)
!
    endsubroutine read_chemistry_init_pars
!***********************************************************************
    subroutine write_chemistry_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit, NML=chemistry_init_pars)
!
    endsubroutine write_chemistry_init_pars
!***********************************************************************
    subroutine read_chemistry_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read (parallel_unit, NML=chemistry_run_pars, IOSTAT=iostat)
!
    endsubroutine read_chemistry_run_pars
!***********************************************************************
    subroutine write_chemistry_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit, NML=chemistry_run_pars)
!
    endsubroutine write_chemistry_run_pars
!***********************************************************************
    subroutine build_stoich_matrix(StartInd,StopInd,k,ChemInpLine,product)
!
!  calculation of the stoichoimetric matrix
!
!  10-mar-08/nils: coded
!
      integer, intent(in) :: StartInd, StopInd,k
      character(len=*), intent(in) :: ChemInpLine
      logical, intent(in) :: product
      integer :: StartSpecie, ind_glob, ind_chem
      real :: stoi
      logical :: found_specie, lreal=.false.
      integer :: stoi_int
!
      stoi = 0.
      if ((ChemInpLine(StartInd:StopInd) /= "M" ) &
          .and. (ChemInpLine(StartInd:StartInd+1) /= "hv" )) then
        StartSpecie = verify(ChemInpLine(StartInd:StopInd),"1234567890")+StartInd-1
!
!  Call to a routine that checks for arbitrary stoiciometric coefficents
!  removes them and shifts the begin of the species index in cheminpline
!
        if (StopInd-StartSpecie >= 3 .and. StartSpecie > 1) &
            call find_remove_real_stoic(ChemInpLine(StartSpecie-1:StopInd),lreal,stoi,StartSpecie)
!
        call find_species_index(ChemInpLine(StartSpecie:StopInd),ind_glob,ind_chem,found_specie)
        if (.not. found_specie) then
          print*,'ChemInpLine=',ChemInpLine
          print*,'StartSpecie,StopInd=',StartSpecie,StopInd
          print*,'ChemInpLine(StartSpecie:StopInd)=',ChemInpLine(StartSpecie:StopInd)
          print*,'ind_glob,ind_chem=',ind_glob,ind_chem
!          if (.not. lpencil_check_small) then
!          if (.not. lpencil_check) then
           call fatal_error("build_stoich_matrix","Did not find species")
!          endif
!          endif
        endif
!
!        if (found_specie) then
        if (StartSpecie == StartInd) then
          stoi = 1.0
        else
          if (.not. lreal) read (unit=ChemInpLine(StartInd:StartInd),fmt='(I1)') stoi_int
        endif
      endif
!
      if (product) then
        Sijm(ind_chem,k) = Sijm(ind_chem,k)+stoi
      else
        Sijp(ind_chem,k) = Sijp(ind_chem,k)+stoi
      endif
!        endif
!
    endsubroutine build_stoich_matrix
!***********************************************************************
    subroutine write_reactions
!
!  write reaction coefficient in the output file
!
!  11-mar-08/nils: coded
!
      use General, only: itoa
!
      integer :: reac, spec
      character(len=80) :: reac_string, product_string, output_string
      character(len=intlen) :: Sijp_string, Sijm_string
      character(len=1) :: separatorp, separatorm
      character(len=fnlen) :: input_file="./data/chem.out"
      integer :: file_id=123
!
      open (file_id,file=input_file,POSITION='APPEND',FORM='FORMATTED')
      write (file_id,*) 'REACTIONS'
      !open(file_id,file=input_file)
!
      do reac = 1,mreactions
        reac_string = ''
        product_string = ''
        separatorp = ''
        separatorm = ''
!
        do spec = 1,nchemspec
          if (Sijp(spec,reac) > 0) then
            Sijp_string = ''
            if (Sijp(spec,reac) /= 1) write (Sijp_string,'(F3.1)') Sijp(spec,reac)
!            if (Sijp(spec,reac)>1) Sijp_string=itoa(Sijp(spec,reac))
            reac_string = trim(reac_string)//trim(separatorp)// &
                          trim(Sijp_string)//trim(varname(ichemspec(spec)))
            separatorp = '+'
          endif
          if (Sijm(spec,reac) > 0) then
            Sijm_string = ''
            if (Sijm(spec,reac) /= 1) write (Sijm_string,'(F3.1)') Sijm(spec,reac)
!            if (Sijm(spec,reac)>1) Sijm_string=itoa(Sijm(spec,reac))
            product_string = trim(product_string)//trim(separatorm)// &
                             trim(Sijm_string)//trim(varname(ichemspec(spec)))
            separatorm = '+'
          endif
        enddo
!
        output_string = trim(reac_string)//'='//trim(product_string)
!
        if (.not. photochem_case(reac)) then
!
!  Note that since the B_n term is in logarithmic form within the code
!  the exponential must be used for output.
!
          write (unit=output_string(30:45),fmt='(E14.4)') exp(B_n(reac))
          write (unit=output_string(47:62),fmt='(E14.4)') alpha_n(reac)
          write (unit=output_string(64:79),fmt='(E14.4)') E_an(reac)
        endif
        write (file_id,*) trim(output_string)
        if (.not. photochem_case(reac)) then
          if (maxval(abs(low_coeff(:,reac))) > 0.) then
            write (file_id,*) 'LOW/',exp(low_coeff(1,reac)),low_coeff(2:,reac)
          elseif (maxval(abs(high_coeff(:,reac))) > 0.) then
            write (file_id,*) 'HIGH/',high_coeff(:,reac)
          endif
          if (maxval(abs(troe_coeff(:,reac))) > 0.) then
            write (file_id,*) 'TROE/',troe_coeff(:,reac)
          endif
          if (minval(a_k4(:,reac)) < impossible) then
            write (file_id,*) "a_k4/", a_k4(:,reac)
          endif
        else
          write (file_id,*) ' min lambda=',lamb_low,' max lambda=',lamb_up
        endif
!
      enddo
!
      write (file_id,*) 'END'
      write (file_id,*) '(M+) case: ',Mplus_case
      write (file_id,*) 'photochemical case: ',photochem_case
!
      close (file_id)
!
    endsubroutine write_reactions
!***********************************************************************
    subroutine chemkin_data(f)
!
!  if the file with chemkin data exists
!  reading the Chemkin data
!
!  21-jul-10/julien: Reading lewis.dat file to collect constant Lewis numbers
!                     for each species
!
      character(len=fnlen) :: input_file
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: stat, k,i
      character(len=fnlen) :: input_file2="./data/stoich.out"
      integer :: file_id=123
      logical :: chemin,cheminp
!
      inquire(file='chem.inp',exist=cheminp)
      inquire(file='chem.in',exist=chemin)
      if (chemin .and. cheminp) call fatal_error('chemistry', &
          'chem.in and chem.inp found, please decide for one')
      if (cheminp) input_file='chem.inp'
      if (chemin) input_file='chem.in'
!
      inquire (file='tran.dat',exist=tran_exist)
      if (.not. tran_exist) inquire (file='tran.in',exist=tran_exist)

      inquire (file='lewis.dat',exist=lew_exist)
!
!  Allocate binary diffusion coefficient array
!
!      if (.not. lreloading) then
!Natalia:
!this does not work for ldiffusion=F
!        if (ldiffusion .and. .not. lfix_Sc) then
        if (.not. lfix_Sc) then
!NILS: Since Bin_diff_coeff is such a huge array we must check if it
!NILS: required to define it for the full domain!!!!!!
          if (.not. lreloading) then
            allocate(Bin_Diff_coef(mx,my,mz,nchemspec,nchemspec),STAT=stat)
            if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate "// &
                                           "binary diffusion coefficients")
!
            allocate(Diff_full(mx,my,mz,nchemspec),STAT=stat)
            if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate "// &
                                           "diffusion coefficients Diff_full")

            allocate(Diff_full_add(mx,my,mz,nchemspec),STAT=stat)
            if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate "// &
                                           "diffusion coefficients Diff_full_add")
          endif
!
        endif
!      endif
!
      if (tran_exist) then
        if (lroot) print*,'tran.in/dat file with transport data is found.'
        call read_transport_data
      endif
!
      if (lew_exist) then
        if (lroot) then
          print*,'lewis.dat file with transport data is found.'
          print*,'Species diffusion coefficients calculated using constant Lewis numbers.'
        endif
        call read_Lewis
      endif
!
      if (.not. lew_exist .and. lDiff_lewis .and. lroot) then
        if (.not. l1step_test) then
          print*, 'No lewis.dat file present, switch to simplified diffusion'
          lDiff_lewis = .false.
          lDiff_simple = .true.
        else
          print*, 'Le=1'
          lDiff_lewis = .true.
        endif
      endif
!
      if (lroot .and. .not. tran_exist .and. .not. lew_exist) then
        if (chem_diff == 0.) call inevitably_fatal_error('chemkin data', 'chem_diff = 0')
        print*,'tran.dat file with transport data is not found.'
        print*,'lewis.dat file with Lewis numbers is not found.'
        print*,'Now diffusion coefficients is ',chem_diff
        print*,'Now species viscosity is ',nu_spec
        Bin_Diff_coef = chem_diff/(unit_length*unit_length/unit_time)
        do k = 1,nchemspec
          species_viscosity(:,:,:,k) = nu_spec(k)/(unit_mass/unit_length/unit_time)
        enddo
      endif
!
!  Find number of ractions
!
      call read_reactions(input_file,NrOfReactions=mreactions)
      if (lroot) print*,'Number of reactions=',mreactions
      if (lroot) print*,'Number of species=',nchemspec
      nreactions = mreactions
!
!  Allocate reaction arrays
!
      if (.not. lreloading) then
        allocate(stoichio(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate stoichio")
        allocate(Sijm(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate Sijm")
        allocate(Sijp(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate Sijp")
        allocate(Sijm_(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate Sijm_")
        allocate(Sijp_(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate Sijp_")
        allocate(Sijm_mod(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate Sijm_mod")
        allocate(Sijp_mod(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate Sijp_mod")
        allocate(reaction_name(mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate reaction_name")
        allocate(back(mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate back")
        allocate(B_n(mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate B_n")
        B_n = 0.
        allocate(alpha_n(mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate alpha_n")
        alpha_n = 0.
        allocate(E_an(mreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate E_an")
        E_an = 0.
!
        allocate(low_coeff(3,nreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate low_coeff")
        allocate(low_coeff_abs_max(nreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate low_coeff_abs_max")
        low_coeff = 0.
        low_coeff_abs_max = 0.
        allocate(high_coeff(3,nreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate high_coeff")
        allocate(high_coeff_abs_max(nreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate high_coeff_abs_max")
        high_coeff = 0.
        high_coeff_abs_max = 0.
        allocate(troe_coeff(3,nreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate troe_coeff")
        allocate(troe_coeff_abs_max(nreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate troe_coeff_abs_max")
        troe_coeff = 0.
        troe_coeff_abs_max = 0.
        allocate(a_k4(nchemspec,nreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate a_k4")
        allocate(a_k4_min(nreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate a_k4_min")
        a_k4 = impossible
        a_k4_min = impossible
        allocate(Mplus_case (nreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate Mplus_case")
        Mplus_case = .false.
        allocate(photochem_case (nreactions),STAT=stat)
        if (stat > 0) call fatal_error('chemkin_data',"Couldn't allocate photochem_case")
        photochem_case = .false.
      endif
!
!  Initialize data
!
      Sijp = 0.
      Sijm = 0.
      Sijp_mod = 0.
      Sijm_mod = 0.
      back = .true.
!
!  read chemistry data
!
      call read_reactions(input_file)
      call write_reactions
!
!  Possibility of modification for the power part in case of global (not detailed) reactions.
!
      if (lchem_detailed) then
        Sijp_=Sijp
        Sijm_=Sijm
      else
        call read_reactions_mod
        Sijp_=Sijp+Sijp_mod
        Sijm_=Sijm+Sijm_mod
        call write_matrices
      endif
!
!  calculate stoichio and nreactions
!
      stoichio = Sijp-Sijm
!
!  print input data for verification
!
      if (lroot .and. nreactions > 0) then
!
        open (file_id,file=input_file2,POSITION='rewind',FORM='FORMATTED')
        write (file_id,*) 'STOICHIOMETRIC MATRIX'
!
        write (file_id,*) 'Species names'
        write (file_id,101) varname(ichemspec(:))
!
        write (file_id,*) 'Sijm'
        do i = 1,nreactions
          write (file_id,100) i,Sijm(:,i)
        enddo
        write (file_id,*) 'Sijp:'
        do i = 1,nreactions
          write (file_id,100) i,Sijp(:,i)
        enddo
        write (file_id,*) 'stoichio='
        do i = 1,nreactions
          write (file_id,100) i,stoichio(:,i)
        enddo
        close (file_id)
      endif
!
100   format(I3,60f6.1)
101   format('    ',60A6)
!
      call keep_compiler_quiet(f)
!
    endsubroutine chemkin_data
!***********************************************************************
    subroutine read_reactions(input_file,NrOfReactions)
!
!  This subroutine reads all reaction information from chem.inp
!  See the chemkin manual for more information on
!  the syntax of chem.inp.
!
!  10-mar-08/nils: coded
!
      logical :: IsReaction=.false., found_new_reaction=.false.
      logical, save :: find_specie, found_specie
      integer, optional :: NrOfReactions
      integer, save :: ind_glob, ind_chem
      integer :: i, k, file_id=123, StartInd, StopInd, StartInd_add
      integer :: StopInd_add, StopInd_add_, StopIndName
      integer :: VarNumber, VarNumber_add, SeparatorInd
      integer :: PlusInd
      integer :: LastLeftCharacter, ParanthesisInd, Mplussind
      integer :: photochemInd, plusind_
      character(len=120) :: ChemInpLine, ChemInpLine_add
      character(len=*) :: input_file

      if (present(NrOfReactions)) NrOfReactions=0
      k=1
      open(file_id,file=input_file)
      dataloop3: do
        if (found_new_reaction) then
          ChemInpLine=ChemInpLine_add
        else
          read(file_id,'(80A)',end=1012) ChemInpLine
        endif
        found_new_reaction=.false.
!
! Check if we are reading a line within the reactions section
!
        if (ChemInpLine(1:9)=="REACTIONS")            IsReaction=.true.
        if (ChemInpLine(1:3)=="END" .and. IsReaction) IsReaction=.false.

        if (present(NrOfReactions)) then
!
! Find number of reactions
!
          if (IsReaction) then
            if (ChemInpLine(1:9) /= "REACTIONS") then
              StartInd=1; StopInd =0
              StopInd=index(ChemInpLine(StartInd:),'=')+StartInd-1
              if (StopInd>0 .and. ChemInpLine(1:1) /= '!') then
                NrOfReactions=NrOfReactions+1
              endif
            endif
          endif
        else
!
! Read in species
!
          if (IsReaction) then
            if (ChemInpLine(1:9) /= "REACTIONS") then
              StartInd=1; StopInd =0
              StopInd=index(ChemInpLine(StartInd:),'=')+StartInd-1
              if (StopInd>0 .and. ChemInpLine(1:1) /= '!') then
!
! Fill in reaction name
!
                StopIndName=index(ChemInpLine(StartInd:),' ')+StartInd-1
                reaction_name(k)=ChemInpLine(StartInd:StopIndName)
!
! Photochemical case
!
               ParanthesisInd=0
               photochemInd=0
               SeparatorInd=0
               photochemInd=index(ChemInpLine(StartInd:),'hv')
               if (photochemInd>0) then
                 photochem_case (k)=.true.

               ParanthesisInd=index(ChemInpLine(photochemInd:),'(') &
                             +photochemInd-1


               if ((ParanthesisInd>0) .and. (photochemInd>0)) then
                StopInd=index(ChemInpLine(StartInd:),'lam')

                SeparatorInd=index(ChemInpLine(ParanthesisInd:StopInd),'<')

                if (SeparatorInd>0) then
                  SeparatorInd=SeparatorInd+ParanthesisInd-1
                  read (unit=ChemInpLine(ParanthesisInd+1:&
                                 SeparatorInd-1),fmt='(E15.8)') lamb_low
                endif
                ParanthesisInd=index(ChemInpLine(StopInd:),')') +StopInd-1
                SeparatorInd=index(ChemInpLine(StopInd:ParanthesisInd),'<') &
                             +StopInd-1
                 if (SeparatorInd>0) then
                   read (unit=ChemInpLine(SeparatorInd+1:&
                                  ParanthesisInd-1),fmt='(E15.8)') lamb_up
                 endif
               endif
               endif
!
! End of the photochemical case
!
! Find reactant side stoichiometric coefficients (if =>, then backward reaction to false)
!
                SeparatorInd=index(ChemInpLine(StartInd:),'<=')
                if (SeparatorInd==0) then
                  SeparatorInd=index(ChemInpLine(StartInd:),'=')
                  if (index(ChemInpLine(StartInd:),'=>')/=0) back(k)=.false.
                endif

                ParanthesisInd=0
                MplussInd=0

                ParanthesisInd=index(ChemInpLine(StartInd:),'(+M)')
                MplussInd=index(ChemInpLine(StartInd:),'+M')

                found_new_reaction=.false.
                if (ParanthesisInd>0 .or. Mplussind>0) then
                  if (ParanthesisInd>0) then
                    LastLeftCharacter=min(ParanthesisInd,SeparatorInd)-1
                  else
                    LastLeftCharacter=min(MplussInd,SeparatorInd)-1
                  endif
!
! reading of the additional data for (+M) case
!
100               read(file_id,'(80A)',end=1012) ChemInpLine_add
                  if (ChemInpLine_add(1:1) == ' ') then


                    if (ParanthesisInd>0) then
                      Mplus_case (k)=.true.
                    endif

                    i=1
                    do while (i<80)
                      if (ChemInpLine_add(i:i)==' ') then
                        i=i+1
                      elseif (ChemInpLine_add(i:i+2)=='LOW') then
                      !  if (lroot) print*,ChemInpLine_add(i:i+2),&
                      !      '   coefficients for reaction ', &
                      !      reaction_name(k),'number ', k
                        VarNumber_add=1; StartInd_add=i+4; StopInd_add=i+4
                        do while (VarNumber_add<4)
                          StopInd_add=index(ChemInpLine_add(StartInd_add:),&
                              ' ')+StartInd_add-2
                          StopInd_add_=index(ChemInpLine_add(StartInd_add:),&
                              '/')+StartInd_add-2
                          StopInd_add=min(StopInd_add,StopInd_add_)
                          if (StopInd_add==StartInd_add) then
                            StartInd_add=StartInd_add+1
                          else
                            if (VarNumber_add==1) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E25.8)') low_coeff(1,k)
                              if (low_coeff(1,k)/=0.) low_coeff(1,k)=log(low_coeff(1,k))
                            elseif (VarNumber_add==2) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') low_coeff(2,k)
                            elseif (VarNumber_add==3) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') low_coeff(3,k)
                            else
                              call fatal_error("read_reactions","no such VarNumber: "// &
                                               trim(itoa(VarNumber_add)))
                            endif
                          endif
                          VarNumber_add=VarNumber_add+1
                          !StartInd_add=StopInd_add
                          StartInd_add=verify(ChemInpLine_add(StopInd_add+1:),&
                              ' ')+StopInd_add
                          StopInd_add=StartInd_add
                        enddo
                        i=80
                      elseif (ChemInpLine_add(i:i+3)=='TROE') then
                      !  if (lroot) print*,ChemInpLine_add(i:i+3),&
                      !      '   coefficients for reaction ', reaction_name(k),&
                      !      'number ', k
                        VarNumber_add=1; StartInd_add=i+5; StopInd_add=i+5
                        do while (VarNumber_add<4)
                          StopInd_add=index(ChemInpLine_add(StartInd_add:),&
                              ' ')+StartInd_add-2
                          StopInd_add_=index(ChemInpLine_add(StartInd_add:),&
                              '/')+StartInd_add-2
                          StopInd_add=min(StopInd_add,StopInd_add_)
                          if (StopInd_add==StartInd_add) then
                            StartInd_add=StartInd_add+1
                          else
                            if (VarNumber_add==1) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') troe_coeff(1,k)
                            elseif (VarNumber_add==2) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') troe_coeff(2,k)
                            elseif (VarNumber_add==3) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') troe_coeff(3,k)
                            else
                              call fatal_error("read_reactions","no such VarNumber: "// &
                                               trim(itoa(VarNumber_add)))
                            endif
                          endif
                          VarNumber_add=VarNumber_add+1
                          !   StartInd_add=StopInd_add
                          StartInd_add=verify(ChemInpLine_add(StopInd_add+1:),&
                              ' ')+StopInd_add
                          StopInd_add=StartInd_add
                        enddo
                        i=80
                      elseif (ChemInpLine_add(i:i+3)=='HIGH') then
                      !  if (lroot) print*,ChemInpLine_add(i:i+3),&
                      !      '   coefficients for reaction ', reaction_name(k),&
                      !      'number ', k
                        VarNumber_add=1; StartInd_add=i+5; StopInd_add=i+5
                        do while (VarNumber_add<4)
                          StopInd_add=index(ChemInpLine_add(StartInd_add:),' ')+StartInd_add-2
                          StopInd_add_=index(ChemInpLine_add(StartInd_add:),'/')+StartInd_add-2
                          StopInd_add=min(StopInd_add,StopInd_add_)
                          if (StopInd_add==StartInd_add) then
                            StartInd_add=StartInd_add+1
                          else
                            if (VarNumber_add==1) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E25.8)') high_coeff(1,k)
                              if (high_coeff(1,k)/=0.) high_coeff(1,k)=log(high_coeff(1,k))
                            elseif (VarNumber_add==2) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') high_coeff(2,k)
                            elseif (VarNumber_add==3) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') high_coeff(3,k)
                            else
                              call fatal_error("read_reactions","no such VarNumber: "// &
                                               trim(itoa(VarNumber_add)))
                            endif
                          endif
                          VarNumber_add=VarNumber_add+1
                          !   StartInd_add=StopInd_add
                          StartInd_add=verify(ChemInpLine_add(StopInd_add+1:),' ')+StopInd_add
                          StopInd_add=StartInd_add
                        enddo
                        i=80
                      else
                        !                a_k4=0.
                        StartInd_add=i; StopInd_add=0; StopInd_add_=0
                        do while (ChemInpLine_add(i:i+1)/='  ')
                          find_specie=.true.
                          do while (StartInd_add/=StopInd_add_)
                            StopInd_add=index(ChemInpLine_add(StartInd_add:),'/')+StartInd_add-2
                            StopInd_add_=index(ChemInpLine_add(StartInd_add:),' ')+StartInd_add-1
                            if (find_specie) then
                              call find_species_index(ChemInpLine_add(StartInd_add:StopInd_add),ind_glob, &
                                  ind_chem,found_specie)
                            else
                              if (found_specie) then
                                read(unit=ChemInpLine_add(StartInd_add:StopInd_add),fmt='(E15.8)') &
                                      a_k4(ind_chem,k)
                              else
                                print*,'ChemInpLine=',ChemInpLine_add
                                print*,'Specie=',ChemInpLine_add(StartInd_add:StopInd_add)
                                print*,'StartInd_add,StopInd_add=',StartInd_add,StopInd_add
                                call fatal_error("read_reactions","did not find specie")
                              endif
                            endif
                            StartInd_add=StopInd_add+2
                            find_specie=.false.
                          enddo
                          i=StopInd_add_
                          StartInd_add=StartInd_add+1
                        enddo
                        i=80
                        !call find_species_index('N2',&
                        !    ind_glob,ind_chem,found_specie)
                        !if (found_specie) a_k4(ind_chem,k)=1.
                        do ind_chem=1,nchemspec
                          if (a_k4(ind_chem,k)==impossible) a_k4(ind_chem,k)=1
                        enddo

                      endif
                    enddo
                    goto 100 ! this is an excellent example of "spaghetti code" [Bourdin.KIS]

                  else
                    found_new_reaction=.true.
                  endif

                else
                  LastLeftCharacter=SeparatorInd-1
                endif

                StartInd=1
                PlusInd=index(ChemInpLine(StartInd:LastLeftCharacter),'+')+StartInd-1
                do while (PlusInd<LastLeftCharacter .AND. PlusInd>0)
                  StopInd=PlusInd-1
                  call build_stoich_matrix(StartInd,StopInd,k,ChemInpLine,.false.)
                  StartInd=StopInd+2
                  plusind_=index(ChemInpLine(StartInd:),'+')
                  if (plusind_ > 0) then
                    PlusInd=plusind_+StartInd-1
                  else
                    PlusInd=10000
                  endif
                enddo
                StopInd=LastLeftCharacter
                call build_stoich_matrix(StartInd,StopInd,k,ChemInpLine,.false.)
!
! Find product side stoichiometric coefficients
!
                StartInd=index(ChemInpLine,'>')+1
                if (StartInd==1) StartInd=index(ChemInpLine,'=')+1
                SeparatorInd=index(ChemInpLine(StartInd:),' ')+StartInd-1

                ParanthesisInd=index(ChemInpLine(StartInd:),'(+M)')+StartInd-1
                MplussInd=index(ChemInpLine(StartInd:),'+M')+StartInd-1

                if (ParanthesisInd>StartInd) then
                  LastLeftCharacter=min(ParanthesisInd,SeparatorInd)-1
                elseif (MplussInd>StartInd) then
                  LastLeftCharacter=min(MplussInd,SeparatorInd)-1
                else
                  LastLeftCharacter=SeparatorInd-1
                endif
                PlusInd=index(ChemInpLine(StartInd:LastLeftCharacter),'+')+StartInd-1
                do while (PlusInd<LastLeftCharacter .AND. PlusInd>StartInd)
                  StopInd=PlusInd-1
                  call build_stoich_matrix(StartInd,StopInd,k,ChemInpLine,.true.)
                  StartInd=StopInd+2
                  PlusInd=index(ChemInpLine(StartInd:),'+')+StartInd-1
                enddo
                StopInd=LastLeftCharacter
                call build_stoich_matrix(StartInd,StopInd,k,ChemInpLine,.true.)
!
! Find Arrhenius coefficients
!
                VarNumber=1; StartInd=1; StopInd =0
                stringloop: do while (VarNumber<4)
                  StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
                  StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
                  StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
                  if (StopInd==StartInd) then
                    StartInd=StartInd+1
                  else
                    if (VarNumber==1) then
                      read (unit=ChemInpLine(StartInd:StopInd),fmt='(E25.8)') B_n(k)
                      if (B_n(k)/=0.) B_n(k)=log(B_n(k))
                    elseif (VarNumber==2) then
                      read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)') alpha_n(k)
                    elseif (VarNumber==3) then
                      read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)') E_an(k)
                    else
                      call fatal_error("read_reactions","no such VarNumber: "//trim(itoa(VarNumber)))
                    endif
                    VarNumber=VarNumber+1
                    StartInd=StopInd
                  endif
                  if (StartInd==80) exit
                enddo stringloop
!
! Increase reaction counter by one
!
                k=k+1
              endif
            endif
          endif
        endif
      enddo dataloop3
1012  continue
      close(file_id)
!
    endsubroutine read_reactions
!***********************************************************************
    subroutine read_reactions_mod
!
!  Read Sijp_mod.dat and Sijm_mod.dat, but so far only Sijp_mod.dat.
!
!  26-apr-24/nils+axel: coded
!
    integer :: lun=23
!
    open(lun,file='Sijp_mod.dat')
    read(lun,*) Sijp_mod
    close(lun)
!
    endsubroutine read_reactions_mod
!***********************************************************************
    subroutine write_matrices
!
!  Print Sijp_mod.dat and others
!
!  26-apr-24/nils+axel: coded
!
    integer :: ireactions
!
    if (lroot) then
      write(*,*) 'Sijp and Sijm based on chem.inp file'
      write(*,*) 'Sijp_mod based on the Sijp_mod.dat file (read only if lchem_detailed=F)'
      write(*,*) 'Sijp_ is the sum of those (used to change power of in reaction equations)'
      write(*,*) 'Sijp_mod ='
      write(*,1000) varname(ichemspec(:))
      do ireactions=1,mreactions
        write(*,1001) ireactions,Sijp_mod(:,ireactions)
      enddo
      write(*,*)
!
      write(*,*) 'Sijp='
      write(*,1000) varname(ichemspec(:))
      do ireactions=1,mreactions
        write(*,1001) ireactions,Sijp(:,ireactions)
      enddo
      write(*,*)
!
      write(*,*) 'Sijp_='
      write(*,1000) varname(ichemspec(:))
      do ireactions=1,mreactions
        write(*,1001) ireactions,Sijp_(:,ireactions)
      enddo
      write(*,*)
!
      write(*,*) 'Sijm_='
      write(*,1000) varname(ichemspec(:))
      do ireactions=1,mreactions
        write(*,1001) ireactions,Sijm_(:,ireactions)
      enddo
      write(*,*)
    endif
!
1000 format(5x,20a6)
1001 format(i2,20f6.1)
    endsubroutine write_matrices
!***********************************************************************
    subroutine get_1step_test_reaction_rate(f,vreact_p,vreact_m,p,reac)
!
!  26-oct-25/TP: carved from get_reaction_rate
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type (pencil_case), intent(in) :: p
      real, dimension(nx,nreactions), intent(out) :: vreact_p, vreact_m
      integer, intent(in) :: reac

      integer :: i

      do i = 1,nx
        !TP: cannot access 1 on GPU
        if (p%TT(i) > Tc) then
          vreact_p(i,reac) = f(l1,m,n,iux)*f(l1,m,n,iux)*p%rho(1)*Cp_const &
              /lambda_const*beta*(beta-1.)*(1.-f(l1-1+i,m,n,ichemspec(ipr)))
        else
          vreact_p(i,reac) = 0.
        endif
      enddo
      vreact_m(:,reac) = 0.
    endsubroutine get_1step_test_reaction_rate
!***********************************************************************
    subroutine get_reaction_rate(f,vreact_p,vreact_m,p)
!
!  This subroutine calculates forward and reverse reaction rates,
!  if chem.inp file exists.
!  For more details see Chemkin Theory Manual
!
!  NILS: This routine is by far the most CPU intencive routine in the code,
!  NILS: so one should maybe put some effort into optimizing it more.
!
!  17-mar-08/natalia: coded
!  11-nov-10/julien: optimized, reaction rates are calculated under a
!                    logarithmic form
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,nreactions), intent(out) :: vreact_p, vreact_m
!
      type (pencil_case) :: p
      real, dimension(nx) :: dSR, dHRT, Kp, Kc
      real, dimension(nx) :: prod1, prod2
      real, dimension(nx) :: kf, kr
      real, dimension(nx) :: rho_cgs, p_atm
      real, dimension(nx) :: mix_conc
      integer :: k, reac, i
      real :: sum_tmp, ddd
      real :: Rcal, Rcal1, lnRgas, l10, lnp_atm
      logical, save :: lwrite_first=.true.
      character(len=fnlen) :: input_file="./data/react.out"
      integer :: file_id=123
      real :: B_n_0, alpha_n_0, E_an_0
      real, dimension(nx) :: kf_0, Pr, sum_sp
      real, dimension(nx) :: Fcent, ccc, nnn, lnPr, FF, tmpF
      real, dimension(nx) :: TT1_loc
!
!  Check which reactions rate method we will use
!
      if (reac_rate_method == 'chemkin') then
!
        TT1_loc = p%TT1
!
        if (lwrite_first)  open (file_id,file=input_file)
!
!  p is in atm units; atm/bar=1./10.13
!  NILS: I think the conversion constant to calories is 4.184 instead 4.14
!
        if (unit_system/="cgs") call fatal_error("get_reaction_rate","Unit system must be cgs.")

        !if (chemkin_units=="SI") then
        !
        !endif

        Rcal = Rgas_unit_sys/4.184*1e-7
        Rcal1 = 1./Rcal
        lnRgas = log(Rgas)
        l10 = log(10.)
        rho_cgs = p%rho*unit_mass/unit_length**3
        lnp_atm = log(1e6*unit_length**3/unit_energy)
        p_atm = 1e6*(unit_length**3)/unit_energy
!
!  calculation of the reaction rate
!
        do reac = 1,nreactions
!
!  Find the product of the species molar consentrations (where
!  each molar consentration is taken to the power of the number)
!
          prod1 = 1.
          prod2 = 1.
          do k = 1,nchemspec
            if (abs(Sijp_(k,reac)) == 1) then
              prod1 = prod1*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)/species_constants(k,imass))
            elseif (abs(Sijp_(k,reac)) == 2) then
              prod1 = prod1*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)/species_constants(k,imass)) &
                           *(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)/species_constants(k,imass))
            elseif (abs(Sijp_(k,reac)) > 0) then
              prod1 = prod1*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)/species_constants(k,imass))**Sijp_(k,reac)
            endif
          enddo
          do k = 1,nchemspec
            if (abs(Sijm_(k,reac)) == 1.0) then
              prod2 = prod2*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)/species_constants(k,imass))
            elseif (abs(Sijm_(k,reac)) == 2.0) then
              prod2 = prod2*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)/species_constants(k,imass)) &
                           *(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)/species_constants(k,imass))
            elseif (abs(Sijm_(k,reac)) > 0.0) then
              prod2 = prod2*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)/species_constants(k,imass))**Sijm_(k,reac)
            endif
          enddo
!
!  Find forward rate constant for reaction 'reac'
!
          if (latmchem) then
            if ((B_n(reac) == 0.) .and. (alpha_n(reac) == 0.) .and. (E_an(reac) == 0.)) then
              do i = 1,nx
                call calc_extra_react(f,reac,kf(i),i,m,n,p)
              enddo
              kf = log(kf)
            else
              kf = B_n(reac)+alpha_n(reac)*p%lnTT-E_an(reac)*TT1_loc
            endif
          else
            kf = B_n(reac)+alpha_n(reac)*p%lnTT-E_an(reac)*Rcal1*TT1_loc
          endif
!
!  Find backward rate constant for reaction 'reac'
!
          dSR = 0.
          dHRT = 0.
          sum_tmp = 0.
          do k = 1,nchemspec
            dSR = dSR+(Sijm(k,reac) -Sijp(k,reac))*p%S0_R(:,k)
            dHRT = dHRT+(Sijm(k,reac)-Sijp(k,reac))*p%H0_RT(:,k)
            sum_tmp = sum_tmp+(Sijm(k,reac)-Sijp(k,reac))
          enddo
          Kp = dSR-dHRT
!
          if (sum_tmp == 0.) then
            Kc = Kp
          else
            Kc = Kp+sum_tmp*(lnp_atm-p%lnTT-lnRgas)
          endif
!
!  Multiply by third body reaction term
!
          if(a_k4_min(reac) < impossible) then
            sum_sp = 0.
            do k = 1,nchemspec
              sum_sp = sum_sp+a_k4(k,reac)*f(l1:l2,m,n,ichemspec(k))  &
                      *rho_cgs(:)/species_constants(k,imass)
            enddo
            mix_conc = sum_sp
          else
            sum_sp = 1.
            mix_conc = rho_cgs(:)*p%mu1(:)/unit_mass
          endif
!
!  The Lindeman approach to the fall of reactions
!
          if (low_coeff_abs_max(reac) > 0.) then
            B_n_0 = low_coeff(1,reac)
            alpha_n_0 = low_coeff(2,reac)
            E_an_0 = low_coeff(3,reac)
            kf_0(:) = B_n_0+alpha_n_0*p%lnTT(:)-E_an_0*Rcal1*TT1_loc(:)
            Pr = exp(kf_0-kf)*mix_conc
            kf = kf+log(Pr/(1.+Pr))
          elseif (high_coeff_abs_max(reac) > 0.) then
            B_n_0 = high_coeff(1,reac)
            alpha_n_0 = high_coeff(2,reac)
            E_an_0 = high_coeff(3,reac)
            kf_0(:) = B_n_0+alpha_n_0*p%lnTT(:)-E_an_0*Rcal1*TT1_loc(:)
            Pr = exp(kf_0-kf)*mix_conc
            kf = kf-log(1.+Pr)
          endif
!
! The Troe approach
!
          if (troe_coeff_abs_max(reac) > 0.) then
            Fcent = (1.-troe_coeff(1,reac))*exp(-p%TT(:)/troe_coeff(2,reac)) &
                    +troe_coeff(1,reac)*exp(-p%TT(:)/troe_coeff(3,reac))
            ccc = -0.4-0.67*log10(Fcent)
            nnn = 0.75-1.27*log10(Fcent)
            ddd = 0.14
            lnPr = log10(Pr)
            tmpF = ((lnPr+ccc)/(nnn-ddd*(lnPr+ccc)))**2
            tmpF = 1./(1.+tmpF)
            FF = tmpF*log10(Fcent)
            FF = FF*l10
            kf = kf+FF
          endif
!
!  Find forward (vreact_p) and backward (vreact_m) rate of
!  progress variable.
!  (vreact_p - vreact_m) is labeled q in the chemkin manual
!
          if (latmchem) then
            kr = kf
          else
            kr = kf-Kc
          endif
          if (lpencil_check_at_work) where (kr > 32) kr = kr / exp(real(nint(alog(kr))))
!
          if (Mplus_case(reac)) then
            where (prod1 > 0.)
              vreact_p(:,reac) = prod1*exp(kf)
            elsewhere
              vreact_p(:,reac) = 0.
            endwhere
            where (prod2 > 0.)
              vreact_m(:,reac) = prod2*exp(kr)
            elsewhere
              vreact_m(:,reac) = 0.
            endwhere
!
          else
            where (prod1 > 0.)
              vreact_p(:,reac) = prod1*exp(kf)*sum_sp
            elsewhere
              vreact_p(:,reac) = 0.
            endwhere
            where (prod2 > 0.)
              vreact_m(:,reac) = prod2*exp(kr)*sum_sp
            elsewhere
              vreact_m(:,reac) = 0.
            endwhere
          endif
!
          if (.not. back(reac)) vreact_m(:,reac) = 0.
        enddo
!
! This part calculates forward and reverse reaction rates
! for the test case R->P
! For more details see Doom, et al., J. Comp. Phys., 226, 2007
!
      elseif (reac_rate_method == '1step_test') then
        call get_1step_test_reaction_rate(f,vreact_p,vreact_m,p,reac)
!
!  Add alternative method for finding the reaction rates based on work by
!  Roux et al. (2009)
!  NILS: Should split the different methods for calculating the reaction
!  NILS: rates into different subroutines at some point.
!
      elseif (reac_rate_method == 'roux') then
        call roux(f,p,vreact_p,vreact_m)
      endif
!
      if (lwrite_first .and. lroot) print*,'get_reaction_rate: writing react.out file'
      lwrite_first = .false.
!
    endsubroutine get_reaction_rate
!***********************************************************************
    subroutine roux(f,p,vreact_p,vreact_m)
!
!  nilshau: 2011.01.11 (coded)
!
!  Calculate reaction rates based on the method of Roux et al. (2009).
!  This is a single step reaction mechanism for propane:
!
!  C3H8 + 5O2 +XN2 -> 3CO2 + 4H2O + XN2
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,nreactions), intent(out) :: vreact_p, vreact_m
      type (pencil_case) :: p
!
      real :: Rcal, f_phi, E_a
      real :: init_C3H8, init_O2
      integer :: i_O2, i_C3H8, j
      real, dimension(nx) :: activation_energy, pre_exp, term1, term2
!
      if (nreactions /= 1) call fatal_error('roux','nreactions should always be 1')
!
!  Check that a global equivalence ratio is given at input
!
      if (global_phi == impossible) call fatal_error('roux','global_phi must be given as input')
!
      Rcal = Rgas_unit_sys/4.14*1e-7
!
!  Check that oxygen and propane exist and find their molar masses
!
      i_O2   = i_O2_glob
      i_C3H8 = i_C3H8_glob
      if (.not. lO2) then
        call fatal_error('roux','O2 is not defined')
      endif
      if (.not. lC3H8) then
        call fatal_error('roux','C3H8 is not defined')
      endif
!
!  Print debugging output
!
      if (headtt) then
        init_O2 = initial_massfractions(ichem_O2)
        init_C3H8 = initial_massfractions(ichem_C3H8)
        print*,'i_O2, i_C3H8, ichem_O2, ichem_C3H8=', i_O2, i_C3H8, ichem_O2, ichem_C3H8
        print*,'lO2, lC3H8=',lO2, lC3H8
        print*,'init_C3H8,init_O2,mO2,mC3H8=',init_C3H8,init_O2,mO2,mC3H8
      endif
!
!  Find Laminar flame speed corrector based on equivalence ratio phi
!
      f_phi = 0.5*(1+tanh((0.8-global_phi)/1.5)) &
             +2.11/4*(1+tanh((global_phi-0.11)/0.2)) &
             *(1+tanh((1.355-global_phi)/0.24))
!
!  Find the classical Arrhenius terms
!
      E_a = 31126
      activation_energy = exp(-E_a*p%TT1/Rcal)
      pre_exp = 3.2916e10
!
!  Find density and mass fraction dependent terms
!
      term1 = (f(l1:l2,m,n,i_C3H8)*p%rho/species_constants(i_C3H8,imass))**0.856
      term2 = (f(l1:l2,m,n,i_O2)*p%rho/species_constants(i_O2,imass))**0.503
!
!  Use the above to find reaction terms
!
      vreact_p(:,1) = f_phi*pre_exp*term1*term2*activation_energy
!
!  Set reaction rate to zero when mass fractions of propane of oxygen is
!  very close to zero.
!
      where (f(l1:l2,m,n,i_C3H8) < 1e-12) vreact_p(:,1) = 0.
      where (f(l1:l2,m,n,i_O2) < 1e-12) vreact_p(:,1) = 0.
      vreact_m(:,1) = 0.
!
    endsubroutine roux
!***********************************************************************
    subroutine calc_reaction_term(f,p)
!
!  Calculation of the reaction term
!
      real :: alpha, eps
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,mreactions) :: vreactions
      real, dimension(nx,nchemspec) :: xdot
      real, dimension(nx) :: rho1
      real, dimension(nx,nchemspec) :: molm
      type (pencil_case) :: p
      integer :: k,j,ii
      integer,parameter :: i1=1, i2=2, i3=3, i4=4, i5=5, i6=6, i7=7, i8=8, i9=9, i10=10
      integer,parameter :: i11=11, i12=12, i13=13, i14=14, i15=15, i16=16, i17=17, i18=18, i19=19
!
      eps = sqrt(epsilon(alpha))
!
      p%DYDt_reac = 0.
      rho1 = 1./p%rho
      if (lcheminp .and. (.not. l1step_test)) then
        do k = 1,nchemspec
          molm(:,k) = rho1*species_constants(k,imass)
        enddo
      else
        molm = 1.
      endif
!
!  if we do reactions, we must calculate the reaction speed vector
!  outside the loop where we multiply it by the stoichiometric matrix
!
      if (.not. lcheminp) then
!
!  Axel' case
!
        do j = 1,nreactions
          if (lkreactions_alpha) then
            alpha = kreactions_alpha(j)
          else
            alpha = 1.
          endif
          vreactions_p(:,j) = alpha*kreactions_p(j)*kreactions_z(n,j)
          vreactions_m(:,j) = alpha*kreactions_m(j)*kreactions_z(n,j)
          do k = 1,nchemspec
            vreactions_p(:,j) = vreactions_p(:,j)*f(l1:l2,m,n,ichemspec(k))**Sijm(k,j)
            vreactions_m(:,j) = vreactions_m(:,j)*f(l1:l2,m,n,ichemspec(k))**Sijp(k,j)
          enddo
        enddo
        vreactions_m = -vreactions_m
        vreactions_p = -vreactions_p
      else
!
!  Chemkin data case
!
        call get_reaction_rate(f,vreactions_p,vreactions_m,p)
      endif
!
!  Calculate rate of reactions (labeled q in the chemkin manual)
!
      vreactions = vreactions_p-vreactions_m
!
!  Calculate production rate for all species k (called \dot(\omega)_k
!  in the chemkin manual)
!
      xdot = 0.
      do k = 1,nchemspec
        do j = 1,nreactions
          xdot(:,k) = xdot(:,k)-stoichio(k,j)*vreactions(:,j)*molm(:,k)
        enddo
      enddo
      p%DYDt_reac = xdot*unit_time


!
! NH:
!
!  Sums for diagnostics
!
      !sum_omega=0.
      !sum_Y=0.
      !do k=1,nchemspec
      !  sum_omega=sum_omega+maxval(p%DYDt_reac(:,k))
      !  sum_Y=sum_Y+maxval(f(l1:l2,m,n,ichemspec(k)))
      !enddo
!
    endsubroutine calc_reaction_term
!***********************************************************************
    subroutine  write_net_reaction
!
!  write net reactions to file
!
      if (lchemistry_diag) then
        open (1,file=trim(datadir)//'/net_reactions.dat',position='append')
        write (1,*) t
        write (1,'(8e10.2)') net_react_p, net_react_m
        close (1)
!
! Reset to zero for next time
!
        net_react_m = 0.
        net_react_p = 0.
      endif
!
    endsubroutine  write_net_reaction
!***********************************************************************
    subroutine  calc_collision_integral(omega,lnTst,Omega_kl)
!
!  Get coefficients for calculating of the collision integral
!  This routine is called from calc_diff_visc_coeff, which again is called from
!  calc_for_chem_mixture, which is why we work on full chunks of arrays here.
!
!  03-apr-08/natalia: coded
!
      character(len=*), intent(in) :: omega
      real, dimension(mx,my,mz), intent(in) :: lnTst
      real, dimension(mx,my,mz), intent(out) :: Omega_kl

      integer :: i
      real, dimension(8) :: aa
!
      select case (omega)
      case ('Omega11')
        aa(1) = 6.96945701E-1
        aa(2) = 3.39628861E-1
        aa(3) = 1.32575555E-2
        aa(4) = -3.41509659E-2
        aa(5) = 7.71359429E-3
        aa(6) = 6.16106168E-4
        aa(7) = -3.27101257E-4
        aa(8) = 2.51567029E-5
      case ('Omega22')
        aa(1) = 6.33225679E-1
        aa(2) = 3.14473541E-1
        aa(3) = 1.78229325E-2
        aa(4) = -3.99489493E-2
        aa(5) = 8.98483088E-3
        aa(6) = 7.00167217E-4
        aa(7) = -3.82733808E-4
        aa(8) = 2.97208112E-5
      case default
        call fatal_error('calc_collision_integral','no such omega: '//trim(omega))
      endselect
!
      Omega_kl = 0.
      do i = 1,8
        Omega_kl = Omega_kl+aa(i)*(lnTst)**(i-1)
      enddo
      Omega_kl = 1./Omega_kl
!
    endsubroutine  calc_collision_integral
!***********************************************************************
    subroutine calc_diff_visc_coef(f)
!
!  Calculation of the binary diffusion coefficients and the species viscosities.
!  This routind is called from calc_for_chem_mixture,
!  which is why we work on full chunks of arrays here.
!
      real, dimension(mx,my,mz,mfarray) :: f
      intent(in) :: f
      real, dimension(mx,my,mz) :: Omega_kl, prefactor
      real, dimension(mx,my,mz) :: lnTjk, lnTk_array
      integer :: k, j, j2, j3
      real :: eps_jk, sigma_jk, m_jk, delta_jk, delta_st
      real :: Na, tmp_local, tmp_local2, delta_jk_star
      character(len=7) :: omega
!
!  Find binary diffusion coefficients
!
      Na=N_avogadro_cgs
      tmp_local = 3./16.*sqrt(2.*k_B_cgs**3/pi)
      prefactor(ll1:ll2,mm1:mm2,nn1:nn2) = tmp_local*sqrt(TT_full(ll1:ll2,mm1:mm2,nn1:nn2)) &
          *unit_length**3/(Rgas_unit_sys*rho_full(ll1:ll2,mm1:mm2,nn1:nn2))
!
      omega = 'Omega11'
!
! Check if we use fixed Schmidt number for speeding up calculations
!
      if (.not. lfix_Sc) then
!
! Check if we use simplified version of the binary diffusion calculation
!
        if (ldiffusion .and. (.not. lDiff_simple) .and. (.not. lDiff_lewis)) then
!
!  Do non-simplified binary diffusion coefficient
!
          do k = 1,nchemspec
            do j = k,nchemspec
!  Account for the difference between eq. 5-4 and 5-31 in the Chemkin theory
!  manual
!
              if (j /= k) then
                eps_jk = sqrt(tran_data(j,2)*tran_data(k,2))
                sigma_jk = 0.5*(tran_data(j,3)+tran_data(k,3))*1e-8
                m_jk = (species_constants(j,imass)*species_constants(k,imass)) &
                    /(species_constants(j,imass)+species_constants(k,imass))/Na
                delta_jk = 0.5*tran_data(j,4)*tran_data(k,4)*1e-18*1e-18
              else
                eps_jk = tran_data(j,2)
                sigma_jk = tran_data(j,3)*1e-8
                m_jk = species_constants(j,imass)/(2*Na)
                delta_jk = 0.5*(tran_data(j,4)*1e-18)*(tran_data(j,4)*1e-18)
              endif
!
! Loop over all grid points
!
              do j3 = nn1,nn2
                do j2 = mm1,mm2
                  if (ltemperature_nolog) then
                    lnTjk(ll1:ll2,j2,j3) = log(f(ll1:ll2,j2,j3,ilnTT)/eps_jk)
                  else
                    lnTjk(ll1:ll2,j2,j3) = f(ll1:ll2,j2,j3,ilnTT)-log(eps_jk)
                  endif
!
                  Omega_kl(ll1:ll2,j2,j3) = &
                      1./(6.96945701E-1   +3.39628861E-1*lnTjk(ll1:ll2,j2,j3) &
                      +1.32575555E-2*lnTjk(ll1:ll2,j2,j3)*lnTjk(ll1:ll2,j2,j3) &
                      -3.41509659E-2*lnTjk(ll1:ll2,j2,j3)**3 &
                      +7.71359429E-3*lnTjk(ll1:ll2,j2,j3)**4 &
                      +6.16106168E-4*lnTjk(ll1:ll2,j2,j3)**5 &
                      -3.27101257E-4*lnTjk(ll1:ll2,j2,j3)**6 &
                      +2.51567029E-5*lnTjk(ll1:ll2,j2,j3)**7)
                  delta_jk_star = delta_jk/(eps_jk*k_B_cgs*sigma_jk**3)
!
                  Omega_kl(ll1:ll2,j2,j3) = Omega_kl(ll1:ll2,j2,j3) &
                      +0.19*delta_jk_star*delta_jk_star/(TT_full(ll1:ll2,j2,j3)/eps_jk)
                  if (j /= k) then
                    Bin_Diff_coef(ll1:ll2,j2,j3,k,j) = prefactor(ll1:ll2,j2,j3)/mu1_full(ll1:ll2,j2,j3) &
                        /(sqrt(m_jk)*sigma_jk*sigma_jk*Omega_kl(ll1:ll2,j2,j3))
                  else
                    Bin_Diff_coef(ll1:ll2,j2,j3,k,j) = prefactor(ll1:ll2,j2,j3) &
                        /(sqrt(m_jk)*sigma_jk*sigma_jk*Omega_kl(ll1:ll2,j2,j3))*species_constants(k,imass)
!
                  endif
                enddo
              enddo
            enddo
          enddo
!
          do k = 1,nchemspec
            do j = 1,k-1
              Bin_Diff_coef(ll1:ll2,mm1:mm2,nn1:nn2,k,j) = Bin_Diff_coef(ll1:ll2,mm1:mm2,nn1:nn2,j,k)
            enddo
          enddo
!
        else
          if (.not. lDiff_simple .and. .not. lDiff_lewis) Bin_Diff_coef = 0.
        endif
      endif
!
!  Calculate viscosity
!
!     if (visc_const==impossible) then
      omega = 'Omega22'
      tmp_local = 5./16.*sqrt(k_B_cgs/(Na*pi))
!
      do k = 1,nchemspec
        tmp_local2 = sqrt(species_constants(k,imass))/  &
            ((tran_data(k,3)*1e-8)*(tran_data(k,3)*1e-8))*tmp_local
!
! 1 Debye = 10**(-18) esu -> (1e-18*tran_data(k,4))
!
        delta_st = (1e-18*tran_data(k,4))*(1e-18*tran_data(k,4))/2./ &
            (tran_data(k,2)*k_B_cgs*(tran_data(k,3)*1e-8)**3)
!
        if (ltemperature_nolog) then
          lnTk_array = log(f(:,:,:,ilnTT)/tran_data(k,2))
        else
          lnTk_array = f(:,:,:,ilnTT)-log(tran_data(k,2))
        endif
        call calc_collision_integral(omega,lnTk_array,Omega_kl)
!
        species_viscosity(ll1:ll2,mm1:mm2,nn1:nn2,k) = &
             sqrt(TT_full(ll1:ll2,mm1:mm2,nn1:nn2))/(Omega_kl(ll1:ll2,mm1:mm2,nn1:nn2) &
            +0.2*delta_st*delta_st/(TT_full(ll1:ll2,mm1:mm2,nn1:nn2)/tran_data(k,2)))*tmp_local2 &
            /(unit_mass/unit_length/unit_time)
      enddo
!      endif
!
    endsubroutine calc_diff_visc_coef
!***********************************************************************
    subroutine calc_therm_diffus_coef
!
!  Calculate the thermal diffusion coefficient based on equation 5-17 in
!  the Chemkin theory manual
!
      real, dimension(mx,my,mz,nchemspec) :: species_cond
      real, dimension(mmx) :: tmp_val, ZZ, FF, tmp_sum, tmp_sum2
      real, dimension(mmx) :: AA, BB, f_tran, f_rot, f_vib
      real, dimension(mmx) :: Cv_vib_R, T_st, pi_1_5, pi_2
      real :: Cv_rot_R, Cv_tran_R
      integer :: j2,j3,k
!
!        call timing('calc_therm_diffus_coef','just entered')
!
      pi_2 = pi*pi
      pi_1_5 = pi*sqrt(pi)
!
!  With lspecies_cond_simplified, species conduction is calculated
!  according to the book "Transport phenomena" by Bird, Warren, & Lightfoot
!  page 276, Eq.(9.3-15).
!
      do j3 = nn1,nn2
        do j2 = mm1,mm2
          tmp_sum = 0.
          tmp_sum2 = 0.
          do k = 1,nchemspec
            if (lspecies_cond_simplified) then
              species_cond(ll1:ll2,j2,j3,k) = (species_viscosity(ll1:ll2,j2,j3,k))*&
                ( cp_spec_glo(ll1:ll2,j2,j3,k) + 1.25*Rgas/species_constants(k,imass) )
            else
!
! Check if the molecule is a single atom (0), linear (1) or non-linear (2).
!
              if (tran_data(k,1) == 0.) then
                Cv_tran_R = 1.5
                Cv_rot_R = 0.
                Cv_vib_R = 0.
              elseif (tran_data(k,1) == 1.) then
                Cv_tran_R = 1.5
                Cv_rot_R = 1.
                Cv_vib_R = cv_R_spec_full(ll1:ll2,j2,j3,k)-2.5
              elseif (tran_data(k,1) == 2.) then
                Cv_tran_R = 1.5
                Cv_rot_R = 1.5
                Cv_vib_R = cv_R_spec_full(ll1:ll2,j2,j3,k)-3.
              else
                Cv_tran_R = 0
                Cv_rot_R = 0
                Cv_vib_R = 0
                call fatal_error('calc_therm_diffus_coef','No such tran_data')
              endif
!
! The rotational and vibrational contributions are zero for the single
! atom molecules but not for the linear or non-linear molecules
!
              if (tran_data(k,1) > 0. .and. (.not. lfix_Sc)) then
                tmp_val = Bin_Diff_coef(ll1:ll2,j2,j3,k,k)*rho_full(ll1:ll2,j2,j3) &
                    /species_viscosity(ll1:ll2,j2,j3,k)
                AA = 2.5-tmp_val
                T_st = tran_data(k,2)/298.
                FF = 1.+pi_1_5/2.*sqrt(T_st)+(pi_2/4.+2.) &
                    *(T_st)+pi_1_5*(T_st)**1.5
                ZZ = tran_data(k,6)*FF
                T_st = tran_data(k,2)/TT_full(ll1:ll2,j2,j3)
                FF = 1.+pi_1_5/2.*sqrt(T_st)+(pi_2/4.+2.) &
                    *(T_st)+pi_1_5*(T_st)**1.5
                ZZ = ZZ/FF
                BB = ZZ+2./pi*(5./3.*Cv_rot_R+tmp_val)
                f_tran = 2.5*(1.- 2./pi*Cv_rot_R/Cv_tran_R*AA/BB)
                f_rot = tmp_val*(1+2./pi*AA/BB)
                f_vib = tmp_val
              else
                f_tran = 2.5
                f_rot = 0.0
                f_vib = 0.0
              endif
              species_cond(ll1:ll2,j2,j3,k) = (species_viscosity(ll1:ll2,j2,j3,k)) &
                  /(species_constants(k,imass)/unit_mass)*Rgas* &
                  (f_tran*Cv_tran_R+f_rot*Cv_rot_R  &
                  +f_vib*Cv_vib_R)
            endif

!
! tmp_sum and tmp_sum2 are used later to find the mixture averaged
! conductivity.
!
            tmp_sum = tmp_sum +XX_full(ll1:ll2,j2,j3,k)*species_cond(ll1:ll2,j2,j3,k)
            tmp_sum2 = tmp_sum2 +XX_full(ll1:ll2,j2,j3,k)/species_cond(ll1:ll2,j2,j3,k)
          enddo
!
! Find the mixture averaged conductivity
!
          where (tmp_sum2 <= 0.)
            lambda_full(ll1:ll2,j2,j3) = 0.
          elsewhere
            lambda_full(ll1:ll2,j2,j3) = 0.5*(tmp_sum+1./tmp_sum2)
          endwhere
        enddo
      enddo
!
      call timing('calc_therm_diffus_coef','just finished')
!
    endsubroutine calc_therm_diffus_coef
!***********************************************************************
    subroutine calc_diffusion_term(f,p)
!
!  Calculate diffusion term, p%DYDt_diff
!
!  22-jun-10/julien: Evaluation of mass diffusion fluxes using Fick's
!                    law for simplified diffusion using constant Lewis numbers
!
      use Sub, only: del2,grad,dot_mn
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p

      real, dimension(nx) :: Xk_Yk
      real, dimension(nx,3) :: gXk_Yk
      real, dimension(nx) :: del2chemspec
      real, dimension(nx) :: diff_op, diff_op1, diff_op2, diff_op3, del2XX, del2lnpp
      real, dimension(nx) :: glnpp_gXkYk, glnrho_glnpp, gD_glnpp, glnpp_glnpp
      real, dimension(nx) :: sum_gdiff, gY_sumdiff, glnmu_glnpp
      real, dimension(nx,3) :: sum_diff, dk_D
      real :: diff_k
      integer :: k,i
!
      intent(in) :: f
!
      p%DYDt_diff = 0.
      diff_k = chem_diff
!
!  Loop over all chemical species.
!
      do k = 1,nchemspec
!
!  Simple chemistry (not using ChemKin formalism).
!  Eliminate need for p%glnrho if density is not included.
!
        if (.not. lcheminp) then
          if (chem_diff /= 0.) then
            diff_k = chem_diff*chem_diff_prefactor(k)
            if (headtt) print*,'dchemistry_dt: k,diff_k=',k,diff_k
            call del2(f,ichemspec(k),del2chemspec)
            !call grad(f,ichemspec(k),gchemspec)
            if (ldensity) then
              call dot_mn(p%glnrho,p%gYYk(:,:,k),diff_op)
              diff_op = diff_op+del2chemspec
            else
              diff_op = del2chemspec
            endif
            p%DYDt_diff(:,k) = diff_k*diff_op
          endif
        else
!
!  Detailed chemistry and transport using CHEMKIN formalism.
!
          if (ldiffusion) then
!
!  Calculate the terms needed by the diffusion fluxes in 3 cases:
!    1) Fickian diffusion law (gradient of species MASS fractions)
!    2) Simplified fluxes (gradient of species MOLAR fractions)
!    3) Detailed transport (with pressure gradients included)
!
            if (lDiff_fick) then
              call del2(f,ichemspec(k),del2chemspec)
!              call grad(f,ichemspec(k),gchemspec)
              call dot_mn(p%glnrho,p%gYYk(:,:,k),diff_op1)
              call dot_mn(p%gdiffk(:,:,k),p%gYYk(:,:,k),diff_op2)
            elseif (lFlux_simple) then
              !call del2(XX_full(:,:,:,k),del2XX)
              call dot_mn(p%glnrho,p%gXXk(:,:,k),diff_op1)
              call dot_mn(p%gdiffk(:,:,k),p%gXXk(:,:,k),diff_op2)
              call dot_mn(p%glnmu,p%gXXk(:,:,k),diff_op3)
            else
              !call del2(XX_full(:,:,:,k),del2XX)
              call dot_mn(p%glnrho,p%gXXk(:,:,k),diff_op1)
              call dot_mn(p%gdiffk(:,:,k),p%gXXk(:,:,k),diff_op2)
              call dot_mn(p%glnmu,p%gXXk(:,:,k),diff_op3)
              call dot_mn(p%glnpp,p%glnpp,glnpp_glnpp)
              do i = 1,3
                gXk_Yk(:,i) = p%gXXk(:,i,k)-p%gYYk(:,i,k)
              enddo
              del2lnpp = p%del2pp/p%pp-glnpp_glnpp
              Xk_Yk = XX_full(l1:l2,m,n,k)-f(l1:l2,m,n,ichemspec(k))
              call dot_mn(p%glnrho,p%glnpp,glnrho_glnpp)
              call dot_mn(p%gdiffk(:,:,k),p%glnpp,gD_glnpp)
              call dot_mn(gXk_Yk,p%glnpp,glnpp_gXkYk)
              call dot_mn(p%glnmu,p%glnpp,glnmu_glnpp)
            endif
!
!  Calculate the diffusion fluxes and dk_D in 3 cases:
!    1) Fickian diffusion law (gradient of species MASS fractions)
!    2) Simplified diffusion fluxes (only gradient of species MOLAR fractions)
!    3) Detailed transport (with pressure gradients included)
!  Note that the ratio Wk/Wm is introduced here and not during the calculation
!  of species diffusion coefficients as before. It indeed depends on the diffusion
!  flux formulation and not on the diffusive properties of each species.
!
            if (lDiff_fick) then
              p%DYDt_diff(:,k) = p%Diff_penc_add(:,k)*(del2chemspec+diff_op1) + diff_op2
              do i = 1,3
                dk_D(:,i) = p%Diff_penc_add(:,k)*p%gYYk(:,i,k)
              enddo
            elseif (lFlux_simple) then
              p%DYDt_diff(:,k) = p%Diff_penc_add(:,k)*p%mukmu1(:,k) &
                  *(p%g2XXk(:,k)+diff_op1-diff_op3) + p%mukmu1(:,k)*diff_op2
              do i = 1,3
                dk_D(:,i) = p%Diff_penc_add(:,k)*p%mukmu1(:,k)*p%gXXk(:,i,k)
              enddo
            else
              p%DYDt_diff(:,k) = p%Diff_penc_add(:,k)*p%mukmu1(:,k) &
                   *(p%g2XXk(:,k)+diff_op1-diff_op3)+ &
                   p%mukmu1(:,k)*diff_op2
              !
              ! Include terms due to pressure gradient
              !
              if (lgradP_terms) then
                p%DYDt_diff(:,k) = p%DYDt_diff(:,k) + &
                     p%Diff_penc_add(:,k)*p%mukmu1(:,k)*Xk_Yk(:) &
                     *(del2lnpp+glnrho_glnpp-glnmu_glnpp)+ &
                     Xk_Yk(:)*p%mukmu1(:,k)*gD_glnpp+p%Diff_penc_add(:,k) &
                     *p%mukmu1(:,k)*glnpp_gXkYk
              endif
              do i = 1,3
                dk_D(:,i) = p%Diff_penc_add(:,k)*p%mukmu1(:,k)*(p%gXXk(:,i,k) + Xk_Yk(:)*p%glnpp(:,i))
              enddo
            endif
!
            if (ldiff_corr) then
              sum_diff = sum_diff+dk_D
              sum_gdiff(:) = sum_gdiff(:)+ p%DYDt_diff(:,k)
            endif
          endif
        endif
      enddo
!
!  Adding correction diffusion velocity to ensure mass balance
!
      if (ldiffusion .and. ldiff_corr) then
        do k = 1,nchemspec
          !call grad(f,ichemspec(k),gchemspec)
          call dot_mn(p%gYYk(:,:,k),sum_diff,gY_sumdiff)
          p%DYDt_diff(:,k) = p%DYDt_diff(:,k) - (gY_sumdiff+f(l1:l2,m,n,ichemspec(k))*sum_gdiff(:))
        enddo
      endif
!
    endsubroutine calc_diffusion_term
!***********************************************************************
    subroutine calc_heatcond_chemistry(f,df,p)
!
!  Calculate gamma*chi*(del2lnT+gradlnTT.grad(lnT+lnrho+lncp+lnchi))
!
!  29-feb-08/natalia: coded
!
      use Sub, only: dot, del2
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: g2TT, g2TTlambda=0., tmp1, del2TT
!
      if (ltemperature_nolog) then
        call dot(p%gTT,p%glambda,g2TTlambda)
        call del2(f,iTT,del2TT)
      else
        call dot(p%glnTT,p%glambda,g2TTlambda)
        call dot(p%glnTT,p%glnTT,g2TT)
      endif
!
!  Add heat conduction to RHS of temperature equation
!
!      if (l1step_test .or. lSmag_heat_transport) then
      if (l1step_test) then
        tmp1 = p%lambda(:)*(p%del2lnTT+g2TT)*p%cv1/p%rho(:)
      else
        if (ltemperature_nolog) then
          tmp1 = (p%lambda(:)*del2TT+g2TTlambda)*p%cv1/p%rho(:)
   !       tmp1 = (p%lambda(:)*p%del2lnTT+g2TTlambda)*p%cv1/p%rho(:)
          df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + tmp1
        else
          tmp1 = (p%lambda(:)*(p%del2lnTT+g2TT)+g2TTlambda)*p%cv1/p%rho(:)
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + tmp1
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_heatcond_chemistry
!***********************************************************************
    subroutine get_RHS_Y_full(RHS_Y)
!
      real, dimension(mx,my,mz,nchemspec) :: RHS_Y
      intent(out) :: RHS_Y
!
      RHS_Y = RHS_Y_full
!
    endsubroutine get_RHS_Y_full
!***********************************************************************
    subroutine get_cs2_full(cs2_full)
!
      real, dimension(mx,my,mz) :: cs2_full
      intent(out) :: cs2_full
      integer :: j,k
!
      do j = 1,my
        do k = 1,mz
          if (minval(cv_full(:,j,k)) <= 0) then
            cs2_full(:,j,k) = 0.
          else
            cs2_full(:,j,k) = cp_full(:,j,k)/cv_full(:,j,k)*mu1_full(:,j,k)*TT_full(:,j,k)*Rgas
          endif
        enddo
      enddo
!
    endsubroutine get_cs2_full
!***********************************************************************
    subroutine get_cs2_slice(f,slice,dir,index)
!
! Find a slice of the speed of sound
!
! 10-dez-09/nils: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(:,:), intent(out) :: slice
      integer, intent(in) :: index, dir
!
      if (dir == 1) then
        slice = cp_full(index,m1:m2,n1:n2)/cv_full(index,m1:m2,n1:n2) &
            *mu1_full(index,m1:m2,n1:n2)*TT_full(index,m1:m2,n1:n2)*Rgas
      elseif (dir == 2) then
        slice = cp_full(l1:l2,index,n1:n2)/cv_full(l1:l2,index,n1:n2) &
            *mu1_full(l1:l2,index,n1:n2)*TT_full(l1:l2,index,n1:n2)*Rgas
      elseif (dir == 3) then
        slice = cp_full(l1:l2,m1:m2,index)/cv_full(l1:l2,m1:m2,index) &
            *mu1_full(l1:l2,m1:m2,index)*TT_full(l1:l2,m1:m2,index)*Rgas
      else
        call fatal_error('get_cs2_slice','No such dir!')
      endif
!
    endsubroutine get_cs2_slice
!***********************************************************************
    subroutine get_gamma_full(gamma_full)
!
      real, dimension(mx,my,mz) :: gamma_full
      intent(out) :: gamma_full
      integer :: j,k
!
      do j = 1,my
        do k = 1,mz
          if (minval(cv_full(:,j,k)) <= 0) then
            gamma_full(:,j,k) = 0.
          else
            gamma_full(:,j,k) = cp_full(:,j,k)/cv_full(:,j,k)
          endif
        enddo
      enddo
!
    endsubroutine get_gamma_full
!***********************************************************************
    subroutine get_gamma_slice(f,slice,dir,index)
!
!  Get a 2D slice of gamma
!
!  10-dez-09/Nils Erland L. Haugen: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(:,:), intent(out) :: slice
      integer, intent(in) :: index, dir
!
      if (dir == 1) then
        slice = cp_full(index,m1:m2,n1:n2)/cv_full(index,m1:m2,n1:n2)
      elseif (dir == 2) then
        slice = cp_full(l1:l2,index,n1:n2)/cv_full(l1:l2,index,n1:n2)
      elseif (dir == 3) then
        slice = cp_full(l1:l2,m1:m2,index)/cv_full(l1:l2,m1:m2,index)
      else
        call fatal_error('get_gamma_slice','No such dir!')
      endif
!
    endsubroutine get_gamma_slice
!***********************************************************************
    subroutine air_field(f,PP)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz) :: sum_Y, tmp
      real :: PP ! (in dynes = 1atm)
!
      logical :: emptyfile=.true.
      logical :: found_specie
      integer :: file_id=123, ind_glob, ind_chem
      character(len=80) :: ChemInpLine,spe
      character(len=10) :: specie_string
      character(len=1) :: tmp_string
      integer :: i, j, k=1
      real :: YY_k, air_mass, TT=300.
      real :: velx=0.
      real, dimension(nchemspec) :: stor2=0.0
      integer, dimension(nchemspec) :: stor1
!
      integer :: StartInd, StopInd, StartInd_1, StopInd_1
      integer :: iostat
      logical :: airdat=.false.,airin=.false.
!
      air_mass = 0.
      StartInd_1 = 1
      StopInd_1 = 0
!
      inquire (file='air.dat',exist=airdat)
      inquire (file='air.in',exist=airin)
      if (airdat .and. airin) &
        call fatal_error('chemistry','both air.in and air.dat found. Please decide for one')

      if (airdat) open(file_id,file='air.dat')
      if (airin) open(file_id,file='air.in')
!
      if (lroot) print*, 'the following parameters and '//&
          'species are found in air.dat (volume fraction fraction in %): '

      dataloop: do

        read(file_id,'(80A)',IOSTAT=iostat) ChemInpLine
        if (iostat < 0) exit dataloop
        emptyFile=.false.
        StartInd_1=1; StopInd_1=0
        StopInd_1=index(ChemInpLine,' ')
        specie_string=trim(ChemInpLine(1:StopInd_1-1))
        tmp_string=trim(ChemInpLine(1:1))

        if (tmp_string == '!' .or. tmp_string == ' ') then
        elseif (tmp_string == 'T') then
          StartInd=1; StopInd =0

          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
          StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1

          read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)') TT
          if (lroot) print*, ' Temperature, K   ', TT

        elseif (tmp_string == 'P') then

          StartInd=1; StopInd =0

          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
          StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1

          read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)') PP
          if (lroot) print*, ' Pressure, Ba   ', PP

        elseif (tmp_string == 'V') then

          StartInd=1; StopInd =0

          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
          StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1

          read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)') velx
          if (lroot) print*, ' Velocity, cm/s   ', velx

        else

          call find_species_index(specie_string,ind_glob,ind_chem,found_specie)

          if (found_specie) then

            StartInd=1; StopInd =0

            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
            StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
            read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)') YY_k
            if (lroot) print*, ' volume fraction, %,    ', YY_k,species_constants(ind_chem,imass)
            if (species_constants(ind_chem,imass)>0.) then
              air_mass=air_mass+YY_k*0.01/species_constants(ind_chem,imass)
              print*,"ind_chem,YY_k,species_constants(ind_chem,imass),air_mass=",&
                   ind_chem,YY_k,species_constants(ind_chem,imass),air_mass
            endif

            if (StartInd==80) exit

            stor1(k)=ind_chem
            stor2(k)=YY_k
            k=k+1
          endif

        endif
      enddo dataloop
!
! Stop if air.dat is empty
!
      if (emptyFile) call fatal_error("air_field",'Input file air.dat was empty')

      if ((sum(stor2) .lt. 99.) .or. (sum(stor2) .gt. 101.)) then
        print*,"sum(stor2)=",sum(stor2)
        call fatal_error("air_field",'The mass fractions should sum to 100.')
      endif
      air_mass=1./air_mass

      do j=1,k-1
        f(:,:,:,ichemspec(stor1(j)))=stor2(j)*0.01
      enddo

      sum_Y=0.

      do j=1,nchemspec
        sum_Y=sum_Y+f(:,:,:,ichemspec(j))
      enddo
      do j=1,nchemspec
        f(:,:,:,ichemspec(j))=f(:,:,:,ichemspec(j))/sum_Y
        initial_massfractions(j)=f(l1,m1,n1,ichemspec(j))
      enddo

      if (mvar < 5) call fatal_error("air_field","can only set existing fields")

      if (.not. reinitialize_chemistry) then
        if (.not.lflame_front .and. .not.ltriple_flame .and. .not.lFlameMaster)  then
          if (ltemperature_nolog) then
            f(:,:,:,iTT)=TT
          else
            f(:,:,:,ilnTT)=alog(TT)!+f(:,:,:,ilnTT)
          endif
          if (ldensity_nolog) then
            f(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)*air_mass/TT)/unit_mass*unit_length**3
          else
            tmp=(PP/(k_B_cgs/m_u_cgs)*air_mass/TT)/unit_mass*unit_length**3
            f(:,:,:,ilnrho)=alog(tmp)
          endif
          if (nxgrid>1) f(:,:,:,iux)=f(:,:,:,iux)+init_ux
        endif
     endif

      if (init_from_file) then
        if (lroot) print*, 'Velocity field read from file, initialization' // &
                           'of density and temperature'
        if (ldensity_nolog) then
          f(:,:,:,ilnrho)=f(:,:,:,ilnrho)*(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
        else
          tmp=f(:,:,:,ilnrho)*(PP/(k_B_cgs/m_u_cgs)*air_mass/TT)/unit_mass*unit_length**3
          f(:,:,:,ilnrho)=alog(tmp)
        endif
        if (ltemperature_nolog) then
          if (ldensity_nolog) then
            f(:,:,:,iTT)=(PP/(k_B_cgs/m_u_cgs)*&
                air_mass/f(:,:,:,ilnrho))/unit_mass*unit_length**3
          else
            f(:,:,:,iTT)=(PP/(k_B_cgs/m_u_cgs)*&
                air_mass/exp(f(:,:,:,ilnrho)))/unit_mass*unit_length**3
          endif
        else
          if (ldensity_nolog) then
            tmp=(PP/(k_B_cgs/m_u_cgs)*&
                air_mass/f(:,:,:,ilnrho))/unit_mass*unit_length**3
            f(:,:,:,ilnTT)=alog(tmp)!+f(:,:,:,ilnTT)
          else
            tmp=(PP/(k_B_cgs/m_u_cgs)*&
                air_mass/exp(f(:,:,:,ilnrho)))/unit_mass*unit_length**3
            f(:,:,:,ilnTT)=alog(tmp)!+f(:,:,:,ilnTT)
          endif
        endif
        if (velx/=0.) f(:,:,:,iux)=f(:,:,:,iux)+velx
      endif

      if (linit_temperature) then
        do i=1,mx
          if (x(i)<=init_x1) f(i,:,:,ilnTT)=alog(init_TT1)
          if (x(i)>=init_x2) f(i,:,:,ilnTT)=alog(init_TT2)
          if (x(i)>init_x1 .and. x(i)<init_x2) then
            if (init_x1 /= init_x2) &
              f(i,:,:,ilnTT)=&
                 alog((x(i)-init_x1)/(init_x2-init_x1)*(init_TT2-init_TT1)+init_TT1)
          endif
        enddo
      endif

      if (linit_density) then
        if (ldensity_nolog) then
          f(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)*air_mass/exp(f(:,:,:,ilnTT)))/unit_mass*unit_length**3
        else
          tmp=(PP/(k_B_cgs/m_u_cgs)*air_mass/exp(f(:,:,:,ilnTT)))/unit_mass*unit_length**3
          f(:,:,:,ilnrho)=alog(tmp)
        endif
      endif

      if (nxgrid>1) f(:,:,:,iux)=f(:,:,:,iux)+velx

      if (lroot) then
        print*, 'Air temperature, K', TT
        print*, 'Air pressure, dyn', PP
        print*, 'Air density, g/cm^3:'
        print '(E10.3)',  PP/(k_B_cgs/m_u_cgs)*air_mass/TT
        print*, 'Air mean weight, g/mol', air_mass
        print*, 'R', k_B_cgs/m_u_cgs

        do j=1,nchemspec
          spe="Y("//trim(varname(ichemspec(j)))//")="
          write (*,'(A10,F10.7)') spe,initial_massfractions(j)
        enddo

      endif

      close(file_id)
!
    endsubroutine air_field

    !***********************************************************************
    subroutine premixed_equiv_ratio(f)
      !      
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nchemspec) :: conc=0.0, initial_massfractions=0.0, dens_spec=0.0
      real :: conc_total, dens_total, mean_molar_mass
      real :: rel_mass_oxidizer, rho, TT
      integer :: i, j, ichem_O2, ichem_N2, ind_chem, ind_glob
      logical :: found_specie
      character(len=30) :: spe, specie_string
!
      !init_premixed_fuel(1)="SiO"
      !init_premixed_fuel(2)="CO"
      !init_fuel_molar_ratio(1)=0.5
      !init_fuel_molar_ratio(2)=0.5
      !init_fuel_O2_demand(1)=0.5
      !init_fuel_O2_demand(2)=0.5
      !phi=1.0

      print*,"init_premixed_fuel=",init_premixed_fuel
      
      call find_species_index("O2",ind_glob,ichem_O2,found_specie)
      conc(ichem_O2)=1.
      call find_species_index("N2",ind_glob,ichem_N2,found_specie)
      conc(ichem_N2)=conc(ichem_O2)*79./21.
      !
      ! Check that the molar ratios of the fuel sums to unity
      !
      if (sum(init_fuel_molar_ratio) .ne. 1.0) then
        !call fatal_error("premixed_equiv_ratio","init_fuel_molar_ratio must sum to unity.")
      endif
      !
      ! Loop over all species and check which that are present
      !
      do i=1,nchemspec
        if (init_fuel_molar_ratio(i) .gt. 0) then
          specie_string=init_premixed_fuel(i)
          call find_species_index(specie_string,ind_glob,ind_chem,found_specie)
          conc(ind_chem)=conc(ichem_O2)*init_phi*init_fuel_molar_ratio(i)*init_fuel_O2_demand(i)
          print*,"i,specie_string,ind_glob_ind_chem,found_specie=",i,specie_string,ind_glob,ind_chem,found_specie,conc(ind_chem)
        endif
      enddo
      !
      ! Find total molar concentration and use this to determine mass fractions
      !
      conc_total=sum(conc)
      dens_spec=conc*species_constants(:,imass)
      dens_total=sum(dens_spec)
      mean_molar_mass=0.
      do i=1,nchemspec
        initial_massfractions(i)=dens_spec(i)/dens_total
        f(:,:,:,ichemspec(i))=initial_massfractions(i)
        mean_molar_mass=mean_molar_mass+conc(i)*species_constants(i,imass)/conc_total
      enddo
      !
      ! Calculate temperature based on temperature of fuel and oxidizer streams
      ! Here we assume that the heat capacities of the two streams are equal
      !
      rel_mass_oxidizer=initial_massfractions(ichem_O2)+initial_massfractions(ichem_N2)
      TT=init_temp_oxidizer*(rel_mass_oxidizer)+init_temp_fuel*(1-rel_mass_oxidizer)
      if (ltemperature_nolog) then
        f(:,:,:,iTT)=TT
      else
        f(:,:,:,ilnTT)=alog(TT)
      endif
      !
      ! Set density
      !
      rho=init_pressure*mean_molar_mass/(Rgas*TT)
      if (ldensity_nolog) then
        f(:,:,:,irho)=rho
      else
        f(:,:,:,ilnrho)=alog(rho)
      endif
      !
      ! Print results
      !
      if (lroot) then
        print*, 'Oxidizer temperature, [K] =', init_temp_oxidizer
        print*, 'Fuel temperature, [K]     =', init_temp_fuel
        print*, 'Mixture temperature, [K]  =', TT
        print*, 'Pressure, [dyn]           =', init_pressure
        print*, 'Density, [g/cm^3]         =', rho
        do j=1,nchemspec
          spe="Y("//trim(varname(ichemspec(j)))//")="
          write (*,'(A10,F10.7)') spe,initial_massfractions(j)
        enddo
        do j=1,nchemspec
          spe="conc("//trim(varname(ichemspec(j)))//")="
          write (*,'(A10,F10.7)') spe,conc(j)
        enddo
      endif
      !
      if (sum(initial_massfractions) .ne. 1.0) then
        call fatal_error("premixed_equiv_ratio","initial_massfractions must sum to unity.")
      endif
!
    endsubroutine premixed_equiv_ratio
!***********************************************************************
!          NSCBC boundary conditions
!***********************************************************************
    subroutine damp_zone_for_NSCBC(f,df)
!
!   16-jul-06/natalia: coded
!   24-jan-11/julien: modified
!    buffer zone to damp the acustic waves!!!!!!!!!!!
!    important for NSCBC
!    Most of this routine is still hard-coded
!    Should contain more tests to detect boundaries and set the reference conditions
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      integer :: sz_x, sz_y, sz_z, ll1, ll2, i
!      integer :: sz_r_x, sz_l_x, sz_r_y, sz_l_y, sz_r_z, sz_l_z
      real :: dt1, func_x !,func_y,func_z
      real :: ux_ref, uy_ref, uz_ref, lnTT_ref, lnrho_ref, gamma, cs
      real :: del
!      logical :: lzone_y=.false.,lzone_z=.false.
      logical :: dir_damp1=.true. !, dir_damp2=.false., dir_damp3=.false.
      logical :: lright=.true., lleft=.true.
!
      ux_ref = 0.
      uy_ref = 0.
      uz_ref = 0.
!      lnTT_ref=6.39693
      lnTT_ref = log(init_TT1)
      lnrho_ref = log(init_pressure)-log(Rgas)-lnTT_ref-log(mu1_full(l1,m1,n1))
      gamma = 1.4
      cs = sqrt(gamma*init_pressure/exp(lnrho_ref))
!
!  The characteristic time scale is taken of the order of a CFL time scale
!
      dt1 = (ux_ref+cs)/(Lxyz(1)/nxgrid)
      del = 0.1
!
      sz_x = int(del*nxgrid)
      sz_y = int(del*nygrid)
      sz_z = int(del*nzgrid)
      ll1 = l1
      ll2 = l2
!
      if (nxgrid /= 1 .and. dir_damp1) then
!
!  On the right side
!
        if (lright) then
          ll1 = nxgrid-sz_x
          ll2 = nxgrid
          do i = l1,l2
            if (x(i) >= xgrid(ll1)) then
              func_x = (x(i)-xgrid(ll1))**3/(xgrid(ll2)-xgrid(ll1))**3
              df(i,m,n,iux) = df(i,m,n,iux)-func_x*(f(i,m,n,iux)-ux_ref)*dt1
              df(i,m,n,iuy) = df(i,m,n,iuy)-func_x*(f(i,m,n,iuy)-uy_ref)*dt1
              df(i,m,n,iuz) = df(i,m,n,iuz)-func_x*(f(i,m,n,iuz)-uz_ref)*dt1
!            df(i,m,n,ilnrho)=df(i,m,n,ilnrho)-func_x*(f(i,m,n,ilnrho)-lnrho_ref)*dt1
!            df(i,m,n,ilnTT)=df(i,m,n,ilnTT)-func_x*(f(i,m,n,ilnTT)-lnTT_ref)*dt1
            endif
          enddo
        endif
!
!  On the left side
!
        if (lleft) then
          ll1 = 1
          ll2 = sz_x
          do i = l1,l2
            if (x(i) <= xgrid(ll2)) then
              func_x = (x(i)-xgrid(ll2))**3/(xgrid(ll1)-xgrid(ll2))**3
              df(i,m,n,iux) = df(i,m,n,iux)-func_x*(f(i,m,n,iux)-ux_ref)*dt1
              df(i,m,n,iuy) = df(i,m,n,iuy)-func_x*(f(i,m,n,iuy)-uy_ref)*dt1
              df(i,m,n,iuz) = df(i,m,n,iuz)-func_x*(f(i,m,n,iuz)-uz_ref)*dt1
!            df(i,m,n,ilnrho)=df(i,m,n,ilnrho)-func_x*(f(i,m,n,ilnrho)-lnrho_ref)*dt1
!            df(i,m,n,ilnTT)=df(i,m,n,ilnTT)-func_x*(f(i,m,n,ilnTT)-lnTT_ref)*dt1
            endif
          enddo
        endif
!
      endif
!
!  The following should be replaced to generalize the formula
!
!       if (nygrid/=1 .and. dir_damp2) then
!
!       if (sz_r_y<=m1) call fatal_error('to use ldamp_zone_NSCBC',&
!                  'you should increase nygrid!')
!
!       if ((m<=sz_l_y) .and. (m>=m1)) then
!        func_y=(y(m)-y(sz_l_y))**3/(y(m1)-y(sz_l_y))**3
!        lzone_y=.true.
!       elseif ((m>=sz_r_y) .and. (m<=m2)) then
!        func_y= (y(m)-y(sz_r_y))**3/(y(m2)-y(sz_r_y))**3
!        lzone_y=.true.
!       endif
!
!       if (lzone_y) then
!        df(sz_l_x:sz_r_x,m,n,iux)=df(sz_l_x:sz_r_x,m,n,iux)&
!           -func_y*(f(sz_l_x:sz_r_x,m,n,iux)-ux_ref)*dt1
!        df(sz_l_x:sz_r_x,m,n,iuy)=df(sz_l_x:sz_r_x,m,n,iuy)&
!           -func_y*(f(sz_l_x:sz_r_x,m,n,iuy)-uy_ref)*dt1
!        df(sz_l_x:sz_r_x,m,n,iuz)=df(sz_l_x:sz_r_x,m,n,iuz)&
!           -func_y*(f(sz_l_x:sz_r_x,m,n,iuz)-uz_ref)*dt1
!        df(sz_l_x:sz_r_x,m,n,ilnrho)=df(sz_l_x:sz_r_x,m,n,ilnrho)&
!           -func_y*(f(sz_l_x:sz_r_x,m,n,ilnrho)-lnrho_ref)*dt1
!        df(sz_l_x:sz_r_x,m,n,ilnTT)=df(sz_l_x:sz_r_x,m,n,ilnTT)&
!           -func_y*(f(sz_l_x:sz_r_x,m,n,ilnTT)-lnTT_ref)*dt1
!        lzone_y=.false.
!       endif
!       endif
!
!      if (nzgrid/=1 .and. dir_damp3) then
!        if (sz_r_z<=n1) call fatal_error('to use ldamp_zone_NSCBC',&
!                  'you should increase nzgrid!')
!
!      if ((n<=sz_l_z) .and. (n>=n1)) then
!        func_z=(z(n)-z(sz_l_z))**3/(z(n1)-z(sz_l_z))**3
!        lzone_z=.true.
!       elseif ((n>=sz_r_z) .and. (n<=n2)) then
!        func_z= (z(n)-z(sz_r_z))**3/(z(n2)-z(sz_r_z))**3
!        lzone_z=.true.
!       endif
!
!       if (lzone_z) then
!        df(sz_l_x:sz_r_x,m,n,iux)=df(sz_l_x:sz_r_x,m,n,iux)&
!           -func_z*(f(sz_l_x:sz_r_x,m,n,iux)-ux_ref)*dt1
!        df(sz_l_x:sz_r_x,m,n,iuy)=df(sz_l_x:sz_r_x,m,n,iuy)&
!           -func_z*(f(sz_l_x:sz_r_x,m,n,iuy)-uy_ref)*dt1
!        df(sz_l_x:sz_r_x,m,n,iuz)=df(sz_l_x:sz_r_x,m,n,iuz)&
!           -func_z*(f(sz_l_x:sz_r_x,m,n,iuz)-uz_ref)*dt1
!        df(sz_l_x:sz_r_x,m,n,ilnrho)=df(sz_l_x:sz_r_x,m,n,ilnrho)&
!           -func_z*(f(sz_l_x:sz_r_x,m,n,ilnrho)-lnrho_ref)*dt1
!        df(sz_l_x:sz_r_x,m,n,ilnTT)=df(sz_l_x:sz_r_x,m,n,ilnTT)&
!           -func_z*(f(sz_l_x:sz_r_x,m,n,ilnTT)-lnTT_ref)*dt1
!        lzone_z=.false.
!       endif
!       endif
!
    endsubroutine damp_zone_for_NSCBC
!***********************************************************************
    subroutine jacobn(f,jacob)
!
! Compute the jacobian, i.e. the matrix  jacob(nchemspec x nchemspec)
! where jacob(i,j)=dv_i/dc_j
! v is the vector of dimension nchemspec of the rates dc_j/dt
! (the c_j being concentrations, stocked in f among other)
!
!  28-may-09/rplasson: coded
!
!   exchange data
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,nchemspec,nchemspec) :: jacob
!
      intent(in) :: f
      intent(out) :: jacob
!   internal data
!
!   indices
      integer :: i, j, k, l, ii
!
!   temporary
      real :: tmp_p, tmp_m
!   Code
!
      jacob = 0.
!
!  identify module
!
      if (headtt .or. ldebug) print*,'jacobn: compute the jacobian matrix'
!
      if (lreactions) then
        do n = n1,n2
          do m = m1,m2
            do l = l1,l2
              do i = 1,nchemspec
                do j = 1,nchemspec
! Compute dv_i/dc_j
!              print*,"dv_",i,"/dc_",j,"="
                  do k = 1,nreactions
! Check if compound i participate in reaction k
                    if (Sijp(i,k) /= Sijm(i,k)) then
! Compute the contribution of reaction k to dv_i/dc_j
!                  print*,"+(",(Sijp(i,k)-Sijm(i,k)),"k_",k,"(+/-)"
                      tmp_p = (Sijp(i,k)-Sijm(i,k))*kreactions_p(k)*kreactions_z(n,k)
                      tmp_m = (Sijm(i,k)-Sijp(i,k))*kreactions_m(k)*kreactions_z(n,k)
                      do  ii = 1,nchemspec
! Compute the contribution of compound ii in reaction k
!                    print*,"**-",Sijm(ii,k),Sijp(ii,k),"-**"
                        if (ii /= j) then
!                      print*,"c_",ii,"^",Sijm(ii,k)," (/) ","c_",ii,"^",Sijp(ii,k)
                          tmp_p = tmp_p*f(l,m,n,ichemspec(ii))**Sijm(ii,k)
                          tmp_m = tmp_m*f(l,m,n,ichemspec(ii))**Sijp(ii,k)
                        else
                          if (Sijm(ii,k) == 0) then
!                        print*,"0*c_",ii
                            tmp_p = 0.
                          elseif (Sijm(ii,k) > 1) then
!                        print*,Sijm(ii,k),"*c_",ii,"^",(Sijm(ii,k)-1)
                            tmp_p = Sijm(ii,k)*tmp_p*f(l,m,n,ichemspec(ii))**(Sijm(ii,k)-1)
!                      else
!                        print*,"c_",ii,"^0"
                          endif
!                      print*," (/) "
                          if (Sijp(ii,k) == 0) then
!                        print*,"0*c_",ii
                            tmp_m = 0.
                          elseif (Sijp(ii,k) > 1) then
!                        print*,Sijp(ii,k),"*c_",ii,"^",(Sijp(ii,k)-1)
                            tmp_m = Sijp(ii,k)*tmp_m*f(l,m,n,ichemspec(ii))**(Sijp(ii,k)-1)
!                      else
!                        print*,"c_",ii,"^0"
                          endif
                        endif
                      enddo
!                  print*,")"
! Add the contribution of reaction k to dv_i/dc_j
                      jacob(l,m,n,i,j) = jacob(l,m,n,i,j)+tmp_p+tmp_m
!                  print*,"(=",tmp_p," (-) ",tmp_m,")"
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
!
    endsubroutine jacobn
!***********************************************************************
    subroutine get_mu1_slice(f,slice,grad_slice,index,sgn,direction)
!
! For the NSCBC boundary conditions the slice of mu1 at the boundary, and
! its gradient, is required.
!
!  10-dez-09/Nils Erland L. Haugen: coded
!
      use Deriv, only: der_onesided_4_slice_other
!
      real, dimension(:,:), intent(out) :: slice
      real, dimension(:,:), intent(out) :: grad_slice
      integer, intent(in) :: index, sgn, direction
      real, dimension(mx,my,mz,mfarray) :: f
!
      if (direction == 1) then
        slice = mu1_full(index,m1:m2,n1:n2)
        call der_onesided_4_slice_other(mu1_full,sgn,grad_slice,index,direction)
      elseif (direction == 2) then
        slice = mu1_full(l1:l2,index,n1:n2)
        call der_onesided_4_slice_other(mu1_full,sgn,grad_slice,index,direction)
      else
        slice = mu1_full(l1:l2,m1:m2,index)
        call der_onesided_4_slice_other(mu1_full,sgn,grad_slice,index,direction)
      endif
!
    endsubroutine get_mu1_slice
!***********************************************************************
    subroutine calc_extra_react(f,reac,kf_loc,i,mm,nn,p)
!
!   character (len=*), intent(in) :: element_name
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case), intent(in) :: p
      integer, intent(in) :: reac, i, mm, nn
      real, intent(out) :: kf_loc
      real :: X_O2, X_N2, X_O2N2, X_M, X_H2O
      real :: K1, K2, K3, K4, KMT06
!
!
      select case (reaction_name(reac))
      case ('O=O3')
        X_O2 = f(l1+i-1,mm,nn,ichemspec(index_O2))*unit_mass &
            /species_constants(index_O2,imass)*p%rho(i)
        X_N2 = f(l1+i-1,mm,nn,ichemspec(index_N2))*unit_mass &
            /species_constants(index_N2,imass)*p%rho(i)
        kf_loc = 5.60D-34*X_O2*X_N2*((1./300.*p%TT(i))**(-2.6)) &
            +6.00D-34*X_O2**2*((1./300.*p%TT(i))**(-2.6))
      case ('O1D=O')
        X_O2 = f(l1+i-1,mm,nn,ichemspec(index_O2))*unit_mass &
            /species_constants(index_O2,imass)*p%rho(i)
        X_N2 = f(l1+i-1,mm,nn,ichemspec(index_N2))*unit_mass &
            /species_constants(index_N2,imass)*p%rho(i)
        kf_loc = 3.20D-11*X_O2*exp(67.*p%TT1(i))+1.80D-11*X_N2*exp(107.*p%TT1(i))
      case ('OH+CO=HO2')
        X_O2N2 = f(l1+i-1,mm,nn,ichemspec(index_O2N2))*unit_mass &
            /species_constants(index_O2N2,imass)*p%rho(i)
        kf_loc = 1.30D-13*(1+((0.6*index_O2N2)/(2.652E+19*(300.*p%TT1(i)))))
      case ('2HO2=H2O2')
        X_M = p%rho(i)*mu1_full(l1+i-1,mm,nn)
        X_H2O = f(l1+i-1,mm,nn,ichemspec(index_H2O))*unit_mass &
            /species_constants(index_H2O,imass)*p%rho(i)
        KMT06 = 1 + (1.4E-21 * EXP(2200.*p%TT1(i)) * X_H2O)
        kf_loc = 2.20D-13*KMT06*EXP(600.*p%TT1(i)) &
            + 1.90D-33*X_M*KMT06*EXP(980.*p%TT1(i))
      case ('OH+HNO3=NO3')
        X_O2N2 = f(l1+i-1,mm,nn,ichemspec(index_O2N2))*unit_mass &
                 /species_constants(index_O2N2,imass)*p%rho(i)
        K1     =  2.4E-14 * EXP(460.*p%TT1(i))
        K3     =  6.5E-34 * EXP(1335.*p%TT1(i))
        K4     =  2.7E-17 * EXP(2199.*p%TT1(i))
        K2     =  (K3 * X_O2N2) / (1 + (K3*X_O2N2/K4))
        kf_loc = K1 + K2
      case default
        call fatal_error('calc_extra_react','no such reaction_name:  "'//trim(reaction_name(reac))//'"')
      endselect
!
    endsubroutine calc_extra_react
!***********************************************************************
    subroutine prerun_1D(f,prerun_dir)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!  11-apr-97/axel: coded
!
      character(len=*) :: prerun_dir
      character(len=100) :: file
      character(len=10) :: processor
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(:,:,:,:), allocatable :: a
      real, dimension(:,:,:), allocatable :: grid
      real, dimension(:), allocatable :: cc
      real :: t, sum
      integer :: j, k, i, ii, imid, ipos
      integer :: mx2, my2, mz2, mv2
      integer :: l22, m22
!
! NILS: For the time beeing the prerun simulation can not have been run in parallell
! NILS: Should fix this soon.
!
! Read dimension of the array stored during the 1D prerun
!
      processor = 'proc0'
      file = trim(prerun_dir)//'/data/'//trim(processor)//'/dim.dat'
      print*,'Reading prerun dim from ',file
      open (1,FILE=file)
      read (1,*) mx2, my2, mz2, mv2
      close (1)
      l22 = mx2-l1+1
      m22 = mx2-2*(l1-1)
      allocate(a(mx2,my2,mz2,mv2), grid(m22,my2-6,mz2-6), cc(m22))
!
! Read the grid used during the prerun
!
      file = trim(prerun_dir)//'/data/'//trim(processor)//'/grid.dat'
      print*,'Reading prerun grid from ',file
      open (1,FILE=file,FORM='unformatted')
      read (1) t,grid(:,1,1),grid(1,:,1),grid(1,1,:)
      close (1)
!
! Read the data stored during the prerun
!
      file = trim(prerun_dir)//'/data/'//trim(processor)//'/var.dat'
      print*,'Reading inlet data from ',file
      open (1,FILE=file,FORM='unformatted')
      read (1) a
      close (1)
!
! Definition of a progress variable
!
      cc(1) = 0.
      ipos = 0
      imid = 0
      do ii = 2, m22-1
        cc(ii) = (exp(a(ii,m1,n1,iuz+2)) - exp(a(l1,m1,n1,iuz+2)))/ &
                 (exp(a(m22-1,m1,n1,iuz+2)) - exp(a(l1,m1,n1,iuz+2)))
        if (cc(ii) > 0.7 .and. cc(ii-1) <= 0.7) imid = ii
        if (grid(ii,1,1) > flame_pos .and. grid(ii-1,1,1) <= flame_pos) ipos = ii
        if (ipos > 0 .and. imid > 0) exit
      enddo
!
!  The center of the flame, arbitrarily identified by cc=0.7, is
!  located at the position flame_pos (chemistry_init_pars)
!  The flame is thus moved in its domain
!
!  Spread the data on the f-array
!
      do j = 1,my
        do k = 1,mz
          do i = 1,mx
!
            do ii = 2, m22-1
              if (x(i) > grid(ii,1,1)-(grid(imid,1,1)-grid(ipos,1,1)) .and. x(i) &
                  <= grid(ii+1,1,1)-(grid(imid,1,1)-grid(ipos,1,1))) then
                if (.not. lperi(1)) then
                  f(i,j,k,iux) = f(i,j,k,iux)+a(ii,m1,n1,iux)+(x(i)-grid(ii,1,1) &
                      +(grid(imid,1,1)-grid(ipos,1,1)))*(a(ii+1,m1,n1,iux) &
                      -a(ii,m1,n1,iux))/(grid(ii+1,1,1)-grid(ii,1,1))
                endif
                f(i,j,k,iuz+1:mvar) = a(ii,m1,n1,iuz+1:mvar)+(x(i)-grid(ii,1,1)+ &
                    (grid(imid,1,1)-grid(ipos,1,1)))*(a(ii+1,m1,n1,iuz+1:mvar) &
                    -a(ii,m1,n1,iuz+1:mvar))/(grid(ii+1,1,1)-grid(ii,1,1))
              elseif (x(i) <= grid(2,1,1)-(grid(imid,1,1)-grid(ipos,1,1))) then
                if (.not. lperi(1)) f(i,j,k,iux) = f(i,j,k,iux)+a(l1,m1,n1,iux)
                f(i,j,k,iuz+1:mvar) = a(l1,m1,n1,iuz+1:mvar)
                exit
              elseif (x(i) >= grid(m22-1,1,1)-(grid(imid,1,1)-grid(ipos,1,1))) then
                if (.not. lperi(1)) f(i,j,k,iux) = f(i,j,k,iux)+a(l22,m1,n1,iux)
                f(i,j,k,iuz+1:mvar) = a(l22,m1,n1,iuz+1:mvar)
                exit
              endif
            enddo
!
!  Renormalize the species mass fractions
!
            sum = 0.
            do ii = 1, nchemspec
              sum = sum + f(i,j,k,ichemspec(ii))
            enddo
            f(i,j,k,ichemspec(1):ichemspec(nchemspec)) = f(i,j,k,ichemspec(1):ichemspec(nchemspec))/sum
!
          enddo
        enddo
      enddo
!
!  Set the y and z velocities to zero in order to avoid random noise
!
!      if (.not. reinitialize_chemistry) then
!        f(:,:,:,iuy:iuz)=0
!      endif
!
      deallocate(a,grid,cc)
!
    endsubroutine prerun_1D
!***********************************************************************
    subroutine prerun_1D_opp(f,prerun_dir)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!  04-aug-10/julien: coded
!
      character(len=*) :: prerun_dir
      character(len=100) :: file
      character(len=10) :: processor
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(:,:,:,:), allocatable :: a
      real, dimension(:,:,:), allocatable :: grid
      real :: t
      real :: x1, x2, xm, xfl
      integer :: j, k, i, ii, ifl
      integer :: mx2, my2, mz2, mv2
      integer :: l22, m22
!
      x1 = xyz0(1)+Lxyz(1) / 4.
      x2 = x1+Lxyz(1) / 2.
      xm = xyz0(1)+Lxyz(1) / 2.
!
! Read dimension of the array stored during the 1D prerun
!
      processor = 'proc0'
      file = trim(prerun_dir)//'/data/'//trim(processor)//'/dim.dat'
      print*,'Reading prerun dim from ',file
      open (1,FILE=file)
      read (1,*) mx2, my2, mz2, mv2
      close (1)
      l22 = mx2-l1+1
      m22 = mx2-2*(l1-1)
      allocate(a(mx2,my2,mz2,mv2), grid(m22,my2-6,mz2-6))
!
! Read the grid used during the prerun
!
      file = trim(prerun_dir)//'/data/'//trim(processor)//'/grid.dat'
      print*,'Reading prerun grid from ',file
      open (1,FILE=file,FORM='unformatted')
      read (1) t,grid(:,1,1),grid(1,:,1),grid(1,1,:)
      close (1)
!
! Read the data stored during the prerun
!
      file = trim(prerun_dir)//'/data/'//trim(processor)//'/var.dat'
      print*,'Reading inlet data from ',file
      open (1,FILE=file,FORM='unformatted')
      read (1) a
      close (1)
!
      do ii = 1, m22-1
        if (exp(a(ii,m1,n1,iuz+2)) < 1200. .and. exp(a(ii+1,m1,n1,iuz+2)) >= 1200.) then
          xfl = grid(ii,1,1)
          ifl = ii
          exit
        endif
      enddo
!
!  Spread the data on the f-array : interpolations
!
      do j = 1,my
        do k = 1,mz
          do i = 1,mx
!
            if (x(i) < xm) then
              do ii = 2, m22-1
                if (x(i)-x1 > grid(ii,1,1)-xfl .and. x(i)-x1 <= grid(ii+1,1,1)-xfl) then
                  f(i,j,k,iuz+1:mvar) = a(ii,m1,n1,iuz+1:mvar)+((x(i)-x1)-(grid(ii,1,1)-xfl))* &
                      (a(ii+1,m1,n1,iuz+1:mvar)-a(ii,m1,n1,iuz+1:mvar))/(grid(ii+1,1,1)-grid(ii,1,1))
                  exit
                elseif (x(i)-x1 <= grid(l1,1,1)-xfl) then
                  f(i,j,k,iuz+1:mvar) = a(l1,m1,n1,iuz+1:mvar)
                  exit
                elseif (x(i)-x1 >= grid(m22,1,1)-xfl) then
                  f(i,j,k,iuz+1:mvar) = a(l22,m1,n1,iuz+1:mvar)
                  exit
                endif
              enddo
!
            elseif (x(i) >= xm) then
              do ii = 2, m22-1
                if (x(i)-x2 > grid(ii,1,1)-xfl .and. x(i)-x2 <= grid(ii+1,1,1)-xfl &
                    .and. 2*ifl > ii+1) then
                  f(i,j,k,iuz+1:mvar) = a(2*ifl-ii,m1,n1,iuz+1:mvar)+((x(i)-x2)-(grid(ii,1,1)-xfl)) &
                      *(a(2*ifl-ii-1,m1,n1,iuz+1:mvar)-a(2*ifl-ii,m1,n1,iuz+1:mvar))/(grid(ii+1,1,1)-grid(ii,1,1))
                  exit
                elseif (x(i)-x2 <= grid(l1,1,1)-xfl) then
                  f(i,j,k,iuz+1:mvar) = a(l22,m1,n1,iuz+1:mvar)
                  exit
                elseif (x(i)-x2 >= grid(m22-1,1,1)-xfl .or. 2*ifl <= ii+1) then
                  f(i,j,k,iuz+1:mvar) = a(l1,m1,n1,iuz+1:mvar)
                  exit
                endif
              enddo
            endif
          enddo
        enddo
      enddo
!
!  Set the y and z velocities to zero in order to avoid random noise
!
!      if (.not. reinitialize_chemistry) f(:,:,:,iuy:iuz)=0
!
      deallocate(a,grid)
!
    endsubroutine prerun_1D_opp
!***********************************************************************
    subroutine FlameMaster_ini(f,file_name)
!
!  read FlameMaster file for flame initialization
!  11-nov-10/julien: coded
!
      character(len=*) :: file_name
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(:,:), allocatable :: a
      real, dimension(:), allocatable :: grid
      real, dimension(:), allocatable :: cc
      integer :: j, k, i, ii, imid, ipos
      integer :: is, js
      integer :: nsp, npts
      character(len=8), dimension(:), allocatable :: name_sp
      character(len=10) :: car10='nothing'
      character(len=12) :: car12
      character(len=20) :: car20='nothing'
      character(len=1) :: car1
!
      lFlameMaster = .true.
!
! Read dimension of the array stored in the FlameMaster initial file
!
      open (1,FILE=trim(file_name))
      if (lroot) print*, 'Reading initial conditions in file ', trim(file_name)
      do while (car10 /= 'FlameThick')
        read (1,*) car10
      enddo
      read (1,*) car12, car1, nsp
      read (1,*) car12, car1, npts
      allocate(a(npts,nsp+3), grid(npts), cc(npts))
      allocate(name_sp(nsp))
!
! Read the data stored during the prerun
!
      do while (car10 /= 'body')
        read (1,*) car10
      enddo
      read (1,*)
      read (1,*) grid
      grid = grid*100.
!
      i = 0
      do while (trim(car20) /= 'trailer')
        read (1,*) car20
        select case (car20(1:12))
        case ('massflowrate')
          read (1,*) a(:,1)
        case ('temperature')
          read (1,*) a(:,3)
        case ('massfraction')
          i = i+1
          read (1,*) a(:,3+i)
          name_sp(i) = car20(14:20)
        case ('density')
          read (1,*) a(:,2)
        endselect
      enddo
      close (1)
      a(:,1) = a(:,1) / a(:,2)
      a(:,1) = a(:,1)*100.
      a(:,2) = a(:,2)/1000.
!
! Definition of a progress variable
!
      cc(1) = 0.
      ipos = 0
      imid = 0
      do ii = 2, npts
        cc(ii) = (a(ii,3) - a(1,3)) / (a(npts,3) - a(1,3))
        if (cc(ii) > 0.7 .and. cc(ii-1) <= 0.7) imid = ii
        if (grid(ii) > flame_pos .and. grid(ii-1) <= flame_pos) ipos = ii
        if (ipos > 0 .and. imid > 0) exit
      enddo
!
!  The center of the flame, arbitrarily identified by cc=0.7, is
!  located at the position flame_pos (chemistry_init_pars)
!  The flame is thus moved in its domain
!
!  Spread the data on the f-array : interpolations
!
      do j = 1,my
        do k = 1,mz
          do i = 1,mx
!
            do ii = 2, npts-1
              if (x(i) > grid(ii)-(grid(imid)-grid(ipos)) .and. x(i) &
                  <= grid(ii+1)-(grid(imid)-grid(ipos))) then
                if (.not. init_from_file) &
                    f(i,j,k,iux) = f(i,j,k,iux)+a(ii,1)+(x(i)-grid(ii)+(grid(imid)-grid(ipos)))* &
                    (a(ii+1,1)-a(ii,1)) / (grid(ii+1)-grid(ii))
                f(i,j,k,iuz+2) = a(ii,3)+(x(i)-grid(ii)+(grid(imid)-grid(ipos)))* &
                    (a(ii+1,3)-a(ii,3)) / (grid(ii+1)-grid(ii))
                f(i,j,k,iuz+1) = a(ii,2)+(x(i)-grid(ii)+(grid(imid)-grid(ipos)))* &
                    (a(ii+1,2)-a(ii,2)) / (grid(ii+1)-grid(ii))
                do is = 1, nsp
                  do js = 1, nchemspec
                    if (trim(name_sp(is)) == trim(varname(ichemspec(js)))) then
                      f(i,j,k,iuz+2+js) = a(ii,is+3)+(x(i)-grid(ii)+(grid(imid)-grid(ipos)))*   &
                          (a(ii+1,is+3)-a(ii,is+3)) / (grid(ii+1)-grid(ii))
                      exit
                    endif
                  enddo
                enddo
!
              elseif (x(i) <= grid(2)-(grid(imid)-grid(ipos))) then
                if (.not. init_from_file) f(i,j,k,iux) = f(i,j,k,iux)+a(1,1)
                f(i,j,k,iuz+2) = a(1,3)
                f(i,j,k,iuz+1) = a(1,2)
                do is = 1, nsp
                  do js = 1, nchemspec
                    if (trim(name_sp(is)) == trim(varname(ichemspec(js)))) then
                      f(i,j,k,iuz+2+js) = a(1,is+3)
                      exit
                    endif
                  enddo
                enddo
                exit
!
              elseif (x(i) >= grid(npts)-(grid(imid)-grid(ipos))) then
                if (.not. init_from_file) f(i,j,k,iux) = f(i,j,k,iux)+a(npts,1)
                f(i,j,k,iuz+2) = a(npts,3)
                f(i,j,k,iuz+1) = a(npts,2)
                do is = 1, nsp
                  do js = 1, nchemspec
                    if (trim(name_sp(is)) == trim(varname(ichemspec(js)))) then
                      f(i,j,k,iuz+2+js) = a(npts,is+3)
                      exit
                    endif
                  enddo
                enddo
                exit
              endif
            enddo
!
!  Renormalize the species mass fractions
!
            f(i,j,k,iuz+3:iuz+2+nchemspec) = f(i,j,k,iuz+3:iuz+2+nchemspec) &
                / sum(f(i,j,k,iuz+3:iuz+2+nchemspec))
!
          enddo
        enddo
      enddo
      if (.not. ldensity_nolog) f(:,:,:,iuz+1) = log(f(:,:,:,iuz+1))
      if (.not. ltemperature_nolog) f(:,:,:,iuz+2) = log(f(:,:,:,iuz+2))
!
!  Set the y and z velocities to zero in order to avoid random noise
!
!      if (.not. reinitialize_chemistry) f(:,:,:,iuy:iuz)=0
!
      deallocate(a,grid,name_sp)
!
    endsubroutine FlameMaster_ini
!***********************************************************************
    subroutine get_reac_rate(wt,f,p)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx) :: wt
      real, dimension(nx,nchemspec) :: ydot
      type (pencil_case) :: p
!
      ydot = p%DYDt_reac
      f(l1:l2,m,n,ireaci(1):ireaci(nchemspec)) = ydot
!         f(l1:l2,m,n,ireaci(1):ireaci(nchemspec))+ydot
!
      if (maux == nchemspec+1) f(l1:l2,m,n,ireaci(nchemspec)+1) = wt
!         f(l1:l2,m,n,ireaci(nchemspec)+1)+wt
!
    endsubroutine get_reac_rate
!***********************************************************************
    subroutine chemistry_clean_up
!
      if (allocated(Bin_diff_coef))  deallocate(Bin_diff_coef)
      if (allocated(stoichio))       deallocate(stoichio)
      if (allocated(vreactions_p))   deallocate(vreactions_p)
      if (allocated(vreactions_m))   deallocate(vreactions_m)
      if (allocated(Sijm))           deallocate(Sijm)
      if (allocated(Sijp))           deallocate(Sijp)
      if (allocated(kreactions_z))   deallocate(kreactions_z)
      if (allocated(kreactions_m))   deallocate(kreactions_m)
      if (allocated(kreactions_p))   deallocate(kreactions_p)
      if (allocated(reaction_name))  deallocate(reaction_name)
      if (allocated(enum_reaction_name))  deallocate(enum_reaction_name)
      if (allocated(B_n))            deallocate(B_n)
      if (allocated(alpha_n))        deallocate(alpha_n)
      if (allocated(E_an))           deallocate(E_an)
      if (allocated(low_coeff))      deallocate(low_coeff)
      if (allocated(high_coeff))     deallocate(high_coeff)
      if (allocated(troe_coeff))     deallocate(troe_coeff)
      if (allocated(low_coeff_abs_max))      deallocate(low_coeff_abs_max)
      if (allocated(high_coeff_abs_max))     deallocate(high_coeff_abs_max)
      if (allocated(troe_coeff_abs_max))     deallocate(troe_coeff_abs_max)
      if (allocated(a_k4))           deallocate(a_k4)
      if (allocated(a_k4_min))       deallocate(a_k4_min)
      if (allocated(Mplus_case))     deallocate(Mplus_case)
      if (allocated(photochem_case)) deallocate(photochem_case)
      if (allocated(net_react_m))    deallocate(net_react_m)
      if (allocated(net_react_p))    deallocate(net_react_p)
      if (allocated(back))           deallocate(back)
      if (allocated(Diff_full))      deallocate(Diff_full)
      if (allocated(Diff_full_add))  deallocate(Diff_full_add)
!
    endsubroutine chemistry_clean_up
!***********************************************************************
    subroutine find_remove_real_stoic(Speciesstring,lreal,stoi,startindex)
!
      character(len=*), intent(in) :: Speciesstring
      logical, intent(inout) :: lreal
      real, intent(inout) :: stoi
      integer, intent(inout) :: startindex
      integer :: lreadable
!
      read (Speciesstring(1:3),'(F3.1)',iostat=lreadable) stoi
      if (lreadable == 0) then
        startindex = startindex+2
        lreal = .True.
      endif
!
    endsubroutine find_remove_real_stoic
!!***********************************************************************
    subroutine chemspec_normalization_N2(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: sum_Y !, sum_Y2
      integer :: k ,isN2, ichemsN2
      logical :: lsN2
!
      call find_species_index('N2', isN2, ichemsN2, lsN2)

      sum_Y=0.0 !; sum_Y2=0.0
      do k=1,nchemspec
        if (k/=ichemsN2) sum_Y=sum_Y+f(:,:,:,ichemspec(k))
      enddo
      f(:,:,:,isN2)=1.0-sum_Y
!
    endsubroutine chemspec_normalization_N2
!***********************************************************************
    subroutine getmu_array(f,mu1_full,linit)
!
!  Calculate mean molecular weight
!
!  16-mar-10/natalia
!  30-jun-17/MR: moved here from eos_chemistry.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: mu1_full
!
      integer :: k
      logical, optional :: linit
      logical :: llinit=.false.
!
      if (.not. present(linit)) then
        llinit=.false.
      else
        llinit=linit
      endif
!
!  Mean molecular weight
!
      mu1_full=0.
      do k=1,nchemspec
        if (species_constants(k,imass)>0.) then
          if (llinit) then
            mu1_full(:,:,:)= &
                 mu1_full(:,:,:)+unit_mass*&
                 f(:,:,:,ichemspec(k)) &
                 /species_constants(k,imass)
          else
            mu1_full(:,mm1:mm2,nn1:nn2)= &
                 mu1_full(:,mm1:mm2,nn1:nn2)+unit_mass*&
                 f(:,mm1:mm2,nn1:nn2,ichemspec(k)) &
                 /species_constants(k,imass)
          endif
        endif
      enddo
!
    endsubroutine getmu_array
!***********************************************************************
   subroutine read_Lewis
!
!  Reading of the species Lewis numbers in an input file
!
!  21-jun-10/julien: coded
!  30-jun-17/MR: moved here from eos_chemistry.
!
!
      logical :: emptyfile
      logical :: found_specie
      integer :: file_id=123, ind_glob, ind_chem, i
      real    :: lewisk
      character (len=10) :: specie_string
!
      emptyFile=.true.
!
      open(file_id,file="lewis.dat")
!
      i=0
      do
        read(file_id,*,end=1000) specie_string, lewisk
        emptyFile=.false.
!
        call find_species_index(specie_string,ind_glob,ind_chem,found_specie)
!
        if (found_specie) then
          if (lroot) print*,specie_string,' ind_glob=',ind_glob,' Lewis=', lewisk
          Lewis_coef(ind_chem) = lewisk
          Lewis_coef1(ind_chem) = 1./lewisk
          i=i+1
        endif

      enddo
!
! Stop if lewis.dat is empty
!
1000  if (emptyFile) call fatal_error('read_Lewis','end of file "lewis.dat"')
!
      if (i == 0) call warning('read_Lewis','File "lewis.dat" empty => Lewis numbers set to unity')
!
      close(file_id)
!
    endsubroutine read_Lewis
!***********************************************************************
    subroutine cond_spec_cond(f,df,p,ad,dustbin_width,mfluxcond)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: dustbin_width
      real, dimension (nx) :: mfluxcond
!
      real, dimension (nx) :: ff_cond_fact, Q_cond
      real, dimension(ndustspec) :: ad
      integer :: ichem, kkk,k
!
!
!  All the bins consume.
!
      p%ff_cond=0.
      p%part_heatcap=0.
      p%latent_heat=0.
      ff_cond_fact=4.*pi*mfluxcond*true_density_cond_spec
      do k=1,ndustspec
        p%ff_cond=p%ff_cond+ff_cond_fact*ad(k)**2*f(l1:l2,m,n,ind(k))*dustbin_width
        p%part_heatcap=p%part_heatcap+4*pi*ad(k)**3*f(l1:l2,m,n,ind(k))*true_density_cond_spec*dustbin_width
      enddo
      do ichem = 1,nchemspec
        kkk=ichemspec(ichem)
        if (kkk==i_cond_spec) then
          df(l1:l2,m,n,kkk) = df(l1:l2,m,n,kkk) + p%ff_cond*(f(l1:l2,m,n,kkk)-1.)/p%rho
        else
          df(l1:l2,m,n,kkk) = df(l1:l2,m,n,kkk) + p%ff_cond*f(l1:l2,m,n,kkk)/p%rho
        endif
      enddo
      if (ldensity_nolog) then
        df(l1:l2,m,n,irho) = df(l1:l2,m,n,irho) - p%ff_cond
      else
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) - p%ff_cond*p%rho1
      endif
      !
      ! Make pencil containing the latent heat due to condesation phase change
      ! This may be added to the energy equation of the gas phase, but this
      ! is currently not implemented for the eulerian approach. See the
      ! lagrangian approach (cond_spec_cond_lagr) below for info of how it
      ! is done.
      !
      if (.not. lnolatentheat) then
        p%latent_heat=p%latent_heat+p%ff_cond*deltaH_cgs/unit_temperature/molar_mass_spec*unit_mass
      endif
! 
    end subroutine cond_spec_cond
!***********************************************************************
    subroutine cond_spec_cond_lagr(f,df,p,rp,ix0,ix,np_swarm,dapdt)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      integer :: ichem, kkk,ix0,ix
      real :: ffcondp, Q_cond
      real, intent(IN) :: rp,dapdt,np_swarm
!
! Modify continuity equation
!
      ffcondp=dapdt*4*pi*rp**2*true_density_cond_spec*np_swarm
      p%ff_cond(ix) = p%ff_cond(ix)+ffcondp
      if (ldensity_nolog) then
        df(ix0,m,n,irho)   = df(ix0,m,n,irho)   - ffcondp
      else
        df(ix0,m,n,ilnrho) = df(ix0,m,n,ilnrho) - ffcondp*p%rho1(ix)
      endif
!
! Loop over all species and modify the species equation
!
      do ichem = 1,nchemspec
        kkk=ichemspec(ichem)
        if (kkk==i_cond_spec) then
          df(ix0,m,n,kkk) = df(ix0,m,n,kkk) + ffcondp*(f(ix0,m,n,kkk)-1.)*p%rho1(ix)
        else
          df(ix0,m,n,kkk) = df(ix0,m,n,kkk) + ffcondp*f(ix0,m,n,kkk)*p%rho1(ix)
        endif
      enddo
      !
      ! Make pencil containing the latent heat due to condesation phase change
      ! This may be added to the energy equation of the gas phase in
      ! particles_temperature.f90.
      !
      if (.not. lnolatentheat) then
        p%latent_heat(ix)=p%latent_heat(ix)+ffcondp*deltaH_cgs/unit_temperature/molar_mass_spec*unit_mass
      endif
!
    end subroutine cond_spec_cond_lagr
!***********************************************************************
    subroutine cond_spec_nucl(f,df,p,kk_vec,ad)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: ichem, kkk,i
      integer, dimension(nx) :: kk_vec
      real, dimension(ndustspec) :: ad
      real :: Q_nucl
!
      do i=1,nx
        !
        ! Check the radius of the nuclii fits within the allocated radius bins
        !
        if (p%nucl_rmin(i)>ad(1) .and. p%nucl_rmin(i)<ad(ndustspec)) then
          !
          !  Generating the nucleii consumes the condensing species
          !
          p%ff_nucl(i)=p%nucl_rate(i)*4.*pi*true_density_cond_spec/3.*ad(kk_vec(i))**3
          do ichem = 1,nchemspec
            kkk=ichemspec(ichem)
            if (kkk==i_cond_spec) then
              df(l1+i-1,m,n,kkk) = df(l1+i-1,m,n,kkk) + p%ff_nucl(i)*(f(l1+i-1,m,n,kkk)-1.)/p%rho(i)
            else
              df(l1+i-1,m,n,kkk) = df(l1+i-1,m,n,kkk) + p%ff_nucl(i)*f(l1+i-1,m,n,kkk)/p%rho(i)
            endif
          enddo
          if (ldensity_nolog) then
            df(l1+i-1,m,n,irho) = df(l1+i-1,m,n,irho) - p%ff_nucl(i)
          else
            df(l1+i-1,m,n,ilnrho) = df(l1+i-1,m,n,ilnrho) - p%ff_nucl(i)*p%rho1(i)
          endif
          !
          ! Add heat due to condesation phase change to the energy equation
          !
          if (.not. lnolatentheat) then
            p%latent_heat(i)=p%latent_heat(i)+p%ff_nucl(i)*deltaH_cgs/unit_temperature/molar_mass_spec*unit_mass
          endif
        else
          p%nucl_rate(i)=0.
        endif
      enddo
!
    end subroutine cond_spec_nucl
!***********************************************************************
    subroutine cond_spec_nucl_lagr(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: ichem, kkk,i,ii
      real, dimension(nx) :: Q_nucl
      !
      if (lnucleation .and. dt>0) then
        !
        ! Calculate mass flux of nucleii
        !
        p%ff_nucl=p%nucl_rate*4.*pi*true_density_cond_spec*p%nucl_rmin**3/3.
        where (p%ff_nucl<0)
          p%ff_nucl=0.0
        end where
        !
        ! The mass of the nucleii is added to the passive scalar equation
        !
        if (lnucl_dynamic) then
          df(l1:l2,m,n,icc+2) = df(l1:l2,m,n,icc+2) &
            + (p%ff_nucl*(1+f(l1:l2,m,n,icc))*p%rho1 - f(l1:l2,m,n,icc+2))/(Ntau*dt)
          df(l1:l2,m,n,icc) = df(l1:l2,m,n,icc) + f(l1:l2,m,n,icc+2)
          !
          ! Check if nucleii have just be added in particles_dust
          !
          do ii=l1,l2
            if (lnucleii_generated(ii-nghost,m-nghost,n-nghost)) then
              df(ii,m,n,icc+2) = df(ii,m,n,icc+2) &
                   - redfrac*f(ii,m,n,icc)/(Ntau*dt**2)
            endif
          enddo
        else
          df(l1:l2,m,n,icc) = df(l1:l2,m,n,icc) + p%ff_nucl*(1+p%cc(:,1))*p%rho1
        endif
        df(l1:l2,m,n,icc+1) = df(l1:l2,m,n,icc+1) + p%nucl_rmin*p%ff_nucl*(1+p%cc(:,2))*p%rho1
        !
        !  Generating the nucleii consumes the condensing species
        !
        do ichem = 1,nchemspec
          kkk=ichemspec(ichem)
          if (kkk==i_cond_spec) then
            df(l1:l2,m,n,kkk) = df(l1:l2,m,n,kkk) + p%ff_nucl*(f(l1:l2,m,n,kkk)-1.)*p%rho1
          else
            df(l1:l2,m,n,kkk) = df(l1:l2,m,n,kkk) + p%ff_nucl*f(l1:l2,m,n,kkk)*p%rho1
          endif
        enddo
        !
        ! The new nucleii also means that the density is reduced
        !
        if (ldensity_nolog) then
          df(l1:l2,m,n,irho) = df(l1:l2,m,n,irho) - p%ff_nucl
        else
          df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) - p%rho1*p%ff_nucl
        endif
        !
        ! Add latent heat due to condesation phase change to the energy equation.
        ! This is done in particles_temperature.f90.
        !
        if (.not. lnolatentheat) then
          p%latent_heat=p%latent_heat+p%ff_nucl*deltaH_cgs/unit_temperature/molar_mass_spec*unit_mass
       endif
        !
      endif
!
    end subroutine cond_spec_nucl_lagr
!***********************************************************************
    subroutine condensing_species_rate(p,mfluxcond)
      !
      real, dimension (nx) :: mfluxcond
      type (pencil_case) :: p
      !
      integer :: ichem, kkk
      !
      mfluxcond=A_spec*(p%chem_conc(:,ichem_cond_spec)-p%conc_sat_spec)*sqrt(p%TT)
      !
      ! Do we allow for evaporation of particles
      !
      if (lnoevap) then
        mfluxcond=max(mfluxcond,0.)
      endif
      !
    end subroutine condensing_species_rate
!***********************************************************************
    subroutine cond_spec_nucl_rate(p,nucleation_rmin,nucleation_rate)
      !
      type (pencil_case), intent(in) :: p
      real, dimension (nx), intent(out) :: nucleation_rate, nucleation_rmin
      !
      integer :: ichem, kkk, i
      real :: volume_spec_cgs=4.5e-23
      real :: volume_spec
      real, dimension (nx) :: sat_ratio_spec
      real, dimension (nx) :: gam_surf_energy, chem_conc
      real :: tmp1, tmp2, tmp3, tmp4, nucleation_rate_coeff
      !
!  compute rmin
!
      if (lnucleation) then
        !
        ! Set minimum concentration of the condensing species to avoid NaNs
        !
        chem_conc=max(1e-20,p%chem_conc(:,ichem_cond_spec))
        sat_ratio_spec=chem_conc/p%conc_sat_spec

        if (isurf_energy == "const") then
          gam_surf_energy=gam_surf_energy_cgs/(unit_mass/unit_time)
        elseif (isurf_energy == "Kingery") then
          gam_surf_energy=(307.+0.031*(p%TT-2073.))/unit_energy*unit_length**2*gam_surf_energy_mul_fac
        else
          call fatal_error("cond_spec_nucl_rate","No such isurf_energy")
        endif
        volume_spec=volume_spec_cgs/unit_length**3
        do i=1,nx
          if (sat_ratio_spec(i) > 1.0) then
            nucleation_rmin(i)=2.*gam_surf_energy(i)*volume_spec/&
                 (k_B*p%TT(i)*alog(sat_ratio_spec(i)))
            !
            !  Compute nucleation rate (of nucleii with r=rmin)
            !
            if (inucl_pre_exp .eq. "const") then
              nucleation_rate_coeff=nucleation_rate_coeff_cgs*unit_length**3
            elseif (inucl_pre_exp .eq. "oxtoby") then
              tmp1=sqrt(2*gam_surf_energy(i)/(pi*atomic_m_spec))
              tmp2=(chem_conc(i)*N_avogadro_cgs)**2
              nucleation_rate_coeff=tmp1*volume_spec*tmp2/sat_ratio_spec(i)
            elseif (inucl_pre_exp .eq. "becker_doring") then
              tmp1=sqrt(2*gam_surf_energy(i)/(pi*atomic_m_spec))
              tmp2=(chem_conc(i)*N_avogadro_cgs)**2
              nucleation_rate_coeff=tmp1*volume_spec*tmp2
            endif
            tmp3=-16.*pi*gam_surf_energy(i)**3*volume_spec**2
            tmp4=3*(k_B*p%TT(i))**3*(alog(sat_ratio_spec(i)))**2
            nucleation_rate(i)=nucleation_rate_coeff*exp(tmp3/tmp4)
          else
            nucleation_rate(i)=0.0
            nucleation_rmin(i)=min_nucl_radius_cgs/unit_length
          endif
        enddo
      else
        nucleation_rate=0.0
        nucleation_rmin=0.0
      endif
!
      end subroutine cond_spec_nucl_rate
!***********************************************************************
      subroutine cond_spec_sat_conc(p,conc_sat_spec)
        !
        ! Calculate the saturation concentration of the condensing species
        !
        type (pencil_case) :: p
        !
        real, dimension (nx) :: conc_sat_spec, tmp1
        real :: P_boil_cgs = 1013250 ! 1bar in Ba (=0.1Pa)
        real :: T_boil_cgs = 2503 ! in K
        real :: P_boil
        !
        if (iconc_sat_spec=="const") then
          conc_sat_spec=conc_sat_spec_cgs*unit_length**3
        elseif (iconc_sat_spec=="Clausius") then
          P_boil=P_boil_cgs/unit_pressure
          tmp1=-deltaH_cgs/(unit_energy*Rgas)*(p%TT1-unit_temperature/T_boil_cgs)
          conc_sat_spec=P_boil*exp(tmp1)/(Rgas*p%TT)
        else
          call fatal_error("cond_spec_sat_conc","no such iconc_sat_spec")
        endif
        !
      end subroutine cond_spec_sat_conc
!***********************************************************************
    subroutine chemistry_init_reduc_pointers
!
! 7-feb-24/TP:  allocates memory needed for reductions
!
      if (allocated(net_react_m)) p_net_react_m => net_react_m
      if (allocated(net_react_p)) p_net_react_p => net_react_p
!
    endsubroutine chemistry_init_reduc_pointers
!***********************************************************************
    subroutine chemistry_diags_reductions
!
! 7-feb-24/TP:  diag_reductions for chemistry
!
      if (allocated(net_react_m)) p_net_react_m = p_net_react_m + net_react_m
      if (allocated(net_react_p)) p_net_react_p = p_net_react_p + net_react_p
!
    endsubroutine chemistry_diags_reductions
!***********************************************************************
    include 'chemistry_common.inc'
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General,  only: string_to_enum

    integer, parameter :: n_pars=150
    integer(KIND=ikind8), dimension(n_pars) :: p_par
    integer :: i

    call copy_addr(rgas,p_par(1))
    call copy_addr(lambda_const,p_par(2))
    call copy_addr(visc_const,p_par(3))
    call copy_addr(lcorr_vel,p_par(4)) ! bool
    call copy_addr(init_tt1,p_par(5))
    call copy_addr(init_pressure,p_par(6))
    call copy_addr(lfilter_strict,p_par(7)) ! bool
    call copy_addr(lreactions,p_par(8)) ! bool
    call copy_addr(ladvection,p_par(9)) ! bool
    call copy_addr(ldiffusion,p_par(10)) ! bool
    call copy_addr(lheatc_chemistry,p_par(11)) ! bool
    call copy_addr(ldiff_simple,p_par(12)) ! bool
    call copy_addr(lthcond_simple,p_par(13)) ! bool
    call copy_addr(lt_const,p_par(14)) ! bool
    call copy_addr(ldiff_fick,p_par(15)) ! bool
    call copy_addr(lflux_simple,p_par(16)) ! bool
    call copy_addr(ldiff_corr,p_par(17)) ! bool
    call copy_addr(lfilter,p_par(18)) ! bool
    call copy_addr(lkreactions_alpha,p_par(19)) ! bool
    call copy_addr(nreactions,p_par(20)) ! int
    call copy_addr(mreactions,p_par(21)) ! int
    call copy_addr(iadv,p_par(22)) ! int
    call copy_addr(ldamp_zone_for_nscbc,p_par(23)) ! bool
    call copy_addr(l1step_test,p_par(24)) ! bool
    call copy_addr(tc,p_par(25))
    call copy_addr(tinf,p_par(26))
    call copy_addr(beta,p_par(27))
    call copy_addr(chem_diff,p_par(28))
    call copy_addr(lcheminp,p_par(29)) ! bool
    call copy_addr(lchem_cdtc,p_par(30)) ! bool
    call copy_addr(lmobility,p_par(31)) ! bool
    call copy_addr(itemp1,p_par(32)) ! int
    call copy_addr(itemp2,p_par(33)) ! int
    call copy_addr(itemp3,p_par(34)) ! int
    call copy_addr(lcloud,p_par(35)) ! bool
    call copy_addr(index_h2o,p_par(36)) ! int
    call copy_addr(lchemistry_diag,p_par(37)) ! bool
    call copy_addr(ireac,p_par(38)) ! int
    call copy_addr(chem_diff_prefactor,p_par(42)) ! (nchemspec)
    call copy_addr(mobility,p_par(44)) ! (nchemspec)
    call copy_addr(species_constants,p_par(47)) ! (nchemspec) (18)
    if (allocated(kreactions_alpha)) call copy_addr(kreactions_alpha,p_par(49)) ! (mreactions)
    call copy_addr(kreactions_p,p_par(50)) ! (mreactions)
    call copy_addr(kreactions_z,p_par(51)) ! (mz) (mreactions)
    call copy_addr(kreactions_m,p_par(52)) ! (mreactions)
    call copy_addr(sijm,p_par(53)) ! (nchemspec) (mreactions)

    call copy_addr(rgas_unit_sys,p_par(54))
    call copy_addr(diff_coef_const,p_par(55))
    call copy_addr(sc_number,p_par(56))
    call copy_addr(lnucleation,p_par(57)) ! bool
    call copy_addr(lgradp_terms,p_par(58)) ! bool
    call copy_addr(global_phi,p_par(59))
    call copy_addr(ldiff_lewis,p_par(60)) ! bool
    call copy_addr(lew_exist,p_par(61)) ! bool
    call copy_addr(lsmag_heat_transport,p_par(62)) ! bool
    call copy_addr(lsmag_diffusion,p_par(63)) ! bool
    call copy_addr(ipr,p_par(64)) ! int
    call copy_addr(z_cloud,p_par(65))
    call copy_addr(pr_turb,p_par(66))
    call copy_addr(ichem_cond_spec,p_par(67)) ! int
    call copy_addr(gam_surf_energy_cgs,p_par(68))
    call copy_addr(nucleation_rate_coeff_cgs,p_par(69))
    call copy_addr(atomic_m_spec,p_par(70))
    call copy_addr(deltah_cgs,p_par(71))
    call copy_addr(gam_surf_energy_mul_fac,p_par(72))
    call copy_addr(conc_sat_spec_cgs,p_par(73))
    call copy_addr(min_nucl_radius_cgs,p_par(74))
    call copy_addr(lnolatentheat,p_par(75)) ! bool
    call copy_addr(latmchem,p_par(76)) ! bool
    call copy_addr(index_o2,p_par(77)) ! int
    call copy_addr(index_n2,p_par(78)) ! int
    call copy_addr(index_o2n2,p_par(79)) ! int
    call copy_addr(stoichio,p_par(84)) ! (nchemspec) (mreactions)
    call copy_addr(sijp,p_par(85)) ! (nchemspec) (mreactions)
    call copy_addr(sijm_,p_par(86)) ! (nchemspec) (mreactions)
    call copy_addr(sijp_,p_par(87)) ! (nchemspec) (mreactions)
    call copy_addr(back,p_par(88)) ! bool (mreactions)
    call copy_addr(initial_massfractions,p_par(89)) ! (nchemspec)
    call copy_addr(b_n,p_par(94)) ! (mreactions)
    call copy_addr(alpha_n,p_par(95)) ! (mreactions)
    call copy_addr(e_an,p_par(96)) ! (mreactions)
    call copy_addr(low_coeff,p_par(97)) ! (3) (nreactions)
    call copy_addr(high_coeff,p_par(98)) ! (3) (nreactions)
    call copy_addr(troe_coeff,p_par(99)) ! (3) (nreactions)
    call copy_addr(a_k4,p_par(100)) ! (nchemspec) (nreactions)
    if(allocated(mplus_case)) call copy_addr(mplus_case,p_par(101)) ! bool (nreactions)
    call copy_addr(lewis_coef1,p_par(102)) ! (nchemspec)
    call string_to_enum(enum_reac_rate_method,reac_rate_method)
    call copy_addr(enum_reac_rate_method,p_par(103)) ! int
    do i = 1,mreactions
        call string_to_enum(enum_reaction_name(i),reaction_name(i))
    enddo
    call copy_addr(enum_reaction_name,p_par(104)) ! int (mreactions)
    call string_to_enum(enum_iconc_sat_spec,iconc_sat_spec)
    call copy_addr(enum_iconc_sat_spec,p_par(105)) ! int
    call string_to_enum(enum_isurf_energy,isurf_energy)
    call copy_addr(enum_isurf_energy,p_par(106)) ! int
    call string_to_enum(enum_inucl_pre_exp,inucl_pre_exp)
    call copy_addr(enum_inucl_pre_exp,p_par(107)) ! int
    call copy_addr(low_coeff_abs_max,p_par(113)) ! (nreactions)
    call copy_addr(high_coeff_abs_max,p_par(114)) ! (nreactions)
    call copy_addr(troe_coeff_abs_max,p_par(115)) ! (nreactions)
    call copy_addr(a_k4_min,p_par(123)) ! (nreactions)
    call copy_addr(i_o2_glob,p_par(124)) ! int
    call copy_addr(ichem_o2,p_par(125)) ! int
    call copy_addr(i_c3h8_glob,p_par(126)) ! int
    call copy_addr(ichem_c3h8,p_par(127)) ! int
    call copy_addr(lo2,p_par(128)) ! bool
    call copy_addr(lc3h8,p_par(129)) ! bool
    call copy_addr(mo2,p_par(130))
    call copy_addr(mc3h8,p_par(131))

  endsubroutine pushpars2c
!***********************************************************************
  subroutine make_flame_index(f)
!
    use Sub, only: grad
!
! Calculate flame index and store in f-array.
!
  real, dimension (mx,my,mz,mfarray) :: f
  integer :: n_loc, m_loc, igrad1, igrad2, ind_chem_spec1, ind_chem_spec2
  integer :: iweight
  logical :: lspec1, lspec2, lweight
  real, dimension(nx) :: FI
  real, dimension(nx,3) :: grad1, grad2
!
  call find_species_index(flameind_spec1,igrad1,ind_chem_spec1,lspec1)
  if (.not. lspec1) then
    print*,"flameind_spec1=",flameind_spec1
    call fatal_error("make_flame_index","did not find spec1")
  endif
  call find_species_index(flameind_spec2,igrad2,ind_chem_spec2,lspec2)
  if (.not. lspec2) then
    print*,"flameind_spec2=",flameind_spec2
    call fatal_error("make_flame_index","did not find spec2")
  endif
  iweight=ireaci(ind_chem_spec2)
  lweight=.true.
!
! Loop over all pencils
!
  do n_loc=n1,n2
    do m_loc=m1,m2
      m=m_loc;n=n_loc
      call grad(f,igrad1,grad1)
      call grad(f,igrad2,grad2)
      f(l1:l2,m,n,iFlameInd)=(grad1(:,1)*grad2(:,1)+grad1(:,2)*grad2(:,2)+grad1(:,3)*grad2(:,3))/sqrt( &
           (grad1(:,1)**2+grad1(:,2)**2+grad1(:,3)**2)* &
           (grad2(:,1)**2+grad2(:,2)**2+grad2(:,3)**2)+tini)
      if (lweight) f(l1:l2,m,n,iFlameInd)=f(l1:l2,m,n,iFlameInd)*abs(f(l1:l2,m,n,iweight))
    enddo
  enddo
!
end subroutine make_flame_index
!***********************************************************************
subroutine make_mixture_fraction(f)
!
! Calculate Bilger mixture fraction and store in f-array. 
!
  real, dimension (mx,my,mz,mfarray) :: f
  integer :: n_loc, m_loc, ichem, imass_spec
  real, dimension(nx) :: var
!
! Select which element to base the mixture fraction on
!
  select case (mixture_fraction_element)
  case ("H")
    imass_spec=imassH
  case ("C")
    imass_spec=imassC
  case ("O")
    imass_spec=imassO
  case ("S")
     imass_spec=imassS
  case ("TI")
     imass_spec=imassTI
   case default
     print*,"mixture_fraction_element=",mixture_fraction_element
     call fatal_error("make_mixture_fraction","No such mixture_fraction_element")
  end select
  !
  ! Loop over all pencils
  !
  do n_loc=n1,n2
    do m_loc=m1,m2
      m=m_loc;n=n_loc
      var=0.
      do ichem=1,nchemspec
        var=var+species_constants(ichem,imass_spec)*f(l1:l2,m,n,ichemspec(ichem))/ &
             species_constants(ichem,imass)
      enddo
      f(l1:l2,m,n,iMixFrac)=var/mix_frac_IH
    enddo
  enddo
!
end subroutine make_mixture_fraction
!***********************************************************************
endmodule Chemistry
