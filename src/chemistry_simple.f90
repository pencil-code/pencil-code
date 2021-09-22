! $Id$
!
!  This modules addes chemical species and reactions.
!  The units used in the chem.in files are cm3,mole,sec,kcal and K
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
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED cv; cv1; cp; cp1; glncp(3)
! PENCILS PROVIDED nu; gradnu(3)
! PENCILS PROVIDED DYDt_reac(nchemspec); DYDt_diff(nchemspec)
! PENCILS PROVIDED lambda; glambda(3)
! PENCILS PROVIDED Diff_penc_add(nchemspec); H0_RT(nchemspec); hhk_full(nchemspec)
! PENCILS PROVIDED ghhk(3,nchemspec); S0_R(nchemspec)

!
!***************************************************************
module Chemistry
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use EquationOfState
  use Messages, only: svn_id, timing, fatal_error, inevitably_fatal_error
  use Mpicomm, only: stop_it
!
  implicit none
!
  include 'chemistry.h'
!
  real :: Rgas, Rgas_unit_sys=1.
  real, pointer :: scale_Rgas
  real, dimension(mx,my,mz,nchemspec) :: cp_R_spec
! parameters for simplified cases
  logical :: init_from_file, reinitialize_chemistry=.false.
  character(len=30) :: reac_rate_method = 'chemkin'
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
!
  logical :: lfilter_strict=.false.
!
!  parameters related to chemical reactions, diffusion and advection
!
  logical :: lreactions=.true.
  logical :: ladvection=.true.
  logical :: ldiffusion=.true.
!
  logical :: lheatc_chemistry=.true.
  logical :: lDiff_simple=.true.
  logical :: lThCond_simple=.false.
  logical :: lT_const=.false.
  logical :: ldiff_corr=.false.
  logical, save :: tran_exist=.false.
  logical, save :: lew_exist=.false.
  logical :: lmech_simple=.false.
!
  logical :: lfilter=.false.
  integer :: nreactions=0, nreactions1=0, nreactions2=0
  integer :: ll1, ll2, mm1, mm2, nn1, nn2
!
  integer :: mreactions
!
!  The stociometric factors need to be reals for arbitrary reaction orders
!
  real, allocatable, dimension(:,:) :: stoichio, Sijm, Sijp
  logical, allocatable, dimension(:) :: back
  character(len=30), allocatable, dimension(:) :: reaction_name
  logical :: lT_tanh=.false.
!  logical :: ldamp_zone_for_NSCBC=.false.
  logical :: linit_temperature=.false., linit_density=.false.
  logical :: lreac_as_aux=.false.
!
  logical :: lflame_front=.false., lflame_front_2D=.false.
!
!  hydro-related parameters
!
  real, dimension(nchemspec) :: amplchemk=0., amplchemk2=0.
  real, dimension(nchemspec) :: chem_diff_prefactor=1., initial_massfractions
  real :: amplchem=1., kx_chem=1., ky_chem=1., kz_chem=1., widthchem=1.
  real :: chem_diff=0.
  character(len=labellen), dimension(ninit) :: initchem='nothing'
  character(len=60) :: prerun_directory='nothing'
  character(len=60) :: file_name='nothing'
!
  real, dimension(mx,my,mz,nchemspec), save :: RHS_Y_full
  real, dimension(nchemspec) :: nu_spec=0., mobility=1.
!
!  Chemkin related parameters
!
  logical :: lcheminp=.false., lchem_cdtc=.false.
  logical :: lmobility=.false.
  ! real, dimension(nchemspec,18) :: species_constants
  integer :: iTemp1=2, iTemp2=3, iTemp3=4
  integer, dimension(7) :: iaa1, iaa2
  real, allocatable, dimension(:) :: B_n, alpha_n, E_an
  real, allocatable, dimension(:,:) :: low_coeff, high_coeff, troe_coeff, a_k4
  logical, allocatable, dimension(:) :: Mplus_case
  logical, allocatable, dimension(:) :: photochem_case
  real, allocatable, dimension(:,:) :: orders_m, orders_p
  real :: lamb_low, lamb_up, Pr_turb=0.7
  real :: intro_time=0
  real :: p_init=1013000 ! cgs !
  integer :: i_CO, ichem_CO, i_CO2, ichem_CO2, ichem_O2
  logical :: lCO, lCO2, lO2
  ! const mass fraction of H2O for lmech_simple when only 4 species are used
  real :: Y_H2O=0.0008
  real :: m_H2O = 2.*1.00794+15.9994
  logical :: lback=.true.
  real :: scale_homo = 0.
!
!   Lewis coefficients
!
 real, dimension(nchemspec) :: Lewis_coef=1., Lewis_coef1=1.
!
!   Transport data
!
 real, dimension(nchemspec,7) :: tran_data
!
!   Diagnostics
!
  real, allocatable, dimension(:,:) :: net_react_m, net_react_p
  logical :: lchemistry_diag=.false.
!
! input parameters
  namelist /chemistry_init_pars/ &
      initchem, amplchem, kx_chem, ky_chem, kz_chem, widthchem, &
      amplchemk,amplchemk2, chem_diff,nu_spec, &
      lThCond_simple,&
      init_x1,init_x2,init_y1,init_y2,init_z1,init_z2,init_TT1,&
      init_TT2,init_rho,&
      init_ux,init_uy,init_uz,init_pressure, &
      str_thick,lT_tanh,lT_const,lheatc_chemistry, &
      prerun_directory, lback, scale_homo, &
      lchemistry_diag,lfilter_strict,linit_temperature, &
      linit_density, init_rho2, intro_time, p_init, Y_H2O, &
      file_name, lreac_as_aux, init_zz1, init_zz2, flame_pos, &
      reac_rate_method,global_phi, Pr_turb, lew_exist, Lewis_coef, lmech_simple 
!
!
! run parameters
  namelist /chemistry_run_pars/ &
      chem_diff,chem_diff_prefactor, nu_spec, ldiffusion, ladvection, &
      lreactions,lchem_cdtc,lheatc_chemistry, lchemistry_diag, &
      lmobility,mobility, lfilter,lT_tanh, &
      lThCond_simple,reinitialize_chemistry,init_from_file, &
      lfilter_strict,init_TT1,init_TT2,init_x1,init_x2, linit_temperature, &
      linit_density, intro_time, &
      ldiff_corr, lreac_as_aux, reac_rate_method,global_phi
!
! diagnostic variables (need to be consistent with reset list below)
!
  integer, dimension(nchemspec) :: idiag_Ym=0      ! DIAG_DOC: $\left<Y_x\right>$
  integer, dimension(nchemspec) :: idiag_dYm=0     ! DIAG_DOC: $\delta\left<Y_x\right>/\delta t$
  integer, dimension(nchemspec) :: idiag_dYmax=0   ! DIAG_DOC: $max\delta\left<Y_x\right>/\delta t$
  integer, dimension(nchemspec) :: idiag_Ymax=0    ! DIAG_DOC: $\left<Y_{x,max}\right>$
  integer, dimension(nchemspec) :: idiag_Ymin=0    ! DIAG_DOC: $\left<Y_{x,min}\right>$
  integer, dimension(nchemspec) :: idiag_hm=0      ! DIAG_DOC: $\left<H_{x,max}\right>$
  integer, dimension(nchemspec) :: idiag_cpm=0     ! DIAG_DOC: $\left<c_{p,x}\right>$
  integer, dimension(nchemspec) :: idiag_diffm=0   ! DIAG_DOC: $\left<D_{x}\right>$
  integer, dimension(nchemspec) :: idiag_Ymz=0     ! DIAG_DOC: $\left<Y_x\right>_{xy}(z)$
  integer :: idiag_dtchem=0     ! DIAG_DOC: $dt_{chem}$
!
!
  integer :: idiag_e_intm=0
!
  integer :: idiag_lambdam=0
  integer :: idiag_num=0
!
!  Auxiliaries.
!
  integer :: ireac=0, ireac_CO2=0, ireac_CO=0, ireac_O2=0
!

  contains
!
!***********************************************************************
    subroutine register_chemistry()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  13-aug-07/steveb: coded
!   8-jan-08/axel: added modifications analogously to dustdensity
!   5-mar-08/nils: Read thermodynamical data from chem.inp
!
      use FArrayManager
!
      integer :: k, ichemspec_tmp
      character(len=fnlen) :: input_file
      logical ::  chemin, cheminp
!
!  Initialize some index pointers
!
      iaa1(1) = 5
      iaa1(2) = 6
      iaa1(3) = 7
      iaa1(4) = 8
      iaa1(5) = 9
      iaa1(6) = 10
      iaa1(7) = 11
!
      iaa2(1) = 12
      iaa2(2) = 13
      iaa2(3) = 14
      iaa2(4) = 15
      iaa2(5) = 16
      iaa2(6) = 17
      iaa2(7) = 18
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
      call farray_register_auxiliary('viscosity',iviscosity)
!
!  Register cp
!
      call farray_register_auxiliary('cp',icp)
!
!  Register individual gas constant
!
      call farray_register_auxiliary('RR',iRR)
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
      if (lcheminp) then
        call read_species(input_file)
      else
!        if (lmech_simple) then
        if (lmech_simple .and. nchemspec==5) then
          varname(ichemspec(1):ichemspec(nchemspec))= (/ 'CO2       ','CO        ','N2        ','O2        ','H2O       '/)
        elseif (lmech_simple) then
          varname(ichemspec(1):ichemspec(nchemspec))= (/ 'CO2       ','CO        ','N2        ','O2        '/)
        else
          varname(ichemspec(1):ichemspec(nchemspec))= (/ 'H2        ','O2        ','H2O       ','H         ','O         ',&
             'OH        ','HO2       ','H2O2      ','AR        ','N2        ','HE        ','CO        ','CO2       '/)
        endif
      endif
!
!  Read data on the thermodynamical properties of the different species.
!  All these data are stored in the array species_constants.
!
      if (lcheminp) then
        call read_thermodyn(input_file)
      else 
        call read_thermodyn_simple()
      endif
!
!  Write all data on species and their thermodynamics to file.
!
      if (lcheminp) call write_thermodyn()
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( "$Id$")
!
    endsubroutine register_chemistry
!***********************************************************************
    subroutine read_thermodyn_simple()

!   Hard-coded H2_flamespeed mechanism
!   SPECIES:
!   1-H2 2-O2 3-H2O 4-H 5-O 6-OH 7-HO2 8-H2O2 9-AR 10-N2 11-HE 12-CO 13-CO2
!   
!   Mech_simple: 2CO+O2 = 2CO2
!   1-CO2 2-CO 3-N2 4-O2 5-H2O 

    integer, dimension(7) :: iaa1,iaa2
    integer :: iTemp1=2,iTemp2=3,iTemp3=4
!
    iaa1(1)=5;iaa1(2)=6;iaa1(3)=7;iaa1(4)=8
    iaa1(5)=9;iaa1(6)=10;iaa1(7)=11
!
    iaa2(1)=12;iaa2(2)=13;iaa2(3)=14;iaa2(4)=15
    iaa2(5)=16;iaa2(6)=17;iaa2(7)=18

    if (lmech_simple) then
      ! O2
      species_constants(4,iaa1(1):iaa2(7)) = (/ 3.69757819E+00, 6.13519689E-04,-1.25884199E-07, 1.77528148E-11, &
                                               -1.13643531E-15,-1.23393018E+03, 3.18916559E+00, 3.21293640E+00, &
                                                1.12748635E-03,-5.75615047E-07, 1.31387723E-09,-8.76855392E-13, &
                                               -1.00524902E+03, 6.03473759E+00 /)
      ! N2
      species_constants(3,iaa1(1):iaa2(7)) = (/ 0.02926640E+02, 0.01487977E-01,-0.05684761E-05, 0.01009704E-08, &
                                               -0.06753351E-13,-0.09227977E+04, 0.05980528E+02, 0.03298677E+02, &
                                                0.01408240E-01,-0.03963222E-04, 0.05641515E-07,-0.02444855E-10, &
                                               -0.01020900E+05, 0.03950372E+02 /)
      ! CO
      species_constants(2,iaa1(1):iaa2(7)) = (/0.03025078E+02, 0.01442689E-01,-0.05630828E-05, 0.01018581E-08, &
                                               -0.06910952E-13,-0.01426835E+06, 0.06108218E+02, 0.03262452E+02, &
                                                0.01511941E-01,-0.03881755E-04, 0.05581944E-07,-0.02474951E-10, &
                                               -0.01431054E+06, 0.04848897E+02 /)
      ! CO2
      species_constants(1,iaa1(1):iaa2(7)) = (/0.04453623E+02, 0.03140169E-01,-0.01278411E-04, 0.02393997E-08, &
                                               -0.01669033E-12,-0.04896696E+06,-0.09553959E+01, 0.02275725E+02, &
                                                0.09922072E-01,-0.01040911E-03, 0.06866687E-07,-0.02117280E-10, &
                                               -0.04837314E+06, 0.01018849E+03 /)
      if (nchemspec == 5) then
        !H2O
        species_constants(5,iaa1(1):iaa2(7)) = (/ 2.67214561E+00, 3.05629289E-03,-8.73026011E-07, 1.20099639E-10, &
                                                 -6.39161787E-15,-2.98992090E+04, 6.86281681E+00, 3.38684249E+00, &
                                                  3.47498246E-03,-6.35469633E-06, 6.96858127E-09,-2.50658847E-12, &
                                                 -3.02081133E+04, 2.59023285E+00 /)
        species_constants(1:5,imass) = (/12.0107+2.*15.9994, 12.0107+15.9994, 2.*14.00674, 2.*15.9994, 2.*1.00794+15.9994/)
      else
        species_constants(1:4,imass) = (/12.0107+2.*15.9994, 12.0107+15.9994, 2.*14.00674, 2.*15.9994/)
      endif
      species_constants(:,iTemp1) = 300.0
      species_constants(:,iTemp3) = 5000.0
      species_constants(:,iTemp2) = 1000.00
!
    else
!
      species_constants(5,iaa1(1):iaa2(7)) = (/ 2.54205966E+00,-2.75506191E-05,-3.10280335E-09, 4.55106742E-12, &
                                               -4.36805150E-16, 2.92308027E+04, 4.92030811E+00, 2.94642878E+00, &
                                               -1.63816649E-03, 2.42103170E-06,-1.60284319E-09, 3.89069636E-13, &
                                                2.91476445E+04, 2.96399498E+00 /)
      species_constants(4,iaa1(1):iaa2(7)) = (/ 2.50000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
                                                0.00000000E+00, 2.54716270E+04,-4.60117638E-01, 2.50000000E+00, &
                                                0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
                                                2.54716270E+04,-4.60117608E-01 /)
      species_constants(6,iaa1(1):iaa2(7)) = (/ 2.86472886E+00, 1.05650448E-03,-2.59082758E-07, 3.05218674E-11, &
                                               -1.33195876E-15, 3.68362875E+03, 5.70164073E+00, 4.12530561E+00, &
                                               -3.22544939E-03, 6.52764691E-06,-5.79853643E-09, 2.06237379E-12, &
                                                3.34630913E+03,-6.90432960E-01 /)
      species_constants(1,iaa1(1):iaa2(7)) = (/ 2.99142337E+00, 7.00064411E-04,-5.63382869E-08,-9.23157818E-12, &
                                                1.58275179E-15,-8.35033997E+02,-1.35511017E+00, 3.29812431E+00, &
                                                8.24944174E-04,-8.14301529E-07,-9.47543433E-11, 4.13487224E-13, &
                                               -1.01252087E+03,-3.29409409E+00 /)
      species_constants(2,iaa1(1):iaa2(7)) = (/ 3.69757819E+00, 6.13519689E-04,-1.25884199E-07, 1.77528148E-11, &
                                               -1.13643531E-15,-1.23393018E+03, 3.18916559E+00, 3.21293640E+00, &
                                                1.12748635E-03,-5.75615047E-07, 1.31387723E-09,-8.76855392E-13, &
                                               -1.00524902E+03, 6.03473759E+00 /)
      species_constants(3,iaa1(1):iaa2(7)) = (/ 2.67214561E+00, 3.05629289E-03,-8.73026011E-07, 1.20099639E-10, &
                                               -6.39161787E-15,-2.98992090E+04, 6.86281681E+00, 3.38684249E+00, &
                                                3.47498246E-03,-6.35469633E-06, 6.96858127E-09,-2.50658847E-12, &
                                               -3.02081133E+04, 2.59023285E+00 /)
      species_constants(7,iaa1(1):iaa2(7)) = (/ 4.01721090E+00, 2.23982013E-03,-6.33658150E-07, 1.14246370E-10, &
                                               -1.07908535E-14, 1.11856713E+02, 3.78510215E+00, 4.30179801E+00, &
                                               -4.74912051E-03, 2.11582891E-05,-2.42763894E-08, 9.29225124E-12, &
                                                2.94808040E+02, 3.71666245E+00 /)
      species_constants(8,iaa1(1):iaa2(7)) = (/ 4.57316685E+00, 4.33613639E-03,-1.47468882E-06, 2.34890357E-10, &
                                               -1.43165356E-14,-1.80069609E+04, 5.01136959E-01, 3.38875365E+00, &
                                                6.56922581E-03,-1.48501258E-07,-4.62580552E-09, 2.47151475E-12, &
                                               -1.76631465E+04, 6.78536320E+00 /)
      species_constants(9,iaa1(1):iaa2(7)) = (/ 0.02500000E+02, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
                                                0.00000000E+00,-0.07453750E+04, 0.04366001E+02, 0.02500000E+02, &
                                                0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
                                               -0.07453750E+04, 0.04366001E+02 /)
      species_constants(10,iaa1(1):iaa2(7)) = (/0.02926640E+02, 0.01487977E-01,-0.05684761E-05, 0.01009704E-08, &
                                               -0.06753351E-13,-0.09227977E+04, 0.05980528E+02, 0.03298677E+02, &
                                                0.01408240E-01,-0.03963222E-04, 0.05641515E-07,-0.02444855E-10, &
                                               -0.01020900E+05, 0.03950372E+02 /)
      species_constants(11,iaa1(1):iaa2(7)) = (/0.02500000E+02, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
                                                0.00000000E+00,-0.07453750E+04, 0.09153489E+01, 0.02500000E+02, &
                                                0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, &
                                               -0.07453750E+04, 0.09153488E+01 /)
      species_constants(12,iaa1(1):iaa2(7)) = (/0.03025078E+02, 0.01442689E-01,-0.05630828E-05, 0.01018581E-08, &
                                               -0.06910952E-13,-0.01426835E+06, 0.06108218E+02, 0.03262452E+02, &
                                                0.01511941E-01,-0.03881755E-04, 0.05581944E-07,-0.02474951E-10, &
                                               -0.01431054E+06, 0.04848897E+02 /)
      species_constants(13,iaa1(1):iaa2(7)) = (/0.04453623E+02, 0.03140169E-01,-0.01278411E-04, 0.02393997E-08, &
                                               -0.01669033E-12,-0.04896696E+06,-0.09553959E+01, 0.02275725E+02, &
                                                0.09922072E-01,-0.01040911E-03, 0.06866687E-07,-0.02117280E-10, &
                                               -0.04837314E+06, 0.01018849E+03 /)

      species_constants(1:13,imass) = (/ 2.*1.00794, 2.*15.9994, 2.*1.00794+15.9994, 1.00794, &
                                         15.9994, 1.00794+15.9994, 1.00794+2.*15.9994, 2.*1.00794+2.*15.9994, &
                                         39.948, 2.*14.00674, 4.0026, 12.0107+15.9994, 12.0107+2.*15.9994 /)

      species_constants(1:13,iTemp1) = (/ 300.00, 300.00, 300.00, 300.00, 300.00, 200.00, 200.00, &
                                          300.00, 300.00, 300.00, 300.00, 300.00, 300.00 /)
      species_constants(1:13,iTemp3) = (/ 5000.00, 5000.00, 5000.00, 5000.00, 5000.00, 6000.00, 3500.00, &
                                        5000.00, 5000.00, 5000.00, 5000.00, 5000.00, 5000.00 /)
      species_constants(1:13,iTemp2) = 1000.00

    endif

    endsubroutine read_thermodyn_simple
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
      use SharedVariables, only: get_shared_variable, put_shared_variable
      use Messages, only: warning
!
      real, dimension(mx,my,mz,mfarray) :: f
      logical :: data_file_exit=.false.
      logical :: exist, exist1, exist2
      character(len=15) :: file1='chemistry_m.dat', file2='chemistry_p.dat'
      integer :: ind_glob, ind_chem, stat
      integer :: i
      logical :: found_specie=.false.
      character(len=2) :: car2
!
!  initialize chemistry
!
        if (unit_temperature /= 1) then
          call fatal_error('initialize_chemistry', &
              'unit_temperature must be unity when using chemistry!')
        endif
!
!  calculate universal gas constant based on Boltzmann constant
!  and the proton mass
!
        call get_shared_variable('scale_Rgas', scale_Rgas)
        if (unit_system == 'cgs') then
          Rgas_unit_sys = k_B_cgs/m_u_cgs
          Rgas = Rgas_unit_sys/unit_energy*scale_Rgas
        endif
!
!  Read in data file in ChemKin format
!
      if (lcheminp) then
        call chemkin_data
        data_file_exit = .true.
      else
        call chemkin_data_simple
        data_file_exit = .true.
      endif
!
! check the existence of a data file
!
      if (.not. data_file_exit) then
        call stop_it('initialize_chemistry: there is no chemistry data file')
      endif
!
      if ((nxgrid == 1) .and. (nygrid == 1) .and. (nzgrid == 1)) then
        ll1 = 1
        ll2 = mx
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
!  allocate memory for net_reaction diagnostics
!
      if (allocated(net_react_p) .and. .not. lreloading) then
        print*, 'this should not be here'
      endif

      if (lchemistry_diag .and. .not. lreloading) then
        allocate(net_react_p(nchemspec,nreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for net_react_p")
        net_react_p = 0.
        allocate(net_react_m(nchemspec,nreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for net_react_m")
        net_react_m = 0.
      endif
!
!  Define the chemical reaction rates as auxiliary variables for output
!
      if (lreac_as_aux) then
        if (ireac == 0) then
          call find_species_index('CO',i_CO,ichem_CO,lCO)
          call find_species_index('CO2',i_CO2,ichem_CO2,lCO2)
          call find_species_index('O2',i_O2,ichem_O2,lO2)
!          call farray_register_auxiliary('reac_CO2',ireac_CO2)
          call farray_register_auxiliary('reac_CO',ireac_CO)
!          call farray_register_auxiliary('reac_O2',ireac_O2)
!          call farray_register_auxiliary('reac',ireac)
!          call farray_register_auxiliary('reac',ireac)
        endif
      endif
!
      if (lew_exist) then 
        Lewis_coef1 = 1./Lewis_coef
      else
        print*,'Lewis numbers need to be read from start.in, no option to read from file'
        print*,'Set all Le = 1'
      endif
!
!  Needed by ogrid_chemistry 
!
   if (lsolid_cells) then
!
      call put_shared_variable('lheatc_chemistry', lheatc_chemistry)
      call put_shared_variable('ldiffusion',ldiffusion)
      call put_shared_variable('ldiff_corr',ldiff_corr)
      call put_shared_variable('lew_exist',lew_exist)
      call put_shared_variable('lcheminp',lcheminp)
      call put_shared_variable('lThCond_simple',lThCond_simple)
      call put_shared_variable('tran_exist',tran_exist)
      call put_shared_variable('tran_data',tran_data)
      call put_shared_variable('lt_const',lt_const)
      call put_shared_variable('ladvection',ladvection)
      call put_shared_variable('lfilter',lfilter)
      call put_shared_variable('lfilter_strict',lfilter_strict)
      call put_shared_variable('Lewis_coef1',Lewis_coef1)
      call put_shared_variable('nreactions',nreactions)
      call put_shared_variable('Sijm',Sijm)
      call put_shared_variable('Sijp',Sijp)
 !     call put_shared_variable('reaction_name',reaction_name)
      call put_shared_variable('back',back)
      call put_shared_variable('B_n',B_n)
      call put_shared_variable('alpha_n',alpha_n)
      call put_shared_variable('E_an',E_an)
      call put_shared_variable('low_coeff',low_coeff)
      call put_shared_variable('high_coeff',high_coeff)
      call put_shared_variable('troe_coeff',troe_coeff)
      call put_shared_variable('a_k4',a_k4)
      call put_shared_variable('Mplus_case',Mplus_case)
      call put_shared_variable('photochem_case',photochem_case)
      call put_shared_variable('stoichio',stoichio)
      call put_shared_variable('iaa1',iaa1)
      call put_shared_variable('iaa2',iaa2)
      call put_shared_variable('lmech_simple',lmech_simple)
      call put_shared_variable('species_constants',species_constants)
      call put_shared_variable('imass',imass)
      call put_shared_variable('lflame_front_2D',lflame_front_2D)
      call put_shared_variable('p_init',p_init)
      call put_shared_variable('lreac_as_aux',lreac_as_aux)
      if (lreac_as_aux) then
        call put_shared_variable('ireac',ireac)
        call put_shared_variable('ireac_CO',ireac_CO)
        call put_shared_variable('ireac_CO2',ireac_CO2)
        call put_shared_variable('ireac_O2',ireac_O2)
      endif
      if (lmech_simple) then
        call put_shared_variable('orders_p',orders_p)
        call put_shared_variable('orders_m',orders_m)
        call put_shared_variable('Y_H2O',Y_H2O)
        call put_shared_variable('m_H2O',m_H2O)
      endif
!
   endif 
!
!  write array dimension to chemistry diagnostics file
!
      open (1,file=trim(datadir)//'/net_reactions.dat',position='append')
      write (1,*) nchemspec,nreactions
      close (1)
!
    endsubroutine initialize_chemistry
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
      real :: PP
      integer :: j,k
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
            f(:,:,:,ichemspec(k)) = amplchemk(k)
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
            call hatwave(amplchem,f,ichemspec(k),widthchem,kx=kx_chem)
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
          if (.not. air_exist) then
            inquire (file='air.in',exist=air_exist)
          endif
          if (air_exist) then
            call air_field(f,PP)
          else
            call stop_it('there is no air.in/dat file')
          endif
        case ('flame_front')
          call flame_front(f)
        case ('flame_front_2D')
          call flame_front_2D(f)
        case default
!
!  Catch unknown values
!
          if (lroot) print*, 'initchem: No such value for initchem: ', &
              trim(initchem(j))
          call stop_it('')
!
        endselect
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_chemistry(f)
!
    endsubroutine init_chemistry
!***********************************************************************
    subroutine pencil_criteria_chemistry()
!
!  All pencils that this chemistry module depends on are specified here.
!
!  13-aug-07/steveb: coded
!
      lpenc_requested(i_ghhk) = .true.
      if (lreactions) lpenc_requested(i_DYDt_reac) = .true.
      lpenc_requested(i_DYDt_diff) = .true.
      lpenc_requested(i_H0_RT) = .true.
      lpenc_requested(i_S0_R) = .true.
      lpenc_requested(i_hhk_full) = .true.
      if (ldiffusion .or. lparticles_chemistry) then
        lpenc_requested(i_Diff_penc_add) = .true.
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
!      lpencil_in(i_lambda) = .true.
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
        lpencil_in(i_RR) = .true.
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
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_chemistry
!***********************************************************************
    subroutine calc_pencils_chemistry(f,p)
!
!  Calculate chemistry pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   13-aug-07/steveb: coded
!   10-jan-11/julien: adapted for the case where chemistry is solved by LSODE
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f
      intent(inout) :: p
      integer :: k,i,j2,j3
      integer :: ii1=1, ii2=2, ii3=3, ii4=4, ii5=5, ii6=6, ii7=7
      real :: T_low, T_up, T_mid
      real, dimension(nx) :: T_loc, TT_2, TT_3, TT_4, D_th
!
      T_loc = p%TT
      TT_2 = T_loc*T_loc
      TT_3 = TT_2*T_loc
      TT_4 = TT_3*T_loc
!
!  Dimensionless Standard-state molar enthalpy H0/RT
!
      if (lpencil(i_H0_RT)) then
        if ((.not. lT_const)) then
          do k = 1,nchemspec
            T_low = species_constants(k,iTemp1)
            T_mid = species_constants(k,iTemp2)
            T_up = species_constants(k,iTemp3)
            where (T_loc <= T_mid)
              p%H0_RT(:,k) = species_constants(k,iaa2(ii1)) &
                           + species_constants(k,iaa2(ii2))*T_loc/2 &
                           + species_constants(k,iaa2(ii3))*TT_2/3  &
                           + species_constants(k,iaa2(ii4))*TT_3/4  &
                           + species_constants(k,iaa2(ii5))*TT_4/5  &
                           + species_constants(k,iaa2(ii6))/T_loc
            elsewhere
              p%H0_RT(:,k) = species_constants(k,iaa1(ii1)) &
                           + species_constants(k,iaa1(ii2))*T_loc/2 &
                           + species_constants(k,iaa1(ii3))*TT_2/3  &
                           + species_constants(k,iaa1(ii4))*TT_3/4  &
                           + species_constants(k,iaa1(ii5))*TT_4/5  &
                           + species_constants(k,iaa1(ii6))/T_loc
            endwhere
          enddo
!
!  Enthalpy flux
!
          if (lpencil(i_hhk_full) ) then
            do k = 1,nchemspec
              if (species_constants(k,imass) > 0.)  then
                p%hhk_full(:,k) = p%H0_RT(:,k)*Rgas*T_loc &
                    /species_constants(k,imass)
              endif
            enddo
          endif
        endif
      endif
!
      if (lpencil(i_ghhk)) then
        do k = 1,nchemspec
          if (species_constants(k,imass) > 0.)  then
            do i = 1,3
!              p%ghhk(:,i,k) = p%cp*p%glnTT(:,i)*p%TT
              p%ghhk(:,i,k) = cp_R_spec(l1:l2,m,n,k)/species_constants(k,imass)*Rgas*p%glnTT(:,i)*p%TT
            enddo
          endif
        enddo
      endif
!
!  Find the entropy by using fifth order temperature fitting function
!
      if (lpencil(i_S0_R)) then
        if (.not. lT_const) then
          do k = 1,nchemspec
            T_low = species_constants(k,iTemp1)
            T_mid = species_constants(k,iTemp2)
            T_up = species_constants(k,iTemp3)
            where (T_loc <= T_mid .and. T_low <= T_loc)
              p%S0_R(:,k) = species_constants(k,iaa2(ii1))*p%lnTT &
                          + species_constants(k,iaa2(ii2))*T_loc  &
                          + species_constants(k,iaa2(ii3))*TT_2/2 &
                          + species_constants(k,iaa2(ii4))*TT_3/3 &
                          + species_constants(k,iaa2(ii5))*TT_4/4 &
                          + species_constants(k,iaa2(ii7))
            elsewhere (T_mid <= T_loc .and. T_loc <= T_up)
              p%S0_R(:,k) = species_constants(k,iaa1(ii1))*p%lnTT &
                          + species_constants(k,iaa1(ii2))*T_loc  &
                          + species_constants(k,iaa1(ii3))*TT_2/2 &
                          + species_constants(k,iaa1(ii4))*TT_3/3 &
                          + species_constants(k,iaa1(ii5))*TT_4/4 &
                          + species_constants(k,iaa1(ii7))
            endwhere
          enddo
        endif
      endif
!
! Calculate the reaction term and the corresponding pencil
!
      if (lreactions .and. lpencil(i_DYDt_reac)) then
        call calc_reaction_term(f,p)
      else
        p%DYDt_reac = 0.
      endif
!
! Calculate the diffusion term and the corresponding pencil
!
      if (lpencil(i_Diff_penc_add)) then
!
! D_th is computed twice in calc_penc in eos and here
!
        D_th = f(l1:l2,m,n,iviscosity)/Pr_number
        do k = 1,nchemspec
          p%Diff_penc_add(:,k) = D_th*Lewis_coef1(k)
        enddo
      endif
!
      if (ldiffusion .and. lpencil(i_DYDt_diff)) then
        call calc_diffusion_term(f,p)
      else
        p%DYDt_diff = 0.
      endif
!
      RHS_Y_full(l1:l2,m,n,:) = p%DYDt_reac+p%DYDt_diff
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
      real, dimension(mx,my,mz) :: mu1_full
      integer :: i, j,k
!
      real :: mO2=0., mH2=0., mN2=0., mH2O=0., mCH4=0., mCO2=0., mCO=0.
      real :: log_inlet_density, del, PP
      integer :: i_H2=0, i_O2=0, i_H2O=0, i_N2=0
      integer :: ichem_H2=0, ichem_O2=0, ichem_N2=0, ichem_H2O=0
      integer :: i_CH4=0, i_CO2=0, i_CO=0, ichem_CH4=0, ichem_CO2=0, ichem_CO=0
      real :: initial_mu1, final_massfrac_O2, final_massfrac_CH4, &
          final_massfrac_H2O, final_massfrac_CO2, final_massfrac_CO
      real :: init_H2, init_O2, init_N2, init_H2O, init_CO2, init_CH4, init_CO
      logical :: lH2=.false., lO2=.false., lN2=.false., lH2O=.false.
      logical :: lCH4=.false., lCO2=.false., lCO=.false.
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
      call find_species_index('CO',i_CO,ichem_CO,lCO)
      if (lCO) then
        mCO = species_constants(ichem_CO,imass)
        init_CO = initial_massfractions(ichem_CO)
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
        final_massfrac_O2 = &
            1. - final_massfrac_CO2 - final_massfrac_H2O  &
            - init_N2
      elseif (lCO .and. .not. lH2 .and. .not. lCH4) then
        final_massfrac_H2O = init_H2O
        final_massfrac_CO2 = mCO2/mCO * init_CO
        final_massfrac_O2 = &
            1. - final_massfrac_CO2 - final_massfrac_H2O  &
            - init_N2
      endif
!
      if (final_massfrac_O2 < 0.) final_massfrac_O2 = 0.
      if (lroot) then
        print*, '          init                      final'
        if (lH2 .and. .not. lCH4) print*, 'H2 :', init_H2, 0.
        if (lCH4) print*, 'CH4 :', init_CH4, 0.
        if (lO2) print*, 'O2 :', init_O2, final_massfrac_O2
        if (lH2O) print*, 'H2O :', 0., final_massfrac_H2O
        if (lCO2)  print*, 'CO2 :', 0., final_massfrac_CO2
        if (lCO .and. .not. lH2 .and. .not. lCH4)  print*, 'CO :', init_CO, 0.
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
              f(k,:,:,ilnTT) = log((x(k)-init_x1)/(init_x2-init_x1) &
              *(init_TT2-init_TT1)+init_TT1)
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
            if (lroot) print*, 'No tanh initial function available for CH4 combustion.'
          endif
!
        else
          if (x(k) > init_x1) then
            if (lH2 .and. .not. lCH4) f(k,:,:,i_H2) = init_H2* &
                (exp(f(k,:,:,ilnTT))-init_TT2)/(init_TT1-init_TT2)
            if (lCH4) f(k,:,:,i_CH4) = init_CH4*(exp(f(k,:,:,ilnTT))-init_TT2) &
                /(init_TT1-init_TT2)
            if (lCO .and. .not. lH2 .and. .not. lCH4) then
              f(k,:,:,i_CO) = init_CO*(exp(f(k,:,:,ilnTT))-init_TT2) &
                /(init_TT1-init_TT2)
            endif
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
              f(k,:,:,i_O2) = (x(k)-init_x1)/(init_x2-init_x1) &
              *(final_massfrac_O2-init_O2)+init_O2
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
            if (lroot) print*, 'No tanh initial function available for CH4 combustion.'
          endif
        enddo
!
      else
        do k = 1,mx
          if (x(k) >= init_x1 .and. x(k) < init_x2) then
            f(k,:,:,i_H2O) = (x(k)-init_x1)/(init_x2-init_x1) &
                *final_massfrac_H2O
            if (lCO .and. .not. lH2 .and. .not. lCH4) then
              f(k,:,:,i_H2O) = final_massfrac_H2O
            endif
            if (lCO2) f(k,:,:,i_CO2) = (x(k)-init_x1)/(init_x2-init_x1) &
                *final_massfrac_CO2
          elseif (x(k) >= init_x2) then
            if (lCO2) f(k,:,:,i_CO2) = final_massfrac_CO2
            if (lH2O) f(k,:,:,i_H2O) = final_massfrac_H2O
          endif
        enddo
      endif
!
      if (unit_system == 'cgs') then
        Rgas_unit_sys = k_B_cgs/m_u_cgs
        Rgas = Rgas_unit_sys/unit_energy*scale_Rgas
      endif
!
!  Find logaritm of density at inlet
!
      initial_mu1 = &
          initial_massfractions(ichem_O2)/(mO2) &
          +initial_massfractions(ichem_H2O)/(mH2O) &
          +initial_massfractions(ichem_N2)/(mN2)
      if (lH2 .and. .not. lCH4) initial_mu1 = initial_mu1+ &
          initial_massfractions(ichem_H2)/(mH2)
      if (lCO2) initial_mu1 = initial_mu1+init_CO2/(mCO2)
      if (lCO .and. .not. lH2 .and. .not. lCH4) initial_mu1 = initial_mu1+init_CO/(mCO)
      if (lCH4) initial_mu1 = initial_mu1+init_CH4/(mCH4)
      log_inlet_density = &
          log(init_pressure)-log(Rgas)-log(init_TT1)-log(initial_mu1)
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
      if (ldensity_nolog) &
          f(l1:l2,m1:m2,n1:n2,irho) = exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ltemperature_nolog) &
          f(:,:,:,iTT) = exp(f(:,:,:,ilnTT))
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
    subroutine flame_front_2D(f)
!
!   Flame_front adjusted to 2D. Sets T, rho and species field but not u.
!   U field is set in the ogrid module.
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz) :: mu1_full
      integer :: i, j,k
!
      real :: mO2=0., mH2=0., mN2=0., mH2O=0., mCH4=0., mCO2=0., mCO=0.
      real :: log_inlet_density, del, PP
      integer :: i_H2=0, i_O2=0, i_H2O=0, i_N2=0
      integer :: ichem_H2=0, ichem_O2=0, ichem_N2=0, ichem_H2O=0
      integer :: i_CH4=0, i_CO2=0, i_CO=0, ichem_CH4=0, ichem_CO2=0, ichem_CO=0
      real :: initial_mu1, final_massfrac_O2, final_massfrac_CH4, &
          final_massfrac_H2O, final_massfrac_CO2, final_massfrac_CO
      real :: init_H2, init_O2, init_N2, init_H2O, init_CO2, init_CH4, init_CO
      logical :: lH2=.false., lO2=.false., lN2=.false., lH2O=.false.
      logical :: lCH4=.false., lCO2=.false., lCO=.false.
!
      lflame_front_2D = .true.
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
      call find_species_index('CO',i_CO,ichem_CO,lCO)
      if (lCO) then
        mCO = species_constants(ichem_CO,imass)
        init_CO = initial_massfractions(ichem_CO)
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
        final_massfrac_O2 = &
            1. - final_massfrac_CO2 - final_massfrac_H2O  &
            - init_N2
      elseif (lCO .and. .not. lH2 .and. .not. lCH4) then
        final_massfrac_H2O = init_H2O
        final_massfrac_CO2 = mCO2/mCO * init_CO
        final_massfrac_O2 = &
            1. - final_massfrac_CO2 - final_massfrac_H2O  &
            - init_N2
      endif
!
      if (final_massfrac_O2 < 0.) final_massfrac_O2 = 0.
      if (lroot) then
        print*, '          init                      final'
        if (lH2 .and. .not. lCH4) print*, 'H2 :', init_H2, 0.
        if (lCH4) print*, 'CH4 :', init_CH4, 0.
        if (lO2) print*, 'O2 :', init_O2, final_massfrac_O2
        if (lH2O) print*, 'H2O :', 0., final_massfrac_H2O
        if (lCO2)  print*, 'CO2 :', 0., final_massfrac_CO2
        if (lCO .and. .not. lH2 .and. .not. lCH4)  print*, 'CO :', init_CO, 0.
      endif
!
!  Initialize temperature and species
!
      do k = 1,my
!
!  Initialize temperature
!
          if (y(k) <= init_x1) f(:,k,:,ilnTT) = log(init_TT1)
          if (y(k) >= init_x2) f(:,k,:,ilnTT) = log(init_TT2)
          if (y(k) > init_x1 .and. y(k) < init_x2) &
              f(:,k,:,ilnTT) = log((y(k)-init_x1)/(init_x2-init_x1) &
              *(init_TT2-init_TT1)+init_TT1)
!
!  Initialize fuel
!
          if (y(k) > init_x1) then
            if (lH2 .and. .not. lCH4) f(:,k,:,i_H2) = init_H2* &
                (exp(f(:,k,:,ilnTT))-init_TT2)/(init_TT1-init_TT2)
            if (lCH4) f(:,k,:,i_CH4) = init_CH4*(exp(f(:,k,:,ilnTT))-init_TT2) &
                /(init_TT1-init_TT2)
            if (lCO .and. .not. lH2 .and. .not. lCH4) then
              f(:,k,:,i_CO) = init_CO*(exp(f(:,k,:,ilnTT))-init_TT2) &
                /(init_TT1-init_TT2)
            endif
          endif
!
!  Initialize oxygen
!
          if (y(k) > init_x2) f(:,k,:,i_O2) = final_massfrac_O2
          if (y(k) > init_x1 .and. y(k) <= init_x2) &
              f(:,k,:,i_O2) = (y(k)-init_x1)/(init_x2-init_x1) &
              *(final_massfrac_O2-init_O2)+init_O2
      enddo
!
! Initialize products
!
        do k = 1,my
          if (y(k) >= init_x1 .and. y(k) < init_x2) then
            f(:,k,:,i_H2O) = (y(k)-init_x1)/(init_x2-init_x1) &
                *final_massfrac_H2O
            if (lCO .and. .not. lH2 .and. .not. lCH4) then
              f(:,k,:,i_H2O) = final_massfrac_H2O
            endif
            if (lCO2) f(:,k,:,i_CO2) = (y(k)-init_x1)/(init_x2-init_x1) &
                *final_massfrac_CO2
          elseif (y(k) >= init_x2) then
            if (lCO2) f(:,k,:,i_CO2) = final_massfrac_CO2
            if (lH2O) f(:,k,:,i_H2O) = final_massfrac_H2O
          endif
        enddo
!
      if (unit_system == 'cgs') then
        Rgas_unit_sys = k_B_cgs/m_u_cgs
        Rgas = Rgas_unit_sys/unit_energy*scale_Rgas
      endif
!
!  Find logaritm of density at inlet
!
      initial_mu1 = &
          initial_massfractions(ichem_O2)/(mO2) &
          +initial_massfractions(ichem_H2O)/(mH2O) &
          +initial_massfractions(ichem_N2)/(mN2)
      if (lH2 .and. .not. lCH4) initial_mu1 = initial_mu1+ &
          initial_massfractions(ichem_H2)/(mH2)
      if (lCO2) initial_mu1 = initial_mu1+init_CO2/(mCO2)
      if (lCO .and. .not. lH2 .and. .not. lCH4) initial_mu1 = initial_mu1+init_CO/(mCO)
      if (lCH4) initial_mu1 = initial_mu1+init_CH4/(mCH4)
      log_inlet_density = &
          log(init_pressure)-log(Rgas)-log(init_TT1)-log(initial_mu1)
!
!  Initialize density
!
      call getmu_array(f,mu1_full)
      f(l1:l2,m1:m2,n1:n2,ilnrho) = log(init_pressure)-log(Rgas)  &
          -f(l1:l2,m1:m2,n1:n2,ilnTT)-log(mu1_full(l1:l2,m1:m2,n1:n2))
!
!  Initialize velocity
!
      f(:,:,:,iuy) = 0
      f(:,:,:,iux) = 0
!
!  Check if we want nolog of density or nolog of temperature
!
      if (ldensity_nolog) &
          f(l1:l2,m1:m2,n1:n2,irho) = exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ltemperature_nolog) &
          f(:,:,:,iTT) = exp(f(:,:,:,ilnTT))
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
    endsubroutine flame_front_2D
!***********************************************************************
    subroutine calc_for_chem_mixture(f)
!
!  Calculate quantities for a mixture
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz) :: mu1_full
      real, dimension(mx) ::  nu_dyn
! lambda_Suth in [g/(cm*s*K^0.5)]
      real :: lambda_Suth = 1.5e-5, Suth_const = 200.
      integer :: j2,j3, k
      real, dimension(mx) :: T_loc, T_loc_2, T_loc_3, T_loc_4
      real :: T_up, T_mid, T_low
!
          call getmu_array(f,mu1_full)
          f(:,:,:,iRR) = mu1_full*Rgas
!
            f(:,:,:,icp) = 0.
            if (Cp_const < impossible) then
!
              f(:,:,:,icp) = Cp_const*mu1_full
!            
            else

!
              do j3 = nn1,nn2
                do j2 = mm1,mm2
                  if (ltemperature_nolog) then 
                    T_loc = (f(:,j2,j3,iTT))
                  else
                    T_loc = exp(f(:,j2,j3,ilnTT))
                  endif
                  T_loc_2 = T_loc*T_loc
                  T_loc_3 = T_loc_2*T_loc
                  T_loc_4 = T_loc_3*T_loc
                  do k = 1,nchemspec
                    if (species_constants(k,imass) > 0.) then
                      T_low = species_constants(k,iTemp1)-1.
                      T_mid = species_constants(k,iTemp2)
                      T_up = species_constants(k,iTemp3)
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
                      where (T_loc >= T_low .and. T_loc <= T_mid)
                        cp_R_spec(:,j2,j3,k) = species_constants(k,iaa2(1)) &
                          +species_constants(k,iaa2(2))*T_loc &
                          +species_constants(k,iaa2(3))*T_loc_2 &
                          +species_constants(k,iaa2(4))*T_loc_3 &
                          +species_constants(k,iaa2(5))*T_loc_4
                      elsewhere (T_loc >= T_mid .and. T_loc <= T_up)
                        cp_R_spec(:,j2,j3,k) = species_constants(k,iaa1(1)) &
                          +species_constants(k,iaa1(2))*T_loc &
                          +species_constants(k,iaa1(3))*T_loc_2 &
                          +species_constants(k,iaa1(4))*T_loc_3 &
                          +species_constants(k,iaa1(5))*T_loc_4
                      endwhere
!
! Check if the temperature are within bounds
!
                      if (maxval(T_loc) > T_up .or. minval(T_loc) < T_low) then
                        print*,'iproc=',iproc
                        print*,'TT_full(:,j2,j3)=',T_loc
                        print*,'j2,j3=',j2,j3
                        call inevitably_fatal_error('calc_for_chem_mixture', &
                        'TT_full(:,j2,j3) is outside range', .true.)
                      endif
!
! Find cp and cv for the mixture for the full domain
!
                      f(:,j2,j3,icp) = f(:,j2,j3,icp)+cp_R_spec(:,j2,j3,k)*Rgas/species_constants(k,imass) &
                                      *f(:,j2,j3,ichemspec(k))
                    endif
                  enddo
                enddo
              enddo
            endif

!  Viscosity of a mixture
!
      do j3 = nn1,nn2
        do j2 = mm1,mm2
          if (ltemperature_nolog) then
            nu_dyn = lambda_Suth*f(:,j2,j3,iTT)**(3./2.)/(Suth_const+f(:,j2,j3,iTT))
          else
            nu_dyn = lambda_Suth*exp(f(:,j2,j3,ilnTT))**(3./2.)/(Suth_const+exp(f(:,j2,j3,ilnTT)))
          endif
          if (ldensity_nolog) then
            f(:,j2,j3,iviscosity) = nu_dyn/f(:,j2,j3,irho)
          else
            f(:,j2,j3,iviscosity) = nu_dyn/exp(f(:,j2,j3,ilnrho))
          endif
        enddo
      enddo
!
    endsubroutine calc_for_chem_mixture
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
      use Diagnostics
      use Sub, only: grad,dot_mn
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx,3) :: gchemspec, dk_D, sum_diff=0.
      real, dimension(nx) :: ugchemspec, sum_DYDT, sum_dhhk=0.
      real, dimension(nx) :: sum_dk_ghk, dk_dhhk, sum_hhk_DYDt_reac
      type (pencil_case) :: p
      real, dimension(nx) :: RHS_T_full, diffus_chem
!
!  indices
!
      integer :: j, k,i
      integer :: i1=1, i2=2, i3=3, i4=4, i5=5, i6=6, i7=7, i8=8, i9=9, i10=10
      integer :: i11=11, i12=12, i13=13, i14=14, i15=15, i16=16, i17=17, i18=18, i19=19
      integer :: iz1=1, iz2=2, iz3=3, iz4=4, iz5=5, iz6=6, iz7=7, iz8=8, iz9=9, iz10=10
      integer :: iz11=11, iz12=12, iz13=13, iz14=14, iz15=15, iz16=16, iz17=17
      integer :: iz18=18, iz19=19
!
      intent(in) :: p,f
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      call timing('dchemistry_dt','entered',mnloop=.true.)
      if (headtt .or. ldebug) print*,'dchemistry_dt: SOLVE dchemistry_dt'
!
!  loop over all chemicals
!
      do k = 1,nchemspec
!
!  advection terms
!
        if (lhydro .and. ladvection) then
          call grad(f,ichemspec(k),gchemspec)
          call dot_mn(p%uu,gchemspec,ugchemspec)
          if (lmobility) ugchemspec = ugchemspec*mobility(k)
          df(l1:l2,m,n,ichemspec(k)) = df(l1:l2,m,n,ichemspec(k))-ugchemspec
        endif
!
!  diffusion operator
!
        if (ldiffusion) then
          df(l1:l2,m,n,ichemspec(k)) = df(l1:l2,m,n,ichemspec(k))+p%DYDt_diff(:,k)
        endif
!
!  chemical reactions:
!  multiply with stoichiometric matrix with reaction speed
!  d/dt(x_i) = S_ij v_j
!
        if (lreactions) then
          df(l1:l2,m,n,ichemspec(k)) = df(l1:l2,m,n,ichemspec(k))+ p%DYDt_reac(:,k)
        endif
!
!  Add filter for negative concentrations
!
        if (lfilter .and. .not. lfilter_strict) then
          do i = 1,mx
            if ((f(i,m,n,ichemspec(k))+df(i,m,n,ichemspec(k))*dt) < -1e-25 ) then
              df(i,m,n,ichemspec(k)) = -1e-25*dt
            endif
            if ((f(i,m,n,ichemspec(k))+df(i,m,n,ichemspec(k))*dt) > 1. ) then
              df(i,m,n,ichemspec(k)) = 1.*dt
            endif
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
            if ((f(i,m,n,ichemspec(k))+df(i,m,n,ichemspec(k))*dt) > 1. ) then
              df(i,m,n,ichemspec(k)) = 1.*dt
            endif
          enddo
        endif
!
      enddo
!
!  Modify RHS of temperature equation
!
      if (ldensity) then 
!
          sum_DYDt = 0.
          sum_hhk_DYDt_reac = 0.
          sum_dk_ghk = 0.
!
          do k = 1,nchemspec
            if (species_constants(k,imass) > 0.) then
              sum_DYDt = sum_DYDt+Rgas/species_constants(k,imass) &
                  *(p%DYDt_reac(:,k)+p%DYDt_diff(:,k))
              if (lreactions) then
                sum_hhk_DYDt_reac = sum_hhk_DYDt_reac-p%hhk_full(:,k)*p%DYDt_reac(:,k)
              endif
!
!  Sum over all species of diffusion terms
!
              if (ldiffusion) then
                call grad(f,ichemspec(k),gchemspec)
                do i = 1,3
                  dk_D(:,i) = gchemspec(:,i)*p%Diff_penc_add(:,k)
                enddo
                call dot_mn(dk_D,p%ghhk(:,:,k),dk_dhhk)
                sum_dk_ghk = sum_dk_ghk+dk_dhhk
                if (ldiff_corr) then
                  do i = 1,3
                    sum_diff(:,i) = sum_diff(:,i)+dk_D(:,i)
                  enddo
                endif
              endif
!
            endif
          enddo
!
! If the correction velocity is added
!
          if (ldiff_corr .and. ldiffusion) then
            do k = 1,nchemspec
              call dot_mn(sum_diff,p%ghhk(:,:,k),sum_dhhk)
              sum_dk_ghk(:) = sum_dk_ghk(:)-f(l1:l2,m,n,ichemspec(k))*sum_dhhk(:)
            enddo
          endif
!
        if (ltemperature_nolog) then
          RHS_T_full = p%cv1*((sum_DYDt(:)-p%RRmix*p%divu)*p%TT &
                  +sum_dk_ghk+sum_hhk_DYDt_reac)
        else
          RHS_T_full = (sum_DYDt(:)-p%RRmix*p%divu)*p%cv1 &
              +sum_dk_ghk*p%TT1(:)*p%cv1+sum_hhk_DYDt_reac*p%TT1(:)*p%cv1
        endif
!
        if (.not. lT_const) then
            df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + RHS_T_full
        endif
!
        if (lheatc_chemistry) call calc_heatcond_chemistry(f,df,p)
      endif
!
      if (lreac_as_aux) then
        if (llast) call get_reac_rate(f,p)
      endif
!
! The part below is commented because ldt = false all the time and it's
! never executed (llsode = F, lfirst = T/F depending on sub-timestep)
! That's because ldt = (dt==0.) in register.f90, and dt is set in run.in
! file. We'd probably want to compute dt dynamically?
!
! this damping zone is needed in a case of NSCBC
!
!      if (ldamp_zone_for_NSCBC) call damp_zone_for_NSCBC(f,df)
!
!  For the timestep calculation, need maximum diffusion
!
!      if (lfirst .and. ldt) then
!print*,'Why never enters here??***********************************************************'
!          diffus_chem=0.
!          do j = 1,nx
!            if (ldiffusion .and. .not. ldiff_simple) then
!
!--------------------------------------
!  This expression should be discussed
!--------------------------------------
!
!              diffus_chem(j) = diffus_chem(j)+ &
!                  maxval(Diff_full_add(l1+j-1,m,n,1:nchemspec))*dxyz_2(j)
!            else
!              diffus_chem(j) = 0.
!            endif
!          enddo
!        maxdiffus=max(maxdiffus,diffus_chem)
!      endif
!
! NB: it should be discussed
!
!      if (lfirst .and. ldt) then
!        if (lreactions .and.(.not. llsode)) then
!
!  calculate maximum of *relative* reaction rate if decaying,
!  or maximum of absolute rate, if growing.
!
!          if (lchem_cdtc) then
!            reac_chem = 0.
!            do k = 1,nchemspec
!              reac_chem = max(reac_chem, &
!                  abs(p%DYDt_reac(:,k)/max(f(l1:l2,m,n,ichemspec(k)),.001)))
!            enddo
!
!          else
!            reac_chem = 0.
!            !sum_reac_rate=0.
!            do k = 1,nchemspec
!              reac_chem = reac_chem+abs(p%DYDt_reac(:,k)/ &
!                  max(f(l1:l2,m,n,ichemspec(k)),0.001))
!              !sum_reac_rate=sum_reac_rate+p%DYDt_reac(:,k)
!            enddo
!            if (maxval(reac_chem) > 1e11) then
!              reac_chem = 1e11
!            endif
!          endif
!        endif
!      endif
!
!  Calculate diagnostic quantities
!
      call timing('dchemistry_dt','before ldiagnos',mnloop=.true.)
      if (ldiagnos) then
        if (idiag_dtchem /= 0) then
          call max_mn_name(reac_chem/cdtc,idiag_dtchem,l_dt=.true.)
        endif
!
        do i = 1,nchemspec
          if (idiag_Ym(i)/= 0) then
            call sum_mn_name(f(l1:l2,m,n,ichemspec(i)),idiag_Ym(i))
          endif
          if (idiag_Ymax(i)/= 0) then
            call max_mn_name(f(l1:l2,m,n,ichemspec(i)),idiag_Ymax(i))
          endif
          if (idiag_Ymin(i)/= 0) then
            call max_mn_name(-f(l1:l2,m,n,ichemspec(i)),idiag_Ymin(i),lneg=.true.)
          endif
        enddo
!
        if (idiag_num /= 0) call sum_mn_name(f(l1:l2,m,n,iviscosity), idiag_num)
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
        do i = 1,nchemspec
          if (idiag_Ymz(i)/= 0) &
            call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(i)),idiag_Ymz(i))
        enddo
      endif
      call timing('dchemistry_dt','finished',mnloop=.true.)
!
    endsubroutine dchemistry_dt
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
      character(len=6) :: diagn_dYm
      character(len=6) :: diagn_hm
      character(len=6) :: diagn_cpm
      character(len=7) :: diagn_diffm
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
        idiag_dYm = 0
        idiag_Ymax = 0
        idiag_Ymin = 0
        idiag_dYmax = 0
        idiag_hm = 0
        idiag_cpm = 0
        idiag_diffm = 0
!
        idiag_e_intm = 0
        idiag_Ymz = 0
!
        idiag_lambdam = 0
        idiag_num = 0
!
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname = 1,nname
        do ii=1,nchemspec
          write (number,'(I2)') ii
          diagn_Ym = 'Y'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_Ym),idiag_Ym(ii))
          diagn_Ymax = 'Y'//trim(adjustl(number))//'max'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_Ymax),idiag_Ymax(ii))
          diagn_Ymin = 'Y'//trim(adjustl(number))//'min'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_Ymin),idiag_Ymin(ii))
          diagn_dYm = 'dY'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_dYm),idiag_dYm(ii))
          diagn_dYmax = 'dY'//trim(adjustl(number))//'max'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_dYmax),idiag_dYmax(ii))
          diagn_hm = 'h'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_hm),idiag_hm(ii))
          diagn_cpm = 'cp'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_cpm),idiag_cpm(ii))
          diagn_diffm = 'diff'//trim(adjustl(number))//'m'
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_diffm),idiag_diffm(ii))
        enddo
        call parse_name(iname,cname(iname),cform(iname),'dtchem',idiag_dtchem)
!
!   Sample for hard-coded heat capacity diagnostics 
!
!        call parse_name(iname,cname(iname),cform(iname),'cp',idiag_cp)
        call parse_name(iname,cname(iname),cform(iname),'e_intm',idiag_e_intm)
        call parse_name(iname,cname(iname),cform(iname),'lambdam',idiag_lambdam)
        call parse_name(iname,cname(iname),cform(iname),'num',idiag_num)
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
          if (get_species_nr(sname,'chemspec',nchemspec,'rprint_chemistry')>0) &
            cformv(iname)='DEFINED'
        endif
        if (sname(1:9)=='viscosity') cformv(iname)='DEFINED'
        if (sname(1:2)=='cp') cformv(iname)='DEFINED'
        if (sname(1:8)=='reac_CO2') cformv(iname)='DEFINED'
        if (sname(1:7)=='reac_CO') cformv(iname)='DEFINED'
        if (sname(1:7)=='reac_O2') cformv(iname)='DEFINED'
        if (sname(1:4)=='reac') cformv(iname)='DEFINED'
      enddo
!
!  Write chemistry index in short notation
!
      if (lwr) then
        call farray_index_append('nchemspec',nchemspec)
        call farray_index_append('ichemspec',ichemspec(1),1,nchemspec)
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
      if (sname .ne. 'viscosity') then
        if (sname(9:)==' ') then    ! 9=len('chemspec')+1
          ispec=1
        else
          read(sname(9:),'(i3)') ispec
        endif
      call assign_slices_scal(slices,f,ichemspec(ispec))
      endif
!
    endsubroutine get_slices_chemistry
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
      if ((ChemInpLine(StartInd:StopInd) /= "M" ) &
          .and. (ChemInpLine(StartInd:StartInd+1) /= "hv" )) then
        StartSpecie = verify(ChemInpLine(StartInd:StopInd), &
            "1234567890")+StartInd-1
!
!  Call to a routine that checks for arbitrary stoiciometric coefficents
!  removes them and shifts the begin of the species index in cheminpline
!
        if (StopInd-StartSpecie >= 3 .and. StartSpecie > 1) &
            call find_remove_real_stoic(ChemInpLine(StartSpecie-1:StopInd),lreal,stoi,StartSpecie)
!
        call find_species_index(ChemInpLine(StartSpecie:StopInd), &
            ind_glob,ind_chem,found_specie)
        if (.not. found_specie) then
          print*,'ChemInpLine=',ChemInpLine
          print*,'StartSpecie,StopInd=',StartSpecie,StopInd
          print*,'ChemInpLine(StartSpecie:StopInd)=',ChemInpLine(StartSpecie:StopInd)
          print*,'ind_glob,ind_chem=',ind_glob,ind_chem
!          if (.not. lpencil_check_small) then
!          if (.not. lpencil_check) then
          call stop_it("build_stoich_matrix: Did not find species!")
!          endif
!          endif
        endif
!
!        if (found_specie) then
        if (StartSpecie == StartInd) then
          stoi = 1.0
        else
          if (.not. lreal) then
            read (unit=ChemInpLine(StartInd:StartInd),fmt='(I1)') stoi_int
          endif
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
    subroutine write_reactions()
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
            if (Sijp(spec,reac) /= 1) then
              write (Sijp_string,'(F3.1)') Sijp(spec,reac)
            endif
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
            write (file_id,*) a_k4(:,reac)
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
    subroutine chemkin_data
!
!  if the file with chemkin data exists
!  reading the Chemkin data
!
      character(len=fnlen) :: input_file
      integer :: stat, k,i
      character(len=fnlen) :: input_file2="./data/stoich.out"
      integer :: file_id=123
      logical :: chemin,cheminp
!
      inquire(file='chem.inp',exist=cheminp)
      inquire(file='chem.in',exist=chemin)
      if (chemin .and. cheminp) call fatal_error('chemistry',&
          'chem.in and chem.inp found, please decide for one')
      if (cheminp) input_file='chem.inp'
      if (chemin) input_file='chem.in'
!
      inquire (file='tran.dat',exist=tran_exist)
      if (.not. tran_exist) then
        inquire (file='tran.in',exist=tran_exist)
      endif
!
      if (tran_exist) then
        if (lroot) then
          print*,'tran.in/dat file with transport data is found.'
        endif
        call read_transport_data
      endif
!
      if (lroot .and. .not. tran_exist .and. .not. lew_exist) then
        if (chem_diff == 0.) &
            call inevitably_fatal_error('chemkin data', 'chem_diff = 0.')
        print*,'tran.dat file with transport data is not found.'
        print*,'lewis.dat file with Lewis numbers is not found.'
        print*,'Now diffusion coefficients is ',chem_diff
        print*,'Now species viscosity is ',nu_spec
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
        if (stat > 0) call stop_it("Couldn't allocate memory for stoichio")
        allocate(Sijm(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for Sijm")
        allocate(Sijp(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for Sijp")
        allocate(reaction_name(mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for reaction_name")
        allocate(back(mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for back")
        allocate(B_n(mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for B_n")
        B_n = 0.
        allocate(alpha_n(mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for alpha_n")
        alpha_n = 0.
        allocate(E_an(mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for E_an")
        E_an = 0.
!
        allocate(low_coeff(3,nreactions),STAT=stat)
        low_coeff = 0.
        if (stat > 0) call stop_it("Couldn't allocate memory for low_coeff")
        allocate(high_coeff(3,nreactions),STAT=stat)
        high_coeff = 0.
        if (stat > 0) call stop_it("Couldn't allocate memory for high_coeff")
        allocate(troe_coeff(3,nreactions),STAT=stat)
        troe_coeff = 0.
        if (stat > 0) call stop_it("Couldn't allocate memory for troe_coeff")
        allocate(a_k4(nchemspec,nreactions),STAT=stat)
        a_k4 = impossible
        if (stat > 0) call stop_it("Couldn't allocate memory for troe_coeff")
        allocate(Mplus_case (nreactions),STAT=stat)
        Mplus_case = .false.
        allocate(photochem_case (nreactions),STAT=stat)
        photochem_case = .false.
        if (stat > 0) call stop_it("Couldn't allocate memory for photochem_case")
        if (lmech_simple) then
          allocate(orders_p(5,nreactions),STAT=stat)
          allocate(orders_m(5,nreactions),STAT=stat)
          orders_p = 0.
          orders_m = 0.
          if (stat > 0) call stop_it("Couldn't allocate memory for orders_p/m")
        endif
      endif
!
!  Initialize data
!
      Sijp = 0.
      Sijm = 0.0
      back = .true.
!
!  read chemistry data
!
      call read_reactions(input_file)
      call write_reactions()
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
100   format(I1,26f6.1)
101   format('    ',26A6)
!
    endsubroutine chemkin_data
!***********************************************************************
    subroutine chemkin_data_simple
!
      character(len=fnlen) :: input_file
      integer :: stat, k,i
      character(len=fnlen) :: input_file2="./data/stoich.out"
      integer :: file_id=123
      logical :: chemin,cheminp
!
! CO2 CO N2 O2 H2O 
!
    if (lmech_simple) then
        ! O2
        tran_data(4,1:7) = (/ 1.0000000000000000        ,107.40000000000001        ,3.4580000000000002        ,&
                              0.0000000000000000E+000   ,1.6000000000000001        ,3.7999999999999998        ,&
                              0.0000000000000000E+000 /)
        ! N2
        tran_data(3,1:7) = (/ 1.0000000000000000        ,97.530000000000001        ,3.6210000000000000        ,&
                               0.0000000000000000E+000   ,1.7600000000000000        ,4.0000000000000000       ,&
                               0.0000000000000000E+000 /)
        ! CO
        tran_data(2,1:7) = (/ 1.0000000000000000        ,98.099999999999994        ,3.6499999999999999        ,&
                               0.0000000000000000E+000   ,1.9500000000000000        ,1.8000000000000000       ,&
                               0.0000000000000000E+000 /)
        ! CO2
        tran_data(1,1:7) = (/ 1.0000000000000000        ,244.00000000000000        ,3.7629999999999999        ,&
                               0.0000000000000000E+000   ,2.6499999999999999        ,2.1000000000000001       ,&
                               0.0000000000000000E+000 /)
        if (nchemspec==5) then
          ! H2O
          tran_data(5,1:7) = (/ 2.0000000000000000        ,572.39999999999998        ,2.6050000000000000        ,&
                                1.8440000000000001        ,0.0000000000000000E+000   ,4.0000000000000000        ,&
                                0.0000000000000000E+000 /)
        endif
!
!  Find number of ractions
!
      mreactions = 1
      if (lroot) print*,'Number of reactions=',mreactions
      if (lroot) print*,'Number of species=',nchemspec
      nreactions = mreactions
!
    else
!
        tran_data(1,1:7) = (/ 1.0000000000000000        ,38.000000000000000        ,2.9199999999999999        ,&
                              0.0000000000000000E+000  ,0.79000000000000004        ,280.00000000000000        ,&
                              0.0000000000000000E+000 /)
        tran_data(2,1:7) = (/ 1.0000000000000000        ,107.40000000000001        ,3.4580000000000002        ,&
                              0.0000000000000000E+000   ,1.6000000000000001        ,3.7999999999999998        ,&
                              0.0000000000000000E+000 /)
        tran_data(3,1:7) = (/ 2.0000000000000000        ,572.39999999999998        ,2.6050000000000000        ,&
                              1.8440000000000001        ,0.0000000000000000E+000   ,4.0000000000000000        ,&
                              0.0000000000000000E+000 /)
        tran_data(4,1:7) = (/ 0.0000000000000000E+000   ,145.00000000000000        ,2.0499999999999998        ,&
                              0.0000000000000000E+000   ,0.0000000000000000E+000   ,0.0000000000000000E+000   ,&
                              0.0000000000000000E+000 /)
        tran_data(5,1:7) = (/ 0.0000000000000000E+000   ,80.000000000000000        ,2.7500000000000000        ,&
                              0.0000000000000000E+000   ,0.0000000000000000E+000   ,0.0000000000000000E+000   ,&
                              0.0000000000000000E+000 /)
        tran_data(6,1:7) = (/ 1.0000000000000000        ,80.000000000000000        ,2.7500000000000000        ,&
                              0.0000000000000000E+000   ,0.0000000000000000E+000   ,0.0000000000000000E+000   ,&
                              0.0000000000000000E+000 /)
        tran_data(7,1:7) = (/ 2.0000000000000000        ,107.40000000000001        ,3.4580000000000002        ,&
                              0.0000000000000000E+000   ,0.0000000000000000E+000   ,1.0000000000000000        ,&
                              0.0000000000000000E+000 /)
        tran_data(8,1:7) = (/ 2.0000000000000000        ,107.40000000000001        ,3.4580000000000002        ,&
                              0.0000000000000000E+000   ,0.0000000000000000E+000   ,3.7999999999999998        ,&
                              0.0000000000000000E+000 /)
        tran_data(9,1:7) = (/ 0.0000000000000000E+000   ,136.50000000000000        ,3.3300000000000001        ,&
                              0.0000000000000000E+000   ,0.0000000000000000E+000   ,0.0000000000000000E+000   ,&
                              0.0000000000000000E+000 /)
        tran_data(10,1:7) = (/ 1.0000000000000000        ,97.530000000000001        ,3.6210000000000000        ,&
                               0.0000000000000000E+000   ,1.7600000000000000        ,4.0000000000000000        ,&
                               0.0000000000000000E+000 /)
        tran_data(11,1:7) = (/ 0.0000000000000000E+000   ,10.199999999999999        ,2.5760000000000001        ,&
                               0.0000000000000000E+000   ,0.0000000000000000E+000   ,0.0000000000000000E+000   ,&
                               0.0000000000000000E+000 /)
        tran_data(12,1:7) = (/ 1.0000000000000000        ,98.099999999999994        ,3.6499999999999999        ,&
                               0.0000000000000000E+000   ,1.9500000000000000        ,1.8000000000000000        ,&
                               0.0000000000000000E+000 /)
        tran_data(13,1:7) = (/ 1.0000000000000000        ,244.00000000000000        ,3.7629999999999999        ,&
                               0.0000000000000000E+000   ,2.6499999999999999        ,2.1000000000000001        ,&
                               0.0000000000000000E+000 /)
!
!  Find number of ractions
!
      mreactions = 25
      if (lroot) print*,'Number of reactions=',mreactions
      if (lroot) print*,'Number of species=',nchemspec
      nreactions = mreactions
!
    endif
!
!  Allocate reaction arrays
!
      if (.not. lreloading) then
        allocate(stoichio(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for stoichio")
        allocate(Sijm(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for Sijm")
        allocate(Sijp(nchemspec,mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for Sijp")
        allocate(reaction_name(mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for reaction_name")
        allocate(back(mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for back")
        allocate(B_n(mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for B_n")
        B_n = 0.
        allocate(alpha_n(mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for alpha_n")
        alpha_n = 0.
        allocate(E_an(mreactions),STAT=stat)
        if (stat > 0) call stop_it("Couldn't allocate memory for E_an")
        E_an = 0.
!
        allocate(low_coeff(3,nreactions),STAT=stat)
        low_coeff = 0.
        if (stat > 0) call stop_it("Couldn't allocate memory for low_coeff")
        allocate(high_coeff(3,nreactions),STAT=stat)
        high_coeff = 0.
        if (stat > 0) call stop_it("Couldn't allocate memory for high_coeff")
        allocate(troe_coeff(3,nreactions),STAT=stat)
        troe_coeff = 0.
        if (stat > 0) call stop_it("Couldn't allocate memory for troe_coeff")
        allocate(a_k4(nchemspec,nreactions),STAT=stat)
        a_k4 = impossible
        if (stat > 0) call stop_it("Couldn't allocate memory for troe_coeff")
        allocate(Mplus_case (nreactions),STAT=stat)
        Mplus_case = .false.
        allocate(photochem_case (nreactions),STAT=stat)
        photochem_case = .false.
        if (stat > 0) call stop_it("Couldn't allocate memory for photochem_case")
        if (lmech_simple) then
          allocate(orders_p(5,nreactions),STAT=stat)
          allocate(orders_m(5,nreactions),STAT=stat)
          orders_p = 0.
          orders_m = 0.
          if (stat > 0) call stop_it("Couldn't allocate memory for orders_p/m")
        endif
      endif
!
!  Initialize data
!
      Sijp = 0.
      Sijm = 0.0
      back = .true.
!
    if (lmech_simple) then
!
! CO2 CO N2 O2 H2O 
!
      Sijp(1:4,1) = (/ &
        0.0000000000000000      ,&
        1.0000000000000000      ,&
        0.0000000000000000E+000 ,&
        0.5000000000000000E+000 /)
!
      Sijm(1:4,1) = (/ &
        1.0000000000000000      ,&
        0.0000000000000000E+000 ,&
        0.0000000000000000E+000 ,&
        0.0000000000000000E+000 /)
!
      if (nchemspec==5) then
        Sijp(nchemspec,1) = 0.0
        Sijm(nchemspec,1) = 0.0
      endif
!
      B_n = 33.6174731212 + scale_homo
     ! B_n = 32.680877  !changed A to lower S_l
      alpha_n = 0.0
      !! E_an in CAL !!
      E_an = 40000
!
      back(:) = lback
      Mplus_case(:) = .false.
      photochem_case(:) = .false.
!
      low_coeff(:,:) = 0.0
      troe_coeff(:,:) = 0.0
      high_coeff(:,:) = 0.0
!
      a_k4(:,:) = 3.9084999999999999E+037
!
      orders_p(:,1) = (/ 0.0, 1.0, 0.0, 0.25, 0.5 /)
      orders_m(:,1) = (/ 1.0, 0.0, 0.0, 0.0, 0.0 /)
!
    else
!
!  read chemistry data
!
 Sijp(:,           1 ) = (/ &
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,           2 ) = (/ &
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,           3 ) = (/ &
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,           4 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,           5 ) = (/ &
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,           6 ) = (/ &
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,           7 ) = (/ &
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,           8 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,           9 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          10 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          11 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          12 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          13 ) = (/ &
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          14 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          15 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          16 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          17 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          18 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          19 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          20 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          21 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          22 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          23 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          24 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijp(:,          25 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,           1 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,           2 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,           3 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,           4 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,           5 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,           6 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,           7 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,           8 ) = (/ &
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,           9 ) = (/ &
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          10 ) = (/ &
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          11 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          12 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          13 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          14 ) = (/ &
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          15 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          16 ) = (/ &
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          17 ) = (/ &
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          18 ) = (/ &
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          19 ) = (/ &
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          20 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          21 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          22 ) = (/ &
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          23 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          24 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 Sijm(:,          25 ) = (/ &
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 B_n = (/&
   35.804878570702172      ,&
   10.835651633566574      ,&
   19.190788965648441      ,&
   14.904072510778882      ,&
   45.270160528558371      ,&
   43.211262470732954      ,&
   43.211262470732954      ,&
   36.357664531526993      ,&
   30.568064393133859      ,&
   30.568064393133859      ,&
   42.997068478406760      ,&
   51.991873112601347      ,&
   28.019679105720332      ,&
   30.440423811291044      ,&
   31.890738863714645      ,&
   31.112261205264240      ,&
   30.994862711046935      ,&
   33.671275827205960      ,&
   25.590800287401994      ,&
   33.318335397877441      ,&
   30.813232956425157      ,&
   31.506380136985104      ,&
   16.072051712456911      ,&
   27.631021115928547      ,&
   33.994049219469012      /)
 alpha_n = (/&
 -0.40600000000000003      ,&
   2.6699999999999999      ,&
   1.5100000000000000      ,&
   2.0200000000000000      ,&
  -1.3999999999999999      ,&
  -1.1000000000000001      ,&
  -1.1000000000000001      ,&
 -0.50000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
  -1.0000000000000000      ,&
  -2.0000000000000000      ,&
  0.59999999999999998      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   2.0000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 E_an = (/&
   16599.000000000000      ,&
   6290.0000000000000      ,&
   3430.0000000000000      ,&
   13400.000000000000      ,&
   104380.00000000000      ,&
   104380.00000000000      ,&
   104380.00000000000      ,&
   0.0000000000000000E+000 ,&
  -1788.0000000000000      ,&
  -1788.0000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   823.00000000000000      ,&
   295.00000000000000      ,&
   0.0000000000000000E+000 ,&
  -497.00000000000000      ,&
   11982.000000000000      ,&
  -1629.3000000000000      ,&
   48430.000000000000      ,&
   3970.0000000000000      ,&
   7950.0000000000000      ,&
   3970.0000000000000      ,&
   0.0000000000000000E+000 ,&
   9557.0000000000000      /)
!
   back(:) = .true.
   Mplus_case(:) = .false.
   Mplus_case(13) = .true.
   Mplus_case(20) = .true.
   photochem_case(:) = .false.

 low_coeff(           1 ,:) = (/&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   47.902673188740813      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   39.327933417011792      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 low_coeff(           2 ,:) = (/&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
  -1.7200000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 low_coeff(           3 ,:) = (/&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   524.79999999999995      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   45500.000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 troe_coeff(           1 ,:) = (/&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
  0.80000000000000004      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
  0.50000000000000000      ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 troe_coeff(           2 ,:) = (/&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   9.9999999999999996E-039 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   9.9999999999999996E-039 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 troe_coeff(           3 ,:) = (/&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000E+022 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   1.0000000000000000E+022 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 high_coeff(           1 ,:) = (/&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 high_coeff(           2 ,:) = (/&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
 high_coeff(           3 ,:) = (/&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 ,&
   0.0000000000000000E+000 /)
!
 a_k4(           1 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   2.5000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   2.5000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   2.5000000000000000      ,&
   2.5000000000000000      ,&
   2.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   2.5000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(           2 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
  0.78000000000000003      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(           3 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   12.000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   12.000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   12.000000000000000      ,&
   12.000000000000000      ,&
   11.000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   12.000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(           4 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(           5 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(           6 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(           7 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(           8 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(           9 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   0.0000000000000000E+000 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   0.0000000000000000E+000 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
  0.75000000000000000      ,&
  0.38000000000000000      ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
  0.64000000000000001      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(          10 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(          11 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   0.0000000000000000E+000 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
  0.75000000000000000      ,&
  0.38000000000000000      ,&
   1.0000000000000000      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
  0.64000000000000001      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(          12 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.8999999999999999      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.8999999999999999      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.8999999999999999      ,&
   1.8999999999999999      ,&
   1.8999999999999999      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   1.8999999999999999      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
 a_k4(          13 ,:) = (/&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.7999999999999998      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.7999999999999998      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.7999999999999998      ,&
   3.7999999999999998      ,&
   3.7999999999999998      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.7999999999999998      ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 ,&
   3.9084999999999999E+037 /)
!
    endif
!
      call write_reactions()
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
100   format(I1,26f6.1)
101   format('    ',26A6)
!
    endsubroutine chemkin_data_simple
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
! Find reactant side stoichiometric coefficients
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
                                  StopInd_add),fmt='(E15.8)') low_coeff(1,k)
                              if (low_coeff(1,k)/=0.) low_coeff(1,k)=log(low_coeff(1,k))
                            elseif (VarNumber_add==2) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') low_coeff(2,k)
                            elseif (VarNumber_add==3) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') low_coeff(3,k)
                            else
                              call stop_it("No such VarNumber!")
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
                              call stop_it("No such VarNumber!")
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
                                  StopInd_add),fmt='(E15.8)') high_coeff(1,k)
                              if (high_coeff(1,k)/=0.) high_coeff(1,k)=log(high_coeff(1,k))
                            elseif (VarNumber_add==2) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') high_coeff(2,k)
                            elseif (VarNumber_add==3) then
                              read (unit=ChemInpLine_add(StartInd_add:&
                                  StopInd_add),fmt='(E15.8)') high_coeff(3,k)
                            else
                              call stop_it("No such VarNumber!")
                            endif
                          endif
                          VarNumber_add=VarNumber_add+1
                          !   StartInd_add=StopInd_add
                          StartInd_add=verify(ChemInpLine_add(StopInd_add+1:),&
                              ' ')+StopInd_add
                          StopInd_add=StartInd_add
                        enddo
                        i=80
                      else
                        !                a_k4=0.
                        StartInd_add=i; StopInd_add=0; StopInd_add_=0
                        do while (ChemInpLine_add(i:i+1)/='  ')
                          find_specie=.true.
                          do while (StartInd_add/=StopInd_add_)
                            StopInd_add=index(ChemInpLine_add(StartInd_add:),&
                                '/')+StartInd_add-2
                            StopInd_add_=index(ChemInpLine_add(StartInd_add:),&
                                ' ')+StartInd_add-1
                            if (find_specie) then
                              call find_species_index(ChemInpLine_add&
                                  (StartInd_add:StopInd_add),ind_glob,&
                                  ind_chem,found_specie)
                            else
                              if (found_specie) then
                                read (unit=ChemInpLine_add(StartInd_add:&
                                    StopInd_add),fmt='(E15.8)') a_k4(ind_chem,k)
                              else
                                print*,'ChemInpLine=',ChemInpLine_add
                                print*,'Specie=',&
                                    ChemInpLine_add(StartInd_add:StopInd_add)
                                print*,'StartInd_add,StopInd_add=',&
                                    StartInd_add,StopInd_add
                                call stop_it("read_reactions: Did not find specie!")
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
                PlusInd=index(ChemInpLine(StartInd:LastLeftCharacter),'+')&
                    +StartInd-1
                do while (PlusInd<LastLeftCharacter .AND. PlusInd>0)
                  StopInd=PlusInd-1
                  call build_stoich_matrix(StartInd,StopInd,k,&
                      ChemInpLine,.false.)
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
                PlusInd=index(ChemInpLine(StartInd:LastLeftCharacter),'+')&
                    +StartInd-1
                do while (PlusInd<LastLeftCharacter .AND. PlusInd>StartInd)
                  StopInd=PlusInd-1
                  call build_stoich_matrix(StartInd,StopInd,k,&
                      ChemInpLine,.true.)
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
                      read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)')  &
                          B_n(k)
                      if (B_n(k)/=0.) B_n(k)=log(B_n(k))
                    elseif (VarNumber==2) then
                      read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)')  &
                          alpha_n(k)
                    elseif (VarNumber==3) then
                      read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)')  &
                          E_an(k)
                    else
                      call stop_it("No such VarNumber!")
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
      real, dimension(nx) :: dSR=0., dHRT=0., Kp, Kc
      real, dimension(nx) :: prod1, prod2
      real, dimension(nx) :: kf=0., kr=0.
      real, dimension(nx) :: rho_cgs, p_atm
      real, dimension(nx) :: mix_conc
      integer :: k, reac, i
      real :: sum_tmp=0., ddd
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
!
        Rcal = Rgas_unit_sys/4.14*1e-7
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
        if (lmech_simple) then
          prod1 = 1.
          prod2 = 1.
          do k = 1,nchemspec
            if (abs(orders_p(k,reac)) > 0.0) then
              prod1 = prod1*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))**orders_p(k,reac)
            endif
            if (abs(orders_m(k,reac)) > 0.0) then
              prod2 = prod2*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))**orders_m(k,reac)
            endif
          enddo
          if (nchemspec==4) prod1=prod1*(Y_H2O*rho_cgs(:)/m_H2O)**orders_p(5,reac)
        else
          prod1 = 1.
          prod2 = 1.
          do k = 1,nchemspec
            if (abs(Sijp(k,reac)) == 1) then
              prod1 = prod1*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))
            elseif (abs(Sijp(k,reac)) == 2) then
              prod1 = prod1*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))*(f(l1:l2,m,n,ichemspec(k)) &
                  *rho_cgs(:)/species_constants(k,imass))
            elseif (abs(Sijp(k,reac)) > 0) then
              prod1 = prod1*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))**Sijp(k,reac)
            endif
          enddo
          do k = 1,nchemspec
            if (abs(Sijm(k,reac)) == 1.0) then
              prod2 = prod2*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))
            elseif (abs(Sijm(k,reac)) == 2.0) then
              prod2 = prod2*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))*(f(l1:l2,m,n,ichemspec(k)) &
                  *rho_cgs(:)/species_constants(k,imass))
            elseif (abs(Sijm(k,reac)) > 0.0) then
              prod2 = prod2*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))**Sijm(k,reac)
            endif
          enddo
        endif
!
!  Find forward rate constant for reaction 'reac'
!
        kf = B_n(reac)+alpha_n(reac)*p%lnTT-E_an(reac)*Rcal1*TT1_loc
!
!  Find backward rate constant for reaction 'reac'
!
        if (lmech_simple) then
          !! 20.0301186564 = ln(5*10e8) !!
          kr = 20.0301186564-E_an(reac)*Rcal1*TT1_loc
          !kr = 19.0833687-E_an(reac)*Rcal1*TT1_loc !changed A to obtain flamespeed as in GRI
          sum_sp = 1.
        else

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
          if (minval(a_k4(:,reac)) < impossible) then
            sum_sp = 0.
            do k = 1,nchemspec
              sum_sp = sum_sp+a_k4(k,reac)*f(l1:l2,m,n,ichemspec(k))  &
                  *rho_cgs(:)/species_constants(k,imass)
            enddo
            mix_conc = sum_sp
          else
            sum_sp = 1.
!            mix_conc = rho_cgs(:)*p%mu1(:)/unit_mass
            mix_conc = rho_cgs(:)*(p%RRmix/Rgas)/unit_mass
          endif
!
!  The Lindeman approach to the fall of reactions
!
          if (maxval(abs(low_coeff(:,reac))) > 0.) then
            B_n_0 = low_coeff(1,reac)
            alpha_n_0 = low_coeff(2,reac)
            E_an_0 = low_coeff(3,reac)
            kf_0(:) = B_n_0+alpha_n_0*p%lnTT(:)-E_an_0*Rcal1*TT1_loc(:)
            Pr = exp(kf_0-kf)*mix_conc
            kf = kf+log(Pr/(1.+Pr))
          elseif (maxval(abs(high_coeff(:,reac))) > 0.) then
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
          if (maxval(abs(troe_coeff(:,reac))) > 0.) then
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
            kr = kf-Kc
!
        endif

          if (Mplus_case (reac)) then
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
      endif
!
      if (lwrite_first .and. lroot) &
          print*,'get_reaction_rate: writing react.out file'
      lwrite_first = .false.
!
    endsubroutine get_reaction_rate
!***********************************************************************
    subroutine calc_reaction_term(f,p)
!
!  Calculation of the reaction term
!
      use Diagnostics, only: sum_mn_name, max_mn_name
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,mreactions) :: vreactions, vreactions_p, vreactions_m
      real, dimension(nx,nchemspec) :: xdot
      real, dimension(nx) :: rho1
      real, dimension(nx,nchemspec) :: molm
      type (pencil_case) :: p
      integer :: k,j,ii
      integer :: i1=1, i2=2, i3=3, i4=4, i5=5, i6=6, i7=7, i8=8, i9=9, i10=10
      integer :: i11=11, i12=12, i13=13, i14=14, i15=15, i16=16, i17=17, i18=18, i19=19
!
      p%DYDt_reac = 0.
      rho1 = 1./p%rho
      do k = 1,nchemspec
        molm(:,k) = rho1*species_constants(k,imass)
      enddo
!
!  Chemkin data case
!
        call get_reaction_rate(f,vreactions_p,vreactions_m,p)
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
      if (t < intro_time) p%DYDt_reac = p%DYDt_reac*t/intro_time
!
!  For diagnostics
!
      if (lchemistry_diag) then
        do k = 1,nchemspec
          do j = 1,nreactions
            net_react_p(k,j) = net_react_p(k,j)+stoichio(k,j) &
                *sum(vreactions_p(:,j))
            net_react_m(k,j) = net_react_m(k,j)+stoichio(k,j) &
                *sum(vreactions_m(:,j))
          enddo
        enddo
      endif
!
!  Calculate diagnostic quantities
!
      if (ldiagnos) then
        do ii=1,nchemspec
          if (idiag_dYm(ii) /= 0) then
            call sum_mn_name(p%DYDt_reac(:,ii),idiag_dYm(ii))
          endif
          if (idiag_dYmax(ii) /= 0) then
            call max_mn_name(abs(p%DYDt_reac(:,ii)),idiag_dYmax(ii))
          endif
          if (idiag_hm(ii) /= 0) then
            call sum_mn_name(p%H0_RT(:,ii)*Rgas* &
                p%TT(:)/species_constants(ii,imass),idiag_hm(ii))
          endif
        enddo
      endif
!
    endsubroutine calc_reaction_term
!***********************************************************************
    subroutine  write_net_reaction
!
!  write net reactions to file

      open (1,file=trim(datadir)//'/net_reactions.dat',position='append')
      write (1,*) t
      write (1,'(8e10.2)') net_react_p, net_react_m
      close (1)
!
    endsubroutine  write_net_reaction
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
      real, dimension(nx,3) :: gchemspec, sum_diff=0., dk_D
      real, dimension(nx,3,nchemspec) :: gDiff_full_add
      real, dimension(nx) :: del2chemspec, diff_op1, diff_op2
      real, dimension(nx) :: sum_gdiff=0., gY_sumdiff
      integer :: k,i
!
      intent(in) :: f
!
      p%DYDt_diff = 0.
!
!  Loop over all chemical species.
!
      do k = 1,nchemspec
!
!  Detailed chemistry and transport using CHEMKIN formalism.
!
        if (ldiffusion) then
!
!  Calculate diffusion coefficient gradients gDiff_full_add for the case:
!    2) Constant Lewis number (= 1 by default or read from start.in) and given Pr
!
            do i = 1,3
              gDiff_full_add(:,i,k) = p%gradnu(:,i)/Pr_number*Lewis_coef1(k)
            enddo
!
!  Calculate the terms needed by the diffusion fluxes in a case:
!    1) Fickian diffusion law (gradient of species MASS fractions)
!
            call del2(f,ichemspec(k),del2chemspec)
            call grad(f,ichemspec(k),gchemspec)
            call dot_mn(p%glnrho,gchemspec,diff_op1)
            call dot_mn(gDiff_full_add(:,:,k),gchemspec,diff_op2)
!
            p%DYDt_diff(:,k) = p%Diff_penc_add(:,k) &
                             *(del2chemspec+diff_op1) + diff_op2
            do i = 1,3
              dk_D(:,i) = p%Diff_penc_add(:,k)*gchemspec(:,i)
            enddo
!
            if (ldiff_corr) then
              do i = 1,3
                sum_diff(:,i) = sum_diff(:,i)+dk_D(:,i)
              enddo
              sum_gdiff(:) = sum_gdiff(:)+ p%DYDt_diff(:,k)
            endif
        endif
      enddo
!
!  Adding correction diffusion velocity to ensure mass balance
!
      if (ldiffusion .and. ldiff_corr) then
        do k = 1,nchemspec
          call grad(f,ichemspec(k),gchemspec)
          call dot_mn(gchemspec(:,:),sum_diff,gY_sumdiff)
          p%DYDt_diff(:,k) = p%DYDt_diff(:,k) &
                - (gY_sumdiff+f(l1:l2,m,n,ichemspec(k))*sum_gdiff(:))
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
!  Add heat conduction to RHS of temperature equation
!
      if (ltemperature_nolog) then
        call dot(p%gTT,p%glambda,g2TTlambda)
        call del2(f,iTT,del2TT)
        tmp1 = (p%lambda(:)*del2TT+g2TTlambda)*p%cv1*p%rho1(:)
        df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + tmp1
      else
        call dot(p%glnTT,p%glambda,g2TTlambda)
        call dot(p%glnTT,p%glnTT,g2TT)
        tmp1 = (p%lambda(:)*(p%del2lnTT+g2TT)+g2TTlambda)*p%cv1*p%rho1(:)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + tmp1
      endif
!
!      call keep_compiler_quiet(f)
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
!
      call keep_compiler_quiet(cs2_full)
      call fatal_error('get_cs2_full',&
        'This function is not working with chemistry_simple since all full arrays are removed')
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
      intent(in) :: f
      integer :: j2, j3, k
      real, dimension(nchemspec) ::  cp_k, cv_k
      real, dimension (:,:), allocatable :: TT_full, cp_full, cv_full
!
      if (Cp_const < impossible) then
!
      do k = 1,nchemspec
        cp_k(k) = Cp_const/species_constants(k,imass)
        cv_k(k) = (Cp_const-Rgas)/species_constants(k,imass)
      enddo
!
      if (dir == 1) then
         allocate (TT_full(my,mz))  
         allocate (cp_full(my,mz))  
         allocate (cv_full(my,mz))  
         cp_full = 0.
         cv_full = 0.
         TT_full = 0.
         if (ltemperature_nolog) then
           TT_full(m1:m2,n1:n2) = f(index,m1:m2,n1:n2,iTT)
         else
           TT_full(m1:m2,n1:n2) = exp(f(index,m1:m2,n1:n2,ilnTT))
         endif
         do k=1,nchemspec
            if (species_constants(k,imass)>0.) then
               do j2=m1,m2
                  do j3=n1,n2
                     cp_full(j2,j3) = &
                          cp_full(j2,j3)+cp_k(k)*f(index,j2,j3,ichemspec(k))
                     cv_full(j2,j3) = &
                          cv_full(j2,j3)+cv_k(k)*f(index,j2,j3,ichemspec(k))
                  enddo
               enddo
            endif
         enddo
         slice = cp_full(m1:m2,n1:n2)/cv_full(m1:m2,n1:n2) &
               * f(index,m1:m2,n1:n2,iRR)*TT_full(m1:m2,n1:n2)
      elseif (dir == 2) then
         allocate (TT_full(mx,mz))  
         allocate (cp_full(mx,mz))  
         allocate (cv_full(mx,mz))  
         cp_full = 0.
         cv_full = 0.
         TT_full = 0.
         if (ltemperature_nolog) then
           TT_full(l1:l2,n1:n2) = f(l1:l2,index,n1:n2,iTT)
         else
           TT_full(l1:l2,n1:n2) = exp(f(l1:l2,index,n1:n2,ilnTT))
         endif
         do k=1,nchemspec
            if (species_constants(k,imass)>0.) then
               do j2=l1,l2
                  do j3=n1,n2
                     cp_full(j2,j3) = &
                          cp_full(j2,j3)+cp_k(k)*f(j2,index,j3,ichemspec(k))
                     cv_full(j2,j3) = &
                          cv_full(j2,j3)+cv_k(k)*f(j2,index,j3,ichemspec(k))
                  enddo
               enddo
            endif
         enddo
         slice = cp_full(l1:l2,n1:n2)/cv_full(l1:l2,n1:n2) &
               * f(l1:l2,index,n1:n2,iRR)*TT_full(l1:l2,n1:n2)
      elseif (dir == 3) then
         allocate (TT_full(mx,my))  
         allocate (cp_full(mx,my))  
         allocate (cv_full(mx,my))  
         cp_full = 0.
         cv_full = 0.
         TT_full = 0.
         if (ltemperature_nolog) then
           TT_full(l1:l2,m1:m2) = f(l1:l2,m1:m2,index,iTT)
         else
           TT_full(l1:l2,m1:m2) = exp(f(l1:l2,m1:m2,index,ilnTT))
         endif
         do k=1,nchemspec
            if (species_constants(k,imass)>0.) then
               do j2=l1,l2
                  do j3=m1,m2
                     cp_full(j2,j3) = &
                          cp_full(j2,j3)+cp_k(k)*f(j2,j3,index,ichemspec(k))
                     cv_full(j2,j3) = &
                          cv_full(j2,j3)+cv_k(k)*f(j2,j3,index,ichemspec(k))
                  enddo
               enddo
            endif
         enddo
         slice = cp_full(l1:l2,m1:m2)/cv_full(l1:l2,m1:m2) &
               * f(l1:l2,m1:m2,index,iRR)*TT_full(l1:l2,m1:m2)
      else
         call fatal_error('get_cs2_slice','No such dir!')
      endif
!
      deallocate (TT_full)  
      deallocate (cp_full)  
      deallocate (cv_full)  
!
      else
!
      if (dir == 1) then
         allocate (TT_full(my,mz))  
         allocate (cp_full(my,mz))  
         allocate (cv_full(my,mz))  
         if (ltemperature_nolog) then
           TT_full(m1:m2,n1:n2) = f(index,m1:m2,n1:n2,iTT)
         else
           TT_full(m1:m2,n1:n2) = exp(f(index,m1:m2,n1:n2,ilnTT))
         endif
         cp_full(m1:m2,n1:n2) = f(index,m1:m2,n1:n2,icp)
         cv_full = cp_full - f(index,:,:,iRR)
         slice = cp_full(m1:m2,n1:n2)/cv_full(m1:m2,n1:n2) &
               * f(index,m1:m2,n1:n2,iRR)*TT_full(m1:m2,n1:n2)
      elseif (dir == 2) then
         allocate (TT_full(mx,mz))  
         allocate (cp_full(mx,mz))  
         allocate (cv_full(mx,mz))  
         if (ltemperature_nolog) then
           TT_full(l1:l2,n1:n2) = f(l1:l2,index,n1:n2,iTT)
         else
           TT_full(l1:l2,n1:n2) = exp(f(l1:l2,index,n1:n2,ilnTT))
         endif
         cp_full(l1:l2,n1:n2) = f(l1:l2,index,n1:n2,icp)
         cv_full = cp_full - f(:,index,:,iRR)
         slice = cp_full(l1:l2,n1:n2)/cv_full(l1:l2,n1:n2) &
               * f(l1:l2,index,n1:n2,iRR)*TT_full(l1:l2,n1:n2)
      elseif (dir == 3) then
         allocate (TT_full(mx,my))  
         allocate (cp_full(mx,my))  
         allocate (cv_full(mx,my))  
         if (ltemperature_nolog) then
           TT_full(l1:l2,m1:m2) = f(l1:l2,m1:m2,index,iTT)
         else
           TT_full(l1:l2,m1:m2) = exp(f(l1:l2,m1:m2,index,ilnTT))
         endif
         cp_full(l1:l2,m1:m2) = f(l1:l2,m1:m2,index,icp)
         cv_full = cp_full - f(:,:,index,iRR)
         slice = cp_full(l1:l2,m1:m2)/cv_full(l1:l2,m1:m2) &
               * f(l1:l2,m1:m2,index,iRR)*TT_full(l1:l2,m1:m2)
      else
         call fatal_error('get_cs2_slice','No such dir!')
      endif
!
      deallocate (TT_full)  
      deallocate (cp_full)  
      deallocate (cv_full)  
!
    endif  
!
    endsubroutine get_cs2_slice
!***********************************************************************
    subroutine get_gamma_full(gamma_full)
!
      real, dimension(mx,my,mz) :: gamma_full
      intent(out) :: gamma_full
!
      call keep_compiler_quiet(gamma_full)
      call fatal_error('get_gamma_full',&
        'This function is not working with chemistry_simple since all full arrays are removed')
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
      intent(in) :: f
      integer :: j2, j3, k
      real, dimension(nchemspec) ::  cp_k, cv_k
      real, dimension (:,:), allocatable :: cp_full, cv_full
!
    if (Cp_const < impossible) then
!
      do k = 1,nchemspec
        cp_k(k) = Cp_const/species_constants(k,imass)
        cv_k(k) = (Cp_const-Rgas)/species_constants(k,imass)
      enddo
!
      if (dir == 1) then
         allocate (cp_full(my,mz))  
         allocate (cv_full(my,mz))  
         cp_full = 0.
         cp_full(m1:m2,n1:n2) = f(index,m1:m2,n1:n2,icp)
         cv_full = 0.
         do k=1,nchemspec
            if (species_constants(k,imass)>0.) then
               do j2=m1,m2
                  do j3=n1,n2
                     cv_full(j2,j3) = &
                          cv_full(j2,j3)+cv_k(k)*f(index,j2,j3,ichemspec(k))
                  enddo
               enddo
            endif
         enddo
        slice = cp_full(m1:m2,n1:n2)/cv_full(m1:m2,n1:n2)
      elseif (dir == 2) then
         allocate (cp_full(mx,mz))  
         allocate (cv_full(mx,mz))  
         cp_full = 0.
         cp_full(l1:l2,n1:n2) = f(l1:l2,index,n1:n2,icp)
         cv_full = 0.
         do k=1,nchemspec
            if (species_constants(k,imass)>0.) then
               do j2=l1,l2
                  do j3=n1,n2
                     cv_full(j2,j3) = &
                          cv_full(j2,j3)+cv_k(k)*f(j2,index,j3,ichemspec(k))
                  enddo
               enddo
            endif
         enddo
        slice = cp_full(l1:l2,n1:n2)/cv_full(l1:l2,n1:n2)
      elseif (dir == 3) then
         allocate (cp_full(mx,my))  
         allocate (cv_full(mx,my))  
         cp_full = 0.
         cp_full(l1:l2,m1:m2) = f(l1:l2,m1:m2,index,icp)
         cv_full = 0.
         do k=1,nchemspec
            if (species_constants(k,imass)>0.) then
               do j2=l1,l2
                  do j3=m1,m2
                     cv_full(j2,j3) = &
                          cv_full(j2,j3)+cv_k(k)*f(j2,j3,index,ichemspec(k))
                  enddo
               enddo
            endif
         enddo
        slice = cp_full(l1:l2,m1:m2)/cv_full(l1:l2,m1:m2)
      else
        call fatal_error('get_gamma_slice','No such dir!')
      endif
!
      deallocate (cp_full)  
      deallocate (cv_full) 
!
    else
!
      if (dir == 1) then
         allocate (cp_full(my,mz))  
         allocate (cv_full(my,mz))  
         cp_full(m1:m2,n1:n2) = f(index,m1:m2,n1:n2,icp)
         cv_full = cp_full - f(index,:,:,iRR)
         slice = cp_full(m1:m2,n1:n2)/cv_full(m1:m2,n1:n2)
      elseif (dir == 2) then
         allocate (cp_full(mx,mz))  
         allocate (cv_full(mx,mz))  
         cp_full(l1:l2,n1:n2) = f(l1:l2,index,n1:n2,icp)
         cv_full = cp_full - f(:,index,:,iRR)
         slice = cp_full(l1:l2,n1:n2)/cv_full(l1:l2,n1:n2)
      elseif (dir == 3) then
         allocate (cp_full(mx,my))  
         allocate (cv_full(mx,my))  
         cp_full(l1:l2,m1:m2) = f(l1:l2,m1:m2,index,icp)
         cv_full = cp_full - f(:,:,index,iRR)
         slice = cp_full(l1:l2,m1:m2)/cv_full(l1:l2,m1:m2)
      else
        call fatal_error('get_gamma_slice','No such dir!')
      endif
!
      deallocate (cp_full)  
      deallocate (cv_full) 
!
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
      character(len=80) :: ChemInpLine
      character(len=10) :: specie_string
      character(len=1) :: tmp_string
      integer :: i, j, k=1
      real :: YY_k, air_mass, TT=300.
      real :: velx=0.
      real, dimension(nchemspec) :: stor2
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
      if (airdat .and. airin) then
        call fatal_error('chemistry',&
            'air.in and air.dat found. Please decide for one')
      endif
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
          if (lroot) print*, ' Pressure, Pa   ', PP

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
            if (lroot) print*, ' volume fraction, %,    ', YY_k, &
                species_constants(ind_chem,imass)

            if (species_constants(ind_chem,imass)>0.) then
             air_mass=air_mass+YY_k*0.01/species_constants(ind_chem,imass)
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
      if (emptyFile)  call stop_it('The input file air.dat was empty!')
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

      if (mvar < 5) then
        call fatal_error("air_field", "I can only set existing fields")
      endif

      if (.not. reinitialize_chemistry) then
      if (.not.lflame_front .or. .not.lflame_front_2D)  then
        if (ltemperature_nolog) then
          f(:,:,:,iTT)=TT
        else
          f(:,:,:,ilnTT)=alog(TT)!+f(:,:,:,ilnTT)
        endif
        if (ldensity_nolog) then
          f(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
        else
          tmp=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
            f(:,:,:,ilnrho)=alog(tmp)
        endif
        if (nxgrid>1) f(:,:,:,iux)=f(:,:,:,iux)+init_ux
      endif
      endif

      if (init_from_file) then
        if (lroot) print*, 'Velocity field read from file, initialization' //&
            'of density and temperature'
        if (ldensity_nolog) then
          f(:,:,:,ilnrho)=f(:,:,:,ilnrho)*(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
        else
          tmp=f(:,:,:,ilnrho)*(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
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
        if (x(i)<=init_x1) then
          f(i,:,:,ilnTT)=alog(init_TT1)
        endif
        if (x(i)>=init_x2) then
          f(i,:,:,ilnTT)=alog(init_TT2)
        endif
        if (x(i)>init_x1 .and. x(i)<init_x2) then
          if (init_x1 /= init_x2) then
            f(i,:,:,ilnTT)=&
               alog((x(i)-init_x1)/(init_x2-init_x1) &
               *(init_TT2-init_TT1)+init_TT1)
          endif
        endif
        enddo
      endif

      if (linit_density) then
        if (ldensity_nolog) then
          f(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/exp(f(:,:,:,ilnTT)))/unit_mass*unit_length**3
        else
          tmp=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/exp(f(:,:,:,ilnTT)))/unit_mass*unit_length**3
          f(:,:,:,ilnrho)=alog(tmp)
        endif
      endif

      if (nxgrid>1) then
        f(:,:,:,iux)=f(:,:,:,iux)+velx
      endif


      if (lroot) print*, 'Air temperature, K', TT
      if (lroot) print*, 'Air pressure, dyn', PP
      if (lroot) print*, 'Air density, g/cm^3:'
      if (lroot) print '(E10.3)',  PP/(k_B_cgs/m_u_cgs)*air_mass/TT
      if (lroot) print*, 'Air mean weight, g/mol', air_mass
      if (lroot) print*, 'R', k_B_cgs/m_u_cgs

      close(file_id)
!
    endsubroutine air_field
!***********************************************************************
!          NSCBC boundary conditions
!***********************************************************************
!    subroutine damp_zone_for_NSCBC(f,df) was here
!***********************************************************************
    subroutine get_mu1_slice(f,slice,grad_slice,index,sgn,direction)
!
! For the NSCBC boudary conditions the slice of mu1 at the boundary, and
! its gradient, is required.
!
!  10-dez-09/Nils Erland L. Haugen: coded
!
      use Deriv, only: der_onesided_4_slice_other
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(:,:), intent(out) :: slice
      real, dimension(:,:), intent(out) :: grad_slice
      integer, intent(in) :: index, sgn, direction
      intent(in) :: f
!
      integer :: j2, j3, k
!
      if (direction == 1) then
        slice = f(index,m1:m2,n1:n2,iRR)/Rgas
      elseif (direction == 2) then
        slice = f(l1:l2,index,n1:n2,iRR)/Rgas
      else
        slice = f(l1:l2,m1:m2,index,iRR)/Rgas
      endif
      call der_onesided_4_slice_other(f(:,:,:,iRR),sgn,grad_slice,index,direction)
      grad_slice = grad_slice/Rgas
!
    endsubroutine get_mu1_slice
!***********************************************************************
    subroutine get_reac_rate(f,p)
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
!      f(l1:l2,m,n,ireac_CO2) = p%DYDt_reac(:,ichem_CO2)
      f(l1:l2,m,n,ireac_CO) = p%DYDt_reac(:,ichem_CO)
!      f(l1:l2,m,n,ireac_O2) = p%DYDt_reac(:,ichem_O2)
!      f(l1:l2,m,n,ireac) = p%DYDt_reac(:,ichem_CO2)+p%DYDt_reac(:,ichem_CO)+p%DYDt_reac(:,ichem_O2)
!
    endsubroutine get_reac_rate
!***********************************************************************
    subroutine chemistry_clean_up()
!
      if (allocated(stoichio))       deallocate(stoichio)
      if (allocated(Sijm))           deallocate(Sijm)
      if (allocated(Sijp))           deallocate(Sijp)
      if (allocated(reaction_name))  deallocate(reaction_name)
      if (allocated(B_n))            deallocate(B_n)
      if (allocated(alpha_n))        deallocate(alpha_n)
      if (allocated(E_an))           deallocate(E_an)
      if (allocated(low_coeff))      deallocate(low_coeff)
      if (allocated(high_coeff))     deallocate(high_coeff)
      if (allocated(troe_coeff))     deallocate(troe_coeff)
      if (allocated(a_k4))           deallocate(a_k4)
      if (allocated(Mplus_case))     deallocate(Mplus_case)
      if (allocated(photochem_case)) deallocate(photochem_case)
      if (allocated(net_react_m))    deallocate(net_react_m)
      if (allocated(net_react_p))    deallocate(net_react_p)
      if (allocated(back))           deallocate(back)
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
! subroutine chemspec_normalization_N2(f) was here
!***********************************************************************
    subroutine getmu_array(f,mu1_full)
!
!  Calculate mean molecular weight
!
!  16-mar-10/natalia
!  30-jun-17/MR: moved here from eos_chemistry.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: mu1_full
      integer :: k,j2,j3
!
!  Mean molecular weight
!
      mu1_full=0.
      do k=1,nchemspec
        if (species_constants(k,imass)>0.) then
          do j2=mm1,mm2
            do j3=nn1,nn2
              mu1_full(:,j2,j3)= &
                  mu1_full(:,j2,j3)+f(:,j2,j3,ichemspec(k)) &
                  /species_constants(k,imass)
            enddo
          enddo
        endif
      enddo
!
    endsubroutine getmu_array
!***********************************************************************
!   subroutine read_Lewis
!
! Removed. Lewis # read from start.in
!
 !   endsubroutine read_Lewis
!***********************************************************************
   subroutine read_transport_data
!
!  Reading of the chemkin transport data
!
!  01-apr-08/natalia: coded
!  30-jun-17/MR: moved here from eos_chemistry.
!
     use Mpicomm, only: stop_it
!
      logical :: emptyfile
      logical :: found_specie
      integer :: file_id=123, ind_glob, ind_chem
      character (len=80) :: ChemInpLine
      character (len=10) :: specie_string
      integer :: VarNumber
      integer :: StartInd,StopInd,StartInd_1,StopInd_1
      logical :: tranin=.false.
      logical :: trandat=.false.
!
      emptyFile=.true.
!
      StartInd_1=1; StopInd_1 =0

      inquire (file='tran.dat',exist=trandat)
      inquire (file='tran.in',exist=tranin)
      if (tranin .and. trandat) then
        call fatal_error('eos_chemistry',&
            'tran.in and tran.dat found. Please decide which one to use.')
      endif

      if (tranin) open(file_id,file='tran.in')
      if (trandat) open(file_id,file='tran.dat')
!
      if (lroot) print*, 'the following species are found '//&
          'in tran.in/dat: beginning of the list:'
!
      dataloop: do
!
        read(file_id,'(80A)',end=1000) ChemInpLine(1:80)
        emptyFile=.false.
!
        StopInd_1=index(ChemInpLine,' ')
        specie_string=trim(ChemInpLine(1:StopInd_1-1))
!
        call find_species_index(specie_string,ind_glob,ind_chem,found_specie)
!
        if (found_specie) then
          if (lroot) print*,specie_string,' ind_glob=',ind_glob,' ind_chem=',ind_chem
!
          VarNumber=1; StartInd=1; StopInd =0
          stringloop: do while (VarNumber<7)
!
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
            StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
!
            if (StopInd==StartInd) then
              StartInd=StartInd+1
            else
              if (VarNumber==1) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E1.0)')  &
                    tran_data(ind_chem,VarNumber)
              elseif (VarNumber==2) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)')  &
                    tran_data(ind_chem,VarNumber)
              elseif (VarNumber==3) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)')  &
                    tran_data(ind_chem,VarNumber)
              elseif (VarNumber==4) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)')  &
                    tran_data(ind_chem,VarNumber)
              elseif (VarNumber==5) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)')  &
                    tran_data(ind_chem,VarNumber)
              elseif (VarNumber==6) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)')  &
                    tran_data(ind_chem,VarNumber)
              else
                call stop_it("No such VarNumber!")
              endif
!
              VarNumber=VarNumber+1
              StartInd=StopInd
            endif
            if (StartInd==80) exit
          enddo stringloop
!
        endif
      enddo dataloop
!
! Stop if tran.dat is empty
!
!
1000  if (emptyFile)  call stop_it('The input file tran.dat was empty!')
!
      if (lroot) print*, 'the following species are found in tran.dat: end of the list:'                    
!
      close(file_id)
!
    endsubroutine read_transport_data
!***********************************************************************
    subroutine jacobn(f,jacob)
!
!   dummy routine
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,nchemspec,nchemspec) :: jacob
!
          call fatal_error('jacobn', &
              'this does not work for simplified chemistry!')
!
!      call keep_compiler_quiet(jacob)
      call keep_compiler_quiet(f)
!
    endsubroutine jacobn
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
    subroutine chemspec_normalization_N2(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: sum_Y !, sum_Y2
      integer :: k ,isN2, ichemsN2
      logical :: lsN2
!
      call find_species_index('N2', isN2, ichemsN2, lsN2 )
      sum_Y=0.0 !; sum_Y2=0.0
      do k=1,nchemspec
        if (k/=ichemsN2) then
          sum_Y=sum_Y+f(:,:,:,ichemspec(k))
        endif
      enddo
      f(:,:,:,isN2)=1.0-sum_Y
!
    endsubroutine chemspec_normalization_N2
!***********************************************************************
endmodule Chemistry
