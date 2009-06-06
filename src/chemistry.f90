! $Id$
!  This modules addes chemical species and reactions.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lchemistry = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED gTT(3); mu1; gamma; gamma1; gamma11; gradcp(3)
! PENCILS PROVIDED cv; cv1; cp; cp1; YY(nchemspec)
! PENCILS PROVIDED cs2; rho1gpp(3); gmu1(3); nu; gradnu(3); nu_art
! PENCILS PROVIDED DYDt_reac(nchemspec); DYDt_diff(nchemspec)
! PENCILS PROVIDED lambda; glnlambda(3); ghYrho(3); cvspec(nchemspec)
! PENCILS PROVIDED DYDt_reac(nchemspec); DYDt_diff(nchemspec)
!
!***************************************************************

module Chemistry

  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'chemistry.h'

  real :: Rgas, Rgas_unit_sys=1.
  real, dimension (mx,my,mz) :: cp_full,cv_full,mu1_full, nu_full, pp_full
  real, dimension (mx,my,mz) :: lambda_full, rho_full, nu_art_full=0.
  real, dimension (mx,my,mz,nchemspec) :: cvspec_full
  real, dimension (mx,my,mz,nchemspec) ::  cp_R_spec
  real, dimension (mx,my,mz) ::  TT_full
  real, dimension (mx,my,mz) :: hYrho_full, e_int_full

  real :: lambda_const=impossible
  real :: visc_const=impossible
  real :: diffus_const=impossible
  real :: Sc_number=0.7, Pr_number=0.7
  real :: Cp_const=impossible
  real :: Cv_const=impossible
  real :: init_x1=-0.2,init_x2=0.2
  real :: init_TT1=400, init_TT2=2400., init_ux
  real :: str_thick=0.02
  real :: init_pressure=10.13e5
!
  logical :: lone_spec=.false.
  logical :: lfix_Sc=.false., lfix_Pr=.false.
 
!
!  parameters related to chemical reactions
!
  logical :: lreactions=.true.
  logical :: ladvection=.true.
  logical :: ldiffusion=.true.

  logical :: lheatc_chemistry=.true.

  logical :: BinDif_simple=.false.
  logical :: visc_simple=.false.
  
  logical :: lfilter=.false.
  logical :: lkreactions_profile=.false.
  integer :: nreactions=0,nreactions1=0,nreactions2=0
  integer :: ll1,ll2,mm1,mm2,nn1,nn2
  real, dimension(2*nchemspec) :: kreactions_profile_width=0.


  integer :: mreactions
  integer, allocatable, dimension(:,:) :: stoichio,Sijm,Sijp
  real,    allocatable, dimension(:,:) :: kreactions_z
  real,    allocatable, dimension(:)   :: kreactions_m,kreactions_p
  character (len=30),allocatable, dimension(:) :: reaction_name
  logical :: lT_tanh=.false.
! 1step_test case 

    logical :: l1step_test=.false.
    integer :: ipr=2
    real :: Tc=440., Tinf=2000., beta=1.09

!
!  hydro-related parameters
!
  real, dimension(nchemspec) :: amplchemk=0.,amplchemk2=0.
  real, dimension(nchemspec) :: chem_diff_prefactor=1.,initial_massfractions
  real :: amplchem=1.,kx_chem=1.,ky_chem=1.,kz_chem=1.,widthchem=1.
  real :: chem_diff=0.
  character (len=labellen), dimension (ninit) :: initchem='nothing'
  character (len=labellen), dimension (2*nchemspec) :: kreactions_profile=''

  real, allocatable, dimension(:,:,:,:,:) :: Bin_Diff_coef
  real, dimension (mx,my,mz,nchemspec) :: Diff_full, Diff_full_add, XX_full
  real, dimension (mx,my,mz,nchemspec) :: species_viscosity
  real, dimension(nchemspec) :: nu_spec=0., mobility=1.
!
!  Chemkin related parameters
!
  logical :: lcheminp=.false., lchem_cdtc=.false.
  logical :: lmobility=.false.
  real, dimension(nchemspec,18) :: species_constants
  integer :: imass=1, iTemp1=2,iTemp2=3,iTemp3=4
  integer, dimension(7) :: iaa1,iaa2
  real, allocatable, dimension(:)  :: B_n, alpha_n, E_an
  real, allocatable, dimension(:,:) :: low_coeff,high_coeff,troe_coeff,a_k4
  logical, allocatable, dimension(:) :: Mplus_case
  real, dimension(nchemspec,7)     :: tran_data
  real, dimension (nx,nchemspec), save  :: S0_R
  real, dimension (mx,my,mz,nchemspec), save :: H0_RT
  real, dimension (mx,my,mz,nchemspec), save :: RHS_Y_full
  real, dimension (mx,my,mz), save :: RHS_T_full

! input parameters
  namelist /chemistry_init_pars/ &
      initchem, amplchem, kx_chem, ky_chem, kz_chem, widthchem, &
      amplchemk,amplchemk2, chem_diff,nu_spec, BinDif_simple, visc_simple, &
      lambda_const, visc_const,Cp_const,Cv_const,diffus_const,init_x1,init_x2, &
      init_TT1,init_TT2,init_ux,l1step_test,Sc_number,init_pressure,lfix_Sc, str_thick, &
      lfix_Pr,lT_tanh


! run parameters
  namelist /chemistry_run_pars/ &
      lkreactions_profile,kreactions_profile,kreactions_profile_width, &
      chem_diff,chem_diff_prefactor, nu_spec, ldiffusion, ladvection, &
      lreactions,lchem_cdtc,lheatc_chemistry, BinDif_simple, visc_simple, &
      lmobility,mobility, lfilter,lT_tanh
!
! diagnostic variables (need to be consistent with reset list below)
!
  integer :: idiag_Y1m=0        ! DIAG_DOC: $\left<Y_1\right>$
  integer :: idiag_Y2m=0        ! DIAG_DOC: $\left<Y_2\right>$
  integer :: idiag_Y3m=0        ! DIAG_DOC: $\left<Y_3\right>$
  integer :: idiag_Y4m=0        ! DIAG_DOC: $\left<Y_4\right>$
  integer :: idiag_Y5m=0        ! DIAG_DOC: $\left<Y_5\right>$
  integer :: idiag_Y6m=0        ! DIAG_DOC: $\left<Y_6\right>$
  integer :: idiag_Y7m=0        ! DIAG_DOC: $\left<Y_7\right>$
  integer :: idiag_Y8m=0        ! DIAG_DOC: $\left<Y_8\right>$
  integer :: idiag_Y9m=0        ! DIAG_DOC: $\left<Y_9\right>$
  integer :: idiag_dY1m=0        ! DIAG_DOC: $\left<dY_1\right>$
  integer :: idiag_dY2m=0        ! DIAG_DOC: $\left<dY_2\right>$
  integer :: idiag_dY3m=0        ! DIAG_DOC: $\left<dY_3\right>$
  integer :: idiag_dY4m=0        ! DIAG_DOC: $\left<dY_4\right>$
  integer :: idiag_dY5m=0        ! DIAG_DOC: $\left<dY_5\right>$
  integer :: idiag_dY6m=0        ! DIAG_DOC: $\left<dY_6\right>$
  integer :: idiag_dY7m=0        ! DIAG_DOC: $\left<dY_7\right>$
  integer :: idiag_dY8m=0        ! DIAG_DOC: $\left<dY_8\right>$
  integer :: idiag_dY9m=0        ! DIAG_DOC: $\left<dY_9\right>$

  integer :: idiag_Y1mz=0        ! DIAG_DOC: $\left<Y_1\right>_{xy}(z)$
  integer :: idiag_Y2mz=0        ! DIAG_DOC: $\left<Y_2\right>_{xy}(z)$
  integer :: idiag_Y3mz=0        ! DIAG_DOC: $\left<Y_3\right>_{xy}(z)$
  integer :: idiag_Y4mz=0        ! DIAG_DOC: $\left<Y_4\right>_{xy}(z)$
  integer :: idiag_Y5mz=0        ! DIAG_DOC: $\left<Y_5\right>_{xy}(z)$
  integer :: idiag_Y6mz=0        ! DIAG_DOC: $\left<Y_6\right>_{xy}(z)$
  integer :: idiag_Y7mz=0        ! DIAG_DOC: $\left<Y_7\right>_{xy}(z)$
  integer :: idiag_Y8mz=0        ! DIAG_DOC: $\left<Y_8\right>_{xy}(z)$
  integer :: idiag_Y9mz=0        ! DIAG_DOC: $\left<Y_9\right>_{xy}(z)$

  integer :: idiag_h1m=0
  integer :: idiag_h2m=0
  integer :: idiag_h3m=0
  integer :: idiag_h4m=0
  integer :: idiag_h5m=0
  integer :: idiag_h6m=0
  integer :: idiag_h7m=0
  integer :: idiag_h8m=0
  integer :: idiag_h9m=0
  
  integer :: idiag_cpfull=0
  integer :: idiag_cvfull=0

  integer :: idiag_cp1m=0
  integer :: idiag_cp2m=0
  integer :: idiag_cp3m=0
  integer :: idiag_cp4m=0
  integer :: idiag_cp5m=0
  integer :: idiag_cp6m=0
  integer :: idiag_cp7m=0
  integer :: idiag_cp8m=0
  integer :: idiag_cp9m=0
  integer :: idiag_e_intm=0
!
  contains

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
      use Cdata
      use FArrayManager
      use General, only: chn
!
      integer :: k, ichemspec_tmp
      character (len=20) :: input_file='chem.inp'
!
!  Initialize some index pointers
!
      iaa1(1)=5;iaa1(2)=6;iaa1(3)=7;iaa1(4)=8
      iaa1(5)=9;iaa1(6)=10;iaa1(7)=11
!
      iaa2(1)=12;iaa2(2)=13;iaa2(3)=14;iaa2(4)=15
      iaa2(5)=16;iaa2(6)=17;iaa2(7)=18
!
!  Set ichemistry to consecutive numbers nvar+1, nvar+2, ..., nvar+nchemspec.
!
      call farray_register_pde('chemspec',ichemspec_tmp,vector=nchemspec)
      do k=1,nchemspec
        ichemspec(k)=ichemspec_tmp+k-1
      enddo
!
!  Read species to be used from chem.inp (if the file exists).
!
      inquire(FILE=input_file, EXIST=lcheminp)
      if (lcheminp) call read_species(input_file)
!
!  Read data on the thermodynamical properties of the different species.
!  All these data are stored in the array species_constants.
!
      if (lcheminp) call read_thermodyn(input_file)
!
!  Write all data on species and their thermodynamics to file.
!
      if (lcheminp) call write_thermodyn()
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_chemistry
!***********************************************************************
    subroutine initialize_chemistry(f)

!  called by run.f90 after reading parameters, but before the time loop
!
!  13-aug-07/steveb: coded
!  19-feb-08/axel: reads in chemistry.dat file
!
      use Mpicomm, only: stop_it

      real, dimension (mx,my,mz,mfarray) :: f
      logical :: data_file_exit=.false.
      logical :: exist,exist1,exist2
      character (len=15) :: file1='chemistry_m.dat',file2='chemistry_p.dat'
!
!  initialize chemistry
!
      if (lcheminp) then
        if (unit_temperature /= 1) then
          call fatal_error('initialize_chemistry', &
              'unit_temperature must be unity when using chemistry!')
        endif
!
!  calculate universal gas constant based on Boltzmann constant
!  and the proton mass
!
        if (unit_system == 'cgs') then
          Rgas_unit_sys = k_B_cgs/m_u_cgs
          Rgas=Rgas_unit_sys/unit_energy
        endif
      endif
!
      if (nchemspec==1) then
        lone_spec=.true.
        lreactions=.false.
        ladvection=.false.
        ldiffusion=.false.
      endif
!
!  check for the existence of chemistry input files
!
      inquire(file=file1,exist=exist1)
      inquire(file=file2,exist=exist2)
      inquire(file='chemistry.dat',exist=exist)
!
!  Read in data file in ChemKin format
!
      if (lcheminp) then
        call chemkin_data(f)
        data_file_exit=.true.
      endif
!
!  Alternatively, read in stoichiometric matrices in explicit format.
!  For historical reasons this is referred to as "astrobiology_data"
!
      if (exist1 .and. exist2) then
        call astrobiology_data(f)
        data_file_exit=.true.
      endif
!
      if (exist) then
        call astrobiology_data(f)
        data_file_exit=.true.
      endif
!
! check the existence of a data file
!
      if (.not. data_file_exit) then
        call stop_it('initialize_chemistry: there is no chemistry data file')
      endif
!

      if ((nxgrid==1) .and. (nygrid==1) .and. (nzgrid==1)) then
       ll1=1
       ll2=mx
       mm1=m1
       mm2=m2
       nn1=n1
       nn2=n2
      else
      if (nxgrid==1) then
       ll1=l1
       ll2=l2
      else
       ll1=1
       ll2=mx
      endif

      if (nygrid==1) then
       mm1=m1
       mm2=m2
      else
       mm1=1
       mm2=my
      endif

      if (nzgrid==1) then
       nn1=n1
       nn2=n2
      else
       nn1=1
       nn2=mz
      endif

     endif


      call keep_compiler_quiet(f)
!
    endsubroutine initialize_chemistry
!***********************************************************************
    subroutine init_chemistry(f)
!
!  initialise chemistry initial condition; called from start.f90
!
!  13-aug-07/steveb: coded
!
      use Cdata
      use Initcond
!AB:  use InitialCondition, only: initial_condition_chemistry
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j,k
      logical :: lnothing, air_exist
!
      intent(inout) :: f
!
!  different initializations of nd (called from start)
!
      lnothing=.false.
      do j=1,ninit
        select case(initchem(j))
!
        case('nothing')
          if (lroot .and. .not. lnothing) print*, 'init_chem: nothing '
          lnothing=.true.
        case('constant')
          do k=1,nchemspec
            f(:,:,:,ichemspec(k))=amplchemk(k)
          enddo
        case('positive-noise')
          do k=1,nchemspec
            call posnoise(amplchemk(k),f,ichemspec(k))
          enddo
        case('innerbox')
          do k=1,nchemspec
            call innerbox(amplchemk(k),amplchemk2(k),f,ichemspec(k),widthchem)
          enddo
        case('cos2x_cos2y_cos2z')
          do k=1,nchemspec
            call cos2x_cos2y_cos2z(amplchemk(k),f,ichemspec(k))
          enddo
        case('coswave-x')
          do k=1,nchemspec
            call coswave(amplchem,f,ichemspec(k),kx=kx_chem)
          enddo
        case('gaussian-x')
          do k=1,nchemspec
            call gaussian(amplchem,f,ichemspec(k),kx=kx_chem)
          enddo
        case('hatwave-x')
          do k=1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),widthchem,kx=kx_chem)
          enddo
        case('hatwave-y')
          do k=1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),widthchem,ky=ky_chem)
          enddo
        case('hatwave-z')
          do k=1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),widthchem,kz=kz_chem)
          enddo
        case('air')
          if (lroot ) print*, 'init_chem: air '
           inquire(file='air.dat',exist=air_exist)
           if (air_exist) then
            call air_field(f)
           else
            call stop_it('there is no air.dat file')
           endif
        case('flame_front')
          call flame_front(f)
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
!AB:  if (linitial_condition) call initial_condition_chemistry(f)
!
      if (lone_spec) then
        f(:,:,:,ichemspec(1))=1.
        if (lroot) print*, 'initchem: this is one specie case'
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_chemistry
!***********************************************************************
    subroutine pencil_criteria_chemistry()
!
!  All pencils that this chemistry module depends on are specified here.
!
!  13-aug-07/steveb: coded
!
      lpenc_requested(i_YY)=.true.
      lpenc_requested(i_cs2)=.true.
      lpenc_requested(i_cvspec)=.true.
!
      lpenc_requested(i_DYDt_reac)=.true.
      lpenc_requested(i_DYDt_diff)=.true.
!
       if (lcheminp) then
         lpenc_requested(i_mu1)=.true.
         lpenc_requested(i_gmu1)=.true.
         lpenc_requested(i_pp)=.true.
         lpenc_requested(i_cp)=.true.
         lpenc_requested(i_cp1)=.true.
         lpenc_requested(i_gradcp)=.true.
         lpenc_requested(i_cv)=.true.
         lpenc_requested(i_cv1)=.true.
         lpenc_requested(i_gamma)=.true.
         lpenc_requested(i_gamma1)=.true.
         lpenc_requested(i_gamma11)=.true.
         lpenc_requested(i_ghYrho)=.true.
!
         lpenc_requested(i_DYDt_reac)=.true.
         lpenc_requested(i_DYDt_diff)=.true.
         if (lheatc_chemistry) then
           lpenc_requested(i_lambda)=.true.
           lpenc_requested(i_glnlambda)=.true.
         endif
!
      endif
!
    endsubroutine pencil_criteria_chemistry
!***********************************************************************
    subroutine pencil_interdep_chemistry(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here
!
!  13-aug-07/steveb: coded
!
      use Sub, only: keep_compiler_quiet
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_chemistry
!***********************************************************************
    subroutine calc_pencils_chemistry(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   13-aug-07/steveb: coded
!
      use Cdata
      use Sub
      use Cparam
      use EquationOfState
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx) :: mu1_cgs, cp_spec
      real, dimension (mx) :: tmp_sum, tmp_sum2
      real, dimension (nchemspec,nchemspec) :: Phi
!
      intent(in) :: f
      intent(inout) :: p
      integer :: k,i,j
      real :: T_local, T_up, T_mid, T_low, tmp,  lnT_local
      real :: mk_mj, nuk_nuj
      logical :: lcheminp_tmp=.false.
!
!  Mass fraction YY
!
      if (lpencil(i_YY)) then
        do k=1,nchemspec;  p%YY(:,k)=f(l1:l2,m,n,ichemspec(k)); enddo
      endif
!
      if (lcheminp) then
!
!  Mean molecular weight
!
        if (lpencil(i_mu1)) then 
          p%mu1=mu1_full(l1:l2,m,n)
        endif
!
        if (lpencil(i_gmu1)) call grad(mu1_full,p%gmu1)
!
!  Mole fraction XX
!
        !if (lpencil(i_XX)) then
        !  do k=1,nchemspec 
        !    p%XX(:,k)=p%YY(:,k)/species_constants(ichemspec(k),imass)/p%mu1
        !  enddo
        !endif
!
!  Temperature
!
        if (lpencil(i_lnTT)) p%lnTT=f(l1:l2,m,n,ilnTT)
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=1./p%TT!
!
!  Temperature laplacian and gradient
!
        if (lpencil(i_glnTT)) call grad(f,ilnTT,p%glnTT)
        if (lpencil(i_del2lnTT)) call del2(f,ilnTT,p%del2lnTT)
!
!  Density
!
        !if (lpencil(i_lnrho)) p%lnrho=f(l1:l2,m,n,ilnrho)
        !if (lpencil(i_rho))   p%rho=exp(p%lnrho)
!
!  Pressure
!
        if (lpencil(i_pp)) p%pp = Rgas*p%TT*p%mu1*p%rho
!
!  Specific heat at constant pressure
!
        if (lpencil(i_cp)) then
          p%cp=cp_full(l1:l2,m,n)
        endif
!
        if (lpencil(i_cp1) .and. maxval(p%cp)>0)   p%cp1 = 1./p%cp
!
!  Gradient of the above
!
        if (lpencil(i_gradcp)) call grad(cp_full,p%gradcp)
!
!  Specific heat at constant volume (i.e. density)
!
        if (lpencil(i_cv)) p%cv = cv_full(l1:l2,m,n)
!
        if (lpencil(i_cv1)) p%cv1=1./p%cv
!
        if (lpencil(i_cvspec)) then
          do k=1,nchemspec
            p%cvspec(:,k)=cvspec_full(l1:l2,m,n,k)
          enddo
        endif
!

!
!  Polytropic index
!
        if (lpencil(i_gamma)) p%gamma = p%cp*p%cv1
        if (lpencil(i_gamma11)) p%gamma11 = p%cv*p%cp1
        if (lpencil(i_gamma1)) p%gamma1 = p%gamma - 1
!
!  Sound speed
!
        !if (lpencil(i_cs2)) p%cs2=p%cp*p%TT*p%gamma1*p%mu1
!     
!------------------------------------------
!  NB:this is wrong for the chemistry case
!------------------------------------------
!  WL: what is wrong?
!
!  Logarithmic pressure gradient
!
        if (lpencil(i_rho1gpp)) then
!
! NILS: rho1gpp should be calculated from gradT, gradlnrho and gradmu 
! NILS: instead. When this is implemented one should remove the 
! NILS: calculation of pp_full
!
          call grad(pp_full,p%rho1gpp)
!
          do i=1,3
            p%rho1gpp(:,i)=p%rho1gpp(:,i)/p%rho(:)
          enddo
!
          !do i=1,3
          !  p%rho1gpp(:,i) = p%gamma11*p%cs2*&
          !      (p%glnrho(:,i)+p%glnTT(:,i)+p%gmu1(:,i)/p%mu1(:))
          !enddo
        endif
!
!  Viscosity of a mixture
!
        if (lpencil(i_nu)) then
          p%nu=nu_full(l1:l2,m,n)
          if (lpencil(i_gradnu)) then
            call grad(nu_full,p%gradnu)
          endif
        endif
!
      endif
!
!  Artificial Viscosity of a mixture
!
      if (lpencil(i_nu_art)) then
        p%nu_art=nu_art_full(l1:l2,m,n)
      endif
! 
! Calculate the reaction term and the corresponding pencil 
!
      if (lreactions .and. lpencil(i_DYDt_reac)) then
        call calc_reaction_term(f,p)
      else
        p%DYDt_reac=0.
      endif
! 
! Calculate the diffusion term and the corresponding pencil 
!
      if (ldiffusion .and. lpencil(i_DYDt_diff)) then
        call calc_diffusion_term(f,p)
      else
        p%DYDt_diff=0.
      endif
!
      RHS_Y_full(l1:l2,m,n,:)=p%DYDt_reac+p%DYDt_diff
!
! Calculate thermal diffusivity
!
      if (lpenc_requested(i_lambda)) then
         p%lambda=lambda_full(l1:l2,m,n)
         if (lpenc_requested(i_glnlambda)) call grad(lambda_full,p%glnlambda)
      endif
!
!  Calculate grad(enthalpy)
!
      if (lpenc_requested(i_ghYrho)) then
        call grad(hYrho_full,p%ghYrho)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_chemistry
!**************************************************************************
subroutine flame_front(f)
!
! 06.05.2009/Nils Erland L. Haugen: adapted from similar 
!                                   routine in special/chem_stream.f90
! This routine set up the initial profiles used in 1D flame speed measurments
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: mu1
      integer :: k,j,i,j1,j2,j3

      real :: mO2, mH2, mN2, mH2O
      real :: log_inlet_density, del
      integer :: i_H2, i_O2, i_H2O, i_N2, ichem_H2, ichem_O2, ichem_N2, ichem_H2O
      real :: initial_mu1, final_massfrac_O2
      logical :: found_specie
!     
! Initialize some indexes
!
      call find_species_index('H2' ,i_H2 ,ichem_H2 ,found_specie)
      call find_species_index('O2' ,i_O2 ,ichem_O2 ,found_specie)
      call find_species_index('N2' ,i_N2 ,ichem_N2 ,found_specie)
      call find_species_index('H2O',i_H2O,ichem_H2O,found_specie)
      mO2 =species_constants(ichem_O2 ,imass)
      mH2 =species_constants(ichem_H2 ,imass)
      mH2O=species_constants(ichem_H2O,imass)
      mN2 =species_constants(ichem_N2 ,imass)
!      
! Find approximate value for the mass fraction of O2 after the flame front
!
      final_massfrac_O2&
          =(initial_massfractions(ichem_O2)/mO2&
          -initial_massfractions(ichem_H2)/(2*mH2))*mO2 
!
!  Initialize temperature and species
!
      do k=1,mx 
!
!  Initialize temperature
!

       if(lT_tanh) then
        del=init_x2-init_x1
         f(k,:,:,ilnTT)=log((init_TT2+init_TT1)*0.5  &
             +((init_TT2-init_TT1)*0.5)  &
             *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del)))
       else
         if (x(k)<init_x1) then
           f(k,:,:,ilnTT)=log(init_TT1)
         endif
         if (x(k)>init_x2) then
           f(k,:,:,ilnTT)=log(init_TT2)
         endif
         if (x(k)>init_x1 .and. x(k)<init_x2) then
           f(k,:,:,ilnTT)=log((x(k)-init_x1)/(init_x2-init_x1) &
               *(init_TT2-init_TT1)+init_TT1)
         endif
        endif
!
!  Initialize steam and hydrogen
!

       if (lT_tanh) then
         del=(init_x2-init_x1)/3.
         f(k,:,:,i_H2)=(0.+f(l1,:,:,i_H2))*0.5  &
           +(0.-f(l1,:,:,i_H2))*0.5  &
           *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
!
         f(k,:,:,i_H2O)=(f(l1,:,:,i_H2)/2.*18.+f(l1,:,:,i_H2O))*0.5  &
             +((f(l1,:,:,i_H2)/2.*18.-f(l1,:,:,i_H2O))*0.5)  &
             *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
!
       else
        if (x(k)>init_x1) then
          f(k,:,:,i_H2O)=initial_massfractions(ichem_H2)/mH2*mH2O &
               *(exp(f(k,:,:,ilnTT))-init_TT1) &
               /(init_TT2-init_TT1)
          f(k,:,:,i_H2)=initial_massfractions(ichem_H2) &
               *(exp(f(k,:,:,ilnTT))-init_TT2) &
               /(init_TT1-init_TT2)
        endif
       endif
!
!  Initialize oxygen
!

       if (lT_tanh) then
        del=(init_x2-init_x1)
        f(k,:,:,i_O2)=(f(l2,:,:,i_O2)+f(l1,:,:,i_O2))*0.5  &
             +((f(l2,:,:,i_O2)-f(l1,:,:,i_O2))*0.5)  &
             *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del))
       else

        if (x(k)>init_x2) then
          f(k,:,:,i_O2)=final_massfrac_O2
        endif
        if (x(k)>init_x1 .and. x(k)<init_x2) then
          f(k,:,:,i_O2)=(x(k)-init_x2)/(init_x1-init_x2) &
               *(initial_massfractions(ichem_O2)-final_massfrac_O2)&
               +final_massfrac_O2
        endif
       endif
      enddo

      call calc_for_chem_mixture(f)


!
!  Find logaritm of density at inlet
!
      initial_mu1&
          =initial_massfractions(ichem_H2)/(mH2)&
          +initial_massfractions(ichem_O2)/(mO2)&
          +initial_massfractions(ichem_H2O)/(mH2O)&
          +initial_massfractions(ichem_N2)/(mN2)
      log_inlet_density=log(init_pressure)-log(Rgas)-log(init_TT1)-log(initial_mu1)



       do j3=nn1,nn2
       do j2=mm1,mm2
       do j1=ll1,ll2

!
!  Initialize density
!
   
       f(j1,j2,j3,ilnrho)=log(init_pressure)-log(Rgas)  &
           -f(j1,j2,j3,ilnTT)-log(mu1_full(j1,j2,j3))
 
!
!  Initialize velocity
!
      f(j1,j2,j3,iux)=f(j1,j2,j3,iux)  &
            +init_ux*exp(log_inlet_density)/exp(f(j1,j2,j3,ilnrho))

       enddo
       enddo
       enddo


!
!  Check if we want nolog of density
!
      if (ldensity_nolog) f(:,:,:,ilnrho)=exp(f(:,:,:,ilnrho))
!
    endsubroutine flame_front
!***********************************************************************
    subroutine calc_for_chem_mixture(f)
!
!  DOCUMENT ME!!!
!
      use Cdata
      use Sub
      use Cparam
      use EquationOfState
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) ::  tmp_sum,tmp_sum2, nuk_nuj, nu_dyn, Phi
   !   real, dimension (mx,my,mz,nchemspec,nchemspec) :: Phi
      real, dimension (mx,my,mz,nchemspec) :: species_cond
!
      intent(in) :: f
      integer :: k,i,j, j1,j2,j3
      integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9
      real :: T_up, T_mid, T_low, T_loc, tmp
      real :: mk_mj
!
      logical :: tran_exist=.false.
      logical,save :: lwrite=.true.
!
      character (len=20) :: output_file="./data/mix_quant.out"
      integer :: file_id=123,lmid

! 
! Density and temperature
!
      rho_full=exp(f(:,:,:,ilnrho))
      TT_full=exp(f(:,:,:,ilnTT))
!
! Now this routine is only for chemkin data !!!
!
      if (lcheminp) then
!
        if (unit_system == 'cgs') then
!
          Rgas_unit_sys = k_B_cgs/m_u_cgs
          Rgas=Rgas_unit_sys/unit_energy
!
!  Mean molecular weight
!
          mu1_full=0.
          do k=1,nchemspec
           do j2=mm1,mm2
           do j3=nn1,nn2
            mu1_full(:,j2,j3)=mu1_full(:,j2,j3)+unit_mass*f(:,j2,j3,ichemspec(k)) &
                /species_constants(k,imass)
           enddo
           enddo
          enddo

          if (l1step_test) then
             species_constants(:,imass)=1.
             mu1_full=1.
          endif
!
!  Mole fraction XX
!
          do k=1,nchemspec 
           do j2=mm1,mm2
            do j3=nn1,nn2
!!$              if (minval(mu1_full(:,j2,j3))<=0) then
!!$                XX_full(:,j2,j3,k)=0.
!!$              else
                XX_full(:,j2,j3,k)=f(:,j2,j3,ichemspec(k))*unit_mass &
                 /(species_constants(k,imass)*mu1_full(:,j2,j3))
!!$              endif
            enddo
           enddo
          enddo
!
! NILS: Is this really necesarry?
!
          if (lpencil(i_rho1gpp)) then
           do j2=mm1,mm2
           do j3=nn1,nn2
            pp_full(:,j2,j3)=Rgas*mu1_full(:,j2,j3) &
                            *rho_full(:,j2,j3)*TT_full(:,j2,j3)
           enddo
           enddo
          endif
!
!  Specific heat at constant pressure
!
          cp_full=0.
          cv_full=0.
          do k=1,nchemspec
            T_low=species_constants(k,iTemp1)
            T_mid=species_constants(k,iTemp2)
            T_up= species_constants(k,iTemp3)
            
!
            do j3=nn1,nn2
              do j2=mm1,mm2
                do j1=1,mx

                 T_loc=TT_full(j1,j2,j3)
!
                  if (j1<=l1 .or. j2>=l2) then
                    T_low=0.
                    T_up=1e10
                  endif
!
                  if (j2<=m1 .or. j2>=m2) then
                    T_low=0.
                    T_up=1e10
                  endif
!
                  if (j3<=n1 .or. j3>=n2) then
                    T_low=0.
                    T_up=1e10
                  endif
!

                  if (T_loc >=T_low .and. T_loc <= T_mid) then
                    tmp=0. 
                    do j=1,5
                      tmp=tmp+species_constants(k,iaa2(j))*T_loc**(j-1) 
                    enddo
                    cp_R_spec(j1,j2,j3,k)=tmp
                    cvspec_full(j1,j2,j3,k)=cp_R_spec(j1,j2,j3,k)-1.
                  elseif (T_loc >=T_mid .and. T_loc<= T_up) then
                !  elseif (TT_full(j1,j2,j3) >=T_mid ) then
                    tmp=0.
                    do j=1,5 
                      tmp=tmp+species_constants(k,iaa1(j))*T_loc**(j-1) 
                    enddo
                    cp_R_spec(j1,j2,j3,k)=tmp
                    cvspec_full(j1,j2,j3,k)=cp_R_spec(j1,j2,j3,k)-1.
                  else
                    print*,'TT_full(j1,j2,j3)=',T_loc
                    print*,'j1,j2,j3=',j1,j2,j3
                    call fatal_error('calc_for_chem_mixture',&
                        'TT_full(j1,j2,j3) is outside range')
                  endif


                 cp_full(j1,j2,j3)=cp_full(j1,j2,j3)+f(j1,j2,j3,ichemspec(k))  &
                  *cp_R_spec(j1,j2,j3,k)/species_constants(k,imass)*Rgas
                 cv_full(j1,j2,j3)=cv_full(j1,j2,j3)+f(j1,j2,j3,ichemspec(k))  &
                  *cvspec_full(j1,j2,j3,k)/species_constants(k,imass)*Rgas


                enddo
              enddo
            enddo
!
          enddo


          if (Cp_const<impossible) then
            Cp_full=Cp_const
            Cv_full=Cp_const-Rgas
          endif


          if (Cv_const<impossible) then
            Cv_full=Cv_const
          endif
!
!  Binary diffusion coefficients
!
          inquire(file='tran.dat',exist=tran_exist)



          if (tran_exist) then
            call calc_diff_visc_coef(f)
          endif
!
!  Viscosity of a mixture
!

        do j3=nn1,nn2
        do j2=mm1,mm2
        do j1=1,mx

          if  (lone_spec) then
            nu_full(j1,j2,j3)=species_viscosity(j1,j2,j3,1)/rho_full(j1,j2,j3)
          elseif (visc_simple) then
            nu_dyn(j1,j2,j3)=0.
            do k=1,nchemspec 
              nu_dyn(j1,j2,j3)=nu_dyn(j1,j2,j3)+XX_full(j1,j2,j3,k) &
                              *species_viscosity(j1,j2,j3,k)
            enddo
            nu_full(j1,j2,j3)=nu_dyn(j1,j2,j3)/rho_full(j1,j2,j3)
          else 
            nu_dyn(j1,j2,j3)=0.
            do k=1,nchemspec
              tmp_sum2(j1,j2,j3)=0.
              do j=1,nchemspec
                mk_mj=species_constants(k,imass) &
                    /species_constants(j,imass)
                nuk_nuj(j1,j2,j3)=species_viscosity(j1,j2,j3,k) &
                    /species_viscosity(j1,j2,j3,j)
                Phi(j1,j2,j3)=1./sqrt(8.)*1./sqrt(1.+mk_mj) &
                    *(1.+sqrt(nuk_nuj(j1,j2,j3))*mk_mj**(-0.25))**2
                tmp_sum2(j1,j2,j3)=tmp_sum2(j1,j2,j3) &
                                  +XX_full(j1,j2,j3,j)*Phi(j1,j2,j3)
              enddo
              nu_dyn(j1,j2,j3)=nu_dyn(j1,j2,j3)+XX_full(j1,j2,j3,k)*&
                  species_viscosity(j1,j2,j3,k)/tmp_sum2(j1,j2,j3)
             enddo

              nu_full(j1,j2,j3)=nu_dyn(j1,j2,j3)/rho_full(j1,j2,j3)      
          endif

        enddo
        enddo
        enddo


          if (visc_const<impossible) then
                nu_full=visc_const
          endif
!
!  Diffusion coeffisient of a mixture
!
         do j3=nn1,nn2
         do j2=mm1,mm2
         do j1=1,mx

            if (.not. lone_spec) then

             if (diffus_const<impossible) then
                Diff_full(j1,j2,j3,:)=diffus_const
             elseif (lfix_Sc) then
              do k=1,nchemspec 
               Diff_full(j1,j2,j3,k)=species_viscosity(j1,j2,j3,k) &
                                    /rho_full(j1,j2,j3)/Sc_number
              enddo
             else
              do k=1,nchemspec
                tmp_sum(j1,j2,j3)=0.
                do j=1,nchemspec
                 if (j== k) then
                 else
                   tmp_sum(j1,j2,j3)=tmp_sum(j1,j2,j3) &
                        +XX_full(j1,j2,j3,j)/Bin_Diff_coef(j1,j2,j3,j,k)
                 endif	
                enddo
                 Diff_full(j1,j2,j3,k)=(1.-f(j1,j2,j3,ichemspec(k)))/tmp_sum(j1,j2,j3)
              enddo
             endif
            endif

            do k=1,nchemspec 
              Diff_full_add(j1,j2,j3,k)=Diff_full(j1,j2,j3,k)*&
                  species_constants(k,imass)/unit_mass*mu1_full(j1,j2,j3)
            enddo
         enddo
         enddo
         enddo
!
!do k=1,nchemspec 
!print*, Diff_full_add(l1,m1,n1,k),Diff_full_add(l2,m1,n1,k),k
!enddo
!
!
!
!  Thermal diffusivity 
!
!
! NB: one should check the coefficient 15/4
!
         do j3=nn1,nn2
         do j2=mm1,mm2
         do j1=1,mx

          tmp_sum(j1,j2,j3)=0.
          tmp_sum2(j1,j2,j3)=0.
!
          do k=1,nchemspec 
            species_cond(j1,j2,j3,k)=(species_viscosity(j1,j2,j3,k)) &
                /(species_constants(k,imass)/unit_mass)*Rgas*15./4.! 15./4.
            tmp_sum(j1,j2,j3)=tmp_sum(j1,j2,j3)  &
                             +XX_full(j1,j2,j3,k)*species_cond(j1,j2,j3,k)
            tmp_sum2(j1,j2,j3)=tmp_sum2(j1,j2,j3) &
                             +XX_full(j1,j2,j3,k)/species_cond(j1,j2,j3,k)

          enddo
         enddo
         enddo
         enddo

          do j1=1,mx
          do j2=mm1,mm2
          do j3=nn1,nn2
           if ((tmp_sum2(j1,j2,j3))<=0.) then
            lambda_full(j1,j2,j3)=0.
           else
            lambda_full(j1,j2,j3)=0.5*(tmp_sum(j1,j2,j3)+1./tmp_sum2(j1,j2,j3))
           endif

           if (lambda_const<impossible) then
            lambda_full(j1,j2,j3)=lambda_const
           endif


          enddo
          enddo
          enddo
 
         
!
!  Dimensionless Standard-state molar enthalpy H0/RT
!
          do k=1,nchemspec
            T_low=species_constants(k,iTemp1)
            T_mid=species_constants(k,iTemp2)
            T_up= species_constants(k,iTemp3)
            do j3=nn1,nn2
              do j2=mm1,mm2
                do j1=1,mx
                  if (TT_full(j1,j2,j3) <= T_mid) then
                    tmp=0. 
                    do j=1,5
                      tmp=tmp+species_constants(k,iaa2(j))*TT_full(j1,j2,j3)**(j-1)/j 
                    enddo
                    H0_RT(j1,j2,j3,k)=tmp+species_constants(k,iaa2(6))/TT_full(j1,j2,j3)
                  else               
                    tmp=0. 
                    do j=1,5
                      tmp=tmp+species_constants(k,iaa1(j))*TT_full(j1,j2,j3)**(j-1)/j 
                    enddo
                    H0_RT(j1,j2,j3,k)=tmp+species_constants(k,iaa1(6))/TT_full(j1,j2,j3)
                  endif
                enddo
              enddo
            enddo
          enddo

!................................................
!
!  Enthalpy flux
!

        do j3=nn1,nn2
        do j2=mm1,mm2
        do j1=1,mx

          hYrho_full(j1,j2,j3)=0.       
          do k=1,nchemspec
            hYrho_full(j1,j2,j3)=hYrho_full(j1,j2,j3)&
               +H0_RT(j1,j2,j3,k)*Rgas*TT_full(j1,j2,j3)&
               *f(j1,j2,j3,ichemspec(k))/species_constants(k,imass)
          enddo
!
!  Internal energy
!
          e_int_full(j1,j2,j3)=hYrho_full(j1,j2,j3) &
                    -Rgas*TT_full(j1,j2,j3)*mu1_full(j1,j2,j3)
          hYrho_full(j1,j2,j3)=hYrho_full(j1,j2,j3)*rho_full(j1,j2,j3)
        enddo
        enddo
        enddo

!
!.......................................................
!................................................
!  Enthalpy flux
!
   !       hYrho_full=0.       
   !       do k=1,nchemspec
   !         hYrho_full(l1:l2,m1:m2,n1:n2)=hYrho_full(l1:l2,m1:m2,n1:n2)&
   !             +H0_RT(l1:l2,m1:m2,n1:n2,k)&
   !             *Rgas*TT_full(l1:l2,m1:m2,n1:n2)&
   !             *f(l1:l2,m1:m2,n1:n2,ichemspec(k))/species_constants(k,imass)
   !       enddo
!
! Also the values at the ghost zones are required for hYrho_full
!
   !       do i=1,nghost
   !         hYrho_full(i,:,:)=hYrho_full(l1,:,:)
   !         hYrho_full(l2+i,:,:)=hYrho_full(l2,:,:)
   !         hYrho_full(:,i,:)=hYrho_full(:,m1,:)
   !         hYrho_full(:,m2+i,:)=hYrho_full(:,m2,:)
   !         hYrho_full(:,:,i)=hYrho_full(:,:,n1)
   !         hYrho_full(:,:,n2+i)=hYrho_full(:,:,n2)
   !       enddo
!
!  Internal energy
!
    !      e_int_full=hYrho_full-Rgas*TT_full*mu1_full
    !      hYrho_full=hYrho_full*rho_full
!
!.......................................................

!!$if (lroot) then
!!$lmid=(l2-l1)/2+l1
!!$print*,'e_int_full=',e_int_full(l1,m1,n1),e_int_full(lmid,m1,n1),e_int_full(l2,m1,n1)
!!$print*,'hYrho_full=',hYrho_full(l1,m1,n1),hYrho_full(lmid,m1,n1),hYrho_full(l2,m1,n1)
!!$print*,'H0_RT=',H0_RT(l1,m1,n1,1),H0_RT(lmid,m1,n1,1),H0_RT(l2,m1,n1,1)
!!$print*,'lambda_full=',lambda_full(l1,m1,n1),lambda_full(lmid,m1,n1),lambda_full(l2,m1,n1)
!!$print*,'Diff_full=',Diff_full(l1,m1,n1,1),Diff_full(lmid,m1,n1,1),Diff_full(l2,m1,n1,1)
!!$print*,'mu1_full=',mu1_full(l1,m1,n1),mu1_full(lmid,m1,n1),mu1_full(l2,m1,n1)
!!$print*,'nu_full=',nu_full(l1,m1,n1),nu_full(lmid,m1,n1),nu_full(l2,m1,n1)
!!$print*,'cp_full=',cp_full(l1,m1,n1),cp_full(lmid,m1,n1),cp_full(l2,m1,n1)
!!$print*,'cv_full=',cv_full(l1,m1,n1),cv_full(lmid,m1,n1),cv_full(l2,m1,n1)
!!$print*,'pp_full=',pp_full(l1,m1,n1),pp_full(lmid,m1,n1),pp_full(l2,m1,n1)
!!$print*,'xx_full=',xx_full(l1,m1,n1,1),xx_full(lmid,m1,n1,1),xx_full(l2,m1,n1,1)
!!$print*,'rho_full=',rho_full(l1,m1,n1),rho_full(lmid,m1,n1),rho_full(l2,m1,n1)
!!$print*,'TT_full=',TT_full(l1,m1,n1),TT_full(lmid,m1,n1),TT_full(l2,m1,n1)
!!$endif

        else
          call stop_it('This case works only for cgs units system!')
        endif
      endif
!
      if (lwrite) then
        open(file_id,file=output_file)
        write(file_id,*) 'Mixture quantities'
        write(file_id,*) '*******************'
        write(file_id,*) ''
        write(file_id,*) 'Mass, g/mole'
        write(file_id,'(7E12.4)') 1./maxval(mu1_full/unit_mass)
        write(file_id,*) ''
        write(file_id,*) 'Density, g/cm^3'
        write(file_id,'(7E12.4)') minval(rho_full)*unit_mass/unit_length**3, &
                                  maxval(rho_full)*unit_mass/unit_length**3
        write(file_id,*) ''
        write(file_id,*) 'Themperature, K'
         ! Commented the next line out because
         ! samples/2d-tests/chemistry_GrayScott apparently has no f(:,:,:,5)
         !write(file_id,'(7E12.4)') exp(maxval(f(:,:,:,5)))*unit_temperature
        write(file_id,*) ''
        write(file_id,*) 'Cp,  erg/mole/K'
        write(file_id,'(7E12.4)') maxval(cp_full)/Rgas*&
            Rgas_unit_sys/maxval(mu1_full/unit_mass)
        write(file_id,*) ''
        write(file_id,*) 'cp, erg/g/K'
        write(file_id,'(7E12.4)') maxval(cp_full)/Rgas*Rgas_unit_sys
        write(file_id,*) ''
        write(file_id,*) 'gamma,max,min'
        write(file_id,'(7E12.4)') maxval(cp_full(l1:l2,m1:m2,n1:n2)) &
                                 /maxval(cv_full(l1:l2,m1:m2,n1:n2)),&
          minval(cp_full(l1:l2,m1:m2,n1:n2))/minval(cv_full(l1:l2,m1:m2,n1:n2))
        write(file_id,*) ''
        write(file_id,*) 'Viscosity, g/cm/s,'
        write(file_id,'(7E12.4)') minval(nu_dyn)*&
            (unit_mass/unit_length/unit_time),maxval(nu_dyn)*&
            (unit_mass/unit_length/unit_time)
        write(file_id,*) ''
        write(file_id,*) 'Species viscosity, g/cm/s,'
        do k=1,nchemspec
        write(file_id,'(7E12.4)') minval(species_viscosity(:,:,:,k)),  &
                                  maxval(species_viscosity(:,:,:,k))
        enddo
        write(file_id,*) ''
        write(file_id,*) 'Thermal cond, erg/(cm K s),'
        write(file_id,'(7E12.4)') (lambda_full(l1,4,4)*&
            unit_energy/unit_time/unit_length/unit_temperature), &
                        (lambda_full(l2,4,4)*&
            unit_energy/unit_time/unit_length/unit_temperature)
        write(file_id,*) ''
        write(file_id,*) 'Species  Diffusion coefficient, cm^2/s'
        do k=1,nchemspec 
        write(file_id,'(7E12.4)')minval(Diff_full(:,:,:,k))*unit_length**2/unit_time, &
                                 maxval(Diff_full(:,:,:,k))*unit_length**2/unit_time
        enddo
        write(file_id,*) ''
      
       
        if (lroot) print*,'calc_for_chem_mixture: writing mix_quant.out file'
        close(file_id)
        lwrite=.false.
      endif
!

    endsubroutine calc_for_chem_mixture
!**********************************************************************
    subroutine astrobiology_data(f)
!
!  Proceedure to read in stoichiometric matrices in explicit format for
!  forward and backward reations. For historical reasons this is referred
!  to as "astrobiology_data".
!
!  28-feb-08/axel: coded
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub, only: keep_compiler_quiet
      use General, only: chn
!
      character (len=80) :: chemicals=''
      character (len=15) :: file1='chemistry_m.dat',file2='chemistry_p.dat'
      character (len=20) :: input_file='chem.inp'
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: exist,exist1,exist2
      integer :: i,j,k,stat,reac,spec
!
!  Find number of ractions
!
      mreactions=2*nchemspec
      if (lroot) print*,'Number of reactions=',mreactions
!
!  Allocate reaction arrays (but not during reloading!)
!
      if (.not.lreloading) then
        allocate(stoichio(nchemspec,mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for stoichio")
        allocate(Sijm(nchemspec,mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for Sijm")
        allocate(Sijp(nchemspec,mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for Sijp")
        allocate(kreactions_z(mz,mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for kreactions_z")
        allocate(kreactions_p(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for kreactions_p")
        allocate(kreactions_m(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for kreactions_m")
        allocate(reaction_name(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for reaction_name")
      endif
!
!  Initialize data
!
      kreactions_z=1.
      Sijp=0
      Sijm=0
!
!  read chemistry data
!
      inquire(file=file1,exist=exist1)
      inquire(file=file2,exist=exist2)
! 
      if (exist1.and.exist2) then
!
!  if both chemistry1.dat and chemistry2.dat are present,
!  then read Sijp and Sijm, and calculate their sum
!
!  file1
!
        open(19,file=file1)
        read(19,*) chemicals
        do j=1,mreactions
          read(19,*,end=994) kreactions_m(j),(Sijm(i,j),i=1,nchemspec)
        enddo
994     close(19)
        nreactions1=j-1
!
!  file2
!
        open(19,file=file2)
        read(19,*) chemicals
        do j=1,mreactions
          read(19,*,end=992) kreactions_p(j),(Sijp(i,j),i=1,nchemspec)
        enddo
992     close(19)
        nreactions2=j-1
!
!  calculate stoichio and nreactions
!
        if (nreactions1==nreactions2) then
          nreactions=nreactions1
          stoichio=Sijp-Sijm
        else
          call stop_it('nreactions1/=nreactions2')
        endif
!
      else 
!
!  old method: read chemistry data, if present
!
        inquire(file='chemistry.dat',exist=exist)
        if (exist) then
          open(19,file='chemistry.dat')
          read(19,*) chemicals
          do j=1,mreactions
            read(19,*,end=990) kreactions_p(j),(stoichio(i,j),i=1,nchemspec)
          enddo
990       close(19)
          nreactions=j-1
          Sijm=-min(stoichio,0)
          Sijp=+max(stoichio,0)
        else
          if (lroot) print*,'no chemistry.dat file to be read.'
          lreactions=.false.
        endif
      endif
!
!  print input data for verification
!
      if (lroot) then
        print*,'chemicals=',chemicals
        print*,'kreactions_m=',kreactions_m(1:nreactions)
        print*,'kreactions_p=',kreactions_p(1:nreactions)
        print*,'Sijm:' ; write(*,100),Sijm(:,1:nreactions)
        print*,'Sijp:' ; write(*,100),Sijp(:,1:nreactions)
        print*,'stoichio=' ; write(*,100),stoichio(:,1:nreactions)
      endif
!
!  possibility of z-dependent kreactions_z profile
!
      if (lkreactions_profile) then
        do j=1,nreactions
          if (kreactions_profile(j)=='cosh') then
            do n=1,mz
              kreactions_z(n,j)=1./cosh(z(n)/kreactions_profile_width(j))**2
            enddo
          elseif (kreactions_profile(j)=='square') then
            do n=1,mz
              if (n < mz/2) then
                kreactions_z(n,j)=kreactions_profile_width(j)
              else 
                kreactions_z(n,j)=0.
              endif
            enddo
          elseif (kreactions_profile(j)=='sin') then
            do n=1,mz
              kreactions_z(n,j)=0.5*(1+cos(pi*z(n)/kreactions_profile_width(j)))
            enddo
          endif
        enddo
      endif
!
100   format(8i4)
!
    endsubroutine astrobiology_data
!**********************************************************************
    subroutine chemkin_data(f)
!
!  if the file with chemkin data exists
! 
!  DOCUMENT ME FURTHER!!!
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub, only: keep_compiler_quiet
      use General, only: chn
!
      character (len=20) :: input_file='chem.inp'
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: stat,k
      logical :: tran_exist

      inquire(file='tran.dat',exist=tran_exist)
!
!  Allocate binary diffusion coefficient array
!
      if (.not.lreloading) then
        if (.not. lfix_Sc) then
          allocate(Bin_Diff_coef(mx,my,mz,nchemspec,nchemspec),STAT=stat)
          if (stat>0) call stop_it("Couldn't allocate memory "//&
              "for binary diffusion coefficients") 
        endif
      endif
!
      if (tran_exist) then 
        if (lroot) then
          print*,'tran.dat file with transport data is found.'
        endif
        call read_transport_data
        !call calc_diff_visc_coef(f)
      else
        if (lroot) then
          print*,'tran.dat file with transport data is not found.'
          print*,'Now diffusion coefficients is ',chem_diff
          print*,'Now species viscosity is ',nu_spec
        endif
        Bin_Diff_coef=chem_diff/(unit_length**2/unit_time)
        do k=1,nchemspec
          species_viscosity(:,:,:,k)=nu_spec(k)/&
              (unit_mass/unit_length/unit_time)
        enddo
      endif
!
!  Find number of ractions
!
      call read_reactions(input_file,NrOfReactions=mreactions)
      if (lroot) print*,'Number of reactions=',mreactions
      nreactions=mreactions
!
!  Allocate reaction arrays
!
      if (.not.lreloading) then
        allocate(stoichio(nchemspec,mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for stoichio")
        allocate(Sijm(nchemspec,mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for Sijm")
        allocate(Sijp(nchemspec,mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for Sijp")
        allocate(kreactions_p(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for kreactions_p")
        allocate(kreactions_m(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for kreactions_m")
        allocate(reaction_name(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for reaction_name")
        
        allocate(B_n(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for B_n")
        allocate(alpha_n(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for alpha_n")
        allocate(E_an(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for E_an")
        
        allocate(low_coeff(3,nreactions),STAT=stat)
        low_coeff=0.
        if (stat>0) call stop_it("Couldn't allocate memory for low_coeff")
        allocate(high_coeff(3,nreactions),STAT=stat)
        high_coeff=0.
        if (stat>0) call stop_it("Couldn't allocate memory for high_coeff")
        allocate(troe_coeff(3,nreactions),STAT=stat)
        troe_coeff=0.
        if (stat>0) call stop_it("Couldn't allocate memory for troe_coeff")
        allocate(a_k4(nchemspec,nreactions),STAT=stat)
        a_k4=impossible 
        if (stat>0) call stop_it("Couldn't allocate memory for troe_coeff")
        allocate(Mplus_case(nreactions),STAT=stat)
        Mplus_case=.false. 
        if (stat>0) call stop_it("Couldn't allocate memory for Mplus_case")
      endif
!
!  Initialize data
!
      Sijp=0
      Sijm=0
!
!  read chemistry data
!
      call read_reactions(input_file)
      call write_reactions()
!
!  calculate stoichio and nreactions
!
      stoichio=Sijp-Sijm
!
!  print input data for verification
!
      if (lroot .and. nreactions>0) then
        print*,'kreactions_m=',kreactions_m(1:nreactions)
        print*,'kreactions_p=',kreactions_p(1:nreactions)
        print*,'Sijm:' ; write(*,100),Sijm(:,1:nreactions)
        print*,'Sijp:' ; write(*,100),Sijp(:,1:nreactions)
        print*,'stoichio=' ; write(*,100),stoichio(:,1:nreactions)
      endif

100   format(8i4)

    endsubroutine chemkin_data
!**********************************************************************
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
!
      use Cdata
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: gchemspec
      real, dimension (nx) :: ugchemspec, sum_DYDT, ghYrho_uu=0.
      type (pencil_case) :: p
!
!  indices
!
      integer :: j,k,i
      integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9
!
      intent(in) :: f,p
      intent(inout) :: df
!
      real,dimension(nchemspec) :: reac_rate=0.
!
      real,dimension(nx) :: sum_reac_rate
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dchemistry_dt: SOLVE dchemistry_dt'
!     !if (headtt) call identify_bcs('ss',iss)
!
!  loop over all chemicals
!

      do k=1,nchemspec
!
!  advection terms
! 
        if (ladvection) then 
          call grad(f,ichemspec(k),gchemspec) 
          call dot_mn(p%uu,gchemspec,ugchemspec)
          if (lmobility) ugchemspec=ugchemspec*mobility(k)
          df(l1:l2,m,n,ichemspec(k))=df(l1:l2,m,n,ichemspec(k))-ugchemspec
        endif
!
!  diffusion operator
!
!  Temporary we check the existence of chem.imp data, 
!  further one should check the existence of a file with 
!  binary diffusion coefficients!
!
        if (ldiffusion) then
          df(l1:l2,m,n,ichemspec(k))=df(l1:l2,m,n,ichemspec(k))+&
              p%DYDt_diff(:,k) 
        endif
!
!  chemical reactions:
!  multiply with stoichiometric matrix with reaction speed
!  d/dt(x_i) = S_ij v_j
!
        if (lreactions) then
          df(l1:l2,m,n,ichemspec(k))=df(l1:l2,m,n,ichemspec(k))+&
              p%DYDt_reac(:,k)
        endif
!
      if (lfilter) then
       do i=1,mx
        if ((f(i,m,n,ichemspec(k))+df(i,m,n,ichemspec(k))*dt)<-1e-25 ) df(i,m,n,ichemspec(k))=-1e-25*dt
        if ((f(i,m,n,ichemspec(k))+df(i,m,n,ichemspec(k))*dt)>1. ) df(i,m,n,ichemspec(k))=1.*dt
       enddo
      endif


     enddo
!
      if (ldensity .and. lcheminp) then

       if (l1step_test) then
        sum_DYDt=0.
        do i=1,nx 
  
        sum_DYDt(i)=-p%rho(1)*(p%TT(i)-Tinf)/p%TT(i) &
           *Cp_const/lambda_const*beta*(beta-1.)*f(l1,m,n,iux)**2
  
      !   if (p%TT(i)>Tc) then
      !     if (x(i)>0.) then

    !        sum_DYDt(i)=f(l1,m,n,iux)**2*(Tinf-p%TT(i))/p%TT(i) & !(-p%TT(i)+Tinf)
    !        *Cp_const/lambda_const*p%rho(1)*beta*(beta-1.)* &
    !        (1.-f(l1-1+i,m,n,ichemspec(ipr)))
    
      
      !   endif
        enddo
       else
        sum_DYDt=0.
        do k=1,nchemspec
          sum_DYDt=sum_DYDt+Rgas/species_constants(k,imass)*&
              (1.-H0_RT(l1:l2,m,n,k))*(p%DYDt_reac(:,k)+p%DYDt_diff(:,k))
        enddo
       endif
!
       if (nxgrid >1) then
        call dot_mn(p%ghYrho,p%uu,ghYrho_uu)
       endif
!

        if (l1step_test) then
          RHS_T_full(l1:l2,m,n)=sum_DYDt(:)
        else
          RHS_T_full(l1:l2,m,n)=(sum_DYDt(:)- Rgas*p%mu1*p%divu)*p%cv1 &
            !/(p%cp-Rgas*p%mu1)&
            -(hYrho_full(l1:l2,m,n)*p%divu(:)+ghYrho_uu(:))/p%TT(:)*p%cv1/p%rho(:)
        endif

        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + RHS_T_full(l1:l2,m,n)

!print*,'nat2',maxval(RHS_T_full),minval(RHS_T_full)

!
        if (lheatc_chemistry) call calc_heatcond_chemistry(f,df,p)
!

      endif
!
!  For the timestep calculation, need maximum diffusion
!
      if (lfirst.and. ldt) then
        if (.not. lcheminp) then
          diffus_chem=chem_diff*maxval(chem_diff_prefactor)*dxyz_2
        else
          do j=1,nx
            if (ldiffusion) then
!
!--------------------------------------
!  This expression should be discussed 
!--------------------------------------
!
              diffus_chem(j)=diffus_chem(j)+&
                  maxval(Diff_full(l1+j-1,m,n,1:nchemspec))*dxyz_2(j)
            else
              diffus_chem(j)=0.
            endif
          enddo
        endif
      endif
!
! NB: it should be discussed
! WL: what should be discussed?
!
      if (lfirst .and. ldt) then
        if (lreactions) then
!
!  calculate maximum of *relative* reaction rate if decaying,
!  or maximum of absolute rate, if growing.
!
          if (lchem_cdtc) then
            reac_chem=0.
            do k=1,nchemspec
              reac_chem=max(reac_chem, &
                  abs(p%DYDt_reac(:,k)/max(p%YY(:,k),.001)))
            enddo
!
          elseif (lcheminp) then
            reac_chem=0.
            !sum_reac_rate=0.
            do k=1,nchemspec
              reac_chem=reac_chem+abs(p%DYDt_reac(:,k)/p%YY(:,k))
              !sum_reac_rate=sum_reac_rate+p%DYDt_reac(:,k)
            enddo
            if (maxval(reac_chem)>1e11) then
              reac_chem=1e11
            endif
          endif
        endif
      endif
!
!  Calculate diagnostic quantities
!
      if (ldiagnos) then
!
!  WL: instead of hardcoding Y1-Y9, wouldn't it be possible
!      to have them all in the same array? The particles_nbody 
!      module, for instance, has idiag_xxspar and idiag_vvspar, which 
!      allows the user to add output the positions and velocities 
!      of as many particle he/she wants. 
!
        if (idiag_Y1m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(i1)),idiag_Y1m)
        if (idiag_Y2m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(i2)),idiag_Y2m)
        if (idiag_Y3m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(i3)),idiag_Y3m)
        if (idiag_Y4m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(i4)),idiag_Y4m)
        if (idiag_Y5m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(i5)),idiag_Y5m)
        if (idiag_Y6m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(i6)),idiag_Y6m)
        if (idiag_Y7m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(i7)),idiag_Y7m)
        if (idiag_Y8m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(i8)),idiag_Y8m)
        if (idiag_Y9m/=0) call sum_mn_name(f(l1:l2,m,n,ichemspec(i9)),idiag_Y9m)
!
        if (idiag_cpfull/=0) call sum_mn_name(cp_full(l1:l2,m,n),idiag_cpfull)
        if (idiag_cvfull/=0) call sum_mn_name(cv_full(l1:l2,m,n),idiag_cvfull)
!
        if (idiag_cp1m/=0) call sum_mn_name(cp_R_spec(l1:l2,m,n,i1)*&
                             Rgas/species_constants(i1,imass),idiag_cp1m)
        if (idiag_cp2m/=0) call sum_mn_name(cp_R_spec(l1:l2,m,n,i2)*&
                             Rgas/species_constants(i2,imass),idiag_cp2m)
        if (idiag_cp3m/=0) call sum_mn_name(cp_R_spec(l1:l2,m,n,i3)*&
                             Rgas/species_constants(i3,imass),idiag_cp3m)
        if (idiag_cp4m/=0) call sum_mn_name(cp_R_spec(l1:l2,m,n,i4)*&
                             Rgas/species_constants(i4,imass),idiag_cp4m)
        if (idiag_cp5m/=0) call sum_mn_name(cp_R_spec(l1:l2,m,n,i5)*&
                             Rgas/species_constants(i5,imass),idiag_cp5m)
        if (idiag_cp6m/=0) call sum_mn_name(cp_R_spec(l1:l2,m,n,i6)*&
                             Rgas/species_constants(i6,imass),idiag_cp6m)
        if (idiag_cp7m/=0) call sum_mn_name(cp_R_spec(l1:l2,m,n,i7)*&
                             Rgas/species_constants(i7,imass),idiag_cp7m)
        if (idiag_cp8m/=0) call sum_mn_name(cp_R_spec(l1:l2,m,n,i8)*&
                             Rgas/species_constants(i8,imass),idiag_cp8m)
        if (idiag_cp9m/=0) call sum_mn_name(cp_R_spec(l1:l2,m,n,i9)*&
                             Rgas/species_constants(i9,imass),idiag_cp9m)
        if (idiag_e_intm/=0) call sum_mn_name(e_int_full(l1:l2,m,n),&
                             idiag_e_intm)
!
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1ddiagnos) then
        if (idiag_Y1mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(i1)),idiag_Y1mz)
        if (idiag_Y2mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(i2)),idiag_Y2mz)
        if (idiag_Y3mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(i3)),idiag_Y3mz)
        if (idiag_Y4mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(i4)),idiag_Y4mz)
        if (idiag_Y5mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(i5)),idiag_Y5mz)
        if (idiag_Y6mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(i6)),idiag_Y6mz)
        if (idiag_Y7mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(i7)),idiag_Y7mz)
        if (idiag_Y8mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(i8)),idiag_Y8mz)
        if (idiag_Y9mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,ichemspec(i9)),idiag_Y9mz)
      endif
!
!  Keep compiler quiet by ensuring every parameter is used
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dchemistry_dt
!***********************************************************************
    subroutine read_chemistry_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=chemistry_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=chemistry_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_chemistry_init_pars
!***********************************************************************
   subroutine write_chemistry_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=chemistry_init_pars)
!
    endsubroutine write_chemistry_init_pars
!***********************************************************************
    subroutine read_chemistry_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=chemistry_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=chemistry_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_chemistry_run_pars
!***********************************************************************
    subroutine write_chemistry_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=chemistry_run_pars)
!
    endsubroutine write_chemistry_run_pars
!***********************************************************************
    subroutine rprint_chemistry(lreset,lwrite)
!
!  reads and registers print parameters relevant to chemistry
!
!  13-aug-07/steveb: coded
!
      use Cdata
      use Diagnostics
      use General, only: chn
!
      integer :: iname,inamez
      logical :: lreset,lwr
      logical, optional :: lwrite
      character (len=5) :: schem,schemspec,snd1,smd1,smi1
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_Y1m=0; idiag_Y2m=0; idiag_Y3m=0; idiag_Y4m=0
        idiag_Y5m=0; idiag_Y6m=0; idiag_Y7m=0; idiag_Y8m=0
        idiag_Y9m=0
        idiag_dY1m=0; idiag_dY2m=0; idiag_dY3m=0; idiag_dY4m=0
        idiag_dY5m=0; idiag_dY6m=0; idiag_dY7m=0; idiag_dY8m=0
        idiag_dY9m=0
        idiag_h1m=0; idiag_h2m=0; idiag_h3m=0; idiag_h4m=0;
        idiag_h5m=0; idiag_h6m=0; idiag_h7m=0; idiag_h8m=0;
        idiag_h9m=0
        idiag_cp1m=0; idiag_cp2m=0; idiag_cp3m=0; idiag_cp4m=0;
        idiag_cp5m=0; idiag_cp6m=0; idiag_cp7m=0; idiag_cp8m=0; 
        idiag_cp9m=0; idiag_cpfull=0; idiag_cvfull=0
        idiag_e_intm=0
        idiag_Y1mz=0; idiag_Y2mz=0; idiag_Y3mz=0; idiag_Y4mz=0
        idiag_Y5mz=0; idiag_Y6mz=0; idiag_Y7mz=0; idiag_Y8mz=0
        idiag_Y9mz=0
      endif
!
      call chn(nchemspec,schemspec)
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Y1m',idiag_Y1m)
        call parse_name(iname,cname(iname),cform(iname),'Y2m',idiag_Y2m)
        call parse_name(iname,cname(iname),cform(iname),'Y3m',idiag_Y3m)
        call parse_name(iname,cname(iname),cform(iname),'Y4m',idiag_Y4m)
        call parse_name(iname,cname(iname),cform(iname),'Y5m',idiag_Y5m)
        call parse_name(iname,cname(iname),cform(iname),'Y6m',idiag_Y6m)
        call parse_name(iname,cname(iname),cform(iname),'Y7m',idiag_Y7m)
        call parse_name(iname,cname(iname),cform(iname),'Y8m',idiag_Y8m)
        call parse_name(iname,cname(iname),cform(iname),'Y9m',idiag_Y9m)
        call parse_name(iname,cname(iname),cform(iname),'dY1m',idiag_dY1m)
        call parse_name(iname,cname(iname),cform(iname),'dY2m',idiag_dY2m)
        call parse_name(iname,cname(iname),cform(iname),'dY3m',idiag_dY3m)
        call parse_name(iname,cname(iname),cform(iname),'dY4m',idiag_dY4m)
        call parse_name(iname,cname(iname),cform(iname),'dY5m',idiag_dY5m)
        call parse_name(iname,cname(iname),cform(iname),'dY6m',idiag_dY6m)
        call parse_name(iname,cname(iname),cform(iname),'dY7m',idiag_dY7m)
        call parse_name(iname,cname(iname),cform(iname),'dY8m',idiag_dY8m)
        call parse_name(iname,cname(iname),cform(iname),'dY9m',idiag_dY9m)
        call parse_name(iname,cname(iname),cform(iname),'h1m',idiag_h1m)
        call parse_name(iname,cname(iname),cform(iname),'h2m',idiag_h2m)
        call parse_name(iname,cname(iname),cform(iname),'h3m',idiag_h3m)
        call parse_name(iname,cname(iname),cform(iname),'h4m',idiag_h4m)
        call parse_name(iname,cname(iname),cform(iname),'h5m',idiag_h5m)
        call parse_name(iname,cname(iname),cform(iname),'h6m',idiag_h6m)
        call parse_name(iname,cname(iname),cform(iname),'h7m',idiag_h7m)
        call parse_name(iname,cname(iname),cform(iname),'h8m',idiag_h8m)
        call parse_name(iname,cname(iname),cform(iname),'h9m',idiag_h9m)
        call parse_name(iname,cname(iname),cform(iname),'cpfull',idiag_cpfull)
        call parse_name(iname,cname(iname),cform(iname),'cvfull',idiag_cvfull)
        call parse_name(iname,cname(iname),cform(iname),'cp1m',idiag_cp1m)
        call parse_name(iname,cname(iname),cform(iname),'cp2m',idiag_cp2m)
        call parse_name(iname,cname(iname),cform(iname),'cp3m',idiag_cp3m)
        call parse_name(iname,cname(iname),cform(iname),'cp4m',idiag_cp4m)
        call parse_name(iname,cname(iname),cform(iname),'cp5m',idiag_cp5m)
        call parse_name(iname,cname(iname),cform(iname),'cp6m',idiag_cp6m)
        call parse_name(iname,cname(iname),cform(iname),'cp7m',idiag_cp7m)
        call parse_name(iname,cname(iname),cform(iname),'cp8m',idiag_cp8m)
        call parse_name(iname,cname(iname),cform(iname),'cp9m',idiag_cp9m)
        call parse_name(iname,cname(iname),cform(iname),'e_intm',idiag_e_intm)
      enddo
!
!  xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Y1mz',idiag_Y1mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Y2mz',idiag_Y2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Y3mz',idiag_Y3mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Y4mz',idiag_Y4mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Y5mz',idiag_Y5mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Y6mz',idiag_Y6mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Y7mz',idiag_Y7mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Y8mz',idiag_Y8mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'Y9mz',idiag_Y9mz)
      enddo
!
!  Write chemistry index in short notation
!
      call chn(ichemspec(1),snd1)
      if (lwr) then
        write(3,*) 'i_Y1m=',idiag_Y1m
        write(3,*) 'i_Y2m=',idiag_Y2m
        write(3,*) 'i_Y3m=',idiag_Y3m
        write(3,*) 'i_Y4m=',idiag_Y4m
        write(3,*) 'i_Y5m=',idiag_Y5m
        write(3,*) 'i_Y6m=',idiag_Y6m
        write(3,*) 'i_Y7m=',idiag_Y7m
        write(3,*) 'i_Y8m=',idiag_Y8m
        write(3,*) 'i_Y9m=',idiag_Y9m
        write(3,*) 'i_dY1m=',idiag_dY1m
        write(3,*) 'i_dY2m=',idiag_dY2m
        write(3,*) 'i_dY3m=',idiag_dY3m
        write(3,*) 'i_dY4m=',idiag_dY4m
        write(3,*) 'i_dY5m=',idiag_dY5m
        write(3,*) 'i_dY6m=',idiag_dY6m
        write(3,*) 'i_dY7m=',idiag_dY7m
        write(3,*) 'i_dY8m=',idiag_dY8m
        write(3,*) 'i_dY9m=',idiag_dY9m
        write(3,*) 'i_h1m=',idiag_h1m
        write(3,*) 'i_h2m=',idiag_h2m
        write(3,*) 'i_h3m=',idiag_h3m
        write(3,*) 'i_h4m=',idiag_h4m
        write(3,*) 'i_h5m=',idiag_h5m
        write(3,*) 'i_h6m=',idiag_h6m
        write(3,*) 'i_h7m=',idiag_h7m
        write(3,*) 'i_h8m=',idiag_h8m
        write(3,*) 'i_h9m=',idiag_h9m
        write(3,*) 'i_cpfull=',idiag_cpfull
        write(3,*) 'i_cvfull=',idiag_cvfull
        write(3,*) 'i_cp1m=',idiag_cp1m
        write(3,*) 'i_cp2m=',idiag_cp2m
        write(3,*) 'i_cp3m=',idiag_cp3m
        write(3,*) 'i_cp4m=',idiag_cp4m
        write(3,*) 'i_cp5m=',idiag_cp5m
        write(3,*) 'i_cp6m=',idiag_cp6m
        write(3,*) 'i_cp7m=',idiag_cp7m
        write(3,*) 'i_cp8m=',idiag_cp8m
        write(3,*) 'i_cp9m=',idiag_cp9m
        write(3,*) 'i_e_intm=',idiag_e_intm
        write(3,*) 'i_Y1mz=',idiag_Y1mz
        write(3,*) 'i_Y2mz=',idiag_Y2mz
        write(3,*) 'i_Y3mz=',idiag_Y3mz
        write(3,*) 'i_Y4mz=',idiag_Y4mz
        write(3,*) 'i_Y5mz=',idiag_Y5mz
        write(3,*) 'i_Y6mz=',idiag_Y6mz
        write(3,*) 'i_Y7mz=',idiag_Y7mz
        write(3,*) 'i_Y8mz=',idiag_Y8mz
        write(3,*) 'i_Y9mz=',idiag_Y9mz
        write(3,*) 'ichemspec=indgen('//trim(schemspec)//') + '//trim(snd1)
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
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Chemical species mass fractions.
!
        case ('chemspec')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ichemspec(i1))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ichemspec(i1))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ichemspec(i1))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ichemspec(i1))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ichemspec(i1))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ichemspec(i1))
          slices%ready=.true.
        case ('chemspec2')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ichemspec(i2))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ichemspec(i2))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ichemspec(i2))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ichemspec(i2))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ichemspec(i2))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ichemspec(i2))
          slices%ready=.true.
        case ('chemspec3')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ichemspec(i3))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ichemspec(i3))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ichemspec(i3))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ichemspec(i3))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ichemspec(i3))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ichemspec(i3))
          slices%ready=.true.
        case ('chemspec4')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ichemspec(i4))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ichemspec(i4))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ichemspec(i4))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ichemspec(i4))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ichemspec(i4))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ichemspec(i4))
          slices%ready=.true.
        case ('chemspec5')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ichemspec(i5))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ichemspec(i5))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ichemspec(i5))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ichemspec(i5))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ichemspec(i5))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ichemspec(i5))
          slices%ready=.true.
        case ('chemspec6')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ichemspec(i6))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ichemspec(i6))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ichemspec(i6))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ichemspec(i6))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ichemspec(i6))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ichemspec(i6))
          slices%ready=.true.
        case ('chemspec7')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ichemspec(i7))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ichemspec(i7))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ichemspec(i7))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ichemspec(i7))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ichemspec(i7))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ichemspec(i7))
          slices%ready=.true.
        case ('chemspec8')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ichemspec(i8))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ichemspec(i8))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ichemspec(i8))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ichemspec(i8))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ichemspec(i8))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ichemspec(i8))
          slices%ready=.true.
        case ('chemspec9')
          slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,ichemspec(i9))
          slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,ichemspec(i9))
          slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,ichemspec(i9))
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,ichemspec(i9))
          if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2 ,m1:m2 ,iz3_loc,ichemspec(i9))
          if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2 ,m1:m2 ,iz4_loc,ichemspec(i9))
          slices%ready=.true.
      endselect
!
    endsubroutine get_slices_chemistry
!***********************************************************************
    subroutine find_species_index(species_name,ind_glob,ind_chem,found_specie)
!
!  Find index in the f array for specie
!
!  05-feb-08/nils: coded
!
      use Cdata
!
      integer, intent(out) :: ind_glob
      integer, intent(inout) :: ind_chem
      character (len=*), intent(in) :: species_name
      integer :: k
      logical, intent(out) :: found_specie
!
      ind_glob=0
    !  ind_chem=0
      do k=1,nchemspec
        if (trim(varname(ichemspec(k)))==species_name) then
          ind_glob=k+ichemspec(1)-1
          ind_chem=k
          exit
        endif
      enddo
!
!  Check if the species was really found
!
      if (ind_glob==0) then
        found_specie=.false.
      else
        found_specie=.true.
        if (lroot) print*,species_name,'   species index= ',ind_chem
      endif
!
    endsubroutine find_species_index
!***********************************************************************
    subroutine find_mass(element_name,MolMass)
!
!  Find mass of element
!
!  05-feb-08/nils: coded
!
      use Mpicomm
!
      character (len=*), intent(in) :: element_name
      real, intent(out) :: MolMass
!
      select case (element_name)
      case('H')  
        MolMass=1.00794
      case('C')  
        MolMass=12.0107
      case('N')  
        MolMass=14.00674
      case('O')  
        MolMass=15.9994
      case('Ar','AR') 
        MolMass=39.948
      case('He','HE') 
        MolMass=4.0026
      case default
        if (lroot) print*,'element_name=',element_name
        call stop_it('find_mass: Element not found!')
      end select
!
    endsubroutine find_mass
!***********************************************************************
    subroutine read_species(input_file)
!
!  This subroutine reads all species information from chem.inp
!  See the chemkin manual for more information on
!  the syntax of chem.inp.
!
!  06-mar-08/nils: coded
!
      use Mpicomm
!
      logical :: IsSpecie=.false., emptyfile
      integer :: k,file_id=123, StartInd, StopInd
      character (len=80) :: ChemInpLine
      character (len=*) :: input_file
!
      emptyFile=.true.
      k=1
      open(file_id,file=input_file)
      dataloop: do
        read(file_id,'(80A)',end=1000) ChemInpLine(1:80)
        emptyFile=.false.
!
!  Check if we are reading a line within the species section
!
        if (ChemInpLine(1:7)=="SPECIES")            IsSpecie=.true.
        if (ChemInpLine(1:3)=="END" .and. IsSpecie) IsSpecie=.false.
!
!  Read in species
!
        if (IsSpecie) then
          if (ChemInpLine(1:7) /= "SPECIES") then
            StartInd=1; StopInd =0
            stringloop: do
              StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
              if (StopInd==StartInd) then
                StartInd=StartInd+1
              else
                varname(ichemspec(k))=trim(ChemInpLine(StartInd:StopInd-1))
                StartInd=StopInd
                k=k+1
                if (k>nvar) then
                  print*,'nchemspec=',nchemspec
                  call stop_it("There were too many species, "//&
                      "please increase nchemspec!")
                endif
              endif
              if (StartInd==80) exit
            enddo stringloop
          endif
        endif
      enddo dataloop
!
!  Stop if chem.inp is empty
!
1000  if (emptyFile)  call stop_it('The input file chem.inp was empty!')
      close(file_id)
!
    endsubroutine read_species
!********************************************************************
    subroutine read_thermodyn(input_file)
!
!  This subroutine reads the thermodynamical data for all species 
!  from chem.inp. See the chemkin manual for more information on
!  the syntax of chem.inp.
!
!  06-mar-08/nils: coded
!
      character (len=*), intent(in) :: input_file
      integer :: file_id=123, ind_glob, ind_chem
      character (len=80) :: ChemInpLine
      integer :: In1,In2,In3,In4,In5,iElement,iTemperature,nn2,StopInd
      integer :: NumberOfElement_i
      logical :: IsThermo=.false., found_specie
      real, dimension(4) :: MolMass
      real, dimension(3) :: tmp_temp
      character (len=5) :: NumberOfElement_string,element_string
      character (len=10) :: specie_string,TemperatureNr_i
      real :: nne

      ind_chem=0
!
      open(file_id,file=input_file)
      dataloop2: do
        read(file_id,'(80A)',end=1001) ChemInpLine(1:80)
!
! Check if we are reading a line within the thermo section
!
        if (ChemInpLine(1:6)=="THERMO") IsThermo=.true.
        if (ChemInpLine(1:3)=="END" .and. IsThermo) IsThermo=.false.
!
! Read in thermo data
!
        if (IsThermo) then
          if (ChemInpLine(1:7) /= "THERMO") then
            StopInd=index(ChemInpLine,' ')
            specie_string=trim(ChemInpLine(1:StopInd-1))

            call find_species_index(specie_string,ind_glob,ind_chem,&
                found_specie)
!

!AB: the following causes an error if ind_chem=0
            if (ind_chem>0 .and. ind_chem<=nchemspec) then

            if (found_specie) then
!
! Find molar mass
!
              MolMass=0
              do iElement=1,4                
                In1=25+(iElement-1)*5
                In2=26+(iElement-1)*5
                In3=27+(iElement-1)*5
                In4=29+(iElement-1)*5
                if (ChemInpLine(In1:In1)==' ') then
                  MolMass(iElement)=0
                else
                  element_string=trim(ChemInpLine(In1:In2))
                  call find_mass(element_string,MolMass(iElement))
                  In5=verify(ChemInpLine(In3:In4),' ')+In3-1
                  NumberOfElement_string=trim(ChemInpLine(In5:In4))
                  read (unit=NumberOfElement_string,fmt='(I5)') &
                      NumberOfElement_i
                  MolMass(iElement)=MolMass(iElement)*NumberOfElement_i
                endif
              enddo
              species_constants(ind_chem,imass)=sum(MolMass)

!
! Find temperature-ranges for low and high temperature fitting
!
              do iTemperature=1,3
                In1=46+(iTemperature-1)*10
                In2=55+(iTemperature-1)*10
                if (iTemperature==3) In2=73
                In3=verify(ChemInpLine(In1:In2),' ')+In1-1
                TemperatureNr_i=trim(ChemInpLine(In3:In2))
                read (unit=TemperatureNr_i,fmt='(F10.1)') nne
                tmp_temp(iTemperature)=nne
              enddo
              species_constants(ind_chem,iTemp1)=tmp_temp(1)
              species_constants(ind_chem,iTemp2)=tmp_temp(3)
              species_constants(ind_chem,iTemp3)=tmp_temp(2)
!
            elseif (ChemInpLine(80:80)=="2") then
              ! Read iaa1(1):iaa1(5)
              read (unit=ChemInpLine(1:75),fmt='(5E15.8)')  &
                  species_constants(ind_chem,iaa1(1):iaa1(5))
!
            elseif (ChemInpLine(80:80)=="3") then
              ! Read iaa1(6):iaa5(3)
              read (unit=ChemInpLine(1:75),fmt='(5E15.8)')  &
                  species_constants(ind_chem,iaa1(6):iaa2(3))
            elseif (ChemInpLine(80:80)=="4") then
              ! Read iaa2(4):iaa2(7)
              read (unit=ChemInpLine(1:75),fmt='(4E15.8)')  &
                  species_constants(ind_chem,iaa2(4):iaa2(7))
            endif

          endif
          endif !(from ind_chem>0 query)
        endif
      enddo dataloop2
1001  continue
      close(file_id)
!

   endsubroutine read_thermodyn
!***********************************************************************
    subroutine read_reactions(input_file,NrOfReactions)
!
!  This subroutine reads all reaction information from chem.inp
!  See the chemkin manual for more information on
!  the syntax of chem.inp.
!
!  10-mar-08/nils: coded
!
      use Mpicomm
!
      logical :: IsReaction=.false.,LastSpecie,found_new_reaction=.false.
      logical, save :: find_specie, found_specie
      integer, optional :: NrOfReactions
      integer, save :: ind_glob, ind_chem
      integer :: i,k,file_id=123, StartInd, StopInd, StartInd_add
      integer :: StopInd_add, StopInd_add_,StopIndName
      integer :: VarNumber, VarNumber_add, SeparatorInd 
      integer :: StartSpecie,stoi, PlusInd
      integer :: LastLeftCharacter,ParanthesisInd,Mplussind
      character (len=80) :: ChemInpLine, ChemInpLine_add
      character (len=*) :: input_file
!
      if (present(NrOfReactions)) NrOfReactions=0
!
      k=1
      open(file_id,file=input_file)
      dataloop3: do
        if (found_new_reaction) then
          ChemInpLine=ChemInpLine_add
        else
          read(file_id,'(80A)',end=1012) ChemInpLine(1:80)
        endif
        found_new_reaction=.false.
!
!  Check if we are reading a line within the reactions section
!
        if (ChemInpLine(1:9)=="REACTIONS")            IsReaction=.true.
        if (ChemInpLine(1:3)=="END" .and. IsReaction) IsReaction=.false.
!
        if (present(NrOfReactions)) then
!
!  Find number of reactions
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
!  Read in species
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
! Find reactant side stoichiometric coefficients
!
                SeparatorInd=index(ChemInpLine(StartInd:),'<=')
                if (SeparatorInd==0) then
                  SeparatorInd=index(ChemInpLine(StartInd:),'=')
                endif
!
                ParanthesisInd=0
                MplussInd=0
!
                ParanthesisInd=index(ChemInpLine(StartInd:),'(+M)')
                MplussInd=index(ChemInpLine(StartInd:),'+M')
!
                found_new_reaction=.false.
                if (ParanthesisInd>0 .or. Mplussind>0) then
                  if (ParanthesisInd>0) then
                    LastLeftCharacter=min(ParanthesisInd,SeparatorInd)-1
                  else
                    LastLeftCharacter=min(MplussInd,SeparatorInd)-1
                  endif
!
!  reading of the additional data for (+M) case
!
100               read(file_id,'(80A)',end=1012) ChemInpLine_add(1:80)
                  if (ChemInpLine_add(1:1) == ' ') then

!
                    if (ParanthesisInd>0) then
                      Mplus_case(k)=.true.
                    endif
!
                    i=1
                    do while (i<80)
                      if (ChemInpLine_add(i:i)==' ') then
                        i=i+1
                      elseif (ChemInpLine_add(i:i+2)=='LOW') then
                        if (lroot) print*,ChemInpLine_add(i:i+2),&
                            '   coefficients for reaction ', &
                            reaction_name(k),'number ', k 
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
                        if (lroot) print*,ChemInpLine_add(i:i+3),&
                            '   coefficients for reaction ', reaction_name(k),&
                            'number ', k 
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
                        if (lroot) print*,ChemInpLine_add(i:i+3),&
                            '   coefficients for reaction ', reaction_name(k),&
                            'number ', k 
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
                        if (lroot) print*,' --------------  a_k4 coefficients-------------'
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
                                if (lroot) print*, 'a_k4(ind_chem,k)=',a_k4(ind_chem,k),&
                                    ind_chem,k
                              else
                                call stop_it("Did not find specie!")
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
!
                      endif
                    enddo
                    goto 100

!
                  else
                    found_new_reaction=.true.
                  endif

                else
                  LastLeftCharacter=SeparatorInd-1
                endif
!
                StartInd=1
                PlusInd=index(ChemInpLine(StartInd:LastLeftCharacter),'+')&
                    +StartInd-1
                do while (PlusInd<LastLeftCharacter .AND. PlusInd>0)
                  StopInd=PlusInd-1
                  call build_stoich_matrix(StartInd,StopInd,k,&
                      ChemInpLine,.false.)
                  StartInd=StopInd+2
                  PlusInd=index(ChemInpLine(StartInd:),'+')+StartInd-1
                enddo
                StopInd=LastLeftCharacter
                call build_stoich_matrix(StartInd,StopInd,k,ChemInpLine,.false.)
!
! Find product side stoichiometric coefficients
!
                StartInd=index(ChemInpLine,'>')+1
                if (StartInd==1) StartInd=index(ChemInpLine,'=')+1
                SeparatorInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
!
                ParanthesisInd=index(ChemInpLine(StartInd:),'(+M)')+StartInd-1
                MplussInd=index(ChemInpLine(StartInd:),'+M')+StartInd-1
!
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
!********************************************************************
    subroutine write_thermodyn()
!
!  This subroutine writes the thermodynamical data for every specie
!  to ./data/chem.out. 
!
!  06-mar-08/nils: coded
!
      use General
!
      character (len=20) :: input_file="./data/chem.out"
      character (len=5) :: ispec
      integer :: file_id=123,k
!
      open(file_id,file=input_file)
      write(file_id,*) 'Specie'
      write(file_id,*) 'MolMass Temp1 Temp2 Temp3'
      write(file_id,*) 'a1(1)  a1(2)  a1(3)  a1(4)  a1(5)  a1(6)  a1(7)'
      write(file_id,*) 'a2(1)  a2(2)  a2(3)  a2(4)  a2(5)  a2(6)  a2(7)'
      write(file_id,*) '***********************************************'
      dataloop2: do k=1,nchemspec
        write(file_id,*) varname(ichemspec(k))
        write(file_id,'(F10.2,3F10.2)') species_constants(k,imass),&
            species_constants(k,iTemp1:iTemp3)
        write(file_id,'(7E12.5)') species_constants(k,iaa1)
        write(file_id,'(7E12.5)') species_constants(k,iaa2)
      enddo dataloop2
!
      close(file_id)
!
      if (lroot) then
        print*,'Write pc_constants.pro in chemistiry.f90'
        open (143,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (143,*) 'specname=strarr(',nchemspec,')'
        write (143,*) 'specmass=fltarr(',nchemspec,')'
        do k=1,nchemspec
          call chn(k-1,ispec)
          write (143,*) 'specname[',trim(ispec),']=',"'",&
              trim(varname(ichemspec(k))),"'"
          write (143,*) 'specmass[',trim(ispec),']=',species_constants(k,imass)
        enddo
        close (143)
      endif
!
    endsubroutine write_thermodyn
!***************************************************************
    subroutine build_stoich_matrix(StartInd,StopInd,k,ChemInpLine,product)
!
!  DOCUMENT ME!
!
!  10-mar-08/nils: coded
!
      use Mpicomm
!
      integer, intent(in) :: StartInd,StopInd,k
      character (len=*), intent(in) :: ChemInpLine
      integer :: StartSpecie,ind_glob,ind_chem,stoi
      logical :: found_specie,product
!
      if (ChemInpLine(StartInd:StopInd) /= "M" ) then
        StartSpecie=verify(ChemInpLine(StartInd:StopInd),&
            "1234567890")+StartInd-1
        call find_species_index(ChemInpLine(StartSpecie:StopInd),&
            ind_glob,ind_chem,found_specie)
!
        if (.not. found_specie) call stop_it("Did not find specie!")
        if (StartSpecie==StartInd) then
          stoi=1
        else
          read (unit=ChemInpLine(StartInd:StartInd),fmt='(I1)') stoi
        endif
        if (product) then
          Sijm(ind_chem,k)=Sijm(ind_chem,k)+stoi
        else
          Sijp(ind_chem,k)=Sijp(ind_chem,k)+stoi
        endif
      endif
!
    endsubroutine build_stoich_matrix
!***************************************************************
    subroutine write_reactions()
!
!  DOCUMENT ME!!
!
!  11-mar-08/nils: coded
!
      use General, only: chn
!
      integer :: reac,spec
      character (len=80) :: reac_string,product_string,output_string
      character (len=5)  :: Sijp_string,Sijm_string
      character (len=1)  :: separatorp,separatorm
      character (len=20) :: input_file="./data/chem.out"
      integer :: file_id=123
!
      open(file_id,file=input_file,POSITION='APPEND',FORM='FORMATTED')
      write(file_id,*) 'REACTIONS'
      !open(file_id,file=input_file)
!
      do reac=1,mreactions
        reac_string=''
        product_string=''
        separatorp=''
        separatorm=''
!
        do spec=1,nchemspec
          if (Sijp(spec,reac)>0) then
            Sijp_string=''
            if (Sijp(spec,reac)>1) call chn(Sijp(spec,reac),Sijp_string)
            reac_string=trim(reac_string)//trim(separatorp)//&
                trim(Sijp_string)//trim(varname(ichemspec(spec)))
            separatorp='+'
          endif
          if (Sijm(spec,reac)>0) then
            Sijm_string=''
            if (Sijm(spec,reac)>1) call chn(Sijm(spec,reac),Sijm_string)
            product_string=trim(product_string)//trim(separatorm)//&
                trim(Sijm_string)//trim(varname(ichemspec(spec)))
            separatorm='+'
          endif
        enddo
!
        output_string=trim(reac_string)//'='//trim(product_string)
        write(unit=output_string(30:45),fmt='(E14.4)') B_n(reac)
        write(unit=output_string(47:62),fmt='(E14.4)') alpha_n(reac)
        write(unit=output_string(64:79),fmt='(E14.4)') E_an(reac)
        write(file_id,*) trim(output_string)
        if (maxval(abs(low_coeff(:,reac))) > 0.) then
          write(file_id,*) 'LOW/',low_coeff(:,reac)
        elseif (maxval(abs(high_coeff(:,reac))) > 0.) then
          write(file_id,*) 'HIGH/',high_coeff(:,reac)
        endif
        if (minval(a_k4(:,reac))<impossible) then
          write(file_id,*) a_k4(:,reac)
        endif
!
      enddo
!
      write(file_id,*) 'END'
      write(file_id,*) '(M+) case: ',Mplus_case
!
      close(file_id) 
!
    endsubroutine write_reactions
!***************************************************************
    subroutine get_reaction_rate(f,vreact_p,vreact_m,p)
!
!  This subroutine calculates forward and reverse reaction rates, 
!  if chem.inp file exists.
!  For more details see Chemkin Theory Manual
!
!  17-mar-08/natalia: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f 
      real, dimension (nx,nreactions), intent(out) :: vreact_p, vreact_m
      real, dimension (nx):: mix_conc
!
      type (pencil_case) :: p
      real, dimension (nx) :: dSR=0.,dHRT=0.,Kp,Kc,prod1,prod2
      real, dimension (nx) :: kf=0., kr=0.
      real, dimension (nx) :: rho_cgs,p_atm
      real :: Rcal
      integer :: k , reac, j, i, v, t
      real  :: sum_tmp=0., T_low, T_mid, T_up, tmp
      logical,save :: lwrite=.true.
      character (len=20) :: input_file="./data/react.out"
      integer :: file_id=123
      real :: B_n_0,alpha_n_0,E_an_0
      real, dimension (nx) ::  kf_0,Kc_0,Pr,sum_sp,prod1_0,prod2_0
      integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9

!
      if (lwrite)  open(file_id,file=input_file)
!
!  p is in atm units; atm/bar=1./10.13
!
      Rcal=Rgas_unit_sys/4.14*1e-7
      rho_cgs=p%rho*unit_mass/unit_length**3
      p_atm=1e6*unit_length**3/unit_energy
      if (lwrite)  write(file_id,*)'T= ',   p%TT
      if (lwrite)  write(file_id,*)'p_atm= ',   p_atm
!
!  Dimensionless Standard-state molar entropy  S0/R
!
      if (lwrite)  write(file_id,*)'**************************'
      if (lwrite) write(file_id,*)'S0_R'
      if (lwrite)  write(file_id,*)'**************************'
!
      do k=1,nchemspec
        T_low=species_constants(k,iTemp1)
        T_mid=species_constants(k,iTemp2)
        T_up= species_constants(k,iTemp3)
        do i=1,nx
          if (p%TT(i) <= T_mid) then
            tmp=0. 
            do j=2,5
              tmp=tmp+species_constants(k,iaa2(j))*p%TT(i)**(j-1)/(j-1) 
            enddo
            S0_R(i,k)=species_constants(k,iaa2(1))*p%lnTT(i)+tmp&
                +species_constants(k,iaa2(7))
          else
            tmp=0. 
            do j=2,5 
              tmp=tmp+species_constants(k,iaa1(j))*p%TT(i)**(j-1)/(j-1) 
            enddo
            S0_R(i,k)=species_constants(k,iaa1(1))*p%lnTT(i)+tmp&
                +species_constants(k,iaa1(7))
          endif
        enddo
!
        if (lwrite)  then
          write(file_id,*)&
              varname(ichemspec(k)),&
              maxval(S0_R(:,k)), &
              minval(S0_R(:,k))
        endif
      enddo
!
!  calculation of the reaction rate
!
      if (lwrite) write(file_id,*)'**************************'
      if (lwrite) write(file_id,*)'Reaction rates'
      if (lwrite) write(file_id,*)'**************************'
!
      do reac=1,nreactions
!
!  Find forward rate constant for reaction 'reac'
!
        kf(:)=B_n(reac)*p%TT(:)**alpha_n(reac)*exp(-E_an(reac)/Rcal/p%TT(:))
        if (lwrite)  write(file_id,*) 'Nreact= ',reac,  'kf=', maxval(kf)
!
!  Find backward rate constant for reaction 'reac'
!
        dSR=0.
        dHRT=0.
        sum_tmp=0.
        do k=1,nchemspec
          dSR(:) =dSR(:)+(Sijm(k,reac) -Sijp(k,reac))*S0_R(:,k)
          dHRT(:)=dHRT(:)+(Sijm(k,reac)-Sijp(k,reac))*H0_RT(l1:l2,m,n,k)  
          sum_tmp=sum_tmp+(Sijm(k,reac)-Sijp(k,reac))
        enddo
        if (lwrite) write(file_id,*) 'Nreact= ',reac,'dSR= ', maxval(dSR)
        if (lwrite) write(file_id,*) 'Nreact= ',reac,'dHRT= ', maxval(dHRT)
        Kp=exp(dSR-dHRT)
        if (sum_tmp==0.) then
          Kc=Kp
        else
          Kc=Kp*(p_atm/(p%TT*Rgas))**sum_tmp
        endif
        kr(:)=kf(:)/Kc
        if (lwrite) write(file_id,*) 'Nreact= ',reac,'Kc= ', maxval(Kc)
        if (lwrite) write(file_id,*) 'Nreact= ',reac,  'kr=', maxval(kr)
        if (lwrite) write(file_id,*)'**************************'
!
!  Find the product of the species molar consentrations (where
!  each molar consentration is taken to the power of the number)
!
        prod1=1.
        prod2=1.
        do k=1,nchemspec
          prod1=prod1*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)&
              /species_constants(k,imass))**Sijp(k,reac)
        enddo
        do k=1,nchemspec
          prod2=prod2*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)&
              /species_constants(k,imass))**Sijm(k,reac)
        enddo
!
!  Finalize writing to file
!
        if (lwrite) write(file_id,*) ''
        if (lwrite) write(file_id,*) '*******************'
        if (lwrite) print*,'get_reaction_rate: writing react.out file'
        if (lwrite) close(file_id)
        lwrite=.false.
!
!  Multiply by third body reaction term
!
        if (minval(a_k4(:,reac))<impossible) then
          sum_sp=0.
          do k=1,nchemspec  
            sum_sp=sum_sp+a_k4(k,reac)*f(l1:l2,m,n,ichemspec(k))  &
                *rho_cgs(:)/species_constants(k,imass)
          enddo
          mix_conc=sum_sp
        else
          sum_sp=1.
          mix_conc=rho_cgs(:)*p%mu1(:)/unit_mass
        endif
!
!  The Lindstrom approach to the fall of reactions
!
        Kc_0=Kc
        if (maxval(abs(low_coeff(:,reac))) > 0.) then
          B_n_0=low_coeff(1,reac) 
          alpha_n_0=low_coeff(2,reac)
          E_an_0=low_coeff(3,reac)                      
          kf_0(:)=B_n_0*p%TT(:)**alpha_n_0*exp(-E_an_0/Rcal/p%TT(:))
          Pr=kf_0/kf*mix_conc
          kf=kf*(Pr/(1.+Pr))
          kr(:)=kf(:)/Kc_0
        elseif (maxval(abs(high_coeff(:,reac))) > 0.) then
          B_n_0=high_coeff(1,reac) 
          alpha_n_0=high_coeff(2,reac)
          E_an_0=high_coeff(3,reac)                      
          kf_0(:)=B_n_0*p%TT(:)**alpha_n_0*exp(-E_an_0/Rcal/p%TT(:))
          Pr=kf_0/kf*mix_conc
          kf=kf*(1./(1.+Pr))
          kr(:)=kf(:)/Kc_0
        endif
!
!  Find forward (vreact_p) and backward (vreact_m) rate of 
!  progress variable. 
!  (vreact_p - vreact_m) is labeled q in the chemkin manual
! 
        if (Mplus_case(reac)) then
          vreact_p(:,reac)=prod1*kf
          vreact_m(:,reac)=prod2*kr
        else
          vreact_p(:,reac)=prod1*kf*sum_sp
          vreact_m(:,reac)=prod2*kr*sum_sp
        endif
!


 ! This part calculates forward and reverse reaction rates
!  for the test case R->P 
!
!  For more details see Doom, et al., J. Comp. Phys., 226, 2007

      if (l1step_test) then
       
        do i=1,nx
         if (p%TT(i) > Tc) then
         vreact_p(i,reac)=f(l1,m,n,iux)**2*p%rho(1)*Cp_const/lambda_const*beta*(beta-1.) &
                          *(1.-f(l1-1+i,m,n,ichemspec(ipr)))
        !  vreact_p(i,reac)=f(l1,m,n,iux)**2*Cp_const/lambda_const*beta*(beta-1.) &
         !                 *(1.-p%TT(i)/Tinf)

        else
         vreact_p(i,reac)=0.
        endif
       enddo 
         vreact_m(:,reac)=0.
      endif
      enddo
!
    endsubroutine get_reaction_rate
!***************************************************************
    subroutine calc_reaction_term(f,p)
!
!  DOCUMENT ME!!
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx,mreactions) :: vreactions,vreactions_p,vreactions_m
      real, dimension (nx) :: xdot
      type (pencil_case) :: p
      integer :: k,j
      real :: sum_omega
      real :: sum_Y
      integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9
!
      p%DYDt_reac=0.
!
!  if we do reactions, we must calculate the reaction speed vector
!  outside the loop where we multiply it by the stoichiometric matrix
!
      if (.not. lcheminp) then
!
!  Axel' case
!
        do j=1,nreactions
          vreactions_p(:,j)=kreactions_p(j)*kreactions_z(n,j)
          vreactions_m(:,j)=kreactions_m(j)*kreactions_z(n,j)
          do k=1,nchemspec
            vreactions_p(:,j)=vreactions_p(:,j)*&
                f(l1:l2,m,n,ichemspec(k))**Sijm(k,j)
            vreactions_m(:,j)=vreactions_m(:,j)*&
                f(l1:l2,m,n,ichemspec(k))**Sijp(k,j)
          enddo
        enddo
      else
!        
!  Chemkin data case
!

         call get_reaction_rate(f,vreactions_p,vreactions_m,p)
      endif
!
!  Calculate rate of progress variable (labeled q in the chemkin manual)
!
      vreactions=vreactions_p-vreactions_m
!
!  Calculate production rate for all species k (called \dot(\omega)_k 
!  in the chemkin manual)
!
      do k=1,nchemspec  
        xdot=0.
        do j=1,nreactions
          xdot=xdot+stoichio(k,j)*vreactions(:,j)
        enddo
        if (lcheminp) then
         if (l1step_test) then
          xdot=-xdot!/p%rho
         else
          xdot=-xdot*species_constants(k,imass)/p%rho
         endif
        endif
        p%DYDt_reac(:,k)=xdot*unit_time

      enddo
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
!  Calculate diagnostic quantities
!
      if (ldiagnos) then
        if (idiag_dY1m/=0) call sum_mn_name(p%DYDt_reac(:,i1),idiag_dY1m)
        if (idiag_dY2m/=0) call sum_mn_name(p%DYDt_reac(:,i2),idiag_dY2m)
        if (idiag_dY3m/=0) call sum_mn_name(p%DYDt_reac(:,i3),idiag_dY3m)
        if (idiag_dY4m/=0) call sum_mn_name(p%DYDt_reac(:,i4),idiag_dY4m)
        if (idiag_dY5m/=0) call sum_mn_name(p%DYDt_reac(:,i5),idiag_dY5m)
        if (idiag_dY6m/=0) call sum_mn_name(p%DYDt_reac(:,i6),idiag_dY6m)
        if (idiag_dY7m/=0) call sum_mn_name(p%DYDt_reac(:,i7),idiag_dY7m)
        if (idiag_dY8m/=0) call sum_mn_name(p%DYDt_reac(:,i8),idiag_dY8m)
        if (idiag_dY9m/=0) call sum_mn_name(p%DYDt_reac(:,i9),idiag_dY9m)
        if (idiag_h1m/=0) call sum_mn_name(H0_RT(l1:l2,m,n,i1)*Rgas*&
            p%TT(:)/species_constants(i1,imass),idiag_h1m)
        if (idiag_h2m/=0) call sum_mn_name(H0_RT(l1:l2,m,n,i2)*Rgas*&
            p%TT(:)/species_constants(i2,imass),idiag_h2m)
        if (idiag_h3m/=0) call sum_mn_name(H0_RT(l1:l2,m,n,i3)*Rgas*&
            p%TT(:)/species_constants(i3,imass),idiag_h3m)
        if (idiag_h4m/=0) call sum_mn_name(H0_RT(l1:l2,m,n,i4)*Rgas*&
            p%TT(:)/species_constants(i4,imass),idiag_h4m)
        if (idiag_h5m/=0) call sum_mn_name(H0_RT(l1:l2,m,n,i5)*Rgas*&
            p%TT(:)/species_constants(i5,imass),idiag_h5m)
        if (idiag_h6m/=0) call sum_mn_name(H0_RT(l1:l2,m,n,i6)*Rgas*&
            p%TT(:)/species_constants(i6,imass),idiag_h6m)
        if (idiag_h7m/=0) call sum_mn_name(H0_RT(l1:l2,m,n,i7)*Rgas*&
            p%TT(:)/species_constants(i7,imass),idiag_h7m)
        if (idiag_h8m/=0) call sum_mn_name(H0_RT(l1:l2,m,n,i8)*Rgas*&
            p%TT(:)/species_constants(i8,imass),idiag_h8m)
        if (idiag_h9m/=0) call sum_mn_name(H0_RT(l1:l2,m,n,i9)*Rgas*&
            p%TT(:)/species_constants(i9,imass),idiag_h9m)
      endif
!
    endsubroutine calc_reaction_term
!***************************************************************
    subroutine read_transport_data  
! 
!  Reading of the chemkin transport data
!
!  01-apr-08/natalia: coded
!
      use Mpicomm
!
      logical :: IsSpecie=.false., emptyfile
      logical ::  found_specie
      integer :: file_id=123, ind_glob, ind_chem
      character (len=80) :: ChemInpLine
      character (len=10) :: specie_string
      integer :: VarNumber
      integer :: StartInd,StopInd,StartInd_1,StopInd_1
!
      emptyFile=.true.
!
      StartInd_1=1; StopInd_1 =0
      open(file_id,file="tran.dat")
!
      if (lroot) print*, 'the following species are found '//&
          'in tran.dat: beginning of the list:'
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

1000  if (emptyFile)  call stop_it('The input file tran.dat was empty!')
!
      if (lroot) print*, 'the following species are found in tran.dat: end of the list:'
!
      close(file_id)
!
    endsubroutine read_transport_data
!***************************************************************
    subroutine  calc_collision_integral(omega,lnTst,Omega_kl)
!
!  Get coefficients for calculating of the collision integral 
!  This routind is called from calc_diff_visc_coeff, which again is called from
!  calc_for_chem_mixture, which is why we work on full chunks of arrays here.  
!  
!  03-apr-08/natalia: coded
!
      use Mpicomm
!
      character (len=*), intent(in) :: omega
      real,  dimension(mx,my,mz), intent(in)  :: lnTst
      real,  dimension(mx,my,mz), intent(out) :: Omega_kl
      integer :: i 
      real, dimension(8) :: aa
!
      select case (omega)
      case('Omega11')
        aa(1)= 6.96945701E-1
        aa(2)= 3.39628861E-1
        aa(3)= 1.32575555E-2
        aa(4)=-3.41509659E-2
        aa(5)= 7.71359429E-3
        aa(6)= 6.16106168E-4
        aa(7)=-3.27101257E-4
        aa(8)= 2.51567029E-5
      case('Omega22')  
        aa(1)= 6.33225679E-1
        aa(2)= 3.14473541E-1
        aa(3)= 1.78229325E-2
        aa(4)=-3.99489493E-2
        aa(5)= 8.98483088E-3
        aa(6)= 7.00167217E-4
        aa(7)=-3.82733808E-4
        aa(8)= 2.97208112E-5
      case default
        call stop_it('Insert Omega_kl')
      end select
!
      Omega_kl=0.
      do i=1,8
        Omega_kl=Omega_kl+aa(i)*(lnTst)**(i-1)
      enddo
      Omega_kl=1./Omega_kl
!
    endsubroutine  calc_collision_integral
!***************************************************************
    subroutine calc_diff_visc_coef(f)
!
!  Calculation of Bin_Diff_coef
!  This routind is called from calc_for_chem_mixture, 
!  which is why we work on full chunks of arrays here.  
!
!  WHO, WHEN?
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
      real, dimension (mx,my,mz) :: Omega_kl, prefactor
      real, dimension (mx,my,mz) :: lnTjk,lnTk_array
      integer :: k,j,j1,j2,j3
      real :: eps_jk, sigma_jk, m_jk, delta_jk, delta_st
      character (len=7) :: omega
      real :: Na=6.022E23,tmp_local,tmp_local2, lnTk

      tmp_local=3./16.*(2.*k_B_cgs**3/pi)**0.5

      do j3=nn1,nn2
      do j2=mm1,mm2
      do j1=ll1,ll2
       prefactor(j1,j2,j3)=tmp_local*(TT_full(j1,j2,j3))**0.5*unit_length**3&
          /(Rgas_unit_sys*rho_full(j1,j2,j3)*mu1_full(j1,j2,j3))
      enddo
      enddo 
      enddo

      omega="Omega11"

    if (.not. lfix_Sc) then

    if (BinDif_simple) then
        do k=1,nchemspec
          do j=k,nchemspec
            eps_jk=(tran_data(j,2)*tran_data(k,2))**0.5
            sigma_jk=0.5*(tran_data(j,3)+tran_data(k,3))*1e-8
            m_jk=(species_constants(j,imass)*species_constants(k,imass)) &
                /(species_constants(j,imass)+species_constants(k,imass))/Na
            tmp_local=(m_jk)**(-0.5)*(sigma_jk*unit_length)**(-2)*unit_time

          do j3=nn1,nn2
          do j2=mm1,mm2
          do j1=ll1,ll2
            lnTjk(j1,j2,j3)=f(j1,j2,j3,ilnTT)-log(eps_jk)
          !  Omega_kl(j1,j2,j3)=(6.96945701E-1 +3.39628861E-1*lnTjk(j1,j2,j3))


           Omega_kl(j1,j2,j3)= &
                 (6.96945701E-1   +3.39628861E-1*lnTjk(j1,j2,j3) &
                  +1.32575555E-2*lnTjk(j1,j2,j3)*lnTjk(j1,j2,j3) &
                  -3.41509659E-2*lnTjk(j1,j2,j3)**3 &
                  +7.71359429E-3*lnTjk(j1,j2,j3)**4 &
                  +6.16106168E-4*lnTjk(j1,j2,j3)**5 &
                  -3.27101257E-4*lnTjk(j1,j2,j3)**6 &
                  +2.51567029E-5*lnTjk(j1,j2,j3)**7)


            Bin_Diff_coef(j1,j2,j3,k,j)=prefactor(j1,j2,j3)  &
                                       *Omega_kl(j1,j2,j3)*tmp_local
          enddo
          enddo
          enddo

          enddo
        enddo

       do j3=nn1,nn2
       do j2=mm1,mm2
       do j1=ll1,ll2
        do k=1,nchemspec
        do j=1,k-1
            Bin_Diff_coef(j1,j2,j3,k,j)=Bin_Diff_coef(j1,j2,j3,j,k)
        enddo
        enddo
       enddo
       enddo
       enddo

      else
         do k=1,nchemspec
          do j=1,nchemspec
!
            eps_jk=(tran_data(j,2)*tran_data(k,2))**0.5
            sigma_jk=0.5*(tran_data(j,3)+tran_data(k,3))*1e-8
            delta_jk=0.5*tran_data(j,4)*tran_data(k,4)
            m_jk=(species_constants(j,imass)*species_constants(k,imass)) &
                /(species_constants(j,imass)+species_constants(k,imass))/Na
!
     !       do j3=nn1,nn2
     !       do j2=mm1,mm2
     !       do j1=ll1,ll2
     !        lnTjk(j1,j2,j3)=f(j1,j2,j3,ilnTT)-log(eps_jk)
     !       enddo
     !      enddo
     !       enddo
!
     !       call calc_collision_integral(omega,lnTjk,Omega_kl)


!
            do j3=nn1,nn2
            do j2=mm1,mm2
            do j1=ll1,ll2

             lnTjk(j1,j2,j3)=f(j1,j2,j3,ilnTT)-log(eps_jk)

            Omega_kl(j1,j2,j3)= &
                 1./(6.96945701E-1   +3.39628861E-1*lnTjk(j1,j2,j3) &
                  +1.32575555E-2*lnTjk(j1,j2,j3)*lnTjk(j1,j2,j3) &
                  -3.41509659E-2*lnTjk(j1,j2,j3)**3 &
                  +7.71359429E-3*lnTjk(j1,j2,j3)**4 &
                  +6.16106168E-4*lnTjk(j1,j2,j3)**5 &
                  -3.27101257E-4*lnTjk(j1,j2,j3)**6 &
                  +2.51567029E-5*lnTjk(j1,j2,j3)**7)


             Bin_Diff_coef(j1,j2,j3,k,j)=prefactor(j1,j2,j3)/sqrt(m_jk)/sigma_jk**2 &
                /(Omega_kl(j1,j2,j3)+0.19*delta_jk/(TT_full(j1,j2,j3)/eps_jk)) &
                /(unit_length**2/unit_time)
            enddo
            enddo
            enddo
! 
          enddo
        enddo

      endif
     endif
      omega="Omega22"

      tmp_local=5./16.*(k_B_cgs/Na/pi)**0.5/(1e-8)**2

      do k=1,nchemspec
!
    !   lnTk=f(:,:,:,ilnTT)-log(tran_data(k,2))

      tmp_local2=(species_constants(k,imass))**0.5/(tran_data(k,3))**2 &
                   *tmp_local

    !    tmp_local2=(species_constants(k,imass))**0.5 &
    !               /(tran_data(k,3))**2*tmp_local
        if (visc_simple) then
          do j3=nn1,nn2
          do j2=mm1,mm2
          do j1=ll1,ll2
           lnTk=f(j1,j2,j3,ilnTT)-log(tran_data(k,2))

           Omega_kl(j1,j2,j3)=(6.33225679E-1 +3.14473541E-1*lnTk+1.78229325E-2*lnTk*lnTk &
                   -3.99489493E-2*lnTk*lnTk*lnTk+8.98483088E-3*lnTk**4 &
                  +7.00167217E-4*lnTk**5-3.82733808E-4*lnTk**6+2.97208112E-5*lnTk**7)

           species_viscosity(j1,j2,j3,k)=TT_full(j1,j2,j3)**0.5 &
                           *(Omega_kl(j1,j2,j3))*tmp_local2 &
                           /(unit_mass/unit_length/unit_time)
          enddo
          enddo
          enddo
        else
         delta_st=tran_data(k,4)**2/2./tran_data(k,2)/&
          (tran_data(k,3))**3*(1e-8**3)
!
         lnTk_array=f(:,:,:,ilnTT)-log(tran_data(k,2))
         call calc_collision_integral(omega,lnTk_array,Omega_kl)
          do j3=nn1,nn2
          do j2=mm1,mm2
          do j1=ll1,ll2
           species_viscosity(j1,j2,j3,k)=TT_full(j1,j2,j3)**0.5/(Omega_kl(j1,j2,j3) &
                 +0.2*delta_st/(TT_full(j1,j2,j3)/tran_data(k,2)))*tmp_local2 &
                 /(unit_mass/unit_length/unit_time)
          enddo
          enddo
          enddo
        endif
      enddo
!
    endsubroutine calc_diff_visc_coef
!***************************************************************
    subroutine calc_diffusion_term(f,p)
!
!  DOCUMENT ME!!!
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: gXX, gDiff_full_add, gchemspec
      real, dimension (nx) :: del2chemspec
      real, dimension (nx) :: diff_op,diff_op1,diff_op2,xdot, del2XX 
      real :: diff_k
      integer :: j,k
!
      intent(in) :: f
!      
      p%DYDt_diff=0.
      diff_k=chem_diff
!
      do k=1,nchemspec
!
        if (.not. lcheminp) then  
          if (chem_diff/=0.) then
            diff_k=chem_diff*chem_diff_prefactor(k)
            if (headtt) print*,'dchemistry_dt: k,diff_k=',k,diff_k
            call del2(f,ichemspec(k),del2chemspec) 
            call grad(f,ichemspec(k),gchemspec) 
            call dot_mn(p%glnrho,gchemspec,diff_op)
            diff_op=diff_op+del2chemspec
            p%DYDt_diff(:,k)=diff_k*diff_op
          endif
        else 
!
          if (ldiffusion) then
!
            call del2(XX_full(:,:,:,k),del2XX)
            call grad(XX_full(:,:,:,k),gXX)
            call dot_mn(p%glnrho,gXX,diff_op1)
!
            call grad(Diff_full_add(:,:,:,k),gDiff_full_add)
            call dot_mn(gDiff_full_add,gXX,diff_op2)
          endif
          p%DYDt_diff(:,k)=Diff_full_add(l1:l2,m,n,k)*(del2XX+diff_op1)+diff_op2
        endif
      enddo
!


    endsubroutine calc_diffusion_term
!***************************************************************
    subroutine calc_heatcond_chemistry(f,df,p)
!
!  Calculate gamma*chi*(del2lnT+gradlnTT.grad(lnT+lnrho+lncp+lnchi))
!
!  29-feb-08/natalia: coded
!
      use Sub
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension(nx) :: g2TT, g2TTlnlambda=0., tmp1
      integer :: i
!

      call dot(p%glnTT,p%glnlambda,g2TTlnlambda)
      call dot(p%glnTT,p%glnTT,g2TT)
!
!  Add heat conduction to RHS of temperature equation
!

      if (l1step_test) then
       tmp1= p%lambda(:)*(p%del2lnTT+g2TT)*p%cv1/p%rho(:)
      else
       tmp1= (p%lambda(:)*(p%del2lnTT+g2TT)+g2TTlnlambda)*p%cv1/p%rho(:)
      endif


      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + tmp1 
!
      RHS_T_full(l1:l2,m,n)=RHS_T_full(l1:l2,m,n) + tmp1 
!

    endsubroutine calc_heatcond_chemistry
!***********************************************************************
    subroutine calc_cs2x(cs2x,topbot,f)
!
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: cs2x
      character (len=3) :: topbot
      integer :: k
!
      select case(topbot)
      case('bot')               ! bottom boundary
        cs2x=cp_full(l1,:,:)/cv_full(l1,:,:)*&
            mu1_full(l1,:,:)*TT_full(l1,:,:)*Rgas
      case('top')               ! top boundary
        cs2x=cp_full(l2,:,:)/cv_full(l2,:,:)*&
            mu1_full(l2,:,:)*TT_full(l2,:,:)*Rgas
      case default
        print*, "calc_cs2x: ", topbot, " should be `top' or `bot'"
      endselect
!
    endsubroutine calc_cs2x
!***********************************************************************
    subroutine get_cs2_full(cs2_full)
!
      use Mpicomm
!
      real, dimension (mx,my,mz) :: cs2_full
!
      intent(out) :: cs2_full

      integer :: j,k

!
     do j=1,my
     do k=1,mz
      if (minval(cv_full(:,j,k))<=0) then
       cs2_full(:,j,k)=0.
      else
       cs2_full(:,j,k)=cp_full(:,j,k)/cv_full(:,j,k)*mu1_full(:,j,k)*TT_full(:,j,k)*Rgas
      endif
     enddo
     enddo
!

    endsubroutine get_cs2_full
!*************************************************************
    subroutine get_gamma_full(gamma_full)
!
      use Mpicomm
!
      real, dimension (mx,my,mz) :: gamma_full
!
      intent(out) :: gamma_full

      integer :: j,k
!      
      do j=1,my
      do k=1,mz
       if (minval(cv_full(:,j,k))<=0) then
        gamma_full(:,j,k)=0.
       else
        gamma_full(:,j,k)=cp_full(:,j,k)/cv_full(:,j,k)
       endif
      enddo
      enddo


    endsubroutine get_gamma_full
!*************************************************************
    subroutine get_p_infx(p_infx,topbot)
!
!  DOCUMENT ME!!
!
      use Mpicomm
!
      real, dimension (ny,nz) :: p_infx
      character (len=3) :: topbot
      logical, save :: read_P=.true.
      real, dimension (2,ny,nz), save :: pp_infx
      real :: PP
!
      logical :: emptyfile=.true., lairinp=.false.
      integer :: file_id=123
      character (len=80) :: ChemInpLine
      character (len=10) :: specie_string
      character (len=1)  :: tmp_string 
      integer :: VarNumber,i,j,k=1
      character (len=20) :: input_file='air.dat'
!
      integer :: StartInd,StopInd,StartInd_1,StopInd_1
      integer :: iostat
!
      StartInd_1=1; StopInd_1 =0
!
      inquire(FILE=input_file, EXIST=lairinp)
!
      if (read_P .and. lairinp) then
!
        StartInd_1=1; StopInd_1 =0
        open(file_id,file="air.dat")
!
        if (lroot) print*, 'the following parameters '//&
            'and species are found in air.dat (volume fraction fraction in %): '
!        
        dataloop: do
!
          read(file_id,'(80A)',IOSTAT=iostat) ChemInpLine(1:80)
!
          if (iostat < 0) exit dataloop
          emptyFile=.false.
          StartInd_1=1; StopInd_1=0
          StopInd_1=index(ChemInpLine,' ') 
          specie_string=trim(ChemInpLine(1:StopInd_1-1)) 
          tmp_string=trim(ChemInpLine(1:1)) 
!
          if (tmp_string == '!' .or. tmp_string == ' ') then
          elseif (tmp_string == 'P') then
!
            StartInd=1; StopInd =0
!
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
            StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
!
            read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)'), PP
            if (lroot) print*, ' Pressure, Pa   ', PP
!
          endif
        enddo dataloop
!
        close(file_id)
!
        if (emptyFile) then
          pp_infx(:,:,:)=0.
          read_P=.false.
          call stop_it('air.dat file is empty')
        else
          pp_infx(:,:,:)=PP*10.
          read_P=.false.
        endif
      endif
!
      select case(topbot)
      case('bot')               ! bottom boundary
        p_infx(1:ny,1:nz)=pp_infx(1,1:ny,1:nz)
      case('top')               ! top boundary
        p_infx(1:ny,1:nz)=pp_infx(2,1:ny,1:nz)
      case default
        print*, "get_p_infx: ", topbot, " should be `top' or `bot'"
      endselect
!
    endsubroutine get_p_infx
!*************************************************************
    subroutine get_p_infy(p_infy,topbot)
!
!  DOCUMENT ME!!
!
      use Mpicomm
      real, dimension (nx,nz) :: p_infy
      character (len=3) :: topbot
      logical, save :: read_P=.true.
      real, dimension (2,nx,nz), save :: pp_infy
      real :: PP
!
      logical :: emptyfile
      integer :: file_id=123
      character (len=80) :: ChemInpLine
      character (len=10) :: specie_string
      character (len=1)  :: tmp_string 
      integer :: VarNumber,i,j,k=1
!
      integer :: StartInd,StopInd,StartInd_1,StopInd_1
      integer :: iostat
!
      StartInd_1=1; StopInd_1 =0
!
      if (read_P) then
!
        StartInd_1=1; StopInd_1 =0
        open(file_id,file="air.dat")

        if (lroot) print*, 'the following parameters and '//&
            'species are found in air.dat (volume fraction fraction in %): '
!
        dataloop: do
!
          read(file_id,'(80A)',IOSTAT=iostat) ChemInpLine(1:80)
          if (iostat < 0) exit dataloop
          emptyFile=.false.
          StartInd_1=1; StopInd_1=0
          StopInd_1=index(ChemInpLine,' ') 
          specie_string=trim(ChemInpLine(1:StopInd_1-1)) 
          tmp_string=trim(ChemInpLine(1:1)) 
!
          if (tmp_string == '!' .or. tmp_string == ' ') then
          elseif (tmp_string == 'P') then
!
            StartInd=1; StopInd =0
!
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
            StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
!
            read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)'), PP
            if (lroot) print*, ' Pressure, Pa   ', PP
!
          endif
        enddo dataloop
        close(file_id)
!
        pp_infy(:,:,:)=PP*10.
        read_P=.false.
      endif
!
      select case(topbot)
      case('bot')               ! bottom boundary
        p_infy(:,:)=pp_infy(1,:,:)
      case('top')               ! top boundary
        p_infy(:,:)=pp_infy(2,:,:)
      case default
        print*, "get_p_infy: ", topbot, " should be `top' or `bot'"
      endselect
!
    endsubroutine get_p_infy
!*************************************************************
    subroutine get_rhs_Y(topbot,j,bound_reac_term)
!
!  DOCUMENT ME!!
!
      use Mpicomm
!
      real, dimension (ny,nz,nchemspec) :: bound_reac_term
      integer :: j
      character (len=3) :: topbot
!
      intent(out) :: bound_reac_term
      intent(in) :: j
!
      if (j==1) then
!
        select case(topbot)
        case('bot')               ! bottom boundary
          bound_reac_term(:,:,1:nchemspec)=&
              RHS_Y_full(l1,m1:m2,n1:n2,1:nchemspec)
        case('top')               ! top boundary
          bound_reac_term(:,:,1:nchemspec)=&
              RHS_Y_full(l2,m1:m2,n1:n2,1:nchemspec)
!
        case default
          print*, "get_rhs_Y: ", topbot, " should be `top' or `bot'"
        endselect
!
      elseif (j==2) then
!
      endif
!
    endsubroutine get_rhs_Y
!*************************************************************
    subroutine get_rhs_T(topbot,j,bound_rhs_T)
!
!  DOCUMENT ME!!
!
      use Mpicomm
!
      real, dimension (ny,nz) :: bound_rhs_T
      integer :: j
      character (len=3) :: topbot
!
      intent(out) :: bound_rhs_T
      intent(in) :: j
!
      if (j==1) then
!
        select case(topbot)
        case('bot')               ! bottom boundary
          bound_rhs_T(:,:)=RHS_T_full(l1,m1:m2,n1:n2)
        case('top')               ! top boundary
          bound_rhs_T(:,:)=RHS_T_full(l2,m1:m2,n1:n2)
!
        case default
          print*, "get_rhs_T: ", topbot, " should be `top' or `bot'"
        endselect
!
      elseif (j==2) then
!
      endif
!
    endsubroutine get_rhs_T
!*************************************************************
    subroutine air_field(f)
!
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: sum_Y
!
      logical :: IsSpecie=.false., emptyfile=.true.
      logical ::  found_specie
      integer :: file_id=123, ind_glob, ind_chem
      character (len=80) :: ChemInpLine
      character (len=10) :: specie_string
      character (len=1)  :: tmp_string 
      integer :: VarNumber,i,j,k=1
      real :: YY_k, air_mass, TT=300., PP=10.13e4
      real, dimension(nchemspec)    :: stor2 
      integer, dimension(nchemspec) :: stor1
!
      integer :: StartInd,StopInd,StartInd_1,StopInd_1
      integer :: iostat
!
      air_mass=0.
      StartInd_1=1; StopInd_1 =0
      open(file_id,file="air.dat")
!
      if (lroot) print*, 'the following parameters and '//&
          'species are found in air.dat (volume fraction fraction in %): '
!
      dataloop: do
!
        read(file_id,'(80A)',IOSTAT=iostat) ChemInpLine(1:80)
        if (iostat < 0) exit dataloop
        emptyFile=.false.
        StartInd_1=1; StopInd_1=0
        StopInd_1=index(ChemInpLine,' ') 
        specie_string=trim(ChemInpLine(1:StopInd_1-1)) 
        tmp_string=trim(ChemInpLine(1:1)) 
!
        if (tmp_string == '!' .or. tmp_string == ' ') then
        elseif (tmp_string == 'T') then
          StartInd=1; StopInd =0
!
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
          StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
!
          read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)'), TT
          if (lroot) print*, ' Temperature, K   ', TT
!
        elseif (tmp_string == 'P') then
!
          StartInd=1; StopInd =0
!
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
          StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
!
          read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)'), PP
          if (lroot) print*, ' Pressure, Pa   ', PP
!
        else
!
          call find_species_index(specie_string,ind_glob,ind_chem,found_specie)
!
          if (found_specie) then
!
            StartInd=1; StopInd =0
!
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
            StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
            read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)'), YY_k  
            if (lroot) print*, ' volume fraction, %,    ', YY_k, &
                species_constants(ind_chem,imass)
!
            air_mass=air_mass+YY_k*0.01/species_constants(ind_chem,imass)
!
            if (StartInd==80) exit
!
            stor1(k)=ind_chem
            stor2(k)=YY_k
            k=k+1
          endif
!
        endif
      enddo dataloop
!
!  Stop if air.dat is empty
!
      if (emptyFile)  call stop_it('The input file tran.dat was empty!')
      air_mass=1./air_mass
!
      do j=1,k-1 
        f(:,:,:,ichemspec(stor1(j)))=stor2(j)*0.01

      enddo
!
      !do j=1,nchemspec
      ! if (maxval(f(:,:,:,ichemspec(j)))<1e-15) then
      !     f(:,:,:,ichemspec(j))=1e-15
      !  endif 
      !enddo
!
      sum_Y=0.
!
      do j=1,nchemspec
        sum_Y=sum_Y+f(:,:,:,ichemspec(j))
      enddo
      do j=1,nchemspec
        f(:,:,:,ichemspec(j))=f(:,:,:,ichemspec(j))/sum_Y
        initial_massfractions(j)=f(l1,m1,n1,ichemspec(j))
      enddo
!
      if (mvar < 5) then
        call fatal_error("air_field", "I can only set existing fields")
      endif
      f(:,:,:,ilnTT)=log(TT)
!
      f(:,:,:,ilnrho)=log((PP*10./(k_B_cgs/m_u_cgs)*&
          air_mass/TT)/unit_mass*unit_length**3)

   !   f(:,:,:,iux)=init_ux

!
      if (lroot) print*, 'Air temperature, K', TT
      if (lroot) print*, 'Air pressure, dyn', PP*10.
      if (lroot) print*, 'Air density, g/cm^3:'
      if (lroot) print '(E10.3)',  PP*10./(k_B_cgs/m_u_cgs)*air_mass/TT
      if (lroot) print*, 'Air mean weight, g/mol', air_mass
      if (lroot) print*, 'R', k_B_cgs/m_u_cgs
!
      close(file_id)
!
    endsubroutine air_field
!********************************************************************
!!!!!!!!!  NSCBC boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine bc_nscbc_subin_x(f,df,topbot,val)
!
!   nscbc case 
!   subsonic inflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    ux,uy,T are fixed, drho  is calculated 
!
!   16-nov-08/natalia: coded.
!
      use MpiComm, only: stop_it
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil
!      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: mom2
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension(ny,nz) :: dux_dx, L_1, L_2, L_5, dpp_dx
      real, dimension(my,mz) :: rho0, gamma0, dmom2_dy, TT0 
      real, dimension (my,mz) :: cs2x,cs0_ar,cs20_ar
      real, dimension (mx,my,mz) :: cs2_full, gamma_full, rho_full, pp
      real, dimension(nchemspec) :: YYi
      real, dimension (mcom), optional :: val
      integer :: lll,i, sgn,k
      real :: Mach_num, u_t, T_t
!
      intent(inout) :: f
      intent(out) :: df
!

 !    if (.not.present(val)) call stop_it(&
 !          'bc_nscbc_subin_x: you must specify fbcx)')

      u_t=val(iux)
      T_t=val(ilnTT)
      do k=1,nchemspec
       YYi(k)=val(ichemspec(k))
      enddo

      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case(topbot)
      case('bot')
        lll = l1
        sgn = 1
      case('top')
        lll = l2
        sgn = -1
      case default
        print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect
!
      if (leos_chemistry) then
         cs20_ar=cs2_full(lll,:,:)
         cs0_ar=cs2_full(lll,:,:)**0.5
         gamma0=gamma_full(lll,:,:)
         TT0=exp(f(lll,:,:,ilnTT))
         rho_full=exp(f(:,:,:,ilnrho))
         rho0(:,:) = rho_full(lll,:,:)
         mom2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)
         do i=1,my
         do k=1,mz
          if (minval(gamma_full(:,i,k))<=0.) then
           pp(:,i,k)=0.
          else
           pp(:,i,k)=cs2_full(:,i,k)*rho_full(:,i,k)/gamma_full(:,i,k)
          endif
         enddo
         enddo
      else
        print*,"bc_nscbc_subin_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented for "//&
            "the chemistry case." 
        print*,"Boundary treatment skipped."
        return
      endif
!
!
      !  call der_onesided_4_slice(f,sgn,ilnrho,dlnrho_dx,lll,1)
        call der_onesided_4_slice(f,sgn,iux,dux_dx,lll,1)
        call der_onesided_4_slice(pp,sgn,dpp_dx,lll,1)
!
        do i=1,mz
          call der_pencil(2,mom2(lll,:,i),dmom2_dy(:,i))
        enddo

        select case(topbot)
        case('bot')
          L_1 = (f(lll,m1:m2,n1:n2,iux) - cs0_ar(m1:m2,n1:n2))*&
              (dpp_dx - rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
          L_5 =L_1-2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*&
              df(lll,m1:m2,n1:n2,iux)
        case('top')
          L_5 = (f(lll,m1:m2,n1:n2,iux) + cs0_ar(m1:m2,n1:n2))*&
              (dpp_dx + rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
          L_1 = L_5+2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*&
              df(lll,m1:m2,n1:n2,iux)
        endselect
        L_2 = 0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1) &
            +rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2)*df(lll,m1:m2,n1:n2,ilnTT)
        if (ldensity_nolog) then
          df(lll,m1:m2,n1:n2,ilnrho) = -1./cs20_ar(m1:m2,n1:n2)*&
              (L_2+0.5*(L_5 + L_1)) ! -dmom2_dy(m1:m2,n1:n2)
        else
          df(lll,m1:m2,n1:n2,ilnrho) = &
              -1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2) &
              *(L_2+0.5*(L_5 + L_1)) !&
          ! -1./rho0(m1:m2,n1:n2)*dmom2_dy(m1:m2,n1:n2)
        endif

!
! falia: this subroutine is still under construction
! Please, do not remove this commented part !!!!!
!
  
!
!  this conditions can be important! 
!  check withour them
!

   !     do k=1,nchemspec
   !      f(lll,m1:m2,n1:n2,ichemspec(k))=YYi(k)  
   !     enddo
   !      f(lll,m1:m2,n1:n2,iux) = u_t
   !     f(lll,m1:m2,n1:n2,ilnTT) = T_t
   !      df(lll,m1:m2,n1:n2,ilnTT)=0.


    endsubroutine bc_nscbc_subin_x
!***********************************************************
    subroutine bc_nscbc_nref_subout_x_(f,df,topbot)
!
!   nscbc case 
!   subsonic non-reflecting outflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    drho. dT, dux, duy  are calculated, p_inf can be 
!   fixed (if nscbc_sigma <>0)
!
!   16-nov-08/natalia: coded.
!
      use MpiComm, only: stop_it
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil
!      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (my,mz) :: rho0,gamma0
      real, dimension (mx,my,mz) :: mom2, rho_ux2, rho_uy2
      real, dimension (mx,my,mz) :: rho_gamma, rhoE_p, pp
      real, dimension (ny,nz) :: dux_dx,duy_dx,duz_dx, drho_dx, dpp_dx,dYk_dx
      real, dimension (ny,nz) :: drho_prefac,p_infx, KK, L_1, L_2, L_3,L_4, L_5
      real, dimension (my,mz) :: cs2x,cs0_ar,cs20_ar,dmom2_dy
      real, dimension (my,mz) :: drhoE_p_dy,dYk_dy, dux_dy
      real, dimension (ny,nz,nchemspec) :: bound_rhs_Y
      real, dimension (ny,nz) :: bound_rhs_T
      real, dimension (mx,my,mz) :: cs2_full, gamma_full, rho_full
!      real, dimension(ny,nz,3,3) :: dui_dxj
!
      integer :: lll, sgn,i,j,k
      real :: Mach_num
!
      intent(inout) :: f
      intent(out) :: df
!
      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case(topbot)
      case('bot')
        lll = l1
        sgn = 1
        if (leos_chemistry) then 
          call get_p_infx(p_infx,'bot')
        endif
      case('top')
        lll = l2
        sgn = -1
        if (leos_chemistry) then
          call get_p_infx(p_infx,'top')
        endif
      case default
        print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect

      if (leos_chemistry) then
         cs20_ar=cs2_full(lll,:,:)
         cs0_ar=cs2_full(lll,:,:)**0.5
         gamma0=gamma_full(lll,:,:)

        if (ldensity_nolog) then
          rho_full = f(:,:,:,ilnrho)
          rho0(:,:) = f(lll,:,:,ilnrho)
          drho_prefac=-1./cs20_ar(m1:m2,n1:n2)
        else
          rho_full = exp(f(:,:,:,ilnrho))
          rho0(:,:) = rho_full(lll,:,:)
          drho_prefac=-1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2)
        endif

         do i=1,my
         do k=1,mz
          if (minval(gamma_full(:,i,k))<=0.) then
           pp(:,i,k)=0.
          else
           pp(:,i,k)=cs2_full(:,i,k)*rho_full(:,i,k)/gamma_full(:,i,k)
          endif
         enddo
         enddo

         mom2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)
         rho_ux2(lll,:,:)=rho0(:,:)*f(lll,:,:,iux)*f(lll,:,:,iux)
         rho_uy2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)*f(lll,:,:,iuy)
         do i=1,my
         do k=1,mz
          if (gamma0(i,k)<=0.) then
           rho_gamma(lll,i,k)=0.
           rhoE_p(lll,i,k)=0.
          else
           rho_gamma(lll,i,k)=rho0(i,k)/gamma0(i,k)
           rhoE_p(lll,i,k)=0.5*rho_ux2(lll,i,k)+0.5*rho_uy2(lll,i,k) &
                         +cs20_ar(i,k)/(gamma0(i,k)-1.)*rho0(i,k)
          endif
         enddo
         enddo
      else
        print*,"bc_nscbc_subin_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented "//&
            "for the chemistry case." 
        print*,"Boundary treatment skipped."
        return
      endif

      call der_onesided_4_slice(rho_full,sgn,drho_dx,lll,1)
      call der_onesided_4_slice(pp,sgn,dpp_dx,lll,1)
      call der_onesided_4_slice(f,sgn,iux,dux_dx,lll,1)
      
      call der_onesided_4_slice(f,sgn,iuy,duy_dx,lll,1) 
      call der_onesided_4_slice(f,sgn,iuz,duz_dx,lll,1)

  !    call der_onesided_4_slice(f,sgn,ilnrho,div_rho(:,:,1),lll,1)
  !    call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,1),lll,1)
  !    call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,1),lll,1)
  !    call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,1),lll,1)

      do i=1,mz
        call der_pencil(2,mom2(lll,:,i),dmom2_dy(:,i))
        call der_pencil(2,f(lll,:,i,iux),dux_dy(:,i))
        call der_pencil(2,rhoE_p(lll,:,i),drhoE_p_dy(:,i))
      enddo
!
      Mach_num=maxval(f(lll,m1:m2,n1:n2,iux)/cs0_ar(m1:m2,n1:n2))
      KK=nscbc_sigma_out*(1.-Mach_num*Mach_num)*cs0_ar(m1:m2,n1:n2)/Lxyz(1)

!
      select case(topbot)
      case('bot')
        L_5=KK*(cs20_ar(m1:m2,n1:n2)/gamma0(m1:m2,n1:n2)*&
            rho0(m1:m2,n1:n2)-p_infx)
        L_1 = (f(lll,m1:m2,n1:n2,iux) - cs0_ar(m1:m2,n1:n2))*&
            (dpp_dx- rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
        call get_rhs_Y('bot',1,bound_rhs_Y)
        call get_rhs_T('bot',1,bound_rhs_T)
     case('top')
        L_1=KK*(cs20_ar(m1:m2,n1:n2)/gamma0(m1:m2,n1:n2)*&
            rho0(m1:m2,n1:n2)-p_infx)
        L_5 = (f(lll,m1:m2,n1:n2,iux) + cs0_ar(m1:m2,n1:n2))*&
            ( dpp_dx+ rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
        call get_rhs_Y('top',1,bound_rhs_Y)
        call get_rhs_T('top',1,bound_rhs_T)
      endselect
!
      L_2 = f(lll,m1:m2,n1:n2,iux)*(cs20_ar(m1:m2,n1:n2)*drho_dx-dpp_dx)
      L_3 = f(lll,m1:m2,n1:n2,iux)*duy_dx(:,:)
      L_4 = f(lll,m1:m2,n1:n2,iux)*duz_dx(:,:)
!
      if (ldensity_nolog) then
        df(lll,m1:m2,n1:n2,ilnrho) = &
            drho_prefac*(L_2+0.5*(L_5 + L_1))!-dmom2_dy(m1:m2,n1:n2)
      else
        df(lll,m1:m2,n1:n2,ilnrho) = &
            drho_prefac*(L_2+0.5*(L_5 + L_1)) !&
     !   -1./rho0(m1:m2,n1:n2)*dmom2_dy(m1:m2,n1:n2)
      endif
      df(lll,m1:m2,n1:n2,iux) = -1./&
          (2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2))*(L_5 - L_1) !&
     !      -f(lll,m1:m2,n1:n2,iux)*dux_dy(m1:m2,n1:n2)
      df(lll,m1:m2,n1:n2,iuy) = -L_3
      df(lll,m1:m2,n1:n2,iuz) = -L_4
      df(lll,m1:m2,n1:n2,ilnTT) = -1./&
          (rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))*(-L_2 &
          +0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1)) 
 !       -1./(rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))* &
 !       (gamma0(m1:m2,n1:n2)-1.) &
 !      *gamma0(m1:m2,n1:n2)*drhoE_p_dy(m1:m2,n1:n2) !&
    !    +bound_rhs_T(:,:)



      if (nchemspec>1) then
       do k=1,nchemspec
          call der_onesided_4_slice(f,sgn,ichemspec(k),dYk_dx,lll,1)
          do i=1,mz
            call der_pencil(2,f(lll,:,i,ichemspec(k)),dYk_dy(:,i))
          enddo
          df(lll,m1:m2,n1:n2,ichemspec(k))=&
              -f(lll,m1:m2,n1:n2,iux)*dYk_dx &
              -f(lll,m1:m2,n1:n2,iuy)*dYk_dy(m1:m2,n1:n2) &
              +bound_rhs_Y(:,:,k)

       ! if (lfilter) then
         do i=m1,m2
         do j=n1,n2
           if ((f(lll,i,j,ichemspec(k))+df(lll,i,j,ichemspec(k))*dt)<-1e-25 ) then
             df(lll,i,j,ichemspec(k))=-1e-25*dt
           endif
           if ((f(lll,i,j,ichemspec(k))+df(lll,i,j,ichemspec(k))*dt)>1.) then
             f(lll,i,j,ichemspec(k))=1.*dt
           endif
         enddo
         enddo
       ! endif
       enddo
      endif
!
      select case(topbot)
      case('bot')
       do i=1,nghost; f(l1-i,:,:,ilnTT)=2*f(l1,:,:,ilnTT)-f(l1+i,:,:,ilnTT); enddo 
       do i=1,nghost; f(l1-i,:,:,ilnrho)=2*f(l1,:,:,ilnrho)-f(l1+i,:,:,ilnrho); enddo 
       do i=1,nghost; f(l1-i,:,:,iux)=2*f(l1,:,:,iux)-f(l1+i,:,:,iux); enddo 
       do k=1,nchemspec
        do i=1,nghost; f(l1-i,:,:,ichemspec(k))=2*f(l1,:,:,ichemspec(k))  &
                            -f(l1+i,:,:,ichemspec(k)); enddo 
       enddo
      case('top')
       do i=1,nghost; f(l2+i,:,:,ilnTT)=2*f(l2,:,:,ilnTT)-f(l2-i,:,:,ilnTT); enddo
       do i=1,nghost; f(l2+i,:,:,ilnrho)=2*f(l2,:,:,ilnrho)-f(l2-i,:,:,ilnrho); enddo
      do i=1,nghost; f(l2+i,:,:,iux)=2*f(l2,:,:,iux)-f(l2-i,:,:,iux); enddo
       do k=1,nchemspec
        do i=1,nghost; f(l2+i,:,:,ichemspec(k))=2*f(l2,:,:,ichemspec(k))  &
                            -f(l2-i,:,:,ichemspec(k)); enddo 
       enddo
     case default
       print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect


    endsubroutine bc_nscbc_nref_subout_x_
!***********************************************************************

!***********************************************************
    subroutine bc_nscbc_nref_subout_x(f,df,topbot)

!
!   nscbc case 
!   subsonic non-reflecting outflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    drho. dT, dux, duy  are calculated, p_inf can be 
!   fixed (if nscbc_sigma <>0)
!
!   16-nov-08/natalia: coded.
!
      use MpiComm, only: stop_it
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil
!      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (my,mz) :: rho0,gamma0
      real, dimension (mx,my,mz) :: mom2, rho_ux2, rho_uy2
      real, dimension (mx,my,mz) :: rho_gamma, rhoE_p, pp
      real, dimension (ny,nz) :: dux_dx,duy_dx,duz_dx, drho_dx, dpp_dx,dYk_dx
      real, dimension (ny,nz) :: drho_prefac,p_infx, KK, L_1, L_2, L_3,L_4, L_5
      real, dimension (my,mz) :: cs2x,cs0_ar,cs20_ar,dmom2_dy
      real, dimension (my,mz) :: drhoE_p_dy,dYk_dy, dux_dy
      real, dimension (ny,nz,nchemspec) :: bound_rhs_Y
      real, dimension (ny,nz) :: bound_rhs_T
      real, dimension (mx,my,mz) :: cs2_full, gamma_full, rho_full
!      real, dimension(ny,nz,3,3) :: dui_dxj
!
      integer :: lll, sgn,i,j,k
      real :: Mach_num
!
      intent(inout) :: f
      intent(out) :: df
!
      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case(topbot)
      case('bot')
        lll = l1
        sgn = 1
        if (leos_chemistry) then 
          call get_p_infx(p_infx,'bot')
        endif
      case('top')
        lll = l2
        sgn = -1
        if (leos_chemistry) then
          call get_p_infx(p_infx,'top')
        endif
      case default
        print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect

      if (leos_chemistry) then
         cs20_ar=cs2_full(lll,:,:)
         cs0_ar=cs2_full(lll,:,:)**0.5
         gamma0=gamma_full(lll,:,:)

        if (ldensity_nolog) then
          rho_full = f(:,:,:,ilnrho)
          rho0(:,:) = f(lll,:,:,ilnrho)
          drho_prefac=-1./cs20_ar(m1:m2,n1:n2)
        else
          rho_full = exp(f(:,:,:,ilnrho))
          rho0(:,:) = rho_full(lll,:,:)
          drho_prefac=-1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2)
        endif

         do i=1,my
         do k=1,mz
          if (minval(gamma_full(:,i,k))<=0.) then
           pp(:,i,k)=0.
          else
           pp(:,i,k)=cs2_full(:,i,k)*rho_full(:,i,k)/gamma_full(:,i,k)
          endif
         enddo
         enddo

         mom2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)
         rho_ux2(lll,:,:)=rho0(:,:)*f(lll,:,:,iux)*f(lll,:,:,iux)
         rho_uy2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)*f(lll,:,:,iuy)
         do i=1,my
         do k=1,mz
          if (gamma0(i,k)<=0.) then
           rho_gamma(lll,i,k)=0.
           rhoE_p(lll,i,k)=0.
          else
           rho_gamma(lll,i,k)=rho0(i,k)/gamma0(i,k)
           rhoE_p(lll,i,k)=0.5*rho_ux2(lll,i,k)+0.5*rho_uy2(lll,i,k) &
                         +cs20_ar(i,k)/(gamma0(i,k)-1.)*rho0(i,k)
          endif
         enddo
         enddo
      else
        print*,"bc_nscbc_subin_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented "//&
            "for the chemistry case." 
        print*,"Boundary treatment skipped."
        return
      endif

      call der_onesided_4_slice(rho_full,sgn,drho_dx,lll,1)
      call der_onesided_4_slice(pp,sgn,dpp_dx,lll,1)
      call der_onesided_4_slice(f,sgn,iux,dux_dx,lll,1)
      
      call der_onesided_4_slice(f,sgn,iuy,duy_dx,lll,1) 
      call der_onesided_4_slice(f,sgn,iuz,duz_dx,lll,1)

  !    call der_onesided_4_slice(f,sgn,ilnrho,div_rho(:,:,1),lll,1)
  !    call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,1),lll,1)
  !    call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,1),lll,1)
  !    call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,1),lll,1)

      do i=1,mz
        call der_pencil(2,mom2(lll,:,i),dmom2_dy(:,i))
        call der_pencil(2,f(lll,:,i,iux),dux_dy(:,i))
        call der_pencil(2,rhoE_p(lll,:,i),drhoE_p_dy(:,i))
      enddo
!
      Mach_num=maxval(f(lll,m1:m2,n1:n2,iux)/cs0_ar(m1:m2,n1:n2))
      KK=nscbc_sigma_out*(1.-Mach_num*Mach_num)*cs0_ar(m1:m2,n1:n2)/Lxyz(1)

!
      select case(topbot)
      case('bot')
        L_5=KK*(cs20_ar(m1:m2,n1:n2)/gamma0(m1:m2,n1:n2)*&
            rho0(m1:m2,n1:n2)-p_infx)
        L_1 = (f(lll,m1:m2,n1:n2,iux) - cs0_ar(m1:m2,n1:n2))*&
            (dpp_dx- rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
        call get_rhs_Y('bot',1,bound_rhs_Y)
        call get_rhs_T('bot',1,bound_rhs_T)
     case('top')
        L_1=KK*(cs20_ar(m1:m2,n1:n2)/gamma0(m1:m2,n1:n2)*&
            rho0(m1:m2,n1:n2)-p_infx)
        L_5 = (f(lll,m1:m2,n1:n2,iux) + cs0_ar(m1:m2,n1:n2))*&
            ( dpp_dx+ rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
        call get_rhs_Y('top',1,bound_rhs_Y)
        call get_rhs_T('top',1,bound_rhs_T)
      endselect
!
      L_2 = f(lll,m1:m2,n1:n2,iux)*(cs20_ar(m1:m2,n1:n2)*drho_dx-dpp_dx)
      L_3 = f(lll,m1:m2,n1:n2,iux)*duy_dx(:,:)
      L_4 = f(lll,m1:m2,n1:n2,iux)*duz_dx(:,:)
!
      if (ldensity_nolog) then
        df(lll,m1:m2,n1:n2,ilnrho) = &
            drho_prefac*(L_2+0.5*(L_5 + L_1))!-dmom2_dy(m1:m2,n1:n2)
      else
        df(lll,m1:m2,n1:n2,ilnrho) = &
            drho_prefac*(L_2+0.5*(L_5 + L_1)) !&
     !   -1./rho0(m1:m2,n1:n2)*dmom2_dy(m1:m2,n1:n2)
      endif
      df(lll,m1:m2,n1:n2,iux) = -1./&
          (2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2))*(L_5 - L_1) !&
     !      -f(lll,m1:m2,n1:n2,iux)*dux_dy(m1:m2,n1:n2)
      df(lll,m1:m2,n1:n2,iuy) = -L_3
      df(lll,m1:m2,n1:n2,iuz) = -L_4
      df(lll,m1:m2,n1:n2,ilnTT) = -1./&
          (rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))*(-L_2 &
          +0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1)) 
 !       -1./(rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))* &
 !       (gamma0(m1:m2,n1:n2)-1.) &
 !      *gamma0(m1:m2,n1:n2)*drhoE_p_dy(m1:m2,n1:n2) !&
    !    +bound_rhs_T(:,:)



      if (nchemspec>1) then
       do k=1,nchemspec
          call der_onesided_4_slice(f,sgn,ichemspec(k),dYk_dx,lll,1)
          do i=1,mz
            call der_pencil(2,f(lll,:,i,ichemspec(k)),dYk_dy(:,i))
          enddo
          df(lll,m1:m2,n1:n2,ichemspec(k))=&
              -f(lll,m1:m2,n1:n2,iux)*dYk_dx &
              -f(lll,m1:m2,n1:n2,iuy)*dYk_dy(m1:m2,n1:n2) &
              +bound_rhs_Y(:,:,k)

       ! if (lfilter) then
         do i=m1,m2
         do j=n1,n2
           if ((f(lll,i,j,ichemspec(k))+df(lll,i,j,ichemspec(k))*dt)<-1e-25 ) then
             df(lll,i,j,ichemspec(k))=-1e-25*dt
           endif
           if ((f(lll,i,j,ichemspec(k))+df(lll,i,j,ichemspec(k))*dt)>1.) then
             f(lll,i,j,ichemspec(k))=1.*dt
           endif
         enddo
         enddo
       ! endif
       enddo
      endif
!
      select case(topbot)
      case('bot')
       do i=1,nghost; f(l1-i,:,:,ilnTT)=2*f(l1,:,:,ilnTT)-f(l1+i,:,:,ilnTT); enddo 
       do i=1,nghost; f(l1-i,:,:,ilnrho)=2*f(l1,:,:,ilnrho)-f(l1+i,:,:,ilnrho); enddo 
       do i=1,nghost; f(l1-i,:,:,iux)=2*f(l1,:,:,iux)-f(l1+i,:,:,iux); enddo 
       do k=1,nchemspec
        do i=1,nghost; f(l1-i,:,:,ichemspec(k))=2*f(l1,:,:,ichemspec(k))  &
                            -f(l1+i,:,:,ichemspec(k)); enddo 
       enddo
      case('top')
       do i=1,nghost; f(l2+i,:,:,ilnTT)=2*f(l2,:,:,ilnTT)-f(l2-i,:,:,ilnTT); enddo
       do i=1,nghost; f(l2+i,:,:,ilnrho)=2*f(l2,:,:,ilnrho)-f(l2-i,:,:,ilnrho); enddo
       do i=1,nghost; f(l2+i,:,:,iux)=2*f(l2,:,:,iux)-f(l2-i,:,:,iux); enddo
       do k=1,nchemspec
        do i=1,nghost; f(l2+i,:,:,ichemspec(k))=2*f(l2,:,:,ichemspec(k))  &
                            -f(l2-i,:,:,ichemspec(k)); enddo 
       enddo
     case default
       print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect


       select case(topbot)
      case('bot')
       do i=1,nghost; f(l1-i,1:m1-1,1:n1-1,ilnTT) &
         =2*f(l1,1:m1-1,1:n1-1,ilnTT)-f(l1+i,1:m1-1,1:n1-1,ilnTT); enddo 
       do i=1,nghost; f(l1-i,m2+1:my,n2+1:mz,ilnTT) &
         =2*f(l1,m2+1:my,n2+1:mz,ilnTT)-f(l1+i,m2+1:my,n2+1:mz,ilnTT); enddo 
       do i=1,nghost; f(l1-i,1:m1-1,1:n1-1,ilnrho) &
          =2*f(l1,1:m1-1,1:n1-1,ilnrho)-f(l1+i,1:m1-1,1:n1-1,ilnrho); enddo 
       do i=1,nghost; f(l1-i,m2+1:my,n2+1:mz,ilnrho) &
          =2*f(l1,m2+1:my,n2+1:mz,ilnrho)-f(l1+i,m2+1:my,n2+1:mz,ilnrho); enddo 
       do i=1,nghost; f(l1-i,1:m1-1,1:n1-1,iux) &
          =2*f(l1,1:m1-1,1:n1-1,iux)-f(l1+i,1:m1-1,1:n1-1,iux); enddo 
       do i=1,nghost; f(l1-i,m2+1:my,n2+1:mz,iux) &
          =2*f(l1,m2+1:my,n2+1:mz,iux)-f(l1+i,m2+1:my,n2+1:mz,iux); enddo 

       do k=1,nchemspec
        do i=1,nghost; f(l1-i,1:m1-1,1:n1-1,ichemspec(k))  &
                            =2*f(l1,1:m1-1,1:n1-1,ichemspec(k))  &
                            -f(l1+i,1:m1-1,1:n1-1,ichemspec(k)); enddo 
       enddo

       do k=1,nchemspec
        do i=1,nghost; f(l1-i,m2+1:my,n2+1:mz,ichemspec(k))  &
                         =2*f(l1,m2+1:my,n2+1:mz,ichemspec(k))  &
                         -f(l1+i,m2+1:my,n2+1:mz,ichemspec(k)); enddo 
       enddo
      case('top')
       do i=1,nghost; f(l2+i,1:m1-1,1:n1-1,ilnTT) &
             =2*f(l2,1:m1-1,1:n1-1,ilnTT)-f(l2-i,1:m1-1,1:n1-1,ilnTT); enddo
       do i=1,nghost; f(l2+i,m2+1:my,n2+1:mz,ilnTT)  &
             =2*f(l2,m2+1:my,n2+1:mz,ilnTT)-f(l2-i,m2+1:my,n2+1:mz,ilnTT); enddo

       do i=1,nghost; f(l2+i,1:m1-1,1:n1-1,ilnrho)  &
             =2*f(l2,1:m1-1,1:n1-1,ilnrho)-f(l2-i,1:m1-1,1:n1-1,ilnrho); enddo
        do i=1,nghost; f(l2+i,m2+1:my,n2+1:mz,ilnrho) &
             =2*f(l2,m2+1:my,n2+1:mz,ilnrho)-f(l2-i,m2+1:my,n2+1:mz,ilnrho); enddo


       do i=1,nghost; f(l2+i,1:m1-1,1:n1-1,iux)  &
                       =2*f(l2,1:m1-1,1:n1-1,iux)-f(l2-i,1:m1-1,1:n1-1,iux); enddo
       do k=1,nchemspec
        do i=1,nghost; f(l2+i,m2+1:my,n2+1:mz,ichemspec(k))  &
                            =2*f(l2,m2+1:my,n2+1:mz,ichemspec(k))  &
                            -f(l2-i,m2+1:my,n2+1:mz,ichemspec(k)); enddo 
     enddo
     case default
       print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect


!
  !    select case(topbot)
  !    case('bot')
  !     do i=1,nghost; f(l1-i,:,:,ilnTT)=2*f(l1,:,:,ilnTT)-f(l1+i,:,:,ilnTT); enddo 
  !     do i=1,nghost; f(l1-i,:,:,ilnrho)=2*f(l1,:,:,ilnrho)-f(l1+i,:,:,ilnrho); enddo 
  !     do i=1,nghost; f(l1-i,:,:,iux)=2*f(l1,:,:,iux)-f(l1+i,:,:,iux); enddo 
  !     do k=1,nchemspec
  !      do i=1,nghost; f(l1-i,:,:,ichemspec(k))=2*f(l1,:,:,ichemspec(k))  &
  !                          -f(l1+i,:,:,ichemspec(k)); enddo 
  !     enddo
  !    case('top')
  !     do i=1,nghost; f(l2+i,:,:,ilnTT)=2*f(l2,:,:,ilnTT)-f(l2-i,:,:,ilnTT); enddo
  !     do i=1,nghost; f(l2+i,:,:,ilnrho)=2*f(l2,:,:,ilnrho)-f(l2-i,:,:,ilnrho); enddo
  !     do i=1,nghost; f(l2+i,:,:,iux)=2*f(l2,:,:,iux)-f(l2-i,:,:,iux); enddo
  !     do k=1,nchemspec
  !      do i=1,nghost; f(l2+i,:,:,ichemspec(k))=2*f(l2,:,:,ichemspec(k))  &
  !                          -f(l2-i,:,:,ichemspec(k)); enddo 
  !     enddo
  !   case default
  !     print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
  !    endselect


    endsubroutine bc_nscbc_nref_subout_x
!***********************************************************************
!***********************************************************************

!********************************************************************
  subroutine chemistry_clean_up()

  use Mpicomm


  if (allocated(Bin_diff_coef))  deallocate(Bin_diff_coef)
  if (allocated(stoichio))       deallocate(stoichio)
  if (allocated(Sijm))           deallocate(Sijm)
  if (allocated(Sijp))           deallocate(Sijp)
  if (allocated(kreactions_z))   deallocate(kreactions_z)
  if (allocated(kreactions_m))   deallocate(kreactions_m)
  if (allocated(kreactions_p))   deallocate(kreactions_p)
  if (allocated(reaction_name))  deallocate(reaction_name)
  if (allocated(B_n))            deallocate(B_n)
  if (allocated(alpha_n))        deallocate(alpha_n)
  if (allocated(E_an))           deallocate(E_an)
  if (allocated(low_coeff))      deallocate(low_coeff)
  if (allocated(high_coeff))     deallocate(high_coeff)
  if (allocated(troe_coeff))     deallocate(troe_coeff)
  if (allocated(a_k4))           deallocate(a_k4)
  if (allocated(Mplus_case))     deallocate(Mplus_case)


!
  endsubroutine chemistry_clean_up
!********************************************************************
endmodule Chemistry

