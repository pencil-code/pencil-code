! $Id$
!
!  This modules addes chemical species and reactions.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lchemistry = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED gTT(3); mu1; gamma
! PENCILS PROVIDED cv; cv1; cp; cp1; YY(nchemspec); gXXk(3,nchemspec); ghhk(3,nchemspec)
! PENCILS PROVIDED cs2; rho1gpp(3); glnpp(3);  nu; gradnu(3),gamma_m1, gradcp(3)
! PENCILS PROVIDED DYDt_reac(nchemspec); DYDt_diff(nchemspec)
! PENCILS PROVIDED lambda; glambda(3)
! PENCILS PROVIDED DYDt_reac(nchemspec); DYDt_diff(nchemspec)
!
!***************************************************************
module Chemistry
!
  use Cparam
  use Cdata
  use Messages
  use Mpicomm, only: stop_it
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'chemistry.h'

  real :: Rgas, Rgas_unit_sys=1.
  real, dimension (mx,my,mz) :: cp_full,cv_full,mu1_full, nu_full, pp_full
  real, dimension (mx,my,mz) :: lambda_full, rho_full!, nu_art_full=0.
  real, dimension (mx,my,mz,nchemspec) :: cv_R_spec_full
  real, dimension (mx,my,mz,nchemspec) ::  cp_R_spec, hhk_full
  real, dimension (mx,my,mz) ::  TT_full
  real, dimension (mx,my,mz) ::  e_int_full

  real :: lambda_const=impossible
  real :: visc_const=impossible
  real :: diffus_const=impossible
  real :: Sc_number=0.7, Pr_number=0.7
  real :: Cp_const=impossible
  real :: Cv_const=impossible
  real :: init_x1=-0.2,init_x2=0.2
  real :: init_y1=-0.2,init_y2=0.2
  real :: init_z1=-0.2,init_z2=0.2
  real :: init_TT1=400, init_TT2=2400., init_ux=0., init_uy=0., init_uz=0.
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
!  logical :: lnospec_eqns=.true.

  logical :: lheatc_chemistry=.true.

  logical :: BinDif_simple=.false.
  logical :: visc_simple=.false.
  logical :: lT_const=.false.

  logical :: lfilter=.false.
  logical :: lkreactions_profile=.false., lkreactions_alpha=.false.
  integer :: nreactions=0,nreactions1=0,nreactions2=0
  integer :: ll1,ll2,mm1,mm2,nn1,nn2
  real, allocatable, dimension(:) :: kreactions_profile_width, kreactions_alpha

  integer :: mreactions
  integer, allocatable, dimension(:,:) :: stoichio,Sijm,Sijp
  real,    allocatable, dimension(:,:) :: kreactions_z
  real,    allocatable, dimension(:)   :: kreactions_m,kreactions_p
  character (len=30),allocatable, dimension(:) :: reaction_name
  logical :: lT_tanh=.false.
  logical :: ldamp_zone_NSCBC=.false.
! 1step_test case

    logical :: l1step_test=.false., lflame_front=.false.
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
  character (len=labellen), allocatable, dimension (:) :: kreactions_profile

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
      init_y1,init_y2,init_z1,init_z2,&
      init_TT1,init_TT2,init_ux,init_uy,init_uz,l1step_test,Sc_number,init_pressure,lfix_Sc, str_thick, &
      lfix_Pr,lT_tanh,ldamp_zone_NSCBC, lT_const,lheatc_chemistry


! run parameters
  namelist /chemistry_run_pars/ &
      lkreactions_profile, lkreactions_alpha, &
      chem_diff,chem_diff_prefactor, nu_spec, ldiffusion, ladvection, &
      lreactions,lchem_cdtc,lheatc_chemistry, BinDif_simple, visc_simple, &
      lmobility,mobility, lfilter,lT_tanh,ldamp_zone_NSCBC
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

  integer :: idiag_lambdam=0
  integer :: idiag_num=0
  integer :: idiag_diff1m=0
  integer :: idiag_diff2m=0
  integer :: idiag_diff3m=0
  integer :: idiag_diff4m=0
  integer :: idiag_diff5m=0
  integer :: idiag_diff6m=0
  integer :: idiag_diff7m=0
  integer :: idiag_diff8m=0  
  integer :: idiag_diff9m=0
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
!  Identify version number (generated automatically by SVN).
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

!     if (lreactions==.false.) 
!     if (ladvection==.false.)
!     if (ldiffusion==.false.) then
!      lnospec_eqns=.true.   
!     endif
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
      use Initcond
!AB:  use InitialCondition, only: initial_condition_chemistry
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
        select case (initchem(j))
!
        case ('nothing')
          if (lroot .and. .not. lnothing) print*, 'init_chem: nothing '
          lnothing=.true.
        case ('constant')
          do k=1,nchemspec
            f(:,:,:,ichemspec(k))=amplchemk(k)
          enddo
        case ('positive-noise')
          do k=1,nchemspec
            call posnoise(amplchemk(k),f,ichemspec(k))
          enddo
        case ('positive-noise-rel')
          do k=1,nchemspec
            call posnoise_rel(amplchemk(k),amplchemk2(k),f,ichemspec(k))
          enddo
        case ('innerbox')
          do k=1,nchemspec
            call innerbox(amplchemk(k),amplchemk2(k),f,ichemspec(k),widthchem)
          enddo
        case ('cos2x_cos2y_cos2z')
          do k=1,nchemspec
            call cos2x_cos2y_cos2z(amplchemk(k),f,ichemspec(k))
          enddo
        case ('coswave-x')
          do k=1,nchemspec
            call coswave(amplchem,f,ichemspec(k),kx=kx_chem)
          enddo
        case ('gaussian-x')
          do k=1,nchemspec
            call gaussian(amplchem,f,ichemspec(k),kx=kx_chem)
          enddo
        case ('gaussian-pos')
          do k=1,nchemspec
            call gaussianpos(amplchemk(k),f,ichemspec(k),widthchem,kx_chem,ky_chem,kz_chem)
            print*,"c(",x(l1),",",y(m1),",",z(n1),")=", f(l1,m1,n1,ichemspec(k))                     
          enddo
        case ('hatwave-x')
          do k=1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),widthchem,kx=kx_chem)
          enddo
        case ('hatwave-y')
          do k=1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),widthchem,ky=ky_chem)
          enddo
        case ('hatwave-z')
          do k=1,nchemspec
            call hatwave(amplchem,f,ichemspec(k),widthchem,kz=kz_chem)
          enddo
        case ('air')
          if (lroot ) print*, 'init_chem: air '
           inquire(file='air.dat',exist=air_exist)
           if (air_exist) then
            call air_field(f)
           else
            call stop_it('there is no air.dat file')
           endif
        case ('flame_front')
          call flame_front(f)
        case ('flame_front_3D')
          call flame_front_3D(f)
        case ('flame_blob')
          call flame_blob(f)
        case ('flame_slab')
          call flame_slab(f)
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
!RP:  I comment the following ... why should one specie case uniformly 
!     initialized to 1 and override the manually entered conditions ????  
!
!      if (lone_spec) then
!        f(:,:,:,ichemspec(1))=1.
!        if (lroot) print*, 'initchem: this is one specie case'
!      endif
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
      lpenc_requested(i_gXXk)=.true.
      lpenc_requested(i_ghhk)=.true.
      lpenc_requested(i_cs2)=.true.
!
      lpenc_requested(i_DYDt_reac)=.true.
      lpenc_requested(i_DYDt_diff)=.true.
!
       if (lcheminp) then
         lpenc_requested(i_rho)=.true.
         lpenc_requested(i_lnrho)=.true.
         lpenc_requested(i_glnpp)=.true.
         lpenc_requested(i_mu1)=.true.
       !  lpenc_requested(i_gmu1)=.true.
         lpenc_requested(i_pp)=.true.
         lpenc_requested(i_cp)=.true.
         lpenc_requested(i_cp1)=.true.
         lpenc_requested(i_cv)=.true.
         lpenc_requested(i_cv1)=.true.
         lpenc_requested(i_gamma)=.true.
    
            
!
         lpenc_requested(i_DYDt_reac)=.true.
         lpenc_requested(i_DYDt_diff)=.true.
         if (lheatc_chemistry) then
           lpenc_requested(i_lambda)=.true.
           lpenc_requested(i_glambda)=.true.
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
      use Sub
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (nx,3) :: gXX_tmp, ghhk_tmp
!
      intent(in) :: f
      intent(inout) :: p
      integer :: k,i
!
!  Mass fraction YY
!
      if (lpencil(i_YY)) then
        do k=1,nchemspec;  p%YY(:,k)=f(l1:l2,m,n,ichemspec(k)); enddo
      endif
      
      if (lpencil(i_gXXk) .and. ldiffusion) then
       do k=1,nchemspec 
         call grad(XX_full(:,:,:,k),gXX_tmp)
         do i=1,3; p%gXXk(:,i,k)=gXX_tmp(:,i); enddo
       enddo
      endif

     if (lpencil(i_ghhk) .and. lreactions) then
       do k=1,nchemspec 
         call grad(hhk_full(:,:,:,k),ghhk_tmp)
         do i=1,3; p%ghhk(:,i,k)=ghhk_tmp(:,i); enddo
       enddo
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
!        if (lpencil(i_gmu1)) call grad(mu1_full,p%gmu1)
!
!  Mole fraction XX
!
        !if (lpencil(i_XX)) then
        !  do k=1,nchemspec
        !    p%XX(:,k)=p%YY(:,k)/species_constants(ichemspec(k),imass)/p%mu1
        !  enddo
        !endif
!_--------------------------------------------------------------
!   moved to eos_chemistry.f90 
!
!  Temperature
!
   !     if (lpencil(i_lnTT)) p%lnTT=f(l1:l2,m,n,ilnTT)
   !     if (lpencil(i_TT)) p%TT=exp(p%lnTT)
   !     if (lpencil(i_TT1)) p%TT1=1./p%TT!
!
!  Temperature laplacian and gradient
!
    !    if (lpencil(i_glnTT)) call grad(f,ilnTT,p%glnTT)
    !    if (lpencil(i_del2lnTT)) call del2(f,ilnTT,p%del2lnTT)

!-----------------------------------------------------------------
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
          !  p%rho1gpp(:,i) = p%gamma_inv*p%cs2*&
          !      (p%glnrho(:,i)+p%glnTT(:,i)+p%gmu1(:,i)/p%mu1(:))
          !enddo
        endif
!
! Gradient of the lnpp
!
       if (lpencil(i_glnpp)) then
            do i=1,3
             p%glnpp(:,i)=p%rho1gpp(:,i)*p%rho(:)/p%pp(:)
            enddo
       endif
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
!        if (lpencil(i_gradcp)) call grad(cp_full,p%gradcp)
!
!  Specific heat at constant volume (i.e. density)
!
        if (lpencil(i_cv)) p%cv = cv_full(l1:l2,m,n)
!
        if (lpencil(i_cv1)) p%cv1=1./p%cv
!

!
!  Polytropic index
!
        if (lpencil(i_gamma)) p%gamma = p%cp*p%cv1
  !      if (lpencil(i_gamma_inv)) p%gamma_inv = p%cv*p%cp1
  !      if (lpencil(i_gamma_m1)) p%gamma_m1 = p%gamma - 1
!
!  Sound speed
!
        !if (lpencil(i_cs2)) p%cs2=p%cp*p%TT*p%gamma_m1*p%mu1
!
!------------------------------------------
!  NB:this is wrong for the chemistry case
!------------------------------------------
!  WL: what is wrong?

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
! Calculate chethermal diffusivity
!
      if (lpenc_requested(i_lambda) .and. lheatc_chemistry) then
         p%lambda=lambda_full(l1:l2,m,n)
         if (lpenc_requested(i_glambda)) call grad(lambda_full,p%glambda)
      endif
!
!  Calculate grad(enthalpy)
!
    !  if (lpenc_requested(i_ghYrho)) then
    !   call grad(hYrho_full,p%ghYrho)
    !  endif

    !   if (lpenc_requested(i_ghYrho_uu)) then
    !    call dot_mn(p%ghYrho,p%uu,p%ghYrho_uu)
    !  endif

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
      integer :: k

      real :: mO2, mH2, mN2, mH2O
      real :: log_inlet_density, del
      integer :: i_H2, i_O2, i_H2O, i_N2, ichem_H2, ichem_O2, ichem_N2, ichem_H2O
      real :: initial_mu1, final_massfrac_O2
      logical :: found_specie
! 
      lflame_front=.true.
!
      call air_field(f)

   if (ltemperature_nolog) f(:,:,:,ilnTT)=log(f(:,:,:,ilnTT))

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
        if (lT_tanh) then
          del=init_x2-init_x1
          f(k,:,:,ilnTT)=f(k,:,:,ilnTT)+log((init_TT2+init_TT1)*0.5  &
              +((init_TT2-init_TT1)*0.5)  &
              *(exp(x(k)/del)-exp(-x(k)/del))/(exp(x(k)/del)+exp(-x(k)/del)))
        else
          if (x(k)<=init_x1) then
            f(k,:,:,ilnTT)=log(init_TT1)
          endif
          if (x(k)>=init_x2) then
            f(k,:,:,ilnTT)=log(init_TT2)
          endif
          if (x(k)>init_x1 .and. x(k)<init_x2) then
            f(k,:,:,ilnTT)=&
                log((x(k)-init_x1)/(init_x2-init_x1) &
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
!
    call calc_for_chem_mixture(f)

!
!  Find logaritm of density at inlet
!
      initial_mu1&
          =initial_massfractions(ichem_H2)/(mH2)&
          +initial_massfractions(ichem_O2)/(mO2)&
          +initial_massfractions(ichem_H2O)/(mH2O)&
          +initial_massfractions(ichem_N2)/(mN2)
      log_inlet_density=&
          log(init_pressure)-log(Rgas)-log(init_TT1)-log(initial_mu1)
print*,'inlet rho=', exp(log_inlet_density),'inlet mu=',1./initial_mu1

!
!
!  Initialize density
!
      f(l1:l2,m1:m2,n1:n2,ilnrho)=log(init_pressure)-log(Rgas)  &
          -f(l1:l2,m1:m2,n1:n2,ilnTT)-log(mu1_full(l1:l2,m1:m2,n1:n2))
!
!  Initialize velocity
!
      f(l1:l2,m1:m2,n1:n2,iux)=f(l1:l2,m1:m2,n1:n2,iux)  &
          +init_ux*exp(log_inlet_density)/exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
!

     call calc_for_chem_mixture(f)
!
!  Check if we want nolog of density or nolog of temperature
!
      if (ldensity_nolog) f(l1:l2,m1:m2,n1:n2,irho)=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))
      if (ltemperature_nolog) f(l1:l2,m1:m2,n1:n2,iTT)=exp(f(l1:l2,m1:m2,n1:n2,ilnTT))
!
    endsubroutine flame_front
!***********************************************************************
   subroutine flame_front_3D(f)
!
! 06.05.2009/Nils Erland L. Haugen: adapted from similar
!                                   routine in special/chem_stream.f90
! This routine set up the initial profiles used in 1D flame speed measurments
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j1,j2,j3

      real :: mO2, mH2, mN2, mH2O
      real :: log_inlet_density, del
      integer :: i_H2, i_O2, i_H2O, i_N2, ichem_H2, ichem_O2, ichem_N2, ichem_H2O
      real :: initial_mu1, final_massfrac_O2,RR,RR1,RR2
      logical :: found_specie
! 
      lflame_front=.true.
!
      call air_field(f)
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
     ! do k=1,mx
       
        do j3=nn1,nn2
        do j2=mm1,mm2
        do j1=ll1,ll2

       RR=0.
       RR1=0.
       Rr2=0.
       del=0.
       if (nxgrid >1) then

       RR=x(j1)**2
       RR1=(init_x1)**2+(init_y1)**2+(init_z1)**2
       RR2=(init_x2)**2+(init_y2)**2+(init_z2)**2
     !  del=((init_x2-init_x1)**2)
        del=init_x2-init_x1

!        Rad=x(j1)**2
       endif
       if (nygrid>1) then
       ! Rad=Rad+y(j2)**2

       RR=RR+y(j2)**2
       RR1=RR1+(init_y1)**2
       RR2=RR2+(init_y2)**2
     !    del=(del+(init_y2-init_y1)**2)
       endif
       if (nzgrid>1) then
      ! Rad=Rad+z(j3)**2
       RR=RR+z(j3)**2
       RR1=RR1+(init_z1)**2
       RR2=RR2+(init_z2)**2
     !     del=(del+(init_z2-init_z1)**2)
       endif


       RR=sqrt(RR)
       RR1=sqrt(RR1)
       RR2=sqrt(RR2)
     !  del=sqrt(del)

!
!  Initialize temperature
!
        if (lT_tanh) then
        !  del=init_x2-init_x1
          f(j1,j2,j3,ilnTT)=f(j1,j2,j3,ilnTT)+log((init_TT2+init_TT1)*0.5  &
              +((init_TT2-init_TT1)*0.5)  &
              *(exp(RR/del)-exp(-RR/del))/(exp(RR/del)+exp(-RR/del)))
        else


      !    if (RR>del) then
      !      f(j1,j2,j3,ilnTT)=f(j1,j2,j3,ilnTT)+log(init_TT1)

      !    else
!            f(j1,j2,j3,ilnTT)=f(j1,j2,j3,ilnTT)+&
!                log((RR-RR1)/(RR2-RR1) &
 !               *(init_TT2-init_TT1)+init_TT1)
     !      f(j1,j2,j3,ilnTT)=f(j1,j2,j3,ilnTT)+log(init_TT2)

     !     endif
     
     
      
        if (RR<del) then
          f(j1,j2,j3,ilnTT)=log(init_TT1)+log(4.)*((del-RR)/del)**2
         else
          f(j1,j2,j3,ilnTT)=log(init_TT1)
         endif

        endif



!
!  Initialize steam and hydrogen
!
        if (lT_tanh) then
     !     del=(init_x2-init_x1)/3.
             del=(RR2-RR1)/3.
          f(j1,j2,j3,i_H2)=(0.+f(l1,m1,n1,i_H2))*0.5  &
              +(0.-f(l1,m1,n1,i_H2))*0.5  &
              *(exp(RR/del)-exp(-RR/del))/(exp(RR/del)+exp(-RR/del))
!
          f(j1,j2,j3,i_H2O)=(f(l1,m1,n1,i_H2)/2.*18.+f(l1,m1,n1,i_H2O))*0.5  &
              +((f(l1,m1,n1,i_H2)/2.*18.-f(l1,m1,n1,i_H2O))*0.5)  &
              *(exp(RR/del)-exp(-RR/del))/(exp(RR/del)+exp(-RR/del))
!
        else
       !   if (RR<del) then
       !     f(j1,j2,j3,i_H2O)=initial_massfractions(ichem_H2)/mH2*mH2O 
       !     f(j1,j2,j3,i_H2)=initial_massfractions(ichem_H2)
       !  endif


          if (RR<del) then
          f(j1,j2,j3,i_H2)=initial_massfractions(ichem_H2)*(1.-((del-RR)/del)**2)
	  f(j1,j2,j3,i_H2O)=initial_massfractions(ichem_H2)/mH2*mH2O*(1.-((del-RR)/del)**2)
         else
          f(j1,j2,j3,i_H2)=initial_massfractions(ichem_H2)
	  f(j1,j2,j3,i_H2O)=initial_massfractions(ichem_H2)/mH2*mH2O
         endif
	  
        endif
!
!  Initialize oxygen
!
        if (lT_tanh) then
        !  del=(init_x2-init_x1)
            del=(RR2-RR1)/3.
          f(j1,j2,j3,i_O2)=(f(l2,m2,n2,i_O2)+f(l1,m1,n1,i_O2))*0.5  &
              +((f(l2,m2,n2,i_O2)-f(l1,m1,n1,i_O2))*0.5)  &
              *(exp(RR/del)-exp(-RR/del))/(exp(RR/del)+exp(-RR/del))
        else

       !   if (RR<del) then
       !     f(j1,j2,j3,i_O2)=final_massfrac_O2
       !  endif
          if (RR<del) then
            f(j1,j2,j3,i_O2)=(final_massfrac_O2-initial_massfractions(ichem_O2))*((del-RR)/del)**2  &
	        +initial_massfractions(ichem_O2)
          else

           f(j1,j2,j3,i_O2)=initial_massfractions(ichem_O2)

 !                 (RR-RR2)/(RR1-RR2) &
 !               *(initial_massfractions(ichem_O2)-final_massfrac_O2)&
 !              +final_massfrac_O2
          endif
        endif
      enddo
      enddo
      enddo
!
      call calc_for_chem_mixture(f)
!
!  Find logaritm of density at inlet
!
      initial_mu1&
          =initial_massfractions(ichem_H2)/(mH2)&
          +initial_massfractions(ichem_O2)/(mO2)&
          +initial_massfractions(ichem_H2O)/(mH2O)&
          +initial_massfractions(ichem_N2)/(mN2)
      log_inlet_density=&
          log(init_pressure)-log(Rgas)-log(init_TT1)-log(initial_mu1)
!
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
                +init_ux!*exp(log_inlet_density)/exp(f(j1,j2,j3,ilnrho))
            f(j1,j2,j3,iuy)=f(j1,j2,j3,iuy)  &
                +init_ux!*exp(log_inlet_density)/exp(f(j1,j2,j3,ilnrho))
            f(j1,j2,j3,iuz)=f(j1,j2,j3,iuz)  &
                +init_ux!*exp(log_inlet_density)/exp(f(j1,j2,j3,ilnrho))
!
          enddo
        enddo
      enddo
!
!  Check if we want nolog of density
!
      if (ldensity_nolog) f(:,:,:,irho)=exp(f(:,:,:,ilnrho))
!
    endsubroutine flame_front_3D
!***********************************************************************
    subroutine flame_blob(f)

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: mu1
      integer :: j1,j2,j3

      real :: mO2, mH2, mN2, mH2O
      integer :: i_H2, i_O2, i_H2O, i_N2, ichem_H2, ichem_O2, ichem_N2, ichem_H2O
      real :: initial_mu1, final_massfrac_O2
      logical :: found_specie

      real :: Rad, sz1,sz2
 
     lflame_front=.true.
!
      call air_field(f)
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
!  it is in air_field(f)
!___________________________________________

      call calc_for_chem_mixture(f)
!
!  Find logaritm of density at inlet
!
      initial_mu1&
          =initial_massfractions(ichem_H2)/(mH2)&
          +initial_massfractions(ichem_O2)/(mO2)&
          +initial_massfractions(ichem_H2O)/(mH2O)&
          +initial_massfractions(ichem_N2)/(mN2)

 
!  
!   log_inlet_density=&
!      log(init_pressure)-log(Rgas)-log(init_TT1)-log(initial_mu1)
    !  do j3=nn1,nn2
    !    do j2=mm1,mm2
    !      do j1=ll1,ll2

      
        do j3=1,mz
      
       do j2=1,my
      do j1=1,mx

  !     do j3=1,nxgrid+8
  !    do j2=1,nygrid+8
  !    do j1=1,nzgrid+8
      
        Rad=0.
       if (nxgrid >1) then
        Rad=x(j1)**2
       endif
       if (nygrid>1) then
        Rad=Rad+y(j2)**2
       endif
       if (nzgrid>1) then
        Rad=Rad+z(j3)**2
       endif

       Rad=(Rad)**0.5

   !  if (j3==4)    print*,Rad,Rad**0.5,x(j1),y(j2),'  ',j1,j2   

  !      Rad=abs(x(j1))
         if (Rad<0.2) then
          f(j1,j2,j3,ilnTT)=log(init_TT1)+log(3.5)*((0.2-Rad)/0.2)**2
         else
          f(j1,j2,j3,ilnTT)=log(init_TT1)
         endif
          mu1(j1,j2,j3)=f(j1,j2,j3,i_H2)/(2.*mH2)+f(j1,j2,j3,i_O2)/(2.*mO2) &
              +f(j1,j2,j3,i_H2O)/(2.*mH2+mO2)+f(j1,j2,j3,i_N2)/(2.*mN2)

         f(j1,j2,j3,ilnrho)=log(init_pressure)-log(Rgas)-f(j1,j2,j3,ilnTT)  &
              -log(mu1_full(j1,j2,j3))

      !  f(j1,j2,j3,ilnrho)=log(init_pressure)-log(Rgas)-f(l1,m1,n1,ilnTT)  &
      !      -log(mu1_full(l1,m1,n1))



 !      f(j1,j2,j3,ilnrho)=log(init_pressure)-log(Rgas)  &
 !          -f(j1,j2,j3,ilnTT)-log(mu1_full(j1,j2,j3))

!
!  Initialize velocity
!
            f(j1,j2,j3,iux)=f(j1,j2,j3,iux)  &
                +init_ux!*exp(log_inlet_density)/exp(f(j1,j2,j3,ilnrho))
            f(j1,j2,j3,iuy)=f(j1,j2,j3,iuy)+ init_uy
            f(j1,j2,j3,iuz)=f(j1,j2,j3,iuz)+ init_uz

           if (nxgrid==1) f(j1,j2,j3,iux)=0. 
           if (nygrid==1) f(j1,j2,j3,iuy)=0. 
           if (nzgrid==1) f(j1,j2,j3,iuz)=0. 

           if (nxgrid/=1) then
            sz1=(xyz0(1)+Lxyz(1)*0.15)
            sz2=(xyz0(1)+Lxyz(1)*(1.-0.15))
            if ((x(j1)<sz1) .or. (sz2<x(j1))) then
         !     f(j1,j2,j3,iux)=0.
            endif   
           endif

          if (nygrid/=1)   then
            sz1=(xyz0(2)+Lxyz(2)*0.15)
            sz2=(xyz0(2)+Lxyz(2)*(1.-0.15))
            if ((y(j2)<sz1) .or. (y(j2)>sz2)) then
          !     f(j1,j2,j3,iuy)=0.
            endif
          endif

          if (nzgrid/=1)  then
           sz1=(xyz0(3)+Lxyz(3)*0.15)
           sz2=(xyz0(3)+Lxyz(3)*(1.-0.15))
           if ((z(j3)<sz1) .or. (z(j3)>sz2)) then
           !   f(j1,j2,j3,iuz)=0.
           endif 
          endif

          enddo
        enddo
      enddo
!
!  Check if we want nolog of density
!
      if (ldensity_nolog) f(:,:,:,irho)=exp(f(:,:,:,ilnrho))
!
    endsubroutine flame_blob
!***********************************************************************
!***********************************************************************
    subroutine flame_slab(f)

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: mu1
      integer :: j1,j2,j3
      real :: mO2, mH2, mN2, mH2O
      integer :: i_H2, i_O2, i_H2O, i_N2, ichem_H2, ichem_O2, ichem_N2, ichem_H2O
      real :: initial_mu1, final_massfrac_O2
      logical :: found_specie
      real :: Rad

     lflame_front=.true.
!
      call air_field(f)
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
!  it is in air_field(f)
!___________________________________________

      call calc_for_chem_mixture(f)
!
!  Find logaritm of density at inlet
!
      initial_mu1&
          =initial_massfractions(ichem_H2)/(mH2)&
          +initial_massfractions(ichem_O2)/(mO2)&
          +initial_massfractions(ichem_H2O)/(mH2O)&
          +initial_massfractions(ichem_N2)/(mN2)

       do j1=1,mx
         Rad=abs(x(j1))
         if (Rad<0.2) then
          f(j1,:,:,ilnTT)=log(init_TT1)+log(3.5)*((0.2-Rad)/0.2)**2
         else
          f(j1,:,:,ilnTT)=log(init_TT1)
         endif
          mu1(j1,:,:)=f(j1,:,:,i_H2)/(2.*mH2)+f(j1,:,:,i_O2)/(2.*mO2) &
              +f(j1,:,:,i_H2O)/(2.*mH2+mO2)+f(j1,:,:,i_N2)/(2.*mN2)

         f(j1,:,:,ilnrho)=log(init_pressure)-log(Rgas)-f(j1,:,:,ilnTT)  &
              -log(mu1_full(j1,:,:))

!
!  Initialize velocity
!
            f(j1,:,:,iux)=f(j1,:,:,iux)  &
                +init_ux!*exp(log_inlet_density)/exp(f(j1,j2,j3,ilnrho))
            f(j1,:,:,iuy)=f(j1,:,:,iuy)+ init_uy
            f(j1,:,:,iuz)=f(j1,:,:,iuz)+ init_uz

           if (nxgrid==1) f(j1,:,:,iux)=0. 
           if (nygrid==1) f(j1,:,:,iuy)=0. 
           if (nzgrid==1) f(j1,:,:,iuz)=0. 

      enddo
!
!  Check if we want nolog of density
!
      if (ldensity_nolog) f(:,:,:,irho)=exp(f(:,:,:,ilnrho))
!
    endsubroutine flame_slab
!***********************************************************************
    subroutine calc_for_chem_mixture(f)
!
!  DOCUMENT ME!!!
!
      use Sub
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) ::  tmp_sum,tmp_sum2, nu_dyn
!
      intent(in) :: f
      integer :: k,j, j1,j2,j3
      real :: T_up, T_mid, T_low, T_loc, tmp
!
      logical :: tran_exist=.false.
      logical,save :: lwrite=.true.
!
      character (len=20) :: output_file="./data/mix_quant.out"
      integer :: file_id=123
!
! Density and temperature
!

     if (ldensity_nolog) then
      rho_full=f(:,:,:,ilnrho)
     else
      rho_full=exp(f(:,:,:,ilnrho))
     endif


     if (ltemperature_nolog) then
         TT_full=f(:,:,:,iTT)
     else
         TT_full=exp(f(:,:,:,ilnTT))
     endif
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
        !  if (lpencil(i_rho1gpp)) then
           do j2=mm1,mm2
           do j3=nn1,nn2
            pp_full(:,j2,j3)=Rgas*mu1_full(:,j2,j3) &
                            *rho_full(:,j2,j3)*TT_full(:,j2,j3)
           enddo
           enddo
        !  endif

!
!  Specific heat at constant pressure
!

        if ((Cp_const<impossible) .or. (Cv_const<impossible)) then

          if (Cp_const<impossible) then
            cp_full=Cp_const*mu1_full
            cv_full=(Cp_const-Rgas)*mu1_full
          endif

          if (Cv_const<impossible) then
            cv_full=Cv_const*mu1_full
          endif

         else
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
                    cv_R_spec_full(j1,j2,j3,k)=cp_R_spec(j1,j2,j3,k)-1.
                  elseif (T_loc >=T_mid .and. T_loc<= T_up) then
                !  elseif (TT_full(j1,j2,j3) >=T_mid ) then
                    tmp=0.
                    do j=1,5
                      tmp=tmp+species_constants(k,iaa1(j))*T_loc**(j-1)
                    enddo
                    cp_R_spec(j1,j2,j3,k)=tmp
                    cv_R_spec_full(j1,j2,j3,k)=cp_R_spec(j1,j2,j3,k)-1.
                  else
                    print*,'TT_full(j1,j2,j3)=',T_loc
                    print*,'j1,j2,j3=',j1,j2,j3
                    call fatal_error('calc_for_chem_mixture',&
                        'TT_full(j1,j2,j3) is outside range')
                  endif
                 cp_full(j1,j2,j3)=cp_full(j1,j2,j3)+f(j1,j2,j3,ichemspec(k))  &
                  *cp_R_spec(j1,j2,j3,k)/species_constants(k,imass)*Rgas
                 cv_full(j1,j2,j3)=cv_full(j1,j2,j3)+f(j1,j2,j3,ichemspec(k))  &
                  *cv_R_spec_full(j1,j2,j3,k)/species_constants(k,imass)*Rgas
                enddo
              enddo
            enddo
!
          enddo
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
      
       if (visc_const<impossible) then
                nu_full=visc_const
       else

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
!              tmp_sum2(j1,j2,j3)=0.
!              do j=1,k
!                mk_mj=species_constants(k,imass) &
!                    /species_constants(j,imass)
!                nuk_nuj(j1,j2,j3)=species_viscosity(j1,j2,j3,k) &
!                    /species_viscosity(j1,j2,j3,j)
!                Phi(j1,j2,j3)=1./sqrt(8.)*1./sqrt(1.+mk_mj) &
!                    *(1.+sqrt(nuk_nuj(j1,j2,j3))*mk_mj**(-0.25))**2
!                tmp_sum2(j1,j2,j3)=tmp_sum2(j1,j2,j3) &
!                                  +XX_full(j1,j2,j3,j)*Phi(j1,j2,j3)
!              enddo
              nu_dyn(j1,j2,j3)=nu_dyn(j1,j2,j3)+XX_full(j1,j2,j3,k)*&
                  species_viscosity(j1,j2,j3,k)!/tmp_sum2(j1,j2,j3)
             enddo

              nu_full(j1,j2,j3)=nu_dyn(j1,j2,j3)/rho_full(j1,j2,j3)
          endif

        enddo
        enddo
        enddo
       endif
!
!  Diffusion coeffisient of a mixture
!
         do j3=nn1,nn2
         do j2=mm1,mm2
         do j1=1,mx

            Diff_full(j1,j2,j3,:)=0.
            if (.not. lone_spec) then

             if (diffus_const<impossible) then
                Diff_full(j1,j2,j3,:)=diffus_const
             elseif (lfix_Sc) then
              do k=1,nchemspec
               Diff_full(j1,j2,j3,k)=species_viscosity(j1,j2,j3,k) &
                                    /rho_full(j1,j2,j3)/Sc_number
              enddo
             elseif (ldiffusion) then
!
! The mixture diffusion coefficient as described in eq. 5-45 of the Chemkin
! manual. Previously eq. 5-44 was used, but due to problems in the limit
! when the mixture becomes a pure specie we changed to the more robust eq. 5-45.
!
              do k=1,nchemspec
                tmp_sum(j1,j2,j3)=0.
                tmp_sum2(j1,j2,j3)=0.
                do j=1,nchemspec
                 if (j .ne. k) then
                   tmp_sum(j1,j2,j3)=tmp_sum(j1,j2,j3) &
                        +XX_full(j1,j2,j3,j)/Bin_Diff_coef(j1,j2,j3,j,k)
                   tmp_sum2(j1,j2,j3)=tmp_sum2(j1,j2,j3) &
                       +XX_full(j1,j2,j3,j)*species_constants(j,imass)

                 endif
                enddo
                 Diff_full(j1,j2,j3,k)=mu1_full(j1,j2,j3)*tmp_sum2(j1,j2,j3)&
                     /tmp_sum(j1,j2,j3)
              enddo
             
            endif
            endif

             do k=1,nchemspec
              Diff_full_add(j1,j2,j3,k)=Diff_full(j1,j2,j3,k)*&
                  species_constants(k,imass)/unit_mass &
                  *mu1_full(j1,j2,j3)!*0.8
             enddo
           
         enddo
         enddo
         enddo

!
!  Thermal diffusivity
!
        if (lheatc_chemistry) call calc_therm_diffus_coef(f)


!
!  Dimensionless Standard-state molar enthalpy H0/RT
!

         if (.not. lT_const) then
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
        if (lreactions) then
        do j3=nn1,nn2
        do j2=mm1,mm2
        do j1=1,mx
          do k=1,nchemspec
            hhk_full(j1,j2,j3,k)=H0_RT(j1,j2,j3,k)*Rgas*TT_full(j1,j2,j3)&
               /species_constants(k,imass)
          enddo
        enddo
        enddo
        enddo
        endif
       endif


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
        write(file_id,'(7E12.4)') rho_full(l1,m1,n1)*unit_mass/unit_length**3, &
                                  rho_full(l2,m2,n2)*unit_mass/unit_length**3
        write(file_id,*) ''
        write(file_id,*) 'Themperature, K'
         ! Commented the next line out because
         ! samples/2d-tests/chemistry_GrayScott apparently has no f(:,:,:,5)
        if (ilnTT>0) write(file_id,'(7E12.4)')  &
        exp(f(l1,m1,n1,ilnTT))*unit_temperature, &
        exp(f(l2,m2,n2,ilnTT))*unit_temperature
        write(file_id,*) ''
        write(file_id,*) 'Cp,  erg/mole/K'
        write(file_id,'(7E12.4)') cp_full(l1,m1,n1)/Rgas*&
            Rgas_unit_sys/mu1_full(l1,m1,n1)/unit_mass,cp_full(l2,m2,n2)/Rgas*&
            Rgas_unit_sys/mu1_full(l2,m2,n2)/unit_mass
        write(file_id,*) ''
        write(file_id,*) 'cp, erg/g/K'
        write(file_id,'(7E12.4)') cp_full(l1,m1,n1)/Rgas*Rgas_unit_sys,cp_full(l2,m2,n2)/Rgas*Rgas_unit_sys
        write(file_id,*) ''
        write(file_id,*) 'gamma,max,min'
        write(file_id,'(7E12.4)') cp_full(l1,m1,n1)/cv_full(l1,m1,n1),&
            cp_full(l2,m2,n2)/cv_full(l2,m2,n2)
        write(file_id,*) ''
        write(file_id,*) 'Viscosity, g/cm/s,'
        write(file_id,'(7E12.4)') nu_dyn(l1,m1,n1)*&
            (unit_mass/unit_length/unit_time),nu_dyn(l2,m1,n1)*&
            (unit_mass/unit_length/unit_time)
        write(file_id,*) ''
        write(file_id,*) 'Species viscosity, g/cm/s,'
        do k=1,nchemspec
        write(file_id,'(7E12.4)') species_viscosity(l1,m1,n1,k),  &
                                  species_viscosity(l2-1,m1,n1,k)
        enddo
        write(file_id,*) ''
        write(file_id,*) 'Thermal cond, erg/(cm K s),'
        write(file_id,'(7E12.4)') (lambda_full(l1,m1,n1)*&
            unit_energy/unit_time/unit_length/unit_temperature), &
                        (lambda_full(l2,m2,n2)*&
            unit_energy/unit_time/unit_length/unit_temperature)
        write(file_id,*) ''
        write(file_id,*) 'Species  Diffusion coefficient, cm^2/s'
        do k=1,nchemspec
        write(file_id,'(7E12.4)')Diff_full(l1,m1,n1,k)*unit_length**2/unit_time, &
                                 Diff_full(l2,m2,n2,k)*unit_length**2/unit_time
        enddo
        write(file_id,*) ''


        if (lroot) print*,'calc_for_chem_mixture: writing mix_quant.out file'
        close(file_id)
        lwrite=.false.
      endif
!!$      if (lwrite) then
!!$        open(file_id,file=output_file)
!!$        write(file_id,*) 'Mixture quantities'
!!$        write(file_id,*) '*******************'
!!$        write(file_id,*) ''
!!$        write(file_id,*) 'Mass, g/mole'
!!$        write(file_id,'(7E12.4)') 1./maxval(mu1_full/unit_mass)
!!$        write(file_id,*) ''
!!$        write(file_id,*) 'Density, g/cm^3'
!!$        write(file_id,'(7E12.4)') minval(rho_full(l1:l2,m1:m2,n1:n2))*unit_mass/unit_length**3, &
!!$                                  maxval(rho_full(l1:l2,m1:m2,n1:n2))*unit_mass/unit_length**3
!!$        write(file_id,*) ''
!!$        write(file_id,*) 'Themperature, K'
!!$         ! Commented the next line out because
!!$         ! samples/2d-tests/chemistry_GrayScott apparently has no f(:,:,:,5)
!!$        if (ilnTT>0) write(file_id,'(7E12.4)')  &
!!$        exp(minval(f(l1:l2,m1:m2,n1:n2,ilnTT)))*unit_temperature, &
!!$        exp(maxval(f(l1:l2,m1:m2,n1:n2,ilnTT)))*unit_temperature
!!$        write(file_id,*) ''
!!$        write(file_id,*) 'Cp,  erg/mole/K'
!!$        write(file_id,'(7E12.4)') maxval(cp_full(l1:l2,m1:m2,n1:n2))/Rgas*&
!!$            Rgas_unit_sys/maxval(mu1_full(l1:l2,m1:m2,n1:n2)/unit_mass)
!!$        write(file_id,*) ''
!!$        write(file_id,*) 'cp, erg/g/K'
!!$        write(file_id,'(7E12.4)') maxval(cp_full(l1:l2,m1:m2,n1:n2))/Rgas*Rgas_unit_sys
!!$        write(file_id,*) ''
!!$        write(file_id,*) 'gamma,max,min'
!!$        write(file_id,'(7E12.4)') maxval(cp_full(l1:l2,m1:m2,n1:n2)) &
!!$                                 /maxval(cv_full(l1:l2,m1:m2,n1:n2)),&
!!$          minval(cp_full(l1:l2,m1:m2,n1:n2))/minval(cv_full(l1:l2,m1:m2,n1:n2))
!!$        write(file_id,*) ''
!!$        write(file_id,*) 'Viscosity, g/cm/s,'
!!$        write(file_id,'(7E12.4)') minval(nu_dyn(l1:l2,m1:m2,n1:n2))*&
!!$            (unit_mass/unit_length/unit_time),maxval(nu_dyn(l1:l2,m1:m2,n1:n2))*&
!!$            (unit_mass/unit_length/unit_time)
!!$        write(file_id,*) ''
!!$        write(file_id,*) 'Species viscosity, g/cm/s,'
!!$        do k=1,nchemspec
!!$        write(file_id,'(7E12.4)') minval(species_viscosity(l1:l2,m1:m2,n1:n2,k)),  &
!!$                                  maxval(species_viscosity(l1:l2,m1:m2,n1:n2,k))
!!$        enddo
!!$        write(file_id,*) ''
!!$        write(file_id,*) 'Thermal cond, erg/(cm K s),'
!!$        write(file_id,'(7E12.4)') (minval(lambda_full(l1:l2,m1:m2,n1:n2))*&
!!$            unit_energy/unit_time/unit_length/unit_temperature), &
!!$                        (maxval(lambda_full(l1:l2,m1:m2,n1:n2))*&
!!$            unit_energy/unit_time/unit_length/unit_temperature)
!!$        write(file_id,*) ''
!!$        write(file_id,*) 'Species  Diffusion coefficient, cm^2/s'
!!$        do k=1,nchemspec
!!$        write(file_id,'(7E12.4)')minval(Diff_full(l1:l2,m1:m2,n1:n2,k))*unit_length**2/unit_time, &
!!$                                 maxval(Diff_full(l1:l2,m1:m2,n1:n2,k))*unit_length**2/unit_time
!!$        enddo
!!$        write(file_id,*) ''
!!$
!!$
!!$        if (lroot) print*,'calc_for_chem_mixture: writing mix_quant.out file'
!!$        close(file_id)
!!$        lwrite=.false.
!!$      endif
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
      use General, only: chn
!
      character (len=80) :: chemicals='' 
      ! Careful, limits the absolut size of the input matrix !!!
      character (len=15) :: file1='chemistry_m.dat',file2='chemistry_p.dat'
!      character (len=20) :: input_file='chem.inp'
      real, dimension (mx,my,mz,mfarray) :: f
      real :: dummy
      logical :: exist,exist1,exist2
      integer :: i,j,stat
      integer :: nchemspectemp
      character :: tmpchar
      logical :: inside

!
!  Find number of reactions by reading how many lines we have in file2
!
      j=1
      open(19,file=file2)
      read(19,*) chemicals
      do while (.true.)
        read(19,*,end=996) dummy
        j=j+1
      enddo
996   close(19)
      mreactions=j-1
      if (lroot) print*,'Number of reactions=',mreactions
!
!  Find number of compounds by reading how many columns we have in file1
!
      open(19,file=file1)
      read(19,fmt="(a80)") chemicals
      nchemspectemp=0
      inside=.true.
      do i=1,len_trim(chemicals)
        tmpchar=chemicals(i:i)
        if (tmpchar == ' ') then
          if (.not. inside) then
            inside = .true.
            nchemspectemp=nchemspectemp+1
          endif
        else
          inside=.false.
        endif
      enddo
      if (inside) nchemspectemp=nchemspectemp-1
      close(19)
      if (lroot) print*,'Number of compounds=',nchemspectemp
      if (nchemspectemp>nchemspec) call &
          stop_it("Too much chemicals! Change NCHEMSPEC in src/cparam.local")
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
        allocate(kreactions_profile(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for kreactions_profile")
        allocate(kreactions_profile_width(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for kreactions_profile_width")
        allocate(kreactions_alpha(mreactions),STAT=stat)
        if (stat>0) call stop_it("Couldn't allocate memory for kreactions_alpha")
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
          read(19,*,end=994) kreactions_m(j),(Sijm(i,j),i=1,nchemspectemp)
        enddo
994     close(19)
        nreactions1=j-1
!
!  file2
!
        open(19,file=file2)
        read(19,*) chemicals
        do j=1,mreactions
          if (lkreactions_profile) then
            if (lkreactions_alpha) then
              read(19,*) kreactions_p(j),(Sijp(i,j),i=1,nchemspectemp),kreactions_profile(j),&
                  kreactions_profile_width(j),kreactions_alpha(j)
            else
              read(19,*) kreactions_p(j),(Sijp(i,j),i=1,nchemspectemp),kreactions_profile(j),&
                  kreactions_profile_width(j)
            endif
          else
            if (lkreactions_alpha) then
              read(19,*) kreactions_p(j),(Sijp(i,j),i=1,nchemspectemp),kreactions_alpha(j)
            else
              read(19,*) kreactions_p(j),(Sijp(i,j),i=1,nchemspectemp)
            endif
          endif
        enddo
        close(19)
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
        if (nreactions /= mreactions) call stop_it('nreactions/=mreactions')
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
!        print*,'chemicals=',chemicals
        print*,'kreactions_m=',kreactions_m(1:nreactions)
        print*,'kreactions_p=',kreactions_p(1:nreactions)
        print*,'Sijm:'
        do i=1,nreactions
          print*,Sijm(:,i)
        enddo
        print*,'Sijp:'
        do i=1,nreactions
          print*,Sijp(:,i)
        enddo
        print*,'stoichio='
        do i=1,nreactions
          print*,stoichio(:,i)
        enddo
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
          elseif (kreactions_profile(j)=='sin-bg') then
            do n=1,mz
              kreactions_z(n,j)=0.5*(1.1+cos(pi*z(n)/kreactions_profile_width(j)))
            enddo
          elseif (kreactions_profile(j)=='spike') then
            do n=1,mz
              if (cos(pi*z(n)/kreactions_profile_width(j)) > 0.99) then
                kreactions_z(n,j)=1
              else
                kreactions_z(n,j)=0
              endif
            enddo
          endif
        enddo
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine astrobiology_data
!**********************************************************************
    subroutine chemkin_data(f)
!
!  if the file with chemkin data exists
!
!  DOCUMENT ME FURTHER!!!
!
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
        allocate(Mplus_case (nreactions),STAT=stat)
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
!
      call keep_compiler_quiet(f)
!
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
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: gchemspec, dk_D
      real, dimension (nx) :: ugchemspec, sum_DYDT
      real, dimension (nx) :: sum_dk_ghk,dk_dhhk,sum_hhk_DYDt_reac
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
      !      real,dimension(nchemspec) :: reac_rate=0.
!
!      real,dimension(nx) :: sum_reac_rate
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
        if (lhydro.and.ladvection) then
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

        sum_DYDt(i)=-p%rho(1)*(p%TT(i)-Tinf)*p%TT1(i) &
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
        
        sum_hhk_DYDt_reac=0.
        sum_dk_ghk=0.
        
!NANANAN

        do k=1,nchemspec
          sum_DYDt=sum_DYDt+Rgas/species_constants(k,imass)&
              *(p%DYDt_reac(:,k)+p%DYDt_diff(:,k))
          if (lreactions) then
           sum_hhk_DYDt_reac=sum_hhk_DYDt_reac-hhk_full(l1:l2,m,n,k)*p%DYDt_reac(:,k)
          endif

         if (ldiffusion) then
          do i=1,3
           dk_D(:,i)=(p%gXXk(:,i,k) &
            +(XX_full(l1:l2,m,n,k)-f(l1:l2,m,n,ichemspec(k)))*p%glnpp(:,i)) &
            *Diff_full_add(l1:l2,m,n,k)
          enddo
           call dot_mn(dk_D,p%ghhk(:,:,k),dk_dhhk)
           sum_dk_ghk=sum_dk_ghk+dk_dhhk
         endif

        enddo
       endif
!
   


        if (l1step_test) then
          RHS_T_full(l1:l2,m,n)=sum_DYDt(:)
        else
          if (ltemperature_nolog) then
        !    RHS_T_full(l1:l2,m,n)=(sum_DYDt(:)- Rgas*p%mu1*p%divu)*p%cv1*p%TT(:)! &
           ! -(hYrho_full(l1:l2,m,n)*p%divu(:)+p%ghYrho_uu(:))*p%cv1/p%rho(:)
          else
            RHS_T_full(l1:l2,m,n)=(sum_DYDt(:)-Rgas*p%mu1*p%divu)*p%cv1 &
        
           +sum_dk_ghk*p%TT1(:)*p%cv1+sum_hhk_DYDt_reac*p%TT1(:)*p%cv1
           
          endif
        endif

        if (.not. lT_const) then
         df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + RHS_T_full(l1:l2,m,n)
        endif     

!
        if (lheatc_chemistry) call calc_heatcond_chemistry(f,df,p)
!
      endif




!
! this damping zone is needed in a case of NSCBC
!

      if (ldamp_zone_NSCBC) call damp_zone_for_NSCBC(f,df)
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

        if (idiag_lambdam/=0) call sum_mn_name(lambda_full(l1:l2,m,n),&
                             idiag_lambdam)
        if (idiag_num/=0) call sum_mn_name(nu_full(l1:l2,m,n),&
                             idiag_num)
        if (idiag_diff1m/=0) call sum_mn_name(diff_full(l1:l2,m,n,i1),&
                             idiag_diff1m)
        if (idiag_diff2m/=0) call sum_mn_name(diff_full(l1:l2,m,n,i2),&
                             idiag_diff2m)
        if (idiag_diff3m/=0) call sum_mn_name(diff_full(l1:l2,m,n,i3),&
                             idiag_diff3m)
        if (idiag_diff4m/=0) call sum_mn_name(diff_full(l1:l2,m,n,i4),&
                             idiag_diff4m)
        if (idiag_diff5m/=0) call sum_mn_name(diff_full(l1:l2,m,n,i5),&
                             idiag_diff5m)
        if (idiag_diff6m/=0) call sum_mn_name(diff_full(l1:l2,m,n,i6),&
                             idiag_diff6m)
        if (idiag_diff7m/=0) call sum_mn_name(diff_full(l1:l2,m,n,i7),&
                             idiag_diff7m)
        if (idiag_diff8m/=0) call sum_mn_name(diff_full(l1:l2,m,n,i8),&
                             idiag_diff8m)
        if (idiag_diff9m/=0) call sum_mn_name(diff_full(l1:l2,m,n,i9),&
                             idiag_diff9m)
!
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
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
      use Diagnostics
      use General, only: chn
!
      integer :: iname,inamez
      logical :: lreset,lwr
      logical, optional :: lwrite
      character (len=5) :: schemspec,snd1
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

        idiag_diff1m=0; idiag_diff2m=0; idiag_diff3m=0; idiag_diff4m=0;
        idiag_diff5m=0; idiag_diff6m=0; idiag_diff7m=0; idiag_diff8m=0; 
        idiag_diff9m=0; idiag_lambdam=0; idiag_num=0

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
        call parse_name(iname,cname(iname),cform(iname),'lambdam',idiag_lambdam)
        call parse_name(iname,cname(iname),cform(iname),'num',idiag_num)
        call parse_name(iname,cname(iname),cform(iname),'diff1m',idiag_diff1m)
        call parse_name(iname,cname(iname),cform(iname),'diff2m',idiag_diff2m)
        call parse_name(iname,cname(iname),cform(iname),'diff3m',idiag_diff3m)
        call parse_name(iname,cname(iname),cform(iname),'diff4m',idiag_diff4m)
        call parse_name(iname,cname(iname),cform(iname),'diff5m',idiag_diff5m)
        call parse_name(iname,cname(iname),cform(iname),'diff6m',idiag_diff6m)
        call parse_name(iname,cname(iname),cform(iname),'diff7m',idiag_diff7m)
        call parse_name(iname,cname(iname),cform(iname),'diff8m',idiag_diff8m)
        call parse_name(iname,cname(iname),cform(iname),'diff9m',idiag_diff9m)

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
        write(3,*) 'i_lambdam=',idiag_lambdam
        write(3,*) 'i_num=',idiag_num
        write(3,*) 'i_diff1m=',idiag_diff1m
        write(3,*) 'i_diff2m=',idiag_diff2m
        write(3,*) 'i_diff3m=',idiag_diff3m
        write(3,*) 'i_diff4m=',idiag_diff4m
        write(3,*) 'i_diff5m=',idiag_diff5m
        write(3,*) 'i_diff6m=',idiag_diff6m
        write(3,*) 'i_diff7m=',idiag_diff7m
        write(3,*) 'i_diff8m=',idiag_diff8m
        write(3,*) 'i_diff9m=',idiag_diff9m
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
      character (len=*), intent(in) :: element_name
      real, intent(out) :: MolMass
!
      select case (element_name)
      case ('H')
        MolMass=1.00794
      case ('C')
        MolMass=12.0107
      case ('N')
        MolMass=14.00674
      case ('O')
        MolMass=15.9994
      case ('Ar','AR')
        MolMass=39.948
      case ('He','HE')
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
      integer :: In1,In2,In3,In4,In5,iElement,iTemperature,StopInd
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
      logical :: IsReaction=.false.,found_new_reaction=.false.
      logical, save :: find_specie, found_specie
      integer, optional :: NrOfReactions
      integer, save :: ind_glob, ind_chem
      integer :: i,k,file_id=123, StartInd, StopInd, StartInd_add
      integer :: StopInd_add, StopInd_add_,StopIndName
      integer :: VarNumber, VarNumber_add, SeparatorInd
      integer :: PlusInd
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
                      Mplus_case (k)=.true.
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
        print*,'Write pc_constants.pro in chemistry.f90'
        open (143,FILE=trim(datadir)//'/pc_constants.pro',POSITION="append")
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
        if (maxval(abs(troe_coeff(:,reac))) > 0.) then
          write(file_id,*) 'TROE/',troe_coeff(:,reac)
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
      integer :: k , reac, j, i
      real  :: sum_tmp=0., T_low, T_mid, T_up, tmp, ddd, T_loc
      logical,save :: lwrite=.true.
      character (len=20) :: input_file="./data/react.out"
      integer :: file_id=123
      real :: B_n_0,alpha_n_0,E_an_0
      real, dimension (nx) :: kf_0,Pr,sum_sp
      real, dimension (nx) :: Fcent, ccc, nnn, lnPr, FF,tmpF
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
        T_low=species_constants(k,iTemp1)-10.
        T_mid=species_constants(k,iTemp2)
        T_up= species_constants(k,iTemp3)
         do i=1,nx
          T_loc=p%TT(i)

         ! if (T_loc<T_low) T_loc=T_low
         ! if (T_loc>T_up) T_loc=T_up

          if (T_loc <= T_mid .and. T_low <= T_loc) then
            tmp=0.
            do j=2,5
              tmp=tmp+species_constants(k,iaa2(j))*p%TT(i)**(j-1)/(j-1)
            enddo
            S0_R(i,k)=species_constants(k,iaa2(1))*p%lnTT(i)+tmp&
                +species_constants(k,iaa2(7))
          elseif (T_mid <= T_loc .and. T_loc <= T_up) then
            tmp=0.
            do j=2,5
              tmp=tmp+species_constants(k,iaa1(j))*p%TT(i)**(j-1)/(j-1)
            enddo
            S0_R(i,k)=species_constants(k,iaa1(1))*p%lnTT(i)+tmp&
                +species_constants(k,iaa1(7))
          else
                  if (i==50)  print*,'p%TT(i)=',p%TT(i)
                  !  print*,'i=',i
                    call fatal_error('get_reaction_rate',&
                        'p%TT(i) is outside range')
                    
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
!  Find the product of the species molar consentrations (where
!  each molar consentration is taken to the power of the number)
!
        prod1=1.
        prod2=1.
        do k=1,nchemspec
         if (abs(Sijp(k,reac))>0) then
          prod1=prod1*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)&
              /species_constants(k,imass))**Sijp(k,reac)
         endif
        enddo
        do k=1,nchemspec
         if (abs(Sijm(k,reac))>0) then
          prod2=prod2*(f(l1:l2,m,n,ichemspec(k))*rho_cgs(:)&
              /species_constants(k,imass))**Sijm(k,reac)
         endif
        enddo

      do i=1,nx
!
!  Find forward rate constant for reaction 'reac'
!
    !    if (prod1(i)==0. .and. prod2(i)==0.) then
    !     kf(i)=0.
    !    else
         kf(i)=B_n(reac)*p%TT(i)**alpha_n(reac)*exp(-E_an(reac)/Rcal*p%TT1(i))
    !    endif
!  
!  Find backward rate constant for reaction 'reac'
!

       if (prod2(i)>0) then
        dSR(i)=0.
        dHRT(i)=0.
        sum_tmp=0.
        do k=1,nchemspec
          dSR(i) =dSR(i)+(Sijm(k,reac) -Sijp(k,reac))*S0_R(i,k)
          dHRT(i)=dHRT(i)+(Sijm(k,reac)-Sijp(k,reac))*H0_RT(l1+i-1,m,n,k)
          sum_tmp=sum_tmp+(Sijm(k,reac)-Sijp(k,reac))
        enddo
     
        Kp(i)=exp(dSR(i)-dHRT(i))
        if (sum_tmp==0.) then
          Kc(i)=Kp(i)
        else
          Kc(i)=Kp(i)*(p_atm(i)/(p%TT(i)*Rgas))**sum_tmp
        endif
        if (Kc(i)==0.) then 
          print*,'Kc(i)=',Kc(i),'i=',i,'dSR(i)-dHRT(i)=',dSR(i)-dHRT(i)
          call fatal_error('get_reaction_rate',&
                        'Kc(i)=0')
        else
         kr(i)=kf(i)/Kc(i)
        endif
       endif
      enddo
     
 
        if (lwrite) write(file_id,*) 'Nreact= ',reac,'dSR= ', maxval(dSR)
        if (lwrite) write(file_id,*) 'Nreact= ',reac,'dHRT= ', maxval(dHRT)
        if (lwrite) write(file_id,*)  'Nreact= ',reac,  'kf= ', maxval(kf)
        if (lwrite) write(file_id,*)  'Nreact= ',reac,  'Kc= ', maxval(Kc)
        if (lwrite) write(file_id,*)  'Nreact= ',reac,  'kr= ', maxval(kr)
        if (lwrite) write(file_id,*)'**************************'

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
     !   Kc_0=Kc
        if (maxval(abs(low_coeff(:,reac))) > 0.) then
          B_n_0=low_coeff(1,reac)
          alpha_n_0=low_coeff(2,reac)
          E_an_0=low_coeff(3,reac)
          kf_0(:)=B_n_0*p%TT(:)**alpha_n_0*exp(-E_an_0/Rcal*p%TT1(:))
          Pr=kf_0/kf*mix_conc
          kf=kf*(Pr/(1.+Pr))
    !      kr(:)=kf(:)/Kc_0

        elseif (maxval(abs(high_coeff(:,reac))) > 0.) then
          B_n_0=high_coeff(1,reac)
          alpha_n_0=high_coeff(2,reac)
          E_an_0=high_coeff(3,reac)
          kf_0(:)=B_n_0*p%TT(:)**alpha_n_0*exp(-E_an_0/Rcal*p%TT1(:))
          Pr=kf_0/kf*mix_conc
          kf=kf*(1./(1.+Pr))
    !      kr(:)=kf(:)/Kc_0
        endif

        if (maxval(abs(troe_coeff(:,reac))) > 0.) then
         Fcent=(1.-troe_coeff(1,reac))*exp(-p%TT(:)/troe_coeff(2,reac)) &
         +troe_coeff(1,reac)*exp(-p%TT(:)/troe_coeff(3,reac))
         ccc=-0.4-0.67*log10(Fcent)
         nnn=0.75-1.27*log10(Fcent)
         ddd=0.14
         lnPr=log10(Pr)
         tmpF=((lnPr+ccc)/(nnn-ddd*(lnPr+ccc)))**2
         tmpF=1./(1.+tmpF)
         FF=tmpF*log10(Fcent)
         FF=10**(FF)
         kf=kf*FF
    !     kr(:)=kf(:)/Kc_0
        endif



!
!  Find forward (vreact_p) and backward (vreact_m) rate of
!  progress variable.
!  (vreact_p - vreact_m) is labeled q in the chemkin manual
!
       do i=1,nx  
       if (prod1(i)>0) then
        if (Mplus_case (reac)) then
          vreact_p(i,reac)=prod1(i)*kf(i)
        else
          vreact_p(i,reac)=prod1(i)*kf(i)*sum_sp(i)
        endif
       else
          vreact_p(i,reac)=0.
       endif
       if (prod2(i)>0) then
        kr(i)=kf(i)/Kc(i)
        if (Mplus_case (reac)) then
         vreact_m(i,reac)=prod2(i)*kr(i)
        else
         vreact_m(i,reac)=prod2(i)*kr(i)*sum_sp(i)
        endif
       else
         vreact_m(i,reac)=0.
       endif
      enddo
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
      real :: alpha
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx,mreactions) :: vreactions,vreactions_p,vreactions_m
      real, dimension (nx) :: xdot
      type (pencil_case) :: p
      integer :: k,j
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
          if (lkreactions_alpha) then
            alpha=kreactions_alpha(j)
          else
            alpha=1.
          endif
          vreactions_p(:,j)=alpha*kreactions_p(j)*kreactions_z(n,j)
          vreactions_m(:,j)=alpha*kreactions_m(j)*kreactions_z(n,j)
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
      logical :: emptyfile
      logical :: found_specie
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
      character (len=*), intent(in) :: omega
      real,  dimension(mx,my,mz), intent(in)  :: lnTst
      real,  dimension(mx,my,mz), intent(out) :: Omega_kl
      integer :: i
      real, dimension(8) :: aa
!
      select case (omega)
      case ('Omega11')
        aa(1)= 6.96945701E-1
        aa(2)= 3.39628861E-1
        aa(3)= 1.32575555E-2
        aa(4)=-3.41509659E-2
        aa(5)= 7.71359429E-3
        aa(6)= 6.16106168E-4
        aa(7)=-3.27101257E-4
        aa(8)= 2.51567029E-5
      case ('Omega22')
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
!  Calculation of the binary diffusion coefficients and the species viscosities.
!  This routind is called from calc_for_chem_mixture,
!  which is why we work on full chunks of arrays here.
!
!  28.08.2009: Nils Erland L. Haugen (Corrected the calculation of the 
!              prefactor.) 
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
      real, dimension (mx,my,mz) :: Omega_kl, prefactor
      real, dimension (mx,my,mz) :: lnTjk,lnTk_array
      integer :: k,j,j1,j2,j3
      real :: eps_jk, sigma_jk, m_jk, delta_jk, delta_st
      character (len=7) :: omega
      real :: Na=6.022E23 ,tmp_local,tmp_local2, delta_jk_star


!
!  Find binary diffusion coefficients
!
      tmp_local=3./16.*(2.*k_B_cgs**3/pi)**0.5
      do j3=nn1,nn2
      do j2=mm1,mm2
      do j1=ll1,ll2
        prefactor(j1,j2,j3)=tmp_local*(TT_full(j1,j2,j3))**0.5*unit_length**3&
            /(Rgas_unit_sys*rho_full(j1,j2,j3))
      enddo
      enddo
      enddo
! NILS: remove omega?
      omega="Omega11"
!
! Check if we use fixed Schmidt number for speeding up calculations
!
      if (.not. lfix_Sc) then
!
! Check if we use simplified version of the binary diffusion calculation
!
        if (BinDif_simple) then
          call stop_it('BinDif_simple case does not work now!')
        elseif (ldiffusion) then
!
!  Do non-simplified binary diffusion coefficient
!
          do k=1,nchemspec
            do j=k,nchemspec
       
!  Account for the difference between eq. 5-4 and 5-31 in the Chemkin theory
!  manual 
!
            if (j/=k) then
                eps_jk=(tran_data(j,2)*tran_data(k,2))**0.5
                sigma_jk=0.5*(tran_data(j,3)+tran_data(k,3))*1e-8
                m_jk=(species_constants(j,imass)*species_constants(k,imass)) &
                    /(species_constants(j,imass)+species_constants(k,imass))/Na
                delta_jk=0.5*tran_data(j,4)*tran_data(k,4)*1e-18*1e-18
              else
                eps_jk=tran_data(j,2)
                sigma_jk=tran_data(j,3)*1e-8
                m_jk=species_constants(j,imass)/(2*Na)
                delta_jk=0.5*(tran_data(j,4)*1e-18)**2
              endif
  
!
! Loop over all grid points
!
              do j3=nn1,nn2
              do j2=mm1,mm2
              do j1=ll1,ll2
                if (ltemperature_nolog) then
                  lnTjk(j1,j2,j3)=log(f(j1,j2,j3,ilnTT)/eps_jk)
                else
                  lnTjk(j1,j2,j3)=f(j1,j2,j3,ilnTT)-log(eps_jk)
                endif
!
                Omega_kl(j1,j2,j3)= &
                    1./(6.96945701E-1   +3.39628861E-1*lnTjk(j1,j2,j3) &
                    +1.32575555E-2*lnTjk(j1,j2,j3)*lnTjk(j1,j2,j3) &
                    -3.41509659E-2*lnTjk(j1,j2,j3)**3 &
                    +7.71359429E-3*lnTjk(j1,j2,j3)**4 &
                    +6.16106168E-4*lnTjk(j1,j2,j3)**5 &
                    -3.27101257E-4*lnTjk(j1,j2,j3)**6 &
                    +2.51567029E-5*lnTjk(j1,j2,j3)**7)
                delta_jk_star=delta_jk/(eps_jk*k_B_cgs*sigma_jk**3)
                
                Omega_kl(j1,j2,j3)=Omega_kl(j1,j2,j3)&
                    +0.19*delta_jk_star*delta_jk_star/(TT_full(j1,j2,j3)/eps_jk)
               if (j/=k) then
                Bin_Diff_coef(j1,j2,j3,k,j)=prefactor(j1,j2,j3)/mu1_full(j1,j2,j3)&
                    /(sqrt(m_jk)*sigma_jk**2*Omega_kl(j1,j2,j3))


               else
                Bin_Diff_coef(j1,j2,j3,k,j)=prefactor(j1,j2,j3)&
                    /(sqrt(m_jk)*sigma_jk**2*Omega_kl(j1,j2,j3))*species_constants(k,imass)
 
               endif
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
       Bin_Diff_coef=0.
      endif
      endif

!
!  Calculate viscosity
!
      omega="Omega22"

      tmp_local=5./16.*(k_B_cgs/(Na*pi))**0.5

      do k=1,nchemspec
        tmp_local2=(species_constants(k,imass))**0.5/(tran_data(k,3)*1e-8)**2 &
            *tmp_local
        if (visc_simple) then
          call stop_it('visc_simple case does not work now!')
        elseif (visc_const==impossible) then
!
! 1 Debye = 10**(-18) esu -> (1e-18*tran_data(k,4))
!
          delta_st=(1e-18*tran_data(k,4))**2/2./ &
             (tran_data(k,2)*k_B_cgs*(tran_data(k,3)*1e-8)**3)
!
          if (ltemperature_nolog) then
            lnTk_array=log(f(:,:,:,ilnTT)/tran_data(k,2))
          else
            lnTk_array=f(:,:,:,ilnTT)-log(tran_data(k,2))
          endif
          call calc_collision_integral(omega,lnTk_array,Omega_kl)


          do j3=nn1,nn2
          do j2=mm1,mm2
          do j1=ll1,ll2
           species_viscosity(j1,j2,j3,k)=TT_full(j1,j2,j3)**0.5&
               /(Omega_kl(j1,j2,j3) &
               +0.2*delta_st*delta_st/(TT_full(j1,j2,j3)/tran_data(k,2)))*tmp_local2 &
               /(unit_mass/unit_length/unit_time)
          enddo
          enddo
          enddo
        endif
      enddo
     !
 
    endsubroutine calc_diff_visc_coef
!***************************************************************
    subroutine calc_therm_diffus_coef(f)
! 
!  Calculate the thermal diffusion coefficient based on equation 5-17 in
!  the Chemkin theory manual
!     
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,nchemspec) :: species_cond
      real, dimension (mx,my,mz) :: tmp_sum, tmp_sum2, tmp_val,ZZ,FF
      real, dimension (mx,my,mz) :: AA,BB, Cv_vib_R, f_tran, f_rot, f_vib
      intent(in) :: f
      integer :: j1,j2,j3,k
      real :: Cv_rot_R, Cv_tran_R,T_st, pi_1_5, pi_2

      pi_1_5=pi**1.5
      pi_2=pi**2.
!
      do j3=nn1,nn2
      do j2=mm1,mm2
      do j1=1,mx
        tmp_sum(j1,j2,j3)=0.
        tmp_sum2(j1,j2,j3)=0.
        do k=1,nchemspec
!
! Check if the molecule is a single atom (0), linear (1) or non-linear (2).
!
          if (tran_data(k,1)==0.) then
            Cv_tran_R=1.5
            Cv_rot_R=0.
            Cv_vib_R=0.
          elseif (tran_data(k,1)==1.) then
            Cv_tran_R=1.5
            Cv_rot_R=1.
            Cv_vib_R(j1,j2,j3)=cv_R_spec_full(j1,j2,j3,k)-2.5
          elseif (tran_data(k,1)==2.) then
            Cv_tran_R=1.5
            Cv_rot_R=1.5
            Cv_vib_R(j1,j2,j3)=cv_R_spec_full(j1,j2,j3,k)-3.
          endif
!
! The rotational and vibrational contributions are zero for the single
! atom molecules but not for the linear or non-linear molecules 
!

          if (tran_data(k,1)>0.) then
            tmp_val(j1,j2,j3)=Bin_Diff_coef(j1,j2,j3,k,k)*rho_full(j1,j2,j3)&
                /species_viscosity(j1,j2,j3,k)
            AA(j1,j2,j3)=2.5-tmp_val(j1,j2,j3)
            T_st=tran_data(k,2)/298.
            FF(j1,j2,j3)=1.+pi_1_5/2.*(T_st)**0.5+(pi_2/4.+2.) &
             *(T_st)+pi_1_5*(T_st)**1.5
            ZZ(j1,j2,j3)=tran_data(k,6)*FF(j1,j2,j3)
            T_st=tran_data(k,2)/TT_full(j1,j2,j3)
            FF(j1,j2,j3)=1.+pi_1_5/2.*(T_st)**0.5+(pi_2/4.+2.) &
             *(T_st)+pi_1_5*(T_st)**1.5
             ZZ(j1,j2,j3)=ZZ(j1,j2,j3)/FF(j1,j2,j3)
            BB(j1,j2,j3)=ZZ(j1,j2,j3)+2./pi*(5./3.*Cv_rot_R&
                +tmp_val(j1,j2,j3))
            f_tran(j1,j2,j3)=2.5*(1.- 2./pi*Cv_rot_R/&
                Cv_tran_R*AA(j1,j2,j3)/BB(j1,j2,j3))
            f_rot(j1,j2,j3)=tmp_val(j1,j2,j3)&
                *(1+2./pi*AA(j1,j2,j3)/BB(j1,j2,j3))
            f_vib(j1,j2,j3)=tmp_val(j1,j2,j3)
          else
            f_tran(j1,j2,j3)=2.5     
            f_rot(j1,j2,j3) =0.0   
            f_vib(j1,j2,j3) =0.0   
          endif
          species_cond(j1,j2,j3,k)=(species_viscosity(j1,j2,j3,k)) &
              /(species_constants(k,imass)/unit_mass)*Rgas* &
                (f_tran(j1,j2,j3)*Cv_tran_R+f_rot(j1,j2,j3)*Cv_rot_R  &
                +f_vib(j1,j2,j3)*Cv_vib_R(j1,j2,j3))
 
!
! tmp_sum and tmp_sum2 are used later to find the mixture averaged
! conductivity.
!
          tmp_sum(j1,j2,j3)=tmp_sum(j1,j2,j3)  &
              +XX_full(j1,j2,j3,k)*species_cond(j1,j2,j3,k)
          tmp_sum2(j1,j2,j3)=tmp_sum2(j1,j2,j3) &
              +XX_full(j1,j2,j3,k)/species_cond(j1,j2,j3,k)
          enddo
        enddo
        enddo
        enddo
!
! Find the mixture averaged conductivity
!
        do j1=1,mx
        do j2=mm1,mm2
        do j3=nn1,nn2
          if ((tmp_sum2(j1,j2,j3))<=0.) then
            lambda_full(j1,j2,j3)=0.
          else
            lambda_full(j1,j2,j3)=0.5*(tmp_sum(j1,j2,j3)+1.&
                /tmp_sum2(j1,j2,j3))
         !   lambda_full(j1,j2,j3)=(tmp_sum(j1,j2,j3))
          endif
          if (lambda_const<impossible) then
            lambda_full(j1,j2,j3)=lambda_const
          endif
        enddo
        enddo
        enddo
!
        call keep_compiler_quiet(f)
!
    endsubroutine calc_therm_diffus_coef
!***************************************************************
    subroutine calc_diffusion_term(f,p)
!
!  Calculate diffusion term, p%DYDt_diff
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz) :: Xk_Yk
      real, dimension (nx,3) :: gDiff_full_add, gchemspec, gXk_Yk
      real, dimension (nx) :: del2chemspec
      real, dimension (nx) :: diff_op,diff_op1,diff_op2,del2XX, del2pp, del2lnpp
      real, dimension (nx) :: glnpp_gXkYk,glnrho_glnpp,gD_glnpp, glnpp_glnpp
      real :: diff_k
      integer :: k
!
      intent(in) :: f
!
      p%DYDt_diff=0.
      diff_k=chem_diff
!
!  Loop over all chemical species.
!
      do k=1,nchemspec
!
!  Simple chemistry (not using ChemKin formalism).
!  Eliminate need for p%glnrho if density is not included.
!
        if (.not. lcheminp) then
          if (chem_diff/=0.) then
            diff_k=chem_diff*chem_diff_prefactor(k)
            if (headtt) print*,'dchemistry_dt: k,diff_k=',k,diff_k
            call del2(f,ichemspec(k),del2chemspec)
            call grad(f,ichemspec(k),gchemspec)
            if (ldensity) then
              call dot_mn(p%glnrho,gchemspec,diff_op)
              diff_op=diff_op+del2chemspec
            else
              diff_op=del2chemspec
            endif
            p%DYDt_diff(:,k)=diff_k*diff_op
          endif
        else
!
          if (ldiffusion) then
!
          
            call del2(XX_full(:,:,:,k),del2XX)
         
            call dot_mn(p%glnrho,p%gXXk(:,:,k),diff_op1)
            call grad(Diff_full_add(:,:,:,k),gDiff_full_add)
            call dot_mn(gDiff_full_add,p%gXXk(:,:,k),diff_op2)
!

       
            call del2(pp_full(:,:,:),del2pp)
            call dot_mn(p%glnpp,p%glnpp,glnpp_glnpp)

            del2lnpp=del2pp/p%pp-glnpp_glnpp

            Xk_Yk=XX_full(:,:,:,k)-f(:,:,:,ichemspec(k))
            call grad(Xk_Yk,gXk_Yk)
            
           

            call dot_mn(p%glnrho,p%glnpp,glnrho_glnpp)
            call dot_mn(gDiff_full_add,p%glnpp,gD_glnpp)
            call dot_mn(gXk_Yk,p%glnpp,glnpp_gXkYk)


          endif
          p%DYDt_diff(:,k)=Diff_full_add(l1:l2,m,n,k)*(del2XX+diff_op1)+diff_op2 &
          +Diff_full_add(l1:l2,m,n,k)*Xk_Yk(l1:l2,m,n)*del2lnpp &
          +Diff_full_add(l1:l2,m,n,k)*Xk_Yk(l1:l2,m,n)*glnrho_glnpp &
          +Xk_Yk(l1:l2,m,n)*gD_glnpp+Diff_full_add(l1:l2,m,n,k)*glnpp_gXkYk
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

      real, dimension(nx) :: g2TT, g2TTlambda=0., tmp1
!

     if (ltemperature_nolog) then
      call dot(p%gTT,p%glambda,g2TTlambda)
     else
      call dot(p%glnTT,p%glambda,g2TTlambda)
      call dot(p%glnTT,p%glnTT,g2TT)
     endif
!
!  Add heat conduction to RHS of temperature equation
!

    !  if (l1step_test) then
    !   tmp1= p%lambda(:)*(p%del2lnTT+g2TT)*p%cv1/p%rho(:)
    !  else
    !   tmp1= (p%lambda(:)*(p%del2lnTT+g2TT)+g2TTlnlambda)*p%cv1/p%rho(:)
    !  endif
      if (l1step_test) then
       tmp1= p%lambda(:)*(p%del2lnTT+g2TT)*p%cv1/p%rho(:)
      else
      if (ltemperature_nolog) then
       tmp1= (p%lambda(:)*p%del2lnTT+g2TTlambda)*p%cv1/p%rho(:)
       df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + tmp1
      else
       tmp1= (p%lambda(:)*(p%del2lnTT+g2TT)+g2TTlambda)*p%cv1/p%rho(:)
       df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + tmp1
      endif
      endif
       

    !  df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + tmp1
!
      RHS_T_full(l1:l2,m,n)=RHS_T_full(l1:l2,m,n) + tmp1
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_heatcond_chemistry
!***********************************************************************
   subroutine get_cs2_full(cs2_full)
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
!***********************************************************************
   subroutine get_cs2_slice(slice,dir,index)
!
! Find a slice of the speed of sound
!
! 2009.12.10: Nils Erland L. Haugen (coded)
!
      real, dimension (:,:), intent(out) :: slice
      integer, intent(in) :: index, dir
!
      if (dir==1) then 
        slice=cp_full(index,:,:)/cv_full(index,:,:)&
            *mu1_full(index,:,:)*TT_full(index,:,:)*Rgas
      elseif (dir==2) then 
        slice=cp_full(:,index,:)/cv_full(:,index,:)&
            *mu1_full(:,index,:)*TT_full(:,index,:)*Rgas
      elseif (dir==3) then 
        slice=cp_full(:,:,index)/cv_full(:,:,index)&
            *mu1_full(:,:,index)*TT_full(:,:,index)*Rgas
      else
        call fatal_error('get_cs2_slice','No such dir!')
      endif
!
    endsubroutine get_cs2_slice
!*************************************************************
    subroutine get_gamma_full(gamma_full)
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
    subroutine get_gamma_slice(slice,dir,index)
!
!  Get a 2D slice of gamma
!
!  2009.12.10: Nils Erland L. Haugen (coded)
!
      real, dimension (:,:), intent(out)  :: slice
      integer, intent(in) :: index, dir
!
      if (dir==1) then 
        slice=cp_full(index,:,:)/cv_full(index,:,:)
      elseif (dir==2) then 
        slice=cp_full(:,index,:)/cv_full(:,index,:)
      elseif (dir==3) then 
        slice=cp_full(:,:,index)/cv_full(:,:,index)
      else
        call fatal_error('get_gamma_slice','No such dir!')
      endif
      !
    endsubroutine get_gamma_slice
!*************************************************************
    subroutine air_field(f)
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: sum_Y
!
      logical :: emptyfile=.true.
      logical :: found_specie
      integer :: file_id=123, ind_glob, ind_chem
      character (len=80) :: ChemInpLine
      character (len=10) :: specie_string
      character (len=1)  :: tmp_string
      integer :: j,k=1
      real :: YY_k, air_mass, TT=300., PP=1.013e6 ! (in dynes = 1atm)
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


     
        if (.not. lflame_front) then

      if (ltemperature_nolog) then
       f(:,:,:,iTT)=TT
      else
       f(:,:,:,ilnTT)=log(TT)!+f(:,:,:,ilnTT)
      endif
!
      if (ldensity_nolog) then
       f(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)*&
          air_mass/TT)/unit_mass*unit_length**3
      else
       f(:,:,:,ilnrho)=log((PP/(k_B_cgs/m_u_cgs)*&
          air_mass/TT)/unit_mass*unit_length**3)
      endif

      endif


   !   f(:,:,:,iux)=init_ux

!
      if (lroot) print*, 'Air temperature, K', TT
      if (lroot) print*, 'Air pressure, dyn', PP
      if (lroot) print*, 'Air density, g/cm^3:'
      if (lroot) print '(E10.3)',  PP/(k_B_cgs/m_u_cgs)*air_mass/TT
      if (lroot) print*, 'Air mean weight, g/mol', air_mass
      if (lroot) print*, 'R', k_B_cgs/m_u_cgs
!
      close(file_id)
!
    endsubroutine air_field
!********************************************************************
!!!!!!!!!  NSCBC boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine bc_nscbc_subin_x_new(f,df,topbot,val)
!
!   nscbc case
!   subsonic inflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    ux,uy,T are fixed, drho  is calculated
!
!   16-nov-09/natalia: coded.
!
!
      use Deriv, only: der_onesided_4_slice,der_pencil, der2_pencil
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (my,mz) :: rho0,TT0,gamma0
      real, dimension (mx,my,mz) :: mom2,mom3!, rho_ux2, rho_uy2
      real, dimension (mx,my,mz) :: rhoE_p
      real, dimension (mx,my,mz,2) ::  rhoE_pU
      real, dimension (my,mz) ::  dYk_dx,dYk_dy,dYk_dz
      real, dimension (ny,nz) :: drho_prefac, KK
      
      real, dimension (ny,nz) :: L_1, L_2, L_3, L_4, L_5
      real, dimension (ny,nz)  :: M_1, M_2, M_3, M_4, M_5
  !    
      real, dimension (my,mz) :: cs0_ar,cs20_ar,dmom2_dy,dmom3_dz
      real, dimension (my,mz,2) :: drhoE_pU!,dYk_dy
   !   real, dimension (mx,my,mz,nchemspec) :: bound_rhs_Y
    !  real, dim(:,nnn)ension (ny,nz) :: bound_rhs_T
      real, dimension (mx,my,mz) :: cs2_full, gamma_full
      real, dimension (nx,ny,nz) :: p_inf
      
      real, dimension (ny,nz,3,3) :: dui_dxj
      real, dimension (my,mz)     :: tmp1,tmp2,tmp3, tmp_rho,tmp_pp       
      real, dimension (ny,nz,3)   :: grad_rho, grad_pp
     ! real, dimension(ny,nz) ::     d2u1_dy2,d2u1_dz2
     ! real, dimension(ny,nz) ::     d2u2_dy2,d2u2_dz2,d2u3_dy2,d2u3_dz2
     ! real, dimension (ny,nz,3) :: dlnT_dxj
      real, dimension (ny,nz) :: T_1_y, T_2_y, T_3_y, T_4_y, T_5_y
      real, dimension (ny,nz) :: T_1_z, T_2_z, T_3_z, T_4_z, T_5_z
      real, dimension (mcom), optional :: val
!
      integer :: lll, sgn,i,j,k, mm, mmm, irho_tmp
      real :: Mach_num
      real :: U0_x,U0_y,U0_z,T0
!
     
      intent(inout) :: f
      intent(out) :: df
      real :: nscbc_sigma=0.5,ita
!
  
      if (present(val)) call keep_compiler_quiet(val)

       U0_x=init_ux;U0_y=init_uy;U0_z=init_uz;T0=init_TT1

      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case (topbot)
      case ('bot')
        lll = l1; sgn = 1
        if (leos_chemistry) then
          p_inf(1,:,:)=pp_full(l1,m1:m2,n1:n2)
        endif
      case ('top')
        lll = l2; sgn = -1
        if (leos_chemistry) then
          p_inf(nx,:,:)=pp_full(l2,m1:m2,n1:n2)
        endif
      case default
        print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect
 

      if (leos_chemistry) then
         cs20_ar=cs2_full(lll,:,:)
         cs0_ar=cs2_full(lll,:,:)**0.5
         gamma0=gamma_full(lll,:,:)

        if (ldensity_nolog) then
         rho_full = f(:,:,:,irho)
          rho0(:,:) = f(lll,:,:,irho)
          drho_prefac=-1./cs20_ar(m1:m2,n1:n2)
          irho_tmp=irho
          call stop_it('bc_nscbc_nref_subout_x: NSCBC works now only for lnrho')
        else
          rho_full = exp(f(:,:,:,ilnrho))
          rho0(:,:) = rho_full(lll,:,:)
          drho_prefac=-1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2)
          irho_tmp=ilnrho
        endif

        if (ltemperature_nolog) then
         TT0(:,:) = f(lll,:,:,iTT)
        else
         TT0(:,:) = exp(f(lll,:,:,ilnTT))
        endif
         mom2(lll,:,:)=rho0(:,:)*f(lll,:,:,iuy)
         mom3(lll,:,:)=rho0(:,:)*f(lll,:,:,iuz)
         rhoE_p(lll,:,:)=0.5*rho_full(lll,:,:) &
            *(f(lll,:,:,iux)**2+f(lll,:,:,iuy)**2+f(lll,:,:,iuz)**2) &
             +gamma_full(lll,:,:)/(gamma_full(lll,:,:)-1)*pp_full(lll,:,:)
         rhoE_pU(lll,:,:,1)=rhoE_p(lll,:,:)*f(lll,:,:,iuy) 
         rhoE_pU(lll,:,:,2)=rhoE_p(lll,:,:)*f(lll,:,:,iuz)


         Mach_num=maxval(f(lll,m1:m2,n1:n2,iux)/cs0_ar(m1:m2,n1:n2))
         KK=nscbc_sigma*(1.-Mach_num*Mach_num)/Lxyz(1)

      else
        print*,"bc_nscbc_subin_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented "//&
            "for the chemistry case."
        print*,"Boundary treatment skipped."
        return
      endif

        call der_onesided_4_slice(rho_full,sgn,grad_rho(:,:,1),lll,1)
        call der_onesided_4_slice(pp_full,sgn,grad_pp(:,:,1),lll,1)
        call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,1),lll,1)
        call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,1),lll,1)
        call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,1),lll,1)
       ! call der_onesided_4_slice(f,sgn,ilnTT,dlnT_dxj(:,:,1),lll,1)


       if (nygrid /= 1) then
    
         do i=n1,n2
          call der_pencil(2,f(lll,:,i,iux),tmp1(:,i))
          call der_pencil(2,f(lll,:,i,iuy),tmp2(:,i))
          call der_pencil(2,f(lll,:,i,iuz),tmp3(:,i))
          call der_pencil(2,rho_full(lll,:,i),tmp_rho(:,i))
          call der_pencil(2,pp_full(lll,:,i),tmp_pp(:,i))
          call der_pencil(2,mom2(lll,:,i),dmom2_dy(:,i))
          call der_pencil(2,rhoE_pU(lll,:,i,1),drhoE_pU(:,i,1))
         enddo

      else
        tmp3=0
        tmp2=0
        tmp1=0
        tmp_rho=0
        tmp_pp=0
        dmom2_dy=0
        drhoE_pU(:,:,1)=0

      endif
        dui_dxj(:,:,1,2)=tmp1(m1:m2,n1:n2)
        dui_dxj(:,:,2,2)=tmp2(m1:m2,n1:n2)
        dui_dxj(:,:,3,2)=tmp3(m1:m2,n1:n2)
        grad_rho(:,:,2)=tmp_rho(m1:m2,n1:n2)
        grad_pp(:,:,2)=tmp_pp(m1:m2,n1:n2)

      if (nzgrid /= 1) then
         do j=m1,m2
          call der_pencil(3,f(lll,j,:,iux),tmp1(j,:))
          call der_pencil(3,f(lll,j,:,iuy),tmp2(j,:))
          call der_pencil(3,f(lll,j,:,iuz),tmp3(j,:))
          call der_pencil(3,rho_full(lll,j,:),tmp_rho(j,:))
          call der_pencil(3,pp_full(lll,j,:),tmp_pp(j,:))
          call der_pencil(3,mom3(lll,j,:),dmom3_dz(j,:))
          call der_pencil(3,rhoE_pU(lll,j,:,2),drhoE_pU(j,:,2))
         enddo
      else
        tmp3=0
        tmp2=0
        tmp1=0
        tmp_rho=0
        tmp_pp=0
        dmom3_dz=0
        drhoE_pU(:,:,2)=0
      endif
        dui_dxj(:,:,1,3)=tmp1(m1:m2,n1:n2)
        dui_dxj(:,:,2,3)=tmp2(m1:m2,n1:n2)
        dui_dxj(:,:,3,3)=tmp3(m1:m2,n1:n2)
        grad_rho(:,:,3)=tmp_rho(m1:m2,n1:n2)
        grad_pp(:,:,3)=tmp_pp(m1:m2,n1:n2)

!


       if ((nygrid /= 1) .or.(nzgrid /= 1))  then

          T_1_y(:,:)=-dmom2_dy(m1:m2,n1:n2)/rho0(m1:m2,n1:n2)
          T_1_z(:,:)=-dmom3_dz(m1:m2,n1:n2)/rho0(m1:m2,n1:n2)
          T_2_y(:,:)=-f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,1,2) 
          T_2_z(:,:)=-f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,1,3) !&
                !   -nu_full(lll,m1:m2,n1:n2)*(d2u1_dy2+d2u1_dz2)
          T_3_y(:,:)=-f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,2,2) &
                     -grad_pp(:,:,2)/rho0(m1:m2,n1:n2) 
         T_3_z(:,:)=-f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,2,3) 
                !   -nu_full(lll,m1:m2,n1:n2)*(d2u2_dy2+d2u2_dz2)
          T_4_y(:,:)=-f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,3,2) 
         T_4_z(:,:)=-f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,3,3) &
                     -grad_pp(:,:,3)/rho0(m1:m2,n1:n2)
          T_5_y(:,:)= drho_prefac(:,:)*(gamma0(m1:m2,n1:n2)-1.) &
                      *gamma0(m1:m2,n1:n2)*drhoE_pU(m1:m2,n1:n2,1)
          T_5_z(:,:)= drho_prefac(:,:)*(gamma0(m1:m2,n1:n2)-1.) &
                      *gamma0(m1:m2,n1:n2)*(drhoE_pU(m1:m2,n1:n2,2))

       endif


      ita=1.

      select case (topbot)
      case ('bot')
        L_5=ita*(1.-Mach_num*Mach_num)/Lxyz(1)*cs20_ar(m1:m2,n1:n2)*rho0(m1:m2,n1:n2)&
           *(f(lll,m1:m2,n1:n2,iux)-U0_x) &
           +(T_5_y +T_5_z) &
           +rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*(T_2_y +T_2_z)
      case ('top')
        L_1=ita*(1.-Mach_num*Mach_num)/Lxyz(1)*cs20_ar(m1:m2,n1:n2)*rho0(m1:m2,n1:n2)&
           *(f(lll,m1:m2,n1:n2,iux)-U0_x)&
           +(T_5_y +T_5_z) &
           -rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*(T_2_y +T_2_z)
     endselect
!
     

       L_2 =ita*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*Rgas/Lxyz(1) &
           *(TT0(m1:m2,n1:n2)-T0)+cs20_ar(m1:m2,n1:n2)*(T_1_y +T_1_z)-(T_5_y +T_5_z)
       L_3 = ita*cs0_ar(m1:m2,n1:n2)/Lxyz(1)*(f(lll,m1:m2,n1:n2,iuy)-U0_y) &
             +(T_3_y +T_3_z)
       L_4 = ita*cs0_ar(m1:m2,n1:n2)/Lxyz(1)*(f(lll,m1:m2,n1:n2,iuz)-U0_z) &
             +(T_4_y +T_4_z)


        if ((nygrid /= 1) .or.(nzgrid /= 1))  then

        M_1 = (f(lll,m1:m2,n1:n2,iuy) - cs0_ar(m1:m2,n1:n2))*&
            (grad_pp(:,:,2)- rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,2,2))

        M_2 = f(lll,m1:m2,n1:n2,iuy)*(cs20_ar(m1:m2,n1:n2)*grad_rho(:,:,2)-grad_pp(:,:,2))
        M_3 = f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,1,2)
        M_4 = f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,3,2)

        M_5 = (f(lll,m1:m2,n1:n2,iuy) + cs0_ar(m1:m2,n1:n2))*&
            (grad_pp(:,:,2)+ rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,2,2))


       do i=1,2

       if (i==1) then
         mm=m1
         mmm=1

       select case (topbot)
         case ('bot')
           M_5(mmm,:)=0.
           L_5(mmm,:)=L_5(mmm,:)-M_1(mmm,:)/2.-rho0(mm,n1:n2)*cs0_ar(mm,n1:n2)*M_2(mmm,:)
           L_3(mmm,:)=L_3(mmm,:)+M_1(mmm,:)/(2.*rho0(mm,n1:n2)*cs0_ar(mm,n1:n2))
         case ('top')
           M_1(mmm,:)=0.
           L_1(mmm,:)=L_1(mmm,:)-M_5(mmm,:)/2.-rho0(mm,n1:n2)*cs0_ar(mm,n1:n2)*M_2(mmm,:)
           L_3(mmm,:)=L_3(mmm,:)-M_5(mmm,:)/(2.*rho0(mm,n1:n2)*cs0_ar(mm,n1:n2))
       endselect
       elseif (i==2) then
         mm=m2
         mmm=ny
        select case (topbot)
         case ('bot')
           M_5(mmm,:)=0.
           L_5(mmm,:)=L_5(mmm,:)-M_1(mmm,:)/2.-rho0(mm,n1:n2)*cs0_ar(mm,n1:n2)*M_2(mmm,:)
           L_3(mmm,:)=L_3(mmm,:)+M_1(mmm,:)/(2.*rho0(mm,n1:n2)*cs0_ar(mm,n1:n2))
          case ('top')
           M_1(mmm,:)=0.
           L_1(mmm,:)=L_1(mmm,:)-M_5(mmm,:)/2.-rho0(mm,n1:n2)*cs0_ar(mm,n1:n2)*M_2(mmm,:)
           L_3(mmm,:)=L_3(mmm,:)-M_5(mmm,:)/(2.*rho0(mm,n1:n2)*cs0_ar(mm,n1:n2))
         endselect
       endif



       L_2(mmm,:) = L_2(mmm,:)-M_3(mmm,:)
       L_4(mmm,:) = L_4(mmm,:)-M_4(mmm,:)


       !  df(lll,m1:m2,n1:n2,irho_tmp) = df(lll,m1:m2,n1:n2,irho_tmp) + T_1_y + T_1_z
       !  df(lll,m1:m2,n1:n2,iux)      = df(lll,m1:m2,n1:n2,iux)      + T_2_y + T_2_z
       !  df(lll,m1:m2,n1:n2,iuy)      = df(lll,m1:m2,n1:n2,iuy)      + T_3_y + T_3_z
       !  df(lll,m1:m2,n1:n2,iuz)      = df(lll,m1:m2,n1:n2,iuz)      + T_4_y + T_4_z
       !  df(lll,m1:m2,n1:n2,ilnTT)    = df(lll,m1:m2,n1:n2,ilnTT)    + T_5_y + T_5_z

       enddo

       endif


!
       df(lll,m1:m2,n1:n2,irho_tmp) = &
         drho_prefac*(L_2+0.5*(L_5 + L_1)) 
       df(lll,m1:m2,n1:n2,iux) =  &
         -1./(2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2))*(L_5 - L_1)
       df(lll,m1:m2,n1:n2,iuy) = -L_3
       df(lll,m1:m2,n1:n2,iuz) = -L_4

       if (ltemperature_nolog) then
        df(lll,m1:m2,n1:n2,iTT) = &
          -1./(rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))*(-L_2 &
          +0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1))*TT0(m1:m2,n1:n2)
       else 
        df(lll,m1:m2,n1:n2,ilnTT) = &
          -1./(rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))*(-L_2 &
          +0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1)) 
       endif



!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nchemspec>1) then
       do k=1,nchemspec
          call der_onesided_4_slice(f,sgn,ichemspec(k),dYk_dx(m1:m2,n1:n2),lll,1)
         do i=n1,n2
          call der_pencil(2,f(lll,:,i,ichemspec(k)),dYk_dy(:,i))
         enddo
         do i=m1,m2
           call der_pencil(3,f(lll,i,:,ichemspec(k)),dYk_dz(i,:))
         enddo
!NMNMNMN
           df(lll,m1:m2,n1:n2,ichemspec(k))=&
              -f(lll,m1:m2,n1:n2,iux)*dYk_dx(m1:m2,n1:n2) &
              -f(lll,m1:m2,n1:n2,iuy)*dYk_dy(m1:m2,n1:n2) &
              -f(lll,m1:m2,n1:n2,iuz)*dYk_dz(m1:m2,n1:n2) &
              + RHS_Y_full(lll,m1:m2,n1:n2,k)

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
    endsubroutine bc_nscbc_subin_x_new
!***********************************************************
   subroutine bc_nscbc_subin_x(f,df,topbot,val)
!
!   nscbc case
!   subsonic inflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    ux,uy,T are fixed, drho  is calculated
!
!   16-nov-08/natalia: coded.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: mom2
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension(ny,nz) :: dux_dx, L_1, L_2, L_5, dpp_dx
      real, dimension(my,mz) :: rho0, gamma0, dmom2_dy, TT0
      real, dimension(my,mz) :: cs0_ar,cs20_ar
    !  real, dimension (my,mz) :: cs2x
      real, dimension (mx,my,mz) :: cs2_full, gamma_full, rho_full, pp
      real, dimension(nchemspec) :: YYi
      real, dimension (mcom), optional :: val
      integer :: lll,i, sgn,k
      real :: u_t, T_t
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
      select case (topbot)
      case ('bot')
        lll = l1
        sgn = 1
      case ('top')
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
       !   TT0=TT_full(lll,:,:)
         if (ltemperature_nolog) then
          TT0=f(lll,:,:,iTT)
         else
          TT0=exp(f(lll,:,:,ilnTT))
         endif
         
         if (ldensity_nolog) then
          rho_full=f(:,:,:,irho)
         else
          rho_full=exp(f(:,:,:,ilnrho))
         endif

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

        select case (topbot)
        case ('bot')
          L_1 = (f(lll,m1:m2,n1:n2,iux) - cs0_ar(m1:m2,n1:n2))*&
              (dpp_dx - rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
          L_5 =L_1-2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*&
              df(lll,m1:m2,n1:n2,iux)
        case ('top')
          L_5 = (f(lll,m1:m2,n1:n2,iux) + cs0_ar(m1:m2,n1:n2))*&
              (dpp_dx + rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dux_dx)
          L_1 = L_5+2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*&
              df(lll,m1:m2,n1:n2,iux)
        endselect

        if (ltemperature_nolog) then
        L_2 = 0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1) &
            +rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2) &
            *df(lll,m1:m2,n1:n2,iTT)/TT0(m1:m2,n1:n2)
        else
        L_2 = 0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1) &
            +rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2)*df(lll,m1:m2,n1:n2,ilnTT)
        endif
        if (ldensity_nolog) then
          df(lll,m1:m2,n1:n2,ilnrho) = -1./cs20_ar(m1:m2,n1:n2)*&
              (L_2+0.5*(L_5 + L_1)) ! -dmom2_dy(m1:m2,n1:n2)
        else
          df(lll,m1:m2,n1:n2,ilnrho) = &
              -1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2) &
              *(L_2+0.5*(L_5 + L_1)) !&
           !-1./rho0(m1:m2,n1:n2)*dmom2_dy(m1:m2,n1:n2)
        endif

!
! natalia: this subroutine is still under construction
! Please, do not remove this commented part !!!!!
!

!
!  this conditions can be important!
!  check without them
!

   !     do k=1,nchemspec
   !      f(lll,m1:m2,n1:n2,ichemspec(k))=YYi(k)
   !     enddo
   !      f(lll,m1:m2,n1:n2,iux) = u_t
   !     f(lll,m1:m2,n1:n2,ilnTT) = T_t
   !      df(lll,m1:m2,n1:n2,ilnTT)=0.
    endsubroutine bc_nscbc_subin_x
!***********************************************************
!***********************************************************
    subroutine bc_nscbc_nref_subout_x(f,df,topbot,nscbc_sigma_out)

!
!   nscbc case
!   subsonic non-reflecting outflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    drho. dT, dux, duy  are calculated, p_inf can be
!   fixed (if nscbc_sigma <>0)
!
!   16-nov-08/natalia: coded.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil, der2_pencil
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (my,mz) :: rho0,gamma0, TT0
      real, dimension (mx,my,mz) :: mom2,mom3!, rho_ux2, rho_uy2
      real, dimension (mx,my,mz) :: rhoE_p
      real, dimension (mx,my,mz,2) ::  rhoE_pU
      real, dimension (my,mz) ::  dYk_dx,dYk_dy,dYk_dz
      real, dimension (ny,nz) :: drho_prefac, KK
      
      real, dimension (ny,nz) :: L_1, L_2, L_3, L_4, L_5
  !    
      real, dimension (my,mz) :: cs0_ar,cs20_ar,dmom2_dy,dmom3_dz
      real, dimension (my,mz,2) :: drhoE_pU!,dYk_dy
   !   real, dimension (mx,my,mz,nchemspec) :: bound_rhs_Y
    !  real, dim(:,nnn)ension (ny,nz) :: bound_rhs_T
      real, dimension (mx,my,mz) :: cs2_full, gamma_full,pp
      real, dimension (nx,ny,nz) :: p_inf
      
      real, dimension (ny,nz,3,3) :: dui_dxj
      real, dimension (my,mz)     :: tmp1,tmp2,tmp3, tmp_rho,tmp_pp       
      real, dimension (ny,nz,3)   :: grad_rho, grad_pp
     ! real, dimension(ny,nz) ::     d2u1_dy2,d2u1_dz2
     ! real, dimension(ny,nz) ::     d2u2_dy2,d2u2_dz2,d2u3_dy2,d2u3_dz2
   !   real, dimension (ny,nz,3) :: dlnT_dxj
      real, dimension (ny,nz) :: T_1_y, T_2_y, T_3_y, T_4_y, T_5_y
      real, dimension (ny,nz) :: T_1_z, T_2_z, T_3_z, T_4_z, T_5_z
!
      integer :: lll, sgn,i,j,k,irho_tmp
      real :: Mach_num, nscbc_sigma_out
!
      intent(inout) :: f
      intent(out) :: df
      intent(in) :: nscbc_sigma_out
!
  
      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case (topbot)
      case ('bot')
        lll = l1; sgn = 1
        if (leos_chemistry) then
          p_inf(1,:,:)=pp_full(l1,m1:m2,n1:n2)
        endif
      case ('top')
        lll = l2; sgn = -1
        if (leos_chemistry) then
          p_inf(nx,:,:)=pp_full(l2,m1:m2,n1:n2)
        endif
      case default
        print*, "bc_nscbc_subin_x: ", topbot, " should be `top' or `bot'"
      endselect
 

      if (leos_chemistry) then
         cs20_ar=cs2_full(lll,:,:)
         cs0_ar=cs2_full(lll,:,:)**0.5
         gamma0=gamma_full(lll,:,:)
         if (ltemperature_nolog) then
          TT0=f(lll,:,:,iTT)
         else
          TT0=exp(f(lll,:,:,ilnTT))
         endif



        if (ldensity_nolog) then
          rho_full = f(:,:,:,irho)
          rho0(:,:) = f(lll,:,:,irho)
          drho_prefac=-1./cs20_ar(m1:m2,n1:n2)
          irho_tmp=irho
        !  call stop_it('bc_nscbc_nref_subout_x: NSCBC works now only for lnrho')
        else
          rho_full = exp(f(:,:,:,ilnrho))
          rho0(:,:) = rho_full(lll,:,:)
          drho_prefac=-1./rho0(m1:m2,n1:n2)/cs20_ar(m1:m2,n1:n2)
          irho_tmp=ilnrho
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
         mom3(lll,:,:)=rho0(:,:)*f(lll,:,:,iuz)
         rhoE_p(lll,:,:)=0.5*rho_full(lll,:,:) &
            *(f(lll,:,:,iux)**2+f(lll,:,:,iuy)**2+f(lll,:,:,iuz)**2) &
             +gamma_full(lll,:,:)/(gamma_full(lll,:,:)-1)*pp(lll,:,:)
         rhoE_pU(lll,:,:,1)=rhoE_p(lll,:,:)*f(lll,:,:,iuy) 
         rhoE_pU(lll,:,:,2)=rhoE_p(lll,:,:)*f(lll,:,:,iuz)


         Mach_num=maxval(f(lll,m1:m2,n1:n2,iux)/cs0_ar(m1:m2,n1:n2))
         KK=nscbc_sigma_out*(1.-Mach_num*Mach_num)*cs0_ar(m1:m2,n1:n2)/Lxyz(1)
      else
        print*,"bc_nscbc_subin_x: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented "//&
            "for the chemistry case."
        print*,"Boundary treatment skipped."
        return
      endif

        call der_onesided_4_slice(rho_full,sgn,grad_rho(:,:,1),lll,1)
        call der_onesided_4_slice(pp,sgn,grad_pp(:,:,1),lll,1)
        call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,1),lll,1)
        call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,1),lll,1)
        call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,1),lll,1)
    !    call der_onesided_4_slice(f,sgn,ilnTT,dlnT_dxj(:,:,1),lll,1)


       if (nygrid /= 1) then
    
         do i=n1,n2
          call der_pencil(2,f(lll,:,i,iux),tmp1(:,i))
          call der_pencil(2,f(lll,:,i,iuy),tmp2(:,i))
          call der_pencil(2,f(lll,:,i,iuz),tmp3(:,i))
          call der_pencil(2,rho_full(lll,:,i),tmp_rho(:,i))
          call der_pencil(2,pp(lll,:,i),tmp_pp(:,i))
          call der_pencil(2,mom2(lll,:,i),dmom2_dy(:,i))
          call der_pencil(2,rhoE_pU(lll,:,i,1),drhoE_pU(:,i,1))
         enddo
  
      !  do i=n1,n2
      !    call der2_pencil(2,f(lll,:,i,iux),tmpy)
      !    d2u1_dy2(:,i-n1+1)=tmpy(:)
      !    call der2_pencil(2,f(lll,:,i,iuy),tmpy)
      !    d2u2_dy2(:,i-n1+1)=tmpy(:)
      !    call der2_pencil(2,f(lll,:,i,iuz),tmpy)
      !    d2u3_dy2(:,i-n1+1)=tmpy(:)
      !    enddo
      else
        tmp3=0
        tmp2=0
        tmp1=0
        tmp_rho=0
        tmp_pp=0
        dmom2_dy=0
        drhoE_pU(:,:,1)=0
     !   d2u1_dy2=0
     !   d2u2_dy2=0
     !   d2u3_dy2=0
      
      endif
        dui_dxj(:,:,1,2)=tmp1(m1:m2,n1:n2)
        dui_dxj(:,:,2,2)=tmp2(m1:m2,n1:n2)
        dui_dxj(:,:,3,2)=tmp3(m1:m2,n1:n2)
        grad_rho(:,:,2)=tmp_rho(m1:m2,n1:n2)
        grad_pp(:,:,2)=tmp_pp(m1:m2,n1:n2)

      if (nzgrid /= 1) then
         do j=m1,m2
          call der_pencil(3,f(lll,j,:,iux),tmp1(j,:))
          call der_pencil(3,f(lll,j,:,iuy),tmp2(j,:))
          call der_pencil(3,f(lll,j,:,iuz),tmp3(j,:))
          call der_pencil(3,rho_full(lll,j,:),tmp_rho(j,:))
          call der_pencil(3,pp(lll,j,:),tmp_pp(j,:))
          call der_pencil(3,mom3(lll,j,:),dmom3_dz(j,:))
          call der_pencil(3,rhoE_pU(lll,j,:,2),drhoE_pU(j,:,2))
         enddo
      else
        tmp3=0
        tmp2=0
        tmp1=0
        tmp_rho=0
        tmp_pp=0
        dmom3_dz=0
        drhoE_pU(:,:,2)=0
      endif
        dui_dxj(:,:,1,3)=tmp1(m1:m2,n1:n2)
        dui_dxj(:,:,2,3)=tmp2(m1:m2,n1:n2)
        dui_dxj(:,:,3,3)=tmp3(m1:m2,n1:n2)
        grad_rho(:,:,3)=tmp_rho(m1:m2,n1:n2)
        grad_pp(:,:,3)=tmp_pp(m1:m2,n1:n2)


!
      select case (topbot)
      case ('bot')
        L_5=KK*(cs20_ar(m1:m2,n1:n2)/gamma0(m1:m2,n1:n2)*&
            rho0(m1:m2,n1:n2)-p_inf(1,1:ny,1:nz)) 
        L_1 = (f(lll,m1:m2,n1:n2,iux) - cs0_ar(m1:m2,n1:n2))*&
           (grad_pp(:,:,1)-rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,1,1))
      case ('top')
        L_1=KK*(cs20_ar(m1:m2,n1:n2)/gamma0(m1:m2,n1:n2)*&
            rho0(m1:m2,n1:n2)-p_inf(nx,1:ny,1:nz))
        L_5 = (f(lll,m1:m2,n1:n2,iux) + cs0_ar(m1:m2,n1:n2))*&
            (grad_pp(:,:,1)+ rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2)*dui_dxj(:,:,1,1))
     endselect
!
       L_2 = f(lll,m1:m2,n1:n2,iux)*(cs20_ar(m1:m2,n1:n2)*grad_rho(:,:,1)-grad_pp(:,:,1))
       L_3 = f(lll,m1:m2,n1:n2,iux)*dui_dxj(:,:,2,1)
       L_4 = f(lll,m1:m2,n1:n2,iux)*dui_dxj(:,:,3,1)
!
       df(lll,m1:m2,n1:n2,irho_tmp) = &
         drho_prefac*(L_2+0.5*(L_5 + L_1)) 
       df(lll,m1:m2,n1:n2,iux) =  &
         -1./(2.*rho0(m1:m2,n1:n2)*cs0_ar(m1:m2,n1:n2))*(L_5 - L_1)
       df(lll,m1:m2,n1:n2,iuy) = -L_3
       df(lll,m1:m2,n1:n2,iuz) = -L_4
       if (ltemperature_nolog) then
        df(lll,m1:m2,n1:n2,ilnTT) = &
          -1./(rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))*(-L_2 &
          +0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1))*TT0(m1:m2,n1:n2) !&
      ! + RHS_T_full(lll,m1:m2,n1:n2)
       else
        df(lll,m1:m2,n1:n2,ilnTT) = &
          -1./(rho0(m1:m2,n1:n2)*cs20_ar(m1:m2,n1:n2))*(-L_2 &
          +0.5*(gamma0(m1:m2,n1:n2)-1.)*(L_5+L_1))
       endif
!NNNNNNNNNN

        if ((nygrid /= 1) .or.(nzgrid /= 1))  then

        
          T_1_y(:,:)=-dmom2_dy(m1:m2,n1:n2)/rho0(m1:m2,n1:n2)
          T_1_z(:,:)=-dmom3_dz(m1:m2,n1:n2)/rho0(m1:m2,n1:n2)
          T_2_y(:,:)=-f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,1,2) 
          T_2_z(:,:)=-f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,1,3) !&
                !   -nu_full(lll,m1:m2,n1:n2)*(d2u1_dy2+d2u1_dz2)
          T_3_y(:,:)=-f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,2,2) &
                     -grad_pp(:,:,2)/rho0(m1:m2,n1:n2) 
         T_3_z(:,:)=-f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,2,3) 
                !   -nu_full(lll,m1:m2,n1:n2)*(d2u2_dy2+d2u2_dz2)
          T_4_y(:,:)=-f(lll,m1:m2,n1:n2,iuy)*dui_dxj(:,:,3,2) 
         T_4_z(:,:)=-f(lll,m1:m2,n1:n2,iuz)*dui_dxj(:,:,3,3) &
                     -grad_pp(:,:,3)/rho0(m1:m2,n1:n2)
          T_5_y(:,:)= drho_prefac(:,:)*(gamma0(m1:m2,n1:n2)-1.) &
                      *gamma0(m1:m2,n1:n2)*drhoE_pU(m1:m2,n1:n2,1)
          T_5_z(:,:)= drho_prefac(:,:)*(gamma0(m1:m2,n1:n2)-1.) &
                      *gamma0(m1:m2,n1:n2)*(drhoE_pU(m1:m2,n1:n2,2))
 
         df(lll,m1:m2,n1:n2,irho_tmp) = df(lll,m1:m2,n1:n2,irho_tmp) + T_1_y + T_1_z
         df(lll,m1:m2,n1:n2,iux)      = df(lll,m1:m2,n1:n2,iux)      + T_2_y + T_2_z
         df(lll,m1:m2,n1:n2,iuy)      = df(lll,m1:m2,n1:n2,iuy)      + T_3_y + T_3_z
         df(lll,m1:m2,n1:n2,iuz)      = df(lll,m1:m2,n1:n2,iuz)      + T_4_y + T_4_z
        if (ltemperature_nolog) then
         df(lll,m1:m2,n1:n2,ilnTT)    = df(lll,m1:m2,n1:n2,ilnTT)   &
                                + (T_5_y + T_5_z)*TT0(m1:m2,n1:n2)
        else
         df(lll,m1:m2,n1:n2,ilnTT)    = df(lll,m1:m2,n1:n2,ilnTT)    + T_5_y + T_5_z
        endif

 !!!
!!!  Corner points are described in bc_nscbc_nref_subout_y and bc_nscbc_nref_subout_z
!!!


       endif







!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nchemspec>1) then
       do k=1,nchemspec
          call der_onesided_4_slice(f,sgn,ichemspec(k),dYk_dx(m1:m2,n1:n2),lll,1)
         do i=n1,n2
          call der_pencil(2,f(lll,:,i,ichemspec(k)),dYk_dy(:,i))
         enddo
         do i=m1,m2
           call der_pencil(3,f(lll,i,:,ichemspec(k)),dYk_dz(i,:))
         enddo
!NMNMNMN
           df(lll,m1:m2,n1:n2,ichemspec(k))=&
              -f(lll,m1:m2,n1:n2,iux)*dYk_dx(m1:m2,n1:n2) &
              -f(lll,m1:m2,n1:n2,iuy)*dYk_dy(m1:m2,n1:n2) &
              -f(lll,m1:m2,n1:n2,iuz)*dYk_dz(m1:m2,n1:n2) &
              + RHS_Y_full(lll,m1:m2,n1:n2,k)

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
  
    endsubroutine bc_nscbc_nref_subout_x
!***********************************************************************

!***********************************************************************************
 !***********************************************************
    subroutine bc_nscbc_nref_subout_y(f,df,topbot,nscbc_sigma_out)
!
!   nscbc case
!   subsonic non-reflecting outflow boundary conditions
!
!   16-jun-09/natalia: coded.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (mx,mz) :: rho0,gamma0
      real, dimension (mx,my,mz) :: mom1, mom3, pp!, rho_ux2, rho_uy2
      real, dimension (mx,mz) ::  dYk_dy,dYk_dx,dYk_dz
      real, dimension (nx,nz) :: drho_prefac, KK
      real, dimension (nx,nz) :: M_1, M_2, M_3, M_4, M_5
      real, dimension (nx,nz)  :: L_1, L_2, L_3, L_4, L_5
      real, dimension (nx,nz)  :: N_1, N_2, N_3, N_4, N_5
      real, dimension (nx,nz) :: T_1_x, T_2_x, T_3_x, T_4_x, T_5_x
      real, dimension (nx,nz) :: T_1_z, T_2_z, T_3_z, T_4_z, T_5_z
      real, dimension (mx,mz) :: cs0_ar,cs20_ar,dmom1_dx,dmom3_dz
      real, dimension (mx,my,mz) :: rhoE_p
      real, dimension (mx,my,mz,2) ::   rhoE_pU
      real, dimension (mx,mz,2) :: drhoE_pU
      !real, dimension (mx,mz) ::  !dYk_dy
  !    real, dimension (nx,nz,nchemspec) :: bound_rhs_Y
    !  real, dimension (nx,nz) :: bound_rhs_T
      real, dimension (mx,my,mz) :: cs2_full, gamma_full
      real, dimension (nx,ny,nz) :: p_inf
      real, dimension (nx,nz,3,3) :: dui_dxj
      real, dimension (mx,mz)     :: tmp11,tmp21,tmp31,tmp13,tmp23,tmp33
      real, dimension (mx,mz)     :: tmp1_rho,tmp3_rho,tmp1_pp,tmp3_pp
      real, dimension (nx,nz,3)   :: grad_rho, grad_pp

!
      integer :: mmm, sgn,i,j,k, irho_tmp, nn, nnn, ll, lll
      real :: Mach_num,nscbc_sigma_out
      logical :: lcorner_x=.false.,lcorner_z=.false.
!
      intent(inout) :: f
      intent(out) :: df
      intent(in) :: nscbc_sigma_out
!
      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!

      select case (topbot)
      case ('bot')
        mmm = m1; sgn = 1
        if (leos_chemistry) then
          p_inf(:,1,:)=pp_full(l1:l2,m1,n1:n2)
        endif
      case ('top')
        mmm = m2; sgn = -1
        if (leos_chemistry) then
          p_inf(:,ny,:)=pp_full(l1:l2,m2,n1:n2)
        endif
      case default
        print*, "bc_nscbc_subout_y: ", topbot, " should be `top' or `bot'"
      endselect


      if (leos_chemistry) then
         cs20_ar=cs2_full(:,mmm,:)
         cs0_ar=cs2_full(:,mmm,:)**0.5
         gamma0=gamma_full(:,mmm,:)

        if (ldensity_nolog) then
          rho_full = f(:,:,:,irho)
          rho0(:,:) = f(:,mmm,:,irho)
          drho_prefac=-1./cs20_ar(l1:l2,n1:n2)
          irho_tmp=irho

         call stop_it('bc_nscbc_nref_subout_y: NSCBC works now only for lnrho')

        else
          rho_full = exp(f(:,:,:,ilnrho))
          rho0(:,:) = rho_full(:,mmm,:)
          drho_prefac=-1./rho0(l1:l2,n1:n2)/cs20_ar(l1:l2,n1:n2)
          irho_tmp=ilnrho
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
        

         rhoE_p(:,mmm,:)=0.5*rho_full(:,mmm,:) &
            *(f(:,mmm,:,iux)**2+f(:,mmm,:,iuy)**2+f(:,mmm,:,iuz)**2) &
             +gamma_full(:,mmm,:)/(gamma_full(:,mmm,:)-1)*pp_full(:,mmm,:)
         rhoE_pU(:,mmm,:,1)=rhoE_p(:,mmm,:)*f(:,mmm,:,iux) 
         rhoE_pU(:,mmm,:,2)=rhoE_p(:,mmm,:)*f(:,mmm,:,iuz)
         mom1(:,mmm,:)=rho0(:,:)*f(:,mmm,:,iux)
         mom3(:,mmm,:)=rho0(:,:)*f(:,mmm,:,iuz)

         Mach_num=maxval(f(l1:l2,mmm,n1:n2,iuy)/cs0_ar(l1:l2,n1:n2))
         KK=nscbc_sigma_out*(1.-Mach_num*Mach_num)*cs0_ar(l1:l2,n1:n2)/Lxyz(2)


      else
        print*,"bc_nscbc_subin_y: leos_idealgas=",leos_idealgas,"."
        print*,"NSCBC subsonic inflos is only implemented "//&
            "for the chemistry case."
        print*,"Boundary treatment skipped."
        return
      endif

      call der_onesided_4_slice(rho_full,sgn,grad_rho(:,:,2),mmm,2)
      call der_onesided_4_slice(pp,sgn,grad_pp(:,:,2),mmm,2)
      call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,2),mmm,2)
      call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,2),mmm,2)
      call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,2),mmm,2)

      if (nxgrid /= 1) then
         do i=n1,n2
          call der_pencil(1,f(:,mmm,i,iux),tmp11(:,i))
          call der_pencil(1,f(:,mmm,i,iuy),tmp21(:,i))
          call der_pencil(1,f(:,mmm,i,iuz),tmp31(:,i))
          call der_pencil(1,rho_full(:,mmm,i),tmp1_rho(:,i))
          call der_pencil(1,pp(:,mmm,i),tmp1_pp(:,i))
          call der_pencil(1,mom1(:,mmm,i),dmom1_dx(:,i))
          call der_pencil(1,rhoE_pU(:,mmm,i,1),drhoE_pU(:,i,1))
         enddo

      else
        tmp31=0
        tmp21=0
        tmp11=0
        tmp1_rho=0
        tmp1_pp=0
        dmom1_dx=0
        drhoE_pU(:,:,1)=0
      endif
        dui_dxj(:,:,1,1)=tmp11(l1:l2,n1:n2)
        dui_dxj(:,:,2,1)=tmp21(l1:l2,n1:n2)
        dui_dxj(:,:,3,1)=tmp31(l1:l2,n1:n2)
        grad_rho(:,:,1)=tmp1_rho(l1:l2,n1:n2)
        grad_pp(:,:,1)=tmp1_pp(l1:l2,n1:n2)

      if (nzgrid /= 1) then
        
         do j=l1,l2
          call der_pencil(3,f(j,mmm,:,iux),tmp13(j,:))
          call der_pencil(3,f(j,mmm,:,iuy),tmp23(j,:))
          call der_pencil(3,f(j,mmm,:,iuz),tmp33(j,:))
          call der_pencil(3,rho_full(j,mmm,:),tmp3_rho(j,:))
          call der_pencil(3,pp(j,mmm,:),tmp3_pp(j,:))
          call der_pencil(3,mom3(j,mmm,:),dmom3_dz(j,:))
          call der_pencil(3,rhoE_pU(j,mmm,:,2),drhoE_pU(j,:,2))

         enddo

      else
        tmp33=0
        tmp23=0
        tmp13=0
        tmp3_rho=0
        tmp3_pp=0
        dmom3_dz=0
        drhoE_pU(:,:,2)=0
      endif
        dui_dxj(:,:,1,3)=tmp13(l1:l2,n1:n2)
        dui_dxj(:,:,2,3)=tmp23(l1:l2,n1:n2)
        dui_dxj(:,:,3,3)=tmp33(l1:l2,n1:n2)
        grad_rho(:,:,3)=tmp3_rho(l1:l2,n1:n2)
        grad_pp(:,:,3)=tmp3_pp(l1:l2,n1:n2)


!
      select case (topbot)
      case ('bot')
        M_5=KK*(cs20_ar(l1:l2,n1:n2)/gamma0(l1:l2,n1:n2)*&
            rho0(l1:l2,n1:n2)-p_inf(1:nx,1,1:nz))
        M_1 = (f(l1:l2,mmm,n1:n2,iuy) - cs0_ar(l1:l2,n1:n2))*&
            (grad_pp(:,:,2)- rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,2,2))
      case ('top')
        M_1=KK*(cs20_ar(l1:l2,n1:n2)/gamma0(l1:l2,n1:n2)*&
            rho0(l1:l2,n1:n2)-p_inf(1:nx,ny,1:nz))
        M_5 = (f(l1:l2,mmm,n1:n2,iuy) + cs0_ar(l1:l2,n1:n2))*&
            (grad_pp(:,:,2)+ rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,2,2))
      endselect
!
      M_2 = f(l1:l2,mmm,n1:n2,iuy)*(cs20_ar(l1:l2,n1:n2)*grad_rho(:,:,2)-grad_pp(:,:,2))
      M_3 = f(l1:l2,mmm,n1:n2,iuy)*dui_dxj(:,:,1,2)
      M_4 = f(l1:l2,mmm,n1:n2,iuy)*dui_dxj(:,:,3,2)

      df(l1:l2,mmm,n1:n2,irho_tmp) = drho_prefac*(M_2+0.5*(M_5 + M_1)) 
      df(l1:l2,mmm,n1:n2,iux) = -M_3 
      df(l1:l2,mmm,n1:n2,iuy) = -1./&
          (2.*rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2))*(M_5 - M_1) !&
      df(l1:l2,mmm,n1:n2,iuz) = -M_4
      df(l1:l2,mmm,n1:n2,ilnTT) = drho_prefac(:,:)*(-M_2 &
          +0.5*(gamma0(l1:l2,n1:n2)-1.)*(M_5+M_1))


!NNNNNNNN
       if ((nxgrid /= 1) .or. (nzgrid /= 1)) then

          T_1_x(:,:)=-dmom1_dx(l1:l2,n1:n2)/rho0(l1:l2,n1:n2)
          T_1_z(:,:)=-dmom3_dz(l1:l2,n1:n2)/rho0(l1:l2,n1:n2)
          T_3_x(:,:)=-f(l1:l2,mmm,n1:n2,iux)*dui_dxj(:,:,2,1) 
          T_3_z(:,:)=-f(l1:l2,mmm,n1:n2,iuz)*dui_dxj(:,:,2,3)
          T_2_x(:,:)=-f(l1:l2,mmm,n1:n2,iux)*dui_dxj(:,:,1,1) &
                    -grad_pp(:,:,1)/rho0(l1:l2,n1:n2)
          T_2_z(:,:)=-f(l1:l2,mmm,n1:n2,iuz)*dui_dxj(:,:,1,3) 
          T_4_x(:,:)=-f(l1:l2,mmm,n1:n2,iux)*dui_dxj(:,:,3,1) 
          T_4_z(:,:)=-f(l1:l2,mmm,n1:n2,iuz)*dui_dxj(:,:,3,3) &
                   -grad_pp(:,:,3)/rho0(l1:l2,n1:n2)
          T_5_x(:,:)=+drho_prefac(:,:)*(gamma0(l1:l2,n1:n2)-1.)*gamma0(l1:l2,n1:n2) &
                   *(drhoE_pU(l1:l2,n1:n2,1))
          T_5_z(:,:)=+drho_prefac(:,:)*(gamma0(l1:l2,n1:n2)-1.)*gamma0(l1:l2,n1:n2) &
                   *(drhoE_pU(l1:l2,n1:n2,2))

  
        df(l1:l2,mmm,n1:n2,irho_tmp) = df(l1:l2,mmm,n1:n2,irho_tmp) + T_1_x+T_1_z
        df(l1:l2,mmm,n1:n2,iux) =      df(l1:l2,mmm,n1:n2,iux)      + T_2_x+T_2_z
        df(l1:l2,mmm,n1:n2,iuy) =      df(l1:l2,mmm,n1:n2,iuy)      + T_3_x+T_3_z
        df(l1:l2,mmm,n1:n2,iuz) =      df(l1:l2,mmm,n1:n2,iuz)      + T_4_x+T_4_z
        df(l1:l2,mmm,n1:n2,ilnTT) =    df(l1:l2,mmm,n1:n2,ilnTT)    + T_5_x+T_5_z

   !    if ((nxgrid /= 1) .and. (nzgrid /= 1)) then
    ! if ((nxgrid /= 1) .and. (nzgrid == 1)) then
   
!!!
!!! Corner points
!!!
!!!

   
        if (nxgrid /= 1) then 
           L_1(:,:)=(f(l1:l2,mmm,n1:n2,iux) - cs0_ar(l1:l2,n1:n2))&
             *(grad_pp(:,:,1)-rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,1,1))
           L_2(:,:)=f(l1:l2,mmm,n1:n2,iux)*(cs20_ar(l1:l2,n1:n2) &
                  *grad_rho(:,:,1)-grad_pp(:,:,1))
           L_3(:,:)=f(l1:l2,mmm,n1:n2,iux)*dui_dxj(:,:,2,1)
           L_4(:,:)=f(l1:l2,mmm,n1:n2,iux)*dui_dxj(:,:,3,1)
           L_5(:,:)=(f(l1:l2,mmm,n1:n2,iux) + cs0_ar(l1:l2,n1:n2))&
            *(grad_pp(:,:,1)+ rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,1,1))


           if (x(l2)==Lxyz(1)) then
            L_1(nx,:)=KK(nx,:)*(cs20_ar(l2,n1:n2)/gamma0(l2,n1:n2) &
               *rho0(l2,n1:n2)-p_inf(nx,mmm-3,:)) 
           endif
           if  (x(l1)==xyz0(1)) then
            L_5(1,:)=KK(1,:)*(cs20_ar(l1,n1:n2)/gamma0(l1,n1:n2)*&
             rho0(l1,n1:n2)-p_inf(1,mmm-3,:))
           endif

        else
        L_1=0; L_2=0; L_3=0; L_4=0; L_5=0 
       endif
!
       if (nzgrid /= 1)  then
         N_1(:,:)=(f(l1:l2,mmm,n1:n2,iuz) - cs0_ar(l1:l2,n1:n2))&
          *(grad_pp(:,:,3)-rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,3,3))

         N_2(:,:)=f(l1:l2,mmm,n1:n2,iuz)*dui_dxj(:,:,1,3)

         N_3(:,:)=f(l1:l2,mmm,n1:n2,iuz)*dui_dxj(:,:,2,3)

         N_4(:,:)=f(l1:l2,mmm,n1:n2,iuz) &
            *(cs20_ar(l1:l2,n1:n2)*grad_rho(:,:,3)-grad_pp(:,:,3))

         N_5(:,:)=(f(l1:l2,mmm,n1:n2,iuz) + cs0_ar(l1:l2,n1:n2))&
            *(grad_pp(:,:,3)+ rho0(l1:l2,n1:n2)*cs0_ar(l1:l2,n1:n2)*dui_dxj(:,:,3,3))

         if (z(n2)==Lxyz(3)) then
          N_1(:,nz)=KK(:,nz)*(cs20_ar(l1:l2,n2)/gamma0(l1:l2,n2)*&
            rho0(l1:l2,n2)-p_inf(:,mmm-3,nz))
         endif
         if (z(n1)==xyz0(3)) then
          N_5(:,1)=KK(:,1)*(cs20_ar(l1:l2,n1)/gamma0(l1:l2,n1)*&
             rho0(l1:l2,n1)-p_inf(:,mmm-3,1))
         endif


       else
        N_1=0; N_2=0; N_3=0; N_4=0; N_5=0
       endif


         do i=1,2
         if (i==1) then
            ll=l1; lll=1
            if (x(l1)==xyz0(1)) lcorner_x=.true.
         elseif (i==2) then
            ll=l2; lll=nx
            if (x(l2)==Lxyz(1)) lcorner_x=.true.
         endif
         if (lcorner_x)  then
          df(ll,mmm,n1:n2,irho_tmp) = &
            drho_prefac(lll,:)*(L_2(lll,:)+0.5*(L_5(lll,:) + L_1(lll,:))) &
          + drho_prefac(lll,:)*(M_2(lll,:) &
          +0.5*(M_5(lll,:) + M_1(lll,:))) +T_1_z(lll,:)

          df(ll,mmm,n1:n2,iux) =  -1./&
           (2.*rho0(ll,n1:n2)*cs0_ar(ll,n1:n2))*(L_5(lll,:) - L_1(lll,:)) &
           -M_3(lll,:)+T_2_z(lll,:)
          df(ll,mmm,n1:n2,iuy) =  - L_3(lll,:) - 1./&
           (2.*rho0(ll,n1:n2)*cs0_ar(ll,n1:n2))*(M_5(lll,:) - M_1(lll,:)) &
           + T_3_z(lll,:)
          df(ll,mmm,n1:n2,iuz) =  - L_4(lll,:) - M_4(lll,:)  +T_4_z(lll,:)

          df(ll,mmm,n1:n2,ilnTT) = drho_prefac(lll,:)*(-L_2(lll,:) &
           +0.5*(gamma0(ll,n1:n2)-1.)*(L_5(lll,:)+L_1(lll,:)))  &
           +drho_prefac(lll,:)*(-M_2(lll,:) &
           +0.5*(gamma0(ll,n1:n2)-1.)*(M_5(lll,:)+M_1(lll,:))) +T_5_z(lll,:)
          lcorner_x=.false.
        endif
        enddo

        do i=1,2 
         if (i==1) then
          nn=n1; nnn=1
          if (z(n1)==xyz0(3)) lcorner_z=.true.
         elseif (i==2) then
          nn=n2; nnn=nz
          if (z(n2)==Lxyz(3))  lcorner_z=.true.
         endif

         if  (lcorner_z) then
          df(l1:l2,mmm,nn,irho_tmp) = T_1_x(:,nnn) &
          + drho_prefac(:,nnn)*(M_2(:,nnn)+0.5*(M_5(:,nnn) + M_1(:,nnn))) &
          + drho_prefac(:,nnn)*(N_4(:,nnn)+0.5*(N_5(:,nnn) + N_1(:,nnn)))
          df(l1:l2,mmm,nn,iux) = T_2_x(:,nnn) -M_3(:,nnn)-N_2(:,nnn)
          df(l1:l2,mmm,nn,iuy) = T_3_x(:,nnn) - 1./&
           (2.*rho0(l1:l2,nn)*cs0_ar(l1:l2,nn))*(M_5(:,nnn) &
                - M_1(:,nnn))-N_3(:,nnn)
          df(l1:l2,mmm,nn,iuz) =  T_4_x(:,nnn) - M_4(:,nnn)  - 1./&
           (2.*rho0(l1:l2,nn)*cs0_ar(l1:l2,nn))*(N_5(:,nnn) - N_1(:,nnn))
          df(l1:l2,mmm,nn,ilnTT) = T_5_x(:,nnn) &
           +drho_prefac(:,nnn)*(-M_2(:,nnn) &
           +0.5*(gamma0(l1:l2,nn)-1.)*(M_5(:,nnn)+M_1(:,nnn))) &
           +drho_prefac(:,nnn)*(-N_4(:,nnn) &
           +0.5*(gamma0(l1:l2,nn)-1.)*(N_5(:,nnn)+N_1(:,nnn)))
           lcorner_z=.false.
         endif
        enddo


        do i=1,2 
         if (i==1) then
          nn=n1; nnn=1
          if (z(n1)==xyz0(3)) lcorner_z=.true.
         elseif (i==2) then
          nn=n2; nnn=nz
          if (z(n2)==Lxyz(3))  lcorner_z=.true.
         endif

        do j=1,2
         if (j==1) then
           ll=l1; lll=1
           if (x(l1)==xyz0(1)) lcorner_x=.true.
         elseif (j==2) then
           ll=l2; lll=nx
           if (x(l2)==Lxyz(1)) lcorner_x=.true.
         endif

!NNNNNNNN

    !   if ((x(l1)==xyz0(1))  .or. (x(l2)==Lxyz(1)))  lcorner_x=.true.
    !   if ((z(n1)==xyz0(3))  .or. (z(n2)==Lxyz(3)))  lcorner_z=.true.


       if ((lcorner_x)  .and. (lcorner_z)) then
        df(ll,mmm,nn,irho_tmp) = &
          drho_prefac(lll,nnn)*(L_2(lll,nnn)+0.5*(L_5(lll,nnn) + L_1(lll,nnn))) &
          + drho_prefac(lll,nnn)*(M_2(lll,nnn)+0.5*(M_5(lll,nnn) + M_1(lll,nnn))) &
          + drho_prefac(lll,nnn)*(N_4(lll,nnn)+0.5*(N_5(lll,nnn) + N_1(lll,nnn)))
        df(ll,mmm,nn,iux) =  -1./&
          (2.*rho0(ll,nn)*cs0_ar(ll,nn))*(L_5(lll,nnn) - L_1(lll,nnn)) &
          -M_3(lll,nnn)-N_2(lll,nnn)
        df(ll,mmm,nn,iuy) =  - L_3(lll,nnn) - 1./&
          (2.*rho0(ll,nn)*cs0_ar(ll,nn))*(M_5(lll,nnn) - M_1(lll,nnn))-N_3(lll,nnn)
        df(ll,mmm,nn,iuz) =  - L_4(lll,nnn) - M_4(lll,nnn)  - 1./&
          (2.*rho0(ll,nn)*cs0_ar(ll,nn))*(N_5(lll,nnn) - N_1(lll,nnn))
        df(ll,mmm,nn,ilnTT) = drho_prefac(lll,nnn)*(-L_2(lll,nnn) &
          +0.5*(gamma0(ll,nn)-1.)*(L_5(lll,nnn)+L_1(lll,nnn)))  &
          +drho_prefac(lll,nnn)*(-M_2(lll,nnn) &
          +0.5*(gamma0(ll,nn)-1.)*(M_5(lll,nnn)+M_1(lll,nnn))) &
           +drho_prefac(lll,nnn)*(-N_4(lll,nnn) &
          +0.5*(gamma0(ll,nn)-1.)*(N_5(lll,nnn)+N_1(lll,nnn)))
           lcorner_x=.false.
           lcorner_z=.false.
       endif
      enddo
      enddo

    endif



      if (nchemspec>1) then
       do k=1,nchemspec
          call der_onesided_4_slice(f,sgn,ichemspec(k),dYk_dy(l1:l2,n1:n2),mmm,2)
         do i=n1,n2
          call der_pencil(1,f(:,mmm,i,ichemspec(k)),dYk_dx(:,i))
         enddo
         do i=l1,l2
          call der_pencil(3,f(i,mmm,:,ichemspec(k)),dYk_dz(i,:))
         enddo


          df(l1:l2,mmm,n1:n2,ichemspec(k))=&
              -f(l1:l2,mmm,n1:n2,iux)*dYk_dx(l1:l2,n1:n2) &
              -f(l1:l2,mmm,n1:n2,iuy)*dYk_dy(l1:l2,n1:n2) &
              -f(l1:l2,mmm,n1:n2,iuz)*dYk_dz(l1:l2,n1:n2) &
              +RHS_Y_full(l1:l2,mmm,n1:n2,k)

       ! if (lfilter) then
         do i=l1,l2
         do j=n1,n2
           if ((f(i,mmm,j,ichemspec(k))+df(i,mmm,j,ichemspec(k))*dt)<-1e-25 ) then
             df(i,mmm,j,ichemspec(k))=-1e-25*dt
           endif
           if ((f(i,mmm,j,ichemspec(k))+df(i,mmm,j,ichemspec(k))*dt)>1.) then
             f(i,mmm,j,ichemspec(k))=1.*dt
           endif
         enddo
         enddo
       ! endif
       enddo
      endif
!
    endsubroutine bc_nscbc_nref_subout_y
!***********************************************************************************
  
!***********************************************************************************
!***********************************************************************************
   subroutine bc_nscbc_nref_subout_z(f,df,topbot,nscbc_sigma_out)
!
!   nscbc case
!   subsonic non-reflecting outflow boundary conditions
!
!   16-jun-09/natalia: coded.
!
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (mx,my) :: rho0,gamma0  
       real, dimension (mx,my,mz) :: mom1,mom2,pp
         real, dimension (mx,my,mz) :: rhoE_p
      real, dimension (mx,my,mz,2) ::   rhoE_pU=0.
      real, dimension (mx,my,2) :: drhoE_pU
      real, dimension (mx,my) ::  dYk_dz=0.,dYk_dx=0.,dYk_dy=0.
      real, dimension (nx,ny) :: drho_prefac, KK, N_1, N_2, N_3,N_4, N_5
      real, dimension (nx,ny)  :: M_1, M_2, M_3, M_4, M_5
      real, dimension (nx,ny)  :: L_1, L_2, L_3, L_4, L_5

      real, dimension (mx,my) :: cs0_ar,cs20_ar,dmom1_dx,dmom2_dy
      real, dimension (mx,my,mz) :: cs2_full, gamma_full
      real, dimension (nx,ny,nz) :: p_inf
      real, dimension(nx,ny,3,3) :: dui_dxj=0.
       real, dimension (mx,my)     :: tmp11,tmp21,tmp31,tmp12,tmp22,tmp32
      real, dimension (mx,my)     :: tmp1_rho,tmp2_rho,tmp1_pp,tmp2_pp
      real, dimension (nx,ny,3)   :: grad_rho=0., grad_pp=0.
      real, dimension (nx,ny) :: T_1_x, T_2_x, T_3_x, T_4_x, T_5_x
      real, dimension (nx,ny) :: T_1_y, T_2_y, T_3_y, T_4_y, T_5_y
!
      integer :: nnn, sgn,i,j,k, mmm,mm,lll,ll, irho_tmp
      real :: Mach_num,nscbc_sigma_out

      logical :: lcorner_x=.false.,lcorner_y=.false.
!
      intent(inout) :: f
      intent(out) :: df
      intent(in) :: nscbc_sigma_out
!
      if (leos_chemistry) then
        call get_cs2_full(cs2_full)
        call get_gamma_full(gamma_full)
      endif
!
      select case (topbot)
      case ('bot')
        nnn = n1; sgn = 1
        if (leos_chemistry) then
          p_inf(:,:,1)=pp_full(l1:l2,m1:m2,n1)
        endif
      case ('top')
        nnn = n2; sgn = -1
        if (leos_chemistry) then
          p_inf(:,:,nz)=pp_full(l1:l2,m1:m2,n2)
        endif
      case default
        print*, "bc_nscbc_subout_z: ", topbot, " should be `top' or `bot'"
      endselect

        if (leos_chemistry) then
         cs20_ar=cs2_full(:,:,nnn)
         cs0_ar=cs2_full(:,:,nnn)**0.5
         gamma0=gamma_full(:,:,nnn)

        if (ldensity_nolog) then
          rho_full = f(:,:,:,irho)
          rho0(:,:) = f(:,:,nnn,irho)
          drho_prefac=-1./cs20_ar(l1:l2,m1:m2)
           irho_tmp=irho
        else
          rho_full = exp(f(:,:,:,ilnrho))
          rho0(:,:) = rho_full(:,:,nnn)
          drho_prefac=-1./rho0(l1:l2,m1:m2)/cs20_ar(l1:l2,m1:m2)
           irho_tmp=ilnrho
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

         mom1(:,:,nnn)=rho0(:,:)*f(:,:,nnn,iux)
         mom2(:,:,nnn)=rho0(:,:)*f(:,:,nnn,iuy)
         rhoE_p(:,:,nnn)=0.5*rho_full(:,:,nnn) &
            *(f(:,:,nnn,iux)**2+f(:,:,nnn,iuy)**2+f(:,:,nnn,iuz)**2) &
             +gamma_full(:,:,nnn)/(gamma_full(:,:,nnn)-1)*pp(:,:,nnn)
         rhoE_pU(:,:,nnn,1)=rhoE_p(:,:,nnn)*f(:,:,nnn,iux) 
         rhoE_pU(:,:,nnn,2)=rhoE_p(:,:,nnn)*f(:,:,nnn,iuy)

        Mach_num=maxval(f(l1:l2,m1:m2,nnn,iuz)/cs0_ar(l1:l2,m1:m2))
        KK=nscbc_sigma_out*(1.-Mach_num*Mach_num)*cs0_ar(l1:l2,m1:m2)/Lxyz(3)

      else
        print*,"bc_nscbc_subin_y: leos_idealgas=",leos_idealgas,"."
       print*,"NSCBC subsonic inflos is only implemented "//&
           "for the chemistry case."
        print*,"Boundary treatment skipped."
        return
      endif

      call der_onesided_4_slice(rho_full,sgn,grad_rho(:,:,3),nnn,3)
      call der_onesided_4_slice(pp,sgn,grad_pp(:,:,3),nnn,3)
      call der_onesided_4_slice(f,sgn,iux,dui_dxj(:,:,1,3),nnn,3)
      call der_onesided_4_slice(f,sgn,iuy,dui_dxj(:,:,2,3),nnn,3)
      call der_onesided_4_slice(f,sgn,iuz,dui_dxj(:,:,3,3),nnn,3)

      if (nxgrid /= 1) then  

          do i=m1,m2
          call der_pencil(1,f(:,i,nnn,iux),tmp11(:,i))
          call der_pencil(1,f(:,i,nnn,iuy),tmp21(:,i))
          call der_pencil(1,f(:,i,nnn,iuz),tmp31(:,i))
          call der_pencil(1,rho_full(:,i,nnn),tmp1_rho(:,i))
          call der_pencil(1,pp(:,i,nnn),tmp1_pp(:,i))
          call der_pencil(1,mom1(:,i,nnn),dmom1_dx(:,i))
          call der_pencil(1,rhoE_pU(:,i,nnn,1),drhoE_pU(:,i,1))
         enddo


      else
        tmp31=0
        tmp21=0
        tmp11=0
        tmp1_rho=0
        tmp1_pp=0
        dmom1_dx=0
        drhoE_pU(:,:,1)=0
      endif
        dui_dxj(:,:,1,1)=tmp11(l1:l2,m1:m2)
        dui_dxj(:,:,2,1)=tmp21(l1:l2,m1:m2)
        dui_dxj(:,:,3,1)=tmp31(l1:l2,m1:m2)
        grad_rho(:,:,1)=tmp1_rho(l1:l2,m1:m2)
        grad_pp(:,:,1)=tmp1_pp(l1:l2,m1:m2)

      if (nygrid /= 1) then

          do i=l1,l2
          call der_pencil(2,f(i,:,nnn,iux),tmp12(i,:))
          call der_pencil(2,f(i,:,nnn,iuy),tmp22(i,:))
          call der_pencil(2,f(i,:,nnn,iuz),tmp32(i,:))
          call der_pencil(2,rho_full(i,:,nnn),tmp2_rho(i,:))
          call der_pencil(2,pp(i,:,nnn),tmp2_pp(i,:))
          call der_pencil(2,mom2(i,:,nnn),dmom2_dy(i,:))
          call der_pencil(2,rhoE_pU(i,:,nnn,2),drhoE_pU(i,:,2))
         enddo
      else
        tmp32=0
        tmp22=0
        tmp12=0
        tmp2_rho=0
        tmp2_pp=0
        dmom2_dy=0
        drhoE_pU(:,:,2)=0
      endif
        dui_dxj(:,:,1,2)=tmp12(l1:l2,m1:m2)
        dui_dxj(:,:,2,2)=tmp22(l1:l2,m1:m2)
        dui_dxj(:,:,3,2)=tmp32(l1:l2,m1:m2)
        grad_rho(:,:,2)=tmp2_rho(l1:l2,m1:m2)
        grad_pp(:,:,2)=tmp2_pp(l1:l2,m1:m2)

!
      select case (topbot)
      case ('bot')
        N_5=KK*(cs20_ar(l1:l2,m1:m2)/gamma0(l1:l2,m1:m2)*&
            rho0(l1:l2,m1:m2)-p_inf(1:nx,1:ny,1))
        N_1 = (f(l1:l2,m1:m2,nnn,iuz) - cs0_ar(l1:l2,m1:m2))*&
            (grad_pp(:,:,3)- rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,3,3))
      case ('top')
        N_1=KK*(cs20_ar(l1:l2,m1:m2)/gamma0(l1:l2,m1:m2)*&
          rho0(l1:l2,m1:m2)-p_inf(1:nx,1:ny,nz))
        N_5 = (f(l1:l2,m1:m2,nnn,iuz) + cs0_ar(l1:l2,m1:m2))*&
            (grad_pp(:,:,3)+ rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,3,3))
      endselect
!
        N_2 = f(l1:l2,m1:m2,nnn,iuz)  &
           *(cs20_ar(l1:l2,m1:m2)*grad_rho(:,:,3)-grad_pp(:,:,3))
        N_3 = f(l1:l2,m1:m2,nnn,iuz)*dui_dxj(:,:,1,3)
        N_4 = f(l1:l2,m1:m2,nnn,iuz)*dui_dxj(:,:,2,3)
!

        df(l1:l2,m1:m2,nnn,irho_tmp) = &
            drho_prefac*(N_2+0.5*(N_5 + N_1))
        df(l1:l2,m1:m2,nnn,iuz) = -1./&
          (2.*rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2))*(N_5 - N_1) 
        df(l1:l2,m1:m2,nnn,iux) = -N_3
        df(l1:l2,m1:m2,nnn,iuy) = -N_4
        df(l1:l2,m1:m2,nnn,ilnTT) = -1./&
         (rho0(l1:l2,m1:m2)*cs20_ar(l1:l2,m1:m2))*(-N_2 &
        +0.5*(gamma0(l1:l2,m1:m2)-1.)*(N_5+N_1))
!NNNNNNNNN

       if ((nxgrid /= 1) .or. (nygrid /= 1)) then

         T_1_x(:,:)=-dmom1_dx(l1:l2,m1:m2)/rho0(l1:l2,m1:m2)
         T_1_y(:,:)=-dmom2_dy(l1:l2,m1:m2)/rho0(l1:l2,m1:m2)
         T_2_x(:,:)=-f(l1:l2,m1:m2,nnn,iux)*dui_dxj(:,:,1,1) &
                  -grad_pp(:,:,1)/rho0(l1:l2,m1:m2)
         T_2_y(:,:)=-f(l1:l2,m1:m2,nnn,iuy)*dui_dxj(:,:,1,2)

         T_3_x(:,:)=-f(l1:l2,m1:m2,nnn,iux)*dui_dxj(:,:,2,1)
         T_3_y(:,:)=-f(l1:l2,m1:m2,nnn,iuy)*dui_dxj(:,:,2,2)  &
                   -grad_pp(:,:,2)/rho0(l1:l2,m1:m2)
         T_4_x(:,:)=-f(l1:l2,m1:m2,nnn,iux)*dui_dxj(:,:,3,1) 
         T_4_y(:,:)=-f(l1:l2,m1:m2,nnn,iuy)*dui_dxj(:,:,3,2)
         T_5_x(:,:)=+drho_prefac(:,:)*(gamma0(l1:l2,m1:m2)-1.)*gamma0(l1:l2,m1:m2) &
                   *(drhoE_pU(l1:l2,m1:m2,1))
         T_5_y(:,:)=+drho_prefac(:,:)*(gamma0(l1:l2,m1:m2)-1.)*gamma0(l1:l2,m1:m2) &
                   *(drhoE_pU(l1:l2,m1:m2,2))
   
         df(l1:l2,m1:m2,nnn,irho_tmp) = df(l1:l2,m1:m2,nnn,irho_tmp) + T_1_x+T_1_y
         df(l1:l2,m1:m2,nnn,iux) =      df(l1:l2,m1:m2,nnn,iux)      + T_2_x+T_2_y
         df(l1:l2,m1:m2,nnn,iuy) =      df(l1:l2,m1:m2,nnn,iuy)      + T_3_x+T_3_y
         df(l1:l2,m1:m2,nnn,iuz) =      df(l1:l2,m1:m2,nnn,iuz)      + T_4_x+T_4_y
         df(l1:l2,m1:m2,nnn,ilnTT) =    df(l1:l2,m1:m2,nnn,ilnTT)    + T_5_x+T_5_y
  
! if ((nxgrid /= 1) .and. (nygrid /= 1)) then
 ! if ((nxgrid /= 1) .and. (nygrid == 1)) then

!
! Corner points
!
     

       if (nxgrid /= 1) then 
         L_1=(f(l1:l2,m1:m2,nnn,iux) - cs0_ar(l1:l2,m1:m2))&
            *(grad_pp(:,:,1)-rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,1,1))
          L_2=f(l1:l2,m1:m2,nnn,iux)*(cs20_ar(l1:l2,m1:m2) &
                  *grad_rho(:,:,1)-grad_pp(:,:,1))
          L_3=f(l1:l2,m1:m2,nnn,iux)*dui_dxj(:,:,2,1)
          L_4=f(l1:l2,m1:m2,nnn,iux)*dui_dxj(:,:,3,1)
          L_5=(f(l1:l2,m1:m2,nnn,iux) + cs0_ar(l1:l2,m1:m2))&
            *(grad_pp(:,:,1)+ rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,1,1))

         if (x(l1)==xyz0(1)) then
          L_1(nx,:)=KK(nx,:)*(cs20_ar(l2,m1:m2)/gamma0(l2,m1:m2)*&
           rho0(l2,m1:m2)-p_inf(nx,:,nnn-3))
         endif
         if (x(l2)==Lxyz(1)) then
          L_5(1,:)=KK(1,:)*(cs20_ar(l1,m1:m2)/gamma0(l1,m1:m2)*&
            rho0(l1,m1:m2)-p_inf(1,:,nnn-3))
         endif

       else
        L_1=0; L_2=0; L_3=0; L_4=0; L_5=0 
       endif

       if (nygrid /= 1)  then
           M_1=(f(l1:l2,m1:m2,nnn,iuy) - cs0_ar(l1:l2,m1:m2))&
            *(grad_pp(:,:,2)-rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,2,2))
          M_3=f(l1:l2,m1:m2,nnn,iuy)*dui_dxj(:,:,1,2)
          M_2=f(l1:l2,m1:m2,nnn,iuy) &
           *(cs20_ar(l1:l2,m1:m2)*grad_rho(:,:,2)-grad_pp(:,:,2))
          M_4=f(l1:l2,m1:m2,nnn,iuy)*dui_dxj(:,:,3,2)
          M_5=(f(l1:l2,m1:m2,nnn,iuy) + cs0_ar(l1:l2,m1:m2))&
            *(grad_pp(:,:,2)+ rho0(l1:l2,m1:m2)*cs0_ar(l1:l2,m1:m2)*dui_dxj(:,:,2,2))

         if (y(m2)==Lxyz(2)) then
          M_1(:,ny)=KK(:,ny)*(cs20_ar(l1:l2,m2)/gamma0(l1:l2,m2)*&
            rho0(l1:l2,m2)-p_inf(:,ny,nnn-3))
         endif

         if (y(m1)==xyz0(2)) then 
          M_5(:,1)=KK(:,1)*(cs20_ar(l1:l2,m1)/gamma0(l1:l2,m1)*&
            rho0(l1:l2,m1)-p_inf(:,1,nnn-3))
         endif

      else
        M_1=0; M_2=0; M_3=0; M_4=0; M_5=0
       endif

        do j=1,2
         if (j==1) then
           ll=l1; lll=1
           if (x(l1)==xyz0(1))   lcorner_x=.true.
         elseif (j==2) then
           ll=l2; lll=nx
           if (x(l2)==Lxyz(1))   lcorner_x=.true.
         endif
        if (lcorner_x)  then
         df(ll,m1:m2,nnn,irho_tmp) = &
           + drho_prefac(lll,:)*(L_2(lll,:)+0.5*(L_5(lll,:) + L_1(lll,:))) &
           + T_1_y(lll,:) &
           + drho_prefac(lll,:)*(N_2(lll,:)+0.5*(N_5(lll,:) + N_1(lll,:)))
         df(ll,m1:m2,nnn,iux) =  -1./&
           (2.*rho0(ll,m1:m2)*cs0_ar(ll,m1:m2))*(L_5(lll,:) - L_1(lll,:)) &
           +T_2_y(lll,:)-N_3(lll,:)
         df(ll,m1:m2,nnn,iuy) =  - L_3(lll,:) +T_3_y(lll,:)-N_4(lll,:)
         df(ll,m1:m2,nnn,iuz) =  - L_4(lll,:)  + T_4_y(lll,:) - 1./&
            (2.*rho0(ll,m1:m2)*cs0_ar(ll,m1:m2))*(N_5(lll,:) - N_1(lll,:))
         df(ll,m1:m2,nnn,ilnTT) = drho_prefac(lll,:)*(-L_2(lll,:) &
           + 0.5*(gamma0(ll,m1:m2)-1.)*(L_5(lll,:)+L_1(lll,:)))  &
           + T_5_y(lll,:) &
           + drho_prefac(lll,:)*(-N_2(lll,:) &
           + 0.5*(gamma0(ll,m1:m2)-1.)*(N_5(lll,:)+N_1(lll,:)))
         lcorner_x=.false.
       endif
       enddo




        do i=1,2 
         if (i==1) then
          mm=m1; mmm=1
          if (y(m1)==xyz0(2))   lcorner_y=.true.
         elseif (i==2) then
          mm=m2; mmm=ny
          if (y(m2)==Lxyz(2))  lcorner_y=.true.
         endif

        if (lcorner_y)  then
         df(l1:l2,mm,nnn,irho_tmp) = T_1_x(:,mmm)&
          + drho_prefac(:,mmm)*(M_2(:,mmm)+0.5*(M_5(:,mmm) + M_1(:,mmm))) &
          + drho_prefac(:,mmm)*(N_2(:,mmm)+0.5*(N_5(:,mmm) &
          + N_1(:,mmm)))
         df(l1:l2,mm,nnn,iux) =  T_2_x(:,mmm) -M_3(:,mmm)-N_3(:,mmm)
         df(l1:l2,mm,nnn,iuy) =  T_3_x(:,mmm) - 1./&
           (2.*rho0(l1:l2,mm)*cs0_ar(l1:l2,mm))*(M_5(:,mmm) - M_1(:,mmm))-N_4(:,mmm)
         df(l1:l2,mm,nnn,iuz) = T_4_x(:,mmm)- M_4(:,mmm)  - 1./&
             (2.*rho0(l1:l2,mm)*cs0_ar(l1:l2,mm))*(N_5(:,mmm) - N_1(:,mmm))
          df(l1:l2,mm,nnn,ilnTT) =T_5_x(:,mmm) &
           +drho_prefac(:,mmm)*(-M_2(:,mmm) &
           +0.5*(gamma0(l1:l2,mm)-1.)*(M_5(:,mmm)+M_1(:,mmm))) &
           +drho_prefac(:,mmm)*(-N_2(:,mmm) &
           +0.5*(gamma0(l1:l2,mm)-1.)*(N_5(:,mmm)+N_1(:,mmm)))
          lcorner_y=.false.
        endif
       enddo


         do i=1,2 
         if (i==1) then
          mm=m1; mmm=1
          if (y(m1)==xyz0(2))   lcorner_y=.true.
         elseif (i==2) then
          mm=m2; mmm=ny
          if (y(m2)==Lxyz(2))  lcorner_y=.true.
         endif

         do j=1,2
         if (j==1) then
           ll=l1; lll=1
           if (x(l1)==xyz0(1))   lcorner_x=.true.
         elseif (j==2) then
           ll=l2; lll=nx
           if (x(l2)==Lxyz(1))   lcorner_x=.true.
         endif

   !     if ((x(l1)==xyz0(1))   .or. (x(l2)==Lxyz(1)))   lcorner_x=.true.
   !     if ((y(m1)==xyz0(2))   .or. (y(m2)==Lxyz(2)))   lcorner_y=.true.

        if ((lcorner_x) .and. (lcorner_y)) then
         df(ll,mm,nnn,irho_tmp) = &
           drho_prefac(lll,mmm)*(L_2(lll,mmm)+0.5*(L_5(lll,mmm) + L_1(lll,mmm))) &
         + drho_prefac(lll,mmm)*(M_2(lll,mmm)+0.5*(M_5(lll,mmm) + M_1(lll,mmm))) &
         + drho_prefac(lll,mmm)*(N_2(lll,mmm)+0.5*(N_5(lll,mmm) + N_1(lll,mmm)))

         df(ll,mm,nnn,iux) =  -1./&
           (2.*rho0(ll,mm)*cs0_ar(ll,mm))*(L_5(lll,mmm) - L_1(lll,mmm)) &
           -M_3(lll,mmm)-N_3(lll,mmm)
         df(ll,mm,nnn,iuy) =  - L_3(lll,mmm) - 1./&
           (2.*rho0(ll,mm)*cs0_ar(ll,mm))*(M_5(lll,mmm) - M_1(lll,mmm))-N_4(lll,mmm)
         df(ll,mm,nnn,iuz) =  - L_4(lll,mmm) - M_4(lll,mmm)  - 1./&
            (2.*rho0(ll,mm)*cs0_ar(ll,mm))*(N_5(lll,mmm) - N_1(lll,mmm))
         df(ll,mm,nnn,ilnTT) = drho_prefac(lll,mmm)*(-L_2(lll,mmm) &
           +0.5*(gamma0(ll,mm)-1.)*(L_5(lll,mmm)+L_1(lll,mmm)))  &
           +drho_prefac(lll,mmm)*(-M_2(lll,mmm) &
           +0.5*(gamma0(ll,mm)-1.)*(M_5(lll,mmm)+M_1(lll,mmm))) &
           +drho_prefac(lll,mmm)*(-N_2(lll,mmm) &
           +0.5*(gamma0(ll,mm)-1.)*(N_5(lll,mmm)+N_1(lll,mmm)))
          lcorner_x=.false. 
          lcorner_y=.false.
       endif
       enddo
       enddo 

     endif

      if (nchemspec>1) then
       do k=1,nchemspec
         call der_onesided_4_slice(f,sgn,ichemspec(k),dYk_dz(l1:l2,m1:m2),nnn,3)

         do i=m1,m2
          call der_pencil(1,f(:,i,nnn,ichemspec(k)),dYk_dx(:,i))
         enddo
         do i=l1,l2
          call der_pencil(2,f(i,:,nnn,ichemspec(k)),dYk_dy(i,:))
         enddo


          df(l1:l2,m1:m2,nnn,ichemspec(k))=&
              -f(l1:l2,m1:m2,nnn,iux)*dYk_dx(l1:l2,m1:m2) &
             -f(l1:l2,m1:m2,nnn,iuy)*dYk_dy(l1:l2,m1:m2) &
              -f(l1:l2,m1:m2,nnn,iuz)*dYk_dz(l1:l2,m1:m2) &
               +RHS_Y_full(l1:l2,m1:m2,nnn,k)

       ! if (lfilter) then
         do i=l1,l2
         do j=m1,m2
           if ((f(i,j,nnn,ichemspec(k))+df(i,j,nnn,ichemspec(k))*dt)<-1e-25 ) then
             df(i,j,nnn,ichemspec(k))=-1e-25*dt
           endif
           if ((f(i,j,nnn,ichemspec(k))+df(i,j,nnn,ichemspec(k))*dt)>1.) then
             f(i,j,nnn,ichemspec(k))=1.*dt
           endif
         enddo
         enddo
       ! endif
       enddo
      endif
!
!

    endsubroutine bc_nscbc_nref_subout_z
!***********************************************************************
!***********************************************************************
  
    subroutine damp_zone_for_NSCBC(f,df)
!
!   16-jul-06/natalia: coded
!    buffer zone to damp the acustic waves!!!!!!!!!!!
!    important for NSCBC  
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, dimension (mx) :: func_x
      integer :: j1
      real :: dt1, func_y,func_z, ux_ref,uy_ref,uz_ref,lnTT_ref,lnrho_ref
      real :: sz1,sz2, sz1_x,sz2_x, del
      logical :: lzone_y=.false.,lzone_z=.false.

       dt1=1./dt
       del=0.1

       ux_ref=0.
       uy_ref=0.
       uz_ref=0.
     !  lnrho_ref=-7.73236
       lnTT_ref=log(init_TT1)

      lnrho_ref=&
         log(init_pressure)-log(Rgas)-lnTT_ref-log(mu1_full(l1,m1,n1)) 

        do j1=l1,l2
                 
       if (nxgrid/=1) then

        
          
         sz1_x=(xyz0(1)+Lxyz(1)*del)
         sz2_x=(xyz0(1)+Lxyz(1)*(1.-del))
        
!         sz1=(xyz0(2)+Lxyz(2)*del)
!         sz2=(xyz0(2)+Lxyz(2)*(1.-del))


     !     if (sz_r_x<=l1) call fatal_error('to use ldamp_zone_NSCBC',&
     !               'you should increase nxgrid!')

          if (x(j1)<sz1_x) then
           func_x(j1)=(sz1_x-x(j1))**3/(del*Lxyz(1))**3
    
          elseif (sz2_x<x(j1)) then
           func_x(j1)=(x(j1)-sz2_x)**3/(del*Lxyz(1))**3
     
          endif

        if ((x(j1)<sz1_x) .or. (sz2_x<x(j1))) then

          df(j1,m,n,iux)=df(j1,m,n,iux)-func_x(j1)*(f(j1,m,n,iux)-ux_ref)*dt1
          df(j1,m,n,iuy)=df(j1,m,n,iuy)-func_x(j1)*(f(j1,m,n,iuy)-uy_ref)*dt1
          df(j1,m,n,iuz)=df(j1,m,n,iuz)-func_x(j1)*(f(j1,m,n,iuz)-uz_ref)*dt1
          df(j1,m,n,ilnrho)=df(j1,m,n,ilnrho)  &
            -func_x(j1)*(f(j1,m,n,ilnrho)-lnrho_ref)*dt1
          df(j1,m,n,ilnTT)=df(j1,m,n,ilnTT)&  
            -func_x(j1)*(f(j1,m,n,ilnTT)-lnTT_ref)*dt1  

        endif
     
       endif
      
       
       if (nygrid/=1) then

      ! if (sz_r_y<=m1) call fatal_error('to use ldamp_zone_NSCBC',&
      !            'you should increase nygrid!')

          

         sz1=(xyz0(2)+Lxyz(2)*del)
         sz2=(xyz0(2)+Lxyz(2)*(1.-del))

         
       if ((y(m)<=sz1) .and. (y(m)>=xyz0(2))) then
        func_y=(sz1-y(m))**3/(Lxyz(2)*del)**3 
        lzone_y=.true.
       elseif ((y(m)>=sz2) .and. (y(m)<=xyz0(2)+Lxyz(2))) then 
        func_y= (y(m)-sz2)**3/(Lxyz(2)*del)**3 
        lzone_y=.true.
       endif     

       if (lzone_y) then
    
!        if ((x(j1)>sz1_x) .and. (x(j1)<sz2_x)) then
        df(j1,m,n,iux)=df(j1,m,n,iux)-func_y*(f(j1,m,n,iux)-ux_ref)*dt1
        df(j1,m,n,iuy)=df(j1,m,n,iuy)-func_y*(f(j1,m,n,iuy)-uy_ref)*dt1
        df(j1,m,n,iuz)=df(j1,m,n,iuz)-func_y*(f(j1,m,n,iuz)-uz_ref)*dt1
        df(j1,m,n,ilnrho)=df(j1,m,n,ilnrho)&  
           -func_y*(f(j1,m,n,ilnrho)-lnrho_ref)*dt1
        df(j1,m,n,ilnTT)=df(j1,m,n,ilnTT)&  
           -func_y*(f(j1,m,n,ilnTT)-lnTT_ref)*dt1
        endif
  
        lzone_y=.false.
     
!       endif
       endif

      if (nzgrid/=1) then

         sz1=(xyz0(3)+Lxyz(3)*del)
         sz2=(xyz0(3)+Lxyz(3)*(1.-del))

         
       if ((z(n)<=sz1) .and. (z(n)>=xyz0(3))) then
        func_z=(sz1-z(n))**3/(Lxyz(3)*del)**3 
        lzone_z=.true.
       elseif ((z(n)>=sz2) .and. (z(n)<=xyz0(3)+Lxyz(3))) then 
        func_z= (z(n)-sz2)**3/(Lxyz(3)*del)**3 
        lzone_z=.true.
       endif  
    
       if (lzone_z) then
!        if ((x(j1)>sz1_x) .and. (x(j1)<sz2_x)) then
         
        df(j1,m,n,iux)=df(j1,m,n,iux)-func_z*(f(j1,m,n,iux)-ux_ref)*dt1
        df(j1,m,n,iuy)=df(j1,m,n,iuy)-func_z*(f(j1,m,n,iuy)-uy_ref)*dt1
        df(j1,m,n,iuz)=df(j1,m,n,iuz)-func_z*(f(j1,m,n,iuz)-uz_ref)*dt1
        df(j1,m,n,ilnrho)=df(j1,m,n,ilnrho) &
                -func_z*(f(j1,m,n,ilnrho)-lnrho_ref)*dt1
        df(j1,m,n,ilnTT)=df(j1,m,n,ilnTT)  &
                -func_z*(f(j1,m,n,ilnTT)-lnTT_ref)*dt1
        lzone_z=.false.
!        endif
    
       endif
       endif

       enddo

    endsubroutine damp_zone_for_NSCBC
!***********************************************************************
    subroutine jacobn(f,jacob)
! Compute the jacobian, i.e. the matrix  jacob(nchemspec x nchemspec)
! where jacob(i,j)=dv_i/dc_j
! v is the vector of dimension nchemspec of the rates dc_j/dt
! (the c_j being concentrations, stocked in f among other)
!
!  28-may-09/rplasson: coded
!

! Check what is necessary...
!
      use Diagnostics
      use Sub
!
      implicit none
!   exchange data
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,nchemspec,nchemspec) :: jacob
      
      intent(in) :: f
      intent(out) :: jacob
!   internal data
!
!   indices
      integer :: i,j,k,l,ii
!
!   temporary 
      real :: tmp_p,tmp_m
!   Code

      jacob=0.
      
!  identify module 
!
      if (headtt.or.ldebug) print*,'jacobn: compute the jacobian matrix'
!
      if (lreactions) then
        do n=n1,n2; do m=m1,m2;do l=l1,l2
          do i=1,nchemspec 
            do j=1,nchemspec
              ! Compute dv_i/dc_j 
!              print*,"dv_",i,"/dc_",j,"="
              do k=1,nreactions
                ! Check if compound i participate in reaction k
                if (Sijp(i,k)/=Sijm(i,k)) then
                  ! Compute the contribution of reaction k to dv_i/dc_j
!                  print*,"+(",(Sijp(i,k)-Sijm(i,k)),"k_",k,"(+/-)"
                  tmp_p=(Sijp(i,k)-Sijm(i,k))*kreactions_p(k)*kreactions_z(n,k)
                  tmp_m=(Sijm(i,k)-Sijp(i,k))*kreactions_m(k)*kreactions_z(n,k)
                  do  ii=1,nchemspec
                    ! Compute the contribution of compound ii in reaction k
!                    print*,"**-",Sijm(ii,k),Sijp(ii,k),"-**"
                    if (ii/=j) then
!                      print*,"c_",ii,"^",Sijm(ii,k)," (/) ","c_",ii,"^",Sijp(ii,k)
                      tmp_p=tmp_p*f(l,m,n,ichemspec(ii))**Sijm(ii,k)
                      tmp_m=tmp_m*f(l,m,n,ichemspec(ii))**Sijp(ii,k)
                    else
                      if (Sijm(ii,k)==0) then 
!                        print*,"0*c_",ii
                        tmp_p=0.
                      elseif (Sijm(ii,k)>1) then
!                        print*,Sijm(ii,k),"*c_",ii,"^",(Sijm(ii,k)-1)
                        tmp_p=Sijm(ii,k)*tmp_p*f(l,m,n,ichemspec(ii))**(Sijm(ii,k)-1) 
!                      else
!                        print*,"c_",ii,"^0"
                      endif
!                      print*," (/) "
                      if (Sijp(ii,k)==0) then
!                        print*,"0*c_",ii
                        tmp_m=0.
                      elseif (Sijp(ii,k)>1) then
!                        print*,Sijp(ii,k),"*c_",ii,"^",(Sijp(ii,k)-1)
                        tmp_m=Sijp(ii,k)*tmp_m*f(l,m,n,ichemspec(ii))**(Sijp(ii,k)-1) 
!                      else
!                        print*,"c_",ii,"^0"
                      endif
                    endif
                  enddo
!                  print*,")"
                  ! Add the contribution of reaction k to dv_i/dc_j
                  jacob(l,m,n,i,j)=jacob(l,m,n,i,j)+tmp_p+tmp_m
!                  print*,"(=",tmp_p," (-) ",tmp_m,")"
                endif
              enddo
            enddo
          enddo
        enddo; enddo; enddo
      endif
      
    endsubroutine jacobn
!***********************************************************************
    subroutine get_mu1_slicex(slice,grad_slice,index,sgn)
!
! For the NSCBC boudary conditions the slice of mu1 at the boundary
! is required.
!
! 2009.12.10: Nils Erland L. Haugen (coded)
!
      use Deriv, only: der_onesided_4_slice_other
!
      real, dimension(ny,nz), intent(out)   :: slice
      real, dimension(ny,nz), intent(out) :: grad_slice
      integer, intent(in) :: index, sgn
!
      slice=mu1_full(index,m1:m2,n1:n2)
      call der_onesided_4_slice_other(mu1_full,sgn,grad_slice,index,1)
!
    end subroutine get_mu1_slicex
!********************************************************************
  subroutine chemistry_clean_up()
!
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

