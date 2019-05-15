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
! MAUX CONTRIBUTION 2
!
! PENCILS PROVIDED cv; cv1; cp; cp1; glncp(3)
! PENCILS PROVIDED nu; gradnu(3)
! PENCILS PROVIDED DYDt_reac(nchemspec); DYDt_diff(nchemspec)
! PENCILS PROVIDED lambda; glambda(3)
! PENCILS PROVIDED Diff_penc_add(nchemspec); H0_RT(nchemspec); hhk_full(nchemspec)
! PENCILS PROVIDED ghhk(3,nchemspec); S0_R(nchemspec)

!
!***************************************************************
module solid_cells_ogrid_chemistry
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use EquationOfState
  use Chemistry, only: Rgas, lreactions, lchemistry_diag
  use Messages, only: svn_id, timing, fatal_error, inevitably_fatal_error
  use Mpicomm, only: stop_it
  use solid_cells_ogrid_cdata
  use solid_cells_ogrid_sub
!
  implicit none
!
public :: initialize_eos_chemistry, calc_for_chem_mixture_ogrid, calc_pencils_eos_ogrid_chem
public :: calc_pencils_chemistry_ogrid, dYk_dt_ogrid, calc_heter_reaction_term
!
  real :: Rgas_unit_sys=1.
  real, dimension(mx_ogrid,my_ogrid,mz_ogrid,nchemspec) :: cp_R_spec_ogrid
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
!  parameters related to chemical reactions, diffusion and advection
!
  integer, pointer :: nreactions
!  real, dimension(:,:), pointer :: species_constants ! available already from eos_simple
!  real, pointer ::cp_const ! available already from eos_simple
!  real, dimension(:), pointer :: Lewis_coef1
  integer, dimension(:), pointer :: iaa1, iaa2
!  character(len=30), pointer, dimension(:) :: reaction_name
  real, dimension(:), pointer :: alpha_n, E_an, B_n
  logical, pointer ::  lcheminp,ldiffusion,ldiff_corr, lmech_simple
  logical, pointer ::  tran_exist, lThCond_simple !,lheatc_chemistry,
  logical, pointer :: lfilter_strict, lfilter, ladvection,lt_const
  real, pointer, dimension(:,:) :: tran_data, Sijm, Sijp, stoichio
  real, pointer, dimension(:,:) :: low_coeff, high_coeff, troe_coeff, a_k4
  logical, pointer, dimension(:) :: photochem_case, Mplus_case, back
!
!  character(len=30), allocatable, dimension(:) :: reaction_name
  logical :: lT_tanh=.false.
  logical :: linit_temperature=.false., linit_density=.false.
  logical :: lreac_as_aux=.false.
!
  logical :: lflame_front=.false.
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
  real, dimension(nchemspec) :: nu_spec=0., mobility=1.
!
!  Chemkin related parameters
!
  logical :: lchem_cdtc=.false.
  logical :: lmobility=.false.
  integer :: iTemp1=2, iTemp2=3, iTemp3=4
  real :: lamb_low, lamb_up, Pr_turb=0.7
!
! carbon molar mass
     real :: M_C=12.0107
!
!   Diagnostics
!
!  real, allocatable, dimension(:,:) :: net_react_m, net_react_p
!  logical :: lchemistry_diag=.false.
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
  integer :: ireac=0
  integer, dimension(nchemspec) :: ireaci=0
!
!  Indices
!
  integer :: i_N2=0, i_CO2=0, i_CO=0
  integer :: ichem_O2=0, ichem_N2=0, ichem_CO2=0, ichem_CO=0
  logical :: lO2=.false., lN2=.false.
  logical :: lCO2=.false., lCO=.false.
!
  contains
!
!***********************************************************************
!    subroutine register_chemistry()
!
! Not used for ogrid
!
!    endsubroutine register_chemistry
!***********************************************************************
!    subroutine read_thermodyn_simple()
!
! species_constants - available from eos, iaa1, iaa2 - copied from chemistry
!
!    endsubroutine read_thermodyn_simple
!***********************************************************************
    subroutine initialize_eos_chemistry
!
! Initialize variable selection code (needed for RELOADing)
!
        ieosvars=-1
        ieosvar_count=0
!
        if ((nxgrid_ogrid==1) .and. (nygrid_ogrid==1) .and. (nzgrid_ogrid==1)) then
          ll1_ogrid=1; ll2_ogrid=mx_ogrid; mm1_ogrid=m1_ogrid; mm2_ogrid=m2_ogrid; nn1_ogrid=n1_ogrid; nn2_ogrid=n2_ogrid
        elseif (nxgrid_ogrid==1) then
          ll1_ogrid=l1_ogrid; ll2_ogrid=l2_ogrid
        else
          ll1_ogrid=1; ll2_ogrid=mx_ogrid
        endif
!
        if (nygrid_ogrid==1) then
          mm1_ogrid=m1_ogrid; mm2_ogrid=m2_ogrid
        else
          mm1_ogrid=1; mm2_ogrid=my_ogrid
        endif
!
        if (nzgrid_ogrid==1) then
          nn1_ogrid=n1_ogrid; nn2_ogrid=n2_ogrid
        else
          nn1_ogrid=1;  nn2_ogrid=mz_ogrid
        endif
!
!      if (.not.ldensity) then
!        call put_shared_variable('rho0',rho0,caller='initialize_eos')
!        call put_shared_variable('lnrho0',lnrho0)
!      endif
!
    endsubroutine initialize_eos_chemistry
!***********************************************************************
    subroutine initialize_chemistry_og(f_og)
!
!  chemistry initialization on the ogrid
!  called from initialize_solid_cells in solid_cells_ogrid
!
      use SharedVariables, only: get_shared_variable, put_shared_variable
!
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f_og
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid) :: mu1_full_og
      integer :: i

!  calculate universal gas constant based on Boltzmann constant
!  and the proton mass
!
        if (unit_system == 'cgs') then
          Rgas_unit_sys = k_B_cgs/m_u_cgs
          Rgas = Rgas_unit_sys/unit_energy
        endif
!
      if ((nxgrid_ogrid == 1) .and. (nygrid_ogrid == 1) .and. (nzgrid_ogrid == 1)) then
        ll1_ogrid = 1
        ll2_ogrid = mx_ogrid
        mm1_ogrid = m1_ogrid
        mm2_ogrid = m2_ogrid
        nn1_ogrid = n1_ogrid
        nn2_ogrid = n2_ogrid
      else
        if (nxgrid_ogrid == 1) then
          ll1_ogrid = l1_ogrid
          ll2_ogrid = l2_ogrid
        else
          ll1_ogrid = 1
          ll2_ogrid = mx_ogrid
        endif
!
        if (nygrid_ogrid == 1) then
          mm1_ogrid = m1_ogrid
          mm2_ogrid = m2_ogrid
        else
          mm1_ogrid = 1
          mm2_ogrid = my_ogrid
        endif
!
        if (nzgrid_ogrid == 1) then
          nn1_ogrid = n1_ogrid
          nn2_ogrid = n2_ogrid
        else
          nn1_ogrid = 1
          nn2_ogrid = mz_ogrid
        endif
      endif
!
!  Reinitialize if required
!
      if (reinitialize_chemistry) then
          call fatal_error('initialize_chemistry_og', &
              'Reinitialize chemistry - not implemented on the ogrid!')
      endif
!
      call getmu_array_ogrid(f_og,mu1_full_og)
      f_og(:,:,:,iRR) = mu1_full_og*Rgas
!
!  allocate memory for net_reaction diagnostics
!
!      if (allocated(net_react_p) .and. .not. lreloading) then
!        print*, 'this should not be here'
!      endif
!
!      if (lchemistry_diag .and. .not. lreloading) then
!        allocate(net_react_p(nchemspec,nreactions),STAT=stat)
!        if (stat > 0) call stop_it("Couldn't allocate memory for net_react_p")
!        net_react_p = 0.
!        allocate(net_react_m(nchemspec,nreactions),STAT=stat)
!        if (stat > 0) call stop_it("Couldn't allocate memory for net_react_m")
!        net_react_m = 0.
!      endif
!
!  Define the chemical reaction rates as auxiliary variables for output
!
!      if (lreac_as_aux) then
!        if (ireac == 0) then
!          call farray_register_auxiliary('reac',ireac,vector=nchemspec)
!          do i = 0, nchemspec-1
!            ireaci(i+1) = ireac+i
!          enddo
!        endif
!        if (ireac /= 0 .and. lroot) then
!          print*, 'initialize_reaction_rates: ireac = ', ireac
!          open (3,file=trim(datadir)//'/index.pro', POSITION='append')
!          write (3,*) 'ireac=',ireac
!          do i = 1, nchemspec
!            write (car2,'(i2)') i
!            write (3,*) 'ireac'//trim(adjustl(car2))//'=', ireaci(i)
!          enddo
!          close (3)
!        endif
!      endif
!
        call get_shared_variable('ldiffusion',ldiffusion)
        call get_shared_variable('ldiff_corr',ldiff_corr)
        call get_shared_variable('lcheminp',lcheminp)
        call get_shared_variable('lThCond_simple',lThCond_simple)
        call get_shared_variable('tran_exist',tran_exist)
        call get_shared_variable('tran_data',tran_data)
        call get_shared_variable('lt_const',lt_const)
        call get_shared_variable('ladvection',ladvection)
        call get_shared_variable('lfilter',lfilter)
        call get_shared_variable('lfilter_strict',lfilter_strict)
        call get_shared_variable('Lewis_coef1',Lewis_coef1)
        call get_shared_variable('nreactions',nreactions)
        call get_shared_variable('Sijm',Sijm)
        call get_shared_variable('Sijp',Sijp)
  !      call get_shared_variable('reaction_name',reaction_name)
        call get_shared_variable('back',back)
        call get_shared_variable('B_n',B_n)
        call get_shared_variable('alpha_n',alpha_n)
        call get_shared_variable('E_an',E_an)
        call get_shared_variable('low_coeff',low_coeff)
        call get_shared_variable('high_coeff',high_coeff)
        call get_shared_variable('troe_coeff',troe_coeff)
        call get_shared_variable('a_k4',a_k4)
        call get_shared_variable('Mplus_case',Mplus_case)
        call get_shared_variable('photochem_case',photochem_case)
        call get_shared_variable('stoichio',stoichio)
        call get_shared_variable('iaa1',iaa1)
        call get_shared_variable('iaa2',iaa2)
        call get_shared_variable('lmech_simple',lmech_simple)
        call get_shared_variable('linterp_pressure',linterp_pressure)
!
        if (lreac_heter) then
          if (lmech_simple) then
            call find_species_index('CO2',i_CO2,ichem_CO2,lCO2)
            call find_species_index('CO',i_CO,ichem_CO,lCO)
            call find_species_index('N2',i_N2,ichem_N2,lN2)
            call find_species_index('O2',i_O2,ichem_O2,lO2)
            allocate(heter_reaction_rate(my_ogrid,mz_ogrid,nchemspec+1))
            heter_reaction_rate(:,:,:)=0
          else
            call fatal_error('initialize_chemistry_og', &
               'Heterogeneous chemistry works only with simplified mechanism (lmech_simple=T)')
          endif
        endif
!
!  write array dimension to chemistry diagnostics file
!
!      open (1,file=trim(datadir)//'/net_reactions.dat',position='append')
!      write (1,*) nchemspec,nreactions
!      close (1)
!
!
    endsubroutine initialize_chemistry_og
!***********************************************************************
!    subroutine init_chemistry_og(f_og)
!
! TODO: decide what to do for the ogrid
! Not used for ogrid
!
!    endsubroutine init_chemistry_og
!***********************************************************************
!    subroutine pencil_criteria_chemistry()
!
! Which pencils to cpmute is specified in solid_cells_ogrid
!
!    endsubroutine pencil_criteria_chemistry
!***********************************************************************
!    subroutine pencil_interdep_chemistry(lpencil_in)
!
! Which pencils to cpmute is specified in solid_cells_ogrid
!
!    endsubroutine pencil_interdep_chemistry
!***********************************************************************
    subroutine calc_pencils_chemistry_ogrid(f_og)
!
!  Calculate chemistry pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
!
      integer :: k,i,j2,j3
      integer :: ii1=1, ii2=2, ii3=3, ii4=4, ii5=5, ii6=6, ii7=7
      real :: T_low, T_up, T_mid
      real, dimension(nx_ogrid) :: T_loc, TT_2, TT_3, TT_4, D_th
!
      T_loc = p_ogrid%TT
      TT_2 = p_ogrid%TT*p_ogrid%TT
      TT_3 = TT_2*p_ogrid%TT
      TT_4 = TT_2*TT_2
!
!  Dimensionless Standard-state molar enthalpy H0/RT
!
      if (lpencil_ogrid(i_og_H0_RT)) then
        if (.not. lT_const) then
          do k = 1,nchemspec
            T_low = species_constants(k,iTemp1)
            T_mid = species_constants(k,iTemp2)
            T_up = species_constants(k,iTemp3)
            where (T_loc <= T_mid)
              p_ogrid%H0_RT(:,k) = species_constants(k,iaa2(ii1))         &
                                 + species_constants(k,iaa2(ii2))*T_loc/2 &
                                 + species_constants(k,iaa2(ii3))*TT_2/3  &
                                 + species_constants(k,iaa2(ii4))*TT_3/4  &
                                 + species_constants(k,iaa2(ii5))*TT_4/5  &
                                 + species_constants(k,iaa2(ii6))/T_loc
            elsewhere
              p_ogrid%H0_RT(:,k) = species_constants(k,iaa1(ii1))         &
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
          if (lpencil_ogrid(i_og_hhk_full) ) then
            do k = 1,nchemspec
              if (species_constants(k,imass) > 0.)  then
                p_ogrid%hhk_full(:,k) = p_ogrid%H0_RT(:,k)*Rgas*T_loc &
                    /species_constants(k,imass)
              endif
            enddo
          endif
        endif
      endif
!
      if (lpencil_ogrid(i_og_ghhk)) then
        do k = 1,nchemspec
          if (species_constants(k,imass) > 0.)  then
            do i = 1,3
              p_ogrid%ghhk(:,i,k) = cp_R_spec_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,k)&
                                  /species_constants(k,imass)*Rgas*p_ogrid%gTT(:,i)   
            enddo
          endif
        enddo
      endif
!
!  Find the entropy by using fifth order temperature fitting function
!
      if (lpencil_ogrid(i_og_S0_R)) then
        if (.not. lT_const) then
          do k = 1,nchemspec
            T_low = species_constants(k,iTemp1)
            T_mid = species_constants(k,iTemp2)
            T_up = species_constants(k,iTemp3)
            where (T_loc <= T_mid .and. T_low <= T_loc)
              p_ogrid%S0_R(:,k) = species_constants(k,iaa2(ii1))*p_ogrid%lnTT &
                                + species_constants(k,iaa2(ii2))*T_loc  &
                                + species_constants(k,iaa2(ii3))*TT_2/2 &
                                + species_constants(k,iaa2(ii4))*TT_3/3 &
                                + species_constants(k,iaa2(ii5))*TT_4/4 &
                                + species_constants(k,iaa2(ii7))
            elsewhere (T_mid <= T_loc .and. T_loc <= T_up)
              p_ogrid%S0_R(:,k) = species_constants(k,iaa1(ii1))*p_ogrid%lnTT &
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
      if (lreactions .and. lpencil_ogrid(i_og_DYDt_reac)) then
        call calc_reaction_term_ogr(f_og,p_ogrid)
      else
        p_ogrid%DYDt_reac = 0.
      endif
!
! Calculate the diffusion term and the corresponding pencil
!
      if (lpencil_ogrid(i_og_Diff_penc_add)) then
!
! D_th is computed twice in calc_penc in eos and here
!
        D_th = f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iviscosity)*Pr_number1
        do k = 1,nchemspec
          p_ogrid%Diff_penc_add(:,k) = D_th*Lewis_coef1(k)
        enddo
      endif
!
      if (ldiffusion .and. lpencil_ogrid(i_og_DYDt_diff)) then
        call calc_diffusion_term_ogrid(f_og,p_ogrid)
      else
        p_ogrid%DYDt_diff = 0.
      endif
!
    endsubroutine calc_pencils_chemistry_ogrid
!***********************************************************************
    subroutine calc_pencils_eos_ogrid_chem(f_og)
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f_og

      real, dimension(nx_ogrid,3) :: glnDiff_full_add, glncp
      real, dimension(nx_ogrid) :: D_th, R_mix
!
      integer :: i,k,j2,j3
!
! Cp/Cv pencils
!
      if (lpencil_ogrid(i_og_cp)) p_ogrid%cp = f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,icp)
!
      if (lpencil_ogrid(i_og_cv))  then
        p_ogrid%cv = 0.
        if (Cp_const < impossible) then
          do k = 1,nchemspec
            p_ogrid%cv = p_ogrid%cv + (Cp_const-Rgas)/species_constants(k,imass)&
                         *f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))
          enddo
        else
          R_mix = 0.
          do k = 1,nchemspec
            R_mix = R_mix + Rgas/species_constants(k,imass)&
                    *f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))
          enddo
          p_ogrid%cv = p_ogrid%cp - R_mix 
        endif
        if (lpencil_ogrid(i_og_cv1)) p_ogrid%cv1 = 1./p_ogrid%cv
        if (lpencil_ogrid(i_og_cp1)) p_ogrid%cp1 = 1./p_ogrid%cp
      endif
!
!  Temperature
!
      if (lpencil_ogrid(i_og_lnTT)) then
        if (ltemperature_nolog) then
          p_ogrid%lnTT=log(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT))
        else
          p_ogrid%lnTT=f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ilnTT)
        endif
      endif
!
      if (lpencil_ogrid(i_og_TT))  then
        if (ltemperature_nolog) then
          p_ogrid%TT=f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT)
        else
          p_ogrid%TT=exp(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ilnTT))
        endif
!
        if (minval(p_ogrid%TT)==0.) then
          call fatal_error('calc_pencils_eos_ogrid','p_ogrid%TT=0!')
        endif         
      endif
!
      if (lpencil_ogrid(i_og_TT1)) then
        if (ltemperature_nolog) then
          p_ogrid%TT1=1./f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT)
        else 
          p_ogrid%TT1=1./exp(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ilnTT))
        endif
      endif
!
!  Temperature laplacian and gradient
!
      if (lpencil_ogrid(i_og_gTT))call grad_ogrid(f_og,iTT,p_ogrid%gTT)
      if (lpencil_ogrid(i_og_glnTT)) then
        if (ltemperature_nolog) then
          p_ogrid%glnTT(:,1)=p_ogrid%gTT(:,1)/p_ogrid%TT(:)
          p_ogrid%glnTT(:,2)=p_ogrid%gTT(:,2)/p_ogrid%TT(:)
          p_ogrid%glnTT(:,3)=p_ogrid%gTT(:,3)/p_ogrid%TT(:)
        else
          call grad_ogrid(f_og,ilnTT,p_ogrid%glnTT)
        endif
      endif
!
      if (lpencil_ogrid(i_og_del2lnTT) .or. lpencil_ogrid(i_og_del2TT)) then
        if (ltemperature_nolog) then
          call del2_ogrid(f_og   ,iTT,p_ogrid%del2TT)
        else
          call del2_ogrid(f_og,ilnTT,p_ogrid%del2lnTT)
        endif
      endif
!
! Viscosity of a mixture
!
      if (lpencil_ogrid(i_og_nu)) then
          p_ogrid%nu = f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iviscosity)
        if (lpencil_ogrid(i_og_gradnu)) then
          call grad_ogrid(f_og,iviscosity,p_ogrid%gradnu)
        endif
      endif
!
! Calculate thermal conductivity & diffusivity
!
      D_th = f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iviscosity)*Pr_number1
      if (lpencil_ogrid(i_og_lambda)) p_ogrid%lambda = p_ogrid%cp*p_ogrid%rho*D_th
      if (lpencil_ogrid(i_og_glambda)) then         
        call grad_ogrid(f_og,icp,p_ogrid%glncp)
        do i = 1,3
          p_ogrid%glncp(:,i) = p_ogrid%glncp(:,i)*p_ogrid%cp1
          glnDiff_full_add(:,i) = p_ogrid%gradnu(:,i)/p_ogrid%nu
          p_ogrid%glambda(:,i) = p_ogrid%lambda*(p_ogrid%glnrho(:,i) &
                               + glnDiff_full_add(:,i) + p_ogrid%glncp(:,i)) 
        enddo
      endif
!
!  Mean molecular weight
!
      if (lpencil_ogrid(i_og_RR)) p_ogrid%RR=f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iRR)
!
      if (lpencil_ogrid(i_og_glnRR)) call grad_ogrid(f_og,iRR,p_ogrid%glnRR)
      do i = 1,3
        p_ogrid%glnRR(:,i) = p_ogrid%glnRR(:,i)/p_ogrid%RR
      enddo
!
!  Pressure
!
      if (lpencil_ogrid(i_og_pp)) p_ogrid%pp = p_ogrid%TT*p_ogrid%RR*p_ogrid%rho
      if (lpencil_ogrid(i_og_cs2)) then
        if (any(p_ogrid%cv1 == 0.0)) then
        else
          p_ogrid%cs2 = p_ogrid%cp*p_ogrid%cv1*p_ogrid%RR*p_ogrid%TT
        endif
      endif
!
!  pressure gradient
!
      if (lpencil_ogrid(i_og_rho1gpp)) then
        do i=1,3
          p_ogrid%fpres(:,i) = -p_ogrid%pp*p_ogrid%rho1 &
              *(p_ogrid%glnrho(:,i)+p_ogrid%glnTT(:,i)+p_ogrid%glnRR(:,i))
!              *(p_ogrid%glnrho(:,i)+p_ogrid%glnTT(:,i)+p_ogrid%gmu1(:,i)/p_ogrid%mu1)
 ! This actually can be computed without log gradients
 !         p_ogrid%rho1gpp(:,i) = Rgas*(p_ogrid%mu1*p_ogrid%TT*p_ogrid%rho1*p_ogrid%grho(:,i)&
 !                              + p_ogrid%mu1*p_ogrid%gTT(:,i)+p_ogrid%TT*p_ogrid%gmu1(:,i))
 ! Or
 !         p_ogrid%rho1gpp(:,i) = p_ogrid%cv*p_ogrid%cp1*p_ogrid%cs2*&
 !               (p_ogrid%glnrho(:,i)+p_ogrid%gTT(:,i)/p_ogrid%TT+p_ogrid%gmu1(:,i)/p_ogrid%mu1)
        enddo
      endif
!
      if (lpres_grad) then
        f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,igpx) = -p_ogrid%fpres(:,1)*p_ogrid%rho
        f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,igpy) = -p_ogrid%fpres(:,2)*p_ogrid%rho
      endif
      if (linterp_pressure) then
        f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ipp ) =  p_ogrid%pp
      endif
!
!  Energy per unit mass 
!
!        if (lpencil_ogrid(i_og_ee)) p_ogrid%ee = p_ogrid%cv*p_ogrid%TT
!
! Gradient of lnpp (removed)
!
    endsubroutine calc_pencils_eos_ogrid_chem
!***********************************************************************
!    subroutine flame_front(f_og)
!
!    endsubroutine flame_front
!***********************************************************************
    subroutine calc_for_chem_mixture_ogrid(f_og)
!
!  Calculate quantities for a mixture
!
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f_og
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid) :: mu1_full_og
      real, dimension(mx_ogrid) ::  nu_dyn
! lambda_Suth in cgs [g/(cm*s*K^0.5)]
      real :: lambda_Suth = 1.5e-5, Suth_const = 200.
      integer :: j2,j3, k
      real, dimension(mx_ogrid) :: T_loc, T_loc_2, T_loc_3, T_loc_4
      real :: T_up, T_mid, T_low
!
            call getmu_array_ogrid(f_og,mu1_full_og)
            f_og(:,:,:,iRR) = mu1_full_og*Rgas
!
            f_og(:,:,:,icp) = 0.
            if (Cp_const < impossible) then
!
              f_og(:,:,:,icp) = Cp_const*mu1_full_og
!            
            else
!
              do j3 = nn1_ogrid,nn2_ogrid
                do j2 = mm1_ogrid,mm2_ogrid
                  T_loc = f_og(:,j2,j3,iTT)
!                  T_loc = exp(f_og(:,j2,j3,ilnTT))
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
                        cp_R_spec_ogrid(:,j2,j3,k) = species_constants(k,iaa2(1)) &
          !              cp_R_spec_ogrid(:) = species_constants(k,iaa2(1)) &
                          +species_constants(k,iaa2(2))*T_loc &
                          +species_constants(k,iaa2(3))*T_loc_2 &
                          +species_constants(k,iaa2(4))*T_loc_3 &
                          +species_constants(k,iaa2(5))*T_loc_4
                      elsewhere (T_loc >= T_mid .and. T_loc <= T_up)
                        cp_R_spec_ogrid(:,j2,j3,k) = species_constants(k,iaa1(1)) &
          !              cp_R_spec_ogrid(:) = species_constants(k,iaa1(1)) &
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
                      !  call inevitably_fatal_error('calc_for_chem_mixture_ogrid', &
                      !  'TT_full(:,j2,j3) is outside range', .true.)
                      endif
!
! Find cp and cv for the mixture for the full domain
!
                      f_og(:,j2,j3,icp) = f_og(:,j2,j3,icp)+cp_R_spec_ogrid(:,j2,j3,k) &
                                        * Rgas/species_constants(k,imass)*f_og(:,j2,j3,ichemspec(k))
                    endif
                  enddo
                enddo
              enddo
            endif
!
!  Viscosity of a mixture
!
      do j3 = nn1_ogrid,nn2_ogrid
        do j2 = mm1_ogrid,mm2_ogrid
          if (ltemperature_nolog) then
            nu_dyn = lambda_Suth*f_og(:,j2,j3,iTT)**(3./2.)/(Suth_const+f_og(:,j2,j3,iTT))
          else
            nu_dyn = lambda_Suth*exp(f_og(:,j2,j3,ilnTT))**(3./2.)/(Suth_const+exp(f_og(:,j2,j3,ilnTT)))
          endif
          if (ldensity_nolog) then
            f_og(:,j2,j3,iviscosity) = nu_dyn/f_og(:,j2,j3,irho)
          else
            f_og(:,j2,j3,iviscosity) = nu_dyn/exp(f_og(:,j2,j3,ilnrho))
          endif
        enddo
      enddo
!
    endsubroutine calc_for_chem_mixture_ogrid
!***********************************************************************
    subroutine dYk_dt_ogrid(f_og,df,dt_ogrid)
!
!  species transport equation
!  calculate dYk/dt = - u.gradYk - div(rho*Yk*Vk)/rho + R/rho
!  add chemistry contribution to temperature equation
!
      use Diagnostics
!
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) ::  f_og
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df
      real, dimension(nx_ogrid,3) :: gchemspec, dk_D, sum_diff=0.
      real, dimension(nx_ogrid) :: ugchemspec, sum_DYDT, sum_dhhk=0.
      real, dimension(nx_ogrid) :: sum_dk_ghk, dk_dhhk, sum_hhk_DYDt_reac
      real, dimension(nx_ogrid) :: RHS_T_full!, diffus_chem

      integer :: j,k,i
      real :: dt_ogrid

      intent(inout) :: df, f_og
!
!  identify module and boundary conditions
!
      call timing('dchemistry_dt','entered',mnloop=.true.)
      if (headtt .or. ldebug) print*,'dchemistry_dt: SOLVE dchemistry_dt'
        do k=1,nchemspec
!
!  Advection 
!
          if (lhydro .and. ladvection) then
            call grad_ogrid(f_og,ichemspec(k),gchemspec)
            call dot_mn_ogrid(p_ogrid%uu,gchemspec,ugchemspec)
!            if (lmobility) ugchemspec = ugchemspec*mobility(k)
            df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k)) = &
                df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k)) - ugchemspec
          endif
!
!  Diffusion
!
          if (ldiffusion) then
            df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k)) = &
                df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k)) + p_ogrid%DYDt_diff(:,k)
          endif
!
!  Reaction
!
          if (lreactions) then
              df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k)) = &
                df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k)) + p_ogrid%DYDt_reac(:,k)
          endif
!
!  Add filter for negative concentrations
! TODO: dt or dt_ogrid or it doesn't matter? 
!
          if (lfilter .and. .not. lfilter_strict) then
            do i = 1,mx_ogrid
              if ((f_og(i,m_ogrid,n_ogrid,ichemspec(k))&
 !                + df(i,m_ogrid,n_ogrid,ichemspec(k))*dt) < -1e-25 ) then
 !               df(i,m_ogrid,n_ogrid,ichemspec(k)) = -1e-25*dt
                 + df(i,m_ogrid,n_ogrid,ichemspec(k))*dt_ogrid) < -1e-25 ) then
                df(i,m_ogrid,n_ogrid,ichemspec(k)) = -1e-25*dt_ogrid
              endif
              if ((f_og(i,m_ogrid,n_ogrid,ichemspec(k))&
      !           + df(i,m_ogrid,n_ogrid,ichemspec(k))*dt) > 1. ) then
      !          df(i,m_ogrid,n_ogrid,ichemspec(k)) = 1.*dt
                 + df(i,m_ogrid,n_ogrid,ichemspec(k))*dt_ogrid) > 1. ) then
                df(i,m_ogrid,n_ogrid,ichemspec(k)) = 1.*dt_ogrid
              endif
            enddo
          endif
!
!  Add strict filter for negative concentrations
!
          if (lfilter_strict) then
            do i = 1,mx_ogrid
              if ((f_og(i,m_ogrid,n_ogrid,ichemspec(k))&
                 + df(i,m_ogrid,n_ogrid,ichemspec(k))*dt) < 0.0 ) then
                if (df(i,m_ogrid,n_ogrid,ichemspec(k)) < 0.)&
                   df(i,m_ogrid,n_ogrid,ichemspec(k)) = 0.
              endif
              if ((f_og(i,m_ogrid,n_ogrid,ichemspec(k))&
                 + df(i,m_ogrid,n_ogrid,ichemspec(k))*dt) > 1. ) then
                df(i,m_ogrid,n_ogrid,ichemspec(k)) = 1.*dt
              endif
            enddo
          endif

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
                  *(p_ogrid%DYDt_reac(:,k)+p_ogrid%DYDt_diff(:,k))
              if (lreactions) then
                sum_hhk_DYDt_reac = sum_hhk_DYDt_reac-p_ogrid%hhk_full(:,k)*p_ogrid%DYDt_reac(:,k)
              endif
!
!  Sum over all species of diffusion terms
!
              if (ldiffusion) then
                do i = 1,3
                  dk_D(:,i) = gchemspec(:,i)*p_ogrid%Diff_penc_add(:,k)
                enddo
                call dot_mn_ogrid(dk_D,p_ogrid%ghhk(:,:,k),dk_dhhk)
                sum_dk_ghk = sum_dk_ghk+dk_dhhk
                if (ldiff_corr) sum_diff(:,k) = sum_diff(:,k)+dk_D(:,k)
              endif
!
            endif
          enddo
!
! If the correction velocity is added
!
          if (ldiff_corr .and. ldiffusion) then
            do k = 1,nchemspec
              call dot_mn_ogrid(sum_diff,p_ogrid%ghhk(:,:,k),sum_dhhk)
              sum_dk_ghk(:) = sum_dk_ghk(:)-f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))*sum_dhhk(:)
            enddo
          endif
!
        if (ltemperature_nolog) then
          RHS_T_full = p_ogrid%cv1*((sum_DYDt(:)-p_ogrid%RR*p_ogrid%divu)*p_ogrid%TT &
                     + sum_dk_ghk+sum_hhk_DYDt_reac)
        else
          RHS_T_full = (sum_DYDt(:)-p_ogrid%RR*p_ogrid%divu)*p_ogrid%cv1 &
            + sum_dk_ghk*p_ogrid%TT1(:)*p_ogrid%cv1+sum_hhk_DYDt_reac*p_ogrid%TT1(:)*p_ogrid%cv1
        endif
!
        if (.not. lT_const) then
            df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ilnTT) = &
               df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ilnTT) + RHS_T_full
        endif
!
!     Heatcond called from denergy_dt so removed here
!        if (lheatc_chemistry) call calc_heatcond_chemistry_ogrid(f_og,df,p_ogrid)
!
      endif
!
!      if (lreactions .and. ireac /= 0 .and. ((.not. llsode))) then
!        if (llast) call get_reac_rate(sum_hhk_DYDt_reac,f_og,p_ogrid)
!      endif
!
!
! The part below is commented because ldt = false all the time and it's
! never executed (llsode = F, lfirst = T/F depending on sub-timestep)
! TODO: uncomment
!
! this damping zone is needed in a case of NSCBC
!
!      if (ldamp_zone_for_NSCBC) call damp_zone_for_NSCBC(f,df)
!
!  For the timestep calculation, need maximum diffusion
!
!      if (lfirst .and. ldt) then
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
    endsubroutine dYk_dt_ogrid
!***********************************************************************
    subroutine rprint_chemistry_og(lreset,lwrite)
!
!  reads and registers print parameters relevant to chemistry
!
!  13-aug-07/steveb: coded
!
      use Diagnostics, only: parse_name
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
!        call parse_name(iname,cname(iname),cform(iname),'cp1m',idiag_cp1m)
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
      enddo
!
!  Write chemistry index in short notation
!
      if (lwr) then
        write (3,*) 'nchemspec=',nchemspec
        write (3,*) 'ichemspec=indgen('//trim(itoa(nchemspec))//') + '//trim(itoa(ichemspec(1)))
      endif
!
    endsubroutine rprint_chemistry_og
!***********************************************************************
!    subroutine get_slices_chemistry(f,slices)
!
! TODO: do I need that for the ogrid?
!
!  Write slices for animation of Chemistry variables.
!
!  13-aug-07/steveb: dummy
!  16-may-09/raphael: added more slices
!
!      use Slices_methods, only: assign_slices_scal
!
!      real, dimension(mx,my,mz,mfarray) :: f
!      type (slice_data) :: slices
!
!      character(len=fmtlen) :: sname
!      integer :: ispec
!
!  Chemical species mass fractions.
!
!      sname=trim(slices%name)
!      if (sname(9:)==' ') then    ! 9=len('chemspec')+1
!        ispec=1
!      else
!        read(sname(9:),'(i3)') ispec
!      endif
! 
!      call assign_slices_scal(slices,f,ichemspec(ispec))
!
!    endsubroutine get_slices_chemistry
!***********************************************************************
!    subroutine build_stoich_matrix(StartInd,StopInd,k,ChemInpLine,product)
!
! Does not involve f array, not computed separately for the ogrid
!
!    endsubroutine build_stoich_matrix
!***********************************************************************
! TODO: do I need that for the ogrid?
!
!    subroutine write_reactions()
!
!  write reaction coefficient in the output file
!
!  11-mar-08/nils: coded
!
!      use General, only: itoa
!
!      integer :: reac, spec
!      character(len=80) :: reac_string, product_string, output_string
!      character(len=intlen) :: Sijp_string, Sijm_string
!      character(len=1) :: separatorp, separatorm
!      character(len=fnlen) :: input_file="./data/chem.out"
!      integer :: file_id=123
!
!      open (file_id,file=input_file,POSITION='APPEND',FORM='FORMATTED')
!      write (file_id,*) 'REACTIONS'
      !open(file_id,file=input_file)
!
!      do reac = 1,mreactions
!        reac_string = ''
!        product_string = ''
!        separatorp = ''
!        separatorm = ''
!
!        do spec = 1,nchemspec
!          if (Sijp(spec,reac) > 0) then
!            Sijp_string = ''
!            if (Sijp(spec,reac) /= 1) then
!              write (Sijp_string,'(F3.1)') Sijp(spec,reac)
!            endif
!            if (Sijp(spec,reac)>1) Sijp_string=itoa(Sijp(spec,reac))
!            reac_string = trim(reac_string)//trim(separatorp)// &
!                trim(Sijp_string)//trim(varname(ichemspec(spec)))
!            separatorp = '+'
!          endif
!          if (Sijm(spec,reac) > 0) then
!            Sijm_string = ''
!            if (Sijm(spec,reac) /= 1) write (Sijm_string,'(F3.1)') Sijm(spec,reac)
!            if (Sijm(spec,reac)>1) Sijm_string=itoa(Sijm(spec,reac))
!            product_string = trim(product_string)//trim(separatorm)// &
!                trim(Sijm_string)//trim(varname(ichemspec(spec)))
!            separatorm = '+'
!          endif
!        enddo
!
!        output_string = trim(reac_string)//'='//trim(product_string)
!
!        if (.not. photochem_case(reac)) then
!
!  Note that since the B_n term is in logarithmic form within the code
!  the exponential must be used for output.
!
!          write (unit=output_string(30:45),fmt='(E14.4)') exp(B_n(reac))
!          write (unit=output_string(47:62),fmt='(E14.4)') alpha_n(reac)
!          write (unit=output_string(64:79),fmt='(E14.4)') E_an(reac)
!        endif
!        write (file_id,*) trim(output_string)
!        if (.not. photochem_case(reac)) then
!          if (maxval(abs(low_coeff(:,reac))) > 0.) then
!            write (file_id,*) 'LOW/',exp(low_coeff(1,reac)),low_coeff(2:,reac)
!          elseif (maxval(abs(high_coeff(:,reac))) > 0.) then
!            write (file_id,*) 'HIGH/',high_coeff(:,reac)
!          endif
!          if (maxval(abs(troe_coeff(:,reac))) > 0.) then
!            write (file_id,*) 'TROE/',troe_coeff(:,reac)
!          endif
!          if (minval(a_k4(:,reac)) < impossible) then
!            write (file_id,*) a_k4(:,reac)
!          endif
!        else
!          write (file_id,*) ' min lambda=',lamb_low,' max lambda=',lamb_up
!        endif
!
!      enddo
!
!      write (file_id,*) 'END'
!      write (file_id,*) '(M+) case: ',Mplus_case
!      write (file_id,*) 'photochemical case: ',photochem_case
!
!      close (file_id)
!
!    endsubroutine write_reactions
!***********************************************************************
!    subroutine chemkin_data
!
! All variables from this function apart from varname shared on the ogrid
!
!    endsubroutine chemkin_data
!***********************************************************************
!    subroutine chemkin_data_simple
!
! Does not involve f array, not computed separately for the ogrid
!
!    endsubroutine chemkin_data_simple
!***********************************************************************
!    subroutine read_reactions(input_file,NrOfReactions)
!
! Does not involve f array, not computed separately for the ogrid
!
!    endsubroutine read_reactions
!***********************************************************************
    subroutine get_reaction_rate_ogr(f_og,vreact_p,vreact_m,p_ogrid)
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
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f_og
      real, dimension(nx_ogrid,nreactions), intent(out) :: vreact_p, vreact_m
!
      type (pencil_case_ogrid) :: p_ogrid
      real, dimension(nx_ogrid) :: dSR=0., dHRT=0., Kp, Kc
      real, dimension(nx_ogrid) :: prod1, prod2
      real, dimension(nx_ogrid) :: kf=0., kr=0.
      real, dimension(nx_ogrid) :: rho_cgs, p_atm
      real, dimension(nx_ogrid) :: mix_conc
      integer :: k, reac, i
      real :: sum_tmp=0., ddd
      real :: Rcal, Rcal1, lnRgas, l10, lnp_atm
      logical, save :: lwrite_first=.true.
      character(len=fnlen) :: input_file="./data/react_ogr.out"
      integer :: file_id=123
      real :: B_n_0, alpha_n_0, E_an_0
      real, dimension(nx_ogrid) :: kf_0, Pr, sum_sp
      real, dimension(nx_ogrid) :: Fcent, ccc, nnn, lnPr, FF, tmpF
      real, dimension(nx_ogrid) :: TT1_loc
      real, dimension(5,nreactions) :: orders_m, orders_p
!
!  Check which reactions rate method we will use
!
      if (reac_rate_method == 'chemkin') then
!
        TT1_loc = p_ogrid%TT1
!
        if (lwrite_first)  open (file_id,file=input_file)
!
!  p is in atm units; atm/bar=1./10.13
!
        Rcal = Rgas_unit_sys/4.14*1e-7
        Rcal1 = 1./Rcal
        lnRgas = log(Rgas)
        l10 = log(10.)
        rho_cgs = p_ogrid%rho*unit_mass/unit_length**3
        lnp_atm = log(1e6*unit_length**3/unit_energy)
        p_atm = 1e6*(unit_length**3)/unit_energy
!
        if (lmech_simple) then
          orders_p(:,1) = (/ 0.0, 1.0, 0.0, 0.25, 0.5 /)
          orders_m(:,1) = (/ 1.0, 0.0, 0.0, 0.0, 0.0 /)
        endif
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
              prod1 = prod1*(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))**orders_p(k,reac)   
            endif
            if (abs(orders_m(k,reac)) > 0.0) then       
              prod2 = prod2*(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))**orders_m(k,reac)
            endif
          enddo
        else
          prod1 = 1.
          prod2 = 1.
          do k = 1,nchemspec
            if (abs(Sijp(k,reac)) == 1) then
              prod1 = prod1*(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))
            elseif (abs(Sijp(k,reac)) == 2) then
              prod1 = prod1*(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))*(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k)) &
                  *rho_cgs(:)/species_constants(k,imass))
            elseif (abs(Sijp(k,reac)) > 0) then
              prod1 = prod1*(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))**Sijp(k,reac)
            endif
          enddo
          do k = 1,nchemspec
            if (abs(Sijm(k,reac)) == 1.0) then
              prod2 = prod2*(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))
            elseif (abs(Sijm(k,reac)) == 2.0) then
              prod2 = prod2*(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))*(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k)) &
                  *rho_cgs(:)/species_constants(k,imass))
            elseif (abs(Sijm(k,reac)) > 0.0) then
              prod2 = prod2*(f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))*rho_cgs(:) &
                  /species_constants(k,imass))**Sijm(k,reac)
            endif
          enddo
        endif
!
!  Find forward rate constant for reaction 'reac'
!
          kf = B_n(reac)+alpha_n(reac)*p_ogrid%lnTT-E_an(reac)*Rcal1*TT1_loc
!
!  Find backward rate constant for reaction 'reac'
!
        if (lmech_simple) then
          !! 20.0301186564 = ln(5*10e8) !!
          kr = 20.0301186564-E_an(reac)*Rcal1*TT1_loc
          sum_sp = 1.
        else

          dSR = 0.
          dHRT = 0.
          sum_tmp = 0.
          do k = 1,nchemspec
            dSR = dSR+(Sijm(k,reac) -Sijp(k,reac))*p_ogrid%S0_R(:,k)
            dHRT = dHRT+(Sijm(k,reac)-Sijp(k,reac))*p_ogrid%H0_RT(:,k)
            sum_tmp = sum_tmp+(Sijm(k,reac)-Sijp(k,reac))
          enddo
          Kp = dSR-dHRT
!
          if (sum_tmp == 0.) then
            Kc = Kp
          else
            Kc = Kp+sum_tmp*(lnp_atm-p_ogrid%lnTT-lnRgas)
          endif
!
!  Multiply by third body reaction term
!
          if (minval(a_k4(:,reac)) < impossible) then
            sum_sp = 0.
            do k = 1,nchemspec
              sum_sp = sum_sp+a_k4(k,reac)*f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))  &
                  *rho_cgs(:)/species_constants(k,imass)
            enddo
            mix_conc = sum_sp
          else
            sum_sp = 1.
            mix_conc = rho_cgs(:)*(p_ogrid%RR(:)/Rgas)/unit_mass
            !mix_conc = rho_cgs(:)*p_ogrid%mu1(:)/unit_mass
          endif
!
!  The Lindeman approach to the fall of reactions
!
          if (maxval(abs(low_coeff(:,reac))) > 0.) then
            B_n_0 = low_coeff(1,reac)
            alpha_n_0 = low_coeff(2,reac)
            E_an_0 = low_coeff(3,reac)
            kf_0(:) = B_n_0+alpha_n_0*p_ogrid%lnTT(:)-E_an_0*Rcal1*TT1_loc(:)
            Pr = exp(kf_0-kf)*mix_conc
            kf = kf+log(Pr/(1.+Pr))
          elseif (maxval(abs(high_coeff(:,reac))) > 0.) then
            B_n_0 = high_coeff(1,reac)
            alpha_n_0 = high_coeff(2,reac)
            E_an_0 = high_coeff(3,reac)
            kf_0(:) = B_n_0+alpha_n_0*p_ogrid%lnTT(:)-E_an_0*Rcal1*TT1_loc(:)
            Pr = exp(kf_0-kf)*mix_conc
            kf = kf-log(1.+Pr)
          endif
!
! The Troe approach
!
          if (maxval(abs(troe_coeff(:,reac))) > 0.) then
            Fcent = (1.-troe_coeff(1,reac))*exp(-p_ogrid%TT(:)/troe_coeff(2,reac)) &
                +troe_coeff(1,reac)*exp(-p_ogrid%TT(:)/troe_coeff(3,reac))
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
!
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
    endsubroutine get_reaction_rate_ogr
!***********************************************************************
    subroutine calc_reaction_term_ogr(f_og,p_ogrid)
!
!  Calculation of the reaction term
!
      use Diagnostics, only: sum_mn_name, max_mn_name
!
      real :: alpha, eps
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) :: f_og
      real, dimension(nx_ogrid,nreactions) :: vreactions, vreactions_p, vreactions_m 
      real, dimension(nx_ogrid,nchemspec) :: xdot
      real, dimension(nx_ogrid) :: rho1
      real, dimension(nx_ogrid,nchemspec) :: molm
      type (pencil_case_ogrid) :: p_ogrid
      integer :: k,j,ii
!
      eps = sqrt(epsilon(alpha))
!
      p_ogrid%DYDt_reac = 0.
      rho1 = 1./p_ogrid%rho
      do k = 1,nchemspec
        molm(:,k) = rho1*species_constants(k,imass)
      enddo
!
!  Chemkin data case
!
      call get_reaction_rate_ogr(f_og,vreactions_p,vreactions_m,p_ogrid)
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
      p_ogrid%DYDt_reac = xdot*unit_time
!
!  For diagnostics
!
!      if (lchemistry_diag) then
!        do k = 1,nchemspec
!          do j = 1,nreactions
!            net_react_p(k,j) = net_react_p(k,j)+stoichio(k,j) &
!                *sum(vreactions_p(:,j))
!            net_react_m(k,j) = net_react_m(k,j)+stoichio(k,j) &
!                *sum(vreactions_m(:,j))
!          enddo
!        enddo
!      endif
!
!  Calculate diagnostic quantities
!
!      if (ldiagnos) then
!        do ii=1,nchemspec
!          if (idiag_dYm(ii) /= 0) then
!            call sum_mn_name(p_ogrid%DYDt_reac(:,ii),idiag_dYm(ii))
!          endif
!          if (idiag_dYmax(ii) /= 0) then
!            call max_mn_name(abs(p_ogrid%DYDt_reac(:,ii)),idiag_dYmax(ii))
!          endif
!          if (idiag_hm(ii) /= 0) then
!            call sum_mn_name(p_ogrid%H0_RT(:,ii)*Rgas* &
!                p_ogrid%TT(:)/species_constants(ii,imass),idiag_hm(ii))
!          endif
!        enddo
!      endif
!
    endsubroutine calc_reaction_term_ogr
!***********************************************************************
    subroutine calc_heter_reaction_term(f_og)
!
!  Calculation of the heterogeneous reaction term
!
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(inout) :: f_og
      real, dimension(2) :: vreact_p,kf,prod
      real, dimension(nchemspec) :: mdot = 0.
      real :: rho1,T_loc1,rho_cgs
      integer :: k,i
      real, dimension(2) :: B_n_het, E_an_het, alpha
      real :: Rcal1, Rcal, u_stefan
!
      rho_cgs = f_og(l1_ogrid,m_ogrid,n_ogrid,irho)
      rho1 = 1./rho_cgs
      T_loc1 = 1./f_og(l1_ogrid,m_ogrid,n_ogrid,iTT)

      Rcal  = Rgas_unit_sys/4.14*1e-7
      Rcal1 = 1./(Rcal)
!
      ! B & A checked in the paper that was cited - all ok
      B_n_het=(/1.97e9,1.291e7/)        ! cm/s
      E_an_het=(/47301.701,45635.267/)  ! cal/mol
!
!  Chemkin data case
!
      prod(1) = f_og(l1_ogrid,m_ogrid,n_ogrid,i_O2)*rho_cgs
      prod(2) = f_og(l1_ogrid,m_ogrid,n_ogrid,i_CO2)*rho_cgs
    !  if (prod(1) .lt. 0.0) prod(1)=0.0
    !  if (prod(2) .lt. 0.0) prod(2)=0.0
      do i=1,2
        kf(i)=log(B_n_het(i))-E_an_het(i)*Rcal1*T_loc1
        vreact_p(i)=prod(i)*exp(kf(i))
      enddo
!
      mdot(ichem_O2)=-vreact_p(1)
      mdot(ichem_CO)=2.0*(species_constants(ichem_CO,imass)/species_constants(ichem_O2,imass)*vreact_p(1)&
                        +species_constants(ichem_CO,imass)/species_constants(ichem_CO2,imass)*vreact_p(2))
      mdot(ichem_CO2)=-vreact_p(2)
!
      mdot=mdot*unit_time
!
      u_stefan = 0.
      do k = 1,nchemspec
        heter_reaction_rate(m_ogrid,n_ogrid,k)= mdot(k)
        u_stefan = u_stefan + mdot(k)
      enddo
      heter_reaction_rate(m_ogrid,n_ogrid,nchemspec+1) = -u_stefan
      ! the above = (2.*M_C/species_constants(ichem_O2,imass)*mdot(ichem_O2) &
      !             + M_C/species_constants(ichem_CO2,imass)*mdot(ichem_CO2))
      u_stefan = u_stefan/rho_cgs
   !   p_ogrid%DYDt_reac(1,:) = 0
      f_og(l1_ogrid,m_ogrid,n_ogrid,iux) = u_stefan
!
    endsubroutine calc_heter_reaction_term
!***********************************************************************
! TODO: do I need that for the ogrid?
!
!    subroutine  write_net_reaction
!
!  write net reactions to file
!
!      open (1,file=trim(datadir)//'/net_reactions.dat',position='append')
!      write (1,*) t
!      write (1,'(8e10.2)') net_react_p, net_react_m
!      close (1)
!
!    endsubroutine  write_net_reaction
!***********************************************************************
    subroutine calc_diffusion_term_ogrid(f_og,p_ogrid)
!
!  Calculate diffusion term, p%DYDt_diff
!
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
      type (pencil_case_ogrid) :: p_ogrid
!
      real, dimension(nx_ogrid) :: sum_gdiff_ogrid=0., gYk_sumdiff_ogrid
      real, dimension(nx_ogrid) :: del2Yk, diff_op1_ogrid, diff_op2_ogrid
      real, dimension(nx_ogrid,3) :: gYk, sum_diff_ogrid=0., dk_D_ogrid
      real, dimension(nx_ogrid,3,nchemspec) :: gDiff_full_add_ogrid
      integer :: k,i
!
      p_ogrid%DYDt_diff = 0.
!
!  Loop over all chemical species.
!
      do k = 1,nchemspec
!
          if (ldiffusion) then
!
!  Calculate diffusion coefficient gradients gDiff_full_add for the case:
!    2) Constant Lewis number (= 1 by default or read from start.in) and given Pr
!
              do i = 1,3
                gDiff_full_add_ogrid(:,i,k) = p_ogrid%gradnu(:,i)*Pr_number1*Lewis_coef1(k)
              enddo
!
!  Calculate the terms needed by the diffusion fluxes for the case:
!    1) Fickian diffusion law (gradient of species MASS fractions)
!
              call del2_ogrid(f_og,ichemspec(k),del2Yk)
              call grad_ogrid(f_og,ichemspec(k),gYk)
              call dot_mn_ogrid(p_ogrid%glnrho,gYk,diff_op1_ogrid)
              call dot_mn_ogrid(gDiff_full_add_ogrid(:,:,k),gYk,diff_op2_ogrid)
!
              p_ogrid%DYDt_diff(:,k) = p_ogrid%Diff_penc_add(:,k) &
                                     *(del2Yk+diff_op1_ogrid) + diff_op2_ogrid
              do i = 1,3
                dk_D_ogrid(:,i) = p_ogrid%Diff_penc_add(:,k)*gYk(:,i)
              enddo
!
              if (ldiff_corr) then
                do i = 1,3
                  sum_diff_ogrid(:,i) = sum_diff_ogrid(:,i)+dk_D_ogrid(:,i)
                enddo
                sum_gdiff_ogrid(:) = sum_gdiff_ogrid(:)+ p_ogrid%DYDt_diff(:,k)
              endif
          endif
      enddo
!
!  Adding correction diffusion velocity to ensure mass balance
!
      if (ldiffusion .and. ldiff_corr) then
        do k = 1,nchemspec
          call grad_ogrid(f_og,ichemspec(k),gYk)
          call dot_mn_ogrid(gYk(:,:),sum_diff_ogrid,gYk_sumdiff_ogrid)
          p_ogrid%DYDt_diff(:,k) = p_ogrid%DYDt_diff(:,k) &
            - (gYk_sumdiff_ogrid+f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ichemspec(k))*sum_gdiff_ogrid(:))
        enddo
      endif

    endsubroutine calc_diffusion_term_ogrid
!***********************************************************************
    subroutine calc_heatcond_chemistry_ogrid(f_og,df)
!
!  Calculate gamma*chi*(del2lnT+gradlnTT.grad(lnT+lnrho+lncp+lnchi))
!
      use Sub, only: dot
!
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) ::  f_og
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df
      real, dimension(nx_ogrid) :: g2TT, g2TTlambda=0., tmp1
!
!  Add heat conduction to RHS of temperature equation
!
        if (ltemperature_nolog) then
          call dot(p_ogrid%gTT,p_ogrid%glambda,g2TTlambda)
          tmp1 = (p_ogrid%lambda(:)*p_ogrid%del2TT+g2TTlambda)*p_ogrid%cv1*p_ogrid%rho1(:)
          df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT) = &
                 df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT) + tmp1
        else
          call dot(p_ogrid%glnTT,p_ogrid%glambda,g2TTlambda)
          call dot(p_ogrid%glnTT,p_ogrid%glnTT,g2TT)
          tmp1 = (p_ogrid%lambda(:)*(p_ogrid%del2lnTT+g2TT)&
               + g2TTlambda)*p_ogrid%cv1*p_ogrid%rho1(:)
          df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ilnTT) = &
                 df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ilnTT) + tmp1
        endif
!
!      call keep_compiler_quiet(f_og)
!
    endsubroutine calc_heatcond_chemistry_ogrid
!***********************************************************************
!    subroutine get_RHS_Y_full_ogrid(RHS_Y_ogrid)
!
!      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,nchemspec) :: RHS_Y_ogrid
!      intent(out) :: RHS_Y_ogrid
!
!      RHS_Y_ogrid = RHS_Y_full_ogrid
!
!    endsubroutine get_RHS_Y_full_ogrid
!***********************************************************************
!    subroutine get_cs2_full_ogr(cs2_full)
!
!    endsubroutine get_cs2_full_ogr
!***********************************************************************
!    subroutine get_cs2_slice(f,slice,dir,index)
!
! TODO: what do I do with BCs?
!
!    endsubroutine get_cs2_slice
!***********************************************************************
!    subroutine get_gamma_slice(f,slice,dir,index)
!
!    endsubroutine get_gamma_slice
!***********************************************************************
    subroutine air_field_ogr(f_og,PP)
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) ::  f_og
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid) :: sum_Y, tmp
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
!
      dataloop: do
!
        read(file_id,'(80A)',IOSTAT=iostat) ChemInpLine
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
          read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)') TT
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
          read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)') PP
          if (lroot) print*, ' Pressure, Pa   ', PP
!
        elseif (tmp_string == 'V') then
!
          StartInd=1; StopInd =0
!
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
          StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
          StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
!
          read (unit=ChemInpLine(StartInd:StopInd),fmt='(E14.7)') velx
          if (lroot) print*, ' Velocity, cm/s   ', velx
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
            read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)') YY_k
            if (lroot) print*, ' volume fraction, %,    ', YY_k, &
               species_constants(ind_chem,imass)
!
            if (species_constants(ind_chem,imass)>0.) then
             air_mass=air_mass+YY_k*0.01/species_constants(ind_chem,imass)
            endif
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
! Stop if air.dat is empty
!
      if (emptyFile)  call stop_it('The input file air.dat was empty!')
       air_mass=1./air_mass
!
      do j=1,k-1
        f_og(:,:,:,ichemspec(stor1(j)))=stor2(j)*0.01
      enddo

      sum_Y=0.

      do j=1,nchemspec
        sum_Y=sum_Y+f_og(:,:,:,ichemspec(j))
      enddo
      do j=1,nchemspec
        f_og(:,:,:,ichemspec(j))=f_og(:,:,:,ichemspec(j))/sum_Y
        initial_massfractions(j)=f_og(l1_ogrid,m1_ogrid,n1_ogrid,ichemspec(j))
      enddo

      if (mvar < 5) then
        call fatal_error("air_field_ogr", "I can only set existing fields")
      endif

      if (.not. reinitialize_chemistry) then
      if (.not.lflame_front)  then
        if (ltemperature_nolog) then
          f_og(:,:,:,iTT)=TT
        else
          f_og(:,:,:,ilnTT)=alog(TT)!+f_og(:,:,:,ilnTT)
        endif
        if (ldensity_nolog) then
          f_og(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
        else
          tmp=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
            f_og(:,:,:,ilnrho)=alog(tmp)
        endif
        if (nxgrid_ogrid>1) f_og(:,:,:,iux)=f_og(:,:,:,iux)+init_ux
      endif
      endif

      if (init_from_file) then
        if (lroot) print*, 'Velocity field read from file, initialization' //&
            'of density and temperature'
        if (ldensity_nolog) then
          f_og(:,:,:,ilnrho)=f_og(:,:,:,ilnrho)*(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
        else
          tmp=f_og(:,:,:,ilnrho)*(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/TT)/unit_mass*unit_length**3
          f_og(:,:,:,ilnrho)=alog(tmp)
        endif
        if (ltemperature_nolog) then
          if (ldensity_nolog) then
            f_og(:,:,:,iTT)=(PP/(k_B_cgs/m_u_cgs)*&
                air_mass/f_og(:,:,:,ilnrho))/unit_mass*unit_length**3
          else
            f_og(:,:,:,iTT)=(PP/(k_B_cgs/m_u_cgs)*&
                air_mass/exp(f_og(:,:,:,ilnrho)))/unit_mass*unit_length**3
          endif
        else
          if (ldensity_nolog) then
            tmp=(PP/(k_B_cgs/m_u_cgs)*&
                air_mass/f_og(:,:,:,ilnrho))/unit_mass*unit_length**3
            f_og(:,:,:,ilnTT)=alog(tmp)!+f_og(:,:,:,ilnTT)
          else
            tmp=(PP/(k_B_cgs/m_u_cgs)*&
                air_mass/exp(f_og(:,:,:,ilnrho)))/unit_mass*unit_length**3
            f_og(:,:,:,ilnTT)=alog(tmp)!+f_og(:,:,:,ilnTT)
          endif
        endif
        if (velx/=0.) f_og(:,:,:,iux)=f_og(:,:,:,iux)+velx
      endif

      if (linit_temperature) then
        do i=1,mx_ogrid
        if (x(i)<=init_x1) then
          f_og(i,:,:,ilnTT)=alog(init_TT1)
        endif
        if (x(i)>=init_x2) then
          f_og(i,:,:,ilnTT)=alog(init_TT2)
        endif
        if (x(i)>init_x1 .and. x(i)<init_x2) then
          if (init_x1 /= init_x2) then
            f_og(i,:,:,ilnTT)=&
               alog((x(i)-init_x1)/(init_x2-init_x1) &
               *(init_TT2-init_TT1)+init_TT1)
          endif
        endif
        enddo
      endif

      if (linit_density) then
        if (ldensity_nolog) then
          f_og(:,:,:,ilnrho)=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/exp(f_og(:,:,:,ilnTT)))/unit_mass*unit_length**3
        else
          tmp=(PP/(k_B_cgs/m_u_cgs)*&
            air_mass/exp(f_og(:,:,:,ilnTT)))/unit_mass*unit_length**3
          f_og(:,:,:,ilnrho)=alog(tmp)
        endif
      endif

      if (nxgrid>1) then
        f_og(:,:,:,iux)=f_og(:,:,:,iux)+velx
      endif


      if (lroot) print*, 'Air temperature, K', TT
      if (lroot) print*, 'Air pressure, dyn', PP
      if (lroot) print*, 'Air density, g/cm^3:'
      if (lroot) print '(E10.3)',  PP/(k_B_cgs/m_u_cgs)*air_mass/TT
      if (lroot) print*, 'Air mean weight, g/mol', air_mass
      if (lroot) print*, 'R', k_B_cgs/m_u_cgs

      close(file_id)
!
    endsubroutine air_field_ogr
!***********************************************************************
!          NSCBC boundary conditions
!***********************************************************************
!    subroutine damp_zone_for_NSCBC(f,df) was here
!***********************************************************************
! Here was the subroutine jacobn(f,jacob) which is called from timestep_stiff.f90
!***********************************************************************
!    subroutine get_mu1_slice(f,slice,grad_slice,index,sgn,direction)
!
!    endsubroutine get_mu1_slice
!***********************************************************************
!    subroutine get_reac_rate(wt,f_og,p_ogrid)
!
!      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) ::  f_og
!      real, dimension(nx_ogrid) :: wt
!      real, dimension(nx_ogrid,nchemspec) :: ydot
!      type (pencil_case_ogrid) :: p_ogrid
!
!      ydot = p_ogrid%DYDt_reac
!      f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ireaci(1):ireaci(nchemspec)) = ydot
!         f(l1:l2,m,n,ireaci(1):ireaci(nchemspec))+ydot
!
!      if (maux == nchemspec+1) f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ireaci(nchemspec)+1) = wt
!         f(l1:l2,m,n,ireaci(nchemspec)+1)+wt
!
!    endsubroutine get_reac_rate
!***********************************************************************
!    subroutine chemistry_clean_up()
!
!    endsubroutine chemistry_clean_up
!***********************************************************************
!    subroutine find_remove_real_stoic(Speciesstring,lreal,stoi,startindex)
!
! Does not involve f array, not executed separately for the ogrid
!
!    endsubroutine find_remove_real_stoic
!!***********************************************************************
! subroutine chemspec_normalization_N2(f) was here
!***********************************************************************
    subroutine getmu_array_ogrid(f_og,mu1_full_og)
!
!  Calculate mean molecular weight
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) ::  f_og
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid) :: mu1_full_og
      integer :: k,j2,j3
!
!  Mean molecular weight
!
      mu1_full_og=0.
      do k=1,nchemspec
        if (species_constants(k,imass)>0.) then
          do j2=mm1_ogrid,mm2_ogrid
            do j3=nn1_ogrid,nn2_ogrid
              mu1_full_og(:,j2,j3)= &
                  mu1_full_og(:,j2,j3)+f_og(:,j2,j3,ichemspec(k))/species_constants(k,imass)
            enddo
          enddo
        endif
      enddo
!
    endsubroutine getmu_array_ogrid
!***********************************************************************
!   subroutine read_Lewis
!
! Removed. Lewis # read from start.in
!
 !   endsubroutine read_Lewis
!***********************************************************************
!   subroutine read_transport_data
!
! Does not involve f array, not executed separately for the ogrid
!
!    endsubroutine read_transport_data
!***********************************************************************
!    subroutine jacobn_ogr(f_og,jacob)
!
!    endsubroutine jacobn_ogr
!***********************************************************************
    subroutine chemspec_normalization_og(f_og)
!
!   20-sep-10/Natalia: coded
!   renormalization of the species
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) ::  f_og
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid) :: sum_Y
      integer :: k
!
      sum_Y = 0.
      do k = 1,nchemspec
        sum_Y = sum_Y+f_og(:,:,:,ichemspec(k))
      enddo
      do k = 1,nchemspec
        f_og(:,:,:,ichemspec(k)) = f_og(:,:,:,ichemspec(k))/sum_Y
      enddo
!
    endsubroutine chemspec_normalization_og
!***********************************************************************
    subroutine chemspec_normalization_N2_og(f_og)
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) ::  f_og
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid) :: sum_Y
      integer :: k ,isN2, ichemsN2
      logical :: lsN2
!
      call find_species_index('N2', isN2, ichemsN2, lsN2 )
      sum_Y=0.0
      do k=1,nchemspec
        if (k/=ichemsN2) then
          sum_Y=sum_Y+f_og(:,:,:,ichemspec(k))
        endif
      enddo
      f_og(:,:,:,isN2)=1.0-sum_Y
!
    endsubroutine chemspec_normalization_N2_og
!***********************************************************************
endmodule solid_cells_ogrid_chemistry
