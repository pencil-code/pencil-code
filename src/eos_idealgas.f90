! $Id$
!
!  Equation of state for an ideal gas without ionization.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cp; cp1; cp1tilde
! PENCILS PROVIDED glnTT(3); TT; TT1; gTT(3); yH; hss(3,3); hlnTT(3,3)
! PENCILS PROVIDED del2ss; del6ss; del2lnTT; cv1; del6lnTT; gamma
! PENCILS PROVIDED del2TT; del6TT; glnmumol(3); ppvap; csvap2
! PENCILS PROVIDED TTb; rho_anel; eth; geth(3); del2eth
!
!***************************************************************
module EquationOfState
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'eos.h'
!
  interface eoscalc ! Overload subroutine `eoscalc' function
    module procedure eoscalc_pencil   ! explicit f implicit m,n
    module procedure eoscalc_point    ! explicit lnrho, ss
    module procedure eoscalc_farray   ! explicit lnrho, ss
  end interface
!
  interface pressure_gradient ! Overload subroutine `pressure_gradient'
    module procedure pressure_gradient_farray  ! explicit f implicit m,n
    module procedure pressure_gradient_point   ! explicit lnrho, ss
  end interface
!
  integer, parameter :: ilnrho_ss=1, ilnrho_ee=2, ilnrho_pp=3
  integer, parameter :: ilnrho_lnTT=4, ilnrho_cs2=5
  integer, parameter :: irho_cs2=6, irho_ss=7, irho_lnTT=8, ilnrho_TT=9
  integer, parameter :: irho_TT=10, ipp_ss=11, ipp_cs2=12, irho_eth=13
  integer, parameter :: ilnrho_eth=14, irho_rhop=15
  integer :: iglobal_cs2, iglobal_glnTT
  real, dimension(mz) :: profz_eos=1.0,dprofz_eos=0.0
  real, dimension(3) :: beta_glnrho_global=0.0, beta_glnrho_scaled=0.0
  real :: lnTT0=impossible, TT0=impossible
  real :: xHe=0.0
  real :: mu=1.0
  real :: cs0=1.0, rho0=1.0, rho01=1.0, pp0=1.0
  real :: cs20=1.0, lnrho0=0.0
  real :: gamma=5.0/3.0
  real :: Rgas_cgs=0.0, Rgas, error_cp=1.0e-6
  real :: gamma_m1    !(=gamma-1)
  real :: gamma1   !(=1/gamma)
  real :: cp=impossible, cp1=impossible, cv=impossible, cv1=impossible
  real :: pres_corr=0.1
  real :: cs2top_ini=impossible, dcs2top_ini=impossible
  real :: cs2bot=1.0, cs2top=1.0
  real :: cs2cool=0.0
  real :: mpoly=1.5, mpoly0=1.5, mpoly1=1.5, mpoly2=1.5
  real :: width_eos_prof=0.2
  real :: sigmaSBt=1.0
  real :: TT_floor=impossible
  integer :: isothtop=0
  integer :: ieosvars=-1, ieosvar1=-1, ieosvar2=-1, ieosvar_count=0
  logical :: leos_isothermal=.false., leos_isentropic=.false.
  logical :: leos_isochoric=.false., leos_isobaric=.false.
  logical :: leos_localisothermal=.false.
  logical :: lanelastic_lin=.false.
  character (len=labellen) :: ieos_profile='nothing'
!
  real, dimension(nchemspec,18) :: species_constants
  real, dimension(nchemspec,7)     :: tran_data
  real, dimension(nchemspec)  :: Lewis_coef, Lewis_coef1
!
!  Input parameters.
!
  namelist /eos_init_pars/ &
      xHe, mu, cp, cs0, rho0, gamma, error_cp, cs2top_ini, &
      dcs2top_ini, sigmaSBt, lanelastic_lin
!
!  Run parameters.
!
  namelist /eos_run_pars/ &
      xHe, mu, cp, cs0, rho0, gamma, error_cp, cs2top_ini,           &
      dcs2top_ini, ieos_profile, width_eos_prof,pres_corr, sigmaSBt, &
      lanelastic_lin, TT_floor
!
  contains
!***********************************************************************
    subroutine register_eos()
!
!  Register variables from the EquationOfState module.
!
!  14-jun-03/axel: adapted from register_eos
!
      leos_idealgas=.true.
!
      iyH=0
      ilnTT=0
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_eos: ionization nvar = ', nvar
      endif
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          '$Id$')
!
    endsubroutine register_eos
!***********************************************************************
    subroutine units_eos()
!
!  This routine calculates things related to units and must be called
!  before the rest of the units are being calculated.
!
!  22-jun-06/axel: adapted from initialize_eos
!   4-aug-09/axel: added possibility of vertical profile function
!
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable
      use Sub, only: erfunc
!
      real :: Rgas_unit_sys, cp_reference
      integer :: ierr
!
!  Set gamma_m1, cs20, and lnrho0, and rho01.
!  (used currently for non-dimensional equation of state)
!
      gamma_m1=gamma-1.0
      gamma1=1/gamma
      lnrho0=log(rho0)
      rho01 = 1./rho0
!
!  Avoid floating overflow if cs0 was not set.
!
      cs20=cs0**2
!
!  Initialize variable selection code (needed for RELOADing).
!
      ieosvars=-1
      ieosvar_count=0
!
!  Unless unit_temperature is set, calculate by default with cp=1.
!  If unit_temperature is set, cp must follow from this.
!  Conversely, if cp is set, then unit_temperature must follow from this.
!  If unit_temperature and cp are set, the problem is overdetermined,
!  but it may still be correct, so this will be checked here.
!  When gamma=1.0 (gamma_m1=0.0), write Rgas=mu*cp or cp=Rgas/mu.
!
      if (unit_system=='cgs') then
        Rgas_unit_sys=k_B_cgs/m_u_cgs
      elseif (unit_system=='SI') then
        Rgas_unit_sys=k_B_cgs/m_u_cgs*1.0e-4
      endif
!
      if (unit_temperature==impossible) then
        if (cp==impossible) cp=1.0
        if (gamma_m1==0.0) then
          Rgas=mu*cp
        else
          Rgas=mu*(1.0-gamma1)*cp
        endif
        unit_temperature=unit_velocity**2*Rgas/Rgas_unit_sys
      else
        Rgas=Rgas_unit_sys*unit_temperature/unit_velocity**2
        if (cp==impossible) then
          if (gamma_m1==0.0) then
            cp=Rgas/mu
          else
            cp=Rgas/(mu*gamma_m1*gamma1)
          endif
        else
!
!  Checking whether the units are overdetermined.
!  This is assumed to be the case when the two differ by error_cp.
!
          if (gamma_m1==0.0) then
            cp_reference=Rgas/mu
          else
            cp_reference=Rgas/(mu*gamma_m1*gamma1)
          endif
          if (abs(cp-cp_reference)/cp > error_cp) then
            if (lroot) print*,'initialize_eos: consistency: cp=', cp , &
                'while: cp_reference=', cp_reference
            call fatal_error('units_eos','cp is not correctly calculated')
          endif
        endif
      endif
      cp1=1/cp
      cv=gamma1*cp
      cv1=gamma*cp1
!
!  Need to calculate the equivalent of cs0.
!  Distinguish between gamma=1 case and not.
!
      if (gamma_m1/=0.0) then
        lnTT0=log(cs20/(cp*gamma_m1))  !(general case)
      else
        lnTT0=log(cs20/cp)  !(isothermal/polytropic cases: check!)
      endif
      pp0=Rgas*exp(lnTT0)*rho0
      TT0=exp(lnTT0)
!
! Shared variables
!
      call put_shared_variable('cs20',cs20,ierr)
      if (ierr/=0) call fatal_error('units_eos','problem when putting cs20')
!
      call put_shared_variable('mpoly',mpoly,ierr)
      if (ierr/=0) call fatal_error('units_eos','problem when putting mpoly')
!
      call put_shared_variable('gamma',gamma,ierr)
      if (ierr/=0) call fatal_error('units_eos','problem when putting gamma')
!
!  Check that everything is OK.
!
      if (lroot) then
        print*, 'initialize_eos: unit_temperature=', unit_temperature
        print*, 'initialize_eos: cp, lnTT0, cs0, pp0=', cp, lnTT0, cs0, pp0
      endif
!
!  Calculate profile functions (used as prefactors to turn off pressure
!  gradient term).
!
      if (ieos_profile=='nothing') then
        profz_eos=1.0
        dprofz_eos=0.0
      elseif (ieos_profile=='surface_z') then
        profz_eos=0.5*(1.0-erfunc(z/width_eos_prof))
        dprofz_eos=-exp(-(z/width_eos_prof)**2)/(sqrtpi*width_eos_prof)
      endif
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos()
!
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable
!
      integer :: ierr
!
!  Perform any post-parameter-read initialization
!
!  Initialize variable selection code (needed for RELOADing).
!
      ieosvars=-1
      ieosvar_count=0
!
!  Write constants to disk. In future we may want to deal with this
!  using an include file or another module.
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,'(a,1pd26.16)') 'k_B=',k_B
        write (1,'(a,1pd26.16)') 'm_H=',m_H
        write (1,*) 'lnrho0=',lnrho0
        write (1,*) 'lnTTO=',lnTT0
        write (1,*) 'cs20=',cs20
        write (1,*) 'cp=',cp
        close (1)
      endif
!
      call put_shared_variable('cp',cp,ierr)
        if (ierr/=0) call stop_it("cp: "//&
             "there was a problem when sharing cp")
      call put_shared_variable('cv',cv,ierr)
        if (ierr/=0) call stop_it("cv: "//&
             "there was a problem when sharing cv")
!
     if (lanelastic) then
        call put_shared_variable('lanelastic_lin',lanelastic_lin,ierr)
        if (ierr/=0) call stop_it("lanelastic_lin: "//&
             "there was a problem when sharing lanelastic_lin")
      endif
!
    endsubroutine initialize_eos
!***********************************************************************
    subroutine select_eos_variable(variable,findex)
!
!  Select eos variable.
!
!   02-apr-06/tony: implemented
!
      use FArrayManager, only: farray_register_global
!
      character (len=*), intent(in) :: variable
      integer, intent(in) :: findex
      integer :: this_var=-1
      integer, save :: ieosvar_selected=0
      integer, parameter :: ieosvar_lnrho = 2**0
      integer, parameter :: ieosvar_rho   = 2**1
      integer, parameter :: ieosvar_ss    = 2**2
      integer, parameter :: ieosvar_lnTT  = 2**3
      integer, parameter :: ieosvar_TT    = 2**4
      integer, parameter :: ieosvar_cs2   = 2**5
      integer, parameter :: ieosvar_pp    = 2**6
      integer, parameter :: ieosvar_eth   = 2**7
      integer, parameter :: ieosvar_rhop  = 2**8
!
      if (ieosvar_count==0) ieosvar_selected=0
!
      if (ieosvar_count>=2) &
          call fatal_error('select_eos_variable', &
          '2 thermodynamic quantities have already been defined '// &
          'while attempting to add a 3rd')
!
      ieosvar_count=ieosvar_count+1
!
      if (variable=='ss') then
        this_var=ieosvar_ss
        if (findex<0) then
          leos_isentropic=.true.
        endif
      elseif (variable=='cs2') then
        this_var=ieosvar_cs2
        if (findex==-2) then
          leos_localisothermal=.true.
          call farray_register_global('cs2',iglobal_cs2)
          call farray_register_global('glnTT',iglobal_glnTT,vector=3)
        elseif (findex<0) then
          leos_isothermal=.true.
        endif
      elseif (variable=='lnTT') then
        this_var=ieosvar_lnTT
        if (findex<0) then
          leos_isothermal=.true.
        endif
      elseif (variable=='TT') then
        this_var=ieosvar_TT
      elseif (variable=='lnrho') then
        this_var=ieosvar_lnrho
        if (findex<0) then
          leos_isochoric=.true.
        endif
      elseif (variable=='rho') then
        this_var=ieosvar_rho
        if (findex<0) then
          leos_isochoric=.true.
        endif
      elseif (variable=='pp') then
        this_var=ieosvar_pp
        if (findex<0) then
          leos_isobaric=.true.
        endif
      elseif (variable=='eth') then
        this_var=ieosvar_eth
      elseif (variable=='rhop') then
        this_var=ieosvar_rhop
      else
        call fatal_error('select_eos_variable','unknown thermodynamic variable')
      endif
      if (ieosvar_count==1) then
        ieosvar1=findex
        ieosvar_selected=ieosvar_selected+this_var
        return
      endif
!
!  Ensure the indexes are in the correct order.
!
      if (this_var<ieosvar_selected) then
        ieosvar2=ieosvar1
        ieosvar1=findex
      else
        ieosvar2=findex
      endif
      ieosvar_selected=ieosvar_selected+this_var
      select case (ieosvar_selected)
        case (ieosvar_lnrho+ieosvar_ss)
          if (lroot) print*, 'select_eos_variable: Using lnrho and ss'
          ieosvars=ilnrho_ss
        case (ieosvar_rho+ieosvar_ss)
          if (lroot) print*, 'select_eos_variable: Using rho and ss'
          ieosvars=irho_ss
        case (ieosvar_lnrho+ieosvar_lnTT)
          if (lroot) print*, 'select_eos_variable: Using lnrho and lnTT'
          ieosvars=ilnrho_lnTT
        case (ieosvar_lnrho+ieosvar_TT)
          if (lroot) print*, 'select_eos_variable: Using lnrho and TT'
          ieosvars=ilnrho_TT
        case (ieosvar_rho+ieosvar_lnTT)
          if (lroot) print*, 'select_eos_variable: Using rho and lnTT'
          ieosvars=irho_lnTT
        case (ieosvar_lnrho+ieosvar_cs2)
          if (lroot) print*, 'select_eos_variable: Using lnrho and cs2'
          ieosvars=ilnrho_cs2
        case (ieosvar_rho+ieosvar_cs2)
          if (lroot) print*, 'select_eos_variable: Using rho and cs2'
          ieosvars=irho_cs2
        case (ieosvar_rho+ieosvar_TT)
          if (lroot) print*, 'select_eos_variable: Using rho and TT'
          ieosvars=irho_TT
        case (ieosvar_pp+ieosvar_ss)
          if (lroot) print*, 'select_eos_variable: Using pp and ss'
          ieosvars=ipp_ss
        case (ieosvar_pp+ieosvar_cs2)
          if (lroot) print*, 'select_eos_variable: Using pp and cs2'
          ieosvars=ipp_cs2
        case (ieosvar_rho+ieosvar_eth)
          if (lroot) print*, 'select_eos_variable: Using rho and eth'
          ieosvars=irho_eth
        case (ieosvar_lnrho+ieosvar_eth)
          if (lroot) print*, 'select_eos_variable: Using lnrho and eth'
          ieosvars=ilnrho_eth
        case (ieosvar_rho+ieosvar_rhop)
          if (lroot) print*, 'select_eos_variable: Using rho and rhop'
          ieosvars=irho_rhop
        case default
          if (lroot) print*, 'select_eos_variable: '// &
              'Thermodynamic variable combination, ieosvar_selected =', &
              ieosvar_selected
          call fatal_error('select_eos_variable', &
              'This thermodynamic variable combination is not implemented')
      endselect
!
    endsubroutine select_eos_variable
!***********************************************************************
    subroutine getmu(f,mu_tmp)
!
!  Calculate average particle mass in the gas relative to
!
!   12-aug-03/tony: implemented
!
      real, dimension (mx,my,mz,mfarray), optional :: f
      real, intent(out) :: mu_tmp
!
!  mu = mu_H * (1 - xHe) + mu_He * xHe
!     = mu_H + (mu_He-mu_H) * xHe
!  mu_H = 1.
!  mu_He = 4.0026 / 1.0079  (molar masses from a Periodic Table)
!        = 3.97
!
      if (mu==0.0) then
        mu_tmp=1.0+2.97153*xHe
      else
        mu_tmp=mu
      endif
!
      call keep_compiler_quiet(present(f))
!
    endsubroutine getmu
!***********************************************************************
    subroutine getmu_array(f,mu1_full_tmp)
!
!  dummy routine to calculate mean molecular weight
!
!   16-mar-10/natalia
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: mu1_full_tmp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(mu1_full_tmp)
!
    endsubroutine getmu_array
!***********************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
!  Writes iyH and ilnTT to index.pro file.
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_eos
!***********************************************************************
    subroutine get_slices_eos(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_eos
!***********************************************************************
    subroutine pencil_criteria_eos()
!
!  All pencils that the EquationOfState module depends on are specified here.
!
!  02-04-06/tony: coded
!
    endsubroutine pencil_criteria_eos
!***********************************************************************
    subroutine pencil_interdep_eos(lpencil_in)
!
!  Interdependency among pencils from the EquationOfState module is specified
!  here.
!
!  20-nov-04/anders: coded
!  15-jul-10/axel: added gTT calculation for ilnrho_ss,irho_ss case
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_cp)) lpencil_in(i_cp1)=.true.
!
      select case (ieosvars)
!
!  Pencils for thermodynamic quantities for given lnrho or rho and ss.
!
      case (ilnrho_ss,irho_ss)
        if (leos_isentropic) then
          if (lpencil_in(i_cs2)) lpencil_in(i_lnrho)=.true.
        elseif (leos_isothermal) then
          if (lpencil_in(i_ss)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_gss)) lpencil_in(i_glnrho)=.true.
          if (lpencil_in(i_hss)) lpencil_in(i_hlnrho)=.true.
          if (lpencil_in(i_del2ss)) lpencil_in(i_del2lnrho)=.true.
          if (lpencil_in(i_del6ss)) lpencil_in(i_del6lnrho)=.true.
        elseif (leos_localisothermal) then
        else
          if (lpencil_in(i_cs2)) then
            lpencil_in(i_ss)=.true.
            lpencil_in(i_lnrho)=.true.
          endif
        endif
        if (lpencil_in(i_lnTT)) then
          lpencil_in(i_ss)=.true.
          lpencil_in(i_lnrho)=.true.
        endif
        if (lpencil_in(i_pp)) then
          lpencil_in(i_lnTT)=.true.
          lpencil_in(i_lnrho)=.true.
        endif
        if (lpencil_in(i_ee)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_glnTT)) then
          lpencil_in(i_glnrho)=.true.
          lpencil_in(i_gss)=.true.
        endif
        if (lpencil_in(i_gTT)) then
          lpencil_in(i_glnTT)=.true.
          lpencil_in(i_TT)=.true.
        endif
        if (lpencil_in(i_del2lnTT)) then
          lpencil_in(i_del2lnrho)=.true.
          lpencil_in(i_del2ss)=.true.
        endif
        if (lpencil_in(i_hlnTT)) then
          lpencil_in(i_hlnrho)=.true.
          lpencil_in(i_hss)=.true.
        endif
!
!  Pencils for thermodynamic quantities for given lnrho or rho and lnTT.
!
      case (ilnrho_lnTT,irho_lnTT)
        if (leos_isentropic) then
          if (lpencil_in(i_lnTT)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_glnTT)) lpencil_in(i_glnrho)=.true.
          if (lpencil_in(i_hlnTT)) lpencil_in(i_hlnrho)=.true.
          if (lpencil_in(i_del2lnTT)) lpencil_in(i_del2lnrho)=.true.
          if (lpencil_in(i_cs2)) lpencil_in(i_lnrho)=.true.
        elseif (leos_isothermal) then
        elseif (leos_localisothermal) then
        else
          if (lpencil_in(i_cs2)) lpencil_in(i_lnTT)=.true.
        endif
        if (lpencil_in(i_ss)) then
          lpencil_in(i_lnTT)=.true.
          lpencil_in(i_lnrho)=.true.
        endif
        if (lpencil_in(i_pp)) then
          lpencil_in(i_lnTT)=.true.
          lpencil_in(i_lnrho)=.true.
        endif
        if (lpencil_in(i_ee)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_gss)) then
          lpencil_in(i_glnTT)=.true.
          lpencil_in(i_glnrho)=.true.
        endif
        if (lpencil_in(i_del2ss)) then
          lpencil_in(i_del2lnTT)=.true.
          lpencil_in(i_del2lnrho)=.true.
        endif
        if (lpencil_in(i_hss)) then
          lpencil_in(i_hlnTT)=.true.
          lpencil_in(i_hlnrho)=.true.
        endif
        if (lpencil_in(i_gTT)) then
          lpencil_in(i_glnTT)=.true.
        endif
!
!  Pencils for thermodynamic quantities for given lnrho or rho and TT.
!
      case (ilnrho_TT,irho_TT)
        if (lpencil_in(i_glnTT)) then
          lpencil_in(i_gTT)=.true.
          lpencil_in(i_TT1)=.true.
        endif
        if (lpencil_in(i_ss)) then
          lpencil_in(i_lnTT)=.true.
          lpencil_in(i_lnrho)=.true.
        endif
        if (lpencil_in(i_del2lnTT)) then
          lpencil_in(i_glnTT)=.true.
          lpencil_in(i_TT1)=.true.
        endif
!
!  Pencils for thermodynamic quantities for given lnrho or rho and cs2.
!
      case (ilnrho_cs2,irho_cs2)
        if (leos_isentropic) then
           call fatal_error('eos_isentropic', 'isentropic case not yet coded')
        elseif (leos_isothermal) then
          if (lpencil_in(i_ss)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_del2ss)) lpencil_in(i_del2lnrho)=.true.
          if (lpencil_in(i_gss)) lpencil_in(i_glnrho)=.true.
          if (lpencil_in(i_hss)) lpencil_in(i_hlnrho)=.true.
          if (lpencil_in(i_pp)) lpencil_in(i_rho)=.true.
        endif
!
!  Pencils for thermodynamic quantities for given pp and ss (anelastic case only).
!
      case(ipp_ss)
        if (leos_isentropic) then
           call fatal_error('eos_isentropic', 'isentropic case not yet coded')
        elseif (leos_isothermal) then
          if (lpencil_in(i_lnrho)) then
            lpencil_in(i_pp)=.true.
!            lpencil_in(i_TT)=.true.
          endif
          if (lpencil_in(i_rho)) lpencil_in(i_lnrho)=.true.
        else
          lpencil_in(i_rho)=.true.
          lpencil_in(i_pp)=.true.
          lpencil_in(i_ss)=.true.
          if (lpencil_in(i_lnrho)) lpencil_in(i_rho)=.true.
          if (lpencil_in(i_lnTT)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_lnTT)) lpencil_in(i_ss)=.true.
          if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
          if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
!         if (lpencil_in(i_lnrho)) then
!           lpencil_in(i_pp)=.true.
!           lpencil_in(i_ss)=.true.
!         endif
          if (lpencil_in(i_rho_anel)) then
              lpencil_in(i_pp)=.true.
              lpencil_in(i_ss)=.true.
          endif
        endif
!
      case (ipp_cs2)
        if (leos_isentropic) then
           call fatal_error('eos_isentropic', 'isentropic case not yet coded')
        elseif (leos_isothermal) then
          if (lpencil_in(i_lnrho)) then
            lpencil_in(i_pp)=.true.
          endif
          if (lpencil_in(i_rho)) lpencil_in(i_lnrho)=.true.
        else
          if (lpencil_in(i_rho)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_TT1)) lpencil_in(i_TT)=.true.
          if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
        endif
!
      case (irho_eth,ilnrho_eth)
        if (lpencil_in(i_cs2).or. &
            lpencil_in(i_TT).or. &
            lpencil_in(i_lnTT).or. &
            lpencil_in(i_TT1)) then
          lpencil_in(i_eth)=.true.
          lpencil_in(i_rho1)=.true.
        endif
        if (lpencil_in(i_pp)) lpencil_in(i_eth)=.true.
        if (lpencil_in(i_ee)) then
          lpencil_in(i_rho1)=.true.
          lpencil_in(i_eth)=.true.
        endif
        if (lpencil_in(i_TT)) then
          lpencil_in(i_cv1)=.true.
          lpencil_in(i_rho1)=.true.
          lpencil_in(i_eth)=.true.
        endif
        if (lpencil_in(i_lnTT)) lpencil_in(i_TT)=.true.
        if (lpencil_in(i_TT1)) lpencil_in(i_TT)=.true.
        if (lpencil_in(i_gTT)) then
          lpencil_in(i_rho1)=.true.
          lpencil_in(i_cv1)=.true.
          lpencil_in(i_geth)=.true.
          lpencil_in(i_TT)=.true.
          lpencil_in(i_rho)=.true.
        endif
        if (lpencil_in(i_del2TT)) then
          lpencil_in(i_rho1)=.true.
          lpencil_in(i_cv1)=.true.
          lpencil_in(i_del2eth)=.true.
          lpencil_in(i_TT)=.true.
          lpencil_in(i_del2rho)=.true.
          lpencil_in(i_grho)=.true.
          lpencil_in(i_gTT)=.true.
        endif
      case (irho_rhop) 
        lpencil_in(i_rho)=.true.
        lpencil_in(i_rhop)=.true.
        if (lpencil_in(i_lnTT)) lpencil_in(i_TT)=.true.
        if (lpencil_in(i_pp)) lpencil_in(i_TT)=.true.
        if (lpencil_in(i_cs2)) then 
          lpencil_in(i_cp)=.true.
          lpencil_in(i_TT)=.true.
        endif
        if (lpencil_in(i_ss)) then 
          lpencil_in(i_cp)=.true.
          lpencil_in(i_TT)=.true.
        endif
      case default
        call fatal_error('pencil_interdep_eos','case not implemented yet')
      endselect
!
    endsubroutine pencil_interdep_eos
!***********************************************************************
    subroutine calc_pencils_eos(f,p)
!
!  Calculate EquationOfState pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  02-apr-06/tony: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      real, dimension(nx) :: tmp
      integer :: i,j
!
!  Inverse cv and cp values.
!
      if (lpencil(i_cv1)) p%cv1=cv1
      if (lpencil(i_cp1)) p%cp1=cp1
      if (lpencil(i_cp))  p%cp=1/p%cp1
      if (lpencil(i_cp1tilde)) p%cp1tilde=cp1
!
      if (lpencil(i_glnmumol)) p%glnmumol(:,:)=0.0
!
!  THE FOLLOWING 2 ARE CONCEPTUALLY WRONG
!  FOR pretend_lnTT since iss actually contain lnTT NOT entropy!
!  The code is not wrong however since this is correctly
!  handled by the eos module.
!
      select case (ieosvars)
!
!  Work out thermodynamic quantities for given lnrho or rho and ss.
!
      case (ilnrho_ss,irho_ss)
        if (leos_isentropic) then
          if (lpencil(i_ss)) p%ss=0.0
          if (lpencil(i_gss)) p%gss=0.0
          if (lpencil(i_hss)) p%hss=0.0
          if (lpencil(i_del2ss)) p%del2ss=0.0
          if (lpencil(i_del6ss)) p%del6ss=0.0
          if (lpencil(i_cs2)) p%cs2=cs20*exp(gamma_m1*(p%lnrho-lnrho0))
        elseif (leos_isothermal) then
          if (lpencil(i_ss)) p%ss=-(cp-cv)*(p%lnrho-lnrho0)
          if (lpencil(i_gss)) p%gss=-(cp-cv)*p%glnrho
          if (lpencil(i_hss)) p%hss=-(cp-cv)*p%hlnrho
          if (lpencil(i_del2ss)) p%del2ss=-(cp-cv)*p%del2lnrho
          if (lpencil(i_del6ss)) p%del6ss=-(cp-cv)*p%del6lnrho
          if (lpencil(i_cs2)) p%cs2=cs20
        elseif (leos_localisothermal) then
          call fatal_error('calc_pencils_eos','leos_localisothermal '// &
              'not implemented for ilnrho_ss, try ilnrho_cs2')
        else
          if (lpencil(i_ss)) p%ss=f(l1:l2,m,n,ieosvar2)
          if (lpencil(i_gss)) call grad(f,ieosvar2,p%gss)
          if (lpencil(i_hss)) call g2ij(f,ieosvar2,p%hss)
          if (lpencil(i_del2ss)) call del2(f,ieosvar2,p%del2ss)
          if (lpencil(i_del6ss)) call del6(f,ieosvar2,p%del6ss)
          if (lpencil(i_cs2)) p%cs2=cs20*exp(cv1*p%ss+gamma_m1*(p%lnrho-lnrho0))
        endif
        if (lpencil(i_lnTT)) p%lnTT=lnTT0+cv1*p%ss+gamma_m1*(p%lnrho-lnrho0)
        if (lpencil(i_pp)) p%pp=(cp-cv)*exp(p%lnTT+p%lnrho)
        if (lpencil(i_ee)) p%ee=cv*exp(p%lnTT)
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
        if (lpencil(i_glnTT)) p%glnTT=gamma_m1*p%glnrho+cv1*p%gss
        if (lpencil(i_gTT)) then
          do i=1,3; p%gTT(:,i)=p%glnTT(:,i)*p%TT; enddo
        endif
        if (lpencil(i_del2lnTT)) p%del2lnTT=gamma_m1*p%del2lnrho+cv1*p%del2ss
        if (lpencil(i_hlnTT)) p%hlnTT=gamma_m1*p%hlnrho+cv1*p%hss
!
!  Work out thermodynamic quantities for given lnrho or rho and lnTT.
!
      case (ilnrho_lnTT,irho_lnTT)
        if (leos_isentropic) then
          if (lpencil(i_lnTT)) p%lnTT=gamma_m1*(p%lnrho-lnrho0)+lnTT0
          if (lpencil(i_glnTT)) p%glnTT=gamma_m1*p%glnrho
          if (lpencil(i_hlnTT)) p%hlnTT=gamma_m1*p%hlnrho
          if (lpencil(i_del2lnTT)) p%del2lnTT=gamma_m1*p%del2lnrho
          if (lpencil(i_cs2)) p%cs2=cs20*exp(gamma_m1*(p%lnrho-lnrho0))
        elseif (leos_isothermal) then
          if (lpencil(i_lnTT)) p%lnTT=lnTT0
          if (lpencil(i_glnTT)) p%glnTT=0.0
          if (lpencil(i_hlnTT)) p%hlnTT=0.0
          if (lpencil(i_del2lnTT)) p%del2lnTT=0.0
          if (lpencil(i_cs2)) p%cs2=cs20
        elseif (leos_localisothermal) then
          call fatal_error('calc_pencils_eos','leos_localisothermal '// &
              'not implemented for ilnrho_ss, try ilnrho_cs2')
        else
          if (lpencil(i_lnTT)) p%lnTT=f(l1:l2,m,n,ieosvar2)
          if (lpencil(i_glnTT)) call grad(f,ieosvar2,p%glnTT)
          if (lpencil(i_hlnTT)) call g2ij(f,ieosvar2,p%hlnTT)
          if (lpencil(i_del2lnTT)) call del2(f,ieosvar2,p%del2lnTT)
          if (lpencil(i_del6lnTT)) call del6(f,ieosvar2,p%del6lnTT)
          if (lpencil(i_cs2)) p%cs2=cp*exp(p%lnTT)*gamma_m1
        endif
        if (lpencil(i_ss)) p%ss=cv*(p%lnTT-lnTT0-gamma_m1*(p%lnrho-lnrho0))
        if (lpencil(i_pp)) p%pp=(cp-cv)*exp(p%lnTT+p%lnrho)
        if (lpencil(i_ee)) p%ee=cv*exp(p%lnTT)
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
        if (lpencil(i_gss)) p%gss=cv*(p%glnTT-gamma_m1*p%glnrho)
        if (lpencil(i_del2ss)) p%del2ss=cv*(p%del2lnTT-gamma_m1*p%del2lnrho)
        if (lpencil(i_hss)) p%hss=cv*(p%hlnTT-gamma_m1*p%hlnrho)
        if (lpencil(i_gTT)) then
          do i=1,3; p%gTT(:,i)=p%TT*p%glnTT(:,i); enddo
        endif
        if (lpencil(i_del6ss)) call fatal_error('calc_pencils_eos', &
            'del6ss not available for ilnrho_lnTT')
!
!  Work out thermodynamic quantities for given lnrho or rho and TT.
!
      case (ilnrho_TT,irho_TT)
        if (lpencil(i_TT))   p%TT=f(l1:l2,m,n,ieosvar2)
        if (lpencil(i_TT1).or.lpencil(i_hlnTT))  p%TT1=1/f(l1:l2,m,n,ieosvar2)
        if (lpencil(i_lnTT).or.lpencil(i_ss).or.lpencil(i_ee)) &
            p%lnTT=log(f(l1:l2,m,n,ieosvar2))
        if (lpencil(i_cs2))  p%cs2=cp*gamma_m1*f(l1:l2,m,n,ieosvar2)
        if (lpencil(i_gTT))  call grad(f,ieosvar2,p%gTT)
        if (lpencil(i_glnTT).or.lpencil(i_hlnTT)) then
          do i=1,3; p%glnTT(:,i)=p%gTT(:,i)*p%TT1; enddo
        endif
        if (lpencil(i_del2TT).or.lpencil(i_del2lnTT)) &
            call del2(f,ieosvar2,p%del2TT)
        if (lpencil(i_del2lnTT)) then
          tmp=0.0
          do i=1,3
            tmp=tmp+p%glnTT(:,i)**2
          enddo
          p%del2lnTT=p%del2TT*p%TT1-tmp
        endif
        if (lpencil(i_hlnTT)) then
          call g2ij(f,iTT,p%hlnTT)
          do i=1,3; do j=1,3
            p%hlnTT(:,i,j)=p%hlnTT(:,i,j)*p%TT1-p%glnTT(:,i)*p%glnTT(:,j)
          enddo; enddo
        endif
        if (lpencil(i_del6TT)) call del6(f,ieosvar2,p%del6TT)
        if (lpencil(i_ss)) p%ss=cv*(p%lnTT-lnTT0-gamma_m1*(p%lnrho-lnrho0))
        if (lpencil(i_pp)) p%pp=cv*gamma_m1*p%rho*p%TT
        if (lpencil(i_ee)) p%ee=cv*exp(p%lnTT)
!
!  Work out thermodynamic quantities for given lnrho or rho and cs2.
!
      case (ilnrho_cs2,irho_cs2)
        if (leos_isentropic) then
          call fatal_error('calc_pencils_eos', &
              'leos_isentropic not implemented for ilnrho_cs2, try ilnrho_ss')
        elseif (leos_isothermal) then
          if (lpencil(i_cs2)) p%cs2=cs20
          if (lpencil(i_lnTT)) p%lnTT=lnTT0
          if (lpencil(i_glnTT)) p%glnTT=0.0
          if (lpencil(i_hlnTT)) p%hlnTT=0.0
          if (lpencil(i_del2lnTT)) p%del2lnTT=0.0
          if (lpencil(i_ss)) p%ss=-(cp-cv)*(p%lnrho-lnrho0)
          if (lpencil(i_del2ss)) p%del2ss=-(cp-cv)*p%del2lnrho
          if (lpencil(i_gss)) p%gss=-(cp-cv)*p%glnrho
          if (lpencil(i_hss)) p%hss=-(cp-cv)*p%hlnrho
          if (lpencil(i_pp)) p%pp=gamma1*p%rho*cs20
        elseif (leos_localisothermal) then
          if (lpencil(i_cs2)) p%cs2=f(l1:l2,m,n,iglobal_cs2)
          if (lpencil(i_lnTT)) call fatal_error('calc_pencils_eos', &
              'temperature not needed for localisothermal')
          if (lpencil(i_glnTT)) &
              p%glnTT=f(l1:l2,m,n,iglobal_glnTT:iglobal_glnTT+2)
          if (lpencil(i_hlnTT)) call fatal_error('calc_pencils_eos', &
              'no gradients yet for localisothermal')
          if (lpencil(i_del2lnTT)) call fatal_error('calc_pencils_eos', &
              'no gradients yet for localisothermal')
          if (lpencil(i_ss)) call fatal_error('calc_pencils_eos', &
              'entropy not needed for localisothermal')
          if (lpencil(i_del2ss)) call fatal_error('calc_pencils_eos', &
              'no gradients yet for localisothermal')
          if (lpencil(i_gss)) call fatal_error('calc_pencils_eos', &
              'entropy gradient not needed for localisothermal')
          if (lpencil(i_hss)) call fatal_error('calc_pencils_eos', &
              'no gradients yet for localisothermal')
          if (lpencil(i_pp)) p%pp=p%rho*p%cs2
        else
          call fatal_error('calc_pencils_eos', &
              'Full equation of state not implemented for ilnrho_cs2')
        endif
!
!  Work out thermodynamic quantities for given pp and ss (anelastic case).
!
      case (ipp_ss)
        if (lanelastic) then
          if (lanelastic_lin) then
            p%pp=f(l1:l2,m,n,ipp)
            p%ss=f(l1:l2,m,n,iss)
            p%TTb=cs20*cp1*exp(gamma*f(l1:l2,m,n,iss_b)*cp1+gamma_m1*p%lnrho)/gamma_m1
            p%cs2=cp*p%TTb*gamma_m1
            p%TT1=1./p%TTb
            p%rho_anel=(f(l1:l2,m,n,ipp)/(f(l1:l2,m,n,irho_b)*p%cs2)- &
                 f(l1:l2,m,n,iss)*cp1)
          else
            if (lpencil(i_pp)) p%pp=f(l1:l2,m,n,ipp)
            if (lpencil(i_ss)) p%ss=f(l1:l2,m,n,iss)
            if (lpencil(i_rho)) p%rho=f(l1:l2,m,n,irho)
            !if (lpencil(i_rho)) p%rho=rho0*(gamma*p%pp/(rho0*cs20*exp(cv1*p%ss)))**gamma1
            if (lpencil(i_lnrho)) p%lnrho=alog(p%rho)
            if (lpencil(i_lnTT)) p%lnTT=lnTT0+cv1*p%ss+gamma_m1*(p%lnrho-lnrho0)
            if (lpencil(i_ee)) p%ee=cv*exp(p%lnTT)
            if (lpencil(i_yH)) p%yH=impossible
            if (lpencil(i_TT)) p%TT=exp(p%lnTT)
            if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
          endif
        endif
        if (leos_isentropic) then
          if (lpencil(i_ss)) p%ss=0.0
          if (lpencil(i_lnrho)) p%lnrho=log(gamma*p%pp/(rho0*cs20))/gamma
          if (lpencil(i_rho)) p%rho=exp(log(gamma*p%pp/(rho0*cs20))/gamma)
          if (lpencil(i_TT)) p%TT=(p%pp/pp0)**(1.-gamma1)
          if (lpencil(i_lnTT)) p%lnTT=(1.-gamma1)*log(gamma*p%pp/(rho0*cs0))
          if (lpencil(i_cs2)) p%cs2=cs20*(p%pp/pp0)**(1.-gamma1)
        elseif (leos_isothermal) then
          if (lpencil(i_lnrho)) p%lnrho=log(gamma*p%pp/(cs20*rho0))-p%lnTT
          if (lpencil(i_rho)) p%rho=exp(p%lnrho)
          if (lpencil(i_cs2)) p%cs2=cs20
          if (lpencil(i_lnTT)) p%lnTT=lnTT0
          if (lpencil(i_glnTT)) p%glnTT=0.0
          if (lpencil(i_hlnTT)) p%hlnTT=0.0
          if (lpencil(i_del2lnTT)) p%del2lnTT=0.0
        elseif (leos_localisothermal) then
          call fatal_error('calc_pencils_eos', &
              'Local Isothermal case not implemented for ipp_ss')
        endif
!
      case (ipp_cs2)
        if (leos_isentropic) then
          call fatal_error('calc_pencils_eos', &
              'isentropic not implemented for (pp,lnTT)')
        elseif (leos_isothermal) then
        if (lanelastic) then
          if (lanelastic_lin) then
            p%pp=f(l1:l2,m,n,ipp)
            p%rho_anel=f(l1:l2,m,n,ipp)/(f(l1:l2,m,n,irho_b)*cs20)
          else  ! lanelastic_lin=F means the non-linearized anelastic approx.
            p%pp=f(l1:l2,m,n,ipp)
          endif
        else
          if (lpencil(i_cs2)) p%cs2=cs20
          if (lpencil(i_lnrho)) p%lnrho=log(p%pp/cs20)
          if (lpencil(i_rho)) p%rho=(p%pp/cs20)
          if (lpencil(i_lnTT)) p%lnTT=lnTT0
          if (lpencil(i_glnTT)) p%glnTT=0.0
          if (lpencil(i_hlnTT)) p%hlnTT=0.0
          if (lpencil(i_del2lnTT)) p%del2lnTT=0.0
        endif
        elseif (leos_localisothermal) then
          call fatal_error('calc_pencils_eos', &
              'Local Isothermal case not implemented for ipp_cs2')
        endif
!
!  Internal energy.
!  For gamma=1, we use R/mu = c_p = c_v, thus ee = c_vT = R/mu T = p/rho = cs^2.
!
        if (lpencil(i_ee)) then
          if (gamma_m1/=0.0) then
            p%ee=(gamma1/gamma_m1)*p%cs2
          else
            p%ee=p%cs2
          endif
        endif
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
        if (lpencil(i_del6ss)) call fatal_error('calc_pencils_eos', &
            'del6ss not available for ilnrho_cs2')
!
!  Work out thermodynamic quantities for given lnrho or rho and eth.
!
      case (irho_eth,ilnrho_eth)
        if (lpencil(i_eth)) p%eth=f(l1:l2,m,n,ieth)
        if (lpencil(i_geth)) call grad(f,ieth,p%geth)
        if (lpencil(i_del2eth)) call del2(f,ieth,p%del2eth)
        if (lpencil(i_cs2)) p%cs2=gamma*gamma_m1*p%eth*p%rho1
        if (lpencil(i_pp)) p%pp=gamma_m1*p%eth
        if (lpencil(i_ee)) p%ee=p%rho1*p%eth
        if (lpencil(i_TT)) p%TT=p%cv1*p%rho1*p%eth
        if (lpencil(i_lnTT)) p%lnTT=alog(p%TT)
        if (lpencil(i_TT1)) p%TT1=1/p%TT
        if (lpencil(i_gTT)) then
          do i=1,3
            p%gTT(:,i)=p%rho1*(p%cv1*p%geth(:,i)-p%TT*p%grho(:,i))
          enddo
        endif
        if (lpencil(i_del2TT)) p%del2TT= &
            p%rho1*(p%cv1*p%del2eth-p%TT*p%del2rho-2*sum(p%grho*p%gTT,2))
!
!  Work out thermodynamic quantities for given rho and rhop. 
!  Check Lyra & Kuchner 2012 (arXiv:1204.6322) for details. 
!
      case (irho_rhop)
        if (lpencil(i_TT)) then
          p%TT=TT0*rho01*p%rhop
          if (TT_floor /= impossible) &
               where(p%TT < TT_floor) p%TT = TT_floor
        endif
        if (lpencil(i_lnTT)) p%lnTT=log(p%TT)
        if (lpencil(i_pp))   p%pp=cv*gamma_m1*p%rho*p%TT
        if (lpencil(i_ss))   p%ss=cv*(p%lnTT-lnTT0-gamma_m1*(p%lnrho-lnrho0))
        if (lpencil(i_cs2))  p%cs2=p%TT*p%cp*gamma_m1
      case default
        call fatal_error('calc_pencils_eos','case not implemented yet')
      endselect
!
    endsubroutine calc_pencils_eos
!***********************************************************************
    subroutine ioninit(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
    real, dimension (mx,my,mz,mfarray) :: f
!
    call keep_compiler_quiet(f)
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)
!
      real, intent(in) :: EE,TT,yH
      real, intent(inout) :: rho
!
      rho = EE * cv1 / TT
      call keep_compiler_quiet(yH)
!
    endsubroutine getdensity
!***********************************************************************
  subroutine gettemperature(f,TT_tmp)
!
     real, dimension (mx,my,mz,mfarray) :: f
     real, dimension (mx,my,mz), intent(out) :: TT_tmp
!
     call keep_compiler_quiet(f)
     call keep_compiler_quiet(TT_tmp)
!
   endsubroutine gettemperature
!***********************************************************************
   subroutine getpressure(pp_tmp)
!
     real, dimension (mx,my,mz), intent(out) :: pp_tmp
!
     call keep_compiler_quiet(pp_tmp)
!
   endsubroutine getpressure
!***********************************************************************
    subroutine get_cp1(cp1_)
!
!  04-nov-06/axel: added to alleviate spurious use of pressure_gradient
!
!  return the value of cp1 to outside modules
!
      real, intent(out) :: cp1_
!
      cp1_=cp1
!
    endsubroutine get_cp1
!***********************************************************************
    subroutine get_cv1(cv1_)
!
!  22-dec-10/PJK: adapted from get_cp1
!
!  return the value of cv1 to outside modules
!
      real, intent(out) :: cv1_
!
      cv1_=cv1
!
    endsubroutine get_cv1
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: cs2,cp1tilde
      real, dimension(nx) :: lnrho,ss
!
      if (ldensity_nolog) then
        lnrho=log(f(l1:l2,m,n,irho))
      else
        lnrho=f(l1:l2,m,n,ilnrho)
      endif
      ss=f(l1:l2,m,n,iss)
!
!  pretend_lnTT
!
      if (pretend_lnTT) then
        cs2=gamma_m1*exp(cv1*ss)
      else
        cs2=cs20*exp(cv1*ss+gamma_m1*(lnrho-lnrho0))
      endif
!! Actual pressure gradient calculation:
!!          do j=1,3
!!            ju=j+iuu-1
!!            if (pretend_lnTT) then
!!              df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) - &
!!                  p%cs2*(p%glnrho(:,j)/gamma + p%cp1tilde*p%gss(:,j))
!!            else
!!              df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) - &
!!                  p%cs2*(p%glnrho(:,j) + p%cp1tilde*p%gss(:,j))
!!            endif
!!           enddo
!
!  inverse cp (will be different from 1 when cp is not 1)
!
      cp1tilde=cp1
!
    endsubroutine pressure_gradient_farray
!***********************************************************************
    subroutine pressure_gradient_point(lnrho,ss,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: cs2,cp1tilde
!
!  pretend_lnTT
!
      if (pretend_lnTT) then
        cs2=gamma_m1*exp(gamma*cp1*ss)
      else
        cs2=cs20*exp(cv1*ss+gamma_m1*(lnrho-lnrho0))
      endif
      cp1tilde=cp1
!
    endsubroutine pressure_gradient_point
!***********************************************************************
    subroutine temperature_gradient(f,glnrho,gss,glnTT)
!
!   Calculate thermodynamical quantities
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(in) :: glnrho,gss
      real, dimension(nx,3), intent(out) :: glnTT
!
      if (gamma_m1==0.) call fatal_error('temperature_gradient', &
        'gamma=1 not allowed with entropy turned on!')
!
!  pretend_lnTT
!
      if (pretend_lnTT) then
        glnTT=gss
      else
        glnTT=gamma_m1*glnrho+cv1*gss
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_laplacian(f,p)
!
!   Calculate thermodynamical quantities
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      use Sub, only: dot2
!
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx) :: tmp
!
      if (gamma_m1==0.) &
           call fatal_error('temperature_laplacian', &
               'gamma=1 not allowed w/entropy')
!
!  pretend_lnTT
!
      if (pretend_lnTT) then
        p%del2lnTT=p%del2ss
      else
        if (ldensity_nolog) then
          call dot2(p%grho,tmp)
          p%del2lnTT=gamma_m1*p%rho1*(p%del2rho+p%rho1*tmp)
        else
          p%del2lnTT=gamma_m1*p%del2lnrho+p%cv1*p%del2ss
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine temperature_laplacian
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1
!   and optionally hlnPP and hlnTT
!   hP/rho=cs2*(hlnrho+cp1*hss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3,3), intent(in) :: hlnrho,hss
      real, dimension(nx,3,3), intent(out) :: hlnTT
!
      if (gamma_m1==0.) call fatal_error('temperature_hessian','gamma=1 not allowed w/entropy')
!
!  pretend_lnTT
!
      if (pretend_lnTT) then
        hlnTT=hss
      else
        hlnTT=gamma_m1*hlnrho+cv1*hss
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine thermal_energy_hessian(f,ivar_eth,del2lneth,hlneth)
!
      use Sub, only: g2ij,grad,dot2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: del2lneth,del2eth,geth2,eth_1
      real, dimension (nx,3,3) :: hlneth,heth
      real, dimension (nx,3) :: geth
      integer :: ivar_eth,i,j
!
      intent (in) :: f,ivar_eth
      intent (out) :: del2lneth,hlneth
!
      call g2ij(f,ivar_eth,heth)
      call grad(f,ivar_eth,geth)
!
      call dot2(geth,geth2)
!
      del2eth = heth(:,1,1) + heth(:,2,2) + heth(:,3,3)
!
      eth_1 = 1./f(l1:l2,m,n,ivar_eth)
!
      del2lneth = eth_1*del2eth - eth_1*eth_1*geth2
!
      do i=1,3
        do j=1,3
          hlneth(:,i,j) = eth_1*(heth(:,i,j) - eth_1*geth(:,i)*geth(:,j))
        enddo
      enddo
!
    endsubroutine thermal_energy_hessian
!***********************************************************************
    subroutine eosperturb(f,psize,ee,pp,ss)
!
!  Set f(l1:l2,m,n,iss), depending on the valyes of ee and pp
!  Adding pressure perturbations is not implemented
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(in), optional :: ee, pp, ss
      real, dimension(psize) :: lnrho_
!
      if (psize==nx) then
        if (ldensity_nolog) then
          lnrho_=log(f(l1:l2,m,n,irho))
        else
          lnrho_=f(l1:l2,m,n,ilnrho)
        endif
        if (present(ee)) then
          if (pretend_lnTT) then
            f(l1:l2,m,n,iss)=log(cv1*ee)
          else
            f(l1:l2,m,n,iss)=cv*(log(cv1*ee)-lnTT0-gamma_m1*(lnrho_-lnrho0))
          endif
        elseif (present(pp)) then
          if (pretend_lnTT) then
            f(l1:l2,m,n,iss)=log(gamma*pp/(gamma_m1*lnrho_))
          else
            f(l1:l2,m,n,iss)=cv*(log(gamma*pp/gamma_m1)-gamma*lnrho_-gamma_m1*lnrho0-lnTT0)
          endif
        elseif (present(ss)) then
          if (pretend_lnTT) then
            f(l1:l2,m,n,iss)=lnTT0+cv1*ss+gamma_m1*(lnrho_-lnrho0)
          else
            f(l1:l2,m,n,iss)=ss
          endif
        endif
!
      elseif (psize==mx) then
        if (ldensity_nolog) then
          lnrho_=log(f(:,m,n,irho))
        else
          lnrho_=f(:,m,n,ilnrho)
        endif
        if (present(ee)) then
          if (pretend_lnTT) then
            f(:,m,n,iss)=log(cv1*ee)
          else
            f(:,m,n,iss)=cv*(log(cv1*ee)-lnTT0-gamma_m1*(lnrho_-lnrho0))
          endif
        elseif (present(pp)) then
          if (pretend_lnTT) then
            f(:,m,n,iss)=log(gamma*pp/(gamma_m1*lnrho_))
          else
            f(:,m,n,iss)=cv*(log(gamma*pp/gamma_m1)-gamma*lnrho_-gamma_m1*lnrho0-lnTT0)
          endif
        elseif (present(ss)) then
          if (pretend_lnTT) then
            f(:,m,n,iss)=lnTT0+cv1*ss+gamma_m1*(lnrho_-lnrho0)
          else
            f(:,m,n,iss)=ss
          endif
        endif
!
      else
        call not_implemented("eosperturb")
      endif
    endsubroutine eosperturb
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,yH,lnTT,ee,pp,kapparho)
!
!   Calculate thermodynamical quantities
!
!   02-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1 to
!                   subroutine pressure_gradient
!
      use Diagnostics, only: max_mn_name, sum_mn_name
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(out), optional :: lnrho
      real, dimension(psize), intent(out), optional :: yH,ee,pp,kapparho
      real, dimension(psize), intent(out), optional :: lnTT
      real, dimension(psize) :: lnTT_, cs2_
      real, dimension(psize) :: lnrho_,ss_
!
!ajwm this test should be done at initialization
!      if (gamma_m1==0.) call fatal_error('eoscalc_farray','gamma=1 not allowed w/entropy')
!
      select case (ieosvars)
!
! Log rho and entropy
!
      case (ilnrho_ss,irho_ss)
        select case (psize)
        case (nx)
          if (ieosvars==ilnrho_ss) then
            lnrho_=f(l1:l2,m,n,ieosvar1)
          else
            lnrho_=log(f(l1:l2,m,n,ieosvar1))
          endif
          if (leos_isentropic) then
            ss_=0
          elseif (leos_isothermal) then
            ss_=-cv*gamma_m1*(lnrho_-lnrho0)
          else
            ss_=f(l1:l2,m,n,ieosvar2)
          endif
        case (mx)
          if (ieosvars==ilnrho_ss) then
            lnrho_=f(:,m,n,ieosvar1)
          else
            lnrho_=log(f(:,m,n,ieosvar1))
          endif
          if (leos_isentropic) then
            ss_=0
          elseif (leos_isothermal) then
            ss_=-cv*gamma_m1*(lnrho_-lnrho0)
          else
            ss_=f(:,m,n,ieosvar2)
          endif
        case default
          call fatal_error('eoscalc_farray','no such pencil size')
        end select
!
        lnTT_=lnTT0+cv1*ss_+gamma_m1*(lnrho_-lnrho0)
        if (gamma_m1==0.) &
            call fatal_error('eoscalc_farray','gamma=1 not allowed w/entropy')
        if (present(lnrho)) lnrho=lnrho_
        if (present(lnTT)) lnTT=lnTT_
        if (present(ee)) ee=cv*exp(lnTT_)
        if (present(pp)) pp=(cp-cv)*exp(lnTT_+lnrho_)
!
! Log rho and Log T
!
      case (ilnrho_lnTT,irho_lnTT)
        select case (psize)
        case (nx)
          if (ieosvars==ilnrho_lnTT) then
            lnrho_=f(l1:l2,m,n,ieosvar1)
          else
            lnrho_=log(f(l1:l2,m,n,ieosvar1))
          endif
          if (leos_isentropic) then
            lnTT_=lnTT0+(cp-cv)*(lnrho_-lnrho0)
          elseif (leos_isothermal) then
            lnTT_=lnTT0
          else
            lnTT_=f(l1:l2,m,n,ieosvar2)
          endif
        case (mx)
          if (ieosvars==ilnrho_lnTT) then
            lnrho_=f(:,m,n,ieosvar1)
          else
            lnrho_=log(f(:,m,n,ieosvar1))
          endif
          if (leos_isentropic) then
            lnTT_=lnTT0+(cp-cv)*(lnrho_-lnrho0)
          elseif (leos_isothermal) then
            lnTT_=lnTT0
          else
            lnTT_=f(:,m,n,ieosvar2)
          endif
        case default
          call fatal_error('eoscalc_farray','no such pencil size')
        end select
!
        if (present(lnrho)) lnrho=lnrho_
        if (present(lnTT)) lnTT=lnTT_
        if (present(ee)) ee=cv*exp(lnTT_)
        if (present(pp)) pp=(cp-cv)*exp(lnTT_+lnrho_)
!
! Log rho or rho and T
!
      case (ilnrho_TT,irho_TT)
          call fatal_error('eoscalc_farray','no implemented for lnrho_TT or rho_TT')
!
! Log rho and cs2
!
      case (ilnrho_cs2,irho_cs2)
        select case (psize)
        case (nx)
          if (ieosvars==ilnrho_cs2) then
            lnrho_=f(l1:l2,m,n,ieosvar1)
          else
            lnrho_=log(f(l1:l2,m,n,ieosvar1))
          endif
          if (leos_isentropic) then
            cs2_=exp(gamma_m1*(lnrho_-lnrho0)+log(cs20))
          elseif (leos_isothermal) then
            cs2_=cs20
          elseif (leos_localisothermal) then
            cs2_=f(l1:l2,m,n,iglobal_cs2)
          else
            call fatal_error('eoscalc_farray','full eos for cs2 not implemented')
          endif
        case (mx)
          if (ieosvars==ilnrho_cs2) then
            lnrho_=f(:,m,n,ieosvar1)
          else
            lnrho_=log(f(:,m,n,ieosvar1))
          endif
          if (leos_isentropic) then
            cs2_=exp(gamma_m1*(lnrho_-lnrho0)+log(cs20))
          elseif (leos_isothermal) then
            cs2_=cs20
          elseif (leos_localisothermal) then
            cs2_=f(:,m,n,iglobal_cs2)
          else
            call fatal_error('eoscalc_farray','full eos for cs2 not implemented')
          endif
        case default
          call fatal_error('eoscalc_farray','no such pencil size')
        end select
!
        if (present(lnrho)) lnrho=lnrho_
        if (present(lnTT)) lnTT=lnTT0+log(cs2_)
        if (present(ee)) ee=gamma1*cs2_/gamma_m1
        if (present(pp)) pp=gamma1*cs2_*exp(lnrho_)
!
      case default
        call fatal_error("eoscalc_farray",'Thermodynamic variable combination not implemented!')
      endselect
!
      if (present(yH)) yH=impossible
!
      if (present(kapparho)) then
        kapparho=0
        call fatal_error("eoscalc","sorry, no Hminus opacity with noionization")
      endif
!
    endsubroutine eoscalc_farray
!***********************************************************************
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1 to
!                   subroutine pressure_gradient
!   27-mar-06/tony: Introduces cv, cv1, gamma1 to make faster
!                   + more explicit
!   31-mar-06/tony: I removed messy lcalc_cp stuff completely. cp=1.
!                   is just fine.
!   22-jun-06/axel: reinstated cp,cp1,cv,cv1 in hopefully all the places.
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp,cs2
      real :: lnrho_,ss_,lnTT_,ee_,pp_,cs2_,TT_
!
      if (gamma_m1==0.and..not.lanelastic) call fatal_error &
        ('eoscalc_point','gamma=1 not allowed w/entropy')
!
      select case (ivars)
!
      case (ilnrho_ss,irho_ss)
        if (ivars==ilnrho_ss) then
          lnrho_=var1
        else
          lnrho_=log(var1)
        endif
        ss_=var2
        lnTT_=lnTT0+cv1*ss_+gamma_m1*(lnrho_-lnrho0)
        ee_=cv*exp(lnTT_)
        pp_=(cp-cv)*exp(lnTT_+lnrho_)
        cs2_=gamma*gamma_m1*ee_
!
      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
        lnTT_=log(cv1*ee_)
        ss_=cv*(lnTT_-lnTT0-gamma_m1*(lnrho_-lnrho0))
        pp_=gamma_m1*ee_*exp(lnrho_)
        cs2_=gamma*gamma_m1*ee_
!
      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
        ss_=cv*(log(pp_*exp(-lnrho_)*gamma/cs20)-gamma_m1*(lnrho_-lnrho0))
        ee_=pp_*exp(-lnrho_)/gamma_m1
        lnTT_=log(cv1*ee_)
        cs2_=gamma*gamma_m1*ee_
!
      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
        ss_=cv*(lnTT_-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*exp(lnTT_)
        pp_=ee_*exp(lnrho_)*gamma_m1
        cs2_=gamma*gamma_m1*ee_
!
      case (ilnrho_TT)
        lnrho_=var1
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*exp(lnrho_)*gamma_m1
        cs2_=cp*gamma_m1*TT_
!
      case (irho_TT)
        lnrho_=log(var1)
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*var1*gamma_m1
        cs2_=cp*gamma_m1*TT_
!
      case (ipp_cs2)
        if (lanelastic) then
          if (lanelastic_lin) then
            lnrho_=log(var1)
            TT_=exp(lnTT0)
            pp_=exp(lnrho_)*cs20/gamma
          else
            if (leos_isothermal) then
              pp_=var1
              lnrho_=log(pp_*cs20)
              TT_=exp(lnTT0)
            endif
          endif
        endif
!
      case (ipp_ss)
        if (lanelastic) then
          if (lanelastic_lin) then
            lnrho_=(var1)
            ss_=var2
            cs2_=exp(gamma*ss_*cp1+gamma_m1*(lnrho_-lnrho0))*cs20
            TT_=cs2_/(gamma_m1*cp)
          else
            pp_=var1
            ss_=var2
            cs2_=exp(ss_*cp1+gamma1*gamma_m1*log(pp_/pp0))*cs20
            TT_=cs2_/(gamma_m1*cp)
            lnrho_=log(gamma*pp_/cs2_)
          endif
        endif
      case default
        call not_implemented('eoscalc_point')
      end select
!
      if (present(lnrho)) lnrho=lnrho_
      if (present(ss)) ss=ss_
      if (present(yH)) yH=impossible
      if (present(lnTT)) lnTT=lnTT_
      if (present(ee)) ee=ee_
      if (present(pp)) pp=pp_
      if (present(cs2)) cs2=cs2_
!
    endsubroutine eoscalc_point
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1 to
!                   subroutine pressure_gradient
!   27-mar-06/tony: Introduces cv, cv1, gamma1 to make faster
!                   + more explicit
!   31-mar-06/tony: I removed messy lcalc_cp stuff completely. cp=1.
!                   is just fine.
!   22-jun-06/axel: reinstated cp,cp1,cv,cv1 in hopefully all the places.
!
      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), intent(out), optional :: lnrho,ss
      real, dimension(nx), intent(out), optional :: yH,lnTT
      real, dimension(nx), intent(out), optional :: ee,pp,cs2
      real, dimension(nx) :: lnrho_,ss_,lnTT_,ee_,pp_,cs2_,TT_
!
      if (gamma_m1==0.) call fatal_error('eoscalc_pencil','gamma=1 not allowed w/entropy')
!
      select case (ivars)
!
      case (ilnrho_ss,irho_ss)
        if (ivars==ilnrho_ss) then
          lnrho_=var1
        else
          lnrho_=log(var1)
        endif
        ss_=var2
        lnTT_=lnTT0+cv1*ss_+gamma_m1*(lnrho_-lnrho0)
        ee_=cv*exp(lnTT_)
        pp_=(cp-cv)*exp(lnTT_+lnrho_)
        cs2_=gamma*gamma_m1*ee_
        cs2_=cs20*cv1*ee_
!
      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
        lnTT_=log(cv1*ee_)
        ss_=cv*(lnTT_-lnTT0-gamma_m1*(lnrho_-lnrho0))
        pp_=gamma_m1*ee_*exp(lnrho_)
        cs2_=gamma*gamma_m1*ee_
!
      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
        ss_=cv*(log(pp_*exp(-lnrho_)*gamma/cs20)-gamma_m1*(lnrho_-lnrho0))
        ee_=pp_*exp(-lnrho_)/gamma_m1
        lnTT_=log(cv1*ee_)
        cs2_=gamma*gamma_m1*ee_
!
      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
        ss_=cv*(lnTT_-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*exp(lnTT_)
        pp_=ee_*exp(lnrho_)*gamma_m1
        cs2_=gamma*gamma_m1*ee_
!
      case (ilnrho_TT)
        lnrho_=var1
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*exp(lnrho_)*gamma_m1
        cs2_=cp*gamma_m1*TT_
!
      case (irho_TT)
        lnrho_=log(var1)
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*var1*gamma_m1
        cs2_=cp*gamma_m1*TT_
!DM+PC
      case (ipp_ss)
        pp_=var1
        ss_=var2
        lnrho_=log(pp_)/gamma-ss_/cp
        TT_=pp_/((gamma_m1)*cv*exp(lnrho_))
        cs2_=cp*gamma_m1*TT_
      case default
        call not_implemented('eoscalc_point')
      end select
!
      if (present(lnrho)) lnrho=lnrho_
      if (present(ss)) ss=ss_
      if (present(yH)) yH=impossible
      if (present(lnTT)) lnTT=lnTT_
      if (present(ee)) ee=ee_
      if (present(pp)) pp=pp_
      if (present(cs2)) cs2=cs2_
!
    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine get_soundspeed(TT,cs2)
!
!  Calculate sound speed for given temperature
!
!  20-Oct-03/tobi: Coded
!
      real, intent(in)  :: TT
      real, intent(out) :: cs2
!
      cs2=gamma_m1*cp*TT
!
    endsubroutine get_soundspeed
!***********************************************************************
    subroutine read_eos_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=eos_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=eos_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_eos_init_pars
!***********************************************************************
    subroutine write_eos_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=eos_init_pars)
!
    endsubroutine write_eos_init_pars
!***********************************************************************
    subroutine read_eos_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=eos_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=eos_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_eos_run_pars
!***********************************************************************
    subroutine write_eos_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=eos_run_pars)
!
    endsubroutine write_eos_run_pars
!***********************************************************************
    subroutine isothermal_entropy(f,T0)
!
!  Isothermal stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!  Sound speed (and hence Temperature), is
!  initialised to the reference value:
!           sound speed: cs^2_0            from start.in
!           density: rho0 = exp(lnrho0)
!
!  11-jun-03/tony: extracted from isothermal routine in Density module
!                  to allow isothermal condition for arbitrary density
!  17-oct-03/nils: works also with leos_ionization=T
!  18-oct-03/tobi: distributed across ionization modules
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: T0
      real, dimension(nx) :: lnrho,ss,lnTT
!      real :: ss_offset=0.
!
!  if T0 is different from unity, we interpret
!  ss_offset = ln(T0)/gamma as an additive offset of ss
!
!      if (T0/=1.) ss_offset=log(T0)/gamma
!
      do n=n1,n2
      do m=m1,m2
        if (ldensity_nolog) then
          lnrho=log(f(l1:l2,m,n,irho))
        else
          lnrho=f(l1:l2,m,n,ilnrho)
        endif
        lnTT=log(T0)
          !+ other terms for sound speed not equal to cs_0
        call eoscalc(ilnrho_lnTT,lnrho,lnTT,ss=ss)
        f(l1:l2,m,n,iss)=ss
      enddo
      enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      call get_soundspeed(T0,cs2bot)
      cs2top=cs2bot
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine isothermal_lnrho_ss(f,T0,rho0)
!
!  Isothermal stratification for lnrho and ss (for yH=0!)
!
!  Currently only implemented for ionization_fixed.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: T0,rho0
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(T0)
      call keep_compiler_quiet(rho0)
!
    endsubroutine isothermal_lnrho_ss
!***********************************************************************
    subroutine Hminus_opacity(f,kapparho)
!
!  dummy routine
!
!  03-apr-2004/tobi: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz), intent(out) :: kapparho
!
      call fatal_error('Hminus_opacity',"opacity_type='Hminus' may not be used with noionization")
!
      call keep_compiler_quiet(kapparho)
      call keep_compiler_quiet(f)
!
    endsubroutine Hminus_opacity
!***********************************************************************
    subroutine bc_ss_flux_orig(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      logical, pointer :: lmultilayer, lheatc_chiconst
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy
      integer :: i,ierr
!
      if (ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
!  Get the shared variables
!
      call get_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond0")
      call get_shared_variable('hcond1',hcond1,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond1")
      call get_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting Fbot")
      call get_shared_variable('Ftop',Ftop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting Ftop")
      call get_shared_variable('FbotKbot',FbotKbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting FbotKbot")
      call get_shared_variable('FtopKtop',FtopKtop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting FtopKtop")
      call get_shared_variable('chi',chi,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting chi")
      call get_shared_variable('lmultilayer',lmultilayer,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lmultilayer")
      call get_shared_variable('lheatc_chiconst',lheatc_chiconst,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lheatc_chiconst")
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case ('bot')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n1,ilnrho)-lnrho0)+cv1*f(:,:,n1,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!AB: are here any cp factors?
!
        if (lheatc_chiconst) then
          tmp_xy=Fbot/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FbotKbot/cs2_xy
        endif
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+(cp-cv)* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case ('top')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0
        endif
!
!  calculate Ftop/(K*cs2)
!
        rho_xy=exp(f(:,:,n2,ilnrho))
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n2,ilnrho)-lnrho0)+cv1*f(:,:,n2,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!
        if (lheatc_chiconst) then
          tmp_xy=Ftop/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FtopKtop/cs2_xy
        endif
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+(cp-cv)* &
              (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)-2*i*dz*tmp_xy)
        enddo
      case default
        call fatal_error('bc_ss_flux','invalid argument')
      endselect
!
    endsubroutine bc_ss_flux_orig
!***********************************************************************
    subroutine get_average_pressure(init_average_density,average_density,&
                                    average_pressure)
!
!   01-dec-2009/piyali+dhruba: coded
!
      real, intent(in) :: init_average_density,average_density
      real, intent(inout) :: average_pressure
!
      if (leos_isothermal.or.lfirst) then
        average_pressure = average_density*cs20
      else
        average_pressure = average_pressure+((average_density/&
                           init_average_density)**gamma-1.0)*pp0*pres_corr
        call fatal_error('get_average_pressure','Non isothermal case no coded yet')
      endif
!
    endsubroutine get_average_pressure
!***********************************************************************
    subroutine bc_ss_flux_tmp(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      logical, pointer :: lmultilayer, lheatc_chiconst
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy,lnrho_xy,ss_xy
      real, dimension (mx,my) :: cs2_xy1,cs2_xy2,T_xy,T_xy1,T_xy2
      real :: eps
      integer :: i,ierr,iter,j,k
      integer,parameter :: niter=4
!
      if (ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
!  Get the shared variables
!
      call get_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond0")
      call get_shared_variable('hcond1',hcond1,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond1")
      call get_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting Fbot")
      call get_shared_variable('Ftop',Ftop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting Ftop")
      call get_shared_variable('FbotKbot',FbotKbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting FbotKbot")
      call get_shared_variable('FtopKtop',FtopKtop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting FtopKtop")
      call get_shared_variable('chi',chi,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting chi")
      call get_shared_variable('lmultilayer',lmultilayer,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lmultilayer")
      call get_shared_variable('lheatc_chiconst',lheatc_chiconst,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lheatc_chiconst")
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case ('bot')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n1,ilnrho)-lnrho0)+cv1*f(:,:,n1,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!AB: are here any cp factors?
!
        if (lheatc_chiconst) then
          tmp_xy=Fbot/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FbotKbot/cs2_xy
        endif
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+(cp-cv)* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case ('top')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0
        endif
!
!  Compute Temperature at the first 3 levels inward
!
        lnrho_xy=f(:,:,n2,ilnrho)
        cs2_xy =cs20*exp(gamma_m1*(f(:,:,n2  ,ilnrho)-lnrho0)+cv1*f(:,:,n2  ,iss))
        cs2_xy1=cs20*exp(gamma_m1*(f(:,:,n2-1,ilnrho)-lnrho0)+cv1*f(:,:,n2-1,iss))
        cs2_xy2=cs20*exp(gamma_m1*(f(:,:,n2-2,ilnrho)-lnrho0)+cv1*f(:,:,n2-2,iss))
        T_xy=cs2_xy/(cp*gamma_m1)
        T_xy1=cs2_xy1/(cp*gamma_m1)
        T_xy2=cs2_xy2/(cp*gamma_m1)
!
!  calculate Ftop/(K*cs2)
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!
        if (lheatc_chiconst) then
          tmp_xy=Ftop/(exp(lnrho_xy)*chi*cs2_xy)
        else
          tmp_xy=FtopKtop/cs2_xy
        endif
!
!  iterate
!
        eps=2.*dz/1.5**4/3.
        do iter=1,niter
          T_xy=(4*T_xy1-T_xy2)/3.-eps*(T_xy/1.5)**4
          if (ip<8) print*,'iter, T_xy(l1,m1)=',iter,T_xy(l1,m1)
        enddo
!
!  use EOS to work out ss on the boundary
!
        call eoscalc_pencil(ilnrho_TT,lnrho_xy,T_xy,ss=ss_xy)
        f(:,:,n2,iss)=ss_xy
!
!  apply to ghost zones
!
          j=iss
          k=n2+1
          f(:,:,k,j)=7*f(:,:,k-1,j) &
                   -21*f(:,:,k-2,j) &
                   +35*f(:,:,k-3,j) &
                   -35*f(:,:,k-4,j) &
                   +21*f(:,:,k-5,j) &
                    -7*f(:,:,k-6,j) &
                      +f(:,:,k-7,j)
          k=n2+2
          f(:,:,k,j)=9*f(:,:,k-1,j) &
                   -35*f(:,:,k-2,j) &
                   +77*f(:,:,k-3,j) &
                  -105*f(:,:,k-4,j) &
                   +91*f(:,:,k-5,j) &
                   -49*f(:,:,k-6,j) &
                   +15*f(:,:,k-7,j) &
                    -2*f(:,:,k-8,j)
          k=n2+3
          f(:,:,k,j)=9*f(:,:,k-1,j) &
                   -45*f(:,:,k-2,j) &
                  +147*f(:,:,k-3,j) &
                  -315*f(:,:,k-4,j) &
                  +441*f(:,:,k-5,j) &
                  -399*f(:,:,k-6,j) &
                  +225*f(:,:,k-7,j) &
                   -72*f(:,:,k-8,j) &
                   +10*f(:,:,k-9,j)
!
      case default
        call fatal_error('bc_ss_flux','invalid argument')
      endselect
!
    endsubroutine bc_ss_flux_tmp
!***********************************************************************
    subroutine bc_ss_flux_tmp2(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      logical, pointer :: lmultilayer, lheatc_chiconst
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy
      integer :: i,ierr
!
      if (ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
!  Get the shared variables
!
      call get_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond0")
      call get_shared_variable('hcond1',hcond1,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond1")
      call get_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting Fbot")
      call get_shared_variable('Ftop',Ftop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting Ftop")
      call get_shared_variable('FbotKbot',FbotKbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting FbotKbot")
      call get_shared_variable('FtopKtop',FtopKtop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting FtopKtop")
      call get_shared_variable('chi',chi,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting chi")
      call get_shared_variable('lmultilayer',lmultilayer,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lmultilayer")
      call get_shared_variable('lheatc_chiconst',lheatc_chiconst,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lheatc_chiconst")
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case ('bot')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n1,ilnrho)-lnrho0)+cv1*f(:,:,n1,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!AB: are here any cp factors?
!
        if (lheatc_chiconst) then
          tmp_xy=Fbot/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FbotKbot/cs2_xy
        endif
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+(cp-cv)* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case ('top')
!
!  Set (dcs2/dz) / (dcs2/dz)_ini = (cs2/cs2top_ini)^4
!  Note that (dcs2/dz) = cs20*[(gamma-1)*dlnrho/dz + gamma*d(s/cp)/dz]
!  So, ds/dz = - (cp-cv)*dlnrho/dz + cv*(dcs2/dz)/cs20
!  calculate tmp_xy
!
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n2,ilnrho)-lnrho0)+cv1*f(:,:,n2,iss))
        tmp_xy=cv*dcs2top_ini/cs20*(cs2_xy/cs2top_ini)**4
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n2+i,iss)=f(:,:,n2-i,iss) &
              -(cp-cv)*(f(:,:,n2+i,ilnrho)-f(:,:,n2-i,ilnrho)) &
              +2*i*dz*tmp_xy
        enddo
      case default
        call fatal_error('bc_ss_flux','invalid argument')
      endselect
!
    endsubroutine bc_ss_flux_tmp2
!***********************************************************************
    subroutine bc_ss_flux(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      real, pointer :: hcond0_kramers, nkramers
      logical, pointer :: lmultilayer, lheatc_chiconst, lheatc_kramers
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy
      integer :: i,ierr
!
      if (ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
!  Get the shared variables
!
      call get_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond0")
      call get_shared_variable('hcond1',hcond1,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond1")
      call get_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting Fbot")
      call get_shared_variable('Ftop',Ftop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting Ftop")
      call get_shared_variable('FbotKbot',FbotKbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting FbotKbot")
      call get_shared_variable('FtopKtop',FtopKtop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting FtopKtop")
      call get_shared_variable('chi',chi,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting chi")
      call get_shared_variable('lmultilayer',lmultilayer,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lmultilayer")
      call get_shared_variable('lheatc_chiconst',lheatc_chiconst,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lheatc_chiconst")
      call get_shared_variable('hcond0_kramers',hcond0_kramers,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond0_kramers")
      call get_shared_variable('nkramers',nkramers,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting nkramers")
      call get_shared_variable('lheatc_kramers',lheatc_kramers,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lheatc_kramers")
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case ('bot')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        if (ldensity_nolog) then
          rho_xy=f(:,:,n1,irho)
          cs2_xy=cs20*exp(gamma_m1*(log(f(:,:,n1,irho))-lnrho0)+cv1*f(:,:,n1,iss))
        else
          rho_xy=exp(f(:,:,n1,ilnrho))
          cs2_xy=cs20*exp(gamma_m1*(f(:,:,n1,ilnrho)-lnrho0)+cv1*f(:,:,n1,iss))
        endif
!
!  Check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!  Check also whether Kramers opacity is used, then hcond itself depends
!  on density and temperature.
!
        if (lheatc_chiconst) then
          tmp_xy=Fbot/(rho_xy*chi*cs2_xy)
        else if (lheatc_kramers) then
          tmp_xy=Fbot*rho_xy**(2*nkramers)*(cp*gamma_m1)**(6.5*nkramers)/(hcond0_kramers*cs2_xy**(6.5*nkramers+1.))
        else
          tmp_xy=FbotKbot/cs2_xy
        endif
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          if (ldensity_nolog) then
            f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+(cp-cv)* &
              (log(f(:,:,n1+i,irho)/f(:,:,n1-i,irho))+2*i*dz*tmp_xy)
          else
            f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+(cp-cv)* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
          endif
        enddo
!
!  top boundary
!  ============
!
      case ('top')
!
!  calculate Fbot/(K*cs2)
!
        if (ldensity_nolog) then
          rho_xy=f(:,:,n2,irho)
          cs2_xy=cs20*exp(gamma_m1*(log(f(:,:,n2,irho))-lnrho0)+cv1*f(:,:,n2,iss))
        else
          rho_xy=exp(f(:,:,n2,ilnrho))
          cs2_xy=cs20*exp(gamma_m1*(f(:,:,n2,ilnrho)-lnrho0)+cv1*f(:,:,n2,iss))
        endif
!
!  Check whether we have chi=constant at top, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!  Check also whether Kramers opacity is used, then hcond itself depends
!  on density and temperature.
!
        if (lheatc_chiconst) then
          tmp_xy=Ftop/(rho_xy*chi*cs2_xy)
        else if (lheatc_kramers) then
          tmp_xy=Ftop*rho_xy**(2*nkramers)*(cp*gamma_m1)**(6.5*nkramers)/(hcond0_kramers*cs2_xy**(6.5*nkramers+1.))
        else
          tmp_xy=FtopKtop/cs2_xy
        endif
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          if (ldensity_nolog) then
            f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+(cp-cv)* &
              (log(f(:,:,n2-i,irho)/f(:,:,n2+i,irho))+2*i*dz*tmp_xy)
          else
            f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+(cp-cv)* &
              (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)-2*i*dz*tmp_xy)
          endif
        enddo
!
      case default
        call fatal_error('bc_ss_flux','invalid argument')
      endselect
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_flux_turb(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcz='Fgs')
!
!   04-may-2009/axel: adapted from bc_ss_flux
!   31-may-2010/pete: replaced sigmaSB by a `turbulent' sigmaSBt
!
      use Mpicomm, only: stop_it
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: chi,chi_t,hcondzbot,hcondztop
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: dsdz_xy,cs2_xy,rho_xy,TT_xy,dlnrhodz_xy
      real :: fac
      integer :: i,ierr
!
      if (ldebug) print*,'bc_ss_flux_turb: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
!  Get the shared variables
!
      call get_shared_variable('chi_t',chi_t,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux_turb: "//&
           "there was a problem when getting chi_t")
      call get_shared_variable('chi',chi,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux_turb: "//&
           "there was a problem when getting chi")
      call get_shared_variable('hcondzbot',hcondzbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux_turb: "//&
           "there was a problem when getting hcondzbot")
      call get_shared_variable('hcondztop',hcondztop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux_turb: "//&
           "there was a problem when getting hcondztop")
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case ('bot')
!
!  set ghost zones such that dsdz_xy obeys
!  - chi_t rho T dsdz_xy - hcond gTT = sigmaSBt*TT^4
!
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n1,ilnrho)-lnrho0)+cv1*f(:,:,n1,iss))
        rho_xy=exp(f(:,:,n1,ilnrho))
        TT_xy=cs2_xy/(gamma_m1*cp)
        fac=(1./60)*dz_1(n1)
        dlnrhodz_xy=fac*(+ 45.0*(f(:,:,n1+1,ilnrho)-f(:,:,n1-1,ilnrho)) &
                         -  9.0*(f(:,:,n1+2,ilnrho)-f(:,:,n1-2,ilnrho)) &
                         +      (f(:,:,n1+3,ilnrho)-f(:,:,n1-3,ilnrho)))
        dsdz_xy=-(sigmaSBt*TT_xy**3+hcondzbot*(gamma_m1)*dlnrhodz_xy)/ &
            (chi_t*rho_xy+hcondzbot/cv)
!
!  enforce ds/dz=-(sigmaSBt*T^3 + hcond*(gamma-1)*glnrho)/(chi_t*rho+hcond/cv)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+2*i*dz*dsdz_xy
        enddo
!
!  top boundary
!  ============
!
      case ('top')
!
!  set ghost zones such that dsdz_xy obeys
!  - chi_t rho T dsdz_xy - hcond gTT = sigmaSBt*TT^4
!
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n2,ilnrho)-lnrho0)+cv1*f(:,:,n2,iss))
        rho_xy=exp(f(:,:,n2,ilnrho))
        TT_xy=cs2_xy/(gamma_m1*cp)
        fac=(1./60)*dz_1(n2)
        dlnrhodz_xy=fac*(+ 45.0*(f(:,:,n2+1,ilnrho)-f(:,:,n2-1,ilnrho)) &
                         -  9.0*(f(:,:,n2+2,ilnrho)-f(:,:,n2-2,ilnrho)) &
                         +      (f(:,:,n2+3,ilnrho)-f(:,:,n2-3,ilnrho)))
        if (hcondztop==impossible) then
          dsdz_xy=-(sigmaSBt*TT_xy**3+chi*rho_xy*cp*(gamma_m1)*dlnrhodz_xy)/ &
              (chi_t*rho_xy+chi*rho_xy*cp/cv)
        else
          dsdz_xy=-(sigmaSBt*TT_xy**3+hcondztop*(gamma_m1)*dlnrhodz_xy)/ &
              (chi_t*rho_xy+hcondztop/cv)
        endif
!
!  enforce ds/dz=-(sigmaSBt*T^3 + hcond*(gamma-1)*glnrho)/(chi_t*rho+hcond/cv)
!
        do i=1,nghost
          f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+2*i*dz*dsdz_xy
        enddo
!
!  capture undefined entries
!
      case default
        call fatal_error('bc_ss_flux_turb','invalid argument')
      endselect
!
    endsubroutine bc_ss_flux_turb
!***********************************************************************
    subroutine bc_ss_flux_turb_x(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcz='Fgs')
!
!   31-may-2010/pete: adapted from bc_ss_flux_turb
!   20-jul-2010/pete: expanded to take into account hcond/=0
!
      use Mpicomm, only: stop_it
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: chi_t,hcondxbot,hcondxtop
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: dsdx_yz,cs2_yz,rho_yz,dlnrhodx_yz,TT_yz
      real :: fac
      integer :: i=0,ierr
!
      if (ldebug) print*,'bc_ss_flux_turb: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `Fgs' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
!  Get the shared variables
!
      call get_shared_variable('chi_t',chi_t,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux_turb_x: "//&
           "there was a problem when getting chi_t")
      call get_shared_variable('hcondxbot',hcondxbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux_turb_x: "//&
           "there was a problem when getting hcondxbot")
      call get_shared_variable('hcondxtop',hcondxtop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux_turb_x: "//&
           "there was a problem when getting hcondxtop")
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case ('bot')
!
! For the case of pretend_lnTT=T, set glnTT=-sigma*T^3/hcond
!
        if (pretend_lnTT) then
            f(l1-1,:,:,iss)=f(l1+i,:,:,iss) + &
                2*i*dx*sigmaSBt*exp(f(l1,:,:,iss))**3/hcondxbot
        else
!
!  set ghost zones such that dsdx_yz obeys
!  - chi_t rho T dsdx_yz - hcond gTT = sigmaSBt*TT^4
!
          cs2_yz=cs20*exp(gamma_m1*(f(l1,:,:,ilnrho)-lnrho0)+cv1*f(l1,:,:,iss))
          TT_yz=cs2_yz/(gamma_m1*cp)
          rho_yz=exp(f(l1,:,:,ilnrho))
          fac=(1./60)*dx_1(l1)
          dlnrhodx_yz=fac*(+ 45.0*(f(l1+1,:,:,ilnrho)-f(l1-1,:,:,ilnrho)) &
                           -  9.0*(f(l1+2,:,:,ilnrho)-f(l1-2,:,:,ilnrho)) &
                           +      (f(l1+3,:,:,ilnrho)-f(l1-3,:,:,ilnrho)))
           dsdx_yz=-(sigmaSBt*TT_yz**3+hcondxbot*(gamma_m1)*dlnrhodx_yz)/ &
              (chi_t*rho_yz+hcondxbot/cv)
!
!  enforce ds/dx = - (sigmaSBt*T^3 + hcond*(gamma-1)*glnrho)/(chi_t*rho+hcond/cv)
!
          do i=1,nghost
            f(l1-1,:,:,iss)=f(l1+i,:,:,iss)-2*i*dx*dsdx_yz
          enddo
        endif
!
!  top boundary
!  ============
!
      case ('top')
!
! For the case of pretend_lnTT=T, set glnTT=-sigma*T^3/hcond
!
        if (pretend_lnTT) then
            f(l2-1,:,:,iss)=f(l2+i,:,:,iss) - &
                2*i*dx*sigmaSBt*exp(f(l2,:,:,iss))**3/hcondxtop
        else
!
!  set ghost zones such that dsdx_yz obeys
!  - chi_t rho T dsdx_yz - hcond gTT = sigmaSBt*TT^4
!
          cs2_yz=cs20*exp(gamma_m1*(f(l2,:,:,ilnrho)-lnrho0)+cv1*f(l2,:,:,iss))
          TT_yz=cs2_yz/(gamma_m1*cp)
          rho_yz=exp(f(l2,:,:,ilnrho))
          fac=(1./60)*dx_1(l2)
          dlnrhodx_yz=fac*(+ 45.0*(f(l2+1,:,:,ilnrho)-f(l2-1,:,:,ilnrho)) &
                           -  9.0*(f(l2+2,:,:,ilnrho)-f(l2-2,:,:,ilnrho)) &
                           +      (f(l2+3,:,:,ilnrho)-f(l2-3,:,:,ilnrho)))
          dsdx_yz=-(sigmaSBt*TT_yz**3+hcondxtop*(gamma_m1)*dlnrhodx_yz)/ &
              (chi_t*rho_yz+hcondxtop/cv)
!
!  enforce ds/dx = - (sigmaSBt*T^3 + hcond*(gamma-1)*glnrho)/(chi_t*rho+hcond/cv)
!
          do i=1,nghost
            f(l2+i,:,:,iss)=f(l2-i,:,:,iss)+2*i*dx*dsdx_yz
          enddo
        endif
!
!  capture undefined entries
!
      case default
        call fatal_error('bc_ss_flux_turb_x','invalid argument')
      endselect
!
    endsubroutine bc_ss_flux_turb_x
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  23-jun-2003/tony: implemented for leos_fixed_ionization
!  26-aug-2003/tony: distributed across ionization modules
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp_old: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for entropy.
!  This assumes that the density is already set (ie density must register
!  first!)
!  tmp_xy = s(x,y) on the boundary.
!  gamma*s/cp = [ln(cs2/cs20)-(gamma-1)ln(rho/rho0)]
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if ((bcz1(ilnrho) /= 'a2') .and. (bcz1(ilnrho) /= 'a3')) &
          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3.')
        if (ldebug) print*, &
                'bc_ss_temp_old: set bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) &
              print*,'bc_ss_temp_old: cannot have cs2bot = ', cs2bot, ' <= 0'
        tmp_xy = (-gamma_m1*(f(:,:,n1,ilnrho)-lnrho0) &
             + log(cs2bot/cs20)) / gamma
        f(:,:,n1,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n1-i,iss) = 2*tmp_xy - f(:,:,n1+i,iss)
        enddo
!
!  top boundary
!
      case ('top')
        if ((bcz1(ilnrho) /= 'a2') .and. (bcz1(ilnrho) /= 'a3')) &
          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3.')
        if (ldebug) print*, &
                   'bc_ss_temp_old: set top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                   'bc_ss_temp_old: cannot have cs2top = ',cs2top, ' <= 0'
  !     if (bcz1(ilnrho) /= 'a2') &
  !          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 4.')
        tmp_xy = (-gamma_m1*(f(:,:,n2,ilnrho)-lnrho0) &
                 + log(cs2top/cs20)) / gamma
        f(:,:,n2,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n2+i,iss) = 2*tmp_xy - f(:,:,n2-i,iss)
        enddo
      case default
        call fatal_error('bc_ss_temp_old','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp_old
!***********************************************************************
    subroutine bc_ss_temp_x(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp_x: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_x: set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_x: cannot have cs2bot<=0'
        if (lentropy .and. .not. pretend_lnTT) then
           tmp = 2*cv*log(cs2bot/cs20)
!  Corrected for linear density
           if (ldensity_nolog) then
              f(l1,:,:,iss) = 0.5*tmp - (cp-cv)*(log(f(l1,:,:,ilnrho)) - lnrho0)
              do i=1,nghost
                 f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
                      - (cp-cv)*(log(f(l1+i,:,:,ilnrho)*f(l1-i,:,:,ilnrho)) - 2*lnrho0)
              enddo
           else
              f(l1,:,:,iss) = 0.5*tmp - (cp-cv)*(f(l1,:,:,ilnrho)-lnrho0)
              do i=1,nghost
                 f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
                      - (cp-cv)*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
              enddo
           endif
!
        elseif (lentropy .and. pretend_lnTT) then
           f(l1,:,:,iss) = log(cs2bot/gamma_m1)
           do i=1,nghost; f(l1-i,:,:,iss)=2*f(l1,:,:,iss)-f(l1+i,:,:,iss); enddo
        elseif (ltemperature) then
           f(l1,:,:,ilnTT) = log(cs2bot/gamma_m1)
           do i=1,nghost; f(l1-i,:,:,ilnTT)=2*f(l1,:,:,ilnTT)-f(l1+i,:,:,ilnTT); enddo
        endif
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                       'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                       'bc_ss_temp_x: cannot have cs2top<=0'
        if (lentropy .and. .not. pretend_lnTT) then
!
!  Distinguish cases for linear and logarithmic density
!
          tmp = 2*cv*log(cs2top/cs20)
          if (ldensity_nolog) then
            f(l2,:,:,iss) = 0.5*tmp - (cp-cv)*(log(f(l2,:,:,ilnrho))-lnrho0)
            do i=1,nghost
              f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
                - (cp-cv)*(log(f(l2-i,:,:,ilnrho)*f(l2+i,:,:,ilnrho))-2*lnrho0)
            enddo
          else
            f(l2,:,:,iss) = 0.5*tmp - (cp-cv)*(f(l2,:,:,ilnrho)-lnrho0)
            do i=1,nghost
              f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
                - (cp-cv)*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
            enddo
          endif
        elseif (lentropy .and. pretend_lnTT) then
          f(l2,:,:,iss) = log(cs2top/gamma_m1)
          do i=1,nghost; f(l2+i,:,:,iss)=2*f(l2,:,:,iss)-f(l2-i,:,:,iss); enddo
        elseif (ltemperature) then
          f(l2,:,:,ilnTT) = log(cs2top/gamma_m1)
          do i=1,nghost; f(l2+i,:,:,ilnTT)=2*f(l2,:,:,ilnTT)-f(l2-i,:,:,ilnTT); enddo
        endif
!
      case default
        call fatal_error('bc_ss_temp_x','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_y(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp_y: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_y: set y bottom temperature - cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_y: cannot have cs2bot<=0'
        tmp = 2*cv*log(cs2bot/cs20)
        f(:,m1,:,iss) = 0.5*tmp - (cp-cv)*(f(:,m1,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m1-i,:,iss) = -f(:,m1+i,:,iss) + tmp &
               - (cp-cv)*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                     'bc_ss_temp_y: set y top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                     'bc_ss_temp_y: cannot have cs2top<=0'
        tmp = 2*cv*log(cs2top/cs20)
        f(:,m2,:,iss) = 0.5*tmp - (cp-cv)*(f(:,m2,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
               - (cp-cv)*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
        enddo
!
      case default
        call fatal_error('bc_ss_temp_y','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp_z: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_z: cannot have cs2bot = ', cs2bot, ' <= 0'
        if (lentropy .and. .not. pretend_lnTT) then
!
!  Distinguish cases for linear and logarithmic density
!
           tmp = 2*cv*log(cs2bot/cs20)
           if (ldensity_nolog) then
             f(:,:,n1,iss) = 0.5*tmp - (cp-cv)*(alog(f(:,:,n1,irho))-lnrho0)
             do i=1,nghost
               f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
               !   - (cp-cv)*(log(f(:,:,n1+i,irho)*f(:,:,n1-i,irho))-2*lnrho0)
!AB: this could be better
                  - 2*(cp-cv)*(log(f(:,:,n1,irho))-lnrho0)
             enddo
           else
             f(:,:,n1,iss) = 0.5*tmp - (cp-cv)*(f(:,:,n1,ilnrho)-lnrho0)
             do i=1,nghost
               f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
                   - (cp-cv)*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
             enddo
           endif
        elseif (lentropy .and. pretend_lnTT) then
            f(:,:,n1,iss) = log(cs2bot/gamma_m1)
            do i=1,nghost; f(:,:,n1-i,iss)=2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
        elseif (ltemperature) then
            if (ltemperature_nolog) then
              f(:,:,n1,iTT)   = cs2bot/gamma_m1
            else
              f(:,:,n1,ilnTT) = log(cs2bot/gamma_m1)
            endif
            do i=1,nghost; f(:,:,n1-i,ilnTT)=2*f(:,:,n1,ilnTT)-f(:,:,n1+i,ilnTT); enddo
        endif
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                   'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                   'bc_ss_temp_z: cannot have cs2top = ', cs2top, ' <= 0'
!DM+PC next two lines need to be looked into.
        if (lread_oldsnap) &
          cs2top=cs20*exp(gamma*f(l2,m2,n2,iss)/cp+gamma_m1*f(l2,m2,n2,ilnrho))
        if (lentropy .and. .not. pretend_lnTT) then
!
!  Distinguish cases for linear and logarithmic density
!
          tmp = 2*cv*log(cs2top/cs20)
          if (ldensity_nolog) then
            f(:,:,n2,iss) = 0.5*tmp - (cp-cv)*(alog(f(:,:,n2,irho))-lnrho0)
            do i=1,nghost
              f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
                   !- (cp-cv)*(log(f(:,:,n2-i,irho)*f(:,:,n2+i,irho))-2*lnrho0)
!AB: this could be better
                   - 2*(cp-cv)*(log(f(:,:,n2,irho))-lnrho0)
            enddo
          else
            f(:,:,n2,iss) = 0.5*tmp - (cp-cv)*(f(:,:,n2,ilnrho)-lnrho0)
            do i=1,nghost
              f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
                   - (cp-cv)*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
            enddo
          endif
        elseif (lentropy .and. pretend_lnTT) then
            f(:,:,n2,iss) = log(cs2top/gamma_m1)
            do i=1,nghost; f(:,:,n2+i,iss)=2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
        elseif (ltemperature) then
            if (ltemperature_nolog) then
              f(:,:,n2,iTT)   = cs2top/gamma_m1
            else
              f(:,:,n2,ilnTT) = log(cs2top/gamma_m1)
            endif
            do i=1,nghost; f(:,:,n2+i,ilnTT)=2*f(:,:,n2,ilnTT)-f(:,:,n2-i,ilnTT); enddo
        endif
!
      case default
        call fatal_error('bc_ss_temp_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_lnrho_temp_z(f,topbot)
!
!  boundary condition for lnrho *and* ss: constant temperature
!
!  27-sep-2002/axel: coded
!  19-aug-2005/tobi: distributed across ionization modules
!
      use Gravity, only: gravz
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'bc_lnrho_temp_z: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*, &
                 'bc_lnrho_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*, &
                 'bc_lnrho_temp_z: cannot have cs2bot<=0'
        tmp = 2*cv*log(cs2bot/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n1,iss) = 0.5*tmp - (cp-cv)*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n1-i,iss) = 2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2bot
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=-gravz/cs2bot
        do i=1,nghost
          f(:,:,n1-i,ilnrho)=f(:,:,n1+i,ilnrho)+cp1*f(:,:,n1+i,iss) &
                                               -cp1*f(:,:,n1-i,iss)+2*i*dz*tmp
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                    'bc_lnrho_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*, &
                    'bc_lnrho_temp_z: cannot have cs2top<=0'
        tmp = 2*cv*log(cs2top/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n2,iss) = 0.5*tmp - (cp-cv)*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n2+i,iss) = 2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2top
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=gravz/cs2top
        do i=1,nghost
          f(:,:,n2+i,ilnrho)=f(:,:,n2-i,ilnrho)+cp1*f(:,:,n2-i,iss) &
                                               -cp1*f(:,:,n2+i,iss)+2*i*dz*tmp
        enddo
!
      case default
        call fatal_error('bc_lnrho_temp_z','invalid argument')
      endselect
!
    endsubroutine bc_lnrho_temp_z
!***********************************************************************
    subroutine bc_lnrho_pressure_z(f,topbot)
!
!  boundary condition for lnrho: constant pressure
!
!   4-apr-2003/axel: coded
!   1-may-2003/axel: added the same for top boundary
!  19-aug-2005/tobi: distributed across ionization modules
!
      use Gravity, only: lnrho_bot,lnrho_top,ss_bot,ss_top
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (ldebug) print*,'bc_lnrho_pressure_z: cs20,cs0=',cs20,cs0
!
!  Constant pressure, i.e. antisymmetric
!  This assumes that the entropy is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('top')
        if (ldebug) print*,'bc_lnrho_pressure_z: lnrho_top,ss_top=',lnrho_top,ss_top
!
!  fix entropy if inflow (uz>0); otherwise leave s unchanged
!  afterwards set s antisymmetrically about boundary value
!
        if (lentropy) then
!         do m=m1,m2
!         do l=l1,l2
!           if (f(l,m,n1,iuz)>=0) then
!             f(l,m,n1,iss)=ss_bot
!           else
!             f(l,m,n1,iss)=f(l,m,n1+1,iss)
!           endif
!         enddo
!         enddo
          f(:,:,n2,iss)=ss_top
          do i=1,nghost; f(:,:,n2+i,iss)=2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n2,ilnrho)=lnrho_top+cp1*(ss_top-f(:,:,n2,iss))
        else
          f(:,:,n2,ilnrho)=lnrho_top
        endif
!
!  make density antisymmetric about boundary
!  another possibility might be to enforce hydrostatics
!  ie to set dlnrho/dz=-g/cs^2, assuming zero entropy gradient
!
        do i=1,nghost
          f(:,:,n2+i,ilnrho)=2*f(:,:,n2,ilnrho)-f(:,:,n2-i,ilnrho)
        enddo
!
!  top boundary
!
      case ('bot')
        if (ldebug) print*,'bc_lnrho_pressure_z: lnrho_bot,ss_bot=',lnrho_bot,ss_bot
!
!  fix entropy if inflow (uz>0); otherwise leave s unchanged
!  afterwards set s antisymmetrically about boundary value
!
        if (lentropy) then
!         do m=m1,m2
!         do l=l1,l2
!           if (f(l,m,n1,iuz)>=0) then
!             f(l,m,n1,iss)=ss_bot
!           else
!             f(l,m,n1,iss)=f(l,m,n1+1,iss)
!           endif
!         enddo
!         enddo
          f(:,:,n1,iss)=ss_bot
          do i=1,nghost; f(:,:,n1-i,iss)=2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n1,ilnrho)=lnrho_bot+ss_bot-f(:,:,n1,iss)
        else
          f(:,:,n1,ilnrho)=lnrho_bot
        endif
!
!  make density antisymmetric about boundary
!  another possibility might be to enforce hydrostatics
!  ie to set dlnrho/dz=-g/cs^2, assuming zero entropy gradient
!
        do i=1,nghost
          f(:,:,n1-i,ilnrho)=2*f(:,:,n1,ilnrho)-f(:,:,n1+i,ilnrho)
        enddo
!
      case default
        call fatal_error('bc_lnrho_pressure_z','invalid argument')
      endselect
!
    endsubroutine bc_lnrho_pressure_z
!***********************************************************************
    subroutine bc_ss_temp2_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp2_z: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*, &
                   'bc_ss_temp2_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp2_z: cannot have cs2bot<=0'
        tmp = cv*log(cs2bot/cs20)
        do i=0,nghost
          f(:,:,n1-i,iss) = tmp &
               - (cp-cv)*(f(:,:,n1-i,ilnrho)-lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                     'bc_ss_temp2_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp2_z: cannot have cs2top<=0'
        tmp = cv*log(cs2top/cs20)
        do i=0,nghost
          f(:,:,n2+i,iss) = tmp &
               - (cp-cv)*(f(:,:,n2+i,ilnrho)-lnrho0)
        enddo
      case default
        call fatal_error('bc_ss_temp2_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_stemp_x: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (cs2bot<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2bot<=0'
        do i=1,nghost
          f(l1-i,:,:,iss) = f(l1+i,:,:,iss) &
               + (cp-cv)*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
        enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2top<=0'
        do i=1,nghost
          f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
               + (cp-cv)*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
        enddo
!
      case default
        call fatal_error('bc_ss_stemp_x','invalid argument')
      endselect
!
    endsubroutine bc_ss_stemp_x
!***********************************************************************
    subroutine bc_ss_stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_stemp_y: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (cs2bot<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2bot<=0'
        do i=1,nghost
          f(:,m1-i,:,iss) = f(:,m1+i,:,iss) &
               + (cp-cv)*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
        enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2top<=0'
        do i=1,nghost
          f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
               + (cp-cv)*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
        enddo
!
      case default
        call fatal_error('bc_ss_stemp_y','invalid argument')
      endselect
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_stemp_z: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
          if (cs2bot<=0.) print*, &
                                  'bc_ss_stemp_z: cannot have cs2bot<=0'
          do i=1,nghost
             f(:,:,n1-i,iss) = f(:,:,n1+i,iss) &
                  + (cp-cv)*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
          enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0.) print*, &
                 'bc_ss_stemp_z: cannot have cs2top<=0'
         do i=1,nghost
           f(:,:,n2+i,iss) = f(:,:,n2-i,iss) &
                + (cp-cv)*(f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho))
         enddo
      case default
        call fatal_error('bc_ss_stemp_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_a2stemp_x(f,topbot)
!
!  Boundary condition for entropy: adopt boundary value for temperature in
!  the ghost zone to handle shock profiles in interstellar with steep +ve
!  1st derivative in cooled remnant shells, followed by steep -ve 1st
!  derivative inside remnant.
!  s or a2 for temperature both unstable and unphysical as the unshocked
!  exterior ISM will be comparatively homogeneous, hence allowing the ghost
!  zone to fluctuate matching the boundary values is a reasonable approx
!  of the physical flow, whilst avoiding unphysical spikes to wreck the
!  calculation.
!
!  25-2010/fred: adapted from bc_ss_stemp_z
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_a2stemp_z: cs20,cs0=',cs20,cs0
!
!  Uniform temperature/sound speed condition for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
        case ('bot')
          if (cs2bot<=0.) print*, &
              'bc_ss_a2stemp_x: cannot have cs2bot<=0'
          do i=1,nghost
            f(l1-i,:,:,iss) = f(l1+1-i,:,:,iss)+(cp-cv)* &
                (f(l1+1-i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
          enddo
!
!  top boundary
!
        case ('top')
          if (cs2top<=0.) print*, &
              'bc_ss_a2stemp_x: cannot have cs2top<=0'
          do i=1,nghost
            f(l2+i,:,:,iss) = f(l2-1+i,:,:,iss)+(cp-cv)* &
                (f(l2-1+i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
          enddo
!
        case default
          call fatal_error('bc_ss_a2stemp_x','invalid argument')
      endselect
!
    endsubroutine bc_ss_a2stemp_x
!***********************************************************************
    subroutine bc_ss_a2stemp_y(f,topbot)
!
!  Boundary condition for entropy: adopt boundary value for temperature in
!  the ghost zone to handle shock profiles in interstellar with steep +ve
!  1st derivative in cooled remnant shells, followed by steep -ve 1st
!  derivative inside remnant.
!  s or a2 for temperature both unstable and unphysical as the unshocked
!  exterior ISM will be comparatively homogeneous, hence allowing the ghost
!  zone to fluctuate matching the boundary values is a reasonable approx
!  of the physical flow, whilst avoiding unphysical spikes to wreck the
!  calculation.
!
!  25-2010/fred: adapted from bc_ss_stemp_z
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_a2stemp_z: cs20,cs0=',cs20,cs0
!
!  Uniform temperature/sound speed condition for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
        case ('bot')
          if (cs2bot<=0.) print*, &
              'bc_ss_a2stemp_y: cannot have cs2bot<=0'
          do i=1,nghost
            f(:,m1-i,:,iss) = f(:,m1+1-i,:,iss)+(cp-cv)* &
                (f(:,m1+1-i,:,ilnrho)-f(:,m1-i,:,ilnrho))
          enddo
!
!  top boundary
!
        case ('top')
          if (cs2top<=0.) print*, &
              'bc_ss_a2stemp_y: cannot have cs2top<=0'
          do i=1,nghost
            f(:,m2+i,:,iss) = f(:,m2-1+i,:,iss)+(cp-cv)* &
                (f(:,m2-1+i,:,ilnrho)-f(:,m2+i,:,ilnrho))
          enddo
!
        case default
          call fatal_error('bc_ss_a2stemp_y','invalid argument')
      endselect
!
    endsubroutine bc_ss_a2stemp_y
!***********************************************************************
    subroutine bc_ss_a2stemp_z(f,topbot)
!
!  Boundary condition for entropy: adopt boundary value for temperature in
!  the ghost zone to handle shock profiles in interstellar with steep +ve
!  1st derivative in cooled remnant shells, followed by steep -ve 1st
!  derivative inside remnant.
!  s or a2 for temperature both unstable and unphysical as the unshocked
!  exterior ISM will be comparatively homogeneous, hence allowing the ghost
!  zone to fluctuate matching the boundary values is a reasonable approx
!  of the physical flow, whilst avoiding unphysical spikes to wreck the
!  calculation.
!
!  25-2010/fred: adapted from bc_ss_stemp_z
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_a2stemp_z: cs20,cs0=',cs20,cs0
!
!  Uniform temperature/sound speed condition for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
        case ('bot')
          if (cs2bot<=0.) print*, &
              'bc_ss_a2stemp_z: cannot have cs2bot<=0'
          do i=1,nghost
            f(:,:,n1-i,iss) = f(:,:,n1+1-i,iss) + (cp-cv)* &
                (f(:,:,n1+1-i,ilnrho)-f(:,:,n1-i,ilnrho))
          enddo
!
!  top boundary
!
        case ('top')
          if (cs2top<=0.) print*, &
              'bc_ss_a2stemp_z: cannot have cs2top<=0'
          do i=1,nghost
            f(:,:,n2+i,iss) = f(:,:,n2-1+i,iss) + (cp-cv)* &
                (f(:,:,n2-1+i,ilnrho)-f(:,:,n2+i,ilnrho))
          enddo
        case default
          call fatal_error('bc_ss_a2stemp_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_a2stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
!  boundary condition for entropy
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the entropy module
!  26-aug-2003/tony: distributed across ionization modules
!
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: cs2_2d
      integer :: i
!
!  The 'ce' boundary condition for entropy makes the energy constant at
!  the boundaries.
!  This assumes that the density is already set (ie density must register
!  first!)
!
    select case (topbot)
!
! Bottom boundary
!
    case ('bot')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma_m1*f(:,:,n1,ilnrho)+cv1*f(:,:,n1,iss))
      do i=1,nghost
         f(:,:,n1-i,iss)=cv*(-gamma_m1*f(:,:,n1-i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
!
! Top boundary
!
    case ('top')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma_m1*f(:,:,n2,ilnrho)+cv1*f(:,:,n2,iss))
      do i=1,nghost
         f(:,:,n2+i,iss)=cv*(-gamma_m1*f(:,:,n2+i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
    case default
      call fatal_error('bc_ss_energy','invalid argument')
    endselect
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_stellar_surface: NOT IMPLEMENTED IN EOS_IDEALGAS")
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)
!
!  Boundary condition for radial centrifugal balance
!
!  This sets
!    \partial_{r} \ln\rho
!  such that
!    \partial_{r} p = uphi**2/rad - \partial_{r} Phi
!  where Phi is the gravitational potential
!
!  i.e. it enforces centrifugal balance at the boundary.
!
!  As it is, works only for isobaric, isothermal and cylindrical coordinates
!
!  21-aug-2006/wlad: coded
!
      use Gravity, only: potential
      use Sub, only: div
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
      real, dimension (my,mz) :: cs2,gravterm,centterm,uphi
      real :: potp,potm,rad,step
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary
!
      case ('bot')
        do i=1,nghost
!
          cs2 = cs20
          call potential(R=x(l1-i),pot=potm)
          call potential(R=x(l1+i),pot=potp)
!
          gravterm= -(potm-potp)/cs2
!
          step=-2*i*dx
          rad=x(l1-i)
          uphi=f(l1-i,:,:,iuy)
!
          centterm= uphi**2 * step/(rad*cs2)
          if (ldensity_nolog) then
            f(l1-i,:,:,ilnrho)=f(l1+i,:,:,irho)*exp(gravterm + centterm)
          else
            f(l1-i,:,:,ilnrho)=f(l1+i,:,:,ilnrho) + gravterm + centterm
          endif
!
          !print*,'potentials',potm,potp,-(potm-potp)
          !print*,'centrifugal',f(l1-i,mpoint,npoint,iuy)**2 *step/rad
          !stop
!
        enddo
!
!  Top boundary
!
      case ('top')
        do i=1,nghost
!
          cs2 = cs20
          call potential(R=x(l2+i),pot=potp)
          call potential(R=x(l2-i),pot=potm)
!
          gravterm= -(potp-potm)/cs2
!
          step=2*i*dx
          rad=x(l2+i)
          uphi=f(l2+i,:,:,iuy)
!
          centterm= uphi**2 * step/(rad*cs2)
          if (ldensity_nolog) then
            f(l2+i,:,:,irho)   = f(l2-i,:,:,irho)*exp(gravterm + centterm)
          else
            f(l2+i,:,:,ilnrho) = f(l2-i,:,:,ilnrho) + gravterm + centterm
          endif
!
          !if (i==nghost) then
          !  print*,'potentials',potp,potm,-potp+potm,-(potp-potm)
          !  print*,'centrifugal',f(l2+i,mpoint,npoint,iuy)**2 *step/rad
          !  stop
          !endif
        enddo
!
      case default
!
      endselect
!
    endsubroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso(f,topbot)
!
!  Boundary condition for density *and* entropy.
!
!  This sets
!    \partial_{z} \ln\rho
!  such that
!    \partial_{z} p = \rho g_{z},
!  i.e. it enforces hydrostatic equlibrium at the boundary.
!
!  Currently this is only correct if
!    \partial_{z} lnT = 0
!  at the boundary.
!
!  12-Juil-2006/dintrans: coded
!
      use Gravity, only: potential, gravz
      use Sub, only: div
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
!
      real, dimension (mx,my) :: cs2
      real, dimension (nx) :: shock,divu
      real :: dlnrhodz, dssdz, cs2_point
      real :: potp,potm
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary
!
      case ('bot')
!
        if (lentropy) then
!
!  The following might work for anelastic
!
          if (ldensity) then
            if (bcz1(iss)/='hs') then
              call fatal_error("bc_lnrho_hydrostatic_z", &
                "This boundary condition for density is "// &
                "currently only correct for bcz1(iss)='hs'")
            endif
!
            call eoscalc(ilnrho_ss,f(l1,m1,n1,ilnrho),f(l1,m1,n1,iss), &
                cs2=cs2_point)
!
            dlnrhodz =  gamma *gravz/cs2_point
            dssdz    = -gamma_m1*gravz/cs2_point
!
            do i=1,nghost
              f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) - 2*i*dz*dlnrhodz
              f(:,:,n1-i,iss   ) = f(:,:,n1+i,iss   ) - 2*i*dz*dssdz
            enddo
          else if (lanelastic) then
            if (bcz1(iss_b)/='hs') then
              call fatal_error("bc_lnrho_hydrostatic_z", &
                "This boundary condition for density is "// &
                "currently only correct for bcz1(iss)='hs'")
            endif
            call eoscalc(ipp_ss,log(f(l1,m1,n1,irho_b)),f(l1,m1,n1,iss_b), &
                cs2=cs2_point)
!
            dlnrhodz =  gamma *gravz/cs2_point
            dssdz    = gamma_m1*gravz/cs2_point
!
            do i=1,nghost
              f(:,:,n1-i,irho_b) = f(:,:,n1+i,irho_b) - 2*i*dz*dlnrhodz*f(:,:,n1+1,irho_b)
              f(:,:,n1-i,iss_b   ) = f(:,:,n1+i,iss_b   ) - 2*i*dz*dssdz
            enddo
          endif
!
        elseif (ltemperature) then
!
!  Energy equation formulated in logarithmic temperature.
!
          if (bcz1(ilntt)/='s') then
            call fatal_error("bc_lnrho_hydrostatic_z", &
                "This boundary condition for density is "// &
                "currently only correct for bcz1(ilntt)='s'")
          endif
!
          call eoscalc(ilnrho_lntt,f(l1,m1,n1,ilnrho),f(l1,m1,n1,ilntt), &
              cs2=cs2_point)
!
          dlnrhodz =  gamma *gravz/cs2_point
!
          do i=1,nghost
            f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) - 2*i*dz*dlnrhodz
          enddo
!
        else
!
!  Isothermal equation of state.
!
          do i=1,nghost
            call potential(z=z(n1-i),pot=potm)
            call potential(z=z(n1+i),pot=potp)
            cs2 = cs2bot
!
            if (.false.) then
              ! Note: Since boundconds_x and boundconds_y are called first,
              ! this doesn't set the corners properly. However, this is
              ! not a problem since cross derivatives of density are never
              ! needed.
              n = n1+i
              do m = m1,m2
                shock = f(l1:l2,m,n,ishock)
                call div(f,iuu,divu)
                cs2(l1:l2,m) = cs2bot - shock*divu
              enddo
            endif
!
            if (ldensity_nolog) then
              f(:,:,n1-i,irho)   = f(:,:,n1+i,irho)*exp(-(potm-potp)/cs2)
            else
              f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) - (potm-potp)/cs2
            endif
!
          enddo
!
        endif
!
!  Top boundary
!
      case ('top')
!
        if (lentropy) then
!
          if (bcz2(iss)/='hs') then
            call fatal_error("bc_lnrho_hydrostatic_z", &
                "This boundary condition for density is "//&
                "currently only correct for bcz2(iss)='hs'")
          endif
!
          call eoscalc(ilnrho_ss,f(l2,m2,n2,ilnrho),f(l2,m2,n2,iss), &
              cs2=cs2_point)
!
          dlnrhodz =  gamma *gravz/cs2_point
          dssdz    = -gamma_m1*gravz/cs2_point
!
          do i=1,nghost
            f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) + 2*i*dz*dlnrhodz
            f(:,:,n2+i,iss   ) = f(:,:,n2-i,iss   ) + 2*i*dz*dssdz
          enddo
!
        elseif (ltemperature) then
!
!  Energy equation formulated in logarithmic temperature.
!
          if (bcz2(ilntt)/='s') then
            call fatal_error("bc_lnrho_hydrostatic_z", &
                "This boundary condition for density is "//&
                "currently only correct for bcz2(ilntt)='s'")
          endif
!
          call eoscalc(ilnrho_lntt,f(l2,m2,n2,ilnrho),f(l2,m2,n2,ilntt), &
              cs2=cs2_point)
!
          dlnrhodz =  gamma *gravz/cs2_point
!
          do i=1,nghost
            f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) + 2*i*dz*dlnrhodz
          enddo
!
        else
!
!  Isothermal equation of state.
!
          do i=1,nghost
            call potential(z=z(n2+i),pot=potp)
            call potential(z=z(n2-i),pot=potm)
            cs2 = cs2bot
            if (.false.) then
              ! Note: Since boundconds_x and boundconds_y are called first,
              ! this doesn't set the corners properly. However, this is
              ! not a problem since cross derivatives of density are never
              ! needed.
              n = n2-i
              do m = m1,m2
                shock = f(l1:l2,m,n,ishock)
                call div(f,iuu,divu)
                cs2(l1:l2,m) = cs2top - shock*divu
              enddo
            else
            endif
            if (ldensity_nolog) then
              f(:,:,n2+i,irho)   = f(:,:,n2-i,irho)*exp(-(potp-potm)/cs2)
            else
              f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) - (potp-potm)/cs2
            endif
          enddo
!
        endif
!
      case default
!
      endselect
!
    endsubroutine bc_lnrho_hds_z_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)
!
!  Smooth out density perturbations with respect to hydrostatic
!  stratification in Fourier space.
!
!  Note: Since boundconds_x and boundconds_y are called first,
!  this doesn't set the corners properly. However, this is
!  not a problem since cross derivatives of density are never
!  needed.
!
!  05-jul-07/tobi: Adapted from bc_aa_pot3
!
      use Fourier, only: fourier_transform_xy_xy, fourier_transform_other
      use Gravity, only: potential
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
!
      real, dimension (nx,ny) :: kx,ky,kappa,exp_fact
      real, dimension (nx,ny) :: tmp_re,tmp_im
      real :: pot
      integer :: i
!
!  Get local wave numbers
!
      kx = spread(kx_fft(ipx*nx+1:ipx*nx+nx),2,ny)
      ky = spread(ky_fft(ipy*ny+1:ipy*ny+ny),1,nx)
!
!  Calculate 1/k^2, zero mean
!
      if (lshear) then
        kappa = sqrt((kx+ky*deltay/Lx)**2+ky**2)
      else
        kappa = sqrt(kx**2 + ky**2)
      endif
!
!  Check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  Potential field condition at the bottom
!
      case ('bot')
!
        do i=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          exp_fact = exp(-kappa*(z(n1+i)-z(n1-i)))
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
          call potential(z=z(n1+i),pot=pot)
          if (ldensity_nolog) then
            tmp_re = f(l1:l2,m1:m2,n1+i,irho)*exp(+pot/cs2bot)
          else
            tmp_re = f(l1:l2,m1:m2,n1+i,ilnrho) + pot/cs2bot
          endif
          tmp_im = 0.0
          if (nxgrid>1 .and. nygrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im)
          else
            call fourier_transform_other(tmp_re,tmp_im)
          endif
          tmp_re = tmp_re*exp_fact
          tmp_im = tmp_im*exp_fact
          ! Transform back
          if (nxgrid>1 .and. nygrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
          else
            call fourier_transform_other(tmp_re,tmp_im,linv=.true.)
          endif
          call potential(z=z(n1-i),pot=pot)
          if (ldensity_nolog) then
            f(l1:l2,m1:m2,n1-i,irho)   = tmp_re*exp(-pot/cs2bot)
          else
            f(l1:l2,m1:m2,n1-i,ilnrho) = tmp_re - pot/cs2bot
          endif
!
        enddo
!
!  Potential field condition at the top
!
      case ('top')
!
        do i=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          exp_fact = exp(-kappa*(z(n2+i)-z(n2-i)))
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
          call potential(z=z(n2-i),pot=pot)
          if (ldensity_nolog) then
            tmp_re = f(l1:l2,m1:m2,n2-i,irho)*exp(+pot/cs2top)
          else
            tmp_re = f(l1:l2,m1:m2,n2-i,ilnrho) + pot/cs2top
          endif
          tmp_im = 0.0
          if (nxgrid>1 .and. nygrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im)
          else
            call fourier_transform_other(tmp_re,tmp_im)
          endif
          tmp_re = tmp_re*exp_fact
          tmp_im = tmp_im*exp_fact
          ! Transform back
          if (nxgrid>1 .and. nygrid>1) then
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
          else
            call fourier_transform_other(tmp_re,tmp_im,linv=.true.)
          endif
          call potential(z=z(n2+i),pot=pot)
          if (ldensity_nolog) then
            f(l1:l2,m1:m2,n2+i,irho)   = tmp_re*exp(-pot/cs2top)
          else
            f(l1:l2,m1:m2,n2+i,ilnrho) = tmp_re - pot/cs2top
          endif
!
        enddo
!
      case default
!
        if (lroot) print*,"bc_lnrho_hydrostatic_z_smooth: invalid argument"
!
      endselect
!
    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine read_transport_data
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine read_transport_data
!***********************************************************************
    subroutine write_thermodyn()
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine write_thermodyn
!***********************************************************************
    subroutine read_thermodyn(input_file)
!
      character (len=*), intent(in) :: input_file
!
      call keep_compiler_quiet(input_file)
!
    endsubroutine read_thermodyn
!***********************************************************************
    subroutine read_species(input_file)
!
      character (len=*) :: input_file
!
      call keep_compiler_quiet(input_file)
!
    endsubroutine read_species
!***********************************************************************
    subroutine find_species_index(species_name,ind_glob,ind_chem,found_specie)
!
      integer, intent(out) :: ind_glob
      integer, intent(inout) :: ind_chem
      character (len=*), intent(in) :: species_name
      logical, intent(out) :: found_specie
!
       call keep_compiler_quiet(ind_glob)
       call keep_compiler_quiet(ind_chem)
       call keep_compiler_quiet(species_name)
       call keep_compiler_quiet(found_specie)
!
     endsubroutine find_species_index
!***********************************************************************
     subroutine find_mass(element_name,MolMass)
!
       character (len=*), intent(in) :: element_name
       real, intent(out) :: MolMass
!
       call keep_compiler_quiet(element_name)
       call keep_compiler_quiet(MolMass)
!
     endsubroutine find_mass
!***********************************************************************
    subroutine read_Lewis
!
!  Dummy routine
!
    endsubroutine read_Lewis
!***********************************************************************
endmodule EquationOfState
