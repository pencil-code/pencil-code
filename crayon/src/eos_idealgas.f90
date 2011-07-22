! $Id: eos_idealgas.f90 13598 2010-04-06 10:26:10Z tavo.buk $
!
!  Equation of state for an ideal gas without ionization.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cp; cp1; cp1tilde
! PENCILS PROVIDED glnTT(3); TT; TT1; gTT(3); yH; hss(3,3); hlnTT(3,3)
! PENCILS PROVIDED del2ss; del2lnTT; cv1; gamma
! PENCILS PROVIDED del2TT; glnmumol(3); ppvap; csvap2
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
  integer, parameter :: irho_TT=10, ipp_ss=11, ipp_cs2=12
  real, dimension(mz) :: profz_eos=1.0,dprofz_eos=0.0
  real, dimension(3) :: beta_glnrho_global=0.0, beta_glnrho_scaled=0.0
  real :: lnTT0=impossible
  real :: xHe=0.0
  real :: mu=1.0
  real :: cs0=1.0, rho0=1.0, pp0=1.0
  real :: cs20=1.0, lnrho0=0.0
  real :: ptlaw=0.0
  real :: gamma=5.0/3.0
  real :: Rgas_cgs=0.0, Rgas, error_cp=1.0e-6
  real :: gamma_m1    !(=gamma-1)
  real :: gamma_inv   !(=1/gamma)
  real :: cp=impossible, cp1=impossible, cv=impossible, cv1=impossible
  real :: pres_corr=0.1
  real :: cs2top_ini=impossible, dcs2top_ini=impossible
  real :: cs2bot=1.0, cs2top=1.0
  real :: cs2cool=0.0
  real :: mpoly=1.5, mpoly0=1.5, mpoly1=1.5, mpoly2=1.5
  real :: width_eos_prof=0.2
  integer :: isothtop=0
  integer :: ieosvars=-1, ieosvar1=-1, ieosvar2=-1, ieosvar_count=0
  logical :: leos_isothermal=.false., leos_isentropic=.false.
  logical :: leos_isochoric=.false., leos_isobaric=.false.
  character (len=labellen) :: ieos_profile='nothing'
!
!  Input parameters.
!
  namelist /eos_init_pars/ &
      xHe, mu, cp, cs0, rho0, gamma, error_cp, ptlaw, cs2top_ini, dcs2top_ini
!
!  Run parameters.
!
  namelist /eos_run_pars/ &
      xHe, mu, cp, cs0, rho0, gamma, error_cp, ptlaw, cs2top_ini, &
      dcs2top_ini, ieos_profile, width_eos_prof,pres_corr
!
  contains
!***********************************************************************
    subroutine register_eos()
!
!  Register variables from the EquationOfState module.
!
!  14-jun-03/axel: adapted from register_eos
!
      leos=.true.
      leos_idealgas=.true.
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_eos: ionization nvar = ', nvar
      endif
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          '$Id: eos_idealgas.f90 13598 2010-04-06 10:26:10Z tavo.buk $')
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
      use SharedVariables, only: put_shared_variable
      use Mpicomm, only: stop_it
      use Sub, only: erfunc
!
      real :: Rgas_unit_sys, cp_reference
      integer :: ierr
!
!  Set gamma_m1, cs20, and lnrho0.
!  (used currently for non-dimensional equation of state)
!
      gamma_m1=gamma-1.0
      gamma_inv=1/gamma
!
!  Avoid floating overflow if cs0 was not set.
!
      cs20=cs0**2
      lnrho0=log(rho0)
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
          Rgas=mu*(1.0-gamma_inv)*cp
        endif
        unit_temperature=unit_velocity**2*Rgas/Rgas_unit_sys
      else
        Rgas=Rgas_unit_sys*unit_temperature/unit_velocity**2
        if (cp==impossible) then
          if (gamma_m1==0.0) then
            cp=Rgas/mu
          else
            cp=Rgas/(mu*gamma_m1*gamma_inv)
          endif
        else
!
!  Checking whether the units are overdetermined.
!  This is assumed to be the case when the two differ by error_cp.
!
          if (gamma_m1==0.0) then
            cp_reference=Rgas/mu
          else
            cp_reference=Rgas/(mu*gamma_m1*gamma_inv)
          endif
          if (abs(cp-cp_reference)/cp > error_cp) then
            if (lroot) print*,'initialize_eos: consistency: cp=', cp , &
                'while: cp_reference=', cp_reference
            call fatal_error('initialize_eos','')
          endif
        endif
      endif
      cp1=1/cp
      cv=gamma_inv*cp
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
    endsubroutine initialize_eos
!***********************************************************************
    subroutine select_eos_variable(variable,findex)
!
!  Select eos variable.
!
!   02-apr-06/tony: implemented
!
      use FArrayManager
!
      character (len=*), intent(in) :: variable
      integer, intent(in) :: findex
      integer :: this_var
      integer, save :: ieosvar_selected=0
      integer, parameter :: ieosvar_lnrho = 2**0
      integer, parameter :: ieosvar_rho   = 2**1
      integer, parameter :: ieosvar_ss    = 2**2
      integer, parameter :: ieosvar_lnTT  = 2**3
      integer, parameter :: ieosvar_TT    = 2**4
      integer, parameter :: ieosvar_cs2   = 2**5
      integer, parameter :: ieosvar_pp    = 2**6
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
        if (findex<0) then
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
    endsubroutine getmu
!***********************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
!  Writes to index.pro file.
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
!  20-11-04/anders: coded
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
            lpencil_in(i_TT)=.true.
          endif
          if (lpencil_in(i_rho)) lpencil_in(i_lnrho)=.true.
        else
          if (lpencil_in(i_cs2)) then
            lpencil_in(i_rho)=.true.
            lpencil_in(i_pp)=.true.
          endif
          if (lpencil_in(i_TT1)) lpencil_in(i_TT)=.true.
          if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
          if (lpencil_in(i_lnTT)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_rho)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_lnrho)) then
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
      integer :: i
!
!  THE FOLLOWING 2 ARE CONCEPTUALLY WRONG
!  FOR pretend_lnTT since iss actually contain lnTT NOT entropy!
!  The code is not wrong however since this is correctly
!  handled by the eos module.

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
          if (lpencil(i_cs2)) p%cs2=cs20*exp(gamma_m1*(p%lnrho-lnrho0))
        elseif (leos_isothermal) then
          if (lpencil(i_ss)) p%ss=-(cp-cv)*(p%lnrho-lnrho0)
          if (lpencil(i_gss)) p%gss=-(cp-cv)*p%glnrho
          if (lpencil(i_hss)) p%hss=-(cp-cv)*p%hlnrho
          if (lpencil(i_del2ss)) p%del2ss=-(cp-cv)*p%del2lnrho
          if (lpencil(i_cs2)) p%cs2=cs20
        else
          if (lpencil(i_ss)) p%ss=f(l1:l2,m,n,ieosvar2)
          if (lpencil(i_gss)) call grad(f,ieosvar2,p%gss)
          if (lpencil(i_hss)) call g2ij(f,ieosvar2,p%hss)
          if (lpencil(i_del2ss)) call del2(f,ieosvar2,p%del2ss)
          if (lpencil(i_cs2)) p%cs2=cs20*exp(cv1*p%ss+gamma_m1*(p%lnrho-lnrho0))
        endif
        if (lpencil(i_lnTT)) p%lnTT=lnTT0+cv1*p%ss+gamma_m1*(p%lnrho-lnrho0)
        if (lpencil(i_pp)) p%pp=(cp-cv)*exp(p%lnTT+p%lnrho)
        if (lpencil(i_ee)) p%ee=cv*exp(p%lnTT)
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
        if (lpencil(i_glnTT)) p%glnTT=gamma_m1*p%glnrho+cv1*p%gss
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
        else
          if (lpencil(i_lnTT)) p%lnTT=f(l1:l2,m,n,ieosvar2)
          if (lpencil(i_glnTT)) call grad(f,ieosvar2,p%glnTT)
          if (lpencil(i_hlnTT)) call g2ij(f,ieosvar2,p%hlnTT)
          if (lpencil(i_del2lnTT)) call del2(f,ieosvar2,p%del2lnTT)
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
!
!  Work out thermodynamic quantities for given lnrho or rho and TT.
!
      case (ilnrho_TT,irho_TT)
        if (lpencil(i_TT))   p%TT=f(l1:l2,m,n,ieosvar2)
        if (lpencil(i_TT1))  p%TT1=1/f(l1:l2,m,n,ieosvar2)
        if (lpencil(i_lnTT)) p%lnTT=log(f(l1:l2,m,n,ieosvar2))
        if (lpencil(i_cs2))  p%cs2=cp*gamma_m1*f(l1:l2,m,n,ieosvar2)
        if (lpencil(i_gTT))  call grad(f,ieosvar2,p%gTT)
        if (lpencil(i_glnTT)) then
          do i=1,3; p%glnTT(:,i)=p%gTT(:,i)*p%TT1; enddo
        endif
        if (lpencil(i_del2TT)) call del2(f,ieosvar2,p%del2TT)
        if (lpencil(i_del2lnTT)) then
          tmp=0.0
          do i=1,3
            tmp=tmp+p%glnTT(:,i)**2
          enddo
          p%del2lnTT=p%del2TT*p%TT1-tmp
        endif
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
          if (lpencil(i_pp)) p%pp=gamma_inv*p%rho*cs20
          call fatal_error('calc_pencils_eos', &
              'Full equation of state not implemented for ilnrho_cs2')
        endif
!
!  Internal energy.
!  For gamma=1, we use R/mu = c_p = c_v, thus ee = c_vT = R/mu T = p/rho = cs^2.
!
        if (lpencil(i_ee)) then
          if (gamma_m1/=0.0) then
            p%ee=(gamma_inv/gamma_m1)*p%cs2
          else
            p%ee=p%cs2
          endif
        endif
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
      case default
        call fatal_error('calc_pencils_eos','case not implemented yet')
      endselect
!
!  Inverse cv and cp values.
!
      if (lpencil(i_cv1)) p%cv1=cv1
      if (lpencil(i_cp1)) p%cp1=cp1
      if (lpencil(i_cp))  p%cp=1/p%cp1
      if (lpencil(i_cp1tilde)) p%cp1tilde=cp1
!
      if (lpencil(i_glnmumol)) p%glnmumol(:,:)=0.
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

     real, dimension (mx,my,mz,mfarray) :: f
     real, dimension (mx,my,mz), intent(out) :: TT_tmp
!
     call keep_compiler_quiet(f)
     call keep_compiler_quiet(TT_tmp)
!  
   endsubroutine gettemperature
!***********************************************************************
 subroutine getpressure(pp_tmp)

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
      cs2=cs20*exp(cv1*ss+gamma_m1*(lnrho-lnrho0))
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
      cs2=cs20*exp(cv1*ss+gamma_m1*(lnrho-lnrho0))
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
      glnTT=gamma_m1*glnrho+cv1*gss
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
      if (ldensity_nolog) then
        call dot2(p%grho,tmp)
        p%del2lnTT=gamma_m1*p%rho1*(p%del2rho+p%rho1*tmp)
      else
        p%del2lnTT=gamma_m1*p%del2lnrho+p%cv1*p%del2ss
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
      if (gamma_m1==0.) call fatal_error('temperature_hessian','gamma=1 not allowed w/entrop')
!
      hlnTT=gamma_m1*hlnrho+cv1*hss
!
      call keep_compiler_quiet(f)
!
    endsubroutine temperature_hessian
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

      if (psize==nx) then
        lnrho_=f(l1:l2,m,n,ilnrho)
        if (present(ee)) then
          f(l1:l2,m,n,iss)=cv*(log(cv1*ee)-lnTT0-gamma_m1*(lnrho_-lnrho0))
        elseif (present(pp)) then
          f(l1:l2,m,n,iss)=cv*(log(gamma*pp/gamma_m1)-&
               gamma*lnrho_-gamma_m1*lnrho0-lnTT0)
        elseif (present(ss)) then
          f(l1:l2,m,n,iss)=ss
        endif
      elseif (psize==mx) then
        lnrho_=f(:,m,n,ilnrho)
        if (present(ee)) then
          f(:,m,n,iss)=cv*(log(cv1*ee)-lnTT0-gamma_m1*(lnrho_-lnrho0))
        elseif (present(pp)) then
          f(:,m,n,iss)=cv*(log(gamma*pp/gamma_m1)-&
               gamma*lnrho_-gamma_m1*lnrho0-lnTT0)
        elseif (present(ss)) then
          f(:,m,n,iss)=ss
        endif
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
            lnrho_=alog(f(l1:l2,m,n,ieosvar1))
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
            lnrho_=alog(f(:,m,n,ieosvar1))
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
          lnrho_=f(l1:l2,m,n,ieosvar1)
          if (leos_isentropic) then
            lnTT_=lnTT0+(cp-cv)*(lnrho_-lnrho0)
          elseif (leos_isothermal) then
            lnTT_=lnTT0
          else
            lnTT_=f(l1:l2,m,n,ieosvar2)
          endif
        case (mx)
          lnrho_=f(:,m,n,ieosvar1)
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
        if (ieosvars==ilnrho_lnTT) then
          if (present(pp)) pp=(cp-cv)*exp(lnTT_+lnrho_)
        else
          if (present(pp)) pp=(cp-cv)*exp(lnTT_)*lnrho_
        endif
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
            lnrho_=alog(f(l1:l2,m,n,ieosvar1))
          endif
          if (leos_isentropic) then
            cs2_=exp(gamma_m1*(lnrho_-lnrho0)+log(cs20))
          elseif (leos_isothermal) then
            cs2_=cs20
          else
            call fatal_error('eoscalc_farray','full eos for cs2 not implemented')
          endif
        case (mx)
          lnrho_=f(:,m,n,ieosvar1)
          if (leos_isentropic) then
            cs2_=exp(gamma_m1*(lnrho_-lnrho0)+log(cs20))
          elseif (leos_isothermal) then
            cs2_=cs20
          else
            call fatal_error('eoscalc_farray','full eos for cs2 not implemented')
          endif
        case default
          call fatal_error('eoscalc_farray','no such pencil size')
        end select
!
        if (present(lnrho)) lnrho=lnrho_
        if (present(lnTT)) lnTT=lnTT0+log(cs2_)
        if (present(ee)) ee=gamma_inv*cs2_/gamma_m1
        if (present(pp)) pp=gamma_inv*cs2_*exp(lnrho_)
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
!   27-mar-06/tony: Introduces cv, cv1, gamma_inv to make faster
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
      if (gamma_m1==0) call fatal_error &
        ('eoscalc_point','gamma=1 not allowed w/entropy')
!
      select case (ivars)

      case (ilnrho_ss)
        lnrho_=var1
        ss_=var2
        lnTT_=lnTT0+cv1*ss_+gamma_m1*(lnrho_-lnrho0)
        ee_=cv*exp(lnTT_)
        pp_=(cp-cv)*exp(lnTT_+lnrho_)
        cs2_=gamma*gamma_m1*ee_

      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
        lnTT_=log(cv1*ee_)
        ss_=cv*(lnTT_-lnTT0-gamma_m1*(lnrho_-lnrho0))
        pp_=gamma_m1*ee_*exp(lnrho_)
        cs2_=gamma*gamma_m1*ee_

      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
        ss_=cv*(log(pp_*exp(-lnrho_)*gamma/cs20)-gamma_m1*(lnrho_-lnrho0))
        ee_=pp_*exp(-lnrho_)/gamma_m1
        lnTT_=log(cv1*ee_)
        cs2_=gamma*gamma_m1*ee_

      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
        ss_=cv*(lnTT_-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*exp(lnTT_)
        pp_=ee_*exp(lnrho_)*gamma_m1
        cs2_=gamma*gamma_m1*ee_

      case (ilnrho_TT)
        lnrho_=var1
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*exp(lnrho_)*gamma_m1
        cs2_=cp*gamma_m1*TT_

      case (irho_TT)
        lnrho_=log(var1)
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*var1*gamma_m1
        cs2_=cp*gamma_m1*TT_

      case (ipp_cs2)
        if (leos_isothermal) then
        pp_=var1
        lnrho_=pp_*cs20
        TT_=exp(lnTT0)
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
!   27-mar-06/tony: Introduces cv, cv1, gamma_inv to make faster
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

      case (ilnrho_ss)
        lnrho_=var1
        ss_=var2
        lnTT_=lnTT0+cv1*ss_+gamma_m1*(lnrho_-lnrho0)
        ee_=cv*exp(lnTT_)
        pp_=(cp-cv)*exp(lnTT_+lnrho_)
        cs2_=gamma*gamma_m1*ee_
        cs2_=cs20*cv1*ee_

      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
        lnTT_=log(cv1*ee_)
        ss_=cv*(lnTT_-lnTT0-gamma_m1*(lnrho_-lnrho0))
        pp_=gamma_m1*ee_*exp(lnrho_)
        cs2_=gamma*gamma_m1*ee_

      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
        ss_=cv*(log(pp_*exp(-lnrho_)*gamma/cs20)-gamma_m1*(lnrho_-lnrho0))
        ee_=pp_*exp(-lnrho_)/gamma_m1
        lnTT_=log(cv1*ee_)
        cs2_=gamma*gamma_m1*ee_

      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
        ss_=cv*(lnTT_-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*exp(lnTT_)
        pp_=ee_*exp(lnrho_)*gamma_m1
        cs2_=gamma*gamma_m1*ee_

      case (ilnrho_TT)
        lnrho_=var1
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*exp(lnrho_)*gamma_m1
        cs2_=cp*gamma_m1*TT_

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
        lnrho_=log(pp)/gamma-ss/cp
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
    subroutine get_soundspeed(lnTT,cs2)
!
!  Calculate sound speed for given temperature
!
!  20-Oct-03/tobi: Coded
!
      real, intent(in)  :: lnTT
      real, intent(out) :: cs2
!
      cs2=gamma_m1*cp*exp(lnTT)
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
    subroutine get_average_pressure(init_average_density,average_density,&
                                    average_pressure)
!
!   01-dec-2009/piyali+dhruba: coded
!
      use Cdata
!
      real,intent(in):: init_average_density,average_density
      real,intent(inout):: average_pressure
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
endmodule EquationOfState
