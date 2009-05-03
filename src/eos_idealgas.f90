! $Id$

!  Equation of state for an ideal gas without ionization.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cp; cp1; cp1tilde
! PENCILS PROVIDED glnTT(3); TT; TT1; gTT(3); yH; hss(3,3); hlnTT(3,3)
! PENCILS PROVIDED del2ss; del6ss; del2lnTT; cv1; del6lnTT; gamma; lncp
! PENCILS PROVIDED del2TT; del6TT
!
!***************************************************************

module EquationOfState

  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'eos.h'

  interface eoscalc ! Overload subroutine `eoscalc' function
    module procedure eoscalc_pencil   ! explicit f implicit m,n
    module procedure eoscalc_point    ! explicit lnrho, ss
    module procedure eoscalc_farray   ! explicit lnrho, ss
  end interface

  interface pressure_gradient ! Overload subroutine `pressure_gradient'
    module procedure pressure_gradient_farray  ! explicit f implicit m,n
    module procedure pressure_gradient_point   ! explicit lnrho, ss
  end interface

! integers specifying which independent variables to use in eoscalc
  integer, parameter :: ilnrho_ss=1,ilnrho_ee=2,ilnrho_pp=3
  integer, parameter :: ilnrho_lnTT=4,ilnrho_cs2=5
  integer, parameter :: irho_cs2=6, irho_ss=7, irho_lnTT=8, ilnrho_TT=9
  integer, parameter :: irho_TT=10

  integer :: iglobal_cs2, iglobal_glnTT

  ! secondary parameters calculated in initialize
!  real :: TT_ion=impossible,TT_ion_=impossible
!  real :: ss_ion=impossible,kappa0=impossible
!  real :: lnrho_H=impossible,lnrho_e=impossible,lnrho_e_=impossible
!  real :: lnrho_p=impossible,lnrho_He=impossible
!
  real :: lnTT0=impossible
!
!  initialize the helium fraction (by mass) to 0.
!  and the mean molecular weight mu to unity.
!
  real :: xHe=0.
  real :: mu=1.

  real :: cs0=1., rho0=1.
  real :: cs20=1., lnrho0=0.
  real :: ptlaw=3./4.
  real :: gamma=5./3.
  real :: Rgas_cgs=0., Rgas, error_cp=1e-6
  real :: gamma1    !(=gamma-1)
  real :: gamma11   !(=1/gamma)
  real :: cp=impossible, cp1=impossible, cv=impossible, cv1=impossible
  real :: cs2top_ini=impossible, dcs2top_ini=impossible
  real :: cs2bot=1., cs2top=1.
  real :: cs2cool=0.
  real :: mpoly=1.5, mpoly0=1.5, mpoly1=1.5, mpoly2=1.5
  real, dimension(3) :: beta_glnrho_global=0., beta_glnrho_scaled=0.
  integer :: isothtop=0

  integer :: ieosvars=-1, ieosvar1=-1, ieosvar2=-1, ieosvar_count=0

  logical :: leos_isothermal=.false., leos_isentropic=.false.
  logical :: leos_isochoric=.false., leos_isobaric=.false.
  logical :: leos_localisothermal=.false.

  ! input parameters
  namelist /eos_init_pars/ xHe, mu, cp, cs0, rho0, gamma, error_cp, ptlaw, &
    cs2top_ini, dcs2top_ini

  ! run parameters
  namelist /eos_run_pars/  xHe, mu, cp, cs0, rho0, gamma, error_cp, ptlaw, &
    cs2top_ini, dcs2top_ini

  contains

!***********************************************************************
    subroutine register_eos()
!
!  14-jun-03/axel: adapted from register_eos
!
      use Cdata
      use Sub
!
      leos=.true.
      leos_idealgas=.true.
!
      iyH = 0
      ilnTT = 0
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_eos: ionization nvar = ', nvar
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
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
!
      use Mpicomm, only: stop_it
!
      real :: Rgas_unit_sys=1., cp_reference
!
!  set gamma1, cs20, and lnrho0
!  (used currently for non-dimensional equation of state)
!
      gamma1=gamma-1.
      gamma11=1./gamma
!
!  avoid floating overflow if cs0 was not set:
!
      cs20=cs0**2
      lnrho0=log(rho0)
!
! Initialize variable selection code (needed for RELOADing)
!
      ieosvars=-1
      ieosvar_count=0
!
!  Unless unit_temperature is set, calculate by default with cp=1.
!  If unit_temperature is set, cp must follow from this.
!  Conversely, if cp is set, then unit_temperature must follow from this.
!  If unit_temperature and cp are set, the problem is overdetermined,
!    but it may still be correct, so this will be checked here.
!  When gamma=1. (gamma1=0.), write Rgas=mu*cp or cp=Rgas/mu.
!
      if (unit_system == 'cgs') then
         Rgas_unit_sys = k_B_cgs/m_u_cgs
      elseif (unit_system == 'SI') then
         Rgas_unit_sys = k_B_cgs/m_u_cgs*1.e-4
      endif
!
      if (unit_temperature == impossible) then
        if (cp == impossible) cp=1.
        if (gamma1 == 0.) then
          Rgas=mu*cp
        else
          Rgas=mu*gamma1*gamma11*cp
        endif
        unit_temperature=unit_velocity**2*Rgas/Rgas_unit_sys
      else
        Rgas=Rgas_unit_sys*unit_temperature/unit_velocity**2
        if (cp == impossible) then
          if (gamma1 == 0.) then
            cp=Rgas/mu
          else
            cp=Rgas/(mu*gamma1*gamma11)
          endif
        else
!
!  checking whether the units are overdetermined.
!  This is assumed to be the case when the to differ by error_cp
!
          if (gamma1 == 0.) then
            cp_reference=Rgas/mu
          else
            cp_reference=Rgas/(mu*gamma1*gamma11)
          endif
          if (abs(cp-cp_reference)/cp > error_cp) then
            if (lroot) print*,'initialize_eos: consistency: cp=',cp, &
               'while: cp_reference=',cp_reference
            call stop_it('initialize_eos')
          endif
        endif
      endif
      cp1=1./cp
      cv=gamma11*cp
      cv1=gamma*cp1
!
!  Need to calculate the equivalent of cs0
!  Distinguish between gamma=1 case and not.
!
      if (gamma1 /= 0.) then
        lnTT0=log(cs20/(cp*gamma1))  !(general case)
      else
        lnTT0=log(cs20/cp)  !(isothermal/polytropic cases: check!)
      endif
!
!  check that everything is OK
!
      if (lroot) then
        print*,'initialize_eos: unit_temperature=',unit_temperature
        print*,'initialize_eos: cp,lnTT0,cs0=',cp,lnTT0,cs0
      endif
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos()
!
      use Mpicomm, only: stop_it
!
! Initialize variable selection code (needed for RELOADing)
!
      ieosvars=-1
      ieosvar_count=0
!
!  write constants to disk. In future we may want to deal with this
!  using an include file or another module.
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
!        write (1,*) 'TT_ion=',TT_ion
!        write (1,*) 'TT_ion_=',TT_ion_
!        write (1,*) 'lnrho_e=',lnrho_e
!        write (1,*) 'lnrho_H=',lnrho_H
!        write (1,*) 'lnrho_p=',lnrho_p
!        write (1,*) 'lnrho_He=',lnrho_He
!        write (1,*) 'lnrho_e_=',lnrho_e_
!        write (1,*) 'ss_ion=',ss_ion
!        write (1,*) 'kappa0=',kappa0
        write (1,'(a,1pd26.16)') 'k_B=',k_B
        write (1,'(a,1pd26.16)') 'm_H=',m_H
        write (1,*) 'lnTTO=',lnTT0
        write (1,*) 'cp=',cp
        close (1)
      endif
!
    endsubroutine initialize_eos
!*******************************************************************
    subroutine select_eos_variable(variable,findex)
!
!  Select eos variable
!
!   02-apr-06/tony: implemented
!
      use FArrayManager

      character (len=*), intent(in) :: variable
      integer, intent(in) :: findex
      integer :: this_var=0
      integer, save :: ieosvar_selected=0
      integer, parameter :: ieosvar_lnrho = 2**0
      integer, parameter :: ieosvar_rho   = 2**1
      integer, parameter :: ieosvar_ss    = 2**2
      integer, parameter :: ieosvar_lnTT  = 2**3
      integer, parameter :: ieosvar_TT    = 2**4
      integer, parameter :: ieosvar_cs2   = 2**5
      integer, parameter :: ieosvar_pp    = 2**6
!
      if (ieosvar_count.eq.0) ieosvar_selected=0
!
      if (ieosvar_count.ge.2) &
        call fatal_error("select_eos_variable", &
             "2 thermodynamic quantities have already been defined while attempting to add a 3rd: ") !//variable)

      ieosvar_count=ieosvar_count+1

!      select case (variable)
      if (variable=='ss') then
          this_var=ieosvar_ss
          if (findex.lt.0) then
            leos_isentropic=.true.
          endif
      elseif (variable=='cs2') then
          this_var=ieosvar_cs2
          if (findex==-2) then
            leos_localisothermal=.true.

            call farray_register_global('cs2',iglobal_cs2)
            call farray_register_global('glnTT',iglobal_glnTT,vector=3)

          elseif (findex.lt.0) then
            leos_isothermal=.true.
          endif
      elseif (variable=='lnTT') then
          this_var=ieosvar_lnTT
          if (findex.lt.0) then
            leos_isothermal=.true.
          endif
      elseif (variable=='TT') then
          this_var=ieosvar_TT
      elseif (variable=='lnrho') then
          this_var=ieosvar_lnrho
          if (findex.lt.0) then
            leos_isochoric=.true.
          endif
      elseif (variable=='rho') then
          this_var=ieosvar_rho
          if (findex.lt.0) then
            leos_isochoric=.true.
          endif
      elseif (variable=='pp') then
          this_var=ieosvar_pp
          if (findex.lt.0) then
            leos_isobaric=.true.
          endif
      else
        call fatal_error("select_eos_variable", &
             "unknown thermodynamic variable")
      endif
      if (ieosvar_count==1) then
        ieosvar1=findex
        ieosvar_selected=ieosvar_selected+this_var
        return
      endif
!
! Ensure the indexes are in the correct order.
!
      if (this_var.lt.ieosvar_selected) then
        ieosvar2=ieosvar1
        ieosvar1=findex
      else
        ieosvar2=findex
      endif
      ieosvar_selected=ieosvar_selected+this_var
      select case (ieosvar_selected)
        case (ieosvar_lnrho+ieosvar_ss)
          if (lroot) print*,"select_eos_variable: Using lnrho and ss"
          ieosvars=ilnrho_ss
        case (ieosvar_rho+ieosvar_ss)
          if (lroot) print*,"select_eos_variable: Using rho and ss"
          ieosvars=irho_ss
        case (ieosvar_lnrho+ieosvar_lnTT)
          if (lroot) print*,"select_eos_variable: Using lnrho and lnTT"
          ieosvars=ilnrho_lnTT
        case (ieosvar_lnrho+ieosvar_TT)
          if (lroot) print*,"select_eos_variable: Using lnrho and TT"
          ieosvars=ilnrho_TT
        case (ieosvar_rho+ieosvar_lnTT)
          if (lroot) print*,"select_eos_variable: Using rho and lnTT"
          ieosvars=irho_lnTT
        case (ieosvar_lnrho+ieosvar_cs2)
          if (lroot) print*,"select_eos_variable: Using lnrho and cs2"
          ieosvars=ilnrho_cs2
        case (ieosvar_rho+ieosvar_cs2)
          if (lroot) print*,"select_eos_variable: Using rho and cs2",iproc
          ieosvars=irho_cs2
        case (ieosvar_rho+ieosvar_TT)
          if (lroot) print*,"select_eos_variable: Using rho and TT"
          ieosvars=irho_TT
        case default
          if (lroot) print*,"select_eos_variable: Thermodynamic variable combination, ieosvar_selected= ",ieosvar_selected
          call fatal_error("select_eos_variable", &
             "This thermodynamic variable combination is not implemented: ")
      endselect
!
    endsubroutine select_eos_variable
!*******************************************************************
    subroutine getmu(mu_tmp)
!
!  Calculate average particle mass in the gas relative to
!
!   12-aug-03/tony: implemented
!
      real, intent(out) :: mu_tmp

!  mu = mu_H * (1 - xHe) + mu_He * xHe
!     = mu_H + (mu_He-mu_H) * xHe
!  mu_H = 1.
!  mu_He = 4.0026 / 1.0079  (molar masses from a Periodic Table)
!        = 3.97
!
      if (mu == 0.) then
        mu_tmp=1.+2.97153*xHe
      else
        mu_tmp=mu
      endif
    endsubroutine getmu
!*******************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
!  Writes iyH and ilnTT to index.pro file
!
!  14-jun-03/axel: adapted from rprint_radiation
!  21-11-04/anders: moved diagnostics to entropy
!
      use Cdata
!
      logical :: lreset
      logical, optional :: lwrite
!
      if (NO_WARN) print*, lreset, lwrite  !(keep compiler quiet)
!
    endsubroutine rprint_eos
!***********************************************************************
    subroutine get_slices_eos(f,slices)
!
      use Sub, only: keep_compiler_quiet
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
!  EOS is a pencil provider but evolves nothing so it is unlokely that
!  it will require any pencils for it's own use.
!
    endsubroutine pencil_criteria_eos
!***********************************************************************
    subroutine pencil_interdep_eos(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_cp)) lpencil_in(i_cp1)=.true.
!
      select case (ieosvars)
!
      case (ilnrho_ss,irho_ss)
        if (lpencil_in(i_ee)) then
          lpencil_in(i_lnTT)=.true.
        endif
        if (lpencil_in(i_pp)) then
          lpencil_in(i_lnTT)=.true.
          lpencil_in(i_lnrho)=.true.
        endif
        if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_lnTT)) then
          lpencil_in(i_lnrho)=.true.
          lpencil_in(i_ss)=.true.
        endif
        if (lpencil_in(i_del2lnTT)) then
          if (ldensity_nolog) then
            lpencil_in(i_del2rho)=.true.
          else
            lpencil_in(i_del2lnrho)=.true.
          endif
          lpencil_in(i_del2ss)=.true.
        endif
        if (lpencil_in(i_glnTT)) then
          lpencil_in(i_glnrho)=.true.
          lpencil_in(i_gss)=.true.
        endif
        if (lpencil_in(i_hlnTT)) then
          lpencil_in(i_hss)=.true.
          lpencil_in(i_hlnrho)=.true.
        endif
        if (leos_isentropic) then
          if (lpencil_in(i_cs2)) lpencil_in(i_lnrho)=.true.
        elseif (leos_isothermal) then
          if (lpencil_in(i_ss)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_gss)) lpencil_in(i_glnrho)=.true.
          if (lpencil_in(i_hss)) lpencil_in(i_hlnrho)=.true.
          if (lpencil_in(i_del2ss)) then
            if (ldensity_nolog) then
              lpencil_in(i_del2rho)=.true.
            else
              lpencil_in(i_del2lnrho)=.true.
            endif
          endif
          if (lpencil_in(i_del6ss)) lpencil_in(i_del6lnrho)=.true.
        else
          if (lpencil_in(i_cs2)) then
            lpencil_in(i_lnrho)=.true.
            lpencil_in(i_ss)=.true.
          endif
        endif
!
      case (ilnrho_lnTT,irho_lnTT)
        if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_ss)) then
          lpencil_in(i_lnrho)=.true.
          lpencil_in(i_lnTT)=.true.
        endif
        if (lpencil_in(i_del2ss)) then
          if (ldensity_nolog) then
            lpencil_in(i_del2rho)=.true.
          else
            lpencil_in(i_del2lnrho)=.true.
          endif
          lpencil_in(i_del2lnTT)=.true.
        endif
        if (lpencil_in(i_gss)) then
          lpencil_in(i_glnrho)=.true.
          lpencil_in(i_glnTT)=.true.
        endif
        if (lpencil_in(i_hss)) then
          lpencil_in(i_hlnTT)=.true.
          lpencil_in(i_hlnrho)=.true.
        endif
        if (lpencil_in(i_ee)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_pp)) then
          lpencil_in(i_lnTT)=.true.
          lpencil_in(i_lnrho)=.true.
        endif
        if (leos_isentropic) then
          if (lpencil_in(i_lnTT)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_glnTT)) lpencil_in(i_glnrho)=.true.
          if (lpencil_in(i_hlnTT)) lpencil_in(i_hlnrho)=.true.
          if (lpencil_in(i_del2lnTT)) then
            if (ldensity_nolog) then
              lpencil_in(i_del2rho)=.true.
            else
              lpencil_in(i_del2lnrho)=.true.
            endif
        endif
          if (lpencil_in(i_cs2)) lpencil_in(i_lnrho)=.true.
        else
          if (lpencil_in(i_cs2)) lpencil_in(i_lnTT)=.true.
        endif
!
      case (ilnrho_cs2,irho_cs2)
        if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_ee)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_pp)) lpencil_in(i_rho)=.true.
        if (leos_isothermal) then
          if (lpencil_in(i_ss)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_del2ss)) then 
            if (ldensity_nolog) then
              lpencil_in(i_del2rho)=.true.
            else
              lpencil_in(i_del2lnrho)=.true.
            endif
          endif
          if (lpencil_in(i_gss)) lpencil_in(i_glnrho)=.true.
          if (lpencil_in(i_hss)) lpencil_in(i_hlnrho)=.true.
        else
          if (lpencil_in(i_pp)) lpencil_in(i_cs2)=.true.
        endif
!
      case (ilnrho_TT)
        if (lpencil_in(i_ss)) then
          lpencil_in(i_lnrho)=.true.
          lpencil_in(i_lnTT)=.true.
        endif
        if (lpencil_in(i_glnTT)) then
          lpencil_in(i_gTT)=.true.
          lpencil_in(i_TT1)=.true.
        endif
!
      case (irho_TT)
        if (lpencil_in(i_ss)) then
          lpencil_in(i_lnrho)=.true.
          lpencil_in(i_lnTT)=.true.
        endif
        if (lpencil_in(i_glnTT)) then
          lpencil_in(i_gTT)=.true.
          lpencil_in(i_TT1)=.true.
        endif
        if (lpencil_in(i_del2lnTT)) then
          lpencil_in(i_del2TT)=.true.
          lpencil_in(i_glnTT)=.true.
          lpencil_in(i_TT1)=.true.
        endif
!
      case default
        call fatal_error("pencil_interdep_eos","case not implemented yet")
      endselect
!
    endsubroutine pencil_interdep_eos
!***********************************************************************
    subroutine calc_pencils_eos(f,p)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  02-apr-06/tony: coded
!
      use Cparam
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
      real, dimension(nx) :: grhogrho,d2rho,tmp
      integer :: i
!
!  Convert del2lnrho to rho-only terms; done here to avoid
!  repeated calls to dot2
!
      if (ldensity_nolog) then
        call dot2(p%grho,grhogrho)
        d2rho=p%rho1*(p%del2rho+p%rho1*grhogrho)
      endif
!
! THE FOLLOWING 2 ARE CONCEPTUALLY WRONG
! FOR pretend_lnTT since iss actually contain lnTT NOT entropy!
! The code is not wrong however since this is correctly
! handled by the eos module.

      select case (ieosvars)
!
!  work out thermodynamic quantities for given lnrho or rho and ss
!
      case (ilnrho_ss,irho_ss)
        if (leos_isentropic) then
          if (lpencil(i_ss)) p%ss=0
          if (lpencil(i_gss)) p%gss=0
          if (lpencil(i_hss)) p%hss=0
          if (lpencil(i_del2ss)) p%del2ss=0
          if (lpencil(i_del6ss)) p%del6ss=0
          if (lpencil(i_cs2)) p%cs2=cs20*exp(gamma1*(p%lnrho-lnrho0))
        elseif (leos_isothermal) then
          if (lpencil(i_ss)) p%ss=-(cp-cv)*(p%lnrho-lnrho0)
          if (lpencil(i_gss)) p%gss=-(cp-cv)*p%glnrho
          if (lpencil(i_hss)) p%hss=-(cp-cv)*p%hlnrho
          if (lpencil(i_del2ss)) then
            if (ldensity_nolog) then
              p%del2ss=-(cp-cv)*d2rho
            else
              p%del2ss=-(cp-cv)*p%del2lnrho
            endif
          endif
          if (lpencil(i_del6ss)) p%del6ss=-(cp-cv)*p%del6lnrho
          if (lpencil(i_cs2)) p%cs2=cs20
        elseif (leos_localisothermal) then
          call fatal_error("calc_pencils_eos","leos_localisothermal not implemented for ilnrho_ss, try ilnrho_cs2")
        else
          if (lpencil(i_ss)) p%ss=f(l1:l2,m,n,ieosvar2)
          if (lpencil(i_gss)) call grad(f,ieosvar2,p%gss)
          if (lpencil(i_hss)) call g2ij(f,ieosvar2,p%hss)
          if (lpencil(i_del2ss)) call del2(f,ieosvar2,p%del2ss)
          if (lpencil(i_del6ss)) call del6(f,ieosvar2,p%del6ss)
          if (lpencil(i_cs2)) p%cs2=cs20*exp(cv1*p%ss+gamma1*(p%lnrho-lnrho0))
        endif
        if (lpencil(i_lnTT)) p%lnTT=lnTT0+cv1*p%ss+gamma1*(p%lnrho-lnrho0)
        if (lpencil(i_pp)) p%pp=(cp-cv)*exp(p%lnTT+p%lnrho)
        if (lpencil(i_ee)) p%ee=cv*exp(p%lnTT)
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
        if (lpencil(i_glnTT)) p%glnTT=gamma1*p%glnrho+cv1*p%gss
        if (lpencil(i_del2lnTT)) then
          if (ldensity_nolog) then
            p%del2lnTT=gamma1*d2rho+cv1*p%del2ss
          else
            p%del2lnTT=gamma1*p%del2lnrho+cv1*p%del2ss
          endif
        endif
        if (lpencil(i_hlnTT)) p%hlnTT=gamma1*p%hlnrho+cv1*p%hss
!
!  work out thermodynamic quantities for given lnrho or rho and lnTT
!
      case (ilnrho_lnTT,irho_lnTT)
        if (leos_isentropic) then
          if (lpencil(i_lnTT)) p%lnTT=gamma1*(p%lnrho-lnrho0)+lnTT0
          if (lpencil(i_glnTT)) p%glnTT=gamma1*p%glnrho
          if (lpencil(i_hlnTT)) p%hlnTT=gamma1*p%hlnrho
          if (lpencil(i_del2lnTT)) then
            if (ldensity_nolog) then
              p%del2lnTT=gamma1*d2rho
            else
              p%del2lnTT=gamma1*p%del2lnrho
            endif
          endif
          if (lpencil(i_cs2)) p%cs2=cs20*exp(gamma1*(p%lnrho-lnrho0))
        elseif (leos_isothermal) then
          if (lpencil(i_lnTT)) p%lnTT=lnTT0
          if (lpencil(i_glnTT)) p%glnTT=0
          if (lpencil(i_hlnTT)) p%hlnTT=0
          if (lpencil(i_del2lnTT)) p%del2lnTT=0
          if (lpencil(i_cs2)) p%cs2=cs20
        elseif (leos_localisothermal) then
          call fatal_error("calc_pencils_eos","leos_localisothermal not implemented for ilnrho_ss, try ilnrho_cs2")
        else
          if (lpencil(i_lnTT)) p%lnTT=f(l1:l2,m,n,ieosvar2)
          if (lpencil(i_glnTT)) call grad(f,ieosvar2,p%glnTT)
          if (lpencil(i_hlnTT)) call g2ij(f,ieosvar2,p%hlnTT)
          if (lpencil(i_del2lnTT)) call del2(f,ieosvar2,p%del2lnTT)
          if (lpencil(i_del6lnTT)) call del6(f,ieosvar2,p%del6lnTT)
          if (lpencil(i_cs2)) p%cs2=cp*exp(p%lnTT)*gamma1
        endif
        if (lpencil(i_ss)) p%ss=cv*(p%lnTT-lnTT0-gamma1*(p%lnrho-lnrho0))
        if (lpencil(i_pp)) p%pp=(cp-cv)*exp(p%lnTT+p%lnrho)
        if (lpencil(i_ee)) p%ee=cv*exp(p%lnTT)
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
        if (lpencil(i_gss)) p%gss=cv*(p%glnTT-gamma1*p%glnrho)
        if (lpencil(i_del2ss)) then
          if (ldensity_nolog) then
            p%del2ss=cv*(p%del2lnTT-gamma1*d2rho)
          else
            p%del2ss=cv*(p%del2lnTT-gamma1*p%del2lnrho)
          endif
        endif
        if (lpencil(i_hss)) p%hss=cv*(p%hlnTT-gamma1*p%hlnrho)
        if (lpencil(i_del6ss)) call fatal_error("calc_pencils_eos","del6ss not available for ilnrho_lnTT")
!
!  work out thermodynamic quantities for given lnrho or rho and TT
!
      case (ilnrho_TT,irho_TT)
        if (lpencil(i_TT))   p%TT=f(l1:l2,m,n,ieosvar2)
        if (lpencil(i_TT1))  p%TT1=1./f(l1:l2,m,n,ieosvar2)
        if (lpencil(i_lnTT)) p%lnTT=log(f(l1:l2,m,n,ieosvar2))
        if (lpencil(i_cs2))  p%cs2=cp*gamma1*f(l1:l2,m,n,ieosvar2)
        if (lpencil(i_gTT))  call grad(f,ieosvar2,p%gTT)
        if (lpencil(i_glnTT)) then
          do i=1,3;p%glnTT(:,i)=p%gTT(:,i)*p%TT1;enddo
        endif
        if (lpencil(i_del2TT)) call del2(f,ieosvar2,p%del2TT)
        if (lpencil(i_del2lnTT)) then 
          tmp=0.
          do i=1,3
            tmp=tmp+p%glnTT(:,i)**2
          enddo
          p%del2lnTT=p%del2TT*p%TT1 - tmp
        endif
        if (lpencil(i_del6TT)) call del6(f,ieosvar2,p%del6TT)
        if (lpencil(i_ss)) p%ss=cv*(p%lnTT-lnTT0-gamma1*(p%lnrho-lnrho0))
!
!  work out thermodynamic quantities for given lnrho or rho and cs2
!
      case (ilnrho_cs2,irho_cs2)
        if (leos_isentropic) then
          call fatal_error("calc_pencils_eos","leos_isentropic not implemented for ilnrho_cs2, try ilnrho_ss")
        elseif (leos_isothermal) then
          if (lpencil(i_cs2)) p%cs2=cs20
          if (lpencil(i_lnTT)) p%lnTT=lnTT0
          if (lpencil(i_glnTT)) p%glnTT=0
          if (lpencil(i_hlnTT)) p%hlnTT=0
          if (lpencil(i_del2lnTT)) p%del2lnTT=0
          if (lpencil(i_ss)) p%ss=-(cp-cv)*(p%lnrho-lnrho0)
          if (lpencil(i_del2ss)) then
            if (ldensity_nolog) then
              p%del2ss=-(cp-cv)*d2rho
            else 
              p%del2ss=-(cp-cv)*p%del2lnrho
            endif
          endif
          if (lpencil(i_gss)) p%gss=-(cp-cv)*p%glnrho
          if (lpencil(i_hss)) p%hss=-(cp-cv)*p%hlnrho
          if (lpencil(i_pp)) p%pp=gamma11*p%rho*cs20
        elseif (leos_localisothermal) then
          if (lpencil(i_cs2)) p%cs2=f(l1:l2,m,n,iglobal_cs2)
          if (lpencil(i_lnTT)) call fatal_error("calc_pencils_eos","temperature not needed for localisothermal")
          if (lpencil(i_glnTT)) p%glnTT=f(l1:l2,m,n,iglobal_glnTT:iglobal_glnTT+2)
          if (lpencil(i_hlnTT)) call fatal_error("calc_pencils_eos","no gradients yet for localisothermal")
          if (lpencil(i_del2lnTT)) call fatal_error("calc_pencils_eos","no gradients yet for localisothermal")
          if (lpencil(i_ss)) call fatal_error("calc_pencils_eos","entropy not needed for localisothermal")
          if (lpencil(i_del2ss)) call fatal_error("calc_pencils_eos","no gradients yet for localisothermal")
          if (lpencil(i_gss)) call fatal_error("calc_pencils_eos","entropy gradient not needed for localisothermal")
          if (lpencil(i_hss)) call fatal_error("calc_pencils_eos","no gradients yet for localisothermal")
          if (lpencil(i_pp)) p%pp=p%rho*p%cs2
        else
          call fatal_error("calc_pencils_eos","Full equation of state not implemented for ilnrho_cs2")
        endif
!
!  internal energy
!  For gamma=1, we use R/mu = c_p = c_v, thus ee = c_vT = R/mu T = p/rho = cs^2.
!
        if (lpencil(i_ee)) then
          if (gamma1 /= 0.) then
            p%ee=(gamma11/gamma1)*p%cs2
          else
            p%ee=p%cs2            
          endif
        endif
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
        if (lpencil(i_del6ss)) call fatal_error("calc_pencils_eos","del6ss not available for ilnrho_cs2")
      case default
        call fatal_error("calc_pencils_eos","case not implemented yet")
      endselect
!
!  inverse cv and cp values
!
       if (lpencil(i_cv1)) p%cv1=cv1
       if (lpencil(i_cp1)) p%cp1=cp1
       if (lpencil(i_cp))  p%cp=1/p%cp1
       if (lpencil(i_cp1tilde)) p%cp1tilde=cp1
!
    endsubroutine calc_pencils_eos
!***********************************************************************
    subroutine ioninit(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      if (NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
    real, dimension (mx,my,mz,mfarray) :: f
!
    if (NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)

      real, intent(in) :: EE,TT,yH
      real, intent(inout) :: rho

      rho = EE * cv1 / TT
      if (NO_WARN) print*,yH

    endsubroutine getdensity
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
    subroutine get_ptlaw(ptlaw_)
!
!  04-jul-07/wlad: return the value of ptlaw to outside modules
!                  ptlaw is temperature gradient in accretion disks
!
      real, intent(out) :: ptlaw_
!
      ptlaw_=ptlaw
!
    endsubroutine get_ptlaw
!***********************************************************************
    subroutine isothermal_density_ion(pot,tmp)
!
      real, dimension (nx), intent(in) :: pot
      real, dimension (nx), intent(out) :: tmp
!
      tmp=pot
!
    endsubroutine isothermal_density_ion
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      use Cdata
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: cs2,cp1tilde
      real, dimension(nx) :: lnrho,ss
!
      if (ldensity_nolog) then
        lnrho=log(f(l1:l2,m,n,ilnrho))
      else
        lnrho=f(l1:l2,m,n,ilnrho)
      endif
      ss=f(l1:l2,m,n,iss)
!
!  pretend_lnTT
!
      if (pretend_lnTT) then
        cs2=gamma1*exp(cv1*ss)
      else
        cs2=cs20*exp(cv1*ss+gamma1*(lnrho-lnrho0))
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
      use Cdata
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: cs2,cp1tilde
!
!  pretend_lnTT
!
      if (pretend_lnTT) then
        cs2=gamma1*exp(gamma*cp1*ss)
      else
        cs2=cs20*exp(cv1*ss+gamma1*(lnrho-lnrho0))
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
      use Cdata
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(in) :: glnrho,gss
      real, dimension(nx,3), intent(out) :: glnTT
!
      if (gamma1==0.) call fatal_error('temperature_gradient', &
        'gamma=1 not allowed with entropy turned on!')
!
!  pretend_lnTT
!
      if (pretend_lnTT) then
        glnTT=gss
      else
        glnTT=gamma1*glnrho+cv1*gss
      endif
!
      if (NO_WARN) print*,f !(keep compiler quiet)
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
      use Cdata
      use Sub, only: dot2
!
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx) :: tmp
!
      if (gamma1==0.) &
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
          p%del2lnTT=gamma1*p%rho1*(p%del2rho+p%rho1*tmp)
        else
          p%del2lnTT=gamma1*p%del2lnrho+p%cv1*p%del2ss
        endif
      endif
!
      if (NO_WARN) print*,f !(keep compiler quiet)
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
      use Cdata
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3,3), intent(in) :: hlnrho,hss
      real, dimension(nx,3,3), intent(out) :: hlnTT
!
      if (gamma1==0.) call fatal_error('temperature_hessian','gamma=1 not allowed w/entropy')
!
!  pretend_lnTT
!
      if (pretend_lnTT) then
        hlnTT=hss
      else
        hlnTT=gamma1*hlnrho+cv1*hss
      endif
!
      if (NO_WARN) print*,f !(keep compiler quiet)
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
          if (pretend_lnTT) then
            f(l1:l2,m,n,iss)=log(cv1*ee)
          else
            f(l1:l2,m,n,iss)=cv*(log(cv1*ee)-lnTT0-gamma1*(lnrho_-lnrho0))
          endif
        elseif (present(pp)) then
          if (pretend_lnTT) then
            f(l1:l2,m,n,iss)=log(gamma*pp/(gamma1*lnrho_))
          else
            f(l1:l2,m,n,iss)=cv*(log(gamma*pp/gamma1)-gamma*lnrho_-gamma1*lnrho0-lnTT0)
          endif
        elseif (present(ss)) then
          if (pretend_lnTT) then
            f(l1:l2,m,n,iss)=lnTT0+cv1*ss+gamma1*(lnrho_-lnrho0)
          else
            f(l1:l2,m,n,iss)=ss
          endif
        endif

      elseif (psize==mx) then
        lnrho_=f(:,m,n,ilnrho)
        if (present(ee)) then
          if (pretend_lnTT) then
            f(:,m,n,iss)=log(cv1*ee)
          else
            f(:,m,n,iss)=cv*(log(cv1*ee)-lnTT0-gamma1*(lnrho_-lnrho0))
          endif
        elseif (present(pp)) then
          if (pretend_lnTT) then
            f(:,m,n,iss)=log(gamma*pp/(gamma1*lnrho_))
          else
            f(:,m,n,iss)=cv*(log(gamma*pp/gamma1)-gamma*lnrho_-gamma1*lnrho0-lnTT0)
          endif
        elseif (present(ss)) then
          if (pretend_lnTT) then
            f(:,m,n,iss)=lnTT0+cv1*ss+gamma1*(lnrho_-lnrho0)
          else
            f(:,m,n,iss)=ss
          endif
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
      use Cdata
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
!      if (gamma1==0.) call fatal_error('eoscalc_farray','gamma=1 not allowed w/entropy')
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
            ss_=-cv*gamma1*(lnrho_-lnrho0)
          else
            ss_=f(l1:l2,m,n,ieosvar2)
          endif
        case (mx)
          lnrho_=f(:,m,n,ieosvar1)
          if (leos_isentropic) then
            ss_=0
          elseif (leos_isothermal) then
            ss_=-cv*gamma1*(lnrho_-lnrho0)
          else
            ss_=f(:,m,n,ieosvar2)
          endif
        case default
          call fatal_error('eoscalc_farray','no such pencil size')
        end select

        lnTT_=lnTT0+cv1*ss_+gamma1*(lnrho_-lnrho0)
        if (gamma1==0.) &
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
            cs2_=exp(gamma1*(lnrho_-lnrho0)+log(cs20))
          elseif (leos_isothermal) then
            cs2_=cs20
          elseif (leos_localisothermal) then
            cs2_=f(l1:l2,m,n,iglobal_cs2)
          else
            call fatal_error('eoscalc_farray','full eos for cs2 not implemented')
          endif
        case (mx)
          lnrho_=f(:,m,n,ieosvar1)
          if (leos_isentropic) then
            cs2_=exp(gamma1*(lnrho_-lnrho0)+log(cs20))
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
        if (present(ee)) ee=gamma11*cs2_/gamma1
        if (present(pp)) pp=gamma11*cs2_*exp(lnrho_)
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
!   27-mar-06/tony: Introduces cv, cv1, gamma11 to make faster
!                   + more explicit
!   31-mar-06/tony: I removed messy lcalc_cp stuff completely. cp=1.
!                   is just fine.
!   22-jun-06/axel: reinstated cp,cp1,cv,cv1 in hopefully all the places.
!
      use Cdata
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp,cs2
      real :: lnrho_,ss_,lnTT_,ee_,pp_,cs2_,TT_
!
      if (gamma1==0.) call fatal_error('eoscalc_point','gamma=1 not allowed w/entropy')
!
      select case (ivars)

      case (ilnrho_ss)
        lnrho_=var1
        ss_=var2
        lnTT_=lnTT0+cv1*ss_+gamma1*(lnrho_-lnrho0)
        ee_=cv*exp(lnTT_)
        pp_=(cp-cv)*exp(lnTT_+lnrho_)
        cs2_=gamma*gamma1*ee_

      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
        lnTT_=log(cv1*ee_)
        ss_=cv*(lnTT_-lnTT0-gamma1*(lnrho_-lnrho0))
        pp_=gamma1*ee_*exp(lnrho_)
        cs2_=gamma*gamma1*ee_

      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
        ss_=cv*(log(pp_*exp(-lnrho_)*gamma/cs20)-gamma1*(lnrho_-lnrho0))
        ee_=pp_*exp(-lnrho_)/gamma1
        lnTT_=log(cv1*ee_)
        cs2_=gamma*gamma1*ee_

      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
        ss_=cv*(lnTT_-lnTT0-gamma1*(lnrho_-lnrho0))
        ee_=cv*exp(lnTT_)
        pp_=ee_*exp(lnrho_)*gamma1
        cs2_=gamma*gamma1*ee_

      case (ilnrho_TT)
        lnrho_=var1
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*exp(lnrho_)*gamma1
        cs2_=cp*gamma1*TT_

      case (irho_TT)
        lnrho_=log(var1)
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*var1*gamma1
        cs2_=cp*gamma1*TT_

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
!   27-mar-06/tony: Introduces cv, cv1, gamma11 to make faster
!                   + more explicit
!   31-mar-06/tony: I removed messy lcalc_cp stuff completely. cp=1.
!                   is just fine.
!   22-jun-06/axel: reinstated cp,cp1,cv,cv1 in hopefully all the places.
!
      use Cdata
!
      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), intent(out), optional :: lnrho,ss
      real, dimension(nx), intent(out), optional :: yH,lnTT
      real, dimension(nx), intent(out), optional :: ee,pp,cs2
      real, dimension(nx) :: lnrho_,ss_,lnTT_,ee_,pp_,cs2_,TT_
!
      if (gamma1==0.) call fatal_error('eoscalc_pencil','gamma=1 not allowed w/entropy')
!
      select case (ivars)

      case (ilnrho_ss)
        lnrho_=var1
        ss_=var2
        lnTT_=lnTT0+cv1*ss_+gamma1*(lnrho_-lnrho0)
        ee_=cv*exp(lnTT_)
        pp_=(cp-cv)*exp(lnTT_+lnrho_)
        cs2_=gamma*gamma1*ee_
        cs2_=cs20*cv1*ee_

      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
        lnTT_=log(cv1*ee_)
        ss_=cv*(lnTT_-lnTT0-gamma1*(lnrho_-lnrho0))
        pp_=gamma1*ee_*exp(lnrho_)
        cs2_=gamma*gamma1*ee_

      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
        ss_=cv*(log(pp_*exp(-lnrho_)*gamma/cs20)-gamma1*(lnrho_-lnrho0))
        ee_=pp_*exp(-lnrho_)/gamma1
        lnTT_=log(cv1*ee_)
        cs2_=gamma*gamma1*ee_

      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
        ss_=cv*(lnTT_-lnTT0-gamma1*(lnrho_-lnrho0))
        ee_=cv*exp(lnTT_)
        pp_=ee_*exp(lnrho_)*gamma1
        cs2_=gamma*gamma1*ee_

      case (ilnrho_TT)
        lnrho_=var1
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*exp(lnrho_)*gamma1
        cs2_=cp*gamma1*TT_

      case (irho_TT)
        lnrho_=log(var1)
        TT_=var2
        ss_=cv*(log(TT_)-lnTT0-gamma1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*var1*gamma1
        cs2_=cp*gamma1*TT_

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
      cs2=gamma1*cp*exp(lnTT)
!
    endsubroutine get_soundspeed
!***********************************************************************
    subroutine read_eos_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=eos_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=eos_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_eos_init_pars
!***********************************************************************
    subroutine write_eos_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=eos_init_pars)

    endsubroutine write_eos_init_pars
!***********************************************************************
    subroutine read_eos_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=eos_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=eos_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_eos_run_pars
!***********************************************************************
    subroutine write_eos_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=eos_run_pars)

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
      use Cdata
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
        lnrho=f(l1:l2,m,n,ilnrho)
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
      call get_soundspeed(log(T0),cs2bot)
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
      if (NO_WARN) print*,f,T0,rho0
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

      call fatal_error('Hminus_opacity',"opacity_type='Hminus' may not be used with noionization")

      if (NO_WARN) then
        kapparho=0
        print*,f
      endif

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
      use Cdata
      use Gravity
      use SharedVariables, only: get_shared_variable
      use Mpicomm, only: stop_it
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
      select case(topbot)
!
!  bottom boundary
!  ===============
!
      case('bot')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma1*(f(:,:,n1,ilnrho)-lnrho0)+cv1*f(:,:,n1,iss))
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
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+(cp-cv)* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case('top')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0
        endif
!
!  calculate Ftop/(K*cs2)
!
        rho_xy=exp(f(:,:,n2,ilnrho))
        cs2_xy=cs20*exp(gamma1*(f(:,:,n2,ilnrho)-lnrho0)+cv1*f(:,:,n2,iss))
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
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
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
    subroutine bc_ss_flux_tmp(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!
      use Cdata
      use Gravity
      use SharedVariables, only: get_shared_variable
      use Mpicomm, only: stop_it
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      logical, pointer :: lmultilayer, lheatc_chiconst
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy,lnrho_xy,ss_xy
      real, dimension (mx,my) :: cs2_xy1,cs2_xy2,T_xy,T_xy1,T_xy2,Told4
      real :: eps
      integer :: i,ierr,iter,niter=4,j,k
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
      select case(topbot)
!
!  bottom boundary
!  ===============
!
      case('bot')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma1*(f(:,:,n1,ilnrho)-lnrho0)+cv1*f(:,:,n1,iss))
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
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+(cp-cv)* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case('top')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0
        endif
!
!  Compute Temperature at the first 3 levels inward
!
        lnrho_xy=f(:,:,n2,ilnrho)
        cs2_xy =cs20*exp(gamma1*(f(:,:,n2  ,ilnrho)-lnrho0)+cv1*f(:,:,n2  ,iss))
        cs2_xy1=cs20*exp(gamma1*(f(:,:,n2-1,ilnrho)-lnrho0)+cv1*f(:,:,n2-1,iss))
        cs2_xy2=cs20*exp(gamma1*(f(:,:,n2-2,ilnrho)-lnrho0)+cv1*f(:,:,n2-2,iss))
        T_xy=cs2_xy/(cp*gamma1)
        T_xy1=cs2_xy1/(cp*gamma1)
        T_xy2=cs2_xy2/(cp*gamma1)
!
!  calculate Ftop/(K*cs2)
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!
        if (lheatc_chiconst) then
          tmp_xy=Ftop/(rho_xy*chi*cs2_xy)
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

      case default
        call fatal_error('bc_ss_flux','invalid argument')
      endselect
!
    endsubroutine bc_ss_flux_tmp
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
      use Cdata
      use Gravity
      use SharedVariables, only: get_shared_variable
      use Mpicomm, only: stop_it
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      logical, pointer :: lmultilayer, lheatc_chiconst
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy,lnrho_xy,ss_xy
      real, dimension (mx,my) :: cs2_xy1,cs2_xy2,T_xy,T_xy1,T_xy2,Told4
      real :: eps
      integer :: i,ierr,iter,niter=4,j,k
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
      select case(topbot)
!
!  bottom boundary
!  ===============
!
      case('bot')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma1*(f(:,:,n1,ilnrho)-lnrho0)+cv1*f(:,:,n1,iss))
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
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+(cp-cv)* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case('top')
!
!  Set (dcs2/dz) / (dcs2/dz)_ini = (cs2/cs2top_ini)^4
!  Note that (dcs2/dz) = cs20*[(gamma-1)*dlnrho/dz + gamma*d(s/cp)/dz]
!  So, ds/dz = - (cp-cv)*dlnrho/dz + cv*(dcs2/dz)/cs20
!  calculate tmp_xy
!
        cs2_xy=cs20*exp(gamma1*(f(:,:,n2,ilnrho)-lnrho0)+cv1*f(:,:,n2,iss))
        tmp_xy=cv*dcs2top_ini/cs20*(cs2_xy/cs2top_ini)**4
!
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
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
    endsubroutine bc_ss_flux
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
      use Cdata
      use Gravity
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
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if ((bcz1(ilnrho) /= 'a2') .and. (bcz1(ilnrho) /= 'a3')) &
          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3.')
        if (ldebug) print*, &
                'bc_ss_temp_old: set bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) &
              print*,'bc_ss_temp_old: cannot have cs2bot = ', cs2bot, ' <= 0'
        tmp_xy = (-gamma1*(f(:,:,n1,ilnrho)-lnrho0) &
             + log(cs2bot/cs20)) / gamma
        f(:,:,n1,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n1-i,iss) = 2*tmp_xy - f(:,:,n1+i,iss)
        enddo

!
!  top boundary
!
      case('top')
        if ((bcz1(ilnrho) /= 'a2') .and. (bcz1(ilnrho) /= 'a3')) &
          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3.')
        if (ldebug) print*, &
                   'bc_ss_temp_old: set top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                   'bc_ss_temp_old: cannot have cs2top = ',cs2top, ' <= 0'
  !     if (bcz1(ilnrho) /= 'a2') &
  !          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 4.')
        tmp_xy = (-gamma1*(f(:,:,n2,ilnrho)-lnrho0) &
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
      use Cdata
      use Gravity
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
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_x: set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_x: cannot have cs2bot<=0'
        if (lentropy .and. .not. pretend_lnTT) then 
           tmp = 2*cv*log(cs2bot/cs20)
           f(l1,:,:,iss) = 0.5*tmp - (cp-cv)*(f(l1,:,:,ilnrho)-lnrho0)
           do i=1,nghost
              f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
               - (cp-cv)*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
           enddo
        elseif (lentropy .and. pretend_lnTT) then
           f(l1,:,:,iss) = log(cs2bot/gamma1)
           do i=1,nghost; f(l1-i,:,:,iss)=2*f(l1,:,:,iss)-f(l1+i,:,:,iss); enddo              
        elseif (ltemperature) then
           f(l1,:,:,ilnTT) = log(cs2bot/gamma1)
           do i=1,nghost; f(l1-i,:,:,ilnTT)=2*f(l1,:,:,ilnTT)-f(l1+i,:,:,ilnTT); enddo              
        endif
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                       'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                       'bc_ss_temp_x: cannot have cs2top<=0'
         if (lentropy .and. .not. pretend_lnTT) then 
            tmp = 2*cv*log(cs2top/cs20)
            f(l2,:,:,iss) = 0.5*tmp - (cp-cv)*(f(l2,:,:,ilnrho)-lnrho0)
            do i=1,nghost
               f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
                    - (cp-cv)*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
            enddo
        elseif (lentropy .and. pretend_lnTT) then
           f(l2,:,:,iss) = log(cs2top/gamma1)
           do i=1,nghost; f(l2+i,:,:,iss)=2*f(l2,:,:,iss)-f(l2-i,:,:,iss); enddo
        elseif (ltemperature) then
           f(l2,:,:,ilnTT) = log(cs2top/gamma1)
           do i=1,nghost; f(l2+i,:,:,ilnTT)=2*f(l2,:,:,ilnTT)-f(l2-i,:,:,ilnTT); enddo           
        endif

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
      use Cdata
      use Gravity
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
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
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
      case('top')
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
      use Cdata
      use Gravity
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
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_z: cannot have cs2bot = ', cs2bot, ' <= 0'
        if (lentropy .and. .not. pretend_lnTT) then 
           tmp = 2*cv*log(cs2bot/cs20)
           f(:,:,n1,iss) = 0.5*tmp - (cp-cv)*(f(:,:,n1,ilnrho)-lnrho0)
           do i=1,nghost
              f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
                   - (cp-cv)*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
           enddo
        elseif (lentropy .and. pretend_lnTT) then
            f(:,:,n1,iss) = log(cs2bot/gamma1)
            do i=1,nghost; f(:,:,n1-i,iss)=2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
        elseif (ltemperature) then
            if (ltemperature_nolog) then 
              f(:,:,n1,ilnTT) = cs2bot/gamma1
            else
              f(:,:,n1,ilnTT) = log(cs2bot/gamma1)
            endif
            do i=1,nghost; f(:,:,n1-i,ilnTT)=2*f(:,:,n1,ilnTT)-f(:,:,n1+i,ilnTT); enddo
        endif
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                   'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                   'bc_ss_temp_z: cannot have cs2top = ', cs2top, ' <= 0'
        if (lentropy .and. .not. pretend_lnTT) then 
           tmp = 2*cv*log(cs2top/cs20)
           f(:,:,n2,iss) = 0.5*tmp - (cp-cv)*(f(:,:,n2,ilnrho)-lnrho0)
           do i=1,nghost
              f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
                   - (cp-cv)*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
           enddo
        elseif (lentropy .and. pretend_lnTT) then
            f(:,:,n2,iss) = log(cs2top/gamma1)
            do i=1,nghost; f(:,:,n2+i,iss)=2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
        elseif (ltemperature) then
            if (ltemperature_nolog) then 
              f(:,:,n2,ilnTT) = cs2top/gamma1
            else
              f(:,:,n2,ilnTT) = log(cs2top/gamma1)
            endif
            do i=1,nghost; f(:,:,n2+i,ilnTT)=2*f(:,:,n2,ilnTT)-f(:,:,n2-i,ilnTT); enddo
        endif

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
      use Cdata
      use Gravity
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
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
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
      case('top')
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
      use Cdata
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
      select case(topbot)
!
!  bottom boundary
!
      case('top')
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
      case('bot')
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
      use Cdata
      use Gravity
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
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
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
      case('top')
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
      use Cdata
      use Gravity
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
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (cs2bot<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2bot<=0'
        do i=1,nghost
          f(l1-i,:,:,iss) = f(l1+i,:,:,iss) &
               + (cp-cv)*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2top<=0'
        do i=1,nghost
          f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
               + (cp-cv)*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
        enddo

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
      use Cdata
      use Gravity
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
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (cs2bot<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2bot<=0'
        do i=1,nghost
          f(:,m1-i,:,iss) = f(:,m1+i,:,iss) &
               + (cp-cv)*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2top<=0'
        do i=1,nghost
          f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
               + (cp-cv)*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
        enddo

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
      use Cdata
      use Gravity
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
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
          if (cs2bot<=0.) print*, &
                                  'bc_ss_stemp_z: cannot have cs2bot<=0'
          do i=1,nghost
             f(:,:,n1-i,iss) = f(:,:,n1+i,iss) &
                  + (cp-cv)*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
          enddo
!
!  top boundary
!
      case('top')
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
    subroutine bc_ss_energy(f,topbot)
!
!  boundary condition for entropy
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the entropy module
!  26-aug-2003/tony: distributed across ionization modules
!
      use Cdata
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
    select case(topbot)
!
! Bottom boundary
!
    case('bot')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma1*f(:,:,n1,ilnrho)+cv1*f(:,:,n1,iss))
      do i=1,nghost
         f(:,:,n1-i,iss)=cv*(-gamma1*f(:,:,n1-i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
!
! Top boundary
!
    case('top')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma1*f(:,:,n2,ilnrho)+cv1*f(:,:,n2,iss))
      do i=1,nghost
         f(:,:,n2+i,iss)=cv*(-gamma1*f(:,:,n2+i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
    case default
      call fatal_error('bc_ss_energy','invalid argument')
    endselect

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
      if (NO_WARN) print*,f,topbot
!
    endsubroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot,j)
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
      use Cdata
      use Gravity
      use Sub, only: div

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
      real, dimension (my,mz) :: cs2,gravterm,centterm,uphi
      real :: dlnrhodz, dssdz
      real :: potp,potm,rad,step
      integer :: i,j

      select case (topbot)

!
!  Bottom boundary
!
      case ('bot')
        do i=1,nghost

          cs2 = cs20
          call potential(R=x(l1-i),pot=potm)
          call potential(R=x(l1+i),pot=potp)

          gravterm= -(potm-potp)/cs2

          step=-2*i*dx
          rad=x(l1-i)
          uphi=f(l1-i,:,:,iuy)

          centterm= uphi**2 * step/(rad*cs2)
          if (ldensity_nolog) then
            f(l1-i,:,:,ilnrho)=f(l1+i,:,:,ilnrho)*exp(gravterm + centterm)
          else  
            f(l1-i,:,:,ilnrho)=f(l1+i,:,:,ilnrho) + gravterm + centterm
          endif
          
          !print*,'potentials',potm,potp,-(potm-potp)
          !print*,'centrifugal',f(l1-i,mpoint,npoint,iuy)**2 *step/rad
          !stop

        enddo

!
!  Top boundary
!
      case ('top')
        do i=1,nghost

          cs2 = cs20
          call potential(R=x(l2+i),pot=potp)
          call potential(R=x(l2-i),pot=potm)
 
          gravterm= -(potp-potm)/cs2

          step=2*i*dx
          rad=x(l2+i)
          uphi=f(l2+i,:,:,iuy)

          centterm= uphi**2 * step/(rad*cs2)          
          if (ldensity_nolog) then
            f(l2+i,:,:,ilnrho) = f(l2-i,:,:,ilnrho)*exp(gravterm + centterm)
          else
            f(l2+i,:,:,ilnrho) = f(l2-i,:,:,ilnrho) + gravterm + centterm
          endif

          !if (i==nghost) then
          !  print*,'potentials',potp,potm,-potp+potm,-(potp-potm)
          !  print*,'centrifugal',f(l2+i,mpoint,npoint,iuy)**2 *step/rad
          !  stop
          !endif
        enddo
            
      case default

      endselect

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
      use Cdata
      use Gravity
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

      select case (topbot)
!
!  Bottom boundary
!
      case ('bot')
!
        if (lentropy) then
!
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
          dssdz    = -gamma1*gravz/cs2_point
!
          do i=1,nghost
            f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) - 2*i*dz*dlnrhodz
            f(:,:,n1-i,iss   ) = f(:,:,n1+i,iss   ) - 2*i*dz*dssdz
          enddo
!
        else
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
              f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho)*exp(-(potm-potp)/cs2)
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
          call eoscalc(ilnrho_ss,f(l1,m1,n1,ilnrho),f(l1,m1,n1,iss), &
              cs2=cs2_point)
!
          dlnrhodz =  gamma *gravz/cs2_point
          dssdz    = -gamma1*gravz/cs2_point
!
          do i=1,nghost
            f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) + 2*i*dz*dlnrhodz
            f(:,:,n2+i,iss   ) = f(:,:,n2-i,iss   ) + 2*i*dz*dssdz
          enddo
!
        else
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
              f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho)*exp(-(potp-potm)/cs2)
            else
              f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) - (potp-potm)/cs2
            endif
          enddo

        endif

      case default

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
      use Cdata
      use Fourier, only: fourier_transform_xy_xy, fourier_transform_other
      use Gravity, only: potential

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot

      real, dimension (mx,my,nghost) :: ghost_zones
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
      select case(topbot)
!
!  Potential field condition at the bottom
!
      case('bot')

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
            tmp_re = f(l1:l2,m1:m2,n1+i,ilnrho)*exp(+pot/cs2bot)
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
            f(l1:l2,m1:m2,n1-i,ilnrho) = tmp_re*exp(-pot/cs2bot)
          else
            f(l1:l2,m1:m2,n1-i,ilnrho) = tmp_re - pot/cs2bot
          endif

        enddo
!
!  Potential field condition at the top
!
      case('top')

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
            tmp_re = f(l1:l2,m1:m2,n2-i,ilnrho)*exp(+pot/cs2top)
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
            f(l1:l2,m1:m2,n2+i,ilnrho) = tmp_re*exp(-pot/cs2top)
          else
            f(l1:l2,m1:m2,n2+i,ilnrho) = tmp_re - pot/cs2top
          endif

        enddo

      case default

        if (lroot) print*,"bc_lnrho_hydrostatic_z_smooth: invalid argument"

      endselect

    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_liso(f,topbot)
!
!  Potential field boundary condition
!
!  02-jul-07/wlad: Adapted from Tobi's bc_aa_pot2
!  Does the same thing as bc_lnrho_hdss_z_iso, but for a local isothermal
!  equation of state (as opposed to strictly isothermal).
!
      use Cdata
      use Fourier, only: fourier_transform_xy_xy, fourier_transform_other
      use Gravity, only: potential

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot

      real, dimension (mx,my,nghost) :: ghost_zones
      real, dimension (nx,ny) :: kx,ky,kappa,exp_fact
      real, dimension (nx,ny) :: tmp_re,tmp_im
      real, dimension (nx) :: pot,rr_cyl,rr_sph,cs2,tmp1,tmp2
      integer :: i,mm_noghost
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
!  Check whether we want to do top or bottom (this is processor dependent)
!
      select case(topbot)
!
!  Potential field condition at the bottom
!
      case('bot')

        do i=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          exp_fact = exp(-kappa*(z(n1+i)-z(n1-i)))
         
          do m=m1,m2
            mm_noghost=m-m1+1
            rr_cyl=sqrt(x(l1:l2)**2+y(m)**2)
            rr_sph=sqrt(x(l1:l2)**2+y(m)**2+z(n1+i)**2)
            cs2=cs20*rr_cyl**(-ptlaw)
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
            call potential(x(l1:l2),y(m),z(n1+i),POT=tmp1,RMN=rr_sph)
            call potential(x(l1:l2),y(m),z(n1+i),POT=tmp2,RMN=rr_cyl)
            pot=tmp1-tmp2
          !call potential(z=z(n1+i),pot=pot)
            
            if (ldensity_nolog) then
              tmp_re(:,mm_noghost) = f(l1:l2,m,n1+i,ilnrho)*exp(+pot/cs2)
            else
              tmp_re(:,mm_noghost) = f(l1:l2,m,n1+i,ilnrho) + pot/cs2
            endif
          enddo

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

          do m=m1,m2
            mm_noghost=m-m1+1
!          call potential(z=z(n1-i),pot=pot)
            rr_cyl=sqrt(x(l1:l2)**2+y(m)**2)
            rr_sph=sqrt(x(l1:l2)**2+y(m)**2+z(n1-i)**2)
            call potential(x(l1:l2),y(m),z(n1-i),POT=tmp1,RMN=rr_sph)
            call potential(x(l1:l2),y(m),z(n1-i),POT=tmp2,RMN=rr_cyl)
            pot=tmp1-tmp2
            cs2=cs20*rr_cyl**(-ptlaw)

            if (ldensity_nolog) then
              f(l1:l2,m,n1-i,ilnrho) = tmp_re(:,mm_noghost)*exp(-pot/cs2)
            else
              f(l1:l2,m,n1-i,ilnrho) = tmp_re(:,mm_noghost) - pot/cs2
            endif
          enddo

        enddo
!
!  Potential field condition at the top
!
      case('top')

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
!          call potential(z=z(n2-i),pot=pot)
          do m=m1,m2
            mm_noghost=m-m1+1
            rr_cyl=sqrt(x(l1:l2)**2+y(m)**2)
            rr_sph=sqrt(x(l1:l2)**2+y(m)**2+z(n2-i)**2)
            call potential(x(l1:l2),y(m),z(n2-i),POT=tmp1,RMN=rr_sph)
            call potential(x(l1:l2),y(m),z(n2-i),POT=tmp2,RMN=rr_cyl)
            pot=tmp1-tmp2
            cs2=cs20*rr_cyl**(-ptlaw)

            if (ldensity_nolog) then
              tmp_re(:,mm_noghost) = f(l1:l2,m,n2-i,ilnrho)*exp(+pot/cs2)
            else
              tmp_re(:,mm_noghost) = f(l1:l2,m,n2-i,ilnrho) + pot/cs2
            endif
          enddo
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

          do m=m1,m2
            mm_noghost=m-m1+1
            rr_cyl=sqrt(x(l1:l2)**2+y(m)**2)
            rr_sph=sqrt(x(l1:l2)**2+y(m)**2+z(n2+i)**2)
            call potential(x(l1:l2),y(m),z(n2+i),POT=tmp1,RMN=rr_sph)
            call potential(x(l1:l2),y(m),z(n2+i),POT=tmp2,RMN=rr_cyl)
            pot=tmp1-tmp2
            cs2=cs20*rr_cyl**(-ptlaw)

!          call potential(z=z(n2+i),pot=pot)
            if (ldensity_nolog) then
              f(l1:l2,m,n2+i,ilnrho) = tmp_re(:,mm_noghost)*exp(-pot/cs2)
            else
              f(l1:l2,m,n2+i,ilnrho) = tmp_re(:,mm_noghost) - pot/cs2
            endif
          enddo
        enddo
        
      case default

        if (lroot) print*,"bc_lnrho_hydrostatic_z_smooth: invalid argument"

      endselect

    endsubroutine bc_lnrho_hdss_z_liso
!***********************************************************************
    subroutine bc_lnrho_hds_z_liso(f,topbot)
!
!  Boundary condition for density
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
!
!  12-Jul-2006/dintrans: coded
!  18-Jul-2007/wlad: adapted for local isothermal equation of state
!
      use Cdata
      use Gravity

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (nx) :: potm,potp,tmp1,tmp2,rr_cyl,rr_sph,cs2
      character (len=3), intent (in) :: topbot
      real :: dlnrhodz, dssdz
      integer :: i
      
      select case (topbot)
!
!  Bottom boundary
!
      case ('bot')
!
        do i=1,nghost
          do m=m1,m2 
            rr_cyl=sqrt(x(l1:l2)**2+y(m)**2)
            rr_sph=sqrt(x(l1:l2)**2+y(m)**2+z(n1-i)**2)
!
            !call the arrays with potentials
            call potential(x(l1:l2),y(m),z(n1-i),POT=tmp1,RMN=rr_sph)
            call potential(x(l1:l2),y(m),z(n1-i),POT=tmp2,RMN=rr_cyl)
            potm=tmp1-tmp2
!
            rr_sph=sqrt(x(l1:l2)**2+y(m)**2+z(n1+i)**2)
            call potential(x(l1:l2),y(m),z(n1+i),POT=tmp1,RMN=rr_sph)
            call potential(x(l1:l2),y(m),z(n1+i),POT=tmp2,RMN=rr_cyl)
            potp=tmp1-tmp2
!
            cs2=cs20*rr_cyl**(-ptlaw)
            if (ldensity_nolog) then
              f(l1:l2,m,n1-i,ilnrho) = f(l1:l2,m,n1+i,ilnrho)*exp((potm-potp)/cs2)
            else
              f(l1:l2,m,n1-i,ilnrho) = f(l1:l2,m,n1+i,ilnrho) + (potm-potp)/cs2
            endif
          enddo
        enddo
!
!  Top boundary
!
      case ('top')
!
        do i=1,nghost
          do m=m1,m2 
            rr_cyl=sqrt(x(l1:l2)**2+y(m)**2)
            rr_sph=sqrt(x(l1:l2)**2+y(m)**2+z(n2-i)**2)
!            
            !call the arrays with potentials
            call potential(x(l1:l2),y(m),z(n2-i),POT=tmp1,RMN=rr_sph)
            call potential(x(l1:l2),y(m),z(n2-i),POT=tmp2,RMN=rr_cyl)
            potm=tmp1-tmp2
!            
            rr_sph=sqrt(x(l1:l2)**2+y(m)**2+z(n2+i)**2)
            call potential(x(l1:l2),y(m),z(n2+i),POT=tmp1,RMN=rr_sph)
            call potential(x(l1:l2),y(m),z(n2+i),POT=tmp2,RMN=rr_cyl)
            potp=tmp1-tmp2
!
            cs2=cs20*rr_cyl**(-ptlaw)
            if (ldensity_nolog) then
              f(l1:l2,m,n2+i,ilnrho) = f(l1:l2,m,n2-i,ilnrho)*exp(-(potp-potm)/cs2)
            else
              f(l1:l2,m,n2+i,ilnrho) = f(l1:l2,m,n2-i,ilnrho) - (potp-potm)/cs2
            endif
          enddo
        enddo
      case default
!        
      endselect
!
    endsubroutine bc_lnrho_hds_z_liso
!***********************************************************************
endmodule EquationOfState
