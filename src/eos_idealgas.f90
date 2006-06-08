! $Id: eos_idealgas.f90,v 1.48 2006-06-08 23:16:38 theine Exp $

!  Dummy routine for ideal gas

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss,gss,ee,pp,lnTT,cs2,cp1tilde,glnTT,TT,TT1
! PENCILS PROVIDED yH,hss,hlnTT,del2ss,del6ss,del2lnTT,cv1
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
  integer, parameter :: irho_cs2=6, irho_ss=7, irho_lnTT=8

  ! secondary parameters calculated in initialize
!  real :: TT_ion=impossible,TT_ion_=impossible
!  real :: ss_ion=impossible,kappa0=impossible 
!  real :: lnrho_H=impossible,lnrho_e=impossible,lnrho_e_=impossible
!  real :: lnrho_p=impossible,lnrho_He=impossible
!
  real :: lnTT0=impossible

  real :: xHe=0.   !0.1
  real :: mu=0.

  real :: cs0=1., rho0=1.
  real :: cs20=1., lnrho0=0.
  real :: gamma=5./3.
  real :: cp_cgs=0.
  real :: gamma1    ! gamma - 1.
  real :: gamma11   ! 1. / gamma
  real :: cs2bot=1., cs2top=1. 
  real :: cs2cool=0.
  real :: mpoly=1.5, mpoly0=1.5, mpoly1=1.5, mpoly2=1.5
  real, dimension(3) :: beta_glnrho_global=0., beta_glnrho_scaled=0.
  integer :: isothtop=0

  integer :: ieosvars=-1, ieosvar1=-1, ieosvar2=-1

  logical :: leos_isothermal=.false., leos_isentropic=.false.
  logical :: leos_isochoric=.false., leos_isobaric=.false.
  logical :: leos_localisothermal=.false.

  ! input parameters
  namelist /eos_init_pars/ xHe, mu, cp_cgs, cs0, rho0, gamma


  ! run parameters
  namelist /eos_run_pars/  xHe, mu, cp_cgs, cs0, rho0, gamma

  contains

!***********************************************************************
    subroutine register_eos()
!
!  14-jun-03/axel: adapted from register_eos
!
      use Cdata
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call fatal_error('register_eos','module registration called twice')
      first = .false.
!
      leos=.true.
      leos_idealgas=.true.
!
      iyH = 0
      ilnTT = 0

      if ((ip<=8) .and. lroot) then
        print*, 'register_eos: ionization nvar = ', nvar
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           '$Id: eos_idealgas.f90,v 1.48 2006-06-08 23:16:38 theine Exp $')
!
!  Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call fatal_error('register_eos','naux > maux')
      endif
!
    endsubroutine register_eos
!***********************************************************************
    subroutine initialize_eos()
!
      real :: mu_tmp, save_unit_temperature
!
!  set gamma1, cs20, and lnrho0
!  (used currently for non-dimensional equation of state)
!
      gamma1=gamma-1.
      gamma11=1./gamma
      ! avoid floating overflow if cs0 was not set:
      cs20=cs0**2
      lnrho0=log(rho0)
!
!
      ! Avoid setting unit_temperature=0 to avoid floating point exceptions
      save_unit_temperature=unit_temperature
      if (gamma1 /= 0.) then
        if (mu /= 0. .or. xHe /= 0.) then
          if (lroot) print*,'initialize_eos: Calculating unit_temperature based on mu or xHe'
          call getmu(mu_tmp)
          unit_temperature=unit_velocity**2*gamma1*mu_tmp/R_cgs*gamma11
        endif
      endif
      if (cp_cgs /= 0.) then
          if (lroot) print*,'initialize_eos: Calculating unit_temperature based on cp_cgs'
          unit_temperature=unit_velocity**2/cp_cgs
      endif
      if (lroot) print*,'initialize_eos: unit_temperature=',unit_temperature
!      if (abs(save_unit_temperature-unit_temperature) > 100*epsi) then
!        call fatal_error("initialize_eos", &
!               "unit_temperature specified does not match that calculated!")
!      endif

!
! Need to recalculate some constants
!
      if (gamma1 /= 0.) then
        lnTT0=log(cs20/gamma1)
      else                      ! gamma==1
        lnTT0=log(cs20)      ! Could the ionizers please check!
      endif
!   
!  write constants to disk. In future we may want to deal with this
!  using an include file or another module.
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro')
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
        close (1)
      endif
!
    endsubroutine initialize_eos
!*******************************************************************
    subroutine select_eos_variable(variable,findex)
!
!  Calculate average particle mass in the gas relative to
!
!   02-apr-06/tony: implemented
!
      character (len=*), intent(in) :: variable
      integer, intent(in) :: findex
      integer :: this_var=0
      integer, save :: ieosvar=0
      integer, save :: ieosvar_selected=0
      integer, parameter :: ieosvar_lnrho = 2**0
      integer, parameter :: ieosvar_rho   = 2**1
      integer, parameter :: ieosvar_ss    = 2**2
      integer, parameter :: ieosvar_lnTT  = 2**3
      integer, parameter :: ieosvar_TT    = 2**4
      integer, parameter :: ieosvar_cs2   = 2**5
      integer, parameter :: ieosvar_pp    = 2**6
!
!!      if (ieosvar.ge.2) &
!!        call fatal_error("select_eos_variable", &
!!             "2 thermodynamic quantities have already been defined while attempting to add a 3rd: ") !//variable)

      ieosvar=ieosvar+1

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
          elseif (findex.lt.0) then
            leos_isothermal=.true.
          endif
      elseif (variable=='lnTT') then
          this_var=ieosvar_lnTT
          if (findex.lt.0) then
            leos_isothermal=.true.
          endif
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
      if (ieosvar==1) then
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
        case (ieosvar_rho+ieosvar_lnTT)
          if (lroot) print*,"select_eos_variable: Using rho and lnTT" 
          ieosvars=irho_lnTT 
        case (ieosvar_lnrho+ieosvar_cs2)
          if (lroot) print*,"select_eos_variable: Using lnrho and cs2" 
          ieosvars=ilnrho_cs2
        case (ieosvar_rho+ieosvar_cs2)
          if (lroot) print*,"select_eos_variable: Using rho and cs2",iproc 
          ieosvars=irho_cs2 
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
      if(NO_WARN) print*,lreset  !(keep compiler quiet)   
!
    endsubroutine rprint_eos
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
      select case (ieosvars)
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
          lpencil_in(i_del2lnrho)=.true.
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
          if (lpencil_in(i_del2ss)) lpencil_in(i_del2lnrho)=.true.
          if (lpencil_in(i_del6ss)) lpencil_in(i_del6lnrho)=.true.
        else
          if (lpencil_in(i_cs2)) then
            lpencil_in(i_lnrho)=.true.
            lpencil_in(i_ss)=.true.
          endif
        endif
      case (ilnrho_lnTT,irho_lnTT)
        if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_ss)) then
          lpencil_in(i_lnrho)=.true.
          lpencil_in(i_lnTT)=.true.
        endif
        if (lpencil_in(i_del2ss)) then
          lpencil_in(i_del2lnrho)=.true.
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
          if (lpencil_in(i_del2lnTT)) lpencil_in(i_del2lnrho)=.true.
          if (lpencil_in(i_cs2)) lpencil_in(i_lnrho)=.true.
        else
          if (lpencil_in(i_cs2)) lpencil_in(i_lnTT)=.true.
        endif
      case (ilnrho_cs2,irho_cs2)
        if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_ee)) lpencil_in(i_lnTT)=.true.
        if (lpencil_in(i_pp)) then
          lpencil_in(i_lnTT)=.true.
          lpencil_in(i_lnrho)=.true.
        endif
        if (leos_isothermal) then
          if (lpencil_in(i_ss)) lpencil_in(i_lnrho)=.true.
          if (lpencil_in(i_del2ss)) lpencil_in(i_del2lnrho)=.true.
          if (lpencil_in(i_gss)) lpencil_in(i_glnrho)=.true.
          if (lpencil_in(i_hss)) lpencil_in(i_hlnrho)=.true.
        endif

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
!  02-04-06/tony: coded
!
      use Cparam
      use Global, only: get_global
      use Sub
!      
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!      
      intent(in) :: f
      intent(inout) :: p
!
! THE FOLLOWING 2 ARE CONCEPTUALLY WRONG 
! FOR pretend_lnTT since iss actually contain lnTT NOT entropy!
! The code is not wrong however since this is correctly
! handled by the eos module. 

      select case (ieosvars)
      case (ilnrho_ss,irho_ss)
        if (leos_isentropic) then
          if (lpencil(i_ss)) p%ss=0
          if (lpencil(i_gss)) p%gss=0
          if (lpencil(i_hss)) p%hss=0
          if (lpencil(i_del2ss)) p%del2ss=0
          if (lpencil(i_del6ss)) p%del6ss=0
          if (lpencil(i_cs2)) p%cs2=cs20*exp(gamma1*(p%lnrho-lnrho0)) 
        elseif (leos_isothermal) then
          if (lpencil(i_ss)) p%ss=-gamma1*gamma11*(p%lnrho-lnrho0)
          if (lpencil(i_gss)) p%gss=-gamma1*gamma11*p%glnrho
          if (lpencil(i_hss)) p%hss=-gamma1*gamma11*p%hlnrho
          if (lpencil(i_del2ss)) p%del2ss=-gamma1*gamma11*p%del2lnrho
          if (lpencil(i_del6ss)) p%del6ss=-gamma1*gamma11*p%del6lnrho
          if (lpencil(i_cs2)) p%cs2=cs20
        elseif (leos_localisothermal) then
          call fatal_error("calc_pencils_eos","leos_localisothermal not implemented for ilnrho_ss, try ilnrho_cs2")
        else
          if (lpencil(i_ss)) p%ss=f(l1:l2,m,n,ieosvar2)
          if (lpencil(i_gss)) call grad(f,ieosvar2,p%gss)
          if (lpencil(i_hss)) call g2ij(f,ieosvar2,p%hss)
          if (lpencil(i_del2ss)) call del2(f,ieosvar2,p%del2ss)
          if (lpencil(i_del6ss)) call del6(f,ieosvar2,p%del6ss)
          if (lpencil(i_cs2)) p%cs2=cs20*exp(gamma*p%ss+gamma1*(p%lnrho-lnrho0)) 
        endif
        if (lpencil(i_lnTT)) p%lnTT=lnTT0+gamma*p%ss+gamma1*(p%lnrho-lnrho0)
        if (lpencil(i_pp)) p%pp=gamma11*gamma1*exp(p%lnTT+p%lnrho)
        if (lpencil(i_ee)) p%ee=gamma11*exp(p%lnTT)
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
        if (lpencil(i_glnTT)) p%glnTT=gamma1*p%glnrho+gamma*p%gss
        if (lpencil(i_del2lnTT)) p%del2lnTT=gamma1*p%del2lnrho+gamma*p%del2ss
        if (lpencil(i_hlnTT)) p%hlnTT=gamma1*p%hlnrho+gamma*p%hss
      case (ilnrho_lnTT,irho_lnTT)
        if (leos_isentropic) then
          if (lpencil(i_lnTT)) p%lnTT=gamma1*(p%lnrho-lnrho0)+lnTT0
          if (lpencil(i_glnTT)) p%glnTT=gamma1*p%glnrho
          if (lpencil(i_hlnTT)) p%hlnTT=gamma1*p%hlnrho
          if (lpencil(i_del2lnTT)) p%del2lnTT=gamma1*p%del2lnrho
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
          if (lpencil(i_cs2)) p%cs2=exp(p%lnTT)*gamma1
        endif
        if (lpencil(i_ss)) p%ss=p%lnTT-lnTT0-gamma1*(p%lnrho-lnrho0)
        if (lpencil(i_pp)) p%pp=gamma11*gamma1*exp(p%lnTT+p%lnrho)
        if (lpencil(i_ee)) p%ee=gamma11*exp(p%lnTT)
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
        if (lpencil(i_gss)) p%gss=gamma11*(p%glnTT-gamma1*p%glnrho)
        if (lpencil(i_del2ss)) p%del2ss=gamma11*(p%del2lnTT-gamma1*p%del2lnrho)
        if (lpencil(i_hss)) p%hss=gamma11*(p%hlnTT-gamma1*p%hlnrho)
        if (lpencil(i_del6ss)) call fatal_error("calc_pencils_eos","del6ss not available for ilnrho_lnTT")
!
!
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
          if (lpencil(i_ss)) p%ss=-gamma1*(p%lnrho-lnrho0)*gamma11
          if (lpencil(i_del2ss)) p%del2ss=-gamma1*p%del2lnrho*gamma11
          if (lpencil(i_gss)) p%gss=-gamma1*p%glnrho*gamma11
          if (lpencil(i_hss)) p%hss=-gamma1*p%hlnrho*gamma11
        elseif (leos_localisothermal) then
          if (lpencil(i_cs2)) call get_global(p%cs2,m,n,'cs2')
          if (lpencil(i_lnTT)) p%lnTT=log(p%cs2/gamma1)
          if (lpencil(i_glnTT)) call fatal_error("calc_pencils_eos","no gradients yet for localisothermal") !p%glnTT=0
          if (lpencil(i_hlnTT)) call fatal_error("calc_pencils_eos","no gradients yet for localisothermal") 
          if (lpencil(i_del2lnTT)) call fatal_error("calc_pencils_eos","no gradients yet for localisothermal") 
          if (lpencil(i_ss)) p%ss=(p%lnTT-lnTT0-gamma1*(p%lnrho-lnrho0))*gamma11
          if (lpencil(i_del2ss)) call fatal_error("calc_pencils_eos","no gradients yet for localisothermal") 
          if (lpencil(i_gss)) call fatal_error("calc_pencils_eos","no gradients yet for localisothermal") 
          if (lpencil(i_hss)) call fatal_error("calc_pencils_eos","no gradients yet for localisothermal") 
        else
          call fatal_error("calc_pencils_eos","Full equation of state not implemented for ilnrho_cs2")
        endif
        if (lpencil(i_pp)) p%pp=gamma11*gamma1*exp(p%lnTT+p%lnrho)
        if (lpencil(i_ee)) p%ee=gamma11*p%cs2
        if (lpencil(i_yH)) p%yH=impossible
        if (lpencil(i_TT)) p%TT=exp(p%lnTT)
        if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
        if (lpencil(i_del6ss)) call fatal_error("calc_pencils_eos","del6ss not available for ilnrho_cs2")
        
      case default
        call fatal_error("calc_pencils_eos","case not implemented yet")
      endselect
!
       if (lpencil(i_cv1)) p%cv1=1.
       if (lpencil(i_cp1tilde)) p%cp1tilde=1.
! From noentropy.f90
!
!      if (lpencil(i_ee)) p%ee=p%cs2*gamma11/gamma1
!      if (lpencil(i_lnTT)) then
!        if (gamma==1. .or. cs20==0) then
!          p%TT1=0.
!        else
!          p%TT1=gamma1/p%cs2
!        endif
!      endif 
!

    endsubroutine calc_pencils_eos
!***********************************************************************
    subroutine ioninit(f)
!   
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
!   
      if(NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
    real, dimension (mx,my,mz,mvar+maux) :: f
!
    if(NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)
      
      real, intent(in) :: EE,TT,yH
      real, intent(inout) :: rho

      rho = EE * gamma / TT
      if (NO_WARN) print*,yH

    end subroutine getdensity
!***********************************************************************
    subroutine isothermal_density_ion(pot,tmp)
!
      real, dimension (nx), intent(in) :: pot
      real, dimension (nx), intent(out) :: tmp
!
      tmp=pot
!
    end subroutine isothermal_density_ion
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
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
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
        cs2=gamma1*exp(gamma*ss)
      else
        cs2=cs20*exp(gamma*ss+gamma1*(lnrho-lnrho0))
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
      cp1tilde=1.
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
        cs2=gamma1*exp(gamma*ss)
      else
        cs2=cs20*exp(gamma*ss+gamma1*(lnrho-lnrho0))
      endif
      cp1tilde=1.
!
    endsubroutine pressure_gradient_point
!***********************************************************************
    subroutine temperature_gradient(f,glnrho,gss,glnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
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
        glnTT=gamma1*glnrho+gamma*gss
      endif
!
      if (NO_WARN) print*,f !(keep compiler quiet)
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_laplacian(f,del2lnrho,del2ss,del2lnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx), intent(in) :: del2lnrho,del2ss
      real, dimension(nx), intent(out) :: del2lnTT
!
      if (gamma1==0.) call fatal_error('temperature_laplacian','gamma=1 not allowed w/entropy')
!
!  pretend_lnTT
!
      if (pretend_lnTT) then
        del2lnTT=del2ss
      else
        del2lnTT=gamma1*del2lnrho+gamma*del2ss
      endif
!
      if (NO_WARN) print*,f !(keep compiler quiet)
    endsubroutine temperature_laplacian
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally hlnPP and hlnTT
!   hP/rho=cs2*(hlnrho+cp1tilde*hss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
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
        hlnTT=gamma1*hlnrho+gamma*hss
      endif
!
      if (NO_WARN) print*,f !(keep compiler quiet)
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eosperturb(f,psize,ee,pp)
      
      real, dimension(mx,my,mz,mvar+maux), intent(inout) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(in), optional :: ee, pp
      real, dimension(psize) :: lnrho_

      if (psize==nx) then
        lnrho_=f(l1:l2,m,n,ilnrho)
        if (present(ee)) then
          if (pretend_lnTT) then
            f(l1:l2,m,n,iss)=log(gamma*ee)
          else
            f(l1:l2,m,n,iss)=(log(gamma*ee)-gamma1*(lnrho_-lnrho0)-lnTT0)*gamma11
          endif
        elseif (present(pp)) then
          if (pretend_lnTT) then
            f(l1:l2,m,n,iss)=log(gamma*pp/(gamma1*lnrho_))
          else
            f(l1:l2,m,n,iss)=(log(gamma*pp/gamma1)-gamma*lnrho_-gamma1*lnrho0-lnTT0)*gamma11
          endif
        endif
      else
        call not_implemented("eosperturb")
      endif
    end subroutine eosperturb
!***********************************************************************
    subroutine eoscalc_farray(f,psize,yH,lnTT,ee,pp,kapparho)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!
      use Cdata
      use Sub, only: max_mn_name, sum_mn_name
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      integer, intent(in) :: psize
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
            ss_=-gamma1*gamma11*(lnrho_-lnrho0)
          else
            ss_=f(l1:l2,m,n,ieosvar2)
          endif
        case (mx)
          lnrho_=f(:,m,n,ieosvar1)
          if (leos_isentropic) then
            ss_=0
          elseif (leos_isothermal) then
            ss_=-gamma1*gamma11*(lnrho_-lnrho0)
          else
            ss_=f(:,m,n,ieosvar2)
          endif
        case default
          call fatal_error('eoscalc_farray','no such pencil size')
        end select

        lnTT_=lnTT0+gamma*ss_+gamma1*(lnrho_-lnrho0)
        if (gamma1==0.) &
            call fatal_error('eoscalc_farray','gamma=1 not allowed w/entropy')
        if (present(lnTT)) lnTT=lnTT_
        if (present(ee)) &
            ee=gamma11*exp(lnTT_)
        if (present(pp)) &
            pp=gamma11*gamma1*exp(lnTT_+lnrho_)
!
! Log rho and Log T
!
      case (ilnrho_lnTT,irho_lnTT)
        select case (psize)
        case (nx)
          if (ieosvars==ilnrho_lnTT) then
            lnrho_=f(l1:l2,m,n,ieosvar1)
          else
            lnrho_=alog(f(l1:l2,m,n,ieosvar1))
          endif
          if (leos_isentropic) then
            ss_=0
          elseif (leos_isothermal) then
          else
            ss_=f(l1:l2,m,n,ieosvar2)
          endif
        case (mx)
          lnrho_=f(:,m,n,ieosvar1)
          if (leos_isentropic) then
            lnTT_=lnTT0+gamma1*gamma11*(lnrho_-lnrho0)
          elseif (leos_isothermal) then
            lnTT_=lnTT0
          else
            lnTT_=f(:,m,n,ieosvar2)
          endif
        case default
          call fatal_error('eoscalc_farray','no such pencil size')
        end select
!
        if (present(lnTT)) lnTT=lnTT_
        if (present(ee)) ee=gamma11*exp(lnTT_)
        if (present(pp)) pp=gamma11*gamma1*exp(lnTT_+lnrho_)
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
          else
            call fatal_error('eoscalc_farray','full eos for cs2 not implemented')
          endif
        case (mx)
          lnrho_=f(:,m,n,ieosvar1)
          if (leos_isentropic) then
            cs2_=exp(gamma1*(lnrho_-lnrho0)+log(cs20))
          elseif (leos_isothermal) then
            cs2_=cs20
          else
            call fatal_error('eoscalc_farray','full eos for cs2 not implemented')
          endif
        case default
          call fatal_error('eoscalc_farray','no such pencil size')
        end select
!
        if (present(lnTT)) lnTT=log(cs2_/gamma1)
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
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!   27-mar-06/tony: Introduces cv, cv1, gamma11 to make faster 
!                   + more explicit
!   31-mar-06/tony: I removed messy lcalc_cp stuff completely. cp=1. 
!                   is just fine.
!
      use Cdata
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp
      real :: lnrho_,ss_,lnTT_,ee_,pp_
!
      if (gamma1==0.) call fatal_error('eoscalc_point','gamma=1 not allowed w/entropy')
!
      select case (ivars)

      case (ilnrho_ss)
        lnrho_=var1
        ss_=var2
        lnTT_=lnTT0+gamma*ss_+gamma1*(lnrho_-lnrho0)
        ee_=gamma11*exp(lnTT_)
        pp_=gamma11*gamma1*exp(lnTT_+lnrho_)

      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
        ss_=gamma11*(log(ee_*gamma)-lnTT0-gamma1*(lnrho_-lnrho0))
        lnTT_=log(gamma11*ee_)
        pp_=gamma1*ee_*exp(lnrho_)

      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
        ss_=gamma11*(log(pp_*exp(-lnrho_)*gamma/cs20)-gamma1*(lnrho_-lnrho0))
        ee_=pp_*exp(-lnrho_)/gamma1
        lnTT_=log(gamma1*ee_)

      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
        ss_=gamma1*(lnTT_-lnTT0-gamma1*(lnrho_-lnrho0))
        ee_=gamma1*exp(lnTT_)
        pp_=ee_*exp(lnrho_)*gamma1

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
!
    endsubroutine eoscalc_point
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!   27-mar-06/tony: Introduces cv, cv1, gamma11 to make faster 
!                   + more explicit
!   31-mar-06/tony: I removed messy lcalc_cp stuff completely. cp=1. 
!                   is just fine.
!
      use Cdata
!
      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), intent(out), optional :: lnrho,ss
      real, dimension(nx), intent(out), optional :: yH,lnTT
      real, dimension(nx), intent(out), optional :: ee,pp
      real, dimension(nx) :: lnrho_,ss_,lnTT_,ee_,pp_
!
      if (gamma1==0.) call fatal_error('eoscalc_pencil','gamma=1 not allowed w/entropy')
!
      select case (ivars)

      case (ilnrho_ss)
        lnrho_=var1
        ss_=var2
        lnTT_=lnTT0+gamma*ss_+gamma1*(lnrho_-lnrho0)
        ee_=gamma11*exp(lnTT_)
        pp_=ee_*exp(lnrho_)*gamma1
        pp_=gamma11*gamma1*exp(lnTT_+lnrho_)

      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
        ss_=gamma11*(log(ee_*gamma)-lnTT0-gamma1*(lnrho_-lnrho0))
        lnTT_=log(gamma11*ee_)
        pp_=gamma1*ee_*exp(lnrho_)

      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
        ss_=gamma11*(log(pp_*exp(-lnrho_)*gamma/cs20)-gamma1*(lnrho_-lnrho0))
        ee_=pp_*exp(-lnrho_)/gamma1
        lnTT_=log(gamma*ee_)

      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
        ss_=gamma11*(lnTT_-lnTT0-gamma1*(lnrho_-lnrho0))
        ee_=gamma11*exp(lnTT_)
        pp_=ee_*exp(lnrho_)*gamma1

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
      cs2=gamma1*exp(lnTT)
!
    end subroutine get_soundspeed
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
      real, dimension(mx,my,mz,mvar+maux), intent(inout) :: f
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
      real, dimension(mx,my,mz,mvar+maux), intent(inout) :: f
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
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(mx,my,mz), intent(out) :: kapparho

      call fatal_error('Hminus_opacity',"opacity_type='Hminus' may not be used with noionization")

      if (NO_WARN) then
        kapparho=0
        print*,f
      endif

    endsubroutine Hminus_opacity
!***********************************************************************
    subroutine bc_ss_flux(f,topbot,hcond0,hcond1,Fheat,FheatK,chi, &
                lmultilayer,lcalc_heatcond_constchi)
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
!
      real, intent(in) :: Fheat, FheatK, hcond0, hcond1, chi
      logical, intent(in) :: lmultilayer, lcalc_heatcond_constchi
      
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy
      integer :: i
      
!
      if(ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!  ===============
!
      case('bot')
        if (lmultilayer) then
          if(headtt) print*,'bc_ss_flux: Fbot,hcond=',Fheat,hcond0*hcond1
        else
          if(headtt) print*,'bc_ss_flux: Fbot,hcond=',Fheat,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma1*(f(:,:,n1,ilnrho)-lnrho0)+gamma*f(:,:,n1,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy. 
!
        if(lcalc_heatcond_constchi) then
          tmp_xy=Fheat/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FheatK/cs2_xy
        endif
!
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+gamma1*gamma11* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case('top')
        if (lmultilayer) then
          if(headtt) print*,'bc_ss_flux: Ftop,hcond=',Fheat,hcond0*hcond1
        else
          if(headtt) print*,'bc_ss_flux: Ftop,hcond=',Fheat,hcond0
        endif
!
!  calculate Ftop/(K*cs2)
!
        rho_xy=exp(f(:,:,n2,ilnrho))
        cs2_xy=cs20*exp(gamma1*(f(:,:,n2,ilnrho)-lnrho0)+gamma*f(:,:,n2,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy. 
!
        if(lcalc_heatcond_constchi) then
          tmp_xy=Fheat/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FheatK/cs2_xy
        endif
!
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+gamma1*gamma11* &
              (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)-2*i*dz*tmp_xy)
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i
!
      if(ldebug) print*,'bc_ss_temp_old: ENTER - cs20,cs0=',cs20,cs0
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'bc_ss_temp_x: cs20,cs0=',cs20,cs0
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
                   'bc_ss_temp_x: set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_x: cannot have cs2bot<=0'
        tmp = 2*gamma11*log(cs2bot/cs20)
        f(l1,:,:,iss) = 0.5*tmp - gamma1*gamma11*(f(l1,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
               - gamma1*gamma11*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                       'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                       'bc_ss_temp_x: cannot have cs2top<=0'
        tmp = 2*gamma11*log(cs2top/cs20)
        f(l2,:,:,iss) = 0.5*tmp - gamma1*gamma11*(f(l2,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
               - gamma1*gamma11*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
        enddo

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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'bc_ss_temp_y: cs20,cs0=',cs20,cs0
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
        tmp = 2*gamma11*log(cs2bot/cs20)
        f(:,m1,:,iss) = 0.5*tmp - gamma1*gamma11*(f(:,m1,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m1-i,:,iss) = -f(:,m1+i,:,iss) + tmp &
               - gamma1*gamma11*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                     'bc_ss_temp_y: set y top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                     'bc_ss_temp_y: cannot have cs2top<=0'
        tmp = 2*gamma11*log(cs2top/cs20)
        f(:,m2,:,iss) = 0.5*tmp - gamma1*gamma11*(f(:,m2,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
               - gamma1*gamma11*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'bc_ss_temp_z: cs20,cs0=',cs20,cs0
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
        tmp = 2*gamma11*log(cs2bot/cs20)
        f(:,:,n1,iss) = 0.5*tmp - gamma1*gamma11*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
               - gamma1*gamma11*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                   'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                   'bc_ss_temp_z: cannot have cs2top = ', cs2top, ' <= 0'
        tmp = 2*gamma11*log(cs2top/cs20)
        f(:,:,n2,iss) = 0.5*tmp - gamma1*gamma11*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
               - gamma1*gamma11*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
        enddo
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'bc_lnrho_temp_z: cs20,cs0=',cs20,cs0
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
        tmp = 2*gamma11*log(cs2bot/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n1,iss) = 0.5*tmp - gamma1*gamma11*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n1-i,iss) = 2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2bot
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=-gravz/cs2bot
        do i=1,nghost
          f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) +f(:,:,n1+i,iss) &
                                                  -f(:,:,n1-i,iss) +2*i*dz*tmp
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                    'bc_lnrho_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*, &
                    'bc_lnrho_temp_z: cannot have cs2top<=0'
        tmp = 2*gamma11*log(cs2top/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n2,iss) = 0.5*tmp - gamma1*gamma11*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n2+i,iss) = 2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2top
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=gravz/cs2top
        do i=1,nghost
          f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) +f(:,:,n2-i,iss) &
                                                  -f(:,:,n2+i,iss) +2*i*dz*tmp
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
      if(ldebug) print*,'bc_lnrho_pressure_z: cs20,cs0=',cs20,cs0
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
        if(lentropy) then
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
          f(:,:,n2,ilnrho)=lnrho_top+ss_top-f(:,:,n2,iss)
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
        if(lentropy) then
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      if(ldebug) print*,'bc_ss_temp2_z: cs20,cs0=',cs20,cs0
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
        tmp = 1*gamma11*log(cs2bot/cs20)
        do i=0,nghost
          f(:,:,n1-i,iss) = tmp &
               - gamma1*gamma11*(f(:,:,n1-i,ilnrho)-lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                     'bc_ss_temp2_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp2_z: cannot have cs2top<=0'
        tmp = 1*gamma11*log(cs2top/cs20)
        do i=0,nghost
          f(:,:,n2+i,iss) = tmp &
               - gamma1*gamma11*(f(:,:,n2+i,ilnrho)-lnrho0)
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
      if(ldebug) print*,'bc_ss_stemp_x: cs20,cs0=',cs20,cs0
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
               + gamma1*gamma11*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2top<=0'
        do i=1,nghost
          f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
               + gamma1*gamma11*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
      if(ldebug) print*,'bc_ss_stemp_y: cs20,cs0=',cs20,cs0
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
               + gamma1*gamma11*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2top<=0'
        do i=1,nghost
          f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
               + gamma1*gamma11*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
      if(ldebug) print*,'bc_ss_stemp_z: cs20,cs0=',cs20,cs0
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
                  + gamma1*gamma11*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
          enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                 'bc_ss_stemp_z: cannot have cs2top<=0'
         do i=1,nghost
           f(:,:,n2+i,iss) = f(:,:,n2-i,iss) &
                + gamma1*gamma11*(f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho))
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      cs2_2d=cs20*exp(gamma1*f(:,:,n1,ilnrho)+gamma*f(:,:,n1,iss))
      do i=1,nghost
         f(:,:,n1-i,iss)=1.*gamma11*(-gamma1*f(:,:,n1-i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo

!
! Top boundary
!
    case('top')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma1*f(:,:,n2,ilnrho)+gamma*f(:,:,n2,iss))
      do i=1,nghost
         f(:,:,n2+i,iss)=1.*gamma11*(-gamma1*f(:,:,n2+i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
    case default
      call fatal_error('bc_ss_energy','invalid argument')
    endselect

    end subroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      call stop_it("bc_stellar_surface: NOT IMPLEMENTED IN EOS_IDEALGAS")
      if (NO_WARN) print*,f(1,1,1,1),topbot
!
    end subroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_stellar_surface_2(f,topbot)
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      call stop_it("bc_stellar_surface_2: NOT IMPLEMENTED IN EOS_IDEALGAS")
      if (NO_WARN) print*,f(1,1,1,1),topbot
!
    end subroutine bc_stellar_surface_2
!***********************************************************************
endmodule EquationOfState
