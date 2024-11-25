! $Id$
!
!  Thermodynamics with Fixed ionization fraction
!
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .true., leos_ionization = .true., leos_temperature_ionization=.false.
! CPARAM logical, parameter :: leos_idealgas = .false., leos_chemistry = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cp; cp1; cp1tilde
! PENCILS PROVIDED glnTT(3); TT; TT1; gTT(3); yH; hss(3,3); hlnTT(3,3)
! PENCILS PROVIDED rho_anel
! PENCILS PROVIDED del2ss; del6ss; del2lnTT; cv; cv1; glnmumol(3); ppvap; csvap2
! PENCILS PROVIDED rho1gpp(3)
!
!***************************************************************
module EquationOfState
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'eos.h'
  include 'eos_params.h'
!
! Constants use in calculation of thermodynamic quantities
!
  real :: lnTTss,lnTTlnrho,lnTT0
!
! secondary parameters calculated in initialize
  real :: TT_ion,lnTT_ion,TT_ion_,lnTT_ion_,Rgas,mu1yHxHe
  real :: ss_ion,ee_ion,kappa0,Srad0
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_p,lnrho_He
  real :: xHe_term,yH_term,one_yH_term
!
  real :: yH0=.0,xHe=0.1,xH2=0.,kappa_cst=1.
  character (len=labellen) :: opacity_type='ionized_H'
!
  namelist /eos_init_pars/ yH0,xHe,xH2,opacity_type,kappa_cst
!
  namelist /eos_run_pars/ yH0,xHe,xH2,opacity_type,kappa_cst
!ajwm  can't use impossible else it breaks reading param.nml
!ajwm  SHOULDN'T BE HERE... But can wait till fully unwrapped
  real :: cs0=1., rho0=1. 
  real :: cs20=1., lnrho0=0.
  real :: gamma=5./3., gamma_m1,gamma1, nabla_ad
!ajwm  can't use impossible else it breaks reading param.nml
  real :: cs2bot=1., cs2top=1.
  integer :: imass=0, ivars_mod
!
  real :: Cp_const=impossible
  real :: Pr_number=0.7
  logical :: lpres_grad = .false.
!
  contains
!***********************************************************************
    subroutine register_eos
!
!  14-jun-03/axel: adapted from register_ionization
!
      use SharedVariables, only: put_shared_variable
!
      iyH = 0
      ilnTT = 0
!
      if ((ip<=8) .and. lroot) print*, 'register_eos: ionization nvar = ', nvar
!
!  identify version number
!
      if (lroot) call svn_id( &
          "$Id$")
!
      call put_shared_variable('gamma',gamma,caller='register_eos')
!
      if (.not.ldensity) then
        call put_shared_variable('rho0',rho0)
        call put_shared_variable('lnrho0',lnrho0)
      endif
!
    endsubroutine register_eos
!*******************************************************************
    subroutine getmu(f,mu_tmp)
!
!  Calculate mean molecular weight of the gas
!
!   12-aug-03/tony: implemented
!   30-mar-04/anders: Added molecular hydrogen to ionization_fixed
!
      real, dimension (mx,my,mz,mfarray), optional :: f
      real, optional, intent(out) :: mu_tmp
!
      mu_tmp = (1.+3.97153*xHe)/(1-xH2+xHe)
!
! Complain if xH2 not between 0 and 0.5
!
      if (xH2 < 0. .or. xH2 > 0.5) call fatal_error('get_mu','xH2 must be <= 0.5 and >= 0.0')
!
      call keep_compiler_quiet(present(f))
!
    endsubroutine getmu
!***********************************************************************
    subroutine units_eos
!
!  dummy: here we don't allow for inputting cp.
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos(f)
!
!  Perform any post-parameter-read initialization, e.g. set derived
!  parameters.
!
!   2-feb-03/axel: adapted from Interstellar module
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: mu1yHxHe
!
      if (pretend_lnTT) then
        call warning('initialize_eos','pretend_lnTT is not used with ionization')
        pretend_lnTT=.false.
      endif
!
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
      if (lroot) print*,'initialize_eos: assume cp is not 1, yH0=',yH0

      Rgas=k_B/m_p
      mu1yHxHe=1.+3.97153*xHe
      TT_ion=chiH/k_B
      lnTT_ion=log(chiH/k_B)
      TT_ion_=chiH_/k_B
      lnTT_ion_=log(chiH_/k_B)
      lnrho_e=1.5*log((m_e/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_H=1.5*log((m_H/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_p=1.5*log((m_p/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_He=1.5*log((m_He/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_e_=1.5*log((m_e/hbar)*(chiH_/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      ss_ion=k_B/m_H/mu1yHxHe
      ee_ion=ss_ion*TT_ion
      gamma_m1=gamma-1.
      gamma1=1./gamma
      nabla_ad=gamma_m1/gamma
      kappa0=sigmaH_/m_H/mu1yHxHe/4.0
      Srad0=sigmaSB*TT_ion**4.0D0/pi
!
      if (lroot) then
        print*,'initialize_eos: reference values for ionization'
        print*,'initialize_eos: TT_ion,lnrho_e,ss_ion=',TT_ion,lnrho_e,ss_ion
      endif
!
      if (yH0>0.) then
        yH_term=yH0*(2*log(yH0)-lnrho_e-lnrho_p)
      elseif (yH0<0.) then
        call fatal_error('initialize_eos','yH0 must not be lower than zero')
      else
        yH_term=0.
      endif
!
      if (yH0<1.) then
        one_yH_term=(1.-yH0)*(log(1.-yH0)-lnrho_H)
      elseif (yH0>1.) then
        call fatal_error('initialize_eos','yH0 must not be greater than one')
      else
        one_yH_term=0.
      endif
!
      if (xHe>0.) then
        xHe_term=xHe*(log(xHe)-lnrho_He)
      elseif (xHe<0.) then
        call fatal_error('initialize_eos','xHe lower than zero makes no sense')
      else
        xHe_term=0.
      endif
!
! Set the reference sound speed (used for noionisation to impossible)
!
      lnTTss=(2./3.)/(1.+yH0+xHe-xH2)/ss_ion
      lnTTlnrho=2./3.
!
      lnTT0=lnTT_ion+(2./3.)*((yH_term+one_yH_term+xHe_term)/(1+yH0+xHe-xH2)-2.5)

      if (lroot) then
        print*,'initialize_eos: reference values for ionization'
        print*,'initialize_eos: TT_ion,ss_ion,kappa0=',TT_ion,ss_ion,kappa0
        print*,'initialize_eos: lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_=', &
                lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_
      endif
!
!  write scale non-free constants to file; to be read by idl
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,*) 'TT_ion=',TT_ion
        write (1,*) 'TT_ion_=',TT_ion_
        write (1,*) 'lnrho_e=',lnrho_e
        write (1,*) 'lnrho_H=',lnrho_H
        write (1,*) 'lnrho_p=',lnrho_p
        write (1,*) 'lnrho_He=',lnrho_He
        write (1,*) 'lnrho_e_=',lnrho_e_
        write (1,*) 'ss_ion=',ss_ion
        write (1,*) 'ee_ion=',ee_ion
        write (1,*) 'kappa0=',kappa0
        write (1,*) 'Srad0=',Srad0
        write (1,*) 'lnTTss=',lnTTss
        write (1,*) 'lnTTlnrho=',lnTTlnrho
        write (1,*) 'lnTT0=',lnTT0
        close (1)
      endif
!
    endsubroutine initialize_eos
!*******************************************************************
    subroutine select_eos_variable(variable,findex)
!
!   02-apr-06/tony: implemented
!
      character (len=*), intent(in) :: variable
      integer, intent(in) :: findex
!
      call keep_compiler_quiet(variable)
      call keep_compiler_quiet(findex)
!  DUMMY ideagas version below
!!      integer :: this_var=0
!!      integer, save :: ieosvar=0
!!      integer, save :: ieosvar_selected=0
!!      integer, parameter :: ieosvar_lnrho = 1
!!      integer, parameter :: ieosvar_rho   = 2
!!      integer, parameter :: ieosvar_ss    = 4
!!      integer, parameter :: ieosvar_lnTT  = 8
!!      integer, parameter :: ieosvar_TT    = 16
!!      integer, parameter :: ieosvar_cs2   = 32
!!      integer, parameter :: ieosvar_pp    = 64
!!!
!!      if (ieosvar>=2) &
!!        call fatal_error("select_eos_variable", &
!!             "2 thermodynamic quantities have already been defined while attempting to add a 3rd: ") !//variable)
!!
!!      ieosvar=ieosvar+1
!!
!!!      select case (variable)
!!      if (variable=='ss') then
!!          this_var=ieosvar_ss
!!          if (findex<0) then
!!            leos_isentropic=.true.
!!          endif
!!      elseif (variable=='cs2') then
!!          this_var=ieosvar_cs2
!!          if (findex==-2) then
!!            leos_localisothermal=.true.
!!          elseif (findex<0) then
!!            leos_isothermal=.true.
!!          endif
!!      elseif (variable=='lnTT') then
!!          this_var=ieosvar_lnTT
!!          if (findex<0) then
!!            leos_isothermal=.true.
!!          endif
!!      elseif (variable=='lnrho') then
!!          this_var=ieosvar_lnrho
!!          if (findex<0) then
!!            leos_isochoric=.true.
!!          endif
!!      elseif (variable=='rho') then
!!          this_var=ieosvar_rho
!!          if (findex<0) then
!!            leos_isochoric=.true.
!!          endif
!!      elseif (variable=='pp') then
!!          this_var=ieosvar_pp
!!          if (findex<0) then
!!            leos_isobaric=.true.
!!          endif
!!      else
!!        call fatal_error("select_eos_variable", &
!!             "unknown thermodynamic variable")
!!      endif
!!      if (ieosvar==1) then
!!        ieosvar1=findex
!!        ieosvar_selected=ieosvar_selected+this_var
!!        return
!!      endif
!!!
!!! Ensure the indexes are in the correct order.
!!!
!!      if (this_var<ieosvar_selected) then
!!        ieosvar2=ieosvar1
!!        ieosvar1=findex
!!      else
!!        ieosvar2=findex
!!      endif
!!      ieosvar_selected=ieosvar_selected+this_var
!!      select case (ieosvar_selected)
!!        case (ieosvar_lnrho+ieosvar_ss)
!!          ieosvars=ilnrho_ss
!!        case (ieosvar_lnrho+ieosvar_lnTT)
!!          ieosvars=ilnrho_lnTT
!!        case (ieosvar_lnrho+ieosvar_cs2)
!!          ieosvars=ilnrho_cs2
!!        case default
!!          print*,"select_eos_variable: Thermodynamic variable combination, ieosvar_selected= ",ieosvar_selected
!!          call fatal_error("select_eos_variable", &
!!             "This thermodynamic variable combination is not implemented: ")
!!      endselect
!
    endsubroutine select_eos_variable
!***********************************************************************
    subroutine pencil_criteria_eos
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
      if (lpencil_in(i_gTT)) then
        lpencil_in(i_glnTT)=.true.
        lpencil_in(i_TT)=.true.
      endif
      if (lpencil_in(i_del2lnTT)) then
        lpencil_in(i_del2lnrho)=.true.
        lpencil_in(i_del2ss)=.true.
      endif
      if (lpencil_in(i_glnTT)) then
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_gss)=.true.
      endif
      if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
      if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
!
      if (lpencil_in(i_hlnTT)) then
        lpencil_in(i_hss)=.true.
        if (.not.pretend_lnTT) lpencil_in(i_hlnrho)=.true.
      endif
!
    endsubroutine pencil_interdep_eos
!***********************************************************************
    subroutine calc_pencils_eos_std(f,p)
!
! Envelope adjusting calc_pencils_eos_pencpar to the standard use with
! lpenc_loc=lpencil
!
!  9-oct-15/MR: coded
!
      real, dimension (mx,my,mz,mfarray),intent(INOUT):: f
      type (pencil_case),                intent(OUT)  :: p
!
      call calc_pencils_eos_pencpar(f,p,lpencil)
!
    endsubroutine calc_pencils_eos_std
!***********************************************************************
    subroutine calc_pencils_eos_pencpar(f,p,lpenc_loc)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  02-04-06/tony: coded
!  09-10-15/MR: added mask parameter lpenc_loc.
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      intent(inout):: f
      intent(in)   :: lpenc_loc
      intent(out)  :: p
!
      integer :: i
!
! THE FOLLOWING 2 ARE CONCEPTUALLY WRONG
! FOR pretend_lnTT since iss actually contain lnTT NOT entropy!
! The code is not wrong however since this is correctly
! handled by the eos module.
! ss
      if (lpenc_loc(i_ss)) p%ss=f(l1:l2,m,n,iss)
! gss
      if (lpenc_loc(i_gss)) call grad(f,iss,p%gss)
! pp
      if (lpenc_loc(i_pp)) call eoscalc(f,nx,pp=p%pp)
! ee
      if (lpenc_loc(i_ee)) call eoscalc(f,nx,ee=p%ee)
! lnTT
      if (lpenc_loc(i_lnTT)) call eoscalc(f,nx,lnTT=p%lnTT)
! yH
      if (lpenc_loc(i_yH)) call eoscalc(f,nx,yH=p%yH)
! TT
      if (lpenc_loc(i_TT)) p%TT=exp(p%lnTT)
! TT1
      if (lpenc_loc(i_TT1)) p%TT1=exp(-p%lnTT)
! cs2 and cp1tilde
      if (lpenc_loc(i_cs2) .or. lpenc_loc(i_cp1tilde)) &
          call pressure_gradient(f,p%cs2,p%cp1tilde)
! glnTT
      if (lpenc_loc(i_glnTT)) call temperature_gradient(f,p%glnrho,p%gss,p%glnTT)
! gTT
      if (lpenc_loc(i_gTT)) then
        do i=1,3; p%gTT(:,i)=p%glnTT(:,i)*p%TT; enddo
      endif
! hss
      if (lpenc_loc(i_hss)) call g2ij(f,iss,p%hss)
! del2ss
      if (lpenc_loc(i_del2ss)) call del2(f,iss,p%del2ss)
! del2lnTT
      if (lpenc_loc(i_del2lnTT)) call temperature_laplacian(f,p)
! del6ss
      if (lpenc_loc(i_del6ss)) call del6(f,iss,p%del6ss)
! hlnTT
      if (lpenc_loc(i_hlnTT)) call temperature_hessian(f,p%hlnrho,p%hss,p%hlnTT)
!
      if (lpenc_loc(i_glnmumol)) p%glnmumol(:,:)=0.
!
!  This routine does not yet compute cv or cv1, but since those pencils
!  are supposed to be provided here, we better set them to impossible.
!
      if (lpenc_loc(i_cv1)) p%cv1=impossible
      if (lpenc_loc(i_cp1)) p%cp1=impossible
      if (lpenc_loc(i_cv))  p%cv=impossible
      if (lpenc_loc(i_cp))  p%cp=impossible
!
    endsubroutine calc_pencils_eos_pencpar
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)
!
      real, intent(in) :: EE,TT,yH
      real, intent(out) :: rho
      real :: lnrho

      lnrho = log(EE) - log(1.5*(1.+yH+xHe-xH2)*ss_ion*TT + yH*ee_ion)
!
      rho=exp(max(lnrho,-15.))
!
    endsubroutine getdensity
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
      real, dimension(nx), intent(out) :: cs2
      real, dimension(nx), intent(out), optional :: cp1tilde

      real, dimension(nx) :: lnrho,ss,lnTT
!
      lnrho=f(l1:l2,m,n,ilnrho)
      ss=f(l1:l2,m,n,iss)
      lnTT=lnTTss*ss+lnTTlnrho*lnrho+lnTT0
!
      cs2=gamma*(1+yH0+xHe-xH2)*ss_ion*exp(lnTT)
      if (present(cp1tilde)) cp1tilde=nabla_ad/(1+yH0+xHe-xH2)/ss_ion
!
    endsubroutine pressure_gradient_farray
!***********************************************************************
    subroutine get_gamma_etc(gamma,cp,cv)
!
      real, optional, intent(OUT) :: gamma, cp,cv
!
      if (headt) call warning('get_gamma_etc','gamma, cp, and cv are not constant in eos_fixed_ionization.'// &
                              achar(10)//'The values provided are for one-atomic ideal gas. Use at own risk')
      if (present(gamma)) gamma=5./3.
      if (present(cp)) cp=1.
      if (present(cv)) cv=3./5.

    endsubroutine get_gamma_etc
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
      real :: lnTT
!
      lnTT=lnTTss*ss+lnTTlnrho*lnrho+lnTT0
!
      cs2=gamma*(1+yH0+xHe-xH2)*ss_ion*exp(lnTT)
      cp1tilde=nabla_ad/(1+yH0+xHe-xH2)/ss_ion
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
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(in) :: glnrho,gss
      real, dimension(nx,3), intent(out) :: glnTT
!
      glnTT=(2.0/3.0)*(glnrho+gss/ss_ion/(1+yH0+xHe-xH2))
!
      call keep_compiler_quiet(f)
!
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally hlnPP and hlnTT
!   hP/rho=cs2*(hlnrho+cp1tilde*hss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3,3), intent(in) :: hlnrho,hss
      real, dimension(nx,3,3), intent(out) :: hlnTT
!
      hlnTT=(2.0/3.0)*(hlnrho+hss/ss_ion/(1+yH0+xHe-xH2))
!
      call keep_compiler_quiet(f)
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,yH,lnTT,ee,pp,cs2,kapparho)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(out), optional :: lnrho,yH,lnTT
      real, dimension(psize), intent(out), optional :: ee,pp,kapparho
      real, dimension(psize), optional :: cs2
      real, dimension(psize) :: lnrho_,ss_,lnTT_,TT_,yH_
!
      select case (psize)
!
      case (nx)
        lnrho_=f(l1:l2,m,n,ilnrho)
        ss_=f(l1:l2,m,n,iss)
      case (mx)
        lnrho_=f(:,m,n,ilnrho)
        ss_=f(:,m,n,iss)
      case default
        call fatal_error("eoscalc_farray","no such pencil size")
      end select
!
      lnTT_=lnTTss*ss_+lnTTlnrho*lnrho_+lnTT0
      if (present(ee) .or. present(pp) .or. present(kapparho)) TT_=exp(lnTT_)
      yH_=yH0
!
      if (present(lnrho)) lnrho=lnrho_
      if (present(yH))    yH=yH_
      if (present(lnTT))  lnTT=lnTT_
      if (present(ee))    ee=1.5*(1+yH_+xHe-xH2)*ss_ion*TT_+yH_*ss_ion*TT_ion
      if (present(pp))    pp=(1+yH_+xHe-xH2)*exp(lnrho_)*TT_*ss_ion
      if (present(cs2)) call not_implemented('eoscalc_farray','calculation of cs2')
!
!  Hminus opacity
!
      if (present(kapparho)) &
        kapparho=exp(2*lnrho_-lnrho_e+1.5*(lnTT_ion_-lnTT_)+TT_ion_/TT_)*yH_*(1-yH_)*kappa0
!
    endsubroutine eoscalc_farray
!***********************************************************************
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)

      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp,cs2

      if (.not.any((/ilnrho_ss, ilnrho_ee, ilnrho_pp, ilnrho_lnTT/)==ivars))  &
        call fatal_error('eoscalc_point','invalid combination of thermodynamic variables')

      ivars_mod=ivars

      call eoscalc_elem(var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)

    endsubroutine eoscalc_point
!***********************************************************************
    elemental subroutine eoscalc_elem(var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp,cs2
      real :: lnrho_,ss_,lnTT_,TT_,rho_,ee_,pp_
!
      select case (ivars_mod)
!
      case (ilnrho_ss)
        lnrho_ = var1
        ss_    = var2
        lnTT_  = lnTTss*ss_+lnTTlnrho*lnrho_+lnTT0
        if (present(ee) .or. present(pp)) then
          TT_  = exp(lnTT_)
          if (present(ee)) ee_ = 1.5*(1+yH0+xHe-xH2)*ss_ion*TT_+yH0*ee_ion
          if (present(pp)) pp_ = (1+yH0+xHe-xH2)*exp(lnrho_)*TT_*ss_ion
        endif
!
      case (ilnrho_ee)
        lnrho_ = var1
        ee_    = var2
        TT_    = (2.0/3.0)*TT_ion*(ee_/ee_ion-yH0)/(1+yH0+xHe-xH2)
        if (present(lnTT).or.present(ss)) then
          lnTT_  = log(TT_)
          if (present(ss)) ss_ = (lnTT_-(lnTTlnrho*lnrho_)-lnTT0)/lnTTss
        endif
        if (present(pp)) pp_ = (1+yH0+xHe-xH2)*exp(lnrho_)*TT_*ss_ion
!
      case (ilnrho_pp)
        lnrho_ = var1
        pp_    = var2
        TT_    = pp_/((1+yH0+xHe-xH2)*ss_ion*exp(lnrho_))
        if (present(lnTT).or.present(ss)) then
          lnTT_  = log(TT_)
          if (present(ss)) ss_ = (lnTT_-(lnTTlnrho*lnrho_)-lnTT0)/lnTTss
        endif
        if (present(ee)) ee_ = 1.5*(1+yH0+xHe-xH2)*ss_ion*TT_+yH0*ee_ion
!
      case (ilnrho_lnTT)
        lnrho_ = var1
        lnTT_  = var2
        if (present(ss)) ss_ = (lnTT_-lnTTlnrho*lnrho_-lnTT0)/lnTTss
        if (present(ee) .or. present(pp)) then
          TT_  = exp(lnTT_)
          if (present(ee)) ee_ = 1.5*(1+yH0+xHe-xH2)*ss_ion*TT_+yH0*ee_ion   !tbchecked
          if (present(pp)) pp_ = (1+yH0+xHe-xH2)*exp(lnrho_)*TT_*ss_ion      !tbchecked
        endif
!
     end select
!
      if (present(lnrho)) lnrho=lnrho_
      if (present(ss)) ss=ss_
      if (present(yH)) yH=yH0
      if (present(lnTT)) lnTT=lnTT_
      if (present(ee)) ee=ee_
      if (present(pp)) pp=pp_
      if (present(cs2)) cs2=impossible
!
    endsubroutine eoscalc_elem
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)

      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), intent(out), optional :: lnrho,ss
      real, dimension(nx), intent(out), optional :: yH,lnTT
      real, dimension(nx), intent(out), optional :: ee,pp,cs2

      if (.not.any((/ilnrho_ss, ilnrho_ee, ilnrho_pp, ilnrho_lnTT/)==ivars))  &
        call fatal_error('eoscalc_point','invalid combination of thermodynamic variables')

      ivars_mod=ivars

      call eoscalc_elem(var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)

    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine read_eos_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=eos_init_pars, IOSTAT=iostat)
!
    endsubroutine read_eos_init_pars
!***********************************************************************
    subroutine write_eos_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=eos_init_pars)
!
    endsubroutine write_eos_init_pars
!***********************************************************************
    subroutine read_eos_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=eos_run_pars, IOSTAT=iostat)
!
    endsubroutine read_eos_run_pars
!***********************************************************************
    subroutine write_eos_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=eos_run_pars)
!
    endsubroutine write_eos_run_pars
!***********************************************************************
    subroutine isothermal_entropy(lnrho,T0,ss)
!
!  Isothermal stratification: initilizes ss from lnrho.
!  This routine should be independent of the gravity module used.
!
!  Sound speed (and hence Temperature), is initialised to the reference value:
!           sound speed: cs^2_0            from start.in
!           density: rho0 = exp(lnrho0)
!
!  11-jun-03/tony: extracted from isothermal routine in Density module
!                  to allow isothermal condition for arbitrary density
!  17-oct-03/nils: works also with leos_ionization=T
!  18-oct-03/tobi: distributed across ionization modules
!
      real, dimension(mx,my,mz), intent(in) :: lnrho
      real, dimension(mx,my,mz), intent(out) :: ss
      real, intent(in) :: T0
!
      ss=ss_ion*((1+yH0+xHe-xH2)*(1.5*log(T0/TT_ion)-lnrho+2.5)-yH_term-one_yH_term-xHe_term)
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine isothermal_lnrho_ss(lnrho,T0,rho0,ss)
!
!  Isothermal stratification for lnrho and ss (for yH=0!)
!
!  Uses T=T0 everywhere in the box and rho=rho0 in the mid-plane
!
!  Currently only works with gravz_profile='linear', but can easily be
!  generalised.
!
!  11-feb-04/anders: Programmed more or less from scratch
!
      use Gravity, only: gravz_profile
!
      real, dimension(mx,my,mz), intent(out) :: lnrho, ss
      real, intent(in) :: T0,rho0
      real, dimension(nx) :: lnTT
!
      if (gravz_profile /= 'linear') call not_implemented &
          ('isothermal_lnrho_ss','for other than linear gravity profile')
!
!  First calculate hydrostatic density stratification when T=T0
!
      do n=n1,n2
        lnrho(l1:l2,m1:m2,n) = -(Omega*z(n))**2/(2*(1.+xHe-xH2)*ss_ion*T0)+log(rho0)
      enddo
!
!  Then calculate entropy as a function of T0 and lnrho
!
      do m=m1,m2
        do n=n1,n2
          lnTT=log(T0)
          call eoscalc_pencil(ilnrho_lnTT,lnrho(l1:l2,m,n),lnTT,ss=ss(l1:l2,m,n))
        enddo
      enddo
!
    endsubroutine isothermal_lnrho_ss
!***********************************************************************
    subroutine bc_ss_flux(f,topbot,lone_sided)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!   3-oct-16/MR: added new optional switch lone_sided
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      logical, optional :: lone_sided
!
      call not_implemented("bc_ss_flux","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_flux_turb(f,topbot)
!
!  dummy routine
!
!   4-may-2009/axel: dummy routine
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb
!***********************************************************************
    subroutine bc_ss_flux_turb_x(f,topbot)
!
!  dummy routine
!
!   31-may-2010/pete: dummy routine
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_x(f,topbot)
!
!   23-apr-2014/pete: dummy
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_mean_x(f,topbot)
!
!   07-jan-2015/pete: dummy
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_mean_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_z(f,topbot)
!
!   15-jul-2014/pete: dummy
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_z
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  23-jun-2003/tony: implemented for leos_fixed_ionization
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_temp_old","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_old
!***********************************************************************
    subroutine bc_ss_temp_x(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_temp_x","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_y(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_temp_y","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot,lone_sided)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      logical, optional :: lone_sided
!
      call not_implemented("bc_ss_temp_z","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_lnrho_temp_z(f,topbot)
!
!  boundary condition for density: constant temperature
!
!  19-aug-2005/tobi: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_lnrho_temp_z","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_temp_z
!***********************************************************************
    subroutine bc_lnrho_pressure_z(f,topbot)
!
!  boundary condition for density: constant pressure
!
!  19-aug-2005/tobi: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_lnrho_pressure_z","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_temp2_z","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_temp3_z(f,topbot)
!
!  31-jan-2013/axel: coded to impose cs2bot and dcs2bot at bottom
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_ss_temp3_z','in eos_fixed_ionization')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp3_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_stemp_x","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_x
!***********************************************************************
    subroutine bc_ss_stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_stemp_y","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_stemp_z","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_a2stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_a2stemp_x","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_x
!***********************************************************************
    subroutine bc_ss_a2stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_a2stemp_y","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_y
!***********************************************************************
    subroutine bc_ss_a2stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_a2stemp_z","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
!  boundary condition for entropy
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the entropy module
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
!  The 'ce' boundary condition for entropy makes the energy constant at
!  the boundaries.
!  This assumes that the density is already set (ie density must register
!  first!)
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_stellar_surface","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_lnrho_cfb_r_iso","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_lnrho_hds_z_iso","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hds_z_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_lnrho_hdss_z_iso","in eos_fixed_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine bc_ism(f,topbot,j)
!
!  30-nov-15/fred: Replaced bc_ctz and bc_cdz.
!  Apply observed scale height locally from Reynolds 1991, Manchester & Taylor
!  1981 for warm ionized gas - dominant scale height above 500 parsecs.
!  Apply constant local temperature across boundary for entropy.
!  Motivation to prevent numerical spikes in shock fronts, which cannot be
!  absorbed in only three ghost cells, but boundary thermodynamics still
!  responsive to interior dynamics.
!  06-jun-22/fred update to allow setting scale height in start.in or run.in
!  default is density_scale_factor=impossible so that scale_factor is 0.9, assuming
!  unit_length = 1 kpc and scale is 400 pc. To change scale height add to
!  start_pars or run_pars density_scale_factor=... in dimensionless units
!  Copied from eos_ionization written for entropy - may need revision
!  Currently not correct for energy variable
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,k
      real :: density_scale1, density_scale
!
      if (density_scale_factor==impossible) then
        density_scale=density_scale_cgs/unit_length
      else
        density_scale=density_scale_factor
      endif
      density_scale1=1./density_scale
!
      select case (topbot)
!
      case(BOT)               ! bottom boundary
        do k=1,nghost
          if (j==irho .or. j==ilnrho) then
            if (ldensity_nolog) then
              f(:,:,k,j)=f(:,:,n1,j)*exp(-(z(n1)-z(k))*density_scale1)
            else
              f(:,:,k,j)=f(:,:,n1,j) - (z(n1)-z(k))*density_scale1
            endif
          else if (j==iss) then
            f(:,:,k,j)=f(:,:,n1,j) + (z(n1)-z(k))*density_scale1
            !if (ldensity_nolog) then
            !  f(:,:,n1-k,j)=f(:,:,n1,j)+(cp-cv) * &
            !      (log(f(:,:,n1,j-1))-log(f(:,:,n1-k,j-1))) + &
            !      cv*log((z(n1)-z(n1-k))*density_scale+1.)
            !else
            !  f(:,:,n1-k,j)=f(:,:,n1,j)+(cp-cv)*&
            !      (f(:,:,n1,j-1)-f(:,:,n1-k,j-1))+&
            !      cv*log((z(n1)-z(n1-k))*density_scale+1.)
            !endif
          else
            call fatal_error('bc_ism','only for irho, ilnrho or iss')
          endif
        enddo
!
      case(TOP)               ! top boundary
        do k=1,nghost
          if (j==irho .or. j==ilnrho) then
            if (ldensity_nolog) then
              f(:,:,n2+k,j)=f(:,:,n2,j)*exp(-(z(n2+k)-z(n2))*density_scale1)
            else
              f(:,:,n2+k,j)=f(:,:,n2,j) - (z(n2+k)-z(n2))*density_scale1
            endif
          else if (j==iss) then
            f(:,:,n2+k,j)=f(:,:,n2,j) + (z(n2+k)-z(n2))*density_scale1
            !if (ldensity_nolog) then
            !  f(:,:,n2+k,j)=f(:,:,n2,j)+(cp-cv)*&
            !      (log(f(:,:,n2,j-1))-log(f(:,:,n2+k,j-1)))+&
            !      cv*log((z(n2+k)-z(n2))*density_scale+1.)
            !else
            !  f(:,:,n2+k,j)=f(:,:,n2,j)+(cp-cv)*&
            !      (f(:,:,n2,j-1)-f(:,:,n2+k,j-1))+&
            !      cv*log((z(n2+k)-z(n2))*density_scale+1.)
            !endif
          else
            call fatal_error('bc_ism','only for irho, ilnrho or iss')
          endif
        enddo
!
      case default
        call fatal_error('bc_ss_flux','topbot should be BOT or TOP')
!
      endselect
!
    endsubroutine bc_ism
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(cs20,p_par(1))

    endsubroutine pushpars2c
!***********************************************************************
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'eos_dummies.inc'
!***********************************************************************
endmodule EquationOfState
