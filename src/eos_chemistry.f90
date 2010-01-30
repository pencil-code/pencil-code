! $Id$
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
! PENCILS PROVIDED lnTT; cp1tilde; glnTT(3); TT; TT1; gTT(3)
! PENCILS PROVIDED hss(3,3); hlnTT(3,3); del2ss; del6ss; del2lnTT; del6lnTT
! PENCILS PROVIDED yH; ee; ss; pp; delta; ppvap; csvap2
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
  integer, parameter :: ilnrho_ss=1,ilnrho_ee=2,ilnrho_pp=3
  integer, parameter :: ilnrho_lnTT=4,ilnrho_cs2=5
  integer, parameter :: irho_cs2=6, irho_ss=7, irho_lnTT=8, ilnrho_TT=9
  integer, parameter :: ipp_ss=11,ipp_cs2=12
!
  integer :: iglobal_cs2, iglobal_glnTT
!
  real :: lnTT0=impossible
!
  real :: mu=1.
  real :: cs0=1., rho0=1.
  real :: cs20=1., lnrho0=0.
  real :: ptlaw=3./4.
  real :: gamma=5./3.
  real :: Rgas_cgs=0., Rgas, Rgas_unit_sys=1.,  error_cp=1e-6
  real :: gamma_m1    !(=gamma-1)
  real :: gamma_inv   !(=1/gamma)
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
  character (len=20) :: input_file='chem.inp'
  logical, SAVE ::  lcheminp_eos=.false.
  logical :: l_gamma_m1=.false.
  logical :: l_gamma=.false.
  logical :: l_cp=.false.
!
  character (len=labellen) :: ieos_profile='nothing'
  real, dimension(mz) :: profz_eos=1.
!
  namelist /eos_init_pars/  mu, cp, cs0, rho0, gamma, error_cp, ptlaw
!
  namelist /eos_run_pars/   mu, cp, cs0, rho0, gamma, error_cp, ptlaw
!
  contains
!***********************************************************************
    subroutine register_eos()
!
!  14-jun-03/axel: adapted from register_eos
!
      use Sub
!
      leos=.true.
      leos_chemistry=.true.
!
      ilnTT = 0
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
!
      use Mpicomm, only: stop_it
!
      real ::  cp_reference

!
!  set gamma_m1, cs20, and lnrho0
!  (used currently for non-dimensional equation of state)
!
      gamma_m1=gamma-1.
      gamma_inv=1./gamma
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
!  When gamma=1. (gamma_m1=0.), write Rgas=mu*cp or cp=Rgas/mu.
!
      if (unit_system == 'cgs') then
         Rgas_unit_sys = k_B_cgs/m_u_cgs
      elseif (unit_system == 'SI') then
         Rgas_unit_sys = k_B_cgs/m_u_cgs*1.e-4
      endif
!
      if (unit_temperature == impossible) then
    !    if (cp == impossible) cp=1.
    !    if (gamma_m1 == 0.) then
    !      Rgas=mu*cp
    !    else
    !      Rgas=mu*gamma_m1*gamma_inv*cp
    !    endif
    !    unit_temperature=unit_velocity**2*Rgas/Rgas_unit_sys
        call stop_it('unit_temperature is not found!')
      else
        Rgas=Rgas_unit_sys*unit_temperature/unit_velocity**2
        if (cp == impossible) then
          if (gamma_m1 == 0.) then
            cp=Rgas/mu
          else
            cp=Rgas/(mu*gamma_m1*gamma_inv)
          endif
        else
!
!  checking whether the units are overdetermined.
!  This is assumed to be the case when the to differ by error_cp
!
          if (gamma_m1 == 0.) then
            cp_reference=Rgas/mu
          else
            cp_reference=Rgas/(mu*gamma_m1*gamma_inv)
          endif
          if (abs(cp-cp_reference)/cp > error_cp) then
            if (lroot) print*,'initialize_eos: consistency: cp=',cp, &
               'while: cp_reference=',cp_reference
            call stop_it('initialize_eos')
          endif
        endif
      endif
      cp1=1./cp
      cv=gamma_inv*cp
      cv1=gamma*cp1
!
!  Need to calculate the equivalent of cs0
!  Distinguish between gamma=1 case and not.
!
      if (gamma_m1 /= 0.) then
        lnTT0=log(cs20/(cp*gamma_m1))  !(general case)
      else
        lnTT0=log(cs20/cp)  !(isothermal/polytropic cases: check!)
      endif


      inquire(FILE=input_file, EXIST=lcheminp_eos)

      if (lcheminp_eos) then
       l_gamma_m1=.true.
       l_gamma=.true.
       l_cp=.true.
      endif
!
!  check that everything is OK
!
      if (lroot) then

       if (.not. l_gamma_m1) then
        print*,'initialize_eos: unit_temperature=',unit_temperature
        print*,'initialize_eos: cp,lnTT0,cs0=',cp,lnTT0,cs0
       else
        print*,'initialize_eos: chem.imp is found! Now cp, cv, gamma, mu are pencils ONLY!'
       endif  
      endif
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos()
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
        mu_tmp=1.!+2.97153*xHe
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
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_eos
!***********************************************************************
!    subroutine get_slices_eos(f,slices)
!
!      real, dimension (mx,my,mz,mfarray) :: f
!      type (slice_data) :: slices
!
!      call keep_compiler_quiet(f)
!      call keep_compiler_quiet(slices%ready)
!
!    endsubroutine get_slices_eos
!***********************************************************************
!***********************************************************************
    subroutine get_slices_eos(f,slices)
!
!  Write slices for animation of Eos variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices.
!
      select case (trim(slices%name))
!
!  Temperature.
!
        case ('lnTT')
          slices%yz =f(ix_loc,m1:m2,n1:n2,ilnTT)
          slices%xz =f(l1:l2,iy_loc,n1:n2,ilnTT)
          slices%xy =f(l1:l2,m1:m2,iz_loc,ilnTT)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,ilnTT)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,ilnTT)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,ilnTT)
          slices%ready=.true.
     
      endselect
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

    lpenc_requested(i_lnTT)=.true.
    lpenc_requested(i_TT)=.true.
    lpenc_requested(i_TT1)=.true.
    lpenc_requested(i_glnTT)=.true.
    lpenc_requested(i_del2lnTT)=.true.


    endsubroutine pencil_criteria_eos
!***********************************************************************
    subroutine pencil_interdep_eos(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
! Modified by Natalia. Taken from  eos_temperature_ionization module

   logical, dimension(npencils) :: lpencil_in

      if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
      if (lpencil_in(i_TT1)) lpencil_in(i_TT)=.true.

      call keep_compiler_quiet(lpencil_in)
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
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
! THE FOLLOWING 2 ARE CONCEPTUALLY WRONG
! FOR pretend_lnTT since iss actually contain lnTT NOT entropy!
! The code is not wrong however since this is correctly
! handled by the eos moduleinteger :: tm1=1, tm2=2,tm3=3.

!Natalia: removed all previous staff and included my own


!
!  Temperature
!
       if (lpencil(i_lnTT)) then
         if (ltemperature_nolog) then
          p%lnTT=log(f(l1:l2,m,n,iTT))
         else
          p%lnTT=f(l1:l2,m,n,ilnTT)
         endif
           
       endif
       if (lpencil(i_TT))  then
         if (ltemperature_nolog) then
           p%TT=f(l1:l2,m,n,iTT)
         else
           p%TT=exp(p%lnTT)
         endif
       endif
       if (lpencil(i_TT1)) p%TT1=1./p%TT!
!
!  Temperature laplacian and gradient
!
        if (lpencil(i_glnTT)) then
         if (ltemperature_nolog) then
           call grad(f,iTT,p%glnTT)
           p%glnTT(:,1)=p%glnTT(:,1)/p%TT(:)
           p%glnTT(:,2)=p%glnTT(:,2)/p%TT(:)
           p%glnTT(:,3)=p%glnTT(:,3)/p%TT(:)
         else
           call grad(f,ilnTT,p%glnTT)
         endif
        endif
        
        if (ltemperature_nolog) then
         if (lpencil(i_gTT)) call grad(f,iTT,p%gTT)
         if (lpencil(i_del2lnTT)) call del2(f,iTT,p%del2lnTT)
        else
         if (lpencil(i_del2lnTT)) call del2(f,ilnTT,p%del2lnTT)
        endif


 call keep_compiler_quiet(f)
 call keep_compiler_quiet(p)

!  Natalia: 26.02.2008: calculation of additional penciles

    endsubroutine calc_pencils_eos
!***********************************************************************
   subroutine ioninit(f)
!
!  the ionization fraction has to be set to a value yH0 < yH < yHmax before
!  rtsafe is called for the first time
!
!  12-jul-03/tobi: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
!   calculate degree of ionization and temperature
!   This routine is called from equ.f90 and operates on the full 3-D array.
!
!   13-jun-03/tobi: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine ioncalc
!***********************************************************************

   subroutine getdensity(EE,TT,yH,rho)

     real, intent(in) :: EE,TT,yH
     real, intent(inout) :: rho

     rho = EE * cv1 / TT
      call keep_compiler_quiet(yH)

   endsubroutine getdensity
!***********************************************************************
    subroutine get_cp1(cp1_)
!
!  04-nov-06/axel: added to alleviate spurious use of pressure_gradient
!
!  return the value of cp1 to outside modules
!
    use Mpicomm, only: stop_it 

      real, intent(out) :: cp1_
!
    if (.not. l_cp) then
      cp1_=cp1
    else
      call stop_it('chem.inp is found: pressure_gradient_point can not be used for this moment')
    endif
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
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: cs2,cp1tilde
!
      cs2=impossible
      cp1tilde=impossible
!
      call keep_compiler_quiet(f)
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
      use Mpicomm, only: stop_it 
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: cs2,cp1tilde
!
      call stop_it(' pressure_gradient_point should never be called')
!
      call keep_compiler_quiet(cs2,cp1tilde,ss,lnrho)
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
      use Mpicomm, only: stop_it 
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(in) :: glnrho,gss
      real, dimension(nx,3), intent(out) :: glnTT
!
      call stop_it('temperature_gradient should never be called')
!
!  given that we just stopped, it cannot become worse by setting
!  cp1tilde to impossible, which allows the compiler to compile.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(glnrho,gss,glnTT)
!
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_laplacian(f,del2lnrho,del2ss,del2lnTT)
!
!   Calculate thermodynamical quantities
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(in) :: del2lnrho,del2ss
      real, dimension(nx), intent(out) :: del2lnTT
  !    type (pencil_case) :: p

!
      call fatal_error('temperature_laplacian','SHould not be called!')
!
     call keep_compiler_quiet(f)
     call keep_compiler_quiet(del2lnrho,del2ss,del2lnTT)
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

    !  type (pencil_case) :: p
!
      call fatal_error('temperature_hessian','This routine is not coded for eos_chemistry')
!
      hlnTT(:,:,:)=0
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(hss)
      call keep_compiler_quiet(hlnrho)
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
      elseif (psize==mx) then
        lnrho_=f(:,m,n,ilnrho)
      else
        call not_implemented("eosperturb")
      endif
!
      call keep_compiler_quiet(present(ee))
      call keep_compiler_quiet(present(pp))
      call keep_compiler_quiet(present(ss))
!
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
      use Diagnostics
     
     ! type (pencil_case) :: p
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
           if (.not. l_gamma_m1) then
            ss_=-cv*gamma_m1*(lnrho_-lnrho0)
           else
      !      ss_=-p%cv*p%gamma_m1*(lnrho_-lnrho0)
           endif
          else
            ss_=f(l1:l2,m,n,ieosvar2)
          endif
        case (mx)
          lnrho_=f(:,m,n,ieosvar1)
          if (leos_isentropic) then
            ss_=0
          elseif (leos_isothermal) then
           if (.not. l_gamma_m1) then
            ss_=-cv*gamma_m1*(lnrho_-lnrho0)
           else
       !     ss_=-p%cv*p%gamma_m1*(lnrho_-lnrho0)
           endif
          else
            ss_=f(:,m,n,ieosvar2)
          endif
        case default
          call fatal_error('eoscalc_farray','no such pencil size')
        end select

        if (.not. l_gamma_m1) then
          lnTT_=lnTT0+cv1*ss_+gamma_m1*(lnrho_-lnrho0)
        else
        !  lnTT_=lnTT0+p%cv1*ss_+p%gamma_m1*(lnrho_-lnrho0)
        endif
        if (gamma_m1==0.) &
            call fatal_error('eoscalc_farray','gamma=1 not allowed w/entropy')
        if (present(lnrho)) lnrho=lnrho_
        if (present(lnTT)) lnTT=lnTT_
       if (.not.  l_gamma_m1) then 
        if (present(ee)) ee=cv*exp(lnTT_)
        if (present(pp)) pp=(cp-cv)*exp(lnTT_+lnrho_)
       else
       ! if (present(ee)) ee=p%cv*exp(lnTT_)
       ! if (present(pp)) pp=(p%cp-p%cv)*exp(lnTT_+lnrho_)
       endif  
!
! Log rho and Log T
!
      case (ilnrho_lnTT,irho_lnTT)
        select case (psize)
        case (nx)
          lnrho_=f(l1:l2,m,n,ieosvar1)
          if (leos_isentropic) then
           if (.not. l_cp) then
            lnTT_=lnTT0+(cp-cv)*(lnrho_-lnrho0)
           else
        !    lnTT_=lnTT0+(p%cp-p%cv)*(lnrho_-lnrho0)
           endif
          elseif (leos_isothermal) then
            lnTT_=lnTT0
          else
            lnTT_=f(l1:l2,m,n,ieosvar2)
          endif
        case (mx)
          lnrho_=f(:,m,n,ieosvar1)
          if (leos_isentropic) then
           if (.not. l_cp) then
            lnTT_=lnTT0+(cp-cv)*(lnrho_-lnrho0)
           else
         !   lnTT_=lnTT0+(p%cp-p%cv)*(lnrho_-lnrho0)
           endif
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

       if (.not. l_cp) then
        if (present(ee)) ee=cv*exp(lnTT_)
        if (ieosvars==ilnrho_lnTT) then
          if (present(pp)) pp=(cp-cv)*exp(lnTT_+lnrho_)
        else
          if (present(pp)) pp=(cp-cv)*exp(lnTT_)*lnrho_
        endif
       else
      !  if (present(ee)) ee=p%cv*exp(lnTT_)
        if (ieosvars==ilnrho_lnTT) then
      !    if (present(pp)) pp=(p%cp-p%cv)*exp(lnTT_+lnrho_)
        else
       !   if (present(pp)) pp=(p%cp-p%cv)*exp(lnTT_)*lnrho_
        endif
       endif
!
! Log rho and T
!
      case (ilnrho_TT)
          call fatal_error('eoscalc_farray','no implemented for lnrho_TT')
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
           if (.not. l_gamma_m1) then 
            cs2_=exp(gamma_m1*(lnrho_-lnrho0)+log(cs20))
           else
        !    cs2_=exp(p%gamma_m1*(lnrho_-lnrho0)+log(cs20))
           endif  
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
           if (.not. l_gamma_m1) then 
            cs2_=exp(gamma_m1*(lnrho_-lnrho0)+log(cs20))
           else
         !   cs2_=exp(p%gamma_m1*(lnrho_-lnrho0)+log(cs20))
           endif
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
       if (.not. l_gamma_m1) then 
        if (present(ee)) ee=gamma_inv*cs2_/gamma_m1
        if (present(pp)) pp=gamma_inv*cs2_*exp(lnrho_)
       else
      !  if (present(ee)) ee=p%gamma_inv*cs2_/p%gamma_m1
      !  if (present(pp)) pp=p%gamma_inv*cs2_*exp(lnrho_)
       endif
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
      use Mpicomm, only: stop_it 
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp,cs2
      real :: lnrho_,ss_,lnTT_,ee_,pp_,cs2_,TT_
!
      if (gamma_m1==0.) call fatal_error('eoscalc_point','gamma=1 not allowed w/entropy')
!

      if (l_gamma) call stop_it('now gamma is a pencil: eoscalc_point can not be used for this moment')
 
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
        cs2_=gamma_m1*TT_

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

    !  type (pencil_case) :: p
!
      if (gamma_m1==0.) call fatal_error('eoscalc_pencil','gamma=1 not allowed w/entropy')
!
      select case (ivars)

      case (ilnrho_ss)
        lnrho_=var1
        ss_=var2
       if (.not. l_gamma_m1) then
        lnTT_=lnTT0+cv1*ss_+gamma_m1*(lnrho_-lnrho0)
        ee_=cv*exp(lnTT_)
        pp_=(cp-cv)*exp(lnTT_+lnrho_)
        cs2_=gamma*gamma_m1*ee_
        cs2_=cs20*cv1*ee_
       else
   !     lnTT_=lnTT0+p%cv1*ss_+p%gamma_m1*(lnrho_-lnrho0)
     !   ee_=p%cv*exp(lnTT_)
     !   pp_=(p%cp-p%cv)*exp(lnTT_+lnrho_)
     !  cs2_=p%gamma*p%gamma_m1*ee_
    !    cs2_=cs20*p%cv1*ee_
       endif

      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
       if (.not. l_gamma_m1) then
        lnTT_=log(cv1*ee_)
        ss_=cv*(lnTT_-lnTT0-gamma_m1*(lnrho_-lnrho0))
        pp_=gamma_m1*ee_*exp(lnrho_)
        cs2_=gamma*gamma_m1*ee_
       else
      !  lnTT_=log(p%cv1*ee_)
      !  ss_=p%cv*(lnTT_-lnTT0-p%gamma_m1*(lnrho_-lnrho0))
      !  pp_=p%gamma_m1*ee_*exp(lnrho_)
      !  cs2_=p%gamma*p%gamma_m1*ee_
       endif

      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
       if (.not. l_gamma_m1) then
        ss_=cv*(log(pp_*exp(-lnrho_)*gamma/cs20)-gamma_m1*(lnrho_-lnrho0))
        ee_=pp_*exp(-lnrho_)/gamma_m1
        lnTT_=log(cv1*ee_)
        cs2_=gamma*gamma_m1*ee_
       else
       ! ss_=p%cv*(log(pp_*exp(-lnrho_)*p%gamma/cs20)-p%gamma_m1*(lnrho_-lnrho0))
       ! ee_=pp_*exp(-lnrho_)/p%gamma_m1
        !lnTT_=log(p%cv1*ee_)
       ! cs2_=p%gamma*p%gamma_m1*ee_
       endif

      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
       if (.not. l_gamma_m1) then
        ss_=cv*(lnTT_-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*exp(lnTT_)
        pp_=ee_*exp(lnrho_)*gamma_m1
        cs2_=gamma*gamma_m1*ee_
       else
       ! ss_=p%cv*(lnTT_-lnTT0-p%gamma_m1*(lnrho_-lnrho0))
       ! ee_=p%cv*exp(lnTT_)
       ! pp_=ee_*exp(lnrho_)*p%gamma_m1
       ! cs2_=p%gamma*p%gamma_m1*ee_
       endif 

      case (ilnrho_TT)
        lnrho_=var1
        TT_=var2
       if (.not. l_gamma_m1) then
        ss_=cv*(log(TT_)-lnTT0-gamma_m1*(lnrho_-lnrho0))
        ee_=cv*TT_
        pp_=ee_*exp(lnrho_)*gamma_m1
        cs2_=gamma_m1*TT_
       else
       ! ss_=p%cv*(log(TT_)-lnTT0-p%gamma_m1*(lnrho_-lnrho0))
       ! ee_=p%cv*TT_
       ! pp_=ee_*exp(lnrho_)*p%gamma_m1
       ! cs2_=p%gamma_m1*TT_
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
    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine get_soundspeed(lnTT,cs2)
!
!  Calculate sound speed for given temperature
!
!  20-Oct-03/tobi: Coded
!
     use Mpicomm, only: stop_it 

      real, intent(in)  :: lnTT
      real, intent(out) :: cs2
!
     if (.not.l_gamma_m1) then 
      cs2=gamma_m1*cp*exp(lnTT)
     else
        call stop_it('chem.inp is found: get_soundspeed can not be used for this moment')
     endif 
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(T0)
      call keep_compiler_quiet(rho0)
!
    endsubroutine isothermal_lnrho_ss
!***********************************************************************
 !   subroutine Hminus_opacity(f,kapparho)
!
!  dummy routine
!
!  03-apr-2004/tobi: coded
!
 !     real, dimension(mx,my,mz,mfarray), intent(in) :: f
 !     real, dimension(mx,my,mz), intent(out) :: kapparho

 !     call fatal_error('Hminus_opacity',"opacity_type='Hminus' may not be used with noionization")

 !   endsubroutine Hminus_opacity

!***********************************************************************
     subroutine get_average_pressure(average_density,average_pressure)

!   01-dec-2009/piyali+dhrube: coded
      use Cdata
!      
      real, intent(in):: average_density
      real, intent(out):: average_pressure
      call keep_compiler_quiet(average_density)
      call keep_compiler_quiet(average_pressure)
    endsubroutine get_average_pressure
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
      use Gravity
      use SharedVariables,only:get_shared_variable
      use Mpicomm, only:stop_it
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      logical, pointer :: lmultilayer, lheatc_chiconst
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy
      integer :: i,ierr
!

    if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_flux_turb(f,topbot)
!
!  dummy routine
!
!   4-may-2009/axel: dummy routine
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i
!

     if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: i

      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

!
      if (ldebug) print*,'bc_ss_temp_x: cs20,cs0=',cs20,cs0
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
            tmp = 2*cv*log(cs2top/cs20)
            f(l2,:,:,iss) = 0.5*tmp - (cp-cv)*(f(l2,:,:,ilnrho)-lnrho0)
            do i=1,nghost
               f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
                    - (cp-cv)*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
            enddo
        elseif (lentropy .and. pretend_lnTT) then
           f(l2,:,:,iss) = log(cs2top/gamma_m1)
           do i=1,nghost; f(l2+i,:,:,iss)=2*f(l2,:,:,iss)-f(l2-i,:,:,iss); enddo
        elseif (ltemperature) then
           f(l2,:,:,ilnTT) = log(cs2top/gamma_m1)
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: i
!

      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: i
!

      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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
           tmp = 2*cv*log(cs2bot/cs20)
           f(:,:,n1,iss) = 0.5*tmp - (cp-cv)*(f(:,:,n1,ilnrho)-lnrho0)
           do i=1,nghost
              f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
                   - (cp-cv)*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
           enddo
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
        if (lentropy .and. .not. pretend_lnTT) then 
           tmp = 2*cv*log(cs2top/cs20)
           f(:,:,n2,iss) = 0.5*tmp - (cp-cv)*(f(:,:,n2,ilnrho)-lnrho0)
           do i=1,nghost
              f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
                   - (cp-cv)*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
           enddo
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: i
!

      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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

      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp
      integer :: i
!

      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!

      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
! 
      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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

    if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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

    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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
      use Gravity
      use Sub, only: div

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
      real, dimension (my,mz) :: cs2,gravterm,centterm,uphi
      real :: potp,potm,rad,step
      integer :: i


      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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
            f(l1-i,:,:,irho)  =f(l1+i,:,:,irho)*exp(gravterm + centterm)
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
            f(l2+i,:,:,irho)   = f(l2-i,:,:,irho)*exp(gravterm + centterm)
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
!
!  12-Juil-2006/dintrans: coded
!
      use Gravity
      use Sub, only: div

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
      real, dimension (mx,my) :: cs2
      real, dimension (nx) :: shock,divu
      real :: dlnrhodz, dssdz
      real :: potp,potm
      integer :: i


     if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

      select case (topbot)

!
!  Bottom boundary
!
      case ('bot')

        if (lentropy) then

          if (bcz1(iss)/='hs') then
            call fatal_error("bc_lnrho_hydrostatic_z", &
                             "This boundary condition for density is "// &
                             "currently only correct for bcz1(iss)='hs'")
          endif

          dlnrhodz = gamma*gravz/cs2bot
          dssdz = -gamma_m1*gravz/cs2bot
          do i=1,nghost
            f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) - 2*i*dz*dlnrhodz
            f(:,:,n1-i,iss   ) = f(:,:,n1+i,iss   ) - 2*i*dz*dssdz
          enddo

        else

          do i=1,nghost
            call potential(z=z(n1-i),pot=potm)
            call potential(z=z(n1+i),pot=potp)
            cs2 = cs2bot
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
            if (ldensity_nolog) then
              f(:,:,n1-i,irho)   = f(:,:,n1+i,irho)*exp(-(potm-potp)/cs2)
            else
              f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) - (potm-potp)/cs2
            endif
          enddo

        endif

!
!  Top boundary
!
      case ('top')

        if (lentropy) then

          if (bcz2(iss)/='hs') then
            call fatal_error("bc_lnrho_hydrostatic_z", &
                             "This boundary condition for density is "//&
                             "currently only correct for bcz2(iss)='hs'")
          endif

          dlnrhodz = gamma*gravz/cs2top
          dssdz = -gamma_m1*gravz/cs2top
          do i=1,nghost
            f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) + 2*i*dz*dlnrhodz
            f(:,:,n2+i,iss   ) = f(:,:,n2-i,iss   ) + 2*i*dz*dssdz
          enddo

        else

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

        endif

      case default

      endselect

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

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot

      real, dimension (nx,ny) :: kx,ky,kappa,exp_fact
      real, dimension (nx,ny) :: tmp_re,tmp_im
      real :: pot
      integer :: i

      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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

        enddo
!
!  Potential field condition at the top
!
      case ('top')

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
      use Fourier, only: fourier_transform_xy_xy, fourier_transform_other
      use Gravity, only: potential
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
!
      real, dimension (nx,ny) :: kx,ky,kappa,exp_fact
      real, dimension (nx,ny) :: tmp_re,tmp_im
      real, dimension (nx) :: pot,rr_cyl,rr_sph,cs2,tmp1,tmp2
      integer :: i,mm_noghost

      if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'
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
      select case (topbot)
!
!  Potential field condition at the bottom
!
      case ('bot')

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
              tmp_re(:,mm_noghost) = f(l1:l2,m,n1+i,irho)*exp(+pot/cs2)
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
              f(l1:l2,m,n1-i,irho)   = tmp_re(:,mm_noghost)*exp(-pot/cs2)
            else
              f(l1:l2,m,n1-i,ilnrho) = tmp_re(:,mm_noghost) - pot/cs2
            endif
          enddo

        enddo
!
!  Potential field condition at the top
!
      case ('top')

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
              tmp_re(:,mm_noghost) = f(l1:l2,m,n2-i,irho)*exp(+pot/cs2)
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
              f(l1:l2,m,n2+i,irho)   = tmp_re(:,mm_noghost)*exp(-pot/cs2)
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
      use Gravity

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (nx) :: potm,potp,tmp1,tmp2,rr_cyl,rr_sph,cs2
      character (len=3), intent (in) :: topbot
      integer :: i
 

     if (l_gamma) print*,'gamma is a pencil now! Be careful with using such boundary conditions!!!'

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
              f(l1:l2,m,n1-i,irho)   = f(l1:l2,m,n1+i,irho)*exp((potm-potp)/cs2)
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
              f(l1:l2,m,n2+i,irho)   = f(l1:l2,m,n2-i,irho)*exp(-(potp-potm)/cs2)
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
