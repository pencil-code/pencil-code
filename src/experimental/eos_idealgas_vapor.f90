! $Id: eos_idealgas_vapor.f90,v 1.7 2010/02/03 10:30:37 ajohan Exp $
!
!  Equation of state for an ideal gas with variable water vapour.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .true.
! CPARAM logical, parameter :: leos_ionization = .false., leos_temperature_ionization=.false.
! CPARAM logical, parameter :: leos_idealgas = .false., leos_chemistry = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 3
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cp; cp1; cp1tilde
! PENCILS PROVIDED glnTT(3); TT; TT1; gTT(3); yH; hss(3,3); hlnTT(3,3)
! PENCILS PROVIDED del2ss; del6ss; del2lnTT; cv; cv1; del6lnTT; gamma
! PENCILS PROVIDED del2TT; del6TT; mumol; mumol1; glnmumol(3)
! PENCILS PROVIDED rho_anel; ppvap; csvap2; fvap; gfvap(3)
!
!***************************************************************
module EquationOfState
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../eos.h'
  include '../eos_params.h'
!
  integer :: iglobal_cs2, iglobal_glnTT
  real, dimension(nchemspec,18) :: species_constants
  real :: Rgas_unit_sys=1.0
  real :: lnTT0=impossible
  real :: mudry=1.0, muvap=1.0, mudry1=1.0, muvap1=1.0
  real :: cs0=1.0, rho0=1.0, pp0=1.0
  real :: cs20=1.0, lnrho0=0.0
  real :: ptlaw=0.0
  real :: gamma=5.0/3.0
  real :: Rgas_cgs=0.0, Rgas, error_cp=1.0e-6
  real :: gamma_m1    !(=gamma-1)
  real :: gamma1   !(=1/gamma)
  real :: cpdry=impossible, cpdry1=impossible
  real :: cvdry=impossible, cvdry1=impossible
  real :: cs2bot=1.0, cs2top=1.0
  integer :: imass=1
  integer :: ieosvars=-1, ieosvar1=-1, ieosvar2=-1, ieosvar_count=0
  logical :: leos_isothermal=.false., leos_isentropic=.false.
  logical :: leos_isochoric=.false., leos_isobaric=.false.
  logical :: leos_localisothermal=.false.
!
! Kishore: I have not checked what these are used for; just copied from eos_idealgas to get this module to compile.
!
  real :: Cp_const=impossible
  real :: Pr_number=0.7
  logical :: lpres_grad=.false.
!
!  Shared variables
!
  real :: fac_cs=1.0
  integer :: isothmid=0
!
! Indices for aux variables
!
  integer :: ifvap=0, imumol1=0
!
!  Input parameters.
!
  namelist /eos_init_pars/ &
      mudry, muvap, cpdry, cs0, rho0, gamma, error_cp, ptlaw
!
!  Run parameters.
!
  namelist /eos_run_pars/ &
      mudry, muvap, cpdry, cs0, rho0, gamma, error_cp, ptlaw
!
  contains
!***********************************************************************
    subroutine register_eos()
!
!  Register variables from the EquationOfState module.
!
!  06-jan-10/anders: adapted from eos_idealgas
!  04-dec-2024/Kishore: added shared variables (copied from eos_idealgas)
!
      use SharedVariables, only: put_shared_variable
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
          '$Id: eos_idealgas_vapor.f90,v 1.7 2010/02/03 10:30:37 ajohan Exp $')
!
! Shared variables
!
      call put_shared_variable('cs20',cs20,caller='register_eos')
      call put_shared_variable('isothmid',isothmid)
      call put_shared_variable('fac_cs',fac_cs)
!
    endsubroutine register_eos
!***********************************************************************
    subroutine units_eos()
!
!  This routine calculates things related to units and must be called
!  before the rest of the units are being calculated.
!
!  06-jan-10/anders: adapted from eos_idealgas
!
      real :: cp_reference
!
!  Set gamma_m1, cs20, and lnrho0
!  (used currently for non-dimensional equation of state)
!
      gamma_m1=gamma-1.0
      gamma1=1/gamma
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
        if (cpdry==impossible) cpdry=1.0
        if (gamma_m1==0.0) then
          Rgas=mudry*cpdry
        else
          Rgas=mudry*(1.0-gamma1)*cpdry
        endif
        unit_temperature=unit_velocity**2*Rgas/Rgas_unit_sys
      else
        Rgas=Rgas_unit_sys*unit_temperature/unit_velocity**2
        if (cpdry==impossible) then
          if (gamma_m1==0.0) then
            cpdry=Rgas/mudry
          else
            cpdry=Rgas/(mudry*gamma_m1*gamma1)
            if (headt) print*,'units_eos: cpdry=',cpdry
          endif
        else
!
!  Checking whether the units are overdetermined.
!  This is assumed to be the case when the to differ by error_cp.
!
          if (gamma_m1==0.0) then
            cp_reference=Rgas/mudry
          else
            cp_reference=Rgas/(mudry*gamma_m1*gamma1)
          endif
          if (abs(cpdry-cp_reference)/cpdry > error_cp) then
            if (lroot) print*,'initialize_eos: consistency: cpdry=', cpdry , &
                'while: cp_reference=', cp_reference
            call fatal_error('initialize_eos','')
          endif
        endif
      endif
      cpdry1=1/cpdry
      cvdry=gamma1*cpdry
      cvdry1=gamma*cpdry1
!
!  Need to calculate the equivalent of cs0.
!  Distinguish between gamma=1 case and not.
!
      if (gamma_m1/=0.0) then
        lnTT0=log(cs20/(cpdry*gamma_m1))  !(general case)
      else
        lnTT0=log(cs20/cpdry)  !(isothermal case)
      endif
      pp0=Rgas*exp(lnTT0)*rho0
!
!  Check that everything is OK.
!
      if (lroot) then
        print*, 'initialize_eos: unit_temperature=', unit_temperature
        print*, 'initialize_eos: cpdry, lnTT0, cs0, pp0=', &
            cpdry, lnTT0, cs0, pp0
      endif
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos(f)
!
!  Perform any post-parameter-read initialization
!
!  06-jan-10/anders: adapted from eos_idealgas
!
      use Sub, only: register_report_aux
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      mudry1=1/mudry
      muvap1=1/muvap
!
!  Initialize variable selection code (needed for RELOADing).
!
      ieosvars=-1
      ieosvar_count=0
!
!  fvap, mumol1, and cp as auxiliary variables.
!
      call register_report_aux('fvap',ifvap)
      call register_report_aux('mumol1',imumol1)
      call register_report_aux('cp',icp)
!
!  Write constants to disk. In future we may want to deal with this
!  using an include file or another module.
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,'(a,1pd26.16)') 'k_B=', k_B
        write (1,'(a,1pd26.16)') 'm_H=', m_H
        write (1,*) 'lnrho0=', lnrho0
        write (1,*) 'lnTTO=', lnTT0
        write (1,*) 'cs20=', cs20
        write (1,*) 'cpdry=', cpdry
        close (1)
      endif
!
    endsubroutine initialize_eos
!***********************************************************************
    subroutine select_eos_variable(variable,findex)
!
!  Select eos variable.
!
!  06-jan-10/anders: adapted from eos_idealgas
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
!  06-jan-10/anders: adapted from eos_idealgas
!
      real, dimension (mx,my,mz,mfarray), optional :: f
      real, intent(out) :: mu_tmp
!
      call fatal_error('getmu','not implemented')
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz), intent(out) :: mu1_full_tmp
!
      call fatal_error('getmu_array','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(mu1_full_tmp)
!
    endsubroutine getmu_array
!***********************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_eos
!***********************************************************************
    subroutine get_slices_eos(f,slices)
!
!  06-jan-10/anders: adapted from eos_idealgas
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_eos
!***********************************************************************
    subroutine pencil_interdep_eos(lpencil_in)
!
!  Interdependency among pencils from the EquationOfState module is specified
!  here.
!
!  06-jan-10/anders: adapted from eos_idealgas
!
      logical, dimension(npencils) :: lpencil_in
!
      if (leos_isentropic.or.leos_isothermal.or.leos_localisothermal) &
          call fatal_error('pencil_interdep_eos','leos_... case not implemented')
!
      if (lpencil_in(i_glnmumol)) then
        lpencil_in(i_mumol)=.true.
        !lpencil_in(i_gcc)=.true.
        !lpencil_in(i_gssat)=.true.
        lpencil_in(i_gfvap)=.true.
        lpencil_in(i_fvap)=.true.
      endif
      if (lpencil_in(i_mumol)) lpencil_in(i_mumol1)=.true.
      !if (lpencil_in(i_mumol1)) lpencil_in(i_cc)=.true.
      if (lpencil_in(i_mumol1)) lpencil_in(i_fvap)=.true.
      if (lpencil_in(i_cv1)) lpencil_in(i_cp1)=.true.
      if (lpencil_in(i_cp1)) lpencil_in(i_cp)=.true.
      if (lpencil_in(i_cv)) lpencil_in(i_cp)=.true.
      if (lpencil_in(i_cp)) lpencil_in(i_mumol1)=.true.
      if (lpencil_in(i_cs2)) then
        lpencil_in(i_cp)=.true.
        lpencil_in(i_TT)=.true.
      endif
      if (lpencil_in(i_pp)) then
        lpencil_in(i_cp)=.true.
        lpencil_in(i_cv)=.true.
        lpencil_in(i_rho)=.true.
        lpencil_in(i_TT)=.true.
      endif
      if (lpencil_in(i_ppvap)) then
        lpencil_in(i_fvap)=.true.
        lpencil_in(i_rho)=.true.
        lpencil_in(i_TT)=.true.
      endif
      if (lpencil_in(i_ee)) then
        lpencil_in(i_cv)=.true.
        lpencil_in(i_lnTT)=.true.
      endif
      if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
      if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
!
!  Pencils that depend on the chosen thermodynamical variables.
!
      select case (ieosvars)
!
!  Pencils for thermodynamic quantities for given lnrho or rho and ss.
!
      case (ilnrho_ss,irho_ss)
        if (lpencil_in(i_lnTT)) then
          lpencil_in(i_cv1)=.true.
          lpencil_in(i_ss)=.true.
          lpencil_in(i_lnrho)=.true.
        endif
        if (lpencil_in(i_glnTT)) then
          lpencil_in(i_glnrho)=.true.
          lpencil_in(i_cv1)=.true.
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
!
      case default
        call fatal_error('pencil_interdep_eos','ieosvars case not implemented yet')
      endselect
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
!  Calculate EquationOfState pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  06-jan-10/anders: adapted from eos_idealgas
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(:),             intent(IN)   :: lpenc_loc
!
      intent(inout) :: f
      intent(inout) :: p
!
      integer :: i,j
!
      if (leos_isentropic.or.leos_isothermal.or.leos_localisothermal) &
          call fatal_error('pencil_interdep_eos','case not implemented')
!
!  Pencils that are independent of the chosen thermodynamical variables.
!
!
!  fvap (rho_vapour/rho_total). The quantity acc (called q_V in [https://doi.org/10.5194/acp-19-639-2019]),
!  is fvap/(1-fvap).
!  Kishore: Axel, in commit e495afe2fecd170243aa475024da0f3f5b3fd632 ,
!  Kishore: you made a change to treat p%ssat as the vapour fraction.
!  Kishore: However, I think that should be calculated as below.
!  Kishore: Do you agree? I do not see eos_idealgas_vapor being used
!  Kishore: in any samples, and so have taken the liberty of changing it.
!
      if (lpenc_loc(i_fvap)) then
        !p%fvap = p%acc/(1+p%acc)
        !Already calculated in eos_before_boundary
        p%fvap = f(l1:l2,m,n,ifvap)
      endif
      if (lpenc_loc(i_gfvap)) then
        do i=1,3
          p%gfvap(:,i)= p%gacc(:,i)/(1+p%acc)**2
        enddo
      endif
!
!  mumol1
!
      if (lpenc_loc(i_mumol1)) then
        !p%mumol1=(1-p%cc(:,1))*mudry1+p%cc(:,1)*muvap1
        !p%mumol1=(1-p%ssat)*mudry1+p%ssat*muvap1
        !p%mumol1 = (1-p%fvap)*mudry1 + p%fvap*muvap1
        !Already calculated in eos_before_boundary
        p%mumol1 = f(l1:l2,m,n,imumol1)
      endif
! mumol
      if (lpenc_loc(i_mumol)) p%mumol=1/p%mumol1
! glnmumol
      if (lpenc_loc(i_glnmumol)) then
        do i=1,3
          !p%glnmumol(:,i)=-p%mumol*(p%gcc(:,i,1)*(muvap1-mudry1))
          !p%glnmumol(:,i)=-p%mumol*(p%gssat(:,i)*(muvap1-mudry1))
          p%glnmumol(:,i)=-p%mumol*p%gfvap(:,i)*(muvap1-mudry1)
        enddo
      endif
! cp
      if (lpenc_loc(i_cp)) then
        !p%cp=cpdry*mudry*p%mumol1
        !Already calculated in eos_before_boundary
        p%cp = f(l1:l2,m,n,icp)
      endif
! cp1
      if (lpenc_loc(i_cp1)) p%cp1=1/p%cp
! cv
      if (lpenc_loc(i_cv)) p%cv=gamma1*p%cp
! cv1
      if (lpenc_loc(i_cv1)) p%cv1=gamma*p%cp1
! cp1tilde
      if (lpenc_loc(i_cp1tilde)) p%cp1tilde=p%cp1
!
!  Pencils that depend on the chosen thermodynamical variables.
!
      select case (ieosvars)
!
!  Work out thermodynamic quantities for given lnrho or rho and ss.
!
      case (ilnrho_ss,irho_ss)
! ss
        if (lpenc_loc(i_ss)) p%ss=f(l1:l2,m,n,ieosvar2)
! gss
        if (lpenc_loc(i_gss)) call grad(f,iss,p%gss)
! hss
        if (lpenc_loc(i_hss)) call g2ij(f,iss,p%hss)
! del2ss
        if (lpenc_loc(i_del2ss)) call del2(f,iss,p%del2ss)
! del6ss
        if (lpenc_loc(i_del6ss)) call del6(f,iss,p%del6ss)
! lnTT
        if (lpenc_loc(i_lnTT)) p%lnTT=lnTT0+p%cv1*p%ss+gamma_m1*(p%lnrho-lnrho0)
! glnTT
        if (lpenc_loc(i_glnTT)) then
          do i=1,3
            p%glnTT(:,i)=gamma_m1*p%glnrho(:,i)+p%cv1*p%gss(:,i)
          enddo
        endif
! del2lnTT
        if (lpenc_loc(i_del2lnTT)) p%del2lnTT=gamma_m1*p%del2lnrho+p%cv1*p%del2ss
! hlnTT
        if (lpenc_loc(i_hlnTT)) then
          do j=1,3; do i=1,3
            p%hlnTT(:,i,j)=gamma_m1*p%hlnrho(:,i,j)+p%cv1*p%hss(:,i,j)
          enddo; enddo
        endif
!
!  Work out thermodynamic quantities for given lnrho or rho and lnTT.
!
      case (ilnrho_lnTT,irho_lnTT)
! ss
        if (lpenc_loc(i_ss)) call fatal_error('calc_pencils_eos', &
            'ss pencil not implemented')
! gss
        if (lpenc_loc(i_gss)) call fatal_error('calc_pencils_eos', &
            'gss pencil not implemented')
! hss
        if (lpenc_loc(i_hss)) call fatal_error('calc_pencils_eos', &
            'hss pencil not implemented')
! del2ss
        if (lpenc_loc(i_del2ss)) call fatal_error('calc_pencils_eos', &
            'del2ss pencil not implemented')
! del6ss
        if (lpenc_loc(i_del6ss)) call fatal_error('calc_pencils_eos', &
            'del6ss pencil not implemented')
! lnTT
        if (lpenc_loc(i_lnTT)) p%lnTT=f(l1:l2,m,n,ieosvar2)
! glnTT
        if (lpenc_loc(i_glnTT)) call grad(f,ieosvar2,p%glnTT)
! del2lnTT
        if (lpenc_loc(i_del2lnTT)) call del2(f,ieosvar2,p%del2lnTT)
! hlnTT
        if (lpenc_loc(i_hlnTT)) call del6(f,ieosvar2,p%del6lnTT)
!
      case default
        call fatal_error('calc_pencils_eos','case not implemented yet')
      endselect
! TT
      if (lpenc_loc(i_TT)) p%TT=exp(p%lnTT)
! TT1
      if (lpenc_loc(i_TT1)) p%TT1=exp(-p%lnTT)
! pp
      if (lpenc_loc(i_pp)) p%pp=(p%cp-p%cv)*p%rho*p%TT
! ppvap
      !if (lpenc_loc(i_ppvap)) p%ppvap=muvap1*Rgas_unit_sys*p%cc(:,1)*p%rho*p%TT
!       if (lpenc_loc(i_ppvap)) p%ppvap=muvap1*Rgas_unit_sys*p%ssat*p%rho*p%TT
      if (lpenc_loc(i_ppvap)) p%ppvap=p%fvap*mudry1*muvap1*cpdry*p%rho*p%TT
! cs2
      if (lpenc_loc(i_cs2)) p%cs2=p%cp*p%TT*gamma_m1
! csvap2
      if (lpenc_loc(i_csvap2)) p%csvap2=p%cs2*p%mumol*muvap1
! ee
      if (lpenc_loc(i_ee)) p%ee=p%cv*exp(p%lnTT)
! yH
      if (lpenc_loc(i_yH)) p%yH=impossible
!
    endsubroutine calc_pencils_eos_pencpar
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
      real :: EE, TT, yH, rho
!
      call fatal_error('getdensity','not implemented')
!
      call keep_compiler_quiet(EE)
      call keep_compiler_quiet(TT)
      call keep_compiler_quiet(yH)
      call keep_compiler_quiet(rho)
!
    endsubroutine getdensity
!***********************************************************************
    subroutine gettemperature(f,TT_tmp)
!
     real, dimension (mx,my,mz,mfarray), optional :: f
     real, dimension (mx,my,mz), intent(out) :: TT_tmp
!
     call fatal_error('gettemperature','not implemented')
!
     call keep_compiler_quiet(present(f))
     call keep_compiler_quiet(TT_tmp)
!
    endsubroutine gettemperature
!***********************************************************************
    subroutine getpressure(pp_tmp)
!
     real, dimension (mx,my,mz), intent(out) :: pp_tmp
!
     call fatal_error('getpressure','not implemented')
!
     call keep_compiler_quiet(pp_tmp)
!
    endsubroutine getpressure
!***********************************************************************
    subroutine get_gamma_etc(gamma_,cp,cv,f)
!
      real, optional, intent(OUT) :: gamma_, cp,cv
      real, dimension(mfarray), optional, intent(IN) :: f
!
      if (present(gamma_)) gamma_=gamma
!
      if (present(f)) then
        if (present(cp)) cp = f(icp)
        if (present(cv)) cv = f(icp)/gamma
      else
        if (present(cp)) then
          call warning('get_gamma_etc','cp is not constant in eos_idealgas_vapor.'// &
            achar(10)//'The value provided is for one-atomic ideal gas. Use at own risk')
          cp=1.
        endif
        if (present(cv)) then
          call warning('get_gamma_etc','cv is not constant in eos_idealgas_vapor.'// &
            achar(10)//'The value provided is for one-atomic ideal gas. Use at own risk')
          cv=1/gamma
        endif
      endif

    endsubroutine get_gamma_etc
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx) :: cs2, cp1tilde
!
      call fatal_error('pressure_gradient_farray','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(cs2)
      call keep_compiler_quiet(cp1tilde)
!
    endsubroutine pressure_gradient_farray
!***********************************************************************
    subroutine pressure_gradient_point(lnrho,ss,cs2,cp1tilde)
!
      real :: lnrho, ss, cs2, cp1tilde
!
      call fatal_error('pressure_gradient_point','not implemented')
!
      call keep_compiler_quiet(lnrho)
      call keep_compiler_quiet(ss)
      call keep_compiler_quiet(cs2)
      call keep_compiler_quiet(cp1tilde)
!
    endsubroutine pressure_gradient_point
!***********************************************************************
    subroutine temperature_gradient(f,glnrho,gss,glnTT)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx,3) :: glnrho, gss, glnTT
!
      call fatal_error('temperature_gradient','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(glnrho)
      call keep_compiler_quiet(gss)
      call keep_compiler_quiet(glnTT)
!
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_laplacian(f,p)
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call fatal_error('temperature_laplacian','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine temperature_laplacian
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx,3,3) :: hlnrho, hss, hlnTT
!
      call fatal_error('temperature_hessian','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(hlnrho)
      call keep_compiler_quiet(hss)
      call keep_compiler_quiet(hlnTT)
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eosperturb(f,psize,ee,pp,ss)
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: psize
      real, dimension(psize), optional :: ee, pp, ss
!
      call fatal_error('eosperturb','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(psize)
      call keep_compiler_quiet(present(ee))
      call keep_compiler_quiet(present(pp))
      call keep_compiler_quiet(present(ss))
!
    endsubroutine eosperturb
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,yH,lnTT,ee,pp,cs2,kapparho)
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(out), optional :: lnrho
      real, dimension(psize), intent(out), optional :: yH,ee,pp,kapparho
      real, dimension(psize), intent(out), optional :: lnTT
      real, dimension(psize), intent(out), optional :: cs2
!
      real, dimension(psize) :: lnTT_, cs2_
      real, dimension(psize) :: lnrho_,ss_
      real, dimension(psize) :: rho, eth
!     real, dimension(:,:), pointer :: reference_state
!
      select case (ieosvars)
!
!  Log rho and entropy
!
      case (ilnrho_ss,irho_ss)
        call fatal_error('eoscalc_farray','ilnrho_ss not implemented')
!
!  Log rho and Log T
!
      case (ilnrho_lnTT,irho_lnTT)
        select case (psize)
        case (nx)
          lnrho_=f(l1:l2,m,n,ieosvar1)
          lnTT_ =f(l1:l2,m,n,ieosvar2)
        case (mx)
          lnrho_=f(:,m,n,ieosvar1)
          lnTT_ =f(:,m,n,ieosvar2)
        case default
          call fatal_error('eoscalc_farray','no such pencil size')
        end select
!
        if (present(lnrho)) lnrho=lnrho_
        if (present(lnTT)) lnTT=lnTT_
        if (present(ee)) call fatal_error('eoscalc_farray', &
            'ee not implemented for ilnrho_lnTT')
        if (present(pp)) call fatal_error('eoscalc_farray', &
            'pp not implemented for ilnrho_lnTT')
        if (present(cs2)) call fatal_error('eoscalc_farray', &
            'cs2 not implemented for ilnrho_lnTT')
!
      case default
        call fatal_error('eoscalc_farray', &
            'Thermodynamic variable combination  not implemented!')
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
      integer :: ivars
      real :: var1, var2
      real, optional :: lnrho, ss, yH, lnTT, ee, pp, cs2
!
      call fatal_error('eoscalc_point','not implemented')
!
      call keep_compiler_quiet(ivars)
      call keep_compiler_quiet(var1)
      call keep_compiler_quiet(var2)
      call keep_compiler_quiet(present(lnrho))
      call keep_compiler_quiet(present(ss))
      call keep_compiler_quiet(present(yH))
      call keep_compiler_quiet(present(lnTT))
      call keep_compiler_quiet(present(ee))
      call keep_compiler_quiet(present(pp))
      call keep_compiler_quiet(present(cs2))
!
    endsubroutine eoscalc_point
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
      integer :: ivars
      real, dimension(nx) :: var1, var2
      real, dimension(nx), optional :: lnrho, ss, yH, lnTT, ee, pp, cs2
!
      call fatal_error('eoscalc_pencil','not implemented')
!
      call keep_compiler_quiet(ivars)
      call keep_compiler_quiet(var1)
      call keep_compiler_quiet(var2)
      call keep_compiler_quiet(present(lnrho))
      call keep_compiler_quiet(present(ss))
      call keep_compiler_quiet(present(yH))
      call keep_compiler_quiet(present(lnTT))
      call keep_compiler_quiet(present(ee))
      call keep_compiler_quiet(present(pp))
      call keep_compiler_quiet(present(cs2))
!
    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine get_soundspeed(lnTT,cs2)
!
      real :: lnTT, cs2
!
      call fatal_error('get_soundspeed','not implemented')
!
      call keep_compiler_quiet(lnTT)
      call keep_compiler_quiet(cs2)
!
    endsubroutine get_soundspeed
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
      real, intent(in) :: T0
      real, dimension(mx,my,mz), intent(in) :: lnrho
      real, dimension(mx,my,mz), intent(out):: ss
!
      call fatal_error('isothermal_entropy','not implemented')
!
      call keep_compiler_quiet(T0)
      call keep_compiler_quiet(lnrho)
      call keep_compiler_quiet(ss)
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine bc_ss_flux(f,topbot,lone_sided)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
      logical, optional :: lone_sided
!
      call fatal_error('bc_ss_flux','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(lone_sided)
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
!
      call fatal_error('bc_ss_temp_old','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_old
!***********************************************************************
    subroutine bc_ss_temp_x(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
!
      call fatal_error('bc_ss_temp_x','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_y(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
!
      call fatal_error('bc_ss_temp_y','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot,lone_sided)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
      logical, optional :: lone_sided
!
      call fatal_error('bc_ss_temp_z','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(lone_sided)
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_lnrho_temp_z(f,topbot)
!
!  boundary condition for lnrho *and* ss: constant temperature
!
!  27-sep-2002/axel: coded
!  19-aug-2005/tobi: distributed across ionization modules
!  04-dec-2024/Kishore: implemented for eos_idealgas_vapor
!
      use Gravity, only: gravz
      use DensityMethods, only: getlnrho
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: TTval
      integer :: i,il,im
      real, dimension(mx,my) :: lnrho_xy, cp, cv
!
      if (ldebug) print*,'bc_lnrho_temp_z: cs20,cs0=',cs20,cs0
!
!  Constant temperature for entropy, with the temperature corresponding
!  to cs2{top,bot}/((gamma-1)*cpdry). The entropy is set to be antisymmetric
!  about its boundary value, while the density is set assuming
!  hydrostatic equilibrium.
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        do il=1,mx
          do im =1,my
            call get_gamma_etc(cp=cp(il,im), cv=cv(il,im), f=f(il,im,n1,:))
          enddo
        enddo
        TTval = cs2bot/(gamma_m1*cpdry)
        if (ldebug) print*, 'bc_lnrho_temp_z: set z bottom temperature: cs2bot=',cs2bot,"; TTbot=",TTval
        if (cs2bot<=0.) call fatal_error('bc_lnrho_temp_z','cannot have cs2bot<=0')
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        call getlnrho(f(:,:,n1,ilnrho),lnrho_xy)
!
!  This formula works because cp,cv are independent of rho,TT
!
        f(:,:,n1,iss) = cv*(log(TTval)-lnTT0) - (cp-cv)*(lnrho_xy-lnrho0)
        if (lreference_state) call not_implemented('bc_lnrho_temp_z','for lSmag_heat_transport=T')
!
        do i=1,nghost; f(:,:,n1-i,iss) = 2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + (1/cp)*ds/dz = gz/(gamma*(cp-cv)*TTval)
!
        do i=1,nghost
          f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) + (f(:,:,n1+i,iss)-f(:,:,n1-i,iss))/cp - dz2_bound(-i)*gravz/(gamma*(cp-cv)*TTval)
        enddo
!
!  top boundary
!
      case(TOP)
        do il=1,mx
          do im =1,my
            call get_gamma_etc(cp=cp(il,im), cv=cv(il,im), f=f(il,im,n2,:))
          enddo
        enddo
        TTval = cs2top/(gamma_m1*cpdry)
        if (ldebug) print*, 'bc_lnrho_temp_z: set z top temperature: cs2top=',cs2top,"; TTtop=",TTval
        if (cs2top<=0.) call fatal_error('bc_lnrho_temp_z','cannot have cs2top<=0')
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        call getlnrho(f(:,:,n2,ilnrho),lnrho_xy)
!
!  This formula works because cp,cv are independent of rho,TT
!
        f(:,:,n2,iss) = cv*(log(TTval)-lnTT0) - (cp-cv)*(lnrho_xy-lnrho0)
        if (lreference_state) call not_implemented('bc_lnrho_temp_z','for lSmag_heat_transport=T')
!
        do i=1,nghost; f(:,:,n2+i,iss) = 2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + (1/cp)*ds/dz = gz/(gamma*(cp-cv)*TTval)
!
        do i=1,nghost
          f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) + (f(:,:,n2-i,iss)-f(:,:,n2+i,iss))/cp + dz2_bound(i)*gravz/(gamma*(cp-cv)*TTval)
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
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
!
      call fatal_error('bc_lnrho_pressure_z','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_pressure_z
!***********************************************************************
    subroutine bc_ss_temp2_z(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call fatal_error('bc_ss_temp2_z','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
!
      call fatal_error('bc_ss_stemp_x','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_x
!***********************************************************************
    subroutine bc_ss_stemp_y(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
!
      call fatal_error('bc_ss_stemp_y','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
!
      call fatal_error('bc_ss_stemp_z','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
!
      call fatal_error('bc_ss_energy','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_lnrho_hdss_z_liso(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, intent(IN) :: topbot
!
      call fatal_error('bc_lnrho_hdss_z_liso','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hdss_z_liso
!***********************************************************************
    subroutine bc_ss_a2stemp_x(f,topbot)
!
!  11-mar-2012/anders: dummmy
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call fatal_error('bc_ss_a2stemp_x','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_x
!***********************************************************************
    subroutine bc_ss_a2stemp_y(f,topbot)
!
!  11-mar-2012/anders: dummmy
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call fatal_error('bc_ss_a2stemp_y','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_y
!***********************************************************************
    subroutine bc_ss_a2stemp_z(f,topbot)
!
!  11-mar-2012/anders: dummmy
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call fatal_error('bc_ss_a2stemp_z','not implemented')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_z
!***********************************************************************
    subroutine eos_before_boundary(f)
!
!     05-dec-2024/kishore: added
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!     Subroutine get_gamma_etc requires cp to be in the f-array. Usage of
!     get_gamma_etc in boundary condition routines then requires that cp be
!     present in the f-array before the boundary conditions are called.
!     Calculating cp requires fvap and mumol1; we then keep them in the f-array
!     to avoid unnecessary recomputation if these are needed as pencils later on.
!
      f(l1:l2,m1:m2,n1:n2,ifvap) = f(l1:l2,m1:m2,n1:n2,iacc)/(1+f(l1:l2,m1:m2,n1:n2,iacc))
      f(l1:l2,m1:m2,n1:n2,imumol1) = (1-f(l1:l2,m1:m2,n1:n2,ifvap))*mudry1 &
                                     +   f(l1:l2,m1:m2,n1:n2,ifvap)*muvap1
      f(l1:l2,m1:m2,n1:n2,icp) = cpdry*mudry*f(l1:l2,m1:m2,n1:n2,imumol1)
!
    endsubroutine eos_before_boundary
!***********************************************************************
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noeos.f90 for any Eos routines     **
!**  not implemented in this file                                  **
!**                                                                **
    include '../eos_dummies.inc'
!***********************************************************************
endmodule EquationOfState
