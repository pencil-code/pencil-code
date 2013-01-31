! $Id$
!
!  This module takes care of everything related to equation of state.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cv1; cp1; cp1tilde
! PENCILS PROVIDED glnTT(3); TT; TT1; cp; gTT(3)
! PENCILS PROVIDED yH; hss(3,3); hlnTT(3,3); del2ss; del6ss; del2TT; del2lnTT
! PENCILS PROVIDED glnmumol(3); ppvap; csvap2; rho_anel
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
! integers specifying which independent variables to use in eoscalc
  integer, parameter :: ilnrho_ss=1, ilnrho_ee=2, ilnrho_pp=3, ilnrho_lnTT=4
  integer, parameter :: irho_ss=7, ilnrho_TT=9, irho_TT=10, ipp_ss=11
  integer, parameter :: ipp_cs2=12
!
  real, dimension (mz) :: profz_eos=1.0,dprofz_eos=0.0
  real, dimension (3) :: beta_glnrho_global=0.0, beta_glnrho_scaled=0.0
  real :: cp=impossible, cp1=impossible
  real :: cs0=1.0, rho0=1.0
  real :: cs20=1.0, lnrho0=0.0
  real :: gamma=5.0/3.0, gamma_m1=2.0/3.0, gamma1=3.0/5.0
  real :: cs2top_ini=impossible, dcs2top_ini=impossible
  real :: cs2bot=1.0, cs2top=1.0
  real :: cs2cool=0.0
  real :: mpoly=1.5, mpoly0=1.5, mpoly1=1.5, mpoly2=1.5
  integer :: isothtop=1
  logical :: lcalc_cp=.false.
  character (len=labellen) :: ieos_profile='nothing'
  real, dimension(nchemspec,18) :: species_constants
  real, dimension(nchemspec,7)  :: tran_data
  real, dimension(nchemspec)    :: Lewis_coef, Lewis_coef1
!
  contains
!***********************************************************************
    subroutine register_eos()
!
!  14-jun-03/axel: adapted from register_eos
!
      iyH=0
      ilnTT=0
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
!  Dummy.
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos()
!
!  Dummy.
!
    endsubroutine initialize_eos
!***********************************************************************
    subroutine select_eos_variable(variable,findex)
!
!   02-apr-06/tony: dummy
!
      character (len=*), intent(in) :: variable
      integer, intent(in) :: findex
!
      call keep_compiler_quiet(variable)
      call keep_compiler_quiet(findex)
!
    endsubroutine select_eos_variable
!***********************************************************************
    subroutine getmu(f,mu)
!
!  Calculate mean molecular weight.
!
!   12-aug-03/tony: dummy
!
      real, dimension (mx,my,mz,mfarray), optional :: f
      real, intent(out) :: mu
!
      call fatal_error('getmu','SHOULD NOT BE CALLED WITH NOEOS')
!
      mu=0.0
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz), intent(out) :: mu1_full_tmp
!
      call fatal_error('getmu_array','SHOULD NOT BE CALLED WITH NOEOS')
!
      mu1_full_tmp=0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine getmu_array
!***********************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
!  Writes iyH and ilnTT to index.pro file
!
!  02-apr-03/tony: dummy
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
!  02-04-06/tony: dummy
!
    endsubroutine pencil_criteria_eos
!***********************************************************************
    subroutine pencil_interdep_eos(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: dummy
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_eos
!***********************************************************************
    subroutine calc_pencils_eos(f,p)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  02-apr-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
!  Set default values.
!
! SVEN: results should not depend on the values set here!!!
! AXEL: yes, but magnetic has a problematic line
! AXEL: "if (ltemperature) lpenc_requested(i_cv1)=.true."
!
      if (lpencil(i_cv1)) p%cv1=0.0
      if (lpencil(i_cp1)) p%cp1=0.0
      if (lpencil(i_cs2)) p%cs2=cs20
      if (lpencil(i_gTT)) p%gTT=0.0
!      if (lanelastic) then
!        p%lnrho=lnrho0
!        p%rho=exp(p%lnrho)
!      endif
!
      call keep_compiler_quiet(f)
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
    subroutine getdensity(f,ee,TT,yH,rho)
!
      real, intent(in), optional :: ee,TT,yH
      real, dimension (mx,my,mz), intent(inout) :: rho
      real, dimension (mx,my,mz,mfarray), optional :: f
!
      call fatal_error('getdensity','SHOULD NOT BE CALLED WITH NOEOS')
!
      call keep_compiler_quiet(rho)
      call keep_compiler_quiet(present(yH))
      call keep_compiler_quiet(present(ee))
      call keep_compiler_quiet(present(TT))
      call keep_compiler_quiet(present(f))
!
    endsubroutine getdensity
!***********************************************************************
    subroutine gettemperature(f,TT_tmp)
!
     real, dimension (mx,my,mz,mfarray), optional :: f
     real, dimension (mx,my,mz), intent(out) :: TT_tmp
!
     call fatal_error('gettemperature','Should not be called with noeos.')
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
     call fatal_error('getpressure','Should not be called with noeos.')
!
     call keep_compiler_quiet(pp_tmp)
!
    endsubroutine getpressure
!***********************************************************************
    subroutine get_cp1(cp1_)
!
      real, intent(out) :: cp1_
!
      call fatal_error('get_cp1','cp1 is not defined with noeos.f90')
!
      cp1_=0.0
!
    endsubroutine get_cp1
!***********************************************************************
    subroutine get_cv1(cv1_)
!
!  23-dec-10/bing: dummy routine
!
!  return the value of cv1 to outside modules
!
      real, intent(out) :: cv1_
!
      call fatal_error('get_cv1','cv1 is not defined with noeos.f90')
!
      cv1_=0.
!
    endsubroutine get_cv1
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
!   02-apr-04/tony: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: cs2,cp1tilde
!
      call fatal_error('pressure_gradient_farray', &
          'should not be called with noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(cs2,cp1tilde)
!
    endsubroutine pressure_gradient_farray
!***********************************************************************
    subroutine pressure_gradient_point(lnrho,ss,cs2,cp1tilde)
!
!   02-apr-04/tony: dummy
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: cs2,cp1tilde
!
      call fatal_error('pressure_gradient_point', &
          'should not be called with noeos')
!
      call keep_compiler_quiet(lnrho,ss)
      call keep_compiler_quiet(cs2,cp1tilde)
!
    endsubroutine pressure_gradient_point
!***********************************************************************
    subroutine temperature_gradient(f,glnrho,gss,glnTT)
!
!   02-apr-04/tony: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx,3), intent(in) :: glnrho,gss
      real, dimension (nx,3), intent(out) :: glnTT
!
      call fatal_error('temperature_gradient','should not be called with noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(glnrho,glnTT)
      call keep_compiler_quiet(gss)
!
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_laplacian(f,del2lnrho,del2ss,del2lnTT)
!
!   12-dec-05/tony: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(in) :: del2lnrho,del2ss
      real, dimension (nx), intent(out) :: del2lnTT
!
      call fatal_error('temperature_laplacian', &
          'should not be called with noeos')
!
      del2lnTT=0.0
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(del2lnrho)
      call keep_compiler_quiet(del2ss)
!
    endsubroutine temperature_laplacian
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
!   13-may-04/tony: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx,3), intent(in) :: hlnrho,hss
      real, dimension (nx,3) :: hlnTT
!
      call fatal_error('temperature_hessian','now I do not believe you'// &
          ' intended to call this!')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(hlnrho,hss,hlnTT)
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eosperturb(f,psize,ee,pp)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: psize
      real, dimension (psize), intent(in), optional :: ee,pp
!
      call not_implemented('eosperturb')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(psize)
      call keep_compiler_quiet(present(ee))
      call keep_compiler_quiet(present(pp))
!
    endsubroutine eosperturb
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,ss,yH,mu1,lnTT,ee,pp,kapparho)
!
!   02-apr-04/tony: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension (psize), intent(out), optional :: lnrho,ss
      real, dimension (psize), intent(out), optional :: yH,lnTT,mu1
      real, dimension (psize), intent(out), optional :: ee,pp,kapparho
!
      call fatal_error('eoscalc_farray','should not be called with noeos')
!
! To keep compiler quiet set variables with intent(out) to zero
!
      if (present(lnrho)) lnrho = 0.0
      if (present(ss)) ss = 0.0
      if (present(yH)) yH = 0.0
      if (present(lnTT)) lnTT = 0.0
      if (present(mu1)) mu1 = 0.0
      if (present(ee)) ee = 0.0
      if (present(pp)) pp = 0.0
      if (present(kapparho)) kapparho = 0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine eoscalc_farray
!***********************************************************************
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
!   02-apr-04/tony: dummy
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp
!
      call fatal_error('eoscalc_point','should not be called with noeos')
!
      if (present(lnrho)) lnrho=0.0
      if (present(ss)) ss=0.0
      if (present(yH)) yH=0.0
      if (present(lnTT)) lnTT=0.0
      if (present(ee)) ee=0.0
      if (present(pp)) pp=0.0
!
      call keep_compiler_quiet(ivars)
      call keep_compiler_quiet(var1,var2)
!
    endsubroutine eoscalc_point
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
      integer, intent(in) :: ivars
      real, dimension (nx), intent(in) :: var1,var2
      real, dimension (nx), intent(out), optional :: lnrho,ss
      real, dimension (nx), intent(out), optional :: yH,lnTT
      real, dimension (nx), intent(out), optional :: ee,pp
!
      call fatal_error('eoscalc_pencil','should not be called with noeos')
!
      if (present(lnrho)) lnrho=0.0
      if (present(ss)) ss=0.0
      if (present(yH)) yH=0.0
      if (present(lnTT)) lnTT=0.0
      if (present(ee)) ee=0.0
      if (present(pp)) pp=0.0
!
      call keep_compiler_quiet(ivars)
      call keep_compiler_quiet(var1,var2)
!
    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine get_soundspeed(TT,cs2)
!
!  02-apr-04/tony: dummy
!
      real, intent(in)  :: TT
      real, intent(out) :: cs2
!
      call not_implemented('get_soundspeed')
      cs2=0.0
      call keep_compiler_quiet(TT)
!
    endsubroutine get_soundspeed
!***********************************************************************
    subroutine read_eos_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_eos_init_pars
!***********************************************************************
    subroutine write_eos_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_eos_init_pars
!***********************************************************************
    subroutine read_eos_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_eos_run_pars
!***********************************************************************
    subroutine write_eos_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
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
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: T0
      real, dimension (nx) :: lnrho,ss
      real :: ss_offset=0.0
!
!  if T0 is different from unity, we interpret
!  ss_offset = ln(T0)/gamma as an additive offset of ss
!
      if (T0/=1.) ss_offset=alog(T0)/gamma
!
      do n=n1,n2
      do m=m1,m2
        lnrho=f(l1:l2,m,n,ilnrho)
        ss=-gamma_m1*(lnrho-lnrho0)/gamma
          !+ other terms for sound speed not equal to cs_0
        f(l1:l2,m,n,iss)=ss+ss_offset
      enddo
      enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      cs2bot=cs20
      cs2top=cs20
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine isothermal_lnrho_ss(f,T0,rho0)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: T0,rho0
!
      call fatal_error('isothermal_lnrho_ss', &
          'should not be called with noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(T0,rho0)
!
    endsubroutine isothermal_lnrho_ss
!***********************************************************************
    subroutine get_average_pressure(average_density,average_pressure)
!
!  01-dec-2009/piyali+dhruba: dimmy
!
      real, intent(in) :: average_density
      real, intent(out) :: average_pressure
!
      call keep_compiler_quiet(average_density)
      call keep_compiler_quiet(average_pressure)
!
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
      use SharedVariables,only: get_shared_variable
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      logical, pointer :: lmultilayer, lheatc_chiconst
      real, dimension (:,:), allocatable :: tmp_xy,cs2_xy,rho_xy
      integer :: i,ierr,stat
!
      if (ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
!
!  Allocate memory for large arrays.
!
      allocate(tmp_xy(mx,my),stat=stat)
      if (stat>0) call fatal_error('bc_ss_flux', &
          'Could not allocate memory for tmp_xy')
      allocate(cs2_xy(mx,my),stat=stat)
      if (stat>0) call fatal_error('bc_ss_flux', &
          'Could not allocate memory for cs2_xy')
      allocate(rho_xy(mx,my),stat=stat)
      if (stat>0) call fatal_error('bc_ss_flux', &
          'Could not allocate memory for rho_xy')
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
!
!  Get the shared variables
!
      call get_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call fatal_error('bc_ss_flux', &
          'there was a problem when getting hcond0')
      call get_shared_variable('hcond1',hcond1,ierr)
      if (ierr/=0) call fatal_error('bc_ss_flux', &
           'there was a problem when getting hcond1')
      call get_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call fatal_error('bc_ss_flux', &
           'there was a problem when getting Fbot')
      call get_shared_variable('Ftop',Ftop,ierr)
      if (ierr/=0) call fatal_error('bc_ss_flux', &
           'there was a problem when getting Ftop')
      call get_shared_variable('FbotKbot',FbotKbot,ierr)
      if (ierr/=0) call fatal_error('bc_ss_flux', &
           'there was a problem when getting FbotKbot')
      call get_shared_variable('FtopKtop',FtopKtop,ierr)
      if (ierr/=0) call fatal_error('bc_ss_flux', &
           'there was a problem when getting FtopKtop')
      call get_shared_variable('chi',chi,ierr)
      if (ierr/=0) call fatal_error('bc_ss_flux', &
           'there was a problem when getting chi')
      call get_shared_variable('lmultilayer',lmultilayer,ierr)
      if (ierr/=0) call fatal_error('bc_ss_flux' , &
           'there was a problem when getting lmultilayer')
      call get_shared_variable('lheatc_chiconst',lheatc_chiconst,ierr)
      if (ierr/=0) call fatal_error('bc_ss_flux', &
           'there was a problem when getting lheatc_chiconst')
!
      select case (topbot)
!
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
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n1,ilnrho)-lnrho0)+gamma*f(:,:,n1,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
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
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+gamma_m1/gamma* &
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
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n2,ilnrho)-lnrho0)+gamma*f(:,:,n2,iss))
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
          f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+gamma_m1/gamma* &
              (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)-2*i*dz*tmp_xy)
        enddo
      case default
        call fatal_error('bc_ss_flux','invalid argument')
      endselect
!
!  Deallocate arrays.
!
      if (allocated(tmp_xy)) deallocate(tmp_xy)
      if (allocated(cs2_xy)) deallocate(cs2_xy)
      if (allocated(rho_xy)) deallocate(rho_xy)
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_flux_turb(f,topbot)
!
!   4-may-2009/axel: dummy
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb
!***********************************************************************
    subroutine bc_ss_flux_turb_x(f,topbot)
!
!   31-may-2010/pete: dummy
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
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
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      character (len=3), intent(in) :: topbot
!
      real, dimension (:,:), allocatable :: tmp_xy
      integer :: i, stat
!
      if (ldebug) print*,'bc_ss_temp_old: ENTER - cs20,cs0=',cs20,cs0
!
!  Allocate memory for large arrays.
!
      allocate(tmp_xy(mx,my),stat=stat)
      if (stat>0) call fatal_error('bc_ss_temp_old', &
          'Could not allocate memory for tmp_xy')
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
              print*,'bc_ss_temp_old: cannot have cs2bot<=0'
        tmp_xy = (-gamma_m1*(f(:,:,n1,ilnrho)-lnrho0) &
             + alog(cs2bot/cs20)) / gamma
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
                   'bc_ss_temp_old: cannot have cs2top<=0'
  !     if (bcz1(ilnrho) /= 'a2') &
  !          call fatal_error(bc_ss_temp_old','Inconsistent boundary conditions 4.')
        tmp_xy = (-gamma_m1*(f(:,:,n2,ilnrho)-lnrho0) &
                 + alog(cs2top/cs20)) / gamma
        f(:,:,n2,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n2+i,iss) = 2*tmp_xy - f(:,:,n2-i,iss)
        enddo
      case default
        call fatal_error('bc_ss_temp_old','invalid argument')
      endselect
!
!  Deallocate arrays.
!
      if (allocated(tmp_xy)) deallocate(tmp_xy)
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
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(l1,:,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(l1,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
               - gamma_m1/gamma*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                       'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                       'bc_ss_temp_x: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(l2,:,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(l2,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
               - gamma_m1/gamma*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
        enddo
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
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,m1,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,m1,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m1-i,:,iss) = -f(:,m1+i,:,iss) + tmp &
               - gamma_m1/gamma*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                     'bc_ss_temp_y: set y top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                     'bc_ss_temp_y: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,m2,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,m2,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
               - gamma_m1/gamma*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
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
                   'bc_ss_temp_z: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,:,n1,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
               - gamma_m1/gamma*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                     'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp_z: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,:,n2,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
               - gamma_m1/gamma*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
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
        tmp = 2/gamma*log(cs2bot/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n1,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
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
      case ('top')
        if (ldebug) print*, &
                    'bc_lnrho_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*, &
                    'bc_lnrho_temp_z: cannot have cs2top<=0'
        tmp = 2/gamma*log(cs2top/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n2,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
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
        tmp = 1/gamma*alog(cs2bot/cs20)
        do i=0,nghost
          f(:,:,n1-i,iss) = tmp &
               - gamma_m1/gamma*(f(:,:,n1-i,ilnrho)-lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                     'bc_ss_temp2_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp2_z: cannot have cs2top<=0'
        tmp = 1/gamma*alog(cs2top/cs20)
        do i=0,nghost
          f(:,:,n2+i,iss) = tmp &
               - gamma_m1/gamma*(f(:,:,n2+i,ilnrho)-lnrho0)
        enddo
      case default
        call fatal_error('bc_ss_temp2_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_temp3_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  22-jan-2013/axel: coded to impose cs2bot and dcs2bot at bottom
!
      use Gravity, only: gravz
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real :: tmp,dcs2bot
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp3_z: cs20,cs0=',cs20,cs0
!
      call fatal_error('bc_ss_temp3_z','invalid argument')
!
    endsubroutine bc_ss_temp3_z
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
               + gamma_m1/gamma*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
        enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2top<=0'
        do i=1,nghost
          f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
               + gamma_m1/gamma*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
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
               + gamma_m1/gamma*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
        enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2top<=0'
        do i=1,nghost
          f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
               + gamma_m1/gamma*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
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
!  Boundary condition for entropy: symmetric temperature.
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
                  + gamma_m1/gamma*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
          enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0.) print*, &
                 'bc_ss_stemp_z: cannot have cs2top<=0'
         do i=1,nghost
           f(:,:,n2+i,iss) = f(:,:,n2-i,iss) &
                + gamma_m1/gamma*(f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho))
         enddo
      case default
        call fatal_error('bc_ss_stemp_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_a2stemp_x(f,topbot)
!
!  boundary condition for entropy: asymmetric temperature vanishing 2nd deriv
!
!  22-sep-2010/fred: adapted from bc_ss_stemp_z
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_a2stemp_x: cs20,cs0=',cs20,cs0
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
              'bc_ss_a2stemp_x: cannot have cs2bot<=0'
          do i=1,nghost
            f(l1-i,:,:,iss) = min( &
                2*f(l1+1-i,:,:,iss)-f(l1+2-i,:,:,iss)+gamma_m1/gamma* &
                (2*f(l1+1-i,:,:,ilnrho)-f(l1+2-i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)), &
                f(l1+i,:,:,iss)+gamma_m1/gamma* &
                (f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)))
          enddo
!
!  top boundary
!
        case ('top')
          if (cs2top<=0.) print*, &
              'bc_ss_a2stemp_x: cannot have cs2top<=0'
          do i=1,nghost
            f(l2+i,:,:,iss) = min( &
                2*f(l2-1+i,:,:,iss)-f(l2+2-i,:,:,iss)+gamma_m1/gamma* &
                (2*f(l2-1+i,:,:,ilnrho)-f(l2+2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho)), &
                f(l2+i,:,:,iss)+gamma_m1/gamma* &
                (f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho)))
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
!  boundary condition for entropy: asymmetric temperature vanishing 2nd deriv
!
!  22-sep-2010/fred: adapted from bc_ss_stemp_y
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_a2stemp_y: cs20,cs0=',cs20,cs0
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
              'bc_ss_a2stemp_y: cannot have cs2bot<=0'
          do i=1,nghost
            f(:,m1-i,:,iss) = min( &
                2*f(:,m1+1-i,:,iss)-f(:,m1+2-i,:,iss)+gamma_m1/gamma* &
                (2*f(:,m1+1-i,:,ilnrho)-f(:,m1+2-i,:,ilnrho)-f(:,m1-i,:,ilnrho)), &
                f(:,m1-i,:,iss)+gamma_m1/gamma* &
                (f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho)))
          enddo
!
!  top boundary
!
        case ('top')
          if (cs2top<=0.) print*, &
              'bc_ss_a2stemp_y: cannot have cs2top<=0'
          do i=1,nghost
            f(:,m2+i,:,iss) = min( &
                2*f(:,m2-1+i,:,iss)-f(:,m2-2+i,:,iss)+gamma_m1/gamma* &
                (2*f(:,m2-1+i,:,ilnrho)-f(:,m2-2+i,:,ilnrho)-f(:,m2+i,:,ilnrho)), &
                f(:,m2-i,:,iss)+gamma_m1/gamma* &
                (f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho)))
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
!  boundary condition for entropy: asymmetric temperature vanishing 2nd deriv
!
!  22-sep-2010/fred: adapted from bc_ss_stemp_z
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_a2stemp_z: cs20,cs0=',cs20,cs0
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
              'bc_ss_a2stemp_z: cannot have cs2bot<=0'
          do i=1,nghost
            f(:,:,n1-i,iss) = min( &
                2*f(:,:,n1+1-i,iss)-f(:,:,n1+2-i,iss) + gamma_m1/gamma* &
                (2*f(:,:,n1+1-i,ilnrho)-f(:,:,n1+2-i,ilnrho)-f(:,:,n1-i,ilnrho)), &
                f(:,:,n1+i,iss)+gamma_m1/gamma* &
                (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)))
          enddo
!
!  top boundary
!
        case ('top')
          if (cs2top<=0.) print*, &
              'bc_ss_a2stemp_z: cannot have cs2top<=0'
          do i=1,nghost
            f(:,:,n2+i,iss) = min( &
                2*f(:,:,n2-1+i,iss)-f(:,:,n2-2+i,iss) + gamma_m1/gamma* &
                (2*f(:,:,n2-1+i,ilnrho)-f(:,:,n2-2+i,ilnrho)-f(:,:,n2+i,ilnrho)), &
                f(:,:,n2-i,iss)+gamma_m1/gamma* &
                (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)))
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
!  Bottom boundary
!
    case ('bot')
!
!  Set cs2 (temperature) in the ghost points to the value on
!  the boundary
!
      cs2_2d=cs20*exp(gamma_m1*f(:,:,n1,ilnrho)+gamma*f(:,:,n1,iss))
      do i=1,nghost
         f(:,:,n1-i,iss)=1./gamma*(-gamma_m1*f(:,:,n1-i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
!
!  Top boundary
!
    case ('top')
!
!  Set cs2 (temperature) in the ghost points to the value on
!  the boundary
!
      cs2_2d=cs20*exp(gamma_m1*f(:,:,n2,ilnrho)+gamma*f(:,:,n2,iss))
      do i=1,nghost
         f(:,:,n2+i,iss)=1./gamma*(-gamma_m1*f(:,:,n2+i,ilnrho)-log(cs20)&
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
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call fatal_error('bc_stellar_surface','not implemented in eos_idealgas')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: topbot
!
      call fatal_error('bc_lnrho_cfb_r_iso','not implemented in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: topbot
!
      call fatal_error('bc_lnrho_hds_z_iso','not implemented in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hds_z_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=3) :: topbot
!
      call fatal_error('bc_lnrho_hdss_z_iso','not implemented in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
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
