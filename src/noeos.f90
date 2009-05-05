! $Id$
!
!  This module takes care of everything related to equation of state.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cp1; cp1tilde
! PENCILS PROVIDED glnTT(3); TT; TT1
! PENCILS PROVIDED yH; hss(3,3); hlnTT(3,3); del2ss; del6ss; del2lnTT
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
! integers specifying which independent variables to use in eoscalc
  integer, parameter :: ilnrho_ss=1,ilnrho_ee=2,ilnrho_pp=3,ilnrho_lnTT=4
  integer, parameter :: ilnrho_TT=9
!
  real :: cp=impossible, cp1=impossible
!
  real :: cs0=1., rho0=1.
  real :: cs20=1., lnrho0=0.
  logical :: lcalc_cp=.false.
  real :: gamma=5./3., gamma1=2./3.,gamma11=3./5.
  real :: cs2top_ini=impossible, dcs2top_ini=impossible
  real :: cs2bot=1., cs2top=1.
  real :: cs2cool=0.,ptlaw=impossible
  real :: mpoly=1.5, mpoly0=1.5, mpoly1=1.5, mpoly2=1.5
  integer :: isothtop=1
  real, dimension (3) :: beta_glnrho_global=0.0,beta_glnrho_scaled=0.0
!
  contains
!***********************************************************************
    subroutine register_eos()
!
!  14-jun-03/axel: adapted from register_eos
!
      use Sub
!
      logical, save :: first=.true.
!
      leos=.false.
!
      iyH = 0
      ilnTT = 0
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
!  dummy
!
      gamma1=gamma-1.
      gamma11=1./gamma
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos()
!
!  dummy
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
!
      if (NO_WARN) print*,variable,findex
!
    endsubroutine select_eos_variable
!*******************************************************************
    subroutine getmu(mu)
!
!  Calculate average particle mass in the gas relative to
!
!   12-aug-03/tony: implemented dummy
!
      real, intent(out) :: mu

      mu=0.
    endsubroutine getmu
!*******************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
!  Writes iyH and ilnTT to index.pro file
!
!  02-apr-03/tony: implemented dummy
!
!
      logical :: lreset
      logical, optional :: lwrite
!
      if (NO_WARN) print*,lreset,present(lwrite)  !(keep compiler quiet)
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
!  DUMMY ROUTINE
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
      if (NO_WARN) print*,lpencil_in
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
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      if (lpencil(i_cp1)) p%cp1=0.
!
      if (NO_WARN) print*,f,p
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

      call fatal_error('getdensity','SHOULD NOT BE CALLED WITH NOEOS')
      rho = 0.
      if (NO_WARN) print*,yH,EE,TT  !(keep compiler quiet)

    end subroutine getdensity
!***********************************************************************
    subroutine get_cp1(cp1_)
!
!  return the value of cp1 to outside modules
!
      real, intent(out) :: cp1_
!
      cp1_=impossible
!
    end subroutine get_cp1
!***********************************************************************
    subroutine get_ptlaw(ptlaw_)
!
!  04-jul-07/wlad: return the value of ptlaw to outside modules
!                  ptlaw is temperature gradient in accretion disks
!
      real, intent(out) :: ptlaw_
!
      ptlaw_=impossible
!
    end subroutine get_ptlaw
!***********************************************************************
    subroutine isothermal_density_ion(pot,tmp)
!
      real, dimension (nx), intent(in) :: pot
      real, dimension (nx), intent(out) :: tmp
!
      tmp=0.
      call fatal_error('isothermal_density_ion','SHOULD NOT BE CALLED WITH NOEOS')
      if (NO_WARN) print*,pot  !(keep compiler quiet)
!
    end subroutine isothermal_density_ion
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   02-apr-04/tony: implemented dummy
!
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: cs2,cp1tilde
!
!
      cs2=impossible
      cp1tilde=impossible
      call fatal_error('pressure_gradient_farray','SHOULD NOT BE CALLED WITH NOEOS')
      if (NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine pressure_gradient_farray
!***********************************************************************
    subroutine pressure_gradient_point(lnrho,ss,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   02-apr-04/tony: implemented dummy
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: cs2,cp1tilde
!
      cs2=impossible
      cp1tilde=impossible
      call fatal_error('pressure_gradient_point','SHOULD NOT BE CALLED WITH NOEOS')
      if (NO_WARN) print*,lnrho,ss  !(keep compiler quiet)
!
    endsubroutine pressure_gradient_point
!***********************************************************************
    subroutine temperature_gradient(f,glnrho,gss,glnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   02-apr-04/tony: implemented dummy
!
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(in) :: glnrho,gss
      real, dimension(nx,3), intent(out) :: glnTT
!
!
      glnTT=0.
      call fatal_error('temperature_gradient','SHOULD NOT BE CALLED WITH NOEOS')
      if (NO_WARN) print*,f,glnrho,gss  !(keep compiler quiet)
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_laplacian(f,del2lnrho,del2ss,del2lnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   12-dec-05/tony: adapted from subroutine temperature_gradient
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(in) :: del2lnrho,del2ss
      real, dimension(nx), intent(out) :: del2lnTT
!
      call fatal_error('temperature_laplacian','SHOULD NOT BE CALLED WITH NOEOS')
!
      del2lnTT=0.
      if (NO_WARN) print*,f,del2lnrho,del2ss !(keep compiler quiet)
    endsubroutine temperature_laplacian
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally hlnPP and hlnTT
!   hP/rho=cs2*(hlnrho+cp1tilde*hss)
!
!   13-may-04/tony: adapted from idealgas dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(in) :: hlnrho,hss
      real, dimension(nx,3) :: hlnTT
!
      call fatal_error('temperature_hessian','now I do not believe you'// &
                           ' intended to call this!')
!
      if (NO_WARN) print*,f,hlnrho,hss !(keep compiler quiet)
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eosperturb(f,psize,ee,pp)

      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(in), optional :: ee,pp

      call not_implemented("eosperturb")

    end subroutine eosperturb
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,ss,yH,mu1,lnTT,ee,pp,kapparho)
!
!   Calculate thermodynamical quantities
!
!   02-apr-04/tony: implemented dummy
!
      use Mpicomm, only: stop_it

      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(out), optional :: lnrho,ss
      real, dimension(psize), intent(out), optional :: yH,lnTT,mu1
      real, dimension(psize), intent(out), optional :: ee,pp,kapparho
!
      call fatal_error('eoscalc_farray','SHOULD NOT BE CALLED WITH NOEOS')
!
      lnrho=0.
      ss=0.
      yH=0.
      mu1=0.
      lnTT=0.
      ee=0.
      pp=0.
      kapparho=0.
      if (NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine eoscalc_farray
!***********************************************************************
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
!   Calculate thermodynamical quantities
!
!   02-apr-04/tony: implemented dummy
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp
!
!
      call fatal_error('eoscalc_point','SHOULD NOT BE CALLED WITH NOEOS')

      if (present(lnrho)) lnrho=0.
      if (present(ss)) ss=0.
      if (present(yH)) yH=impossible
      if (present(lnTT)) lnTT=0.
      if (present(ee)) ee=0.
      if (present(pp)) pp=0.
      if (NO_WARN) print*,ivars,var1,var2  !(keep compiler quiet)
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
!
      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), intent(out), optional :: lnrho,ss
      real, dimension(nx), intent(out), optional :: yH,lnTT
      real, dimension(nx), intent(out), optional :: ee,pp
!
!
      call fatal_error('eoscalc_pencil','SHOULD NOT BE CALLED WITH NOEOS')

      if (present(lnrho)) lnrho=0.
      if (present(ss)) ss=0.
      if (present(yH)) yH=impossible
      if (present(lnTT)) lnTT=0.
      if (present(ee)) ee=0.
      if (present(pp)) pp=0.
      if (NO_WARN) print*,ivars,var1,var2  !(keep compiler quiet)
!
    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine get_soundspeed(lnTT,cs2)
!
!  Calculate sound speed for given temperature
!
!  02-apr-04/tony: dummy coded
!
      real, intent(in)  :: lnTT
      real, intent(out) :: cs2
!
      call not_implemented('get_soundspeed')
      cs2=0.
      if (NO_WARN) print*,lnTT
!
    end subroutine get_soundspeed
!***********************************************************************
    subroutine read_eos_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_eos_init_pars
!***********************************************************************
    subroutine write_eos_init_pars(unit)
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine write_eos_init_pars
!***********************************************************************
    subroutine read_eos_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_eos_run_pars
!***********************************************************************
    subroutine write_eos_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit !(keep compiler quiet)

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
      real, dimension(nx) :: lnrho,ss
      real :: ss_offset=0.
!
!  if T0 is different from unity, we interpret
!  ss_offset = ln(T0)/gamma as an additive offset of ss
!
      if (T0/=1.) ss_offset=alog(T0)/gamma
!
      do n=n1,n2
      do m=m1,m2
        lnrho=f(l1:l2,m,n,ilnrho)
        ss=-gamma1*(lnrho-lnrho0)/gamma
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
      use Mpicomm,        only:stop_it
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
        cs2_xy=cs20*exp(gamma1*(f(:,:,n1,ilnrho)-lnrho0)+gamma*f(:,:,n1,iss))
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
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+gamma1/gamma* &
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
        cs2_xy=cs20*exp(gamma1*(f(:,:,n2,ilnrho)-lnrho0)+gamma*f(:,:,n2,iss))
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
          f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+gamma1/gamma* &
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
      if (NO_WARN) print*,f,topbot  !(keep compiler quiet)
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
              print*,'bc_ss_temp_old: cannot have cs2bot<=0'
        tmp_xy = (-gamma1*(f(:,:,n1,ilnrho)-lnrho0) &
             + alog(cs2bot/cs20)) / gamma
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
                   'bc_ss_temp_old: cannot have cs2top<=0'
  !     if (bcz1(ilnrho) /= 'a2') &
  !          call fatal_error(bc_ss_temp_old','Inconsistent boundary conditions 4.')
        tmp_xy = (-gamma1*(f(:,:,n2,ilnrho)-lnrho0) &
                 + alog(cs2top/cs20)) / gamma
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
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_x: set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_x: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(l1,:,:,iss) = 0.5*tmp - gamma1/gamma*(f(l1,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
               - gamma1/gamma*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                       'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                       'bc_ss_temp_x: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(l2,:,:,iss) = 0.5*tmp - gamma1/gamma*(f(l2,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
               - gamma1/gamma*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
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
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,m1,:,iss) = 0.5*tmp - gamma1/gamma*(f(:,m1,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m1-i,:,iss) = -f(:,m1+i,:,iss) + tmp &
               - gamma1/gamma*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                     'bc_ss_temp_y: set y top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                     'bc_ss_temp_y: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,m2,:,iss) = 0.5*tmp - gamma1/gamma*(f(:,m2,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
               - gamma1/gamma*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
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
                   'bc_ss_temp_z: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,:,n1,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
               - gamma1/gamma*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                     'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp_z: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,:,n2,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
               - gamma1/gamma*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
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
        tmp = 2/gamma*log(cs2bot/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n1,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
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
        tmp = 2/gamma*log(cs2top/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n2,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
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
        tmp = 1/gamma*alog(cs2bot/cs20)
        do i=0,nghost
          f(:,:,n1-i,iss) = tmp &
               - gamma1/gamma*(f(:,:,n1-i,ilnrho)-lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                     'bc_ss_temp2_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp2_z: cannot have cs2top<=0'
        tmp = 1/gamma*alog(cs2top/cs20)
        do i=0,nghost
          f(:,:,n2+i,iss) = tmp &
               - gamma1/gamma*(f(:,:,n2+i,ilnrho)-lnrho0)
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
               + gamma1/gamma*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2top<=0'
        do i=1,nghost
          f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
               + gamma1/gamma*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
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
               + gamma1/gamma*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2top<=0'
        do i=1,nghost
          f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
               + gamma1/gamma*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
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
                  + gamma1/gamma*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
          enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                 'bc_ss_stemp_z: cannot have cs2top<=0'
         do i=1,nghost
           f(:,:,n2+i,iss) = f(:,:,n2-i,iss) &
                + gamma1/gamma*(f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho))
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
      cs2_2d=cs20*exp(gamma1*f(:,:,n1,ilnrho)+gamma*f(:,:,n1,iss))
      do i=1,nghost
         f(:,:,n1-i,iss)=1./gamma*(-gamma1*f(:,:,n1-i,ilnrho)-log(cs20)&
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
         f(:,:,n2+i,iss)=1./gamma*(-gamma1*f(:,:,n2+i,ilnrho)-log(cs20)&
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
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_stellar_surface: NOT IMPLEMENTED IN EOS_IDEALGAS")
      if (NO_WARN) print*,f,topbot
!
    end subroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot,j)
      use Mpicomm, only: stop_it
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      call stop_it("bc_lnrho_cfb_r_iso: NOT IMPLEMENTED IN NOEOS")
      if (NO_WARN) print*,f,topbot,j
    end subroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso(f,topbot)
      use Mpicomm, only: stop_it
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      call stop_it("bc_lnrho_hds_z_iso: NOT IMPLEMENTED IN NOEOS")
      if (NO_WARN) print*,f,topbot
    end subroutine bc_lnrho_hds_z_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_liso(f,topbot)
      use Mpicomm, only: stop_it
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      call stop_it("bc_lnrho_hds_z_liso: NOT IMPLEMENTED IN NOEOS")
      if (NO_WARN) print*,f,topbot
    end subroutine bc_lnrho_hds_z_liso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)
      use Mpicomm, only: stop_it
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      call stop_it("bc_lnrho_hdss_z_iso: NOT IMPLEMENTED IN NOEOS")
      if (NO_WARN) print*,f,topbot
    end subroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_liso(f,topbot)
      use Mpicomm, only: stop_it
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      call stop_it("bc_lnrho_hdss_z_liso: NOT IMPLEMENTED IN NOEOS")
      if (NO_WARN) print*,f,topbot
    end subroutine bc_lnrho_hdss_z_liso
!***********************************************************************
endmodule EquationOfState
