! $Id$
!
!  This module takes care of everything related to equation of state.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .false., leos_ionization=.false., leos_temperature_ionization=.false.
! CPARAM logical, parameter :: leos_idealgas = .false., leos_chemistry = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cv1; cp1; cp1tilde
! PENCILS PROVIDED glnTT(3); TT; TT1; cp; cv; gTT(3); mu1; glnmu(3); gmu1(3)
! PENCILS PROVIDED yH; hss(3,3); hlnTT(3,3); del2ss; del6ss; del2TT; del2lnTT; del6TT; del6lnTT
! PENCILS PROVIDED glnmumol(3); csvap2; rho_anel
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
  integer :: imass=1
!
  real :: cs0=1.0, rho0=1.0, rho02
  real :: cs20=1.0, lnrho0=0.0 
  real :: gamma=5.0/3.0, gamma_m1, gamma1
  real :: cs2bot=1.0, cs2top=1.0
  real :: Cp_const=impossible
  real :: Pr_number=0.7
  real :: lnTT0=impossible, TT0=impossible
  logical :: lpres_grad=.false.

  contains
!***********************************************************************
    subroutine register_eos
!
!  14-jun-03/axel: adapted from register_eos
!
      use SharedVariables, only: put_shared_variable
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          '$Id$')

      call put_shared_variable('gamma',gamma,caller='register_eos')
      if (.not.ldensity) then
        call put_shared_variable('rho0',rho0)
        call put_shared_variable('lnrho0',lnrho0)
      endif
!
    endsubroutine register_eos
!***********************************************************************
    subroutine units_eos
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos(f)
!
      real, dimension (mx,my,mz,mfarray),intent(INOUT):: f
!
      gamma_m1=gamma-1.; gamma1=1./gamma
      rho02 = rho0**2
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
    subroutine getmu(f,mu_tmp)
!
!  Calculate mean molecular weight.
!
!   12-aug-03/tony: dummy
!
      real, dimension (mx,my,mz,mfarray), optional :: f
      real, intent(out) :: mu_tmp
!
      call not_implemented('getmu','for this EOS')
!
      mu_tmp=0.0
!
      call keep_compiler_quiet(present(f))
!
    endsubroutine getmu
!***********************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
!  02-apr-03/tony: dummy
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
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
    subroutine pencil_criteria_eos
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
!  02-apr-06/tony: dummy
!  09-oct-15/MR: added mask parameter lpenc_loc.
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      intent(in) :: f,lpenc_loc
      intent(inout) :: p
!
!  Set default values.
!
! SVEN: results should not depend on the values set here!!!
! AXEL: yes, but magnetic has a problematic line
! AXEL: "if (ltemperature) lpenc_requested(i_cv1)=.true."
!
      if (lpenc_loc(i_cv)) p%cv=0.0
      if (lpenc_loc(i_cp)) p%cp=0.0
      if (lpenc_loc(i_cv1)) p%cv1=0.0
      if (lpenc_loc(i_cp1)) p%cp1=0.0
      if (lpenc_loc(i_cs2)) p%cs2=cs20
      if (lpenc_loc(i_gTT)) p%gTT=0.0
      if (lpenc_loc(i_mu1)) p%mu1=0.0
      if (lpenc_loc(i_glnmu)) p%glnmu=1.
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_eos_pencpar
!***********************************************************************
    subroutine init_eos(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_eos
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
      call not_implemented('getdensity','for this EOS')
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
     call not_implemented('gettemperature','for this EOS')
!
     call keep_compiler_quiet(present(f))
     call keep_compiler_quiet(TT_tmp)
!
    endsubroutine gettemperature
!***********************************************************************
    subroutine getpressure(pp_tmp,TT_tmp,rho_tmp,mu1_tmp)
!
     real, dimension (nx), intent(out) :: pp_tmp
     real, dimension (nx), intent(in)  :: TT_tmp,rho_tmp,mu1_tmp
!
     call not_implemented('getpresssure','for this EOS')
!
     call keep_compiler_quiet(pp_tmp)
     call keep_compiler_quiet(TT_tmp)
     call keep_compiler_quiet(rho_tmp)
     call keep_compiler_quiet(mu1_tmp)
!
    endsubroutine getpressure
!***********************************************************************
    subroutine get_gamma_etc(gamma,cp,cv,f)
!
      real, optional, intent(OUT) :: gamma, cp,cv
      real, dimension(mfarray), optional, intent(IN) :: f
!
      call keep_compiler_quiet(f)
!
      if (headt) call warning('get_gamma_etc','gamma, cp, and cv are undefined in noeos.'// &
                              achar(10)//'The values provided are for one-atomic ideal gas. Use at own risk')
      if (present(gamma)) gamma=5./3.
      if (present(cp)) cp=1.
      if (present(cv)) cv=3./5.

    endsubroutine get_gamma_etc
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
!   02-apr-04/tony: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: cs2
      real, dimension (nx), intent(out), optional :: cp1tilde
!
      call not_implemented('presssure_gradient_farray','for this EOS')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(cs2)
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
      call not_implemented('presssure_gradient_point','for this EOS')
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
      call not_implemented('temperature_gradient','for this EOS')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(glnrho,glnTT)
      call keep_compiler_quiet(gss)
!
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_laplacian(f,p)
!
!   12-dec-05/tony: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (pencil_case),                intent(OUT)  :: p
!
      call not_implemented('temperature_laplacian','for this EOS')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
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
      call not_implemented('temperature_hessian','for this EOS')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(hlnrho,hss,hlnTT)
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,yH,mu1,lnTT,ee,pp,cs2,kapparho)
!
!   02-apr-04/tony: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension (psize), intent(out), optional :: lnrho
      real, dimension (psize), intent(out), optional :: yH,lnTT,mu1
      real, dimension (psize), intent(out), optional :: ee,pp,cs2,kapparho
!
      call not_implemented('eoscalc_farray','for this EOS')
!
! To keep compiler quiet set variables with intent(out) to zero
!
      if (present(lnrho)) lnrho = 0.0
      if (present(yH)) yH = 0.0
      if (present(lnTT)) lnTT = 0.0
      if (present(mu1)) mu1 = 0.0
      if (present(ee)) ee = 0.0
      if (present(pp)) pp = 0.0
      if (present(kapparho)) kapparho = 0.0
      if (present(cs2)) cs2 = 0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine eoscalc_farray
!***********************************************************************
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
!   02-apr-04/tony: dummy
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp,cs2
!
      call not_implemented('eoscalc_point','for this EOS')
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
    subroutine eoscalc_point_f(ivars,f,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
!     07-dec-2024/Kishore: dummy
!
      integer, intent(in) :: ivars
      real, dimension(mfarray), intent(in) :: f
      real, optional, intent(out) :: lnrho, ss, yH, lnTT, ee, pp, cs2
!
      call not_implemented('eoscalc_point_f','for this EOS')
!
      if (present(lnrho)) lnrho=0.0
      if (present(ss)) ss=0.0
      if (present(yH)) yH=0.0
      if (present(lnTT)) lnTT=0.0
      if (present(ee)) ee=0.0
      if (present(pp)) pp=0.0
!
      call keep_compiler_quiet(ivars)
      call keep_compiler_quiet(f)
!
    endsubroutine eoscalc_point_f
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
      integer, intent(in) :: ivars
      real, dimension (nx), intent(in) :: var1,var2
      real, dimension (nx), intent(out), optional :: lnrho,ss
      real, dimension (nx), intent(out), optional :: yH,lnTT
      real, dimension (nx), intent(out), optional :: ee,pp,cs2
!
      call not_implemented('eoscalc_pencil','for this EOS')
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
    elemental subroutine get_soundspeed(TT,cs2)
!
!  02-apr-04/tony: dummy
!
      real, intent(in)  :: TT
      real, intent(out) :: cs2
!
    endsubroutine get_soundspeed
!***********************************************************************
    subroutine read_eos_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
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
    subroutine read_eos_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
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
    subroutine isothermal_entropy(lnrho,T0,ss)
!
!  Isothermal stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
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
      real, intent(in) :: T0
      real, dimension(mx,my,mz), intent(in) :: lnrho
      real, dimension(mx,my,mz), intent(out):: ss

      real :: ss_offset
!
!  if T0 is different from unity, we interpret
!  ss_offset = ln(T0)/gamma as an additive offset of ss
!
      ss_offset=0.
      if (T0/=1.) ss_offset=alog(T0)/gamma
!
      ss=-gamma_m1*(lnrho-lnrho0)/gamma+ss_offset
          !+ other terms for sound speed not equal to cs_0
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      cs2bot=cs20
      cs2top=cs20
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine isothermal_lnrho_ss(lnrho,T0,rho0,ss)
!
      real, dimension (mx,my,mz), intent(out) :: lnrho,ss
      real, intent(in) :: T0,rho0
!
      call not_implemented('isothermal_lnrho_ss','for this EOS')
!
      call keep_compiler_quiet(lnrho,ss)
      call keep_compiler_quiet(T0,rho0)
!
    endsubroutine isothermal_lnrho_ss
!***********************************************************************
    subroutine get_average_pressure(average_density,average_pressure)
!
!  01-dec-2009/piyali+dhruba: dummy
!
      real, intent(in) :: average_density
      real, intent(out) :: average_pressure
!
      call not_implemented('get_average_pressure','for this EOS')
!
      call keep_compiler_quiet(average_density)
      call keep_compiler_quiet(average_pressure)
!
    endsubroutine get_average_pressure
!***********************************************************************
    subroutine bc_ss_flux(f,topbot,lone_sided)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!  This implementation is intended to work for all supported EOS.
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!  13-mar-2011/pete: c1 condition for z-boundaries with Kramers' opacity
!   4-jun-2015/MR: factor cp added in front of tmp_xy
!  30-sep-2016/MR: changes for use of one-sided BC formulation (chosen by setting new optional switch lone_sided)
!  21-jun-2024/Kishore: account for bounds on Kramers conductivity
!  07-dec-2024/Kishore: copied from eos_idealgas to eos_idealgas_vapor and modified to be EOS-agnostic. Also reduced duplicated code.
!  09-dec-2024/Kishore: moved to noeos.f90 to replace an earlier implementation that looked correct only for eos_idealgas.
!  TODO: after verifying that this works correctly for all EOS, move it to entropy.f90
!  NOTE: non-dummy implementations (excluding this) are in eos_ionization and eos_temperature_ionization
!
      use DensityMethods, only: getdlnrho_z, getderlnrho_z, getrho
      use Deriv, only: bval_from_neumann, set_ghosts_for_onesided_ders
      use General, only: loptest
      use SharedVariables, only: get_shared_variable
!
      real, pointer :: Flux,FbyK,chi
      real, pointer :: hcond0_kramers, nkramers, chimax_kramers, chimin_kramers
      logical, pointer :: lmultilayer, lheatc_chiconst, lheatc_kramers, lheatc_Kprof, lheatc_Kconst
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      logical, optional :: lone_sided
!
      real, dimension (mx,my) :: FbyKT_xy, TT_xy, rho_xy, Krho1kr_xy, cp, cv, pp_xy
      integer :: i, il, im, ivars, n, ig1, ig2, dir
!
      if (ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
!
      call get_shared_variable('lheatc_chiconst',lheatc_chiconst)
      call get_shared_variable('lheatc_kramers',lheatc_kramers)
      call get_shared_variable('lheatc_Kprof',lheatc_Kprof)
      call get_shared_variable('lheatc_Kconst',lheatc_Kconst)
!
!     TODO: is there a neat way to guard against simultaneous use of heat conductivities (which is currently not supported here)?
!
      if (lheatc_kramers) then
        call get_shared_variable('hcond0_kramers',hcond0_kramers)
        call get_shared_variable('nkramers',nkramers)
        call get_shared_variable('chimax_kramers',chimax_kramers)
        call get_shared_variable('chimin_kramers',chimin_kramers)
      endif
!
      if (lheatc_chiconst) then
        call get_shared_variable('chi',chi)
      endif
!
      select case (topbot)
!
      case(BOT)
!
        if (lheatc_chiconst.or.lheatc_kramers) then
          call get_shared_variable('Fbot',Flux)
        else
          call get_shared_variable('FbotKbot',FbyK)
        endif
        n=n1
        ig1=-1
        ig2=-nghost
        dir=-1
!
      case(TOP)
!
        if (lheatc_chiconst.or.lheatc_kramers) then
          call get_shared_variable('Ftop',Flux)
        else
          call get_shared_variable('FtopKtop',FbyK)
        endif
        n=n2
        ig1=1
        ig2=nghost
        dir=1
!
      case default
        call fatal_error('bc_ss_flux','invalid argument')
      endselect
!
!AXEL: not needed for noeos. Can remove if auto-test ok.
!KG: restored since this subroutine is used other EOS modules (even
!KG: eos_idealgas) through the eos_dummies.inc mechanism. Commenting this out
!KG: breaks most of the convection auto-tests (e.g. conv-slab)
      do il=1,mx
        do im =1,my
          call get_gamma_etc(cp=cp(il,im), cv=cv(il,im), f=f(il,im,n,:))
        enddo
      enddo
!
      if (pretend_lnTT) then
!
        if (.not.(lheatc_Kprof.or.lheatc_Kconst)) call not_implemented('bc_ss_flux', &
          'for this heat conductivity when lpretend_lnTT=T')
!
        FbyKT_xy=-FbyK/exp(f(:,:,n,iss))
        do i=ig1,ig2,dir
          f(:,:,n+i,iss) = f(:,:,n-i,iss) - dir*dz2_bound(i)*FbyKT_xy
        enddo
      else
!
        call getrho(f(:,:,n,ilnrho),rho_xy)
!
        if (ldensity_nolog) then
          ivars=irho_ss
        else
          ivars=ilnrho_ss
        endif
        do il=1,mx
          do im=1,my
            call eoscalc(ivars,f(il,im,n,:),lnTT=TT_xy(il,im), pp=pp_xy(il,im))
          enddo
        enddo
        TT_xy = exp(TT_xy)
!
!  Check whether we have chi=constant at the boundary, in which case
!  we have the nonconstant rho_xy*chi in FbyKT_xy.
!  Check also whether Kramers opacity is used, then hcond itself depends
!  on density and temperature.
!
        if (lheatc_chiconst) then
          !FbyKT_xy=Flux/(rho_xy*chi*TT_xy)
          !Kishore: I think the above misses a factor of cp.
          FbyKT_xy=Flux/(rho_xy*chi*cp*TT_xy)
        else if (lheatc_kramers) then
          Krho1kr_xy = hcond0_kramers*rho_xy**(-2*nkramers-1)*TT_xy**(6.5*nkramers)
!
          if (chimin_kramers>0) Krho1kr_xy = max(Krho1kr_xy, chimin_kramers*cp)
          if (chimax_kramers>0) Krho1kr_xy = min(Krho1kr_xy, chimax_kramers*cp)
!
          FbyKT_xy=Flux/(rho_xy*Krho1kr_xy*TT_xy)
        else
          if (headtt.and..not.(lheatc_Kprof.or.lheatc_Kconst)) then
            call warning('bc_ss_flux', 'FtopKtop and FbotKbot may not be correctly set.')
          endif
          FbyKT_xy=FbyK/TT_xy
        endif
!
!  enforce ds/dz + (p/(rho*TT))*dlnrho/dz = - cv*F/(K*TT)
!
        if (loptest(lone_sided)) then
          call not_implemented('bc_ss_flux', 'one-sided BC')
!           call getderlnrho_z(f,n,rho_xy)                           ! rho_xy=d_z ln(rho)
!           call bval_from_neumann(f,topbot,iss,3,rho_xy)
!           call set_ghosts_for_onesided_ders(f,topbot,iss,3,.true.)
        else
          pp_xy = pp_xy/(rho_xy*TT_xy) !pp_xy is now P/(rho*T)
          FbyKT_xy = dir*cv*FbyKT_xy !now (+-)cv*F/(K*T)
          do i=ig1,ig2,dir
            call getdlnrho_z(f(:,:,:,ilnrho),n,i,rho_xy) !rho_xy is now dlnrho
            f(:,:,n+i,iss) =   f(:,:,n-i,iss) - pp_xy*rho_xy &
                             - dz2_bound(i)*FbyKT_xy
          enddo
        endif
      endif
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_flux_turb(f,topbot)
!
!   4-may-2009/axel: dummy
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call fatal_error('bc_ss_flux_turb','not implemented')
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call fatal_error('bc_ss_flux_turb_x','not implemented')
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
      call fatal_error('bc_ss_flux_condturb_x','not implemented')
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
      call fatal_error('bc_ss_flux_condturb_mean_x','not implemented')
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
      call fatal_error('bc_ss_flux_condturb_z','not implemented')
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
!  26-aug-2003/tony: distributed across ionization modules
!
      real, dimension (:,:,:,:), intent(inout) :: f
      integer, intent(IN) :: topbot
!
      real, dimension (:,:), allocatable :: tmp_xy
      integer :: i, stat
!
      if (ldebug) print*,'bc_ss_temp_old: ENTER - cs20,cs0=',cs20,cs0
!
!  Allocate memory for large arrays.
!
      allocate(tmp_xy(size(f,1),size(f,2)),stat=stat)
      if (stat>0) call fatal_error('bc_ss_temp_old','could not allocate tmp_xy')
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
      case(BOT)
        if ((bcz12(ilnrho,1) /= 'a2') .and. (bcz12(ilnrho,1) /= 'a3')) &
          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3')
        if (ldebug) print*, 'bc_ss_temp_old: set bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*,'bc_ss_temp_old: cannot have cs2bot<=0'
!
        if (ldensity_nolog) then
          tmp_xy = (-gamma_m1*log(f(:,:,n1,irho)/rho0) + log(cs2bot/cs20))*gamma1
        else
          tmp_xy = (-gamma_m1*(f(:,:,n1,ilnrho)-lnrho0) + log(cs2bot/cs20))*gamma1
        endif
        f(:,:,n1,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n1-i,iss) = 2*tmp_xy - f(:,:,n1+i,iss)
        enddo
!
!  top boundary
!
      case(TOP)
        if ((bcz12(ilnrho,2) /= 'a2') .and. (bcz12(ilnrho,2) /= 'a3')) &
          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3')
        if (ldebug) print*, 'bc_ss_temp_old: set top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, 'bc_ss_temp_old: cannot have cs2top<=0'
  !     if (bcz12(ilnrho,1) /= 'a2') &
  !          call fatal_error(bc_ss_temp_old','Inconsistent boundary conditions 4')
!
        if (ldensity_nolog) then
          tmp_xy = (-gamma_m1*log(f(:,:,n2,irho)/rho0) + log(cs2top/cs20))*gamma1
        else
          tmp_xy = (-gamma_m1*(f(:,:,n2,ilnrho)-lnrho0) + log(cs2top/cs20))*gamma1
        endif
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
      case(BOT)
        if (ldebug) print*, 'bc_ss_temp_x: set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, 'bc_ss_temp_x: cannot have cs2bot<=0'
        tmp = 2*gamma1*alog(cs2bot/cs20)
!
        if (ldensity_nolog) then
          f(l1,:,:,iss) = 0.5*tmp - gamma_m1/gamma*log(f(l1,:,:,irho)/rho0)
        else
          f(l1,:,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(l1,:,:,ilnrho)-lnrho0)
        endif
!
        do i=1,nghost
          if (ldensity_nolog) then
            f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
                 - gamma_m1/gamma*log(f(l1+i,:,:,irho)*f(l1-i,:,:,irho)/rho02)
          else
            f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
                 - gamma_m1/gamma*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
          endif
        enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, 'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, 'bc_ss_temp_x: cannot have cs2top<=0'
        tmp = 2*gamma1*alog(cs2top/cs20)
!
        if (ldensity_nolog) then
          f(l2,:,:,iss) = 0.5*tmp - gamma_m1/gamma*log(f(l2,:,:,irho)/rho0)
        else
          f(l2,:,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(l2,:,:,ilnrho)-lnrho0)
        endif
!
        do i=1,nghost
          if (ldensity_nolog) then
            f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
                 - gamma_m1/gamma*log(f(l2-i,:,:,irho)*f(l2+i,:,:,irho)/rho02)
          else
            f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
                 - gamma_m1/gamma*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
          endif
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
      case(BOT)
        if (ldebug) print*, 'bc_ss_temp_y: set y bottom temperature - cs2bot=',cs2bot
        if (cs2bot<=0.) print*, 'bc_ss_temp_y: cannot have cs2bot<=0'
        tmp = 2*gamma1*alog(cs2bot/cs20)
        if (ldensity_nolog) then
          f(:,m1,:,iss) = 0.5*tmp - gamma_m1/gamma*log(f(:,m1,:,irho)/rho0)
        else
          f(:,m1,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,m1,:,ilnrho)-lnrho0)
        endif
        do i=1,nghost
          if (ldensity_nolog) then
            f(:,m1-i,:,iss) = -f(:,m1+i,:,iss) + tmp &
                 - gamma_m1/gamma*log(f(:,m1+i,:,irho)*f(:,m1-i,:,irho)/rho02)
          else
            f(:,m1-i,:,iss) = -f(:,m1+i,:,iss) + tmp &
                 - gamma_m1/gamma*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
          endif
        enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, 'bc_ss_temp_y: set y top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, 'bc_ss_temp_y: cannot have cs2top<=0'
        tmp = 2*gamma1*alog(cs2top/cs20)
        if (ldensity_nolog) then
          f(:,m2,:,iss) = 0.5*tmp - gamma_m1/gamma*log(f(:,m2,:,irho)/rho0)
        else
          f(:,m2,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,m2,:,ilnrho)-lnrho0)
        endif
        do i=1,nghost
          if (ldensity_nolog) then
            f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
                 - gamma_m1/gamma*log(f(:,m2-i,:,irho)*f(:,m2+i,:,irho)/rho02)
          else
            f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
                 - gamma_m1/gamma*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
          endif
        enddo
!
      case default
        call fatal_error('bc_ss_temp_y','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot,lone_sided)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!   3-oct-16/MR: added new optional switch lone_sided
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      logical, optional :: lone_sided
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
      case(BOT)
        if (ldebug) print*, 'bc_ss_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, 'bc_ss_temp_z: cannot have cs2bot<=0'
        tmp = 2*gamma1*alog(cs2bot/cs20)
        if (ldensity_nolog) then
          f(:,:,n1,iss) = 0.5*tmp - gamma_m1/gamma*log(f(:,:,n1,irho)/rho0)
        else
          f(:,:,n1,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
        endif
        do i=1,nghost
          if (ldensity_nolog) then
            f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
                 - gamma_m1/gamma*log(f(:,:,n1+i,irho)*f(:,:,n1-i,irho)/rho02)
          else
            f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
                 - gamma_m1/gamma*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
          endif
        enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, 'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp_z: cannot have cs2top<=0'
        tmp = 2*gamma1*alog(cs2top/cs20)
        if (ldensity_nolog) then
          f(:,:,n2,iss) = 0.5*tmp - gamma_m1/gamma*log(f(:,:,n2,irho)/rho0)
        else
          f(:,:,n2,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        endif
        do i=1,nghost
          if (ldensity_nolog) then
            f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
                 - gamma_m1/gamma*log(f(:,:,n2-i,irho)*f(:,:,n2+i,irho)/rho02)
          else
            f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
                 - gamma_m1/gamma*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
          endif
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
      case(BOT)
        if (ldebug) print*, 'bc_lnrho_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*, 'bc_lnrho_temp_z: cannot have cs2bot<=0'
        tmp = 2*gamma1*log(cs2bot/cs20)
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
                                                  -f(:,:,n1-i,iss)+dz2_bound(-i)*tmp
        enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*,'bc_lnrho_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*,'bc_lnrho_temp_z: cannot have cs2top<=0'
        tmp = 2*gamma1*log(cs2top/cs20)
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
                                                  -f(:,:,n2+i,iss)+dz2_bound(i)*tmp
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
      case(TOP)
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
      case(BOT)
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
      case(BOT)
        if (ldebug) print*, 'bc_ss_temp2_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, 'bc_ss_temp2_z: cannot have cs2bot<=0'
        tmp = gamma1*alog(cs2bot/cs20)
        do i=0,nghost
          if (ldensity_nolog) then
            f(:,:,n1-i,iss) = tmp &
                 - gamma_m1/gamma*log(f(:,:,n1-i,irho)/rho0)
          else
            f(:,:,n1-i,iss) = tmp &
                 - gamma_m1/gamma*(f(:,:,n1-i,ilnrho)-lnrho0)
          endif
        enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, 'bc_ss_temp2_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp2_z: cannot have cs2top<=0'
        tmp = gamma1*alog(cs2top/cs20)
        do i=0,nghost
          if (ldensity_nolog) then
            f(:,:,n2+i,iss) = tmp - gamma_m1/gamma*log(f(:,:,n2+i,irho)/rho0)
          else
            f(:,:,n2+i,iss) = tmp - gamma_m1/gamma*(f(:,:,n2+i,ilnrho)-lnrho0)
          endif
        enddo
      case default
        call fatal_error('bc_ss_temp2_z','invalid argument')
      endselect
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
      call not_implemented('bc_ss_temp3_z','in noeos')
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
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
      case(BOT)
        if (cs2bot<=0.) print*,'bc_ss_stemp_x: cannot have cs2bot<=0'
        do i=1,nghost
          if (ldensity_nolog) then
            f(l1-i,:,:,iss) = f(l1+i,:,:,iss) &
               + gamma_m1/gamma*log(f(l1+i,:,:,irho)/f(l1-i,:,:,irho))
          else
            f(l1-i,:,:,iss) = f(l1+i,:,:,iss) &
               + gamma_m1/gamma*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
          endif
        enddo
!
!  top boundary
!
      case(TOP)
        if (cs2top<=0.) print*,'bc_ss_stemp_x: cannot have cs2top<=0'
        do i=1,nghost
          if (ldensity_nolog) then
            f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
                 + gamma_m1/gamma*log(f(l2-i,:,:,irho)/f(l2+i,:,:,irho))
          else
            f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
                 + gamma_m1/gamma*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
          endif
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
      case(BOT)
        if (cs2bot<=0.) print*,'bc_ss_stemp_y: cannot have cs2bot<=0'
        do i=1,nghost
          if (ldensity_nolog) then
            f(:,m1-i,:,iss) = f(:,m1+i,:,iss) &
                 + gamma_m1/gamma*log(f(:,m1+i,:,irho)/f(:,m1-i,:,irho))
          else
            f(:,m1-i,:,iss) = f(:,m1+i,:,iss) &
                 + gamma_m1/gamma*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
          endif
        enddo
!
!  top boundary
!
      case(TOP)
        if (cs2top<=0.) print*,'bc_ss_stemp_y: cannot have cs2top<=0'
        do i=1,nghost
          if (ldensity_nolog) then
            f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
                 + gamma_m1/gamma*log(f(:,m2-i,:,irho)/f(:,m2+i,:,irho))
          else
            f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
                 + gamma_m1/gamma*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
          endif
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
      case(BOT)
          if (cs2bot<=0.) print*, &
                                  'bc_ss_stemp_z: cannot have cs2bot<=0'
          do i=1,nghost
            if (ldensity_nolog) then
               f(:,:,n1-i,iss) = f(:,:,n1+i,iss) &
                    + gamma_m1/gamma*log(f(:,:,n1+i,irho)/f(:,:,n1-i,irho))
            else
               f(:,:,n1-i,iss) = f(:,:,n1+i,iss) &
                    + gamma_m1/gamma*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
            endif
          enddo
!
!  top boundary
!
      case(TOP)
        if (cs2top<=0.) print*, &
                 'bc_ss_stemp_z: cannot have cs2top<=0'
          do i=1,nghost
            if (ldensity_nolog) then
              f(:,:,n2+i,iss) = f(:,:,n2-i,iss) &
                   + gamma_m1/gamma*log(f(:,:,n2-i,irho)/f(:,:,n2+i,irho))
            else
              f(:,:,n2+i,iss) = f(:,:,n2-i,iss) &
                   + gamma_m1/gamma*(f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho))
            endif
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
        case(BOT)
          if (cs2bot<=0.) call fatal_error('bc_ss_a2stemp_x','cannot have cs2bot<=0')
          do i=1,nghost
            if (ldensity_nolog) then
              f(l1-i,:,:,iss) = min( &
                  2*f(l1+1-i,:,:,iss)-f(l1+2-i,:,:,iss)+gamma_m1/gamma* &
                  log(f(l1+1-i,:,:,irho)**2/f(l1+2-i,:,:,irho)/f(l1-i,:,:,irho)), &
                  f(l1+i,:,:,iss)+gamma_m1/gamma* &
                  log(f(l1+i,:,:,irho)/f(l1-i,:,:,irho)))
            else
              f(l1-i,:,:,iss) = min( &
                  2*f(l1+1-i,:,:,iss)-f(l1+2-i,:,:,iss)+gamma_m1/gamma* &
                  (2*f(l1+1-i,:,:,ilnrho)-f(l1+2-i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)), &
                  f(l1+i,:,:,iss)+gamma_m1/gamma* &
                  (f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)))
            endif
          enddo
!
!  top boundary
!
        case(TOP)
          if (cs2top<=0.) call fatal_error('bc_ss_a2stemp_x','cannot have cs2top<=0')
          do i=1,nghost
            if (ldensity_nolog) then
              f(l2+i,:,:,iss) = min( &
                  2*f(l2-1+i,:,:,iss)-f(l2+2-i,:,:,iss)+gamma_m1/gamma* &
                  log(f(l2-1+i,:,:,irho)**2/f(l2+2-i,:,:,irho)/f(l2+i,:,:,irho)), &
                  f(l2+i,:,:,iss)+gamma_m1/gamma* &
                  log(f(l2-i,:,:,irho)/f(l2+i,:,:,irho)))
            else
              f(l2+i,:,:,iss) = min( &
                  2*f(l2-1+i,:,:,iss)-f(l2+2-i,:,:,iss)+gamma_m1/gamma* &
                  (2*f(l2-1+i,:,:,ilnrho)-f(l2+2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho)), &
                  f(l2+i,:,:,iss)+gamma_m1/gamma* &
                  (f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho)))
            endif
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
        case(BOT)
          if (cs2bot<=0.) call fatal_error('bc_ss_a2stemp_y','cannot have cs2bot<=0')
          do i=1,nghost
            if (ldensity_nolog) then
              f(:,m1-i,:,iss) = min( &
                  2*f(:,m1+1-i,:,iss)-f(:,m1+2-i,:,iss)+gamma_m1/gamma* &
                  log(f(:,m1+1-i,:,irho)**2/f(:,m1+2-i,:,irho)/f(:,m1-i,:,irho)), &
                  f(:,m1-i,:,iss)+gamma_m1/gamma* &
                  log(f(:,m1+i,:,irho)/f(:,m1-i,:,irho)))
            else
              f(:,m1-i,:,iss) = min( &
                  2*f(:,m1+1-i,:,iss)-f(:,m1+2-i,:,iss)+gamma_m1/gamma* &
                  (2*f(:,m1+1-i,:,ilnrho)-f(:,m1+2-i,:,ilnrho)-f(:,m1-i,:,ilnrho)), &
                  f(:,m1-i,:,iss)+gamma_m1/gamma* &
                  (f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho)))
            endif
          enddo
!
!  top boundary
!
        case(TOP)
          if (cs2top<=0.) call fatal_error('bc_ss_a2stemp_y','cannot have cs2top<=0')
          do i=1,nghost
            if (ldensity_nolog) then
              f(:,m2+i,:,iss) = min( &
                  2*f(:,m2-1+i,:,iss)-f(:,m2-2+i,:,iss)+gamma_m1/gamma* &
                  log(f(:,m2-1+i,:,irho)**2/f(:,m2-2+i,:,irho)/f(:,m2+i,:,irho)), &
                  f(:,m2-i,:,iss)+gamma_m1/gamma* &
                  log(f(:,m2-i,:,irho)/f(:,m2+i,:,irho)))
            else
              f(:,m2+i,:,iss) = min( &
                  2*f(:,m2-1+i,:,iss)-f(:,m2-2+i,:,iss)+gamma_m1/gamma* &
                  (2*f(:,m2-1+i,:,ilnrho)-f(:,m2-2+i,:,ilnrho)-f(:,m2+i,:,ilnrho)), &
                  f(:,m2-i,:,iss)+gamma_m1/gamma* &
                  (f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho)))
            endif
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
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
        case(BOT)
          if (cs2bot<=0.) call fatal_error('bc_ss_a2stemp_z','cannot have cs2bot<=0')
          do i=1,nghost
            if (ldensity_nolog) then
              f(:,:,n1-i,iss) = min( &
                  2*f(:,:,n1+1-i,iss)-f(:,:,n1+2-i,iss) + gamma_m1/gamma* &
                  log(f(:,:,n1+1-i,irho)**2/f(:,:,n1+2-i,irho)/f(:,:,n1-i,irho)), &
                  f(:,:,n1+i,iss)+gamma_m1/gamma* &
                  log(f(:,:,n1+i,irho)/f(:,:,n1-i,irho)))
            else
              f(:,:,n1-i,iss) = min( &
                  2*f(:,:,n1+1-i,iss)-f(:,:,n1+2-i,iss) + gamma_m1/gamma* &
                  (2*f(:,:,n1+1-i,ilnrho)-f(:,:,n1+2-i,ilnrho)-f(:,:,n1-i,ilnrho)), &
                  f(:,:,n1+i,iss)+gamma_m1/gamma* &
                  (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)))
            endif
          enddo
!
!  top boundary
!
        case(TOP)
          if (cs2top<=0.) call fatal_error('bc_ss_a2stemp_z','cannot have cs2top<=0')
          do i=1,nghost
            if (ldensity_nolog) then
              f(:,:,n2+i,iss) = min( &
                  2*f(:,:,n2-1+i,iss)-f(:,:,n2-2+i,iss) + gamma_m1/gamma* &
                  log(f(:,:,n2-1+i,irho)**2/f(:,:,n2-2+i,irho)/f(:,:,n2+i,irho)), &
                  f(:,:,n2-i,iss)+gamma_m1/gamma* &
                  log(f(:,:,n2-i,irho)/f(:,:,n2+i,irho)))
            else
              f(:,:,n2+i,iss) = min( &
                  2*f(:,:,n2-1+i,iss)-f(:,:,n2-2+i,iss) + gamma_m1/gamma* &
                  (2*f(:,:,n2-1+i,ilnrho)-f(:,:,n2-2+i,ilnrho)-f(:,:,n2+i,ilnrho)), &
                  f(:,:,n2-i,iss)+gamma_m1/gamma* &
                  (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)))
            endif
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,1),size(f,2)) :: cs2_2d
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
    case(BOT)
!
!  Set cs2 (temperature) in the ghost points to the value on
!  the boundary
!
      if (ldensity_nolog) then
        cs2_2d=cs20*exp(gamma_m1*log(f(:,:,n1,irho))+gamma*f(:,:,n1,iss))
      else
        cs2_2d=cs20*exp(gamma_m1*f(:,:,n1,ilnrho)+gamma*f(:,:,n1,iss))
      endif
      do i=1,nghost
        if (ldensity_nolog) then
          f(:,:,n1-i,iss)= gamma1*(-gamma_m1*log(f(:,:,n1-i,irho))-log(cs20)+log(cs2_2d))
        else
          f(:,:,n1-i,iss)= gamma1*(-gamma_m1*f(:,:,n1-i,ilnrho)-log(cs20)+log(cs2_2d))
        endif
      enddo
!
!  Top boundary
!
    case(TOP)
!
!  Set cs2 (temperature) in the ghost points to the value on
!  the boundary
!
      if (ldensity_nolog) then
        cs2_2d=cs20*exp(gamma_m1*log(f(:,:,n2,irho))+gamma*f(:,:,n2,iss))
      else
        cs2_2d=cs20*exp(gamma_m1*f(:,:,n2,ilnrho)+gamma*f(:,:,n2,iss))
      endif

      do i=1,nghost
        if (ldensity_nolog) then
          f(:,:,n2+i,iss)= gamma1*(-gamma_m1*log(f(:,:,n2+i,irho))-log(cs20)+log(cs2_2d))
        else
          f(:,:,n2+i,iss)= gamma1*(-gamma_m1*f(:,:,n2+i,ilnrho)-log(cs20)+log(cs2_2d))
        endif
      enddo
    case default
      call fatal_error('bc_ss_energy','invalid argument')
    endselect
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_stellar_surface','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call not_implemented('bc_lnrho_cfb_r_iso','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso(f,topbot)
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call not_implemented('bc_lnrho_hds_z_iso','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hds_z_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call not_implemented('bc_lnrho_hdss_z_iso','in noeos')
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
!  unit_length = 1 kpc and scale is 900 pc. To change scale height add to
!  start_pars or run_pars density_scale_factor=... in dimensionless units
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
          else
            call fatal_error('bc_ism','only for irho, ilnrho')
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
          else
            call fatal_error('bc_ism','only for irho, ilnrho')
          endif
        enddo
!
      case default
        print*, "bc_ism: topbot should be BOT or TOP"
!
      endselect
!
    endsubroutine bc_ism
!***********************************************************************
    subroutine write_thermodyn
!
      call not_implemented('write_thermodyn','for this EOS')
!
    endsubroutine write_thermodyn
!***********************************************************************
    subroutine read_thermodyn(input_file)
!
      character (len=*), intent(in) :: input_file
!
      call not_implemented('read_thermodyn','for this EOS')
!
      call keep_compiler_quiet(input_file)
!
    endsubroutine read_thermodyn
!***********************************************************************
    subroutine read_species(input_file)
!
      character (len=*) :: input_file
!
      call not_implemented('read_species','for this EOS')
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
      call not_implemented('find_species_index','for this EOS')
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
      call not_implemented('find_mass','for this EOS')
!
      call keep_compiler_quiet(element_name)
      call keep_compiler_quiet(MolMass)
!
    endsubroutine find_mass
!***********************************************************************
    subroutine get_stratz(z, rho0z, dlnrho0dz, eth0z)
!
!  Get background stratification in z direction.
!
!  13-oct-14/ccyang: dummy
!
      real, dimension(:), intent(in) :: z
      real, dimension(:), intent(out), optional :: rho0z, dlnrho0dz, eth0z
!
      call not_implemented('get_stratz','stratification for this EOS')
!
      call keep_compiler_quiet(z)
      if (present(rho0z)) call keep_compiler_quiet(rho0z)
      if (present(dlnrho0dz)) call keep_compiler_quiet(dlnrho0dz)
      if (present(eth0z)) call keep_compiler_quiet(eth0z)
!
    endsubroutine get_stratz
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=20
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    real, save :: cv, cp, lnTT0

    cv=0.; cp=0.; lnTT0=0.

    call copy_addr(cs20,p_par(1))
    call copy_addr(cv,p_par(2))
    call copy_addr(cp,p_par(3))
    call copy_addr(gamma,p_par(4))
    call copy_addr(lnrho0,p_par(5))
    call copy_addr(lnTT0,p_par(6))
    call copy_addr(rho0,p_par(7))
    call copy_addr(imass,p_par(8)) ! int
    call copy_addr(cp_const,p_par(9))

    endsubroutine pushpars2c
!***********************************************************************
    subroutine eos_before_boundary(f)
!
!     05-dec-2024/Kishore: added dummy routine
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine eos_before_boundary
!***********************************************************************
endmodule EquationOfState
