! $Id$
!
!  This module takes care of entropy (initial condition
!  and time advance) for a fluid consisting of gas and perfectly
!  coupled pressureless dust.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .true.
! CPARAM logical, parameter :: ltemperature = .false.
! CPARAM logical, parameter :: lthermal_energy = .false.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cp1tilde; glnTT(3)
! PENCILS PROVIDED TT; TT1; Ma2; ugss; hss(3,3); hlnTT(3,3)
! PENCILS PROVIDED del2ss; del6ss; del2lnTT; cv1
!
!***************************************************************
module Entropy
!
  use Cdata
  use Cparam
  use EquationOfState, only: gamma, gamma_m1, cs20, beta_glnrho_global
  use Interstellar
  use Messages
  use Sub, only: keep_compiler_quiet
  use Viscosity
!
  implicit none
!
  include 'entropy.h'
!
  real :: chi=0.0, chi_t=0.0, chi_shock=0.0, chi_hyper3=0.0
  real :: ss_const=0.0
  real :: T0=1.0
  real :: kx_ss=1.0
  real :: hcond0=impossible
  real :: Kbot=impossible
  integer, parameter :: nheatc_max=4
  logical :: lheatc_Kconst=.false., lheatc_simple=.false.
  logical :: lheatc_chiconst=.false.
  logical :: lheatc_shock=.false., lheatc_hyper3ss=.false.
  logical :: lupw_ss=.false.
  logical, pointer :: lpressuregradient_gas
  logical :: ladvection_entropy=.true.
  character (len=labellen), dimension(ninit) :: initss='nothing'
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  character (len=5) :: iinit_str
!
  namelist /entropy_init_pars/ &
      initss, grads0, ss_const, T0, kx_ss, beta_glnrho_global, &
      ladvection_entropy
!
  namelist /entropy_run_pars/ &
      hcond0, chi_t, chi_shock, chi, iheatcond, Kbot, lupw_ss,chi_hyper3, &
      beta_glnrho_global, ladvection_entropy
!
  integer :: idiag_dtc=0,idiag_ethm=0,idiag_ethdivum=0,idiag_ssm=0
  integer :: idiag_ugradpm=0,idiag_ethtot=0,idiag_dtchi=0,idiag_ssmphi=0
  integer :: idiag_yHm=0,idiag_yHmax=0,idiag_TTm=0,idiag_TTmax=0,idiag_TTmin=0
  integer :: idiag_fconvz=0,idiag_dcoolz=0,idiag_fradz=0,idiag_fturbz=0
  integer :: idiag_ssmz=0,idiag_ssmy=0,idiag_ssmx=0,idiag_TTmz=0
!
  contains
!***********************************************************************
    subroutine register_entropy()
!
!  Initialise variables which should know that we solve an entropy
!  equation: iss, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use FArrayManager
      use SharedVariables, only: get_shared_variables
!
      integer :: ierr
!
      call farray_register_pde('ss',iss)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
      if (lhydro) then
        call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas,ierr)
        if (ierr/=0) call fatal_error('register_entropy','lpressuregradient_gas')
      endif
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  21-jul-2002/wolf: coded
!
      use Gravity, only: gravz,g0
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      integer :: i
      logical :: lnothing
!
!  check any module dependencies.
!
      if (.not. leos) then
        call fatal_error('initialize_entropy','EOS=noeos but entropy requires an EQUATION OF STATE for the fluid')
      endif
      call select_eos_variable('ss',iss)
!
!  For global density gradient beta=H/r*dlnrho/dlnr, calculate actual
!  gradient dlnrho/dr = beta/H.
!
      if (maxval(abs(beta_glnrho_global))/=0.0) then
        beta_glnrho_scaled=beta_glnrho_global*Omega/cs0
        if (lroot) print*, 'initialize_entropy: Global density gradient '// &
            'with beta_glnrho_global=', beta_glnrho_global
      endif
!
!  Turn off pressure gradient term and advection for 0-D runs.
!
      if (nxgrid*nygrid*nzgrid==1) then
        lpressuregradient_gas=.false.
        ladvection_entropy=.false.
        print*, 'initialize_entropy: 0-D run, turned off pressure gradient term'
        print*, 'initialize_entropy: 0-D run, turned off advection of entropy'
      endif
!
!  Initialize heat conduction.
!
      lheatc_Kconst=.false.
      lheatc_simple=.false.
      lheatc_chiconst=.false.
      lheatc_shock=.false.
      lheatc_hyper3ss=.false.
!
      lnothing=.false.
!
!  Select which radiative heating we are using.
!
      if (lroot) print*,'initialize_entropy: nheatc_max,iheatcond=',nheatc_max,iheatcond(1:nheatc_max)
      do i=1,nheatc_max
        select case (iheatcond(i))
        case ('simple')
          lheatc_simple=.true.
          if (lroot) print*, 'heat conduction: simple'
        case ('chi-const')
          lheatc_chiconst=.true.
          if (lroot) print*, 'heat conduction: constant chi'
        case ('shock')
          lheatc_shock=.true.
          if (lroot) print*, 'heat conduction: shock'
        case ('hyper3_ss')
          lheatc_hyper3ss=.true.
          if (lroot) print*, 'heat conduction: hyperdiffusivity of ss'
        case ('nothing')
          if (lroot .and. (.not. lnothing)) print*,'heat conduction: nothing'
        case default
          if (lroot) then
            write(unit=errormsg,fmt=*)  &
                'No such value iheatcond = ', trim(iheatcond(i))
            call fatal_error('initialize_entropy',errormsg)
          endif
        endselect
        lnothing=.true.
      enddo
!
!  A word of warning...
!
      if (lheatc_Kconst .and. hcond0==0.0) then
        call warning('initialize_entropy', 'hcond0 is zero!')
      endif
      if (lheatc_chiconst .and. chi==0.0) then
        call warning('initialize_entropy','chi is zero!')
      endif
      if (all(iheatcond=='nothing') .and. hcond0/=0.0) then
        call warning('initialize_entropy', 'No heat conduction, but hcond0 /= 0')
      endif
      if (lheatc_simple .and. Kbot==0.0) then
        call warning('initialize_entropy','Kbot is zero!')
      endif
      if (lheatc_hyper3ss .and. chi_hyper3==0.0) then
        call warning('initialize_entropy','chi_hyper3 is zero!')
      endif
      if (lheatc_shock .and. chi_shock==0.0) then
        call warning('initialize_entropy','chi_shock is zero!')
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
      endsubroutine initialize_entropy
!***********************************************************************
    subroutine read_entropy_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=entropy_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=entropy_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_entropy_init_pars
!***********************************************************************
    subroutine write_entropy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_init_pars)
!
    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=entropy_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=entropy_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_run_pars)
!
    endsubroutine write_entropy_run_pars
!***********************************************************************
    subroutine init_ss(f)
!
!  Initialise entropy; called from start.f90.
!
      use Sub
      use Gravity
      use General, only: chn
      use Initcond
      use InitialCondition, only: initial_condition_ss
      use EquationOfState,  only: mpoly, isothtop, &
                                rho0, lnrho0, isothermal_entropy, &
                                isothermal_lnrho_ss, eoscalc, ilnrho_pp
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: j
      logical :: lnothing
!
      intent(inout) :: f
!
      lnothing=.true.
      do j=1,ninit
!
        if (initss(j)/='nothing') then
!
          lnothing=.false.
          call chn(j,iinit_str)
!
!  select different initial conditions
!
          select case (initss(j))

            case ('zero', '0'); f(:,:,:,iss) = 0.
            case ('const_ss'); f(:,:,:,iss)=f(:,:,:,iss)+ss_const
            case ('isothermal'); call isothermal_entropy(f,T0)
            case ('isothermal_lnrho_ss')
              print*, 'init_ss: Isothermal density and entropy stratification'
              call isothermal_lnrho_ss(f,T0,rho0)
            case default
!
!  Catch unknown values
!
              write(unit=errormsg,fmt=*) 'No such value for initss(' &
                               //trim(iinit_str)//'): ',trim(initss(j))
              call fatal_error('init_ss',errormsg)
          endselect
!
        endif
!
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_ss(f)
!
    endsubroutine init_ss
!***********************************************************************
    subroutine pencil_criteria_entropy()
!
!  All pencils that the Entropy module depends on are specified here.
!
!  20-11-04/anders: coded
!
      use EquationOfState, only: beta_glnrho_scaled
!
      if (ldt) lpenc_requested(i_cs2)=.true.
      if (lpressuregradient_gas) then
        lpenc_requested(i_cs2)=.true.
        lpenc_requested(i_cp1tilde)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_epsd)=.true.
      endif
      if (ladvection_entropy) lpenc_requested(i_ugss)=.true.
      if (pretend_lnTT) lpenc_requested(i_divu)=.true.
      if (lheatc_simple) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2lnrho)=.true.
        lpenc_requested(i_del2ss)=.true.
      endif
      if (lheatc_Kconst) then
        if (hcond0/=0) then
          lpenc_requested(i_rho1)=.true.
          lpenc_requested(i_glnrho)=.true.
          lpenc_requested(i_gss)=.true.
          lpenc_requested(i_del2lnrho)=.true.
          lpenc_requested(i_del2ss)=.true.
        endif
        if (chi_t/=0) then
          lpenc_requested(i_del2ss)=.true.
        endif
      endif
      if (lheatc_chiconst) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2lnrho)=.true.
        lpenc_requested(i_del2ss)=.true.
      endif
      if (lheatc_shock) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
        lpenc_requested(i_del2ss)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_glnTT)=.true.
      endif
      if (lheatc_hyper3ss) lpenc_requested(i_del6ss)=.true.
      if (lpressuregradient_gas) lpenc_requested(i_cp1tilde)=.true.
!
      if (maxval(abs(beta_glnrho_scaled))/=0.0) lpenc_requested(i_cs2)=.true.
!
      lpenc_diagnos2d(i_ss)=.true.
!
      if (idiag_dtchi/=0) lpenc_diagnos(i_rho1)=.true.
      if (idiag_ethdivum/=0) lpenc_diagnos(i_divu)=.true.
      if (idiag_ssm/=0 .or. idiag_ssmz/=0 .or. idiag_ssmy/=0.or.idiag_ssmx/=0) &
          lpenc_diagnos(i_ss)=.true.
      if (idiag_ethm/=0 .or. idiag_ethtot/=0 .or. idiag_ethdivum/=0) then
          lpenc_diagnos(i_rho)=.true.
          lpenc_diagnos(i_ee)=.true.
      endif
      if (idiag_fconvz/=0 .or. idiag_fturbz/=0 ) then
          lpenc_diagnos(i_rho)=.true.
          lpenc_diagnos(i_TT)=.true.  !(to be replaced by enthalpy)
      endif
      if (idiag_TTm/=0 .or. idiag_TTmz/=0 .or. idiag_TTmax/=0 &
        .or. idiag_TTmin/=0) &
          lpenc_diagnos(i_TT)=.true.
      if (idiag_yHm/=0 .or. idiag_yHmax/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_dtc/=0) lpenc_diagnos(i_cs2)=.true.
!
    endsubroutine pencil_criteria_entropy
!***********************************************************************
    subroutine pencil_interdep_entropy(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
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
      if (lpencil_in(i_ugss)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gss)=.true.
      endif
      if (pretend_lnTT .and. lpencil_in(i_glnTT)) lpencil_in(i_gss)=.true.
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif
      if (lpencil_in(i_hlnTT)) then
        if (pretend_lnTT) then
          lpencil_in(i_hss)=.true.
        else
          lpencil_in(i_hlnrho)=.true.
          lpencil_in(i_hss)=.true.
        endif
      endif
!  The pencils cs2 and cp1tilde come in a bundle, so enough to request one.
      if (lpencil_in(i_cs2) .and. lpencil_in(i_cp1tilde)) &
          lpencil_in(i_cp1tilde)=.false.
!
    endsubroutine pencil_interdep_entropy
!***********************************************************************
    subroutine calc_pencils_entropy(f,p)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use EquationOfState
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: glnrhod, glnrhotot
      real, dimension (nx) :: eps
      integer :: i
!
      intent(in) :: f
      intent(inout) :: p
! ss
      if (lpencil(i_ss)) p%ss=f(l1:l2,m,n,iss)
! gss
      if (lpencil(i_gss)) call grad(f,iss,p%gss)
! pp
      if (lpencil(i_pp)) call eoscalc(f,nx,pp=p%pp)
! ee
      if (lpencil(i_ee)) call eoscalc(f,nx,ee=p%ee)
! lnTT
      if (lpencil(i_lnTT)) call eoscalc(f,nx,lnTT=p%lnTT)
! TT
      if (lpencil(i_TT)) p%TT=exp(p%lnTT)
! TT1
      if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
! cs2 and cp1tilde
      if (lpencil(i_cs2) .or. lpencil(i_cp1tilde)) &
          call pressure_gradient(f,p%cs2,p%cp1tilde)
! Ma2
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
! glnTT
      if (lpencil(i_glnTT)) then
        if (pretend_lnTT) then
           p%glnTT=p%gss
        else
          call temperature_gradient(f,p%glnrho,p%gss,p%glnTT)
        endif
      endif
! ugss
      if (lpencil(i_ugss)) &
          call u_dot_grad(f,iss,p%gss,p%uu,p%ugss,UPWIND=lupw_ss)
!ajwm Should probably combine the following two somehow.
! hss
      if (lpencil(i_hss)) then
        call g2ij(f,iss,p%hss)
      endif
! del2ss
      if (lpencil(i_del2ss)) then
        call del2(f,iss,p%del2ss)
      endif
! del2lnTT
      if (lpencil(i_del2lnTT)) then
          call temperature_laplacian(f,p)
      endif
! del6ss
      if (lpencil(i_del6ss)) then
        call del6(f,iss,p%del6ss)
      endif
! hlnTT
      if (lpencil(i_hlnTT)) then
         if (pretend_lnTT) then
           p%hlnTT=p%hss
         else
           call temperature_hessian(f,p%hlnrho,p%hss,p%hlnTT)
         endif
       endif
! uud (for dust continuity equation in one fluid approximation)
      if (lpencil(i_uud)) p%uud(:,:,1)=p%uu
! divud (for dust continuity equation in one fluid approximation)
      if (lpencil(i_divud)) p%divud(:,1)=p%divu
!
    endsubroutine calc_pencils_entropy
!**********************************************************************
    subroutine dss_dt(f,df,p)
!
!  Calculate right hand side of entropy equation.
!
      use Diagnostics
      use EquationOfState, only: beta_glnrho_global, beta_glnrho_scaled
      use Special, only: special_calc_entropy
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax
      integer :: j,ju
!
      intent(inout)  :: f,p
      intent(out) :: df
!
      Hmax = 0.0
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dss_dt: SOLVE dss_dt'
      if (headtt) call identify_bcs('ss',iss)
!
      if (lhydro) then
!
!  Pressure term in momentum equation (setting lpressuregradient_gas to
!  .false. allows suppressing pressure term for test purposes).
!
        if (lpressuregradient_gas) then
          do j=1,3
            ju=j+iuu-1
            df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) - &
                1/(1+p%epsd(:,1))*p%cs2*(p%glnrho(:,j) + p%cp1tilde*p%gss(:,j))
          enddo
!
!  Add pressure force from global density gradient.
!  WARNING (AJ): This may be implemented inconsistently, since we have
!  here linearised rho and P independently.
!
          if (maxval(abs(beta_glnrho_global))/=0.0) then
            if (headtt) print*, 'dss_dt: adding global pressure gradient force'
              do j=1,3
                df(l1:l2,m,n,(iux-1)+j) = df(l1:l2,m,n,(iux-1)+j) &
                    - 1/(1+p%epsd(:,1))*p%cs2*beta_glnrho_scaled(j)
              enddo
            endif
          endif
        endif
!
!  Advection term.
!
      if (ladvection_entropy) df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%ugss
!
!  Calculate viscous contribution to entropy.
!
      if (lviscosity) call calc_viscous_heat(f,df,p,Hmax)
!
!  Thermal conduction.
!
      if (lheatc_Kconst)   call calc_heatcond(f,df,p)
      if (lheatc_simple)   call calc_heatcond_simple(f,df,p)
      if (lheatc_chiconst) call calc_heatcond_constchi(f,df,p)
      if (lheatc_hyper3ss) call calc_heatcond_hyper3(f,df,p)
!
!  Entry possibility for "personal" entries.
!
      if (lspecial) call special_calc_entropy(f,df,p)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  Calculate entropy related diagnostics.
!
      if (ldiagnos) then
        if (idiag_TTmax/=0) call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        if (idiag_TTm/=0) call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_dtc/=0) &
            call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        if/(idiag_ethm/=0) call sum_mn_name(p%rho*p%ee,idiag_ethm)
        if (idiag_ethtot/=0) call integrate_mn_name(p%rho*p%ee,idiag_ethtot)
        if (idiag_ethdivum/=0) &
            call sum_mn_name(p%rho*p%ee*p%divu,idiag_ethdivum)
        if (idiag_ssm/=0) call sum_mn_name(p%ss,idiag_ssm)
     endif
!
!  1D averages. Happens at every it1d timesteps, NOT at every it1
!  idiag_fradz is done in the calc_headcond routine
!
     if (l1davgfirst) then
        if (idiag_fconvz/=0) call xysum_mn_name_z(p%rho*p%uu(:,3)*p%TT,idiag_fconvz)
        if (idiag_ssmz/=0) call xysum_mn_name_z(p%ss,idiag_ssmz)
        if (idiag_ssmy/=0) call xzsum_mn_name_y(p%ss,idiag_ssmy)
        if (idiag_ssmx/=0) call yzsum_mn_name_x(p%ss,idiag_ssmx)
        if (idiag_TTmz/=0) call xysum_mn_name_z(p%TT,idiag_TTmz)
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_lentropy_pars(f)
!
!  dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lentropy_pars
!***********************************************************************
    subroutine calc_heatcond_constchi(f,df,p)
!
!  Heat conduction for constant value of chi=K/(rho*cp)
!  This routine also adds in turbulent diffusion, if chi_t /= 0.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi*rho*gradT + chi_t*rho*T*grads
!  This routine is currently not correct when ionization is used.
!
!  29-sep-02/axel: adapted from calc_heatcond_simple
!
      use Diagnostics
      use Gravity
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: glnT,glnP
      real, dimension (nx) :: thdiff,g2
!
      intent(in) :: f
      intent(out) :: df
!
!  check that chi is ok
!
      if (headtt) print*,'calc_heatcond_constchi: chi=', chi
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!  The variable g2 is reused to calculate glnP.gss a few lines below.
!
      glnT = gamma*p%gss + gamma_m1*p%glnrho
      glnP = gamma*p%gss + gamma*p%glnrho
      call dot(glnP,glnT,g2)
      thdiff = chi * (gamma*p%del2ss+gamma_m1*p%del2lnrho + g2)
      if (chi_t/=0.) then
        call dot(glnP,p%gss,g2)
        thdiff = thdiff + chi_t*(p%del2ss+g2)
      endif
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_constchi: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+(gamma*chi+chi_t)*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_hyper3(f,df,p)
!
!  Naive hyperdiffusivity of entropy.
!
!  17-jun-05/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: thdiff
!
      intent(in) :: f
      intent(out) :: df
!
!  check that chi_hyper3 is ok
!
      if (headtt) print*, 'calc_heatcond_hyper3: chi_hyper3=', chi_hyper3
!
!  Heat conduction
!
      thdiff = chi_hyper3 * p%del6ss
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_hyper3: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) diffus_chi3=diffus_chi3+chi_hyper3*dxyz_6
!
    endsubroutine calc_heatcond_hyper3
!***********************************************************************
    subroutine calc_heatcond_shock(f,df,p)
!
!  Adds in shock entropy diffusion. There is potential for
!  recycling some quantities from previous calculations.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi_shock*rho*T*grads
!  (in comments we say chi_shock, but in the code this is "chi_shock*shock")
!  This routine should be ok with ionization.
!
!  20-jul-03/axel: adapted from calc_heatcond_constchi
!  19-nov-03/axel: added chi_t also here.
!
      use Diagnostics
      use Gravity
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2,gshockgss
!
      intent(in) :: f,p
      intent(out) :: df
!
!  check that chi is ok
!
      if (headtt) print*,'calc_heatcond_shock: chi_t,chi_shock=',chi_t,chi_shock
!
!  calculate terms for shock diffusion
!  Ds/Dt = ... + chi_shock*[del2ss + (glnchi_shock+glnpp).gss]
!
      call dot(p%gshock,p%gss,gshockgss)
      call dot(p%glnTT+p%glnrho,p%gss,g2)
!
!  shock entropy diffusivity
!  Write: chi_shock = chi_shock0*shock, and gshock=grad(shock), so
!  Ds/Dt = ... + chi_shock0*[shock*(del2ss+glnpp.gss) + gshock.gss]
!
      if (headtt) print*,'calc_heatcond_shock: use shock diffusion'
      thdiff=(chi_shock*p%shock+chi_t)*(p%del2ss+g2)+chi_shock*gshockgss
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_shock: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chi*del2ss
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+(chi_t+chi_shock*p%shock)*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_shock
!***********************************************************************
    subroutine calc_heatcond_simple(f,df,p)
!
!  heat conduction
!
!   8-jul-02/axel: adapted from Wolfgang's more complex version
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: glnT,glnThcond !,glhc
      real, dimension (nx) :: chix
      real, dimension (nx) :: thdiff,g2
      real, dimension (nx) :: hcond
!
      intent(in) :: f,p
      intent(out) :: df
!
!  This particular version assumes a simple polytrope, so mpoly is known
!
      hcond=Kbot
      if (headtt) print*,'calc_heatcond_simple: hcond=', maxval(hcond)
!
!  Heat conduction
!  Note: these routines require revision when ionization turned on
!
      chix = p%rho1*hcond
      glnT = gamma*p%gss + gamma_m1*p%glnrho ! grad ln(T)
      glnThcond = glnT !... + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
      call dot(glnT,glnThcond,g2)
      thdiff = chix * (gamma*p%del2ss+gamma_m1*p%del2lnrho + g2)
!
!  add heat conduction to entropy equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
      if (headtt) print*,'calc_heatcond_simple: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  With heat conduction, the second-order term for entropy is
!  gamma*chix*del2ss
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_simple
!***********************************************************************
    subroutine calc_heatcond(f,df,p)
!
!  heat conduction
!
!  17-sep-01/axel: coded
!  14-jul-05/axel: corrected expression for chi_t diffusion.
!
      use Diagnostics
      use Gravity
      use IO
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: glnT,glnThcond,glhc,glnP,glchit_prof
      real, dimension (nx) :: chix
      real, dimension (nx) :: thdiff,g2
      real, dimension (nx) :: hcond,chit_prof
      real, save :: z_prev=-1.23e20
!
      save :: hcond,glhc
!
      intent(in) :: f,p
      intent(out) :: df
!
!  Heat conduction / entropy diffusion
!
      if (hcond0/=0.0) then
        if (headtt) print*,'calc_heatcond: hcond0=', hcond0
!        if (lgravz) then
!          if (headtt) print*,'calc_heatcond: lgravz=',lgravz
          ! For vertical geometry, we only need to calculate this for each
          ! new value of z -> speedup by about 8% at 32x32x64
!          if (z_mn(1) /= z_prev) then
!            call heatcond(hcond)
!            call gradloghcond(glhc)
!            call chit_profile(chit_prof)
!            call gradlogchit_profile(glchit_prof)
!            z_prev = z_mn(1)
!          endif
!        else
!          call heatcond(hcond)       ! returns hcond=hcond0
!          call gradloghcond(glhc)    ! returns glhc=0
!          call chit_profile(chit_prof)
!          call gradlogchit_profile(glchit_prof)
!        endif
        chix = p%rho1*hcond
        glnT = gamma*p%gss + gamma_m1*p%glnrho             ! grad ln(T)
        glnThcond = glnT + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
        call dot(glnT,glnThcond,g2)
        thdiff = chix * (gamma*p%del2ss+gamma_m1*p%del2lnrho + g2)
      else
        chix   = 0.0
        thdiff = 0.0
        hcond  = 0.0
        glhc   = 0.0
      endif
!
!  write z-profile (for post-processing)
!
      call write_zprof('hcond',hcond)
!
!  Write radiative flux array
!
      if (l1davgfirst) then
        if (idiag_fradz/=0) call xysum_mn_name_z(-hcond*p%TT*glnT(:,3),idiag_fradz)
        if (idiag_fturbz/=0) call xysum_mn_name_z(-chi_t*p%rho*p%TT*p%gss(:,3),idiag_fturbz)
      endif
!
!  "turbulent" entropy diffusion
!  should only be present if g.gradss > 0 (unstable stratification)
!
      if (chi_t/=0.) then
        if (headtt) then
          print*,'calc_headcond: "turbulent" entropy diffusion: chi_t=',chi_t
          if (hcond0 /= 0) then
            call warning('calc_heatcond',"hcond0 and chi_t combined don't seem to make sense")
          endif
        endif
        glnP=gamma*(p%gss+p%glnrho)
        call dot(glnP+glchit_prof,p%gss,g2)
        !thdiff=thdiff+chi_t*(p%del2ss+g2)
        thdiff=thdiff+chi_t*chit_prof*(p%del2ss+g2)
      endif
!
!  check for NaNs initially
!
      if (headt .and. (hcond0 /= 0)) then
        if (notanumber(glhc))      print*,'calc_heatcond: NaNs in glhc'
        if (notanumber(p%rho1))    print*,'calc_heatcond: NaNs in rho1'
        if (notanumber(hcond))     print*,'calc_heatcond: NaNs in hcond'
        if (notanumber(chix))      print*,'calc_heatcond: NaNs in chix'
        if (notanumber(p%del2ss))    print*,'calc_heatcond: NaNs in del2ss'
        if (notanumber(p%del2lnrho)) print*,'calc_heatcond: NaNs in del2lnrho'
        if (notanumber(glhc))      print*,'calc_heatcond: NaNs in glhc'
        if (notanumber(1/hcond))   print*,'calc_heatcond: NaNs in 1/hcond'
        if (notanumber(glnT))      print*,'calc_heatcond: NaNs in glnT'
        if (notanumber(glnThcond)) print*,'calc_heatcond: NaNs in glnThcond'
        if (notanumber(g2))        print*,'calc_heatcond: NaNs in g2'
        if (notanumber(thdiff))    print*,'calc_heatcond: NaNs in thdiff'
!
!  Most of these should trigger the following trap.
!
        if (notanumber(thdiff)) then
          print*, 'calc_heatcond: m,n,y(m),z(n)=',m,n,y(m),z(n)
          call fatal_error('calc_heatcond','NaNs in thdiff')
        endif
      endif

      if (headt .and. lfirst .and. ip<=9) then
        call output_pencil(trim(directory)//'/chi.dat',chix,1)
        call output_pencil(trim(directory)//'/hcond.dat',hcond,1)
        call output_pencil(trim(directory)//'/glhc.dat',glhc,3)
      endif
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + thdiff
!
      if (headtt) print*,'calc_heatcond: added thdiff'
!
!  check maximum diffusion from thermal diffusion
!  NB: With heat conduction, the second-order term for entropy is
!    gamma*chix*del2ss
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+(gamma*chix+chi_t)*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond
!***********************************************************************
    subroutine calc_heatcond_ADI(f)
!
!  Dummy subroutine.
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_heatcond_ADI
!***********************************************************************
    subroutine get_slices_entropy(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_entropy
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics
!
      integer :: iname,inamez,inamey,inamex,irz
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtc=0; idiag_ethm=0; idiag_ethdivum=0; idiag_ssm=0
        idiag_ugradpm=0; idiag_ethtot=0; idiag_dtchi=0; idiag_ssmphi=0
        idiag_yHmax=0; idiag_yHm=0; idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_fconvz=0; idiag_dcoolz=0; idiag_fradz=0; idiag_fturbz=0
        idiag_ssmz=0; idiag_ssmy=0; idiag_ssmx=0; idiag_TTmz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'ethtot',idiag_ethtot)
        call parse_name(iname,cname(iname),cform(iname),'ethdivum',idiag_ethdivum)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',idiag_ugradpm)
        call parse_name(iname,cname(iname),cform(iname),'yHm',idiag_yHm)
        call parse_name(iname,cname(iname),cform(iname),'yHmax',idiag_yHmax)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fturbz',idiag_fturbz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fconvz',idiag_fconvz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'dcoolz',idiag_dcoolz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradz',idiag_fradz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssmz',idiag_ssmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTmz',idiag_TTmz)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'ssmy',idiag_ssmy)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ssmx',idiag_ssmx)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ssmphi',idiag_ssmphi)
      enddo
!
!  Write column where which entropy variable is stored.
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'iss=',iss
        write(3,*) 'iyH=',iyH
        write(3,*) 'ilnTT=',ilnTT
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine fill_farray_pressure(f)
!
!  Fill f array with the pressure, to be able to calculate pressure gradient
!  directly from the pressure.
!
!  18-feb-10/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine fill_farray_pressure
!***********************************************************************
endmodule Entropy
