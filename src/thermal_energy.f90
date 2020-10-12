! $Id$
!
!  This module can replace the energy module by using the thermal energy
!  eth as dependent variable. For a perfect gas we have
!
!    deth/dt + div(u*eth) = -P*div(u)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .false.
! CPARAM logical, parameter :: ltemperature = .false.
! CPARAM logical, parameter :: lthermal_energy = .true.
!
! MVAR CONTRIBUTION 1
!
! PENCILS PROVIDED Ma2; fpres(3); ugeths; transpeth; sglnTT(3)
!
!***************************************************************
module Energy
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'energy.h'
!
  real :: eth_left, eth_right, widtheth, eth_const=1.0
  real :: chi=0.0, chi_shock=0.0, chi_hyper3_mesh=0.
  real :: kappa_rosseland=1.0
  real :: energy_floor = 0.0, temperature_floor = 0.0
  real :: rk_eps = 1.e-3
  real :: nu_z=0., cs_z=0.
  real :: TTref = 0., tau_cool = 1.
  integer :: rk_nmax = 100
  integer :: njeans = 4
  logical :: lsplit_update=.false.
  logical :: lviscosity_heat=.true.
  logical :: lupw_eth=.false.
  logical :: lcheck_negative_energy=.false.
  logical :: lKI02 = .false.
  logical :: lSD93 = .false.
  logical :: lconst_cooling_time = .false.
  logical :: ljeans_floor=.false.
  logical :: lchi_rosseland=.false.
  logical, pointer :: lpressuregradient_gas
  logical :: lenergy_slope_limited=.false.
  character (len=labellen), dimension(ninit) :: initeth='nothing'
  character(len=labellen) :: feedback = 'linear'
!
!  Input parameters.
!
  namelist /entropy_init_pars/ &
      initeth, eth_left, eth_right, widtheth, eth_const
!
!  Run parameters.
!
  namelist /entropy_run_pars/ &
      lviscosity_heat, lupw_eth, &
      chi, chi_shock, chi_hyper3_mesh, &
      lchi_rosseland, kappa_rosseland, &
      energy_floor, temperature_floor, lcheck_negative_energy, &
      rk_eps, rk_nmax, lKI02, lSD93, nu_z, cs_z, &
      lconst_cooling_time, TTref, tau_cool, &
      ljeans_floor, njeans, w_sldchar_ene
!
!  Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_TTmax=0    ! DIAG_DOC: $\max (T)$
  integer :: idiag_TTmin=0    ! DIAG_DOC: $\min (T)$
  integer :: idiag_ppm=0      ! DIAG_DOC: $\left< p \right>$
  integer :: idiag_TTm=0      ! DIAG_DOC: $\left<T\right>$
  integer :: idiag_pdivum=0
  integer :: idiag_ethm=0     ! DIAG_DOC: $\left< e_{\text{th}}\right> =
                              ! DIAG_DOC:  \left< c_v \rho T \right> $
                              ! DIAG_DOC: \quad(mean thermal energy)
  integer :: idiag_ethtot=0   ! DIAG_DOC: $\int_V e_{\text{th}}\,dV$
                              ! DIAG_DOC:   \quad(total thermal energy)
  integer :: idiag_ethmin=0   ! DIAG_DOC: $\mathrm{min} e_\text{th}$
  integer :: idiag_ethmax=0   ! DIAG_DOC: $\mathrm{max} e_\text{th}$
  integer :: idiag_eem=0      ! DIAG_DOC: $\left< e \right> =
                              ! DIAG_DOC:  \left< c_v T \right>$
                              ! DIAG_DOC: \quad(mean internal energy)
  integer :: idiag_etot=0     ! DIAG_DOC: $\langle e_\textrm{th} + \rho u^2 / 2\rangle$
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_ppmz=0     ! XYAVG_DOC: $\left<p\right>_{xy}$
  integer :: idiag_TTmz=0     ! XYAVG_DOC: $\left<T\right>_{xy}$
!
! xz averaged diagnostics given in xzaver.in
!
  integer :: idiag_ppmy=0     ! XZAVG_DOC: $\left<p\right>_{xz}$
  integer :: idiag_TTmy=0     ! XZAVG_DOC: $\left<T\right>_{xz}$
!
! yz averaged diagnostics given in yzaver.in
!
  integer :: idiag_ppmx=0     ! YZAVG_DOC: $\left<p\right>_{yz}$
  integer :: idiag_TTmx=0     ! YZAVG_DOC: $\left<T\right>_{yz}$
!
! y averaged diagnostics given in yaver.in
!
  integer :: idiag_TTmxz=0    ! YAVG_DOC: $\left<T\right>_{y}$
!
! z averaged diagnostics given in zaver.in
!
  integer :: idiag_TTmxy=0    ! ZAVG_DOC: $\left<T\right>_{z}$
!
! variables for slices given in video.in
!
  integer :: ivid_pp=0
  real, dimension(:,:), allocatable :: pp_xz,pp_yz,pp_xy,pp_xy2,pp_xy3,pp_xy4,pp_xz2
!
! General variables for operator split terms.
!
  real :: cv1 = 0.0, cv1_temp = 0.0
!
! Coefficients for the KI02 terms
!
  real, parameter :: KI_heat = 2.0e-26                             ! heating rate [erg/s]
  real, parameter :: KI_v1 = 1.e7, KI_v2 = 1.4e-2                  ! coefficients of the two terms in Lambda / Gamma [cm^3]
  real, parameter :: KI_T1 = 1.184e5, KI_T2 = 1.e3, KI_T3 = 92.    ! coefficients for the temperatures [K]
  real :: KI_a0 = 0., KI_a1 = 0., KI_a2 = 0.
!
! Sutherland & Dopita (1993) cooling function.
!
  real, dimension(:), allocatable :: SD_logTT, SD_logLambda
  integer :: SD_nt = 0
  real :: SD_a0 = 0.
!
! Coefficient for Jeans energy floor
!
  real :: Jeans_c0 = 0.
!
! Background stratification
!
  real, dimension(mz) :: eth0z = 0.0, dlneth0dz = 0.0
  real, dimension(mz) :: rho0z = 0.0
  real, dimension(nx) :: diffus_chi, diffus_chi3
!
  contains
!***********************************************************************
    subroutine register_energy()
!
!  Initialise variables which should know that we solve an energy equation.
!
!  04-nov-10/anders+evghenii: adapted
!
      use FArrayManager, only: farray_register_pde, farray_register_auxiliary
      use SharedVariables, only: get_shared_variable
!
      integer :: ierr
!
      call farray_register_pde('eth',ieth)
!
!  logical variable lpressuregradient_gas shared with hydro modules
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas,ierr)
      if (ierr/=0) call fatal_error('register_energy','lpressuregradient_gas')
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_energy
!***********************************************************************
    subroutine initialize_energy(f)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  04-nov-10/anders+evghenii: adapted
!  03-oct-11/ccyang: add initialization for KI02
!
      !use Slice_methods, only: alloc_slice_buffers
      use EquationOfState, only: select_eos_variable, get_stratz, get_cv1, getmu, gamma, gamma_m1, cs0, cs20
      use SharedVariables
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      real, pointer :: tsg
!
      integer :: istat
      integer :: i, j, k
      real :: mu
      real :: c0, c1
!
      call select_eos_variable('eth',ieth)
!
      call put_shared_variable('lviscosity_heat',lviscosity_heat)
!
!  Decide if operator splitting is required.
!
      lsplit_update = lconst_cooling_time .or. lKI02 .or. lSD93
      if (lsplit_update .and. .not. ldensity) call fatal_error('initialize_energy', 'Density is required for split_update_energy.')
!
!  General variables required by split_update_energy.
!
      ideal_gas: if (lsplit_update) then
        if (.not. leos_idealgas) call fatal_error('initialize_energy', 'currently assumes eos_idealgas')
        call get_cv1(cv1)
        cv1_temp = cv1 * unit_temperature
        if (lstratz) cv1_temp = cv1_temp * cs20 / (gamma * gamma_m1)
      endif ideal_gas
!
!  Initialize the KI02 terms.
!
      KI02: if (lKI02) then
        call getmu(mu_tmp=mu)
        c1 = log10(mu * m_u)
        KI_a0 = log10(KI_heat) + log10(unit_time) - log10(unit_energy)
        c0 = KI_a0 - 3. * log10(unit_length) - 2. * c1
        KI_a0 = KI_a0 - c1
        SI: if (unit_system == 'SI') then
          KI_a0 = KI_a0 - 7.
          c0 = c0 - 13.
        endif SI
        KI_a0 = 10.**KI_a0
        c0 = 10.**c0
        disk: if (nzgrid == 1) then
          c0 = c0 * nu_z / (2. * sqrtpi)
          if (cs_z > 0.) c0 = c0 / cs_z
        endif disk
        KI_a1 = c0 * KI_v1
        KI_a2 = c0 * KI_v2
      endif KI02
!
!  Initialize the SD93 terms.
!
      if (lSD93) call init_cooling_SD93
!
!  Initialize the coefficient of the Jeans energy floor.
!
      Jeans: if (ljeans_floor) then
        if (nzgrid /= 1) then
          call fatal_error('initialize_energy', '3D Jeans floor under construction')
        else
          Jeans_c0 = real(njeans) * G_Newton * dxmax / (gamma * gamma_m1)
        endif
      endif Jeans
!
!  Get background stratification, if any.
!
      if (lstratz) call get_stratz(z, rho0z, dlneth0dz, eth0z)  ! dlnrho0dz = dlneth0dz
!
      if (llocal_iso) &
          call fatal_error('initialize_energy', &
          'llocal_iso switches on the local isothermal approximation. ' // &
          'Use ENERGY=noenergy in src/Makefile.local')
!
      if (ivid_pp/=0) then
        !call alloc_slice_buffers(pp_xy,pp_xz,pp_yz,pp_xy2,pp_xy3,pp_xy4)
        if (lwrite_slice_xy .and..not.allocated(pp_xy) ) allocate(pp_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(pp_xz) ) allocate(pp_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(pp_yz) ) allocate(pp_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(pp_xy2)) allocate(pp_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(pp_xy3)) allocate(pp_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(pp_xy4)) allocate(pp_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(pp_xz2)) allocate(pp_xz2(nx,nz))
      endif

      call keep_compiler_quiet(f)
!
    endsubroutine initialize_energy
!***********************************************************************
    subroutine init_energy(f)
!
!  Initialise thermal energy.
!
!  04-nov-10/anders+evghenii: adapted
!
      use General, only: itoa
      use Initcond, only: jump
      use EquationOfState, only: rho0, cs20, gamma, gamma_m1
      use InitialCondition, only: initial_condition_ss
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      integer :: j
      logical :: lnothing=.true.
      character (len=intlen) :: iinit_str
!
      do j=1,ninit
!
        if (initeth(j)/='nothing') then
!
          lnothing=.false.
!
          iinit_str=itoa(j)
!
!  Select between various initial conditions.
!
          select case (initeth(j))
!
          case ('zero', '0'); f(:,:,:,ieth) = 0.0
!
          case ('xjump'); call jump(f,ieth,eth_left,eth_right,widtheth,'x')
!
          case ('const_eth'); f(:,:,:,ieth)=f(:,:,:,ieth)+eth_const
!
          case ('basic_state')
            if (lstratz) then
              f(:,:,:,ieth) = 0.0
            else
              f(:,:,:,ieth) = rho0 * cs20 / (gamma * gamma_m1)
            endif
!
          case default
!
!  Catch unknown values.
!
            write(unit=errormsg,fmt=*) 'No such value for initeth(' &
                //trim(iinit_str)//'): ',trim(initeth(j))
            call fatal_error('init_energy',errormsg)
!
          endselect
        endif
      enddo
!
!  Call for the InitialCondition facility.
!
      if (linitial_condition) call initial_condition_ss(f)
!
    endsubroutine init_energy
!***********************************************************************
    subroutine pencil_criteria_energy()
!
!  All pencils that the Energy module depends on are specified here.
!
!  04-nov-10/anders+evghenii: adapted
!
      lpenc_requested(i_divu) = .true.
      lpenc_requested(i_pp) = .true.
      stratz: if (lstratz) then
        lpenc_requested(i_eths) = .true.
        lpenc_requested(i_ugeths) = .true.
        lpenc_requested(i_uu) = .true.
      else stratz
        weno: if (lweno_transport) then
          lpenc_requested(i_transpeth) = .true.
        else weno
          lpenc_requested(i_geth) = .true.
          lpenc_requested(i_eth) = .true.
          lpenc_requested(i_uu) = .true.
        endif weno
        lpenc_requested(i_pp) = .true.
      endif stratz
      if (lhydro) lpenc_requested(i_fpres)=.true.
      if (ldt) lpenc_requested(i_cs2)=.true.
      if (lviscosity.and.lviscosity_heat) lpenc_requested(i_visc_heat)=.true.
      if (chi/=0.0) then
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_cp)=.true.
        lpenc_requested(i_del2TT)=.true.
        lpenc_requested(i_grho)=.true.
        lpenc_requested(i_gTT)=.true.
      endif
      if (lchi_rosseland) then
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_del2TT)=.true.
        lpenc_requested(i_grho)=.true.
        lpenc_requested(i_gTT)=.true.
      endif
!
      shock: if (chi_shock /= 0.0) then
        lpenc_requested(i_shock) = .true.
        lpenc_requested(i_gshock) = .true.
        stratz1: if (lstratz) then
          lpenc_requested(i_geths) = .true.
        else stratz1
          lpenc_requested(i_geth) = .true.
          lpenc_requested(i_del2eth) = .true.
        endif stratz1
      endif shock
!
!  Diagnostic pencils.
!
      if (idiag_ethm/=0 .or. idiag_ethmin/=0 .or. idiag_ethmax/=0 .or. idiag_ethtot/=0) lpenc_diagnos(i_eth)=.true.
      if (idiag_eem/=0) lpenc_diagnos(i_ee)=.true.
      etot: if (idiag_etot /= 0) then
        lpenc_diagnos(i_eth) = .true.
        lpenc_diagnos(i_ekin) = .true.
      endif etot
      if (idiag_ppm/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_pdivum/=0) then
        lpenc_diagnos(i_pp)=.true.
        lpenc_diagnos(i_divu)=.true.
      endif
!
      if (idiag_TTm/=0   .or. idiag_TTmin/=0 .or. idiag_TTmax/=0 .or. &
          idiag_TTmxy/=0 .or. idiag_TTmxz/=0 .or. idiag_TTmx/=0  .or. &
          idiag_TTmy/=0  .or. idiag_TTmz/=0 ) lpenc_diagnos(i_TT)=.true.
!
    endsubroutine pencil_criteria_energy
!***********************************************************************
    subroutine pencil_interdep_energy(lpencil_in)
!
!  Interdependency among pencils from the Energy module is specified here.
!
!  04-nov-10/anders+evghenii: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      fpres: if (lpencil_in(i_fpres)) then
        lpencil_in(i_rho1) = .true.
        stratz: if (lstratz) then
          lpencil_in(i_geths) = .true.
        else stratz
          lpencil_in(i_geth) = .true.
        endif stratz
      endif fpres
!
      ugeths: if (lpencil_in(i_ugeths)) then
        lpencil_in(i_geths) = .true.
        lpencil_in(i_uu) = .true.
      endif ugeths
!
    endsubroutine pencil_interdep_energy
!***********************************************************************
    subroutine calc_pencils_energy(f,p)
!
!  Calculate Energy pencils.
!  This routine is called after  calc_pencils_eos
!  Most basic pencils should come first, as others may depend on them.
!
!  04-nov-10/anders+evghenii: adapted
!  14-feb-11/bing: moved eth dependend pecnils to eos_idealgas
!
      use EquationOfState, only: gamma, gamma_m1, cs20
      use Sub, only: u_dot_grad
      use WENO_transport, only: weno_transp
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension(nx) :: penc
      integer :: j
!
! fpres
!
      fpres: if (lpencil(i_fpres)) then
        stratz: if (lstratz) then
          penc = gamma_m1 * eth0z(n) * p%rho1
          p%fpres = -p%geths
          p%fpres(:,3) = p%fpres(:,3) + (f(l1:l2,m,n,irho) - f(l1:l2,m,n,ieth)) * dlneth0dz(n)
          forall(j = 1:3) p%fpres(:,j) = penc * p%fpres(:,j)
        else stratz
          forall(j = 1:3) p%fpres(:,j) = -gamma_m1 * p%rho1 * p%geth(:,j)
        endif stratz
      endif fpres
!
! ugeths
!
      if (lpencil(i_ugeths)) call u_dot_grad(f, ieth, p%geths, p%uu, p%ugeths)
!
! transpeth
!
      if (lpencil(i_transpeth)) &
          call weno_transp(f,m,n,ieth,-1,iux,iuy,iuz,p%transpeth,dx_1,dy_1,dz_1)
! sglnTT 
      if (lpencil(i_sglnTT)) &
        call fatal_error('calc_pencils_energy', &
            'Pencil sglnTT not yet implemented for thermal_energy')
!
    endsubroutine calc_pencils_energy
!***********************************************************************
    subroutine denergy_dt(f,df,p)
!
!  Calculate right hand side of energy equation.
!
!    deth/dt + div(u*eth) = -P*div(u)
!
!  04-nov-10/anders+evghenii: coded
!  02-aug-11/ccyang: add mesh hyper-diffusion
!
      use EquationOfState, only: gamma
      use Special, only: special_calc_energy
      use Sub, only: identify_bcs, u_dot_grad, del2
      use Viscosity, only: calc_viscous_heat
      use Deriv, only: der6
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(nx) :: Hmax=0.0, ugeth, d6eth, penc
      integer :: j
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*, 'denergy_dt: solve deth_dt'
      if (headtt) call identify_bcs('eth',ieth)
!
!  Sound speed squared.
!
      if (headtt) print*, 'denergy_dt: cs20=', p%cs2(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (ldensity.and.lhydro.and.lfirst.and.ldt) then
        advec_cs2=p%cs2*dxyz_2
        if (headtt.or.ldebug) print*,'denergy_dt: max(advec_cs2) =',maxval(advec_cs2)
      endif
!
!  Add pressure gradient term in momentum equation.
!
      if (lhydro) df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%fpres
!
!  Entry possibility for "personal" entries.
!  In that case you'd need to provide your own "special" routine.
!
      if (lspecial) call special_calc_energy(f,df,p)
!
      stratz: if (lstratz) then
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) - p%ugeths - p%eths * (gamma * p%divu + p%uu(:,3) * dlneth0dz(n))
      else stratz
!
!  Add energy transport term.
!
        if (lweno_transport) then
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) - p%transpeth
        else
          call u_dot_grad(f, ieth, p%geth, p%uu, ugeth, UPWIND=lupw_eth)
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) - p%eth*p%divu - ugeth
        endif
!
!  Add P*dV work.
!
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) - p%pp*p%divu
!
      endif stratz
!
!  Calculate viscous contribution to temperature.
!
      if (lviscosity.and.lviscosity_heat) call calc_viscous_heat(df,p,Hmax)
!
!  Thermal energy diffusion.
!
      diffus_chi=0; diffus_chi3=0.
      if (chi/=0.0) then
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + chi * p%cp * (p%rho*p%del2TT + sum(p%grho*p%gTT, 2))
        if (lfirst .and. ldt) diffus_chi = diffus_chi + gamma*chi*dxyz_2
      endif
!
!  Shock diffusion
!
      shock: if (chi_shock /= 0.0) then
        stratz1: if (lstratz) then
          call del2(f, ieth, penc)
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + chi_shock * (p%shock * penc + sum(p%gshock*p%geths, 2))
        else stratz1
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + chi_shock * (p%shock * p%del2eth + sum(p%gshock*p%geth, 2))
        endif stratz1
        if (lfirst .and. ldt) diffus_chi = diffus_chi + chi_shock * p%shock * dxyz_2
      endif shock
!
!  Mesh hyper-diffusion
!
      if (chi_hyper3_mesh/=0.0) then
        do j=1,3
          call der6(f, ieth, d6eth, j, IGNOREDX=.true.)
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + &
              chi_hyper3_mesh * d6eth * dline_1(:,j)
        enddo
        if (lfirst .and. ldt) diffus_chi3 = diffus_chi3 + chi_hyper3_mesh*sum(dline_1,2)                  
      endif
!
!  Radiative diffusion (Rosseland approximation) through thermal energy diffusion.
!
      if (lchi_rosseland) then
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + ( 16.0*sigmaSB/(3.0*kappa_rosseland))*p%TT*p%TT*p%rho1*(  &
                                 3.0*sum(p%gTT*p%gTT, 2) &
                                - p%TT*p%rho1*sum(p%grho*p%gTT, 2) &
                                + p%TT*p%del2TT )
!       This expression has an extra TT/eth = cv1*rho1  to convert to the right units - needed because of the sigmaSB
!       The timestep gets the factor TT/eth = cv1*rho1
        if (lfirst .and. ldt) diffus_chi = diffus_chi &
            +(16.0*sigmaSB/(3.0*kappa_rosseland))*p%TT*p%TT*p%TT/p%rho *cv1*p%rho1 * dxyz_2
     endif
     if (lfirst .and. ldt) then
       maxdiffus=max(maxdiffus,diffus_chi)
       maxdiffus3=max(maxdiffus3,diffus_chi3)
     endif
!
!  Diagnostics.
!
     call calc_diagnostics_energy(f,p)
!
    endsubroutine denergy_dt
!***********************************************************************
    subroutine calc_diagnostics_energy(f,p)

      use Diagnostics
      use Slices_methods, only: store_slices

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      call keep_compiler_quiet(f)

      if (ldiagnos) then
        call sum_mn_name(p%TT,idiag_TTm)
        call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0)  call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        call sum_mn_name(p%eth,idiag_ethm)
        call integrate_mn_name(p%eth,idiag_ethtot)
        if (idiag_ethmin/=0) call max_mn_name(-p%eth,idiag_ethmin,lneg=.true.)
        call max_mn_name(p%eth,idiag_ethmax)
        call sum_mn_name(p%ee,idiag_eem)
        call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_etot /= 0) call sum_mn_name(p%eth + p%ekin, idiag_etot)
        if (idiag_pdivum/=0) call sum_mn_name(p%pp*p%divu,idiag_pdivum)
      endif
!
      if (l1davgfirst) then
        call yzsum_mn_name_x(p%pp,idiag_ppmx)
        call xzsum_mn_name_y(p%pp,idiag_ppmy)
        call xysum_mn_name_z(p%pp,idiag_ppmz)
        call yzsum_mn_name_x(p%TT,idiag_TTmx)
        call xzsum_mn_name_y(p%TT,idiag_TTmy)
        call xysum_mn_name_z(p%TT,idiag_TTmz)
      endif
!
!  2-D averages.
!
      if (l2davgfirst) then
        call zsum_mn_name_xy(p%TT,idiag_TTmxy)
        call ysum_mn_name_xz(p%TT,idiag_TTmxz)
      endif
!
      if (lvideo.and.lfirst) then
        if (ivid_pp/=0) call store_slices(p%pp,pp_xy,pp_xz,pp_yz,pp_xy2,pp_xy3,pp_xy4,pp_xz2)
      endif

    endsubroutine calc_diagnostics_energy
!***********************************************************************
    subroutine energy_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!   1-apr-20/joern: coded
!
      use EquationOfState, only : gamma_m1, gamma
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: cs2
!
!    Slope limited diffusion: update characteristic speed
!    Not staggered yet
!
     if (lslope_limit_diff .and. llast) then
       cs2=0.
       do m=1,my
       do n=1,mz
         if (ldensity_nolog) then
           cs2 = gamma*gamma_m1*f(:,m,n,ieth)/f(:,m,n,irho)
         else
           cs2 = gamma*gamma_m1*f(:,m,n,ieth)*exp(-f(:,m,n,ilnrho))
         endif
         f(:,m,n,isld_char)=f(:,m,n,isld_char)+w_sldchar_ene*sqrt(cs2)
       enddo
       enddo
     endif
!
    endsubroutine energy_before_boundary
!***********************************************************************
    subroutine energy_after_boundary(f)
!
!  Dummy routine.
!
!  04-nov-10/anders+evghenii: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine energy_after_boundary
!***********************************************************************
    subroutine energy_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine energy_after_timestep
!***********************************************************************
    subroutine read_energy_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_init_pars, IOSTAT=iostat)
!
    endsubroutine read_energy_init_pars
!***********************************************************************
    subroutine write_energy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_init_pars)
!
    endsubroutine write_energy_init_pars
!***********************************************************************
    subroutine read_energy_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_run_pars, IOSTAT=iostat)
!
    endsubroutine read_energy_run_pars
!***********************************************************************
    subroutine write_energy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_run_pars)
!
    endsubroutine write_energy_run_pars
!***********************************************************************
    subroutine rprint_energy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
!  04-nov-10/anders+evghenii: adapted from temperature_idealgas.f90
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez, inamexy, inamexz
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_TTm=0; idiag_TTmax=0; idiag_TTmin=0
        idiag_ethm=0; idiag_ethmin=0; idiag_ethmax=0; idiag_eem=0
        idiag_etot = 0; idiag_ethtot=0
        idiag_pdivum=0; idiag_ppm=0
        ivid_pp=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ethmin',idiag_ethmin)
        call parse_name(iname,cname(iname),cform(iname),'ethmax',idiag_ethmax)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'etot',idiag_etot)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'pdivum',idiag_pdivum)
        call parse_name(iname,cname(iname),cform(iname),'ethtot',idiag_ethtot)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ppmx',idiag_ppmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'TTmx',idiag_TTmx)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'ppmy',idiag_ppmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'TTmy',idiag_TTmy)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ppmz',idiag_ppmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTmz',idiag_TTmz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'TTmxy', &
            idiag_TTmxy)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'TTmxz', &
            idiag_TTmxz)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then 
        where(cnamev=='TT'.or.cnamev=='lnTT') cformv='DEFINED'
      endif
      do iname=1,nnamev
        call parse_name(iname,cnamev(iname),cformv(iname),'pp',ivid_pp)
      enddo
!
    endsubroutine rprint_energy
!***********************************************************************
    subroutine get_slices_energy(f,slices)
!
!  04-nov-10/anders+evghenii: adapted
!
      use Slices_methods, only: assign_slices_scal, process_slices, exp2d
      use EquationOfState, only: eoscalc, irho_eth, ilnrho_eth
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      real, dimension(nx) :: penc
      integer :: ieosvars
      integer :: m, n
!
!  Loop over slices
!
      slice_name: select case (trim(slices%name))
!
!  Pressure
!
        case ('pp') slice_name
          call assign_slices_scal(slices,pp_xy,pp_xz,pp_yz,pp_xy2,pp_xy3,pp_xy4,pp_xz2)
!
!  Temperature
!
        case ('TT', 'lnTT') slice_name
          density: if (ldensity_nolog) then
            ieosvars = irho_eth
          else density
            ieosvars = ilnrho_eth
          endif density
!
          do m = m1, m2
            if (lwrite_slice_xy) &
              call eoscalc(ieosvars, f(l1:l2,m,iz_loc,ilnrho), f(l1:l2,m,iz_loc,ieth), iz=iz_loc, lnTT=slices%xy(:,m-nghost))
            if (lwrite_slice_xy2) &
              call eoscalc(ieosvars, f(l1:l2,m,iz2_loc,ilnrho), f(l1:l2,m,iz2_loc,ieth), iz=iz2_loc, lnTT=slices%xy(:,m-nghost))
            if (lwrite_slice_xy3) &
              call eoscalc(ieosvars, f(l1:l2,m,iz3_loc,ilnrho), f(l1:l2,m,iz3_loc,ieth), iz=iz3_loc, lnTT=slices%xy(:,m-nghost))
            if (lwrite_slice_xy4) &
              call eoscalc(ieosvars, f(l1:l2,m,iz4_loc,ilnrho), f(l1:l2,m,iz4_loc,ieth), iz=iz4_loc, lnTT=slices%xy(:,m-nghost))
          enddo
!
          if (lwrite_slice_xz.or.lwrite_slice_xz2) then
            do n = n1, n2
              if (lwrite_slice_xz) &
                call eoscalc(ieosvars, f(l1:l2,iy_loc,n,ilnrho), f(l1:l2,iy_loc,n,ieth), iz=n, lnTT=slices%xz(:,n-nghost))
              if (lwrite_slice_xz2) &
                call eoscalc(ieosvars, f(l1:l2,iy2_loc,n,ilnrho), f(l1:l2,iy2_loc,n,ieth), iz=n, lnTT=slices%xz2(:,n-nghost))
            enddo
          endif
!
          if (lwrite_slice_yz) then
            do n = n1, n2
              do m = m1, m2
                call eoscalc(ieosvars, f(ix_loc,m,n,ilnrho), f(ix_loc,m,n,ieth), iz=n, lnTT=slices%yz(m-nghost,n-nghost))
              enddo
            enddo
          endif
!
          if (trim(slices%name) == 'TT') call process_slices(slices,exp2d)
!
          slices%ready = .true.
!
      endselect slice_name
!
    endsubroutine get_slices_energy
!***********************************************************************
    subroutine fill_farray_pressure(f)
!
!  04-nov-10/anders+evghenii: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine fill_farray_pressure
!***********************************************************************
    subroutine impose_energy_floor(f)
!
!  Trap any negative energy or impose a floor in minimum thermal energy.
!
!  29-aug-11/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      real, dimension(nx) :: eth, rho, eth1
      integer :: i, j, k
      real :: deth
!
!  Impose the energy floor.
!
      efloor: if (energy_floor > 0.0) then
        stratz: if (lstratz) then
          zscan: do k = 1, mz
            deth = energy_floor / eth0z(k) - 1.0
            where(f(:,:,k,ieth) < deth) f(:,:,k,ieth) = deth
          enddo zscan
        else stratz
          where(f(:,:,:,ieth) < energy_floor) f(:,:,:,ieth) = energy_floor
        endif stratz
      endif efloor
!
!  Apply Jeans energy floor.
!
      jfloor: if (ljeans_floor) then
        zscan1: do k = n1, n2
          yscan1: do j = m1, m2
            stratz1: if (lstratz) then
              eth = eth0z(k) * (1.0 + f(l1:l2,j,k,ieth))
              rho = rho0z(k) * (1.0 + f(l1:l2,j,k,irho))
              eth1 = eth
              call jeans_floor(eth, rho)
              where(eth1 /= eth) f(l1:l2,j,k,ieth) = eth / eth0z(k) - 1.0
            else stratz1
              eth = f(l1:l2,j,k,ieth)
              rho = f(l1:l2,j,k,irho)
              if (.not. ldensity_nolog) rho = exp(rho)
              call jeans_floor(eth, rho)
              f(l1:l2,j,k,ieth) = eth
            endif stratz1
          enddo yscan1
        enddo zscan1
      endif jfloor
!
!  Stop the code if negative energy exists.
!
      chkneg: if (lcheck_negative_energy) then
        negeth: if (any(f(l1:l2,m1:m2,n1:n2,ieth) <= merge(-1.0, 0.0, lstratz))) then
          zscan2: do k = n1, n2
            yscan2: do j = m1, m2
              xscan2: do i = l1, l2
                if (lstratz) then
                  if (f(i,j,k,ieth) <= -1.0) print 10, f(i,j,k,ieth), x(i), y(j), z(k)
                else
                  if (f(i,j,k,ieth) <= 0.0) print 20, f(i,j,k,ieth), x(i), y(j), z(k)
                endif
                10 format (1x, 'deth = ', es13.6, ' at x = ', es13.6, ', y = ', es13.6, ', z = ', es13.6)
                20 format (1x, 'eth = ', es13.6, ' at x = ', es13.6, ', y = ', es13.6, ', z = ', es13.6)
              enddo xscan2
            enddo yscan2
          enddo zscan2
          call fatal_error('impose_energy_floor', 'negative energy detected')
        endif negeth
      endif chkneg
!
    endsubroutine impose_energy_floor
!***********************************************************************
    subroutine dynamical_thermal_diffusion(uc)
!
!  Dynamically set thermal diffusion coefficient given fixed mesh Reynolds number.
!
!  02-aug-11/ccyang: coded
!
!  Input Argument
!      uc
!          Characteristic velocity of the system.
!
      real, intent(in) :: uc
!
!  Hyper-diffusion coefficient
!
      if (chi_hyper3_mesh /= 0.0) chi_hyper3_mesh = pi5_1 * uc / re_mesh / sqrt(real(dimensionality))
!
    endsubroutine dynamical_thermal_diffusion
!***********************************************************************
    subroutine split_update_energy(f)
!
!  Update the thermal energy by integrating the operator split energy
!  terms.
!
!  19-jan-13/ccyang: coded
!
      use Boundcond, only: zero_ghosts, update_ghosts
      use Density, only: impose_density_floor
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      integer, dimension(mx,my,mz) :: status
      real, dimension(mx,my,mz) :: delta_eth
      character(len=256) :: message
      integer :: i, j, k
!
      split_update: if (lsplit_update) then
!
!  Impose density and energy floors.
!
        call impose_density_floor(f)
        call impose_energy_floor(f)
!
!  Update the ghost cells.
!
        call zero_ghosts(f, ieth)
        call update_ghosts(f, ieth)
        gcrho: if (ldensity_nolog) then
          call zero_ghosts(f, irho)
          call update_ghosts(f, irho)
        else gcrho
          call zero_ghosts(f, ilnrho)
          call update_ghosts(f, ilnrho)
        endif gcrho
!
!  Update the energy.
!
        deth: if (lstratz) then
          zscan: do k = 1, mz
            call get_delta_eth(real(t), dt, eth0z(k) * (1.0 + f(:,:,k,ieth)), rho0z(k) * (1.0 + f(:,:,k,irho)), &
                               delta_eth(:,:,k), status(:,:,k))
            delta_eth(:,:,k) = delta_eth(:,:,k) / eth0z(k)
          enddo zscan
        elseif (ldensity_nolog) then deth
          call get_delta_eth(real(t), dt, f(:,:,:,ieth), f(:,:,:,irho), delta_eth, status)
        else deth
          call get_delta_eth(real(t), dt, f(:,:,:,ieth), exp(f(:,:,:,ilnrho)), delta_eth, status)
        endif deth
!
        if (any(status < 0)) then
          write(message,10) count(status == -1), ' cells had underflows and ', &
                            count(status == -2), ' cells reached maximum iterations.'
          10 format(1x, i3, a, i3, a)
          call warning('split_update_energy', trim(message))
          do k = n1, n2; do j = m1, m2; do i = l1, l2
            if (status(i,j,k) == -2) print 20, f(i,j,k,irho), get_temperature(f(i,j,k,ieth), f(i,j,k,irho))
            20 format (1x, 'rho = ', es13.6, ' & TT = ', es13.6)
          enddo; enddo; enddo
        endif
!
        if (ldebug) then
          where(status < 0) status = rk_nmax
          print *, 'Minimum, maximum, and average numbers of iterations = ', &
            minval(status), maxval(status), real(sum(status)) / real(nw)
        endif
!
        f(l1:l2,m1:m2,n1:n2,ieth) = f(l1:l2,m1:m2,n1:n2,ieth) + delta_eth(l1:l2,m1:m2,n1:n2)
      endif split_update
!
    endsubroutine
!***********************************************************************
    pure subroutine const_cooling_time(eth, rho, dt)
!
!  Exactly integrate the cooling term dTT/dt = (TTref - TT) / tau_cool.
!  Ideal equation of state is assumed.
!
!  31-jan-13/ccyang: coded.
!
      real, intent(inout) :: eth
      real, intent(in) :: rho, dt
!
      real :: temp, a, b, c, x
!
      c = cv1 / rho
      temp = c * eth
      x = dt / tau_cool
      a = exp(-x)
      b = 1.0 - a
      if (b == 0.0) b = x * (1.0 - 0.5 * x)
      temp = temp * a + TTref * b
      eth = temp / c
!
    endsubroutine const_cooling_time
!***********************************************************************
    elemental subroutine get_delta_eth(t, dt, eth, rho, delta_eth, status)
!
!  Integrate the operator split energy terms that are local:
!
!    de/dt = H(e(t),rho)
!
!  08-aug-11/ccyang: coded
!
      real, intent(in) :: t, dt, eth, rho
      real, intent(out) :: delta_eth
      integer, intent(out), optional :: status
!
!  Cash-Karp parameters for embedded Runga-Kutta method
!
      real, parameter :: safety = 0.9, errcon = (5. / safety)**5
      real, parameter :: a2 = .2, a3 = .3, a4 = .6, a5 = 1., a6 = .875
      real, parameter :: b2 = .2
      real, dimension(2), parameter :: b3 = (/ .075, .225 /)
      real, dimension(3), parameter :: b4 = (/ .3, -.9, 1.2 /)
      real, dimension(4), parameter :: b5 = (/ -11./54., 2.5, -70./27., 35./27. /)
      real, dimension(5), parameter :: b6 = (/ 1631./55296., 175./512., 575./13824., 44275./110592., 253./4096. /)
      real, dimension(6), parameter :: c = (/ 37./378., 0., 250./621., 125./594., 0., 512./1771. /)
      real, dimension(6), parameter :: d = c - (/ 2825./27648., 0., 18575./48384., 13525./55296., 277./14336., .25 /)
!
      real, dimension(6) :: de
      logical :: last, lovershoot
      integer :: nok, nbad, status1
      real :: tf, t1, dt1, eth1, ethj, deth, error
!
!  Initialization
!
      delta_eth = 0.
!
      tf = abs(t + dt)
      t1 = t
      dt1 = dt
      eth1 = eth
      nok = 0
      nbad = 0
      last = .true.
!
      rkck: do
!
!  Embedded Runga-Kutta step
!
        de(1) = dt1 * calc_heat_split(t1, eth1, rho)
        de(2) = dt1 * calc_heat_split(t1 + a2 * dt1, eth1 + b2 * de(1), rho)
        de(3) = dt1 * calc_heat_split(t1 + a3 * dt1, eth1 + sum(b3 * de(:2)), rho)
        de(4) = dt1 * calc_heat_split(t1 + a4 * dt1, eth1 + sum(b4 * de(:3)), rho)
        de(5) = dt1 * calc_heat_split(t1 + a5 * dt1, eth1 + sum(b5 * de(:4)), rho)
        de(6) = dt1 * calc_heat_split(t1 + a6 * dt1, eth1 + sum(b6 * de(:5)), rho)
!
        deth = sum(c * de)
!
!  Check overshoot towards negative energy.
!
        lovershoot = .false.
        floor: if (ljeans_floor) then
          ethj = eth1 + deth
          call jeans_floor(ethj, rho)
          if (ethj > eth1 + deth) then
            lovershoot = .true.
            deth = ethj - eth1
            error = 0.03125
          endif
        elseif (eth1 + deth <= 0.) then floor
          lovershoot = .true.
          error = 32.
        endif floor
!
!  Estimate the error.
!
        if (.not. lovershoot) error = abs(sum(d * de)) / (rk_eps * (eth1 + abs(de(1))))
!
!  Time step control
!
        passed: if (error <= 1.) then
          nok = nok + 1
          delta_eth = delta_eth + deth
!
          done: if (last) then
            if (lovershoot) then
              status1 = 1
            else
              status1 = 0
            endif
            exit rkck
          else if (t1 + dt1 == t1) then done
            status1 = -1
            exit rkck
          endif done
!
          t1 = t1 + dt1
          eth1 = eth + delta_eth
!
          if (error > errcon) then
            dt1 = safety * dt1 * error**(-0.2)
          else
            dt1 = 5. * dt1
          endif
          if (abs(t1 + dt1) >= tf) then
            dt1 = t + dt - t1
            last = .true.
          endif
        else passed
          nbad = nbad + 1
          dt1 = sign(max(abs(safety*dt1*error**(-0.25)), 0.1*abs(dt1)), dt1)
          last = .false.
        endif passed
!
        too_long: if (nok + nbad > rk_nmax) then
          status1 = -2
          exit rkck
        endif too_long
!
      enddo rkck
!
      status_code: if (present(status)) then
        if (ldebug) then
          status = nok + nbad
        else
          status = status1
        endif
      endif status_code
!
!  Cases for which there exists exact solutions.
!
      const: if (lconst_cooling_time) then
        eth1 = eth + delta_eth
        call const_cooling_time(eth1, rho, dt)
        delta_eth = eth1 - eth
      endif const
!
    endsubroutine get_delta_eth
!***********************************************************************
    elemental function calc_heat_split(t, eth, rho) result(heat)
!
!  Calculate the total heat generated from the operator split energy
!  terms.
!
!  06-aug-11/ccyang: coded
!
      real, intent(in) :: t, eth, rho
      real :: heat
!
      real :: eth1
!
      if (t == t) &    ! keep the compiler quiet; can be later removed if t is ever used.
!
      eth1 = eth
      if (ljeans_floor) call jeans_floor(eth1, rho)
!
      heat = 0.
      if (lKI02 .and. lSD93) then
        heat = heat - cool_KI02_SD93(eth1, rho)
      elseif (lKI02) then
        heat = heat + heat_KI02(eth1, rho)
      elseif (lSD93) then
        heat = heat - cool_SD93(eth1, rho)
      endif
!
    endfunction calc_heat_split
!***********************************************************************
    elemental real function heat_KI02(eth, rho)
!
!  Add Koyama & Inutsuka (2002) heating and cooling terms.
!
!  03-oct-11/ccyang: coded
!
      use EquationOfState, only: get_soundspeed
!
      real, intent(in) :: eth, rho
!
      real :: temp, cs2
      real :: c0, c1, c2
!
!  Find the temperature and its derived values.
!
      heat_KI02 = 0.
      temp = get_temperature(eth, rho)
      if (temp <= 0.) return
      c1 = exp(-KI_T1 / (temp + KI_T2))
      c2 = exp(-KI_T3 / temp)
!
!  Calculate the net heat.
!    For 3D run: rho is volume mass density
!    For 2D run: rho is surface mass density
!    cs_z, if > 0, fix the disk vertical scale height
!
      c0 = rho
      if (nzgrid == 1 .and. cs_z <= 0.) then
        call get_soundspeed(real(temp / unit_temperature), cs2)
        c0 = rho / sqrt(cs2)
      endif
      heat_KI02 = rho * (KI_a0 - c0 * (KI_a1 * c1 + KI_a2 * sqrt(temp) * c2))
!
    endfunction
!***********************************************************************
    subroutine init_cooling_SD93()
!
!  Initializes the tabulated cooling function of Sutherland and Dopita (1993).
!
!  02-feb-13/ccyang: coded.
!   8-nov-13/axel: changed some calls to f95 standard
!
      use Mpicomm
      use Syscalls, only: get_env_var
      use EquationOfState, only: getmu
!
      integer, parameter :: lun = 1
      integer, parameter :: ncol = 12
      real, dimension(ncol) :: col
      character(len=256) :: msg, src
      integer :: stat
      real :: mu
      integer :: i
!
!  Find the number of table entries.
!
      get_nt: if (lroot) then
!-      call get_environment_variable('PENCIL_HOME', src)
        call get_env_var('PENCIL_HOME', src)
        src = trim(src) // '/src/cooling_SD93_Z00.dat'
!-      open (unit=lun, file=src, action='read', iostat=stat, iomsg=msg)
        open (unit=lun, file=src, action='read', iostat=stat)
        if (stat /= 0) call fatal_error('init_cooling_SD93', 'cannot open the cooling table; ' // trim(msg), force=.true.)
        nline: do
!-        read (lun,*,iostat=stat,iomsg=msg) col
          read (lun,*,iostat=stat) col
          if (stat < 0) exit nline
          if (stat > 0) call fatal_error('init_cooling_SD93', 'error in reading the cooling table; ' // trim(msg), force=.true.)
          SD_nt = SD_nt + 1
        enddo nline
!-      close (unit=lun, iostat=stat, iomsg=msg)
        close (unit=lun, iostat=stat)
        if (stat /= 0) call fatal_error('init_cooling_SD93', 'cannot close the cooling table; ' // trim(msg), force=.true.)
      endif get_nt
      call mpibcast_int(SD_nt,comm=MPI_COMM_WORLD)
!
!  Allocate the cooling table.
!
      allocate (SD_logTT(SD_nt), SD_logLambda(SD_nt), stat=stat)
      if (stat /= 0) call fatal_error('init_cooling_SD93', 'cannot allocate the cooling table. ')
!
!  Read the cooling table.
!
      table: if (lroot) then
        open (unit=lun, file=src, action='read')
        reading: do i = 1, SD_nt
          read (lun,*) col
          SD_logTT(i) = col(1)
          SD_logLambda(i) = col(5)
        enddo reading
        close (unit=lun)
      endif table
      call mpibcast_real(SD_logTT, SD_nt, comm=MPI_COMM_WORLD)
      call mpibcast_real(SD_logLambda, SD_nt,comm=MPI_COMM_WORLD)
!
      if (lroot) print *, 'init_cooling_SD93: read the Sutherland & Dopita (1993) cooling table. '
!
!  Save the conversion factor.
!
      call getmu(mu_tmp=mu)
      SD_a0 = exp(log(unit_time) - log(unit_energy) - 3. * log(unit_length) - 2. * log(mu * m_u))
      if (unit_system == 'SI') SD_a0 = 1.e-13 * SD_a0
      disk: if (nzgrid == 1) then
        SD_a0 = 0.5 * SD_a0 * nu_z / sqrtpi
        if (cs_z > 0.) SD_a0 = SD_a0 / cs_z
      endif disk
!
    endsubroutine init_cooling_SD93
!***********************************************************************
    elemental real function cool_SD93(eth, rho)
!
!  Add Sutherland & Dopita (1993) cooling function.
!
!  02-feb-13/ccyang: coded.
!
      use General, only: spline
      use EquationOfState, only: get_soundspeed
!
      real, intent(in) :: eth, rho
!
      logical :: err
      real, dimension(1) :: logLambda
      real :: temp, cs2
!
!  Find the temperature.
!
      cool_SD93 = 0.
      temp = get_temperature(eth, rho)
      if (temp <= 0.) return
!
!  Interpolate the cooling table.
!
      call spline(SD_logTT, SD_logLambda, (/log10(temp)/), logLambda, SD_nt, 1, err)
      if (err) return
!
!  Calculate the net cooling rate.
!    For 3D run: rho is volume mass density
!    For 2D run: rho is surface mass density
!    cs_z, if > 0, fix the disk vertical scale height
!
      cool_SD93 = SD_a0 * rho**2 * 10.**logLambda(1)
      disk: if (nzgrid == 1 .and. cs_z <= 0.) then
        call get_soundspeed(real(temp / unit_temperature), cs2)
        cool_SD93 = cool_SD93 / sqrt(cs2)
      endif disk
!
    endfunction cool_SD93
!***********************************************************************
    elemental real function cool_KI02_SD93(eth, rho) result(cool)
!
!  Smoothly connects the KI02 and the SD93 heating/cooling functions.
!
!  04-feb-13/ccyang: coded.
!
      real, intent(in) :: eth, rho
!
      real, parameter :: logTTmin = 4., logTTmax = 4.2
      real, parameter :: dlogTT = logTTmax - logTTmin
      real :: temp, logTT, cool1, cool2
!
!  Find the temperature.
!
      temp = get_temperature(eth, rho)
      error: if (temp <= 0.) then
        cool = 0.
        return
      endif error
      logTT = log10(temp)
!
!  Linearly combine the two heating/cooling functions.
!
      step: if (logTT < logTTmin) then
        cool = -heat_KI02(eth, rho)
      elseif (logTT > logTTmax) then step
        cool = cool_SD93(eth, rho)
      else step
        cool1 = -heat_KI02(eth, rho)
        cool2 = cool_SD93(eth, rho)
        temp = cos(pi * (logTT - logTTmin) / dlogTT)
        cool = 0.5 * ((1. + temp) * cool1 + (1. - temp) * cool2)
      endif step
!
    endfunction cool_KI02_SD93
!***********************************************************************
    elemental subroutine jeans_floor(eth, rho)
!
!  Apply a floor to energy where Jeans length is unresolved.
!
!  29-aug-11/ccyang: coded.
!
      real, intent(inout) :: eth
      real, intent(in) :: rho
!
      real :: eth_min
!
!  Find the Jeans energy floor.
!
      if (nzgrid /= 1) then
!       3D run: rho is volume mass density
      else
!       2D run: rho is surface mass density
        eth_min = Jeans_c0 * rho**2
      endif
!
!  Impose the floor where needed.
!
      if (eth < eth_min) eth = eth_min
!
    endsubroutine
!***********************************************************************
    subroutine expand_shands_energy()
!
!  Presently dummy, for possible use
!
    endsubroutine expand_shands_energy
!***********************************************************************
    pure real function get_temperature(eth, rho) result(temp)
!
!  Finds the absolute temperature in Kelvins assuming the ideal-gas EOS.
!
!  09-sep-14/ccyang: coded.
!
      real, intent(in) :: eth, rho
!
!  Find the temperature.
!
      if (lstratz) then
        temp = cv1_temp * (1.0 + eth) / (1.0 + rho)
      else
        temp = cv1_temp * eth / rho
      endif
!
    endfunction get_temperature
!***********************************************************************
    subroutine update_char_vel_energy(f)
!
! TB implemented.
!
!   25-sep-15/MR+joern: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f

      call keep_compiler_quiet(f)

      call warning('update_char_vel_energy', &
           'characteristic velocity not yet implemented for thermal energy')

    endsubroutine update_char_vel_energy
!***********************************************************************
    subroutine pushdiags2c(p_diag)

    integer, parameter :: n_diags=0
    integer(KIND=ikind8), dimension(:) :: p_diag

    call keep_compiler_quiet(p_diag)

    endsubroutine pushdiags2c
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
 
    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(chi,p_par(1))

    endsubroutine pushpars2c
!***********************************************************************
endmodule Energy
