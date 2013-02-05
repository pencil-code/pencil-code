! $Id$
!
!  This module can replace the entropy module by using the thermal energy
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
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED Ma2; fpres(3); transpeth
!
!***************************************************************
module Entropy
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'entropy.h'
!
  real :: eth_left, eth_right, widtheth, eth_const=1.0
  real :: chi=0.0, chi_shock=0.0, chi_shock_gradTT=0., chi_hyper3_mesh=0.
  real :: energy_floor = 0., temperature_floor = 0.
  real :: rk_eps = 1.e-3
  real :: nu_z=0., cs_z=0.
  real :: detonation_factor = 1.
  real :: TTref = 0., tau_cool = 1.
  integer :: rk_nmax = 100
  integer :: njeans = 4
  integer :: idet = 0
  integer :: detonation_power = 0
  logical :: lsplit_update=.false.
  logical :: lviscosity_heat=.true.
  logical :: lupw_eth=.false.
  logical :: lcheck_negative_energy=.false.
  logical :: lKI02 = .false.
  logical :: lSD93 = .false.
  logical :: lconst_cooling_time = .false.
  logical :: ljeans_floor=.false.
  logical :: ldetonate=.false.
  logical, pointer :: lpressuregradient_gas
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
      chi, chi_shock, chi_shock_gradTT, chi_hyper3_mesh, &
      energy_floor, temperature_floor, lcheck_negative_energy, &
      rk_eps, rk_nmax, lKI02, lSD93, nu_z, cs_z, &
      lconst_cooling_time, TTref, tau_cool, &
      ljeans_floor, njeans, &
      ldetonate, detonation_factor, feedback
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
  integer :: idiag_ethmin=0   ! DIAG_DOC: $\mathrm{min} e_\text{th}$
  integer :: idiag_ethmax=0   ! DIAG_DOC: $\mathrm{max} e_\text{th}$
  integer :: idiag_eem=0      ! DIAG_DOC: $\left< e \right> =
                              ! DIAG_DOC:  \left< c_v T \right>$
                              ! DIAG_DOC: \quad(mean internal energy)
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
  real, dimension(nx,nz) :: pp_xz
  real, dimension(ny,nz) :: pp_yz
  real, dimension(nx,ny) :: pp_xy,pp_xy2,pp_xy3,pp_xy4
!
! General variables for operator split terms.
!
  real :: cv1 = 0., cv1_temp = 0.
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
! Variables for the detonations.
!
  logical, dimension(-nghost:nghost,-nghost:nghost,-nghost:nghost) :: mask_sphere = .false.
  real, dimension(-nghost:nghost,-nghost:nghost,-nghost:nghost) :: smooth = 0.
  integer :: nxs = 0, nys = 0, nzs = 0, nyz = ny * nz
  real :: deposit = 0.
!
  contains
!***********************************************************************
    subroutine register_entropy()
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
      call farray_register_auxiliary('det', idet, communicated=.true.)
!
!  logical variable lpressuregradient_gas shared with hydro modules
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas,ierr)
      if (ierr/=0) call fatal_error('register_entropy','lpressuregradient_gas')
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  04-nov-10/anders+evghenii: adapted
!  03-oct-11/ccyang: add initialization for KI02
!
      use EquationOfState, only: select_eos_variable, get_cv1, getmu, gamma, gamma_m1, cs0, cs20
      use SharedVariables
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      logical, intent (in) :: lstarting
!
      integer :: i, j, k
      real :: mu
      real :: c0
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
      call select_eos_variable('eth',ieth)
!
      call put_shared_variable('lviscosity_heat',lviscosity_heat)
!
!  Decide if operator splitting is required.
!
      lsplit_update = lconst_cooling_time .or. lKI02 .or. lSD93
      if (lsplit_update .and. .not. ldensity) call fatal_error('initialize_entropy', 'Density is required for split_update_energy.')
!
!  General variables required by split_update_energy.
!
      ideal_gas: if (lsplit_update) then
        if (.not. leos_idealgas) call fatal_error('initialize_entropy', 'currently assumes eos_idealgas')
        call get_cv1(cv1)
        cv1_temp = cv1 * unit_temperature
      endif ideal_gas
!
!  Initialize the KI02 terms.
!
      KI02: if (lKI02) then
        call getmu(mu_tmp=mu)
        KI_a0 = (KI_heat * unit_time / unit_energy) / (mu * m_u)
        c0 = KI_a0 / (mu * m_u * unit_length**3)
        SI: if (unit_system == 'SI') then
          KI_a0 = 1.e-7 * KI_a0
          c0 = 1.e-13 * c0
        endif SI
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
          call fatal_error('initialize_entropy', '3D Jeans floor under construction')
        else
          Jeans_c0 = real(njeans) * G_Newton * dxmax / (gamma * gamma_m1)
        endif
      endif Jeans
!
!  Initialize the variables for detonations.
!
      detonate: if (ldetonate) then
!       Determine the amount of deposit into each cell to be detonated.
        method: select case (feedback)
        case ('linear') method
          deposit = detonation_factor * cs20
          detonation_power = 1
        case ('binding') method
          deposit = detonation_factor * Jeans_c0
          detonation_power = 2
        case default method
          call fatal_error('detonate', 'unknown feedback method')
        endselect method
!       Smoothing kernal
        smooth(0,0,0) = 1.
        xdir: if (nxgrid > 1) then
          nxs = nghost
          forall (i=-nghost:nghost, i/=0) smooth(i,0,0) = smooth(0,0,0) * exp(-0.5 * real(i * i))
        endif xdir
        ydir: if (nygrid > 1) then
          nys = nghost
          forall (i=-nghost:nghost, j=-nghost:nghost, j/=0) smooth(i,j,0) = smooth(i,0,0) * exp(-0.5 * real(j * j))
        endif ydir
        zdir: if (nzgrid > 1) then
          nzs = nghost
          forall (i=-nghost:nghost, j=-nghost:nghost, k=-nghost:nghost, k/=0) &
            smooth(i,j,k) = smooth(i,j,0) * exp(-0.5 * real(k * k))
        endif zdir
        smooth = smooth / sum(smooth)
!       Spherical mask
        forall (i=-nxs:nxs, j=-nys:nys, k=-nzs:nzs, i*i+j*j+k*k <= 2**2) mask_sphere(i,j,k) = .true.
      endif detonate
!
      if (llocal_iso) &
           call fatal_error('initialize_entropy', &
           'llocal_iso switches on the local isothermal approximation. ' // &
           'Use ENTROPY=noentropy in src/Makefile.local')
!
    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f)
!
!  Initialise thermal energy.
!
!  04-nov-10/anders+evghenii: adapted
!
      use General, only: itoa
      use Initcond, only: jump
      use EquationOfState, only: rho0, cs20, gamma, gamma_m1
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
            if (gamma_m1 /= 0.) then
              f(:,:,:,ieth) = f(:,:,:,ieth) + rho0 * cs20 / (gamma * gamma_m1)
            else
              f(:,:,:,ieth) = f(:,:,:,ieth) + rho0 * cs20
            endif
!
          case default
!
!  Catch unknown values.
!
            write(unit=errormsg,fmt=*) 'No such value for initeth(' &
                //trim(iinit_str)//'): ',trim(initeth(j))
            call fatal_error('init_ss',errormsg)
!
          endselect
        endif
      enddo
!
    endsubroutine init_ss
!***********************************************************************
    subroutine pencil_criteria_entropy()
!
!  All pencils that the Entropy module depends on are specified here.
!
!  04-nov-10/anders+evghenii: adapted
!
      lpenc_requested(i_pp)=.true.
      lpenc_requested(i_divu)=.true.
      if (lweno_transport) then
        lpenc_requested(i_transpeth)=.true.
      else
        lpenc_requested(i_geth)=.true.
        lpenc_requested(i_eth)=.true.
      endif
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
      if (chi_shock /= 0.) then
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_geth)=.true.
        lpenc_requested(i_del2eth)=.true.
      endif
      if (chi_shock_gradTT/=0.0) then
        lpenc_requested(i_cp)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_del2TT)=.true.
        lpenc_requested(i_grho)=.true.
        lpenc_requested(i_gTT)=.true.
        lpenc_requested(i_gshock)=.true.
      endif
!
!  Diagnostic pencils.
!
      if (idiag_ethm/=0 .or. idiag_ethmin/=0 .or. idiag_ethmax/=0) lpenc_diagnos(i_eth)=.true.
      if (idiag_eem/=0) lpenc_diagnos(i_ee)=.true.
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
    endsubroutine pencil_criteria_entropy
!***********************************************************************
    subroutine pencil_interdep_entropy(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  04-nov-10/anders+evghenii: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_fpres)) then
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_geth)=.true.
      endif
!
    endsubroutine pencil_interdep_entropy
!***********************************************************************
    subroutine calc_pencils_entropy(f,p)
!
!  Calculate Entropy pencils.
!  This routine is called after  calc_pencils_eos
!  Most basic pencils should come first, as others may depend on them.
!
!  04-nov-10/anders+evghenii: adapted
!  14-feb-11/bing: moved eth dependend pecnils to eos_idealgas
!
      use EquationOfState, only: gamma_m1
      use Sub, only: del2
      use WENO_transport, only: weno_transp
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      integer :: j
! fpres
      if (lpencil(i_fpres)) then
        do j=1,3
          p%fpres(:,j)=-p%rho1*gamma_m1*p%geth(:,j)
        enddo
      endif
! transpeth
      if (lpencil(i_transpeth)) &
          call weno_transp(f,m,n,ieth,-1,iux,iuy,iuz,p%transpeth,dx_1,dy_1,dz_1)
!
    endsubroutine calc_pencils_entropy
!***********************************************************************
    subroutine dss_dt(f,df,p)
!
!  Calculate right hand side of energy equation.
!
!    deth/dt + div(u*eth) = -P*div(u)
!
!  04-nov-10/anders+evghenii: coded
!  02-aug-11/ccyang: add mesh hyper-diffusion
!
      use Diagnostics
      use EquationOfState, only: gamma
      use Special, only: special_calc_entropy
      use Sub, only: identify_bcs, u_dot_grad
      use Viscosity, only: calc_viscous_heat
      use Deriv, only: der6
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax=0.0,ugeth
      real, dimension(nx,3) :: uu
      real, dimension(nx) :: d6eth
      integer :: j
!
      intent(inout) :: f,p
      intent(out) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*, 'dss_dt: solve deth_dt'
      if (headtt) call identify_bcs('eth',ieth)
!
!  Sound speed squared.
!
      if (headtt) print*, 'dss_dt: cs20=', p%cs2(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  Add pressure gradient term in momentum equation.
!
      if (lhydro) df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%fpres
!
!  Entry possibility for "personal" entries.
!  In that case you'd need to provide your own "special" routine.
!
      if (lspecial) call special_calc_entropy(f,df,p)
!
!  Add energy transport term.
!
      if (lweno_transport) then
        if (lconst_advection) call fatal_error('dss_dt', 'constant background advection is not implemented with WENO transport.')
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) - p%transpeth
      else
        uu = p%uu
        if (lconst_advection) uu = uu + spread(u0_advec,1,nx)
        call u_dot_grad(f, ieth, p%geth, uu, ugeth, UPWIND=lupw_eth)
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) - p%eth*p%divu - ugeth
      endif
!
!  Add P*dV work.
!
      df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) - p%pp*p%divu
!
!  Calculate viscous contribution to temperature.
!
      if (lviscosity.and.lviscosity_heat) call calc_viscous_heat(f,df,p,Hmax)
!
!  Thermal energy diffusion.
!
      if (chi/=0.0) then
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + chi * p%cp * (p%rho*p%del2TT + sum(p%grho*p%gTT, 2))
        if (lfirst .and. ldt) diffus_chi = diffus_chi + gamma*chi*dxyz_2
      endif
!
!  Shock diffusion
!
      if (chi_shock /= 0.) then
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + chi_shock * (p%shock*p%del2eth + sum(p%gshock*p%geth, 2))
        if (lfirst .and. ldt) diffus_chi = diffus_chi + chi_shock*p%shock*dxyz_2
      endif
!
      if (chi_shock_gradTT/=0.0) then
        df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + chi_shock_gradTT * p%cp * ( &
            p%shock * (p%rho*p%del2TT + sum(p%grho*p%gTT, 2)) + p%rho * sum(p%gshock*p%gTT, 2))
        if (lfirst .and. ldt) diffus_chi = diffus_chi + gamma*chi_shock_gradTT*p%shock*dxyz_2
      endif
!
!  Mesh hyper-diffusion
!
      if (chi_hyper3_mesh/=0.0) then
        do j=1,3
          call der6(f, ieth, d6eth, j, IGNOREDX=.true.)
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + &
              chi_hyper3_mesh * d6eth * dline_1(:,j)
        enddo
        if (lfirst .and. ldt) diffus_chi3 = diffus_chi3 + chi_hyper3_mesh* &
            (abs(dline_1(:,1))+abs(dline_1(:,2))+abs(dline_1(:,3)))
      endif
!
!  Diagnostics.
!
      if (ldiagnos) then
        if (idiag_TTm/=0)    call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_TTmax/=0)  call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0)  call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        if (idiag_ethm/=0)   call sum_mn_name(p%eth,idiag_ethm)
        if (idiag_ethmin/=0) call max_mn_name(-p%eth,idiag_ethmin,lneg=.true.)
        if (idiag_ethmax/=0) call max_mn_name(p%eth,idiag_ethmax)
        if (idiag_eem/=0)    call sum_mn_name(p%ee,idiag_eem)
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
        if (idiag_TTmxy/=0) call zsum_mn_name_xy(p%TT,idiag_TTmxy)
        if (idiag_TTmxz/=0) call ysum_mn_name_xz(p%TT,idiag_TTmxz)
      endif
!
      if (lvideo.and.lfirst) then
        pp_yz(m-m1+1,n-n1+1)=p%pp(ix_loc-l1+1)
        if (m==iy_loc)  pp_xz(:,n-n1+1)=p%pp
        if (n==iz_loc)  pp_xy(:,m-m1+1)=p%pp
        if (n==iz2_loc) pp_xy2(:,m-m1+1)=p%pp
        if (n==iz3_loc) pp_xy3(:,m-m1+1)=p%pp
        if (n==iz4_loc) pp_xy4(:,m-m1+1)=p%pp
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_lentropy_pars(f)
!
!  Dummy routine.
!
!  04-nov-10/anders+evghenii: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lentropy_pars
!***********************************************************************
    subroutine read_entropy_init_pars(unit,iostat)
!
!  04-nov-10/anders+evghenii: coded
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
!  04-nov-10/anders+evghenii: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_init_pars)
!
    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(unit,iostat)
!
!  04-nov-10/anders+evghenii: coded
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
!  04-nov-10/anders+evghenii: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_run_pars)
!
    endsubroutine write_entropy_run_pars
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
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
        idiag_pdivum=0; idiag_ppm=0
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
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'pdivum',idiag_pdivum)
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
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine get_slices_entropy(f,slices)
!
!  04-nov-10/anders+evghenii: adapted
!
      use EquationOfState, only: eoscalc, irho_eth, ilnrho_eth
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      real, dimension(nx) :: penc
      integer :: ieosvars, idensity
      integer :: m, n
      real :: a
!
!  Loop over slices
!
      slice_name: select case (trim(slices%name))
!
!  Pressure
!
        case ('pp') slice_name
          slices%yz =pp_yz
          slices%xz =pp_xz
          slices%xy =pp_xy
          slices%xy2=pp_xy2
          if (lwrite_slice_xy3) slices%xy3=pp_xy3
          if (lwrite_slice_xy4) slices%xy4=pp_xy4
          slices%ready=.true.
!
!  Temperature
!
        case ('TT', 'lnTT') slice_name
          density: if (ldensity_nolog) then
            ieosvars = irho_eth
            idensity = irho
          else density
            ieosvars = ilnrho_eth
            idensity = ilnrho
          endif density
!
          xy: do m = m1, m2
            call eoscalc(ieosvars, f(l1:l2,m,iz_loc,idensity), f(l1:l2,m,iz_loc,ieth), lnTT=penc)
            slices%xy(:,m-nghost) = penc
            call eoscalc(ieosvars, f(l1:l2,m,iz2_loc,idensity), f(l1:l2,m,iz2_loc,ieth), lnTT=penc)
            slices%xy2(:,m-nghost) = penc
            call eoscalc(ieosvars, f(l1:l2,m,iz3_loc,idensity), f(l1:l2,m,iz3_loc,ieth), lnTT=penc)
            slices%xy3(:,m-nghost) = penc
            call eoscalc(ieosvars, f(l1:l2,m,iz4_loc,idensity), f(l1:l2,m,iz4_loc,ieth), lnTT=penc)
            slices%xy4(:,m-nghost) = penc
          enddo xy
!
          xz: do n = n1, n2
            call eoscalc(ieosvars, f(l1:l2,iy_loc,n,idensity), f(l1:l2,iy_loc,n,ieth), lnTT=penc)
            slices%xz(:,n-nghost) = penc
          enddo xz
!
          yz_z: do n = n1, n2
            yz_y: do m = m1, m2
              call eoscalc(ieosvars, f(ix_loc,m,n,idensity), f(ix_loc,m,n,ieth), lnTT=a)
              slices%yz(m-nghost,n-nghost) = a
            enddo yz_y
          enddo yz_z
!
          nolog: if (trim(slices%name) == 'TT') then
            slices%xy = exp(slices%xy)
            slices%xy2 = exp(slices%xy2)
            slices%xy3 = exp(slices%xy3)
            slices%xy4 = exp(slices%xy4)
            slices%xz = exp(slices%xz)
            slices%yz = exp(slices%yz)
          endif nolog
!
          slices%ready = .true.
!
      endselect slice_name
!
    endsubroutine get_slices_entropy
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
      integer :: i, j, k
      real, dimension(nx) :: eth, rho
!
!  Impose the energy floor.
!
      if (energy_floor > 0.) where(f(:,:,:,ieth) < energy_floor) f(:,:,:,ieth) = energy_floor
!
!  Apply Jeans energy floor.
!
      if (ljeans_floor) then
        do k = n1, n2
          do j = m1, m2
            eth = f(l1:l2,j,k,ieth)
            rho = f(l1:l2,j,k,irho)
            if (.not. ldensity_nolog) rho = exp(rho)
            call jeans_floor(eth, rho)
            f(l1:l2,j,k,ieth) = eth
          enddo
        enddo
      endif
!
!  Stop the code if negative energy exists.
!
      if (lcheck_negative_energy) then
        if (any(f(l1:l2,m1:m2,n1:n2,ieth) <= 0.)) then
          do k = n1, n2
            do j = m1, m2
              do i = l1, l2
                if (f(i,j,k,ieth) <= 0.) print 10, f(i,j,k,ieth), x(i), y(j), z(k)
                10 format (1x, 'eth = ', es13.6, ' at x = ', es13.6, ', y = ', es13.6, ', z = ', es13.6)
              enddo
            enddo
          enddo
          call fatal_error('impose_energy_floor', 'negative energy detected')
        endif
      endif
!
    endsubroutine impose_energy_floor
!***********************************************************************
    subroutine dynamical_thermal_diffusion(umax)
!
!  Dynamically set thermal diffusion coefficient given fixed mesh Reynolds number.
!
!  02-aug-11/ccyang: coded
!
      real, intent(in) :: umax
!
!  Hyper-diffusion coefficient
!
      if (chi_hyper3_mesh /= 0.) chi_hyper3_mesh = pi5_1 * umax / re_mesh
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
      use Boundcond, only: update_ghosts
      use Density, only: impose_density_floor
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      integer, dimension(mx,my,mz) :: status
      real, dimension(mx,my,mz) :: delta_eth
      character(len=256) :: message
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
        call update_ghosts(f, ieth, ieth)
        if (ldensity_nolog) then
            call update_ghosts(f, irho, irho)
        else
            call update_ghosts(f, ilnrho, ilnrho)
        endif
!
!  Update the energy.
!
        if (ldensity_nolog) then
          call get_delta_eth(real(t), dt, f(:,:,:,ieth), f(:,:,:,irho), delta_eth, status)
        else
          call get_delta_eth(real(t), dt, f(:,:,:,ieth), exp(f(:,:,:,ilnrho)), delta_eth, status)
        endif
!
        if (any(status < 0)) then
          write(message,10) count(status == -1), ' cells had underflows and ', &
                            count(status == -2), ' cells reached maximum iterations.'
          10 format(1x, i3, a, i3, a)
          call warning('split_update_energy', trim(message))
        endif
!
        if (ldebug) then
          where(status < 0) status = rk_nmax
          print *, 'Minimum, maximum, and average numbers of iterations = ', &
            minval(status), maxval(status), real(sum(status)) / real(nw)
        endif
!
        f(l1:l2,m1:m2,n1:n2,ieth) = f(l1:l2,m1:m2,n1:n2,ieth) + delta_eth(l1:l2,m1:m2,n1:n2)
!
!  Detonate those cells violating the Jeans criterion.
!
        if (ldetonate) call detonate(f, status == 1)
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
      hot: if (temp > TTref) then
        x = dt / tau_cool
        a = exp(-x)
        b = 1. - a
        if (b == 0.) b = x * (1. - 0.5 * x)
        temp = temp * a + TTref * b
        eth = temp / c
      endif hot
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
      real, dimension(2), parameter :: b3 = [.075, .225]
      real, dimension(3), parameter :: b4 = [.3, -.9, 1.2]
      real, dimension(4), parameter :: b5 = [-11./54., 2.5, -70./27., 35./27.]
      real, dimension(5), parameter :: b6 = [1631./55296., 175./512., 575./13824., 44275./110592., 253./4096.]
      real, dimension(6), parameter :: c = [37./378., 0., 250./621., 125./594., 0., 512./1771.]
      real, dimension(6), parameter :: d = c - [2825./27648., 0., 18575./48384., 13525./55296., 277./14336., .25]
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
!
      use Mpicomm
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
        call get_environment_variable('PENCIL_HOME', src)
        src = trim(src) // '/src/cooling_SD93_Z00.dat'
        open (unit=lun, file=src, action='read', iostat=stat, iomsg=msg)
        if (stat /= 0) call fatal_error('init_cooling_SD93', 'cannot open the cooling table; ' // trim(msg), force=.true.)
        nline: do
          read (lun,*,iostat=stat,iomsg=msg) col
          if (stat < 0) exit nline
          if (stat > 0) call fatal_error('init_cooling_SD93', 'error in reading the cooling table; ' // trim(msg), force=.true.)
          SD_nt = SD_nt + 1
        enddo nline
        close (unit=lun, iostat=stat, iomsg=msg)
        if (stat /= 0) call fatal_error('init_cooling_SD93', 'cannot close the cooling table; ' // trim(msg), force=.true.)
      endif get_nt
      call mpibcast_int(SD_nt)
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
      call mpibcast_real(SD_logTT, SD_nt)
      call mpibcast_real(SD_logLambda, SD_nt)
!
      if (lroot) print *, 'init_cooling_SD93: read the Sutherland & Dopita (1993) cooling table. '
!
!  Save the conversion factor.
!
      call getmu(mu_tmp=mu)
      SD_a0 = unit_time / (unit_energy * unit_length**3 * (mu * m_u)**2)
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
        cool = exp(0.5 * ((1. + temp) * log(cool1) + (1. - temp) * log(cool2)))
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
    subroutine detonate(f, unstable)
!
!  Detonate specified cells by specified method.
!
!  20-jan-13/ccyang: coded.
!
      use Mpicomm
      use Boundcond
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      logical, dimension(mx,my,mz), intent(in) :: unstable
!
      real, dimension(nx,ny,nz) :: delta
      real, dimension(nx) :: ps
      real, dimension(-3:3) :: uu
      integer :: ndet
      integer :: i, j, k, imn, m, n
      integer :: ll1, ll2, mm1, mm2, nn1, nn2
      real :: divu, r
!
!  Prepare for communicating the cells that should be detonated.
!
      ndet = 0
      f(:,:,:,idet) = 0.
      zscan: do k = n1, n2
        nn1 = k - nzs
        nn2 = k + nzs
        yscan: do j = m1, m2
          mm1 = j - nys
          mm2 = j + nys
          xscan: do i = l1, l2
            ll1 = i - nxs
            ll2 = i + nxs
            trigger: if (maxval(f(ll1:ll2,mm1:mm2,nn1:nn2,irho), &
                                mask=mask_sphere(-nxs:nxs,-nys:nys,-nzs:nzs)) == f(i,j,k,irho) .and. &
                         all(unstable(ll1:ll2,mm1:mm2,nn1:nn2) .and. mask_sphere(-nxs:nxs,-nys:nys,-nzs:nzs) &
                                                               .eqv. mask_sphere(-nxs:nxs,-nys:nys,-nzs:nzs))) then
              divu = 0.
              duxdx: if (nxs > 0) then
                uu = f(ll1:ll2,j,k,iux)
                divu = divu + derivative(uu, dx_1(i))
              endif duxdx
              duydy: if (nys > 0) then
                uu = f(i,mm1:mm2,k,iuy)
                divu = divu + derivative(uu, dy_1(j))
              endif duydy
              duzdz: if (nzs > 0) then
                uu = f(i,j,nn1:nn2,iuz)
                divu = divu + derivative(uu, dz_1(k))
              endif duzdz
              converge: if (divu < 0.) then
                ndet = ndet + 1
                f(i,j,k,idet) = deposit * f(i,j,k,irho)**detonation_power
              endif converge
            endif trigger
          enddo xscan
        enddo yscan
      enddo zscan
      if (ldebug .and. ndet > 0) print *, 'Detonated ', ndet, ' cells. '
!
!  Smooth out the detonation energy.
!
      shear: if (lshear) then
        call boundconds_y(f, idet, idet)
        call initiate_isendrcv_bdry(f, idet, idet)
        call finalize_isendrcv_bdry(f, idet, idet)
      endif shear
      call boundconds_x(f, idet, idet)
      call initiate_isendrcv_bdry(f, idet, idet)
!
      pencil: do imn = 1, nyz
        n = nn(imn)
        m = mm(imn)
        final: if (necessary(imn)) then
          call finalize_isendrcv_bdry(f, idet, idet)
          call boundconds_y(f, idet, idet)
          call boundconds_z(f, idet, idet)
        endif final
!
        ps = 0.
        zdir: do k = -nzs, nzs
          ydir: do j = -nys, nys
            xdir: do i = -nxs, nxs
              ps = ps + smooth(i,j,k) * f(l1+i:l2+i, m+j, n+k, idet)
            enddo xdir
          enddo ydir
        enddo zdir
        delta(:,m-nghost,n-nghost) = ps
      enddo pencil
!
!  Make the deposit into the cells.
!
      f(l1:l2,m1:m2,n1:n2,ieth) = f(l1:l2,m1:m2,n1:n2,ieth) + delta(:,:,:)
!
    endsubroutine detonate
!***********************************************************************
    subroutine expand_shands_entropy()
!
!  Presently dummy, for possible use
!
    endsubroutine expand_shands_entropy
!***********************************************************************
    real function derivative(a, dx1)
!
!  Returns the first derivative at the center of a seven-point data
!  array.
!
!  19-jan-13/ccyang: coded.
!
      real, dimension(-3:3), intent(in) :: a
      real, intent(in) :: dx1
!
      real, parameter :: sixtieth = 1. / 60.
!
      derivative = sixtieth * dx1 * ((a(3) - a(-3)) - 9.0 * (a(2) - a(-2)) + 45.0 * (a(1) - a(-1)))
!
    endfunction
!***********************************************************************
    pure real function get_temperature(eth, rho) result(temp)
!
!  Finds the absolute temperature in Kelvins assuming the ideal-gas EOS.
!
!  04-feb-13/ccyang: coded.
!
      real, intent(in) :: eth, rho
!
!  Find the temperature.
!
      if (cv1_temp > 0. .and. eth > 0. .and. rho > 0.) then
        temp = cv1_temp * eth / rho
      else
        temp = 0.
      endif
!
    endfunction get_temperature
!***********************************************************************
endmodule Entropy
