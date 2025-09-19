! $Id$
!
!  Electric field, dE/dt = curlB, originally only for the special case
!  of no fluid induction, but now fluid motions are also included.
!
!  25-feb-07/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 4
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED e2; edot2; el(3); a0; ga0(3); del2ee(3); curlE(3); BcurlE
! PENCILS PROVIDED rhoe, divJ, divE, gGamma(3); sigE, sigB; eb; count_eb0
! PENCILS PROVIDED boost; gam_EB; eprime; bprime; jprime; GammaY
! PENCILS PROVIDED jj_higgsY(3); rhoe_higgsY
! PENCILS EXPECTED phi, infl_phi, dphi, infl_dphi, gphi(3); cov_der(4,4)
! PENCILS EXPECTED curlb(3), jj_ohm(3), phi_doublet(3)
! PENCILS EXPECTED gpsi(3), dpsi
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
! input parameters
!
  real, dimension (ninit) :: amplee=0.0 !, kx_aa=1.0, ky_aa=1.0, kz_aa=1.0
  real :: alpf=0., alpfpsi=0.
  real :: ampl_ex=0.0, ampl_ey=0.0, ampl_ez=0.0, ampl_a0=0.0
  real :: kx_ex=0.0, kx_ey=0.0, kx_ez=0.0
  real :: ky_ex=0.0, ky_ey=0.0, ky_ez=0.0
  real :: kz_ex=0.0, kz_ey=0.0, kz_ez=0.0
  real :: kx_a0=0.0, ky_a0=0.0, kz_a0=0.0
  real :: phase_ex=0.0, phase_ey=0.0, phase_ez=0.0, phase_a0=0.0
  real :: initpower_ee=0.0, initpower2_ee=0.0
  real :: cutoff_ee=0.0, ncutoff_ee=0.0, kpeak_ee=0.0
  real :: relhel_ee=0.0, kgaussian_ee=0.0
  real :: ampla0=0.0, initpower_a0=0.0, initpower2_a0=0.0
  real :: cutoff_a0=0.0, ncutoff_a0=0.0, kpeak_a0=0.0
  real :: relhel_a0=0.0, kgaussian_a0=0.0, eta_ee=0.0
  real :: sigE_prefactor=1., sigB_prefactor=1.
  real :: weight_longitudinalE=2.0, mass_chi=0.
  real :: coupl_gy=.345 ! electroweak SU(2) x U(1) coupling of Higgs to U(1)
  logical :: luse_scale_factor_in_sigma=.false., lapply_Gamma_corr=.true.
  logical, pointer :: lohm_evolve, lphi_doublet, lphi_hypercharge
  logical, pointer :: lwaterfall
  real, pointer :: eta, Hscript, echarge, sigEm_all, sigBm_all
  integer :: iGamma=0, ia0=0, idiva_name=0, ieedot=0, iedotx=0, iedoty=0, iedotz=0
  integer :: idivE=0, isigE=0, isigB=0
  logical :: llongitudinalE=.true., llorenz_gauge_disp=.false., lskip_projection_ee=.false.
  logical :: lscale_tobox=.true., lskip_projection_a0=.false.
  logical :: lpower_profile_file=.false.
  logical :: lvectorpotential=.false., lphi_hom=.false., lphi_linear_regime=.false.
  logical :: lpsi_hom=.false.
  logical :: lno_noise_ee=.false., lnoncollinear_EB=.false., lnoncollinear_EB_aver=.false.
  logical :: lcollinear_EB=.false., lcollinear_EB_aver=.false.
  logical :: leedot_as_aux=.false., lcurlyA=.true., lsolve_chargedensity=.false.
  logical :: ldivE_as_aux=.false., lsigE_as_aux=.false., lsigB_as_aux=.false.
  logical :: lrandom_ampl_ee=.false., lfixed_phase_ee=.false., lallow_bprime_zero=.true.
  logical :: lswitch_off_divJ=.false., lswitch_off_Gamma=.false., lmass_suppression=.false.
  logical :: loverride_c_light=.false.
  character(len=labellen) :: inita0='zero'
  character (len=labellen), dimension(ninit) :: initee='nothing'
  character (len=labellen) :: power_filename='power_profile.dat'
!
  namelist /special_init_pars/ &
    initee, inita0, alpf, &
    ampl_ex, ampl_ey, ampl_ez, ampl_a0, &
    kx_ex, kx_ey, kx_ez, &
    ky_ex, ky_ey, ky_ez, &
    kz_ex, kz_ey, kz_ez, &
    kx_a0, ky_a0, kz_a0, &
    phase_ex, phase_ey, phase_ez, phase_a0, &
    llongitudinalE, llorenz_gauge_disp, lphi_hom, lphi_linear_regime, &
    amplee, initpower_ee, initpower2_ee, lscale_tobox, &
    cutoff_ee, ncutoff_ee, kpeak_ee, relhel_ee, kgaussian_ee, &
    ampla0, initpower_a0, initpower2_a0, lno_noise_ee, &
    cutoff_a0, ncutoff_a0, kpeak_a0, relhel_a0, kgaussian_a0, &
    leedot_as_aux, ldivE_as_aux, lsigE_as_aux, lsigB_as_aux, &
    lsolve_chargedensity, weight_longitudinalE, lswitch_off_Gamma, &
    lrandom_ampl_ee, lfixed_phase_ee, lskip_projection_ee, &
    luse_scale_factor_in_sigma, lpower_profile_file, power_filename, &
    coupl_gy, alpfpsi, lpsi_hom, loverride_c_light
!
  ! run parameters
  real :: beta_inflation=0., rescale_ee=1.
  logical :: reinitialize_ee=.false.
  namelist /special_run_pars/ &
    alpf, llongitudinalE, llorenz_gauge_disp, lphi_hom, lphi_linear_regime, &
    leedot_as_aux, ldivE_as_aux, lsigE_as_aux, lsigB_as_aux, &
    eta_ee, lcurlyA, beta_inflation, &
    weight_longitudinalE, lswitch_off_divJ, lswitch_off_Gamma, &
    lnoncollinear_EB, lnoncollinear_EB_aver, luse_scale_factor_in_sigma, &
    lcollinear_EB, lcollinear_EB_aver, sigE_prefactor, sigB_prefactor, &
    reinitialize_ee, initee, rescale_ee, lmass_suppression, mass_chi, &
    lallow_bprime_zero, lapply_Gamma_corr, coupl_gy, lpsi_hom, alpfpsi, &
    loverride_c_light
!
! Declare any index variables necessary for main or
!
  real :: c_light2
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_EEEM=0       ! DIAG_DOC: $\left<\Ev^2+\Bv^2\right>/2$
  integer :: idiag_erms=0       ! DIAG_DOC: $\left<\Ev^2\right>^{1/2}$
  integer :: idiag_eprimerms=0  ! DIAG_DOC: $\left<(E')^2\right>^{1/2}$
  integer :: idiag_bprimerms=0  ! DIAG_DOC: $\left<(B')^2\right>^{1/2}$
  integer :: idiag_jprimerms=0  ! DIAG_DOC: $\left<(J')^2\right>^{1/2}$
  integer :: idiag_gam_EBrms=0  ! DIAG_DOC: $\left<(\gamma')^2\right>^{1/2}$
  integer :: idiag_boostprms=0  ! DIAG_DOC: $\left<\mbox{boost}^2\right>^{1/2}$
  integer :: idiag_edotrms=0    ! DIAG_DOC: $\left<\dot{\Ev}^2\right>^{1/2}$
  integer :: idiag_emax=0       ! DIAG_DOC: $\max(|\Ev|)$
  integer :: idiag_a0rms=0      ! DIAG_DOC: $\left<A_0^2\right>^{1/2}$
  integer :: idiag_grms=0       ! DIAG_DOC: $\left<C-\nabla\cdot\Av\right>^{1/2}$
  integer :: idiag_da0rms=0     ! DIAG_DOC: $\left<C-\nabla\cdot\Av\right>^{1/2}$
  integer :: idiag_BcurlEm=0    ! DIAG_DOC: $\left<\Bv\cdot\nabla\times\Ev\right>$
  integer :: idiag_divJrms=0    ! DIAG_DOC: $\left<\nab\Jv^2\right>^{1/2}$
  integer :: idiag_divErms=0    ! DIAG_DOC: $\left<\nab\Ev^2\right>^{1/2}$
  integer :: idiag_rhoerms=0    ! DIAG_DOC: $\left<\rho_e^2\right>^{1/2}$
  integer :: idiag_divJm=0      ! DIAG_DOC: $\left<\nab\Jv\right>$
  integer :: idiag_divEm=0      ! DIAG_DOC: $\left<\nab\Ev\right>$
  integer :: idiag_rhoem=0      ! DIAG_DOC: $\left<\rho_e\right>$
  integer :: idiag_count_eb0=0  ! DIAG_DOC: $f_\mathrm{EB0}$
  integer :: idiag_mfpf=0       ! DIAG_DOC: $-f'/f$
  integer :: idiag_fppf=0       ! DIAG_DOC: $f''/f$
  integer :: idiag_afact=0      ! DIAG_DOC: $a$ (scale factor)
  integer :: idiag_constrainteqn=0  ! DIAG_DOC: $<deldotE+>$
  integer :: idiag_exm=0        ! DIAG_DOC: $\left<E_x\right>$
  integer :: idiag_eym=0        ! DIAG_DOC: $\left<E_y\right>$
  integer :: idiag_ezm=0        ! DIAG_DOC: $\left<E_z\right>$
  integer :: idiag_sigEm=0      ! DIAG_DOC: $\left<\sigma_\mathrm{E}\right>$
  integer :: idiag_sigBm=0      ! DIAG_DOC: $\left<\sigma_\mathrm{B}\right>$
  integer :: idiag_sigErms=0    ! DIAG_DOC: $\left<\sigma_\mathrm{E}^2\right>^{1/2}$
  integer :: idiag_sigBrms=0    ! DIAG_DOC: $\left<\sigma_\mathrm{B}^2\right>^{1/2}$
  integer :: idiag_sigEE2m=0    ! DIAG_DOC: $\left<\sigma_\mathrm{E}\Ev^2\right>$
  integer :: idiag_sigBBEm=0    ! DIAG_DOC: $\left<\sigma_\mathrm{E}\Bv\cdot\Ev\right>$
  integer :: idiag_adphiBm=0    ! DIAG_DOC: $\left<(\alpha/f)<\phi'\Bv\cdot\Ev\right>$
  integer :: idiag_Johmrms=0    ! DIAG_DOC: $\left<\Jv^2\right>^{1/2}$
  integer :: idiag_echarge=0    ! DIAG_DOC: $\left<e_\mathrm{eff}\right>$
  integer :: idiag_ebm=0        ! DIAG_DOC: $\left<\Ev\cdot\Bv\right>$
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_exmz=0       ! XYAVG_DOC: $\left<{\cal E}_x\right>_{xy}$
  integer :: idiag_eymz=0       ! XYAVG_DOC: $\left<{\cal E}_y\right>_{xy}$
  integer :: idiag_ezmz=0       ! XYAVG_DOC: $\left<{\cal E}_z\right>_{xy}$
!
! yz averaged diagnostics given in yzaver.in
!
  integer :: idiag_e2mx = 0     ! YZAVG_DOC: $\langle E^2\rangle_{yz}$
!
  contains
!
!***********************************************************************
    subroutine register_special
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  18-mar-21/axel: coded Faraday displacement current
!
      use FArrayManager
      use Sub, only: register_report_aux
      use SharedVariables, only: put_shared_variable

      ldisp_current =.true.
!
!  It would have been more consistent to call the indices to the
!  three components iex, iey, and iez
!
      call farray_register_pde('ee',iee,vector=3)
      iex=iee; iey=iee+1; iez=iee+2
!
      if (leedot_as_aux) &
        call register_report_aux('eedot', ieedot, iedotx, iedoty, iedotz)
!
      if (ldivE_as_aux) &
        call register_report_aux('divE', idivE)
!
      if (lsigE_as_aux) &
        call register_report_aux('sigE', isigE)
!
      if (lsigB_as_aux) &
        call register_report_aux('sigB', isigB)
!
      if (lsolve_chargedensity) &
        call farray_register_pde('rhoe',irhoe)
!
      if (llorenz_gauge_disp) then
        call farray_register_pde('a0',ia0)
        call farray_register_pde('diva_name',idiva_name)
      endif
!
!  For llongitudinalE=T, we replace graddiv in curlb by the gradient of Gamma.
!  According to later work, this does not seem advantageous, however.
!
      if (llongitudinalE) &
        call farray_register_pde('Gamma',iGamma)
!
!  The following variables are also used in special/backreact_infl.f90
!
      call put_shared_variable('alpf',alpf,caller='register_disp_current')
      call put_shared_variable('alpfpsi',alpfpsi)
      call put_shared_variable('lphi_hom',lphi_hom)
      call put_shared_variable('lpsi_hom',lpsi_hom)
      call put_shared_variable('lphi_linear_regime',lphi_linear_regime)
      call put_shared_variable('sigE_prefactor',sigE_prefactor)
      call put_shared_variable('sigB_prefactor',sigB_prefactor)
      call put_shared_variable('lcollinear_EB',lcollinear_EB)
      call put_shared_variable('lcollinear_EB_aver',lcollinear_EB_aver)
      call put_shared_variable('lnoncollinear_EB',lnoncollinear_EB)
      call put_shared_variable('lnoncollinear_EB_aver',lnoncollinear_EB_aver)
      call put_shared_variable('lmass_suppression',lmass_suppression)
      call put_shared_variable('lallow_bprime_zero',lallow_bprime_zero)
      call put_shared_variable('mass_chi',mass_chi)
      call put_shared_variable('llongitudinalE',llongitudinalE)
      call put_shared_variable('coupl_gy',coupl_gy)
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  20-mar-21/axel: coded
!
      use FArrayManager
      use SharedVariables, only: get_shared_variable
      use Initcond, only: gaunoise
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
!  Initialize module variables which are parameter dependent
!  If one really wants to work with c_light /= 1,
!  then one needs to override this (loverride_c_light=T).
!
      if (c_light/=1. .and. .not. loverride_c_light) call fatal_error('disp_current', &
          "use unit_system='set' or put loverride_c_light=T")
      c_light2=c_light**2
!
      if (lmagnetic .and. .not.lswitch_off_divJ) &
        call get_shared_variable('eta',eta, caller='initialize_magnetic')
!
!  The following are only obtained when luse_scale_factor_in_sigma=T
!  (luse_scale_factor_in_sigma=F by default, because they are defined
!  in special/backreact_infl.f90, which may not be always be used).
!
      if (luse_scale_factor_in_sigma) then
        call get_shared_variable('Hscript', Hscript)
        call get_shared_variable('echarge', echarge)
        call get_shared_variable('sigEm_all', sigEm_all)
        call get_shared_variable('sigBm_all', sigBm_all)
        call get_shared_variable('lohm_evolve', lohm_evolve)
      else
        if (.not.associated(Hscript)) allocate(Hscript,echarge,sigEm_all,sigBm_all)
        Hscript=0.
        echarge=0.
        sigEm_all=0.
        sigBm_all=0.
        allocate(lohm_evolve)
        lohm_evolve=.false.
      endif
!
!  Reinitialize magnetic field using a small selection of perturbations
!  that were mostly also available as initial conditions.
!
      if (reinitialize_ee) then
        do j=1,ninit
          select case (initee(j))
          case ('rescale')
            f(:,:,:,iex:iez)=rescale_ee*f(:,:,:,iex:iez)
            if (llongitudinalE) f(:,:,:,iGamma)=rescale_ee*f(:,:,:,iGamma)
          case ('gaussian-noise'); call gaunoise(amplee(j),f,iex,iez)
          case default
          endselect
        enddo
      endif

      if (lklein_gordon) then
        call get_shared_variable('lphi_doublet',lphi_doublet, caller='initialize_disp_current')
        call get_shared_variable('lphi_hypercharge',lphi_hypercharge, &
          caller='initialize_disp_current')
        call get_shared_variable('lwaterfall',lwaterfall, caller='initialize_disp_current')
      else
        if (.not.associated(lphi_doublet)) allocate(lphi_doublet,lphi_hypercharge,lwaterfall)
        lphi_doublet=.false.
        lphi_hypercharge=.false.
        lwaterfall=.false.
      endif
!
      if (lphi_hom) weight_longitudinalE=0.
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Initcond
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: divA, divE
      integer :: j
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      do j=1,ninit
        select case (initee(j))
          case ('nothing'); if (lroot) print*,'initee: nothing'
          case ('zero'); f(:,:,:,iex:iez)=0.
          case ('coswave-phase')
            call coswave_phase(f,iex,ampl_ex,kx_ex,ky_ex,kz_ex,phase_ex)
            call coswave_phase(f,iey,ampl_ey,kx_ey,ky_ey,kz_ey,phase_ey)
            call coswave_phase(f,iez,ampl_ez,kx_ez,ky_ez,kz_ez,phase_ez)
          case ('sinwave-phase')
            call sinwave_phase(f,iex,ampl_ex,kx_ex,ky_ex,kz_ex,phase_ex)
            call sinwave_phase(f,iey,ampl_ey,kx_ey,ky_ey,kz_ey,phase_ey)
            call sinwave_phase(f,iez,ampl_ez,kx_ez,ky_ez,kz_ez,phase_ez)
          case ('power_randomphase_hel')
            call power_randomphase_hel(amplee(j),initpower_ee,initpower2_ee, &
              cutoff_ee,ncutoff_ee,kpeak_ee,f,iex,iez,relhel_ee,kgaussian_ee, &
              lskip_projection_ee, lvectorpotential, lscale_tobox=lscale_tobox, &
              lpower_profile_file=lpower_profile_file, power_filename=power_filename, &
              lno_noise=lno_noise_ee, lrandom_ampl=lrandom_ampl_ee, &
              lfixed_phase=lfixed_phase_ee)
          case default
            call fatal_error('init_ee','no such init_ee: "'//trim(initee(j))//'"')
        endselect
      enddo
!
!  Initial condition for Gamma = divA
!
      if (llongitudinalE) then
        do n=n1,n2; do m=m1,m2
          call div(f,iaa,diva)
          f(l1:l2,m,n,iGamma)=diva
        enddo; enddo
      endif
!
!  Initial condition for rhoe = divE
!
      if (lsolve_chargedensity) then
        do n=n1,n2; do m=m1,m2
          call div(f,iee,divE)
          f(l1:l2,m,n,irhoe)=divE
        enddo; enddo
      endif
!
!  Initialize diva_name if llorenz_gauge_disp=T
!
      if (llorenz_gauge_disp) then
        do n=n1,n2; do m=m1,m2
          call div(f,iaa,diva)
          f(l1:l2,m,n,idiva_name)=diva
        enddo; enddo
!
!  initial conditions for A0 (provided llorenz_gauge_disp=T)
!
        select case (inita0)
          case ('coswave-phase')
            call coswave_phase(f,ia0,ampl_a0,kx_a0,ky_a0,kz_a0,phase_a0)
          case ('zero'); f(:,:,:,ia0)=0.
          case ('power_randomphase')
            call power_randomphase_hel(ampla0,initpower_a0,initpower2_a0, &
              cutoff_a0,ncutoff_a0,kpeak_a0,f,ia0,ia0, &
              relhel_a0,kgaussian_a0, lskip_projection_a0, lvectorpotential, &
              lscale_tobox, lpower_profile_file=.false.)
          case default
            call fatal_error("init_special","no such inita0: "//trim(inita0))
        endselect
      endif
!
      if (leedot_as_aux) f(:,:,:,iedotx:iedotz)=0.
      if (ldivE_as_aux) f(:,:,:,idivE)=0.
      if (lsigE_as_aux) f(:,:,:,isigE)=0.
      if (lsigB_as_aux) f(:,:,:,isigB)=0.
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  25-feb-07/axel: adapted
!
!
!  compulsory pencils
!
      lpenc_requested(i_aa)=.true.
      lpenc_requested(i_el)=.true.
      ! alberto: should be this pencil only requested if llorenz_gauge_disp=T?
      !lpenc_requested(i_ga0)=.true.
      lpenc_requested(i_curlb)=.true.
      lpenc_requested(i_gGamma)=.true.
      lpenc_requested(i_jj_ohm)=.true.
!
!  Pencils for axion-like coupling alpf * phi F Fdual
!
      if (alpf/=0.) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_phi)=.true.
        lpenc_requested(i_dphi)=.true.
        lpenc_requested(i_gphi)=.true.
      endif
      if (alpfpsi/=0. .and. lwaterfall) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_dpsi)=.true.
        lpenc_requested(i_gpsi)=.true.
      endif
!
      if (llorenz_gauge_disp) lpenc_requested(i_ga0)=.true.
!
! Pencils for lnoncollinear_EB and lcollinear_EB cases.
!
      if (lnoncollinear_EB .or. lnoncollinear_EB_aver &
        .or. lcollinear_EB .or. lcollinear_EB_aver) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_b2)=.true.
        lpenc_requested(i_e2)=.true.
      endif
!
   !  if (lnoncollinear_EB) then
   !    lpenc_requested(i_eb)=.true.
   !  endif
      lpenc_requested(i_eb)=.true.
!
      !if (llorenz_gauge_disp) then
      !  lpenc_requested(i_diva)=.true.
      !endif
!
!  Terms for Gamma evolution.
!
      if (llongitudinalE) then
        lpenc_requested(i_divE)=.true.
        ! lpenc_requested(i_gGamma)=.true.
        ! alberto: gGamma is always requested, as curlb depends on it
      endif
!
      if (idiag_divEm/=0. .or. idiag_divErms/=0.) then
        lpenc_requested(i_divE)=.true.
      endif
!
!  charge density
!
      if (lsolve_chargedensity) then
        lpenc_requested(i_divJ)=.true.
        lpenc_requested(i_uij)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_bb)=.true.
      endif
!
!  diffusion term.
!
      if (eta_ee/=0.) lpenc_requested(i_del2ee)=.true.
!
!  Higgs pencils
!
      if (lklein_gordon .and. lphi_doublet .and. lphi_hypercharge) then
        lpenc_requested(i_phi)=.true.
        lpenc_requested(i_phi_doublet)=.true.
        lpenc_requested(i_cov_der)=.true.
        lpenc_requested(i_jj_higgsY)=.true.
      endif
!
!  Diagnostics pencils:
!
      ! if (eta_ee/=0.) lpenc_requested(i_del2ee)=.true.
      ! if (eta_ee/=0.) lpenc_requested(i_del2ee)=.true.

      if (idiag_BcurlEm/=0) then
        lpenc_diagnos(i_curlE)=.true.
        lpenc_diagnos(i_BcurlE)=.true.
      endif
!
      if (idiag_ebm/=0) lpenc_diagnos(i_eb)=.true.
      if (idiag_a0rms/=0) lpenc_diagnos(i_a0)=.true.
      if (idiag_grms/=0) lpenc_diagnos(i_diva)=.true.
      if (idiag_edotrms/=0) lpenc_diagnos(i_edot2)=.true.
      if (idiag_EEEM/=0 .or. idiag_erms/=0 .or. idiag_emax/=0 & 
        .or. idiag_e2mx/=0) lpenc_diagnos(i_e2)=.true.
      ! if (idiag_exmz/=0 .or. idiag_eymz/=0 .or. idiag_ezmz/=0 ) lpenc_diagnos(i_el)=.true.
      ! if (idiag_exm/=0 .or. idiag_eym/=0 .or. idiag_ezm/=0 ) lpenc_diagnos(i_el)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      use Sub, only: grad, div, curl, del2v, dot2_mn, dot, levi_civita, del2v_etc
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx) :: tmp, mass_suppression_fact
      integer :: i,j,k
!
      intent(inout) :: f
      intent(inout) :: p
!
!  Pencil for charge density.
!
      if (lsolve_chargedensity) p%rhoe=f(l1:l2,m,n,irhoe)
!
!  Terms for Gamma evolution.
!
      if (lpenc_requested(i_divE)) call div(f,iee,p%divE)
      if (lpenc_requested(i_gGamma)) then
        if (llongitudinalE) then
          call grad(f,iGamma,p%gGamma)
        ! alberto: when llongitudinalE=F, we should compute
        ! grad div a from f-array
        else
          !call del2v_etc(f,iaa,GRADDIV=p%gGamma)
          !16-sep-25/axel: but in this case, pGamma should not be needed.
          call fatal_error("calc_pencils_special","Gamma is not defined")
        endif
      endif
!
!  Replace p%curlb by the combination -p%del2a+p%gGamma.
!
      ! alberto: as curlb is always requested, this allows to
      ! compute curlb also if llongitudinalE=F
      !16-sep-25/axel: but curlb was already computed in calc_pencils_magnetic_pencpar (in magnetic).
      if (lpenc_requested(i_curlb)) then
        if (lapply_Gamma_corr) then
          if (llongitudinalE) then
            p%curlb=-p%del2a+p%gGamma
          else
            call fatal_error("calc_pencils_special","Gamma is not defined")
          endif
        endif
      endif
      ! if (llongitudinalE) then
      !   ! call div(f,iee,p%divE)
      !   ! call grad(f,iGamma,p%gGamma)
      !   if (lapply_Gamma_corr) p%curlb=-p%del2a+p%gGamma
      !   ! if (lsolve_chargedensity) p%rhoe=f(l1:l2,m,n,irhoe)
      ! endif
!
! el and e2 (note that this is called after magnetic, where sigma is computed)
!
      p%el=f(l1:l2,m,n,iex:iez)
      call dot2_mn(p%el,p%e2)
!
!  eb pencil
!
      if (lpenc_requested(i_eb)) then
        call dot(p%el,p%bb,p%eb)
      endif
!
!  Compute fully non-collinear expression for the current density.
!  This is for the *spatially dependent* sigE and sigB.
!  The averaged ones are computed in backreact_infl.f90.
!  Any change to the code must be the same both here and there.
!  Location 1 for conductivity.
!
      if (lnoncollinear_EB .or. lnoncollinear_EB_aver &
        .or. lcollinear_EB .or. lcollinear_EB_aver) then
        if (lnoncollinear_EB) then
          p%boost=sqrt((p%e2-p%b2)**2+4.*p%eb**2)
          p%gam_EB=sqrt21*sqrt(1.+(p%e2+p%b2)/p%boost)
          p%eprime=sqrt21*sqrt(p%e2-p%b2+p%boost)
          p%bprime=sqrt21*sqrt(p%b2-p%e2+p%boost)*sign(1.,p%eb)
          if (lallow_bprime_zero) then
            where (p%eprime/=0.)
              where (p%bprime/=0.)
                p%jprime=Chypercharge*echarge**3/(6.*pi**2*Hscript)*p%eprime*abs(p%bprime)/tanh(pi*abs(p%bprime)/p%eprime)
              elsewhere
                p%jprime=Chypercharge*echarge**3/(6.*pi**3*Hscript)*p%eprime**2
              endwhere
              p%sigE=sigE_prefactor*abs(p%jprime)*p%eprime/(p%gam_EB*p%boost)
              p%sigB=sigB_prefactor*abs(p%jprime)*p%eb/(p%eprime*p%gam_EB*p%boost)
              p%count_eb0=0.
            elsewhere
              p%jprime=0.
              p%sigE=0.
              p%sigB=0.
              p%count_eb0=1.
            endwhere
          else
            where (p%eprime/=0. .and. p%bprime/=0.)
              p%jprime=Chypercharge*echarge**3/(6.*pi**2*Hscript)*p%eprime*abs(p%bprime)/tanh(pi*abs(p%bprime)/p%eprime)
              p%sigE=sigE_prefactor*abs(p%jprime)*p%eprime/(p%gam_EB*p%boost)
              p%sigB=sigB_prefactor*abs(p%jprime)*p%eb/(p%eprime*p%gam_EB*p%boost)
              p%count_eb0=0.
            elsewhere
              p%jprime=0.
              p%sigE=0.
              p%sigB=0.
              p%count_eb0=1.
            endwhere
          endif

        elseif (lcollinear_EB) then
          p%eprime=sqrt(p%e2)
          p%bprime=sqrt(p%b2)
          where (p%eprime/=0. .and. p%bprime/=0.)
            p%sigE=sigE_prefactor*Chypercharge*echarge**3/(6.*pi**2*Hscript)*p%bprime/tanh(pi*abs(p%bprime)/p%eprime)
            p%count_eb0=0.
          elsewhere
            p%sigE=0.
            p%count_eb0=1.
          endwhere
          if (lmass_suppression) then
            mass_suppression_fact=exp(-pi*mass_chi**2/(Chypercharge**onethird*echarge*sqrt(p%e2)))
            p%sigE=p%sigE*mass_suppression_fact
          endif
          p%sigB=0.
        elseif (lnoncollinear_EB_aver .or. lcollinear_EB_aver) then
!
!  This is based on <E> and <B>. Later, when sigEm and sigBm diagonstics are being
!  computed, those are then based on the same p%sigE and p%sigB values, and therefore
!  also the same as sigEma and sigBma.
!
          p%sigE=sigEm_all
          p%sigB=sigBm_all
        endif
!
!  Now compute current, using any of the 4 expressions above.
!
        if (lohm_evolve) then
          p%jj_ohm=f(l1:l2,m,n,ijx:ijz)
        else
          do j=1,3
            p%jj_ohm(:,j)=p%sigE*p%el(:,j)+p%sigB*p%bb(:,j)
          enddo
          if (ijx/=0) f(l1:l2,m,n,ijx:ijz) = p%jj_ohm
        endif
!
      endif
!
! edot2
!
      if (leedot_as_aux) then
        call dot2_mn(f(l1:l2,m,n,iedotx:iedotz),p%edot2)
      else
        p%edot2=0.
      endif
!
!  del2ee  !(AB: not needed)
!
      if (eta_ee/=0.) call del2v(f,iex,p%del2ee)
!
!  curle
!
      if (idiag_BcurlEm/=0) then
        call curl(f,iex,p%curle)
        call dot(p%bb,p%curle,p%BcurlE)
      endif
! !
! !  del2ee
! !
!       if (eta_ee/=0.) call del2v(f,iex,p%del2ee)
!
! a0 & ga0
!
      if (ia0>0) then
        p%a0=f(l1:l2,m,n,ia0)
        call grad(f,ia0,p%ga0)
      endif
!
!  divJ (using Ohm's law)
!  divJ=sigma*[divE+eps_ijk*(u_j,i * b_k + u_j * b_k,i)]
!  The use if eta may be suspect and should be checked.
!
      if (lpenc_requested(i_divJ)) then
        tmp=0.
        do i=1,3
        do j=1,3
        do k=1,3
          tmp=tmp+levi_civita(i,j,k)*(p%uij(:,j,i)*p%bb(:,k)+p%uu(:,j)*p%bij(:,k,i))
        enddo
        enddo
        enddo
        if (lswitch_off_divJ) then
          p%divJ=0.
        else
          if (eta==0.) then
!
!  The following expression ignores gradients of p%sigE
!
            !call fatal_error('disp_current/calc_pencils_special', "eta=0 not ok here")
            p%divJ=(p%divE+tmp)*p%sigE
          else
            p%divJ=(p%divE+tmp)/(mu0*eta)
          endif
        endif
      endif

!  pencils for klein_gordon module
      if (lpenc_requested(i_GammaY)) then
        if (llongitudinalE) then
          p%GammaY=f(l1:l2,m,n,iGamma)
        else
          call div(f,iaa,p%GammaY)
        endif
      endif
      if (alpf/=0.and..not.lklein_gordon) p%dphi=p%infl_dphi
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine calc_constrainteqn(p,tmp,constrainteqn)
      type(pencil_case) :: p
      real, dimension(nx),intent(IN) :: tmp
      real, dimension(nx), intent(OUT) :: constrainteqn
      real, dimension(nx) :: constrainteqn1
      !constrainteqn1=sqrt(p%divE**2+tmp**2)
!
!  in the following, should use "where"
!
      constrainteqn1=sqrt(p%divE**2+tmp**2)
      if (any(constrainteqn1 == 0.)) then
        constrainteqn=0.
      else
        constrainteqn=(p%divE-tmp)/constrainteqn1
      endif
     endsubroutine calc_constrainteqn
!***********************************************************************
      real function get_mfpf()
              get_mfpf = beta_inflation*Hp_target
      end function get_mfpf
!***********************************************************************
      real function get_fppf()
              get_fppf=beta_inflation*((beta_inflation+1.)*Hp_target**2-appa_target)
      end function get_fppf
!***********************************************************************
      subroutine calc_axion_term(p,dst,gphi,alpff,lphihom)
!
!  Compute -(alpha/f)*B.gradphi axion term (when alpha/f/=0).
!
      use Sub
!
      type(pencil_case) :: p
      real, intent(in) :: alpff
      logical, intent(in) :: lphihom
      real, dimension(nx), intent(out) :: dst
      real, dimension(nx,3), intent(in) :: gphi
!
!   14-sep-25/alberto: added axion coupling for a second scalar field psi
!
      if (lphihom) then
        dst=0.
      else
        if (alpff/=0.) then
        ! call dot(p%bb,p%gphi,dst2)
          call dot(p%bb,gphi,dst)
          dst=-alpff*dst
        else
         dst=0.
        endif
      endif

      ! if (.not. lpsi_hom .and. alpfpsi/=0. .and. lwaterfall) then
      !   call dot(p%bb,p%gpsi,dst2)
      !   dst=dst-alpfpsi*dst2
      ! endif

      endsubroutine calc_axion_term
!***********************************************************************
      subroutine calc_helical_term(p,gtmp,dphi,gphi,lphihom)
          use Sub
          type(pencil_case), intent(IN) :: p
          real, dimension(nx,3), intent(out) :: gtmp
          logical, intent(in) :: lphihom
          real, dimension(nx), intent(in) :: dphi
          real, dimension(nx,3), intent(in) :: gphi
!
          if (lphihom) then
            call multsv(dphi,p%bb,gtmp)
          else
            call cross(gphi,p%el,gtmp)
            call multsv_add(gtmp,dphi,p%bb,gtmp)
          endif

          ! if (alpfpsi/=0. .and. lwaterfall) then
          !   if (lpsi_hom) then
          !     call multsv(p%dpsi,p%bb,gtmp2)
          !   else
          !     call cross(p%gpsi,p%el,gtmp2)
          !     call multsv_add(gtmp2,p%dpsi,p%bb,gtmp2)
          !   endif
          !   gtmp=gtmp+gtmp2*alpfpsi
          ! endif

      endsubroutine
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   18-mar-21/axel: coded Faraday displacement current
!   07-sep-25/alberto: coded charges from Higgs doublet
!
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: gtmp, dJdt, del2JJ
      real, dimension (nx) :: tmp, tmp2, del2a0, constrainteqn, constrainteqn1
      real :: inflation_factor=0., mfpf=0., fppf=0.
      integer :: j
!
      intent(inout) :: p
      intent(inout) :: f, df
      integer :: i
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
      if (headtt) call identify_bcs('ee',iee)
!
!  Calculate rhs of Gamma equation and update curl
!  Initialize tmp with axion term.
!
      ! call calc_axion_term(p,tmp)
      call calc_axion_term(p,tmp,p%gphi,alpf,lphi_hom)
      if (lwaterfall) then
        call calc_axion_term(p,tmp2,p%gpsi,alpfpsi,lpsi_hom)
        tmp=tmp+tmp2
      endif
!
!  Solve for Gamma (unless lswitch_off_Gamma) and possibly for charge density.
!  Add to the existing tmp and update df(l1:l2,m,n,irhoe).
!
      if (llongitudinalE) then
        if (lsolve_chargedensity) tmp=tmp+f(l1:l2,m,n,irhoe)
        ! add charge from Higgs field to Gauss constraint
        if (lphi_doublet .and. lphi_hypercharge .and. coupl_gy /= 0) then
          p%rhoe_higgsY=-coupl_gy*(p%phi*p%cov_der(:,1,2) - &
                  p%phi_doublet(:,1)*p%cov_der(:,1,1) + &
                  p%phi_doublet(:,2)*p%cov_der(:,1,4) - &
                  p%phi_doublet(:,3)*p%cov_der(:,1,3))
          tmp = tmp + p%rhoe_higgsY
        endif
        if (.not.lswitch_off_Gamma) df(l1:l2,m,n,iGamma)=df(l1:l2,m,n,iGamma) &
          -(1.-weight_longitudinalE)*p%divE-weight_longitudinalE*tmp
      endif
!
!  solve: dE/dt = curlB - ...
!  Calculate curlB as -del2a, because curlB leads to instability.
!  Solve dA/dt = -E.
!
      if (lmagnetic) then
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-p%el
!
!  Maxwell equation otherwise the same in both gauges.
!
        df(l1:l2,m,n,iex:iez)=df(l1:l2,m,n,iex:iez)+c_light2*(p%curlb-mu0*p%jj_ohm)
!
!  Solve for charge density
!
        if (lsolve_chargedensity) df(l1:l2,m,n,irhoe)=df(l1:l2,m,n,irhoe)-p%divJ
!
!  Magneto-genesis from reheating. In the papers by Subramanian (2010) and Sharma+17,
!  as well as BS21, the calliographic variable curly-A=f*A was introduced to get
!  rid of the first derivative of A. But the disadvantage is that the generation
!  term, (f"/f)*<A.E> is then gauge-dependent. Because of this and other reasons,
!  it is better to work with the original 2(f'/f)*A' = -2(f'/f)*E term, which is
!  gauge-independent.
!
        if (beta_inflation/=0.) then
          if (ip<14.and.lroot) print*,'scl_factor_target, Hp_target, appa_target, wweos_target=', &
                                       scl_factor_target, Hp_target, appa_target, wweos_target
          mfpf=beta_inflation*Hp_target
          fppf=beta_inflation*((beta_inflation+1.)*Hp_target**2-appa_target)
          if (lcurlyA) then
            inflation_factor=fppf
            df(l1:l2,m,n,iex:iez)=df(l1:l2,m,n,iex:iez)-c_light2*inflation_factor*p%aa
          else
            inflation_factor=-2.*mfpf
            df(l1:l2,m,n,iex:iez)=df(l1:l2,m,n,iex:iez)-c_light2*inflation_factor*p%el
          endif
          if (ip<15.and.lroot.and.lfirst) print*,'t, inflation_factor=',t, inflation_factor
        endif
!
!  if particles, would add J=sum(qi*Vi*ni)
!
!  A0 equation: 3 terms
!  dA0/dt = divA
!  dAA/dt = ... + gradA0
!
        ! alberto: llorenz_gauge_disp can also be considered when
        ! alpf=0, moved addition of del2a0 from conditions below
        if (llorenz_gauge_disp) then
          call del2(f,ia0,del2a0)
          df(l1:l2,m,n,idiva_name)=df(l1:l2,m,n,idiva_name)+del2a0
        endif

!  helical term:
!  dEE/dt = ... -alp/f (dphi*BB + gradphi x E)
!  Use the combined routine multsv_add if both terms are included.
!
        if (alpf/=0.) then
          call calc_helical_term(p,gtmp,p%dphi,p%gphi,lphi_hom)
!          print*,"p%infl_phi",p%infl_phi
!          print*,"p%infl_dphi",p%infl_dphi
          df(l1:l2,m,n,iex:iez)=df(l1:l2,m,n,iex:iez)-alpf*gtmp
          if (llorenz_gauge_disp) then
            ! if (lphi_hom) then
            !   df(l1:l2,m,n,idiva_name)=df(l1:l2,m,n,idiva_name)+del2a0
            ! else
            if (.not. lphi_hom) then
              call dot_mn(p%gphi,p%bb,tmp)
              df(l1:l2,m,n,idiva_name)=df(l1:l2,m,n,idiva_name)+alpf*tmp
            endif
          endif
        endif
        if (lwaterfall .and. alpfpsi/=0.) then
          call calc_helical_term(p,gtmp,p%dpsi,p%gpsi,lpsi_hom)
          df(l1:l2,m,n,iex:iez)=df(l1:l2,m,n,iex:iez)-alpfpsi*gtmp
          if (llorenz_gauge_disp) then
            if (.not. lpsi_hom) then
              call dot_mn(p%gpsi,p%bb,tmp)
              df(l1:l2,m,n,idiva_name)=df(l1:l2,m,n,idiva_name)+alpfpsi*tmp
            endif
          endif
        endif
!
            !endif
!
!  Evolution of the equation for the scalar potential, moved below
!
            !df(l1:l2,m,n,ia0)=df(l1:l2,m,n,ia0)+p%diva
            !df(l1:l2,m,n,ia0)=df(l1:l2,m,n,ia0)+f(l1:l2,m,n,idiva_name)
            !df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+p%ga0
        !   endif
        ! endif
        ! alberto: If llorenz_gauge_disp, add the two terms to the A0 and A equations
        ! here, so that they are always added, regardless of alpf.
        if (llorenz_gauge_disp) then
          df(l1:l2,m,n,ia0)=df(l1:l2,m,n,ia0)+f(l1:l2,m,n,idiva_name)
          df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+p%ga0
        endif
        if (eta_ee/=0.) df(l1:l2,m,n,iex:iez)=df(l1:l2,m,n,iex:iez)+c_light2*eta_ee*p%del2ee
!
!  If Higgs doublet, add current from Higgs U(1) hypercharge.
!
        if (lphi_doublet .and. lphi_hypercharge .and. coupl_gy /= 0) then
          ! compute jj_higgsY pencil (not done in calc_pencils_special as
          ! pencils from other special modules might not be available there yet)
          do i=1,3
            p%jj_higgsY(:,i)=coupl_gy*(p%phi*p%cov_der(:,i+1,2) - &
                      p%phi_doublet(:,1)*p%cov_der(:,i+1,1) + &
                      p%phi_doublet(:,2)*p%cov_der(:,i+1,4) - &
                      p%phi_doublet(:,3)*p%cov_der(:,i+1,3))
          enddo
          ! do i=1,3
          !   df(l1:l2,m,n,iex+i-1)=df(l1:l2,m,n,iex+i-1) - &
          !         coupl_gy*(p%phi*p%cov_der(:,i+1,2) - &
          !         p%phi_doublet(:,1)*p%cov_der(:,i+1,1) + &
          !         p%phi_doublet(:,2)*p%cov_der(:,i+1,4) - &
          !         p%phi_doublet(:,3)*p%cov_der(:,i+1,3))
          ! enddo
          df(l1:l2,m,n,iex:iez)=df(l1:l2,m,n,iex:iez) - p%jj_higgsY
        endif
      endif
!
!  Compute eedot_as_aux; currently ignore alpf/=0.
!
      if (leedot_as_aux) f(l1:l2,m,n,iedotx:iedotz)=c_light2*(p%curlb-mu0*p%jj_ohm)
!
!  Evolving current
!
      if (lnoncollinear_EB .or. lnoncollinear_EB_aver &
        .or. lcollinear_EB .or. lcollinear_EB_aver) then
        if (lohm_evolve) then
    !     if (tau_jj>0) then
    !       tau1_jj=1./tau_jj
    !       do j=1,3
!
!  Here we would need to add tau*sigmaB*B
!
    !         dJdt(:,j)=tau1_jj*(p%el(:,j)+p%uxb(:,j))*mu01/eta_total
    !       enddo
    !       if (ell_jj/=0.) then
    !         call del2v(f,ijx,del2jj)
    !         dJdt=dJdt+(ell_jj**2*tau1_jj)*del2jj
    !       endif
    !       df(l1:l2,m,n,ijx:ijz)=df(l1:l2,m,n,ijx:ijz)+dJdt
    !     else
!if (f(l1,m,n,ijx)==0.) then
!  print*,'AXEL start j/=0 '
!            do j=1,3
!  f(l1:l2,m,n,ijx-1+j)=p%sigE*p%el(:,j)
!            enddo
!endif
          do j=1,3
            dJdt(:,j)=3.*Hscript*p%sigE*p%el(:,j)
          enddo
          df(l1:l2,m,n,ijx:ijz)=df(l1:l2,m,n,ijx:ijz)+dJdt
        endif
      endif
!
!  If requested, put sigE and sigB into f array as auxiliaries.
!
      if (ldivE_as_aux) f(l1:l2,m,n,idivE)=p%divE
      if (lsigE_as_aux) f(l1:l2,m,n,isigE)=p%sigE
      if (lsigB_as_aux) f(l1:l2,m,n,isigB)=p%sigB
!
!  timestep constraint
!
      if (lfirst.and.ldt) advec_cs2=max(advec_cs2,c_light2*dxyz_2)
!
!  diagnostics
!
      if (ldiagnos) then
        call calc_diagnostics_special(f,p)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine calc_diagnostics_special(f,p)
      use Sub
      use Diagnostics
      real, dimension(mx,my,mz,mfarray) :: f
      type(pencil_case) :: p
      real, dimension(nx) :: tmp,constrainteqn
      real :: mfpf=0.,fppf=0.
      real, dimension(nx,3) :: gtmp

      call save_name(get_mfpf(),idiag_mfpf)
      call save_name(get_fppf(),idiag_fppf)
      call save_name(scl_factor_target,idiag_afact)
      if (idiag_EEEM/=0) call sum_mn_name(.5*(p%e2+p%b2),idiag_EEEM)
      call sum_mn_name(p%el(:,1),idiag_exm)
      call sum_mn_name(p%el(:,2),idiag_eym)
      call sum_mn_name(p%el(:,3),idiag_ezm)
      call sum_mn_name(p%sigE,idiag_sigEm)
      call sum_mn_name(p%sigB,idiag_sigBm)
      call sum_mn_name(p%eb,idiag_ebm)
      if (idiag_sigErms/=0) call sum_mn_name(p%sigE**2,idiag_sigErms,lsqrt=.true.)
      if (idiag_sigBrms/=0) call sum_mn_name(p%sigB**2,idiag_sigBrms,lsqrt=.true.)
      call sum_mn_name(p%sigE*p%e2,idiag_sigEE2m)
      call sum_mn_name(p%sigB*p%eb,idiag_sigBBEm)
      if (idiag_adphiBm/=0) then
        if (alpf/=0.) call calc_helical_term(p,gtmp,p%dphi,p%gphi,lphi_hom)
        call dot(alpf*gtmp,p%el,tmp)
        call sum_mn_name(tmp,idiag_adphiBm)
      endif
      if (idiag_Johmrms/=0) then
        call dot2_mn(p%jj_ohm,tmp)
        call sum_mn_name(tmp,idiag_Johmrms,lsqrt=.true.)
      endif
      call save_name(echarge,idiag_echarge)
      call sum_mn_name(p%e2,idiag_erms,lsqrt=.true.)
      call sum_mn_name(p%edot2,idiag_edotrms,lsqrt=.true.)
      call max_mn_name(p%e2,idiag_emax,lsqrt=.true.)
      call sum_mn_name(p%eprime**2,idiag_eprimerms,lsqrt=.true.)
      call sum_mn_name(p%bprime**2,idiag_bprimerms,lsqrt=.true.)
      call sum_mn_name(p%jprime**2,idiag_jprimerms,lsqrt=.true.)
      call sum_mn_name(p%gam_EB**2,idiag_gam_EBrms,lsqrt=.true.)
      call sum_mn_name(p%boost**2 ,idiag_boostprms,lsqrt=.true.)
      if (idiag_a0rms/=0) call sum_mn_name(p%a0**2,idiag_a0rms,lsqrt=.true.)
      call sum_mn_name(p%BcurlE,idiag_BcurlEm)
  !   if (lsolve_chargedensity) then
      call sum_mn_name(p%rhoe,idiag_rhoem)
      call sum_mn_name(p%count_eb0,idiag_count_eb0)
      call sum_mn_name(p%rhoe**2,idiag_rhoerms,lsqrt=.true.)
  !   endif
      if (idiag_divErms/=0) call sum_mn_name(p%divE**2,idiag_divErms,lsqrt=.true.)
      if(idiag_constrainteqn > 0) then
        call calc_axion_term(p,tmp,p%gphi,alpf,lphi_hom)
        call calc_constrainteqn(p,tmp,constrainteqn)
        call sum_mn_name(constrainteqn,idiag_constrainteqn)
      endif
!
      call calc_1d_diagnostics_special(p)
!
    endsubroutine calc_diagnostics_special
!******************************************************************************
    subroutine calc_1d_diagnostics_special(p)
!
!  2-D averages.
!  Note that this does not necessarily happen with ldiagnos=.true.
!
!  13-sep-25/axel: adapted from magnetic
!
      use Diagnostics
!
      type(pencil_case) :: p
!
      real, dimension(nx) :: fres2, tmp1, Rmmz, bdel2a, jdel2a
      real, dimension(nx,3) :: tmp2
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1.
!
      if (l1davgfirst .or. (ldiagnos .and. ldiagnos_need_zaverages)) then
        call yzsum_mn_name_x(p%e2, idiag_e2mx)
      endif
!
    endsubroutine calc_1d_diagnostics_special
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
      use Sub
!
!  define counters
!
      integer :: iname,inamex,inamez
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_EEEM=0; idiag_erms=0; idiag_exm=0;idiag_eym=0;  idiag_ezm=0; idiag_emax=0
        idiag_edotrms=0; idiag_a0rms=0; idiag_grms=0; idiag_da0rms=0; idiag_BcurlEm=0
        idiag_mfpf=0; idiag_fppf=0; idiag_afact=0
        idiag_rhoerms=0; idiag_divErms=0; idiag_divJrms=0
        idiag_rhoem=0; idiag_count_eb0=0; idiag_divEm=0; idiag_divJm=0; idiag_constrainteqn=0
        idiag_ebm=0; idiag_sigEm=0; idiag_sigBm=0; idiag_sigErms=0; idiag_sigBrms=0
        idiag_Johmrms=0; idiag_adphiBm=0; idiag_sigEE2m=0; idiag_sigBBEm=0
        idiag_eprimerms=0; idiag_bprimerms=0; idiag_jprimerms=0; idiag_gam_EBrms=0; 
        idiag_boostprms=0; idiag_echarge=0; idiag_e2mx=0
        cformv=''
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'EEEM',idiag_EEEM)
        call parse_name(iname,cname(iname),cform(iname),'erms',idiag_erms)
        call parse_name(iname,cname(iname),cform(iname),'exm',idiag_exm)
        call parse_name(iname,cname(iname),cform(iname),'eym',idiag_eym)
        call parse_name(iname,cname(iname),cform(iname),'ezm',idiag_ezm)
        call parse_name(iname,cname(iname),cform(iname),'eprimerms',idiag_eprimerms)
        call parse_name(iname,cname(iname),cform(iname),'bprimerms',idiag_bprimerms)
        call parse_name(iname,cname(iname),cform(iname),'jprimerms',idiag_jprimerms)
        call parse_name(iname,cname(iname),cform(iname),'gam_EBrms',idiag_gam_EBrms)
        call parse_name(iname,cname(iname),cform(iname),'boostprms',idiag_boostprms)
        call parse_name(iname,cname(iname),cform(iname),'edotrms',idiag_edotrms)
        call parse_name(iname,cname(iname),cform(iname),'emax',idiag_emax)
        call parse_name(iname,cname(iname),cform(iname),'a0rms',idiag_a0rms)
        call parse_name(iname,cname(iname),cform(iname),'grms',idiag_grms)
        call parse_name(iname,cname(iname),cform(iname),'da0rms',idiag_da0rms)
        call parse_name(iname,cname(iname),cform(iname),'BcurlEm',idiag_BcurlEm)
        call parse_name(iname,cname(iname),cform(iname),'divErms',idiag_divErms)
        call parse_name(iname,cname(iname),cform(iname),'divJrms',idiag_divJrms)
        call parse_name(iname,cname(iname),cform(iname),'rhoerms',idiag_rhoerms)
        call parse_name(iname,cname(iname),cform(iname),'divEm',idiag_divEm)
        call parse_name(iname,cname(iname),cform(iname),'divJm',idiag_divJm)
        call parse_name(iname,cname(iname),cform(iname),'rhoem',idiag_rhoem)
        call parse_name(iname,cname(iname),cform(iname),'count_eb0',idiag_count_eb0)
        call parse_name(iname,cname(iname),cform(iname),'sigEm',idiag_sigEm)
        call parse_name(iname,cname(iname),cform(iname),'sigBm',idiag_sigBm)
        call parse_name(iname,cname(iname),cform(iname),'ebm',idiag_ebm)
        call parse_name(iname,cname(iname),cform(iname),'sigErms',idiag_sigErms)
        call parse_name(iname,cname(iname),cform(iname),'sigBrms',idiag_sigBrms)
        call parse_name(iname,cname(iname),cform(iname),'sigEE2m',idiag_sigEE2m)
        call parse_name(iname,cname(iname),cform(iname),'sigBBEm',idiag_sigBBEm)
        call parse_name(iname,cname(iname),cform(iname),'adphiBm',idiag_adphiBm)
        call parse_name(iname,cname(iname),cform(iname),'Johmrms',idiag_Johmrms)
        call parse_name(iname,cname(iname),cform(iname),'echarge',idiag_echarge)
        call parse_name(iname,cname(iname),cform(iname),'mfpf',idiag_mfpf)
        call parse_name(iname,cname(iname),cform(iname),'fppf',idiag_fppf)
        call parse_name(iname,cname(iname),cform(iname),'afact',idiag_afact)
        call parse_name(iname,cname(iname),cform(iname),'constrainteqn',idiag_constrainteqn)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'e2mx',idiag_e2mx)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'exmz',idiag_exmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'eymz',idiag_eymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ezmz',idiag_ezmz)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='ee') cformv='DEFINED'
        !where(cnamev=='ee' .or. cnamev=='sigE' .or. cnamev=='sigB') cformv='DEFINED'
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!     use Poisson
!, only: inverse_laplacian
!     use Sub, only: div
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: tmp
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of electric potential
!
!  26-feb-07/axel: adapted from gross_pitaevskii
!
      use Slices_methods, only: assign_slices_vec
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (slice_data) :: slices
!
      integer :: inamev
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Electric field.
!
      case ('ee');   call assign_slices_vec (slices,f,iee)
      !case ('sigE'); call assign_slices_scal(slices,f,isigE)
      !case ('sigB'); call assign_slices_scal(slices,f,isigB)
!
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=100
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(alpf,p_par(1))
    call copy_addr(eta_ee,p_par(2))
    call copy_addr(sige_prefactor,p_par(3))
    call copy_addr(sigb_prefactor,p_par(4))
    call copy_addr(mass_chi,p_par(5))
    call copy_addr(igamma,p_par(6)) ! int
    call copy_addr(ia0,p_par(7)) ! int
    call copy_addr(idiva_name,p_par(8)) ! int
    call copy_addr(llongitudinale,p_par(9)) ! bool
    call copy_addr(llorenz_gauge_disp,p_par(10)) ! bool
    call copy_addr(lphi_hom,p_par(11)) ! bool
    call copy_addr(lnoncollinear_eb,p_par(12)) ! bool
    call copy_addr(lnoncollinear_eb_aver,p_par(13)) ! bool
    call copy_addr(lcollinear_eb,p_par(14)) ! bool
    call copy_addr(lcollinear_eb_aver,p_par(15)) ! bool
    call copy_addr(leedot_as_aux,p_par(16)) ! bool
    call copy_addr(lcurlya,p_par(17)) ! bool
    call copy_addr(lsolve_chargedensity,p_par(18)) ! bool
    call copy_addr(ldive_as_aux,p_par(19)) ! bool
    call copy_addr(lsige_as_aux,p_par(20)) ! bool
    call copy_addr(lsigb_as_aux,p_par(21)) ! bool
    call copy_addr(lallow_bprime_zero,p_par(22)) ! bool
    call copy_addr(lswitch_off_divj,p_par(23)) ! bool
    call copy_addr(lswitch_off_gamma,p_par(24)) ! bool
    call copy_addr(lmass_suppression,p_par(25)) ! bool
    call copy_addr(beta_inflation,p_par(26))
    call copy_addr(c_light2,p_par(27))
    call copy_addr(idiag_bcurlem,p_par(28)) ! int
    call copy_addr(idiag_adphibm,p_par(29)) ! int
    call copy_addr(idiag_johmrms,p_par(30)) ! int
    call copy_addr(lapply_gamma_corr,p_par(31)) ! bool
    call copy_addr(lphi_linear_regime,p_par(32)) ! bool
    call copy_addr(weight_longitudinale,p_par(33))
    call copy_addr(alpfpsi,p_par(34))
    call copy_addr(coupl_gy,p_par(35))
    call copy_addr(lpsi_hom,p_par(36)) ! bool
    call copy_addr(iex,p_par(37)) ! int
    call copy_addr(iey,p_par(38)) ! int
    call copy_addr(iez,p_par(39)) ! int
    call copy_addr(iedotx,p_par(40)) ! int
    call copy_addr(iedoty,p_par(41)) ! int
    call copy_addr(iedotz,p_par(42)) ! int
    call copy_addr(idive,p_par(43)) ! int
    call copy_addr(isige,p_par(44)) ! int
    call copy_addr(isigb,p_par(45)) ! int
    call copy_addr(irhoe,p_par(46)) ! int



    endsubroutine pushpars2c
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
!
endmodule Special
