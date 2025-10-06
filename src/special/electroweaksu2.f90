! $Id$
!
!  Dynamical equations for electroweak SU(2) non-Abelian gauge fields
!
!  12-aug-25/alberto: adapted from disp_current.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 21
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED W1(3); W2(3); W3(3)
! PENCILS PROVIDED dW1(3); dW2(3); dW3(3)
! PENCILS PROVIDED GammaW1; GammaW2; GammaW3
! PENCILS PROVIDED W1sq; W2sq; W3sq
! PENCILS PROVIDED dW1sq; dW2sq; dW3sq
! PENCILS PROVIDED rhoW1; rhoW2; rhoW3
! PENCILS PROVIDED W1ddot_gw(3); W2ddot_gw(3); W3ddot_gw(3)
! PENCILS PROVIDED W1ddot_gw2(3); W2ddot_gw2(3); W3ddot_gw2(3)
! PENCILS PROVIDED divdotW1; divdotW2, divdotW3
! PENCILS PROVIDED del2W1(3); del2W2(3); del2W3(3)
! PENCILS PROVIDED W1ij(3,3); W2ij(3,3); W3ij(3,3)
! PENCILS PROVIDED curlBW1(3); curlBW2(3); curlBW3(3)
! PENCILS PROVIDED gGammaW1(3); gGammaW2(3); gGammaW3(3)
! PENCILS PROVIDED jj_higgsW1(3); jj_higgsW2(3); jj_higgsW3(3)
! PENCILS PROVIDED rhoe_higgsW1; rhoe_higgsW2; rhoe_higgsW3
! PENCILS PROVIDED W1dotW1; W1dotW2; W1dotW3
! PENCILS PROVIDED W2dotW1; W2dotW2; W2dotW3
! PENCILS PROVIDED W3dotW1; W3dotW2; W3dotW3
! PENCILS PROVIDED W1ddotsq, W2ddotsq, W3ddotsq
! PENCILS EXPECTED phi, dphi, cov_der(4,4), phi_doublet(3)
!***************************************************************
!
!
!  compulsory pencils
!
!
module Special
!
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
  integer :: iGammaW=0, iW0=0, irhoW=0
  integer :: idWW1=0, idWW2=0, idWW3=0, iWW1ddot=0, iWW2ddot=0
  integer :: iWW3ddot=0, idWW=0
  ! Wddot = d^2 W / dt^2 ~ dt E
  integer :: idivdotW=0 ! -div dt W (=div E) in Weyl gauge
  real, dimension (ninit, 3) :: amplWW=0.0, ampldWW=0.0
  real, dimension(3) :: ampl_Wx=0.0, ampl_Wy=0.0, ampl_Wz=0.0, ampl_W0=0.0
  real, dimension(3) :: ampl_dWx=0.0, ampl_dWy=0.0, ampl_dWz=0.0, ampl_dW0=0.0
  real, dimension(3) :: kx_Wx=0.0, kx_Wy=0.0, kx_Wz=0.0
  real, dimension(3) :: ky_Wx=0.0, ky_Wy=0.0, ky_Wz=0.0
  real, dimension(3) :: kz_Wx=0.0, kz_Wy=0.0, kz_Wz=0.0
  real, dimension(3) :: kx_dWx=0.0, kx_dWy=0.0, kx_dWz=0.0
  real, dimension(3) :: ky_dWx=0.0, ky_dWy=0.0, ky_dWz=0.0
  real, dimension(3) :: kz_dWx=0.0, kz_dWy=0.0, kz_dWz=0.0
  real, dimension(3) :: phase_Wx=0.0, phase_Wy=0.0, phase_Wz=0.0
  real, dimension(3) :: phase_dWx=0.0, phase_dWy=0.0, phase_dWz=0.0
  real, dimension(3) :: initpower_WW=0.0, initpower2_WW=0.0
  real, dimension(3) :: initpower_dWW=0.0, initpower2_dWW=0.0
  real, dimension(3) :: cutoff_WW=0.0, ncutoff_WW=0.0, kpeak_WW=0.0
  real, dimension(3) :: cutoff_dWW=0.0, ncutoff_dWW=0.0, kpeak_dWW=0.0
  real, dimension(3) :: relhel_WW=0.0, kgaussian_WW=0.0
  real, dimension(3) :: relhel_dWW=0.0, kgaussian_dWW=0.0
  real :: weight_longitudinalW=2.0, eta_WW=0.0
!   logical :: luse_scale_factor_in_sigma=.false., lapply_Gamma_corr=.true.
!   logical, pointer :: lohm_evolve
  real :: coupl_gw=.65    ! electroweak SU(2) x U(1) coupling of Higgs to SU(2)
  logical, pointer :: lphi_doublet, lphi_weakcharge, lphi_hom
  logical :: llongitudinalW=.true., lskip_projection_WW=.false.
  logical :: lWddot_as_aux=.false.,ldivdotW_as_aux=.false. !, lsolve_chargedensityW=.false.
  logical :: lscale_tobox=.true., lpower_profile_file=.false.
  logical :: lno_noise_WW=.false., lno_noise_dWW=.false.
  logical :: lfixed_phase_WW=.false., lfixed_phase_dWW=.false.
  logical :: lskip_projection=.false., lvectorpotential=.false.
  character(len=labellen) :: initWW0='nothing', initdWW0='nothing'
  character (len=labellen), dimension(ninit) :: initWW='nothing'
  character (len=labellen), dimension(ninit) :: initdWW='nothing'
  character (len=labellen) :: power_filename='power_profile.dat'
!
  namelist /special_init_pars/ &
    initWW, initWW0, ampl_Wx, ampl_Wy, ampl_Wz, ampl_W0, &
    initdWW, initdWW0, ampl_dWx, ampl_dWy, ampl_dWz, ampl_dW0, &
    kx_Wx, kx_Wy, kx_Wz, &
    ky_Wx, ky_Wy, ky_Wz, &
    kz_Wx, kz_Wy, kz_Wz, &
    kx_dWx, kx_dWy, kx_dWz, &
    ky_dWx, ky_dWy, ky_dWz, &
    kz_dWx, kz_dWy, kz_dWz, &
    phase_Wx, phase_Wy, phase_Wz, phase_dWx, phase_dWy, phase_dWz, &
    ! lsolve_chargedensityW, &
    llongitudinalW, amplWW, initpower_WW, initpower2_WW, lscale_tobox, &
    initpower_dWW, initpower2_dWW, &
    cutoff_dWW, ncutoff_dWW, kpeak_dWW, relhel_dWW, kgaussian_dWW, &
    lno_noise_WW, weight_longitudinalW, coupl_gw, lWddot_as_aux, &
    lfixed_phase_WW, lskip_projection_WW, &
    lpower_profile_file, power_filename, ldivdotW_as_aux
!
  ! run parameters
  logical :: reinitialize_WW=.false., reinitialize_dWW=.false.
  real, dimension(3) :: rescale_WW=1.0, rescale_dWW=1.0
  namelist /special_run_pars/ &
    llongitudinalW, eta_WW, weight_longitudinalW, &
    lWddot_as_aux, ldivdotW_as_aux, reinitialize_WW, initWW, &
    rescale_WW, rescale_dWW, initdWW, reinitialize_dWW
!
! Declare any index variables necessary for main or
!
  real :: c_light2
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_W1rms=0       ! DIAG_DOC: $\left<{\Wv^1}^2\right>^{1/2}$
  integer :: idiag_W2rms=0       ! DIAG_DOC: $\left<{\Wv^2}^2\right>^{1/2}$
  integer :: idiag_W3rms=0       ! DIAG_DOC: $\left<{\Wv^3}^2\right>^{1/2}$
  integer :: idiag_W1max=0       ! DIAG_DOC: $\max(|\Wv^1|)$
  integer :: idiag_W2max=0       ! DIAG_DOC: $\max(|\Wv^2|)$
  integer :: idiag_W3max=0       ! DIAG_DOC: $\max(|\Wv^3|)$
  integer :: idiag_dW1rms=0       ! DIAG_DOC: $\left<{{\dot{\Wv}}_1}^2\right>^{1/2}$
  integer :: idiag_dW2rms=0       ! DIAG_DOC: $\left<{{\dot{\Wv}}_2}^2\right>^{1/2}$
  integer :: idiag_dW3rms=0       ! DIAG_DOC: $\left<{{\dot{\Wv}}_3}^2\right>^{1/2}$
  integer :: idiag_dW1max=0       ! DIAG_DOC: $\max(|\dot{\Wv}_1|)$
  integer :: idiag_dW2max=0       ! DIAG_DOC: $\max(|\dot{\Wv}_2|)$
  integer :: idiag_dW3max=0       ! DIAG_DOC: $\max(|\dot{\Wv}_3|)$
  integer :: idiag_W1ddotrms=0.  ! DIAG_DOC: $\left<{\ddot{\Wv}_1}^2\right>^{1/2}$
  integer :: idiag_W2ddotrms=0.  ! DIAG_DOC: $\left<{\ddot{\Wv}_2}^2\right>^{1/2}$
  integer :: idiag_W3ddotrms=0.  ! DIAG_DOC: $\left<{\ddot{\Wv}_3}^2\right>^{1/2}$
  integer :: idiag_divW1rms=0    ! DIAG_DOC: $\left<{\nab{\Wv}_1}^2\right>^{1/2}$
  integer :: idiag_divW2rms=0    ! DIAG_DOC: $\left<{\nab{\Wv}_2}^2\right>^{1/2}$
  integer :: idiag_divW3rms=0    ! DIAG_DOC: $\left<{\nab{\Wv}_3}^2\right>^{1/2}$
  integer :: idiag_divW1m=0      ! DIAG_DOC: $\left<\nab{\Wv}_1\right>$
  integer :: idiag_divW2m=0      ! DIAG_DOC: $\left<\nab{\Wv}_2\right>$
  integer :: idiag_divW3m=0      ! DIAG_DOC: $\left<\nab{\Wv}_3\right>$
  integer :: idiag_divdotW1rms=0 ! DIAG_DOC: $\left<{\nab\dot{\Wv}_1}^2\right>^{1/2}$
  integer :: idiag_divdotW2rms=0 ! DIAG_DOC: $\left<{\nab\dot{\Wv}_2}^2\right>^{1/2}$
  integer :: idiag_divdotW3rms=0 ! DIAG_DOC: $\left<{\nab\dot{\Wv}_3}^2\right>^{1/2}$
  integer :: idiag_rhoW1rms=0   ! DIAG_DOC: $\left<{\rho_{\Wv_1}}^2\right>^{1/2}$
  integer :: idiag_rhoW2rms=0   ! DIAG_DOC: $\left<{\rho_{\Wv_2}}^2\right>^{1/2}$
  integer :: idiag_rhoW3rms=0   ! DIAG_DOC: $\left<{\rho_{\Wv_3}}^2\right>^{1/2}$
  integer :: idiag_divdotW1m=0   ! DIAG_DOC: $\left<\nab\dot\Wv_1\right>$
  integer :: idiag_divdotW2m=0   ! DIAG_DOC: $\left<\nab\dot\Wv_2\right>$
  integer :: idiag_divdotW3m=0   ! DIAG_DOC: $\left<\nab\dot\Wv_3\right>$
  integer :: idiag_rhoW1m=0     ! DIAG_DOC: $\left<\rho_e\Wv_1\right>$
  integer :: idiag_rhoW2m=0     ! DIAG_DOC: $\left<\rho_e\Wv_2\right>$
  integer :: idiag_rhoW3m=0     ! DIAG_DOC: $\left<\rho_e\Wv_3\right>$
  integer :: idiag_constrainteqnW=0  ! DIAG_DOC: $<deldotW+>$
  integer :: idiag_W1xm=0        ! DIAG_DOC: $\left<W_x^1\right>$
  integer :: idiag_W1ym=0        ! DIAG_DOC: $\left<W_y^1\right>$
  integer :: idiag_W1zm=0        ! DIAG_DOC: $\left<W_z^1\right>$
  integer :: idiag_W2xm=0        ! DIAG_DOC: $\left<W_x^2\right>$
  integer :: idiag_W2ym=0        ! DIAG_DOC: $\left<W_y^2\right>$
  integer :: idiag_W2zm=0        ! DIAG_DOC: $\left<W_z^2\right>$
  integer :: idiag_W3xm=0        ! DIAG_DOC: $\left<W_x^3\right>$
  integer :: idiag_W3ym=0        ! DIAG_DOC: $\left<W_y^3\right>$
  integer :: idiag_W3zm=0        ! DIAG_DOC: $\left<W_z^3\right>$
  integer :: idiag_dW1xm=0       ! DIAG_DOC: $\left<\dot{W}^1_x\right>$
  integer :: idiag_dW1ym=0       ! DIAG_DOC: $\left<\dot{W}^1_y\right>$
  integer :: idiag_dW1zm=0       ! DIAG_DOC: $\left<\dot{W}^1_z\right>$
  integer :: idiag_dW2xm=0       ! DIAG_DOC: $\left<\dot{W}^2_x\right>$
  integer :: idiag_dW2ym=0       ! DIAG_DOC: $\left<\dot{W}^2_y\right>$
  integer :: idiag_dW2zm=0       ! DIAG_DOC: $\left<\dot{W}^2_z\right>$
  integer :: idiag_dW3xm=0       ! DIAG_DOC: $\left<\dot{W}^3_x\right>$
  integer :: idiag_dW3ym=0       ! DIAG_DOC: $\left<\dot{W}^3_y\right>$
  integer :: idiag_dW3zm=0       ! DIAG_DOC: $\left<\dot{W}^3_z\right>$
  integer :: idiag_W1dotW1m=0     ! DIAG_DOC: $\left<\dot{\Wv}_1\cdot{\Wv}_1\right>$
  integer :: idiag_W2dotW2m=0     ! DIAG_DOC: $\left<\dot{\Wv}_2\cdot{\Wv}_2\right>$
  integer :: idiag_W3dotW3m=0     ! DIAG_DOC: $\left<\dot{\Wv}_3\cdot{\Wv}_3\right>$
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_W1xmz=0       ! XYAVG_DOC: $\left<{\cal W}_x^1\right>_{xy}$
  integer :: idiag_W1ymz=0       ! XYAVG_DOC: $\left<{\cal W}_y^1\right>_{xy}$
  integer :: idiag_W1zmz=0       ! XYAVG_DOC: $\left<{\cal W}_z^1\right>_{xy}$
  integer :: idiag_W2xmz=0       ! XYAVG_DOC: $\left<{\cal W}_x^2\right>_{xy}$
  integer :: idiag_W2ymz=0       ! XYAVG_DOC: $\left<{\cal W}_y^2\right>_{xy}$
  integer :: idiag_W2zmz=0       ! XYAVG_DOC: $\left<{\cal W}_z^2\right>_{xy}$
  integer :: idiag_W3xmz=0       ! XYAVG_DOC: $\left<{\cal W}_x^3\right>_{xy}$
  integer :: idiag_W3ymz=0       ! XYAVG_DOC: $\left<{\cal W}_y^3\right>_{xy}$
  integer :: idiag_W3zmz=0       ! XYAVG_DOC: $\left<{\cal W}_z^3\right>_{xy}$
  integer :: idiag_dW1xmz=0      ! XYAVG_DOC: $\left<\dot{\cal W}_x^1\right>_{xy}$
  integer :: idiag_dW1ymz=0      ! XYAVG_DOC: $\left<\dot{\cal W}_y^1\right>_{xy}$
  integer :: idiag_dW1zmz=0      ! XYAVG_DOC: $\left<\dot{\cal W}_z^1\right>_{xy}$
  integer :: idiag_dW2xmz=0      ! XYAVG_DOC: $\left<\dot{\cal W}_x^2\right>_{xy}$
  integer :: idiag_dW2ymz=0      ! XYAVG_DOC: $\left<\dot{\cal W}_y^2\right>_{xy}$
  integer :: idiag_dW2zmz=0      ! XYAVG_DOC: $\left<\dot{\cal W}_z^2\right>_{xy}$
  integer :: idiag_dW3xmz=0      ! XYAVG_DOC: $\left<\dot{\cal W}_x^3\right>_{xy}$
  integer :: idiag_dW3ymz=0      ! XYAVG_DOC: $\left<\dot{\cal W}_y^3\right>_{xy}$
  integer :: idiag_dW3zmz=0      ! XYAVG_DOC: $\left<\dot{\cal W}_z^3\right>_{xy}$
!
  contains
!
!***********************************************************************
    subroutine register_special
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  12-aug-25/alberto: coded equations of motion for SU(2) x U(1)
!                     electroweak gauge fields
!
      use FArrayManager
      use Sub, only: register_report_aux
      use SharedVariables, only: put_shared_variable
      integer :: j
!
! alberto: non-Abelian SU(2) (Wia) gauge fields in the W_0^i = 0 gauge
!
      call farray_register_pde('WW1',iWW1,vector=3)
      iWW=iWW1
      call farray_register_pde('WW2',iWW2,vector=3)
      call farray_register_pde('WW3',iWW3,vector=3)
      call farray_register_pde('dWW1',idWW1,vector=3)
      idWW=idWW1
      call farray_register_pde('dWW2',idWW2,vector=3)
      call farray_register_pde('dWW3',idWW3,vector=3)

      if (lWddot_as_aux) &
        call farray_register_auxiliary('WW1ddot', iWW1ddot, array=3)
        call farray_register_auxiliary('WW2ddot', iWW2ddot, array=3)
        call farray_register_auxiliary('WW3ddot', iWW3ddot, array=3)

      if (ldivdotW_as_aux) &
        call farray_register_auxiliary('divdotW', idivdotW, array=3)

      ! if(lsolve_chargedensityW) &
      !   call farray_register_pde('rhoW',irhoW,vector=3)
      !
      if (llongitudinalW) &
        call farray_register_pde('GammaW',iGammaW,vector=3)
!
!  The following variables are also used in special/klein_gordon.f90
!
      call put_shared_variable('llongitudinalW',llongitudinalW)
      call put_shared_variable('coupl_gw',coupl_gw,caller='register_electroweaksu2')
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
      integer :: i, j
!
!  Initialize module variables which are parameter dependent
!  If one really wants to work with c_light /= 1,
!  then one needs to override this.
!
      if (c_light/=1.) call fatal_error('electroweaksu2', "use unit_system='set'")
      c_light2=c_light**2

      ! if (lmagnetic .and. .not.lswitch_off_divJ) &
      !   call get_shared_variable('eta',eta, caller='initialize_magnetic')

!  The following are only obtained when luse_scale_factor_in_sigma=T
!  (luse_scale_factor_in_sigma=F by default, because they are defined
!  in special/backreact_infl.f90, which may not be always be used).
!
      ! if (luse_scale_factor_in_sigma) then
      !   call get_shared_variable('Hscript', Hscript)
      !   call get_shared_variable('echarge', echarge)
      !   call get_shared_variable('sigEm_all', sigEm_all)
      !   call get_shared_variable('sigBm_all', sigBm_all)
      !   call get_shared_variable('lohm_evolve', lohm_evolve)
      ! else
      !   if (.not.associated(Hscript)) allocate(Hscript,echarge,sigEm_all,sigBm_all)
      !   Hscript=0.
      !   echarge=0.
      !   sigEm_all=0.
      !   sigBm_all=0.
      !   allocate(lohm_evolve)
      !   lohm_evolve=.false.
      ! endif
!
!  Reinitialize gauge field W and time derivative using a small selection of perturbations
!  that were mostly also available as initial conditions.
!
      if (reinitialize_WW) then
        do j=1,ninit
          select case (initWW(j))
          case ('rescale')
            do i=0,2
              f(:,:,:,iWW+i*3:iWW+i*3+2)=rescale_WW(i+1)*f(:,:,:,iWW+i*3:iWW+i*3+2)
              if (llongitudinalW) &
                  f(:,:,:,iGammaW+i)=rescale_WW(i+1)*f(:,:,:,iGammaW+i)
            enddo
          case ('gaussian-noise')
            do i=0,2
              call gaunoise(amplWW(j,i+1),f,iWW+i*3,iWW+i*3+2)
            enddo
          case default
          endselect
        enddo
      endif
      if (reinitialize_dWW) then
        do j=1,ninit
          select case (initdWW(j))
          case ('rescale')
            do i=0,2
              f(:,:,:,idWW+i*3:idWW+i*3+2)=rescale_dWW(i+1)*f(:,:,:,idWW+i*3:idWW+i*3+2)
            enddo
          case ('gaussian-noise')
            do i=0,2
              call gaunoise(ampldWW(j, i+1),f,idWW+i*3,idWW+i*3+2)
            enddo
          case default
          endselect
        enddo
      endif

      if (lklein_gordon) then
        call get_shared_variable('lphi_doublet',lphi_doublet, caller='initialize_electroweaksu2')
        call get_shared_variable('lphi_weakcharge',lphi_weakcharge, &
          caller='initialize_electroweaksu2')
      else
        if (.not.associated(lphi_doublet)) allocate(lphi_doublet,lphi_weakcharge)
        lphi_doublet=.false.
        lphi_weakcharge=.false.
      endif

      if (iex /= 0) then
        call get_shared_variable('lphi_hom',lphi_hom, caller='initialize_electroweaksu2')
      else
        if(.not.associated(lphi_hom)) allocate(lphi_hom)
        lphi_hom=.false.
      endif

      if (lphi_hom) weight_longitudinalW=0.0
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
      real, dimension (nx) :: divW, divE
      integer :: j, i
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      do j=1,ninit
        select case (initWW(j))
          case ('nothing'); if (lroot) print*,'initWW: nothing'
          case ('zero'); f(:,:,:,iWW1:iWW3+2)=0.
          case ('coswave-phase')
            do i=0,2
              call coswave_phase(f,iWW+3*i,ampl_Wx(i+1),kx_Wx(i+1),ky_Wx(i+1),kz_Wx(i+1),phase_Wx(i+1))
              call coswave_phase(f,iWW+3*i+1,ampl_Wy(i+1),kx_Wy(i+1),ky_Wy(i+1),kz_Wy(i+1),phase_Wy(i+1))
              call coswave_phase(f,iWW+3*i+2,ampl_Wz(i+1),kx_Wz(i+1),ky_Wz(i+1),kz_Wz(i+1),phase_Wz(i+1))
            enddo
          case ('sinwave-phase')
            do i=0,2
              call sinwave_phase(f,iWW+3*i,ampl_Wx(i+1),kx_Wx(i+1),ky_Wx(i+1),kz_Wx(i+1),phase_Wx(i+1))
              call sinwave_phase(f,iWW+3*i+1,ampl_Wy(i+1),kx_Wy(i+1),ky_Wy(i+1),kz_Wy(i+1),phase_Wy(i+1))
              call sinwave_phase(f,iWW+3*i+2,ampl_Wz(i+1),kx_Wz(i+1),ky_Wz(i+1),kz_Wz(i+1),phase_Wz(i+1))
            enddo
          case ('power_randomphase_hel')
            do i=1,3
              call power_randomphase_hel(amplWW(j,i+1),initpower_WW(i),initpower2_WW(i), &
                  cutoff_WW(i),ncutoff_WW(i),kpeak_WW(i),f,iWW+3*(i-1),iWW+3*i-1,relhel_WW(i),kgaussian_WW(i), &
                  lskip_projection, lvectorpotential, &
                  lscale_tobox=lscale_tobox, lpower_profile_file=lpower_profile_file, &
                  power_filename=power_filename, &
                  lno_noise=lno_noise_WW, lfixed_phase=lfixed_phase_WW)
            enddo
          case default
            call fatal_error('init_WW','no such init_WW: "'//trim(initWW(j))//'"')
        endselect
      enddo

      do j=1,ninit
        select case (initdWW(j))
          case ('nothing'); if (lroot) print*,'initdWW: nothing'
          case ('zero'); f(:,:,:,idWW1:idWW3+2)=0.
          case ('coswave-phase')
            do i=0,2
              call coswave_phase(f,idWW+3*i,ampl_dWx(i+1),kx_dWx(i+1),ky_dWx(i+1),kz_dWx(i+1),phase_dWx(i+1))
              call coswave_phase(f,idWW+3*i+1,ampl_dWy(i+1),kx_dWy(i+1),ky_dWy(i+1),kz_dWy(i+1),phase_dWy(i+1))
              call coswave_phase(f,idWW+3*i+2,ampl_dWz(i+1),kx_dWz(i+1),ky_dWz(i+1),kz_dWz(i+1),phase_dWz(i+1))
            enddo
          case ('sinwave-phase')
            do i=0,2
              call sinwave_phase(f,idWW+3*i,ampl_dWx(i+1),kx_dWx(i+1),ky_dWx(i+1),kz_dWx(i+1),phase_dWx(i+1))
              call sinwave_phase(f,idWW+3*i+1,ampl_dWy(i+1),kx_dWy(i+1),ky_dWy(i+1),kz_dWy(i+1),phase_dWy(i+1))
              call sinwave_phase(f,idWW+3*i+2,ampl_dWz(i+1),kx_dWz(i+1),ky_dWz(i+1),kz_dWz(i+1),phase_dWz(i+1))
            enddo
          case ('power_randomphase_hel')
            do i=1,3
              call power_randomphase_hel(ampldWW(j,i+1),initpower_dWW(i),initpower2_dWW(i), &
                  cutoff_dWW(i),ncutoff_dWW(i),kpeak_dWW(i),f,idWW+3*(i-1),idWW+3*i-1,relhel_dWW(i),kgaussian_dWW(i), &
                  lskip_projection, lvectorpotential, &
                  lscale_tobox=lscale_tobox, lpower_profile_file=lpower_profile_file, &
                  power_filename=power_filename, &
                  lno_noise=lno_noise_dWW, lfixed_phase=lfixed_phase_dWW)
            enddo
          case default
            call fatal_error('init_WW','no such init_WW: "'//trim(initWW(j))//'"')
        endselect
      enddo
!
!  Initial condition for GammaW^a = div W^a
!
      if (llongitudinalW) then
        do n=n1,n2; do m=m1,m2; do i=0,2
          call div(f,iWW+3*i,divW)
          f(l1:l2,m,n,iGammaW+i)=divW
        enddo; enddo; enddo
      endif
!
!  Initial condition for rhoW = div dt W
!
      ! if (lsolve_chargedensityW) then
      !   do n=n1,n2; do m=m1,m2; do i=0,2
      !     call div(f,idWW+3*i,divE)
      !     f(l1:l2,m,n,irhoW+i)=divE
      !   enddo; enddo; enddo
      ! endif
!
!  Initialize diva_name if llorenz_gauge_disp=T
!
!       if (llorenz_gauge_disp) then
!         do n=n1,n2; do m=m1,m2
!           call div(f,iaa,diva)
!           f(l1:l2,m,n,idiva_name)=diva
!         enddo; enddo
! !
! !  initial conditions for A0 (provided llorenz_gauge_disp=T)
! !
!         select case (inita0)
!           case ('coswave-phase')
!             call coswave_phase(f,ia0,ampl_a0,kx_a0,ky_a0,kz_a0,phase_a0)
!           case ('zero'); f(:,:,:,ia0)=0.
!           case ('power_randomphase')
!             call power_randomphase_hel(ampla0,initpower_a0,initpower2_a0, &
!               cutoff_a0,ncutoff_a0,kpeak_a0,f,ia0,ia0, &
!               relhel_a0,kgaussian_a0, lskip_projection_a0, lvectorpotential, &
!               lscale_tobox, lpower_profile_file=.false.)
!           case default
!             call fatal_error("init_special","no such inita0: "//trim(inita0))
!         endselect
!       endif

      if (lWddot_as_aux) f(:,:,:,iWW1ddot:iWW3ddot+2)=0.
      if (ldivdotW_as_aux) f(:,:,:,idivdotW:idivdotW+2)=0.
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  25-feb-07/axel: adapted
!  12-aug-25/alberto: adapted for SU(2)xU(1) electroweak gauge fields
!
!  compulsory pencils
!
!   Mandatory pencils
      ! lpenc_requested(i_W1)=.true.
      ! lpenc_requested(i_W2)=.true.
      ! lpenc_requested(i_W3)=.true.
      ! lpenc_requested(i_dW1)=.true.
      ! lpenc_requested(i_dW2)=.true.
      ! lpenc_requested(i_dW3)=.true.
      ! lpenc_requested(i_gGammaW1)=.true.
      ! lpenc_requested(i_gGammaW2)=.true.
      ! lpenc_requested(i_gGammaW3)=.true.
      ! lpenc_requested(i_curlBW1)=.true.
      ! lpenc_requested(i_curlBW2)=.true.
      ! lpenc_requested(i_curlBW3)=.true.
      ! lpenc_requested(i_del2W1)=.true.
      ! lpenc_requested(i_del2W2)=.true.
      ! lpenc_requested(i_del2W3)=.true.
      ! lpenc_requested(i_W1sq)=.true.
      ! lpenc_requested(i_W2sq)=.true.
      ! lpenc_requested(i_W3sq)=.true.
!
!  Terms for Gamma evolution.
!
      if (llongitudinalW) then
        lpenc_requested(i_divdotW1)=.true.
        lpenc_requested(i_divdotW2)=.true.
        lpenc_requested(i_divdotW3)=.true.
        lpenc_diagnos(i_W1dotW1)=.true.
        lpenc_diagnos(i_W1dotW2)=.true.
        lpenc_diagnos(i_W1dotW3)=.true.
        lpenc_diagnos(i_W2dotW1)=.true.
        lpenc_diagnos(i_W2dotW2)=.true.
        lpenc_diagnos(i_W2dotW3)=.true.
        lpenc_diagnos(i_W3dotW1)=.true.
        lpenc_diagnos(i_W3dotW2)=.true.
        lpenc_diagnos(i_W3dotW3)=.true.
      endif
!
      if (idiag_divW1m/=0. .or. idiag_divW1rms/=0.) then
        lpenc_requested(i_GammaW1)=.true.
      endif
      if (idiag_divW2m/=0. .or. idiag_divW2rms/=0.) then
        lpenc_requested(i_GammaW2)=.true.
      endif
      if (idiag_divW3m/=0. .or. idiag_divW3rms/=0.) then
        lpenc_requested(i_GammaW3)=.true.
      endif

      if (idiag_divdotW1m/=0. .or. idiag_divdotW1rms/=0.) then
        lpenc_requested(i_divdotW1)=.true.
      endif
      if (idiag_divdotW2m/=0. .or. idiag_divdotW2rms/=0.) then
        lpenc_requested(i_divdotW2)=.true.
      endif
      if (idiag_divdotW3m/=0. .or. idiag_divdotW3rms/=0.) then
        lpenc_requested(i_divdotW3)=.true.
      endif
!
!  charge density
!
      ! if (lsolve_chargedensityW) then
      !   lpenc_requested(i_divJ)=.true.
      !   lpenc_requested(i_uij)=.true.
      !   lpenc_requested(i_bij)=.true.
      !   lpenc_requested(i_uu)=.true.
      !   lpenc_requested(i_bb)=.true.
      ! endif
!
!  diffusion term.
!
      if (eta_WW/=0.) then
        lpenc_requested(i_del2W1)=.true.
        lpenc_requested(i_del2W2)=.true.
        lpenc_requested(i_del2W3)=.true.
      endif
!
!  Higgs pencils
!
      if (lklein_gordon .and. lphi_doublet .and. lphi_weakcharge) then
        lpenc_requested(i_phi)=.true.
        lpenc_requested(i_phi_doublet)=.true.
        lpenc_requested(i_cov_der)=.true.
      endif
!
!  Diagnostics pencils:
!

      ! if (idiag_BcurlEm/=0) then
      !   lpenc_diagnos(i_curlE)=.true.
      !   lpenc_diagnos(i_BcurlE)=.true.
      ! endif
!
      ! if (idiag_ebm/=0) lpenc_diagnos(i_eb)=.true.
      ! if (idiag_a0rms/=0) lpenc_diagnos(i_a0)=.true.
      ! if (idiag_grms/=0) lpenc_diagnos(i_diva)=.true.
      ! if (idiag_edotrms/=0) lpenc_diagnos(i_edot2)=.true.
      if (idiag_W1rms/=0 .or. idiag_W1max/=0) lpenc_diagnos(i_W1sq)=.true.
      if (idiag_W2rms/=0 .or. idiag_W2max/=0) lpenc_diagnos(i_W2sq)=.true.
      if (idiag_W3rms/=0 .or. idiag_W3max/=0) lpenc_diagnos(i_W3sq)=.true.
      if (idiag_dW1rms/=0 .or. idiag_dW1max/=0) lpenc_diagnos(i_dW1sq)=.true.
      if (idiag_dW2rms/=0 .or. idiag_dW2max/=0) lpenc_diagnos(i_dW2sq)=.true.
      if (idiag_dW3rms/=0 .or. idiag_dW3max/=0) lpenc_diagnos(i_dW3sq)=.true.
      if (idiag_W1ddotrms/=0) lpenc_diagnos(i_W1ddotsq)=.true.
      if (idiag_W2ddotrms/=0) lpenc_diagnos(i_W2ddotsq)=.true.
      if (idiag_W3ddotrms/=0) lpenc_diagnos(i_W3ddotsq)=.true.
      if (idiag_W1dotW1m /= 0) lpenc_diagnos(i_W1dotW1)=.true.
      if (idiag_W2dotW2m /= 0) lpenc_diagnos(i_W2dotW2)=.true.
      if (idiag_W3dotW3m /= 0) lpenc_diagnos(i_W3dotW3)=.true.
      ! if (idiag_Wxmz/=0 .or. idiag_Wymz/=0 .or. idiag_Wzmz/=0 ) lpenc_diagnos(i_WW)=.true.
      ! if (idiag_Wxm/=0 .or. idiag_Wym/=0 .or. idiag_Wzm/=0 ) lpenc_diagnos(i_WW)=.true.
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
      ! use Sub, only: grad, div, curl, del2v, dot2_mn, dot, levi_civita, gij, &
      !                del2v_etc, dot_mn_sv_pencil
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p

      real, dimension (nx,3) :: tmp1, tmp2
      real, dimension (nx) :: tmp0
!
      ! real, dimension (nx) :: tmp
      integer :: i,j,k
!
      intent(inout) :: f
      intent(inout) :: p
!
!  Pencil for charge density.
!
      ! if (lsolve_chargedensityW) then
      !   p%rhoW1=f(l1:l2,m,n,irhoW)
      !   p%rhoW2=f(l1:l2,m,n,irhoW+1)
      !   p%rhoW3=f(l1:l2,m,n,irhoW+2)
      ! endif
!
!  Terms for Gamma evolution (required if llongitudinalW=T)
!
      if (lpenc_requested(i_divdotW1)) call div(f,idWW1,p%divdotW1)
      if (lpenc_requested(i_divdotW2)) call div(f,idWW2,p%divdotW2)
      if (lpenc_requested(i_divdotW3)) call div(f,idWW3,p%divdotW3)
!
!  GammaW1, GammaW2, GammaW3 pencils (required from klein_gordon if lphi_weakcharge=T)
!
      if (lpenc_requested(i_GammaW1)) then
        if (llongitudinalW) then
          p%GammaW1=f(l1:l2,m,n,iGammaW)
        else
          call div(f,iWW1,p%GammaW1)
        endif
      endif
      if (lpenc_requested(i_GammaW2)) then
        if (llongitudinalW) then
          p%GammaW2=f(l1:l2,m,n,iGammaW+1)
        else
          call div(f,iWW2,p%GammaW2)
        endif
      endif
      if (lpenc_requested(i_GammaW3)) then
        if (llongitudinalW) then
          p%GammaW3=f(l1:l2,m,n,iGammaW+2)
        else
          call div(f,iWW3,p%GammaW3)
        endif
      endif
!
!  Pencils del2W1, del2W2, del2W3 and W1ij, W2ij, W3ij
!  mandatory for curlBW1, curlBW2, curlBW3
!
      call gij(f,iWW1,p%W1ij,1)
      call gij(f,iWW2,p%W2ij,1)
      call gij(f,iWW3,p%W3ij,1)

      call del2v(f,iWW1,p%del2W1)
      call del2v(f,iWW2,p%del2W2)
      call del2v(f,iWW3,p%del2W3)
!
! Compute gradient of GammaW (mandatory for curlBW1, curlBW2, curlBW3)
!
      if (llongitudinalW) then
        call grad(f,iGammaW,p%gGammaW1)
        call grad(f,iGammaW+1,p%gGammaW2)
        call grad(f,iGammaW+2,p%gGammaW3)
      else
        call del2v_etc(f,iWW1,GRADDIV=p%gGammaW1)
        call del2v_etc(f,iWW2,GRADDIV=p%gGammaW2)
        call del2v_etc(f,iWW3,GRADDIV=p%gGammaW3)
      endif
!
! Compute curlBW1, curlBW2, curlBW3 as - del2W + grad(divW)
! (mandatory for EOM)
!
      p%curlBW1=-p%del2W1+p%gGammaW1
      p%curlBW2=-p%del2W2+p%gGammaW2
      p%curlBW3=-p%del2W3+p%gGammaW3
!
! pencils dW1, dW2, dW3
!
      if (lpenc_requested(i_dW1)) p%dW1=f(l1:l2,m,n,idWW1:idWW1+2)
      if (lpenc_requested(i_dW1sq)) call dot2_mn(p%dW1,p%dW1sq)
      if (lpenc_requested(i_dW2)) p%dW2=f(l1:l2,m,n,idWW2:idWW2+2)
      if (lpenc_requested(i_dW2sq)) call dot2_mn(p%dW2,p%dW2sq)
      if (lpenc_requested(i_dW3)) p%dW3=f(l1:l2,m,n,idWW3:idWW3+2)
      if (lpenc_requested(i_dW3sq)) call dot2_mn(p%dW3,p%dW3sq)
!
! pencils W1, W2, W3, W1sq, W2sq, W3sq (mandatory for EOM)
!
      p%W1=f(l1:l2,m,n,iWW1:iWW1+2)
      call dot2_mn(p%W1,p%W1sq)
      p%W2=f(l1:l2,m,n,iWW2:iWW2+2)
      call dot2_mn(p%W2,p%W2sq)
      p%W3=f(l1:l2,m,n,iWW3:iWW3+2)
      call dot2_mn(p%W3,p%W3sq)
!
!  WdotW pencils (required if llongitudinalW=T)
!
      if (lpenc_requested(i_W1dotW1)) call dot(p%W1,p%dW1,p%W1dotW1)
      if (lpenc_requested(i_W1dotW2)) call dot(p%W1,p%dW2,p%W1dotW2)
      if (lpenc_requested(i_W1dotW3)) call dot(p%W1,p%dW3,p%W1dotW3)
      if (lpenc_requested(i_W2dotW1)) call dot(p%W2,p%dW1,p%W2dotW1)
      if (lpenc_requested(i_W2dotW2)) call dot(p%W2,p%dW2,p%W2dotW2)
      if (lpenc_requested(i_W2dotW3)) call dot(p%W2,p%dW3,p%W2dotW3)
      if (lpenc_requested(i_W3dotW1)) call dot(p%W3,p%dW1,p%W3dotW1)
      if (lpenc_requested(i_W3dotW2)) call dot(p%W3,p%dW2,p%W3dotW2)
      if (lpenc_requested(i_W3dotW3)) call dot(p%W3,p%dW3,p%W3dotW3)
!
! WWddot2
!
      if (lWddot_as_aux) then
        call dot2_mn(f(l1:l2,m,n,iWW1ddot:iWW1ddot+2),p%W1ddotsq)
        call dot2_mn(f(l1:l2,m,n,iWW2ddot:iWW2ddot+2),p%W2ddotsq)
        call dot2_mn(f(l1:l2,m,n,iWW3ddot:iWW3ddot+2),p%W3ddotsq)
      else
        p%W1ddotsq=0.
        p%W2ddotsq=0.
        p%W3ddotsq=0.
      endif
!
!  Pencils W1ddot_gw, W2ddot_gw, W3ddot_gw (mandatory for EOM)
!  They correspond to the terms in the EOM of dW1, dW2, dW3 that
!  are proportional to coupl_gw.
!  These terms are characteristic of Yang-Mills fields and appear
!  due to the non-Abelian nature of the gauge group SU(2).
!
      if (coupl_gw /= 0) then
        ! Add terms from coupl_gw x f^{abc} W^b GammaW^c
        ! Antonino: Tanmay's gauge fixing terms (-gw eps^abc W_c^i GammaW_b)
        call dot_mn_sv_pencil(p%W2,p%GammaW3,tmp1)
        call dot_mn_sv_pencil(p%W3,p%GammaW2,tmp2)
        p%W1ddot_gw=-coupl_gw*(tmp1-tmp2)
        call dot_mn_sv_pencil(p%W3,p%GammaW1,tmp1)
        call dot_mn_sv_pencil(p%W1,p%GammaW3,tmp2)
        p%W2ddot_gw=-coupl_gw*(tmp1-tmp2)
        call dot_mn_sv_pencil(p%W1,p%GammaW2,tmp1)
        call dot_mn_sv_pencil(p%W2,p%GammaW1,tmp2)
        p%W3ddot_gw=-coupl_gw*(tmp1-tmp2)
        ! Add terms from coupl_gw x f^{abc} d_k W^b W_k^c
        ! compute terms W3 . grad W2 and W2 . grad W3
        call u_dot_grad(f,iWW2,p%W2ij,p%W3,tmp1)
        call u_dot_grad(f,iWW3,p%W3ij,p%W2,tmp2)
        p%W1ddot_gw=p%W1ddot_gw-2.*coupl_gw*(tmp1-tmp2)
        ! compute terms W1 . grad W3 and W3 . grad W1
        call u_dot_grad(f,iWW3,p%W3ij,p%W1,tmp1)
        call u_dot_grad(f,iWW1,p%W1ij,p%W3,tmp2)
        p%W2ddot_gw=p%W2ddot_gw-2.*coupl_gw*(tmp1-tmp2)
        ! compute terms W2 . grad W1 and W1 . grad W2
        call u_dot_grad(f,iWW1,p%W1ij,p%W2,tmp1)
        call u_dot_grad(f,iWW2,p%W2ij,p%W1,tmp2)
        p%W3ddot_gw=p%W3ddot_gw-2.*coupl_gw*(tmp1-tmp2)

        ! two more additional terms come from coupl_gw x f^{abc} W_k^b W_{ik}^c
        ! the first one is coupl_gw x f^{abc} W_k^b nabla W_k^c
        do i=1,3; do j=1,3
        p%W1ddot_gw(:,j)=p%W1ddot_gw(:,j) - &
            coupl_gw*(p%W2(:,i)*p%W3ij(:,i,j) - p%W3(:,i)*p%W2ij(:,i,j))
        p%W2ddot_gw(:,j)=p%W2ddot_gw(:,j) - &
            coupl_gw*(p%W3(:,i)*p%W1ij(:,i,j) - p%W1(:,i)*p%W3ij(:,i,j))
        p%W3ddot_gw(:,j)=p%W3ddot_gw(:,j) - &
            coupl_gw*(p%W1(:,i)*p%W2ij(:,i,j) - p%W2(:,i)*p%W1ij(:,i,j))
        enddo; enddo

        ! the second one is coupl_gw^2 f^{abc} f^{cde} W_k^b W_k^e W_i^d
        ! which can be split into two parts:
        ! 1) - coupl_gw^2 (W_b . W_b) W^a
        ! 2) + coupl_gw^2 (W_a . W_b) W^b
        tmp0=p%W1sq+p%W2sq+p%W3sq
        call dot_mn_sv_pencil(p%W1,tmp0,tmp1)
        p%W1ddot_gw2=-coupl_gw**2*tmp1
        call dot_mn_sv_pencil(p%W2,tmp0,tmp1)
        p%W2ddot_gw2=-coupl_gw**2*tmp1
        call dot_mn_sv_pencil(p%W3,tmp0,tmp1)
        p%W3ddot_gw2=-coupl_gw**2*tmp1
        ! df(l1:l2,m,n,idWW1:idWW1+2)=df(l1:l2,m,n,idWW1:idWW1+2) - &
        !     coupl_gw**2*(p%W1sq+p%W2sq+p%W3sq)*f(l1:l2,m,n,iWW1:iWW1+2)
        ! df(l1:l2,m,n,idWW2:idWW2+2)=df(l1:l2,m,n,idWW2:idWW2+2) - &
        !     coupl_gw**2*(p%W1sq+p%W2sq+p%W3sq)*f(l1:l2,m,n,iWW2:iWW2+2)
        ! df(l1:l2,m,n,idWW3:idWW3+2)=df(l1:l2,m,n,idWW3:idWW3+2) - &
        !     coupl_gw**2*(p%W1sq+p%W2sq+p%W3sq)*f(l1:l2,m,n,iWW3:iWW3+2)
        do i=1,3
          p%W1ddot_gw2(:,i)=p%W1ddot_gw2(:,i) + &
              coupl_gw**2*p%W1(:,i)*p%W1sq
          p%W2ddot_gw2(:,i)=p%W2ddot_gw2(:,i) + &
              coupl_gw**2*p%W2(:,i)*p%W2sq
          p%W3ddot_gw2(:,i)=p%W3ddot_gw2(:,i) + &
              coupl_gw**2*p%W3(:,i)*p%W3sq
          ! df(l1:l2,m,n,idWW1+i)=df(l1:l2,m,n,idWW1+i) + &
          !     coupl_gw**2*p%W1(:,i+1)*p%W1sq
          ! df(l1:l2,m,n,idWW2+i)=df(l1:l2,m,n,idWW2+i) + &
          !     coupl_gw**2*p%W2(:,i+1)*p%W2sq
          ! df(l1:l2,m,n,idWW3+i)=df(l1:l2,m,n,idWW3+i) + &
          !     coupl_gw**2*p%W3(:,i+1)*p%W3sq
        enddo
        call dot_mn(p%W1,p%W2,tmp0)
        do i=1,3
          p%W1ddot_gw2(:,i)=p%W1ddot_gw2(:,i) + &
              coupl_gw**2*p%W2(:,i)*tmp0
          p%W2ddot_gw2(:,i)=p%W2ddot_gw2(:,i) + &
              coupl_gw**2*p%W1(:,i)*tmp0
          ! df(l1:l2,m,n,idWW1+i)=df(l1:l2,m,n,idWW1+i) + &
          !     coupl_gw**2*p%W2(:,i+1)*tmp
          ! df(l1:l2,m,n,idWW2+i)=df(l1:l2,m,n,idWW2+i) + &
          !     coupl_gw**2*p%W1(:,i+1)*tmp
        enddo
        call dot_mn(p%W1,p%W3,tmp0)
        do i=1,3
          p%W1ddot_gw2(:,i)=p%W1ddot_gw2(:,i) + &
              coupl_gw**2*p%W3(:,i)*tmp0
          p%W3ddot_gw2(:,i)=p%W3ddot_gw2(:,i) + &
              coupl_gw**2*p%W1(:,i)*tmp0
          ! df(l1:l2,m,n,idWW1+i)=df(l1:l2,m,n,idWW1+i) + &
          !     coupl_gw**2*p%W3(:,i+1)*tmp
          ! df(l1:l2,m,n,idWW3+i)=df(l1:l2,m,n,idWW3+i) + &
          !     coupl_gw**2*p%W3(:,i+1)*tmp
        enddo
        call dot_mn(p%W2,p%W3,tmp0)
        do i=1,3
          p%W2ddot_gw2(:,i)=p%W2ddot_gw2(:,i) + &
              coupl_gw**2*p%W3(:,i)*tmp0
          p%W3ddot_gw2(:,i)=p%W3ddot_gw2(:,i) + &
              coupl_gw**2*p%W2(:,i)*tmp0
          ! df(l1:l2,m,n,idWW2+i)=df(l1:l2,m,n,idWW2+i) + &
          !     coupl_gw**2*p%W3(:,i+1)*tmp
          ! df(l1:l2,m,n,idWW3+i)=df(l1:l2,m,n,idWW3+i) + &
          !     coupl_gw**2*p%W2(:,i+1)*tmp
        enddo
      endif
      !

      ! call u_dot_grad(f,iuu,p%uij,p%oo,p%ogu,UPWIND=lupw_uu)
      ! compute terms W3 . grad W2 and W2 . grad W3
      ! call u_dot_grad(f,iWW2,p%W2ij,p%W3,tmp1)
      ! call u_dot_grad(f,iWW3,p%W3ij,p%W2,tmp2)
      ! df(l1:l2,m,n,idWW1:idWW1+2)=df(l1:l2,m,n,idWW1:idWW1+2) - &
      !    2.*coupl_gw*(tmp1-tmp2)
      ! ! compute terms W1 . grad W3 and W3 . grad W1
      ! call u_dot_grad(f,iWW3,p%W3ij,p%W1,tmp1)
      ! call u_dot_grad(f,iWW1,p%W1ij,p%W3,tmp2)
      ! df(l1:l2,m,n,idWW2:idWW2+2)=df(l1:l2,m,n,idWW2:idWW2+2) - &
      !    2.*coupl_gw*(tmp1-tmp2)
      ! ! compute terms W2 . grad W1 and W1 . grad W2
      ! call u_dot_grad(f,iWW1,p%W1ij,p%W2,tmp1)
      ! call u_dot_grad(f,iWW2,p%W2ij,p%W1,tmp2)
      ! df(l1:l2,m,n,idWW3:idWW3+2)=df(l1:l2,m,n,idWW3:idWW3+2) - &
      !    2.*coupl_gw*(tmp1-tmp2)
      ! two more additional terms come from coupl_gw x f^{abc} W_k^b W_{ik}^c
      ! the first one is coupl_gw x f^{abc} W_k^b d_i W_k^c
      ! do i=1,3; do j=1,3
      !   df(l1:l2,m,n,idWW1+j-1)=df(l1:l2,m,n,idWW1+j-1) - &
      !     coupl_gw*(p%W2(:,i)*p%W3ij(:,i,j) - p%W3(:,i)*p%W2ij(:,i,j))
      !   df(l1:l2,m,n,idWW2+j-1)=df(l1:l2,m,n,idWW2+j-1) - &
      !     coupl_gw*(p%W3(:,i)*p%W1ij(:,i,j) - p%W1(:,i)*p%W3ij(:,i,j))
      !   df(l1:l2,m,n,idWW3+j-1)=df(l1:l2,m,n,idWW3+j-1) - &
      !     coupl_gw*(p%W1(:,i)*p%W2ij(:,i,j) - p%W2(:,i)*p%W1ij(:,i,j))
      ! enddo; enddo
!
! !  curle
! !
!       if (idiag_BcurlEm/=0) then
!         call curl(f,iex,p%curle)
!         call dot(p%bb,p%curle,p%BcurlE)
!       endif
!
! !  del2ee
! !
!       if (eta_ee/=0.) call del2v(f,iex,p%del2ee)
!
! a0 & ga0
!
      ! if (iW0>0) then
      !   p%a0=f(l1:l2,m,n,iW0)
      !   call grad(f,iW0,p%ga0)
      ! endif
!
!  divJ (using Ohm's law)
!  divJ=sigma*[divE+eps_ijk*(u_j,i * b_k + u_j * b_k,i)]
!  The use if eta may be suspect and should be checked.
! !
!       if (lpenc_requested(i_divJ)) then
!         tmp=0.
!         do i=1,3
!         do j=1,3
!         do k=1,3
!           tmp=tmp+levi_civita(i,j,k)*(p%uij(:,j,i)*p%bb(:,k)+p%uu(:,j)*p%bij(:,k,i))
!         enddo
!         enddo
!         enddo
!         if (lswitch_off_divJ) then
!           p%divJ=0.
!         else
!           if (eta==0.) then
! !
! !  The following expression ignores gradients of p%sigE
! !
!             !call fatal_error('disp_current/calc_pencils_special', "eta=0 not ok here")
!             p%divJ=(p%divE+tmp)*p%sigE
!           else
!             p%divJ=(p%divE+tmp)/(mu0*eta)
!           endif
!         endif
!       endif
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine calc_constrainteqn(p,tmp,constrainteqn)

      type(pencil_case) :: p
      real, dimension(nx,3),intent(IN) :: tmp
      real, dimension(nx,3), intent(OUT) :: constrainteqn
      real, dimension(nx,3) :: constrainteqn1

      ! constraint equation to be adapted for WW
      ! constrainteqn1=sqrt(p%divE**2+tmp**2)
      constrainteqn1(:,1)=sqrt(p%divdotW1**2 + tmp(:,1)**2)
      constrainteqn1(:,2)=sqrt(p%divdotW2**2 + tmp(:,2)**2)
      constrainteqn1(:,3)=sqrt(p%divdotW3**2 + tmp(:,3)**2)
!
!  in the following, should use "where"
!
      if (any(constrainteqn1 == 0.)) then
        constrainteqn=0.
      else
        !constrainteqn=(p%divE-tmp)/constrainteqn1
        !constrainteqn=(p%divdotW1+p%divdotW2+p%divdotW3-tmp)/constrainteqn1
        constrainteqn(:,1)=(p%divdotW1 + tmp(:,1))/constrainteqn1(:,1)
        constrainteqn(:,2)=(p%divdotW2 + tmp(:,2))/constrainteqn1(:,2)
        constrainteqn(:,3)=(p%divdotW3 + tmp(:,3))/constrainteqn1(:,3)
      endif
     endsubroutine calc_constrainteqn
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
!
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      ! real, dimension (nx,3,3) :: gtmp, dJdt, del2JJ
      real, dimension (nx,3) :: tmp=0. ! del2a0 !constrainteqn, constrainteqn1
      ! real :: inflation_factor=0. ! mfpf=0., fppf=0.
      integer :: j
!
      intent(inout) :: p
      intent(inout) :: f, df
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
!
!  Solve for Gamma (unless lswitch_off_Gamma) and possibly for charge density.
!  Add to the existing tmp and update df(l1:l2,m,n,irhoe).
!
      if (llongitudinalW) then
        ! if (lsolve_chargedensityW) then
        !   tmp(:,1)=tmp(:,1)+p%rhoW1
        !   tmp(:,2)=tmp(:,2)+p%rhoW2
        !   tmp(:,3)=tmp(:,3)+p%rhoW3
        ! endif
        if (lphi_doublet .and. lphi_weakcharge .and. coupl_gw /= 0) then
          p%rhoe_higgsW1=-coupl_gw*(p%phi*p%cov_der(:,1,4) - &
                                p%phi_doublet(:,1)*p%cov_der(:,1,3) + &
                                p%phi_doublet(:,2)*p%cov_der(:,1,2) - &
                                p%phi_doublet(:,3)*p%cov_der(:,1,1))
          p%rhoe_higgsW2=-coupl_gw*(-p%phi*p%cov_der(:,1,3) - &
                                p%phi_doublet(:,1)*p%cov_der(:,1,4) + &
                                p%phi_doublet(:,2)*p%cov_der(:,1,1) + &
                                p%phi_doublet(:,3)*p%cov_der(:,1,2))
          p%rhoe_higgsW3=-coupl_gw*(p%phi*p%cov_der(:,1,2) - &
                                p%phi_doublet(:,1)*p%cov_der(:,1,1) - &
                                p%phi_doublet(:,2)*p%cov_der(:,1,4) + &
                                p%phi_doublet(:,3)*p%cov_der(:,1,3))
          !
          ! add Higgs charge density to tmp
          !
          tmp(:,1) = tmp(:,1) - p%rhoe_higgsW1
          tmp(:,2) = tmp(:,2) - p%rhoe_higgsW2
          tmp(:,3) = tmp(:,3) - p%rhoe_higgsW3
        endif
        ! add coupl_gw x f^{abc} x W_b . dW_b contraction to Gauss constraint
        if (coupl_gw /= 0) then
          tmp(:,1) = tmp(:,1) - coupl_gw*(p%W2dotW3 - p%W3dotW2)
          tmp(:,2) = tmp(:,2) - coupl_gw*(p%W3dotW1 - p%W1dotW3)
          tmp(:,3) = tmp(:,3) - coupl_gw*(p%W1dotW2 - p%W2dotW1)
        endif
        df(l1:l2,m,n,iGammaW)=df(l1:l2,m,n,iGammaW) &
            +(1.+weight_longitudinalW)*p%divdotW1+weight_longitudinalW*tmp(:,1)
        df(l1:l2,m,n,iGammaW+1)=df(l1:l2,m,n,iGammaW+1) &
            +(1.+weight_longitudinalW)*p%divdotW2+weight_longitudinalW*tmp(:,2)
        df(l1:l2,m,n,iGammaW+2)=df(l1:l2,m,n,iGammaW+2) &
            +(1.+weight_longitudinalW)*p%divdotW3+weight_longitudinalW*tmp(:,3)
      endif
!
!  solve: dE/dt = curlB - ...
!  Calculate curlB as -del2a, because curlB leads to instability.
!  Solve dA/dt = -E.
!
      ! if (lmagnetic) then
      !   df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-p%el
      df(l1:l2,m,n,iWW1:iWW1+2)= &
            df(l1:l2,m,n,iWW1:iWW1+2)+p%dW1
      df(l1:l2,m,n,iWW2:iWW2+2)= &
            df(l1:l2,m,n,iWW2:iWW2+2)+p%dW2
      df(l1:l2,m,n,iWW3:iWW3+2)= &
            df(l1:l2,m,n,iWW3:iWW3+2)+p%dW3

      !  Maxwell equation otherwise the same in both gauges.
      df(l1:l2,m,n,idWW1:idWW1+2)=&
            df(l1:l2,m,n,idWW1:idWW1+2)-c_light2*p%curlBW1
      df(l1:l2,m,n,idWW2:idWW2+2)=&
            df(l1:l2,m,n,idWW2:idWW2+2)-c_light2*p%curlBW2
      df(l1:l2,m,n,idWW3:idWW3+2)=&
            df(l1:l2,m,n,idWW3:idWW3+2)-c_light2*p%curlBW3
!
!  Add Yang-Mills terms proportional to coupl_gw and coupl_gw^2
!
      if (coupl_gw /= 0) then
        df(l1:l2,m,n,idWW1:idWW1+2)=df(l1:l2,m,n,idWW1:idWW1+2) + &
              p%W1ddot_gw + p%W1ddot_gw2
        df(l1:l2,m,n,idWW2:idWW2+2)=df(l1:l2,m,n,idWW2:idWW2+2) + &
              p%W2ddot_gw + p%W2ddot_gw2
        df(l1:l2,m,n,idWW3:idWW3+2)=df(l1:l2,m,n,idWW3:idWW3+2) + &
              p%W3ddot_gw + p%W3ddot_gw2
      endif
!
      !   df(l1:l2,m,n,iex:iez)=df(l1:l2,m,n,iex:iez)+c_light2*(p%curlb-mu0*p%jj_ohm)
!
!  Solve for charge density
!
      ! to be updated with Higgs charges
      ! if (lsolve_chargedensityW) then
      !   df(l1:l2,m,n,irhow)=df(l1:l2,m,n,irhow)-p%divJ
      ! endif
!
! Add diffusion term if eta_WW /= 0
      if (eta_WW/=0.) then
        df(l1:l2,m,n,iWW1:iWW1+2)=df(l1:l2,m,n,iWW1:iWW1+2)+c_light2*eta_WW*p%del2W1
        df(l1:l2,m,n,iWW2:iWW2+2)=df(l1:l2,m,n,iWW2:iWW2+2)+c_light2*eta_WW*p%del2W2
        df(l1:l2,m,n,iWW3:iWW3+2)=df(l1:l2,m,n,iWW3:iWW3+2)+c_light2*eta_WW*p%del2W3
      endif
!
!  If Higgs doublet, add current from Higgs SU(2) weak charge,
!  compute pencils p%jj_higgsW1, p%jj_higgsW2, p%jj_higgsW3
!
      if (lphi_doublet .and. lphi_weakcharge .and. coupl_gw /= 0) then
        ! currents from charged Higgs doublet
        do j=1,3
          p%jj_higgsW1(:,j)=coupl_gw*(p%phi*p%cov_der(:,j+1,4) - &
                                p%phi_doublet(:,1)*p%cov_der(:,j+1,3) + &
                                p%phi_doublet(:,2)*p%cov_der(:,j+1,2) - &
                                p%phi_doublet(:,3)*p%cov_der(:,j+1,1))
          p%jj_higgsW2(:,j)=coupl_gw*(-p%phi*p%cov_der(:,j+1,3) - &
                                p%phi_doublet(:,1)*p%cov_der(:,j+1,4) + &
                                p%phi_doublet(:,2)*p%cov_der(:,j+1,1) + &
                                p%phi_doublet(:,3)*p%cov_der(:,j+1,2))
          p%jj_higgsW3(:,j)=coupl_gw*(p%phi*p%cov_der(:,j+1,2) - &
                                p%phi_doublet(:,1)*p%cov_der(:,j+1,1) - &
                                p%phi_doublet(:,2)*p%cov_der(:,j+1,4) + &
                                p%phi_doublet(:,3)*p%cov_der(:,j+1,3))
        enddo
        df(l1:l2,m,n,idWW1:idWW1+2)=df(l1:l2,m,n,idWW1:idWW1+2) + &
            p%jj_higgsW1
        df(l1:l2,m,n,idWW2:idWW2+2)=df(l1:l2,m,n,idWW2:idWW2+2) + &
            p%jj_higgsW2
        df(l1:l2,m,n,idWW3:idWW3+2)=df(l1:l2,m,n,idWW3:idWW3+2) + &
            p%jj_higgsW3
      endif
!
!  Compute WWddot_as_aux
!
      if (lWddot_as_aux) then
        f(l1:l2,m,n,iWW1ddot:iWW1ddot+2) = - &
            c_light2*p%curlBW1
        f(l1:l2,m,n,iWW2ddot:iWW2ddot+2) = - &
            c_light2*p%curlBW2
        f(l1:l2,m,n,iWW3ddot:iWW3ddot+2) = - &
            c_light2*p%curlBW3
        if (coupl_gw /= 0) then
          f(l1:l2,m,n,iWW1ddot:iWW1ddot+2)=f(l1:l2,m,n,iWW1ddot:iWW1ddot+2) + &
                p%W1ddot_gw + p%W1ddot_gw2
          f(l1:l2,m,n,iWW2ddot:iWW2ddot+2)=f(l1:l2,m,n,iWW2ddot:iWW2ddot+2) + &
                p%W2ddot_gw + p%W2ddot_gw2
          f(l1:l2,m,n,iWW3ddot:iWW3ddot+2)=f(l1:l2,m,n,iWW3ddot:iWW3ddot+2) + &
                p%W3ddot_gw + p%W3ddot_gw2
          if (lphi_doublet .and. lphi_weakcharge) then
            f(l1:l2,m,n,iWW1ddot:iWW1ddot+2)=f(l1:l2,m,n,iWW1ddot:iWW1ddot+2) - &
                p%jj_higgsW1
            f(l1:l2,m,n,iWW2ddot:iWW2ddot+2)=f(l1:l2,m,n,iWW2ddot:iWW2ddot+2) - &
                p%jj_higgsW2
            f(l1:l2,m,n,iWW3ddot:iWW3ddot+2)=f(l1:l2,m,n,iWW3ddot:iWW3ddot+2) - &
                p%jj_higgsW3
          endif
        endif
      endif
!
!  If requested, put divE
!
      if (ldivdotW_as_aux) then
        f(l1:l2,m,n,idivdotW)=p%divdotW1
        f(l1:l2,m,n,idivdotW+1)=p%divdotW2
        f(l1:l2,m,n,idivdotW+2)=p%divdotW3
      endif
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
      real, dimension(nx,3) :: tmp, constrainteqn
      ! real :: mfpf=0.,fppf=0.
      ! real, dimension(nx,3) :: gtmp
      integer :: i,j

      ! call save_name(get_mfpf(),idiag_mfpf)
      ! call save_name(get_fppf(),idiag_fppf)
      ! call save_name(scl_factor_target,idiag_afact)
      ! if (idiag_EEEM/=0) call sum_mn_name(.5*(p%e2+p%b2),idiag_EEEM)
      if (idiag_W1dotW1m /= 0) call sum_mn_name(p%W1dotW1,idiag_W1dotW1m)
      if (idiag_W2dotW2m /= 0) call sum_mn_name(p%W2dotW2,idiag_W2dotW2m)
      if (idiag_W3dotW3m /= 0) call sum_mn_name(p%W3dotW3,idiag_W3dotW3m)
      if (idiag_W1rms/=0) call sum_mn_name(p%W1**2,idiag_W1rms,lsqrt=.true.)
      if (idiag_W2rms/=0) call sum_mn_name(p%W2**2,idiag_W2rms,lsqrt=.true.)
      if (idiag_W3rms/=0) call sum_mn_name(p%W3**2,idiag_W3rms,lsqrt=.true.)
      if (idiag_dW1rms/=0) call sum_mn_name(p%dW1**2,idiag_dW1rms,lsqrt=.true.)
      if (idiag_dW2rms/=0) call sum_mn_name(p%dW2**2,idiag_dW2rms,lsqrt=.true.)
      if (idiag_dW3rms/=0) call sum_mn_name(p%dW3**2,idiag_dW3rms,lsqrt=.true.)
      call sum_mn_name(p%dW1(:,1),idiag_dW1xm)
      call sum_mn_name(p%dW1(:,2),idiag_dW1ym)
      call sum_mn_name(p%dW1(:,3),idiag_dW1zm)
      call sum_mn_name(p%dW2(:,1),idiag_dW2xm)
      call sum_mn_name(p%dW2(:,2),idiag_dW2ym)
      call sum_mn_name(p%dW2(:,3),idiag_dW2zm)
      call sum_mn_name(p%dW3(:,1),idiag_dW3xm)
      call sum_mn_name(p%dW3(:,2),idiag_dW3ym)
      call sum_mn_name(p%dW3(:,3),idiag_dW3zm)
      call sum_mn_name(p%W1(:,1),idiag_W1xm)
      call sum_mn_name(p%W1(:,2),idiag_W1ym)
      call sum_mn_name(p%W1(:,3),idiag_W1zm)
      call sum_mn_name(p%W2(:,1),idiag_W2xm)
      call sum_mn_name(p%W2(:,2),idiag_W2ym)
      call sum_mn_name(p%W2(:,3),idiag_W2zm)
      call sum_mn_name(p%W3(:,1),idiag_W3xm)
      call sum_mn_name(p%W3(:,2),idiag_W3ym)
      call sum_mn_name(p%W3(:,3),idiag_W3zm)
      ! call sum_mn_name(p%sigE,idiag_sigEm)
      ! call sum_mn_name(p%sigB,idiag_sigBm)
      call max_mn_name(p%W1sq,idiag_W1max,lsqrt=.true.)
      call max_mn_name(p%W2sq,idiag_W2max,lsqrt=.true.)
      call max_mn_name(p%W3sq,idiag_W3max,lsqrt=.true.)
      call max_mn_name(p%dW1sq,idiag_dW1max,lsqrt=.true.)
      call max_mn_name(p%dW2sq,idiag_dW2max,lsqrt=.true.)
      call max_mn_name(p%dW3sq,idiag_dW3max,lsqrt=.true.)

      if (idiag_W1ddotrms /= 0) &
          call sum_mn_name(p%W1ddotsq,idiag_W1ddotrms,lsqrt=.true.)
      if (idiag_W2ddotrms /= 0) &
          call sum_mn_name(p%W2ddotsq,idiag_W2ddotrms,lsqrt=.true.)
      if (idiag_W3ddotrms /= 0) &
          call sum_mn_name(p%W3ddotsq,idiag_W3ddotrms,lsqrt=.true.)

      if (idiag_divW1rms /= 0) call sum_mn_name(p%GammaW1**2,idiag_divW1rms,lsqrt=.true.)
      if (idiag_divW2rms /= 0) call sum_mn_name(p%GammaW2**2,idiag_divW2rms,lsqrt=.true.)
      if (idiag_divW3rms /= 0) call sum_mn_name(p%GammaW3**2,idiag_divW3rms,lsqrt=.true.)
      if (idiag_divW1m /= 0) call sum_mn_name(p%GammaW1,idiag_divW1m)
      if (idiag_divW2m /= 0) call sum_mn_name(p%GammaW2,idiag_divW2m)
      if (idiag_divW3m /= 0) call sum_mn_name(p%GammaW3,idiag_divW3m)

      if (idiag_divdotW1rms /= 0) &
          call sum_mn_name(p%divdotW1**2,idiag_divdotW1rms,lsqrt=.true.)
      if (idiag_divdotW2rms /= 0) &
          call sum_mn_name(p%divdotW2**2,idiag_divdotW2rms,lsqrt=.true.)
      if (idiag_divdotW3rms /= 0) &
          call sum_mn_name(p%divdotW3**2,idiag_divdotW3rms,lsqrt=.true.)

      if (idiag_divdotW1m /= 0) call sum_mn_name(p%divdotW1,idiag_divdotW1m)
      if (idiag_divdotW2m /= 0) call sum_mn_name(p%divdotW2,idiag_divdotW2m)
      if (idiag_divdotW3m /= 0) call sum_mn_name(p%divdotW3,idiag_divdotW3m)

      call sum_mn_name(p%rhoW1,idiag_rhoW1m)
      call sum_mn_name(p%rhoW2,idiag_rhoW2m)
      call sum_mn_name(p%rhoW3,idiag_rhoW3m)

      call sum_mn_name(p%rhoW1**2,idiag_rhoW1rms,lsqrt=.true.)
      call sum_mn_name(p%rhoW2**2,idiag_rhoW2rms,lsqrt=.true.)
      call sum_mn_name(p%rhoW3**2,idiag_rhoW3rms,lsqrt=.true.)

  !   endif
      if(idiag_constrainteqnW > 0) then
        call calc_constrainteqn(p,tmp,constrainteqn)
        call sum_mn_name(constrainteqn,idiag_constrainteqnW)
      endif
    endsubroutine calc_diagnostics_special
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
      integer :: iname,inamez
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
            idiag_W1rms=0; idiag_W2rms=0; idiag_W3rms=0;
            idiag_W1max=0; idiag_W2max=0; idiag_W3max=0;
            idiag_dW1rms=0; idiag_dW2rms=0; idiag_dW3rms=0;
            idiag_dW1max=0; idiag_dW2max=0; idiag_dW3max=0;
            idiag_W1ddotrms=0; idiag_W2ddotrms=0; idiag_W3ddotrms=0;
            idiag_divW1rms=0; idiag_divW2rms=0; idiag_divW3rms=0;
            idiag_divW1m=0; idiag_divW2m=0; idiag_divW3m=0;
            idiag_divdotW1rms=0; idiag_divdotW2rms=0; idiag_divdotW3rms=0;
            idiag_divdotW1m=0; idiag_divdotW2m=0; idiag_divdotW3m=0;
            idiag_rhoW1rms=0; idiag_rhoW2rms=0; idiag_rhoW3rms=0;
            idiag_rhoW1m=0; idiag_rhoW2m=0; idiag_rhoW3m=0;
            idiag_constrainteqnW=0;
            idiag_W1dotW1m=0; idiag_W2dotW2m=0; idiag_W3dotW3m=0;
            idiag_W1xm=0; idiag_W1ym=0; idiag_W1zm=0;
            idiag_W2xm=0; idiag_W2ym=0; idiag_W2zm=0;
            idiag_W3xm=0; idiag_W3ym=0; idiag_W3zm=0;
            idiag_dW1xm=0; idiag_dW1ym=0; idiag_dW1zm=0;
            idiag_dW2xm=0; idiag_dW2ym=0; idiag_dW2zm=0;
            idiag_dW3xm=0; idiag_dW3ym=0; idiag_dW3zm=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'W1rms',idiag_W1rms)
        call parse_name(iname,cname(iname),cform(iname),'W2rms',idiag_W2rms)
        call parse_name(iname,cname(iname),cform(iname),'W3rms',idiag_W3rms)
        call parse_name(iname,cname(iname),cform(iname),'W1max',idiag_W1max)
        call parse_name(iname,cname(iname),cform(iname),'W2max',idiag_W2max)
        call parse_name(iname,cname(iname),cform(iname),'W3max',idiag_W3max)
        call parse_name(iname,cname(iname),cform(iname),'dW1rms',idiag_dW1rms)
        call parse_name(iname,cname(iname),cform(iname),'dW2rms',idiag_dW2rms)
        call parse_name(iname,cname(iname),cform(iname),'dW3rms',idiag_dW3rms)
        call parse_name(iname,cname(iname),cform(iname),'dW1max',idiag_dW1max)
        call parse_name(iname,cname(iname),cform(iname),'dW2max',idiag_dW2max)
        call parse_name(iname,cname(iname),cform(iname),'dW3max',idiag_dW3max)
        call parse_name(iname,cname(iname),cform(iname),'W1ddotrms',idiag_W1ddotrms)
        call parse_name(iname,cname(iname),cform(iname),'W2ddotrms',idiag_W2ddotrms)
        call parse_name(iname,cname(iname),cform(iname),'W3ddotrms',idiag_W3ddotrms)
      enddo
!
      ! do inamez=1,nnamez
      !   call parse_name(inamez,cnamez(inamez),cformz(inamez),'exmz',idiag_exmz)
      !   call parse_name(inamez,cnamez(inamez),cformz(inamez),'eymz',idiag_eymz)
      !   call parse_name(inamez,cnamez(inamez),cformz(inamez),'ezmz',idiag_ezmz)
      ! enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='WW1') cformv='DEFINED'
        where(cnamev=='WW2') cformv='DEFINED'
        where(cnamev=='WW3') cformv='DEFINED'
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
      case ('WW1');   call assign_slices_vec (slices,f,iWW1)
      case ('WW2');   call assign_slices_vec (slices,f,iWW2)
      case ('WW3');   call assign_slices_vec (slices,f,iWW3)
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

    integer, parameter :: n_pars=40
    integer(KIND=ikind8), dimension(n_pars) :: p_par

!     call copy_addr(alpf,p_par(1))
!     call copy_addr(eta_ee,p_par(2))
!     call copy_addr(sige_prefactor,p_par(3))
!     call copy_addr(sigb_prefactor,p_par(4))
!     call copy_addr(mass_chi,p_par(5))
!     call copy_addr(igamma,p_par(6)) ! int
!     call copy_addr(ia0,p_par(7)) ! int
!     call copy_addr(idiva_name,p_par(8)) ! int
!     call copy_addr(llongitudinale,p_par(9)) ! bool
!     call copy_addr(llorenz_gauge_disp,p_par(10)) ! bool
!     call copy_addr(lphi_hom,p_par(11)) ! bool
!     call copy_addr(lnoncollinear_eb,p_par(12)) ! bool
!     call copy_addr(lnoncollinear_eb_aver,p_par(13)) ! bool
!     call copy_addr(lcollinear_eb,p_par(14)) ! bool
!     call copy_addr(lcollinear_eb_aver,p_par(15)) ! bool
!     call copy_addr(leedot_as_aux,p_par(16)) ! bool
!     call copy_addr(lcurlya,p_par(17)) ! bool
!     call copy_addr(lsolve_chargedensity,p_par(18)) ! bool
!     call copy_addr(ldive_as_aux,p_par(19)) ! bool
!     call copy_addr(lsige_as_aux,p_par(20)) ! bool
!     call copy_addr(lsigb_as_aux,p_par(21)) ! bool
!     call copy_addr(lallow_bprime_zero,p_par(22)) ! bool
!     call copy_addr(lswitch_off_divj,p_par(23)) ! bool
!     call copy_addr(lswitch_off_gamma,p_par(24)) ! bool
!     call copy_addr(lmass_suppression,p_par(25)) ! bool
!     call copy_addr(beta_inflation,p_par(26))
!     call copy_addr(c_light2,p_par(27))
!     call copy_addr(idiag_bcurlem,p_par(28)) ! int
!     call copy_addr(idiag_adphibm,p_par(29)) ! int
!     call copy_addr(idiag_johmrms,p_par(30)) ! int
!     call copy_addr(lapply_gamma_corr,p_par(31)) ! bool
!     call copy_addr(lphi_linear_regime,p_par(32)) ! bool
!     call copy_addr(weight_longitudinale,p_par(33))


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
