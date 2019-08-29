! $Id$
!
!  This module directly evolves the magnetic field instead of
!  the vector potential; special care needs to be taken to
!  guarantee the divergence of the field remains zero.
!
!  REFERENCE:
!    Maron, J. L., Mac Low, M.-M., & Oishi, J. S. 2008, ApJ, 677, 520
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagnetic = .true.
! CPARAM logical, parameter :: lbfield = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 6
! COMMUNICATED AUXILIARIES 6
!
! PENCILS PROVIDED bb(3); bbb(3); b2
! PENCILS PROVIDED bij(3,3); jj(3); j2; divb
! PENCILS PROVIDED curle(3); jxbr(3)
! PENCILS PROVIDED beta; va2
! PENCILS PROVIDED bunit(3)
!
! PENCILS PROVIDED aa(3); ss12
!
!***************************************************************
module Magnetic
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include 'magnetic.h'
  integer :: pushpars2c, pushdiags2c  ! should be procedure pointer (F2003)
!
!  Initialization parameters
!
  real, dimension(3) :: b_ext = 0.0, b0_ext = 0.0
  real :: t_bext = 0.0, t0_bext = 0.0
!
  namelist /magnetic_init_pars/ b_ext, b0_ext, t_bext, t0_bext
!
!  Runtime parameters
!
  real, dimension(2) :: eta_zdep_coeff = 0.0
  character(len=8) :: eta_zdep_prof = ''
  logical :: lbext = .false.
  logical :: limplicit_resistivity = .false.
  logical :: lcovariant_magnetic= .false.
  logical :: lresis_const = .false.
  logical :: lresis_zdep = .false.
  logical :: lresis_shock = .false.
  logical :: lresis_hyper3_mesh = .false.
  logical :: lohmic_heat = .true.
  real :: eta = 0.0
  real :: eta_shock = 0.0
  real :: eta_hyper3_mesh = 0.0
!
  namelist /magnetic_run_pars/ &
    b_ext, b0_ext, t_bext, t0_bext, &
    eta, eta_zdep_prof, eta_zdep_coeff, eta_shock, eta_hyper3_mesh, &
    limplicit_resistivity, lohmic_heat
!
!  Diagnostic variables
!
  integer :: idiag_bmax = 0     ! DIAG_DOC: $\max B$
  integer :: idiag_bmin = 0     ! DIAG_DOC: $\min B$
  integer :: idiag_brms = 0     ! DIAG_DOC: $\langle B^2\rangle^{1/2}$
  integer :: idiag_bm = 0       ! DIAG_DOC: $\langle B\rangle$
  integer :: idiag_b2m = 0      ! DIAG_DOC: $\langle B^2\rangle$
  integer :: idiag_bxmax = 0    ! DIAG_DOC: $\max|B_x|$
  integer :: idiag_bymax = 0    ! DIAG_DOC: $\max|B_y|$
  integer :: idiag_bzmax = 0    ! DIAG_DOC: $\max|B_z|$
  integer :: idiag_bxm = 0      ! DIAG_DOC: $\langle B_x\rangle$
  integer :: idiag_bym = 0      ! DIAG_DOC: $\langle B_y\rangle$
  integer :: idiag_bzm = 0      ! DIAG_DOC: $\langle B_z\rangle$
  integer :: idiag_bx2m = 0     ! DIAG_DOC: $\langle B_x^2\rangle$
  integer :: idiag_by2m = 0     ! DIAG_DOC: $\langle B_y^2\rangle$
  integer :: idiag_bz2m = 0     ! DIAG_DOC: $\langle B_z^2\rangle$
  integer :: idiag_bxbym = 0    ! DIAG_DOC: $\langle B_x B_y\rangle$
  integer :: idiag_bxbzm = 0    ! DIAG_DOC: $\langle B_x B_z\rangle$
  integer :: idiag_bybzm = 0    ! DIAG_DOC: $\langle B_y B_z\rangle$
  integer :: idiag_dbxmax = 0   ! DIAG_DOC: $\max|B_x - B_{\mathrm{ext,}x}|$
  integer :: idiag_dbymax = 0   ! DIAG_DOC: $\max|B_y - B_{\mathrm{ext,}y}|$
  integer :: idiag_dbzmax = 0   ! DIAG_DOC: $\max|B_z - B_{\mathrm{ext,}z}|$
  integer :: idiag_dbxm = 0     ! DIAG_DOC: $\langle B_x - B_{\mathrm{ext,}x}\rangle$
  integer :: idiag_dbym = 0     ! DIAG_DOC: $\langle B_y - B_{\mathrm{ext,}y}\rangle$
  integer :: idiag_dbzm = 0     ! DIAG_DOC: $\langle B_z - B_{\mathrm{ext,}z}\rangle$
  integer :: idiag_dbx2m = 0    ! DIAG_DOC: $\langle\left(B_x - B_{\mathrm{ext,}x}\right)^2\rangle$
  integer :: idiag_dby2m = 0    ! DIAG_DOC: $\langle\left(B_y - B_{\mathrm{ext,}y}\right)^2\rangle$
  integer :: idiag_dbz2m = 0    ! DIAG_DOC: $\langle\left(B_z - B_{\mathrm{ext,}z}\right)^2\rangle$
  integer :: idiag_jmax = 0     ! DIAG_DOC: $\max J$
  integer :: idiag_jmin = 0     ! DIAG_DOC: $\min J$
  integer :: idiag_jrms = 0     ! DIAG_DOC: $\langle J^2\rangle^{1/2}$
  integer :: idiag_jm = 0       ! DIAG_DOC: $\langle J\rangle$
  integer :: idiag_j2m = 0      ! DIAG_DOC: $\langle J^2\rangle$
  integer :: idiag_jxmax = 0    ! DIAG_DOC: $\max|J_x|$
  integer :: idiag_jymax = 0    ! DIAG_DOC: $\max|J_y|$
  integer :: idiag_jzmax = 0    ! DIAG_DOC: $\max|J_z|$
  integer :: idiag_jxm = 0      ! DIAG_DOC: $\langle J_x\rangle$
  integer :: idiag_jym = 0      ! DIAG_DOC: $\langle J_y\rangle$
  integer :: idiag_jzm = 0      ! DIAG_DOC: $\langle J_z\rangle$
  integer :: idiag_jx2m = 0     ! DIAG_DOC: $\langle J_x^2\rangle$
  integer :: idiag_jy2m = 0     ! DIAG_DOC: $\langle J_y^2\rangle$
  integer :: idiag_jz2m = 0     ! DIAG_DOC: $\langle J_z^2\rangle$
  integer :: idiag_divbmax = 0  ! DIAG_DOC: $\max|\nabla\cdot\vec{B}|$
  integer :: idiag_divbrms = 0  ! DIAG_DOC: $\langle\left(\nabla\cdot\vec{B}\right)^2\rangle^{1/2}$
  integer :: idiag_betamax = 0  ! DIAG_DOC: $\max\beta$
  integer :: idiag_betamin = 0  ! DIAG_DOC: $\min\beta$
  integer :: idiag_betam = 0    ! DIAG_DOC: $\langle\beta\rangle$
  integer :: idiag_vAmax = 0    ! DIAG_DOC: $\max v_A$
  integer :: idiag_vAmin = 0    ! DIAG_DOC: $\min v_A$
  integer :: idiag_vAm = 0      ! DIAG_DOC: $\langle v_A\rangle$
!
!  xy averages
!
  integer :: idiag_bmz = 0      ! XYAVG_DOC: $\langle B\rangle_{xy}$
  integer :: idiag_b2mz = 0     ! XYAVG_DOC: $\langle B^2\rangle_{xy}$
  integer :: idiag_bxmz = 0     ! XYAVG_DOC: $\langle B_x\rangle_{xy}$
  integer :: idiag_bymz = 0     ! XYAVG_DOC: $\langle B_y\rangle_{xy}$
  integer :: idiag_bzmz = 0     ! XYAVG_DOC: $\langle B_z\rangle_{xy}$
  integer :: idiag_bx2mz = 0    ! XYAVG_DOC: $\langle B_x^2\rangle_{xy}$
  integer :: idiag_by2mz = 0    ! XYAVG_DOC: $\langle B_y^2\rangle_{xy}$
  integer :: idiag_bz2mz = 0    ! XYAVG_DOC: $\langle B_z^2\rangle_{xy}$
  integer :: idiag_bxbymz = 0   ! XYAVG_DOC: $\langle B_x B_y\rangle_{xy}$
  integer :: idiag_bxbzmz = 0   ! XYAVG_DOC: $\langle B_x B_z\rangle_{xy}$
  integer :: idiag_bybzmz = 0   ! XYAVG_DOC: $\langle B_y B_z\rangle_{xy}$
  integer :: idiag_betamz = 0   ! XYAVG_DOC: $\langle\beta\rangle_{xy}$
  integer :: idiag_beta2mz = 0  ! XYAVG_DOC: $\langle\beta^2\rangle_{xy}$
!
!  yz averages
!
  integer :: idiag_bmx = 0      ! YZAVG_DOC: $\langle B\rangle_{yz}$
  integer :: idiag_b2mx = 0     ! YZAVG_DOC: $\langle B^2\rangle_{yz}$
  integer :: idiag_bxmx = 0     ! YZAVG_DOC: $\langle B_x\rangle_{yz}$
  integer :: idiag_bymx = 0     ! YZAVG_DOC: $\langle B_y\rangle_{yz}$
  integer :: idiag_bzmx = 0     ! YZAVG_DOC: $\langle B_z\rangle_{yz}$
  integer :: idiag_bx2mx = 0    ! YZAVG_DOC: $\langle B_x^2\rangle_{yz}$
  integer :: idiag_by2mx = 0    ! YZAVG_DOC: $\langle B_y^2\rangle_{yz}$
  integer :: idiag_bz2mx = 0    ! YZAVG_DOC: $\langle B_z^2\rangle_{yz}$
  integer :: idiag_bxbymx = 0   ! YZAVG_DOC: $\langle B_x B_y\rangle_{yz}$
  integer :: idiag_bxbzmx = 0   ! YZAVG_DOC: $\langle B_x B_z\rangle_{yz}$
  integer :: idiag_bybzmx = 0   ! YZAVG_DOC: $\langle B_y B_z\rangle_{yz}$
  integer :: idiag_betamx = 0   ! YZAVG_DOC: $\langle\beta\rangle_{yz}$
  integer :: idiag_beta2mx = 0  ! YZAVG_DOC: $\langle\beta^2\rangle_{yz}$
!
!  Module variables
!
  real, dimension(nx) :: maxdiffus_eta = 0.0
  real, dimension(nx) :: maxdiffus_eta3 = 0.0
  real, dimension(mx) :: uy0 = 0.0
  real, dimension(mz) :: eta_zdep, detadz
  logical :: lresistivity = .false.
  logical :: lexplicit_resistivity = .false.
  logical :: lhyper3x_mesh = .false.
  real :: diff_hyper3x_mesh = 0.0
!
!  Dummy but public variables (unfortunately)
!
  real, dimension(mz,3), parameter :: aamz = 0.0
  real, dimension(nz,3), parameter :: bbmz = 0.0, jjmz = 0.0
  logical, dimension(7) :: lresi_dep = .false.
  real, dimension(3) :: b_ext_inv = 0.0
  logical, parameter :: lcalc_aamean = .false.
  logical, parameter :: lcalc_aameanz = .false.
  logical, parameter :: lelectron_inertia = .false.
  integer, parameter :: idiag_bcosphz = 0, idiag_bsinphz = 0
  integer, parameter :: idiag_axmz = 0, idiag_aymz = 0
  real, parameter :: inertial_length = 0.0, linertial_2 = 0.0
!
  contains
!***********************************************************************
!***********************************************************************
!
!  PUBLIC AND ACTIVE ROUTINES GO BELOW HERE.
!
!***********************************************************************
!***********************************************************************
    subroutine register_magnetic()
!
!  Register the variables evolved by this magnetic module.
!
!  25-oct-13/ccyang: coded.
!
      use FArrayManager, only: farray_register_pde, farray_register_auxiliary
!
      integer :: istat
!
!  Identify version number.
!
      if (lroot) call svn_id("$Id$")
!
!  Request variable for the magnetic field.
!
      call farray_register_pde('bb', ibb, vector=3, ierr=istat)
      if (istat /= 0) call fatal_error('register_magnetic', 'cannot register the variable bb. ')
      ibx = ibb
      iby = ibx + 1
      ibz = iby + 1
!
!  Request auxiliary variable for the effective electric field.
!
      call farray_register_auxiliary('ee', iee, vector=3, communicated=.true., ierr=istat)
      if (istat /= 0) call fatal_error('register_magnetic', 'cannot register the variable ee. ')
      ieex = iee
      ieey = ieex + 1
      ieez = ieey + 1
!
!  Request auxiliary variable for the current density.
!
      call farray_register_auxiliary('jj', ijj, vector=3, communicated=.true., ierr=istat)
      if (istat /= 0) call fatal_error('register_magnetic', 'cannot register the variable jj. ')
      ijx = ijj
      ijy = ijx + 1
      ijz = ijy + 1
!
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic(f)
!
!  Conducts post-parameter-read initialization for Magnetic.
!
!  29-aug-13/ccyang: coded.
!
      use SharedVariables, only: put_shared_variable
      use Shear, only: get_hyper3x_mesh
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: ierr
!
      call keep_compiler_quiet(f)
!
!  Check the existence of external field.
!
      lbext = any(b_ext /= 0.0) .or. any(b0_ext /= 0.0)
!
!  Share it.
!
      call put_shared_variable('B_ext', B_ext, ierr)
      if (ierr /= 0) call fatal_error('initialize_magnetic', 'unable to share variable B_ext')
!
!  Calculates variables required by Magnetic and possibly other module(s).
!
      mu01 = 1.0 / mu0
!
!  B_ext_inv: currently required by test_methods, unfortunately a public variable.
!
      if (lbext) b_ext_inv = b_ext / sum(b_ext**2)
!
!  Check the switches for resistivities.
!
      if (eta /= 0.0) lresis_const = .true.
      zdep: if (eta_zdep_prof /= '') then
        call get_eta_zdep(z, eta_zdep, detadz)
        if (all(eta_zdep == 0.0)) call fatal_error('initialize_magnetic', 'unknown eta_zdep_prof ' // eta_zdep_prof)
        lresis_zdep = .true.
      endif zdep
      if (eta_shock /= 0.0) lresis_shock = .true.
      if (eta_hyper3_mesh /= 0.0) lresis_hyper3_mesh = .true.
      lresistivity = lresis_const .or. lresis_zdep .or. lresis_shock
!
!  Get hyper3x_mesh if shear is present.
!
      if (lshear) call get_hyper3x_mesh(lhyper3x_mesh, diff_hyper3x_mesh)
!
!  Sanity check
!
      if (lresis_shock .and. .not. lshock) &
        call fatal_error('initialize_magnetic', 'Shock module is required for shock resistivity. ')
!
!  Determine if any resistivity by explicit solver is present.
!
      lexplicit_resistivity = lresis_shock .or. (.not. limplicit_resistivity .and. (lresis_const .or. lresis_zdep))
!
!  Information output
!
      resis: if (lroot) then
        if (lresis_const) print *, 'initialize_magnetic: constant resistivity, eta = ', eta
        if (lresis_zdep) print *, 'initialize_magnetic: z-dependent resistivity profile ', eta_zdep_prof
        if (lresis_shock) print *, 'initialize_magnetic: shock resistivity, eta_shock = ', eta_shock
        if (lresis_hyper3_mesh) print *, 'initialize_magnetic: mesh hyper-resistivity, eta_hyper3_mesh = ', eta_hyper3_mesh
      endif resis
!
    endsubroutine initialize_magnetic
!***********************************************************************
    subroutine init_aa(f)
!
!  Sets the initial conditions relevant to Magnetic.
!
!  25-jun-13/ccyang: coded
!
      use InitialCondition, only: initial_condition_aa
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      if (.not. linitial_condition) &
        call fatal_error('init_aa', 'Please use the InitialCondition module to set up your initial conditions. ')
!
      call initial_condition_aa(f)
!
    endsubroutine init_aa
!***********************************************************************
    subroutine pencil_criteria_magnetic()
!
!  Specifies all pencils that Magnetic requires.
!
!  26-jun-13/ccyang: coded.
!
!  PDE related
!
      lpenc_requested(i_curle) = .true.
      lpenc_requested(i_bunit) = .true.
!
      ohmic: if (lenergy .and. lresistivity .and. lohmic_heat) then
        lpenc_requested(i_j2) = .true.
        entropy: if (lentropy) then
          lnTT: if (pretend_lnTT) then
            lpenc_requested(i_cv1) = .true.
            lpenc_requested(i_rho1) = .true.
            lpenc_requested(i_TT1) = .true.
          else lnTT
            lpenc_requested(i_rho1) = .true.
            lpenc_requested(i_TT1) = .true.
          endif lnTT
        else if (ltemperature) then entropy
          nolog: if (ltemperature_nolog) then
            lpenc_requested(i_cv1) = .true.
            lpenc_requested(i_rho1) = .true.
          else nolog
            lpenc_requested(i_cv1) = .true.
            lpenc_requested(i_rho1) = .true.
            lpenc_requested(i_TT1) = .true.
          endif nolog
        endif entropy
      endif ohmic
!
      if (lhydro) lpenc_requested(i_jxbr) = .true.
!
      time_step: if (ldt) then
        lpenc_requested(i_bb) = .true.
        lpenc_requested(i_rho1) = .true.
      endif time_step
!
!  Diagnostic related
!
      if (idiag_bmax /= 0 .or. idiag_bmin /= 0 .or. idiag_brms /= 0 .or. &
          idiag_bm /= 0 .or. idiag_b2m /= 0) lpenc_diagnos(i_b2) = .true.
      if (idiag_bxmax /= 0 .or. idiag_bymax /= 0 .or. idiag_bzmax /= 0 .or. &
          idiag_bxm /= 0 .or. idiag_bym /= 0 .or. idiag_bzm /= 0 .or. &
          idiag_bx2m /= 0 .or. idiag_by2m /= 0 .or. idiag_bz2m /= 0 .or. &
          idiag_bxbym /= 0 .or. idiag_bxbzm /= 0 .or. idiag_bybzm /= 0) lpenc_diagnos(i_bb) = .true.
      if (idiag_dbxmax /= 0 .or. idiag_dbymax /= 0 .or. idiag_dbzmax /= 0 .or. &
          idiag_dbxm /= 0 .or. idiag_dbym /= 0 .or. idiag_dbzm /= 0 .or. &
          idiag_dbx2m /= 0 .or. idiag_dby2m /= 0 .or. idiag_dbz2m /= 0) lpenc_diagnos(i_bbb) = .true.
      if (idiag_jmax /= 0 .or. idiag_jmin /= 0 .or. idiag_jrms /= 0 .or. &
          idiag_jm /= 0 .or. idiag_j2m /= 0) lpenc_diagnos(i_j2) = .true.
      if (idiag_jxmax /= 0 .or. idiag_jymax /= 0 .or. idiag_jzmax /= 0 .or. &
          idiag_jxm /= 0 .or. idiag_jym /= 0 .or. idiag_jzm /= 0 .or. &
          idiag_jx2m /= 0 .or. idiag_jy2m /= 0 .or. idiag_jz2m /= 0) lpenc_diagnos(i_jj) = .true.
      if (idiag_divbmax /= 0 .or. idiag_divbrms /= 0) lpenc_diagnos(i_divb) = .true.
      if (idiag_betamax /= 0 .or. idiag_betamin /= 0 .or. idiag_betam /= 0) lpenc_diagnos(i_beta) = .true.
      if (idiag_vAmax /= 0 .or. idiag_vAmin /= 0 .or. idiag_vAm /= 0) lpenc_diagnos(i_va2) = .true.
!
!  xy-averages related
!
      if (idiag_bmz /= 0 .or. idiag_b2mz /= 0) lpenc_diagnos(i_b2) = .true.
      if (idiag_bxmz /= 0 .or. idiag_bymz /= 0 .or. idiag_bzmz /= 0 .or. &
          idiag_bx2mz /= 0 .or. idiag_by2mz /= 0 .or. idiag_bz2mz /= 0 .or. &
          idiag_bxbymz /= 0 .or. idiag_bxbzmz /= 0 .or. idiag_bybzmz /= 0) lpenc_diagnos(i_bb) = .true.
      if (idiag_betamz /= 0 .or. idiag_beta2mz /= 0) lpenc_diagnos(i_beta) = .true.
!
!  yz-averages related
!
      if (idiag_bmx /= 0 .or. idiag_b2mx /= 0) lpenc_diagnos(i_b2) = .true.
      if (idiag_bxmx /= 0 .or. idiag_bymx /= 0 .or. idiag_bzmx /= 0 .or. &
          idiag_bx2mx /= 0 .or. idiag_by2mx /= 0 .or. idiag_bz2mx /= 0 .or. &
          idiag_bxbymx /= 0 .or. idiag_bxbzmx /= 0 .or. idiag_bybzmx /= 0) lpenc_diagnos(i_bb) = .true.
      if (idiag_betamx /= 0 .or. idiag_beta2mx /= 0) lpenc_diagnos(i_beta) = .true.
!
    endsubroutine pencil_criteria_magnetic
!***********************************************************************
    subroutine pencil_interdep_magnetic(lpencil_in)
!
!  Speicifies the dependency of the pencils provided by Magnetic.
!
!  25-oct-13/ccyang: coded.
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      if (lpencil_in(i_b2)) lpencil_in(i_bb) = .true.
!
      if (lpencil_in(i_bunit)) then
        lpencil_in(i_bb)=.true.
        lpencil_in(i_b2)=.true.
      endif
!
      divb: if (lpencil_in(i_divb)) then
        lpencil_in(i_bij) = .true.
        if (.not. lcartesian_coords) lpencil_in(i_bb) = .true.
      endif divb
!
      if (lpencil_in(i_j2)) lpencil_in(i_jj) = .true.
!
      jxbr: if (lpencil_in(i_jxbr)) then
        lpencil_in(i_jj) = .true.
        lpencil_in(i_bb) = .true.
        lpencil_in(i_rho1) = .true.
      endif jxbr
!
      beta1: if (lpencil_in(i_beta)) then
        lpencil_in(i_b2) = .true.
        lpencil_in(i_pp) = .true.
      endif beta1
!
      va2: if (lpencil_in(i_va2)) then
        lpencil_in(i_b2) = .true.
        lpencil_in(i_rho1) = .true.
      endif va2
!
    endsubroutine pencil_interdep_magnetic
!***********************************************************************
    subroutine magnetic_after_boundary(f)
!
!  Conducts any preprocessing required before the pencil calculations.
!
!  25-oct-13/ccyang: coded.
!
      use Boundcond, only: update_ghosts, zero_ghosts
      use Grid, only: get_grid_mn
      use Shear, only: get_uy0_shear
      use Sub, only: gij, curl_mn
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(nx,3,3) :: bij
      real, dimension(mx,3) :: uum, bbm
      real, dimension(nx,3) :: jj, bb
      real, dimension(mx) :: eta_penc
      real, dimension(3) :: b_ext
!
!  Replace B_ext locally to accommodate its time dependence.
!
      if (lbext) call get_bext(b_ext)
!
!  Reset maxdiffus_eta for time step constraint.
!
      if (lfirst .and. ldt) maxdiffus_eta = 0.0
!
!  Zero E field or initialize it with mesh hyper-resistivity.
!
      if (lresis_hyper3_mesh .or. (lshear .and. lhyper3x_mesh)) then
        call mesh_hyper_resistivity(f)
!       Note: The jj field contains some garbage after this call.
      else
        f(:,:,:,ieex:ieez) = 0.0
      endif
!
!  Find the current density J.
!
      mn_loop: do imn = 1, ny * nz
        n = nn(imn)
        m = mm(imn)
        call gij(f, ibb, bij, 1)
        cartesian: if (lcartesian_coords) then
          call curl_mn(bij, jj)
        else cartesian
          if (lbext) then
            bb = f(l1:l2,m,n,ibx:ibz) + spread(b_ext,1,nx)
          else
            bb = f(l1:l2,m,n,ibx:ibz)
          endif
          call curl_mn(bij, jj, bb)
        endif cartesian
        f(l1:l2,m,n,ijx:ijz) = jj
      enddo mn_loop
!
!  Communicate the J field.
!
      call zero_ghosts(f, ijx, ijz)
      call update_ghosts(f, ijx, ijz)
!
      resis: if (lexplicit_resistivity) then
!
!  Add normal resistivity.
!
        zscan1: do n = 1, mz
          yscan1: do m = 1, my
            call get_resistivity(f, eta_penc)
            f(:,m,n,ieex:ieez) = f(:,m,n,ieex:ieez) + spread(eta_penc, 2, 3) * f(:,m,n,ijx:ijz)
!           Time-step constraint
            timestep: if (lfirst .and. ldt) then
              if (.not. lcartesian_coords .or. .not. all(lequidist)) call get_grid_mn
              maxdiffus_eta = max(maxdiffus_eta, eta_penc(l1:l2) * dxyz_2)
            endif timestep
          enddo yscan1
        enddo zscan1
      endif resis
!
!  Give the J field correct units.
!
      f(:,:,:,ijx:ijz) = mu01 * f(:,:,:,ijx:ijz)
!
!  Get the shear velocity if it exists.
!
      if (lshear) call get_uy0_shear(uy0, x=x)
!
!  Add uu cross bb, including ghost cells.
!
      zscan2: do n = 1, mz
        yscan2: do m = 1, my
          if (lbext) then
            bbm = f(:,m,n,ibx:ibz) + spread(b_ext,1,mx)
          else
            bbm = f(:,m,n,ibx:ibz)
          endif
          if (iuu /= 0) then
            uum = f(:,m,n,iux:iuz)
          else
            uum = 0.0
          endif
          if (lshear) uum(:,2) = uum(:,2) + uy0
          f(:,m,n,ieex) = f(:,m,n,ieex) - (uum(:,2) * bbm(:,3) - uum(:,3) * bbm(:,2))
          f(:,m,n,ieey) = f(:,m,n,ieey) - (uum(:,3) * bbm(:,1) - uum(:,1) * bbm(:,3))
          f(:,m,n,ieez) = f(:,m,n,ieez) - (uum(:,1) * bbm(:,2) - uum(:,2) * bbm(:,1))
        enddo yscan2
      enddo zscan2
!
    endsubroutine magnetic_after_boundary
!***********************************************************************
    subroutine calc_pencils_magnetic_std(f,p)
!
!  Standard version (_std): global variable lpencil contains information about needed pencils.
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      type (pencil_case),                 intent(out)  :: p
!
      call calc_pencils_magnetic_pencpar(f,p,lpencil)
!
    endsubroutine calc_pencils_magnetic_std
!***********************************************************************
    subroutine calc_pencils_magnetic_pencpar(f,p,lpenc_loc)
!
!  Calculates Magnetic pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  25-jun-13/ccyang: coded.
!
      use Sub, only: gij, div_mn, curl_mn, cross, dot2_mn
!
      real, dimension(mx,my,mz,mfarray), intent(in)   :: f
      type(pencil_case),                 intent(inout):: p
      logical, dimension(:),             intent(in)   :: lpenc_loc
!
      real, dimension(nx,3,3) :: eij
      real, dimension(3) :: b_ext
!
      bb: if (lpenc_loc(i_bb)) then
        if (lbext) then
          call get_bext(b_ext)
          p%bb = f(l1:l2,m,n,ibx:ibz) + spread(b_ext,1,nx)
        else
          p%bb = f(l1:l2,m,n,ibx:ibz)
        endif
        p%bunit = p%bb
      endif bb
!
      if (lpenc_loc(i_bbb)) p%bbb = f(l1:l2,m,n,ibx:ibz)
!
      if (lpenc_loc(i_b2)) p%b2 = sum(p%bb**2, dim=2)
!
      if (lpenc_loc(i_bij)) call gij(f, ibb, p%bij, 1)
!
      if (lpenc_loc(i_jj)) p%jj = f(l1:l2,m,n,ijx:ijz)
!
      if (lpenc_loc(i_j2)) call dot2_mn(p%jj, p%j2)
!
      if (lpenc_loc(i_divb)) call div_mn(p%bij, p%divb, p%bb)
!
      jxbr: if (lpenc_loc(i_jxbr)) then
        call cross(p%jj, p%bb, p%jxbr)
        p%jxbr = spread(p%rho1,2,3) * p%jxbr
      endif jxbr
!
      curle: if (lpenc_loc(i_curle)) then
        call gij(f, iee, eij, 1)
        call curl_mn(eij, p%curle, f(:,m,n,ieex:ieez))
      endif curle
!
      if (lpenc_loc(i_beta)) p%beta = 2.0 * mu0 * p%pp / max(p%b2, tiny(1.0))
!
      if (lpenc_loc(i_va2)) p%va2 = mu01 * p%b2 * p%rho1
!
!  Dummy pencils
!
      if (lpenc_loc(i_aa)) call fatal_error('calc_pencils_magnetic', 'pencil aa is not implemented. ')
!
      if (lpenc_loc(i_ss12)) call fatal_error('calc_pencils_magnetic', 'pencil ss12 is not implemented. ')
!
    endsubroutine calc_pencils_magnetic_pencpar
!***********************************************************************
    subroutine daa_dt(f, df, p)
!
!  Evaluates the time derivative of the magnetic field.
!
!  25-oct-13/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(mx) :: eta_penc
!
!  dB/dt = -curl(E).
!
      df(l1:l2,m,n,ibx:ibz) =  df(l1:l2,m,n,ibx:ibz) - p%curle
!
!  Lorentz force
!
      if (lhydro) df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%jxbr
!
!  Ohmic heating
!
      ohmic: if (lresistivity .and. lenergy .and. lohmic_heat) then
        call get_resistivity(f, eta_penc, getall=.true.)
        eta_penc = mu0 * eta_penc
        eth: if (lthermal_energy) then
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + eta_penc(l1:l2) * p%j2
        else if (lentropy) then eth
          if (pretend_lnTT) then
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + eta_penc(l1:l2) * p%cv1 * p%j2 *p%rho1 * p%TT1
          else
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + eta_penc(l1:l2) * p%j2 * p%rho1 *p%TT1
          endif
        else if (ltemperature) then eth
          if (ltemperature_nolog) then
            df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + eta_penc(l1:l2) * p%cv1 * p%j2 * p%rho1
          else
            df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + eta_penc(l1:l2) * p%cv1 * p%j2 * p%rho1 * p%TT1
          endif
        endif eth
      endif ohmic
!
!  Constrain the time step.
!
      timestep: if (lfirst .and. ldt) then
        call set_advec_va2(p)
        maxdiffus = max(maxdiffus,maxdiffus_eta)
        maxdiffus3 = max(maxdiffus3,maxdiffus_eta3)
      endif timestep
!
!  Evaluate Magnetic diagnostics.
!
      if (ldiagnos) call diagnostic_magnetic(p)
!
      avg1d: if (l1davgfirst) then
        call xyaverages_magnetic(p)
        call xzaverages_magnetic(p)
        call yzaverages_magnetic(p)
      endif avg1d
!
    endsubroutine daa_dt
!******************************************************************************
    subroutine calc_diagnostics_magnetic(f,p)

      real, dimension(:,:,:,:) :: f
      type(pencil_case) :: p
!
      intent(in) :: f, p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_diagnostics_magnetic
!***********************************************************************
   subroutine magnetic_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine magnetic_after_timestep
!***********************************************************************
    subroutine split_update_magnetic(f)
!
! Calls for ImplicitDiffusion to implicitly evolve the resistivity
! term(s).
!
! 21-aug-13/ccyang: coded
!
      use ImplicitDiffusion, only: integrate_diffusion
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      if (limplicit_resistivity) call integrate_diffusion(get_resistivity_implicit, f, ibx, ibz)
!
    endsubroutine split_update_magnetic
!***********************************************************************
    subroutine read_magnetic_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=magnetic_init_pars, IOSTAT=iostat)
!
    endsubroutine read_magnetic_init_pars
!***********************************************************************
    subroutine write_magnetic_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=magnetic_init_pars)
!
    endsubroutine write_magnetic_init_pars
!***********************************************************************
    subroutine read_magnetic_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=magnetic_run_pars, IOSTAT=iostat)
!
    endsubroutine read_magnetic_run_pars
!***********************************************************************
    subroutine write_magnetic_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=magnetic_run_pars)
!
    endsubroutine write_magnetic_run_pars
!***********************************************************************
    subroutine rprint_magnetic(lreset, lwrite)
!
!  Reads and registers print parameters relevant for magnetic fields.
!
!  25-jun-13/ccyang: coded.
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
!
      logical, intent(in) :: lreset
      logical, intent(in), optional :: lwrite
!
      integer :: iname
!
!  Reset everything in case of RELOAD.
!
      reset: if (lreset) then
!
!       Diagnostics
!
        idiag_bmax = 0
        idiag_bmin = 0
        idiag_brms = 0
        idiag_bm = 0
        idiag_b2m = 0
        idiag_bxmax = 0
        idiag_bymax = 0
        idiag_bzmax = 0
        idiag_bxm = 0
        idiag_bym = 0
        idiag_bzm = 0
        idiag_bx2m = 0
        idiag_by2m = 0
        idiag_bz2m = 0
        idiag_bxbym = 0
        idiag_bxbzm = 0
        idiag_bybzm = 0
        idiag_dbxmax = 0
        idiag_dbymax = 0
        idiag_dbzmax = 0
        idiag_dbxm = 0
        idiag_dbym = 0
        idiag_dbzm = 0
        idiag_dbx2m = 0
        idiag_dby2m = 0
        idiag_dbz2m = 0
        idiag_jmax = 0
        idiag_jmin = 0
        idiag_jrms = 0
        idiag_jm = 0
        idiag_j2m = 0
        idiag_jxmax = 0
        idiag_jymax = 0
        idiag_jzmax = 0
        idiag_jxm = 0
        idiag_jym = 0
        idiag_jzm = 0
        idiag_jx2m = 0
        idiag_jy2m = 0
        idiag_jz2m = 0
        idiag_divbmax = 0
        idiag_divbrms = 0
        idiag_betamax = 0
        idiag_betamin = 0
        idiag_betam = 0
        idiag_vAmax = 0
        idiag_vAmin = 0
        idiag_vAm = 0
!
!       xy averages
!
        idiag_bmz = 0
        idiag_b2mz = 0
        idiag_bxmz = 0
        idiag_bymz = 0
        idiag_bzmz = 0
        idiag_bx2mz = 0
        idiag_by2mz = 0
        idiag_bz2mz = 0
        idiag_bxbymz = 0
        idiag_bxbzmz = 0
        idiag_bybzmz = 0
        idiag_betamz = 0
        idiag_beta2mz = 0
!
!       yz averages
!
        idiag_bmx = 0
        idiag_b2mx = 0
        idiag_bxmx = 0
        idiag_bymx = 0
        idiag_bzmx = 0
        idiag_bx2mx = 0
        idiag_by2mx = 0
        idiag_bz2mx = 0
        idiag_bxbymx = 0
        idiag_bxbzmx = 0
        idiag_bybzmx = 0
        idiag_betamx = 0
        idiag_beta2mx = 0
!
      endif reset
!
!  Parse the names from print.in.
!
      diag: do iname = 1, nname
        call parse_name(iname, cname(iname), cform(iname), 'bmax', idiag_bmax)
        call parse_name(iname, cname(iname), cform(iname), 'bmin', idiag_bmin)
        call parse_name(iname, cname(iname), cform(iname), 'brms', idiag_brms)
        call parse_name(iname, cname(iname), cform(iname), 'bm', idiag_bm)
        call parse_name(iname, cname(iname), cform(iname), 'b2m', idiag_b2m)
        call parse_name(iname, cname(iname), cform(iname), 'bxmax', idiag_bxmax)
        call parse_name(iname, cname(iname), cform(iname), 'bymax', idiag_bymax)
        call parse_name(iname, cname(iname), cform(iname), 'bzmax', idiag_bzmax)
        call parse_name(iname, cname(iname), cform(iname), 'bxm', idiag_bxm)
        call parse_name(iname, cname(iname), cform(iname), 'bym', idiag_bym)
        call parse_name(iname, cname(iname), cform(iname), 'bzm', idiag_bzm)
        call parse_name(iname, cname(iname), cform(iname), 'bx2m', idiag_bx2m)
        call parse_name(iname, cname(iname), cform(iname), 'by2m', idiag_by2m)
        call parse_name(iname, cname(iname), cform(iname), 'bz2m', idiag_bz2m)
        call parse_name(iname, cname(iname), cform(iname), 'bxbym', idiag_bxbym)
        call parse_name(iname, cname(iname), cform(iname), 'bxbzm', idiag_bxbzm)
        call parse_name(iname, cname(iname), cform(iname), 'bybzm', idiag_bybzm)
        call parse_name(iname, cname(iname), cform(iname), 'dbxmax', idiag_dbxmax)
        call parse_name(iname, cname(iname), cform(iname), 'dbymax', idiag_dbymax)
        call parse_name(iname, cname(iname), cform(iname), 'dbzmax', idiag_dbzmax)
        call parse_name(iname, cname(iname), cform(iname), 'dbxm', idiag_dbxm)
        call parse_name(iname, cname(iname), cform(iname), 'dbym', idiag_dbym)
        call parse_name(iname, cname(iname), cform(iname), 'dbzm', idiag_dbzm)
        call parse_name(iname, cname(iname), cform(iname), 'dbx2m', idiag_dbx2m)
        call parse_name(iname, cname(iname), cform(iname), 'dby2m', idiag_dby2m)
        call parse_name(iname, cname(iname), cform(iname), 'dbz2m', idiag_dbz2m)
        call parse_name(iname, cname(iname), cform(iname), 'jmax', idiag_jmax)
        call parse_name(iname, cname(iname), cform(iname), 'jmin', idiag_jmin)
        call parse_name(iname, cname(iname), cform(iname), 'jrms', idiag_jrms)
        call parse_name(iname, cname(iname), cform(iname), 'jm', idiag_jm)
        call parse_name(iname, cname(iname), cform(iname), 'j2m', idiag_j2m)
        call parse_name(iname, cname(iname), cform(iname), 'jxmax', idiag_jxmax)
        call parse_name(iname, cname(iname), cform(iname), 'jymax', idiag_jymax)
        call parse_name(iname, cname(iname), cform(iname), 'jzmax', idiag_jzmax)
        call parse_name(iname, cname(iname), cform(iname), 'jxm', idiag_jxm)
        call parse_name(iname, cname(iname), cform(iname), 'jym', idiag_jym)
        call parse_name(iname, cname(iname), cform(iname), 'jzm', idiag_jzm)
        call parse_name(iname, cname(iname), cform(iname), 'jx2m', idiag_jx2m)
        call parse_name(iname, cname(iname), cform(iname), 'jy2m', idiag_jy2m)
        call parse_name(iname, cname(iname), cform(iname), 'jz2m', idiag_jz2m)
        call parse_name(iname, cname(iname), cform(iname), 'divbmax', idiag_divbmax)
        call parse_name(iname, cname(iname), cform(iname), 'divbrms', idiag_divbrms)
        call parse_name(iname, cname(iname), cform(iname), 'betamax', idiag_betamax)
        call parse_name(iname, cname(iname), cform(iname), 'betamin', idiag_betamin)
        call parse_name(iname, cname(iname), cform(iname), 'betam', idiag_betam)
        call parse_name(iname, cname(iname), cform(iname), 'vAmax', idiag_vAmax)
        call parse_name(iname, cname(iname), cform(iname), 'vAmin', idiag_vAmin)
        call parse_name(iname, cname(iname), cform(iname), 'vAm', idiag_vAm)
      enddo diag
!
!  Parse the names from xy-averages.
!
      xyavg: do iname = 1, nnamez
        call parse_name(iname, cnamez(iname), cformz(iname), 'bmz', idiag_bmz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'b2mz', idiag_b2mz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'bxmz', idiag_bxmz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'bymz', idiag_bymz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'bzmz', idiag_bzmz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'bx2mz', idiag_bx2mz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'by2mz', idiag_by2mz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'bz2mz', idiag_bz2mz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'bxbymz', idiag_bxbymz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'bxbzmz', idiag_bxbzmz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'bybzmz', idiag_bybzmz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'betamz', idiag_betamz)
        call parse_name(iname, cnamez(iname), cformz(iname), 'beta2mz', idiag_beta2mz)
      enddo xyavg
!
!  Parse the names from yz-averages.
!
      yzavg: do iname = 1, nnamex
        call parse_name(iname, cnamex(iname), cformx(iname), 'bmx', idiag_bmx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'b2mx', idiag_b2mx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'bxmx', idiag_bxmx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'bymx', idiag_bymx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'bzmx', idiag_bzmx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'bx2mx', idiag_bx2mx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'by2mx', idiag_by2mx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'bz2mx', idiag_bz2mx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'bxbymx', idiag_bxbymx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'bxbzmx', idiag_bxbzmx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'bybzmx', idiag_bybzmx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'betamx', idiag_betamx)
        call parse_name(iname, cnamex(iname), cformx(iname), 'beta2mx', idiag_beta2mx)
      enddo yzavg
!
!  check for those quantities for which we want video slices
!       
      if (lwrite_slices) then 
        where(cnamev=='bb'.or.cnamev=='b2'.or.cnamev=='db2') cformv='DEFINED'
      endif
!       
!  Write variable indices for IDL.
!
      option: if (present(lwrite)) then
        indices: if (lwrite) then
          call farray_index_append('ibb',ibb)
          call farray_index_append('ibx',ibx)
          call farray_index_append('iby',iby)
          call farray_index_append('ibz',ibz)
          call farray_index_append('iEE',iEE)
          call farray_index_append('iEEx',iEEx)
          call farray_index_append('iEEy',iEEy)
          call farray_index_append('iEEz',iEEz)
          call farray_index_append('ijj',ijj)
          call farray_index_append('ijx',ijx)
          call farray_index_append('ijy',ijy)
          call farray_index_append('ijz',ijz)
        endif indices
      endif option
!
    endsubroutine rprint_magnetic
!***********************************************************************
    subroutine get_slices_magnetic(f, slices)
!
!  Prepares the slices requested by video.in.
!
!  26-jun-13/ccyang: coded.
!
      use Slices_methods, only: assign_slices_vec
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type(slice_data), intent(inout) :: slices
!
      integer :: ivar
      real, dimension(3) :: b_ext
!
!  Replace B_ext locally to accommodate its time dependence.
!
      if (lbext) call get_bext(b_ext)
!
      slice: select case (trim(slices%name))
!
!  Magnetic field
!
      case ('bb') slice
        call assign_slices_vec(slices,f,ibb)
!
!  Magnetic energy or fluctuating part of the magnetic energy
!
      case ('b2','db2') slice
        if (trim(slices%name)=='b2'.and.any(b_ext /= 0.0)) then
          if (lwrite_slice_yz)  slices%yz = sum((f(ix_loc,m1:m2,n1:n2,ibx:ibz) + spread(spread(b_ext,1,ny),2,nz))**2, dim=3)
          if (lwrite_slice_xz)  slices%xz = sum((f(l1:l2,iy_loc,n1:n2,ibx:ibz) + spread(spread(b_ext,1,nx),2,nz))**2, dim=3)
          if (lwrite_slice_xz2) slices%xz2 = sum((f(l1:l2,iy2_loc,n1:n2,ibx:ibz) + spread(spread(b_ext,1,nx),2,nz))**2, dim=3)
          if (lwrite_slice_xy)  slices%xy = sum((f(l1:l2,m1:m2,iz_loc,ibx:ibz) + spread(spread(b_ext,1,nx),2,ny))**2, dim=3)
          if (lwrite_slice_xy2) slices%xy2 = sum((f(l1:l2,m1:m2,iz2_loc,ibx:ibz) + spread(spread(b_ext,1,nx),2,ny))**2, dim=3)
          if (lwrite_slice_xy3) slices%xy3 = sum((f(l1:l2,m1:m2,iz3_loc,ibx:ibz) + spread(spread(b_ext,1,nx),2,ny))**2, dim=3)
          if (lwrite_slice_xy4) slices%xy4 = sum((f(l1:l2,m1:m2,iz4_loc,ibx:ibz) + spread(spread(b_ext,1,nx),2,ny))**2, dim=3)
        else
          if (lwrite_slice_yz)  slices%yz = sum(f(ix_loc,m1:m2,n1:n2,ibx:ibz)**2, dim=3)
          if (lwrite_slice_xz)  slices%xz = sum(f(l1:l2,iy_loc,n1:n2,ibx:ibz)**2, dim=3)
          if (lwrite_slice_xz2) slices%xz2 = sum(f(l1:l2,iy2_loc,n1:n2,ibx:ibz)**2, dim=3)
          if (lwrite_slice_xy)  slices%xy = sum(f(l1:l2,m1:m2,iz_loc,ibx:ibz)**2, dim=3)
          if (lwrite_slice_xy2) slices%xy2 = sum(f(l1:l2,m1:m2,iz2_loc,ibx:ibz)**2, dim=3)
          if (lwrite_slice_xy3) slices%xy3 = sum(f(l1:l2,m1:m2,iz3_loc,ibx:ibz)**2, dim=3)
          if (lwrite_slice_xy4) slices%xy4 = sum(f(l1:l2,m1:m2,iz4_loc,ibx:ibz)**2, dim=3)
        endif
        slices%ready = .true.
!
      endselect slice
!
    endsubroutine get_slices_magnetic
!***********************************************************************
    subroutine dynamical_resistivity(uc)
!
!  Dynamically set resistivity coefficient given fixed mesh Reynolds number.
!
!  19-may-14/ccyang: coded
!
!  Input Argument
!      uc
!          Characteristic velocity of the system.
!
      real, intent(in) :: uc
!
!  Mesh hyper-resistivity coefficient
!
      if (lresis_hyper3_mesh) eta_hyper3_mesh = pi5_1 * uc / re_mesh / sqrt(real(dimensionality))
!
    endsubroutine dynamical_resistivity
!***********************************************************************
    subroutine update_char_vel_magnetic(f)
!
!   Add the vector potential to the characteristic velocity
!   for slope limited diffusion.
!
!   25-sep-15/MR+joern: for slope limited diffusion
!
      real, dimension(mx,my,mz,mfarray), intent(inout):: f
!
      if (lslope_limit_diff) &
        f(:,:,:,iFF_diff2)=f(:,:,:,iFF_diff2)+sum(f(:,:,:,ibx:ibz)**2,4)
      
    endsubroutine update_char_vel_magnetic
!***********************************************************************
!***********************************************************************
!
!  LOCAL ROUTINES GO BELOW HERE.
!
!***********************************************************************
!***********************************************************************
    subroutine diagnostic_magnetic(p)
!
!  Evaluates the diagnostic variables.
!
!  25-jun-13/ccyang: coded.
!
      use Diagnostics, only: max_mn_name, sum_mn_name
!
      type(pencil_case), intent(in) :: p
!
      if (idiag_bmax /= 0) call max_mn_name(p%b2, idiag_bmax, lsqrt=.true.)
      if (idiag_bmin /= 0) call max_mn_name(-sqrt(p%b2), idiag_bmin, lneg=.true.)
      if (idiag_brms /= 0) call sum_mn_name(p%b2, idiag_brms, lsqrt=.true.)
      if (idiag_bm /= 0) call sum_mn_name(sqrt(p%b2), idiag_bm)
      if (idiag_b2m /= 0) call sum_mn_name(p%b2, idiag_b2m)
      if (idiag_bxmax /= 0) call max_mn_name(abs(p%bb(:,1)), idiag_bxmax)
      if (idiag_bymax /= 0) call max_mn_name(abs(p%bb(:,2)), idiag_bymax)
      if (idiag_bzmax /= 0) call max_mn_name(abs(p%bb(:,3)), idiag_bzmax)
      if (idiag_bxm /= 0) call sum_mn_name(p%bb(:,1), idiag_bxm)
      if (idiag_bym /= 0) call sum_mn_name(p%bb(:,2), idiag_bym)
      if (idiag_bzm /= 0) call sum_mn_name(p%bb(:,3), idiag_bzm)
      if (idiag_bx2m /= 0) call sum_mn_name(p%bb(:,1)**2, idiag_bx2m)
      if (idiag_by2m /= 0) call sum_mn_name(p%bb(:,2)**2, idiag_by2m)
      if (idiag_bz2m /= 0) call sum_mn_name(p%bb(:,3)**2, idiag_bz2m)
      if (idiag_bxbym /= 0) call sum_mn_name(p%bb(:,1) * p%bb(:,2), idiag_bxbym)
      if (idiag_bxbzm /= 0) call sum_mn_name(p%bb(:,1) * p%bb(:,3), idiag_bxbzm)
      if (idiag_bybzm /= 0) call sum_mn_name(p%bb(:,2) * p%bb(:,3), idiag_bybzm)
      if (idiag_dbxmax /= 0) call max_mn_name(abs(p%bbb(:,1)), idiag_dbxmax)
      if (idiag_dbymax /= 0) call max_mn_name(abs(p%bbb(:,2)), idiag_dbymax)
      if (idiag_dbzmax /= 0) call max_mn_name(abs(p%bbb(:,3)), idiag_dbzmax)
      if (idiag_dbxm /= 0) call sum_mn_name(p%bbb(:,1), idiag_dbxm)
      if (idiag_dbym /= 0) call sum_mn_name(p%bbb(:,2), idiag_dbym)
      if (idiag_dbzm /= 0) call sum_mn_name(p%bbb(:,3), idiag_dbzm)
      if (idiag_dbx2m /= 0) call sum_mn_name(p%bbb(:,1)**2, idiag_dbx2m)
      if (idiag_dby2m /= 0) call sum_mn_name(p%bbb(:,2)**2, idiag_dby2m)
      if (idiag_dbz2m /= 0) call sum_mn_name(p%bbb(:,3)**2, idiag_dbz2m)
      if (idiag_jmax /= 0) call max_mn_name(p%j2, idiag_jmax, lsqrt=.true.)
      if (idiag_jmin /= 0) call max_mn_name(-sqrt(p%j2), idiag_jmin, lneg=.true.)
      if (idiag_jrms /= 0) call sum_mn_name(p%j2, idiag_jrms, lsqrt=.true.)
      if (idiag_jm /= 0) call sum_mn_name(sqrt(p%j2), idiag_jm)
      if (idiag_j2m /= 0) call sum_mn_name(p%j2, idiag_j2m)
      if (idiag_jxmax /= 0) call max_mn_name(abs(p%jj(:,1)), idiag_jxmax)
      if (idiag_jymax /= 0) call max_mn_name(abs(p%jj(:,2)), idiag_jymax)
      if (idiag_jzmax /= 0) call max_mn_name(abs(p%jj(:,3)), idiag_jzmax)
      if (idiag_jxm /= 0) call sum_mn_name(p%jj(:,1), idiag_jxm)
      if (idiag_jym /= 0) call sum_mn_name(p%jj(:,2), idiag_jym)
      if (idiag_jzm /= 0) call sum_mn_name(p%jj(:,3), idiag_jzm)
      if (idiag_jx2m /= 0) call sum_mn_name(p%jj(:,1)**2, idiag_jx2m)
      if (idiag_jy2m /= 0) call sum_mn_name(p%jj(:,2)**2, idiag_jy2m)
      if (idiag_jz2m /= 0) call sum_mn_name(p%jj(:,3)**2, idiag_jz2m)
      if (idiag_divbmax /= 0) call max_mn_name(abs(p%divb), idiag_divbmax)
      if (idiag_divbrms /= 0) call sum_mn_name(p%divb**2, idiag_divbrms, lsqrt=.true.)
      if (idiag_betamax /= 0) call max_mn_name(p%beta, idiag_betamax)
      if (idiag_betamin /= 0) call max_mn_name(-p%beta, idiag_betamin, lneg=.true.)
      if (idiag_betam /= 0) call sum_mn_name(p%beta, idiag_betam)
      if (idiag_vAmax /= 0) call max_mn_name(p%va2, idiag_vAmax, lsqrt=.true.)
      if (idiag_vAmin /= 0) call max_mn_name(-sqrt(p%va2), idiag_vAmin, lneg=.true.)
      if (idiag_vAm /= 0) call sum_mn_name(sqrt(p%va2), idiag_vAm)
!
    endsubroutine diagnostic_magnetic
!***********************************************************************
    subroutine xyaverages_magnetic(p)
!
!  Accumulates the xy-averaged diagnostics.
!
!  28-oct-13/ccyang: coded.
!
      use Diagnostics, only: xysum_mn_name_z
!
      type(pencil_case), intent(in) :: p
!
      if (idiag_bmz /= 0) call xysum_mn_name_z(sqrt(p%b2), idiag_bmz)
      if (idiag_b2mz /= 0) call xysum_mn_name_z(p%b2, idiag_b2mz)
      if (idiag_bxmz /= 0) call xysum_mn_name_z(p%bb(:,1), idiag_bxmz)
      if (idiag_bymz /= 0) call xysum_mn_name_z(p%bb(:,2), idiag_bymz)
      if (idiag_bzmz /= 0) call xysum_mn_name_z(p%bb(:,3), idiag_bzmz)
      if (idiag_bx2mz /= 0) call xysum_mn_name_z(p%bb(:,1)**2, idiag_bx2mz)
      if (idiag_by2mz /= 0) call xysum_mn_name_z(p%bb(:,2)**2, idiag_by2mz)
      if (idiag_bz2mz /= 0) call xysum_mn_name_z(p%bb(:,3)**2, idiag_bz2mz)
      if (idiag_bxbymz /= 0) call xysum_mn_name_z(p%bb(:,1) * p%bb(:,2), idiag_bxbymz)
      if (idiag_bxbzmz /= 0) call xysum_mn_name_z(p%bb(:,1) * p%bb(:,3), idiag_bxbzmz)
      if (idiag_bybzmz /= 0) call xysum_mn_name_z(p%bb(:,2) * p%bb(:,3), idiag_bybzmz)
      if (idiag_betamz /= 0) call xysum_mn_name_z(p%beta, idiag_betamz)
      if (idiag_beta2mz /= 0) call xysum_mn_name_z(p%beta**2, idiag_beta2mz)
!
    endsubroutine xyaverages_magnetic
!***********************************************************************
    subroutine xzaverages_magnetic(p)
!
!  Accumulates the xz-averaged diagnostics.
!
!  28-oct-13/ccyang: coded.
!
      type(pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(p)
!
    endsubroutine xzaverages_magnetic
!***********************************************************************
    subroutine yzaverages_magnetic(p)
!
!  Accumulates the yz-averaged diagnostics.
!
!  28-oct-13/ccyang: coded.
!
      use Diagnostics, only: yzsum_mn_name_x
!
      type(pencil_case), intent(in) :: p
!
      if (idiag_bmx /= 0) call yzsum_mn_name_x(sqrt(p%b2), idiag_bmx)
      if (idiag_b2mx /= 0) call yzsum_mn_name_x(p%b2, idiag_b2mx)
      if (idiag_bxmx /= 0) call yzsum_mn_name_x(p%bb(:,1), idiag_bxmx)
      if (idiag_bymx /= 0) call yzsum_mn_name_x(p%bb(:,2), idiag_bymx)
      if (idiag_bzmx /= 0) call yzsum_mn_name_x(p%bb(:,3), idiag_bzmx)
      if (idiag_bx2mx /= 0) call yzsum_mn_name_x(p%bb(:,1)**2, idiag_bx2mx)
      if (idiag_by2mx /= 0) call yzsum_mn_name_x(p%bb(:,2)**2, idiag_by2mx)
      if (idiag_bz2mx /= 0) call yzsum_mn_name_x(p%bb(:,3)**2, idiag_bz2mx)
      if (idiag_bxbymx /= 0) call yzsum_mn_name_x(p%bb(:,1) * p%bb(:,2), idiag_bxbymx)
      if (idiag_bxbzmx /= 0) call yzsum_mn_name_x(p%bb(:,1) * p%bb(:,3), idiag_bxbzmx)
      if (idiag_bybzmx /= 0) call yzsum_mn_name_x(p%bb(:,2) * p%bb(:,3), idiag_bybzmx)
      if (idiag_betamx /= 0) call yzsum_mn_name_x(p%beta, idiag_betamx)
      if (idiag_beta2mx /= 0) call yzsum_mn_name_x(p%beta**2, idiag_beta2mx)
!
    endsubroutine yzaverages_magnetic
!***********************************************************************
    subroutine set_advec_va2(p)
!
!  Evaluates advec_va2 to constrain the time step.
!
!  25-jun-13/ccyang: coded.
!
      type(pencil_case), intent(in) :: p
      real, dimension(nx) :: advec_va2
!
      if (lspherical_coords) then
        advec_va2 = ((p%bb(:,1) * dx_1(l1:l2))**2 + &
                     (p%bb(:,2) * dy_1(m) * r1_mn)**2 + &
                     (p%bb(:,3) * dz_1(n) * r1_mn * sin1th(m))**2) * mu01 * p%rho1
      elseif (lcylindrical_coords) then
        advec_va2 = ((p%bb(:,1) * dx_1(l1:l2))**2 + &
                     (p%bb(:,2) * dy_1(m) * rcyl_mn1)**2 + &
                     (p%bb(:,3) * dz_1(n))**2) * mu01 * p%rho1
      else
        advec_va2 = ((p%bb(:,1) * dx_1(l1:l2))**2 + &
                     (p%bb(:,2) * dy_1(m))**2 + &
                     (p%bb(:,3) * dz_1(n))**2) * mu01 * p%rho1
      endif
      advec2=advec2+advec_va2
!
    endsubroutine set_advec_va2
!***********************************************************************
    subroutine get_resistivity(f, eta_penc, getall)
!
!  Gets the total normal resistivity along one pencil.
!  The unit of eta_penc is unit_length^2 / unit_time.
!
!  25-oct-13/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx), intent(out) :: eta_penc
      logical, intent(in), optional :: getall
!
      logical :: lall
!
!  Check if all resistivities are requested.
!
      if (present(getall)) then
          lall = getall .or. .not. limplicit_resistivity
      else
          lall = .not. limplicit_resistivity
      endif
!
!  Shock resistivity
!
      if (lresis_shock) then
        eta_penc = eta_shock * f(:,m,n,ishock)
      else
        eta_penc = 0.0
      endif
!
      explicit: if (lall) then
!
!  Constant resistivity
!
        if (lresis_const) eta_penc = eta_penc + eta
!
!  z-dependent resistivity
!
        if (lresis_zdep) eta_penc = eta_penc + eta_zdep(n)
!
      endif explicit
!
    endsubroutine get_resistivity
!***********************************************************************
    subroutine get_resistivity_implicit(ndc, diffus_coeff)
!
!  Gets the diffusion coefficient along a given pencil for the implicit algorithm.
!
!  21-aug-13/ccyang: coded.
!
      integer, intent(in) :: ndc
      real, dimension(ndc), intent(out) :: diffus_coeff
!
!  Constant resistivity
!
      if (lresis_const) then
        diffus_coeff = eta
      else
        diffus_coeff = 0.0
      endif
!
!  z-dependent resistivity
!
      if (lresis_zdep) call fatal_error('get_resistivity_implicit', &
                                        'z-dependent resistivity with implicit solver is under construction.')
!
    endsubroutine get_resistivity_implicit
!***********************************************************************
    subroutine mesh_hyper_resistivity(f)
!
!  Adds the mesh hyper-resistivity in divergence conserving form.
!
!  18-jul-14/ccyang: coded
!
      use Boundcond, only: update_ghosts, zero_ghosts
      use Grid, only: get_grid_mn
      use Sub, only: curl, del4, find_xyrms_fvec
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(nx,3) :: pv
      real, dimension(nx) :: eta3
      integer :: j
!
!  Reset maxdiffus_eta3 for time-step constraint.
!
      if (lfirst .and. ldt) maxdiffus_eta3 = 0.0
!
!  Calculate the mesh curl of B (assuming the boundary conditions have been applied).
!  Note: The auxiliary J field is temporarily used as working array here.
!
      curlbb: do imn = 1, ny * nz
        m = mm(imn)
        n = nn(imn)
        call curl(f, ibb, pv, ignoredx=.true.)
        f(l1:l2,m,n,ijx:ijz) = pv
      enddo curlbb
      call zero_ghosts(f, ijx, ijz)
      call update_ghosts(f, ijx, ijz)
!
!  Get the hyper resistivity.
!
      if (lshear .and. lhyper3x_mesh) then
        eta3 = eta_hyper3_mesh + diff_hyper3x_mesh * abs(Sshear) * xprim(l1:l2)
      else
        eta3 = eta_hyper3_mesh
      endif
!
!  Calculate the fourth-order mesh derivative of J.
!
      getd4jj: do imn = 1, ny * nz
        m = mm(imn)
        n = nn(imn)
        comp: do j = 1, 3
          call del4(f, ijj+j-1, pv(:,j), ignoredx=.true.)
          pv(:,j) = eta3 * pv(:,j)
        enddo comp
        f(l1:l2,m,n,ieex:ieez) = pv
!       Time-step constraint
        timestep: if (lfirst .and. ldt) then
          if (.not. lcartesian_coords .or. .not. all(lequidist)) call get_grid_mn
          maxdiffus_eta3 = max(maxdiffus_eta3, eta3 * (abs(dline_1(:,1)) + abs(dline_1(:,2)) + abs(dline_1(:,3))))
        endif timestep
      enddo getd4jj
!
!  Communicate the E field.
!
      call zero_ghosts(f, ieex, ieez)
      call update_ghosts(f, ieex, ieez)
!
    endsubroutine mesh_hyper_resistivity
!***********************************************************************
    elemental subroutine get_eta_zdep(z, eta, detadz)
!
! Finds the resistivity at given height z and specified profile
! zdep_prof.  If an unknown profile is given, zero resistivity is
! returned.
!
! 26-oct-13/ccyang: coded.
!
      use Sub, only: erfunc
      use EquationOfState, only: cs0
!
      real, intent(in) :: z
      real, intent(out), optional :: eta
      real, intent(out), optional :: detadz
!
      real :: h, zoh, a
!
      profile: select case (eta_zdep_prof)
!
      case ('FS03', 'FS2003') profile
!       Fleming, T., & Stone, J. M. 2003, ApJ, 585, 908
!         eta_zdep_coeff(1): eta at z = 0
!         eta_zdep_coeff(2): surface density ratio Sigma_0 / Sigma_CR, where CR stands for cosmic rays
        h = max(sqrt(2.0) * abs(cs0) / max(abs(Omega), tiny(1.0)), tiny(1.0))
        zoh = abs(z) / h
        a = eta_zdep_coeff(1) * exp(-0.5 * zoh**2 - 0.25 * eta_zdep_coeff(2) * erfunc(zoh))
        if (present(eta)) eta = a
        if (present(detadz)) detadz = -a * (zoh + 0.5 * exp(-zoh**2) / sqrtpi) / h
!
      case default profile
!       Unknown profile
        if (present(eta)) eta = 0.0
        if (present(detadz)) detadz = 0.0
!
      endselect profile
!
    endsubroutine get_eta_zdep
!***********************************************************************
!***********************************************************************
!
!  DUMMY BUT PUBLIC ROUTINES GO BELOW HERE.
!
!***********************************************************************
!***********************************************************************
    subroutine magnetic_before_boundary(f)
!
!  Conduct pre-processing required before boundary conditions and pencil
!  calculations.
!
!  29-may-14/ccyang: dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in):: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine magnetic_before_boundary
!***********************************************************************
    subroutine time_integrals_magnetic(f,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine time_integrals_magnetic
!***********************************************************************
    subroutine df_diagnos_magnetic(df,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  ::  df, p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine df_diagnos_magnetic
!***********************************************************************
    subroutine rescaling_magnetic(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine rescaling_magnetic
!***********************************************************************
    subroutine calc_mfield
!
!  Dummy routine
!
    endsubroutine calc_mfield
!***********************************************************************
    subroutine bb_unitvec_shock(f,bb_hat)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,3), intent (out) :: bb_hat
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bb_hat)
!
    endsubroutine bb_unitvec_shock
!***********************************************************************
    subroutine input_persistent_magnetic(id,done)
!
!  Dummy routine
!
      integer :: id
      logical :: done
!
      call keep_compiler_quiet(id)
      call keep_compiler_quiet(done)
!
    endsubroutine input_persistent_magnetic
!***********************************************************************
    logical function output_persistent_magnetic()
!
!  Dummy routine
!
      output_persistent_magnetic = .false.
!
    endfunction output_persistent_magnetic
!***********************************************************************
    subroutine expand_shands_magnetic()
!
!  Dummy
!
    endsubroutine expand_shands_magnetic
!***********************************************************************
    subroutine get_bext(B_ext_out)
!
!  Get the external magnetic field at current time step.
!
!  11-jun-14/ccyang: coded
!
      real, dimension(3), intent(out) :: B_ext_out
!
      bext: if (t_bext > 0.0 .and. t < t_bext) then
        if (t <= t0_bext) then
          B_ext_out = B0_ext
        else
          B_ext_out = B0_ext + 0.5 * (1.0 - cos(pi * (t - t0_bext) / (t_bext - t0_bext))) * (B_ext - B0_ext)
        endif
      else bext
        B_ext_out = B_ext
      endif bext
!
    endsubroutine get_bext
!***********************************************************************
endmodule Magnetic
