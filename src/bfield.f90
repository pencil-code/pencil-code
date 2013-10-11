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
! MAUX CONTRIBUTION 3
! COMMUNICATED AUXILIARIES 3
!
! PENCILS PROVIDED bb(3); bbb(3); b2
! PENCILS PROVIDED bij(3,3); jj(3); j2; divb
! PENCILS PROVIDED curle(3); jxbr(3)
! PENCILS PROVIDED beta; va2
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
!
!  Initialization parameters
!
  real, dimension(3) :: b_ext = 0.0
!
  namelist /magnetic_init_pars/ b_ext
!
!  Runtime parameters
!
  logical :: lbext = .false.
  logical :: limplicit_resistivity = .false.
  logical :: lresis_const = .false.
  logical :: lresis_shock = .false.
  logical :: lresis_hyper3_mesh = .false.
  logical :: lohmic_heat = .true.
  logical, dimension(7) :: lresi_dep=.false.
  real :: eta = 0.0
  real :: eta_shock = 0.0
  real :: eta_hyper3_mesh = 0.0
!
  namelist /magnetic_run_pars/ b_ext, eta, eta_shock, eta_hyper3_mesh, limplicit_resistivity, lohmic_heat
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
!  Module variables
!
  real, dimension(nx) :: maxdiffus_eta = 0.0
  real, dimension(nx) :: maxdiffus_eta3 = 0.0
  real, dimension(mx) :: uy0 = 0.0
  logical :: lresistivity = .false.
  logical :: lexplicit_resistivity = .false.
!
!  Dummy but public variables (unfortunately)
!
  real, dimension(mz,3), parameter :: aamz = 0.0
  real, dimension(nz,3), parameter :: bbmz = 0.0, jjmz = 0.0
  real, dimension(3) :: b_ext_inv = 0.0
  logical, parameter :: lcalc_aameanz = .false.
  logical, parameter :: lelectron_inertia = .false.
  integer, parameter :: idiag_bcosphz = 0, idiag_bsinphz = 0
  integer, parameter :: idiag_axmz = 0, idiag_aymz = 0
  integer, parameter :: idiag_bxmz = 0, idiag_bymz = 0
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
!  24-jun-13/ccyang: coded.
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
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic(f, lstarting)
!
!  Conducts post-parameter-read initialization for Magnetic.
!
!  29-aug-13/ccyang: coded.
!
      use SharedVariables, only: put_shared_variable
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      logical, intent(in) :: lstarting
!
      integer :: ierr
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
!  Check the existence of external field.
!
      lbext = any(b_ext /= 0.0)
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
      if (eta_shock /= 0.0) lresis_shock = .true.
      if (eta_hyper3_mesh /= 0.0) lresis_hyper3_mesh = .true.
      lresistivity = lresis_const .and. lresis_shock
!
!  Sanity check
!
      if (lresis_shock .and. .not. lshock) &
        call fatal_error('initialize_magnetic', 'Shock module is required for shock resistivity. ')
!
!  Determine if any resistivity by explicit solver is present.
!
      lexplicit_resistivity = (.not. limplicit_resistivity .and. lresis_const) .or. lresis_shock
!
!  Information output
!
      resis: if (lroot) then
        if (lresis_const) print *, 'initialize_magnetic: constant resistivity, eta = ', eta
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
    endsubroutine pencil_criteria_magnetic
!***********************************************************************
    subroutine pencil_interdep_magnetic(lpencil_in)
!
!  Speicifies the dependency of the pencils provided by Magnetic.
!
!  24-jun-13/ccyang: coded.
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      if (lpencil_in(i_b2)) lpencil_in(i_bb) = .true.
!
      bder: if (lpencil_in(i_divb) .or. lpencil_in(i_jj)) then
        lpencil_in(i_bij) = .true.
        if (.not. lcartesian_coords) lpencil_in(i_bb) = .true.
      endif bder
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
    subroutine calc_lmagnetic_pars(f)
!
!  Conducts any preprocessing required before the pencil calculations.
!
!  08-oct-13/ccyang: coded.
!
      use Boundcond, only: update_ghosts
      use Grid, only: get_grid_mn
      use Shear, only: get_uy0_shear
      use Sub, only: gij, curl_mn
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(nx,3,3) :: bij
      real, dimension(mx,3) :: uum, bbm
      real, dimension(nx,3) :: jjn, bbn
      real, dimension(nx) :: eta_penc
      integer :: j, k
!
!  Reset maxdiffus_eta for time step constraint.
!
      if (lfirst .and. ldt) maxdiffus_eta = 0.0
!
!  Zero E field or initialize it with mesh hyper-resistivity.
!
      if (lresis_hyper3_mesh) then
        call mesh_hyper_resistivity(f)
      else
        f(:,:,:,ieex:ieez) = 0.0
      endif
!
!  Add normal resistivities.
!
      resis: if (lexplicit_resistivity) then
        mn_loop: do imn = 1, ny * nz
          n = nn(imn)
          m = mm(imn)
          if (lbext) then
            bbn = f(l1:l2,m,n,ibx:ibz) + spread(b_ext,1,nx)
          else
            bbn = f(l1:l2,m,n,ibx:ibz)
          endif
          call get_resistivity(f, eta_penc)
          call gij(f, ibb, bij, 1)
          call curl_mn(bij, jjn, bbn)
          f(l1:l2,m,n,ieex:ieez) = f(l1:l2,m,n,ieex:ieez) - spread(mu01 * eta_penc, 2, 3) * jjn
!         Time-step constraint
          timestep: if (lfirst .and. ldt) then
            if (.not. lcartesian_coords .or. .not. all(lequidist)) call get_grid_mn
            maxdiffus_eta = max(maxdiffus_eta, eta_penc * dxyz_2)
          endif timestep
        enddo mn_loop
      endif resis
!
!  Communicate the E field.
!
      if (lresis_hyper3_mesh .or. lexplicit_resistivity) call update_ghosts(f, ieex, ieez)
!
!  Get the shear velocity if existed.
!
      if (lshear) call get_uy0_shear(uy0, x=x)
!
!  Add uu cross bb, including ghost cells.
!
      uxb: if (iuu /= 0) then
        zscan: do j = 1, mz
          yscan: do k = 1, my
            if (lbext) then
              bbm = f(:,j,k,ibx:ibz) + spread(b_ext,1,mx)
            else
              bbm = f(:,j,k,ibx:ibz)
            endif
            uum = f(:,j,k,iux:iuz)
            if (lshear) uum(:,2) = uum(:,2) + uy0
            f(:,j,k,ieex) = f(:,j,k,ieex) + (uum(:,2) * bbm(:,3) - uum(:,3) * bbm(:,2))
            f(:,j,k,ieey) = f(:,j,k,ieey) + (uum(:,3) * bbm(:,1) - uum(:,1) * bbm(:,3))
            f(:,j,k,ieez) = f(:,j,k,ieez) + (uum(:,1) * bbm(:,2) - uum(:,2) * bbm(:,1))
          enddo yscan
        enddo zscan
      endif uxb
!
    endsubroutine calc_lmagnetic_pars
!***********************************************************************
    subroutine calc_pencils_magnetic(f, p)
!
!  Calculates Magnetic pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  25-jun-13/ccyang: coded.
!
      use Sub, only: gij, div_mn, curl_mn, cross, dot2_mn
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type(pencil_case), intent(inout) :: p
!
      real, dimension(nx,3,3) :: eij
!
      bb: if (lpencil(i_bb)) then
        if (lbext) then
          p%bb = f(l1:l2,m,n,ibx:ibz) + spread(b_ext,1,nx)
        else
          p%bb = f(l1:l2,m,n,ibx:ibz)
        endif
      endif bb
!
      if (lpencil(i_bbb)) p%bbb = f(l1:l2,m,n,ibx:ibz)
!
      if (lpencil(i_b2)) p%b2 = sum(p%bb**2, dim=2)
!
      if (lpencil(i_bij)) call gij(f, ibb, p%bij, 1)
!
      jj: if (lpencil(i_jj)) then
        call curl_mn(p%bij, p%jj, p%bb)
        p%jj = mu01 * p%jj
      endif jj
!
      if (lpencil(i_j2)) call dot2_mn(p%jj, p%j2)
!
      if (lpencil(i_divb)) call div_mn(p%bij, p%divb, p%bb)
!
      jxbr: if (lpencil(i_jxbr)) then
        call cross(p%jj, p%bb, p%jxbr)
        p%jxbr = spread(p%rho1,2,3) * p%jxbr
      endif jxbr
!
      curle: if (lpencil(i_curle)) then
        call gij(f, iee, eij, 1)
        call curl_mn(eij, p%curle, f(l1:l2,m,n,ieex:ieez))
      endif curle
!
      if (lpencil(i_beta)) p%beta = 2.0 * mu0 * p%pp / max(p%b2, tiny(1.0))
!
      if (lpencil(i_va2)) p%va2 = mu01 * p%b2 * p%rho1
!
!  Dummy pencils
!
      if (lpencil(i_aa)) call fatal_error('calc_pencils_magnetic', 'pencil aa is not implemented. ')
!
      if (lpencil(i_ss12)) call fatal_error('calc_pencils_magnetic', 'pencil ss12 is not implemented. ')
!
    endsubroutine calc_pencils_magnetic
!***********************************************************************
    subroutine daa_dt(f, df, p)
!
!  Evaluates the time derivative of the magnetic field.
!
!  09-jul-13/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type(pencil_case), intent(in) :: p
!
      real, dimension(nx) :: eta_penc
!
!  dB/dt = curl(E).
!
      df(l1:l2,m,n,ibx:ibz) =  df(l1:l2,m,n,ibx:ibz) + p%curle
!
!  Lorentz force
!
      if (lhydro) df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%jxbr
!
!  Ohmic heating
!
      ohmic: if (lresistivity .and. lenergy .and. lohmic_heat) then
        call get_resistivity(f, eta_penc)
        eth: if (lthermal_energy) then
          df(l1:l2,m,n,ieth) = df(l1:l2,m,n,ieth) + mu0 * eta_penc * p%j2
        else if (lentropy) then eth
          if (pretend_lnTT) then
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + mu0 * eta_penc * p%cv1 * p%j2 *p%rho1 * p%TT1
          else
            df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + mu0 * eta_penc * p%j2 * p%rho1 *p%TT1
          endif
        else if (ltemperature) then eth
          if (ltemperature_nolog) then
            df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + mu0 * eta_penc * p%cv1 * p%j2 * p%rho1
          else
            df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + mu0 * eta_penc * p%cv1 * p%j2 * p%rho1 * p%TT1
          endif
        endif eth
      endif ohmic
!
!  Constrain the time step.
!
      timestep: if (lfirst .and. ldt) then
        call set_advec_va2(p)
        if (lshear) advec_shear = abs(uy0(l1:l2) * dy_1(m))
        diffus_eta = maxdiffus_eta
        diffus_eta3 = maxdiffus_eta3
      endif timestep
!
!  Evaluate Magnetic diagnostics.
!
      if (ldiagnos) call diagnostic_magnetic(p)
!
    endsubroutine daa_dt
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
    subroutine read_magnetic_init_pars(unit, iostat)
!
! Reads the initialization parameters for Magnetic.
!
! 19-jun-13/ccyang: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      integer :: stat
!
      read(unit, NML=magnetic_init_pars, IOSTAT=stat)
      if (present(iostat)) then
        iostat = stat
      else if (stat /= 0) then
        call fatal_error('read_magnetic_init_pars', 'cannot read magnetic_init_pars. ')
      endif
!
    endsubroutine read_magnetic_init_pars
!***********************************************************************
    subroutine write_magnetic_init_pars(unit)
!
! Writes the initialization parameters for Magnetic.
!
! 19-jun-13/ccyang: coded
!
      integer, intent(in) :: unit
!
      integer :: stat
!
      write(unit, NML=magnetic_init_pars, IOSTAT=stat)
      if (stat /= 0) call fatal_error('write_magnetic_init_pars', 'cannot write magnetic_init_pars. ')
!
    endsubroutine write_magnetic_init_pars
!***********************************************************************
    subroutine read_magnetic_run_pars(unit, iostat)
!
! Reads the runtime parameters for Magnetic.
!
! 19-jun-13/ccyang: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      integer :: stat
!
      read(unit, NML=magnetic_run_pars, IOSTAT=stat)
      if (present(iostat)) then
        iostat = stat
      else if (stat /= 0) then
        call fatal_error('read_magnetic_run_pars', 'cannot read magnetic_run_pars. ')
      endif
!
    endsubroutine read_magnetic_run_pars
!***********************************************************************
    subroutine write_magnetic_run_pars(unit)
!
! Writes the runtime parameters for Magnetic.
!
! 19-jun-13/ccyang: coded
!
      integer, intent(in) :: unit
!
      integer :: stat
!
      write(unit, NML=magnetic_run_pars, IOSTAT=stat)
      if (stat /= 0) call fatal_error('write_magnetic_run_pars', 'cannot write magnetic_run_pars. ')
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
!
      logical, intent(in) :: lreset
      logical, intent(in), optional :: lwrite
!
      integer :: iname
!
!  Reset everything in case of RELOAD.
!
      reset: if (lreset) then
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
!  Write variable indices for IDL.
!
      option: if (present(lwrite)) then
        indices: if (lwrite) then
          write(3,*) 'ibb = ', ibb
          write(3,*) 'ibx = ', ibx
          write(3,*) 'iby = ', iby
          write(3,*) 'ibz = ', ibz
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
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type(slice_data), intent(inout) :: slices
!
      integer :: ivar
!
      slice: select case (trim(slices%name))
!
!  Magnetic field
!
      case ('bb') slice
        bcomp: if (slices%index >= 3) then
          slices%ready=.false.
        else bcomp
          ivar = ibb + slices%index
          slices%index = slices%index + 1
          slices%yz = f(ix_loc,m1:m2,n1:n2,ivar)
          slices%xz = f(l1:l2,iy_loc,n1:n2,ivar)
          slices%xy = f(l1:l2,m1:m2,iz_loc,ivar)
          slices%xy2 = f(l1:l2,m1:m2,iz2_loc,ivar)
          if (lwrite_slice_xy3) slices%xy3 = f(l1:l2,m1:m2,iz3_loc,ivar)
          if (lwrite_slice_xy4) slices%xy4 = f(l1:l2,m1:m2,iz4_loc,ivar)
          if (slices%index <= 3) slices%ready = .true.
        endif bcomp
!
!  Magnetic energy
!
      case ('b2') slice
        slices%yz = sum((f(ix_loc,m1:m2,n1:n2,ibx:ibz) + spread(spread(b_ext,1,ny),2,nz))**2, dim=3)
        slices%xz = sum((f(l1:l2,iy_loc,n1:n2,ibx:ibz) + spread(spread(b_ext,1,nx),2,nz))**2, dim=3)
        slices%xy = sum((f(l1:l2,m1:m2,iz_loc,ibx:ibz) + spread(spread(b_ext,1,nx),2,ny))**2, dim=3)
        slices%xy2 = sum((f(l1:l2,m1:m2,iz2_loc,ibx:ibz) + spread(spread(b_ext,1,nx),2,ny))**2, dim=3)
        if (lwrite_slice_xy3) slices%xy3 = sum((f(l1:l2,m1:m2,iz3_loc,ibx:ibz) + spread(spread(b_ext,1,nx),2,ny))**2, dim=3)
        if (lwrite_slice_xy4) slices%xy4 = sum((f(l1:l2,m1:m2,iz4_loc,ibx:ibz) + spread(spread(b_ext,1,nx),2,ny))**2, dim=3)
        slices%ready = .true.
!
!  Fluctuating part of the magnetic energy
!
      case ('db2') slice
        slices%yz = sum(f(ix_loc,m1:m2,n1:n2,ibx:ibz)**2, dim=3)
        slices%xz = sum(f(l1:l2,iy_loc,n1:n2,ibx:ibz)**2, dim=3)
        slices%xy = sum(f(l1:l2,m1:m2,iz_loc,ibx:ibz)**2, dim=3)
        slices%xy2 = sum(f(l1:l2,m1:m2,iz2_loc,ibx:ibz)**2, dim=3)
        if (lwrite_slice_xy3) slices%xy3 = sum(f(l1:l2,m1:m2,iz3_loc,ibx:ibz)**2, dim=3)
        if (lwrite_slice_xy4) slices%xy4 = sum(f(l1:l2,m1:m2,iz4_loc,ibx:ibz)**2, dim=3)
        slices%ready = .true.
!
      endselect slice
!
    endsubroutine get_slices_magnetic
!***********************************************************************
    subroutine dynamical_resistivity(umax)
!
!  Dynamically set resistivity coefficient given fixed mesh Reynolds number.
!
!  09-jul-13/ccyang: coded
!
      real, intent(in) :: umax
!
!  Mesh hyper-resistivity coefficient
!
      if (lresis_hyper3_mesh) eta_hyper3_mesh = pi5_1 * umax / re_mesh
!
    endsubroutine dynamical_resistivity
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
    subroutine set_advec_va2(p)
!
!  Evaluates advec_va2 to constrain the time step.
!
!  25-jun-13/ccyang: coded.
!
      type(pencil_case), intent(in) :: p
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
!
    endsubroutine set_advec_va2
!***********************************************************************
    subroutine get_resistivity(f, eta_penc)
!
!  Gets the total normal resistivity along one pencil.
!
!  21-aug-13/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: eta_penc
!
!  Constant resistivity
!
      if (lresis_const .and. .not. limplicit_resistivity) then
        eta_penc = eta
      else
        eta_penc = 0.0
      endif
!
!  Shock resistivity
!
      if (lresis_shock) eta_penc = eta_penc + eta_shock * f(l1:l2,m,n,ishock)
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
      if (lresis_const) then
        diffus_coeff = mu01 * eta
      else
        diffus_coeff = 0.0
      endif
!
    endsubroutine get_resistivity_implicit
!***********************************************************************
    subroutine mesh_hyper_resistivity(f)
!
!  Adds the mesh hyper-resistivity in divergence conserving form.
!
!  20-sep-13/ccyang: coded
!
      use Boundcond, only: update_ghosts
      use Deriv, only: der4
      use Grid, only: get_grid_mn
      use Sub, only: curl
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(nx,ny,nz,3) :: d4jj
      real, dimension(nx,3) :: pv
      integer :: j, k, mc, nc
!
!  Reset maxdiffus_eta3 for time-step constraint.
!
      if (lfirst .and. ldt) maxdiffus_eta3 = 0.0
!
!  Calculate the mesh curl of B (assuming the boundary conditions have been applied).
!
      curlbb: do imn = 1, ny * nz
        m = mm(imn)
        n = nn(imn)
        call curl(f, ibb, pv, ignoredx=.true.)
        f(l1:l2,m,n,ieex:ieez) = pv
      enddo curlbb
      call update_ghosts(f, ieex, ieez)
!
!  Calculate its fourth-order mesh derivative.
!
      d4jj = 0.0
      d4jjp: do imn = 1, ny * nz
        m = mm(imn)
        n = nn(imn)
        mc = m - nghost
        nc = n - nghost
        dir: do k = 1, 3
          comp: do j = 1, 3
            call der4(f, iee+j-1, pv(:,j), k, ignoredx=.true.)
          enddo comp
          d4jj(:,mc,nc,:) = d4jj(:,mc,nc,:) + pv
        enddo dir
!       Time-step constraint
        timestep: if (lfirst .and. ldt) then
          if (.not. lcartesian_coords .or. .not. all(lequidist)) call get_grid_mn
          maxdiffus_eta3 = max(maxdiffus_eta3, eta_hyper3_mesh * (abs(dline_1(:,1)) + abs(dline_1(:,2)) + abs(dline_1(:,3))))
        endif timestep
      enddo d4jjp
!
!  Assign it to the E field.
!
      f(l1:l2,m1:m2,n1:n2,ieex:ieez) = eta_hyper3_mesh * d4jj
!
    endsubroutine mesh_hyper_resistivity
!***********************************************************************
!***********************************************************************
!
!  DUMMY BUT PUBLIC ROUTINES GO BELOW HERE.
!
!***********************************************************************
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
endmodule Magnetic
