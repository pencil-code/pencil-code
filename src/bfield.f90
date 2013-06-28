! $Id$
!
!  This module directly evolves the magnetic field instead of
!  the vector potential; special care needs to be taken to
!  guarantee the divergence of the field remains zero.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagnetic = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 3
! COMMUNICATED AUXILIARIES 3
!
! PENCILS PROVIDED bb(3); bbb(3); b2
! PENCILS PROVIDED bij(3,3); jj(3); divb
! PENCILS PROVIDED curle(3); jxbr(3)
!
! PENCILS PROVIDED ss12
!
!***************************************************************
module Magnetic
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include 'magnetic.h'
!
!  Module variables
!
  real, dimension(3) :: b_ext = 0.0
  logical :: lbext = .false.
!
!  Initialization parameters
!
  namelist /magnetic_init_pars/ b_ext
!
!  Runtime parameters
!
  namelist /magnetic_run_pars/ b_ext
!
!  Diagnostic variables
!
  integer :: idiag_bmax = 0     ! DIAG_DOC: $\max|\mathbf(B)|$
  integer :: idiag_brms = 0     ! DIAG_DOC: $\left<B^2\right>^{1/2}$
  integer :: idiag_bxmax = 0    ! DIAG_DOC: $\max|B_x|$
  integer :: idiag_bymax = 0    ! DIAG_DOC: $\max|B_y|$
  integer :: idiag_bzmax = 0    ! DIAG_DOC: $\max|B_z|$
  integer :: idiag_dbxmax = 0   ! DIAG_DOC: $\max|\Delta B_x|$
  integer :: idiag_dbymax = 0   ! DIAG_DOC: $\max|\Delta B_y|$
  integer :: idiag_dbzmax = 0   ! DIAG_DOC: $\max|\Delta B_z|$
  integer :: idiag_divbmax = 0  ! DIAG_DOC: $\max|\nabla\cdot\mathbf{B}|$
!
!  Dummy but public variables (unfortunately)
!
  real, dimension(mz,3), parameter :: aamz = 0.0
  real, dimension(nz,3), parameter :: bbmz = 0.0, jjmz = 0.0
  real, dimension(3) :: b_ext_inv = 0.0
  logical, parameter :: lcalc_aamean = .false.
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
!  20-jun-13/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      logical, intent(in) :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
!  Check the existence of external field.
!
      lbext = any(b_ext /= 0.0)
!
!  Calculates variables required by Magnetic and possibly other module(s).
!
      mu01 = 1.0 / mu0
!
!  Calculates variables that are public and currently used by other module(s).
!
      if (lbext) b_ext_inv = b_ext / sum(b_ext**2)
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
      if (lhydro) lpenc_requested(i_jxbr) = .true.
!
      time_step: if (ldt) then
        lpenc_requested(i_bb) = .true.
        lpenc_requested(i_rho1) = .true.
      endif time_step
!
!  Diagnostic related
!
      if (idiag_bmax /= 0 .or. idiag_brms /= 0) lpenc_diagnos(i_b2) = .true.
      if (idiag_bxmax /= 0 .or. idiag_bymax /= 0 .or. idiag_bzmax /= 0) lpenc_diagnos(i_bb) = .true.
      if (idiag_dbxmax /= 0 .or. idiag_dbymax /= 0 .or. idiag_dbzmax /= 0) lpenc_diagnos(i_bbb) = .true.
      if (idiag_divbmax /= 0) lpenc_diagnos(i_divb) = .true.
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
    endsubroutine pencil_interdep_magnetic
!***********************************************************************
    subroutine calc_lmagnetic_pars(f)
!
!  Conducts any preprocessing required before the pencil calculations.
!
!  27-jun-13/ccyang: coded.
!
      use Boundcond, only: update_ghosts
      use Sub, only: cross
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(nx,3) :: bb, ee
      integer :: imn, m, n
!
!  Update ghost cells of the velocity and the magnetic fields.
!
      if (iuu /= 0) call update_ghosts(f, iux, iuz)
      call update_ghosts(f, ibx, ibz)
!
!  Calculate the effective electric field along each pencil.
!
      mn_loop: do imn = 1, ny * nz
        n = nn(imn)
        m = mm(imn)
!
!  Add uu cross bb.
!
        uxb: if (iuu /= 0) then
          if (lbext) then
            bb = f(l1:l2,m,n,ibx:ibz) + spread(b_ext,1,nx)
          else
            bb = f(l1:l2,m,n,ibx:ibz)
          endif
          call cross(f(l1:l2,m,n,iux:iuz), bb, ee)
        else uxb
          ee = 0.0
        endif uxb
!
        f(l1:l2,m,n,ieex:ieez) = ee
      enddo mn_loop
!
!  Communicate the E field.
!
      call update_ghosts(f, ieex, ieez)
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
      use Sub, only: gij, div_mn, curl_mn, cross
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
!  Dummy pencils
!
      if (lpencil(i_ss12)) call fatal_error('calc_pencils_magnetic', 'pencil ss12 is not implemented. ')
!
    endsubroutine calc_pencils_magnetic
!***********************************************************************
    subroutine daa_dt(f, df, p)
!
!  Evaluates the time derivative of the magnetic field.
!
!  25-jun-13/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type(pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f)
!
!  dB/dt = curl(E).
!
      df(l1:l2,m,n,ibx:ibz) =  df(l1:l2,m,n,ibx:ibz) + p%curle
!
!  Lorentz force
!
      if (lhydro) df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%jxbr
!
!  Constrain the time step.
!
      if (lfirst .and. ldt) call set_advec_va2(p)
!
!  Evaluate Magnetic diagnostics.
!
      if (ldiagnos) call diagnostic_magnetic(p)
!
    endsubroutine daa_dt
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
        idiag_brms = 0
        idiag_bxmax = 0
        idiag_bymax = 0
        idiag_bzmax = 0
        idiag_dbxmax = 0
        idiag_dbymax = 0
        idiag_dbzmax = 0
        idiag_divbmax = 0
      endif reset
!
!  Parse the names from print.in.
!
      diag: do iname = 1, nname
        call parse_name(iname, cname(iname), cform(iname), 'bmax', idiag_bmax)
        call parse_name(iname, cname(iname), cform(iname), 'brms', idiag_brms)
        call parse_name(iname, cname(iname), cform(iname), 'bxmax', idiag_bxmax)
        call parse_name(iname, cname(iname), cform(iname), 'bymax', idiag_bymax)
        call parse_name(iname, cname(iname), cform(iname), 'bzmax', idiag_bzmax)
        call parse_name(iname, cname(iname), cform(iname), 'dbxmax', idiag_dbxmax)
        call parse_name(iname, cname(iname), cform(iname), 'dbymax', idiag_dbymax)
        call parse_name(iname, cname(iname), cform(iname), 'dbzmax', idiag_dbzmax)
        call parse_name(iname, cname(iname), cform(iname), 'divbmax', idiag_divbmax)
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
      if (idiag_brms /= 0) call sum_mn_name(p%b2, idiag_brms, lsqrt=.true.)
      if (idiag_bxmax /= 0) call max_mn_name(abs(p%bb(:,1)), idiag_bxmax)
      if (idiag_bymax /= 0) call max_mn_name(abs(p%bb(:,2)), idiag_bymax)
      if (idiag_bzmax /= 0) call max_mn_name(abs(p%bb(:,3)), idiag_bzmax)
      if (idiag_dbxmax /= 0) call max_mn_name(abs(p%bbb(:,1)), idiag_dbxmax)
      if (idiag_dbymax /= 0) call max_mn_name(abs(p%bbb(:,2)), idiag_dbymax)
      if (idiag_dbzmax /= 0) call max_mn_name(abs(p%bbb(:,3)), idiag_dbzmax)
      if (idiag_divbmax /= 0) call max_mn_name(abs(p%divb), idiag_divbmax)
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
    subroutine dynamical_resistivity(umax)
!
!  dummy
!
      real, intent(in) :: umax
!     
      call keep_compiler_quiet(umax)
!
    endsubroutine dynamical_resistivity
!***********************************************************************
    subroutine split_update_magnetic(f)
!
!  dummy
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine split_update_magnetic
!***********************************************************************
    subroutine expand_shands_magnetic()
!
!  Dummy
!
    endsubroutine expand_shands_magnetic
!***********************************************************************
endmodule Magnetic
