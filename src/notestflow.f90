! $Id$

!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestflow = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Testflow
!
  use Cparam
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'testflow.h'
!
  real :: dummy
  namelist /testflow_init_pars/ dummy
!
  namelist /testflow_run_pars/ dummy
!
  contains
!***********************************************************************
    subroutine register_testflow()
!
!  Dummy routine
!
    endsubroutine register_testflow
!***********************************************************************
    subroutine initialize_testflow(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_testflow
!***********************************************************************
    subroutine init_uutest(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_uutest
!***********************************************************************
    subroutine pencil_criteria_testflow()
!
!  Dummy routine
!
    endsubroutine pencil_criteria_testflow
!***********************************************************************
    subroutine pencil_interdep_testflow(lpencil_in)
!
!  Dummy routine
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_testflow
!***********************************************************************
    subroutine read_testflow_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_testflow_init_pars
!***********************************************************************
    subroutine write_testflow_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_testflow_init_pars
!***********************************************************************
    subroutine read_testflow_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_testflow_run_pars
!***********************************************************************
    subroutine write_testflow_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_testflow_run_pars
!***********************************************************************
    subroutine duutest_dt(f,df,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)     :: f,p
      intent(inout)  :: df
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine duutest_dt
!***********************************************************************
    subroutine get_slices_testflow(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_testflow
!***********************************************************************
    subroutine testflow_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!   15-dec-10/MR: adapted from density
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine testflow_before_boundary
!***********************************************************************
    subroutine calc_ltestflow_nonlin_terms(f,df)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      intent(in)    :: f
      intent(inout) :: df
!
      call keep_compiler_quiet(f,df)
!
    endsubroutine calc_ltestflow_nonlin_terms
!***********************************************************************
    subroutine rprint_testflow(lreset,lwrite)
!
!  Dummy routine
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_testflow
!***********************************************************************
endmodule Testflow 
