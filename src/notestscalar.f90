! $Id$
!
!  This modules deals with all aspects of testscalar fields; if no
!  testscalar fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testscalar relevant subroutines listed in here.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestscalar = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Testscalar
!
  use Cparam
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'testscalar.h'
!
  real :: dummy=0.
!
  namelist /testscalar_init_pars/ &
       dummy
  namelist /testscalar_run_pars/ &
       dummy
!
  contains
!***********************************************************************
    subroutine register_testscalar()
!
!  Dummy routine
!
    endsubroutine register_testscalar
!***********************************************************************
    subroutine initialize_testscalar(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_testscalar
!***********************************************************************
    subroutine init_cctest(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_cctest
!***********************************************************************
    subroutine pencil_criteria_testscalar()
!
!  Dummy routine
!
    endsubroutine pencil_criteria_testscalar
!***********************************************************************
    subroutine pencil_interdep_testscalar(lpencil_in)
!
!  Dummy routine
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_testscalar
!***********************************************************************
    subroutine read_testscalar_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_testscalar_init_pars
!***********************************************************************
    subroutine write_testscalar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_testscalar_init_pars
!***********************************************************************
    subroutine read_testscalar_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_testscalar_run_pars
!***********************************************************************
    subroutine write_testscalar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_testscalar_run_pars
!***********************************************************************
    subroutine dcctest_dt(f,df,p)
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
    endsubroutine dcctest_dt
!***********************************************************************
    subroutine get_slices_testscalar(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_testscalar
!***********************************************************************
    subroutine testscalar_after_boundary(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine testscalar_after_boundary
!***********************************************************************
    subroutine rescaling_testscalar(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine rescaling_testscalar
!***********************************************************************
    subroutine rprint_testscalar(lreset,lwrite)
!
!  Dummy routine
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_testscalar
!***********************************************************************
endmodule Testscalar
