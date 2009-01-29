! $Id$

!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Testfield

  use Cparam

  implicit none

  include 'testfield.h'

  real :: dummy=0.

  namelist /testfield_init_pars/ &
       dummy
  namelist /testfield_run_pars/ &
       dummy

  contains

!***********************************************************************
    subroutine register_testfield()
!
!  Dummy routine
!
    endsubroutine register_testfield
!***********************************************************************
    subroutine initialize_testfield(f)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (ip==0) print*,f  !(to keep compiler quiet)
    endsubroutine initialize_testfield
!***********************************************************************
    subroutine init_aatest(f)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (ip==0) print*,f  !(to keep compiler quiet)
!
    endsubroutine init_aatest
!***********************************************************************
    subroutine pencil_criteria_testfield()
!
!  Dummy routine
!
    endsubroutine pencil_criteria_testfield
!***********************************************************************
    subroutine pencil_interdep_testfield(lpencil_in)
!
!  Dummy routine
!
      use Cdata
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in
!
    endsubroutine pencil_interdep_testfield
!***********************************************************************
    subroutine read_testfield_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_testfield_init_pars
!***********************************************************************
    subroutine write_testfield_init_pars(unit)
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine write_testfield_init_pars
!***********************************************************************
    subroutine read_testfield_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_testfield_run_pars
!***********************************************************************
    subroutine write_testfield_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit !(keep compiler quiet)

    endsubroutine write_testfield_run_pars
!***********************************************************************
    subroutine daatest_dt(f,df,p)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)     :: f, p
      intent(inout)  :: df
!
      if (ip==0) print*, f, df, p  !(to keep compiler quiet)
!
    endsubroutine daatest_dt
!***********************************************************************
    subroutine get_slices_testfield(f,slices)
!
      use Sub, only: keep_compiler_quiet
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_testfield
!***********************************************************************
    subroutine calc_ltestfield_pars(f)
!
!  29-jan-06/axel: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in)     :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine calc_ltestfield_pars
!***********************************************************************
    subroutine rescaling_testfield(f)
!
!  18-may-08/axel: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine rescaling_testfield
!***********************************************************************
    subroutine rprint_testfield(lreset,lwrite)
!
!  Dummy routine
!
      use Cdata
!
      logical :: lreset
      logical, optional :: lwrite
!
      if (ip==0) print*,lreset,lwrite  !(to keep compiler quiet)
    endsubroutine rprint_testfield
!***********************************************************************

endmodule Testfield
