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

module Testflow

  use Cparam
  use Messages

  implicit none

  include 'testflow.h'

  real :: dummy
  namelist /testflow_init_pars/ dummy

  namelist /testflow_run_pars/ dummy

  contains

!***********************************************************************
    subroutine register_testflow()
!
!  Dummy routine
!
      use Cdata
!
    endsubroutine register_testflow
!***********************************************************************
    subroutine initialize_testflow(f)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (ip==0) print*,f  !(to keep compiler quiet)
    endsubroutine initialize_testflow
!***********************************************************************
    subroutine init_uutest(f,xx,yy,zz)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      if (ip==0) print*,f,xx,yy,zz  !(to keep compiler quiet)
    endsubroutine init_uutest
!***********************************************************************
    subroutine pencil_criteria_testflow()
!
!  Dummy routine
!
      use Cdata
!
    endsubroutine pencil_criteria_testflow
!***********************************************************************
    subroutine pencil_interdep_testflow(lpencil_in)
!
!  Dummy routine
!
      use Cdata
!
      logical, dimension(npencils) :: lpencil_in
!
    endsubroutine pencil_interdep_testflow
!***********************************************************************
    subroutine read_testflow_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_testflow_init_pars
!***********************************************************************
    subroutine write_testflow_init_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit !(keep compiler quiet)

    endsubroutine write_testflow_init_pars
!***********************************************************************
    subroutine read_testflow_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)

    endsubroutine read_testflow_run_pars
!***********************************************************************
    subroutine write_testflow_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit !(keep compiler quiet)

    endsubroutine write_testflow_run_pars
!***********************************************************************
    subroutine duutest_dt(f,df,p)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      intent(in)     :: f,p
      intent(inout)  :: df
!
      if (ip==0) print*, f, df, p  !(to keep compiler quiet)
!
    endsubroutine duutest_dt
!***********************************************************************
    subroutine get_slices_testflow(f,slices)
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
    endsubroutine get_slices_testflow
!***********************************************************************
    subroutine calc_ltestflow_nonlin_terms(f,df)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      intent(inout) :: f,df
!
      if (NO_WARN) print*, f
!
    endsubroutine calc_ltestflow_nonlin_terms
!***********************************************************************
    subroutine rprint_testflow(lreset,lwrite)
!
!  Dummy routine
!
      use Cdata
!
      logical :: lreset
      logical, optional :: lwrite
!
      if (ip==0) print*,lreset,lwrite  !(to keep compiler quiet)
    endsubroutine rprint_testflow

endmodule Testflow
