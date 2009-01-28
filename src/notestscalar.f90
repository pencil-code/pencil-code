! $Id$

!  This modules deals with all aspects of testscalar fields; if no
!  testscalar fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testscalar relevant subroutines listed in here.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Testscalar

  use Cparam
  use Messages

  implicit none

  include 'testscalar.h'

  real :: dummy=0.

  namelist /testscalar_init_pars/ &
       dummy
  namelist /testscalar_run_pars/ &
       dummy

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
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (ip==0) print*,f  !(to keep compiler quiet)
!
    endsubroutine initialize_testscalar
!***********************************************************************
    subroutine init_cctest(f,xx,yy,zz)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if (ip==0) print*,f,xx,yy,zz  !(to keep compiler quiet)
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
      use Cdata
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in
!
    endsubroutine pencil_interdep_testscalar
!***********************************************************************
    subroutine read_testscalar_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_testscalar_init_pars
!***********************************************************************
    subroutine write_testscalar_init_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit !(keep compiler quiet)

    endsubroutine write_testscalar_init_pars
!***********************************************************************
    subroutine read_testscalar_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_testscalar_run_pars
!***********************************************************************
    subroutine write_testscalar_run_pars(unit)
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine write_testscalar_run_pars
!***********************************************************************
    subroutine dcctest_dt(f,df,p)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)     :: f,p
      intent(inout)  :: df
!
      if (ip==0) print*, f, df, p  !(to keep compiler quiet)
!
    endsubroutine dcctest_dt
!***********************************************************************
    subroutine get_slices_testscalar(f,slices)
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
    endsubroutine get_slices_testscalar
!***********************************************************************
    subroutine calc_ltestscalar_pars(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: c,s
!
      if (NO_WARN) print*, f
!
    endsubroutine calc_ltestscalar_pars
!***********************************************************************
    subroutine rescaling_testscalar(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine rescaling_testscalar
!***********************************************************************
    subroutine rprint_testscalar(lreset,lwrite)
!
!  Dummy routine
!
      use Cdata
!
      logical :: lreset
      logical, optional :: lwrite
!
      if (ip==0) print*,lreset,lwrite  !(to keep compiler quiet)
    endsubroutine rprint_testscalar

endmodule Testscalar
