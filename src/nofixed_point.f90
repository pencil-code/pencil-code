! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
module Fixed_point
!
  use Cdata
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
! a few constants
  integer :: FINISHED_FIXED = 98
! the arrays with the values for x, y, and z for all cores (tracer xyz)
  real, pointer, dimension (:) :: xt, yt, zt
! fixed points array
  real, dimension (1000,3) :: fixed_points
! total number of fixed points of this core
  integer :: fidx
  real :: fidx_read
!
  real, public :: tfixed_points
  integer, public :: nfixed_points
! Time value to be written together with the fixed points.
  real :: tfixed_points_write
!
  contains
!
!***********************************************************************
  subroutine fixed_points_prepare()
!
!  dummy
!
  endsubroutine fixed_points_prepare
!***********************************************************************
  subroutine wfixed_points(f,path)
!
!  dummy
!
    real, dimension (mx,my,mz,mfarray) :: f
    character(len=*) :: path
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(path)
!
  endsubroutine wfixed_points
!***********************************************************************

endmodule Fixed_point
