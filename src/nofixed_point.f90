! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
module Fixed_point
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  public :: fixed_points_prepare, get_fixed_point, edge, pindex
  public :: get_fixed_points, wfixed_points
!
  contains
!
!***********************************************************************
    subroutine fixed_points_prepare()
!
!    dummy
!
    endsubroutine fixed_points_prepare
!***********************************************************************
    subroutine get_fixed_point(f,point,fixed_point,q,vv)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: point(2), fixed_point(2)
      real :: q
      real, pointer, dimension (:,:,:,:) :: vv
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(point)
      call keep_compiler_quiet(fixed_point)
      call keep_compiler_quiet(q)
      call keep_compiler_quiet(vv)
!
    endsubroutine get_fixed_point
!***********************************************************************
    subroutine edge(f,sx,sy,diff1,diff2,phi_min,vv,rec)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: sx(2), sy(2)
      real diff1(2), diff2(2)
      real :: phi_min
      real, pointer, dimension (:,:,:,:) :: vv
      integer :: rec
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sx)
      call keep_compiler_quiet(sy)
      call keep_compiler_quiet(diff1)
      call keep_compiler_quiet(diff2)
      call keep_compiler_quiet(phi_min)
      call keep_compiler_quiet(vv)
      call keep_compiler_quiet(rec)
!
    endsubroutine edge
!***********************************************************************
    subroutine pindex(f,sx,sy,diff,phi_min,vv,poincare)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: sx(2), sy(2)
      real diff(4,2)
      real :: phi_min
      real, pointer, dimension (:,:,:,:) :: vv
      real :: poincare
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(sx)
      call keep_compiler_quiet(sy)
      call keep_compiler_quiet(diff)
      call keep_compiler_quiet(phi_min)
      call keep_compiler_quiet(vv)
      call keep_compiler_quiet(poincare)
!
    endsubroutine pindex
!***********************************************************************
    subroutine get_fixed_points(f,tracers,vv)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, pointer, dimension (:,:) :: tracers
      real, pointer, dimension (:,:,:,:) :: vv
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(tracers)
      call keep_compiler_quiet(vv)
!
    endsubroutine get_fixed_points
!***********************************************************************
    subroutine wfixed_points(f,path)
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
