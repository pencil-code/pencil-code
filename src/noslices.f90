! $Id$
!
!  This module produces slices for animation purposes.
!
module Slices
!
  use Quiet
!
  implicit none
!
  private
!
  public :: wvid, wvid_prepare, setup_slices, wslice
!
  contains
!***********************************************************************
    subroutine wvid_prepare
!
    endsubroutine wvid_prepare
!***********************************************************************
    subroutine wvid(f)
!
!  23-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine wvid
!***********************************************************************
    subroutine wslice(filename,a,pos,ndim1,ndim2)
!
!  23-nov-09/anders: dummy
!
      integer :: ndim1,ndim2
      character (len=*) :: filename

      real, dimension (ndim1,ndim2) :: a
      real, intent(in) :: pos
!
      call keep_compiler_quiet(filename)
      call keep_compiler_quiet(a)
      call keep_compiler_quiet(pos)
      call keep_compiler_quiet(ndim1)
      call keep_compiler_quiet(ndim2)
!
    endsubroutine wslice
!***********************************************************************
    subroutine setup_slices
!
    endsubroutine setup_slices
!***********************************************************************
endmodule Slices
