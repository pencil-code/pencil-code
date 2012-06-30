! $Id$
!
!  This module produces slices for animation purposes.
!
module Slices
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  private
!
  public :: wvid, wvid_prepare, setup_slices, wslice
!
  real, public :: tvid=0.0
  integer, public :: nvid=0
  real :: tslice=0.0
!
  contains
!***********************************************************************
    subroutine wvid_prepare
!
!  23-nov-09/anders: dummy
!
    endsubroutine wvid_prepare
!***********************************************************************
    subroutine wvid(f,path)
!
!  23-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      character(len=*) :: path
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(path)
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
    subroutine setup_slices()
!
!  23-nov-09/anders: dummy
!
    endsubroutine setup_slices
!***********************************************************************
endmodule Slices
