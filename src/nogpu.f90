! $Id$
!
! MODULE_DOC: This module contains GPU related dummy types and functions.
!
! CPARAM logical, parameter :: lgpu = .false.
!
!**************************************************************************
!
module GPU
!
  use Cdata
  use General, only: keep_compiler_quiet

  implicit none

  include 'gpu.h'

contains

!***********************************************************************
    subroutine initialize_GPU
!
    endsubroutine initialize_GPU
!**************************************************************************
    subroutine finalize_GPU
!
    endsubroutine finalize_GPU
!**************************************************************************
    subroutine rhs_GPU(f,itsub,lsnap)
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: itsub
      logical :: lsnap
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(itsub)
      call keep_compiler_quiet(lsnap)
!
    endsubroutine rhs_GPU
!**************************************************************************
endmodule  GPU
