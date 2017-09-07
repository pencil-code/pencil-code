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
    subroutine rhs_GPU(f,itsub)
!
      real, dimension (:,:,:,:) :: f
      integer :: itsub
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(itsub)
!
    endsubroutine rhs_GPU
!**************************************************************************
endmodule  GPU
