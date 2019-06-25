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
    subroutine gpu_init
!
    endsubroutine gpu_init
!**************************************************************************
    subroutine register_GPU(f)
!
      real, dimension(:,:,:,:), intent(IN) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine register_GPU
!**************************************************************************
    subroutine finalize_GPU
!
    endsubroutine finalize_GPU
!**************************************************************************
    subroutine rhs_GPU(f,itsub,early_finalize)
!
      real, dimension (:,:,:,:) :: f
      integer :: itsub
      logical :: early_finalize
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(itsub)
      call keep_compiler_quiet(early_finalize)
!
    endsubroutine rhs_GPU
!**************************************************************************
    subroutine copy_farray_from_GPU(f)

      real, dimension (:,:,:,:), intent(OUT) :: f

      call keep_compiler_quiet(f)

    endsubroutine copy_farray_from_GPU
!**************************************************************************
endmodule  GPU
