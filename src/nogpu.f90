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
  use General, only: keep_compiler_quiet, lpointer

  implicit none
  type(lpointer), dimension(1) :: lsnap_flags_to_wait_on
  type(lpointer), dimension(1) :: ldiag_flags_to_wait_on
  logical, target :: always_false_ng = .false.
  logical, target :: always_true_ng = .true.

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
      lsnap_flags_to_wait_on(1)%p => always_true_ng
      ldiag_flags_to_wait_on(1)%p => always_true_ng
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
    subroutine copy_farray_from_GPU(f,lflags_to_wait_on)

      real, dimension (:,:,:,:), intent(OUT) :: f
      type(lpointer), dimension(:) :: lflags_to_wait_on

      call keep_compiler_quiet(f)

    endsubroutine copy_farray_from_GPU
!**************************************************************************
    subroutine load_farray_to_GPU(f)

      real, dimension (:,:,:,:), intent(OUT) :: f

      call keep_compiler_quiet(f)

    endsubroutine  load_farray_to_GPU
!**************************************************************************
    subroutine test_rhs_gpu(f,df,p,mass_per_proc,early_finalize,cpu_version)
!
!  Used to test different implementations of rhs_cpu.
!
!  13-nov-23/TP: Written
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mfarray) :: df
      type (pencil_case) :: p
      real, dimension(1), intent(inout) :: mass_per_proc
      logical ,intent(in) :: early_finalize

      interface
        subroutine cpu_version(f,df,p,mass_per_proc,early_finalize)
          import mx
          import my
          import mz
          import mfarray
          import pencil_case
          real, dimension (mx,my,mz,mfarray) :: f
          real, dimension (mx,my,mz,mfarray) :: df
          type (pencil_case) :: p
          real, dimension(1), intent(inout) :: mass_per_proc
          logical ,intent(in) :: early_finalize

          intent(inout) :: f
          intent(inout) :: p
          intent(out) :: df

        endsubroutine cpu_version
      endinterface

    endsubroutine test_rhs_gpu
!**************************************************************************
endmodule  GPU
