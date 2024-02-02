! $Id$
!
! MODULE_DOC: This module contains GPU related dummy types and functions.
!
! CPARAM logical, parameter :: lgpu = .false.  !
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
 subroutine test_rhs_gpu(f,df,p,mass_per_proc,early_finalize,cpu_version)
!  Used to test different implementations of rhs_cpu.
!
!  13-nov-23/TP: Written
!
      ! use MPIcomm
      real, dimension (mx,my,mz,mfarray) :: f,f_copy
      real, dimension (mx,my,mz,mfarray) :: df,df_copy
      type (pencil_case) :: p,p_copy
      real, dimension(1), intent(inout) :: mass_per_proc
      logical ,intent(in) :: early_finalize
      integer :: i,j,k,n
      logical :: passed
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
  end subroutine  test_rhs_gpu
!**************************************************************************
endmodule  GPU
