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

  include 'gpu.h'

contains
!***********************************************************************
    subroutine initialize_GPU(f)
!
      real, dimension(:,:,:,:), intent(IN) :: f

      call keep_compiler_quiet(f)

    endsubroutine initialize_GPU
!**************************************************************************
    subroutine gpu_init
    endsubroutine gpu_init
!**************************************************************************
    subroutine register_GPU
    endsubroutine register_GPU
!**************************************************************************
    subroutine finalize_GPU
    endsubroutine finalize_GPU
!**************************************************************************
    subroutine get_farray_ptr_gpu
    endsubroutine get_farray_ptr_gpu
!**************************************************************************
    subroutine rhs_GPU(f,itsub)
!
      real, dimension (:,:,:,:) :: f
      integer :: itsub
      logical :: early_finalize
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(itsub)
!
    endsubroutine rhs_GPU
!**************************************************************************
    function get_ptr_GPU(ind1,ind2,lout) result(pFarr)

      integer :: ind1
      integer, optional :: ind2
      logical, optional :: lout

      real, dimension(:,:,:,:), pointer :: pFarr

      call keep_compiler_quiet(ind1,ind2)
      call keep_compiler_quiet(lout)
      call keep_compiler_quiet(pFarr)

    endfunction get_ptr_GPU
!**************************************************************************
    function get_ptr_GPU_for_training(ind1,ind2,lout) result(pFarr)

      integer :: ind1
      integer, optional :: ind2
      logical, optional :: lout

      real, dimension(:,:,:,:,:), pointer :: pFarr

      call keep_compiler_quiet(ind1,ind2)
      call keep_compiler_quiet(lout)
      call keep_compiler_quiet(pFarr)

    endfunction get_ptr_GPU_for_training
!**************************************************************************
    subroutine copy_farray_from_GPU(f)

      real, dimension (:,:,:,:), intent(OUT) :: f

      call keep_compiler_quiet(f)

    endsubroutine copy_farray_from_GPU
!**************************************************************************
    subroutine load_farray_to_GPU(f)

      real, dimension (:,:,:,:), intent(OUT) :: f

      call keep_compiler_quiet(f)

    endsubroutine  load_farray_to_GPU
!**************************************************************************
    subroutine reload_GPU_config

    endsubroutine reload_GPU_config
!**************************************************************************
    subroutine update_on_gpu(index,varname,value)

      integer :: index
      character(LEN=*), optional :: varname
      real, optional :: value

      call keep_compiler_quiet(index)
      if (present(varname)) call keep_compiler_quiet(varname)
      if (present(value)) call keep_compiler_quiet(value)

    endsubroutine update_on_gpu
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

      external cpu_version

    endsubroutine test_rhs_gpu
!**************************************************************************
    subroutine gpu_set_dt()
    endsubroutine gpu_set_dt
!**************************************************************************
    subroutine infer_c(flag)
    integer :: flag
    endsubroutine infer_c
!**************************************************************************
    subroutine train_c(f)
    real :: f
    endsubroutine train_c
!**************************************************************************

endmodule  GPU
