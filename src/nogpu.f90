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
  use Cparam
  use General, only: keep_compiler_quiet, lpointer, keep_compiler_quiet_dble

  implicit none

  logical :: ltest_bcs
  include 'gpu.h'

contains
!***********************************************************************
    subroutine initialize_GPU(f)
!
      real, dimension(:,:,:,:), intent(IN) :: f

      call keep_compiler_quiet(f)

    endsubroutine initialize_GPU
!**************************************************************************
    subroutine read_gpu_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_gpu_run_pars
!***********************************************************************
    subroutine write_gpu_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_gpu_run_pars
!***********************************************************************
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
    subroutine before_boundary_gpu(f,lrmv,itsub,t)
!
      real, dimension (:,:,:,:) :: f
      integer :: itsub
      logical :: early_finalize
      logical :: lrmv
      real(KIND=rkind8), intent(IN) :: t
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lrmv)
      call keep_compiler_quiet(itsub)
      call keep_compiler_quiet_dble(t)
!
    endsubroutine before_boundary_gpu
!**************************************************************************
    subroutine update_after_substep_gpu
    endsubroutine update_after_substep_gpu
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
    function get_ptr_GPU_training(ind1,ind2,lout) result(pFarr)

      integer :: ind1
      integer, optional :: ind2
      logical, optional :: lout

      real, dimension(:,:,:,:,:), pointer :: pFarr

      call keep_compiler_quiet(ind1,ind2)
      call keep_compiler_quiet(lout)
      call keep_compiler_quiet(pFarr)

    endfunction get_ptr_GPU_training
!**************************************************************************
    subroutine copy_farray_from_GPU(f,nowait_)

      real, dimension (:,:,:,:), intent(OUT) :: f
      logical, optional :: nowait_

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
    subroutine gpu_set_dt()

    endsubroutine gpu_set_dt
!**************************************************************************
    subroutine infer_gpu(flag)

    integer :: flag

    call keep_compiler_quiet(flag)

    endsubroutine infer_gpu
!**************************************************************************
    subroutine train_gpu(f, itsub)

    real :: f
    real :: itsub
    call keep_compiler_quiet(f)

    endsubroutine train_gpu
!**************************************************************************
    subroutine radtransfer_gpu
    endsubroutine
!**************************************************************************
    subroutine get_gpu_reduced_vars(dst)
      real, dimension(:) :: dst
      call keep_compiler_quiet(dst)
    endsubroutine get_gpu_reduced_vars
!**************************************************************************
    subroutine test_gpu_bcs
    endsubroutine test_gpu_bcs
!**************************************************************************
    subroutine split_update_gpu(f)
      real, dimension(:,:,:,:) :: f
    endsubroutine split_update_gpu
!**************************************************************************
    subroutine pushpars2c(p_par)
    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par
    endsubroutine pushpars2c
!**************************************************************************
endmodule GPU
