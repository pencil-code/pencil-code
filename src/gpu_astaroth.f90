! $Id$
!
! MODULE_DOC: This module contains GPU related types and functions to be used with the ASTAROTH nucleus.
!
! CPARAM logical, parameter :: lgpu = .true.
!
!**************************************************************************
!
module GPU
!
  use Cdata
  use General, only: keep_compiler_quiet, lpointer, ioptest, loptest
  use Messages
!$ use, intrinsic :: iso_c_binding
  use iso_c_binding

  implicit none

  include 'gpu.h'

!$  interface
!$    subroutine random_initial_condition() bind(C)
!$    endsubroutine random_initial_condition
!$  end interface
  
  external initialize_gpu_c
  external register_gpu_c
  external finalize_gpu_c
  external get_farray_ptr_gpu_c
  external rhs_gpu_c
  external before_boundary_gpu_c
  external after_timestep_gpu_c 
  external source_function_and_opacity_gpu_c
  external load_farray_c
  external reload_gpu_config_c
  external test_rhs_c
  external copy_farray_c
  external update_on_gpu_arr_by_ind_c
  external update_on_gpu_scal_by_ind_c
  external pos_real_ptr_c
  external gpu_set_dt_c
  external torchtrain_c 
  external torchinfer_c
  external calcq_gpu_c
  external get_gpu_reduced_vars_c
  external test_bcs_c

  integer, external :: update_on_gpu_arr_by_name_c
  integer, external :: update_on_gpu_scal_by_name_c

  !integer(KIND=ikind8) :: pFarr_GPU_in, pFarr_GPU_out
  type(C_PTR) :: pFarr_GPU_in, pFarr_GPU_out

  namelist /gpu_run_pars/ &
     ltest_bcs,lac_sparse_autotuning,lcpu_timestep_on_gpu,lcumulative_df_on_gpu,lread_all_vars_from_device,lcuda_aware_mpi

contains
!***********************************************************************
  subroutine train_gpu(loss)

    real :: loss
    call torchtrain_c(loss)

  endsubroutine train_gpu 
!***********************************************************************
  subroutine infer_gpu(flag)

    integer :: flag
    call torchinfer_c(flag)

  endsubroutine infer_gpu
!***********************************************************************
    subroutine initialize_GPU(f)
!
      use Mpicomm, only: MPI_COMM_PENCIL

      real, dimension(:,:,:,:), intent(IN) :: f

      character(LEN=512) :: str
!
      str=''
      if (lanelastic) str=trim(str)//', '//'anelastic'
      if (lboussinesq) str=trim(str)//', '//'boussinesq'
      if (lhyperresistivity_strict) str=trim(str)//', '//'hyperresi_strict'
      if (lhyperviscosity_strict) str=trim(str)//', '//'hypervisc_strict'
      if (lADI) str=trim(str)//', '//'implicit_physics'
      if (llorenz_gauge) str=trim(str)//', '//'lorenz_gauge'
      if (ltestscalar) str=trim(str)//', '//'testscalar'
      if (ltestfield) str=trim(str)//', '//'testfield'
      if (ltestflow) str=trim(str)//', '//'testflow'
      if (ldetonate) str=trim(str)//', '//'detonate'
      if (lopacity) str=trim(str)//', '//'opacity'
      if (lpointmasses) str=trim(str)//', '//'pointmasses'
      if (lsolid_cells) str=trim(str)//', '//'solid_cells'
      if (lparticles) str=trim(str)//', '//'particles'

      if (str/='') call fatal_error('initialize_GPU','no GPU implementation available for module(s) "'// &
                                    trim(str(3:))//'"')
!
      if (dt<=0.) dt = dtmin
      call initialize_gpu_c(f,MPI_COMM_PENCIL,t,nt)
!
! Load farray to gpu
!
      if (nt>0) call load_farray_to_GPU(f)

  !print'(a,1x,Z0,1x,Z0)', 'pFarr_GPU_in,pFarr_GPU_out=', pFarr_GPU_in,pFarr_GPU_out
    endsubroutine initialize_GPU
!**************************************************************************
    subroutine read_gpu_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=gpu_run_pars, IOSTAT=iostat)
!
    endsubroutine read_gpu_run_pars 
!***********************************************************************
    subroutine write_gpu_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=gpu_run_pars)
!
    endsubroutine write_gpu_run_pars
!**************************************************************************
    subroutine register_GPU
!
      call register_gpu_c
!
    endsubroutine register_GPU
!**************************************************************************
    subroutine finalize_gpu
!
      call finalize_gpu_c
!
    endsubroutine finalize_GPU
!**************************************************************************
    subroutine get_farray_ptr_gpu

      call get_farray_ptr_gpu_c(pFarr_GPU_in)

    endsubroutine get_farray_ptr_gpu
!**************************************************************************
    subroutine rhs_GPU(f,isubstep)
!
      use General, only: notanumber

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f
      integer,                            intent(IN)    :: isubstep
!
      call rhs_gpu_c(isubstep,t)
!
    endsubroutine rhs_GPU
!**************************************************************************
    subroutine before_boundary_gpu(f,lrmv,isubstep,t)
!
      use General, only: notanumber

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f
      integer,                            intent(IN)    :: isubstep
      real(KIND=rkind8), intent(IN) :: t
      logical :: lrmv
      integer :: lrmv_int
!
      !TP: pass int since integers are more compatible with C than logical to booleans
      if (lrmv) then
        lrmv_int = 1
      else
        lrmv_int = 0
      endif
      call before_boundary_gpu_c(lrmv_int,isubstep,t)
!
    endsubroutine before_boundary_gpu
!**************************************************************************
    subroutine after_timestep_gpu
      call after_timestep_gpu_c
    endsubroutine after_timestep_gpu
!**************************************************************************
    subroutine gpu_set_dt
!
      call gpu_set_dt_c(t)
!
    endsubroutine gpu_set_dt
!**************************************************************************
    function get_ptr_GPU(ind1,ind2,lout) result(pFarr)
!
!  Fetches the address of the f-array counterpart on the GPU for slots from ind1 to ind2
!  and transforms it to a Fortran pointer.
!
      integer :: ind1
      integer, optional :: ind2
      logical, optional :: lout

      real, dimension(:,:,:,:), pointer :: pFarr

      integer :: i2

      interface
        type(c_ptr) function pos_real_ptr_c(ptr,ind)
          import :: c_ptr, ikind8
          type(c_ptr) :: ptr
          integer :: ind
        endfunction
      endinterface

      i2 = ioptest(ind2,ind1)
      if (loptest(lout)) then
        call c_f_pointer(pos_real_ptr_c(pFarr_GPU_out,ind1-1),pFarr,(/mx,my,mz,i2-ind1+1/))
      else
        call c_f_pointer(pos_real_ptr_c(pFarr_GPU_in,ind1-1),pFarr,(/mx,my,mz,i2-ind1+1/))
      endif

    endfunction get_ptr_GPU
!**************************************************************************
    function get_ptr_GPU_training(ind1,ind2,lout) result(pFarr)

! For training, pointed-to array needs to be 5-dimensional (last dimension being the batch size)

      use Cparam
      use iso_c_binding

      integer :: ind1
      integer, optional :: ind2
      logical, optional :: lout

      real, dimension(:,:,:,:,:), pointer :: pFarr

      integer :: i2

      interface
        type(c_ptr) function pos_real_ptr_c(ptr,ind)
          import :: c_ptr, ikind8
          type(c_ptr) :: ptr
          integer :: ind
        endfunction
      endinterface

      i2 = ioptest(ind2,ind1)
      if (loptest(lout)) then
        call c_f_pointer(pos_real_ptr_c(pFarr_GPU_out,ind1-1),pFarr,(/mx,my,mz,i2-ind1+1,1/))
      else
        call c_f_pointer(pos_real_ptr_c(pFarr_GPU_in,ind1-1),pFarr,(/mx,my,mz,i2-ind1+1,1/))
      endif

    endfunction get_ptr_GPU_training
!**************************************************************************
    subroutine copy_farray_from_GPU(f,nowait_)

!$    use General, only: signal_wait

      real, dimension (mx,my,mz,mfarray), intent(OUT) :: f
      logical, optional :: nowait_
      logical :: nowait
      integer :: i

      nowait = loptest(nowait_)
      if (nowait) then
        call copy_farray_c(f)
        return
      endif
!
!$    if (lfarray_copied) return
!
! Have to wait since if doing diagnostics don't want to overwrite f.
!
!$    call signal_wait(lhelper_perf, .false.)
      call copy_farray_c(f)
!$    lfarray_copied = .true.

    endsubroutine copy_farray_from_GPU
!**************************************************************************
    subroutine load_farray_to_GPU(f)

      real, dimension (mx,my,mz,mfarray), intent(IN) :: f

      call load_farray_c

    endsubroutine load_farray_to_GPU
!**************************************************************************
    subroutine reload_GPU_config

      call reload_gpu_config_c

    endsubroutine reload_GPU_config
!**************************************************************************
    subroutine update_on_gpu(index, varname, value)
!
!  Updates an element of the Astaroth configuration, identified by name or index, on the GPU.
!
      integer, intent(inout) :: index
      character(LEN=*),optional :: varname
      real, optional :: value

      if (index>=0) then
        if (present(value)) then
          call update_on_gpu_scal_by_ind_c(index,value)
        else
          call update_on_gpu_arr_by_ind_c(index)
        endif
      else
        if (present(value)) then
          index = update_on_gpu_scal_by_name_c(varname//char(0),value)
        else
          index = update_on_gpu_arr_by_name_c(varname//char(0))
        endif
        if (index<0) call fatal_error('update_on_gpu','variable "'//trim(varname)//'" not found')
      endif

    endsubroutine update_on_gpu
!**************************************************************************
    subroutine calcQ_gpu(idir, dir, stop, unit_vec, lperiodic)

      integer :: idir
      integer, dimension(3) :: dir, stop
      real, dimension(3) :: unit_vec
      logical :: lperiodic

      call calcQ_gpu_c(idir, dir, stop, unit_vec, lperiodic)
 
    endsubroutine calcQ_gpu
!**************************************************************************
    subroutine source_function_and_opacity_gpu(inu)
      integer :: inu
      call source_function_and_opacity_gpu_c(inu)
    endsubroutine
!**************************************************************************
    subroutine get_gpu_reduced_vars(dst)
      real, dimension(10) :: dst
      call get_gpu_reduced_vars_c(dst)
    endsubroutine get_gpu_reduced_vars
!**************************************************************************
    subroutine test_gpu_bcs
            call test_bcs_c
    endsubroutine test_gpu_bcs
!**************************************************************************
endmodule GPU
