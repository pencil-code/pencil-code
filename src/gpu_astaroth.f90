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
  use General, only: lpointer, ioptest, loptest
  use Quiet
  use Messages
!$ use, intrinsic :: iso_c_binding
  use iso_c_binding

  implicit none

  include 'gpu.h'
  
  external initialize_gpu_c
  external register_gpu_c
  external finalize_gpu_c
  external get_farray_ptr_gpu_c
  external rhs_gpu_c
  external before_boundary_gpu_c
  external update_after_substep_gpu_c 
  external radtransfer_gpu_c
  external load_farray_c
  external reload_gpu_config_c
  external copy_farray_c
  external update_on_gpu_arr_by_ind_c
  external update_on_gpu_scal_by_ind_c
  external update_on_gpu_vec_by_ind_c
  external pos_real_ptr_c
  external gpu_prepare_for_first_substep_c
  external get_gpu_reduced_vars_c
  external test_bcs_c
  external split_update_gpu_c

  ! Torchfort
  external tf_train_c 
  external tf_infer_c
  external tf_create_model_c
  external print_snapshot_c
  external tf_load_model_c
  external tf_load_model_checkpoint_c
  external tf_save_model_c
  external tf_save_checkpoint_c

  integer, external :: update_on_gpu_arr_by_name_c
  integer, external :: update_on_gpu_scal_by_name_c
  integer, external :: update_on_gpu_vec_by_name_c

  type(C_PTR) :: pFarr_GPU_in, pFarr_GPU_out

  logical :: lcpu_timestep_on_gpu=.false.
  !PAR-DOC: Since on the GPU, calculation of df and update of f happen without synchronization
  !PAR-DOC: of the decomposed portions of the subdomains (or different processes) we can't compute
  !PAR-DOC: dt 'on the fly' but take it from the previous timestep (calculated at the moment on
  !PAR-DOC: the first substep but the last would make more sense). 
  !PAR-DOC: lcpu_timestep_on_gpu=T means to calculate the rhs an extra time
  !PAR-DOC: on the first substep to compute dt as is done normally on the CPUs

  logical :: lsingle_precision_timestep =.true.
  !PAR-DOC: Whether we compute the timestep in single precision regardless of the precision of real.
  !PAR-DOC: We are a bit brave and do it by default in single since the exact value of the timestep
  !PAR-DOC: Should never be that important: safety factors and other controls are surely more important
  !PAR-DOC: than the negligible amount of round-off.
  !PAR-DOC: Anyways one is supposed to be able to turn it off for testing.
  !MR: Why is single for dt relevant?

  logical :: lac_sparse_autotuning=.true.
  !PAR-DOC: Astaroth empirically finds the best configuration parameters of the kernels by running them with different 
  !PAR-DOC: options and picking the best one. Sparse autotuning means to prune the parameter space more than 
  !PAR-DOC: usual by picking only the most likely ones (empirically gives for large grids the same as the 
  !PAR-DOC: larger search, but is considerably faster).

  logical :: lac_sparse_autotuning_always=.false.

  logical :: lcumulative_df_on_gpu=.false.
  !PAR-DOC: Whether df is 0 or the accumulated value at the start of the kernel.  MR: don't understand this: which kernel?
  !PAR-DOC: For simplicity df is zero and then accumulated to the buffer, but sometimes you need to be careful like with
  !PAR-DOC: short stopping time approximation in dustvelocity (the only case we are aware at the moment where the difference matters).

  logical :: lskip_rtime_compilation=.false.    !PAR-DOC: Placeholder

  logical :: lread_all_vars_from_device = .false.
  !PAR-DOC: By default only pde variables and those aux variables that are registered to be always read are read from the device.
  !PAR-DOC: If this is true all (Field) variables are always read.

  logical :: lcuda_aware_mpi=.true.
  !PAR-DOC: Whether to use CUDA-aware MPI (faster). If you have it you should always want to use it, but sometimes you do not have it or using
  !PAR-DOC: it is more unstable than routing the communication via the host yourself.

  logical :: ltest_bcs =.false.    !PAR-DOC: Whether to test the agreement of bcs on GPU and CPU.

  logical :: ltest_rhs =.false.
  !PAR_DOC: Whether to test the agreement of RHS on GPU and CPU
  !PAR-DOC: Will run both versions and compare the differences.

  integer :: it_test_rhs = 1
  !PAR-DOC: At which timestep to perform the comparison. It is useful to be able to vary this in case some samples need some
  !PAR-DOC: integration time for some of the fields to have meaningful values.

  integer, dimension(3) :: thread_block_loop_factors = (/1,1,1/)
  !PAR-DOC: How should reduction kernels decompose the subdomain into smaller chunks.
  !PAR-DOC: The decomposition allows to do partial reductions before the global reduction. Saves memory and work.

  logical :: lonly_default_stream_for_taskgraphs = .false.
  !PAR-DOC: Saves memory at the potential cost of performance by not creating streams,
  !PAR-DOC: which take surprisingly large amount of memory. Mainly needed for the autotest in norlx51

  namelist /gpu_run_pars/ &
     ltest_bcs,lac_sparse_autotuning,lac_sparse_autotuning_always,lcuda_aware_mpi, &
     lcpu_timestep_on_gpu,lsingle_precision_timestep,lcumulative_df_on_gpu, &
     lread_all_vars_from_device,ltest_rhs,it_test_rhs,thread_block_loop_factors,lonly_default_stream_for_taskgraphs

contains
!***********************************************************************
    subroutine TF_train_gpu(model_name, loss, itsub, t)
  
      real :: loss
      integer :: itsub 
      real(KIND=rkind8), intent(IN) :: t
      character(len=*), intent(in) :: model_name
  
      call tf_train_c(trim(model_name)//c_null_char, loss, itsub, t)
  
    endsubroutine TF_train_gpu 
!***********************************************************************
    subroutine TF_infer_gpu(model_name, flag)
  
      integer :: flag
      character(len=*), intent(in) :: model_name

      call tf_infer_c(trim(model_name)//c_null_char, flag)
  
    endsubroutine TF_infer_gpu
!***********************************************************************
    subroutine TF_create_model(model_name, config_file_path, lmpicomm)
  
      use Mpicomm, only: MPI_COMM_PENCIL
  
      integer :: lmpicomm_int
      logical :: lmpicomm
      character(len=*), intent(in) :: model_name, config_file_path
  
      lmpicomm_int = merge(1,0,lmpicomm)
      call tf_create_model_c(trim(model_name)//c_null_char, trim(config_file_path)//c_null_char, MPI_COMM_PENCIL, lmpicomm_int)
  
    endsubroutine TF_create_model
!***********************************************************************
    subroutine tau_snapshots()
  
      call print_snapshot_c()
  
    endsubroutine tau_snapshots
!***********************************************************************
    subroutine TF_load_model(model_name, fname)
  
      character(len=*), intent(in) :: model_name, fname
  
      call tf_load_model_c(trim(model_name)//c_null_char, trim(fname)//c_null_char)
  
    endsubroutine TF_load_model
!***********************************************************************
    subroutine TF_load_model_checkpoint(model_name, checkpoint_dir)
  
      character(len=*), intent(in) :: model_name, checkpoint_dir
  
      call tf_load_model_checkpoint_c(trim(model_name)//c_null_char, trim(checkpoint_dir)//c_null_char)
  
    endsubroutine TF_load_model_checkpoint
!***********************************************************************
    subroutine TF_save_model(model_name, fname)
  
      character(len=*), intent(in) :: model_name, fname
  
      call tf_save_model_c(trim(model_name)//c_null_char, trim(fname)//c_null_char)
  
    endsubroutine TF_save_model
!***********************************************************************
    subroutine TF_save_checkpoint(model_name, checkpoint_dir)
  
      character(len=*), intent(in) :: model_name, checkpoint_dir
  
      call tf_save_checkpoint_c(trim(model_name)//c_null_char, trim(checkpoint_dir)//c_null_char)
  
    endsubroutine TF_save_checkpoint
!***********************************************************************
    subroutine initialize_GPU(f)
!
      use Mpicomm, only: MPI_COMM_PENCIL

      real, dimension(:,:,:,:), intent(IN) :: f
      integer :: lread_all_vars_from_device_int
      integer :: lcpu_timestep_on_gpu_int
      integer :: lac_sparse_autotuning_int

      character(LEN=512) :: str
!
      if (ltest_rhs) lread_all_vars_from_device = .true.
      !If there are enough GPUs we can distribute the autotuning between them
      !and it won't take too long
      if (ncpus >= 8) lac_sparse_autotuning = .false.

      if (lac_sparse_autotuning_always) lac_sparse_autotuning = .true.

      str=''
      !List of unsupported modules
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
      if (dt<=0.) dt = dtmin

      lread_all_vars_from_device_int = merge(1,0,lread_all_vars_from_device)
      lcpu_timestep_on_gpu_int       = merge(1,0,lcpu_timestep_on_gpu)
      lac_sparse_autotuning_int      = merge(1,0,lac_sparse_autotuning)
      call initialize_gpu_c(f,MPI_COMM_PENCIL,t,nt,lread_all_vars_from_device_int,&
                            lcpu_timestep_on_gpu_int,lac_sparse_autotuning_int)
!
! Load farray to gpu
!
      if (nt>0) call load_farray_to_GPU(f)

      call get_farray_ptr_gpu
!print'(a,1x,Z0,1x,Z0)', 'pFarr_GPU_in,pFarr_GPU_out=', pFarr_GPU_in,pFarr_GPU_out
!flush(6)

    endsubroutine initialize_GPU
!**************************************************************************
    subroutine read_gpu_run_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(LEN=*), intent(out) :: iomsg
      integer :: iostat
!
      read(parallel_unit, NML=gpu_run_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
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

      call get_farray_ptr_gpu_c(pFarr_GPU_in,pFarr_GPU_out)

    endsubroutine get_farray_ptr_gpu
!**************************************************************************
    subroutine rhs_GPU(f,isubstep)
!
      use General, only: notanumber

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f
      integer,                            intent(IN)    :: isubstep
!
      call keep_compiler_quiet(f)
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
      integer :: lrmv_int,lsubstepping_in_time_int
!
      call keep_compiler_quiet(f)
      !transform to integers since they are more compatible with C than logical to booleans
      lrmv_int = merge(1,0,lrmv)
      lsubstepping_in_time_int = merge(1,0,lsubstepping_in_time)
      call before_boundary_gpu_c(lrmv_int,isubstep,t,lsubstepping_in_time_int)
!
    endsubroutine before_boundary_gpu
!**************************************************************************
    subroutine update_after_substep_gpu

      call update_after_substep_gpu_c

    endsubroutine update_after_substep_gpu
!**************************************************************************
    subroutine gpu_prepare_for_first_substep
!
      call gpu_prepare_for_first_substep_c(t)
!
    endsubroutine gpu_prepare_for_first_substep 
!**************************************************************************
    function get_ptr_GPU(ind1,ind2,lout) result(pFarr)
!
!  Fetches the address of the f-array counterpart on the GPU for slots from ind1 to ind2
!  and transforms it to a Fortran pointer.
!
      integer, optional :: ind1,ind2
      logical, optional :: lout

      real, dimension(:,:,:,:), pointer :: pFarr

      integer :: i1,i2

      interface
        type(c_ptr) function pos_real_ptr_c(ptr,ind)
          import :: c_ptr, ikind8
          type(c_ptr) :: ptr
          integer :: ind
        endfunction
      endinterface

      i1 = ioptest(ind1)
      i2 = ioptest(ind2,ind1)

      if (loptest(lout)) then
        call c_f_pointer(pos_real_ptr_c(pFarr_GPU_out,i1-1),pFarr,(/mx,my,mz,i2-i1+1/))
      else
        call c_f_pointer(pos_real_ptr_c(pFarr_GPU_in,i1-1),pFarr,(/mx,my,mz,i2-i1+1/))
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

      call keep_compiler_quiet(f)
      call load_farray_c

    endsubroutine load_farray_to_GPU
!**************************************************************************
    subroutine reload_GPU_config

      call reload_gpu_config_c

    endsubroutine reload_GPU_config
!**************************************************************************
    subroutine update_on_gpu_vec(index, varname, value)
!
!  Updates a 3-vector element of the Astaroth configuration, identified by name or index, on the GPU.
!
      integer, intent(inout) :: index
      character(LEN=*),optional :: varname
      real, dimension(3), optional :: value

      if (index>=0) then
        if (present(value)) then
          call update_on_gpu_vec_by_ind_c(index,value)
        endif
      else
        if (present(value)) then
          index = update_on_gpu_vec_by_name_c(varname//char(0),value)
        endif
        if (index<0) call fatal_error('update_on_gpu','variable "'//trim(varname)//'" not found')
      endif

    endsubroutine update_on_gpu_vec
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
    subroutine radtransfer_gpu

      call radtransfer_gpu_c

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
    subroutine split_update_gpu(f)

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f

      call split_update_gpu_c(f)

    endsubroutine split_update_gpu
!**************************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=20
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(lskip_rtime_compilation,p_par(3)) ! bool
    call copy_addr(lcumulative_df_on_gpu,p_par(4)) ! bool
    call copy_addr(ltest_bcs,p_par(7)) ! bool
    call copy_addr(lsingle_precision_timestep,p_par(8)) ! bool
    call copy_addr(thread_block_loop_factors,p_par(9)) ! int3 dconst
    call copy_addr(lonly_default_stream_for_taskgraphs,p_par(10)) ! bool dconst
    call copy_addr(lcuda_aware_mpi,p_par(11)) ! bool dconst

    endsubroutine pushpars2c
!**************************************************************************
endmodule GPU
