  public :: register_GPU, initialize_GPU, finalize_GPU, get_farray_ptr_gpu, rhs_GPU, &
            copy_farray_from_GPU, &
            read_gpu_run_pars, write_gpu_run_pars, &
            load_farray_to_GPU, reload_GPU_config, update_on_gpu, get_ptr_GPU, get_ptr_GPU_training, &
            before_boundary_gpu, &
            update_after_substep_gpu, &
            gpu_set_dt, train_gpu, infer_gpu,radtransfer_gpu, &
            get_gpu_reduced_vars,test_gpu_bcs, split_update_gpu, &
            pushpars2c,ltest_bcs

  private
