  public :: register_GPU, initialize_GPU, finalize_GPU, get_farray_ptr_gpu, rhs_GPU, &
            copy_farray_from_GPU, &
	    read_gpu_run_pars, write_gpu_run_pars, &
            test_rhs_gpu, load_farray_to_GPU, reload_GPU_config, update_on_gpu, get_ptr_GPU, get_ptr_GPU_training, &
            calcQ_gpu, before_boundary_gpu, &
	    gpu_set_dt, train_gpu, infer_gpu,source_function_and_opacity_gpu

  private
