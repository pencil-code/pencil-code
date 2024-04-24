  call omp_set_max_active_levels(2)
  num_helper_threads = omp_get_max_threads()-1
  lmultithread=.true.
  call signal_init

