if (lstart) then
  num_helper_threads = omp_get_max_threads()
else
  call omp_set_max_active_levels(2)
  num_helper_threads = omp_get_max_threads()-1
  if (num_helper_threads==0) call fatal_error('run','zero helper threads in multithreaded version')
  lmultithread=.true.
  call signal_init
  loffload = omp_get_num_devices() /= 0
endif
