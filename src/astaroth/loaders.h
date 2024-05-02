
#if SINGLEPASS
  auto single_loader0= [](ParamLoadingInfo p)
  {
	  p.params -> singlepass_solve.step_num = 0;
	  p.params -> singlepass_solve.dt = p.device->local_config.real_params[AC_dt];
  };
  auto single_loader1= [](ParamLoadingInfo p)
  {
	  p.params -> singlepass_solve.step_num = 1;
	  p.params -> singlepass_solve.dt = p.device->local_config.real_params[AC_dt];
  };
  auto single_loader2= [](ParamLoadingInfo p)
  {
	  p.params -> singlepass_solve.step_num = 2;
	  p.params -> singlepass_solve.dt = p.device->local_config.real_params[AC_dt];
  };
#else
  auto intermediate_loader_0= [](ParamLoadingInfo p)
  {
	  p.params -> twopass_solve_intermediate.step_num = 0;
	  p.params -> twopass_solve_intermediate.dt = p.device->local_config.real_params[AC_dt];
  };
  auto final_loader_0 = [](ParamLoadingInfo p)
  {
	  p.params -> twopass_solve_final.step_num = 0;
  };
  auto intermediate_loader_1= [](ParamLoadingInfo p)
  {
	  p.params -> twopass_solve_intermediate.step_num = 1;
	  p.params -> twopass_solve_intermediate.dt = p.device->local_config.real_params[AC_dt];
  };
  auto final_loader_1 = [](ParamLoadingInfo p)
  {
	  p.params -> twopass_solve_final.step_num = 1;
  };
  auto intermediate_loader_2= [](ParamLoadingInfo p)
  {
	  p.params -> twopass_solve_intermediate.step_num = 2;
	  p.params -> twopass_solve_intermediate.dt = p.device->local_config.real_params[AC_dt];
  };
  auto final_loader_2= [](ParamLoadingInfo p)
  {
	  p.params -> twopass_solve_final.step_num = 2;
  };
#endif
