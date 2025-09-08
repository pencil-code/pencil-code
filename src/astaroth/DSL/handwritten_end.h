 if(AC_lcourant_dt__mod__cdata)
 {
  	if (AC_iuu__mod__cdata != 0)    write( F_UU,  rk_intermediate(F_UU , DF_UU,  step_num, AC_dt__mod__cdata) )
  	if ((AC_ilnrho__mod__cdata + AC_irho__mod__cdata) != 0)  write(F_RHO, rk_intermediate(F_RHO, DF_RHO, step_num, AC_dt__mod__cdata) )
  	if (AC_iss__mod__cdata != 0)  write( F_SS,  rk_intermediate(F_SS, DF_SS,  step_num, AC_dt__mod__cdata) )
  	if (AC_iaa__mod__cdata != 0) write( F_AA,  rk_intermediate(F_AA , DF_AA,  step_num, AC_dt__mod__cdata) )
#if LGRAVITATIONAL_WAVES_HTXK
	if (AC_lfirst__mod__cdata)
	{
		write(F_STRESS_0,DF_STRESS_0)
		write(F_STRESS_1,DF_STRESS_1)
		write(F_STRESS_2,DF_STRESS_2)
		write(F_STRESS_3,DF_STRESS_3)
		write(F_STRESS_4,DF_STRESS_4)
		write(F_STRESS_5,DF_STRESS_5)
	}
	if (AC_lfirst__mod__cdata)
	{
		reduce_max(dt1_max_loc,AC_dt1_max)
	}
#endif
 }
 else
 {
 	maximum_error = 0.0 
 	if (AC_iuu__mod__cdata != 0) maximum_error = rkf4_update(DF_UU,step_num,AC_dt__mod__cdata,ERROR_UU,BETA_UU,F_UU,uux_initial_max,uuy_initial_max,uuz_initial_max,maximum_error,AC_dt_ratio__mod__cdata,AC_dt_epsi__mod__cdata)
 	if (AC_iaa__mod__cdata != 0) maximum_error = rkf4_update(DF_AA,step_num,AC_dt__mod__cdata,ERROR_AA,BETA_AA,F_AA,aax_initial_max,aay_initial_max,aaz_initial_max,maximum_error,AC_dt_ratio__mod__cdata,AC_dt_epsi__mod__cdata)
 	if ((AC_ilnrho__mod__cdata + AC_irho__mod__cdata) != 0) maximum_error = rkf4_update(DF_RHO,step_num,AC_dt__mod__cdata,ERROR_RHO,BETA_RHO,F_RHO,rho_initial_max,maximum_error,AC_dt_ratio__mod__cdata,AC_dt_epsi__mod__cdata)
 	if (AC_iss__mod__cdata != 0) maximum_error = rkf4_update(DF_SS,step_num,AC_dt__mod__cdata,ERROR_SS,BETA_SS,F_SS,ss_initial_max,maximum_error,AC_dt_ratio__mod__cdata,AC_dt_epsi__mod__cdata)
 	if(step_num == 4) reduce_max(maximum_error,AC_maximum_error)
 }

