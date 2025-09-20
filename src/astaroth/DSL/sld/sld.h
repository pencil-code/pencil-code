Kernel sld_calc_char_speed(PC_SUB_STEP_NUMBER step_num)
{
#if LMAGNETIC
#if LEOS
#if LHYDRO
	if(AC_lslope_limit_diff__mod__cdata && step_num == AC_num_substeps__mod__cdata-1)
	{
		res = calculate_characteristic_speed(AC_w_sldchar_hyd__mod__hydro, UU, AC_w_sldchar_ene__mod__energy, AC_cs0__mod__equationofstate, AC_w_sldchar_mag__mod__magnetic, curl(AA), LNRHO, AC_mu0__mod__cdata)
		write(SLD_CHAR_SPEED,res)
	}
#endif
#endif
#endif
}
