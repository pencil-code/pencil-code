if(AC_lfreeze_var_all__mod__cdata)
{
	if(AC_wfreeze_int__mod__cdata != 0.0)
	{
                pfreeze=quintic_step(ac_transformed_pencil_r_mn,AC_rfreeze_int__mod__cdata,AC_wfreeze_int__mod__cdata,AC_fshift_int__mod__cdata)
		DF_UVEC *= pfreeze
		DF_AVEC *= pfreeze
		DF_RHO  *= pfreeze
		DF_SS   *= pfreeze
	}
}
