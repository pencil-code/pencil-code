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
if(AC_lfrozen_bcs_x__mod__cdata)
{
	if(AC_iux__mod__cdata != 0 && AC_lfrozen_bot_var_x__mod__cdata[AC_iux__mod__cdata-1])
	{
		DF_UX = (AC_lfirst_proc_x__mod__cdata) ? 0.0 : DF_UX;
	}
	if(AC_iuy__mod__cdata != 0 && AC_lfrozen_bot_var_x__mod__cdata[AC_iuy__mod__cdata-1])
	{
		DF_UY = (AC_lfirst_proc_x__mod__cdata) ? 0.0 : DF_UY;
	}
	if(AC_iuz__mod__cdata != 0 && AC_lfrozen_bot_var_x__mod__cdata[AC_iuz__mod__cdata-1])
	{
		DF_UZ = (AC_lfirst_proc_x__mod__cdata) ? 0.0 : DF_UZ;
	}
	if(AC_ilnrho__mod__cdata != 0 && AC_lfrozen_bot_var_x__mod__cdata[AC_ilnrho__mod__cdata-1])
	{
		DF_RHO = (AC_lfirst_proc_x__mod__cdata) ? 0.0 : DF_RHO;
	}

	if(AC_iux__mod__cdata != 0 && AC_lfrozen_top_var_x__mod__cdata[AC_iux__mod__cdata-1])
	{
		DF_UX = (AC_llast_proc_x__mod__cdata) ? 0.0 : DF_UX;
	}
	if(AC_iuy__mod__cdata != 0 && AC_lfrozen_top_var_x__mod__cdata[AC_iuy__mod__cdata-1])
	{
		DF_UY = (AC_llast_proc_x__mod__cdata) ? 0.0 : DF_UY;
	}
	if(AC_iuz__mod__cdata != 0 && AC_lfrozen_top_var_x__mod__cdata[AC_iuz__mod__cdata-1])
	{
		DF_UZ = (AC_llast_proc_x__mod__cdata) ? 0.0 : DF_UZ;
	}
	if(AC_ilnrho__mod__cdata != 0 && AC_lfrozen_top_var_x__mod__cdata[AC_ilnrho__mod__cdata-1])
	{
		DF_RHO = (AC_llast_proc_x__mod__cdata) ? 0.0 : DF_RHO;
	}
}
