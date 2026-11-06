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
if(AC_lfrozen_bcs_z__mod__cdata)
{
	const bool lower_boundary = AC_lfirst_proc_z__mod__cdata && vertexIdx.z == AC_n1-1
	const bool upper_boundary = AC_llast_proc_z__mod__cdata && vertexIdx.z == AC_n2-1
	if(AC_iux__mod__cdata != 0 && AC_lfrozen_bot_var_z__mod__cdata[max(AC_iux__mod__cdata-1,0)])
	{
		DF_UX = lower_boundary ? 0.0 : DF_UX;
	}
	if(AC_iux__mod__cdata != 0 && AC_lfrozen_top_var_z__mod__cdata[max(AC_iux__mod__cdata-1,0)])
	{
		DF_UX = upper_boundary ? 0.0 : DF_UX;
	}
	if(AC_iuy__mod__cdata != 0 && AC_lfrozen_bot_var_z__mod__cdata[max(AC_iuy__mod__cdata-1,0)])
	{
		DF_UY = lower_boundary ? 0.0 : DF_UY;
	}
	if(AC_iuy__mod__cdata != 0 && AC_lfrozen_top_var_z__mod__cdata[max(AC_iuy__mod__cdata-1,0)])
	{
		DF_UY = upper_boundary ? 0.0 : DF_UY;
	}
	if(AC_iuz__mod__cdata != 0 && AC_lfrozen_bot_var_z__mod__cdata[max(AC_iuz__mod__cdata-1,0)])
	{
		DF_UZ = lower_boundary ? 0.0 : DF_UZ;
	}
	if(AC_iuz__mod__cdata != 0 && AC_lfrozen_top_var_z__mod__cdata[max(AC_iuz__mod__cdata-1,0)])
	{
		DF_UZ = upper_boundary ? 0.0 : DF_UZ;
	}
	if(AC_ilnrho__mod__cdata != 0 && AC_lfrozen_bot_var_z__mod__cdata[max(AC_ilnrho__mod__cdata-1,0)])
	{
		DF_RHO = lower_boundary ? 0.0 : DF_RHO;
	}
	if(AC_ilnrho__mod__cdata != 0 && AC_lfrozen_top_var_z__mod__cdata[max(AC_ilnrho__mod__cdata-1,0)])
	{
		DF_RHO = upper_boundary ? 0.0 : DF_RHO;
	}
	if(AC_ilntt__mod__cdata != 0 && AC_lfrozen_bot_var_z__mod__cdata[max(AC_ilntt__mod__cdata-1,0)])
	{
		DF_TT = lower_boundary ? 0.0 : DF_TT;
	}
	if(AC_ilntt__mod__cdata != 0 && AC_lfrozen_top_var_z__mod__cdata[max(AC_ilntt__mod__cdata-1,0)])
	{
		DF_TT = upper_boundary ? 0.0 : DF_TT;
	}
	if(AC_iss__mod__cdata != 0 && AC_lfrozen_bot_var_z__mod__cdata[max(AC_iss__mod__cdata-1,0)])
	{
		DF_SS = lower_boundary ? 0.0 : DF_SS;
	}
	if(AC_iss__mod__cdata != 0 && AC_lfrozen_top_var_z__mod__cdata[max(AC_iss__mod__cdata-1,0)])
	{
		DF_SS = upper_boundary ? 0.0 : DF_SS;
	}
	if(AC_iax__mod__cdata != 0 && AC_lfrozen_bot_var_z__mod__cdata[max(AC_iax__mod__cdata-1,0)])
	{
		DF_AX = lower_boundary ? 0.0 : DF_AX;
	}
	if(AC_iax__mod__cdata != 0 && AC_lfrozen_top_var_z__mod__cdata[max(AC_iax__mod__cdata-1,0)])
	{
		DF_AX = upper_boundary ? 0.0 : DF_AX;
	}
	if(AC_iay__mod__cdata != 0 && AC_lfrozen_bot_var_z__mod__cdata[max(AC_iay__mod__cdata-1,0)])
	{
		DF_AY = lower_boundary ? 0.0 : DF_AY;
	}
	if(AC_iay__mod__cdata != 0 && AC_lfrozen_top_var_z__mod__cdata[max(AC_iay__mod__cdata-1,0)])
	{
		DF_AY = upper_boundary ? 0.0 : DF_AY;
	}
	if(AC_iaz__mod__cdata != 0 && AC_lfrozen_bot_var_z__mod__cdata[max(AC_iaz__mod__cdata-1,0)])
	{
		DF_AZ = lower_boundary ? 0.0 : DF_AZ;
	}
	if(AC_iaz__mod__cdata != 0 && AC_lfrozen_top_var_z__mod__cdata[max(AC_iaz__mod__cdata-1,0)])
	{
		DF_AZ = upper_boundary ? 0.0 : DF_AZ;
	}
}
