#if LHYDRO

Profile<XZ> AC_uu_average_cyl__mod__hydro 
Profile<XY> AC_uu_average_sph__mod__hydro

Kernel hydro_before_boundary(PC_SUB_STEP_NUMBER step_num)
{
	if((AC_lfargo_advection__mod__cdata && step_num == 0) || AC_lcalc_uuavg__mod__hydro)
	{
		if(AC_lcylindrical_coords__mod__cdata)
		{
	          reduce_sum(AC_inv_ngrid.y*UUY,AC_uu_average_cyl__mod__hydro)
		}
		else if(AC_lspherical_coords__mod__cdata)
		{
		  reduce_sum(AC_inv_ngrid.z*UUZ,AC_uu_average_sph__mod__hydro)
		}
	}
}
#else
Kernel hydro_before_boundary(PC_SUB_STEP_NUMBER step_num)
{
	suppress_unused_warning(step_num)
}
#endif
