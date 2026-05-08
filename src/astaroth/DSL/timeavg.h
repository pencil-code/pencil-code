#if LTIMEAVG
Field TIME_AVERAGES[AC_mtavg__mod__timeavg]
Kernel update_timeavgs_kernel(int timestep, real dt)
{
	weigth = min(dt/AC_tavg__mod__timeavg,1.0)
	for i in 0:mtavg
	{
		Field f = Field(AC_idx_tavg__mod_timeavg[i])
		if(timestep == 0)
		{
			write(TIME_AVERAGES[i],f)
		}
		else
		{
			write(TIME_AVERAGES[i], TIME_AVERAGES[i] + weight*(f - TIME_AVERAGES[i]))
		}
	}
}
ComputeSteps
AC_update_timeavgs(boundconds)
{
	update_timeavgs_kernel(AC_timestep_number,AC_dt)
}
#endif
