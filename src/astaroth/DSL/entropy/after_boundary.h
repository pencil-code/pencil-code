#if Lentropy_MODULE
Kernel entropy_reductions()
{
	if(AC_lcalc_ssmean__mod__energy)
	{
		reduce_sum(SS/nxygrid,AC_ssmz__mod__energy)
	}
	if(AC_lcalc_ssmeanxy__mod__energy)
	{
		reduce_sum(SS/nyzgrid,AC_ssmx__mod__energy)
	}
}

Kernel entropy_update_running_average(real t, real dt)
{
  if(AC_lss_running_aver__mod__energy)
  {
	real fss_run_aver = value(F_SS_RUN_AVER)
  	if(t == 0. || t<dt)
  	{
		fss_run_aver = value(F_SS)
  	}
  	frac = dt/AC_tau_aver1__mod__cdata
	F_SS_RUN_AVER[vertexIdx.x][vertexIdx.y][vertexIdx.z] = (1.0-frac)*fss_run_aver + frac*F_SS
  }
}

Kernel entropy_mean_derivs_z()
{
	if(AC_lcalc_ssmean__mod__energy)
	{
		AC_gssmz__mod__energy[vertexIdx.z][0] = 0.0
		AC_gssmz__mod__energy[vertexIdx.z][1] = 0.0
		AC_gssmz__mod__energy[vertexIdx.z][2] = derz(AC_ssmz__mod__energy)
		write(AC_del2ssmz__mod__energy,derzz(AC_ssmz__mod__energy))
	}
}
Kernel entropy_mean_derivs_x()
{
	if(AC_lcalc_ssmeanxy__mod__energy)
	{
	  derssmx = derx(AC_ssmx__mod__energy)
	  AC_gssmx__mod__energy[vertexIdx.x-NGHOST][0] = derssmx
	  AC_gssmx__mod__energy[vertexIdx.x-NGHOST][1] = 0.0
	  AC_gssmx__mod__energy[vertexIdx.x-NGHOST][2] = 0.0
	  del2ssmx = derxx(AC_ssmx__mod__energy)
	  if(AC_lspherical_coords__mod__cdata)
	  {
		  del2ssmx += 2.*derssmx/AC_x__mod__cdata[vertexIdx.x]
	  }
	  write(AC_del2ssmx__mod__energy,del2ssmx)
	}
}
Kernel
entropy_smooth(bool lrmv)
{
  if(AC_lss_running_aver__mod__energy && AC_lsmooth_ss_run_aver__mod__energy && lrmv)
  {
  	write(F_SS_RUN_AVER,binomial_smooth(F_SS_RUN_AVER))
  }
}
#else
Kernel entropy_reductions()
{
}
Kernel entropy_update_running_average(real t, real dt)
{
	suppress_unused_warning(t)
	suppress_unused_warning(dt)
}
Kernel entropy_mean_derivs_x()
{
}
Kernel entropy_mean_derivs_z()
{
}
Kernel entropy_smooth(bool lrmv)
{
	suppress_unused_warning(lrmv)
}
#endif
