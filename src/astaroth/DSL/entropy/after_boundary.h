#if Lentropy_MODULE
Kernel entropy_reductions()
{
	if(AC_lcalc_ssmean__mod__energy)
	{
		reduce_sum(SS/nxygrid,AC_ssmz__mod__energy)
	}
}
Kernel entropy_mean_derivs()
{
	if(AC_lcalc_ssmean__mod__energy)
	{
		AC_gssmz__mod__energy[vertexIdx.z][0] = 0.0
		AC_gssmz__mod__energy[vertexIdx.z][1] = 0.0
		AC_gssmz__mod__energy[vertexIdx.z][2] = derz(AC_ssmz__mod__energy)
		write(AC_del2ssmz__mod__energy,derzz(AC_ssmz__mod__energy))
	}
}
#else
Kernel entropy_reductions()
{
}
Kernel entropy_mean_derivs()
{
}
#endif
