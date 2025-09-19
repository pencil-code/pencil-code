#if LDENSITY
global output real AC_current_total_mass

Kernel get_current_total_mass(bool lrmv)
{
	if (lrmv && AC_lconserve_total_mass__mod__density && AC_total_mass__mod__density > 0.0)
	{
		rho = (AC_ldensity_nolog__mod__cdata) ? value(RHO) : exp(LNRHO);
		integration_weight =  AC_dvol_x__mod__cdata[vertexIdx.x]
			            * AC_dvol_y__mod__cdata[vertexIdx.y]
			            * AC_dvol_z__mod__cdata[vertexIdx.z]
		reduce_sum(rho*integration_weight,AC_current_total_mass)
	}
}
//Using the array syntax on the lhs signals to Astaroth that the writes are done inplace and no swapping is thus need
//This is important since otherwise the buffers will be swapped incorrectly for substeps > 1 (with normal rk3 timestepping)
Kernel fix_mass_drift(bool lrmv)
{
    if (lrmv && AC_lconserve_total_mass__mod__density && AC_total_mass__mod__density > 0.0)
    {
     	real fact=AC_total_mass__mod__density/AC_current_total_mass
     	if (AC_ldensity_nolog__mod__cdata)
     	{
		RHO[vertexIdx.x][vertexIdx.y][vertexIdx.z] = fact*RHO
     	}
     	else
     	{
		LNRHO[vertexIdx.x][vertexIdx.y][vertexIdx.z] = LNRHO + log(fact)
     	}
#if LHYDRO
        UU[vertexIdx.x][vertexIdx.y][vertexIdx.z] = UU/fact
#endif
    }
}
#else
Kernel get_current_total_mass(bool lrmv)
{
}
Kernel fix_mass_drift(bool lrmv)
{
}
#endif
