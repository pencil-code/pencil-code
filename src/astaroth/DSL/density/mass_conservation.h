#if LDENSITY
global output real AC_current_total_mass

Kernel get_current_total_mass(bool lrmv)
{
	if(lrmv && AC_lconserve_total_mass__mod__density)
	{
		rho = (AC_ldensity_nolog__mod__cdata) ? value(RHO) : exp(LNRHO);
		integration_weight = AC_ds.x*AC_ds.y*AC_ds.z;
		reduce_sum(rho*integration_weight,AC_current_total_mass)
	}
}
Kernel fix_mass_drift(bool lrmv)
{

    if(lrmv && AC_lconserve_total_mass__mod__density)
    {
     	real fact=AC_total_mass__mod__density/AC_current_total_mass
     	if(AC_ldensity_nolog__mod__cdata)
     	{
     	        write(RHO,fact*RHO)
     	}
     	else
     	{
     	        write(LNRHO,LNRHO + log(fact))
     	}
#if LHYDRO
     	write(UU, UU/fact)
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
