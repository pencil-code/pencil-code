#if LDENSITY
Kernel get_current_total_mass(bool lrmv)
{
	if(lrmv && lconserve_total_mass)
	{
		rho = (ldensity_nolog) ? value(RHO) : exp(LNRHO);
		integration_weight = AC_ds.x*AC_ds.y*AC_ds.z;
		reduce_sum(rho*integration_weight,AC_current_total_mass)
	}
}
Kernel fix_mass_drift(bool lrmv)
{

    if(lrmv && lconserve_total_mass)
    {
     	real fact=total_mass/AC_current_total_mass
     	if(ldensity_nolog)
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
Kernel get_current_total_mass()
{
}
Kernel fix_mass_drift()
{
}
#endif
