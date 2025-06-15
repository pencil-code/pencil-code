#if LSELFGRAVITY
#include "../stdlib/poisson.h"
Kernel selfgravity_calc_rhs()
{
	rhs_poisson = 0.0
	if(ldensity && AC_lselfgravity_gas__mod__selfgravity)
	{
		if(AC_lstratz__mod__cdata)
		{
			rhs_poisson += AC_rho0z__mod__selfgravity[vertexIdx.z] * (1.0 + RHO)
		}
		else if(AC_ldensity_nolog__mod__cdata)
		{
			rhs_poisson += value(RHO)
		}
		else
		{
			rhs_poisson += exp(LNRHO)
		}
	}
	if(ldustdensity && AC_lselfgravity_dust__mod__selfgravity)
	{
		if(AC_ldustdensity_log__mod__cdata)
		{
			rhs_poisson += exp(Field(AC_ind__mod__cdata[0]))
		}
		else
		{
			rhs_poisson += value(Field(AC_ind__mod__cdata[0]))
		}
	}
	if(lneutraldensity && AC_lselfgravity_neutrals__mod__selfgravity)
	{
		if(AC_lneutraldensity_nolog__mod__cdata)
		{
			rhs_poisson += value(RHO_NEUTRAL)
		}
		else
		{
			rhs_poisson += exp(RHO_NEUTRAL)
		}
	}
	write(RHS_POISSON,rhs_poisson)
}
Kernel calc_final_potential(real t)
{
	if(AC_tselfgrav_gentle__mod__selfgravity > 0.0 && t < AC_tstart_selfgrav__mod__selfgravity + AC_tselfgrav_gentle__mod__selfgravity)
	{
		write(
			SELFGRAVITY_POTENTIAL, 
		        0.5*AC_rhs_poisson_const__mod__selfgravity*(1.0-cos(AC_REAL_PI * (t - AC_tstart_selfgrav__mod__selfgravity) / AC_tselfgrav_gentle__mod__selfgravity)) * SELFGRAVITY_POTENTIAL
		)
	}
	else
	{
		write(
			SELFGRAVITY_POTENTIAL,
			AC_rhs_poisson_const__mod__selfgravity*SELFGRAVITY_POTENTIAL
		)
	}
}
Kernel selfgravity_sor_step(int color)
{
	poisson_sor_red_black(color,RHS_POISSON,SELFGRAVITY_POTENTIAL,1.0)		
}
#else
Kernel selfgravity_calc_rhs(){}
Kernel calc_final_potential(real t){suppress_unused_warning(t)}
Kernel selfgravity_sor_step(int color){suppress_unused_warning(color)}
#endif

