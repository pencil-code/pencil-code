#if LHYDRO
#if Lhydro_MODULE

Profile<XZ> AC_uu_average_cyl__mod__hydro 
Profile<XY> AC_uu_average_sph__mod__hydro
global output real AC_angmom_global
global output real AC_rhosint_global

Kernel hydro_before_boundary_reductions(bool lrmv)
{
	if(AC_lremove_mean_angmom__mod__hydro && lrmv)
	{
      		fac = 1./(AC_lxyz__mod__cdata.x*(cos(AC_y0__mod__cdata)-cos(AC_y0__mod__cdata+AC_lxyz__mod__cdata.y))*AC_lxyz__mod__cdata.z)
		rho = AC_ldensity_nolog__mod__cdata ? value(RHO) : exp(LNRHO)
		wx = AC_x__mod__cdata[vertexIdx.x]*AC_dvol_x__mod__cdata[vertexIdx.x]
		tmp = rho*wx	
		wmn = AC_sinth__mod__cdata[vertexIdx.y]*AC_dvol_y__mod__cdata[vertexIdx.y]*AC_dvol_z__mod__cdata[vertexIdx.z]
		angmom  = tmp*UU.z*wmn
		rhosint = tmp*wmn

		reduce_sum(fac*angmom,AC_angmom_global)
		reduce_sum(fac*rhosint,AC_rhosint_global)
	}
}

Kernel hydro_before_boundary(bool lrmv, PC_SUB_STEP_NUMBER step_num)
{
	if((AC_lfargo_advection__mod__cdata && step_num == 0) || AC_lcalc_uuavg__mod__hydro)
	{
		if(AC_lcylindrical_coords__mod__cdata)
		{
	          reduce_sum(AC_ngrid_inv.y*UUY,AC_uu_average_cyl__mod__hydro)
		}
		else if(AC_lspherical_coords__mod__cdata)
		{
		  reduce_sum(AC_ngrid_inv.z*UUZ,AC_uu_average_sph__mod__hydro)
		}
	}
	if(AC_lremove_mean_angmom__mod__hydro && lrmv)
	{
		um = AC_angmom_global/AC_rhosint_global
		UUZ[vertexIdx.x][vertexIdx.y][vertexIdx.z] = UU.z-um
	}
}

#else
Kernel hydro_before_boundary_reductions(bool lrmv)
{
	suppress_unused_warning(lrmv)
}
Kernel hydro_before_boundary(bool lrmv, PC_SUB_STEP_NUMBER step_num)
{
	suppress_unused_warning(lrmv)
	suppress_unused_warning(step_num)
}
#endif
#else
Kernel hydro_before_boundary(bool lrmv, PC_SUB_STEP_NUMBER step_num)
{
	suppress_unused_warning(lrmv)
	suppress_unused_warning(step_num)
}
#endif
