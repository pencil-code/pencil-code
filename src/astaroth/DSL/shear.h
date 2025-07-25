#if LSHEAR

if(lshear)
{
  if(!AC_lshearadvection_as_shift__mod__shear)
  {
        const real uy0 = AC_uy0__mod__shear[vertexIdx.x-NGHOST]
	DF_UU  -= uy0*dery(UU)
	DF_AA  -= uy0*dery(AA)
	DF_SS  -= uy0*dery(SS)
	DF_RHO -= uy0*dery(RHO)
	#if LDUSTDENSITY
	for k in 0:ndustspec
	{
	      DF_DUST_VELOCITY[k]  -= AC_sshear1__mod__shear*uy0*dery(F_DUST_VELOCITY[k])
	      DF_DUST_DENSITY[k]   -= AC_sshear1__mod__shear*uy0*dery(F_DUST_DENSITY[k])
	}
	#endif
  }
  if (lhydro && AC_lshear_acceleration__mod__shear) DF_UU.y  -= AC_sshear1__mod__shear * UU.x
  if(AC_lhyper3x_mesh__mod__shear)
  {
  	const real d = AC_diff_hyper3x_mesh__mod__shear*abs(AC_sshear__mod__cdata)
	DF_UU  += d*der6x_ignore_spacing(UU)
	DF_AA  += d*der6x_ignore_spacing(AA)
	DF_RHO += d*der6x_ignore_spacing(RHO)
	DF_SS  += d*der6x_ignore_spacing(SS)
	#if LDUSTDENSITY
	for k in 0:ndustspec
	{
	      DF_DUST_VELOCITY[k]  += d*der6x_ignore_spacing(F_DUST_VELOCITY[k])
	      DF_DUST_DENSITY[k]   += d*der6x_ignore_spacing(F_DUST_DENSITY[k])
	}
        #endif
	if(AC_lupdate_courant_dt__mod__cdata)
	{
		diffus_shear3 = d
		maxdiffus3__mod__cdata=max(maxdiffus3__mod__cdata,diffus_shear3)
	}
  }
  if (lmagnetic && !lbfield && AC_lmagnetic_stretching__mod__shear)
  {
	DF_AVEC.x      -= AC_sshear__mod__cdata*F_AVEC.y
  	if(AC_lmagnetic_tilt__mod__shear)
  	{
	  DF_AVEC.x      -= AC_sshear_sini__mod__shear*F_AVEC.x
	  DF_AVEC.y      += AC_sshear_sini__mod__shear*F_AVEC.y
  	}
  }
  #if LDUSTVELOCITY
  if(ldustvelocity)
  {
	  for k in 0:ndustspec
	  {
		DF_DUST_VELOCITY[k].y  -= AC_sshear1__mod__shear*DF_DUST_VELOCITY[k].x
	  }
  }
  #endif
  /**
  if(ltestfield)
  {
	int j = AC_iaatest__mod__cdata-1
	while(j < AC_iaztestpq__mod__cdata)
	{
		shear_rhs[j] -= AC_sshear__mod__cdata*Field(j+1)
		j += 3
	}
	if(AC_iuutest__mod__cdata != 0)
	{
		j = AC_iuutest__mod__cdata-1
		while(j < AC_iuztestpq__mod__cdata)
		{
			shear_rhs[j+1] -= AC_sshear1__mod__shear*Field(j)
			j += 3
		}
	}
	
  }
  if (AC_iam__mod__cdata!=0) shear_rhs[AC_iamx__mod__cdata-1] -=AC_sshear__mod__cdata*Field(AC_iamy__mod__cdata-1)
  **/
}
#endif
