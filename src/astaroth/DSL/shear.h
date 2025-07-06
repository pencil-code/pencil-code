#if LSHEAR

if(lshear)
{
  if(!AC_lshearadvection_as_shift__mod__shear)
  {
	DF_UU  -= AC_uy0__mod__shear[vertexIdx.x-NGHOST]*dery(UU)
	DF_AA  -= AC_uy0__mod__shear[vertexIdx.x-NGHOST]*dery(AA)
	DF_SS  -= AC_uy0__mod__shear[vertexIdx.x-NGHOST]*dery(SS)
	DF_RHO -= AC_uy0__mod__shear[vertexIdx.x-NGHOST]*dery(RHO)
  }
  if (lhydro && AC_lshear_acceleration__mod__shear) DF_UU.y  -= AC_sshear1__mod__shear * UU.x
  if(AC_lhyper3x_mesh__mod__shear)
  {
  	d = AC_diff_hyper3x_mesh__mod__shear*abs(AC_sshear__mod__cdata)
	DF_UU  += d*der6x_ignore_spacing(UU)
	DF_AA  += d*der6x_ignore_spacing(AA)
	DF_RHO += d*der6x_ignore_spacing(RHO)
	DF_SS  += d*der6x_ignore_spacing(SS)
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
  /**
  if(ldustvelocity)
  {
	  for k in 0:ndustspec
	  {
		int dust_spec_y = AC_iudy__mod__cdata[k]
		int dust_spec_x = AC_iudx__mod__cdata[k]
		shear_rhs[dust_spec_y] -= AC_sshear1__mod__shear*Field(dust_spec_x)
	  }
  }
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
