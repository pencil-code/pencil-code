real shear_rhs[AC_nvar__mod__cdata]
for i in 0:AC_nvar__mod__cdata
{
  shear_rhs[i] = 0.0
}

if(lshear)
{
  if(!AC_lshearadvection_as_shift__mod__shear)
  {
  	for j in 0:AC_nvar__mod__cdata
  	{
  		if(lbfield && (j >= AC_ibx__mod__cdata) && (j <= AC_ibz__mod__cdata)) continue
            	if (ltestflow && (j >= AC_iuutest__mod__cdata) && (j <= AC_iuutest__mod__cdata+AC_ntestflow__mod__cdata-1)) continue
		shear_rhs[j] -= AC_uy0__mod__shear[vertexIdx.x-NGHOST]*deryy(Field(j))
  	}
  }
  if (lhydro && AC_lshear_acceleration__mod__shear) shear_rhs[UU.y]  -= AC_sshear1__mod__shear * UU.x
  if(AC_lhyper3x_mesh__mod__shear)
  {
  	d = AC_diff_hyper3x_mesh__mod__shear*abs(AC_sshear__mod__cdata)
  	for j in  0:AC_nvar__mod__cdata
  	{
            if ((lbfield && AC_ibx__mod__cdata <= j && j <= AC_ibz__mod__cdata) ||
                (lpscalar && AC_icc__mod__cdata <= j && j <= AC_icc__mod__cdata+npscalar-1)) continue
           shear_rhs[j] += d*der6x_ignore_spacing(Field(j))
  	}
  	//TP: not workable in handwritten code since no global maxdiffus3
          //if (lupdate_courant_dt) then
          //  diffus_shear3 = d
          //  maxdiffus3=max(maxdiffus3,diffus_shear3)
          //endif
  }
  if(ldustvelocity)
  {
	  for k in 0:ndustspec
	  {
		int dust_spec_y = AC_iudy__mod__cdata[k]
		int dust_spec_x = AC_iudx__mod__cdata[k]
		shear_rhs[dust_spec_y] -= AC_sshear1__mod__shear*Field(dust_spec_x)
	  }
  }
  if (lmagnetic && !lbfield && AC_lmagnetic_stretching__mod__shear)
  {
        shear_rhs[AAX] -= AC_sshear__mod__cdata*AAY
  	if(AC_lmagnetic_tilt__mod__shear)
  	{
	  shear_rhs[AAX] -= AC_sshear_sini__mod__shear*AAX
	  shear_rhs[AAY] += AC_sshear_sini__mod__shear*AAY
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
}
