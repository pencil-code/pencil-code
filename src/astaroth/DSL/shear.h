real shear_rhs[nvar]
for i in 0:nvar
{
  shear_rhs[i] = 0.0
}

if(lshear)
{
  if(!lshearadvection_as_shift)
  {
  	for j in 0:nvar
  	{
  		if(lbfield && (j >= ibx) && (j <= ibz)) continue
            	if (ltestflow && (j >= iuutest) && (j <= iuutest+ntestflow-1)) continue
		shear_rhs[j] -= uy0[vertexIdx.x-NGHOST]*deryy(Field(j))
  	}
  }
  if (lhydro && lshear_acceleration) shear_rhs[UU.y]  -= sshear1 * UU.x
  if(lhyper3x_mesh)
  {
  	d = diff_hyper3x_mesh*abs(sshear)
  	for j in  0:nvar
  	{
            if ((lbfield && ibx <= j && j <= ibz) ||
                (lpscalar && icc <= j && j <= icc+npscalar-1)) continue
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
		int dust_spec_y = iudy[k]
		int dust_spec_x = iudx[k]
		shear_rhs[dust_spec_y] -= sshear1*Field(dust_spec_x)
	  }
  }
  if (lmagnetic && !lbfield && lmagnetic_stretching)
  {
        shear_rhs[AAX] -= sshear*AAY
  	if(lmagnetic_tilt)
  	{
	  shear_rhs[AAX] -= sshear_sini*AAX
	  shear_rhs[AAY] += sshear_sini*AAY
  	}
  }
  if(ltestfield)
  {
	int j = iaatest-1
	while(j < iaztestpq)
	{
		shear_rhs[j] -= sshear*Field(j+1)
		j += 3
	}
	if(iuutest != 0)
	{
		j = iuutest-1
		while(j < iuztestpq)
		{
			shear_rhs[j+1] -= sshear1*Field(j)
			j += 3
		}
	}
	
  }
  if (iam!=0) shear_rhs[iamx-1] -=sshear*Field(iamy-1)
}
