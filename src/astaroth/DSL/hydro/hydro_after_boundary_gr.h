//TP: hydro: hydro_after_boundary_gr
//Conserved-to-primitive inversion in curved spacetime (lgr_metric=T), see hydro.f90 for the equations.
//The metric is taken from special/z4c, bb=curl(aa)+B_ext=sqrt(gamma)*B^i is computed here since it is
//not stored on the GPU. Only the interior is computed, the halos are communicated when needed.
#if Lhydro_MODULE && Ldensity_MODULE && LZ4C
gr_matvec(Matrix m, real3 v)
{
	return real3(m[0][0]*v.x+m[0][1]*v.y+m[0][2]*v.z,
	             m[1][0]*v.x+m[1][1]*v.y+m[1][2]*v.z,
	             m[2][0]*v.x+m[2][1]*v.y+m[2][2]*v.z)
}
gr_zfunc(real zz, real EE, real S2, real B2, real SB, real w)
{
	v2 = min((S2+SB*SB*(2.*zz+B2)/(zz*zz))/((zz+B2)*(zz+B2)),1.-1e-12)
	return zz-w*zz*(1.-v2)/(1.+w)+0.5*B2*(1.+v2)-0.5*SB*SB/(zz*zz)-EE
}
//S_ij = S^k_i gamma_kj
gr_lower_smix(Matrix smix, Matrix gam, int i, int j)
{
	return smix[0][i]*gam[0][j]+smix[1][i]*gam[1][j]+smix[2][i]*gam[2][j]
}
Kernel hydro_after_boundary_gr()
{
	if(AC_lconservative__mod__hydro && AC_lgr_metric__mod__z4c)
	{
		cs20 = AC_cs20__mod__equationofstate
//
//  Metric: gamma_ij = gt_ij/chi, gamma^ij = chi gt^ij, sqrt(gamma) = sqrt(det(gt)/chi^3).
//
		chi = value(Field(AC_ichi__mod__z4c-1))
		alp = value(Field(AC_ialp__mod__z4c-1))
		g11 = value(Field(AC_igt__mod__z4c[0]-1))
		g22 = value(Field(AC_igt__mod__z4c[1]-1))
		g33 = value(Field(AC_igt__mod__z4c[2]-1))
		g12 = value(Field(AC_igt__mod__z4c[3]-1))
		g23 = value(Field(AC_igt__mod__z4c[4]-1))
		g13 = value(Field(AC_igt__mod__z4c[5]-1))
		beta = real3(value(Field(AC_ibeta__mod__z4c.x-1)),value(Field(AC_ibeta__mod__z4c.y-1)),value(Field(AC_ibeta__mod__z4c.z-1)))
		u11 = g22*g33-g23*g23
		u22 = g11*g33-g13*g13
		u33 = g11*g22-g12*g12
		u12 = g13*g23-g12*g33
		u23 = g12*g13-g11*g23
		u13 = g12*g23-g22*g13
		detgt = g11*u11+g12*u12+g13*u13
		sg = sqrt(detgt/(chi*chi*chi))
		gam  = (1./chi)*Matrix(real3(g11,g12,g13),real3(g12,g22,g23),real3(g13,g23,g33))
		gamu = (chi/detgt)*Matrix(real3(u11,u12,u13),real3(u12,u22,u23),real3(u13,u23,u33))

		EE = value(Field(AC_irho__mod__cdata-1))/sg
		Sd = value(F_UVEC)/sg
		Bu = real3(0.,0.,0.)
#if LMAGNETIC
		Bu = (curl(F_AVEC)+AC_b_ext__mod__magnetic)/sg
#endif
		Su = gr_matvec(gamu,Sd)
		Bd = gr_matvec(gam,Bu)
		S2 = dot(Su,Sd)
		B2 = dot(Bu,Bd)
		SB = dot(Sd,Bu)
//
//  Start from the exact solution without magnetic field for the fluid
//  energy E-B^2/2, then do Newton iterations.
//
		e = max(EE-0.5*B2,AC_tini__mod__cparam)
		zz = 0.5*((1.+cs20)*e+sqrt(max(((1.+cs20)*e)*((1.+cs20)*e)-4.*cs20*S2,0.)))
		zz = max(zz,AC_tini__mod__cparam)
		converged = false
		for iter in 1:AC_niter_gr__mod__hydro+1 {
			if(!converged)
			{
				ff = gr_zfunc(zz,EE,S2,B2,SB,cs20)
				if(abs(ff)/EE < 1e-14)
				{
					converged = true
				}
				else
				{
					dz = 1e-7*zz
					ffp = (gr_zfunc(zz+dz,EE,S2,B2,SB,cs20)-ff)/dz
					zz = max(zz-ff/ffp,0.5*zz)
				}
			}
		}
//
//  Primitives: v_i, v^i, W^2, e and the ideal electric field E_i = -sqrt(gamma) [ijk] v^j B^k.
//
		vd = (Sd+(SB/zz)*Bd)/(zz+B2)
		vu = gr_matvec(gamu,vd)
		v2 = min(dot(vu,vd),1.-1e-12)
		e = zz*(1.-v2)/(1.+cs20)
		Ed = -sg*cross(vu,Bu)
		Eu = gr_matvec(gamu,Ed)
		E2 = dot(Eu,Ed)
		ptot = cs20*e+0.5*(E2+B2)
//
//  Stress tensor S^j_i (smix[j][i]) and fluxes.
//
		smix = Matrix(zz*vu.x*vd-Eu.x*Ed-Bu.x*Bd+real3(ptot,0.,0.),
		              zz*vu.y*vd-Eu.y*Ed-Bu.y*Bd+real3(0.,ptot,0.),
		              zz*vu.z*vd-Eu.z*Ed-Bu.z*Bd+real3(0.,0.,ptot))
		fe = sg*(alp*Su-EE*beta)
		F_FGREX[vertexIdx.x][vertexIdx.y][vertexIdx.z] = fe.x
		F_FGREY[vertexIdx.x][vertexIdx.y][vertexIdx.z] = fe.y
		F_FGREZ[vertexIdx.x][vertexIdx.y][vertexIdx.z] = fe.z
		//F^j_i at iFgrS+3*(i-1)+j-1
                F_FGRS_1[vertexIdx.x][vertexIdx.y][vertexIdx.z] = sg*(alp*smix[0][0]-beta.x*Sd.x)
                F_FGRS_2[vertexIdx.x][vertexIdx.y][vertexIdx.z] = sg*(alp*smix[1][0]-beta.y*Sd.x)
                F_FGRS_3[vertexIdx.x][vertexIdx.y][vertexIdx.z] = sg*(alp*smix[2][0]-beta.z*Sd.x)
                F_FGRS_4[vertexIdx.x][vertexIdx.y][vertexIdx.z] = sg*(alp*smix[0][1]-beta.x*Sd.y)
                F_FGRS_5[vertexIdx.x][vertexIdx.y][vertexIdx.z] = sg*(alp*smix[1][1]-beta.y*Sd.y)
                F_FGRS_6[vertexIdx.x][vertexIdx.y][vertexIdx.z] = sg*(alp*smix[2][1]-beta.z*Sd.y)
                F_FGRS_7[vertexIdx.x][vertexIdx.y][vertexIdx.z] = sg*(alp*smix[0][2]-beta.x*Sd.z)
                F_FGRS_8[vertexIdx.x][vertexIdx.y][vertexIdx.z] = sg*(alp*smix[1][2]-beta.y*Sd.z)
                F_FGRS_9[vertexIdx.x][vertexIdx.y][vertexIdx.z] = sg*(alp*smix[2][2]-beta.z*Sd.z)
		VX[vertexIdx.x][vertexIdx.y][vertexIdx.z] = alp*vu.x-beta.x
		VY[vertexIdx.x][vertexIdx.y][vertexIdx.z] = alp*vu.y-beta.y
		VZ[vertexIdx.x][vertexIdx.y][vertexIdx.z] = alp*vu.z-beta.z
		F_LORENTZ[vertexIdx.x][vertexIdx.y][vertexIdx.z] = 1./(1.-v2)
		F_EGR[vertexIdx.x][vertexIdx.y][vertexIdx.z] = e
//
//  ADM sources (if the Eadm etc. slots exist), stored in the order xx,yy,zz,xy,yz,xz.
//
		if(AC_ieadm__mod__z4c != 0)
		{
			F_EADM[vertexIdx.x][vertexIdx.y][vertexIdx.z] = EE
			F_SADMX[vertexIdx.x][vertexIdx.y][vertexIdx.z] = Sd.x
			F_SADMY[vertexIdx.x][vertexIdx.y][vertexIdx.z] = Sd.y
			F_SADMZ[vertexIdx.x][vertexIdx.y][vertexIdx.z] = Sd.z
			F_SADMIJ_1[vertexIdx.x][vertexIdx.y][vertexIdx.z] = gr_lower_smix(smix,gam,0,0)
			F_SADMIJ_2[vertexIdx.x][vertexIdx.y][vertexIdx.z] = gr_lower_smix(smix,gam,1,1)
			F_SADMIJ_3[vertexIdx.x][vertexIdx.y][vertexIdx.z] = gr_lower_smix(smix,gam,2,2)
			F_SADMIJ_4[vertexIdx.x][vertexIdx.y][vertexIdx.z] = gr_lower_smix(smix,gam,0,1)
			F_SADMIJ_5[vertexIdx.x][vertexIdx.y][vertexIdx.z] = gr_lower_smix(smix,gam,1,2)
			F_SADMIJ_6[vertexIdx.x][vertexIdx.y][vertexIdx.z] = gr_lower_smix(smix,gam,0,2)
		}
	}
}
#else
Kernel hydro_after_boundary_gr()
{
}
#endif
