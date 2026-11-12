//TP: special/z4c: special_before_boundary
//Enforce the algebraic constraints det(gt)=1 and gt^ij At_ij=0 and keep chi above chi_floor.
//The symmetric tensors are stored in the order xx,yy,zz,xy,yz,xz (csym).
#if LZ4C
Kernel z4c_before_boundary()
{
	if(AC_lenforce_algebraic__mod__z4c)
	{
		g11 = value(F_GT_1)
		g22 = value(F_GT_2)
		g33 = value(F_GT_3)
		g12 = value(F_GT_4)
		g23 = value(F_GT_5)
		g13 = value(F_GT_6)
		u11 = g22*g33-g23*g23
		u22 = g11*g33-g13*g13
		u33 = g11*g22-g12*g12
		u12 = g13*g23-g12*g33
		u23 = g12*g13-g11*g23
		u13 = g12*g23-g22*g13
		detg = g11*u11+g12*u12+g13*u13
		fac = pow(detg,-1./3.)
		//gt^ij of the rescaled gt_ij
		ufac = 1./(detg*fac)
		trA = ufac*(u11*value(F_AT_1)+u22*value(F_AT_2)+u33*value(F_AT_3)
		      +2.*(u12*value(F_AT_4)+u23*value(F_AT_5)+u13*value(F_AT_6)))
		F_GT_1[vertexIdx.x][vertexIdx.y][vertexIdx.z] = fac*g11
		F_GT_2[vertexIdx.x][vertexIdx.y][vertexIdx.z] = fac*g22
		F_GT_3[vertexIdx.x][vertexIdx.y][vertexIdx.z] = fac*g33
		F_GT_4[vertexIdx.x][vertexIdx.y][vertexIdx.z] = fac*g12
		F_GT_5[vertexIdx.x][vertexIdx.y][vertexIdx.z] = fac*g23
		F_GT_6[vertexIdx.x][vertexIdx.y][vertexIdx.z] = fac*g13

		F_AT_1[vertexIdx.x][vertexIdx.y][vertexIdx.z] = F_AT_1-fac*g11*trA/3.
		F_AT_2[vertexIdx.x][vertexIdx.y][vertexIdx.z] = F_AT_2-fac*g22*trA/3.
		F_AT_3[vertexIdx.x][vertexIdx.y][vertexIdx.z] = F_AT_3-fac*g33*trA/3.
		F_AT_4[vertexIdx.x][vertexIdx.y][vertexIdx.z] = F_AT_4-fac*g12*trA/3.
		F_AT_5[vertexIdx.x][vertexIdx.y][vertexIdx.z] = F_AT_5-fac*g23*trA/3.
		F_AT_6[vertexIdx.x][vertexIdx.y][vertexIdx.z] = F_AT_6-fac*g13*trA/3.

		F_CHI[vertexIdx.x][vertexIdx.y][vertexIdx.z] = max(value(F_CHI),AC_chi_floor__mod__z4c)
	}
}
#else
Kernel z4c_before_boundary()
{
}
#endif
