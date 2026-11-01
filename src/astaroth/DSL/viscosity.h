#if LVISCOSITY && LHYDRO
Kernel viscosity_after_boundary()
{
	if(AC_lnusmag_as_aux__mod__viscosity)
	{
		fatal_error_message(AC_lshear_rateofstrain__mod__hydro,"viscosity_after_boundary: lnusmag_as_aux with lshear_rateofstrain not implemented on the GPU")
		Matrix uij_smag
		if(AC_ivv__mod__cdata != 0)
		{
			uij_smag = gradient_tensor(F_VVEC)
		}
		else
		{
			uij_smag = gradient_tensor(F_UVEC)
		}
		sij_smag = traceless_strain(uij_smag)
		sij2_smag = 0.0
		for i in 0:3 {
			for j in 0:3 {
				sij2_smag += sij_smag[i][j]*sij_smag[i][j]
			}
		}
		nu_smag = (AC_c_smag__mod__viscosity*AC_dxmax__mod__cdata)*(AC_c_smag__mod__viscosity*AC_dxmax__mod__cdata)*sqrt(2.*sij2_smag)
		if(AC_gamma_smag__mod__viscosity != 0.)
		{
			nu_smag = nu_smag/sqrt(1.+AC_gamma_smag__mod__viscosity*sij2_smag)
		}
		F_NUSMAG[vertexIdx.x][vertexIdx.y][vertexIdx.z] = nu_smag
	}
	if(AC_lrate_of_strain_as_aux__mod__viscosity)
	{
		if(AC_lconservative__mod__hydro)
		{
			uij = gradient_tensor(MOM,F_RHO)
			Sij = traceless_strain(uij)
			write(F_SIJ_0,Sij[0][0])
			write(F_SIJ_1,Sij[1][1])
			write(F_SIJ_2,Sij[2][2])
			write(F_SIJ_3,Sij[0][1])
			write(F_SIJ_4,Sij[0][2])
			write(F_SIJ_5,Sij[1][2])
		}
	}
}
#else
Kernel viscosity_after_boundary(){}
#endif
