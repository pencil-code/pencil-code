#if LVISCOSITY && LHYDRO
Kernel viscosity_after_boundary()
{
	if(AC_lrate_of_strain_as_aux__mod__viscosity)
	{
		if(AC_lconservative__mod__hydro)
		{
			uij = gradient_tensor(MOM,RHO)
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
