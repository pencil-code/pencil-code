#if LIMPLICIT_DIFFUSION
//TP: Hacky way to make sure the buffers to hold the Fourier space values are allocated
Kernel initialize_implicit_diffusion_buffers()
{
	if(AC_limplicit_diffusion_with_fft__mod__implicitdiffusion)
	{ write(SPLIT_DIFFUSION_UPDATE_BUFFER_REAL,0.0)
		write(SPLIT_DIFFUSION_UPDATE_BUFFER_IMAG,0.0)
	}
}

//Hard-coded for nwo
Field CG_R[6]
Field CG_P[6]
Field CG_RHS[6]
const Field CG_LHS = 
[
	F_BX,
	F_BY,
	F_BZ,
	UUX,
	UUY,
	UUZ
]



enum CG_FIELD
{
	CG_FIELD_0,
	CG_FIELD_1,
	CG_FIELD_2,
	CG_FIELD_3,
	CG_FIELD_4,
	CG_FIELD_5
}

//TP: using implicit euler the update becomes (I-dt*nu*nabla^2)u^{n+1} = u^n
implicit_euler_lhs(Field f, real dt, real diffusion_coeff)
{
	return f - dt*diffusion_coeff*laplace(f)
}

Kernel ID_init_cg_residual_kernel(real dt)
{
	if(AC_limplicit_diffusion_with_cg__mod__implicitdiffusion)
	{
		if(AC_limplicit_resistivity__mod__magnetic)
		{
			write(CG_R[0],F_BX-implicit_euler_lhs(F_BX,dt,AC_eta__mod__magnetic))
			write(CG_R[1],F_BY-implicit_euler_lhs(F_BY,dt,AC_eta__mod__magnetic))
			write(CG_R[2],F_BZ-implicit_euler_lhs(F_BZ,dt,AC_eta__mod__magnetic))

			write(CG_P[0],F_BX-implicit_euler_lhs(F_BX,dt,AC_eta__mod__magnetic))
			write(CG_P[1],F_BY-implicit_euler_lhs(F_BY,dt,AC_eta__mod__magnetic))
			write(CG_P[2],F_BZ-implicit_euler_lhs(F_BZ,dt,AC_eta__mod__magnetic))

			write(CG_RHS[0],F_BX)
			write(CG_RHS[1],F_BY)
			write(CG_RHS[2],F_BZ)
		}
		if(AC_limplicit_viscosity__mod__viscosity)
		{
			write(CG_R[3],UUX-implicit_euler_lhs(UUX,dt,AC_nu__mod__viscosity))
			write(CG_R[4],UUY-implicit_euler_lhs(UUY,dt,AC_nu__mod__viscosity))
			write(CG_R[5],UUZ-implicit_euler_lhs(UUZ,dt,AC_nu__mod__viscosity))

			write(CG_P[3],UUX-implicit_euler_lhs(UUX,dt,AC_nu__mod__viscosity))
			write(CG_P[4],UUY-implicit_euler_lhs(UUY,dt,AC_nu__mod__viscosity))
			write(CG_P[5],UUZ-implicit_euler_lhs(UUZ,dt,AC_nu__mod__viscosity))

			write(CG_RHS[3],UUX)
			write(CG_RHS[4],UUY)
			write(CG_RHS[5],UUZ)
		}
	}
}

Kernel ID_cg_compute_alpha(real dt, CG_FIELD f, real diffusion_coeff)
{

	if(AC_limplicit_diffusion_with_cg__mod__implicitdiffusion)
	{
		Ap = implicit_euler_lhs(CG_P[f],dt,diffusion_coeff)
		r = value(CG_R[f])
		p = value(CG_P[f])
		cg_compute_inner_products(Ap,p,r)
	}
}

Kernel ID_cg_advance(real dt, CG_FIELD f, real diffusion_coeff)
{
	if(AC_limplicit_diffusion_with_cg__mod__implicitdiffusion)
	{
		Ap = implicit_euler_lhs(CG_P[f],dt,diffusion_coeff)
		p = value(CG_P[f])
		cg_advance_solution(CG_LHS[f],Ap,p,CG_R[f])
	}
}

Kernel ID_cg_advance_p(CG_FIELD f)
{
	if(AC_limplicit_diffusion_with_cg__mod__implicitdiffusion)
	{
		cg_next_direction(CG_P[f],CG_R[f])
	}
}

global output real AC_implicit_diffusion_residual

Kernel ID_get_residual(CG_FIELD f, real dt, real diffusion_coeff)
{
	//It is well known that one should not use the implicit residual inside CG to determine determination since it is not as accurate
	//as computing the residual outside the solver
	if(AC_limplicit_diffusion_with_cg__mod__implicitdiffusion)
	{
		r = CG_RHS[f] - implicit_euler_lhs(CG_LHS[f],dt,diffusion_coeff)
		reduce_sum(r*r,AC_implicit_diffusion_residual)
	}
}

input CG_FIELD AC_CG_FIELD
input real AC_implicit_diffusion_coefficient 

ComputeSteps
implicit_diffusion_init_cg(boundconds)
{
	ID_init_cg_residual_kernel(AC_dt)
}
ComputeSteps
implicit_diffusion_cg_step(boundconds)
{
	ID_cg_compute_alpha(AC_dt,AC_CG_FIELD,AC_implicit_diffusion_coefficient)
	ID_cg_advance(AC_dt,AC_CG_FIELD,AC_implicit_diffusion_coefficient)
	ID_cg_advance_p(AC_CG_FIELD)
}

ComputeSteps
implicit_diffusion_get_residual(boundconds)
{
	ID_get_residual(AC_CG_FIELD,AC_dt,AC_implicit_diffusion_coefficient)
}



#endif
