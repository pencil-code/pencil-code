//Radiative transfer along the rays of radiation_ray.f90 (rays along the coordinate axes only).
//Each ray direction is integrated separately (radTransfer in gpu_astaroth.cc loops over the rays in the
//order of initialize_radiation); idir (AC_ray_idir in the ComputeSteps) is the index of the current ray for dlength and weight.
//Vertical rays use Q_UP/TAU_UP and Q_DOWN/TAU_DOWN, horizontal ones share Q_HOR/TAU_HOR.
Field Q_UP
Field Q_DOWN

Field TAU_UP
Field TAU_DOWN

Field Q_HOR
Field TAU_HOR

//Line element 0.5*(dlength(n-nrad)+dlength(n)) and 0.5*(dlength(n)+dlength(n+nrad)) of Qintrinsic
rad_ray_func(int3 direction, Field Q, Field TAU, int idir)
{
	const int nrad = direction.z
	ac_ray_func(direction,F_KAPPARHO,SRAD,Q,TAU,
			0.5*(      AC_dlength__mod__radiation[vertexIdx.z-nrad][idir-1]
				+  AC_dlength__mod__radiation[vertexIdx.z][idir-1]
			    ),
			0.5*(      AC_dlength__mod__radiation[vertexIdx.z+nrad][idir-1]
				+  AC_dlength__mod__radiation[vertexIdx.z][idir-1]
			    )
		)
}

Kernel integrate_Q_up(int idir)
{
	const int3 direction = (int3){0,0,1};
	rad_ray_func(direction,Q_UP,TAU_UP,idir)
}

Kernel integrate_Q_down(int idir)
{
	const int3 direction = (int3){0,0,-1};
	rad_ray_func(direction,Q_DOWN,TAU_DOWN,idir)
}

Kernel integrate_Q_xp(int idir)
{
	const int3 direction = (int3){1,0,0};
	rad_ray_func(direction,Q_HOR,TAU_HOR,idir)
}

Kernel integrate_Q_xm(int idir)
{
	const int3 direction = (int3){-1,0,0};
	rad_ray_func(direction,Q_HOR,TAU_HOR,idir)
}

Kernel integrate_Q_yp(int idir)
{
	const int3 direction = (int3){0,1,0};
	rad_ray_func(direction,Q_HOR,TAU_HOR,idir)
}

Kernel integrate_Q_ym(int idir)
{
	const int3 direction = (int3){0,-1,0};
	rad_ray_func(direction,Q_HOR,TAU_HOR,idir)
}

revise_and_accumulate(int idir, int3 direction, Field Q, Field TAU)
{
	Q_new = extrinsic_ray_update(direction,Q,TAU)
	write(Q,Q_new)
	write(F_QRAD,F_QRAD+AC_weight__mod__radiation[idir-1]*Q_new*F_KAPPARHO)
}

//Vertical rays: the boundary intensity is set in the upstream ghost zone by rad_intensity_bc_z
Kernel revise_and_accumulate_Q_up(int idir)
{
	revise_and_accumulate(idir,(int3){0,0,1},Q_UP,TAU_UP)
}

Kernel revise_and_accumulate_Q_down(int idir)
{
	revise_and_accumulate(idir,(int3){0,0,-1},Q_DOWN,TAU_DOWN)
}

Kernel revise_and_accumulate_Q_xp(int idir)
{
	revise_and_accumulate(idir,(int3){+1,0,0},Q_HOR,TAU_HOR)
}

Kernel revise_and_accumulate_Q_xm(int idir)
{
	revise_and_accumulate(idir,(int3){-1,0,0},Q_HOR,TAU_HOR)
}

Kernel revise_and_accumulate_Q_yp(int idir)
{
	revise_and_accumulate(idir,(int3){0,+1,0},Q_HOR,TAU_HOR)
}

Kernel revise_and_accumulate_Q_ym(int idir)
{
	revise_and_accumulate(idir,(int3){0,-1,0},Q_HOR,TAU_HOR)
}

Kernel reset_Qrad()
{
	write(F_QRAD,0.0)
}

