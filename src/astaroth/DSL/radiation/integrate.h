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

//Vertical rays: the boundary intensity is set in the upstream ghost zone by rad_intensity_bc_z
Kernel revise_and_accumulate_Q_up(int idir)
{
	const int3 direction = {0,0,1}
	const real Q0 = extrinsic_ray_update(direction,Q_UP,TAU_UP)
	write(F_QRAD,F_QRAD+AC_weight__mod__radiation[idir-1]*Q0*F_KAPPARHO)
}

Kernel revise_and_accumulate_Q_down(int idir)
{
	const int3 direction = {0,0,-1}
	const real Q0 = extrinsic_ray_update(direction,Q_DOWN,TAU_DOWN)
	write(F_QRAD,F_QRAD+AC_weight__mod__radiation[idir-1]*Q0*F_KAPPARHO)
}

//Periodic horizontal rays on a single process in the ray direction (Qperiodic):
//the boundary heating rate is Q0 = Q(end)/(1-exp(-tau(end))) with the intrinsic values
//at the downstream end of the line, and then Q = Q + Q0*exp(-tau) (Qrevision), which is
//accumulated into Qrad directly. Q_HOR is not overwritten with the revised value since all
//points of the line read the intrinsic Q_HOR at ray_end: an in-place update would race with that.
revise_and_accumulate_periodic_ray(int3 direction, int idir)
{
	int3 ray_end = vertexIdx
	if (direction.x > 0) ray_end.x = AC_nlocal_max.x-1
	if (direction.x < 0) ray_end.x = NGHOST
	if (direction.y > 0) ray_end.y = AC_nlocal_max.y-1
	if (direction.y < 0) ray_end.y = NGHOST
	const real Q0 = Q_HOR[ray_end.x][ray_end.y][ray_end.z]/get_one_m_negexp(TAU_HOR[ray_end.x][ray_end.y][ray_end.z])
	write(F_QRAD,F_QRAD+AC_weight__mod__radiation[idir-1]*(Q_HOR+Q0*exp(-TAU_HOR))*F_KAPPARHO)
}

Kernel revise_and_accumulate_Q_xp(int idir)
{
	revise_and_accumulate_periodic_ray((int3){1,0,0},idir)
}

Kernel revise_and_accumulate_Q_xm(int idir)
{
	revise_and_accumulate_periodic_ray((int3){-1,0,0},idir)
}

Kernel revise_and_accumulate_Q_yp(int idir)
{
	revise_and_accumulate_periodic_ray((int3){0,1,0},idir)
}

Kernel revise_and_accumulate_Q_ym(int idir)
{
	revise_and_accumulate_periodic_ray((int3){0,-1,0},idir)
}

//f(:,:,:,iQrad) = f(:,:,:,iQrad) + weight(idir)*Qrad*f(:,:,:,ikapparho) of calcQ
Kernel reset_Qrad()
{
	write(F_QRAD,0.0)
}

