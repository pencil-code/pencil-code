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

accumulate(int idir, int3 direction, Field Q, Field TAU)
{
	Q_new = extrinsic_ray_update(direction,Q,TAU)
	write(F_QRAD,F_QRAD+AC_weight__mod__radiation[idir-1]*Q_new*F_KAPPARHO)
}

accumulate_vertical(int idir, int3 direction, Field Q, Field TAU)
{
	const real Q_intrinsic = value(Q)
	const real Q_extrinsic = get_incoming_boundary_ray(Q,direction)*exp(-value(TAU))
	//On the downstream boundary we have already revised so don't apply the extrinsic value.
	const bool downstream_boundary = (direction.z > 0) ? vertexIdx.z == AC_nlocal_max.z-1 : vertexIdx.z == NGHOST
	const real Q_revised = downstream_boundary ? Q_intrinsic : Q_intrinsic + Q_extrinsic
	write(F_QRAD,F_QRAD+AC_weight__mod__radiation[idir-1]*Q_revised*F_KAPPARHO)
}

Kernel accumulate_Q_up(int idir)
{
	accumulate_vertical(idir,(int3){0,0,1},Q_UP,TAU_UP)
}

Kernel accumulate_Q_down(int idir)
{
	accumulate_vertical(idir,(int3){0,0,-1},Q_DOWN,TAU_DOWN)
}

//Revision of only the downstream boundary layer (launched there by acRayUpdate), in place so that the received
//upstream ghost zone stays valid: this is the boundary heating rate sent to the next process along the ray
Kernel propagate_Q_up_boundary()
{
	Q_UP[vertexIdx.x][vertexIdx.y][vertexIdx.z] = Q_UP + get_incoming_boundary_ray(Q_UP,(int3){0,0,1})*exp(-TAU_UP)
}

Kernel propagate_Q_down_boundary()
{
	Q_DOWN[vertexIdx.x][vertexIdx.y][vertexIdx.z] = Q_DOWN + get_incoming_boundary_ray(Q_DOWN,(int3){0,0,-1})*exp(-TAU_DOWN)
}

Kernel accumulate_Q_xp(int idir)
{
	accumulate(idir,(int3){+1,0,0},Q_HOR,TAU_HOR)
}

Kernel accumulate_Q_xm(int idir)
{
	accumulate(idir,(int3){-1,0,0},Q_HOR,TAU_HOR)
}

Kernel accumulate_Q_yp(int idir)
{
	accumulate(idir,(int3){0,+1,0},Q_HOR,TAU_HOR)
}

Kernel accumulate_Q_ym(int idir)
{
	accumulate(idir,(int3){0,-1,0},Q_HOR,TAU_HOR)
}

Kernel reset_Qrad()
{
	write(F_QRAD,0.0)
}

