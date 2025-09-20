Field Q_UP
Field Q_DOWN

Field TAU_UP
Field TAU_DOWN
Kernel integrate_Q_up()
{
	const int nrad=1;
	const int3 direction = (int3){0,0,nrad};
	ac_ray_func(direction,F_KAPPARHO,SRAD,Q_UP,TAU_UP,
			0.5*(      AC_dlength__mod__radiation[vertexIdx.z-nrad][1]
				+  AC_dlength__mod__radiation[vertexIdx.z][1]
			    ),
			0.5*(      AC_dlength__mod__radiation[vertexIdx.z+nrad][1]
				+  AC_dlength__mod__radiation[vertexIdx.z][1]
			    )
		)
}

Kernel integrate_Q_down()
{
	const int nrad=-1;
	const int3 direction = (int3){0,0,nrad};
	ac_ray_func(direction,F_KAPPARHO,SRAD,Q_DOWN,TAU_DOWN,
			0.5*(      AC_dlength__mod__radiation[vertexIdx.z-nrad][0]
				+  AC_dlength__mod__radiation[vertexIdx.z][0]
			    ),
			0.5*(      AC_dlength__mod__radiation[vertexIdx.z+nrad][0]
				+  AC_dlength__mod__radiation[vertexIdx.z][0]
			    )
		)
}

Kernel revise_Q_up()
{
	const int3 direction = {0,0,1}
	write(Q_UP,extrinsic_ray_update(direction,Q_UP,TAU_UP))
}

Kernel revise_Q_down()
{
	const int3 direction = {0,0,-1}
	write(Q_DOWN,extrinsic_ray_update(direction,Q_DOWN,TAU_DOWN))
}

Kernel sum_up_rays()
{
	down_idir = 1
	up_idir = 2
	res  = AC_weight__mod__radiation[up_idir-1]*Q_UP*F_KAPPARHO
	res += AC_weight__mod__radiation[down_idir-1]*Q_DOWN*F_KAPPARHO
	write(F_QRAD,res)
}
