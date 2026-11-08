//Tracer particles (particles_tracers.f90) on the GPU.
//
//The particle state lives in gmem arrays with capacity npar (>= mpar_loc), stored as
//separate components (structure of arrays). They are loaded from / stored to the
//host fp array by gpu_astaroth.cc. The kernels are launched directly over the volume
//(0,0,0)-(n,1,1), i.e. vertexIdx.x is the particle (or list) index. There are deliberately no
//"if (index < n)" guards: only threads inside the launch volume run, and with such guards the kernel
//analysis (run with zero kernel parameters) would skip the bodies, so acc would treat the particle
//arrays as unused and remove them.
//
//particles_tracer_step is launched after the halos of the velocity have been updated.
//Per RK substep it does what dxxp_dt_pencil + particles_timestep_second +
//boundconds_particles do on the CPU for Cartesian, periodic, non-shearing boxes:
//  dfp = alpha*dfp + u(xp)   (alpha, beta: low-storage RK coefficients of AC_rk_order,
//  xp  = xp + dt*beta*dfp     u interpolated linearly as in interpolate_linear)
//  periodic wrap of xp into the global box [xyz0,xyz1)
//  AC_par_leaving = 1 if xp is outside of the local subdomain (as in migrate_particles)
//
//Migration between processes (gpu_astaroth.cc, particles_migrate): the host reads
//AC_par_leaving, uploads index lists to AC_par_idx/AC_par_idx2 and uses
//  particles_pack:   records of the leaving particles   -> AC_par_buf (sent by MPI)
//  particles_move:   particles from the tail into the holes left by them
//  particles_unpack: received records in AC_par_buf      -> appended after npar_loc
//Index lists are stored as reals to need only the real-array transfer API.

#if LPARTICLES
gmem real AC_fp_xp[AC_npar]
gmem real AC_fp_yp[AC_npar]
gmem real AC_fp_zp[AC_npar]
gmem real AC_fp_vpx[AC_npar]
gmem real AC_fp_vpy[AC_npar]
gmem real AC_fp_vpz[AC_npar]
gmem real AC_dfp_xp[AC_npar]
gmem real AC_dfp_yp[AC_npar]
gmem real AC_dfp_zp[AC_npar]
gmem real AC_fp_ipar[AC_npar]

gmem real AC_par_leaving[AC_npar]
gmem real AC_par_idx[AC_npar]
gmem real AC_par_idx2[AC_npar]
//Particle records: column-major, so record i occupies AC_par_buf[0:PAR_NREC-1][i] contiguously
//Order: xp,yp,zp,vpx,vpy,vpz,dfp_xp,dfp_yp,dfp_zp,ipar (must agree with gpu_astaroth.cc)
const int PAR_NREC = 10
gmem real AC_par_buf[PAR_NREC][AC_npar]

//Index of the grid point at or below pos in one direction, as in map_nearest_grid
//(with lcheck_exact_frontier) followed by the correction in interpolate_linear.
//Clamped so that i0 and i0+1 are always valid indices, also for garbage positions.
particle_lower_index(real pos, real coord_first, real dcoord_1, int nloc)
{
	int i0 = int(round((pos - coord_first)*dcoord_1)) + NGHOST
	if (i0 < NGHOST) i0 = NGHOST
	if (i0 > NGHOST+nloc-1) i0 = NGHOST+nloc-1
	return i0
}

particle_interpolate_linear(Field fld, int ix0, int iy0, int iz0, real xp0, real yp0, real zp0, real dx1, real dy1, real dz1)
{
	g1 = fld[ix0  ][iy0  ][iz0  ]
	g2 = fld[ix0+1][iy0  ][iz0  ]
	g3 = fld[ix0  ][iy0+1][iz0  ]
	g4 = fld[ix0+1][iy0+1][iz0  ]
	g5 = fld[ix0  ][iy0  ][iz0+1]
	g6 = fld[ix0+1][iy0  ][iz0+1]
	g7 = fld[ix0  ][iy0+1][iz0+1]
	g8 = fld[ix0+1][iy0+1][iz0+1]
	return g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) +
	       xp0*yp0*dx1*dy1*(g1-g2-g3+g4) + xp0*zp0*dx1*dz1*(g1-g2-g5+g6) +
	       yp0*zp0*dy1*dz1*(g1-g3-g5+g7) +
	       xp0*yp0*zp0*dx1*dy1*dz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
}

particle_periodic_wrap(real pos, real pos0, real pos1, real period)
{
	real res = pos
	if (res <  pos0) res = res + period
	if (res >= pos1) res = res - period
	return res
}

//Low-storage RK coefficients of substep step_num (0-based), the same tables the fluid uses
particle_rk_alpha(int step_num)
{
	if (AC_rk_order == 1) return rk1_alpha[step_num]
	if (AC_rk_order == 2) return rk2_alpha[step_num]
	return rk3_alpha[step_num]
}
particle_rk_beta(int step_num)
{
	if (AC_rk_order == 1) return rk1_beta[step_num]
	if (AC_rk_order == 2) return rk2_beta[step_num]
	return rk3_beta[step_num]
}

//xx0_*, xx1_*: bounds of the local subdomain relative to xyz0 (xyz0_loc-xyz0, xyz1_loc-xyz0)
//as used in migrate_particles; the host passes -/+huge in directions without migration.
Kernel particles_tracer_step(int step_num, real dt_step,
                             real xx0_x, real xx0_y, real xx0_z, real xx1_x, real xx1_y, real xx1_z)
{
	const real alpha   = particle_rk_alpha(step_num)
	const real dt_beta = dt_step*particle_rk_beta(step_num)
	const int k = vertexIdx.x
	if (k >= 0)   //always true; bare blocks are not allowed in the DSL
	{
		const real xp = AC_fp_xp[k]
		const real yp = AC_fp_yp[k]
		const real zp = AC_fp_zp[k]

		int ix0 = NGHOST
		int iy0 = NGHOST
		int iz0 = NGHOST
		real xp0 = 0.0
		real yp0 = 0.0
		real zp0 = 0.0
		if (AC_nxgrid != 1)
		{
			ix0 = particle_lower_index(xp, AC_x__mod__cdata[NGHOST], AC_dx_1__mod__cdata[NGHOST], AC_nlocal.x)
			if (AC_x__mod__cdata[ix0] > xp) ix0 = ix0-1
			xp0 = xp - AC_x__mod__cdata[ix0]
		}
		if (AC_nygrid != 1)
		{
			iy0 = particle_lower_index(yp, AC_y__mod__cdata[NGHOST], AC_dy_1__mod__cdata[NGHOST], AC_nlocal.y)
			if (AC_y__mod__cdata[iy0] > yp) iy0 = iy0-1
			yp0 = yp - AC_y__mod__cdata[iy0]
		}
		if (AC_nzgrid != 1)
		{
			iz0 = particle_lower_index(zp, AC_z__mod__cdata[NGHOST], AC_dz_1__mod__cdata[NGHOST], AC_nlocal.z)
			if (AC_z__mod__cdata[iz0] > zp) iz0 = iz0-1
			zp0 = zp - AC_z__mod__cdata[iz0]
		}
		const real dx1 = AC_dx_1__mod__cdata[ix0]
		const real dy1 = AC_dy_1__mod__cdata[iy0]
		const real dz1 = AC_dz_1__mod__cdata[iz0]

		const real ux = particle_interpolate_linear(F_UX, ix0, iy0, iz0, xp0, yp0, zp0, dx1, dy1, dz1)
		const real uy = particle_interpolate_linear(F_UY, ix0, iy0, iz0, xp0, yp0, zp0, dx1, dy1, dz1)
		const real uz = particle_interpolate_linear(F_UZ, ix0, iy0, iz0, xp0, yp0, zp0, dx1, dy1, dz1)

		//vvp is an auxiliary equal to the gas velocity at the particle (only for output)
		AC_fp_vpx[k] = ux
		AC_fp_vpy[k] = uy
		AC_fp_vpz[k] = uz

		real dxp = alpha*AC_dfp_xp[k]
		real dyp = alpha*AC_dfp_yp[k]
		real dzp = alpha*AC_dfp_zp[k]
		if (AC_nxgrid != 1) dxp = dxp + ux
		if (AC_nygrid != 1) dyp = dyp + uy
		if (AC_nzgrid != 1) dzp = dzp + uz
		AC_dfp_xp[k] = dxp
		AC_dfp_yp[k] = dyp
		AC_dfp_zp[k] = dzp

		const real xn = particle_periodic_wrap(xp + dt_beta*dxp, AC_xyz0__mod__cdata.x, AC_xyz1__mod__cdata.x, AC_lxyz__mod__cdata.x)
		const real yn = particle_periodic_wrap(yp + dt_beta*dyp, AC_xyz0__mod__cdata.y, AC_xyz1__mod__cdata.y, AC_lxyz__mod__cdata.y)
		const real zn = particle_periodic_wrap(zp + dt_beta*dzp, AC_xyz0__mod__cdata.z, AC_xyz1__mod__cdata.z, AC_lxyz__mod__cdata.z)
		AC_fp_xp[k] = xn
		AC_fp_yp[k] = yn
		AC_fp_zp[k] = zn

		real leaving = 0.0
		if ((xn-AC_xyz0__mod__cdata.x) >= xx1_x || (xn-AC_xyz0__mod__cdata.x) < xx0_x) leaving = 1.0
		if ((yn-AC_xyz0__mod__cdata.y) >= xx1_y || (yn-AC_xyz0__mod__cdata.y) < xx0_y) leaving = 1.0
		if ((zn-AC_xyz0__mod__cdata.z) >= xx1_z || (zn-AC_xyz0__mod__cdata.z) < xx0_z) leaving = 1.0
		AC_par_leaving[k] = leaving
	}
}

//Clamp list entries to valid indices (they are garbage during the kernel analysis)
particle_list_index(real val)
{
	int i = int(val)
	if (i < 0) i = 0
	if (i > AC_npar-1) i = AC_npar-1
	return i
}

Kernel particles_pack()
{
	const int i = vertexIdx.x
	if (i >= 0)
	{
		const int k = particle_list_index(AC_par_idx[i])
		//Only particles flagged as leaving are listed; the check also makes AC_par_leaving
		//read by a kernel, otherwise acc can consider it as dead since it is only read by the host
		if (AC_par_leaving[k] != 0.0)
		{
		AC_par_buf[0][i] = AC_fp_xp[k]
		AC_par_buf[1][i] = AC_fp_yp[k]
		AC_par_buf[2][i] = AC_fp_zp[k]
		AC_par_buf[3][i] = AC_fp_vpx[k]
		AC_par_buf[4][i] = AC_fp_vpy[k]
		AC_par_buf[5][i] = AC_fp_vpz[k]
		AC_par_buf[6][i] = AC_dfp_xp[k]
		AC_par_buf[7][i] = AC_dfp_yp[k]
		AC_par_buf[8][i] = AC_dfp_zp[k]
		AC_par_buf[9][i] = AC_fp_ipar[k]
		}
	}
}

//Sources are all >= the new npar_loc and destinations all below it, so there are no conflicts
Kernel particles_move()
{
	const int i = vertexIdx.x
	if (i >= 0)
	{
		const int src = particle_list_index(AC_par_idx[i])
		const int dst = particle_list_index(AC_par_idx2[i])
		AC_fp_xp[dst]   = AC_fp_xp[src]
		AC_fp_yp[dst]   = AC_fp_yp[src]
		AC_fp_zp[dst]   = AC_fp_zp[src]
		AC_fp_vpx[dst]  = AC_fp_vpx[src]
		AC_fp_vpy[dst]  = AC_fp_vpy[src]
		AC_fp_vpz[dst]  = AC_fp_vpz[src]
		AC_dfp_xp[dst]  = AC_dfp_xp[src]
		AC_dfp_yp[dst]  = AC_dfp_yp[src]
		AC_dfp_zp[dst]  = AC_dfp_zp[src]
		AC_fp_ipar[dst] = AC_fp_ipar[src]
	}
}

Kernel particles_unpack(int offset)
{
	const int i = vertexIdx.x
	if (i >= 0)
	{
		int k = offset+i
		if (k > AC_npar-1) k = AC_npar-1
		AC_fp_xp[k]   = AC_par_buf[0][i]
		AC_fp_yp[k]   = AC_par_buf[1][i]
		AC_fp_zp[k]   = AC_par_buf[2][i]
		AC_fp_vpx[k]  = AC_par_buf[3][i]
		AC_fp_vpy[k]  = AC_par_buf[4][i]
		AC_fp_vpz[k]  = AC_par_buf[5][i]
		AC_dfp_xp[k]  = AC_par_buf[6][i]
		AC_dfp_yp[k]  = AC_par_buf[7][i]
		AC_dfp_zp[k]  = AC_par_buf[8][i]
		AC_fp_ipar[k] = AC_par_buf[9][i]
	}
}
#endif
