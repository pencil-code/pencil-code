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
//The right-hand sides of the particle equations (kernel particles_rhs) are generated from dxxp_dt and dvvp_dt at
//compilation time, see transpile.py and fortran-parser/particles_rhs.py, for tracers and for the self-gravitating
//particles. particles_tracer_step is launched after the halos of the velocity have been updated.
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

//Right-hand sides from particles_rhs: xp,yp,zp,vpx,vpy,vpz
gmem real AC_par_rhs[6][AC_npar]
gmem real AC_par_leaving[AC_npar]
gmem real AC_par_idx[AC_npar]
gmem real AC_par_idx2[AC_npar]
//Particle records: column-major, so record i occupies AC_par_buf[0:PAR_NREC-1][i] contiguously
//Order: xp,yp,zp,vpx,vpy,vpz,dfp_xp,dfp_yp,dfp_zp,ipar (must agree with gpu_astaroth.cc)
#if LPARTICLES_DUST
//Dust records continue with: dfp_vpx,dfp_vpy,dfp_vpz,rhopswarm,aps,cell (only for the sink particles, not unpacked)
const int PAR_NREC = 16
gmem real AC_dfp_vpx[AC_npar]
gmem real AC_dfp_vpy[AC_npar]
gmem real AC_dfp_vpz[AC_npar]
gmem real AC_fp_rhopswarm[AC_npar]
gmem real AC_fp_aps[AC_npar]
gmem real AC_par_dt1[AC_npar]      //inverse time step of each particle, read by the host
gmem real AC_par_drag_x[AC_npar]
gmem real AC_par_drag_y[AC_npar]
gmem real AC_par_drag_z[AC_npar]
gmem real AC_par_tausp1_k[AC_npar]
gmem real AC_par_icell[AC_npar]    //cell of each particle, read by the host
gmem real AC_par_order[AC_npar]    //particles in cell order
//Sink particles: the list of all of them (from all processes), what the particles found out about them
const int PAR_SINK_CAP = 4096
gmem real AC_sink_x[PAR_SINK_CAP]
gmem real AC_sink_y[PAR_SINK_CAP]
gmem real AC_sink_z[PAR_SINK_CAP]
gmem real AC_sink_rad[PAR_SINK_CAP]
gmem real AC_par_accreted[AC_npar]   //1+index in the sink list of the sink accreting this particle, or 0
gmem real AC_par_interp[AC_npar]     //rhop at the particle, -1 for sink particles
gmem real AC_dep_interp[AC_npar]
gmem real AC_par_cell_best[AC_mlocal.x][AC_mlocal.y][AC_mlocal.z]
//Copies of the particle data in cell order
gmem real AC_dep_x[AC_npar]
gmem real AC_dep_y[AC_npar]
gmem real AC_dep_z[AC_npar]
gmem real AC_dep_fx[AC_npar]
gmem real AC_dep_fy[AC_npar]
gmem real AC_dep_fz[AC_npar]
gmem real AC_dep_rhop[AC_npar]
gmem real AC_dep_tausp1[AC_npar]
gmem real AC_par_cell_start[AC_mlocal.x][AC_mlocal.y][AC_mlocal.z]
gmem real AC_par_cell_count[AC_mlocal.x][AC_mlocal.y][AC_mlocal.z]
gmem real AC_par_dt1_cell[AC_mlocal.x][AC_mlocal.y][AC_mlocal.z]
//Only the surface of the subdomain goes through the host: the first ghost layer and the first interior layer
const int PAR_FOLD_CAP = (nxgrid/nprocx+2)*(nygrid/nprocy+2)*(nzgrid/nprocz+2) - (nxgrid/nprocx)*(nygrid/nprocy)*(nzgrid/nprocz)
gmem real AC_par_fold_idx[3][PAR_FOLD_CAP]
gmem real AC_par_fold_idx2[3][PAR_FOLD_CAP]
gmem real AC_par_fold_buf[PAR_NDEP][PAR_FOLD_CAP]
#else
const int PAR_NREC = 10
#endif
gmem real AC_par_buf[PAR_NREC][AC_npar]
#if LPARTICLES_NBODY
//Gradient of the potential, with the first ghost layer filled from the interior for the interpolation
gmem real AC_par_gpot[3][AC_mlocal.x][AC_mlocal.y][AC_mlocal.z]
#endif

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

//z index of the cell of a particle; a direction with one grid point has no ghost zones on the device
particle_cell_z(real zp)
{
	if (AC_nzgrid == 1) return AC_nmin.z
	return particle_lower_index(zp, AC_z__mod__cdata[NGHOST], AC_dz_1__mod__cdata[NGHOST], AC_nlocal.z)
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

//Interpolation primitives of the generated right-hand sides (interpolate_linear): as dxp0 etc. in that routine.
//The Fortran ineargrid is 1-based, as the grid arrays in the translated code.
particle_ineargrid_x(int k)
{
	if (AC_nxgrid == 1) return NGHOST+1
	return particle_lower_index(AC_fp_xp[k], AC_x__mod__cdata[NGHOST], AC_dx_1__mod__cdata[NGHOST], AC_nlocal.x)+1
}
particle_ineargrid_y(int k)
{
	if (AC_nygrid == 1) return NGHOST+1
	return particle_lower_index(AC_fp_yp[k], AC_y__mod__cdata[NGHOST], AC_dy_1__mod__cdata[NGHOST], AC_nlocal.y)+1
}
particle_ineargrid_z(int k)
{
	if (AC_nzgrid == 1) return NGHOST+1
	return particle_lower_index(AC_fp_zp[k], AC_z__mod__cdata[NGHOST], AC_dz_1__mod__cdata[NGHOST], AC_nlocal.z)+1
}
particle_ix0(real xp)
{
	int i = AC_nmin.x
	if (AC_nxgrid != 1)
	{
		i = particle_lower_index(xp, AC_x__mod__cdata[NGHOST], AC_dx_1__mod__cdata[NGHOST], AC_nlocal.x)
		if (AC_x__mod__cdata[i] > xp) i = i-1
	}
	return i
}
particle_iy0(real yp)
{
	int i = AC_nmin.y
	if (AC_nygrid != 1)
	{
		i = particle_lower_index(yp, AC_y__mod__cdata[NGHOST], AC_dy_1__mod__cdata[NGHOST], AC_nlocal.y)
		if (AC_y__mod__cdata[i] > yp) i = i-1
	}
	return i
}
particle_iz0(real zp)
{
	int i = AC_nmin.z
	if (AC_nzgrid != 1)
	{
		i = particle_lower_index(zp, AC_z__mod__cdata[NGHOST], AC_dz_1__mod__cdata[NGHOST], AC_nlocal.z)
		if (AC_z__mod__cdata[i] > zp) i = i-1
	}
	return i
}
particle_offset(real pos, real gridpos, int nglobal)
{
	if (nglobal != 1) return pos-gridpos
	return 0.0
}
#if LPARTICLES_TRACERS
particle_linear_gas_u(int k)
{
	const int ix0 = particle_ix0(AC_fp_xp[k])
	const int iy0 = particle_iy0(AC_fp_yp[k])
	const int iz0 = particle_iz0(AC_fp_zp[k])
	const real xp0 = particle_offset(AC_fp_xp[k], AC_x__mod__cdata[ix0], AC_nxgrid)
	const real yp0 = particle_offset(AC_fp_yp[k], AC_y__mod__cdata[iy0], AC_nygrid)
	const real zp0 = particle_offset(AC_fp_zp[k], AC_z__mod__cdata[iz0], AC_nzgrid)
	const real dx1 = AC_dx_1__mod__cdata[ix0]
	const real dy1 = AC_dy_1__mod__cdata[iy0]
	const real dz1 = AC_dz_1__mod__cdata[iz0]
	return real3(particle_interpolate_linear(F_UX, ix0, iy0, iz0, xp0, yp0, zp0, dx1, dy1, dz1),
	             particle_interpolate_linear(F_UY, ix0, iy0, iz0, xp0, yp0, zp0, dx1, dy1, dz1),
	             particle_interpolate_linear(F_UZ, ix0, iy0, iz0, xp0, yp0, zp0, dx1, dy1, dz1))
}
#include "particles_rhs.h"
#endif

//Advances the particles by one RK substep with the right-hand sides of particles_rhs and applies the periodic
//wrap and the migration flags, as particles_timestep_second + boundconds_particles do on the CPU.
//xx0_*, xx1_*: bounds of the local subdomain relative to xyz0 (xyz0_loc-xyz0, xyz1_loc-xyz0)
//as used in migrate_particles; the host passes -/+huge in directions without migration.
Kernel particles_advance(int step_num, real dt_step,
                         real xx0_x, real xx0_y, real xx0_z, real xx1_x, real xx1_y, real xx1_z)
{
	const real alpha   = particle_rk_alpha(step_num)
	const real dt_beta = dt_step*particle_rk_beta(step_num)
	const int k = vertexIdx.x
	if (k >= 0)   //always true; bare blocks are not allowed in the DSL
	{
		const real dxp = alpha*AC_dfp_xp[k] + AC_par_rhs[0][k]
		const real dyp = alpha*AC_dfp_yp[k] + AC_par_rhs[1][k]
		const real dzp = alpha*AC_dfp_zp[k] + AC_par_rhs[2][k]
		AC_dfp_xp[k] = dxp
		AC_dfp_yp[k] = dyp
		AC_dfp_zp[k] = dzp
#if LPARTICLES_NBODY
		const real dvpx = alpha*AC_dfp_vpx[k] + AC_par_rhs[3][k]
		const real dvpy = alpha*AC_dfp_vpy[k] + AC_par_rhs[4][k]
		const real dvpz = alpha*AC_dfp_vpz[k] + AC_par_rhs[5][k]
		AC_dfp_vpx[k] = dvpx
		AC_dfp_vpy[k] = dvpy
		AC_dfp_vpz[k] = dvpz
		AC_fp_vpx[k]  = AC_fp_vpx[k] + dt_beta*dvpx
		AC_fp_vpy[k]  = AC_fp_vpy[k] + dt_beta*dvpy
		AC_fp_vpz[k]  = AC_fp_vpz[k] + dt_beta*dvpz
#endif
		const real xn = particle_periodic_wrap(AC_fp_xp[k] + dt_beta*dxp, AC_xyz0__mod__cdata.x, AC_xyz1__mod__cdata.x, AC_lxyz__mod__cdata.x)
		const real yn = particle_periodic_wrap(AC_fp_yp[k] + dt_beta*dyp, AC_xyz0__mod__cdata.y, AC_xyz1__mod__cdata.y, AC_lxyz__mod__cdata.y)
		const real zn = particle_periodic_wrap(AC_fp_zp[k] + dt_beta*dzp, AC_xyz0__mod__cdata.z, AC_xyz1__mod__cdata.z, AC_lxyz__mod__cdata.z)
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

#if LPARTICLES_DUST
#if LPARTICLES_NBODY
//Self-gravitating particles without gas (particles_dust + particles_selfgravity), one process only.
//The kernels are named as for the dust so that gpu_astaroth.cc can drive both in the same way.

//Fills the first ghost layer of the gradient from the interior (periodic, one process)
Kernel particles_gpot_wrap()
{
	const int ix = vertexIdx.x
	const int iy = vertexIdx.y
	const int iz = vertexIdx.z
	int wx = ix
	int wy = iy
	int wz = iz
	if (wx <  AC_nmin.x)             wx = wx + AC_nlocal.x
	if (wx >= AC_nmin.x+AC_nlocal.x) wx = wx - AC_nlocal.x
	if (wy <  AC_nmin.y)             wy = wy + AC_nlocal.y
	if (wy >= AC_nmin.y+AC_nlocal.y) wy = wy - AC_nlocal.y
	if (wz <  AC_nmin.z)             wz = wz + AC_nlocal.z
	if (wz >= AC_nmin.z+AC_nlocal.z) wz = wz - AC_nlocal.z
	if (wx != ix || wy != iy || wz != iz)
	{
		AC_par_gpot[0][ix][iy][iz] = AC_par_gpot[0][wx][wy][wz]
		AC_par_gpot[1][ix][iy][iz] = AC_par_gpot[1][wx][wy][wz]
		AC_par_gpot[2][ix][iy][iz] = AC_par_gpot[2][wx][wy][wz]
	}
}

//interpolate_linear of a component of the gradient
particle_gpot_linear(int c, int ix0, int iy0, int iz0, real xp0, real yp0, real zp0, real dx1, real dy1, real dz1)
{
	int iz1 = iz0
	if (AC_nzgrid != 1) iz1 = iz0+1
	g1 = AC_par_gpot[c][ix0  ][iy0  ][iz0]
	g2 = AC_par_gpot[c][ix0+1][iy0  ][iz0]
	g3 = AC_par_gpot[c][ix0  ][iy0+1][iz0]
	g4 = AC_par_gpot[c][ix0+1][iy0+1][iz0]
	g5 = AC_par_gpot[c][ix0  ][iy0  ][iz1]
	g6 = AC_par_gpot[c][ix0+1][iy0  ][iz1]
	g7 = AC_par_gpot[c][ix0  ][iy0+1][iz1]
	g8 = AC_par_gpot[c][ix0+1][iy0+1][iz1]
	return g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) +
	       xp0*yp0*dx1*dy1*(g1-g2-g3+g4) + xp0*zp0*dx1*dz1*(g1-g2-g5+g6) +
	       yp0*zp0*dy1*dz1*(g1-g3-g5+g7) +
	       xp0*yp0*zp0*dx1*dy1*dz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
}

//interpolate_linear of the gradient of the potential at a particle, primitive of the generated right-hand side
particle_linear_gpot(int k)
{
	const int ix0 = particle_ix0(AC_fp_xp[k])
	const int iy0 = particle_iy0(AC_fp_yp[k])
	const int iz0 = particle_iz0(AC_fp_zp[k])
	const real xp0 = particle_offset(AC_fp_xp[k], AC_x__mod__cdata[ix0], AC_nxgrid)
	const real yp0 = particle_offset(AC_fp_yp[k], AC_y__mod__cdata[iy0], AC_nygrid)
	const real zp0 = particle_offset(AC_fp_zp[k], AC_z__mod__cdata[iz0], AC_nzgrid)
	const real dx1 = AC_dx_1__mod__cdata[ix0]
	const real dy1 = AC_dy_1__mod__cdata[iy0]
	const real dz1 = AC_dz_1__mod__cdata[iz0]
	return real3(particle_gpot_linear(0, ix0, iy0, iz0, xp0, yp0, zp0, dx1, dy1, dz1),
	             particle_gpot_linear(1, ix0, iy0, iz0, xp0, yp0, zp0, dx1, dy1, dz1),
	             particle_gpot_linear(2, ix0, iy0, iz0, xp0, yp0, zp0, dx1, dy1, dz1))
}
#include "particles_rhs.h"
#else
//Quadratic spline (TSC) interpolation as in interpolate_quadratic_spline, 3D only.
//The arguments are the nearest grid point and the weights at offsets -1, 0, +1 in each direction.
particle_interpolate_tsc(Field fld, int ix0, int iy0, int iz0,
                         real wxm, real wx0, real wxp, real wym, real wy0, real wyp, real wzm, real wz0, real wzp)
{
	real res = 0.0
	for c in -1:2
	{
		real wz = wz0
		if (c == -1) wz = wzm
		if (c ==  1) wz = wzp
		for b in -1:2
		{
			real wy = wy0
			if (b == -1) wy = wym
			if (b ==  1) wy = wyp
			res = res + wz*wy*(wxm*fld[ix0-1][iy0+b][iz0+c] + wx0*fld[ix0][iy0+b][iz0+c] + wxp*fld[ix0+1][iy0+b][iz0+c])
		}
	}
	return res
}

//0-based species of a particle, as jspec in get_frictiontime
particle_species(real ipar)
{
	int j = (npar_species*(int(ipar)-1))/AC_npar
	if (j < 0) j = 0
	if (j > npar_species-1) j = npar_species-1
	return j
}

//Inverse friction time of a species, as get_frictiontime
particle_tausp1(int jspec)
{
	real res = AC_tausp1__mod__particles
	if (npar_species > 1) res = AC_tausp1_species__mod__particles[jspec]
	return res
}

//Quadratic spline interpolation of the gas velocity at a particle, primitive of the generated right-hand side
particle_tsc_gas_u(int k)
{
	const real xp  = AC_fp_xp[k]
	const real yp  = AC_fp_yp[k]
	const real zp  = AC_fp_zp[k]
	const int ix0 = particle_lower_index(xp, AC_x__mod__cdata[NGHOST], AC_dx_1__mod__cdata[NGHOST], AC_nlocal.x)
	const int iy0 = particle_lower_index(yp, AC_y__mod__cdata[NGHOST], AC_dy_1__mod__cdata[NGHOST], AC_nlocal.y)
	const int iz0 = particle_lower_index(zp, AC_z__mod__cdata[NGHOST], AC_dz_1__mod__cdata[NGHOST], AC_nlocal.z)
	const real dxp0 = (xp-AC_x__mod__cdata[ix0])*AC_dx_1__mod__cdata[ix0]
	const real dyp0 = (yp-AC_y__mod__cdata[iy0])*AC_dy_1__mod__cdata[iy0]
	const real dzp0 = (zp-AC_z__mod__cdata[iz0])*AC_dz_1__mod__cdata[iz0]
	const real wxm = 0.5*(0.5-dxp0)*(0.5-dxp0)
	const real wx0 = 0.75-dxp0*dxp0
	const real wxp = 0.5*(0.5+dxp0)*(0.5+dxp0)
	const real wym = 0.5*(0.5-dyp0)*(0.5-dyp0)
	const real wy0 = 0.75-dyp0*dyp0
	const real wyp = 0.5*(0.5+dyp0)*(0.5+dyp0)
	const real wzm = 0.5*(0.5-dzp0)*(0.5-dzp0)
	const real wz0 = 0.75-dzp0*dzp0
	const real wzp = 0.5*(0.5+dzp0)*(0.5+dzp0)
	return real3(particle_interpolate_tsc(F_UX, ix0, iy0, iz0, wxm, wx0, wxp, wym, wy0, wyp, wzm, wz0, wzp),
	             particle_interpolate_tsc(F_UY, ix0, iy0, iz0, wxm, wx0, wxp, wym, wy0, wyp, wzm, wz0, wzp),
	             particle_interpolate_tsc(F_UZ, ix0, iy0, iz0, wxm, wx0, wxp, wym, wy0, wyp, wzm, wz0, wzp))
}
#include "particles_rhs.h"

//RK update of a substep with the generated right-hand sides, and the wrap and migration flags as boundconds_particles
//for a 3D Cartesian shearing box. The shear offset used in the wrap is the one after this substep: the CPU wraps at the start of the next one.
Kernel particles_dust_step(int step_num, real dt_step,
                           real xx0_x, real xx0_y, real xx0_z, real xx1_x, real xx1_y, real xx1_z,
                           real delta_y_now)
{
	const real delta_y_rate = -AC_sshear__mod__cdata*AC_lxyz__mod__cdata.x
	const real alpha   = particle_rk_alpha(step_num)
	const real dt_beta = dt_step*particle_rk_beta(step_num)
	const int k = vertexIdx.x
	if (k >= 0)   //always true; bare blocks are not allowed in the DSL
	{
		const real xp  = AC_fp_xp[k]
		const real yp  = AC_fp_yp[k]
		const real zp  = AC_fp_zp[k]
		const real dxp  = alpha*AC_dfp_xp[k] + AC_par_rhs[0][k]
		const real dyp  = alpha*AC_dfp_yp[k] + AC_par_rhs[1][k]
		const real dzp  = alpha*AC_dfp_zp[k] + AC_par_rhs[2][k]
		const real dvpx = alpha*AC_dfp_vpx[k] + AC_par_rhs[3][k]
		const real dvpy = alpha*AC_dfp_vpy[k] + AC_par_rhs[4][k]
		const real dvpz = alpha*AC_dfp_vpz[k] + AC_par_rhs[5][k]
		AC_dfp_xp[k]  = dxp
		AC_dfp_yp[k]  = dyp
		AC_dfp_zp[k]  = dzp
		AC_dfp_vpx[k] = dvpx
		AC_dfp_vpy[k] = dvpy
		AC_dfp_vpz[k] = dvpz
		AC_fp_vpx[k]  = AC_fp_vpx[k] + dt_beta*dvpx
		AC_fp_vpy[k]  = AC_fp_vpy[k] + dt_beta*dvpy
		AC_fp_vpz[k]  = AC_fp_vpz[k] + dt_beta*dvpz

		//Shear-periodic wrap in x, as in boundconds_particles (delta_y advanced as in advance_shear)
		real delta_y_end = delta_y_now + delta_y_rate*dt_beta
		delta_y_end = delta_y_end - int(delta_y_end/AC_lxyz__mod__cdata.y)*AC_lxyz__mod__cdata.y
		real xn = xp + dt_beta*dxp
		real yn = yp + dt_beta*dyp
		if (xn <  AC_xyz0__mod__cdata.x)
		{
			xn = xn + AC_lxyz__mod__cdata.x
			yn = yn - delta_y_end
		}
		if (xn >= AC_xyz1__mod__cdata.x)
		{
			xn = xn - AC_lxyz__mod__cdata.x
			yn = yn + delta_y_end
		}
		yn = particle_periodic_wrap(yn, AC_xyz0__mod__cdata.y, AC_xyz1__mod__cdata.y, AC_lxyz__mod__cdata.y)
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
#endif

//Particle-to-grid deposition. A particle cannot write to the 27 grid points it influences (several particles
//would write to the same point), so instead every grid point gathers from the particles in the surrounding
//cells and writes only to itself. The particles are copied in cell order to AC_dep_*
//and AC_par_cell_start/count give the range of each cell. The order and the ranges come from a counting sort on the host.
//Reads inside loops are not seen by the kernel analysis when the loops are empty there, hence the particle_dep_touch calls.

//Flags the particles outside of the local subdomain, after the host has changed the particles
Kernel particles_dust_flag_leaving(real xx0_x, real xx0_y, real xx0_z, real xx1_x, real xx1_y, real xx1_z)
{
	const int k = vertexIdx.x
	if (k >= 0)
	{
		const real xn = AC_fp_xp[k]-AC_xyz0__mod__cdata.x
		const real yn = AC_fp_yp[k]-AC_xyz0__mod__cdata.y
		const real zn = AC_fp_zp[k]-AC_xyz0__mod__cdata.z
		real leaving = 0.0
		if (xn >= xx1_x || xn < xx0_x) leaving = 1.0
		if (yn >= xx1_y || yn < xx0_y) leaving = 1.0
		if (zn >= xx1_z || zn < xx0_z) leaving = 1.0
		AC_par_leaving[k] = leaving
	}
}

//Cell (nearest grid point) of each particle, the sort key
Kernel particles_cell_index()
{
	const int k = vertexIdx.x
	if (k >= 0)
	{
		const int ix0 = particle_lower_index(AC_fp_xp[k], AC_x__mod__cdata[NGHOST], AC_dx_1__mod__cdata[NGHOST], AC_nlocal.x)
		const int iy0 = particle_lower_index(AC_fp_yp[k], AC_y__mod__cdata[NGHOST], AC_dy_1__mod__cdata[NGHOST], AC_nlocal.y)
		const int iz0 = particle_cell_z(AC_fp_zp[k])
		AC_par_icell[k] = ix0 + AC_mlocal.x*(iy0 + AC_mlocal.y*iz0)
	}
}

Kernel particles_dep_copy()
{
	const int i = vertexIdx.x
	if (i >= 0)
	{
		const int k = particle_list_index(AC_par_order[i])
#if LPARTICLES_NBODY
		const real rhop = AC_rhop_swarm__mod__particles
#else
		const real rhop = AC_fp_rhopswarm[k]
#endif
		AC_dep_x[i]      = AC_fp_xp[k]
		AC_dep_y[i]      = AC_fp_yp[k]
		AC_dep_z[i]      = AC_fp_zp[k]
		AC_dep_fx[i]     = rhop*AC_par_drag_x[k]
		AC_dep_fy[i]     = rhop*AC_par_drag_y[k]
		AC_dep_fz[i]     = rhop*AC_par_drag_z[k]
		real rhop_map = rhop
		if (AC_fp_aps[k] > 0.0) rhop_map = 0.0   //sink particles are not part of rhop
		AC_dep_rhop[i]   = rhop_map
		AC_dep_tausp1[i] = AC_par_tausp1_k[k]
	}
}

//TSC weight of a particle at normalized distance d from a grid point, off: is it not the nearest grid point
particle_tsc_weight(real d, bool off)
{
	if (off) return 1.125 - 1.5*abs(d) + 0.5*d*d
	return 0.75 - d*d
}

particle_dep_touch()
{
	return AC_dep_x[0]+AC_dep_y[0]+AC_dep_z[0]+AC_dep_fx[0]+AC_dep_fy[0]+AC_dep_fz[0]+AC_dep_rhop[0]+AC_dep_tausp1[0]
}

#if LPARTICLES_NBODY
//CIC weight of a particle at normalized distance d from a grid point, in a direction with one grid point only the
//particle counts fully and only for the point of its own cell
particle_cic_weight(real d, int offset, bool active)
{
	if (!active)
	{
		if (offset == 0) return 1.0
		return 0.0
	}
	return max(0.0, 1.0-abs(d))
}

//rhop (CIC, folded afterwards) and np (nearest grid point), as map_xxp_grid
Kernel particles_gather_maps()
{
	const int ix = vertexIdx.x
	const int iy = vertexIdx.y
	const int iz = vertexIdx.z
	const real touch = particle_dep_touch()
	const real xv = AC_x__mod__cdata[ix]
	const real yv = AC_y__mod__cdata[iy]
	const real zv = AC_z__mod__cdata[iz]
	const real dx1 = AC_dx_1__mod__cdata[ix]
	const real dy1 = AC_dy_1__mod__cdata[iy]
	const real dz1 = AC_dz_1__mod__cdata[iz]
	real rhop = 0.0
	real np = 0.0
	for c in -1:2
	{
		for b in -1:2
		{
			for a in -1:2
			{
				const int cx = ix+a
				const int cy = iy+b
				const int cz = iz+c
				if (cx >= AC_nmin.x && cx < AC_nmin.x+AC_nlocal.x && cy >= AC_nmin.y && cy < AC_nmin.y+AC_nlocal.y && cz >= AC_nmin.z && cz < AC_nmin.z+AC_nlocal.z)
				{
					const int first = int(AC_par_cell_start[cx][cy][cz])
					const int count = int(AC_par_cell_count[cx][cy][cz])
					if (a == 0 && b == 0 && c == 0) np = count
					for i in first:first+count
					{
						rhop = rhop + AC_dep_rhop[i]*
						       particle_cic_weight((AC_dep_x[i]-xv)*dx1, a, AC_nxgrid != 1)*
						       particle_cic_weight((AC_dep_y[i]-yv)*dy1, b, AC_nygrid != 1)*
						       particle_cic_weight((AC_dep_z[i]-zv)*dz1, c, AC_nzgrid != 1)
					}
				}
			}
		}
	}
	AC_par_dep[3][ix][iy][iz] = rhop
	AC_par_dep[4][ix][iy][iz] = np
	if (touch < -1.0e300) AC_par_dep[3][ix][iy][iz] = 0.0   //never true, see particle_dep_touch
}

#else
//Drag force on the gas, launched over the interior and the first ghost layer (folded afterwards)
Kernel particles_gather_force()
{
	const int ix = vertexIdx.x
	const int iy = vertexIdx.y
	const int iz = vertexIdx.z
	const real touch = particle_dep_touch()
	const real xv = AC_x__mod__cdata[ix]
	const real yv = AC_y__mod__cdata[iy]
	const real zv = AC_z__mod__cdata[iz]
	const real dx1 = AC_dx_1__mod__cdata[ix]
	const real dy1 = AC_dy_1__mod__cdata[iy]
	const real dz1 = AC_dz_1__mod__cdata[iz]
	real fx = 0.0
	real fy = 0.0
	real fz = 0.0
	for c in -1:2
	{
		for b in -1:2
		{
			for a in -1:2
			{
				const int cx = ix+a
				const int cy = iy+b
				const int cz = iz+c
				if (cx >= NGHOST && cx < NGHOST+AC_nlocal.x && cy >= NGHOST && cy < NGHOST+AC_nlocal.y && cz >= NGHOST && cz < NGHOST+AC_nlocal.z)
				{
					const int first = int(AC_par_cell_start[cx][cy][cz])
					const int last  = first + int(AC_par_cell_count[cx][cy][cz])
					for i in first:last
					{
						const real w = particle_tsc_weight((AC_dep_x[i]-xv)*dx1, a != 0)*
						               particle_tsc_weight((AC_dep_y[i]-yv)*dy1, b != 0)*
						               particle_tsc_weight((AC_dep_z[i]-zv)*dz1, c != 0)
						fx = fx + w*AC_dep_fx[i]
						fy = fy + w*AC_dep_fy[i]
						fz = fz + w*AC_dep_fz[i]
					}
				}
			}
		}
	}
	const real rho1 = 1.0/F_RHO[ix][iy][iz]
	AC_par_dep[0][ix][iy][iz] = -rho1*fx
	AC_par_dep[1][ix][iy][iz] = -rho1*fy
	AC_par_dep[2][ix][iy][iz] = -rho1*fz
	if (touch < -1.0e300) AC_par_dep[0][ix][iy][iz] = 0.0   //never true, see particle_dep_touch
}

//rhop (TSC, folded afterwards) and np (nearest grid point), as map_xxp_grid. Only needed for output.
Kernel particles_gather_maps()
{
	const int ix = vertexIdx.x
	const int iy = vertexIdx.y
	const int iz = vertexIdx.z
	const real touch = particle_dep_touch()
	const real xv = AC_x__mod__cdata[ix]
	const real yv = AC_y__mod__cdata[iy]
	const real zv = AC_z__mod__cdata[iz]
	const real dx1 = AC_dx_1__mod__cdata[ix]
	const real dy1 = AC_dy_1__mod__cdata[iy]
	const real dz1 = AC_dz_1__mod__cdata[iz]
	real rhop = 0.0
	real np = 0.0
	for c in -1:2
	{
		for b in -1:2
		{
			for a in -1:2
			{
				const int cx = ix+a
				const int cy = iy+b
				const int cz = iz+c
				if (cx >= NGHOST && cx < NGHOST+AC_nlocal.x && cy >= NGHOST && cy < NGHOST+AC_nlocal.y && cz >= NGHOST && cz < NGHOST+AC_nlocal.z)
				{
					const int first = int(AC_par_cell_start[cx][cy][cz])
					const int count = int(AC_par_cell_count[cx][cy][cz])
					if (a == 0 && b == 0 && c == 0) np = count
					for i in first:first+count
					{
						rhop = rhop + AC_dep_rhop[i]*
						       particle_tsc_weight((AC_dep_x[i]-xv)*dx1, a != 0)*
						       particle_tsc_weight((AC_dep_y[i]-yv)*dy1, b != 0)*
						       particle_tsc_weight((AC_dep_z[i]-zv)*dz1, c != 0)
					}
				}
			}
		}
	}
	AC_par_dep[3][ix][iy][iz] = rhop
	AC_par_dep[4][ix][iy][iz] = np
	if (touch < -1.0e300) AC_par_dep[3][ix][iy][iz] = 0.0   //never true, see particle_dep_touch
}

//Inverse friction time step of each cell: dt1_drag in dvvp_dt_blocks. Launched over the interior.
Kernel particles_cell_dt1()
{
	const int ix = vertexIdx.x
	const int iy = vertexIdx.y
	const int iz = vertexIdx.z
	const real touch = particle_dep_touch()
	const int first = int(AC_par_cell_start[ix][iy][iz])
	const int last  = first + int(AC_par_cell_count[ix][iy][iz])
	real dust = 0.0
	real gas = 0.0
	for i in first:last
	{
		dust = max(dust, AC_dep_tausp1[i])
		gas  = gas + AC_dep_rhop[i]*AC_dep_tausp1[i]
	}
	real gas_term = 0.0
	if (AC_ldragforce_gas_par__mod__particles) gas_term = gas/F_RHO[ix][iy][iz]
	AC_par_dt1_cell[ix][iy][iz] = (dust + gas_term)/AC_cdtp__mod__particles
	if (touch < -1.0e300) AC_par_dt1_cell[ix][iy][iz] = 0.0   //never true, see particle_dep_touch
}

#endif

//Copies the listed cells of AC_par_dep to AC_par_fold_buf, for the host. No loops over the components here:
//they would be empty in the kernel analysis, which runs with zero parameters, and the arrays would be removed.
//Clears the first ghost layer of rhop after the fold, as fold_f: the interpolation of rhop at the particles sees
//zero there. Component 3 only: the others are not interpolated.
Kernel particles_fold_clear()
{
	const int i = vertexIdx.x
	if (i >= 0)
	{
		const int x = int(AC_par_fold_idx[0][i])
		const int y = int(AC_par_fold_idx[1][i])
		const int z = int(AC_par_fold_idx[2][i])
		AC_par_dep[3][x][y][z] = 0.0
	}
}

Kernel particles_fold_pack()
{
	const int i = vertexIdx.x
	if (i >= 0)
	{
		const int x = int(AC_par_fold_idx[0][i])
		const int y = int(AC_par_fold_idx[1][i])
		const int z = int(AC_par_fold_idx[2][i])
		AC_par_fold_buf[0][i] = AC_par_dep[0][x][y][z]
		AC_par_fold_buf[1][i] = AC_par_dep[1][x][y][z]
		AC_par_fold_buf[2][i] = AC_par_dep[2][x][y][z]
		AC_par_fold_buf[3][i] = AC_par_dep[3][x][y][z]
		AC_par_fold_buf[4][i] = AC_par_dep[4][x][y][z]
	}
}

particle_fold_mask(int c, int comp0, int ncomp)
{
	if (c >= comp0 && c < comp0+ncomp) return 1.0
	return 0.0
}

//Adds the contributions from the host to the components comp0:comp0+ncomp-1 of the listed cells
Kernel particles_fold_add(int comp0, int ncomp)
{
	const int i = vertexIdx.x
	if (i >= 0)
	{
		const int x = int(AC_par_fold_idx2[0][i])
		const int y = int(AC_par_fold_idx2[1][i])
		const int z = int(AC_par_fold_idx2[2][i])
		AC_par_dep[0][x][y][z] = AC_par_dep[0][x][y][z] + particle_fold_mask(0,comp0,ncomp)*AC_par_fold_buf[0][i]
		AC_par_dep[1][x][y][z] = AC_par_dep[1][x][y][z] + particle_fold_mask(1,comp0,ncomp)*AC_par_fold_buf[1][i]
		AC_par_dep[2][x][y][z] = AC_par_dep[2][x][y][z] + particle_fold_mask(2,comp0,ncomp)*AC_par_fold_buf[2][i]
		AC_par_dep[3][x][y][z] = AC_par_dep[3][x][y][z] + particle_fold_mask(3,comp0,ncomp)*AC_par_fold_buf[3][i]
		AC_par_dep[4][x][y][z] = AC_par_dep[4][x][y][z] + particle_fold_mask(4,comp0,ncomp)*AC_par_fold_buf[4][i]
	}
}

//Sink particles (particles_sink.f90). Flags the particles near a sink (within AC_sink_rad, which the host sets to
//twice the sink radius); the few flagged particles are then treated by the host as in sink_particle_accretion.
Kernel particles_sink_detect(int nsink, real delta_y_now)
{
	const int k = vertexIdx.x
	if (k >= 0)
	{
		const real touch = AC_sink_x[0]+AC_sink_y[0]+AC_sink_z[0]+AC_sink_rad[0]   //see particle_dep_touch
		const real lenx = AC_lxyz__mod__cdata.x
		const real leny = AC_lxyz__mod__cdata.y
		const real lenz = AC_lxyz__mod__cdata.z
		const real xp = AC_fp_xp[k]
		const real yp = AC_fp_yp[k]
		const real zp = AC_fp_zp[k]
		const bool is_sink = AC_fp_aps[k] != 0.0
		real hit = 0.0
		for j in 0:nsink
		{
			const real rad = AC_sink_rad[j]
			//Nearest image by directional splitting, ties as minloc over dis=(-1,0,1)
			real sepx = xp-AC_sink_x[j]+lenx
			int disx = -1
			if (abs(xp-AC_sink_x[j]) < abs(sepx))
			{
				sepx = xp-AC_sink_x[j]
				disx = 0
			}
			if (abs(xp-AC_sink_x[j]-lenx) < abs(sepx))
			{
				sepx = xp-AC_sink_x[j]-lenx
				disx = 1
			}
			real sy = AC_sink_y[j]
			if (disx == -1) sy = sy + delta_y_now
			if (disx ==  1) sy = sy - delta_y_now
			real sepy = yp-sy+leny
			if (abs(yp-sy) < abs(sepy)) sepy = yp-sy
			if (abs(yp-sy-leny) < abs(sepy)) sepy = yp-sy-leny
			real sepz = zp-AC_sink_z[j]+lenz
			if (abs(zp-AC_sink_z[j]) < abs(sepz)) sepz = zp-AC_sink_z[j]
			if (abs(zp-AC_sink_z[j]-lenz) < abs(sepz)) sepz = zp-AC_sink_z[j]-lenz
			if (!is_sink && rad > 0.0 && abs(sepx) <= rad && abs(sepy) <= rad && abs(sepz) <= rad && sepx*sepx+sepy*sepy+sepz*sepz <= rad*rad) hit = j+1
		}
		if (touch < -1.0e300) hit = 0.0
		AC_par_accreted[k] = hit
	}
}

//rhop at each particle as in create_particles_sink
Kernel particles_sink_interp()
{
	const int k = vertexIdx.x
	if (k >= 0)
	{
		const real xp = AC_fp_xp[k]
		const real yp = AC_fp_yp[k]
		const real zp = AC_fp_zp[k]
		const int ix0 = particle_lower_index(xp, AC_x__mod__cdata[NGHOST], AC_dx_1__mod__cdata[NGHOST], AC_nlocal.x)
		const int iy0 = particle_lower_index(yp, AC_y__mod__cdata[NGHOST], AC_dy_1__mod__cdata[NGHOST], AC_nlocal.y)
		const int iz0 = particle_lower_index(zp, AC_z__mod__cdata[NGHOST], AC_dz_1__mod__cdata[NGHOST], AC_nlocal.z)
		const real dxp0 = (xp-AC_x__mod__cdata[ix0])*AC_dx_1__mod__cdata[ix0]
		const real dyp0 = (yp-AC_y__mod__cdata[iy0])*AC_dy_1__mod__cdata[iy0]
		const real dzp0 = (zp-AC_z__mod__cdata[iz0])*AC_dz_1__mod__cdata[iz0]
		const real wxm = 0.5*(0.5-dxp0)*(0.5-dxp0)
		const real wx0 = 0.75-dxp0*dxp0
		const real wxp = 0.5*(0.5+dxp0)*(0.5+dxp0)
		const real wym = 0.5*(0.5-dyp0)*(0.5-dyp0)
		const real wy0 = 0.75-dyp0*dyp0
		const real wyp = 0.5*(0.5+dyp0)*(0.5+dyp0)
		const real wzm = 0.5*(0.5-dzp0)*(0.5-dzp0)
		const real wz0 = 0.75-dzp0*dzp0
		const real wzp = 0.5*(0.5+dzp0)*(0.5+dzp0)
		real res = 0.0
		for c in -1:2
		{
			real wz = wz0
			if (c == -1) wz = wzm
			if (c ==  1) wz = wzp
			for b in -1:2
			{
				real wy = wy0
				if (b == -1) wy = wym
				if (b ==  1) wy = wyp
				res = res + wz*wy*(wxm*AC_par_dep[3][ix0-1][iy0+b][iz0+c] + wx0*AC_par_dep[3][ix0][iy0+b][iz0+c] + wxp*AC_par_dep[3][ix0+1][iy0+b][iz0+c])
			}
		}
		if (AC_fp_aps[k] != 0.0) res = -1.0
		AC_par_interp[k] = res
	}
}

Kernel particles_sink_copy_interp()
{
	const int i = vertexIdx.x
	if (i >= 0)
	{
		const int k = particle_list_index(AC_par_order[i])
		AC_dep_interp[i] = AC_par_interp[k]
	}
}

//Largest rhop at the particles of each cell, for lsink_create_one_per_cell
Kernel particles_sink_cell_best()
{
	const int ix = vertexIdx.x
	const int iy = vertexIdx.y
	const int iz = vertexIdx.z
	const real touch = AC_dep_interp[0]   //see particle_dep_touch
	const int first = int(AC_par_cell_start[ix][iy][iz])
	const int last  = first + int(AC_par_cell_count[ix][iy][iz])
	real best = -1.0
	for i in first:last
	{
		best = max(best, AC_dep_interp[i])
	}
	if (touch < -1.0e300) best = 0.0
	AC_par_cell_best[ix][iy][iz] = best
}

#if LPARTICLES_SINK
//A particle becomes a sink if rhop reaches rhop_sink_create, with lsink_create_one_per_cell only the one with the
//largest rhop in its cell
Kernel particles_sink_create()
{
	const int k = vertexIdx.x
	if (k >= 0)
	{
		const int ix0 = particle_lower_index(AC_fp_xp[k], AC_x__mod__cdata[NGHOST], AC_dx_1__mod__cdata[NGHOST], AC_nlocal.x)
		const int iy0 = particle_lower_index(AC_fp_yp[k], AC_y__mod__cdata[NGHOST], AC_dy_1__mod__cdata[NGHOST], AC_nlocal.y)
		const int iz0 = particle_lower_index(AC_fp_zp[k], AC_z__mod__cdata[NGHOST], AC_dz_1__mod__cdata[NGHOST], AC_nlocal.z)
		const real interp = AC_par_interp[k]
		const real best = AC_par_cell_best[ix0][iy0][iz0]
		real aps = AC_fp_aps[k]
		if (aps == 0.0 && interp >= AC_rhop_sink_create__mod__particles_sink &&
		    (!AC_lsink_create_one_per_cell__mod__particles_sink || interp >= best)) aps = AC_sink_radius__mod__particles_sink
		AC_fp_aps[k] = aps
	}
}
#endif

//New state of the local sink particles listed in AC_par_idx, from the records in AC_par_buf
Kernel particles_sink_store()
{
	const int i = vertexIdx.x
	if (i >= 0)
	{
		const int k = particle_list_index(AC_par_idx[i])
		AC_fp_xp[k]  = AC_par_buf[0][i]
		AC_fp_yp[k]  = AC_par_buf[1][i]
		AC_fp_zp[k]  = AC_par_buf[2][i]
		AC_fp_vpx[k] = AC_par_buf[3][i]
		AC_fp_vpy[k] = AC_par_buf[4][i]
		AC_fp_vpz[k] = AC_par_buf[5][i]
		AC_fp_rhopswarm[k] = AC_par_buf[13][i]
	}
}

#if LPARTICLES_NBODY
#else
//Inverse time step of each particle: the advective limit of particles_rhs and the friction limit of its cell (particles_cell_dt1)
Kernel particles_dust_dt1()
{
	const int k = vertexIdx.x
	if (k >= 0)
	{
		const int ix0 = particle_lower_index(AC_fp_xp[k], AC_x__mod__cdata[NGHOST], AC_dx_1__mod__cdata[NGHOST], AC_nlocal.x)
		const int iy0 = particle_lower_index(AC_fp_yp[k], AC_y__mod__cdata[NGHOST], AC_dy_1__mod__cdata[NGHOST], AC_nlocal.y)
		const int iz0 = particle_lower_index(AC_fp_zp[k], AC_z__mod__cdata[NGHOST], AC_dz_1__mod__cdata[NGHOST], AC_nlocal.z)
		AC_par_dt1[k] = max(AC_par_dt1[k], AC_par_dt1_cell[ix0][iy0][iz0])
	}
}
#endif
#endif

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
		//acc can consider AC_par_leaving as dead since otherwise it is only read by the host
		if (i >= 0)   //always true; bare blocks are not allowed in the DSL
		{
		if (AC_par_leaving[k] < -1.0) AC_par_buf[0][i] = 0.0   //never true: makes AC_par_leaving read by a kernel
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
#if LPARTICLES_DUST
		AC_par_buf[10][i] = AC_dfp_vpx[k]
		AC_par_buf[11][i] = AC_dfp_vpy[k]
		AC_par_buf[12][i] = AC_dfp_vpz[k]
		AC_par_buf[13][i] = AC_fp_rhopswarm[k]
		AC_par_buf[14][i] = AC_fp_aps[k]
		AC_par_buf[15][i] = AC_par_icell[k]
		//Never true: makes AC_par_dt1 read by a kernel, otherwise acc considers it dead since only the host reads it
		if (AC_par_dt1[k] < -1.0) AC_par_buf[14][i] = 0.0
#endif
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
#if LPARTICLES_DUST
		AC_dfp_vpx[dst]      = AC_dfp_vpx[src]
		AC_dfp_vpy[dst]      = AC_dfp_vpy[src]
		AC_dfp_vpz[dst]      = AC_dfp_vpz[src]
		AC_fp_rhopswarm[dst] = AC_fp_rhopswarm[src]
		AC_fp_aps[dst]       = AC_fp_aps[src]
#endif
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
#if LPARTICLES_DUST
		AC_dfp_vpx[k]      = AC_par_buf[10][i]
		AC_dfp_vpy[k]      = AC_par_buf[11][i]
		AC_dfp_vpz[k]      = AC_par_buf[12][i]
		AC_fp_rhopswarm[k] = AC_par_buf[13][i]
		AC_fp_aps[k]       = AC_par_buf[14][i]
#endif
	}
}
#endif
