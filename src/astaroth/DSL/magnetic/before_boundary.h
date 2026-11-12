#if Lmagnetic_MODULE
Profile<XY> AX_mean_z
Profile<XY> AY_mean_z
Profile<XY> AZ_mean_z
Kernel magnetic_before_boundary_reductions(bool lrmv)
{
      if (AC_lquench_eta_aniso__mod__magnetic)
      {
	 reduce_rms(F_AVEC,AC_Arms)
      }
      if (AC_lremove_meanaxy__mod__magnetic && lrmv)
      {
 	 fact=1./AC_nzgrid_eff__mod__cdata
	 reduce_sum(fact*F_AVEC.x,AX_mean_z)
	 reduce_sum(fact*F_AVEC.y,AY_mean_z)
	 reduce_sum(fact*F_AVEC.z,AZ_mean_z)
      }
}
Kernel magnetic_before_boundary(bool lrmv)
{
      if (AC_lremove_meanaxy__mod__magnetic && lrmv)
      {
	 F_AX[vertexIdx.x][vertexIdx.y][vertexIdx.z] = F_AX - AC_tau_remove_meanaxy__mod__magnetic*AX_mean_z
	 F_AY[vertexIdx.x][vertexIdx.y][vertexIdx.z] = F_AY - AC_tau_remove_meanaxy__mod__magnetic*AY_mean_z
	 F_AZ[vertexIdx.x][vertexIdx.y][vertexIdx.z] = F_AZ - AC_tau_remove_meanaxy__mod__magnetic*AZ_mean_z
      }
}
//TP: bb as communicated auxiliary (magnetic_before_boundary), read by the rhs through value(F_BVEC)
Kernel magnetic_calc_bb_comaux()
{
      if (AC_lbb_as_comaux__mod__magnetic)
      {
	 bb = curl(F_AVEC)
	 if (AC_lb_ext_in_comaux__mod__magnetic)
	 {
	    bb += AC_b_ext__mod__magnetic
	 }
	 F_BX[vertexIdx.x][vertexIdx.y][vertexIdx.z] = bb.x
	 F_BY[vertexIdx.x][vertexIdx.y][vertexIdx.z] = bb.y
	 F_BZ[vertexIdx.x][vertexIdx.y][vertexIdx.z] = bb.z
      }
}
#else
Kernel magnetic_calc_bb_comaux()
{
}
Kernel magnetic_before_boundary_reductions(bool lrmv)
{
	 suppress_unused_warning(lrmv)
}
Kernel magnetic_before_boundary(bool lrmv)
{
	 suppress_unused_warning(lrmv)
}
#endif
#if Lbfield_MODULE
Kernel bfield_get_j()
{
	jj = AC_mu01__mod__cdata*curl(F_BVEC)
	write(F_JVEC, jj)
}
Kernel bfield_get_e()
{
	ee = -cross(UU,F_BVEC + AC_b_ext__mod__magnetic)
	write(F_EVEC, ee)
}
#else
Kernel bfield_get_j()
{
}
Kernel bfield_get_e()
{
}
#endif
#if LSOLAR_CORONA && Lmagnetic_MODULE
Kernel solar_corona_mag_driver(real dt_driver)
{
	const int ix = vertexIdx.x-NGHOST_VAL
	const int iy = vertexIdx.y-NGHOST_VAL
	const real ainit_x = AC_a_init__mod__solar_corona[ix][iy][0]
	const real ainit_y = AC_a_init__mod__solar_corona[ix][iy][1]
	real ax = F_AX[vertexIdx.x][vertexIdx.y][vertexIdx.z]
	real ay = F_AY[vertexIdx.x][vertexIdx.y][vertexIdx.z]
	if (AC_domain_coordinates.z == 0)
	{
		const real fac = dt_driver*AC_b_tau__mod__solar_corona
		ax = ax*(1.0-fac) + ainit_x*fac
		ay = ay*(1.0-fac) + ainit_y*fac
	}
	F_AX[vertexIdx.x][vertexIdx.y][vertexIdx.z] = ax
	F_AY[vertexIdx.x][vertexIdx.y][vertexIdx.z] = ay
}
#endif
#if LCORONAE
//Bottom-plane drivers of special/coronae: the granular velocities (lgranulation, set at the first substep) and the
//relaxation of A towards the magnetogram (lmag_time_bound). The host computes them in special_before_boundary_shared.
//All processes launch the kernel on their lowest local plane so that they have the same task graphs; the kernel
//changes the values only on the lower boundary and elsewhere writes its input back.
Kernel coronae_bottom_driver(PC_SUB_STEP_NUMBER step_num, real dt_driver)
{
	const int ix = vertexIdx.x-NGHOST_VAL
	const int iy = vertexIdx.y-NGHOST_VAL
	const bool lbottom = AC_domain_coordinates.z == 0
	const real ux_gran = AC_ux__mod__coronae[ix][iy]*AC_u_amplifier__mod__coronae
	const real uy_gran = AC_uy__mod__coronae[ix][iy]*AC_u_amplifier__mod__coronae
	real ux = F_UX[vertexIdx.x][vertexIdx.y][vertexIdx.z]
	real uy = F_UY[vertexIdx.x][vertexIdx.y][vertexIdx.z]
	real uz = F_UZ[vertexIdx.x][vertexIdx.y][vertexIdx.z]
	if (lbottom && AC_lgranulation__mod__coronae && step_num == PC_FIRST_SUB_STEP)
	{
		ux = ux_gran
		uy = uy_gran
		uz = 0.0
	}
	F_UX[vertexIdx.x][vertexIdx.y][vertexIdx.z] = ux
	F_UY[vertexIdx.x][vertexIdx.y][vertexIdx.z] = uy
	F_UZ[vertexIdx.x][vertexIdx.y][vertexIdx.z] = uz
#if Lmagnetic_MODULE
	const real ax_mag = AC_ax_bound__mod__coronae[ix][iy]
	const real ay_mag = AC_ay_bound__mod__coronae[ix][iy]
	const real az_mag = AC_az_bound__mod__coronae[ix][iy]
	real ax = F_AX[vertexIdx.x][vertexIdx.y][vertexIdx.z]
	real ay = F_AY[vertexIdx.x][vertexIdx.y][vertexIdx.z]
	real az = F_AZ[vertexIdx.x][vertexIdx.y][vertexIdx.z]
	if (lbottom && AC_lmag_time_bound__mod__coronae)
	{
		const real b_tau = AC_b_tau__mod__coronae
		if (b_tau > 0.0)
		{
			ax = ax*(1.0-dt_driver*b_tau) + ax_mag*dt_driver*b_tau
			ay = ay*(1.0-dt_driver*b_tau) + ay_mag*dt_driver*b_tau
			if (AC_lmag_bound_az__mod__coronae) az = az*(1.0-dt_driver*b_tau) + az_mag*dt_driver*b_tau
		}
		else
		{
			ax = ax_mag
			ay = ay_mag
			if (AC_lmag_bound_az__mod__coronae) az = az_mag
		}
	}
	F_AX[vertexIdx.x][vertexIdx.y][vertexIdx.z] = ax
	F_AY[vertexIdx.x][vertexIdx.y][vertexIdx.z] = ay
	F_AZ[vertexIdx.x][vertexIdx.y][vertexIdx.z] = az
#endif
}
#endif
