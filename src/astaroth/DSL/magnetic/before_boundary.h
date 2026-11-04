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
