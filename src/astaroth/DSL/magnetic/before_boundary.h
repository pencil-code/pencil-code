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
#else
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
