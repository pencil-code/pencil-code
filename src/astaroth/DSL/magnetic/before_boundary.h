#if Lmagnetic_MODULE
Profile<XY> AX_mean_z
Profile<XY> AY_mean_z
Profile<XY> AZ_mean_z
Kernel magnetic_before_boundary_reductions(bool lrmv)
{
      if (AC_lquench_eta_aniso__mod__magnetic)
      {
	 reduce_rms(AA,AC_Arms)
      }
      if (AC_lremove_meanaxy__mod__magnetic && lrmv)
      {
 	 fact=1./AC_nzgrid_eff__mod__cdata
	 reduce_sum(fact*AA.x,AX_mean_z)
	 reduce_sum(fact*AA.y,AY_mean_z)
	 reduce_sum(fact*AA.z,AZ_mean_z)
      }
}
Kernel magnetic_before_boundary(bool lrmv)
{
      if (AC_lremove_meanaxy__mod__magnetic && lrmv)
      {
	 AAX[vertexIdx.x][vertexIdx.y][vertexIdx.z] = AAX - AC_tau_remove_meanaxy__mod__magnetic*AX_mean_z
	 AAY[vertexIdx.x][vertexIdx.y][vertexIdx.z] = AAY - AC_tau_remove_meanaxy__mod__magnetic*AY_mean_z
	 AAZ[vertexIdx.x][vertexIdx.y][vertexIdx.z] = AAZ - AC_tau_remove_meanaxy__mod__magnetic*AZ_mean_z
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
