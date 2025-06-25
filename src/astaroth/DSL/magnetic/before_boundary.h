#if LMAGNETIC
Kernel magnetic_before_boundary_reductions()
{
      if (AC_lquench_eta_aniso__mod__magnetic)
      {
	reduce_rms(AA,AC_Arms)
      }
}
#else
Kernel magnetic_before_boundary_reductions()
{
}
#endif
