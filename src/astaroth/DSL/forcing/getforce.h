get_deltacorr_force(int step_num, real t) {
  real3 force = real3(0.,0.,0.)
  suppress_unused_warning(step_num)
#ifdef LFORCING
  if ( !AC_lforcing_cont__mod__cdata ) {
    if (step_num==AC_num_substeps__mod__cdata-1) {force = forcing()}
  }
#endif
  if (AC_lforce_ramp_down__mod__forcing)
  {
    tmp= max(0.0,1.0+min(0.0,(t-AC_t__mod__cdata)/AC_tauforce_ramp_down__mod__forcing))
    force *= tmp
  }
  return force
}
