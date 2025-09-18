get_deltacorr_force(int step_num) {
  real3 force = real3(0.,0.,0.)
#ifdef LFORCING
  if ( !AC_lforcing_cont__mod__cdata ) {
    if (step_num==2) {force = forcing()}
  }
#endif
  return force
}
