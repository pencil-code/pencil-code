real get_deltacorr_force(int step_num) {

  real3 force = real3(0.,0.,0.)
#ifdef LFORCING
  if ( !lforcing_cont ) {
    if (step_num==2) {force = forcing()}
  }
#endif
  return force
}
