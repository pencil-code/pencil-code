real get_deltacorr_force(int step_num) {

  real3 force=real3(0.,0.,0.)

  if ( !lforcing_cont ) {
    if (step_num==2) {force=forcing()}
  }
  return force
}
