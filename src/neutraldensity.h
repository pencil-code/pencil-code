!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_neutraldensity, initialize_neutraldensity
  public :: read_neutraldensity_init_pars, write_neutraldensity_init_pars
  public :: read_neutraldensity_run_pars,  write_neutraldensity_run_pars
  public :: rprint_neutraldensity
  public :: init_lnrhon, dlnrhon_dt

  public :: pencil_criteria_neutraldensity, pencil_interdep_neutraldensity
  public :: calc_pencils_neutraldensity, calc_diagnostics_neutraldens
  public :: neutraldensity_after_boundary
