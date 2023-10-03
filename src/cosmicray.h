!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_cosmicray, initialize_cosmicray
  public :: read_cosmicray_init_pars, write_cosmicray_init_pars
  public :: read_cosmicray_run_pars,  write_cosmicray_run_pars
  public :: rprint_cosmicray, get_slices_cosmicray
  public :: init_ecr, decr_dt
  public :: pencil_criteria_cosmicray, pencil_interdep_cosmicray
  public :: calc_pencils_cosmicray, calc_diagnostics_cosmicray
  public :: impose_ecr_floor

