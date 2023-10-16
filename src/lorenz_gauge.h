!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_lorenz_gauge, initialize_lorenz_gauge
  public :: read_lorenz_gauge_init_pars, write_lorenz_gauge_init_pars
  public :: read_lorenz_gauge_run_pars,  write_lorenz_gauge_run_pars
  public :: rprint_lorenz_gauge
  public :: get_slices_lorenz_gauge
  public :: init_lorenz_gauge

  public :: dlorenz_gauge_dt

  public :: calc_pencils_lorenz_gauge, calc_diagnostics_lorenz_gauge
  public :: pencil_criteria_lorenz_gauge
  public :: pencil_interdep_lorenz_gauge
