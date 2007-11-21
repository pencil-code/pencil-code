!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_density, initialize_density
  public :: read_density_init_pars, write_density_init_pars
  public :: read_density_run_pars,  write_density_run_pars
  public :: rprint_density
  public :: init_lnrho, dlnrho_dt, impose_density_floor

  public :: pencil_criteria_density, pencil_interdep_density
  public :: calc_pencils_density

