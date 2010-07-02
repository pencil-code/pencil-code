!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_hydro, initialize_hydro
  public :: read_hydro_init_pars, write_hydro_init_pars
  public :: read_hydro_run_pars,  write_hydro_run_pars
  public :: rprint_hydro
  public :: get_slices_hydro
  public :: init_uu, duu_dt, calc_pencils_hydro
  public :: pencil_criteria_hydro, pencil_interdep_hydro
  public :: coriolis_cartesian
