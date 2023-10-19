!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_heatflux, initialize_heatflux
  public :: read_heatflux_run_pars,  write_heatflux_run_pars
  public :: rprint_heatflux
  public :: get_slices_heatflux
  public :: init_heatflux
  public :: dheatflux_dt
  public :: calc_pencils_heatflux, calc_diagnostics_heatflux
  public :: pencil_criteria_heatflux
  public :: pencil_interdep_heatflux
