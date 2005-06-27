!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_cosmicrayflux, initialize_cosmicrayflux
  public :: read_cosmicrayflux_init_pars, write_cosmicrayflux_init_pars
  public :: read_cosmicrayflux_run_pars,  write_cosmicrayflux_run_pars
  public :: rprint_cosmicrayflux
  public :: init_fcr, dfcr_dt

  public :: pencil_criteria_cosmicrayflux, pencil_interdep_cosmicrayflux
  public :: calc_pencils_cosmicrayflux

