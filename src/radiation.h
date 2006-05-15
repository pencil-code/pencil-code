!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_radiation, initialize_radiation
  public :: read_radiation_init_pars, write_radiation_init_pars
  public :: read_radiation_run_pars,  write_radiation_run_pars
  public :: rprint_radiation
  public :: pencil_criteria_radiation, pencil_interdep_radiation
  public :: calc_pencils_radiation

  public :: init_rad, radtransfer, de_dt, radiative_cooling, radiative_pressure

  public :: bc_ee_inflow_x, bc_ee_outflow_x
!ajwm SHOULDN'T BE EXPORTED
  public :: DFF_new
                                                                                                       


