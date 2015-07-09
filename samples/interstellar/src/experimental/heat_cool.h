!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_heat_cool, initialize_heat_cool
  public :: rprint_heat_cool
  public :: init_heat_cool
  public :: get_heat_coolstellar
  public :: read_heat_cool_init_pars, write_heat_cool_init_pars
  public :: read_heat_cool_run_pars,  write_heat_cool_run_pars
  public :: pencil_criteria_heat_cool
  public :: interstellar_before_boundary
  public :: calc_heat_cool_interstellar, check_SN
  public :: calc_snr_damping
  public :: calc_snr_damp_int
  public :: calc_snr_unshock
  public :: input_persistent_interstellar, output_persistent_interstellar
  public :: addmassflux

