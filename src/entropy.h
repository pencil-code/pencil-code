!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_entropy, initialize_entropy
  public :: read_entropy_init_pars, write_entropy_init_pars
  public :: read_entropy_run_pars, write_entropy_run_pars
  public :: rprint_entropy, get_slices_entropy
  public :: init_ss, dss_dt, calc_lentropy_pars
  public :: pencil_criteria_entropy, pencil_interdep_entropy
  public :: calc_pencils_entropy, fill_farray_pressure
  public :: dynamical_thermal_diffusion
