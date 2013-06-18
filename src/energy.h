!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_energy, initialize_energy
  public :: read_energy_init_pars, write_energy_init_pars
  public :: read_energy_run_pars, write_energy_run_pars
  public :: rprint_energy, get_slices_energy
  public :: init_energy, denergy_dt, calc_lenergy_pars
  public :: pencil_criteria_energy, pencil_interdep_energy
  public :: calc_pencils_energy, fill_farray_pressure
  public :: impose_energy_floor
  public :: dynamical_thermal_diffusion
  public :: split_update_energy
  public :: expand_shands_energy
