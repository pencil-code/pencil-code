!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
  real :: w_sldchar_ene=0.1

  public :: register_energy, initialize_energy
  public :: read_energy_init_pars, write_energy_init_pars
  public :: read_energy_run_pars, write_energy_run_pars
  public :: rprint_energy, get_slices_energy
  public :: init_energy, denergy_dt, energy_after_boundary
  public :: pencil_criteria_energy, pencil_interdep_energy
  public :: calc_pencils_energy, fill_farray_pressure
  public :: impose_energy_floor, energy_before_boundary
  public :: dynamical_thermal_diffusion
  public :: split_update_energy
  public :: expand_shands_energy
  public :: update_char_vel_energy
  public :: energy_after_timestep
  public :: pushpars2c
  public :: calc_diagnostics_energy
