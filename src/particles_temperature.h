!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_particles_temperature, initialize_particles_temperature
  public :: init_particles_temperature
  public :: particles_temperature_prepencil_calc
  public :: pencil_criteria_par_temperature
  public :: dpTT_dt, dpTT_dt_pencil
  public :: read_particles_temperature_init_pars, write_particles_temperature_init_pars
  public :: read_particles_temperature_run_pars, write_particles_temperature_run_pars
  public :: rprint_particles_temperature

