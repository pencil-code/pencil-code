!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_particles_radius, initialize_particles_radius
  public :: pencil_criteria_par_radius
  public :: dap_dt, dap_dt_pencil
  public :: read_particles_rad_init_pars, write_particles_rad_init_pars
  public :: read_particles_rad_run_pars, write_particles_rad_run_pars
  public :: rprint_particles_radius, set_particle_radius
  public :: get_stbin,get_mass_from_radius,get_maxrad
  public :: calc_diagnostics_particles_rad
