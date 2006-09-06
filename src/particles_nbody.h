!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_particles_nbody, initialize_particles_nbody, init_particles_nbody
  public :: pencil_criteria_par_nbody, pencil_interdep_par_nbody
  public :: calc_pencils_par_nbody
  public :: dvvp_dt_nbody_pencil,dvvp_dt_nbody
  public :: rprint_particles_nbody
  public :: read_particles_nbody_init_pars, write_particles_nbody_init_pars
  public :: read_particles_nbody_run_pars, write_particles_nbody_run_pars
  public :: get_particles_interdistances,share_sinkparticles
  
