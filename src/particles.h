!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_particles, initialize_particles, init_particles
  public :: pencil_criteria_particles, pencil_interdep_particles
  public :: calc_pencils_particles
  public :: dxxp_dt_pencil, dvvp_dt_pencil, dxxp_dt, dvvp_dt
  public :: remove_particles_sink, create_sink_particles
  public :: powersnap_particles, rprint_particles
  public :: read_particles_init_pars, write_particles_init_pars
  public :: read_particles_run_pars, write_particles_run_pars
  public :: insert_particles
