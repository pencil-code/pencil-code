!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_particles, initialize_particles, init_particles
  public :: pencil_criteria_particles, pencil_interdep_particles
  public :: calc_pencils_particles, particles_dragforce_stiff
  public :: dxxp_dt, dvvp_dt, dxxp_dt_pencil, dvvp_dt_pencil
  public :: dxxp_dt_blocks, dvvp_dt_blocks
  public :: remove_particles_sink_simple, create_particles_sink_simple
  public :: powersnap_particles, rprint_particles
  public :: read_particles_init_pars, write_particles_init_pars
  public :: read_particles_run_pars, write_particles_run_pars
  public :: insert_particles
