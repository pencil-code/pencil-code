!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_particles_caustics
  public :: initialize_particles_caustics,init_particles_caustics,reinitialize_caustics
  public :: read_pcaustics_init_pars, write_pcaustics_init_pars
  public :: read_pcaustics_run_pars, write_pcaustics_run_pars
  public :: rprint_particles_caustics
  public :: dcaustics_dt,dcaustics_dt_pencil
  public :: reset_caustics
