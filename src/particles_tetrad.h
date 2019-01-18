!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_particles_tetrad
  public :: initialize_particles_tetrad,init_particles_tetrad,reinitialize_tetrad
  public :: read_ptetrad_init_pars, write_ptetrad_init_pars
  public :: read_ptetrad_run_pars, write_ptetrad_run_pars
  public :: rprint_particles_tetrad
  public :: dtetrad_dt,dtetrad_dt_pencil
