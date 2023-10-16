!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: initialize_particles_surf
  public :: register_particles_surfspec
  public :: read_particles_surf_init_pars
  public :: write_particles_surf_init_pars
  public :: read_particles_surf_run_pars
  public :: write_particles_surf_run_pars
  public :: rprint_particles_surf
  public :: init_particles_surf
  public :: dpsurf_dt, dpsurf_dt_pencil
  public :: calc_psurf_pencils
  public :: cleanup_surf_pencils
  public :: particles_surfspec_clean_up
  public :: calc_diagnostics_particles_surf
