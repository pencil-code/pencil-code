!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: particles_register_modules, particles_rprint_list
  public :: particles_initialize_modules, particles_init
  public :: particles_boundconds, particles_special
  public :: particles_timestep_first, particles_timestep_second
  public :: particles_pencil_criteria, particles_pencil_interdep
  public :: particles_calc_pencils, particles_calc_selfpotential
  public :: particles_read_snapshot, particles_calc_nbodygravity
  public :: particles_write_snapshot,particles_write_dsnapshot
  public :: particles_pde, particles_pde_pencil, particles_write_pdim
  public :: read_particles_init_pars_wrap
  public :: write_particles_init_pars_wrap
  public :: read_particles_run_pars_wrap
  public :: write_particles_run_pars_wrap
  public :: particles_powersnap, get_slices_particles
  public :: particles_doprepencil_calc
