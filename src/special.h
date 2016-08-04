!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_special
  public :: register_particles_special
  public :: initialize_special
  public :: finalize_special
  public :: read_special_init_pars
  public :: write_special_init_pars
  public :: read_special_run_pars
  public :: write_special_run_pars
  public :: rprint_special
  public :: get_slices_special
  public :: init_special

  public :: dspecial_dt

  public :: calc_pencils_special
  public :: pencil_criteria_special
  public :: pencil_interdep_special

  public :: calc_lspecial_pars

  public :: special_calc_hydro
  public :: special_calc_density
  public :: special_calc_dustdensity
  public :: special_calc_energy
  public :: special_calc_magnetic
  public :: special_calc_pscalar
  public :: special_calc_particles
  public :: special_calc_particles_nbody
  public :: special_calc_chemistry

  public :: special_boundconds
  public :: special_before_boundary
  public :: special_after_boundary
  public :: special_after_timestep
    
  public :: set_init_parameters   

