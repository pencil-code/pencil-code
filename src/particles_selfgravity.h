!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_particles_selfgrav
  public :: initialize_particles_selfgrav
  public :: calc_selfpotential_particles
  public :: pencil_criteria_par_selfgrav, pencil_interdep_par_selfgrav
  public :: calc_pencils_par_selfgrav
  public :: dvvp_dt_selfgrav_pencil, dvvp_dt_selfgrav
  public :: read_particles_selfg_init_pars
  public :: write_particles_selfg_init_pars
  public :: read_particles_selfg_run_pars
  public :: write_particles_selfg_run_pars
  public :: rprint_particles_selfgrav
  public :: calc_diagnostics_particles_selg
