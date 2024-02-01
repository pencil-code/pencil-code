!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_selfgravity, initialize_selfgravity
  public :: pencil_criteria_selfgravity, pencil_interdep_selfgravity
  public :: calc_pencils_selfgravity, calc_selfpotential
  public :: addselfgrav, calc_diagnostics_selfgrav
  public :: read_selfgravity_init_pars, write_selfgravity_init_pars
  public :: read_selfgravity_run_pars, write_selfgravity_run_pars
  public :: rprint_selfgravity
  public :: rhs_poisson_const
