!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
  private
!
  public :: register_shock, initialize_shock, rprint_shock
  public :: get_slices_shock
  public :: calc_shock_profile, calc_shock_profile_simple
  public :: shock_before_boundary
  public :: read_shock_run_pars,  write_shock_run_pars
  public :: pencil_criteria_shock, pencil_interdep_shock
  public :: calc_pencils_shock, calc_diagnostics_shock
  public :: pushpars2c
!
