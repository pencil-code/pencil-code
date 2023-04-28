!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_neutralvelocity, initialize_neutralvelocity
  public :: read_neutralvelocity_init_pars, write_neutralvelocity_init_pars
  public :: read_neutralvelocity_run_pars,  write_neutralvelocity_run_pars
  public :: rprint_neutralvelocity
  public :: init_uun, calc_pencils_neutralvelocity, duun_dt, calc_diagnostics_neutralvel
  public :: pencil_criteria_neutralvelocity, pencil_interdep_neutralvelocity
