
!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_ascalar, initialize_ascalar

  public :: read_ascalar_init_pars, write_ascalar_init_pars
  public :: read_ascalar_run_pars,  write_ascalar_run_pars
  public :: rprint_ascalar 
  public :: dacc_dt
  public :: init_acc
  public :: pencil_criteria_ascalar, pencil_interdep_ascalar
  public :: calc_pencils_ascalar, calc_diagnostics_ascalar
