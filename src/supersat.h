
!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_supersat, initialize_supersat

  public :: read_supersat_init_pars, write_supersat_init_pars
  public :: read_supersat_run_pars,  write_supersat_run_pars
  public :: rprint_supersat 
  public :: dssat_dt
  public :: init_ssat
  public :: pencil_criteria_supersat, pencil_interdep_supersat
  public :: calc_pencils_supersat
