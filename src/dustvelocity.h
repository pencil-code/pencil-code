!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_dustvelocity, initialize_dustvelocity
  public :: read_dustvelocity_init_pars, write_dustvelocity_init_pars
  public :: read_dustvelocity_run_pars,  write_dustvelocity_run_pars
  public :: rprint_dustvelocity, get_slices_dustvelocity
  public :: init_uud, calc_pencils_dustvelocity, duud_dt,calc_diagnostics_dustvelocity
  public :: pencil_criteria_dustvelocity, pencil_interdep_dustvelocity

  public :: copy_bcs_dust

! Shouldn't really be here
  public :: nd0, rhod0, ldustcoagulation, ldustcondensation
  public :: dust_binning
