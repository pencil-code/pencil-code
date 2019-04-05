!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_testfield, initialize_testfield
  public :: read_testfield_init_pars, write_testfield_init_pars
  public :: read_testfield_run_pars,  write_testfield_run_pars
  public :: rprint_testfield
  public :: get_slices_testfield
  public :: init_aatest, daatest_dt
  public :: pencil_criteria_testfield, pencil_interdep_testfield
  public :: testfield_before_boundary, testfield_after_boundary
  public :: rescaling_testfield
