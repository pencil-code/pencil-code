!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_testscalar, initialize_testscalar
  public :: read_testscalar_init_pars, write_testscalar_init_pars
  public :: read_testscalar_run_pars,  write_testscalar_run_pars
  public :: rprint_testscalar
  public :: get_slices_testscalar
  public :: init_cctest, dcctest_dt
  public :: pencil_criteria_testscalar, pencil_interdep_testscalar
  public :: testscalar_after_boundary
  public :: rescaling_testscalar

