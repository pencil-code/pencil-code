!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_testflow, initialize_testflow
  public :: read_testflow_init_pars, write_testflow_init_pars
  public :: read_testflow_run_pars,  write_testflow_run_pars
  public :: rprint_testflow
  public :: get_slices_testflow
  public :: init_uutest, duutest_dt
  public :: pencil_criteria_testflow, pencil_interdep_testflow
  public :: calc_ltestflow_nonlin_terms, testflow_before_boundary

