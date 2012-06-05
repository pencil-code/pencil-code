!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_shear, initialize_shear
  public :: read_shear_init_pars, write_shear_init_pars
  public :: read_shear_run_pars, write_shear_run_pars
  public :: boundcond_shear, rprint_shear
  public :: shear_before_boundary, shearing, advance_shear
  public :: pencil_criteria_shear, pencil_interdep_shear
  public :: calc_pencils_shear
  public :: shear_variables
  public :: sheared_advection_fft
