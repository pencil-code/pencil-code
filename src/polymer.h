!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_polymer, initialize_polymer
  public :: read_polymer_init_pars, write_polymer_init_pars
  public :: read_polymer_run_pars,  write_polymer_run_pars
  public :: rprint_polymer
  public :: get_slices_polymer
  public :: init_poly
  public :: dpoly_dt
  public :: calc_pencils_polymer
  public :: pencil_criteria_polymer
  public :: pencil_interdep_polymer
  public :: calc_polymer_after_boundary
