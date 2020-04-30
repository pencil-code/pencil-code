!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_chiral, initialize_chiral

  public :: read_chiral_init_pars, write_chiral_init_pars
  public :: read_chiral_run_pars,  write_chiral_run_pars
  public :: rprint_chiral, get_slices_chiral
  public :: init_chiral, dXY_chiral_dt
  public :: pencil_criteria_chiral, pencil_interdep_chiral
  public :: calc_pencils_chiral
  public :: chiral_before_boundary
  public :: iXX_chiral,iYY_chiral
