!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_dustdensity, initialize_dustdensity
  public :: read_dustdensity_init_pars, write_dustdensity_init_pars
  public :: read_dustdensity_run_pars,  write_dustdensity_run_pars
  public :: rprint_dustdensity, get_slices_dustdensity
  public :: init_nd, calc_pencils_dustdensity, dndmd_dt
  public :: dustdensity_after_boundary, dustdensity_before_boundary, calc_diagnostics_dustdensity
  public :: pencil_criteria_dustdensity, pencil_interdep_dustdensity

  public :: impose_dustdensity_floor
!  public :: Ntot_i
