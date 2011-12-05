!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_pscalar, initialize_pscalar

  public :: read_pscalar_init_pars, write_pscalar_init_pars
  public :: read_pscalar_run_pars,  write_pscalar_run_pars
  public :: rprint_pscalar, get_slices_pscalar
  public :: init_lncc, dlncc_dt
  public :: pencil_criteria_pscalar, pencil_interdep_pscalar
  public :: calc_pencils_pscalar
  public :: pscalar_after_boundary
  public :: calc_mpscalar

!ajwm SHOULD NOT BE PUBLIC
!axel: yes, but at the moment we need this to calculate PDFs.
  public :: rhoccm, idiag_rhoccm
  public :: cc2m, idiag_cc2m
  public :: gcc2m, idiag_gcc2m

