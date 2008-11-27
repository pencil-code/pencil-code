!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_chemistry, initialize_chemistry
  public :: read_chemistry_init_pars, write_chemistry_init_pars
  public :: read_chemistry_run_pars,  write_chemistry_run_pars
  public :: rprint_chemistry
  public :: get_slices_chemistry
  public :: init_chemistry

  public :: dchemistry_dt

  public :: calc_pencils_chemistry
  public :: pencil_criteria_chemistry
  public :: pencil_interdep_chemistry


  public :: calc_for_chem_mixture

  public :: get_p_infx
  public :: get_p_infy
  public :: get_rhs_Y
  public :: get_rhs_T
  public :: get_cs2_full
  public :: get_gamma_full

! public :: chemistry_calc_density
! public :: chemistry_calc_hydro
! public :: chemistry_calc_entropy
! public :: chemistry_calc_magnetic

! public :: chemistry_before_boundary
