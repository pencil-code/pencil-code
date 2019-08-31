!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_density, initialize_density
  public :: read_density_init_pars, write_density_init_pars
  public :: read_density_run_pars,  write_density_run_pars
  public :: rprint_density, get_slices_density, get_slices_pressure
  public :: init_lnrho, dlnrho_dt, impose_density_floor, impose_density_ceiling
  public :: density_before_boundary
  public :: split_update_density

  public :: pencil_criteria_density, pencil_interdep_density
  public :: calc_pencils_density
  public :: get_init_average_density, anelastic_after_mn,density_after_boundary
  public :: dynamical_diffusion, boussinesq
  public :: mean_density
  public :: update_char_vel_density
  public :: density_after_timestep
  public :: pushpars2c, pushdiags2c
  public :: calc_diagnostics_density
!
! WL: ONLY SUBROUTINES SHOULD BE PUBLIC. THESE DO NOT QUALIFY!!!!!
!
  public :: lnrhomz,glnrhomz,lcalc_lnrhomean,lcalc_glnrhomean,lupw_lnrho
  public :: beta_glnrho_global,beta_glnrho_scaled
