!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_density, initialize_density
  public :: read_density_init_pars, write_density_init_pars
  public :: read_density_run_pars,  write_density_run_pars
  public :: rprint_density
  public :: init_lnrho, dlnrho_dt

  public :: bc_lnrho_temp_z, bc_lnrho_pressure_z
  public :: pencil_criteria_density, pencil_interdep_density
  public :: calc_pencils_density
                                                                                               
!ajwm SHOULD NOT BE PUBLIC
  public :: cs2cool,mpoly, mpoly0, mpoly1, mpoly2, isothtop
  public :: beta_dlnrhodr, beta_dlnrhodr_scaled
  public :: initlnrho
  public :: idiag_ekin, idiag_ekintot, idiag_rhom, idiag_rhomin, idiag_rhomax
  public :: b_ell, rbound

