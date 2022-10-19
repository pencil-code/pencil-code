!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
  private
!
  public :: register_viscosity, initialize_viscosity, rprint_viscosity
  public :: calc_viscosity, calc_viscous_heat, calc_viscous_force
  public :: read_viscosity_run_pars,  write_viscosity_run_pars
  public :: pencil_criteria_viscosity, pencil_interdep_viscosity
  public :: calc_pencils_viscosity
  public :: viscosity_after_boundary
  public :: calc_visc_heat_ppd, getnu
  public :: dynamical_viscosity
  public :: split_update_viscosity
  public :: pushpars2c
!
!ajwm SHOULDN'T BE SHARED
  public :: lvisc_first
!
