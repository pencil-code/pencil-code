!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)

  private

  public :: register_gravity, initialize_gravity
  public :: read_gravity_init_pars, write_gravity_init_pars
  public :: read_gravity_run_pars,  write_gravity_run_pars
  public :: rprint_gravity
  public :: init_gg, calc_pencils_gravity, duu_dt_grav
  public :: pencil_criteria_gravity,pencil_interdep_gravity
  public :: compute_gravity_star
  public :: get_xgravity
  public :: potential,acceleration

!ajwm SHOULDN'T BE SHARED
  public :: gravz,nu_epicycle,g0,gravz_const
  public :: gravz_profile
  public :: zref,z1,z2,zinfty,zgrav,reduced_top
  public :: lnrho_bot,lnrho_top
  public :: ss_bot,ss_top
  public :: lnumerical_equilibrium

