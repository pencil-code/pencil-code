!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_hydro, initialize_hydro
  public :: read_hydro_init_pars, write_hydro_init_pars
  public :: read_hydro_run_pars,  write_hydro_run_pars
  public :: rprint_hydro
  public :: init_uu, duu_dt, calc_pencils_hydro
  public :: pencil_criteria_hydro, pencil_interdep_hydro

  public :: calc_mflow, calc_turbulence_pars

!ajwm SHOULDN'T BE EXPORTED
  public :: idiag_Marms, idiag_Mamax, idiag_dtv
  public :: orms, idiag_orms
!!  public :: nu_turb
  public :: nu_turb0, nu_turb1, tau_nuturb, lcalc_turbulence_pars
  public :: theta
  public :: ul0,tl0,teta,ueta,tl01,teta1
  public :: kep_cutoff_pos_ext,kep_cutoff_width_ext
  public :: kep_cutoff_pos_int,kep_cutoff_width_int,u_out_kep
