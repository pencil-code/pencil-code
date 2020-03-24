!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
!
! functions
  public :: register_hydro, initialize_hydro
  public :: read_hydro_init_pars, write_hydro_init_pars
  public :: read_hydro_run_pars,  write_hydro_run_pars
  public :: input_persistent_hydro, output_persistent_hydro
  public :: rprint_hydro
  public :: get_slices_hydro
  public :: init_uu, duu_dt, hydro_after_boundary, calc_pencils_hydro
  public :: time_integrals_hydro
  public :: pencil_criteria_hydro, pencil_interdep_hydro
  public :: calc_mflow, remove_mean_momenta, remove_mean_flow
  public :: impose_velocity_ceiling
  public :: hydro_clean_up
  public :: coriolis_cartesian
  public :: kinematic_random_phase, kinematic_random_ampl
  public :: kinematic_random_wavenumber
  public :: hydro_before_boundary
  public :: expand_shands_hydro
  public :: calc_means_hydro
  public :: update_char_vel_hydro
  public :: hydro_after_timestep
  public :: calc_gradu
  public :: pushpars2c, pushdiags2c
  public :: calc_diagnostics_hydro, df_diagnos_hydro
!
! WL: SHOULDN'T BE PUBLIC!
!
  public :: uumx,        lcalc_uumeanx
  public :: uumz, guumz, lcalc_uumeanz, lupw_uu
  public :: uumxy, uumxz,lcalc_uumeanxy, lcalc_uumeanxz
  public :: ampl_fcont_uu

  interface calc_pencils_hydro
    module procedure calc_pencils_hydro_pencpar
    module procedure calc_pencils_hydro_std
  endinterface calc_pencils_hydro
