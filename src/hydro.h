!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_hydro, initialize_hydro
  public :: read_hydro_init_pars, write_hydro_init_pars
  public :: read_hydro_run_pars,  write_hydro_run_pars
  public :: output_persistent_hydro, input_persistent_hydro
  public :: rprint_hydro
  public :: get_slices_hydro
  public :: init_uu, duu_dt, calc_lhydro_pars, calc_pencils_hydro
  public :: time_integrals_hydro
  public :: pencil_criteria_hydro, pencil_interdep_hydro
  public :: calc_mflow, remove_mean_momenta, impose_velocity_ceiling
  public :: uumz,guumz,lcalc_uumean,lupw_uu
  public :: lforcing_cont_uu, ampl_fcont_uu
  public :: hydro_clean_up
  public :: traceless_strain, coriolis_cartesian
  public:: kinematic_random_phase
!ajwm SHOULDN'T BE EXPORTED
!
! Keplerian velocity boundary condition parameters
!    (needed by boundcond.f90)
!
!WL: commented out the routine that needed them
!
!  public :: kep_cutoff_pos_ext,kep_cutoff_width_ext
!  public :: kep_cutoff_pos_int,kep_cutoff_width_int
!  public :: u_out_kep
