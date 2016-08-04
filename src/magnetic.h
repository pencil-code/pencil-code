! $Id$
!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
  private

  public :: register_magnetic, initialize_magnetic
  public :: read_magnetic_init_pars, write_magnetic_init_pars
  public :: read_magnetic_run_pars,  write_magnetic_run_pars
  public :: output_persistent_magnetic, input_persistent_magnetic
  public :: rprint_magnetic
  public :: get_slices_magnetic
  public :: init_aa, daa_dt
  public :: magnetic_before_boundary, calc_lmagnetic_pars
  public :: time_integrals_magnetic
  public :: df_diagnos_magnetic

  public :: pencil_criteria_magnetic, pencil_interdep_magnetic
  public :: calc_pencils_magnetic
  public :: get_bext

  public :: calc_mfield, idiag_bcosphz, idiag_bsinphz
  public :: lcalc_aameanz, aamz, bbmz, jjmz
  public :: rescaling_magnetic
  public :: B_ext_inv
  public :: bb_unitvec_shock
  public :: lelectron_inertia, inertial_length, linertial_2
  public :: idiag_axmz,idiag_aymz
  public :: idiag_bxmz,idiag_bymz
  public :: dynamical_resistivity
  public :: split_update_magnetic
  public :: expand_shands_magnetic
  public :: update_char_vel_magnetic

  public :: lresi_dep
