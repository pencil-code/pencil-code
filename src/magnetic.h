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
  public :: init_aa, daa_dt, calc_lmagnetic_pars
  public :: time_integrals_magnetic
  public :: df_diagnos_magnetic

  public :: pencil_criteria_magnetic, pencil_interdep_magnetic
  public :: calc_pencils_magnetic

  public :: calc_mfield, idiag_bcosphz, idiag_bsinphz
  public :: lcalc_aamean, aamz, bbmz, jjmz
! public :: pert_aa, rescaling_magnetic
  public :: rescaling_magnetic
  public :: bb_unitvec_shock, remove_mean_emf
  public :: lelectron_inertia, inertial_length, linertial_2
  public :: idiag_axmz,idiag_aymz
  public :: idiag_bxmz,idiag_bymz

!ajwm SHOULDN'T BE SHARED
!
! Used to get parameters into nohydro for kinematic dynamo simulations!
!
  public :: ABC_A, KZ_AA, ABC_C, KY_AA, ABC_B, KX_AA
!
!ajwm  Are these totally dead now? [29-03-06]
!
  !public :: eta !(needed for alpm [20-11-04/axel])
  !public :: meanfield_EMFdotB !(needed for alpm [20-11-04/axel])
