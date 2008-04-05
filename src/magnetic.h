!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_magnetic, initialize_magnetic
  public :: read_magnetic_init_pars, write_magnetic_init_pars
  public :: read_magnetic_run_pars,  write_magnetic_run_pars
  public :: output_persistent_magnetic, input_persistent_magnetic
  public :: rprint_magnetic
  public :: get_slices_magnetic
  public :: init_aa, daa_dt
  public :: df_diagnos_magnetic

  public :: pencil_criteria_magnetic, pencil_interdep_magnetic
  public :: calc_pencils_magnetic

  public :: calc_mfield
  public :: bc_aa_pot, bc_aa_pot2, bc_aa_pot3, bc_frozen_in_bb
  public :: pert_aa, rescaling
  public :: bb_unitvec_shock, remove_mean_emf
  public :: lelectron_inertia, inertial_length, linertial_2

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
