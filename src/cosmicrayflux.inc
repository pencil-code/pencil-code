!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_cosmicrayflux, initialize_cosmicrayflux
  public :: read_cosmicrayflux_init_pars, write_cosmicrayflux_init_pars
  public :: read_cosmicrayflux_run_pars,  write_cosmicrayflux_run_pars
  public :: rprint_cosmicrayflux
  public :: init_fcr, dfcr_dt

  public :: pencil_criteria_cosmicrayflux, pencil_interdep_cosmicrayflux
  public :: calc_pencils_cosmicrayflux

!  public :: calc_mfield
!  public :: bc_aa_pot, bc_frozen_in_bb_z
!  public :: pert_aa, rescaling

!ajwm SHOULDN'T BE SHARED
!  public :: brms, idiag_brms, idiag_bmx, idiag_bmy
!  public :: kinflow
!  public :: ABC_A, KZ_AA, ABC_C, KY_AA, ABC_B, KX_AA
  !public :: eta !(needed for alpm [20-11-04/axel])
  !public :: meanfield_EMFdotB !(needed for alpm [20-11-04/axel])
