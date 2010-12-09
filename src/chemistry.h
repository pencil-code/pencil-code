!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_chemistry, initialize_chemistry
  public :: read_chemistry_init_pars, write_chemistry_init_pars
  public :: read_chemistry_run_pars,  write_chemistry_run_pars
  public :: rprint_chemistry, get_slices_chemistry
  public :: init_chemistry

  public :: dchemistry_dt

  public :: calc_pencils_chemistry
  public :: pencil_criteria_chemistry
  public :: pencil_interdep_chemistry

  public :: calc_for_chem_mixture
  public :: chemspec_normalization
!  public :: bc_nscbc_subin_x
!  public :: bc_nscbc_nref_subout_x
!  public :: bc_nscbc_nref_subout_y
!  public :: bc_nscbc_nref_subout_z

  public :: jacobn
  public :: get_mu1_slice
  public :: get_gamma_slice
  public :: get_cs2_slice
  public :: get_cs2_full
  public :: get_gamma_full
  public :: get_RHS_Y_full
  public :: get_reac_rate
 

  public :: Rgas

!  public :: get_p_infx
!  public :: get_p_infy
!  public :: get_rhs_Y
!  public :: get_rhs_T
 

  public :: chemistry_clean_up

! public :: chemistry_calc_density
! public :: chemistry_calc_hydro
! public :: chemistry_calc_entropy
! public :: chemistry_calc_magnetic

! public :: chemistry_before_boundary
  public :: write_net_reaction
  public :: lchemistry_diag
  public :: lreactions
