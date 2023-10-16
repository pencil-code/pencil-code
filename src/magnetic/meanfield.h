! $Id$
!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
  private

  public :: register_magn_mf
  public :: initialize_magn_mf, init_aa_mf
  public :: pencil_criteria_magn_mf
  public :: pencil_interdep_magn_mf
  public :: calc_pencils_magn_mf, calc_diagnostics_meanfield
  public :: daa_dt_meanfield
  public :: read_magn_mf_init_pars, write_magn_mf_init_pars
  public :: read_magn_mf_run_pars, write_magn_mf_run_pars
  public :: rprint_magn_mf
!- public :: meanfield_chitB
  public :: pc_aasb_const_alpha
  public :: meanfield_after_boundary
