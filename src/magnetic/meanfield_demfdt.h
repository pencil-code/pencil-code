! $Id$
!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
  private

  public :: register_magnetic_mf_demfdt
  public :: initialize_magnetic_mf_demfdt, init_aa_mf_demfdt
  public :: pencil_criteria_magnetic_mf_demfdt
  public :: pencil_interdep_magnetic_mf_demfdt
  public :: calc_pencils_magnetic_mf_demfdt
  public :: demf_dt_meanfield
  public :: read_magnetic_mf_demfdt_init_pars, write_magnetic_mf_demfdt_init_pars
  public :: read_magnetic_mf_demfdt_run_pars, write_magnetic_mf_demfdt_run_pars
  public :: rprint_magnetic_mf_demfdt
