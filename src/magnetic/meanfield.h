! $Id$
!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
  private

  public :: register_magnetic_mf
  public :: initialize_magnetic_mf, init_aa_mf
  public :: pencil_criteria_magnetic_mf
  public :: pencil_interdep_magnetic_mf
  public :: calc_pencils_magnetic_mf
  public :: daa_dt_meanfield
  public :: read_magnetic_mf_init_pars, write_magnetic_mf_init_pars
  public :: read_magnetic_mf_run_pars, write_magnetic_mf_run_pars
  public :: rprint_magnetic_mf
