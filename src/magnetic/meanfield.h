! $Id: magnetic.h 14416 2010-07-23 15:13:00Z Bourdin.KIS $
!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
  private

  public :: initialize_magnetic_mf
  public :: read_magnetic_mf_init_pars, write_magnetic_mf_init_pars
  public :: read_magnetic_mf_run_pars, write_magnetic_mf_run_pars
  public :: pencil_criteria_magnetic_mf
  public :: pencil_interdep_magnetic_mf
  public :: calc_pencils_magnetic_mf
  public :: rprint_magnetic_mf
  public :: daa_dt_meanfield
