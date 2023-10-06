!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_pointmasses, initialize_pointmasses
  public :: init_pointmasses
  public :: pencil_criteria_pointmasses, pencil_interdep_pointmasses
  public :: calc_pencils_pointmasses, calc_diagnostics_pointmasses
  public :: pointmasses_pde_pencil,pointmasses_pde
  public :: rprint_pointmasses
  public :: read_pointmasses_init_pars, write_pointmasses_init_pars
  public :: read_pointmasses_run_pars, write_pointmasses_run_pars
  public :: get_totalmass
  public :: pointmasses_read_snapshot
  public :: pointmasses_write_snapshot,pointmasses_write_qdim
  public :: pointmasses_timestep_first,pointmasses_timestep_second
  public :: boundconds_pointmasses

