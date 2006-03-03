!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_planet,initialize_planet
  public :: rprint_planet
  public :: read_planet_init_pars,  write_planet_init_pars
  public :: read_planet_run_pars,  write_planet_run_pars
  public :: pencil_criteria_planet
  
  public :: gravity_companion,local_isothermal
  public :: gravity_star,wave_damping !,calc_torque,gravity_companion

!public stuff to use "abroad"
  public :: llocal_iso,gc,b,lwavedamp,lramp
  
