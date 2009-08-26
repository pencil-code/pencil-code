!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim) 

  public :: initialize_border_profiles
  public :: request_border_driving
  public :: pencil_criteria_borderprofiles
  public :: calc_pencils_borderprofiles
  public :: border_quenching, border_driving
  public :: set_border_initcond

  integer, parameter, public :: i_BORDER_ZERO    = 1
  integer, parameter, public :: i_BORDER_SPECIAL = 2
  integer, parameter, public :: i_BORDER_        = 3
