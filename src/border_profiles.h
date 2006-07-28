!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim) 

  public :: initialize_border_profiles

  public :: border_quenching, border_driving

  integer, parameter, public :: i_BORDER_ZERO    = 1
  integer, parameter, public :: i_BORDER_SPECIAL = 2
  integer, parameter, public :: i_BORDER_        = 3
