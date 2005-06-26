!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: initialize_timeavg
  public :: update_timeavgs
  public :: wsnap_timeavgs

! Variables
  public :: ltavg
  
!ajwm SHOULDN'T BE SHARED
  public :: tavg, idx_tavg

