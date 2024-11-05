!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: initialize_timeavg
  public :: update_timeavgs
  public :: wsnap_timeavgs

! Variables
! [ajwm] SHOULDN'T BE SHARED
! [PABourdin] should be in cdata.f90, because these are namelist parameters
  public :: tavg, idx_tavg

