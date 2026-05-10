!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: initialize_timeavg
  public :: update_timeavgs
  public :: wsnap_timeavgs
  public :: pushpars2c

! Variables
! [ajwm] SHOULDN'T BE SHARED
! [PABourdin] Should be in cdata.f90, because these are namelist parameters.
! But "idx_tavg" is of dimension (mtavg) that is not yet known in cdata.
  public :: idx_tavg

