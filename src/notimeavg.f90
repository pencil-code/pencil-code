! $Id: notimeavg.f90,v 1.10 2005-06-26 17:34:13 eos_merger_tony Exp $ 

module Timeavg

!
!   Dummy module
!
  use Cdata

  implicit none

  include 'timeavg.inc' 
!
!  real, dimension(mx,my,mz,mtavg) :: f_tavg
!  integer, dimension(mtavg) :: idx_tavg=0
!  real :: tavg=0.

  integer :: idx_tavg=0         ! just scalar, since unused and no mtavg known
  real :: tavg=0.

  logical :: ltavg=.false.

  contains

!***********************************************************************
    subroutine initialize_timeavg(a)
!
      real, dimension(mx,my,mz,mvar+maux) :: a
      
      intent (in) :: a
!
      if (ip < 0) print*, a(1,1,1,1)
!
    endsubroutine initialize_timeavg
!***********************************************************************
    subroutine update_timeavgs(a,dt,init)
!
!
      real, dimension(mx,my,mz,mvar+maux) :: a
      real :: dt
      logical, optional :: init

      intent (in) :: a
!
      if (ip < 0) print*, a(1,1,1,1),dt,present(init)
!
    endsubroutine update_timeavgs
!***********************************************************************
    subroutine wsnap_timeavgs(chsnap,enum,flist)
!
      character (len=*) :: chsnap,flist
      logical :: enum
      optional :: flist
!
      if(chsnap=='X') enum=.false. !(to keep compiler quiet)
      if(NO_WARN) print*,flist !(to keep compiler quiet)
    endsubroutine wsnap_timeavgs
!***********************************************************************

endmodule Timeavg

!!! End of file timeavg.f90
