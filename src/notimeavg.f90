! $Id: notimeavg.f90,v 1.5 2003-06-16 04:41:11 brandenb Exp $ 

module Timeavg

!
!   Dummy module
!
  use Cdata

  implicit none
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
    subroutine wsnap_timeavgs(chsnap,llabel)
!
      character (len=*) :: chsnap
      logical llabel
!
      if(chsnap=='X') llabel=.false. !(to keep compiler quiet)
    endsubroutine wsnap_timeavgs
!***********************************************************************

endmodule Timeavg

!!! End of file timeavg.f90
