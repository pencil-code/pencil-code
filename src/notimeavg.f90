! $Id: notimeavg.f90,v 1.1 2002-10-07 20:39:22 dobler Exp $ 

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
    subroutine timeavg_run_hook(a)
!
      real, dimension(mx,my,mz,mvar) :: a
!
      if (ip < 0) print*, a(1,1,1,1)
!
    endsubroutine timeavg_run_hook
!***********************************************************************
    subroutine update_timeavgs(a,dt,init)
!
!
      real, dimension(mx,my,mz,mvar) :: a
      real :: dt
      logical, optional :: init
!
      if (ip < 0) print*, a(1,1,1,1),dt,present(init)
!
    endsubroutine update_timeavgs
!***********************************************************************
    subroutine write_timeavgs()
!
    endsubroutine write_timeavgs
!***********************************************************************

endmodule Timeavg

!!! End of file timeavg.f90
