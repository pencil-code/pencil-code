! $Id: timeavg.f90,v 1.1 2002-10-07 20:39:22 dobler Exp $ 

module Timeavg

!
!   Data array and infrastructure for time averages of fundamental
!   variables.
!
  use Cdata

  implicit none
!
  include 'ctimeavg.local'
!
  real, dimension(mx,my,mz,mtavg) :: f_tavg
  integer, dimension(mtavg) :: idx_tavg=0
  real :: tavg=0.

  logical :: ltavg=.true.

  contains

!***********************************************************************
    subroutine timeavg_run_hook(a)
!
!  Initialize time averages to the corresponding values
!
!  7-oct-02/wolf: coded
!
      real, dimension(mx,my,mz,mvar) :: a
      integer :: i
!
!  initialize values
!
      call update_timeavgs(a,0.,INIT=.true.)
!
!  Build list of indexes for f; defaults to indices 1 to mtavg, so
!  setting mtavg=mvar will automatically average all variables.
!  Currently we do not check for (inefficient, but not terrible) multiple
!  occurences of an index, which can occur when combining the defaulting
!  mechanism with setting of some indices.
!
      do i=1,mtavg
        if (idx_tavg(i) == 0) then
          print*, 'TIMEAVG: defaulting idx no.', i, ' to ', i
          idx_tavg(i)=i
        endif
      enddo
!
    endsubroutine timeavg_run_hook
!***********************************************************************
    subroutine update_timeavgs(a,dt,init)
!
!  Called after a time step, this routine updates the time averages from
!  the variable vector a.
!  With optional argument INIT=.true., initialize to a.
!
!  7-oct-02/wolf: coded
!
      real, dimension(mx,my,mz,mvar) :: a
      real :: dt,weight
      integer :: i,idx
      logical, optional :: init
      logical :: init1=.false.
!
      init1=.false.             ! somehow the initialization above doen
                                ! not seem to work on Cincinnatus
      if (tavg <= 0) return
print*, 'PRESENT 1: ', present(init),init1
      if (present(init)) init1=init
print*, 'PRESENT 1: ', present(init),init1
!
      weight = min(dt/tavg,1.)
      do i=1,mtavg
        idx = idx_tavg(i)
        if (idx > 0) then
          if (init1) then
print*, 'Initializing...'
            f_tavg(:,:,:,i) = a(:,:,:,idx)
print*, '  minmax = ', minval(f_tavg), maxval(f_tavg)
          else
print*, 'Averaging...'
            f_tavg(:,:,:,i) = weight*a(:,:,:,idx) + (1-weight)*f_tavg(:,:,:,i)
print*, '  minmax = ', minval(f_tavg), maxval(f_tavg)
          endif
        endif
      enddo
!
    endsubroutine update_timeavgs
!***********************************************************************
    subroutine write_timeavgs()
!
      use IO
!
      call output(trim(directory)//'/timeavg.dat',f_tavg,mtavg)
!
    endsubroutine write_timeavgs
!***********************************************************************

endmodule Timeavg

!!! End of file timeavg.f90
