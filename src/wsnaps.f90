! $Id: wsnaps.f90,v 1.11 2002-10-09 14:05:31 mee Exp $

!!!!!!!!!!!!!!!!!!!!!!!
!!!   wsnaps.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Distributed IO (i.e. each process writes its own file data/procX)
!!!  07-Nov-2001/wd: Put into separate module, so one can choose
!!!  alternative IO mechanism.

module Wsnaps

  implicit none

contains

!***********************************************************************
    subroutine wsnap(chsnap,a,llabel)
!
!  Write snapshot file, labelled consecutively if llabel==.true.
!  Otherwise just write a snapshot without label (used for var.dat)
!
!  30-sep-97/axel: coded
!  08-oct-02/tony: expanded file to handle 120 character datadir // '/tsnap.dat'
!
      use Cdata
      use Mpicomm
      use Boundcond
      use Sub
      use Io
!
      real, dimension (mx,my,mz,mvar) :: a
      character (len=4) :: ch
      character (len=135) :: file
      character (len=*) :: chsnap
      logical lsnap,llabel
      integer, save :: ifirst,nsnap
      real, save :: tsnap
!
!  Output snapshot with label in 'tsnap' time intervals
!  file keeps the information about number and time of last snapshot
!
      if (llabel) then
        file=trim(datadir)//'/tsnap.dat'
!
!  at first call, need to initialize tsnap
!  tsnap calculated in out1, but only available to root processor
!
        if (ifirst==0) then
          call out1 (trim(file),tsnap,nsnap,dsnap,t)
          ifirst=1
        endif
!
!  Check whether we want to output snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently)
!
        call out2 (trim(file),tsnap,nsnap,dsnap,t,lsnap,ch,.true.)
        if (lsnap) then
          call update_ghosts(a)
          call output(chsnap//ch,a,mvar)
        endif
!
      else
!
!  write snapshot without label (typically, var.dat)
!
        call update_ghosts(a)
        call output(chsnap,a,mvar)
      endif
!
    endsubroutine wsnap
!***********************************************************************

endmodule Wsnaps

