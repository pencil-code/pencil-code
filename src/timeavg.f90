! $Id: timeavg.f90,v 1.9 2003-08-07 17:06:56 dobler Exp $ 

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
    subroutine initialize_timeavg(a)
!
!  Initialize time averages to the corresponding values
!
!  7-oct-02/wolf: coded
!
      real, dimension(mx,my,mz,mvar+maux) :: a
      integer :: i

      intent (in) :: a
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
          if (lroot) print*, 'TIMEAVG: defaulting idx no.', i, ' to ', i
          idx_tavg(i)=i
        endif
      enddo
!
    endsubroutine initialize_timeavg
!***********************************************************************
    subroutine update_timeavgs(a,dt,init)
!
!  Called after a time step, this routine updates the time averages from
!  the variable vector a.
!  With optional argument INIT=.true., initialize to a.
!
!  7-oct-02/wolf: coded
!
      real, dimension(mx,my,mz,mvar+maux) :: a
      real :: dt,weight
      integer :: i,idx
      logical, optional :: init
      logical :: init1=.false.

      intent (in) :: a
!
      init1=.false.             ! somehow the initialization above doen
                                ! not seem to work on Cincinnatus
      if (tavg <= 0) return
      if (present(init)) init1=init
!
      weight = min(dt/tavg,1.)
      do i=1,mtavg
        idx = idx_tavg(i)
        if (idx > 0) then       ! should always be the case
          if (init1) then
            f_tavg(:,:,:,i) = a(:,:,:,idx)
          else
            f_tavg(:,:,:,i) = weight*a(:,:,:,idx) + (1-weight)*f_tavg(:,:,:,i)
          endif
        endif
      enddo
!
    endsubroutine update_timeavgs
!***********************************************************************
    subroutine wsnap_timeavgs(chsnap,enumerate)
!
!  Write snapshot file for time averages, labelled consecutively if
!  enumerate==.true. 
!  Otherwise write without label (used for timeavg.dat).
!   Currently uses the same parameter dsnap as wsnaps
!
!   9-oct-02/wolf: adapted from wsnaps
!
      use Cdata
      use General
      use Mpicomm
      use Sub
      use Io
!
      character (len=4) :: ch
      character (len=80) :: file
      character (len=*) :: chsnap
      logical lsnap,enumerate
      integer, save :: ifirst,nsnap
      real, save :: tsnap
!
!  Output snapshot with label in 'tsnap' time intervals
!  file keeps the information about number and time of last snapshot
!
      if (enumerate) then
        ! file=trim(datadir)//'/tavgsnap.dat'
        call safe_character_assign(file, trim(datadir)//'/tavgsnap.dat')
!
!  at first call, need to initialize tsnap
!  tsnap calculated in read_snaptime, but only available to root processor
!
        if (ifirst==0) then
          call read_snaptime(trim(file),tsnap,nsnap,dsnap,t)
          ifirst=1
        endif
!
!  Check whether we want to output snapshot.
!
        call update_snaptime(trim(file),tsnap,nsnap,dsnap,t,lsnap,ch, &
                             ENUMERATE=.true.)
        if (lsnap) then
          call output(chsnap//ch,f_tavg,mtavg)
        endif
!
      else
!
!  write snapshot without label (typically, timeavg.dat)
!
        call output(chsnap,f_tavg,mtavg)
      endif
!
    endsubroutine wsnap_timeavgs
!***********************************************************************

endmodule Timeavg

!!! End of file timeavg.f90
