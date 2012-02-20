! $Id$

module Timeavg

!
!   Data array and infrastructure for time averages of fundamental
!   variables.
!
  use Cdata

  implicit none

  include 'timeavg.h'
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
      real, dimension(mx,my,mz,mfarray) :: a
      integer :: i

      intent (in) :: a
!
!  initialize values
!
      call update_timeavgs(a,0.,INIT=.true.)
!
!  Build list of indices for f; defaults to indices 1 to mtavg, so
!  setting mtavg=mvar will automatically average all variables.
!  Currently we do not check (inefficient, but not terrible) for multiple
!  occurences of an index, which can occur when combining the default
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
      real, dimension(mx,my,mz,mfarray) :: a
      real :: dt,weight
      integer :: i,idx
      logical, optional :: init
      logical :: init1=.false.

      intent (in) :: a
!
      init1=.false.             ! somehow the initialization above does
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
!            f_tavg(:,:,:,i) = weight*a(:,:,:,idx) + (1-weight)*f_tavg(:,:,:,i)
            ! numerically slightly better (probably irrelevant):
            f_tavg(:,:,:,i) = f_tavg(:,:,:,i) &
                              + weight*(a(:,:,:,idx)-f_tavg(:,:,:,i) )
          endif
        endif
      enddo
!
    endsubroutine update_timeavgs
!***********************************************************************
    subroutine wsnap_timeavgs(chsnap,enum,flist)
!
!  Write snapshot file for time averages, labelled consecutively if
!  enum==.true.
!  Otherwise write without label (used for timeavg.dat).
!  If `flist' is given, add name of snapshot to file flist.
!
!   9-oct-02/wolf: adapted from wsnaps
!
      use Cdata
      use General
      use IO, only: log_filename_to_file, output_globals
      use Mpicomm
      use Sub
!
      character (len=intlen) :: ch
      character (len=fnlen) :: file
      character (len=*) :: chsnap,flist
      logical :: lsnap,enum
      integer, save :: ifirst=0,nsnap
      real, save :: tsnap
      optional :: flist
!
!  Output snapshot with label in 'tsnap' time intervals
!  file keeps the information about number and time of last snapshot
!
      if (enum) then
        ! file=trim(datadir)//'/tavgsnap.dat'
        call safe_character_assign(file, trim(datadir)//'/tavgsnap.dat')
!
!  at first call, need to initialize tsnap
!  tsnap calculated in read_snaptime, but only available to root processor
!
        if (ifirst==0) then
          call read_snaptime(file,tsnap,nsnap,tavg,t)
          ifirst=1
        endif
!
!  Check whether we want to output snapshot.
!
        call update_snaptime(file,tsnap,nsnap,tavg,t,lsnap,ch)
        if (lsnap) then
          call output_globals(chsnap//ch,f_tavg,mtavg)
          if (present(flist)) call log_filename_to_file(chsnap//ch,flist)
        endif
!
      else
!
!  write snapshot without label (typically, timeavg.dat)
!
        call output_globals(chsnap,f_tavg,mtavg)
      endif
!
    endsubroutine wsnap_timeavgs
!***********************************************************************

endmodule Timeavg

!!! End of file timeavg.f90
