! $Id: wsnaps.f90,v 1.13 2002-11-09 08:31:56 brandenb Exp $

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
          if(ip<=10.and.lroot) print*,'wsnap: written snapshot ',chsnap//ch
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
   subroutine powersnap(a)
!
!  Write a snapshot of power spectrum 
!
!  30-sep-97/axel: coded
!  07-oct-02/nils: adapted from wsnap
!  08-oct-02/tony: expanded file to handle 120 character datadir // '/tspec.dat'
!
      use Cdata
      use Mpicomm
      use Boundcond
      use Sub
      use Io
      use Power_spectrum
!
      real, dimension (mx,my,mz,mvar) :: a
      character (len=135) :: file
      character (len=4) :: ch
      logical lspec
      integer, save :: ifirst,nspec
      real, save :: tspec
!
!  Output snapshot in 'tpower' time intervals
!  file keeps the information about time of last snapshot
!
      file=trim(datadir)//'/tspec.dat'
!
!  at first call, need to initialize tspec
!  tspec calculated in out1, but only available to root processor
!
      if (ifirst==0) then
         call out1 (trim(file),tspec,nspec,dspec,t)
         ifirst=1
      endif
!
!  Check whether we want to output power snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently)
!
      call out2 (trim(file),tspec,nspec,dspec,t,lspec,ch,.false.)
      if (lspec) then
         call update_ghosts(a)
         if (vel_spec) call power(a,'u')
         if (mag_spec) call power(a,'b')
         if (vec_spec) call power(a,'a')
         if (ab_spec)  call powerhel(a,'mag')
         if (ou_spec)  call powerhel(a,'kin')
      endif
!
    endsubroutine powersnap
!***********************************************************************
endmodule Wsnaps

