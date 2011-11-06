! $Id$
!
!  Module to handle variables whose state should persist between executions of
!  run.x, e.g. the random number seeds and some other forcing state information.
!
!  25-Apr-2005/tony: Implemented initial try at backwards compatible
!                    additions to var.dat files.
!
!  The idea is to use integer block and record type tags to store arbitrary
!  extra information in the var files along with the actual field information.
!
!  The integers representing the various block/record types are defined in a
!  separate file, record_types.inc.  These numbers MUST remain unique and MUST
!  not be altered, though adding new types is acceptable (else old var.dat
!  files may become unreadable).
!
module Persist
!
  use Cdata
!
  implicit none
!
  private
!
  public :: input_persistent, output_persistent
!
  include 'record_types.h'
!
  contains
!***********************************************************************
    subroutine input_persistent(lun)
!
!  Read auxiliary information from snapshot file.
!  lun should be set to the same lun as that of the snapshot.
!
!  26-may-03/axel: adapted from output_vect
!   6-apr-08/axel: added input_persistent_magnetic
!
      use Interstellar, only: input_persistent_interstellar
      use Forcing, only: input_persistent_forcing
      use General, only: input_persistent_general
      use Magnetic, only: input_persistent_magnetic
      use Hydro, only: input_persistent_hydro
!
      integer :: lun
      integer :: id, dummy,ierr
      logical :: done =.false.
!
      if ((ip<=8).and.lroot) print*,'input_persistent: '
!
      read(lun,end=1000) id
      if (id/=id_block_PERSISTENT) then
        if ((ip<=8).and.lroot) &
            print*,'input_persistent: No persistent data to read'
        return
      endif
!
dataloop: do
        read(lun,iostat=ierr,end=1000) id
        done=.false.
        if (id==id_block_PERSISTENT) then
          exit dataloop
        endif
        if (ierr<0) exit dataloop
        if (.not.done) call input_persistent_general(id,lun,done)
        if (.not.done) call input_persistent_interstellar(id,lun,done)
        if (.not.done) call input_persistent_forcing(id,lun,done)
        if (.not.done) call input_persistent_magnetic(id,lun,done)
        if (.not.done) call input_persistent_hydro(id,lun,done)
        if (.not.done) read(lun,end=1000) dummy
      enddo dataloop
!
      if ((ip<=8).and.lroot) print*,'input_persistent: DONE'
      return
1000  if ((ip<=8).and.lroot) print*,'input_persistent: EOF termination'
!
    endsubroutine input_persistent
!***********************************************************************
    subroutine output_persistent(lun_output)
!
!  Write auxiliary information into snapshot file.
!  lun should be set to the same lun as that of the snapshot
!
!  26-may-03/axel: adapted from output_vect
!   6-apr-08/axel: added output_persistent_magnetic
!   5-nov-11/MR: IOSTAT handling added
!
      use Interstellar, only: output_persistent_interstellar
      use Forcing, only: output_persistent_forcing
      use General, only: output_persistent_general
      use Magnetic, only: output_persistent_magnetic
      use Hydro, only: output_persistent_hydro
      use Messages, only: outlog
!
      integer :: lun_output
!
      integer :: iostat
!
      if ((ip<=8).and.lroot) print*,'output_persistent: '
!
      write(lun_output,iostat=IOSTAT) id_block_PERSISTENT
      call outlog(iostat,'write id_block_PERSISTENT')
      call output_persistent_general(lun_output)
      call output_persistent_interstellar(lun_output)
      call output_persistent_forcing(lun_output)
      call output_persistent_magnetic(lun_output)
      call output_persistent_hydro(lun_output)
      write(lun_output,iostat=IOSTAT) id_block_PERSISTENT
      call outlog(iostat,'write id_block_PERSISTENT')
!
    endsubroutine output_persistent
!***********************************************************************
endmodule Persist
