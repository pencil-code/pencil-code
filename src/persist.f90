! $Id: persist.f90,v 1.4 2005-12-21 16:45:13 mee Exp $

!!!!!!!!!!!!!!!!!!!!!!!
!!!   persist.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Module to handle variables whose state should persist
!!!  between executions of run.x, eg. the randome number seeds
!!!  and some other forcing state information. 
!!!
!!!  25-Apr-2005/tony: Implemented initial try at backwards compatible
!!!                    additions to var.dat files.
!!!
!!!  The idea is to use integer block and record type tags to store 
!!!  arbitrary extra information in the var files along with the
!!!  actual field information.
!!!
!!!  The integers representing the various block/record types are defined
!!!  in a separate file, record_types.inc.  These numbers MUST remain unique
!!!  and MUST not be altered, though adding new types is acceptable.
!!!  (Else old var.dat files may become unreadable)
!!!

module Persist

  implicit none

  private

  public :: input_persistent, output_persistent

  include 'record_types.h'

contains
!***********************************************************************
    subroutine input_persistent(lun)
!
!  write auxiliary information into snapshot file
!
!  26-may-03/axel: adapted from output_vect
!
      use Cdata
      Use Interstellar, only: input_persistent_interstellar
      Use Forcing, only: input_persistent_forcing
      use General, only: input_persistent_general
!
      integer :: lun
      integer :: id, dummy,ierr
      logical :: done =.false.
!
      if ((ip<=8).and.lroot) print*,'input_persistent: '
!
      read(lun,end=1000) id
      if (id/=id_block_PERSISTANT) then
        if (lroot) print*,'input_persistent: No persistent data to read'
        return
      endif
      
dataloop: do 
        read(lun,iostat=ierr,end=1000) id
        done=.false.
        if (id==id_block_PERSISTANT) then
          exit dataloop
        endif
        if (ierr<0) exit dataloop
        if (.not.done) call input_persistent_general(id,lun,done)
        if (.not.done) call input_persistent_interstellar(id,lun,done)
        if (.not.done) call input_persistent_forcing(id,lun,done)
        if (.not.done) read(lun,end=1000) dummy
      enddo dataloop

      if (lroot) print*,'input_persistent: DONE'
      return
1000  if (lroot) print*,'input_persistent: EOF termination'
!
    endsubroutine input_persistent
!***********************************************************************
    subroutine output_persistent(lun_output)
!
!  write auxiliary information into snapshot file
!
!  26-may-03/axel: adapted from output_vect
!
      use Cdata
      use Interstellar, only: output_persistent_interstellar
      use Forcing, only: output_persistent_forcing
      use General, only: output_persistent_general
!
      integer :: lun_output
!
      if ((ip<=8).and.lroot) print*,'output_persistent: '
!
      write(lun_output) id_block_PERSISTANT
      call output_persistent_general(lun_output)
      call output_persistent_interstellar(lun_output)
      call output_persistent_forcing(lun_output)
      write(lun_output) id_block_PERSISTANT
!
    endsubroutine output_persistent

endmodule Persist
