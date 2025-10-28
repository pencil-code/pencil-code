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
!  separate file, record_types.h.  These numbers MUST remain unique and MUST
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
  interface input_persistent_general
     module procedure input_persist_general_by_id
     module procedure input_persist_general_by_label
  endinterface
!
  contains
!***********************************************************************
    subroutine input_persistent(file)
!
!  Read auxiliary information from snapshot file.
!  lun should be set to the same lun as that of the snapshot.
!
!  26-may-03/axel: adapted from output_vect
!   6-apr-08/axel: added input_persistent_magnetic
!
      use IO, only: init_read_persist, read_persist_id, IO_strategy
      use Hydro, only: input_persistent_hydro
      use Interstellar, only: input_persistent_interstellar
      use Forcing, only: input_persistent_forcing
      use Magnetic, only: input_persistent_magnetic
      use Special, only: input_persistent_special
!
      character (len=*), intent(in), optional :: file
!
      integer :: id
      logical :: done
!
      if (lroot .and. (ip <= 8)) print *, 'input_persistent: START '//trim (file)
!
      if (IO_strategy == 'HDF5') then
        if (.not. read_persist_id ('INITIAL_BLOCK_ID', id, .true.)) return
        call input_persistent_general
        call input_persistent_hydro
        call input_persistent_interstellar
        call input_persistent_forcing
        call input_persistent_magnetic
        call input_persistent_special
        return
      endif
!
      if (read_persist_id ('INITIAL_BLOCK_ID', id, .true.)) then
        if (.not. present (file)) return
        if (file == 'var.dat') then
          if (init_read_persist ('pers_'//file)) return
        elseif (index (file, 'VAR') == 1) then
          if (init_read_persist ('PERS_'//file(4:))) return
        else
          return
        endif
        if (read_persist_id ('INITIAL_BLOCK_ID', id)) return
      else
        if (init_read_persist ()) return
!print*,'nach read_persist_id, INITIAL3: id,id_block_PERSISTENT=', id,id_block_PERSISTENT
      endif
!
      if (id /= id_block_PERSISTENT) then
        if (lroot .and. (ip <= 8)) print *, 'input_persistent: Missing initial persistent block ID'
        return
      endif
!
      if (read_persist_id ('FIRST_BLOCK_ID', id)) return
      do while (id /= id_block_PERSISTENT)
        done = .false.
        if (.not. done) call input_persistent_general (id, done)
        if (.not. done) call input_persistent_hydro (id, done)
        if (.not. done) call input_persistent_interstellar (id, done)
        if (.not. done) call input_persistent_forcing (id, done)
        if (.not. done) call input_persistent_magnetic (id, done)
        if (.not. done) call input_persistent_special (id, done)
        if (read_persist_id ('NEXT_BLOCK_ID', id)) return
      enddo
!
      if (lroot .and. (ip <= 8)) print *, 'input_persistent: DONE'
!
    endsubroutine input_persistent
!***********************************************************************
    subroutine output_persistent(file)
!
!  Write auxiliary information into snapshot file.
!  lun should be set to the same lun as that of the snapshot
!
!  26-may-03/axel: adapted from output_vect
!   6-apr-08/axel: added output_persistent_magnetic
!  16-nov-11/MR: calls adapted
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use IO, only: init_write_persist
      use Hydro, only: output_persistent_hydro
      use Interstellar, only: output_persistent_interstellar
      use Forcing, only: output_persistent_forcing
      use Magnetic, only: output_persistent_magnetic
      use Special, only: output_persistent_special
!
      character (len=*), intent(in) :: file
!
      if (lroot .and. (ip <= 8)) print *, 'output_persistent: START '//trim (file)
!
      if (lseparate_persist) then
        if ((file == 'var.dat') .or. (file == 'crash.dat')) then
          if (init_write_persist('pers_'//file)) return
        elseif (index (file, 'VAR') == 1) then
          if (init_write_persist('PERS_'//file(4:))) return
        endif
      endif
!
      if (output_persistent_general()) return
      if (output_persistent_hydro()) return
      if (output_persistent_interstellar()) return
      if (output_persistent_forcing()) return
      if (output_persistent_magnetic()) return
      if (output_persistent_special()) return
!
    endsubroutine output_persistent
!***********************************************************************
    subroutine input_persist_general_by_id(id,done)
!
!  Reads seed from a snapshot.
!  A read is performed depending on id
!
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use Cdata, only: seed, nseed
      use General, only: random_seed_wrapper
      use IO, only: read_persist
!
      integer, intent(in) :: id
      logical, intent(inout) :: done
!
      real :: dely, dtmp, deps
!
      select case (id)
        case (id_record_RANDOM_SEEDS)
          call random_seed_wrapper (GET=seed,CHANNEL=1)
          done=read_persist ('RANDOM_SEEDS', seed(1:nseed))
          if (.not.done) call random_seed_wrapper (PUT=seed,CHANNEL=1)
        case (id_record_RANDOM_SEEDS2)
          call random_seed_wrapper (GET=seed2,CHANNEL=2)
          done=read_persist ('RANDOM_SEEDS2', seed2(1:nseed))
          if (.not.done) call random_seed_wrapper (PUT=seed2,CHANNEL=2)
        case (id_record_SHEAR_DELTA_Y)
          done=read_persist ('SHEAR_DELTA_Y', dely)
          if (.not.done) deltay=dely
        case (id_record_TIME_STEP)
          if (dt0==0) then
            done=read_persist ('TIME_STEP', dtmp)
            if (.not.done) dt0=dtmp
          endif
        case (id_record_EPS_RKF)
          if (eps_rkf0==0) then
            done=read_persist ('EPS_RKF', deps)
            if (.not.done) eps_rkf=deps
          endif
      endselect
!
    endsubroutine input_persist_general_by_id
!***********************************************************************
    subroutine input_persist_general_by_label
!
!  Reads seed from a snapshot.
!
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use Cdata, only: seed, nseed
      use General, only: random_seed_wrapper
      use IO, only: persist_exists, read_persist
      use Messages, only: warning
!
      real :: dtmp, deps
!
      if (persist_exists ('RANDOM_SEEDS')) then
        call random_seed_wrapper(GET=seed,CHANNEL=1)
        if (.not. read_persist('RANDOM_SEEDS', seed(1:nseed))) then
           call random_seed_wrapper(PUT=seed,CHANNEL=1)
        else
          call warning('input_persist_general_by_label','no persistent value of seed1 found')
        endif
      endif
!
      if (persist_exists ('RANDOM_SEEDS2')) then
        call random_seed_wrapper(GET=seed2,CHANNEL=2)
        if (.not. read_persist ('RANDOM_SEEDS2', seed2(1:nseed))) then
           call random_seed_wrapper(PUT=seed2,CHANNEL=2)
        else
          call warning('input_persist_general_by_label','no persistent value of seed2 found')
        endif
      endif
!
      if (lshear) then
        if (read_persist ('SHEAR_DELTA_Y', deltay)) &
          call warning('input_persist_general_by_label','no persistent value of deltay found')
      endif
!
      if (dt0==0) then
        if (read_persist ('TIME_STEP', dtmp)) then
          call warning('input_persist_general_by_label','no persistent value of dt found')
        else
          dt0=dtmp
        endif
      endif
!
      if (eps_rkf0==0.) then
        deps=eps_rkf
        if (read_persist ('EPS_RKF', deps)) then
          call warning('input_persist_general_by_label','no persistent value of eps_rkf found')
        else
          eps_rkf=deps
        endif
      endif
!
    endsubroutine input_persist_general_by_label
!***********************************************************************
    logical function output_persistent_general()
!
!  Writes general persistent variables to a snapshot.
!
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use Cdata, only: seed, nseed, lshear, deltay
      use General, only: random_seed_wrapper
      use IO, only: write_persist

      real :: dt_save
!
      output_persistent_general = .false.
!
      ! Don't write the seeds, if they are unchanged from their default value.
      call random_seed_wrapper (GET=seed,CHANNEL=1)
      if (any (seed(1:nseed) /= seed0)) then
        if (write_persist ('RANDOM_SEEDS', id_record_RANDOM_SEEDS, seed(1:nseed))) &
            output_persistent_general = .true.
      endif
!
      call random_seed_wrapper (GET=seed2,CHANNEL=2)
      if (any (seed2(1:nseed) /= seed0)) then
        if (write_persist ('RANDOM_SEEDS2', id_record_RANDOM_SEEDS2, seed2(1:nseed))) &
            output_persistent_general = .true.
      endif
!
      if (lshear) then
        if (write_persist ('SHEAR_DELTA_Y', id_record_SHEAR_DELTA_Y, deltay)) &
            output_persistent_general = .true.
      endif
!
      if (.not.lstart.and.(.not.lcourant_dt)) then
        if(ldt .and. .not. lcourant_dt) then
          dt_save = dt_next
        else
          dt_save = dt
        endif
        if (write_persist ('TIME_STEP', id_record_TIME_STEP, dt_save)) then
          output_persistent_general = .true.
        endif
      endif
!
      if (.not.ldt) then
        if (write_persist ('EPS_RKF', id_record_EPS_RKF, eps_rkf)) &
          output_persistent_general = .true.
      endif
!
    endfunction output_persistent_general
!***********************************************************************
endmodule Persist
