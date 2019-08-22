! $Id$
!
!  This module goes straight and implements distributed file Input/Output.
!  We use here only F95 features and are explicitly not HPC-friendly.
!
module File_io
!
  implicit none
!
  integer, parameter :: parallel_unit=14, parallel_unit_vec=14
!
  include 'file_io.h'

  private

  contains
!***********************************************************************
    subroutine parallel_open(file,form,remove_comments,nitems)
!
!  Read a global file.
!
!  18-mar-10/Bourdin.KIS: implemented
!
      use Messages, only: fatal_error
      use Cparam, only: ALWAYS_FALSE
!
      character(len=*)          , intent(in) :: file
      character(len=*), optional, intent(in) :: form
      logical,          optional, intent(in) :: remove_comments
      integer,          optional, intent(out):: nitems
!
      logical :: exists
!
!  Test if file exists.
!
      if (present(nitems)) nitems=0
      if (ALWAYS_FALSE) print*, present(remove_comments)
!
      inquire(file=file,exist=exists)
      if (.not. exists) call fatal_error('parallel_open', 'file "'//trim(file)//'" not found')
!
!  Open file.
!
      if (present(form)) then
        open(parallel_unit, file=file, FORM=form, status='old')
      else
        open(parallel_unit, file=file, status='old')
      endif
!
    endsubroutine parallel_open
!***********************************************************************
    subroutine parallel_rewind
!
!  Rewind a file unit opened by parallel_open.
!
!  23-May-2014/Bourdin.KIS: implemented
!
      rewind(parallel_unit)
!
    endsubroutine parallel_rewind
!***********************************************************************
    subroutine parallel_close
!
!  Close a file unit opened by parallel_open and remove temporary file.
!
!  17-mar-10/Bourdin.KIS: implemented
!
      close(parallel_unit)
!
    endsubroutine parallel_close
!**************************************************************************
    !function find_namelist(name) result (lfound)
    subroutine find_namelist(name,lfound)
!
!  Tests if the namelist is present and reports a missing namelist.
!
!  26-Sep-2015/PABourdin: coded
!   6-oct-2015/MR: turned into subroutine because of CRAY compiler bug;
!                  easily revertable by shifting comment char at beginning and end.
!
      use Cdata, only: comment_char
      use General, only: lower_case
      use Mpicomm, only: lroot, mpibcast, MPI_COMM_WORLD
      use Messages, only: warning
!
      character(len=*), intent(in) :: name
      logical :: lfound
!
      integer :: ierr, pos, state, max_len, line_len
      character(len=36000) :: line
      character :: ch
!
      if (lroot) then
        lfound = .false.
!
        max_len = len (name)
        ierr = 0
        do while (ierr == 0)
          state = 0
          read (parallel_unit,'(A)',iostat=ierr) line
          if (ierr /= 0) cycle
          line_len = len_trim (line)
          do pos = 1, line_len
            ch = lower_case (line(pos:pos))
            if (ch .eq. char(10)) then
              state = 0
              cycle
            endif
            if ((ch == '!') .or. (ch == comment_char)) state = -2
            if ((ch == ' ') .or. (ch == char(9)) .or. (state < 0)) cycle
            if ((state == 0) .and. (ch == '&')) then
              state = 1
              cycle
            endif
            if (state >= 1) then
              if (ch == lower_case (name(state:state))) then
                if (state == max_len) then
                  if (pos == line_len) then
                    lfound = .true.
                    exit
                  endif
                  ch = lower_case (line(pos+1:pos+1))
                  if ((ch == ' ') .or. (ch == char(9)) .or. (ch == '!') .or. (ch == comment_char)) then
                    lfound = .true.
                  endif
                  if (lfound) exit
                  state = -1
                  cycle
                endif
                state = state + 1
                cycle
              endif
            endif
            state = -1
          enddo
        enddo
        call parallel_rewind
        if (.not. lfound) call warning ('find_namelist', 'namelist "'//trim(name)//'" is missing!')
      endif
!
      call mpibcast (lfound,comm=MPI_COMM_WORLD)
!
    endsubroutine find_namelist
    !endfunction find_namelist
!***********************************************************************
    subroutine flush_file(unit)

      use General, only: keep_compiler_quiet

      integer, intent(IN) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine flush_file
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that allows   **
!**  to store replicated code for any File-IO routines not         **
!**  implemented in this file                                      **
!**                                                                **
    include 'file_io_common.inc'
!********************************************************************
endmodule File_io
