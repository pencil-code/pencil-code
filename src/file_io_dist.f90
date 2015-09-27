! $Id$
!
!  This module goes straight and implements distributed file Input/Output.
!  We use here only F95 features and are explicitly not HPC-friendly.
!
module File_io
!
  implicit none
!
  integer, parameter :: parallel_unit = 14
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
    function find_namelist(name)
!
!  Tests if the namelist is present and reports a missing namelist.
!
!  26-Sep-2015/PABourdin: coded
!
      use Cdata, only: iproc, comment_char
      use Cparam, only: linelen
      use General, only: lower_case
      use Mpicomm, only: lroot, mpibcast
      use Messages, only: warning
!
      logical :: find_namelist
      character(len=*), intent(in) :: name
!
      integer :: ierr, pos, state, max_len, line_len
      character(len=linelen) :: line
      character :: ch
!
      find_namelist = .false.
!
      if (lroot) then
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
                    find_namelist = .true.
                    exit
                  endif
                  ch = lower_case (line(pos+1:pos+1))
                  if ((ch == ' ') .or. (ch == char(9)) .or. (ch == '!') .or. (ch == comment_char)) then
                    find_namelist = .true.
                  endif
                  if (find_namelist) exit
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
        call parallel_rewind()
        if (.not. find_namelist) call warning ('find_namelist', 'namelist "'//trim(name)//'" is missing!')
      endif
!
      call mpibcast (find_namelist)
!
    endfunction find_namelist
!***********************************************************************
endmodule File_io
