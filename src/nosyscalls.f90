! $Id: $
!
!  This module takes care of system calls and provides ANSI-C functionality.
!
module Syscalls
!
  implicit none
!
  contains

!***********************************************************************
    function file_exists(file, delete)
!
!  Determines if a file exists.
!  If delete is true, deletes the file.
!
!  Returns:
!  * Logical containing the existence of a given file
!
!  23-mar-10/Bourdin.KIS: implemented
!
      use Cparam, only : ip
!
      implicit none
!
      logical :: file_exists
      character(len=*) :: file
      logical, optional :: delete
!
      integer :: unit=1
!
      inquire(file=file, exist=file_exists)
!
      if (file_exists .and. present(delete)) then
        if (delete) then
          if (ip <= 6) print *, 'remove_file: Removing file <'//trim(file)//'>'
          open(unit,FILE=file)
          close(unit,STATUS='DELETE')
        endif
      endif
!
    endfunction
!***********************************************************************
    subroutine touch_file(file)
!
!  Touches a given file (used for code locking)
!  25-may-03/axel: coded
!  24-mar-10/Bourdin.KIS: moved here from sub.f90 and mpicomm.f90
!
      character (len=*) :: file
!
      integer :: unit=1
!
      open(unit,FILE=file)
      close(unit)
!
    endsubroutine touch_file
!***********************************************************************
    function file_size(file)
!
!  Determines the size of a given file.
!
!  Returns:
!  * positive integer containing the file size of a given file
!  * -2 if the file could not be found or opened
!  * -1 if retrieving the file size failed
!
!  18-mar-10/Bourdin.KIS: implemented
!
      implicit none
!
      character(len=*) :: file
      integer :: file_size
!
      integer :: ierr, unit=1
      logical :: exists
      integer, parameter :: buf_len=128
      character (len=buf_len) :: chunk
      integer :: n_chunks, trim_len
!
      ierr=0
      exists=.false.
!
      ! file must exist
      file_size=-2
      inquire(file=file, exist=exists)
      if (.not. exists) return
!
      ! open file and determine its size by reading chunks until EOF
      file_size=-1
      open(unit, FILE=file, FORM='unformatted', RECL=buf_len, ACCESS='direct', STATUS='old', IOSTAT=ierr)
      if (ierr /= 0) return
!
      n_chunks=0
      file_size=0
      trim_len=0
      do while(ierr == 0)
        chunk = char(0)
        n_chunks=n_chunks+1
        read(unit, REC=n_chunks, IOSTAT=ierr) chunk
        if (ierr == 0) then
          file_size=file_size+buf_len
          trim_len=len(trim(chunk))
        endif
      enddo
      close(unit)
!
      ! calculate file size and allocate a buffer
      file_size=file_size-buf_len+trim_len
      if (file_size < 0) file_size=-1
!
    endfunction file_size
!***********************************************************************
    function count_lines(file)
!
!  Determines the number of lines in a file
!
!  Returns:
!  * Integer containing the number of lines in a given file
!  * -1 on error
!
!  23-mar-10/Bourdin.KIS: implemented
!
      implicit none
!
      character(len=*) :: file
      integer :: count_lines
!
      integer :: unit=1, ierr
!
      count_lines=-1
      if (.not. file_exists(file)) return
!
      count_lines=-1
      open(unit, FILE=file, STATUS='old', IOSTAT=ierr)
      if (ierr /= 0) return
      do while (ierr == 0)
        read(unit,*,iostat=ierr)
        if (ierr == 0) count_lines=count_lines+1
      enddo
      close(unit)
!
    endfunction
!***********************************************************************
endmodule Syscalls

