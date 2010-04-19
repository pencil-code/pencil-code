! $Id$
!
!  This module takes care of system calls and provides ANSI-C functionality.
!
module Syscalls
!
  implicit none
!
  external file_size_c
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
    endfunction file_exists
!***********************************************************************
    subroutine touch_file(file)
!
!  Touches a given file (used for code locking).
!
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
!  19-mar-10/Bourdin.KIS: coded
!
      implicit none
!
      integer file_size
      character(len=*) :: file
!
      file_size=-1
      call file_size_c(trim(file)//char(0), file_size)
!
    endfunction file_size
!***********************************************************************
    function count_lines(file)
!
!  Determines the number of lines in a file.
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
      count_lines=0
      open(unit, FILE=file, STATUS='old', IOSTAT=ierr)
      if (ierr/=0) return
      do while (ierr==0)
        read(unit,*,iostat=ierr)
        if (ierr==0) count_lines=count_lines+1
      enddo
      close(unit)
!
    endfunction count_lines
!***********************************************************************
endmodule Syscalls
