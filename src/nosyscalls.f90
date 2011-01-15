! $Id$
!
!  This module takes care of system calls and provides ANSI-C functionality.
!
module Syscalls
!
  implicit none
!
  interface is_nan
     module procedure is_nan_0D
  endinterface
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
      use Cparam, only: ip
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
!  18-mar-10/Bourdin.KIS: implemented
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
      do while(ierr==0)
        chunk = char(0)
        n_chunks=n_chunks+1
        read(unit, REC=n_chunks, IOSTAT=ierr) chunk
        if (ierr==0) then
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
!  Determines the number of lines in a file.
!
!  Returns:
!  * Integer containing the number of lines in a given file
!  * -1 on error
!
!  23-mar-10/Bourdin.KIS: implemented
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
      do while (ierr == 0)
        read(unit,*,iostat=ierr)
        if (ierr==0) count_lines=count_lines+1
      enddo
      close(unit)
!
    endfunction count_lines
!***********************************************************************
    function get_PID()
!
!  The Fortran95 standard has no means to fetch the real PID.
!  If one need the real PID, please use the 'syscalls' module.
!
!   4-aug-10/Bourdin.KIS: coded
!
      integer :: get_PID
!
      get_PID = -1
!
      print *, 'get_PID: Not available in the F95 standard, please use the syscalls module.'
      stop
!
    endfunction get_PID
!***********************************************************************
    subroutine get_env_var(name,value)
!
!  Reads in an environment variable.
!
!  Returns:
!  * String containing the content of a given environment variable name
!  * Empty string, if the variable doesn't exist
!
!   4-aug-10/Bourdin.KIS: implemented
!
      character(len=*) :: name
      character(len=*) :: value
!
      value = char(0)
!
      print *, 'get_env_var('//name//'): Not available in the F95 standard, please use the syscalls module.'
      stop
!
    endsubroutine get_env_var
!***********************************************************************
    function get_tmp_prefix()
!
!  Determines the proper temp directory and adds a unique prefix.
!
!  Returns:
!  * String containing the location of a usable temp directory
!  * Default is '/tmp'
!
!   4-aug-10/Bourdin.KIS: coded
!
      use Cparam, only: fnlen
!
      character(len=fnlen) :: get_tmp_prefix
!
      get_tmp_prefix = '/tmp/pencil-'
!
    endfunction get_tmp_prefix
!***********************************************************************
    function is_nan_0D(value)
!
!  Determines if value is not a number (NaN).
!  This function is a trick to circumvent the lack of isnan in F95.
!
!  Usage of the syscalls module is highly recommended, because there,
!  the implementation if is_nan is made with ANSI standard routines.
!
!  Returns:
!  * true, if value is not a number (NaN)
!  * false, otherwise
!
!  14-jan-2011/Bourdin.KIS: coded
!
      logical :: is_nan_0D
      real, intent(in) :: value
!
      is_nan_0D = .not. (abs (value) <= huge (value))
!
    endfunction is_nan_0D
!***********************************************************************
endmodule Syscalls
