! $Id$
!
!  This module takes care of system calls and provides ANSI-C functionality.
!
module Syscalls
!
  implicit none
!
  external file_size_c
  external get_pid_c
  external get_env_var_c
  external is_nan_c
  external system_c
!
  interface is_nan
     module procedure is_nan_0D
     module procedure is_nan_1D
     module procedure is_nan_2D
     module procedure is_nan_3D
     module procedure is_nan_4D
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
      use Cparam, only : ip
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
    function get_PID()
!
!  Determines the PID of the current process.
!
!  Returns:
!  * Integer containing the PID of the current process
!  * -1 if retrieving of the PID failed
!
!   4-aug-10/Bourdin.KIS: coded
!
      integer :: get_PID
!
      integer, save :: my_PID = -1
!
      if (my_PID == -1) call get_PID_c(my_PID)
      get_PID = my_PID
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
!   4-aug-10/Bourdin.KIS: coded
!
      character(len=*) :: name
      character(len=*) :: value
!
      value = ' '
      call get_env_var_c(trim(name)//char(0), value)
      value = trim(value)
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
      character(len=fnlen) :: tmp_dir
!
      call get_env_var('TMPDIR', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) call get_env_var('TEMP', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) call get_env_var('TMP', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) call get_env_var('TMP_DIR', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) call get_env_var('PBS_TEMP', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) call get_env_var('PBS_O_LOCAL', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) tmp_dir = '/tmp'
!
      write (get_tmp_prefix,'(A,A,I0,A)') trim(tmp_dir), '/pencil-', get_PID(), '-'
!
    endfunction get_tmp_prefix
!***********************************************************************
    subroutine system(command)
!
!  launches system command.
!
!  3-nov-11/MR: coded
!
      character(len=*) :: command
!
      call system_c(trim(command)//char(0))
!
    endsubroutine system
!***********************************************************************
    function is_nan_0D(value)
!
!  Determines if value is not a number (NaN).
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
      integer :: result
!
      call is_nan_c (value, result)
      if (result == -1) then
        print *, 'is_nan_0D: serious failure in is_nan_c'
        stop
      endif
      is_nan_0D = (result == 1)
!
    endfunction is_nan_0D
!***********************************************************************
    function is_nan_1D(value)
!
!  Determines if value is not a number (NaN).
!
!  Returns:
!  * true, if value is not a number (NaN)
!  * false, otherwise
!
!  15-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:), intent(in) :: value
      logical, dimension(size (value, 1)) :: is_nan_1D
!
      integer, dimension(size (value, 1)) :: result
      integer :: pos
!
      do pos = 1, size (value, 1)
        call is_nan_c (value(pos), result(pos))
      enddo
      if (any (result == -1)) then
        print *, 'is_nan_1D: serious failure in is_nan_c'
        stop
      endif
      is_nan_1D = (result == 1)
!
    endfunction is_nan_1D
!***********************************************************************
    function is_nan_2D(value)
!
!  Determines if value is not a number (NaN).
!
!  Returns:
!  * true, if value is not a number (NaN)
!  * false, otherwise
!
!  15-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: value
      logical, dimension(size (value, 1),size (value, 2)) :: is_nan_2D
!
      integer, dimension(size (value, 1),size (value, 2)) :: result
      integer :: pos_x, pos_y
!
      do pos_x = 1, size (value, 1)
        do pos_y = 1, size (value, 2)
          call is_nan_c (value(pos_x,pos_y), result(pos_x,pos_y))
        enddo
      enddo
      if (any (result == -1)) then
        print *, 'is_nan_2D: serious failure in is_nan_c'
        stop
      endif
      is_nan_2D = (result == 1)
!
    endfunction is_nan_2D
!***********************************************************************
    function is_nan_3D(value)
!
!  Determines if value is not a number (NaN).
!
!  Returns:
!  * true, if value is not a number (NaN)
!  * false, otherwise
!
!  15-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: value
      logical, dimension(size (value, 1),size (value, 2),size (value, 3)) :: is_nan_3D
!
      integer, dimension(size (value, 1),size (value, 2),size (value, 3)) :: result
      integer :: pos_x, pos_y, pos_z
!
      do pos_x = 1, size (value, 1)
        do pos_y = 1, size (value, 2)
          do pos_z = 1, size (value, 3)
            call is_nan_c (value(pos_x,pos_y,pos_z), result(pos_x,pos_y,pos_z))
          enddo
        enddo
      enddo
      if (any (result == -1)) then
        print *, 'is_nan_3D: serious failure in is_nan_c'
        stop
      endif
      is_nan_3D = (result == 1)
!
    endfunction is_nan_3D
!***********************************************************************
    function is_nan_4D(value)
!
!  Determines if value is not a number (NaN).
!
!  Returns:
!  * true, if value is not a number (NaN)
!  * false, otherwise
!
!  15-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: value
      logical, dimension(size (value, 1),size (value, 2),size (value, 3),size (value, 4)) :: is_nan_4D
!
      integer, dimension(size (value, 1),size (value, 2),size (value, 3),size (value, 4)) :: result
      integer :: pos_x, pos_y, pos_z, pos_a
!
      do pos_x = 1, size (value, 1)
        do pos_y = 1, size (value, 2)
          do pos_z = 1, size (value, 3)
            do pos_a = 1, size (value, 4)
              call is_nan_c (value(pos_x,pos_y,pos_z,pos_a), result(pos_x,pos_y,pos_z,pos_a))
            enddo
          enddo
        enddo
      enddo
      if (any (result == -1)) then
        print *, 'is_nan_4D: serious failure in is_nan_c'
        stop
      endif
      is_nan_4D = (result == 1)
!
    endfunction is_nan_4D
!***********************************************************************
endmodule Syscalls
