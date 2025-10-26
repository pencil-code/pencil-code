! $Id$
!
!  This module takes care of system calls and provides ANSI-C functionality.
!
module Syscalls
!
  use Cparam, only: ikind8
  use Geometrical_types, only: torus_rect

  implicit none
!
  !integer, parameter :: ikind8=selected_int_kind(14)  ! 8-byte integer kind
  integer, parameter :: rkind8=selected_real_kind(12) ! 8-byte real kind

  external directory_exists_c
  external get_PID_c
  external get_env_var_c
  external is_nan_c
  external system_c
  external sizeof_real_c
  external copy_addr_c
  external copy_addr_c_bool
  external extract_string_c
  external mem_usage_c
!
  interface is_nan
    module procedure is_nan_0D
    module procedure is_nan_1D
    module procedure is_nan_2D
    module procedure is_nan_3D
    module procedure is_nan_4D
  endinterface
!
  interface copy_addr
    module procedure copy_addr_int
    module procedure copy_addr_int_1D
    module procedure copy_addr_int_2D
    module procedure copy_addr_log
    module procedure copy_addr_log_1D
    module procedure copy_addr_real
    module procedure copy_addr_real_1D
    module procedure copy_addr_real_2D
    module procedure copy_addr_real_3D
    module procedure copy_addr_real_4D
    module procedure copy_addr_torus_rect
    !module procedure copy_addr_dble
    !module procedure copy_addr_dble_1D
  endinterface
!
  interface
    function sleep_c(ms) bind (c, name="sleep")
      use iso_c_binding, only : c_int
      integer(c_int), value :: ms
      integer(c_int) :: sleep_c
    endfunction sleep_c
  endinterface
!
  contains
!***********************************************************************
    subroutine system_cmd(command)
!
!  Executes a system command.
!
!  3-nov-11/MR: coded
!
      character(len=*) :: command
!
      integer :: ret
      integer :: system_c

      ret=system_c(trim(command)//char(0))
!
    endsubroutine system_cmd
!***********************************************************************
    subroutine sleep(seconds)
!
!  Sleep a pre-defined time.
!
!  26-Oct-2025/PABourdin: coded
!
      use iso_c_binding, only: c_int
!
      integer, intent(in) :: seconds
!
      integer(c_int) :: ms_c
      integer :: result
!
      ms_c = seconds * 1000
      result = sleep_c (ms_c)
!
    endsubroutine sleep
!***********************************************************************
    function sizeof_real()
!
!  Determines the size of a real in bytes.
!
!  Returns:
!  * The number of bytes used for a real.
!
!  16-Feb-2012/Bourdin.KIS: coded
!
      integer :: sizeof_real
!
      call sizeof_real_c(1.0, sizeof_real)
!
    endfunction sizeof_real
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
      call get_environment_variable (trim (name), value)
      value = trim (value)
!
    endsubroutine get_env_var
!***********************************************************************
    function directory_exists(path)
!
!  Checks for existence of a directory.
!
!  Returns:
!  * True, if 'path' points to a directory
!  * False, otherwise
!
!   2-sep-15/PABourdin: coded
!
      logical :: directory_exists
      character(len=*) :: path
!
      integer :: exists
!
      exists = -1

      call directory_exists_c(trim(path)//char(0), exists)
      if (exists == -1) then
        write (*,*) 'WARNING: failure while checking if "'//trim(path)//'" exists!'
        ! This line is temporary code to allow Nils debugging further.
        !exists = 1
      endif
!
      directory_exists = (exists == 1)
!
    endfunction directory_exists
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
    logical function readlink(filename,link)
!
!  Returns the file pointed to by symbolic link filename in link.
!  Returns .true. if successful.
!
!  21-mar-20/MR: coded
!
      character(LEN=*), intent(IN) :: filename
      character(LEN=*), intent(OUT):: link

      integer :: readlink_c

      integer :: retlen, exists

      readlink=.false.

      retlen=readlink_c(trim(filename)//char(0),link,len(link))
      if (retlen==-1) return

      link(retlen+1:len(link))=''

      call directory_exists_c(trim(link)//char(0), exists)
      if (exists==-1) return

      readlink=.true.

    endfunction readlink 
!***********************************************************************
    logical function islink(filename)
!
!  Tests whether filename is a symbolic link.
!
!  21-mar-20/MR: coded
!
      character(LEN=*), intent(IN) :: filename
      integer :: islink_c

      islink = islink_c(trim(filename)//char(0))==1

    endfunction islink
!***********************************************************************
    subroutine extract_str(cmd,result)
!
!  Extracts a string by cmd, e.g., from a file (would be included in cmd)
!  and returns it in result.
!
!  21-mar-20/MR: coded
!
    character(LEN=*), intent(IN) :: cmd
    character(LEN=*), intent(OUT):: result

    result(len(result):)=char(0)

    call extract_string_c(trim(cmd)//char(0),result,len(result))

    result(index(result,char(0))-1:) = ''  ! as last char is newline

    endsubroutine extract_str
!***********************************************************************
    subroutine copy_addr_int(var, caddr)

    integer, intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_int
!***********************************************************************
    subroutine copy_addr_log(var, caddr)

    logical, intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_log
!***********************************************************************
    subroutine copy_addr_real_1D(var, caddr)

    real, dimension(:), intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_real_1D
!***********************************************************************
    subroutine copy_addr_log_1D(var, caddr)

    logical, dimension(:), intent(IN) :: var
    integer, dimension(size(var)) :: tmp
    integer(KIND=ikind8), intent(OUT) :: caddr
    
    integer :: i

    tmp = 0
    do i = 1,size(var)
        if(var(i)) tmp(i) = 1
    enddo

    call copy_addr_c_bool(tmp,caddr,size(var))

    endsubroutine copy_addr_log_1D
!***********************************************************************
    subroutine copy_addr_int_1D(var, caddr)

    integer, dimension(:), intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_int_1D
!***********************************************************************
    subroutine copy_addr_int_2D(var, caddr)

    integer, dimension(:,:), intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_int_2D
!***********************************************************************
    subroutine copy_addr_real_2D(var, caddr)

    real, dimension(:,:), intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_real_2D
!***********************************************************************
    subroutine copy_addr_real_3D(var, caddr)

    real, dimension(:,:,:), intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_real_3D
!***********************************************************************
    subroutine copy_addr_real_4D(var, caddr)

    real, dimension(:,:,:,:), intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_real_4D
!***********************************************************************
    subroutine copy_addr_dble(var, caddr)

    real(KIND=rkind8), intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_dble
!***********************************************************************
    subroutine copy_addr_dble_1D(var, caddr)

    real(KIND=rkind8), dimension(:), intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_dble_1D
!***********************************************************************
    subroutine copy_addr_real(var, caddr)

    real, intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_real
!***********************************************************************
    subroutine copy_addr_torus_rect(var, caddr)

    type(torus_rect), intent(IN) :: var
    integer(KIND=ikind8), intent(OUT) :: caddr

    call copy_addr_c(var,caddr)

    endsubroutine copy_addr_torus_rect
!***********************************************************************
    subroutine get_char_arr(strarr,strlen,arrlen)

    integer :: strlen,arrlen
    character(LEN=strlen), dimension(arrlen), intent(IN) :: strarr

    !call get_char_arr_c(strarr,strlen,arrlen)

    endsubroutine get_char_arr
!***********************************************************************
    integer function memusage()

      integer :: mem_usage_c

      memusage=mem_usage_c()

    endfunction memusage
!***********************************************************************
endmodule Syscalls
