! $Id$
!
!  This module takes care of massive parallel file Input/Output.
!  We use here only F95 features for HPC-friendly behaviour.
!
module File_io
!
  use Cdata
  use Cparam
  use Mpicomm
!
  implicit none
!
  external file_size_c
  external write_binary_file_c
!
  interface write_binary_file
    module procedure write_binary_file_str
    module procedure write_binary_file_carr
  endinterface
!
  integer, parameter :: parallel_unit = 14
  logical :: lignore_namelist_errors, lnamelist_error
!
  contains
!
!***********************************************************************
    function get_unit()
!
!  Returns the buffer for parallel reading.
!
!  29-May-2015/Bourdin.KIS: implemented
!
      integer :: get_unit
!
      get_unit = parallel_unit
!
    endfunction get_unit
!***********************************************************************
    function parallel_read(file,buffer,remove_comments)
!
!  Returns a buffer with the content of a global file read in parallel.
!
!  28-May-2015/Bourdin.KIS: implemented
!
      integer :: parallel_read
      character, dimension(:), allocatable :: buffer
      character (len=*), intent(in) :: file
      logical, optional :: remove_comments
!
      integer :: bytes, unit = 11
!
      if (lroot) then
        if (.not. file_exists(file)) call stop_it_if_any(.true., &
            'parallel_read: file not found "'//trim(file)//'"')
        bytes=file_size(file)
        if (bytes < 0) call stop_it_if_any(.true., &
            'parallel_read: could not determine file size "'//trim(file)//'"')
        if (bytes == 0) call stop_it_if_any(.true., &
            'parallel_read: file is empty "'//trim(file)//'"')
      endif
!
      ! Catch conditional errors of the MPI root rank.
      call stop_it_if_any(.false.,'')
!
      ! Broadcast the file size.
      call mpibcast_int(bytes)
      parallel_read = bytes
!
      ! Allocate temporary memory.
      if (allocated (buffer)) deallocate (buffer)
      allocate (buffer(bytes))
      buffer(1:) = char(0)
!
      if (lroot) then
        ! Read file content into buffer.
        open(unit, file=file, FORM='unformatted', RECL=bytes, ACCESS='direct', status='old')
        read(unit, REC=1) buffer
        close(unit)
!
        if (present (remove_comments)) then
          ! Strip comments.
          if (remove_comments) call strip_comments(buffer)
        endif
      endif
!
      ! Broadcast buffer to all MPI ranks.
      call mpibcast_char(buffer, bytes)
!
    endfunction parallel_read
!***********************************************************************
    subroutine parallel_open(unit,file,form,remove_comments)
!
!  Open a global file in parallel.
!
!  17-mar-10/Bourdin.KIS: implemented
!  28-May-2015/Bourdin.KIS: reworked
!
      use Cparam, only: fnlen
!
      integer, intent(out) :: unit
      character (len=*) :: file
      character (len=*), optional, intent(in) :: form
      logical, optional, intent(in) :: remove_comments
!
      integer :: bytes, pos
      character (len=fnlen) :: filename
      character, dimension(:), allocatable :: buffer
      character (len=fnlen) :: tmp_prefix
!
      bytes = parallel_read(file, buffer, remove_comments)
!
      ! Create unique temporary filename.
      pos=scan(file, '/')
      do while(pos /= 0)
        file(pos:pos)='_'
        pos=scan(file, '/')
      enddo
      tmp_prefix = get_tmp_prefix()
      write(filename,'(A,A,A,I0)') trim(tmp_prefix), file, '-', iproc
!
      ! Write temporary file into local RAM disk (/tmp).
      call write_binary_file(filename, bytes, buffer)
      deallocate(buffer)
!
      ! Open temporary file.
      if (present(form)) then
        open(parallel_unit, file=filename, FORM=form, status='old')
      else
        open(parallel_unit, file=filename, status='old')
      endif
      unit = parallel_unit
      ! unit is now reading from a local filesystem (probably /tmp)
      ! and is ready to be used on all ranks in parallel.
!
    endsubroutine parallel_open
!***********************************************************************
    subroutine parallel_rewind(unit)
!
!  Rewind a file unit opened by parallel_open.
!
!  23-May-2014/Bourdin.KIS: implemented
!
      integer, intent(in) :: unit
!
      rewind(unit)
!
    endsubroutine parallel_rewind
!***********************************************************************
    subroutine parallel_close(unit)
!
!  Close a file unit opened by parallel_open and remove temporary file.
!
!  17-mar-10/Bourdin.KIS: implemented
!
      integer, intent(in) :: unit
!
      close(unit,status='delete')
!
    endsubroutine parallel_close
!***********************************************************************
    function parallel_count_lines(file,ignore_comments)
!
!  Determines in parallel the number of lines in a file.
!
!  Returns:
!  * Integer containing the number of lines in a given file
!  * -1 on error
!
!  23-mar-10/Bourdin.KIS: implemented
!  26-aug-13/MR: optional parameter ignore_comments added
!  28-May-2015/Bourdin.KIS: reworked
!
      integer :: parallel_count_lines
      character (len=*), intent(in) :: file
      logical, optional, intent(in) :: ignore_comments
!
      if (lroot) parallel_count_lines = count_lines(file,ignore_comments)
      call mpibcast_int(parallel_count_lines)
!
    endfunction
!***********************************************************************
    function parallel_file_exists(file, delete)
!
!  Determines in parallel if a given file exists.
!  If delete is true, deletes the file.
!
!  Returns:
!  * Integer containing the number of lines in a given file
!  * -1 on error
!
!  23-mar-10/Bourdin.KIS: implemented
!
      character (len=*) :: file
      logical :: parallel_file_exists
      logical, optional :: delete
!
      logical :: ldelete
!
      if (present(delete)) then
        ldelete = delete
      else
        ldelete = .false.
      endif
!
      ! Let the root node do the dirty work
      if (lroot) parallel_file_exists = file_exists(file,ldelete)
!
      call mpibcast_logical(parallel_file_exists)
!
    endfunction
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
      character (len=*) :: file
      logical, optional :: delete
!
      integer :: unit = 1
!
      inquire (file=file, exist=file_exists)
!
      if (file_exists .and. present(delete)) then
        if (delete) then
          if (ip <= 6) print *, 'remove_file: Removing file <'//trim(file)//'>'
          open (unit, file=file)
          close (unit, status='delete')
        endif
      endif
!
    endfunction file_exists
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
!  23-may-2015/Bourdin.KIS: coded
!
      integer :: file_size
      character (len=*) :: file
!
      file_size = -2
      if (file_exists (trim(file)//char(0))) then
        file_size = -1
        call file_size_c(trim(file)//char(0), file_size)
      endif
!
    endfunction file_size
!***********************************************************************
    subroutine write_binary_file_carr(file,bytes,buffer)
!
!  Writes a given buffer (vector of single characters) to a binary file.
!
!  06-Apr-2014/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: bytes
      character, dimension(:), intent(in) :: buffer
!
      integer :: result
!
      call write_binary_file_c(trim(file)//char(0), bytes, buffer, result)

      if (result /= bytes) then
        if (result < 0) then
          print *, 'write_binary_file: could not open file for writing "'//trim(file)//'"'
        elseif (result == 0) then
          print *, 'write_binary_file: could not start writing "'//trim(file)//'"'
        else
          print *, 'write_binary_file: could not finish writing "'//trim(file)//'"', result
        endif
        stop
      endif
!
    endsubroutine write_binary_file_carr
!***********************************************************************
    subroutine write_binary_file_str(file,bytes,buffer)
!
!  Writes a given buffer (string) to a binary file.
!
!  21-jan-2015/MR: copied from write_binary_file_carr 
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: bytes
      character (len=*), intent(in) :: buffer
!
      integer :: result
!
      call write_binary_file_c(trim(file)//char(0), bytes, buffer, result)
      
      if (result /= bytes) then
        if (result < 0) then
          print *, 'write_binary_file: could not open file for writing "'//trim(file)//'"'
        elseif (result == 0) then
          print *, 'write_binary_file: could not start writing "'//trim(file)//'"'
        else
          print *, 'write_binary_file: could not finish writing "'//trim(file)//'"', result
        endif
        stop
      endif
!
    endsubroutine write_binary_file_str
!***********************************************************************
    function count_lines(file,ignore_comments)
!
!  Determines the number of lines in a file.
!
!  Returns:
!  * Integer containing the number of lines in a given file
!  * -1 on error
!
!  23-mar-10/Bourdin.KIS: implemented
!  26-aug-13/MR: optional parameter ignore_comments added
!  28-May-2015/Bourdin.KIS: reworked
!
      use General, only: operator(.in.)
!
      integer :: count_lines
      character (len=*), intent(in) :: file
      logical, optional, intent(in) :: ignore_comments
!
      integer :: ierr, unit = 1
      character :: ch
!
      count_lines = -1
      if (.not. file_exists(file)) return
!
      open (unit, file=file, status='old', iostat=ierr)
      if (ierr /= 0) return
      count_lines = 0
      do while (ierr == 0)
        read (unit,'(a)',iostat=ierr) ch
        if (ierr == 0) then
          if (present(ignore_comments)) then
            if (ignore_comments .and. (ch .in. (/ '!', comment_char /))) cycle
          endif
          count_lines = count_lines + 1
        endif
      enddo
      close (unit)
!
    endfunction count_lines
!**************************************************************************
    subroutine strip_comments(buffer)
!
!  Strip comments from a *.in-file
!
!  28-May-2015/Bourdin.KIS: inspired by MR's read_infile
!
      character, dimension(:), allocatable :: buffer
!
      integer :: num_bytes, pos
      logical :: lcomment
!
      if (.not. allocated (buffer)) return
!
      num_bytes = size (buffer)
      lcomment = .false.
      do pos = 1, num_bytes
        if (buffer(pos) == char(10)) then
          lcomment = .false.
        elseif ((buffer(pos) == '!') .or. (buffer(pos) == comment_char)) then
          lcomment = .true.
        endif
        if (lcomment) buffer(pos) = ' '
      enddo
!
    endsubroutine strip_comments
!**************************************************************************
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
      use Syscalls, only: get_PID, get_env_var
!
      character (len=fnlen) :: get_tmp_prefix
      character (len=fnlen) :: tmp_dir
!
      call get_env_var('TMPDIR', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) call get_env_var('TEMP', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) call get_env_var('TMP', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) call get_env_var('TMP_DIR', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) call get_env_var('PBS_TEMP', tmp_dir)
      if (len(trim(tmp_dir)) <= 0) tmp_dir = '/tmp'
!
      write (get_tmp_prefix,'(A,A,I0,A)') trim(tmp_dir), '/pencil-', get_PID(), '-'
!
    endfunction get_tmp_prefix
!***********************************************************************
    subroutine read_namelist(reader,name,lactive)
!
!  encapsulates reading of pars + error handling
!
!  31-oct-13/MR: coded
!  16-dec-13/MR: handling of optional ierr corrected
!  18-dec-13/MR: changed handling of ierr to avoid compiler trouble
!  19-dec-13/MR: changed ierr into logical lierr to avoid compiler trouble
!  11-jul-14/MR: end-of-file caught to avoid program crash when a namelist is missing
!  18-aug-15/PABourdin: reworked to simplify code and display all errors at once
!  19-aug-15/PABourdin: renamed from 'read_pars' to 'read_namelist'
!
      use Messages, only: warning
!
      interface
        subroutine reader(iostat)
          integer, intent(out) :: iostat
        endsubroutine reader
      endinterface
!
      character(len=*), intent(in) :: name
      logical, intent(in), optional :: lactive
!
      integer :: ierr
      character(len=linelen) :: type
      include "parallel_unit.h"
!
      if (present (lactive)) then
        if (.not. lactive) return
      endif
!
      call reader(ierr)
!
      if (lstart) then
        type = '_init'
        if (name == 'init') type = ''
      else
        type = '_run'
        if (name == 'run') type = ''
      endif
!
      if (.not. lignore_namelist_errors) then
        if (ierr == -1) then
          if (lroot) call warning ('read_namelist', 'namelist "'//trim(name)//trim(type)//'_pars" missing!')
          lnamelist_error = .true.
        elseif (ierr /= 0) then
          if (lroot) call warning ('read_namelist', 'namelist "'//trim(name)//trim(type)//'_pars" has an error!')
          lnamelist_error = .true.
        endif
      endif
!
      call parallel_rewind(parallel_unit)
!
    endsubroutine read_namelist
!**************************************************************************
endmodule File_io
