! $Id$
!
!  This module takes care of massive parallel file Input/Output.
!  We use here only F95 features for HPC-friendly behaviour.
!
!  5-Aug-15/MR: shifted parallel_file_exists and parallel_count_lines to Sub,
!               file_exists, file_size and count_lines to General
!               removed unit parameters from parallel_open, parallel_rewind, 
!               parallel_close: it is always parallel_unit.
!               made parallel_unit public, hence removed get_unit.
!               moved lnamelist_error to cdata.
!
module File_io
!
  use Mpicomm
!
  implicit none
!
  external write_binary_file_c
!
  interface write_binary_file
    module procedure write_binary_file_str
    module procedure write_binary_file_carr
  endinterface
!
  integer, parameter :: parallel_unit = 14

  include 'file_io.h'
!
  private
  contains
!***********************************************************************
    function parallel_read(file,buffer,remove_comments)
!
!  Fills buffer with the content of a global file. Buffer is read in parallel.
!  Returns length of buffer.
!
!  28-May-2015/Bourdin.KIS: implemented
!
      use General, only: loptest, file_size, file_exists
      use Messages, only: fatal_error

      integer :: parallel_read
      character, dimension(:), allocatable :: buffer
      character (len=*), intent(in) :: file
      logical,           intent(in),  optional :: remove_comments
!
      integer :: bytes
      integer, parameter :: unit = 11
!
      if (lroot) then
        if (.not. file_exists(file)) call fatal_error( &
            'parallel_read', 'file "'//trim(file)//'" not found', force=.true.)
        bytes=file_size(file)
        if (bytes < 0) call fatal_error( &
            'parallel_read', 'could not determine size of file "'//trim(file)//'"', force=.true.)
        if (bytes == 0) call fatal_error( &
            'parallel_read', 'file "'//trim(file)//'" is empty', force=.true.)
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
      buffer = char(0)
!
      if (lroot) then
        ! Read file content into buffer.
        open(unit, file=file, FORM='unformatted', RECL=bytes, ACCESS='direct', status='old')
        read(unit, REC=1) buffer
        close(unit)
!
        ! Strip comments.
        if (loptest(remove_comments)) call strip_comments(buffer)
      endif
!
      ! Broadcast buffer to all MPI ranks.
      call mpibcast_char(buffer, bytes)
!
    endfunction parallel_read
!***********************************************************************
    subroutine parallel_open(file,form,remove_comments,nitems)
!
!  Open a global file in parallel.
!
!  17-mar-10/Bourdin.KIS: implemented
!  28-May-2015/Bourdin.KIS: reworked
!   5-Aug-2015/MR: added (dummy) parameter nitems

      use Cparam, only: fnlen
!
      character (len=*)                       :: file
      character (len=*), optional, intent(in) :: form
      logical,           optional, intent(in) :: remove_comments
      integer,           optional, intent(out):: nitems
!
      integer :: bytes, pos
      character (len=fnlen) :: filename
      character, dimension(:), allocatable :: buffer
      character (len=fnlen) :: tmp_prefix
!
      if (present(nitems)) nitems=0

      bytes = parallel_read(file, buffer, remove_comments)
!
      ! Create unique temporary filename.
      pos=scan(file, '/')
      do while(pos /= 0)
        file(pos:pos)='_'
        pos=scan(file, '/')
      enddo
      tmp_prefix = get_tmp_prefix()
      write(filename,'(A,I0)') trim(tmp_prefix)//trim(file)//'-', iproc
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
      ! parallel_unit is now reading from a local filesystem (probably /tmp)
      ! and is ready to be used on all ranks in parallel.
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
      close(parallel_unit,status='delete')
!
    endsubroutine parallel_close
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
    subroutine strip_comments(buffer)
!
!  Strip comments from a *.in-file
!
!  28-May-2015/Bourdin.KIS: inspired by MR's read_infile
!
      use Cdata, only: comment_char

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
!**************************************************************************
endmodule File_io
