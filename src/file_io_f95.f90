! $Id$
!
!  This module takes care of massive parallel file Input/Output.
!  We use here only F95 features for HPC-friendly behaviour.
!
!  5-Aug-15/MR: removed unit parameters from parallel_open, parallel_rewind, 
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
    module procedure write_binary_file_char
    module procedure write_binary_file_str
  endinterface
!
  integer, parameter :: parallel_unit=14, parallel_unit_vec=14

  include 'file_io.h'
!
  private

  character, dimension(:), allocatable :: buffer

  contains
!***********************************************************************
    function parallel_read(file,remove_comments)
!
!  Fills buffer with the content of a global file. Buffer is read in parallel.
!  Returns length of buffer.
!
!  28-May-2015/Bourdin.KIS: implemented
!
      use General, only: loptest
      use Messages, only: fatal_error

      integer :: parallel_read
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
      call mpibcast_int(bytes,comm=MPI_COMM_WORLD)
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
        if (loptest(remove_comments)) call strip_comments
      endif
!
      ! Broadcast buffer to all MPI ranks.
      call mpibcast_char(buffer, bytes, comm=MPI_COMM_WORLD)
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
      use Cdata, only: iproc_world
!
      character (len=*),           intent(in) :: file
      character (len=*), optional, intent(in) :: form
      logical,           optional, intent(in) :: remove_comments
      integer,           optional, intent(out):: nitems
!
      integer :: bytes, pos
      character (len=fnlen) :: filename, tmp_prefix
!
      if (present(nitems)) nitems=0

      bytes = parallel_read(file, remove_comments)
!
      ! Create unique temporary filename.
      filename = file
      pos = scan(filename, '/')
      do while(pos /= 0)
        filename(pos:pos) = '_'
        pos = scan(filename, '/')
      enddo
      tmp_prefix = get_tmp_prefix()
      tmp_prefix = 'data/'
      write(filename,'(A,I0)') trim(tmp_prefix)//trim(filename)//'-', iproc_world
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
    subroutine write_binary_file_char(file,bytes,buffer)
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
          print *, 'write_binary_file_char: could not open file for writing "'//trim(file)//'"'
        elseif (result == 0) then
          print *, 'write_binary_file_char: could not start writing "'//trim(file)//'"'
        else
          print *, 'write_binary_file_char: could not finish writing "'//trim(file)//'"', result
        endif
        stop
      endif
!
    endsubroutine write_binary_file_char
!***********************************************************************
    subroutine write_binary_file_str(file,bytes,buffer)
!
!  Writes a given buffer (string) to a binary file.
!
!  21-jan-2015/MR: copied from write_binary_file_char 
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
          print *, 'write_binary_file_str: could not open file for writing "'//trim(file)//'"'
        elseif (result == 0) then
          print *, 'write_binary_file_str: could not start writing "'//trim(file)//'"'
        else
          print *, 'write_binary_file_str: could not finish writing "'//trim(file)//'"', result
        endif
        stop
      endif
!
    endsubroutine write_binary_file_str
!***********************************************************************
    subroutine strip_comments
!
!  Strip comments from a *.in-file
!
!  28-May-2015/Bourdin.KIS: inspired by MR's read_infile
!
      use Cdata, only: comment_char

      integer :: num_bytes, pos
      logical :: lcomment
!
      if (.not. allocated(buffer)) return

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
    !function find_namelist(name) result(lfound)
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
      use Mpicomm, only: lroot, mpibcast
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
!
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
