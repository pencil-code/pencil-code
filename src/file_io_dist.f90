! $Id$
!
!  This module goes straight and implements distributed file Input/Output.
!  We use here only F95 features and are explicitly not HPC-friendly.
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
!
  integer, parameter :: parallel_unit = 14
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
        bytes = file_size(file)
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
      allocate(buffer(bytes))
      buffer=char(0)
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
!  Read a global file.
!
!  18-mar-10/Bourdin.KIS: implemented
!
      integer, intent(out) :: unit
      character (len=*) :: file
      character (len=*), optional :: form
      logical, optional :: remove_comments
!
      logical :: exists
!
!  Test if file exists.
!
      inquire(file=file,exist=exists)
      if (.not. exists) call stop_it('parallel_open: file not found "'//trim(file)//'"')
!
!  Open file.
!
      if (present(form)) then
        open(unit, file=file, FORM=form, status='old')
      else
        open(parallel_unit, file=file, status='old')
      endif
      unit = parallel_unit
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
      close(unit)
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
      character(len=*), intent(in) :: file
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
      logical :: parallel_file_exists
      character(len=*) :: file
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
      character(len=*) :: file
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
      character(len=*) :: file
!
      file_size = -2
      if (file_exists (trim(file)//char(0))) then
        file_size = -1
        call file_size_c(trim(file)//char(0), file_size)
      endif
!
    endfunction file_size
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
write (*,*) 'NEW: ==>', buffer, '<=='
!
    endsubroutine strip_comments
!***********************************************************************
endmodule File_io
