! $Id$
!
!  This module takes care of massive parallel file Input/Output.
!  We use here F2003 features for HPC-friendly behaviour.
!
module File_io
!
  use Cdata
  use Cparam
  use Mpicomm
!
  implicit none
!
  character (len=:), allocatable, protected :: parallel_unit
!
  include 'file_io.h'
!
  private
!
  contains
!
!***********************************************************************
    function get_parallel_unit()
!
!  Returns the internal file buffer for parallel reading.
!
!  27-Sep-2015/PABourdin: implemented
!
      character (len=:), allocatable :: parallel_unit
!
      get_parallel_unit = parallel_unit
!
    endfunction parallel_unit
!***********************************************************************
    function parallel_read(file,buffer,remove_comments)
!
!  Returns a buffer with the content of a global file read in parallel.
!
!  28-May-2015/PABourdin: implemented
!
      integer :: parallel_read
      character (len=:), allocatable :: buffer
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
 write (*,*) 'PRA: ', iproc
 flush (6)
      ! Catch conditional errors of the MPI root rank.
      call stop_it_if_any(.false.,'')
!
 write (*,*) 'PRB: ', iproc
 flush (6)
      ! Broadcast the file size.
      call mpibcast_int(bytes)
      parallel_read = bytes
!
 write (*,*) 'PRC: ', iproc
 flush (6)
      ! Allocate temporary memory.
      if (allocated (buffer)) deallocate (buffer)
      allocate (character (len=bytes) :: buffer)
      buffer(1:) = char(0)
!
      if (lroot) then
        ! Read file content into buffer.
 write (*,*) 'PRD: ', iproc
 flush (6)
        open(unit, file=file, form='unformatted', access='stream', status='old')
        read(unit, REC=1) buffer
        close(unit)
!
 write (*,*) 'PRE: ', iproc
 flush (6)
        if (present (remove_comments)) then
          ! Strip comments.
          if (remove_comments) call strip_comments(buffer)
        endif
      endif
!
 write (*,*) 'PRF: ', iproc
 flush (6)
      ! Broadcast buffer to all MPI ranks.
      call mpibcast_char_alloc(buffer, bytes)
 write (*,*) 'PRG: ', iproc
 flush (6)
 write (*,*) iproc,file,'===>',buffer,'<==='
 flush (6)
 write (*,*) 'PRH: ', iproc
 flush (6)
!
    endfunction parallel_read
!***********************************************************************
    subroutine parallel_open(file,form,remove_comments,nitems)
!
!  Open a global file in parallel.
!
!  17-mar-10/PABourdin: implemented
!  28-May-2015/PABourdin: reworked
!
      use Cparam, only: fnlen
!
      character (len=*) :: file
      character (len=*), optional, intent(in) :: form
      logical, optional, intent(in) :: remove_comments
      integer, optional, intent(out) :: nitems
!
      integer :: bytes
!
 write (*,*) 'POA: ', iproc
 flush (6)
      bytes = parallel_read(file, parallel_unit, remove_comments)
      if (present (nitems)) nitems = 1
 write (*,*) 'POB: ', iproc
 flush (6)
      ! parallel_unit is now reading from an internal file (from RAM)
      ! and is ready to be used on all ranks in parallel.
!
    endsubroutine parallel_open
!***********************************************************************
    subroutine parallel_rewind()
!
!  Rewind a parallel file unit.
!
!  26-Sep-2015/PABourdin: implemented
!
    endsubroutine parallel_rewind
!***********************************************************************
    subroutine parallel_close()
!
!  Close a parallel file unit.
!
!  17-mar-10/PABourdin: implemented
!
      character (len=:), allocatable :: buffer
!
 write (*,*) 'PCA: ', iproc
 flush (6)
      parallel_buffer(1:) = char(0)
 write (*,*) 'PCB: ', iproc
 flush (6)
      deallocate (parallel_buffer)
 write (*,*) 'PCC: ', iproc
 flush (6)
!      if (allocated (buffer)) then
! write (*,*) "PC: DOUBLE DEALLOCATION!"
! flush (6)
!        deallocate (buffer)
!      endif
! write (*,*) 'PCD: ', iproc
! flush (6)
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
!  23-mar-10/PABourdin: implemented
!  26-aug-13/MR: optional parameter ignore_comments added
!  28-May-2015/PABourdin: reworked
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
!  23-mar-10/PABourdin: implemented
!
      character(len=*) :: file
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
!  23-mar-10/PABourdin: implemented
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
!  23-may-2015/PABourdin: coded
!
      integer :: file_size
      character(len=*) :: file
!
      file_size = -2
      if (file_exists (trim(file)//char(0))) then
        file_size = -1
        inquire (file=trim(file)//char(0), size=file_size)
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
!  23-mar-10/PABourdin: implemented
!  26-aug-13/MR: optional parameter ignore_comments added
!  28-May-2015/PABourdin: reworked
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
!  28-May-2015/PABourdin: inspired by MR's read_infile
!
      use Cdata, only: comment_char
!
      character (len=:), allocatable :: buffer
!
      integer :: num_bytes, pos
      character :: ch
      logical :: lcomment
!
      num_bytes = len_trim (buffer)
      lcomment = .false.
      do pos = 1, num_bytes
        ch = buffer(pos:pos)
        if (ch == char(10)) then
          lcomment = .false.
        elseif ((ch == '!') .or. (ch == comment_char)) then
          lcomment = .true.
        endif
        if (lcomment) buffer(pos:pos) = ' '
      enddo
!
    endsubroutine strip_comments
!**************************************************************************
    function find_namelist(name)
!
!  Tests if the namelist is present and reports a missing namelist.
!
!  26-Sep-2015/PABourdin: coded
!
      use Cdata, only: iproc, comment_char
      use General, only: operator(.in.)
      use Messages, only: warning
      use Mpicomm, only: lroot, mpibcast
!
      logical :: find_namelist
      character(len=*), intent(in) :: name
!
      integer :: ierr, state, pos, len, max_len
      character :: ch
!
      find_namelist = .false.
!
      if (lroot) then
        state = 0
        len = len_trim (name)
 write (10+iproc,*) 'FIND: ', '&'//name
 write (10+iproc,*) 'MAX: ', len_trim (parallel_unit)
 flush (10+iproc)
        ! need to substract one begin (&) and one end marker (/) chars
        max_len = len_trim (parallel_unit) - len + 1 - 2
        do pos = 1, max_len
          ch = parallel_unit(pos:pos)
          if (ch .eq. char(10)) state = 0
          if (ch .eq. '!') state = -2
          if ((ch .eq. ' ') .or. (ch .eq. char(9)) .or. (state < 0)) cycle
          if (ch .ne. '&') then
            state = -1
            cycle
          endif
          if (trim (name) == parallel_unit(pos+1:pos+len)) then
 write (10+iproc,*) 'FOUND POS: ', pos
 flush (10+iproc)
            find_namelist = .true.
            exit
          endif
        enddo
        if (.not. find_namelist) call warning ('find_namelist', 'namelist "'//trim(name)//'" is missing!')
      endif
!
      call mpibcast (find_namelist)
!
    endfunction find_namelist
!***********************************************************************
endmodule File_io
