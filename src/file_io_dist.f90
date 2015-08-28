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
  integer, parameter :: parallel_unit = 14
!
  include 'file_io.h'

  private

  contains
!***********************************************************************
    subroutine parallel_open(file,form,remove_comments,nitems)
!
!  Read a global file.
!
!  18-mar-10/Bourdin.KIS: implemented
!
      character (len=*) :: file
      character (len=*), optional :: form
      logical, optional :: remove_comments
      integer, optional :: nitems
!
      logical :: exists
!
!  Test if file exists.
!
      inquire(file=file,exist=exists)
      if (.not. exists) call stop_it('parallel_open: file "'//trim(file)//'" not found')
!
!  Open file.
!
      if (present(form)) then
        open(parallel_unit, file=file, FORM=form, status='old')
      else
        open(parallel_unit, file=file, status='old')
      endif
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
      close(parallel_unit)
!
    endsubroutine parallel_close
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
!***********************************************************************
endmodule File_io
