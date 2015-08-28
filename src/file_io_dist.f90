! $Id$
!
!  This module goes straight and implements distributed file Input/Output.
!  We use here only F95 features and are explicitly not HPC-friendly.
!
module File_io
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
      use Messages, only: fatal_error
!
      character(len=*)          , intent(in) :: file
      character(len=*), optional, intent(in) :: form
      logical,          optional, intent(in) :: remove_comments
      integer,          optional, intent(out):: nitems
!
      logical :: exists
!
!  Test if file exists.
!
      if (present(nitems)) nitems=0
!
      inquire(file=file,exist=exists)
      if (.not. exists) call fatal_error('parallel_open', 'file "'//trim(file)//'" not found')
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
endmodule File_io
