! $Id$
!
!  This module goes straight and implements distributed file Input/Output.
!  We use here only F2003 features (HPC-friendly).
!
!  5-Aug-15/MR: shifted parallel_file_exists and parallel_count_lines to Sub,
!               file_exists, file_size and count_lines to General
!               removed unit parameters from parallel_open, parallel_rewind, 
!               parallel_close: it is always parallel_unit, now an allocatable
!               string array.
!               Made parallel_unit public (but protected) hence removed get_unit.
!
module File_io
!
  implicit none
!
! Fixed string length necessary as gfortran is compiling incorrectly otherwise.
! For future debugged gfortran versions the commented lines should be used.
!
  character(LEN=:), dimension(:), allocatable, protected :: parallel_unit
  !character(LEN=36000), dimension(:), allocatable, protected :: parallel_unit
!
  include 'file_io.h'

  private

  contains
!***********************************************************************
    subroutine parallel_read(file,remove_comments,nitems)
!
!  Provides the (comment-purged) contents if an input file in parallel_unit.
!  Reading from file system is done by root only.
!  Returns number of valid records in nitems, if requested.
!
!  28-May-2015/Bourdin.KIS: implemented
!  31-Jul-2015/MR: adapted to allocatable parallel_unit, 
!                  reinstated comment removal as in former read_infile.
!                  (seems to be necessary)
!
      use Mpicomm, only: lroot, mpibcast_int, mpibcast_char
      use General, only: loptest, parser, file_exists, file_size
      use Messages, only: fatal_error
      use Cdata, only: comment_char, root

      character (len=*), intent(in)            :: file
      logical,           intent(in),  optional :: remove_comments
      integer,           intent(out), optional :: nitems
!
      integer :: bytes, ios, ind, indc, inda, inda2, lenbuf, indmax, ni
      integer, parameter :: unit = 11
      character(LEN=14000) :: linebuf             ! string length overdimensioned, but needed so for some compilers.
      character(LEN=:), allocatable :: buffer
      character :: sepchar
      logical :: l0
!
      if (lroot) then
        if (.not. file_exists(file)) call fatal_error(&
            'parallel_read', 'file "'//trim(file)//'" not found', force=.true.)
        bytes = file_size(file)
        if (bytes < 0) call fatal_error(&
            'parallel_read', 'could not determine size of file "'//trim(file)//'"', force=.true.)
        if (bytes == 0) call fatal_error(&
            'parallel_read', 'file "'//trim(file)//'" is empty', force=.true.)

        ! Allocate temporary memory.
        allocate(character(len=bytes) :: parallel_unit(1))
        !allocate(parallel_unit(1))
!
        ! Read file content into buffer.
        open(unit, file=file, status='old')

        ! Namelist-read and non-namelist-read files need to be treated differently:
        ! For the former a blank, for the latter, LF is inserted between the records
        ! and the number of (non-comment) records is counted in nitems.
        
        if (present(nitems)) then
          sepchar=char(10)
        else
          sepchar=' '
        endif

        l0=.true.; ni=0; indmax=0
        do
          read(unit,'(a)',iostat=ios) linebuf
          if (ios<0) exit

          linebuf=adjustl(linebuf)
!
          if (loptest(remove_comments)) then
            inda=index(linebuf,"'")                                  ! position of an opening string bracket
            ind=index(linebuf,'!'); indc=index(linebuf,comment_char) ! positions of comment indicators

            ! check if comment indicators are within a string, hence irrelevant.

            if (inda>0) then 
              inda2=index(linebuf(inda+1:),"'")+inda    ! position of a closing string bracket
              if (inda2==inda) inda2=len(linebuf)+1     ! if closing string bracket missing, assume it at end of record
              if (ind>inda.and.ind<inda2) ind=0
              if (indc>inda.and.indc<inda2) indc=0
            endif

            if (indc>0) ind=min(max(ind,1),indc)        ! determine smaller of the two comment indicators
          else
            ind=0
          endif
!
          if (ind==0) then
            ind=len(trim(linebuf))
          else
            ind=ind-1
            if (ind>0) ind=len(trim(linebuf(1:ind)))  ! if comment appended, remove it
          endif
          indmax = max(indmax,ind)                    ! update maximum length of record
!
          if (ind==0) then                            ! line is a comment or empty -> skip
            cycle
          elseif (l0) then                            ! otherwise append it to parallel_unit with 
            parallel_unit(1)=linebuf(1:ind)           ! separating character sepchar
            lenbuf=ind
            l0=.false.
          else
            parallel_unit(1)=parallel_unit(1)(1:lenbuf)//sepchar//linebuf(1:ind)
            lenbuf=lenbuf+ind+1
          endif
          ni=ni+1                                     ! update number of valid records

        enddo
        close(unit)
!
      endif
!
      ! Broadcast the size of parallel_unit.
      call mpibcast_int(lenbuf)
!
      if (present(nitems)) then

        ! broadcast number of valid records and maximum record length
        if (lroot) nitems=ni
        call mpibcast_int(nitems)
        if (nitems==0) return
        ni=nitems
        call mpibcast_int(indmax)

      else
        ! let ni and indmax have valid values for namelist-read files
        ni=1; indmax=lenbuf
      endif

      ! prepare broadcasting of parallel_unit.
      if (lroot) then
        if (present(nitems)) then

          ! for non-namelist-read files: organize parallel_unit as array 
          allocate(character(len=lenbuf) :: buffer)
          buffer=parallel_unit(1)(1:lenbuf)
          deallocate(parallel_unit)
          allocate(character(len=indmax) :: parallel_unit(nitems))
          !allocate(parallel_unit(nitems))
          ! decompose former parallel_unit into records guided by sepchar
          if (parser(buffer,parallel_unit,sepchar)/=nitems) &
            call fatal_error('parallel_read', 'too less elements found when parsing buffer')

        endif 
      else  
        allocate(character(len=indmax) :: parallel_unit(ni))
        !allocate(parallel_unit(ni))
        parallel_unit=' '
      endif
!
      ! broadcast parallel_unit
      if (present(nitems)) then
        call mpibcast_char(parallel_unit, nitems, root)
      else  
        call mpibcast_char(parallel_unit(1)(1:lenbuf), root)
      endif
!
    endsubroutine parallel_read
!***********************************************************************
    subroutine parallel_open(file,form,remove_comments,nitems)
!
!  Read a global file.
!
!  18-mar-10/Bourdin.KIS: implemented
!  30-jul-15/MR: reworked
!
      character (len=*), intent(in)            :: file
      character (len=*), intent(in),  optional :: form
      logical,           intent(in),  optional :: remove_comments
      integer,           intent(out), optional :: nitems
!
!  Parameter form is ignored as parallel_read is at present implemented for formatted reading only.
!
      call parallel_read(file, remove_comments, nitems)
!
    endsubroutine parallel_open
!***********************************************************************
    subroutine parallel_rewind
!
!  Dummy.
!
!  30-Jul-2015/MR: implemented
!
    endsubroutine parallel_rewind
!***********************************************************************
    subroutine parallel_close
!
!  Deallocates the internal file parallel_unit opened by parallel_open.
!  Each call to parallel_open must be followed by a call to parallel_close
!  before parallel_open can be called again.

!  30-Jul-2015/MR: implemented
!
      deallocate(parallel_unit)
!
    endsubroutine parallel_close
!***********************************************************************
endmodule File_io
