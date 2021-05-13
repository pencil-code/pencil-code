! $Id$
!
!  This module goes straight and implements distributed file Input/Output.
!  We use here only F2003 features (HPC-friendly).
!
!  5-Aug-15/MR: removed unit parameters from parallel_open, parallel_rewind, 
!               parallel_close: it is always parallel_unit, now an allocatable
!               string array.
!               Made parallel_unit public (but protected) hence removed get_unit.
!  5-Oct-15/MR: intro'd two parallel units: scalar/vector to avoid ugly big constant string length;
!               old code still in commented lines
!
module File_io
!
  implicit none
!
  ! Fixed length or scalar string necessary as gfortran is compiling incorrectly otherwise.
  ! For future debugged gfortran versions the commented lines should be used.
!
  !character(len=:), dimension(:), allocatable, protected :: parallel_unit                ! gfortran v4.9.2 will not compile correctly with this line
  character(len=:),                          allocatable, protected :: parallel_unit      ! gfortran v4.8.4 will not compile this line
  ! Temporary replacement code for the preceding line(has some other consequences):
  !character(len=36000), protected :: parallel_unit
  integer, parameter :: fixed_buflen=128
  character(len=fixed_buflen), dimension(:), allocatable, protected :: parallel_unit_vec
!
  include 'file_io.h'

  private

  contains
!***********************************************************************
    subroutine parallel_read(file,remove_comments,nitems,lbinary)
!
!  Provides the (comment-purged) contents if an input file in parallel_unit.
!  Reading from file system is done by root only.
!  Returns number of valid records in nitems, if requested.
!
!  28-May-2015/Bourdin.KIS: implemented
!  31-Jul-2015/MR: adapted to allocatable parallel_unit, 
!                  reinstated comment removal as in former read_infile.
!                  (seems to be necessary)
!   6-oct-15/MR: parameter lbinary added for reading data as a byte stream
!
      use Mpicomm, only: lroot, mpibcast_int, mpibcast_char, MPI_COMM_WORLD
      use General, only: loptest, parser
      use Messages, only: fatal_error
      use Cdata, only: comment_char

      character (len=*), intent(in)            :: file
      logical,           intent(in),  optional :: remove_comments
      integer,           intent(out), optional :: nitems
      logical,           intent(in),  optional :: lbinary
!
      integer :: bytes, ios, ind, indc, inda, inda2, lenbuf, indmax, ni
      integer, parameter :: unit = 11
      character(len=14000) :: linebuf             ! string length overdimensioned, but needed so for some compilers.
      !character(len=:), allocatable :: buffer    ! g95 v0.92 will not compile this line
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

        ! Allocate internal file parallel_unit.
        allocate(character(LEN=bytes) :: parallel_unit)
!
        ! Read file content into buffer.
        if (loptest(lbinary)) then
          open(unit, file=file, status='old',access='stream')
          read(unit) parallel_unit
          lenbuf=bytes
        else
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
              parallel_unit=linebuf(1:ind)           ! separating character sepchar
              lenbuf=ind
              l0=.false.
            else
              parallel_unit=parallel_unit(1:lenbuf)//sepchar//linebuf(1:ind)
              lenbuf=lenbuf+ind+1
            endif
            ni=ni+1                                     ! update number of valid records

          enddo
        endif
        close(unit)
      endif
!
      if (present(nitems)) then

        ! broadcast number of valid records and maximum record length
        if (lroot) nitems=ni
        call mpibcast_int(nitems,comm=MPI_COMM_WORLD)
        if (nitems==0) return
        !call mpibcast_int(indmax)
        allocate(parallel_unit_vec(nitems))

      else
      ! Broadcast the size of parallel_unit.
        call mpibcast_int(lenbuf,comm=MPI_COMM_WORLD)
      endif

      ! prepare broadcasting of parallel_unit.
      if (lroot) then
        if (present(nitems)) then

          ! for non-namelist-read files: organize parallel_unit as array 
          !allocate(character(len=lenbuf) :: buffer)   ! gfortran v4.6.3 will not compile this line, v4.8.4 works
          ! Temporary replacement code for the above line:
          !buffer = (repeat (char (0), lenbuf))

          !buffer=parallel_unit(1)(1:lenbuf)
          !deallocate(parallel_unit)
          !allocate(character(len=indmax) :: parallel_unit(nitems))
          ! decompose former parallel_unit into records guided by sepchar
          !if (parser(buffer,parallel_unit_vec,sepchar)/=nitems) &
          if (parser(parallel_unit,parallel_unit_vec,sepchar)/=nitems) &
            call fatal_error('parallel_read', 'too less elements found when parsing buffer')

        endif 
      else
        if (.not.present(nitems)) allocate(character(LEN=lenbuf) :: parallel_unit)
      endif
!
      ! broadcast parallel_unit
      if (present(nitems)) then
        call mpibcast_char(parallel_unit_vec, nitems, comm=MPI_COMM_WORLD)
      else  
        call mpibcast_char(parallel_unit(1:lenbuf), comm=MPI_COMM_WORLD)
      endif
!
    endsubroutine parallel_read
!***********************************************************************
    subroutine parallel_open(file,form,remove_comments,nitems,lbinary)
!
!  Read a global file.
!
!  18-mar-10/Bourdin.KIS: implemented
!  30-jul-15/MR: reworked
!   6-oct-15/MR: parameter lbinary added for reading data as a byte stream
!
      character (len=*), intent(in)            :: file
      character (len=*), intent(in),  optional :: form
      logical,           intent(in),  optional :: remove_comments
      integer,           intent(out), optional :: nitems
      logical,           intent(in),  optional :: lbinary
!
!  Parameter form is ignored as parallel_read is at present implemented for formatted reading only.
!
      call parallel_read(file, remove_comments, nitems, lbinary)
!
    endsubroutine parallel_open
!***********************************************************************
    subroutine parallel_close
!
!  Deallocates the internal file parallel_unit opened by parallel_open.
!  Each call to parallel_open must be followed by a call to parallel_close
!  before parallel_open can be called again.

!  30-Jul-2015/MR: implemented
!
      if (allocated(parallel_unit)) deallocate(parallel_unit)
      if (allocated(parallel_unit_vec)) deallocate(parallel_unit_vec)
!
    endsubroutine parallel_close
!***********************************************************************
    !function find_namelist(name) result(lfound)
    subroutine find_namelist(name,lfound)
!
!  Tests if the namelist is present and reports a missing namelist.
!
!  26-Sep-2015/PABourdin: coded
!   6-oct-2015/MR: turned into subroutine because of CRAY compiler bug;
!                  easily revertable by shifting comment char at beginning and end.

      use Cdata, only: comment_char
      use General, only: lower_case, operator(.in.)
      use Messages, only: warning
      use Mpicomm, only: lroot, mpibcast,MPI_COMM_WORLD
!
      character(len=*), intent(in) :: name
      logical :: lfound
!
      integer :: pos, len, max_len
!
      if (lroot) then
!print*, 'name=', name
        lfound = .false.
        len = len_trim (name) + 1
        ! need to subtract two chars for the end marker of an empty namelist
        max_len = len_trim (parallel_unit) - len + 1 - 2
        do pos = 1, max_len
          if ('&'//lower_case (trim (name)) == lower_case (parallel_unit(pos:pos+len-1))) then
!print*, 'line=',parallel_unit(pos:pos+len-1) 
            if (parallel_unit(pos+len:pos+len) .in. (/ ' ', '!', comment_char /)) then
              if (pos == 1) then
!print*, 'line=',pos,'#'//parallel_unit(pos+len:pos+len)//'#' 
                lfound = .true.
                exit
              elseif (parallel_unit(pos-1:pos-1) .eq. ' ') then
!print*, 'line=',pos,'#'//parallel_unit(pos-1:pos-1)//'#' 
                lfound = .true.
                exit
              endif
            endif
          endif
        enddo
        if (.not. lfound) call warning ('find_namelist', 'namelist "'//trim(name)//'" is missing!')
      endif
!
      call mpibcast (lfound,comm=MPI_COMM_WORLD)
!
    !endfunction find_namelist
    endsubroutine find_namelist
!***********************************************************************
    subroutine flush_file(unit)

      integer, intent(IN) :: unit

      flush(unit)

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
