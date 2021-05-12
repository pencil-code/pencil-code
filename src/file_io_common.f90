! $Id$
!
!  This module holds only common standard code for all File-IO modules.
!  It is NOT meant to be used as a replacement for the other File-IO modules!
!
!  26-Nov-2016/PABourdin: coded
!
module File_io
!
  implicit none
!
  private
!
  contains
!***********************************************************************
    subroutine parallel_rewind
!
    endsubroutine parallel_rewind
!***********************************************************************
!***********************************************************************
! The following routines are for later usage (transferred from General and Sub):
!***********************************************************************
!***********************************************************************
    subroutine fseek_pos(unit, rec_len, num_rec, reference)
!
!  Seeks to a given position in an opened file relative to a reference point.
!  If reference=0, this is relative to the beginning of the file,
!  if reference=1, this is relative to the current position in the file,
!  and if reference=2, this is relative to the end of the file.
!  'rec_len' and 'num_rec' are referring to a record length and a given number
!  of records that one likes to seek, boths must be representable in integer.
!  If 'num_rec' is negative, seeking is done backwards.
!
!  20-Feb-2012/PABourdin: coded
!
      use General, only: itoa, keep_compiler_quiet
      use Messages, only: fatal_error
!
      integer, intent(in) :: unit
      integer(kind=8) :: rec_len, num_rec
      integer, intent(in) :: reference
!
      integer :: i, num, len
!
      if (num_rec < 0) then
        num_rec = -num_rec
        rec_len = -rec_len
      endif
!
      ! all numbers must be representable as integer(kind=4)
      len = rec_len
      num = num_rec
      if (len /= rec_len) call fatal_error ('fseek_pos on unit '//trim (itoa (unit)), &
          "rec_len is not representable as integer(kind=4).", .true.)
      if (num /= num_rec) call fatal_error ('fseek_pos on unit '//trim (itoa (unit)), &
          "num_rec is not representable as integer(kind=4).", .true.)
!
! WORKAROUND:
! Even though the ifort manual states that ifort would be able to fseek
! with an 64-bit integer argument, this is NOT working!
! Therefore, we have to iterate the fseek with a 32-bit integer to be save.
! Note: gfortran would be able to seek with a 64-bit integer value, though.
! (20-Feb-2012, PABourdin)
!
      !!!call fseek (unit, rec_len, reference)
      if (num >= 2) then
        do i = 2, num
          !!!call fseek (unit, rec_len, 1)
          call keep_compiler_quiet(reference)
        enddo
      endif
!
    endsubroutine fseek_pos
!***********************************************************************
    subroutine backskip_to_time(lun,lroot)
!
!  Skips over possible persistent data blocks from end of snapshot to time record.
!
!  9-mar-15/MR: coded
!  8-sep-15/MR: excluded false detection of id_block_PERSISTENT for double precision version.
!               (in single precision false detection is impossible as id_block_PERSISTENT=2000
!                corresponds to the non-normalized real 2.80259693E-42)
!
      use Cparam, only: rkind8
      use General, only: loptest
!
      integer,           intent(in) :: lun
      logical, optional, intent(in) :: lroot
!
      integer :: i,id,ios
      real :: x
!
      backspace(lun)
      read(lun) id
!
      ios=1
      if (id==id_block_PERSISTENT) then
!
        backspace(lun)
        if (kind(x)==rkind8) then      ! if PC is in double precision version
          read(lun,iostat=ios) x       ! try to read a double precision number from the same position as id
          backspace(lun)
        endif
!
        if (ios/=0) then               ! if read try not done or unsuccessful: id_block_PERSISTENT was properly found
          do
            do i=1,3; backspace(lun); enddo
            read(lun) id
            if (id==id_block_PERSISTENT) exit
          enddo
          backspace(lun)
        endif
      endif
!
      if (ios/=0) backspace(lun)         ! if read try successful (ios==0), i.e., id_block_PERSISTENT was falsely detected,
                                         ! one backspace already done
      if (loptest(lroot)) backspace(lun)
!
    endsubroutine backskip_to_time
!***********************************************************************
    subroutine delete_file(file)
!
!  Deletes a file. Needed on CRAYs as status='replace' in open is not sufficient
!  to avoid unwanted file growth.
!
! 11-jan-15/MR: coded
!
      character(len=*), intent(in) :: file
!
      integer, parameter :: lun=111
      logical :: exists
!
      inquire(FILE=file, EXIST=exists)
      if (exists) then
        open (lun, FILE=file)
        close(lun, status='delete')
      endif
!
    endsubroutine delete_file
!****************************************************************************
    subroutine file_remove(file)
!
!  Removes a file if it exists.
!
!  05-Nov-2018/PABourdin: coded
!
      character(len=*), intent(in) :: file
!
      logical :: removed
!
      removed = file_exists(file, delete=.true.)
!
    endsubroutine file_remove
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
      use Cdata, only: ip
      use General, only: loptest
!
      logical :: file_exists
      character(len=*), intent(in) :: file
      logical, optional, intent(in) :: delete
!
      integer, parameter :: unit = 1
!
      inquire (file=file, exist=file_exists)
!
      if (file_exists .and. loptest(delete)) then
        if (ip <= 6) print *, 'file_exists: Removing file <'//trim(file)//'>'
        open (unit, file=file)
        close (unit, status='delete')
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
      character (len=*), intent(in) :: file
!
      file_size = -2
      if (file_exists(file)) then
        file_size = -1
        call file_size_c(trim(file)//char(0), file_size)
      endif
!
    endfunction file_size
!***********************************************************************
    function count_lines(file,ignore_comments,lbinary)
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
!   1-Dec-2017/MR: option for binary files added
!
      use Cdata, only: comment_char
      use General, only: loptest, operator (.in.)
!
      integer :: count_lines
      character (len=*), intent(in) :: file
      logical, optional, intent(in) :: ignore_comments, lbinary
!
      integer :: ierr, idum
      integer, parameter :: unit = 1
      character :: ch
      logical :: lbin
!
      count_lines = -1
      if (.not. file_exists(file)) return
!
      lbin=loptest(lbinary)
      if (lbin) then
        open (unit, file=file, status='old', iostat=ierr, form='unformatted')
      else
        open (unit, file=file, status='old', iostat=ierr)
      endif
      if (ierr /= 0) return
      count_lines = 0
      do while (ierr == 0)
        if (lbin) then
          read (unit,iostat=ierr) idum
        else
          read (unit,'(a)',iostat=ierr) ch
        endif
        if (ierr == 0) then
          if (loptest(ignore_comments) .and. (ch .in. (/ '!', comment_char /))) cycle
          count_lines = count_lines + 1
        endif
      enddo
      close (unit)
!
    endfunction count_lines
!***********************************************************************
    function list_files(name,options,only_number) result (num)
! 
! Wrapper for UNIX command ls. At present returns only number of found files.
!
! 20-may-18/MR: coded
!
      use General, only: coptest
      use Syscalls, only: system_cmd

      integer :: num
      
      character(LEN=*),           intent(IN) :: name
      character(LEN=*), optional, intent(IN) :: options 
      logical,          optional, intent(IN) :: only_number
    
      call system_cmd('ls '//coptest(options)//name//' > tmplsout 2> /dev/null')
      num=count_lines('tmplsout')
      call delete_file('tmplsout')
    
    endfunction list_files
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
      use Cdata, only: lroot
      use Mpicomm, only: mpibcast_int, MPI_COMM_WORLD
!
      integer :: parallel_count_lines
      character (len=*), intent(in) :: file
      logical, optional, intent(in) :: ignore_comments
!
      if (lroot) parallel_count_lines = count_lines(file,ignore_comments)
      call mpibcast_int(parallel_count_lines,comm=MPI_COMM_WORLD)
!
    endfunction parallel_count_lines
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
      use Cdata, only: lroot
      use General, only: loptest
      use Mpicomm, only: mpibcast_logical, MPI_COMM_WORLD
!
      logical :: parallel_file_exists
      character (len=*) :: file
      logical, optional :: delete
!
      ! Let the root node do the dirty work
      if (lroot) parallel_file_exists = file_exists(file,loptest(delete))
!
      call mpibcast_logical(parallel_file_exists,comm=MPI_COMM_WORLD)
!
    endfunction parallel_file_exists
!***********************************************************************
    subroutine read_namelist(reader,name,lactive)
!
!  Encapsulates reading of pars + error handling.
!
!  31-oct-13/MR: coded
!  16-dec-13/MR: handling of optional ierr corrected
!  18-dec-13/MR: changed handling of ierr to avoid compiler trouble
!  19-dec-13/MR: changed ierr into logical lierr to avoid compiler trouble
!  11-jul-14/MR: end-of-file caught to avoid program crash when a namelist is missing
!  18-aug-15/PABourdin: reworked to simplify code and display all errors at once
!  19-aug-15/PABourdin: renamed from 'read_pars' to 'read_namelist'
!
      use Cdata, only: lnamelist_error, lparam_nml, lstart, lroot
      use General, only: loptest, itoa
      use Messages, only: warning
!
      interface
        subroutine reader(iostat)
          integer, intent(out) :: iostat
        endsubroutine reader
      endinterface
!
      character(len=*), intent(in) :: name
      logical, optional, intent(in) :: lactive
!
      integer :: ierr
      logical :: found
      character(len=5) :: type, suffix
!
      if (.not. loptest (lactive, .true.)) return
!
      if (lstart .or. lparam_nml) then
        type = 'init'
      else
        type = 'run'
      endif
      if (name /= '') type = '_'//type
      suffix = '_pars'
      if (name == 'initial_condition_pars') then
        type = ''
        suffix = ''
      endif
!
      !if (.not. find_namelist (trim(name)//trim(type)//trim(suffix))) then
      call find_namelist (trim(name)//trim(type)//trim(suffix),found)
!
      ierr = 0 ! G95 complains 'ierr' is used but not set, even though 'reader' has intent(out).
      call reader(ierr)
!
      if (ierr /= 0) then
      
        if (.not.found) then
          if (.not. lparam_nml) lnamelist_error = .true.
          call parallel_rewind
          return
        endif
!
        lnamelist_error = .true.
        if (lroot) then
          if (ierr == -1) then
            call warning ('read_namelist', 'namelist "'//trim(name)//trim(type)//trim(suffix)//'" is missing!')
          else
            call warning ('read_namelist', 'namelist "'//trim(name)//trim(type)//trim(suffix)// &
                          '" has an error ('//trim(itoa(ierr))//')!')
          endif
        endif
      endif
!
      call parallel_rewind
!
    endsubroutine read_namelist
!***********************************************************************
    subroutine read_zaver(f,k1,k2,source,nav,indav,nstart_,ltaver)
    
      use Cparam, only: nx,ny,nz,l1,l2,m1,m2,n1,n2,lactive_dimension
      use Cdata, only: directory_snap
      use General, only: directory_names_std, loptest, ioptest
      use Messages, only: warning

      real, dimension(:,:,:,:), intent(OUT) :: f
      integer, intent(IN) :: k1,k2,nav,indav
      character(LEN=*) :: source
      integer, optional, intent(IN) :: nstart_
      logical, optional, intent(IN) :: ltaver

      integer :: k,nt,it,nstart,ios,klen
      logical :: s0
      integer, parameter :: unit=111
      real, dimension(nx,ny,nav) :: read_arr
      real :: tav, tav0
      real, dimension(:,:,:), allocatable :: buffer

      klen=k2-k1+1
      call directory_names_std(.true.)

      !if (file_exists(trim(directory_snap)//'/zaverages0.dat')) then
      !  open(unit,FILE=trim(directory_snap)//'/zaverages0.dat',FORM='unformatted', STATUS='old')
      if (file_exists(trim(source)//trim(directory_snap)//'/zaverages.dat')) then
        nstart=ioptest(nstart_,-1)
        allocate(buffer(nx,ny,klen))
        if (loptest(ltaver)) then
          open (unit, FILE=trim(source)//trim(directory_snap)//'/zaverages.dat', &
                FORM='unformatted', status='old')
          buffer=0.
        else
          open (unit, FILE=trim(source)//trim(directory_snap)//'/zaverages.dat', &
                FORM='unformatted', status='old', position='append')
          backspace(unit)
          backspace(unit)
        endif

        ios=0; s0=.true.; nt=0; it=1
        do while(ios==0)

          read(unit,iostat=ios) tav
          if (ios/=0) exit
          if (loptest(ltaver) .and. it<nstart) then
            read(unit,iostat=ios) tav
          else
            if (s0) then
              tav0=tav
              s0=.false.
            endif
            read(unit,iostat=ios) read_arr
            if (ios==0) then
              if (loptest(ltaver)) then
                nt=nt+1
                buffer=buffer+read_arr(:,:,indav:indav+klen-1)
              else
                buffer=read_arr(:,:,indav:indav+klen-1)
              endif
            endif
          endif
          it=it+1

        enddo
        close(unit)

        if (loptest(ltaver)) buffer=buffer/nt

        if (.not.lactive_dimension(3)) then
          f(l1:l2,m1:m2,n1,k1:k2) = buffer
        else
          do k=k1,k2
            f(l1:l2,m1:m2,n1:n2,k) = spread(buffer(:,:,k),3,nz)
          enddo
        endif
      else
        call warning('read_zaver','file "'// &
                     trim(source)//trim(directory_snap)//'/zaverages.dat" does not exist')
      endif

    endsubroutine read_zaver
!***********************************************************************
endmodule File_io
