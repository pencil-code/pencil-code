! $Id$
!
!  This module takes care of code messages.
!
module Messages
!
  use Cdata
  use Mpicomm
!
  implicit none
!
  private
!
  public :: svn_id, timing
  public :: initialize_messages
  public :: information, warning, error
  public :: fatal_error, inevitably_fatal_error, not_implemented
  public :: fatal_error_local, fatal_error_local_collect
  public :: life_support_on, life_support_off
  public :: outlog
!
  integer, public, parameter :: iterm_DEFAULT   = 0
  integer, public, parameter :: iterm_BRIGHT    = 1
  integer, public, parameter :: iterm_UNDERLINE = 4
  integer, public, parameter :: iterm_FLASH     = 5
  integer, public, parameter :: iterm_FG_BLACK  = 30
  integer, public, parameter :: iterm_FG_RED    = 31
  integer, public, parameter :: iterm_FG_GREEN  = 32
  integer, public, parameter :: iterm_FG_YELLOW = 33
  integer, public, parameter :: iterm_FG_BLUE   = 34
  integer, public, parameter :: iterm_FG_MAGENTA= 35
  integer, public, parameter :: iterm_FG_CYAN   = 36
  integer, public, parameter :: iterm_FG_WHITE  = 37
  integer, public, parameter :: iterm_BG_BLACK  = 40
  integer, public, parameter :: iterm_BG_RED    = 41
  integer, public, parameter :: iterm_BG_GREEN  = 42
  integer, public, parameter :: iterm_BG_YELLOW = 43
  integer, public, parameter :: iterm_BG_BLUE   = 44
  integer, public, parameter :: iterm_BG_MAGENTA= 45
  integer, public, parameter :: iterm_BG_CYAN   = 46
  integer, public, parameter :: iip_EVERYTHING  = 0
  integer, public, parameter :: iip_DEFAULT     = 0
  integer, parameter :: iinformation_ip = 1000
  integer :: warnings=0
  integer :: errors=0
  integer :: fatal_errors=0, fatal_errors_total=0
  logical :: ldie_onwarning=.false.
  logical :: ldie_onerror=.true.
  logical :: ldie_onfatalerror=.true.
  logical :: llife_support=.false.
!
  logical :: ltermcap_color=.false.
!
  contains
!***********************************************************************
    subroutine initialize_messages
!
! Set a flag if colored output has been requested.
! Also set a flag if fake_parallel_io is requested.
!
      inquire(FILE="COLOR", EXIST=ltermcap_color)
      inquire(FILE="FAKE_PARALLEL_IO", EXIST=lfake_parallel_io)
!
    endsubroutine initialize_messages
!***********************************************************************
    subroutine not_implemented(location,message)
!
      character(len=*)           :: location
      character(len=*), optional :: message
!
      if (.not.llife_support) then
        errors=errors+1
!
        if (lroot .or. (ncpus<=16 .and. (message/=''))) then
          call terminal_highlight_error()
          write (*,'(A18)',ADVANCE='NO') "NOT IMPLEMENTED: "
          call terminal_defaultcolor()
          if (present(message)) then
            write(*,*) trim(location) // ": " // trim(message)
          else
            write(*,*) trim(location) // ": " // &
                "Some feature waits to get implemented -- by you?"
          endif
        endif
!
        if (ldie_onfatalerror) call die_gracefully
!
      endif
!
    endsubroutine not_implemented
!***********************************************************************
   subroutine fatal_error(location,message,force)
!
      character(len=*) :: location
      character(len=*) :: message
      logical, optional :: force
!
      logical :: fatal = .false.
!
      if (.not.llife_support) then
        errors=errors+1
        if (present(force)) fatal=force
!
        if (lroot .or. (ncpus<=16 .and. (message/='')) .or. fatal) then
          call terminal_highlight_fatal_error()
          write (*,'(A13)',ADVANCE='NO') "FATAL ERROR: "
          call terminal_defaultcolor()
          write (*,*) trim(location) // ": " // trim(message)
        endif
!
        if (ldie_onfatalerror .and. fatal) call die_immediately
        if (ldie_onfatalerror) call die_gracefully
!
      endif
!
    endsubroutine fatal_error
!***********************************************************************
    subroutine inevitably_fatal_error(location,message,force)
!
!  A fatal error that doesn't care for llife_support
!  Use (sparingly) in those cases where things should fail even during
!  pencil_consistency_test
!  07-26-2011: Julien\ Added forced exit if "force" is set to .true.
!
      character(len=*) :: location
      character(len=*) :: message
!
      logical, optional :: force
!
      logical :: fatal = .false.
!
      if (present(force)) fatal=force
      errors=errors+1
!
      if (lroot .or. (ncpus<=16 .and. (message/='')) .or. fatal) then
        call terminal_highlight_fatal_error()
        write (*,'(A13)',ADVANCE='NO') "FATAL ERROR: "
        call terminal_defaultcolor()
        write (*,*) trim(location) // ": " // trim(message)
      endif
!
      if (fatal) call die_immediately()
      call die_gracefully()
!
    endsubroutine inevitably_fatal_error
!***********************************************************************
    subroutine fatal_error_local(location,message)
!
!  Register a fatal error happening at one processor. The code will die
!  at the end of the time-step.
!
!  17-may-2006/anders: coded
!
      character(len=*) :: location
      character(len=*) :: message
!
      if (.not.llife_support) then
        fatal_errors=fatal_errors+1
!
        if (lroot .or. (ncpus<=16 .and. (message/=''))) then
          call terminal_highlight_fatal_error()
          write (*,'(A13)',ADVANCE='NO') "FATAL ERROR: "
          call terminal_defaultcolor()
          write (*,*) trim(location) // ": " // trim(message)
        endif
!
      endif
!
    endsubroutine fatal_error_local
!***********************************************************************
    subroutine fatal_error_local_collect()
!
!  Collect fatal errors from processors and die if there are any.
!
!  17-may-2006/anders: coded
!
      if (.not.llife_support) then
        call mpireduce_sum_int(fatal_errors,fatal_errors_total)
        call mpibcast_int(fatal_errors_total,1)
!
        if (fatal_errors_total/=0) then
          if (lroot) then
            print*, 'DYING - there were', fatal_errors_total, 'errors.'
            print*, 'This is probably due to one or more fatal errors that'
            print*, 'have occurred only on a single processor.'
          endif
          if (ldie_onfatalerror) call die_gracefully
        endif
!
        fatal_errors=0
        fatal_errors_total=0
!
      endif
!
    endsubroutine fatal_error_local_collect
!***********************************************************************
    subroutine error(location,message)
!
      character(len=*) :: location
      character(len=*) :: message
!
      if (.not.llife_support) then
        errors=errors+1
!
        if (lroot .or. (ncpus<=16 .and. (message/=''))) then
          call terminal_highlight_error()
          write (*,'(A7)',ADVANCE='NO') "ERROR: "
          call terminal_defaultcolor()
          write (*,*) trim(location) // ": " // trim(message)
        endif
!
        if (ldie_onerror) call die_gracefully
!
      endif
!
    endsubroutine error
!***********************************************************************
    subroutine warning(location,message,ip)
!
!  Print out colored warning.
!
!  30-jun-05/tony: coded
!
      character (len=*) :: message,location
      integer, optional :: ip
!
      if (.not.llife_support) then
        if (lroot .or. (ncpus<=16 .and. (message/=''))) then
          call terminal_highlight_warning()
          write (*,'(A9)',ADVANCE='NO') "WARNING:"
          call terminal_defaultcolor()
          write (*,*) trim(location) // ": " // trim(message)
!          call flush() ! has to wait until F2003
        endif
!
        if (ldie_onwarning) call die_gracefully
!
      endif
!
      if (NO_WARN) print*, ip
!
    endsubroutine warning
!***********************************************************************
    subroutine information(location,message,level)
!
!  Print out colored warning.
!
!  30-jun-05/tony: coded
!
      character (len=*) :: message,location
      integer, optional :: level
      integer :: level_ = iinformation_ip
!
      if (present(level)) level_=level
!
      if (ip<=level_) write (*,*) trim(location) // ": " // trim(message)
!
    endsubroutine information
!***********************************************************************
    subroutine svn_id(svnid)
!
!  Print SVN Revision info in a compact, yet structured form.
!  Expects the standard "SVN Id:" line as argument.
!
!  25-jun-02/wolf: coded
!
      character (len=*) :: svnid
!
      character (len=20) :: filename, revision, author, date
      character (len=200) :: fmt
      character (len=20) :: tmp1,tmp2,tmp3,tmp4
      integer :: if0,if1,iv0,iv1,iy0,iy1,it0,it1,ia0,ia1,iat
      integer :: wf=18, wv=7, wd=19 ! width of individual fields
      integer :: wd1=0
      logical, save :: lfirstcall=.true.
!
!  Write string to screen and to 'svnid.dat' file.
!
      if (lfirstcall) then
        open(1, file=trim(datadir)//'/svnid.dat', status='replace')
        lfirstcall=.false.
      else
        open(1, file=trim(datadir)//'/svnid.dat', status='old', position='append')
      endif
!
!  Construct format
!  Need to set explicit format below, to avoid problems when the
!  -i8 compiler option is invoked. Hope that the format i5 is sufficient.
!
      write(tmp1,'(i5)') wf
      write(tmp2,'(i5)') 6+wf
      write(tmp3,'(i5)') 6+wf+4+wv
      write(tmp4,'(i5)') 6+wf+4+wv+2+wd
!      fmt = '(A, A' // trim(adjustl(tmp1)) &
      fmt = '(A, A' &
           // ', T' // trim(adjustl(tmp2)) &
           // ', " v. ", A, T' // trim(adjustl(tmp3)) &
           // ', " (", A, T' // trim(adjustl(tmp4)) &
           // ', ") ", A)'
      if ((svnid(1:1) == "$") .and. (svnid(2:4) == "Id:")) then
      ! starts with `$...' --> SVN Id: line
!
!  file name
!
        if0 = index(svnid, ": ") + 2
        if1 = if0 + index(svnid(if0+1:), " ") - 1
        call extract_substring(svnid, if0, if1, filename)
!
!  Revision number
!
        iv0 = if1 + 2
        iv1 = iv0 + index(svnid(iv0+1:), " ") - 1
        call extract_substring(svnid, iv0, iv1, revision)
!
!  Date
!
        iy0 = iv1 + 2             ! first char of year
        iy1 = iy0 + 10            ! last char of year
        it0 = iy1 + 2             ! first char of time-of-day
        it1 = it0 + index(svnid(it0+1:), " ") - 1
        if (svnid(it1:it1) == "Z") then
          call extract_substring(svnid, iy0, it1-1, date) ! strip trailing `Z'
        else
          call extract_substring(svnid, iy0, it1, date)
        endif
!
!  Author
!
        ia0 = it1 + 2
        ! strip @some.where part off some user names
        iat = index(svnid(ia0+1:), "@")
        if (iat > 0) then
          ia1 = ia0 + iat - 1
        else
          ia1 = ia0 + index(svnid(ia0+1:), " ") - 1
        endif
        call extract_substring(svnid, ia0, ia1, author)
!
        write(*,fmt) "SVN: ", &
            trim(filename), &
            revision(1:wv), &
            date(1:wd), &
            trim(author)
!
        write(1,fmt) "SVN: ", &
            trim(filename), &
            revision(1:wv), &
            date(1:wd), &
            trim(author)
      else                      ! not a SVN line; maybe `[No ID given]'
        wd1 = min(wd, len(svnid))
        write(*,fmt) "SVN: ", &
            '-------', &
            '', &
            '', &
            svnid(1:wd1)
        write(1,fmt) "SVN: ", &
            '-------', &
            '', &
            '', &
            svnid(1:wd1)
      endif
      !write(*,'(A)') '123456789|123456789|123456789|123456789|123456789|12345'
      !write(*,'(A)') '         1         2         3         4         5'
!
      close(1)
!
    endsubroutine svn_id
!***********************************************************************
    subroutine timing(location,message,instruct,mnloop)
!
!  Timer: write the current systems time to standart output
!  provided it=it_timing.
!
      integer :: lun=9
      character(len=*) :: location
      character(len=*) :: message
      double precision :: time
      double precision, save :: time_initial
      character(len=*), optional :: instruct
      logical, optional :: mnloop
      integer :: mul_fac
!
!  work on the timing only when it==it_timing
!
      if (it==it_timing) then
        if (lroot) then
!
!  initialize
!
          if (present(instruct)) then
            if (trim(instruct)=='initialize') then
              open(lun,file=trim(datadir)//'/timing.dat', status='replace')
              time_initial=mpiwtime()
            endif
          endif
!
!  write current timing to the timing file
!
          if (lfirst) then
            if ((present(mnloop).and.lfirstpoint).or. .not.present(mnloop)) then
              time=mpiwtime()-time_initial
              if (present(mnloop)) then
                mul_fac=ny*nz
              else
                mul_fac=1
              endif
              write (lun,*) time,mul_fac,trim(location)//": "//trim(message)
            endif
          endif
!
!  finalize
!
          if (present(instruct)) then
            if (trim(instruct)=='finalize') close(lun)
          endif
!
        endif
      endif
!
    endsubroutine timing
!***********************************************************************
    subroutine extract_substring(string, idx0, idx1, substring)
!
!  Extract a substring after sanity check
!
      intent(in) :: string, idx0, idx1
      intent(out) :: substring
      character(len=*) :: string
      integer :: idx0, idx1
      character(len=*) substring
!
      if (1 <= idx0 .and. idx0 <= idx1 .and. idx1 <= len(string)) then
        substring = string(idx0:idx1)
      else
        substring = "???"
      endif
!
    endsubroutine extract_substring
!***********************************************************************
    subroutine life_support_off(message)
!
!  Allow code to die on errors
!
!  30-jun-05/tony: coded
!
      character(len=*) :: message
!
!  set llife_support
!
      llife_support=.false.
      call information('life_support_off',message,level=12)
!
    endsubroutine life_support_off
!***********************************************************************
    subroutine life_support_on(message)
!
!  Prevent the code from dying on errors
!
!  30-jun-05/tony: coded
!
      character(len=*) :: message
!
      call information('life_support_on',message,level=12)
      llife_support=.true.
!
    endsubroutine life_support_on
!***********************************************************************
    subroutine terminal_setfgcolor(col)
!
!  Set foreground color of terminal text
!
!  08-jun-05/tony: coded
!
      integer :: col
!
      if (ltermcap_color) then
        write(*,fmt='(A1,A1,I2,A1)',ADVANCE='no') CHAR(27), '[', col, 'm'
      endif
!
    endsubroutine terminal_setfgcolor
!***********************************************************************
    subroutine terminal_setfgbrightcolor(col)
!
!  Set bright terminal colors
!
!  08-jun-05/tony: coded
!
      integer :: col
!
      if (ltermcap_color) then
        write(*,fmt='(A1,A1,I1,A1,I2,A1)',ADVANCE='no') &
            CHAR(27), '[', iterm_BRIGHT, ';', col, 'm'
      endif
!
    endsubroutine terminal_setfgbrightcolor
!***********************************************************************
    subroutine terminal_defaultcolor
!
!  Set terminal color to default value
!
!  08-jun-05/tony: coded
!
      if (ltermcap_color) then
        write(*,fmt='(A1,A1,I1,A1)',ADVANCE='no') &
            CHAR(27), '[', iterm_DEFAULT, 'm'
      endif
!
    endsubroutine terminal_defaultcolor
!***********************************************************************
    subroutine terminal_highlight_warning
!
!  Change to warning color
!
!  08-jun-05/tony: coded
!
      if (ltermcap_color) then
        write(*,fmt='(A1,A1,I1,A1,I2,A1)',ADVANCE='no') &
            CHAR(27), '[', iterm_BRIGHT, ';', iterm_FG_MAGENTA, 'm'
      endif
!
    endsubroutine terminal_highlight_warning
!***********************************************************************
    subroutine terminal_highlight_error
!
!  Change to warning color
!
!  08-jun-05/tony: coded
!
      if (ltermcap_color) then
        write(*,fmt='(A1,A1,I1,A1,I2,A1)',ADVANCE='no') &
            CHAR(27), '[', iterm_BRIGHT, ';', iterm_FG_RED, 'm'
      endif
!
    endsubroutine terminal_highlight_error
!***********************************************************************
    subroutine terminal_highlight_fatal_error
!
!  Change to warning color
!
!  08-jun-05/tony: coded
!
      if (ltermcap_color) then
        write(*,fmt='(A1,A1,I1,A1,I2,A1)',ADVANCE='no') &
            CHAR(27), '[', iterm_BRIGHT, ';', iterm_FG_RED, 'm'
      endif
!
    endsubroutine terminal_highlight_fatal_error
!***********************************************************************
  subroutine outlog(code,msg,file)
!
!  Creates log entries for I/O errors in ioerrors.log.
!  Notifies user via e-mail if address mailaddress is given unless it
!  stops program if lstop_on_ioerror is set.
!
!  code(IN): errorcode from IOSTAT
!  msg (IN): describes failed action, starts with 'open', 'read', 'write' or 'close'
!  file(IN): file with which operation failed,
!            if omitted assumed to be the one saved in curfile
!
!  3-nov-11/MR: coded
!
    use Syscalls, only: system
    use General, only: itoa,date_time_string
!
    integer,                     intent(IN) :: code
    character (LEN=*),           intent(IN) :: msg
    character (LEN=*), optional, intent(IN) :: file
!
    character (LEN=80), save :: curfile=''
    integer :: unit=87, iostat
    character (LEN=20) :: date, codestr, strarr(2)
!
    if (present(file)) curfile = file
!
    if (code == 0) return
!
    errormsg = 'ERROR'
!
    if ( msg(1:5)=='open '  .or. msg(1:5)=='read ' .or. &
         msg(1:6)=='write ' .or. msg(1:6)=='close ' ) then
!
      errormsg = trim(errormsg)//' when '//msg(1:4)//'ing'
      if ( msg(1:4)=='read' ) errormsg = trim(errormsg)//trim(msg(6:))//' from'
      errormsg = trim(errormsg)//' file "'
!
    else
      errormsg = trim(errormsg)//': file "'
    endif
!
    codestr = itoa(code)
    errormsg = trim(errormsg)//trim(curfile)//'". Code: '//trim(codestr)
 
    if ( trim(mailaddress) == '' ) then
      lstop_on_ioerror = .true.
    else if (lroot) then 
!    
! scan of ioerrors.log to avoid multiple entries for the same file with same error code.
! When user eliminates cause of error, (s)he should also remove the corresponding line(s) in ioerror.log.
!
           strarr(1) = trim(curfile)
           strarr(2) = trim(codestr)
!
           if ( .not.scanfile('ioerrors.log',2,strarr,'all') ) then
!
             open(unit,file='ioerrors.log',position='append',iostat=IOSTAT)
             if (iostat==0) then
               call date_time_string(date)
               write(unit,'(a)',iostat=IOSTAT) date//' '//trim(errormsg)
               close(unit,iostat=IOSTAT)
             endif
!
             call system( &
              'echo '//trim(errormsg)//'| mail -s PencilCode Message '//trim(mailaddress) )
!
           endif
!
         endif
!
    if (lstop_on_ioerror) call stop_it('Stopped due to: '//trim(errormsg))
!
  endsubroutine outlog
!***********************************************************************
  logical function scanfile(file,nstr,strings,mode)
!
!  Scans a file for a line in which one or all strings in list strings occur.
!  Returns on first hit.
!
!  3-nov-11/MR: coded
!
    character (LEN=*),                           intent(IN) :: file
    integer,                                     intent(IN) :: nstr
    character (LEN=*), dimension(nstr),          intent(IN) :: strings
    character (LEN=3),                 optional, intent(IN) :: mode
!
    character (LEN=3) :: model
    character (LEN=120) :: line
    integer :: lun=111,i,count,iostat
   
    if ( .not.present(mode) ) then
      model='any'
    else
      model=mode
    endif
!
    scanfile = .false.
    open(lun,file=file,ERR=99,IOSTAT=iostat) 
!
    do
      read(lun,'(a)',ERR=97,END=98) line
!
      count=0
      do i=1,nstr
        if ( index(line,trim(strings(i)))/=0 ) then
          if (mode=='any') then
            scanfile = .true.
            goto 98
          else
            count = count+1
          endif
        endif
      enddo
!
      if ( count==nstr ) then
        scanfile = .true.
        goto 98 
      endif 
      
97  enddo
 
98  close(lun,ERR=99)
99  continue

  end function scanfile
!***********************************************************************
    subroutine input_array(file,a,dimx,dimy,dimz,dimv)
!
!  Generalized form of input, allows specifying dimension.
!
!  27-sep-03/axel: coded
!  04-nov-11/MR: moved here from General
!
      character (len=*) :: file
      integer :: dimx,dimy,dimz,dimv
      real, dimension (dimx,dimy,dimz,dimv) :: a
!
      integer :: iostat
!
      open(1,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (iostat /= 0) call stop_it("Cannot open "//trim(file)//" for reading",iostat)
      read(1,IOSTAT=iostat) a
      if (iostat /= 0) call stop_it("Cannot read a from "//trim(file),iostat)
      close(1,IOSTAT=iostat)
      if (iostat /= 0) call stop_it("Cannot close "//trim(file),iostat)
!
    endsubroutine input_array
!***********************************************************************
endmodule Messages
