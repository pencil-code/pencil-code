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
      integer :: wd1=0, dist, iostat
      logical, save :: lfirstcall=.true.
      logical :: writ
!
!  Write string to screen and to 'svnid.dat' file.
!
      writ = .true.
      if (lfirstcall) then
        open(1, file=trim(datadir)//'/svnid.dat', status='replace',IOSTAT=iostat)
        dist=0
        lfirstcall=.false.
      else
        open(1, file=trim(datadir)//'/svnid.dat', status='old', position='append',IOSTAT=iostat)
        dist=-1
      endif
      if (outlog(iostat,'open',trim(datadir)//'/svnid.dat',dist=dist)) writ=.false.
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
        if (writ) then
          write(1,fmt,IOSTAT=iostat) "SVN: ", &
              trim(filename), &
              revision(1:wv), &
              date(1:wd), &
              trim(author)
          if (outlog(iostat,'write filename, revision etc.')) return
        endif
      else                      ! not a SVN line; maybe `[No ID given]'
        wd1 = min(wd, len(svnid))
        write(*,fmt) "SVN: ", &
            '-------', &
            '', &
            '', &
            svnid(1:wd1)
        if (writ) then
          write(1,fmt,IOSTAT=iostat) "SVN: ", &
              '-------', &
              '', &
              '', &
              svnid(1:wd1)
          if (outlog(iostat,'write svnid')) return
        endif
      endif
      !write(*,'(A)') '123456789|123456789|123456789|123456789|123456789|12345'
      !write(*,'(A)') '         1         2         3         4         5'
!
      if (writ) then
        close(1,IOSTAT=iostat)
        if (outlog(iostat,'close')) continue
      endif
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
      integer :: mul_fac,iostat
      logical, save :: opened=.true.
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
              open(lun,file=trim(datadir)//'/timing.dat', status='replace',IOSTAT=iostat)
              if (outlog(iostat,'open',trim(datadir)//'/timing.dat')) opened=.false.
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
              if (.not.opened) then
                open(lun,file=trim(datadir)//'/timing.dat',position='append',IOSTAT=iostat)
                opened = .not.(outlog(iostat,'open',trim(datadir)//'/timing.dat'))
              endif
              if (opened) then
                write(lun,*,IOSTAT=iostat) time,mul_fac,trim(location)//": "//trim(message)
                if (outlog(iostat,'write time,mul_fac etc.')) opened=.false.
              endif
            endif
          endif
!
!  finalize
!
          if (present(instruct)) then
            if (trim(instruct)=='finalize' .and. opened) then
              close(lun,IOSTAT=iostat)
              if (outlog(iostat,'close')) continue
            endif
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
  logical function outlog(code,mode,file,dist,msg)
!
!  Creates log entries for I/O errors in ioerrors.log.
!  Notifies user via e-mail if address mailaddress is given.
!  stops program if lstop_on_ioerror is set.
!  reverts incompletely written files to a defined state, in particular
!  distributed files in data/procN/ -> all sub-files in a coherent state
!
!  code(IN): errorcode from IOSTAT
!  mode(IN): describes failed action, starts with 'open', 'read', 'write' or 'close'
!            for 'read' and 'write': should contain the name of the relevant variable(s) 
!  file(IN): file at which operation failed, 
!            if omitted assumed to be the one saved in curfile
!            usually set by the call with mode='open'
!  dist(IN): indicator for distributed files (>0) and need of synchronization of file states across nodes
!            or simple backskipping (<0);|dist| = logical unit number
!            only considered in calls with mode='open'
!  msg(IN) : additional message text
!
!  return value: flag for 'I/O error has occurred'. If so execution should jump immediately after the 'close'
!                statement ending the present group of I/O operations as outlog closes (tries to close) the file.
!                It is in the responsibility of the programmer that by this jump no relevant statements are missed.
!                          
!  3-nov-11/MR: coded;
! 16-nov-11/MR: modified; experimental version which always stops program on I/O error
! 13-Dec-2011/Bourdin.KIS: added EOF sensing, which is not an error.
!
    use Syscalls, only: system_cmd
    use General, only: itoa,date_time_string,safe_character_append,safe_character_prepend,backskip
    use Mpicomm, only: report_clean_output
!
    integer,                     intent(IN) :: code
    character (LEN=*),           intent(IN) :: mode
    character (LEN=*), optional, intent(IN) :: file,msg
    integer,           optional, intent(IN) :: dist
!
    character (LEN=fnlen), save :: curfile=''
    ! default: file not distributed, no backskipping
    integer, save :: curdist=0, curback=0
!
    integer :: unit=90, iostat, ind
    character (LEN=intlen) :: date, codestr
    character (LEN=fnlen), dimension(2) :: strarr
    character (LEN=fnlen) :: filename, message, submsg
    logical :: lopen, lclose, lread, lwrite, lsync
!
    outlog = .false.
!
    filename = ""
    submsg = ""
    if (present (file)) filename = ' "'//trim (file)//'"'
    if (present (msg)) message = ': '//trim (msg)
!
    ! Set the following expression to .false. to activate the experimental code
    if (.true.) then
      if (code < 0) then
        outlog = .true.
        call warning(mode, 'reached End-Of-File'//trim (filename)//trim (message))
      else if (code > 0) then
        outlog = .true.
        call fatal_error(msg, 'I/O error'//trim (filename)//trim (message), .true.)
      endif
      return
    endif
!
    ! EXPERIMENTAL CODE:
    lopen  = mode(1:4)=='open'
    lread  = mode(1:5)=='read '
    lwrite = mode(1:6)=='write '
    lclose = mode(1:5)=='close'
!
    if (present(file)) curfile = file
!
    if (lopen) then
!
      if (present(dist)) then
        curdist = dist
      else
        curdist = 0
      endif
!
      curback = 0                                      ! counter for successfully read records set back  
!
    endif
!
    if (lwrite.and.code==0) curback = curback+1        ! number of succesfully written records after open
!
    if ( lroot ) then
      errormsg = ''; submsg = 'File'
    endif
!
    lsync = .false.
!
    if (curdist/=0) then                                       ! backskipping enabled 
!
      if ( ncpus==1 .and. curdist>0 ) curdist = -curdist
!
      if ( ncpus>1 .and. curdist>0 .and. (lwrite.or.lclose) ) then      
                                                               ! write/close on distributed file failed (somewhere)
        lsync = report_clean_output(code/=0, errormsg) 
        if (lsync) then                                        ! synchronization necessary as at least 
                                                               ! one processor failed in writing
          if (lserial_io) then                                 ! no backskipping, needs to be checked
            submsg = ' not synchronized (lserial_io=T)!'
            lsync = .false.                                    ! should file be closed??? Is return value then correct?
                                                               ! better to stop immediately?
          else
!
            if (lclose.and.code==0) open(curdist,file=curfile,position='append') ! re-open successfully written and closed files
!
            if ( backskip(curdist,curback) ) then                                ! try to set back file pointer by curback records
              if (lroot) submsg = trim(submsg)//' not synchronized!' 
            else
              if (lroot) submsg = trim(submsg)//' synchronized.'
            endif  
!
            close(curdist,IOSTAT=iostat)               ! try to close file
            if (iostat/=0) call safe_character_append(submsg,'. File not closed!')
!
          endif
!
        endif
!
      else if ( curdist<0 .and. code/=0 ) then         ! undistributed file, operation failed
!
        lsync = .false.
        if (lwrite.or.lclose) then 
!
          if ( backskip(abs(curdist),curback) ) submsg = trim(submsg)//' pointer not set back!'
!
          close(abs(curdist),IOSTAT=iostat)            ! try to close file
          if (iostat/=0) call safe_character_append(submsg,'. File not closed!')
!
        endif
!
      endif
    endif
!
    if ( lroot ) then
!
      if ( code==0 .and. .not.lsync ) return
!
      call safe_character_prepend( errormsg, 'ERROR ' )
!
      if ( lopen.or.lread.or.lwrite.or.lclose ) then
!
        call safe_character_append(errormsg,' when '//mode(1:4)//'ing ')
        if ( lread .or. lwrite ) then
          call safe_character_append(errormsg,trim(mode(6:)))
          if (lread) then
            call safe_character_append(errormsg,' from')
          else
            call safe_character_append(errormsg,' to')
          endif
        endif
        call safe_character_append(errormsg,' file "')
!
      else
        call safe_character_append(errormsg,': file "')
      endif
!
      filename = curfile
      if ( ncpus>1 .and. curdist>=0 ) then         ! for synchronized file replace 'procN' by 'proc*'
        ind = index(curfile,'proc')+4
        filename(ind:ind) = '*'
      endif
!
      codestr = itoa(code)
      call safe_character_append(errormsg,trim(filename)//'". Code: '//trim(codestr))
      if ( submsg/='File' ) call safe_character_append(errormsg,'. '//trim(submsg))
      if ( present(msg) ) call safe_character_append(errormsg,'. '//trim(msg))
!    
! scan of ioerrors.log to avoid multiple entries for the same file with same error code.
! When user eliminates cause of error, (s)he should also remove the corresponding line(s) in ioerrors.log.
!
      strarr(1) = filename
      strarr(2) = codestr
!
      if ( .not.scanfile('ioerrors.log',2,strarr,'all') ) then
!
        open(unit,file='ioerrors.log',position='append',iostat=IOSTAT)
!
        if (iostat==0) then
          call date_time_string(date)
          write(unit,'(a)',IOSTAT=iostat) date//' '//trim(errormsg)
          close(unit,IOSTAT=iostat)
        endif
!
        if (iostat/=0) write(*,'(a)',iostat=IOSTAT) date//' '//trim(errormsg)
!
        if ( index(mailaddress,'@') > 0 ) &        ! send mail to user if address could be valid
          call system_cmd( &
               'echo '//trim(errormsg)//'| mail -s PencilCode Message '//trim(mailaddress) )
      endif
    endif
!
    if (code/=0.or.lsync) then
      outlog = .true.
      if (lstop_on_ioerror) call stop_it('due to '//trim(errormsg))     ! stop on error if requested by user
    endif
!
  end function outlog
!***********************************************************************
  logical function scanfile(file,nstr,strings,mode)
!
!  Scans a file for a line in which one or all strings in list strings occur.
!  Returns on first hit.
!
!  3-nov-11/MR: coded
!  15-Feb-2012/Bourdin: removed deprecated features
!
    character (LEN=*),                           intent(IN) :: file
    integer,                                     intent(IN) :: nstr
    character (LEN=*), dimension(nstr),          intent(IN) :: strings
    character (LEN=3),                 optional, intent(IN) :: mode
!
    character (LEN=3) :: model
    character (LEN=fnlen) :: line
    integer :: lun=90,i,count,io_err
   
    if ( .not.present(mode) ) then
      model='any'
    else
      model=mode
    endif
!
    scanfile = .false.
!
    open(lun,file=file,IOSTAT=io_err)
    if (io_err /= 0) return
!
    do
      read(lun,'(a)',IOSTAT=io_err) line
      if (io_err /= 0) exit
!
      count=0
      do i=1,nstr
        if (index(line,trim(strings(i))) /= 0) then
          if (mode=='any') then
            scanfile = .true.
            exit
          endif
          count = count+1
        endif
      enddo
!
      if (count == nstr) scanfile = .true.
      if (scanfile) exit 
!
    enddo
!
    close(lun)
!
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
      if (outlog(iostat,'close',file)) continue
!
    endsubroutine input_array
!***********************************************************************
endmodule Messages
