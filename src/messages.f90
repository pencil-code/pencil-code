! $Id: messages.f90,v 1.1 2005-07-01 02:56:08 mee Exp $

!  This module takes care of entropy (initial condition
!  and time advance)
module Messages
  use Mpicomm, only: stop_it, die_gracefully

  implicit none

  private

  public :: initialize_messages
  public :: information, warning, error
  public :: fatal_error, not_implemented
  public :: life_support_on, life_support_off
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
!
  integer, public, parameter :: iip_EVERYTHING  = 0
  integer, public, parameter :: iip_DEFAULT     = 0
!
  integer, parameter :: iwarning_ip     = 1000
  integer, parameter :: iinformation_ip = 1000
!
  integer :: warnings=0
  integer :: errors=0
  integer :: fatal_errors=0
!
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
!
        inquire(FILE="COLOR", EXIST=ltermcap_color)

    endsubroutine initialize_messages
!***********************************************************************
    subroutine not_implemented(location)
!
      character(len=*) :: location

      if (.not.llife_support) then
        errors=errors+1
!
        call terminal_highlight_error()
        write (*,'(A7)',ADVANCE='NO') "NOT IMPLEMENTED: "
        call terminal_defaultcolor()
        write (*,*) "Attempted to use a routine not capable of handling the"// &
                    "current parameters at "//trim(location)
!
        if (ldie_onfatalerror) call die_gracefully
      endif
    endsubroutine not_implemented
!***********************************************************************
    subroutine fatal_error(location,message)
!
      character(len=*) :: location
      character(len=*) :: message
!
      if (.not.llife_support) then
        errors=errors+1
!  
        call terminal_highlight_fatal_error()
        write (*,'(A13)',ADVANCE='NO') "FATAL ERROR: "
        call terminal_defaultcolor()
        write (*,*) trim(message)//" occured at "//trim(location)
!  
        if (ldie_onfatalerror) call die_gracefully
      endif
    endsubroutine fatal_error
!***********************************************************************
    subroutine error(location,message)
!
      character(len=*) :: location
      character(len=*) :: message

      if (.not.llife_support) then
        errors=errors+1
!
        call terminal_highlight_error()
        write (*,'(A7)',ADVANCE='NO') "ERROR: "
        call terminal_defaultcolor()
        write (*,*) trim(message)//" occured at "//trim(location)
!
        if (ldie_onerror) call die_gracefully
      endif
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
        call terminal_highlight_warning()
        write (*,'(A9)',ADVANCE='NO') "WARNING:"
        call terminal_defaultcolor()
        write (*,*) trim(message)//" occured at "//trim(location)
!
        if (ldie_onwarning) call die_gracefully
      endif
    endsubroutine warning
!***********************************************************************
    subroutine information(location,message,level)
!    
!  Print out colored warning.
!
!  30-jun-05/tony: coded
!
      use Cdata, only: ip
!
      character (len=*) :: message,location
      integer, optional :: level
      integer :: level_ = iinformation_ip
!
      if (present(level)) level_=level 
!
      if (ip<=level_) write (*,*) trim(location)//":"//trim(message)
!
    endsubroutine information
!***********************************************************************
    subroutine life_support_off
!    
!  Allow code to die on errors 
!
!  30-jun-05/tony: coded
!
      llife_support=.false.
      call warning('life_support_off','death on error restored')
!
    endsubroutine life_support_off
!***********************************************************************
    subroutine life_support_on
!    
!  Prevent the code from dying on errors 
!
!  30-jun-05/tony: coded
!
      llife_support=.true.
      call warning('life_support_on','death on error has been suspended')
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
!
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
endmodule Messages
