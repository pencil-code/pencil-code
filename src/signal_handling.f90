! $Id$

module Signal_handling

!  Used for signal handling in run.f90

  implicit none

  private

  public :: signal_prepare
!
  contains
!***********************************************************************
subroutine regexit()

! signal handling routine, touches STOP if root node, does nothing otherwise

!  12-nov-09/MR: coded

  use Cdata, only: lroot

  integer :: stat,system

  if (lroot) then
    print *,'End of the program run.x due to signal'
    stat=system('touch STOP')
  endif

endsubroutine regexit
!*****************************************************************************
subroutine signal_prepare()
!
! example signal catching for SIGINT and SIGUSR1
!
! declarations for signal handling
!
  integer, parameter :: SIGFPE=8, SIGINT=2, SIGHUP=1, SIGTERM=15, SIGUSR1=30
  integer, parameter :: USER=-1
  integer :: sigret,signal
!
  sigret = signal( SIGINT , regexit, USER )
  sigret = signal( SIGUSR1, regexit, USER )
!
endsubroutine signal_prepare
!*****************************************************************************
endmodule Signal_handling
