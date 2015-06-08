! $Id$


!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lsignal = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Signal_handling

!  Used for signal handling in run.f90

  implicit none

  private

  public :: signal_prepare
  public :: read_signal_init_pars
  public :: write_signal_init_pars
  public :: emergency_stop

!
  logical :: emergency_stop=.false.
  integer, dimension(2) :: sigval=-1  ! 2 is the max number of signal to catch
! input parameters
  namelist /signal_init_pars/ &
      sigval

  contains
!***********************************************************************
subroutine regexit()

! signal handling routine, touches STOP if root node, does nothing otherwise

!  12-nov-09/MR: coded

  use Cdata, only: lroot

  integer :: stat,system

  if (lroot) then
    print *,'End of the program run.x due to signal'
    !stat=system('touch STOP')
  endif
  emergency_stop=.true.

endsubroutine regexit
!*****************************************************************************
subroutine signal_prepare()
!
! example signal catching for SIGINT and SIGUSR1
!
! declarations for signal handling
!
!  integer, parameter :: SIGFPE=8, SIGINT=2, SIGHUP=1, SIGTERM=15, SIGUSR1=10
!  Signal numbers are arch dependent.
! Instead, should be declared by user in start.in, signal_init_pars section
!  integer, parameter :: USER=-1
  integer :: i,sigret,signal
!
  do i=1, 2
    if (sigval(i) /= -1) then
      sigret = signal( sigval(i) , regexit) !, USER ) ! Do not compile here (gfortran) with third param. Raphael.
    endif
  enddo
!
endsubroutine signal_prepare
!***********************************************************************
    subroutine read_signal_init_pars(iostat)
!
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include "parallel_unit.h"
!
      read(parallel_unit, NML=signal_init_pars, IOSTAT=iostat)
!
    endsubroutine read_signal_init_pars
!***********************************************************************
    subroutine write_signal_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=signal_init_pars)
!
    endsubroutine write_signal_init_pars
!***********************************************************************
endmodule Signal_handling
