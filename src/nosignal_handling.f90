! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lsignal = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Signal_handling
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!  Signal handling in run.f90
!
  implicit none
  private
  public :: signal_prepare
  public :: read_signal_init_pars
  public :: write_signal_init_pars
  public :: emergency_stop
!
  logical :: emergency_stop = .false.
  contains
!***********************************************************************
    subroutine signal_prepare()
!
!  dummy routine
!
    endsubroutine signal_prepare
!***********************************************************************
    subroutine read_signal_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_signal_init_pars
!***********************************************************************
    subroutine write_signal_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_signal_init_pars
!***********************************************************************
endmodule Signal_handling
