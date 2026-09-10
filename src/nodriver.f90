! $Id$
!
!  Dummy module for boundary driving from external files.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldriver = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED
!
!***************************************************************
module Driver
!
  use Quiet
!
  implicit none
!
  private
!
  public :: initialize_driver, finalize_driver
  public :: read_driver_run_pars, write_driver_run_pars
  public :: driver_apply
!
  contains
!***********************************************************************
    subroutine initialize_driver
!
! Dummy routine.
!
    endsubroutine initialize_driver
!***********************************************************************
    subroutine finalize_driver
!
! Dummy routine.
!
    endsubroutine finalize_driver
!***********************************************************************
    subroutine read_driver_run_pars(iomsg)
!
      character(LEN=*), intent(out) :: iomsg
!
      iomsg = ""
!
    endsubroutine read_driver_run_pars
!***********************************************************************
    subroutine write_driver_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_driver_run_pars
!***********************************************************************
    subroutine driver_apply(f, df)
!
! Dummy routine.
!
! 07-Sep-2026/PABourdin: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
!
      call keep_compiler_quiet(f,df)
!
    endsubroutine driver_apply
!***********************************************************************
endmodule Driver
