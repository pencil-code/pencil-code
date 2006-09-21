! $Id: noplanet.f90,v 1.36 2006-09-21 23:19:17 wlyra Exp $
!
!  Dummy module
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lplanet = .false.
!
!***************************************************************

module Planet
  use Cparam
  use Cdata
  use Messages
  implicit none
contains
!*************************************
  subroutine pencil_criteria_planet()
  endsubroutine pencil_criteria_planet
!*************************************
  subroutine runtime_phiavg(p)
    type (pencil_case) :: p
    if (NO_WARN) print*, p
  endsubroutine runtime_phiavg
!*************************************
endmodule Planet
