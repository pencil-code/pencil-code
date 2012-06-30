! $Id$

!
!  This module applies a sixth order hyperviscosity to the equation
!  of motion (following Haugen & Brandenburg 2004).
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lhyperviscosity_strict=.false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Hypervisc_strict
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  contains
!***********************************************************************
    subroutine register_hypervisc_strict()
!
!  Set up indices for hyperviscosity auxiliary slots.
!
!  20-aug-07/anders: coded
!
    endsubroutine register_hypervisc_strict
!***********************************************************************
    subroutine hyperviscosity_strict(f)
!
!  Apply momentum-conserving, symmetric, sixth order hyperviscosity with
!  positive define heating rate (see Haugen & Brandenburg 2004).
!
!  20-aug-2007/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine hyperviscosity_strict
!***********************************************************************
endmodule Hypervisc_strict

