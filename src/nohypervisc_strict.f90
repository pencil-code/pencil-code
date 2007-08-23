! $Id: nohypervisc_strict.f90,v 1.4 2007-08-23 21:05:41 wlyra Exp $

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

  use Cdata
  use Cparam
  use Fourier
  use Messages

  implicit none

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
      if (NO_WARN) print*, f !(keep compiler quiet)
!
    endsubroutine hyperviscosity_strict
!***********************************************************************

endmodule Hypervisc_strict

