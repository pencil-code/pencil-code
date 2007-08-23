! $Id: nohypervisc_strict.f90,v 1.2 2007-08-23 02:43:58 brandenb Exp $

!
!  This module applies a sixth order hyperviscosity to the equation
!  of motion (following Haugen & Brandenburg 2004).
!
!  The rate of strain tensor
!    S^(3) = (-nab^2)^2*S
!  is a high order generalisation of the first order operator
!    2*S_ij = u_i,j + u_j,i - 2/3*delta_ij*div(u)
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lhyperviscosity_strict=.true.
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
! 
    endsubroutine register_hypervisc_strict
!***********************************************************************
    subroutine hyperviscosity_strict(f,j)
!
!  Apply momentum-conserving, symmetric, sixth order hyperviscosity with
!  positive define heating rate (see Haugen & Brandenburg 2004).
!
!  20-aug-2007/anders: coded
!
      use Fourier
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      if (NO_WARN) print*,f,j !(keep compiler quiet)
!
    endsubroutine hyperviscosity_strict
!***********************************************************************

endmodule Hypervisc_strict

