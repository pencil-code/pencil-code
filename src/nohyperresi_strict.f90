! $Id$
!
!  This module applies a sixth order hyperresistivity to the induction
!  equation (following Brandenburg & Sarson 2002). This hyperresistivity
!  ensures that the energy dissipation rate is positive define everywhere.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lhyperresistivity_strict=.false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Hyperresi_strict
!
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'hyperresi_strict.h'
!
  contains
!***********************************************************************
    subroutine register_hyperresi_strict()
!
!  Set up indices for hyperresistivity auxiliary slots.
!
!  23-aug-07/anders: adapted from register_hypervisc_strict
!
    endsubroutine register_hyperresi_strict
!***********************************************************************
    subroutine hyperresistivity_strict(f)
!
!  Apply sixth order hyperresistivity with positive definite heating rate
!  (see Brandenburg & Sarson 2002).
!
!  23-aug-07/anders: adapted from hyperviscosity_strict
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine hyperresistivity_strict
!***********************************************************************
endmodule Hyperresi_strict
