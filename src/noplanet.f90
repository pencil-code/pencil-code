! $Id: noplanet.f90,v 1.39 2007-02-05 22:09:30 wlyra Exp $
!
!  Dummy module
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************

module Planet
  use Cparam
  use Cdata
  use Messages

  implicit none
!
  include 'planet.h'

contains
!***********************************************************************
    subroutine pencil_criteria_planet()
!
!  All pencils that the Planet module depends on are specified here.
!
!  16-nov-06/tony: coded
!
    endsubroutine pencil_criteria_planet
!***********************************************************************
    subroutine pencil_interdep_planet(lpencil_in)
!
!  Interdependency among pencils from the Planet module is specified here.
!
!  16-nov-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if(NO_WARN) print*,lpencil_in(1)
!
    endsubroutine pencil_interdep_planet
!***********************************************************************
    subroutine calc_pencils_planet(f,p)
!
!  Calculate Planet pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  16-nov-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      if (NO_WARN) print*, f
!
    endsubroutine calc_pencils_planet
!*************************************
endmodule Planet
