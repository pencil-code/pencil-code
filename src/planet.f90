! $Id: planet.f90,v 1.82 2007-02-05 22:09:29 wlyra Exp $
!
!  This modules contains the routines for accretion disc and planet
!  building simulations.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!
!***************************************************************
!
! This module takes care of (mostly) everything related to the
! planet module
!
module Planet
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
!
  include 'planet.h'
!
  contains
!
!***********************************************************************
    subroutine pencil_criteria_planet()
!
!  All pencils that the Planet module depends on are specified here.
!
!  06-nov-05/wlad: coded
!  16-nov-06/tony: pencilised coordinates and averages
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
    endsubroutine pencil_interdep_planet
!*******************************************************************
    subroutine calc_pencils_planet(f,p)
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
    endsubroutine calc_pencils_planet
!***************************************************************
  endmodule Planet

