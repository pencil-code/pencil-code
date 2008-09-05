! $Id$

!  Dummy routine for ideal gas

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED
! PENCILS PROVIDED
!
!***************************************************************

module SpecialFunctions

  use Cparam

  implicit none

  contains
  
!***********************************************************************
    subroutine sp_besselj_l(r,l,res)
!
!   6-aug-07/axel: dummy routine
! 
      real :: r,res
      integer :: l

      res=0.

    end subroutine sp_besselj_l
!***********************************************************************
    subroutine sp_bessely_l(res,l,r)
!
!   6-aug-07/axel: dummy routine
! 
      real :: r,res
      integer :: l

      res=0.

    end subroutine sp_bessely_l
!***********************************************************************
    subroutine sp_harm_real(plm,l,m,theta,phi)
!
!   6-aug-07/axel: dummy routine
! 
      real :: plm,theta,phi
      integer :: l,m

      plm=0.

    end subroutine sp_harm_real
!***********************************************************************
    subroutine sp_harm_imag(plm,l,m,theta,phi)
!
!   6-aug-07/axel: dummy routine
! 
      real :: plm,theta,phi
      integer :: l,m

      plm=0.

    end subroutine sp_harm_imag
!***********************************************************************
endmodule SpecialFunctions
