! $Id: nogsl.f90,v 1.2 2007-09-14 05:47:25 ajohan Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Gsl

  use Cparam

  implicit none

  contains

!***********************************************************************
    subroutine sp_besselj_l(y,l,x)
!
!  dummy routine
!
      real :: y,x
      integer :: l
!
      if (NO_WARN) print*,y,l,x
!
    endsubroutine sp_besselj_l
!***********************************************************************
    subroutine sp_bessely_l_(y,l,x)
!
!  dummy routine
!
      real :: y,x
      integer :: l
!
      if (NO_WARN) print*,y,l,x
!
    endsubroutine sp_bessely_l_
!***********************************************************************
    subroutine sp_harm_real_(y,l,m,theta,phi)
!
!  dummy routine
!
      real :: y,theta,phi
      integer :: l,m
!
      if (NO_WARN) print*,y,l,m,theta,phi
!
    endsubroutine sp_harm_real_
!***********************************************************************
    subroutine sp_harm_imag_(y,l,m,theta,phi)
!
!  dummy routine
!
      real :: y,theta,phi
      integer :: l,m
!
      if (NO_WARN) print*,y,l,m,theta,phi
!
    endsubroutine sp_harm_imag_
!***********************************************************************

endmodule Gsl

