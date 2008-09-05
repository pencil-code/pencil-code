! $Id$

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
!  regular spherical bessel function - dummy routine
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
!  irregular spherical bessel function - dummy routine
!
      real :: y,x
      integer :: l
!
      if (NO_WARN) print*,y,l,x
!
    endsubroutine sp_bessely_l_
!***********************************************************************
    subroutine sp_bessel_jnu_(y,nu,x)
!
!  regular cylindrical bessel function - dummy routine
!
      real :: y,x,nu
!
      if (NO_WARN) print*,y,nu,x
!
    endsubroutine sp_bessel_jnu_
!***********************************************************************
    subroutine sp_bessely_l_(y,nu,x)
!
!  irregular cylindrical bessel function - dummy routine
!
      real :: y,x,nu
!
      if (NO_WARN) print*,y,nu,x
!
    endsubroutine sp_bessel_ynu_
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

