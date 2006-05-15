! $Id: nopoisson.f90,v 1.1 2006-05-15 21:30:50 ajohan Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Poisson

  use Cdata
  use Cparam
  use Messages

  implicit none

  include 'poisson.h'

  contains

!***********************************************************************
    subroutine poisson_solver_fft(f,irhs,ipoisson,rhs_const)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: rhs_const
      integer :: irhs, ipoisson
!
    endsubroutine poisson_solver_fft
!***********************************************************************

endmodule Poisson

