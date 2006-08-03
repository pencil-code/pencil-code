! $Id: nopoisson.f90,v 1.4 2006-08-03 07:07:28 ajohan Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpoisson=.false.
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
    subroutine poisson_solver_fft(a1,lklimit,kmax)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension (nx,ny,nz) :: a1
      real :: kmax
      logical :: lklimit
!
      if (NO_WARN) print*, a1, lklimit, kmax
!
    endsubroutine poisson_solver_fft
!***********************************************************************

endmodule Poisson

