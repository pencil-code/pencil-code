! $Id: nopoisson.f90,v 1.8 2007-02-20 17:50:30 dobler Exp $

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
    subroutine poisson_solver_fft(a1,kmax)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension (nx,ny,nz) :: a1
      real, optional :: kmax
!
      if (NO_WARN) print*, a1, kmax
!
    endsubroutine poisson_solver_fft
!***********************************************************************
    subroutine poisson_solver_fftxy_discretez(a1)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  20-dec-2006/anders: dummy
!
      real, dimension (nx,ny,nz) :: a1
!
      if (NO_WARN) print*, a1
!
    endsubroutine poisson_solver_fftxy_discretez
!***********************************************************************

endmodule Poisson

