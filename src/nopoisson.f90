! $Id: nopoisson.f90,v 1.10 2007-08-16 22:06:45 dobler Exp $

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
    subroutine inverse_laplacian(phi,kmax,c)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension (nx,ny,nz) :: phi
      real, optional             :: kmax,c
!
      if (NO_WARN) print*, phi, kmax
!
    endsubroutine inverse_laplacian
!***********************************************************************

endmodule Poisson

