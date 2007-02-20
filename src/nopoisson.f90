! $Id: nopoisson.f90,v 1.7 2007-02-20 17:46:22 dobler Exp $

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
    subroutine inverse_laplacian_fft(phi,kmax,c)
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
    endsubroutine inverse_laplacian_fft
!***********************************************************************
    subroutine inverse_laplacian_semispectral(phi,c)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  20-dec-2006/anders: dummy
!
      real, dimension (nx,ny,nz) :: phi
      real, optional             :: c
!
      if (NO_WARN) print*, phi
!
    endsubroutine inverse_laplacian_semispectral
!***********************************************************************

endmodule Poisson

