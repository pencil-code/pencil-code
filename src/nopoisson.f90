! $Id$

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
    subroutine initialize_poisson()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-oct-07/anders: dummy
!
    endsubroutine initialize_poisson
!***********************************************************************
    subroutine inverse_laplacian(phi)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension (nx,ny,nz) :: phi
!
      if (NO_WARN) print*, phi
!
    endsubroutine inverse_laplacian
!***********************************************************************
    subroutine read_poisson_init_pars(unit,iostat)
!
!  Read Poisson init parameters.
!
!  18-oct-2007/anders: dummy
!
      integer :: unit
      integer, optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_poisson_init_pars
!***********************************************************************
    subroutine write_poisson_init_pars(unit)
!
!  Write Poisson init parameters.
!
!  18-oct-2007/anders: dummy
!
      integer :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_poisson_init_pars
!***********************************************************************
    subroutine read_poisson_run_pars(unit,iostat)
!
!  Read Poisson run parameters.
!
!  18-oct-2007/anders: dummy
!
      integer :: unit
      integer, optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_Poisson_run_pars
!***********************************************************************
    subroutine write_poisson_run_pars(unit)
!
!  Write Poisson run parameters.
!
!  18-oct-2007/anders: dummy
!
      integer :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_poisson_run_pars
!***********************************************************************

endmodule Poisson

