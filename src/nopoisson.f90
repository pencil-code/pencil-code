! $Id$
!
!  This module solves the Poisson equation.
!
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
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'poisson.h'
!
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
    subroutine inverse_laplacian(f,phi)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: phi
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(phi)
      call keep_compiler_quiet(f)
!
    endsubroutine inverse_laplacian
!***********************************************************************
    subroutine inverse_laplacian_semispectral(f,phi)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: phi
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(phi)
      call keep_compiler_quiet(f)
!
    endsubroutine inverse_laplacian_semispectral
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
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
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
      call keep_compiler_quiet(unit)
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
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_poisson_run_pars
!***********************************************************************
    subroutine write_poisson_run_pars(unit)
!
!  Write Poisson run parameters.
!
!  18-oct-2007/anders: dummy
!
      integer :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_poisson_run_pars
!***********************************************************************
endmodule Poisson
