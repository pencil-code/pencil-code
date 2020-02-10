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
    subroutine inverse_laplacian(phi)
!
!  Solve the Poisson equation by Fourier transforming on a periodic grid.
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension(nx,ny,nz), intent(in) :: phi
!
      call keep_compiler_quiet(phi)
!
    endsubroutine inverse_laplacian
!***********************************************************************
    subroutine inverse_laplacian_fft_z(phi)
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension(nx,ny,nz), intent(in) :: phi
!
      call keep_compiler_quiet(phi)
!
    endsubroutine inverse_laplacian_fft_z
!***********************************************************************
    subroutine inverse_laplacian_z_2nd_neumann(f)
!
!  15-may-2006/anders+jeff: dummy
!
      real, dimension(:,:,:,:), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine inverse_laplacian_z_2nd_neumann
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
    subroutine read_poisson_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_poisson_init_pars
!***********************************************************************
    subroutine write_poisson_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_poisson_init_pars
!***********************************************************************
    subroutine read_poisson_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_poisson_run_pars
!***********************************************************************
    subroutine write_poisson_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_poisson_run_pars
!***********************************************************************
    subroutine get_acceleration(acceleration)
!
      real, dimension(nx,ny,nz,3), intent(out) :: acceleration           !should I (CAN I?) make this allocatable?
!
      call keep_compiler_quiet(acceleration)
!
    endsubroutine get_acceleration
!***********************************************************************
endmodule Poisson
