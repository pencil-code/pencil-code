! $Id$
!
!  This module contains FFT wrapper subroutines.
!
module Fourier
!
  use Cdata
  use Messages, only: fatal_error
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'fourier.h'
!
  interface fourier_transform_other
    module procedure fourier_transform_other_1
    module procedure fourier_transform_other_2
  endinterface
!
  interface fft_x_parallel
    module procedure fft_x_parallel_1D
    module procedure fft_x_parallel_2D
    module procedure fft_x_parallel_3D
    module procedure fft_x_parallel_4D
  endinterface
!
  interface fft_y_parallel
    module procedure fft_y_parallel_1D
    module procedure fft_y_parallel_2D
    module procedure fft_y_parallel_3D
    module procedure fft_y_parallel_4D
  endinterface
!
  interface fft_z_parallel
    module procedure fft_z_parallel_1D
    module procedure fft_z_parallel_2D
    module procedure fft_z_parallel_3D
    module procedure fft_z_parallel_4D
  endinterface
!
  interface fft_xy_parallel
    module procedure fft_xy_parallel_2D
    module procedure fft_xy_parallel_3D
    module procedure fft_xy_parallel_4D
  endinterface
!
  interface fft_xyz_parallel
    module procedure fft_xyz_parallel_3D
    module procedure fft_xyz_parallel_4D
  endinterface
!
  contains
!***********************************************************************
    subroutine fourier_transform(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform.
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform
!***********************************************************************
    subroutine fourier_transform_xy(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x- and y-directions.
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_xy', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform_xy
!***********************************************************************
    subroutine fourier_transform_xz(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x- and z-directions.
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_xz', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform_xz
!***********************************************************************
    subroutine fourier_transform_x(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x-direction.
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_x', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform_x
!***********************************************************************
    subroutine fourier_transform_y(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x-direction.
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_y', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform_y
!***********************************************************************
    subroutine fourier_transform_shear(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in shearing coordinates.
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_shear', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform_shear
!***********************************************************************
    subroutine fourier_transform_shear_xy(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in shearing coordinates in x- and
!  y-directions.
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_shear_xy', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform_shear_xy
!***********************************************************************
    subroutine fourier_transform_other_1(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform on a 1-D array of arbitrary size.
!
      real, dimension(:) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_other_1', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform_other_1
!***********************************************************************
    subroutine fourier_transform_other_2(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform of a 2-D array of arbitrary size.
!
      real, dimension(:,:) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_other_2', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform_other_2
!***********************************************************************
    subroutine fourier_transform_xy_xy(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do Fourier transform of a 2-D array of arbitrary size.
!
      real, dimension(nx,ny), intent(inout) :: a_re,a_im
      logical, optional, intent(in) :: linv,lneed_im
!
      call fatal_error('fourier_transform_xy_xy', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fourier_transform_xy_xy
!***********************************************************************
    subroutine fourier_transform_xy_xy_other(a_re,a_im,linv)
!
! Subroutine to do Fourier transform of a 2-D array of arbitrary size.
!
      real, dimension(:,:) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_xy_xy', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform_xy_xy_other
!***********************************************************************
    subroutine fft_x_parallel_1D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 1D data in the x-direction.
!
      real, dimension (nx), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_x_parallel_1D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_x_parallel_1D
!***********************************************************************
    subroutine fft_x_parallel_2D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 2D data in the x-direction.
!
      real, dimension (nx,ny), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_x_parallel_2D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_x_parallel_2D
!***********************************************************************
    subroutine fft_x_parallel_3D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 3D data in the x-direction.
!
      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_x_parallel_3D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_x_parallel_3D
!***********************************************************************
    subroutine fft_x_parallel_4D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 4D data in the x-direction.
!
      real, dimension (:,:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_x_parallel_4D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_x_parallel_4D
!***********************************************************************
    subroutine fft_y_parallel_1D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 1D data in the y-direction.
!
      real, dimension (ny), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_y_parallel_1D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_y_parallel_1D
!***********************************************************************
    subroutine fft_y_parallel_2D(a_re,a_im,linv,lneed_im,shift_y)
!
!  Subroutine to do FFT of distributed 2D data in the y-direction.
!
      real, dimension (nx,ny), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
      real, dimension(nxgrid), optional :: shift_y
!
      call fatal_error('fft_y_parallel_2D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
      call keep_compiler_quiet(present(shift_y))
!
    endsubroutine fft_y_parallel_2D
!***********************************************************************
    subroutine fft_y_parallel_3D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 3D data in the y-direction.
!
      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_y_parallel_3D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_y_parallel_3D
!***********************************************************************
    subroutine fft_y_parallel_4D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 4D data in the y-direction.
!
      real, dimension (:,:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_y_parallel_4D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_y_parallel_4D
!***********************************************************************
    subroutine fft_z_parallel_1D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 1D data in the z-direction.
!
      real, dimension (nz), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_z_parallel_1D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_z_parallel_1D
!***********************************************************************
    subroutine fft_z_parallel_2D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 2D data in the z-direction.
!
      real, dimension (:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_z_parallel_2D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_z_parallel_2D
!***********************************************************************
    subroutine fft_z_parallel_3D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 3D data in the z-direction.
!
      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_z_parallel_3D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_z_parallel_3D
!***********************************************************************
    subroutine fft_z_parallel_4D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 4D data in the z-direction.
!
      real, dimension (:,:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_z_parallel_4D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_z_parallel_4D
!***********************************************************************
    subroutine fft_xy_parallel_2D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 2D data in the x- and y-direction.
!
      real, dimension (nx,ny), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_xy_parallel_2D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_xy_parallel_2D
!***********************************************************************
    subroutine fft_xy_parallel_3D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 3D data in the x- and y-direction.
!
      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_xy_parallel_3D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_xy_parallel_3D
!***********************************************************************
    subroutine fft_xy_parallel_4D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 4D data in the x- and y-direction.
!
      real, dimension (:,:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_xy_parallel_4D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_xy_parallel_4D
!***********************************************************************
    subroutine fft_xyz_parallel_3D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 3D data in the x-, y- and z-direction.
!
      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_xyz_parallel_3D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_xyz_parallel_3D
!***********************************************************************
    subroutine fft_xyz_parallel_4D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 4D data in the x-, y- and z-direction.
!
      real, dimension (:,:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      call fatal_error('fft_xyz_parallel_4D', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(lneed_im))
!
    endsubroutine fft_xyz_parallel_4D
!***********************************************************************
    subroutine vect_pot_extrapol_z_parallel(in,out,factor)
!
!  Subroutine to do a z-extrapolation of a vector potential using
!  'factor' as a multiplication factor to the Fourier coefficients.
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
      real, dimension(:,:,:), intent(in) :: factor
!
      call fatal_error('vect_pot_extrapol_z_parallel', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(in)
      call keep_compiler_quiet(out)
      call keep_compiler_quiet(factor)
!
    endsubroutine vect_pot_extrapol_z_parallel
!***********************************************************************
    subroutine field_extrapol_z_parallel(in,out,factor)
!
!  Subroutine to do a z-extrapolation of a fields z-component using
!  'factor' as a multiplication factor to the Fourier coefficients.
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
      real, dimension(:,:,:), intent(in) :: factor
!
      call fatal_error('field_extrapol_z_parallel', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(in)
      call keep_compiler_quiet(out)
      call keep_compiler_quiet(factor)
!
    endsubroutine field_extrapol_z_parallel
!***********************************************************************
    subroutine fourier_transform_y_y(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform of a 1-D array under MPI.
!
      real, dimension(ny) :: a_re, a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_y_y', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(present(linv))
!
    endsubroutine fourier_transform_y_y
!***********************************************************************
    subroutine fourier_shift_yz_y(a_re,shift_y)
!
!  Performs a periodic shift in the y-direction of an entire y-z plane by
!  the amount shift_y.
!
!  02-oct-07/anders: dummy
!
      real, dimension (ny,nz) :: a_re
      real :: shift_y
!
      call fatal_error('fourier_shift_yz_y', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(shift_y)
!
    endsubroutine fourier_shift_yz_y
!***********************************************************************
    subroutine fourier_shift_y(a_re,shift_y)
!
!  Performs a periodic shift in the y-direction by the amount shift_y.
!
!  04-oct-07/anders: dummy
!
      real, dimension (nx,ny,nz) :: a_re
      real, dimension (nx) :: shift_y
!
      call fatal_error('fourier_transform_y', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(shift_y)
!
    endsubroutine fourier_shift_y
!***********************************************************************
    subroutine fourier_transform_real_1(a,na,ifirst_fft,wsavex_temp,linv)
!
!   1-jul-08/axel: dummy routine
!
      integer, intent(in) :: na
      real, dimension(na) :: a
      integer, intent(in) :: ifirst_fft
      logical, optional :: linv
      real, dimension(2*na+15), optional :: wsavex_temp
!
      call fatal_error('fourier_transform_real_1', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(a)
      call keep_compiler_quiet(na)
      call keep_compiler_quiet(ifirst_fft)
      call keep_compiler_quiet(present(linv))
      call keep_compiler_quiet(present(wsavex_temp))
!
    endsubroutine fourier_transform_real_1
!***********************************************************************
endmodule Fourier
