! $Id$
!
!  This module contains FFT wrapper subroutines.
!
module Fourier
!
  use Cdata
  use Cparam
  use Messages, only: fatal_error
  use Mpicomm, only: transp
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
  contains
!
!***********************************************************************
    subroutine fourier_transform(a_re,a_im,linv)
!
!  Subroutine to do fourier transform
!  The routine overwrites the input data
!
!  22-oct-02/axel+tarek: adapted from transform
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      if (present(linv)) then
        if (linv) then
          if (lroot) print*, 'fourier_transform: only implemented for '// &
              'forwards transform!'
          call fatal_error('fourier_transform','')
        endif
      endif
!
      if (lroot .and. ip<10) print*,'transform_i: doing three FFTs'
      call fft(a_re,a_im, nx*ny*nz, nx, nx,-1)
      call transp(a_re,'y')
      call transp(a_im,'y')
      call fft(a_re,a_im, nx*ny*nz, nx, nx,-1)
      call transp(a_re,'z')
      call transp(a_im,'z')
      call fft(a_re,a_im, nx*ny*nz, nx, nx,-1)
!
!  Normalize
!
      a_re=a_re/nwgrid
      a_im=a_im/nwgrid
      if (lroot .and. ip<10) print*,'transform_i: fft has finished'
!
    endsubroutine fourier_transform
!***********************************************************************
    subroutine fourier_transform_xz(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x- and z-directions
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_xz', &
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
!
    endsubroutine fourier_transform_xz
!***********************************************************************
    subroutine fourier_transform_xy(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x- and y-directions
!
      real, dimension(nx,ny) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_xy', &
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
!
    endsubroutine fourier_transform_xy
!***********************************************************************
    subroutine fourier_transform_shear_xy(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x- and y-directions
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_xy_shear', &
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
!
    endsubroutine fourier_transform_shear_xy
!***********************************************************************
    subroutine fourier_transform_xy_xy(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x- and y-directions
!
      real, dimension(nx,ny) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_xy_xy', &
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
!
    endsubroutine fourier_transform_xy_xy
!***********************************************************************
    subroutine fourier_transform_y_y(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform of a 1-D array under MPI.
!
      real, dimension(ny) :: a_re, a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_y_y', &
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
!
    endsubroutine fourier_transform_y_y
!***********************************************************************
    subroutine fourier_transform_xy_xy_other(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x- and y-directions
!
      real, dimension(nx,ny) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_xy_xy_other', &
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
!
    endsubroutine fourier_transform_xy_xy_other
!***********************************************************************
    subroutine fourier_transform_x(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x-direction.
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_x', &
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
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
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
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
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
!
    endsubroutine fourier_transform_shear
!***********************************************************************
    subroutine fourier_transform_other_1(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform on a 1-D array of arbitrary size.
!
      real, dimension(:) :: a_re,a_im
      logical, optional :: linv
!
      call fatal_error('fourier_transform_other_1', &
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
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
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(a_im)
      call keep_compiler_quiet(linv)
!
    endsubroutine fourier_transform_other_2
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
          'this sub is not available in fourier_fft.f90!')
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
      call fatal_error('fourier_shif_y', &
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a_re)
      call keep_compiler_quiet(shift_y)
!
    endsubroutine fourier_shift_y
!***********************************************************************
    subroutine fourier_transform_real_1(a,na,ifirst_fft,wsavex_temp,linv)
!
!   3-jul-08/dhruba: dummy routine
!
      real, dimension(na) :: a
      integer, intent(in) :: na,ifirst_fft
      logical, optional :: linv
      real, dimension(2*na+15),optional :: wsavex_temp
!
      call fatal_error('fourier_transform_real', &
          'this sub is not available in fourier_fft.f90!')
!
      call keep_compiler_quiet(a)
      call keep_compiler_quiet(na)
      call keep_compiler_quiet(ifirst_fft)
      call keep_compiler_quiet(linv)
      call keep_compiler_quiet(wsavex_temp)
!
    endsubroutine fourier_transform_real_1
!***********************************************************************
    subroutine fourier_transform_xy_xy_parallel(in_re,out_re,factor)
!
! Subroutine to do multi functional Fourier transform of a 2-D
! array under MPI in parallel for ghost cells.
!
      real, dimension(:,:,:) :: in_re,factor
      real, dimension(:,:,:,:) :: out_re
!
      call fatal_error('fourier_transform_xy_xy_parallel', &
          'this sub is not available in nofourier.f90!')
!
      call keep_compiler_quiet(in_re)
      call keep_compiler_quiet(out_re)
      call keep_compiler_quiet(factor)
!
    endsubroutine fourier_transform_xy_xy_parallel
!***********************************************************************
endmodule Fourier
