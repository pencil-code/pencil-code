! $Id: nogeneral_fft.f90,v 1.2 2006-07-28 13:29:41 ajohan Exp $
!
!  This module contains FFT wrapper subroutines.
!
module General_FFT

  use Cdata
  use Cparam
  use Messages
  use Mpicomm, only: transp

  implicit none

  include 'general_fft.h'

  contains

!***********************************************************************
    subroutine transform_i(a_re,a_im)
!
!  Subroutine to do fourier transform
!  The routine overwrites the input data
!
!  28-jul-06/anders: dummy
!
      real, dimension(nx,ny,nz) :: a_re,a_im
!
      call fatal_error('transform_i','sub not available in nogeneral_fft.f90!')
!
      if (NO_WARN) print*, a_re, a_im
!
    endsubroutine transform_i
!***********************************************************************
    subroutine transform_fftpack(a_re,a_im,direction)
!
!  Subroutine to do Fourier transform
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  28-jul-06/anders: dummy
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      integer, optional :: direction
!
      call fatal_error('transform_fftpack', &
          'sub not available in nogeneral_fft.f90!')
!
      if (NO_WARN) print*, a_re, a_im, direction
!
    endsubroutine transform_fftpack
!***********************************************************************
    subroutine transform_fftpack_2d(a_re,a_im,direction)
!
!  Subroutine to do Fourier transform
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  28-jul-06/anders: dummy
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      integer, optional :: direction
!
      call fatal_error('transform_fftpack_2d', &
          'sub not available in nogeneral_fft.f90!')
!
      if (NO_WARN) print*, a_re, a_im, direction
!
    endsubroutine transform_fftpack_2d
!***********************************************************************
    subroutine transform_fftpack_1d(a_re,a_im)
!
!  Subroutine to do Fourier transform
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  28-jul-06/anders: dummy
!
      real, dimension(nx,ny,nz) :: a_re,a_im
!
      call fatal_error('transform_fftpack_1d', &
          'sub not available in nogeneral_fft.f90!')
!
      if (NO_WARN) print*, a_re, a_im      
!
    endsubroutine transform_fftpack_1d
!***********************************************************************
    subroutine transform_fftpack_shear(a_re,a_im,direction)
!
!  Subroutine to do Fourier transform in shearing coordinates.
!  The routine overwrites the input data
!
!  28-jul-06/anders: dummy
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      integer :: direction
!
      call fatal_error('transform_fftpack_shear', &
          'sub not available in nogeneral_fft.f90!')
!
      if (NO_WARN) print*, a_re, a_im, direction
!
    endsubroutine transform_fftpack_shear
!***********************************************************************
endmodule General_FFT
