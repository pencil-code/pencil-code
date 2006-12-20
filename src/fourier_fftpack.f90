! $Id: fourier_fftpack.f90,v 1.12 2006-12-20 13:10:50 ajohan Exp $
!
!  This module contains FFT wrapper subroutines.
!
module Fourier

  use Cdata
  use Cparam
  use Messages
  use Mpicomm, only: transp

  implicit none

  include 'fourier.h'

  interface fourier_transform_other
    module procedure fourier_transform_other_1
    module procedure fourier_transform_other_2
  endinterface

  contains

!***********************************************************************
    subroutine fourier_transform(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in 3-D.
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  27-oct-02/axel: adapted from transform_i, for fftpack
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension(nx) :: ax
      real, dimension(4*nx+15) :: wsavex
      integer :: l,m,n
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif
!
!  Need to initialize cfft only once, because we require nxgrid=nygrid=nzgrid
!
      call cffti(nx,wsavex)
!
      if (lforward) then
        if (lroot .and. ip<10) print*, 'fourier_transform: doing FFTpack in x'
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nx,ax,wsavex)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
!
!  Transform y-direction.
!
        if (nygrid/=1) then
          if (nygrid/=nxgrid) then
            if (lroot) print*, 'fourier_transform: must have nygrid=nxgrid!'
            call fatal_error('fourier_transform','')
          endif
          call transp(a_re,'y')
          call transp(a_im,'y')
!
!  The length of the array in the y-direction is nx.
!
          if (lroot .and. ip<10) print*, 'fourier_transform: doing FFTpack in y'
          do n=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftf(nx,ax,wsavex)
            a_re(:,l,n)=real(ax)
            a_im(:,l,n)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform z-direction.
!
        if (nzgrid/=1) then
          if (nzgrid/=nxgrid) then
            if (lroot) print*, 'fourier_transform: must have nzgrid=nxgrid!'
            call fatal_error('fourier_transform','')
          endif
          call transp(a_re,'z')
          call transp(a_im,'z')
!
!  The length of the array in the z-direction is also nx.
!
          if (lroot .and. ip<10) print*, 'fourier_transform: doing FFTpack in z'
          do m=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,m),a_im(:,l,m))
            call cfftf(nx,ax,wsavex)
            a_re(:,l,m)=real(ax)
            a_im(:,l,m)=aimag(ax)
          enddo; enddo
        endif
      else
!
!  Transform from k-space to real space. Works on arrays that have been put in
!  (z,x,y)-order by the parallel Fourier transform.
!
        if (nzgrid/=1) then
!
!  Transform z-direction back.
!
          if (nzgrid/=nxgrid) then
            if (lroot) print*, 'fourier_transform: must have nzgrid=nxgrid!'
            call fatal_error('fourier_transform','')
          endif
!
          if (lroot .and. ip<10) print*, 'fourier_transform: doing FFTpack in z'
          do m=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,m),a_im(:,l,m))
            call cfftb(nx,ax,wsavex)
            a_re(:,l,m)=real(ax)
            a_im(:,l,m)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform y-direction back. Must transpose to go from (z,x,y) to (y,x,z).
!
        if (nygrid/=1) then
          if (nygrid/=nxgrid) then
            if (lroot) print*, 'fourier_transform: must have nygrid=nxgrid!'
            call fatal_error('fourier_transform','')
          endif
!
          if (nzgrid/=1) then
            call transp(a_re,'z')
            call transp(a_im,'z')
          endif
!
          if (lroot .and. ip<10) print*, 'fourier_transform: doing FFTpack in y'
          do n=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftb(nx,ax,wsavex)
            a_re(:,l,n)=real(ax)
            a_im(:,l,n)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform x-direction back. Transpose to go from (y,x,z) to (x,y,z).
!
        if (lroot .and. ip<10) print*, 'fourier_transform: doing FFTpack in x'
        if (nygrid==1) then
          call transp(a_re,'z')
          call transp(a_im,'z')
        else
          call transp(a_re,'y')
          call transp(a_im,'y')
        endif
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nx,ax,wsavex)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/nwgrid
        a_im=a_im/nwgrid
      endif
!
      if (lroot .and. ip<10) print*, 'fourier_transform: fft has finished'
!
    endsubroutine fourier_transform
!***********************************************************************
    subroutine fourier_transform_xy(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in x and y.
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  19-12t-06/anders: adapted from fourier_transform
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension(nx) :: ax
      real, dimension(4*nx+15) :: wsavex
      integer :: l,m,n
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif
!
!  Need to initialize cfft only once, because we require nxgrid=nygrid.
!
      call cffti(nx,wsavex)
!
      if (lforward) then
        if (lroot .and. ip<10) print*, 'fourier_transform_xy: doing FFTpack in x'
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nx,ax,wsavex)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
!
!  Transform y-direction.
!
        if (nygrid/=1) then
          if (nygrid/=nxgrid) then
            if (lroot) print*, 'fourier_transform_xy: must have nygrid=nxgrid!'
            call fatal_error('fourier_transform_xy','')
          endif
          call transp(a_re,'y')
          call transp(a_im,'y')
!
!  The length of the array in the y-direction is nx.
!
          if (lroot .and. ip<10) print*, 'fourier_transform_xy: doing FFTpack in y'
          do n=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftf(nx,ax,wsavex)
            a_re(:,l,n)=real(ax)
            a_im(:,l,n)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform y-direction back.
!
        if (nygrid/=1) then
          if (nygrid/=nxgrid) then
            if (lroot) print*, 'fourier_transform: must have nygrid=nxgrid!'
            call fatal_error('fourier_transform','')
          endif
!
          if (lroot .and. ip<10) print*, 'fourier_transform_xy: doing FFTpack in y'
          do n=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftb(nx,ax,wsavex)
            a_re(:,l,n)=real(ax)
            a_im(:,l,n)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform x-direction back. Transpose to go from (y,x,z) to (x,y,z).
!
        if (lroot .and. ip<10) print*, 'fourier_transform_xy: doing FFTpack in x'
        if (nygrid/=1) then
          call transp(a_re,'y')
          call transp(a_im,'y')
        endif
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nx,ax,wsavex)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/(nxgrid*nygrid)
        a_im=a_im/(nxgrid*nygrid)
      endif
!
      if (lroot .and. ip<10) print*, 'fourier_transform_xy: fft has finished'
!
    endsubroutine fourier_transform_xy
!***********************************************************************
    subroutine fourier_transform_xz(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x- and z-directions.
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  27-oct-02/axel: adapted from transform_fftpack
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension(nx) :: ax
      real, dimension(4*nx+15) :: wsavex
      integer :: l,m,n
!
      if (present(linv)) then
        if (linv) then
          if (lroot) print*, 'fourier_transform_xz: only implemented for '// &
              'forwards transform!'
          call fatal_error('fourier_transform_xz','')
        endif
      endif
!
!  check whether nxgrid=nygrid=nzgrid
!
      if (nxgrid/=nygrid .or. (nxgrid/=nzgrid .and. nzgrid/=1)) then
        if (lroot) &
            print*, 'fourier_transform_xz: must have nxgrid=nygrid=nzgrid!'
        call fatal_error('fourier_transform_xz','')
      endif
!
!  need to initialize cfft only once, because nxgrid=nygrid=nzgrid
!
      call cffti(nx,wsavex)
!
      if (lroot .and. ip<10) print*, 'fourier_transform_xz: doing FFTpack in x'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        call cfftf(nx,ax,wsavex)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
      call transp(a_re,'z')
      call transp(a_im,'z')
!
!  The length of the array in the z-direction is also nx
!
      if (lroot .and. ip<10) print*, 'fourier_transform_xz: doing FFTpack in z'
      do l=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,l),a_im(:,m,l))
        call cfftf(nx,ax,wsavex)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
!
!  Normalize
!
      a_re=a_re/nwgrid
      a_im=a_im/nwgrid
      if (lroot .and. ip<10) print*, 'fourier_transform_xz: fft has finished'
!
    endsubroutine fourier_transform_xz
!***********************************************************************
    subroutine fourier_transform_x(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in the x-direction.
!  WARNING: It is not cache efficient to Fourier transform in any other
!  direction, so better transpose the array and only transform in x.
!
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  06-feb-03/nils: adapted from transform_fftpack
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension(nx) :: ax
      real, dimension(4*nx+15) :: wsavex
      integer :: m,n
!
      if (present(linv)) then
        if (linv) then
          if (lroot) print*, 'fourier_transform_x: only implemented for '// &
              'forwards transform!'
          call fatal_error('fourier_transform_x','')
        endif
      endif
!
!  check whether nxgrid=nygrid=nzgrid
!
      if ( (nygrid/=1.and.nygrid/=nxgrid) .or. &
           (nzgrid/=1.and.nzgrid/=nxgrid) ) then
        if (lroot) &
            print*, 'fourier_transform_x: must have nxgrid=nygrid=nzgrid!'
        call fatal_error('fourier_transform_x','')
      endif
!
!  need to initialize cfft only once, because nxgrid=nygrid=nzgrid
!
      call cffti(nx,wsavex)
!
      if (lroot .and. ip<10) print*, 'fourier_transform_x: doing FFTpack in x'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        call cfftf(nx,ax,wsavex)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
!
!  Normalize
!
      a_re=a_re/nxgrid
      a_im=a_im/nxgrid
      if (lroot .and. ip<10) print*, 'fourier_transform_x: fft has finished'
!
    endsubroutine fourier_transform_x
!***********************************************************************
    subroutine fourier_transform_shear(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in shearing coordinates.
!  The routine overwrites the input data
!
!  The transform from real space to Fourier space goes as follows:
!
!  (x,y,z)
!     |
!  (y,x,z) - (ky,x,z) - (ky',x,z)
!                            |
!                       (x,ky',z) - (kx,ky',z)
!                                        |
!                                   (z,ky',kx) - (kz,ky',kx)
!
!  Vertical lines refer to a transposition operation, horizontal lines to any
!  other operation. Here ky' refers to a coordinate frame where y has been
!  transformed so that the x-direction is purely periodic (not shear periodic).
!  The transformation from Fourier space to real space takes place like sketched
!  in the diagram, only in the opposite direction (starting lower right).
!
!  For 2-D runs two cases of degeneracy have to be taken into account:
!  (- refers to a missing direction).
!
!  (x,y,-)
!     |
!  (y,x,-) - (ky,x,-) - (ky',x,-)
!                            |
!                       (x, ky',-) - (kx,ky',-)
!
!  (x,-,z) - (kx,-,z)
!                |
!            (z,-,kx) - (kz,-,kx)
!
!  25-may-06/anders: adapted from transform_fftpack
!
      use Cdata, only: pi, Lx, x0, x, deltay, ky_fft
!
      real, dimension (nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (nxgrid) :: ax
      complex, dimension (nxgrid) :: ay
      complex, dimension (nxgrid) :: az
      real, dimension (4*nxgrid+15) :: wsave
      real :: deltay_x
      integer :: l,m,n
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif
!
!  if nxgrid/=nygrid/=nzgrid, stop.
!
      if (nygrid/=nxgrid .and. nygrid /= 1) then
        print*, 'fourier_transform_shear: '// &
            'need to have nygrid=nxgrid if nygrid/=1.'
        call fatal_error('fourier_transform_shear','')
      endif
      if (nzgrid/=nxgrid .and. nzgrid /= 1) then
        print*,'fourier_transform_shear: '// &
            'need to have nzgrid=nxgrid if nzgrid/=1.'
        call fatal_error('fourier_transform_shear','')
      endif
!
!  Need to initialize cfft only once, because we require nxgrid=nygrid=nzgrid.
!
      call cffti(nxgrid,wsave)
!
      if (lforward) then
!
!  Transform y-direction. Must start with y, because x is not periodic (yet).
!
        if (nygrid/=1) then
          if (lroot.and.ip<10) &
              print*, 'fourier_transform_shear: doing FFTpack in y'
          call transp(a_re,'y')
          call transp(a_im,'y')
          do n=1,nz; do l=1,ny
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftf(nxgrid,ay,wsave)
!  Shift y-coordinate so that x-direction is periodic. This is best done in
!  k-space, by multiplying the complex amplitude by exp[i*ky*deltay(x)].
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            ay(2:nxgrid)=ay(2:nxgrid)*exp(cmplx(0.0, ky_fft(2:nxgrid)*deltay_x))
            a_re(:,l,n)=real(ay)
            a_im(:,l,n)=aimag(ay)
          enddo; enddo
        endif
!
!  Transform x-direction.
!
        if (lroot.and.ip<10) &
            print*, 'fourier_transform_shear: doing FFTpack in x'
        if (nygrid/=1) then
          call transp(a_re,'y')
          call transp(a_im,'y')
        endif
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nxgrid,ax,wsave)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
!
!  Transform z-direction.
!
        if (nzgrid/=1) then
          if (lroot.and.ip<10) &
              print*, 'fourier_transform_shear: doing FFTpack in z'
          call transp(a_re,'z')
          call transp(a_im,'z')
          do l=1,nz; do m=1,ny
            az=cmplx(a_re(:,m,l),a_im(:,m,l))
            call cfftf(nxgrid,az,wsave)
            a_re(:,m,l)=real(az)
            a_im(:,m,l)=aimag(az)
          enddo; enddo
        endif
      else
!
!  Transform z-direction back.
!
        if (nzgrid/=1) then
          if (lroot.and.ip<10) &
              print*, 'fourier_transform_shear: doing FFTpack in z'
          do l=1,nz; do m=1,ny
            az=cmplx(a_re(:,m,l),a_im(:,m,l))
            call cfftb(nxgrid,az,wsave)
            a_re(:,m,l)=real(az)
            a_im(:,m,l)=aimag(az)
          enddo; enddo
        endif
!
!  Transform x-direction back.
!
        if (lroot.and.ip<10) &
            print*, 'fourier_transform_shear: doing FFTpack in x'
        if (nzgrid/=1) then
          call transp(a_re,'z')
          call transp(a_im,'z')
        endif
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nxgrid,ax,wsave)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
!
!  Transform y-direction back. Would be more practical to end with the
!  x-direction, but the transformation from y' to y means that we must transform
!  the y-direction last.
!
        if (nygrid/=1) then
          call transp(a_re,'y')
          call transp(a_im,'y')
          if (lroot.and.ip<10) &
              print*, 'fourier_transform_shear: doing FFTpack in y'
          do n=1,nz; do l=1,ny
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
!  Shift y-coordinate back to regular frame (see above).
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            ay(2:nxgrid)=ay(2:nxgrid)*exp(cmplx(0.0,-ky_fft(2:nxgrid)*deltay_x))
            call cfftb(nxgrid,ay,wsave)
            a_re(:,l,n)=real(ay)
            a_im(:,l,n)=aimag(ay)
          enddo; enddo
          call transp(a_re,'y')  ! Deliver array back in (x,y,z) order.
          call transp(a_im,'y')
        endif
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/nwgrid
        a_im=a_im/nwgrid
      endif
!
    endsubroutine fourier_transform_shear
!***********************************************************************
    subroutine fourier_transform_shear_xy(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in shearing coordinates.
!  The routine overwrites the input data
!
!  The transform from real space to Fourier space goes as follows:
!
!  (x,y,z)
!     |
!  (y,x,z) - (ky,x,z) - (ky',x,z)
!                            |
!                       (x,ky',z) - (kx,ky',z)
!
!  Vertical lines refer to a transposition operation, horizontal lines to any
!  other operation. Here ky' refers to a coordinate frame where y has been
!  transformed so that the x-direction is purely periodic (not shear periodic).
!  The transformation from Fourier space to real space takes place like sketched
!  in the diagram, only in the opposite direction (starting lower right).
!
!  19-dec-06/anders: adapted from fourier_transform_shear
!
      use Cdata, only: pi, Lx, x0, x, deltay, ky_fft
!
      real, dimension (nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (nxgrid) :: ax
      complex, dimension (nxgrid) :: ay
      complex, dimension (nxgrid) :: az
      real, dimension (4*nxgrid+15) :: wsave
      real :: deltay_x
      integer :: l,m,n
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif
!
!  if nxgrid/=nygrid/=nzgrid, stop.
!
      if (nygrid/=nxgrid .and. nygrid /= 1) then
        print*, 'fourier_transform_shear_xy: '// &
            'need to have nygrid=nxgrid if nygrid/=1.'
        call fatal_error('fourier_transform_shear','')
      endif
!
!  Need to initialize cfft only once, because we require nxgrid=nygrid=nzgrid.
!
      call cffti(nxgrid,wsave)
!
      if (lforward) then
!
!  Transform y-direction. Must start with y, because x is not periodic (yet).
!
        if (nygrid/=1) then
          if (lroot.and.ip<10) &
              print*, 'fourier_transform_shear: doing FFTpack in y'
          call transp(a_re,'y')
          call transp(a_im,'y')
          do n=1,nz; do l=1,ny
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftf(nxgrid,ay,wsave)
!  Shift y-coordinate so that x-direction is periodic. This is best done in
!  k-space, by multiplying the complex amplitude by exp[i*ky*deltay(x)].
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            ay(2:nxgrid)=ay(2:nxgrid)*exp(cmplx(0.0, ky_fft(2:nxgrid)*deltay_x))
            a_re(:,l,n)=real(ay)
            a_im(:,l,n)=aimag(ay)
          enddo; enddo
        endif
!
!  Transform x-direction.
!
        if (lroot.and.ip<10) &
            print*, 'fourier_transform_shear: doing FFTpack in x'
        if (nygrid/=1) then
          call transp(a_re,'y')
          call transp(a_im,'y')
        endif
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nxgrid,ax,wsave)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
      else
!
!  Transform x-direction back.
!
        if (lroot.and.ip<10) &
            print*, 'fourier_transform_shear: doing FFTpack in x'
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nxgrid,ax,wsave)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
!
!  Transform y-direction back. Would be more practical to end with the
!  x-direction, but the transformation from y' to y means that we must transform
!  the y-direction last.
!
        if (nygrid/=1) then
          call transp(a_re,'y')
          call transp(a_im,'y')
          if (lroot.and.ip<10) &
              print*, 'fourier_transform_shear: doing FFTpack in y'
          do n=1,nz; do l=1,ny
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
!  Shift y-coordinate back to regular frame (see above).
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            ay(2:nxgrid)=ay(2:nxgrid)*exp(cmplx(0.0,-ky_fft(2:nxgrid)*deltay_x))
            call cfftb(nxgrid,ay,wsave)
            a_re(:,l,n)=real(ay)
            a_im(:,l,n)=aimag(ay)
          enddo; enddo
          call transp(a_re,'y')  ! Deliver array back in (x,y,z) order.
          call transp(a_im,'y')
        endif
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/(nxgrid*nygrid)
        a_im=a_im/(nxgrid*nygrid)
      endif
!
    endsubroutine fourier_transform_shear_xy
!***********************************************************************
    subroutine fourier_transform_other_1(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform on a 1-D array of arbitrary size.
!  The routine overwrites the input data.
!
!  28-jul-2006/anders: adapted from fourier_transform
!
      real, dimension(:) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension(size(a_re,1)) :: ax
      real, dimension(4*size(a_re,1)+15) :: wsavex
      integer :: nx_other
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif
!
      nx_other=size(a_re,1)
!
!  Initialize fftpack.
!
      call cffti(nx_other,wsavex)
!
!  Transform x-direction.
!
      if (lforward) then
        if (lroot .and. ip<10) &
            print*, 'fourier_transform_other_1: doing FFTpack in x'
        ax=cmplx(a_re,a_im)
        call cfftf(nx_other,ax,wsavex)
        a_re=real(ax)
        a_im=aimag(ax)
      else
!
!  Transform x-direction back.
!
        if (lroot .and. ip<10) &
            print*, 'fourier_transform_other_1: doing FFTpack in x'
        ax=cmplx(a_re,a_im)
        call cfftb(nx_other,ax,wsavex)
        a_re=real(ax)
        a_im=aimag(ax)
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/nx_other
        a_im=a_im/nx_other
      endif
!
      if (lroot .and. ip<10) &
          print*, 'fourier_transform_other_1: fft has finished'
!
    endsubroutine fourier_transform_other_1
!***********************************************************************
    subroutine fourier_transform_other_2(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform of a 2-D array of arbitrary size.
!  The routine overwrites the input data.
!
!  28-jul-2006/anders: adapted from fourier_transform_1
!
      real, dimension(:,:) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension(size(a_re,1)) :: ax
      complex, dimension(size(a_re,2)) :: ay
      real, dimension(4*size(a_re,1)+15) :: wsavex
      real, dimension(4*size(a_re,2)+15) :: wsavey
      integer :: l, m, nx_other, ny_other
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif
!
      nx_other=size(a_re,1); ny_other=size(a_re,2)
!
      if (lforward) then
!
!  Transform x-direction.
!
        call cffti(nx_other,wsavex)
!
        if (lroot .and. ip<10) &
            print*, 'fourier_transform_other_2: doing FFTpack in x'
        do m=1,ny_other
          ax=cmplx(a_re(:,m),a_im(:,m))
          call cfftf(nx_other,ax,wsavex)
          a_re(:,m)=real(ax)
          a_im(:,m)=aimag(ax)
        enddo
!
!  Transform y-direction.
!
        call cffti(ny_other,wsavey)
!
        if (lroot .and. ip<10) &
            print*, 'fourier_transform_other_2: doing FFTpack in y'
        do l=1,nx_other
          ay=cmplx(a_re(l,:),a_im(l,:))
          call cfftf(ny_other,ay,wsavey)
          a_re(l,:)=real(ay)
          a_im(l,:)=aimag(ay)
        enddo
      else
!
!  Transform x-direction back.
!
        call cffti(nx_other,wsavex)
!
        if (lroot .and. ip<10) &
            print*, 'fourier_transform_other_2: doing FFTpack in x'
        do m=1,ny_other
          ax=cmplx(a_re(:,m),a_im(:,m))
          call cfftb(nx_other,ax,wsavex)
          a_re(:,m)=real(ax)
          a_im(:,m)=aimag(ax)
        enddo
!
!  Transform y-direction back.
!
        call cffti(ny_other,wsavey)
!
        if (lroot .and. ip<10) &
            print*, 'fourier_transform_other_2: doing FFTpack in y'
        do l=1,nx_other
          ay=cmplx(a_re(l,:),a_im(l,:))
          call cfftb(ny_other,ay,wsavey)
          a_re(l,:)=real(ay)
          a_im(l,:)=aimag(ay)
        enddo
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/(nx_other*ny_other)
        a_im=a_im/(nx_other*ny_other)
      endif
!
      if (lroot .and. ip<10) &
          print*, 'fourier_transform_other_2: fft has finished'
!
    endsubroutine fourier_transform_other_2
!***********************************************************************
    subroutine fourier_transform_xy_parallel(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform of a 2-D array under MPI.
!  nxgrid is restricted to be an integer multiple of nygrid.
!  The routine overwrites the input data.
!
!   6-oct-2006/tobi: adapted from fourier_transform_other_2
!
      use Mpicomm, only: transp_xy

      real, dimension(nx,ny) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension(nxgrid) :: ax
      complex, dimension(nygrid) :: ay
      real, dimension(4*nxgrid+15) :: wsavex
      real, dimension(4*nygrid+15) :: wsavey
      integer :: l,m,ibox
      logical :: lforward

      if (mod(nxgrid,nygrid)/=0) then
        call fatal_error('fourier_transform_xy_parallel', &
                         'nxgrid needs to be an integer multiple of nygrid.')
      endif

      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif

      if (lforward) then

        if (nxgrid>1) then
!
!  Transform x-direction.
!
          call cffti(nxgrid,wsavex)

          do m=1,ny
            ax=cmplx(a_re(:,m),a_im(:,m))
            call cfftf(nxgrid,ax,wsavex)
            a_re(:,m)=real(ax)
            a_im(:,m)=aimag(ax)
          enddo

        endif

        if (nygrid>1) then
!
!  Transform y-direction.
!
          call transp_xy(a_re)
          call transp_xy(a_im)

          call cffti(nygrid,wsavey)

          do ibox=0,nxgrid/nygrid-1
            iy=ibox*nygrid
            do l=1,ny
              ay=cmplx(a_re(iy+1:iy+nygrid,l),a_im(iy+1:iy+nygrid,l))
              call cfftf(nygrid,ay,wsavey)
              a_re(iy+1:iy+nygrid,l)=real(ay)
              a_im(iy+1:iy+nygrid,l)=aimag(ay)
            enddo
          enddo

          call transp_xy(a_re)
          call transp_xy(a_im)

        endif

      else

        if (nygrid>1) then
!
!  Transform y-direction back.
!
          call transp_xy(a_re)
          call transp_xy(a_im)

          call cffti(nygrid,wsavey)

          do ibox=0,nxgrid/nygrid-1
            iy=ibox*nygrid
            do l=1,ny
              ay=cmplx(a_re(iy+1:iy+nygrid,l),a_im(iy+1:iy+nygrid,l))
              call cfftb(nygrid,ay,wsavey)
              a_re(iy+1:iy+nygrid,l)=real(ay)
              a_im(iy+1:iy+nygrid,l)=aimag(ay)
            enddo
          enddo

          call transp_xy(a_re)
          call transp_xy(a_im)

        endif

        if (nxgrid>1) then
!
!  Transform x-direction back.
!
          call cffti(nxgrid,wsavex)

          do m=1,ny
            ax=cmplx(a_re(:,m),a_im(:,m))
            call cfftb(nxgrid,ax,wsavex)
            a_re(:,m)=real(ax)
            a_im(:,m)=aimag(ax)
          enddo

        endif

      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/(nxgrid*nygrid)
        a_im=a_im/(nxgrid*nygrid)
      endif

    endsubroutine fourier_transform_xy_parallel
!***********************************************************************
endmodule Fourier
