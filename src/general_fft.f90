! $Id: general_fft.f90,v 1.1 2006-07-28 11:49:12 ajohan Exp $
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
    subroutine transform(a1,a2,a3,b1,b2,b3)
!
!  Subroutine to do fourier transform
!  The routine overwrites the input data
!
!  03-sep-02/nils: coded
!  05-nov-02/axel: added normalization factor
!
      real,dimension(nx,ny,nz) :: a1,b1,a2,b2,a3,b3

! Doing the x field
      if (lroot .and. ip<10) print*,'transform: doing fft of x-component'
      call fft(a1,b1, nx*ny*nz, nx, nx,-1) ! x-direction
      call transp(a1,'y')
      call transp(b1,'y')
      call fft(a1,b1, nx*ny*nz, nx, nx,-1) ! y-direction
      call transp(a1,'z')
      call transp(b1,'z')
      call fft(a1,b1, nx*ny*nz, nx, nx,-1) ! z-direction

! Doing the y field
      if (lroot .and. ip<10) print*,'transform: doing fft of y-component'
      call fft(a2,b2, nx*ny*nz, nx, nx,-1) ! x-direction
      call transp(a2,'y')
      call transp(b2,'y')
      call fft(a2,b2, nx*ny*nz, nx, nx,-1) ! y-direction
      call transp(a2,'z')
      call transp(b2,'z')
      call fft(a2,b2, nx*ny*nz, nx, nx,-1) ! z-direction

! Doing the z field
      if (lroot .and. ip<10) print*,'transform: doing fft of z-component'
      call fft(a3,b3, nx*ny*nz, nx, nx,-1) ! x-direction
      call transp(a3,'y')
      call transp(b3,'y')
      call fft(a3,b3, nx*ny*nz, nx, nx,-1) ! y-direction
      call transp(a3,'z')
      call transp(b3,'z')
      call fft(a3,b3, nx*ny*nz, nx, nx,-1) ! z-direction
!
!  Normalize
!
      a1=a1/nwgrid; a2=a2/nwgrid; a3=a3/nwgrid
      b1=b1/nwgrid; b2=b2/nwgrid; b3=b3/nwgrid
      if (lroot .and. ip<10) print*,'transform: fft has finished'
!
    endsubroutine transform
!***********************************************************************
    subroutine transform_i(a_re,a_im)
!
!  Subroutine to do fourier transform
!  The routine overwrites the input data
!
!  22-oct-02/axel+tarek: adapted from transform
!
      real,dimension(nx,ny,nz) :: a_re,a_im

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
    endsubroutine transform_i
!***********************************************************************
    subroutine transform_fftpack(a_re,a_im,direction)
!
!  Subroutine to do Fourier transform
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  27-oct-02/axel: adapted from transform_i, for fftpack
!
      real,dimension(nx,ny,nz) :: a_re,a_im
      complex,dimension(nx) :: ax
      real,dimension(4*nx+15) :: wsavex
      integer :: m,n
      integer,optional :: direction
      logical :: lforward=.true.
!
      if (present(direction)) then
        if (direction==-1) then
          lforward=.false.
        else
          lforward=.true.
        endif
      endif
!
!  Need to initialize cfft only once, because we require nxgrid=nygrid=nzgrid
!
      call cffti(nx,wsavex)
!
      if (lforward) then
        if (lroot .and. ip<10) print*, 'transform_fftpack: doing FFTpack in x'
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
            if (lroot) print*, 'transform_fftpack: must have nygrid=nxgrid!'
            call fatal_error('transform_fftpack','')
          endif
          call transp(a_re,'y')
          call transp(a_im,'y')
!
!  The length of the array in the y-direction is nx.
!
          if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in y'
          do n=1,nz; do m=1,ny
            ax=cmplx(a_re(:,m,n),a_im(:,m,n))
            call cfftf(nx,ax,wsavex)
            a_re(:,m,n)=real(ax)
            a_im(:,m,n)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform z-direction.
!      
        if (nzgrid/=1) then
          if (nzgrid/=nxgrid) then
            if (lroot) &
                print*, 'transform_fftpack: must have nzgrid=nxgrid!'
            call fatal_error('transform_fftpack','')
          endif
          call transp(a_re,'z')
          call transp(a_im,'z')
!
!  The length of the array in the z-direction is also nx.
!
          if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in z'
          do n=1,nz; do m=1,ny
            ax=cmplx(a_re(:,m,n),a_im(:,m,n))
            call cfftf(nx,ax,wsavex)
            a_re(:,m,n)=real(ax)
            a_im(:,m,n)=aimag(ax)
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
            if (lroot) print*, 'transform_fftpack: must have nzgrid=nxgrid!'
            call fatal_error('transform_fftpack','')
          endif
!
          if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in z'
          do n=1,nz; do m=1,ny
            ax=cmplx(a_re(:,m,n),a_im(:,m,n))
            call cfftb(nx,ax,wsavex)
            a_re(:,m,n)=real(ax)
            a_im(:,m,n)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform y-direction back. Must transpose to go from (z,x,y) to (y,x,z).
!
        if (nygrid/=1) then
          if (nygrid/=nxgrid) then
            if (lroot) print*, 'transform_fftpack: must have nygrid=nxgrid!'
            call fatal_error('transform_fftpack','')
          endif
!
          if (nzgrid/=1) then
            call transp(a_re,'z')
            call transp(a_im,'z')
          endif
!
          if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in y'
          do n=1,nz; do m=1,ny
            ax=cmplx(a_re(:,m,n),a_im(:,m,n))
            call cfftb(nx,ax,wsavex)
            a_re(:,m,n)=real(ax)
            a_im(:,m,n)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform x-direction back. Transpose to go from (y,x,z) to (x,y,z).
!
        if (lroot .and. ip<10) print*, 'transform_fftpack: doing FFTpack in x'
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
      if (lroot .and. ip<10) print*,'transform_fftpack: fft has finished'
!
    endsubroutine transform_fftpack
!***********************************************************************
    subroutine transform_fftpack_2d(a_re,a_im,dummy)
!
!  Subroutine to do Fourier transform
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  27-oct-02/axel: adapted from transform_fftpack
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      complex, dimension(nx) :: ax
      real, dimension(4*nx+15) :: wsavex
      integer :: m,n
      integer, optional :: dummy
!
!  check whether nxgrid=nygrid=nzgrid
!
      if (nxgrid/=nygrid .or. (nxgrid/=nzgrid .and. nzgrid/=1)) then
        print*,'transform_fftpack: must have nxgrid=nygrid=nzgrid!'
        call fatal_error('must have nxgrid=nygrid=nzgrid!','')
      endif
!
!  need to initialize cfft only once, because nxgrid=nygrid=nzgrid
!
      call cffti(nx,wsavex)
!
      if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in x'
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
      if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in z'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        call cfftf(nx,ax,wsavex)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
!
!  Normalize
!
      a_re=a_re/nwgrid
      a_im=a_im/nwgrid
      if (lroot .and. ip<10) print*,'transform_fftpack: fft has finished'
!
      if (NO_WARN) print*,dummy  !(keep compiler quiet)
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
!  06-feb-03/nils: adapted from transform_fftpack
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      complex, dimension(nx) :: ax
      real, dimension(4*nx+15) :: wsavex
      integer :: m,n
!
!  check whether nxgrid=nygrid=nzgrid
!
      if ( (nygrid/=1.and.nygrid/=nxgrid) .or. &
           (nzgrid/=1.and.nzgrid/=nxgrid) ) then
        if (lroot) &
            print*, 'transform_fftpack_1d: must have nxgrid=nygrid=nzgrid!'
        call fatal_error( &
            'transform_fftpack_1d: must have nxgrid=nygrid=nzgrid!','')
      endif
!
!  need to initialize cfft only once, because nxgrid=nygrid=nzgrid
!
      call cffti(nx,wsavex)
!
      if (lroot .and. ip<10) print*,'transform_fftpack_1d: doing FFTpack in x'
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
      if (lroot .and. ip<10) print*,'transform_fftpack_1d: fft has finished'
!
    endsubroutine transform_fftpack_1d
!***********************************************************************
    subroutine transform_fftpack_shear(a_re,a_im,direction)
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
      integer, optional :: direction
!
      complex, dimension (nxgrid) :: ax
      complex, dimension (nxgrid) :: ay
      complex, dimension (nxgrid) :: az
      real, dimension (4*nxgrid+15) :: wsave
      real :: deltay_x
      integer :: l,m,n
      logical :: lforward=.true.
!
!  if nxgrid/=nygrid/=nzgrid, stop.

      if (nygrid/=nxgrid .and. nygrid /= 1) then
        print*,'transform_fftpack_shear: need to have nygrid=nxgrid if nygrid /=1.'
        call fatal_error('Inconsistency: nygrid/=nxgrid','')
      endif
      if (nzgrid/=nxgrid .and. nzgrid /= 1) then
        print*,'transform_fftpack_shear: need to have nzgrid=nxgrid if nzgrid /=1.'
        call fatal_error('Inconsistency: nzgrid/=nxgrid','')
      endif
!
      if (present(direction)) then
        if (direction==-1) then
          lforward=.false.
        else
          lforward=.true.
        endif
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
              print*, 'doing FFTpack in y, direction=',direction
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
        if (lroot.and.ip<10) print*, 'doing FFTpack in x, direction=',direction
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
              print*, 'doing FFTpack in z, direction=',direction
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
              print*, 'doing FFTpack in z, direction=',direction
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
        if (lroot.and.ip<10) print*, 'doing FFTpack in x, direction=',direction
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
              print*, 'doing FFTpack in y, direction=',direction
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
    endsubroutine transform_fftpack_shear
!***********************************************************************
    subroutine transform_nr(a_re,a_im)
!
!  Subroutine to do Fourier transform using Numerical Recipes routine.
!  Note that this routine requires that nx, ny, and nz are powers of 2.
!  The routine overwrites the input data.
!
!  30-oct-02/axel: adapted from transform_fftpack for Numerical Recipes
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      complex, dimension(nx) :: ax
      integer :: m,n
!
!  This Fourier transform would work, but it's very slow!
!  Even the compilation is very slow, so we better get rid of it!
!
      call fatal_error('transform_nr: currently disabled!','')
!
      if (lroot .and. ip<10) print*,'transform_nr: doing FFT_nr in x'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        !call four1(ax,nx,-1)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
      call transp(a_re,'y')
      call transp(a_im,'y')
!
!  The length of the array in the y-direction is nx
!  (remember: nxgrid=nygrid=nzgrid!)
!
      if (lroot .and. ip<10) print*,'transform_nr: doing FFT_nr in y'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        !call four1(ax,nx,-1)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
      call transp(a_re,'z')
      call transp(a_im,'z')
!
!  The length of the array in the z-direction is also nx
!
      if (lroot .and. ip<10) print*,'transform_nr: doing FFT_nr in z'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        !call four1(ax,nx,-1)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
!
!  Normalize
!
      a_re=a_re/nwgrid
      a_im=a_im/nwgrid
      if (lroot .and. ip<10) print*,'transform_nr: fft has finished'
!
    endsubroutine transform_nr
!***********************************************************************
endmodule General_FFT
