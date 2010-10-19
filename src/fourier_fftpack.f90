! $Id$
!
!  This module contains FFT wrapper subroutines.
!
module Fourier
!
  use Cdata
  use Cparam
  use Messages
  use Mpicomm, only: transp,transp_other
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
  interface fft_xy_parallel
    module procedure fft_xy_parallel_2D
    module procedure fft_xy_parallel_3D
    module procedure fft_xy_parallel_4D
  endinterface
!
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
      real, dimension (nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (nx) :: ax
      real, dimension (4*nx+15) :: wsavex
      integer :: l,m,n
      logical :: lforward
!
      if (nprocx>1) call fatal_error('fourier_transform','Must have nprocx=1!')
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
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
      real, dimension (nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (nx) :: ax
      real, dimension (4*nx+15) :: wsavex
      integer :: l,m,n
      logical :: lforward
!
      if (nprocx>1) &
          call fatal_error('fourier_transform_xy','Must have nprocx=1!')
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
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
      else
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
      real, dimension (nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (nx) :: ax
      real, dimension (4*nx+15) :: wsavex
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
      if (nprocx>1) &
          call fatal_error('fourier_transform_xz','Must have nprocx=1!')
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
      real, dimension (nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (nx) :: ax
      real, dimension (4*nx+15) :: wsavex
      integer :: m,n
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      if (nprocx>1) &
          call fatal_error('fourier_transform_x','Must have nprocx=1!')
!
!  Check whether nxgrid=nygrid=nzgrid.
!
      if ( (nygrid/=1.and.nygrid/=nxgrid) .or. &
           (nzgrid/=1.and.nzgrid/=nxgrid) ) then
        if (lroot) &
            print*, 'fourier_transform_x: must have nxgrid=nygrid=nzgrid!'
        call fatal_error('fourier_transform_x','')
      endif
!
!  need to initialize cfft only once, because nxgrid=nygrid
!
      call cffti(nx,wsavex)
!
      if (lforward) then
!
!  Transform x-direction to fourier space
!
        if (lroot.and.ip<10) print*, 'fourier_transform_x: doing FFTpack in x'
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nx,ax,wsavex)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
!
      else
!
!  Transform x-direction back to real space
!
        if (lroot.and.ip<10) print*, 'fourier_transform_x: doing FFTpack in x'
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nx,ax,wsavex)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
!
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/nxgrid
        a_im=a_im/nxgrid
      endif
!
      if (lroot .and. ip<10) print*, 'fourier_transform_x: fft has finished'
!
    endsubroutine fourier_transform_x
!***********************************************************************
    subroutine fourier_transform_y(a_re,a_im,linv,lnorm)
!
!  Subroutine to do Fourier transform in the y-direction.
!  As it is not cache efficient to Fourier transform in other
!  direction than x, we transpose the array. The original array
!  and the transposed one are shown below in the left and right
!  sides, respectively
!
!
!    ---------------  ^ nygrid     ----------------  ^ nygrid
!   |x z Q W E R T Y| |           |r s m W |i h c Y| |
!   |j k l m n b v c| |           |e a l Q |u g v T| |
!   |---------------| | ^         |----------------| | ^
!   |o p a s d f g h| | | ny      |w p k z |y f b R| | | ny
!   |q w e r t y u i| | |         |q o j x |t d n E| | |
!    ---------------               ----------------
!    --------------->              <nygrid><nygrid>
!                   nx             ----------------> nx
!
!  The transposed array has the elements in the correct order, but
!  for obvious book-keeping reasons, the dimension is still (nx,ny),
!  instead of (nygrid,nx). The fourier transform then can be done, but
!  the calculation has to be split into boxes of dimension (nygrid,ny).
!  Therefore, it only works when nx is a multiple of nygrid.
!
!  For the case nx<nygrid, interpolate the x-values to the dimension of
!  nygrid prior to transposing. Then transform and interpolate back
!
!  07-may-08/wlad: coded
!
      use General, only: spline
!
      real, dimension (nx,ny,nz)     :: a_re,a_im
      real, dimension (nygrid,ny,nz) :: tmp_re,tmp_im
      real, dimension (nygrid),save  :: xnyg
      real, dimension (4*nygrid+15)  :: wsave
      complex, dimension (nygrid)    :: ay
      real    :: dnx
      integer :: l,n,iarr,ix,ido,iup,i
      logical :: lforward,lnormalize,err
      logical, save :: lfirstcall=.true.
      logical, optional :: linv
      logical, optional :: lnorm
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      lnormalize=.true.
      if (present(lnorm)) lnormalize=lnorm
!
      if (nprocx>1) &
          call fatal_error('fourier_transform_y','Must have nprocx=1!')
!
! Separate the problem in two cases. nxgrid>= nygrid and its
! opposite
!
      if (nxgrid >= nygrid) then
        if (mod(nxgrid,nygrid)/=0) then
          print*,'fourier_transform_y: when nxgrid>= nygrid, '//&
               'nxgrid needs to be an integer multiple of nygrid.'
          call fatal_error('fourier_transform_y','mod(nxgrid,nygrid)/=0')
        endif
      endif
!
!  initialize cfft (coefficients for fft?)
!
      call cffti(nygrid,wsave)
!
      if (lroot.and.ip<10) print*, 'fourier_transform_y: doing FFTpack in y'
!
!  Transform differently according to sizes of x and y
!
      if (nxgrid>=nygrid) then
!
!  Transpose, transform and transpose back
!
        if (lroot .and. ip<10) print*, &
             'fourier_transform_y: nxgrid>=nygrid'
!
        call transp(a_re,'y') ; call transp(a_im,'y')
        do n=1,nz; do l=1,ny
!  Divide a_re into arrays of size nygrid to fit ay
          do iarr=0,nxgrid/nygrid-1
            ix=iarr*nygrid ; ido=ix+1 ; iup=ix+nygrid
            ay=cmplx(a_re(ido:iup,l,n),a_im(ido:iup,l,n))
            if (lforward) then
              call cfftf(nygrid,ay,wsave)
            else
              call cfftb(nygrid,ay,wsave)
            endif
            a_re(ido:iup,l,n)=real(ay)
            a_im(ido:iup,l,n)=aimag(ay)
          enddo
        enddo;enddo
        call transp(a_re,'y') ; call transp(a_im,'y')
!
! Normalize if forward
!
        if (lforward) then
          if (lnormalize) then
            a_re=a_re/nygrid
            a_im=a_im/nygrid
          endif
        endif
!
      else !case nxgrid<nygrid
!
        if (lroot .and. ip<10) print*, &
             'fourier_transform_y: nxgrid<nygrid'
!
! Save interpolated values of x to dimension nygrid
!
        if (lfirstcall) then
          dnx=Lxyz(1)/(nygrid-1)
          do i=0,nygrid-1
            xnyg(i+1)=x(l1)+i*dnx
          enddo
          lfirstcall=.false.
        endif
!
! Interpolate nx to the same dimension as nygrid, so we
! can transpose
!
        do n=1,nz;do m=1,ny
          call spline(x(l1:l2),a_re(:,m,n),xnyg,tmp_re(:,m,n),nx,nygrid,err)
          call spline(x(l1:l2),a_im(:,m,n),xnyg,tmp_im(:,m,n),nx,nygrid,err)
        enddo;enddo
!
! Transpose, transform, transpose back
!
        call transp_other(tmp_re,'y') ; call transp_other(tmp_im,'y')
        do n=1,nz;do l=1,ny
          ay=cmplx(tmp_re(:,l,n),tmp_im(:,l,n))
          if (lforward) then
            call cfftf(nygrid,ay,wsave)
          else
            call cfftb(nygrid,ay,wsave)
          endif
          tmp_re(:,l,n)=real(ay)
          tmp_im(:,l,n)=aimag(ay)
        enddo; enddo
        call transp_other(tmp_re,'y') ; call transp_other(tmp_im,'y')
!
! Normalize if forward
!
        if (lforward) then
          if (lnormalize) then
            tmp_re=tmp_re/nygrid
            tmp_im=tmp_im/nygrid
          endif
        endif
!
! Interpolate (coarsen) back to dimension nx
!
        do n=1,nz;do m=1,ny
          call spline(xnyg,tmp_re(:,m,n),x(l1:l2),a_re(:,m,n),nygrid,nx,err)
          call spline(xnyg,tmp_im(:,m,n),x(l1:l2),a_im(:,m,n),nygrid,nx,err)
        enddo;enddo
      endif
!
      if (lroot .and. ip<10) print*, 'fourier_transform_y: fft has finished'
!
    endsubroutine fourier_transform_y
!***********************************************************************
    subroutine fourier_transform_shear(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform in shearing coordinates.
!  The routine overwrites the input data.
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
      real, dimension (nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (nx) :: ax
      complex, dimension (nx) :: ay
      complex, dimension (nx) :: az
      real, dimension (4*nx+15) :: wsave
      real :: deltay_x
      integer :: l,m,n,two
      logical :: lforward
!
      two = 2         ! avoid 'array out of bounds' below for nygrid=1
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      if (nprocx>1) &
          call fatal_error('fourier_transform_shear','Must have nprocx=1!')
!
!  If nxgrid/=nygrid/=nzgrid, stop.
!
      if (nygrid/=nxgrid .and. nygrid/=1) then
        print*, 'fourier_transform_shear: '// &
            'need to have nygrid=nxgrid if nygrid/=1.'
        call fatal_error('fourier_transform_shear','')
      endif
      if (nzgrid/=nxgrid .and. nzgrid/=1) then
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
            ay(two:nxgrid)=ay(two:nxgrid)*exp(cmplx(0.0, ky_fft(two:nxgrid)*deltay_x))
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
            ay(two:nxgrid)=ay(two:nxgrid)*exp(cmplx(0.0,-ky_fft(two:nxgrid)*deltay_x))
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
!  The routine overwrites the input data.
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
      real, dimension (nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (nx) :: ax
      complex, dimension (nx) :: ay
      real, dimension (4*nx+15) :: wsave
      real :: deltay_x
      integer :: l,m,n,two
      logical :: lforward
!
      two = 2         ! avoid 'array out of bounds' below for nygrid=1
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      if (nprocx>1) &
          call fatal_error('fourier_transform_shear_xy','Must have nprocx=1!')
!
!  If nxgrid/=nygrid/=nzgrid, stop.
!
      if (nygrid/=nxgrid .and. nygrid/=1) then
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
        if (nygrid>1) then
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
            ay(two:nxgrid)=ay(two:nxgrid)*exp(cmplx(0.0, ky_fft(two:nxgrid)*deltay_x))
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
        if (nygrid>1) then
          call transp(a_re,'y')
          call transp(a_im,'y')
          if (lroot.and.ip<10) &
              print*, 'fourier_transform_shear: doing FFTpack in y'
          do n=1,nz; do l=1,ny
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
!  Shift y-coordinate back to regular frame (see above).
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            ay(two:nxgrid)=ay(two:nxgrid)*exp(cmplx(0.0,-ky_fft(two:nxgrid)*deltay_x))
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
!  This routine does not operate in parallel, but should be used to Fourier
!  transform an array present in its entirety on the local processor.
!  The routine overwrites the input data.
!
!  28-jul-2006/anders: adapted from fourier_transform
!
      real, dimension (:) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (size(a_re,1)) :: ax
      real, dimension (4*size(a_re,1)+15) :: wsavex
      integer :: nx_other
      logical :: lforward
!
      nx_other=size(a_re,1)
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
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
!  This routine does not operate in parallel, but should be used to Fourier
!  transform an array present in its entirety on the local processor.
!  The routine overwrites the input data.
!
!  28-jul-2006/anders: adapted from fourier_transform_1
!
      real, dimension (:,:) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (size(a_re,1)) :: ax
      complex, dimension (size(a_re,2)) :: ay
      real, dimension (4*size(a_re,1)+15) :: wsavex
      real, dimension (4*size(a_re,2)+15) :: wsavey
      integer :: l, m, nx_other, ny_other
      logical :: lforward
!
      nx_other=size(a_re,1); ny_other=size(a_re,2)
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
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
    subroutine fourier_transform_xy_xy(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do Fourier transform of a 2-D array under MPI.
!  nxgrid is restricted to be an integer multiple of nygrid.
!  If it's not, 'fft_xy_parallel' is called, which has fewer restrictions.
!  The routine overwrites the input data.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!
!   6-oct-2006/tobi: adapted from fourier_transform_other_2
!
      use Mpicomm, only: transp_xy
!
      real, dimension (nx,ny), intent(inout) :: a_re,a_im
      logical, optional, intent(in) :: linv,lneed_im
!
      complex, dimension (nx) :: ax
      complex, dimension (nygrid) :: ay
      real, dimension (4*nx+15) :: wsavex
      real, dimension (4*nygrid+15) :: wsavey
      real, dimension (ny) :: deltay_x
      integer :: l,m,ibox
      logical :: lforward,lcompute_im
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      lcompute_im=.true.
      if (present(lneed_im)) lcompute_im=lneed_im
!
      if ((nprocx>1) .or. (mod(nxgrid,nygrid)/=0)) then
        call fft_xy_parallel(a_re,a_im,linv)
        return
      endif
!
      if (lshear) deltay_x=-deltay*(x(m1+ipy*ny:m2+ipy*ny)-(x0+Lx/2))/Lx
!
      if (lforward) then
!
        if (nygrid>1) then
!
!  Transform y-direction.
!
          call transp_xy(a_re)
          if (lcompute_im) then
            call transp_xy(a_im)
          else
            a_im=0.0
          endif
!
          call cffti(nygrid,wsavey)
!
          do ibox=0,nxgrid/nygrid-1
            iy=ibox*nygrid
            do l=1,ny
              ay=cmplx(a_re(iy+1:iy+nygrid,l),a_im(iy+1:iy+nygrid,l))
              call cfftf(nygrid,ay,wsavey)
              if (lshear) ay = ay*exp(cmplx(0.,+ky_fft*deltay_x(l)))
              a_re(iy+1:iy+nygrid,l)=real(ay)
              a_im(iy+1:iy+nygrid,l)=aimag(ay)
            enddo
          enddo
!
          call transp_xy(a_re)
          call transp_xy(a_im)
!
        endif
!
        if (nxgrid>1) then
!
!  Transform x-direction.
!
          call cffti(nxgrid,wsavex)
!
          do m=1,ny
            ax=cmplx(a_re(:,m),a_im(:,m))
            call cfftf(nxgrid,ax,wsavex)
            a_re(:,m)=real(ax)
            a_im(:,m)=aimag(ax)
          enddo
!
        endif
!
      else
!
        if (nxgrid>1) then
!
!  Transform x-direction back.
!
          call cffti(nxgrid,wsavex)
!
          do m=1,ny
            ax=cmplx(a_re(:,m),a_im(:,m))
            call cfftb(nxgrid,ax,wsavex)
            a_re(:,m)=real(ax)
            a_im(:,m)=aimag(ax)
          enddo
!
        endif
!
        if (nygrid>1) then
!
!  Transform y-direction back.
!
          call transp_xy(a_re)
          call transp_xy(a_im)
!
          call cffti(nygrid,wsavey)
!
          do ibox=0,nxgrid/nygrid-1
            iy=ibox*nygrid
            do l=1,ny
              ay=cmplx(a_re(iy+1:iy+nygrid,l),a_im(iy+1:iy+nygrid,l))
              if (lshear) ay = ay*exp(cmplx(0.,-ky_fft*deltay_x(l)))
              call cfftb(nygrid,ay,wsavey)
              a_re(iy+1:iy+nygrid,l)=real(ay)
              if (lcompute_im) a_im(iy+1:iy+nygrid,l)=aimag(ay)
            enddo
          enddo
!
          call transp_xy(a_re)
          if (lcompute_im) call transp_xy(a_im)
!
        endif
!
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/(nxgrid*nygrid)
        a_im=a_im/(nxgrid*nygrid)
      endif
!
    endsubroutine fourier_transform_xy_xy
!***********************************************************************
    subroutine fourier_transform_xy_xy_other(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform of a 2-D array
!  of arbitrary size under MPI.
!  The routine overwrites the input data.
!
!  15-jan-2008/wlad: adapted from fourier_transform_xy_xy
!                    and fourier_transform_other_2
!
      use Mpicomm, only: transp_xy_other
!
      real, dimension (:,:) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (size(a_re,1)) :: ax
      complex, dimension (nprocy*size(a_re,2)) :: ay
      real, dimension (4*size(a_re,1)+15) :: wsavex
      real, dimension (4*nprocy*size(a_re,2)+15) :: wsavey
      integer :: l,m,nx_other,ny_other
      integer :: nxgrid_other,nygrid_other
      logical :: lforward
!
      nx_other=size(a_re,1); ny_other=size(a_re,2)
      nxgrid_other=nx_other
      nygrid_other=ny_other*nprocy
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      if (nxgrid_other/=nygrid_other) &
        call fatal_error('fourier_transform_xy_xy_other', &
             'nxgrid_other needs to be equal to nygrid_other.',lfirst_proc_xy)
!
      if (lforward) then
        if (nygrid_other > 1) then
!
!  Transform y-direction.
!
          call transp_xy_other(a_re)
          call transp_xy_other(a_im)
!
          call cffti(nygrid_other,wsavey)
!
          do l=1,ny_other
            ay=cmplx(a_re(:,l),a_im(:,l))
            call cfftf(nygrid_other,ay,wsavey)
            a_re(:,l)=real(ay)
            a_im(:,l)=aimag(ay)
          enddo
!
          call transp_xy_other(a_re)
          call transp_xy_other(a_im)
!
      endif
!
      if (nxgrid > 1) then
!
!  Transform x-direction
!
        call cffti(nxgrid_other,wsavex)
!
        do m=1,ny_other
          ax=cmplx(a_re(:,m),a_im(:,m))
          call cfftf(nxgrid_other,ax,wsavex)
          a_re(:,m)=real(ax)
          a_im(:,m)=aimag(ax)
        enddo
!
      endif
!
    else
!
        if (nxgrid_other>1) then
!
!  Transform x-direction back.
!
          call cffti(nxgrid_other,wsavex)
!
          do m=1,ny_other
            ax=cmplx(a_re(:,m),a_im(:,m))
            call cfftb(nxgrid_other,ax,wsavex)
            a_re(:,m)=real(ax)
            a_im(:,m)=aimag(ax)
          enddo
!
        endif
!
        if (nygrid_other>1) then
!
!  Transform y-direction back.
!
          call transp_xy_other(a_re)
          call transp_xy_other(a_im)
!
          call cffti(nygrid_other,wsavey)
!
          do l=1,ny_other
            ay=cmplx(a_re(:,l),a_im(:,l))
            call cfftb(nygrid_other,ay,wsavey)
            a_re(:,l)=real(ay)
            a_im(:,l)=aimag(ay)
          enddo
!
          call transp_xy_other(a_re)
          call transp_xy_other(a_im)
!
        endif
!
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/(nxgrid_other*nygrid_other)
        a_im=a_im/(nxgrid_other*nygrid_other)
      endif
!
    endsubroutine fourier_transform_xy_xy_other
!***********************************************************************
    subroutine fft_xy_parallel_2D(a_re,a_im,linv,lneed_im,shift_y)
!
!  Subroutine to do FFT of distributed 2D data in the x- and y-direction.
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  17-aug-2010/Bourdin.KIS: adapted from fft_xy_parallel_4D
!  09-sep-2010/wlad: added possibility of shift
!
      use Mpicomm, only: remap_to_pencil_xy, transp_pencil_xy, unmap_from_pencil_xy
!
      real, dimension (nx,ny), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
      real, dimension (nx), optional :: shift_y
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
      real, dimension (:,:), allocatable :: p_re, p_im   ! data in pencil shape
      real, dimension (:,:), allocatable :: t_re, t_im   ! data in transposed pencil shape
      complex, dimension (nxgrid) :: ax
      complex, dimension (nygrid) :: ay
      real, dimension (4*nxgrid+15) :: wsavex
      real, dimension (4*nygrid+15) :: wsavey
      real, dimension (tny) :: deltay_x
      real, dimension (tny) :: dshift_y
!
      integer :: l, m, stat, x_offset
      logical :: lforward, lcompute_im
      logical :: lshift
!
      lforward = .true.
      if (present (linv)) lforward = .not.linv
!
      lshift=.false.
      if (present (shift_y)) lshift=.true.
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error('fft_xy_parallel_2D', &
          'nxgrid needs to be an integer multiple of nprocx*nprocy', &
          lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error('fft_xy_parallel_2D', &
          'nygrid needs to be an integer multiple of nprocx*nprocy', &
          lfirst_proc_xy)
!
!  Allocate memory for large arrays.
!
      allocate (p_re(pnx,pny), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_2D', &
          'Could not allocate memory for p_re', .true.)
      allocate (p_im(pnx,pny), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_2D', &
          'Could not allocate memory for p_im', .true.)
      allocate (t_re(tnx,tny), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_2D', &
          'Could not allocate memory for t_re', .true.)
      allocate (t_im(tnx,tny), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_2D', &
          'Could not allocate memory for t_im', .true.)
!
      if (lshear) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        deltay_x = -deltay * (xgrid(x_offset:x_offset+tny-1) - (x0+Lx/2))/Lx
      endif
!
      if (lshift) then
        x_offset = 1+ipy*tny
        dshift_y = shift_y(x_offset:x_offset+tny-1)
      endif
!
      call cffti(nxgrid, wsavex)
      call cffti(nygrid, wsavey)
!
      if (lforward) then
!
!  Forward FFT.
!
!  Remap the data we need into pencil shape.
!
        call remap_to_pencil_xy(a_re, p_re)
        if (lcompute_im) then
          call remap_to_pencil_xy(a_im, p_im)
        else
          p_im = 0.0
        endif
!
        if (nxgrid > 1) then
          do m = 1, pny
!
!  Transform x-direction.
!
            ax = cmplx (p_re(:,m), p_im(:,m))
            call cfftf(nxgrid, ax, wsavex)
            p_re(:,m) = real (ax)
            p_im(:,m) = aimag (ax)
          enddo
        endif
!
        if (nygrid > 1) then
          call transp_pencil_xy(p_re, t_re)
          call transp_pencil_xy(p_im, t_im)
!
          do l = 1, tny
!
!  Transform y-direction.
!
            ay = cmplx (t_re(:,l), t_im(:,l))
            call cfftf(nygrid, ay, wsavey)
            if (lshear) ay = ay * exp (cmplx (0, ky_fft * deltay_x(l)))
            if (lshift) ay = ay * exp (cmplx (0,-ky_fft * dshift_y(l)))
!
            t_re(:,l) = real (ay)
            t_im(:,l) = aimag (ay)
          enddo
!
          call transp_pencil_xy(t_re, p_re)
          call transp_pencil_xy(t_im, p_im)
        endif
!
!  Unmap the results back to normal shape.
!
        call unmap_from_pencil_xy(p_re, a_re)
        call unmap_from_pencil_xy(p_im, a_im)
!
!  Apply normalization factor to fourier coefficients.
!
        a_re = a_re / (nxgrid*nygrid)
        a_im = a_im / (nxgrid*nygrid)
!
      else
!
!  Inverse FFT.
!
!  Remap the data we need into transposed pencil shape.
!
        call remap_to_pencil_xy(a_re, p_re)
        call remap_to_pencil_xy(a_im, p_im)
!
        if (nygrid > 1) then
          call transp_pencil_xy(p_re, t_re)
          call transp_pencil_xy(p_im, t_im)
!
          do l = 1, tny
!
!  Transform y-direction back.
!
            ay = cmplx (t_re(:,l), t_im(:,l))
            if (lshear) ay = ay * exp (cmplx (0, -ky_fft*deltay_x(l)))
            call cfftb (nygrid, ay, wsavey)
            t_re(:,l) = real (ay)
            t_im(:,l) = aimag (ay)
          enddo
!
          call transp_pencil_xy(t_re, p_re)
          call transp_pencil_xy(t_im, p_im)
        endif
!
        if (nxgrid > 1) then
          do m = 1, pny
!
!  Transform x-direction back.
!
            ax = cmplx (p_re(:,m), p_im(:,m))
            call cfftb (nxgrid, ax, wsavex)
            p_re(:,m) = real (ax)
            if (lcompute_im) p_im(:,m) = aimag (ax)
          enddo
        endif
!
!  Unmap the results back to normal shape.
!
        call unmap_from_pencil_xy(p_re, a_re)
        if (lcompute_im) call unmap_from_pencil_xy(p_im, a_im)
!
      endif
!
!  Deallocate large arrays.
!
      if (allocated (p_re)) deallocate (p_re)
      if (allocated (p_im)) deallocate (p_im)
      if (allocated (t_re)) deallocate (t_re)
      if (allocated (t_im)) deallocate (t_im)
!
    endsubroutine fft_xy_parallel_2D
!***********************************************************************
    subroutine fft_xy_parallel_3D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 3D data in the x- and y-direction.
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  17-aug-2010/Bourdin.KIS: adapted from fft_xy_parallel_4D
!
      use Mpicomm, only: remap_to_pencil_xy, transp_pencil_xy, unmap_from_pencil_xy
!
      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
      real, dimension (:,:,:), allocatable :: p_re, p_im   ! data in pencil shape
      real, dimension (:,:,:), allocatable :: t_re, t_im   ! data in transposed pencil shape
      complex, dimension (nxgrid) :: ax
      complex, dimension (nygrid) :: ay
      real, dimension (4*nxgrid+15) :: wsavex
      real, dimension (4*nygrid+15) :: wsavey
      real, dimension (tny) :: deltay_x
      integer :: inz ! size of the third dimension
      integer :: l, m, stat, x_offset, pos_z
      logical :: lforward, lcompute_im
!
      lforward = .true.
      if (present (linv)) lforward = .not.linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      inz = size (a_re, 3)
!
      if (inz /= size (a_im, 3)) &
          call fatal_error('fft_xy_parallel_3D', &
          'third dimension differs for real and imaginary part', &
          lfirst_proc_xy)
!
      if ((size (a_re, 1) /= nx) .or. (size (a_re, 2) /= ny)) &
          call fatal_error('fft_xy_parallel_3D', &
          'real array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (a_im, 1) /= nx) .or. (size (a_im, 2) /= ny)) &
          call fatal_error('fft_xy_parallel_3D', &
          'imaginary array size mismatch /= nx,ny', lfirst_proc_xy)
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error('fft_xy_parallel_3D', &
          'nxgrid needs to be an integer multiple of nprocx*nprocy', &
          lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error('fft_xy_parallel_3D', &
          'nygrid needs to be an integer multiple of nprocx*nprocy', &
          lfirst_proc_xy)
!
!  Allocate memory for large arrays.
!
      allocate (p_re(pnx,pny,inz), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_3D', &
          'Could not allocate memory for p_re', .true.)
      allocate (p_im(pnx,pny,inz), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_3D', &
          'Could not allocate memory for p_im', .true.)
      allocate (t_re(tnx,tny,inz), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_3D', &
          'Could not allocate memory for t_re', .true.)
      allocate (t_im(tnx,tny,inz), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_3D', &
          'Could not allocate memory for t_im', .true.)
!
      if (lshear) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        deltay_x = -deltay * (xgrid(x_offset:x_offset+tny-1) - (x0+Lx/2))/Lx
      endif
!
      call cffti(nxgrid, wsavex)
      call cffti(nygrid, wsavey)
!
      if (lforward) then
!
!  Forward FFT.
!
!  Remap the data we need into pencil shape.
!
        call remap_to_pencil_xy(a_re, p_re)
        if (lcompute_im) then
          call remap_to_pencil_xy(a_im, p_im)
        else
          p_im = 0.0
        endif
!
        if (nxgrid > 1) then
          do pos_z = 1, inz
            do m = 1, pny
!
!  Transform x-direction.
!
              ax = cmplx (p_re(:,m,pos_z), p_im(:,m,pos_z))
              call cfftf(nxgrid, ax, wsavex)
              p_re(:,m,pos_z) = real (ax)
              p_im(:,m,pos_z) = aimag (ax)
            enddo
          enddo
        endif
!
        if (nygrid > 1) then
          call transp_pencil_xy(p_re, t_re)
          call transp_pencil_xy(p_im, t_im)
!
          do pos_z = 1, inz
            do l = 1, tny
!
!  Transform y-direction.
!
              ay = cmplx (t_re(:,l,pos_z), t_im(:,l,pos_z))
              call cfftf(nygrid, ay, wsavey)
              if (lshear) ay = ay * exp (cmplx (0, ky_fft * deltay_x(l)))
              t_re(:,l,pos_z) = real (ay)
              t_im(:,l,pos_z) = aimag (ay)
            enddo
          enddo
!
          call transp_pencil_xy(t_re, p_re)
          call transp_pencil_xy(t_im, p_im)
        endif
!
!  Unmap the results back to normal shape.
!
        call unmap_from_pencil_xy(p_re, a_re)
        call unmap_from_pencil_xy(p_im, a_im)
!
!  Apply normalization factor to fourier coefficients.
!
        a_re = a_re / (nxgrid*nygrid)
        a_im = a_im / (nxgrid*nygrid)
!
      else
!
!  Inverse FFT.
!
!  Remap the data we need into transposed pencil shape.
!
        call remap_to_pencil_xy(a_re, p_re)
        call remap_to_pencil_xy(a_im, p_im)
!
        if (nygrid > 1) then
          call transp_pencil_xy(p_re, t_re)
          call transp_pencil_xy(p_im, t_im)
!
          do pos_z = 1, inz
            do l = 1, tny
!
!  Transform y-direction back.
!
              ay = cmplx (t_re(:,l,pos_z), t_im(:,l,pos_z))
              if (lshear) ay = ay * exp (cmplx (0, -ky_fft*deltay_x(l)))
              call cfftb (nygrid, ay, wsavey)
              t_re(:,l,pos_z) = real (ay)
              t_im(:,l,pos_z) = aimag (ay)
            enddo
          enddo
!
          call transp_pencil_xy(t_re, p_re)
          call transp_pencil_xy(t_im, p_im)
        endif
!
        if (nxgrid > 1) then
          do pos_z = 1, inz
            do m = 1, pny
!
!  Transform x-direction back.
!
              ax = cmplx (p_re(:,m,pos_z), p_im(:,m,pos_z))
              call cfftb (nxgrid, ax, wsavex)
              p_re(:,m,pos_z) = real (ax)
              if (lcompute_im) p_im(:,m,pos_z) = aimag (ax)
            enddo
          enddo
        endif
!
!  Unmap the results back to normal shape.
!
        call unmap_from_pencil_xy(p_re, a_re)
        if (lcompute_im) call unmap_from_pencil_xy(p_im, a_im)
!
      endif
!
!  Deallocate large arrays.
!
      if (allocated (p_re)) deallocate (p_re)
      if (allocated (p_im)) deallocate (p_im)
      if (allocated (t_re)) deallocate (t_re)
      if (allocated (t_im)) deallocate (t_im)
!
    endsubroutine fft_xy_parallel_3D
!***********************************************************************
    subroutine fft_xy_parallel_4D(a_re,a_im,linv,lneed_im)
!
!  Subroutine to do FFT of distributed 4D data in the x- and y-direction.
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  28-jul-2010/Bourdin.KIS: backport from vect_pot_extrapol_z_parallel
!
      use Mpicomm, only: remap_to_pencil_xy, transp_pencil_xy, unmap_from_pencil_xy
!
      real, dimension (:,:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
      real, dimension (:,:,:,:), allocatable :: p_re, p_im ! data in pencil shape
      real, dimension (:,:,:,:), allocatable :: t_re, t_im ! data in transposed pencil shape
      complex, dimension (nxgrid) :: ax
      complex, dimension (nygrid) :: ay
      real, dimension (4*nxgrid+15) :: wsavex
      real, dimension (4*nygrid+15) :: wsavey
      real, dimension (tny) :: deltay_x
      integer :: inz, ina ! size of the third dimension
      integer :: l, m, stat, x_offset, pos_z, pos_a
      logical :: lforward, lcompute_im
!
      lforward = .true.
      if (present (linv)) lforward = .not.linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      inz = size (a_re, 3)
      ina = size (a_re, 4)
!
      if (inz /= size (a_im, 3)) &
          call fatal_error('fft_xy_parallel_4D', &
          'third dimension differs for real and imaginary part', &
          lfirst_proc_xy)
      if (ina /= size (a_im, 4)) &
          call fatal_error('fft_xy_parallel_4D', &
          'fourth dimension differs for real and imaginary part', &
          lfirst_proc_xy)
!
      if ((size (a_re, 1) /= nx) .or. (size (a_re, 2) /= ny)) &
          call fatal_error('fft_xy_parallel_4D', &
          'real array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (a_im, 1) /= nx) .or. (size (a_im, 2) /= ny)) &
          call fatal_error('fft_xy_parallel_4D', &
          'imaginary array size mismatch /= nx,ny', lfirst_proc_xy)
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error('fft_xy_parallel_4D', &
          'nxgrid needs to be an integer multiple of nprocx*nprocy', &
          lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error('fft_xy_parallel_4D', &
          'nygrid needs to be an integer multiple of nprocx*nprocy', &
          lfirst_proc_xy)
!
!  Allocate memory for large arrays.
!
      allocate (p_re(pnx,pny,inz,ina), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_4D', &
          'Could not allocate memory for p_re', .true.)
      allocate (p_im(pnx,pny,inz,ina), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_4D', &
          'Could not allocate memory for p_im', .true.)
      allocate (t_re(tnx,tny,inz,ina), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_4D', &
          'Could not allocate memory for t_re', .true.)
      allocate (t_im(tnx,tny,inz,ina), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_4D', &
          'Could not allocate memory for t_im', .true.)
!
      if (lshear) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        deltay_x = -deltay * (xgrid(x_offset:x_offset+tny-1) - (x0+Lx/2))/Lx
      endif
!
      call cffti(nxgrid, wsavex)
      call cffti(nygrid, wsavey)
!
      if (lforward) then
!
!  Forward FFT.
!
!  Remap the data we need into pencil shape.
!
        call remap_to_pencil_xy(a_re, p_re)
        if (lcompute_im) then
          call remap_to_pencil_xy(a_im, p_im)
        else
          p_im = 0.0
        endif
!
        if (nxgrid > 1) then
          do pos_a = 1, ina
            do pos_z = 1, inz
              do m = 1, pny
!
!  Transform x-direction.
!
                ax = cmplx (p_re(:,m,pos_z,pos_a), p_im(:,m,pos_z,pos_a))
                call cfftf(nxgrid, ax, wsavex)
                p_re(:,m,pos_z,pos_a) = real (ax)
                p_im(:,m,pos_z,pos_a) = aimag (ax)
              enddo
            enddo
          enddo
        endif
!
        if (nygrid > 1) then
          call transp_pencil_xy(p_re, t_re)
          call transp_pencil_xy(p_im, t_im)
!
          do pos_a = 1, ina
            do pos_z = 1, inz
              do l = 1, tny
!
!  Transform y-direction.
!
                ay = cmplx (t_re(:,l,pos_z,pos_a), t_im(:,l,pos_z,pos_a))
                call cfftf(nygrid, ay, wsavey)
                if (lshear) ay = ay * exp (cmplx (0, ky_fft * deltay_x(l)))
                t_re(:,l,pos_z,pos_a) = real (ay)
                t_im(:,l,pos_z,pos_a) = aimag (ay)
              enddo
            enddo
          enddo
!
          call transp_pencil_xy(t_re, p_re)
          call transp_pencil_xy(t_im, p_im)
        endif
!
!  Unmap the results back to normal shape.
!
        call unmap_from_pencil_xy(p_re, a_re)
        call unmap_from_pencil_xy(p_im, a_im)
!
!  Apply normalization factor to fourier coefficients.
!
        a_re = a_re / (nxgrid*nygrid)
        a_im = a_im / (nxgrid*nygrid)
!
      else
!
!  Inverse FFT.
!
!  Remap the data we need into transposed pencil shape.
!
        call remap_to_pencil_xy(a_re, p_re)
        call remap_to_pencil_xy(a_im, p_im)
!
        if (nygrid > 1) then
          call transp_pencil_xy(p_re, t_re)
          call transp_pencil_xy(p_im, t_im)
!
          do pos_a = 1, ina
            do pos_z = 1, inz
              do l = 1, tny
!
!  Transform y-direction back.
!
                ay = cmplx (t_re(:,l,pos_z,pos_a), t_im(:,l,pos_z,pos_a))
                if (lshear) ay = ay * exp (cmplx (0, -ky_fft*deltay_x(l)))
                call cfftb (nygrid, ay, wsavey)
                t_re(:,l,pos_z,pos_a) = real (ay)
                t_im(:,l,pos_z,pos_a) = aimag (ay)
              enddo
            enddo
          enddo
!
          call transp_pencil_xy(t_re, p_re)
          call transp_pencil_xy(t_im, p_im)
        endif
!
        if (nxgrid > 1) then
          do pos_a = 1, ina
            do pos_z = 1, inz
              do m = 1, pny
!
!  Transform x-direction back.
!
                ax = cmplx (p_re(:,m,pos_z,pos_a), p_im(:,m,pos_z,pos_a))
                call cfftb (nxgrid, ax, wsavex)
                p_re(:,m,pos_z,pos_a) = real (ax)
                if (lcompute_im) p_im(:,m,pos_z,pos_a) = aimag (ax)
              enddo
            enddo
          enddo
        endif
!
!  Unmap the results back to normal shape.
!
        call unmap_from_pencil_xy(p_re, a_re)
        if (lcompute_im) call unmap_from_pencil_xy(p_im, a_im)
!
      endif
!
!  Deallocate large arrays.
!
      if (allocated (p_re)) deallocate (p_re)
      if (allocated (p_im)) deallocate (p_im)
      if (allocated (t_re)) deallocate (t_re)
      if (allocated (t_im)) deallocate (t_im)
!
    endsubroutine fft_xy_parallel_4D
!***********************************************************************
    subroutine vect_pot_extrapol_z_parallel(in,out,factor)
!
!  Subroutine to do a z-extrapolation of a vector potential using
!  'factor' as a multiplication factor to the Fourier coefficients.
!  The normalization needs to be already included in 'factor'.
!  Backwards and forwards transforms are done efficiently in one go.
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!
!   7-jul-2010/Bourdin.KIS: coded, adapted parts of bc_aa_pot2 and mdi_init
!
      use Mpicomm, only: remap_to_pencil_xy, transp_pencil_xy, unmap_from_pencil_xy
!
      real, dimension (:,:,:), intent(in) :: in
      real, dimension (:,:,:,:), intent(out) :: out
      real, dimension (:,:,:), intent(in) :: factor
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
      real, dimension (:,:,:), allocatable :: p_re, p_im   ! data in pencil shape
      real, dimension (:,:,:), allocatable :: t_re, t_im   ! data in transposed pencil shape
      real, dimension (:,:,:,:), allocatable :: e_re, e_im ! extrapolated data in transposed pencil shape
      real, dimension (:,:,:,:), allocatable :: b_re, b_im ! backtransformed data in pencil shape
      complex, dimension (nxgrid) :: ax
      complex, dimension (nygrid) :: ay, ay_extra
      real, dimension (4*nxgrid+15) :: wsavex
      real, dimension (4*nygrid+15) :: wsavey
      integer :: ina ! number of components in the output data (usually 3)
      integer :: onz, ona ! number of ghost cells and components in the output data (usually 3)
      integer :: l, m, stat, pos_a, pos_z
!
!
      ina = size (in, 3)
      onz = size (out, 3)
      ona = size (out, 4)
!
      if (ina /= ona) &
          call fatal_error('vect_pot_extrapol_z_parallel', &
          'number of components is different for input and ouput arrays', &
          lfirst_proc_xy)
!
      if ((size (in, 1) /= nx) .or. (size (in, 2) /= ny)) &
          call fatal_error('vect_pot_extrapol_z_parallel', &
          'input array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (out, 1) /= nx) .or. (size (out, 2) /= ny)) &
          call fatal_error('vect_pot_extrapol_z_parallel', &
          'output array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (factor, 1) /= tnx) .or. (size (factor, 2) /= tny)) &
          call fatal_error('vect_pot_extrapol_z_parallel', &
          'factor array size mismatch /= tnx,tny', lfirst_proc_xy)
      if (size (factor, 3) /= onz) &
          call fatal_error('vect_pot_extrapol_z_parallel', &
          'number of ghost cells differs between multiplication factor '// &
          'and ouput array', lfirst_proc_xy)
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error('vect_pot_extrapol_z_parallel', &
          'nxgrid needs to be an integer multiple of nprocx*nprocy', &
          lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error('vect_pot_extrapol_z_parallel', &
          'nygrid needs to be an integer multiple of nprocx*nprocy', &
          lfirst_proc_xy)
!
      if (lshear) call fatal_error('vect_pot_extrapol_z_parallel', &
          'shearing is not implemented in this routine!', lfirst_proc_xy)
!
!  Allocate memory for large arrays.
!
      allocate (p_re(pnx,pny,ona), stat=stat)
      if (stat > 0) call fatal_error('vect_pot_extrapol_z_parallel', &
          'Could not allocate memory for p_re', .true.)
      allocate (p_im(pnx,pny,ona), stat=stat)
      if (stat > 0) call fatal_error('vect_pot_extrapol_z_parallel', &
          'Could not allocate memory for p_im', .true.)
      allocate (t_re(tnx,tny,ona), stat=stat)
      if (stat > 0) call fatal_error('vect_pot_extrapol_z_parallel', &
          'Could not allocate memory for t_re', .true.)
      allocate (t_im(tnx,tny,ona), stat=stat)
      if (stat > 0) call fatal_error('vect_pot_extrapol_z_parallel', &
          'Could not allocate memory for t_im', .true.)
!
      call cffti(nxgrid, wsavex)
      call cffti(nygrid, wsavey)
!
!  Collect the data we need.
!
      call remap_to_pencil_xy(in, p_re)
      p_im = 0.0
!
      do pos_a = 1, ona
        do m = 1, pny
!
!  Transform x-direction.
!
          ax = cmplx (p_re(:,m,pos_a), p_im(:,m,pos_a))
          call cfftf(nxgrid, ax, wsavex)
          p_re(:,m,pos_a) = real (ax)
          p_im(:,m,pos_a) = aimag (ax)
        enddo
      enddo
!
      call transp_pencil_xy(p_re, t_re)
      call transp_pencil_xy(p_im, t_im)
!
      if (allocated (p_re)) deallocate (p_re)
      if (allocated (p_im)) deallocate (p_im)
!
      allocate (e_re(tnx,tny,onz,ona), stat=stat)
      if (stat > 0) call fatal_error('vect_pot_extrapol_z_parallel', &
          'Could not allocate memory for e_re', .true.)
!
      allocate (e_im(tnx,tny,onz,ona), stat=stat)
      if (stat > 0) call fatal_error('vect_pot_extrapol_z_parallel', &
          'Could not allocate memory for e_im', .true.)
!
      do pos_a = 1, ona
        do l = 1, tny
!
!  Transform y-direction.
!
          ay = cmplx (t_re(:,l,pos_a), t_im(:,l,pos_a))
          call cfftf(nygrid, ay, wsavey)
!
!  Transform y-direction back in each z layer.
!
          do pos_z = 1, onz
!
!  Apply factor to fourier coefficients.
!
            ay_extra = ay * factor(:,l,pos_z)
            call cfftb (nygrid, ay_extra, wsavey)
            e_re(:,l,pos_z,pos_a) = real (ay_extra)
            e_im(:,l,pos_z,pos_a) = aimag (ay_extra)
          enddo
        enddo
      enddo
!
      if (allocated (t_re)) deallocate (t_re)
      if (allocated (t_im)) deallocate (t_im)
!
      allocate (b_re(pnx,pny,onz,ona), stat=stat)
      if (stat > 0) call fatal_error('vect_pot_extrapol_z_parallel', &
          'Could not allocate memory for b_re', .true.)
!
      allocate (b_im(pnx,pny,onz,ona), stat=stat)
      if (stat > 0) call fatal_error('vect_pot_extrapol_z_parallel', &
          'Could not allocate memory for b_im', .true.)
!
      call transp_pencil_xy(e_re, b_re)
      call transp_pencil_xy(e_im, b_im)
!
      if (allocated (e_re)) deallocate (e_re)
      if (allocated (e_im)) deallocate (e_im)
!
      do pos_a = 1, ona
        do pos_z = 1, onz
          do m = 1, pny
!
!  Transform x-direction back in each z layer.
!
            ax = cmplx (b_re(:,m,pos_z,pos_a), b_im(:,m,pos_z,pos_a))
            call cfftb (nxgrid, ax, wsavex)
            b_re(:,m,pos_z,pos_a) = real (ax)
          enddo
        enddo
      enddo
!
!  Distribute the results back in normal shape.
!
      call unmap_from_pencil_xy(b_re, out)
!
      if (allocated (b_re)) deallocate (b_re)
      if (allocated (b_im)) deallocate (b_im)
!
    endsubroutine vect_pot_extrapol_z_parallel
!***********************************************************************
    subroutine field_extrapol_z_parallel(in,out,factor)
!
!  Subroutine to do a z-extrapolation of a fields z-component using
!  'factor' as a multiplication factor to the Fourier coefficients.
!  The normalization needs to be already included in 'factor'.
!  'in' and 'factor' are assumed to be already in pencil shape.
!  Backwards and forwards transforms are done efficiently in one go.
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!
!   7-jul-2010/Bourdin.KIS: coded, adapted parts of bc_aa_pot2 and mdi_init
!
      use Mpicomm, only: remap_to_pencil_xy, transp_pencil_xy, unmap_from_pencil_xy
!
      real, dimension (:,:), intent(in) :: in
      real, dimension (:,:,:,:), intent(out) :: out
      real, dimension (:,:,:), intent(in) :: factor
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
      real, dimension (:,:), allocatable :: p_re, p_im     ! data in pencil shape
      real, dimension (:,:), allocatable :: t_re, t_im     ! data in transposed pencil shape
      real, dimension (:,:,:,:), allocatable :: e_re, e_im ! extrapolated data in transposed pencil shape
      real, dimension (:,:,:,:), allocatable :: b_re, b_im ! backtransformed data in pencil shape
      complex, dimension (nxgrid) :: ax
      complex, dimension (nygrid) :: ay, ay_extra, ay_extra_x, ay_extra_y
      real, dimension (4*nxgrid+15) :: wsavex
      real, dimension (4*nygrid+15) :: wsavey
      integer :: onz ! number of layers in the output data
      integer :: l, m, stat, pos_z
!
      onz = size (out, 3)
!
      if ((size (in, 1) /= pnx) .or. (size (in, 2) /= pny)) &
          call fatal_error('field_extrapol_z_parallel', &
          'input array size mismatch /= pnx,pny', lfirst_proc_xy)
      if ((size (out, 1) /= nx) .or. (size (out, 2) /= ny)) &
          call fatal_error('field_extrapol_z_parallel', &
          'output array size mismatch /= nx,ny', lfirst_proc_xy)
      if (size (out, 4) /= 2) &
          call fatal_error('field_extrapol_z_parallel', &
          'output array must have two components', lfirst_proc_xy)
      if ((size (factor, 1) /= tnx) .or. (size (factor, 2) /= tny)) &
          call fatal_error('field_extrapol_z_parallel', &
          'factor array size mismatch /= pnx,pny', lfirst_proc_xy)
      if (size (factor, 3) /= onz) &
          call fatal_error('field_extrapol_z_parallel', &
          'number of ghost cells differs between multiplication factor '// &
          'and ouput array', lfirst_proc_xy)
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error('field_extrapol_z_parallel', &
          'nxgrid needs to be an integer multiple of nprocx*nprocy', &
          lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error('field_extrapol_z_parallel', &
          'nygrid needs to be an integer multiple of nprocx*nprocy', &
          lfirst_proc_xy)
!
      if (lshear) call fatal_error('field_extrapol_z_parallel', &
          'shearing is not implemented in this routine!', lfirst_proc_xy)
!
!  Allocate memory for large arrays.
!
      allocate (p_re(pnx,pny), stat=stat)
      if (stat > 0) call fatal_error('field_extrapol_z_parallel', &
          'Could not allocate memory for p_re', .true.)
      allocate (p_im(pnx,pny), stat=stat)
      if (stat > 0) call fatal_error('field_extrapol_z_parallel', &
          'Could not allocate memory for p_im', .true.)
      allocate (t_re(tnx,tny), stat=stat)
      if (stat > 0) call fatal_error('field_extrapol_z_parallel', &
          'Could not allocate memory for t_re', .true.)
      allocate (t_im(tnx,tny), stat=stat)
      if (stat > 0) call fatal_error('field_extrapol_z_parallel', &
          'Could not allocate memory for t_im', .true.)
!
      call cffti(nxgrid, wsavex)
      call cffti(nygrid, wsavey)
!
!  Collect the data we need.
!
      p_re = in
      p_im = 0.0
!
      do m = 1, pny
!
!  Transform x-direction.
!
        ax = cmplx (p_re(:,m), p_im(:,m))
        call cfftf(nxgrid, ax, wsavex)
        p_re(:,m) = real (ax)
        p_im(:,m) = aimag (ax)
      enddo
!
      call transp_pencil_xy(p_re, t_re)
      call transp_pencil_xy(p_im, t_im)
!
      if (allocated (p_re)) deallocate (p_re)
      if (allocated (p_im)) deallocate (p_im)
!
      allocate (e_re(tnx,tny,onz,2), stat=stat)
      if (stat > 0) call fatal_error('field_extrapol_z_parallel', &
          'Could not allocate memory for e_re', .true.)
!
      allocate (e_im(tnx,tny,onz,2), stat=stat)
      if (stat > 0) call fatal_error('field_extrapol_z_parallel', &
          'Could not allocate memory for e_im', .true.)
!
      do l = 1, tny
!
!  Transform y-direction.
!
        ay = cmplx (t_re(:,l), t_im(:,l))
        call cfftf(nygrid, ay, wsavey)
!
!  Transform y-direction back in each z layer.
!
        do pos_z = 1, onz
!
!  Apply factor to fourier coefficients.
!
          ay_extra = ay * factor(:,l,pos_z)
!
!  x-component of A.
!
          ay_extra_x = cmplx (-aimag (ay_extra), real (ay_extra)) * ky_fft
          call cfftb (nygrid, ay_extra_x, wsavey)
          e_re(:,l,pos_z,1) = real (ay_extra_x)
          e_im(:,l,pos_z,1) = aimag (ay_extra_x)
!
!  y-component of A.
!
          ay_extra_y = cmplx (aimag (ay_extra), -real (ay_extra)) * kx_fft(l+(iproc-ipz*nprocxy)*tny)
          call cfftb (nygrid, ay_extra_y, wsavey)
          e_re(:,l,pos_z,2) = real (ay_extra_y)
          e_im(:,l,pos_z,2) = aimag (ay_extra_y)
        enddo
      enddo
!
      if (allocated (t_re)) deallocate (t_re)
      if (allocated (t_im)) deallocate (t_im)
!
      allocate (b_re(pnx,pny,onz,2), stat=stat)
      if (stat > 0) call fatal_error('field_extrapol_z_parallel', &
          'Could not allocate memory for b_re', .true.)
!
      allocate (b_im(pnx,pny,onz,2), stat=stat)
      if (stat > 0) call fatal_error('field_extrapol_z_parallel', &
          'Could not allocate memory for b_im', .true.)
!
      call transp_pencil_xy(e_re, b_re)
      call transp_pencil_xy(e_im, b_im)
!
      if (allocated (e_re)) deallocate (e_re)
      if (allocated (e_im)) deallocate (e_im)
!
!  Transform x-direction back.
!
      do pos_z = 1, onz
        do m = 1, pny
!
!  x-component of A.
!
          ax = cmplx (b_re(:,m,pos_z,1), b_im(:,m,pos_z,1))
          call cfftb (nxgrid, ax, wsavex)
          b_re(:,m,pos_z,1) = real (ax)
!
!  y-component of A.
!
          ax = cmplx (b_re(:,m,pos_z,2), b_im(:,m,pos_z,2))
          call cfftb (nxgrid, ax, wsavex)
          b_re(:,m,pos_z,2) = real (ax)
        enddo
      enddo
!
!  Distribute the results.
!
      call unmap_from_pencil_xy(b_re, out)
!
      if (allocated (b_re)) deallocate (b_re)
      if (allocated (b_im)) deallocate (b_im)
!
    endsubroutine field_extrapol_z_parallel
!***********************************************************************
    subroutine fourier_transform_y_y(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform of a 1-D array under MPI. Not very
!  efficient since the ipy=0 processors do all the work.
!
!   3-sep-2008/anders: adapted from fourier_transform_xy_xy
!
      use Mpicomm, only: mpirecv_real, mpisend_real
!
      real, dimension (ny) :: a_re, a_im
      logical, optional :: linv
!
      real, dimension (nygrid) :: a_re_full, a_im_full
      complex, dimension (nygrid) :: ay
      real, dimension (4*nygrid+15) :: wsavey
      integer :: ipy_send
      integer, parameter :: itag1=100, itag2=200
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      if (lforward) then
!
        if (nygrid>1) then
!
!  Transform y-direction.
!
          if (lfirst_proc_y) then
            a_re_full(1:ny)=a_re
            a_im_full(1:ny)=a_im
            do ipy_send=1,nprocy-1
              call mpirecv_real(a_re_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,itag1)
              call mpirecv_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,itag2)
            enddo
          else
            call mpisend_real(a_re,ny,iproc-ipy,itag1)
            call mpisend_real(a_im,ny,iproc-ipy,itag2)
          endif
!
          if (lfirst_proc_y) then
!
            call cffti(nygrid,wsavey)
!
            ay=cmplx(a_re_full,a_im_full)
            call cfftf(nygrid,ay,wsavey)
            a_re_full=real(ay)
            a_im_full=aimag(ay)
!
            a_re=a_re_full(1:ny)
            a_im=a_im_full(1:ny)
!
            do ipy_send=1,nprocy-1
              call mpisend_real(a_re_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,itag1)
              call mpisend_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,itag2)
            enddo
          else
            call mpirecv_real(a_re,ny,iproc-ipy,itag1)
            call mpirecv_real(a_im,ny,iproc-ipy,itag2)
          endif
!
        endif
!
      else
!
!  Transform y-direction back.
!
        if (nygrid>1) then
          if (lfirst_proc_y) then
            a_re_full(1:ny)=a_re
            a_im_full(1:ny)=a_im
            do ipy_send=1,nprocy-1
              call mpirecv_real(a_re_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,itag1)
              call mpirecv_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,itag2)
            enddo
          else
            call mpisend_real(a_re,ny,iproc-ipy,itag1)
            call mpisend_real(a_im,ny,iproc-ipy,itag2)
          endif
!
          if (lfirst_proc_y) then
!
            call cffti(nygrid,wsavey)
!
            ay=cmplx(a_re_full,a_im_full)
            call cfftb(nygrid,ay,wsavey)
            a_re_full=real(ay)
            a_im_full=aimag(ay)
!
            a_re=a_re_full(1:ny)
            a_im=a_im_full(1:ny)
!
            do ipy_send=1,nprocy-1
              call mpisend_real(a_re_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,itag1)
              call mpisend_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,itag2)
            enddo
          else
            call mpirecv_real(a_re,ny,iproc-ipy,itag1)
            call mpirecv_real(a_im,ny,iproc-ipy,itag2)
          endif
!
        endif
!
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/nygrid
        a_im=a_im/nygrid
      endif
!
    endsubroutine fourier_transform_y_y
!***********************************************************************
    subroutine fourier_shift_yz_y(a_re,shift_y)
!
!  Performs a periodic shift in the y-direction of an entire y-z plane by
!  the amount shift_y. The shift is done in Fourier space for maximum
!  interpolation accuracy.
!
!  19-jul-06/anders: coded
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      real, dimension (ny,nz) :: a_re
      real :: shift_y
!
      complex, dimension (nygrid) :: a_cmplx, cmplx_shift
      real, dimension (nygrid,max(nz/nprocy,1)) :: a_re_new, a_im_new
      integer :: n, nz_new, ipy_from, ipy_to, iproc_from, iproc_to
      integer :: nprocy_used
      integer, parameter :: itag=666
!
!  Fourier transform of the subdivided y-interval is done by collecting
!  pencils of length nygrid at each processor. Consider the processor space
!  for a given ipz for nygrid=32, nz=8:
!
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    |               |               |               |               | 8
!    |               |               |               |               | 7
!    |               |               |               |               | 6
!    |               |               |               |               | 5
!  z |     ipy=0     |     ipy=1     |     ipy=2     |     ipy=3     | 4
!    |               |               |               |               | 3
!    |               |               |               |               | 2
!    |               |               |               |               | 1
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                                    y
!
!  This is the resulting processor division that can be used for the Fourier
!  transform:
!
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    |                             ipy=3                             | 8
!    |_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _| 7
!    |                             ipy=2                             | 6
!    |_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _| 5
!  z |                             ipy=1                             | 4
!    |_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _| 3
!    |                             ipy=0                             | 2
!    |                                                               | 1
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                                    y
!
!  The height at each processor is nz/nprocy.
!
      nz_new=max(nz/nprocy,1)
      if (nprocy/=1) then
        if (nzgrid==1) then
!
!  Degenerate z-direction. Let y-root processor do the shift of the single
!  nygrid pencil.
!
          do ipy_from=1,nprocy-1
            if (lfirst_proc_y) then
              call mpirecv_real( &
                  a_re_new(ipy_from*ny+1:(ipy_from+1)*ny,1), &
                  ny,ipy_from*nprocx+ipx,itag)
            else
              ipy_to=0
              if (ipy==ipy_from) &
                  call mpisend_real(a_re(:,1),ny,ipy_to*nprocx+ipx,itag)
            endif
          enddo
          if (lfirst_proc_y) a_re_new(1:ny,1)=a_re(:,1)
        else
!
!  Present z-direction. If nz<nprocy we have less y-pencils to shift than
!  we have processors. In that case we give one pencil to the first
!  nprocy_used processors, while the other processors get zero pencils.
!
          if (nz>=nprocy) then
            nprocy_used=nprocy
          else
            nprocy_used=nz
          endif
!
          do ipy_from=0,nprocy-1
            iproc_from=ipz*nprocy*nprocx+ipy_from*nprocx+ipx
            if (ipy/=ipy_from) then
              if (ipy<nprocy_used) call mpirecv_real( &
                  a_re_new(ipy_from*ny+1:(ipy_from+1)*ny,:), &
                  (/ny,nz_new/),iproc_from,itag)
            else
              if (ipy<nprocy_used) a_re_new(ipy*ny+1:(ipy+1)*ny,:) = &
                  a_re(:,ipy*nz_new+1:(ipy+1)*nz_new)
              do ipy_to=0,nprocy_used-1
                iproc_to=ipz*nprocy*nprocx+ipy_to*nprocx+ipx
                if (ipy/=ipy_to) call mpisend_real( &
                    a_re(:,ipy_to*nz_new+1:(ipy_to+1)*nz_new), &
                    (/ny,nz_new/),iproc_to,itag)
              enddo
            endif
          enddo
        endif
      else
!
!  Only parallelization along z (or not at all).
!
        a_re_new(1:ny,1:nz_new)=a_re(1:ny,1:nz_new)
      endif
!
      a_im_new=0.0
      cmplx_shift=exp(cmplx(0.0,-ky_fft*shift_y))
!
!  Transform to Fourier space.
!
      do n=1,nz_new
        call fourier_transform_other(a_re_new(:,n),a_im_new(:,n))
        a_cmplx=cmplx(a_re_new(:,n),a_im_new(:,n))
        a_cmplx=a_cmplx*cmplx_shift
        a_re_new(:,n)=real(a_cmplx)
        a_im_new(:,n)=aimag(a_cmplx)
      enddo
!
!  Back to real space.
!
      do n=1,nz_new
        call fourier_transform_other(a_re_new(:,n),a_im_new(:,n),linv=.true.)
      enddo
!
!  Reinstate original division of yz-plane.
!
      if (nprocy/=1) then
        if (nzgrid==1) then
!  No z-direction.
          if (lfirst_proc_y) then
            do ipy_to=1,nprocy-1
              call mpisend_real( &
                  a_re_new(ipy_to*ny+1:(ipy_to+1)*ny,1), &
                  ny,ipy_to*nprocx+ipx,itag)
            enddo
          else
            call mpirecv_real(a_re(:,1),ny,ipx,itag)
          endif
          if (lfirst_proc_y) a_re(:,1)=a_re_new(1:ny,1)
        else
!  Present z-direction.
          do ipy_from=0,nprocy_used-1
            iproc_from=ipz*nprocy*nprocx+ipy_from*nprocx+ipx
            if (ipy/=ipy_from) then
              call mpirecv_real( &
                  a_re(:,ipy_from*nz_new+1:(ipy_from+1)*nz_new), &
                  (/ny,nz_new/),iproc_from,itag+100)
            else
              if (ipy<nprocy_used) a_re(:,ipy*nz_new+1:(ipy+1)*nz_new)= &
                  a_re_new(ipy*ny+1:(ipy+1)*ny,:)
              do ipy_to=0,nprocy-1
                iproc_to=ipz*nprocy*nprocx+ipy_to*nprocx+ipx
                if (ipy/=ipy_to) call mpisend_real( &
                    a_re_new(ipy_to*ny+1:(ipy_to+1)*ny,:), &
                    (/ny,nz_new/),iproc_to,itag+100)
              enddo
            endif
          enddo
        endif
      else
!  Only parallelization along z (or not at all).
        a_re(1:ny,1:nz_new)=a_re_new(1:ny,1:nz_new)
      endif
!
    endsubroutine fourier_shift_yz_y
!***********************************************************************
    subroutine fourier_shift_y(a_re,shift_y)
!
!  Performs a periodic shift in the y-direction by the amount shift_y(x).
!  The shift is done in Fourier space for maximum interpolation accuracy.
!
!  04-oct-07/anders: adapted from fourier_transform_shear
!
      real, dimension (nx,ny,nz) :: a_re
      real, dimension (nx) :: shift_y
!
      real, dimension (nx,ny,nz) :: a_im
      complex, dimension (nx) :: ay
      real, dimension (4*nx+15) :: wsave
      integer :: l,n,two
!
      two = 2         ! avoid 'array out of bounds' below for nygrid=1
!
!  if nxgrid/=nygrid, then stop.
!
      if (nygrid/=nxgrid .and. nygrid/=1) then
        print*, 'fourier_shift_y: need to have nygrid=nxgrid if nygrid/=1'
        call fatal_error('fourier_transform_shear','')
      endif
!
!  Initialize cfft.
!
      call cffti(nygrid,wsave)
!
!  Transform y-direction.
!
      if (nygrid/=1) then
        a_im=0.0
        if (lroot.and.ip<10) print*, 'fourier_shift_y: doing FFTpack in y'
        call transp(a_re,'y')
        do n=1,nz; do l=1,ny
          ay=cmplx(a_re(:,l,n),a_im(:,l,n))
          call cfftf(nygrid,ay,wsave)
!
!  Shift all modes by the amount shift_y(x).
!
          ay(two:nxgrid)=ay(two:nxgrid)*exp(cmplx(0.0,-ky_fft(two:nxgrid)*shift_y(l+ipy*ny)))
          a_re(:,l,n)=real(ay)
          a_im(:,l,n)=aimag(ay)
        enddo; enddo
!
!  Transform y-direction back.
!
        if (lroot.and.ip<10) print*, 'fourier_shift_y: doing FFTpack in y'
        do n=1,nz; do l=1,ny
          ay=cmplx(a_re(:,l,n),a_im(:,l,n))
          call cfftb(nygrid,ay,wsave)
          a_re(:,l,n)=real(ay)
          a_im(:,l,n)=aimag(ay)
        enddo; enddo
        call transp(a_re,'y')
        call transp(a_im,'y')
!
!  Normalize
!
        a_re=a_re/nygrid
      endif
!
    endsubroutine fourier_shift_y
!***********************************************************************
    subroutine fourier_transform_real_1(a,na,ifirst_fft,wsavex_temp,linv)
!
!  Subroutine to do Fourier transform on a 1-D *real* array of arbitrary size.
!  This routine does not operate in parallel, but should be used to Fourier
!  transform an array present in its entirety on the local processor.
!  The routine overwrites the input data.
!
!  1-jul-2006/dhruba: Adapted from fourier_transform_other_1
!
      integer, intent (in) :: na,ifirst_fft
      real, dimension (na) :: a
      logical, optional :: linv
      real, dimension (2*na+15) :: wsavex_temp
!
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      if (ifirst_fft==1) then
! Initialize fftpack
        call rffti(na,wsavex_temp)
      else
      endif
!
!  Transform x-direction.
!
      if (lforward) then
        if (lroot .and. ip<10) &
            print*, 'fourier_transform_real_1: doing FFTpack in x'
        call rfftf(na,a,wsavex_temp)
      else
!
!  Transform x-direction back.
!
        if (lroot .and. ip<10) &
            print*, 'fourier_transform_real_1: doing FFTpack in x'
        call rfftb(na,a,wsavex_temp)
      endif
!
!  Normalize
!
      if (lforward) then
        a=a/na
      endif
!
      if (lroot .and. ip<10) &
          print*, 'fourier_transform_real_1: fft has finished'
!
    endsubroutine fourier_transform_real_1
!***********************************************************************
endmodule Fourier
