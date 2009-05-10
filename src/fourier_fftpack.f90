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
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif
!
!  need to initialize cfft only once, because nxgrid=nygrid
!
      call cffti(nx,wsavex)
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
      real,dimension(nx,ny,nz)     :: a_re,a_im
      real,dimension(nygrid,ny,nz) :: tmp_re,tmp_im
      real,dimension(nygrid),save  :: xnyg
      real,dimension(4*nygrid+15)  :: wsave
      complex,dimension(nygrid)    :: ay
      real    :: dnx
      integer :: l,n,iarr,ix,ido,iup,i
      logical :: lforward,lnormalize,err,lfirstcall=.true.
      logical, optional :: linv
      logical, optional :: lnorm
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
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif
!
      lnormalize=.true.
      if (present(lnorm)) then
        if (.not.lnorm) lnormalize=.false.
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
      real, dimension (nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (nxgrid) :: ax
      complex, dimension (nxgrid) :: ay
      complex, dimension (nxgrid) :: az
      real, dimension (4*nxgrid+15) :: wsave
      real :: deltay_x
      integer :: l,m,n,two
      logical :: lforward
!
      two = 2         ! avoid `array out of bounds' below for nygrid=1
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
      real, dimension (nx,ny,nz) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension (nxgrid) :: ax
      complex, dimension (nxgrid) :: ay
      complex, dimension (nxgrid) :: az
      real, dimension (4*nxgrid+15) :: wsave
      real :: deltay_x
      integer :: l,m,n,two
      logical :: lforward
!
      two = 2         ! avoid `array out of bounds' below for nygrid=1
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
      use Mpicomm, only: stop_it
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
!  This routine does not operate in parallel, but should be used to Fourier
!  transform an array present in its entirety on the local processor.
!  The routine overwrites the input data.
!
!  28-jul-2006/anders: adapted from fourier_transform_1
!
      use Mpicomm, only: stop_it
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
    subroutine fourier_transform_xy_xy(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform of a 2-D array under MPI.
!  nxgrid is restricted to be an integer multiple of nygrid.
!  The routine overwrites the input data.
!
!   6-oct-2006/tobi: adapted from fourier_transform_other_2
!
      use Mpicomm, only: transp_xy
!
      real, dimension(nx,ny) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension(nxgrid) :: ax
      complex, dimension(nygrid) :: ay
      real, dimension(4*nxgrid+15) :: wsavex
      real, dimension(4*nygrid+15) :: wsavey
      real, dimension(ny) :: deltay_x
      integer :: l,m,ibox
      logical :: lforward
!
      if (mod(nxgrid,nygrid)/=0) then
        call fatal_error('fourier_transform_xy_xy', &
                         'nxgrid needs to be an integer multiple of nygrid.')
      endif
!
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
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
          call transp_xy(a_im)
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
              a_im(iy+1:iy+nygrid,l)=aimag(ay)
            enddo
          enddo
!
          call transp_xy(a_re)
          call transp_xy(a_im)
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
      use Mpicomm, only: transp_xy_other,transp_xy
!
      real, dimension(:,:) :: a_re,a_im
      logical, optional :: linv
!
      complex, dimension(size(a_re,1)) :: ax
      complex, dimension(nprocy*size(a_re,2)) :: ay
      real, dimension(4*size(a_re,1)+15) :: wsavex
      real, dimension(4*nprocy*size(a_re,2)+15) :: wsavey
      integer :: l,m,ibox,nx_other,ny_other
      integer :: nxgrid_other,nygrid_other
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif
!
      nx_other=size(a_re,1); ny_other=size(a_re,2)
      nxgrid_other=nx_other
      nygrid_other=ny_other*nprocy     
!
      if (nxgrid_other/=nygrid_other) then
        call fatal_error('fourier_transform_xy_xy_other', &
             'nxgrid_other needs to be equal to nygrid_other.')
      endif
!
      if (lforward) then
        if (nygrid_other > 1) then
!
!  Transform y-direction.
!
          call transp_xy_other(a_re)
          call transp_xy_other(a_im)

          call cffti(nygrid_other,wsavey)

          do l=1,ny_other
            ay=cmplx(a_re(:,l),a_im(:,l))
            call cfftf(nygrid_other,ay,wsavey)
            a_re(:,l)=real(ay)
            a_im(:,l)=aimag(ay)
          enddo

          call transp_xy_other(a_re)
          call transp_xy_other(a_im)

      endif
      
      if (nxgrid > 1) then
!
!  Transform x-direction
!
        call cffti(nxgrid_other,wsavex)

        do m=1,ny_other
          ax=cmplx(a_re(:,m),a_im(:,m))
          call cfftf(nxgrid_other,ax,wsavex)
          a_re(:,m)=real(ax)
          a_im(:,m)=aimag(ax)
        enddo
        
      endif

    else

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

          call cffti(nygrid_other,wsavey)

          do l=1,ny_other
            ay=cmplx(a_re(:,l),a_im(:,l))
            call cfftb(nygrid_other,ay,wsavey)
            a_re(:,l)=real(ay)
            a_im(:,l)=aimag(ay)
          enddo

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
    subroutine fourier_transform_y_y(a_re,a_im,linv)
!
!  Subroutine to do Fourier transform of a 1-D array under MPI. Not very
!  efficient since the ipy=0 processors do all the work.
!
!   3-sep-2008/anders: adapted from fourier_transform_xy_xy
!
      use Mpicomm, only: mpirecv_real, mpisend_real
!
      real, dimension(ny) :: a_re, a_im
      logical, optional :: linv
!
      real, dimension(nygrid) :: a_re_full, a_im_full
      complex, dimension(nygrid) :: ay
      real, dimension(4*nygrid+15) :: wsavey
      integer :: ipy_send
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
      endif
!
      if (lforward) then
!
        if (nygrid>1) then
!
!  Transform y-direction.
!
          if (ipy==0) then
            a_re_full(1:ny)=a_re
            a_im_full(1:ny)=a_im
            do ipy_send=1,nprocy-1
              call mpirecv_real(a_re_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,111)
              call mpirecv_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,222)
            enddo
          else
            call mpisend_real(a_re,ny,iproc-ipy,111)
            call mpisend_real(a_im,ny,iproc-ipy,222)
          endif
!
          if (ipy==0) then
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
                  ny,iproc+ipy_send,111)
              call mpisend_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,222)
            enddo
          else
            call mpirecv_real(a_re,ny,iproc-ipy,111)
            call mpirecv_real(a_im,ny,iproc-ipy,222)
          endif
!
        endif
!
      else
!
!  Transform y-direction back.
!
        if (nygrid>1) then
          if (ipy==0) then
            a_re_full(1:ny)=a_re
            a_im_full(1:ny)=a_im
            do ipy_send=1,nprocy-1
              call mpirecv_real(a_re_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,111)
              call mpirecv_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,222)
            enddo
          else
            call mpisend_real(a_re,ny,iproc-ipy,111)
            call mpisend_real(a_im,ny,iproc-ipy,222)
          endif
!
          if (ipy==0) then
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
                  ny,iproc+ipy_send,111)
              call mpisend_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,iproc+ipy_send,222)
            enddo
          else
            call mpirecv_real(a_re,ny,iproc-ipy,111)
            call mpirecv_real(a_im,ny,iproc-ipy,222)
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
      use Mpicomm
!
      real, dimension (ny,nz) :: a_re
      real :: shift_y
!
      complex, dimension(nygrid) :: a_cmplx, cmplx_shift
      real, dimension(nygrid,max(nz/nprocy,1)) :: a_re_new, a_im_new
      integer :: n, nz_new, ipy_from, ipy_to, iproc_from, iproc_to
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
!  Degenerate z-direction. Let root processor do the shift of the single nygrid
!  pencil.
!
          do iproc_from=1,ncpus-1
            if (lroot) then
              call mpirecv_real( &
                  a_re_new(iproc_from*ny+1:(iproc_from+1)*ny,1), &
                  ny, iproc_from, 666)
            else
              iproc_to=0
              if (iproc==iproc_from) &
                  call mpisend_real(a_re(:,1), ny, iproc_to, 666)
            endif
          enddo
          if (lroot) a_re_new(1:ny,1)=a_re(:,1)
        else
!
!  Present z-direction. Here nz must be a whole multiple of nprocy (e.g. nz=8,
!  nprocy=4). This constraint could be removed with a bit of work if it should
!  become necessary.
!
          if (modulo(nz,nprocy)/=0) then
            if (lroot) print*, 'fourier_shift_yz_y: nz must be a whole '// &
                'multiple of nprocy!'
            call fatal_error('fourier_shift_yz_y','')
          endif
!
          do ipy_from=0,nprocy-1
            iproc_from=ipz*nprocy+ipy_from
            if (ipy/=ipy_from) then
              call mpirecv_real( &
                  a_re_new(ipy_from*ny+1:(ipy_from+1)*ny,:), &
                  (/ny,nz_new/), iproc_from, 666)
            else
              a_re_new(ipy*ny+1:(ipy+1)*ny,:) = &
                  a_re(:,ipy*nz_new+1:(ipy+1)*nz_new)
              do ipy_to=0,nprocy-1
                iproc_to=ipz*nprocy+ipy_to
                if (ipy/=ipy_to) call mpisend_real( &
                    a_re(:,ipy_to*nz_new+1:(ipy_to+1)*nz_new), &
                    (/ny,nz_new/), iproc_to, 666)
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
          if (.not. lroot) then
            iproc_from=0
            call mpirecv_real(a_re(:,1), ny, iproc_from, 666)
          else
            do iproc_to=1,ncpus-1
              call mpisend_real( &
                  a_re_new(iproc_to*ny+1:(iproc_to+1)*ny,1), &
                  ny, iproc_to, 666)
            enddo
          endif
          if (lroot) a_re(:,1)=a_re_new(1:ny,1)
        else
!  Present z-direction.
          do ipy_from=0,nprocy-1
            iproc_from=ipz*nprocy+ipy_from
            if (ipy/=ipy_from) then
              call mpirecv_real( &
                  a_re(:,ipy_from*nz_new+1:(ipy_from+1)*nz_new), &
                  (/ny,nz_new/), iproc_from, 666)
            else
              a_re(:,ipy*nz_new+1:(ipy+1)*nz_new)= &
                  a_re_new(ipy*ny+1:(ipy+1)*ny,:)
              do ipy_to=0,nprocy-1
                iproc_to=ipz*nprocy+ipy_to
                if (ipy/=ipy_to) call mpisend_real( &
                    a_re_new(ipy_to*ny+1:(ipy_to+1)*ny,:), &
                    (/ny,nz_new/), iproc_to, 666)
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
      complex, dimension (nxgrid) :: ay
      real, dimension (4*nxgrid+15) :: wsave
      integer :: l,m,n,two
!
      two = 2         ! avoid `array out of bounds' below for nygrid=1
!
!  if nxgrid/=nygrid, then stop.
!
      if (nygrid/=nxgrid .and. nygrid /= 1) then
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
      use Mpicomm, only: stop_it
!
      real, dimension(na) :: a
      integer, intent(in) :: na,ifirst_fft
      logical, optional :: linv
      real, dimension(2*na+15) :: wsavex_temp  
!
      integer :: nx_other
      logical :: lforward
!----------
      if (ifirst_fft==1) then
! Initialize fftpack
        call rffti(na,wsavex_temp)
      else
      endif
!---------
      lforward=.true.
      if (present(linv)) then
        if (linv) lforward=.false.
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
