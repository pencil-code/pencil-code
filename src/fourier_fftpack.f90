! $Id$
!
!  This module contains FFT wrapper subroutines.
!
module Fourier
!
  use Cdata
  use Messages
  use Mpicomm, only: transp,transp_other 
  use General, only: ioptest
!$ use OMP_LIB
!
  implicit none
!
  include 'fourier.h'
!
  complex, dimension(nx) :: ax,ay,az
  !$omp threadprivate(ax,ay,az)
  
  real, dimension(:,:), allocatable :: wsavex
  real, dimension(:,:), allocatable :: wsavey
  real, dimension(:,:), allocatable :: wsavez

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
  real, dimension(:,:), allocatable :: p_re_2d, p_im_2d   ! data in pencil shape
  real, dimension(:,:), allocatable :: t_re_2d, t_im_2d   ! data in transposed pencil shape

  real, dimension(:,:,:), allocatable :: p_re_3d, p_im_3d   ! data in pencil shape
  real, dimension(:,:,:), allocatable :: t_re_3d, t_im_3d   ! data in transposed pencil shape

  real, dimension(:,:,:,:), allocatable :: p_re_4d, p_im_4d   ! data in pencil shape
  real, dimension(:,:,:,:), allocatable :: t_re_4d, t_im_4d   ! data in transposed pencil shape
!
  contains
!***********************************************************************
    subroutine initialize_fourier

      include 'fourier_common.h'
!
!  Initializations of module auxiliaries.
!
      if (.not.lreloading) then 
        if (lactive_dimension(1)) allocate(wsavex(4*nxgrid+15,num_helper_threads))
        if (lactive_dimension(2)) allocate(wsavey(4*nygrid+15,num_helper_threads))
        if (lactive_dimension(3)) allocate(wsavez(4*nzgrid+15,num_helper_threads))
      endif

!$omp parallel num_threads(num_helper_threads)
!$      thread_id = omp_get_thread_num()+1
        if (lactive_dimension(1)) call cffti(nxgrid,wsavex(:,thread_id))
        if (lactive_dimension(2)) call cffti(nygrid,wsavey(:,thread_id))
        if (lactive_dimension(3)) call cffti(nzgrid,wsavez(:,thread_id))
!$omp end parallel
!
    endsubroutine initialize_fourier
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
      real, dimension (nx,ny,nz) :: a_re,a_im   ! effectively, nx=nxgrid due to nprocx=1 required
      logical, optional :: linv
!
      integer :: l,m,n
      logical :: lforward
!
      if (nprocx>1) call fatal_error('fourier_transform','must have nprocx=1')
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      if (lforward) then
        if (ip<10) call information('fourier_transform','doing FFTpack in x')
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nxgrid,ax,wsavex(:,thread_id))
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
        !!$omp end parallel
!
!  Transform y-direction.
!
        if (nygrid/=1) then
          if (nygrid/=nxgrid) call fatal_error('fourier_transform','must have nygrid=nxgrid')
          
          call transp(a_re,'y' )
          call transp(a_im,'y' )
!
!  The length of the array in the y-direction is nx.
!
          if (ip<10) call information('fourier_transform','doing FFTpack in y')
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do n=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftf(nxgrid,ax,wsavex(:,thread_id))
            a_re(:,l,n)=real(ax)
            a_im(:,l,n)=aimag(ax)
          enddo; enddo
          !!$omp end parallel
        endif
!
!  Transform z-direction.
!
        if (nzgrid/=1) then
          if (nzgrid/=nxgrid) call fatal_error('fourier_transform','must have nzgrid=nxgrid')
          
          call transp(a_re,'z' )
          call transp(a_im,'z' )
!
!  The length of the array in the z-direction is also nx.
!
          if (ip<10) call information('fourier_transform','doing FFTpack in z')
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do m=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,m),a_im(:,l,m))
            call cfftf(nxgrid,ax,wsavex(:,thread_id))
            a_re(:,l,m)=real(ax)
            a_im(:,l,m)=aimag(ax)
          enddo; enddo
          !!$omp end parallel
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
          if (nzgrid/=nxgrid) call fatal_error('fourier_transform','must have nzgrid=nxgrid')
!
          if (ip<10) call information('fourier_transform','doing FFTpack in z')
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do m=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,m),a_im(:,l,m))
            call cfftb(nxgrid,ax,wsavex(:,thread_id))
            a_re(:,l,m)=real(ax)
            a_im(:,l,m)=aimag(ax)
          enddo; enddo
          !!$omp end parallel
        endif
!
!  Transform y-direction back. Must transpose to go from (z,x,y) to (y,x,z).
!
        if (nygrid/=1) then
          if (nygrid/=nxgrid) call fatal_error('fourier_transform','must have nygrid=nxgrid')
!
          if (nzgrid/=1) then
            call transp(a_re,'z' )
            call transp(a_im,'z' )
          endif
!
          if (ip<10) call information('fourier_transform','doing FFTpack in y')
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do n=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftb(nxgrid,ax,wsavex(:,thread_id))
            a_re(:,l,n)=real(ax)
            a_im(:,l,n)=aimag(ax)
          enddo; enddo
          !!$omp end parallel
        endif
!
!  Transform x-direction back. Transpose to go from (y,x,z) to (x,y,z).
!
        if (ip<10) call information('fourier_transform','doing FFTpack in x')
        
        if (nygrid==1) then
          call transp(a_re,'z' )
          call transp(a_im,'z' )
        else
          call transp(a_re,'y' )
          call transp(a_im,'y' )
        endif
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nxgrid,ax,wsavex(:,thread_id))
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
        !!$omp end parallel
      endif
!
!  Normalize
!
      if (lforward) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp workshare
        a_re=a_re/nwgrid
        a_im=a_im/nwgrid
        !$omp end workshare 
        !!$omp end parallel
      endif
!
      if (ip<10) call information('fourier_transform','fft has finished')
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
!  19-dec-06/anders: adapted from fourier_transform
!
      real, dimension (:,:,:) :: a_re,a_im
      logical, optional :: linv
!
      integer :: l,m,n
      logical :: lforward
!
      if (nprocx>1) call fatal_error('fourier_transform_xy','must have nprocx=1')
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      if (lforward) then
        if (ip<10) call information('fourier_transform_xy','doing FFTpack in x')
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nxgrid,ax,wsavex(:,thread_id))
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
        !!$omp end parallel
!
!  Transform y-direction.
!
        if (nygrid/=1) then
          if (nygrid/=nxgrid) call fatal_error('fourier_transform_xy','must have nygrid=nxgrid')
          
          call transp(a_re,'y' )
          call transp(a_im,'y' )
!
!  The length of the array in the y-direction is nx.
!
          if (ip<10) call information('fourier_transform_xy','doing FFTpack in y')
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do n=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftf(nxgrid,ax,wsavex(:,thread_id))
            a_re(:,l,n)=real(ax)
            a_im(:,l,n)=aimag(ax)
          enddo; enddo
          !!$omp end parallel
        endif
      else
!
!  Transform y-direction back.
!
        if (nygrid/=1) then
          if (nygrid/=nxgrid) call fatal_error('fourier_transform','must have nygrid=nxgrid')
!
          if (ip<10) call information('fourier_transform_xy','doing FFTpack in y')
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do n=1,nz; do l=1,ny
            ax=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftb(nxgrid,ax,wsavex(:,thread_id))
            a_re(:,l,n)=real(ax)
            a_im(:,l,n)=aimag(ax)
          enddo; enddo
          !!$omp end parallel
        endif
!
!  Transform x-direction back. Transpose to go from (y,x,z) to (x,y,z).
!
        if (ip<10) call information('fourier_transform_xy','doing FFTpack in x')
        if (nygrid/=1) then
          call transp(a_re,'y' )
          call transp(a_im,'y' )
        endif
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nxgrid,ax,wsavex(:,thread_id))
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
        !!$omp end parallel
      endif
!
!  Normalize
!
      if (lforward) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp workshare
        a_re=a_re/nxygrid
        a_im=a_im/nxygrid
        !$omp end workshare 
        !!$omp end parallel
      endif
!
      if (ip<10) call information('fourier_transform_xy','fft has finished')
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
      real, dimension (nx,ny,nz) :: a_re,a_im   ! effectively, nx=nxgrid due to nprocx=1 required
      logical, optional :: linv
!
      integer :: l,m,n
!
      if (present(linv)) then
        if (linv) call not_implemented('fourier_transform_xz','for backwards transform')
      endif
!
      if (nprocx>1) call fatal_error('fourier_transform_xz','must have nprocx=1')
!
!  check whether nxgrid=nygrid=nzgrid
!
      if (nxgrid/=nygrid .or. (nxgrid/=nzgrid .and. nzgrid/=1)) &
        call fatal_error('fourier_transform_xz','must have nxgrid=nygrid=nzgrid')
!
      if (ip<10) call information('fourier_transform_xz','doing FFTpack in x')
      !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
      !$omp do collapse(2)
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        call cfftf(nxgrid,ax,wsavex(:,thread_id))
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
      !!$omp end parallel
      
      call transp(a_re,'z')
      call transp(a_im,'z')
!
!  The length of the array in the z-direction is also nx. Normalization is included.
!
      if (ip<10) call information('fourier_transform_xz','doing FFTpack in z')

      !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
      !$omp do collapse(2)
      do l=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,l),a_im(:,m,l))
        call cfftf(nxgrid,ax,wsavex(:,thread_id))
        a_re(:,m,n)=real(ax)/nwgrid     ! normalize
        a_im(:,m,n)=aimag(ax)/nwgrid
      enddo; enddo
      !!$omp end parallel
!
      if (ip<10) call information('fourier_transform_xz','fft has finished')
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
      real, dimension (nx,ny,nz) :: a_re,a_im   ! effectively, nx=nxgrid due to nprocx=1 required
      logical, optional :: linv
!
      integer :: m,n
      logical :: lforward
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      if (nprocx>1) call fatal_error('fourier_transform_x','must have nprocx=1')
!
!  Check whether nxgrid=nygrid=nzgrid.
!
      if ( (nygrid/=1.and.nygrid/=nxgrid) .or. &
           (nzgrid/=1.and.nzgrid/=nxgrid) ) &
        call fatal_error('fourier_transform_x','must have nxgrid=nygrid=nzgrid')
!
      if (lforward) then
!
!  Transform x-direction to fourier space. Normalization is included.
!
        if (ip<10) call information('fourier_transform_x','doing FFTpack in x')
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nxgrid,ax,wsavex(:,thread_id))
          a_re(:,m,n)=real(ax)/nxgrid
          a_im(:,m,n)=aimag(ax)/nxgrid
        enddo; enddo
        !!$omp end parallel
!
      else
!
!  Transform x-direction back to real space
!
        if (ip<10) call information('fourier_transform_x','doing FFTpack in x')

        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nxgrid,ax,wsavex(:,thread_id))
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
        !!$omp end parallel
!
      endif
!
      if (ip<10) call information('fourier_transform_x','fft has finished')
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
      real    :: dnx
      integer :: l,n,m,iarr,ix,ido,iup,i
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
      if (nprocx>1) call fatal_error('fourier_transform_y','must have nprocx=1')
!
! Separate the problem in two cases. nxgrid>= nygrid and its
! opposite
!
      if (nxgrid >= nygrid) then
        if (mod(nxgrid,nygrid)/=0) call fatal_error('fourier_transform_y', &
            'when nxgrid>= nygrid, nxgrid needs to be an integer multiple of nygrid')
      endif
!
!  initialize cfft (coefficients for fft?)
!
      if (ip<10) call information('fourier_transform_y','doing FFTpack in y')
!
!  Transform differently according to sizes of x and y
!
      if (nxgrid>=nygrid) then
!
!  Transpose, transform and transpose back
!
        if (ip<10) call information('fourier_transform_y','nxgrid>=nygrid')
!
        call transp(a_re,'y' ) ; call transp(a_im,'y' )
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do n=1,nz; do l=1,ny
!  Divide a_re into arrays of size nygrid to fit ay
          do iarr=0,nxgrid/nygrid-1
            ix=iarr*nygrid ; ido=ix+1 ; iup=ix+nygrid
            ay=cmplx(a_re(ido:iup,l,n),a_im(ido:iup,l,n))
            if (lforward) then
              call cfftf(nygrid,ay,wsavey(:,thread_id))
            else
              call cfftb(nygrid,ay,wsavey(:,thread_id))
            endif
            a_re(ido:iup,l,n)=real(ay)
            a_im(ido:iup,l,n)=aimag(ay)
          enddo
        enddo;enddo
        !!$omp end parallel

        call transp(a_re,'y' ) ; call transp(a_im,'y' )
!
! Normalize if forward
!
        if (lforward) then
          if (lnormalize) then
            !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
            !$omp workshare
            a_re=a_re/nygrid
            a_im=a_im/nygrid
            !$omp end workshare 
            !!$omp end parallel
          endif
        endif
!
      else !case nxgrid<nygrid
!
        if (ip<10) call information('fourier_transform_y','nxgrid<nygrid')
!
! Save interpolated values of x to dimension nygrid
!
        if (lfirstcall) then
          dnx=Lxyz(1)/(nygrid-1)
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do
          do i=0,nygrid-1
            xnyg(i+1)=x(l1)+i*dnx
          enddo
          !!$omp end parallel
          lfirstcall=.false.
        endif
!
! Interpolate nx to the same dimension as nygrid, so we
! can transpose
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz;do m=1,ny
          call spline(x(l1:l2),a_re(:,m,n),xnyg,tmp_re(:,m,n),nx,nygrid,err)
          call spline(x(l1:l2),a_im(:,m,n),xnyg,tmp_im(:,m,n),nx,nygrid,err)
        enddo;enddo
        !!$omp end parallel
!
! Transpose, transform, transpose back
!
        call transp_other(tmp_re,'y' ) ; call transp_other(tmp_im,'y' )
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz;do l=1,ny
          ay=cmplx(tmp_re(:,l,n),tmp_im(:,l,n))
          if (lforward) then
            call cfftf(nygrid,ay,wsavey(:,thread_id))
          else
            call cfftb(nygrid,ay,wsavey(:,thread_id))
          endif
          tmp_re(:,l,n)=real(ay)
          tmp_im(:,l,n)=aimag(ay)
        enddo; enddo
        !!$omp end parallel
        call transp_other(tmp_re,'y' ) ; call transp_other(tmp_im,'y' )
!
! Normalize if forward
!
        if (lforward) then
          if (lnormalize) then
            !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
            !$omp workshare
            tmp_re=tmp_re/nygrid
            tmp_im=tmp_im/nygrid
            !$omp end workshare 
            !!$omp end parallel
          endif
        endif
!
! Interpolate (coarsen) back to dimension nx
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz;do m=1,ny
          call spline(xnyg,tmp_re(:,m,n),x(l1:l2),a_re(:,m,n),nygrid,nx,err)
          call spline(xnyg,tmp_im(:,m,n),x(l1:l2),a_im(:,m,n),nygrid,nx,err)
        enddo;enddo
        !!$omp end parallel
      endif
!
      if (ip<10) call information('fourier_transform_y','fft has finished')
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
      real, dimension (nx,ny,nz) :: a_re,a_im   ! effectively, nx=nxgrid due to nprocx=1 required
      logical, optional :: linv
!
      complex, dimension (nxgrid) :: ay
      complex, dimension (nxgrid) :: az
      real :: deltay_x
      integer :: l,m,n
      logical :: lforward
      integer, parameter :: two=2  ! avoid 'array out of bounds' below for nygrid=1
!
      lforward=.true.
      if (present(linv)) lforward=.not.linv
!
      if (nprocx>1) call fatal_error('fourier_transform_shear','must have nprocx=1')
!
!  If nxgrid/=nygrid/=nzgrid, stop.
!
      if (nygrid/=nxgrid .and. nygrid/=1) &
        call fatal_error('fourier_transform_shear','need to have nygrid=nxgrid if nygrid/=1')
      if (nzgrid/=nxgrid .and. nzgrid/=1) &
        call fatal_error('fourier_transform_shear','need to have nzgrid=nxgrid if nzgrid/=1')
!
      if (lforward) then
!
!  Transform y-direction. Must start with y, because x is not periodic (yet).
!
        if (nygrid/=1) then
          if (ip<10) call information('fourier_transform_shear','doing FFTpack in y')
          
          call transp(a_re,'y' )
          call transp(a_im,'y' )
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do n=1,nz; do l=1,ny
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftf(nxgrid,ay,wsavex(:,thread_id))
!  Shift y-coordinate so that x-direction is periodic. This is best done in
!  k-space, by multiplying the complex amplitude by exp[i*ky*deltay(x)].
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            ay(two:nxgrid)=ay(two:nxgrid)*exp(cmplx(0.0, ky_fft(two:nxgrid)*deltay_x))
            a_re(:,l,n)=real(ay)
            a_im(:,l,n)=aimag(ay)
          enddo; enddo
          !!$omp end parallel
        endif
!
!  Transform x-direction and normalize.
!
        if (ip<10) call information('fourier_transform_shear','doing FFTpack in x')
        if (nygrid/=1) then
          call transp(a_re,'y' )
          call transp(a_im,'y' )
        endif
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nxgrid,ax,wsavex(:,thread_id))
          a_re(:,m,n)=real(ax)/nwgrid
          a_im(:,m,n)=aimag(ax)/nwgrid
        enddo; enddo
        !!$omp end parallel
!
!  Transform z-direction.
!
        if (nzgrid/=1) then
          if (ip<10) call information('fourier_transform_shear','doing FFTpack in z')
          
          call transp(a_re,'z' )
          call transp(a_im,'z' )
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do l=1,nz; do m=1,ny
            az=cmplx(a_re(:,m,l),a_im(:,m,l))
            call cfftf(nxgrid,az,wsavex(:,thread_id))
            a_re(:,m,l)=real(az)
            a_im(:,m,l)=aimag(az)
          enddo; enddo
          !!$omp end parallel
        endif
      else
!
!  Transform z-direction back.
!
        if (nzgrid/=1) then
          if (ip<10) call information('fourier_transform_shear','doing FFTpack in z')
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do l=1,nz; do m=1,ny
            az=cmplx(a_re(:,m,l),a_im(:,m,l))
            call cfftb(nxgrid,az,wsavex(:,thread_id))
            a_re(:,m,l)=real(az)
            a_im(:,m,l)=aimag(az)
          enddo; enddo
          !!$omp end parallel
        endif
!
!  Transform x-direction back.
!
        if (ip<10) call information('fourier_transform_shear','doing FFTpack in x')
        if (nzgrid/=1) then
          
          call transp(a_re,'z' )
          call transp(a_im,'z' )
        endif
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nxgrid,ax,wsavex(:,thread_id))
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
        !!$omp end parallel
!
!  Transform y-direction back. Would be more practical to end with the
!  x-direction, but the transformation from y' to y means that we must transform
!  the y-direction last.
!
        if (nygrid/=1) then
          
          call transp(a_re,'y' )
          call transp(a_im,'y' )
          if (ip<10) call information('fourier_transform_shear','doing FFTpack in y')
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do n=1,nz; do l=1,ny
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
!  Shift y-coordinate back to regular frame (see above).
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            ay(two:nxgrid)=ay(two:nxgrid)*exp(cmplx(0.0,-ky_fft(two:nxgrid)*deltay_x))
            call cfftb(nxgrid,ay,wsavex(:,thread_id))
            a_re(:,l,n)=real(ay)
            a_im(:,l,n)=aimag(ay)
          enddo; enddo
          !!$omp end parallel
          call transp(a_re,'y' )  ! Deliver array back in (x,y,z) order.
          call transp(a_im,'y' )
        endif
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
      real, dimension (nx,ny,nz) :: a_re,a_im   ! effectively, nx=nxgrid due to nprocx=1 required
      logical, optional :: linv
!
      complex, dimension (nxgrid) :: ay
      real :: deltay_x
      integer :: l,m,n,two
      logical :: lforward
!
      lforward=.true.
      ! avoid 'array out of bounds' below for nygrid=1
      two = 2
      if (present(linv)) lforward=.not.linv
!
      if (nprocx>1) call fatal_error('fourier_transform_shear_xy','must have nprocx=1')
!
!  If nxgrid/=nygrid/=nzgrid, stop.
!
      if (nygrid/=nxgrid .and. nygrid/=1) &
        call fatal_error('fourier_transform_shear_xy','need to have nygrid=nxgrid if nygrid/=1')
!
      if (lforward) then
!
!  Transform y-direction. Must start with y, because x is not periodic (yet).
!
        if (nygrid>1) then
          if (ip<10) call information('fourier_transform_shear','doing FFTpack in y')
          
          call transp(a_re,'y' )
          call transp(a_im,'y' )
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do n=1,nz; do l=1,ny
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftf(nxgrid,ay,wsavex(:,thread_id))
!  Shift y-coordinate so that x-direction is periodic. This is best done in
!  k-space, by multiplying the complex amplitude by exp[i*ky*deltay(x)].
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            !This seems strange, breaks if two is a parameter 
            ay(two:nx)=ay(two:nx)*exp(cmplx(0.0, ky_fft(two:nxgrid)*deltay_x))
            a_re(:,l,n)=real(ay)
            a_im(:,l,n)=aimag(ay)
          enddo; enddo
          !!$omp end parallel
        endif
!
!  Transform x-direction. Normalization is included.
!
        if (ip<10) call information('fourier_transform_shear','doing FFTpack in x')
        
        if (nygrid/=1) then
          call transp(a_re,'y' )
          call transp(a_im,'y' )
        endif
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nxgrid,ax,wsavex(:,thread_id))
          a_re(:,m,n)=real(ax)/nxygrid
          a_im(:,m,n)=aimag(ax)/nxygrid
        enddo; enddo
        !!$omp end parallel
      else
!
!  Transform x-direction back.
!
        if (ip<10) call information('fourier_transform_shear','doing FFTpack in x')
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nxgrid,ax,wsavex(:,thread_id))
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
        !!$omp end parallel
!
!  Transform y-direction back. Would be more practical to end with the
!  x-direction, but the transformation from y' to y means that we must transform
!  the y-direction last.
!
        if (nygrid>1) then
          
          call transp(a_re,'y' )
          call transp(a_im,'y' )
          if (ip<10) call information('fourier_transform_shear','doing FFTpack in y')
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do n=1,nz; do l=1,ny
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
!  Shift y-coordinate back to regular frame (see above).
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            ay(two:nxgrid)=ay(two:nxgrid)*exp(cmplx(0.0,-ky_fft(two:nxgrid)*deltay_x))
            call cfftb(nxgrid,ay,wsavex(:,thread_id))
            a_re(:,l,n)=real(ay)
            a_im(:,l,n)=aimag(ay)
          enddo; enddo
          !!$omp end parallel
          call transp(a_re,'y' )  ! Deliver array back in (x,y,z) order.
          call transp(a_im,'y' )
        endif
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
        if (ip<10) call information('fourier_transform_other_1','doing FFTpack in x')
        ax=cmplx(a_re,a_im)
        call cfftf(nx_other,ax,wsavex)
        a_re=real(ax)/nx_other
        a_im=aimag(ax)/nx_other
      else
!
!  Transform x-direction back.
!
        if (ip<10) call information('fourier_transform_other_1','doing FFTpack in x')
        ax=cmplx(a_re,a_im)
        call cfftb(nx_other,ax,wsavex)
        a_re=real(ax)
        a_im=aimag(ax)
      endif
!
      if (ip<10) call information('fourier_transform_other_1','fft has finished')
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
      call cffti(nx_other,wsavex)
      call cffti(ny_other,wsavey)
!
      if (lforward) then
!
!  Transform x-direction.
!
        if (ip<10) call information('fourier_transform_other_2','doing FFTpack in x')
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do m=1,ny_other
          ax=cmplx(a_re(:,m),a_im(:,m))
          call cfftf(nx_other,ax,wsavex)
          a_re(:,m)=real(ax)
          a_im(:,m)=aimag(ax)
        enddo
        !!$omp end parallel
!
!  Transform y-direction.
!
        if (ip<10) call information('fourier_transform_other_2','doing FFTpack in y')
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do l=1,nx_other
          ay=cmplx(a_re(l,:),a_im(l,:))
          call cfftf(ny_other,ay,wsavey)
          a_re(l,:)=real(ay)/(nx_other*ny_other)
          a_im(l,:)=aimag(ay)/(nx_other*ny_other)
        enddo
        !!$omp end parallel
      else
!
!  Transform x-direction back.
!
        if (ip<10) call information('fourier_transform_other_2','doing FFTpack in x')
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do m=1,ny_other
          ax=cmplx(a_re(:,m),a_im(:,m))
          call cfftb(nx_other,ax,wsavex)
          a_re(:,m)=real(ax)
          a_im(:,m)=aimag(ax)
        enddo
        !!$omp end parallel
!
!  Transform y-direction back.
!
        if (ip<10) call information('fourier_transform_other_2','doing FFTpack in y')
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do l=1,nx_other
          ay=cmplx(a_re(l,:),a_im(l,:))
          call cfftb(ny_other,ay,wsavey)
          a_re(l,:)=real(ay)
          a_im(l,:)=aimag(ay)
        enddo
        !!$omp end parallel
      endif
!
      if (ip<10) call information('fourier_transform_other_2','fft has finished')
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
      real, dimension (:,:), intent(inout) :: a_re,a_im
      logical, optional, intent(in) :: linv,lneed_im
!
      complex, dimension (size(a_re,1)) :: ax
      real, dimension (size(a_re,2)) :: deltay_x
      integer :: l,m,ibox,nyl
      logical :: lforward,lcompute_im
!
      nyl=size(a_re,2)
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
!  From here on, nx=nxgrid can be assumed.
!
      if (lshear) deltay_x=-deltay*(x(m1+ipy*nyl:m2+ipy*nyl)-(x0+Lx/2))/Lx
!
      if (lforward) then
!
        if (nygrid>1) then
!
!  Transform y-direction.
!
          call transp_xy(a_re )
          if (lcompute_im) then
            call transp_xy(a_im )
          else
            !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
            !$omp workshare
            a_im=0.0
            !$omp end workshare
            !!$omp end parallel
          endif
!
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do ibox=0,nxgrid/nygrid-1
            do l=1,nyl
              iy=ibox*nygrid
              ay=cmplx(a_re(iy+1:iy+nygrid,l),a_im(iy+1:iy+nygrid,l))
              call cfftf(nygrid,ay,wsavey(:,thread_id))
              if (lshear) ay = ay*exp(cmplx(0.,+ky_fft*deltay_x(l)))
              a_re(iy+1:iy+nygrid,l)=real(ay)
              a_im(iy+1:iy+nygrid,l)=aimag(ay)
            enddo
          enddo
          !!$omp end parallel
!
          call transp_xy(a_re )
          call transp_xy(a_im )
!
        endif
!
        if (nxgrid>1) then
!
!  Transform x-direction.
!
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do
          do m=1,nyl
            ax=cmplx(a_re(:,m),a_im(:,m))
            call cfftf(nxgrid,ax,wsavex(:,thread_id))    ! requires ax to have dimension nxgrid, but see declaration
            a_re(:,m)=real(ax)
            a_im(:,m)=aimag(ax)
          enddo
          !!$omp end parallel
!
        endif
!
      else
!
        if (nxgrid>1) then
!
!  Transform x-direction back.
!
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do
          do m=1,nyl
            ax=cmplx(a_re(:,m),a_im(:,m))
            call cfftb(nxgrid,ax,wsavex(:,thread_id))    ! requires ax to have dimension nxgrid, but see declaration
            a_re(:,m)=real(ax)
            a_im(:,m)=aimag(ax)
          enddo
          !!$omp end parallel
!
        endif
!
        if (nygrid>1) then
!
!  Transform y-direction back.
!
          call transp_xy(a_re )
          call transp_xy(a_im )
!
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do collapse(2)
          do ibox=0,nxgrid/nygrid-1
            do l=1,nyl
              iy=ibox*nygrid
              ay=cmplx(a_re(iy+1:iy+nygrid,l),a_im(iy+1:iy+nygrid,l))
              if (lshear) ay = ay*exp(cmplx(0.,-ky_fft*deltay_x(l)))
              call cfftb(nygrid,ay,wsavey(:,thread_id))
              a_re(iy+1:iy+nygrid,l)=real(ay)
              if (lcompute_im) a_im(iy+1:iy+nygrid,l)=aimag(ay)
            enddo
          enddo
          !!$omp end parallel
!
          call transp_xy(a_re )
          if (lcompute_im) call transp_xy(a_im )
!
        endif
!
      endif
!
!  Normalize
!
      if (lforward) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp workshare
        a_re=a_re/nxygrid
        a_im=a_im/nxygrid
        !$omp end workshare 
        !!$omp end parallel
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
      if (nxgrid_other/=nygrid_other) call fatal_error('fourier_transform_xy_xy_other', &
             'nxgrid_other needs to be equal to nygrid_other',lfirst_proc_xy)
!
      if (lforward) then
        if (nygrid_other > 1) then
!
!  Transform y-direction.
!
          call transp_xy_other(a_re )
          call transp_xy_other(a_im )
!
          call cffti(nygrid_other,wsavey)
!
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do
          do l=1,ny_other
            ay=cmplx(a_re(:,l),a_im(:,l))
            call cfftf(nygrid_other,ay,wsavey)
            a_re(:,l)=real(ay)
            a_im(:,l)=aimag(ay)
          enddo
          !!$omp end parallel
!
          call transp_xy_other(a_re )
          call transp_xy_other(a_im )
!
        endif
!
        if (nxgrid_other > 1) then
!
!  Transform x-direction
!
          call cffti(nxgrid_other,wsavex)
!
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do
          do m=1,ny_other
            ax=cmplx(a_re(:,m),a_im(:,m))
            call cfftf(nxgrid_other,ax,wsavex)
            a_re(:,m)=real(ax)
            a_im(:,m)=aimag(ax)
          enddo
          !!$omp end parallel
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
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do
          do m=1,ny_other
            ax=cmplx(a_re(:,m),a_im(:,m))
            call cfftb(nxgrid_other,ax,wsavex)
            a_re(:,m)=real(ax)
            a_im(:,m)=aimag(ax)
          enddo
          !!$omp end parallel
!
        endif
!
        if (nygrid_other>1) then
!
!  Transform y-direction back.
!
          call transp_xy_other(a_re )
          call transp_xy_other(a_im )
!
          call cffti(nygrid_other,wsavey)
!
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp do
          do l=1,ny_other
            ay=cmplx(a_re(:,l),a_im(:,l))
            call cfftb(nygrid_other,ay,wsavey)
            a_re(:,l)=real(ay)
            a_im(:,l)=aimag(ay)
          enddo
          !!$omp end parallel
!
          call transp_xy_other(a_re )
          call transp_xy_other(a_im )
!
        endif
!
      endif
!
!  Normalize
!
      if (lforward) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp workshare
        a_re=a_re/(nxgrid_other*nygrid_other)
        a_im=a_im/(nxgrid_other*nygrid_other)
        !$omp end workshare 
        !!$omp end parallel
      endif
!
    endsubroutine fourier_transform_xy_xy_other
!***********************************************************************
    subroutine fft_x_parallel_1D(a_re,a_im,linv,lneed_im,lignore_shear)
!
!  Subroutine to do FFT of distributed 1D data in the x-direction.
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  08-dec-2010/Bourdin.KIS: adapted from fft_xy_parallel_3D
!
      use Mpicomm, only: remap_to_pencil_x, unmap_from_pencil_x
!
      real, dimension (:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
!
      real, dimension(nxgrid) :: p_re, p_im
      logical :: lforward, lcompute_im, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      if (size (a_re, 1) /= nx) &
          call fatal_error ('fft_x_parallel_1D', 'size of input must be nx', lfirst_proc_x)
      if (size (a_re, 1) /= size (a_im, 1)) &
          call fatal_error ('fft_x_parallel_1D', 'size differs for real and imaginary part', lfirst_proc_x)
!
      if (lshear_loc) call not_implemented('fft_x_parallel_1D', 'Shearing', lfirst_proc_x)
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_x (a_re, p_re )
        if (lcompute_im) then
          call remap_to_pencil_x (a_im, p_im )
        else
          p_im = 0.0
        endif
!
        ! Transform x-direction and normalize.
        ax = cmplx (p_re, p_im)
        call cfftf (nxgrid, ax, wsavex(:,thread_id))
        p_re = real (ax)/nxgrid
        p_im = aimag (ax)/nxgrid
!
        ! Unmap the results back to normal shape.
        
        call unmap_from_pencil_x (p_re, a_re)
        call unmap_from_pencil_x (p_im, a_im)
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_x (a_re, p_re )
        call remap_to_pencil_x (a_im, p_im )
!
        ! Transform x-direction back.
        ax = cmplx (p_re, p_im)
        call cfftb (nxgrid, ax, wsavex(:,thread_id))
        p_re = real (ax)
        if (lcompute_im) p_im = aimag (ax)
!
        ! Unmap the results back to normal shape.
        
        call unmap_from_pencil_x (p_re, a_re)
        if (lcompute_im) call unmap_from_pencil_x (p_im, a_im)
!
      endif
!
    endsubroutine fft_x_parallel_1D
!***********************************************************************
    subroutine fft_x_parallel_2D(a_re,a_im,linv,lneed_im, lignore_shear)
!
!  Subroutine to do FFT of distributed 2D data in the x-direction.
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  08-dec-2010/Bourdin.KIS: adapted from fft_xy_parallel_3D
!
      use Mpicomm, only: remap_to_pencil_xy, unmap_from_pencil_xy
!
      real, dimension (nx,ny), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
!
      integer :: m, stat
      logical :: lforward, lcompute_im, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error ('fft_x_parallel_2D', 'nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_x)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error ('fft_x_parallel_2D', 'nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_x)
!
      if (lshear_loc) call not_implemented('fft_x_parallel_2D', 'Shearing', lfirst_proc_x)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_2d(pnx,pny), p_im_2d(pnx,pny),stat=stat)
      if (stat > 0) call fatal_error ('fft_x_parallel_2D', 'Could not allocate p' , .true.)
      !$omp end single
      !$omp barrier
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_2d )
        if (lcompute_im) then
          call remap_to_pencil_xy (a_im, p_im_2d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_2d = 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        ! Transform x-direction and normalize.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do m = 1, pny
          ax = cmplx (p_re_2d(:,m), p_im_2d(:,m))
          call cfftf (nxgrid, ax, wsavex(:,thread_id))
          p_re_2d(:,m) = real (ax)/nxgrid
          p_im_2d(:,m) = aimag (ax)/nxgrid
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re_2d, a_re )
        call unmap_from_pencil_xy (p_im_2d, a_im )
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_2d )
        call remap_to_pencil_xy (a_im, p_im_2d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do m = 1, pny
          ! Transform x-direction back.
          ax = cmplx (p_re_2d(:,m), p_im_2d(:,m))
          call cfftb (nxgrid, ax, wsavex(:,thread_id))
          p_re_2d(:,m) = real (ax)
          if (lcompute_im) p_im_2d(:,m) = aimag (ax)
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re_2d, a_re )
        if (lcompute_im) call unmap_from_pencil_xy (p_im_2d, a_im )
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_2d, p_im_2d)
      !$omp end single
!
    endsubroutine fft_x_parallel_2D
!***********************************************************************
    subroutine fft_x_parallel_3D(a_re,a_im,linv,lneed_im,lignore_shear)
!
!  Subroutine to do FFT of distributed 3D data in the x-direction.
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  08-dec-2010/Bourdin.KIS: adapted from fft_x_parallel_4D
!
      use Mpicomm, only: remap_to_pencil_xy, unmap_from_pencil_xy
!
      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
      integer :: inz ! size of the third dimension
      integer :: m, stat, pos_z
      logical :: lforward, lcompute_im, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      inz = size (a_re, 3)
!
      if (inz /= size (a_im, 3)) &
          call fatal_error ('fft_x_parallel_3D', 'third dimension differs for real and imaginary part', lfirst_proc_x)
!
      if ((size (a_re, 1) /= nx) .or. (size (a_re, 2) /= ny)) &
          call fatal_error ('fft_x_parallel_3D', 'real array size mismatch /= nx,ny', lfirst_proc_x)
      if ((size (a_im, 1) /= nx) .or. (size (a_im, 2) /= ny)) &
          call fatal_error ('fft_x_parallel_3D', 'imaginary array size mismatch /= nx,ny', lfirst_proc_x)
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error ('fft_x_parallel_3D', 'nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_x)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error ('fft_x_parallel_3D', 'nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_x)
!
      if (lshear_loc) call not_implemented('fft_x_parallel_3D', 'Shearing', lfirst_proc_x)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_3d(pnx,pny,inz), p_im_3d(pnx,pny,inz), stat=stat)
      if (stat > 0) call fatal_error ('fft_x_parallel_3D', 'Could not allocate p', .true.)
      !$omp end single
      !$omp barrier
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_3d )
        if (lcompute_im) then
          call remap_to_pencil_xy (a_im, p_im_3d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_3d = 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do pos_z = 1, inz
          do m = 1, pny
            ! Transform x-direction and normalize.
            ax = cmplx (p_re_3d(:,m,pos_z), p_im_3d(:,m,pos_z))
            call cfftf (nxgrid, ax, wsavex(:,thread_id))
            p_re_3d(:,m,pos_z) = real (ax)/nxgrid
            p_im_3d(:,m,pos_z) = aimag (ax)/nxgrid
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re_3d, a_re )
        call unmap_from_pencil_xy (p_im_3d, a_im )
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_3d )
        call remap_to_pencil_xy (a_im, p_im_3d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do pos_z = 1, inz
          do m = 1, pny
            ! Transform x-direction back.
            ax = cmplx (p_re_3d(:,m,pos_z), p_im_3d(:,m,pos_z))
            call cfftb (nxgrid, ax, wsavex(:,thread_id))
            p_re_3d(:,m,pos_z) = real (ax)
            if (lcompute_im) p_im_3d(:,m,pos_z) = aimag (ax)
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re_3d, a_re )
        if (lcompute_im) call unmap_from_pencil_xy (p_im_3d, a_im )
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_3d, p_im_3d)
      !$omp end single
!
    endsubroutine fft_x_parallel_3D
!***********************************************************************
    subroutine fft_x_parallel_4D(a_re,a_im,linv,lneed_im,lignore_shear)
!
!  Subroutine to do FFT of distributed 4D data in the x-direction.
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  08-dec-2010/Bourdin.KIS: adapted from fft_xy_parallel_4D
!
      use Mpicomm, only: remap_to_pencil_xy, unmap_from_pencil_xy
!
      real, dimension (:,:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
      integer :: inz, ina ! size of the third and fourth dimension
      integer :: m, stat, pos_z, pos_a
      logical :: lforward, lcompute_im, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not.linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      inz = size (a_re, 3)
      ina = size (a_re, 4)
!
      if (inz /= size (a_im, 3)) &
          call fatal_error ('fft_x_parallel_4D','third dimension differs for real and imaginary part', lfirst_proc_x)
      if (ina /= size (a_im, 4)) &
          call fatal_error ('fft_x_parallel_4D','fourth dimension differs for real and imaginary part', lfirst_proc_x)
!
      if ((size (a_re, 1) /= nx) .or. (size (a_re, 2) /= ny)) &
          call fatal_error ('fft_x_parallel_4D','real array size mismatch /= nx,ny', lfirst_proc_x)
      if ((size (a_im, 1) /= nx) .or. (size (a_im, 2) /= ny)) &
          call fatal_error ('fft_x_parallel_4D','imaginary array size mismatch /= nx,ny', lfirst_proc_x)
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error ('fft_x_parallel_4D','nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_x)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error ('fft_x_parallel_4D','nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_x)
!
      if (lshear_loc) call not_implemented('fft_x_parallel_4D', 'Shearing', lfirst_proc_x)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_4d(pnx,pny,inz,ina), p_im_4d(pnx,pny,inz,ina), stat=stat)
      if (stat > 0) call fatal_error ('fft_x_parallel_4D', 'Could not allocate p', .true.)
      !$omp end single
      !$omp barrier
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_4d )
        if (lcompute_im) then
          call remap_to_pencil_xy (a_im, p_im_4d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_4d = 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do pos_z = 1, inz
            do m = 1, pny
              ! Transform x-direction and normalize.
              ax = cmplx (p_re_4d(:,m,pos_z,pos_a), p_im_4d(:,m,pos_z,pos_a))
              call cfftf (nxgrid, ax, wsavex(:,thread_id))
              p_re_4d(:,m,pos_z,pos_a) = real (ax)/nxgrid
              p_im_4d(:,m,pos_z,pos_a) = aimag (ax)/nxgrid
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re_4d, a_re )
        call unmap_from_pencil_xy (p_im_4d, a_im )
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_4d )
        call remap_to_pencil_xy (a_im, p_im_4d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do pos_z = 1, inz
            do m = 1, pny
              ! Transform x-direction back.
              ax = cmplx (p_re_4d(:,m,pos_z,pos_a), p_im_4d(:,m,pos_z,pos_a))
              call cfftb (nxgrid, ax, wsavex(:,thread_id))
              p_re_4d(:,m,pos_z,pos_a) = real (ax)
              if (lcompute_im) p_im_4d(:,m,pos_z,pos_a) = aimag (ax)
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re_4d, a_re )
        if (lcompute_im) call unmap_from_pencil_xy (p_im_4d, a_im )
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_4d, p_im_4d)
      !$omp end single
!
    endsubroutine fft_x_parallel_4D
!***********************************************************************
    subroutine fft_y_parallel_1D(a_re,a_im,linv,lneed_im,shift_y,lignore_shear)
!
!  Subroutine to do FFT of distributed 1D data in the y-direction.
!  For y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding y-coordinate.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  08-dec-2010/Bourdin.KIS: adapted from fft_y_parallel_2D
!
      use Mpicomm, only: remap_to_pencil_y, unmap_from_pencil_y
!
      real, dimension (:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      real, optional :: shift_y
!
      real, dimension(nygrid) :: p_re, p_im
      logical :: lforward, lcompute_im, lshift, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_y)) lshift = .true.
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      if (size (a_re, 1) /= ny) &
          call fatal_error ('fft_y_parallel_1D', 'size of input must be ny', lfirst_proc_y)
      if (size (a_re, 1) /= size (a_im, 1)) &
          call fatal_error ('fft_y_parallel_1D', 'size differs for real and imaginary part', lfirst_proc_y)
!
      if (lshear_loc) call not_implemented('fft_y_parallel_1D', 'Shearing for 1D data', lfirst_proc_y)
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_y (a_re, p_re )
        if (lcompute_im) then
          call remap_to_pencil_y (a_im, p_im )
        else
          p_im = 0.0
        endif
!
        ! Transform y-direction and normalize.
        ay = cmplx (p_re, p_im)
        call cfftf (nygrid, ay, wsavey(:,thread_id))
        if (lshift) ay = ay * exp (cmplx (0,-ky_fft * shift_y))
        p_re = real (ay)/nygrid
        p_im = aimag (ay)/nygrid
!
        ! Unmap the results back to normal shape.
        
        call unmap_from_pencil_y (p_re, a_re)
        call unmap_from_pencil_y (p_im, a_im)
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_y (a_re, p_re )
        call remap_to_pencil_y (a_im, p_im )
!
        ! Transform y-direction back.
        ay = cmplx (p_re, p_im)
        call cfftb (nygrid, ay, wsavey(:,thread_id))
        p_re = real (ay)
        p_im = aimag (ay)
!
        ! Unmap the results back to normal shape.
        
        call unmap_from_pencil_y (p_re, a_re)
        if (lcompute_im) call unmap_from_pencil_y (p_im, a_im)
!
      endif
!
    endsubroutine fft_y_parallel_1D
!***********************************************************************
    subroutine fft_y_parallel_2D(a_re,a_im,linv,lneed_im,shift_y,lignore_shear)
!
!  Subroutine to do FFT of distributed 2D data in the y-direction.
!  For y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding y-coordinate.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  08-dec-2010/Bourdin.KIS: adapted from fft_y_parallel_4D
!
      use Mpicomm, only: remap_to_pencil_y, unmap_from_pencil_y
!
      real, dimension (nx,ny), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      real, dimension (nx), optional :: shift_y
!
      real, dimension (nx) :: deltay_x
!
      integer :: l, stat
      logical :: lforward, lcompute_im, lshift, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_y)) lshift = .true.
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_2d(nx,nygrid), p_im_2d(nx,nygrid), stat=stat)
      if (stat > 0) call fatal_error ('fft_y_parallel_2D', 'Could not allocate p', .true.)
      !$omp end single
      !$omp barrier
!
      if (lshear_loc) deltay_x = -deltay * (x(l1:l2) - (x0+Lx/2))/Lx
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_y (a_re, p_re_2d )
        if (lcompute_im) then
          call remap_to_pencil_y (a_im, p_im_2d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_2d = 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        ! Transform y-direction and normalize.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do l = 1, nx
          ay = cmplx (p_re_2d(l,:), p_im_2d(l,:))
          call cfftf (nygrid, ay, wsavey(:,thread_id))
          if (lshear_loc) ay = ay * exp (cmplx (0, ky_fft * deltay_x(l)))
          if (lshift) ay = ay * exp (cmplx (0,-ky_fft * shift_y(l)))
          p_re_2d(l,:) = real (ay)/nygrid
          p_im_2d(l,:) = aimag (ay)/nygrid
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_y (p_re_2d, a_re)
        call unmap_from_pencil_y (p_im_2d, a_im)
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_y (a_re, p_re_2d )
        call remap_to_pencil_y (a_im, p_im_2d )
!
        ! Transform y-direction back.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do l = 1, nx
          ay = cmplx (p_re_2d(l,:), p_im_2d(l,:))
          if (lshear_loc) ay = ay * exp (cmplx (0, -ky_fft*deltay_x(l)))
          call cfftb (nygrid, ay, wsavey(:,thread_id))
          p_re_2d(l,:) = real (ay)
          p_im_2d(l,:) = aimag (ay)
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_y (p_re_2d, a_re)
        if (lcompute_im) call unmap_from_pencil_y (p_im_2d, a_im)
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_2d, p_im_2d)
      !$omp end single
!
    endsubroutine fft_y_parallel_2D
!***********************************************************************
    subroutine fft_y_parallel_3D(a_re,a_im,linv,lneed_im,shift_y,lignore_shear)
!
!  Subroutine to do FFT of distributed 3D data in the y-direction.
!  For y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding y-coordinate.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  08-dec-2010/Bourdin.KIS: adapted from fft_y_parallel_4D
!
      use Mpicomm, only: remap_to_pencil_y, unmap_from_pencil_y
!
      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      real, dimension (nx), optional :: shift_y
!
      integer :: inz ! size of the third dimension
      real, dimension (nx) :: deltay_x
!
      integer :: l, stat, pos_z
      logical :: lforward, lcompute_im, lshift, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_y)) lshift = .true.
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      inz = size (a_re, 3)
!
      if (inz /= size (a_im, 3)) &
          call fatal_error ('fft_y_parallel_3D', 'third dimension differs for real and imaginary part', lfirst_proc_y)
!
      if ((size (a_re, 1) /= nx) .or. (size (a_re, 2) /= ny)) &
          call fatal_error ('fft_y_parallel_3D', 'real array size mismatch /= nx,ny', lfirst_proc_y)
      if ((size (a_im, 1) /= nx) .or. (size (a_im, 2) /= ny)) &
          call fatal_error ('fft_y_parallel_3D', 'imaginary array size mismatch /= nx,ny', lfirst_proc_y)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_3d(nx,nygrid,inz), p_im_3d(nx,nygrid,inz), stat=stat)
      if (stat > 0) call fatal_error ('fft_y_parallel_3D', 'Could not allocate p', .true.)
      !$omp end single
      !$omp barrier
!
      if (lshear_loc) deltay_x = -deltay * (x(l1:l2) - (x0+Lx/2))/Lx
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        call remap_to_pencil_y (a_re, p_re_3d )
        if (lcompute_im) then
          call remap_to_pencil_y (a_im, p_im_3d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_3d = 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        ! Transform y-direction and normalize.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do pos_z = 1, inz
          do l = 1, nx
            ay = cmplx (p_re_3d(l,:,pos_z), p_im_3d(l,:,pos_z))
            call cfftf (nygrid, ay, wsavey(:,thread_id))
            if (lshear_loc) ay = ay * exp (cmplx (0, ky_fft * deltay_x(l)))
            if (lshift) ay = ay * exp (cmplx (0,-ky_fft * shift_y(l)))
            p_re_3d(l,:,pos_z) = real (ay)/nygrid
            p_im_3d(l,:,pos_z) = aimag (ay)/nygrid
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_y (p_re_3d, a_re)
        call unmap_from_pencil_y (p_im_3d, a_im)
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_y (a_re, p_re_3d )
        call remap_to_pencil_y (a_im, p_im_3d )
!
        ! Transform y-direction back.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do pos_z = 1, inz
          do l = 1, nx
            ay = cmplx (p_re_3d(l,:,pos_z), p_im_3d(l,:,pos_z))
            if (lshear_loc) ay = ay * exp (cmplx (0, -ky_fft*deltay_x(l)))
            call cfftb (nygrid, ay, wsavey(:,thread_id))
            p_re_3d(l,:,pos_z) = real (ay)
            p_im_3d(l,:,pos_z) = aimag (ay)
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_y (p_re_3d, a_re)
        if (lcompute_im) call unmap_from_pencil_y (p_im_3d, a_im)
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_3d, p_im_3d)
      !$omp end single
!
    endsubroutine fft_y_parallel_3D
!***********************************************************************
    subroutine fft_y_parallel_4D(a_re,a_im,linv,lneed_im,shift_y,lignore_shear)
!
!  Subroutine to do FFT of distributed 4D data in the y-direction.
!  For y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding y-coordinate.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  08-dec-2010/Bourdin.KIS: adapted from fft_xy_parallel_4D
!
      use Mpicomm, only: remap_to_pencil_y, unmap_from_pencil_y
!
      real, dimension (:,:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      
      real, dimension (nx), optional :: shift_y
!
      integer :: inz, ina ! size of the third and fourth dimension
      real, dimension (nx) :: deltay_x
!
      integer :: l, stat, pos_z, pos_a
      logical :: lforward, lcompute_im, lshift, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not.linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_y)) lshift = .true.
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      inz = size (a_re, 3)
      ina = size (a_re, 4)
!
      if (inz /= size (a_im, 3)) &
          call fatal_error ('fft_y_parallel_4D', 'third dimension differs for real and imaginary part', lfirst_proc_y)
      if (ina /= size (a_im, 4)) &
          call fatal_error ('fft_y_parallel_4D', 'fourth dimension differs for real and imaginary part', lfirst_proc_y)
!
      if ((size (a_re, 1) /= nx) .or. (size (a_re, 2) /= ny)) &
          call fatal_error ('fft_y_parallel_4D', 'real array size mismatch /= nx,ny', lfirst_proc_y)
      if ((size (a_im, 1) /= nx) .or. (size (a_im, 2) /= ny)) &
          call fatal_error ('fft_y_parallel_4D', 'imaginary array size mismatch /= nx,ny', lfirst_proc_y)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_4d(nx,nygrid,inz,ina), p_im_4d(nx,nygrid,inz,ina), stat=stat)
      if (stat > 0) call fatal_error ('fft_y_parallel_4D', 'Could not allocate p', .true.)
      !$omp end single
      !$omp barrier
!
      if (lshear_loc) deltay_x = -deltay * (x(l1:l2) - (x0+Lx/2))/Lx
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_y (a_re, p_re_4d )
        if (lcompute_im) then
          call remap_to_pencil_y (a_im, p_im_4d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_4d= 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        ! Transform y-direction and normalize.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do pos_z = 1, inz
            do l = 1, nx
              ay = cmplx (p_re_4d(l,:,pos_z,pos_a), p_im_4d(l,:,pos_z,pos_a))
              call cfftf (nygrid, ay, wsavey(:,thread_id))
              if (lshear_loc) ay = ay * exp (cmplx (0, ky_fft * deltay_x(l)))
              if (lshift) ay = ay * exp (cmplx (0,-ky_fft * shift_y(l)))
              p_re_4d(l,:,pos_z,pos_a) = real (ay)/nygrid
              p_im_4d(l,:,pos_z,pos_a) = aimag (ay)/nygrid
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_y (p_re_4d, a_re)
        call unmap_from_pencil_y (p_im_4d, a_im)
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        call remap_to_pencil_y (a_re, p_re_4d )
        call remap_to_pencil_y (a_im, p_im_4d )
!
        ! Transform y-direction back.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do pos_z = 1, inz
            do l = 1, nx
              ay = cmplx (p_re_4d(l,:,pos_z,pos_a), p_im_4d(l,:,pos_z,pos_a))
              if (lshear_loc) ay = ay * exp (cmplx (0, -ky_fft*deltay_x(l)))
              call cfftb (nygrid, ay, wsavey(:,thread_id))
              p_re_4d(l,:,pos_z,pos_a) = real (ay)
              p_im_4d(l,:,pos_z,pos_a) = aimag (ay)
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_y (p_re_4d, a_re)
        if (lcompute_im) call unmap_from_pencil_y (p_im_4d, a_im)
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_4d, p_im_4d)
      !$omp end single
!
    endsubroutine fft_y_parallel_4D
!***********************************************************************
    subroutine fft_z_parallel_1D(a_re,a_im,linv,lneed_im,shift_z,lignore_shear)
!
!  Subroutine to do FFT of distributed 1D data in the z-direction.
!  For z-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding z-coordinate.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  13-dec-2010/Bourdin.KIS: adapted from fft_z_parallel_2D
!
      use Mpicomm, only: remap_to_pencil_z, unmap_from_pencil_z
!
      real, dimension (:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      real, optional :: shift_z
!
      real, dimension(nzgrid) :: p_re, p_im
      logical :: lforward, lcompute_im, lshift, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_z)) lshift = .true.
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      if (size (a_re, 1) /= nz) &
          call fatal_error ('fft_z_parallel_1D', 'size of input must be nz!', lfirst_proc_z)
      if (size (a_re, 1) /= size (a_im, 1)) &
          call fatal_error ('fft_z_parallel_1D', 'size differs for real and imaginary part', lfirst_proc_z)
!
      if (lshear_loc) call not_implemented('fft_z_parallel_1D', 'Shearing', lfirst_proc_z)
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_z (a_re, p_re )
        if (lcompute_im) then
          call remap_to_pencil_z (a_im, p_im )
        else
          p_im= 0.0
        endif
!
        ! Transform z-direction and normalize.
        az = cmplx (p_re, p_im)
        call cfftf (nzgrid, az, wsavez(:,thread_id))
        if (lshift) az = az * exp (cmplx (0,-kz_fft * shift_z))
        p_re = real (az)/nzgrid
        p_im = aimag (az)/nzgrid
!
        ! Unmap the results back to normal shape.
        
        call unmap_from_pencil_z (p_re, a_re)
        call unmap_from_pencil_z (p_im, a_im)
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_z (a_re, p_re )
        call remap_to_pencil_z (a_im, p_im )
!
        ! Transform z-direction back.
        az = cmplx (p_re, p_im)
        call cfftb (nzgrid, az, wsavez(:,thread_id))
        p_re = real (az)
        p_im = aimag (az)
!
        ! Unmap the results back to normal shape.
        
        call unmap_from_pencil_z (p_re, a_re)
        if (lcompute_im) call unmap_from_pencil_z (p_im, a_im)
!
      endif
!
    endsubroutine fft_z_parallel_1D
!***********************************************************************
    subroutine fft_z_parallel_2D(a_re,a_im,linv,lneed_im,shift_z,lignore_shear)
!
!  Subroutine to do FFT of distributed 2D data in the z-direction.
!  The z-component is expected to be in the first dimension.
!  For z-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding z-coordinate.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  13-dec-2010/Bourdin.KIS: adapted from fft_z_parallel_4D
!
      use Mpicomm, only: remap_to_pencil_z, unmap_from_pencil_z
!
      real, dimension (:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      real, dimension (:), optional :: shift_z
!
      integer :: inz, ina ! size of the first and second dimension
      integer :: pos_a, stat
      logical :: lforward, lcompute_im, lshift, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_z)) lshift = .true.
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      inz = size (a_re, 1)
      ina = size (a_re, 2)
!
      if (inz /= nz) &
          call fatal_error ('fft_z_parallel_2D', 'first dimension must be the z-coordinate (nz)', lfirst_proc_z)
      if (inz /= size (a_im, 1)) &
          call fatal_error ('fft_z_parallel_2D', 'first dimension differs for real and imaginary part', lfirst_proc_z)
      if (ina /= size (a_im, 2)) &
          call fatal_error ('fft_z_parallel_2D', 'second dimension differs for real and imaginary part', lfirst_proc_z)
!
      if (lshear_loc) call not_implemented('fft_z_parallel_2D', 'Shearing', lfirst_proc_z)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_2d(nzgrid,ina), p_im_2d(nzgrid,ina), stat=stat)
      if (stat > 0) call fatal_error ('fft_z_parallel_2D', 'Could not allocate p', .true.)
      !$omp end single
      !$omp barrier
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_z (a_re, p_re_2d )
        if (lcompute_im) then
          call remap_to_pencil_z (a_im, p_im_2d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_2d= 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        ! Transform z-direction and normalize.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do pos_a = 1, ina
          az = cmplx (p_re_2d(:,pos_a), p_im_2d(:,pos_a))
          call cfftf (nzgrid, az, wsavez(:,thread_id))
          if (lshift) az = az * exp (cmplx (0,-kz_fft * shift_z(pos_a)))
          p_re_2d(:,pos_a) = real (az)/nzgrid
          p_im_2d(:,pos_a) = aimag (az)/nzgrid
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_z (p_re_2d, a_re)
        call unmap_from_pencil_z (p_im_2d, a_im)
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_z (a_re, p_re_2d )
        call remap_to_pencil_z (a_im, p_im_2d )
!
        ! Transform z-direction back.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do pos_a = 1, ina
          az = cmplx (p_re_2d(:,pos_a), p_im_2d(:,pos_a))
          call cfftb (nzgrid, az, wsavez(:,thread_id))
          p_re_2d(:,pos_a) = real (az)
          p_im_2d(:,pos_a) = aimag (az)
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_z (p_re_2d, a_re)
        if (lcompute_im) call unmap_from_pencil_z (p_im_2d, a_im)
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_2d, p_im_2d)
      !$omp end single
!
    endsubroutine fft_z_parallel_2D
!***********************************************************************
    subroutine fft_z_parallel_3D(a_re,a_im,linv,lneed_im,lignore_shear)
!
!  Subroutine to do FFT of distributed 3D data in the z-direction.
!  The z-component is expected to be in the third dimension.
!  For z-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding z-coordinate.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  08-dec-2010/Bourdin.KIS: adapted from fft_z_parallel_4D
!
      use Mpicomm, only: remap_to_pencil_z, unmap_from_pencil_z
!
      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
!
      integer :: inx, iny ! size of the third dimension
      integer :: l, m, stat
      logical :: lforward, lcompute_im, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      inx = size (a_re, 1)
      iny = size (a_re, 2)
!
      if (inx /= size (a_im, 1)) &
          call fatal_error ('fft_z_parallel_3D', 'first dimension differs for real and imaginary part', lfirst_proc_z)
      if (iny /= size (a_im, 2)) &
          call fatal_error ('fft_z_parallel_3D', 'second dimension differs for real and imaginary part', lfirst_proc_z)
      if (size (a_re, 3) /= nz) &
          call fatal_error ('fft_z_parallel_3D', 'third dimension must be the z-coordinate (nz)', lfirst_proc_z)
!
      if (lshear_loc) call not_implemented('fft_z_parallel_3D', 'Shearing', lfirst_proc_z)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_3d(inx,iny,nzgrid), p_im_3d(inx,iny,nzgrid), stat=stat)
      if (stat > 0) call fatal_error ('fft_z_parallel_3D', 'Could not allocate p', .true.)
      !$omp end single
      !$omp barrier
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_z (a_re, p_re_3d )
        if (lcompute_im) then
          call remap_to_pencil_z (a_im, p_im_3d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_3d= 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        ! Transform z-direction and normalize.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do m = 1, iny
          do l = 1, inx
            az = cmplx (p_re_3d(l,m,:), p_im_3d(l,m,:))
            call cfftf (nzgrid, az, wsavez(:,thread_id))
            p_re_3d(l,m,:) = real (az)/nzgrid
            p_im_3d(l,m,:) = aimag (az)/nzgrid
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_z (p_re_3d, a_re)
        call unmap_from_pencil_z (p_im_3d, a_im)
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_z (a_re, p_re_3d )
        call remap_to_pencil_z (a_im, p_im_3d )
!
        ! Transform z-direction back.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do m = 1, iny
          do l = 1, inx
            az = cmplx (p_re_3d(l,m,:), p_im_3d(l,m,:))
            call cfftb (nzgrid, az, wsavez(:,thread_id))
            p_re_3d(l,m,:) = real (az)
            p_im_3d(l,m,:) = aimag (az)
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_z (p_re_3d, a_re)
        if (lcompute_im) call unmap_from_pencil_z (p_im_3d, a_im)
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_3d, p_im_3d)
      !$omp end single
!
    endsubroutine fft_z_parallel_3D
!***********************************************************************
    subroutine fft_z_parallel_4D(a_re,a_im,linv,lneed_im,lignore_shear)
!
!  Subroutine to do FFT of distributed 4D data in the z-direction.
!  The z-component is expected to be in the third dimension.
!  For z-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding z-coordinate.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  08-dec-2010/Bourdin.KIS: adapted from fft_xy_parallel_4D
!
      use Mpicomm, only: remap_to_pencil_z, unmap_from_pencil_z
!
      real, dimension (:,:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
!
      integer :: inx, iny, ina ! size of the third and fourth dimension
      integer :: l, m, stat, pos_a
      logical :: lforward, lcompute_im, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not.linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      inx = size (a_re, 1)
      iny = size (a_re, 2)
      ina = size (a_re, 4)
!
      if (inx /= size (a_im, 1)) &
          call fatal_error ('fft_z_parallel_4D', 'first dimension differs for real and imaginary part', lfirst_proc_z)
      if (iny /= size (a_im, 2)) &
          call fatal_error ('fft_z_parallel_4D', 'second dimension differs for real and imaginary part', lfirst_proc_z)
      if (size (a_re, 3) /= nz) &
          call fatal_error ('fft_z_parallel_4D', 'third dimension must be the z-coordinate (nz)', lfirst_proc_z)
      if (ina /= size (a_im, 4)) &
          call fatal_error ('fft_z_parallel_4D', 'fourth dimension differs for real and imaginary part', lfirst_proc_z)
!
      if (lshear_loc) call not_implemented('fft_z_parallel_3D', 'Shearing', lfirst_proc_z)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_4d(inx,iny,nzgrid,ina), p_im_4d(inx,iny,nzgrid,ina), stat=stat)
      if (stat > 0) call fatal_error ('fft_z_parallel_4D', 'Could not allocate p', .true.)
      !$omp end single
      !$omp barrier
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_z (a_re, p_re_4d )
        if (lcompute_im) then
          call remap_to_pencil_z (a_im, p_im_4d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_4d= 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        ! Transform z-direction and normalize.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do m = 1, iny
            do l = 1, inx
              az = cmplx (p_re_4d(l,m,:,pos_a), p_im_4d(l,m,:,pos_a))
              call cfftf (nzgrid, az, wsavez(:,thread_id))
              p_re_4d(l,m,:,pos_a) = real (az)/nzgrid
              p_im_4d(l,m,:,pos_a) = aimag (az)/nzgrid
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_z (p_re_4d, a_re)
        call unmap_from_pencil_z (p_im_4d, a_im)
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_z (a_re, p_re_4d )
        call remap_to_pencil_z (a_im, p_im_4d )
!
        ! Transform z-direction back.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do m = 1, iny
            do l = 1, inx
              az = cmplx (p_re_4d(l,m,:,pos_a), p_im_4d(l,m,:,pos_a))
              call cfftb (nzgrid, az, wsavez(:,thread_id))
              p_re_4d(l,m,:,pos_a) = real (az)
              p_im_4d(l,m,:,pos_a) = aimag (az)
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_z (p_re_4d, a_re)
        if (lcompute_im) call unmap_from_pencil_z (p_im_4d, a_im)
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_4d, p_im_4d)
      !$omp end single
!
    endsubroutine fft_z_parallel_4D
!***********************************************************************
    subroutine fft_xy_parallel_2D(a_re,a_im,linv,lneed_im,shift_y,lignore_shear)
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
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      
      real, dimension (nxgrid), optional :: shift_y
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
      real, dimension (tny) :: deltay_x
      real, dimension (tny) :: dshift_y
!
      integer :: l, m, stat, x_offset
      logical :: lforward, lcompute_im, lshift, lnoshear, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_y)) lshift = .true.
!
      lnoshear = .false.
      if (present (lignore_shear)) lnoshear = lignore_shear
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      ! Check for degenerate cases.
      if (nxgrid == 1) then
        !$omp single         ! as a_* is shared and no worksharing inside fft_y_parallel_1D
        if (lshift) then
          call fft_y_parallel (a_re(1,:), a_im(1,:), .not. lforward, lcompute_im, shift_y(1), lignore_shear=lnoshear)
        else
          call fft_y_parallel (a_re(1,:), a_im(1,:), .not. lforward, lcompute_im, lignore_shear=lnoshear)
        endif
        !$omp end single
        return
      endif
      if (nygrid == 1) then
        !$omp single         ! as a_* is shared and no worksharing inside fft_y_parallel_1D
        call fft_x_parallel (a_re(:,1), a_im(:,1), .not. lforward, lcompute_im, lignore_shear=lnoshear)
        !$omp end single
        return
      endif
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error('fft_xy_parallel_2D','nxgrid needs to be an integer multiple of nprocx*nprocy',lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error('fft_xy_parallel_2D','nygrid needs to be an integer multiple of nprocx*nprocy',lfirst_proc_xy)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_2d(pnx,pny), p_im_2d(pnx,pny), t_re_2d(tnx,tny), t_im_2d(tnx,tny), stat=stat)
      if (stat > 0) call fatal_error ('fft_xy_parallel_2D', 'Could not allocate p and t', .true.)
      !$omp end single
      !$omp barrier
!
      if (lshear_loc) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        deltay_x = -deltay * (xgrid(x_offset:x_offset+tny-1) - (x0+Lx/2))/Lx
      endif
!
      if (lshift) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        dshift_y = shift_y(x_offset:x_offset+tny-1)
      endif
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_2d )
        if (lcompute_im) then
          call remap_to_pencil_xy (a_im, p_im_2d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_2d= 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        call transp_pencil_xy (p_re_2d, t_re_2d )
        call transp_pencil_xy (p_im_2d, t_im_2d )
!
        ! Transform y-direction.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do l = 1, tny
          ay = cmplx (t_re_2d(:,l), t_im_2d(:,l))
          call cfftf(nygrid, ay, wsavey(:,thread_id))
          if (lshear_loc) ay = ay * exp (cmplx (0, ky_fft * deltay_x(l)))
          if (lshift) ay = ay * exp (cmplx (0,-ky_fft * dshift_y(l)))
          t_re_2d(:,l) = real (ay)
          t_im_2d(:,l) = aimag (ay)
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy (t_re_2d, p_re_2d )
        call transp_pencil_xy (t_im_2d, p_im_2d )
!
        ! Transform x-direction.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do m = 1, pny
          ax = cmplx (p_re_2d(:,m), p_im_2d(:,m))
          call cfftf (nxgrid, ax, wsavex(:,thread_id))
          p_re_2d(:,m) = real (ax)/nxygrid
          p_im_2d(:,m) = aimag (ax)/nxygrid
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re_2d, a_re )
        call unmap_from_pencil_xy (p_im_2d, a_im )
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_2d )
        call remap_to_pencil_xy (a_im, p_im_2d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do m = 1, pny
          ! Transform x-direction back.
          ax = cmplx (p_re_2d(:,m), p_im_2d(:,m))
          call cfftb (nxgrid, ax, wsavex(:,thread_id))
          p_re_2d(:,m) = real (ax)
          p_im_2d(:,m) = aimag (ax)
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy (p_re_2d, t_re_2d )
        call transp_pencil_xy (p_im_2d, t_im_2d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do l = 1, tny
          ! Transform y-direction back.
          ay = cmplx (t_re_2d(:,l), t_im_2d(:,l))
          if (lshear_loc) ay = ay * exp (cmplx (0, -ky_fft*deltay_x(l)))
          call cfftb (nygrid, ay, wsavey(:,thread_id))
          t_re_2d(:,l) = real (ay)
          if (lcompute_im) t_im_2d(:,l) = aimag (ay)
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy (t_re_2d, p_re_2d )
        if (lcompute_im) call transp_pencil_xy (t_im_2d, p_im_2d )
!
        ! Unmap the results back to normal shape.
        
        call unmap_from_pencil_xy (p_re_2d, a_re )
        if (lcompute_im) call unmap_from_pencil_xy (p_im_2d, a_im )
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_2d, p_im_2d, t_re_2d, t_im_2d)
      !$omp end single
!
    endsubroutine fft_xy_parallel_2D
!***********************************************************************
    subroutine fft_xy_parallel_2D_other(a_re,a_im,linv,lignore_shear)
!
!  For x- and/or y-parallelization the calculation will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  You can set lneed_im=false if the imaginary data can be disregarded.
!  Attention: input data will be overwritten.
!
!  09-jul-2013/wlad: adpated from fft_xy_parallel
!
      use Mpicomm, only: remap_to_pencil_xy_2D_other, transp_pencil_xy, unmap_from_pencil_xy_2D_other
!
      real, dimension (:,:) :: a_re, a_im
      logical, optional :: linv, lignore_shear
!
      integer :: pnx,pny,tnx,tny,nx_other,ny_other,nxgrid_other,nygrid_other
!
      complex, dimension (nprocx*size(a_re,1)) :: ax_other
      complex, dimension (nprocy*size(a_re,2)) :: ay_other
      real, dimension (4*nprocx*size(a_re,1)+15) :: wsavex_other
      real, dimension (4*nprocy*size(a_re,2)+15) :: wsavey_other

      integer :: l, m, stat
      logical :: lforward, lnoshear
!
! pencil shaped data sizes
!
      nx_other=size(a_re,1); ny_other=size(a_re,2)
      nxgrid_other=nx_other*nprocx
      nygrid_other=ny_other*nprocy
!
      pnx=nxgrid_other
      pny=nygrid_other/nprocxy
!
! pencil shaped transposed data sizes
!
      tnx=nprocy*size(a_re,2)
      tny=nprocx*size(a_re,1)/nprocxy
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lnoshear = .false.
      if (present (lignore_shear)) lnoshear = lignore_shear
!
      ! Check for degenerate cases.
      if (nxgrid == 1) then
        !$omp single     ! as a_* is shared and no worksharing inside fft_y_parallel_1D
        call fft_y_parallel(a_re(1,:),a_im(1,:),.not.lforward,lignore_shear=lnoshear)
        !$omp end single
        return
      endif
      if (nygrid == 1) then
        !$omp single     ! as a_* is shared and no worksharing inside fft_y_parallel_1D
        call fft_x_parallel(a_re(:,1), a_im(:,1),.not.lforward,lignore_shear=lnoshear)
        !$omp end single
        return
      endif
!
      if (mod (nxgrid_other, nprocxy) /= 0) call fatal_error('fft_xy_parallel_2D_other', &
          'nxgrid_other needs to be an integer multiple of nprocx*nprocy',lfirst_proc_xy)
      if (mod (nygrid_other, nprocxy) /= 0) call fatal_error('fft_xy_parallel_2D_other', &
          'nygrid_other needs to be an integer multiple of nprocx*nprocy',lfirst_proc_xy)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_2d(pnx,pny), p_im_2d(pnx,pny), t_re_2d(tnx,tny), t_im_2d(tnx,tny), stat=stat)
      if (stat > 0) call fatal_error('fft_xy_parallel_2D','Could not allocate p and t', .true.)
      !$omp end single
      !$omp barrier
!
      call cffti(nxgrid_other,wsavex_other)
      call cffti(nygrid_other,wsavey_other)
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_xy_2D_other(a_re,p_re_2d )
        call remap_to_pencil_xy_2D_other(a_im,p_im_2d )
!
        ! Transform x-direction.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do m=1,pny
          ax_other = cmplx(p_re_2d(:,m),p_im_2d(:,m))
          call cfftf(nxgrid_other,ax_other,wsavex_other)
          p_re_2d(:,m) = real(ax_other)
          p_im_2d(:,m) = aimag(ax_other)
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy(p_re_2d,t_re_2d )
        call transp_pencil_xy(p_im_2d,t_im_2d )
!
        ! Transform y-direction and normalize.
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do l=1,tny
          ay_other = cmplx(t_re_2d(:,l), t_im_2d(:,l))
          call cfftf(nygrid_other,ay_other,wsavey_other)
          t_re_2d(:,l) = real(ay_other)/(nxgrid_other*nygrid_other)
          t_im_2d(:,l) = aimag(ay_other)/(nxgrid_other*nygrid_other)
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy(t_re_2d,p_re_2d )
        call transp_pencil_xy(t_im_2d,p_im_2d )
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy_2D_other(p_re_2d,a_re )
        call unmap_from_pencil_xy_2D_other(p_im_2d,a_im )
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_xy_2D_other(a_re, p_re_2d )
        call remap_to_pencil_xy_2D_other(a_im, p_im_2d )
!
        call transp_pencil_xy(p_re_2d, t_re_2d )
        call transp_pencil_xy(p_im_2d, t_im_2d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do l=1,tny
          ! Transform y-direction back.
          ay_other = cmplx(t_re_2d(:,l),t_im_2d(:,l))
          call cfftb(nygrid_other,ay_other,wsavey_other)
          t_re_2d(:,l) = real(ay_other)
          t_im_2d(:,l) = aimag(ay_other)
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy(t_re_2d,p_re_2d )
        call transp_pencil_xy(t_im_2d,p_im_2d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do m=1,pny
          ! Transform x-direction back.
          ax_other = cmplx(p_re_2d(:,m),p_im_2d(:,m))
          call cfftb(nxgrid_other,ax_other,wsavex_other)
          p_re_2d(:,m) = real(ax_other)
          p_im_2d(:,m) = aimag(ax_other)
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy_2D_other(p_re_2d,a_re )
        call unmap_from_pencil_xy_2D_other(p_im_2d,a_im )
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_2d, p_im_2d, t_re_2d, t_im_2d)
      !$omp end single
!
    endsubroutine fft_xy_parallel_2D_other
!***********************************************************************
    subroutine fft_xy_parallel_3D(a_re,a_im,linv,lneed_im,shift_y,lignore_shear)
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

      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      
      real, dimension (nxgrid), optional :: shift_y
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
      
      integer :: inz ! size of the third dimension
      real, dimension (tny) :: deltay_x
      real, dimension (tny) :: dshift_y
!
      integer :: l, m, stat, pos_z, x_offset, i
      logical :: lforward, lcompute_im, lshift, lnoshear, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_y)) lshift = .true.
!
      lnoshear = .false.
      if (present (lignore_shear)) lnoshear = lignore_shear
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      inz = size (a_re, 3)
!     Check for degenerate cases.
!     KG: (27-sep-2023) it seems to me that fft_y_parallel_2D does the FFT along the second axis; here, we want FFT along the first axis, and moreover have to handle the case where nx!=ny!=nz. This seems the simplest way to accomplish that.
      if (nxgrid == 1) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do i=1,inz
          call fft_y_parallel (a_re(1,:,i), a_im(1,:,i), .not. lforward, lcompute_im, lignore_shear=lnoshear)
        enddo
        !!$omp end parallel
        return
      endif
      if (nygrid == 1) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do i=1,inz
          call fft_x_parallel (a_re(:,1,i), a_im(:,1,i), .not. lforward, lcompute_im, lignore_shear=lnoshear)
        enddo
        !!$omp end parallel
        return
      endif
!
      if (inz /= size (a_im, 3)) &
          call fatal_error ('fft_xy_parallel_3D', 'third dimension differs for real and imaginary part', lfirst_proc_xy)
!
      if ((size (a_re, 1) /= nx) .or. (size (a_re, 2) /= ny)) &
          call fatal_error ('fft_xy_parallel_3D', 'real array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (a_im, 1) /= nx) .or. (size (a_im, 2) /= ny)) &
          call fatal_error ('fft_xy_parallel_3D', 'imaginary array size mismatch /= nx,ny', lfirst_proc_xy)
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error ('fft_xy_parallel_3D', 'nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error ('fft_xy_parallel_3D', 'nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_3d(pnx,pny,inz), p_im_3d(pnx,pny,inz), t_re_3d(tnx,tny,inz), t_im_3d(tnx,tny,inz), stat=stat)
      if (stat > 0) call fatal_error ('fft_xy_parallel_3D', 'Could not allocate p and t', .true.)
      !$omp end single
      !$omp barrier
!
      if (lshear_loc) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        deltay_x = -deltay * (xgrid(x_offset:x_offset+tny-1) - (x0+Lx/2))/Lx
      endif
!
      if (lshift) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        dshift_y = shift_y(x_offset:x_offset+tny-1)
      endif
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_3d )
        if (lcompute_im) then
          call remap_to_pencil_xy (a_im, p_im_3d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_3d = 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        call transp_pencil_xy (p_re_3d, t_re_3d )
        call transp_pencil_xy (p_im_3d, t_im_3d )
!§
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do pos_z = 1, inz
          do l = 1, tny
            ! Transform y-direction.
            ay = cmplx (t_re_3d(:,l,pos_z), t_im_3d(:,l,pos_z))
            call cfftf (nygrid, ay, wsavey(:,thread_id))
            if (lshear_loc) ay = ay * exp (cmplx (0, ky_fft * deltay_x(l)))
            if (lshift) ay = ay * exp (cmplx (0,-ky_fft * dshift_y(l)))
            t_re_3d(:,l,pos_z) = real (ay)
            t_im_3d(:,l,pos_z) = aimag (ay)
          enddo
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy (t_re_3d, p_re_3d )
        call transp_pencil_xy (t_im_3d, p_im_3d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do pos_z = 1, inz
          do m = 1, pny
            ! Transform x-direction and normalize.
            ax = cmplx (p_re_3d(:,m,pos_z), p_im_3d(:,m,pos_z))
            call cfftf (nxgrid, ax, wsavex(:,thread_id))
            p_re_3d(:,m,pos_z) = real (ax)/nxygrid
            p_im_3d(:,m,pos_z) = aimag (ax)/nxygrid
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re_3d, a_re )
        call unmap_from_pencil_xy (p_im_3d, a_im )
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_3d )
        call remap_to_pencil_xy (a_im, p_im_3d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do pos_z = 1, inz
          do m = 1, pny
            ! Transform x-direction back.
            ax = cmplx (p_re_3d(:,m,pos_z), p_im_3d(:,m,pos_z))
            call cfftb (nxgrid, ax, wsavex(:,thread_id))
            p_re_3d(:,m,pos_z) = real (ax)
            p_im_3d(:,m,pos_z) = aimag (ax)
          enddo
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy (p_re_3d, t_re_3d )
        call transp_pencil_xy (p_im_3d, t_im_3d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do pos_z = 1, inz
          do l = 1, tny
            ! Transform y-direction back.
            ay = cmplx (t_re_3d(:,l,pos_z), t_im_3d(:,l,pos_z))
            if (lshear_loc) ay = ay * exp (cmplx (0, -ky_fft*deltay_x(l)))
            call cfftb (nygrid, ay, wsavey(:,thread_id))
            t_re_3d(:,l,pos_z) = real (ay)
            if (lcompute_im) t_im_3d(:,l,pos_z) = aimag (ay)
          enddo
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy (t_re_3d, p_re_3d )
        if (lcompute_im) call transp_pencil_xy (t_im_3d, p_im_3d )
!
        ! Unmap the results back to normal shape.
        
        call unmap_from_pencil_xy (p_re_3d, a_re )
        if (lcompute_im) call unmap_from_pencil_xy (p_im_3d, a_im )
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_3d, p_im_3d, t_re_3d, t_im_3d)
      !$omp end single
!
    endsubroutine fft_xy_parallel_3D
!***********************************************************************
    subroutine fft_xy_parallel_4D(a_re,a_im,linv,lneed_im,shift_y,lignore_shear)
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
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      
      real, dimension (nxgrid), optional :: shift_y
!
      integer, parameter :: pnx=nxgrid, pny=nygrid/nprocxy ! pencil shaped data sizes
      integer, parameter :: tnx=nygrid, tny=nxgrid/nprocxy ! pencil shaped transposed data sizes
      integer :: inz, ina ! size of the third and fourth dimension
      real, dimension (tny) :: deltay_x
      real, dimension (tny) :: dshift_y
!
      integer :: l, m, stat, pos_z, pos_a, x_offset, iz, ia
      logical :: lforward, lcompute_im, lshift, lnoshear, lshear_loc
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_y)) lshift = .true.
!
      lnoshear = .false.
      if (present (lignore_shear)) lnoshear = lignore_shear
!
      lshear_loc = lshear
      if (present (lignore_shear)) lshear_loc = (.not.lignore_shear).and.lshear
!
      inz = size (a_re, 3)
      ina = size (a_re, 4)
!     Check for degenerate cases.
!     KG: (27-sep-2023) it seems to me that fft_y_parallel_3D does the FFT along the second axis; here, we want FFT along the first axis, and moreover have to handle the case where nx!=ny!=nz. This seems the simplest way to accomplish that.
      if (nxgrid == 1) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do iz=1,inz
          do ia = 1,ina
            call fft_y_parallel (a_re(1,:,iz,ia), a_im(1,:,iz,ia), .not. lforward, lcompute_im, lignore_shear=lnoshear)
          enddo
        enddo
        !!$omp end parallel
        return
      endif
      if (nygrid == 1) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do iz=1,inz
          do ia = 1,ina
            call fft_x_parallel (a_re(:,1,iz,ia), a_im(:,1,iz,ia), .not. lforward, lcompute_im, lignore_shear=lnoshear)
          enddo
        enddo
        !!$omp end parallel
        return
      endif
!
      if (inz /= size (a_im, 3)) &
          call fatal_error ('fft_xy_parallel_4D', 'third dimension differs for real and imaginary part', lfirst_proc_xy)
      if (ina /= size (a_im, 4)) &
          call fatal_error ('fft_xy_parallel_4D', 'fourth dimension differs for real and imaginary part', lfirst_proc_xy)
!
      if ((size (a_re, 1) /= nx) .or. (size (a_re, 2) /= ny)) &
          call fatal_error ('fft_xy_parallel_4D', 'real array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (a_im, 1) /= nx) .or. (size (a_im, 2) /= ny)) &
          call fatal_error ('fft_xy_parallel_4D', 'imaginary array size mismatch /= nx,ny', lfirst_proc_xy)
!
      if (mod (nxgrid, nprocxy) /= 0) &
          call fatal_error ('fft_xy_parallel_4D', 'nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) &
          call fatal_error ('fft_xy_parallel_4D', 'nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_4d(pnx,pny,inz,ina), p_im_4d(pnx,pny,inz,ina), t_re_4d(tnx,tny,inz,ina), t_im_4d(tnx,tny,inz,ina), stat=stat)
      if (stat > 0) call fatal_error ('fft_xy_parallel_4D', 'Could not allocate p and t', .true.)
      !$omp end single
      !$omp barrier
!
      if (lshear_loc) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        deltay_x = -deltay * (xgrid(x_offset:x_offset+tny-1) - (x0+Lx/2))/Lx
      endif
!
      if (lshift) then
        x_offset = 1 + (ipx+ipy*nprocx)*tny
        dshift_y = shift_y(x_offset:x_offset+tny-1)
      endif
!
      if (lforward) then
!
!  Forward FFT:
!
        ! Remap the data we need into pencil shape.
        
        call remap_to_pencil_xy (a_re, p_re_4d )
        if (lcompute_im) then
          call remap_to_pencil_xy (a_im, p_im_4d )
        else
          !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
          !$omp workshare
          p_im_4d= 0.0
          !$omp end workshare
          !!$omp end parallel
        endif
!
        call transp_pencil_xy (p_re_4d, t_re_4d)
        call transp_pencil_xy (p_im_4d, t_im_4d)
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do pos_z = 1, inz
            do l = 1, tny
              ! Transform y-direction.
              ay = cmplx (t_re_4d(:,l,pos_z,pos_a), t_im_4d(:,l,pos_z,pos_a))
              call cfftf (nygrid, ay, wsavey(:,thread_id))
              if (lshear_loc) ay = ay * exp (cmplx (0, ky_fft * deltay_x(l)))
              if (lshift) ay = ay * exp (cmplx (0,-ky_fft * dshift_y(l)))
              t_re_4d(:,l,pos_z,pos_a) = real (ay)
              t_im_4d(:,l,pos_z,pos_a) = aimag (ay)
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy (t_re_4d, p_re_4d)
        call transp_pencil_xy (t_im_4d, p_im_4d)
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do pos_z = 1, inz
            do m = 1, pny
              ! Transform x-direction and normalize.
              ax = cmplx (p_re_4d(:,m,pos_z,pos_a), p_im_4d(:,m,pos_z,pos_a))
              call cfftf (nxgrid, ax, wsavex(:,thread_id))
              p_re_4d(:,m,pos_z,pos_a) = real (ax)/nxygrid
              p_im_4d(:,m,pos_z,pos_a) = aimag (ax)/nxygrid
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re_4d, a_re )
        call unmap_from_pencil_xy (p_im_4d, a_im )
!
      else
!
!  Inverse FFT:
!
        ! Remap the data we need into transposed pencil shape.
        call remap_to_pencil_xy (a_re, p_re_4d )
        call remap_to_pencil_xy (a_im, p_im_4d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do pos_z = 1, inz
            do m = 1, pny
              ! Transform x-direction back.
              ax = cmplx (p_re_4d(:,m,pos_z,pos_a), p_im_4d(:,m,pos_z,pos_a))
              call cfftb (nxgrid, ax, wsavex(:,thread_id))
              p_re_4d(:,m,pos_z,pos_a) = real (ax)
              p_im_4d(:,m,pos_z,pos_a) = aimag (ax)
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy (p_re_4d, t_re_4d)
        call transp_pencil_xy (p_im_4d, t_im_4d)
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do pos_z = 1, inz
            do l = 1, tny
              ! Transform y-direction back.
              ay = cmplx (t_re_4d(:,l,pos_z,pos_a), t_im_4d(:,l,pos_z,pos_a))
              if (lshear_loc) ay = ay * exp (cmplx (0, -ky_fft*deltay_x(l)))
              call cfftb (nygrid, ay, wsavey(:,thread_id))
              t_re_4d(:,l,pos_z,pos_a) = real (ay)
              if (lcompute_im) t_im_4d(:,l,pos_z,pos_a) = aimag (ay)
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        call transp_pencil_xy (t_re_4d, p_re_4d)
        if (lcompute_im) call transp_pencil_xy (t_im_4d, p_im_4d)
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_xy (p_re_4d, a_re)
        if (lcompute_im) call unmap_from_pencil_xy (p_im_4d, a_im)
!
      endif
!
      !$omp barrier
      !$omp single
      deallocate (p_re_4d, p_im_4d, t_re_4d, t_im_4d)
      !$omp end single
!
    endsubroutine fft_xy_parallel_4D
!***********************************************************************
    subroutine fft_xyz_parallel_3D(a_re,a_im,linv,lneed_im,shift_y,lignore_shear)
!
!  Subroutine to do FFT of distributed 3D data in the x-, y- and z-direction.
!  For the x- and y-direction, the subroutine 'fft_xy_parallel' is used.
!  For z-parallelization the transform in the z-direction will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  In addition to the restrictions of 'fft_xy_parallel' (see there),
!  ny is restricted to be an integer multiple of nprocz.
!  linv indicates forward (=false, default) or backward (=true) transform.
!
!  27-oct-2010/Bourdin.KIS: adapted from fft_xyz_parallel_4D
!
      use Mpicomm, only: remap_to_pencil_yz, unmap_from_pencil_yz
!
      real, dimension (:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      
      real, dimension (nxgrid), optional :: shift_y
!
      integer, parameter :: pny=ny/nprocz, pnz=nzgrid    ! z-pencil shaped data sizes
!
      integer :: l, m, stat
      logical :: lforward, lcompute_im, lshift, lnoshear
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_y)) lshift = .true.
!
      lnoshear = .false.
      if (present (lignore_shear)) lnoshear = lignore_shear
!
      ! Check for degenerate cases.
      if (nzgrid == 1) then
        if (lshift) then
          call fft_xy_parallel (a_re(:,:,1), a_im(:,:,1), .not. lforward, lcompute_im, shift_y, lignore_shear=lnoshear)
        else
          call fft_xy_parallel (a_re(:,:,1), a_im(:,:,1), .not. lforward, lcompute_im, lignore_shear=lnoshear)
        endif
        return
      endif
      if ((nxgrid == 1) .and. (nygrid == 1)) then
        !$omp single    ! as a_* is shared and no worksharing inside fft_y_parallel_1D
        call fft_z_parallel (a_re(1,1,:), a_im(1,1,:), .not. lforward, lcompute_im, lignore_shear=lnoshear)
        !$omp end single
        return
      endif
!
      if (size (a_re, 3) /= nz) &
          call fatal_error ('fft_xyz_parallel_3D', 'real array size mismatch /= nx,ny,nz', lroot)
      if (size (a_im, 3) /= nz) &
          call fatal_error ('fft_xyz_parallel_3D', 'imaginary array size mismatch /= nx,ny,nz', lroot)
      if (mod (nygrid, nprocyz) /= 0) &
          call fatal_error ('fft_xyz_parallel_3D', 'nygrid needs to be an integer multiple of nprocy*nprocz', lroot)
!
      if (lforward) then
!
!  Forward FFT:
!
        if (lshift) then
          call fft_xy_parallel (a_re, a_im, .not. lforward, lcompute_im, shift_y, lignore_shear=lnoshear)
        else
          call fft_xy_parallel (a_re, a_im, .not. lforward, lcompute_im, lignore_shear=lnoshear)
        endif
!
        ! Allocate memory for large arrays.
        !$omp single
        allocate (p_re_3d(nx,pny,pnz), p_im_3d(nx,pny,pnz), stat=stat)
        if (stat > 0) call fatal_error ('fft_xyz_parallel_3D', 'Could not allocate p', .true.)
        !$omp end single
        !$omp barrier
!
        ! Remap the data we need into z-pencil shape.
!
        call remap_to_pencil_yz (a_re, p_re_3d )
        call remap_to_pencil_yz (a_im, p_im_3d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do l = 1, nx
          do m = 1, pny
            ! Transform z-direction and normalize.
            az = cmplx (p_re_3d(l,m,:), p_im_3d(l,m,:))
            call cfftf (nzgrid, az, wsavez(:,thread_id))
            p_re_3d(l,m,:) = real (az)/nzgrid
            p_im_3d(l,m,:) = aimag (az)/nzgrid
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_yz (p_re_3d, a_re )
        call unmap_from_pencil_yz (p_im_3d, a_im )
!
        !$omp barrier
        !$omp single
        deallocate (p_re_3d, p_im_3d)
        !$omp end single
!
      else
!
!  Inverse FFT:
!
        ! Allocate memory for large arrays.
        !$omp single
        allocate (p_re_3d(nx,pny,pnz), p_im_3d(nx,pny,pnz), stat=stat)
        if (stat > 0) call fatal_error ('fft_xyz_parallel_3D', 'Could not allocate p', .true.)
        !$omp end single
        !$omp barrier
!
        ! Remap the data we need into transposed z-pencil shape.
        
        call remap_to_pencil_yz (a_re, p_re_3d )
        call remap_to_pencil_yz (a_im, p_im_3d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2) 
        do l = 1, nx
          do m = 1, pny
            ! Transform z-direction back.
            az = cmplx (p_re_3d(l,m,:), p_im_3d(l,m,:))
            call cfftb (nzgrid, az, wsavez(:,thread_id))
            p_re_3d(l,m,:) = real (az)
            p_im_3d(l,m,:) = aimag (az)
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_yz (p_re_3d, a_re )
        call unmap_from_pencil_yz (p_im_3d, a_im )
!
        !$omp barrier
        !$omp single
        deallocate (p_re_3d, p_im_3d)
        !$omp end single
!
        if (lshift) then
          call fft_xy_parallel (a_re, a_im, .not. lforward, lcompute_im, shift_y, lignore_shear=lnoshear)
        else
          call fft_xy_parallel (a_re, a_im, .not. lforward, lcompute_im, lignore_shear=lnoshear)
        endif
!
      endif
!
    endsubroutine fft_xyz_parallel_3D
!***********************************************************************
    subroutine fft_xyz_parallel_4D(a_re,a_im,linv,lneed_im,shift_y,lignore_shear)
!
!  Subroutine to do FFT of distributed 4D data in the x- and y-direction.
!  For the x- and y-direction, the subroutine 'fft_xy_parallel' is used.
!  For z-parallelization the transform in the z-direction will be done under
!  MPI in parallel on all processors of the corresponding xy-plane.
!  In addition to the restrictions of 'fft_xy_parallel' (see there),
!  ny is restricted to be an integer multiple of nprocz.
!  linv indicates forward (=false, default) or backward (=true) transform.
!  Attention: input data will be overwritten.
!
!  27-oct-2010/Bourdin.KIS: coded
!
      use Mpicomm, only: remap_to_pencil_yz, unmap_from_pencil_yz
!
      real, dimension (:,:,:,:), intent(inout) :: a_re, a_im
      logical, optional, intent(in) :: linv, lneed_im, lignore_shear
      
      real, dimension (nxgrid), optional :: shift_y
!
      integer, parameter :: pny=ny/nprocz, pnz=nzgrid      ! z-pencil shaped data sizes
     
      integer :: ina ! size of the fourth dimension
!
      integer :: l, m, stat, pos_a
      logical :: lforward, lcompute_im, lshift, lnoshear
!
      lforward = .true.
      if (present (linv)) lforward = .not. linv
!
      lcompute_im = .true.
      if (present (lneed_im)) lcompute_im = lneed_im
!
      lshift = .false.
      if (present (shift_y)) lshift = .true.
!
      lnoshear = .false.
      if (present (lignore_shear)) lnoshear = lignore_shear
!
      ! Check for degenerate cases.
      if (nzgrid == 1) then
        if (lshift) then
          call fft_xy_parallel (a_re(:,:,1,:), a_im(:,:,1,:), .not. lforward, lcompute_im, shift_y, lignore_shear=lnoshear)
        else
          call fft_xy_parallel (a_re(:,:,1,:), a_im(:,:,1,:), .not. lforward, lcompute_im, lignore_shear=lnoshear)
        endif
        return
      endif
      if ((nxgrid == 1) .and. (nygrid == 1)) then
        call fft_z_parallel (a_re(1,1,:,:), a_im(1,1,:,:), .not. lforward, lcompute_im, lignore_shear=lnoshear)
        return
      endif
!
      ina = size (a_re, 4)
!
      if (size (a_re, 3) /= nz) &
          call fatal_error ('fft_xyz_parallel_4D', 'real array size mismatch /= nx,ny,nz', lroot)
      if (size (a_im, 3) /= nz) &
          call fatal_error ('fft_xyz_parallel_4D', 'imaginary array size mismatch /= nx,ny,nz', lroot)
!
      if (mod (nygrid, nprocyz) /= 0) &
          call fatal_error ('fft_xyz_parallel_4D', 'nygrid needs to be an integer multiple of nprocy*nprocz', lroot)
!
      if (lforward) then
!
!  Forward FFT:
!
        if (lshift) then
          call fft_xy_parallel (a_re, a_im, .not. lforward, lcompute_im, shift_y, lignore_shear=lnoshear)
        else
          call fft_xy_parallel (a_re, a_im, .not. lforward, lcompute_im, lignore_shear=lnoshear)
        endif
!
        ! Allocate memory for large arrays.
        !$omp single
        allocate (p_re_4d(nx,pny,pnz,ina), p_im_4d(nx,pny,pnz,ina), stat=stat)
        if (stat > 0) call fatal_error ('fft_xyz_parallel_4D', 'Could not allocate p', .true.)
        !$omp end single
        !$omp barrier
!
        ! Remap the data we need into z-pencil shape.
        
        call remap_to_pencil_yz (a_re, p_re_4d )
        call remap_to_pencil_yz (a_im, p_im_4d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do l = 1, nx
            do m = 1, pny
              ! Transform z-direction and normalize.
              az = cmplx (p_re_4d(l,m,:,pos_a), p_im_4d(l,m,:,pos_a))
              call cfftf (nzgrid, az, wsavez(:,thread_id))
              p_re_4d(l,m,:,pos_a) = real (az)/nzgrid
              p_im_4d(l,m,:,pos_a) = aimag (az)/nzgrid
            enddo
          enddo
        enddo
        !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_yz (p_re_4d, a_re )
        call unmap_from_pencil_yz (p_im_4d, a_im )
!
        !$omp barrier
        !$omp single
        deallocate (p_re_4d, p_im_4d)
        !$omp end single
!
      else
!
!  Inverse FFT:
!
        ! Allocate memory for large arrays.
        !$omp single
        allocate (p_re_4d(nx,pny,pnz,ina), p_im_4d(nx,pny,pnz,ina), stat=stat)
        if (stat > 0) call fatal_error ('fft_xyz_parallel_4D', 'Could not allocate p', .true.)
        !$omp end single
        !$omp barrier
!
        ! Remap the data we need into transposed pencil shape.
        
        call remap_to_pencil_yz (a_re, p_re_4d )
        call remap_to_pencil_yz (a_im, p_im_4d )
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(3)
        do pos_a = 1, ina
          do l = 1, nx
            do m = 1, pny
              ! Transform z-direction back.
              az = cmplx (p_re_4d(l,m,:,pos_a), p_im_4d(l,m,:,pos_a))
              call cfftb (nzgrid, az, wsavez(:,thread_id))
              p_re_4d(l,m,:,pos_a) = real (az)
              p_im_4d(l,m,:,pos_a) = aimag (az)
            enddo
          enddo
        enddo
       !!$omp end parallel
!
        ! Unmap the results back to normal shape.
        call unmap_from_pencil_yz (p_re_4d, a_re )
        call unmap_from_pencil_yz (p_im_4d, a_im )
!
        !$omp barrier
        !$omp single
        deallocate (p_re_4d, p_im_4d)
        !$omp end single
!
        if (lshift) then
          call fft_xy_parallel (a_re, a_im, .not. lforward, lcompute_im, shift_y, lignore_shear=lnoshear)
        else
          call fft_xy_parallel (a_re, a_im, .not. lforward, lcompute_im, lignore_shear=lnoshear)
        endif
!
      endif
!
    endsubroutine fft_xyz_parallel_4D
!***********************************************************************
    subroutine setup_extrapol_fact(z,ref_z,factor,reduce)
!
!  Subroutine to setup 'factor' for z-extrapolation of a vector potential.
!  'factor' is the multiplication factor for the Fourier coefficients,
!  including the normalization.
!  'z' gives the z-coordinates.
!  'ref_z' gives the reference z-coordinate for the extrapolation.
!  'reduce' is used to reduce the eventual enhancement of contrast
!  in cases where a z is smaller than ref_z (="intrapolation").
!  nx is restricted to be an integer multiple of nprocy.
!  ny is restricted to be an integer multiple of nprocx.
!
!  25-jan-2011/Bourdin.KIS: coded
!
      real, dimension (:), intent(in) :: z
      real, intent(in) :: ref_z
      real, dimension (:,:,:), intent(out) :: factor
      real, intent(in), optional :: reduce
!
      integer, parameter :: enx=nygrid, eny=nx/nprocy ! transposed data in pencil shape
      integer :: onz, pos_z, kx_start, alloc_err
      real :: delta_z
      real, dimension (:,:), allocatable :: k_2
!
      onz = size (z, 1)
!
      if (size (factor, 3) /= onz) &
          call fatal_error ('setup_extrapol_fact', 'z dimension differs between z and factor', lfirst_proc_xy)
      if ((size (factor, 1) /= enx) .or. (size (factor, 2) /= eny)) &
          call fatal_error ('setup_extrapol_fact', 'factor x/y-dimension is invalid', lfirst_proc_xy)
!
      !$omp single
      allocate (k_2(enx,eny), stat=alloc_err)
      if (alloc_err > 0) call fatal_error ('setup_extrapol_fact', 'Could not allocate k_2', .true.)
      !$omp end single
      !$omp barrier
!
      ! Get wave numbers in transposed pencil shape and calculate exp(|k|)
      kx_start = (ipx + ipy*nprocx)*eny
      k_2 = spread (ky_fft, 2, eny)**2 + spread (kx_fft(kx_start+1:kx_start+eny), 1, enx)**2
      ! Set dummy value to avoid later division by zero:
      if (kx_start == 0) k_2(1,1) = 1.0
!
      ! Setup fourier coefficients factor for extrapolation/intrapolation
      factor(:,:,onz) = exp (sqrt (k_2))
      if (kx_start == 0) then
        ! Special case for kx=0 and ky=0
        factor(1,1,onz) = 1.0
      endif
      !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
      !$omp do
      do pos_z = 1, onz
        delta_z = ref_z - z(pos_z)
        ! delta_z negative => decay of contrast
        ! delta_z positive => enhancement of contrast (can be reduced)
        if (present (reduce) .and. (delta_z > 0.0)) delta_z = delta_z * reduce
        ! Include normalization of FFT: 1/(nxgrid*nygrid)
        factor(:,:,pos_z) = factor(:,:,onz) ** delta_z / (k_2 * nxgrid*nygrid)
      enddo
      !!$omp end parallel
!
      !$omp barrier
      !$omp single
      deallocate (k_2)
      !$omp end single
!
    endsubroutine setup_extrapol_fact
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
      
      real, dimension (:,:,:,:), allocatable :: e_re, e_im ! extrapolated data in transposed pencil shape
      real, dimension (:,:,:,:), allocatable :: b_re, b_im ! backtransformed data in pencil shape
      complex, dimension (nygrid) :: ay_extra
      integer :: ina ! number of components in the output data (usually 3)
      integer :: onz, ona ! number of ghost cells and components in the output data (usually 3)
      integer :: l, m, stat, pos_a, pos_z
!
      ina = size (in, 3)
      onz = size (out, 3)
      ona = size (out, 4)
!
      if (ina /= ona) call fatal_error ('vect_pot_extrapol_z_parallel', &
              'number of components is different for input and ouput arrays', lfirst_proc_xy)
!
      if ((size (in, 1) /= nx) .or. (size (in, 2) /= ny)) &
          call fatal_error ('vect_pot_extrapol_z_parallel', 'input array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (out, 1) /= nx) .or. (size (out, 2) /= ny)) &
          call fatal_error ('vect_pot_extrapol_z_parallel', 'output array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (factor, 1) /= tnx) .or. (size (factor, 2) /= tny)) &
          call fatal_error ('vect_pot_extrapol_z_parallel', 'factor array size mismatch /= tnx,tny', lfirst_proc_xy)
      if (size (factor, 3) /= onz) &
          call fatal_error ('vect_pot_extrapol_z_parallel', &
              'number of ghost cells differs between multiplication factor and ouput array', lfirst_proc_xy)
!
      if (mod (nxgrid, nprocxy) /= 0) call fatal_error ('vect_pot_extrapol_z_parallel', &
              'nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) call fatal_error ('vect_pot_extrapol_z_parallel', &
              'nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
!
      if (lshear) call not_implemented('vect_pot_extrapol_z_parallel','shearing',lfirst_proc_xy) 
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_3d(pnx,pny,ona), p_im_3d(pnx,pny,ona), t_re_3d(tnx,tny,ona), t_im_3d(tnx,tny,ona), stat=stat)
      if (stat > 0) call fatal_error ('vect_pot_extrapol_z_parallel', 'Could not allocate p and t', .true.)
      !$omp end single
      !$omp barrier
!
      ! Collect the data we need.
      
      call remap_to_pencil_xy (in, p_re_3d )
      p_im_3d= 0.
!
      !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
      !$omp do collapse(2)
      do pos_a = 1, ona
        do m = 1, pny
          ! Transform x-direction.
          ax = cmplx (p_re_3d(:,m,pos_a), p_im_3d(:,m,pos_a))
          call cfftf (nxgrid, ax, wsavex(:,thread_id))
          p_re_3d(:,m,pos_a) = real (ax)
          p_im_3d(:,m,pos_a) = aimag (ax)
        enddo
      enddo
      !!$omp end parallel
!
      call transp_pencil_xy (p_re_3d, t_re_3d )
      call transp_pencil_xy (p_im_3d, t_im_3d )
!
      !$omp barrier
      !$omp single
      deallocate (p_re_3d, p_im_3d)
      allocate (e_re(tnx,tny,onz,ona), e_im(tnx,tny,onz,ona), stat=stat)
      if (stat > 0) call fatal_error ('vect_pot_extrapol_z_parallel', 'Could not allocate e', .true.)
      !$omp end single
      !$omp barrier
!
      !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
      !$omp do collapse(2)
      do pos_a = 1, ona
        do l = 1, tny
          ! Transform y-direction.
          ay = cmplx (t_re_3d(:,l,pos_a), t_im_3d(:,l,pos_a))
          call cfftf (nygrid, ay, wsavey(:,thread_id))
          do pos_z = 1, onz
            ! Apply factor to fourier coefficients.
            ay_extra = ay * factor(:,l,pos_z)
            ! Transform y-direction back in each z layer.
            call cfftb (nygrid, ay_extra, wsavey(:,thread_id))
            e_re(:,l,pos_z,pos_a) = real (ay_extra)
            e_im(:,l,pos_z,pos_a) = aimag (ay_extra)
          enddo
        enddo
      enddo
      !!$omp end parallel
!
      !$omp barrier
      !$omp single
      deallocate (t_re_3d, t_im_3d)
      allocate (b_re(pnx,pny,onz,ona), b_im(pnx,pny,onz,ona), stat=stat)
      if (stat > 0) call fatal_error ('vect_pot_extrapol_z_parallel', 'Could not allocate b', .true.)
      !$omp end single
      !$omp barrier
!
      call transp_pencil_xy (e_re, b_re)
      call transp_pencil_xy (e_im, b_im)
!
      !$omp barrier
      !$omp single
      deallocate (e_re, e_im)
      !$omp end single
!
      ! Transform x-direction back in each z layer.
      !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
      !$omp do collapse(3)
      do pos_a = 1, ona
        do pos_z = 1, onz
          do m = 1, pny
            ax = cmplx (b_re(:,m,pos_z,pos_a), b_im(:,m,pos_z,pos_a))
            call cfftb (nxgrid, ax, wsavex(:,thread_id))
            b_re(:,m,pos_z,pos_a) = real (ax)
          enddo
        enddo
      enddo
      !!$omp end parallel
!
      ! Distribute the results back in normal shape.
      call unmap_from_pencil_xy (b_re, out )
!
      !$omp barrier
      !$omp single
      deallocate (b_re, b_im)
      !$omp end single
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
      
      real, dimension (:,:,:,:), allocatable :: e_re, e_im ! extrapolated data in transposed pencil shape
      real, dimension (:,:,:,:), allocatable :: b_re, b_im ! backtransformed data in pencil shape
      complex, dimension (nygrid) :: ay_extra, ay_extra_x, ay_extra_y
      integer :: onz ! number of layers in the output data
      integer :: l, m, stat, pos_z
!
      onz = size (out, 3)
!
      if ((size (in, 1) /= pnx) .or. (size (in, 2) /= pny)) &
          call fatal_error ('field_extrapol_z_parallel', &
              'input array size mismatch /= pnx,pny', lfirst_proc_xy)
      if ((size (out, 1) /= nx) .or. (size (out, 2) /= ny)) &
          call fatal_error ('field_extrapol_z_parallel', &
               'output array size mismatch /= nx,ny', lfirst_proc_xy)
      if (size (out, 4) /= 2) call fatal_error ('field_extrapol_z_parallel', &
              'output array must have two components', lfirst_proc_xy)
      if ((size (factor, 1) /= tnx) .or. (size (factor, 2) /= tny)) &
          call fatal_error ('field_extrapol_z_parallel', &
              'factor array size mismatch /= pnx,pny', lfirst_proc_xy)
      if (size (factor, 3) /= onz) call fatal_error ('field_extrapol_z_parallel', &
              'number of ghost cells differs between multiplication factor and ouput array', lfirst_proc_xy)
!
      if (mod (nxgrid, nprocxy) /= 0) call fatal_error ('field_extrapol_z_parallel', &
              'nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
      if (mod (nygrid, nprocxy) /= 0) call fatal_error ('field_extrapol_z_parallel', &
              'nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
!
      if (lshear) call not_implemented('field_extrapol_z_parallel','shearing', lfirst_proc_xy) 
!
      ! Allocate memory for large arrays.
      !$omp single
      allocate (p_re_2d(pnx,pny), p_im_2d(pnx,pny), t_re_2d(tnx,tny), t_im_2d(tnx,tny), stat=stat)
      if (stat > 0) call fatal_error ('field_extrapol_z_parallel', 'Could not allocate p and t', .true.)
      !$omp end single
      !$omp barrier
!
      ! Collect the data we need.
      !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
      !$omp workshare
      p_re_2d= in
      p_im_2d= 0.
      !$omp end workshare
      !!$omp end parallel
!
      !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
      !$omp do
      do m = 1, pny
        ! Transform x-direction.
        ax = cmplx (p_re_2d(:,m), p_im_2d(:,m))
        call cfftf (nxgrid, ax, wsavex(:,thread_id))
        p_re_2d(:,m) = real (ax)
        p_im_2d(:,m) = aimag (ax)
      enddo
      !!$omp end parallel
!
      call transp_pencil_xy (p_re_2d, t_re_2d )
      call transp_pencil_xy (p_im_2d, t_im_2d )
!
      !$omp barrier
      !$omp single
      deallocate (p_re_2d, p_im_2d)
      allocate (e_re(tnx,tny,onz,2), e_im(tnx,tny,onz,2), stat=stat)
      if (stat > 0) call fatal_error ('field_extrapol_z_parallel', 'Could not allocate e', .true.)
      !$omp end single
      !$omp barrier
!
      !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
      !$omp do
      do l = 1, tny
        ! Transform y-direction.
        ay = cmplx (t_re_2d(:,l), t_im_2d(:,l))
        call cfftf (nygrid, ay, wsavey(:,thread_id))
        ! Transform y-direction back in each z layer.
        do pos_z = 1, onz
          ! Apply factor to fourier coefficients.
          ay_extra = ay * factor(:,l,pos_z)
          ! x-component of A:
          ay_extra_x = cmplx (-aimag (ay_extra), real (ay_extra)) * ky_fft
          call cfftb (nygrid, ay_extra_x, wsavey(:,thread_id))
          e_re(:,l,pos_z,1) = real (ay_extra_x)
          e_im(:,l,pos_z,1) = aimag (ay_extra_x)
          ! y-component of A:
          ay_extra_y = cmplx (aimag (ay_extra), -real (ay_extra)) * kx_fft(l+(iproc-ipz*nprocxy)*tny)
          call cfftb (nygrid, ay_extra_y, wsavey(:,thread_id))
          e_re(:,l,pos_z,2) = real (ay_extra_y)
          e_im(:,l,pos_z,2) = aimag (ay_extra_y)
        enddo
      enddo
      !!$omp end parallel
!
      !$omp barrier
      !$omp single
      deallocate (t_re_2d, t_im_2d)
      allocate (b_re(pnx,pny,onz,2), b_im(pnx,pny,onz,2), stat=stat)
      if (stat > 0) call fatal_error ('field_extrapol_z_parallel', 'Could not allocate b', .true.)
      !$omp end single
      !$omp barrier
      
      call transp_pencil_xy (e_re, b_re )
      call transp_pencil_xy (e_im, b_im )
!
      !$omp barrier
      !$omp single
      deallocate (e_re, e_im)
      !$omp end single
!
      ! Transform x-direction back.
      !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
      !$omp do collapse(2)
      do pos_z = 1, onz
        do m = 1, pny
          ! x-component of A:
          ax = cmplx (b_re(:,m,pos_z,1), b_im(:,m,pos_z,1))
          call cfftb (nxgrid, ax, wsavex(:,thread_id))
          b_re(:,m,pos_z,1) = real (ax)
          ! y-component of A:
          ax = cmplx (b_re(:,m,pos_z,2), b_im(:,m,pos_z,2))
          call cfftb (nxgrid, ax, wsavex(:,thread_id))
          b_re(:,m,pos_z,2) = real (ax)
        enddo
      enddo
      !!$omp end parallel
!
      ! Distribute the results.
      call unmap_from_pencil_xy (b_re, out )
!
      !$omp barrier
      !$omp single
      deallocate (b_re, b_im)
      !$omp end single
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
      use General, only: find_proc
!
      real, dimension (ny) :: a_re, a_im
      logical, optional :: linv
!
      real, dimension (nygrid) :: a_re_full, a_im_full
      integer :: ipy_send,partner
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
              partner=find_proc(0,ipy+ipy_send,ipz)
              call mpirecv_real(a_re_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,partner,itag1 )   !iproc+ipy_send,itag1)
              call mpirecv_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,partner,itag2 )   !iproc+ipy_send,itag2)
            enddo
          else
            partner=find_proc(0,0,ipz)
            call mpisend_real(a_re,ny,partner,itag1 )
            call mpisend_real(a_im,ny,partner,itag2 )
          endif
!
          if (lfirst_proc_y) then
!
            ay=cmplx(a_re_full,a_im_full)
            call cfftf(nygrid,ay,wsavey(:,thread_id))
            a_re_full=real(ay)
            a_im_full=aimag(ay)
!
            a_re=a_re_full(1:ny)
            a_im=a_im_full(1:ny)
!
            do ipy_send=1,nprocy-1
              partner=find_proc(0,ipy+ipy_send,ipz)
              call mpisend_real(a_re_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,partner,itag1 )    !iproc+ipy_send,itag1)
              call mpisend_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,partner,itag2 )    !iproc+ipy_send,itag2)
            enddo
          else
            call mpirecv_real(a_re,ny,partner,itag1 )   !iproc-ipy,itag1)
            call mpirecv_real(a_im,ny,partner,itag2 )   !iproc-ipy,itag2)
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
              partner=find_proc(0,ipy+ipy_send,ipz)
              call mpirecv_real(a_re_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,partner,itag1 )    !iproc+ipy_send,itag1)
              call mpirecv_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,partner,itag2 )    !iproc+ipy_send,itag2)
            enddo
          else
            partner=find_proc(0,0,ipz)
            call mpisend_real(a_re,ny,partner,itag1 )   !iproc-ipy,itag1)
            call mpisend_real(a_im,ny,partner,itag2 )   !iproc-ipy,itag2)
          endif
!
          if (lfirst_proc_y) then
!
            ay=cmplx(a_re_full,a_im_full)
            call cfftb(nygrid,ay,wsavey(:,thread_id))
            a_re_full=real(ay)
            a_im_full=aimag(ay)
!
            a_re=a_re_full(1:ny)
            a_im=a_im_full(1:ny)
!
            do ipy_send=1,nprocy-1
              partner=find_proc(0,ipy+ipy_send,ipz)
              call mpisend_real(a_re_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,partner,itag1 )    !iproc+ipy_send,itag1)
              call mpisend_real(a_im_full(ipy_send*ny+1:(ipy_send+1)*ny), &
                  ny,partner,itag2 )    !iproc+ipy_send,itag2)
            enddo
          else
            call mpirecv_real(a_re,ny,partner,itag1 )    !iproc-ipy,itag1)
            call mpirecv_real(a_im,ny,partner,itag2 )    !iproc-ipy,itag2)
          endif
!
        endif
!
      endif
!
!  Normalize
!
      if (lforward) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp workshare
        a_re=a_re/nygrid
        a_im=a_im/nygrid
        !$omp end workshare
        !!$omp end parallel
      endif
!
    endsubroutine fourier_transform_y_y
!***********************************************************************
    subroutine fourier_shift_yz_y(a_re,shift_y)
!
!  Performs a periodic shift in the y-direction of an entire y-z plane by
!  the amount shift_y. The shift is done in Fourier space for maximum
!  interpolation accuracy. Data is re-distributed in elongated shape,
!  so that the FFT will be performed in parallel on all processors.
!
!  19-jul-06/anders: coded
!
      use Mpicomm, only: mpisend_real, mpirecv_real
      use General, only: find_proc
!
      real, dimension (ny,nz) :: a_re
      real :: shift_y
!
      complex, dimension (nygrid) :: a_cmplx, cmplx_shift
      real, dimension (nygrid,max(nz/nprocy,1)) :: a_re_new, a_im_new
      integer :: n, nz_new, ipy_from, ipy_to, iproc_from, iproc_to
      integer :: nprocy_used
      real, dimension (:,:), allocatable :: buffer
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
          nprocy_used = 1
          do ipy_from=1,nprocy-1
            if (lfirst_proc_y) then
              call mpirecv_real( a_re_new(ipy_from*ny+1:(ipy_from+1)*ny,1), &
                  ny,find_proc(ipx,ipy_from,0),itag )     !ipy_from*nprocx+ipx,itag)
            else
              if (ipy==ipy_from) &
                  call mpisend_real(a_re(:,1),ny,find_proc(ipx,0,0),itag )   !ipy_to*nprocx+ipx,itag)
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
          ! re-distribute data in elongated pencil shape
          ! all processors need to send their data portion
          ! each participating processor receives data
          !$omp single
          allocate (buffer(ny,nz_new))
          !$omp end single
          !$omp barrier
          do ipy_from=0,nprocy-1
            iproc_from=find_proc(ipx,ipy_from,ipz)   !ipz*nprocy*nprocx+ipy_from*nprocx+ipx
            if (ipy/=ipy_from) then
              if (ipy<nprocy_used) then
                call mpirecv_real(buffer,(/ny,nz_new/),iproc_from,itag )
                a_re_new(ipy_from*ny+1:(ipy_from+1)*ny,:) = buffer
              endif
            else
              if (ipy<nprocy_used) a_re_new(ipy*ny+1:(ipy+1)*ny,:) = &
                  a_re(:,ipy*nz_new+1:(ipy+1)*nz_new)
              do ipy_to=0,nprocy_used-1
                iproc_to=find_proc(ipx,ipy_to,ipz)   !ipz*nprocy*nprocx+ipy_to*nprocx+ipx
                if (ipy/=ipy_to) call mpisend_real( &
                    a_re(:,ipy_to*nz_new+1:(ipy_to+1)*nz_new),(/ny,nz_new/),iproc_to,itag )
              enddo
            endif
          enddo
        endif
      else
!
!  Only parallelization along z (or not at all).
!
        nprocy_used=1
        a_re_new(1:ny,1:nz_new)=a_re(1:ny,1:nz_new)
      endif
!
      ! perform FFT on all participating procesors
      if (ipy < nprocy_used) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp workshare
        a_im_new=0.0
        cmplx_shift=exp(cmplx(0.0,-ky_fft*shift_y))
        !$omp end workshare
        !!$omp end parallel
!
!  Transform to Fourier space.
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do
        do n=1,nz_new
          call fourier_transform_other(a_re_new(:,n),a_im_new(:,n))
          a_cmplx=cmplx(a_re_new(:,n),a_im_new(:,n))*cmplx_shift
          a_re_new(:,n)=real(a_cmplx)
          a_im_new(:,n)=aimag(a_cmplx)
        enddo
        !!$omp end parallel
!
!  Back to real space.
!
        do n=1,nz_new
          call fourier_transform_other(a_re_new(:,n),a_im_new(:,n),linv=.true.)
        enddo
      endif
!
!  Reinstate original division of yz-plane.
!
      if (nprocy/=1) then
        if (nzgrid==1) then
!  No z-direction.
          if (lfirst_proc_y) then
            do ipy_to=1,nprocy-1
              call mpisend_real( a_re_new(ipy_to*ny+1:(ipy_to+1)*ny,1), &
                  ny,find_proc(ipx,ipy_to,0),itag )    !ipy_to*nprocx+ipx,itag)
            enddo
          else
            call mpirecv_real(a_re(:,1),ny,find_proc(ipx,0,0),itag )  !ipx,itag)
          endif
          if (lfirst_proc_y) a_re(:,1)=a_re_new(1:ny,1)
        else
!  Present z-direction.
          ! distribute data back to their regular subdomains
          ! each participating processor needs to send data
          ! all processors receive their data portion
          do ipy_from=0,nprocy_used-1
            iproc_from=find_proc(ipx,ipy_from,ipz)    !ipz*nprocy*nprocx+ipy_from*nprocx+ipx
            if (ipy/=ipy_from) then
              call mpirecv_real( a_re(:,ipy_from*nz_new+1:(ipy_from+1)*nz_new), &
                  (/ny,nz_new/),iproc_from,itag+100 )
            else
              if (ipy<nprocy_used) a_re(:,ipy*nz_new+1:(ipy+1)*nz_new)= &
                  a_re_new(ipy*ny+1:(ipy+1)*ny,:)
              do ipy_to=0,nprocy-1
                iproc_to=find_proc(ipx,ipy_to,ipz)    !ipz*nprocy*nprocx+ipy_to*nprocx+ipx
                if (ipy/=ipy_to) then
                  buffer = a_re_new(ipy_to*ny+1:(ipy_to+1)*ny,:)
                  call mpisend_real(buffer,(/ny,nz_new/),iproc_to,itag+100 )
                endif
              enddo
            endif
          enddo
          !$omp barrier
          !$omp single
          deallocate (buffer)
          !$omp end single
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
!  MR: NOT USED; POTENTIALLY INCORRECT
!
!  04-oct-07/anders: adapted from fourier_transform_shear
!
      real, dimension (nx,ny,nz) :: a_re
      real, dimension (nx) :: shift_y
!
      real, dimension (nx,ny,nz) :: a_im
      complex, dimension (nx) :: ay
      real, dimension (4*nx+15) :: wsave
      integer :: l,n
      integer, parameter :: two=2
!
!  if nxgrid/=nygrid, then stop.
!
      if (nygrid/=nxgrid .and. nygrid/=1) &
        call fatal_error('fourier_shft_y','need to have nygrid=nxgrid if nygrid/=1')
!
!  Initialize cfft.
!
      call cffti(nygrid,wsave)    !MR: whatif nygrid/=nx?
!
!  Transform y-direction.
!
      if (ip<10) call information('fourier_shift_y','doing FFTpack in y')
      if (nygrid/=1) then
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp workshare
        a_im=0.0
        !$omp end workshare
        !!$omp end parallel
        call transp(a_re,'y' )
        
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do l=1,ny
          ay=cmplx(a_re(:,l,n),a_im(:,l,n))
          call cfftf(nygrid,ay,wsave)   !MR: whatif nygrid/=nx? See also next statement!
!
!  Shift all modes by the amount shift_y(x).
!
          !!!ay(two:nxgrid)=ay(two:nxgrid)*exp(cmplx(0.0,-ky_fft(two:nxgrid)*shift_y(l+ipy*ny)))
          a_re(:,l,n)=real(ay)
          a_im(:,l,n)=aimag(ay)
        enddo; enddo
        !!$omp end parallel
!
!  Transform y-direction back and normalize.
!
        !!$omp parallel num_threads(num_helper_threads) if(.not.omp_in_parallel())
        !$omp do collapse(2)
        do n=1,nz; do l=1,ny
          ay=cmplx(a_re(:,l,n),a_im(:,l,n))
          call cfftb(nygrid,ay,wsave)    !MR: whatif nygrid/=nx?
          a_re(:,l,n)=real(ay)/nygrid
          a_im(:,l,n)=aimag(ay)/nygrid
        enddo; enddo
        !!$omp end parallel

        call transp(a_re,'y' )
        call transp(a_im,'y' )
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
      endif
      if (ip<10) call information('fourier_transform_real_1','doing FFTpack in x')
!
!  Transform x-direction and normalize.
!
      if (lforward) then
        call rfftf(na,a,wsavex_temp)
        a=a/na
      else
!
!  Transform x-direction back.
!
        call rfftb(na,a,wsavex_temp)
      endif
!
      if (ip<10) call information('fourier_transform_real_1','fft has finished')
!
    endsubroutine fourier_transform_real_1
!***********************************************************************
endmodule Fourier
