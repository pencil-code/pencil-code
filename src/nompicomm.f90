! Id: nompicomm.f90,v 1.35 2002/08/16 21:23:48 brandenb Exp $

!!!!!!!!!!!!!!!!!!!!!!!
!!!  nompicomm.f90  !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Module with dummy MPI stuff.
!!!  This allows code to be run on single cpu machine

module Mpicomm

  use Cparam
  use Cdata, only: iproc,ipx,ipy,ipz,root,lroot

  implicit none

  interface mpibcast_real      ! Overload the `mpibcast_real' function
    module procedure mpibcast_real_scl
    module procedure mpibcast_real_arr
  endinterface

  interface mpibcast_logical   ! Overload
    module procedure mpibcast_logical_scl
  endinterface

  integer :: ierr

  contains

!***********************************************************************
    subroutine mpicomm_init()
!
!  Before the communication has been completed, the nghost=3 layers next
!  to the processor boundary (m1, m2, n1, or n2) cannot be used yet.
!  In the mean time we can calculate the interior points sufficiently far
!  away from the boundary points. Here we calculate the order in which
!  m and n are executed. At one point, necessary(imn)=.true., which is
!  the moment when all communication must be completed.
!
!   6-jun-02/axel: generalized to allow for ny=1
!  23-nov-02/axel: corrected problem with ny=4 or less
!
      use General
      use Cdata, only: lmpicomm
!
!  sets iproc in order that we write in the correct directory
!
!  consistency check
!
      if (ncpus > 1) then
        call stop_it("Inconsistency: MPICOMM=nompicomm, but ncpus >= 2")
      endif
!
!  for single cpu machine, set processor to zero
!
      lmpicomm = .false.
      iproc = 0
      lroot = .true.
      ipx = 0
      ipy = 0
      ipz = 0
!
      call setup_mm_nn()
!
    endsubroutine mpicomm_init
!***********************************************************************
    subroutine initiate_isendrcv_bdry(f)
!
      use Cdata
!
!  for one processor, use periodic boundary conditions
!  but in this dummy routine this is done in finalise_isendrcv_bdry
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: dummy
!
      dummy=f(1,1,1,1)  !(prevent compiler warning "unused variable ...")
!
    endsubroutine initiate_isendrcv_bdry
!***********************************************************************
    subroutine finalise_isendrcv_bdry(f)
!
      use Cparam
!
!  apply boundary conditions
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (ip==0) print*,'FINALIZE_ISENDRCV_BDRY: f=',f
    endsubroutine finalise_isendrcv_bdry
!***********************************************************************
    subroutine initiate_shearing(f)
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: dummy
!    
      dummy=f(1,1,1,1)
    endsubroutine initiate_shearing
!***********************************************************************
    subroutine finalise_shearing(f)
!
  use Cdata
!
!  for one processor, use periodic boundary conditions
!  but in this dummy routine this is done in finalise_isendrcv_bdry
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      double precision    :: deltay_dy, frak, c1, c2, c3, c4, c5, c6
      integer :: displs
!
!  Periodic boundary conditions in x, with shearing sheat
!
      if (nygrid==1) then !If 2D
         f( 1:l1-1,:,:,:) = f(l2i:l2,:,:,:)
         f(l2+1:mx,:,:,:) = f(l1:l1i,:,:,:)
      else
         deltay_dy=deltay/dy
         displs=int(deltay_dy)
         frak=deltay_dy-displs
         c1 = -          (frak+1.)*frak*(frak-1.)*(frak-2.)*(frak-3.)/120.
         c2 = +(frak+2.)          *frak*(frak-1.)*(frak-2.)*(frak-3.)/24.
         c3 = -(frak+2.)*(frak+1.)     *(frak-1.)*(frak-2.)*(frak-3.)/12.
         c4 = +(frak+2.)*(frak+1.)*frak          *(frak-2.)*(frak-3.)/12.
         c5 = -(frak+2.)*(frak+1.)*frak*(frak-1.)          *(frak-3.)/24.
         c6 = +(frak+2.)*(frak+1.)*frak*(frak-1.)*(frak-2.)          /120.
         f( 1:l1-1,m1:m2,:,:)=c1*cshift(f(l2i:l2,m1:m2,:,:),-displs+2,2) &
                             +c2*cshift(f(l2i:l2,m1:m2,:,:),-displs+1,2) &
                             +c3*cshift(f(l2i:l2,m1:m2,:,:),-displs,2) &
                             +c4*cshift(f(l2i:l2,m1:m2,:,:),-displs-1,2) &
                             +c5*cshift(f(l2i:l2,m1:m2,:,:),-displs-2,2) &
                             +c6*cshift(f(l2i:l2,m1:m2,:,:),-displs-3,2)  
         f(l2+1:mx,m1:m2,:,:)=c1*cshift(f(l1:l1i,m1:m2,:,:),displs-2,2) &
                             +c2*cshift(f(l1:l1i,m1:m2,:,:),displs-1,2) &
                             +c3*cshift(f(l1:l1i,m1:m2,:,:),displs,2) &
                             +c4*cshift(f(l1:l1i,m1:m2,:,:),displs+1,2) &
                             +c5*cshift(f(l1:l1i,m1:m2,:,:),displs+2,2) &
                             +c6*cshift(f(l1:l1i,m1:m2,:,:),displs+3,2) 
      end if
    end subroutine finalise_shearing
!***********************************************************************
    subroutine radboundary_zx_recv(rady0,mrad,Ibuf_zx,taubuf_zx)
!
!   2-jul-03/tony: dummy created
!
      integer :: rady0,mrad
      real, dimension(mx,rady0,mz) :: Ibuf_zx
      real, dimension(mx,rady0,mz), optional :: taubuf_zx
!
      if (ip==0) then
         print*,rady0,mrad,Ibuf_zx
         if (present(taubuf_zx)) print*,taubuf_zx
      endif
!
    endsubroutine radboundary_zx_recv
!***********************************************************************
    subroutine radboundary_xy_recv(radz0,nrad,Ibuf_xy,taubuf_xy)
!
!   2-jul-03/tony: dummy created
!
      integer :: radz0,nrad
      real, dimension(mx,my,radz0) :: Ibuf_xy
      real, dimension(mx,my,radz0), optional :: taubuf_xy
!
      if (ip==0) then
         print*,radz0,nrad,Ibuf_xy
         if (present(taubuf_xy)) print*,taubuf_xy
      endif
!
    endsubroutine radboundary_xy_recv
!***********************************************************************
    subroutine radboundary_zx_send(rady0,mrad,Ibuf_zx,taubuf_zx)
!
!   2-jul-03/tony: dummy created
!
      integer :: rady0,mrad
      real, dimension(mx,rady0,mz) :: Ibuf_zx
      real, dimension(mx,rady0,mz), optional :: taubuf_zx
!
      if (ip==0) then
         print*,rady0,mrad,Ibuf_zx
         if (present(taubuf_zx)) print*,taubuf_zx
      endif
!
    endsubroutine radboundary_zx_send
!***********************************************************************
    subroutine radboundary_xy_send(radz0,nrad,Ibuf_xy,taubuf_xy)
!
!   2-jul-03/tony: dummy created
!
      integer :: radz0,nrad
      real, dimension(mx,my,radz0) :: Ibuf_xy
      real, dimension(mx,my,radz0), optional :: taubuf_xy
!
      if (ip==0) then
         print*,radz0,nrad,Ibuf_xy
         if (present(taubuf_xy)) print*,taubuf_xy
      endif
!
    endsubroutine radboundary_xy_send
!***********************************************************************
    subroutine mpibcast_int(ibcast_array,nbcast_array)
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: ibcast_array,dummy
!    
      dummy=ibcast_array
    endsubroutine mpibcast_int
!***********************************************************************
    subroutine mpibcast_logical_scl(ibcast_array,nbcast_array)
!
      integer :: nbcast_array
      logical :: ibcast_array,dummy
!    
      if (nbcast_array/=1) stop "problem in mpibcast_real_scl"
      dummy=ibcast_array
    endsubroutine mpibcast_logical_scl
!***********************************************************************
    subroutine mpibcast_real_scl(bcast_array,nbcast_array)
!
      integer :: nbcast_array
      real :: bcast_array,dummy
!
      if (nbcast_array/=1) stop "problem in mpibcast_real_scl"
      dummy=bcast_array
    endsubroutine mpibcast_real_scl
!***********************************************************************
    subroutine mpibcast_real_arr(bcast_array,nbcast_array)
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array,dummy
!
      dummy=bcast_array
    endsubroutine mpibcast_real_arr
!***********************************************************************
    subroutine mpibcast_real_nonroot(bcast_array,nbcast_array,ibcast_proc)
!
      integer :: nbcast_array,ibcast_proc
      real, dimension(nbcast_array) :: bcast_array,dummy
!
      dummy=bcast_array
      dummy=ibcast_proc
    endsubroutine mpibcast_real_nonroot
!***********************************************************************
    subroutine mpireduce_max(fmax_tmp,fmax,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp, fmax
!
      fmax=fmax_tmp
    endsubroutine mpireduce_max
!***********************************************************************
    subroutine mpireduce_sum(fsum_tmp,fsum,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fsum_tmp,fsum
!
      fsum=fsum_tmp
    endsubroutine mpireduce_sum
!***********************************************************************
    subroutine start_serialize()
    endsubroutine start_serialize
!***********************************************************************
    subroutine end_serialize()
    endsubroutine end_serialize
!***********************************************************************
    subroutine mpibarrier()
    endsubroutine mpibarrier
!***********************************************************************
    subroutine mpifinalize()
    endsubroutine mpifinalize
!***********************************************************************
    function mpiwtime()
!
!  Mimic the MPI_WTIME() timer function. On many machines, the
!  implementation through system_clock() will overflow after about 50
!  minutes, so MPI_WTIME() is better.
!
!   5-oct-2002/wolf: coded
!
      double precision :: mpiwtime
      integer :: count_rate,time
!
      call system_clock(COUNT_RATE=count_rate)
      call system_clock(COUNT=time)

      if (count_rate /= 0) then
        mpiwtime = (time*1.)/count_rate
      else                      ! occurs with ifc 6.0 after long (> 2h) runs
        mpiwtime = 0
      endif
!
    endfunction mpiwtime
!***********************************************************************
    function mpiwtick()
!
!  Mimic the MPI_WTICK() function for measuring timer resolution.
!
!   5-oct-2002/wolf: coded
!
      double precision :: mpiwtick
      integer :: count_rate
!
      call system_clock(COUNT_RATE=count_rate)
      if (count_rate /= 0) then
        mpiwtick = 1./count_rate
      else                      ! occurs with ifc 6.0 after long (> 2h) runs
        mpiwtick = 0
      endif
!
    endfunction mpiwtick
!***********************************************************************
    subroutine stop_it(msg)
!
!  Print message and stop
!  6-nov-01/wolf: coded
!
      character (len=*) :: msg
!      
      if (lroot) write(0,'(A,A)') 'STOPPED: ', msg
      call mpifinalize
      STOP
    endsubroutine stop_it
!***********************************************************************
    subroutine transp(a,var)
!
!  Doing a transpose (dummy version for single processor)
!
!   5-sep-02/axel: adapted from version in mpicomm.f90
!
      real, dimension(nx,ny,nz) :: a
      real, dimension(nz) :: tmp_z
      real, dimension(ny) :: tmp_y
      integer :: i,j
      character :: var
!
      if (ip<10) print*,'transp for single processor'
!
!  Doing x-y transpose if var='y'
!
if (var=='y') then
!
      do i=1,ny
        do j=i+1,ny
          tmp_z=a(i,j,:)
          a(i,j,:)=a(j,i,:)
          a(j,i,:)=tmp_z
        enddo
      enddo
!
!  Doing x-z transpose if var='z'
!
elseif (var=='z') then
!
      do i=1,nz
        do j=i+1,nz
          tmp_y=a(i,:,j)
          a(i,:,j)=a(j,:,i)
          a(j,:,i)=tmp_y
        enddo
      enddo
!
endif
!
 end subroutine transp
!***********************************************************************
subroutine transform(a1,a2,a3,b1,b2,b3)
!
!  Subroutine to do fourier transform
!  The routine overwrites the input data
!
!  03-nov-02/nils: coded
!  05-nov-02/axel: added normalization factor
!
  real,dimension(nx,ny,nz) :: a1,a2,a3,b1,b2,b3
!
  if(lroot .AND. ip<10) print*,'doing fft of x-component'
  ! Doing the x field
  call fft(a1,b1, nx*ny*nz, nx, nx      ,-1) ! x-direction
  call fft(a1,b1, nx*ny*nz, ny, nx*ny   ,-1) ! y-direction
  call fft(a1,b1, nx*ny*nz, nz, nx*ny*nz,-1) ! z-direction
  
  ! Doing the y field
  if(lroot .AND. ip<10) print*,'doing fft of y-component'
  call fft(a2,b2, nx*ny*nz, nx, nx      ,-1) ! x-direction
  call fft(a2,b2, nx*ny*nz, ny, nx*ny   ,-1) ! y-direction
  call fft(a2,b2, nx*ny*nz, nz, nx*ny*nz,-1) ! z-direction
  
  ! Doing the z field
  if(lroot .AND. ip<10) print*,'doing fft of z-component'
  call fft(a3,b3, nx*ny*nz, nx, nx      ,-1) ! x-direction
  call fft(a3,b3, nx*ny*nz, ny, nx*ny   ,-1) ! y-direction
  call fft(a3,b3, nx*ny*nz, nz, nx*ny*nz,-1) ! z-direction

  ! Normalize
  a1=a1/nwgrid; a2=a2/nwgrid; a3=a3/nwgrid
  b1=b1/nwgrid; b2=b2/nwgrid; b3=b3/nwgrid

end subroutine transform
!***********************************************************************
subroutine transform_i(a_re,a_im)
!
!  Subroutine to do fourier transform
!  The routine overwrites the input data
!
!  22-oct-02/axel+tarek: adapted from transform
!
  real,dimension(nx,ny,nz) :: a_re,a_im
!
  if(lroot .AND. ip<10) print*,'doing three FFTs'
  call fft(a_re,a_im, nx*ny*nz, nx, nx      ,-1)
  call fft(a_re,a_im, nx*ny*nz, ny, nx*ny   ,-1)
  call fft(a_re,a_im, nx*ny*nz, nz, nx*ny*nz,-1)
!
!  Normalize
!
  a_re=a_re/nwgrid
  a_im=a_im/nwgrid
!
end subroutine transform_i
!***********************************************************************
subroutine transform_fftpack(a_re,a_im,direction)
!
!  Subroutine to do Fourier transform
!  The routine overwrites the input data
!
!  27-oct-02/axel: adapted from transform_i, for fftpack
!
  real,dimension(nx,ny,nz) :: a_re,a_im
  complex,dimension(nx) :: ax
  complex,dimension(ny) :: ay
  complex,dimension(nz) :: az
  real,dimension(4*nx+15) :: wsavex
  real,dimension(4*ny+15) :: wsavey
  real,dimension(4*nz+15) :: wsavez
  logical :: lforward=.true.
  integer,optional :: direction
  integer :: l,m,n

  if (present(direction).and. (direction.eq.-1)) lforward=.false.
!
  if(lroot .AND. ip<10) print*,'doing FFTpack in x, direction =',direction
  call cffti(nx,wsavex)
  do m=1,ny
  do n=1,nz
    ax=cmplx(a_re(:,m,n),a_im(:,m,n))
    if (lforward) then 
        call cfftf(nx,ax,wsavex)
    else 
        call cfftb(nx,ax,wsavex)
    endif
    a_re(:,m,n)=real(ax)
    a_im(:,m,n)=aimag(ax)
  enddo
  enddo
!
  if(lroot .AND. ip<10) print*,'doing FFTpack in y, direction =',direction
  call cffti(ny,wsavey)
  do l=1,nx
  do n=1,nz
    ay=cmplx(a_re(l,:,n),a_im(l,:,n))
    if (lforward) then 
        call cfftf(ny,ay,wsavey)
    else 
        call cfftb(ny,ay,wsavey)
    endif
    a_re(l,:,n)=real(ay)
    a_im(l,:,n)=aimag(ay)
  enddo
  enddo
!
  if(lroot .AND. ip<10) print*,'doing FFTpack in z, direction =',direction
  call cffti(nz,wsavez)
  do l=1,nx
  do m=1,ny
    az=cmplx(a_re(l,m,:),a_im(l,m,:))
    if (lforward) then 
       call cfftf(nz,az,wsavez)
    else 
       call cfftb(nz,az,wsavez)
    endif
    a_re(l,m,:)=real(az)
    a_im(l,m,:)=aimag(az)
  enddo
  enddo
!
!  Normalize
!

  if (lforward) then 
    a_re=a_re/nwgrid
    a_im=a_im/nwgrid
  endif
!
end subroutine transform_fftpack
!***********************************************************************
subroutine transform_nr(a_re,a_im)
!
!  Subroutine to do Fourier transform using Numerical Recipes routine.
!  Note that this routine requires that nx, ny, and nz are powers of 2.
!  The routine overwrites the input data
!
!  30-oct-02/axel: adapted from transform_fftpack for Numerical Recipes
!
  real,dimension(nx,ny,nz) :: a_re,a_im
  complex,dimension(nx) :: ax
  complex,dimension(ny) :: ay
  complex,dimension(nz) :: az
  integer :: l,m,n
!
!  This Fourier transform would work, but it's very slow!
!  Even the compilation is very slow, so we better get rid of it!  
!
  print*,'fft_nr currently disabled!'
  call stop_it("")
!
  if(lroot .AND. ip<10) print*,'doing FFT_nr in x'
  do m=1,ny
  do n=1,nz
    ax=cmplx(a_re(:,m,n),a_im(:,m,n))
    !call four1(ax,nx,-1)
    a_re(:,m,n)=real(ax)
    a_im(:,m,n)=aimag(ax)
  enddo
  enddo
!
  if(lroot .AND. ip<10) print*,'doing FFT_nr in y'
  do l=1,nx
  do n=1,nz
    ay=cmplx(a_re(l,:,n),a_im(l,:,n))
    !call four1(ay,ny,-1)
    a_re(l,:,n)=real(ay)
    a_im(l,:,n)=aimag(ay)
  enddo
  enddo
!
  if(lroot .AND. ip<10) print*,'doing FFT_nr in z'
  do l=1,nx
  do m=1,ny
    az=cmplx(a_re(l,m,:),a_im(l,m,:))
    !call four1(az,nz,-1)
    a_re(l,m,:)=real(az)
    a_im(l,m,:)=aimag(az)
  enddo
  enddo
!
!  Normalize
!
  a_re=a_re/nwgrid
  a_im=a_im/nwgrid
!
end subroutine transform_nr
!***********************************************************************
subroutine transform_fftpack_1d(a_re,a_im)
!
!  Subroutine to do Fourier transform
!  The routine overwrites the input data
!
!  06-feb-03/nils: adapted from transform_fftpack
!
  real,dimension(nx,ny,nz) :: a_re,a_im
  complex,dimension(nx) :: ax
  real,dimension(4*nx+15) :: wsavex
  integer :: m,n
!
  if(lroot .AND. ip<10) print*,'doing FFTpack in x'
  call cffti(nx,wsavex)
  do m=1,ny
  do n=1,nz
    ax=cmplx(a_re(:,m,n),a_im(:,m,n))
    call cfftf(nx,ax,wsavex)
    a_re(:,m,n)=real(ax)
    a_im(:,m,n)=aimag(ax)
  enddo
  enddo
!
!  Normalize
!
  a_re=a_re/nxgrid
  a_im=a_im/nxgrid
!
end subroutine transform_fftpack_1d
!***********************************************************************
endmodule Mpicomm
