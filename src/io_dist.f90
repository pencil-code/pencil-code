!!!!!!!!!!!!!!!!!!!!!!!
!!!   io_dist.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Distributed IO (i.e. each process writes its own file tmp/procX)
!!!  07-Nov-2001/wd: Put into separate module, so one can choose
!!!  alternative IO mechanism.

module Io

  implicit none

  interface output              ! Overload the `output' function
    module procedure output_vect
    module procedure output_scal
  endinterface

  interface output_stenc        ! Overload the `output_stenc' function
    module procedure output_stenc_vect
    module procedure output_stenc_scal
  endinterface

  ! Interface to external C function
  ! Does not work, since stenc can be either a 1-d or a 2-d array and the
  ! C function does not care.
  !   interface output_stenciled_c
  !     subroutine output_stenciled_c(filename,stenc,&
  !                                   ndim,i,iy,iz,t, &
  !                                   nx,ny,nz,nghost,fnlen)
  !       use Cdata, only: mx
  
  !       real,dimension(mx,*) :: stenc
  !       real,dimension(mx) :: stenc
  !       real :: t
  !       integer :: ndim,i,iy,iz,nx,ny,nz,nghost,fnlen
  !       character (LEN=*) :: filename
  !     endsubroutine output_stenciled_c
  !   endinterface

  external output_stenciled_c   ! Note really needed, but self-documenting

contains

!***********************************************************************
    subroutine input(file,a,nn,mode)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!  11-apr-97/axel: coded
!
      use Cdata
!
      character (LEN=*) :: file
      integer :: nn,mode
      real, dimension (mx,my,mz,nn) :: a
!
      open(1,file=file,form='unformatted')
      if (ip<=8) print*,'open, mx,my,mz,nn=',mx,my,mz,nn
      read(1) a
      if (ip<=8) print*,'read ',file
      if (mode==1) then
        read(1) t,x,y,z
        if (ip<=3) print*,'ip,x',ip,x
        if (ip<=3) print*,'y',y
        if (ip<=3) print*,'z',z
!
!  assume uniform mesh; use the first two *interior* points
!  to calculate mesh spacing
!
        if (mx.gt.1) then; dx=x(5)-x(4); else; dx=0.; endif
        if (my.gt.1) then; dy=y(5)-y(4); else; dy=0.; endif
        if (mz.gt.1) then; dz=z(5)-z(4); else; dz=0.; endif
        dxmax=max(dx,dy,dz)
        if (mx.eq.1) dx=dxmax
        if (my.eq.1) dy=dxmax
        if (mz.eq.1) dz=dxmax
        dxmin=min(dx,dy,dz)
        Lx=dx*mx
        Ly=dy*my
        Lz=dz*mz
!
        if (ip<=4) print*
        if (ip<=4) print*,'dt,dx,dy,dz=',dt,dx,dy,dz
        if (ip<=8) print*,'pi,nu=',pi,nu
      endif
!
      close(1)
    endsubroutine input
!***********************************************************************
!    subroutine output(file,a,nn)
!!
!!  write snapshot file, always write mesh and time, could add other things
!!  11-apr-97/axel: coded
!!
!      use Cdata
!
!      integer :: nn
!      real, dimension (mx,my,mz,nn) :: a
!      character (LEN=*) :: file
!!
!      open(1,file=file,form='unformatted')
!      write(1) a
!      write(1) t,x,y,z
!      close(1)
!    endsubroutine output
!!***********************************************************************
!    subroutine output_mvarvect(file,a,nn)
!!
!!  write snapshot file, always write time and mesh, could add other things
!!  version for vector field
!!  11-apr-97/axel: coded
!!
!      use Cdata
!!
!      integer :: nn
!      real, dimension (mx,my,mz,mvar) :: a
!      character (LEN=*) :: file
!!
!!      print*,'OUTPUT_VECTOR'
!      if (nn /= mvar) STOP "OUTPUT_3vect called with nn/=mvar"
!      open(91,file=file,form='unformatted')
!      write(91) a
!      write(91) t,x,y,z
!      close(91)
!    endsubroutine output_mvarvect
!!***********************************************************************
!    subroutine output_3vect(file,a,nn)
!!
!!  write snapshot file, always write time and mesh, could add other things
!!  version for vector field
!!  11-apr-97/axel: coded
!!
!      use Cdata
!!
!      integer :: nn
!      real, dimension (mx,my,mz,3) :: a
!      character (LEN=*) :: file
!!
!!      print*,'OUTPUT_VECTOR'
!      if (nn /= 3) STOP "OUTPUT_3vect called with nn/=3"
!      open(91,file=file,form='unformatted')
!      write(91) a
!      write(91) t,x,y,z
!      close(91)
!    endsubroutine output_3vect
!***********************************************************************
    subroutine output_vect(file,a,nn)
!
!  write snapshot file, always write time and mesh, could add other things
!  version for vector field
!  11-apr-97/axel: coded
!
      use Cdata
!
      integer :: nn
      real, dimension (mx,my,mz,nn) :: a
      character (LEN=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'OUTPUT_VECTOR: nn =', nn
      open(91,file=file,form='unformatted')
      write(91) a
      write(91) t,x,y,z
      close(91)
    endsubroutine output_vect
!***********************************************************************
    subroutine output_scal(file,a,nn)
!
!  write snapshot file, always write time and mesh, could add other things
!  version for scalar field
!  11-apr-97/axel: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
!
      integer :: nn
      real, dimension (mx,my,mz) :: a
      character (LEN=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'OUTPUT_SCALAR'
      if (nn /= 1) call stop_it("OUTPUT called with scalar field, but nn/=1")
      open(91,file=file,form='unformatted')
      write(91) a
      write(91) t,x,y,z
      close(91)
    endsubroutine output_scal
!***********************************************************************
    subroutine output_stenc_vect(file,a,ndim)
!
!  Write snapshot file of stenciled vector data (for debugging).
!  Wrapper to the C routine output_stenciled_c.
!
!  15-feb-02/wolf: coded
!
      use Cdata
      use Mpicomm, only: mm,nn
!
      integer :: ndim,imn
      real, dimension (nx,ndim) :: a
      character (LEN=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'OUTPUT_STENC_VECT: nn =', nn
!
      call output_stenciled_c(file, a, ndim, &
                              imn, mm(imn), nn(imn), t, &
                              nx, ny, nz, nghost, len(file))
!
    endsubroutine output_stenc_vect
!***********************************************************************
    subroutine output_stenc_scal(file,a,ndim)
!
!  Write snapshot file of stenciled scalar data (for debugging).
!  Wrapper to the C routine output_stenciled_c.
!
!  15-feb-02/wolf: coded
!
      use Cdata
      use Mpicomm, only: imn,mm,nn,lroot,stop_it

!
      integer :: ndim
      real, dimension (nx) :: a
      character (LEN=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'OUTPUT_SCALAR'
      if (ndim /= 1) &
           call stop_it("OUTPUT called with scalar field, but ndim/=1")
!
      call output_stenciled_c(file, a, ndim, &
                              imn, mm(imn), nn(imn), t, &
                              nx, ny, nz, nghost, len(file))
!
    endsubroutine output_stenc_scal
!***********************************************************************
    subroutine outpus(file,a,nn)
!
!  write snapshot file, always write mesh and time, could add other things
!  11-oct-98/axel: adapted
!
      use Cdata
!
      integer :: nn
      character (LEN=*) :: file
      real, dimension (mx,my,mz,nn) :: a
!
!nn=mx/2
!l1=(mx-nn)/2+1; l2=l1+nn-1
!m1=(mx-nn)/2+1; m2=m1+nn-1
!n1=(mx-nn)/2+1; n2=n1+nn-1
      open(1,file=file,form='unformatted')
      write(1) a(l1:l2,m1:m2,n1:n2,:)
      write(1) t,x,y,z
      close(1)
    endsubroutine outpus
!***********************************************************************
    subroutine wgrid (file)
!
!  Write processor-local part of grid coordinates.
!  21-jan-02/wolf: coded
!
      use Cdata, only: t,x,y,z,dx,dy,dz
!      use Mpicomm
!
      character (LEN=*) :: file
!
      open(1,FILE=file,FORM='unformatted')
      write(1) t,x,y,z
      write(1) dx,dy,dz
!
    endsubroutine wgrid
!***********************************************************************
    subroutine rgrid (file)
!
!  Read processor-local part of grid coordinates.
!  21-jan-02/wolf: coded
!
      use Cdata, only: x,y,z,dx,dy,dz
!      use Mpicomm
!
      real :: tdummy
      character (LEN=*) :: file
!
      open(1,FILE=file,FORM='unformatted')
      read(1) tdummy,x,y,z
      read(1) dx,dy,dz
!
    endsubroutine rgrid
!***********************************************************************

endmodule Io
