! $Id: io_dist.f90,v 1.22 2002-06-01 02:56:21 brandenb Exp $

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

  interface output_pencil        ! Overload the `output_pencil' function
    module procedure output_pencil_vect
    module procedure output_pencil_scal
  endinterface

  ! Interface to external C function
  ! Does not work, since a pencil can be either a 1-d or a 2-d array and the
  ! C function does not care.
  !   interface output_penciled_c
  !     subroutine output_penciled_c(filename,pencil,&
  !                                   ndim,i,iy,iz,t, &
  !                                   nx,ny,nz,nghost,fnlen)
  !       use Cdata, only: mx
  
  !       real,dimension(mx,*) :: pencil
  !       real,dimension(mx) :: pencil
  !       real :: t
  !       integer :: ndim,i,iy,iz,nx,ny,nz,nghost,fnlen
  !       character (len=*) :: filename
  !     endsubroutine output_penciled_c
  !   endinterface

  external output_penciled_c   ! Note really needed, but self-documenting

contains

!***********************************************************************
    subroutine input(file,a,nn,mode)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!  11-apr-97/axel: coded
!
      use Cdata
!
      character (len=*) :: file
      integer :: nn,mode
      real, dimension (mx,my,mz,nn) :: a
!
      open(1,file=file,form='unformatted')
      if (ip<=8) print*,'open, mx,my,mz,nn=',mx,my,mz,nn
      read(1) a
      if (ip<=8) print*,'read ',file
      if (mode==1) then
        read(1) t,x,y,z,dx,dy,dz
        if (ip<=3) print*,'ip,x',ip,x
        if (ip<=3) print*,'y',y
        if (ip<=3) print*,'z',z
!
!  assume uniform mesh; use the first two *interior* points
!  to calculate mesh spacing
!
        dxmax=max(dx,dy,dz)
        dxmin=min(dx,dy,dz)
        Lx=dx*nx*nprocx
        Ly=dy*ny*nprocy
        Lz=dz*nz*nprocz
!
        if (ip<=4) print*
        if (ip<=4) print*,'dt,dx,dy,dz=',dt,dx,dy,dz
        if (ip<=8) print*,'pi,nu=',pi,nu
      endif
!
      close(1)
    endsubroutine input
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
      character (len=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'OUTPUT_VECTOR: nn =', nn
      open(91,file=file,form='unformatted')
      write(91) a
      write(91) t,x,y,z,dx,dy,dz
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
      character (len=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'OUTPUT_SCALAR'
      if (nn /= 1) call stop_it("OUTPUT called with scalar field, but nn/=1")
      open(91,file=file,form='unformatted')
      write(91) a
      write(91) t,x,y,z,dx,dy,dz
      close(91)
    endsubroutine output_scal
!***********************************************************************
    subroutine output_pencil_vect(file,a,ndim)
!
!  Write snapshot file of penciled vector data (for debugging).
!  Wrapper to the C routine output_penciled_c.
!
!  15-feb-02/wolf: coded
!
      use Cdata
      use Mpicomm, only: imn,mm,nn
!
      integer :: ndim
      real, dimension (nx,ndim) :: a
      character (len=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'OUTPUT_PENCIL_VECT: ndim =', ndim
!
       call output_penciled_c(file, a, ndim, &
                               imn, mm(imn), nn(imn), t, &
                               nx, ny, nz, nghost, len(file))
!
    endsubroutine output_pencil_vect
!***********************************************************************
    subroutine output_pencil_scal(file,a,ndim)
!
!  Write snapshot file of penciled scalar data (for debugging).
!  Wrapper to the C routine output_penciled_c.
!
!  15-feb-02/wolf: coded
!
      use Cdata
      use Mpicomm, only: imn,mm,nn,lroot,stop_it

!
      integer :: ndim
      real, dimension (nx) :: a
      character (len=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'OUTPUT_PENCIL_SCAL'
      if (ndim /= 1) &
           call stop_it("OUTPUT called with scalar field, but ndim/=1")
!
      call output_penciled_c(file, a, ndim, &
                              imn, mm(imn), nn(imn), t, &
                              nx, ny, nz, nghost, len(file))
!
    endsubroutine output_pencil_scal
!***********************************************************************
    subroutine outpus(file,a,nn)
!
!  write snapshot file, always write mesh and time, could add other things
!  11-oct-98/axel: adapted
!
      use Cdata
!
      integer :: nn
      character (len=*) :: file
      real, dimension (mx,my,mz,nn) :: a
!
!nn=mx/2
!l1=(mx-nn)/2+1; l2=l1+nn-1
!m1=(mx-nn)/2+1; m2=m1+nn-1
!n1=(mx-nn)/2+1; n2=n1+nn-1
      open(1,file=file,form='unformatted')
      write(1) a(l1:l2,m1:m2,n1:n2,:)
      write(1) t,x,y,z,dx,dy,dz
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
      character (len=*) :: file
!
      open(1,FILE=file,FORM='unformatted')
      write(1) t,x,y,z,dx,dy,dz
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
      character (len=*) :: file
!
      open(1,FILE=file,FORM='unformatted')
      read(1) tdummy,x,y,z,dx,dy,dz
      read(1) dx,dy,dz
!
    endsubroutine rgrid
!***********************************************************************

endmodule Io
