!!!!!!!!!!!!!!!!!!!!!!!
!!!   io_dist.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Distributed IO (i.e. each process writes its own file tmp/procX)
!!!  07-Nov-2001/wd: Put into separate module, so one can choose
!!!  alternative IO mechanism.

module Io

  interface output              ! Overload the `output' function
    module procedure output_vect
    module procedure output_scal
  endinterface

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
        if (mx.eq.1) dx=dmax
        if (my.eq.1) dy=dmax
        if (mz.eq.1) dz=dmax
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
!!***********************************************************************
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
      use Cdata, only: x,y,z,dx,dy,dz
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
      character (LEN=*) :: file
!
      open(1,FILE=file,FORM='unformatted')
      read(1) t,x,y,z
      read(1) dx,dy,dz
!
    endsubroutine rgrid
!***********************************************************************

endmodule Io
