! $Id: io_mpio.f90,v 1.1 2002-09-20 10:16:27 dobler Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   io_mpi-io.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  Parallel IO via MPI2 (i.e. all process write to the same file, e.g.
!!!  tmp/var.dat
!!!  19-sep-2002/wolf: started

module Io

  use Cparam

  implicit none
 
  interface output              ! Overload the `output' function
    module procedure output_vect
    module procedure output_scal
  endinterface

  interface output_pencil        ! Overload the `output_pencil' function
    module procedure output_pencil_vect
    module procedure output_pencil_scal
  endinterface

  include 'mpiof.h'

  integer, dimension(3) :: globalsize=(/nx,ny,nz/)
  integer, dimension(3) :: localsize=(/ nx/nprocx,ny/nprocy,nz/nprocz /)

  !
  ! Interface to external C function(s).
  ! Need to have two different C functions in order to have F90
  ! interfaces, since a pencil can be either a 1-d or a 2-d array.
  !
!   interface output_penciled_vect_c
!     subroutine output_penciled_vect_c(filename,pencil,&
!                                       ndim,i,iy,iz,t, &
!                                       nx,ny,nz,nghost,fnlen)
!       use Cdata, only: mx
!       real,dimension(mx,*) :: pencil
!       real :: t
!       integer :: ndim,i,iy,iz,nx,ny,nz,nghost,fnlen
!       character (len=*) :: filename
!     endsubroutine output_penciled_vect_c
!   endinterface
!   !
!   interface output_penciled_scal_c
!     subroutine output_penciled_scal_c(filename,pencil,&
!                                       ndim,i,iy,iz,t, &
!                                       nx,ny,nz,nghost,fnlen)
!       use Cdata, only: mx
!       real,dimension(mx) :: pencil
!       real :: t
!       integer :: ndim,i,iy,iz,nx,ny,nz,nghost,fnlen
!       character (len=*) :: filename
!     endsubroutine output_penciled_scal_c
!   endinterface
  !
  !  Still not possible with the NAG compiler (`No specific match for
  !  reference to generic OUTPUT_PENCILED_SCAL_C')
  !
  external output_penciled_scal_c
  external output_penciled_vect_c

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
!
!  check whether we want to read deltay from snapshot
!
        if (lshear) then
          read(1) t,x,y,z,dx,dy,dz,deltay
        else
          read(1) t,x,y,z,dx,dy,dz
        endif
!
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
      write(91) t,x,y,z,dx,dy,dz,deltay
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
      write(91) t,x,y,z,dx,dy,dz,deltay
      close(91)
    endsubroutine output_scal
!***********************************************************************
    subroutine output_scal_mpio(file,a,nn)
!
!  write snapshot file; currently without ghost zones and any meta data
!  like time, etc.
!  19-sep-02/wolf: coded
!
      use Cdata
!
      integer :: nn
      real, dimension (mx,my,mz) :: a
      character (len=*) :: file
      integer, dimension(3) :: start_index
      integer :: filetype,fhandle
      integer, dimension(MPI_STATUS_SIZE) :: status

      if ((ip<=8) .and. lroot) print*,'OUTPUT_SCALAR'
      if (nn /= 1) call stop_it("OUTPUT called with scalar field, but nn/=1")
      !
      !  global indices of first element of iproc's data in the file
      !
      start_index(1) = ipx*localsize(1)
      start_index(2) = ipy*localsize(2)
      start_index(3) = ipz*localsize(3)
      !
      !  construct the new datatype `filetype' 
      !
      call MPI_TYPE_CREATE_SUBARRAY(3, &
               globalsize, localsize, start_index, &
               MPI_ORDER_FORTRAN, MPI_REAL, &
               filetype, ierr)
      call MPI_TYPE_COMMIT(filetype,ierr)
      !
      !  open file and set view (specify which file positions we can access)
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, file, &
               ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
               MPI_INFO_NULL, fhandle)
      call MPI_FILE_SET_VIEW(fhandle, 0, MPI_REAL, filetype, "native", &
               MPI_INFO_NULL, ierr)
      !
      !  write data
      !
      call MPI_FILE_WRITE_ALL(fhandle, a, product(localsize), MPI_REAL, &
               status, ierr)
      call MPI_FILE_CLOSE(fhandle)

!      open(91,file=file,form='unformatted')
!      write(91) a
!      write(91) t,x,y,z,dx,dy,dz,deltay
!      close(91)

    endsubroutine output_scal_mpio
!***********************************************************************
    subroutine output_pencil_vect(file,a,ndim)
!
!  Write snapshot file of penciled vector data (for debugging).
!  Wrapper to the C routine output_penciled_vect_c.
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
      if (ip<9.and.lroot.and.imn==1) &
           print*,'output_pencil_vect('//file//'): ndim=',ndim
!
      if (headt .and. (imn==1)) print*, &
           'OUTPUT_PENCIL: Writing to ', trim(file), &
           ' for debugging -- this may slow things down'
!
       call output_penciled_vect_c(file, a, ndim, &
                                   imn, mm(imn), nn(imn), t, &
                                   nx, ny, nz, nghost, len(file))
!
    endsubroutine output_pencil_vect
!***********************************************************************
    subroutine output_pencil_scal(file,a,ndim)
!
!  Write snapshot file of penciled scalar data (for debugging).
!  Wrapper to the C routine output_penciled_scal_c.
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
      if ((ip<=8) .and. lroot .and. imn==1) &
           print*,'output_pencil_scal('//file//')'
!
      if (ndim /= 1) &
           call stop_it("OUTPUT called with scalar field, but ndim/=1")
!
      if (headt .and. (imn==1)) print*, &
           'OUTPUT_PENCIL: Writing to ', trim(file), &
           ' for debugging -- this may slow things down'
!
      call output_penciled_scal_c(file, a, ndim, &
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
      open(1,file=file,form='unformatted')
      write(1) a(l1:l2,m1:m2,n1:n2,:)
      write(1) t,x,y,z,dx,dy,dz,deltay
      close(1)
    endsubroutine outpus
!***********************************************************************
    subroutine wgrid (file)
!
!  Write processor-local part of grid coordinates.
!  21-jan-02/wolf: coded
!
      use Cdata, only: t,x,y,z,dx,dy,dz
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
