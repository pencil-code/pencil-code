! $Id$
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   debug_io_mpi.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Parallel debug-IO via MPI2 (i.e. write to a single file in data/allprocs)
!
!  The file format written by output() (and used, e.g. in var.dat)
!  consists of the followinig Fortran records:
!    1. data(mxgrid,mygrid,mzgrid,nvar)
!    2. t(1), x(mxgrid), y(mygrid), z(mzgrid), dx(1), dy(1), dz(1)
!    3. deltay(1) [optional: if lshear==.true.]
!  Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
!  for one vector field, 8 for var.dat in the case of MHD with entropy.
!
!  11-Dec-2011/Bourdin.KIS: moved debug-IO from 'io_mpio' into separate module
!
module Debug_IO
!
  use Cdata
!
  implicit none
!
  public :: lun_input, lun_output
  public :: output, output_pencil
!
  interface output              ! Overload the `output' function
    module procedure output_vect
    module procedure output_scal
  endinterface
!
  interface output_pencil        ! Overload the `output_pencil' function
    module procedure output_pencil_vect
    module procedure output_pencil_scal
  endinterface
!
  ! define unique logical unit number for input and output calls
  integer :: lun_input=89
  integer :: lun_output=92
!
  external output_penciled_scal_c
  external output_penciled_vect_c
!
  include 'mpif.h'
!
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: data_start=4
  integer :: io_filetype,io_memtype,io_filetype_v,io_memtype_v
  integer :: fhandle,ierr
  logical :: io_initialized=.false.
!
contains
!***********************************************************************
    subroutine commit_io_type_vect(nv)
!
!  For a new value of nv, commit MPI types needed for output_vect(). If
!  called with the same value of nv as the previous time, do nothing.
!  20-sep-02/wolf: coded
!
      integer, intent(in) :: nv
!
      integer, dimension(4) :: globalsize_v,localsize_v,memsize_v
      integer, dimension(4) :: start_index_v,mem_start_index_v
      integer,save :: lastnv=-1 ! value of nv at previous call
!
      if (nv /= lastnv) then
        if (lastnv > 0) then
          ! free old types, so we can re-use them
          call MPI_TYPE_FREE(io_filetype_v, ierr)
          call MPI_TYPE_FREE(io_memtype_v, ierr)
        endif
        globalsize_v=(/nxgrid,nygrid,nzgrid,nv/)
        localsize_v =(/nx    ,ny    ,nz    ,nv/)
        memsize_v   =(/mx    ,my    ,mz    ,nv/)
!
!  global indices of first element of iproc's data in the file
!
        start_index_v(1) = ipx*nx
        start_index_v(2) = ipy*ny
        start_index_v(3) = ipz*nz
        start_index_v(4) = 0
!
!  construct the new datatype `io_filetype_n'
!
        call MPI_TYPE_CREATE_SUBARRAY(4, &
                 globalsize_v, localsize_v, start_index_v, &
                 MPI_ORDER_FORTRAN, MPI_REAL, &
                 io_filetype_v, ierr)
        call MPI_TYPE_COMMIT(io_filetype_v,ierr)
!
!  create a derived datatype describing the data layout in memory
!  (i.e. including the ghost zones)
!
        mem_start_index_v(1:3) = nghost
        mem_start_index_v(4)   = 0
        call MPI_TYPE_CREATE_SUBARRAY(4, &
                 memsize_v, localsize_v, mem_start_index_v, &
                 MPI_ORDER_FORTRAN, MPI_REAL, &
                 io_memtype_v, ierr)
        call MPI_TYPE_COMMIT(io_memtype_v,ierr)
      endif
!
    endsubroutine commit_io_type_vect
!***********************************************************************
    subroutine commit_io_type_vect_1D(nr,nv)
!
!  For a new value of nv, commit MPI types needed for output_vect(). If
!  called with the same value of nv as the previous time, do nothing.
!  20-sep-02/wolf: coded
!
      integer, intent(in) :: nr,nv
!
      integer, dimension(2) :: globalsize_v,localsize_v,memsize_v
      integer, dimension(2) :: start_index_v,mem_start_index_v
      integer,save :: lastnv=-1 ! value of nv at previous call
!
      if (nv /= lastnv) then
        if (lastnv > 0) then
          ! free old types, so we can re-use them
          call MPI_TYPE_FREE(io_filetype_v, ierr)
          call MPI_TYPE_FREE(io_memtype_v, ierr)
        endif
        globalsize_v=(/nr,nv/)
        localsize_v =(/nr,nv/)
        memsize_v   =(/nr,nv/)
!
!  global indices of first element of iproc's data in the file
!
        start_index_v(1) = ipx*nr
        start_index_v(2) = 0
!
!  construct the new datatype `io_filetype_n'
!
        call MPI_TYPE_CREATE_SUBARRAY(2, &
                 globalsize_v, localsize_v, start_index_v, &
                 MPI_ORDER_FORTRAN, MPI_REAL, &
                 io_filetype_v, ierr)
        call MPI_TYPE_COMMIT(io_filetype_v,ierr)
!
!  create a derived datatype describing the data layout in memory
!  (i.e. including the ghost zones)
!
        mem_start_index_v(1:2) = 0
        call MPI_TYPE_CREATE_SUBARRAY(2, &
                 memsize_v, localsize_v, mem_start_index_v, &
                 MPI_ORDER_FORTRAN, MPI_REAL, &
                 io_memtype_v, ierr)
        call MPI_TYPE_COMMIT(io_memtype_v,ierr)
      endif
!
    endsubroutine commit_io_type_vect_1D
!***********************************************************************
    subroutine output_vect(file,a,nv)
!
!  Write snapshot file; currently without ghost zones and grid.
!    Looks like we need to commit the MPI type anew each time we are called,
!  since nv may vary.
!
!  20-sep-02/wolf: coded
!
      use Mpicomm, only: stop_it
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: a
!
      character (len=fnlen) :: filename
!
      if ((ip<=8) .and. lroot) print*,'output_vect: nv =', nv
      if (.not. io_initialized) &
           call stop_it("output_vect: Need to call register_io first")
      filename = trim (datadir_snap)//'/'//trim (file)
!
      call commit_io_type_vect(nv) ! will free old type if new one is needed
      !
      !  open file and set view (specify which file positions we can access)
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
               ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
               MPI_INFO_NULL, fhandle, ierr)
      call MPI_FILE_SET_VIEW(fhandle, data_start, MPI_REAL, io_filetype_v, &
               "native", MPI_INFO_NULL, ierr)
      !
      !  write data
      !
      call MPI_FILE_WRITE_ALL(fhandle, a, 1, io_memtype_v, status, ierr)
      call MPI_FILE_CLOSE(fhandle, ierr)
      !
      !  write meta data (to make var.dat as identical as possible to
      !  what a single-processor job would write with io_dist.f90)
      !
      call write_record_info(filename, nv)
!
    endsubroutine output_vect
!***********************************************************************
    subroutine output_scal(file,a,nv)
!
!  Write snapshot file; currently without ghost zones and grid
!
!  20-sep-02/wolf: coded
!
      use Mpicomm, only: stop_it
!
      character (len=*), intent(in) :: file
      real, dimension (mx,my,mz), intent(in) :: a
      integer, intent(in) :: nv
!
      character (len=fnlen) :: filename
!
      if ((ip<=8) .and. lroot) print*,'output_scal: ENTER'
      if (.not. io_initialized) &
           call stop_it("output_scal: Need to call register_io first")
      if (nv /= 1) call stop_it("output_scal: called with scalar field, but nv/=1")
      filename = trim (datadir_snap)//'/'//trim (file)
      !
      !  open file and set view (specify which file positions we can access)
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
               ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
               MPI_INFO_NULL, fhandle, ierr)
      call MPI_FILE_SET_VIEW(fhandle, data_start, MPI_REAL, io_filetype, &
               "native", MPI_INFO_NULL, ierr)
      !
      !  write data
      !
      call MPI_FILE_WRITE_ALL(fhandle, a, 1, io_memtype, status, ierr)
      call MPI_FILE_CLOSE(fhandle, ierr)
      !
      !  write meta data (to make var.dat as identical as possible to
      !  what a single-processor job would write with io_dist.f90)
      !
      call write_record_info(filename, 1)
!
    endsubroutine output_scal
!***********************************************************************
    subroutine output_pencil_vect(file,a,ndim)
!
!  Write snapshot file of penciled vector data (for debugging).
!  Wrapper to the C routine output_penciled_vect_c.
!
!  15-feb-02/wolf: coded
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: ndim
      real, dimension (nx,ndim), intent(in) :: a
!
      real :: t_sp   ! t in single precision for backwards compatibility
      character (len=fnlen) :: filename
!
      t_sp = t
      if (ip<9.and.lroot.and.imn==1) &
           print*,'output_pencil_vect('//file//'): ndim=',ndim
      filename = trim (datadir_snap)//'/'//trim (file)
!
      if (headt .and. (imn==1)) write(*,'(A)') &
           'output_pencil_vect: Writing to ' // trim(filename) // &
           ' for debugging -- this may slow things down'
!
       call output_penciled_vect_c(filename, a, ndim, &
                                   imn, mm(imn), nn(imn), t_sp, &
                                   nx, ny, nz, nghost, len(filename))
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
      use Mpicomm, only: stop_it
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: ndim
      real, dimension (nx), intent(in) :: a
!
      real :: t_sp   ! t in single precision for backwards compatibility
      character (len=fnlen) :: filename
!
      t_sp = t
      if ((ip<=8) .and. lroot .and. imn==1) &
           print*,'output_pencil_scal('//file//')'
      filename = trim (datadir_snap)//'/'//trim (file)
!
      if (ndim /= 1) &
           call stop_it("OUTPUT called with scalar field, but ndim/=1")
!
      if (headt .and. (imn==1)) print*, &
           'output_pencil_scal: Writing to ', trim(filename), &
           ' for debugging -- this may slow things down'
!
      call output_penciled_scal_c(filename, a, ndim, &
                                  imn, mm(imn), nn(imn), t_sp, &
                                  nx, ny, nz, nghost, len(filename))
!
    endsubroutine output_pencil_scal
!***********************************************************************
    subroutine write_record_info(file, nv)
!
!  Add record markers and time to file, so it looks as similar as
!  possible/necessary to a file written by io_dist.f90. Currently, we
!  don't handle writing the grid, but that does not seem to be used
!  anyway.
!
      use Mpicomm, only: stop_it
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv
!
      integer(KIND=ikind8) :: reclen
      integer(kind=MPI_OFFSET_KIND) :: fpos
      real :: t_sp   ! t in single precision for backwards compatibility
      integer, parameter :: test_int = 0
      real, parameter :: test_float = 0
      integer :: byte_per_int, byte_per_float
!
      t_sp = t
      inquire (IOLENGTH=byte_per_int) test_int
      inquire (IOLENGTH=byte_per_float) test_float
!
      call MPI_FILE_OPEN(MPI_COMM_WORLD, file, & ! MPI_FILE_OPEN is collective
                         ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
                         MPI_INFO_NULL, fhandle, ierr)
      if (lroot) then           ! only root writes
        !
        ! record markers for (already written) data block
        !
        fpos = 0                ! open-record marker
        reclen = nxgrid*nygrid*nzgrid*nv*byte_per_float
        if (reclen>max_int) call stop_it('write_record_info: reclen>max_int')
        call MPI_FILE_WRITE_AT(fhandle,fpos,int(reclen),1,MPI_INTEGER,status,ierr)
        fpos = fpos + byte_per_int
                                ! the data itself has already been written by ouput_vect
        fpos = fpos + reclen    ! close-record marker
        call MPI_FILE_WRITE_AT(fhandle,fpos,int(reclen),1,MPI_INTEGER,status,ierr)
        fpos = fpos + byte_per_int
        !
        ! time in a new record
        !
        call MPI_FILE_WRITE_AT(fhandle,fpos,byte_per_float,1,MPI_INTEGER,status,ierr)
        fpos = fpos + byte_per_int
        call MPI_FILE_WRITE_AT(fhandle,fpos,t_sp,1,MPI_REAL,status,ierr)
        fpos = fpos + byte_per_float
        call MPI_FILE_WRITE_AT(fhandle,fpos,byte_per_float,1,MPI_INTEGER,status,ierr)
        fpos = fpos + byte_per_int
      endif
      !
      call MPI_FILE_CLOSE(fhandle,ierr)
!
    endsubroutine write_record_info
!***********************************************************************
    subroutine commit_gridio_types(nglobal,nlocal,mlocal,ipvar,filetype,memtype)
!
!  define MPI data types needed for MPI-IO of grid files
!  20-sep-02/wolf: coded
!
      integer :: nglobal,nlocal,mlocal,ipvar
      integer :: filetype,memtype
!
!  construct new datatype `filetype'
!
      call MPI_TYPE_CREATE_SUBARRAY(1, &
               nglobal, nlocal, ipvar*nlocal, &
               MPI_ORDER_FORTRAN, MPI_REAL, &
               filetype, ierr)
      call MPI_TYPE_COMMIT(filetype,ierr)
!
!  create a derived datatype describing the data layout in memory
!  (i.e. including the ghost zones)
!
      call MPI_TYPE_CREATE_SUBARRAY(1, &
               mlocal, nlocal, nghost, &
               MPI_ORDER_FORTRAN, MPI_REAL, &
               memtype, ierr)
      call MPI_TYPE_COMMIT(memtype,ierr)
!
    endsubroutine commit_gridio_types
!***********************************************************************
    subroutine write_grid_data(file,nglobal,nlocal,mlocal,ipvar,var)
!
!  write grid positions for var to file
!  20-sep-02/wolf: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension(*) :: var ! x, y or z
      integer :: nglobal,nlocal,mlocal,ipvar
      integer :: filetype,memtype
      character (len=*) :: file
!
      call commit_gridio_types(nglobal,nlocal,mlocal,ipvar,filetype,memtype)
!
!  open file and set view (specify which file positions we can access)
!
      call MPI_FILE_OPEN(MPI_COMM_WORLD, file, &
               ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
               MPI_INFO_NULL, fhandle, ierr)
!
      if (ierr /= 0) call stop_it( &
           "Cannot MPI_FILE_OPEN " // trim(file) // &
           " (or similar) for writing" // &
           " -- is data/ visible from all nodes?")
!
      call MPI_FILE_SET_VIEW(fhandle, data_start, MPI_REAL, filetype, &
               "native", MPI_INFO_NULL, ierr)
!
!  write data and free type handle
!
      call MPI_FILE_WRITE_ALL(fhandle, var, 1, memtype, status, ierr)
      call MPI_FILE_CLOSE(fhandle, ierr)
!
      call MPI_TYPE_FREE(filetype, ierr)
      call MPI_TYPE_FREE(memtype, ierr)
!
    endsubroutine write_grid_data
!***********************************************************************
    subroutine read_grid_data(file,nglobal,nlocal,mlocal,ipvar,var)
!
!  read grid positions for var to file
!  20-sep-02/wolf: coded
!
      real, dimension(*) :: var ! x, y or z
      integer :: nglobal,nlocal,mlocal,ipvar
      integer :: filetype,memtype
      character (len=*) :: file
!
      call commit_gridio_types(nglobal,nlocal,mlocal,ipvar,filetype,memtype)
!
!  open file and set view (specify which file positions we can access)
!
      call MPI_FILE_OPEN(MPI_COMM_WORLD, file, &
               MPI_MODE_RDONLY, &
               MPI_INFO_NULL, fhandle, ierr)
      call MPI_FILE_SET_VIEW(fhandle, data_start, MPI_REAL, filetype, &
               "native", MPI_INFO_NULL, ierr)
!
!  read data and free type handle
!
      call MPI_FILE_READ_ALL(fhandle, var, 1, memtype, status, ierr)
      call MPI_FILE_CLOSE(fhandle, ierr)
!
      call MPI_TYPE_FREE(filetype, ierr)
      call MPI_TYPE_FREE(memtype, ierr)
!
    endsubroutine read_grid_data
!***********************************************************************
    subroutine wgrid(file)
!
!  Write processor-local part of grid coordinates.
!  21-jan-02/wolf: coded
!
      use Cdata, only: directory,dx,dy,dz
!
      character (len=*) :: file ! not used
      integer :: iostat
!
      call write_grid_data(trim(directory)//'/gridx.dat',nxgrid,nx,mx,ipx,x)
      call write_grid_data(trim(directory)//'/gridy.dat',nygrid,ny,my,ipy,y)
      call write_grid_data(trim(directory)//'/gridz.dat',nzgrid,nz,mz,ipz,z)
!
! write dx,dy,dz to their own file
!
      if (lroot) then
        open(lun_output,FILE=trim(directory)//'/dxyz.dat',FORM='unformatted',IOSTAT=iostat)
        if (outlog(iostat,'open',trim(directory)//'/dxyz.dat')) return
!
        write(lun_output,IOSTAT=iostat) dx,dy,dz
        if (outlog(iostat,'write dx,dy,dz')) return
!
        close(lun_output,IOSTAT=iostat)
        if (outlog(iostat,'close')) continue
      endif
!
      call keep_compiler_quiet(file)
!
    endsubroutine wgrid
!***********************************************************************
    subroutine rgrid(file)
!
!  Read processor-local part of grid coordinates.
!  21-jan-02/wolf: coded
!
      use Cdata, only: directory,dx,dy,dz,ip
!
      integer :: i,iostat
      character (len=*) :: file ! not used
!
      call read_grid_data(trim(directory)//'/gridx.dat',nxgrid,nx,mx,ipx,x)
      call read_grid_data(trim(directory)//'/gridy.dat',nygrid,ny,my,ipy,y)
      call read_grid_data(trim(directory)//'/gridz.dat',nzgrid,nz,mz,ipz,z)
!
!  read dx,dy,dz from file
!  We can't just compute them, since different
!  procs may obtain different results at machine precision, which causes
!  nasty communication failures with shearing sheets.
!
      open(lun_input,FILE=trim(directory)//'/dxyz.dat',FORM='unformatted',IOSTAT=iostat)
      if (iostat /= 0) call stop_it("Cannot open "//trim(directory)//'/dxyz.dat'//" for reading",iostat)
!
      read(lun_input,IOSTAT=iostat) dx,dy,dz
      if (iostat /= 0) call stop_it("Cannot read dx,dy,dz from "//trim(directory)//'/dxyz.dat',iostat)
!
      close(lun_input,IOSTAT=iostat)
      if (outlog(iostat,'close',(directory)//'/dxyz.dat')) continue
!
!  reconstruct ghost values
!
      do i=1,nghost
        x(l1-i) = x(l1) - i*dx
        y(m1-i) = y(m1) - i*dy
        z(n1-i) = z(n1) - i*dz
        !
        x(l2+i) = x(l2) + i*dx
        y(m2+i) = y(m2) + i*dy
        z(n2+i) = z(n2) + i*dz
      enddo
!
!  inherit Lx, Ly, Lz from start, and assume uniform mesh
!
      Lx=dx*nx*nprocx
      Ly=dy*ny*nprocy
      Lz=dz*nz*nprocz
!
      if (ip<=4) print*,'rgrid: dt,dx,dy,dz=',dt,dx,dy,dz
!
      call keep_compiler_quiet(file)
!
    endsubroutine rgrid
!***********************************************************************
endmodule Io
