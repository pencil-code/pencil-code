! $Id$
!
!!!!!!!!!!!!!!!!!!!!!!!
!!!   io_mpio.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!
!
!  Parallel IO via MPI2 (i.e. all process write to the same file, like
!  data/allprocs/var.dat)
!
!  19-sep-02/wolf: started
!
!  The file format written by output() (and used, e.g. in var.dat)
!  consists of the followinig Fortran records:
!    1. data(mx,my,mz,nvar)
!    2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
!  Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
!  for one vector field, 8 for var.dat in the case of MHD with entropy.
!
module Io
!
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'io.h'
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
  integer :: lun_input=88
  integer :: lun_output=91
!
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
!       double precision :: t
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
!       double precision :: t
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
!
  include 'mpif.h'
!  include 'mpiof.h'
!
  integer, dimension(MPI_STATUS_SIZE) :: status
! LAM-MPI does not know the MPI2 constant MPI_OFFSET_KIND, but LAM
! doesn't work with this module anyway
!  integer(kind=8) :: data_start=4
  integer(kind=MPI_OFFSET_KIND) :: data_start=4
  integer :: io_filetype,io_memtype,io_filetype_v,io_memtype_v
  integer :: fhandle,ierr
  logical :: io_initialized=.false.
!
contains
!
!***********************************************************************
    subroutine register_io()
!
!  define MPI data types needed for MPI-IO
!  closely follwing Gropp et al. `Using MPI-2'
!
!  20-sep-02/wolf: coded
!  25-oct-02/axel: removed assignment of datadir; now set in cdata.f90
!
      use Cdata, only: datadir,directory_snap
      use Mpicomm, only: stop_it
!
      integer, dimension(3) :: globalsize=(/nxgrid,nygrid,nzgrid/)
      integer, dimension(3) :: localsize =(/nx    ,ny    ,nz    /)
      integer, dimension(3) :: memsize   =(/mx    ,my    ,mz    /)
      integer, dimension(3) :: start_index,mem_start_index
!
      lmonolithic_io = .true.   ! we write f to one single file
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  consistency check
!
      if (.not. lmpicomm) call stop_it('Need mpicomm for io_mpio')
!
!  global indices of first element of iproc's data in the file
!
      start_index(1) = ipx*localsize(1)
      start_index(2) = ipy*localsize(2)
      start_index(3) = ipz*localsize(3)
!
!  construct the new datatype `io_filetype'
!
      call MPI_TYPE_CREATE_SUBARRAY(3, &
               globalsize, localsize, start_index, &
               MPI_ORDER_FORTRAN, MPI_REAL, &
               io_filetype, ierr)
      call MPI_TYPE_COMMIT(io_filetype,ierr)
!
!  create a derived datatype describing the data layout in memory
!  (i.e. including the ghost zones)
!
      mem_start_index = nghost
      call MPI_TYPE_CREATE_SUBARRAY(3, &
               memsize, localsize, mem_start_index, &
               MPI_ORDER_FORTRAN, MPI_REAL, &
               io_memtype, ierr)
      call MPI_TYPE_COMMIT(io_memtype,ierr)
!
      io_initialized=.true.
!
!  initialize datadir and directory_snap (where var.dat and VAR# go)
!  -- may be overwritten in *.in parameter file
!
      directory_snap = ''
!
    endsubroutine register_io
!***********************************************************************
    subroutine directory_names()
!
!  Set up the directory names:
!  initialize datadir to `data' (would have been `tmp' with older runs)
!  set directory name for the output (one subdirectory for each processor)
!  if  datadir_snap (where var.dat, VAR# go) is empty, initialize to datadir
!
!  02-oct-2002/wolf: coded
!
      use Cdata, only: datadir,directory,datadir_snap,directory_snap
      use General, only: safe_character_assign
!
      if ((datadir_snap == '') .or. (index(datadir_snap,'allprocs')>0)) then
        datadir_snap = datadir
      endif
!
      call safe_character_assign(directory, trim(datadir)//'/allprocs')
      call safe_character_assign(directory_snap,trim(datadir_snap)//'/allprocs')
!
    endsubroutine directory_names
!***********************************************************************
    subroutine commit_io_type_vect(nv)
!
!  For a new value of nv, commit MPI types needed for output_vect(). If
!  called with the same value of nv as the previous time, do nothing.
!  20-sep-02/wolf: coded
!
      integer, dimension(4) :: globalsize_v,localsize_v,memsize_v
      integer, dimension(4) :: start_index_v,mem_start_index_v
      integer,save :: lastnv=-1 ! value of nv at previous call
      integer :: nv
!
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
      integer :: nr,nv
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
    subroutine input(file,a,nv,mode)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!  11-apr-97/axel: coded
!
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
      integer :: nv,mode
      real, dimension (mx,my,mz,nv) :: a
!
      if (ip<=8) print*,'input: mx,my,mz,nv=',mx,my,mz,nv
      if (.not. io_initialized) &
           call stop_it("input: Need to call register_io first")
!
      call commit_io_type_vect(nv)
!
!  open file and set view (specify which file positions we can access)
!
      call MPI_FILE_OPEN(MPI_COMM_WORLD, file, &
               MPI_MODE_RDONLY, &
               MPI_INFO_NULL, fhandle, ierr)
      call MPI_FILE_SET_VIEW(fhandle, data_start, MPI_REAL, io_filetype_v, &
               "native", MPI_INFO_NULL, ierr)
!
!  read data
!
      call MPI_FILE_READ_ALL(fhandle, a, 1, io_memtype_v, status, ierr)
      call MPI_FILE_CLOSE(fhandle, ierr)
!
      call keep_compiler_quiet(mode)
!
    endsubroutine input
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
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'output_vect: nv =', nv
      if (.not. io_initialized) &
           call stop_it("output_vect: Need to call register_io first")
!
      call commit_io_type_vect(nv) ! will free old type if new one is needed
      !
      !  open file and set view (specify which file positions we can access)
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, file, &
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
      call write_record_info(file, nv)
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
      real, dimension (mx,my,mz) :: a
      integer :: nv
      character (len=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'output_scal: ENTER'
      if (.not. io_initialized) &
           call stop_it("output_scal: Need to call register_io first")
      if (nv /= 1) call stop_it("output_scal: called with scalar field, but nv/=1")
      !
      !  open file and set view (specify which file positions we can access)
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, file, &
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
      call write_record_info(file, 1)
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
      integer :: ndim
      real, dimension (nx,ndim) :: a
      character (len=*) :: file
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = t
      if (ip<9.and.lroot.and.imn==1) &
           print*,'output_pencil_vect('//file//'): ndim=',ndim
!
      if (headt .and. (imn==1)) write(*,'(A)') &
           'output_pencil_vect: Writing to ' // trim(file) // &
           ' for debugging -- this may slow things down'
!
       call output_penciled_vect_c(file, a, ndim, &
                                   imn, mm(imn), nn(imn), t_sp, &
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
      use Mpicomm, only: stop_it
!
!
      integer :: ndim
      real, dimension (nx) :: a
      character (len=*) :: file
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = t
      if ((ip<=8) .and. lroot .and. imn==1) &
           print*,'output_pencil_scal('//file//')'
!
      if (ndim /= 1) &
           call stop_it("OUTPUT called with scalar field, but ndim/=1")
!
      if (headt .and. (imn==1)) print*, &
           'output_pencil_scal: Writing to ', trim(file), &
           ' for debugging -- this may slow things down'
!
      call output_penciled_scal_c(file, a, ndim, &
                                  imn, mm(imn), nn(imn), t_sp, &
                                  nx, ny, nz, nghost, len(file))
!
    endsubroutine output_pencil_scal
!***********************************************************************
    subroutine outpus(file,a,nv)
!
!  write snapshot file, always write mesh and time, could add other things
!  11-oct-98/axel: adapted
!
      use Mpicomm, only: stop_it
!
      integer :: nv
      character (len=*) :: file
      real, dimension (mx,my,mz,nv) :: a
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = t
      call stop_it("outpus: doesn't work with io_mpio yet -- but wasn't used anyway")
!
      open(1,FILE=file,FORM='unformatted')
      write(1) a(l1:l2,m1:m2,n1:n2,:)
      write(1) t_sp,x,y,z,dx,dy,dz,deltay
      close(1)
    endsubroutine outpus
!***********************************************************************
    subroutine write_record_info(file, nv)
!
!  Add record markers and time to file, so it looks as similar as
!  possible/necessary to a file written by io_dist.f90. Currently, we
!  don't handle writing the grid, but that does not seem to be used
!  anyway.
!
      integer :: nv,reclen
      integer(kind=MPI_OFFSET_KIND) :: fpos
      character (len=*) :: file
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
        call MPI_FILE_WRITE_AT(fhandle,fpos,reclen,1,MPI_INTEGER,status,ierr)
        fpos = fpos + byte_per_int
                                ! the data itself has already been written by ouput_vect
        fpos = fpos + reclen    ! close-record marker
        call MPI_FILE_WRITE_AT(fhandle,fpos,reclen,1,MPI_INTEGER,status,ierr)
        fpos = fpos + byte_per_int
        !
        ! time in a new record
        !
        reclen = byte_per_float
        call MPI_FILE_WRITE_AT(fhandle,fpos,reclen,1,MPI_INTEGER,status,ierr)
        fpos = fpos + byte_per_int
        call MPI_FILE_WRITE_AT(fhandle,fpos,t_sp,1,MPI_REAL,status,ierr)
        fpos = fpos + reclen
        call MPI_FILE_WRITE_AT(fhandle,fpos,reclen,1,MPI_INTEGER,status,ierr)
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
    subroutine log_filename_to_file(filename,flist)
!
!  In the directory containing `filename', append one line to file
!  `flist' containing the file part of filename
!
!  13-may-08/wolf: adapted from io_dist.f90
!
      use Cdata, only: lroot,lcopysnapshots_exp,datadir
      use Cparam, only: fnlen
      use General, only: parse_filename
!
      character (len=*) :: filename,flist
      character (len=fnlen) :: dir,fpart
!
      if (lroot) then
        call parse_filename(filename,dir,fpart)
        open(lun_output,FILE=trim(dir)//'/'//trim(flist),POSITION='append')
        write(lun_output,'(A)') trim(fpart)
        close(lun_output)
      endif
!
    endsubroutine log_filename_to_file
!***********************************************************************
    subroutine wgrid(file)
!
!  Write processor-local part of grid coordinates.
!  21-jan-02/wolf: coded
!
      use Cdata, only: directory,dx,dy,dz
!
      character (len=*) :: file ! not used
!
      call write_grid_data(trim(directory)//'/gridx.dat',nxgrid,nx,mx,ipx,x)
      call write_grid_data(trim(directory)//'/gridy.dat',nygrid,ny,my,ipy,y)
      call write_grid_data(trim(directory)//'/gridz.dat',nzgrid,nz,mz,ipz,z)
!
! write dx,dy,dz to their own file
!
      if (lroot) then
        open(lun_output,FILE=trim(directory)//'/dxyz.dat',FORM='unformatted')
        write(lun_output) dx,dy,dz
        close(lun_output)
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
      use Cdata, only: directory,dx,dy,dz
!
      integer :: i
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
      open(lun_input,FILE=trim(directory)//'/dxyz.dat',FORM='unformatted')
      read(lun_input) dx,dy,dz
      close(lun_input)
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
!  Find minimum/maximum grid spacing. Note that
!    minval( (/dx,dy,dz/), MASK=((/nxgrid,nygrid,nzgrid/) > 1) )
!  will be undefined if all n[x-z]grid=1, so we have to add the fourth
!  component with a test that is always true
!
      dxmin = minval( (/dx,dy,dz,huge(dx)/), &
                MASK=((/nxgrid,nygrid,nzgrid,2/) > 1) )
      dxmax = maxval( (/dx,dy,dz,epsilon(dx)/), &
                MASK=((/nxgrid,nygrid,nzgrid,2/) > 1) )
!
!  Fill pencil with maximum gridspacing. Will be overwritten
!  during the mn loop in the non equiditant case
!
      dxmax_pencil(:) = dxmax
      dxmin_pencil(:) = dxmin
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
    subroutine wproc_bounds(file)
!
!   Export processor boundaries to file.
!
!   20-aug-09/bourdin: adapted
!
      use Cdata, only: procy_bounds,procz_bounds
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
      integer :: ierr
!
      if (lroot) then
        open(lun_output,FILE=file,FORM='unformatted',IOSTAT=ierr)
        if (ierr /= 0) call stop_it( &
            "Cannot open " // trim(file) // " (or similar) for writing" // &
            " -- is data/ visible from all nodes?")
        write(lun_output) procy_bounds
        write(lun_output) procz_bounds
        close(lun_output)
      endif
!
    endsubroutine wproc_bounds
!***********************************************************************
    subroutine rproc_bounds(file)
!
!   Import processor boundaries from file.
!
!   20-aug-09/bourdin: adapted
!
      use Cdata, only: procy_bounds,procz_bounds
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
      integer :: ierr
!
      open(lun_input,FILE=file,FORM='unformatted',IOSTAT=ierr)
      if (ierr /= 0) call stop_it( &
          "Cannot open " // trim(file) // " (or similar) for reading" // &
          " -- is data/ visible from all nodes?")
      read(lun_input) procy_bounds
      read(lun_input) procz_bounds
      close(lun_input)
!
    endsubroutine rproc_bounds
!***********************************************************************
    subroutine wtime(file,tau)
!
!  Write t to file
!  21-sep-02/wolf: coded
!
      use Mpicomm, only: stop_it
!
      double precision :: tau
      character (len=*) :: file
!
      integer :: ierr
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = tau
      if (lroot) then
        open(lun_output,FILE=file,IOSTAT=ierr)
        if (ierr /= 0) call stop_it( &
            "Cannot open " // trim(file) // " for writing")
        write(lun_output,*) t_sp
        close(lun_output)
      endif
!
    endsubroutine wtime
!***********************************************************************
    subroutine rtime(file,tau)
!
!  Read t from file
!  21-sep-02/wolf: coded
!
      use Mpicomm, only: stop_it
!
      double precision :: tau
      character (len=*) :: file
!
      integer :: ierr
      real:: t_sp   ! t in single precision for backwards compatibility
!
      open(lun_input,FILE=file,IOSTAT=ierr)
      if (ierr /= 0) call stop_it( &
          "Cannot open " // trim(file) // " for reading")
      read(lun_input,*) t_sp
      close(lun_input)
      tau = t_sp
!
    endsubroutine rtime
!***********************************************************************
endmodule Io
