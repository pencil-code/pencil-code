! $Id: io_mpidist.f90,v 1.11 2004-11-22 21:13:31 dobler Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   io_mpidist.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  IO via MPI2 (i.e. all process write to separate files, e.g.
!!!  data/proc0/var.dat
!!!  WITH NON-BLOCKING WRITES
!!!  12-jun-03/tony: started  (based on Wolfgang's io-mpio.f90 code)

module Io

  use Cdata

  implicit none
 
  interface output              ! Overload the `output' function
    module procedure output_vect
    module procedure output_scal
  endinterface

  interface output_pencil        ! Overload the `output_pencil' function
    module procedure output_pencil_vect
    module procedure output_pencil_scal
  endinterface

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
!
  include 'mpif.h'
!  include 'mpiof.h'
!
  integer :: lun_output=91

  integer, dimension(MPI_STATUS_SIZE) :: status
! LAM-MPI does not know the MPI2 constant MPI_OFFSET_KIND, 8 is hopfully OK
!  integer(kind=MPI_OFFSET_KIND) :: dist_zero=0
  integer(kind=8) :: dist_zero=0
  integer :: io_filetype_s,io_memtype_s,io_filetype_v,io_memtype_v
  integer :: fhandle,ierr
  logical :: io_initialized=.false.

  logical :: lbusy_var = .false.

contains

!***********************************************************************
    subroutine register_io()
!
!  define MPI data types needed for MPI-IO
!  closely follwing Gropp et al. `Using MPI-2'
!
!  20-sep-02/wolf: coded
!  25-oct-02/axel: removed assignment of datadir; now set in cdata.f90
!
      use General
      use Cdata
      use Sub
      use Mpicomm, only: lroot,stop_it
!
      logical, save :: first=.true.
!
!      if (.not. first) call stop_it('register_io called twice')
      first = .false.
!
!  identify version number
!
      if (lroot) call cvs_id("$Id: io_mpidist.f90,v 1.11 2004-11-22 21:13:31 dobler Exp $")
!
      io_initialized=.true.
!
!  initialize datadir and directory_snap (where var.dat and VAR# go)
!  -- may be overwritten in *.in parameter file
!
      directory_snap = ''
!
    endsubroutine register_io
!
!ajwm need and initialize_io ???
!
!
!***********************************************************************
    subroutine directory_names()
!
!  Set up the directory names:
!  set directory name for the output (one subdirectory for each processor)
!  if datadir_snap (where var.dat, VAR# go) is empty, initialize to datadir
!
!  02-oct-2002/wolf: coded
!
      use Cdata, only: datadir,directory,datadir_snap,directory_snap
      use Mpicomm, only: iproc,stop_it
      use General, only: chn, safe_character_assign
!
      character (len=5) :: chproc=''
!
!  check whether directory_snap contains `/proc0' -- if so, revert to the
!  default name.
!  Rationale: if directory_snap was not explicitly set in start.in, it
!  will be written to param.nml as 'data/proc0', but this should in fact
!  be data/procN on processor N.
!
      if ((datadir_snap == '') .or. (index(datadir_snap,'proc0')>0)) then
        datadir_snap = datadir
      endif
!
      call chn(iproc,chproc,'directory_names')
      call safe_character_assign(directory, trim(datadir)//'/proc'//chproc)
      call safe_character_assign(directory_snap, &
                                            trim(datadir_snap)//'/proc'//chproc)
!
    endsubroutine directory_names
!***********************************************************************
    subroutine commit_io_type_vect(nv,a,mode)
!
!  For a new value of nv, commit MPI types needed for output_vect(). If
!  called with the same value of nv as the previous time, do nothing.
!  12-jun-03/tony: coded
!
!      use Cdata
!
      integer, dimension(9) :: blen, datatype
      integer (kind=MPI_ADDRESS_KIND), dimension(9) :: disp
      integer,save :: lastnv=-1, lastmode ! value of nv at previous call
      real, dimension (mx,my,mz,nv) :: a
      integer :: nv, mode, mode_vars
      integer (kind=MPI_ADDRESS_KIND) :: a_addr, var_addr
      integer (kind=MPI_ADDRESS_KIND), save :: lasta_addr=-1

!               a    t   x   y   z dx dy dz deltay

!
      call MPI_GET_ADDRESS(a,a_addr,ierr)

      if (nv /= lastnv .or. a_addr/=lasta_addr .or. mode.ne.lastmode) then
        if (lastnv > 0) then
          ! free old types, so we can re-use them
          call MPI_TYPE_FREE(io_memtype_v, ierr)
          call MPI_TYPE_FREE(io_filetype_v, ierr)
        endif
  
        mode_vars=8
        if (lshear) mode_vars=9
        if (mode.eq.0) mode_vars=1

        blen     = (/ mw*nv, 1, mx, my, mz, 1, 1, 1,     1/)
        datatype(:) = MPI_REAL
        disp(1)=0
        call MPI_GET_ADDRESS(t,var_addr,ierr)
        disp(2)=var_addr - a_addr
        call MPI_GET_ADDRESS(x,var_addr,ierr)
        disp(3)=var_addr - a_addr
        call MPI_GET_ADDRESS(y,var_addr,ierr)
        disp(4)=var_addr - a_addr
        call MPI_GET_ADDRESS(z,var_addr,ierr)
        disp(5)=var_addr - a_addr
        call MPI_GET_ADDRESS(dx,var_addr,ierr)
        disp(6)=var_addr - a_addr
        call MPI_GET_ADDRESS(dy,var_addr,ierr)
        disp(7)=var_addr - a_addr
        call MPI_GET_ADDRESS(dz,var_addr,ierr)
        disp(8)=var_addr - a_addr
        call MPI_GET_ADDRESS(deltay,var_addr,ierr)
        disp(9)=var_addr - a_addr
        
!
        call MPI_TYPE_CREATE_STRUCT(mode_vars, blen(1:mode_vars), disp(1:mode_vars), &
                                    datatype(1:mode_vars), io_memtype_v, ierr)
        call MPI_TYPE_COMMIT(io_memtype_v,ierr)

        call MPI_TYPE_CONTIGUOUS(sum(blen(1:mode_vars)),MPI_REAL,io_filetype_v,ierr)
        call MPI_TYPE_COMMIT(io_filetype_v,ierr)

        lasta_addr=a_addr
      endif
!
    endsubroutine commit_io_type_vect
!***********************************************************************
    subroutine commit_io_type_scal(nv,a,mode)
!
!  For a new value of nv, commit MPI types needed for output_vect(). If
!  called with the same value of nv as the previous time, do nothing.
!  12-jun-03/tony: coded
!
!      use Cdata
!
      integer, dimension(9) :: blen, datatype
      integer (kind=MPI_ADDRESS_KIND), dimension(9) :: disp
      integer,save :: lastnv=-1, lastmode ! value of nv at previous call
      real, dimension (mx,my,mz,nv) :: a
      integer :: nv, mode, mode_vars
      integer (kind=MPI_ADDRESS_KIND) :: a_addr, var_addr
      integer (kind=MPI_ADDRESS_KIND), save :: lasta_addr=-1

!               a    t   x   y   z dx dy dz deltay

!
      call MPI_GET_ADDRESS(a,a_addr,ierr)

      if (nv /= lastnv .or. a_addr/=lasta_addr .or. mode.ne.lastmode) then
        if (lastnv > 0) then
          ! free old types, so we can re-use them
          call MPI_TYPE_FREE(io_memtype_s, ierr)
          call MPI_TYPE_FREE(io_filetype_s, ierr)
        endif
  
        mode_vars=8
        if (lshear) mode_vars=9
        if (mode.eq.0) mode_vars=1

        blen     = (/ mw*nv, 1, mx, my, mz, 1, 1, 1,     1/)
        datatype(:) = MPI_REAL
        disp(1)=0
        call MPI_GET_ADDRESS(t,var_addr,ierr)
        disp(2)=var_addr - a_addr
        call MPI_GET_ADDRESS(x,var_addr,ierr)
        disp(3)=var_addr - a_addr
        call MPI_GET_ADDRESS(y,var_addr,ierr)
        disp(4)=var_addr - a_addr
        call MPI_GET_ADDRESS(z,var_addr,ierr)
        disp(5)=var_addr - a_addr
        call MPI_GET_ADDRESS(dx,var_addr,ierr)
        disp(6)=var_addr - a_addr
        call MPI_GET_ADDRESS(dy,var_addr,ierr)
        disp(7)=var_addr - a_addr
        call MPI_GET_ADDRESS(dz,var_addr,ierr)
        disp(8)=var_addr - a_addr
        call MPI_GET_ADDRESS(deltay,var_addr,ierr)
        disp(9)=var_addr - a_addr
        
!
        call MPI_TYPE_CREATE_STRUCT(mode_vars, blen(1:mode_vars), disp(1:mode_vars), &
                                    datatype(1:mode_vars), io_memtype_s, ierr)
        call MPI_TYPE_COMMIT(io_memtype_s,ierr)

        call MPI_TYPE_CONTIGUOUS(sum(blen(1:mode_vars)),MPI_REAL,io_filetype_s,ierr)
        call MPI_TYPE_COMMIT(io_filetype_s,ierr)

        lasta_addr=a_addr
      endif
!
    endsubroutine commit_io_type_scal
!***********************************************************************
    subroutine input(file,a,nv,mode)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!  11-apr-97/axel: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
!
      character (len=*) :: file
      integer :: nv,mode                  ,i
      real, dimension (mx,my,mz,nv) :: a
!
      if (ip<=8) print*,'input: mx,my,mz,nv=',mx,my,mz,nv
      if (.not. io_initialized) &
           call stop_it("input: Need to call init_io first")
!
      call commit_io_type_vect(nv,a,mode)
!
!  open file and set view (specify which file positions we can access)
!
      call MPI_FILE_OPEN(MPI_COMM_SELF, file, &
               MPI_MODE_RDONLY, &
               MPI_INFO_NULL, fhandle, ierr)
      call MPI_FILE_SET_VIEW(fhandle, dist_zero, io_filetype_v, io_filetype_v, &
               "native", MPI_INFO_NULL, ierr)
!
!  read data
!
      call MPI_FILE_READ(fhandle, a, 1, io_memtype_v, status, ierr)
      call MPI_FILE_CLOSE(fhandle, ierr)
      
!
    endsubroutine input
!***********************************************************************
    subroutine output_vect(file,a,nv,noclose)
!
!  write snapshot file; currently without ghost zones and any meta data
!  like time, etc.
!    Looks like we nee to commit the MPI type anew each time we are called,
!  since nv may vary.
!  20-sep-02/wolf: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
!
      logical, intent(in), optional :: noclose
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'output_vect: nv =', nv
      if (.not. io_initialized) &
           call stop_it("output_vect: Need to call init_io first")
!
      call commit_io_type_vect(nv,a,1)
      !
      !  open file and set view (specify which file positions we can access)
      !
      call MPI_FILE_OPEN(MPI_COMM_SELF, file, &
               ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
               MPI_INFO_NULL, fhandle, ierr)
      call MPI_FILE_SET_VIEW(fhandle, dist_zero, MPI_REAL, io_filetype_v, &
               "native", MPI_INFO_NULL, ierr)
      !
      !  write data
      !
      call MPI_FILE_WRITE(fhandle, a, 1, io_memtype_v, status, ierr)
      call MPI_FILE_CLOSE(fhandle, ierr)
!
    endsubroutine output_vect
!***********************************************************************
    subroutine output_scal(file,a,nv,noclose)
!
!  write snapshot file; currently without ghost zones and any meta data
!  like time, etc.
!  20-sep-02/wolf: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
!
      logical, intent(in), optional :: noclose
      real, dimension (mx,my,mz) :: a
      integer :: nv
      character (len=*) :: file

      if ((ip<=8) .and. lroot) print*,'output_scal'
      if (.not. io_initialized) &
           call stop_it("output_scal: Need to call init_io first")
      if (nv /= 1) call stop_it("output_scal: called with scalar field, but nv/=1")
!
!  open file and set view (specify which file positions we can access)
!
      call commit_io_type_scal(nv,a,1)

      call MPI_FILE_OPEN(MPI_COMM_SELF, file, &
               ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
               MPI_INFO_NULL, fhandle, ierr)
      call MPI_FILE_SET_VIEW(fhandle, dist_zero, io_filetype_s, io_filetype_s, &
               "native", MPI_INFO_NULL, ierr)
!
!  write data
!
      call MPI_FILE_WRITE(fhandle, a, 1, io_memtype_s, status, ierr)
      call MPI_FILE_CLOSE(fhandle, ierr)
!
    endsubroutine output_scal
!***********************************************************************
    subroutine output_auxiliary(lun_output,nn1,nn2,a)
!
!  write auxiliary information into snapshot file
!  26-may-03/axel: adapted from output_vect
!
      use Cdata
!
      integer :: lun_output,nn1,nn2
      real, dimension (mx,my,mz,nn1+nn2) :: a
      logical :: lauxiliary
!
      if (lroot) print*, &
       'output_auxiliary: ERROR - OUTPUT AUXILIARY NOT IMPLEMENTED FOR IO_MPIDIST.F90'
!
!  determine whether we want to write auxiliary output
!  (currently we always do this provided maux>0)
!
!      lauxiliary=(maux>0)
!      if(lauxiliary) write(lun_output) a(:,:,:,nn1+1:nn1+1+nn2)
!
    endsubroutine output_auxiliary
!***********************************************************************
    subroutine output_pencil_vect(file,a,ndim)
!
!  Write snapshot file of penciled vector data (for debugging).
!  Wrapper to the C routine output_penciled_vect_c.
!
!  15-feb-02/wolf: coded
!
      use Cdata
!
      integer :: ndim
      real, dimension (nx,ndim) :: a
      character (len=*) :: file
!
      if (ip<9.and.lroot.and.imn==1) &
           print*,'output_pencil_vect('//file//'): ndim=',ndim
!
      if (headt .and. (imn==1)) write(*,'(A)') &
           'output_pencil_vect: Writing to ' // trim(file) // &
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
      use Mpicomm, only: lroot,stop_it

!
      integer :: ndim
      real, dimension (nx) :: a
      character (len=*) :: file
!
      if ((ip<=8) .and. lroot .and. imn==1) &
           print*,'output_pencil_scal('//file//')'
!
      if (ndim /= 1) &
           call stop_it("output_pencil_scal: called with scalar field, but ndim/=1")
!
      if (headt .and. (imn==1)) print*, &
           'output_pencil_scal: Writing to ', trim(file), &
           ' for debugging -- this may slow things down'
!
      call output_penciled_scal_c(file, a, ndim, &
                                  imn, mm(imn), nn(imn), t, &
                                  nx, ny, nz, nghost, len(file))
!
    endsubroutine output_pencil_scal
!***********************************************************************
    subroutine outpus(file,a,nv)
!
!  write snapshot file, always write mesh and time, could add other things
!  11-oct-98/axel: adapted
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
!
      integer :: nv
      character (len=*) :: file
      real, dimension (mx,my,mz,nv) :: a
!
      call stop_it("output: doesn't work with io_mpio yet -- but wasn't used anyway")
!
      open(1,FILE=file,FORM='unformatted')
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
      use Mpicomm, only: stop_it_if_any
!
      character (len=*) :: file
      logical :: ioerr
!
      ioerr = .true.            ! will be overridden unless we go 911
      open(1,FILE=file,FORM='unformatted')
      write(1) t,x,y,z,dx,dy,dz
      write(1) dx,dy,dz
      close(1)
      ioerr = .false.
!
!  Something went wrong. Catches cases that would make mpich 1.x hang,
!  provided that this is the first collective write call
!
911   call stop_it_if_any(ioerr, &
          "Cannot open " // trim(file) // " (or similar) for writing" // &
          " -- is data/ visible from all nodes?")
!
    endsubroutine wgrid
!***********************************************************************
    subroutine rgrid (file)
!
!  Read processor-local part of grid coordinates.
!  21-jan-02/wolf: coded
!
      use Cdata
!
      real :: tdummy
      character (len=*) :: file
!
      open(1,FILE=file,FORM='unformatted')
      read(1) tdummy,x,y,z,dx,dy,dz
      read(1) dx,dy,dz
      close(1)
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

      Lx=dx*nx*nprocx
      Ly=dy*ny*nprocy
      Lz=dz*nz*nprocz
!
      if (ip<=4) print*,'rgrid: dx,dy,dz=',dx,dy,dz
      if (ip<=4) print*,'rgrid: dxmin,dxmax=',dxmin,dxmax
!
    endsubroutine rgrid
!***********************************************************************
    subroutine wtime(file,tau)
!
      real :: tau
      character (len=*) :: file
!
      if (.false.) print*,tau,file
    endsubroutine wtime
!***********************************************************************
    subroutine rtime(file,tau)
!
      real :: tau
      character (len=*) :: file
!
      if (.false.) print*,tau,file
    endsubroutine rtime
!***********************************************************************

endmodule Io
