! $Id$
!
!  Distributed IO (i.e. each process writes its own file data/procX)
!
!  The file format written by output_snap() (and used, e.g. in var.dat)
!  consists of the followinig Fortran records:
!    1. data(mx,my,mz,nvar)
!    2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
!  Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
!  for one vector field, 8 for var.dat in the case of MHD with entropy.
!
!  04-nov-11/MR: IOSTAT handling generally introduced
!  16-nov-11/MR: calls to outlog adapted
!  10-Dec-2011/Bourdin.KIS: major cleanup
!  01-Jun-2015/Bourdin.KIS: outlog removed to see clear-text error messages
!
module Io
!
  use Cdata
  use Cparam, only: intlen, fnlen, max_int
  use File_io, only: delete_file, file_exists
  use Messages, only: fatal_error, warning, svn_id
!
  implicit none
!
  include 'io.h'
  include 'record_types.h'
!
  interface read_snap
    module procedure read_snap_double
    module procedure read_snap_single
  endinterface
!
  interface read_part
    module procedure read_part_double
    module procedure read_part_single
  endinterface
!
  interface read_globals
    module procedure read_globals_double
    module procedure read_globals_single
  endinterface
!
  interface input_grid
    module procedure input_grid_double
    module procedure input_grid_single
  endinterface
!
  interface input_proc_bounds
    module procedure input_proc_bounds_double
    module procedure input_proc_bounds_single
  endinterface
!
  ! set ireset_tstart to 1 or 2 to coordinate divergent timestamp
  integer, parameter :: MINT=1, MAXT=2
!
! Indicates if IO is done distributed (each proc writes into a procdir)
! or collectively (eg. by specialized IO-nodes or by MPI-IO).
!
  logical :: lcollective_IO=.false.
  character (len=labellen) :: IO_strategy="dist"
!
  integer :: persist_last_id=-max_int
!
  contains
!***********************************************************************
    subroutine register_io
!
!  dummy routine, generates separate directory for each processor.
!  VAR#-files are written to the directory directory_snap which will
!  be the same as directory, unless specified otherwise.
!
!  20-sep-02/wolf: coded
!
      use Mpicomm, only: lroot
!
!  identify version number
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_io
!***********************************************************************
    subroutine finalize_io
!
    endsubroutine finalize_io
!***********************************************************************
    subroutine directory_names
!
!  Set up the directory names:
!  set directory name for the output (one subdirectory for each processor)
!  if datadir_snap (where var.dat, VAR# go) is empty, initialize to datadir
!
!  02-oct-2002/wolf: coded
!
      use General, only: directory_names_std
!
!  check whether directory_snap contains `/proc0' -- if so, revert to the
!  default name.
!  Rationale: if directory_snap was not explicitly set in start.in, it
!  will be written to param.nml as 'data/proc0', but this should in fact
!  be data/procN on processor N.
!
      if ((datadir_snap == '') .or. (index(datadir_snap,'proc0')>0)) &
        datadir_snap = datadir
!
      call directory_names_std(.true.)      
!
    endsubroutine directory_names
!***********************************************************************
    subroutine output_snap(a,nv1,nv2,file)
!
!  Write snapshot file, always write time and mesh, could add other things.
!
!  11-apr-97/axel: coded
!  13-Dec-2011/Bourdin.KIS: reworked
!  13-feb-2014/MR: made 'file' optional, 'a' assumed-shape (for downsampled output);
!                  moved donwsampling stuff to snapshot
!
      use Mpicomm, only: start_serialize, end_serialize
      use General, only: get_range_no, ioptest
!
      real, dimension (:,:,:,:),  intent(IN) :: a
      integer,           optional,intent(IN) :: nv1,nv2
      character (len=*), optional,intent(IN) :: file
!
      real :: t_sp   ! t in single precision for backwards compatibility
      integer :: na, ne
!
      t_sp = real (t)
!
      if (lserial_io) call start_serialize
      if (present(file)) then
        call delete_file(trim(directory_snap)//'/'//file)
        open (lun_output, FILE=trim(directory_snap)//'/'//file, FORM='unformatted', status='new')
      endif

      na=ioptest(nv1,1)
      ne=ioptest(nv2,mvar_io)
      if (lroot .and. (ip <= 8)) print *, 'output_snap: nv1,nv2 =', na,ne
!
      if (lwrite_2d) then
        if (nx == 1) then
          write (lun_output) a(l1,:,:,na:ne)
        elseif (ny == 1) then
          write (lun_output) a(:,m1,:,na:ne)
        elseif (nz == 1) then
          write (lun_output) a(:,:,n1,na:ne)
        else
          call fatal_error ('output_snap', 'lwrite_2d used for 3D simulation!')
        endif
      else
        write (lun_output) a(:,:,:,na:ne)
      endif
!
!  Write shear at the end of x,y,z,dx,dy,dz.
!  At some good moment we may want to treat deltay like with
!  other modules and call a corresponding i/o parameter module.
!
      if (lshear) then
        write (lun_output) t_sp, x(1:size(a,1)), y(1:size(a,2)), z(1:size(a,3)), dx, dy, dz, deltay
      else
        write (lun_output) t_sp, x(1:size(a,1)), y(1:size(a,2)), z(1:size(a,3)), dx, dy, dz
      endif
!
      if (lserial_io) call end_serialize
!
    endsubroutine output_snap
!***********************************************************************
    subroutine output_snap_finalize
!
!  Close snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: adapted from output_snap
!
      use Mpicomm, only: end_serialize
!
      if (persist_initialized) then
        write (lun_output) id_block_PERSISTENT
        persist_initialized = .false.
        persist_last_id = -max_int
      endif
!
      close (lun_output)
!
      if (lserial_io) call end_serialize
!
    endsubroutine output_snap_finalize
!***********************************************************************
    subroutine output_average_2D(label, nc, name, data, time, lwrite, header)
!
!   Output average to a file.
!
!   16-Nov-2018/PABourdin: coded
!   08-nov-2022/ccyang: remove plain text output (never used)
!
      use Cdata, only: directory_dist
      use General, only: keep_compiler_quiet
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: nc
      character (len=fmtlen), dimension(nc), intent(in) :: name
      real, dimension(:,:,:), intent(in) :: data
      real, intent(in) :: time
      logical, intent(in) :: lwrite
      real, dimension(:), optional, intent(in) :: header
!
      character (len=fnlen) :: filename
      integer :: ia
!
      call keep_compiler_quiet(name)
!
      if (.not. lwrite .or. (nc <= 0)) return
!
      filename = trim(directory_dist) // '/' // trim(label) // 'averages.dat'
      open(lun_output, file=filename, form='unformatted', position='append')
      if (present (header)) write(lun_output) header
      write(lun_output) time
      if (label(1:1) == 'z') then
        write(lun_output) ( data(ia,:,:), ia=1, nc )
      else
        write(lun_output) data(:,:,1:nc)
      endif
      close(lun_output)
!
    endsubroutine output_average_2D
!***********************************************************************
    subroutine output_slice_position()
!
!  Record slice positions.
!
!  13-nov-20/ccyang: wrapper
!
      use HDF5_IO, only: hdf5_output_slice_position
!
      call hdf5_output_slice_position()
!
    endsubroutine output_slice_position
!***********************************************************************
    subroutine output_slice(lwrite, time, label, suffix, pos, grid_pos, data)
!
!  Append to a slice file
!
!  13-nov-20/ccyang: wrapper
!
      use HDF5_IO, only: hdf5_output_slice
!
      logical, intent(in) :: lwrite
      real, intent(in) :: time
      character(len=*), intent(in) :: label, suffix
      real, intent(in) :: pos
      integer, intent(in) :: grid_pos
      real, dimension(:,:), pointer :: data
!
      call hdf5_output_slice(lwrite, time, label, suffix, pos, grid_pos, data)
!
    endsubroutine output_slice
!***********************************************************************
    subroutine output_part_snap(ipar, ap, mv, nv, file, label, ltruncate)
!
!  Write particle snapshot file, always write mesh and time.
!
!  22-Oct-2018/PABourdin: adapted from output_snap
!
      integer, intent(in) :: mv, nv
      integer, dimension (mv), intent(in) :: ipar
      real, dimension (mv,mparray), intent(in) :: ap
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: ltruncate
!
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = real (t)
      open(lun_output,FILE=trim(directory_dist)//'/'//file,FORM='unformatted')
!
      ! write the number of particles present at the processor
      write(lun_output) nv
      if (nv /= 0) then
        ! write particle index
        write(lun_output) ipar(1:nv)
        ! write particle data
        write(lun_output) ap(1:nv,:)
      endif
!
      ! write time and grid parameters
      write(lun_output) t_sp, x, y, z, dx, dy, dz
!
      close(lun_output)
!
    endsubroutine output_part_snap
!***********************************************************************
    subroutine output_stalker_init(num, nv, snap, ID)
!
!  Open stalker particle snapshot file and initialize with snapshot time.
!
!  03-May-2019/PABourdin: coded
!
      integer, intent(in) :: num, nv, snap
      integer, dimension(nv), intent(in) :: ID
!
      real :: t_sp
!
      call output_stalker ('', 0, nv, (/ 0.0 /), num)
!
      open (lun_output, file=trim(directory_dist)//'/particles_stalker.dat', form='unformatted', position='append')
!
      ! write the time, number, and ID of stalked particles at this processor
      t_sp = t
      write (lun_output) t_sp, nv
      if (nv >= 1) write (lun_output) ID
!
    endsubroutine output_stalker_init
!***********************************************************************
    subroutine output_stalker(label, mv, nv, data, nvar, lfinalize)
!
!  Write stalker particle quantity to snapshot file.
!
!  03-May-2019/PABourdin: coded
!
      use General, only: loptest
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: mv, nv
      real, dimension (mv), intent(in) :: data
      integer, intent(in), optional :: nvar
      logical, intent(in), optional :: lfinalize
!
      character (len=fnlen) :: dataset
      real, dimension(:,:), allocatable, save :: stalk_data
      integer, save :: pos = 0
!
      if (loptest (lfinalize)) then
        ! deallocate temporary stalker particle space
        if (pos > 0) then
          write (lun_output) stalk_data
          pos = 0
        endif
        if (allocated (stalk_data)) deallocate (stalk_data)
        return
      endif
!
      if (present (nvar)) then
        ! allocate temporary stalker particle space
        if (nv > 0) allocate (stalk_data (nvar, nv))
        return
      endif
!
      if (nv > 0) then
        ! write stalker particles to temporary space
        pos = pos + 1
        stalk_data(pos,:) = data(1:nv)
      endif
!
    endsubroutine output_stalker
!***********************************************************************
    subroutine output_part_finalize
!
!  Close particle snapshot file.
!
!  03-May-2019/PABourdin: coded
!
      call output_stalker ('', 0, 1, (/ 0.0 /), 0, lfinalize=.true.)
      close (lun_output)
!
    endsubroutine output_part_finalize
!***********************************************************************
    subroutine output_pointmass(file, labels, fq, mv, nc)
!
!  Write pointmass snapshot file with time.
!
!  26-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mv, nc
      character (len=*), dimension (nc), intent(in) :: labels
      real, dimension (mv,nc), intent(in) :: fq
!
      if (.not. lroot) return
!
      open(lun_output,FILE=trim(directory_snap)//'/'//trim(file),FORM='unformatted')
      write(lun_output) mv
      if (mv > 0) write(lun_output) fq
      write(lun_output) t
      close(lun_output)
      if (ip<=10) print*,'written pointmass snapshot ', trim (file)
!
    endsubroutine output_pointmass
!***********************************************************************
    logical function init_write_persist(file)
!
!  Initialize writing of persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in), optional :: file
!
      character (len=fnlen), save :: filename=""
!
      persist_last_id = -max_int
      init_write_persist = .false.
!
      if (present (file)) then
        filename = file
        persist_initialized = .false.
        return
      endif
!
      if (filename /= "") then
        close (lun_output)
        call delete_file(trim (directory_snap)//'/'//filename)
        open (lun_output, FILE=trim (directory_snap)//'/'//filename, FORM='unformatted', status='new')
        filename = ""
      endif
!
      if (lroot .and. (ip <= 9)) write (*,*) 'begin persistent block'
      write (lun_output) id_block_PERSISTENT
      init_write_persist = .false.
!
    endfunction init_write_persist
!***********************************************************************
    logical function write_persist_id(label, id)
!
!  Write persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
!
!write(20+iproc,*) 'write_persist_id, label, persist_initialized=',label,persist_initialized 
      write_persist_id = .true.
      if (.not. persist_initialized) write_persist_id = init_write_persist()
      if (.not. persist_initialized) return
!
      if (persist_last_id /= id) then
        write (lun_output) id
        persist_last_id = id
      endif
!
      write_persist_id = .false.
!
    endfunction write_persist_id
!***********************************************************************
    logical function write_persist_logical_0D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, intent(in) :: value
!
      write_persist_logical_0D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output) value
      write_persist_logical_0D = .false.
!
    endfunction write_persist_logical_0D
!***********************************************************************
    logical function write_persist_logical_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, dimension(:), intent(in) :: value
!
      write_persist_logical_1D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output) value
      write_persist_logical_1D = .false.
!
    endfunction write_persist_logical_1D
!***********************************************************************
    logical function write_persist_int_0D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, intent(in) :: value
!
      write_persist_int_0D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output) value
      write_persist_int_0D = .false.
!
    endfunction write_persist_int_0D
!***********************************************************************
    logical function write_persist_int_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, dimension(:), intent(in) :: value
!
      write_persist_int_1D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output) value
      write_persist_int_1D = .false.
!
    endfunction write_persist_int_1D
!***********************************************************************
    logical function write_persist_real_0D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, intent(in) :: value
!
      write_persist_real_0D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output) value
      write_persist_real_0D = .false.
!
    endfunction write_persist_real_0D
!***********************************************************************
    logical function write_persist_real_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, dimension(:), intent(in) :: value
!
      write_persist_real_1D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output) value
      write_persist_real_1D = .false.
!
    endfunction write_persist_real_1D
!***********************************************************************
    logical function write_persist_torus_rect(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  16-May-2020/MR: coded
!
      use Geometrical_types

      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      type(torus_rect), intent(in) :: value
!
      write_persist_torus_rect = .true.
      if (write_persist_id (label, id)) return
!
      !write (lun_output) value
      write_persist_torus_rect = .false.
!
    endfunction write_persist_torus_rect
!***********************************************************************
    subroutine input_snap(file,a,nv,mode)
!
!  manages reading of snapshot from different precision
!
!  24-oct-13/MR: coded
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real, dimension (mx,my,mz,nv), intent(out) :: a

      real(KIND=rkind8), dimension(:,:,:,:), allocatable :: adb
      real(KIND=rkind4), dimension(:,:,:,:), allocatable :: asg

      real(KIND=rkind8), dimension(:), allocatable :: xdb,ydb,zdb
      real(KIND=rkind4), dimension(:), allocatable :: xsg,ysg,zsg

      real(KIND=rkind8) :: dxdb,dydb,dzdb,deltaydb
      real(KIND=rkind4) :: dxsg,dysg,dzsg,deltaysg

      if (lread_from_other_prec) then
        if (kind(a)==rkind4) then
          allocate(adb(mx,my,mz,nv),xdb(mx),ydb(my),zdb(mz))
          call read_snap(file,adb,xdb,ydb,zdb,dxdb,dydb,dzdb,deltaydb,nv,mode)
          a=adb; x=xdb; y=ydb; z=zdb; dx=dxdb; dy=dydb; dz=dzdb; deltay=deltaydb
        elseif (kind(a)==rkind8) then
          allocate(asg(mx,my,mz,nv),xsg(mx),ysg(my),zsg(mz))
          call read_snap(file,asg,xsg,ysg,zsg,dxsg,dysg,dzsg,deltaysg,nv,mode)
          a=asg; x=xsg; y=ysg; z=zsg; dx=dxsg; dy=dysg; dz=dzsg; deltay=deltaysg
        endif
      else
        call read_snap(file,a,x,y,z,dx,dy,dz,deltay,nv,mode)
      endif

    endsubroutine input_snap
!***********************************************************************
    subroutine read_snap_single(file,a,x,y,z,dx,dy,dz,deltay,nv,mode)
!
!  Read snapshot file in single precision, possibly with mesh and time (if mode=1).
!
!  24-oct-13/MR: derived from input_snap
!  28-oct-13/MR: consistency check for t_sp relaxed for restart from different precision
!   6-mar-14/MR: if timestamp of snapshot inconsistent, now three choices:
!                if ireset_tstart= 0: cancel program
!                                = MINT, tstart unspecified: use minimum time of all var.dat 
!                                                        for start
!                                = MAXT, tstart unspecified: use maximum time of all var.dat 
!                                                        for start
!                                > 0, tstart specified: use this value
!  13-feb-20/MR: added possibility to omit a range of variables when reading the snapshot,
!                range is defined as [ivar_omi1,ivar_omi2]
!  22-oct-20/MR: outsourced code shared with read_snap_double into io_dist.h
!                Added "rescue mechanism": if snapshot is not properly readable,
!                reding is tried from the (min. 11, max. 26) neighboring
!                processors' snaphots.
!
      use Mpicomm, only: start_serialize, end_serialize, mpibcast, mpibcast_real, mpiallreduce_or, &
                         stop_it, mpiallreduce_min, mpiallreduce_max, MPI_COMM_PENCIL, &
                         mpisend_char, mpirecv_char, mpiallreduce_and
      use Syscalls, only: islink, system_cmd, readlink, extract_str
      use General, only: itoa, find_proc, atoi
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real(KIND=rkind4), dimension (mx,my,mz,nv), intent(out) :: a
!
      real(KIND=rkind4) :: t_sp, t_red

      real(KIND=rkind4),                 intent(out) :: dx, dy, dz, deltay
      real(KIND=rkind4), dimension (mx), intent(out) :: x
      real(KIND=rkind4), dimension (my), intent(out) :: y
      real(KIND=rkind4), dimension (mz), intent(out) :: z

      real(KIND=rkind4), dimension(:,:,:,:), allocatable :: tmp_omit,tmp
!
      include 'io_dist.h'
!
    endsubroutine read_snap_single
!***********************************************************************
    subroutine read_snap_double(file,a,x,y,z,dx,dy,dz,deltay,nv,mode)
!
!  Read snapshot file in double precision, possibly with mesh and time (if mode=1).
!
!  24-oct-13/MR: derived from input_snap
!  28-oct-13/MR: consistency check for t_sp relaxed for restart from different precision
!   6-mar-14/MR: if timestamp of snapshot inconsistent, now three choices:
!                if ireset_tstart= 0: cancel program
!                                = MINT, tstart unspecified: use minimum time of all var.dat 
!                                                        for start
!                                = MAXT, tstart unspecified: use maximum time of all var.dat 
!                                                        for start
!                                > 0, tstart specified: use this value
!                             
      use Mpicomm, only: start_serialize, end_serialize, mpibcast, mpibcast_real, mpiallreduce_or, &
                         stop_it, mpiallreduce_min, mpiallreduce_max, MPI_COMM_PENCIL, &
                         mpisend_char, mpirecv_char,mpiallreduce_and
      use Syscalls, only: islink, system_cmd, readlink, extract_str
      use General, only: itoa, find_proc, atoi
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real(KIND=rkind8), dimension (mx,my,mz,nv), intent(out) :: a
!
      real(KIND=rkind8) :: t_sp, t_red

      real(KIND=rkind8), intent(out) :: dx, dy, dz, deltay
      real(KIND=rkind8), dimension (mx), intent(out) :: x
      real(KIND=rkind8), dimension (my), intent(out) :: y
      real(KIND=rkind8), dimension (mz), intent(out) :: z

      real(KIND=rkind8), dimension(:,:,:,:), allocatable :: tmp_omit,tmp

      include 'io_dist.h'
!
    endsubroutine read_snap_double
!***********************************************************************
    subroutine input_snap_finalize
!
!  Close snapshot file.
!
!  11-apr-97/axel: coded
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use Mpicomm, only: end_serialize
      use Syscalls, only: system_cmd
!
      close (lun_input)
      if (snaplink/='') then
        call system_cmd('rm -f '//snaplink)
        snaplink=''
      endif

      if (lserial_io) call end_serialize
!
    endsubroutine input_snap_finalize
!***********************************************************************
    subroutine input_slice_real_arr(file, time, pos, data)
!
!  Read a slice file.
!
!  13-nov-20/ccyang: wrapper
!
      use HDF5_IO, only: hdf5_input_slice
!
      character(len=*), intent(in) :: file
      real, intent(out):: time, pos
      real, dimension(:,:,:), intent(out):: data
!
      call hdf5_input_slice(file, time, pos, data)
!
    endsubroutine input_slice_real_arr
!***********************************************************************
    subroutine input_slice_scat_arr(file, pos, data, ivar, nt)
!
!  Read a slice file.
!
!  13-nov-20/ccyang: wrapper
!
      use General, only: scattered_array
      use HDF5_IO, only: hdf5_input_slice
!
      character(len=*), intent(in) :: file
      real, intent(out):: pos
      type(scattered_array), pointer :: data   !intent(inout)
      integer, intent(in) :: ivar
      integer, intent(out):: nt
!
      call hdf5_input_slice(file, pos, data, ivar, nt)
!
    endsubroutine input_slice_scat_arr
!***********************************************************************
    subroutine input_part_snap(ipar, ap, mv, nv, npar_total, file, label)
!
!  Read particle snapshot file, mesh and time are read in 'input_snap'.
!
!  24-Oct-2018/PABourdin: apadpted and moved to IO module
!
      use Mpicomm, only: mpireduce_sum_int
!
      integer :: mparray_tmp=1
      integer, intent(in) :: mv
      integer, dimension (mv), intent(out) :: ipar
      real, dimension (mv,mparray), intent(out) :: ap
      integer, intent(out) :: nv, npar_total
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label
!
      open (lun_input, FILE=trim(directory_dist)//'/'//file, FORM='unformatted')
      ! Read number of particles for this processor.
      read (lun_input) nv
      if (nv > 0) then
        ! Read particle number identifier and then particle data.
        read (lun_input) ipar(1:nv)
        if (lread_oldsnap_nosink) then
          if (lroot) print*,'read old snapshot file (but without sink particles)'
          read (lun_input) ap(1:nv,1:mparray-1)
          mparray_tmp=mparray !(trick to prevent compiler problem for mparray=0)
          ap(1:nv,mparray_tmp)=0.0
        else
          read (lun_input) ap(1:nv,:)
        endif
      endif
      close (lun_input)
      if (ip<=8 .and. lroot) print*, 'input_particles: read ', file
!
      ! Sum the total number of all particles on the root processor.
      call mpireduce_sum_int (nv, npar_total)
!
    endsubroutine input_part_snap
!***********************************************************************
    subroutine input_pointmass(file, labels, fq, mv, nc)
!
!  Read pointmass snapshot file.
!
!  26-Oct-2018/PABourdin: coded
!
      use Mpicomm, only: mpibcast_real
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mv, nc
      character (len=*), dimension (nc), intent(in) :: labels
      real, dimension (mv,nc), intent(out) :: fq
!
      integer :: mv_in
!
      if (lroot) then
        if (ip<=8) print*, 'read pointmass snapshot', trim (file)
        open(lun_input,FILE=trim(directory_snap)//'/'//trim(file),FORM='unformatted')
        read(lun_input) mv_in
        if (mv_in /= mv) call fatal_error("input_pointmass","dimensions differ between file and 'fq' array")
        if (mv_in > 0) read(lun_input) fq
        close(lun_input)
      endif
!
      call mpibcast_real(fq,(/mv,nc/))
!
    endsubroutine input_pointmass
!***********************************************************************
    logical function init_read_persist(file)
!
!  Initialize writing of persistent data to persistent file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      use File_io, only: file_exists
      use Mpicomm, only: mpibcast_logical, MPI_COMM_WORLD
!
      character (len=*), intent(in), optional :: file
!
      init_read_persist = .true.
!
      if (present (file)) then
        if (lroot) init_read_persist = .not. file_exists (trim (directory_snap)//'/'//file)
        call mpibcast_logical(init_read_persist,comm=MPI_COMM_WORLD)
        if (init_read_persist) return
      endif
!
      if (present (file)) then
        close (lun_input)
        open (lun_input, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', status='old')
      endif
!
      if (lroot .and. (ip <= 9)) write (*,*) 'begin persistent block'
      init_read_persist = .false.
!
    endfunction init_read_persist
!***********************************************************************
    logical function persist_exists(label)
!
!  Dummy routine
!
!  12-Oct-2019/PABourdin: coded
!
      character (len=*), intent(in) :: label
!
      persist_exists = .false.
!
    endfunction persist_exists
!***********************************************************************
    logical function read_persist_id(label, id, lerror_prone)
!
!  Read persistent block ID from snapshot file.
!  Each persist data block starts with a unique ID followed by the data.
!  Both are Fortran records. Returns a positive ID number on success.
!
!  17-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: mpibcast_int, MPI_COMM_WORLD
!
      character (len=*), intent(in) :: label
      integer, intent(out) :: id
      logical, intent(in), optional :: lerror_prone
!
      logical :: lcatch_error
      integer :: io_err
!
      lcatch_error = .false.
      if (present (lerror_prone)) lcatch_error = lerror_prone
!
      if (lroot .and. (ip <= 9)) write (*,*) 'read persistent ID '//trim (label)
      if (lcatch_error) then
        read (lun_input, iostat=io_err) id
        if (io_err /= 0) id = -max_int
      else
        read (lun_input) id
      endif
!
      call mpibcast_int(id,comm=MPI_COMM_WORLD)
!
      read_persist_id = .false.
      if (id == -max_int) read_persist_id = .true.
!
    endfunction read_persist_id
!***********************************************************************
    logical function read_persist_logical_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      logical, intent(out) :: value
!
      read (lun_input) value
      read_persist_logical_0D = .false.
!
    endfunction read_persist_logical_0D
!***********************************************************************
    logical function read_persist_logical_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      logical, dimension(:), intent(out) :: value
!
      read (lun_input) value
      read_persist_logical_1D = .false.
!
    endfunction read_persist_logical_1D
!***********************************************************************
    logical function read_persist_int_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      integer, intent(out) :: value
!
      read (lun_input) value
      read_persist_int_0D = .false.
!
    endfunction read_persist_int_0D
!***********************************************************************
    logical function read_persist_int_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      integer, dimension(:), intent(out) :: value
!
      read (lun_input) value
      read_persist_int_1D = .false.
!
    endfunction read_persist_int_1D
!***********************************************************************
    logical function read_persist_real_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!  23-oct-2013/MR: modified for reading of different precision
!
      character (len=*), intent(in) :: label
      real, intent(out) :: value
!
      real(KIND=rkind8) :: vdb
      real(KIND=rkind4) :: vsg
!
      if (lread_from_other_prec) then
        if (kind(value)==rkind4) then
          read (lun_input) vdb
          value=vdb
        elseif (kind(value)==rkind8) then
          read (lun_input) vsg
          value=vsg
        endif
      else
        read (lun_input) value
      endif

      read_persist_real_0D = .false.
!
    endfunction read_persist_real_0D
!***********************************************************************
    logical function read_persist_real_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!  23-oct-2013/MR: modified for reading of different precision
!
      character (len=*), intent(in) :: label
      real, dimension(:), intent(out) :: value
!
      real(KIND=rkind8), dimension(:), allocatable :: vdb
      real(KIND=rkind4), dimension(:), allocatable :: vsg
!
      if (lread_from_other_prec) then
        if (kind(value)==rkind4) then
          allocate(vdb(size(value)))
          read (lun_input) vdb
          value=vdb
        elseif (kind(value)==rkind8) then
          allocate(vsg(size(value)))
          read (lun_input) vsg
          value=vsg
        endif
      else
        read (lun_input) value
      endif
!
      read_persist_real_1D = .false.
!
    endfunction read_persist_real_1D
!***********************************************************************
    logical function read_persist_torus_rect(label,value)
!
!  Read persistent data from snapshot file.
!
!  16-May-2020/MR: coded
!
      use Geometrical_types

      character (len=*), intent(in) :: label
      type(torus_rect), intent(out) :: value
!
      !read (lun_input) value
      read_persist_torus_rect = .false.
!
    endfunction read_persist_torus_rect
!***********************************************************************
    subroutine output_globals(file, a, nv, label)
!
!  Write snapshot file of globals, ignoring mesh.
!
!  10-nov-06/tony: coded
!
      use Mpicomm, only: start_serialize, end_serialize
!
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
      character (len=*), intent(in), optional :: label
!
      if (lserial_io) call start_serialize
      open(lun_output,FILE=trim(directory_snap)//'/'//file,FORM='unformatted',status='replace')
!
      if (lwrite_2d) then
        if (nx==1) then
          write(lun_output) a(nghost+1,:,:,:)
        elseif (ny==1) then
          write(lun_output) a(:,nghost+1,:,:)
        elseif (nz==1) then
          write(lun_output) a(:,:,nghost+1,:)
        else
          call fatal_error('output_globals','lwrite_2d used for 3-D simulation!')
        endif
      else
        write(lun_output) a
      endif
      close(lun_output)
!
      if (lserial_io) call end_serialize
!
    endsubroutine output_globals
!***********************************************************************
    subroutine input_globals(file, a, nv)
!
!  Read globals snapshot file, ignoring mesh.
!
!  10-nov-06/tony: coded
!
      use Mpicomm, only: start_serialize,end_serialize
!
      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
!
      real(KIND=rkind8), dimension(:,:,:,:), allocatable :: adb
      real(KIND=rkind4), dimension(:,:,:,:), allocatable :: asg

!
      if (lserial_io) call start_serialize
!
      open(lun_input,FILE=trim(directory_snap)//'/'//file,FORM='unformatted',status='old')

      if (lread_from_other_prec) then
        if (kind(a)==rkind4) then
          allocate(adb(mx,my,mz,nv))
          call read_globals(adb)
          a=adb
        elseif (kind(a)==rkind8) then
          allocate(asg(mx,my,mz,nv))
          call read_globals(asg)
          a=asg
        endif
      else
        call read_globals(a)
      endif

      close(lun_input)
!
      if (lserial_io) call end_serialize
!
    endsubroutine input_globals
!***********************************************************************
    subroutine read_globals_double(a)
!
!  Read globals snapshot file in double precision
!
!  23-oct-13/MR  : derived from input_globals
!
      real(KIND=rkind8), dimension (:,:,:,:) :: a
!
      if (lwrite_2d) then
        if (nx==1) then
          read(lun_input) a(nghost+1,:,:,:)
        elseif (ny==1) then
          read(lun_input) a(:,nghost+1,:,:)
        elseif (nz==1) then
          read(lun_input) a(:,:,nghost+1,:)
        else
          call fatal_error('input_globals','lwrite_2d used for 3-D simulation!')
        endif
      else
        read(lun_input) a
      endif
!
    endsubroutine read_globals_double
!***********************************************************************
    subroutine read_globals_single(a)
!
!  Read globals snapshot file in single precision
!
!  23-oct-13/MR  : derived from input_globals
!
      real(KIND=rkind4), dimension (:,:,:,:) :: a
!
      if (lwrite_2d) then
        if (nx==1) then
          read(lun_input) a(nghost+1,:,:,:)
        elseif (ny==1) then
          read(lun_input) a(:,nghost+1,:,:)
        elseif (nz==1) then
          read(lun_input) a(:,:,nghost+1,:)
        else
          call fatal_error('input_globals','lwrite_2d used for 3-D simulation!')
        endif
      else
        read(lun_input) a
      endif
!
    endsubroutine read_globals_single
!***********************************************************************
    subroutine log_filename_to_file(file,flist)
!
!  In the directory containing `filename', append one line to file
!  `flist' containing the file part of filename
!
      use General, only: parse_filename, safe_character_assign
!
      character (len=*) :: file,flist
!
      character (len=fnlen) :: dir,fpart
!
      call parse_filename(file,dir,fpart)
      if (dir == '.') call safe_character_assign(dir,directory_snap)
!
      open(lun_output,FILE=trim(dir)//'/'//flist,POSITION='append')
      write(lun_output,'(A,1x,e16.8)') trim(fpart), t
      close(lun_output)
!
      if (lcopysnapshots_exp) then
        if (lroot) then
          open(lun_output,FILE=trim(datadir)//'/move-me.list',POSITION='append')
          write(lun_output,'(A)') trim(fpart)
          close(lun_output)
        endif
      endif
!
    endsubroutine log_filename_to_file
!***********************************************************************
    subroutine wgrid(file,mxout,myout,mzout,lwrite)
!
!  Write processor-local part of grid coordinates.
!
!  21-jan-02/wolf: coded
!  15-jun-03/axel: Lx,Ly,Lz are now written to file (Tony noticed the mistake)
!  30-sep-13/MR  : optional parameters mxout,myout,mzout added
!                  to be able to output coordinate vectors with lengths differing from
!                  mx,my,mz
!  25-Aug-2016/PABourdin: now a global "grid.dat" is always written from all IO modules
!   5-oct-2016/MR: modifications for output of downsampled grid
!
      use Mpicomm, only: collect_grid
      use General, only: ioptest, loptest
!
      character (len=*) :: file
      integer, optional :: mxout,myout,mzout
      logical, optional :: lwrite
!
      integer :: mxout1,myout1,mzout1
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (loptest(lwrite,.not.luse_oldgrid)) then
        mxout1=ioptest(mxout,mx)
        myout1=ioptest(myout,my)
        mzout1=ioptest(mzout,mz)
!
        t_sp = real (t)

        open (lun_output,FILE=trim(directory)//'/'//file,FORM='unformatted',status='replace')
        write(lun_output) t_sp,x(:mxout1),y(:myout1),z(:mzout1),dx,dy,dz
        write(lun_output) dx,dy,dz
        write(lun_output) Lx,Ly,Lz
        write(lun_output) dx_1(:mxout1),dy_1(:myout1),dz_1(:mzout1)
        write(lun_output) dx_tilde(:mxout1),dy_tilde(:myout1),dz_tilde(:mzout1)
        close(lun_output)
!
        if (lyang) return      ! grid collection only needed on Yin grid, as grids are identical

        ! write also a global "data/allprocs/grid.dat"
        if (lroot) then
          if (ldownsampling) then
            allocate (gx(ceiling(real(nxgrid)/downsampl(1))+2*nghost), &
                      gy(ceiling(real(nygrid)/downsampl(2))+2*nghost), &
                      gz(ceiling(real(nzgrid)/downsampl(3))+2*nghost), stat=alloc_err)
          else
            allocate (gx(nxgrid+2*nghost), gy(nygrid+2*nghost), gz(nzgrid+2*nghost), stat=alloc_err)
          endif
          if (alloc_err > 0) call fatal_error ('wgrid', 'Could not allocate memory for gx,gy,gz', .true.)
!
          open (lun_output, FILE=trim(directory_collect)//'/'//file, FORM='unformatted', status='replace')
          t_sp = real (t)
        endif

        call collect_grid(x(:mxout1), y(:myout1), z(:mzout1), gx, gy, gz)
        if (lroot) then
          write (lun_output) t_sp, gx, gy, gz, dx, dy, dz
          write (lun_output) dx, dy, dz
          write (lun_output) Lx, Ly, Lz
        endif
        
        call collect_grid(dx_1(:mxout1), dy_1(:myout1), dz_1(:mzout1), gx, gy, gz)
        if (lroot) write (lun_output) gx, gy, gz
        
        call collect_grid(dx_tilde(:mxout1), dy_tilde(:myout1), dz_tilde(:mzout1), gx, gy, gz)
        if (lroot) then
          write (lun_output) gx, gy, gz
          close (lun_output)
       endif
     endif
!
    endsubroutine wgrid
!***********************************************************************
    subroutine input_grid_single(x,y,z,dx,dy,dz,Lx,Ly,Lz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde)
!
!  Read grid in single precision
!
!  23-oct-13/MR: derived from input_grid
!
      real(KIND=rkind4),                intent(OUT) :: dx,dy,dz,Lx,Ly,Lz
      real(KIND=rkind4), dimension (mx),intent(OUT) :: x,dx_1,dx_tilde
      real(KIND=rkind4), dimension (my),intent(OUT) :: y,dy_1,dy_tilde
      real(KIND=rkind4), dimension (mz),intent(OUT) :: z,dz_1,dz_tilde
!
      integer :: io_err
      real(KIND=rkind4) :: t_sp   ! t in single precision for backwards compatibility
!
      read(lun_input) t_sp,x,y,z,dx,dy,dz
      read(lun_input) dx,dy,dz
      read(lun_input, iostat=io_err) Lx,Ly,Lz
      if (io_err < 0) then
        ! End-Of-File: give notification that box dimensions are not read.
        ! This should only happen when reading old files.
        ! We should allow this for the time being.
        call warning ('input_grid', "Lx,Ly,Lz are not yet in grid.dat")
      else
        read(lun_input) dx_1,dy_1,dz_1
        read(lun_input) dx_tilde,dy_tilde,dz_tilde
      endif
!
    endsubroutine input_grid_single
!***********************************************************************
    subroutine input_grid_double(x,y,z,dx,dy,dz,Lx,Ly,Lz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde)
!
!  Read grid in double precision
!
!  23-oct-13/MR: derived from input_grid
!
      real(KIND=rkind8),                intent(OUT) :: dx,dy,dz,Lx,Ly,Lz
      real(KIND=rkind8), dimension (mx),intent(OUT) :: x,dx_1,dx_tilde
      real(KIND=rkind8), dimension (my),intent(OUT) :: y,dy_1,dy_tilde
      real(KIND=rkind8), dimension (mz),intent(OUT) :: z,dz_1,dz_tilde
!
      integer :: io_err
      real(KIND=rkind8) :: t_sp   ! t in double precision for backwards compatibility
!
      read(lun_input) t_sp,x,y,z,dx,dy,dz
      read(lun_input) dx,dy,dz
      read(lun_input, iostat=io_err) Lx,Ly,Lz
      if (io_err < 0) then
        ! End-Of-File: give notification that box dimensions are not read.
        ! This should only happen when reading old files.
        ! We should allow this for the time being.
        call warning ('input_grid_double', "Lx,Ly,Lz are not yet in grid.dat")
      else
        read(lun_input) dx_1,dy_1,dz_1
        read(lun_input) dx_tilde,dy_tilde,dz_tilde
      endif
!
    endsubroutine input_grid_double
!***********************************************************************
    subroutine rgrid (file)
!
!  Read processor-local part of grid coordinates.
!
!  21-jan-02/wolf: coded
!  15-jun-03/axel: Lx,Ly,Lz are now read in from file (Tony noticed the mistake)
!  24-oct-13/MR  : handling of reading from different precision introduced
!  28-oct-13/MR  : added overwriting of grid.dat if restart from different precision
!   3-mar-15/MR  : calculation of d[xyz]2_bound added: contain twice the distances of
!                  three neighbouring points from the boundary point
!  15-apr-15/MR  : automatic detection of precision added; doesn't work for old files 
!                  which don't contain Lx, Ly, Lz etc.
!
      use File_io, only: file_size
!
      character (len=*) :: file
!
      integer :: datasize, filesize
      integer, parameter :: nrec=5
      logical :: lotherprec
!
      real(KIND=rkind8), dimension(:), allocatable :: xdb,ydb,zdb,dx_1db,dy_1db,dz_1db,dx_tildedb,dy_tildedb,dz_tildedb
      real(KIND=rkind4), dimension(:), allocatable :: xsg,ysg,zsg,dx_1sg,dy_1sg,dz_1sg,dx_tildesg,dy_tildesg,dz_tildesg

      real(KIND=rkind8) :: dxdb,dydb,dzdb,Lxdb,Lydb,Lzdb
      real(KIND=rkind4) :: dxsg,dysg,dzsg,Lxsg,Lysg,Lzsg

      open(lun_input,FILE=trim(directory)//'/'//file,FORM='unformatted',status='old')

      lotherprec=.false.

      if (lread_from_other_prec) then

        datasize = 3*(mx+my+mz) + 10
        filesize = file_size(trim(directory)//'/'//file) - 8*nrec    ! length without record marker
!
        if (kind(x)==rkind4) then
          lotherprec = filesize/=4*datasize
          if (lotherprec) then
            allocate(xdb(mx),ydb(my),zdb(mz),dx_1db(mx),dy_1db(my),dz_1db(mz),dx_tildedb(mx),dy_tildedb(my),dz_tildedb(mz))
            call input_grid(xdb,ydb,zdb,dxdb,dydb,dzdb,Lxdb,Lydb,Lzdb, &
                            dx_1db,dy_1db,dz_1db,dx_tildedb,dy_tildedb,dz_tildedb)
            x=xdb; y=ydb; z=zdb; dx=dxdb; dy=dydb; dz=dzdb
            Lx=Lxdb; Ly=Lydb; Lz=Lzdb; dx_1=dx_1db; dy_1=dy_1db; dz_1=dz_1db
            dx_tilde=dx_tildedb; dy_tilde=dy_tildedb; dz_tilde=dz_tildedb
          endif
        elseif (kind(x)==rkind8) then
          lotherprec = filesize/=8*datasize
          if (lotherprec) then
            allocate(xsg(mx),ysg(my),zsg(mz),dx_1sg(mx),dy_1sg(my),dz_1sg(mz),dx_tildesg(mx),dy_tildesg(my),dz_tildesg(mz))
            call input_grid(xsg,ysg,zsg,dxsg,dysg,dzsg,Lxsg,Lysg,Lzsg, &
                            dx_1sg,dy_1sg,dz_1sg,dx_tildesg,dy_tildesg,dz_tildesg)
            x=xsg; y=ysg; z=zsg; dx=dxsg; dy=dysg; dz=dzsg
            Lx=Lxsg; Ly=Lysg; Lz=Lzsg; dx_1=dx_1sg; dy_1=dy_1sg; dz_1=dz_1sg;
            dx_tilde=dx_tildesg; dy_tilde=dy_tildesg; dz_tilde=dz_tildesg
          endif
        endif
      endif
      if (.not.lotherprec) &
        call input_grid(x,y,z,dx,dy,dz,Lx,Ly,Lz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde)
      
      close(lun_input)
!
      if (lotherprec) call wgrid(file,lwrite=.true.)         ! perhaps not necessary
!
!  debug output
!
      if (ip<=4.and.lroot) then
        print*,'rgrid: Lx,Ly,Lz=',Lx,Ly,Lz
        print*,'rgrid: dx,dy,dz=',dx,dy,dz
      endif
!
    endsubroutine rgrid
!***********************************************************************
    subroutine wproc_bounds(file)
!
!  Export processor boundaries to file.
!
!  22-Feb-2012/PABourdin: adapted from io_dist
!  27-nov-2020/ccyang: make the file single
!
      character(len=*), intent(in) :: file
!
      integer :: ierr
!
!  Only one process is needed.
!
      if (.not. lroot) return
!
!  Write proc[xyz]_bounds.
!
      open (lun_output, FILE=file, FORM='unformatted', IOSTAT=ierr, status='replace')
      if (ierr /= 0) call fatal_error("wproc_bounds", "Cannot open " // trim(file))
      write (lun_output) procx_bounds
      write (lun_output) procy_bounds
      write (lun_output) procz_bounds
      close (lun_output)
!
    endsubroutine wproc_bounds
!***********************************************************************
    subroutine rproc_bounds(file)
!
!   Import processor boundaries from file.
!
!   10-jul-08/kapelrud: coded
!   16-Feb-2012/Bourdin.KIS: rewritten
!
      character (len=*) :: file
!
      real(KIND=rkind4), dimension(0:nprocx):: procx_boundssg
      real(KIND=rkind4), dimension(0:nprocy):: procy_boundssg
      real(KIND=rkind4), dimension(0:nprocz):: procz_boundssg
!
      real(KIND=rkind8), dimension(0:nprocx):: procx_boundsdb
      real(KIND=rkind8), dimension(0:nprocy):: procy_boundsdb
      real(KIND=rkind8), dimension(0:nprocz):: procz_boundsdb
!
      open(lun_input,FILE=file,FORM='unformatted',status='old')

      if (lread_from_other_prec) then
        if (kind(x)==rkind4) then
          call input_proc_bounds(procx_boundsdb,procy_boundsdb,procz_boundsdb)
          procx_bounds=procx_boundsdb; procy_bounds=procy_boundsdb; procz_bounds=procz_boundsdb
        elseif (kind(x)==rkind8) then
          call input_proc_bounds(procx_boundssg,procy_boundssg,procz_boundssg)
          procx_bounds=procx_boundssg; procy_bounds=procy_boundssg; procz_bounds=procz_boundssg
        endif
      else
        call input_proc_bounds(procx_bounds,procy_bounds,procz_bounds)
      endif

      close(lun_output)
!
    endsubroutine rproc_bounds
!***********************************************************************
    subroutine input_proc_bounds_double(procx_bounds,procy_bounds,procz_bounds)
!
!   Import processor boundaries from file.in double precision
!
!   23-oct-13/MR: derivced from rproc_bounds
!
      real(KIND=rkind8), dimension(0:nprocx), intent(OUT):: procx_bounds
      real(KIND=rkind8), dimension(0:nprocy), intent(OUT):: procy_bounds
      real(KIND=rkind8), dimension(0:nprocz), intent(OUT):: procz_bounds
!
      read(lun_input) procx_bounds
      read(lun_input) procy_bounds
      read(lun_input) procz_bounds
!
    endsubroutine input_proc_bounds_double
!***********************************************************************
    subroutine input_proc_bounds_single(procx_bounds,procy_bounds,procz_bounds)
!
!   Import processor boundaries from file.in single precision
!
!   23-oct-13/MR: derived from rproc_bounds
!
      real(KIND=rkind4), dimension(0:nprocx), intent(OUT):: procx_bounds
      real(KIND=rkind4), dimension(0:nprocy), intent(OUT):: procy_bounds
      real(KIND=rkind4), dimension(0:nprocz), intent(OUT):: procz_bounds
!
      read(lun_input) procx_bounds
      read(lun_input) procy_bounds
      read(lun_input) procz_bounds
!
    endsubroutine input_proc_bounds_single
!***********************************************************************
    function read_part_double(lun,a,omit) result (ios)
!
!  Reads part of array, stored on disk, defined by range of variable indices ivar_omit(1):ivar_omit(2)
!  into array a. Uses dummy buffer omit of size (:,:,:,1:ivar_omit(2)-ivar_omit(1)+1).
!
!   31-mar-21/MR: coded
!
      integer, intent(IN) :: lun
      real(KIND=rkind8), dimension(:,:,:,:), intent(OUT) :: a, omit
      integer :: ios

      if (ivar_omit(1)==1) then
        read (lun,iostat=ios) omit,a
      elseif (ivar_omit(1)>1) then
        read (lun,iostat=ios) a(:,:,:,:ivar_omit(1)-1),omit,a(:,:,:,ivar_omit(1):)
      else
        read (lun,iostat=ios) a
      endif

    endfunction read_part_double
!***********************************************************************
    function read_part_single(lun,a,omit) result (ios)
!
!  Reads part of array, stored on disk, defined by range of variable indices ivar_omit(1):ivar_omit(2)
!  into array a. Uses dummy buffer omit of size (:,:,:,1:ivar_omit(2)-ivar_omit(1)+1).
!
!   31-mar-21/MR: coded
!
      integer, intent(IN) :: lun
      real(KIND=rkind4), dimension(:,:,:,:), intent(OUT) :: a, omit
      integer :: ios

      if (ivar_omit(1)==1) then
        read (lun,iostat=ios) omit,a
      elseif (ivar_omit(1)>1) then
        read (lun,iostat=ios) a(:,:,:,:ivar_omit(1)-1),omit,a(:,:,:,ivar_omit(1):)
      else
        read (lun,iostat=ios) a
      endif

    endfunction read_part_single
!***********************************************************************
endmodule Io
