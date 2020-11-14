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
!
module Io
!
  use Cdata
  use Cparam, only: intlen, fnlen, max_int
  use File_io, only: file_exists
  use HDF5_IO, only: output_dim
  use Messages, only: fatal_error, outlog, warning, svn_id
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
      if (lread_from_other_prec) &
        call warning('register_io','Reading from other precision not implemented')

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
    subroutine output_snap(a,nv,file)
!
!  Write snapshot file, always write time and mesh, could add other things.
!
!  11-apr-97/axel: coded
!  13-Dec-2011/Bourdin.KIS: reworked
!  13-feb-2014/MR: made 'file' optional, 'a' assumed-shape (for downsampled output);
!                  moved donwsampling stuff to snapshot
!
      use Mpicomm, only: start_serialize, end_serialize
      use General, only: get_range_no
!
      character (len=*), intent(in), optional :: file
      integer, intent(in) :: nv
      real, dimension (:,:,:,:), intent(in) :: a
!
      real :: t_sp   ! t in single precision for backwards compatibility
      integer :: io_err
      logical :: lerror
!
      t_sp = real (t)
      if (lroot .and. (ip <= 8)) print *, 'output_vect: nv =', nv
!
      if (lserial_io) call start_serialize
      if (present(file)) then
        open (lun_output, FILE=trim(directory_snap)//'/'//file, FORM='unformatted', IOSTAT=io_err, status='replace')
        lerror = outlog (io_err, 'openw', file, dist=lun_output, location='output_snap')
      endif
!
      if (lwrite_2d) then
        if (nx == 1) then
          write (lun_output, IOSTAT=io_err) a(l1,:,:,:)
        elseif (ny == 1) then
          write (lun_output, IOSTAT=io_err) a(:,m1,:,:)
        elseif (nz == 1) then
          write (lun_output, IOSTAT=io_err) a(:,:,n1,:)
        else
          io_err = 0
          call fatal_error ('output_snap', 'lwrite_2d used for 3D simulation!')
        endif
      elseif (ldownsampl) then
        write (lun_output, IOSTAT=io_err) a(firstind(1):l2:downsampl(1), &
                                            firstind(2):m2:downsampl(2), &
                                            firstind(3):n2:downsampl(3), :)
      else
        write (lun_output, IOSTAT=io_err) a
      endif
!
      lerror = outlog(io_err, 'main data')
!
!  Write shear at the end of x,y,z,dx,dy,dz.
!  At some good moment we may want to treat deltay like with
!  other modules and call a corresponding i/o parameter module.
!
      if (lshear) then
        write (lun_output, IOSTAT=io_err) t_sp, x, y, z, dx, dy, dz, deltay
        lerror = outlog(io_err, 'additional data and deltay')
      else
        write (lun_output, IOSTAT=io_err) t_sp, x, y, z, dx, dy, dz
        lerror = outlog(io_err, 'additional data')
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
      integer :: io_err
      logical :: lerror
!
      if (persist_initialized) then
        write (lun_output, iostat=io_err) id_block_PERSISTENT
        lerror = outlog (io_err, 'id_block_PERSISTENT')
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
    subroutine output_part_snap(ipar, ap, mv, nv, file, label, ltruncate)
!
!  Write particle snapshot file, always write mesh and time.
!
!  23-Oct-2018/PABourdin: adapted from output_snap
!
      integer, intent(in) :: mv, nv
      integer, dimension (mv), intent(in) :: ipar
      real, dimension (mv,mparray), intent(in) :: ap
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: ltruncate
!
      integer :: io_err
      logical :: lerror
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = real (t)
      open(lun_output,FILE=trim(directory_dist)//'/'//file,FORM='unformatted', IOSTAT=io_err, status='new')
      lerror = outlog (io_err, 'openw', file, dist=lun_output, location='output_part_snap')
!
      ! write the number of particles present at the processor
      write(lun_output, IOSTAT=io_err) nv
      lerror = outlog (io_err, 'number of particles per processor')
      if (nv/=0) then
        ! write particle index
        write(lun_output, IOSTAT=io_err) ipar(1:nv)
        ! write particle data
        write(lun_output, IOSTAT=io_err) ap(1:nv,:)
      endif
      lerror = outlog(io_err, 'main data')
!
      ! write time and grid parameters
      write(lun_output, IOSTAT=io_err) t_sp, x, y, z, dx, dy, dz
      lerror = outlog(io_err, 'additional data')
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
      character (len=fnlen) :: filename
      integer :: io_err
      logical :: lerror
      real :: t_sp
!
      call output_stalker ('', 0, nv, (/ 0.0 /), num)
!
      filename = trim(directory_dist)//'/particles_stalker.dat'
      open (lun_output, file=filename, form='unformatted', IOSTAT=io_err, position='append')
      lerror = outlog (io_err, 'openw', filename, dist=lun_output, location='output_stalker_init')
!
      ! write the time, number, and ID of stalked particles at this processor
      t_sp = t
      write (lun_output, IOSTAT=io_err) t_sp, nv
      lerror = outlog (io_err, 'number of stalker particles per processor')
      if (nv >= 1) then
        write (lun_output, IOSTAT=io_err) ID
        lerror = outlog (io_err, 'stalker particles ID')
      endif
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
      integer :: io_err
      logical :: lerror
!
      if (loptest (lfinalize)) then
        ! deallocate temporary stalker particle space
        if (pos > 0) then
          write (lun_output, IOSTAT=io_err) stalk_data
          lerror = outlog (io_err, 'stalker particles data')
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
!  03-May-2019/PABourdin: adapted from output_snap_finalize
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
      integer :: io_err
      logical :: lerror
!
      if (.not. lroot) return
!
      open(lun_output,FILE=trim(directory_snap)//'/'//trim(file),FORM='unformatted', IOSTAT=io_err, status='new')
      lerror = outlog (io_err, 'openw', file, dist=lun_output, location='output_pointmass')
      write(lun_output, IOSTAT=io_err) mv
      lerror = outlog (io_err, 'number of pointmasses')
      if (mv > 0) then
        write(lun_output, IOSTAT=io_err) fq
        lerror = outlog (io_err, 'main data')
      endif
      write(lun_output, IOSTAT=io_err) t
      lerror = outlog (io_err, 'additional data')
      close(lun_output)
      if (ip<=10) print*,'written pointmass snapshot ', trim (file)
!
    endsubroutine output_pointmass
!***********************************************************************
    subroutine output_slice_position
!
!  'data/procN/slice_position.dat' is distributed, but may not be synchronized
!  on I/O error (-> dist=0) as this would make it disfunctional; correct a posteriori if necessary.
!
!  27-Oct-2018/PABourdin: cleaned up
!
      integer :: io_err
!
      open (lun_output, file=trim(directory)//'/slice_position.dat', STATUS='unknown', IOSTAT=io_err)
!
      if (outlog (io_err, 'open', trim(directory)//'/slice_position.dat')) return
!
      write (lun_output, '(l5,i5," XY")', IOSTAT=io_err) lwrite_slice_xy, iz_loc
      if (outlog (io_err, 'write lwrite_slice_xy,iz_loc')) return
!
      write (lun_output, '(l5,i5," XY2")', IOSTAT=io_err) lwrite_slice_xy2, iz2_loc
      if (outlog (io_err, 'write lwrite_slice_xy2,iz2_loc')) return
!
      write (lun_output, '(l5,i5," XY3")', IOSTAT=io_err) lwrite_slice_xy3, iz3_loc
      if (outlog (io_err, 'write lwrite_slice_xy3,iz3_loc')) return
!
      write (lun_output, '(l5,i5," XY4")', IOSTAT=io_err) lwrite_slice_xy4, iz4_loc
      if (outlog (io_err, 'write lwrite_slice_xy4,iz4_loc')) return
!
      write (lun_output, '(l5,i5," XZ")', IOSTAT=io_err) lwrite_slice_xz, iy_loc
      if (outlog (io_err, 'write lwrite_slice_xz,iy_loc')) return
!
      write (lun_output, '(l5,i5," XZ2")', IOSTAT=io_err) lwrite_slice_xz2, iy2_loc
      if (outlog (io_err, 'write lwrite_slice_xz2,iy2_loc')) return
!
      write (lun_output, '(l5,i5," YZ")', IOSTAT=io_err) lwrite_slice_yz, ix_loc
      if (outlog (io_err, 'write lwrite_slice_yz,ix_loc')) return
!
      close (lun_output, IOSTAT=io_err)
      if (outlog (io_err,'close')) return
!
    endsubroutine output_slice_position
!***********************************************************************
    subroutine output_slice(lwrite, time, label, suffix, pos, grid_pos, data)
!
!  append to a slice file
!
!  12-nov-02/axel: coded
!  26-jun-06/anders: moved to slices
!  22-sep-07/axel: changed Xy to xy2, to be compatible with Mac
!  28-Oct-2018/PABourdin: cleaned up, moved to IO module
!
      logical, intent(in) :: lwrite
      real, intent(in) :: time
      character (len=*), intent(in) :: label, suffix
      real, intent(in) :: pos
      integer, intent(in) :: grid_pos
      real, dimension (:,:), pointer :: data
!
      character (len=fnlen) :: filename
      integer :: io_err
!
      if (.not. lwrite .or. .not. associated(data)) return
!
!  files data/procN/slice*.* are distributed and will be synchronized a-posteriori on I/O error
!
      filename = trim(directory)//'/slice_'//trim(label)//'.'//trim(suffix)
      open (lun_output, file=filename, form='unformatted', position='append', &
          IOSTAT=io_err)
      if (outlog (io_err, 'open', filename, dist=1)) return
!
      write (lun_output, IOSTAT=io_err) data, time, pos
      if (outlog (io_err, 'write data,time,pos')) return
!
      close (lun_output, IOSTAT=io_err)
      if (outlog (io_err, 'close')) continue
!
    endsubroutine output_slice
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
      integer :: io_err
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
        open (lun_output, FILE=trim (directory_snap)//'/'//filename, FORM='unformatted', &
              IOSTAT=io_err, status='replace')
        init_write_persist = outlog (io_err, 'openw persistent file', &
                             trim (directory_snap)//'/'//filename, location='init_write_persist' )
        filename = ""
      endif
!
      if (lroot .and. (ip <= 9)) write (*,*) 'begin persistent block'
      write (lun_output, iostat=io_err) id_block_PERSISTENT
      init_write_persist = outlog (io_err, 'id_block_PERSISTENT')
      persist_initialized = .not. init_write_persist
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
      integer :: io_err
!
      write_persist_id = .true.
      if (.not. persist_initialized) write_persist_id = init_write_persist()
      if (.not. persist_initialized) return
!
      if (persist_last_id /= id) then
        write (lun_output, iostat=io_err) id
        write_persist_id = outlog (io_err, 'persistent ID '//label)
        persist_last_id = id
      else
        write_persist_id = .false.
      endif
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
      integer :: io_err
!
      write_persist_logical_0D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=io_err) value
      write_persist_logical_0D = outlog (io_err, 'persistent '//label)
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
      integer :: io_err
!
      write_persist_logical_1D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=io_err) value
      write_persist_logical_1D = outlog (io_err, 'persistent '//label)
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
      integer :: io_err
!
      write_persist_int_0D = .true.
      if (.not. persist_initialized) return
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=io_err) value
      write_persist_int_0D = outlog (io_err, 'persistent '//label)
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
      integer :: io_err
!
      write_persist_int_1D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=io_err) value
      write_persist_int_1D = outlog (io_err, 'persistent '//label)
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
      integer :: io_err
!
      write_persist_real_0D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=io_err) value
      write_persist_real_0D = outlog (io_err, 'persistent '//label)
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
      integer :: io_err
!
      write_persist_real_1D = .true.
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=io_err) value
      write_persist_real_1D = outlog (io_err, 'persistent '//label)
!
    endfunction write_persist_real_1D
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
!                if lreset_tstart=F: cancel program
!                                =T, tstart unspecified: use minimum time of all var.dat
!                                =T, tstart specified: use this value
!
      use Mpicomm, only: start_serialize, end_serialize, mpibcast_real, mpiallreduce_or, &
                         stop_it, mpiallreduce_min, MPI_COMM_WORLD
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real(KIND=rkind4), dimension (mx,my,mz,nv), intent(out) :: a
!
      real(KIND=rkind4),                 intent(out) :: dx, dy, dz, deltay
      real(KIND=rkind4), dimension (mx), intent(out) :: x
      real(KIND=rkind4), dimension (my), intent(out) :: y
      real(KIND=rkind4), dimension (mz), intent(out) :: z
!
      real(KIND=rkind4) :: t_sp, t_sgl
      real :: t_test   ! t in single precision for backwards compatibility
!
      integer :: io_err
      logical :: lerror, ltest
      logical :: lreset_tstart
!
      lreset_tstart = .false.
!
      if (lserial_io) call start_serialize
      open (lun_input, FILE=trim(directory_snap)//'/'//file, FORM='unformatted', &
            IOSTAT=io_err, status='old')
      lerror = outlog (io_err, "openr snapshot data", trim(directory_snap)//'/'//file, &
                       location='read_snap_single')
!      if (ip<=8) print *, 'read_snap_single: open, mx,my,mz,nv=', mx, my, mz, nv
      if (lwrite_2d) then
        if (nx == 1) then
          read (lun_input, IOSTAT=io_err) a(4,:,:,:)
        elseif (ny == 1) then
          read (lun_input, IOSTAT=io_err) a(:,4,:,:)
        elseif (nz == 1) then
          read (lun_input, IOSTAT=io_err) a(:,:,4,:)
        else
          io_err = 0
          call fatal_error ('read_snap_single', 'lwrite_2d used for 3-D simulation!')
        endif
      else
!
!  Possibility of reading data with different numbers of ghost zones.
!  In that case, one must regenerate the mesh with luse_oldgrid=T.
!
        if (nghost_read_fewer==0) then
          read (lun_input, IOSTAT=io_err) a
        elseif (nghost_read_fewer>0) then
          read (lun_input, IOSTAT=io_err) &
              a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                1+nghost_read_fewer:my-nghost_read_fewer, &
                1+nghost_read_fewer:mz-nghost_read_fewer,:)
!
!  The following 3 possibilities allow us to replicate 1-D data input
!  in x (nghost_read_fewer=-1), y (-2), or z (-3) correspondingly.
!
        elseif (nghost_read_fewer==-1) then
          read (lun_input, IOSTAT=io_err) a(:,1:1+nghost*2,1:1+nghost*2,:)
          a=spread(spread(a(:,m1,n1,:),2,my),3,mz)
        elseif (nghost_read_fewer==-2) then
          read (lun_input, IOSTAT=io_err) a(1:1+nghost*2,:,1:1+nghost*2,:)
          a=spread(spread(a(l1,:,n1,:),1,mx),3,mz)
        elseif (nghost_read_fewer==-3) then
          read (lun_input, IOSTAT=io_err) a(1:1+nghost*2,1:1+nghost*2,:,:)
          a=spread(spread(a(l1,m1,:,:),1,mx),2,my)
        else
          call fatal_error('read_snap_single','nghost_read_fewer must be >=0')
        endif
      endif
      lerror = outlog (io_err, 'main data')

      if (ip <= 8) print *, 'read_snap_single: read ', file
      if (mode == 1) then
!
!  Check whether we want to read deltay from snapshot.
!
        if (lshear) then
          read (lun_input, IOSTAT=io_err) t_sp, x, y, z, dx, dy, dz, deltay
          lerror = outlog (io_err, 'additional data + deltay')
        else
          if (nghost_read_fewer==0) then
            read (lun_input, IOSTAT=io_err) t_sp, x, y, z, dx, dy, dz
          elseif (nghost_read_fewer>0) then
            read (lun_input, IOSTAT=io_err) t_sp
          endif
          lerror = outlog (io_err, 'additional data')
        endif
!
!  Verify consistency of the snapshots regarding their timestamp,
!  unless lreset_tstart=T, in which case we reset all times to tstart.
!
        if (.not.lreset_tstart.or.tstart==impossible) then
!
          t_test = t_sp
          call mpibcast_real(t_test,comm=MPI_COMM_WORLD)
          call mpiallreduce_or(t_test /= t_sp .and. .not.lread_from_other_prec &
                               .or. abs(t_test-t_sp)>1.e-6,ltest,MPI_COMM_WORLD)
!
!  If timestamps deviate at any processor
!
          if (ltest) then
            if (lreset_tstart) then
!
!  If reset of tstart enabled and tstart unspecified, use minimum of all t_sp
!
              call mpiallreduce_min(t_sp,t_sgl,MPI_COMM_WORLD)
              tstart = t_sgl
              if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT. Using t=', tstart,'.'
            else
              write (*,*) 'ERROR: '//trim(directory_snap)//'/'//trim(file)// &
                          ' IS INCONSISTENT: t=', t_sp
              call stop_it('read_snap_single')
            endif
          else
            tstart = t_sp
          endif
!
        endif
!
!  Set time or overwrite it by a given value.
!
        if (lreset_tstart) then
          t = tstart
        else
          t = t_sp
        endif
!
!  Verify the read values for x, y, z, and t.
!
        if (ip <= 3) print *, 'read_snap_single: x=', x
        if (ip <= 3) print *, 'read_snap_single: y=', y
        if (ip <= 3) print *, 'read_snap_single: z=', z
        if (ip <= 3) print *, 'read_snap_single: t=', t
!
      endif
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
!                if lreset_tstart=F: cancel program
!                                =T, tstart unspecified: use minimum time of all var.dat 
!                                                        for start
!                                =T, tstart specified: use this value
!                             
      use Mpicomm, only: start_serialize, end_serialize, mpibcast_real, mpiallreduce_or, &
                         stop_it, mpiallreduce_min, MPI_COMM_WORLD
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real(KIND=rkind8), dimension (mx,my,mz,nv), intent(out) :: a
!
      real(KIND=rkind8) :: t_sp, t_dbl

      real(KIND=rkind8), intent(out) :: dx, dy, dz, deltay
      real(KIND=rkind8), dimension (mx), intent(out) :: x
      real(KIND=rkind8), dimension (my), intent(out) :: y
      real(KIND=rkind8), dimension (mz), intent(out) :: z

      real :: t_test   ! t in single precision for backwards compatibility
      integer :: io_err
      logical :: lerror,ltest
      logical :: lreset_tstart
!
      lreset_tstart = .false.
!
      if (lserial_io) call start_serialize
      open (lun_input, FILE=trim(directory_snap)//'/'//file, FORM='unformatted', &
            IOSTAT=io_err, status='old')
      lerror = outlog (io_err, "openr snapshot data", trim(directory_snap)//'/'//file, &
                       location='read_snap_double')
!      if (ip<=8) print *, 'read_snap_double: open, mx,my,mz,nv=', mx, my, mz, nv
      if (lwrite_2d) then
        if (nx == 1) then
          read (lun_input, IOSTAT=io_err) a(4,:,:,:)
        elseif (ny == 1) then
          read (lun_input, IOSTAT=io_err) a(:,4,:,:)
        elseif (nz == 1) then
          read (lun_input, IOSTAT=io_err) a(:,:,4,:)
        else
          io_err = 0
          call fatal_error ('read_snap_double', 'lwrite_2d used for 3-D simulation!')
        endif
      else
!
!  Possibility of reading data with different numbers of ghost zones.
!  In that case, one must regenerate the mesh with luse_oldgrid=T.
!
        if (nghost_read_fewer==0) then
          read (lun_input, IOSTAT=io_err) a
        elseif (nghost_read_fewer>0) then
          read (lun_input, IOSTAT=io_err) &
              a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                1+nghost_read_fewer:my-nghost_read_fewer, &
                1+nghost_read_fewer:mz-nghost_read_fewer,:)
!
!  The following 3 possibilities allow us to replicate 1-D data input
!  in x (nghost_read_fewer=-1), y (-2), or z (-3) correspondingly.
!
        elseif (nghost_read_fewer==-1) then
          read (lun_input, IOSTAT=io_err) a(:,1:1+nghost*2,1:1+nghost*2,:)
          a=spread(spread(a(:,m1,n1,:),2,my),3,mz)
        elseif (nghost_read_fewer==-2) then
          read (lun_input, IOSTAT=io_err) a(1:1+nghost*2,:,1:1+nghost*2,:)
          a=spread(spread(a(l1,:,n1,:),1,mx),3,mz)
        elseif (nghost_read_fewer==-3) then
          read (lun_input, IOSTAT=io_err) a(1:1+nghost*2,1:1+nghost*2,:,:)
          a=spread(spread(a(l1,m1,:,:),1,mx),2,my)
        else
          call fatal_error('read_snap_double','nghost_read_fewer must be >=0')
        endif
      endif
      lerror = outlog (io_err, 'main data')

      if (ip <= 8) print *, 'read_snap: read ', file
      if (mode == 1) then
!
!  Check whether we want to read deltay from snapshot.
!
        if (lshear) then
          read (lun_input, IOSTAT=io_err) t_sp, x, y, z, dx, dy, dz, deltay
          lerror = outlog (io_err, 'additional data + deltay')
        else
          if (nghost_read_fewer==0) then
            read (lun_input, IOSTAT=io_err) t_sp, x, y, z, dx, dy, dz
          elseif (nghost_read_fewer>0) then
            read (lun_input, IOSTAT=io_err) t_sp
          endif
          lerror = outlog (io_err, 'additional data')
        endif
!
!  Verify consistency of the snapshots regarding their timestamp,
!  unless lreset_tstart=T, in which case we reset all times to tstart.
!
        if (.not.lreset_tstart.or.tstart==impossible) then
!
          t_test = t_sp
          call mpibcast_real(t_test,comm=MPI_COMM_WORLD)
          call mpiallreduce_or(t_test /= t_sp .and. .not.lread_from_other_prec &
                               .or. abs(t_test-t_sp)>1.e-6,ltest,MPI_COMM_WORLD)
!
!  If timestamp deviates at any processor
!
          if (ltest) then
            if (lreset_tstart) then
!
!  If reset of tstart enabled and tstart unspecified, use minimum of all t_sp
!
              call mpiallreduce_min(t_sp,t_dbl,MPI_COMM_WORLD)
              tstart = t_dbl
              if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT. Using t=', tstart, '.'
            else
              write (*,*) 'ERROR: '//trim(directory_snap)//'/'//trim(file)// &
                          ' IS INCONSISTENT: t=', t_sp
              call stop_it('read_snap_double')
            endif
          else
            tstart = t_sp
          endif
!
        endif
!
!  Set time or overwrite it by a given value.
!
        if (lreset_tstart) then
          t = tstart
        else
          t = t_sp
        endif
!
!  Verify the read values for x, y, z, and t.
!
        if (ip <= 3) print *, 'read_snap_double: x=', x
        if (ip <= 3) print *, 'read_snap_double: y=', y
        if (ip <= 3) print *, 'read_snap_double: z=', z
        if (ip <= 3) print *, 'read_snap_double: t=', t
!
      endif
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
!
      close (lun_input)
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
!  25-Oct-2018/PABourdin: apadpted and moved to IO module
!
      use Mpicomm, only: mpireduce_sum_int
!
      integer, intent(in) :: mv
      integer, dimension (mv), intent(out) :: ipar
      real, dimension (mv,mparray), intent(out) :: ap
      integer, intent(out) :: nv, npar_total
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label
!
      integer :: io_err
      logical :: lerror
!
      open (1, FILE=trim(directory_dist)//'/'//file, FORM='unformatted', IOSTAT=io_err, status='old')
      lerror = outlog (io_err, "openr snapshot data", trim(directory_snap)//'/'//file, location='read_snap_double')
!
      ! Read number of particles for this processor.
      read (1, IOSTAT=io_err) nv
      lerror = outlog (io_err, 'number of particles per processor')
      if (nv > 0) then
        ! Read particle number identifier and then particle data.
        read (1, IOSTAT=io_err) ipar(1:nv)
        read (1, IOSTAT=io_err) ap(1:nv,:)
        lerror = outlog (io_err, 'main data')
      endif
      close (1)
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
      integer :: io_err
      logical :: lerror
!
      if (lroot) then
        open(lun_input,FILE=trim(directory_snap)//'/'//trim(file),FORM='unformatted')
        lerror = outlog (io_err, 'openr', file, dist=lun_input, location='input_pointmass')
        read(lun_input) mv_in
        lerror = outlog (io_err, 'number of pointmasses')
        if (mv_in /= mv) call fatal_error("","")
        if (mv_in /= 0) read(lun_input) fq
        lerror = outlog (io_err, 'main data')
        close(lun_input)
        if (ip<=8) print*, 'read pointmass snapshot', trim (file)
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
      integer :: io_err
!
      init_read_persist = .true.
!
      if (present (file)) then
        if (lroot) init_read_persist = .not. file_exists (trim (directory_snap)//'/'//file)
        call mpibcast_logical (init_read_persist,comm=MPI_COMM_WORLD)
        if (init_read_persist) return
      endif
!
      if (present (file)) then
        close (lun_input)
        open (lun_input, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', IOSTAT=io_err, status='old')
        init_read_persist = outlog (io_err, 'openr persistent data',file,location='init_read_persist')
      endif
!
      if (lroot .and. (ip <= 9)) write (*,*) 'begin persistent block'
      init_read_persist = .false.
      persist_initialized = .true.
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
!
!  17-Feb-2012/Bourdin.KIS: coded
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
      read (lun_input, iostat=io_err) id
      if (lcatch_error) then
        if (io_err /= 0) then
          id = -max_int
          read_persist_id = .true.
        else
          read_persist_id = .false.
        endif
      else
        read_persist_id = outlog (io_err, 'persistent ID '//label,lcont=.true.)
      endif
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
      integer :: io_err
!
      read (lun_input, iostat=io_err) value
      read_persist_logical_0D = outlog(io_err, 'persistent '//label,lcont=.true.)
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
      integer :: io_err
!
      read (lun_input, iostat=io_err) value
      read_persist_logical_1D = outlog(io_err, 'persistent '//label,lcont=.true.)
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
      integer :: io_err
!
      read (lun_input, iostat=io_err) value
      read_persist_int_0D = outlog(io_err, 'persistent '//label,lcont=.true.)
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
      integer :: io_err
!
      read (lun_input, iostat=io_err) value
      read_persist_int_1D = outlog(io_err, 'persistent '//label,lcont=.true.)
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
      integer :: io_err
      real(KIND=rkind8) :: vdb
      real(KIND=rkind4) :: vsg
!
      if (lread_from_other_prec) then
        if (kind(value)==rkind4) then
          read (lun_input, iostat=io_err) vdb
          value=vdb
        elseif (kind(value)==rkind8) then
          read (lun_input, iostat=io_err) vsg
          value=vsg
        endif
      else
        read (lun_input, iostat=io_err) value
      endif

      read_persist_real_0D = outlog(io_err, 'persistent '//label,lcont=.true.)
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
      integer :: io_err
      real(KIND=rkind8), dimension(:), allocatable :: vdb
      real(KIND=rkind4), dimension(:), allocatable :: vsg
!
      if (lread_from_other_prec) then
        if (kind(value)==rkind4) then
          allocate(vdb(size(value)))
          read (lun_input, iostat=io_err) vdb
          value=vdb
        elseif (kind(value)==rkind8) then
          allocate(vsg(size(value)))
          read (lun_input, iostat=io_err) vsg
          value=vsg
        endif
      else
        read (lun_input, iostat=io_err) value
      endif
!
      read_persist_real_1D = outlog(io_err, 'persistent '//label, lcont=.true.)
!
    endfunction read_persist_real_1D
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
      integer :: io_err
      logical :: lerror
!
      if (lserial_io) call start_serialize
      open(lun_output,FILE=trim(directory_snap)//'/'//file,FORM='unformatted',IOSTAT=io_err,status='replace')
      lerror = outlog(io_err,"openw",file,location='output_globals')
!
      if (lwrite_2d) then
        if (nx==1) then
          write(lun_output,IOSTAT=io_err) a(4,:,:,:)
        elseif (ny==1) then
          write(lun_output,IOSTAT=io_err) a(:,4,:,:)
        elseif (nz==1) then
          write(lun_output,IOSTAT=io_err) a(:,:,4,:)
        else
          io_err=0
          call fatal_error('output_globals','lwrite_2d used for 3-D simulation!')
        endif
      else
        write(lun_output,IOSTAT=io_err) a
      endif
      lerror = outlog(io_err,"data block")
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
      integer :: io_err
      logical :: lerror
      real(KIND=rkind8), dimension(:,:,:,:), allocatable :: adb
      real(KIND=rkind4), dimension(:,:,:,:), allocatable :: asg

!
      if (lserial_io) call start_serialize
!
      open(lun_input,FILE=trim(directory_snap)//'/'//file,FORM='unformatted',IOSTAT=io_err,status='old')
      lerror = outlog(io_err,"openr globals",file,location='input_globals')

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
      integer :: io_err
      logical :: lerror

      if (lwrite_2d) then
        if (nx==1) then
          read(lun_input,IOSTAT=io_err) a(4,:,:,:)
        elseif (ny==1) then
          read(lun_input,IOSTAT=io_err) a(:,4,:,:)
        elseif (nz==1) then
          read(lun_input,IOSTAT=io_err) a(:,:,4,:)
        else
          io_err=0
          call fatal_error('input_globals','lwrite_2d used for 3-D simulation!')
        endif
      else
        read(lun_input,IOSTAT=io_err) a
      endif
      lerror = outlog(io_err,"data block",location='read_globals_double')
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
      integer :: io_err
      logical :: lerror

      if (lwrite_2d) then
        if (nx==1) then
          read(lun_input,IOSTAT=io_err) a(4,:,:,:)
        elseif (ny==1) then
          read(lun_input,IOSTAT=io_err) a(:,4,:,:)
        elseif (nz==1) then
          read(lun_input,IOSTAT=io_err) a(:,:,4,:)
        else
          io_err=0
          call fatal_error('input_globals','lwrite_2d used for 3-D simulation!')
        endif
      else
        read(lun_input,IOSTAT=io_err) a
      endif
      lerror = outlog(io_err,"data block",location='read_globals_single')
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
      integer :: io_err
      logical :: lerror
!
      call parse_filename(file,dir,fpart)
      if (dir == '.') call safe_character_assign(dir,directory_snap)
!
      open(lun_output,FILE=trim(dir)//'/'//flist,POSITION='append',IOSTAT=io_err)
      ! file not distributed?, backskipping enabled
      lerror = outlog(io_err,"openw",trim(dir)//'/'//flist,dist=-lun_output, &
                      location='log_filename_to_file')
      write(lun_output,'(A,1x,e16.8)',IOSTAT=io_err) trim(fpart), t
      lerror = outlog(io_err,"fpart", trim(dir)//'/'//flist)
      close(lun_output)
!
      if (lcopysnapshots_exp) then
        if (lroot) then
          open(lun_output,FILE=trim(datadir)//'/move-me.list',POSITION='append',IOSTAT=io_err)
          ! file not distributed?, backskipping enabled
          lerror = outlog(io_err,"openw",trim(datadir)//'/move-me.list',dist=-lun_output, &
                          location='log_filename_to_file')
          write(lun_output,'(A)',IOSTAT=io_err) trim(fpart)
          lerror = outlog(io_err,"fpart")
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
!                  to be able to output coordinate vectors with coordinates differing from
!                  mx,my,mz
!  25-Aug-2016/PABourdin: now a global "grid.dat" is always written from all IO modules
!
      use Mpicomm, only: collect_grid
      use General, only: loptest
!
      character (len=*) :: file
      integer, optional :: mxout,myout,mzout
      logical, optional :: lwrite
!
      integer :: mxout1,myout1,mzout1
      integer :: io_err
      logical :: lerror
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (loptest(lwrite,.not.luse_oldgrid)) then
        if (present(mzout)) then
          mxout1=mxout
          myout1=myout
          mzout1=mzout
        else
          mxout1=mx
          myout1=my
          mzout1=mz
        endif
!
        t_sp = real (t)

        open(lun_output,FILE=trim(directory)//'/'//file,FORM='unformatted',IOSTAT=io_err,status='replace')
        if (io_err /= 0) call fatal_error('wgrid', &
            "Cannot open " // trim(file) // " (or similar) for writing" // &
            " -- is data/ visible from all nodes?", .true.)
        lerror = outlog(io_err,"openw",trim(directory)//'/'//file,location='wgrid')
        write(lun_output,IOSTAT=io_err) t_sp,x(1:mxout1),y(1:myout1),z(1:mzout1),dx,dy,dz
        lerror = outlog(io_err,"main data block")
        write(lun_output,IOSTAT=io_err) dx,dy,dz
        lerror = outlog(io_err,"dx,dy,dz")
        write(lun_output,IOSTAT=io_err) Lx,Ly,Lz
        lerror = outlog(io_err,"Lx,Ly,Lz")
        write(lun_output,IOSTAT=io_err) dx_1(1:mxout1),dy_1(1:myout1),dz_1(1:mzout1)
        lerror = outlog(io_err,"dx_1,dy_1,dz_1")
        write(lun_output,IOSTAT=io_err) dx_tilde(1:mxout1),dy_tilde(1:myout1),dz_tilde(1:mzout1)
        lerror = outlog(io_err,"dx_tilde,dy_tilde,dz_tilde")
        close(lun_output,IOSTAT=io_err)
        lerror = outlog(io_err,'close')
!
        if (lyang) return      ! grid collection only needed on Yin grid, as grids are identical

        ! write also a global "data/allprocs/grid.dat"
        if (lroot) then
          allocate (gx(nxgrid+2*nghost), gy(nygrid+2*nghost), gz(nzgrid+2*nghost), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('wgrid', 'Could not allocate memory for gx,gy,gz', .true.)
!
          open (lun_output, FILE=trim(directory_collect)//'/'//file, FORM='unformatted', status='replace')
          t_sp = real (t)
        endif

        call collect_grid (x, y, z, gx, gy, gz)
        if (lroot) then
          write (lun_output) t_sp, gx, gy, gz, dx, dy, dz
          write (lun_output) dx, dy, dz
          write (lun_output) Lx, Ly, Lz
        endif

        call collect_grid (dx_1, dy_1, dz_1, gx, gy, gz)
        if (lroot) write (lun_output) gx, gy, gz

        call collect_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
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
!  28-oct-13/MR: added lcont and location parameters to calls of outlog where appropriate
!
      real(KIND=rkind4),                intent(OUT) :: dx,dy,dz,Lx,Ly,Lz
      real(KIND=rkind4), dimension (mx),intent(OUT) :: x,dx_1,dx_tilde
      real(KIND=rkind4), dimension (my),intent(OUT) :: y,dy_1,dy_tilde
      real(KIND=rkind4), dimension (mz),intent(OUT) :: z,dz_1,dz_tilde

      integer :: io_err
      logical :: lerror
      real(KIND=rkind4) :: t_sp   ! t in single precision for backwards compatibility
!
      read(lun_input,IOSTAT=io_err) t_sp,x,y,z,dx,dy,dz
      lerror = outlog(io_err,"main data block",location='input_grid_single', lcont=.true.)
      read(lun_input,IOSTAT=io_err) dx,dy,dz
      lerror = outlog(io_err,"dx,dy,dz",lcont=.true.)
      read(lun_input,IOSTAT=io_err) Lx,Ly,Lz
      if (io_err < 0) then
        ! End-Of-File: give notification that box dimensions are not read.
        ! This should only happen when reading old files.
        ! We should allow this for the time being.
        call warning ('input_grid', "Lx,Ly,Lz are not yet in grid.dat")
      else
        lerror = outlog(io_err,"Lx,Ly,Lz")
        read(lun_input,IOSTAT=io_err) dx_1,dy_1,dz_1
        lerror = outlog(io_err,"dx_1,dy_1,dz_1")
        read(lun_input,IOSTAT=io_err) dx_tilde,dy_tilde,dz_tilde
        lerror = outlog(io_err,"dx_tilde,dy_tilde,dz_tilde")
      endif
!
    endsubroutine input_grid_single
!***********************************************************************
    subroutine input_grid_double(x,y,z,dx,dy,dz,Lx,Ly,Lz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde)
!
!  Read grid in double precision
!
!  23-oct-13/MR: derived from input_grid
!  28-oct-13/MR: added lcont and location parameters to calls of outlog where appropriate
!
      real(KIND=rkind8),                intent(OUT) :: dx,dy,dz,Lx,Ly,Lz
      real(KIND=rkind8), dimension (mx),intent(OUT) :: x,dx_1,dx_tilde
      real(KIND=rkind8), dimension (my),intent(OUT) :: y,dy_1,dy_tilde
      real(KIND=rkind8), dimension (mz),intent(OUT) :: z,dz_1,dz_tilde

      integer :: io_err
      logical :: lerror
      real(KIND=rkind8) :: t_sp   ! t in single precision for backwards compatibility
!
      read(lun_input,IOSTAT=io_err) t_sp,x,y,z,dx,dy,dz
      lerror = outlog(io_err,"main data block",location='input_grid_double',lcont=.true.)
      read(lun_input,IOSTAT=io_err) dx,dy,dz
      lerror = outlog(io_err,"dx,dy,dz",lcont=.true.)
      read(lun_input,IOSTAT=io_err) Lx,Ly,Lz
      if (io_err < 0) then
        ! End-Of-File: give notification that box dimensions are not read.
        ! This should only happen when reading old files.
        ! We should allow this for the time being.
        call warning ('input_grid_double', "Lx,Ly,Lz are not yet in grid.dat")
      else
        lerror = outlog(io_err,"Lx,Ly,Lz",lcont=.true.)
        read(lun_input,IOSTAT=io_err) dx_1,dy_1,dz_1
        lerror = outlog(io_err,"dx_1,dy_1,dz_1",lcont=.true.)
        read(lun_input,IOSTAT=io_err) dx_tilde,dy_tilde,dz_tilde
        lerror = outlog(io_err,"dx_tilde,dy_tilde,dz_tilde",lcont=.true.)
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
!  15-apr-15/MR  : automatic detection of precision added
!
      use File_io, only: file_size
!
      character (len=*) :: file
!
      integer :: io_err, datasize, filesize
      integer, parameter :: nrec=5
      logical :: lerror, lotherprec
!
      real(KIND=rkind8), dimension(:), allocatable :: xdb,ydb,zdb,dx_1db,dy_1db,dz_1db,dx_tildedb,dy_tildedb,dz_tildedb
      real(KIND=rkind4), dimension(:), allocatable :: xsg,ysg,zsg,dx_1sg,dy_1sg,dz_1sg,dx_tildesg,dy_tildesg,dz_tildesg

      real(KIND=rkind8) :: dxdb,dydb,dzdb,Lxdb,Lydb,Lzdb
      real(KIND=rkind4) :: dxsg,dysg,dzsg,Lxsg,Lysg,Lzsg

      open(lun_input,FILE=trim(directory)//'/'//file,FORM='unformatted',IOSTAT=io_err,status='old')
      if (io_err /= 0) call fatal_error('rgrid', &
          "Cannot open " // trim(file) // " (or similar) for reading" // &
          " -- is data/ visible from all nodes?",.true.)
      lerror = outlog(io_err,'openr',file,location='rgrid')

      if (lread_from_other_prec) then

        datasize = 3*(mx+my+mz) + 10
        filesize = file_size(trim(directory)//'/'//file) - 8*nrec
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
      else
        call input_grid(x,y,z,dx,dy,dz,Lx,Ly,Lz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde)
      endif
      close(lun_input,IOSTAT=io_err)
      lerror = outlog(io_err,'close')
!
      if (lread_from_other_prec.and.lotherprec) call wgrid(file,lwrite=.true.)         ! perhaps not necessary
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
    subroutine output_profile(fname, coord, a, type, lsave_name, lhas_ghost)
!
!  Writes a profile to a file.
!
!  10-jul-05/axel: coded
!  07-Nov-2018/PABourdin: moved to IO modules
!
      use General, only: loptest
!
      real, dimension(:) :: coord, a
      character (len=*) :: fname
      character :: type
      logical, optional :: lsave_name, lhas_ghost
!
      integer :: pos, np, io_err
      logical lerror
!
!  If within a loop, do this only for the first step (indicated by lwrite_prof).
!
      if (lwrite_prof) then
        ! Write profile.
        open(lun_output,file=trim(directory)//'/'//type//'prof_'//trim(fname)//'.dat',position='append',IOSTAT=io_err)
        lerror = outlog(io_err,"openw",trim(directory)//type//'/prof_'//trim(fname)//'.dat',location='output_profile')
        np = size(coord)
        do pos=1, np
          write(lun_output,*,IOSTAT=io_err) coord(pos), a(pos)
          lerror = outlog(io_err,'output_profile')
        enddo
        close(lun_output,IOSTAT=io_err)
        lerror = outlog(io_err,'close')
!
        ! Add file name to list of profiles if lsave_name is true.
        if (loptest(lsave_name)) then
          open(lun_output,file=trim(directory)//'/'//type//'prof_list.dat',position='append',IOSTAT=io_err)
          lerror = outlog(io_err,"openw",trim(directory)//type//'/prof_list.dat',location='output_profile')
          write(lun_output,*,IOSTAT=io_err) fname
          lerror = outlog(io_err,'output_profile')
          close(lun_output,IOSTAT=io_err)
          lerror = outlog(io_err,'close')
        endif
      endif
!
    endsubroutine output_profile
!***********************************************************************
    subroutine input_profile(fname, type, a, np, lhas_ghost)
!
!  Read a profile from a file (taken from stratification, for example).
!
!  15-may-15/axel: coded
!  05-Nov-2018/PABourdin: generalized to any direction and moved to IO module
!
      character (len=*), intent(in) :: fname
      character, intent(in) :: type
      integer, intent(in) :: np
      real, dimension(np), intent(out) :: a
      logical, optional :: lhas_ghost
!
      integer :: pos, io_err
      logical lerror
      real, dimension(np) :: coord
!
      ! Read profile.
      open(lun_input,file=trim(directory)//type//'/prof_'//trim(fname)//'.dat',IOSTAT=io_err)
      lerror = outlog(io_err,"openr",trim(directory)//type//'/prof_'//trim(fname)//'.dat',location='input_profile')
      do pos=1, np
        read(lun_input,*,IOSTAT=io_err) coord(pos), a(pos)
        lerror = outlog(io_err,'input_profile')
      enddo
      close(lun_input,IOSTAT=io_err)
      lerror = outlog(io_err,'close')
!
!  Should we check that coord == z for type == 'z'?
!
    endsubroutine input_profile
!***********************************************************************
    subroutine output_average_1D(path, label, nc, name, data, time, lbinary, lwrite, header)
!
!   Output 1D average to a file.
!
!   16-Nov-2018/PABourdin: coded
!
      character (len=*), intent(in) :: path, label
      integer, intent(in) :: nc
      character (len=fmtlen), dimension(nc), intent(in) :: name
      real, dimension(:,:), intent(in) :: data
      real, intent(in) :: time
      logical, intent(in) :: lbinary, lwrite
      real, dimension(:), optional, intent(in) :: header
!
      character (len=fnlen) :: filename
      integer :: io_err
      logical :: lerror
!
      if (.not. lwrite .or. (nc <= 0)) return
!
      filename = trim(path) // '/' // trim(label) // 'averages.dat'
      if (trim (label) == 'phi_z') filename = trim(path) // '/phizaverages.dat'
      if (lbinary) then
        open(lun_output, file=filename, form='unformatted', position='append', IOSTAT=io_err)
        lerror = outlog(io_err,"openw",trim(filename),location='output_average_1D')
        if (present (header)) then
          write(lun_output,IOSTAT=io_err) header
          lerror = outlog(io_err,'1D-average header')
        endif
        write(lun_output,IOSTAT=io_err) time
        lerror = outlog(io_err,'1D-average time')
        write(lun_output,IOSTAT=io_err) data(:,1:nc)
        lerror = outlog(io_err,'1D-average data')
        close(lun_output)
      else
        open(lun_output, file=filename, position='append', IOSTAT=io_err)
        lerror = outlog(io_err,"openw",trim(filename),location='output_average_1D')
        if (present (header)) then
          write(lun_output,'(1p,8e14.5e3)',IOSTAT=io_err) header
          lerror = outlog(io_err,'1D-average header')
        endif
        write(lun_output,'(1pe12.5)',IOSTAT=io_err) time
        lerror = outlog(io_err,'1D-average time')
        write(lun_output,'(1p,8e14.5e3)',IOSTAT=io_err) data(:,1:nc)
        lerror = outlog(io_err,'1D-average data')
        close(lun_output)
      endif
!
    endsubroutine output_average_1D
!***********************************************************************
    subroutine output_average_2D(path, label, nc, name, data, time, lbinary, lwrite, header)
!
!   Output average to a file.
!
!   16-Nov-2018/PABourdin: coded
!
      character (len=*), intent(in) :: path, label
      integer, intent(in) :: nc
      character (len=fmtlen), dimension(nc), intent(in) :: name
      real, dimension(:,:,:), intent(in) :: data
      real, intent(in) :: time
      logical, intent(in) :: lbinary, lwrite
      real, dimension(:), optional, intent(in) :: header
!
      character (len=fnlen) :: filename
      integer :: ia
      integer :: io_err
      logical :: lerror
!
      if (.not. lwrite .or. (nc <= 0)) return
!
      filename = trim(path) // '/' // trim(label) // 'averages.dat'
      if (lbinary) then
        open(lun_output, file=filename, form='unformatted', position='append', IOSTAT=io_err)
        lerror = outlog(io_err,"openw",trim(filename),location='output_average_2D')
        if (present (header)) then
          write(lun_output,IOSTAT=io_err) header
          lerror = outlog(io_err,'2D-average header')
        endif
        write(lun_output,IOSTAT=io_err) time
        lerror = outlog(io_err,'2D-average time')
        if (label == 'z') then
          write(lun_output,IOSTAT=io_err) ( data(ia,:,:), ia=1, nc )
          lerror = outlog(io_err,'2D-average z-data')
        else
          write(lun_output,IOSTAT=io_err) data(:,:,1:nc)
          lerror = outlog(io_err,'2D-average data')
        endif
        close(lun_output)
      else
        open(lun_output, file=filename, position='append', IOSTAT=io_err)
        lerror = outlog(io_err,"openw",trim(filename),location='output_average_2D')
        if (present (header)) then
          write(lun_output,'(1p,8e14.5e3)',IOSTAT=io_err) header
          lerror = outlog(io_err,'2D-average header')
        endif
        write(lun_output,'(1pe12.5)',IOSTAT=io_err) time
        lerror = outlog(io_err,'2D-average time')
        if (label == 'z') then
          write(lun_output,'(1p,8e14.5e3)',IOSTAT=io_err) ( data(ia,:,:), ia=1, nc )
          lerror = outlog(io_err,'2D-average z-data')
        else
          write(lun_output,'(1p,8e14.5e3)',IOSTAT=io_err) data(:,:,1:nc)
          lerror = outlog(io_err,'2D-average data')
        endif
        close(lun_output)
      endif
!
    endsubroutine output_average_2D
!***********************************************************************
    subroutine output_average_phi(path, number, nr, nc, name, data, time, r, dr)
!
!   Output phi average to a file with these records:
!   1) nr_phiavg, nz_phiavg, nvars, nprocz
!   2) t, r_phiavg, z_phiavg, dr, dz
!   3) data
!   4) len(labels),labels
!
!   27-Nov-2014/PABourdin: cleaned up code from write_phiaverages
!   25-Nov-2018/PABourdin: coded
!
      use General, only: safe_character_append
!
      character (len=*), intent(in) :: path, number
      integer, intent(in) :: nr, nc
      character (len=fmtlen), dimension(nc), intent(in) :: name
      real, dimension(:,:,:,:), intent(in) :: data
      real, intent(in) :: time
      real, dimension(nr), intent(in) :: r
      real, intent(in) :: dr
!
      character (len=fnlen) :: filename
      character (len=1024) :: labels
      integer :: pos
      integer :: io_err
      logical :: lerror
!
      if (.not. lroot .or. (nc <= 0)) return
!
      filename = 'PHIAVG' // trim(number)
      open(lun_output, file=trim(path)//'/averages/'//trim(filename), form='unformatted', position='append', IOSTAT=io_err)
      lerror = outlog(io_err,"openw",trim(path)//'/averages/'//trim(filename),location='output_average_phi')
      write(lun_output,IOSTAT=io_err) nr, nzgrid, nc, nprocz
      lerror = outlog(io_err,'phi-average header')
      write(lun_output,IOSTAT=io_err) time, r, z(n1)+(/(pos*dz, pos=0, nzgrid-1)/), dr, dz
      lerror = outlog(io_err,'phi-average time')
      ! note: due to passing data as implicit-size array,
      ! the indices (0:nz) are shifted to (1:nz+1),
      ! so that we have to write only the portion (2:nz+1).
      ! data has to be repacked to avoid writing an array temporary
      ! ngrs: writing an array temporary outputs corrupted data on copson
      write(lun_output,IOSTAT=io_err) pack(data(:,2:nz+1,:,1:nc), .true.)
      lerror = outlog(io_err,'phi-average data')
      labels = trim(name(1))
      if (nc >= 2) then
        do pos = 2, nc
          call safe_character_append (labels, ",", trim(name(pos)))
        enddo
      endif
      write(lun_output,IOSTAT=io_err) len(labels), labels
      lerror = outlog(io_err,'phi-average labels')
      close(lun_output)
!
      ! append filename to list
      open (lun_output, FILE=trim(datadir)//'/averages/phiavg.files', POSITION='append', IOSTAT=io_err)
      lerror = outlog(io_err,"openw",trim(datadir)//'/averages/phiavg.files',location='output_average_phi')
      write (lun_output,'(A)',IOSTAT=io_err) trim(filename)
      lerror = outlog(io_err,'phi-average list')
      close (lun_output)
!
    endsubroutine output_average_phi
!***********************************************************************
    subroutine wproc_bounds(file)
!
!   Export processor boundaries to file.
!
!   10-jul-08/kapelrud: coded
!   16-Feb-2012/Bourdin.KIS: rewritten
!
      character (len=*) :: file
!
      integer :: io_err
      logical :: lerror
!
      open(lun_output,FILE=file,FORM='unformatted',IOSTAT=io_err,status='replace')
      lerror = outlog(io_err,"openw",file,location='wproc_bounds')
      write(lun_output,IOSTAT=io_err) procx_bounds
      lerror = outlog(io_err,'procx_bounds')
      write(lun_output,IOSTAT=io_err) procy_bounds
      lerror = outlog(io_err,'procy_bounds')
      write(lun_output,IOSTAT=io_err) procz_bounds
      lerror = outlog(io_err,'procz_bounds')
      close(lun_output,IOSTAT=io_err)
      lerror = outlog(io_err,'close')
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
      integer :: io_err
      logical :: lerror
      real(KIND=rkind4), dimension(0:nprocx):: procx_boundssg
      real(KIND=rkind4), dimension(0:nprocy):: procy_boundssg
      real(KIND=rkind4), dimension(0:nprocz):: procz_boundssg
!
      real(KIND=rkind8), dimension(0:nprocx):: procx_boundsdb
      real(KIND=rkind8), dimension(0:nprocy):: procy_boundsdb
      real(KIND=rkind8), dimension(0:nprocz):: procz_boundsdb
!
      open(lun_input,FILE=file,FORM='unformatted',IOSTAT=io_err,status='old')

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

      close(lun_output,IOSTAT=io_err)
      lerror = outlog(io_err,'close')
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

      integer :: io_err
      logical :: lerror
!
      read(lun_input,IOSTAT=io_err) procx_bounds
      lerror = outlog(io_err,'procx_bounds')
      read(lun_input,IOSTAT=io_err) procy_bounds
      lerror = outlog(io_err,'procy_bounds')
      read(lun_input,IOSTAT=io_err) procz_bounds
      lerror = outlog(io_err,'procz_bounds')
!
    endsubroutine input_proc_bounds_double
!***********************************************************************
    subroutine input_proc_bounds_single(procx_bounds,procy_bounds,procz_bounds)
!
!   Import processor boundaries from file.in single precision
!
!   23-oct-13/MR: derivced from rproc_bounds
!
      real(KIND=rkind4), dimension(0:nprocx), intent(OUT):: procx_bounds
      real(KIND=rkind4), dimension(0:nprocy), intent(OUT):: procy_bounds
      real(KIND=rkind4), dimension(0:nprocz), intent(OUT):: procz_bounds

      integer :: io_err
      logical :: lerror
!
      read(lun_input,IOSTAT=io_err) procx_bounds
      lerror = outlog(io_err,'procx_bounds')
      read(lun_input,IOSTAT=io_err) procy_bounds
      lerror = outlog(io_err,'procy_bounds')
      read(lun_input,IOSTAT=io_err) procz_bounds
      lerror = outlog(io_err,'procz_bounds')
!
    endsubroutine input_proc_bounds_single
!***********************************************************************
endmodule Io
