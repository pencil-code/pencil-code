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
  use HDF5_IO, only: output_dim
  use Messages, only: fatal_error, warning, svn_id
!
  implicit none
!
  include 'io.h'
  include 'record_types.h'
!
  interface write_persist
    module procedure write_persist_logical_0D
    module procedure write_persist_logical_1D
    module procedure write_persist_int_0D
    module procedure write_persist_int_1D
    module procedure write_persist_real_0D
    module procedure write_persist_real_1D
  endinterface
!
  interface read_persist
    module procedure read_persist_logical_0D
    module procedure read_persist_logical_1D
    module procedure read_persist_int_0D
    module procedure read_persist_int_1D
    module procedure read_persist_real_0D
    module procedure read_persist_real_1D
  endinterface
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
  ! define unique logical unit number for input and output calls
  integer :: lun_input=88
  integer :: lun_output=91
  ! set ireset_tstart to 1 or 2 to coordinate divergent timestamp
  integer :: MINT=1
  integer :: MAXT=2
!
  ! Indicates if IO is done distributed (each proc writes into a procdir)
  ! or collectively (eg. by specialized IO-nodes or by MPI-IO).
  logical :: lcollective_IO=.false.
  character (len=labellen) :: IO_strategy="dist"
!
  logical :: persist_initialized=.false.
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
      ldistribute_persist = .true.
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
      use General, only: get_range_no, ioptest
!
      real, dimension (:,:,:,:),  intent(IN) :: a
      integer,                    intent(IN) :: nv
      character (len=*), optional,intent(IN) :: file
!
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = real (t)
      if (lroot .and. (ip <= 8)) print *, 'output_vect: nv =', nv
!
      if (lserial_io) call start_serialize
      if (present(file)) then
        call delete_file(trim(directory_snap)//'/'//file)
        open (lun_output, FILE=trim(directory_snap)//'/'//file, FORM='unformatted', status='new')
      endif
!
      if (lwrite_2d) then
        if (nx == 1) then
          write (lun_output) a(l1,:,:,:)
        elseif (ny == 1) then
          write (lun_output) a(:,m1,:,:)
        elseif (nz == 1) then
          write (lun_output) a(:,:,n1,:)
        else
          call fatal_error ('output_snap', 'lwrite_2d used for 3D simulation!')
        endif
      else
        write (lun_output) a
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
    subroutine output_slice_position
!
!  'data/procN/slice_position.dat' is distributed, but may not be synchronized
!  on I/O error (-> dist=0) as this would make it disfunctional; correct a posteriori if necessary.
!
!  27-Oct-2018/PABourdin: cleaned up
!
      open (lun_output, file=trim(directory)//'/slice_position.dat', STATUS='unknown')
      write (lun_output, '(l5,i5," XY")') lwrite_slice_xy, iz_loc
      write (lun_output, '(l5,i5," XY2")') lwrite_slice_xy2, iz2_loc
      write (lun_output, '(l5,i5," XY3")') lwrite_slice_xy3, iz3_loc
      write (lun_output, '(l5,i5," XY4")') lwrite_slice_xy4, iz4_loc
      write (lun_output, '(l5,i5," XZ")') lwrite_slice_xz, iy_loc
      write (lun_output, '(l5,i5," XZ2")') lwrite_slice_xz2, iy2_loc
      write (lun_output, '(l5,i5," YZ")') lwrite_slice_yz, ix_loc
      close (lun_output)
!
    endsubroutine output_slice_position
!***********************************************************************
    subroutine output_slice(lwrite, time, label, suffix, grid, pos, grid_pos, data, ndim1, ndim2)
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
      real, dimension (:) :: grid
      integer, intent(in) :: pos, grid_pos
      integer, intent(in) :: ndim1, ndim2
      real, dimension (:,:), pointer :: data
!
      if (.not. lwrite .or. .not. associated(data)) return
!
!  files data/procN/slice*.* are distributed and will be synchronized a-posteriori on I/O error
!
      open (lun_output, file=trim(directory)//'/slice_'//trim(label)//'.'//trim(suffix), form='unformatted', position='append')
      write (lun_output) data, time, grid(pos)
      close (lun_output)
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
      persist_initialized = .true.
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
!
      use Mpicomm, only: start_serialize, end_serialize, mpibcast_real, mpiallreduce_or, &
                         stop_it, mpiallreduce_min, mpiallreduce_max, MPI_COMM_WORLD
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real(KIND=rkind4), dimension (mx,my,mz,nv), intent(out) :: a
!
      real(KIND=rkind4) :: t_sp, t_sgl

      real(KIND=rkind4),                 intent(out) :: dx, dy, dz, deltay
      real(KIND=rkind4), dimension (mx), intent(out) :: x
      real(KIND=rkind4), dimension (my), intent(out) :: y
      real(KIND=rkind4), dimension (mz), intent(out) :: z

      real :: t_test   ! t in single precision for backwards compatibility

      logical :: ltest
!
      if (lserial_io) call start_serialize
      open (lun_input, FILE=trim(directory_snap)//'/'//file, FORM='unformatted', status='old')
!      if (ip<=8) print *, 'read_snap_single: open, mx,my,mz,nv=', mx, my, mz, nv
      if (lwrite_2d) then
        if (nx == 1) then
          read (lun_input) a(4,:,:,:)
        elseif (ny == 1) then
          read (lun_input) a(:,4,:,:)
        elseif (nz == 1) then
          read (lun_input) a(:,:,4,:)
        else
          call fatal_error ('read_snap_single', 'lwrite_2d used for 3-D simulation!')
        endif
      else
!
!  Possibility of reading data with different numbers of ghost zones.
!  In that case, one must regenerate the mesh with luse_oldgrid=T.
!
        if (nghost_read_fewer==0) then
          read (lun_input) a
        elseif (nghost_read_fewer>0) then
          read (lun_input) &
              a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                1+nghost_read_fewer:my-nghost_read_fewer, &
                1+nghost_read_fewer:mz-nghost_read_fewer,:)
!
!  The following 3 possibilities allow us to replicate 1-D data input
!  in x (nghost_read_fewer=-1), y (-2), or z (-3) correspondingly.
!
        elseif (nghost_read_fewer==-1) then
          read (lun_input) a(:,1:1+nghost*2,1:1+nghost*2,:)
          a=spread(spread(a(:,m1,n1,:),2,my),3,mz)
        elseif (nghost_read_fewer==-2) then
          read (lun_input) a(1:1+nghost*2,:,1:1+nghost*2,:)
          a=spread(spread(a(l1,:,n1,:),1,mx),3,mz)
        elseif (nghost_read_fewer==-3) then
          read (lun_input) a(1:1+nghost*2,1:1+nghost*2,:,:)
          a=spread(spread(a(l1,m1,:,:),1,mx),2,my)
        else
          call fatal_error('read_snap_single','nghost_read_fewer must be >=0')
        endif
      endif

      if (ip <= 8) print *, 'read_snap_single: read ', file
      if (mode == 1) then
!
!  Check whether we want to read deltay from snapshot.
!
        if (lshear.and..not.lread_oldsnap_noshear) then
          read (lun_input) t_sp, x, y, z, dx, dy, dz, deltay
        else
          if (nghost_read_fewer==0) then
            read (lun_input) t_sp, x, y, z, dx, dy, dz
          elseif (nghost_read_fewer>0) then
            read (lun_input) t_sp
          endif
        endif
!
!  Verify consistency of the snapshots regarding their timestamp,
!  unless ireset_tstart=T, in which case we reset all times to tstart.
!
        if ((ireset_tstart == 0) .or. (tstart == impossible)) then
!
          t_test = t_sp
          call mpibcast_real(t_test,comm=MPI_COMM_WORLD)
          call mpiallreduce_or((t_test /= t_sp) .and. .not. lread_from_other_prec &
                                .or. (abs(t_test-t_sp) > 1.e-6),ltest,MPI_COMM_WORLD)
!
!  If timestamps deviate at any processor
!
          if (ltest) then
            if (ireset_tstart > 0) then
!
!  If reset of tstart enabled and tstart unspecified, use minimum of all t_sp
!
              if (ireset_tstart == MINT) then
                call mpiallreduce_min(t_sp,t_sgl,MPI_COMM_WORLD)
                if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                       ' Using (min) t=', t_sgl,'with ireset_tstart=', MINT,'.'
              elseif (ireset_tstart >= MAXT) then
                call mpiallreduce_max(t_sp,t_sgl,MPI_COMM_WORLD)
                if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                       ' Using (max) t=', t_sgl,'with ireset_tstart=', MAXT,'.'
              endif
              tstart = t_sgl
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
        if (ireset_tstart > 0) then
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
!                if ireset_tstart= 0: cancel program
!                                = MINT, tstart unspecified: use minimum time of all var.dat 
!                                                        for start
!                                = MAXT, tstart unspecified: use maximum time of all var.dat 
!                                                        for start
!                                > 0, tstart specified: use this value
!                             
      use Mpicomm, only: start_serialize, end_serialize, mpibcast_real, mpiallreduce_or, &
                         stop_it, mpiallreduce_min, mpiallreduce_max, MPI_COMM_WORLD
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
      logical :: ltest
!
      if (lserial_io) call start_serialize
      open (lun_input, FILE=trim(directory_snap)//'/'//file, FORM='unformatted', status='old')
!      if (ip<=8) print *, 'read_snap_double: open, mx,my,mz,nv=', mx, my, mz, nv
      if (lwrite_2d) then
        if (nx == 1) then
          read (lun_input) a(4,:,:,:)
        elseif (ny == 1) then
          read (lun_input) a(:,4,:,:)
        elseif (nz == 1) then
          read (lun_input) a(:,:,4,:)
        else
          call fatal_error ('read_snap_double', 'lwrite_2d used for 3-D simulation!')
        endif
      else
!
!  Possibility of reading data with different numbers of ghost zones.
!  In that case, one must regenerate the mesh with luse_oldgrid=T.
!
        if (nghost_read_fewer==0) then
          read (lun_input) a
        elseif (nghost_read_fewer>0) then
          read (lun_input) &
              a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                1+nghost_read_fewer:my-nghost_read_fewer, &
                1+nghost_read_fewer:mz-nghost_read_fewer,:)
!
!  The following 3 possibilities allow us to replicate 1-D data input
!  in x (nghost_read_fewer=-1), y (-2), or z (-3) correspondingly.
!
        elseif (nghost_read_fewer==-1) then
          read (lun_input) a(:,1:1+nghost*2,1:1+nghost*2,:)
          a=spread(spread(a(:,m1,n1,:),2,my),3,mz)
        elseif (nghost_read_fewer==-2) then
          read (lun_input) a(1:1+nghost*2,:,1:1+nghost*2,:)
          a=spread(spread(a(l1,:,n1,:),1,mx),3,mz)
        elseif (nghost_read_fewer==-3) then
          read (lun_input) a(1:1+nghost*2,1:1+nghost*2,:,:)
          a=spread(spread(a(l1,m1,:,:),1,mx),2,my)
        else
          call fatal_error('read_snap_double','nghost_read_fewer must be >=0')
        endif
      endif

      if (ip <= 8) print *, 'read_snap: read ', file
      if (mode == 1) then
!
!  Check whether we want to read deltay from snapshot.
!
        if (lshear.and..not.lread_oldsnap_noshear) then
          read (lun_input) t_sp, x, y, z, dx, dy, dz, deltay
        else
          if (nghost_read_fewer==0) then
            read (lun_input) t_sp, x, y, z, dx, dy, dz
          elseif (nghost_read_fewer>0) then
            read (lun_input) t_sp
          endif
        endif
!
!  Verify consistency of the snapshots regarding their timestamp,
!  unless ireset_tstart=T, in which case we reset all times to tstart.
!
        if ((ireset_tstart == 0) .or. (tstart == impossible)) then
!
          t_test = t_sp
          call mpibcast_real(t_test,comm=MPI_COMM_WORLD)
          call mpiallreduce_or((t_test /= t_sp) .and. .not. lread_from_other_prec &
                               .or. (abs(t_test-t_sp) > 1.e-6),ltest, MPI_COMM_WORLD)
!
!  If timestamp deviates at any processor
!
          if (ltest) then
            if (ireset_tstart > 0) then
!
!  If reset of tstart enabled and tstart unspecified, use minimum of all t_sp
!
              if (ireset_tstart == MINT) then
                call mpiallreduce_min(t_sp,t_dbl,MPI_COMM_WORLD)
                if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                       ' Using (min) t=', t_dbl,'with ireset_tstart=', MINT,'.'
              elseif (ireset_tstart >= MAXT) then
                call mpiallreduce_max(t_sp,t_dbl,MPI_COMM_WORLD)
                if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                       ' Using (max) t=', t_dbl,'with ireset_tstart=', MAXT,'.'
              endif
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
        if (ireset_tstart > 0) then
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
    subroutine input_part_snap(ipar, ap, mv, nv, npar_total, file, label)
!
!  Read particle snapshot file, mesh and time are read in 'input_snap'.
!
!  24-Oct-2018/PABourdin: apadpted and moved to IO module
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
      open (lun_input, FILE=trim(directory_dist)//'/'//file, FORM='unformatted')
      ! Read number of particles for this processor.
      read (lun_input) nv
      if (nv > 0) then
        ! Read particle number identifier and then particle data.
        read (lun_input) ipar(1:nv)
        read (lun_input) ap(1:nv,:)
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
      persist_initialized = .true.
!
    endfunction init_read_persist
!***********************************************************************
    logical function read_persist_id(label, id, lerror_prone)
!
!  Read persistent block ID from snapshot file.
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
    subroutine output_timeseries(data, data_im)
!
!  Append diagnostic data to a binary file.
!
!  01-Apr-2019/PABourdin: coded
!
      real, dimension(2*nname), intent(in) :: data, data_im
!
      ! dummy routine
!
    endsubroutine output_timeseries
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
          write(lun_output) a(4,:,:,:)
        elseif (ny==1) then
          write(lun_output) a(:,4,:,:)
        elseif (nz==1) then
          write(lun_output) a(:,:,4,:)
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
          read(lun_input) a(4,:,:,:)
        elseif (ny==1) then
          read(lun_input) a(:,4,:,:)
        elseif (nz==1) then
          read(lun_input) a(:,:,4,:)
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
          read(lun_input) a(4,:,:,:)
        elseif (ny==1) then
          read(lun_input) a(:,4,:,:)
        elseif (nz==1) then
          read(lun_input) a(:,:,4,:)
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
      write(lun_output,'(A)') trim(fpart)
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
    subroutine wgrid(file,mxout,myout,mzout)
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
      use General, only: ioptest
!
      character (len=*) :: file
      integer, optional :: mxout,myout,mzout
!
      integer :: mxout1,myout1,mzout1
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
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
        deallocate(gx,gy,gz) 
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
      if (lotherprec) call wgrid(file)         ! perhaps not necessary
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
    subroutine wdim_default_grid(file)
!
!  Write dimension to file.
!
!  02-Nov-2018/PABourdin: redesigned
!
      character (len=*), intent(in) :: file
!
      call output_dim (file, mx, my, mz, mxgrid, mygrid, mzgrid, mvar, maux, mglobal)
!
    endsubroutine wdim_default_grid
!***********************************************************************
    subroutine wdim_default(file, mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out)
!
!  Write dimension to file.
!
!  02-Nov-2018/PABourdin: redesigned
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out
!
      call output_dim (file, mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar, maux, mglobal)
!
    endsubroutine wdim_default
!***********************************************************************
    subroutine wdim(file, mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out)
!
!  Write dimension to file.
!
!   8-sep-01/axel: adapted to take my_out,mz_out
!   4-oct-16/MR: added optional parameters mvar_out,maux_out
!  02-Nov-2018/PABourdin: redesigned, moved to IO modules
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out
!
      call output_dim (file, mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out, mglobal)
!
    endsubroutine wdim
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
      integer :: pos, np
!
!  If within a loop, do this only for the first step (indicated by lwrite_prof).
!
      if (lwrite_prof) then
        ! Write profile.
        open(lun_output,file=trim(directory)//'/'//type//'prof_'//trim(fname)//'.dat',position='append')
        np = size(coord)
        do pos=1, np
          write(lun_output,*) coord(pos), a(pos)
        enddo
        close(lun_output)
!
        ! Add file name to list of profiles if lsave_name is true.
        if (loptest(lsave_name)) then
          open(lun_output,file=trim(directory)//'/'//type//'prof_list.dat',position='append')
          write(lun_output,*) fname
          close(lun_output)
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
      integer :: pos
      real, dimension(np) :: coord
!
      ! Read profile.
      open(lun_input,file=trim(directory)//type//'/prof_'//trim(fname)//'.dat')
      do pos=1, np
        read(lun_input,*) coord(pos), a(pos)
      enddo
      close(lun_input)
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
!
      if (.not. lwrite .or. (nc <= 0)) return
!
      filename = trim(path) // '/' // trim(label) // 'averages.dat'
      if (trim (label) == 'phi_z') filename = trim(path) // '/phizaverages.dat'
      if (lbinary) then
        open(lun_output, file=filename, form='unformatted', position='append')
        if (present (header)) write(lun_output) header
        write(lun_output) time
        write(lun_output) data(:,1:nc)
        close(lun_output)
      else
        open(lun_output, file=filename, position='append')
        if (present (header)) write(lun_output,'(1p,8e14.5e3)') header
        write(lun_output,'(1pe12.5)') time
        write(lun_output,'(1p,8e14.5e3)') data(:,1:nc)
        close(lun_output)
      endif
!
    endsubroutine output_average_1D
!***********************************************************************
    subroutine output_average_1D_chunked(path, label, nc, name, data, full, time, lbinary, lwrite, header)
!
!   Output 1D chunked average to a file.
!
!   16-Nov-2018/PABourdin: coded
!
      character (len=*), intent(in) :: path, label
      integer, intent(in) :: nc
      character (len=fmtlen), dimension(nc), intent(in) :: name
      real, dimension(:,:,:), intent(in) :: data
      integer, intent(in) :: full
      real, intent(in) :: time
      logical, intent(in) :: lbinary, lwrite
      real, dimension(:), optional, intent(in) :: header
!
      if (present (header)) then
        call output_average_2D(path, label, nc, name, data, time, lbinary, lwrite, header)
      else
        call output_average_2D(path, label, nc, name, data, time, lbinary, lwrite)
      endif
!
    endsubroutine output_average_1D_chunked
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
!
      if (.not. lwrite .or. (nc <= 0)) return
!
      filename = trim(path) // '/' // trim(label) // 'averages.dat'
      if (lbinary) then
        open(lun_output, file=filename, form='unformatted', position='append')
        if (present (header)) write(lun_output) header
        write(lun_output) time
        if (label == 'z') then
          write(lun_output) ( data(ia,:,:), ia=1, nc )
        else
          write(lun_output) data(:,:,1:nc)
        endif
        close(lun_output)
      else
        open(lun_output, file=filename, position='append')
        if (present (header)) write(lun_output,'(1p,8e14.5e3)') header
        write(lun_output,'(1pe12.5)') time
        if (label == 'z') then
          write(lun_output,'(1p,8e14.5e3)') ( data(ia,:,:), ia=1, nc )
        else
          write(lun_output,'(1p,8e14.5e3)') data(:,:,1:nc)
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
!
      if (.not. lroot .or. (nc <= 0)) return
!
      filename = 'PHIAVG' // trim(number)
      open(lun_output, file=trim(path)//'/averages/'//trim(filename), form='unformatted', position='append')
      write(lun_output) nr, nzgrid, nc, nprocz
      write(lun_output) time, r, z(n1)+(/(pos*dz, pos=0, nzgrid-1)/), dr, dz
      ! note: due to passing data as implicit-size array,
      ! the indices (0:nz) are shifted to (1:nz+1),
      ! so that we have to write only the portion (2:nz+1).
      ! data has to be repacked to avoid writing an array temporary
      ! ngrs: writing an array temporary outputs corrupted data on copson
      write(lun_output) pack(data(:,2:nz+1,:,1:nc), .true.)
      labels = trim(name(1))
      if (nc >= 2) then
        do pos = 2, nc
          call safe_character_append (labels, ",", trim(name(pos)))
        enddo
      endif
      write(lun_output) len(labels), labels
      close(lun_output)
!
      ! append filename to list
      open (lun_output, FILE=trim(datadir)//'/averages/phiavg.files', POSITION='append')
      write (lun_output,'(A)') trim(filename)
      close (lun_output)
!
    endsubroutine output_average_phi
!***********************************************************************
    subroutine trim_average(path, plane, ngrid, nname)
!
!  Trim a 1D-average file for times past the current time.
!
!  25-apr-16/ccyang: coded
!  23-Nov-2018/PABourdin: moved to IO module
!
      character (len=*), intent(in) :: path, plane
      integer, intent(in) :: ngrid, nname
!
      real, dimension(:), allocatable :: tmp
      character(len=fnlen) :: filename
      integer :: pos, num_rec, ioerr, alloc_err
      integer :: lun_input = 84
      real :: time
!
      if (.not. lroot) return
      if ((ngrid <= 0) .or. (nname <= 0)) return
      filename = trim(path) // '/' // trim(plane) // 'averages.dat'
      if (.not. file_exists (filename)) return
!
      allocate (tmp(ngrid * nname), stat = alloc_err)
      if (alloc_err > 0) call fatal_error('trim_average', 'Could not allocate memory for averages')
!
      if (lwrite_avg1d_binary) then
        open(lun_input, file=filename, form='unformatted', action='readwrite')
      else
        open(lun_input, file=filename, action='readwrite')
      endif
!
      ! count number of records
      num_rec = 0
      ioerr = 0
      time = 0.0
      do while ((t > time) .and. (ioerr == 0))
        num_rec = num_rec + 1
        if (lwrite_avg1d_binary) then
          read(lun_input, iostat=ioerr) time
          read(lun_input, iostat=ioerr) tmp
        else
          read(lun_input, *, iostat=ioerr) time
          read(lun_input, *, iostat=ioerr) tmp
        endif
      enddo
      if (time == t) num_rec = num_rec - 1
!
      if ((time >= t) .and. (ioerr == 0)) then
        ! trim excess data at the end of the average file
        rewind(lun_input)
        if (num_rec > 0) then
          do pos = 1, num_rec
            if (lwrite_avg1d_binary) then
              read(lun_input, iostat=ioerr) time
              read(lun_input, iostat=ioerr) tmp
            else
              read(lun_input, *, iostat=ioerr) time
              read(lun_input, *, iostat=ioerr) tmp
            endif
          enddo
        endif
        endfile (lun_input)
        if (ip <= 10) print *, 'trim_average: trimmed '//trim(plane)//'-averages for t >= ', t
      endif
!
      close (lun_input)
      deallocate (tmp)
!
    endsubroutine trim_average
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
      call delete_file(file) 
      open(lun_output,FILE=file,FORM='unformatted',status='new')
      write(lun_output) procx_bounds
      write(lun_output) procy_bounds
      write(lun_output) procz_bounds
      close(lun_output)
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
!   23-oct-13/MR: derivced from rproc_bounds
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
endmodule Io
