! $Id$
!
!  I/O via the MPI v2 standard IO routines.
!  (storing data into one file, e.g. data/allprocs/var.dat)
!
!  The file written by output_snap() (and used e.g. for 'var.dat')
!  consists of the followinig records (not using record markers):
!    1. data(mxgrid,mygrid,mzgrid,nvar)
!    2. t(1), x(mxgrid), y(mygrid), z(mzgrid), dx(1), dy(1), dz(1)
!  Where nvar denotes the number of variables to be saved.
!  In the case of MHD with entropy, nvar is 8 for a 'var.dat' file.
!  Only outer ghost-layers are written, so mzlocal is between nz and mz,
!  depending on the corresponding ipz-layer.
!
!  To read these snapshots in IDL, the parameter allprocs needs to be set:
!  IDL> pc_read_var, obj=vars, /allprocs
!  or in a much more efficient way by reading into an array:
!  IDL> pc_read_var_raw, obj=data, tags=tags, grid=grid, /allprocs
!
!  20-Mar-2012/PABourdin: adapted from io_collect.f90 and io_collect_xy.f90
!  06-Oct-2015/PABourdin: reworked, should work now
!
module Io
!
  use Cdata
  use Cparam, only: intlen, fnlen, max_int
  use File_io, only: file_exists
  use Messages, only: fatal_error, svn_id, warning
!
  implicit none
!
  include 'io.h'
  include 'record_types.h'
  include 'mpif.h'
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
  ! define unique logical unit number for input and output calls
  integer :: lun_input=88
  integer :: lun_output=91
!
  ! Indicates if IO is done distributed (each proc writes into a procdir)
  ! or collectively (eg. by specialized IO-nodes or by MPI-IO).
  logical :: lcollective_IO=.true.
  character (len=labellen) :: IO_strategy="MPI-IO"
!
  logical :: lread_add=.true., lwrite_add=.true.
  logical :: persist_initialized=.false.
  integer :: persist_last_id=-max_int
!
  integer :: local_type, global_type, mpi_err, io_dims
  integer, parameter :: mpi_comm=MPI_COMM_WORLD
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer (kind=MPI_OFFSET_KIND), parameter :: displacement=0
  integer, parameter :: order=MPI_ORDER_FORTRAN, io_info=MPI_INFO_NULL
!  integer, dimension(io_dims) :: local_size, local_start, global_size, global_start, subsize
  integer, dimension (:), allocatable :: local_size, local_start, global_size, global_start, subsize

!
  contains
!***********************************************************************
    subroutine register_io
!
!  dummy routine, generates separate directory for each processor.
!  VAR#-files are written to the directory directory_snap which will
!  be the same as directory, unless specified otherwise.
!
!  06-Oct-2015/PABourdin: reworked, should work now
!
      integer :: alloc_err
!
      if (lroot) call svn_id ("$Id$")
      if (ldistribute_persist) call fatal_error ('io_mpi2', "Distibuted persistent variables are not allowed with MPI-IO.")
!
      lmonolithic_io = .true.
!
      if (lwrite_2d) then
        io_dims=3
      else
        io_dims=4
      endif
      allocate (local_size(io_dims), local_start(io_dims), global_size(io_dims), global_start(io_dims), subsize(io_dims), &
          stat=alloc_err)
!
      if (.not. lwrite_2D) then
!
! Define local full data size.
!
        local_size(1) = mx
        local_size(2) = my
        local_size(3) = mz
        local_size(4) = -1
!
! Define the subsize of the local data portion to be saved in the global file.
!
        subsize(1) = nx
        subsize(2) = ny
        subsize(3) = nz
        subsize(4) = -1
!
! We need to save also the outermost ghost layers if we are on either boundary.
! This data subsize is identical for the local portion and the global file.
!
        if (lfirst_proc_x) subsize(1) = subsize(1) + nghost
        if (lfirst_proc_y) subsize(2) = subsize(2) + nghost
        if (lfirst_proc_z) subsize(3) = subsize(3) + nghost
        if (llast_proc_x) subsize(1) = subsize(1) + nghost
        if (llast_proc_y) subsize(2) = subsize(2) + nghost
        if (llast_proc_z) subsize(3) = subsize(3) + nghost
!
! The displacements in 'local_start' use C-like format, ie. start from zero.
!
        local_start(1) = nghost
        local_start(2) = nghost
        local_start(3) = nghost
        local_start(4) = 0
!
! We need to include lower ghost cells if we are on a lower edge;
! inclusion of upper ghost cells is taken care of by increased subsize.
!
        if (lfirst_proc_x) local_start(1) = local_start(1) - nghost
        if (lfirst_proc_y) local_start(2) = local_start(2) - nghost
        if (lfirst_proc_z) local_start(3) = local_start(3) - nghost
!
! Define the size of this processor's data in the global file.
!
        global_size(1) = mxgrid
        global_size(2) = mygrid
        global_size(3) = mzgrid
        global_size(4) = 0
!
! Define starting position of this processor's data portion in global file.
!
        global_start(1) = ipx*nx + local_start(1)
        global_start(2) = ipy*ny + local_start(2)
        global_start(3) = ipz*nz + local_start(3)
        global_start(4) = 0
      else
!
! 2D version for ny=1 or nz=1 when setting lwrite_2D=True
!
        local_size(1) = mx
        local_size(3) = -1
        subsize(1) = nx
        subsize(3) = -1
        if (lfirst_proc_x) subsize(1) = subsize(1) + nghost
        if (llast_proc_x) subsize(1) = subsize(1) + nghost
        local_start(1) = nghost
        local_start(2) = nghost
        local_start(3) = 0
        if (lfirst_proc_x) local_start(1) = local_start(1) - nghost
        global_size(1) = mxgrid
        global_size(3) = 0
        global_start(1) = ipx*nx + local_start(1)
        global_start(3) = 0
        if (ny == 1) then
          local_size(2) = mz
          subsize(2) = nz
          if (lfirst_proc_z) subsize(2) = subsize(2) + nghost
          if (llast_proc_z) subsize(2) = subsize(2) + nghost
          if (lfirst_proc_z) local_start(2) = local_start(2) - nghost
          global_size(2) = mzgrid
          global_start(2) = ipz*nz + local_start(2)
        else
          local_size(2) = my
          subsize(2) = ny
          if (lfirst_proc_y) subsize(2) = subsize(2) + nghost
          if (llast_proc_y) subsize(2) = subsize(2) + nghost
          if (lfirst_proc_y) local_start(2) = local_start(2) - nghost
          global_size(2) = mygrid
          global_start(2) = ipy*ny + local_start(2)
        endif
      endif
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
!  check whether directory_snap contains `/allprocs' -- if so, revert to the
!  default name.
!  Rationale: if directory_snap was not explicitly set in start.in, it
!  will be written to param.nml as 'data/allprocs'.
!
      if ((datadir_snap == '') .or. (index(datadir_snap,'allprocs')>0)) &
        datadir_snap = datadir

      call directory_names_std
!
    endsubroutine directory_names
!***********************************************************************
    subroutine distribute_grid(x, y, z, gx, gy, gz)
!
!  This routine distributes the global grid to all processors.
!
!  11-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      real, dimension(mx), intent(out) :: x
      real, dimension(my), intent(out) :: y
      real, dimension(mz), intent(out) :: z
      real, dimension(nxgrid+2*nghost), intent(in), optional :: gx
      real, dimension(nygrid+2*nghost), intent(in), optional :: gy
      real, dimension(nzgrid+2*nghost), intent(in), optional :: gz
!
      integer :: px, py, pz, partner
      integer, parameter :: tag_gx=680, tag_gy=681, tag_gz=682
!
      if (lroot) then
        ! send local x-data to all leading yz-processors along the x-direction
        x = gx(1:mx)
        do px = 0, nprocx-1
          if (px == 0) cycle
          call mpisend_real (gx(px*nx+1:px*nx+mx), mx, px, tag_gx)
        enddo
        ! send local y-data to all leading xz-processors along the y-direction
        y = gy(1:my)
        do py = 0, nprocy-1
          if (py == 0) cycle
          call mpisend_real (gy(py*ny+1:py*ny+my), my, py*nprocx, tag_gy)
        enddo
        ! send local z-data to all leading xy-processors along the z-direction
        z = gz(1:mz)
        do pz = 0, nprocz-1
          if (pz == 0) cycle
          call mpisend_real (gz(pz*nz+1:pz*nz+mz), mz, pz*nprocxy, tag_gz)
        enddo
      endif
      if (lfirst_proc_yz) then
        ! receive local x-data from root processor
        if (.not. lroot) call mpirecv_real (x, mx, 0, tag_gx)
        ! send local x-data to all other processors in the same yz-plane
        do py = 0, nprocy-1
          do pz = 0, nprocz-1
            partner = ipx + py*nprocx + pz*nprocxy
            if (partner == iproc) cycle
            call mpisend_real (x, mx, partner, tag_gx)
          enddo
        enddo
      else
        ! receive local x-data from leading yz-processor
        call mpirecv_real (x, mx, ipx, tag_gx)
      endif
      if (lfirst_proc_xz) then
        ! receive local y-data from root processor
        if (.not. lroot) call mpirecv_real (y, my, 0, tag_gy)
        ! send local y-data to all other processors in the same xz-plane
        do px = 0, nprocx-1
          do pz = 0, nprocz-1
            partner = px + ipy*nprocx + pz*nprocxy
            if (partner == iproc) cycle
            call mpisend_real (y, my, partner, tag_gy)
          enddo
        enddo
      else
        ! receive local y-data from leading xz-processor
        call mpirecv_real (y, my, ipy*nprocx, tag_gy)
      endif
      if (lfirst_proc_xy) then
        ! receive local z-data from root processor
        if (.not. lroot) call mpirecv_real (z, mz, 0, tag_gz)
        ! send local z-data to all other processors in the same xy-plane
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (partner == iproc) cycle
            call mpisend_real (z, mz, partner, tag_gz)
          enddo
        enddo
      else
        ! receive local z-data from leading xy-processor
        call mpirecv_real (z, mz, ipz*nprocxy, tag_gz)
      endif
!
    endsubroutine distribute_grid
!***********************************************************************
    subroutine check_success(routine, message, file)
!
!  Check success of MPI2 file access and issue error if necessary.
!
!  03-Oct-2015/PABourdin: coded
!
      character (len=*), intent(in) :: routine, message, file
!
      if (mpi_err /= MPI_SUCCESS) call fatal_error (routine//'_snap', 'Could not '//message//' "'//file//'"')
!
    endsubroutine check_success
!***********************************************************************
    subroutine output_snap(a, nv, file, mode)
!
!  Write snapshot file, always write mesh and time, could add other things.
!
!  10-Feb-2012/PABourdin: coded
!  13-feb-2014/MR: made file optional (prep for downsampled output)
!
      use Mpicomm, only: globalize_xy, collect_grid, mpi_precision
!
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: a
      character (len=*), optional, intent(in) :: file
      integer, optional, intent(in) :: mode
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: handle, alloc_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (.not. present (file)) call fatal_error ('output_snap', 'downsampled output not implemented for IO_mpi2')
!
      lwrite_add = .true.
      if (present (mode)) lwrite_add = (mode == 1)
!
      local_size(io_dims) = nv
      global_size(io_dims) = nv
      subsize(io_dims) = nv
!
! Create 'local_type' to be the local data portion that is being saved.
!
      call MPI_TYPE_CREATE_SUBARRAY (io_dims, local_size, subsize, local_start, order, mpi_precision, local_type, mpi_err)
      call check_success ('output', 'create local subarray', file)
      call MPI_TYPE_COMMIT (local_type, mpi_err)
      call check_success ('output', 'commit local type', file)
!
! Create 'global_type' to indicate the local data portion in the global file.
!
      call MPI_TYPE_CREATE_SUBARRAY (io_dims, global_size, subsize, global_start, order, mpi_precision, global_type, mpi_err)
      call check_success ('output', 'create global subarray', file)
      call MPI_TYPE_COMMIT (global_type, mpi_err)
      call check_success ('output', 'commit global type', file)
!
      call MPI_FILE_OPEN (mpi_comm, trim (directory_snap)//'/'//file, MPI_MODE_CREATE+MPI_MODE_WRONLY, io_info, handle, mpi_err)
      call check_success ('output', 'open', trim (directory_snap)//'/'//file)
!
! Setting file view and write raw binary data, ie. 'native'.
!
      call MPI_FILE_SET_VIEW (handle, displacement, mpi_precision, global_type, 'native', io_info, mpi_err)
      call check_success ('output', 'create view', file)
!
      if (lwrite_2D) then
        if (ny == 1) then
          call MPI_FILE_WRITE_ALL (handle, a(:,m1,:,:), 1, local_type, status, mpi_err)
        else
          call MPI_FILE_WRITE_ALL (handle, a(:,:,n1,:), 1, local_type, status, mpi_err)
        endif
      else
        call MPI_FILE_WRITE_ALL (handle, a, 1, local_type, status, mpi_err)
      endif
      call check_success ('output', 'write', file)
!
      call MPI_FILE_CLOSE (handle, mpi_err)
      call check_success ('output', 'close', file)
!
      ! write additional data:
      if (lwrite_add) then
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('output_snap', 'Could not allocate memory for gx,gy,gz', .true.)
        endif
!
        call collect_grid (x, y, z, gx, gy, gz)
        if (lroot) then
          open (lun_output, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', position='append', status='old')
          t_sp = real (t)
          write (lun_output) t_sp, gx, gy, gz, dx, dy, dz
        endif
      endif
!
    endsubroutine output_snap
!***********************************************************************
    subroutine output_snap_finalize
!
!  Close snapshot file.
!
!  11-Feb-2012/PABourdin: coded
!
      if (persist_initialized) then
        if (lroot .and. (ip <= 9)) write (*,*) 'finish persistent block'
        if (lroot) then
          write (lun_output) id_block_PERSISTENT
          close (lun_output)
        endif
        persist_initialized = .false.
        persist_last_id = -max_int
      elseif (lwrite_add .and. lroot) then
        close (lun_output)
      endif
!
    endsubroutine output_snap_finalize
!***********************************************************************
    subroutine output_part_snap(ipar, a, mv, nv, file, label, ltruncate)
!
!  Write particle snapshot file, always write mesh and time.
!
!  23-Oct-2018/PABourdin: adapted from output_snap
!
      integer, intent(in) :: mv, nv
      integer, dimension (mv), intent(in) :: ipar
      real, dimension (mv,mparray), intent(in) :: a
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: ltruncate
!
      call fatal_error ('output_part_snap', 'not implemented for "io_mpi2"', .true.)
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
      call fatal_error ('output_stalker_init', 'not implemented for "io_mpi2"', .true.)
!
    endsubroutine output_stalker_init
!***********************************************************************
    subroutine output_stalker(label, mv, nv, data)
!
!  Write stalker particle quantity to snapshot file.
!
!  03-May-2019/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: mv, nv
      real, dimension (mv), intent(in) :: data
!
      call fatal_error ('output_stalker', 'not implemented for "io_mpi2"', .true.)
!
    endsubroutine output_stalker
!***********************************************************************
    subroutine output_part_finalize
!
!  Close particle snapshot file.
!
!  03-May-2019/PABourdin: coded
!
      call fatal_error ('output_part_finalize', 'not implemented for "io_mpi2"', .true.)
!
    endsubroutine output_part_finalize
!***********************************************************************
    subroutine output_pointmass(file, labels, fq, mv, nc)
!
!  Write pointmass snapshot file with time.
!
!  26-Oct-2018/PABourdin: adapted from output_snap
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mv, nc
      character (len=*), dimension (mqarray), intent(in) :: labels
      real, dimension (mv,mparray), intent(in) :: fq
!
      call fatal_error ('output_pointmass', 'not implemented for "io_mpi2"', .true.)
!
    endsubroutine output_pointmass
!***********************************************************************
    subroutine input_snap(file, a, nv, mode)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!  10-Feb-2012/PABourdin: coded
!  10-mar-2015/MR: avoided use of fseek
!
      use File_io, only: backskip_to_time
      use Mpicomm, only: localize_xy, mpibcast_real, mpi_precision
!
      character (len=*) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(out) :: a
      integer, optional, intent(in) :: mode
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: handle, alloc_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
      lread_add = .true.
      if (present (mode)) lread_add = (mode == 1)
!
      local_size(io_dims) = nv
      global_size(io_dims) = nv
      subsize(io_dims) = nv
!
! Create 'local_type' to be the local data portion that is being saved.
!
      call MPI_TYPE_CREATE_SUBARRAY (io_dims, local_size, subsize, local_start, order, mpi_precision, local_type, mpi_err)
      call check_success ('input', 'create local subarray', file)
      call MPI_TYPE_COMMIT (local_type, mpi_err)
      call check_success ('input', 'commit local subarray', file)
!
! Create 'global_type' to indicate the local data portion in the global file.
!
      call MPI_TYPE_CREATE_SUBARRAY (io_dims, global_size, subsize, global_start, order, mpi_precision, global_type, mpi_err)
      call check_success ('input', 'create global subarray', file)
      call MPI_TYPE_COMMIT (global_type, mpi_err)
      call check_success ('input', 'commit global subarray', file)
!
      call MPI_FILE_OPEN (mpi_comm, trim (directory_snap)//'/'//file, MPI_MODE_RDONLY, io_info, handle, mpi_err)
      call check_success ('input', 'open', trim (directory_snap)//'/'//file)
!
! Setting file view and read raw binary data, ie. 'native'.
!
      call MPI_FILE_SET_VIEW (handle, displacement, mpi_precision, global_type, 'native', io_info, mpi_err)
      call check_success ('input', 'create view', file)
!
      if (lwrite_2D) then
        if (ny == 1) then
          call MPI_FILE_READ_ALL (handle, a(:,m1,:,:), 1, local_type, status, mpi_err)
        else
          call MPI_FILE_READ_ALL (handle, a(:,:,n1,:), 1, local_type, status, mpi_err)
        endif
      else
        call MPI_FILE_READ_ALL (handle, a, 1, local_type, status, mpi_err)
      endif
      call check_success ('input', 'read', file)
!
      call MPI_FILE_CLOSE (handle, mpi_err)
      call check_success ('input', 'close', file)
!
      ! read additional data
      if (lread_add) then
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('input_snap', 'Could not allocate memory for gx,gy,gz', .true.)
!
          open (lun_input, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', position='append', status='old')
          call backskip_to_time(lun_input)
          read (lun_input) t_sp, gx, gy, gz, dx, dy, dz
          call distribute_grid (x, y, z, gx, gy, gz)
          deallocate (gx, gy, gz)
        else
          call distribute_grid (x, y, z)
        endif
        call mpibcast_real (t_sp,comm=MPI_COMM_WORLD)
        t = t_sp
      endif
!
    endsubroutine input_snap
!***********************************************************************
    subroutine input_snap_finalize
!
!  Close snapshot file.
!
!  11-Feb-2012/PABourdin: coded
!
      if (persist_initialized) then
        if (ldistribute_persist .or. lroot) close (lun_input)
        persist_initialized = .false.
        persist_last_id = -max_int
      elseif (lread_add .and. lroot) then
        close (lun_input)
      endif
!
    endsubroutine input_snap_finalize
!***********************************************************************
    subroutine input_part_snap(ipar, ap, mv, nv, npar_total, file, label)
!
!  Read particle snapshot file, mesh and time are read in 'input_snap'.
!
!  24-Oct-2018/PABourdin: apadpted and moved to IO module
!
      integer, intent(in) :: mv
      integer, dimension (mv), intent(out) :: ipar
      real, dimension (mv,mparray), intent(out) :: ap
      integer, intent(out) :: nv, npar_total
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label
!
      call fatal_error ('input_part_snap', 'not implemented for "io_mpi2"', .true.)
!
    endsubroutine input_part_snap
!***********************************************************************
    subroutine input_pointmass(file, labels, fq, mv, nc)
!
!  Read pointmass snapshot file.
!
!  26-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mv, nc
      character (len=*), dimension (nc), intent(in) :: labels
      real, dimension (mv,nc), intent(out) :: fq
!
      call fatal_error ('input_pointmass', 'not implemented for "io_mpi2"', .true.)
!
    endsubroutine input_pointmass
!***********************************************************************
    logical function init_write_persist(file)
!
!  Initialize writing of persistent data to persistent file.
!
!  13-Dec-2011/PABourdin: coded
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
      if (lroot) then
        if (filename /= "") then
          if (lroot .and. (ip <= 9)) write (*,*) 'begin write persistent block'
          close (lun_output)
          open (lun_output, FILE=trim (directory_snap)//'/'//filename, FORM='unformatted', status='replace')
          filename = ""
        endif
        write (lun_output) id_block_PERSISTENT
      endif
!
      init_write_persist = .false.
      persist_initialized = .true.
!
    endfunction init_write_persist
!***********************************************************************
    logical function write_persist_id(label, id)
!
!  Write persistent data to snapshot file.
!
!  13-Dec-2011/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
!
      write_persist_id = .true.
      if (.not. persist_initialized) write_persist_id = init_write_persist ()
      if (.not. persist_initialized) return
!
      if (persist_last_id /= id) then
        if (lroot) then
          if (ip <= 9) write (*,*) 'write persistent ID '//trim (label)
          write (lun_output) id
        endif
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
!  12-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, intent(in) :: value
!
      integer :: px, py, pz, partner, alloc_err
      integer, parameter :: tag_log_0D = 700
      logical, dimension (:,:,:), allocatable :: global
      logical :: buffer
!
      write_persist_logical_0D = .true.
      if (write_persist_id (label, id)) return
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_logical_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1,ipz+1) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpirecv_logical (buffer, partner, tag_log_0D)
              global(px+1,py+1,pz+1) = buffer
            enddo
          enddo
        enddo
        if (ip <= 9) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global)
      else
        call mpisend_logical (value, 0, tag_log_0D)
      endif
!
      write_persist_logical_0D = .false.
!
    endfunction write_persist_logical_0D
!***********************************************************************
    logical function write_persist_logical_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  12-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, dimension(:), intent(in) :: value
!
      integer :: px, py, pz, partner, nv, alloc_err
      integer, parameter :: tag_log_1D = 701
      logical, dimension (:,:,:,:), allocatable :: global
      logical, dimension (:), allocatable :: buffer
!
      write_persist_logical_1D = .true.
      if (write_persist_id (label, id)) return
!
      nv = size (value)
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz,nv), buffer(nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_logical_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1,ipz+1,:) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpirecv_logical (buffer, nv, partner, tag_log_1D)
              global(px+1,py+1,pz+1,:) = buffer
            enddo
          enddo
        enddo
        if (ip <= 9) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global, buffer)
      else
        call mpisend_logical (value, nv, 0, tag_log_1D)
      endif
!
      write_persist_logical_1D = .false.
!
    endfunction write_persist_logical_1D
!***********************************************************************
    logical function write_persist_int_0D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  12-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, intent(in) :: value
!
      integer :: px, py, pz, partner, alloc_err
      integer, parameter :: tag_int_0D = 702
      integer, dimension (:,:,:), allocatable :: global
      integer :: buffer
!
      write_persist_int_0D = .true.
      if (write_persist_id (label, id)) return
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_int_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1,ipz+1) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpirecv_int (buffer, partner, tag_int_0D)
              global(px+1,py+1,pz+1) = buffer
            enddo
          enddo
        enddo
        if (ip <= 9) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global)
      else
        call mpisend_int (value, 0, tag_int_0D)
      endif
!
      write_persist_int_0D = .false.
!
    endfunction write_persist_int_0D
!***********************************************************************
    logical function write_persist_int_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  12-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, dimension (:), intent(in) :: value
!
      integer :: px, py, pz, partner, nv, alloc_err
      integer, parameter :: tag_int_1D = 703
      integer, dimension (:,:,:,:), allocatable :: global
      integer, dimension (:), allocatable :: buffer
!
      write_persist_int_1D = .true.
      if (write_persist_id (label, id)) return
!
      nv = size (value)
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz,nv), buffer(nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_int_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1,ipz+1,:) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpirecv_int (buffer, nv, partner, tag_int_1D)
              global(px+1,py+1,pz+1,:) = buffer
            enddo
          enddo
        enddo
        if (ip <= 9) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global, buffer)
      else
        call mpisend_int (value, nv, 0, tag_int_1D)
      endif
!
      write_persist_int_1D = .false.
!
    endfunction write_persist_int_1D
!***********************************************************************
    logical function write_persist_real_0D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  12-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, intent(in) :: value
!
      integer :: px, py, pz, partner, alloc_err
      integer, parameter :: tag_real_0D = 704
      real, dimension (:,:,:), allocatable :: global
      real :: buffer
!
      write_persist_real_0D = .true.
      if (write_persist_id (label, id)) return
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_real_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1,ipz+1) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpirecv_real (buffer, partner, tag_real_0D)
              global(px+1,py+1,pz+1) = buffer
            enddo
          enddo
        enddo
        if (ip <= 9) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global)
      else
        call mpisend_real (value, 0, tag_real_0D)
      endif
!
      write_persist_real_0D = .false.
!
    endfunction write_persist_real_0D
!***********************************************************************
    logical function write_persist_real_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  12-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, dimension (:), intent(in) :: value
!
      integer :: px, py, pz, partner, nv, alloc_err
      integer, parameter :: tag_real_1D = 705
      real, dimension (:,:,:,:), allocatable :: global
      real, dimension (:), allocatable :: buffer
!
      write_persist_real_1D = .true.
      if (write_persist_id (label, id)) return
!
      nv = size (value)
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz,nv), buffer(nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_real_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1,ipz+1,:) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpirecv_real (buffer, nv, partner, tag_real_1D)
              global(px+1,py+1,pz+1,:) = buffer
            enddo
          enddo
        enddo
        if (ip <= 9) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global, buffer)
      else
        call mpisend_real (value, nv, 0, tag_real_1D)
      endif
!
      write_persist_real_1D = .false.
!
    endfunction write_persist_real_1D
!***********************************************************************
    logical function init_read_persist(file)
!
!  Initialize reading of persistent data from persistent file.
!
!  13-Dec-2011/PABourdin: coded
!
      use File_io, only: file_exists
      use Mpicomm, only: mpibcast_logical
!
      character (len=*), intent(in), optional :: file
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
        if (lroot .and. (ip <= 9)) write (*,*) 'begin read persistent block'
        open (lun_input, FILE=trim (directory_dist)//'/'//file, FORM='unformatted', status='old')
      endif
!
      init_read_persist = .false.
      persist_initialized = .true.
!
    endfunction init_read_persist
!***********************************************************************
    logical function read_persist_id(label, id, lerror_prone)
!
!  Read persistent block ID from snapshot file.
!
!  17-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpibcast_int
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
      if (lroot) then
        if (ip <= 9) write (*,*) 'read persistent ID '//trim (label)
        if (lcatch_error) then
          read (lun_input, iostat=io_err) id
          if (io_err /= 0) id = -max_int
        else
          read (lun_input) id
        endif
      endif
!
      call mpibcast_int (id,comm=MPI_COMM_WORLD)
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
!  11-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      logical, intent(out) :: value
!
      integer :: px, py, pz, partner, alloc_err
      integer, parameter :: tag_log_0D = 706
      logical, dimension (:,:,:), allocatable :: global
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_logical_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (ip <= 9) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1,ipz+1)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpisend_logical (global(px+1,py+1,pz+1), partner, tag_log_0D)
            enddo
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_logical (value, 0, tag_log_0D)
      endif
!
      read_persist_logical_0D = .false.
!
    endfunction read_persist_logical_0D
!***********************************************************************
    logical function read_persist_logical_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      logical, dimension(:), intent(out) :: value
!
      integer :: px, py, pz, partner, nv, alloc_err
      integer, parameter :: tag_log_1D = 707
      logical, dimension (:,:,:,:), allocatable :: global
!
      nv = size (value)
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_logical_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (ip <= 9) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1,ipz+1,:)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpisend_logical (global(px+1,py+1,pz+1,:), nv, partner, tag_log_1D)
            enddo
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_logical (value, nv, 0, tag_log_1D)
      endif
!
      read_persist_logical_1D = .false.
!
    endfunction read_persist_logical_1D
!***********************************************************************
    logical function read_persist_int_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, intent(out) :: value
!
      integer :: px, py, pz, partner, alloc_err
      integer, parameter :: tag_int_0D = 708
      integer, dimension (:,:,:), allocatable :: global
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_int_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (ip <= 9) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1,ipz+1)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpisend_int (global(px+1,py+1,pz+1), partner, tag_int_0D)
            enddo
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_int (value, 0, tag_int_0D)
      endif
!
      read_persist_int_0D = .false.
!
    endfunction read_persist_int_0D
!***********************************************************************
    logical function read_persist_int_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, dimension(:), intent(out) :: value
!
      integer :: px, py, pz, partner, nv, alloc_err
      integer, parameter :: tag_int_1D = 709
      integer, dimension (:,:,:,:), allocatable :: global
!
      nv = size (value)
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_int_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (ip <= 9) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1,ipz+1,:)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpisend_int (global(px+1,py+1,pz+1,:), nv, partner, tag_int_1D)
            enddo
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_int (value, nv, 0, tag_int_1D)
      endif
!
      read_persist_int_1D = .false.
!
    endfunction read_persist_int_1D
!***********************************************************************
    logical function read_persist_real_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      real, intent(out) :: value
!
      integer :: px, py, pz, partner, alloc_err
      integer, parameter :: tag_real_0D = 710
      real, dimension (:,:,:), allocatable :: global
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_real_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (ip <= 9) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1,ipz+1)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpisend_real (global(px+1,py+1,pz+1), partner, tag_real_0D)
            enddo
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_real (value, 0, tag_real_0D)
      endif
!
      read_persist_real_0D = .false.
!
    endfunction read_persist_real_0D
!***********************************************************************
    logical function read_persist_real_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/PABourdin: coded
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      real, dimension(:), intent(out) :: value
!
      integer :: px, py, pz, partner, nv, alloc_err
      integer, parameter :: tag_real_1D = 711
      real, dimension (:,:,:,:), allocatable :: global
!
      nv = size (value)
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_real_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (ip <= 9) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1,ipz+1,:)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpisend_real (global(px+1,py+1,pz+1,:), nv, partner, tag_real_1D)
            enddo
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_real (value, nv, 0, tag_real_1D)
      endif
!
      read_persist_real_1D = .false.
!
    endfunction read_persist_real_1D
!***********************************************************************
    subroutine output_globals(file, a, nv, label)
!
!  Write snapshot file of globals, ignore time and mesh.
!
!  10-Feb-2012/PABourdin: coded
!
      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*), intent(in), optional :: label
!
      call output_snap (a, nv, file, 0)
      call output_snap_finalize
!
    endsubroutine output_globals
!***********************************************************************
    subroutine input_globals(file, a, nv)
!
!  Read globals snapshot file, ignore time and mesh.
!
!  10-Feb-2012/PABourdin: coded
!
      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
!
      call input_snap (file, a, nv, 0)
      call input_snap_finalize
!
    endsubroutine input_globals
!***********************************************************************
    subroutine log_filename_to_file(filename, flist)
!
!  In the directory containing 'filename', append one line to file
!  'flist' containing the file part of filename
!
      use General, only: parse_filename, safe_character_assign
      use Mpicomm, only: mpibarrier
!
      character (len=*) :: filename, flist
!
      character (len=fnlen) :: dir, fpart
!
      call parse_filename (filename, dir, fpart)
      if (dir == '.') call safe_character_assign (dir, directory_snap)
!
      if (lroot) then
        open (lun_output, FILE=trim (dir)//'/'//trim (flist), POSITION='append')
        write (lun_output, '(A)') trim (fpart)
        close (lun_output)
      endif
!
      if (lcopysnapshots_exp) then
        call mpibarrier
        if (lroot) then
          open (lun_output,FILE=trim (datadir)//'/move-me.list', POSITION='append')
          write (lun_output,'(A)') trim (fpart)
          close (lun_output)
        endif
      endif
!
    endsubroutine log_filename_to_file
!***********************************************************************
    subroutine wgrid(file,mxout,myout,mzout)
!
!  Write grid coordinates.
!
!  10-Feb-2012/PABourdin: adapted for collective IO
!
      use Mpicomm, only: collect_grid
!
      character (len=*) :: file
      integer, optional :: mxout,myout,mzout
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (lyang) return      ! grid collection only needed on Yin grid, as grids are identical

      if (lroot) then
        allocate (gx(nxgrid+2*nghost), gy(nygrid+2*nghost), gz(nzgrid+2*nghost), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('wgrid', 'Could not allocate memory for gx,gy,gz', .true.)
!
        open (lun_output, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', status='replace')
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
!
    endsubroutine wgrid
!***********************************************************************
    subroutine rgrid(file)
!
!  Read grid coordinates.
!
!  21-jan-02/wolf: coded
!  15-jun-03/axel: Lx,Ly,Lz are now read in from file (Tony noticed the mistake)
!  10-Feb-2012/PABourdin: adapted for collective IO
!
      use Mpicomm, only: mpibcast_int, mpibcast_real
!
      character (len=*) :: file
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (lroot) then
        allocate (gx(nxgrid+2*nghost), gy(nygrid+2*nghost), gz(nzgrid+2*nghost), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('rgrid', 'Could not allocate memory for gx,gy,gz', .true.)
!
        open (lun_input, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', status='old')
        read (lun_input) t_sp, gx, gy, gz, dx, dy, dz
        call distribute_grid (x, y, z, gx, gy, gz)
        read (lun_input) dx, dy, dz
        read (lun_input) Lx, Ly, Lz
        read (lun_input) gx, gy, gz
        call distribute_grid (dx_1, dy_1, dz_1, gx, gy, gz)
        read (lun_input) gx, gy, gz
        call distribute_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
        close (lun_input)
!
        deallocate (gx, gy, gz)
      else
        call distribute_grid (x, y, z)
        call distribute_grid (dx_1, dy_1, dz_1)
        call distribute_grid (dx_tilde, dy_tilde, dz_tilde)
      endif
!
      call mpibcast_real (dx,comm=MPI_COMM_WORLD)
      call mpibcast_real (dy,comm=MPI_COMM_WORLD)
      call mpibcast_real (dz,comm=MPI_COMM_WORLD)
      call mpibcast_real (Lx,comm=MPI_COMM_WORLD)
      call mpibcast_real (Ly,comm=MPI_COMM_WORLD)
      call mpibcast_real (Lz,comm=MPI_COMM_WORLD)
!
      if (lroot.and.ip <= 4) then
        print *, 'rgrid: Lx,Ly,Lz=', Lx, Ly, Lz
        print *, 'rgrid: dx,dy,dz=', dx, dy, dz
      endif
!
    endsubroutine rgrid
!***********************************************************************
    subroutine wproc_bounds(file)
!
!   Export processor boundaries to file.
!
!   22-Feb-2012/PABourdin: adapted from io_dist
!
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
!
      integer :: ierr
!
      open (lun_output, FILE=file, FORM='unformatted', IOSTAT=ierr, status='replace')
      if (ierr /= 0) call stop_it ( &
          'Cannot open "' // trim(file) // '" (or similar) for reading' // &
          ' -- is "data/" visible from all nodes?')
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
!   22-Feb-2012/PABourdin: adapted from io_dist
!
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
!
      integer :: ierr
!
      open (lun_input, FILE=file, FORM='unformatted', IOSTAT=ierr, status='old')
      if (ierr /= 0) call stop_it ( &
          'Cannot open "' // trim(file) // '" (or similar) for reading' // &
          ' -- is "data/" visible from all nodes?')
      read (lun_input) procy_bounds
      read (lun_input) procz_bounds
      close (lun_input)
!
    endsubroutine rproc_bounds
!***********************************************************************
endmodule Io
