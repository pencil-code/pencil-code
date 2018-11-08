! $Id$
!
!  I/O via the HDF5 hyperslab-by-chunk IO routines.
!  (storing data into one file, e.g. data/allprocs/VAR#.h5)
!
!  The data format is self-contained. Only outer ghost-layers are stored.
!
!  19-Sep-2012/PABourdin: adapted from io_mpi2.f90
!  28-Oct-2016/PABourdin: first fully working version
!
module Io
!
  use Cdata
  use Cparam, only: intlen, fnlen, max_int
  use File_io, only: delete_file
  use HDF5_IO
  use Messages, only: fatal_error, svn_id, warning
!
  implicit none
!
  include 'hdf5_io.h'
  include 'io.h'
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
  integer :: lun_input = 88
  integer :: lun_output = 91
!
  ! Indicates if IO is done distributed (each proc writes into a procdir)
  ! or collectively (eg. by specialized IO-nodes or by MPI-IO).
  logical :: lcollective_IO = .true.
  character (len=labellen) :: IO_strategy = "HDF5"
!
  character (len=fnlen) :: last_snapshot = ""
  logical :: persist_initialized = .false.
!
  contains
!***********************************************************************
    subroutine register_io
!
!  dummy routine, generates separate directory for each processor.
!  VAR#-files are written to the directory directory_snap which will
!  be the same as directory, unless specified otherwise.
!
!  04-Jul-2011/Boudin.KIS: coded
!
      if (lroot) call svn_id ("$Id$")
!
      if (lread_from_other_prec) &
        call warning('register_io','Reading from other precision not implemented')
!
      lmonolithic_io = .true.
!
    endsubroutine register_io
!***********************************************************************
    subroutine finalize_io
!
      ! close the HDF5 library
      call finalize_hdf5
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
!
      call directory_names_std

    endsubroutine directory_names
!***********************************************************************
    subroutine distribute_grid(x, y, z, gx, gy, gz)
!
!  This routine distributes the global grid to all processors.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
    subroutine output_snap(a, nv, file, mode, ltruncate, label)
!
!  Write snapshot file, always write mesh and time, could add other things.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!  13-feb-2014/MR: made file optional (prep for downsampled output)
!
      use File_io, only: parallel_file_exists
      use Mpicomm, only: collect_grid, mpi_precision
!
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: a
      character (len=*), optional, intent(in) :: file
      integer, optional, intent(in) :: mode
      logical, optional, intent(in) :: ltruncate
      character (len=*), optional, intent(in) :: label
!
      logical :: ltrunc, lexists, lwrite_add
      character (len=fnlen) :: filename, dataset
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
!
      if (.not. present (file)) call fatal_error ('output_snap', 'downsampled output not implemented for IO_hdf5')
      filename = trim(directory_snap)//'/'//trim(file)//'.h5'
      dataset = 'f'
      if (present (label)) dataset = label
      lexists = parallel_file_exists(filename)
      ltrunc = .true.
      if (present (ltruncate)) ltrunc = ltruncate
      if (.not. lexists) ltrunc = .true.
!
      lwrite_add = .true.
      if (present (mode)) lwrite_add = (mode == 1)
!
      ! open global HDF5 file and write main data
      call file_open_hdf5 (filename, truncate=ltrunc)
      call output_hdf5 (dataset, a, nv)
!
      ! write additional data:
      if (lwrite_add) then
        call file_close_hdf5
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
        else
          allocate (gx(1), gy(1), gz(1), stat=alloc_err)
        endif
        if (alloc_err > 0) call fatal_error ('output_snap', 'allocate memory for gx,gy,gz', .true.)
!
        call collect_grid (x, y, z, gx, gy, gz)
        call file_open_hdf5 (filename, truncate=.false., global=.false.)
        call output_hdf5 ('t', t)
        call create_group_hdf5 ('grid')
        call output_hdf5 ('grid/x', gx, mxgrid)
        call output_hdf5 ('grid/y', gy, mygrid)
        call output_hdf5 ('grid/z', gz, mzgrid)
        call output_hdf5 ('grid/dx', dx)
        call output_hdf5 ('grid/dy', dy)
        call output_hdf5 ('grid/dz', dz)
        call output_hdf5 ('grid/Lx', Lx)
        call output_hdf5 ('grid/Ly', Ly)
        call output_hdf5 ('grid/Lz', Lz)
        call collect_grid (dx_1, dy_1, dz_1, gx, gy, gz)
        call output_hdf5 ('grid/dx_1', gx, mxgrid)
        call output_hdf5 ('grid/dy_1', gy, mygrid)
        call output_hdf5 ('grid/dz_1', gz, mzgrid)
        call collect_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
        call output_hdf5 ('grid/dx_tilde', gx, mxgrid)
        call output_hdf5 ('grid/dy_tilde', gy, mygrid)
        call output_hdf5 ('grid/dz_tilde', gz, mzgrid)
        call create_group_hdf5 ('dim')
        call output_hdf5 ('dim/mx', nxgrid+2*nghost)
        call output_hdf5 ('dim/my', nygrid+2*nghost)
        call output_hdf5 ('dim/mz', nzgrid+2*nghost)
        call output_hdf5 ('dim/nx', nxgrid)
        call output_hdf5 ('dim/ny', nygrid)
        call output_hdf5 ('dim/nz', nzgrid)
        call output_hdf5 ('dim/l1', nghost)
        call output_hdf5 ('dim/m1', nghost)
        call output_hdf5 ('dim/n1', nghost)
        call output_hdf5 ('dim/l2', nghost+nxgrid-1)
        call output_hdf5 ('dim/m2', nghost+nygrid-1)
        call output_hdf5 ('dim/n2', nghost+nzgrid-1)
        call output_hdf5 ('dim/nghost', nghost)
        call output_hdf5 ('dim/mvar', mvar)
        call output_hdf5 ('dim/maux', maux)
        call output_hdf5 ('dim/mglobal', mglobal)
        call output_hdf5 ('dim/nprocx', nprocx)
        call output_hdf5 ('dim/nprocy', nprocy)
        call output_hdf5 ('dim/nprocz', nprocz)
        if (mpi_precision == MPI_REAL) then
          call output_hdf5 ('dim/precision', 'S')
        else
          call output_hdf5 ('dim/precision', 'D')
        endif
        call file_close_hdf5
        deallocate (gx, gy, gz)
        call file_open_hdf5 (filename, truncate=.false.)
      endif
!
    endsubroutine output_snap
!***********************************************************************
    subroutine output_snap_finalize
!
!  Close snapshot file.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      call file_close_hdf5
      if (persist_initialized) persist_initialized = .false.
!
    endsubroutine output_snap_finalize
!***********************************************************************
    subroutine output_part_snap(ipar, ap, mv, nv, file, label, ltruncate)
!
!  Write particle snapshot file, always write mesh and time.
!
!  22-Oct-2018/PABourdin: adapted from output_snap
!
      use File_io, only: parallel_file_exists
      use Mpicomm, only: collect_grid, mpi_precision
!
      integer, intent(in) :: mv, nv
      integer, dimension (mv), intent(in) :: ipar
      real, dimension (mv,mparray), intent(in) :: ap
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label
      logical, optional, intent(in) :: ltruncate
!
      logical :: ltrunc, lexists
      character (len=fnlen) :: filename, dataset
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
!
      dataset = 'fp'
      if (present (label)) dataset = label
      filename = trim(directory_snap)//'/'//trim(file)//'.h5'
      lexists = parallel_file_exists(filename)
      ltrunc = .true.
      if (present (ltruncate)) ltrunc = ltruncate
      if (.not. lexists) ltrunc = .true.
!
      ! open global HDF5 file and write particle data
      call file_open_hdf5 (filename, truncate=ltrunc)
      call output_hdf5 (dataset, ap, mv, mparray, nv)
      call output_hdf5 ('part/ID', ipar(1:nv), nv)
      call file_close_hdf5
!
      ! write additional data:
      if (ltrunc .and. (trim (dataset) == 'fp')) then
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
        else
          allocate (gx(1), gy(1), gz(1), stat=alloc_err)
        endif
        if (alloc_err > 0) call fatal_error ('output_part_snap', 'allocate memory for gx,gy,gz', .true.)
!
        call collect_grid (x, y, z, gx, gy, gz)
        call file_open_hdf5 (filename, truncate=.false., global=.false.)
        call output_hdf5 ('t', t)
        call create_group_hdf5 ('grid')
        call output_hdf5 ('grid/x', gx, mxgrid)
        call output_hdf5 ('grid/y', gy, mygrid)
        call output_hdf5 ('grid/z', gz, mzgrid)
        call output_hdf5 ('grid/dx', dx)
        call output_hdf5 ('grid/dy', dy)
        call output_hdf5 ('grid/dz', dz)
        call output_hdf5 ('grid/Lx', Lx)
        call output_hdf5 ('grid/Ly', Ly)
        call output_hdf5 ('grid/Lz', Lz)
        call collect_grid (dx_1, dy_1, dz_1, gx, gy, gz)
        call output_hdf5 ('grid/dx_1', gx, mxgrid)
        call output_hdf5 ('grid/dy_1', gy, mygrid)
        call output_hdf5 ('grid/dz_1', gz, mzgrid)
        call collect_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
        call output_hdf5 ('grid/dx_tilde', gx, mxgrid)
        call output_hdf5 ('grid/dy_tilde', gy, mygrid)
        call output_hdf5 ('grid/dz_tilde', gz, mzgrid)
        call create_group_hdf5 ('dim')
        call output_hdf5 ('dim/mx', nxgrid+2*nghost)
        call output_hdf5 ('dim/my', nygrid+2*nghost)
        call output_hdf5 ('dim/mz', nzgrid+2*nghost)
        call output_hdf5 ('dim/nx', nxgrid)
        call output_hdf5 ('dim/ny', nygrid)
        call output_hdf5 ('dim/nz', nzgrid)
        call output_hdf5 ('dim/l1', nghost)
        call output_hdf5 ('dim/m1', nghost)
        call output_hdf5 ('dim/n1', nghost)
        call output_hdf5 ('dim/l2', nghost+nxgrid-1)
        call output_hdf5 ('dim/m2', nghost+nygrid-1)
        call output_hdf5 ('dim/n2', nghost+nzgrid-1)
        call output_hdf5 ('dim/nghost', nghost)
        call output_hdf5 ('dim/mvar', mvar)
        call output_hdf5 ('dim/maux', maux)
        call output_hdf5 ('dim/mglobal', mglobal)
        call output_hdf5 ('dim/nprocx', nprocx)
        call output_hdf5 ('dim/nprocy', nprocy)
        call output_hdf5 ('dim/nprocz', nprocz)
        if (mpi_precision == MPI_REAL) then
          call output_hdf5 ('dim/precision', 'S')
        else
          call output_hdf5 ('dim/precision', 'D')
        endif
        call file_close_hdf5
        deallocate (gx, gy, gz)
      endif
!
      ! write processor boundaries
      call file_open_hdf5 (filename, truncate=.false., global=.false.)
      call output_hdf5 ('proc/bounds_x', procx_bounds(0:nprocx), nprocx+1)
      call output_hdf5 ('proc/bounds_y', procy_bounds(0:nprocy), nprocy+1)
      call output_hdf5 ('proc/bounds_z', procz_bounds(0:nprocz), nprocz+1)
      call file_close_hdf5
!
    endsubroutine output_part_snap
!***********************************************************************
    subroutine output_pointmass(file, labels, fq, mv, nc)
!
!  Write pointmass snapshot file with time.
!
!  26-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mv, nc
      character (len=*), dimension (mqarray), intent(in) :: labels
      real, dimension (mv,mparray), intent(in) :: fq
!
      integer :: pos
      character (len=fnlen) :: filename
!
      if (.not. lroot) return
!
      filename = trim (directory_snap)//'/'//trim(file)//'.h5'
      call file_open_hdf5 (filename, global=.false.)
      call output_hdf5 ('number', mv)
      if (mv > 0) then
        call create_group_hdf5 ('points')
        do pos=1, nc
          call output_hdf5 ('points/'//trim(labels(pos)), fq(:,pos), mv)
        enddo
      endif
      call output_hdf5 ('t', t)
      call file_close_hdf5
!
    endsubroutine output_pointmass
!***********************************************************************
    subroutine initialize_slice(label, pos)
!
!  27-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: pos
!
!!      if (pos >= 0) call create_group_hdf5 (label)
!
    endsubroutine initialize_slice
!***********************************************************************
    subroutine output_slice_position
!
!  27-Oct-2018/PABourdin: no action required for HDF5 output
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
!  30-Oct-2018/PABourdin: coded
!
      use File_io, only: parallel_file_exists
      use General, only: itoa
      use Mpicomm, only: mpibcast_int, mpiallreduce_min, MPI_COMM_WORLD
!
      logical, intent(in) :: lwrite
      real, intent(in) :: time
      character (len=*), intent(in) :: label, suffix
      real, dimension (:) :: grid
      integer, intent(in) :: pos, grid_pos
      integer, intent(in) :: ndim1, ndim2
      real, dimension (:,:), pointer :: data
!
      character (len=fnlen) :: filename, group
      integer :: last, this_proc, slice_proc
      real :: time_last, slice_pos
      logical :: lexists, lhas_data
!
      if (grid_pos < 0) return
!
      last = 0
      filename = trim (directory_snap)//'/slices/'//trim(label)//'_'//trim(suffix)//'.h5'
      lexists = parallel_file_exists (filename)
      if (lroot .and. lexists) then
        call file_open_hdf5 (filename, global=.false., read_only=.true.)
        if (exists_in_hdf5 ('last')) call input_hdf5 ('last', last)
        do while (last >= 1)
          call input_hdf5 (trim(itoa(last))//'/time', time_last)
          if (time > time_last) exit
          last = last - 1
        enddo
        call file_close_hdf5
      endif
      last = last + 1
      call mpibcast_int (last)
      group = trim(itoa(last))//'/'
!
      this_proc = iproc
      if (.not. lwrite) this_proc = ncpus + 1
      call mpiallreduce_min (this_proc, slice_proc, MPI_COMM_WORLD)
      if (slice_proc == iproc) then
        call file_open_hdf5 (filename, global=.false., truncate=(.not. lexists))
        call output_hdf5 ('last', last)
        call create_group_hdf5 (group)
        call output_hdf5 (trim(group)//'time', time)
        call output_hdf5 (trim(group)//'position', grid(pos))
        call output_hdf5 (trim(group)//'coordinate', grid_pos)
        call file_close_hdf5
      endif
!
      call file_open_hdf5 (filename, truncate=.false.)
      ! collect data along 'xy', 'xz', or 'yz'
      lhas_data = lwrite .and. associated(data)
      select case (suffix)
      case ('xy')
        call output_hdf5 (trim(group)//'data', data, ndim1, ndim2, nx, ny, ipx, ipy, lhas_data)
      case ('xy2')
        call output_hdf5 (trim(group)//'data', data, ndim1, ndim2, nx, ny, ipx, ipy, lhas_data)
      case ('xy3')
        call output_hdf5 (trim(group)//'data', data, ndim1, ndim2, nx, ny, ipx, ipy, lhas_data)
      case ('xy4')
        call output_hdf5 (trim(group)//'data', data, ndim1, ndim2, nx, ny, ipx, ipy, lhas_data)
      case ('xz')
        call output_hdf5 (trim(group)//'data', data, ndim1, ndim2, nx, nz, ipx, ipz, lhas_data)
      case ('xz2')
        call output_hdf5 (trim(group)//'data', data, ndim1, ndim2, nx, nz, ipx, ipz, lhas_data)
      case ('yz')
        call output_hdf5 (trim(group)//'data', data, ndim1, ndim2, ny, nz, ipy, ipz, lhas_data)
      case default
        call fatal_error ('output_hdf5', 'unknown 2D slice "'//trim (suffix)//'"', .true.)
      endselect
      call file_close_hdf5
!
    endsubroutine output_slice
!***********************************************************************
    subroutine input_snap(file, a, nv, mode, label)
!
!  Read snapshot file. Also read mesh and time, if mode==1.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!  10-Mar-2015/MR: avoided use of fseek;
!                  this subroutine seems not yet to be adapted to HDF5
!
      use File_io, only: backskip_to_time
      use Mpicomm, only: mpibcast_real, MPI_COMM_WORLD
!
      character (len=*) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(out) :: a
      integer, optional, intent(in) :: mode
      character (len=*), optional :: label
!
      character (len=fnlen) :: filename, dataset
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
      logical :: lread_add
!
      filename = trim(directory_snap)//'/'//trim(file)//'.h5'
      dataset = 'f'
      if (present (label)) dataset = label
!
      lread_add = .true.
      if (present (mode)) lread_add = (mode == 1)
!
      ! open existing HDF5 file and read data
      call file_open_hdf5 (filename, read_only=.true.)
      call input_hdf5 (dataset, a, nv)
      call file_close_hdf5
!
      ! read additional data
      if (lread_add) then
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
        else
          allocate (gx(1), gy(1), gz(1), stat=alloc_err)
        endif
        if (alloc_err > 0) call fatal_error ('input_snap', 'Could not allocate memory for gx,gy,gz', .true.)
!
        call file_open_hdf5 (filename, global=.false., read_only=.true.)
        call input_hdf5 ('t', t)
        call input_hdf5 ('grid/x', gx, mxgrid)
        call input_hdf5 ('grid/y', gy, mygrid)
        call input_hdf5 ('grid/z', gz, mzgrid)
        call distribute_grid (x, y, z, gx, gy, gz)
        call input_hdf5 ('grid/dx', dx)
        call input_hdf5 ('grid/dy', dy)
        call input_hdf5 ('grid/dz', dz)
        call input_hdf5 ('grid/Lx', Lx)
        call input_hdf5 ('grid/Ly', Ly)
        call input_hdf5 ('grid/Lz', Lz)
        call input_hdf5 ('grid/dx_1', gx, mxgrid)
        call input_hdf5 ('grid/dy_1', gy, mygrid)
        call input_hdf5 ('grid/dz_1', gz, mzgrid)
        call distribute_grid (dx_1, dy_1, dz_1, gx, gy, gz)
        call input_hdf5 ('grid/dx_tilde', gx, mxgrid)
        call input_hdf5 ('grid/dy_tilde', gy, mygrid)
        call input_hdf5 ('grid/dz_tilde', gz, mzgrid)
        call distribute_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
        call file_close_hdf5
        deallocate (gx, gy, gz)
!
        call mpibcast_real (t, comm=MPI_COMM_WORLD)
        call mpibcast_real (dx, comm=MPI_COMM_WORLD)
        call mpibcast_real (dy, comm=MPI_COMM_WORLD)
        call mpibcast_real (dz, comm=MPI_COMM_WORLD)
        call mpibcast_real (Lx, comm=MPI_COMM_WORLD)
        call mpibcast_real (Ly, comm=MPI_COMM_WORLD)
        call mpibcast_real (Lz, comm=MPI_COMM_WORLD)
      endif
!
    endsubroutine input_snap
!***********************************************************************
    subroutine input_snap_finalize
!
!  Close snapshot file.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      if (persist_initialized) then
        call file_close_hdf5
        persist_initialized = .false.
      endif
!
    endsubroutine input_snap_finalize
!***********************************************************************
    subroutine input_part_snap(ipar, ap, mv, nv, npar_total, file, label)
!
!  Read particle snapshot file, mesh and time are read in 'input_snap'.
!
!  24-Oct-2018/PABourdin: coded
!
      use Mpicomm, only: mpireduce_sum_int, mpibcast
!
      integer, intent(in) :: mv
      integer, dimension (mv), intent(out) :: ipar
      real, dimension (mv,mparray), intent(out) :: ap
      integer, intent(out) :: nv, npar_total
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label
!
      character (len=fnlen) :: filename, dataset
!
      dataset = 'fp'
      if (present (label)) dataset = label
      filename = trim(directory_snap)//'/'//trim(file)//'.h5'
!
      ! open global HDF5 file and read particle data
      call file_open_hdf5 (filename, read_only=.true.)
      call input_hdf5 (dataset, ap, mv, mparray, nv)
      call input_hdf5 ('part/ID', ipar(1:nv), nv)
      call file_close_hdf5
!
      ! Sum the total number of all particles on the root processor.
      call mpireduce_sum_int (nv, npar_total)
!
      ! read processor boundaries
      call file_open_hdf5 (filename, global=.false., read_only=.true.)
      call input_hdf5 ('proc/bounds_x', procx_bounds(0:nprocx), nprocx+1)
      call input_hdf5 ('proc/bounds_y', procy_bounds(0:nprocy), nprocy+1)
      call input_hdf5 ('proc/bounds_z', procz_bounds(0:nprocz), nprocz+1)
      call mpibcast (procx_bounds, nprocx+1)
      call mpibcast (procy_bounds, nprocy+1)
      call mpibcast (procz_bounds, nprocz+1)
      call file_close_hdf5
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
      integer :: pos, mv_in
      character (len=fnlen) :: filename
!
      filename = trim (directory_snap)//'/'//trim(file)//'.h5'
      if (lroot) then
        call file_open_hdf5 (filename, read_only=.true., global=.false.)
        call input_hdf5 ('number', mv_in)
        if (mv_in /= mv) call fatal_error("","")
        if (mv_in /= 0) then
          do pos=1, nc
            call input_hdf5 ('points/'//trim(labels(pos)), fq(:,pos), mv)
          enddo
        endif
        call file_close_hdf5
      endif
!
      call mpibcast_real(fq,(/nqpar,mqarray/))
!
    endsubroutine input_pointmass
!***********************************************************************
    logical function init_write_persist(file)
!
!  Initialize writing of persistent data to persistent file.
!
!  26-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in), optional :: file
!
      character (len=fnlen), save :: filename=""
!
      init_write_persist = .false.
!
      if (present (file)) then
        filename = trim(directory_snap)//'/'//trim(file)//'.h5'
        persist_initialized = .false.
        return
      endif
!
      if (filename /= "") then
        call file_close_hdf5
        call file_open_hdf5 (filename)
        filename = ""
      endif
!
      if (.not. exists_in_hdf5 ('persist')) call create_group_hdf5 ('persist')
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
!
      write_persist_id = .true.
      if (.not. persist_initialized) write_persist_id = init_write_persist ()
      if (.not. persist_initialized) return
      write_persist_id = .false.
!
    endfunction write_persist_id
!***********************************************************************
    logical function write_persist_logical_0D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  26-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, intent(in) :: value
!
      logical, dimension(1) :: out
!
      out(1) = value
      write_persist_logical_0D = write_persist_logical_1D (label, id, out)
!
    endfunction write_persist_logical_0D
!***********************************************************************
    logical function write_persist_logical_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  26-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, dimension(:), intent(in) :: value
!
      integer, dimension(size (value)) :: value_int
!
      write_persist_logical_1D = .true.
      if (write_persist_id (label, id)) return
!
      value_int = 0
      where (value) value_int = 1
      call output_hdf5 ('persist/'//label, value_int, size (value), same_size=.true.)
!
      write_persist_logical_1D = .false.
!
    endfunction write_persist_logical_1D
!***********************************************************************
    logical function write_persist_int_0D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  26-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, intent(in) :: value
!
      integer, dimension(1) :: out
!
      out(1) = value
      write_persist_int_0D = write_persist_int_1D (label, id, out)
!
    endfunction write_persist_int_0D
!***********************************************************************
    logical function write_persist_int_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  26-Oct-2018/PABourdin: coded
!
      use Mpicomm, only: mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, dimension (:), intent(in) :: value
!
      write_persist_int_1D = .true.
      if (write_persist_id (label, id)) return
!
      call output_hdf5 ('persist/'//label, value, size (value), same_size=.true.)
!
      write_persist_int_1D = .false.
!
    endfunction write_persist_int_1D
!***********************************************************************
    logical function write_persist_real_0D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  26-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, intent(in) :: value
!
      real, dimension(1) :: out
!
      out(1) = value
      write_persist_real_0D = write_persist_real_1D (label, id, out)
!
    endfunction write_persist_real_0D
!***********************************************************************
    logical function write_persist_real_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  26-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, dimension (:), intent(in) :: value
!
      write_persist_real_1D = .true.
      if (write_persist_id (label, id)) return
!
      call output_hdf5 ('persist/'//label, value, size (value), same_size=.true.)
!
      write_persist_real_1D = .false.
!
    endfunction write_persist_real_1D
!***********************************************************************
    logical function init_read_persist(file)
!
!  Initialize reading of persistent data from persistent file.
!
!  27-Oct-2018/PABourdin: coded
!
      use File_io, only: parallel_file_exists
!
      character (len=*), intent(in), optional :: file
!
      character (len=fnlen) :: filename
!
      init_read_persist = .true.
!
      if (present (file)) then
        filename = trim (directory_snap)//'/'//trim (file)//'.h5'
        init_read_persist = .not. parallel_file_exists (filename)
        if (init_read_persist) return
        call file_open_hdf5 (filename, read_only=.true.)
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
!  27-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(out) :: id
      logical, intent(in), optional :: lerror_prone
!
      logical :: lcatch_error, lexists
!
      lcatch_error = .false.
      if (present (lerror_prone)) lcatch_error = lerror_prone
!
      lexists = exists_in_hdf5('persist')
      if (lexists) lexists = exists_in_hdf5('persist/'//trim (label))
      if (lcatch_error .and. .not. lexists) id = -max_int
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
!  27-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      logical, intent(out) :: value
!
      logical, dimension(1) :: read
!
      read_persist_logical_0D = read_persist_logical_1D(label, read)
      value = read(1)
!
    endfunction read_persist_logical_0D
!***********************************************************************
    logical function read_persist_logical_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  27-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      logical, dimension(:), intent(out) :: value
!
      integer, dimension(size (value)) :: value_int
!
      call input_hdf5 ('persist/'//label, value_int, size (value), same_size=.true.)
      value = .false.
      where (value_int > 0) value = .true.
!
      read_persist_logical_1D = .false.
!
    endfunction read_persist_logical_1D
!***********************************************************************
    logical function read_persist_int_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  27-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(out) :: value
!
      integer, dimension(1) :: read
!
      read_persist_int_0D = read_persist_int_1D(label, read)
      value = read(1)
!
    endfunction read_persist_int_0D
!***********************************************************************
    logical function read_persist_int_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  27-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, dimension(:), intent(out) :: value
!
      call input_hdf5 ('persist/'//label, value, size (value), same_size=.true.)
!
      read_persist_int_1D = .false.
!
    endfunction read_persist_int_1D
!***********************************************************************
    logical function read_persist_real_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  27-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      real, intent(out) :: value
!
      real, dimension(1) :: read
!
      read_persist_real_0D = read_persist_real_1D(label, read)
      value = read(1)
!
    endfunction read_persist_real_0D
!***********************************************************************
    logical function read_persist_real_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  27-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      real, dimension(:), intent(out) :: value
!
      call input_hdf5 ('persist/'//label, value, size (value), same_size=.true.)
!
      read_persist_real_1D = .false.
!
    endfunction read_persist_real_1D
!***********************************************************************
    subroutine output_globals(file,a,nv)
!
!  Write snapshot file of globals, ignore time and mesh.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
!
      call output_snap (a, nv, file, mode=0, label='globals')
      call output_snap_finalize
!
    endsubroutine output_globals
!***********************************************************************
    subroutine input_globals(file,a,nv)
!
!  Read globals snapshot file, ignore time and mesh.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
!
      call input_snap (file, a, nv, mode=0, label='globals')
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
!  27-Oct-2018/PABourdin: coded
!
      use Mpicomm, only: collect_grid, mpi_precision
!
      character (len=*), intent(in) :: file
      integer, intent(in), optional :: mxout,myout,mzout
!
      character (len=fnlen) :: filename
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
!
      if (lyang) return      ! grid collection only needed on Yin grid, as grids are identical
!
      filename = trim (datadir)//'/grid.h5'
      if (lroot) then
        allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
      else
        allocate (gx(1), gy(1), gz(1), stat=alloc_err)
      endif
      if (alloc_err > 0) call fatal_error ('wgrid', 'allocate memory for gx,gy,gz', .true.)
!
      call collect_grid (x, y, z, gx, gy, gz)
      call file_open_hdf5 (filename, global=.false.)
      call create_group_hdf5 ('grid')
      call output_hdf5 ('grid/x', gx, mxgrid)
      call output_hdf5 ('grid/y', gy, mygrid)
      call output_hdf5 ('grid/z', gz, mzgrid)
      call output_hdf5 ('grid/dx', dx)
      call output_hdf5 ('grid/dy', dy)
      call output_hdf5 ('grid/dz', dz)
      call output_hdf5 ('grid/Lx', Lx)
      call output_hdf5 ('grid/Ly', Ly)
      call output_hdf5 ('grid/Lz', Lz)
      call collect_grid (dx_1, dy_1, dz_1, gx, gy, gz)
      call output_hdf5 ('grid/dx_1', gx, mxgrid)
      call output_hdf5 ('grid/dy_1', gy, mygrid)
      call output_hdf5 ('grid/dz_1', gz, mzgrid)
      call collect_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
      call output_hdf5 ('grid/dx_tilde', gx, mxgrid)
      call output_hdf5 ('grid/dy_tilde', gy, mygrid)
      call output_hdf5 ('grid/dz_tilde', gz, mzgrid)
      call create_group_hdf5 ('dim')
      call output_hdf5 ('dim/mx', nxgrid+2*nghost)
      call output_hdf5 ('dim/my', nygrid+2*nghost)
      call output_hdf5 ('dim/mz', nzgrid+2*nghost)
      call output_hdf5 ('dim/nx', nxgrid)
      call output_hdf5 ('dim/ny', nygrid)
      call output_hdf5 ('dim/nz', nzgrid)
      call output_hdf5 ('dim/l1', nghost)
      call output_hdf5 ('dim/m1', nghost)
      call output_hdf5 ('dim/n1', nghost)
      call output_hdf5 ('dim/l2', nghost+nxgrid-1)
      call output_hdf5 ('dim/m2', nghost+nygrid-1)
      call output_hdf5 ('dim/n2', nghost+nzgrid-1)
      call output_hdf5 ('dim/nghost', nghost)
      call output_hdf5 ('dim/mvar', mvar)
      call output_hdf5 ('dim/maux', maux)
      call output_hdf5 ('dim/mglobal', mglobal)
      call output_hdf5 ('dim/nprocx', nprocx)
      call output_hdf5 ('dim/nprocy', nprocy)
      call output_hdf5 ('dim/nprocz', nprocz)
      if (mpi_precision == MPI_REAL) then
        call output_hdf5 ('dim/precision', 'S')
      else
        call output_hdf5 ('dim/precision', 'D')
      endif
      call file_close_hdf5
      deallocate (gx, gy, gz)
!
    endsubroutine wgrid
!***********************************************************************
    subroutine rgrid(file)
!
!  Read grid coordinates.
!
!  27-Oct-2018/PABourdin: coded
!
      use Mpicomm, only: mpibcast_real, MPI_COMM_WORLD
!
      character (len=*) :: file
!
      character (len=fnlen) :: filename
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
!
      filename = trim (datadir)//'/grid.h5'
      if (lroot) then
        allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
      else
        allocate (gx(1), gy(1), gz(1), stat=alloc_err)
      endif
      if (alloc_err > 0) call fatal_error ('rgrid', 'Could not allocate memory for gx,gy,gz', .true.)
!
      call file_open_hdf5 (filename, global=.false., read_only=.true.)
      call input_hdf5 ('grid/x', gx, mxgrid)
      call input_hdf5 ('grid/y', gy, mygrid)
      call input_hdf5 ('grid/z', gz, mzgrid)
      call distribute_grid (x, y, z, gx, gy, gz)
      call input_hdf5 ('grid/dx', dx)
      call input_hdf5 ('grid/dy', dy)
      call input_hdf5 ('grid/dz', dz)
      call input_hdf5 ('grid/Lx', Lx)
      call input_hdf5 ('grid/Ly', Ly)
      call input_hdf5 ('grid/Lz', Lz)
      call input_hdf5 ('grid/dx_1', gx, mxgrid)
      call input_hdf5 ('grid/dy_1', gy, mygrid)
      call input_hdf5 ('grid/dz_1', gz, mzgrid)
      call distribute_grid (dx_1, dy_1, dz_1, gx, gy, gz)
      call input_hdf5 ('grid/dx_tilde', gx, mxgrid)
      call input_hdf5 ('grid/dy_tilde', gy, mygrid)
      call input_hdf5 ('grid/dz_tilde', gz, mzgrid)
      call distribute_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
      call file_close_hdf5
      deallocate (gx, gy, gz)
!
      call mpibcast_real (t,comm=MPI_COMM_WORLD)
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
    subroutine wdim_default_grid(file)
!
!  Write dimension to file.
!
!  02-Nov-2018/PABourdin: redesigned
!
      character (len=*), intent(in) :: file
!
      if (file == 'dim.dat') return
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
    subroutine output_profile(name, coord, a, type, lsave_name, lhas_ghost)
!
!  Writes a profile along any direction.
!
!  07-Nov-2018/PABourdin: coded
!
      use General, only: loptest
      use File_io, only: parallel_file_exists
!
      real, dimension(:) :: coord, a
      character (len=*) :: name
      character :: type
      logical, optional :: lsave_name, lhas_ghost
!
      character (len=fnlen) :: filename
      integer :: np, ng, ip, np_global, np1, np2
      logical :: lexists, lwrite, lp1, lp2
!
      if (.not. lwrite_prof) return
!
      np = size(coord)
      ng = 0
      if (loptest (lhas_ghost)) ng = 3
      select case (type)
      case ('x')
        np_global = (np - 2*ng) * nprocx + 2*ng
        ip = ipx
        lp1 = lfirst_proc_x
        lp2 = llast_proc_x
        lwrite = lfirst_proc_yz
      case ('y')
        np_global = (np - 2*ng) * nprocy + 2*ng
        ip = ipy
        lp1 = lfirst_proc_y
        lp2 = llast_proc_y
        lwrite = lfirst_proc_xz
      case ('z')
        np_global = (np - 2*ng) * nprocz + 2*ng
        ip = ipz
        lp1 = lfirst_proc_z
        lp2 = llast_proc_z
        lwrite = lfirst_proc_xy
      case default
        call fatal_error ('output_profile', 'unknown direction "'//type//'"')
      endselect
      np1 = 1
      np2 = np
      if (.not. lp1) np1 = np1 + ng
      if (.not. lp2) np2 = np2 - ng
!
      ! write profile
      filename = trim(datadir)//'/'//'profile_'//type//'.h5'
      lexists = parallel_file_exists (filename)
      call file_open_hdf5 (filename, truncate=(.not. lexists))
      call create_group_hdf5 (trim(name))
      call output_hdf5 (trim(name)//'/data', a, np, np_global, ip, np1, np2, ng, lwrite)
      call output_hdf5 (trim(name)//'/position', coord, np, np_global, ip, np1, np2, ng, lwrite)
      call file_close_hdf5
!
    endsubroutine output_profile
!***********************************************************************
    subroutine input_profile(name, type, a, np, lhas_ghost)
!
!  Reads a profile from a file along any direction.
!
!  07-Nov-2018/PABourdin: coded
!
      use General, only: loptest
      use Mpicomm, only: mpibcast_real, MPI_COMM_WORLD
!
      character (len=*), intent(in) :: name
      character, intent(in) :: type
      integer, intent(in) :: np
      real, dimension(np), intent(out) :: a
      logical, optional :: lhas_ghost
!
      character (len=fnlen) :: filename
      integer :: np_global, np1, np2, ng, alloc_err
      real, dimension(:), allocatable :: profile
!
      ng = 0
      if (loptest (lhas_ghost)) ng = 3
      select case (type)
      case ('x')
        np_global = (np - 2*ng) * nprocx + 2*ng
        np1 = 1 + ipx * (np - 2*ng)
        np2 = (ipx + 1) * (np - 2*ng) + 2*ng
      case ('y')
        np_global = (np - 2*ng) * nprocy + 2*ng
        np1 = 1 + ipy * (np - 2*ng)
        np2 = (ipy + 1) * (np - 2*ng) + 2*ng
      case ('z')
        np_global = (np - 2*ng) * nprocz + 2*ng
        np1 = 1 + ipy * (np - 2*ng)
        np2 = (ipy + 1) * (np - 2*ng) + 2*ng
      case default
        call fatal_error ('input_profile', 'unknown direction "'//type//'"')
      endselect
      allocate (profile(np_global), stat=alloc_err)
      if (alloc_err > 0) call fatal_error ('input_profile', 'allocate memory for profile', .true.)
!
      ! read profile
      filename = trim(datadir)//'/'//'profile_'//type//'.h5'
      call file_open_hdf5 (filename, read_only=.true., global=.false.)
      call input_hdf5 (trim(name)//'/data', profile, np_global, lhas_ghost)
      call file_close_hdf5
!
      call mpibcast_real (profile, np_global, comm=MPI_COMM_WORLD)
      a = profile(np1:np2)
      deallocate (profile)
!
!  Should we check that coord == z for type == 'z' etc.?
!
    endsubroutine input_profile
!***********************************************************************
    subroutine wproc_bounds(file)
!
!   Export processor boundaries to file.
!
!   27-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: file
!
      ! already written in particle snapshot
!
    endsubroutine wproc_bounds
!***********************************************************************
    subroutine rproc_bounds(file)
!
!   Import processor boundaries from file.
!
!   27-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: file
!
      ! already read with particle snapshot
!
    endsubroutine rproc_bounds
!***********************************************************************
endmodule Io
