! $Id$
!
!  I/O via the HDF5 hyperslab-by-chunk IO routines.
!  (storing data into one file, e.g. data/allprocs/VAR#.h5)
!
!  The data format is self-contained. Only outer ghost-layers are stored.
!
!  19-Sep-2012/PABourdin: adapted from io_mpi2.f90
!  28-Oct-2016/PABourdin: first fully working version
!  28-Nov-2018/PABourdin: first beta-test version
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
  logical :: lread_add
  character (len=fnlen) :: varfile_name
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
!  28-Oct-2016/PABourdin: redesigned
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
    subroutine output_snap(a, nv1, nv2, file, mode, ltruncate, label)
!
!  Write snapshot file, always write mesh and time, could add other things.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!  13-feb-2014/MR: made file optional (prep for downsampled output)
!  28-Oct-2016/PABourdin: redesigned
!
      use General, only: ioptest
      use File_io, only: parallel_file_exists
!
      integer, optional, intent(in) :: nv1,nv2
      real, dimension (:,:,:,:), intent(in) :: a
      character (len=*), optional, intent(in) :: file
      integer, optional, intent(in) :: mode
      logical, optional, intent(in) :: ltruncate
      character (len=*), optional, intent(in) :: label
!
      integer :: pos,na,ne
      logical :: ltrunc, lexists, lwrite_add
      character (len=fnlen) :: filename, dataset
!
      if (.not. present (file)) call fatal_error ('output_snap', 'downsampled output not implemented for IO_hdf5')
      dataset = 'f'
      if (present (label)) dataset = label
      if (dataset == 'globals') then
        filename = trim(datadir_snap)//'/'//trim(file)//'.h5'
      elseif (dataset == 'timeavg') then
        filename = trim(datadir)//'/averages/'//trim(file)//'.h5'
      else
        filename = trim(directory_snap)//'/'//trim(file)//'.h5'
      endif
      lexists = parallel_file_exists(filename)
      ltrunc = .true.
      if (present (ltruncate)) ltrunc = ltruncate
      if (.not. lexists) ltrunc = .true.
!
      lwrite_add = .true.
      if (present (mode)) lwrite_add = (mode == 1)
!
      na=ioptest(nv1,1)
      ne=ioptest(nv2,mvar_io)
!
      ! open global HDF5 file and write main data
      call file_open_hdf5 (filename, truncate=ltrunc)
      if ((dataset == 'f') .or. (dataset == 'timeavg')) then
        call create_group_hdf5 ('data')
        ! write components of f-array
        do pos=na,ne
          if (index_get(pos) == '') cycle
          call output_hdf5 ('data/'//index_get(pos), a(:,:,:,pos))
        enddo
      elseif (dataset == 'globals') then
        call create_group_hdf5 ('data')
        ! write components of global array
        do pos=na,ne
          if (index_get(mvar_io + pos) == '') cycle
          call output_hdf5 ('data/'//index_get(mvar_io + pos), a(:,:,:,pos))
        enddo
      else
        ! write other type of data array
        call output_hdf5 (dataset, a(:,:,:,na:ne), ne-na+1)
      endif
      call file_close_hdf5
!
      ! write additional settings
      call file_open_hdf5 (filename, global=.false., truncate=.false.)
      call output_settings (real (t), time_only=((.not. lwrite_add) .or. lomit_add_data))
      call file_close_hdf5
!
      call file_open_hdf5 (filename, truncate=.false.)
!
    endsubroutine output_snap
!***********************************************************************
    subroutine output_snap_finalize
!
!  Close snapshot file.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!  28-Oct-2016/PABourdin: redesigned
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
      call create_group_hdf5 ('part')
      call output_hdf5 (dataset, ap, mv, mparray, nv)
      call output_hdf5 ('part/ID', ipar(1:nv), nv)
      call create_group_hdf5 ('proc')
      call output_hdf5 ('proc/distribution', nv)
      call file_close_hdf5
!
      call file_open_hdf5 (filename, global=.false., truncate=.false.)
      ! write additional data
      call output_settings (real (t), time_only=((.not. (ltrunc .and. (trim (dataset) == 'fp'))) .or. lomit_add_data))
      ! write processor boundaries
      call output_hdf5 ('proc/bounds_x', procx_bounds(0:nprocx), nprocx+1)
      call output_hdf5 ('proc/bounds_y', procy_bounds(0:nprocy), nprocy+1)
      call output_hdf5 ('proc/bounds_z', procz_bounds(0:nprocz), nprocz+1)
      call file_close_hdf5
!
    endsubroutine output_part_snap
!***********************************************************************
    subroutine output_stalker_init(num, nv, snap, ID)
!
!  Open stalker particle snapshot file and initialize with snapshot time.
!
!  02-May-2019/PABourdin: coded
!
      use General, only: itoa
!
      integer, intent(in) :: num, nv, snap
      integer, dimension(nv), intent(in) :: ID
!
      character (len=fnlen) :: filename
      real :: t_sp
!
      filename = trim(directory_snap)//'/PSTALK'//trim(itoa(snap))//'.h5'
!
      ! open global HDF5 file and write particle data
      call file_open_hdf5 (filename, global=.true.)
      call create_group_hdf5 ('stalker')
      call output_hdf5 ('stalker/ID', ID, nv)
      call create_group_hdf5 ('proc')
      call output_hdf5 ('proc/distribution', nv)
      call file_close_hdf5
!
      if (lroot) then
        call file_open_hdf5 (filename, global=.false., truncate=.false.)
        t_sp = t
        call output_hdf5 ('time', t_sp)
        call file_close_hdf5
      endif
!
      call file_open_hdf5 (filename, global=.true., truncate=.false.)
!
    endsubroutine output_stalker_init
!***********************************************************************
    subroutine output_stalker(label, mv, nv, data)
!
!  Write stalker particle quantity to snapshot file.
!
!  02-May-2019/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: mv, nv
      real, dimension (mv), intent(in) :: data
!
      character (len=fnlen) :: dataset
!
      dataset = 'stalker/'//trim(label)
      call output_hdf5 (dataset, data(1:nv), nv)
!
    endsubroutine output_stalker
!***********************************************************************
    subroutine output_part_finalize
!
!  Close particle snapshot file.
!
!  02-May-2019/PABourdin: coded
!
      call file_close_hdf5
!
    endsubroutine output_part_finalize
!***********************************************************************
    subroutine output_settings(time, time_only)
!
!  Write additional settings and grid.
!
!  13-Nov-2018/PABourdin: moved from other functions
!
      use General, only: loptest
      use Mpicomm, only: collect_grid, mpi_precision
!
      real, optional, intent(in) :: time
      logical, optional, intent(in) :: time_only
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
!
      if (present (time)) call output_hdf5 ('time', time)
      if (loptest (time_only)) return
!
      if (lroot) then
        allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
      else
        allocate (gx(1), gy(1), gz(1), stat=alloc_err)
      endif
      if (alloc_err > 0) call fatal_error ('output_settings', 'allocate memory for gx,gy,gz', .true.)
!
      call create_group_hdf5 ('grid')
      call collect_grid (x, y, z, gx, gy, gz)
      call output_hdf5 ('grid/x', gx, mxgrid)
      call output_hdf5 ('grid/y', gy, mygrid)
      call output_hdf5 ('grid/z', gz, mzgrid)
      call output_hdf5 ('grid/dx', dx)
      call output_hdf5 ('grid/dy', dy)
      call output_hdf5 ('grid/dz', dz)
      call output_hdf5 ('grid/Lx', Lx)
      call output_hdf5 ('grid/Ly', Ly)
      call output_hdf5 ('grid/Lz', Lz)
      call output_hdf5 ('grid/Ox', x0)
      call output_hdf5 ('grid/Oy', y0)
      call output_hdf5 ('grid/Oz', z0)
      call collect_grid (dx_1, dy_1, dz_1, gx, gy, gz)
      call output_hdf5 ('grid/dx_1', gx, mxgrid)
      call output_hdf5 ('grid/dy_1', gy, mygrid)
      call output_hdf5 ('grid/dz_1', gz, mzgrid)
      call collect_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
      call output_hdf5 ('grid/dx_tilde', gx, mxgrid)
      call output_hdf5 ('grid/dy_tilde', gy, mygrid)
      call output_hdf5 ('grid/dz_tilde', gz, mzgrid)
      call create_group_hdf5 ('unit')
      call output_hdf5 ('unit/system', unit_system)
      call output_hdf5_double ('unit/density', unit_density)
      call output_hdf5_double ('unit/length', unit_length)
      call output_hdf5_double ('unit/velocity', unit_velocity)
      call output_hdf5_double ('unit/magnetic', unit_magnetic)
      call output_hdf5_double ('unit/temperature', unit_temperature)
      call output_hdf5_double ('unit/mass', unit_mass)
      call output_hdf5_double ('unit/energy', unit_energy)
      call output_hdf5_double ('unit/time', unit_time)
      call output_hdf5_double ('unit/flux', unit_flux)
      call create_group_hdf5 ('settings')
      call output_hdf5 ('settings/mx', nxgrid+2*nghost)
      call output_hdf5 ('settings/my', nygrid+2*nghost)
      call output_hdf5 ('settings/mz', nzgrid+2*nghost)
      call output_hdf5 ('settings/nx', nxgrid)
      call output_hdf5 ('settings/ny', nygrid)
      call output_hdf5 ('settings/nz', nzgrid)
      call output_hdf5 ('settings/l1', nghost)
      call output_hdf5 ('settings/m1', nghost)
      call output_hdf5 ('settings/n1', nghost)
      call output_hdf5 ('settings/l2', nghost+nxgrid-1)
      call output_hdf5 ('settings/m2', nghost+nygrid-1)
      call output_hdf5 ('settings/n2', nghost+nzgrid-1)
      call output_hdf5 ('settings/nghost', nghost)
      call output_hdf5 ('settings/mvar', mvar)
      call output_hdf5 ('settings/maux', maux)
      call output_hdf5 ('settings/mglobal', mglobal)
      call output_hdf5 ('settings/nprocx', nprocx)
      call output_hdf5 ('settings/nprocy', nprocy)
      call output_hdf5 ('settings/nprocz', nprocz)
      if (mpi_precision == MPI_REAL) then
        call output_hdf5 ('settings/precision', 'S')
      else
        call output_hdf5 ('settings/precision', 'D')
      endif
      ! versions represent only non-compatible file formats
      ! 0 : experimental
      ! 1 : first public release
      call output_hdf5 ('settings/version', 0)
!
      deallocate (gx, gy, gz)
!
    endsubroutine output_settings
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
      integer :: pos, last
      character (len=fnlen) :: filename, label
!
      if (.not. lroot) return
!
      filename = trim (directory_snap)//'/'//trim(file)//'.h5'
      call file_open_hdf5 (filename, global=.false.)
      call output_hdf5 ('number', mv)
      if (mv > 0) then
        call create_group_hdf5 ('points')
        do pos=1, nc
          label = trim (labels(pos))
          if (label(1:1) == 'i') then
            label = trim (label(2:))
            last = len (trim (label))
            if (label(last:last) == 'q') label = trim (label(1:last-1))
          endif
          call output_hdf5 ('points/'//trim(label), fq(:,pos), mv)
        enddo
      endif
      call output_hdf5 ('time', real (t))
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
    subroutine input_snap(file, a, nv, mode, label)
!
!  Read snapshot file. Also read mesh and time, if mode==1.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!  10-Mar-2015/MR: avoided use of fseek;
!                  this subroutine seems not yet to be adapted to HDF5
!  28-Oct-2016/PABourdin: redesigned
!
      use File_io, only: backskip_to_time
!
      character (len=*) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(out) :: a
      integer, optional, intent(in) :: mode
      character (len=*), optional :: label
!
      character (len=fnlen) :: dataset
      integer :: pos
!
      varfile_name = trim(directory_snap)//'/'//trim(file)//'.h5'
      dataset = 'f'
      if (present (label)) dataset = label
      if (dataset == 'globals') then
        if ((file(1:7) == 'timeavg') .or. (file(1:4) == 'TAVG')) then
          varfile_name = trim(datadir)//'/averages/'//trim(file)//'.h5'
        else
          varfile_name = trim(datadir_snap)//'/'//trim(file)//'.h5'
        endif
      endif
!
      lread_add = .true.
      if (present (mode)) lread_add = (mode == 1)
!
      ! open existing HDF5 file and read data
      call file_open_hdf5 (varfile_name, read_only=.true.)
      if (dataset == 'f') then
        ! read components of f-array
        do pos=1, nv
          if (index_get(pos) == '') cycle
          call input_hdf5 ('data/'//index_get(pos), a(:,:,:,pos))
        enddo
      elseif (dataset == 'globals') then
        ! read components of globals array
        do pos=1, nv
          if (index_get(mvar_io + pos) == '') cycle
          call input_hdf5 ('data/'//index_get(mvar_io + pos), a(:,:,:,pos))
        enddo
      else
        ! read other type of data array
        call input_hdf5 (dataset, a, nv)
      endif
!
    endsubroutine input_snap
!***********************************************************************
    subroutine input_snap_finalize
!
!  Close snapshot file.
!
!  12-Oct-2019/PABourdin: moved code from 'input_snap'
!
      use Mpicomm, only: mpibcast_real, MPI_COMM_WORLD
!
      real :: time
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
!
      call file_close_hdf5
      persist_initialized = .false.
!
      ! read additional data
      if (lread_add .and. lomit_add_data) then
        call file_open_hdf5 (varfile_name, global=.false., read_only=.true.)
        call input_hdf5 ('time', time)
        call file_close_hdf5
!
        call mpibcast_real (time, comm=MPI_COMM_WORLD)
        t = time
      elseif (lread_add) then
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
        else
          allocate (gx(1), gy(1), gz(1), stat=alloc_err)
        endif
        if (alloc_err > 0) call fatal_error ('input_snap', 'Could not allocate memory for gx,gy,gz', .true.)
!
        call file_open_hdf5 (varfile_name, global=.false., read_only=.true.)
        call input_hdf5 ('time', time)
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
        call mpibcast_real (time, comm=MPI_COMM_WORLD)
        t = time
        call mpibcast_real (dx, comm=MPI_COMM_WORLD)
        call mpibcast_real (dy, comm=MPI_COMM_WORLD)
        call mpibcast_real (dz, comm=MPI_COMM_WORLD)
        call mpibcast_real (Lx, comm=MPI_COMM_WORLD)
        call mpibcast_real (Ly, comm=MPI_COMM_WORLD)
        call mpibcast_real (Lz, comm=MPI_COMM_WORLD)
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
      call input_hdf5 ('proc/distribution', nv)
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
      integer :: pos, mv_in, last
      character (len=fnlen) :: filename, label
!
      filename = trim (directory_snap)//'/'//trim(file)//'.h5'
      if (lroot) then
        call file_open_hdf5 (filename, read_only=.true., global=.false.)
        call input_hdf5 ('number', mv_in)
        if (mv_in /= mv) call fatal_error ("input_pointmass", "Number of points seems incorrect.", .true.)
        if (mv_in > 0) then
          do pos=1, nc
            label = trim (labels(pos))
            if (label(1:1) == 'i') then
              label = trim (label(2:))
              last = len (trim (label))
              if (label(last:last) == 'q') label = trim (label(1:last-1))
            endif
            call input_hdf5 ('points/'//trim(label), fq(:,pos), mv)
          enddo
        endif
        call file_close_hdf5
      endif
!
      if (mv > 0) call mpibcast_real (fq, (/ nqpar, mqarray /))
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
      call create_group_hdf5 ('persist')
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
      use General, only: lower_case
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
      call output_hdf5 ('persist/'//lower_case (label), value_int, size (value), same_size=.true.)
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
      use General, only: lower_case
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, dimension (:), intent(in) :: value
!
      write_persist_int_1D = .true.
      if (write_persist_id (label, id)) return
!
      call output_hdf5 ('persist/'//lower_case (label), value, size (value), same_size=.true.)
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
      use General, only: lower_case
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, dimension (:), intent(in) :: value
!
      write_persist_real_1D = .true.
      if (write_persist_id (label, id)) return
!
      call output_hdf5 ('persist/'//lower_case (label), value, size (value), same_size=.true.)
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
        call file_close_hdf5
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
    logical function persist_exists(label)
!
!  Check if a persistent variable exists.
!
!  12-Oct-2019/PABourdin: coded
!
      use General, only: lower_case
!
      character (len=*), intent(in) :: label
!
      persist_exists = exists_in_hdf5('persist')
      if (persist_exists) persist_exists = exists_in_hdf5('persist/'//lower_case (label))
!
    endfunction persist_exists
!***********************************************************************
    logical function read_persist_id(label, id, lerror_prone)
!
!  Read persistent block ID from snapshot file.
!
!  27-Oct-2018/PABourdin: coded
!
      use General, only: lower_case
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
      id = -1
      lexists = exists_in_hdf5('persist')
      if (lexists) lexists = exists_in_hdf5('persist/'//lower_case (label))
      if (lexists) id = 0
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
      use General, only: lower_case
!
      character (len=*), intent(in) :: label
      logical, dimension(:), intent(out) :: value
!
      integer, dimension(size (value)) :: value_int
!
      read_persist_logical_1D = .true.
      if (.not. persist_exists (label)) return
!
      call input_hdf5 ('persist/'//lower_case (label), value_int, size (value), same_size=.true.)
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
      use General, only: lower_case
!
      character (len=*), intent(in) :: label
      integer, dimension(:), intent(out) :: value
!
      read_persist_int_1D = .true.
      if (.not. persist_exists (label)) return
!
      call input_hdf5 ('persist/'//lower_case (label), value, size (value), same_size=.true.)
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
      use General, only: lower_case
!
      character (len=*), intent(in) :: label
      real, dimension(:), intent(out) :: value
!
      read_persist_real_1D = .true.
      if (.not. persist_exists (label)) return
!
      call input_hdf5 ('persist/'//lower_case (label), value, size (value), same_size=.true.)
!
      read_persist_real_1D = .false.
!
    endfunction read_persist_real_1D
!***********************************************************************
    subroutine output_globals(file, a, nv, label)
!
!  Write snapshot file of globals, ignore time and mesh.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!  29-Nov-2018/PABourdin: extended for output of time averages
!
      use General, only : coptest
!
      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*), intent(in), optional :: label
!
      call output_snap (a, nv, file=file, mode=0, label=coptest(label,'globals'))
      call output_snap_finalize
!
    endsubroutine output_globals
!***********************************************************************
    subroutine input_globals(file, a, nv)
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
      if (dir == '.') then
        if (flist(1:4) == 'tavg') then
          call safe_character_assign (dir, trim(datadir)//'/averages')
        else
          call safe_character_assign (dir, directory_snap)
        endif
      endif
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
      use File_io, only: file_exists
!
      character (len=*), intent(in) :: file
      integer, intent(in), optional :: mxout,myout,mzout
!
      character (len=fnlen) :: filename
      logical :: lexists
!
      if (lyang) return      ! grid collection only needed on Yin grid, as grids are identical
!
      filename = trim (datadir)//'/grid.h5'
      lexists = file_exists (filename)
      call file_open_hdf5 (filename, global=.false., truncate=.not. lexists)
      call output_settings
      call file_close_hdf5
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
      call mpibcast_real (dx, comm=MPI_COMM_WORLD)
      call mpibcast_real (dy, comm=MPI_COMM_WORLD)
      call mpibcast_real (dz, comm=MPI_COMM_WORLD)
      call mpibcast_real (Lx, comm=MPI_COMM_WORLD)
      call mpibcast_real (Ly, comm=MPI_COMM_WORLD)
      call mpibcast_real (Lz, comm=MPI_COMM_WORLD)
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
