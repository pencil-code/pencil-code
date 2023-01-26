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
!
  interface input_proc_bounds
    module procedure input_proc_bounds_double
    module procedure input_proc_bounds_single
  endinterface
!
! Indicates if IO is done distributed (each proc writes into a procdir)
! or collectively (eg. by specialized IO-nodes or by MPI-IO).
!
  logical :: lcollective_IO=.true.
  character (len=labellen) :: IO_strategy="HDF5"
!
  character (len=fnlen) :: last_snapshot = ""
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
      character(LEN=5) :: locking

      if (lroot) call svn_id ("$Id$")
!
      if (lread_from_other_prec) &
        call warning('register_io','Reading from other precision not implemented')
!
      lmonolithic_io = .true.

      call getenv('HDF5_USE_FILE_LOCKING',locking)
      if (trim(locking)/='FALSE') &
        call warning('register_io','HDF5 files are possibly locked; writing may fail.')
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
      use General, only: coptest, ioptest, itoa
      use File_io, only: parallel_file_exists
!
      integer, optional, intent(in) :: nv1,nv2
      real, dimension (:,:,:,:), intent(in) :: a
      character (len=*), intent(in) :: file
      integer, optional, intent(in) :: mode
      logical, optional, intent(in) :: ltruncate
      character (len=*), optional, intent(in) :: label
!
      integer :: pos,na,ne
      logical :: ltrunc, lexists, lwrite_add
      character (len=fnlen) :: filename, dataset, group
!
      dataset = coptest(label,'f')
      if (dataset == 'globals') then
        filename = trim(datadir_snap)//'/'//trim(file)//'.h5'
      elseif (dataset == 'timeavg') then
        filename = trim(datadir)//'/averages/'//trim(file)//'.h5'
      else
        filename = trim(directory_snap)//'/'//trim(file)//'.h5'
      endif
      if (present (ltruncate)) then
        ltrunc = ltruncate
      else
        lexists = parallel_file_exists(filename)
        ltrunc = .not. lexists
      endif
!
      lwrite_add = .true.
      if (present (mode)) lwrite_add = (mode == 1)
!
!  nv1 and nv2 must always be within the range of registered variables, which are
!  eligible for storage in snapshots
!
      if (ioptest(nv1)>mvar+maux) &
          call fatal_error ('output_snap', "nv1 outside legal range mvar+maux")
      if (ioptest(nv2)>mvar+maux) &
          call fatal_error ('output_snap', "nv2 outside legal range mvar+maux")
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
          group=index_get(pos)
          if (group == '') cycle
          call output_hdf5 ('data/'//trim(group), a(:,:,:,pos))
        enddo
      elseif (dataset == 'globals') then
        if (.not. present (nv1)) &
          call fatal_error ('output_snap', "for dataset == 'globals' nv1 needs to be set")
        call create_group_hdf5 ('data')
        ! write components of global array
        do pos=1,nv1
          group=index_get(mvar_io + pos)
          if (group == '') cycle
          call output_hdf5 ('data/'//trim(group), a(:,:,:,pos))
        enddo
      else
        ! write other type of data array
        call output_hdf5 (dataset, a(:,:,:,na:ne), ne-na+1)
      endif
      call file_close_hdf5
!
      ! write additional settings
      call file_open_hdf5 (filename, global=.false., truncate=.false.)
      call output_settings (real(t), time_only=((.not. lwrite_add) .or. lomit_add_data))
      call file_close_hdf5
!
!      call file_open_hdf5 (filename, truncate=.false.)
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
    subroutine output_average_2D(label, nc, name, data, time, lwrite, header)
!
!   Output average to a file.
!
!   16-Nov-2018/PABourdin: coded
!
      use File_io, only: parallel_file_exists
      use General, only: itoa
      use Mpicomm, only: mpibcast_int
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: nc
      character (len=fmtlen), dimension(nc), intent(in) :: name
      real, dimension(:,:,:), intent(in), target :: data
      real, intent(in) :: time
      logical, intent(in) :: lwrite
      real, dimension(:), optional, intent(in) :: header
!
      character (len=fnlen) :: filename
      character (len=intlen) :: group
      integer :: last, ia, alloc_err
      logical :: lexists, lglobal
      real, dimension (:,:), allocatable, target :: component
      real, pointer :: local_data(:,:)
!
      if (nc <= 0) return
!
      filename = trim(datadir)//'/averages/'//trim(label)//'.h5'
      lexists = parallel_file_exists (filename)
      lglobal = ((label == 'z') .or. (label == 'y'))
!
      last = 0
      if (lexists) then
        call file_open_hdf5 (filename, global=.false., read_only=.true.)
        if (exists_in_hdf5 ('last')) then
          call input_hdf5 ('last', last)
          last = last + 1
        endif
        call file_close_hdf5
        call mpibcast_int (last)
      endif
!
      call file_open_hdf5 (filename, global=lglobal, truncate=(.not. lexists), write=lwrite)
      if ((last == 0) .and. present (header)) call output_hdf5 ('r', header, size(header))
      group = trim(itoa(last))
      call create_group_hdf5 (group)
      if (label == 'z') then
        allocate (component(size(data,2),size(data,3)), stat = alloc_err)
        if (alloc_err > 0) call fatal_error('output_average_2D', 'Could not allocate memory for component')
        do ia = 1, nc
          component = data(ia,:,:)
          local_data => component
          call output_hdf5 (trim(group)//'/'//trim(name(ia)), local_data, nx, ny, nxgrid, nygrid, ipx, ipy, lwrite)
        enddo
        deallocate (component)
      elseif (label == 'y') then
        do ia = 1, nc
          local_data => data(:,:,ia)
          call output_hdf5 (trim(group)//'/'//trim(name(ia)), local_data, nx, nz, nxgrid, nzgrid, ipx, ipz, lwrite)
        enddo
      else
        do ia = 1, nc
          call output_hdf5 (trim(group)//'/'//trim(name(ia)), data(:,:,ia), size(data,1), size(data,2))
        enddo
      endif
      if (lglobal) then
        call file_close_hdf5
        call file_open_hdf5 (filename, global=.false., truncate=.false.)
      endif
      call output_hdf5 (trim(group)//'/time', time)
      call output_hdf5 ('last', last)
      call file_close_hdf5
!
    endsubroutine output_average_2D
!***********************************************************************
    subroutine output_slice_position
!
!  Record slice positions.
!
!  13-nov-20/ccyang: wrapper
!
      use HDF5_IO, only: hdf5_output_slice_position
!
      call hdf5_output_slice_position
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
      if (present (ltruncate)) then
        ltrunc = ltruncate
      else
        lexists = parallel_file_exists(filename)
        ltrunc = .not. lexists
      endif
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
    subroutine output_stalker(label, mv, nv, data, nvar, lfinalize)
!
!  Write stalker particle quantity to snapshot file.
!
!  02-May-2019/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: mv, nv
      real, dimension (mv), intent(in) :: data
      logical, intent(in), optional :: lfinalize
      integer, intent(in), optional :: nvar
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
        if (nc > 0) then
          pos = 1
          do while (pos <= nc)
            label = trim (labels(pos))
            if (label(1:1) == 'i') then
              label = trim (label(2:))
              last = len (trim (label))
              if (label(last:last) == 'q') label = trim (label(1:last-1))
            endif
            call output_hdf5 ('points/'//trim(label), fq(:,pos), mv)
            pos = pos + 1
          enddo
        endif
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
      use Syscalls, only: islink
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
      if (islink(varfile_name)) snaplink=varfile_name
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
      use Mpicomm, only: mpibcast_real, mpibcast_logical, MPI_COMM_WORLD
      use Syscalls, only: system_cmd
!
      real :: time
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
      logical :: lerrcont
!
      call file_close_hdf5
      persist_initialized = .false.
!
! Read additional data.
!
      if (lread_add) then
!
!  Time is always read.
!
        if (lroot) then
          call file_open_hdf5 (varfile_name, global=.false., read_only=.true.)
          lerrcont=.true.
          call input_hdf5 ('time', time, lerrcont)
          call file_close_hdf5
          if (lerrcont) call recover_time_from_series(time)
        endif
!
        call mpibcast_real (time, comm=MPI_COMM_WORLD)
        t = time
!
!  Read further data if not to be omitted.
!
        if (.not.lomit_add_data) then
!
          if (lroot) then
            allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
            if (alloc_err > 0) call fatal_error ('input_snap', 'Could not allocate memory for gx,gy,gz', .true.)
          endif

          lerrcont=.false.
          if (lroot) then
            call file_open_hdf5 (varfile_name, global=.false., read_only=.true.)
            lerrcont=.true.
            call input_hdf5 ('grid/x', gx, mxgrid, lerrcont=lerrcont)
            if (lerrcont) goto 100
            lerrcont=.true.
            call input_hdf5 ('grid/y', gy, mygrid, lerrcont=lerrcont)
            if (lerrcont) goto 100
            lerrcont=.true.
            call input_hdf5 ('grid/z', gz, mzgrid, lerrcont=lerrcont)
            if (lerrcont) goto 100
          endif
          call distribute_grid (x, y, z, gx, gy, gz)
          if (lroot) then
            lerrcont=.true.
            call input_hdf5 ('grid/dx', dx, lerrcont)
            if (lerrcont) goto 100
            lerrcont=.true.
            call input_hdf5 ('grid/dy', dy, lerrcont)
            if (lerrcont) goto 100
            lerrcont=.true.
            call input_hdf5 ('grid/dz', dz, lerrcont)
            if (lerrcont) goto 100
            lerrcont=.true.
            call input_hdf5 ('grid/Lx', Lx)
            lerrcont=.true.
            call input_hdf5 ('grid/Ly', Ly)
            lerrcont=.true.
            call input_hdf5 ('grid/Lz', Lz)
            lerrcont=.true.
            call input_hdf5 ('grid/dx_1', gx, mxgrid, lerrcont=lerrcont)
            if (lerrcont) goto 100
            lerrcont=.true.
            call input_hdf5 ('grid/dy_1', gy, mygrid, lerrcont=lerrcont)
            if (lerrcont) goto 100
            lerrcont=.true.
            call input_hdf5 ('grid/dz_1', gz, mzgrid, lerrcont=lerrcont)
            if (lerrcont) goto 100
          endif
          call distribute_grid (dx_1, dy_1, dz_1, gx, gy, gz)
          if (lroot) then
            lerrcont=.true.
            call input_hdf5 ('grid/dx_tilde', gx, mxgrid, lerrcont=lerrcont)
            if (lerrcont) goto 100
            lerrcont=.true.
            call input_hdf5 ('grid/dy_tilde', gy, mygrid, lerrcont=lerrcont)
            if (lerrcont) goto 100
            lerrcont=.true.
            call input_hdf5 ('grid/dz_tilde', gz, mzgrid, lerrcont=lerrcont)
            if (lerrcont) goto 100
          endif
          call distribute_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)

          call mpibcast_real (dx, comm=MPI_COMM_WORLD)
          call mpibcast_real (dy, comm=MPI_COMM_WORLD)
          call mpibcast_real (dz, comm=MPI_COMM_WORLD)
          call mpibcast_real (Lx, comm=MPI_COMM_WORLD)
          call mpibcast_real (Ly, comm=MPI_COMM_WORLD)
          call mpibcast_real (Lz, comm=MPI_COMM_WORLD)

100       if (lroot) call file_close_hdf5
          call mpibcast_logical(lerrcont, comm=MPI_COMM_WORLD)
          if (lerrcont) then
            call warning('input_snap_finalize','grid data corrupted, reading grid from grid.h5')
            call rgrid('')
          endif
!
          if (snaplink/='') then
            call system_cmd('rm -f '//snaplink)
            snaplink=''
          endif
!
        endif
      endif
!
contains
!
!------------------------------------------------------------------------------------------------
      subroutine recover_time_from_series(time)
!
!  Requires the start and end positions of t and dt in a record of time_series.dat to be stored in
!  the second record of data/legend.dat in the order start(t), end(t), start(dt), end(dt).
!  If nothing can be read, time is set to zero.
!
        real, intent(OUT) :: time

        integer, dimension(2) :: itpos=(/0,0/), idtpos=(/0,0/)
        !character(LEN=:), allocatable :: ctime
        character(LEN=256) :: ctime
        real :: dtime

        if (lroot) then
          time=0.; dtime=0.
          open(11,file='data/legend.dat', position='append')
          backspace 11
          read(11,*,err=100,end=100) itpos, idtpos
          close(11)
          open(11,file='data/time_series.dat', position='append')
         ! if (max(itpos(2), idtpos(2))>0) allocate(character(LEN=max(itpos(2), idtpos(2))) :: ctime)
          if (itpos(1)/=0) then
            backspace 11
            read(11,'(a)') ctime(1:itpos(2))
            read(ctime(itpos(1):itpos(2)),*) time
          endif
          if (idtpos(1)/=0) then
            backspace 11
            read(11,'(a)') ctime(1:idtpos(2))
            read(ctime(idtpos(1):idtpos(2)),*) dtime
          endif
100       close(11)
          time=time+dtime
          write(ctime,'(e15.7)') time
          call warning('input_snap_finalize', 'snapshot corrupted; time was set to '//trim(ctime))

          !call system_cmd("tail -n 1 data/time_series.dat | tac | "// &
          !                "sed -e's/^ *[0-9][0-9]*  *\([0-9][0-9]*\.[0-9E+-][0-9E+-]*  *[0-9]\.[0-9E+-][0-9E+-]*\) *.*$/\1/'"// &
          !                " > time.tmp")
        endif

       endsubroutine recover_time_from_series

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
      if (lroot) then
        call file_open_hdf5 (filename, global=.false., read_only=.true.)
        call input_hdf5 ('proc/bounds_x', procx_bounds, nprocx+1)
        call input_hdf5 ('proc/bounds_y', procy_bounds, nprocy+1)
        call input_hdf5 ('proc/bounds_z', procz_bounds, nprocz+1)
        call file_close_hdf5
      endif
      call mpibcast (procx_bounds, nprocx+1)
      call mpibcast (procy_bounds, nprocy+1)
      call mpibcast (procz_bounds, nprocz+1)
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
    logical function write_persist_torus_rect(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  13-May-2020/MR: coded
!
      use Geometrical_types
      use General, only: lower_case
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      type(torus_rect), intent(in) :: value
!
      write_persist_torus_rect = .true.
      if (write_persist_id (label, id)) return
!
      call output_hdf5 ('persist/'//lower_case (label), value)
      write_persist_torus_rect = .false.
!
    endfunction write_persist_torus_rect
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
      if (persist_exists) persist_exists = exists_in_hdf5('persist/'//lower_case(label))
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
    logical function read_persist_torus_rect(label,value)
!
!  Read persistent data from snapshot file.
!
!  16-May-2020/MR: coded
!
      use Geometrical_types

      character (len=*), intent(in) :: label
      type(torus_rect) :: value    !, intent(out) :: value
!
      read_persist_torus_rect = .false.
!
    endfunction read_persist_torus_rect
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
        write (lun_output, '(A,1x,e16.8)') trim (fpart), t
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
    subroutine wgrid(file,mxout,myout,mzout,lwrite)
!
!  Write grid coordinates.
!
!  27-Oct-2018/PABourdin: coded
!
      use File_io, only: file_exists
      use General, only: loptest
!
      character (len=*), intent(in) :: file
      integer, intent(in), optional :: mxout,myout,mzout
      logical, optional :: lwrite
!
      character (len=fnlen) :: filename
      logical :: lexists
!
      if (lyang) return      ! grid collection only needed on Yin grid, as grids are identical
!
      if (loptest(lwrite,.not.luse_oldgrid)) then
        filename = trim (datadir)//'/grid.h5'
        lexists = file_exists (filename)
        call file_open_hdf5 (filename, global=.false., truncate=.not. lexists)
        call output_settings
        call file_close_hdf5
      endif
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
      character (len=*) :: file         ! not used
!
      character (len=fnlen) :: filename
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
!
      if (lroot) then
        allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('rgrid', 'Could not allocate memory for gx,gy,gz', .true.)

        filename = trim (datadir)//'/grid.h5'
        call file_open_hdf5 (filename, global=.false., read_only=.true.)
        call input_hdf5 ('grid/x', gx, mxgrid)
        call input_hdf5 ('grid/y', gy, mygrid)
        call input_hdf5 ('grid/z', gz, mzgrid)
        call input_hdf5 ('grid/dx', dx)
        call input_hdf5 ('grid/dy', dy)
        call input_hdf5 ('grid/dz', dz)
        call input_hdf5 ('grid/Lx', Lx)
        call input_hdf5 ('grid/Ly', Ly)
        call input_hdf5 ('grid/Lz', Lz)
      endif
      call distribute_grid (x, y, z, gx, gy, gz)
      if (lroot) then
        call input_hdf5 ('grid/dx_1', gx, mxgrid)
        call input_hdf5 ('grid/dy_1', gy, mygrid)
        call input_hdf5 ('grid/dz_1', gz, mzgrid)
      endif
      call distribute_grid (dx_1, dy_1, dz_1, gx, gy, gz)
      if (lroot) then
        call input_hdf5 ('grid/dx_tilde', gx, mxgrid)
        call input_hdf5 ('grid/dy_tilde', gy, mygrid)
        call input_hdf5 ('grid/dz_tilde', gz, mzgrid)
        call file_close_hdf5
      endif
      call distribute_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
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
