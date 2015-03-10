! $Id$
!
!  I/O via the HDF5 hyperslab-by-chunk IO routines.
!  (storing data into one file, e.g. data/allprocs/var.h5)
!
!  The data format is self-contained. Only outer ghost-layers are stored.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_hdf5.f90
!
module Io
!
  use Cdata
  use Cparam, only: intlen, fnlen, max_int
  use Messages, only: fatal_error, svn_id
!
  implicit none
!
  include 'io.h'
  include 'record_types.h'
!
  interface output_form
    module procedure output_form_int_0D
  endinterface
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
  logical :: persist_initialized=.false.
  integer :: persist_last_id=-max_int
!
  integer :: local_type, global_type, h5_err
  integer(HID_T) :: h5_file, h5_dset, h5_plist, h5_fspace, h5_mspace
  integer, parameter :: n_dims=3
  integer, dimension(n_dims+1) :: local_size, local_subsize, local_start
  integer, dimension(n_dims+1) :: global_size, global_subsize, global_start
!
  contains
!***********************************************************************
    subroutine register_io()
!
!  dummy routine, generates separate directory for each processor.
!  VAR#-files are written to the directory directory_snap which will
!  be the same as directory, unless specified otherwise.
!
!  04-Jul-2011/Boudin.KIS: coded
!
      use Mpicomm, only: stop_it_if_any
!
!  identify version number
!
      if (lroot) call svn_id ("$Id$")
      if (.not. lseparate_persist) call fatal_error ('io_HDF5', &
          "This module only works with the setting lseparate_persist=.true.")
!
! Create datatype to describe internal elements of data, ie. the core data
! excluding the halos, unless we are on an edge and have to include them.
!
      local_size(1) = mx
      local_size(2) = my
      local_size(3) = mz
!
      local_subsize(1) = nx
      local_subsize(2) = ny
      local_subsize(3) = nz
!
! We need to save the outermost halos if we are on either edge.
!
      if (lfirst_proc_x) local_subsize(1) = local_subsize(1) + nghost
      if (lfirst_proc_y) local_subsize(2) = local_subsize(2) + nghost
      if (lfirst_proc_z) local_subsize(3) = local_subsize(3) + nghost
!
      if (llast_proc_x)  local_subsize(1) = local_subsize(1) + nghost
      if (llast_proc_y)  local_subsize(2) = local_subsize(2) + nghost
      if (llast_proc_z)  local_subsize(3) = local_subsize(3) + nghost
!
! The displacements in 'local_start' uses C-like format, ie. start from zero.
!
      local_start(1) = l1 - 1
      local_start(2) = m1 - 1
      local_start(3) = n1 - 1
      local_start(4) = 0

! We need to include lower ghost cells if we are on a lower edge
! inclusion of upper ghost cells is taken care of by increased subsize.

      if (lfirst_proc_x) local_start(1) = local_start(1) - nghost
      if (lfirst_proc_y) local_start(2) = local_start(2) - nghost
      if (lfirst_proc_z) local_start(3) = local_start(3) - nghost
!
! Now define the position of this processors data in the global file:
!
! Create datatype to describe the section of the global file that
! is owned by this process (the ftype). 'global_size' has now the global
! sizes, but 'global_subsize' stays the same. 'global_start' must again count
! in C-manner, ie. from 0.
!
! Global size of array·
!
      global_size(1) = mxgrid
      global_size(2) = mygrid
      global_size(3) = mzgrid
!
! Starting position of this processors data portion.
!
      global_start(1) = nghost + ipx*nx
      global_start(2) = nghost + ipy*ny
      global_start(3) = nghost + ipz*nz
      global_start(4) = 0
!
! Take account of inclusion of lower halos on lower edges.
!
      if (lfirst_proc_x) global_start(1) = global_start(1) - nghost
      if (lfirst_proc_y) global_start(2) = global_start(2) - nghost
      if (lfirst_proc_z) global_start(3) = global_start(3) - nghost
!
    endsubroutine register_io
!***********************************************************************
    subroutine directory_names()
!
!  Set up the directory names:
!  set directory name for the output (one subdirectory for each processor)
!  if datadir_snap (where var.dat, VAR# go) is empty, initialize to datadir
!
!  02-oct-2002/wolf: coded
!
      use General, only: safe_character_assign, itoa
!
      character (len=intlen) :: chproc
!
!  check whether directory_snap contains `/allprocs' -- if so, revert to the
!  default name.
!  Rationale: if directory_snap was not explicitly set in start.in, it
!  will be written to param.nml as 'data/allprocs'.
!
      if ((datadir_snap == '') .or. (index(datadir_snap,'allprocs')>0)) then
        datadir_snap = datadir
      endif
!
      chproc = itoa (iproc)
      call safe_character_assign (directory, trim (datadir)//'/proc'//chproc)
      call safe_character_assign (directory_dist, &
                                            trim (datadir_snap)//'/proc'//chproc)
      call safe_character_assign (directory_snap, trim (datadir_snap)//'/allprocs')
      call safe_character_assign (directory_collect, trim (datadir_snap)//'/allprocs')
!
    endsubroutine directory_names
!***********************************************************************
    subroutine collect_grid(x, y, z, gx, gy, gz)
!
!  This routine collects the global grid on the root processor.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      real, dimension(mx), intent(in) :: x
      real, dimension(my), intent(in) :: y
      real, dimension(mz), intent(in) :: z
      real, dimension(nxgrid+2*nghost), intent(out), optional :: gx
      real, dimension(nygrid+2*nghost), intent(out), optional :: gy
      real, dimension(nzgrid+2*nghost), intent(out), optional :: gz
!
      real, dimension(l2) :: buf_x
      real, dimension(m2) :: buf_y
      real, dimension(n2) :: buf_z
      integer :: px, py, pz
      integer, parameter :: tag_gx=677, tag_gy=678, tag_gz=679
!
      if (lroot) then
        ! collect the global x-data from all leading processors in the yz-plane
        gx(1:l2) = x(1:l2)
        if (nprocx > 1) then
          do px = 1, nprocx-1
            call mpirecv_real (buf_x, l2, px, tag_gx)
            gx(px*nx+l1:px*nx+mx) = buf_x
          enddo
        endif
        ! collect the global y-data from all leading processors in the xz-plane
        gy(1:m2) = y(1:m2)
        if (nprocy > 1) then
          do py = 1, nprocy-1
            call mpirecv_real (buf_y, m2, py*nprocx, tag_gy)
            gy(py*ny+m1:py*ny+my) = buf_y
          enddo
        endif
        ! collect the global z-data from all leading processors in the xy-plane
        gz(1:n2) = z(1:n2)
        if (nprocz > 1) then
          do pz = 1, nprocz-1
            call mpirecv_real (buf_z, n2, pz*nprocxy, tag_gz)
            gz(pz*nz+n1:pz*nz+mz) = buf_z
          enddo
        endif
      else
        ! leading processors send their local coordinates
        if (lfirst_proc_yz) call mpisend_real (x(l1:mx), l2, 0, tag_gx)
        if (lfirst_proc_xz) call mpisend_real (y(m1:my), m2, 0, tag_gy)
        if (lfirst_proc_xy) call mpisend_real (z(n1:mz), n2, 0, tag_gz)
      endif
!
    endsubroutine collect_grid
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
    subroutine output_snap(a, nv, file, mode)
!
!  Write snapshot file, always write mesh and time, could add other things.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!  13-feb-2014/MR: made file optional (prep for downsampled output)
!
      use Mpicomm, only: globalize_xy, mpirecv_real, mpi_precision
!
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: a
      character (len=*), optional, intent(in) :: file
      integer, optional, intent(in) :: mode
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
      integer, dimension (n_dims+1) :: stride, count
      integer, dimension(MPI_STATUS_SIZE) :: status
      logical :: lwrite_add
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (.not.present(file)) call fatal_error('output_snap', &
          'downsampled output not implemented for IO_hdf5')
!
      lwrite_add = .true.
      if (present (mode)) lwrite_add = (mode == 1)
!
! Initialize parallel HDF5 Fortran libaray.
!
      call h5open_f (h5_err)
      call stop_it_if_any (h5_err /= 0, "io_HDF5: can't initialize parallel HDF5 library")
!
! Setup file access property list.
!
      call h5pcreate_f (H5P_FILE_ACCESS_F, h5_plist, h5_err)
      call stop_it_if_any (h5_err /= 0, "io_HDF5: can't create file access property list")
      call h5pset_fapl_mpio_f (h5_plist, MPI_COMM_WORLD, MPI_INFO_NULL, h5_err)
      call stop_it_if_any (h5_err /= 0, "io_HDF5: can't initialize file access property list")
!
! Create (and truncate) HDF5 file.
!
      call h5fcreate_f (trim (directory_snap)//'/'//file, H5F_ACC_TRUNC_F, h5_file, h5_err, access_prp=h5_plist)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not open: "'//trim (directory_snap)//'/'//file//'"')
      call h5pclose_f (h5_plist, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not close parameter list')
!
! Define 'file-space' to indicate the data portion in the global file.
!
      global_size(4) = nv
      call h5screate_simple_f (n_dims+1, global_size, h5_fspace, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not create global file space')
!
! Define 'memory-space' to indicate the local data portion in memory.
!
      local_size(4) = nv
      call h5screate_simple_f (n_dims+1, local_size, h5_mspace, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not create local memory space')
!
! Create the dataset.
!
      call h5pcreate_f (H5P_DATASET_CREATE_F, h5_plist, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not prepare property list')
      call h5pset_chunk_f (h5_plist, n_dims+1, local_size, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not set chunk properties')
      call h5dcreate_f (h5_file, file, H5T_NATIVE_INTEGER, h5_fspace, h5_dset, h5_err, h5_plist)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not create dataset')
      call h5sclose_f (h5_fspace, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not close global file space')
!
! Define local 'hyper-slab' in the global file.
!
      local_subsize(4) = nv
      stride(*) = 1
      count(*) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not get dataset for file space')
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start, count, h5_err, h5_stride, local_size)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not select hyperslab within file')
!
! Prepare data transfer.
!
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not set data transfer properties')
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not select collective IO')
!
! Collectively write the data.
!
      call h5dwrite_f (h5_dset, H5T_NATIVE_INTEGER, a, global_size, h5_err, &
          file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not write the data')
!
      ! write additional data:
      if (lwrite_add) then
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('output_snap', 'Could not allocate memory for gx,gy,gz', .true.)
          call collect_grid (x, y, z, gx, gy, gz)
!
          open (lun_output, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', position='append', status='old')
          t_sp = t
          write (lun_output) t_sp, gx, gy, gz, dx, dy, dz
          close (lun_output)
          deallocate (gx, gy, gz)
        else
          call collect_grid (x, y, z)
        endif
      endif
!
    endsubroutine output_snap
!***********************************************************************
    subroutine output_snap_finalize()
!
!  Close snapshot file.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      if (persist_initialized) then
        if (lroot .and. (ip <= 9)) write (*,*) 'finish persistent block'
        write (lun_output) id_block_PERSISTENT
        persist_initialized = .false.
        close (lun_input)
      endif
!
! Close data spaces, dataset, property list, and the file itself.
!
      call h5sclose_f (h5_fspace, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not close the file space')
      call h5sclose_f (h5_mspace, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not close the memory space')
      call h5dclose_f (h5_dset, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not close the dataset')
      call h5pclose_f (h5_plist, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not close the parameter list')
      call h5fclose_f (h5_file, h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not close the file')
      call hh5close_f (h5_err)
      call stop_it_if_any (h5_err /= 0, 'output_snap: Could not close the parallel HDF5 library')
!
    endsubroutine output_snap_finalize
!***********************************************************************
    subroutine input_snap(file, a, nv, mode)
!
!  Read snapshot file. Also read mesh and time, if mode==1.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!  10-Mar-2015/MR: avoided use of fseek;
!                  this subroutine seems not yet to be adapted to HDF5
!
      use Mpicomm, only: localize_xy, mpibcast_real, mpi_precision
      use General, only: backskip_to_time
!
      character (len=*) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(out) :: a
      integer, optional, intent(in) :: mode
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: comm, handle, alloc_err, io_info=MPI_INFO_NULL
      integer, dimension(MPI_STATUS_SIZE) :: status
      logical :: lread_add
      real :: t_sp   ! t in single precision for backwards compatibility
!
      lread_add = .true.
      if (present (mode)) lread_add = (mode == 1)
!
! Define 'local_type' to be the local data portion that is being saved.
!
      local_size(4) = nv
      local_subsize(4) = nv
      call MPI_TYPE_CREATE_SUBARRAY(n_dims+1, local_size, local_subsize, &
          local_start, order, mpi_precision, local_type, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error ('input_snap', &
          'Could not create local subarray: "'//file//'"')
      call MPI_TYPE_COMMIT(local_type, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error ('input_snap', &
          'Could not commit local subarray: "'//file//'"')
!
! Define 'global_type' to indicate the local data portion in the global file.
!
      global_size(4) = nv
      call MPI_TYPE_CREATE_SUBARRAY(n_dims+1, global_size, local_subsize, &
          global_start, order, mpi_precision, global_type, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error ('input_snap', &
          'Could not create global subarray: "'//file//'"')
      call MPI_TYPE_COMMIT(global_type, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error ('input_snap', &
          'Could not commit global subarray: "'//file//'"')
!
      call MPI_FILE_OPEN(comm, trim (directory_snap)//'/'//file, &
          MPI_MODE_CREATE+MPI_MODE_WRONLY, io_info, handle, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error ('input_snap', &
          'Could not open: "'//trim (directory_snap)//'/'//file//'"')
!
! Setting file view and write raw binary data, ie. 'native'.
!
      call MPI_FILE_SET_VIEW(handle, 0, mpi_precision, global_type, &
          'native', io_info, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error ('input_snap', &
          'Could not create view for "'//file//'"')
!
      call MPI_FILE_READ_ALL(handle, a, 1, local_type, status, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error ('input_snap', &
          'Could not write "'//file//'"')
!
      call MPI_FILE_CLOSE(handle, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error ('input_snap', &
          'Could not close "'//file//'"')
!
      call MPI_FINALIZE(mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error ('input_snap', &
          'Could not finalize "'//file//'"')
!
      ! read additional data
      if (lread_add) then
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('input_snap', 'Could not allocate memory for gx,gy,gz', .true.)
!
          open (lun_input, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', status='old',position='append')
          call backskip_to_time(lun_input)
          read (lun_input) t_sp, gx, gy, gz, dx, dy, dz
          call distribute_grid (x, y, z, gx, gy, gz)
          deallocate (gx, gy, gz)
        else
          call distribute_grid (x, y, z)
        endif
        call mpibcast_real (t_sp)
        t = t_sp
      endif
!
    endsubroutine input_snap
!***********************************************************************
    subroutine input_snap_finalize()
!
!  Close snapshot file.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      if (persist_initialized) then
        close (lun_input)
        persist_initialized = .false.
        persist_last_id = -max_int
      endif
!
    endsubroutine input_snap_finalize
!***********************************************************************
    logical function init_write_persist(file)
!
!  Initialize writing of persistent data to persistent file.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
        open (lun_output, FILE=trim (directory_dist)//'/'//filename, FORM='unformatted', status='replace')
        if (ip <= 9) write (*,*) 'begin persistent block'
        write (lun_output) id_block_PERSISTENT
        filename = ""
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      use Mpicomm, only: mpibcast_logical
      use Syscalls, only: file_exists
!
      character (len=*), intent(in), optional :: file
!
      init_read_persist = .true.
!
      if (present (file)) then
        if (lroot) init_read_persist = .not. file_exists (trim (directory_snap)//'/'//file)
        call mpibcast_logical (init_read_persist)
        if (init_read_persist) return
      endif
!
      if (present (file)) then
        if (lroot .and. (ip <= 9)) write (*,*) 'begin persistent block'
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
      call mpibcast_int (id)
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
      call output_snap (a, nv, file, 0)
      call output_snap_finalize ()
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
      call input_snap (file, a, nv, 0)
      call input_snap_finalize ()
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
        call mpibarrier ()
        if (lroot) then
          open (lun_output,FILE=trim (datadir)//'/move-me.list', POSITION='append')
          write (lun_output,'(A)') trim (fpart)
          close (lun_output)
        endif
      endif
!
    endsubroutine log_filename_to_file
!***********************************************************************
    subroutine wgrid(file)
!
!  Write grid coordinates.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      character (len=*) :: file
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (lroot) then
        allocate (gx(nxgrid+2*nghost), gy(nygrid+2*nghost), gz(nzgrid+2*nghost), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('wgrid', 'Could not allocate memory for gx,gy,gz', .true.)
!
        open (lun_output, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', status='replace')
        t_sp = t
        call collect_grid (x, y, z, gx, gy, gz)
        write (lun_output) t_sp, gx, gy, gz, dx, dy, dz
        write (lun_output) dx, dy, dz
        write (lun_output) Lx, Ly, Lz
        call collect_grid (dx_1, dy_1, dz_1, gx, gy, gz)
        write (lun_output) gx, gy, gz
        call collect_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
        write (lun_output) gx, gy, gz
        close (lun_output)
!
        deallocate (gx, gy, gz)
      else
        call collect_grid (x, y, z)
        call collect_grid (dx_1, dy_1, dz_1)
        call collect_grid (dx_tilde, dy_tilde, dz_tilde)
      endif
!
    endsubroutine wgrid
!***********************************************************************
    subroutine rgrid(file)
!
!  Read grid coordinates.
!
!  19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
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
      call mpibcast_real (dx)
      call mpibcast_real (dy)
      call mpibcast_real (dz)
      call mpibcast_real (Lx)
      call mpibcast_real (Ly)
      call mpibcast_real (Lz)
!
!  Find minimum/maximum grid spacing. Note that
!    minval( (/dx,dy,dz/), MASK=((/nxgrid,nygrid,nzgrid/) > 1) )
!  will be undefined if all n[xyz]grid==1, so we have to add the fourth
!  component with a test that is always true
!
      dxmin = minval ((/ dx, dy, dz,    huge (dx) /), MASK=((/ nxgrid, nygrid, nzgrid, 2 /) > 1))
      dxmax = maxval ((/ dx, dy, dz, epsilon (dx) /), MASK=((/ nxgrid, nygrid, nzgrid, 2 /) > 1))
!
!  Fill pencil with maximum gridspacing. Will be overwritten
!  during the mn loop in the non equiditant case
!
      dxmax_pencil(:) = dxmax
!
      if (lroot) then
        if (ip <= 4) then
          print *, 'rgrid: Lx,Ly,Lz=', Lx, Ly, Lz
          print *, 'rgrid: dx,dy,dz=', dx, dy, dz
          print *, 'rgrid: dxmin,dxmax=', dxmin, dxmax
        endif
        if (dxmin == 0) call fatal_error ("rgrid", "check Lx,Ly,Lz: is one of them 0?", .true.)
      endif
!
    endsubroutine rgrid
!***********************************************************************
    subroutine wproc_bounds(file)
!
!   Export processor boundaries to file.
!
!   19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
!
      integer :: ierr
!
      open (lun_output, FILE=file, FORM='unformatted', IOSTAT=ierr, status='replace')
      if (ierr /= 0) call stop_it ( &
          "Cannot open " // trim(file) // " (or similar) for writing" // &
          " -- is data/ visible from all nodes?")
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
!   19-Sep-2012/Bourdin.KIS: adapted from io_mpi2
!
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
!
      integer :: ierr
!
      open (lun_input, FILE=file, FORM='unformatted', IOSTAT=ierr, status='old')
      if (ierr /= 0) call stop_it ( &
          "Cannot open " // trim(file) // " (or similar) for reading" // &
          " -- is data/ visible from all nodes?")
      read (lun_input) procy_bounds
      read (lun_input) procz_bounds
      close (lun_input)
!
    endsubroutine rproc_bounds
!***********************************************************************
endmodule Io
