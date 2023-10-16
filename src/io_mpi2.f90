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
  use Cparam, only: fnlen, max_int
  use Messages, only: fatal_error, fatal_error_local, fatal_error_local_collect, svn_id, warning
  use Mpicomm, only: mpi_precision
!
  implicit none
!
  include "mpif.h"
  include 'io.h'
  include 'record_types.h'
!
! Indicates if IO is done distributed (each proc writes into a procdir)
! or collectively (eg. by specialized IO-nodes or by MPI-IO).
!
  logical :: lcollective_IO=.true.
  character (len=labellen) :: IO_strategy="MPI-IO"
!
  logical :: lread_add=.true., lwrite_add=.true.
  integer :: persist_last_id=-max_int
!
  integer :: local_type, global_type, mpi_err, io_dims
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer (kind=MPI_OFFSET_KIND), parameter :: displacement=0
  integer, parameter :: order=MPI_ORDER_FORTRAN, io_info=MPI_INFO_NULL
!  integer, dimension(io_dims) :: local_size, local_start, global_size, global_start, subsize
  integer, dimension (:), allocatable :: local_size, local_start, global_size, global_start, subsize
!
  interface input_proc_bounds
    module procedure input_proc_bounds_double
    module procedure input_proc_bounds_single
  endinterface
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
      if (ldistribute_persist) call fatal_error ('io_mpi2', "Distributed persistent variables not allowed with MPI-IO.")
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
!  04-Sep-2015/PABourdin: coded
!
      use Mpicomm, only: mpibcast_real
!
      real, dimension(mx), intent(out) :: x
      real, dimension(my), intent(out) :: y
      real, dimension(mz), intent(out) :: z
      real, dimension(mxgrid), intent(in), optional :: gx
      real, dimension(mygrid), intent(in), optional :: gy
      real, dimension(mzgrid), intent(in), optional :: gz
!
      real, dimension(mxgrid+mygrid+mzgrid) :: tmp_grid
      integer, parameter :: tag_gx=680, tag_gy=681, tag_gz=682
!
      if (lroot) then
        tmp_grid(1:mxgrid) = gx
        tmp_grid(mxgrid+1:mxgrid+mygrid) = gy
        tmp_grid(mxgrid+mygrid+1:mxgrid+mygrid+mzgrid) = gz
      endif
      call mpibcast_real (tmp_grid, mxgrid+mygrid+mzgrid)
      x = tmp_grid(ipx*nx+1:ipx*nx+mx)
      y = tmp_grid(mxgrid+ipy*ny+1:mxgrid+ipy*ny+my)
      z = tmp_grid(mxgrid+mygrid+ipz*nz+1:mxgrid+mygrid+ipz*nz+mz)
!
    endsubroutine distribute_grid
!***********************************************************************
    subroutine check_consistency(x, mask)
!
!  Check if a scalar value is consistent between active processes.
!
!  15-nov-20/ccyang: coded
!
!  Input Arguments
!      xlocal
!          Local value.
!  Optional Input Arguments
!      mask
!          If present, processes to be included in the check.
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, intent(in) :: x
      logical, dimension(ncpus), intent(in), optional :: mask
!
      real, dimension(ncpus) :: buf
      real :: xmin, xmax
!
!  Collect the local values.
!
      buf = 0.0
      buf(iproc+1) = x
      call mpiallreduce_sum(buf, ncpus)
!
!  Find the range of the values.
!
      val: if (present(mask)) then
        if (.not. any(mask)) return
        xmin = minval(buf, mask=mask)
        xmax = maxval(buf, mask=mask)
      else val
        xmin = minval(buf)
        xmax = maxval(buf)
      endif val
!
!  Determine if the local values are consistent.
!
      cons: if (xmin /= xmax) then
        print *, "check_consistency: xmin, xmax = ", xmin, xmax
        call fatal_error("check_consistency", "local values are not consistent")
      endif cons
!
    endsubroutine check_consistency
!***********************************************************************
    subroutine check_success(routine, message, file)
!
!  Check success of MPI2 file access and issue error if necessary.
!
!  03-Oct-2015/PABourdin: coded
!
      character (len=*), intent(in) :: routine, message, file
!
      if (mpi_err /= MPI_SUCCESS) call fatal_error (trim(routine)//'_snap', 'Could not '//trim(message)//' "'//file//'"')
!
    endsubroutine check_success
!***********************************************************************
    subroutine check_success_local(routine, message)
!
!  Check success locally of MPI2 file access and issue error if necessary.
!
!  11-nov-20/ccyang: coded
!
      character(len=*), intent(in) :: routine, message
!
      if (mpi_err /= MPI_SUCCESS) call fatal_error_local(trim(routine) // "_snap", "could not " // trim(message))
!
    endsubroutine check_success_local
!***********************************************************************
    subroutine get_dimensions(dir, nxyzgrid, nxyz, ipxyz)
!
!  Get n[xyz]grid, n[xyz], and ip[xyz] according to dir.
!
!  14-nov-20/ccyang: coded
!
      character, intent(in) :: dir
      integer, intent(out) :: nxyzgrid, nxyz, ipxyz
!
      chdir: if (dir == 'x') then
        nxyzgrid = nxgrid
        nxyz = nx
        ipxyz = ipx
!
      elseif (dir == 'y') then chdir
        nxyzgrid = nygrid
        nxyz = ny
        ipxyz = ipy
!
      elseif (dir == 'z') then chdir
        nxyzgrid = nzgrid
        nxyz = nz
        ipxyz = ipz
!
      else chdir
        call fatal_error("get_dimensions", "unknown direction " // dir)
!
      endif chdir
!
    endsubroutine get_dimensions
!***********************************************************************
    pure integer(KIND=MPI_OFFSET_KIND) function get_disp_to_par_real(npar_tot) result(disp)
!
!  Gets the displacement in bytes to the beginning of particle real
!  data.
!
!  12-nov-20/ccyang: coded
!
      use Mpicomm, only: size_of_int
!
      integer, intent(in) :: npar_tot
!
      disp = int(1 + npar_tot, KIND=MPI_OFFSET_KIND) * size_of_int
!
    endfunction get_disp_to_par_real
!***********************************************************************
    subroutine output_snap(a, nv1, nv2, file, mode)
!
!  Write snapshot file, always write mesh and time, could add other things.
!
!  10-Feb-2012/PABourdin: coded
!  13-feb-2014/MR: made file optional (prep for downsampled output)
!
      use General, only: ioptest
      use Mpicomm, only: globalize_xy, collect_grid
!
      integer, optional, intent(in) :: nv1,nv2
      real, dimension (:,:,:,:), intent(in) :: a
      character (len=*), optional, intent(in) :: file
      integer, optional, intent(in) :: mode
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer(kind=MPI_OFFSET_KIND) total_size
      integer :: handle, alloc_err, na, ne, nv
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (.not. present (file)) call fatal_error ('output_snap', 'downsampled output not implemented for IO_mpi2')
!
      lwrite_add = .true.
      if (present (mode)) lwrite_add = (mode == 1)
!
      na=ioptest(nv1,1)
      ne=ioptest(nv2,mvar_io)
      nv=ne-na+1
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
      call MPI_FILE_OPEN (MPI_COMM_WORLD, trim (directory_snap)//'/'//file, MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                          io_info, handle, mpi_err)
      call check_success ('output', 'open', trim (directory_snap)//'/'//file)
!
! Truncate the file to fit the global data.
!
      total_size = product(int(global_size, KIND=MPI_OFFSET_KIND))
      call MPI_FILE_SET_SIZE(handle, total_size, mpi_err)
      call check_success ("output", "set size of", file)
!
! Setting file view and write raw binary data, ie. 'native'.
!
      call MPI_FILE_SET_VIEW (handle, displacement, mpi_precision, global_type, 'native', io_info, mpi_err)
      call check_success ('output', 'create view', file)
!
      if (lwrite_2D) then
        if (ny == 1) then
          call MPI_FILE_WRITE_ALL (handle, a(:,m1,:,na:ne), 1, local_type, status, mpi_err)
        else
          call MPI_FILE_WRITE_ALL (handle, a(:,:,n1,na:ne), 1, local_type, status, mpi_err)
        endif
      else
        call MPI_FILE_WRITE_ALL (handle, a(:,:,:,na:ne), 1, local_type, status, mpi_err)
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
          if (lshear) then
            write (lun_output) t_sp, gx, gy, gz, dx, dy, dz, deltay
          else
            write (lun_output) t_sp, gx, gy, gz, dx, dy, dz
          endif
        endif

        if (lode) call output_ode(file)

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
    subroutine output_average_2D(label, navg, avgname, avgdata, time, lwrite, header)
!
!  Output 2D average to a file.
!
!  18-dec-2022/ccyang: coded
!
      use Cparam, only: root
      use General, only: keep_compiler_quiet
      use Mpicomm, only: size_of_real
!
      integer, intent(in) :: navg
      character(len=fmtlen), dimension(navg), intent(in) :: avgname
      character(len=*), intent(in) :: label
      real, dimension(:,:,:), intent(in) :: avgdata
      real, intent(in) :: time
      logical, intent(in) :: lwrite
      real, dimension(:), optional, intent(in) :: header
!
      character(len=*), parameter :: rname = "output_average_2D"
      integer, dimension(2) :: asizes, asubs, astarts
      character(len=fnlen) :: fpath
      integer(KIND=MPI_OFFSET_KIND) :: fsize, dsize, disp
      integer :: handle, dtype
      integer :: i, n
      real :: tread
!
      call keep_compiler_quiet(avgname)
      if (present(header)) call keep_compiler_quiet(header)
!
!  Find the data size.
!
      dims: if (label == "x") then
        asizes = (/ nygrid, nzgrid /)
        asubs = (/ ny, nz /)
        astarts = (/ ipy * ny, ipz * nz /)
      elseif (label == "y") then dims
        asizes = (/ nxgrid, nzgrid /)
        asubs = (/ nx, nz /)
        astarts = (/ ipx * nx, ipz * nz /)
      elseif (label == "z") then dims
        asizes = (/ nxgrid, nygrid /)
        asubs = (/ nx, ny /)
        astarts = (/ ipx * nx, ipy * ny /)
      else dims
        call fatal_error(rname, "unknown label " // trim(label))
      endif dims
!
      dsize = int(navg * product(asizes) + 1, kind=MPI_OFFSET_KIND) * size_of_real
!
!  Open average file.
!
      fpath = trim(directory_snap) // '/' // trim(label) // "averages.dat"
      call MPI_FILE_OPEN(MPI_COMM_WORLD, fpath, ior(MPI_MODE_CREATE, MPI_MODE_RDWR), io_info, handle, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error(rname, "unable to open file '" // trim(fpath) // "'")
!
!  Get the file size.
!
      call MPI_FILE_GET_SIZE(handle, fsize, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error(rname, "unable to get file size")
!
      ckfsize: if (mod(fsize, dsize) /= 0) then
        if (lroot) print *, rname, ": fsize, dsize = ", fsize, dsize
        call fatal_error(rname, "file size is not a multiple of data size. ")
      endif ckfsize
!
!  Back scan.
!
      bscan: if (lroot) then
        disp = fsize
        tread = time
        step: do while (tread >= time .and. disp > 0_MPI_OFFSET_KIND)
          disp = disp - dsize
          call MPI_FILE_READ_AT(handle, disp, tread, 1, mpi_precision, status, mpi_err)
          if (mpi_err /= MPI_SUCCESS) call fatal_error_local(rname, "unable to read time")
        enddo step
        if (tread < time) disp = disp + dsize
      endif bscan
      call fatal_error_local_collect()
!
!  Truncate file if needed.
!
      call MPI_BCAST(disp, 1, MPI_OFFSET, root, MPI_COMM_WORLD, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error(rname, "unable to broadcast displacement")
!
      trunc: if (disp < fsize) then
        n = int((fsize - disp) / dsize)
        fsize = disp
        call MPI_FILE_SET_SIZE(handle, fsize, mpi_err)
        if (mpi_err /= MPI_SUCCESS) call fatal_error(rname, "unable to truncate file. ")
        if (lroot) print *, rname, ": trunacted ", n, " chunks of 2D averages. "
      endif trunc
!
!  Write time.
!
      wtime: if (lroot) then
        call MPI_FILE_SEEK(handle, 0_MPI_OFFSET_KIND, MPI_SEEK_END, mpi_err)
        if (mpi_err /= MPI_SUCCESS) call fatal_error_local(rname, "unable to move handle")
        call MPI_FILE_WRITE(handle, time, 1, mpi_precision, status, mpi_err)
        if (mpi_err /= MPI_SUCCESS) call fatal_error_local(rname, "unable to write time")
      endif wtime
      call fatal_error_local_collect()
      disp = fsize + size_of_real
!
!  Decompose the write by processes.
!
      call MPI_TYPE_CREATE_SUBARRAY(2, asizes, asubs, astarts, order, mpi_precision, dtype, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error(rname, "unable to create subarray")
!
      call MPI_TYPE_COMMIT(dtype, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error(rname, "unable to commit data type")
!
      call MPI_FILE_SET_VIEW(handle, disp, mpi_precision, dtype, "native", io_info, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error(rname, "unable to set local view")
!
!  Write the averages.
!
      n = merge(product(asubs), 0, lwrite)
      comp: do i = 1, navg
        if (label == "z") then
          call MPI_FILE_WRITE_ALL(handle, avgdata(i,:,:), n, mpi_precision, status, mpi_err)
        else
          call MPI_FILE_WRITE_ALL(handle, avgdata(:,:,i), n, mpi_precision, status, mpi_err)
        endif
        if (mpi_err /= MPI_SUCCESS) call fatal_error(rname, "unable to write average")
      enddo comp
!
!  Clean up and close average file.
!
      call MPI_TYPE_FREE(dtype, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error(rname, "unable to free data type")
!
      call MPI_FILE_CLOSE(handle, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error(rname, "unable to close file")
!
    endsubroutine output_average_2D
!***********************************************************************
    subroutine output_slice_position()
!
!  Dummy subroutine; no need to record slice positions.
!
!  14-nov-20/ccyang: dummy
!
    endsubroutine output_slice_position
!***********************************************************************
    subroutine output_slice(lwrite, time, label, suffix, pos, grid_pos, data)
!
!  Append to a slice file.
!
!  15-nov-20/ccyang: coded
!
      use General, only: keep_compiler_quiet
      use Mpicomm, only: mpiallreduce_or, size_of_real
!
      real, dimension(:,:), pointer :: data
      character(len=*), intent(in) :: label, suffix
      logical, intent(in) :: lwrite
      integer, intent(in) :: grid_pos
      real, intent(in) :: time, pos
!
      integer, parameter :: nadd = 2  ! number of additional real data
!
      logical, dimension(ncpus) :: lwrite_proc
      integer, dimension(2) :: sizes, subsizes, starts
      character(len=fnlen) :: fpath
      integer(KIND=MPI_OFFSET_KIND) :: dsize, fsize, offset
      integer :: handle, dtype, k
      real :: tprev, tcut
!
      call keep_compiler_quiet(grid_pos)
!
!  Communicate the lwrite status.
!
      lwrite_proc = .false.
      lwrite_proc(iproc+1) = lwrite
      call mpiallreduce_or(lwrite_proc, ncpus)
      if (.not. any(lwrite_proc)) return
!
!  Check the consistency of time and position.
!
      call check_consistency(time, mask=lwrite_proc)
      call check_consistency(pos, mask=lwrite_proc)
!
!  Define the data structure.
!
      do k = 1, 2
        call get_dimensions(suffix(k:k), sizes(k), subsizes(k), starts(k))
      enddo
      starts = starts * subsizes
!
      call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, order, mpi_precision, dtype, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error_local("output_slice", "cannot create subarray type")
      call fatal_error_local_collect()
!
      dsize = int(product(sizes), KIND=MPI_OFFSET_KIND) * size_of_real
      call MPI_TYPE_CREATE_STRUCT(2, (/ 1, nadd /), (/ 0_MPI_OFFSET_KIND, dsize /), (/ dtype, mpi_precision /), dtype, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error_local("output_slice", "cannot create struct type")
      call fatal_error_local_collect()
      dsize = dsize + int(nadd, KIND=MPI_OFFSET_KIND) * size_of_real
!
!  Open the slice file for write.
!
      fpath = trim(directory_snap) // "/slice_" // trim(label) // '.' // trim(suffix)
      call MPI_FILE_OPEN(MPI_COMM_WORLD, fpath, ior(MPI_MODE_APPEND, ior(MPI_MODE_CREATE, MPI_MODE_RDWR)), io_info, handle, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("output_slice", "cannot open file '" // trim(fpath) // "'")
!
!  Check the file size.
!
      call MPI_FILE_GET_SIZE(handle, fsize, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error_local("output_slice", "cannot get file size")
      call fatal_error_local_collect()
!
      ckfsize: if (mod(fsize, dsize) /= 0) then
        if (lroot) print *, "output_slice: fsize, dsize = ", fsize, dsize
        call fatal_error("output_slice", "The file size is not a multiple of the data size. ")
      endif ckfsize
!
!  Back scan the existing slices.
!
      master: if (lroot) then
        tcut = dvid * real(nint(time / dvid))
        offset = fsize - int(nadd, KIND=MPI_OFFSET_KIND) * size_of_real
        bscan: do while (offset > 0_MPI_OFFSET_KIND)
          call MPI_FILE_READ_AT(handle, offset, tprev, 1, mpi_precision, status, mpi_err)
          if (mpi_err /= MPI_SUCCESS) call fatal_error_local("output_slice", "cannot read time")
          if (tprev < tcut) exit
          offset = offset - dsize
        enddo bscan
        offset = offset + int(nadd, KIND=MPI_OFFSET_KIND) * size_of_real
      endif master
      call fatal_error_local_collect()
      call MPI_BCAST(offset, 1, MPI_OFFSET, root, MPI_COMM_WORLD, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("output_slice", "unable to broadcast offset")
!
!  Truncate the slices with later times.
!
      trunc: if (fsize > offset) then
        k = int((fsize - offset) / dsize)
        fsize = offset
        call MPI_FILE_SET_SIZE(handle, fsize, mpi_err)
        if (mpi_err /= MPI_SUCCESS) call fatal_error("output_slice", "cannot set file size")
        if (lroot) print *, "output_slice: truncated from ", trim(fpath), k, " slices with t >= ", tcut
      endif trunc
!
!  Write the slice.
!
      call MPI_TYPE_COMMIT(dtype, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error_local("output_slice", "cannot commit data type")
      call fatal_error_local_collect()
!
      call MPI_FILE_SET_VIEW(handle, fsize, mpi_precision, dtype, "native", io_info, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("output_slice", "cannot set local view")
!
      if (lwrite) then
        call MPI_FILE_WRITE_ALL(handle, data, product(subsizes), mpi_precision, status, mpi_err)
      else
        call MPI_FILE_WRITE_ALL(handle, huge(0.0), 0, mpi_precision, status, mpi_err)
      endif
      if (mpi_err /= MPI_SUCCESS) call fatal_error("output_slice", "cannot write data")
!
      call MPI_TYPE_FREE(dtype, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error_local("output_slice", "cannot free MPI data type")
      call fatal_error_local_collect()
!
!  Write time and slice position.
!
      if (lwrite) then
        call MPI_FILE_WRITE_ALL(handle, (/ time, pos /), nadd, mpi_precision, status, mpi_err)
      else
        call MPI_FILE_WRITE_ALL(handle, huge(0.0), 0, mpi_precision, status, mpi_err)
      endif
      if (mpi_err /= MPI_SUCCESS) call fatal_error("output_slice", "cannot write additional data")
!
!  Close the file.
!
      call MPI_FILE_CLOSE(handle, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("output_slice", "cannot close file '" // trim(fpath) // "'")
!
    endsubroutine output_slice
!***********************************************************************
    subroutine output_part_snap(ipar, a, mv, nv, file, label, ltruncate)
!
!  Write particle snapshot file.
!
!  23-Oct-2018/PABourdin: stub
!  12-nov-2020/ccyang: coded
!
!  NOTE: The optional argument ltruncate is required by IO=io_hdf5.
!
      use General, only: keep_compiler_quiet
      use Mpicomm, only: mpiallreduce_sum_int, size_of_int, size_of_real
!
      integer, intent(in) :: mv, nv
      integer, dimension(mv), intent(in) :: ipar
      real, dimension(mv,mparray), intent(in) :: a
      character(len=*), intent(in) :: file
      character(len=*), intent(in), optional :: label
      logical, intent(in), optional :: ltruncate
!
      integer, dimension(ncpus) :: npar_proc
      character(len=fnlen) :: fpath
      integer :: handle, ftype, npar_tot, ip0, nv1
      integer(KIND=MPI_OFFSET_KIND) :: offset
!
      if (present(ltruncate)) call keep_compiler_quiet(ltruncate)
      if (present(label)) call warning("output_part_snap", "The argument label has no effects. ")
!
!  Broadcast number of particles from each process.
!
      npar_proc = 0
      npar_proc(iproc+1) = nv
      call mpiallreduce_sum_int(npar_proc, ncpus)
      npar_tot = sum(npar_proc)
!
!  Open snapshot file for write.
!
      fpath = trim(directory_snap) // '/' // file
      call MPI_FILE_OPEN(MPI_COMM_WORLD, fpath, ior(MPI_MODE_CREATE, MPI_MODE_WRONLY), io_info, handle, mpi_err)
      call check_success("output_part", "open", fpath)
!
!  Write total number of particles.
!
      wnpar: if (lroot) then
        call MPI_FILE_WRITE(handle, npar_tot, 1, MPI_INTEGER, status, mpi_err)
        call check_success_local("output_part", "write number of particles")
      endif wnpar
      call fatal_error_local_collect()
!
!  Write particle IDs.
!
      ip0 = sum(npar_proc(:iproc))
      nv1 = max(nv, 1)
      call MPI_TYPE_CREATE_SUBARRAY(1, (/ npar_tot /), (/ nv1 /), (/ ip0 /), order, MPI_INTEGER, ftype, mpi_err)
      call check_success_local("output_part", "create MPI subarray")
      call fatal_error_local_collect()
!
      call MPI_TYPE_COMMIT(ftype, mpi_err)
      call check_success_local("output_part", "commit MPI data type")
      call fatal_error_local_collect()
!
      call MPI_FILE_SET_VIEW(handle, size_of_int, MPI_INTEGER, ftype, "native", io_info, mpi_err)
      call check_success("output_part", "set view of", fpath)
!
      call MPI_FILE_WRITE_ALL(handle, ipar, nv, MPI_INTEGER, status, mpi_err)
      call check_success("output_part", "write particle IDs", fpath)
!
      call MPI_TYPE_FREE(ftype, mpi_err)
      call check_success_local("output_part", "free MPI data type")
      call fatal_error_local_collect()
!
!  Write particle data.
!
      call MPI_TYPE_CREATE_SUBARRAY(2, (/ npar_tot, mparray /), (/ nv1, mparray /), (/ ip0, 0 /), &
                                    order, mpi_precision, ftype, mpi_err)
      call check_success_local("output_part", "create MPI subarray")
      call fatal_error_local_collect()
!
      call MPI_TYPE_COMMIT(ftype, mpi_err)
      call check_success_local("output_part", "commit MPI data type")
      call fatal_error_local_collect()
!
      offset = get_disp_to_par_real(npar_tot)
      call MPI_FILE_SET_VIEW(handle, offset, mpi_precision, ftype, "native", io_info, mpi_err)
      call check_success("output_part", "set global view of", fpath)
!
      call MPI_FILE_WRITE_ALL(handle, a(1:nv,:), nv * mparray, mpi_precision, status, mpi_err)
      call check_success("output_part", "write particle data at", fpath)
!
      call MPI_TYPE_FREE(ftype, mpi_err)
      call check_success_local("output_part", "free subarray type")
      call fatal_error_local_collect()
!
!  Write additional data.
!
      offset = offset + int(npar_tot * mparray, KIND=MPI_OFFSET_KIND) * size_of_real
      call MPI_FILE_SET_VIEW(handle, offset, mpi_precision, mpi_precision, "native", io_info, mpi_err)
      call check_success("output_part", "set global view of", fpath)
      if (lroot) then
        call MPI_FILE_WRITE_ALL(handle, real(t), 1, mpi_precision, status, mpi_err)
      else
        call MPI_FILE_WRITE_ALL(handle, huge(0.0), 0, mpi_precision, status, mpi_err)
      endif
      call check_success("output_part", "write additional data", fpath)
!
!  Close snapshot file.
!
      call MPI_FILE_CLOSE(handle, mpi_err)
      call check_success("output_part", "close", fpath)
!
    endsubroutine output_part_snap
!***********************************************************************
    subroutine output_stalker_init(num, nv, snap, ID)
!
!  Open stalker particle snapshot file and initialize with snapshot time.
!
!  03-May-2019/PABourdin: coded
!
      use General, only: keep_compiler_quiet
!
      integer, intent(in) :: num, nv, snap
      integer, dimension(nv), intent(in) :: ID
!
      call fatal_error ('output_stalker_init', 'not implemented for "io_mpi2"', .true.)
      call keep_compiler_quiet(num, nv, snap)
      call keep_compiler_quiet(ID)
!
    endsubroutine output_stalker_init
!***********************************************************************
    subroutine output_stalker(label, mv, nv, data, nvar, lfinalize)
!
!  Write stalker particle quantity to snapshot file.
!
!  03-May-2019/PABourdin: coded
!
      use General, only: keep_compiler_quiet
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: mv, nv
      real, dimension (mv), intent(in) :: data
      logical, intent(in), optional :: lfinalize
      integer, intent(in), optional :: nvar
!
      call fatal_error ('output_stalker', 'not implemented for "io_mpi2"', .true.)
      call keep_compiler_quiet(label)
      call keep_compiler_quiet(data)
      call keep_compiler_quiet(mv, nv)
      if (present(lfinalize)) call warning("output_stalker", "The argument lfinalize has no effects. ")
      if (present(nvar)) call warning("output_stalker", "The argument nvar has no effects. ")
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
      use General, only: keep_compiler_quiet
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mv, nc
      character (len=*), dimension (mqarray), intent(in) :: labels
      real, dimension (mv,mparray), intent(in) :: fq
!
      call fatal_error ('output_pointmass', 'not implemented for "io_mpi2"', .true.)
      call keep_compiler_quiet(file)
      call keep_compiler_quiet(labels)
      call keep_compiler_quiet(fq)
      call keep_compiler_quiet(mv, nc)
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
      use Mpicomm, only: localize_xy, mpibcast_real
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
      call MPI_FILE_OPEN (MPI_COMM_WORLD, trim (directory_snap)//'/'//file, MPI_MODE_RDONLY, io_info, handle, mpi_err)
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
!
        if (lode) call input_ode(file)

      endif

    endsubroutine input_snap
!***********************************************************************
    subroutine input_snap_finalize
!
!  Close snapshot file.
!
!  11-Feb-2012/PABourdin: coded
!
      use Messages, only: not_implemented
!
      if (persist_initialized) then
        if (ldistribute_persist .or. lroot) close (lun_input)
        persist_initialized = .false.
        persist_last_id = -max_int
      elseif (lread_add .and. lroot) then
        close (lun_input)
      endif
!
      if (snaplink /= "") call not_implemented("input_snap_finalize", "module variable snaplink")
!
    endsubroutine input_snap_finalize
!***********************************************************************
    subroutine input_slice_real_arr(file, time, pos, data)
!
!  Read a slice file.
!
!  14-nov-20/ccyang: stub
!
      use General, only: keep_compiler_quiet
      use Messages, only: not_implemented
!
      real, dimension(:,:,:), intent(out):: data
      character(len=*), intent(in) :: file
      real, intent(out):: time, pos
!
      call not_implemented("output_slice_position")
      call keep_compiler_quiet(data)
      call keep_compiler_quiet(file)
      call keep_compiler_quiet(time, pos)
!
    endsubroutine input_slice_real_arr
!***********************************************************************
    subroutine input_slice_scat_arr(file, pos, data, ivar, nt)
!
!  Read a slice file.
!
!  14-nov-20/ccyang: stub
!
      use General, only: keep_compiler_quiet, scattered_array
      use Messages, only: not_implemented
!
      type(scattered_array), pointer :: data   !intent(inout)
      character(len=*), intent(in) :: file
      integer, intent(in) :: ivar
      integer, intent(out):: nt
      real, intent(out):: pos
!
      call not_implemented("output_slice_position")
      call keep_compiler_quiet(data%dim1)
      call keep_compiler_quiet(file)
      call keep_compiler_quiet(ivar, nt)
      call keep_compiler_quiet(pos)
!
    endsubroutine input_slice_scat_arr
!***********************************************************************
    subroutine input_part_snap(ipar, ap, mv, nv, npar_tot, file, label)
!
!  Read particle snapshot file.
!
!  24-Oct-2018/PABourdin: stub
!  12-nov-20/ccyang: coded
!
      use General, only: keep_compiler_quiet
      use Particles_cdata, only: ixp, iyp, izp
      use Mpicomm, only: size_of_int, size_of_real
!
      integer, intent(in) :: mv
      integer, dimension(mv), intent(out) :: ipar
      real, dimension(mv,mparray), intent(out) :: ap
      integer, intent(out) :: nv, npar_tot
      character(len=*), intent(in) :: file
      character(len=*), intent(in), optional :: label
!
      logical, allocatable, dimension(:) :: lpar_loc
      integer, allocatable, dimension(:) :: indices
!
      real, dimension(:), allocatable :: rbuf
      character(len=fnlen) :: fpath
      integer(KIND=MPI_OFFSET_KIND) :: offset, dsize
      integer :: istat, handle, ftype, i, k
!
      if (present(label)) call warning("input_part_snap", "The argument label has no effects. ")
!
!  Open snapshot file for read.
!
      fpath = trim(directory_snap) // '/' // file
      call MPI_FILE_OPEN(MPI_COMM_WORLD, fpath, MPI_MODE_RDONLY, io_info, handle, mpi_err)
      call check_success("input_part", "open", fpath)
!
!  Read total number of particles.
!
      nptot: if (lroot) then
        call MPI_FILE_READ(handle, npar_tot, 1, MPI_INTEGER, status, mpi_err)
        call check_success_local("input_part", "read total number of particles")
      endif nptot
      call fatal_error_local_collect()
      call MPI_BCAST(npar_tot, 1, MPI_INTEGER, root, MPI_COMM_WORLD, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("input_part", "unable to broadcast total number of particles")
!
      cknp: if (ceiling(real(npar_tot) / real(ncpus)) > mv) then
        if (lroot) print *, "input_part_snap: npar_tot, mv = ", npar_tot, mv
        call fatal_error("input_part_snap", "too many particles")
      endif cknp
!
      allocate(rbuf(npar_tot), stat=istat)
      if (istat /= 0) call fatal_error_local("input_part_snap", "cannot allocate buffer")
      call fatal_error_local_collect()
!
!  Identify local particles.
!
      allocate(lpar_loc(npar_tot), stat=mpi_err)
      call check_success_local("input_part", "allocate lpar_loc")
      call fatal_error_local_collect()
      lpar_loc = .true.
!
      offset = get_disp_to_par_real(npar_tot)
      call MPI_FILE_SET_VIEW(handle, offset + (ixp - 1) * int(npar_tot, KIND=MPI_OFFSET_KIND), &
                             mpi_precision, mpi_precision, "native", io_info, mpi_err)
      call check_success("input_part", "set view of", fpath)
!
      inx: if (lactive_dimension(1)) then
        xp: if (lroot) then
          call MPI_FILE_READ(handle, rbuf, npar_tot, mpi_precision, status, mpi_err)
          call check_success_local("input_part", "read xp of")
        endif xp
        call fatal_error_local_collect()
        call MPI_BCAST(rbuf, npar_tot, mpi_precision, root, MPI_COMM_WORLD, mpi_err)
        if (mpi_err /= MPI_SUCCESS) call fatal_error("input_part", "unable to broadcast xp. ")
        lpar_loc = lpar_loc .and. procx_bounds(ipx) <= rbuf(1:npar_tot) .and. rbuf(1:npar_tot) < procx_bounds(ipx+1)
      endif inx
!
      iny: if (lactive_dimension(2)) then
        yp: if (lroot) then
          call MPI_FILE_READ(handle, rbuf, npar_tot, mpi_precision, status, mpi_err)
          call check_success_local("input_part", "read yp of")
        endif yp
        call fatal_error_local_collect()
        call MPI_BCAST(rbuf, npar_tot, mpi_precision, root, MPI_COMM_WORLD, mpi_err)
        if (mpi_err /= MPI_SUCCESS) call fatal_error("input_part", "unable to broadcast yp. ")
        lpar_loc = lpar_loc .and. procy_bounds(ipy) <= rbuf(1:npar_tot) .and. rbuf(1:npar_tot) < procy_bounds(ipy+1)
      endif iny
!
      inz: if (lactive_dimension(3)) then
        zp: if (lroot) then
          call MPI_FILE_READ(handle, rbuf, npar_tot, mpi_precision, status, mpi_err)
          call check_success_local("input_part", "read zp of")
        endif zp
        call fatal_error_local_collect()
        call MPI_BCAST(rbuf, npar_tot, mpi_precision, root, MPI_COMM_WORLD, mpi_err)
        if (mpi_err /= MPI_SUCCESS) call fatal_error("input_part", "unable to broadcast zp. ")
        lpar_loc = lpar_loc .and. procz_bounds(ipz) <= rbuf(1:npar_tot) .and. rbuf(1:npar_tot) < procz_bounds(ipz+1)
      endif inz
!
!  Count local particles.
!
      nv = count(lpar_loc)
      cknv: if (nv > mv) then
        print *, "input_part_snap: iproc, nv, mv = ", iproc, nv, mv
        call fatal_error_local("input_part_snap", "too many local particles")
      endif cknv
      call fatal_error_local_collect()
!
      allocate(indices(nv), stat=mpi_err)
      call check_success_local("input_part", "allocate indices")
      call fatal_error_local_collect()
      k = 0
      par: do i = 1, npar_tot
        loc: if (lpar_loc(i)) then
          k = k + 1
          indices(k) = i - 1
        endif loc
      enddo par
      ! The following was inefficient for large npar_tot.
      !indices = pack((/ (i, i = 0, npar_tot - 1) /), lpar_loc)
!
!  Decompose the integer data domain and read.
!
      call MPI_TYPE_INDEXED(nv, spread(1,1,nv), indices, MPI_INTEGER, ftype, mpi_err)
      call check_success_local("input_part", "create MPI data type")
      call fatal_error_local_collect()
!
      call MPI_TYPE_COMMIT(ftype, mpi_err)
      call check_success_local("input_part", "commit MPI data type")
      call fatal_error_local_collect()
!
      call MPI_FILE_SET_VIEW(handle, size_of_int, MPI_INTEGER, ftype, "native", io_info, mpi_err)
      call check_success("input_part", "set view of", fpath)
!
      call MPI_FILE_READ_ALL(handle, ipar, nv, MPI_INTEGER, status, mpi_err)
      call check_success("input_part", "read particle IDs", fpath)
!
      call MPI_TYPE_FREE(ftype, mpi_err)
      call check_success_local("input_part", "free MPI data type")
      call fatal_error_local_collect()
!
!  Decompose the real data domain and read.
!
      call MPI_TYPE_INDEXED(nv, spread(1,1,nv), indices, mpi_precision, ftype, mpi_err)
      call check_success_local("input_part", "create MPI data type")
      call fatal_error_local_collect()
!
      call MPI_TYPE_COMMIT(ftype, mpi_err)
      call check_success_local("input_part", "commit MPI data type")
      call fatal_error_local_collect()
!
      dsize = int(npar_tot, KIND=MPI_OFFSET_KIND) * size_of_real
      real1: do i = 1, mparray
        call MPI_FILE_SET_VIEW(handle, offset, mpi_precision, ftype, "native", io_info, mpi_err)
        call check_success("input_part", "set view of", fpath)
!
        call MPI_FILE_READ_ALL(handle, rbuf, nv, mpi_precision, status, mpi_err)
        call check_success("input_part", "read particle data", fpath)
!
        ap(1:nv,i) = rbuf(1:nv)
        offset = offset + dsize
      enddo real1
!
      call MPI_TYPE_FREE(ftype, mpi_err)
      call check_success_local("input_part", "free MPI data type")
      call fatal_error_local_collect()
!
      deallocate(lpar_loc, indices)
!
!  Close snapshot file.
!
      call MPI_FILE_CLOSE(handle, mpi_err)
      call check_success("input_part", "close", fpath)
!
      deallocate(rbuf)
!
    endsubroutine input_part_snap
!***********************************************************************
    subroutine input_pointmass(file, labels, fq, mv, nc)
!
!  Read pointmass snapshot file.
!
!  26-Oct-2018/PABourdin: coded
!
      use General, only: keep_compiler_quiet
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mv, nc
      character (len=*), dimension (nc), intent(in) :: labels
      real, dimension (mv,nc), intent(out) :: fq
!
      call fatal_error ('input_pointmass', 'not implemented for "io_mpi2"', .true.)
      call keep_compiler_quiet(file)
      call keep_compiler_quiet(labels(1))
      call keep_compiler_quiet(fq)
      call keep_compiler_quiet(mv, nc)
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
      logical, dimension(size(value)) :: buffer
!
      write_persist_logical_1D = .true.
      if (write_persist_id (label, id)) return
!
      nv = size (value)
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_logical_1D', &
            'Could not allocate memory for global', .true.)
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
        deallocate (global)
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
      integer, dimension(size(value)) :: buffer
!
      write_persist_int_1D = .true.
      if (write_persist_id (label, id)) return
!
      nv = size (value)
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_int_1D', &
            'Could not allocate memory for global', .true.)
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
        deallocate (global)
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
      real, dimension(size(value)) :: buffer
!
      write_persist_real_1D = .true.
      if (write_persist_id (label, id)) return
!
      nv = size (value)
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_real_1D', &
            'Could not allocate memory for global', .true.)
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
        deallocate (global)
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
!
    endfunction init_read_persist
!***********************************************************************
    logical function persist_exists(label)
!
!  Dummy routine
!
!  12-Oct-2019/PABourdin: coded
!
      character (len=*), intent(in), optional :: label
!
      if (present(label)) call warning("persist_exists", "The argument label has no effects. ")
!
      persist_exists = .false.
!
    endfunction persist_exists
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
    logical function write_persist_torus_rect(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  16-May-2020/MR: coded
!
      use Geometrical_types

      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      type(torus_rect), intent(in), optional :: value
!
      if (present(value)) call warning("write_persist_torus_rect", "The argument value has no effects. ")
!
      write_persist_torus_rect = .true.
      if (write_persist_id (label, id)) return
!
      !write (lun_output) value
      write_persist_torus_rect = .false.
!
    endfunction write_persist_torus_rect
!***********************************************************************
    logical function read_persist_torus_rect(label,value)
!
!  Read persistent data from snapshot file.
!
!  16-May-2020/MR: coded
!
      use Geometrical_types

      character (len=*), intent(in), optional :: label
      type(torus_rect), intent(out), optional :: value
!
      if (present(label)) call warning("read_persist_torus_rect", "The argument label has no effects. ")
      if (present(value)) call warning("read_persist_torus_rect", "The argument value has no effects. ")
!
      !read (lun_input) value
      read_persist_torus_rect = .false.
!
    endfunction read_persist_torus_rect
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
      if (present(label)) call warning("output_globals", "The argument label has no effects. ")
!
      call output_snap (a, nv2=nv, file=file, mode=0)
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
!  10-Feb-2012/PABourdin: adapted for collective IO
!
      use Mpicomm, only: collect_grid
      use General, only: loptest
!
      character (len=*) :: file
      integer, optional :: mxout,myout,mzout
      logical, optional :: lwrite
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (present(mxout)) call warning("wgrid", "The argument mxout has no effects. ")
      if (present(myout)) call warning("wgrid", "The argument myout has no effects. ")
      if (present(mzout)) call warning("wgrid", "The argument mzout has no effects. ")
!
      if (lyang) return      ! grid collection only needed on Yin grid, as grids are identical

      if (loptest(lwrite,.not.luse_oldgrid)) then
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
    subroutine wproc_bounds_mpi(file)
!
! Export processor boundaries to file.
!
! 12-nov-20/ccyang: coded
!
      character(len=*), intent(in) :: file
!
      integer :: handle
!
! Open file.
!
      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(file), ior(MPI_MODE_CREATE, MPI_MODE_WRONLY), io_info, handle, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("wproc_bounds_mpi", "could not open file "//trim(file)//" for writing")
!
! Write proc[xyz]_bounds.
!
      wproc: if (lroot) then
!
        call MPI_FILE_WRITE(handle, procx_bounds, nprocx + 1, mpi_precision, status, mpi_err)
        if (mpi_err /= MPI_SUCCESS) call fatal_error_local("wproc_bounds_mpi", "could not write procx_bounds")
!
        call MPI_FILE_WRITE(handle, procy_bounds, nprocy + 1, mpi_precision, status, mpi_err)
        if (mpi_err /= MPI_SUCCESS) call fatal_error_local("wproc_bounds_mpi", "could not write procy_bounds")
!
        call MPI_FILE_WRITE(handle, procz_bounds, nprocz + 1, mpi_precision, status, mpi_err)
        if (mpi_err /= MPI_SUCCESS) call fatal_error_local("wproc_bounds_mpi", "could not write procz_bounds")
!
      endif wproc
      call fatal_error_local_collect()
!
! Close file.
!
      call MPI_FILE_CLOSE(handle, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("wproc_bounds_mpi", "could not close file " // trim(file))
!
    endsubroutine wproc_bounds_mpi
!***********************************************************************
    subroutine rproc_bounds_mpi(file)
!
! Import processor boundaries from file.
!
! 12-nov-20/ccyang: coded
!
      character(len=*), intent(in) :: file
!
      integer :: handle
!
! Open file.
!
      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(file), MPI_MODE_RDONLY, io_info, handle, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("rproc_bounds_mpi", "could not open file "//trim(file)//" for reading")
!
      call MPI_FILE_SET_VIEW(handle, 0_MPI_OFFSET_KIND, mpi_precision, mpi_precision, "native", io_info, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("rproc_bounds_mpi", "could not set view")
!
! Read proc[xyz]_bounds.
!
      call MPI_FILE_READ_ALL(handle, procx_bounds, nprocx + 1, mpi_precision, status, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("rproc_bounds_mpi", "could not read procx_bounds")
!
      call MPI_FILE_READ_ALL(handle, procy_bounds, nprocy + 1, mpi_precision, status, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("rproc_bounds_mpi", "could not read procy_bounds")
!
      call MPI_FILE_READ_ALL(handle, procz_bounds, nprocz + 1, mpi_precision, status, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("rproc_bounds_mpi", "could not read procz_bounds")
!
! Close file.
!
      call MPI_FILE_CLOSE(handle, mpi_err)
      if (mpi_err /= MPI_SUCCESS) call fatal_error("rproc_bounds_mpi", "could not close file " // trim(file))
!
    endsubroutine rproc_bounds_mpi
!***********************************************************************
    include 'io_common.inc'
!***********************************************************************
endmodule Io
