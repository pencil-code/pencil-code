! $Id$
!
!  I/O via MPI by collecting data from all processors in the xy-plane.
!  (storing data into files, e.g. data/proc(0,7,15,...)/var.dat)
!
!  The file written by output_snap() (and used e.g. for 'var.dat')
!  consists of the followinig records (not using F77 record markers):
!    1. data(mxgrid,mygrid,mzlocal,nvar)
!    2. marker(4 bytes), t(1)
!    3. marker(4 bytes), x(mxgrid), y(mygrid), z(mzgrid),
!       dx_1(mxgrid), dy_1(mygrid), dz_1(mzgrid),
!       dx_tilde(mxgrid), dy_tilde(mygrid), dz_tilde(mzgrid) [only in /proc0/]
!  Where nvar denotes the number of variables to be saved,
!  In the case of MHD with entropy, nvar is 8 for a 'var.dat' file.
!  Only outer ghost-layers are written, so mzlocal is between nz and mz,
!  depending on the corresponding ipz-layer.
!
!  To read these snapshots in IDL efficiently:
!  IDL> pc_read_var_raw, obj=data, tags=tags, grid=grid
!  or in the old-fashioned way, which is NOT recommended for big data:
!  IDL> pc_read_var, obj=vars
!
!  04-Sep-2015/PABourdin: adapted from io_collect_xy.f90
!
!  ================================================
!  ================================================
!  ===                                          ===
!  ===  PLEASE NEVER CHANGE THIS FILE YOURSELF  ===
!  ===                                          ===
!  ===  > In case of problems, please report <  ===
!  ===                                          ===
!  ================================================
!  ================================================
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
!
! Indicates if IO is done distributed (each proc writes into a procdir)
! or collectively (eg. by specialized IO-nodes or by MPI-IO).
!
  logical :: lcollective_IO=.true.
  character (len=labellen) :: IO_strategy="collect_xy_stream"
!
  character (len=8) :: record_marker, dead_beef = 'DEADBEEF'
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
!  identify version number
!
      if (lroot) call svn_id ("$Id$")
      if (ldistribute_persist) &
          call fatal_error ('io_collect_xy_stream', "Distibuted persistent variables are fatal with this IO method!")
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use General, only: directory_names_std
!
!  check whether directory_snap contains `/proc0' -- if so, revert to the
!  default name.
!  Rationale: if directory_snap was not explicitly set in start.in, it
!  will be written to param.nml as 'data/proc0', but this should in fact
!  be data/procN on processor N.
!
      if ((datadir_snap == '') .or. (index (datadir_snap,'proc0') > 0)) &
        datadir_snap = datadir
!
      call directory_names_std(.true.)
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
      integer :: px, py, pz, partner
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
    subroutine output_snap(a, nv1, nv2, file, mode)
!
!  write snapshot file, always write mesh and time, could add other things.
!
!  10-Feb-2012/Bourdin.KIS: coded
!  13-feb-2014/MR: made file optional (prep for downsampled output)
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: globalize_xy, collect_grid
      use General, only: ioptest
!
      integer,           optional,intent(IN) :: nv1,nv2
      real, dimension (:,:,:,:), intent(in) :: a
      character (len=*), optional, intent(in) :: file
      integer, optional, intent(in) :: mode
!
      real, dimension (:,:,:,:), allocatable :: ga
      real, dimension (:), allocatable :: gx, gy, gz
      integer, parameter :: tag_ga = 676
      integer :: alloc_err,na,ne
      logical :: lwrite_add
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (.not.present(file)) call fatal_error('output_snap', &
          'downsampled output not implemented for IO_collect_xy_stream')
!
      lwrite_add = .true.
      if (present (mode)) lwrite_add = (mode == 1)
!
      na=ioptest(nv1,1)
      ne=ioptest(nv2,mvar_io)
!
      if (lfirst_proc_xy) then
        allocate (ga(mxgrid,mygrid,mz,ne-na+1), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('output_snap', 'Could not allocate memory for ga', .true.)
!
        ! receive data from the xy-plane of the pz-layer
        call globalize_xy (a(:,:,:,na:ne), ga)

        open (lun_output, FILE=trim (directory_snap)//'/'//file, status='replace', access='stream', form='unformatted')
        write (lun_output) ga
        deallocate (ga)
!
      else
        ! send data to root processor
        call globalize_xy (a(:,:,:,na:ne))
      endif
!
      ! write additional data:
      if (lwrite_add) then
        if (lfirst_proc_xy) then
          if (lroot) then
            allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
            if (alloc_err > 0) call fatal_error ('output_snap', 'Could not allocate memory for gx,gy,gz', .true.)
          endif
!
          t_sp = real (t)
          write (lun_output) dead_beef, t_sp
        endif

        call collect_grid (x, y, z, gx, gy, gz)
        if (lroot) write (lun_output) dead_beef, gx, gy, gz, dx, dy, dz
        call collect_grid (dx_1, dy_1, dz_1, gx, gy, gz)
        if (lroot) write (lun_output) gx, gy, gz
        call collect_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
        if (lroot) write (lun_output) gx, gy, gz
      endif
!
    endsubroutine output_snap
!***********************************************************************
    subroutine output_snap_finalize
!
!  Close snapshot file.
!
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      if (lfirst_proc_xy) then
        if (persist_initialized) then
          if (lroot .and. (ip <= 9)) write (*,*) 'finish persistent block'
          write (lun_output) id_block_PERSISTENT
          persist_initialized = .false.
        endif
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
      call fatal_error ('output_part_snap', 'not implemented for "io_collect_xy_f2003"', .true.)
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
      call fatal_error ('output_stalker_init', 'not implemented for "io_collect_xy_f2003"', .true.)
!
    endsubroutine output_stalker_init
!***********************************************************************
    subroutine output_stalker(label, mv, nv, data, nvar, lfinalize)
!
!  Write stalker particle quantity to snapshot file.
!
!  03-May-2019/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: mv, nv
      real, dimension (mv), intent(in) :: data
      logical, intent(in), optional :: lfinalize
      integer, intent(in), optional :: nvar
!
      call fatal_error ('output_stalker', 'not implemented for "io_collect_xy_f2003"', .true.)
!
    endsubroutine output_stalker
!***********************************************************************
    subroutine output_part_finalize
!
!  Close particle snapshot file.
!
!  03-May-2019/PABourdin: coded
!
      call fatal_error ('output_part_finalize', 'not implemented for "io_collect_xy_f2003"', .true.)
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
      call fatal_error ('output_pointmass', 'not implemented for "io_collect_xy_f2003"', .true.)
!
    endsubroutine output_pointmass
!***********************************************************************
    subroutine input_snap(file, a, nv, mode)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use File_io, only: backskip_to_time
      use Mpicomm, only: localize_xy, mpibcast_real, stop_it_if_any, MPI_COMM_WORLD
      use Syscalls, only: sizeof_real
!
      character (len=*) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(out) :: a
      integer, optional, intent(in) :: mode
!
      real, dimension (:,:,:,:), allocatable :: ga
      real, dimension (:), allocatable :: gx, gy, gz
      integer, parameter :: tag_ga = 675
      integer :: alloc_err
      logical :: lread_add
      real :: t_sp, t_test   ! t in single precision for backwards compatibility
!
      lread_add = .true.
      if (present (mode)) lread_add = (mode == 1)
!
      if (lfirst_proc_xy) then
        allocate (ga(mxgrid,mygrid,mz,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('input_snap', 'Could not allocate memory for ga', .true.)
!
        if (ip <= 8) print *, 'input_snap: open ', file
!
        open (lun_input, FILE=trim (directory_snap)//'/'//file, access='stream', form='unformatted', status='old')
        read (lun_input) ga
!
        ! distribute data in the xy-plane of the pz-layer
        call localize_xy (a, ga)
        deallocate (ga)
!
      else
        ! receive data from root processor of pz-layer
        call localize_xy (a)
      endif
!
      ! read additional data
      if (lread_add) then
        if (lfirst_proc_xy) then
          read (lun_input) record_marker, t_sp
          t_test = t_sp
        endif
!
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('input_snap', 'Could not allocate memory for gx,gy,gz', .true.)
          read (lun_input) record_marker, gx, gy, gz, dx, dy, dz
          call distribute_grid (x, y, z, gx, gy, gz)
          read (lun_input) gx, gy, gz
          call distribute_grid (dx_1, dy_1, dz_1, gx, gy, gz)
          read (lun_input) gx, gy, gz
          call distribute_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
          deallocate (gx, gy, gz)
        else
          call distribute_grid (x, y, z)
          call distribute_grid (dx_1, dy_1, dz_1)
          call distribute_grid (dx_tilde, dy_tilde, dz_tilde)
        endif
!
        call mpibcast_real (t_sp,comm=MPI_COMM_WORLD)
        if (.not. lfirst_proc_xy) t_test = t_sp
        if (t_test /= t_sp) &
            write (*,*) 'ERROR: '//trim(directory_snap)//'/'//trim(file)//' IS INCONSISTENT: t=', t_sp
        call stop_it_if_any ((t_test /= t_sp), '')
        t = t_sp
      endif
!
    endsubroutine input_snap
!***********************************************************************
    subroutine input_snap_finalize
!
!  Close snapshot file.
!
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      if (persist_initialized) then
        persist_initialized = .false.
        persist_last_id = -max_int
      endif
!
      if (lfirst_proc_xy) close (lun_input)
!
    endsubroutine input_snap_finalize
!***********************************************************************
    subroutine input_part_snap(ipar, ap, mv, nv, npar_total, file, label)
!
!  Read particle snapshot file, mesh and time are read in 'input_snap'.
!
!  25-Oct-2018/PABourdin: apadpted and moved to IO module
!
      integer, intent(in) :: mv
      integer, dimension (mv), intent(out) :: ipar
      real, dimension (mv,mparray), intent(out) :: ap
      integer, intent(out) :: nv, npar_total
      character (len=*), intent(in) :: file
      character (len=*), optional, intent(in) :: label
!
      call fatal_error ('input_part_snap', 'not implemented for "io_collect_xy_f2003"', .true.)
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
      call fatal_error ('input_pointmass', 'not implemented for "io_collect_xy_f2003"', .true.)
!
    endsubroutine input_pointmass
!***********************************************************************
    logical function init_write_persist(file)
!
!  Initialize writing of persistent data to persistent file.
!
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
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
      if (lfirst_proc_xy) then
        if (lroot .and. (ip <= 9)) write (*,*) 'begin persistent block'
        if (filename /= "") then
          if (lfirst_proc_xy) close (lun_output)
          open (lun_output, FILE=trim(directory_snap)//'/'//filename, FORM='unformatted', status='replace')
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
!
      write_persist_id = .true.
      if (.not. persist_initialized) write_persist_id = init_write_persist()
      if (.not. persist_initialized) return
!
      if (lfirst_proc_xy .and. (persist_last_id /= id)) then
        if (lroot .and. (ip <= 9)) write (*,*) 'write persistent ID '//trim (label)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, intent(in) :: value
!
      integer :: px, py, partner, alloc_err
      integer, parameter :: tag_log_0D = 700
      logical, dimension (:,:), allocatable :: global
      logical :: buffer
!
      write_persist_logical_0D = .true.
      if (write_persist_id (label, id)) return
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_logical_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpirecv_logical (buffer, partner, tag_log_0D)
            global(px+1,py+1) = buffer
          enddo
        enddo
        if (lroot .and. (ip <= 9)) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global)
      else
        call mpisend_logical (value, ipz*nprocxy, tag_log_0D)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, dimension(:), intent(in) :: value
!
      integer :: px, py, partner, nv, alloc_err
      integer, parameter :: tag_log_1D = 701
      logical, dimension (:,:,:), allocatable :: global
      logical, dimension (:), allocatable :: buffer
!
      write_persist_logical_1D = .true.
      if (write_persist_id (label, id)) return
!
      nv = size (value)
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy,nv), buffer(nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_logical_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1,:) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpirecv_logical (buffer, nv, partner, tag_log_1D)
            global(px+1,py+1,:) = buffer
          enddo
        enddo
        if (lroot .and. (ip <= 9)) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global, buffer)
      else
        call mpisend_logical (value, nv, ipz*nprocxy, tag_log_1D)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, intent(in) :: value
!
      integer :: px, py, partner, alloc_err
      integer, parameter :: tag_int_0D = 702
      integer, dimension (:,:), allocatable :: global
      integer :: buffer
!
      write_persist_int_0D = .true.
      if (write_persist_id (label, id)) return
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_int_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpirecv_int (buffer, partner, tag_int_0D)
            global(px+1,py+1) = buffer
          enddo
        enddo
        if (lroot .and. (ip <= 9)) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global)
      else
        call mpisend_int (value, ipz*nprocxy, tag_int_0D)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, dimension (:), intent(in) :: value
!
      integer :: px, py, partner, nv, alloc_err
      integer, parameter :: tag_int_1D = 703
      integer, dimension (:,:,:), allocatable :: global
      integer, dimension (:), allocatable :: buffer
!
      write_persist_int_1D = .true.
      if (write_persist_id (label, id)) return
!
      nv = size (value)
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy,nv), buffer(nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_int_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1,:) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpirecv_int (buffer, nv, partner, tag_int_1D)
            global(px+1,py+1,:) = buffer
          enddo
        enddo
        if (lroot .and. (ip <= 9)) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global, buffer)
      else
        call mpisend_int (value, nv, ipz*nprocxy, tag_int_1D)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: collect_xy
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, intent(in) :: value
!
      integer :: alloc_err
      integer, parameter :: tag_real_0D = 704
      real, dimension (:,:), allocatable :: global
!
      write_persist_real_0D = .true.
      if (write_persist_id (label, id)) return
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_real_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        call collect_xy (value, global)
        if (lroot .and. (ip <= 9)) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global)
      else
        call collect_xy (value)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, dimension (:), intent(in) :: value
!
      integer :: px, py, partner, nv, alloc_err
      integer, parameter :: tag_real_1D = 705
      real, dimension (:,:,:), allocatable :: global
      real, dimension (:), allocatable :: buffer
!
      write_persist_real_1D = .true.
      if (write_persist_id (label, id)) return
!
      nv = size (value)
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy,nv), buffer(nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('write_persist_real_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        global(ipx+1,ipy+1,:) = value
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpirecv_real (buffer, nv, partner, tag_real_1D)
            global(px+1,py+1,:) = buffer
          enddo
        enddo
        if (lroot .and. (ip <= 9)) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global, buffer)
      else
        call mpisend_real (value, nv, ipz*nprocxy, tag_real_1D)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpibcast_logical, MPI_COMM_WORLD
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
      if (lfirst_proc_xy) then
        if (lroot .and. (ip <= 9)) write (*,*) 'begin persistent block'
        if (present (file)) then
          if (lfirst_proc_xy) close (lun_input)
          open (lun_input, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', status='old')
        endif
      endif
!
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
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
      if (lfirst_proc_xy) then
        if (lroot .and. (ip <= 9)) write (*,*) 'read persistent ID '//trim (label)
        if (lcatch_error) then
          if (lfirst_proc_xy) then
            read (lun_input, iostat=io_err) id
            if (io_err /= 0) id = -max_int
          endif
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      logical, intent(out) :: value
!
      integer :: px, py, partner, alloc_err
      integer, parameter :: tag_log_0D = 706
      logical, dimension (:,:), allocatable :: global
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_logical_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (lroot .and. (ip <= 9)) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpisend_logical (global(px+1,py+1), partner, tag_log_0D)
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_logical (value, ipz*nprocxy, tag_log_0D)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      logical, dimension(:), intent(out) :: value
!
      integer :: px, py, partner, nv, alloc_err
      integer, parameter :: tag_log_1D = 707
      logical, dimension (:,:,:), allocatable :: global
!
      nv = size (value)
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_logical_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (lroot .and. (ip <= 9)) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1,:)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpisend_logical (global(px+1,py+1,:), nv, partner, tag_log_1D)
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_logical (value, nv, ipz*nprocxy, tag_log_1D)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, intent(out) :: value
!
      integer :: px, py, partner, alloc_err
      integer, parameter :: tag_int_0D = 708
      integer, dimension (:,:), allocatable :: global
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_int_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (lroot .and. (ip <= 9)) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpisend_int (global(px+1,py+1), partner, tag_int_0D)
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_int (value, ipz*nprocxy, tag_int_0D)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, dimension(:), intent(out) :: value
!
      integer :: px, py, partner, nv, alloc_err
      integer, parameter :: tag_int_1D = 709
      integer, dimension (:,:,:), allocatable :: global
!
      nv = size (value)
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_int_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (lroot .and. (ip <= 9)) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1,:)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpisend_int (global(px+1,py+1,:), nv, partner, tag_int_1D)
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_int (value, nv, ipz*nprocxy, tag_int_1D)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      real, intent(out) :: value
!
      integer :: px, py, partner, alloc_err
      integer, parameter :: tag_real_0D = 710
      real, dimension (:,:), allocatable :: global
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_real_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (lroot .and. (ip <= 9)) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpisend_real (global(px+1,py+1), partner, tag_real_0D)
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_real (value, ipz*nprocxy, tag_real_0D)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      real, dimension(:), intent(out) :: value
!
      integer :: px, py, partner, nv, alloc_err
      integer, parameter :: tag_real_1D = 711
      real, dimension (:,:,:), allocatable :: global
!
      nv = size (value)
!
      if (lfirst_proc_xy) then
        allocate (global(nprocx,nprocy,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_real_1D', &
            'Could not allocate memory for global buffer', .true.)
!
        if (lroot .and. (ip <= 9)) write (*,*) 'read persistent '//trim (label)
        read (lun_input) global
        value = global(ipx+1,ipy+1,:)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpisend_real (global(px+1,py+1,:), nv, partner, tag_real_1D)
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_real (value, nv, ipz*nprocxy, tag_real_1D)
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
!  Write snapshot file of globals, ignore time and mesh.
!
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*), intent(in), optional :: label
!
      call output_snap (a, nv2=nv, file=file, mode=1)
      call output_snap_finalize
!
    endsubroutine output_globals
!***********************************************************************
    subroutine input_globals(file, a, nv)
!
!  Read globals snapshot file, ignore time and mesh.
!
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
!
      call input_snap (file, a, nv, 1)
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
      if (dir == '.') call safe_character_assign (dir, directory_collect)
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
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
      if (lyang) return      ! grid collection only needed on Yin grid, as grids are identical

      if (loptest(lwrite,.not.luse_oldgrid)) then
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
        if (lroot)  write (lun_output) gx, gy, gz

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
!  10-Feb-2012/Bourdin.KIS: adapted for collective IO
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
!
      use Mpicomm, only: mpibcast_int, mpibcast_real, MPI_COMM_WORLD
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
        open (lun_input, FILE=trim (directory_collect)//'/'//file, FORM='unformatted', status='old')
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
      if (lroot .and. (ip <= 4)) then
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
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
!  04-Sep-2015/PABourdin: adapted from 'io_collect_xy'
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
