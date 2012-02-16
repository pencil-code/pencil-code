! $Id$
!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   io_collect.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  I/O via MPI root rank by collecting data from all processors.
!  (storing data into one file, e.g. data/allprocs/var.dat)
!
!  The file written by output() (and used e.g. for 'var.dat')
!  consists of the followinig records (not using record markers):
!    1. data(mxgrid,mygrid,mzgrid,nvar)
!    2. t(1), x(mxgrid), y(mygrid), z(mzgrid), dx(1), dy(1), dz(1)
!    3. deltay(1) [optional: if lshear==.true.]
!  Where mxgrid is nxgrid+2*nghost, the same applies for mygrid and mzgrid,
!  and nvar denotes the number of variables to be saved.
!  In the case of MHD with entropy, nvar is 8 for a 'var.dat' file.
!
!  13-Jan-2012/Bourdin.KIS: adapted from io_dist.f90
!
module Io
!
  use Cdata
  use Cparam, only: intlen, fnlen, max_int
  use Messages, only: fatal_error, outlog, svn_id
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
  ! define unique logical unit number for input and output calls
  integer :: lun_input=88
  integer :: lun_output=91
!
  ! Indicates if IO is done distributed (each proc writes into a procdir)
  ! or collectively (eg. by specialized IO-nodes or by MPI-IO).
  logical :: lcollective_IO=.true.
!
  logical :: persist_initialized=.false.
  integer :: persist_last_id=-max_int
!
contains
!***********************************************************************
    subroutine register_io()
!
!  dummy routine, generates separate directory for each processor.
!  VAR#-files are written to the directory directory_snap which will
!  be the same as directory, unless specified otherwise.
!
!  04-jul-2011/Boudin.KIS: coded
!
      use Mpicomm, only: lroot
!
!  identify version number
!
      if (lroot) call svn_id ("$Id$")
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
      use Cdata, only: datadir, directory, datadir_snap, directory_snap
      use General, only: safe_character_assign
      use Mpicomm, only: lroot
!
      if ((datadir_snap == '') .or. (index (datadir_snap, 'allprocs') > 0)) then
        datadir_snap = datadir
      endif
!
      call safe_character_assign (directory, trim (datadir)//'/allprocs')
      call safe_character_assign (directory_snap, trim (datadir_snap)//'/allprocs')
!
    endsubroutine directory_names
!***********************************************************************
    subroutine collect_grid(x, y, z, gx, gy, gz)
!
!  This routine collects the global grid on the root processor.
!
!  05-jul-2011/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_real, mpirecv_real
!
      real, dimension(mx), intent(in) :: x
      real, dimension(my), intent(in) :: y
      real, dimension(mz), intent(in) :: z
      real, dimension(nxgrid+2*nghost), intent(out) :: gx
      real, dimension(nygrid+2*nghost), intent(out) :: gy
      real, dimension(nzgrid+2*nghost), intent(out) :: gz
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
    subroutine distribute_grid(gx, gy, gz, x, y, z)
!
!  This routine distributes the global grid to all processors.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_real, mpirecv_real
!
      real, dimension(nxgrid+2*nghost), intent(in) :: gx
      real, dimension(nygrid+2*nghost), intent(in) :: gy
      real, dimension(nzgrid+2*nghost), intent(in) :: gz
      real, dimension(mx), intent(out) :: x
      real, dimension(my), intent(out) :: y
      real, dimension(mz), intent(out) :: z
!
      integer :: px, py, pz, partner
      integer, parameter :: tag_gx=677, tag_gy=678, tag_gz=679
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
    subroutine output_snap(file, a, nv, mode)
!
!  write snapshot file, always write mesh and time, could add other things.
!
!  10-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, globalize_xy, mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: a
      integer, optional, intent(in) :: mode
!
      real, dimension (:,:,:,:), allocatable :: ga
      real, dimension (:), allocatable :: gx, gy, gz
      integer, parameter :: ngx=nxgrid+2*nghost, ngy=nygrid+2*nghost, ngz=nzgrid+2*nghost
      integer, parameter :: tag_ga=676
      integer :: pz, pa, io_len, alloc_err, z_start, z_end
      logical :: lwrite_add
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (lfirst_proc_xy) then
        allocate (ga(ngx,ngy,mz,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('output_snap', 'Could not allocate memory for ga', .true.)
        if (lroot) then
          allocate (gx(ngy), gy(ngy), gz(ngz), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('output_snap', 'Could not allocate memory for gx,gy,gz', .true.)
        endif
      endif
!
      lwrite_add = .true.
      if (present (mode)) lwrite_add = (mode == 1)
      if (lwrite_add) call collect_grid (x, y, z, gx, gy, gz)
!
      if (lroot) then
        inquire (IOLENGTH=io_len) t_sp
        open (lun_output, FILE=file, status='replace', access='direct', recl=ngx*ngy*io_len)
!
        ! iterate through xy-leading processors in the z-direction
        do pz = 0, nprocz-1
          if (pz == 0) then
            ! receive data from the local xy-plane
            call globalize_xy (a, ga)
          else
            ! receive collected data from leading processors in the xy-planes
            call mpirecv_real (ga, (/ ngx, ngy, mz, nv /), pz*nprocxy, tag_ga)
          endif
          z_start = n1
          z_end = n2
          if (pz == 0) z_start = 1
          if (pz == nprocz-1) z_end = n2 + nghost
          ! iterate through variables
          do pa = 1, mvar_io
            ! iterate through xy-planes and write each plane separately
            do iz = z_start, z_end
              write (lun_output, rec=iz+pz*nz+(pa-1)*ngz) ga(:,:,iz,pa)
            enddo
          enddo
        enddo
!
        ! write additional data:
        if (lwrite_add) then
          close (lun_output)
          open (lun_output, FILE=file, FORM='unformatted', position='append')
          t_sp = t
          write (lun_output) t_sp, gx, gy, gz, dx, dy, dz
          if (lshear) write (lun_output) deltay
        endif
!
      else
        if (lfirst_proc_xy) then
          ! receive data from the local xy-plane
          call globalize_xy (a, ga)
          ! send collected data to root processor
          call mpisend_real (ga, (/ ngx, ngy, mz, nv /), 0, tag_ga)
        else
          ! send data to leading processor in the local xy-plane
          call globalize_xy (a)
        endif
      endif
!
    endsubroutine output_snap
!***********************************************************************
    subroutine output_snap_finalize()
!
!  Close snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot
!
      integer :: io_err
      logical :: lerror
!
      if (persist_initialized) then
        if (lroot) then
          write (lun_output, iostat=io_err) id_block_PERSISTENT
          lerror = outlog (io_err, 'write id_block_PERSISTENT')
        endif
        persist_initialized = .false.
        persist_last_id = -max_int
      endif
!
      if (lroot) close (lun_output)
!
    endsubroutine output_snap_finalize
!***********************************************************************
    subroutine input_snap(file, a, nv, mode)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!  10-Feb-2012/Bourdin.KIS: coded
!
      use Cdata
      use Mpicomm, only: lroot, localize_xy, mpisend_real, mpirecv_real, mpibcast_real
      use Syscalls, only: sizeof_real
!
      character (len=*) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(out) :: a
      integer, optional, intent(in) :: mode
!
      real, dimension (:,:,:,:), allocatable :: ga
      real, dimension (:), allocatable :: gx, gy, gz
      integer, parameter :: ngx=nxgrid+2*nghost, ngy=nygrid+2*nghost, ngz=nzgrid+2*nghost
      integer, parameter :: tag_ga=675
      integer :: pz, pa, z_start, z_end, io_len, alloc_err, io_err
      logical :: lread_add
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (lfirst_proc_xy) then
        allocate (ga(ngx,ngy,mz,nv), gx(ngx), gy(ngy), gz(ngz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('input_snap', 'Could not allocate memory for ga,gx,gy,gz', .true.)
      endif
!
      lread_add = .true.
      if (present (mode)) lread_add = (mode == 1)
!
      if (lroot) then
        if (ip <= 8) print *, 'input_snap: open ', file
        inquire (IOLENGTH=io_len) t_sp
        open (lun_input, FILE=file, access='direct', recl=ngx*ngy*io_len)
!
        if (ip <= 8) print *, 'input_snap: read dim=', ngx, ngy, ngz, nv
        ! iterate through xy-leading processors in the z-direction
        do pz = 0, nprocz-1
          z_start = n1
          z_end = n2
          if (pz == 0) z_start = 1
          if (pz == nprocz-1) z_end = n2 + nghost
          ! iterate through variables
          do pa = 1, mvar_io
            ! iterate through xy-planes and read each plane separately
            do iz = z_start, z_end
              read (lun_input, rec=iz+pz*nz+(pa-1)*ngz) ga(:,:,iz,pa)
            enddo
          enddo
          if (pz == 0) then
            ! distribute data in the local xy-plane
            call localize_xy (ga, a)
          else
            ! send collective data to the leading processors in the xy-planes
            call mpisend_real (ga, (/ ngx, ngy, mz, nv /), pz*nprocxy, tag_ga)
          endif
        enddo
!
        if (lread_add) then
          ! read additional data
          close (lun_input)
          open (lun_input, FILE=file, FORM='unformatted')
          call fseek (lun_input, ngx*ngy*ngz*nv*sizeof_real(), 0)
          read (lun_input, IOSTAT=io_err) t_sp, gx, gy, gz, dx, dy, dz
          if (io_err > 0) call fatal_error ('input_snap', 'Could not read additional data', .true.)
          if (lshear) read (lun_input) deltay
        endif
!
      else
        if (lfirst_proc_xy) then
          call mpirecv_real (ga, (/ ngx, ngy, mz, nv /), 0, tag_ga)
        endif
        call localize_xy (ga, a)
      endif
!
      if (lread_add) then
        call mpibcast_real (t_sp)
        t = t_sp
        if (lshear) call mpibcast_real (deltay)
      endif
!
    endsubroutine input_snap
!***********************************************************************
    subroutine input_snap_finalize()
!
!  Close snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      if (lroot) close (lun_input)
!
    endsubroutine input_snap_finalize
!***********************************************************************
    logical function init_write_persist()
!
!  Initialize writing of persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical
!
      integer :: io_err
!
      if (lroot) then
        write (lun_output, iostat=io_err) id_block_PERSISTENT
        init_write_persist = outlog (io_err, 'write id_block_PERSISTENT')
      endif
      call mpibcast_logical (init_write_persist)
      persist_initialized = .not. init_write_persist
      persist_last_id = -max_int
!
    endfunction init_write_persist
!***********************************************************************
    logical function write_persist_id(label, id)
!
!  Write persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
!
      integer :: io_err
!
      write_persist_id = .true.
      if (.not. persist_initialized) return
!
      if (persist_last_id /= id) then
        if (lroot) then
          write (lun_output, iostat=io_err) id
          write_persist_id = outlog (io_err, 'write persistent ID '//label)
        endif
        call mpibcast_logical (write_persist_id)
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
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, intent(in) :: value
!
      integer :: px, py, pz, partner, io_err, alloc_err
      integer, parameter :: tag_log_0D = 700
      logical, dimension (:,:,:), allocatable :: global
      logical :: buffer
!
      write_persist_logical_0D = .true.
      if (.not. persist_initialized) return
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
              call mpirecv_logical (buffer, 1, partner, tag_log_0D)
              global(px+1,py+1,pz+1) = buffer
            enddo
          enddo
        enddo
        write (lun_output, iostat=io_err) global
        write_persist_logical_0D = outlog (io_err, 'write persistent '//label)
!
        deallocate (global)
      else
        call mpisend_logical (value, 1, 0, tag_log_0D)
      endif
!
      call mpibcast_logical (write_persist_logical_0D)
!
    endfunction write_persist_logical_0D
!***********************************************************************
    logical function write_persist_logical_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, dimension(:), intent(in) :: value
!
      integer :: px, py, pz, partner, nv, io_err, alloc_err
      integer, parameter :: tag_log_1D = 701
      logical, dimension (:,:,:,:), allocatable :: global
      logical, dimension (:), allocatable :: buffer
!
      write_persist_logical_1D = .true.
      if (.not. persist_initialized) return
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
        write (lun_output, iostat=io_err) global
        write_persist_logical_1D = outlog (io_err, 'write persistent '//label)
!
        deallocate (global, buffer)
      else
        call mpisend_logical (value, nv, 0, tag_log_1D)
      endif
!
      call mpibcast_logical (write_persist_logical_1D)
!
    endfunction write_persist_logical_1D
!***********************************************************************
    logical function write_persist_int_0D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, intent(in) :: value
!
      integer :: px, py, pz, partner, io_err, alloc_err
      integer, parameter :: tag_int_0D = 702
      integer, dimension (:,:,:), allocatable :: global
      integer :: buffer
!
      write_persist_int_0D = .true.
      if (.not. persist_initialized) return
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
              call mpirecv_int (buffer, 1, partner, tag_int_0D)
              global(px+1,py+1,pz+1) = buffer
            enddo
          enddo
        enddo
        write (lun_output, iostat=io_err) global
        write_persist_int_0D = outlog (io_err, 'write persistent '//label)
!
        deallocate (global)
      else
        call mpisend_int (value, 1, 0, tag_int_0D)
      endif
!
      call mpibcast_logical (write_persist_int_0D)
!
    endfunction write_persist_int_0D
!***********************************************************************
    logical function write_persist_int_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      integer, dimension (:), intent(in) :: value
!
      integer :: px, py, pz, partner, nv, io_err, alloc_err
      integer, parameter :: tag_int_1D = 703
      integer, dimension (:,:,:,:), allocatable :: global
      integer, dimension (:), allocatable :: buffer
!
      write_persist_int_1D = .true.
      if (.not. persist_initialized) return
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
        write (lun_output, iostat=io_err) global
        write_persist_int_1D = outlog (io_err, 'write persistent '//label)
!
        deallocate (global, buffer)
      else
        call mpisend_int (value, nv, 0, tag_int_1D)
      endif
!
      call mpibcast_logical (write_persist_int_1D)
!
    endfunction write_persist_int_1D
!***********************************************************************
    logical function write_persist_real_0D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, intent(in) :: value
!
      integer :: px, py, pz, partner, io_err, alloc_err
      integer, parameter :: tag_real_0D = 704
      real, dimension (:,:,:), allocatable :: global
      real :: buffer
!
      write_persist_real_0D = .true.
      if (.not. persist_initialized) return
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
              call mpirecv_real (buffer, 1, partner, tag_real_0D)
              global(px+1,py+1,pz+1) = buffer
            enddo
          enddo
        enddo
        write (lun_output, iostat=io_err) global
        write_persist_real_0D = outlog (io_err, 'write persistent '//label)
!
        deallocate (global)
      else
        call mpisend_real (value, 1, 0, tag_real_0D)
      endif
!
      call mpibcast_logical (write_persist_real_0D)
!
    endfunction write_persist_real_0D
!***********************************************************************
    logical function write_persist_real_1D(label, id, value)
!
!  Write persistent data to snapshot file.
!
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      real, dimension (:), intent(in) :: value
!
      integer :: px, py, pz, partner, nv, io_err, alloc_err
      integer, parameter :: tag_real_1D = 705
      real, dimension (:,:,:,:), allocatable :: global
      real, dimension (:), allocatable :: buffer
!
      write_persist_real_1D = .true.
      if (.not. persist_initialized) return
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
        write (lun_output, iostat=io_err) global
        write_persist_real_1D = outlog (io_err, 'write persistent '//label)
!
        deallocate (global, buffer)
      else
        call mpisend_real (value, nv, 0, tag_real_1D)
      endif
!
      call mpibcast_logical (write_persist_real_1D)
!
    endfunction write_persist_real_1D
!***********************************************************************
    logical function read_persist_logical_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      logical, intent(out) :: value
!
      integer :: px, py, pz, partner, io_err, alloc_err
      integer, parameter :: tag_log_0D = 706
      logical, dimension (:,:,:), allocatable :: global
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_logical_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        read (lun_input, iostat=io_err) global
        read_persist_logical_0D = outlog (io_err, 'read persistent '//label)
        value = global(ipx+1,ipy+1,ipz+1)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpisend_logical (global(px+1,py+1,pz+1), 1, partner, tag_log_0D)
            enddo
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_logical (value, 1, 0, tag_log_0D)
      endif
!
      call mpibcast_logical (read_persist_logical_0D)
!
    endfunction read_persist_logical_0D
!***********************************************************************
    logical function read_persist_logical_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_logical, mpirecv_logical
!
      character (len=*), intent(in) :: label
      logical, dimension(:), intent(out) :: value
!
      integer :: px, py, pz, partner, nv, io_err, alloc_err
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
        read (lun_input, iostat=io_err) global
        read_persist_logical_1D = outlog (io_err, 'read persistent '//label)
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
      call mpibcast_logical (read_persist_logical_1D)
!
    endfunction read_persist_logical_1D
!***********************************************************************
    logical function read_persist_int_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, intent(out) :: value
!
      integer :: px, py, pz, partner, io_err, alloc_err
      integer, parameter :: tag_int_0D = 708
      integer, dimension (:,:,:), allocatable :: global
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_int_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        read (lun_input, iostat=io_err) global
        read_persist_int_0D = outlog (io_err, 'read persistent '//label)
        value = global(ipx+1,ipy+1,ipz+1)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpisend_int (global(px+1,py+1,pz+1), 1, partner, tag_int_0D)
            enddo
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_int (value, 1, 0, tag_int_0D)
      endif
!
      call mpibcast_logical (read_persist_int_0D)
!
    endfunction read_persist_int_0D
!***********************************************************************
    logical function read_persist_int_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: label
      integer, dimension(:), intent(out) :: value
!
      integer :: px, py, pz, partner, nv, io_err, alloc_err
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
        read (lun_input, iostat=io_err) global
        read_persist_int_1D = outlog (io_err, 'read persistent '//label)
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
      call mpibcast_logical (read_persist_int_1D)
!
    endfunction read_persist_int_1D
!***********************************************************************
    logical function read_persist_real_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      real, intent(out) :: value
!
      integer :: px, py, pz, partner, io_err, alloc_err
      integer, parameter :: tag_real_0D = 710
      real, dimension (:,:,:), allocatable :: global
!
      if (lroot) then
        allocate (global(nprocx,nprocy,nprocz), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_persist_real_0D', &
            'Could not allocate memory for global buffer', .true.)
!
        read (lun_input, iostat=io_err) global
        read_persist_real_0D = outlog (io_err, 'read persistent '//label)
        value = global(ipx+1,ipy+1,ipz+1)
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            do pz = 0, nprocz-1
              partner = px + py*nprocx + pz*nprocxy
              if (iproc == partner) cycle
              call mpisend_real (global(px+1,py+1,pz+1), 1, partner, tag_real_0D)
            enddo
          enddo
        enddo
!
        deallocate (global)
      else
        call mpirecv_real (value, 1, 0, tag_real_0D)
      endif
!
      call mpibcast_logical (read_persist_real_0D)
!
    endfunction read_persist_real_0D
!***********************************************************************
    logical function read_persist_real_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_logical, mpisend_real, mpirecv_real
!
      character (len=*), intent(in) :: label
      real, dimension(:), intent(out) :: value
!
      integer :: px, py, pz, partner, nv, io_err, alloc_err
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
        read (lun_input, iostat=io_err) global
        read_persist_real_1D = outlog (io_err, 'read persistent '//label)
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
      call mpibcast_logical (read_persist_real_1D)
!
    endfunction read_persist_real_1D
!***********************************************************************
    subroutine output_globals(file,a,nv)
!
!  Write snapshot file of globals, ignore time and mesh.
!
!  10-Feb-2012/Bourdin.KIS: coded
!
      character (len=*) :: file
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
!
      call output_snap (file, a, nv, 0)
      call output_snap_finalize ()
!
    endsubroutine output_globals
!***********************************************************************
    subroutine input_globals(file,a,nv)
!
!  Read globals snapshot file, ignore time and mesh.
!
!  10-Feb-2012/Bourdin.KIS: coded
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
      use Cdata, only: lroot, lcopysnapshots_exp, datadir
      use Cparam, only: fnlen
      use General, only: parse_filename
      use Mpicomm, only: mpibarrier
!
      character (len=*) :: filename, flist
      character (len=fnlen) :: dir, fpart
!
      call parse_filename (filename, dir, fpart)
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
!  Write processor-local part of grid coordinates.
!
!  10-Feb-2012/Bourdin.KIS: adapted for collective IO
!
      use Cdata
      use Mpicomm, only: lroot
!
      character (len=*) :: file
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err, io_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
      if (lroot) then
        allocate (gx(nxgrid+2*nghost), gy(nygrid+2*nghost), gz(nzgrid+2*nghost), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('wgrid', 'Could not allocate memory for gx,gy,gz', .true.)
      endif
!
      if (lroot) then
        open (lun_output, FILE=file, FORM='unformatted', IOSTAT=io_err)
        if (io_err /= 0) call fatal_error ('wgrid', &
            "Cannot open " // trim (file) // " (or similar) for writing" // &
            " -- is data/ visible from all nodes?", .true.)
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
      else
        call collect_grid (x, y, z, gx, gy, gz)
        call collect_grid (dx_1, dy_1, dz_1, gx, gy, gz)
        call collect_grid (dx_tilde, dy_tilde, dz_tilde, gx, gy, gz)
      endif
!
    endsubroutine wgrid
!***********************************************************************
    subroutine rgrid(file)
!
!  Read processor-local part of grid coordinates.
!
!  21-jan-02/wolf: coded
!  15-jun-03/axel: Lx,Ly,Lz are now read in from file (Tony noticed the mistake)
!  10-Feb-2012/Bourdin.KIS: adapted for collective IO
!
      use Cdata
      use Mpicomm, only: lroot, mpibcast_int, mpibcast_real
!
      character (len=*) :: file
!
      real, dimension (:), allocatable :: gx, gy, gz
      integer :: alloc_err, io_err
      real :: t_sp   ! t in single precision for backwards compatibility
!
      allocate (gx(nxgrid+2*nghost), gy(nygrid+2*nghost), gz(nzgrid+2*nghost), stat=alloc_err)
      if (alloc_err > 0) call fatal_error ('rgrid', 'Could not allocate memory for gx,gy,gz', .true.)
!
      if (lroot) then
        open (lun_input, FILE=file, FORM='unformatted', IOSTAT=io_err)
        if (io_err /= 0) call fatal_error ('rgrid', &
            "Cannot open " // trim(file) // " (or similar) for reading" // &
            " -- is data/ visible from all nodes?", .true.)
      endif
      if (lroot) read (lun_input) t_sp, gx, gy, gz, dx, dy, dz
      call distribute_grid (gx, gy, gz, x, y, z)
      if (lroot) read (lun_input) dx, dy, dz
      call mpibcast_real (dx)
      call mpibcast_real (dy)
      call mpibcast_real (dz)
      if (lroot) read (lun_input) Lx, Ly, Lz
      call mpibcast_real (Lx)
      call mpibcast_real (Ly)
      call mpibcast_real (Lz)
      if (lroot) read (lun_input) gx, gy, gz
      call distribute_grid (gx, gy, gz, dx_1, dy_1, dz_1)
      if (lroot) read (lun_input) gx, gy, gz
      call distribute_grid (gx, gy, gz, dx_tilde, dy_tilde, dz_tilde)
      if (lroot) close (lun_input)
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
!   20-aug-09/bourdin: adapted
!
      use Cdata, only: procy_bounds,procz_bounds
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
      integer :: ierr
!
      if (lroot) then
        open(lun_output,FILE=file,FORM='unformatted',IOSTAT=ierr)
        if (ierr /= 0) call stop_it( &
            "Cannot open " // trim(file) // " (or similar) for writing" // &
            " -- is data/ visible from all nodes?")
        write(lun_output) procy_bounds
        write(lun_output) procz_bounds
        close(lun_output)
      endif
!
    endsubroutine wproc_bounds
!***********************************************************************
    subroutine rproc_bounds(file)
!
!   Import processor boundaries from file.
!
!   20-aug-09/bourdin: adapted
!
      use Cdata, only: procy_bounds,procz_bounds
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
      integer :: ierr
!
      open(lun_input,FILE=file,FORM='unformatted',IOSTAT=ierr)
      if (ierr /= 0) call stop_it( &
          "Cannot open " // trim(file) // " (or similar) for reading" // &
          " -- is data/ visible from all nodes?")
      read(lun_input) procy_bounds
      read(lun_input) procz_bounds
      close(lun_input)
!
    endsubroutine rproc_bounds
!***********************************************************************
endmodule Io
