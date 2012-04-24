! $Id$
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   io_collect_xy.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  I/O via MPI by collecting data from all processors in the xy-plane.
!  (storing data into files, e.g. data/proc(0,7,15,...)/var.dat)
!
!  The file written by output_snap() (and used e.g. for 'var.dat')
!  consists of the followinig records (not using record markers):
!    1. data(mxgrid,mygrid,mzlocal,nvar)
!    2. t(1)
!    3. x(mxgrid), y(mygrid), z(mzgrid), dx(1), dy(1), dz(1) [only in /proc0/]
!  Where nvar denotes the number of variables to be saved,
!  In the case of MHD with entropy, nvar is 8 for a 'var.dat' file.
!  Only outer ghost-layers are written, so mzlocal is between nz and mz,
!  depending on the corresponding ipz-layer.
!
!  To read these snapshots in IDL, the parameter allprocs=2 needs to be given:
!  IDL> pc_read_var, obj=vars, allprocs=2
!  or in a much more efficient way by reading into an array:
!  IDL> pc_read_var_raw, obj=data, tags=tags, grid=grid, allprocs=2
!
!  05-Mar-2012/Bourdin.KIS: adapted from io_collect.f90
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
  character (len=labellen) :: IO_strategy="collect_xy"
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
      if (.not. lseparate_persist) call fatal_error ('io_collect_xy', &
          "This module only works with the setting lseparate_persist=.true.")
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
      use Mpicomm, only: lroot
!
      character (len=intlen) :: chproc
!
!  check whether directory_snap contains `/proc0' -- if so, revert to the
!  default name.
!  Rationale: if directory_snap was not explicitly set in start.in, it
!  will be written to param.nml as 'data/proc0', but this should in fact
!  be data/procN on processor N.
!
      if ((datadir_snap == '') .or. (index(datadir_snap,'proc0')>0)) then
        datadir_snap = datadir
      endif
!
      chproc = itoa (iproc)
      call safe_character_assign (directory, trim (datadir)//'/proc'//chproc)
      call safe_character_assign (directory_dist, &
                                            trim (datadir_snap)//'/proc'//chproc)
      call safe_character_assign (directory_snap, trim (datadir_snap)//'/proc'//chproc)
      call safe_character_assign (directory_collect, trim (datadir_snap)//'/allprocs')
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
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_real, mpirecv_real
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
      integer, parameter :: tag_ga=676
      integer :: io_len, alloc_err
      integer(kind=8) :: rec_len
      logical :: lwrite_add
      real :: t_sp   ! t in single precision for backwards compatibility
!
      lwrite_add = .true.
      if (present (mode)) lwrite_add = (mode == 1)
!
      if (lfirst_proc_xy) then
        allocate (ga(mxgrid,mygrid,mz,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('output_snap', 'Could not allocate memory for ga', .true.)
!
        inquire (IOLENGTH=io_len) t_sp
        rec_len = int (mxgrid, kind=8) * int (mygrid, kind=8) * mz
        rec_len = rec_len * nv * io_len
        open (lun_output, FILE=trim (directory_snap)//'/'//file, status='replace', access='direct', recl=rec_len)
!
        ! receive data from the xy-plane of the pz-layer
        call globalize_xy (a, ga)
        write (lun_output, rec=1) ga
        deallocate (ga)
!
      else
        ! send data to root processor
        call globalize_xy (a)
      endif
!
      ! write additional data:
      if (lwrite_add) then
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('output_snap', 'Could not allocate memory for gx,gy,gz', .true.)
          call collect_grid (x, y, z, gx, gy, gz)
!
          close (lun_output)
          open (lun_output, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', position='append', status='old')
          t_sp = t
          write (lun_output) t_sp
          write (lun_output) gx, gy, gz, dx, dy, dz
          deallocate (gx, gy, gz)
        else
          call collect_grid (x, y, z)
!
          if (lfirst_proc_xy) then
            close (lun_output)
            open (lun_output, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', position='append', status='old')
            t_sp = t
            write (lun_output) t_sp
          endif
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
      if (persist_initialized) then
        if (lroot .and. (ip <= 9)) write (*,*) 'finish persistent block'
        write (lun_output) id_block_PERSISTENT
        persist_initialized = .false.
        if (.not. lfirst_proc_xy) close (lun_input)
      endif
!
      if (lfirst_proc_xy) close (lun_input)
!
    endsubroutine output_snap_finalize
!***********************************************************************
    subroutine output_form_int_0D(file,data,lappend)
!
!  Write formatted integer data to a file.
!  Set lappend to false to overwrite the file, default is to append the data.
!
!  19-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_int, mpirecv_int
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: data
      logical, intent(in), optional :: lappend
!
      integer, dimension(ncpus) :: data_global
      integer :: partner, buffer
      integer, parameter :: tag_data = 683
      logical :: lappend_opt
!
      lappend_opt = .false.
      if (present(lappend)) lappend_opt = lappend
!
      if (lroot) then
        ! collect data into global array
        data_global(iproc) = data
        do partner = 0, ncpus-1
          if (partner == iproc) cycle
          call mpirecv_int(buffer, 1, partner, tag_data)
          data_global(partner) = buffer
        enddo
        ! write global data to file
        if (lappend_opt) then
          open (lun_output, file=trim (directory_collect)//'/'//file, form='formatted', position='append', status='old')
        else
          open (lun_output, file=trim (directory_collect)//'/'//file, form='formatted', status='replace')
        endif
        write (lun_output,*) data_global
        close (lun_output)
      else
        ! send local data
        call mpisend_int(data, 1, 0, tag_data)
      endif
!
    endsubroutine output_form_int_0D
!***********************************************************************
    subroutine fseek_pos(unit, rec_len, num_rec, reference)
!
!  Seeks to a given position in an opened file relative to a reference point.
!  If reference=0, this is relative to the beginning of the file,
!  if reference=1, this is relative to the current position in the file,
!  and if reference=2, this is relative to the end of the file.
!  'rec_len' and 'num_rec' are referring to a record length and a given number
!  of records that one likes to seek, boths must be representable in integer.
!  If 'num_rec' is negative, seeking is done backwards.
!
!  20-Feb-2012/Bourdin.KIS: coded
!
      use General, only: itoa
!
      integer, intent(in) :: unit
      integer(kind=8) :: rec_len, num_rec
      integer, intent(in) :: reference
!
      integer :: i, num, len
!
      if (num_rec < 0) then
        num_rec = -num_rec
        rec_len = -rec_len
      endif
!
      ! all numbers must be representable as integer(kind=4)
      len = rec_len
      num = num_rec
      if (len /= rec_len) call fatal_error ('fseek_pos on unit '//trim (itoa (unit)), &
          "rec_len is not representable as integer(kind=4).", .true.)
      if (num /= num_rec) call fatal_error ('fseek_pos on unit '//trim (itoa (unit)), &
          "num_rec is not representable as integer(kind=4).", .true.)
!
! WORKAROUND:
! Even though the ifort manual states that ifort would be able to fseek
! with an 64-bit integer argument, this is NOT working!
! Therefore, we have to iterate the fseek with a 32-bit integer to be save.
! Note: gfortran would be able to seek with a 64-bit integer value, though.
! (20-Feb-2012, Bourdin.KIS)
!
      call fseek (unit, rec_len, reference)
      if (num >= 2) then
        do i = 2, num
          call fseek (unit, rec_len, 1)
        enddo
      endif
!
    endsubroutine fseek_pos
!***********************************************************************
    subroutine input_snap(file, a, nv, mode)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!  10-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, localize_xy, mpisend_real, mpirecv_real, mpibcast_real, stop_it_if_any
      use Syscalls, only: sizeof_real
!
      character (len=*) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(out) :: a
      integer, optional, intent(in) :: mode
!
      real, dimension (:,:,:,:), allocatable :: ga
      real, dimension (:), allocatable :: gx, gy, gz
      integer, parameter :: tag_ga=675
      integer :: io_len, alloc_err
      integer(kind=8) :: rec_len, num_rec
      logical :: lread_add
      real :: t_sp, t_test   ! t in single precision for backwards compatibility
!
      lread_add = .true.
      if (present (mode)) lread_add = (mode == 1)
!
      if (lfirst_proc_xy) then
        allocate (ga(mxgrid,mygrid,mz,nv), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('input_snap', 'Could not allocate memory for ga', .true.)
        if (ip <= 8) print *, 'input_snap: open ', file
!
        inquire (IOLENGTH=io_len) t_sp
        rec_len = int (mxgrid, kind=8) * int (mygrid, kind=8) * mz
        rec_len = rec_len * nv * io_len
        open (lun_input, FILE=trim (directory_snap)//'/'//file, access='direct', recl=rec_len, status='old')
!
        if (ip <= 8) print *, 'input_snap: read dim=', mxgrid, mygrid, mz, nv
        ! iterate through xy-leading processors in the z-direction
        read (lun_input, rec=1) ga
        ! distribute data in the xy-plane of the pz-layer
        call localize_xy (a, ga)
        deallocate (ga)
!
      else
        ! receive data from root processor
        call localize_xy (a)
      endif
!
      ! read additional data
      if (lread_add) then
        if (lroot) then
          allocate (gx(mxgrid), gy(mygrid), gz(mzgrid), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('input_snap', 'Could not allocate memory for gx,gy,gz', .true.)
!
          rec_len = int (mxgrid, kind=8) * int (mygrid, kind=8)
          num_rec = int (mz, kind=8) * int (nv*sizeof_real(), kind=8)
          close (lun_input)
          open (lun_input, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', status='old')
          call fseek_pos (lun_input, rec_len, num_rec, 0)
          read (lun_input) t_sp
          read (lun_input) gx, gy, gz, dx, dy, dz
          call distribute_grid (x, y, z, gx, gy, gz)
          deallocate (gx, gy, gz)
        else
          call distribute_grid (x, y, z)
!
          if (lfirst_proc_xy) then
            rec_len = int (mxgrid, kind=8) * int (mygrid, kind=8)
            num_rec = int (mz, kind=8) * int (nv*sizeof_real(), kind=8)
            close (lun_input)
            open (lun_input, FILE=trim (directory_snap)//'/'//file, FORM='unformatted', status='old')
            call fseek_pos (lun_input, rec_len, num_rec, 0)
            read (lun_input) t_sp
          endif
        endif
        if (lfirst_proc_xy) t_test = t_sp
        call mpibcast_real (t_sp)
        if (.not. lfirst_proc_xy) t_test = t_sp
        if (t_test /= t_sp) &
            write (*,*) 'ERROR: '//trim(directory_snap)//'/'//trim(file)//' IS INCONSISTENT: t=', t_sp
        call stop_it_if_any ((t_test /= t_sp), '')
        t = t_sp
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
      if (persist_initialized) then
        if (.not. lfirst_proc_xy) close (lun_input)
        persist_initialized = .false.
        persist_last_id = -max_int
      endif
!
      if (lfirst_proc_xy) close (lun_input)
!
    endsubroutine input_snap_finalize
!***********************************************************************
    logical function init_write_persist(file)
!
!  Initialize writing of persistent data to persistent file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot
!
      character (len=*), intent(in), optional :: file
!
      persist_last_id = -max_int
!
      if (present (file)) then
        if (lfirst_proc_xy) close (lun_output)
        open (lun_output, FILE=trim (directory_dist)//'/'//file, FORM='unformatted', status='replace')
        if (ip <= 9) write (*,*) 'begin persistent block'
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
!  13-Dec-2011/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
!
      write_persist_id = .true.
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
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_logical, mpirecv_logical
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
        if (ip <= 9) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global)
      else
        call mpisend_logical (value, 1, 0, tag_log_0D)
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
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_logical, mpirecv_logical
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
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_int, mpirecv_int
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
        if (ip <= 9) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global)
      else
        call mpisend_int (value, 1, 0, tag_int_0D)
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
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_int, mpirecv_int
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
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_real, mpirecv_real
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
        if (ip <= 9) write (*,*) 'write persistent '//trim (label)
        write (lun_output) global
!
        deallocate (global)
      else
        call mpisend_real (value, 1, 0, tag_real_0D)
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
!  12-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_real, mpirecv_real
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
!  13-Dec-2011/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot
!
      character (len=*), intent(in), optional :: file
!
      if (lroot .and. (ip <= 9)) write (*,*) 'begin persistent block'
      if (present (file)) then
        if (lfirst_proc_xy) close (lun_input)
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
!  17-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpibcast_int
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
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_logical, mpirecv_logical
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
      read_persist_logical_0D = .false.
!
    endfunction read_persist_logical_0D
!***********************************************************************
    logical function read_persist_logical_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_logical, mpirecv_logical
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
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_int, mpirecv_int
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
      read_persist_int_0D = .false.
!
    endfunction read_persist_int_0D
!***********************************************************************
    logical function read_persist_int_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_int, mpirecv_int
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
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_real, mpirecv_real
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
      read_persist_real_0D = .false.
!
    endfunction read_persist_real_0D
!***********************************************************************
    logical function read_persist_real_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      use Mpicomm, only: lroot, mpisend_real, mpirecv_real
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
!  10-Feb-2012/Bourdin.KIS: adapted for collective IO
!
      use Mpicomm, only: lroot
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
        open (lun_output, FILE=trim (directory_collect)//'/'//file, FORM='unformatted', status='replace')
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
!  21-jan-02/wolf: coded
!  15-jun-03/axel: Lx,Ly,Lz are now read in from file (Tony noticed the mistake)
!  10-Feb-2012/Bourdin.KIS: adapted for collective IO
!
      use Mpicomm, only: lroot, mpibcast_int, mpibcast_real
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
!   22-Feb-2012/Bourdin.KIS: adapted from io_dist
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
!   22-Feb-2012/Bourdin.KIS: adapted from io_dist
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
