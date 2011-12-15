! $Id$
!
!!!!!!!!!!!!!!!!!!!!!!
!!   io_dist.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!
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
!
module Io
!
  use Cdata
  use Cparam, only: intlen, fnlen, max_int
  use Messages
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
  logical :: lcollective_IO=.false.
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
!  20-sep-02/wolf: coded
!
      use Mpicomm, only: lroot
!
!  identify version number
!
      if (lroot) call svn_id("$Id$")
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
      use Mpicomm, only: iproc
      use General, only: itoa, safe_character_assign
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
      chproc=itoa(iproc)
      call safe_character_assign(directory, trim(datadir)//'/proc'//chproc)
      call safe_character_assign(directory_snap, &
                                            trim(datadir_snap)//'/proc'//chproc)
!
    endsubroutine directory_names
!***********************************************************************
    subroutine output_snap(file,a,nv)
!
!  Write snapshot file, always write time and mesh, could add other things
!  version for vector field.
!
!  11-apr-97/axel: coded
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use Mpicomm, only: start_serialize, end_serialize
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: a
!
      real :: t_sp   ! t in single precision for backwards compatibility
      integer :: ierr
!
      t_sp = t
      if (lroot .and. (ip <= 8)) print *, 'output_vect: nv =', nv
!
      if (lserial_io) call start_serialize()
      open (lun_output, FILE=file, FORM='unformatted', IOSTAT=ierr)
      if (outlog (ierr, 'open', file=file, dist=lun_output)) then
        if (lserial_io) call end_serialize()
        return
      endif
!
      if (lwrite_2d) then
        if (nx == 1) then
          write (lun_output, IOSTAT=ierr) a(l1,:,:,:)
        elseif (ny == 1) then
          write (lun_output, IOSTAT=ierr) a(:,m1,:,:)
        elseif (nz == 1) then
          write (lun_output, IOSTAT=ierr) a(:,:,n1,:)
        else
          ierr = 0
          call fatal_error ('output_snap', 'lwrite_2d used for 3D simulation!')
        endif
      else
        write (lun_output, IOSTAT=ierr) a
      endif
      ! problematic if fatal_error doesn't stop program:
      if (outlog (ierr, 'write a')) then
        if (lserial_io) call end_serialize()
        return
      endif
!
!  Write shear at the end of x,y,z,dx,dy,dz.
!  At some good moment we may want to treat deltay like with
!  other modules and call a corresponding i/o parameter module.
!
      if (lshear) then
        write (lun_output, IOSTAT=ierr) t_sp, x, y, z, dx, dy, dz, deltay
        if (outlog (ierr, 'write t_sp,x,y,z,dx,dy,dz,deltay')) then
          if (lserial_io) call end_serialize()
          return
        endif
      else
        write (lun_output, IOSTAT=ierr) t_sp, x, y, z, dx, dy, dz
        if (outlog (ierr, 'write t_sp,x,y,z,dx,dy,dz')) then
          if (lserial_io) call end_serialize()
          return
        endif
      endif
!
    endsubroutine output_snap
!***********************************************************************
    subroutine output_snap_finalize()
!
!  Close snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: adapted from output_snap
!
      use Mpicomm, only: end_serialize
!
      integer :: ierr
      logical :: lerror
!
      if (persist_initialized) then
        write (lun_output, iostat=ierr) id_block_PERSISTENT
        lerror = outlog (ierr, 'write id_block_PERSISTENT')
        persist_initialized = .false.
        persist_last_id = -max_int
      endif
!
      close (lun_output, IOSTAT=ierr)
      lerror = outlog (ierr, 'close')
!
      if (lserial_io) call end_serialize()
!
    endsubroutine output_snap_finalize
!***********************************************************************
    logical function init_write_persist()
!
!  Initialize writing of persistent data to snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      integer :: ierr
!
      write (lun_output, iostat=ierr) id_block_PERSISTENT
      init_write_persist = outlog (ierr, 'write id_block_PERSISTENT')
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
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
!
      integer :: ierr
!
      write_persist_id = .true.
      if (.not. persist_initialized) return
!
      if (persist_last_id /= id) then
        write (lun_output, iostat=ierr) id
        write_persist_id = outlog (ierr, 'write persistent ID '//label)
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
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: id
      logical, intent(in) :: value
!
      integer :: ierr
!
      write_persist_logical_0D = .true.
      if (.not. persist_initialized) return
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=ierr) value
      write_persist_logical_0D = outlog (ierr, 'write persistent '//label)
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
      integer :: ierr
!
      write_persist_logical_1D = .true.
      if (.not. persist_initialized) return
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=ierr) value
      write_persist_logical_1D = outlog (ierr, 'write persistent '//label)
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
      integer :: ierr
!
      write_persist_int_0D = .true.
      if (.not. persist_initialized) return
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=ierr) value
      write_persist_int_0D = outlog (ierr, 'write persistent '//label)
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
      integer :: ierr
!
      write_persist_int_1D = .true.
      if (.not. persist_initialized) return
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=ierr) value
      write_persist_int_1D = outlog (ierr, 'write persistent '//label)
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
      integer :: ierr
!
      write_persist_real_0D = .true.
      if (.not. persist_initialized) return
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=ierr) value
      write_persist_real_0D = outlog (ierr, 'write persistent '//label)
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
      integer :: ierr
!
      write_persist_real_1D = .true.
      if (.not. persist_initialized) return
      if (write_persist_id (label, id)) return
!
      write (lun_output, iostat=ierr) value
      write_persist_real_1D = outlog (ierr, 'write persistent '//label)
!
    endfunction write_persist_real_1D
!***********************************************************************
    subroutine input_snap(file,a,nv,mode)
!
!  Read snapshot file, possibly with mesh and time (if mode=1).
!
!  11-apr-97/axel: coded
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use Mpicomm, only: start_serialize, end_serialize, stop_it
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real, dimension (mx,my,mz,nv), intent(out) :: a
!
      real :: t_sp   ! t in single precision for backwards compatibility
      integer :: ierr
!
      if (lserial_io) call start_serialize()
      open (lun_input, FILE=file, FORM='unformatted', IOSTAT=ierr)
      if (ierr /= 0) call stop_it ("Cannot open "//trim(file)//" for reading", ierr)
!      if (ip<=8) print *, 'input_snap: open, mx,my,mz,nv=', mx, my, mz, nv
      if (lwrite_2d) then
        if (nx == 1) then
          read (lun_input, IOSTAT=ierr) a(4,:,:,:)
        elseif (ny == 1) then
          read (lun_input, IOSTAT=ierr) a(:,4,:,:)
        elseif (nz == 1) then
          read (lun_input, IOSTAT=ierr) a(:,:,4,:)
        else
          ierr = 0
          call fatal_error ('input_snap', 'lwrite_2d used for 3-D simulation!')
        endif
      else
        read (lun_input, IOSTAT=ierr) a
      endif
      if (ierr /= 0) call stop_it ("Cannot read a from "//trim(file), ierr)

      if (ip <= 8) print *, 'input_snap: read ', file
      if (mode == 1) then
!
!  Check whether we want to read deltay from snapshot.
!
        if (lshear) then
          read (lun_input, IOSTAT=ierr) t_sp, x, y, z, dx, dy, dz, deltay
          if (ierr /= 0) call stop_it ("Cannot read t_sp,x,y,z,dx,dy,dz,deltay from "//trim(file), ierr)
        else
          read (lun_input, IOSTAT=ierr) t_sp, x, y, z, dx, dy, dz
          if (ierr /= 0) call stop_it ("Cannot read t_sp,x,y,z,dx,dy,dz from "//trim(file), ierr)
        endif
!
!  set initial time to that of snapshot, unless
!  this is overridden
!
        if (lreset_tstart) then
          t = tstart
        else
          t = t_sp
        endif
!
!  verify the ip, x, y, and z readings
!
        if (ip <= 3) print *, 'input_snap: ip,x=', ip, x
        if (ip <= 3) print *, 'input_snap: y=', y
        if (ip <= 3) print *, 'input_snap: z=', z
!
      endif
!
    endsubroutine input_snap
!***********************************************************************
    subroutine input_snap_finalize(file)
!
!  Close snapshot file.
!
!  11-apr-97/axel: coded
!  13-Dec-2011/Bourdin.KIS: reworked
!
      use Mpicomm, only: end_serialize
!
      character (len=*), intent(in) :: file
!
      integer :: ierr
!
      close (lun_input, IOSTAT=ierr)
      if (outlog (ierr, 'close', file)) return
!
      if (lserial_io) call end_serialize()
!
    endsubroutine input_snap_finalize
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
      integer :: ierr
!
      read (lun_input, iostat=ierr) value
      read_persist_logical_0D = outlog (ierr, 'read persistent '//label)
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
      integer :: ierr
!
      read (lun_input, iostat=ierr) value
      read_persist_logical_1D = outlog (ierr, 'read persistent '//label)
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
      integer :: ierr
!
      read (lun_input, iostat=ierr) value
      read_persist_int_0D = outlog (ierr, 'read persistent '//label)
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
      integer :: ierr
!
      read (lun_input, iostat=ierr) value
      read_persist_int_1D = outlog (ierr, 'read persistent '//label)
!
    endfunction read_persist_int_1D
!***********************************************************************
    logical function read_persist_real_0D(label, value)
!
!  Read persistent data from snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      real, intent(out) :: value
!
      integer :: ierr
!
      read (lun_input, iostat=ierr) value
      read_persist_real_0D = outlog (ierr, 'read persistent '//label)
!
    endfunction read_persist_real_0D
!***********************************************************************
    logical function read_persist_real_1D(label, value)
!
!  Read persistent data from snapshot file.
!
!  13-Dec-2011/Bourdin.KIS: coded
!
      character (len=*), intent(in) :: label
      real, dimension(:), intent(out) :: value
!
      integer :: ierr
!
      read (lun_input, iostat=ierr) value
      read_persist_real_1D = outlog (ierr, 'read persistent '//label)
!
    endfunction read_persist_real_1D
!***********************************************************************
    subroutine output_globals(file,a,nv)
!
!  Write snapshot file of globals, always write mesh.
!
!  10-nov-06/tony: coded
!
      use Mpicomm, only: start_serialize,end_serialize
!
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
!
      integer :: iostat

      if (ip<=8.and.lroot) print*,'output_vect: nv =', nv
!
      if (lserial_io) call start_serialize()
      open(lun_output,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (outlog(iostat,'open',file)) goto 99
!
      if (lwrite_2d) then
        if (nx==1) then
          write(lun_output,IOSTAT=iostat) a(4,:,:,:)
        elseif (ny==1) then
          write(lun_output,IOSTAT=iostat) a(:,4,:,:)
        elseif (nz==1) then
          write(lun_output,IOSTAT=iostat) a(:,:,4,:)
        else
          iostat=0
          call fatal_error('output_globals','lwrite_2d used for 3-D simulation!')
        endif
      else
        write(lun_output,IOSTAT=iostat) a
      endif
      if (outlog(iostat,'write a')) goto 99                   !problematic if fatal_error doesn't stop program
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,'close')) continue
!
99    if (lserial_io) call end_serialize()
!
    endsubroutine output_globals
!***********************************************************************
    subroutine input_globals(filename,a,nv)
!
!  Read globals snapshot file, ignoring mesh.
!
!  10-nov-06/tony: coded
!
      use Mpicomm, only: start_serialize,end_serialize,stop_it
!
      character (len=*) :: filename
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
!
      integer :: iostat
!
      if (lserial_io) call start_serialize()
!
      open(lun_input,FILE=filename,FORM='unformatted',IOSTAT=iostat)
      if (iostat /= 0) call stop_it("Cannot open "//trim(filename)//" for reading",iostat)

      if (ip<=8) print*,'input_globals: open, mx,my,mz,nv=',mx,my,mz,nv
      if (lwrite_2d) then
        if (nx==1) then
          read(lun_input,IOSTAT=iostat) a(4,:,:,:)
        elseif (ny==1) then
          read(lun_input,IOSTAT=iostat) a(:,4,:,:)
        elseif (nz==1) then
          read(lun_input,IOSTAT=iostat) a(:,:,4,:)
        else
          iostat=0
          call fatal_error('input_globals','lwrite_2d used for 3-D simulation!')
        endif
      else
        read(lun_input,IOSTAT=iostat) a
      endif
      if (iostat /= 0) call stop_it("Cannot read a from "//trim(filename),iostat)
      if (ip<=8) print*,'input_globals: read ',filename
      close(lun_input,IOSTAT=iostat)
      if (outlog(iostat,'close',filename)) continue
!
      if (lserial_io) call end_serialize()
!
    endsubroutine input_globals
!***********************************************************************
    subroutine log_filename_to_file(filename,flist)
!
!  In the directory containing `filename', append one line to file
!  `flist' containing the file part of filename
!
      use General, only: parse_filename
      use Mpicomm, only: mpibarrier
!
      character (len=*) :: filename,flist
      character (len=fnlen) :: dir,fpart
      integer :: iostat
!
      call parse_filename(filename,dir,fpart)
      open(lun_output,FILE=trim(dir)//'/'//trim(flist),POSITION='append',IOSTAT=iostat)
! file not distributed???, backskipping enabled 
      if (outlog(iostat,"open",trim(dir)//'/'//trim(flist),dist=-lun_output)) goto 99
!
      write(lun_output,'(A)',IOSTAT=iostat) trim(fpart)
      if (outlog(iostat,"write fpart")) goto 99
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,"close")) continue
!
99    if (lcopysnapshots_exp) then
        call mpibarrier()
        if (lroot) then
          open(lun_output,FILE=trim(datadir)//'/move-me.list',POSITION='append',IOSTAT=iostat)
! file not distributed, backskipping enabled
          if (outlog(iostat,"open",trim(datadir)//'/move-me.list',dist=-lun_output)) return
!
          write(lun_output,'(A)',IOSTAT=iostat) trim(fpart)
          if (outlog(iostat,"write fpart")) return

          close(lun_output,IOSTAT=iostat)
          if (outlog(iostat,"close")) continue

        endif
      endif
!
    endsubroutine log_filename_to_file
!***********************************************************************
    subroutine wgrid(file)
!
!  Write processor-local part of grid coordinates.
!
!  21-jan-02/wolf: coded
!  15-jun-03/axel: Lx,Ly,Lz are now written to file (Tony noticed the mistake)
!
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
      integer :: iostat
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = t

      open(lun_output,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (iostat /= 0) call stop_it( &
          "Cannot open " // trim(file) // " (or similar) for writing" // &
          " -- is data/ visible from all nodes?", iostat)
      write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz
      write(lun_output,IOSTAT=iostat) dx,dy,dz
      write(lun_output,IOSTAT=iostat) Lx,Ly,Lz
      write(lun_output,IOSTAT=iostat) dx_1,dy_1,dz_1
      write(lun_output,IOSTAT=iostat) dx_tilde,dy_tilde,dz_tilde
      if (iostat /= 0) call stop_it( &
          "Cannot write "//trim(file)//" properly", iostat)
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,'close',file)) continue
!
    endsubroutine wgrid
!***********************************************************************
    subroutine rgrid (file)
!
!  Read processor-local part of grid coordinates.
!
!  21-jan-02/wolf: coded
!  15-jun-03/axel: Lx,Ly,Lz are now read in from file (Tony noticed the mistake)
!   3-jun-04/bing: added xiprim, psiprim ,zetaprim, etc.
!
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
!
      integer :: iostat, ierr
      real :: t_sp   ! t in single precision for backwards compatibility
!
!  if xiprim etc is not written, just ignore it
!
      open(lun_input,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (iostat /= 0) call stop_it( &
          "Cannot open " // trim(file) // " (or similar) for reading" // &
          " -- is data/ visible from all nodes?",iostat)
!
      read(lun_input,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz
      if (iostat/=0) call stop_it("Error when reading t_sp,x,y,z,dx,dy,dz from "//trim(file),iostat) 
!
      read(lun_input,IOSTAT=iostat) dx,dy,dz
      if (iostat/=0) call stop_it("Error when reading dx,dy,dz from "//trim(file),iostat) 
!
      read(lun_input,IOSTAT=ierr) Lx,Ly,Lz     
!      print*, 'Lx,Ly,Lz=', Lx,Ly,Lz
!
      read(lun_input,end=990,IOSTAT=iostat) dx_1,dy_1,dz_1
      if (outlog(iostat,"read dx_1,dy_1,dz_1",file)) continue
!
      read(lun_input,end=990,IOSTAT=iostat) dx_tilde,dy_tilde,dz_tilde
      if (outlog(iostat,"read dx_tilde,dy_tilde,dz_tilde")) continue
!
990   close(lun_input,IOSTAT=iostat)
      if (outlog(iostat,'close')) continue
!
!  give notification if Lx is not read in
!  This should only happen when reading in old data files
!  We should keep this for the time being
!
      if (ierr /= 0) then
        if (ierr < 0) then
          print*,'rgrid: Lx,Ly,Lz are not yet in grid.dat'
        else
          print*, 'rgrid: IOSTAT=', ierr
          call stop_it("rgrid: error when reading Lx,Ly,Lz from "//trim(file),ierr)
        endif
      endif
!
!  Find minimum/maximum grid spacing. Note that
!    minval( (/dx,dy,dz/), MASK=((/nxgrid,nygrid,nzgrid/) > 1) )
!  will be undefined if all n[x-z]grid=1, so we have to add the fourth
!  component with a test that is always true
!
      dxmin = minval( (/dx,dy,dz,huge(dx)/), &
                MASK=((/nxgrid,nygrid,nzgrid,2/) > 1) )
      dxmax = maxval( (/dx,dy,dz,epsilon(dx)/), &
                MASK=((/nxgrid,nygrid,nzgrid,2/) > 1) )
!
!  Fill pencil with maximum gridspacing. Will be overwritten
!  during the mn loop in the non equiditant case
!
      dxmax_pencil(:) = dxmax
      dxmin_pencil(:) = dxmin
!
!  debug output
!
      if (ip<=4.and.lroot) then
        print*,'rgrid: Lx,Ly,Lz=',Lx,Ly,Lz
        print*,'rgrid: dx,dy,dz=',dx,dy,dz
        print*,'rgrid: dxmin,dxmax=',dxmin,dxmax
      endif
!
!  should stop if dxmin=0
!
      if (dxmin==0) call stop_it("rgrid: check Lx,Ly,Lz: is one of them 0?")
!
    endsubroutine rgrid
!***********************************************************************
    subroutine wproc_bounds(file)
!
!   Export processor boundaries to file.
!
!   10-jul-08/kapelrud: coded
!
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
      integer :: iostat
!
      open(lun_output,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (outlog(iostat,'open',file)) return
!
      write(lun_output,IOSTAT=iostat) procx_bounds
      if (outlog(iostat,'write procx_bounds')) return
!
      write(lun_output,IOSTAT=iostat) procy_bounds
      if (outlog(iostat,'write procy_bounds')) return
!
      write(lun_output,IOSTAT=iostat) procz_bounds
      if (outlog(iostat,'write procz_bounds')) return
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,'close' )) continue
!
    endsubroutine wproc_bounds
!***********************************************************************
    subroutine rproc_bounds(file)
!
!   Import processor boundaries from file.
!
!   10-jul-08/kapelrud: coded
!
      use Mpicomm, only: stop_it
!
      character (len=*) :: file
!
      integer :: iostat
!
      open(lun_input,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (iostat/=0) call stop_it("Cannot open "//trim(file)//" for reading",iostat) 
!     
      read(lun_input,IOSTAT=iostat) procx_bounds
      if (iostat/=0) call stop_it("Error when reading procx_bounds from "//trim(file),iostat) 

      read(lun_input,IOSTAT=iostat) procy_bounds
      if (iostat/=0) call stop_it("Error when reading procy_bounds from "//trim(file),iostat) 

      read(lun_input,IOSTAT=iostat) procz_bounds
      if (iostat/=0) call stop_it("Error when reading procz_bounds from "//trim(file),iostat) 
!
      close(lun_input,IOSTAT=iostat)
      if (outlog(iostat,'close',file)) continue
!
    endsubroutine rproc_bounds
!***********************************************************************
    subroutine wtime(file,tau)
!
      double precision :: tau,tmp
      character (len=*) :: file
!
!     nothing needs to be done here
!
! temporary work around to keep the compiler quiet
      tmp = tau
      file = trim (file)
!
    endsubroutine wtime
!***********************************************************************
    subroutine rtime(file,tau)
!
      double precision :: tau,tmp
      character (len=*) :: file
!
!     nothing needs to be done here
!
! temporary work around to keep the compiler quiet
      tmp = tau
      file = trim (file)
!
    endsubroutine rtime
!***********************************************************************
endmodule Io
