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
  use Cparam, only: intlen, fnlen
  use Messages
!
  implicit none
!
  include 'io.h'
!
  ! define unique logical unit number for input and output calls
  integer :: lun_input=88
  integer :: lun_output=91
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
!  28-jun-10/julien: Added different file formats
!
      use Mpicomm, only: start_serialize, end_serialize
      use Persist, only: output_persistent
!
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
      real :: t_sp   ! t in single precision for backwards compatibility
      integer :: iostat
!
      t_sp = t
      if (ip<=8.and.lroot) print*,'output_vect: nv =', nv
!
      if (lserial_io) call start_serialize()
      open(lun_output,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (outlog(iostat,'open',file=file,dist=lun_output)) goto 99 
!
      if (lwrite_2d) then
        if (nx==1) then
          write(lun_output,IOSTAT=iostat) a(l1,:,:,:)
        elseif (ny==1) then
          write(lun_output,IOSTAT=iostat) a(:,m1,:,:)
        elseif (nz==1) then
          write(lun_output,IOSTAT=iostat) a(:,:,n1,:)
        else
          iostat=0
          call fatal_error('output_snap','lwrite_2d used for 3-D simulation!')
        endif
      else
        write(lun_output,IOSTAT=iostat) a
      endif
      if (outlog(iostat,'write a')) goto 99          ! problematic if fatal_error doesn't stop program
!
!  Write shear at the end of x,y,z,dx,dy,dz.
!  At some good moment we may want to treat deltay like with
!  other modules and call a corresponding i/o parameter module.
!
      if (lshear) then
        write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz,deltay
        if (outlog(iostat,'write t_sp,x,y,z,dx,dy,dz,deltay')) goto 99
      else
        write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz
        if (outlog(iostat,'write t_sp,x,y,z,dx,dy,dz')) goto 99
      endif
      call output_persistent(lun_output)
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,'close')) continue
!
99    if (lformat) call output_snap_form (file,a,nv)
!
      if (ltec) call output_snap_tec (file,a,nv)
!
      if (lserial_io) call end_serialize()
!
    endsubroutine output_snap
!***********************************************************************
    subroutine input_snap(file,a,nv,mode)
!
!  Read snapshot file, possibly with mesh and time (if mode=1).
!
!  11-apr-97/axel: coded
!
      use Mpicomm, only: start_serialize, end_serialize, stop_it
      use Persist, only: input_persistent
!
      character (len=*) :: file
      integer :: nv,mode
      real, dimension (mx,my,mz,nv) :: a
      real :: t_sp   ! t in single precision for backwards compatibility
      integer :: iostat
!
      if (lserial_io) call start_serialize()
      open(lun_input,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (iostat /= 0) call stop_it("Cannot open "//trim(file)//" for reading",iostat)
!      if (ip<=8) print*,'input_snap: open, mx,my,mz,nv=',mx,my,mz,nv
      if (lwrite_2d) then
        if (nx==1) then
          read(lun_input,IOSTAT=iostat) a(4,:,:,:)
        elseif (ny==1) then
          read(lun_input,IOSTAT=iostat) a(:,4,:,:)
        elseif (nz==1) then
          read(lun_input,IOSTAT=iostat) a(:,:,4,:)
        else
          iostat=0
          call fatal_error('input_snap','lwrite_2d used for 3-D simulation!')
        endif
      else
        read(lun_input,IOSTAT=iostat) a
      endif
      if (iostat /= 0) call stop_it("Cannot read a from "//trim(file),iostat)

      if (ip<=8) print*,'input_snap: read ',file
      if (mode==1) then
!
!  Check whether we want to read deltay from snapshot.
!
        if (lshear) then
          read(lun_input,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz,deltay
          if (iostat /= 0) call stop_it("Cannot read t_sp,x,y,z,dx,dy,dz,deltay from "//trim(file),iostat)
        else
          read(lun_input,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz
          if (iostat /= 0) call stop_it("Cannot read t_sp,x,y,z,dx,dy,dz from "//trim(file),iostat)
        endif
!
!  set initial time to that of snapshot, unless
!  this is overridden
!
        if (lreset_tstart) then
          t=tstart
        else
          t=t_sp
        endif
!
!  verify the ip, x, y, and z readings
!
        if (ip<=3) print*,'input_snap: ip,x=',ip,x
        if (ip<=3) print*,'input_snap: y=',y
        if (ip<=3) print*,'input_snap: z=',z
!
      endif
!
      call input_persistent(lun_input)
!
      close(lun_input,IOSTAT=iostat)
      if (outlog(iostat,'close',file)) return
!
      if (lserial_io) call end_serialize()
!
    endsubroutine input_snap
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
    subroutine output_snap_form(file,a,nv)
!
!  Write FORMATTED snapshot file
!
!  28-june-10/julien: coded (copy from output_snap)
!
      integer :: nv
      integer :: i, j, k
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
!
      integer :: iostat

      open(lun_output,FILE=trim(file)//'.form',IOSTAT=iostat)
      if (outlog(iostat,'open',trim(file)//'.form',dist=lun_output)) return 
!
      if (lwrite_2d) then
!
        if (nx==1) then
          do i = m1, m2
            do j = n1, n2
              write(lun_output,'(40(f12.5))',IOSTAT=iostat) x(l1),y(i),z(j),dx,dy,dz,a(l1,i,j,:)
              if (outlog(iostat,'write x, y, z, dx, dy, dz, a(l1,i,j,:)')) return
            enddo
          enddo
        elseif (ny==1) then
          do i = l1, l2
            do j = n1, n2
              write(lun_output,'(40(f12.5))',IOSTAT=iostat) x(i),y(m1),z(j),dx,dy,dz,a(i,m1,j,:)
              if (outlog(iostat,'write x, y, z, dx, dy, dz, a(i,m1,j,:)')) return
            enddo
          enddo
        elseif (nz==1) then
          do i = l1, l2
            do j = m1, m2
              write(lun_output,'(40(f12.5))',IOSTAT=iostat) x(i),y(j),z(n1),dx,dy,dz,a(i,j,n1,:)
              if (outlog(iostat,'write x, y, z, dx, dy, dz, a(i,j,n1,:)')) return
            enddo
          enddo
        else
          call fatal_error('output_snap','lwrite_2d used for 3-D simulation!')
        endif
!
      else if (ny==1.and.nz==1) then
!
        do i = l1, l2
          write(lun_output,'(40(f12.5))',IOSTAT=iostat) x(i),a(i,m1,n1,:)
          if (outlog(iostat,'write x, a')) return
        enddo
!
      else
!
        do i = l1, l2
          do j = m1, m2
            do k = n1, n2
              write(lun_output,'(40(f12.5))',IOSTAT=iostat) x(i),y(j),z(k),dx,dy,dz,a(i,j,k,:)
              if (outlog(iostat,'write x, y, z, dx, dy, dz, a')) return
            enddo
          enddo
        enddo
!
      endif
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,'close')) continue
!
    endsubroutine output_snap_form
!***********************************************************************
    subroutine output_snap_tec(file,a,nv)
!
!  Write TECPLOT output files (binary)
!
!  28-june-10/julien: coded
!
      integer :: nv
      integer :: i, j, k, kk, iostat
      real, dimension (mx,my,mz,nv) :: a
      real, dimension (nx*ny*nz) :: xx, yy, zz
      character (len=*) :: file
      character(len=2) :: car
      character (len=8), dimension (nv) :: name
      character (len=120) :: filel
!
      filel=trim(file)//'.tec'

      open(lun_output,FILE=filel,IOSTAT=iostat)
      if (outlog(iostat,'open',filel,dist=lun_output)) return 
!
      kk = 0
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            xx(kk+i) = x(i)
            yy(kk+i) = y(j)
            zz(kk+i) = z(k)
          enddo
          kk = kk + nx
        enddo
      enddo
!
!  Write header
!
      write(lun_output,*,IOSTAT=iostat) 'TITLE     = "output"'
      if (outlog(iostat,'write TITLE')) return
!
      if (lwrite_2d) then
!
        if (nx==1) then
          write(lun_output,*,IOSTAT=iostat) 'VARIABLES = "y"'
          if (outlog(iostat,'write "VARIABLES = y"')) return
          write(lun_output,*,IOSTAT=iostat) '"z"'
          if (outlog(iostat,'write "z"')) return
        elseif (ny==1) then
          write(lun_output,*,IOSTAT=iostat) 'VARIABLES = "x"'
          if (outlog(iostat,'write "VARIABLES = x"')) return
          write(lun_output,*,IOSTAT=iostat) '"z"'
          if (outlog(iostat,'write "z"')) return
        elseif (nz==1) then
          write(lun_output,*,IOSTAT=iostat) 'VARIABLES = "x"'
          if (outlog(iostat,'write "VARIABLES = x"')) return
          write(lun_output,*,IOSTAT=iostat) '"y"'
          if (outlog(iostat,'write "y"')) return
        endif
!
      else
!
        if (ny==1.and.nz==1) then
          write(lun_output,*,IOSTAT=iostat) 'VARIABLES = "x"'
          if (outlog(iostat,'write "VARIABLES = x"')) return
        else
          write(lun_output,*,IOSTAT=iostat) 'VARIABLES = "x"'
          if (outlog(iostat,'write "VARIABLES = x"')) return
          write(lun_output,*,IOSTAT=iostat) '"y"'
          if (outlog(iostat,'write "VARIABLES = y"')) return
          write(lun_output,*,IOSTAT=iostat) '"z"'
          if (outlog(iostat,'write "VARIABLES = z"')) return
        endif
!
      endif
      do i = 1, nv
        write(car,'(i2)') i
        name(i) = 'VAR_'//adjustl(car)
        write(lun_output,*,IOSTAT=iostat) '"'//trim(name(i))//'"'
        if (outlog(iostat,'write name(i)')) return
      enddo
!
      write(lun_output,*,IOSTAT=iostat) 'ZONE T="Zone"'
      if (outlog(iostat,'write "ZONE T=Zone"')) return
!
      if (lwrite_2d) then
        if (nx==1) then
          write(lun_output,*,IOSTAT=iostat) ' I=1, J=',ny, ', K=',nz
          if (outlog(iostat,'write ny, nz')) return
        endif
        if (ny==1) then
          write(lun_output,*,IOSTAT=iostat) ' I=',nx, ', J=1, K=',nz
          if (outlog(iostat,'write nx, nz')) return
        endif
        if (nz==1) then
          write(lun_output,*,IOSTAT=iostat) ' I=',nx, ', J=',ny, ', K=1'
          if (outlog(iostat,'write nx, ny')) return
        endif
      else
        if (ny==1.and.nz==1) then
          write(lun_output,*,IOSTAT=iostat) ' I=',nx, ', J=  1, K=  1'
          if (outlog(iostat,'write nx')) return
        else
          write(lun_output,*,IOSTAT=iostat) ' I=',nx, ', J=',ny, ', K=',nz
          if (outlog(iostat,'write nx, ny, nz')) return
        endif
      endif
!
      write(lun_output,*,IOSTAT=iostat) ' DATAPACKING=BLOCK'
      if (outlog(iostat,'write "DATAPACKING=BLOCK"')) return
!
!
!  Write data
!
      if (lwrite_2d) then
        if (nx==1) then
!
          write(lun_output,*,IOSTAT=iostat) yy
          if (outlog(iostat,'write yy')) return
          write(lun_output,*,IOSTAT=iostat) zz
          if (outlog(iostat,'write zz')) return
!
          do j = 1, nv
            write(lun_output,*,IOSTAT=iostat) a(l1,m1:m2,n1:n2,j)
            if (outlog(iostat,'write a')) return
          enddo
!
        elseif (ny==1) then
!
          write(lun_output,*,IOSTAT=iostat) xx
          if (outlog(iostat,'write xx')) return
!
          write(lun_output,*,IOSTAT=iostat) zz
          if (outlog(iostat,'write zz')) return
!
          do j = 1, nv
            write(lun_output,*,IOSTAT=iostat) a(l1:l2,m1,n1:n2,j)
            if (outlog(iostat,'write a')) return
          enddo
!
        elseif (nz==1) then
          write(lun_output,*,IOSTAT=iostat) xx
          if (outlog(iostat,'write xx')) return
!
          write(lun_output,*,IOSTAT=iostat) yy
          if (outlog(iostat,'write yy')) return
!
          do j = 1, nv
            write(lun_output,*,IOSTAT=iostat) a(l1:l2,m1:m2,n1,j)
            if (outlog(iostat,'write a')) return
          enddo
!
        else
          call fatal_error('output_snap','lwrite_2d used for 3-D simulation!')
        endif
      else if (ny==1.and.nz==1) then
!
             write(lun_output,*,IOSTAT=iostat) xx
             if (outlog(iostat,'write xx')) return
!
             do j = 1, nv
               write(lun_output,*,IOSTAT=iostat) a(l1:l2,m1,n1,j)
               if (outlog(iostat,'write a')) return
             enddo
!
           else
!
             write(lun_output,*,IOSTAT=iostat) xx
             if (outlog(iostat,'write xx')) return
!
             write(lun_output,*,IOSTAT=iostat) yy
             if (outlog(iostat,'write yy')) return
!
             write(lun_output,*,IOSTAT=iostat) zz
             if (outlog(iostat,'write zz')) return
!
             do j = 1, nv
               write(lun_output,*,IOSTAT=iostat) a(l1:l2,m1:m2,n1:n2,j)
               if (outlog(iostat,'write a')) return
             enddo
!
           endif
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,'close')) continue
!
    endsubroutine output_snap_tec
!***********************************************************************
endmodule Io
