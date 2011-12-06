! $Id$
!
!!!!!!!!!!!!!!!!!!!!!!
!!   io_dist.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!
!
!  Distributed IO (i.e. each process writes its own file data/procX)
!  07-Nov-2001/wd: Put into separate module, so one can choose
!  alternative IO mechanism.
!
!  The file format written by output() (and used, e.g. in var.dat)
!  consists of the followinig Fortran records:
!    1. data(mx,my,mz,nvar)
!    2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
!  Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
!  for one vector field, 8 for var.dat in the case of MHD with entropy.
!
!  04-nov-11/MR: IOSTAT handling generally introduced
!  16-nov-11/MR: calls to outlog adapted
!
module Io
!
  use Cparam, only: intlen
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'io.h'
!
  interface output              ! Overload the `output' function
    module procedure output_vect
    module procedure output_scal
  endinterface
!
  interface output_pencil        ! Overload the `output_pencil' function
    module procedure output_pencil_vect
    module procedure output_pencil_scal
  endinterface
!
  ! define unique logical unit number for input and output calls
  integer :: lun_input=88
  integer :: lun_output=91
!
  !
  ! Interface to external C function(s).
  ! Need to have two different C functions in order to have F90
  ! interfaces, since a pencil can be either a 1-d or a 2-d array.
  !
!   interface output_penciled_vect_c
!     subroutine output_penciled_vect_c(filename,pencil,&
!                                       ndim,i,iy,iz,t, &
!                                       nx,ny,nz,nghost,fnlen)
!       use Cdata, only: mx
!       real,dimension(mx,*) :: pencil
!       double precision :: t
!       integer :: ndim,i,iy,iz,nx,ny,nz,nghost,fnlen
!       character (len=*) :: filename
!     endsubroutine output_penciled_vect_c
!   endinterface
!   !
!   interface output_penciled_scal_c
!     subroutine output_penciled_scal_c(filename,pencil,&
!                                       ndim,i,iy,iz,t, &
!                                       nx,ny,nz,nghost,fnlen)
!       use Cdata, only: mx
!       real,dimension(mx) :: pencil
!       double precision :: t
!       integer :: ndim,i,iy,iz,nx,ny,nz,nghost,fnlen
!       character (len=*) :: filename
!     endsubroutine output_penciled_scal_c
!   endinterface
  !
  !  Still not possible with the NAG compiler (`No specific match for
  !  reference to generic OUTPUT_PENCILED_SCAL_C')
  !
  external output_penciled_scal_c
  external output_penciled_vect_c
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
      use Cdata, only: datadir,directory,datadir_snap,directory_snap
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
    subroutine input(file,a,nv,mode)
!
!  read snapshot file, possibly with mesh and time (if mode=1)
!  11-apr-97/axel: coded
!
      use Cdata
      use Mpicomm, only: start_serialize,end_serialize,stop_it
!
      character (len=*) :: file
      integer :: nv,mode
      real, dimension (mx,my,mz,nv) :: a
      real :: t_sp   ! t in single precision for backwards compatibility
!
      integer :: iostat
!
      if (lserial_io) call start_serialize()
      open(lun_input,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (iostat /= 0) call stop_it("Cannot open "//trim(file)//" for reading.", iostat)
!
      if (ip<=8) print*,'input: open, mx,my,mz,nv=',mx,my,mz,nv
      read(lun_input,IOSTAT=iostat) a
      if (iostat /= 0) call stop_it( "Cannot read a from "//trim(file), iostat)
!
      if (ip<=8) print*,'input: read ',file
      if (mode==1) then
!
!  check whether we want to read deltay from snapshot
!
        if (lshear) then
          read(lun_input,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz,deltay
          if (iostat /= 0) call stop_it( &
            "Cannot read t_sp,x,y,z,dx,dy,dz,deltay from "//trim(file), iostat)
        else
          read(lun_input,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz
          if (iostat /= 0) call stop_it( &
            "Cannot read t_sp,x,y,z,dx,dy,dz from "//trim(file), iostat)
        endif
        t = t_sp
!
        if (ip<=3) print*,'input: ip,x=',ip,x
        if (ip<=3) print*,'input: y=',y
        if (ip<=3) print*,'input: z=',z
!
      endif
!
      close(lun_input,IOSTAT=iostat)
      if (outlog(iostat,'close',file)) continue
!
      if (lserial_io) call end_serialize()
!
    endsubroutine input
!***********************************************************************
    subroutine output_vect(file,a,nv)
!
!  write snapshot file, always write time and mesh, could add other things
!  version for vector field
!
!  11-apr-97/axel: coded
!
      use Cdata
      use Mpicomm, only: start_serialize,end_serialize
!
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
      real :: t_sp   ! t in single precision for backwards compatibility
!
      integer :: iostat
!
      t_sp = t
      if (ip<=8.and.lroot) print*,'output_vect: nv =', nv
!
      if (lserial_io) call start_serialize()
!
      open(lun_output,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (outlog(iostat,'open',file)) goto 99
!
      write(lun_output,IOSTAT=iostat) a
      if (outlog(iostat,'write a')) goto 99
!
      if (lshear) then
        write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz,deltay
        if (outlog(iostat,'write t_sp,x,y,z,dx,dy,dz,deltay')) goto 99
      else
        write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz
        if (outlog(iostat,'write t_sp,x,y,z,dx,dy,dz')) goto 99
      endif
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,'close')) continue
!
99    if (lserial_io) call end_serialize()
!
    endsubroutine output_vect
!***********************************************************************
    subroutine output_scal(file,a,nv)
!
!  write snapshot file, always write time and mesh, could add other things
!  version for scalar field
!
!  11-apr-97/axel: coded
!
      use Cdata
      use Mpicomm, only: stop_it,start_serialize,end_serialize
!
      integer :: nv
      real, dimension (mx,my,mz) :: a
      character (len=*) :: file
      real :: t_sp   ! t in single precision for backwards compatibility
!
      integer :: iostat
!
      t_sp = t
      if ((ip<=8) .and. lroot) print*,'output_scal'
      if (nv /= 1) call stop_it("output_scal: called with scalar field, but nv/=1")
!
      if (lserial_io) call start_serialize()
!
      open(lun_output,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (outlog(iostat,'open',file)) goto 99
!
      write(lun_output,IOSTAT=iostat) a
      if (outlog(iostat,'write a')) goto 99
!
      if (lshear) then
        write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz,deltay
        if (outlog(iostat,'write t_sp,x,y,z,dx,dy,dz,deltay')) goto 99
      else
        write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz
        if (outlog(iostat,'write t_sp,x,y,z,dx,dy,dz')) goto 99
      endif
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,'close')) continue
!
99    if (lserial_io) call end_serialize()
!
    endsubroutine output_scal
!***********************************************************************
    subroutine output_pencil_vect(file,a,ndim)
!
!  Write snapshot file of penciled vector data (for debugging).
!  Wrapper to the C routine output_penciled_vect_c.
!
!  15-feb-02/wolf: coded
!
      use Cdata
!
      integer :: ndim
      real, dimension (nx,ndim) :: a
      character (len=*) :: file
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = t
      if (ip<9.and.lroot.and.imn==1) &
           print*,'output_pencil_vect('//file//'): ndim=',ndim
!
      if (headt .and. (imn==1)) write(*,'(A)') &
           'output_pencil: Writing to ' // trim(file) // &
           ' for debugging -- this may slow things down'
!
       call output_penciled_vect_c(file, a, ndim, &
                                   imn, mm(imn), nn(imn), t_sp, &
                                   nx, ny, nz, nghost, len(file))
!
    endsubroutine output_pencil_vect
!***********************************************************************
    subroutine output_pencil_scal(file,a,ndim)
!
!  Write snapshot file of penciled scalar data (for debugging).
!  Wrapper to the C routine output_penciled_scal_c.
!
!  15-feb-02/wolf: coded
!
      use Cdata
      use Mpicomm, only: stop_it
!
      integer :: ndim
      real, dimension (nx) :: a
      character (len=*) :: file
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = t
      if ((ip<=8) .and. lroot .and. imn==1) &
           print*,'output_pencil_scal('//file//')'
!
      if (ndim /= 1) &
           call stop_it("output_pencil_scal: called with scalar field, but ndim/=1")
!
      if (headt .and. (imn==1)) print*, &
           'output_pencil_scal: Writing to ', trim(file), &
           ' for debugging -- this may slow things down'
!
      call output_penciled_scal_c(file, a, ndim, &
                                  imn, mm(imn), nn(imn), t_sp, &
                                  nx, ny, nz, nghost, len(file))
!
    endsubroutine output_pencil_scal
!***********************************************************************
    subroutine outpus(file,a,nv)
!
!  write snapshot file, always write mesh and time, could add other things
!
!  11-oct-98/axel: adapted
!
      use Cdata
!
      integer :: nv
      character (len=*) :: file
      real, dimension (mx,my,mz,nv) :: a
      real :: t_sp   ! t in single precision for backwards compatibility
!
      integer :: iostat
!
      t_sp = t
      open(lun_output,FILE=file,FORM='unformatted',IOSTAT=iostat)
      if (outlog(iostat,'open',file)) return
!
      write(lun_output,IOSTAT=iostat) a(l1:l2,m1:m2,n1:n2,:)
      if (outlog(iostat,"write a")) return
!
      write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz,deltay
      if (outlog(iostat,"write t_sp,x,y,z,dx,dy,dz,deltay")) return
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,"close")) return
!
    endsubroutine outpus
!***********************************************************************
    subroutine log_filename_to_file(filename,flist)
!
!  In the directory containing `filename', append one line to file
!  `flist' containing the file part of filename
!
      use Cdata, only: lroot,lcopysnapshots_exp,datadir
      use Cparam, only: fnlen
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
      use Cdata
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
      use Cdata
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
      call keep_compiler_quiet(t_sp)
!
    endsubroutine rgrid
!***********************************************************************
    subroutine wproc_bounds(file)
!
!   Export processor boundaries to file.
!
!   10-jul-08/kapelrud: coded
!
      use Cdata, only: procx_bounds,procy_bounds,procz_bounds
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
      use Cdata, only: procx_bounds,procy_bounds,procz_bounds
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
      call keep_compiler_quiet(file)
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
      call keep_compiler_quiet(file)
!
    endsubroutine rtime
!***********************************************************************
endmodule Io
