! $Id$
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   debug_io_dist.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Distributed debug-IO (i.e. each process writes its own file data/procX)
!
!  The file format written by output() (and used, e.g. in var.dat)
!  consists of the followinig Fortran records:
!    1. data(mx,my,mz,nvar)
!    2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1)
! or 2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
!  Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
!  for one vector field, 8 for var.dat in the case of MHD with entropy.
!
!  11-Dec-2011/Bourdin.KIS: moved debug-IO from 'io_dist' into separate module
!
module Debug_IO
!
  use Cdata
!
  implicit none
!
  public :: lun_input, lun_output
  public :: output, output_pencil
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
  integer :: lun_input=89
  integer :: lun_output=92
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
    subroutine output_vect(file,a,nv)
!
!  write snapshot file, always write time and mesh, could add other things
!  version for vector field
!
!  11-apr-97/axel: coded
!
      use Messages, only: outlog
      use Mpicomm, only: start_serialize, end_serialize
!
      integer :: nv
      real, dimension (mx,my,mz,nv) :: a
      character (len=*) :: file
      real :: t_sp   ! t in single precision for backwards compatibility
!
      integer :: iostat
      character (len=fnlen) :: filename
!
      t_sp = t
      if (ip<=8.and.lroot) print*,'output_vect: nv =', nv
      filename = trim (datadir_snap)//'/'//trim (file)
!
      if (lserial_io) call start_serialize()
!
      open(lun_output,FILE=filename,FORM='unformatted',IOSTAT=iostat)
      if (outlog(iostat,'open',filename)) then
        if (lserial_io) call end_serialize()
        return
      endif
!
      write(lun_output,IOSTAT=iostat) a
      if (outlog(iostat,'write a')) then
        if (lserial_io) call end_serialize()
        return
      endif
!
      if (lshear) then
        write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz,deltay
        if (outlog(iostat,'write t_sp,x,y,z,dx,dy,dz,deltay')) then
          if (lserial_io) call end_serialize()
          return
        endif
      else
        write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz
        if (outlog(iostat,'write t_sp,x,y,z,dx,dy,dz')) then
          if (lserial_io) call end_serialize()
          return
        endif
      endif
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,'close')) continue
!
      if (lserial_io) call end_serialize()
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
      use Messages, only: outlog
      use Mpicomm, only: stop_it, start_serialize, end_serialize
!
      integer :: nv
      real, dimension (mx,my,mz) :: a
      character (len=*) :: file
      real :: t_sp   ! t in single precision for backwards compatibility
!
      integer :: iostat
      character (len=fnlen) :: filename
!
      t_sp = t
      if ((ip<=8) .and. lroot) print*,'output_scal'
      if (nv /= 1) call stop_it("output_scal: called with scalar field, but nv/=1")
      filename = trim (datadir_snap)//'/'//trim (file)
!
      if (lserial_io) call start_serialize()
!
      open(lun_output,FILE=filename,FORM='unformatted',IOSTAT=iostat)
      if (outlog(iostat,'open',filename)) then
        if (lserial_io) call end_serialize()
        return
      endif
!
      write(lun_output,IOSTAT=iostat) a
      if (outlog(iostat,'write a')) then
        if (lserial_io) call end_serialize()
        return
      endif
!
      if (lshear) then
        write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz,deltay
        if (outlog(iostat,'write t_sp,x,y,z,dx,dy,dz,deltay')) then
          if (lserial_io) call end_serialize()
          return
        endif
      else
        write(lun_output,IOSTAT=iostat) t_sp,x,y,z,dx,dy,dz
        if (outlog(iostat,'write t_sp,x,y,z,dx,dy,dz')) then
          if (lserial_io) call end_serialize()
          return
        endif
      endif
!
      close(lun_output,IOSTAT=iostat)
      if (outlog(iostat,'close')) continue
!
      if (lserial_io) call end_serialize()
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
      character (len=*), intent(in) :: file
      integer, intent(in) :: ndim
      real, dimension (nx,ndim), intent(in) :: a
!
      real :: t_sp   ! t in single precision for backwards compatibility
      character (len=fnlen) :: filename
!
      filename = trim (datadir_snap)//'/'//trim (file)
      t_sp = t
      if (ip<9.and.lroot.and.imn==1) &
           print*,'output_pencil_vect('//filename//'): ndim=',ndim
!
      if (headt .and. (imn==1)) write(*,'(A)') &
           'output_pencil: Writing to ' // trim(filename) // &
           ' for debugging -- this may slow things down'
!
       call output_penciled_vect_c(filename, a, ndim, &
                                   imn, mm(imn), nn(imn), t_sp, &
                                   nx, ny, nz, nghost, len(filename))
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
      use Mpicomm, only: stop_it
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: ndim
      real, dimension (nx), intent(in) :: a
!
      real :: t_sp   ! t in single precision for backwards compatibility
      character (len=fnlen) :: filename
!
      t_sp = t
      filename = trim (datadir_snap)//'/'//trim (file)
      if ((ip<=8) .and. lroot .and. imn==1) &
           print*,'output_pencil_scal('//trim(filename)//')'
!
      if (ndim /= 1) &
           call stop_it("output_pencil_scal: called with scalar field, but ndim/=1")
!
      if (headt .and. (imn==1)) print*, &
           'output_pencil_scal: Writing to ', trim(filename), &
           ' for debugging -- this may slow things down'
!
      call output_penciled_scal_c(filename, a, ndim, &
                                  imn, mm(imn), nn(imn), t_sp, &
                                  nx, ny, nz, nghost, len(filename))
!
    endsubroutine output_pencil_scal
!***********************************************************************
endmodule Debug_IO
