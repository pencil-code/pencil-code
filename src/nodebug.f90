! $Id$

!!!!!!!!!!!!!!!!!!!!!!!
!!!   nodebug.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
!!!  Date:   24-Jun-2002
!!!
!!!  Description:
!!!   Two dummy debugging routines (in case the C stuff doesn't work).

!***********************************************************************
subroutine output_penciled_vect_c(filename,pencil,&
                                  ndim,i,iy,iz,t, &
                                  nx,ny,nz,nghost,fnlen)
  use Cparam, only: NO_WARN
  use Cdata, only: ip,mx,headt,imn

  real,dimension(mx,*) :: pencil
  real :: t
  integer :: ndim,i,iy,iz,nx,ny,nz,nghost,fnlen
  character (len=*) :: filename

  if (headt .and. (imn==1)) print*, &
       'OUTPUT_PENCIL: Not writing to ', trim(filename), &
       ' since DEBUG=nodebug'

  ! to keep compiler quiet...
  if (NO_WARN) print*,pencil(1,1),ndim,i,iy,iz,t,nx,ny,nz,nghost,fnlen
!
endsubroutine output_penciled_vect_c
!***********************************************************************
subroutine output_penciled_scal_c(filename,pencil,&
                                  ndim,i,iy,iz,t, &
                                  nx,ny,nz,nghost,fnlen)
  use Cparam, only: NO_WARN
  use Cdata, only: ip,mx,headt,imn

  real,dimension(mx) :: pencil
  real :: t
  integer :: ndim,i,iy,iz,nx,ny,nz,nghost,fnlen
  character (len=*) :: filename

  if (headt .and. (imn==1)) print*, &
       'OUTPUT_PENCIL: Not writing to ', trim(filename), &
       ' since DEBUG=nodebug'

  ! to keep compiler quiet...
  if (NO_WARN) print*,pencil(1),ndim,i,iy,iz,t,nx,ny,nz,nghost,fnlen
endsubroutine output_penciled_scal_c


!!! End of file nodebug.f90
