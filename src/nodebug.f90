! $Id$
!
!  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
!  Date:   24-Jun-2002
!
!  Description:
!   Two dummy debugging routines (in case the C stuff doesn't work).
!
!***********************************************************************
subroutine output_penciled_vect_c(filename,pencil,ndim,i,iy,iz,t, &
                                  nx,ny,nz,nghost,fnlen)
!
  use Cdata, only: mx,headt,imn
  use General, only: keep_compiler_quiet
!
  implicit none
!
  real,dimension(mx,*) :: pencil
  real :: t
  integer :: ndim,i,iy,iz,nx,ny,nz,nghost,fnlen
  character (len=*) :: filename
!
  if (headt .and. (imn==1)) print*, &
       'OUTPUT_PENCIL: Not writing to ', trim(filename), &
       ' since DEBUG=nodebug'
!
  call keep_compiler_quiet(pencil(1,1))
  call keep_compiler_quiet(ndim)
  call keep_compiler_quiet(i)
  call keep_compiler_quiet(iy)
  call keep_compiler_quiet(iz)
  call keep_compiler_quiet(t)
  call keep_compiler_quiet(nx)
  call keep_compiler_quiet(ny)
  call keep_compiler_quiet(nz)
  call keep_compiler_quiet(nghost)
  call keep_compiler_quiet(fnlen)
!
endsubroutine output_penciled_vect_c
!***********************************************************************
subroutine output_penciled_scal_c(filename,pencil,ndim,i,iy,iz,t, &
                                  nx,ny,nz,nghost,fnlen)
!
  use Cdata, only: mx,headt,imn
  use General, only: keep_compiler_quiet
!
  real,dimension(mx) :: pencil
  real :: t
  integer :: ndim,i,iy,iz,nx,ny,nz,nghost,fnlen
  character (len=*) :: filename
!
  if (headt .and. (imn==1)) print*, &
       'OUTPUT_PENCIL: Not writing to ', trim(filename), &
       ' since DEBUG=nodebug'
!
  call keep_compiler_quiet(pencil(1))
  call keep_compiler_quiet(ndim)
  call keep_compiler_quiet(i)
  call keep_compiler_quiet(iy)
  call keep_compiler_quiet(iz)
  call keep_compiler_quiet(t)
  call keep_compiler_quiet(nx)
  call keep_compiler_quiet(ny)
  call keep_compiler_quiet(nz)
  call keep_compiler_quiet(nghost)
  call keep_compiler_quiet(fnlen)
!
endsubroutine output_penciled_scal_c
