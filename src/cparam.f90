module Cparam

!!!  Parameters

!  (nx,ny,nz) is the size of the computational mesh
!  The total nmumber of meshpoints is (nx,ny,nz*ncpus).
!  The number of ghost zones is NOT counted.
!
  integer, parameter :: ncpus=1,nprocz=1,nprocy=ncpus/nprocz,nprocx=1
  integer, parameter :: nxgrid=32,nygrid=nxgrid,nzgrid=nxgrid
  integer, parameter :: nghost=3,mk=3000,mvar=8
  integer, parameter :: bclen=3
!
!  derived parameters
!
  integer, parameter :: nx=nxgrid,ny=nygrid/nprocy,nz=nzgrid/nprocz,nw=nx*ny*nz
  integer, parameter :: mx=nx+2*nghost,l1=1+nghost,l2=mx-nghost
  integer, parameter :: my=ny+2*nghost,m1=1+nghost,m2=my-nghost
  integer, parameter :: mz=nz+2*nghost,n1=1+nghost,n2=mz-nghost
  integer, parameter :: mw=mx*my*mz
!
  integer, parameter :: l1i=l1+nghost-1,l2i=l2-nghost+1
  integer, parameter :: m1i=m1+nghost-1,m2i=m2-nghost+1
  integer, parameter :: n1i=n1+nghost-1,n2i=n2-nghost+1
!
!  array dimension for reduce operation (maxima and sums)
!  use here symbol mreduce, use nreduce in call
!  u, divu (gives a total of 2)
!
  integer, parameter :: mreduce=6
  integer :: ip=14
!

endmodule Cparam
