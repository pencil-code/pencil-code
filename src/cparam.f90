! $Id: cparam.f90,v 1.27 2002-10-04 18:02:27 dobler Exp $

module Cparam

!!!  Parameters

!  (nx,ny,nz) is the size of the computational mesh
!  The total nmumber of meshpoints is (nx*nprocx,ny*nprocy,nz*nprocz).
!  The number of ghost zones is NOT counted.
!
!  In practice, the user will change the number of cpus (in y and z)
!  and the number of mesh points, and recompile.
!  Dependening on what is invoked under Makefile.local,
!  one needs to adjust nvar.
!  This part is now isolated in a separate cparam.local file.
!
  include 'cparam.local'
  include 'cparam.inc'
!
!  derived and fixed parameters
!
  integer, parameter :: nghost=3
  integer, parameter :: nx=nxgrid,ny=nygrid/nprocy,nz=nzgrid/nprocz,nw=nx*ny*nz
  integer, parameter :: mx=nx+2*nghost,l1=1+nghost,l2=mx-nghost
  integer, parameter :: my=ny+2*nghost,m1=1+nghost,m2=my-nghost
  integer, parameter :: mz=nz+2*nghost,n1=1+nghost,n2=mz-nghost
  integer, parameter :: mw=mx*my*mz,nwgrid=nxgrid*nygrid*nzgrid
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
!  length of strings for boundary condition,
!            labels a la initss, initaa,
!            lines to be read in
!            date-and-time string
!
  integer, parameter :: bclen=3,labellen=25,linelen=256,datelen=30
!
!  significant length of random number generator state
!  Different compilers have different lengths:
!    SGI: 64, Intel: 47, NAG: 1, Compaq: 2
  integer, parameter :: mseed=100
!
!  a marker value that is highly unlikely (``impossible'') to ever occur
!  during a meaningful run.
!  Maybe using NaN (how do you set this in F90?) would be better..
!
  real, parameter :: impossible=3.9085e37
!
endmodule Cparam

