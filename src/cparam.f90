! $Id: cparam.f90,v 1.39 2004-05-27 20:07:30 mee Exp $

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
  integer, parameter :: ikind8=selected_int_kind(14) ! 8-byte integer kind
  integer, parameter :: nghost=3
  integer, parameter :: nx=nxgrid,ny=nygrid/nprocy,nz=nzgrid/nprocz
  integer(KIND=ikind8), parameter :: nw=nx*ny*nz
  integer, parameter :: mx=nx+2*nghost,l1=1+nghost,l2=mx-nghost
  integer, parameter :: my=ny+2*nghost,m1=1+nghost,m2=my-nghost
  integer, parameter :: mz=nz+2*nghost,n1=1+nghost,n2=mz-nghost
  integer, parameter :: mw=mx*my*mz,nwgrid=nxgrid*nygrid*nzgrid
!
  integer, parameter :: l1i=l1+nghost-1,l2i=l2-nghost+1
  integer, parameter :: m1i=m1+nghost-1,m2i=m2-nghost+1
  integer, parameter :: n1i=n1+nghost-1,n2i=n2-nghost+1
!
  integer, parameter :: nrcyl=nx/2 ! used for azimuthal averages
!
!  array dimension for reduce operation (maxima and sums)
!  use here symbol mreduce, use nreduce in call
!
  integer, parameter :: mreduce=6
  integer :: ip=14
!
!  length of file names
!            strings for boundary condition,
!            labels a la initss, initaa,
!            lines to be read in
!            date-and-time string
!
  integer, parameter :: fnlen=128,bclen=3,labellen=25,linelen=256,datelen=30
!
!  number of slots in initlnrho etc.
!
  integer, parameter :: ninit=4
!
!  significant length of random number generator state
!  Different compilers have different lengths:
!    NAG: 1, Compaq: 2, Intel: 47, SGI: 64, NEC: 256
  integer, parameter :: mseed=256
!
!  a marker value that is highly unlikely (``impossible'') to ever occur
!  during a meaningful run.
!  Maybe using NaN (how do you set this in F90?) would be better..
!
  real, parameter :: impossible=3.9085e37
!
!
! Diagnostic variable types 
!
!     values greater than 0 get maxed across all  processors before any transformation using mpi_reduce_max
!     values less than 0 get summed over all processors before any transformation using mpi_reduce_sum
!     the value 0 causes the value simply to be used from the root processor
!
  integer, parameter :: ilabel_max=-1,ilabel_sum=1,ilabel_save=0,ilabel_max_sqrt=-2,ilabel_sum_sqrt=2
  integer, parameter :: ilabel_max_dt=-3,ilabel_max_neg=-4, ilabel_max_reciprocal=-5
  integer, parameter :: ilabel_integrate=3,ilabel_surf=4
!
endmodule Cparam

