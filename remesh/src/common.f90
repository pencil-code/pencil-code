module Common

!!! Global variables

  use Cdata
  !
  include 'common.local'
  !
  !  New grid parameters
  !
  integer, parameter :: nxcoll=divx*nx, nycoll=divy*ny, nzcoll=divz*nz
  integer, parameter :: mxcoll=nxcoll+2*nghost, mycoll=nycoll+2*nghost, mzcoll=nzcoll+2*nghost

  integer, parameter :: mprocs=mulz*muly*mulx
  
  integer, parameter :: nnx=remesh_parx*nxcoll/mulx
  integer, parameter :: nny=remesh_pary*nycoll/muly
  integer, parameter :: nnz=remesh_parz*nzcoll/mulz
  integer, parameter :: mmx=nnx+2*nghost
  integer, parameter :: mmy=nny+2*nghost
  integer, parameter :: mmz=nnz+2*nghost
  integer, parameter :: nnx_grid=remesh_parx*nxcoll
  integer, parameter :: nny_grid=remesh_pary*nycoll
  integer, parameter :: nnz_grid=remesh_parz*nzcoll
  integer, parameter :: mmx_grid=nnx_grid+2*nghost
  integer, parameter :: mmy_grid=nny_grid+2*nghost
  integer, parameter :: mmz_grid=nnz_grid+2*nghost
  !
  integer, parameter :: ll1=l1
  integer, parameter :: ll2=nnx_grid+nghost
  integer, parameter :: mm1=m1
  integer, parameter :: mm2=nny_grid+nghost
  integer, parameter :: nn1=n1
  integer, parameter :: nn2=nnz_grid+nghost
  !
  integer, parameter :: lll1=l1
  integer, parameter :: lll2=nnx+nghost
  integer, parameter :: mmm1=m1
  integer, parameter :: mmm2=nny+nghost
  integer, parameter :: nnn1=n1
  integer, parameter :: nnn2=nnz+nghost

end module Common
