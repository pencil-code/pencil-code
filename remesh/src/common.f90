module Common

!!! Global variables

  use Cdata
  !
  include 'common.local'
  !
  !  New grid parameters
  !
  integer, parameter :: nnx=remesh_parx*nx/mulx
  integer, parameter :: nny=remesh_pary*ny/muly
  integer, parameter :: nnz=remesh_parz*nz/mulz
  integer, parameter :: mmx=remesh_parx*nx/mulx+2*nghost
  integer, parameter :: mmy=remesh_pary*ny/muly+2*nghost
  integer, parameter :: mmz=remesh_parz*nz/mulz+2*nghost
  integer, parameter :: nnx_grid=remesh_parx*nx
  integer, parameter :: nny_grid=remesh_pary*ny
  integer, parameter :: nnz_grid=remesh_parz*nz
  integer, parameter :: mmx_grid=remesh_parx*nx+2*nghost
  integer, parameter :: mmy_grid=remesh_pary*ny+2*nghost
  integer, parameter :: mmz_grid=remesh_parz*nz+2*nghost
  integer :: ll1=l1
  integer :: ll2=remesh_parx*nx+nghost
  integer :: mm1=m1
  integer :: mm2=remesh_pary*ny+nghost
  integer :: nn1=n1
  integer :: nn2=remesh_parz*nz+nghost
  integer :: lll1=l1
  integer :: lll2=remesh_parx*nx/mulx+nghost
  integer :: mmm1=m1
  integer :: mmm2=remesh_pary*ny/muly+nghost
  integer :: nnn1=n1
  integer :: nnn2=remesh_parz*nz/mulz+nghost

end module Common
