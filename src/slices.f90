! $Id: slices.f90,v 1.3 2002-06-01 09:36:38 brandenb Exp $

module Slices

  use Cdata

!  generate slices for animation purposes

  real, dimension (nx,ny,3) :: uu_xy
  real, dimension (nx,ny) :: lnrho_xy,divu_xy

  real, dimension (nx,nz,3) :: uu_xz
  real, dimension (nx,nz) :: lnrho_xz,divu_xz

endmodule Slices
