! $Id: slices.f90,v 1.5 2002-08-17 08:40:22 brandenb Exp $

module Slices

  use Cdata

!  generate slices for animation purposes

  real, dimension (nx,ny,3) :: uu_xy,bb_xy
  real, dimension (nx,ny) :: lnrho_xy,divu_xy

  real, dimension (nx,nz,3) :: uu_xz,bb_xz
  real, dimension (nx,nz) :: lnrho_xz,divu_xz,ss_xz

endmodule Slices
