! $Id: slices.f90,v 1.2 2002-06-01 02:56:21 brandenb Exp $

module Slices

  use Cdata

  real, dimension (nx,ny,3) :: uu_slice  !(slices for animation purposes)
  real, dimension (nx,ny) :: divu_slice

endmodule Slices
