module Slices

  use Cdata

  real, dimension (nx,ny,3) :: uu_slice  !(slices for animation purposes)
  real, dimension (nx,ny) :: divu_slice

endmodule Slices
