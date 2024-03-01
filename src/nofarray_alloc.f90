  module farray_alloc

  use Cparam

  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (mx,my,mz,mvar) :: df

  public :: f,df
  public :: initialize, finalize

  contains
!***********************************************
  subroutine initialize

  endsubroutine initialize
!***********************************************
  subroutine finalize

  endsubroutine finalize
!***********************************************
  endmodule
