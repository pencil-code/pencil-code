  module farray_alloc

  use Cparam

  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (mx,my,mz,mvar) :: df
  type (pencil_case) :: p

  public :: f,df,p
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
