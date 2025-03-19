  module farray_alloc

  use Cparam

  implicit none

  real, dimension(mx,my,mz,mfarray) :: f
  real, dimension(:,:,:,:), allocatable :: df

  public :: f,df
  public :: initialize, finalize

  contains
!***********************************************
  subroutine initialize

    use Cdata, only: nt
    use Messages, only: fatal_error

    integer :: stat

    if (nt>0.and..not.lgpu) then
      allocate(df(mx,my,mz,mvar),STAT=stat)
      if (stat>0) call fatal_error('farray_alloc','Could not allocate df')
    else
      allocate(df(1,1,1,1))
    endif

  endsubroutine initialize
!***********************************************
  subroutine finalize

    if (allocated(df)) deallocate(df)

  endsubroutine finalize
!***********************************************
  endmodule
