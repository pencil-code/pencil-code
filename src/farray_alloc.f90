  module farray_alloc
!
! Dynamical allocation of f and df.
!
  real, dimension(:,:,:,:), allocatable :: f,df

  public :: f,df
  public :: alloc

  contains
!******************************************************************************
  subroutine alloc 

  use Cdata
  use Messages, only: fatal_error

  integer :: stat

    allocate( f(mx,my,mz,nvar+naux+nscratch+nglobal),STAT=stat)
    if (stat>0) call fatal_error('farray_alloc','Could not allocate memory for f')
    allocate(df(mx,my,mz,nvar),STAT=stat)
    if (stat>0) call fatal_error('farray_alloc','Could not allocate memory for df')
    print*, 'stat,mfarry,nfarray=',stat,mvar,maux,mscratch,mglobal,nvar,naux,nscratch,nglobal

  endsubroutine alloc 
!******************************************************************************
  endmodule
