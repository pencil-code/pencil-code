  module Farray_alloc
!
! Dynamical allocation of f and df.
! 
  implicit none

  real, dimension(:,:,:,:), allocatable :: df
  real, dimension(:,:,:,:), pointer :: f
  real, dimension(:,:,:,:), allocatable, target :: f_arr

  public :: f, df
  public :: initialize, finalize

  external :: allocate_shm

  contains
!******************************************************************************
  subroutine initialize

  use Cdata
  use Messages, only: fatal_error
  use Syscalls, only: memusage
  use iso_c_binding

  integer :: stat
  integer(KIND=ikind8) :: nelems=mx*my*mz*mfarray
  type(C_PTR) :: fp

  interface
    type(C_PTR) function allocate_shm(num,name)
      import :: c_ptr, ikind8
      integer(KIND=ikind8) :: num
      character(LEN=*) :: name
    endfunction
  end interface

    !allocate( f(mx,my,mz,nvar+naux+nscratch+nglobal),STAT=stat)
    !if (stat>0) call fatal_error('farray_alloc','Could not allocate memory for f')
    !allocate(df(mx,my,mz,nvar),STAT=stat)
    !if (stat>0) call fatal_error('farray_alloc','Could not allocate memory for df')
    !!print*, 'stat,mfarry,nfarray=',stat,mvar,maux,mscratch,mglobal,nvar,naux,nscratch,nglobal

    !mvar=nvar; maux=naux; maux_com=naux_com; mscratch=nscratch; mglobal=nglobal

    if (shared_mem_name/='') then
      fp = allocate_shm(nelems,shared_mem_name//char(0))
      call c_f_pointer(fp,f,(/mx,my,mz,mfarray/))
    else
      allocate(f_arr(mx,my,mz,mfarray),STAT=stat)
      f => f_arr
    endif

    if (stat>0) call fatal_error('farray_alloc','Could not allocate f')
    if (nt>0.and..not.lgpu) then
      allocate(df(mx,my,mz,mvar),STAT=stat)
      if (stat>0) call fatal_error('farray_alloc','Could not allocate df')
    else
      allocate(df(1,1,1,1))
    endif

  endsubroutine initialize
!******************************************************************************
  subroutine finalize

    if (allocated(f_arr)) deallocate(f_arr)
    !if (shared_mem_name=='') deallocate fp   !?
    if (allocated(df)) deallocate(df) 

  endsubroutine finalize
!******************************************************************************
  endmodule Farray_alloc
