! $Id$
!
!  Dummy module.
!
module Slices_methods
!
  use Cdata
  use General, only: keep_compiler_quiet

  implicit none
!
  include 'slices_methods.h'
!
  interface assign_slices_scal
    module procedure assign_slices_f_scal
    module procedure assign_slices_sep_scal
  endinterface

  interface assign_slices_vec
    module procedure assign_slices_f_vec
    module procedure assign_slices_sep_vec
  endinterface

  interface alloc_slice_buffers
    module procedure alloc_slice_buffers_scal
    module procedure alloc_slice_buffers_vec
  endinterface
!
contains
!***********************************************************************
    subroutine assign_slices_sep_scal(slices,xy,xz,yz,xy2,xy3,xy4,xz2)

      type (slice_data) :: slices
      real, dimension(:,:), target :: xy,xz,yz,xy2,xy3,xy4,xz2

      call keep_compiler_quiet(xy,xz,yz,xy2)
      call keep_compiler_quiet(xy3,xy4,xz2)
      slices%ready=.false.

    endsubroutine assign_slices_sep_scal
!***********************************************************************
    subroutine assign_slices_sep_vec(slices,xy,xz,yz,xy2,xy3,xy4,xz2)

      type (slice_data) :: slices
      real, dimension(:,:,:), target :: xy,xz,yz,xy2,xy3,xy4,xz2

      call keep_compiler_quiet(xy,xz,yz,xy2)
      call keep_compiler_quiet(xy3,xy4,xz2)
      slices%ready=.false.

    endsubroutine assign_slices_sep_vec
!***********************************************************************
    subroutine assign_slices_f_scal(slices,f,ind)

      type (slice_data) :: slices
      real, dimension(:,:,:,:) :: f
      integer :: ind

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ind)
      slices%ready=.false.

    endsubroutine assign_slices_f_scal
!***********************************************************************
    subroutine assign_slices_f_vec(slices,f,ind)

      type (slice_data) :: slices
      real, dimension(:,:,:,:) :: f
      integer :: ind

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ind)
      slices%ready=.false.

    endsubroutine assign_slices_f_vec
!***********************************************************************
    subroutine store_slices(pencil,xy,xz,yz,xy2,xy3,xy4,xz2)

      real, dimension(nx) :: pencil
      real, dimension(:,:) :: xy,xz,yz,xy2,xy3,xy4,xz2

      call keep_compiler_quiet(xy,xz,yz,xy2)
      call keep_compiler_quiet(pencil,xy3,xy4,xz2)

    endsubroutine store_slices
!***********************************************************************
    subroutine alloc_slice_buffers_scal(xy,xz,yz,xy2,xy3,xy4,xz2)

      real, dimension(:,:), allocatable :: xy,xy2,xy3,xy4
      real, dimension(:,:), allocatable :: xz,xz2
      real, dimension(:,:), allocatable :: yz

      call keep_compiler_quiet(xy,xz,yz,xy2)
      call keep_compiler_quiet(xy3,xy4,xz2)

    endsubroutine alloc_slice_buffers_scal
!***********************************************************************
    subroutine alloc_slice_buffers_vec(xy,xz,yz,xy2,xy3,xy4,xz2)

      real, dimension(:,:,:), allocatable :: xy,xy2,xy3,xy4
      real, dimension(:,:,:), allocatable :: xz,xz2
      real, dimension(:,:,:), allocatable :: yz

      call keep_compiler_quiet(xy,xz,yz,xy2)
      call keep_compiler_quiet(xy3,xy4,xz2)

    endsubroutine alloc_slice_buffers_vec
!***********************************************************************
endmodule Slices_methods
