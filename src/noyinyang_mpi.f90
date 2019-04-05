! $Id$
!
! MODULE_DOC: This module contains Yin-Yang related dummy types and functions.
!
!**************************************************************************
!
module Yinyang_mpi
!
  use Cparam

  include 'yinyang_mpi.h'

contains

!**************************************************************************
    subroutine reduce_zsum(arrm,temp)
!
!  Calculate z-sums (still depending on x and y).
!
!  25-apr-16/MR: outsourced from zaverages_xy
!
      use Mpicomm, only: mpireduce_sum

      real, dimension(:,:,:), intent(IN) :: arrm
      real, dimension(:,:,:), intent(OUT):: temp
!
      integer, dimension(3) :: sz
!
!  Communicate over all processors along z beams.
!  The result is only present on the z-root processors.
!
        sz=(/size(arrm,1),size(arrm,2),size(arrm,3)/)
!
!  Summation of local arrays arrm into temp of z root processors.
!
        call mpireduce_sum(arrm,temp,sz,idir=3)
!
    endsubroutine reduce_zsum
!***********************************************************************
    subroutine zsum_yy(fnamexy,iname,m,n,a)

      use General, only: keep_compiler_quiet

      real, dimension(:,:,:), intent(INOUT) :: fnamexy
      integer,                intent(IN)    :: iname,m,n
      real, dimension(:),     intent(IN)    :: a

      call keep_compiler_quiet(fnamexy)
      call keep_compiler_quiet(iname,m,n)
      call keep_compiler_quiet(a)

    endsubroutine zsum_yy
!***********************************************************************
    subroutine initialize_zaver_yy(nlines,nycap_)

      use General, only: keep_compiler_quiet

      integer :: nlines,nycap_

      call keep_compiler_quiet(nlines)
      call keep_compiler_quiet(nycap_)

    endsubroutine initialize_zaver_yy
!**************************************************************************
endmodule Yinyang_mpi
