module Debugging

!
! Debugging utilities.
!

  implicit none

  interface output_stenc        ! Overload the `output_stenc' function
    module procedure output_stenc_vect
    module procedure output_stenc_scal
  endinterface

  contains

!***********************************************************************
    subroutine output_stenc_vect(file,a,nn,imn)
!
!  Write snapshot file of stenciled vector data.
!  Wrapper to the C routine output_stenciled_c.
!
!  15-feb-02/wolf: coded
!
      use Cdata
!
      integer :: nn,imn
      real, dimension (mx,nn) :: a
      character (LEN=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'OUTPUT_STENC_VECT: nn =', nn
!
      call output_stenciled_c(file,a,nn,imn,t,nx,ny,nz,len(file))
!
    endsubroutine output_stenc_vect
!***********************************************************************
    subroutine output_stenc_scal(file,a,nn,imn)
!
!  Write snapshot file of stenciled scalar data.
!  Wrapper to the C routine output_stenciled_c.
!
!  15-feb-02/wolf: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
!
      integer :: nn,imn
      real, dimension (mx) :: a
      character (LEN=*) :: file
!
      if ((ip<=8) .and. lroot) print*,'OUTPUT_SCALAR'
      if (nn /= 1) call stop_it("OUTPUT called with scalar field, but nn/=1")
!
      call output_stenciled_c(file,a,nn,imn,t,nx,ny,nz,len(file))
!
    endsubroutine output_stenc_scal
!***********************************************************************

endmodule Debugging

!!! End of file debugging.f90
