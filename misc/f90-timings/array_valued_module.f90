!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   array_valued.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  Author: wd (wdobler [at] gmail [dot] com)
!!!  Date:   18-Apr-2010
!!!
!!!  Description:
!!!   Provide array_valued functions for common operations on pencils.

module Array_valued

  use Timings

  implicit none

  interface operator (*)
    module procedure pencil_multiply1
  endinterface

  interface operator (/)
    module procedure pencil_divide1
  endinterface

contains

!***********************************************************************
  subroutine multsv_mn(a,b,c)
!
!  vector multiplied with scalar, gives vector
!   22-nov-01/nils erland: coded
!   10-oct-03/axel: a is now the scalar (now consistent with old routines)
!
!    use Cdata
!
    intent(in)  :: a,b
    intent(out) :: c
!
    real, dimension (nx,3) :: b,c
    real, dimension (nx) :: a
    integer :: i
!
    do i=1,3
      c(:,i)=a*b(:,i)
    enddo
!
  endsubroutine multsv_mn
!***********************************************************************
  function pencil_multiply1(s,v)
!
!  The `*' operator may be extended through this function to allow
!  elementwise multiplication of a `pencil-scalar' with a `pencil-vector'
!
!   6-Sep-05/tobi: coded
!
!    use Cdata

    real, dimension(nx), intent(in) :: s
    real, dimension(nx,3), intent(in) :: v
    real, dimension(nx,3) :: pencil_multiply1

    integer :: i

    do i=1,3
      pencil_multiply1(:,i) = s(:) * v(:,i)
    enddo

  endfunction pencil_multiply1
!***********************************************************************
  function pencil_divide1(v,s)
!
!  The `*' operator may be extended through this function to allow
!  elementwise multiplication of a `pencil-scalar' with a `pencil-vector'
!
!   18-Sep-05/tony: coded
!
!    use Cdata

    real, dimension(nx), intent(in) :: s
    real, dimension(nx,3), intent(in) :: v
    real, dimension(nx,3) :: pencil_divide1

    integer :: i

    do i=1,3
      pencil_divide1(:,i) = v(:,i) / s(:)
    enddo

  endfunction pencil_divide1
!***********************************************************************

endmodule

!!! End of file array_valued.f90
