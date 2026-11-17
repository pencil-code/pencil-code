! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
module InitialCondition
!
  use Cdata
  use Quiet
!
  implicit none
!
  include '../initial_condition.h'
!
  contains
!***********************************************************************
    subroutine initial_condition_lncc(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real :: alpha=10
!
      do n=n1,n2
        do m=m1,m2
          f(l1:l2,m,n,icc) = (1 + tanh((1 - 2*abs(x(l1:l2)))*alpha))*(1 + tanh((1-1.4*y(m))*alpha/4))/2
        enddo
      enddo
!
    endsubroutine initial_condition_lncc
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
