! $Id$
!
!  This module provide a way for users to specify custom initial
!  conditions.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'initial_condition.h'
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      use Messages, only: svn_id
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize each component of the magnetic vector potential with
!  a delta function at the origin.
!
!  22-may-12/ccyang: coded.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: i, j, k
!
      f(:,:,:,iax:iaz) = 0.
!
      forall (i = l1:l2, j = m1:m2, k = n1:n2, &
              abs(x(i)) * dx_1(i) < 1. .and. abs(y(j)) * dy_1(j) < 1. .and. abs(z(k)) * dz_1(k) < 1.) &
        f(i,j,k,iax:iaz) = dx_1(i) * dy_1(j) * dz_1(k)
!
    endsubroutine initial_condition_aa
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include 'initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
