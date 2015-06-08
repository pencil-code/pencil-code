! $Id$
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
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
contains
!***********************************************************************
  subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  04-sep-10/bing: coded
!
    if (lroot) call svn_id( &
        "$Id$")
!
  endsubroutine register_initial_condition
!***********************************************************************
  subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  14-dec-10/bing: coded
!
    real, dimension (mx,my,mz,mfarray) :: f
!
    call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
  subroutine initial_condition_all(f,profiles)
!
!  Initialize logarithmic density.
!
!  04-sep-10/bing: coded
!  10-feb-15/MR: added optional parameter 'profiles' (intended to replace f)
!
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
!
print*,f(l1:l2,m1:m2,n1:n2,:)
print*
!print*,df(l1:l2,m1:m2,n1:n2,:)
!print*
!    f=f+.0001*df
print*,f(l1:l2,m1:m2,n1:n2,:)
print*
!
     call keep_compiler_quiet(profiles)
!
  endsubroutine initial_condition_all
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
