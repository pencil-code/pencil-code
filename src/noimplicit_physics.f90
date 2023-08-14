! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lADI = .false.
!
!***************************************************************
module ImplicitPhysics
!
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'implicit_physics.h'
!
  contains
!***********************************************************************
    subroutine register_implicit_physics
!
    endsubroutine register_implicit_physics
!***********************************************************************
    subroutine initialize_implicit_physics
!
    endsubroutine initialize_implicit_physics
!***********************************************************************
    subroutine init_implicit_physics(f)
!
      real, dimension(mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine init_implicit_physics
!***********************************************************************
    subroutine calc_heatcond_ADI(f)
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_heatcond_ADI
!***********************************************************************
endmodule ImplicitPhysics
