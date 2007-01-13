! $Id: noborder_profiles.f90,v 1.7 2007-01-13 21:46:44 dobler Exp $
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lborder_profiles = .false.
!
!***************************************************************

module BorderProfiles

  use Cparam
  use Cdata

  implicit none

  private

  include 'border_profiles.h'
!

  contains

!***********************************************************************
    subroutine initialize_border_profiles()
!
!  Position-dependent quenching factor that multiplies rhs of pde
!  by a factor that goes gradually to zero near the boundaries.
!  border_frac is a 3-D array, separately for all three directions.
!  border_frac=1 would affect everything between center and border.
!
    ! DUMMY ROUTINE
!
    endsubroutine initialize_border_profiles
!***********************************************************************
    subroutine pencil_criteria_borderprofiles()
!
!  All pencils that this module depends on are specified here.
!
!  DUMMY ROUTINE
!
    endsubroutine pencil_criteria_borderprofiles
!***********************************************************************
    subroutine border_driving(f,df,p,f_target,j)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: f_target
      integer :: j
!
!  Dummy routine
!
      if (NO_WARN) print*,j,f,p,df,f_target
!
    endsubroutine border_driving
!***********************************************************************
    subroutine border_quenching(df,j)
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: j
!
!  Dummy routine
!
      if (NO_WARN) print*,j,df
!
    endsubroutine border_quenching
!***********************************************************************
endmodule BorderProfiles
