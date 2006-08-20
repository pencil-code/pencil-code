! $Id: noborder_profiles.f90,v 1.3 2006-08-20 22:19:56 wlyra Exp $ 
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
    subroutine border_driving(f,df,f_target,j)
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: f_target
      integer :: j
! 
!  Dummy routine
!
      if (NO_WARN) print*,j,f,df,f_target
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
