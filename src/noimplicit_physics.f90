! $Id$
!***************************************************************
module ImplicitPhysics
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'implicit_physics.h'
!
!
  interface heatcond_TT ! Overload subroutine `hcond_TT' function
    module procedure heatcond_TT_0d  ! get one value (hcond, dhcond)
    module procedure heatcond_TT_1d  ! get 1d-arrays (hcond, dhcond)
    module procedure heatcond_TT_2d  ! get 2d-arrays (hcond, dhcond)
  end interface
!
  real, pointer :: hcond0, Fbot, Tbump, Kmax, hole_slope, hole_width, &
                   hole_alpha
!
  contains
!***********************************************************************
    subroutine init_param_ADI()
!
    endsubroutine init_param_ADI
!***********************************************************************
    subroutine calc_heatcond_ADI(finit,f)
!
!  10-sep-07/gastine+dintrans: wrapper to the two possible ADI subroutines
!  ADI_Kconst: constant radiative conductivity
!  ADI_Kprof: radiative conductivity depends on T, i.e. hcond(T)
!
      real, dimension(mx,my,mz,mfarray) :: finit, f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(finit)
!
    endsubroutine calc_heatcond_ADI
!***********************************************************************
    subroutine heatcond_TT_2d(TT, hcond, dhcond)
!
      implicit none
!
      real, dimension(mx,mz) :: TT, arg, hcond
      real, dimension(mx,mz), optional :: dhcond
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
    endsubroutine heatcond_TT_2d
!***********************************************************************
    subroutine heatcond_TT_1d(TT, hcond, dhcond)
!
      implicit none
!
      real, dimension(:)           :: TT, hcond
      real, dimension(:), optional :: dhcond
      real, dimension(size(TT,1))  :: arg
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
    endsubroutine heatcond_TT_1d
!***********************************************************************
    subroutine heatcond_TT_0d(TT, hcond, dhcond)
!
!
      implicit none
!
      real :: TT, arg, hcond
      real, optional :: dhcond
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
    endsubroutine heatcond_TT_0d
!***********************************************************************
endmodule ImplicitPhysics
