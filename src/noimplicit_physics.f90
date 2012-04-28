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
  use Cparam
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
  real :: Tbump, Kmax, Kmin, hole_slope, hole_width, hole_alpha
!
  contains
!***********************************************************************
    subroutine register_implicit_physics()
!
    endsubroutine register_implicit_physics
!***********************************************************************
    subroutine initialize_implicit_physics(f)
!
      use SharedVariables, only: get_shared_variable
      use MpiComm, only: stop_it
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(:), pointer :: hole_params
      integer :: ierr
!
!  Get the hole parameters if we want to run a kappa-mechanism simulation
!  in the fully-explicit case (mainly for testing purposes)
!
      if (ltemperature .and. leos_idealgas) then
        call get_shared_variable('hole_params', hole_params, ierr)
        if (ierr/=0) call stop_it("implicit_physics: "//&
            "there was a problem when getting hole_params")
        Tbump=hole_params(1)
        Kmax=hole_params(2)
        Kmin=hole_params(3)
        hole_slope=hole_params(4)
        hole_width=hole_params(5)
        hole_alpha=(Kmax-Kmin)/(pi/2.+atan(hole_slope*hole_width**2))
        if (lroot .and. ldebug) then
          print*, '************ hole parameters ************'
          print*,'Tbump, Kmax, Kmin, hole_slope, hole_width, hole_alpha=', &
              Tbump, Kmax, Kmin, hole_slope, hole_width, hole_alpha
          print*, '*****************************************'
        endif
!
        if (lrun) then
! hcondADI is dynamically shared with boundcond() for the 'c3' BC
          call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
        else
          hcondADI=spread(Kmax, 1, mx)
        endif
      endif
!
    endsubroutine initialize_implicit_physics
!***********************************************************************
    subroutine calc_heatcond_ADI(f)
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_heatcond_ADI
!***********************************************************************
    subroutine heatcond_TT_2d(TT, hcond, dhcond)
!
      real, dimension(mx,mz) :: TT, arg, hcond
      real, dimension(mx,mz), optional :: dhcond
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) &
        dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
    endsubroutine heatcond_TT_2d
!***********************************************************************
    subroutine heatcond_TT_1d(TT, hcond, dhcond)
!
      real, dimension(:)           :: TT, hcond
      real, dimension(:), optional :: dhcond
      real, dimension(size(TT,1))  :: arg
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) &
        dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
    endsubroutine heatcond_TT_1d
!***********************************************************************
    subroutine heatcond_TT_0d(TT, hcond, dhcond)
!
      real :: TT, arg, hcond
      real, optional :: dhcond
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) &
        dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
    endsubroutine heatcond_TT_0d
!***********************************************************************
endmodule ImplicitPhysics
