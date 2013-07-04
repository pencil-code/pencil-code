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
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'initial_condition.h'
!
  real, dimension(3) :: amp0 = 0.0, k0 = 0.0
  real :: phase0 = 0.0
!
  namelist /initial_condition_pars/ amp0, k0, phase0
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Registers variables associated with this module; likely none.
!
!  25-jun-13/ccyang: coded.
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initializes any module variables which are parameter dependent.
!
!  25-jun-13/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initializes the velocity field.
!
!  25-jun-13/ccyang: coded.
!
      use Initcond, only: sinwave_phase
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: i
!
      do i = 1, 3
        call sinwave_phase(f, iuu+i-1, amp0(i), k0(1), k0(2), k0(3), phase0)
      enddo
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initializes the magnetic field.
!
!  28-jun-13/ccyang: coded.
!
      use EquationOfState, only: rho0
      use Initcond, only: sinwave_phase
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: i
      real :: c
!
      c = sqrt(mu0 * rho0)
      do i = 1, 3
        call sinwave_phase(f, ibb+i-1, c * amp0(i), k0(1), k0(2), k0(3), phase0)
      enddo
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(unit, iostat)
!
! Reads the initialization parameters for the initial conditions.
!
!  25-jun-13/ccyang: coded.
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      integer :: istat
!
      read(unit, NML=initial_condition_pars, IOSTAT=istat)
!
      if (present(iostat)) then
        iostat = istat
      else if (istat /= 0) then
        call fatal_error('read_initial_condition_pars', 'cannot read the namelist initial_condition_pars. ')
      endif
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
! Writes the initialization parameters for the initial conditions.
!
!  25-jun-13/ccyang: coded.
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
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
