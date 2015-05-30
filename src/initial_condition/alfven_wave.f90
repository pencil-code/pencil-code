! $Id$
!
!  This module sets up an Alfvéén wave of amplitude init_amp0
!  with wave vector init_k0 and an optional phase init_phase0.
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
  include '../initial_condition.h'
!
  real, dimension(3) :: init_amp0 = 0.0, init_k0 = 0.0
  real :: init_phase0 = 0.0
!
  namelist /initial_condition_pars/ init_amp0, init_k0, init_phase0
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
!  26-sep-14/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
!  Oscillations must be perpendicular to the wave vector.
!
      if (sum(init_amp0 * init_k0) /= 0.0) &
        call fatal_error('initialize_initial_condition', 'init_amp0 and init_k0 are not perpendicular. ')
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
        call sinwave_phase(f, iuu+i-1, init_amp0(i), init_k0(1), init_k0(2), init_k0(3), init_phase0)
      enddo
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initializes the magnetic field.
!
!  26-sep-14/ccyang: coded.
!
      use EquationOfState, only: rho0
      use Initcond, only: sinwave_phase, coswave_phase
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: i
      real :: c
!
      c = sqrt(mu0 * rho0)
      bfield: if (lbfield) then
        do i = 1, 3
          call sinwave_phase(f, ibb+i-1, c * init_amp0(i), init_k0(1), init_k0(2), init_k0(3), init_phase0)
        enddo
      else bfield
        c = c / init_k0(3)
        call coswave_phase(f, iax, -c * init_amp0(2), init_k0(1), init_k0(2), init_k0(3), init_phase0)
        call coswave_phase(f, iay, c * init_amp0(1), init_k0(1), init_k0(2), init_k0(3), init_phase0)
        f(l1:l2,m1:m2,n1:n2,iaz) = 0.0
      endif bfield
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include "../parallel_unit.h"
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
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
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
