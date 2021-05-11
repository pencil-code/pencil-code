! $Id$
!
!  This module is the dummy for the viscosity module
!  in which e.g the viscous force or viscous heat is
!  computed.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lviscosity = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED fvisc(3); diffus_total
! PENCILS PROVIDED visc_heat; nu; gradnu(3); nu_smag
!
!***************************************************************
module Viscosity
!
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'viscosity.h'
!
  logical :: lvisc_first=.false.
  integer :: pushpars2c, pushdiags2c  ! should be procedure pointer (F2003)
!
  contains
!***********************************************************************
    subroutine register_viscosity()
!
!  19-nov-02/tony: coded
!
      use Messages, only: svn_id
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_viscosity
!***********************************************************************
    subroutine initialize_viscosity()
!
    endsubroutine initialize_viscosity
!***********************************************************************
    subroutine read_viscosity_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_viscosity_run_pars
!***********************************************************************
    subroutine write_viscosity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_viscosity_run_pars
!***********************************************************************
    subroutine rprint_viscosity(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_viscosity
!***********************************************************************
    subroutine pencil_criteria_viscosity()
!
!  All pencils that the Viscosity module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_viscosity
!***********************************************************************
    subroutine pencil_interdep_viscosity(lpencil_in)
!
!  Interdependency among pencils from the Viscosity module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_viscosity
!***********************************************************************
    subroutine calc_pencils_viscosity(f,p)
!
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
!  define visc_heat
!
      if (lpencil(i_fvisc)) p%fvisc=0.0
      if (lpencil(i_visc_heat)) p%visc_heat=0.0
      if (lpencil(i_nu)) p%nu=0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_viscosity
!***********************************************************************
    subroutine viscosity_after_boundary(f)

      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine viscosity_after_boundary
!***********************************************************************
    subroutine calc_viscosity(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_viscosity
!***********************************************************************
    subroutine calc_viscous_heat(df,p,Hmax)
!
      real, dimension (mx,my,mz,mvar)    :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax
!
      intent(in) :: df,p,Hmax
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(Hmax)
!
    endsubroutine calc_viscous_heat
!***********************************************************************
    subroutine calc_viscous_force(df,p)
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent (in) :: df,p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_viscous_force
!***********************************************************************
    subroutine calc_visc_heat_ppd(df,p)
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      intent (in) :: df,p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_visc_heat_ppd
!***********************************************************************
    subroutine getnu(nu_input,nu_pencil,ivis)
!
      real, optional, intent(out) :: nu_input
      real, dimension(nx), optional, intent(out) :: nu_pencil
      character (len=labellen), optional :: ivis
!
!  use ivis='nu-const' and put nu=1 to make
!  the particle module work (for now).
!
      if (present(nu_input))  nu_input=1.0
      if (present(nu_pencil)) nu_pencil=1.0
      if (present(ivis))      ivis='nu-const'
!
    endsubroutine getnu
!***********************************************************************
    subroutine dynamical_viscosity(uc)
!
!  Dummy routine
!
      real, intent(in) :: uc
!
      call keep_compiler_quiet(uc)
!
    endsubroutine dynamical_viscosity
!***********************************************************************
    subroutine split_update_viscosity(f)
!
!  Dummy
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine split_update_viscosity
!***********************************************************************
endmodule Viscosity
