! $Id$
!
!  This module calculates a divergence of u based shock finding
!  profile used by shock viscosities and diffusion terms.
!    eg. the total voscosity is taken as:
!           nu_total = nu + nu_shock*dx*smooth(max5(-(div u))))
!    where dx*smooth(max5(-(div u)))) is the profile calculated
!    here in.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lshock = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED shock; gshock(3); shock_perp; gshock_perp(3)
!
!***************************************************************
module Shock
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'shock.h'
!
  contains
!***********************************************************************
    subroutine register_shock()
!
!  19-nov-02/tony: coded
!  24-jan-05/tony: modified from visc_shock.f90
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_shock
!***********************************************************************
    subroutine initialize_shock(f,lstarting)
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_shock
!***********************************************************************
    subroutine read_shock_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_shock_init_pars
!***********************************************************************
    subroutine write_shock_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_shock_init_pars
!***********************************************************************
    subroutine read_shock_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_shock_run_pars
!***********************************************************************
    subroutine write_shock_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_shock_run_pars
!***********************************************************************
    subroutine rprint_shock(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_shock
!***********************************************************************
    subroutine get_slices_shock(f,slices)
!
!  Write slices for animation of Shock variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_shock
!***********************************************************************
    subroutine pencil_criteria_shock()
!
!  All pencils that the Shock module depends on are specified here.
!
!  20-11-04/anders: coded
!  24-jan-05/tony: modified from visc_shock.f90
!
!
    endsubroutine pencil_criteria_shock
!***********************************************************************
    subroutine pencil_interdep_shock(lpencil_in)
!
!  Interdependency among pencils from the Viscosity module is specified here.
!
!  20-11-04/anders: coded
!  24-jan-05/tony: modified from visc_shock.f90
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_shock
!***********************************************************************
    subroutine calc_pencils_shock(f,p)
!
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!  24-jan-05/tony: modified from visc_shock.f90
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! shock
      if (lpencil(i_shock)) p%shock=0.
! gshock
      if (lpencil(i_gshock)) p%gshock=0.
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_shock
!!***********************************************************************
    subroutine calc_shock_profile_simple(f)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  23-nov-02/tony: coded
!  24-jan-05/tony: modified from visc_shock.f90
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_shock_profile_simple
!!***********************************************************************
    subroutine calc_shock_profile(f)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  23-nov-02/tony: coded
!  24-jan-05/tony: modified from visc_shock.f90
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_shock_profile
!***********************************************************************
endmodule Shock
