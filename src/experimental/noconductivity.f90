! $Id$
!
!  Dummy module for conductivity.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lconductivity = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED hcond, diffus_chi
!
!***************************************************************
module Conductivity
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'conductivity.h'
!
  contains
!***********************************************************************
    subroutine register_conductivity()
!
!  Identify version number. 
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_conductivity
!***********************************************************************
    subroutine initialize_conductivity(f,lstarting)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      logical, intent(in) :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_conductivity
!***********************************************************************
    subroutine read_conductivity_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_conductivity_init_pars
!***********************************************************************
    subroutine write_conductivity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_conductivity_init_pars
!***********************************************************************
    subroutine read_conductivity_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_conductivity_run_pars
!***********************************************************************
    subroutine write_conductivity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_conductivity_run_pars
!***********************************************************************
    subroutine pencil_criteria_conductivity()
!
!  All pencils that the conductivity module depends on are specified here.
!
    endsubroutine pencil_criteria_conductivity
!***********************************************************************
    subroutine pencil_interdep_conductivity(lpencil_in)
!
!  Interdependency among pencils from the conductivity module is specified here.
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_conductivity
!***********************************************************************
    subroutine calc_pencils_conductivity(f,p)
!
!  Calculate conductivity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (pencil_case), intent(inout) :: p
!
      call keep_compiler_quiet(f)
      p%hcond = 0.0
      p%diffus_chi = 0.0
!
    endsubroutine calc_pencils_conductivity
!**********************************************************************
    subroutine heat_conductivity(f,df,p)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(in) :: df
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine heat_conductivity
!***********************************************************************
    subroutine dynamical_thermal_diffusion(umax)
!
!  Dummy 
!
      real, intent(in) :: umax
!
      call keep_compiler_quiet(umax)
      call fatal_error('dynamical_thermal_diffusion', 'not implemented yet')
!
    endsubroutine dynamical_thermal_diffusion
!***********************************************************************
    subroutine rprint_conductivity(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_conductivity
!***********************************************************************
  endmodule Conductivity
