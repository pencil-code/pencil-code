! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lheatflux = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Heatflux
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include 'heatflux.h'
!
  contains
!***********************************************************************
    subroutine register_heatflux()
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_heatflux
!***********************************************************************
    subroutine initialize_heatflux(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_heatflux
!***********************************************************************
    subroutine finalize_heatflux(f)
!
!  Called right before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine finalize_heatflux
!***********************************************************************
    subroutine init_heatflux(f)
!
!  initialise heatflux condition; called from start.f90
!
!  06-oct-2003/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_heatflux
!***********************************************************************
    subroutine pencil_criteria_heatflux()
!
!  All pencils that this heatflux module depends on are specified here.
!
!  18-07-06/tony: coded
!
    endsubroutine pencil_criteria_heatflux
!***********************************************************************
    subroutine pencil_interdep_heatflux(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_heatflux
!***********************************************************************
    subroutine calc_pencils_heatflux(f,p)
!
!  Calculate Heatflux pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_heatflux
!***********************************************************************
    subroutine dheatflux_dt(f,df,p)
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dheatflux_dt
!***********************************************************************
    subroutine read_heatflux_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_heatflux_init_pars
!***********************************************************************
    subroutine write_heatflux_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_heatflux_init_pars
!***********************************************************************
    subroutine read_heatflux_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_heatflux_run_pars
!***********************************************************************
    subroutine write_heatflux_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_heatflux_run_pars
!***********************************************************************
    subroutine rprint_heatflux(lreset,lwrite)
!
!  Reads and registers print parameters relevant to heatflux.
!
!  06-oct-03/tony: coded
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_heatflux
!***********************************************************************
    subroutine get_slices_heatflux(f,slices)
!
!  Write slices for animation of Heatflux variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_heatflux
!***********************************************************************
endmodule Heatflux
