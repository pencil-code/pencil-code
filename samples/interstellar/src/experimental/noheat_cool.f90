! $Id$
!
!  Dummy module.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lheat_cool = .false.
!
!***************************************************************
module Heat_Cool
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'heat_cool.h'
!
  contains
!***********************************************************************
    subroutine register_heat_cool()
!
!  08-jul-15/fred: coded
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_heat_cool
!***********************************************************************
    subroutine initialize_heat_cool(f)
!
!  Perform any post-parameter-read initialization eg. set derived
!  parameters
!
!  08-jul-15/fred: coded - dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_heat_cool
!***********************************************************************
    subroutine rprint_heat_cool(lreset,lwrite)
!
!  reads and registers print parameters relevant to heat_cool
!
!   08-jul-15/fred: adapted from magnetic fields
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_heat_cool
!***********************************************************************
    subroutine get_slices_heat_cool(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_heat_cool
!***********************************************************************
    subroutine read_heat_cool_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_heat_cool_init_pars
!***********************************************************************
    subroutine write_heat_cool_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_heat_cool_init_pars
!***********************************************************************
    subroutine read_heat_cool_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_heat_cool_run_pars
!***********************************************************************
    subroutine write_heat_cool_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_heat_cool_run_pars
!***********************************************************************
    subroutine init_heat_cool(f)
!
!  initialise magnetic field; called from start.f90
!  30-jul-2006/tony: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_heat_cool
!***********************************************************************
    subroutine pencil_criteria_heat_cool()
!
!  All pencils that the Heat_Cool module depends on are specified here.
!
!  26-03-05/tony: coded
!
    endsubroutine pencil_criteria_heat_cool
!***********************************************************************
    subroutine heat_cool_before_boundary(f)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  01-aug-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine heat_cool_before_boundary
!***********************************************************************
    subroutine calc_heat_cool_heat_cool(f,df,p,Hmax)
!
!  adapted from calc_heat_cool
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(in) :: df
      type (pencil_case), intent(in) :: p
      real, dimension(nx), intent(in) :: Hmax
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(Hmax)
!
    endsubroutine calc_heat_cool_heat_cool
!***********************************************************************
endmodule heat_cool
