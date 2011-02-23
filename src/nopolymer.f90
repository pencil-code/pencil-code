! $Id$
!
!  This modules deals with all aspects of polymers.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpolymer = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!
!***************************************************************
module Polymer
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'record_types.h'
  include 'polymer.h'
!
  contains
!***********************************************************************
    subroutine register_polymer()
!
!  Initialise variables which should know that we solve for the vector
!  potential: ipoly, etc; increase nvar accordingly.
!
!  14-Aug-08 : Dhruba
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_polymer
!***********************************************************************
    subroutine initialize_polymer(f,lstarting)
!
!  Perform any post-parameter-read initialization.
!
!  14-aug-08/dhruba: initialize polymer field (dummy at present)
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_polymer
!***********************************************************************
    subroutine init_poly(f)
!
!  Initialise polymer field.
!
!   14-aug-2008/dhruba: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_poly
!***********************************************************************
    subroutine pencil_criteria_polymer()
!
!   All pencils that the Polymer module depends on are specified here.
!
    endsubroutine pencil_criteria_polymer
!***********************************************************************
    subroutine pencil_interdep_polymer(lpencil_in)
!
!  Interdependency among pencils from the Polymer module is specified here.
!
!  18-aug-2008/dhruba: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_polymer
!***********************************************************************
    subroutine calc_pencils_polymer(f,p)
!
!  Calculate Magnetic pencils.
!  Most basic pencils should come first, as others may depend on them.
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
    endsubroutine calc_pencils_polymer
!***********************************************************************
    subroutine dpoly_dt(f,df,p)
!
!  Polymer evolution.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dpoly_dt
!***********************************************************************
    subroutine calc_polymer_after_boundary(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_polymer_after_boundary(f)
!***********************************************************************
    subroutine read_polymer_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_polymer_init_pars
!***********************************************************************
    subroutine write_polymer_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_polymer_init_pars
!***********************************************************************
    subroutine read_polymer_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_polymer_run_pars
!***********************************************************************
    subroutine write_polymer_run_pars(unit)
!
     integer, intent(in) :: unit
!
     call keep_compiler_quiet(unit)
!
    endsubroutine write_polymer_run_pars
!***********************************************************************
    subroutine get_slices_polymer(f,slices)
!
!  Write slices for animation of polymeric variables.
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_polymer
!***********************************************************************
    subroutine rprint_polymer(lreset,lwrite)
!
!  Reads and registers print parameters relevant for polymer.
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_polymer
!***********************************************************************
endmodule Polymer
