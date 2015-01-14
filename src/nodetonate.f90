! $Id$
!
!  This module is dummy.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldetonate = .false.
!
!***************************************************************
module Detonate
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id
!
  implicit none
!
  include 'detonate.h'
!
  contains
!***********************************************************************
    subroutine register_detonate()
!
!  Set up indices for variables.
!
!  13-feb-14/ccyang: coded
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_detonate
!***********************************************************************
    subroutine initialize_detonate(f)
!
!  Initializes module specific variables for later use.
!
!  13-feb-14/ccyang: dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_detonate
!***********************************************************************
    subroutine read_detonate_run_pars(unit, iostat)
!
!  Reads the runtime parameters.
!
!  13-feb-14/ccyang: dummy
!
      include 'unit.h'
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_detonate_run_pars
!***********************************************************************
    subroutine write_detonate_run_pars(unit)
!
!  Writes the runtime parameters.
!
!  13-feb-14/ccyang: dummy
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_detonate_run_pars
!***********************************************************************
    subroutine rprint_detonate(lreset, lwrite)
!
!  Reads and registers print parameters.
!
!  13-feb-14/ccyang: dummy
!
      logical, intent(in) :: lreset
      logical, intent(in), optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_detonate
!***********************************************************************
    subroutine detonate_before_boundary(f)
!
!  Possibility to modify the f before boundary communications.
!
!  14-feb-14/ccyang: dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine detonate_before_boundary
!***********************************************************************
endmodule Detonate
