! $Id$
!
!  This modules solves the passive scalar advection equation
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lcosmicray = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Cosmicray
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'cosmicray.h'
!
  contains
!***********************************************************************
    subroutine register_cosmicray
!
!  Initialise variables which should know that we solve for active
!  scalar: iecr - the cosmic ray energy density; increase nvar accordingly
!
!  09-oct-03/tony: coded
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_cosmicray
!***********************************************************************
    subroutine initialize_cosmicray(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  09-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      cr_spec=.false.
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_cosmicray
!***********************************************************************
    subroutine read_cosmicray_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_cosmicray_init_pars
!***********************************************************************
    subroutine write_cosmicray_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_cosmicray_init_pars
!***********************************************************************
    subroutine read_cosmicray_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_cosmicray_run_pars
!***********************************************************************
    subroutine write_cosmicray_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_cosmicray_run_pars
!***********************************************************************
    subroutine init_ecr(f)
!
!  initialise magnetic field; called from start.f90
!  We have an init parameter (initlncc) to stear magnetic i.c. independently.
!
!   6-jul-02/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_ecr
!***********************************************************************
    subroutine pencil_criteria_cosmicray
!
!  All pencils that the Cosmicray module depends on are specified here.
!
!  20-11-04/anders: coded
!
!
    endsubroutine pencil_criteria_cosmicray
!***********************************************************************
    subroutine pencil_interdep_cosmicray(lpencil_in)
!
!  Interdependency among pencils provided by the Cosmicray module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_cosmicray
!***********************************************************************
    subroutine calc_pencils_cosmicray(f,p)
!
!  Calculate Cosmicray pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_cosmicray
!***********************************************************************
    subroutine decr_dt(f,df,p)
!
!  Cosmic ray density evolution.
!
!   09-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,df,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine decr_dt
!***********************************************************************
    subroutine rprint_cosmicray(lreset,lwrite)
!
!  Reads and registers print parameters relevant for cosmic rays.
!
!   09-oct-03/tony: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_cosmicray
!***********************************************************************
    subroutine get_slices_cosmicray(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_cosmicray
!***********************************************************************
    subroutine impose_ecr_floor(f)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_ecr_floor
!***********************************************************************
endmodule Cosmicray
