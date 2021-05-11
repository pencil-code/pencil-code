! $Id$
!
!  This modules solves the passive scalar advection equation.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lascalar = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Ascalar
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'ascalar.h'
!
  contains
!***********************************************************************
    subroutine register_ascalar
!
!  Initialise variables which should know that we solve for passive
!  scalar: iacc; increase nvar accordingly.
!
!  6-jul-02/axel: coded
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_ascalar
!***********************************************************************
    subroutine initialize_ascalar(f)
!
!  Perform any necessary post-parameter read initialization.
!
!  24-nov-02/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_ascalar
!***********************************************************************
    subroutine pencil_criteria_ascalar
!
!  All pencils that the Pscalar module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_ascalar
!***********************************************************************
    subroutine pencil_interdep_ascalar(lpencil_in)
!
!  Interdependency among pencils provided by the Pscalar module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_ascalar
!***********************************************************************
    subroutine calc_pencils_ascalar(f,p)
!
!  Calculate Pscalar pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-nov-04/anders: coded
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
    endsubroutine calc_pencils_ascalar
!***********************************************************************
    subroutine init_acc(f)
!
!  Initialise energy; called from start.f90.
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_acc
!***********************************************************************

    subroutine dacc_dt(f,df,p)
!
!  Passive scalar evolution.
!
!   6-jul-02/axel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  :: f,df,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dacc_dt
!***********************************************************************
    subroutine read_ascalar_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_ascalar_init_pars
!***********************************************************************
    subroutine write_ascalar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_ascalar_init_pars
!***********************************************************************
    subroutine read_ascalar_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_ascalar_run_pars
!***********************************************************************
    subroutine write_ascalar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_ascalar_run_pars
!***********************************************************************
    subroutine rprint_ascalar(lreset,lwrite)
!
!  Reads and registers print parameters relevant for passive scalar.
!
!   6-jul-02/axel: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_ascalar
!***********************************************************************
endmodule Ascalar
