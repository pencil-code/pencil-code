! $Id$
!
!  This module is used both for the initial condition and during run time.
!  It contains dlnrhon_dt and init_lnrhon, among other auxiliary routines.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lneutraldensity = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Neutraldensity
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'neutraldensity.h'
!
  integer :: idiag_rhonm=0
!
  contains
!***********************************************************************
    subroutine register_neutraldensity()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrhon; increase nvar accordingly.
!
!  18-mar-03/axel: adapted from neutraldensity
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_neutraldensity
!***********************************************************************
    subroutine initialize_neutraldensity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-mar-03/axel: adapted from neutraldensity
!
!  do nothing
!
    endsubroutine initialize_neutraldensity
!***********************************************************************
    subroutine init_lnrhon(f)
!
!  Initialise lnrhon; called from start.f90.
!
!  18-mar-03/axel: adapted from neutraldensity
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_lnrhon
!***********************************************************************
    subroutine pencil_criteria_neutraldensity()
!
!  All pencils that the Neutraldensity module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_neutraldensity
!***********************************************************************
    subroutine pencil_interdep_neutraldensity(lpencil_in)
!
!  Interdependency among pencils provided by the Neutraldensity module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_neutraldensity
!***********************************************************************
    subroutine calc_pencils_neutraldensity(f,p)
!
!  Calculate Neutraldensity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  13-nov-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f, p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_neutraldensity
!***********************************************************************
    subroutine dlnrhon_dt(f,df,p)
!
!  Continuity equation.
!  Calculate dlnrhon/dt = - u.gradlnrhon - divud
!
!  18-mar-03/axel: adapted from neutraldensity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dlnrhon_dt
!***********************************************************************
    subroutine read_neutraldensity_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_neutraldensity_init_pars
!***********************************************************************
    subroutine write_neutraldensity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_neutraldensity_init_pars
!***********************************************************************
    subroutine read_neutraldensity_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_neutraldensity_run_pars
!***********************************************************************
    subroutine write_neutraldensity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_neutraldensity_run_pars
!***********************************************************************
    subroutine rprint_neutraldensity(lreset,lwrite)
!
!  Reads and registers print parameters relevant for neutral density.
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write column where which neutral density variable is stored.
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'ilnrhon=',ilnrhon
      endif
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_neutraldensity
!***********************************************************************
endmodule Neutraldensity
