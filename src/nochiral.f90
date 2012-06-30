! $Id$
!
!  This modules solves two reactive scalar advection equations
!  This is used for modeling the spatial evolution of left and
!  right handed aminoacids.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lchiral = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Chiral
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'chiral.h'
!
  integer :: dummy           ! We cannot define empty namelists
!
  namelist /chiral_init_pars/ dummy
  namelist /chiral_run_pars/  dummy
!
  integer :: idiag_XX_chiralmax=0, idiag_YY_chiralmax=0
!
  contains
!***********************************************************************
    subroutine register_chiral()
!
!  Initialise variables which should know that we solve for passive
!  scalar: iXX_chiral and iYY_chiral; increase nvar accordingly
!
!  28-may-04/axel: adapted from pscalar
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_chiral
!***********************************************************************
    subroutine initialize_chiral(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  28-may-04/axel: adapted from pscalar
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_chiral
!***********************************************************************
    subroutine init_chiral(f)
!
!  initialise passive scalar field; called from start.f90
!
!  28-may-04/axel: adapted from pscalar
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_chiral
!***********************************************************************
    subroutine pencil_criteria_chiral()
!
!  All pencils that the Chiral module depends on are specified here.
!
!  21-11-04/anders: coded
!
    endsubroutine pencil_criteria_chiral
!***********************************************************************
    subroutine pencil_interdep_chiral(lpencil_in)
!
!  Interdependency among pencils provided by the Chiral module
!  is specified here.
!
!  21-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_chiral
!***********************************************************************
    subroutine calc_pencils_chiral(f,p)
!
!  Calculate Chiral pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  21-11-04/anders: coded
!
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_chiral
!***********************************************************************
    subroutine dXY_chiral_dt(f,df,p)
!
!  passive scalar evolution
!  calculate dc/dt=-uu.gcc + pscaler_diff*[del2cc + glnrho.gcc]
!
!  28-may-04/axel: adapted from pscalar
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
    endsubroutine dXY_chiral_dt
!***********************************************************************
    subroutine chiral_before_boundary(f)
!
!  dummy routine
!
!   4-jul-09/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
    endsubroutine chiral_before_boundary
!***********************************************************************
    subroutine read_chiral_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_chiral_init_pars
!***********************************************************************
    subroutine write_chiral_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_chiral_init_pars
!***********************************************************************
    subroutine read_chiral_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_chiral_run_pars
!***********************************************************************
    subroutine write_chiral_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_chiral_run_pars
!***********************************************************************
    subroutine rprint_chiral(lreset,lwrite)
!
!  reads and registers print parameters relevant for chirality
!
!  28-may-04/axel: adapted from pscalar
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_XX_chiralmax=0
        idiag_YY_chiralmax=0
      endif
!
!  write column where which chiral variable is stored
!
      if (lwr) then
        write(3,*) 'iXX_chiral=0'
        write(3,*) 'iYY_chiral=0'
      endif
!
    endsubroutine rprint_chiral
!***********************************************************************
    subroutine get_slices_chiral(f,slices)
!
!  Write slices for animation of chiral variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_chiral
!***********************************************************************
endmodule Chiral
