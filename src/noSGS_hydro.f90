! $Id$
!
!  This module is the dummy for the SGS_hydro module
!  in which e.g the SGS viscous force or viscous heat is
!  computed.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lSGS_hydro = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!! PENCILS PROVIDED fSGS(3)
! PENCILS PROVIDED SGS_heat   
!; diffus_total
!
!***************************************************************
module SGS_hydro
!
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'SGS_hydro.h'
!
  contains
!***********************************************************************
    subroutine register_SGS_hydro
!
!  19-nov-02/tony: coded
!
      use Messages, only: svn_id
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_SGS_hydro
!***********************************************************************
    subroutine initialize_SGS_hydro
!
    endsubroutine initialize_SGS_hydro
!***********************************************************************
    subroutine read_SGS_hydro_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_SGS_hydro_run_pars
!***********************************************************************
    subroutine write_SGS_hydro_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_SGS_hydro_run_pars
!***********************************************************************
    subroutine rprint_SGS_hydro(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_SGS_hydro
!***********************************************************************
    subroutine pencil_criteria_SGS_hydro
!
!  All pencils that the SGS_hydro module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_SGS_hydro
!***********************************************************************
    subroutine pencil_interdep_SGS_hydro(lpencil_in)
!
!  Interdependency among pencils from the SGS_hydro module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_SGS_hydro
!***********************************************************************
    subroutine calc_pencils_SGS_hydro(f,p)
!
!  Calculate SGS_hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
!  define SGS heat
!
      if (lpencil(i_SGS_heat)) p%SGS_heat=0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_SGS_hydro
!***********************************************************************
    subroutine SGS_hydro_after_boundary(f)

      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine SGS_hydro_after_boundary
!***********************************************************************
    subroutine calc_SGS_hydro_heat(df,p,Hmax)
!
      real, dimension (mx,my,mz,mvar)    :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax
!
      intent(in) :: df,p,Hmax
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(Hmax)
!
    endsubroutine calc_SGS_hydro_heat
!***********************************************************************
    subroutine calc_SGS_hydro_force(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar)    :: df
      type (pencil_case) :: p
!
      intent (in) :: df,p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_SGS_hydro_force
!***********************************************************************
endmodule SGS_hydro
