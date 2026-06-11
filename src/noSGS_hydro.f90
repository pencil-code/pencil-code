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
    subroutine read_SGS_hydro_run_pars(iomsg)
!
      character(LEN=*), intent(out) :: iomsg
!
      iomsg=""
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
!TP: calc pencils of this module is not used anymore!
!    subroutine calc_pencils_SGS_hydro(f,p)
!!
!!  Calculate SGS_hydro pencils.
!!  Most basic pencils should come first, as others may depend on them.
!!
!!  20-11-04/anders: coded
!!
!      real, dimension (mx,my,mz,mfarray) :: f
!      type (pencil_case) :: p
!!
!      intent(in) :: f
!      intent(inout) :: p
!!
!!  define SGS heat
!!
!      if (lpencil(i_SGS_heat)) p%SGS_heat=0.0
!!
!      call keep_compiler_quiet(f)
!!
!    endsubroutine calc_pencils_SGS_hydro
!!***********************************************************************
    subroutine SGS_hydro_after_boundary(f)

      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine SGS_hydro_after_boundary
!***********************************************************************
    subroutine calc_SGS_hydro_force(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar)    :: df
      type (pencil_case) :: p
!
      intent (in) :: df,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_SGS_hydro_force
!***********************************************************************
    subroutine calc_diagnostics_SGS_hydro(p)

      type (pencil_case) :: p

      call keep_compiler_quiet(p)
!
    endsubroutine calc_diagnostics_SGS_hydro
!***********************************************************************
endmodule SGS_hydro
