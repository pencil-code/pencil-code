! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lcosmicrayflux = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!
!***************************************************************
module Cosmicrayflux
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'cosmicrayflux.h'
!
  contains
!***********************************************************************
    subroutine register_cosmicrayflux()
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_cosmicrayflux
!***********************************************************************
    subroutine initialize_cosmicrayflux(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_cosmicrayflux
!***********************************************************************
    subroutine read_cosmicrayflux_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_cosmicrayflux_init_pars
!***********************************************************************
    subroutine write_cosmicrayflux_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_cosmicrayflux_init_pars
!***********************************************************************
    subroutine read_cosmicrayflux_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_cosmicrayflux_run_pars
!***********************************************************************
    subroutine write_cosmicrayflux_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_cosmicrayflux_run_pars
!***********************************************************************
    subroutine init_fcr(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_fcr
!***********************************************************************
    subroutine pencil_criteria_cosmicrayflux()
!
!   All pencils that the CosmicrayFlux module depends on are specified here.
!
    endsubroutine pencil_criteria_cosmicrayflux
!***********************************************************************
    subroutine pencil_interdep_cosmicrayflux(lpencil_in)
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_cosmicrayflux
!***********************************************************************
    subroutine calc_pencils_cosmicrayflux(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_cosmicrayflux
!***********************************************************************
    subroutine dfcr_dt(f,df,p)
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
    endsubroutine dfcr_dt
!***********************************************************************
    subroutine rprint_cosmicrayflux(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_cosmicrayflux
!***********************************************************************
endmodule Cosmicrayflux
