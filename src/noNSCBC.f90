! $Id$
!
!  Module for NSCBC boundary conditions.
!  To be included from boundcond.f90.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lnscbc = .false.
!
!***************************************************************
module NSCBC
!
  use Cparam
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'NSCBC.h'
!
  contains
!***********************************************************************
    subroutine nscbc_boundtreat(f,df)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df

      intent(inout) :: f
      intent(inout) :: df
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
!
    endsubroutine nscbc_boundtreat
!***********************************************************************
    subroutine read_NSCBC_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_NSCBC_init_pars
!***********************************************************************
    subroutine write_NSCBC_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_NSCBC_init_pars
!***********************************************************************
    subroutine read_NSCBC_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_NSCBC_run_pars
!***********************************************************************
    subroutine write_NSCBC_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_NSCBC_run_pars
!***********************************************************************
    subroutine NSCBC_clean_up
!
    endsubroutine NSCBC_clean_up
!***********************************************************************
endmodule NSCBC
