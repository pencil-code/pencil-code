! $Id$

!  test perturbation method

module TestPerturb

  implicit none

  real :: dummy

  include 'testperturb.h'

  namelist /testperturb_init_pars/ &
      dummy

  namelist /testperturb_run_pars/ &
      dummy

  contains

!***********************************************************************
    subroutine register_testperturb()
!
!  Dummy routine
!
    endsubroutine register_testperturb
!***********************************************************************
    subroutine initialize_testperturb()
!
!  Dummy routine
!
    endsubroutine initialize_testperturb
!***********************************************************************
    subroutine read_testperturb_init_pars(unit,iostat)
!
!  Dummy routine
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
    endsubroutine read_testperturb_init_pars
!***********************************************************************
    subroutine write_testperturb_init_pars(unit)
!
!  Dummy routine
!
      integer, intent(in) :: unit
!
    endsubroutine write_testperturb_init_pars
!***********************************************************************
    subroutine read_testperturb_run_pars(unit,iostat)
!
!  Dummy routine
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
    endsubroutine read_testperturb_run_pars
!***********************************************************************
    subroutine write_testperturb_run_pars(unit)
!
!  Dummy routine
!
      integer, intent(in) :: unit
!
    endsubroutine write_testperturb_run_pars
!***********************************************************************
    subroutine testperturb_begin(f,df)
!
!  Dummy routine
!
      use Cparam, only: mx,my,mz,mfarray,mvar
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      call keep_compiler_quiet(f,df)
    endsubroutine testperturb_begin
!***********************************************************************
    subroutine testperturb_finalize(f)
!
!  Dummy routine
!
      use Cparam, only: mx,my,mz,mfarray,mvar
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
    endsubroutine testperturb_finalize
!***********************************************************************
    subroutine rprint_testperturb(lreset,lwrite)
!
!  Dummy routine
!
      logical :: lreset
      logical, optional :: lwrite
!
    endsubroutine rprint_testperturb
!***********************************************************************
  endmodule TestPerturb
