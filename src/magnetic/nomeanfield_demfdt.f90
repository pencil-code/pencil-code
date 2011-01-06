! $Id$
!
!  Dummy module for mean-field contributions
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagnetic_mf_demfdt = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED
!
!***************************************************************
module Magnetic_meanfield_demfdt
!
  use Cdata
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'meanfield_demfdt.h'
!
  contains
!***********************************************************************
    subroutine register_magnetic_mf_demfdt()
!
!  Dummy routine
!
    endsubroutine register_magnetic_mf_demfdt
!***********************************************************************
    subroutine initialize_magnetic_mf_demfdt(f,lstarting)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_magnetic_mf_demfdt
!***********************************************************************
    subroutine init_aa_mf_demfdt(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_aa_mf_demfdt
!***********************************************************************
    subroutine pencil_criteria_magnetic_mf_demfdt()
!
!  Dummy routine
!
    endsubroutine pencil_criteria_magnetic_mf_demfdt
!***********************************************************************
    subroutine pencil_interdep_magnetic_mf_demfdt(lpencil_in)
!
!  Dummy routine
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_magnetic_mf_demfdt
!***********************************************************************
    subroutine calc_pencils_magnetic_mf_demfdt(f,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_magnetic_mf_demfdt
!***********************************************************************
    subroutine demf_dt_meanfield(f,df,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine demf_dt_meanfield
!***********************************************************************
    subroutine read_magnetic_mf_demfdt_init_pars(unit,iostat)
!
!  Dummy routine
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_magnetic_mf_demfdt_init_pars
!***********************************************************************
    subroutine write_magnetic_mf_demfdt_init_pars(unit)
!
!  Dummy routine
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_magnetic_mf_demfdt_init_pars
!***********************************************************************
    subroutine read_magnetic_mf_demfdt_run_pars(unit,iostat)
!
!  Dummy routine
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_magnetic_mf_demfdt_run_pars
!***********************************************************************
    subroutine write_magnetic_mf_demfdt_run_pars(unit)
!
!  Dummy routine
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_magnetic_mf_demfdt_run_pars
!***********************************************************************
    subroutine rprint_magnetic_mf_demfdt(lreset,lwrite)
!
!  Dummy routine
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_magnetic_mf_demfdt
!***********************************************************************
endmodule Magnetic_meanfield_demfdt
