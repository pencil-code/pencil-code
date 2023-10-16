! $Id$
!
!  Dummy module for mean-field contributions
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagn_mf_demfdt = .false.
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
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'meanfield_demfdt.h'
!
  contains
!***********************************************************************
    subroutine register_magn_mf_demfdt()
!
!  Dummy routine
!
    endsubroutine register_magn_mf_demfdt
!***********************************************************************
    subroutine initialize_magn_mf_demfdt(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_magn_mf_demfdt
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
    subroutine pencil_criteria_magn_mf_demfdt()
!
!  Dummy routine
!
    endsubroutine pencil_criteria_magn_mf_demfdt
!***********************************************************************
    subroutine pencil_interdep_magn_mf_demfdt(lpencil_in)
!
!  Dummy routine
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_magn_mf_demfdt
!***********************************************************************
    subroutine calc_pencils_magn_mf_demfdt(f,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_magn_mf_demfdt
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
    subroutine calc_diagnostics_dt_meanfield(f)
!
      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)

    endsubroutine calc_diagnostics_dt_meanfield
!***********************************************************************
    subroutine read_magn_mf_demfdt_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_magn_mf_demfdt_init_pars
!***********************************************************************
    subroutine write_magn_mf_demfdt_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_magn_mf_demfdt_init_pars
!***********************************************************************
    subroutine read_magn_mf_demfdt_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_magn_mf_demfdt_run_pars
!***********************************************************************
    subroutine write_magn_mf_demfdt_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_magn_mf_demfdt_run_pars
!***********************************************************************
    subroutine rprint_magn_mf_demfdt(lreset,lwrite)
!
!  Dummy routine
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_magn_mf_demfdt
!***********************************************************************
endmodule Magnetic_meanfield_demfdt
