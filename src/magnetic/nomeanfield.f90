! $Id$
!
!  This modules solves mean-field contributions to both the
!  induction and the momentum equations.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagn_mf = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED mf_EMF(3); mf_EMFdotB
!
!***************************************************************
module Magnetic_meanfield
!
  use Cparam
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'meanfield.h'
!
  contains
!***********************************************************************
    subroutine register_magn_mf()
!
!  Dummy routine
!
    endsubroutine register_magn_mf
!***********************************************************************
    subroutine initialize_magn_mf(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_magn_mf
!***********************************************************************
    subroutine init_aa_mf(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_aa_mf
!***********************************************************************
    subroutine pencil_criteria_magn_mf()
!
!  Dummy routine
!
    endsubroutine pencil_criteria_magn_mf
!***********************************************************************
    subroutine pencil_interdep_magn_mf(lpencil_in)
!
!  Dummy routine
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_magn_mf
!***********************************************************************
    subroutine calc_pencils_magn_mf(f,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lpencil(i_mf_EMF)) p%mf_EMF=0.0
      if (lpencil(i_mf_EMFdotB)) p%mf_EMFdotB=0.0
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_magn_mf
!***********************************************************************
    subroutine daa_dt_meanfield(f,df,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  :: p
      intent(inout)  :: f,df
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine daa_dt_meanfield
!***********************************************************************
    subroutine  calc_diagnostics_meanfield(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine  calc_diagnostics_meanfield
!***********************************************************************
    subroutine read_magn_mf_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_magn_mf_init_pars
!***********************************************************************
    subroutine write_magn_mf_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_magn_mf_init_pars
!***********************************************************************
    subroutine read_magn_mf_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_magn_mf_run_pars
!***********************************************************************
    subroutine write_magn_mf_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_magn_mf_run_pars
!***********************************************************************
    subroutine rprint_magn_mf(lreset,lwrite)
!
!  Dummy routine
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_magn_mf
!***********************************************************************
    subroutine pc_aasb_const_alpha(f,topbot,j)

    real, dimension(:,:,:,:) :: f
    integer :: j, topbot

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)

    endsubroutine pc_aasb_const_alpha
!***********************************************************************
    subroutine meanfield_after_boundary(f)
!
!  precalculate parameters that are new at each timestep,
!  but the same for all pencils
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine meanfield_after_boundary
!***********************************************************************
endmodule Magnetic_meanfield

