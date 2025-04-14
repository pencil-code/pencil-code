! $Id$
!
!  This module takes care of everything related to dust density.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldustdensity = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED rhod(ndustspec); udropav(3), rhodsum, glnrhodsum(3)
!
!***************************************************************
module Dustdensity
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'dustdensity.h'
!
  contains
!***********************************************************************
    subroutine register_dustdensity
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_dustdensity
!***********************************************************************
    subroutine initialize_dustdensity(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      nd_spec=.false.
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_dustdensity
!***********************************************************************
    subroutine init_nd(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_nd
!***********************************************************************
    subroutine pencil_criteria_dustdensity
!
    endsubroutine pencil_criteria_dustdensity
!***********************************************************************
    subroutine pencil_interdep_dustdensity(lpencil_in)
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_dustdensity
!***********************************************************************
    subroutine dustdensity_after_boundary(f)
!
      use General, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine dustdensity_after_boundary
!***********************************************************************
    subroutine dustdensity_before_boundary(f)
!
      use General, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine dustdensity_before_boundary
!***********************************************************************
    subroutine calc_pencils_dustdensity(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f, p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_dustdensity
!***********************************************************************
    subroutine dndmd_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dndmd_dt
!***********************************************************************
    subroutine calc_diagnostics_dustdensity(f,p)
!
      type (pencil_case) :: p
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_diagnostics_dustdensity
!***********************************************************************
    subroutine read_dustdensity_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_dustdensity_init_pars
!***********************************************************************
    subroutine write_dustdensity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_dustdensity_init_pars
!***********************************************************************
    subroutine read_dustdensity_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_dustdensity_run_pars
!***********************************************************************
    subroutine write_dustdensity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_dustdensity_run_pars
!***********************************************************************
    subroutine rprint_dustdensity(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_dustdensity
!***********************************************************************
    subroutine get_slices_dustdensity(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_dustdensity
!***********************************************************************
   subroutine impose_dustdensity_floor(f)
!
     real, dimension (mx,my,mz,mfarray) :: f
!
     call keep_compiler_quiet(f)
!
   endsubroutine impose_dustdensity_floor
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par
    endsubroutine pushpars2c
!***********************************************************************
endmodule Dustdensity
