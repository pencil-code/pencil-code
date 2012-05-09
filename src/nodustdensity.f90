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
! PENCILS PROVIDED rhod(ndustspec); udropav(3), rhodsum, grhodsum(3)
!
!***************************************************************
module Dustdensity
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'dustdensity.h'
!
  contains
!***********************************************************************
    subroutine register_dustdensity()
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
      call keep_compiler_quiet(f)
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
    subroutine pencil_criteria_dustdensity()
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
    subroutine read_dustdensity_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
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
    subroutine read_dustdensity_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
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
    subroutine redist_mdbins(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine redist_mdbins
!***********************************************************************
    subroutine null_dust_vars(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine null_dust_vars
!***********************************************************************
    subroutine rprint_dustdensity(lreset,lwrite)
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lwr) then
        write(3,*) 'ind=',ind
      endif
!
      call keep_compiler_quiet(lreset)
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
  subroutine dustspec_normalization(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
   endsubroutine dustspec_normalization
!***********************************************************************
   subroutine impose_dustdensity_floor(f)
!
     real, dimension (mx,my,mz,mfarray) :: f
!
     call keep_compiler_quiet(f)
!
   endsubroutine impose_dustdensity_floor
!***********************************************************************
endmodule Dustdensity
