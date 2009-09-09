! $Id$
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lchemistry = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Chemistry
!
  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'chemistry.h'
!
  contains
!***********************************************************************
    subroutine register_chemistry()
!
    endsubroutine register_chemistry
!***********************************************************************
    subroutine initialize_chemistry(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_chemistry
!***********************************************************************
    subroutine init_chemistry(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_chemistry
!***********************************************************************
    subroutine pencil_criteria_chemistry()
!
    endsubroutine pencil_criteria_chemistry
!***********************************************************************
    subroutine pencil_interdep_chemistry(lpencil_in)
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_chemistry
!***********************************************************************
    subroutine calc_pencils_chemistry(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_chemistry
!***********************************************************************
    subroutine calc_for_chem_mixture(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_for_chem_mixture
!**********************************************************************
    subroutine dchemistry_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dchemistry_dt
!***********************************************************************
    subroutine read_chemistry_init_pars(unit,iostat)
!
      integer :: unit
      integer, optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_chemistry_init_pars
!***********************************************************************
   subroutine write_chemistry_init_pars(unit)
!
      integer :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_chemistry_init_pars
!***********************************************************************
    subroutine read_chemistry_run_pars(unit,iostat)
!
      integer :: unit
      integer, optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_chemistry_run_pars
!***********************************************************************
    subroutine write_chemistry_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_chemistry_run_pars
!***********************************************************************
    subroutine rprint_chemistry(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_chemistry
!***********************************************************************
    subroutine get_slices_chemistry(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices)
!
    endsubroutine get_slices_chemistry
!********************************************************************
    subroutine bc_nscbc_subin_x(f,df,topbot,val)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (mcom), optional :: val
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(val)
!
    endsubroutine bc_nscbc_subin_x
!********************************************************************
    subroutine bc_nscbc_nref_subout_x(f,df,topbot,nscbc_sigma_out)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real :: nscbc_sigma_out
!     
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(nscbc_sigma_out)
!
    endsubroutine bc_nscbc_nref_subout_x
!***********************************************************************
    subroutine bc_nscbc_nref_subout_y(f,df,topbot,nscbc_sigma_out)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real :: nscbc_sigma_out
!     
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(nscbc_sigma_out)
!
    endsubroutine bc_nscbc_nref_subout_y
!***********************************************************************
    subroutine bc_nscbc_nref_subout_z(f,df,topbot,nscbc_sigma_out)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real :: nscbc_sigma_out
!     
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(nscbc_sigma_out)
!
    endsubroutine bc_nscbc_nref_subout_z
!***********************************************************************
    subroutine chemistry_clean_up()
!
    endsubroutine chemistry_clean_up
!***********************************************************************
    subroutine jacobn(f,jacob)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,nchemspec,nchemspec) :: jacob
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(jacob(1,1,1,1,1))
!
    endsubroutine jacobn
!***********************************************************************
endmodule Chemistry
