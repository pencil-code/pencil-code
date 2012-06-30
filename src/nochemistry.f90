! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lchemistry = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED Ywater
!
!***************************************************************
module Chemistry
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use EquationOfState
!
  implicit none
!
  real :: Rgas
  logical :: lchemistry_diag=.false.
  logical :: lreactions=.false.
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
!***********************************************************************
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
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_chemistry
!***********************************************************************
    subroutine chemspec_normalization(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine chemspec_normalization
!***********************************************************************
    subroutine get_slices_chemistry(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_chemistry
!***********************************************************************
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
      call keep_compiler_quiet(present(val))
!
    endsubroutine bc_nscbc_subin_x
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
    subroutine get_mu1_slice(slice,grad_slice,index,sgn,direction)
!
! For the NSCBC boudary conditions the slice of mu1 at the boundary
! is required.
!
! 2009.12.10: Nils Erland L. Haugen (coded)
!
      real, dimension(ny,nz), intent(out) :: slice, grad_slice
      integer, intent(in) :: index, sgn,direction
!
      call keep_compiler_quiet(slice)
      call keep_compiler_quiet(grad_slice)
      call keep_compiler_quiet(index,sgn,direction)
!
    end subroutine get_mu1_slice
!***********************************************************************
    subroutine get_gamma_slice(slice,index,dir)
!
!  Get a 2D slice of gamma
!
!  2009.12.10: Nils Erland L. Haugen (coded)
!
      real, dimension (:,:), intent(out)  :: slice
      integer, intent(in) :: index,dir
!
      call keep_compiler_quiet(slice)
      call keep_compiler_quiet(index)
      call keep_compiler_quiet(dir)
      !
    endsubroutine get_gamma_slice
!***********************************************************************
    subroutine get_cs2_slice(slice,index,dir)
!
!  Get a 2D slice of cs2
!
!  2009.12.10: Nils Erland L. Haugen (coded)
!
      real, dimension (:,:), intent(out)  :: slice
      integer, intent(in) :: index,dir
!
      call keep_compiler_quiet(slice)
      call keep_compiler_quiet(index)
      call keep_compiler_quiet(dir)
      !
    endsubroutine get_cs2_slice
!***********************************************************************
   subroutine get_cs2_full(cs2_full)
!
      real, dimension (mx,my,mz) :: cs2_full
!
      intent(out) :: cs2_full
!
      call keep_compiler_quiet(cs2_full)
!
    endsubroutine get_cs2_full
!***********************************************************************
    subroutine get_gamma_full(gamma_full)
!
      real, dimension (mx,my,mz) :: gamma_full
!
      intent(out) :: gamma_full
!
      call keep_compiler_quiet(gamma_full)
!
    endsubroutine get_gamma_full
!***********************************************************************
    subroutine get_RHS_Y_full(RHS_Y)
!
      real, dimension (mx,my,mz,nchemspec) :: RHS_Y
!
      intent(out) :: RHS_Y
!
      call keep_compiler_quiet(RHS_Y)
!
    endsubroutine get_RHS_Y_full
!***********************************************************************
    subroutine  write_net_reaction
    endsubroutine  write_net_reaction
!***********************************************************************
    subroutine get_reac_rate(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine get_reac_rate
!***********************************************************************
endmodule Chemistry
