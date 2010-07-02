! $Id: nodensity.f90 12731 2009-12-30 06:54:37Z ajohan@strw.leidenuniv.nl $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .false.
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED rho; lnrho; rho1; glnrho(3); del2rho; del2lnrho
! PENCILS PROVIDED hlnrho(3,3); grho(3); glnrho2
! PENCILS PROVIDED del6lnrho; uij5glnrho(3); uglnrho; ugrho; sglnrho(3)
!
!***************************************************************
module Density
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  logical :: lcalc_glnrhomean=.false.,lupw_lnrho=.false.
  real, dimension (nz,3) :: glnrhomz

  include 'density.h'
!
  contains
!***********************************************************************
    subroutine register_density()
!
      use Sub
!
      if (lroot) call svn_id( &
          "$Id: nodensity.f90 12731 2009-12-30 06:54:37Z ajohan@strw.leidenuniv.nl $")
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      use EquationOfState, only: select_eos_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Tell the equation of state that we're here and we don't have a
!  variable => isochoric (constant density).
!
      call select_eos_variable('lnrho',-1)
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_density
!***********************************************************************
    subroutine init_lnrho(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_lnrho
!***********************************************************************
    subroutine calc_ldensity_pars(f)

      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
  endsubroutine calc_ldensity_pars
!***********************************************************************
    subroutine pencil_criteria_density()
!
!  All pencils that the Density module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_density
!***********************************************************************
    subroutine pencil_interdep_density(lpencil_in)
!
!  Interdependency among pencils from the Density module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_density
!***********************************************************************
    subroutine calc_pencils_density(f,p)
!
!  Calculate Density pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use EquationOfState, only: lnrho0, rho0
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! rho
      if (lpencil(i_rho)) p%rho=rho0
! lnrho
      if (lpencil(i_lnrho)) p%lnrho=lnrho0
! rho1
      if (lpencil(i_rho1)) p%rho1=1/rho0
! glnrho
      if (lpencil(i_glnrho)) p%glnrho=0.0
! grho
      if (lpencil(i_grho)) p%grho=0.0
! del6lnrho
      if (lpencil(i_del6lnrho)) p%del6lnrho=0.0
! hlnrho
      if (lpencil(i_hlnrho)) p%hlnrho=0.0
! sglnrho
      if (lpencil(i_sglnrho)) p%sglnrho=0.0
! uglnrho
      if (lpencil(i_uglnrho)) p%uglnrho=0.0
! ugrho
      if (lpencil(i_ugrho)) p%ugrho=0.0
! uij5glnrho
      if (lpencil(i_uij5glnrho)) p%uij5glnrho=0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_density
!***********************************************************************
    subroutine density_before_boundary(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine density_before_boundary
!***********************************************************************
    subroutine dlnrho_dt(df,p)
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: df,p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine read_density_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_density_init_pars
!***********************************************************************
    subroutine write_density_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_density_init_pars
!***********************************************************************
    subroutine read_density_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_density_run_pars
!***********************************************************************
    subroutine write_density_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_density_run_pars
!***********************************************************************
    subroutine rprint_density(lreset,lwrite)
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_density
!***********************************************************************
    subroutine get_slices_density(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_density
!***********************************************************************
endmodule Density
