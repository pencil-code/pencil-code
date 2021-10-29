! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .false.
! CPARAM logical, parameter :: lanelastic = .false.
! CPARAM logical, parameter :: lboussinesq = .false.
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED rho; lnrho; rho1; glnrho(3); del2rho; del2lnrho
! PENCILS PROVIDED hlnrho(3,3); grho(3); glnrho2
! PENCILS PROVIDED del6lnrho; uij5glnrho(3); uglnrho; ugrho; sglnrho(3)
! PENCILS PROVIDED ekin; transprho
! PENCILS PROVIDED glnrhos(3)
! PENCILS PROVIDED totenergy_rel
!
!***************************************************************
module Density
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  logical :: lcalc_lnrhomean=.false., lcalc_glnrhomean=.false.,lupw_lnrho=.false.
  real, dimension (mz) :: lnrhomz
  real, dimension(3) :: beta_glnrho_global=0.0, beta_glnrho_scaled=0.0
!
  include 'density.h'
  integer :: pushpars2c, pushdiags2c  ! should be procedure pointer (F2003)

  interface calc_pencils_density
    module procedure calc_pencils_density_pnc
    module procedure calc_pencils_density_std
  endinterface calc_pencils_density
!
  contains
!***********************************************************************
    subroutine register_density
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      use EquationOfState, only: select_eos_variable
      use DensityMethods, only: initialize_density_methods
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Tell the equation of state that we're here and we don't have a
!  variable => isochoric (constant density).
!
      call select_eos_variable('lnrho',-1)
      call initialize_density_methods
!
      call keep_compiler_quiet(f)
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
    subroutine density_after_boundary(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
  endsubroutine density_after_boundary
!***********************************************************************
    subroutine pencil_criteria_density
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
    subroutine calc_pencils_density_std(f,p)
!
! Envelope adjusting calc_pencils_density_pnc to the standard use with
! lpenc_loc=lpencil
!
! 21-sep-13/MR: coded
!
      real, dimension (mx,my,mz,mfarray),intent(IN) :: f
      type (pencil_case),                intent(OUT):: p
!
      call calc_pencils_density_pnc(f,p,lpencil)
!
      endsubroutine calc_pencils_density_std
!***********************************************************************
    subroutine calc_pencils_density_pnc(f,p,lpenc_loc)
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
      logical, dimension(:) :: lpenc_loc
!
      intent(in) :: f, lpenc_loc
      intent(inout) :: p
! rho
      if (lpenc_loc(i_rho)) p%rho=rho0
! lnrho
      if (lpenc_loc(i_lnrho)) p%lnrho=lnrho0
! rho1
      if (lpenc_loc(i_rho1)) p%rho1=1/rho0
! glnrho
      if (lpenc_loc(i_glnrho)) p%glnrho=0.0
! grho
      if (lpenc_loc(i_grho)) p%grho=0.0
! del6lnrho
      if (lpenc_loc(i_del6lnrho)) p%del6lnrho=0.0
! hlnrho
      if (lpenc_loc(i_hlnrho)) p%hlnrho=0.0
! sglnrho
      if (lpenc_loc(i_sglnrho)) p%sglnrho=0.0
! uglnrho
      if (lpenc_loc(i_uglnrho)) p%uglnrho=0.0
! ugrho
      if (lpenc_loc(i_ugrho)) p%ugrho=0.0
! uij5glnrho
      if (lpenc_loc(i_uij5glnrho)) p%uij5glnrho=0.0
! ekin
      if (lpenc_loc(i_ekin)) p%ekin=0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_density_pnc
!***********************************************************************
    subroutine density_before_boundary(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine density_before_boundary
!***********************************************************************
    subroutine dlnrho_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,df,p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine calc_diagnostics_density(f,p)

      real, dimension(mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      call keep_compiler_quiet(p)
      call keep_compiler_quiet(f)

    endsubroutine calc_diagnostics_density
!***********************************************************************
    subroutine density_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine density_after_timestep
!***********************************************************************
    subroutine split_update_density(f)
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine split_update_density
!***********************************************************************
    subroutine impose_density_floor(f)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_density_floor
!***********************************************************************
    subroutine read_density_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
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
    subroutine read_density_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
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
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
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
    subroutine get_slices_pressure(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_pressure
!***********************************************************************
    subroutine get_init_average_density(f,init_average_density)
!
!  10-dec-09/piyali: added to pass initial average density
!
    real, dimension (mx,my,mz,mfarray):: f
    real:: init_average_density
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(init_average_density)
!
    endsubroutine get_init_average_density
!***********************************************************************
    subroutine anelastic_after_mn(f, p, df, mass_per_proc)
!
!  14-dec-09/dintrans: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(1) :: mass_per_proc
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(mass_per_proc)
!
    endsubroutine anelastic_after_mn
!***********************************************************************
    subroutine dynamical_diffusion(uc)
!
!  dummy routine
!
      real, intent(in) :: uc
!
      call keep_compiler_quiet(uc)
!
    endsubroutine dynamical_diffusion
!***********************************************************************
    subroutine boussinesq(f)
!
!  23-mar-2012/dintrans: coded
!  dummy routine for the Boussinesq approximation
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine boussinesq
!***********************************************************************
    function mean_density(f)
!
!  Return mean density as rho0 from eos
!
!  1-mar-15/MR: derived from mean_density in density
!
      use EquationOfState, only: rho0
!
      real :: mean_density
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      mean_density=rho0
!
      call keep_compiler_quiet(f)
!
    endfunction mean_density
!***********************************************************************
    subroutine update_char_vel_density(f)
!
!  Updates characteristic velocity for slope-limited diffusion.
!  Most likely not yet a good method
!
!  21-oct-15/MR: coded
!
      use EquationOfState, only: rho0
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      if (lslope_limit_diff) f(2:mx-2,2:my-2,2:mz-2,iFF_char_c) &
                            =f(2:mx-2,2:my-2,2:mz-2,iFF_char_c) + rho0**2
!
    endsubroutine update_char_vel_density
!***********************************************************************
    subroutine impose_density_ceiling(f)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f

      call keep_compiler_quiet(f)

    endsubroutine impose_density_ceiling
!***********************************************************************
endmodule Density
