! $Id$
!
!  Dummy module for energy equation.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!! CPARAM logical, parameter :: lentropy = .false.
!! CPARAM logical, parameter :: ltemperature = .false.
!! CPARAM logical, parameter :: lthermal_energy = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED Ma2; fpres(3); tcond
!
!***************************************************************
module Energy
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
  integer :: pushpars2c  ! should be procedure pointer (F2003)
!
  include 'energy.h'
!
  contains
!***********************************************************************
    subroutine register_energy
!
!  Identify version number.
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_energy
!***********************************************************************
    subroutine initialize_energy(f)
!
!  Called after reading parameters, but before the time loop.
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_energy
!***********************************************************************
    subroutine read_energy_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_energy_init_pars
!***********************************************************************
    subroutine write_energy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_energy_init_pars
!***********************************************************************
    subroutine read_energy_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_energy_run_pars
!***********************************************************************
    subroutine write_energy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_energy_run_pars
!***********************************************************************
    subroutine init_ee(f)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_ee
!***********************************************************************
    subroutine pencil_criteria_energy
!
!  All pencils that the Entropy module depends on are specified here.
!
    endsubroutine pencil_criteria_energy
!***********************************************************************
    subroutine pencil_interdep_energy(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
      logical, dimension(npencils), intent(in) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_energy
!***********************************************************************
    subroutine calc_pencils_energy(f,p)
!
!  Calculate Energy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_energy
!**********************************************************************
    subroutine dee_dt(f,df,p)
!
!  Calculate right hand side of energy equation.
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(in) :: df
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dee_dt
!***********************************************************************
    subroutine energy_after_boundary(f)
!
!  dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine energy_after_boundary
!***********************************************************************
    subroutine read_energy_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_energy_init_pars
!***********************************************************************
    subroutine write_energy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_energy_init_pars
!***********************************************************************
    subroutine read_energy_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_energy_run_pars
!***********************************************************************
    subroutine write_energy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_energy_run_pars
!***********************************************************************
    subroutine get_slices_energy(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_energy
!***********************************************************************
    subroutine fill_farray_pressure(f)
!
!  Fill f array with the pressure, to be able to calculate pressure gradient
!  directly from the pressure.
!
!  18-feb-10/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine fill_farray_pressure
!***********************************************************************
    subroutine impose_energy_floor(f)
!
!  Dummy subroutine; may not be necessary for lnTT
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_energy_floor
!***********************************************************************
    subroutine split_update_energy(f)
!
!  Dummy subroutine
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine
!***********************************************************************
    subroutine expand_shands_energy
!
!  Presently dummy, for possible use
!
    endsubroutine expand_shands_energy
!***********************************************************************
    subroutine rprint_energy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_energy
!***********************************************************************
endmodule Energy
