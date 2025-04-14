! $Id$
!
!  This module takes care of everything related to dust velocity.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldustvelocity = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED uud(3,ndustspec); divud(ndustspec); sdij(3,3,ndustspec)
!
!***************************************************************
module Dustvelocity
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'dustvelocity.h'
!
  public :: dust_geometry, dimd1, rhods, surfd, mdplus, mdminus
  public :: ad, scolld, ustcst, tausd1, tausd
  public :: unit_md, dust_chemistry, mumon, mmon, md
  character (len=labellen) :: dust_geometry='sphere', dust_chemistry='nothing'
  real, dimension(ndustspec,ndustspec) :: scolld
  real, dimension(nx,ndustspec) :: tausd1
  real, dimension(ndustspec) :: md=1.0, mdplus=0.0, mdminus=0.0, surfd=0.0
  real, dimension(ndustspec) :: ad=1.0, tausd=0.0
  character (len=labellen) :: dust_binning='log_mass'
  real :: dimd1=0.0, rhods=0.0, ustcst=0.0, unit_md=0.0
  real :: mumon=0.0, mmon=0.0
  !! SHOULDN'T REALLY BE PUBLIC!!
  real :: nd0=1.,rhod0=1.
  logical :: ldustcoagulation=.false., ldustcondensation=.false.
!
  contains
!***********************************************************************
    subroutine register_dustvelocity()
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_dustvelocity
!***********************************************************************
    subroutine initialize_dustvelocity(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      ud_spec=.false.
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_dustvelocity
!***********************************************************************
    subroutine copy_bcs_dust
!
    endsubroutine copy_bcs_dust
!***********************************************************************
    subroutine init_uud(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_uud
!***********************************************************************
    subroutine pencil_criteria_dustvelocity()
!
    endsubroutine pencil_criteria_dustvelocity
!***********************************************************************
    subroutine pencil_interdep_dustvelocity(lpencil_in)
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_dustvelocity
!***********************************************************************
    subroutine calc_pencils_dustvelocity(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_dustvelocity
!***********************************************************************
    subroutine duud_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine duud_dt
!***********************************************************************
    subroutine calc_diagnostics_dustvelocity(p)
!
      type (pencil_case) :: p
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_diagnostics_dustvelocity
!***********************************************************************
    subroutine read_dustvelocity_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_dustvelocity_init_pars
!***********************************************************************
    subroutine write_dustvelocity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_dustvelocity_init_pars
!***********************************************************************
    subroutine read_dustvelocity_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_dustvelocity_run_pars
!***********************************************************************
    subroutine write_dustvelocity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_dustvelocity_run_pars
!***********************************************************************
    subroutine rprint_dustvelocity(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_dustvelocity
!***********************************************************************
    subroutine get_slices_dustvelocity(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_dustvelocity
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    endsubroutine pushpars2c
!***********************************************************************
endmodule Dustvelocity
