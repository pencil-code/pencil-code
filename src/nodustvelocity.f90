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
  public :: unit_md, dust_chemistry, mumon, mmon, mi, md
  character (len=labellen) :: dust_geometry='sphere', dust_chemistry='nothing'
  real, dimension(ndustspec,ndustspec) :: scolld
  real, dimension(nx,ndustspec) :: tausd1
  real, dimension(ndustspec) :: md=1.0, mdplus=0.0, mdminus=0.0, surfd=0.0
  real, dimension(ndustspec) :: mi=0.0, ad=1.0, tausd=0.0
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
    subroutine read_dustvelocity_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
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
    subroutine read_dustvelocity_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
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
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lwr) then
       if (lmdvar .and. lmice ) then
        write(3,*) 'iuud=',iuud
        write(3,*) 'iudx=',iudx
        write(3,*) 'iudy=',iudy
        write(3,*) 'iudz=',iudz
       endif
      endif
!
      call keep_compiler_quiet(lreset)
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
endmodule Dustvelocity
