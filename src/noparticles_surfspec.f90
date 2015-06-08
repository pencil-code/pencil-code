! $Id: noparticles_surfspec.f90 19193 2012-06-30 12:55:46Z wdobler $
!
!  This module takes care of everything related to particle
!  surface fractions.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_surfspec=.false.
!
!***************************************************************
module Particles_surfspec

  use Cdata
  use General, only: keep_compiler_quiet
  use Particles_cdata

  implicit none

  include 'particles_surfspec.h'

  contains
! ******************************************************************************
!  09.09.14/jonas : coded

  subroutine register_particles_surfspec()
  endsubroutine register_particles_surfspec
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine register_dep_psurfchem()
  endsubroutine register_dep_psurfchem
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine register_indep_psurfchem()
  endsubroutine register_indep_psurfchem
! ******************************************************************************

  subroutine init_particles_surf(f,fp)
    real, dimension(mx,my,mz,mfarray) :: f
    real, dimension(mpar_loc,mparray) :: fp

    ! 19.09.2014/Jonas:coded
    call keep_compiler_quiet(f)
    call keep_compiler_quiet(fp)
  endsubroutine init_particles_surf
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine read_particles_surf_init_pars(iostat)
!
    integer, intent(out) :: iostat

    iostat = 0
  endsubroutine read_particles_surf_init_pars
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine write_particles_surf_init_pars(unit)
    integer, intent(in) :: unit

    call keep_compiler_quiet(unit)
  endsubroutine write_particles_surf_init_pars
! ******************************************************************************
!  19.09.2014/Jonas:coded
!

  subroutine read_particles_surf_run_pars(iostat)
!
    integer, intent(out) :: iostat

    iostat = 0
  endsubroutine read_particles_surf_run_pars
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine write_particles_surf_run_pars(unit)
    integer, intent(in) :: unit

    call keep_compiler_quiet(unit)
  endsubroutine write_particles_surf_run_pars
! ******************************************************************************
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-sep-14/jonas coded
!  JONAS: needs to be filled with life

  subroutine initialize_particles_surf(f)
    real, dimension(mx,my,mz,mfarray) :: f

    call keep_compiler_quiet(f)
  endsubroutine initialize_particles_surf
! ******************************************************************************
!  Evolution of particle near field composition
!
!  28-aug-14/jonas+nils: coded

  subroutine dpsurf_dt(f,df,fp,dfp,ineargrid)
    real, dimension(mx,my,mz,mfarray) :: f
    real, dimension(mx,my,mz,mvar) :: df
    real, dimension(mpar_loc,mparray) :: fp
    real, dimension(mpar_loc,mpvar) :: dfp
    integer, dimension(mpar_loc,3) :: ineargrid

    call keep_compiler_quiet(f)
    call keep_compiler_quiet(df)
    call keep_compiler_quiet(fp)
    call keep_compiler_quiet(dfp)
    call keep_compiler_quiet(ineargrid)
  endsubroutine dpsurf_dt
! ******************************************************************************
!  Evolution of particle surface fraction
!
!  01-sep-14/jonas: coded

  subroutine dpsurf_dt_pencil(f,df,fp,dfp,p,ineargrid)
    real, dimension(mx,my,mz,mfarray) :: f
    real, dimension(mx,my,mz,mvar) :: df
    real, dimension(mpar_loc,mparray) :: fp
    real, dimension(mpar_loc,mpvar) :: dfp
    type (pencil_case) :: p
    integer, dimension(mpar_loc,3) :: ineargrid

    intent(in) :: f, df, fp, ineargrid
    intent(inout) :: dfp

    call keep_compiler_quiet(f)
    call keep_compiler_quiet(df)
    call keep_compiler_quiet(fp)
    call keep_compiler_quiet(dfp)
    call keep_compiler_quiet(p)
    call keep_compiler_quiet(ineargrid)
  endsubroutine dpsurf_dt_pencil
! ******************************************************************************
!  Read and register print parameters relevant for particles
!  surface fractions.
!
!  28-aug-14/jonas+nils: coded

  subroutine rprint_particles_surf(lreset,lwrite)
    logical :: lreset
    logical, optional :: lwrite

    if (present(lwrite)) call keep_compiler_quiet(lwrite)

    call keep_compiler_quiet(lreset)
  endsubroutine rprint_particles_surf
! ******************************************************************************
!
!  nov-14/jonas: coded

  subroutine calc_psurf_pencils(f,fp,p,ineargrid)
    real, dimension(mpar_loc,mparray) :: fp
    real, dimension(mx,my,mz,mfarray) :: f
    integer, dimension(mpar_loc,3) :: ineargrid
    type (pencil_case) :: p
    call keep_compiler_quiet(f)
    call keep_compiler_quiet(fp)
    call keep_compiler_quiet(p)
    call keep_compiler_quiet(ineargrid)
  endsubroutine calc_psurf_pencils
! ******************************************************************************
!  nov-14/jonas: coded

  subroutine cleanup_surf_pencils()
  endsubroutine cleanup_surf_pencils
! ******************************************************************************
  subroutine particles_surfspec_clean_up()
  endsubroutine particles_surfspec_clean_up
! ******************************************************************************
endmodule Particles_surfspec
