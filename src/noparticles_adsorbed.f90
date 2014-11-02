! $Id: noparticles_adsorbed.f90 19193 2012-06-30 12:55:46Z wdobler $
!
!  This module takes care of everything related to particle
!  surface fractions.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_adsorbed=.false.
!
!***************************************************************
module Particles_adsorbed
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Particles_cdata
!
  implicit none
!
  include 'particles_adsorbed.h'
!
  contains
!***********************************************************************
    subroutine register_particles_ads()
!
!  Set up indices for access to the fp and dfp arrays
!
!  28-aug-14/jonas+nils: coded
!
    endsubroutine register_particles_ads
!***********************************************************************
    subroutine initialize_particles_ads(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  28-aug-14/jonas+nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_ads
!***********************************************************************
    subroutine init_particles_ads(f,fp)
!
!  Initial surface fractions of particles.
!
!  28-aug-14/jonas+nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine init_particles_ads
!***********************************************************************
    subroutine pencil_criteria_par_ads()
!
!  All pencils that the Particles_adsorbed
!  module depends on are specified here.
!
!  28-aug-14/jonas+nils: coded
!
    endsubroutine pencil_criteria_par_ads
!***********************************************************************
    subroutine dpads_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle surface fraction
!
!  01-sep-14/jonas: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, df, fp, ineargrid
      intent (inout) :: dfp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpads_dt_pencil
!***********************************************************************
    subroutine dpads_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle surface fractions
!
!  28-aug-14/jonas+nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpads_dt
!***********************************************************************
    subroutine read_particles_ads_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_ads_init_pars
!***********************************************************************
    subroutine write_particles_ads_init_pars(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_ads_init_pars
!***********************************************************************
    subroutine read_particles_ads_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_ads_run_pars
!***********************************************************************
    subroutine write_particles_ads_run_pars(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_ads_run_pars
!***********************************************************************
    subroutine rprint_particles_ads(lreset,lwrite)
!
!  Read and register print parameters relevant for particles
!  surface fractions.
!
!  28-aug-14/jonas+nils: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_particles_ads
!***********************************************************************
    subroutine particles_ads_prepencil_calc(f)
!
!  28-aug-14/jonas+nils: coded
!
      real,dimension(mx,my,mz,mfarray),intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_ads_prepencil_calc
!***********************************************************************
endmodule Particles_adsorbed
