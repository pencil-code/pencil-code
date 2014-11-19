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
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Particles_cdata
!
  implicit none
!
  include 'particles_surfspec.h'
!
  contains
!***********************************************************************
    subroutine register_particles_surfspec()
!
!  09.09.14/jonas : coded
!
    end subroutine register_particles_surfspec
!***********************************************************************
    subroutine register_dep_psurfchem()
!
!  19.09.2014/Jonas:coded
!
    end subroutine register_dep_psurfchem
!***********************************************************************
    subroutine register_indep_psurfchem()
!
!  19.09.2014/Jonas:coded
!
    end subroutine register_indep_psurfchem
!**********************************************************************
    subroutine init_particles_surf(f,fp)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
!
!  19.09.2014/Jonas:coded
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    end subroutine init_particles_surf
!***********************************************************************
    subroutine read_particles_surf_run_pars(unit,iostat)
!
!  19.09.2014/Jonas:coded
!
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    end subroutine read_particles_surf_run_pars
!***********************************************************************
    subroutine read_particles_surf_init_pars(unit,iostat)
!
!  19.09.2014/Jonas:coded
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    end subroutine read_particles_surf_init_pars
!***********************************************************************
    subroutine write_particles_surf_init_pars(unit)
!
!  19.09.2014/Jonas:coded
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    end subroutine write_particles_surf_init_pars
!***********************************************************************
    subroutine write_particles_surf_run_pars(unit)
!
!  19.09.2014/Jonas:coded
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    end subroutine write_particles_surf_run_pars
!***********************************************************************
    subroutine initialize_particles_surf(f,fp)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-sep-14/jonas coded
!  JONAS: needs to be filled with life
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    end subroutine initialize_particles_surf
!**************************************************************
    subroutine dpsurf_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle near field composition
!
!  28-aug-14/jonas+nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpsurf_dt
!**************************************************************
    subroutine dpsurf_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle surface fraction
!
!  01-sep-14/jonas: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
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
    endsubroutine dpsurf_dt_pencil
!***********************************************************************
    subroutine rprint_particles_surf(lreset,lwrite)
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
    endsubroutine rprint_particles_surf
!***********************************************************************
end module Particles_surfspec
