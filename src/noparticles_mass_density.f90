! $Id: particles_number.f90 14421 2010-07-23 23:55:02Z Bourdin.KIS $
!
!  This module takes care of everything related to the mass density
!  represented by each (super)particle.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_mass_density=.false.
!
!***************************************************************
module Particles_mass_density
!
  use Cdata
  use Messages
  use Particles_cdata
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_mass_density.h'
!
  contains
!***********************************************************************
    subroutine register_particles_mass_density()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  22-nov-10/anders+michiel: dummy
!
    endsubroutine register_particles_mass_density
!***********************************************************************
    subroutine initialize_particles_mass_density(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  22-nov-10/anders+michiel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_mass_density
!***********************************************************************
    subroutine init_particles_mass_density(f,fp)
!
!  Initial particle mass density.
!
!  22-nov-10/anders+michiel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine init_particles_mass_density
!***********************************************************************
    subroutine pencil_criteria_par_mass_density()
!
!  All pencils that the Particles_mass_density module depends on are specified
!  here.
!
!  22-nov-10/anders+michiel: dummy
!
    endsubroutine pencil_criteria_par_mass_density
!***********************************************************************
    subroutine drhopswarm_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle mass density.
!
!  22-nov-10/anders+michiel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, fp
      intent (out) :: dfp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine drhopswarm_dt_pencil
!***********************************************************************
    subroutine drhopswarm_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of internal particle number.
!
!  22-nov-10/anders+michiel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine drhopswarm_dt
!***********************************************************************
    subroutine read_particles_dens_init_pars(unit,iostat)
!
!  22-nov-10/anders+michiel: dummy
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_dens_init_pars
!***********************************************************************
    subroutine write_particles_dens_init_pars(unit)
!
!  22-nov-10/anders+michiel: dummy
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_dens_init_pars
!***********************************************************************
    subroutine read_particles_dens_run_pars(unit,iostat)
!
!  22-nov-10/anders+michiel: dummy
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_dens_run_pars
!***********************************************************************
    subroutine write_particles_dens_run_pars(unit)
!
!  22-nov-10/anders+michiel: dummy
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_dens_run_pars
!***********************************************************************
    subroutine rprint_particles_mass_density(lreset,lwrite)
!
!  Read and register print parameters relevant for particle mass density.
!
!  22-nov-10/anders+michiel: dummy
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_mass_density
!***********************************************************************
endmodule Particles_mass_density
