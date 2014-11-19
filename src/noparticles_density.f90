! $Id: noparticles_density.f90 19193 2012-06-30 12:55:46Z wdobler $
!
!  This module takes care of everything related to the density represented by
!  each (super)particle.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_density=.false.
!
!***************************************************************
module Particles_density
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
!
  implicit none
!
  include 'particles_density.h'
!
  contains
!***********************************************************************
    subroutine register_particles_density()
!
!  22-nov-10/anders+michiel: dummy
!
    endsubroutine register_particles_density
!***********************************************************************
    subroutine initialize_particles_density(f)
!
!  22-nov-10/anders+michiel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_density
!***********************************************************************
    subroutine init_particles_density(f,fp)
!
!  22-nov-10/anders+michiel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine init_particles_density
!***********************************************************************
    subroutine pencil_criteria_par_density()
!
!  22-nov-10/anders+michiel: dummy
!
    endsubroutine pencil_criteria_par_density
!***********************************************************************
    subroutine drhopswarm_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  22-nov-10/anders+michiel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      	real, dimension (mpar_loc,mparray) :: fp
	real, dimension (mpar_loc,mpvar) :: dfp
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
!  22-nov-10/anders+michiel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      	real, dimension (mpar_loc,mparray) :: fp
	real, dimension (mpar_loc,mpvar) :: dfp
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
    subroutine rprint_particles_density(lreset,lwrite)
!
!  22-nov-10/anders+michiel: dummy
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_density
!***********************************************************************
endmodule Particles_density
