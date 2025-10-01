! $Id: particles_temperature.f90 21950 2014-07-08 08:53:00Z michiel.lambrechts $
!
!  This module takes care of everything related to the mass of the particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
! CPARAM logical, parameter :: lparticles_mass=.false.
!
!! PENCILS PROVIDED TTp
!
!***************************************************************
module Particles_mass
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_mpicomm
  use Particles_sub
!
  implicit none
!
  include 'particles_mass.h'
!
  contains
!***********************************************************************
    subroutine register_particles_mass()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  23-sep-14/Nils: dummy
!
    endsubroutine register_particles_mass
!***********************************************************************
    subroutine initialize_particles_mass(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  23-sep-14/Nils: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_mass
!***********************************************************************
    subroutine init_particles_mass(f,fp)
!
!  Initial particle mass.
!
!  23-sep-14/Nils: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine init_particles_mass
!***********************************************************************
    subroutine pencil_criteria_par_mass()
!
!  All pencils that the Particles_mass module depends on are specified
!  here.
!
!  23-sep-14/Nils: adapted
!
    endsubroutine pencil_criteria_par_mass
!***********************************************************************
    subroutine dpmass_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle temperature.
!
!  23-sep-14/Nils: coded
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
    endsubroutine dpmass_dt
!***********************************************************************
    subroutine dpmass_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle temperature.
!
!  23-sep-14/Nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(fp,dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpmass_dt_pencil
!***********************************************************************
    subroutine read_particles_mass_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_mass_init_pars
!***********************************************************************
    subroutine write_particles_mass_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_mass_init_pars
!***********************************************************************
    subroutine read_particles_mass_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_mass_run_pars
!***********************************************************************
    subroutine write_particles_mass_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_mass_run_pars
!***********************************************************************
    subroutine rprint_particles_mass(lreset,lwrite)
!
!  Read and register print parameters relevant for particle mass.
!
!  23-sep-14/Nils: adapted
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_mass
!***********************************************************************
endmodule Particles_mass
