!
!  This modules takes care of condensation / evaporation or deposition /
!  sublimation of superparticles
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_condensation= .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_condensation
!
  use Cdata
  use Cparam
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_condensation.h'
!
  contains
!***********************************************************************
    subroutine initialize_particles_cond(f)
!
!  18-jun-17/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Fatal error if Particle_radius module not used.
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_cond
!***********************************************************************
    subroutine particles_condensation_timestep(fp,ineargrid)
!
!  17-jun-17/anders: coded
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_condensation_timestep
!***********************************************************************
    subroutine particles_condensation_pencils(fp,ineargrid)
!
!  17-jun-17/anders: coded
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: ineargrid
      intent (inout) :: fp
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_condensation_pencils
!***********************************************************************
    subroutine particles_condensation_blocks(fp,ineargrid)
!
!  18-jun-17/anders: coded
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_condensation_blocks
!***********************************************************************
    subroutine read_particles_cond_init_pars(iostat)
!
!  18-jun-17/anders: coded
!
      integer, intent(out) :: iostat
!
      call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_cond_init_pars
!***********************************************************************
    subroutine write_particles_cond_init_pars(unit)
!
!  18-jun-17/anders: coded
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_cond_init_pars
!***********************************************************************
    subroutine read_particles_cond_run_pars(iostat)
!
!  18-jun-17/anders: coded
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_cond_run_pars
!***********************************************************************
    subroutine write_particles_cond_run_pars(unit)
!
!  18-jun-17/anders: coded
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_cond_run_pars
!***********************************************************************
    subroutine rprint_particles_condensation(lreset,lwrite)
!
!  18-jun-17/anders: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_condensation
!***********************************************************************
endmodule Particles_condensation
