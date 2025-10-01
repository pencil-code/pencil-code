! $Id$
!
!  This modules takes care of instantaneous collisions between
!  superparticles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_collisions = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_collisions
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
!
  implicit none
!
  include 'particles_collisions.h'
!
  contains
!***********************************************************************
    subroutine initialize_particles_collisions(f)
!
!  23-mar-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_collisions
!***********************************************************************
    subroutine particles_collisions_timestep(fp,ineargrid)
!
!  30-nov-10/anders: dummy
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_collisions_timestep
!***********************************************************************
    subroutine particles_collisions_pencils(fp,ineargrid)
!
!  23-mar-09/anders: dummy
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_collisions_pencils
!***********************************************************************
    subroutine particles_collisions_blocks(fp,ineargrid)
!
!  17-nov-09/anders: dummy
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_collisions_blocks
!***********************************************************************
    subroutine read_particles_coll_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_coll_run_pars
!***********************************************************************
    subroutine write_particles_coll_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_coll_run_pars
!***********************************************************************
    subroutine rprint_particles_collisions(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_particles_collisions
!***********************************************************************
endmodule Particles_collisions
