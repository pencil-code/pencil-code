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
  use Cparam
  use Cdata
  use Messages
  use Particles_cdata
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_collisions.h'
!
  contains
!***********************************************************************
    subroutine initialize_particles_collisions(f,lstarting)
!
!  23-mar-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_collisions
!***********************************************************************
    subroutine particles_collisions_timestep(fp,ineargrid)
!
!  30-nov-10/anders: dummy
!
      real, dimension (mpar_loc,mpvar) :: fp
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
      real, dimension (mpar_loc,mpvar) :: fp
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
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_collisions_blocks
!***********************************************************************
    subroutine read_particles_coll_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
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
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_particles_collisions
!***********************************************************************
endmodule Particles_collisions
