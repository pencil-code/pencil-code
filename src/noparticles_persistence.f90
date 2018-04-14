! $Id$
!
! This module calculates the probability distribution function
! the particles that have moved a certain distance away from their
! initial position. 
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_persistence=.false.
! MPAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_persistence
!
  use Cdata
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
!
  implicit none
!
  include 'particles_persistence.h'
!
contains
!***********************************************************************
    subroutine register_particles_persistence
!
!  Dummy.
!
    endsubroutine register_particles_persistence
!***********************************************************************
    subroutine initialize_particles_persist(fp)
!
!  Dummy.
!
      use General, only: keep_compiler_quiet

      real, dimension (mpar_loc,mparray), intent (in) :: fp
!
      call keep_compiler_quiet(fp)

    endsubroutine initialize_particles_persist
!***********************************************************************
    subroutine init_particles_persistence(fp)
!
      use General, only: keep_compiler_quiet
!
      real, dimension (mpar_loc,mparray), intent (out) :: fp
!
      call keep_compiler_quiet(fp)

    endsubroutine init_particles_persistence
!***********************************************************************
    subroutine dpersist_dt(f,df,fp,dfp,ineargrid)
!
      use General, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,my,mz,mvar), intent (inout) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar), intent (inout) :: dfp
      integer, dimension (mpar_loc,3), intent (in) :: ineargrid
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(fp,dfp)
      call keep_compiler_quiet(ineargrid)
    endsubroutine dpersist_dt
!***********************************************************************
    subroutine read_ppersist_init_pars(iostat)
!
      use General, only: keep_compiler_quiet
!
      integer, intent(out) :: iostat
!
      call keep_compiler_quiet(iostat)
    endsubroutine read_ppersist_init_pars
!***********************************************************************
    subroutine write_ppersist_init_pars(unit)
!
      use General, only: keep_compiler_quiet
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
    endsubroutine write_ppersist_init_pars
!***********************************************************************
    subroutine read_ppersist_run_pars(iostat)
!
      use General, only: keep_compiler_quiet
!
      integer, intent(out) :: iostat
!
      call keep_compiler_quiet(iostat)
!
    endsubroutine read_ppersist_run_pars
!***********************************************************************
    subroutine write_ppersist_run_pars(unit)
!
      use General, only: keep_compiler_quiet
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
    endsubroutine write_ppersist_run_pars
!***********************************************************************
    subroutine rprint_particles_persist(lreset,lwrite)
!
      use General, only: keep_compiler_quiet
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_particles_persist
!***********************************************************************
endmodule Particles_persistence 
