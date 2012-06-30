! $Id$
!
!  This modules takes care of instantaneous collisions between
!  superparticles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_diagnos_state = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_diagnos_state
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
!
  implicit none
!
  include 'particles_diagnos_state.h'
!
  contains
!***********************************************************************
    subroutine register_pars_diagnos_state()
!
!  Dummy
!
    endsubroutine register_pars_diagnos_state
!***********************************************************************
    subroutine initialize_pars_diagnos_state(f,lstarting)
!
!  Dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_pars_diagnos_state
!***********************************************************************
    subroutine init_particles_diagnos_state(fp)
!
!  Dummy
!
      real, dimension (mpar_loc,mpvar) :: fp
!
      call keep_compiler_quiet(fp)
!
    endsubroutine init_particles_diagnos_state
!***********************************************************************
    subroutine insert_particles_diagnos_state(fp, npar_loc_old)
!
!  Dummy
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: npar_loc_old
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(npar_loc_old)
!
    endsubroutine insert_particles_diagnos_state
!***********************************************************************
    subroutine read_pars_diag_state_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_pars_diag_state_run_pars
!***********************************************************************
    subroutine write_pars_diag_state_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_pars_diag_state_run_pars
!*******************************************************************
    subroutine rprint_particles_diagnos_state(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_particles_diagnos_state
!***********************************************************************
    subroutine persistence_check(fp, k, uup)
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: k
      real, dimension(3) :: uup
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(uup)
!
    endsubroutine persistence_check
!***********************************************************************
endmodule Particles_diagnos_state
