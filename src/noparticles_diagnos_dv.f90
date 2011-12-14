! $Id: particles_collisions.f90 17701 2011-10-07 15:14:06Z michiel.lambrechts $
!
!  This modules takes care of instantaneous collisions between
!  superparticles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_diagnos_dv = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_diagnos_dv
!
  use Cparam
  use Cdata
  use Messages
  use Particles_cdata
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_diagnos_dv.h'
!
  contains
!***********************************************************************
    subroutine initialize_particles_diagnos_dv(f,lstarting)
!
!  Dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_diagnos_dv
!***********************************************************************
    subroutine read_pars_diagnos_dv_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_pars_diagnos_dv_run_pars
!***********************************************************************
    subroutine write_pars_diagnos_dv_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_pars_diagnos_dv_run_pars
!*******************************************************************
    subroutine rprint_particles_diagnos_dv(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_particles_diagnos_dv
!***********************************************************************
    subroutine collisions(fp)
!
      real, dimension (mpar_loc,mpvar) :: fp
!
      call keep_compiler_quiet(fp)
!
    endsubroutine collisions
!***********************************************************************
    subroutine repeated_init(fp,init_repeat)
!
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: init_repeat
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(init_repeat)
!
    endsubroutine repeated_init
!***********************************************************************
endmodule Particles_diagnos_dv
