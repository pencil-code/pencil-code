! $Id$
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
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
!
  implicit none
!
  include 'particles_diagnos_dv.h'
!
  contains
!***********************************************************************
    subroutine initialize_particles_diagnos_dv(f)
!
!  Dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_diagnos_dv
!***********************************************************************
    subroutine read_pars_diagnos_dv_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
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
!***********************************************************************
    subroutine rprint_particles_diagnos_dv(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_particles_diagnos_dv
!***********************************************************************
    subroutine collisions(fp)
!
      real, dimension (mpar_loc,mparray) :: fp
!
      call keep_compiler_quiet(fp)
!
    endsubroutine collisions
!***********************************************************************
    subroutine repeated_init(fp,init_repeat)
!
!
      real, dimension (mpar_loc,mparray) :: fp
      integer :: init_repeat
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(init_repeat)
!
    endsubroutine repeated_init
!***********************************************************************
endmodule Particles_diagnos_dv
