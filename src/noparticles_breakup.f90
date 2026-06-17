! $Id$
!
!  Dummy module for particle breakup.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! CPARAM logical, parameter :: lparticles_breakup = .false.
!
! MPVAR CONTRIBUTION 0
! MPAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_breakup
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Particles_cdata
!
  implicit none
!
  include 'particles_breakup.h'
!
contains
!***********************************************************************
    subroutine register_particles_breakup()
    endsubroutine register_particles_breakup
!***********************************************************************
    subroutine initialize_particles_breakup(f)
      real, dimension(mx,my,mz,mfarray) :: f
      call keep_compiler_quiet(f)
    endsubroutine initialize_particles_breakup
!***********************************************************************
    subroutine read_particles_breakup_init_pars(iomsg)
      character(LEN=iomsglen), intent(out) :: iomsg
      iomsg=""
    endsubroutine read_particles_breakup_init_pars
!***********************************************************************
    subroutine write_particles_breakup_init_pars(unit)
      integer, intent(in) :: unit
      call keep_compiler_quiet(unit)
    endsubroutine write_particles_breakup_init_pars
!***********************************************************************
    subroutine read_particles_breakup_run_pars(iomsg)
      character(LEN=iomsglen), intent(out) :: iomsg
      iomsg=""
    endsubroutine read_particles_breakup_run_pars
!***********************************************************************
    subroutine write_particles_breakup_run_pars(unit)
      integer, intent(in) :: unit
      call keep_compiler_quiet(unit)
    endsubroutine write_particles_breakup_run_pars
!***********************************************************************
    subroutine rprint_particles_breakup(lreset,lwrite)
      logical, intent(in) :: lreset
      logical, optional, intent(in) :: lwrite
      call keep_compiler_quiet(lreset,lwrite)
    endsubroutine rprint_particles_breakup
!***********************************************************************
    subroutine dbreakup_dt(f,df,fp,dfp,ineargrid)
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mpar_loc,mpvar), intent(inout) :: dfp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
    endsubroutine dbreakup_dt
!***********************************************************************
    subroutine particles_breakup_pencils(f,fp,ineargrid)
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
      integer, dimension(mpar_loc,3), intent(inout) :: ineargrid
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
    endsubroutine particles_breakup_pencils
!***********************************************************************
endmodule Particles_breakup
