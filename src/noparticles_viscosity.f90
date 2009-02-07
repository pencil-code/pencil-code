! $Id$
!
!  This dummy module takes care of the viscosity of inertial particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_viscosity = .false.
!
!***************************************************************
module Particles_viscosity
!
  use Cparam
  use Cdata
  use Messages
  use Particles_cdata
!
  implicit none
!
  include 'particles_viscosity.h'
!
  contains
!***********************************************************************
    subroutine register_particles_viscosity()
!
!  07-oct-08/anders: dummy
!
    endsubroutine register_particles_viscosity
!***********************************************************************
    subroutine initialize_particles_viscosity(lstarting)
!
!  07-oct-08/anders: dummy
!
      logical, intent(in) :: lstarting
!
      if (NO_WARN) print*, lstarting
!
    endsubroutine initialize_particles_viscosity
!***********************************************************************
    subroutine calc_particles_viscosity(f,fp,ineargrid)
!
!  07-oct-08/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f 
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      if (NO_WARN) print*, f, fp, ineargrid
!
    endsubroutine calc_particles_viscosity
!***********************************************************************
    subroutine dvvp_dt_viscosity_pencil(f,df,fp,dfp,ineargrid)
!
!  07-oct-08/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: df, f, fp, dfp, ineargrid
!
      if (NO_WARN) print*, f, df, fp, dfp, ineargrid
!
    endsubroutine dvvp_dt_viscosity_pencil
!***********************************************************************
    subroutine read_particles_visc_init_pars(unit,iostat)
!
!  07-oct-08/anders: dummy
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*, unit
!
    endsubroutine read_particles_visc_init_pars
!***********************************************************************
    subroutine write_particles_visc_init_pars(unit)
!
!  07-oct-08/anders: dummy
!
      integer, intent(in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_visc_init_pars
!***********************************************************************
    subroutine read_particles_visc_run_pars(unit,iostat)
!
!  07-oct-08/anders: dummy
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_particles_visc_run_pars
!***********************************************************************
    subroutine write_particles_visc_run_pars(unit)
!
!  07-oct-08/anders: dummy
!
      integer, intent(in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_visc_run_pars
!***********************************************************************
    subroutine calc_viscosity(f)
!
!  07-oct-08/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine calc_viscosity
!*******************************************************************
    subroutine rprint_particles_viscosity(lreset,lwrite)
!
!  07-oct-08/anders: dummy
!
      logical :: lreset
      logical, optional :: lwrite
!
      if (NO_WARN) print*, lreset, lwrite
!
    endsubroutine rprint_particles_viscosity
!***********************************************************************
endmodule Particles_viscosity
