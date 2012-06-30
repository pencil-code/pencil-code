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
  use Cdata
  use General, only: keep_compiler_quiet
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
    subroutine initialize_particles_viscosity(f,lstarting)
!
!  07-oct-08/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
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
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_visc_init_pars
!***********************************************************************
    subroutine write_particles_visc_init_pars(unit)
!
!  07-oct-08/anders: dummy
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
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
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_visc_run_pars
!***********************************************************************
    subroutine write_particles_visc_run_pars(unit)
!
!  07-oct-08/anders: dummy
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_visc_run_pars
!***********************************************************************
    subroutine calc_viscosity(f)
!
!  07-oct-08/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_viscosity
!***********************************************************************
    subroutine rprint_particles_viscosity(lreset,lwrite)
!
!  07-oct-08/anders: dummy
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_particles_viscosity
!***********************************************************************
endmodule Particles_viscosity
