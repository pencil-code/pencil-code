! $Id: viscosity.f90 9840 2008-09-05 07:29:37Z ajohan $
!
!  This modules takes care of viscosity of inertial particles.
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
!  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_viscosity.h'
!
  contains
!***********************************************************************
    subroutine register_particles_viscosity()
!
!  07-oct-08/anders: coded
!
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
!  Identify version number.
!
      if (lroot) call cvs_id( &
           "$Id: viscosity.f90 9840 2008-09-05 07:29:37Z ajohan $")
!
    endsubroutine register_particles_viscosity
!***********************************************************************
    subroutine initialize_particles_viscosity(lstarting)
!
!  07-oct-08/anders: coded
!
      logical, intent(in) :: lstarting
!
!      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_viscosity
!***********************************************************************
    subroutine calc_particles_viscosity(f,fp,ineargrid)
!
!
!
      use Particles_sub, only: map_vvp_grid
      use Sub, only: del2v
!
      real, dimension (mx,my,mz,mfarray) :: f 
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
!      call keep_compiler_quiet(f,fp,ineargrid)
!
    endsubroutine calc_particles_viscosity
!***********************************************************************
    subroutine calc_particles_viscous_force(df,p)
!
!  calculate viscous force term for right hand side of  equation
!
!  20-nov-02/tony: coded
!   9-jul-04/nils: added Smagorinsky viscosity

      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: nu_smag
      real, dimension (nx,3) :: nuD2uxb
      type (pencil_case) :: p
!
      intent (in) :: p
      intent (inout) :: df 
!
!      call keep_compiler_quiet(df,p)
!
    endsubroutine calc_particles_viscous_force
!***********************************************************************
    subroutine read_particles_visc_init_pars(unit,iostat)
!
!
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit
!
    endsubroutine read_particles_visc_init_pars
!***********************************************************************
    subroutine write_particles_visc_init_pars(unit)
!
!
!
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit
!
    endsubroutine write_particles_visc_init_pars
!***********************************************************************
    subroutine read_particles_visc_run_pars(unit,iostat)
!
!
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
    endsubroutine read_particles_visc_run_pars
!***********************************************************************
    subroutine write_particles_visc_run_pars(unit)
!
!
!
      integer, intent(in) :: unit
!
    endsubroutine write_particles_visc_run_pars
!***********************************************************************
    subroutine calc_viscosity(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
    endsubroutine calc_viscosity
!*******************************************************************
    subroutine rprint_particles_viscosity(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  07-oct-08/anders: adapted
!
      use Sub
!
      logical :: lreset
      logical, optional :: lwrite
!
      if (NO_WARN) print*, lreset, lwrite
!
    endsubroutine rprint_particles_viscosity
!***********************************************************************
endmodule Particles_viscosity
