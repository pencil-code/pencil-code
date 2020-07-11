! $Id$
!
!  This module takes care of everything related to particle radius.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_radius=.false.
!
!***************************************************************
module Particles_radius
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Particles_cdata
!
  implicit none
!
  include 'particles_radius.h'
!
  contains
!***********************************************************************
    subroutine register_particles_radius()
!
!  Set up indices for access to the fp and dfp arrays
!
!  22-aug-05/anders: dummy
!
    endsubroutine register_particles_radius
!***********************************************************************
    subroutine initialize_particles_radius(f,fp)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  25-nov-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initialize_particles_radius
!***********************************************************************
    subroutine set_particle_radius(f,fp,npar_low,npar_high,init)
!
!  Set radius of new particles.
!
!  18-sep-09/nils: adapted from init_particles_radius
!
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer :: npar_low,npar_high
      logical, optional :: init
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(npar_low)
      call keep_compiler_quiet(npar_high)
      call keep_compiler_quiet(present(init))
!
    endsubroutine set_particle_radius
!***********************************************************************
    subroutine pencil_criteria_par_radius()
!
!  All pencils that the Particles_radius module depends on are specified here.
!
!  21-nov-06/anders: dummy
!
    endsubroutine pencil_criteria_par_radius
!***********************************************************************
    subroutine dap_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle radius.
!
!  21-nov-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(p)
!
    endsubroutine dap_dt_pencil
!***********************************************************************
    subroutine dap_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle radius.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dap_dt
!***********************************************************************
    subroutine read_particles_rad_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_rad_init_pars
!***********************************************************************
    subroutine write_particles_rad_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_rad_init_pars
!***********************************************************************
    subroutine read_particles_rad_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_rad_run_pars
!***********************************************************************
    subroutine write_particles_rad_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_rad_run_pars
!***********************************************************************
    subroutine rprint_particles_radius(lreset,lwrite)
!
!  Read and register print parameters relevant for particles radius.
!
!  22-aug-05/anders: dummy
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_radius
!***********************************************************************
    subroutine get_stbin(iStbin,fp,ip)
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, intent(out) :: iStbin
      integer, intent(in) :: ip
!
      iStbin = -1
      call keep_compiler_quiet(ip)
      call keep_compiler_quiet(fp)
!
    endsubroutine get_stbin
!***********************************************************************
    subroutine get_mass_from_radius(mpi,fp,ip)
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, intent(in) :: ip
      real, intent(out) :: mpi
      real :: api
!
      mpi = -1.0
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ip)
!
    endsubroutine get_mass_from_radius
!***********************************************************************
    subroutine get_maxrad(apmax)
      real :: apmax

    endsubroutine get_maxrad
!***********************************************************************
endmodule Particles_radius
