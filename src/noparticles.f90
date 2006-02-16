! $Id: noparticles.f90,v 1.8 2006-02-16 12:51:45 ajohan Exp $
!
!  This module takes care of everything related to no particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles=.false.
! CPARAM logical, parameter :: lparticles_planet=.false.
!
!***************************************************************
module Particles

  use Cdata
  use Particles_cdata

  implicit none

  include 'particles.h'

  contains

!***********************************************************************
    subroutine register_particles()
!
!  Set up indices for access to the fp and dfp arrays
!
!  22-aug-05/anders: dummy
!
    endsubroutine register_particles
!***********************************************************************
    subroutine initialize_particles(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  22-aug-05/anders: dummy
!
      logical :: lstarting
!
      if (NO_WARN) print*, lstarting
!
    endsubroutine initialize_particles
!***********************************************************************
    subroutine init_particles(f,fp)
!
!  Initial positions and velocities of particles.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp
!
      if (NO_WARN) print*, f, fp
!
    endsubroutine init_particles
!***********************************************************************
    subroutine pencil_interdep_particles(lpencil_in)
!   
!  Interdependency among pencils provided by the Particles module
!  is specified here.
!         
!  16-feb-06/anders: dummy
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in
!
    endsubroutine pencil_interdep_particles
!***********************************************************************
    subroutine calc_pencils_particles(f,p)
!   
!  Calculate particle pencils.
!
!  16-feb-06/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!
      if (NO_WARN) print*, f, p
!
    endsubroutine calc_pencils_particles
!***********************************************************************
    subroutine dxxp_dt(f,fp,dfp)
!
!  Evolution of particle position.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mpar_loc,mpvar) :: fp, dfp
!
      if (NO_WARN) print*, f, fp, dfp
!
    endsubroutine dxxp_dt
!***********************************************************************
    subroutine dvvp_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle velocity.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      if (NO_WARN) print*, f, df, fp, dfp, ineargrid
!
    endsubroutine dvvp_dt
!***********************************************************************
    subroutine read_particles_init_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_particles_init_pars
!***********************************************************************
    subroutine write_particles_init_pars(unit)
!    
      integer, intent (in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(unit,iostat)
!    
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_particles_run_pars
!***********************************************************************
    subroutine write_particles_run_pars(unit)
!    
      integer, intent (in) :: unit
!
       if (NO_WARN) print*, unit
!
    endsubroutine write_particles_run_pars
!***********************************************************************
    subroutine powersnap_particles(f)
!
!  Calculate power spectra of particle variables.
!
!  01-jan-06/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine powersnap_particles
!***********************************************************************
    subroutine rprint_particles(lreset,lwrite)
!   
!  Read and register print parameters relevant for particles
!
!  22-aug-05/anders: dummy
!
      logical :: lreset, lwr
      logical, optional :: lwrite
! 
!  Write information to index.pro
! 
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      
      if (lwr) then
        write(3,*) 'ixp=', ixp
        write(3,*) 'iyp=', iyp
        write(3,*) 'izp=', izp
        write(3,*) 'ivpx=', ivpx
        write(3,*) 'ivpy=', ivpy
        write(3,*) 'ivpz=', ivpz
      endif
!
    endsubroutine rprint_particles
!***********************************************************************

endmodule Particles
