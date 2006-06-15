! $Id: noparticles_selfgravity.f90,v 1.1 2006-06-15 19:34:43 ajohan Exp $
!
!  This module takes care of everything related to particle self-gravity.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MAUX CONTRIBUTION 4
! CPARAM logical, parameter :: lparticles_selfgravity=.true.
!
!***************************************************************
module particles_selfgravity

  use Cdata
  use Particles_cdata
  use Particles_sub

  implicit none

  include 'particles_selfgravity.h'

  contains

!***********************************************************************
    subroutine register_particles_selfgravity()
!
!  Set up indices for access to the fp and dfp arrays.
!
!  14-jun-06/anders: dummy
!
    endsubroutine register_particles_selfgravity
!***********************************************************************
    subroutine initialize_particles_selfgravity(lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  14-jun-06/anders: dummy
!
      logical :: lstarting
!
      if (NO_WARN) print*, lstarting
!
    endsubroutine initialize_particles_selfgravity
!***********************************************************************
    subroutine calc_selfpotential_particles(f,rhs_poisson,rhs_poisson_const,    lcontinued)
!
!  Calculate the potential of the dust particles.
!
!  13-jun-06/anders: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx,ny,nz) :: rhs_poisson
      real :: rhs_poisson_const
      logical :: lcontinued
!
      if (NO_WARN) print*, f, rhs_poisson, rhs_poisson_const, lcontinued
!
    endsubroutine calc_selfpotential_particles
!***********************************************************************
    subroutine dvvp_dt_selfgravity(f,df,fp,dfp,ineargrid)
!
!  Add self-gravity to particle equation of motion.
!
!  14-jun-06/anders: dummy
!
      use Messages, only: fatal_error
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, fp, dfp
!
      if (NO_WARN) print*, f, df, fp, dfp, ineargrid
!
    endsubroutine dvvp_dt_selfgravity
!***********************************************************************
    subroutine read_particles_selfg_init_pars(unit,iostat)
!
!  14-jun-06/anders: dummy
!
      integer, intent (in) :: unit, iostat 
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_particles_selfg_init_pars
!***********************************************************************
    subroutine write_particles_selfg_init_pars(unit)
!
!  14-jun-06/anders: dummy
!
      integer, intent (in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_selfg_init_pars
!***********************************************************************
    subroutine read_particles_selfg_run_pars(unit,iostat)
!
!  14-jun-06/anders: dummy
!
      integer, intent (in) :: unit, iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_particles_selfg_run_pars
!***********************************************************************
    subroutine write_particles_selfg_run_pars(unit)
!
!  14-jun-06/anders: dummy
!    
      integer, intent (in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_selfg_run_pars
!***********************************************************************
    subroutine rprint_particles_selfgravity(lreset,lwrite)
!   
!  Read and register print parameters relevant for particle self-gravity.
!
!  14-jun-06/anders: dummy
!
      use Cdata
!
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write information to index.pro
! 
      if (lwr) then
        write(3,*) 'igpotselfx=0'
        write(3,*) 'igpotselfy=0'
        write(3,*) 'igpotselfz=0'
        write(3,*) 'irhop     =0'
      endif
!
      if (NO_WARN) print*, lreset
!
    endsubroutine rprint_particles_selfgravity
!***********************************************************************

endmodule particles_selfgravity
