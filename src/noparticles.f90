! $Id: noparticles.f90,v 1.3 2005-06-27 00:14:19 mee Exp $
!
!  This module takes care of everything related to particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles=.false.
!
!***************************************************************
module Particles

  use Cdata

  implicit none

  include 'particles.h'

  contains

!***********************************************************************
    subroutine particles_register_modules()
!
!  Register particle modules.
!
!  07-jan-05/anders: dummy coded
!
    endsubroutine particles_register_modules
!***********************************************************************
    subroutine particles_rprint_list(lreset)
!
!  Read names of diagnostic particle variables to print out during run.
!
!  07-jan-05/anders: dummy coded
!
      logical :: lreset
!
      if (NO_WARN) print*, lreset
!
    endsubroutine particles_rprint_list
!***********************************************************************
    subroutine particles_initialize_modules(lstarting)
!
!  Initialize particle modules.
!
!  07-jan-05/anders: dummy coded
!
      logical :: lstarting
!
      if (NO_WARN) print*, lstarting
!
    endsubroutine particles_initialize_modules
!***********************************************************************
    subroutine particles_init(f)
!
!  Set up initial conditios for particle modules.
!
!  07-jan-05/anders: dummy coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      intent (in) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine particles_init
!***********************************************************************
    subroutine particles_read_snapshot(filename)
!
!  Read particle snapshot from file.
!
!  07-jan-05/anders: dummy coded
!
      character (len=*) :: filename
!
      if (NO_WARN) print*, filename
!
    endsubroutine particles_read_snapshot
!***********************************************************************
    subroutine particles_write_snapshot(chsnap,msnap,enum,flist)
!
!  Write particle snapshot to file.
!
!  07-jan-05/anders: dummy coded
!
      integer :: msnap
      logical :: enum
      character (len=*) :: chsnap, flist
      optional :: flist
!
      if (NO_WARN) print*, chsnap, msnap, enum, flist
!
    endsubroutine particles_write_snapshot
!***********************************************************************
    subroutine particles_write_pdim(filename)
!   
!  Write npar and mpvar to file.
!
!  09-jan-05/anders: dummy coded
!
      character (len=*) :: filename
!
      if (NO_WARN) print*, filename
!
    endsubroutine particles_write_pdim
!***********************************************************************
    subroutine particles_timestep_first()
!
!  Setup dfp in the beginning of each itsub.
!
!  07-jan-05/anders: dummy coded
!
    endsubroutine particles_timestep_first
!***********************************************************************
    subroutine particles_timestep_second()
!
!  Time evolution of particle variables.
!
!  07-jan-05/anders: dummy coded
!
    endsubroutine particles_timestep_second
!***********************************************************************
    subroutine particles_pde(f,df)
!
!  07-jan-05/anders: dummy coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      intent (in) :: f, df
!
      if (NO_WARN) print*, f, df
!
    endsubroutine particles_pde
!***********************************************************************
    subroutine read_particles_init_pars(unit,iostat)
!
!  07-jan-05/anders: dummy coded
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
!  07-jan-05/anders: dummy coded
!
      integer, intent (in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_init_pars
!***********************************************************************
    subroutine read_particles_run_pars(unit,iostat)
!
!  07-jan-05/anders: dummy coded
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
!  07-jan-05/anders: dummy coded
!
      integer, intent (in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_run_pars
!***********************************************************************

endmodule Particles
