! $Id$
!
!  This module contains all the main structure needed for particles.
!
!***********************************************************************
module Particles_main

  use Cdata
  use Particles_cdata

  implicit none

  include 'particles_main.h'

  contains

!***********************************************************************
    subroutine particles_register_modules()
!
!  Register particle modules.
!
!  22-aug-05/anders: dummy
!
    endsubroutine particles_register_modules
!***********************************************************************
    subroutine particles_rprint_list(lreset)
!
!  Read names of diagnostic particle variables to print out during run.
!
!  22-aug-05/anders: dummy
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
!  22-aug-05/anders: dummy
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
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine particles_init
!***********************************************************************
    subroutine particles_read_snapshot(filename)
!
!  Read particle snapshot from file.
!
!  22-aug-05/anders: dummy
!
      character (len=*) :: filename
!
      if (NO_WARN) print*, filename
!
    endsubroutine particles_read_snapshot
!***********************************************************************
    subroutine particles_write_snapshot(chsnap,enum,flist)
!
!  Write particle snapshot to file.
!
!  22-aug-05/anders: dummy
!
      logical :: enum
      character (len=*) :: chsnap,flist
      optional :: flist
!
      if (NO_WARN) print*, chsnap, enum, flist
!
    endsubroutine particles_write_snapshot
!***********************************************************************
    subroutine particles_write_dsnapshot(chsnap,enum,flist)
!
!  Write particle derivative snapshot to file.
!
!  22-aug-05/anders: dummy
!
      logical :: enum
      character (len=*) :: chsnap,flist
      optional :: flist
!
      if (NO_WARN) print*, chsnap, enum, flist
!
    endsubroutine particles_write_dsnapshot
!***********************************************************************
    subroutine particles_write_pdim(filename)
!
!  Write npar and mpvar to file.
!
!  22-aug-05/anders: dummy
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
!  22-aug-05/anders: dummy
!
    endsubroutine particles_timestep_first
!***********************************************************************
    subroutine particles_timestep_second()
!
!  Time evolution of particle variables.
!
!  22-aug-05/anders: dummy
!
    endsubroutine particles_timestep_second
!***********************************************************************
    subroutine particles_boundconds(f)
!
!  Particle boundary conditions and parallel communication.
!
!  16-feb-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine particles_boundconds
!***********************************************************************
    subroutine particles_pencil_criteria()
!
!  Request pencils for particles.
!
!  20-apr-06/anders: dummy
!
    endsubroutine particles_pencil_criteria
!***********************************************************************
    subroutine particles_pencil_interdep(lpencil_in)
!
!  Calculate particle pencils.
!
!  16-feb-06/anders: dummy
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in
!
    endsubroutine particles_pencil_interdep
!***********************************************************************
    subroutine particles_calc_pencils(f,p)
!
!  Calculate particle pencils.
!
!  16-feb-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (NO_WARN) print*, f, p
!
    endsubroutine particles_calc_pencils
!***********************************************************************
    subroutine particles_calc_selfpotential(f,rhs_poisson,rhs_poisson_const,lcontinued)
!
!  Calculate the potential of the dust particles.
!
!  13-jun-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: rhs_poisson
      real :: rhs_poisson_const
      logical :: lcontinued
!
      if (NO_WARN) print*, f, rhs_poisson, rhs_poisson_const, lcontinued
!
    endsubroutine particles_calc_selfpotential
!***********************************************************************
    subroutine particles_before_boundary(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine particles_before_boundary
!***********************************************************************
    subroutine particles_special
!
      real :: dummy_=0.
!
      if (NO_WARN) print*,dummy_
!
    endsubroutine particles_special
!***********************************************************************
    subroutine particles_pde(f,df)
!
!  Dynamical evolution of particle variables.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      if (NO_WARN) print*, f, df
!
    endsubroutine particles_pde
!***********************************************************************
    subroutine particles_pde_pencil(f,df,p)
!
!  Dynamical evolution of particle variables.
!
!  20-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      if (NO_WARN) print*, f, df, p
!
    endsubroutine particles_pde_pencil
!***********************************************************************
    subroutine particles_create_sinks(f)
!
! Fetch fp and ineargrid to create_sinks 
! 
! 14-mar-08/wlad: dummy 
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*,f
!
    endsubroutine particles_create_sinks
!***********************************************************************
    subroutine read_particles_init_pars_wrap(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_particles_init_pars_wrap
!***********************************************************************
    subroutine write_particles_init_pars_wrap(unit)
!
      integer, intent (in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_init_pars_wrap
!***********************************************************************
    subroutine read_particles_run_pars_wrap(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_particles_run_pars_wrap
!***********************************************************************
    subroutine write_particles_run_pars_wrap(unit)
!
      integer, intent (in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_particles_run_pars_wrap
!***********************************************************************
    subroutine particles_powersnap(f)
!
!  Calculate power spectra of particle variables.
!
!  01-jan-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine particles_powersnap
!***********************************************************************
    subroutine get_slices_particles(f,slices)
!
!  Write slices for animation of particle variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      if (NO_WARN) print*, f, slices%ready
!
    endsubroutine get_slices_particles
!***********************************************************************
    subroutine particles_doprepencil_calc(f,ivar1,ivar2)
!
!  12-aug-08/kapelrud: dummy coded
!
      real,dimension(mx,my,mz,mfarray),intent(inout) :: f
      integer, intent(out) :: ivar1, ivar2
!
      ivar1 = 0
      ivar2 = 0
!
    endsubroutine particles_doprepencil_calc
!***********************************************************************

endmodule Particles_main
