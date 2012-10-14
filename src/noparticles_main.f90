! $Id$
!
!  This module contains all the main structure needed for particles.
!
!***********************************************************************
module Particles_main
!
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_main.h'
!
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
      call keep_compiler_quiet(lreset)
!
    endsubroutine particles_rprint_list
!***********************************************************************
    subroutine particles_initialize_modules(f,lstarting)
!
!  Initialize particle modules.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
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
      call keep_compiler_quiet(f)
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
      call keep_compiler_quiet(filename)
!
    endsubroutine particles_read_snapshot
!***********************************************************************
    subroutine particles_write_snapshot(chsnap,f,enum,flist)
!
!  Write particle snapshot to file.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: enum
      character (len=*) :: chsnap,flist
      optional :: flist
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(chsnap)
      call keep_compiler_quiet(enum)
      call keep_compiler_quiet(present(flist))
!
    endsubroutine particles_write_snapshot
!***********************************************************************
    subroutine particles_write_dsnapshot(chsnap,f)
!
!  Write particle derivative snapshot to file.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=*) :: chsnap
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(chsnap)
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
      call keep_compiler_quiet(filename)
!
    endsubroutine particles_write_pdim
!***********************************************************************
    subroutine particles_write_block(filename)
!
!  Write block domain decomposition parameters to file.
!
!  16-nov-09/anders: dummy
!
      character (len=*) :: filename
!
      call keep_compiler_quiet(filename)
!
    endsubroutine particles_write_block
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
    subroutine particles_load_balance(f)
!
!  Redistribute particles among the processors for better load balancing.
!
!  04-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_load_balance
!***********************************************************************
    subroutine particles_boundconds(f)
!
!  Particle boundary conditions and parallel communication.
!
!  16-feb-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
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
      call keep_compiler_quiet(lpencil_in)
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine particles_calc_pencils
!***********************************************************************
    subroutine particles_calc_selfpotential(f,rhs_poisson,lcontinued)
!
!  Calculate the potential of the dust particles.
!
!  13-jun-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: rhs_poisson
      logical :: lcontinued
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rhs_poisson)
      call keep_compiler_quiet(lcontinued)
!
    endsubroutine particles_calc_selfpotential
!***********************************************************************
    subroutine particles_before_boundary(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_before_boundary
!***********************************************************************
    subroutine particles_special
!
      real :: dummy_=0.
!
      call keep_compiler_quiet(dummy_)
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
!
    endsubroutine particles_pde
!***********************************************************************
    subroutine particles_pde_pencil(f,df,p)
!
!  Dynamical evolution of particle variables in pencils.
!
!  20-apr-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine particles_pde_pencil
!***********************************************************************
    subroutine particles_pde_blocks(f,df)
!
!  Dynamical evolution of particle variables in blocks.
!
!  30-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine particles_pde_blocks
!***********************************************************************
    subroutine particles_create_sinks(f)
!
! Fetch fp and ineargrid to create_sinks
!
! 14-mar-08/wlad: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_create_sinks
!***********************************************************************
    subroutine particles_read_startpars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine particles_read_startpars
!***********************************************************************
    subroutine particles_rparam(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine particles_rparam
!***********************************************************************
    subroutine particles_wparam(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine particles_wparam
!***********************************************************************
    subroutine particles_read_runpars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine particles_read_runpars
!***********************************************************************
    subroutine particles_wparam2(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine particles_wparam2
!***********************************************************************
    subroutine particles_powersnap(f)
!
!  Calculate power spectra of particle variables.
!
!  01-jan-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_powersnap
!***********************************************************************
    subroutine get_slices_particles(f,slices)
!
!  Write slices for animation of Particle variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
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
      call keep_compiler_quiet(f)
!
    endsubroutine particles_doprepencil_calc
!***********************************************************************
    subroutine insert_particles_now(f)
!
! dummy
!
      !
      real, dimension (mx,my,mz,mfarray) :: f
      !
      !
      call keep_compiler_quiet(f)
      !
!
!
    endsubroutine insert_particles_now
!***********************************************************************
    subroutine particles_insert_continuously(f)
!
!  Insert particles continuously, i.e. add particles in
!  the beginning of a time step.
!
!  sep-09/kragset: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_insert_continuously
!***********************************************************************
    subroutine write_dim_particles(datadir)
!
      character (len=*) :: datadir
!
      call keep_compiler_quiet(datadir)
!
    endsubroutine write_dim_particles
!***********************************************************************
    subroutine write_snapshot_particles(snap_directory,f,enum,snapnum)
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: enum
      character (len=*) :: snap_directory
      integer, optional :: snapnum
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(snap_directory)
      call keep_compiler_quiet(enum)
      call keep_compiler_quiet(present(snapnum))
!
    endsubroutine write_snapshot_particles
!***********************************************************************
    subroutine read_snapshot_particles(snap_directory)
!
      character (len=*) :: snap_directory
!
     call keep_compiler_quiet(snap_directory)
!
    endsubroutine read_snapshot_particles
!***********************************************************************
    subroutine particles_cleanup
!
! dummy subroutine
!
    endsubroutine particles_cleanup
!***********************************************************************
endmodule Particles_main
