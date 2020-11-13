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
    subroutine particles_initialize_modules(f)
!
!  Initialize particle modules.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
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
    subroutine particles_finalize()
!
!  Finalize particle modules.
!
!  01-May-2019/PABourdin: dummy
!
    endsubroutine particles_finalize
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
    subroutine particles_timestep_first(f)
!
!  Setup dfp in the beginning of each itsub.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_timestep_first
!***********************************************************************
    subroutine particles_timestep_second(f)
!
!  Time evolution of particle variables.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_timestep_second
!***********************************************************************
    subroutine split_update_particles(f, dt)
!
!  Wrapper for operator split terms for particle dynamics.
!
!  14-dec-14/ccyang: dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, intent(in) :: dt
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(dt)
!
    endsubroutine split_update_particles
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
    subroutine particles_special_bfre_bdary(f)
!
!  Particle special before boundary.
!  Fetch fp (and fsp) array to special module.
!
!  11-jul-16/wlad: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine particles_special_bfre_bdary
!***********************************************************************
    subroutine particles_special_after_dtsub(f, dtsub)
!
!  Send fp and dfp to Special for processing in the end of a
!  sub-time-step.
!
!  27-feb-18/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, intent(in) :: dtsub
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine particles_special_after_dtsub
!***********************************************************************
    subroutine particles_pde(f,df,p)
!
!  Dynamical evolution of particle variables.
!
!  22-aug-05/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
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
    subroutine read_all_particles_init_pars()
!
    endsubroutine read_all_particles_init_pars
!***********************************************************************
    subroutine write_all_particles_init_pars(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_all_particles_init_pars
!***********************************************************************
    subroutine read_all_particles_run_pars()
!
    endsubroutine read_all_particles_run_pars
!***********************************************************************
    subroutine write_all_particles_run_pars(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_all_particles_run_pars
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
    subroutine write_dim_particles(datadir)
!
      character (len=*) :: datadir
!
      call keep_compiler_quiet(datadir)
!
    endsubroutine write_dim_particles
!***********************************************************************
    subroutine write_snapshot_particles(f,enum,snapnum)
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: enum
      integer, optional :: snapnum
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(enum)
      call keep_compiler_quiet(present(snapnum))
!
    endsubroutine write_snapshot_particles
!***********************************************************************
    subroutine read_snapshot_particles()
!
! dummy subroutine
!
    endsubroutine read_snapshot_particles
!***********************************************************************
    subroutine particles_cleanup
!
! dummy subroutine
!
    endsubroutine particles_cleanup
!***********************************************************************
    subroutine particles_stochastic
!
! dummy subroutine
!
    endsubroutine particles_stochastic
!***********************************************************************
    subroutine fetch_nparloc(dummy)
! dummy subroutine
      integer :: dummy
      call keep_compiler_quiet(dummy)
    endsubroutine fetch_nparloc
!***********************************************************************
    subroutine append_particle_index(label,ilabel)
!
      character (len=*), intent(in) :: label
      integer, intent(out) :: ilabel
!
      call keep_compiler_quiet(label)
      call keep_compiler_quiet(ilabel)
!
    endsubroutine append_particle_index
!*********************************************************************** 
    subroutine fetch_fp_array(fp_aux,dfp_aux,ixw,iyw,izw,ivxw,ivyw,ivzw)
!
! dummy subroutine
!
      real, dimension (mpar_loc,mparray), intent(out) :: fp_aux
      real, dimension (mpar_loc,mpvar), intent(out) :: dfp_aux
      integer, intent(out) :: ixw,iyw,izw,ivxw,ivyw,ivzw
!
      call keep_compiler_quiet(fp_aux)
      call keep_compiler_quiet(dfp_aux)
      call keep_compiler_quiet(ixw,iyw,izw)
      call keep_compiler_quiet(ivxw,ivyw,ivzw)
!
    endsubroutine fetch_fp_array
!***********************************************************************
    subroutine return_fp_array(fp_aux,dfp_aux,flag)
!
! dummy subroutine
!
      real, dimension (mpar_loc,mparray), intent(in) :: fp_aux
      real, dimension (mpar_loc,mpvar), intent(in) :: dfp_aux
      logical, dimension (mpar_loc), optional :: flag
!
      call keep_compiler_quiet(fp_aux)
      call keep_compiler_quiet(dfp_aux)
      call keep_compiler_quiet(flag)
!
    endsubroutine return_fp_array
!***********************************************************************
  endmodule Particles_main
