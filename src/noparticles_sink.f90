! $Id: particles_dust.f90 19206 2012-06-30 21:40:24Z sven.bingert $
!
!  This module takes care of everything related to sink particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! CPARAM logical, parameter :: lparticles_sink=.false.
!
!***************************************************************
module Particles_sink
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Particles_cdata
!
  implicit none
!
  include 'particles_sink.h'
!
  contains
!***********************************************************************
    subroutine register_particles_sink
!
!  07-aug-12/anders: dummy
!
    endsubroutine register_particles_sink
!***********************************************************************
    subroutine initialize_particles_sink(f)
!
!  07-aug-12/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_sink
!***********************************************************************
    subroutine init_particles_sink(f,fp)
!
!  07-aug-12/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine init_particles_sink
!***********************************************************************
    subroutine calc_selfpot_sinkparticles(f,rhs_poisson,fp,ineargrid)
!
!  14-oct-12/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: rhs_poisson
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(rhs_poisson)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine calc_selfpot_sinkparticles
!***********************************************************************
    subroutine create_particles_sink(f,fp,dfp,ineargrid)
!
!  07-aug-12/anders: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine create_particles_sink
!***********************************************************************
    subroutine remove_particles_sink(f,fp,dfp,ineargrid)
!
!  07-aug-12/anders: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension(mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine remove_particles_sink
!***********************************************************************
    subroutine read_particles_sink_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_sink_init_pars
!***********************************************************************
    subroutine write_particles_sink_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_sink_init_pars
!***********************************************************************
    subroutine read_particles_sink_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_sink_run_pars
!***********************************************************************
    subroutine write_particles_sink_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_sink_run_pars
!***********************************************************************
    subroutine rprint_particles_sink(lreset,lwrite)
!
!  11-aug-12/anders: dummy
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_particles_sink
!***********************************************************************
endmodule Particles_sink
