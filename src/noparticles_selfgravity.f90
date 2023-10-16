! $Id$
!
!  This module takes care of everything related to particle self-gravity.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MAUX CONTRIBUTION 0
! CPARAM logical, parameter :: lparticles_selfgravity=.false.
!
!***************************************************************
module Particles_selfgravity
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Particles_cdata
!
  implicit none
!
  include 'particles_selfgravity.h'
!
  contains
!***********************************************************************
    subroutine register_particles_selfgrav
!
!  Set up indices for access to the fp and dfp arrays.
!
!  14-jun-06/anders: dummy
!
    endsubroutine register_particles_selfgrav
!***********************************************************************
    subroutine initialize_particles_selfgrav(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  14-jun-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_selfgrav
!***********************************************************************
    subroutine calc_selfpotential_particles(f,rhs_poisson,lcontinued)
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
    endsubroutine calc_selfpotential_particles
!***********************************************************************
    subroutine pencil_criteria_par_selfgrav
!
!  All pencils that the Particles_selfgrav module depends on are specified here.
!
!  02-jul-06/anders: dummy
!
    endsubroutine pencil_criteria_par_selfgrav
!***********************************************************************
    subroutine pencil_interdep_par_selfgrav(lpencil_in)
!
!  Interdependency among pencils provided by the Particles_selfgrav module
!  is specified here.
!
!  02-jul-06/anders: dummy
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_par_selfgrav
!***********************************************************************
    subroutine calc_pencils_par_selfgrav(f,p)
!
!  Calculate particle pencils.
!
!  02-jul-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_par_selfgrav
!***********************************************************************
    subroutine dvvp_dt_selfgrav_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Add self-gravity to particle equation of motion.
!
!  14-jun-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, df, p, fp, dfp, ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dvvp_dt_selfgrav_pencil
!***********************************************************************
    subroutine calc_diagnostics_particles_selg(p)

      type (pencil_case) :: p
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_diagnostics_particles_selg
!***********************************************************************
    subroutine dvvp_dt_selfgrav(f,df,fp,dfp,ineargrid)
!
!  Add self-gravity to particle equation of motion.
!
!  14-jun-06/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f, fp, dfp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dvvp_dt_selfgrav
!***********************************************************************
    subroutine read_particles_selfg_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_selfg_init_pars
!***********************************************************************
    subroutine write_particles_selfg_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_selfg_init_pars
!***********************************************************************
    subroutine read_particles_selfg_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_selfg_run_pars
!***********************************************************************
    subroutine write_particles_selfg_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_selfg_run_pars
!***********************************************************************
    subroutine rprint_particles_selfgrav(lreset,lwrite)
!
!  Read and register print parameters relevant for particle self-gravity.
!
!  14-jun-06/anders: dummy
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_selfgrav
!***********************************************************************
endmodule Particles_selfgravity
