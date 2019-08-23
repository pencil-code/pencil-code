! $Id$
!
!  This module tries to solve for gradient matrix of particle velocities
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 0
! MPAUX CONTRIBUTION 0
!
! CPARAM logical, parameter :: lparticles_grad=.false.
!
!***************************************************************
module Particles_grad
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Particles_cdata
!
  implicit none
!
  include 'particles_grad.h'
!
  contains
!***********************************************************************
    subroutine register_particles_grad
!
!  Set up indices for access to the fp and dfp arrays.
!
!  17-sep-15/dhruba: dummy
!
    endsubroutine register_particles_grad
!***********************************************************************
    subroutine initialize_particles_grad(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  15-sep-15/dhruba: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_grad
!***********************************************************************
    subroutine pencil_criteria_par_grad
!
!  All pencils that the Particles_grad module depends on are specified here.
!
!  17-sep-15/dhruba: dummy
!
    endsubroutine pencil_criteria_par_grad
!***********************************************************************
    subroutine set_particle_grad(f,fp,npar_low,npar_high,ineargrid,init)
!
!  Set radius of new particles.
!
!  18-sep-15/dhruba: dummy
!
      use General, only: random_number_wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer :: npar_low,npar_high
      integer, dimension(mpar_loc,3) :: ineargrid
      logical, optional :: init
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(npar_low)
      call keep_compiler_quiet(npar_high)
      call keep_compiler_quiet(init)
!
      endsubroutine set_particle_grad
!***********************************************************************
    subroutine dsigmap_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of the gradient of particle velocities.
!
!  17-sep-15/dhruba: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
!
      intent (in) :: f
      intent (out) :: dfp
      intent (inout) :: fp

      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dsigmap_dt_pencil
!***********************************************************************
    subroutine dsigmap_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle gradient (actually merely the calculation of diagnostic 
! variables on particles because the actual evolution is calculated in a pencilized manner )
!
!  17-sep-15/dhruba: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc,3) :: ineargrid

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dsigmap_dt
!***********************************************************************
    subroutine read_particles_grad_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat=0
!
    endsubroutine read_particles_grad_init_pars
!***********************************************************************
    subroutine write_particles_grad_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_grad_init_pars
!***********************************************************************
    subroutine read_particles_grad_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_grad_run_pars
!***********************************************************************
    subroutine write_particles_grad_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_grad_run_pars
!***********************************************************************
    subroutine rprint_particles_grad(lreset,lwrite)
!
!  Read and register print parameters relevant for particles grad.
!
!  22-aug-05/anders: coded
!
      use Diagnostics, only: parse_name
!
      logical :: lreset
      logical, optional :: lwrite
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_grad
!***********************************************************************
endmodule Particles_grad
