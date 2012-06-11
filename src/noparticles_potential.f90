! $Id: noparticles_potential.f90  $
!
! The no module for particles potential
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_potential=.false.
!
!***************************************************************
module Particles_potential
!
  use Cdata
  use Particles_cdata
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'particles_potential.h'
!
  contains
!***********************************************************************
    subroutine register_particles_potential()
!
!  Set up indices for access to the fp and dfp arrays
!
!  22-aug-05/anders: dummy
!
    endsubroutine register_particles_potential
!***********************************************************************
    subroutine initialize_particles_potential(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  25-nov-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_particles_potential
!***********************************************************************
    subroutine pencil_criteria_par_potential()
!
!  All pencils that the Particles_radius module depends on are specified here.
!
!  21-nov-06/anders: dummy
!
    endsubroutine pencil_criteria_par_potential
!***********************************************************************
    subroutine get_interparticle_accn(fp,k,interparticle_acceleration)

!
!  dhruba: 
!
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: k
      real, dimension(3) :: interparticle_acceleration
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(k)
      call keep_compiler_quiet(interparticle_acceleration)
!
    endsubroutine get_interparticle_accn
!***********************************************************************
    subroutine read_particles_pot_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_pot_init_pars
!***********************************************************************
    subroutine write_particles_pot_init_pars(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_pot_init_pars
!***********************************************************************
    subroutine read_particles_pot_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_pot_run_pars
!***********************************************************************
    subroutine write_particles_pot_run_pars(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_pot_run_pars
!***********************************************************************
    subroutine rprint_particles_potential(lreset,lwrite)
!
!  Read and register print parameters relevant for particles potential.
!
!  22-aug-05/anders: dummy
!
      logical :: lreset, lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (lwr) then
        write(3,*) 'iap=', iap
      endif
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_particles_potential
!***********************************************************************
endmodule Particles_potential
