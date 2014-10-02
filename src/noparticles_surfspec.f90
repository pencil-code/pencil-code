! $Id: noparticles_surfspec.f90 19193 2012-06-30 12:55:46Z wdobler $
!
!  This module takes care of everything related to particle
!  surface fractions.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_surfspec=.false.
!
!***************************************************************
module Particles_surfspec
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Particles_cdata
!
  implicit none
!
  include 'particles_surfspec.h'
!
  contains
!***********************************************************************
  subroutine register_particles_surfspec()
!
!  09.09.14/jonas : coded
!
  end subroutine register_particles_surfspec
!***********************************************************************
  subroutine register_dep_psurfchem()
!
!  19.09.2014/Jonas:coded
!
  end subroutine register_dep_psurfchem
!***********************************************************************
  subroutine register_indep_psurfchem()
!
!  19.09.2014/Jonas:coded
!
  end subroutine register_indep_psurfchem
!**********************************************************************
  subroutine init_particles_surf(f,fp)

    real, dimension(mx,my,mz,mpvar) :: f
    real, dimension(mpar_loc,mpvar) :: fp
!
!  19.09.2014/Jonas:coded
!
    call keep_compiler_quiet(f)
    call keep_compiler_quiet(fp)
!
  end subroutine init_particles_surf
!***********************************************************************
  subroutine read_particles_surf_run_pars()
!
!  19.09.2014/Jonas:coded
!
  end subroutine read_particles_surf_run_pars
!***********************************************************************
  subroutine read_particles_surf_init_pars()
!
!  19.09.2014/Jonas:coded
!
  end subroutine read_particles_surf_init_pars
!***********************************************************************
  subroutine rprint_particles_surf()
!
!  19.09.2014/Jonas:coded
!
  end subroutine rPrint_particles_surf
!***********************************************************************
  subroutine write_particles_surf_init_pars()
!
!  19.09.2014/Jonas:coded
!
  end subroutine write_particles_surf_init_pars
!***********************************************************************
  subroutine write_particles_surf_run_pars()
!
!  19.09.2014/Jonas:coded
!
  end subroutine write_particles_surf_run_pars
!***********************************************************************
  subroutine calc_pchem_factors(fp)
!
!  19.09.2014/jonas:coded
!
    real, dimension(mpar_loc,mpvar) :: fp
!
    call keep_compiler_quiet(fp)
!
  end subroutine calc_pchem_factors
!***********************************************************************
  subroutine initialize_particles_surf(fp,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-sep-14/jonas coded
!  JONAS: needs to be filled with life
!
      real, dimension (mpar_loc,mpvar) :: fp
      logical :: lstarting
!
!  
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(lstarting)
!
  end subroutine initialize_particles_surf
!**************************************************************
end module Particles_surfspec
