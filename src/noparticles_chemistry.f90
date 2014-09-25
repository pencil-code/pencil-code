! $Id: noparticles_chemistry.f90 19193 2012-06-30 12:55:46Z wdobler $
!
!  This module takes care of everything related to particle
!  chemistry
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_chemistry=.false.
!
!***************************************************************
module Particles_chemistry
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
  use Particles_mpicomm
  use Particles_radius
!
  implicit none
!
  include 'particles_chemistry.h'
!
  integer :: N_adsorbed_species=0
  real :: mol_mass_carbon=12.0
!
  contains
!***********************************************************************
  subroutine register_particles_surfchem()
!
!  09.09.14/jonas : coded
!
  end subroutine register_particles_surfchem
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
!***********************************************************************
    subroutine get_pchem_info()
!
!  09.09.14/jonas : coded
!
    end subroutine get_pchem_info
!***********************************************************************
    subroutine get_R_c_hat(var,fp)
!
!  09.09.14/jonas : coded
!
      real, dimension(mpar_loc), intent(out) :: var
      real, dimension(mpar_loc,mpvar) :: fp
!
      call keep_compiler_quiet(var)
      call keep_compiler_quiet(fp)
!
    end subroutine get_R_c_hat
!***********************************************************************
    subroutine get_R_j_hat(var)
!
!  09.09.14/jonas : coded
!
      real, dimension(mpar_loc,N_adsorbed_species), intent(out) :: var
!
      call keep_compiler_quiet(var)
!
    end subroutine get_R_j_hat
!***********************************************************************
    subroutine get_mod_surf_area(var,fp,irhopswarm,iap)
!
!  09.09.14/jonas : coded
!
      real, dimension(mpar_loc), intent(out) :: var
      real, dimension(mpar_loc,mpvar) :: fp
      integer :: irhopswarm,iap
!
      call keep_compiler_quiet(var)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(iap)
      call keep_compiler_quiet(irhopswarm)
!
    end subroutine get_mod_surf_area
!***********************************************************************
    subroutine get_St(var,fp)
!
!  09.09.14/jonas : coded
!
      real, dimension(mpar_loc),intent(out) :: var
      real, dimension(mpar_loc,mpvar) :: fp
!
      call keep_compiler_quiet(var)
      call keep_compiler_quiet(fp)
!
   end subroutine get_St
!**************************************************
  integer function count_reactions()
!
!  09.09.14/jonas : coded
!
     count_reactions=0
!
  end function count_reactions
!***********************************************************************
integer function find_species()
!
   implicit none
!
    find_species=0
!
  end function find_species
!**********************************************************************
  subroutine get_species_list(string,list)
!
    character(*) :: string
    character(10) :: list
!
    call keep_compiler_quiet(string)
    call keep_compiler_quiet(list)
!
  end subroutine get_species_list
!**********************************************************************
  subroutine get_conversion()
!
!  19.09.2014/Jonas:coded
!
  end subroutine get_conversion
!***********************************************************************
  subroutine calc_St()
!
!  19.09.2014/Jonas:coded
!
  end subroutine calc_St
!***********************************************************************
  subroutine calc_surf_enthalpy()
!
!  19.09.2014/Jonas:coded
!
  end subroutine calc_surf_enthalpy
!***********************************************************************
  subroutine calc_surf_entropy()
!
!  19.09.2014/Jonas:coded
!
  end subroutine calc_surf_entropy
!***********************************************************************
  subroutine calc_mod_surf_area()
!
!  19.09.2014/Jonas:coded
!
  end subroutine calc_mod_surf_area
!***********************************************************************
  subroutine find_enthalpy_of_reaction()
!
!  19.09.2014/Jonas:coded
!
  end subroutine find_enthalpy_of_reaction
!***********************************************************************
  subroutine find_entropy_of_reaction()
!
!  19.09.2014/Jonas:coded
!
  end subroutine find_entropy_of_reaction
!***********************************************************************
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
  subroutine calc_conversion()
!
!  19.09.2014/Jonas:coded
!
  end subroutine calc_conversion
!***********************************************************************
  subroutine calc_R_c_hat()
!
!  19.09.2014/Jonas:coded
!
  end subroutine calc_R_c_hat
!***********************************************************************
  end module Particles_chemistry
