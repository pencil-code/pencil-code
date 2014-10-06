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
    subroutine get_pchem_info()
!
!  09.09.14/jonas : coded
!
    end subroutine get_pchem_info
!***********************************************************************
    subroutine get_R_c_hat(var,start,end)
!
!  09.09.14/jonas : coded
!
      real, dimension(:), intent(out) :: var
      integer :: start,end
!
      call keep_compiler_quiet(var)
      call keep_compiler_quiet(start)
      call keep_compiler_quiet(end)
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
    subroutine get_St(var,start,end)
!
!  09.09.14/jonas : coded
!
      real, dimension(:),intent(out) :: var
      integer :: start, end
!
      call keep_compiler_quiet(var)
      call keep_compiler_quiet(start)
      call keep_compiler_quiet(end)
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
    find_species=0
!
  end function find_species
!*********************************************************************
  integer function count_max_elements()
!
    count_max_elements=0
!
  end function count_max_elements
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
  subroutine calc_chemistry_pencils(f,fp)
!
!  06-oct-14/jonas:coded
!
    real, dimension(mpar_loc,mpvar)  :: fp
    real, dimension(mx,my,mz,mfarray)  :: f
!
  call keep_compiler_quiet(f)
  call keep_compiler_quiet(fp)
!
  end subroutine calc_chemistry_pencils
!***********************************************************************
  subroutine calc_mod_surf_area()
!
!  19.09.2014/Jonas:coded
!
  end subroutine calc_mod_surf_area
!***********************************************************************
  subroutine calc_enthalpy_of_reaction()
!
!  19.09.2014/Jonas:coded
!
  end subroutine calc_enthalpy_of_reaction
!***********************************************************************
  subroutine calc_entropy_of_reaction()
!
!  19.09.2014/Jonas:coded
!
  end subroutine calc_entropy_of_reaction
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
 subroutine create_dependency()
!
! 30.09.2014/Jonas:coded
!
  end subroutine create_dependency
!***********************************************************************
 subroutine create_ad_sol_lists()
!
! 30.09.2014/Jonas:coded
!
  end subroutine create_ad_sol_lists
!***********************************************************************
 subroutine create_occupancy()
!
! 30.09.2014/Jonas:coded
!
  end subroutine create_occupancy
!***********************************************************************
 subroutine create_dngas()
!
! 30.09.2014/Jonas:coded
!
  end subroutine create_dngas
!***********************************************************************
 subroutine create_stoc()
!
! 30.09.2014/Jonas:coded
!
  end subroutine create_stoc
!***********************************************************************
 subroutine get_ac()
!
! 30.09.2014/Jonas:coded
!
  end subroutine get_ac
!***********************************************************************
 subroutine get_part()
!
! 30.09.2014/Jonas:coded
!
  end subroutine get_part
!***********************************************************************
 subroutine get_reactants()
!
! 30.09.2014/Jonas:coded
!
  end subroutine get_reactants
!***********************************************************************
 subroutine get_RR_hat()
!
! 30.09.2014/Jonas:coded
!
  end subroutine get_RR_hat
!***********************************************************************
 subroutine get_total_carbon_sites()
!
! 30.09.2014/Jonas:coded
!
  end subroutine get_total_carbon_sites
!***********************************************************************
 subroutine sort_compounds()
!
! 30.09.2014/Jonas:coded
!
  end subroutine sort_compounds
!***********************************************************************
 subroutine register_particles_chem()
!
! 30.09.2014/Jonas:coded
!
  end subroutine register_particles_chem
!***********************************************************************
 subroutine calc_RR_hat()
!
! 30.09.2014/Jonas:coded
!
 end subroutine calc_RR_hat
!***********************************************************************
 subroutine calc_ndot_mdot_R_j_hat()
!
! 30.09.2014/Jonas:coded
!
 end subroutine calc_ndot_mdot_R_j_hat
!***********************************************************************
  subroutine calc_surf_enthalpy()
!
!  30.09.2014/jonas:coded
!
  end subroutine calc_surf_enthalpy
!**************************************************************
  subroutine calc_surf_entropy()
!
!  30.09.2014/jonas:coded
!
  end subroutine calc_surf_entropy
!**************************************************************
  subroutine calc_ads_enthalpy(fp)
!
!  30.09.2014/jonas:coded
!
    real, dimension(:,:) :: fp
!
    call keep_compiler_quiet(fp)
!
  end subroutine calc_ads_enthalpy
!***********************************************************************
  subroutine calc_ads_entropy(fp)
!
!  30.09.2014/jonas:coded
!
    real, dimension(:,:) :: fp
!
    call keep_compiler_quiet(fp)
!
  end subroutine calc_ads_entropy
!***********************************************************************
  subroutine cleanup_chemistry_pencils()
!
!  06-oct-14/jonas: coded
!
  end subroutine cleanup_chemistry_pencils
!***********************************************************************
!***************************************************
    subroutine read_particles_chem_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if(present(iostat)) call keep_compiler_quiet(iostat) 
!
    endsubroutine read_particles_chem_init_pars
!***********************************************************************
    subroutine write_particles_chem_init_pars(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_chem_init_pars
!***********************************************************************
    subroutine read_particles_chem_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if(present(iostat)) call keep_compiler_quiet(iostat) 
!
    endsubroutine read_particles_chem_run_pars
!***********************************************************************
    subroutine write_particles_chem_run_pars(unit)
!
      integer, intent (in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_chem_run_pars
!***********************************************************************
  end module Particles_chemistry
