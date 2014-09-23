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
    subroutine register_indep_pchem()
!
!  09.09.14/jonas : coded
!
    end subroutine register_indep_pchem
!***********************************************************************
    subroutine register_dep_pchem()
!
!  09.09.14/jonas : coded
!
    end subroutine register_dep_pchem
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
!***********************************************************************
    subroutine get_conversion(var,fp)
 !
!  09.09.14/jonas : coded
!
      real, dimension(mpar_loc), intent(out) :: var
      real, dimension(mpar_loc,mpvar) :: fp
!
      call keep_compiler_quiet(var)
      call keep_compiler_quiet(fp)
!
   end subroutine get_conversion
!***********************************************************************
integer function count_max_elements(inputfile)
  character(*) :: inputfile
!
!  09.09.14/jonas : coded
!
     call keep_compiler_quiet(inputfile)
      count_max_elements=0
!
end function count_max_elements
!**************************************************
  integer function count_reactions(inputfile)
  character(*) :: inputfile
!
!  09.09.14/jonas : coded
!
      call keep_compiler_quiet(inputfile)
      count_reactions=0
!
  end function count_reactions
!***********************************************************************
integer function find_species(species,unique_species,nlist)
!
   implicit none
!
    integer :: i,nlist
    character(len=*) :: species
    character(len=*) :: unique_species
!
    call keep_compiler_quiet(species)
    call keep_compiler_quiet(unique_species)
    call keep_compiler_quiet(nlist)
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
  end module Particles_chemistry
