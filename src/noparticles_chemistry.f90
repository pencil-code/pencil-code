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
!***************************************************************!
!  Particle independent variables below here                    !
!***************************************************************!
  real, dimension(:), allocatable :: reaction_order
  real, dimension(:), allocatable :: uscale,fscale,constr
  real, dimension(:), allocatable :: qk_reac
  real, dimension(:), allocatable :: effectiveness_factor_reaction
  real, dimension(:), allocatable :: effectiveness_factor_old
  real, dimension(:), allocatable :: mass_trans_coeff_reactants
  real, dimension(:), allocatable :: diff_coeff_reactants
  real, dimension(:,:), allocatable :: nu, nu_prime, mu, mu_prime
  real, dimension(:), allocatable :: B_k, ER_k, ac, aac, RR_method
  real, dimension(:), allocatable :: site_occupancy, sigma_k,dngas
  real, allocatable, dimension(:) :: omega_pg_dbl
  real, dimension(50) :: reaction_enhancement=1
  character(3), dimension(:), allocatable, save :: reaction_direction
  character(3), dimension(:), allocatable ::flags
  character(10), dimension(:,:), allocatable, save :: part
  character(10), dimension(:), allocatable :: solid_species 
  character(10), dimension(50) :: species_name,adsorbed_species_names
  integer, dimension(:), allocatable :: dependent_reactant
  integer, dimension(:), allocatable :: j_of_inu
  integer :: N_species, N_reactions 
  integer :: N_surface_reactions, N_adsorbed_species
  integer :: N_surface_reactants, N_surface_species
  integer :: inuH2,inuCO2,inuH2O,inuCO,inuCH4,inuO2,inuCH,inuHCO,inuCH2,inuCH3
  integer :: imufree, imuadsO, imuadsO2, imuadsOH, imuadsH, imuadsCO,iads_end
  real :: f_RPM, St_save, Sg_save
  real :: x_s_total, rho_part
  real :: effectiveness_factor
  real :: effectiveness_factor_timeaver=1.
  real :: eta_int=0.,delta_rho_surface=0.
  logical :: first_pchem=.true.
!
!  Some physical constants
!
  real :: mol_mass_carbon=12.0
  real :: Sgc_init=3e5
  real :: struct_par=4.6
!
!  is already in the code (R_CGS), with ergs as unit!!!
!
  real :: gas_constant=8314.0 ![J/kmol/K]
!
!*********************************************************************!
!               Particle dependent variables below here               !
!*********************************************************************!
!  
  real, dimension(:), allocatable :: St
  real, dimension(:,:), allocatable :: heating_k,entropy_k,R_c_hat
  real, dimension(:,:), allocatable :: mdot_ck,RR_hat,K_k,x_surface
  real, dimension(:,:), allocatable :: thiele
  real, dimension(:,:), allocatable :: ndot, Cs
  real, dimension(:,:), allocatable :: X_infty_reactants, St_array
  real, dimension(:), allocatable :: initial_mass, initial_radius
  real, dimension(:), allocatable :: initial_density
  real, dimension(:,:), allocatable :: adsorbed_species_enthalpy
  real, dimension(:,:), allocatable :: surface_species_enthalpy
  real, dimension(:,:), allocatable :: adsorbed_species_entropy
  real, dimension(:,:), allocatable :: surface_species_entropy
  real, dimension(:), allocatable :: Qh,Qc,Qreac,Qrad
  real, dimension(:), allocatable :: A_p,ndot_total  
  real, dimension(:), allocatable :: Particle_temperature  
  real, dimension(:), allocatable :: mod_surf_area
  real, dimension(:), allocatable :: St_init,rho_p_init,m_p_init
!
  private
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
    subroutine get_R_c_hat(var,fp,iap)
!
!  09.09.14/jonas : coded 
!    
      real, dimension(mpar_loc), intent(out) :: var 
      real, dimension(mpar_loc,mpvar) :: fp
      integer :: iap
!
      call keep_compiler_quiet(var)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(iap)

    end subroutine get_R_c_hat
!***********************************************************************
    subroutine get_R_j_hat(var)
!
!  09.09.14/jonas : coded 
!    
      real, dimension(mpar_loc,N_adsorbed_species), intent(out) :: var 
!
!      call keep_compiler_quiet(var)
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
    subroutine get_St(var,fp,iap)
!
!  09.09.14/jonas : coded 
!        
      real, dimension(mpar_loc),intent(out) :: var
      real, dimension(mpar_loc,mpvar) :: fp
      integer :: irhopswarm,iap
!
      call keep_compiler_quiet(var)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(iap)
!
   end subroutine get_St
!***********************************************************************
    subroutine get_conversion(var,fp,irhopswarm,iap)
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
!
  end function count_reactions
!***************************************************

!*********************************************************************** 
integer function find_species(species,unique_species,nlist)
!
   implicit none
!
    integer :: i,nlist
    character(len=*) :: species
    character(len=*), dimension(:) :: unique_species
!
!      call keep_compiler_quiet(species)
!      call keep_compiler_quiet(unique_species)
!      call keep_compiler_quiet(nlist)
!
  end function find_species
!**********************************************************************
!**********************************************************************
  end module Particles_chemistry
