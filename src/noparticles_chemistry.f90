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

  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
  use Particles_mpicomm

  implicit none

  include 'particles_chemistry.h'

  real :: mol_mass_carbon=12.0
  integer :: N_surface_reactions=0, N_adsorbed_species=0
  integer :: N_species=0, N_surface_reactants = 0
  integer :: N_surface_species=0

  integer :: inuH2=0, inuCO2=0, inuH2O=0, inuCO=0, inuCH4=0, inuO2=0
  integer :: inuCH=0, inuHCO=0, inuCH2=0, inuCH3=0
  integer :: imufree=0, imuadsO=0, imuadsO2=0, imuadsOH=0, imuadsH=0, imuadsCO=0
  integer :: mu=0, mu_prime=0, ac=0, aac=0, nu=0, nu_prime=0
  integer :: jmap=0
  integer :: dependent_reactant = 0
  logical :: lbaum_and_street = .false.
  logical :: lsurface_nopores
  logical :: lpreactions=.false.
  real, dimension(2) :: mass_loss

  contains
! ******************************************************************************
!  09.09.14/jonas : coded

  subroutine get_pchem_info()
  endsubroutine get_pchem_info
!***********************************************************************
    subroutine pencil_criteria_par_chem()
!
!  16.09.2015/jonas + nils: coded
!
    endsubroutine pencil_criteria_par_chem
!***********************************************************************
    subroutine calc_pencils_par_chem(f,p)
!
!  16.09.2015/jonas + nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_par_chem
! ******************************************************************************
!  09.09.14/jonas : coded

  subroutine get_R_c_hat(var,start,end)
    real, dimension(:), intent(out) :: var
    integer :: start, end

    call keep_compiler_quiet(var)
    call keep_compiler_quiet(start)
    call keep_compiler_quiet(end)
  endsubroutine get_R_c_hat
! ******************************************************************************
!  09.09.14/jonas : coded

  subroutine get_R_j_hat(var)
    real, dimension(mpar_loc,N_adsorbed_species), intent(out) :: var

    call keep_compiler_quiet(var)
  endsubroutine get_R_j_hat
! ******************************************************************************
!  09.09.14/jonas : coded

  subroutine get_mod_surf_area(var,fp,irhopswarm,iap)
    real, dimension(mpar_loc), intent(out) :: var
    real, dimension(mpar_loc,mparray) :: fp
    integer :: irhopswarm, iap

    call keep_compiler_quiet(var)
    call keep_compiler_quiet(fp)
    call keep_compiler_quiet(iap)
    call keep_compiler_quiet(irhopswarm)
  endsubroutine get_mod_surf_area
! ******************************************************************************
!  09.09.14/jonas : coded

  subroutine get_St(var,start,end)
    real, dimension(:), intent(out) :: var
    integer :: start, end

    call keep_compiler_quiet(var)
    call keep_compiler_quiet(start)
    call keep_compiler_quiet(end)
  endsubroutine get_St
! ******************************************************************************
! 09.09.14/jonas : coded

  function count_reactions()
    integer :: count_reactions
    count_reactions = 0
  endfunction count_reactions
! ******************************************************************************
  function find_species()
    integer :: find_species
    find_species = 0
  endfunction find_species
! ******************************************************************************
  function count_max_elements()
    integer :: count_max_elements
    count_max_elements = 0
  endfunction count_max_elements
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine get_conversion()
  endsubroutine get_conversion
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine calc_St()
  endsubroutine calc_St
! ******************************************************************************
!  06-oct-14/jonas:coded

  subroutine calc_pchemistry_pencils(f,fp,p,ineargrid)
    real, dimension(mpar_loc,mparray) :: fp
    real, dimension(mx,my,mz,mfarray)  :: f
    integer, dimension(:,:) :: ineargrid
    type (pencil_case) :: p

    call keep_compiler_quiet(f)
    call keep_compiler_quiet(fp)
    call keep_compiler_quiet(ineargrid)
    call keep_compiler_quiet(p)
  endsubroutine calc_pchemistry_pencils
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine calc_mod_surf_area()
  endsubroutine calc_mod_surf_area
! ******************************************************************************

!  19.09.2014/Jonas:coded

  subroutine calc_get_mod_surf_area()
  endsubroutine calc_get_mod_surf_area
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine calc_enthalpy_of_reaction()
  endsubroutine calc_enthalpy_of_reaction
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine calc_entropy_of_reaction()
  endsubroutine calc_entropy_of_reaction
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine calc_conversion()
  endsubroutine calc_conversion
! ******************************************************************************
!  19.09.2014/Jonas:coded

  subroutine calc_R_c_hat()
  endsubroutine calc_R_c_hat
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine create_dependency()
  endsubroutine create_dependency
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine create_ad_sol_lists()
  endsubroutine create_ad_sol_lists
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine create_occupancy()
  endsubroutine create_occupancy
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine create_dngas()
  endsubroutine create_dngas
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine create_stoc()
  endsubroutine create_stoc
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine get_ac()
  endsubroutine get_ac
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine get_part()
  endsubroutine get_part
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine get_reactants()
  endsubroutine get_reactants
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine get_RR_hat()
  endsubroutine get_RR_hat
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine sort_compounds()
  endsubroutine sort_compounds
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine register_particles_chem()
  endsubroutine register_particles_chem
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine calc_RR_hat()
  endsubroutine calc_RR_hat
! ******************************************************************************
! 30.09.2014/Jonas:coded

  subroutine calc_ndot_mdot_R_j_hat()
  endsubroutine calc_ndot_mdot_R_j_hat
! ******************************************************************************
!  30.09.2014/jonas:coded

  subroutine calc_surf_enthalpy()
  endsubroutine calc_surf_enthalpy
! ******************************************************************************
!  30.09.2014/jonas:coded

  subroutine calc_surf_entropy()
  endsubroutine calc_surf_entropy
! ******************************************************************************
!  30.09.2014/jonas:coded

  subroutine calc_ads_enthalpy(fp)
    real, dimension(:,:) :: fp

    call keep_compiler_quiet(fp)
  endsubroutine calc_ads_enthalpy
! ******************************************************************************
!  30.09.2014/jonas:coded

  subroutine calc_ads_entropy(fp)
    real, dimension(:,:) :: fp

    call keep_compiler_quiet(fp)
  endsubroutine calc_ads_entropy
! ******************************************************************************
!  06-oct-14/jonas: coded

  subroutine cleanup_chemistry_pencils()
  endsubroutine cleanup_chemistry_pencils
! ******************************************************************************
  subroutine read_particles_chem_init_pars(iostat)
    integer, intent(out) :: iostat

    iostat = 0
  endsubroutine read_particles_chem_init_pars
! ******************************************************************************
  subroutine write_particles_chem_init_pars(unit)
    integer, intent(in) :: unit

    call keep_compiler_quiet(unit)
  endsubroutine write_particles_chem_init_pars
! ******************************************************************************
  subroutine read_particles_chem_run_pars(iostat)
    integer, intent(out) :: iostat

    iostat = 0
  endsubroutine read_particles_chem_run_pars
! ******************************************************************************
  subroutine write_particles_chem_run_pars(unit)
    integer, intent(in) :: unit

    call keep_compiler_quiet(unit)
  endsubroutine write_particles_chem_run_pars
! ******************************************************************************
!  07-oct-14/jonas: coded

  subroutine calc_rho_p()
  endsubroutine calc_rho_p
! ******************************************************************************
!  11-nov-2014/jonas: coded

  subroutine get_q_reac(var)
    real, dimension(:) :: var

    call keep_compiler_quiet(var)
  endsubroutine get_q_reac
! ******************************************************************************
!  11-nov-2014/jonas: coded

  subroutine get_Nusselt(var)
    real, dimension(:) :: var

    call keep_compiler_quiet(var)
  endsubroutine get_Nusselt
! ******************************************************************************
!  11-nov-2014/jonas: coded

  subroutine get_mass_chemistry(var1,var2,var3)
    real, dimension(:) :: var1,var2
    real, dimension(:,:) :: var3
!
    call keep_compiler_quiet(var1)
    call keep_compiler_quiet(var2)
    call keep_compiler_quiet(var3)
!
  endsubroutine get_mass_chemistry
! ******************************************************************************
!  11-nov-2014/jonas: coded

  subroutine get_surface_chemistry()
  endsubroutine get_surface_chemistry
! ******************************************************************************
!  11-nov-2014/jonas: coded

  subroutine get_adsorbed_chemistry()
  endsubroutine get_adsorbed_chemistry
! ******************************************************************************
!  11-nov-2014/jonas: coded

  subroutine get_radius_chemistry(var1,var2)
    real, dimension(:) :: var1,var2
!
    call keep_compiler_quiet(var1)
    call keep_compiler_quiet(var2)
!
  endsubroutine get_radius_chemistry
! ******************************************************************************
!  11-nov-2014/jonas: coded

  subroutine get_temperature_chemistry(var1,var2)
    real, dimension(:) :: var1,var2
!
    call keep_compiler_quiet(var1)
    call keep_compiler_quiet(var2)
!
  endsubroutine get_temperature_chemistry
! ******************************************************************************
  subroutine particles_chemistry_clean_up()
  endsubroutine particles_chemistry_clean_up
! ******************************************************************************
endmodule Particles_chemistry
