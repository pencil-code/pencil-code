!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
!
!  These are routines
!
 
  public :: register_particles_chem
  public :: get_St
  public :: get_conversion
  public :: get_R_c_hat
!  public :: get_mod_surf_area
  public :: get_pchem_info
  public :: get_reactants
  public :: get_part
!  public :: get_total_carbon_sites
  public :: get_RR_hat
  public :: get_ac
  public :: get_q_reac
  public :: get_nusselt

  public :: find_species

  public :: calc_pchemistry_pencils
  public :: calc_entropy_of_reaction
  public :: calc_enthalpy_of_reaction
  public :: calc_RR_hat
  public :: calc_ndot_mdot_R_j_hat
  public :: calc_mod_surf_area
  public :: calc_R_c_hat
  public :: calc_conversion
  public :: calc_St
  public :: calc_ads_entropy
  public :: calc_ads_enthalpy
  public :: calc_surf_enthalpy
  public :: calc_surf_entropy
  public :: create_stoc
  public :: create_dependency
  public :: create_dngas
  public :: create_occupancy
  public :: create_ad_sol_lists

  
  public :: cleanup_chemistry_pencils
  
  public :: sort_compounds

  public :: count_reactions, count_max_elements

  public :: mol_mass_carbon

  public :: read_particles_chem_init_pars, write_particles_chem_init_pars
  public :: read_particles_chem_run_pars, write_particles_chem_run_pars
!
!  These are variables
!
  public :: N_surface_reactions, N_adsorbed_species,N_species
  public :: N_surface_reactants, N_surface_species
  public :: inuH2, inuCO2, inuH2O, inuCO, inuCH4, inuO2
  public :: inuCH, inuHCO, inuCH2, inuCH3
  public :: imufree, imuadsO, imuadsO2, imuadsOH, imuadsH, imuadsCO
  public :: mu, mu_prime, ac, aac, nu, nu_prime
  public :: jmap


  

