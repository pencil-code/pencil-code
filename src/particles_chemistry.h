!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
!
!  These are routines that are used once or to transport variables to 
!  _surfspec or _adsorbed
!
 
  public :: register_particles_chem
  public :: calc_get_mod_surf_area
  public :: get_pchem_info
  public :: get_reactants
  public :: get_part
!  public :: get_total_carbon_sites
  public :: get_ac
  public :: find_species
!
!  Routines used in the mn loop for each pencil
!
  public :: calc_pchemistry_pencils
  public :: cleanup_chemistry_pencils
!
!  Transport routines used in _radius, _mass,_surfspec and _adsorbed
!
  public :: get_mass_chemistry
  public :: get_radius_chemistry
  public :: get_adsorbed_chemistry
  public :: get_surface_chemistry
  public :: get_temperature_chemistry
!
!  Routines that are used once in startup or per processor
!
  public :: calc_ads_entropy
  public :: calc_ads_enthalpy
  public :: calc_surf_enthalpy
  public :: calc_surf_entropy
  public :: create_stoc
  public :: create_dependency
  public :: create_dngas
  public :: create_occupancy
  public :: create_ad_sol_lists
  public :: sort_compounds
  public :: count_reactions, count_max_elements
  public :: particles_chemistry_clean_up
  public :: calc_pencils_par_chem, calc_diagnostics_particles_chem
  public :: pencil_criteria_par_chem
  public :: rprint_particles_chem
!
!  Obligatory routines for reading in of the start and run namelists
!
  public :: read_particles_chem_init_pars, write_particles_chem_init_pars
  public :: read_particles_chem_run_pars, write_particles_chem_run_pars  
!
!  These are variables that are needed in particles_mass, _radius, 
!  _surfspec and _adsorbed in their pencil variable evolution 
!
  public :: mol_mass_carbon
!
!  These are variables that need to be communicated to _surfspec and _adsorbed
!
  public :: N_surface_reactions, N_adsorbed_species,N_species
  public :: N_surface_reactants, N_surface_species
  public :: inuH2, inuCO2, inuH2O, inuCO, inuCH4, inuO2
  public :: inuCH, inuHCO, inuCH2, inuCH3
  public :: imufree, imuadsO, imuadsO2, imuadsOH, imuadsH, imuadsCO
  public :: mu, mu_prime, ac, aac, nu, nu_prime
  public :: jmap
  public :: dependent_reactant
  public :: lbaum_and_street
  public :: lsurface_nopores
  public :: lpreactions


  

