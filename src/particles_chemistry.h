!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: register_particles_surfchem
  public :: register_indep_psurfchem
  public :: register_dep_psurfchem
  public :: get_St, calc_St
  public :: get_conversion, calc_conversion
  public :: get_R_c_hat, calc_R_c_hat
  public :: get_mod_surf_area, calc_mod_surf_area
  public :: get_species_list
  public :: get_pchem_info

  public :: find_species
  public :: find_entropy_of_reaction
  public :: find_enthalpy_of_reaction
  public :: calc_surf_enthalpy
  public :: calc_surf_entropy
  public :: mol_mass_carbon

  public :: read_particles_surf_init_pars
  public :: write_particles_surf_init_pars
  public :: read_particles_surf_run_pars
  public :: write_particles_surf_run_pars
  public :: rprint_particles_surf
  public :: init_particles_surf

