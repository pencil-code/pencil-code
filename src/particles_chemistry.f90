! $Id: particles_sub.f90 20578 2013-06-19 07:46:27Z kdittrich85 $
!
!  This module contains useful subroutines for the particle modules.
!  Subroutines that depend on domain decomposition should be put in
!  the Particles_map module.
!
module Particles_chemistry
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_mpicomm
!
  implicit none
!
  real, dimension(:), allocatable :: reaction_order
  real, dimension(:), allocatable :: heating_k,entropy_k
  real, dimension(:), allocatable :: uscale,fscale,constr
  real, dimension(:), allocatable :: mdot_ck,RR_hat,K_k,x_surface
  real, dimension(:), allocatable :: thiele,qk_reac
  real, dimension(:), allocatable :: effectiveness_factor_reaction
  real, dimension(:), allocatable :: effectiveness_factor_old
  real, dimension(:), allocatable :: mass_trans_coeff_reactants, ndot, Cs
  real, dimension(:), allocatable :: X_infty_reactants, St_array
  real, dimension(:), allocatable :: diff_coeff_reactants
  integer, dimension(:), allocatable :: dependent_reactant
  real :: Qh,Qc,Qreac,Qrad,f_RPM, St_save, Sg_save
  real :: x_s_total, rho_part
  real :: effectiveness_factor
  real :: effectiveness_factor_timeaver=1.
  real :: eta_int=0.,delta_rho_surface=0.
  real :: A_p,ndot_total  
  real :: Particle_temperature


!  private
!
!  public :: input_particles, output_particles, boundconds_particles
!
!
  contains
!***********************************************************************
    subroutine register_particles_chemistry()
!
!  Set up required variables etc. for reactive particles
!
!  02-sep-14/nils: coded
!
!  Count number of heterogeneous reactions and elements
!
      N_surface_reactions = count_reactions('mechanics.in')
      N_max_elements = count_max_elements('mechanics.in') + 1 
!
!  Allocate some arrays
!
      allocate(part(N_max_elements,N_surface_reactions))
      allocate(reaction_direction(N_surface_reactions))
!
!  Read the heterogeneous chemical kinetics mechanism
!
      call read_mechanics_file('mechanics.in',part,n_max_elements,&
          reaction_direction)
!
!  Count heterogeneous species
!
      call count_species(part,species,reactants,products,ns,nr,np)
      call count_species_type(species,N_adsorbed_species,&
          N_surface_species,ns)
      call count_species_type(reactants,N_trash,N_surface_reactants,nr)
!
!  Create list of adsorbed and solid species
!
      call create_ad_sol_lists(species,adsorbed_species_names,'ad',ns)
      call sort_compounds(reactants,adsorbed_species_names,N_adsorbed_species,nr)
!
!  Set some indeces (this is hard-coded for now)
!    
      imuadsO =find_species('C(O)',adsorbed_species_names,N_adsorbed_species)
      imuadsO2=find_species('C2(O2)',adsorbed_species_names,N_adsorbed_species)
      imuadsOH=find_species('C(OH)',adsorbed_species_names,N_adsorbed_species)
      imuadsH =find_species('C(H)',adsorbed_species_names,N_adsorbed_species)
      imuadsCO=find_species('C(CO)',adsorbed_species_names,N_adsorbed_species)
      imufree =find_species('Cf',adsorbed_species_names,N_adsorbed_species)
!
! Check if any of the reactions are enhanced
!
      lenhance=.false.
      do i=1,N_surface_reactions
        if (reaction_enhancement(i) .ne. 1) then
          print*,'**************** WARNING! ****************************'
          write(*,'(A5,I2,A25,F10.2)') &
              'Reac ',i,' is enhanced by a factor ',reaction_enhancement(i)
          lenhance=.true.
        endif
      enddo
      if (lenhance) call sleep(4)      
!
! Allocate memory for a number of arrays
!
    allocate(omega_pg_dbl(N_species),STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for omega_pg_dbl')
    allocate(nu(N_surface_species,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for nu')
    allocate(nu_prime(N_surface_species,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for nu_prime')
    allocate(reaction_order(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for reaction_order')
    allocate(B_k(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for B_k')
    allocate(Er_k(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for Er_k')
    allocate(sigma_k(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for sigma_k')
    allocate(RR_method(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for RR_method')
    allocate(heating_k(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for heating_k')
    allocate(entropy_k(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for heating_k')
    allocate(mu(N_adsorbed_species,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for mu')
    allocate(mu_prime(N_adsorbed_species,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for mu_prime')
    allocate(ac(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for ac')
    allocate(aac(N_adsorbed_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for aac')
    allocate(j_of_inu(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for j_of_inu')
    allocate(site_occupancy(N_adsorbed_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for site_occupancy')
    allocate(solid_species(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for solid_species')
    allocate(adsorbed_species_enthalpy(N_adsorbed_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for adsorbed_species_enthalpy')
    allocate(surface_species_enthalpy(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for surface_species_enthalpy')
    allocate(adsorbed_species_entropy(N_adsorbed_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for adsorbed_species_entropy')
    allocate(surface_species_entropy(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for surface_species_entropy')
    allocate(x_s(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for x_s')
    allocate(x_surface(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for x_surface')
    allocate(RR_hat(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for RR_hat')
    allocate(k_k(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for k_k')
    allocate(mass_trans_coeff_reactants(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for mass_trans_coeff_reactants')
    allocate(diff_coeff_reactants(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for diff_coeff_reactants')
    allocate(x_infty_reactants(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for x_infty_reactants')
    allocate(ndot(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for ndot')
    allocate(Cs(N_adsorbed_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for Cs')
    allocate(St_array(N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for St_array')
    allocate(mdot_ck(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for mdot_ck')
    allocate(R_j_hat(N_adsorbed_species-1)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for R_j_hat')
    allocate(dngas(N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
         'Could not allocate memory for dngas')
    allocate(uscale(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for uscale')
    allocate(fscale(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for fscale')
    allocate(constr(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for constr')
    allocate(dependent_reactant(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for dependent_reactant')
    allocate(effectiveness_factor_reaction(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for effectiveness_factor_reaction')
    effectiveness_factor_reaction=1.0
    allocate(effectiveness_factor_old(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for effectiveness_factor_old')
    effectiveness_factor_old=1.
    allocate(thiele(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for thiele')
    allocate(qk_reac(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_solid_phase',&
        'Could not allocate memory for qk_reac')


!
! Define the aac array which gives the amount of carbon in the
! adsorbed species.
!
    aac=0
    if (imuadsCO>0) aac(imuadsCO)=1
!
! Define the Stoichiometric matrixes. Here i=1 correspond to H2O, i=2 is CO2,
! i=3 is H2 and finally i=4 is O2
!
    call create_ad_sol_lists(species,solid_species,'sol',ns)
    call sort_compounds(reactants,solid_species,n_surface_species,nr)

    inuH2O=find_species('H2O',solid_species,n_surface_species)
    inuCO2=find_species('CO2',solid_species,n_surface_species)
    inuH2 =find_species('H2',solid_species,n_surface_species)
    inuO2 =find_species('O2',solid_species,n_surface_species)
    inuCO =find_species('CO',solid_species,n_surface_species)
    inuCH =find_species('CH',solid_species,n_surface_species)
    inuHCO=find_species('HCO',solid_species,n_surface_species)
    inuCH2=find_species('CH2',solid_species,n_surface_species)
    inuCH3=find_species('CH3',solid_species,n_surface_species)
!
! Set number of carbon atoms for each surface species
!
    call produce_ac(ac,solid_species,N_surface_species)
!
!  Set the stoichiometric matrixes
!
    call create_stoc(part,solid_species,nu,.true.,N_surface_species)
    call create_stoc(part,solid_species,nu_prime,.false.,N_surface_species)
    call create_stoc(part,adsorbed_species_names,mu,.true.,N_adsorbed_species)
    call create_stoc(part,adsorbed_species_names,mu_prime,.false.,&
        N_adsorbed_species)
    call create_occupancy(adsorbed_species_names,site_occupancy)
!
! Define which gas phase reactants the given reaction depends on
!
    call create_dependency(nu,dependent_reactant,&
        n_surface_reactions,n_surface_reactants)
!
! Find the mole production of the forward reaction
!
    call create_dngas(nu,nu_prime,dngas)
!
! Define the Arrhenius coefficients. The term ER_k is the activation energy
! divided by the gas constant (R)
!
    call create_arh_param(part,B_k,ER_k,sigma_k)
!
! Define which method to use for calculating RR (and also St)
!  1: the method of Qiao
!  2: the method of Mithcell
!
    call set_RR(part,RR_method)
!
! Some debuging output
!
    if (ldebug) then
      do k=1,N_surface_reactions
        write(*,'(A4," ",I4,3F4.0)') 'nu=',k,nu(1:3,k)
      end do
 !
      do k=1,N_surface_reactions
        write(*,'(A4," ",I4,4F4.0)') 'mu=',k,mu(1:3,k)
      end do
 !
      do k=1,N_surface_reactions
        write(*,'(A4," ",I4,3F4.0)') 'nu_prime=',k,nu_prime(1:3,k)
      end do
 !
      do k=1,N_surface_reactions
        write(*,'(A4," ",I4,4F4.0)') 'mu_prime=',k,mu_prime(1:3,k)
      end do
!
      do k=1,N_surface_reactions
        write(*,'(A12,I4,2E12.5)') 'ER_k, B_k=',&
            k,B_k(k),ER_k(k)/(1e6/gas_constant)
      enddo
      do k=1,N_surface_reactions
        write(*,'(A12,I4,E12.5)') 'Dngas',k,dngas(k)
      end do
!
      print*,'Adsorbed_species_names'
      print*, adsorbed_species_names
      print*,'Solid_species'
      print*, solid_species
!        
      write(*,'(A20," ",10I4)') 'j_of_inu=', j_of_inu
      write(*,'(A20," ",10F4.0)') 'ac=',ac
      write(*,'(A20," ",10F4.0)') 'site_occupancy=',site_occupancy
      write(*,'(A20," ",30I4)') 'dependent_reactant=',dependent_reactant
      write(*,'(A20," ",10E12.5)') 'sigma_k=',sigma_k
        
    endif
!
    end subroutine register_particles_chemistry
!***********************************************************************
  end module Particles_chemistry
