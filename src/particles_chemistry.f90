! $Id: Particles_chemistry.f90 21950 2014-07-08 08:53:00Z michiel.lambrechts $
!
!  This module takes care of everything related to reactive particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 0
!
! CPARAM logical, parameter :: lparticles_chemistry=.true.
!
!***************************************************************
module Particles_chemistry
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata, only: iap, irhopswarm, iTp, imp
  use Particles_sub
  use Particles_radius
!
  implicit none
!
  include 'particles_chemistry.h'
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
  real, dimension(:), allocatable :: omega_pg_dbl
  real, dimension(50) :: reaction_enhancement=1
  character(3), dimension(:), allocatable, save :: reaction_direction
  character(3), dimension(:), allocatable ::flags
  character(10), dimension(:,:), allocatable, save :: part
  character(10), dimension(:), allocatable :: solid_species 
  character(10), dimension(50) :: species_name,adsorbed_species_names
  integer, dimension(:), allocatable :: dependent_reactant
  integer, dimension(:), allocatable :: j_of_inu
  integer :: N_species, N_reactions 
  integer :: N_surface_reactions, N_adsorbed_species=1
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
  real, dimension(:,:), allocatable :: heating_k,entropy_k
  real, dimension(:), allocatable :: R_c_hat
  real, dimension(:,:), allocatable :: mdot_ck,RR_hat,K_k,x_surface
  real, dimension(:,:), allocatable :: thiele
  real, dimension(:,:), allocatable :: ndot, Cs
  real, dimension(:,:), allocatable :: X_infty_reactants, St_array
  real, dimension(:), allocatable :: initial_radius
  real, dimension(:), allocatable :: initial_density
  real, dimension(:,:), allocatable :: adsorbed_species_enthalpy
  real, dimension(:,:), allocatable :: surface_species_enthalpy
  real, dimension(:,:), allocatable :: adsorbed_species_entropy
  real, dimension(:,:), allocatable :: surface_species_entropy
  real, dimension(:), allocatable :: Qh,Qc,Qreac,Qrad
  real, dimension(:), allocatable :: A_p_init,ndot_total  
  real, dimension(:), allocatable :: Particle_temperature  
  real, dimension(:), allocatable :: mod_surf_area
  real, dimension(:), allocatable :: St_init
  real, dimension(mpar_loc) :: init_mass,rho_p_init
  real, dimension(mpar_loc) :: conversion
!
  contains
!***********************************************************************
  subroutine register_indep_pchem()
!
      integer :: i,k,nr,ns,stat,dummy,N_surface_species
      logical :: lenhance
      character(10), dimension(40) :: species,reactants
!
      call get_pchem_info(species,'N_surface_species',N_surface_species,'quiet')
!      
      if (nsurfreacspec/=N_surface_species) then
         print*,'N_surface_species: ', N_surface_species
         call fatal_error('register_particles_ads', &
              'wrong size of storage for surface species allocated')
         else
      endif
!
      call get_species_list('solid_species',solid_species)
!      
      call create_ad_sol_lists(species(:ns),adsorbed_species_names,'ad',ns)
      call sort_compounds(reactants,adsorbed_species_names,N_adsorbed_species,nr)
!
!  Increase of npvar according to N_surface_species, which is
!  the concentration of gas phase species at the particle surface
!
      if (N_surface_species>1) then
         isurf = npvar+1
         do i=1,N_surface_species
! JONAS: where do we save this
            pvarname(isurf+i-1)=solid_species(i)
         enddo
         npvar=npvar+N_surface_species-1
         isurf_end=isurf+N_surface_species-1
      else
         call fatal_error('register_particles_ads', &
              'N_surface_species must be > 1')
      endif
!
!  Set some indeces (this is hard-coded for now)
!    
      if (nadsspec>0) then
         imuadsO =find_species('C(O)',adsorbed_species_names,N_adsorbed_species)
         imuadsO2=find_species('C2(O2)',adsorbed_species_names,N_adsorbed_species)
         imuadsOH=find_species('C(OH)',adsorbed_species_names,N_adsorbed_species)
         imuadsH =find_species('C(H)',adsorbed_species_names,N_adsorbed_species)
         imuadsCO=find_species('C(CO)',adsorbed_species_names,N_adsorbed_species)
         imufree =find_species('Cf',adsorbed_species_names,N_adsorbed_species)
      else
      end if
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
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for omega_pg_dbl')
    allocate(nu(N_surface_species,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for nu')
    allocate(nu_prime(N_surface_species,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for nu_prime')
    allocate(reaction_order(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for reaction_order')
    allocate(B_k(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for B_k')
    allocate(Er_k(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for Er_k')
    allocate(sigma_k(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for sigma_k')
    allocate(RR_method(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for RR_method')
    allocate(heating_k(mpar_loc,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for heating_k')
    allocate(entropy_k(mpar_loc,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for heating_k')
    allocate(ac(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for ac')
    allocate(j_of_inu(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for j_of_inu')
    allocate(solid_species(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for solid_species')
    allocate(mass_trans_coeff_reactants(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for mass_trans_coeff_reactants')
    allocate(diff_coeff_reactants(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for diff_coeff_reactants')
    allocate(St_array(mpar_loc,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for St_array')
    allocate(dngas(N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
         'Could not allocate memory for dngas')
    allocate(uscale(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for uscale')
    allocate(fscale(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for fscale')
    allocate(constr(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for constr')
    allocate(dependent_reactant(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for dependent_reactant')
    allocate(effectiveness_factor_reaction(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for effectiveness_factor_reaction')
    effectiveness_factor_reaction=1.0
    allocate(effectiveness_factor_old(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for effectiveness_factor_old')
    effectiveness_factor_old=1.
    allocate(thiele(mpar_loc,N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for thiele')
    allocate(qk_reac(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
         'Could not allocate memory for qk_reac')
    allocate(St(mpar_loc)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
         'Could not allocate memory for St')
    allocate(R_c_hat(mpar_loc)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
         'Could not allocate memory for R_c_hat')
!
!  Allocate only if adsorbed species in mechanism
!
    if(Nadsspec>0) then
       allocate(mu(N_adsorbed_species,N_surface_reactions),STAT=stat)
       if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for mu')
       allocate(mu_prime(N_adsorbed_species,N_surface_reactions)   ,STAT=stat)
       if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for mu_prime')
       allocate(aac(N_adsorbed_species)   ,STAT=stat)
       if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for aac')
       allocate(site_occupancy(N_adsorbed_species)   ,STAT=stat)
       if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for site_occupancy')
    else
    end if
       
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
    call create_ad_sol_lists(species(:ns),solid_species,'sol',ns)
    call sort_compounds(reactants,solid_species,n_surface_species,nr)
!
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
    call get_ac(ac,solid_species,N_surface_species)
!
!  Set the stoichiometric matrixes
!
    call create_stoc(part,solid_species,nu,.true.,N_surface_species)
    call create_stoc(part,solid_species,nu_prime,.false.,N_surface_species)
    call create_stoc(part,adsorbed_species_names,mu,.true.,N_adsorbed_species)
    call create_stoc(part,adsorbed_species_names,mu_prime,.false.,&
        N_adsorbed_species)
    if(nadsspec>0) call create_occupancy(adsorbed_species_names,site_occupancy)
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
! Some debugging output
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
!        
    endif
!
    end subroutine register_indep_pchem
!***********************************************************************
    subroutine register_dep_pchem()
!      
      integer :: stat,dummy
      character(10), dimension(40) :: trash   
!      
      if (first_pchem) then
         call get_pchem_info(trash,'dummy',dummy,'verbose')
         first_pchem = .false.
      else
      endif
!
    allocate(init_mass(mpar_loc),STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for init_mass')
    allocate(x_surface(mpar_loc,N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for x_surface')
    allocate(RR_hat(mpar_loc,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for RR_hat')
    allocate(ndot(mpar_loc,N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for ndot')
    allocate(mdot_ck(mpar_loc,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for mdot_ck')
    allocate(x_infty_reactants(mpar_loc,N_surface_reactants) &
         ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for x_infty_reactants')
    allocate(surface_species_enthalpy(mpar_loc,N_surface_species) &
         ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for surface_species_enthalpy')
    allocate(surface_species_entropy(mpar_loc,N_surface_species) &
         ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for surface_species_entropy')
    allocate(k_k(mpar_loc,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for k_k')
    allocate(A_p_init(mpar_loc)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for A_p_init')
    allocate(St(mpar_loc)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for St')
!
    if (nadsspec>0) then
           allocate(adsorbed_species_entropy(mpar_loc,N_adsorbed_species) &
                ,STAT=stat)    
           if (stat>0) call fatal_error('register_dep_pchem',&
                'Could not allocate memory for adsorbed_species_entropy')
           allocate(adsorbed_species_enthalpy(mpar_loc,N_adsorbed_species) &
                ,STAT=stat)
           if (stat>0) call fatal_error('register_dep_pchem',&
                'Could not allocate memory for adsorbed_species_enthalpy')
           allocate(Cs(mpar_loc,N_adsorbed_species)   ,STAT=stat)
           if (stat>0) call fatal_error('register_dep_pchem',&
                'Could not allocate memory for Cs')
    else
    end if
!
    end subroutine register_dep_pchem
!***********************************************************************
    subroutine get_pchem_info(species,string,variable,talk)
!
!  Set up required variables etc. for reactive particles
!
!  02-sep-14/nils: coded
!
!  Count number of heterogeneous reactions and elements
!
      integer :: ns,nr,np,N_trash,N_max_elements,variable,stat
      character(10), dimension(40) :: species,products,reactants
      character(10), dimension(:,:), allocatable, save :: part
      character(*) :: string,talk
!
      N_surface_reactions = count_reactions('mechanics.in')
      N_max_elements = count_max_elements('mechanics.in') + 1 
!
!  Allocate some arrays
!
      if (.not. allocated(part)) then
         allocate(part(N_max_elements,N_surface_reactions))
      else
      end if
      if (.not. allocated(reaction_direction)) then
         allocate(reaction_direction(N_surface_reactions))
      else
      end if
      if (.not. allocated(flags)) then
         allocate(flags(N_surface_reactions),STAT=stat)
      else
      end if
!      if (stat>0) call fatal_error('register_indep_pchem',&
!           'Could not allocate memory for flags')
!
!  Read the heterogeneous chemical kinetics mechanism
!
      call read_mechanics_file('mechanics.in',part,n_max_elements,&
          reaction_direction,talk)
!
!  Count heterogeneous species
!
      call count_species(part,species,reactants,products,ns,nr,np)
      call count_species_type(species(:ns),N_adsorbed_species,&
          N_surface_species,ns)
      call count_species_type(reactants,N_trash,N_surface_reactants,nr)
!
      if (trim(string)=='N_surface_species') variable=N_surface_species
      if (trim(string)=='N_adsorbed_species') variable=N_adsorbed_species
      if (trim(string)=='N_surface_reactants') variable=N_surface_reactants
!
    end subroutine get_pchem_info
!***********************************************************************
    subroutine get_R_c_hat(var,fp)
!
      real, dimension(mpar_loc), intent(out) :: var 
      real, dimension(mpar_loc,mpvar) :: fp
!
!  JONAS talk to nils how to implement things from equ.f90
!  JONAS implement mdot_ck
!
      R_c_hat = 0.0
      call get_St(St,fp)
      R_c_hat(:) = -sum(mdot_ck,DIM=2)/(St(:)*mol_mass_carbon)
!
      var = R_c_hat
!
    end subroutine get_R_c_hat
!***********************************************************************
    subroutine get_R_j_hat(var)
!
      real, dimension(mpar_loc,N_adsorbed_species), intent(out) :: var 
      real, dimension(mpar_loc,N_adsorbed_species) :: R_j_hat
      integer :: j,k
!
!  Calculation of R_j_hat according to eq.50 in 8th US combustion 
!  meeting, coal and  biomass combustion and gasification.
!      
      R_j_hat = 0.0
      do k=1,N_surface_reactions
         do j=1,N_adsorbed_species-1
            R_j_hat(:,j)=R_j_hat(:,j)+(mu_prime(j,k)-mu(j,k))*RR_hat(j,k)
         enddo
      enddo
      var = R_j_hat
!
    end subroutine get_R_j_hat
!***********************************************************************
    subroutine get_mod_surf_area(var,fp)
!
      real, dimension(mpar_loc), intent(out) :: var
      real, dimension(mpar_loc) :: mod_all,Sgc,mod_surf_area
      real, dimension(mpar_loc,mpvar) :: fp
      integer :: end
      
!
      call get_St(St,fp)
!
!  mod_all: Middle term in eq. 40 of 8th US combustion meeting, coal and
!  biomass combustion and gasification.
!
      mod_all(:) = St_init(:)*St_init(:)*struct_par* & 
           (1-conversion(:)) * (1-conversion(:)) / &
           (2*St(:)**2)
      Sgc(:) = St(:)/fp(:,irhopswarm)
!
      mod_surf_area(:) = (1-mod_all(:))*Sgc(:)*mol_mass_carbon
!
      var = mod_surf_area
!
    end subroutine get_mod_surf_area
!***********************************************************************
    subroutine get_St(var,fp)
!      
      real, dimension(mpar_loc),intent(out) :: var
      real, dimension(mpar_loc,mpvar) :: fp
      integer :: end
!
!
!  Evolution of the total surface area according to 8th US combustion
!  meeting, coal and biomass combustion and gasification eq. 19  
!   
      St = 0.0
      St(:)=(1-conversion(:))*St_init(:)* & 
           sqrt(1.0 - struct_par*log(fp(:,irhopswarm)/rho_p_init(:)))
      var = St
!
   end subroutine get_St
!***********************************************************************
integer function count_max_elements(inputfile)
!
  character(*) :: inputfile
  integer :: stat,i,k,j,maxelement
  character(150) :: line
  character :: tab = char(9)
  character :: spc = char(32)
!      
 open(30, file=inputfile,iostat=stat)
    if(stat==0) then
          maxelement=1
!       write(*,*) 'Counting elements'
       do
          read(30,'(A150)', end=530) line
          if (line(:1)/='!') then
          i=1
          k=0
          do while (i <len(line))
             if (line(i:i)/='+'.and. &
                  line(i:i)/=spc.and.line(i:i)/=tab) then
                j=1
                do while (line(i+j:i+j)/='+'.and.&
                 line(i+j:i+j)/=spc.and.line(i+j:i+j)/=tab)
                   j=j+1
                enddo
                k=k+1
                i=i+j
             else
                i=i+1
             endif
          if (k>maxelement) then
             maxelement=k
          else
          end if
          enddo
          else
          end if
       end do
530 close(30)
      count_max_elements=maxelement
    else 
       maxelement=0
       write(*,*) 'Problem with the single elements of the mechanics.in file'
    endif
!
end function count_max_elements
!**************************************************
  integer function count_reactions(inputfile)
!
!  this function counts the uncommented lines in the mechanics.in
!
    integer :: stat,reactions
    character(150) :: line
    character(*) :: inputfile
!    
 open(20, file=inputfile,iostat=stat)
    if(stat==0) then
!       write(*,*) 'Counting reactions'
       reactions = 0
       do
          read(20,fmt =610, end=520) line
          if (line(:1)/='!') then
             reactions = reactions + 1
             if (index(line,'<>') > 0) then
                reactions = reactions + 1
             else
             end if
          else
          end if
       end do
520    close(20)
       count_reactions=reactions
610 format(A150)
    else 
       count_reactions=0
       write(*,*) 'Could not open mechanics file'
    endif
!
  end function count_reactions
!*********************************************************************** 
integer function find_species(species,unique_species,nlist)
!
!  function to replace the imu/inuX variables
!
   implicit none
!
    integer :: i,nlist
    character(len=*) :: species
    character(len=*), dimension(:) :: unique_species
!   
    find_species = 0
!
    do i=1,nlist
       if (trim(species) == trim(unique_species(i))) then
          find_species = i
       else
       end if
    end do
!
  end function find_species
!**********************************************************************
  subroutine set_RR(part,RR_method)
!
    integer :: i,j,stat,RR
    real, dimension(:) :: RR_method
    character(10), dimension(:,:) :: part
    character(10) :: element
!    
    do i=1,size(part,2)
       do j=1,size(part,1)
          element = part(j,i)
          if (element(:2) == 'RR') then
             read(element(2:),'(I1.1)',iostat=stat) RR
             RR_method = RR
          else
          end if
       end do
    end do
!
  end subroutine set_RR
!**********************************************************************
  subroutine create_arh_param(part,B_k,ER_k,sigma_k)
!
!takes the first numerical in part and writes it to b_k
!
    character(10), dimension(:,:) :: part
    real, dimension(:) :: B_k,ER_k,sigma_k
    character(10) :: el_B_k,el_ER_k,el_sigma
    real :: B_k_single,ER_k_single,sigma_single
    logical :: done
    integer :: i,k,stat
!
    B_k = 0.0
    ER_k = 0.0
    sigma_k = 0.0
!    
    do i=1, size(part,2)
       if (part(size(part,1),i) == 'rev') then
          B_k(i) = 1e1
          ER_k(i) = 1.
          sigma_k(i) = 1e1
       else
       done = .false.
       do k=1, size(part,1)-2
          el_B_k = part(k,i)
          el_ER_k = part(k+1,i)
          el_sigma = part(k+2,i)
             read(el_B_k,*,iostat=stat) B_k_single
             if (stat == 0 .and. (done .eqv. .false.)) then
                B_k(i) = B_k_single
                done = .true.
                read(el_ER_k,*,iostat=stat) ER_k_single
                if (stat == 0) then
                   ER_k(i) = ER_k_single
                else
                end if
                read(el_sigma,*,iostat=stat) sigma_single
                if (stat == 0) then
                   sigma_k(i) = sigma_single
                else
                end if
             else
             end if
       end do
       end if
    end do
!
    ER_k = ER_k*1e6/gas_constant
!          
  end subroutine create_arh_param
!**********************************************************************
  subroutine create_dependency(nu,dependent_reactant,&
     n_surface_reactions,n_surface_reactants)
!   
    integer :: i,k,n_surface_reactions,n_surface_reactants
    real, dimension(:,:) :: nu
    integer, dimension(:) :: dependent_reactant
!   
    dependent_reactant = 0
!   
    do i=1,n_surface_reactants
       do k=1,n_surface_reactions
          if (nu(i,k) > 0) then
             dependent_reactant(k) = i
          else
          end if
       end do
    end do
!   
  end subroutine create_dependency
!**********************************************************************
  subroutine create_occupancy(adsorbed_species_names,site_occupancy)
!    
    integer :: i
    character(10), dimension(:) :: adsorbed_species_names
    real, dimension(:) :: site_occupancy
!
    site_occupancy = 1
    do i=1, size(site_occupancy,1)
       if  (index(adsorbed_species_names(i), '(O2)') > 0) then
          site_occupancy(i) = 2
       else
       end if
    end do
!
  end subroutine create_occupancy
!**********************************************************************
  subroutine create_stoc(part,list,targ,lhs,nlist)
!
  integer :: i,j,k,stat,nlist
  real :: multi
  character(10), dimension(:,:) :: part
  character(10), dimension(:) :: list
  character(10) :: element
  real, dimension(:,:) :: targ
  logical :: lhs,fwd
!
!  list where the stochiometry is saved in
!
  targ = 0
  do i=1,size(part,2)
     fwd = lhs
     do j=1,size(part,1)
        do k=1,nlist
           if (part(j,i) == '->' .or. &
                part(j,i) == '<>') then
              fwd =  .not. lhs
           else
           end if
           element = part(j,i)
!
!  check if character is numeric
!
           read(element(:1),*,iostat=stat) multi
           if (stat==0) then
              element = element(2:)
           else
              multi = 1.0
           end if
!
!  if string is numeric, change stochiometric factor accordingly
!
           if (element==list(k) .and. &
                fwd .eqv. .true.) then
              targ(k,i) =real(multi)
           else
           end if
        end do
      end do
   end do
   targ(:,:) = int(targ(:,:))
!
  end subroutine create_stoc
!**********************************************************************
  subroutine get_ac(ac,list,nlist)
!
!  gets how many c atoms are on the surface species
!
    integer :: i,c_place,stat,nc,nlist
    character(len=10) :: species_in_q
    logical :: numeric
    real, dimension(:) :: ac
    character(10), dimension(:) :: list
    ac = 0
!
    do i = 1,nlist
       if (scan(list(i),'C') > 0) then
          c_place = scan(list(i),'C')
          species_in_q = list(i)
          read(species_in_q(c_place+1:c_place+1),120,iostat=stat) nc
120 format (I1.1)
          numeric = (stat == 0)
          if (numeric) then
          ac(i) = nc
          else
          ac(i) = 1
          end if
        else
        end if
    end do
!    
  end subroutine get_ac
!**********************************************************************
 subroutine sort_compounds(lhslist,species_list,nlist,n_big)
!
!  this file reads in the order of the 
!  compounds as prescribed
!
   integer :: nlist,i,n_big,j
   character(10), dimension(:) :: lhslist
   integer :: front, end
   character(10), dimension(:) :: species_list
   character(10), dimension(:), allocatable :: temp_list
   character(10) :: temp
   logical :: is_reactant
!
   allocate(temp_list(nlist))
!
   end = 0
   front = 0
!
   do i=1,nlist
      is_reactant = .false.
      do j=1,n_big
         if (species_list(i) == lhslist(j)) then
            is_reactant = .true.
         else
         end if
      end do
!
      if (.not. is_reactant) then
         temp_list(nlist-end) = species_list(i)
         end = end + 1
      else 
         temp_list(1+front) = species_list(i)
         front = front + 1
      end if
   end do
!   
   do i=1,nlist-1
      if (temp_list(i) == 'Cf') then
         temp = temp_list(nlist)
         temp_list(nlist) = 'Cf'
         temp_list(i)=temp
      else
      end if
   end do
!
   species_list(:nlist) = temp_list(:nlist)
!
  end subroutine sort_compounds
!**********************************************************************
  subroutine create_ad_sol_lists(list,target_list,ad_sol,nlist)
!
!  create lists of adsorbed and solid species
!
  character(10), dimension(:) :: list,target_list
  character(*) :: ad_sol
  integer :: i ,nlist
  integer :: place
  place = 1
!    
    do i = 1,nlist
       if (ad_sol == 'ad') then
          if (scan(list(i),'()') > 0 .or.&
               list(i) == 'Cf') then
             target_list(place) = list(i)
             place = place + 1
          else
          end if
       else
       end if
       if (ad_sol == 'sol') then
          if (scan(list(i),'()') == 0 .and.&
               list(i)/='Cb' .and.&
               list(i)/='Cf') then
             target_list(place) = list(i)
             place = place + 1
          else
          end if
       else
       end if
    end do
!
  end subroutine create_ad_sol_lists
!**********************************************************************
  subroutine count_species_type(list,n_ad,n_sol,nlist)
!
    integer :: i,parenthes,nlist
    integer :: n_ad, n_sol
    character(10), dimension(:) :: list
!
!  count adsorbed and surface species
!
    n_ad = 0
    n_sol = 0
    do i = 1,nlist
       parenthes = 0
       parenthes = scan(list(i),'()')
       if (parenthes > 0 .or. &
            list(i) == 'Cf') then
          n_ad = n_ad + 1
       else
          if (list(i)/='Cb') then
             n_sol = n_sol + 1
          else
          end if
       end if
    end do
   end subroutine count_species_type
!**********************************************************************
  logical function is_not_in_list(list, entry)
!
  integer :: k
  character(10), dimension(:) :: list
  character(10) :: entry
!  
  is_not_in_list = .true.
!  
  do k=1,size(list,1)
     if (entry == list(k)) then
        is_not_in_list = .false.
     else
     end if
  end do
!
  end function is_not_in_list
!**********************************************************************
  subroutine count_species(part,species,reactants,products,ns,nr,np)
!
    character(10), dimension(:,:) :: part
    character(10) :: element
    real :: numeric
    integer :: i,j,jmax,stat,place,ns,nr,np,number
    integer :: place_reac, place_prod
    logical :: lhs, to_append, to_append_prod,to_append_reac
    character(10), dimension(40) :: species,reactants,products
    character(10), dimension(40) :: temp,temp_reac,temp_prod
!  
    temp = 'nothing'
    jmax = size(part,1)
    place = 1
    place_reac = 1
    place_prod = 1

    do i=1,n_surface_reactions
       lhs = .true.
       do j=1,jmax
          element = part(j,i)
!
!  switch when the reaction arrow is read
!
          if (element == '->' .or. &
               element == '<>') then
             lhs = .false.
          else
          end if
!
!  if element can be read as real, disregard
!
          read(element,*,iostat=stat) numeric
          if (stat /= 0 .and. &
               element /= '->' .and. &
               element /= '<>' .and. &
               element(:2) /= 'RR') then
             read(element(:1),*,iostat=stat) number
             if (stat==0)  element = element(2:)
!
!  appending the components to the list 
!  of global, reactand and product uniques          
!
          to_append = is_not_in_list(temp,element)
          to_append_reac = is_not_in_list(temp_reac,element)
          to_append_prod = is_not_in_list(temp_prod,element)
!
          if (to_append) then
          temp(place) = element
          place = place+1
          else
          end if  
!
          if (to_append_reac .and. lhs) then
          temp_reac(place_reac) = element
          place_reac = place_reac+1
          else
          end if  
!
          if (to_append_prod .and. (lhs .eqv. .false.)) then
          temp_prod(place_prod) = element
          place_prod = place_prod+1
          else
          end if  
!
          else
          end if
       end do
    end do
!
!creating the lists
!
    species(:place-1)=temp(:place-1)
    reactants(:place_reac-1)=temp_reac(:place_reac-1)
    products(:place_prod-1)=temp_prod(:place_prod-1)
    ns = place-1
    nr = place_reac-1
    np = place_prod-1
!
    if(ldebug) then
       print*,'Number of unique species:', ns
       print*,species(:ns)
       print*,'Number of unique reactands:', nr
       print*,reactants(:nr)
       print*,'Number of unique products:', np
       print*,products(:np)
    else
    endif
!
  end subroutine count_species
!**********************************************************************
subroutine flip_and_parse(string,ireaction,target_list,direction)
!
    character(150) :: string,flipped_string
    character(50) :: lhs,sign,rhs,end
    character(3), dimension(:) :: direction
    integer :: ireaction
    integer :: i,numerical,marker
    real :: real_number 
    character(10), dimension(:,:) :: target_list
!    
    marker = index(string,'<>')
    numerical = 1
    i = marker
!
    do while (numerical  /= 0 )
       i = i + 1
       read(string(i:i+7),*,iostat=numerical) real_number
       if (real_number < 10) then
          numerical = 1
       else
       end if
       if (i > len(string)-10) then
          print*,'no numericals found after sign!'
          numerical = 0
       else
       end if
    end do
!    
    lhs = trim(string(:marker-1))
    sign = trim(string(marker:marker+1))
    rhs = trim(string(marker+2:i))
    end = trim(string(i:))
!
    flipped_string = trim(rhs)//'  '//trim(sign)//'  '//trim(lhs)// trim(end)
    flags(ireaction) = 'rev'
    call parse(flipped_string,ireaction,target_list,'rev',direction)
!
  end subroutine flip_and_parse
!**********************************************************************
  subroutine parse(string,ireaction,target_list,flag,direction)
!
    character(150) :: string
    character :: tab = char(9)
    character :: spc = char(32)
    character(3) :: flag
    character(3), dimension(:) :: direction
    integer :: i,j,k,ireaction
    character(10) :: formatting
    character(10), dimension(:,:) :: target_list
!
    i=1
    k=1
    j=1
!
    do while (i<=len(string))
       if (string(i:i)==tab) then
          string(i:i) = spc
       else
       i=i+1
       endif
    end do
    string = trim(string)
    i=1
    do while (i <len(string))
       if (string(i:i)/='+'.and. &
           string(i:i)/=spc) then
           j=1
          do while (string(i+j:i+j)/='+'.and.&
                 string(i+j:i+j)/=spc)
             j=j+1
          enddo
          formatting = adjustl(string(i:i+j))
          target_list(k,ireaction)=formatting
          k=k+1
          i=i+j
       else
       i=i+1
       endif
   enddo
   if (i==len(string)) then
      target_list(k-1:,ireaction) = '0.0'
   else
   end if
!
   direction(ireaction) = flag
!
  end subroutine parse
!**********************************************************************
  subroutine read_mechanics_file(inputfile,target_list,n_max_elements,&
    reaction_direction,talk)
!  
    integer :: stat,ireaction,i,n_max_elements
    character(150) :: string
    character(*) :: inputfile,talk
    character(10) :: writeformat
    character(3), dimension(:) :: reaction_direction
    character(10), dimension(:,:) :: target_list
!
    writeformat = '(  A10,A3)'
    write(writeformat(2:3),'(I2)') n_max_elements
    open(20, file=inputfile,iostat=stat)
    open(29, file='mech_outputfile.dat',iostat=stat)
!    
    if(stat==0) then
       if(talk=='verbose') write(*,*) 'Opened mechanics file'
       ireaction = 1
       do
          read(20,fmt =510, end=500) string
          if ( (string(:1)) /='!') then
             flags(ireaction) = 'fwd'
          call parse(string,ireaction,target_list,'fwd',reaction_direction)
          ireaction = ireaction + 1
          if (index(string,'<>') > 0) then
            call flip_and_parse(string,ireaction,target_list,reaction_direction)
            ireaction = ireaction + 1
          else
          end if
          else
          end if
       end do
500 if(talk=='verbose') print*,'Done parsing mechanics file'
       close(20)
       if(talk=='verbose') then
          do i=1,N_surface_reactions
             write(*,writeformat) target_list(:,i),reaction_direction(i)
             write(29,writeformat) target_list(:,i),reaction_direction(i)
          enddo
       else
       end if
       close(29)
510 format (A150)
    else 
       write(*,*) 'Could not open mechanics file'
    endif
!
  end subroutine read_mechanics_file
!**********************************************************************
  subroutine get_species_list(string,list)
!  
    character(*) :: string
    character(10), dimension(40) :: list,species,reactants
    integer :: nr,ns,stat,dummy,N_adsorbed_species
    integer :: N_surface_species
!
    call get_pchem_info(species,'N_adsorbed_species',N_adsorbed_species,'quiet')
    call get_pchem_info(species,'N_surface_species',N_surface_species,'quiet')

    call create_ad_sol_lists(species,adsorbed_species_names,'ad',ns)
    call sort_compounds(reactants,adsorbed_species_names,N_adsorbed_species,nr)
!
    allocate(solid_species(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for solid_species')
!
    call create_ad_sol_lists(species(:ns),solid_species,'sol',ns)
    call sort_compounds(reactants,solid_species,N_surface_species,nr)
!
    if (trim(string)=='solid_species') & 
         list(1:N_surface_species)=solid_species(:)
    if (trim(string)=='adsorbed_species_names') &
         list(1:N_adsorbed_species)=adsorbed_species_names(:N_adsorbed_species)
!
    deallocate(solid_species)
!
  end subroutine get_species_list
!**********************************************************************
  subroutine create_dngas(nu,nu_prime,dngas)
!
!  Find the mole production of the forward reaction. This will later
!  be used for the calculation of the reverse reaction rate.
!
  real, dimension(:) :: dngas
  real, dimension(:,:) :: nu,nu_prime
  integer :: k
!
  do k=1,N_surface_reactions
    dngas(k) = sum(nu_prime(:,k)) - sum(nu(:,k))
  end do
!
end subroutine create_dngas
!**********************************************************************
  subroutine find_entropy_of_reaction(fp)
!  
    integer :: k,j,i
    real, dimension(mpar_loc,mpvar) :: fp
!
! JONAS: units are in j/(kmol*K)
!
    if (inuH2O>0) surface_species_entropy(:,inuH2O)= 189.00e3+(0.0425e3*fp(:,iTp))
    if (inuO2>0)  surface_species_entropy(:,inuO2) = 222.55e3+(0.0219e3*fp(:,iTp))
    if (inuCO2>0) surface_species_entropy(:,inuCO2)= 212.19e3+(0.0556e3*fp(:,iTp))
    if (inuH2>0)  surface_species_entropy(:,inuH2 )= 133.80e3+(0.0319e3*fp(:,iTp))
    if (inuCO>0)  surface_species_entropy(:,inuCO )= 199.35e3+(0.0342e3*fp(:,iTp))
!
!  taken from chemistry  webbook (1bar)
!
    if (inuCH>0)  surface_species_entropy(:,inuCH )=183.00e3
    if (inuHCO>0) surface_species_entropy(:,inuHCO)=223.114e3+(0.0491e3*fp(:,iTp))
    if (inuCH2>0) surface_species_entropy(:,inuCH2)=193.297e3+(0.0467e3*fp(:,iTp))
!
!  taken from chemistry webbook (1bar)
!
    if (inuCH4>0) surface_species_entropy(:,inuCH4)= 189.00e3
    if (inuCH3>0) surface_species_entropy(:,inuCH3)= 190.18e3+(0.0601e3*fp(:,iTp))
!
    if (imuadsO>0)    then
       adsorbed_species_entropy(:,imuadsO) = &
    (164.19e3+(0.0218e3*fp(:,iTp)))*0.72 - (3.3*gas_constant)
    else
    end if
!
!  this is guessed
!
    if (imuadsO2>0)   then
       adsorbed_species_entropy(:,imuadsO2) =  &
            2*adsorbed_species_entropy(:,imuadsO)
    else
    end if
    if (imuadsOH>0)   then
       adsorbed_species_entropy(:,imuadsOH) = &
         ((0.0319e3*fp(:,iTp)) + 186.88e3) * 0.7 - (3.3*gas_constant)
    else
    end if
    if (imuadsH>0)    then 
       adsorbed_species_entropy(:,imuadsH) = &
           (117.49e3+(0.0217e3*fp(:,iTp)))*0.54 - (3.3*gas_constant)
    else
    end if
    if (imuadsCO>0)   then
       adsorbed_species_entropy(:,imuadsCO) = &
            surface_species_entropy(:,inuCO)* &
            0.6*(1+(1.44e-4*fp(:,iTp))) - (3.3*gas_constant)
    else
    end if
!
!  taken from nist
!
    if (imufree>0)    adsorbed_species_entropy(:,imufree) = 0

      entropy_k=0
!
    do k=1,N_surface_reactions
      do i=1,N_surface_species
        entropy_k(:,k)=entropy_k(:,k)&
            +nu_prime(i,k)*surface_species_entropy(:,i)&
            -nu(i,k)*surface_species_entropy(:,i)
      enddo
      if (N_adsorbed_species > 0) then
      do j=1,N_adsorbed_species
        entropy_k(:,k)=entropy_k(:,k)&
            +mu_prime(j,k)*adsorbed_species_entropy(:,j)&
            -mu(j,k)*adsorbed_species_entropy(:,j)
      enddo
      else
      end if
    enddo      
!
  end subroutine find_entropy_of_reaction
!**********************************************************************
  subroutine get_reverse_K_k(k,K_k,fp)
!  
    use Particles_cdata, only: interp_pp
!
  integer :: k
  real, dimension(mpar_loc) :: k_p,k_c
  real, dimension(mpar_loc) :: denominator, exponent
  real, dimension(mpar_loc,mpvar) :: fp
  real, dimension(:,:) :: K_k
!
  denominator(:) = heating_k(:,k-1) - (entropy_k(:,k-1)*fp(:,iTp))
  exponent(:) = denominator(:)/(gas_constant*fp(:,iTp))
  k_p(:) = exp(-exponent(:))
  k_c(:) = k_p(:) / (((gas_constant)*fp(:,iTp)/interp_pp)**(dngas(k-1)))
  K_k(:,k) = (K_k(:,k-1) / k_c(:))
!
end subroutine get_reverse_K_k
!**********************************************************************
  subroutine save_current_conversion(fp)
!
    real, dimension(mpar_loc,mpvar) :: fp
!
    conversion(:) = 4/3*pi*fp(:,irhop)*fp(:iap)*fp(:iap)*fp(:iap) &
         / initial_mass(:)
!
  end subroutine save_current_conversion
!**********************************************************************
  end module Particles_chemistry
