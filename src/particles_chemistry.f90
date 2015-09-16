! $Id: Particles_chemistry.f90 21950 2014-07-08 08:53:00Z jonas.kruger $
!
!  This module takes care of everything related to reactive particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_chemistry=.true.
!
!***************************************************************
!
!  The assumptions and equations implemented in this routine are
!  deduced partially from following papers:
!  8th US Combustion Meeting – Paper # 070CO-0312
!  Transient simulations of char gasification     by
!  Nils Erland L. Haugen
!  Reginald E. Mitchell
!  Matt Tilghman
!
!***************************************************************
!
module Particles_chemistry
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_mpicomm
  use Particles_sub
  use EquationOfState
  use Chemistry
  use SharedVariables, only: put_shared_variable
!
  implicit none
!
  include 'particles_chemistry.h'
!
!***************************************************************!
!  Particle independent variables below here                    !
!***************************************************************!
!
  real, dimension(:), allocatable :: reaction_order
  real, dimension(:), allocatable :: effectiveness_factor_old
  real, dimension(:), allocatable :: B_k, ER_k, ac
  real, dimension(:), allocatable :: sigma_k
  real, dimension(:), allocatable :: omega_pg_dbl
  real, dimension(:,:), allocatable :: mu, mu_prime
  real, dimension(:,:), allocatable :: nu, nu_prime
  real, dimension(:), allocatable :: aac, T_k
  real, dimension(:), allocatable :: dngas
  integer, dimension(:), allocatable :: dependent_reactant
  real, dimension(:,:), allocatable, save :: part_power
!
  real, dimension(nchemspec) :: molar_mass=1.0
  real, dimension(50) :: reaction_enhancement=1
  character(len=3), dimension(:), allocatable, save :: reaction_direction
  character(len=3), dimension(:), allocatable ::flags
  character(len=10), dimension(:,:), allocatable, save :: part
  character(len=10), dimension(50) :: species_name
  character(len=10), dimension(40) :: reactants, products
  character(len=20) :: element, writeformat
  character(len=10), dimension(40) :: species
!
  integer :: N_species, N_reactions
  integer :: placeholder=1
  integer :: N_surface_reactions, N_adsorbed_species
  integer :: N_surface_reactants, N_surface_species
  integer :: nr, ns, np, N_max_elements
  integer :: inuH2, inuCO2, inuH2O, inuCO, inuCH4, inuO2
  integer :: inuCH, inuHCO, inuCH2, inuCH3
  integer :: imufree, imuadsO, imuadsO2, imuadsOH, imuadsH, imuadsCO
  integer :: iter=0
  integer, dimension(:), allocatable :: jmap
  real :: x_s_total
  real :: effectiveness_factor_timeaver=1.
  real :: eta_int=0., delta_rho_surface=0.
  real :: A_p_first, rho_p_first
  real :: diffusivity = 0.0
  real, target :: total_carbon_sites=1.08e-8 ! [mol/cm^2]
  real :: chemplaceholder=0.0
  real :: tortuosity=3.
  real, target :: true_density_carbon=1.800 ! g/cm^3
  real :: structural_parameter=8. ! [-]
  real :: startup_time=0.
  real :: startup_quench
  real :: pre_energy=1.0
!
  logical :: first_pchem=.true.
  logical :: lpchem_debug = .false.
  logical :: lthiele=.false.
  logical :: lreactive_heating=.false.
!
!*********************************************************************!
!             Particle dependent variables below here                 !
!*********************************************************************!
!
  real, dimension(:), allocatable :: conversion
  real, dimension(:,:), allocatable :: mdot_ck, RR_hat
  real, dimension(:), allocatable :: St, rho_p
  real, dimension(:), allocatable :: R_c_hat
  real, dimension(:,:), allocatable :: heating_k, entropy_k
  real, dimension(:), allocatable :: Particle_temperature
  real, dimension(:), allocatable :: effectiveness_factor
  real, dimension(:,:), allocatable :: effectiveness_factor_species
  real, dimension(:,:), allocatable :: surface_species_enthalpy
  real, dimension(:,:), allocatable :: surface_species_entropy
  real, dimension(:,:), allocatable :: adsorbed_species_enthalpy
  real, dimension(:,:), allocatable :: adsorbed_species_entropy
  real, dimension(:), allocatable :: ndot_total
  real, dimension(:,:), allocatable :: Rck, Rck_max
  real, dimension(:), allocatable :: f_RPM
  real, dimension(:,:), allocatable :: effectiveness_factor_reaction
  real, dimension(:,:), allocatable :: thiele
  real, dimension(:,:), allocatable :: ndot, K_k
  real, dimension(:,:), allocatable :: R_j_hat
  real, dimension(:,:), allocatable :: Cs
  real, dimension(:), allocatable :: initial_density
  real, dimension(:), allocatable :: Cg, A_p
  real, dimension(:), allocatable :: q_reac
  real, dimension(:), allocatable :: Nu_p
  real, dimension(:), allocatable :: mass_loss
!
  real :: mol_mass_carbon=12.0 !g/mol
  real :: Sgc_init=3e6 ! cm^2/g
!
!  is already in the code (R_CGS), with ergs as unit!!!
!  this one is used for
!
  real :: gas_constant=8.3144727 ![J/mol/K] MUST ALWAYS BE IN SI UNITS!!!
!
  namelist /particles_chem_init_pars/ &
      reaction_enhancement, &
      total_carbon_sites, &
      diffusivity, &
      tortuosity, &
      structural_parameter, &
      Sgc_init, &
      lpchem_debug, &
      true_density_carbon, &
      startup_time
!
  namelist /particles_chem_run_pars/ chemplaceholder, lthiele, lreactive_heating
!
  contains
! ******************************************************************************
!  Wrapper routine for particle dependent and independent chemistry
!  variables
!
!  oct-14/Jonas: coded
!
    subroutine register_particles_chem()
      integer :: ierr
      if (lroot) call svn_id( &
          "$Id: particles_chemistry.f90 20843 2014-10-06 18:45:43Z jonas.kruger $")
!
      call register_unit_system()
!
      call get_pchem_info(species,'N_species',N_species,'quiet')
!    print*, 'Number of species in mechanics file: ', N_species
!    print*, 'Species found: ', species(:N_species)
      call register_indep_pchem()
      call register_dep_pchem()
!
      call put_shared_variable('total_carbon_sites',total_carbon_sites,ierr)
      if (ierr /= 0) call fatal_error('register_particles_chem', 'unable to share total_carbon')
      call put_shared_variable('true_density_carbon',true_density_carbon,ierr)
      if (ierr /= 0) call fatal_error('register_particles_chem', 'unable to share true_density')
!
    endsubroutine register_particles_chem
! ******************************************************************************
!  Allocation of variables that are independent of the local particle number
!
!  oct-14/Jonas: coded
!
    subroutine register_indep_pchem()
      integer :: stat, i
      logical :: lenhance
!
      allocate(dngas(N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_psurfchem', &
          'Could not allocate memory for dngas')
      allocate(omega_pg_dbl(N_species),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_pchem', &
          'Could not allocate memory for omega_pg_dbl')
      allocate(B_k(N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_pchem', &
          'Could not allocate memory for B_k')
      allocate(Er_k(N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_pchem', &
          'Could not allocate memory for Er_k')
      allocate(sigma_k(N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_pchem', &
          'Could not allocate memory for sigma_k')
      allocate(reaction_order(N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_psurfchem', &
          'Could not allocate memory for reaction_order')
      allocate(effectiveness_factor_old(N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_pchem', &
          'Could not allocate memory for effectiveness_factor_old')
!
      effectiveness_factor_old = 1.
!
! Check if any of the reactions are enhanced
      lenhance = .false.
      do i = 1,N_surface_reactions
        if (reaction_enhancement(i)  /=  1) then
          print*,'**************** WARNING! ****************************'
          write (*,'(A5,I2,A25,F10.2)') &
              'Reac ',i,' is enhanced by a factor ',reaction_enhancement(i)
          lenhance = .true.
        endif
      enddo
!if (lenhance) call sleep(4)
!
! Define the Arrhenius coefficients. The term ER_k is the activation energy
! divided by the gas constant (R)
!
      call create_arh_param(part,B_k,ER_k,sigma_k)
    endsubroutine register_indep_pchem
! ******************************************************************************
!  Allocate memory for chemical variables that change from particle to particle
!
!  oct-14/Jonas: coded
!
    subroutine register_dep_pchem()
      integer :: stat
!
      allocate(R_c_hat(mpar_loc),STAT=stat)
      if (stat > 0) call fatal_error('register_dep_pchem', &
          'Could not allocate memory for R_c_hat')
      allocate(heating_k(mpar_loc,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_dep_pchem', &
          'Could not allocate memory for heating_k')
      allocate(entropy_k(mpar_loc,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_dep_pchem', &
          'Could not allocate memory for heating_k')
    endsubroutine register_dep_pchem
!***********************************************************************
    subroutine pencil_criteria_par_chem()
!
!  All pencils that the Particles_chemistry module depends on are specified here.
!
!  16.09.2015/jonas + nils: coded
!
      if (lthiele) then
        lpenc_requested(i_Diff_penc_add)=.true.
      endif
!
    endsubroutine pencil_criteria_par_chem
!***********************************************************************
    subroutine calc_pencils_par_chem(f,p)
!
!  Calculate Particles pencils.
!  Most basic pencils should come first, as others may depend on them.
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
!  Set up required variables etc. for reactive particles
!
!  02-sep-14/nils: coded
!
!  Count number of heterogeneous reactions and elements
!
    subroutine get_pchem_info(species,string,variable,talk)
      integer :: N_trash, variable, stat
      character(len=10), dimension(:) :: species
      character(len=*) :: string, talk
!
      N_surface_reactions = count_reactions('mechanics.in')
      N_max_elements = count_max_elements('mechanics.in') + 1
!
! Allocate some arrays
      if (.not. allocated(part)) then
        allocate(part(N_max_elements,N_surface_reactions))
      endif
      if (.not. allocated(part_power)) then
        allocate(part_power(N_max_elements,N_surface_reactions))
      endif
      if (.not. allocated(reaction_direction)) then
        allocate(reaction_direction(N_surface_reactions))
      endif
      if (.not. allocated(flags)) then
        allocate(flags(N_surface_reactions),STAT=stat)
      endif
      if (.not. allocated(T_k)) then
        allocate(T_k(N_surface_reactions),STAT=stat)
      endif
!      if (stat>0) call fatal_error('register_indep_pchem',&
!           'Could not allocate memory for flags')
!
!  Read the heterogeneous chemical kinetics mechanism
!
      call read_mechanics_file('mechanics.in',part,n_max_elements, &
          reaction_direction,talk)
!
! Count heterogeneous species
      call count_species(part,species,reactants,products)
      call count_species_type(species(:ns),N_adsorbed_species, &
          N_surface_species,ns)
      call count_species_type(reactants,N_trash,N_surface_reactants,nr)
!
      if (trim(string) == 'N_surface_species') variable = N_surface_species
      if (trim(string) == 'N_adsorbed_species') variable = N_adsorbed_species
      if (trim(string) == 'N_surface_reactants') variable = N_surface_reactants
      if (trim(string) == 'N_species') variable = ns
      if (trim(string) == 'N_reactants') variable = nr
    endsubroutine get_pchem_info
! ******************************************************************************
!  Calculate the the molar flux of carbon from the particle surface
!
!  oct-14/Jonas: coded
!
    subroutine calc_R_c_hat()
!
      integer :: k,k1,k2      
!      
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      R_c_hat = 0.0
!
      do k =k1,k2
        R_c_hat(k) = -sum(mdot_ck(k,:))/(St(k)*mol_mass_carbon)
      enddo
!
    endsubroutine calc_R_c_hat
! ******************************************************************************
!  Calculation of the modified surface area of a particle due to changes
!  in total site concentration and conversion
!
!  oct-14/Jonas: coded
!
    subroutine calc_get_mod_surf_area(mod_surf_area,fp)
      real, dimension(:), allocatable :: mod_all, Sgc
      real, dimension(k1_imn(imn):k2_imn(imn)) :: mod_surf_area
      real, dimension(:,:) :: fp
      integer :: k, k1, k2
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      allocate(mod_all(k1:k2))
      allocate(Sgc(k1:k2))
!
!  mod_all: Middle term in eq. 40 of 8th US combustion meeting, coal and
!  biomass combustion and gasification.
!
      do k = k1,k2
        mod_all(k) = (Sgc_init*fp(k,impinit))**2*structural_parameter* &
            (1-conversion(k)) * (1-conversion(k)) / &
            (2*St(k)**2)
        Sgc(k) = St(k)/ rho_p(k)
        mod_surf_area(k) = (1-mod_all(k))*Sgc(k)*mol_mass_carbon
      enddo
!
      deallocate(mod_all)
      deallocate(Sgc)
    endsubroutine calc_get_mod_surf_area
! ******************************************************************************
!  Evolution of the total surface area of each coal particle
!
!  oct-14/Jonas: coded
!
    subroutine calc_St(fp)
      real, dimension(:,:) :: fp
      real :: rho_p_init
      integer :: k, k1, k2
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
!  Evolution of the total surface area according to 8th US combustion
!  meeting, coal and biomass combustion and gasification eq. 19
!
      St = 0.0
!
      do k = k1,k2
        rho_p_init = fp(k,impinit)/(4.*pi/3.*fp(k,iapinit)**3)
        St(k) = fp(k,imp)*Sgc_init* &
            sqrt(1.0 - structural_parameter*log(rho_p(k)/rho_p_init))
      enddo
    endsubroutine calc_St
! ******************************************************************************
!  Reading in of the Arrhenius parameters from the parsed mechanics.in file
!
!  oct-14/Jonas: coded
!
!takes the first numerical in part and writes it to b_k
!
    subroutine create_arh_param(part,B_k,ER_k,sigma_k)
      character(len=10), dimension(:,:) :: part
      real, dimension(:) :: B_k, ER_k, sigma_k
      character(len=10) :: el_B_k, el_ER_k, el_sigma
      real :: B_k_single, ER_k_single, sigma_single
      logical :: done
      integer :: i, k, stat, stat1, stat2
!
      B_k = 0.0
      ER_k = 0.0
      sigma_k = 0.0
!
      do i = 1, size(part,2)
        if (part(size(part,1),i) == 'rev') then
          B_k(i) = 1e1
          ER_k(i) = 1.
          sigma_k(i) = 1e1
        else
          done = .false.
          do k = 1, size(part,1)-2
            el_B_k = part(k,i)
            el_ER_k = part(k+1,i)
            el_sigma = part(k+2,i)
            read (el_B_k,*,iostat=stat) B_k_single
            if (stat == 0 .and. (done .eqv. .false.)) then
              B_k(i) = B_k_single
              done = .true.
              read (el_ER_k,*,iostat=stat) ER_k_single
              if (stat == 0) then
                ER_k(i) = ER_k_single
              endif
              read (el_sigma,*,iostat=stat) sigma_single
              if (stat == 0) then
                sigma_k(i) = sigma_single
              endif
            endif
          enddo
        endif
      enddo
!
      ER_k = ER_k/gas_constant
      print*, 'ER_k'
    endsubroutine create_arh_param
! ******************************************************************************
    subroutine create_dependency(nu,dependent_reactant,n_surface_reactions, &
        n_surface_reactants)
!
!  Creation of a list that contains which reaction has which gas phase species
!  in its left hand side
!
!  oct-14/Jonas: coded
!
      integer :: i, k, n_surface_reactions, n_surface_reactants
      real, dimension(:,:) :: nu
      integer, dimension(:) :: dependent_reactant
!
      dependent_reactant = 0
!
      do i = 1,n_surface_reactants
        do k = 1,n_surface_reactions
          if (nu(i,k) > 0) then
            dependent_reactant(k) = i
          endif
        enddo
      enddo
    endsubroutine create_dependency
! ******************************************************************************
!  Creation of a list of how many C-atoms an adsorbed species contains
!  (right now: 1 for all but C2O2)
!
!  oct-14/Jonas: coded
!
    subroutine create_occupancy(adsorbed_species_names,site_occupancy)
      integer :: i
      character(len=10), dimension(:) :: adsorbed_species_names
      real, dimension(:) :: site_occupancy
!
      site_occupancy = 1
      do i = 1, size(site_occupancy,1)
        if  (index(adsorbed_species_names(i), '(O2)') > 0) then
          site_occupancy(i) = 2
        endif
      enddo
    endsubroutine create_occupancy
! ******************************************************************************
!  Creation of mu, nu and their counterparts, mu_prime and nu_prime
!  nu is what gas phase species are used in the reaction, mu, which adsorbed species
!  nu_prime is what gas phase species are produced, mu_prime which adsorbed species
!  the power term is written if the reaction has a arbitrary order for one or more
!  species concentrations
!
!  oct-14/Jonas: coded
!
    subroutine create_stoc(list,targ,lhs,nlist,power)
      integer :: i, j, k, stat, nlist
      real :: multi
!    character(len=10), dimension(:,:) :: part
      character(len=10), dimension(:) :: list
      character(len=10) :: element
      real, dimension(:,:) :: targ
      real, dimension(:,:), optional :: power
      logical :: lhs, fwd, forpower
!
! list where the stochiometry is saved in
      if (lhs) then
        forpower = .true.
      else
        forpower = .false.
      endif
      if (present(power)) power = 0.0
      targ = 0
      do i = 1,size(part,2)
        fwd = lhs
        do j = 1,size(part,1)
          do k = 1,nlist
            if (part(j,i) == '->' .or. part(j,i) == '<>') then
              fwd =  .not. lhs
            endif
            element = part(j,i)
!
! check if character is numeric
            read (element(:1),*,iostat=stat) multi
            if (stat == 0) then
              element = element(2:)
            else
              multi = 1.0
            endif
!
! if string is numeric, change stochiometric factor accordingly
            if (element == list(k) .and. fwd .eqv. .true.) then
              targ(k,i) = real(multi)
              if (forpower) then
                if (present(power)) power(k,i) = part_power(j,i)
              endif
            endif
          enddo
        enddo
      enddo
!
    endsubroutine create_stoc
! ******************************************************************************
!  Gets how many c atoms are on the surface species
!
!  oct-14/Jonas: coded
!
    subroutine get_ac(ac,list,nlist)
      integer :: i, c_place, stat, nc, nlist
      character(len=10) :: species_in_q
      logical :: numeric
      real, dimension(:) :: ac
      character(len=10), dimension(:) :: list
      ac = 0
!
      do i = 1,nlist
        if (scan(list(i),'C') > 0) then
          c_place = scan(list(i),'C')
          species_in_q = list(i)
          read (species_in_q(c_place+1:c_place+1),'(I1.1)',iostat=stat) nc
          numeric = (stat == 0)
          if (numeric) then
            ac(i) = nc
          else
            ac(i) = 1
          endif
        endif
      enddo
    endsubroutine get_ac
! ******************************************************************************
!  This script reads in the order of the compounds as prescribed
!
!  oct-14/Jonas: coded
!
    subroutine sort_compounds(lhslist,species_list,nlist)
      integer :: nlist, i,j
      character(len=10), dimension(:) :: lhslist
      integer :: front, ende
      character(len=10), dimension(:) :: species_list
      character(len=10), dimension(:), allocatable :: temp_list
      character(len=10) :: temp
      logical :: is_reactant
!
      allocate(temp_list(nlist))
!
      ende = 0
      front = 0
!
      do i = 1,nlist
        is_reactant = .false.
        do j = 1,nr
          if (species_list(i) == lhslist(j)) then
            is_reactant = .true.
          endif
        enddo
!
        if (.not. is_reactant) then
          temp_list(nlist-ende) = species_list(i)
          ende = ende + 1
        else
          temp_list(1+front) = species_list(i)
          front = front + 1
        endif
      enddo
!
      do i = 1,nlist-1
        if (temp_list(i) == 'Cf') then
          temp = temp_list(nlist)
          temp_list(nlist) = 'Cf'
          temp_list(i) = temp
        endif
      enddo
!
      species_list(:nlist) = temp_list(:nlist)
!      print*,'species_list=',species_list
    endsubroutine sort_compounds
! ******************************************************************************
!  Create lists of adsorbed and solid species
!
!  oct-14/Jonas: coded
!
    subroutine create_ad_sol_lists(target_list,ad_sol)
      character(len=10), dimension(:) :: target_list
      character(len=*) :: ad_sol
      integer :: i
      integer :: place
      place = 1
!
      do i = 1,ns
        if (ad_sol == 'ad') then
          if (scan(species(i),'()') > 0 .or. species(i) == 'Cf') then
            target_list(place) = species(i)
            place = place + 1
          endif
        endif
        if (ad_sol == 'sol') then
          if (scan(species(i),'()') == 0 .and. &
              species(i) /= 'Cb' .and. &
              species(i) /= 'Cf') then
            target_list(place) = species(i)
            place = place + 1
          endif
        endif
      enddo
    endsubroutine create_ad_sol_lists
! ******************************************************************************
!  Counts the number of unique adsorbed or gas phase species in a list
!
!  oct-14/Jonas: coded
!
    subroutine count_species_type(list,n_ad,n_sol,nlist)
      integer :: i, parenthes, nlist
      integer :: n_ad, n_sol
      character(len=10), dimension(:) :: list
!
! count adsorbed and surface species
      n_ad = 0
      n_sol = 0
      do i = 1,nlist
        parenthes = 0
        parenthes = scan(list(i),'()')
        if (parenthes > 0 .or. list(i) == 'Cf') then
          n_ad = n_ad + 1
        else
          if (list(i) /= 'Cb') then
            n_sol = n_sol + 1
          endif
        endif
      enddo
    endsubroutine count_species_type
! ******************************************************************************
!  Count the number of unique species in the parsed mechanics.in file
!  Given the parsed mechanics.in file and a list of species (adsorbed or gas phase)
!  creates lists having sorted the given list in a reactands and products.
!  If a species is reactand and product, it will appear on both
!
!  oct-14/Jonas: coded
!
    subroutine count_species(part,species,reactants,products)
      character(len=10), dimension(:,:) :: part
      character(len=10) :: element
      real :: numeric
      integer :: i, j, jmax, stat, place, number
      integer :: place_reac, place_prod
      logical :: lhs
      character(len=10), dimension(40) :: species, reactants, products
      character(len=10), dimension(40) :: temp, temp_reac, temp_prod
!
      temp = 'nothing'
      jmax = size(part,1)
      place = 1
      place_reac = 1
      place_prod = 1
!
      do i = 1,n_surface_reactions
        lhs = .true.
        do j = 1,jmax
          element = part(j,i)
!
! switch when the reaction arrow is read
          if (element == '->' .or. element == '<>') then
            lhs = .false.
          endif
!
! if element can be read as real, disregard
          read (element,*,iostat=stat) numeric
          if (stat /= 0 .and. &
              element /= '->' .and. &
              element /= '<>' .and. &
              element(:2) /= 'RR') then
            read (element(:1),*,iostat=stat) number
            if (stat == 0)  element = element(2:)
!
!  appending the components to the list
!  of global, reactand and product uniques
!
            if (.not. any(temp  ==  element)) then
              temp(place) = element
              place = place+1
            endif
!
            if ((.not. any(temp_reac  ==  element)) .and. lhs) then
              temp_reac(place_reac) = element
              place_reac = place_reac+1
            endif
!
            if ((.not. any(temp_prod == element)) .and. (lhs .eqv. .false.)) then
              temp_prod(place_prod) = element
              place_prod = place_prod+1
            endif
!
          endif
        enddo
      enddo
!
! creating the lists
      species(:place-1) = temp(:place-1)
      reactants(:place_reac-1) = temp_reac(:place_reac-1)
      products(:place_prod-1) = temp_prod(:place_prod-1)
      ns = place-1
      nr = place_reac-1
      np = place_prod-1
    endsubroutine count_species
! ******************************************************************************
!  If the reaction contains '<>', the string is passed to this subroutine
!  After the normal read in. The species are flipped around the '<>' element
!  and parsed again, this way, the backwards reaction is created
!
!  oct-14/Jonas: coded
!
    subroutine flip_and_parse(string,ireaction,target_list,direction)
      character(len=150) :: string, flipped_string
      character(len=50) :: lhs, sign, rhs, ende
      character(len=3), dimension(:) :: direction
      integer :: ireaction
      integer :: i, numerical, marker
      real :: real_number
      character(len=10), dimension(:,:) :: target_list
!
      marker = index(string,'<>')
      numerical = 1
      i = marker
!
      do while (numerical  /= 0 )
        i = i + 1
!
!  NILS: real_number is not necessarily a number when done like this, this
!  NILS: may cause problems for some compilers (e.g.
!  NILS: hosts/nordita/norlx51-daily-test.conf)
!
        read (string(i:i+7),*,iostat=numerical) real_number
        if (numerical == 0) then
          if (real_number < 10.) then
            numerical = 1
          endif
          if (i > len(string)-10) then
            print*,'no numericals found after sign!'
            numerical = 0
          endif
        endif
      enddo
!
      lhs = trim(string(:marker-1))
      sign = trim(string(marker:marker+1))
      rhs = trim(string(marker+2:i))
      ende = trim(string(i:))
!
      flipped_string = trim(rhs)//'  '//trim(sign)//'  '//trim(lhs)// trim(ende)
      flags(ireaction) = 'rev'
      call parse(flipped_string,ireaction,target_list,'rev',direction)
    endsubroutine flip_and_parse
! ******************************************************************************
!  Is given line by line of the mechanics.in file and divides the line by all
!  elements that are divided by any number of whitespaces in between. The '+'
!  common in the mechanics file is overlooked. If the last element is read in, all
!  Column elements afterwards are filled up with 0.0
!
!  oct-14/Jonas: coded
!
    subroutine parse(string,ireaction,target_list,flag,direction)
      character(len=150) :: string
      character :: tab = char(9)
      character :: spc = char(32)
      character(len=3) :: flag
      character(len=3), dimension(:) :: direction
      integer :: i, j, k, ireaction
      character(len=10) :: formatting
      character(len=10), dimension(:,:) :: target_list
!
      i = 1
      k = 1
      j = 1
!
      do while (i <= len(string))
        if (string(i:i) == tab) then
          string(i:i) = spc
        else
          i = i+1
        endif
      enddo
      string = trim(string)
      i = 1
      do while (i < len(string))
        if (string(i:i) /= '+' .and. string(i:i) /= spc) then
          j = 1
          do while (string(i+j:i+j) /= '+' .and. string(i+j:i+j) /= spc .and. string(i+j:i+j) /= tab)
            j = j+1
          enddo
          formatting = adjustl(string(i:i+j))
          target_list(k,ireaction) = formatting
          k = k+1
          i = i+j
        else
          i = i+1
        endif
      enddo
      if (i == len(string)) then
        target_list(k:,ireaction) = '0.0'
      endif
!   print*, target_list(:,ireaction)
!
      direction(ireaction) = flag
    endsubroutine parse
! ******************************************************************************
    subroutine read_mechanics_file(inputfile,target_list,n_max_elements, &
        reaction_direction,talk)
!
!  Open the mechanics.in file and if the line begins with something else than !
!  It is given to the parsing routine above. If a '<>' is encoutered, the line
!  is given to the flipping routine afterwards for creation of the backwards
!  reaction.
!
!  oct-14/Jonas: coded
!
      integer :: stat, ireaction, i, n_max_elements
      character(len=150) :: string
      character(len=*) :: inputfile, talk
      character(len=3), dimension(:) :: reaction_direction
      character(len=10), dimension(:,:) :: target_list
!
      writeformat = '(  A10,A3)'
      write (writeformat(2:3),'(I2)') n_max_elements
      open (20, file=inputfile,iostat=stat)
!
      if (stat == 0) then
        if (talk == 'verbose') then
          write (*,*) 'Opened mechanics file'
        endif
        ireaction = 1
        do while (stat == 0)
          read (20,fmt='(A150)', iostat=stat) string
          if (stat == 0) then
            if ((string(:1)) /= '!') then
              flags(ireaction) = 'fwd'
              call parse(string,ireaction,target_list,'fwd',reaction_direction)
              ireaction = ireaction + 1
              if (index(string,'<>') > 0) then
                call flip_and_parse(string,ireaction,target_list,reaction_direction)
                ireaction = ireaction + 1
              endif
            endif
          endif
        enddo
        if (talk == 'verbose') print*,'Done parsing mechanics file'
        close (20)
!
        call remove_save_T_k(target_list)
        call remove_save_powers(target_list)
!
        if (talk == 'verbose') then
          open (29, file='mech_outputfile.dat',iostat=stat)
          do i = 1,N_surface_reactions
            write (*,writeformat) target_list(:,i),reaction_direction(i)
            write (29,writeformat) target_list(:,i),reaction_direction(i)
          enddo
          close (29)
        endif
      else
        write (*,*) 'Could not open mechanics file'
      endif
    endsubroutine read_mechanics_file
! ******************************************************************************
!  Find the mole production of the forward reaction. This will later
!  be used for the calculation of the reverse reaction rate.
!
!  oct-14/Jonas: coded
!
    subroutine create_dngas()
      integer :: k
!
      do k = 1,N_surface_reactions
        dngas(k) = sum(nu_prime(:,k)) - sum(nu(:,k))
      enddo
    endsubroutine create_dngas
! ******************************************************************************
!  Calculate the change in entropy over a specific reaction. This is done by taking
!  the difference between the specific entropies of each product and species of one
!  reaction
!
!  oct-14/Jonas: coded
!
    subroutine calc_entropy_of_reaction()
      integer :: i, j, l, k, k1, k2
!
      entropy_k = 0
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      do k = k1,k2
        do l = 1,N_surface_reactions
          do i = 1,N_surface_species
            entropy_k(k,l) = entropy_k(k,l) &
                +nu_prime(i,l)*surface_species_entropy(k,i) &
                -nu(i,l)*surface_species_entropy(k,i)
          enddo
          if (N_adsorbed_species > 1) then
            do j = 1,N_adsorbed_species
              entropy_k(k,l) = entropy_k(k,l) &
                  +mu_prime(j,l)*adsorbed_species_entropy(k,j) &
                  -mu(j,l)*adsorbed_species_entropy(k,j)
            enddo
          endif
!  This is the cgs-SI switch
          entropy_k(k,l) = entropy_k(k,l) * pre_energy
        enddo
      enddo
    endsubroutine calc_entropy_of_reaction
! ******************************************************************************
!  Starting from the K_k of the forward reaction, calculate the K_k
!  of the backwards reaction by first obtaining the equilibrium constant
!  from the change in gibbs free energy over the forward reaction
!
!  oct-14/Jonas: coded
!
    subroutine get_reverse_K_k(l,fp,k)
      integer :: l, k
      real :: pre_pressure
      real :: k_p, k_c
      real :: denominator, expo
      real, dimension(:,:) :: fp
!
!  This should be kept in SI units, therefore a prefactor to the pressure is
!  introduced to keep interp_pp in SI units
!
      if (unit_system == 'cgs') then
        pre_pressure = 0.1
      else
        pre_pressure = 1.
      endif
!
      denominator = heating_k(k,l-1) - (entropy_k(k,l-1)*fp(k,iTp))
      expo = denominator/(Rgas*fp(k,iTp))
      k_p = exp(-expo)
!
! k_p == pressure independent equilibrium constant
!
      k_c = k_p / (((gas_constant)*fp(k,iTp)/ &
          (pre_pressure * interp_pp(k)))**(dngas(l-1)))
!
! k_c == pressure dependent equilibrium constant
!
      K_k(k,l) = (K_k(k,l-1) / k_c)
    endsubroutine get_reverse_K_k
! ******************************************************************************
!  Calculate the conversion of the particle
!
!  oct-14/Jonas: coded
!
    subroutine calc_conversion(fp)
      real, dimension(:,:) :: fp
      integer :: k
!
      do k = k1_imn(imn),k2_imn(imn)
        conversion(k) = 1 - fp(k,imp) / fp(k,impinit)
      enddo
    endsubroutine calc_conversion
! ******************************************************************************
!  Calculate the heating of each reaction by taking the difference of
!  specific enthalpies of products and reactands for each reaction
!
!  oct-14/Jonas: coded
!
    subroutine calc_enthalpy_of_reaction()
      integer :: i, k, k1, k2, j,l
!
      heating_k = 0
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      do k = k1,k2
        do l = 1,N_surface_reactions
          do i = 1,N_surface_species
            heating_k(k,l) = heating_k(k,l) &
                +nu_prime(i,l)*surface_species_enthalpy(k,i) &
                -nu(i,l)*surface_species_enthalpy(k,i)
          enddo
          if (N_adsorbed_species > 0) then
            do j = 1,N_adsorbed_species
              heating_k(k,l) = heating_k(k,l) &
                  +mu_prime(j,l)*adsorbed_species_enthalpy(k,j) &
                  -mu(j,l)*adsorbed_species_enthalpy(k,j)
            enddo
          endif
!  This is the cgs-SI switch
          heating_k(k,l) = heating_k(k,l) * pre_energy
        enddo
      enddo
    endsubroutine calc_enthalpy_of_reaction
! ******************************************************************************
!  01-oct-2014/jonas:coded
!
!  Calculates the area specific molar reaction rate from K_k
!
    subroutine calc_RR_hat(f,fp,ineargrid,p)
      type (pencil_case) :: p
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mx,my,mz,mfarray) :: f
      real :: pre_Cg, pre_Cs, pre_RR_hat
      integer :: i, j, k, k1, k2,l
      integer, dimension(mpar_loc,3) :: ineargrid
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
!  The heterogeneous kinetics in the mechanism file is always given in SI units.
!  Here we intorduce correction factors if cgs units are used.
!
      if (unit_system == 'cgs') then
        pre_Cg = 1e6
        pre_Cs = 1e4
        pre_RR_hat = 1e-4
      else
        pre_Cg = 1.
        pre_Cs = 1.
        pre_RR_hat = 1.
      endif
!
      do k = k1,k2
        do j = 1,N_surface_reactions
          RR_hat(k,j) = K_k(k,j)*reaction_enhancement(j)
          do i = 1,N_surface_reactants
            if (nu(i,j) > 0) RR_hat(k,j) = RR_hat(k,j)* &
                (pre_Cg * Cg(k)*fp(k,isurf-1+i))**nu(i,j)
          enddo
          if (N_adsorbed_species > 1) then
            do i = 1,N_adsorbed_species
              if (mu(i,j) > 0) RR_hat(k,j) = RR_hat(k,j)*(pre_Cs*Cs(k,i))**mu(i,j)
            enddo
          endif
          RR_hat(k,j) = RR_hat(k,j)*(fp(k,iTp)**T_k(j))
        enddo
      enddo
!
!  Make sure RR_hat is given in the right unit system (above it is always
!  in SI units).
!
      RR_hat = pre_RR_hat*RR_hat
!
! Find the maximum possible carbon consumption rate (i.e. the rate
! experienced in the case where the effectiveness factor is unity).
!
      Rck_max = 0.
      do k = k1,k2
        do l = 1,N_surface_reactions
          do i = 1,N_surface_species
            Rck_max(k,l) = Rck_max(k,l)+mol_mass_carbon*RR_hat(k,l) &
                *(nu_prime(i,l)-nu(i,l))*ac(i)
          enddo
        enddo
      enddo
!
! Find molar reaction rate of adsorbed surface species
      if (N_adsorbed_species > 1) then
        do k = k1,k2
          do l = 1,N_surface_reactions
            do j = 1,N_adsorbed_species-1
              Rck_max(k,l) = Rck_max(k,l)+mol_mass_carbon*RR_hat(k,l) &
                  *(mu_prime(j,l)-mu(j,l))*aac(j)
            enddo
          enddo
        enddo
      endif
!
!  Adapt the reaction rate according to the internal gradients,
!  after thiele. (8th US combustion Meeting, Paper #070CO-0312)
!  equation 56 ff.
!
      if (lthiele) then
        call calc_effectiveness_factor(fp,ineargrid,p)
        do j = 1,N_surface_reactions
          RR_hat(:,j) = RR_hat(:,j) * effectiveness_factor_reaction(:,j)
        enddo
      endif
!
    endsubroutine calc_RR_hat
! ******************************************************************************
! Then find the carbon consumption rate and molar flux of each gaseous
! species at the particle surface
!
!  oct-14/jonas (coded)
!
    subroutine calc_ndot_mdot_R_j_hat(fp)
      real, dimension(:,:) :: fp
      integer :: n, i, j, k, l, k1, k2
!
      ndot = 0.
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
      Rck = 0.
!
      do k = k1,k2
        do l = 1,N_surface_reactions
          do i = 1,N_surface_species
            Rck(k,l) = Rck(k,l)+mol_mass_carbon*RR_hat(k,l) &
                *(nu_prime(i,l)-nu(i,l))*ac(i)
            ndot(k,i) = ndot(k,i)+RR_hat(k,l)*(nu_prime(i,l)-nu(i,l))*St(k)/ &
                (fp(k,iap)**2*4.*pi)
          enddo
        enddo
      enddo
!
! Find molar reaction rate of adsorbed surface species
      if (N_adsorbed_species > 1) then
        R_j_hat = 0.
        do k = k1,k2
          do l = 1,N_surface_reactions
            do j = 1,N_adsorbed_species-1
              R_j_hat(k,j) = R_j_hat(k,j)+(mu_prime(j,l)-mu(j,l))*RR_hat(k,l)
              Rck(k,l) = Rck(k,l)+mol_mass_carbon*RR_hat(k,l) &
                  *(mu_prime(j,l)-mu(j,l))*aac(j)
            enddo
          enddo
        enddo
      endif
!
! Find mdot_ck
      do k = k1,k2
        do l = 1,N_surface_reactions
          mdot_ck(k,l) = -St(k)*Rck(k,l)
        enddo
        mass_loss(k) = sum(mdot_ck(k,:))
      enddo
!
! Find the sum of all molar fluxes at the surface
      do k = k1,k2
        ndot_total(k) = sum(ndot(k,:))
      enddo
    endsubroutine calc_ndot_mdot_R_j_hat
! ******************************************************************************
!  Get area specific molar reaction rate
!
!  oct-14/Jonas: coded
!
    subroutine get_RR_hat(var,start,end)
      real, dimension(:,:) :: var
      integer :: start, end
!
      var(start:end,:) = RR_hat(start:end,:)
    endsubroutine get_RR_hat
! ******************************************************************************
!  Get the parsed mechanics.in file
!
!  oct-14/Jonas: coded
!
    subroutine get_part(var)
      character(len=10), dimension(:,:) :: var
!
      var = part
    endsubroutine get_part
! ******************************************************************************
!  Get a list of all reactands
!
!  oct-14/Jonas: coded
!
    subroutine get_reactants(var)
      character(len=10), dimension(40) :: var
!
      intent(out) :: var
!
      var = reactants
    endsubroutine get_reactants
! ******************************************************************************
!  Calculates the efficiency factor for each reaction and particle by
!  comparing the maximum possible reaction rate to the one taking in account
!  the internal gradients of species concentration and the porosity of
!  the particle
!
!  01-Oct-2014/Jonas: coded
!  taken from solid_reac L. 149 and equations.pdf eq 35ff
!
    subroutine calc_effectiveness_factor(fp,ineargrid,p)
      real, dimension(:,:) :: fp
      type (pencil_case) :: p
!
      real, dimension(:,:), allocatable :: R_i_hat, D_eff
      real :: Knudsen, tmp3, tmp2, tmp1
      real, dimension(:), allocatable :: pore_radius
      real, dimension(:), allocatable ::  phi, sum_eta_i_R_i_hat_max
      real, dimension(:), allocatable :: sum_R_i_hat_max
      real, dimension(:), allocatable :: volume, porosity
      integer :: i, dep, k, k1, k2,l
      integer, dimension(mpar_loc,3) :: ineargrid
      real :: r_f,test
!
      test = 2.0
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      allocate(R_i_hat(k1:k2,N_surface_reactants))
      allocate(D_eff(k1:k2,N_surface_reactants))
      allocate(volume(k1:k2))
      allocate(pore_radius(k1:k2))
      allocate(phi(k1:k2))
      allocate(sum_R_i_hat_max(k1:k2))
      allocate(sum_eta_i_R_i_hat_max(k1:k2))
      allocate(porosity(k1:k2))
!
! Find reaction rate for each of the reactants
      R_i_hat = 0
      do k = k1,k2
        do l = 1,N_surface_reactions
          do i = 1,N_surface_reactants
            R_i_hat(k,i) = R_i_hat(k,i)+(nu_prime(i,l)-nu(i,l))*RR_hat(k,l)
          enddo
        enddo
      enddo
!
! Find particle volume
      do k = k1,k2
        volume(k) = 4.*pi*(fp(k,iap)**3)/3.
        porosity(k) = 1- (fp(k,imp)/volume(k))/true_density_carbon
      enddo
!
! Find pore radius
      r_f = 2.
      pore_radius(:) = 2*r_f*porosity(:)*volume(:)/St(k1:k2)
!
!  Find effective diffusion coefficient (based on bulk and Knudsen diffusion)
!
      do k = k1,k2
        do i = 1,N_surface_reactants
          tmp3 = 8*Rgas*fp(k,iTp)/(pi*(species_constants(jmap(i),imass)))
          Knudsen = 2*pore_radius(k)*porosity(k)*sqrt(tmp3)/(3*tortuosity)
          tmp1 = 1./p%Diff_penc_add(ineargrid(k,1)-nghost,jmap(i))
          tmp2 = 1./Knudsen
          D_eff(k,i) = 1./(tmp1+tmp2)
        enddo
      enddo
!      print*,'D_eff', D_eff
!
!  Find thiele modulus and effectiveness factor
!  JONAS: was with volume(k) before,
!
      do k = k1,k2
        do i = 1,N_surface_reactants
          if (R_i_hat(k,i) < 0.0) then
            thiele(k,i) = fp(k,iap)*sqrt(-R_i_hat(k,i)*St(k)/ &
                (fp(k,imp)*Cg(k)*fp(k,isurf-1+i)*D_eff(k,i)))
          else
            thiele(k,i) = 1.e-2
          endif
!          print*, 'thiele: ',thiele
!
! calculate the effectivenes factor (all particles, all reactants)
          effectiveness_factor_species(k,i) = 3/thiele(k,i)* &
              (1./tanh(thiele(k,i))-1./thiele(k,i))
!          print*,'effectiveness', effectiveness_factor_species
        enddo
      enddo
!
!  The mean effectiveness factor can be found from the above by using R_i_max
!  (all particles, all reactants)
!
      sum_R_i_hat_max = 0.0
!
      do k = k1,k2
        do i = 1,N_surface_reactants
          if (R_i_hat(k,i) < 0) then
            sum_R_i_hat_max(k) = sum_R_i_hat_max(k)+R_i_hat(k,i)
          endif
        enddo
      enddo
!
      do k = k1,k2
        sum_eta_i_R_i_hat_max(k) = 0
        do i = 1,N_surface_reactants
          if (R_i_hat(k,i) < 0) then
            sum_eta_i_R_i_hat_max(k) = sum_eta_i_R_i_hat_max(k) &
                +R_i_hat(k,i)*effectiveness_factor_species(k,i)
          endif
        enddo
      enddo
!      print*,'sum_eta_i_R_i_hat_max',sum_eta_i_R_i_hat_max
!      print*,'sum_R_i_hat_max',sum_R_i_hat_max
!
      do k = k1,k2
        if (sum_R_i_hat_max(k)/=0.0) then
          effectiveness_factor(k) = sum_eta_i_R_i_hat_max(k)/sum_R_i_hat_max(k)
        endif
      enddo
!
!  The efficiency factor for each reaction can now be found from the thiele modulus
!
!  JONAS: the iteration parameters have switched from solid_reac in order to
!  be consistent! be careful when checking!
!
      do k = k1,k2
        do l = 1,N_surface_reactions
          dep = dependent_reactant(l)
          if (dep > 0) then
            if (R_i_hat(k,dep) > 0.) then
              effectiveness_factor_reaction(k,l) = 1.
            else
              phi(k) = thiele(k,dependent_reactant(l))
              effectiveness_factor_reaction(k,l) = 3* &
                  (1./tanh(phi(k))-1./phi(k))/phi(k)
            endif
          else
            effectiveness_factor_reaction(k,l) = 1
          endif
        enddo
      enddo
!
      deallocate(R_i_hat)
      deallocate(D_eff)
      deallocate(volume)
      deallocate(pore_radius)
      deallocate(sum_R_i_hat_max)
      deallocate(sum_eta_i_R_i_hat_max)
      deallocate(phi)
      deallocate(porosity)
    endsubroutine calc_effectiveness_factor
! ******************************************************************************
!  Calculate specific enthalpy of each gas phase species on the
!  particle surface
!
!  oct-14/Jonas: coded
!
    subroutine calc_surf_enthalpy(fp)
      real, dimension(mpar_loc,mparray) :: fp
      integer :: k, k1, k2
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
! Unit: J/(mol)
      do k = k1,k2
        if (inuH2O > 0) surface_species_enthalpy(k,inuH2O) = &
            -242.18e3-(5.47*fp(k,iTp))
        if (inuO2 > 0)  surface_species_enthalpy(k,inuO2) = 0.
        if (inuCO2 > 0) surface_species_enthalpy(k,inuCO2) = &
            -392.52e3-(2.109*fp(k,iTp))
        if (inuH2 > 0)  surface_species_enthalpy(k,inuH2 ) = 0.
        if (inuCO > 0)  surface_species_enthalpy(k,inuCO ) = &
            -105.95e3-6.143*fp(k,iTp)
        if (inuCH > 0)  surface_species_enthalpy(k,inuCH ) = 594.13e3
        if (inuHCO > 0) surface_species_enthalpy(k,inuHCO) = &
            45.31e3-(5.94*fp(k,iTp))
        if (inuCH2 > 0) surface_species_enthalpy(k,inuCH2) = &
            387.93e3-(5.8*fp(k,iTp))
        if (inuCH4 > 0) surface_species_enthalpy(k,inuCH4) = -75e3
        if (inuCH3 > 0) surface_species_enthalpy(k,inuCH3) = &
            144.65e3-(6.79*fp(k,iTp))
      enddo
    endsubroutine calc_surf_enthalpy
! ******************************************************************************
!  Calculate the specific entropy of the gas phase species at the
!  particle surface
!
!  oct-14/Jonas: coded
!
    subroutine calc_surf_entropy(fp)
      real, dimension(:,:) :: fp
      integer :: k, k1, k2
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
! JONAS: units are in j/(mol*K)
      do k = k1,k2
        if (inuH2O > 0) surface_species_entropy(k,inuH2O) = &
            189.00+(0.0425*fp(k,iTp))
        if (inuO2 > 0)  surface_species_entropy(k,inuO2) = &
            222.55+(0.0219*fp(k,iTp))
        if (inuCO2 > 0) surface_species_entropy(k,inuCO2) = &
            212.19+(0.0556*fp(k,iTp))
        if (inuH2 > 0)  surface_species_entropy(k,inuH2 ) = &
            133.80+(0.0319*fp(k,iTp))
        if (inuCO > 0)  surface_species_entropy(k,inuCO ) = &
            199.35+(0.0342*fp(k,iTp))
!
! taken from chemistry  webbook (1bar)
        if (inuCH > 0)  surface_species_entropy(k,inuCH ) = 183.00
        if (inuHCO > 0) surface_species_entropy(k,inuHCO) = &
            223.114+(0.0491*fp(k,iTp))
        if (inuCH2 > 0) surface_species_entropy(k,inuCH2) = &
            193.297+(0.0467*fp(k,iTp))
!
! taken from chemistry webbook (1bar)
        if (inuCH4 > 0) surface_species_entropy(k,inuCH4) = 189.00
        if (inuCH3 > 0) surface_species_entropy(k,inuCH3) = &
            190.18+(0.0601*fp(k,iTp))
!
      enddo
    endsubroutine calc_surf_entropy
! ******************************************************************************
!  Calculate the specific entropy of adsorbed species for each particle
!
!  oct-14/Jonas: coded
!
    subroutine calc_ads_entropy(fp)
      real, dimension(mpar_loc,mparray) :: fp
      integer :: k, k1, k2
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
! Unit: J/(mol*K)
      do k = k1,k2
!
        if (imuadsO > 0)    then
          adsorbed_species_entropy(k,imuadsO) = &
              (164.19+(0.0218*fp(k,iTp)))*0.72 - (3.3*gas_constant)
        endif
!
! this is guessed
        if (imuadsO2 > 0)   then
          adsorbed_species_entropy(k,imuadsO2) =  &
              2*adsorbed_species_entropy(k,imuadsO)
        endif
        if (imuadsOH > 0)   then
          adsorbed_species_entropy(k,imuadsOH) = &
              ((0.0319*fp(k,iTp)) + 186.88) * 0.7 - (3.3*gas_constant)
        endif
        if (imuadsH > 0)    then
          adsorbed_species_entropy(k,imuadsH) = &
              (117.49+(0.0217*fp(k,iTp)))*0.54 - (3.3*gas_constant)
        endif
        if (imuadsCO > 0)   then
          adsorbed_species_entropy(k,imuadsCO) = &
              (199.35+(0.0342*fp(k,iTp))) * &
              0.6*(1+(1.44e-4*fp(k,iTp))) - (3.3*gas_constant)
        endif
!
! taken from nist
        if (imufree > 0)    adsorbed_species_entropy(k,imufree) = 0
!
      enddo
    endsubroutine calc_ads_entropy
! ******************************************************************************
!  Calculate the specific enthalpy of adsorbed species for each particle
!
!  oct-14/Jonas: coded
!
    subroutine calc_ads_enthalpy(fp)
      real, dimension(:,:) :: fp
      integer :: k, k1, k2
!
!  JONAS: values are from nist and solid_phase.f90 of the stanford
!  code Units: J/mol
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      if (imuadsO > 0)    adsorbed_species_enthalpy(k1:k2,imuadsO) = &
          -148.14e3 + (0.0024e3*(fp(k1:k2,iTp)-273.15))
      if (imuadsO2 > 0)   adsorbed_species_enthalpy(k1:k2,imuadsO2) = &
          2 *  (-148.14e3 + (0.0024e3*(fp(k1:k2,iTp)-273.15)))
      if (imuadsOH > 0) adsorbed_species_enthalpy(k1:k2,imuadsOH) = -148e3
      if (imuadsH > 0) adsorbed_species_enthalpy(k1:k2,imuadsH) = -19.5e3
      if (imuadsCO > 0)   adsorbed_species_enthalpy(k1:k2,imuadsCO) = &
          -199.94e3 - (0.0167e3*(fp(k1:k2,iTp)-273.15))
      if (imufree > 0)    adsorbed_species_enthalpy(k1:k2,imufree) = 0.
    endsubroutine calc_ads_enthalpy
! ******************************************************************************
!  Wrapper routine to calculate all terms needed to be used in the variables to
!  be evolved on each particle on the current pencil
!  These are also used in the particle mass, temperature, surfspec and adsorbed
!  modules
!
!  oct-14/Jonas: coded
!
    subroutine calc_pchemistry_pencils(f,fp,p,ineargrid)
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mx,my,mz,mfarray) :: f
      integer, dimension(mpar_loc,3) :: ineargrid
      type (pencil_case) :: p
      integer :: k1,k2      
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      if (allocated(part)) deallocate(part)
      if (allocated(part_power)) deallocate(part_power)
!
! Routine to calcute quantities used for reactive particles
!
      if (k1 /= 0 .and. k2 /= 0) then
        call allocate_variable_pencils()
!
        call calc_rho_p(fp)
        call calc_conversion(fp)
        call calc_St(fp)
!
        call calc_ads_entropy(fp)
        call calc_ads_enthalpy(fp)
        call calc_surf_entropy(fp)
        call calc_surf_enthalpy(fp)
!
        call calc_enthalpy_of_reaction()
        call calc_entropy_of_reaction()
!
        call calc_K_k(f,fp,p,ineargrid)
        call calc_Cg(fp,p,ineargrid)
        if (N_adsorbed_species > 1) call calc_Cs(fp)
        call calc_A_p(fp)
        call calc_RR_hat(f,fp,ineargrid,p)
!
        if (lreactive_heating) then
          call calc_q_reac()
!JONAS:  Commented out the calculation of the Nusselt number for now
!JONAS:  since this has to be fixed before it can be used.
!        call calc_Nusselt()
          Nu_p = 2.
!JONAS:  We also need to make a subroutine calculcating the Sherwood
!JONAS:  number.
        else
          Nu_p = 2.
          q_reac = 0.0
        endif
!
        call calc_ndot_mdot_R_j_hat(fp)
        call calc_R_c_hat()
    endif
!
    endsubroutine calc_pchemistry_pencils
! ******************************************************************************
!  allocate variables used in the pencil variable calculation
!
!  Allocate terms needed for the attributes of the particles to be evolved
!  Allocation is for all particles on the current pencil
!
!  oct-14/Jonas: coded
!
    subroutine allocate_variable_pencils()
      integer :: k1, k2, stat, ierr
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      allocate(rho_p(k1:k2),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for rho_p')
      allocate(q_reac(k1:k2),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for q_reac')
      allocate(St(k1:k2),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for St')
      allocate(mass_loss(k1:k2),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for mass_loss')
!
      allocate(conversion(k1:k2),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for conversion')
!
      allocate(surface_species_enthalpy(k1:k2,N_surface_species), STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for surface_species_enthalpy')
      allocate(surface_species_entropy(k1:k2,N_surface_species), STAT=stat)
      allocate(adsorbed_species_entropy(k1:k2,N_adsorbed_species), STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for adsorbed_species_entropy')
      allocate(adsorbed_species_enthalpy(k1:k2,N_adsorbed_species), STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for adsorbed_species_enthalpy')
!
      allocate(thiele(k1:k2,N_surface_reactants),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for thiele')
!
      allocate(effectiveness_factor(k1:k2),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for thiele')
      allocate(effectiveness_factor_species(k1:k2,N_surface_reactants),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for effectiveness_factor_species')
      allocate(effectiveness_factor_reaction(k1:k2,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for effectiveness_factor_reaction')
!
      allocate(ndot(k1:k2,N_surface_species),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for ndot')
      allocate(ndot_total(k1:k2),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for ndot_total')
!
      allocate(Rck(k1:k2,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for Rck')
      allocate(Rck_max(k1:k2,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for Rck_max')
!
      allocate(R_j_hat(k1:k2,N_adsorbed_species-1),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for R_j_hat')
      allocate(Cs(k1:k2,N_adsorbed_species),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for Cs')
      allocate(RR_hat(k1:k2,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for RR_hat')
      allocate(K_k(k1:k2,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for K_k')
      allocate(Cg(k1:k2), STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for Cg')
      allocate(A_p(k1:k2), STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for A_p')
      allocate(Nu_p(k1:k2), STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for Nu_p')
      allocate(mdot_ck(k1:k2,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_dep_pchem', &
          'Could not allocate memory for mdot_ck')
    endsubroutine allocate_variable_pencils
! ******************************************************************************
!  Deallocate variables after use
!
!  oct-14/Jonas: coded
!
    subroutine cleanup_chemistry_pencils()
!
!  Do only if particles are present on the pencil
!
      if (npar_imn(imn) /= 0) then
        deallocate(surface_species_enthalpy)
        deallocate(surface_species_entropy)
        deallocate(adsorbed_species_entropy)
        deallocate(adsorbed_species_enthalpy)
        deallocate(thiele)
        deallocate(effectiveness_factor)
        deallocate(effectiveness_factor_species)
        deallocate(effectiveness_factor_reaction)
        deallocate(ndot)
        deallocate(conversion)
        deallocate(ndot_total)
        deallocate(Rck)
        deallocate(RR_hat)
        deallocate(Rck_max)
        deallocate(R_j_hat)
        deallocate(Cs)
        deallocate(rho_p)
        deallocate(St)
        deallocate(K_k)
        deallocate(Cg)
        deallocate(A_p)
        deallocate(q_reac)
        deallocate(Nu_p)
        deallocate(mdot_ck)
        deallocate(mass_loss)
      endif
    endsubroutine cleanup_chemistry_pencils
! ******************************************************************************
!  Read input parameters for the particles_chemistry module in start.x
!
!  oct-14/Jonas: coded
!
    subroutine read_particles_chem_init_pars(iostat)
      use File_io, only: parallel_unit
      integer, intent(out) :: iostat
!
      read (parallel_unit, NML=particles_chem_init_pars, IOSTAT=iostat)
    endsubroutine read_particles_chem_init_pars
! ******************************************************************************
!  Write input parameters for the particles_chemistry module in start.x
!
!  oct-14/Jonas: coded
!
    subroutine write_particles_chem_init_pars(unit)
      integer, intent(in) :: unit
!
      write (unit, NML=particles_chem_init_pars)
    endsubroutine write_particles_chem_init_pars
! ******************************************************************************
!  Read run parameters for particles_chemistry module, run.x
!
!  oct-14/Jonas: coded
!
    subroutine read_particles_chem_run_pars(iostat)
      use File_io, only: parallel_unit
      integer, intent(out) :: iostat
!
      read (parallel_unit, NML=particles_chem_run_pars, IOSTAT=iostat)
    endsubroutine read_particles_chem_run_pars
! ******************************************************************************
!   Write run parameters for particles_chemistry module, run.x
!
!  oct-14/Jonas: coded
!
    subroutine write_particles_chem_run_pars(unit)
      integer, intent(in) :: unit
!
      write (unit, NML=particles_chem_run_pars)
    endsubroutine write_particles_chem_run_pars
! ******************************************************************************
!  Calculate the current density of the particle
!
!  07-oct-14/Jonas: coded
!
    subroutine calc_rho_p(fp)
      real, dimension(:,:) :: fp
      integer :: k, k1, k2
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      do k = k1,k2
        rho_p(k) = fp(k,imp) / (4./3. *pi * fp(k,iap)**3)
      enddo
    endsubroutine calc_rho_p
! ******************************************************************************
!  Looks for terms of ^x.yz after the species for arbitrary reaction order
!  of species concentrations. The power terms are removed from the list and
!  stored in an own part_power variable used in create_stoc
!
!  10-oct-14/Jonas: coded
!
    subroutine remove_save_powers(target_list)
      character(len=10), dimension(:,:) :: target_list
      character(len=10) :: element, writeformat
      real :: power
      integer :: i, j, pow_pla
      integer :: isreal
!
      writeformat = '(  F6.2)'
      write (writeformat(2:3),'(I2)') N_max_elements
!
      i = 1
      do while (i <= N_surface_reactions)
        j = 1
        do while (j <= N_max_elements)
          element = target_list(j,i)
          if (index(element,'^') > 0) then
            pow_pla = index(element,'^')
            read (element(pow_pla+1:pow_pla+4),*,iostat=isreal) power
            target_list(j,i) = element(1:pow_pla-1)
            if (isreal == 0) then
              part_power(j,i) = power
            else
              print*, 'wrong format of power, needs to be x.yz'
            endif
          else
            part_power(j,i) = 1.
          endif
          j = j+1
        enddo
        i = i+1
      enddo
    endsubroutine remove_save_powers
! ******************************************************************************
!  Some printouts to look if the mechanics.in file has been read in correctly
!
!  oct-14/Jonas: coded
!
    subroutine print_debug_info()
      integer :: k
!
      writeformat = '(A12," ",I4,  F7.2)'
      write (writeformat(13:14),'(I2)') N_surface_species
!
      do k = 1,N_surface_reactions
        write (*,writeformat) 'nu=',k,nu(:,k)
      enddo
!
      do k = 1,N_surface_reactions
        write (*,writeformat) 'nu_prime=',k,nu_prime(:,k)
      enddo
!
      write (writeformat(13:14),'(I2)') N_adsorbed_species
!
      do k = 1,N_surface_reactions
        write (*,writeformat) 'mu=',k,mu(:,k)
      enddo
!
      do k = 1,N_surface_reactions
        write (*,writeformat) 'mu_prime=',k,mu_prime(:,k)
      enddo
!
      do k = 1,N_surface_reactions
        write (*,'(A12,I4,2E12.5)') 'ER_k, B_k=',k,B_k(k),ER_k(k)/ (gas_constant)
      enddo
      do k = 1,N_surface_reactions
        write (*,'(A12,I4,E12.5)') 'Dngas',k,dngas(k)
      enddo
      do k = 1,N_surface_reactions
        write (*,'(A12,I4,E12.5)') 'sigma',k,sigma_k(k)
      enddo
      do k = 1,N_surface_reactions
        write (*,'(A12,I4,I4)') 'dep',k,dependent_reactant(k)
      enddo
!
      write (*,'(A20," ",10I4)') 'jmap=', jmap
      write (*,'(A20," ",10F4.0)') 'ac=',ac
      write (*,'(A20," ",10F4.0)') 'site_occupancy=',aac
    endsubroutine print_debug_info
! ******************************************************************************
!  Calculation of K_k from the Arrhenius factors
!
!  27-10-2014/jonas:coded
!
    subroutine calc_K_k(f,fp,p,ineargrid)
      real, dimension(:,:,:,:) :: f
      real, dimension(:,:) :: fp
      type (pencil_case) :: p
      integer, dimension(:,:) :: ineargrid
      integer :: k, k1, k2, i, j,l
      integer :: N_iter = 20
      real :: int_k, gg_old, gg, ff, energy, dE, delta_E
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
! Calculation of the 'plain' K_k (reaction rate)
      do k = k1,k2
        do l = 1,N_surface_reactions
          if (sigma_k(l) == 0) then
            K_k(k,l) = B_k(l)*exp(-ER_k(l)/fp(k,iTp))
          else
            delta_E = sigma_k(l)*6
            dE = 2*delta_E/(N_iter-1)
            energy = ER_k(l)*gas_constant-delta_E
            if (energy < 0) then
              print*,'k,delta_E,sigma_k(k),ER_k(k)=', l,delta_E,sigma_k(l),ER_k(l)
              call fatal_error('get_reaction_rates_solid', &
                  'delta_E is too large!')
            endif
            ff = exp(-0.5*((energy-ER_K(l)*gas_constant)/sigma_k(l))**2)/ &
                (sigma_k(l)*sqrt(2*pi))
            gg_old = ff*B_k(l)*exp(-energy/(gas_constant*fp(k,iTp)))
            int_k = 0
            do j = 2,N_iter
              energy = (j-1)*2*delta_E/(N_iter-1)+ER_k(l)*gas_constant-delta_E
              ff = exp(-0.5*((energy-ER_K(l)*gas_constant)/sigma_k(l))**2)/ &
                  (sigma_k(l)*sqrt(2*pi))
              gg = ff*B_k(l)*exp(-energy/(gas_constant*fp(k,iTp)))
              int_k = int_k+dE*(gg_old+0.5*(gg-gg_old))
              gg_old = gg
            enddo
            K_k(k,l) = int_k
          endif
          if (reaction_direction(l) == 'rev') then
            call get_reverse_K_k(l,fp,k)
          endif
        enddo
      enddo
!
!  Application of the startup-quenching to facilitate convergence on
!  startup
!
      if (t < startup_time) then
        startup_quench = (tanh(6*(t-startup_time/2)/(startup_time))+1.0)/2.
        if (iter < 1) startup_quench = 0.
        K_k = K_k*startup_quench
      endif
!      print*,'k_k   ', k_k
!      print*,'B_k   ', B_k
!      print*,'ER_k  ', ER_k
!      print*,'fp(T) ',fp(:,iTp)
!
    endsubroutine calc_K_k
! ******************************************************************************
!  Gas concentration in the cell of the particle
!
!  oct-14/Jonas: coded
!
    subroutine calc_Cg(fp,p,ineargrid)
      real, dimension(:,:) :: fp
      integer :: k, k1, k2
      integer :: ix0
      real :: Tfilm
      type (pencil_case) :: p
      integer, dimension(:,:) :: ineargrid
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      do k = k1,k2
        Tfilm=fp(k,iTp)+(interp_TT(k)-fp(k,iTp))/3.
        Cg(k) = interp_pp(k)/(Rgas*Tfilm)
      enddo
    endsubroutine calc_Cg
! ******************************************************************************
!  Calculate particle surface area
!
!  oct-14/Jonas: coded
!
    subroutine calc_A_p(fp)
      real, dimension(:,:) :: fp
      integer :: k1, k2
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      A_p(k1:k2) = 4 * pi * fp(k1:k2,iap) * fp(k1:k2,iap)
    endsubroutine calc_A_p
! ******************************************************************************
!  calculate site concentration
!
!  oct-14/Jonas: coded
!
    subroutine calc_Cs(fp)
      real, dimension(:,:) :: fp
      integer :: k, k1, k2,i
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      do k = k1,k2
        do i = 1,N_adsorbed_species-1
          Cs(k,i) = fp(k,iads+i-1)
        enddo
        Cs(k,N_adsorbed_species) = 1 - sum(Cs(k,1:N_adsorbed_species-1))
      enddo
!
      Cs = Cs * total_carbon_sites
    endsubroutine calc_Cs
! ******************************************************************************
!  calculate the reactive heating of the particle
!
!  oct-14/Jonas: coded
!
    subroutine calc_q_reac()
!
      integer :: k, k1, k2, i, ierr
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      q_reac = 0.0
      do k = k1,k2
        do i = 1,N_surface_reactions
          q_reac(k) = q_reac(k)+RR_hat(k,i)*St(k)*(-heating_k(k,i))
        enddo
      enddo
!
    endsubroutine calc_q_reac
! ******************************************************************************
!
!  Calculate the Nusselt number of each particle. This is still incomplete
!
!  oct-14/Jonas: coded
!
    subroutine calc_Nusselt()
!
      use SharedVariables, only: put_shared_variable
!
      integer :: k, k1, k2
      real :: Pr_g, Sc_g, rep, ierr
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
!  JONAS: access to some chemistry variables as well to particles_dust for rep(k)
!  for Pr and Sc is needed, these are only placeholders!!!
!
      call fatal_error('calc_Nusselt', &
          'Prandtl, Schmidt and Reynolds number are wrong!')
      !NILS: We must do this correctly now, otherwise we are likely to
      !NILS: forget and we will get wrong results.
      Pr_g = 0.7
      Sc_g = 0.65
      rep = 0.0002
!
      do k = k1,k2
!  case from the paper "Dimensionless heat-mass transfer coefficients
!  for forced convection around a sphere: a general low
!  Reynolds number correlation
        if (rep < 3.5) then
          Nu_p(k) = 0.5+0.5*(1+2*rep*Pr_g)**(1./3.)
!  First case in "The effect of particle packing on the Nusselt
!  number in a cluster of particles
        elseif (rep >= 3.5 .and. rep <= 7.6e4 .and. Pr_g >= 0.7 .and. Pr_g <= 380.) then
          Nu_p(k) = 2.0+(0.4*rep**0.5+0.06*rep**0.67)* Pr_g**0.4
        else
          Nu_p = 2.0
        endif
      enddo
!
!    call put_shared_variable('Nu_p',Nu_p,ierr)
!    if (ierr/=0) call fatal_error('particles_chemistry:',&
!        'failed to put Nu_p')
!
    endsubroutine calc_Nusselt
! ******************************************************************************
!  This subroutine parses trough the list of elements of the mechanics.in
!  and looks for terms of the form T^x.yz where x.yz is a real number.
!  This is saved in T_k and is for temperature dependent reactions
!
!  oct-14/Jonas: coded
!
    subroutine remove_save_T_k(target_list)
      character(len=10), dimension(:,:) :: target_list
      character(len=10) :: el_T_k
      real :: T_k_single
      integer :: i, k, stat
!
      T_k = 0.0
      do i = 1,size(target_list,2)
        do k = 1,size(target_list,1)
          if (target_list(size(target_list,1),i) == 'rev') then
          else
            el_T_k = target_list(k,i)
            if (el_T_k(1:2) == 'T^') then
              read (el_T_k(3:),*,iostat=stat) T_k_single
              if (stat == 0) then
                T_k(i) = T_k_single
                target_list(k:size(target_list,1)-2,i) = &
                    target_list(k+1:size(target_list,1)-1,i)
                target_list(size(target_list,1)-1,i) = '0.0'
              else
                call fatal_error('create_arh_param', &
                    'T^ was given, but exponent could not be read!')
              endif
            endif
          endif
        enddo
      enddo
    endsubroutine remove_save_T_k
! ******************************************************************************
!  This function counts the maximum number of parsable elements
!  in the mechanics.in file, namely the species, reaction direction
!  activation energy, preexponential factor and optional, distribution of
!  activation energy and  temperature dependency
!
!  oct-14/Jonas: coded
!
    function count_max_elements(inputfile)
      integer :: count_max_elements
      character(len=*) :: inputfile
      integer :: stat, i, k, j, maxelement
      character(len=150) :: line
      character :: tab = char(9)
      character :: spc = char(32)
!
      open (30, file=inputfile,iostat=stat)
      if (stat == 0) then
        maxelement = 1
        do while (stat == 0)
          read (30,'(A150)', iostat=stat) line
          if (stat == 0) then
            if (line(:1) /= '!') then
              i = 1
              k = 0
              do while (i < len(line))
                if (line(i:i) /= '+' .and. line(i:i) /= spc .and. line(i:i) /= tab) then
                  j = 1
                  do while (line(i+j:i+j) /= '+' .and. &
                        line(i+j:i+j) /= spc .and. line(i+j:i+j) /= tab)
                    j = j+1
                  enddo
                  k = k+1
                  i = i+j
                else
                  i = i+1
                endif
                if (k > maxelement) then
                  maxelement = k
                endif
              enddo
            endif
          endif
        enddo
        close (30)
        count_max_elements = maxelement
      else
        maxelement = 0
        write (*,*) 'Problem with the single elements of the mechanics.in file'
      endif
    endfunction count_max_elements
! ******************************************************************************
!  this function counts the uncommented lines in the mechanics.in
!
!  oct-14/Jonas: coded
!
    function count_reactions(inputfile)
      integer :: count_reactions
      integer :: stat, reactions
      character(len=150) :: line
      character(len=*) :: inputfile
!
      open (20, file=inputfile,iostat=stat)
      if (stat == 0) then
        reactions = 0
        do while (stat == 0)
          read (20,fmt ='(A150)', iostat=stat) line
          if (stat == 0) then
            if (line(:1) /= '!') then
              reactions = reactions + 1
              if (index(line,'<>') > 0) then
                reactions = reactions + 1
              endif
            endif
          endif
        enddo
        close (20)
        count_reactions = reactions
      else
        count_reactions = 0
        write (*,*) 'Could not open mechanics file'
      endif
    endfunction count_reactions
! ******************************************************************************
!  Function to replace the imu/inuX variables
!
!  oct-14/Jonas: coded
!
    function find_species(species,unique_species,nlist)
      implicit none
!
      integer :: find_species
      integer :: i, nlist
      character(len=*) :: species
      character(len=*), dimension(:) :: unique_species
!
      find_species = 0
!
      do i = 1,nlist
        if (trim(species) == trim(unique_species(i))) then
          find_species = i
        endif
      enddo
    endfunction find_species
! ******************************************************************************
!  Get mass-chemistry dependent variables!
!  oct-14/Jonas: coded
!
    subroutine get_mass_chemistry(mass_loss_targ,St_targ,Rck_max_targ)
      real, dimension(:) :: mass_loss_targ, St_targ
      real, dimension(:,:) :: Rck_max_targ
!
      mass_loss_targ = mass_loss
      St_targ = St
      Rck_max_targ = Rck_max
    endsubroutine get_mass_chemistry
! ******************************************************************************
!  Get radius-chemistry dependent variables!
!  oct-14/Jonas: coded
!
    subroutine get_radius_chemistry(mass_loss_targ,effectiveness_targ)
      real, dimension(:) :: mass_loss_targ, effectiveness_targ
!
      mass_loss_targ = mass_loss
      effectiveness_targ = effectiveness_factor
    endsubroutine get_radius_chemistry
! ******************************************************************************
!  Get adsorbed-chemistry dependent variables!
!  oct-14/Jonas: coded
!
    subroutine get_adsorbed_chemistry(R_j_targ,R_c_targ)
      real, dimension(:) :: R_c_targ
      real, dimension(:,:) :: R_j_targ
!
      R_j_targ = R_j_hat
      R_c_targ = R_c_hat
!
    endsubroutine get_adsorbed_chemistry
! ******************************************************************************
!  Get surface-chemistry dependent variables!
!  oct-14/Jonas: coded
!
    subroutine get_surface_chemistry(Cg_targ,ndot_targ)
      real, dimension(:,:), optional :: ndot_targ
      real, dimension(:) :: Cg_targ
!
      if (present(ndot_targ)) ndot_targ = ndot
      Cg_targ = Cg
!
    endsubroutine get_surface_chemistry
! ******************************************************************************
!  Get temperature-chemistry dependent variables!
!  dec-11/Jonas: coded
!
    subroutine get_temperature_chemistry(q_reac_targ,Nu_p_targ)
      real, dimension(:) :: q_reac_targ, Nu_p_targ
!
      Nu_p_targ = Nu_p
      q_reac_targ = q_reac
!      q_reac_targ = q_reac
!
    endsubroutine get_temperature_chemistry
! ******************************************************************************
    subroutine particles_chemistry_clean_up
!
      if (allocated(dngas)) deallocate(dngas)
      if (allocated(omega_pg_dbl)) deallocate(omega_pg_dbl)
      if (allocated(B_k)) deallocate(B_k)
      if (allocated(Er_k)) deallocate(Er_k)
      if (allocated(sigma_k)) deallocate(sigma_k)
      if (allocated(reaction_order)) deallocate(reaction_order)
      if (allocated(effectiveness_factor_old)) deallocate(effectiveness_factor_old)
      if (allocated(reaction_direction)) deallocate(reaction_direction)
      if (allocated(flags)) deallocate(flags)
      if (allocated(T_k)) deallocate(T_k)
      if (allocated(R_c_hat)) deallocate(R_c_hat)
      if (allocated(heating_k)) deallocate(heating_k)
      if (allocated(entropy_k)) deallocate(entropy_k)
      if (allocated(part)) deallocate(part)
      if (allocated(part_power)) deallocate(part_power)
!
    endsubroutine particles_chemistry_clean_up
! ******************************************************************************
    subroutine register_unit_system()
!
!  Switch some variables that are affected by the switch from SI to cgs.
!  They are declared in cgs and need to be switched to SI. Do nothing if
!  cgs
!  Enthalpy and Entropie have to get a pre factor if cgs
!
      if (unit_system == 'cgs') then
        ! this must be used since the enthalpies and entropies are always
        ! given in SI units
        pre_energy = pre_energy * 1.0e7 !J to ergs
      else
        mol_mass_carbon = mol_mass_carbon/1000.
      endif
!
    endsubroutine register_unit_system
! ******************************************************************************
endmodule Particles_chemistry
