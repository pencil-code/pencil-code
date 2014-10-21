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
module Particles_chemistry
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
  use Particles_radius
!
  implicit none
!
  include 'particles_chemistry.h'
!
  public :: mu,mu_prime,ac,aac,nu,nu_prime, adsorbed_species_entropy
  public :: adsorbed_species_enthalpy
  public :: N_surface_reactions, N_adsorbed_species
  public :: N_surface_reactants, N_surface_species
  public :: inuH2,inuCO2,inuH2O,inuCO,inuCH4,inuO2
  public :: inuCH,inuHCO,inuCH2,inuCH3
  public :: imufree, imuadsO, imuadsO2, imuadsOH, imuadsH, imuadsCO
  public :: iads, iads_end
  public :: isurf,isurf_end
  public :: R_j_hat
  public :: j_of_inu
  public :: solid_species
  public :: mass_trans_coeff_reactants
  public :: part
  public :: diff_coeff_reactants
  public :: adsorbed_species_names
  public :: reactants
  public :: nr,ns,np
  public :: nu_power, nu_pr_power
  public :: mu_power, mu_pr_power
 !
!***************************************************************!
!  Particle independent variables below here                    !
!***************************************************************!
!
  real, dimension(:), allocatable :: reaction_order
  real, dimension(:), allocatable :: effectiveness_factor_old
  real, dimension(:), allocatable :: B_k, ER_k, ac, RR_method
  real, dimension(:), allocatable :: sigma_k,rho_p
  real, dimension(:), allocatable :: omega_pg_dbl
  real, dimension(:,:), allocatable :: mu, mu_prime
  real, dimension(:,:), allocatable :: nu, nu_prime
  real, dimension(:,:), allocatable :: nu_power, mu_power
  real, dimension(:,:), allocatable :: nu_pr_power, mu_pr_power
  real, dimension(:), allocatable :: aac
  real, dimension(:), allocatable :: dngas
  real, dimension(:), allocatable :: diff_coeff_reactants
  real, dimension(:,:), allocatable, save :: part_power


!
!  JONAS: implement something to calculate molar_mass of 
!  gas phase species
!
  real, dimension(nchemspec) :: molar_mass=1.0
  real, dimension(50) :: reaction_enhancement=1
  character(3), dimension(:), allocatable, save :: reaction_direction
  character(3), dimension(:), allocatable ::flags
  character(10), dimension(:,:), allocatable, save :: part
  character(10), dimension(:), allocatable :: solid_species
  character(10), dimension(50) :: species_name,adsorbed_species_names
  character(10), dimension(40) :: reactants,products
  character(20) :: element, writeformat
!
  integer :: N_species, N_reactions
  integer :: placeholder=1
  integer :: N_surface_reactions, N_adsorbed_species
  integer :: N_surface_reactants, N_surface_species
  integer :: nr, ns, np,N_max_elements
  integer :: iads=0, iads_end=0
  integer :: isurf=0,isurf_end=0
  integer :: inuH2,inuCO2,inuH2O,inuCO,inuCH4,inuO2
  integer :: inuCH,inuHCO,inuCH2,inuCH3
  integer :: imufree, imuadsO, imuadsO2, imuadsOH, imuadsH, imuadsCO
  integer, dimension(:), allocatable :: dependent_reactant
  integer, dimension(:), allocatable :: j_of_inu
  real :: St_save, Sg_save
  real :: x_s_total, rho_part
  real :: effectiveness_factor_timeaver=1.
  real :: eta_int=0.,delta_rho_surface=0.
  real :: St_first, A_p_first, rho_p_first
  real :: diffusivity = 0.0
  real :: molar_mass_carbon=12.
  real :: total_carbon_sites=1.08e-7
  real :: chemplaceholder=0.0
  real :: porosity=2.
  real :: tortuosity=3.
  real :: Init_density_part=1.300 ! g/cm^3
  real :: structural_parameter=0.05
  real :: pi_loc=3.14152
!
!  JONAS: implement something to calculate molar_mass of 
!  gas phase species
!
  logical :: first_pchem=.true.
  logical :: lpchem_debug = .false.
!
!*********************************************************************!
!             Particle dependent variables below here                 !
!*********************************************************************!
!
  real, dimension(:,:), allocatable ::  St_array
  real, dimension(:), allocatable :: conversion
  real, dimension(:), allocatable :: init_mass,rho_p_init
  real, dimension(:,:), allocatable :: mdot_ck,RR_hat
  real, dimension(:,:), allocatable :: qk_reac
  real, dimension(:), allocatable :: St
  real, dimension(:), allocatable :: A_p_init
  real, dimension(:), allocatable :: R_c_hat
  real, dimension(:,:), allocatable :: heating_k,entropy_k
  real, dimension(:), allocatable :: St_init
  real, dimension(:), allocatable :: mod_surf_area
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
  real, dimension(:,:), allocatable :: ndot,K_k
  real, dimension(:,:), allocatable :: R_j_hat
  real, dimension(:,:), allocatable :: Cs
  real, dimension(:), allocatable :: mass_trans_coeff_reactants
  real, dimension(:), allocatable :: initial_density

!
!  Some physical constants
!
  real :: mol_mass_carbon=12.0
  !real :: Sgc_init=3e5 ! m^2/kg
  real :: Sgc_init=3e4 ! cm^2/g NILS: Must set up a system for these dimensional parameters
!
!  is already in the code (R_CGS), with ergs as unit!!!
!
  real :: gas_constant=8314.0 ![J/kmol/K]
!
  namelist /particles_chem_init_pars/ &
       reaction_enhancement, &
       total_carbon_sites, &
       diffusivity,&
       porosity,&
       tortuosity,&
       structural_parameter, &
       Sgc_init, &
       lpchem_debug
!
  namelist /particles_chem_run_pars/ &
       chemplaceholder
!
  contains
!***********************************************************************
  subroutine register_particles_chem()
!
    character(10), dimension(40) :: species
!
!  wrapper routine for particle dependent and independent chemistry 
!  variables
!
          if (lroot) call svn_id( &
          "$Id: particles_chemistry.f90 20843 2014-10-06 18:45:43Z jonas.kruger $")

    call get_pchem_info(species,'N_species',N_species,'quiet')
    print*, 'Number of species in mechanics file: ', N_species
    print*, 'Species found: ', species(:N_species)
    call register_indep_pchem()    
    call register_dep_pchem()
!
  end subroutine register_particles_chem
!***********************************************************************
  subroutine register_indep_pchem()
!
    integer :: stat, i
    logical :: lenhance
!
    allocate(dngas(N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
         'Could not allocate memory for dngas')
    allocate(omega_pg_dbl(N_species),STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for omega_pg_dbl')
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
    allocate(reaction_order(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for reaction_order')
    allocate(effectiveness_factor_old(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for effectiveness_factor_old')
   allocate(nu_power(N_surface_species,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_chem',&
        'Could not allocate memory for nu_power')
   allocate(nu_pr_power(N_surface_species,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_chem',&
        'Could not allocate memory for nu_pr_power')
    allocate(mu_power(N_adsorbed_species,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_chem',&
        'Could not allocate memory for mu_power')
   allocate(mu_pr_power(N_adsorbed_species,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_chem',&
        'Could not allocate memory for mu_pr_power')
!
    effectiveness_factor_old=1.
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
      !if (lenhance) call sleep(4)
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
  end subroutine register_indep_pchem
!***********************************************************************
  subroutine register_dep_pchem()
!
    integer :: stat
!
    allocate(k_k(mpar_loc,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for k_k')
    allocate(A_p_init(mpar_loc)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for A_p_init')
    allocate(St_array(mpar_loc,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for St_array')
    allocate(St_init(mpar_loc),STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for St_init')
    allocate(init_mass(mpar_loc),STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for init_mass')
    allocate(rho_p_init(mpar_loc),STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for rho_p_init')
    allocate(mdot_ck(mpar_loc,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for mdot_ck')
!!$    allocate(f_RPM(mpar_loc)   ,STAT=stat)
!!$    if (stat>0) call fatal_error('register_dep_pchem',&
!!$        'Could not allocate memory for f_RPM')
    allocate(qk_reac(mpar_loc,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
         'Could not allocate memory for qk_reac')
    allocate(R_c_hat(mpar_loc)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
         'Could not allocate memory for R_c_hat')
    allocate(heating_k(mpar_loc,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for heating_k')
    allocate(entropy_k(mpar_loc,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for heating_k')
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
      integer :: N_trash,variable,stat
      character(10), dimension(:) :: species
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
      if (.not. allocated(part_power)) then
         allocate(part_power(N_max_elements,N_surface_reactions))
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
      call count_species(part,species,reactants,products)
      call count_species_type(species(:ns),N_adsorbed_species,&
          N_surface_species,ns)
      call count_species_type(reactants,N_trash,N_surface_reactants,nr)
!
      if (trim(string)=='N_surface_species') variable=N_surface_species
      if (trim(string)=='N_adsorbed_species') variable=N_adsorbed_species
      if (trim(string)=='N_surface_reactants') variable=N_surface_reactants
      if (trim(string)=='N_species') variable=ns
      if (trim(string)=='N_reactants') variable=nr
!
    end subroutine get_pchem_info
!***********************************************************************
    subroutine calc_R_c_hat()
! 
!
!  JONAS talk to nils how to implement things from equ.f90
!  JONAS implement right sequence so that mdot_ck is always calculated
!  first
!
      R_c_hat = 0.0
      R_c_hat(:) = -sum(mdot_ck,DIM=2)/(St(:)*mol_mass_carbon)
!
    end subroutine calc_R_c_hat
!***********************************************************************
    subroutine calc_mod_surf_area(fp)
!
      real, dimension(:), allocatable :: mod_all,Sgc
      real, dimension(:,:) :: fp
      integer :: k,k1,k2
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
    do k=k1,k2
       mod_all(k) = St_init(k)*St_init(k)*structural_parameter* &
            (1-conversion(k)) * (1-conversion(k)) / &
            (2*St(k)**2)
       Sgc(k) = St(k)/ rho_p(k)      
       mod_surf_area(k) = (1-mod_all(k))*Sgc(k)*mol_mass_carbon
    enddo
!
    deallocate(mod_all)
    deallocate(Sgc)
!
    end subroutine calc_mod_surf_area
!***********************************************************************
    subroutine calc_St(fp)
!
      real, dimension(:,:) :: fp
      integer :: k,k1,k2
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
!  Evolution of the total surface area according to 8th US combustion
!  meeting, coal and biomass combustion and gasification eq. 19
!
      St = 0.0
!      
      do k=k1,k2
         St(k)=fp(k,imp)*Sgc_init* &
              sqrt(1.0 - structural_parameter*log(rho_p(k)/rho_p_init(k)))
      enddo
!
   end subroutine calc_St
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
  subroutine create_stoc(part,list,targ,lhs,nlist,power)
!
  integer :: i,j,k,stat,nlist
  real :: multi
  character(10), dimension(:,:) :: part
  character(10), dimension(:) :: list
  character(10) :: element
  real, dimension(:,:) :: targ,power
  logical :: lhs,fwd,forpower
!
!  list where the stochiometry is saved in
!
  if (lhs) then 
     forpower=.true. 
  else
     forpower=.false.
  endif
  power = 0.0
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
              if (forpower) then 
                 power(k,i) = part_power(j,i)
              else
                 power(k,i) = 1.0
              end if
           else
           end if
        end do
      end do
   end do
!   targ(:,:) = int(targ(:,:))
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
  subroutine count_species(part,species,reactants,products)
!
    character(10), dimension(:,:) :: part
    character(10) :: element
    real :: numeric
    integer :: i,j,jmax,stat,place,number
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
       call remove_save_powers(target_list)
!
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
    character(10), dimension(40) :: list
!!$    character(10), dimension(40) :: species,reactants
!!$    integer :: stat
!!$!
!!$    call get_pchem_info(species,'N_adsorbed_species',N_adsorbed_species,'quiet')
!!$    call get_pchem_info(species,'N_surface_species',N_surface_species,'quiet')
!!$
!!$    call create_ad_sol_lists(species,adsorbed_species_names,'ad',ns)
!!$    call sort_compounds(reactants,adsorbed_species_names,N_adsorbed_species,nr)
!!$!
!!$    allocate(solid_species(N_surface_species)   ,STAT=stat)
!!$    if (stat>0) call fatal_error('register_indep_pchem',&
!!$        'Could not allocate memory for solid_species')
!!$!
!!$    call create_ad_sol_lists(species,solid_species,'sol',ns)
!!$    call sort_compounds(reactants,solid_species,N_surface_species,nr)
!!$!
!!$    if (trim(string)=='solid_species') &
!!$         list(1:N_surface_species)=solid_species(:)
!!$    if (trim(string)=='adsorbed_species_names') &
!!$         list(1:N_adsorbed_species)=adsorbed_species_names(:N_adsorbed_species)
!!$    if (trim(string)=='all') then
!!$       list(1:N_surface_species)=solid_species(:)
!!$       list(N_surface_species+1:N_surface_species+N_adsorbed_species)= &
!!$            adsorbed_species_names(:N_adsorbed_species)
!!$    else
!!$    endif
!!$!
!!$    deallocate(solid_species)
!!$!
  end subroutine get_species_list
!**********************************************************************
  subroutine create_dngas()
!
!  Find the mole production of the forward reaction. This will later
!  be used for the calculation of the reverse reaction rate.
!
    integer :: k
!
    do k=1,N_surface_reactions
       dngas(k) = sum(nu_prime(:,k)) - sum(nu(:,k))
    enddo
!
  end subroutine create_dngas
!**********************************************************************
 subroutine calc_entropy_of_reaction()

   integer :: i,j,l,k,k1,k2
!
    entropy_k=0
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    do k=k1,k2
    do l=1,N_surface_reactions
      do i=1,N_surface_species
        entropy_k(k,l)=entropy_k(k,l)&
            +nu_prime(i,l)*surface_species_entropy(k,i)&
            -nu(i,l)*surface_species_entropy(k,i)
      enddo
      if (N_adsorbed_species > 0) then
      do j=1,N_adsorbed_species
        entropy_k(k,l)=entropy_k(k,l)&
            +mu_prime(j,l)*adsorbed_species_entropy(k,j)&
            -mu(j,l)*adsorbed_species_entropy(k,j)
      enddo
      else
      end if
    enddo
    enddo
!
  end subroutine calc_entropy_of_reaction
!**********************************************************************
  subroutine get_reverse_K_k(l,K_k,fp)
!
    use Particles_cdata, only: interp_pp
!
    integer :: l,k,k1,k2
    real, dimension(:), allocatable :: k_p,k_c
    real, dimension(:), allocatable :: denominator, exponent
    real, dimension(:,:) :: fp
    real, dimension(:,:) :: K_k
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    do k=k1,k2
       denominator(l) = heating_k(k,l-1) - (entropy_k(k,l-1)*fp(k,iTp))
       exponent(k) = denominator(k)/(gas_constant*fp(k,iTp))
       k_p(k) = exp(-exponent(k))
       k_c(k) = k_p(k) / (((gas_constant)*fp(k,iTp)/interp_pp(k))**(dngas(l-1)))
       K_k(k,l) = (K_k(k,l-1) / k_c(k))
    enddo
!
  end subroutine get_reverse_K_k
!**********************************************************************
  subroutine calc_conversion(fp)
!
    real, dimension(:,:) :: fp
    integer :: k
!
    do k=k1_imn(imn),k2_imn(imn)
      conversion(k) = fp(k,imp) / init_mass(k)
    enddo
!
  end subroutine calc_conversion
!***********************************************************************
  subroutine get_St(var,start,end)
!
    real, dimension(:) :: var
    integer :: start,end
!
    var = St(start:end)
!
  end subroutine get_St
!***********************************************************************
  subroutine get_R_c_hat(var,start,end)
!
    real, dimension(:) :: var
    integer :: start, end
!
    var = R_c_hat(start:end)
!
  end subroutine get_R_c_hat
!***********************************************************************
  subroutine get_mod_surf_area(var,start,end)
!
    real, dimension(:) :: var
    integer :: start, end
!
    var = mod_surf_area(start:end)
!
  end subroutine get_mod_surf_area
!***********************************************************************
  subroutine get_conversion(var,start,end)
!
    real, dimension(:) :: var
    integer :: start,end
!
    var = conversion(start:end)
!
  end subroutine get_conversion
!***********************************************************************  
  subroutine calc_enthalpy_of_reaction()
!
    integer :: i,k,k1,k2,j,l
!
    heating_k=0
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    do k=k1,k2
    do l=1,N_surface_reactions
      do i=1,N_surface_species
        heating_k(k,l)=heating_k(k,l)&
            +nu_prime(i,l)*surface_species_enthalpy(k,i)&
            -nu(i,l)*surface_species_enthalpy(k,i)
      enddo
      if (N_adsorbed_species > 0) then
      do j=1,N_adsorbed_species
        heating_k(k,l)=heating_k(k,l)&
            +mu_prime(j,l)*adsorbed_species_enthalpy(k,j)&
            -mu(j,l)*adsorbed_species_enthalpy(k,j)
      enddo
      else
      end if
    enddo
    enddo
!
  end subroutine calc_enthalpy_of_reaction
!**********************************************************************
  subroutine calc_RR_hat(f,fp)
!
!  01-oct-2014/jonas:coded
!
!  JONAS: needs to be filled with life
!  solid_reac L:30
!  needed: local pressure UNITS!!!!
!  needed: Cg, Cs, thiele
!
    real, dimension(mpar_loc,mpvar) :: fp
    real, dimension(mx,my,mz,mfarray) :: f
    integer :: i,j,k,k1,k2

    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
!!$    do k=k1,k2
!!$!
!!$          interp_pp(k)
!!$    do j=1,N_surface_reactions
!!$    do i=1,N_surface_reactants
!!$          if (nu(i,j) > 0) RR_hat(k,j)=RR_hat(k,j)*&
!!$!
!!$!  old formula            (Press*x_surface(i))*
!!$!  new                  (Cg*fp(j,isurf-1+i))*
!!$!
!!$          1**nu(i,j)/molar_mass_carbon
!!$    enddo
!!$       adsorbed: if (ladsorbed_species) then
!!$       do i=1,N_adsorbed_species
!!$          if(mu(i,j)> 0) RR_hat(k,j)=RR_hat(k,j)*(Cs(i))**mu(i,j)
!!$       enddo
!!$       else adsorbed
!!$       endif adsorbed
!!$    enddo
!!$    enddo
!!$!
!!$    if (lthiele) then
!!$       calc_effectiveness_factor(effectiveness_factor,fp)
!!$       do j=1,N_surface_reactions
!!$       RR_hat(:,j) = RR_hat(:,j) * effectiveness_factor(:)
!!$       enddo
!!$    else
!!$    endif
!
    RR_hat = 1.0
!
  end subroutine calc_RR_hat
!********************************************************************
  subroutine calc_ndot_mdot_R_j_hat(fp)
!
!  taken from solid_reac L. 108 ff 
!
! Then find the carbon consumption rate and molar flux of each gaseous 
! species at the particle surface
!
    real, dimension(:,:) :: fp
    integer :: n,i,j,k,l,k1,k2
!

    ndot=0    
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
    Rck=0
!
    do k=k1,k2
       do l=1,N_surface_reactions
          do i=1,N_surface_species
             Rck(k,l)=Rck(k,l)+molar_mass_carbon*RR_hat(k,l)*(nu_prime(i,l)-nu(i,l))*ac(i)
             ndot(k,i)=ndot(k,i)+RR_hat(k,l)*(nu_prime(i,l)-nu(i,l))*St_array(k,l)/ &
                  (fp(k,iap)*fp(k,iap)*4.*pi_loc)
          enddo
       enddo
    enddo
!
! Find molar reaction rate of adsorbed surface species
!
  if (N_adsorbed_species>1) then
    R_j_hat=0
    do k=k1,k2
    do l=1,N_surface_reactions
      do j=1,N_adsorbed_species-1
        R_j_hat(k,j)=R_j_hat(k,j)+(mu_prime(j,l)-mu(j,l))*RR_hat(k,l)
        Rck(k,l)=Rck(k,l)+molar_mass_carbon*RR_hat(k,l)*(mu_prime(j,l)-mu(j,l))*aac(j)
      enddo
    enddo
    enddo
  endif
!
! Find mdot_ck
!
  do k=k1,k2
  do l=1,N_surface_reactions
    mdot_ck(k,l)=-St_array(k,l)*Rck(k,l)
  enddo
  enddo
!
! Find the sum of all molar fluxes at the surface
!
  do k=k1,k2
     ndot_total(k)=sum(ndot(k,:))
  enddo
!
  end subroutine calc_ndot_mdot_R_j_hat
!**********************************************************************
  subroutine get_RR_hat(var,start,end)
!
    real, dimension(:,:) :: var
    integer :: start, end
!
    var(start:end,:) = RR_hat(start:end,:)
!
  end subroutine get_RR_hat
!**********************************************************************
  subroutine get_part(var)
!
    character(10), dimension(:,:) :: var
!
    var = part
!
  end subroutine get_part
!**********************************************************************
  subroutine get_reactants(var)
!
    character(10), dimension(40) :: var
!
    intent(out) :: var
!
    var = reactants
!
  end subroutine get_reactants
!**********************************************************************
  subroutine get_total_carbon_sites(var)
!
    real :: var
!
    var = total_carbon_sites
!
  end subroutine get_total_carbon_sites
!**********************************************************************
  subroutine calc_effectiveness_factor(fp,var)
!
!  01-Oct-2014/Jonas: coded
!  taken from solid_reac L. 149 and equations.pdf eq 35ff
!
    real, dimension(:) :: var
    real, dimension(:,:) :: fp
!
!    var = 0.0
!
    real, dimension(:,:), allocatable :: R_i_hat,D_eff
    real, dimension(:,:), allocatable :: effectiveness_factor_species
    real, dimension(:), allocatable :: Knudsen, pore_radius
    real, dimension(:), allocatable :: tmp1,tmp2,tmp3
    real, dimension(:), allocatable ::  phi,sum_eta_i_R_i_hat_max
    real, dimension(:), allocatable :: sum_R_i_hat_max
    real, dimension(:), allocatable :: volume
    integer :: i,dep,n,k,k1,k2,l
    real :: r_f,Cg
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    allocate(R_i_hat(k1:k2,N_surface_reactants))
    allocate(D_eff(k1:k2,N_surface_reactants))
    allocate(volume(k1:k2))
    allocate(pore_radius(k1:k2))
    allocate(Knudsen(k1:k2))
    allocate(tmp1(k1:k2))
    allocate(phi(k1:k2))
    allocate(tmp2(k1:k2))
    allocate(tmp3(k1:k2))
    allocate(sum_R_i_hat_max(k1:k2))
!
!  Find reaction rate for each of the reactants
!
    R_i_hat=0
    do k=k1,k2
    do l=1,N_surface_reactions
       do i=1,N_surface_reactants
          R_i_hat(k,i)=R_i_hat(k,i)+(nu_prime(i,l)-nu(i,l))*RR_hat(k,l)
       enddo
    enddo
    enddo
!
!  Find particle volume
!
    volume(:)=4.*pi_loc*(fp(k1:k2,iap)**3)/3.
!
! Find pore radius
!
    r_f=2.
    pore_radius(:)=2*r_f*porosity*volume(:)/St(k1:k2)
!
!  Find effective diffusion coefficient (based on bulk and Knudsen diffusion)
!
!  JONAS: check which attributes change from particle to particle!!
!
    do k=k1,k2
    do i=1,N_surface_reactants
       tmp3(k)=8*gas_constant*fp(k,iTp)/(pi_loc*molar_mass(j_of_inu(i)))
       Knudsen(k)=2*pore_radius(k)*porosity*sqrt(tmp3(k))/(3*tortuosity)
       tmp1(i)=1./diff_coeff_reactants(i)
       tmp2(k)=1./Knudsen(k)
       D_eff(k,i)=1./(tmp1(i)+tmp2(k))
    enddo
    enddo
!
!  Find thiele modulus and effectiveness factor
!  JONAS: was with volume(k) before,
!
    do k=k1,k2
       do i=1,N_surface_reactants
          if (R_i_hat(k,i)<0) then
             thiele(k,i)=fp(k,iap)*sqrt(-R_i_hat(k,i)*St(k)/&
                  (fp(k,imp)*Cg*fp(k,iads-1+i)*D_eff(k,i)))
          else
             thiele(k,i)=1e-2
          endif
!
!  calculate the effectivenes factor (all particles, all reactants)
!
          effectiveness_factor_species(k,i)=3/thiele(k,i)* &
               (1./tanh(thiele(k,i))-1./thiele(k,i))
       enddo
    enddo
!
!  The mean effectiveness factor can be found from the above by using R_i_max
!  (all particles, all reactants) 
!
    sum_R_i_hat_max=0.0
!             
    do k=k1,k2
       do i=1,N_surface_reactants
          if (R_i_hat(k,i)<0) then
             sum_R_i_hat_max(k)=sum_R_i_hat_max(k)+R_i_hat(k,i)
          endif
       enddo
    enddo
!
    do k=k1,k2
       sum_eta_i_R_i_hat_max(n)=0
       do i=1,N_surface_reactants
          if (R_i_hat(k,i)<0) then
             sum_eta_i_R_i_hat_max(k)=sum_eta_i_R_i_hat_max(k)&
                  +R_i_hat(k,i)*effectiveness_factor_species(k,i)
          endif
       enddo
    enddo
!
    effectiveness_factor(:)=sum_eta_i_R_i_hat_max(:)/sum_R_i_hat_max(:)
!
!  The efficiency factor for each reaction can now be found from the thiele modulus
!
!  JONAS: the iteration parameters have switched from solid_reac in order to 
!  be consistent! be careful when checking!
!
    do k=k1,k2
    do l=1,N_surface_reactions
       dep=dependent_reactant(l)
       if (dep>0) then
          if (R_i_hat(k,dep)>0.) then
             effectiveness_factor_reaction(k,l)=1.
          else
             phi(k)=thiele(k,dependent_reactant(l))
             effectiveness_factor_reaction(k,l)=3* &
                  (1./tanh(phi(k))-1./phi(k))/phi(k)
          endif
       else
          effectiveness_factor_reaction(k,l)=1
       endif
    enddo
    enddo
!
    deallocate(R_i_hat)
    deallocate(D_eff)
    deallocate(volume)
    deallocate(pore_radius)
    deallocate(Knudsen)
    deallocate(tmp1)
    deallocate(tmp2)
    deallocate(tmp3)
    deallocate(sum_R_i_hat_max)
    deallocate(phi)
!
  end subroutine calc_effectiveness_factor
!**********************************************************************
    subroutine calc_surf_enthalpy(fp)
!
    real, dimension(mpar_loc,mpvar) :: fp    
    integer :: k,k1,k2
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
!  Unit: J/kmolK
!
    do k=k1,k2
       if (inuH2O>0) surface_species_enthalpy(k,inuH2O)=-242.18e6-(5.47e3*fp(k,iTp))
       if (inuO2>0)  surface_species_enthalpy(k,inuO2) = 0.
       if (inuCO2>0) surface_species_enthalpy(k,inuCO2)=-392.52e6-(2.109e3*fp(k,iTp))
       if (inuH2>0)  surface_species_enthalpy(k,inuH2 )= 0.
       if (inuCO>0)  surface_species_enthalpy(k,inuCO )=-105.95e6-6.143e3*fp(k,iTp)
       if (inuCH>0)  surface_species_enthalpy(k,inuCH )= 594.13e6
       if (inuHCO>0) surface_species_enthalpy(k,inuHCO)=  45.31e6-(5.94e3*fp(k,iTp))
       if (inuCH2>0) surface_species_enthalpy(k,inuCH2)= 387.93e6-(5.8e3*fp(k,iTp))
       if (inuCH4>0) surface_species_enthalpy(k,inuCH4)= -75e6
       if (inuCH3>0) surface_species_enthalpy(k,inuCH3)= 144.65e6-(6.79e3*fp(k,iTp))
    enddo
!
    end subroutine calc_surf_enthalpy
!**********************************************************************
      subroutine calc_surf_entropy(fp)
!
    real, dimension(:,:) :: fp
    integer :: k,k1,k2
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
! JONAS: units are in j/(kmol*K)
!
    do k=k1,k2
    if (inuH2O>0) surface_species_entropy(k,inuH2O)= & 
         189.00e3+(0.0425e3*fp(k,iTp))
    if (inuO2>0)  surface_species_entropy(k,inuO2) = &
         222.55e3+(0.0219e3*fp(k,iTp))
    if (inuCO2>0) surface_species_entropy(k,inuCO2)= &
         212.19e3+(0.0556e3*fp(k,iTp))
    if (inuH2>0)  surface_species_entropy(k,inuH2 )= &
         133.80e3+(0.0319e3*fp(k,iTp))
    if (inuCO>0)  surface_species_entropy(k,inuCO )= &
         199.35e3+(0.0342e3*fp(k,iTp))
!
!  taken from chemistry  webbook (1bar)
!
    if (inuCH>0)  surface_species_entropy(k,inuCH )=183.00e3
    if (inuHCO>0) surface_species_entropy(k,inuHCO)=223.114e3+(0.0491e3*fp(k,iTp))
    if (inuCH2>0) surface_species_entropy(k,inuCH2)=193.297e3+(0.0467e3*fp(k,iTp))
!
!  taken from chemistry webbook (1bar)
!
    if (inuCH4>0) surface_species_entropy(k,inuCH4)= 189.00e3
    if (inuCH3>0) surface_species_entropy(k,inuCH3)= 190.18e3+(0.0601e3*fp(k,iTp))
!
    enddo
!
    end subroutine calc_surf_entropy
!*********************************************************************
  subroutine calc_ads_entropy(fp)
!
    real, dimension(mpar_loc,mpvar) :: fp
    integer :: k,k1,k2
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
!  Unit: J/kmolK
!
    adsloop: do k=k1,k2
!
    if (imuadsO>0)    then
       adsorbed_species_entropy(k,imuadsO) = &
    (164.19e3+(0.0218e3*fp(k,iTp)))*0.72 - (3.3*gas_constant)
    else
    end if
!
!  this is guessed
!
    if (imuadsO2>0)   then
       adsorbed_species_entropy(k,imuadsO2) =  &
            2*adsorbed_species_entropy(k,imuadsO)
    else
    end if
    if (imuadsOH>0)   then
       adsorbed_species_entropy(k,imuadsOH) = &
         ((0.0319e3*fp(k,iTp)) + 186.88e3) * 0.7 - (3.3*gas_constant)
    else
    end if
    if (imuadsH>0)    then 
       adsorbed_species_entropy(k,imuadsH) = &
           (117.49e3+(0.0217e3*fp(k,iTp)))*0.54 - (3.3*gas_constant)
    else
    end if
    if (imuadsCO>0)   then
       adsorbed_species_entropy(k,imuadsCO) = &
            surface_species_entropy(k,inuCO)* &
            0.6*(1+(1.44e-4*fp(k,iTp))) - (3.3*gas_constant)
    else
    end if
!
!  taken from nist
!
    if (imufree>0)    adsorbed_species_entropy(k,imufree) = 0
!
    enddo adsloop
!
  end subroutine calc_ads_entropy
!*********************************************************************
  subroutine calc_ads_enthalpy(fp)
!
    real, dimension(:,:) :: fp
    integer :: k,k1,k2
!
!  JONAS: values are from nist and solid_phase.f90 of the stanford
!  code Units: J/kmol
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    if (imuadsO>0)    adsorbed_species_enthalpy(k1:k2,imuadsO) = &
         -148.14e6 + (0.0024e6*(fp(k1:k2,iTp)-273.15))
    if (imuadsO2>0)   adsorbed_species_enthalpy(k1:k2,imuadsO2)= &
         2 *  (-148.14e6 + (0.0024e6*(fp(k1:k2,iTp)-273.15)))
    if (imuadsOH>0)   adsorbed_species_enthalpy(k1:k2,imuadsOH)= &
         -148e6
    if (imuadsH>0)    adsorbed_species_enthalpy(k1:k2,imuadsH) = &
         -19.5e6
    if (imuadsCO>0)   adsorbed_species_enthalpy(k1:k2,imuadsCO)= &
         -199.94e6 - (0.0167e6*(fp(k1:k2,iTp)-273.15))
    if (imufree>0)    adsorbed_species_enthalpy(k1:k2,imufree) = 0.
!
  end subroutine calc_ads_enthalpy
!*********************************************************************
  subroutine calc_mass_init(fp)
!
!  21-oct-2014/nils: coded
!
    real, dimension(:,:) :: fp
    integer :: k
!
    do k=1,mpar_loc
      init_mass(k)=fp(k,imp)
    enddo
!
  end subroutine calc_mass_init
!*********************************************************************
  subroutine calc_St_init(fp)
!
!  3-oct-2014/jonas:coded
!  this subroutine calculates the initial total surface area
!
    real, dimension(:,:) :: fp
    integer :: k
!    real :: rel_part_dens
!
!    k1=k1_imn(imn)
!    k2=k2_imn(imn)
!
!    allocate(f_RPM(k1:k2))
!
!  JONAS: there was a complicated formula in the code. i leave it
!  just commented out for now.
!
   do k=1,mpar_loc
!
!  "easy" way:
      St_init(k) = Sgc_init * fp(k,imp)
!
!!$      rel_part_dens = rho_p_init(k) / Init_density_part
!!$      f_RPM(k) = sqrt(1-structural_parameter*log(rel_part_dens))
!!$      St_init(k) = 4* fp(k,iap)**2 *pi_loc * f_RPM(k)
   enddo
!
!   deallocate(f_RPM)
!
!
!!$    if(lpchem_debug) then
!!$       call print_debug_info()
!!$    else
!!$    endif
!
  end subroutine calc_St_init
!*********************************************************************
  subroutine calc_chemistry_pencils(f,fp)
!
    real, dimension(mpar_loc,mpvar) :: fp
    real, dimension(mx,my,mz,mfarray) :: f
    integer :: stat
!
!  Routine to calcute quantities used for reactive particles
!
    call allocate_variable_pencils()
!   
    call calc_rho_p(fp)
    call calc_conversion(fp)
    call calc_St(fp)
    call calc_mod_surf_area(fp)
!
    call calc_ads_entropy(fp)
    call calc_ads_enthalpy(fp)
    call calc_surf_entropy(fp)
    call calc_surf_enthalpy(fp)
!
    call calc_enthalpy_of_reaction()
    call calc_entropy_of_reaction()
!
    call calc_RR_hat(f,fp)
    call calc_ndot_mdot_R_j_hat(fp)
    call calc_R_c_hat()
!
  end subroutine calc_chemistry_pencils
!***************************************************
  subroutine allocate_variable_pencils()
!
!
!
    integer :: k1,k2,stat
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    allocate(rho_p(k1:k2) ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for rho_p')
    allocate(St(k1:k2) ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for St')
    allocate(mod_surf_area(k1:k2) ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for St')
!
    allocate(conversion(k1:k2) ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for conversion')
!
    allocate(surface_species_enthalpy(k1:k2,N_surface_species) &
         ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for surface_species_enthalpy')
    allocate(surface_species_entropy(k1:k2,N_surface_species) &
         ,STAT=stat)
    allocate(adsorbed_species_entropy(k1:k2,N_adsorbed_species) &
         ,STAT=stat)    
    if (stat>0) call fatal_error('register_dep_pchem',&
         'Could not allocate memory for adsorbed_species_entropy')
    allocate(adsorbed_species_enthalpy(k1:k2,N_adsorbed_species) &
         ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
         'Could not allocate memory for adsorbed_species_enthalpy')
!
    allocate(thiele(k1:k2,N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for thiele')
!
    allocate(effectiveness_factor(k1:k2))
    allocate(effectiveness_factor_species(k1:k2,N_surface_reactants))
    allocate(effectiveness_factor_reaction(k1:k2,N_surface_reactions))
!
    allocate(ndot(k1:k2,N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for ndot')
    allocate(ndot_total(k1:k2)   ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for ndot_total')
!
    allocate(Rck(k1:k2,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for Rck')
    allocate(Rck_max(k1:k2,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for Rck_max')
!
    allocate(R_j_hat(k1:k2,N_adsorbed_species-1)   ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
           'Could not allocate memory for R_j_hat')
      allocate(Cs(k1:k2,N_adsorbed_species)   ,STAT=stat)
      if (stat>0) call fatal_error('allocate_variable_pencils',&
           'Could not allocate memory for Cs')
      allocate(RR_hat(k1:k2,N_surface_reactions)   ,STAT=stat)
      if (stat>0) call fatal_error('allocate_variable_pencils',&
           'Could not allocate memory for RR_hat')

!
  end subroutine allocate_variable_pencils
!***************************************************
  subroutine cleanup_chemistry_pencils()
!
!  06-oct-14/jonas:coded
!
    if (lparticles_chemistry) then
    deallocate(mod_surf_area)
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
    endif
!
  end subroutine cleanup_chemistry_pencils
!***************************************************
    subroutine read_particles_chem_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_chem_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_chem_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_chem_init_pars
!***********************************************************************
    subroutine write_particles_chem_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_chem_init_pars)
!
    endsubroutine write_particles_chem_init_pars
!***********************************************************************
    subroutine read_particles_chem_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_chem_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_chem_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_chem_run_pars
!***********************************************************************
    subroutine write_particles_chem_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_chem_run_pars)
!
    endsubroutine write_particles_chem_run_pars
!***********************************************************************
    subroutine calc_rho_p_init(fp)
!
!  07-oct-14/jonas:coded
!
!  now in g/cm^3!!!!
!
      real, dimension(:,:) :: fp
!
      rho_p_init(:) = fp(:,imp) / (4.*pi_loc/3. * fp(:,iap)**3)
!
    end subroutine calc_rho_p_init
!***********************************************************************
    subroutine calc_rho_p(fp)
!
!  07-oct-14/jonas:coded
!
      real, dimension(:,:) :: fp
      integer :: k,k1,k2
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      do k=k1,k2
         rho_p(k) = fp(k,imp) / (4./3. *pi_loc * fp(k,iap)**3)
      enddo
!
    end subroutine calc_rho_p
!***********************************************************************
  subroutine remove_save_powers(target_list)
!
!  10-oct-14/jonas:coded
!
    character(10), dimension(:,:) :: target_list
    character(10) :: element, writeformat
    real :: power
    integer :: i,j,pow_pla
    integer :: isreal
!
    writeformat = '(  F6.2)'
    write(writeformat(2:3),'(I2)') N_max_elements
!
    i = 1
    do while (i<=N_surface_reactions)
       j = 1
       do while (j <= N_max_elements)
          element = target_list(j,i)
          if (index(element,'^') > 0) then
             pow_pla = index(element,'^')
             read(element(pow_pla+1:pow_pla+4),*,iostat=isreal) power
             target_list(j,i) = element(1:pow_pla-1)
             if (isreal==0) then
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
!
  end subroutine remove_save_powers
!**********************************************************************
  subroutine print_debug_info()
!
!
!
    integer :: k
!
    print*,'Solid_species'
    print*, solid_species
    writeformat = '(A12," ",I4,  F7.2)'
    write(writeformat(13:14),'(I2)') N_surface_species
!
    do k=1,N_surface_reactions
       write(*,writeformat) 'nu=',k,nu(:,k)
    end do

    do k=1,N_surface_reactions
       write(*,writeformat) 'nu_power=',k,nu_power(:,k)
    end do
!
    do k=1,N_surface_reactions
       write(*,writeformat) 'nu_prime=',k,nu_prime(:,k)
    end do
!
    write(writeformat(13:14),'(I2)') N_adsorbed_species
!
    print*,'Adsorbed_species_names'
    print*, adsorbed_species_names
    do k=1,N_surface_reactions
       write(*,writeformat) 'mu=',k,mu(:,k)
    end do 
!
    do k=1,N_surface_reactions
       write(*,writeformat) 'mu_power=',k,mu_power(:,k)
    end do
 !
    do k=1,N_surface_reactions
       write(*,writeformat) 'mu_prime=',k,mu_prime(:,k)
    end do
!
    do k=1,N_surface_reactions
       write(*,'(A12,I4,2E12.5)') 'ER_k, B_k=',k,B_k(k),ER_k(k)/ &
               (1e6/gas_constant)
    enddo
    do k=1,N_surface_reactions
       write(*,'(A12,I4,E12.5)') 'Dngas',k,dngas(k)
    end do
!        
        write(*,'(A20," ",10I4)') 'j_of_inu=', j_of_inu
        write(*,'(A20," ",10F4.0)') 'ac=',ac
        write(*,'(A20," ",10F4.0)') 'site_occupancy=',aac
        write(*,'(A20," ",30I4)') 'dependent_reactant=',dependent_reactant
        write(*,'(A20," ",10E12.5)') 'sigma_k=',sigma_k
!
  end subroutine print_debug_info
!*******************************************************************
  end module Particles_chemistry
