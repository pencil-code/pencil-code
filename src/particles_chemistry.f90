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
  use Particles_radius
  use Particles_map
  use EquationOfState
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
  public :: jmap
  public :: solid_species
  public :: mass_trans_coeff_reactants
  public :: mass_trans_coeff_species
  public :: part
  public :: diff_coeff_reactants
  public :: adsorbed_species_names
  public :: reactants
  public :: nr,ns,np
  public :: nu_power
  public :: mu_power
  public :: species
  public :: dummy
  public :: lboundary_explicit
  public :: linfinite_diffusion
  public :: gas_constant
  public :: total_carbon_sites
  public :: R_c_hat
  public :: mod_surf_area
  public :: K_k, init_mass
  public :: mdot_ck
  public :: ndot_total
  public :: ndot
  public :: porosity
  public :: Cg
  public :: lreactive_heating
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
  real, dimension(:,:), allocatable :: nu_power, mu_power
  real, dimension(:), allocatable :: aac, T_k
  real, dimension(:), allocatable :: dngas
  real, dimension(:), allocatable :: diff_coeff_reactants
  real, dimension(:,:), allocatable, save :: part_power
  real, dimension(2,2) :: dummy

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
  character(10), dimension(40) :: species

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
  integer :: iter
  integer, dimension(:), allocatable :: dependent_reactant
  integer, dimension(:), allocatable :: jmap
  real :: St_save, Sg_save
  real :: x_s_total, rho_part
  real :: effectiveness_factor_timeaver=1.
  real :: eta_int=0.,delta_rho_surface=0.
  real :: St_first, A_p_first, rho_p_first
  real :: diffusivity = 0.0
  real :: total_carbon_sites=1.08e-8 ! [mol/cm^2]
  real :: chemplaceholder=0.0
  real :: tortuosity=3.
  real :: Init_density_part=1.300 ! g/cm^3
  real :: true_density_carbon=1.800 ! g/cm^3
  real :: structural_parameter=8. ! [-]
  real :: startup_time=0.
  real :: startup_quench
  real :: init_mass
!
!  JONAS: implement something to calculate molar_mass of 
!  gas phase species
!
  logical :: first_pchem=.true.
  logical :: lpchem_debug = .false.
  logical :: lthiele=.false.
  logical :: lboundary_explicit=.false.
  logical :: linfinite_diffusion=.false.
  logical :: lreactive_heating=.false.
!
!*********************************************************************!
!             Particle dependent variables below here                 !
!*********************************************************************!
!
  real, dimension(:), allocatable :: conversion
  real, dimension(:), allocatable :: rho_p_init
  real, dimension(:,:), allocatable :: mdot_ck,RR_hat
  real, dimension(:,:), allocatable :: qk_reac
  real, dimension(:), allocatable :: St,rho_p,porosity
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
  real, dimension(:,:), allocatable :: mass_trans_coeff_reactants
  real, dimension(:,:), allocatable :: mass_trans_coeff_species
  real, dimension(:), allocatable :: initial_density
  real, dimension(:), allocatable :: diff_coeffs_species,Cg,A_p
  real, dimension(:,:), allocatable :: x_surf
  real, dimension(:), allocatable :: q_reac
  real, dimension(:), allocatable :: Nu_p

!
!  Some physical constants
!
  real :: mol_mass_carbon=12.0
  !real :: Sgc_init=3e5 ! m^2/kg
  real :: Sgc_init=3e4 ! cm^2/g NILS: Must set up a system for these dimensional parameters
!
!  is already in the code (R_CGS), with ergs as unit!!!
!
  real :: gas_constant=8.314 ![J/mol/K]
!
  namelist /particles_chem_init_pars/ &
       reaction_enhancement, &
       total_carbon_sites, &
       diffusivity,&
       tortuosity,&
       structural_parameter, &
       Sgc_init, &
       lpchem_debug, &
       true_density_carbon, &
       startup_time
!
  namelist /particles_chem_run_pars/ &
       chemplaceholder, &
       lthiele, &
       lboundary_explicit, &
       linfinite_diffusion, &
       lreactive_heating
!
  contains
!***********************************************************************
  subroutine register_particles_chem()
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
    allocate(reaction_order(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for reaction_order')
    allocate(effectiveness_factor_old(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_pchem',&
        'Could not allocate memory for effectiveness_factor_old')
   allocate(nu_power(N_surface_species,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_chem',&
        'Could not allocate memory for nu_power')
    allocate(mu_power(N_adsorbed_species,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_chem',&
        'Could not allocate memory for mu_power')
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
  end subroutine register_indep_pchem
!***********************************************************************
  subroutine register_dep_pchem()
!
    integer :: stat
!
    allocate(A_p_init(mpar_loc)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for A_p_init')
    allocate(St_init(mpar_loc),STAT=stat)
    if (stat>0) call fatal_error('register_dep_pchem',&
        'Could not allocate memory for St_init')
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
      endif
      if (.not. allocated(part_power)) then
         allocate(part_power(N_max_elements,N_surface_reactions))
      end if
      if (.not. allocated(reaction_direction)) then
         allocate(reaction_direction(N_surface_reactions))
      end if
      if (.not. allocated(flags)) then
         allocate(flags(N_surface_reactions),STAT=stat)
      end if
      if (.not. allocated(T_k)) then
         allocate(T_k(N_surface_reactions)   ,STAT=stat)
      endif
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
  subroutine create_arh_param(part,B_k,ER_k,sigma_k)
!
!takes the first numerical in part and writes it to b_k
!
    character(10), dimension(:,:) :: part
    real, dimension(:) :: B_k,ER_k,sigma_k
    character(10) :: el_B_k,el_ER_k,el_sigma
    real :: B_k_single,ER_k_single,sigma_single
    logical :: done
    integer :: i,k,stat,stat1,stat2
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
!!$          k=3
!!$          print*,part(:,i)
!!$          do while (k < size(part,1)-2)
!!$             read(part(k,i),*,iostat=stat) B_k_single
!!$             read(part(k+1,i),*,iostat=stat1) ER_k_single
!!$             read(part(k+2,i),*,iostat=stat2) sigma_single
!!$             if (stat==0 .and. stat1==0 .and. stat2==0) then
!!$                B_k(i) = B_k_single
!!$                ER_k(i) = ER_k_single
!!$                sigma_k(i) = sigma_single
!!$                k=size(part,1)
!!$             else
!!$                k=k+1
!!$             endif
!!$          enddo
!!$          if (B_k(i)==0.0) call fatal_error('create_arh_param','error in line')
!!$       endif
!!$    enddo
                
       done = .false.
       do k=1, size(part,1)-2
          el_B_k = part(k,i)
          el_ER_k = part(k+1,i)
          el_sigma = part(k+2,i)
             read(el_B_k,*,iostat=stat) B_k_single
             !print*,B_k_single
             if (stat == 0 .and. (done .eqv. .false.)) then
                B_k(i) = B_k_single
                done = .true.
                read(el_ER_k,*,iostat=stat) ER_k_single
               ! print*,ER_k_single
                if (stat == 0) then
                   ER_k(i) = ER_k_single
                else
                end if
                read(el_sigma,*,iostat=stat) sigma_single
                !print*,sigma_single
                if (stat == 0) then
                   sigma_k(i) = sigma_single
                else
                end if
             else
             end if
          enddo
       end if
    enddo


!
    ER_k = ER_k/gas_constant
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
       enddo
    enddo
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
    enddo
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
              end if
           else
           end if
        enddo
      enddo
   enddo
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
    enddo
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
      enddo
!
      if (.not. is_reactant) then
         temp_list(nlist-end) = species_list(i)
         end = end + 1
      else
         temp_list(1+front) = species_list(i)
         front = front + 1
      end if
   enddo
!
   do i=1,nlist-1
      if (temp_list(i) == 'Cf') then
         temp = temp_list(nlist)
         temp_list(nlist) = 'Cf'
         temp_list(i)=temp
      else
      end if
   enddo
!
   species_list(:nlist) = temp_list(:nlist)

print*,'species_list=',species_list

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
    enddo
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
    enddo
   end subroutine count_species_type
!**********************************************************************
  subroutine count_species(part,species,reactants,products)
!
    character(10), dimension(:,:) :: part
    character(10) :: element
    real :: numeric
    integer :: i,j,jmax,stat,place,number
    integer :: place_reac, place_prod
    logical :: lhs
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
          if (.not. any(temp .eq. element)) then
          temp(place) = element
          place = place+1
          else
          end if
!
          if ((.not. any(temp_reac .eq. element)) .and. lhs) then
          temp_reac(place_reac) = element
          place_reac = place_reac+1
          else
          end if
!
          if ((.not. any(temp_prod .eq. element)) .and. &
               (lhs .eqv. .false.)) then
          temp_prod(place_prod) = element
          place_prod = place_prod+1
          else
          end if
!
          else
          end if
       enddo
    enddo
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
    enddo
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
    enddo
    string = trim(string)
    i=1
    do while (i <len(string))
       if (string(i:i)/='+'.and. &
           string(i:i)/=spc) then
           j=1
          do while (string(i+j:i+j)/='+'.and.&
                 string(i+j:i+j)/=spc .and. &
                 string(i+j:i+j)/=tab)
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
      target_list(k:,ireaction) = '0.0'
   else
   end if
   print*, target_list(:,ireaction)
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
       enddo
500 if(talk=='verbose') print*,'Done parsing mechanics file'
       close(20)
!
       call remove_save_T_k(target_list)
       call remove_save_powers(target_list)
!
       if(talk=='verbose') then
          open(29, file='mech_outputfile.dat',iostat=stat)
          do i=1,N_surface_reactions
             write(*,writeformat) target_list(:,i),reaction_direction(i)
             write(29,writeformat) target_list(:,i),reaction_direction(i)
          enddo
          close(29)
       else
       end if
510 format (A150)
    else
       write(*,*) 'Could not open mechanics file'
    endif
!
  end subroutine read_mechanics_file
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
      if (N_adsorbed_species > 1) then
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
  subroutine get_reverse_K_k(l,fp)
!
    integer :: l,k,k1,k2
    real, dimension(:), allocatable :: k_p,k_c
    real, dimension(:), allocatable :: denominator, exponent_
    real, dimension(:,:) :: fp
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!    
    allocate(k_p(k1:k2))
    allocate(k_c(k1:k2))
    allocate(denominator(k1:k2))
    allocate(exponent_(k1:k2))
!
    do k=k1,k2
       denominator(k) = heating_k(k,l-1) - (entropy_k(k,l-1)*fp(k,iTp))
       exponent_(k) = denominator(k)/(gas_constant*fp(k,iTp))
       k_p(k) = exp(-exponent_(k))
       k_c(k) = k_p(k) / (((gas_constant)*fp(k,iTp)/interp_pp(k))**(dngas(l-1)))
       K_k(k,l) = (K_k(k,l-1) / k_c(k))
    enddo
!
    deallocate(k_p)
    deallocate(k_c)
    deallocate(denominator)
    deallocate(exponent_)
!
  end subroutine get_reverse_K_k
!**********************************************************************
  subroutine calc_conversion(fp)
!
    real, dimension(:,:) :: fp
    integer :: k
!
    do k=k1_imn(imn),k2_imn(imn)
      conversion(k) = fp(k,imp) / init_mass
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
!!$  subroutine get_mod_surf_area(var,start,end)
!!$!
!!$    real, dimension(:) :: var
!!$    integer :: start, end
!!$!
!!$    var = mod_surf_area(start:end)
!!$!
!!$  end subroutine get_mod_surf_area
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
!
    real, dimension(mpar_loc,mpvar) :: fp
    real, dimension(mx,my,mz,mfarray) :: f
    real :: pre_Cg, pre_Cs, pre_RR_hat
    integer :: i,j,k,k1,k2

    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
!  The heterogeneous kinetics in the mechanism file is always given in SI units. 
!  Here we intorduce correction factors if cgs units are used.
!
    if (unit_system=='cgs') then
      pre_Cg=1e6
      pre_Cs=1e4
      pre_RR_hat=1e-4
    else
      pre_Cg=1.
      pre_Cs=1.  
      pre_RR_hat=1.
    endif
!  
    do k=k1,k2
      do j=1,N_surface_reactions
        RR_hat(k,j)=K_k(k,j)*reaction_enhancement(j)
        do i=1,N_surface_reactants
          if (nu(i,j) > 0) RR_hat(k,j)=RR_hat(k,j)*&
              (pre_Cg*Cg(k)*fp(k,isurf-1+i))**nu(i,j)
        enddo
        if (N_adsorbed_species>1) then
          do i=1,N_adsorbed_species
            if (mu(i,j) > 0) RR_hat(k,j)=RR_hat(k,j)*(pre_Cs*Cs(k,i))**mu(i,j)
          enddo
        endif
        RR_hat(k,j) = RR_hat(k,j)*(fp(k,iTp)**T_k(j))
      enddo
    enddo


!
!!$    print*,'RR_hattwo'
!!$    write(*,'(12E12.4)') RR_hat(k1,:)
!!$    print*, 'Cs' 
!!$    write(*,'(4E12.4)')Cs(1,1:N_adsorbed_species)
!!$    print*, 'xsurf' 
!!$    write(*,'(3E12.4)')fp(1,isurf:isurf_end)
!
!  Make sure RR_hat is given in the right unit system (above it is always
!  in SI units).
!
    RR_hat=pre_RR_hat*RR_hat

  write(*,'(A13,12E10.3)') 'RR_hat',RR_hat(k1,:)
  print*, 'Cg',Cg(k1)
  print*,'fp_surf',fp(k1,isurf:isurf_end)
  print*,'Cs',Cs(k1,:)
!
!  Adapt the reaction rate according to the internal gradients, 
!  after thiele. (8th US combustion Meeting, Paper #070CO-0312)
!  equation 56 ff.
!
    if (lthiele) then
       call calc_effectiveness_factor(fp)
       do j=1,N_surface_reactions
         RR_hat(:,j) = RR_hat(:,j) * effectiveness_factor_reaction(:,j)
       enddo
    endif
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

    ndot=0.    
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
    Rck=0.
!
    do k=k1,k2
       do l=1,N_surface_reactions
          do i=1,N_surface_species
             Rck(k,l)=Rck(k,l)+mol_mass_carbon*RR_hat(k,l)&
                 *(nu_prime(i,l)-nu(i,l))*ac(i)
             ndot(k,i)=ndot(k,i)+RR_hat(k,l)*(nu_prime(i,l)-nu(i,l))*St(k)/ &
                  (fp(k,iap)*fp(k,iap)*4.*pi)
          enddo
       enddo
    enddo
!
! Find molar reaction rate of adsorbed surface species
!
  if (N_adsorbed_species>1) then
    R_j_hat=0.
    do k=k1,k2
       do l=1,N_surface_reactions
          do j=1,N_adsorbed_species-1
             R_j_hat(k,j)=R_j_hat(k,j)+(mu_prime(j,l)-mu(j,l))*RR_hat(k,l)
             Rck(k,l)=Rck(k,l)+mol_mass_carbon*RR_hat(k,l)&
                  *(mu_prime(j,l)-mu(j,l))*aac(j)
          enddo
       enddo
    enddo
 endif
!
! Find mdot_ck
!
  do k=k1,k2
     do l=1,N_surface_reactions
        mdot_ck(k,l)=-St(k)*Rck(k,l)
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
!!$  subroutine get_total_carbon_sites(var)
!!$!
!!$    real :: var
!!$!
!!$    var = total_carbon_sites
!!$!
!!$  end subroutine get_total_carbon_sites
!**********************************************************************
  subroutine calc_effectiveness_factor(fp)
!
!  01-Oct-2014/Jonas: coded
!  taken from solid_reac L. 149 and equations.pdf eq 35ff
!
    real, dimension(:,:) :: fp
!
    real, dimension(:,:), allocatable :: R_i_hat,D_eff
    real, dimension(:), allocatable :: Knudsen, pore_radius
    real, dimension(:), allocatable :: tmp1,tmp2,tmp3
    real, dimension(:), allocatable ::  phi,sum_eta_i_R_i_hat_max
    real, dimension(:), allocatable :: sum_R_i_hat_max
    real, dimension(:), allocatable :: volume
    integer :: i,dep,n,k,k1,k2,l
    real :: r_f
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
    volume(:)=4.*pi*(fp(k1:k2,iap)**3)/3.
!
! Find pore radius
!
    r_f=2.
    pore_radius(:)=2*r_f*porosity(:)*volume(:)/St(k1:k2)
!
!  Find effective diffusion coefficient (based on bulk and Knudsen diffusion)
!
!  JONAS: check which attributes change from particle to particle!!
!
    do k=k1,k2
    do i=1,N_surface_reactants
       tmp3(k)=8*gas_constant*fp(k,iTp)/(pi*molar_mass(jmap(i)))
       Knudsen(k)=2*pore_radius(k)*porosity(k)*sqrt(tmp3(k))/(3*tortuosity)
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
                  (fp(k,imp)*Cg(i)*fp(k,iads-1+i)*D_eff(k,i)))
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
!  Unit: J/(mol)
!
    do k=k1,k2
       if (inuH2O>0) surface_species_enthalpy(k,inuH2O)=-242.18e3-(5.47*fp(k,iTp))
       if (inuO2>0)  surface_species_enthalpy(k,inuO2) = 0.
       if (inuCO2>0) surface_species_enthalpy(k,inuCO2)=-392.52e3-(2.109*fp(k,iTp))
       if (inuH2>0)  surface_species_enthalpy(k,inuH2 )= 0.
       if (inuCO>0)  surface_species_enthalpy(k,inuCO )=-105.95e3-6.143*fp(k,iTp)
       if (inuCH>0)  surface_species_enthalpy(k,inuCH )= 594.13e3
       if (inuHCO>0) surface_species_enthalpy(k,inuHCO)=  45.31e3-(5.94*fp(k,iTp))
       if (inuCH2>0) surface_species_enthalpy(k,inuCH2)= 387.93e3-(5.8*fp(k,iTp))
       if (inuCH4>0) surface_species_enthalpy(k,inuCH4)= -75e3
       if (inuCH3>0) surface_species_enthalpy(k,inuCH3)= 144.65e3-(6.79*fp(k,iTp))
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
! JONAS: units are in j/(mol*K)
!
    do k=k1,k2
    if (inuH2O>0) surface_species_entropy(k,inuH2O)= & 
         189.00+(0.0425*fp(k,iTp))
    if (inuO2>0)  surface_species_entropy(k,inuO2) = &
         222.55+(0.0219*fp(k,iTp))
    if (inuCO2>0) surface_species_entropy(k,inuCO2)= &
         212.19+(0.0556*fp(k,iTp))
    if (inuH2>0)  surface_species_entropy(k,inuH2 )= &
         133.80+(0.0319*fp(k,iTp))
    if (inuCO>0)  surface_species_entropy(k,inuCO )= &
         199.35+(0.0342*fp(k,iTp))
!
!  taken from chemistry  webbook (1bar)
!
    if (inuCH>0)  surface_species_entropy(k,inuCH )=183.00
    if (inuHCO>0) surface_species_entropy(k,inuHCO)=223.114+(0.0491*fp(k,iTp))
    if (inuCH2>0) surface_species_entropy(k,inuCH2)=193.297+(0.0467*fp(k,iTp))
!
!  taken from chemistry webbook (1bar)
!
    if (inuCH4>0) surface_species_entropy(k,inuCH4)= 189.00
    if (inuCH3>0) surface_species_entropy(k,inuCH3)= 190.18+(0.0601*fp(k,iTp))
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
!  Unit: J/(mol*K)
!
    adsloop: do k=k1,k2
!
    if (imuadsO>0)    then
       adsorbed_species_entropy(k,imuadsO) = &
    (164.19+(0.0218*fp(k,iTp)))*0.72 - (3.3*gas_constant)
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
         ((0.0319*fp(k,iTp)) + 186.88) * 0.7 - (3.3*gas_constant)
    else
    end if
    if (imuadsH>0)    then 
       adsorbed_species_entropy(k,imuadsH) = &
           (117.49+(0.0217*fp(k,iTp)))*0.54 - (3.3*gas_constant)
    else
    end if
    if (imuadsCO>0)   then
       adsorbed_species_entropy(k,imuadsCO) = &
           (199.35+(0.0342*fp(k,iTp))) * &
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
!  code Units: J/mol
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    if (imuadsO>0)    adsorbed_species_enthalpy(k1:k2,imuadsO) = &
         -148.14e3 + (0.0024e3*(fp(k1:k2,iTp)-273.15))
    if (imuadsO2>0)   adsorbed_species_enthalpy(k1:k2,imuadsO2)= &
         2 *  (-148.14e3 + (0.0024e3*(fp(k1:k2,iTp)-273.15)))
    if (imuadsOH>0)   adsorbed_species_enthalpy(k1:k2,imuadsOH)= &
         -148e3
    if (imuadsH>0)    adsorbed_species_enthalpy(k1:k2,imuadsH) = &
         -19.5e3
    if (imuadsCO>0)   adsorbed_species_enthalpy(k1:k2,imuadsCO)= &
         -199.94e3 - (0.0167e3*(fp(k1:k2,iTp)-273.15))
    if (imufree>0)    adsorbed_species_enthalpy(k1:k2,imufree) = 0.
!
  end subroutine calc_ads_enthalpy
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
    if(lpchem_debug) then
       call print_debug_info()
    else
    endif
!
  end subroutine calc_St_init
!*********************************************************************
  subroutine calc_pchemistry_pencils(f,fp,p,ineargrid)
!
    real, dimension(mpar_loc,mpvar) :: fp
    real, dimension(mx,my,mz,mfarray) :: f
    integer, dimension(mpar_loc,3) :: ineargrid
    type (pencil_case) :: p
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
    call calc_K_k(f,fp,p,ineargrid)
    call calc_Cg(p,ineargrid)
    if (N_adsorbed_species>1) call calc_Cs(fp)
    call calc_A_p(fp)
    call calc_mass_trans_coeff(f,fp,p,ineargrid)
    call calc_x_surf(f,fp,ineargrid)
    call calc_RR_hat(f,fp)
!
    if (lreactive_heating) then
       call calc_q_reac()
       call calc_Nusselt()
    else
       q_reac=0.0
    endif
!
    call calc_ndot_mdot_R_j_hat(fp)
    call calc_R_c_hat()
!
  end subroutine calc_pchemistry_pencils
!********************************************************************
  subroutine allocate_variable_pencils()
!
!  allocate variables used in the pencil variable calculation
!
    integer :: k1,k2,stat
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    allocate(rho_p(k1:k2) ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for rho_p')
    allocate(q_reac(k1:k2) ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for q_reac')
    allocate(porosity(k1:k2) ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for porosity')
    allocate(St(k1:k2) ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for St')
    allocate(mod_surf_area(k1:k2) ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for mod_surf_area')
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
    if (stat>0) call fatal_error('allocate_variable_pencils',&
         'Could not allocate memory for adsorbed_species_entropy')
    allocate(adsorbed_species_enthalpy(k1:k2,N_adsorbed_species) &
         ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
         'Could not allocate memory for adsorbed_species_enthalpy')
!
    allocate(thiele(k1:k2,N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for thiele')
!
    allocate(effectiveness_factor(k1:k2),STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for thiele')
    allocate(effectiveness_factor_species(k1:k2,N_surface_reactants),STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for effectiveness_factor_species')
    allocate(effectiveness_factor_reaction(k1:k2,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('allocate_variable_pencils',&
        'Could not allocate memory for effectiveness_factor_reaction')
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
      allocate(K_k(k1:k2,N_surface_reactions)   ,STAT=stat)
      if (stat>0) call fatal_error('allocate_variable_pencils',&
           'Could not allocate memory for K_k')
      allocate(mass_trans_coeff_species(k1:k2,N_species), STAT=stat)
      if (stat>0) call fatal_error('allocate_variable_pencils',&
           'Could not allocate memory for mass_trans_coeff_species')
    allocate(mass_trans_coeff_reactants(k1:k2,N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for mass_trans_coeff_reactants')
      allocate(Cg(k1:k2), STAT=stat)
      if (stat>0) call fatal_error('allocate_variable_pencils',&
           'Could not allocate memory for Cg')
      allocate(A_p(k1:k2), STAT=stat)
      if (stat>0) call fatal_error('allocate_variable_pencils',&
           'Could not allocate memory for A_p')
      allocate(x_surf(k1:k2,N_surface_species), STAT=stat)
      if (stat>0) call fatal_error('allocate_variable_pencils',&
           'Could not allocate memory for x_surf')
     allocate(Nu_p(k1:k2), STAT=stat)
      if (stat>0) call fatal_error('allocate_variable_pencils',&
           'Could not allocate memory for Nu_p')
!
  end subroutine allocate_variable_pencils
!***************************************************
  subroutine cleanup_chemistry_pencils()
!
!  06-oct-14/jonas:coded
!
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
    deallocate(porosity)
    deallocate(St)
    deallocate(K_k)
    deallocate(mass_trans_coeff_species)
    deallocate(mass_trans_coeff_reactants)
    deallocate(Cg)
    deallocate(A_p)   
    deallocate(x_surf)
    deallocate(q_reac)
    deallocate(Nu_p)
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
      rho_p_init(:) = fp(:,imp) / (4.*pi/3. * fp(:,iap)**3)
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
         rho_p(k) = fp(k,imp) / (4./3. *pi * fp(k,iap)**3)
         porosity(k)=1-rho_p(k)/true_density_carbon
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
    enddo

    do k=1,N_surface_reactions
       write(*,writeformat) 'nu_power=',k,nu_power(:,k)
    enddo
!
    do k=1,N_surface_reactions
       write(*,writeformat) 'nu_prime=',k,nu_prime(:,k)
    enddo
!
    write(writeformat(13:14),'(I2)') N_adsorbed_species
!
    print*,'Adsorbed_species_names'
    print*, adsorbed_species_names
    do k=1,N_surface_reactions
       write(*,writeformat) 'mu=',k,mu(:,k)
    enddo 
!
    do k=1,N_surface_reactions
       write(*,writeformat) 'mu_power=',k,mu_power(:,k)
    enddo
 !
    do k=1,N_surface_reactions
       write(*,writeformat) 'mu_prime=',k,mu_prime(:,k)
    enddo
!
    do k=1,N_surface_reactions
       write(*,'(A12,I4,2E12.5)') 'ER_k, B_k=',k,B_k(k),ER_k(k)/ &
               (gas_constant)
    enddo
    do k=1,N_surface_reactions
       write(*,'(A12,I4,E12.5)') 'Dngas',k,dngas(k)
    enddo
    do k=1,N_surface_reactions
       write(*,'(A12,I4,E12.5)') 'sigma',k,sigma_k(k)
    enddo
    do k=1,N_surface_reactions
       write(*,'(A12,I4,I4)') 'dep',k,dependent_reactant(k)
    enddo
!        
        write(*,'(A20," ",10I4)') 'jmap=', jmap
        write(*,'(A20," ",10F4.0)') 'ac=',ac
        write(*,'(A20," ",10F4.0)') 'site_occupancy=',aac
!
  end subroutine print_debug_info
!*******************************************************************
  subroutine calc_K_k(f,fp,p,ineargrid)
!
!  27-10-2014/jonas:coded
!  calculation of K_k from the Arrhenius factors
!
    real, dimension(:,:,:,:) :: f
    real, dimension(:,:) :: fp
    type (pencil_case) :: p
    integer, dimension(:,:) :: ineargrid
    integer :: k,k1,k2,i,j,l
    integer :: N_iter = 20
    real :: int_k, gg_old,gg,ff,energy,dE,delta_E
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
!  Calculation of the 'plain' K_k (reaction rate)
!
    particles: do k=k1,k2
        do l=1,N_surface_reactions
           if (sigma_k(l)==0) then
              K_k(k,l)=B_k(l)*exp(-ER_k(l)/fp(k,iTp))
           else
              delta_E=sigma_k(l)*6
              dE=2*delta_E/(N_iter-1)
              energy=ER_k(l)*gas_constant-delta_E
              if (energy<0) then
                print*,'k,delta_E,sigma_k(k),ER_k(k)=',&
                    l,delta_E,sigma_k(l),ER_k(l)
                call fatal_error('get_reaction_rates_solid',&
                    'delta_E is too large!')
              endif
              ff=exp(-0.5*((energy-ER_K(l)*gas_constant)/sigma_k(l))**2)/&
                   (sigma_k(l)*sqrt(2*pi))
              gg_old=ff*B_k(l)*exp(-energy/(gas_constant*fp(k,iTp)))
              int_k=0
              eff_kk: do j=2,N_iter
                 energy=(j-1)*2*delta_E/(N_iter-1)+ER_k(l)*gas_constant-delta_E
                 ff=exp(-0.5*((energy-ER_K(l)*gas_constant)/sigma_k(l))**2)/&
                      (sigma_k(l)*sqrt(2*pi))
                 gg=ff*B_k(l)*exp(-energy/(gas_constant*fp(k,iTp)))
                 int_k=int_k+dE*(gg_old+0.5*(gg-gg_old))
                 gg_old=gg
              enddo eff_kk
              K_k(k,l)=int_k
           endif
              if (reaction_direction(l) == 'rev') then
                 call get_reverse_K_k(l,fp)
              else
              end if
        enddo
    enddo particles
!
!  Application of the startup-quenching to facilitate convergence on
!  startup
!
    if (t < startup_time) then
      startup_quench=(tanh(6*(t-startup_time/2)/(startup_time))+1.0)/2.
      if (iter<1) startup_quench=0.
      K_k=K_k*startup_quench
    endif
!    
  end subroutine calc_K_k
!*******************************************************************
  subroutine calc_x_surf(f,fp,ineargrid)
!
!  Routine to calculate the gas composition in the particle
!  near field
!
    real, dimension(:,:,:,:) :: f
    real, dimension(:,:) :: fp
    integer :: k,k1,k2,i,j,l
    integer, dimension(:,:) :: ineargrid
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
!  allocation and calculation of factors used in determining the
!  gas surface composition
!
      do k=k1,k2
         do i=1,N_surface_species
            x_surf(k,i) = fp(k,isurf-1+i)
         enddo
      enddo
!
  end subroutine calc_x_surf
!*******************************************************************
  subroutine calc_Cg(p,ineargrid)
!
!  Gas concentration in the cell of the particle
!
   integer :: k,k1,k2
   integer :: ix0
   type (pencil_case) :: p
   integer, dimension(:,:) :: ineargrid
!
   k1 = k1_imn(imn)
   k2 = k2_imn(imn)
!
   do k=k1,k2
         Cg(k) = interp_pp(k)/(R_cgs*interp_TT(k))
   enddo
!
  end subroutine calc_Cg
!************************************************************************
  subroutine calc_mass_trans_coeff(f,fp,p,ineargrid)
!
!  diffusion coefficients of gas species at nearest grid point
!
  real, dimension(:,:,:,:) :: f
  real, dimension(:,:) :: fp
  type (pencil_case) :: p
  integer, dimension(:,:) :: ineargrid
  integer :: k,k1,k2,i
  integer :: ix0,iy0,iz0
  integer::spec_glob,spec_chem
  real, dimension(:,:), allocatable :: diff_coeff_species
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      allocate(diff_coeff_species(k1:k2,N_species))
!
      do k=k1,k2
         ix0 = ineargrid(k,1)
         iy0 = ineargrid(k,2)
         iz0 = ineargrid(k,3)
         do i=1,N_surface_species
            diff_coeff_species(k,i) = p%Diff_penc_add(ix0-nghost,jmap(i))
         enddo
      enddo
!
      do k=k1,k2
         do i=1,N_surface_species
          mass_trans_coeff_species(k,i)=Cg(k)*diff_coeff_species(k,i)/ &
               fp(k,iap)
         enddo
      enddo
!
      deallocate(diff_coeff_species)
!
  end subroutine calc_mass_trans_coeff
!*************************************************************************
  subroutine calc_A_p(fp)
!
!  particle surface
!
    real, dimension(:,:) :: fp
    integer :: k1,k2
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    A_p(k1:k2) = 4 * pi * fp(k1:k2,iap) * fp(k1:k2,iap)    
!
  end subroutine calc_A_p
!*************************************************************************
  subroutine calc_Cs(fp)
!
!  calculate site concentration
!
    real, dimension(:,:) :: fp
    integer :: k,k1,k2,i
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    do k=k1,k2
       do i=1,N_adsorbed_species-1
          Cs(k,i) = fp(k,iads+i-1)
       enddo
       Cs(k,N_adsorbed_species) = 1 - sum(Cs(k,1:N_adsorbed_species-1))
    enddo
!
    Cs = Cs * total_carbon_sites
!
  end subroutine calc_Cs
!*************************************************************************
  subroutine calc_q_reac()
!
!  calculate the reactive heating of the particle
!
        integer :: k,k1,k2,i
        real, dimension(:,:), allocatable :: qk_reac
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)
!
    allocate(qk_reac(k1:k2,N_surface_reactions))
    do k=k1,k2
       qk_reac(k,:) = RR_hat(k,:)*St(k)*(-heating_k(k,:))
       q_reac(k) = sum(qk_reac(k,:))
    enddo
!
!  JONAS: UNITS!!!!! qreac is in Joule, which is 10^7ergs
!
    if (allocated(q_reac)) deallocate(q_reac)
!
  end subroutine calc_q_reac
!*************************************************************************
  subroutine get_q_reac(var)
!
!  transport routine to particles_temperature module
!
!
  real, dimension(:) :: var
!
  var = q_reac
!
  end subroutine get_q_reac
!************************************************************************
    subroutine calc_Nusselt()
!
!
!
      integer :: k,k1,k2
      real :: Pr_g, Sc_g,rep
!
      k1=k1_imn(imn)
      k2=k2_imn(imn)
!
!  JONAS: access to some chemistry variables as well to particles_dust for rep(k) 
!  for Pr and Sc is needed, these are only placeholders!!!
!
         Pr_g = 0.7 
         Sc_g = 0.65
         rep = 0.0002
!
      do k=k1,k2
!  case from the paper "Dimensionless heat-mass transfer coefficients 
!  for forced convection around a sphere: a general low 
!  Reynolds number correlation
         if (rep<3.5) then
            Nu_p(k) = 0.5+0.5*(1+2*rep*Pr_g)**(1./3.)
!  First case in "The effect of particle packing on the Nusselt 
!  number in a cluster of particles
         else if (rep>=3.5 .and. rep<=7.6e4 .and. &
              Pr_g>=0.7 .and. Pr_g<=380.) then
            Nu_p(k)=2.0+(0.4*rep**0.5+0.06*rep**0.67)*&
                 Pr_g**0.4
         else
            Nu_p=2.0
         endif
      enddo
!
    end subroutine calc_Nusselt
!***********************************************************************
    subroutine get_Nusselt(var)
!
!
!
      real, dimension(:) :: var
!
      var = Nu_p
!
    end subroutine get_Nusselt
!**********************************************************************
  subroutine remove_save_T_k(target_list)
!
    character(10), dimension(:,:) :: target_list
    character(10) :: el_T_k
    real :: T_k_single
    integer :: i,k,stat
!
    T_k=0.0
    do i=1,size(target_list,2)
       do k=1,size(target_list,1)
          if (target_list(size(target_list,1),i) == 'rev') then
          else
             el_T_k = target_list(k,i)
             if (el_T_k(1:2) == 'T^') then
                read(el_T_k(3:),*,iostat=stat) T_k_single
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
!
  end subroutine remove_save_T_k
!********************************************************************
  end module Particles_chemistry
