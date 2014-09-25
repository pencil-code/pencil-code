! $Id: particles_adsorbed.f90 21950 2014-07-08 08:53:00Z jonas.kruger $
!
!  This module takes care of everything related to reactive particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
! NADSSPEC CONTRIBUTION 5
!
! CPARAM logical, parameter :: lparticles_adsorbed=.true.
!
!! PENCILS PROVIDED adsp
!  PENCILS PROVIDED theta(nadsspec)
!  PENCILS PROVIDED x_surf(nsurfreacspec)
!
!***************************************************************
module Particles_adsorbed
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_radius
  use Particles_chemistry
  use Particles_temperature
!
  implicit none
!
  include 'particles_adsorbed.h'
!
  character (len=labellen), dimension (ninit) :: init_adsorbed='nothing'
  character(10), dimension(15)  :: species_name,adsorbed_species_names
  real, dimension(10) :: init_surf_ads_frac,init_surf_gas_frac
  real, dimension(:,:), allocatable :: R_j_hat
  real :: diffusivity=0.0
  real :: init_thCO=0.0
  real :: total_carbon_sites=1.08e-7
  integer :: iads=0, iads_end=0
  real :: gas_constant=8314.0 ![J/kmol/K]
!
!*****************************************************************
! Particle number independent variables
!*****************************************************************
!
  real, dimension(:,:), allocatable :: mu, mu_prime
  real, dimension(:), allocatable :: site_occupancy
  real, dimension(:), allocatable :: aac
  integer :: imufree, imuadsO, imuadsO2, imuadsOH, imuadsH, imuadsCO
  integer :: N_surface_species, N_surface_reactions, nr,ns,np
  integer :: N_adsorbed_species
!
!*********************************************************************!
!               Particle dependent variables below here               !
!*********************************************************************!
!  
  real, dimension(:,:), allocatable :: adsorbed_species_enthalpy
  real, dimension(:,:), allocatable :: adsorbed_species_entropy
  real, dimension(:,:), allocatable :: Cs
!
  namelist /particles_ads_init_pars/ &
      init_adsorbed, &
      init_surf_ads_frac, &
      init_surf_gas_frac, diffusivity, &
      total_carbon_sites
!
  namelist /particles_ads_run_pars/ &
      diffusivity
!
  contains
!***********************************************************************
  subroutine register_particles_ads()
!
    character(10), dimension(40) :: species
!  this is a wrapper function for registering particle number
!  in- and independent variables
!
    call get_pchem_info(species,'N_adsorbed_species',N_adsorbed_species,'quiet')
    call register_indep_ads()
    call register_dep_ads()
 !   
  end subroutine register_particles_ads 
!***********************************************************************
  subroutine register_indep_ads()
!
!  Set up indices for access to the fp and dfp arrays
!
!  29-aug-14/jonas: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      integer :: i,N_surface_species,stat,j
      character(10), dimension(40) :: species
      character(10), dimension(40) :: solid_species, adsorbed_species_names
!
      if (lroot) call svn_id( &
           "$Id: particles_adsorbed.f90 21950 2014-07-08 08:53:00Z jonas.kruger $")
!
!  check if enough storage was reserved in fp and dfp to store
!  the adsorbed species surface and gas phase surface concentrations
!
      call get_pchem_info(species,'N_adsorbed_species',N_adsorbed_species,'quiet')
!
      if (nadsspec/=(N_adsorbed_species-1) .and. &
           N_adsorbed_species>0) then
         print*,'N_adsorbed_species: ',N_adsorbed_species
         call fatal_error('register_particles_ads', &
              'wrong size of storage for adsorbed species allocated.')
         else
      endif
!
      call get_species_list('adsorbed_species_names',adsorbed_species_names)
      call sort_compounds(reactants,adsorbed_species_names,N_adsorbed_species,nr)
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
!  Indices of adsorbed species mole fraction
!
      if (N_adsorbed_species>1) then
         j=1
         iads = npvar+1
         i=1
!  JONAS: commented for now, wait for own pvarnamespace
!         do while (i<=(N_adsorbed_species-1))
!            if (adsorbed_species_names(i)/='Cf') then
!            pvarname(iads+j-1) = adsorbed_species_names(i)
!            i = i+1
!            j = j+1
!            else
!            i = i+1
!            end if 
!         enddo
!
!  Increase of npvar according to N_adsorbed_species
!  The -2 is there to account for Cf being in adsorbed_species_names
!  but not in the calculation
!
         npvar=npvar+N_adsorbed_species-1
         iads_end=iads+N_adsorbed_species-2
      else
         call fatal_error('register_particles_ads', &
              'N_adsorbed_species must be > 1')
      endif
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_ads','npvar > mpvar')
      endif
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
    call create_stoc(part,adsorbed_species_names,mu,.true.,N_adsorbed_species)
    call create_stoc(part,adsorbed_species_names,mu_prime,.false.,&
        N_adsorbed_species)
    call create_occupancy(adsorbed_species_names,site_occupancy)
!
    end subroutine register_indep_ads
!***********************************************************************
  subroutine register_dep_ads()
!
      integer :: stat
      character(10), dimension(40) :: trash   
!
      allocate(R_j_hat(mpar_loc,N_adsorbed_species-1)   ,STAT=stat)
      if (stat>0) call fatal_error('register_particles_ads',&
           'Could not allocate memory for R_j_hat')
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
!
  end subroutine register_dep_ads
!***********************************************************************

    subroutine initialize_particles_ads(fp,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-aug-14/jonas coded
!
      real, dimension (mpar_loc,mpvar) :: fp
      logical :: lstarting
!
      if(lstarting) then
         call get_initial_values(fp)
      else
      endif
!
    end subroutine initialize_particles_ads
!***********************************************************************
    subroutine init_particles_ads(f,fp)
!
!  Initial particle surface fractions
!
!  01-sep-14/jonas: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer :: j,i
!
      intent (out) :: f, fp
!
      call keep_compiler_quiet(f)
!
      fp(:,iads:iads_end)=0.
      do j=1,ninit
!
!  Writing the initial adsorbed species fractions
!
        select case (init_adsorbed(j))
        case ('nothing')
          if (lroot .and. j==1) print*, 'init_particles_ads,adsorbed: nothing'
        case ('constant')
          if (lroot) print*, 'init_particles_ads: Initial Adsorbed Fractions'
          do i=1,mpar_loc
             fp(i,iads:iads_end)=fp(i,iads:iads_end)+ &
                  init_surf_ads_frac(1:(iads_end-iads+1))
          enddo
        case default
          if (lroot) &
              print*, 'init_particles_ads: No such such value for init_adsorbed: ', &
              trim(init_adsorbed(j))
          call fatal_error('init_particles_ads','')
        endselect
      enddo
!
    endsubroutine init_particles_ads
!***********************************************************************
subroutine pencil_criteria_par_ads()
!
!  All pencils that the Particles_adsorbed
!  module depends on are specified here.
!
!  01-sep-14/jonas: coded
!
    endsubroutine pencil_criteria_par_ads
!***********************************************************************
    subroutine dpads_dt(f,df,fp,dfp,ineargrid)
!
!  Evolution of particle surface fractions.
!
!  01-sep-14/jonas: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpads_dt
!***********************************************************************
    subroutine dpads_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle surface fractions
!
!  03-sep-14/jonas: coded
!
      real, dimension (mx,my,my,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      real, dimension (mpar_loc) :: mod_surf_area,R_c_hat
      type (pencil_case) :: p
      real, dimension(mpar_loc,3) :: ineargrid
      integer :: i, n_ads
!
      intent (inout) :: dfp
      intent (in) :: fp
!
!  JONAS incosistency, should i implement another get()?
!   or read it from memory
!
!      call _R_j_hat(R_j_hat)
      call get_mod_surf_area(mod_surf_area)
      call get_R_c_hat(R_c_hat)
!
      n_ads = iads_end - iads +1
!
      do i=1,mpar_loc
         dfp(i,iads:iads_end)=R_j_hat(i,1:n_ads)/ &
              total_carbon_sites + mod_surf_area(i)* &
              R_c_hat(i)*fp(i,iads:iads_end)
      enddo
!
    endsubroutine dpads_dt_pencil
!***********************************************************************
    subroutine read_particles_ads_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_ads_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_ads_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_ads_init_pars
!***********************************************************************
    subroutine write_particles_ads_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_ads_init_pars)
!
    endsubroutine write_particles_ads_init_pars
!***********************************************************************
    subroutine read_particles_ads_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_ads_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_ads_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_ads_run_pars
!***********************************************************************
    subroutine write_particles_ads_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_ads_run_pars)
!
    endsubroutine write_particles_ads_run_pars
!***********************************************************************
    subroutine rprint_particles_ads(lreset,lwrite)
!
!  Read and register print parameters relevant for
! vparticles coverage fraction.
!
!  29-aug-14/jonas: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'iox=', iox
!
      call keep_compiler_quiet(lreset)
!
    end subroutine rprint_particles_ads
!***********************************************************************
    subroutine particles_ads_prepencil_calc(f)
!
!  28-aug-14/jonas+nils: coded
!
      real,dimension(mx,my,mz,mfarray),intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    end subroutine particles_ads_prepencil_calc
!***********************************************************************
   subroutine calc_R_j_hat()
!
      integer :: j,k
!
!  Calculation of R_j_hat according to eq.50 in 8th US combustion 
!  meeting, coal and  biomass combustion and gasification.
!
!  JONAS: RR_hat still needed, get() from chemistry ?
!      
      R_j_hat = 0.0
      do k=1,N_surface_reactions
         do j=1,N_adsorbed_species-1
            R_j_hat(:,j)=R_j_hat(:,j)+(mu_prime(j,k)-mu(j,k))*RR_hat(:,k)
         enddo
      enddo
!
    end subroutine calc_R_j_hat
!***********************************************************************
  subroutine calc_pchem_factors()
!
    call calc_conversion(fp)
    call calc_St(fp)
    call calc_mod_surf_area(fp)
    call calc_R_c_hat()
!
       call calc_R_j_hat()
!
  end subroutine calc_pchem_factors
!***********************************************************************
  subroutine calc_ads_entropy(fp)
!
    real, dimension(mpar_loc,mpvar) :: fp
!
!  Unit: J/kmolK
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
  end subroutine calc_ads_entropy
!*********************************************************************
  subroutine calc_ads_enthalpy(fp)
!
    real, dimension(mpar_loc,mpvar) :: fp
    real, dimension(mpar_loc) :: enth_ads_O, enth_ads_O2
    real, dimension(mpar_loc) :: enth_ads_OH,enth_ads_H
    real, dimension(mpar_loc) :: enth_ads_CO
!
!  JONAS: values are from nist and solid_phase.f90 of the stanford
!  code Units: J/kmol
!
    

    enth_ads_O = -148.14e6 + (0.0024e6*(fp(:,iTp)-273.15))
    enth_ads_CO = -199.94e6 - (0.0167e6*(fp(:,iTp)-273.15))
    enth_ads_H = -19.5e6
    enth_ads_OH = -148e6
    enth_ads_O2(:) =2 * enth_ads_O(:)
!
    if (imuadsO>0)    adsorbed_species_enthalpy(:,imuadsO) = enth_ads_O
    if (imuadsO2>0)   adsorbed_species_enthalpy(:,imuadsO2)= enth_ads_O2
    if (imuadsOH>0)   adsorbed_species_enthalpy(:,imuadsOH)= enth_ads_OH
    if (imuadsH>0)    adsorbed_species_enthalpy(:,imuadsH) = enth_ads_H
    if (imuadsCO>0)   adsorbed_species_enthalpy(:,imuadsCO)= enth_ads_CO
    if (imufree>0)    adsorbed_species_enthalpy(:,imufree) = 0.
!
  end subroutine calc_ads_enthalpy
!*********************************************************************
  end module Particles_adsorbed
