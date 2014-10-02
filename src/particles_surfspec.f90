! $Id: particles_surfspec.f90 21950 2014-07-08 08:53:00Z jonas.kruger $
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
! NADSSPEC CONTRIBUTION 
!
! CPARAM logical, parameter :: lparticles_surfspec=.true.
!
!***************************************************************
module Particles_surfspec
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_map
  use Particles_mpicomm
  use Particles_sub
  use Particles_chemistry
  use Particles_adsorbed
!
!
  implicit none
!
  include 'particles_surfspec.h'
!
!*********************************************************************!
!               Particle independent variables below here             !
!*********************************************************************!
!
  character (len=labellen), dimension (ninit) :: init_surf='nothing'
  character(10), dimension(40) :: reactants
  real, dimension(10) :: init_surf_gas_frac
  real :: diffusivity=0.0
  real :: surfplaceholder=0.0
  real, dimension(:), allocatable :: ac,dngas
  real, dimension(:), allocatable :: mass_trans_coeff_reactants
  real, dimension(:), allocatable :: diff_coeff_reactants
  real, dimension(:), allocatable :: uscale,fscale,constr
  integer :: inuH2,inuCO2,inuH2O,inuCO,inuCH4,inuO2
  integer :: inuCH,inuHCO,inuCH2,inuCH3
  integer :: nr,ns,np
  integer :: N_max_elements
  integer, dimension(:), allocatable :: j_of_inu
  integer, dimension(:), allocatable :: dependent_reactant
  integer :: iads,iads_end
!
!  JONAS: implement j_of_inu to communicate with
!  gas phase
!
  integer :: jH2=1, jO2=2, jCO2=3, jCO=4, jCH4=5, jN2=6, jH2O=7
  integer :: jOH=8, jAR=9, jO=10, jH=11, jCH=12, jCH2=13, jHCO=14, jCH3=15
!
!*********************************************************************!
!               Particle dependent variables below here               !
!*********************************************************************!
!
  real, dimension(:,:), allocatable :: thiele
  real, dimension(:,:), allocatable :: ndot,K_k
  real, dimension(:,:), allocatable :: X_infty_reactants
  real, dimension(:), allocatable :: initial_density
  real, dimension(:), allocatable :: Qh,Qc,Qreac,Qrad

!
  namelist /particles_surf_init_pars/ &
       init_surf, &
       init_surf_gas_frac
!
  namelist /particles_surf_run_pars/ &
       surfplaceholder
!
  contains
!***********************************************************************
  subroutine register_particles_surfspec()
!
    character(10), dimension(40) :: species
!
!  This is a wrapper routine for particle dependent and particle
!  independent variables
!  JONAS: Back to standalone via mpar_loc=1?
!
    N_surface_reactions = count_reactions('mechanics.in')
!
    call get_pchem_info(species,'N_surface_species',N_surface_species,'quiet')
    call get_pchem_info(species,'N_surface_reactants',N_surface_reactants,'quiet')
    call get_pchem_info(species,'N_species',ns,'quiet')
    call get_pchem_info(species,'N_reactants',nr,'quiet')
    call get_pchem_info(species,'N_adsorbed_species',N_adsorbed_species,'quiet')
    call register_indep_psurfspec()
    call register_dep_psurfspec()
!!$!
  end subroutine register_particles_surfspec
!************************************************************************
  subroutine register_indep_psurfspec()
!
      integer :: i,k,stat
      character(10), dimension(40) :: all_species
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
!  Increase of npvar according to N_surface_species, which is
!  the concentration of gas phase species at the particle surface
!
      if (N_surface_species>1) then
         isurf = npvar+1
!         do i=1,N_surface_species
! JONAS: where do we save this
! JONAS: commented for now
!            pvarname(isurf+i-1)=solid_species(i)
!         enddo
         npvar=npvar+N_surface_species-1
         isurf_end=isurf+N_surface_species-1
      else
         call fatal_error('register_particles_', &
              'N_surface_species must be > 1')
      endif
!
!  Check that the fp and dfp arrays are big enough
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_ads','npvar > mpvar')
      endif
!
! Allocate memory for a number of arrays
!
    allocate(dngas(N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
         'Could not allocate memory for dngas')
    allocate(dependent_reactant(N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for dependent_reactant')
    allocate(nu(N_surface_species,N_surface_reactions),STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for nu')
    allocate(nu_prime(N_surface_species,N_surface_reactions)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for nu_prime')
    allocate(ac(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for ac')
    allocate(j_of_inu(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for j_of_inu')
    allocate(solid_species(N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for solid_species')
    allocate(mass_trans_coeff_reactants(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for mass_trans_coeff_reactants')
    allocate(diff_coeff_reactants(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for diff_coeff_reactants')
    allocate(uscale(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for uscale')
    allocate(fscale(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for fscale')
    allocate(constr(N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for constr')
!
! Define the Stoichiometric matrixes. Here i=1 correspond to H2O, i=2 is CO2,
! i=3 is H2 and finally i=4 is O2
!
    call get_species_list('all',all_species)
    call create_ad_sol_lists(all_species,solid_species,'sol',ns)
    call get_reactants(reactants)
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
    if(inuH2O > 0)     j_of_inu(inuH2O)=jH2O
    if(inuCO2 > 0)     j_of_inu(inuCO2)=jCO2
    if(inuH2 > 0)      j_of_inu(inuH2) =jH2
    if(inuO2 > 0)      j_of_inu(inuO2) =jO2
    if(inuCO > 0)      j_of_inu(inuCO) =jCO
    if(inuCH > 0)      j_of_inu(inuCH) =jCH
    if(inuHCO > 0)     j_of_inu(inuHCO)=jHCO
    if(inuCH2 > 0)     j_of_inu(inuCH2)=jCH2
    if(inuCH3 > 0)     j_of_inu(inuCH3)=jCH3
!
! Set number of carbon atoms for each surface species
!
    call get_ac(ac,solid_species,N_surface_species)
!
!  JONAS: VERY roundabout way of doing things
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
!
      call get_part(part)
!
!  Set the stoichiometric matrixes
!
    call create_stoc(part,solid_species,nu,.true.,N_surface_species)
    call create_stoc(part,solid_species,nu_prime,.false.,N_surface_species)
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
    end subroutine register_indep_psurfspec
!***********************************************************************
    subroutine register_dep_psurfspec()
!
      integer :: stat
!
    allocate(x_surface(mpar_loc,N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_psurfchem',&
        'Could not allocate memory for x_surface')
    allocate(ndot(mpar_loc,N_surface_species)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_psurfchem',&
        'Could not allocate memory for ndot')
    allocate(ndot_total(mpar_loc)   ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_psurfchem',&
        'Could not allocate memory for ndot_total')
    allocate(x_infty_reactants(mpar_loc,N_surface_reactants) &
         ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_psurfchem',&
        'Could not allocate memory for x_infty_reactants')
    allocate(surface_species_enthalpy(mpar_loc,N_surface_species) &
         ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_psurfchem',&
        'Could not allocate memory for surface_species_enthalpy')
    allocate(surface_species_entropy(mpar_loc,N_surface_species) &
         ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_psurfchem',&
        'Could not allocate memory for surface_species_entropy')
    allocate(thiele(mpar_loc,N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for thiele')
    allocate(Rck(mpar_loc,N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for Rck')
    allocate(Rck_max(mpar_loc,N_surface_reactants)   ,STAT=stat)
    if (stat>0) call fatal_error('register_indep_psurfchem',&
        'Could not allocate memory for Rck_max')

!
    end subroutine register_dep_psurfspec
!***********************************************************************
    subroutine read_particles_surf_init_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_surf_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_surf_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_surf_init_pars
!***********************************************************************
    subroutine write_particles_surf_init_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_surf_init_pars)
!
    endsubroutine write_particles_surf_init_pars
!***********************************************************************
    subroutine read_particles_surf_run_pars(unit,iostat)
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=particles_surf_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=particles_surf_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_particles_surf_run_pars
!***********************************************************************
    subroutine write_particles_surf_run_pars(unit)
!
      integer, intent (in) :: unit
!
      write(unit,NML=particles_surf_run_pars)
!
    endsubroutine write_particles_surf_run_pars
!***********************************************************************
    subroutine rprint_particles_surf(lreset,lwrite)
!
!  Read and register print parameters relevant for 
!  solid species volume fraction in the near field.
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
    end subroutine rprint_particles_surf
!***********************************************************************
    subroutine init_particles_surf(f,fp)
!
!  Initial particle surface fractions
!
!  01-sep-14/jonas: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      real :: sum_surf_spec
      integer :: j,i
!
      intent (in) :: f 
      intent (out) :: fp
!
      call keep_compiler_quiet(f)
!
      fp(:,isurf:isurf_end)=0.
      do j=1,ninit
!
!  Writing the initial surface species fractions
!
       init: select case (init_surf(j))
        case ('nothing') init
          if (lroot .and. j==1) print*, 'init_particles_surf,gas phase: nothing'
        case ('constant') init
          if (lroot) print*, 'init_particles_surf: Initial Surface Fractions'
!
!  This ensures that we don't have unphysical values as init
!
          sum_surf_spec = sum(init_surf_gas_frac(1:N_surface_species))
          if (sum_surf_spec > 1) then
             print*, 'Sum of all surface fractions >1, normalizing...'
             init_surf_gas_frac(1:N_surface_species) = &
                  init_surf_gas_frac(1:N_surface_species) / sum_surf_spec
          else
          endif
!
          do i=1,mpar_loc
             fp(i,isurf:isurf_end)=fp(i,isurf:isurf_end) + &
                 init_surf_gas_frac(1:N_surface_species)
          enddo
          case default init
          if (lroot) &
              print*, 'init_particles_ads: No such such value for init_surf: ', &
              trim(init_surf(j))
          call fatal_error('init_particles_surf','')
        endselect init
      enddo
!
    endsubroutine init_particles_surf
!**********************************************************************
  subroutine calc_pchem_factors(f,fp)
!
    real, dimension(mpar_loc,mpvar) :: fp
    real, dimension(mx,my,mz,mfarray) :: f
    integer :: stat,k1,k2
!
!  Routine to calcute quantities used for reactive particles
!
    k1 = k1_imn(imn)
    k2 = k2_imn(imn)

    call calc_conversion(fp(k1:k2,:))
    call calc_St(fp(k1:k2,:))
    call calc_mod_surf_area(fp(k1:k2,:))
!
    call calc_ads_entropy(fp(k1:k2,:))
    call calc_ads_enthalpy(fp(k1:k2,:))
    call calc_surf_entropy(fp(k1:k2,:))
    call calc_surf_enthalpy(fp(k1:k2,:))
!
    call calc_enthalpy_of_reaction()
    call calc_entropy_of_reaction()
!
    call calc_RR_hat(f,fp(k1:k2,:))
    call calc_ndot_mdot_R_j_hat(fp(k1:k2,:))
    call calc_R_c_hat(fp(k1:k2,:))
!
  end subroutine calc_pchem_factors
!***************************************************    
subroutine initialize_particles_surf(fp,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-sep-14/jonas coded
!  JONAS: needs to be filled with life
!
      real, dimension (mpar_loc,mpvar) :: fp
      logical :: lstarting
!
!  
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(lstarting)
!
    end subroutine initialize_particles_surf
!**************************************************************
    subroutine dpsurf_dt(f,df,fp,dfp,ineargrid)
!
!  evolution of particle surface fractions
!  (all particles on one node)
!
!  1-oct-14/Jonas: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      integer, dimension (mpar_loc,3) :: ineargrid
!
!  JONAS: equations.tex eq 37   
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpsurf_dt
!*****************************************************************
    subroutine dpsurf_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle temperature.
!  (all particles on one pencil)
!
!  23-sep-14/Nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp, dfp
      type (pencil_case) :: p
      integer, dimension (mpar_loc,3) :: ineargrid
      integer :: k

      intent (in) :: f, fp, ineargrid
      intent (inout) :: dfp, df
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
!  JONAS: equations.pdf eq. 37, look in 
!  particles_mass for looping and such
!
   do k=k1_imn(imn),k2_imn(imn)
!
!  JONAS: implicit?/explicit?
!  communicating with gas phase?
!
   end do
    endsubroutine dpsurf_dt_pencil
!***********************************************************************
  end module Particles_surfspec
