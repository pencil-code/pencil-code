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
  use Particles_cdata, only: iap, irhopswarm,pvarname,npvar
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
  character (len=labellen), dimension (ninit) :: init_surf='nothing'
  character(10), dimension(15)  :: species_name,adsorbed_species_names
  real, dimension(10) :: init_surf_ads_frac,init_surf_gas_frac
  real, dimension(:,:), allocatable :: R_j_hat
  real, dimension(:), allocatable :: R_c_hat
  real :: diffusivity=0.0
  real :: init_thCO=0.0
  real :: total_carbon_sites=1.08e-7
  integer :: iads=0, iads_end=0,isurf=0,isurf_end=0
!
  namelist /particles_ads_init_pars/ &
      init_adsorbed, init_surf, &
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
!  Set up indices for access to the fp and dfp arrays
!
!  29-aug-14/jonas: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      integer :: i,N_surface_species,N_adsorbed_species,stat,j
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
      call get_pchem_info(species,'N_surface_species',N_surface_species,'quiet')
!      
      if (nadsspec/=(N_adsorbed_species-1) .and. &
           nsurfreacspec/=N_surface_species .and. &
           N_adsorbed_species>0) then
         print*,'N_adsorbed_species: ',N_adsorbed_species
         print*,'N_surface_species: ', N_surface_species
         call fatal_error('register_particles_ads', & 
              'wrong size of storage for adsorbed and surface species allocated.')
         elseif (nadsspec/=(N_adsorbed_species-1) .and. &
           N_adsorbed_species>0) then
         print*,'N_adsorbed_species: ',N_adsorbed_species
         call fatal_error('register_particles_ads', & 
              'wrong size of storage for adsorbed species allocated.')
         elseif (nsurfreacspec/=N_surface_species) then
         print*,'N_surface_species: ', N_surface_species
         call fatal_error('register_particles_ads', &
              'wrong size of storage for surface species allocated')
         else
      endif
!      
      call get_species_list('solid_species',solid_species)
      call get_species_list('adsorbed_species_names',adsorbed_species_names)
!
!  Indices of adsorbed species mole fraction
!
      if (N_adsorbed_species>1) then
         j=1
         iads = npvar+1
         i=1
         do while (i<=(N_adsorbed_species-1))
            if (adsorbed_species_names(i)/='Cf') then
            pvarname(iads+j-1) = adsorbed_species_names(i)
            i = i+1
            j = j+1
            else
            i = i+1
            end if 
         enddo
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
!  Increase of npvar according to N_surface_species, which is
!  the concentration of gas phase species at the particle surface
!
      if (N_surface_species>1) then
         isurf = npvar+1
         do i=1,N_surface_species
            pvarname(isurf+i-1)=solid_species(i)
         enddo
         npvar=npvar+N_surface_species-1
         isurf_end=isurf+N_surface_species-1
      else
         call fatal_error('register_particles_ads', &
              'N_surface_species must be > 1')
      endif
!
!  Check that the fp and dfp arrays are big enough.
!
      if (npvar > mpvar) then
        if (lroot) write(0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
        call fatal_error('register_ads','npvar > mpvar')
      endif
!
    allocate(R_j_hat(mpar_loc,N_adsorbed_species-1)   ,STAT=stat)
    if (stat>0) call fatal_error('register_particles_ads',&
        'Could not allocate memory for R_j_hat')
    allocate(R_c_hat(mpar_loc)   ,STAT=stat)
    if (stat>0) call fatal_error('register_particles_ads',&
        'Could not allocate memory for R_c_hat')
!
    call register_indep_pchem()
    call register_dep_pchem()
!
    end subroutine register_particles_ads
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
         call transfer_initial_values(fp)
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
!
        endselect
!
!  Writing the initial adsorbed species fractions
!
        select case (init_surf(j))
        case ('nothing')
          if (lroot .and. j==1) print*, 'init_particles_ads,gas phase: nothing'
        case ('constant')
          if (lroot) print*, 'init_particles_ads: Initial Surface Fractions'
          do i=1,mpar_loc
             fp(i,isurf:isurf_end)=fp(i,isurf:isurf_end)+init_surf_gas_frac(1:)
          enddo
          case default
          if (lroot) &
              print*, 'init_particles_ads: No such such value for init_surf: ', &
              trim(init_surf(j))
          call fatal_error('init_particles_ads','')
!
        endselect
!
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
      real, dimension (mpar_loc) :: mod_surf_area
      type (pencil_case) :: p
      real, dimension(mpar_loc,3) :: ineargrid
      integer :: i, n_ads
!
      intent (inout) :: dfp
      intent (in) :: fp 
!
!
      call get_R_j_hat(R_j_hat)
      call get_mod_surf_area(mod_surf_area,fp)
      call get_R_c_hat(R_c_hat,fp)
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
  end module Particles_adsorbed
