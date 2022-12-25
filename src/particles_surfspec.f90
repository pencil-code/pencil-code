!$Id: particles_surfspec.f90 21950 2014-07-08 08:53:00Z jonas.kruger $
!
!  MOUDLE_DOC: This module takes care the gas phase species in the
!  MODULE_DOC: immediate vicinity of reactive particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
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
  use Particles_mpicomm
  use Particles_sub
  use Particles_chemistry
  use Chemistry
  use EquationOfState
  use SharedVariables
!
  implicit none
!
  include 'particles_surfspec.h'
!
!*********************************************************************!
!               Particle independent variables below here             !
!*********************************************************************!
!
  character(len=labellen), dimension(ninit) :: init_surf='nothing'
  character(len=10), dimension(:), allocatable :: solid_species
  real, dimension(10) :: init_surf_mol_frac ! INIT_DOC: Initial surface fraction
  real, dimension(:,:), allocatable :: nu_power
  integer, dimension(nchemspec) :: ispecaux=0
  integer :: ispecenth
  integer :: ndiffsteps=3 ! RUN_DOC: Number of mass transfer diffusion steps
  integer :: Ysurf
  logical :: lspecies_transfer=.true. ! RUN_DOC: species transfer between particle and gas
  logical :: linfinite_diffusion=.true. ! RUN_DOC: infinitely fast diffusion between particle and gas
  logical :: lboundary_explicit=.true. ! RUN_DOC: explicit evolution of particle surface species
  logical :: lpchem_cdtc=.false.
  logical :: lpchem_mass_enth=.true.
  logical :: lpfilter =.true.
  logical :: lwrite=.true.
  real :: rdiffconsts=0.1178
  logical :: ldiffuse_backspec=.false., ldiffs=.false.
  logical :: ldiffuse_backenth=.false., ldiffenth=.false.
!
  integer :: jH2, jO2, jCO2, jCO, jCH4, jN2, jH2O
  integer :: jOH, jAR, jO, jH, jCH, jCH2, jHCO, jCH3
  integer, dimension(:), allocatable :: idiag_surf
!
!*********************************************************************!
!               Particle dependent variables below here               !
!*********************************************************************!
!
  real, dimension(:,:), allocatable :: X_infty_reactants
  real, dimension(:,:), allocatable :: mass_trans_coeff_reactants
  real, dimension(:,:), allocatable :: mass_trans_coeff_species
  real, dimension(:,:,:), allocatable :: weight_array, dmass_frac_dt
!
  namelist /particles_surf_init_pars/ init_surf, init_surf_mol_frac
!
  namelist /particles_surf_run_pars/ &
      lboundary_explicit, linfinite_diffusion, lspecies_transfer, &
      lpchem_mass_enth, ldiffuse_backspec, ldiffs,rdiffconsts, &
      ndiffsteps,ldiffuse_backenth,ldiffenth,lpfilter
!
  integer :: idiag_dtpchem=0   ! DIAG_DOC: $dt_{particle,chemistry}$
!
  contains
! ******************************************************************************
!  This is a wrapper routine for particle dependent and particle
!  independent variables
!  JONAS: Back to standalone via mpar_loc=1?
!
    subroutine register_particles_surfspec()
!
      use FArrayManager, only: farray_register_auxiliary
!
!
      character(len=11) :: chemspecaux
      integer :: i
!
!      if (lroot) call svn_id( &
!          "$Id: particles_surfspec.f90 20849 2014-10-06 18:45:43Z jonas.kruger $")
!
      call register_indep_psurfspec()
      call register_dep_psurfspec()
!
!  We need to register an auxiliary array to diffuse the species transfer
!
      if (ldiffuse_backspec .and. lspecies_transfer) then
        chemspecaux = 'ichemspec  '
        do i = 1,nchemspec
          write (chemspecaux(10:11),'(I2)') i
          call farray_register_auxiliary(chemspecaux,ispecaux(i),communicated=.true.)
        enddo
      elseif (ldiffuse_backspec .and. .not. lspecies_transfer) then
        call fatal_error('particles_surfspec:', &
            'diffusion of the species transfer needs lspecies_transfer')
      endif
!
!  Diffusion of mass bound enthalpy
!
      if (ldiffuse_backenth .and. lspecies_transfer) then
        call farray_register_auxiliary('ispecenth',ispecenth,communicated=.true.)
      elseif (ldiffuse_backenth .and. .not. lspecies_transfer) then
        call fatal_error('particles_surfspec:', &
            'diffusion of the mass bound enthalpy needs lspecies_transfer')
      endif
!
    endsubroutine register_particles_surfspec
! ******************************************************************************
!
    subroutine register_indep_psurfspec()
      integer :: i, k, stat, k1, k2
      character(len=10), dimension(40) :: reactants
!
      if (nsurfreacspec /= N_surface_species) then
        print*,'N_surface_species: ', N_surface_species
        print*,'NSURFREACSPEC :', nsurfreacspec
        call fatal_error('register_particles_surf', &
            'wrong size of storage for surface species allocated')
      endif
!
!  Increase of npvar according to N_surface_species, which is
!  the concentration of gas phase species at the particle surface
!
      if (N_surface_species <= 1) then
        call fatal_error('register_particles_', 'N_surface_species must be > 1')
      endif
!
! Allocate memory for a number of arrays
      allocate(dependent_reactant(N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_psurfchem', &
          'Could not allocate memory for dependent_reactant')
      allocate(nu(N_surface_species,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_psurfchem', &
          'Could not allocate memory for nu')
      allocate(nu_prime(N_surface_species,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_psurfchem', &
          'Could not allocate memory for nu_prime')
      allocate(ac(N_surface_species),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_psurfchem', &
          'Could not allocate memory for ac')
      allocate(jmap(N_surface_species),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_psurfchem', &
          'Could not allocate memory for jmap')
      allocate(solid_species(N_surface_species),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_psurfchem', &
          'Could not allocate memory for solid_species')
      allocate(nu_power(N_surface_species,N_surface_reactions),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_chem', &
          'Could not allocate memory for nu_power')
      allocate(idiag_surf(N_surface_species),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_chem', &
          'Could not allocate memory for idiag_surf')
      idiag_surf = 0
!
      call create_ad_sol_lists(solid_species,'sol')
      call get_reactants(reactants)
      call sort_compounds(reactants,solid_species,n_surface_species)
!
      do i = 1,N_surface_species
        call append_npvar('i'//solid_species(i),isurf)
      enddo
      isurf_end = isurf
!
!  create binding between gas phase and near field gas species
!  chemistry and particles_chemistry module
!
      call create_jmap()
!
!  print*, jH2O,JCO2,jH2,jO2,jCO,jCH,jHCO,jCH2,jCH3
!
! Set number of carbon atoms for each surface species
!
      call get_ac(ac,solid_species,N_surface_species)
!
! Set the stoichiometric matrixes
!      print*, 'Number of surface species: ',N_surface_species
!
      call create_stoc(solid_species,nu,.true., N_surface_species,nu_power)
      call create_stoc(solid_species,nu_prime,.false., N_surface_species)
!
! Define which gas phase reactants the given reaction depends on
      call create_dependency(nu,dependent_reactant, &
          n_surface_reactions,n_surface_reactants)
!
! Find the mole production of the forward reaction
      call create_dngas()
!
      if (lwrite) then
        call write_outputfile()
      endif
      lwrite = .false.
!
    endsubroutine register_indep_psurfspec
! ******************************************************************************
    subroutine register_dep_psurfspec()
    endsubroutine register_dep_psurfspec
! ******************************************************************************
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-sep-14/jonas coded
!
    subroutine initialize_particles_surf(f)
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: dimx, dimy, dimz, ncells=1
!
!      print*, weight_array
!
      if (lparticlemesh_gab) ncells = 7
      if (lparticlemesh_tsc) ncells = 3
      if (lparticlemesh_cic) ncells = 2
!
      if (nxgrid /= 1) then
        dimx = ncells
      else
        dimx = 1
      endif
      if (nygrid /= 1) then
        dimy = ncells
      else
        dimy = 1
      endif
      if (nzgrid /= 1) then
        dimz = ncells
      else
        dimz = 1
      endif
!
      if (allocated(weight_array)) deallocate(weight_array)
      if (lparticlemesh_gab .or. lparticlemesh_tsc .or. lparticlemesh_cic) then
        allocate(weight_array(dimx,dimy,dimz))
      endif
      if (.not. allocated(weight_array)) then
        allocate(weight_array(1,1,1))
      endif
      call precalc_weights(weight_array)
!
      call keep_compiler_quiet(f)
    endsubroutine initialize_particles_surf
! ******************************************************************************
    subroutine read_particles_surf_init_pars(iostat)
      use File_io, only: parallel_unit
      integer, intent(out) :: iostat
!
      read (parallel_unit, NML=particles_surf_init_pars, IOSTAT=iostat)
    endsubroutine read_particles_surf_init_pars
! ******************************************************************************
    subroutine write_particles_surf_init_pars(unit)
      integer, intent(in) :: unit
!
      write (unit, NML=particles_surf_init_pars)
    endsubroutine write_particles_surf_init_pars
! ******************************************************************************
    subroutine read_particles_surf_run_pars(iostat)
      use File_io, only: parallel_unit
      integer, intent(out) :: iostat
!
      read (parallel_unit, NML=particles_surf_run_pars, IOSTAT=iostat)
    endsubroutine read_particles_surf_run_pars
! ******************************************************************************
    subroutine write_particles_surf_run_pars(unit)
      integer, intent(in) :: unit
!
      write (unit, NML=particles_surf_run_pars)
    endsubroutine write_particles_surf_run_pars
! ******************************************************************************
!  Initial particle surface fractions
!
!  01-sep-14/jonas: coded
!
    subroutine init_particles_surf(f,fp,ineargrid)
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
      real :: sum_surf_spec
      real :: mean_molar_mass
      integer :: i, j, k, ix0, igas, iy0, iz0
!
      intent(in) :: f
      intent(out) :: fp
!
      fp(:,isurf:isurf_end) = 0.
      do j = 1,ninit
!
! Writing the initial surface species fractions
        select case (init_surf(j))
!
        case ('nothing')
          if (lroot .and. j == 1) print*, 'init_particles_surf,gas phase: nothing'
!
        case ('constant')
          if (lroot) print*, 'init_particles_surf: Initial Surface Fractions'
!
! This ensures that we don't have unphysical values as init
          sum_surf_spec = sum(init_surf_mol_frac(1:N_surface_species))
          if (sum_surf_spec > 1) then
            if (lroot)  print*, 'Sum of all surface fractions >1, normalizing...'
            init_surf_mol_frac(1:N_surface_species) = &
                init_surf_mol_frac(1:N_surface_species) / sum_surf_spec
          endif
!
          if (sum_surf_spec == 0.0) then
            call fatal_error('particles_surfspec', &
                'initial molefraction was given without value. '// &
                'please specify the initial gas fraction')
          endif
!
          do k = 1,mpar_loc
            fp(k,isurf:isurf_end) = init_surf_mol_frac(1:N_surface_species)
          enddo
        case ('gas')
!
!  The starting particle surface mole fraction is equal to the
!  gas phase composition at the particles position.
!  This is not functional, and this would need the initialization of
!  The gas field before
          do k = 1,mpar_loc
            mean_molar_mass = 0.0
!
            do i = 1,nchemspec
              mean_molar_mass = mean_molar_mass + &
                  species_constants(i,imass) * f(4,4,4,ichemspec(i))
            enddo
!
!
            do i = 1, N_surface_species
              igas = ichemspec(jmap(i))
              fp(k,isurf+i-1) = f(4,4,4,igas) / &
                  species_constants(jmap(i),imass)*mean_molar_mass
            enddo
          enddo
!
        case default
          if (lroot) &
              print*, 'init_particles_ads: No such such value for init_surf: ', &
              trim(init_surf(j))
          call fatal_error('init_particles_surf','')
        endselect
      enddo
    endsubroutine init_particles_surf
! ******************************************************************************
!  evolution of particle surface fractions
!  (all particles on one node)
!
!  1-oct-14/Jonas: coded
!
    subroutine dpsurf_dt(f,df,fp,dfp,ineargrid)
!
      use Boundcond
      use Mpicomm
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      real :: dt1
      integer, dimension(mpar_loc,3) :: ineargrid
      integer :: i, j,g
      integer :: li,mi,ni
!
!  equations.tex eq 37
!      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
      dt1 = 1./dt
!
!  Diffuse the transfer of species from particle to fluid, and adapt
!  the bulk species (ichemspec(nchemspec)) to ensure species conservation
!
      if (lspecies_transfer .and. ldiffuse_backspec) then
        if (ldensity_nolog) call fatal_error('particles_surf', &
            'not implemented for ldensity_nolog')
        do j = 1,nchemspec-1
          do i = 1,ndiffsteps
            call boundconds_x(f,ispecaux(j),ispecaux(j))
            call initiate_isendrcv_bdry(f,ispecaux(j),ispecaux(j))
            call finalize_isendrcv_bdry(f,ispecaux(j),ispecaux(j))
            call boundconds_y(f,ispecaux(j),ispecaux(j))
            call boundconds_z(f,ispecaux(j),ispecaux(j))
!
            call diffuse_interaction(f(:,:,:,ispecaux(j)),ldiffs,.False.,rdiffconsts)
          enddo
!
!  Control that we don't influence points that would bring the mass fraction into the negative 
!
          if (.not. lpfilter) then
            df(l1:l2,m1:m2,n1:n2,ichemspec(j)) =  df(l1:l2,m1:m2,n1:n2,ichemspec(j)) + &
                f(l1:l2,m1:m2,n1:n2,ispecaux(j))
          else
            where ((f(l1:l2,m1:m2,n1:n2,ichemspec(j))+ &
                f(l1:l2,m1:m2,n1:n2,ispecaux(j))*dt)< 0.0)
              f(l1:l2,m1:m2,n1:n2,ispecaux(j)) = f(l1:l2,m1:m2,n1:n2,ichemspec(j))*dt1*(-0.9)
              where ((f(l1:l2,m1:m2,n1:n2,ichemspec(j)) <= 1E-25))
                f(l1:l2,m1:m2,n1:n2,ispecaux(j)) = 0.0
              end where
            end where
            df(l1:l2,m1:m2,n1:n2,ichemspec(j)) =  df(l1:l2,m1:m2,n1:n2,ichemspec(j)) + &
                f(l1:l2,m1:m2,n1:n2,ispecaux(j))
          endif
        enddo
        df(l1:l2,m1:m2,n1:n2,ichemspec(nchemspec)) = &
            df(l1:l2,m1:m2,n1:n2,ichemspec(nchemspec))- &
            sum(df(l1:l2,m1:m2,n1:n2,ichemspec(:)),DIM=4)
      endif
!
!  Diffusion of the mass bound enthalpy
!
      if (lspecies_transfer .and. ldiffuse_backenth) then
        if (ldensity_nolog .and. ltemperature_nolog) &
            call fatal_error('particles_surf', 'diffusion of mass bound enthalpy &
            & not implemented for ldensity_nolog')
        do i = 1,ndiffsteps
          call boundconds_x(f,ispecenth,ispecenth)
          call initiate_isendrcv_bdry(f,ispecenth,ispecenth)
          call finalize_isendrcv_bdry(f,ispecenth,ispecenth)
          call boundconds_y(f,ispecenth,ispecenth)
          call boundconds_z(f,ispecenth,ispecenth)
          call diffuse_interaction(f(:,:,:,ispecenth),ldiffenth,.False.,rdiffconsts)
!
        enddo
        df(l1:l2,m1:m2,n1:n2,ilnTT) =  df(l1:l2,m1:m2,n1:n2,ilnTT) + &
            f(l1:l2,m1:m2,n1:n2,ispecenth)
      endif
!
      if (ldiagnos) then
        do i = 1,N_surface_species
          if (idiag_surf(i) /= 0) then
            call sum_par_name(fp(1:npar_loc,isurf+i-1),idiag_surf(i))
          endif
        enddo
      endif
!
    endsubroutine dpsurf_dt
! ******************************************************************************
    subroutine dpsurf_dt_pencil(f,df,fp,dfp,p,ineargrid)
!
!  Evolution of particle surface fraction.
!  (all particles on one pencil)
!
!  23-sep-14/Nils: coded
!
      use Diagnostics
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(mpar_loc,mparray) :: fp
      real, dimension(mpar_loc,mpvar) :: dfp
      real, dimension(nx,nchemspec) :: chem_reac
      real, dimension(:,:), allocatable :: term, ndot
      real, dimension(:), allocatable :: Cg_surf, mass_loss
      real :: porosity
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3) :: ineargrid
      integer :: k, k1, k2, i, ix0, iy0, iz0
      real :: weight, volume_cell, rho1_point
      real :: mean_molar_mass, dmass, A_p, denth
      real :: reac_pchem_weight, max_reac_pchem
      real :: summan, m1_cell
      real :: dmass_ndot
      integer :: ix1, iy1, iz1, ierr
      integer :: ixx, iyy, izz
      integer :: ixx0, iyy0, izz0
      integer :: ixx1, iyy1, izz1
      integer :: index1, index2
      real, pointer :: true_density_carbon
      integer :: density_index
      real :: dmass_frac_dt=0.0
      real :: diffusion_transfer=0.0
!
      intent(in) :: ineargrid
      intent(inout) :: dfp, df, fp,f
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
!  Set the particle reaction time to 0 if no particles are present
!
      reac_pchem = 0.0
!
      if (lpreactions) then
!  initializing the auxiliary pencils for the mass and mass bound enthalpy diffusi
!
        volume_cell = (lxyz(1)*lxyz(2)*lxyz(3))/(nxgrid*nygrid*nzgrid)
        if (ldiffuse_backenth) f(l1:l2,m,n,ispecenth) = 0.0
        if (ldiffuse_backspec) then
          do i = 1,nchemspec
            f(l1:l2,m,n,ispecaux(i)) = 0.0
          enddo
        endif
!
!  Do only if particles are present on the current pencil
!
        if (npar_imn(imn) /= 0) then
!
          k1 = k1_imn(imn)
          k2 = k2_imn(imn)
!
          allocate(Cg_surf(k1:k2))
          allocate(mass_loss(k1:k2))
          allocate(ndot(k1:k2,N_surface_species))
          allocate(term(k1:k2,N_surface_reactants))
!
          call get_shared_variable('true_density_carbon',true_density_carbon,caller='dpsurf_dt_pencil')
!
          call calc_mass_trans_reactants()
          call get_surface_chemistry(Cg_surf,ndot,mass_loss)
!
!  set surface gas species composition
!  (infinite diffusion, set as far field and convert from
!  mass fractions to mole fractions)
!
          do k = k1,k2
!
!  particle on pencil loop
!
            porosity = 1.0 - (fp(k,imp)/(fp(k,iap)**3*4./3.*pi))/true_density_carbon
            A_p = 4.*pi*(fp(k,iap)**2)
            mean_molar_mass = (interp_rho(k)*Rgas*interp_TT(k)/ interp_pp(k))
!
            if (lboundary_explicit) then
              if (lparticles_adsorbed) print*, 'values in surfspec begin',fp(k,isurf:isurf_end)
              do i = 1,N_surface_reactants
                diffusion_transfer = mass_trans_coeff_reactants(k,i) * (interp_species(k,jmap(i)) / &
                    species_constants(jmap(i),imass) * mean_molar_mass-fp(k,isurf+i-1))
                term(k,i) = ndot(k,i) - fp(k,isurf+i-1) * sum(ndot(k,:)) + diffusion_transfer
                if (lparticles_adsorbed) then
                  print*, '---------------------------'
                  print*, ' mass_trans_', mass_trans_coeff_reactants(k,i)
                  print*, ' ndot(k,i)  ', ndot(k,i)
                  print*, 'fp(k,isurf+)',fp(k,isurf+i-1)
                  print*, 'sum(ndot(k,:',sum(ndot(k,:))
                  print*, 'x_spec      ',interp_species(k,jmap(i)) / species_constants(jmap(i),imass) * &
                      mean_molar_mass
                  print*, '---------------------------'
                endif
!
! the term 3/fp(k,iap) is ratio of the surface of a sphere to its volume
!
                dfp(k,isurf+i-1) = dfp(k,isurf+i-1) + 3*term(k,i)/(porosity*Cg_surf(k)*fp(k,iap))
                if (lparticles_adsorbed) then
                  print*, 'term(k,i)',term(k,i)
                  print*, 'dfp(k,isurf+i-1): ',k,i,dfp(k,isurf+i-1)
                endif
              enddo
!              if (lparticles_adsorbed) stop
            else
              if (linfinite_diffusion .or. lbaum_and_street) then
                do i = 1,N_surface_reactants
                  fp(k,isurf+i-1) = interp_species(k,jmap(i)) / &
                      species_constants(jmap(i),imass) * mean_molar_mass
                  dfp(k,isurf+i-1) = 0.
                enddo
              else
                print*,'Must set linfinite_diffusion=T if lboundary_explicit=F.'
                call fatal_error('dpsurf_dt_pencil', &
                    'Implicit solver for surface consentrations is not implemented.')
              endif
            endif
!
!  the following block is thoroughly commented in particles_temperature
!  find values for transfer of variables from particle to fluid
            if (lspecies_transfer) then
              ix0 = ineargrid(k,1)
              iy0 = ineargrid(k,2)
              iz0 = ineargrid(k,3)
              call find_interpolation_indeces(ixx0,ixx1,iyy0,iyy1,izz0,izz1, &
                  fp,k,ix0,iy0,iz0)
!
! positive dmass means particle is losing mass
! jmap: gives the ichemspec of a surface specie
              dmass = 0.
              denth = 0.
!
              do i = 1,N_surface_species
                index1 = jmap(i)
                dmass = dmass-ndot(k,i)*species_constants(index1,imass)*A_p
                if (lpchem_mass_enth) then
                  if (ndot(k,i) > 0.0) then
!
!  Species enthalpy at the particle phase temperature
!  Qenth = ndot(mol/s/m2)*A(m2)*enth(J/mol)=J/s
!  to convert to erg, the value has to be multiplied with 1e7 for values
!  fetched from particles_chemistry
!
!  Enthalpy from Pencil code at particles temperature
!
                    if (lpencil(i_H0_RT)) then
                      denth = denth+ndot(k,i)*A_p*p%cv(ix0-nghost)*(fp(k,iTp)-298.15)
                    else
                      call fatal_error('particles_surfspec','mass bound enthalpy transfer needs p%H0_RT')
                    endif
                  else
!
! Species enthalpy at the gas phase temperature
!
                    if (lpencil(i_H0_RT)) then
                      denth = denth+ndot(k,i)*A_p*p%cv(ix0-nghost)*(p%TT(ix0-nghost)-298.15)
                    else
                      call fatal_error('particles_surfspec', &
                          'mass bound enthalpy transfer needs p%H0_RT')
                    endif
                  endif
                endif
              enddo
!
!  Prepare the max_reac_pchem
!
              if (lfirst .and. ldt) max_reac_pchem = 0.0
!
!  Sum for testing
!
              summan = 0.0
!
!  The whole following with find_index is a workaround
!  to enable the accessing of nonreacting gas phase species by this routine
!  since they are as well affected when the particle is adding species to the
!  gas phase.
!  find_index maps the i from the gas phase chemistry to the surface_species
!  from particles_chemistry which is needed to access ndot, since this
!  only contains surface species
!
              if (ldensity_nolog) then
                if (.not. ldiffuse_backspec) then
                  do i = 1,nchemspec-1
                    reac_pchem_weight = 0.0
                    if (find_index(i,jmap,N_surface_species) > 0 ) then
!NILS: Isn't index1=i?
!                      index1 = jmap(find_index(i,jmap,N_surface_species))
                      index2 = ichemspec(i)
                      df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,index2) = &
                          df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,index2)+(A_p* &
                          ndot(k,find_index(i,jmap,N_surface_species))* &
                          species_constants(i,imass)+ &
                          dmass*interp_species(k,i))*weight_array/ &
                          (f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho)*volume_cell)
                      if (lfirst .and. ldt) then
                        reac_pchem_weight = max(reac_pchem_weight, &
                            abs(maxval(df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,index2))/ &
                            max(maxval(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,index2)),1e-10)))
                      endif
                    else
                      df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ichemspec(i)) = &
                          df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ichemspec(i)) + dmass * &
                          interp_species(k,i)*weight_array/ &
                          (f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho)*volume_cell)
                    endif
                  enddo
                else
                  call fatal_error('particles_surf','backdiffusion not yet implemented &
                      & for ldensity_nolog')
                endif
              else
!
!  For density_log, there is the possibility to diffuse the species influence
!  this happens in the else block of the next if clause
!
                if (.not. ldiffuse_backspec) then
                  do i = 1,nchemspec-1
                    reac_pchem_weight = 0.0
                    if (find_index(i,jmap,N_surface_species) > 0 ) then
!NILS: Isn't index1=i?
!                      index1 = jmap(find_index(i,jmap,N_surface_species))
                      index2 = ichemspec(i)
                      df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,index2) = &
                          df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,index2)+(A_p* &
                          ndot(k,find_index(i,jmap,N_surface_species))* &
                          species_constants(i,imass)+ &
                          dmass*interp_species(k,i))*weight_array/ &
                          (exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))*volume_cell)
                      if (lfirst .and. ldt) then
                        reac_pchem_weight = max(reac_pchem_weight, &
                            abs(max(maxval(df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,index2)/ &
                            f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,index2)),1e-10)))
                      endif
                    else
                      df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ichemspec(i)) = &
                          df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ichemspec(i)) + &
                          dmass*interp_species(k,i)*weight_array/ &
                          (exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))*volume_cell)
                    endif
!                  summan = summan+dmass_frac_dt
                  enddo
                else
                  do i = 1,nchemspec-1
                    if (find_index(i,jmap,N_surface_species) > 0 ) then
!  NILS: Isn't index1=i?
!                      index1 = jmap(find_index(i,jmap,N_surface_species))
                      f(ix0,iy0,iz0,ispecaux(i)) = &
                          f(ix0,iy0,iz0,ispecaux(i))+(A_p* &
                          ndot(k,find_index(i,jmap,N_surface_species))* &
                          species_constants(i,imass)+ &
                          dmass*interp_species(k,i))/ &
                          (exp(f(ix0,iy0,iz0,ilnrho))*volume_cell)
                    else
                      f(ix0,iy0,iz0,ispecaux(i)) = &
                          f(ix0,iy0,iz0,ispecaux(i)) + &
                          dmass*interp_species(k,i)/ &
                          (exp(f(ix0,iy0,iz0,ilnrho))*volume_cell)
                    endif
                  enddo
                endif
              endif
!
!  Solving for all but the other values, setting the last one to the
!  negative values of all others.
!
              if (.not. ldiffuse_backspec) then
                df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ichemspec(nchemspec)) = &
                    df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ichemspec(nchemspec))- &
                    sum(df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ichemspec(:)),DIM=4)
              endif
!
!  Compare the current maximum reaction rate to the previous one
!
              if (lfirst .and. ldt) max_reac_pchem = &
                  max(max_reac_pchem, reac_pchem_weight)
!
!  Enthalpy transfer via mass transfer!
!
              if (lpchem_mass_enth .and. .not. ldiffuse_backenth) then
                if (ldensity_nolog) then
                  if (ltemperature_nolog) then
                    df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,iTT) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,iTT) &
                        +denth*p%cv1(ix0-nghost)*weight_array/ &
                        (f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho)*volume_cell)
                  else
                    df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnTT) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnTT) &
                        +denth*p%cv1(ix0-nghost)*p%TT1(ix0-nghost)*weight_array/ &
                        (f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,irho)*volume_cell)
                  endif
                else
                  if (ltemperature_nolog) then
                    df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,iTT) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,iTT) &
                        +denth*p%cv1(ix0-nghost)*weight_array/ &
                        (exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))*volume_cell)
                  else
                    df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnTT) = df(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnTT) &
                        +denth*p%cv1(ix0-nghost)*p%TT1(ix0-nghost)*weight_array/ &
                        (exp(f(ixx0:ixx1,iyy0:iyy1,izz0:izz1,ilnrho))*volume_cell)
                  endif
                endif
!
!  Diffusion of mass bound enthalpy transfer
!
              elseif (lpchem_mass_enth .and. ldiffuse_backenth) then
                if (ldensity_nolog .and. ltemperature_nolog) then
                  call fatal_error('particles_surf','diffusion not implemented for &
                      & ldensity_nolog')
                elseif (ldensity_nolog .and. .not. ltemperature_nolog) then
                  call fatal_error('particles_surf','diffusion not implemented for &
                      & ldensity_nolog')
                elseif (.not. ldensity_nolog .and. ltemperature_nolog) then
                  call fatal_error('particles_surf','diffusion not implemented for &
                      & ltemperature_nolog')
                else
                  f(ix0,iy0,iz0,ispecenth) =  f(ix0,iy0,iz0,ispecenth) &
                      +denth*p%cv1(ix0-nghost)*p%TT1(ix0-nghost) / &
                      (exp(f(ix0,iy0,iz0,ilnrho))*volume_cell)
                endif
              endif
!
!  Compare the current maximum reaction rate to the current maximum
!  reaction rate in the current pencil
!
              if (lfirst .and. ldt) reac_pchem = max(reac_pchem,max_reac_pchem)
!
            endif
!
!  end of particle on pencil loop
!
!            if (lparticles_adsorbed) print*, 'values in surfspec end',fp(k,isurf:isurf_end)
          enddo
!
!
          if (ldiagnos) then
            if (idiag_dtpchem /= 0 ) call max_name(reac_pchem/cdtc,idiag_dtpchem,l_dt=.true.)
          endif
!
          if (allocated(term)) deallocate(term)
          if (allocated(ndot)) deallocate(ndot)
          if (allocated(Cg_surf))   deallocate(Cg_surf)
          if (allocated(mass_loss))   deallocate(mass_loss)
        endif
!
      endif
!
    endsubroutine dpsurf_dt_pencil
! ******************************************************************************
!  Read and register print parameters relevant for
!  particles near field gas composition
!
!  06-oct-14/jonas: adapted
!
    subroutine rprint_particles_surf(lreset,lwrite)
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
      integer :: iname,i
      character(len=7) :: diagn_surf, number
!
! Write information to index.pro
      lwr = .false.
      if (present(lwrite)) lwr = lwrite
      if (lreset) then
        idiag_dtpchem = 0
        idiag_surf = 0
      endif
!
      if (lroot .and. ip < 14) print*,'rprint_particles_surf: run through parse list'
!
      do iname = 1,nname
        do i = 1,N_surface_species
          write (number,'(I2)') i
          diagn_surf = 'Ysurf'//trim(adjustl(number))
          call parse_name(iname,cname(iname),cform(iname),trim(diagn_surf),idiag_surf(i))
        enddo
        call parse_name(iname,cname(iname),cform(iname),'dtpchem',idiag_dtpchem)
      enddo
!
    endsubroutine rprint_particles_surf
! ******************************************************************************
!  07-oct-2014/jonas: coded
!  29-oct-2014/jonas: moved parts of the code into routine, implemented
!                     error when gas phase species is in
!                     particles_chemistry, but not in chemistry
!
!  find the indexes used in chemistry.f90 of species present
!  in the near field of the particle
!
    subroutine create_jmap()
      use EquationOfState
!
      integer :: index_glob, index_chem
      logical :: found_species=.false.
!
      inuH2O = find_species('H2O',solid_species,n_surface_species)
      inuCO2 = find_species('CO2',solid_species,n_surface_species)
      inuH2 = find_species('H2',solid_species,n_surface_species)
      inuO2 = find_species('O2',solid_species,n_surface_species)
      inuCO = find_species('CO',solid_species,n_surface_species)
      inuCH = find_species('CH',solid_species,n_surface_species)
      inuHCO = find_species('HCO',solid_species,n_surface_species)
      inuCH2 = find_species('CH2',solid_species,n_surface_species)
      inuCH3 = find_species('CH3',solid_species,n_surface_species)
!
      call find_species_index('H2',index_glob,index_chem,found_species)
      if (found_species) then
        jH2 = index_chem
      else
        if (inuH2 > 0) call fatal_error('create_jmap','no H2 found')
      endif
!
      call find_species_index('O2',index_glob,index_chem,found_species)
      if (found_species) then
        jO2 = index_chem
      else
        if (inuO2 > 0) call fatal_error('create_jmap','no O2 found')
      endif
!
      call find_species_index('CO2',index_glob,index_chem,found_species)
      if (found_species) then
        jCO2 = index_chem
      else
        if (inuCO2 > 0) call fatal_error('create_jmap','no CO2 found')
      endif
!
      call find_species_index('CO',index_glob,index_chem,found_species)
      if (found_species) then
        jCO = index_chem
      else
        if (inuCO > 0) call fatal_error('create_jmap','no CO found')
      endif
!
      call find_species_index('CH4',index_glob,index_chem,found_species)
      if (found_species) then
        jCH4 = index_chem
      else
        if (inuCH4 > 0) call fatal_error('create_jmap','no CH4 found')
      endif
!
      call find_species_index('H2O',index_glob,index_chem,found_species)
      if (found_species) then
        jH2O = index_chem
      else
        if (inuH2O > 0) call fatal_error('create_jmap','no H2O found')
      endif
!
      call find_species_index('CH3',index_glob,index_chem,found_species)
      if (found_species) then
        jCH3 = index_chem
      else
        if (inuCH3 > 0) call fatal_error('create_jmap','no CH3 found')
      endif
!
      call find_species_index('CH',index_glob,index_chem,found_species)
      if (found_species) then
        jCH = index_chem
      else
        if (inuCH > 0) call fatal_error('create_jmap','no CH found')
      endif
!
      call find_species_index('CH2',index_glob,index_chem,found_species)
      if (found_species) then
        jCH2 = index_chem
      else
        if (inuCH2 > 0) call fatal_error('create_jmap','no CH2 found')
      endif
!      call find_species_index('HCO',index_glob,index_chem,found_species)
      if (found_species) then
        jHCO = index_chem
      else
        if (inuHCO > 0)  call fatal_error('create_jmap','no HCO found')
      endif
!
      if (inuH2O > 0)     jmap(inuH2O) = jH2O
      if (inuCO2 > 0)     jmap(inuCO2) = jCO2
      if (inuH2 > 0)      jmap(inuH2) = jH2
      if (inuO2 > 0)      jmap(inuO2) = jO2
      if (inuCO > 0)      jmap(inuCO) = jCO
      if (inuCH > 0)      jmap(inuCH) = jCH
      if (inuHCO > 0)     jmap(inuHCO) = jHCO
      if (inuCH2 > 0)     jmap(inuCH2) = jCH2
      if (inuCH3 > 0)     jmap(inuCH3) = jCH3
    endsubroutine create_jmap
! ******************************************************************************
!  restrict mass trans coefficients to surface reactants only
!
    subroutine calc_mass_trans_reactants()
      integer :: k, i, k1, k2
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      do k = k1,k2
        do i = 1, N_surface_reactants
          mass_trans_coeff_reactants(k,i) = mass_trans_coeff_species(k,i)
        enddo
      enddo
    endsubroutine calc_mass_trans_reactants
! ******************************************************************************
!  Allocate pencils for mass transport coefficiens
!
!  nov-14/jonas: coded
!
    subroutine allocate_surface_pencils()
      integer :: k1, k2, stat
!
      k1 = k1_imn(imn)
      k2 = k2_imn(imn)
!
      allocate(mass_trans_coeff_species(k1:k2,N_species), STAT=stat)
      if (stat > 0) call fatal_error('allocate_variable_pencils', &
          'Could not allocate memory for mass_trans_coeff_species')
      allocate(mass_trans_coeff_reactants(k1:k2,N_surface_reactants),STAT=stat)
      if (stat > 0) call fatal_error('register_indep_psurfchem', &
          'Could not allocate memory for mass_trans_coeff_reactants')
!
    endsubroutine allocate_surface_pencils
! ******************************************************************************
!
!  Clean up surface pencils
!
!  nov-14/jonas: coded
!
    subroutine cleanup_surf_pencils()
!
      deallocate(mass_trans_coeff_species)
      deallocate(mass_trans_coeff_reactants)
!
    endsubroutine cleanup_surf_pencils
! ******************************************************************************
!
!  Calculate pencils used for surface gas phase fraction evolution
!
!  nov-14/jonas: coded
!
    subroutine calc_psurf_pencils(f,fp,p,ineargrid)
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      type (pencil_case) :: p
!
      call allocate_surface_pencils()
      call calc_mass_trans_coeff(f,fp,p,ineargrid)
!
    endsubroutine calc_psurf_pencils
! ******************************************************************************
!  diffusion coefficients of gas species at nearest grid point
!
!  oct-14/Jonas: coded
!
    subroutine calc_mass_trans_coeff(f,fp,p,ineargrid)
      real, dimension(mx,my,mz,mparray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      type (pencil_case) :: p
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      integer :: k, k1, k2,i
      integer :: ix0
      integer :: spec_glob, spec_chem
      real :: diff_coeff_species
      real, dimension(:), allocatable :: Cg_diff
!
      if (npar_imn(imn) /= 0) then
        k1 = k1_imn(imn)
        k2 = k2_imn(imn)
!
!        allocate(diff_coeff_species(k1:k2,N_species))
        allocate(Cg_diff(k1:k2))
!
!        do k = k1,k2
!          ix0 = ineargrid(k,1)
!          do i = 1,N_surface_species
!            diff_coeff_species(k,i) = p%Diff_penc_add(ix0-nghost,jmap(i))
!          enddo
!        enddo
!
        do k = k1,k2
          Cg_diff(k) = k
        enddo

        call get_surface_chemistry(Cg_diff(:))
        do k = k1,k2
          ix0 = ineargrid(k,1)
          do i = 1,N_surface_species
            diff_coeff_species = p%Diff_penc_add(ix0-nghost,jmap(i))
            mass_trans_coeff_species(k,i) = Cg_diff(k)*diff_coeff_species/ fp(k,iap)
          enddo
        enddo
!
        if (allocated(Cg_diff)) deallocate(Cg_diff)
      endif
!
    endsubroutine calc_mass_trans_coeff
! ******************************************************************************
    subroutine particles_surfspec_clean_up()
!
      if (allocated(dependent_reactant)) deallocate(dependent_reactant)
      if (allocated(nu)) deallocate(nu)
      if (allocated(nu_prime)) deallocate(nu_prime)
      if (allocated(ac)) deallocate(ac)
      if (allocated(jmap)) deallocate(jmap)
      if (allocated(solid_species)) deallocate(solid_species)
      if (allocated(nu_power)) deallocate(nu_power)
      if (allocated(idiag_surf)) deallocate(idiag_surf)
!
    endsubroutine particles_surfspec_clean_up
! ******************************************************************************
    function find_index(element, list,lengthlist)
!
      implicit none
!
      integer :: find_index
      integer :: element,i
      integer, dimension(:) :: list
      integer :: lengthlist
!
      find_index = 0
!
      do i = 1,lengthlist
        if (element == list(i)) find_index = i
      enddo
!
    endfunction find_index
!***********************************************************************
    subroutine write_outputfile()
!
!  Write particle chemistry info to ./data/particle_chemistry.out
!
      integer :: file_id=123
      integer :: k
      character(len=20) :: writeformat
      character(len=30) :: output='./data/particle_chemistry.out'
!
      open (file_id,file=output,position='append')
!
      writeformat = '(  A8)'
      write (writeformat(2:3),'(I2)') N_surface_species
      write (file_id,*) 'Gas phase surface species'
      write (file_id,writeformat) solid_species
      write (file_id,*) ''
!
      writeformat = '(I2,  F5.2)'
      write (writeformat(5:6),'(I2)') N_surface_species
      write (file_id,*) 'Nu'
      do k = 1,N_surface_reactions
        write (file_id,writeformat) k,nu(:,k)
      enddo
      write (file_id,*) ''
!
      write (file_id,*) 'Nu*'
      do k = 1,N_surface_reactions
        write (file_id,writeformat) k,nu_prime(:,k)
      enddo
      write (file_id,*) ''
!
      writeformat = '(  I3)'
      write (writeformat(2:3),'(I2)') N_surface_species
      write (file_id,*) 'J_map, mapping to nchemspec'
      write (file_id,writeformat) jmap
      write (file_id,*) ''
!
    endsubroutine write_outputfile
!***********************************************************************
endmodule Particles_surfspec
