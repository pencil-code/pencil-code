! $Id: particles_surfspec.f90 21950 2014-07-08 08:53:00Z jonas.kruger $
!
!  This module takes care of everything related to reactive particles.
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

  implicit none

  include 'particles_surfspec.h'

!*********************************************************************!
!               Particle independent variables below here             !
!*********************************************************************!

  character(len=labellen), dimension(ninit) :: init_surf='nothing'
  character(len=10), dimension(:), allocatable :: solid_species
  real, dimension(10) :: init_surf_mol_frac
  real :: surfplaceholder=0.0
  real, dimension(:,:), allocatable :: nu_power
  real, dimension(:), allocatable :: uscale, fscale, constr
  integer :: Ysurf
  logical :: lspecies_transfer=.true.
  logical :: linfinite_diffusion=.true.
  logical :: lboundary_explicit=.true.

  integer :: jH2, jO2, jCO2, jCO, jCH4, jN2, jH2O
  integer :: jOH, jAR, jO, jH, jCH, jCH2, jHCO, jCH3
  integer :: idiag_surf=0

!*********************************************************************!
!               Particle dependent variables below here               !
!*********************************************************************!

  real, dimension(:,:), allocatable :: X_infty_reactants
  real, dimension(:,:), allocatable :: mass_trans_coeff_reactants
  real, dimension(:,:), allocatable :: mass_trans_coeff_species

  namelist /particles_surf_init_pars/ init_surf, init_surf_mol_frac

  namelist /particles_surf_run_pars/ &
      surfplaceholder, &
      lboundary_explicit, &
      linfinite_diffusion, &
      lspecies_transfer

  contains
! ******************************************************************************
!  This is a wrapper routine for particle dependent and particle
!  independent variables
!  JONAS: Back to standalone via mpar_loc=1?

  subroutine register_particles_surfspec()
    if (lroot) call svn_id( &
        "$Id: particles_surfspec.f90 20849 2014-10-06 18:45:43Z jonas.kruger $")

    call register_indep_psurfspec()
    call register_dep_psurfspec()
  endsubroutine register_particles_surfspec
! ******************************************************************************

  subroutine register_indep_psurfspec()
    integer :: i, k, stat, k1, k2
    character(len=10), dimension(40) :: reactants

    if (nsurfreacspec /= N_surface_species) then
      print*,'N_surface_species: ', N_surface_species
      print*,'NSURFREACSPEC :', nsurfreacspec
      call fatal_error('register_particles_surf', &
          'wrong size of storage for surface species allocated')
    endif

!  Increase of npvar according to N_surface_species, which is
!  the concentration of gas phase species at the particle surface

    if (N_surface_species > 1) then
      isurf = npvar+1
      npvar = npvar+N_surface_species
      isurf_end = isurf+N_surface_species-1
    else
      call fatal_error('register_particles_', 'N_surface_species must be > 1')
    endif

    ! Check that the fp and dfp arrays are big enough
    if (npvar > mpvar) then
      if (lroot) write (0,*) 'npvar = ', npvar, ', mpvar = ', mpvar
      call fatal_error('register_ads','npvar > mpvar')
    endif

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
    allocate(uscale(N_surface_reactants),STAT=stat)
    if (stat > 0) call fatal_error('register_indep_psurfchem', &
        'Could not allocate memory for uscale')
    allocate(fscale(N_surface_reactants),STAT=stat)
    if (stat > 0) call fatal_error('register_indep_psurfchem', &
        'Could not allocate memory for fscale')
    allocate(constr(N_surface_reactants),STAT=stat)
    if (stat > 0) call fatal_error('register_indep_psurfchem', 'Could not allocate memory for constr')
    allocate(nu_power(N_surface_species,N_surface_reactions),STAT=stat)
    if (stat > 0) call fatal_error('register_indep_chem', &
        'Could not allocate memory for nu_power')

    call create_ad_sol_lists(solid_species,'sol')
    call get_reactants(reactants)
    call sort_compounds(reactants,solid_species,n_surface_species)

    if (N_surface_species > 1) then
      do i = 1,N_surface_species
        pvarname(isurf+i-1) = solid_species(i)
      enddo
    endif

!  create binding between gas phase and near field gas species
!  chemistry and particles_chemistry module

    call create_jmap()

!  print*, jH2O,JCO2,jH2,jO2,jCO,jCH,jHCO,jCH2,jCH3
!
! Set number of carbon atoms for each surface species

    call get_ac(ac,solid_species,N_surface_species)

    ! Set the stoichiometric matrixes
    print*, 'Number of surface species: ',N_surface_species

    call create_stoc(solid_species,nu,.true., N_surface_species,nu_power)
    call create_stoc(solid_species,nu_prime,.false., N_surface_species)

    ! Define which gas phase reactants the given reaction depends on
    call create_dependency(nu,dependent_reactant, &
        n_surface_reactions,n_surface_reactants)

    ! Find the mole production of the forward reaction
    call create_dngas()

  endsubroutine register_indep_psurfspec
! ******************************************************************************

  subroutine register_dep_psurfspec()
    integer :: stat

    allocate(x_infty_reactants(mpar_loc,N_surface_reactants), STAT = stat)
    if (stat > 0) call fatal_error('register_dep_psurfchem', &
        'Could not allocate memory for x_infty_reactants')
  endsubroutine register_dep_psurfspec
! ******************************************************************************
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-sep-14/jonas coded
!  JONAS: needs to be filled with life

  subroutine initialize_particles_surf(f)
    real, dimension(mx,my,mz,mfarray) :: f

    call keep_compiler_quiet(f)
  endsubroutine initialize_particles_surf
! ******************************************************************************

  subroutine read_particles_surf_init_pars(unit,iostat)
    integer, intent(in) :: unit
    integer, intent(inout), optional :: iostat

    if (present(iostat)) then
      read (unit,NML=particles_surf_init_pars,ERR=99, IOSTAT=iostat)
    else
      read (unit,NML=particles_surf_init_pars,ERR=99)
    endif

99  return
  endsubroutine read_particles_surf_init_pars
! ******************************************************************************

  subroutine write_particles_surf_init_pars(unit)
    integer, intent(in) :: unit

    write (unit,NML=particles_surf_init_pars)
  endsubroutine write_particles_surf_init_pars
! ******************************************************************************

  subroutine read_particles_surf_run_pars(unit,iostat)
    integer, intent(in) :: unit
    integer, intent(inout), optional :: iostat

    if (present(iostat)) then
      read (unit,NML=particles_surf_run_pars,ERR=99, IOSTAT=iostat)
    else
      read (unit,NML=particles_surf_run_pars,ERR=99)
    endif

99  return
  endsubroutine read_particles_surf_run_pars
! ******************************************************************************

  subroutine write_particles_surf_run_pars(unit)
    integer, intent(in) :: unit

    write (unit,NML=particles_surf_run_pars)
  endsubroutine write_particles_surf_run_pars
! ******************************************************************************
!  Initial particle surface fractions
!
!  01-sep-14/jonas: coded

  subroutine init_particles_surf(f,fp)
    real, dimension(mx,my,mz,mfarray) :: f
    real, dimension(mpar_loc,mparray) :: fp
    real :: sum_surf_spec
    integer :: i, j,k

    intent(in) :: f
    intent(out) :: fp

    call keep_compiler_quiet(f)

    fp(:,isurf:isurf_end) = 0.
    do j = 1,ninit

      ! Writing the initial surface species fractions
      select case (init_surf(j))
      case ('nothing')
        if (lroot .and. j == 1) print*, 'init_particles_surf,gas phase: nothing'
      case ('constant')
        if (lroot) print*, 'init_particles_surf: Initial Surface Fractions'

        ! This ensures that we don't have unphysical values as init
        sum_surf_spec = sum(init_surf_mol_frac(1:N_surface_species))
        if (sum_surf_spec > 1) then
          print*, 'Sum of all surface fractions >1, normalizing...'
          init_surf_mol_frac(1:N_surface_species) = &
              init_surf_mol_frac(1:N_surface_species) / sum_surf_spec
        endif

        do k = 1,mpar_loc
          fp(k,isurf:isurf_end) = fp(k,isurf:isurf_end) + &
              init_surf_mol_frac(1:N_surface_species)
        enddo
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

  subroutine dpsurf_dt(f,df,fp,dfp,ineargrid)
    real, dimension(mx,my,mz,mfarray) :: f
    real, dimension(mx,my,mz,mvar) :: df
    real, dimension(mpar_loc,mparray) :: fp
    real, dimension(mpar_loc,mpvar) :: dfp
    integer, dimension(mpar_loc,3) :: ineargrid
    integer :: i

    ! JONAS: equations.tex eq 37
    call keep_compiler_quiet(f)
    call keep_compiler_quiet(df)
!      call keep_compiler_quiet(fp)
    call keep_compiler_quiet(dfp)
    call keep_compiler_quiet(ineargrid)

    if (ldiagnos) then
      if (idiag_surf /= 0) then
        do i=1,isurf_end-isurf
          call sum_par_name(fp(1:mpar_loc,isurf+i-1),idiag_surf+i-1)
        enddo
      endif
    endif
  endsubroutine dpsurf_dt
! ******************************************************************************
!  Evolution of particle surface fraction.
!  (all particles on one pencil)
!
!  23-sep-14/Nils: coded

  subroutine dpsurf_dt_pencil(f,df,fp,dfp,p,ineargrid)
    real, dimension(mx,my,mz,mfarray) :: f
    real, dimension(mx,my,mz,mvar) :: df
    real, dimension(mpar_loc,mparray) :: fp
    real, dimension(mpar_loc,mpvar) :: dfp
    real, dimension(:,:), allocatable :: term,ndot
    real :: porosity
    type (pencil_case) :: p
    integer, dimension(mpar_loc,3) :: ineargrid
    integer :: k, k1, k2, i, ix0, iy0, iz0
    real :: weight, volume_cell, rho1_point
    real :: mean_molar_mass, dmass, A_p, Cg
    integer :: ix1, iy1, iz1,ierr
    integer :: ixx, iyy, izz
    integer :: ixx0, iyy0, izz0
    integer :: ixx1, iyy1, izz1
    integer :: index1, index2
    real, pointer :: true_density_carbon

    intent(in) :: f, ineargrid
    intent(inout) :: dfp, df, fp

    call keep_compiler_quiet(df)
    call keep_compiler_quiet(p)
    call keep_compiler_quiet(ineargrid)

    k1 = k1_imn(imn)
    k2 = k2_imn(imn)

    allocate(term(k1:k2,1:N_surface_reactants))
    allocate(ndot(k1:k2,1:N_surface_species))
!    
    call get_shared_variable('true_density_carbon',true_density_carbon,ierr)
!
    call calc_mass_trans_reactants()
    call get_surface_chemistry(ndot)
!
!  set surface gas species composition
!  (infinite diffusion, set as far field and convert from
!  mass fractions to mole fractions)
!
    do k = k1,k2
      porosity = 1 - (fp(k,imp)/(fp(k,iap)**3*4./3.*pi))/true_density_carbon
      Cg = interp_pp(k)/(R_cgs*interp_TT(k))
      A_p = fp(k1,iap)**2*4.*pi
      mean_molar_mass = (interp_rho(k)*R_cgs*interp_TT(k)/ interp_pp(k))

      if (lboundary_explicit) then
!SOLVE explicitly
      else
        if (linfinite_diffusion) then
          do i = 1,N_surface_reactants
            fp(k,isurf+i-1) = interp_species(k,jmap(i)) / &
                species_constants(jmap(i),imass) * mean_molar_mass
            dfp(k,isurf+i-1) = 0.
          enddo
        else
          do i = 1,N_surface_reactants
            term(k,i) = ndot(k,i)-fp(k,isurf+i-1)*sum(ndot(k,:))+ &
                mass_trans_coeff_reactants(k,i)* &
                (interp_species(k,jmap(i)) / &
                species_constants(jmap(i),imass) * &
                mean_molar_mass-fp(k,isurf+i-1))

            ! the term 3/fp(k,iap) is ratio of the surface of a sphere to its volume
            dfp(k,isurf+i-1) = 3*term(k,i)/(porosity*Cg*fp(k,iap))

            ! JONAS: what was implemented was infinite diffusion, but without created species
          enddo
!SOLVE implicit
        endif
      endif

!  the following block is thoroughly commented in particles_temperature
!  find values for transfer of variables from particle to fluid
      if (lspecies_transfer) then
        ix0 = ineargrid(k,1)
        iy0 = ineargrid(k,2)
        iz0 = ineargrid(k,3)
        call find_interpolation_indeces(ixx0,ixx1,iyy0,iyy1,izz0,izz1, &
            fp,k,ix0,iy0,iz0)
        do izz = izz0,izz1
          do iyy = iyy0,iyy1
            do ixx = ixx0,ixx1
              call find_interpolation_weight(weight,fp,k,ixx,iyy,izz,ix0,iy0,iz0)
              call find_grid_volume(ixx,iyy,izz,volume_cell)
              if ( (iyy /= m).or.(izz /= n).or.(ixx < l1).or.(ixx > l2) ) then
                rho1_point = 1.0 / get_gas_density(f,ixx,iyy,izz)
              else
                rho1_point = p%rho1(ixx-nghost)
              endif

              ! negative dmass means particle is losing mass
              dmass = 0.
              do i = 1,N_surface_species
                index1 = jmap(i)
                dmass = dmass+ndot(k,i)*species_constants(index1,imass)*A_p
              enddo

              do i = 1,N_surface_species
                index1 = jmap(i)
                index2 = ichemspec(index1)
                df(ixx,iyy,izz,index2) = &
                    df(ixx,iyy,izz,index2)+ &
                    (A_p*ndot(k,i)*species_constants(index1,imass)+ &
                    dmass*interp_species(k,index1))* &
                    rho1_point*weight/volume_cell
              enddo
            enddo
          enddo
        enddo
      endif
    enddo

    deallocate(term)
    deallocate(ndot)
  endsubroutine dpsurf_dt_pencil
! ******************************************************************************
!  Read and register print parameters relevant for
!  particles near field gas composition
!
!  06-oct-14/jonas: adapted

  subroutine rprint_particles_surf(lreset,lwrite)
    use Diagnostics

    logical :: lreset
    logical, optional :: lwrite

    logical :: lwr
    integer :: iname

    ! Write information to index.pro
    lwr = .false.
    if (present(lwrite)) lwr = lwrite
    if (lwr) write (3,*) 'Ysurf=', isurf

    if (lreset) then
      idiag_surf = 0
    endif

    if (lroot .and. ip < 14) print*,'rprint_particles_surf: run through parse list'

    do iname = 1,nname
      call parse_name(iname,cname(iname),cform(iname),'Ysurf',idiag_surf)
    enddo
  endsubroutine rprint_particles_surf
! ******************************************************************************
!  07-oct-2014/jonas: coded
!  29-oct-2014/jonas: moved parts of the code into routine, implemented
!                     error when gas phase species is in
!                     particles_chemistry, but not in chemistry
!
!  find the indexes used in chemistry.f90 of species present
!  in the near field of the particle

  subroutine create_jmap()
    use EquationOfState

    integer :: index_glob, index_chem
    logical :: found_species=.false.

    inuH2O = find_species('H2O',solid_species,n_surface_species)
    inuCO2 = find_species('CO2',solid_species,n_surface_species)
    inuH2 = find_species('H2',solid_species,n_surface_species)
    inuO2 = find_species('O2',solid_species,n_surface_species)
    inuCO = find_species('CO',solid_species,n_surface_species)
    inuCH = find_species('CH',solid_species,n_surface_species)
    inuHCO = find_species('HCO',solid_species,n_surface_species)
    inuCH2 = find_species('CH2',solid_species,n_surface_species)
    inuCH3 = find_species('CH3',solid_species,n_surface_species)

    call find_species_index('H2',index_glob,index_chem,found_species)
    if (found_species) then
      jH2 = index_chem
    else
      if (inuH2 > 0) call fatal_error('create_jmap','no H2 found')
    endif

    call find_species_index('O2',index_glob,index_chem,found_species)
    if (found_species) then
      jO2 = index_chem
    else
      if (inuO2 > 0) call fatal_error('create_jmap','no O2 found')
    endif

    call find_species_index('CO2',index_glob,index_chem,found_species)
    if (found_species) then
      jCO2 = index_chem
    else
      if (inuCO2 > 0) call fatal_error('create_jmap','no CO2 found')
    endif

    call find_species_index('CO',index_glob,index_chem,found_species)
    if (found_species) then
      jCO = index_chem
    else
      if (inuCO > 0) call fatal_error('create_jmap','no CO found')
    endif

    call find_species_index('CH4',index_glob,index_chem,found_species)
    if (found_species) then
      jCH4 = index_chem
    else
      if (inuCH4 > 0) call fatal_error('create_jmap','no CH4 found')
    endif

    call find_species_index('H2O',index_glob,index_chem,found_species)
    if (found_species) then
      jH2O = index_chem
    else
      if (inuH2O > 0) call fatal_error('create_jmap','no H2O found')
    endif

    call find_species_index('CH3',index_glob,index_chem,found_species)
    if (found_species) then
      jCH3 = index_chem
    else
      if (inuCH3 > 0) call fatal_error('create_jmap','no CH3 found')
    endif

    call find_species_index('CH',index_glob,index_chem,found_species)
    if (found_species) then
      jCH = index_chem
    else
      if (inuCH > 0) call fatal_error('create_jmap','no CH found')
    endif

    call find_species_index('CH2',index_glob,index_chem,found_species)
    if (found_species) then
      jCH2 = index_chem
    else
      if (inuCH2 > 0) call fatal_error('create_jmap','no CH2 found')
    endif

    call find_species_index('HCO',index_glob,index_chem,found_species)
    if (found_species) then
      jHCO = index_chem
    else
      if (inuHCO > 0)  call fatal_error('create_jmap','no HCO found')
    endif

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

  subroutine calc_mass_trans_reactants()
    integer :: k, i, k1, k2

    k1 = k1_imn(imn)
    k2 = k2_imn(imn)

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

  subroutine allocate_surface_pencils()
    integer :: k1, k2, stat

    k1 = k1_imn(imn)
    k2 = k2_imn(imn)

    allocate(mass_trans_coeff_species(k1:k2,N_species), STAT=stat)
    if (stat > 0) call fatal_error('allocate_variable_pencils', &
        'Could not allocate memory for mass_trans_coeff_species')
    allocate(mass_trans_coeff_reactants(k1:k2,N_surface_reactants),STAT=stat)
    if (stat > 0) call fatal_error('register_indep_psurfchem', &
        'Could not allocate memory for mass_trans_coeff_reactants')
  endsubroutine allocate_surface_pencils
! ******************************************************************************
!
!  Clean up surface pencils
!
!  nov-14/jonas: coded

  subroutine cleanup_surf_pencils()
    deallocate(mass_trans_coeff_species)
    deallocate(mass_trans_coeff_reactants)
  endsubroutine cleanup_surf_pencils
! ******************************************************************************
!
!  Calculate pencils used for surface gas phase fraction evolution
!
!  nov-14/jonas: coded

  subroutine calc_psurf_pencils(f,fp,p,ineargrid)
    real, dimension(mpar_loc,mparray) :: fp
    real, dimension(mx,my,mz,mfarray) :: f
    integer, dimension(mpar_loc,3) :: ineargrid
    type (pencil_case) :: p

    call allocate_surface_pencils()
    call calc_mass_trans_coeff(f,fp,p,ineargrid)
  endsubroutine calc_psurf_pencils
! ******************************************************************************
!  diffusion coefficients of gas species at nearest grid point
!
!  oct-14/Jonas: coded

  subroutine calc_mass_trans_coeff(f,fp,p,ineargrid)
    real, dimension(:,:,:,:) :: f
    real, dimension(:,:) :: fp
    type (pencil_case) :: p
    integer, dimension(:,:) :: ineargrid
    integer :: k, k1, k2,i
    integer :: ix0, iy0, iz0
    integer::spec_glob, spec_chem
    real, dimension(:,:), allocatable :: diff_coeff_species
    real :: Cg

    k1 = k1_imn(imn)
    k2 = k2_imn(imn)

    allocate(diff_coeff_species(k1:k2,N_species))

    do k = k1,k2
      ix0 = ineargrid(k,1)
      iy0 = ineargrid(k,2)
      iz0 = ineargrid(k,3)
      do i = 1,N_surface_species
        diff_coeff_species(k,i) = p%Diff_penc_add(ix0-nghost,jmap(i))
      enddo
    enddo

    do k = k1,k2
      Cg = interp_pp(k)/(R_cgs*interp_TT(k))
      do i = 1,N_surface_species
        mass_trans_coeff_species(k,i) = Cg*diff_coeff_species(k,i)/ fp(k,iap)
      enddo
    enddo

    deallocate(diff_coeff_species)
  endsubroutine calc_mass_trans_coeff
! ******************************************************************************
endmodule Particles_surfspec
