! $Id: particles_surfspec.f90 21950 2014-07-08 08:53:00Z jonas.kruger $
!
!  This module takes care of everything related to reactive particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 8
! MAUX CONTRIBUTION 0
! NADSSPEC CONTRIBUTION 0
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
  real, dimension(10) :: init_surf_gas_frac
  real :: surfplaceholder=0.0
  real, dimension(:), allocatable :: ac
  real, dimension(:), allocatable :: uscale,fscale,constr
  integer, dimension(:), allocatable :: dependent_reactant
!
!  JONAS: implement j_of_inu to communicate with
!  gas phase
!
  integer :: jH2=1, jO2=2, jCO2=3, jCO=4, jCH4=5, jN2=6, jH2O=7
  integer :: jOH=8, jAR=9, jO=10, jH=11, jCH=12, jCH2=13, jHCO=14, jCH3=15
  integer :: idiag_isurf
!
!*********************************************************************!
!               Particle dependent variables below here               !
!*********************************************************************!
!
  real, dimension(:,:), allocatable :: X_infty_reactants
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
!!
!  This is a wrapper routine for particle dependent and particle
!  independent variables
!  JONAS: Back to standalone via mpar_loc=1?
!
          if (lroot) call svn_id( &
          "$Id: particles_surfspec.f90 20849 2014-10-06 18:45:43Z jonas.kruger $")
!
    call register_indep_psurfspec()
    call register_dep_psurfspec()
!!$!
  end subroutine register_particles_surfspec
!************************************************************************
  subroutine register_indep_psurfspec()
!
      integer :: i,k,stat
!
      if (nsurfreacspec/=N_surface_species) then
         print*,'N_surface_species: ', N_surface_species
         call fatal_error('register_particles_ads', &
              'wrong size of storage for surface species allocated')
         else
      endif
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
    call create_dngas()
!
    end subroutine register_indep_psurfspec
!***********************************************************************
    subroutine register_dep_psurfspec()
!
      integer :: stat
!
    allocate(x_infty_reactants(mpar_loc,N_surface_reactants) &
         ,STAT=stat)
    if (stat>0) call fatal_error('register_dep_psurfchem',&
        'Could not allocate memory for x_infty_reactants')
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
subroutine initialize_particles_surf(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  29-sep-14/jonas coded
!  JONAS: needs to be filled with life
!
      real, dimension (mx,my,z,mfarray) :: f
      logical :: lstarting
!
!  
      call keep_compiler_quiet(f)
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
    subroutine rprint_particles_surf(lreset,lwrite)
!
!  Read and register print parameters relevant for
! vparticles near field gas composition
!
!  06-oct-14/jonas: adapted
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      logical :: lwr
      integer :: iname
!
!  Write information to index.pro
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
      if (lwr) write(3,*) 'surf=', isurf
!
      if (lreset) then
         idiag_isurf=0;
      endif
!
      if (lroot .and. ip<14) print*,'rprint_particles_ads: run through parse list'
!
      do iname=1,nname
         call parse_name(iname,cname(iname),cform(iname),'isurf',idiag_isurf)
      enddo
!
    end subroutine rprint_particles_surf
!***********************************************************************
  end module Particles_surfspec
