! $Id$
!
!  A module for setting up the f-array and related variables (`register' the
!  entropy, magnetic, etc modules).
!
module Register
!
  use Messages
!
  implicit none
!
  private
!
  public :: register_modules, initialize_modules, rprint_list
  public :: choose_pencils
!
  contains
!***********************************************************************
    subroutine register_modules()
!
!  Call all registration routines, i.e. initialise MPI and register
!  physics modules. Registration implies getting slices of the f-array
!  and setting logicals like lentropy to .true. This routine is called by
!  both, start.x and run.x .
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use General,          only: setup_mm_nn
      use Io,               only: register_io
      use Mpicomm,          only: stop_it
      use Param_Io,         only: get_datadir,get_snapdir
      use Sub
      use Chemistry,        only: register_chemistry
      use Chiral,           only: register_chiral
      use CosmicrayFlux,    only: register_cosmicrayflux
      use Cosmicray,        only: register_cosmicray
      use Density,          only: register_density
      use Dustdensity,      only: register_dustdensity
      use Dustvelocity,     only: register_dustvelocity
      use Entropy,          only: register_entropy
      use EquationOfState,  only: register_eos
      use Forcing,          only: register_forcing
      use Gravity,          only: register_gravity
      use Hydro,            only: register_hydro
      use Hyperresi_strict, only: register_hyperresi_strict
      use Hypervisc_strict, only: register_hypervisc_strict
      use InitialCondition, only: register_initial_condition
      use Interstellar,     only: register_interstellar
      use Lorenz_gauge,     only: register_lorenz_gauge
      use Magnetic,         only: register_magnetic
      use NeutralDensity,   only: register_neutraldensity
      use NeutralVelocity,  only: register_neutralvelocity
      use Polymer,          only: register_polymer
      use Pscalar,          only: register_pscalar
      use Radiation,        only: register_radiation
      use Selfgravity,      only: register_selfgravity
      use Shear,            only: register_shear
      use Shock,            only: register_shock
      use Special,          only: register_special
      use Testfield,        only: register_testfield
      use Testflow,         only: register_testflow
      use TestPerturb,      only: register_testperturb
      use Testscalar,       only: register_testscalar
      use Viscosity,        only: register_viscosity
      use ImplicitPhysics,  only: register_implicit_physics
!
      integer :: ierr
!
!  Overwrite datadir from datadir.in, if that exists.
!
      call get_datadir(datadir)
      call get_snapdir(datadir_snap)
!
!  Initialize index.pro file.
!
      if (lroot) then
        open(3,file=trim(datadir)//'/index.pro',status='replace')
        close(3)
      endif
!
!  Set up the ordering of the pencils.
!
      call setup_mm_nn
!
!  Initialize nvar; is increased by the following routines.
!
      nvar     = 0
      naux     = 0
      naux_com = 0
!
!  Writing files for use with IDL.
!
      if (lroot) then
        if (ldebug) print *, 'Creating ' // trim(datadir) // '/def_var.pro and variables.pro'
        open(15,FILE=trim(datadir)//'/def_var.pro',IOSTAT=ierr)
        if (ierr /= 0) call stop_it("Cannot open "//trim(datadir)// &
            "/def_var.pro for writing -- is "//trim(datadir)//" visible from root node?")
        if (ldebug) print *, 'Creating ' // trim(datadir) // '/variables.pro'
        open(4,FILE=trim(datadir)//'/variables.pro',IOSTAT=ierr)
        if (ierr /= 0) call stop_it("Cannot open "//trim(datadir)// &
            "/variables.pro for writing -- is "//trim(datadir)//" visible from root node?")
        write(4,*) 'close,1'
        write(4,*) "openr,1, datadir+'/'+varfile, /F77"
        write(4,*) 'readu,1 $'
      endif
!
!  Initialize file for writing constants to be read by IDL.
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro')
        write (1,*) '; This file contain pc constants of interest to IDL'
        close (1)
      endif
!
!  Register variables in the f-array.
!
      call register_io
      call register_initial_condition
      call register_eos
      call register_shock
      call register_viscosity             !(should go under hydro)
      call register_hydro
      call register_gravity
      call register_selfgravity           !(why not connected with gravity)
      call register_density
      call register_forcing
      call register_entropy
      call register_magnetic
      call register_lorenz_gauge          !(should go under magnetic)
      call register_polymer
      call register_testscalar
      call register_testfield
      call register_testflow
      call register_radiation
      call register_pscalar
      call register_chiral
      call register_chemistry
      call register_dustvelocity
      call register_dustdensity
      call register_neutralvelocity
      call register_neutraldensity
      call register_cosmicray
      call register_cosmicrayflux
      call register_interstellar
      call register_shear
      call register_hypervisc_strict
      call register_hyperresi_strict
      call register_implicit_physics
      call register_special
!
!  Writing files for use with IDL.
!
      if (lroot) then
        do aux_count=1,maux
          write(4,'(a10)') aux_var(aux_count)
        enddo
        close(4)
        close(15)
      endif
!
      if (nvar /= mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_modules: nvar /= mvar. '// &
           'Check your MVAR and/or MAUX CONTRIBUTION in cparam.local')
      endif
!
!  Initialize headt for root processor only.
!
      if (lroot) headt=.true.
!
    endsubroutine register_modules
!***********************************************************************
    subroutine initialize_modules(f,lstarting)
!
!  Call initialization routines, i.e. initialize physics and technical
!  modules. This implies some preparation of auxiliary quantities, often
!  based on input parameters. This routine is called by run.x (but not by
!  start.x) initially and each time the run parameters have been reread.
!
!  6-nov-01/wolf: coded
! 23-feb-03/axel: added physical constants conversion
!  7-oct-03/david: initialize_gravity before density, etc (its needed there)
! 11-sep-04/axel: began adding spherical coordinates
!
      use Cdata
      use Param_IO
      use Mpicomm,          only: mpireduce_sum,mpibcast_real,&
                                  mpisend_real,mpirecv_real
      use Sub,              only: remove_zprof
      use BorderProfiles,   only: initialize_border_profiles
      use Chemistry,        only: initialize_chemistry
      use Chiral,           only: initialize_chiral
      use CosmicrayFlux,    only: initialize_cosmicrayflux
      use Cosmicray,        only: initialize_cosmicray
      use Density,          only: initialize_density
      use Deriv,            only: initialize_deriv
      use Diagnostics,      only: initialize_prints
      use Dustdensity,      only: initialize_dustdensity
      use Dustvelocity,     only: initialize_dustvelocity
      use Entropy,          only: initialize_entropy
      use EquationOfState,  only: initialize_eos, units_eos
      use Forcing,          only: initialize_forcing
      use Gravity,          only: initialize_gravity
      use Hydro,            only: initialize_hydro
      use InitialCondition, only: initialize_initial_condition
      use Interstellar,     only: initialize_interstellar
      use Magnetic,         only: initialize_magnetic
      use Lorenz_gauge,     only: initialize_lorenz_gauge
      use Polymer,          only: initialize_polymer
      use NeutralDensity,   only: initialize_neutraldensity
      use NeutralVelocity,  only: initialize_neutralvelocity
      use Poisson,          only: initialize_poisson
      use Pscalar,          only: initialize_pscalar
      use Radiation,        only: initialize_radiation
      use Selfgravity,      only: initialize_selfgravity
      use Shear,            only: initialize_shear
      use Shock,            only: initialize_shock
      use Solid_Cells,      only: initialize_solid_cells
      use Special,          only: initialize_special
      use Testfield,        only: initialize_testfield
      use Testflow,         only: initialize_testflow
      use TestPerturb,      only: initialize_testperturb
      use Testscalar,       only: initialize_testscalar
      use Timeavg,          only: initialize_timeavg
      use Viscosity,        only: initialize_viscosity
      use ImplicitPhysics,  only: initialize_implicit_physics
      use Grid,             only: initialize_grid
!
      real, dimension(mx,my,mz,mfarray) :: f
      logical :: lstarting
      integer :: ivar
!
!  Defaults for some logicals; will later be set to true if needed.
!
      lpenc_requested(:) = .false.
!
!  Evaluate physical units. Used currently only in eos, but later also
!  in the interstellar and radiation modules, for example.
!
      call units_general()
      call units_eos()
!
!  Calculated derived units.
!
      unit_mass=unit_density*unit_length**3
      unit_energy=unit_mass*unit_velocity**2
      unit_time=unit_length/unit_velocity
      unit_flux=unit_energy/(unit_length**2*unit_time)
!
!  Convert physical constants to code units.
!
      if (unit_system=='cgs') then
        if (lroot.and.leos_ionization.and.ip<14) print*,'initialize_modules: ' &
          //'unit_velocity, unit_density, etc, are in cgs'
        hbar=hbar_cgs/(unit_energy*unit_time)
        mu0=mu0_cgs*unit_density*(unit_velocity/unit_magnetic)**2
        if (unit_temperature/=impossible) then
          sigmaSB=sigmaSB_cgs/(unit_flux/unit_temperature**4)
          k_B=k_B_cgs/(unit_energy/unit_temperature)
        else
          sigmaSB=impossible
          k_B=impossible
        endif
        m_u=m_u_cgs/unit_mass
        m_p=m_p_cgs/unit_mass
        m_e=m_e_cgs/unit_mass
        eV=eV_cgs/unit_energy
        sigmaH_=sigmaH_cgs/unit_length**2
        kappa_es=kappa_es_cgs/(unit_length**2/unit_mass)
        c_light=c_light_cgs/unit_velocity
        G_Newton=G_Newton_cgs*unit_length**2*unit_density/unit_velocity**2
      elseif (unit_system=='SI') then
        if (lroot.and.leos_ionization) print*,&
            'initialize_modules: unit_velocity, unit_density, etc, are in SI'
        hbar=hbar_cgs*1e-7/(unit_energy*unit_time)
        mu0=1e-7*mu0_cgs*unit_density*(unit_velocity/unit_magnetic)**2
        if (unit_temperature/=impossible) then
          sigmaSB=sigmaSB_cgs*1e-3/(unit_flux/unit_temperature**4)
          k_B=1e-7*k_B_cgs/(unit_energy/unit_temperature)
        else
          sigmaSB=impossible
          k_B=impossible
        endif
        m_u=m_u_cgs*1e-3/unit_mass
        m_p=m_p_cgs*1e-3/unit_mass
        m_e=m_e_cgs*1e-3/unit_mass
        eV=eV_cgs*1e-7/unit_energy
        sigmaH_=sigmaH_cgs*1e-4/unit_length**2
        kappa_es=kappa_es_cgs*1e-1/(unit_length**2/unit_mass)
        c_light=c_light_cgs*1e-2/unit_velocity
        G_Newton=G_Newton_cgs*1e-3*unit_length**2*unit_density/unit_velocity**2
      endif
!
!  Calculate additional constants (now all in code units).
!
      m_H=m_p+m_e
      m_He=3.97153*m_H
      chiH=13.6*eV
      chiH_=0.754*eV
!
!  Print parameters in code units, but only when used.
!
      if (lroot.and.ip<14) then
         if (leos_ionization.or.lradiation.or.lradiation_ray.or.linterstellar) then
            write(*,'(a,1p,4e14.6)') ' register: k_B,m_p,m_e,eV=',k_B,m_p,m_e,eV
         endif
      endif
!
!  initialize time integrals
!  (leads to linker problem)
!
!---  call initialize_time_integrals(f)
!
! Store the value of impossible for use in IDL
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,*) 'impossible=',impossible
        close (1)
      endif
!
!  print summary of variable names
!
      if (lroot) then
        open(3,file=trim(datadir)//'/varname.dat',status='replace')
        do ivar=1,nvar
          write(3,"(i4,2x,a)") ivar,varname(ivar)
        enddo
        close(3)
      endif
!
!  Coordinate-related issues, initialize specific grid variables
!
      call initialize_grid()
!
!  timestep: distinguish two cases,
!  (a) dt explicitly given in run.in -> ldt=.false.
!  (b) dt not given in run.in        -> ldt=.true.  -> calculate dt dynamically
!  Note that ldt will not change unless you RELOAD parameters.
!
!  Note that this should not be moved to timestep.f90 as
!  run_hooks_timestep(), because maybe not, because initialize_modules
!  can also be run from start.f90, which has no knowledge of timestep.f90
!
!  The calculation of ldt needs to be done calculate dt dynamically,
!  because the time step can be changed after a reload.
!
      ldt = (dt==0.)
      if (lroot .and. ip<14) then
        if (ldt) then
          print*,'timestep based on CFL cond; cdt=',cdt
        else
          print*, 'absolute timestep dt=', dt
        endif
      endif
!
!  Run rest of initialization of individual modules.
!
      call initialize_deriv()
      call initialize_prints()
      call initialize_timeavg(f)
      call initialize_initial_condition(f)
      call initialize_eos()
      call initialize_gravity(f,lstarting)
      call initialize_selfgravity(f)
      call initialize_poisson()
      call initialize_density(f,lstarting)
      call initialize_hydro(f,lstarting)
      call initialize_forcing(lstarting)
      call initialize_entropy(f,lstarting)
      call initialize_magnetic(f,lstarting)
      call initialize_lorenz_gauge(f)
      call initialize_polymer(f,lstarting)
      call initialize_testscalar(f)
      call initialize_testfield(f,lstarting)
      call initialize_testflow(f)
      call initialize_radiation()
      call initialize_pscalar(f)
      call initialize_chiral(f)
      call initialize_chemistry(f)
      call initialize_dustvelocity(f)
      call initialize_dustdensity(f)
      call initialize_neutraldensity()
      call initialize_neutralvelocity()
      call initialize_cosmicray(f)
      call initialize_cosmicrayflux(f)
      call initialize_interstellar(f,lstarting)
      call initialize_shear()
      call initialize_testperturb()
      call initialize_shock(f,lstarting)
      call initialize_viscosity(lstarting)
      call initialize_special(f,lstarting)
      call initialize_border_profiles()
      call initialize_solid_cells(f)
      call initialize_implicit_physics(f)
!
    endsubroutine initialize_modules
!***********************************************************************
    subroutine units_general()
!
!  This routine calculates things related to units and must be called
!  before the rest of the units are being calculated.
!
!  12-jul-06/axel: adapted from units_eos
!
      use Cdata, only: unit_system,G_Newton,c_light,hbar,lroot, &
        unit_length,unit_velocity,unit_density,unit_magnetic
      use Cparam, only: G_Newton_cgs,c_light_cgs,hbar_cgs,impossible
!
!  Unless G_Newton,c_light,hbar are all set,
!  unit_velocity,unit_density,unit_length will be calculated
!
      if (G_Newton == impossible .or. &
           c_light == impossible .or. &
              hbar == impossible) then
        if (unit_velocity == impossible) unit_velocity=1.
        if (unit_density == impossible) unit_density=1.
        if (unit_length == impossible) unit_length=1.
      else
        if (unit_system == 'cgs') then
           unit_velocity=c_light_cgs/c_light
           unit_density=unit_velocity**5/((G_Newton_cgs/G_Newton)**2 &
                *(hbar_cgs/hbar))
           unit_length=sqrt((G_Newton_cgs/G_Newton) &
                *(hbar_cgs/hbar)/unit_velocity**3)
        elseif (unit_system == 'SI') then
           unit_velocity=c_light_cgs*1e-2/c_light
           unit_density=unit_velocity**5/((G_Newton_cgs*1e-3/G_Newton)**2 &
                *(hbar_cgs*1e-7/hbar))
           unit_length=sqrt((G_Newton_cgs*1e-3/G_Newton) &
                *(hbar_cgs*1e-7/hbar)/unit_velocity**3)
        endif
      endif
!
!  Set unit_magnetic=3.5449077018110318=sqrt(4*pi), unless it was set already.
!  Note that unit_magnetic determines the value of mu_0 in the rest of the code.
!
      if (unit_magnetic == impossible) unit_magnetic=3.5449077018110318
!
!  Check that everything is OK.
!
      if (lroot) print*,'units_general: unit_velocity=',unit_velocity
      if (lroot) print*,'units_general: unit_density=',unit_density
      if (lroot) print*,'units_general: unit_length=',unit_length
      if (lroot) print*,'units_general: unit_magnetic=',unit_magnetic
!
    endsubroutine units_general
!***********************************************************************
    subroutine choose_pencils()
!
!  Find out which pencils are needed for all time-steps and also for
!  diagnostics only. Also takes care of interdependent pencils.
!
!  20-nov-04/anders: coded
!
      use Cdata
!
      integer :: i
!
      if (lroot) print*, 'choose_pencils: finding out which pencils '// &
          'are needed for the pencil case'
!
!  Must set all pencil arrays to false in case of reload.
!
      lpenc_requested=.false.
      lpenc_diagnos=.false.
      lpenc_diagnos2d=.false.
!
!  Find out which pencils are needed for the pencil case.
!
      call pencil_criteria()
!
!  Set interdependent pencils.
!
      do i=1,3
        call pencil_interdep(lpenc_requested)
        call pencil_interdep(lpenc_diagnos)
        call pencil_interdep(lpenc_diagnos2d)
        call pencil_interdep(lpenc_video)
      enddo
!
!  Swap logical content of pencil number ipencil_swap (for testing).
!
      if (ipencil_swap/=0) then
        if (lpencil_requested_swap) then
          lpenc_requested(ipencil_swap) = (.not. lpenc_requested(ipencil_swap))
          print*, 'choose_pencils: Swapped requested pencil number ', &
              ipencil_swap, ' to ', lpenc_requested(ipencil_swap)
        endif
        if (lpencil_diagnos_swap) then
          lpenc_diagnos(ipencil_swap) = (.not. lpenc_diagnos(ipencil_swap))
          print*, 'choose_pencils: Swapped diagnostic pencil number ', &
              ipencil_swap, ' to ', lpenc_diagnos(ipencil_swap)
        endif
      endif
!
    endsubroutine choose_pencils
!***********************************************************************
    subroutine pencil_criteria()
!
!  Find out which pencils are needed for all the modules. In each call
!  the called module will inform about the pencils that it needs locally.
!  Interdependency among pencils is not solved here.
!
!  Note: No pencils can exist for the forcing module, because it is
!  used outside the pencil mn loop, so rho and 1/rho needs to be
!  calculated separately.
!
!  20-11-04/anders: coded
!
      use Cdata
      use Grid,            only: pencil_criteria_grid
      use BorderProfiles,  only: pencil_criteria_borderprofiles
      use EquationOfState, only: pencil_criteria_eos
      use Hydro,           only: pencil_criteria_hydro
      use Density,         only: pencil_criteria_density
      use Forcing,         only: pencil_criteria_forcing
      use Shock,           only: pencil_criteria_shock
      use Viscosity,       only: pencil_criteria_viscosity
      use Entropy,         only: pencil_criteria_entropy
      use Gravity,         only: pencil_criteria_gravity
      use Selfgravity,     only: pencil_criteria_selfgravity
      use Pscalar,         only: pencil_criteria_pscalar
      use Chemistry,       only: pencil_criteria_chemistry
      use Dustvelocity,    only: pencil_criteria_dustvelocity
      use Dustdensity,     only: pencil_criteria_dustdensity
      use NeutralVelocity, only: pencil_criteria_neutralvelocity
      use NeutralDensity,  only: pencil_criteria_neutraldensity
      use Magnetic,        only: pencil_criteria_magnetic
      use Lorenz_gauge,    only: pencil_criteria_lorenz_gauge
      use Polymer,         only: pencil_criteria_polymer
      use Testscalar,      only: pencil_criteria_testscalar
      use Testfield,       only: pencil_criteria_testfield
      use Testflow,        only: pencil_criteria_testflow
      use Cosmicray,       only: pencil_criteria_cosmicray
      use Cosmicrayflux,   only: pencil_criteria_cosmicrayflux
      use Chiral,          only: pencil_criteria_chiral
      use Radiation,       only: pencil_criteria_radiation
      use Interstellar,    only: pencil_criteria_interstellar
      use Shear,           only: pencil_criteria_shear
      use Special,         only: pencil_criteria_special
      use Particles_main,  only: particles_pencil_criteria
      use Solid_cells,     only: pencil_criteria_solid_cells
!
      call pencil_criteria_grid()
      call pencil_criteria_borderprofiles()
      call pencil_criteria_density()
      call pencil_criteria_forcing()
      call pencil_criteria_eos()
      call pencil_criteria_hydro()
      call pencil_criteria_shock()
      call pencil_criteria_viscosity()
      call pencil_criteria_entropy()
      call pencil_criteria_gravity()
      call pencil_criteria_selfgravity()
      call pencil_criteria_pscalar()
      call pencil_criteria_interstellar()
      call pencil_criteria_chemistry()
      call pencil_criteria_dustvelocity()
      call pencil_criteria_dustdensity()
      call pencil_criteria_neutralvelocity()
      call pencil_criteria_neutraldensity()
      call pencil_criteria_magnetic()
      call pencil_criteria_lorenz_gauge()
      call pencil_criteria_polymer()
      call pencil_criteria_testscalar()
      call pencil_criteria_testfield()
      call pencil_criteria_testflow()
      call pencil_criteria_cosmicray()
      call pencil_criteria_cosmicrayflux()
      call pencil_criteria_chiral()
      call pencil_criteria_radiation()
      call pencil_criteria_shear()
      call pencil_criteria_special()
      call pencil_criteria_solid_cells()
      if (lparticles) call particles_pencil_criteria()
!
    endsubroutine pencil_criteria
!***********************************************************************
    subroutine pencil_interdep(lpencil_in)
!
!  Find out about interdependent pencils. Each module knows what its own
!  pencils depend on. The dependency only needs to be specified one level
!  up, since this subroutine is called several times (currently three).
!
!  Note: No pencils can exist for the forcing module, because it is
!  used outside the pencil mn loop, so rho and 1/rho needs to be
!  calculated separately.
!
!  20-11-04/anders: coded
!
      use Cdata
      use EquationOfState, only: pencil_interdep_eos
      use Hydro, only: pencil_interdep_hydro
      use Density, only: pencil_interdep_density
      use Forcing, only: pencil_interdep_forcing
      use Shock, only: pencil_interdep_shock
      use Viscosity, only: pencil_interdep_viscosity
      use Entropy, only: pencil_interdep_entropy
      use Gravity, only: pencil_interdep_gravity
      use Selfgravity, only: pencil_interdep_selfgravity
      use Magnetic, only: pencil_interdep_magnetic
      use Lorenz_gauge, only: pencil_interdep_lorenz_gauge
      use Polymer, only: pencil_interdep_polymer
      use Testscalar, only: pencil_interdep_testscalar
      use Testfield, only: pencil_interdep_testfield
      use Testflow, only: pencil_interdep_testflow
      use Pscalar, only: pencil_interdep_pscalar
      use Chemistry, only: pencil_interdep_chemistry
      use Dustvelocity, only: pencil_interdep_dustvelocity
      use Dustdensity, only: pencil_interdep_dustdensity
      use NeutralVelocity, only: pencil_interdep_neutralvelocity
      use NeutralDensity, only: pencil_interdep_neutraldensity
      use Cosmicray, only: pencil_interdep_cosmicray
      use Cosmicrayflux, only: pencil_interdep_cosmicrayflux
      use Chiral, only: pencil_interdep_chiral
      use Radiation, only: pencil_interdep_radiation
      use Shear, only: pencil_interdep_shear
      use Special, only: pencil_interdep_special
      use Grid, only: pencil_interdep_grid
      use Particles_main, only: particles_pencil_interdep
!
      logical, dimension (npencils) :: lpencil_in
!
      call pencil_interdep_grid(lpencil_in)
      call pencil_interdep_density(lpencil_in)
      call pencil_interdep_forcing(lpencil_in)
      call pencil_interdep_eos(lpencil_in)
      call pencil_interdep_hydro(lpencil_in)
      call pencil_interdep_shock(lpencil_in)
      call pencil_interdep_viscosity(lpencil_in)
      call pencil_interdep_entropy(lpencil_in)
      call pencil_interdep_gravity(lpencil_in)
      call pencil_interdep_selfgravity(lpencil_in)
      call pencil_interdep_chemistry(lpencil_in)
      call pencil_interdep_dustvelocity(lpencil_in)
      call pencil_interdep_dustdensity(lpencil_in)
      call pencil_interdep_neutralvelocity(lpencil_in)
      call pencil_interdep_neutraldensity(lpencil_in)
      call pencil_interdep_pscalar(lpencil_in)
      call pencil_interdep_magnetic(lpencil_in)
      call pencil_interdep_lorenz_gauge(lpencil_in)
      call pencil_interdep_polymer(lpencil_in)
      call pencil_interdep_testscalar(lpencil_in)
      call pencil_interdep_testfield(lpencil_in)
      call pencil_interdep_testflow(lpencil_in)
      call pencil_interdep_cosmicray(lpencil_in)
      call pencil_interdep_cosmicrayflux(lpencil_in)
      call pencil_interdep_chiral(lpencil_in)
      call pencil_interdep_radiation(lpencil_in)
      call pencil_interdep_shear(lpencil_in)
      call pencil_interdep_special(lpencil_in)
      if (lparticles) call particles_pencil_interdep(lpencil_in)
!
    endsubroutine pencil_interdep
!***********************************************************************
    logical function read_name_format(in_file,cnamel,nnamel)
!
!  Unifies reading of *.in files which contain requests for diagnostic
!  output in the form <name>(<format>); returns number of items read !
!  properly from file 'in_file' (comment lines excluded) in
!  'nnamel'; returns items in cnamel, which is allocated, if necessary, with
!  the length <number of lines in 'in_file' + initial value of
!  'nnamel'>; further allocations done in subroutine argument 'allocator' which
!  takes same length as its parameter;
!
!  Return value is nnamel>0.
!
!   11-jan-11/MR: coded
!   23-jan-11/MR: pointer based handling of cname-like arrays instead of allocatable
!                 dummy parameter (the latter standardized only since FORTRAN 2000)
!   24-jan-11/MR: removed allocation
!
      use Mpicomm, only: parallel_open, parallel_close
!
      character (len=*) ,               intent(in)    :: in_file
      character (len=30), dimension(*), intent(out)   :: cnamel
      integer,                          intent(inout) :: nnamel
!
      character (len=30) :: cname_tmp
      integer            :: iname, mname
      integer, parameter :: unit=1
!
      call parallel_open(unit,FILE=trim(in_file))
!
!  Read names and formats.
!
      mname  = nnamel
      nnamel = 0
      do iname=1,mname
        read(unit,*,end=99) cname_tmp
        if ((cname_tmp(1:1)/='!') .and. (cname_tmp(1:1)/='#')) then
          nnamel=nnamel+1
          cnamel(nnamel)=cname_tmp
        endif
      enddo
!
99    call parallel_close(unit)
!
      read_name_format = nnamel>0
!
    endfunction read_name_format
!***********************************************************************
    subroutine rprint_list(lreset)
!
!  Read variables to print and to calculate averages of from control files.
!
!  3-may-01/axel: coded
!  11-jan-11/MR: introduced read_name_format calls for each of the lists
!                for homogeneity
!
!  All numbers like nname etc. need to be initialized to zero in cdata!
!
      use Cdata
      use Sub,             only: numeric_precision
      use Diagnostics
      use Hydro,           only: rprint_hydro
      use Density,         only: rprint_density
      use Forcing,         only: rprint_forcing
      use Entropy,         only: rprint_entropy
      use Magnetic,        only: rprint_magnetic
      use Lorenz_gauge,    only: rprint_lorenz_gauge
      use Polymer,         only: rprint_polymer
      use Testscalar,      only: rprint_testscalar
      use Testfield,       only: rprint_testfield
      use Testflow,        only: rprint_testflow
      use Radiation,       only: rprint_radiation
      use EquationOfState, only: rprint_eos
      use Pscalar,         only: rprint_pscalar
      use Chiral,          only: rprint_chiral
      use Interstellar,    only: rprint_interstellar
      use Chemistry,       only: rprint_chemistry
      use Dustvelocity,    only: rprint_dustvelocity
      use Dustdensity,     only: rprint_dustdensity
      use NeutralVelocity, only: rprint_neutralvelocity
      use NeutralDensity,  only: rprint_neutraldensity
      use Cosmicray,       only: rprint_cosmicray
      use CosmicrayFlux,   only: rprint_cosmicrayflux
      use Gravity,         only: rprint_gravity
      use Selfgravity,     only: rprint_selfgravity
      use Special,         only: rprint_special
      use Shock,           only: rprint_shock
      use solid_cells,     only: rprint_solid_cells
      use Viscosity,       only: rprint_viscosity
      use Shear,           only: rprint_shear
      use TestPerturb,     only: rprint_testperturb
      use Mpicomm,         only: parallel_file_exists, parallel_count_lines, &
                                 parallel_open, parallel_close
!
      integer :: unit=1
      integer :: ierr,iadd,ios
      logical :: lreset, ldummy
      character (len=30) :: cname_tmp
!
      character (LEN=15)           :: print_in_file
      character (LEN=*), parameter :: video_in_file    = 'video.in'
      character (LEN=*), parameter :: sound_in_file    = 'sound.in'
      character (LEN=*), parameter :: xyaver_in_file   = 'xyaver.in'
      character (LEN=*), parameter :: xzaver_in_file   = 'xzaver.in'
      character (LEN=*), parameter :: yzaver_in_file   = 'yzaver.in'
      character (LEN=*), parameter :: phizaver_in_file = 'phizaver.in'
      character (LEN=*), parameter :: yaver_in_file    = 'yaver.in'
      character (LEN=*), parameter :: zaver_in_file    = 'zaver.in'
      character (LEN=*), parameter :: phiaver_in_file  = 'phiaver.in'
!
!  In all lists:  recognize "!" and "#" as comments.
!
!  Read print.in.double if applicable, else print.in.
!  Read in the list of variables to be printed.
!
      print_in_file = 'print.in'
      if (numeric_precision() == 'D') then
        if (parallel_file_exists(trim(print_in_file)//'.double')) &
            print_in_file = trim(print_in_file)//'.double'
      endif
!
      if (lroot) print*, 'Reading print formats from '//trim(print_in_file)
!
      nname = parallel_count_lines(print_in_file)
!
      if (nname>0) then
        call allocate_fnames(100)
        ldummy = read_name_format(print_in_file,cname,nname)
      endif
      if ( nname==0 ) &
          call fatal_error('rprint_list','You must have a "'// &
          trim(print_in_file)// &
          '" file in the run directory with valid print requests!')
!
      if (lroot .and. (ip<14)) print*, 'rprint_list: nname=', nname
!
!  Read in the list of variables for video slices.
!
      if ( dvid/=0.0 ) then

        nnamev = parallel_count_lines(video_in_file)
!
        if (nnamev>0) then
          call allocate_vnames(nnamev)
          if ( .not.read_name_format(video_in_file,cnamev,nnamev) ) &
            dvid=0.0
        endif
      endif
      if (lroot .and. (ip<14)) &
          print*, 'rprint_list: ix,iy,iz,iz2=', ix,iy,iz,iz2
      if (lroot .and. (ip<14)) print*, 'rprint_list: nnamev=', nnamev
!
!  Read in the list of variables for "sound".
!
!  In allocate_sound the relevant arrays are allocated and the list of
!  coordinates sound_coords_list is read in.
!  nsound_location and lwrite_sound are set there, too.
!
      if ( dimensionality>0 .and. dsound/=0.0 ) then
        
        nname_sound = parallel_count_lines(sound_in_file)
!       
        if (nname_sound>0) then
          call allocate_sound(nname_sound)

          if ( read_name_format(sound_in_file,cname_sound,nname_sound) &
               .and. lwrite_sound ) then
!
!  Read the last sound output time from a soundfile, will be set to
!  starttime otherwise
!
!            tsound=rnan
            tsound=-1.0
            open(1,file=trim(directory)//'/sound.dat',position='append', &
                 status='old',iostat=ios)
            if (ios==0) then
              backspace(1)
              read(1,*) tsound
            endif
            close(1)
!
          endif
!
        else
          nname_sound=0
        endif
      endif
      if (lroot .and. (ip<14)) &
          print*, 'sound_print_list: nname_sound=', nname_sound
!
!  Read in the list of variables for xy-averages.
!
      nnamez = parallel_count_lines(xyaver_in_file)
!
      if (nnamez>0) then
        call allocate_xyaverages(nnamez)
        ldummy = read_name_format(xyaver_in_file,cnamez,nnamez)
      endif
      if (lroot .and. (ip<14)) print*, 'rprint_list: nnamez=', nnamez
!
!  Read in the list of variables for xz-averages.
!
      nnamey = parallel_count_lines(xzaver_in_file)
!
      if (nnamey>0) then
        call allocate_xzaverages(nnamey)
        ldummy = read_name_format(xyaver_in_file,cnamey,nnamey)
      endif
      if (lroot .and. (ip<14)) print*, 'rprint_list: nnamey=', nnamey
!
!  Read in the list of variables for yz-averages.
!
      nnamex = parallel_count_lines(yzaver_in_file)
!
      if (nnamex>0) then
        call allocate_yzaverages(nnamex)
        ldummy = read_name_format(yzaver_in_file,cnamex,nnamex)
      endif

      if (lroot .and. (ip<14)) print*, 'rprint_list: nnamex=', nnamex
!
!  Read in the list of variables for phi-z-averages.
!
      nnamer = parallel_count_lines(phizaver_in_file)
!
      if (nnamer>0) then
        call allocate_phizaverages(nnamer)
        ldummy = read_name_format(phizaver_in_file,cnamer,nnamer)
      endif
      if (lroot .and. (ip<14)) print*, 'rprint_list: nnamer=', nnamer
!
!  2-D averages:
!
!  Read in the list of variables for y-averages.
!
      nnamexz = parallel_count_lines(yaver_in_file)
!
      if (nnamexz>0) then
        call allocate_yaverages(nnamexz)
        lwrite_yaverages = read_name_format(yaver_in_file,cnamexz,nnamexz)
      endif

      if (lroot .and. (ip<14)) print*, 'rprint_list: nnamexz=', nnamexz
!
!  Read in the list of variables for z-averages.
!
      nnamexy = parallel_count_lines(zaver_in_file)
!     
      if (nnamexy>0) then                  
        call allocate_zaverages(nnamexy)
        lwrite_zaverages = read_name_format(zaver_in_file,cnamexy,nnamexy)
      endif

      if (lroot .and. (ip<14)) print*, 'rprint_list: nnamexy=', nnamexy
!
!  Read in the list of variables for phi-averages.
!
      if (parallel_file_exists(phiaver_in_file)) then
        call parallel_open(unit,file=phiaver_in_file)
        iadd=0
!
        do
          read(unit,*,iostat=ierr) cname_tmp
          if (ierr==0) then
            nnamerz=nnamerz+1
            if ( cname_tmp=='uumphi' .or. cname_tmp=='uusphmphi' .or.  &
                 cname_tmp=='bbmphi' .or. cname_tmp=='bbsphmphi' .or.  &
                 cname_tmp=='uxbmphi'.or. cname_tmp=='jxbmphi'       ) &
              iadd=iadd+2
          else
            exit
          endif
        enddo
        call parallel_close(unit)
!
        if (nnamerz>0) then                  
          call allocate_phiaverages(nnamerz+iadd)
          lwrite_phiaverages = read_name_format(phiaver_in_file,cnamerz,nnamerz)
        endif
      endif
      if (lroot .and. (ip<14)) print*, 'rprint_list: nnamerz=', nnamerz
!
!  Set logical for 2-D averages.
!
      lwrite_2daverages= &
          lwrite_yaverages.or.lwrite_zaverages.or.lwrite_phiaverages
!
!  Check which variables are set.
!  For the convenience of idl users, the indices of variables in
!  the f-array and the time_series.dat files are written to data/index.pro.
!
      if (lroot) open(3,file=trim(datadir)//'/index.pro',position='append')
!
      call rprint_general         (lreset,LWRITE=lroot)
      call rprint_hydro           (lreset,LWRITE=lroot)
      call rprint_density         (lreset,LWRITE=lroot)
      call rprint_forcing         (lreset,LWRITE=lroot)
      call rprint_entropy         (lreset,LWRITE=lroot)
      call rprint_magnetic        (lreset,LWRITE=lroot)
      call rprint_lorenz_gauge    (lreset,LWRITE=lroot)
      call rprint_polymer         (lreset,LWRITE=lroot)
      call rprint_testscalar      (lreset,LWRITE=lroot)
      call rprint_testfield       (lreset,LWRITE=lroot)
      call rprint_testflow        (lreset,LWRITE=lroot)
      call rprint_radiation       (lreset,LWRITE=lroot)
      call rprint_eos             (lreset,LWRITE=lroot)
      call rprint_pscalar         (lreset,LWRITE=lroot)
      call rprint_chiral          (lreset,LWRITE=lroot)
      call rprint_interstellar    (lreset,LWRITE=lroot)
      call rprint_chemistry       (lreset,LWRITE=lroot)
      call rprint_dustvelocity    (lreset,LWRITE=lroot)
      call rprint_dustdensity     (lreset,LWRITE=lroot)
      call rprint_neutralvelocity (lreset,LWRITE=lroot)
      call rprint_neutraldensity  (lreset,LWRITE=lroot)
      call rprint_cosmicray       (lreset,LWRITE=lroot)
      call rprint_cosmicrayflux   (lreset,LWRITE=lroot)
      call rprint_gravity         (lreset,LWRITE=lroot)
      call rprint_selfgravity     (lreset,LWRITE=lroot)
      call rprint_special         (lreset,LWRITE=lroot)
      call rprint_shock           (lreset,LWRITE=lroot)
      call rprint_solid_cells     (lreset,LWRITE=lroot)
      call rprint_viscosity       (lreset,LWRITE=lroot)
      call rprint_shear           (lreset,LWRITE=lroot)
      call rprint_testperturb     (lreset,LWRITE=lroot)
      if (lroot) close(3)
!
    endsubroutine rprint_list
!***********************************************************************
    subroutine rprint_general(lreset,lwrite)
!
!  Reads and registers *general* print parameters.
!
!   8-jun-02/axel: adapted from hydro
!
      use Cdata
      use Diagnostics
!
      integer :: iname,irz
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of RELOAD.
!  (general variables that are defined in Cdata)
!
      if (lreset) then
        idiag_t=0; idiag_it=0; idiag_dt=0; idiag_walltime=0
        idiag_timeperstep=0
        idiag_rcylmphi=0; idiag_phimphi=0; idiag_zmphi=0; idiag_rmphi=0
        idiag_dtv=0; idiag_dtdiffus=0; idiag_Rmesh=0; idiag_Rmesh3=0
        idiag_maxadvec=0
      endif
!
!  iname runs through all possible names that may be listed in print.in.
!
      if (lroot.and.ip<14) print*,'rprint_register: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'t',idiag_t)
        call parse_name(iname,cname(iname),cform(iname),'it',idiag_it)
        call parse_name(iname,cname(iname),cform(iname),'dt',idiag_dt)
        call parse_name(iname,cname(iname),cform(iname),'dtv',idiag_dtv)
        call parse_name(iname,cname(iname),cform(iname),'dtdiffus',idiag_dtdiffus)
        call parse_name(iname,cname(iname),cform(iname),'Rmesh',idiag_Rmesh)
        call parse_name(iname,cname(iname),cform(iname),'Rmesh3',idiag_Rmesh3)
        call parse_name(iname,cname(iname),cform(iname),'maxadvec',idiag_maxadvec)
        call parse_name(iname,cname(iname),cform(iname),&
            'walltime',idiag_walltime)
        call parse_name(iname,cname(iname),cform(iname),&
            'timeperstep',idiag_timeperstep)
      enddo
!
!  phi-averages
!
      if (nnamerz>0) then
!
!  Expand some shorthand labels.
!
        call expand_cname(cnamerz,nnamerz,'uumphi','urmphi','upmphi','uzmphi')
        if (lname_is_present(cnamerz,'upmphi')) then
          call expand_cname(cnamerz,nnamerz,'uusphmphi','ursphmphi','uthmphi')
        else
          call expand_cname(cnamerz,nnamerz,'uusphmphi','ursphmphi','uthmphi','upmphi')
        endif
        call expand_cname(cnamerz,nnamerz,'bbmphi','brmphi','bpmphi','bzmphi')
        if (lname_is_present(cnamerz,'bpmphi')) then
          call expand_cname(cnamerz,nnamerz,'bbsphmphi','brsphmphi','bthmphi')
        else
          call expand_cname(cnamerz,nnamerz,'bbsphmphi','brsphmphi','bthmphi','bpmphi')
        endif
        call expand_cname(cnamerz,nnamerz,'uxbmphi','uxbrmphi','uxbpmphi','uxbzmphi')
        call expand_cname(cnamerz,nnamerz,'jxbmphi','jxbrmphi','jxbpmphi','jxbzmphi')
!
!  Some generic quantities (mostly coordinates for debugging).
!
        do irz=1,nnamerz
          call parse_name(irz,cnamerz(irz),cformrz(irz),'rcylmphi',idiag_rcylmphi)
          call parse_name(irz,cnamerz(irz),cformrz(irz),'phimphi',idiag_phimphi)
          call parse_name(irz,cnamerz(irz),cformrz(irz),'zmphi',idiag_zmphi)
          call parse_name(irz,cnamerz(irz),cformrz(irz),'rmphi',idiag_rmphi)
        enddo
!
!  Output in phiavg.list the list of fields after the taking into
!  account of possible shorthands in phiaver.in
!
        if (lroot) then
          open(11,file=trim(datadir)//'/averages/phiavg.list',status='unknown')
          do irz=1,nnamerz
            write(11,'(A30)') cnamerz(irz)
          enddo
          close(11)
        endif
      endif
!
    endsubroutine rprint_general
!***********************************************************************
endmodule Register
