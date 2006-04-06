! $Id: register.f90,v 1.166 2006-04-06 19:02:11 theine Exp $

!!!  A module for setting up the f-array and related variables (`register' the
!!!  entropy, magnetic, etc modules).


module Register

  implicit none 

  private
  
  public :: register_modules, initialize_modules, rprint_list
  public :: choose_pencils

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
      use Mpicomm,         only: mpicomm_init,stop_it,stop_it_if_any
      use Sub 
      use Param_IO,        only: get_datadir,get_snapdir
      use IO,              only: register_io
      use Global,          only: register_global
      use EquationOfState, only: register_eos
      use Shock,           only: register_shock
      use Gravity,         only: register_gravity
      use Hydro,           only: register_hydro
      use Density,         only: register_density
      use Forcing,         only: register_forcing
      use Entropy,         only: register_entropy
      use Magnetic,        only: register_magnetic
      use Testfield,       only: register_testfield
      use Radiation,       only: register_radiation
      use Pscalar,         only: register_pscalar
      use Chiral,          only: register_chiral
      use Dustdensity,     only: register_dustdensity
      use Dustvelocity,    only: register_dustvelocity
      use CosmicRay,       only: register_cosmicray
      use CosmicrayFlux,   only: register_cosmicrayflux
      use Interstellar,    only: register_interstellar
      use Shear,           only: register_shear
      use Viscosity,       only: register_viscosity
      use Special,         only: register_special
      use Planet,          only: register_planet
!
      logical :: ioerr
!
!  initialize all mpi stuff
!
      call mpicomm_init
!
!  initialize nvar; is increased by the following routines
!
      nvar     = 0 
      naux     = 0
      naux_com = 0
!
!  Writing files for use with IDL
!
      ioerr = .true.            ! will be overridden unless we go 911
      if (lroot) then
        open(15,FILE=trim(datadir)//'/def_var.pro',ERR=911)
        open(4,FILE=trim(datadir)//'/variables.pro',ERR=911)
        write(4,*) 'close,1'
        write(4,*) "openr,1, datadir+'/'+varfile, /F77"
        write(4,*) 'readu,1 $'
      endif
      ioerr = .false.
!
      call register_io
      call register_global
      call register_eos

      call register_shock
      call register_viscosity
      call register_hydro
      call register_gravity
      call register_density
      call register_forcing
      call register_entropy
      call register_magnetic
      call register_testfield
      call register_radiation
      call register_pscalar
      call register_chiral
      call register_dustvelocity
      call register_dustdensity
      call register_cosmicray
      call register_cosmicrayflux
      call register_interstellar
      call register_shear
      call register_special
      call register_planet
!
!  Writing files for use with IDL
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
        call stop_it('Initialize: nvar /= mvar. Fix mvar in cparam.local')
      endif
!
!  initialize headt for root processor only
!
      if (lroot) headt=.true.
!
!  overwrite datadir from datadir.in, if that exists
!
      call get_datadir(datadir)
      call get_snapdir(datadir_snap)
!
!  Something went wrong. Catches cases that would make mpich 1.x hang,
!  provided that this is the first attempt to write a file
!

911   call stop_it_if_any(ioerr, &
          "Cannot open data/def_var.pro for writing" // &
          " -- is data/ visible from root node?")
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
      use Sub, only: remove_zprof
      use Param_IO
      use Print
!      use Hydro
!      use Density
      use Timeavg,         only: initialize_timeavg
      use EquationOfState, only: initialize_eos
      use CosmicrayFlux,   only: initialize_cosmicrayflux
      use Hydro,           only: initialize_hydro
      use Density,         only: initialize_density
      use Shock,           only: initialize_shock
      use Gravity,         only: initialize_gravity
      use Forcing,         only: initialize_forcing
      use Entropy,         only: initialize_entropy
      use Magnetic,        only: initialize_magnetic
      use Testfield,       only: initialize_testfield
      use Radiation,       only: initialize_radiation
      use Pscalar,         only: initialize_pscalar
      use Chiral,          only: initialize_chiral
      use Dustvelocity,    only: initialize_dustvelocity
      use Dustdensity,     only: initialize_dustdensity
      use CosmicRay,       only: initialize_cosmicray
      use Interstellar,    only: initialize_interstellar
      use Shear,           only: initialize_shear
      use Viscosity,       only: initialize_viscosity
      use Special,         only: initialize_special
      use Planet,          only: initialize_planet
!     use Timestep,        only: border_profiles

      real, dimension(mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!
!  Defaults for some logicals; will later be set to true if needed
      lpenc_requested(:) = .false.
!
!  evaluate physical units
!  used currently only in eos, but later also in
!  the interstellar and radiation modules, for example
!
      unit_mass=unit_density*unit_length**3
      unit_energy=unit_mass*unit_velocity**2
      unit_time=unit_length/unit_velocity
      unit_flux=unit_energy/(unit_length**2*unit_time)
!
!  convert physical constants
!
      if (unit_system=='cgs') then
        if(lroot.and.leos_ionization.and.ip<14) print*,'initialize_modules: ' &
          //'unit_velocity, unit_density, etc, are in cgs'
        hbar=hbar_cgs/(unit_energy*unit_time)
        k_B=k_B_cgs/(unit_energy/unit_temperature)
        sigmaSB=sigmaSB_cgs/(unit_flux/unit_temperature**4)
        m_u=m_u_cgs/unit_mass
        m_p=m_p_cgs/unit_mass
        m_e=m_e_cgs/unit_mass
        eV=eV_cgs/unit_energy
        sigmaH_=sigmaH_cgs/unit_length**2
        kappa_es=kappa_es_cgs/(unit_length**2/unit_mass)
      elseif (unit_system=='SI') then
        if(lroot.and.leos_ionization) print*,&
            'initialize_modules: unit_velocity, unit_density, etc, are in SI'
        hbar=hbar_cgs*1e-7/(unit_energy*unit_time)
        k_B=1e-7*k_B_cgs/(unit_energy/unit_temperature)
        m_u=m_u_cgs*1e-3/unit_mass
        m_p=m_p_cgs*1e-3/unit_mass
        m_e=m_e_cgs*1e-3/unit_mass
        eV=eV_cgs*1e-7/unit_energy
        sigmaH_=sigmaH_cgs*1e-4/unit_length**2
        sigmaSB=sigmaSB_cgs*1e-3/(unit_flux/unit_temperature**4)
        kappa_es=kappa_es_cgs*1e-1/(unit_length**2/unit_mass)
      endif
!
!  calculate additional constants
!
      m_H=m_p+m_e
      m_He=3.97153*m_H
      chiH=13.6*eV
      chiH_=0.75*eV        
!
!  run initialization of individual modules
!   allow initialize_eos to go early so that it may change the unit temperature.
!     IFF it does it must then be careful not to use any of the consants above
!     that depend upon unit_temperature without extreme care!
!
!      call initialize_io
      call initialize_eos()
!
!  recalculate anything that depends upon unit_temperature
!
      if (unit_system=='cgs') then
        k_B=k_B_cgs/(unit_energy/unit_temperature)
        sigmaSB=sigmaSB_cgs/(unit_flux/unit_temperature**4)
      elseif (unit_system=='SI') then
        k_B=1e-7*k_B_cgs/(unit_energy/unit_temperature)
        sigmaSB=sigmaSB_cgs*1e-3/(unit_flux/unit_temperature**4)
      endif
!
!  print parameters in code units, but only when used
!
      if (lroot.and.ip<14) then
         if (leos_ionization.or.lradiation.or.lradiation_ray.or.linterstellar) then
            write(*,'(a,1p,4e14.6)') ' register: k_B,m_p,m_e,eV=',k_B,m_p,m_e,eV
         endif
      endif

!
!  run rest of initialization of individual modules
!
      call initialize_prints()
      call initialize_timeavg(f) ! initialize time averages
!
      call initialize_gravity()
      call initialize_density(f,lstarting)
      call initialize_hydro(f,lstarting)
      call initialize_forcing(lstarting)   ! get random seed from file, ..
      call initialize_entropy(f,lstarting) ! calculate radiative conductivity,..
      call initialize_magnetic(f,lstarting)
      call initialize_testfield(f)
      call initialize_radiation()
      call initialize_pscalar(f)
      call initialize_chiral(f)
      call initialize_dustvelocity()
      call initialize_dustdensity()
      call initialize_cosmicray(f)
      call initialize_cosmicrayflux(f)
      call initialize_interstellar(lstarting)
      call initialize_shear()
      call initialize_shock(lstarting)
      call initialize_viscosity(lstarting)
      call initialize_special(f)
      call initialize_planet(f,lstarting) !will need f for torque
!
! Store the value of impossible for use in IDL
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,*) 'impossible=',impossible
        close (1)
      endif
!
!----------------------------------------------------------------------------
!  Coordinate-related issues: nonuniform meshes, different corrdinate systems
!
!  The following is only kept for backwards compatibility with
!  an old grid.dat.
!
!AB:  Removed for now, because
!AB:  (i) it should be taken care of by grid.f90
!AB:  (ii) it doesn't take care of 1-D cases with dx=0.
!
!     if (lequidist(1)) dx_1=1./dx
!     if (lequidist(2)) dy_1=1./dy
!     if (lequidist(3)) dz_1=1./dz
!
!  For spherical coordinate system, calculate 1/r, cot(theta)/r, etc
!
      if (coord_system=='cartesian') then
        lspherical=.false.
        lcylindric=.false.
      elseif (coord_system=='spherical') then
        lspherical=.true.
        lcylindric=.false.
        if (x(l1)==0.) then
          r1_mn(2:)=1./x(l1+1:l2)
          r1_mn(1)=0.
        else
          r1_mn=1./x(l1:l2)
        endif
      endif
!
!  print the value for which output is being produced
!  (Have so far only bothered about single processor output.)
!
      if (lroot) then
        print*,'x(lpoint),y(mpoint),z(npoint)=',x(lpoint),y(mpoint),z(npoint)
      endif
!
!  DOCUMENT ME
!  AB: should check whether this can come under initialize_modules
!
!       call border_profiles()
!
!  cleanup profile files
!
      call remove_zprof()
      lwrite_prof=.true.
!
!----------------------------------------------------------------------------
!  timestep: distinguish two cases,
!  (a) dt explicitly given in run.in -> ldt=.false.
!  (b) dt not given in run.in        -> ldt=.true.  -> calculate dt dynamically
!  Note that ldt will not change unless you RELOAD parameters.
!
! Why is this here?...
!   ajwm should this be moved to timestep.f90 as run_hooks_timestep() ??
!   AB   maybe not, because initialize_modules can also be run from start.f90,
!   AB   which has no knowledge of timestep.f90
!
      ldt = (dt==0.)            ! need to calculate dt dynamically?
      if (lroot .and. ip<14) then
        if (ldt) then
          print*,'timestep based on CFL cond; cdt=',cdt
        else
          print*, 'absolute timestep dt=', dt
        endif
      endif
!
    endsubroutine initialize_modules
!***********************************************************************
    subroutine choose_pencils()
!
!  Find out which pencils are needed for all time-steps and also for
!  diagnostics only. Also takes care of interdependant pencils.
!
!  20-11-04/anders: coded
!
      use Cdata
!
      integer :: i
!
      if (lroot) print*, 'choose_pencils: Finding out which pencils '// &
          'are needed for the pencil_case'
!
!  Must set all pencil arrays to false in case of reload or reinit
!
      lpenc_requested=.false.
      lpenc_diagnos=.false.
      lpenc_diagnos2d=.false.
!
!  Find out which pencils are needed for the pencil case.
!
      call pencil_criteria()
!
!  Set interdependent pencils
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
!  20-11-04/anders: coded
!
      use EquationOfState, only: pencil_criteria_eos
      use Hydro, only: pencil_criteria_hydro
      use Density, only: pencil_criteria_density
      use Shock, only: pencil_criteria_shock
      use Viscosity, only: pencil_criteria_viscosity
      use Entropy, only: pencil_criteria_entropy
      use Gravity, only: pencil_criteria_gravity
      use Pscalar, only: pencil_criteria_pscalar
      use Dustvelocity, only: pencil_criteria_dustvelocity
      use Dustdensity, only: pencil_criteria_dustdensity
      use Magnetic, only: pencil_criteria_magnetic
      use Testfield, only: pencil_criteria_testfield
      use Cosmicray, only: pencil_criteria_cosmicray
      use Cosmicrayflux, only: pencil_criteria_cosmicrayflux
      use Chiral, only: pencil_criteria_chiral
      use Radiation, only: pencil_criteria_radiation
      use Interstellar, only: pencil_criteria_interstellar
      use Planet, only: pencil_criteria_planet
!
      call pencil_criteria_density()
      call pencil_criteria_eos()
      call pencil_criteria_hydro()
      call pencil_criteria_shock()
      call pencil_criteria_viscosity()
      call pencil_criteria_entropy()
      call pencil_criteria_gravity()
      call pencil_criteria_pscalar()
      call pencil_criteria_interstellar()
      call pencil_criteria_dustvelocity()
      call pencil_criteria_dustdensity()
      call pencil_criteria_magnetic()
      call pencil_criteria_testfield()
      call pencil_criteria_cosmicray()
      call pencil_criteria_cosmicrayflux()
      call pencil_criteria_chiral()
      call pencil_criteria_radiation()
      call pencil_criteria_planet()
!    
    endsubroutine pencil_criteria
!***********************************************************************
    subroutine pencil_interdep(lpencil_in)
!
!  Find out about interdependant pencils. Each module knows what its own
!  pencils depend on. The dependency only needs to be specified one level
!  up, since this subroutine is called several times (currently three).
!
!
!  20-11-04/anders: coded
!
      use Cdata
      use EquationOfState, only: pencil_interdep_eos
      use Hydro, only: pencil_interdep_hydro
      use Density, only: pencil_interdep_density
      use Shock, only: pencil_interdep_shock
      use Viscosity, only: pencil_interdep_viscosity
      use Entropy, only: pencil_interdep_entropy
      use Gravity, only: pencil_interdep_gravity
      use Magnetic, only: pencil_interdep_magnetic
      use Testfield, only: pencil_interdep_testfield
      use Pscalar, only: pencil_interdep_pscalar
      use Dustvelocity, only: pencil_interdep_dustvelocity
      use Dustdensity, only: pencil_interdep_dustdensity
      use Cosmicray, only: pencil_interdep_cosmicray
      use Cosmicrayflux, only: pencil_interdep_cosmicrayflux
      use Chiral, only: pencil_interdep_chiral
      use Radiation, only: pencil_interdep_radiation
      use Particles_main, only: particles_pencil_interdep
!      
      logical, dimension (npencils) :: lpencil_in
!
      call pencil_interdep_density(lpencil_in)
      call pencil_interdep_eos(lpencil_in)
      call pencil_interdep_hydro(lpencil_in)
      call pencil_interdep_shock(lpencil_in)
      call pencil_interdep_viscosity(lpencil_in)
      call pencil_interdep_entropy(lpencil_in)
      call pencil_interdep_gravity(lpencil_in)
      call pencil_interdep_dustvelocity(lpencil_in)
      call pencil_interdep_dustdensity(lpencil_in)
      call pencil_interdep_pscalar(lpencil_in)
      call pencil_interdep_magnetic(lpencil_in)
      call pencil_interdep_testfield(lpencil_in)
      call pencil_interdep_cosmicray(lpencil_in)
      call pencil_interdep_cosmicrayflux(lpencil_in)
      call pencil_interdep_chiral(lpencil_in)
      call pencil_interdep_radiation(lpencil_in)
      if (lparticles) call particles_pencil_interdep(lpencil_in)
!    
    endsubroutine pencil_interdep
!***********************************************************************
    subroutine rprint_list(lreset)
!
!  read variables to print and to calculate averages of from control files
!
!   3-may-01/axel: coded
!
      use Cdata
      use Param_IO
      use Hydro,           only: rprint_hydro
      use Density,         only: rprint_density
      use Forcing,         only: rprint_forcing
      use Entropy,         only: rprint_entropy
      use Magnetic,        only: rprint_magnetic
      use Testfield,       only: rprint_testfield
      use Radiation,       only: rprint_radiation
      use EquationOfState, only: rprint_eos
      use Pscalar,         only: rprint_pscalar
      use Chiral,          only: rprint_chiral
      use Interstellar,    only: rprint_interstellar
      use Dustvelocity,    only: rprint_dustvelocity
      use Dustdensity,     only: rprint_dustdensity
      use CosmicRay,       only: rprint_cosmicray
      use CosmicRayFlux,   only: rprint_cosmicrayflux
      use Gravity,         only: rprint_gravity
      use Special,         only: rprint_special
      use Shock,           only: rprint_shock
      use Viscosity,       only: rprint_viscosity
      use Shear,           only: rprint_shear
      use Planet,          only: rprint_planet
!
      integer :: iname,inamev,inamez,inamey,inamex,inamexy,inamexz,inamerz
      integer :: ix_,iy_,iz_,iz2_,io_stat,iname_tmp
      integer :: isubstract
      logical :: lreset,exist
      character (LEN=30) :: cname_tmp
!
!  read in the list of variables to be printed
!  recognize "!" and "#" as comments
!
      open(1,file='print.in')
      iname=0
      do iname_tmp=1,mname
        read(1,*,end=99) cname_tmp
        if (cname_tmp(1:1)/='!'.and.cname_tmp(1:1)/='#') then
          iname=iname+1
          cname(iname)=cname_tmp
        endif
      enddo
99    nname=iname
      if (lroot.and.ip<14) print*,'rprint_list: nname=',nname
      close(1)
!
!  read in the list of variables for video slices
!
      inquire(file='video.in',exist=exist)
      if (exist .and. dvid/=0.0) then
        lwrite_slices=.true.
        isubstract=0
        open(1,file='video.in')
        do inamev=1,mnamev
          read(1,*,end=98,iostat=io_stat) ix_,iy_,iz_,iz2_
          if (io_stat/=0) then
            backspace(1)
            read(1,*,end=98) cnamev(inamev-isubstract)
          else
            ix=ix_; iy=iy_; iz=iz_; iz2=iz2_
            isubstract=isubstract+1
          endif
        enddo
98      nnamev=inamev-1-isubstract
        close(1)
      endif
      if (lroot.and.ip<14) print*,'rprint_list: ix,iy,iz,iz2=',ix,iy,iz,iz2
      if (lroot.and.ip<14) print*,'rprint_list: nnamev=',nnamev
!
!  read in the list of variables for xy-averages
!
      inquire(file='xyaver.in',exist=exist)
      if (exist) then
        open(1,file='xyaver.in')
        do inamez=1,mnamez
          read(1,*,end=97) cnamez(inamez)
        enddo
97      nnamez=inamez-1
        close(1)
      endif
      if (lroot.and.ip<14) print*,'rprint_list: nnamez=',nnamez
!
!  read in the list of variables for xz-averages
!
      inquire(file='xzaver.in',exist=exist)
      if (exist) then
        open(1,file='xzaver.in')
        do inamey=1,mnamey
          read(1,*,end=92) cnamey(inamey)
        enddo
92      nnamey=inamey-1
        close(1)
      endif
      if (lroot.and.ip<14) print*,'rprint_list: nnamex=',nnamex
!
!  read in the list of variables for yz-averages
!
      inquire(file='yzaver.in',exist=exist)
      if (exist) then
        open(1,file='yzaver.in')
        do inamex=1,mnamex
          read(1,*,end=93) cnamex(inamex)
        enddo
93      nnamex=inamex-1
        close(1)
      endif
      if (lroot.and.ip<14) print*,'rprint_list: nnamex=',nnamex
!
!  read in the list of variables for y-averages
!
      inquire(file='yaver.in',exist=exist)
      if (exist) then
        open(1,file='yaver.in')
        do inamexz=1,mnamexz
          read(1,*,end=94) cnamexz(inamexz)
        enddo
94      nnamexz=inamexz-1
        close(1)
      else
        lwrite_yaverages = .false. ! switch yaverages off
      endif
      if (lroot.and.ip<14) print*,'rprint_list: nnamexz=',nnamexz
!
!  read in the list of variables for z-averages
!
      inquire(file='zaver.in',exist=exist)
      if (exist) then
        open(1,file='zaver.in')
        do inamexy=1,mnamexy
          read(1,*,end=96) cnamexy(inamexy)
        enddo
96      nnamexy=inamexy-1
        close(1)
      else
        lwrite_zaverages = .false. ! switch zaverages off
      endif
      if (lroot.and.ip<14) print*,'rprint_list: nnamexy=',nnamexy
!
!  read in the list of variables for phi-averages
!
      inquire(file='phiaver.in',exist=exist)
      if (exist) then
        open(1,file='phiaver.in')
        do inamerz=1,mnamerz
          read(1,*,end=95) cnamerz(inamerz)
        enddo
95      nnamerz=inamerz-1
        close(1)
!
      else
        lwrite_phiaverages = .false. ! switch phiaverages off
      endif
      if (lroot.and.ip<14) print*,'rprint_list: nnamerz=',nnamerz
!
!  set logical for 2-D averages
!
      lwrite_2daverages=lwrite_yaverages&
                    .or.lwrite_zaverages&
                    .or.lwrite_phiaverages
!
!  check which variables are set
!  For the convenience of idl users, the indices of variables in
!  the f-array and the time_series.dat files are written to data/index.pro
!
      if (lroot) open(3,file=trim(datadir)//'/index.pro')
      call rprint_general      (lreset,LWRITE=lroot)
      call rprint_hydro        (lreset,LWRITE=lroot)
      call rprint_density      (lreset,LWRITE=lroot)
      call rprint_forcing      (lreset,LWRITE=lroot)
      call rprint_entropy      (lreset,LWRITE=lroot)
      call rprint_magnetic     (lreset,LWRITE=lroot)
      call rprint_testfield    (lreset,LWRITE=lroot)
      call rprint_radiation    (lreset,LWRITE=lroot)
      call rprint_eos          (lreset,LWRITE=lroot)
      call rprint_pscalar      (lreset,LWRITE=lroot)
      call rprint_chiral       (lreset,LWRITE=lroot)
      call rprint_interstellar (lreset,LWRITE=lroot)
      call rprint_dustvelocity (lreset,LWRITE=lroot)
      call rprint_dustdensity  (lreset,LWRITE=lroot)
      call rprint_cosmicray    (lreset,LWRITE=lroot)
      call rprint_cosmicrayflux(lreset,LWRITE=lroot)
      call rprint_gravity      (lreset,LWRITE=lroot)
      call rprint_special      (lreset,LWRITE=lroot)
      call rprint_shock        (lreset,LWRITE=lroot)
      call rprint_viscosity    (lreset,LWRITE=lroot)
      call rprint_shear        (lreset,LWRITE=lroot)
      call rprint_planet       (lreset,LWRITE=lroot)
      if (lroot) close(3)
!
    endsubroutine rprint_list
!***********************************************************************
    subroutine rprint_general(lreset,lwrite)
!
!  reads and registers *general* print parameters
!
!   8-jun-02/axel: adapted from hydro
!
      use Cdata
      use Sub
!
      integer :: iname,irz
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of RELOAD
!  (general variables that are defined in Cdata)
!
      if (lreset) then
        idiag_t=0; idiag_it=0; idiag_dt=0; idiag_walltime=0
        idiag_timeperstep=0
        idiag_rcylmphi=0; idiag_phimphi=0; idiag_zmphi=0; idiag_rmphi=0
        idiag_dtv=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'t',idiag_t)
        call parse_name(iname,cname(iname),cform(iname),'it',idiag_it)
        call parse_name(iname,cname(iname),cform(iname),'dt',idiag_dt)
        call parse_name(iname,cname(iname),cform(iname),'dtv',idiag_dtv)
        call parse_name(iname,cname(iname),cform(iname),&
            'walltime',idiag_walltime)
        call parse_name(iname,cname(iname),cform(iname),&
            'timeperstep',idiag_timeperstep)
      enddo
!
!  phi-averages
!
      !
      !  expand some shorthand labels 
      !
      call expand_cname(cnamerz,nnamerz,'uumphi','urmphi','upmphi','uzmphi')
      call expand_cname(cnamerz,nnamerz,'bbmphi','brmphi','bpmphi','bzmphi')
      call expand_cname(cnamerz,nnamerz,'uxbmphi','uxbrmphi','uxbpmphi','uxbzmphi')
      !
      !  some generic quantities (mostly coordinates for debugging)
      !
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rcylmphi',idiag_rcylmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'phimphi', idiag_phimphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'zmphi',   idiag_zmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rmphi',   idiag_rmphi)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_t=',idiag_t
        write(3,*) 'i_it=',idiag_it
        write(3,*) 'i_dt=',idiag_dt
        write(3,*) 'i_walltime=',idiag_walltime
        write(3,*) 'i_timeperstep=',idiag_timeperstep
        write(3,*) 'i_rcylmphi=',idiag_rcylmphi
        write(3,*) 'i_phimphi=',idiag_phimphi
        write(3,*) 'i_zmphi=',idiag_zmphi
        write(3,*) 'i_rmphi=',idiag_rmphi
        write(3,*) 'i_dtv=',idiag_dtv
        write(3,*) 'nname=',nname
      endif
!
    endsubroutine rprint_general
!***********************************************************************

endmodule Register

!!! End of file register.f90
