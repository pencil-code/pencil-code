! $Id: register.f90,v 1.109 2003-11-14 15:28:11 theine Exp $

!!!  A module for setting up the f-array and related variables (`register' the
!!!  entropy, magnetic, etc modules).


module Register

  implicit none 

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
      use Mpicomm,      only: mpicomm_init, stop_it
      use Sub 
      use Param_IO,     only: get_datadir, get_snapdir
      use IO,           only: register_io
      use Gravity,      only: register_gravity
      use Hydro,        only: register_hydro
      use Density,      only: register_density
      use Forcing,      only: register_forcing
      use Entropy,      only: register_entropy
      use Magnetic,     only: register_magnetic
      use Radiation,    only: register_radiation
      use Ionization,   only: register_ionization
      use Pscalar,      only: register_pscalar
      use Dustdensity,  only: register_dustdensity
      use Dustvelocity, only: register_dustvelocity
      use CosmicRay,    only: register_cosmicray
      use Interstellar, only: register_interstellar
      use Shear,        only: register_shear
      use Viscosity,    only: register_viscosity
      use Special,      only: register_special
!
!  initialize all mpi stuff
!
      call mpicomm_init
!
!  initialize nvar; is increased by the following routines
!
      nvar = 0 
      naux = 0 
!
!  Writing files for use with IDL
!
      if (lroot) then
        open(15,file=trim(datadir)//'/def_var.pro')
        open(4,file=trim(datadir)//'/variables.pro')
        write(4,*) 'close,1'
        write(4,*) "openr,1, datadir+'/'+varfile, /F77"
        write(4,*) 'readu,1 $'
      endif
!
      call register_io
!
      call register_hydro
      call register_density
      call register_viscosity
      call register_forcing
      call register_entropy
      call register_magnetic
      call register_radiation
      call register_ionization
      call register_pscalar
      call register_dustvelocity
      call register_dustdensity
      call register_gravity
      call register_cosmicray
      call register_interstellar
      call register_shear
      call register_special
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
    endsubroutine register_modules
!***********************************************************************
    subroutine initialize_modules(f,lstart)
!
!  Call initialization routines, i.e. initialize physics and technical
!  modules. This implies some preparation of auxiliary quantities, often
!  based on input parameters. This routine is called by run.x (but not by
!  start.x) initially and each time the run parameters have been reread.
!
!  6-nov-01/wolf: coded
! 23-feb-03/axel: added physical constants conversion
!  7-oct-03/david: initialize_gravity before density, etc (its needed there)
!
      use Cdata
      use Param_IO
      use Print
!      use Hydro
!      use Density
      use Timeavg,      only: initialize_timeavg
      use Gravity,      only: initialize_gravity
      use Forcing,      only: initialize_forcing
      use Entropy,      only: initialize_entropy
!      use Magnetic
      use Radiation,    only: initialize_radiation
      use Ionization,   only: initialize_ionization
      use Pscalar,      only: initialize_pscalar
!      use Dustdensity,  only: initialize_dustdensity
!      use Dustvelocity, only: initialize_dustvelocity
      use CosmicRay,    only: initialize_cosmicray
      use Interstellar, only: initialize_interstellar
      use Shear,        only: initialize_shear
      use Viscosity,    only: initialize_viscosity
      use Special,      only: initialize_special

      real, dimension(mx,my,mz,mvar+maux) :: f
      logical :: lstart
!
!  Defaults for some logicals; will later be set to true if needed
      lneed_sij = .false.
      lneed_glnrho = .false.
!
!  evaluate physical units
!  used currently only in ionization, but later also in
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
        if(lroot.and.lionization) print*,'initialize_modules: unit_velocity, unit_density, etc, are in cgs'
        hbar=hbar_cgs/(unit_energy*unit_time)
        k_B=k_B_cgs/(unit_energy/unit_temperature)
        m_p=m_p_cgs/unit_mass
        m_e=m_e_cgs/unit_mass
        eV=eV_cgs/unit_energy
        sigmaH_=sigmaH_cgs/unit_length**2
        sigmaSB=sigmaSB_cgs/(unit_flux/unit_temperature**4)
        kappa_es=kappa_es_cgs/(unit_length**2/unit_mass)
      elseif (unit_system=='SI') then
        if(lroot.and.lionization) print*,'initialize_modules: unit_velocity, unit_density, etc, are in SI'
        hbar=hbar_cgs*1e-7/(unit_energy*unit_time)
        k_B=1e-7*k_B_cgs/(unit_energy/unit_temperature)
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
!  print parameters in code units, but only when used
!
      if (lroot) then
         if (lionization.or.lradiation.or.lradiation_ray.or.linterstellar) then
            print'(a,1p,4e14.6)',' register: k_B,m_p,m_e,eV=',k_B,m_p,m_e,eV
         endif
      endif
!
!  set gamma1, cs20, and lnrho0
!  (used currently for non-dimensional equation of state)
!
      gamma1=gamma-1.
      cs20=cs0**2
      lnrho0=alog(rho0)
!
!  run initialization of individual modules
!
!      call initialize_io
      call initialize_prints
!ajwm timeavg needs tidying to be similar structure to other modules
      call initialize_timeavg(f) ! initialize time averages
!
      call initialize_gravity
      call initialize_hydro
      call initialize_density
      call initialize_forcing(lstart)  ! get random seed from file, ..
      call initialize_entropy          ! calculate radiative conductivity, etc.
!      call initialize_magnetic
      call initialize_radiation
      call initialize_ionization
      call initialize_pscalar(f)
!      call initialize_dustvelocity
!      call initialize_dustdensity
      call initialize_cosmicray(f)
      call initialize_interstellar(lstart)
      call initialize_shear
      call initialize_viscosity
      call initialize_special
!
!  timestep: if dt=0 (ie not initialized), ldt=.true.
!
!ajwm should this be moved to timestep.f90 as run_hooks_timestep() ??
!AB: maybe not, because initialize_modules can also be run from start.f90,
!AB: which has no knowledge of timestep.f90
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
    subroutine rprint_list(lreset)
!
!  read variables to print and to calculate averages of from control files
!
!   3-may-01/axel: coded
!
      use Cdata
      use Param_IO
      use Hydro,        only: rprint_hydro
      use Entropy,      only: rprint_entropy
      use Magnetic,     only: rprint_magnetic
      use Radiation,    only: rprint_radiation
      use Ionization,   only: rprint_ionization
      use Pscalar,      only: rprint_pscalar
      use Dustvelocity, only: rprint_dustvelocity
      use Dustdensity,  only: rprint_dustdensity
      use CosmicRay,    only: rprint_cosmicray
      use Gravity,      only: rprint_gravity
      use Special,      only: rprint_special
!
      integer :: iname,inamev,inamez,inamexy,inamerz
      integer :: ix_,iy_,iz_,iz2_,io_stat
      integer :: isubstract
      character(len=30) :: cnamev_
      logical :: lreset,exist
!
!  read in the list of variables to be printed
!
      open(1,file='print.in')
      do iname=1,mname
        read(1,*,end=99) cname(iname)
      enddo
99    nname=iname-1
      if (lroot.and.ip<14) print*,'rprint_list: nname=',nname
      close(1)
!
!  read in the list of variables for video slices
!
      inquire(file='video.in',exist=exist)
      if (exist) then
        lwrite_slices=.true.
        isubstract=0
        ix=n1; iy=m1; iz=n1; iz2=n2
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
      if (lroot.and.ip<14) print*,'rprint_list: cnamev=',cnamev(1:nnamev)
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
      else
        lwrite_phiaverages = .false. ! switch phiaverages off
      endif
      if (lroot.and.ip<14) print*,'rprint_list: nnamerz=',nnamerz
!
!  check which variables are set
!  For the convenience of idl users, the indices of variables in
!  the f-array and the time_series.dat files are written to data/index.pro
!
      if (lroot) open(3,file=trim(datadir)//'/index.pro')
      call rprint_general     (lreset,LWRITE=lroot)
      call rprint_hydro       (lreset,LWRITE=lroot)
      call rprint_density     (lreset,LWRITE=lroot)
      call rprint_entropy     (lreset,LWRITE=lroot)
      call rprint_magnetic    (lreset,LWRITE=lroot)
      call rprint_radiation   (lreset,LWRITE=lroot)
      call rprint_ionization  (lreset,LWRITE=lroot)
      call rprint_pscalar     (lreset,LWRITE=lroot)
      call rprint_dustvelocity(lreset,LWRITE=lroot)
      call rprint_dustdensity (lreset,LWRITE=lroot)
      call rprint_cosmicray   (lreset,LWRITE=lroot)
      call rprint_gravity     (lreset,LWRITE=lroot)
      call rprint_special     (lreset,LWRITE=lroot)
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
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_t=0;i_it=0;i_dt=0;i_dtc=0;i_walltime=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'t',i_t)
        call parse_name(iname,cname(iname),cform(iname),'it',i_it)
        call parse_name(iname,cname(iname),cform(iname),'dt',i_dt)
        call parse_name(iname,cname(iname),cform(iname),'dtc',i_dtc)
        call parse_name(iname,cname(iname),cform(iname),'walltime',i_walltime)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_t=',i_t
        write(3,*) 'i_it=',i_it
        write(3,*) 'i_dt=',i_dt
        write(3,*) 'i_dtc=',i_dtc
        write(3,*) 'i_walltime=',i_walltime
        write(3,*) 'nname=',nname
!ajwm Not really the correct place to put this...?
        write(3,*) 'ishock=',ishock
      endif
!
    endsubroutine rprint_general
!***********************************************************************

endmodule Register

!!! End of file register.f90
