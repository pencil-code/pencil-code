! $Id$
!
!  This module contains all the main structure needed for particles.
!
module Particles_main
!
  use Cdata
  use Messages
  use Particles
  use Particles_cdata
  use Particles_collisions
  use Particles_nbody
  use Particles_number
  use Particles_radius
  use Particles_spin
  use Particles_selfgravity
  use Particles_stalker
  use Particles_sub
  use Particles_viscosity
!
  implicit none
!
  include 'particles_main.h'
!
  real, dimension (mpar_loc,mpvar) :: fp, dfp
  integer, dimension (mpar_loc,3) :: ineargrid
!
  contains
!***********************************************************************
    subroutine particles_register_modules()
!
!  Register particle modules.
!
!  07-jan-05/anders: coded
!
      call register_particles          ()
      call register_particles_radius   ()
      call register_particles_spin     ()
      call register_particles_number   ()
      call register_particles_selfgrav ()
      call register_particles_nbody    ()
      call register_particles_viscosity()
!
    endsubroutine particles_register_modules
!***********************************************************************
    subroutine particles_rprint_list(lreset)
!
!  Read names of diagnostic particle variables to print out during run.
!
!  07-jan-05/anders: coded
!
      logical :: lreset
!
      if (lroot) open(3, file=trim(datadir)//'/index.pro', &
          STATUS='old', POSITION='append')
      call rprint_particles           (lreset,LWRITE=lroot)
      call rprint_particles_radius    (lreset,LWRITE=lroot)
      call rprint_particles_spin      (lreset,LWRITE=lroot)
      call rprint_particles_number    (lreset,LWRITE=lroot)
      call rprint_particles_selfgrav  (lreset,LWRITE=lroot)
      call rprint_particles_nbody     (lreset,LWRITE=lroot)
      call rprint_particles_viscosity (lreset,LWRITE=lroot)
      call rprint_particles_collisions(lreset,LWRITE=lroot)
      if (lroot) close(3)
!
    endsubroutine particles_rprint_list
!***********************************************************************
    subroutine particles_initialize_modules(f,lstarting)
!
!  Initialize particle modules.
!
!  07-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Check if there is enough total space allocated for particles.
!
      if (ncpus*mpar_loc<npar) then
        if (lroot) then
          print*, 'particles_initialize_modules: '// &
          'total number of particle slots available at the processors '// &
          'is smaller than the number of particles!'
          print*, 'particles_initialize_modules: npar/ncpus=', npar/ncpus
          print*, 'particles_initialize_modules: mpar_loc-ncpus*npar_mig=', &
              mpar_loc-ncpus*npar_mig
        endif
        call fatal_error('particles_initialize_modules','')
      endif
!
      call initialize_particles           (f,lstarting)
      call initialize_particles_radius    (f,lstarting)
      call initialize_particles_spin      (f,lstarting)
      call initialize_particles_number    (f,lstarting)
      call initialize_particles_selfgrav  (f,lstarting)
      call initialize_particles_nbody     (f,lstarting)
      call initialize_particles_viscosity (f,lstarting)
      call initialize_particles_collisions(f,lstarting)
      call initialize_particles_stalker   (f,lstarting)
!
!  Make sure all requested interpolation variables are available.
!
      call interpolation_consistency_check()
!
!  Set internal and external radii of particles.
!
      if (rp_int==-impossible .and. r_int>epsi) rp_int=r_int
      if (rp_ext==-impossible)                  rp_ext=r_ext
!
    endsubroutine particles_initialize_modules
!***********************************************************************
    subroutine particles_init(f)
!
!  Set up initial condition for particle modules.
!
!  07-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent (out) :: f
!
      call init_particles(f,fp,ineargrid)
      if (lparticles_radius) call init_particles_radius(f,fp)
      if (lparticles_spin)   call init_particles_spin(f,fp)
      if (lparticles_number) call init_particles_number(f,fp)
      if (lparticles_nbody)  call init_particles_nbody(f,fp)
!
    endsubroutine particles_init
!***********************************************************************
    subroutine particles_read_snapshot(filename)
!
!  Read particle snapshot from file.
!
!  07-jan-05/anders: coded
!
      character (len=*) :: filename
!
      call input_particles(filename,fp,npar_loc,ipar)
!
    endsubroutine particles_read_snapshot
!***********************************************************************
    subroutine particles_write_snapshot(chsnap,enum,flist)
!
!  Write particle snapshot to file.
!
!  07-jan-05/anders: coded
!
      logical :: enum
      character (len=*) :: chsnap,flist
      optional :: flist
!
      logical :: lsnap
!
      if (present(flist)) then
        call wsnap_particles(chsnap,fp,enum,lsnap,dsnap_par_minor,dsnap_par,npar_loc,ipar,flist)
      else
        call wsnap_particles(chsnap,fp,enum,lsnap,dsnap_par_minor,dsnap_par,npar_loc,ipar)
      endif
!
    endsubroutine particles_write_snapshot
!***********************************************************************
    subroutine particles_write_dsnapshot(chsnap)
!
!  Write particle derivative snapshot to file.
!
!  07-jan-05/anders: coded
!
      character (len=*) :: chsnap
!
      call wsnap_particles(chsnap,dfp,.false.,.false.,0.0,0.0,npar_loc,ipar,nobound=.true.)
!
    endsubroutine particles_write_dsnapshot
!***********************************************************************
    subroutine particles_write_pdim(filename)
!
!  Write npar and mpvar to file.
!
!  09-jan-05/anders: coded
!
      character (len=*) :: filename
!
      open(1,file=filename)
        write(1,'(3i9)') npar, mpvar, npar_stalk
      close(1)
!
    endsubroutine particles_write_pdim
!***********************************************************************
    subroutine particles_timestep_first()
!
!  Setup dfp in the beginning of each itsub.
!
!  07-jan-05/anders: coded
!
      if (itsub==1) then
        dfp(1:npar_loc,:)=0.
      else
        dfp(1:npar_loc,:)=alpha_ts(itsub)*dfp(1:npar_loc,:)
      endif
!
    endsubroutine particles_timestep_first
!***********************************************************************
    subroutine particles_timestep_second()
!
!  Time evolution of particle variables.
!
!  07-jan-05/anders: coded
!
      fp(1:npar_loc,:) = fp(1:npar_loc,:) + dt_beta_ts(itsub)*dfp(1:npar_loc,:)
!
      if (lparticles_collisions .and. itsub==itorder) &
          call calc_particles_collisions(fp,ineargrid)
!
    endsubroutine particles_timestep_second
!***********************************************************************
    subroutine particles_boundconds(f)
!
!  Particle boundary conditions and parallel communication.
!
!  16-feb-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  First apply boundary conditions to the newly updated particle positions.
!
      call boundconds_particles(fp,npar_loc,ipar,dfp=dfp)
!
!  Remove particles that are too close to sink particles or sink points.
!  WARNING: ineargrid and the mapped particle density have not been updated
!  yet, and the sink particle subroutine must not rely on those arrays.
!
      if (lparticles)       call remove_particles_sink(f,fp,dfp,ineargrid)
      if (lparticles_nbody) call remove_particles_sink_nbody(f,fp,dfp,ineargrid)
!
!  Create new sink particles or sink points.
!
      if (lparticles)       call create_sink_particles(f,fp,dfp,ineargrid)
      if (lparticles_nbody) call create_sink_particles_nbody(f,fp,dfp,ineargrid)
!
!  Map the particle positions and velocities on the grid for later use.
!
      call map_nearest_grid(fp,ineargrid)
      call map_xxp_grid(f,fp,ineargrid)
      call map_vvp_grid(f,fp,ineargrid)
!
!  Sort particles so that they can be accessed contiguously in the memory.
!
      call sort_particles_imn(fp,ineargrid,ipar,dfp=dfp)
!
!  Distribute the n-body particles across processors
!
      if (lparticles_nbody) call bcast_nbodyarray(fp)
!
    endsubroutine particles_boundconds
!***********************************************************************
    subroutine particles_doprepencil_calc(f,ivar1,ivar2)
!
!  Do some pre-pencil-loop calculation on the f array.
!  The returned indices should be used for recommication of ghost zones.
!
!  11-aug-08/kapelrud: coded
!
      real, dimension(mx,my,mz,mfarray),intent(inout) :: f
      integer, intent(out) :: ivar1, ivar2
!
      if (lparticles_spin) then
        call particles_spin_prepencil_calc(f)
        ivar1=iox
        ivar2=ioz
      else
        ivar1=-1
        ivar2=-1
      endif
!
    endsubroutine particles_doprepencil_calc
!***********************************************************************
    subroutine particles_calc_selfpotential(f,rhs_poisson,rhs_poisson_const, &
        lcontinued)
!
!  Calculate the potential of the dust particles (wrapper).
!
!  13-jun-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: rhs_poisson
      real :: rhs_poisson_const
      logical :: lcontinued
!
      call calc_selfpotential_particles(f,rhs_poisson,rhs_poisson_const, &
          lcontinued)
!
    endsubroutine particles_calc_selfpotential
!***********************************************************************
    subroutine particles_before_boundary(f)
!
!  Calculate particle-related properties before boundary conditions are
!  set.
!
!  07-feb-09/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Calculate the summed gravity due to the massive particles. 
!  Needed because it is too slow to calculate the gravity at the 
!  position of all dust particles. So we calculate the gravity 
!  for the grid and interpolate to the position of the particles.
!
      if (lparticles_nbody)     call calc_nbodygravity_particles(f)
      if (lparticles_viscosity) call calc_particles_viscosity(f,fp,ineargrid)
!
    endsubroutine particles_before_boundary
!***********************************************************************
    subroutine particles_special
!
!  Fetch fp (and fsp) array to special module.
!
!  01-mar-08/wlad: coded
!
      use Special, only: special_calc_particles
!
      call special_calc_particles(fp)
      if (lparticles_nbody) call particles_nbody_special
!
    endsubroutine particles_special
!***********************************************************************
    subroutine particles_pencil_criteria()
!
!  Request pencils for particles.
!
!  20-apr-06/anders: coded
!
      if (lparticles)             call pencil_criteria_particles()
      if (lparticles_radius)      call pencil_criteria_par_radius()
      if (lparticles_spin)        call pencil_criteria_par_spin()
      if (lparticles_number)      call pencil_criteria_par_number()
      if (lparticles_selfgravity) call pencil_criteria_par_selfgrav()
      if (lparticles_nbody)       call pencil_criteria_par_nbody()
!
    endsubroutine particles_pencil_criteria
!***********************************************************************
    subroutine particles_pencil_interdep(lpencil_in)
!
!  Calculate particle pencils.
!
!  15-feb-06/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lparticles)             call pencil_interdep_particles(lpencil_in)
      if (lparticles_selfgravity) call pencil_interdep_par_selfgrav(lpencil_in)
      if (lparticles_nbody)       call pencil_interdep_par_nbody(lpencil_in)
!
    endsubroutine particles_pencil_interdep
!***********************************************************************
    subroutine particles_calc_pencils(f,p)
!
!  Calculate particle pencils.
!
!  14-feb-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lparticles)             call calc_pencils_particles(f,p)
      if (lparticles_selfgravity) call calc_pencils_par_selfgrav(f,p)
      if (lparticles_nbody)       call calc_pencils_par_nbody(f,p)
!
    endsubroutine particles_calc_pencils
!***********************************************************************
    subroutine particles_pde_pencil(f,df,p)
!
!  Dynamical evolution of particle variables.
!
!  20-apr-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent (in) :: p
      intent (inout) :: f, df
!
!  Create shepherd/neighbour list of required.
!
      if (lshepherd_neighbour) &
          call shepherd_neighbour(fp,ineargrid,kshepherd,kneighbour)
!
!  Interpolate required quantities using the predefined policies. Variables
!  are found in interp.
!  (Clean-up should be performed at end of this subroutine!)
!
     call interpolate_quantities(f,fp,ineargrid)
!
!  Dynamical equations.
!
      if (lparticles)           call dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles)           call dvvp_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_radius)    call dap_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_spin)      call dps_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_number)    call dnptilde_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_selfgravity) &
          call dvvp_dt_selfgrav_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_nbody) &
          call dvvp_dt_nbody_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_viscosity) &
          call dvvp_dt_viscosity_pencil(f,df,fp,dfp,ineargrid)
!
      call cleanup_interpolated_quantities()
!
    endsubroutine particles_pde_pencil
!***********************************************************************
    subroutine particles_pde(f,df)
!
!  Dynamical evolution of particle variables.
!
!  07-jan-05/anders: coded
!
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      intent (in)  :: f
      intent (out) :: df
!
!  Write information about local particle environment to file.
!
      if (itsub==1)               call particles_stalker_sub(f,fp,ineargrid)
!
!  Dynamical equations.
!
      if (lparticles)             call dxxp_dt(f,df,fp,dfp,ineargrid)
      if (lparticles)             call dvvp_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_radius)      call dap_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_spin)        call dps_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_number)      call dnptilde_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_selfgravity) call dvvp_dt_selfgrav(f,df,fp,dfp,ineargrid)
      if (lparticles_nbody)       call dxxp_dt_nbody(dfp)
      if (lparticles_nbody)       call dvvp_dt_nbody(f,df,fp,dfp,ineargrid)
!
!  Correct for curvilinear geometry.
!
      call correct_curvilinear
!
    endsubroutine particles_pde
!***********************************************************************    
    subroutine correct_curvilinear
!
!  Curvilinear corrections to acceleration only.
!  Corrections to velocity were already taken into account
!  in the dxx_dp of particles_dust.f90
!
!  15-sep-07/wlad: coded
!
      real :: rad,raddot,phidot,thtdot,sintht,costht
      integer :: k


      do k=1,npar_loc
!
!  Correct acceleration.
!
        if (lcylindrical_coords) then
          rad=fp(k,ixp);raddot=fp(k,ivpx);phidot=fp(k,ivpy)/max(rad,tini)
          dfp(k,ivpx) = dfp(k,ivpx) + rad*phidot**2
          dfp(k,ivpy) = dfp(k,ivpy) - 2*raddot*phidot
        elseif (lspherical_coords) then
          rad=fp(k,ixp)
          sintht=sin(fp(k,iyp));costht=cos(fp(k,iyp))
          raddot=fp(k,ivpx);thtdot=fp(k,ivpy)/max(rad,tini)
          phidot=fp(k,ivpz)/(max(rad,tini)*sintht)
!
          dfp(k,ivpx) = dfp(k,ivpx) &
               + rad*(thtdot**2 + (sintht*phidot)**2)
          dfp(k,ivpy) = dfp(k,ivpy) &
               - 2*raddot*thtdot + rad*sintht*costht*phidot**2
          dfp(k,ivpz) = dfp(k,ivpz) &
               - 2*phidot*(sintht*raddot + rad*costht*thtdot)
        endif
      enddo
!
    endsubroutine correct_curvilinear
!***********************************************************************
    subroutine particles_read_startpars(unit,iostat)
!
!  Read particle parameters from start.in.
!
!  01-sep-05/anders: coded
!  17-aug-08/wlad: added individual check for the modules inside the wrap
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call read_particles_init_pars(unit,iostat)
      if (present(iostat)) then
        if (iostat/=0) then
          call samplepar_startpars('particles_init_pars',iostat); return
        endif
      endif
!
      if (lparticles_radius) then 
        call read_particles_rad_init_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_startpars('particles_rad_init_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_spin) then
        call read_particles_spin_init_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_startpars('particles_spin_init_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_number) then 
        call read_particles_num_init_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_startpars('particles_num_init_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_selfgravity) then
        call read_particles_selfg_init_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_startpars('particles_selg_init_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_nbody) then
        call read_particles_nbody_init_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_startpars('particles_nbody_init_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_viscosity) then
        call read_particles_visc_init_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_startpars('particles_visc_init_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_stalker) then
        call read_pstalker_init_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_startpars('particles_stalker_init_pars',iostat)
            return
          endif
        endif
      endif
!
    endsubroutine particles_read_startpars
!***********************************************************************
    subroutine particles_rparam(unit)
!
!  Read particle parameters from start.in.
!
!  13-may-09/anders: coded
!
      integer, intent (in) :: unit
!
      call read_particles_init_pars(unit)
      if (lparticles_radius)      call read_particles_rad_init_pars(unit)
      if (lparticles_spin)        call read_particles_spin_init_pars(unit)
      if (lparticles_number)      call read_particles_num_init_pars(unit)
      if (lparticles_selfgravity) call read_particles_selfg_init_pars(unit)
      if (lparticles_nbody)       call read_particles_nbody_init_pars(unit)
      if (lparticles_viscosity)   call read_particles_visc_init_pars(unit)
      if (lparticles_stalker)     call read_pstalker_init_pars(unit)
!
    endsubroutine particles_rparam
!***********************************************************************
    subroutine samplepar_startpars(label,iostat)
!
!  Print sample of particle part of start.in.
!
!  17-aug-08/wlad: copied from param_io
!
      character (len=*), optional :: label
      integer, optional :: iostat
!
      if (lroot) then
        print*
        print*,'-----BEGIN sample particles namelist ------'
        if (lparticles) &
            print*,'&particles_init_pars         /'
        if (lparticles_radius) &
            print*,'&particles_radius_init_pars  /'
        if (lparticles_spin) &
            print*,'&particles_spin_init_pars    /'
        if (lparticles_number) &
            print*,'&particles_number_init_pars  /'
        if (lparticles_selfgravity) &
            print*,'&particles_selfgrav_init_pars/'
        if (lparticles_nbody) &
            print*,'&particles_nbody_init_pars   /'
        if (lparticles_viscosity) &
            print*,'&particles_visc_init_pars    /'
        if (lparticles_stalker) &
            print*,'&particles_stalker_init_pars/'
        print*,'------END sample particles namelist -------'
        print*
        if (present(label)) &
            print*, 'Found error in input namelist "' // trim(label)
        if (present(iostat)) print*, 'iostat = ', iostat
        if (present(iostat).or.present(label)) &
            print*,  '-- use sample above.'
      endif
!
    endsubroutine samplepar_startpars
!***********************************************************************
    subroutine particles_wparam(unit)
!
!  Write particle start parameters to file.
!
      integer, intent (in) :: unit
!
      call write_particles_init_pars(unit)
      if (lparticles_radius)      call write_particles_rad_init_pars(unit)
      if (lparticles_spin)        call write_particles_spin_init_pars(unit)
      if (lparticles_number)      call write_particles_num_init_pars(unit)
      if (lparticles_selfgravity) call write_particles_selfg_init_pars(unit)
      if (lparticles_nbody)       call write_particles_nbody_init_pars(unit)
      if (lparticles_viscosity)   call write_particles_visc_init_pars(unit)
      if (lparticles_stalker)     call write_pstalker_init_pars(unit)
!
    endsubroutine particles_wparam
!***********************************************************************
    subroutine particles_read_runpars(unit,iostat)
!
!  Read particle run parameters from run.in.
!
      integer, intent (in) :: unit
      integer, intent (inout), optional :: iostat
!
      call read_particles_run_pars(unit,iostat)
      if (present(iostat)) then
        if (iostat/=0) then
          call samplepar_runpars('particles_run_pars',iostat); return
        endif
      endif
!      
      if (lparticles_radius) then 
        call read_particles_rad_run_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_runpars('particles_rad_run_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_spin) then
        call read_particles_spin_run_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_runpars('particles_spin_run_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_number) then
        call read_particles_num_run_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_runpars('particles_num_run_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_selfgravity) then
        call read_particles_selfg_run_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_runpars('particles_selfg_run_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_nbody) then
        call read_particles_nbody_run_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_runpars('particles_nbody_run_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_viscosity) then
        call read_particles_visc_run_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_runpars('particles_visc_run_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_collisions) then
        call read_particles_coll_run_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_runpars('particles_coll_run_pars',iostat); return
          endif
        endif
      endif
!
      if (lparticles_stalker) then
        call read_pstalker_run_pars(unit,iostat)
        if (present(iostat)) then
          if (iostat/=0) then
            call samplepar_runpars('particles_stalker_run_pars',iostat); return
          endif
        endif
      endif
!
    endsubroutine particles_read_runpars
!***********************************************************************
    subroutine samplepar_runpars(label,iostat)
!
!  Print sample of particle part of run.in.
!
      character (len=*), optional :: label
      integer, optional :: iostat
!
      if (lroot) then
        print*
        print*,'-----BEGIN sample particle namelist ------'
        if (lparticles)             print*,'&particles_run_pars         /'
        if (lparticles_radius)      print*,'&particles_radius_run_pars  /'
        if (lparticles_spin)        print*,'&particles_spin_run_pars    /'
        if (lparticles_number)      print*,'&particles_number_run_pars  /'
        if (lparticles_selfgravity) print*,'&particles_selfgrav_run_pars/'
        if (lparticles_nbody)       print*,'&particles_nbody_run_pars   /'
        if (lparticles_viscosity)   print*,'&particles_visc_run_pars    /'
        if (lparticles_collisions)  print*,'&particles_coll_run_pars    /'
        if (lparticles_stalker)     print*,'&particles_stalker_run_pars/'
        print*,'------END sample particle namelist -------'
        print*
        if (present(label)) &
            print*, 'Found error in input namelist "' // trim(label)
        if (present(iostat)) print*, 'iostat = ', iostat
        if (present(iostat).or.present(label)) &
            print*,  '-- use sample above.'
      endif
!
    endsubroutine samplepar_runpars
!***********************************************************************
    subroutine particles_wparam2(unit)
!
!  Write particle run parameters to file.
!
      integer, intent (in) :: unit
!
      if (lparticles)             call write_particles_run_pars(unit)
      if (lparticles_radius)      call write_particles_rad_run_pars(unit)
      if (lparticles_spin)        call write_particles_spin_run_pars(unit)
      if (lparticles_number)      call write_particles_num_run_pars(unit)
      if (lparticles_selfgravity) call write_particles_selfg_run_pars(unit)
      if (lparticles_nbody)       call write_particles_nbody_run_pars(unit)
      if (lparticles_viscosity)   call write_particles_visc_run_pars(unit)
      if (lparticles_collisions)  call write_particles_coll_run_pars(unit)
      if (lparticles_stalker)     call write_pstalker_run_pars(unit)
!
    endsubroutine particles_wparam2
!***********************************************************************
    subroutine particles_powersnap(f)
!
!  Calculate power spectra of particle variables.
!
!  01-jan-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call powersnap_particles(f)
!
    endsubroutine particles_powersnap
!***********************************************************************
    subroutine get_slices_particles(f,slices)
!
!  Write slices for animation of Particle variables.
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Particle number density
!
        case ('np')
          slices%yz= f(slices%ix,m1:m2    ,n1:n2     ,inp)
          slices%xz= f(l1:l2    ,slices%iy,n1:n2     ,inp)
          slices%xy= f(l1:l2    ,m1:m2    ,slices%iz ,inp)
          slices%xy2=f(l1:l2    ,m1:m2    ,slices%iz2,inp)
          slices%ready = .true.
!
!  Particle mass density
!
        case ('rhop')
          if (irhop/=0) then
            slices%yz= f(slices%ix,m1:m2    ,n1:n2     ,irhop)
            slices%xz= f(l1:l2    ,slices%iy,n1:n2     ,irhop)
            slices%xy= f(l1:l2    ,m1:m2    ,slices%iz ,irhop)
            slices%xy2=f(l1:l2    ,m1:m2    ,slices%iz2,irhop)
            slices%ready = .true.
          else
            slices%yz= rhop_tilde * f(slices%ix,m1:m2    ,n1:n2     ,inp)
            slices%xz= rhop_tilde * f(l1:l2    ,slices%iy,n1:n2     ,inp)
            slices%xy= rhop_tilde * f(l1:l2    ,m1:m2    ,slices%iz ,inp)
            slices%xy2=rhop_tilde * f(l1:l2    ,m1:m2    ,slices%iz2,inp)
           slices%ready = .true.
          endif
!
      endselect
!
    endsubroutine get_slices_particles
!***********************************************************************
endmodule Particles_main
