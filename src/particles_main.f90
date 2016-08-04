! $Id$
!
!  This module contains all the main structure needed for particles.
!
module Particles_main
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles
  use Particles_adaptation
  use Particles_adsorbed
  use Particles_cdata
  use Particles_chemistry
  use Particles_coagulation
  use Particles_collisions
  use Particles_density
  use Particles_diagnos_dv
  use Particles_diagnos_state
  use Particles_drag
  use Particles_map
  use Particles_mass
  use Particles_mpicomm
  use Particles_nbody
  use Particles_number
  use Particles_radius
  use Particles_potential
  use Particles_grad
  use Particles_selfgravity
  use Particles_sink
  use Particles_spin
  use Particles_stalker
  use Particles_stirring
  use Particles_sub
  use Particles_surfspec
  use Particles_temperature
  use Particles_lyapunov
!
  implicit none
!
  include 'particles_main.h'
!
  real, dimension(mpar_loc,mparray) :: fp = 0.0
  real, dimension(mpar_loc,mpvar) :: dfp = 0.0
  integer, dimension(mpar_loc,3) :: ineargrid = 0
!
  contains
!***********************************************************************
    subroutine particles_register_modules()
!
!  Register particle modules.
!
!  07-jan-05/anders: coded
!
      use Special, only: register_particles_special
!
      integer :: ipvar
!
      call register_particles              ()
     call register_particles_potential    ()
      call register_particles_radius       ()
      call register_particles_grad         ()
      call register_particles_spin         ()
      call register_particles_number       ()
      call register_particles_density      ()
      call register_particles_selfgrav     ()
      call register_particles_nbody        ()
      call register_particles_sink         ()
      call register_particles_TT           ()
      call register_particles_mass         ()
      call register_particles_drag         ()
      call register_particles_chem         ()
      call register_particles_ads          ()
      call register_particles_surfspec     ()
      call register_pars_diagnos_state     ()
      call register_particles_special      (npvar)
!
!  Print summary of variable names.
!
      if (lroot) then
        open(3,file=trim(datadir)//'/pvarname.dat',status='replace')
        do ipvar=1,mparray
          write(3,"(i4,2x,a)") ipvar, pvarname(ipvar)
        enddo
        close(3)
      endif
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
      call rprint_particles              (lreset,LWRITE=lroot)
      call rprint_particles_radius       (lreset,LWRITE=lroot)
      call rprint_particles_grad         (lreset,LWRITE=lroot)
      call rprint_particles_lyapunov     (lreset,LWRITE=lroot)
      call rprint_particles_sink         (lreset,LWRITE=lroot)
      call rprint_particles_spin         (lreset,LWRITE=lroot)
      call rprint_particles_number       (lreset,LWRITE=lroot)
      call rprint_particles_density      (lreset,LWRITE=lroot)
      call rprint_particles_selfgrav     (lreset,LWRITE=lroot)
      call rprint_particles_nbody        (lreset,LWRITE=lroot)
      call rprint_particles_TT           (lreset,LWRITE=lroot)
      call rprint_particles_mass         (lreset,LWRITE=lroot)
      call rprint_particles_ads          (lreset,LWRITE=lroot)
      call rprint_particles_surf         (lreset,LWRITE=lroot)
      call rprint_particles_chem         (lreset,LWRITE=lroot)
      call rprint_particles_coagulation  (lreset,LWRITE=lroot)
!      call rprint_particles_potential    (lreset,LWRITE=lroot)
      call rprint_particles_collisions   (lreset,LWRITE=lroot)
      call rprint_particles_diagnos_dv   (lreset,LWRITE=lroot)
      call rprint_particles_diagnos_state(lreset,LWRITE=lroot)
      if (lroot) close(3)
!
    endsubroutine particles_rprint_list
!***********************************************************************
    subroutine particles_initialize_modules(f)
!
!  Initialize particle modules.
!
!  07-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lyinyang) &
        call fatal_error('particles_initialize_modules','Particles not implemented on Yin-Yang grid')
!
!  Check if there is enough total space allocated for particles.
!
      if (npar/mpar_loc>ncpus) then
        if (lroot) then
          print*, 'particles_initialize_modules: '// &
          'total number of particle slots available at the processors '// &
          'is smaller than the number of particles!'
          print*, 'particles_initialize_modules: npar/ncpus=', npar/ncpus
          print*, 'particles_initialize_modules: mpar_loc-ncpus*npar_mig=', &
              mpar_loc-ncpus*npar_mig
        endif
        call fatal_error('particles_initialize_modules','not enough space')
      endif
!
!  Set mass and number density of individual particles inside each
!  superparticle. The mass of a superparticle is defined through
!
!    *      mpmat : mass of constituent particles
!    *   np_swarm : number density of constituent particles
!    * rhop_swarm : mass density of superparticle
!
!  One can either input these quantities by hand or set the wanted number
!  density and mass density per grid cell (np_const or rhop_const).
!
      if (rhop_const /= 0.0) then
        rhop_swarm = rhop_const / (real(npar) / nwgrid)
        if (all(lequidist)) mp_swarm = rhop_swarm * dvol
        if (mpmat /= 0.0 .or. np_swarm /= 0.0) then
          if (lparticles_radius) &
              call fatal_error('particles_initialize_modules', 'may not set mpmat or np_swarm when setting rhop_const')
          if (mpmat /= 0.0 .and. np_swarm /= 0.0) &
              call fatal_error('particles_initialize_modules', 'must set only mpmat or np_swarm when using rhop_const')
          if (mpmat == 0.0) mpmat = rhop_swarm / np_swarm
          if (np_swarm == 0.0) np_swarm = rhop_swarm / mpmat
        endif
      elseif (np_const/=0.0) then
        if (lparticles_number) then
          if (lroot) print*, 'particles_initialize_modules: '// &
              'can not use np_const together with Particles_number module'
          call fatal_error('particles_initialize_modules','')
        endif
        if (.not.lparticles_radius) then
          if (mpmat==0.0) then
            if (lroot) print*, 'particles_initialize_modules: '// &
                'must have mpmat non zero when setting np_const'
            call fatal_error('particles_initialize_modules','')
          endif
        endif
        np_swarm=np_const/(real(npar)/nwgrid)
        rhop_swarm=np_swarm*mpmat
      else
        if (rhop_swarm==0.0) rhop_swarm=mpmat*np_swarm
      endif
!
      if (lparticles_radius .and. rhopmat>0.0 .and. &
        (np_swarm>0.0 .or. lparticles_number)) lignore_rhop_swarm=.true.

!
!  Initialize individual modules.
!
      call initialize_particles_mpicomm      (f)
      call initialize_particles              (f)
      call initialize_particles_map
      call initialize_particles_adaptation   (f)
      call initialize_particles_density      (f)
      call initialize_particles_nbody        (f)
      call initialize_particles_number       (f)
      call initialize_particles_radius       (f)
      call initialize_particles_grad         (f)
      call initialize_particles_selfgrav     (f)
      call initialize_particles_sink         (f)
      call initialize_particles_spin         (f)
      call initialize_particles_stalker      (f)
      call initialize_particles_TT           (f)
      call initialize_particles_mass         (f)
      call initialize_particles_drag
      call initialize_particles_ads          (f)
      call initialize_particles_surf         (f)
      call initialize_particles_coag         (f)
      call initialize_particles_collisions   (f)
      call initialize_pars_diagnos_state     (f)
      call initialize_particles_diagnos_dv   (f)
      call initialize_particles_potential    (f)
!
      if (lparticles_blocks) then
        if (lrun) then
          if (lroot) print*, 'particles_initialize_modules: reblocking particles'
          call map_nearest_grid(fp,ineargrid)
          call boundconds_particles(fp,ipar)
          call map_nearest_grid(fp,ineargrid)
          call sort_particles_iblock(fp,ineargrid,ipar)
          call map_xxp_grid(f,fp,ineargrid)
          call load_balance_particles(f,fp,ipar)
        else
          inearblock=0
        endif
      endif
!
!  Stop if rhop_swarm is zero.
!
      if (irhop/=0 .and. rhop_swarm==0.0 .and. (.not.lparticles_density) &
      .and. (.not.lignore_rhop_swarm)) then
        if (lroot) then
          print*, 'particles_initialize_modules: rhop_swarm is zero'
          print*, 'particles_initialize_modules: '// &
              'np_swarm, mpmat, rhop_swarm=', np_swarm, mpmat, rhop_swarm
        endif
        call fatal_error('particles_initialize_modules','rhop_swarm is zero')
      endif
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
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      if (lparticles_radius) call set_particle_radius(f,fp,1,npar_loc,init=.true.)
      if (lparticles_grad) call set_particle_grad(f,fp,1,npar_loc,init=.true.)
      if (lparticles_number)        call init_particles_number(f,fp)
      if (lparticles_density)       call init_particles_density(f,fp)
      call init_particles(f,fp,ineargrid)
      if (lparticles_nbody)         call init_particles_nbody(f,fp)
      if (lparticles_sink)          call init_particles_sink(f,fp)
      if (lparticles_spin)          call init_particles_spin(f,fp)
      if (lparticles_temperature)   call init_particles_TT(f,fp)
      if (lparticles_mass)          call init_particles_mass(f,fp)
      if (lparticles_drag)          call init_particles_drag(f,fp)
      if (lparticles_adsorbed)      call init_particles_ads(f,fp)
      if (lparticles_surfspec)      call init_particles_surf(f,fp,ineargrid)
      if (lparticles_diagnos_state) call init_particles_diagnos_state(fp)
      if (lparticles_lyapunov)      call init_particles_lyapunov(f,fp)
!
    endsubroutine particles_init
!***********************************************************************
    subroutine read_snapshot_particles(snap_directory)
!
      character (len=*) :: snap_directory
!
      call particles_read_snapshot(trim(snap_directory)//'/pvar.dat')
      if (lparticles_nbody) &
           call particles_nbody_read_snapshot(&
           trim(snap_directory)//'/spvar.dat')
!
    endsubroutine read_snapshot_particles
!***********************************************************************
    subroutine write_dim_particles(datadir)
!
      character (len=*) :: datadir
!
      call particles_write_pdim(trim(datadir)//'/pdim.dat')
      call particles_write_block(trim(datadir)//'/bdim.dat')
      if (lparticles_nbody) &
          call particles_nbody_write_spdim(trim(datadir)//'/spdim.dat')
!
    endsubroutine write_dim_particles
!***********************************************************************
    subroutine read_all_particles_init_pars()
!
      use Sub, only: read_namelist
!
      call read_namelist(read_particles_init_pars      ,'particles'         ,lparticles)
      call read_namelist(read_particles_rad_init_pars  ,'particles_radius'  ,lparticles_radius)
!     call read_namelist(read_particles_pot_init_pars  ,'particles_potential',lparticles_potential)
      call read_namelist(read_particles_spin_init_pars ,'particles_spin'    ,lparticles_spin)
      call read_namelist(read_particles_sink_init_pars ,'particles_sink'    ,lparticles_sink)
      call read_namelist(read_particles_num_init_pars  ,'particles_number'  ,lparticles_number)
      call read_namelist(read_particles_dens_init_pars ,'particles_dens'    ,lparticles_density)
      call read_namelist(read_particles_selfg_init_pars,'particles_selfgrav',lparticles_selfgravity)
      call read_namelist(read_particles_nbody_init_pars,'particles_nbody'   ,lparticles_nbody)
      call read_namelist(read_particles_mass_init_pars ,'particles_mass'    ,lparticles_mass)
      call read_namelist(read_particles_drag_init_pars ,'particles_drag'    ,lparticles_drag)
      call read_namelist(read_particles_TT_init_pars   ,'particles_TT'      ,lparticles_temperature)
      call read_namelist(read_particles_ads_init_pars  ,'particles_ads'     ,lparticles_adsorbed)
      call read_namelist(read_particles_surf_init_pars ,'particles_surf'    ,lparticles_surfspec)
      call read_namelist(read_particles_chem_init_pars ,'particles_chem'    ,lparticles_chemistry)
      call read_namelist(read_pstalker_init_pars       ,'particles_stalker' ,lparticles_stalker)
 call read_namelist(read_plyapunov_init_pars       ,'particles_lyapunov' ,lparticles_lyapunov)
!
    endsubroutine read_all_particles_init_pars
!***********************************************************************
    subroutine read_all_particles_run_pars()
!
      use Sub, only: read_namelist
!
      call read_namelist(read_particles_run_pars          ,'particles'              ,lparticles)
      call read_namelist(read_particles_adapt_run_pars    ,'particles_adapt'        ,lparticles_adaptation)
      call read_namelist(read_particles_rad_run_pars      ,'particles_radius'       ,lparticles_radius)
!      call read_namelist(read_particles_potential_run_pars,'particles_potential',lparticles_potential)
      call read_namelist(read_particles_spin_run_pars     ,'particles_spin'         ,lparticles_spin)
      call read_namelist(read_particles_sink_run_pars     ,'particles_sink'         ,lparticles_sink)
      call read_namelist(read_particles_num_run_pars      ,'particles_number'       ,lparticles_number)
      call read_namelist(read_particles_selfg_run_pars    ,'particles_selfgrav'     ,lparticles_selfgravity)
      call read_namelist(read_particles_nbody_run_pars    ,'particles_nbody'        ,lparticles_nbody)
      call read_namelist(read_particles_coag_run_pars     ,'particles_coag'         ,lparticles_coagulation)
      call read_namelist(read_particles_coll_run_pars     ,'particles_coll'         ,lparticles_collisions)
      call read_namelist(read_particles_stir_run_pars     ,'particles_stirring'     ,lparticles_stirring)
      call read_namelist(read_pstalker_run_pars           ,'particles_stalker'      ,lparticles_stalker)
      call read_namelist(read_pars_diagnos_dv_run_pars    ,'particles_diagnos_dv'   ,lparticles_diagnos_dv)
      call read_namelist(read_pars_diag_state_run_pars    ,'particles_diagnos_state',lparticles_diagnos_state)
      call read_namelist(read_particles_mass_run_pars     ,'particles_mass'         ,lparticles_mass)
      call read_namelist(read_particles_drag_run_pars     ,'particles_drag'         ,lparticles_drag)
      call read_namelist(read_particles_TT_run_pars       ,'particles_TT'           ,lparticles_temperature)
      call read_namelist(read_particles_ads_run_pars      ,'particles_ads'          ,lparticles_adsorbed)
      call read_namelist(read_particles_surf_run_pars     ,'particles_surf'         ,lparticles_surfspec)
      call read_namelist(read_particles_chem_run_pars     ,'particles_chem'         ,lparticles_chemistry)
      call read_namelist(read_plyapunov_run_pars     ,'particles_lyapunov'         ,lparticles_lyapunov)
!
    endsubroutine read_all_particles_run_pars
!***********************************************************************
    subroutine particles_read_snapshot(filename)
!
!  Read particle snapshot from file.
!
!  07-jan-05/anders: coded
!
      character (len=*) :: filename
!
      call input_particles(filename,fp,ipar)
!
    endsubroutine particles_read_snapshot
!***********************************************************************
    subroutine write_snapshot_particles(snap_directory,f,enum,snapnum)
!
      use General, only: itoa
!
      character (len=*) :: snap_directory
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: enum
      integer, optional :: snapnum
!
      if (present(snapnum)) then
        call particles_write_snapshot(trim(snap_directory)//'/PVAR'//itoa(snapnum),f, &
            enum=.false.,FLIST='pvarN.list')
!
        if (lparticles_nbody) call particles_nbody_write_snapshot(&
            trim(snap_directory)//'/SPVAR'//itoa(snapnum),enum=.false.)
!
      elseif (enum) then
        call particles_write_snapshot(trim(snap_directory)//'/PVAR',f, &
            ENUM=.true.,FLIST='pvarN.list')
!
        if (lparticles_nbody) call particles_nbody_write_snapshot(&
            trim(snap_directory)//'/SPVAR',enum=.true.,flist='spvarN.list')
      else
        call particles_write_snapshot( &
            trim(snap_directory)//'/pvar.dat',f,enum=.false.)
!
        if (lparticles_nbody.and.lroot) then
          call particles_nbody_write_snapshot( &
              trim(snap_directory)//'/spvar.dat',enum=.false.)
        endif
      endif
!
    endsubroutine write_snapshot_particles
!***********************************************************************
    subroutine particles_write_snapshot(chsnap,f,enum,flist)
!
!  Write particle snapshot to file.
!
!  07-jan-05/anders: coded
!
      character (len=*) :: chsnap
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: enum
      character (len=*), optional :: flist
!
      logical :: lsnap
!
      if (present(flist)) then
        call wsnap_particles(chsnap,f,fp,enum,lsnap,dsnap_par_minor,dsnap_par,ipar,mparray,flist)
      else
        call wsnap_particles(chsnap,f,fp,enum,lsnap,dsnap_par_minor,dsnap_par,ipar,mparray)
      endif
!
    endsubroutine particles_write_snapshot
!***********************************************************************
    subroutine particles_write_dsnapshot(chsnap,f)
!
!  Write particle derivative snapshot to file.
!
!  07-jan-05/anders: coded
!
      character (len=*) :: chsnap
      real, dimension (mx,my,mz,mfarray) :: f
!
      call wsnap_particles(chsnap,f,dfp,.false.,.false.,0.0,0.0,ipar,mpvar,nobound=.true.)
!
    endsubroutine particles_write_dsnapshot
!***********************************************************************
    subroutine particles_write_pdim(filename)
!
!  Write npar and mpvar to file.
!
!  09-jan-05/anders: coded
!  10-oct-14/jonas,nils: added the number of particle auxiliaries
!
      character (len=*) :: filename
!
      open(1,file=filename)
        write(1,'(4i9)') npar, mpvar, npar_stalk, mpaux
      close(1)
!
    endsubroutine particles_write_pdim
!***********************************************************************
    subroutine particles_write_block(filename)
!
!  Write block domain decomposition parameters to file.
!
!  05-nov-09/anders: coded
!
      character (len=*) :: filename
!
      if (lparticles_blocks) then
        open(1,file=filename)
          write(1,'(4i9)') nbrickx, nbricky, nbrickz, nblockmax
          write(1,'(4i9)') mxb, myb, mzb, nghostb
          write(1,'(3i9)') nxb, nyb, nzb
          write(1,'(6i9)') l1b, l2b, m1b, m2b, n1b, n2b
        close(1)
      endif
!
    endsubroutine particles_write_block
!***********************************************************************
    subroutine particles_timestep_first(f)
!
!  Setup dfp in the beginning of each itsub.
!
!  07-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
      if (lfirst) then
        dfp(1:npar_loc,:)=0.0
      else
        dfp(1:npar_loc,:)=alpha_ts(itsub)*dfp(1:npar_loc,:)
      endif
!
!  Insert new particles continuously during the run
!
      if (lfirst) call particles_insert_continuously(f)
!
    endsubroutine particles_timestep_first
!***********************************************************************
    subroutine particles_timestep_second(f)
!
!  Time evolution of particle variables.
!
!  07-jan-05/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (.not.(lparticles_nbody.and.(lcylindrical_coords.or.lspherical_coords))) then
        fp(1:npar_loc,1:mpvar) = fp(1:npar_loc,1:mpvar) + dt_beta_ts(itsub)*dfp(1:npar_loc,1:mpvar)
      else
        call advance_particles_in_cartesian(fp,dfp)
      endif
!
!  Discrete particle collisions. Must be done at the end of the time-step.
!   This call also sorts the particles into mn 
!
      call particles_discrete_collisions()
!
!  Adapt the number of particles in each grid cell to a desired number
!
      if (lparticles_adaptation .and. llast) then
        call sort_particles_imn(fp,ineargrid,ipar,f=f)
        call map_xxp_grid(f,fp,ineargrid)
        call map_vvp_grid(f,fp,ineargrid)
        call particles_adaptation_pencils(f,fp,dfp,ipar,ineargrid)
      endif
!
!  Insert "old" (removed) particles.
!
      call insert_particles_now(f)
!
    endsubroutine particles_timestep_second
!***********************************************************************
    subroutine split_update_particles(f, dt)
!
!  Wrapper for operator split terms for particle dynamics.
!
!  08-may-16/ccyang: coded.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: dt
!
      drag: if (lparticles_drag) then
        call particles_boundconds(f)
        call integrate_drag(f, fp, dt)
      endif drag
!
    endsubroutine split_update_particles
!***********************************************************************
    subroutine particles_discrete_collisions()
!
!  Discrete particle collisions.
!
!  13-nov-09/anders: coded
!
      if ( lparticles_stirring .and. llast ) then
        call particle_stirring(fp,ineargrid)
      endif
!
      if ( (lparticles_collisions.or.lparticles_coagulation) .and. llast ) then
!
        call boundconds_particles(fp,ipar)
        call map_nearest_grid(fp,ineargrid)
!
        if (lparticles_blocks) then
          call sort_particles_iblock(fp,ineargrid,ipar)
          if (lparticles_collisions) then
            call particles_collisions_blocks(fp,ineargrid)
          endif
          if (lparticles_coagulation) then
            call particles_coagulation_blocks(fp,ineargrid)
          endif
        else
          call sort_particles_imn(fp,ineargrid,ipar)
          if (lparticles_collisions) then
            call particles_collisions_pencils(fp,ineargrid)
          endif
          if (lparticles_coagulation) then
            call particles_coagulation_pencils(fp,ineargrid)
          endif
        endif
!
      endif
!
    endsubroutine particles_discrete_collisions
!***********************************************************************
    subroutine particles_load_balance(f)
!
!  Redistribute particles among the processors for better load balancing.
!
!  04-nov-09/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lparticles_blocks .and. mod(it,it1_loadbalance)==0) then
        call particles_boundconds(f)
        call load_balance_particles(f,fp,ipar)
      endif
!
    endsubroutine particles_load_balance
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
      call boundconds_particles(fp,ipar,dfp=dfp)
!
!  Remove particles that are too close to sink particles or sink points.
!  WARNING: ineargrid and the mapped particle density have not been updated
!  yet, and the sink particle subroutine must not rely on those arrays.
!
      if  (.not.lpencil_check_at_work) then
        if (lparticles) &
            call remove_particles_sink_simple(f,fp,dfp,ineargrid)
        if (lparticles_nbody) &
            call remove_particles_sink_nbody(f,fp,dfp,ineargrid)
        if (lparticles_sink)  call remove_particles_sink(f,fp,dfp,ineargrid)
      endif
!
!  Find nearest grid point for each particle.
!
      call map_nearest_grid(fp,ineargrid)
!
!  Sort particles so that they can be accessed contiguously in the memory.
!
      if (lparticles_blocks) then
        call sort_particles_iblock(fp,ineargrid,ipar,dfp=dfp)
      else
        call sort_particles_imn(fp,ineargrid,ipar,dfp=dfp,f=f)
      endif
!
!  Map the particle positions and velocities on the grid.
!
      call map_xxp_grid(f,fp,ineargrid)
      call map_vvp_grid(f,fp,ineargrid)
!
!  Distribute the n-body particles across processors
!
      if (lparticles_nbody) call bcast_nbodyarray(fp)
!
    endsubroutine particles_boundconds
!***********************************************************************
    subroutine particles_calc_selfpotential(f,rhs_poisson,lcontinued)
!
!  Calculate the potential of the dust particles (wrapper).
!
!  13-jun-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: rhs_poisson
      logical :: lcontinued
!
      call calc_selfpotential_particles(f,rhs_poisson,lcontinued)
      call calc_selfpot_sinkparticles(f,rhs_poisson,fp,ineargrid)
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
      if (lparticles) then
        call particles_dragforce_stiff(f,fp,ineargrid)
        call periodic_boundcond_on_aux(f)
      endif
!
!  Calculate the summed gravity due to the massive particles.
!  Needed because it is too slow to calculate the gravity at the
!  position of all dust particles. So we calculate the gravity
!  for the grid and interpolate to the position of the particles.
!
      if (lparticles_nbody)     call calc_nbodygravity_particles(f)
!
    endsubroutine particles_before_boundary
!***********************************************************************
    subroutine particles_special(f)
!
!  Fetch fp (and fsp) array to special module.
!
!  01-mar-08/wlad: coded
!
      use Special, only: special_calc_particles
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call special_calc_particles(f,fp,ineargrid)
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
      if (lparticles_density)     call pencil_criteria_par_density()
      if (lparticles_selfgravity) call pencil_criteria_par_selfgrav()
      if (lparticles_nbody)       call pencil_criteria_par_nbody()
      if (lparticles_temperature) call pencil_criteria_par_TT()
      if (lparticles_mass)        call pencil_criteria_par_mass()
      if (lparticles_adsorbed)    call pencil_criteria_par_ads()
      if (lparticles_chemistry)   call pencil_criteria_par_chem()
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
      if (lparticles_chemistry)   call calc_pencils_par_chem(f,p)
!
    endsubroutine particles_calc_pencils
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
      if (lfirst .and. (.not.lpencil_check_at_work)) &
          call particles_stalker_sub(f,fp,ineargrid)
!
!  Dynamical equations.
!
      if (lparticles)             call dxxp_dt(f,df,fp,dfp,ineargrid)
      if (lparticles)             call dvvp_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_potential)   call dvvp_dt_potential(f,df,fp,dfp,ineargrid)
      if (lparticles_lyapunov)    call dlyapunov_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_radius)      call dap_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_spin)        call dps_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_mass)        call dpmass_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_temperature) call dpTT_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_adsorbed)    call dpads_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_surfspec)    call dpsurf_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_number)      call dnpswarm_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_density)     call drhopswarm_dt(f,df,fp,dfp,ineargrid)
      if (lparticles_selfgravity) call dvvp_dt_selfgrav(f,df,fp,dfp,ineargrid)
      if (lparticles_nbody)       call dxxp_dt_nbody(dfp)
      if (lparticles_nbody)       call dvvp_dt_nbody(f,df,fp,dfp,ineargrid)
!
!  Create new sink particles or sink points.
!
      if  (.not.lpencil_check_at_work) then
        if (lparticles) &
            call create_particles_sink_simple(f,fp,dfp,ineargrid)
        if (lparticles_nbody) call create_particles_sink_nbody(f,fp,dfp,ineargrid)
        if (lparticles_sink)  call create_particles_sink(f,fp,dfp,ineargrid)
      endif
!
!  Correct for curvilinear geometry.
!
      if (.not.lparticles_nbody) call correct_curvilinear
!
!  Output particle size distribution to file.
!
      if (lparticles_radius .and. loutput_psize_dist .and. ldiagnos .and. &
          .not. lpencil_check_at_work) call output_particle_size_dist(fp)
!
    endsubroutine particles_pde
!***********************************************************************
    subroutine particles_pde_pencil(f,df,p)
!
!  Dynamical evolution of particle variables in pencils.
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
      call timing('particles_pde_pencil','entered',mnloop=.true.)
      if (lshepherd_neighbour) &
          call shepherd_neighbour_pencil(fp,ineargrid,kshepherd,kneighbour)
!
!  Interpolate required quantities using the predefined policies. Variables
!  are found in interp.
!  (Clean-up should be performed at end of this subroutine!)
!
     call interpolate_quantities(f,fp,p,ineargrid)
!
!  If reactive particles are enabled, needed quantities are calculated
!
     if (lparticles_chemistry) then
        call calc_pchemistry_pencils(f,fp,p,ineargrid)
     endif
     if (lparticles_surfspec) then
       call calc_psurf_pencils(f,fp,p,ineargrid)
     endif
       
!
!  Dynamical equations.
!
      if (lparticles)        call dxxp_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles)        then
!        write(*,*) 'DM: calling dvvp_dt_pencil'
        call dvvp_dt_pencil(f,df,fp,dfp,p,ineargrid)
      endif
      if (lparticles_lyapunov) call dlyapunov_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_radius)   call dap_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_spin)     call dps_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_mass)     call dpmass_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_temperature) call dpTT_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_adsorbed) call dpads_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_surfspec) call dpsurf_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_number) call dnpswarm_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_density)   call drhopswarm_dt_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_selfgravity) &
          call dvvp_dt_selfgrav_pencil(f,df,fp,dfp,p,ineargrid)
      if (lparticles_nbody) &
          call dvvp_dt_nbody_pencil(f,df,fp,dfp,p,ineargrid)
!      if (lparticles_potential) &
!          call dvvp_dt_potential_pencil(f,df,fp,dfp,ineargrid)
!      if (lparticles_polymer) &
!          call dRR_dt_pencil(f,df,fp,dfp,ineargrid)
!
!  Time-step contribution from discrete particle collisions.
!
      if (lparticles_collisions) &
          call particles_collisions_timestep(fp,ineargrid)
      if (lparticles_coagulation) &
          call particles_coagulation_timestep(fp,ineargrid)
!
      call cleanup_chemistry_pencils()
      call cleanup_surf_pencils()
      call cleanup_interpolated_quantities()
      call timing('particles_pde_pencil','finished',mnloop=.true.)
!
    endsubroutine particles_pde_pencil
!***********************************************************************
    subroutine particles_pde_blocks(f,df)
!
!  Dynamical evolution of particle variables in blocks.
!
!  29-nov-09/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      intent (inout) :: f, df
!
      if (lparticles_blocks) then
!
!  Fill adopted blocks with gas density and gas velocity field, for calculating
!  drag forces.
!
        if (lfill_blocks_density) then
          if (ldensity_nolog) then
            call fill_blocks_with_bricks(f,fb,mfarray,irho)
          else
            call fill_blocks_with_bricks(f,fb,mfarray,ilnrho)
          endif
        endif
!
        if (lfill_blocks_velocity) then
          call fill_blocks_with_bricks(f,fb,mfarray,iux)
          call fill_blocks_with_bricks(f,fb,mfarray,iuy)
          call fill_blocks_with_bricks(f,fb,mfarray,iuz)
        endif
!
!  Fill adopted blocks with gravitational acceleration.
!
        if (lfill_blocks_gpotself) then
          call fill_blocks_with_bricks(f,fb,mfarray,igpotselfx)
          call fill_blocks_with_bricks(f,fb,mfarray,igpotselfy)
          call fill_blocks_with_bricks(f,fb,mfarray,igpotselfz)
        endif
!
!  Zero the block contribution to time evolution of gas velocity.
!
        if (lhydro) dfb(:,:,:,iux:iuz,0:nblock_loc-1)=0.0
!
!  Dynamical equations.
!
        if (lparticles) call dxxp_dt_blocks(f,df,fp,dfp,ineargrid)
        if (lparticles) call dvvp_dt_blocks(f,df,fp,dfp,ineargrid)
!
!  Add block contribution to df to the main grid.
!
        if (lfill_bricks_velocity) then
          call fill_bricks_with_blocks(df,dfb,mvar,iux)
          call fill_bricks_with_blocks(df,dfb,mvar,iuy)
          call fill_bricks_with_blocks(df,dfb,mvar,iuz)
        endif
!
      endif
!
    endsubroutine particles_pde_blocks
!***********************************************************************
    subroutine correct_curvilinear
!
!  Curvilinear corrections to acceleration only.
!  Corrections to velocity were already taken into account
!  in the dxx_dt of particles_dust.f90
!
!  In the case that the N-body code is used, the update in polar grids
!  in done by transforming the variables first to Cartesian, to achieve a
!  better conservation of the Jacobi constant, and this code is not called.
!
!  15-sep-07/wlad: coded
!
      real :: rad,raddot,phidot,thtdot,sintht,costht
      integer :: k
!
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
    subroutine write_all_particles_init_pars(unit)
!
!  Write particle start parameters to file.
!
      integer, intent (in) :: unit
!
      call write_particles_init_pars(unit)
      if (lparticles_radius)      call write_particles_rad_init_pars(unit)
!      if (lparticles_potential)   call write_particles_pot_init_pars(unit)
      if (lparticles_spin)        call write_particles_spin_init_pars(unit)
      if (lparticles_sink)        call write_particles_sink_init_pars(unit)
      if (lparticles_number)      call write_particles_num_init_pars(unit)
      if (lparticles_density)     call write_particles_dens_init_pars(unit)
      if (lparticles_selfgravity) call write_particles_selfg_init_pars(unit)
      if (lparticles_nbody)       call write_particles_nbody_init_pars(unit)
      if (lparticles_stalker)     call write_pstalker_init_pars(unit)
      if (lparticles_mass)        call write_particles_mass_init_pars(unit)
      if (lparticles_drag)        call write_particles_drag_init_pars(unit)
      if (lparticles_temperature) call write_particles_TT_init_pars(unit)
      if (lparticles_adsorbed)    call write_particles_ads_init_pars(unit)
      if (lparticles_surfspec)    call write_particles_surf_init_pars(unit)
      if (lparticles_chemistry)   call write_particles_chem_init_pars(unit)
      if (lparticles_lyapunov)   call write_plyapunov_init_pars(unit)
!
    endsubroutine write_all_particles_init_pars
!***********************************************************************
    subroutine write_all_particles_run_pars(unit)
!
!  Write particle run parameters to file.
!
      integer, intent (in) :: unit
!
      if (lparticles)                call write_particles_run_pars(unit)
      if (lparticles_radius)         call write_particles_rad_run_pars(unit)
!      if (lparticles_potential)      call write_particles_pot_run_pars(unit)
      if (lparticles_spin)           call write_particles_spin_run_pars(unit)
      if (lparticles_sink)           call write_particles_sink_run_pars(unit)
      if (lparticles_number)         call write_particles_num_run_pars(unit)
      if (lparticles_selfgravity)    call write_particles_selfg_run_pars(unit)
      if (lparticles_nbody)          call write_particles_nbody_run_pars(unit)
      if (lparticles_coagulation)    call write_particles_coag_run_pars(unit)
      if (lparticles_collisions)     call write_particles_coll_run_pars(unit)
      if (lparticles_stirring)       call write_particles_stir_run_pars(unit)
      if (lparticles_stalker)        call write_pstalker_run_pars(unit)
      if (lparticles_diagnos_dv)     call write_pars_diagnos_dv_run_pars(unit)
      if (lparticles_diagnos_state)  call write_pars_diag_state_run_pars(unit)
      if (lparticles_mass)           call write_particles_mass_run_pars(unit)
      if (lparticles_drag)           call write_particles_drag_run_pars(unit)
      if (lparticles_temperature)    call write_particles_TT_run_pars(unit)
      if (lparticles_adsorbed)       call write_particles_ads_run_pars(unit)
      if (lparticles_surfspec)       call write_particles_surf_run_pars(unit)
      if (lparticles_chemistry)      call write_particles_chem_run_pars(unit)
      if (lparticles_lyapunov)       call write_plyapunov_run_pars(unit)
!
    endsubroutine write_all_particles_run_pars
!***********************************************************************
    subroutine wsnap_particles(snapbase,f,fp,enum,lsnap,dsnap_par_minor,dsnap_par,ipar,varsize,flist,nobound)
!
!  Write particle snapshot file, labelled consecutively if enum==.true.
!  Otherwise just write a snapshot without label (used e.g. for pvar.dat)
!
!  29-dec-04/anders: adapted from wsnap
!  04-oct-08/ccyang: use a separate log file for minor snapshots
!  26-nov-08/ccyang: add independent sequence for particle snapshots
!
      use General, only: itoa, safe_character_assign
      use IO, only: log_filename_to_file
      use Particles_mpicomm, only: output_blocks
      use Sub
!
      integer :: varsize
      character (len=*) :: snapbase, flist
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,varsize) :: fp
      logical :: enum, lsnap, nobound
      real :: dsnap_par_minor, dsnap_par
      integer, dimension (mpar_loc) :: ipar
!
      logical, save :: lfirst_call=.true.
      integer, save :: nsnap, nsnap_minor, nsnap_par
      real, save :: tsnap, tsnap_minor, tsnap_par
      character (len=fnlen), save :: fmajor, fminor, fpar
      logical :: lsnap_minor=.false., lsnap_par=.false.
      character (len=fnlen) :: snapname
      character (len=intlen) :: nsnap_ch,nsnap_minor_ch,nsnap_par_ch
!
      optional :: flist, nobound
!
!  Output snapshot with label in 'tsnap' time intervals
!  file keeps the information about number and time of last snapshot
!
      if (enum) then
!
!  At first call, need to initialize tsnap.
!  tsnap calculated in read_snaptime, but only available to root processor.
!
        if (lfirst_call) then
          call safe_character_assign(fmajor,trim(datadir)//'/tsnap.dat')
          call read_snaptime(fmajor,tsnap,nsnap,dsnap,t)
          if (dsnap_par_minor>0.0) then
            call safe_character_assign(fminor,trim(datadir)//'/tsnap_minor.dat')
            call read_snaptime(fminor,tsnap_minor,nsnap_minor,dsnap_par_minor,t)
          endif
          if (dsnap_par>0.0) then
            call safe_character_assign(fpar,trim(datadir)//'/tsnap_par.dat')
            call read_snaptime(fpar,tsnap_par,nsnap_par,dsnap_par,t)
          endif
          lfirst_call=.false.
        endif
!
!  Output independent sequence of particle snapshots.
!
        if (dsnap_par>0.0) then
          call update_snaptime(fpar,tsnap_par,nsnap_par,dsnap_par,t,lsnap_par,nsnap_par_ch)
          if (lsnap_par) then
            snapname=trim(snapbase)//'_'//nsnap_par_ch
            call particles_boundconds(f)
            call output_particles(snapname,fp,ipar)
            if (ip<=10 .and. lroot) &
                print*,'wsnap_particles: written snapshot ', snapname
            if (present(flist)) call log_filename_to_file(snapname,flist)
          endif
        endif
!
!  Possible to output minor particle snapshots (e.g. for a movie).
!
        if (dsnap_par_minor>0.0) then
          call update_snaptime(fminor,tsnap_minor,nsnap_minor,dsnap_par_minor,t,lsnap_minor,nsnap_minor_ch)
          if (lsnap_minor) then
            snapname=trim(snapbase)//trim(itoa(nsnap-1))//'.'//nsnap_minor_ch
            call particles_boundconds(f)
            call output_particles(snapname,fp,ipar)
            if (ip<=10 .and. lroot) &
                print*,'wsnap_particles: written snapshot ', snapname
            if (present(flist)) call log_filename_to_file(snapname,flist)
          endif
        endif
!
!  Regular data snapshots must come synchronized with the fluid snapshots.
!
        call update_snaptime(fmajor,tsnap,nsnap,dsnap,t,lsnap,nsnap_ch,nowrite=.true.)
        if (lsnap) then
          snapname=trim(snapbase)//nsnap_ch
          call particles_boundconds(f)
          call output_particles(snapname,fp,ipar)
          if (lparticles_blocks) &
              call output_blocks(trim(directory_dist)//'/BLOCKS'//nsnap_ch)
          if (ip<=10 .and. lroot) &
              print*,'wsnap_particles: written snapshot ', snapname
          if (present(flist)) call log_filename_to_file(snapname,flist)
          nsnap_minor=1
        endif
!
      else
!
!  Write snapshot without label
!
        snapname=snapbase
        if (present(nobound)) then
          if (.not. nobound) call particles_boundconds(f)
        else
          call particles_boundconds(f)
        endif
        call output_particles(snapname,fp,ipar)
        if (lparticles_blocks) &
            call output_blocks(trim(directory_dist)//'/blocks.dat')
        if (ip<=10 .and. lroot) &
             print*,'wsnap_particles: written snapshot ', snapname
        if (present(flist)) call log_filename_to_file(snapname,flist)
      endif
!
    endsubroutine wsnap_particles
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
      integer :: l
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Particle number density
!
        case ('np')
          slices%yz= f(ix_loc,m1:m2 ,n1:n2  ,inp)
          slices%xz= f(l1:l2 ,iy_loc,n1:n2  ,inp)
          slices%xy= f(l1:l2 ,m1:m2 ,iz_loc ,inp)
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,inp)
          slices%ready = .true.
!
!  Particle mass density
!
        case ('rhop')
          if (irhop/=0) then
            slices%yz= f(ix_loc,m1:m2 ,n1:n2  ,irhop)
            slices%xz= f(l1:l2 ,iy_loc,n1:n2  ,irhop)
            slices%xy= f(l1:l2 ,m1:m2 ,iz_loc ,irhop)
            slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,irhop)
            slices%ready = .true.
          else
            if (lcartesian_coords.and.(all(lequidist))) then
              slices%yz= rhop_swarm*f(ix_loc,m1:m2 ,n1:n2   ,inp)
              slices%xz= rhop_swarm*f(l1:l2 ,iy_loc,n1:n2  ,inp)
              slices%xy= rhop_swarm*f(l1:l2 ,m1:m2 ,iz_loc ,inp)
              slices%xy2=rhop_swarm*f(l1:l2 ,m1:m2 ,iz2_loc,inp)
            else
              do m=m1,m2 ; do n=n1,n2
                slices%yz(m,n) = rhop_swarm*f(ix_loc,m,n,inp)
              enddo;enddo
              do l=l1,l2 ; do n=n1,n2
                slices%xz(l,n) = rhop_swarm*f(l,iy_loc,n,inp)
              enddo;enddo
              do l=l1,l2 ; do m=m1,m2
                slices%xy(l,m) = rhop_swarm*f(l,m,iz_loc,inp)
                slices%xy2(l,m) = rhop_swarm*f(l,m,iz2_loc,inp)
              enddo;enddo
            endif
            slices%ready = .true.
          endif
!
      endselect
!
    endsubroutine get_slices_particles
!***********************************************************************
    subroutine insert_particles_now(f)
!
!  Insert particles which has been removed from the system as the
!  system has evolved.
!
!  2012-oct-19/dhruba: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (linsert_particle) call insert_lost_particles(f,fp,ineargrid)
!
    endsubroutine insert_particles_now
!***********************************************************************
    subroutine particles_stochastic()
      if (lparticles_lyapunov) call particles_stochastic_lyapunov(fp)
    endsubroutine particles_stochastic
!***********************************************************************
    subroutine particles_insert_continuously(f)
!
!  Insert particles continuously, i.e. add particles in
!  the beginning of a time step.
!
!  sep-09/kragset: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (linsert_particles_continuously) then
        call insert_particles(f,fp,ineargrid)
      endif
!
    endsubroutine particles_insert_continuously
!***********************************************************************
    subroutine particles_cleanup
!
!      call particles_final_clean_up()
      if (lparticles_chemistry) then 
        call particles_chemistry_clean_up()
        call particles_surfspec_clean_up()
        call particles_adsorbed_clean_up()
      endif
!
    endsubroutine particles_cleanup
!***********************************************************************
endmodule Particles_main
