! $Id$
!
!***********************************************************************
!
!  The Pencil Code is a high-order finite-difference code for compressible
!  hydrodynamic flows with magnetic fields and particles. It is highly
!  modular and can easily be adapted to different types of problems.
!
!      MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!      MMMMMMM7MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!      MMMMMMMIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!      MMMMMMMIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!      MMMMMMM7IMMMMMMMMMMMMMMMMMMM7IMMMMMMMMMMMMMMMMMMMMDMMMMMMM
!      MMMMMMZIIMMMMMMMMMMMMMMMMMMMIIMMMMMMMMMMMMMMMMMMMMIMMMMMMM
!      MMMMMMIIIZMMMMMMMMMMMMMMMMMMIIMMMMMMMMMMMMMMMMMMMMIMMMMMMM
!      MMMMMMIIIIMMMMMMMMMMMMMMMMMNII$MMMMMMMMMMMMMMMMMM$IMMMMMMM
!      MMMMM8IIIIMMMMMMMMMMMMMMMMM$IIIMMMMMMMMMMMMMMMMMMII7MMMMMM
!      MMMMMD7II7MMMMMMMMMMMMMMMMMIIIIMMMMMMMMMMMMMMMMMMIIIMMMMMM
!      MMMMN..:~=ZMMMMMMMMMMMMMMMMIIIIDMMMMMMMMMMMMMMMMDIIIMMMMMM
!      MMMM8.,:~=?MMMMMMMMMMMMMMMMOII7NMMMMMMMMMMMMMMMMZIIIDMMMMM
!      MMMM. ,::=+MMMMMMMMMMMMMMM..,~=?MMMMMMMMMMMMMMMMIIII$MMMMM
!      MMMM..,:~=+MMMMMMMMMMMMMM8 .,~=+DMMMMMMMMMMMMMMM8II78MMMMM
!      MMMM .,:~==?MMMMMMMMMMMMMN .,~=+NMMMMMMMMMMMMMM8 ,~~+MMMMM
!      MMM7  ,:~+=?MMMMMMMMMMMMM  .,~==?MMMMMMMMMMMMMM..,~~+DMMMM
!      MMM.  ,:~==?MMMMMMMMMMMMM  .,~~=?MMMMMMMMMMMMMM. ,~~+?MMMM
!      MMM.  ,:~~=??MMMMMMMMMMMM  ,,~~=?MMMMMMMMMMMMM8 .,~~=?NMMM
!      MMM.  ,:~~=+?MMMMMMMMMMMI  ,,~~=+?MMMMMMMMMMMM. .,~~=?MMMM
!      MM~. .,:~:==?MMMMMMMMMMM   ,,~~==?MMMMMMMMMMMM  .,~~=+?MMM
!      MMN8 .D,D+=M=8MMMMMMMMMM   ,,~~==?MMMMMMMMMMMM. .,~~==?MMM
!      MM==8.   .8===MMMMMMMMMM  .,,~~==?$MMMMMMMMMM=  .,~~==?MMM
!      MM==D.    +===MMMMMMMMMO$ .I?7$=7=NMMMMMMMMMM.  ,,~~==?$MM
!      MM==D.    +===MMMMMMMMM==O?..  ====MMMMMMMMMM.  ,,~~==??MM
!      MM==D.    +===MMMMMMMMM===.    .===MMMMMMMMMM+$.?=,7==?7MM
!      MM==D.    +===MMMMMMMMM===.    .===MMMMMMMMMZ==8    .8==MM
!      MM==D.    +===MMMMMMMMM===.    .===MMMMMMMMMZ==I    .Z==MM
!      MM==D. .  +===MMMMMMMMM===.... .===MMMMMMMMMZ==I    .Z==MM
!      MM==D.    +===MMMMMMMMM===.    .===MMMMMMMMMZ==I    .Z==MM
!      MM==D.    +===MMMMMMMMM===.    .===MMMMMMMMMZ==I    .Z==MM
!      MM==D.    +===MMMMMMMMM===..   .===MMMMMMMMMZ==I    .Z==MM
!      MM==D. .  +===MMMMMMMMM===.... .===MMMMMMMMMZ==I    .Z==MM
!
!  More information can be found in the Pencil Code manual and at the
!  website http://www.nordita.org/software/pencil-code/.
!
!***********************************************************************
program run
!
  use Boundcond,       only: update_ghosts
  use Cdata
  use Chemistry,       only: chemistry_clean_up, write_net_reaction, lchemistry_diag
  use Diagnostics
  use Dustdensity,     only: init_nd
  use Dustvelocity,    only: init_uud
  use Equ,             only: debug_imn_arrays,initialize_pencils
  use EquationOfState, only: ioninit,ioncalc
  use FArrayManager,   only: farray_clean_up
  use Filter
  use Forcing,         only: forcing_clean_up,addforce
  use Hydro,           only: hydro_clean_up,kinematic_random_phase
  use ImplicitPhysics, only: calc_heatcond_ADI
  use Interstellar,    only: check_SN,addmassflux
  use IO
  use Magnetic,        only: rescaling_magnetic
  use Messages
  use Mpicomm
  use NSCBC,           only: NSCBC_clean_up
  use Param_IO
  use Particles_main
  use Pencil_check,    only: pencil_consistency_check
  use Register
  use SharedVariables, only: sharedvars_clean_up
  use Signal_handling, only: signal_prepare, emergency_stop
  use Slices
  use Snapshot
  use Solid_Cells,     only: solid_cells_clean_up
  use Sub
  use Testscalar,      only: rescaling_testscalar
  use Testfield,       only: rescaling_testfield
  use TestPerturb,     only: testperturb_begin, testperturb_finalize
  use Timeavg
  use Timestep,        only: rk_2n
!
  implicit none
!
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (mx,my,mz,mvar) :: df
  type (pencil_case) :: p
  double precision :: time1, time2
  double precision :: time_last_diagnostic, time_this_diagnostic
  real :: wall_clock_time=0.0, time_per_step=0.0
  integer :: icount, i, mvar_in
  integer :: it_last_diagnostic, it_this_diagnostic
  logical :: lstop=.false., timeover=.false., resubmit=.false.
  logical :: suppress_pencil_check=.false.
  logical :: lreload_file=.false., lreload_always_file=.false.
!
  lrun=.true.
!
!  Get processor numbers and define whether we are root.
!
  call mpicomm_init
!
!  Identify version.
!
  if (lroot) call svn_id( &
      '$Id$')
!
!  Initialize the message subsystem, eg. color setting etc.
!
  call initialize_messages()
!
!  Read parameters from start.x (default values; may be overwritten by
!  read_runpars).
!
  call rparam()
!
!  Read parameters and output parameter list.
!
  call read_runpars()
!
!  Initialise MPI communication.
!
  call initialize_mpicomm()
!
!  Derived parameters (that may still be overwritten).
!  [might better be put into another routine, possibly even in rparam or
!  read_runpars]
!
  x0 = xyz0(1) ; y0 = xyz0(2) ; z0 = xyz0(3)
  Lx = Lxyz(1) ; Ly = Lxyz(2) ; Lz = Lxyz(3)
!
!  Size of box at local processor.
!
  Lxyz_loc(1)=Lxyz(1)/nprocx
  Lxyz_loc(2)=Lxyz(2)/nprocy
  Lxyz_loc(3)=Lxyz(3)/nprocz
  xyz0_loc(1)=xyz0(1)+ipx*Lxyz_loc(1)
  xyz0_loc(2)=xyz0(2)+ipy*Lxyz_loc(2)
  xyz0_loc(3)=xyz0(3)+ipz*Lxyz_loc(3)
  xyz1_loc(1)=xyz0_loc(1)+Lxyz_loc(1)
  xyz1_loc(2)=xyz0_loc(2)+Lxyz_loc(2)
  xyz1_loc(3)=xyz0_loc(3)+Lxyz_loc(3)
!
!  Register physics modules.
!
  call register_modules()
  if (lparticles) call particles_register_modules()
!
!  Inform about verbose level.
!
  if (lroot) print*, 'The verbose level is ip=', ip, ' (ldebug=', ldebug, ')'
!
!  Call rprint_list to initialize diagnostics and write indices to file.
!
  call rprint_list(LRESET=.false.)
  if (lparticles) call particles_rprint_list(.false.)
!
!  Populate wavenumber arrays for fft and calculate Nyquist wavenumber.
!
  if (nxgrid/=1) then
    kx_fft=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*2*pi/Lx
    kx_fft2=kx_fft**2
    kx_ny =nxgrid/2 * 2*pi/Lx
  else
    kx_fft=0.0
    kx_ny =0.0
  endif
!
  if (nygrid/=1) then
    ky_fft=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*2*pi/Ly
    ky_fft2=ky_fft**2
    ky_ny =nygrid/2 * 2*pi/Ly
  else
    ky_fft=0.0
    ky_ny =0.0
  endif
!
  if (nzgrid/=1) then
    kz_fft=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*2*pi/Lz
    kz_fft2=kz_fft**2
    kz_ny =nzgrid/2 * 2*pi/Lz
  else
    kz_fft=0.0
    kz_ny =0.0
  endif
!
!  Position of equator (if any).
!
  if (lequatory) yequator=xyz0(2)+0.5*Lxyz(2)
  if (lequatorz) zequator=xyz0(3)+0.5*Lxyz(3)
!
!  Limits to xaveraging.
!
  if (lav_smallx) call init_xaver
!
!  Inner radius for freezing variables defaults to r_min.
!  Note: currently (July 2005), hydro.f90 uses a different approach:
!  r_int will override rdampint, which doesn't seem to make much sense (if
!  you want rdampint to be overridden, then don't specify it in the first
!  place).
!
  if (rfreeze_int==-impossible .and. r_int>epsi) rfreeze_int=r_int
  if (rfreeze_ext==-impossible) rfreeze_ext=r_ext
!
!  Will we write all slots of f?
!
  if (lwrite_aux) then
    mvar_io=mvar+maux
  else
    mvar_io=mvar
  endif
!
! Shall we read also auxiliary variables or fewer variables (ex: turbulence
! field with 0 species as an input file for a chemistry problem)?
!
  if (lread_aux) then
    mvar_in=mvar+maux
  else if (lread_less) then
    mvar_in=4
  else
    mvar_in=mvar
  endif
!
!  Print resolution and dimension of the simulation.
!
  dimensionality=min(nxgrid-1,1)+min(nygrid-1,1)+min(nzgrid-1,1)
  if (lroot) write(*,'(a,i1,a)') ' This is a ', dimensionality, '-D run'
  if (lroot) print*, 'nxgrid, nygrid, nzgrid=', nxgrid, nygrid, nzgrid
  if (lroot) print*, 'Lx, Ly, Lz=', Lxyz
  if (lroot) print*, '      Vbox=', Lxyz(1)*Lxyz(2)*Lxyz(3)
!
!  Check if we want to divide cdtv by dimensionality.
!  (old_cdtv defaults to .false.)
!  [This is obsolete now that we calculate the time step in a different
!   manner -- could somebody please adjust visc_var and remove cdtvDim?]
!
  if (old_cdtv) then
    cdtvDim=cdtv
  else
    cdtvDim=cdtv/max(dimensionality,1)
  endif
!
!  Set up directory names `directory' and `directory_snap'.
!
  call directory_names()
!
!  Get state length of random number generator.
!
  call get_nseed(nseed)
!
!  Write particle block dimensions to file (may have been changed for better
!  load balancing).
!
  if (lroot) then
    if (lparticles) call particles_write_block(trim(datadir)//'/bdim.dat')
  endif
!
!  Read data.
!  Snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!  NOTE: for io_dist, rtime doesn't read the time, only for io_mpio.
!
  call rsnap(trim(directory_snap)//'/var.dat',f,mvar_in)
!
  if (lparticles) call read_snapshot_particles(directory_snap)
!
  call get_nseed(nseed)
!
!  Read time and global variables (if any).
!
  call rtime(trim(directory)//'/time.dat',t)
  if (mglobal/=0)  &
      call input_globals(trim(directory_snap)//'/global.dat', &
      f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
!
!  Set initial time to zero if requested.
!
  if (lini_t_eq_zero) t=0.0
!
!  Read coordinates.
!
  if (ip<=6.and.lroot) print*, 'reading grid coordinates'
  call rgrid(trim(directory)//'/grid.dat')
!
!  Read processor boundaries.
!
  if (lparticles) then
    if (ip<=6.and.lroot) print*, 'reading processor boundaries'
    call rproc_bounds(trim(directory)//'/proc_bounds.dat')
  endif
!
!  The following is here to avoid division in sub.f90 for diagnostic
!  outputs of integrated values in the non equidistant case.
!  Do this even for uniform meshes, in which case xprim=dx, etc.
!  Remember that dx_1=0 for runs without extent in that direction.
!
  if (nxgrid==1) then; xprim=1.0; else; xprim=1/dx_1; endif
  if (nygrid==1) then; yprim=1.0; else; yprim=1/dy_1; endif
  if (nzgrid==1) then; zprim=1.0; else; zprim=1/dz_1; endif
!
!  Determine slice positions and whether slices are to be written on this
!  processor. This can only be done after the grid has been established.
!
  call setup_slices()
!
!  Write parameters to log file (done after reading var.dat, since we
!  want to output time t.
!
  call print_runpars(FILE=trim(datadir)//'/params.log',ANNOTATION='Running')
!
!  Allow modules to do any physics modules do parameter dependent
!  initialization. And final pre-timestepping setup.
!  (must be done before need_XXXX can be used, for example)
!
  call initialize_modules(f,LSTARTING=.false.)
!
!  Initialize ionization array.
!
  if (leos_ionization) call ioninit(f)
  if (leos_temperature_ionization) call ioncalc(f)
!
!  Prepare particles.
!
  if (lparticles) then
    call particles_rprint_list(.false.)
    call particles_initialize_modules(f,lstarting=.false.)
  endif
!
!  Write data to file for IDL.
!
  call wparam2()
!
!  Possible debug output (can only be done after "directory" is set).
!  Check whether mn array is correct.
!
  if (ip<=3) call debug_imn_arrays
!
!  Find out which pencils are needed and write information about required,
!  requested and diagnostic pencils to disc.
!
  call choose_pencils()
  call write_pencil_info()
!
  if (mglobal/=0)  &
      call output_globals(trim(directory_snap)//'/global.dat', &
      f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
!
!  Update ghost zones, so rprint works corrected for at the first
!  time step even if we didn't read ghost zones.
!
  call update_ghosts(f)
!
!  Save spectrum snapshot.
!
  if (dspec/=impossible) call powersnap(f)
!
!  Save soundfile
!
!  lout_sound = (t-tsnap) >= dsound
!
!  Initialize pencils in the pencil_case.
!
  if (lpencil_init) call initialize_pencils(p,0.0)
!
!  Perform pencil_case consistency check if requested.
!
  suppress_pencil_check = control_file_exists("NO-PENCIL-CHECK")
  if ( (lpencil_check .and. .not. suppress_pencil_check) .or. &
       ((.not.lpencil_check).and.lpencil_check_small) ) then
    call pencil_consistency_check(f,df,p)
  endif
!
!  Start timing for final timing statistics.
!  Initialize timestep diagnostics during the run (whether used or not,
!  see idiag_timeperstep).
!
  if (lroot) then
    time1=mpiwtime()
    time_last_diagnostic=time1
    icount=0
    it_last_diagnostic=icount
  endif
!
  if (it1d==impossible_int) then
    it1d=it1
  else
    if (it1d<it1) then
      if (lroot) call stop_it_if_any(.true.,'run: it1d smaller than it1')
    endif
  endif
  ! globally catch eventual 'stop_it_if_any' call from single MPI ranks
  call stop_it_if_any(.false.,'')
!
!  Prepare signal catching
!
  call signal_prepare()
!
!  Do loop in time.
!
  Time_loop: do while (it<=nt)
!
    lout  =mod(it-1,it1) ==0
    l1davg=mod(it-1,it1d)==0
!
    if (lout .or. emergency_stop) then
!
!  Exit do loop if file `STOP' exists.
!
      lstop=control_file_exists('STOP',DELETE=.true.)
      if (lstop .or. t>tmax .or. emergency_stop) then
        if (lroot) then
          print*
          if (emergency_stop) print*, 'Emergency stop requested'
          if (lstop) print*, 'Found STOP file'
          if (t>tmax) print*, 'Maximum simulation time exceeded'
        endif
        resubmit=control_file_exists('RESUBMIT',DELETE=.true.)
        if (resubmit) print*, 'Cannot be resubmitted'
        exit Time_loop
      endif
!
!  initialize timer
!
      call timing('run','entered Time_loop',INSTRUCT='initialize')
!
!  Re-read parameters if file `RELOAD' exists; then remove the file
!  (this allows us to change parameters on the fly).
!  Re-read parameters if file `RELOAD_ALWAYS' exists; don't remove file
!  (only useful for debugging RELOAD issues).
!
      lreload_file       =control_file_exists('RELOAD')
      lreload_always_file=control_file_exists('RELOAD_ALWAYS')
      lreloading         =lreload_file .or. lreload_always_file
!
      if (lreloading) then
        if (lroot) write(*,*) 'Found RELOAD file -- reloading parameters'
!  Re-read configuration
        dt=0.0
        call read_runpars(PRINT=.true.,FILE=.true.,ANNOTATION='Reloading')
!
!  Before reading the rprint_list deallocate the arrays allocated for
!  1-D and 2-D diagnostics.
!
                                 call xyaverages_clean_up()
                                 call xzaverages_clean_up()
                                 call yzaverages_clean_up()
        if (lwrite_phizaverages) call phizaverages_clean_up()
        if (lwrite_yaverages)    call yaverages_clean_up()
        if (lwrite_zaverages)    call zaverages_clean_up()
        if (lwrite_phiaverages)  call phiaverages_clean_up()
        if (lforcing)            call forcing_clean_up()
        if (lhydro_kinematic)    call hydro_clean_up()
        if (lsolid_cells)        call solid_cells_clean_up()
        call rprint_list(LRESET=.true.) !(Re-read output list)
        call initialize_modules(f,LSTARTING=.false.)
        if (lparticles) then
          call particles_rprint_list(.false.)
          call particles_initialize_modules(f,lstarting=.false.)
        endif
        call choose_pencils()
        call wparam2()
!
        lreload_file=control_file_exists('RELOAD', DELETE=.true.)
        lreload_file        = .false.
        lreload_always_file = .false.
        lreloading          = .false.
      endif
    endif
!
! Insert particles (if linsert_particles_continuously==T)
!
      call particles_insert_continuously(f)
!
!  Remove wiggles in lnrho in sporadic time intervals.
!  Necessary on moderate-sized grids. When this happens,
!  this is often an indication of bad boundary conditions!
!  iwig=500 is a typical value. For iwig=0 no action is taken.
!  (The two queries below must come separately on compaq machines.)
!
!  N.B.: We now (July 2003) know that boundary conditions
!  change practically nothing, apart from avoiding stationary
!  stagnation points/surfaces in the first place.
!    rmwig is superseeded by the switches lupw_lnrho, lupw_ss,
!  which should provide a far better solution to the problem.
!
    if (iwig/=0) then
      if (mod(it,iwig)==0) then
        if (lrmwig_xyaverage) call rmwig_xyaverage(f,ilnrho)
        if (lrmwig_full) call rmwig(f,df,ilnrho,ilnrho,awig)
        if (lrmwig_rho) call rmwig(f,df,ilnrho,ilnrho,awig,explog=.true.)
      endif
    endif
!
!  If we want to write out video data, wvid sets lvideo=.true.
!  This allows pde to prepare some of the data.
!
    if (lwrite_slices) call wvid_prepare()
    if (lwrite_2daverages) call write_2daverages_prepare()
!
!  Find out which pencils to calculate at current time-step.
!
    lpencil = lpenc_requested
    if (lout)   lpencil=lpencil .or. lpenc_diagnos
    if (l2davg) lpencil=lpencil .or. lpenc_diagnos2d
    if (lvideo) lpencil=lpencil .or. lpenc_video
!
!  Save state vector prior to update for the (implicit) ADI scheme.
!
    if (lADI) f(:,:,:,iTTold)=f(:,:,:,ilnTT)
!
    if (ltestperturb) call testperturb_begin(f,df)
!
!  A random phase for the hydro_kinematic module
!
    if (lhydro_kinematic) call kinematic_random_phase()
!
!  Time advance.
!
    call rk_2n(f,df,p)
!
!  Print diagnostic averages to screen and file.
!
    if (lout) then
      call prints()
      if (lchemistry_diag) call write_net_reaction
     endif
    if (l1davg) call write_1daverages()
    if (l2davg) call write_2daverages()
    if (lwrite_sound) call write_sound()
!
!  Ensure better load balancing of particles by giving equal number of
!  particles to each CPU. This only works when block domain decomposition of
!  particles is activated.
!
    if (lparticles) call particles_load_balance(f)
!
!  07-Sep-07/dintrans+gastine: Implicit advance of the radiative diffusion
!  in the temperature equation (using the implicit_physics module).
!
    if (lADI) call calc_heatcond_ADI(f)
    if (ltestperturb) call testperturb_finalize(f)
!
    if (lroot) icount=icount+1  !  reliable loop count even for premature exit
!
!  Update time averages and time integrals.
!
    if (ltavg) call update_timeavgs(f,dt)
!
!  Add forcing and/or do rescaling (if applicable).
!
    if (lforcing) call addforce(f)
    if (lrescaling_magnetic)  call rescaling_magnetic(f)
    if (lrescaling_testfield) call rescaling_testfield(f)
    if (lrescaling_testscalar) call rescaling_testscalar(f)
!
!  Check for SNe, and update f if necessary (see interstellar.f90).
!
    if (linterstellar) call check_SN(f)
!
!  Check if mass flux replacement required fred test
!
    if (linterstellar) call addmassflux(f)
!
!  Check wall clock time, for diagnostics and for user supplied simulation time
!  limit.
!
    if (lroot.and.(idiag_walltime/=0.or.max_walltime/=0.0)) then
      time2=mpiwtime()
      wall_clock_time=(time2-time1)
      if (lout.and.idiag_walltime/=0) &
          call save_name(wall_clock_time,idiag_walltime)
    endif
!
    if (lout.and.lroot.and.idiag_timeperstep/=0) then
      it_this_diagnostic   = it
      time_this_diagnostic = mpiwtime()
      time_per_step = (time_this_diagnostic - time_last_diagnostic) &
                     /(  it_this_diagnostic -   it_last_diagnostic)
      it_last_diagnostic   =   it_this_diagnostic
      time_last_diagnostic = time_this_diagnostic
      call save_name(time_per_step,idiag_timeperstep)
    endif
!
!  Setting ialive=1 can be useful on flaky machines!
!  Each processor writes it's processor number (if it is alive!)
!  Set ialive=0 to fully switch this off.
!
    if (ialive /= 0) then
      if (mod(it,ialive)==0) &
          call outpui(trim(directory)//'/alive.info', &
          spread(it,1,1) ,1) ! (all procs alive?)
    endif
    if (lparticles) &
        call write_snapshot_particles(directory_snap,f,ENUM=.true.)
!
    call wsnap(trim(directory_snap)//'/VAR',f, &
        mvar_io,ENUM=.true.,FLIST='varN.list')
    call wsnap_timeavgs(trim(directory_snap)//'/TAVG',ENUM=.true., &
         FLIST='tavgN.list')
!
!  Write slices (for animation purposes).
!
    if (lvideo.and.lwrite_slices) call wvid(f,trim(directory)//'/slice_')
!
!  Save snapshot every isnap steps in case the run gets interrupted.
!  The time needs also to be written.
!
    if (isave/=0.and..not.lnowrite) then
      if (mod(it,isave)==0) then
        if (lparticles) &
            call write_snapshot_particles(directory_snap,f,ENUM=.false.)
!
        call wsnap(trim(directory_snap)//'/var.dat', &
                   f,mvar_io,ENUM=.false.,noghost=noghost_for_isave)
        call wsnap_timeavgs(trim(directory_snap)//'/timeavg.dat', &
                            ENUM=.false.)
        call wtime(trim(directory)//'/time.dat',t)
      endif
    endif
!
!  Save spectrum snapshot.
!
    if (dspec/=impossible) call powersnap(f)
!
!  Save global variables.
!
    if (isaveglobal/=0) then
      if (mod(it,isaveglobal)==0) then
        if (mglobal/=0)  &
          call output_globals(trim(directory_snap)//'/global.dat', &
              f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
      endif
    endif
!
!  Do exit when timestep has become too short.
!  This may indicate an MPI communication problem, so the data are useless
!  and won't be saved!
!
    if ((it<nt) .and. (dt<dtmin)) then
      if (lroot) &
          write(*,*) ' Time step has become too short: dt = ', dt
      save_lastsnap=.false.
      exit Time_loop
    endif
!
!  Exit do loop if wall_clock_time has exceeded max_walltime.
!
    if (max_walltime>0.0) then
      if (lroot.and.(wall_clock_time>max_walltime)) timeover=.true.
      call mpibcast_logical(timeover,1)
      if (timeover) then
        if (lroot) then
          print*
          print*, 'Maximum walltime exceeded'
        endif
        exit Time_loop
      endif
    endif
!
!  Fatal errors sometimes occur only on a specific processor. In that case all
!  processors must be informed about the problem before the code can stop.
!
    call fatal_error_local_collect()
    call timing('run','at the end of Time_loop',INSTRUCT='finalize')
!
    it=it+1
    headt=.false.
  enddo Time_loop
!
  if (lroot) then
    print*
    print*, 'Simulation finished after ', icount, ' time-steps'
  endif
!
  if (lroot) time2=mpiwtime()
!
!  Write data at end of run for restart.
!
  if (lroot) then
    print*
    print*, 'Writing final snapshot at time t =', t
  endif
!
  call wtime(trim(directory)//'/time.dat',t)
!
  if (.not.lnowrite) then
    if (save_lastsnap) then
      if (lparticles) &
          call write_snapshot_particles(directory_snap,f,ENUM=.false.)
!
      call wsnap(trim(directory_snap)//'/var.dat',f,mvar_io,ENUM=.false.)
      call wsnap_timeavgs(trim(directory_snap)//'/timeavg.dat',ENUM=.false.)
!
!  dvar is written for analysis and debugging purposes only.
!
      if (ip<=11 .or. lwrite_dvar) then
        call wsnap(trim(directory)//'/dvar.dat',df,mvar,ENUM=.false., &
            noghost=.true.)
        call particles_write_dsnapshot(trim(directory)//'/dpvar.dat',f)
      endif
!
!  Write crash files before exiting if we haven't written var.dat already
!
    else
      call wsnap(trim(directory_snap)//'/crash.dat',f,mvar_io,ENUM=.false.)
      if (ip<=11) call wsnap(trim(directory)//'/dcrash.dat',df,mvar,ENUM=.false.)
    endif
  endif
!
!  Save spectrum snapshot.
!
  if (save_lastsnap) then
    if (dspec/=impossible) call powersnap(f,.true.)
  endif
!
!  Print wall clock time and time per step and processor for diagnostic
!  purposes.
!
  if (lroot) then
    wall_clock_time=time2-time1
    print*
    write(*,'(A,1pG10.3,A,1pG9.2,A)') &
        ' Wall clock time [hours] = ', wall_clock_time/3600.0, &
        ' (+/- ', real(mpiwtick())/3600.0, ')'
    if (it>1) then
      if (lparticles) then
        write(*,'(A,1pG10.3)') &
            ' Wall clock time/timestep/(meshpoint+particle) [microsec] =', &
            wall_clock_time/icount/(nw+npar/ncpus)/ncpus/1.0e-6
      else
        write(*,'(A,1pG10.3)') &
            ' Wall clock time/timestep/meshpoint [microsec] =', &
            wall_clock_time/icount/nw/ncpus/1.0e-6
      endif
    endif
    print*
  endif
!
!  Stop MPI.
!
  call mpifinalize
!
!  Free any allocated memory.
!
  call farray_clean_up()
  call sharedvars_clean_up()
  call chemistry_clean_up()
  call NSCBC_clean_up()
  call vnames_clean_up()
  call xyaverages_clean_up()
  call xzaverages_clean_up()
  call yzaverages_clean_up()
  if (lwrite_phizaverages) call phizaverages_clean_up()
  if (lwrite_yaverages)    call yaverages_clean_up()
  if (lwrite_zaverages)    call zaverages_clean_up()
  if (lwrite_phiaverages)  call phiaverages_clean_up()
  if (lwrite_sound)        call sound_clean_up()
!
endprogram run
