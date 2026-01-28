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
module Run_module

    use Cdata
!
    implicit none

    integer :: icount, it_last_diagnostic
    real(KIND=rkind8) :: time1, time_last_diagnostic
    real :: time_doing_diagnostics = 0., time_in_timestep = 0.

contains
!***********************************************************************
subroutine helper_loop(f,p)
!
  use Boundcond, only: update_ghosts
  use Equ, only: perform_diagnostics
  use Diagnostics, only:  restore_diagnostic_controls,allocate_fnames
!$ use General, only: signal_wait, signal_send
  use Snapshot, only: perform_powersnap, perform_wsnap_ext, perform_wsnap_down
  use Mpicomm, only: mpiwtime
  use Sub, only: check_for_nans_globally
!
  real, dimension (mx,my,mz,mfarray) :: f
  type (pencil_case) :: p

  real :: start_time,end_time
!
! 7-feb-24/TP: coded
    call allocate_fnames(nname)
!
!$  do while(lhelper_run)
!$    call signal_wait(lhelper_perf,lhelper_run)
      call check_for_nans_globally(f,'after reading from the GPU')
      start_time = mpiwtime()
!$    if (lhelper_run) call restore_diagnostic_controls

!$    if (lhelper_run) call update_ghosts(f)
!$    if (lhelper_run .and. lhelperflags(PERF_DIAGS)) then
        call perform_diagnostics(f,p)
!$    else
!$      lhelperflags(PERF_DIAGS) = .false.
!$    endif
!$    if (lhelper_run .and. lhelperflags(PERF_WSNAP)) then
        call perform_wsnap_ext(f)
!$    else
!$      lhelperflags(PERF_WSNAP) = .false.
!$    endif
!$    if (lhelper_run .and. lhelperflags(PERF_WSNAP_DOWN)) then
        call perform_wsnap_down(f)
!$    else
!$      lhelperflags(PERF_WSNAP_DOWN) = .false.
!$    endif
!$    if (lhelper_run .and. lhelperflags(PERF_POWERSNAP)) then
        call perform_powersnap(f)
!$    else
!$      lhelperflags(PERF_POWERSNAP) = .false.
!$    endif
!$    call signal_send(lhelper_perf,.false.)
      end_time = mpiwtime()
      time_doing_diagnostics = time_doing_diagnostics + end_time-start_time

!$  enddo

endsubroutine helper_loop
!***********************************************************************
  subroutine reload(f, lreload_file, lreload_always_file)

    use GPU,         only: copy_farray_from_GPU, reload_GPU_config, load_farray_to_GPU
    use Register,    only: initialize_modules, rprint_list, choose_pencils
    use Sub,         only: control_file_exists
    use Param_IO,    only: read_all_init_pars, read_all_run_pars, write_all_run_pars
    use Boundcond,   only: initialize_boundcond
    use Forcing,     only: forcing_clean_up
    use Hydro,       only: hydro_clean_up
    use Solid_cells, only: solid_cells_clean_up
    use Timestep,    only: initialize_timestep
    use HDF5_IO,     only: initialize_hdf5
    use Diagnostics, only: report_undefined_diagnostics,diagnostics_clean_up
    use Particles_main,   only: particles_rprint_list, particles_initialize_modules
!$ use General, only: signal_wait, signal_send

    real, dimension (mx,my,mz,mfarray) :: f
    logical :: lreload_file, lreload_always_file
    real :: dtmp

    if (lroot) write(*,*) 'Found RELOAD file -- reloading parameters'
!
!  Re-read configuration
!
! If rkf timestep retain the current dt for continuity rather than reset
! with option to initialize_timestep to dt0/=0 from run.in if intended.
!
!TP: important for synchronization purposes that this happens before anything else
    if (lgpu) call copy_farray_from_GPU(f)
! Free the lock as if reloading never happened
!$  if (lfarray_copied) call signal_send(lhelper_perf,.false.)
!$  lfarray_copied = .false.

    if (ldt) then
      dt0=dt
      dt=0.0
    else
      dtmp=dt
    endif
    call read_all_run_pars
    if (.not.ldt) dt=dtmp
    if (ldt .and. .not. lcourant_dt) dt0 = dt_next
!
!  Before reading the rprint_list deallocate the arrays allocated for
!  1-D and 2-D diagnostics.
!
    call diagnostics_clean_up
    if (lforcing)            call forcing_clean_up
    if (lhydro_kinematic)    call hydro_clean_up
    if (lsolid_cells)        call solid_cells_clean_up

    call rprint_list(LRESET=.true.) !(Re-read output list)
    if (lparticles) call particles_rprint_list(.false.) !MR: shouldn't this be called with lreset=.true.?
    call report_undefined_diagnostics

    call initialize_hdf5
    call initialize_timestep
    call initialize_modules(f)
    call initialize_boundcond(f)
    if (lparticles) call particles_initialize_modules(f)
    call choose_pencils
    if (lgpu) then
      lpencil = lpenc_requested
      call reload_GPU_config
      call load_farray_to_GPU(f)
    endif
    call write_all_run_pars('IDL')       ! data to param2.nml
    call write_all_run_pars              ! diff data to params.log
!
    lreload_file=control_file_exists('RELOAD', DELETE=.true.)
    lreload_file        = .false.
    lreload_always_file = .false.
    lreloading          = .false.

  endsubroutine reload
!***********************************************************************
  subroutine gen_output(f)
!
! 5-sep-2024/TP: extracted from timeloop
!
    use Equ,             only: write_diagnostics
    use Snapshot,        only: powersnap, powersnap_prepare, wsnap, wsnap_down, output_form
    use Particles_main,  only: write_snapshot_particles
    use PointMasses,     only: pointmasses_write_snapshot
    use Mpicomm,         only: mpiwtime
    use Sub,             only: control_file_exists
    use Solid_Cells,     only: wsnap_ogrid
    use Timeavg,         only: wsnap_timeavgs
    use Fixed_point,     only: wfixed_points
    use Streamlines,     only: wtracers
    use Io,              only: output_globals

    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
    integer :: isave_shift=0
    real :: tvar1
!
!  Setting ialive=1 can be useful on flaky machines!
!  The iteration number is written into the file "data/proc*/alive.info".
!  Set ialive=0 to fully switch this off.
!
    if (ialive /= 0) then
      if (mod(it,ialive)==0) call output_form('alive.info',it,.false.)
    endif
    if (lparticles) call write_snapshot_particles(f,ENUM=.true.)  !MR: no *.list file for particles?
    if (lpointmasses) call pointmasses_write_snapshot('QVAR',ENUM=.true.,FLIST='qvarN.list')
!
!  Added possibility of outputting only the chunks starting
!  from nv1_capitalvar in the capitalvar file.
!
    call wsnap('VAR',f,mvar_io,ENUM=.true.,FLIST='varN.list',nv1=nv1_capitalvar)
    if (ldownsampl) call wsnap_down(f)
    call wsnap_timeavgs('TAVG',ENUM=.true.,FLIST='tavgN.list')
    !MR: what about ogrid data here?
!
!   Diagnostic output in concurrent thread.
!
    if (lmultithread) then
      if (lout.or.l1davg.or.l1dphiavg.or.l2davg) then
!!$      last_pushed_task= push_task(c_funloc(write_diagnostics_wrapper), last_pushed_task,&
!!$      1, default_task_type, 1, depend_on_all, f, mx, my, mz, mfarray)
      endif
    else
      call write_diagnostics(f)
    endif
!
!  Write tracers (for animation purposes).
!
    if (ltracers.and.lwrite_tracers) call wtracers(f,trim(directory)//'/tracers_')
!
!  Write fixed points (for animation purposes).
!
    if (lfixed_points.and.lwrite_fixed_points) call wfixed_points(f,trim(directory)//'/fixed_points_')
!
!  Save snapshot every isnap steps in case the run gets interrupted or when SAVE file is found.
!  The time needs also to be written.
!
    if (lout) lsave = control_file_exists('SAVE', DELETE=.true.)

    if (lsave .or. ((isave /= 0) .and. .not. lnowrite)) then
      if (lsave .or. (mod(it-isave_shift, isave) == 0)) then
        lsave = .false.
        if (ip<=12.and.lroot) tvar1=mpiwtime()
        call wsnap('var.dat',f, mvar_io,ENUM=.false.,noghost=noghost_for_isave)
        if (ip<=12.and.lroot) print*,'wsnap: written snapshot var.dat in ', &
                                     mpiwtime()-tvar1,' seconds'
        call wsnap_timeavgs('timeavg.dat',ENUM=.false.)
        if (lparticles) call write_snapshot_particles(f,ENUM=.false.)
        if (lpointmasses) call pointmasses_write_snapshot('qvar.dat',ENUM=.false.)
        if (lsave) isave_shift = mod(it+isave-isave_shift, isave) + isave_shift
        if (lsolid_cells) call wsnap_ogrid('ogvar.dat',ENUM=.false.)
      endif
    endif
!
!  Allow here for the possibility to have spectral output
!  *after* the first time. As explained in the comment to
!  the GW module, the stress spectrum is only available
!  when the GW solver has advanced by one step, but the time
!  of the stress spectrum is said to be t+dt, even though
!  it really belongs to the time t. The GW spectra, on the
!  other hand, are indeed at the correct d+dt. Therefore,
!  when lspec_first=T, we output spectra for both t and t+dt.
!
    if (lspec_start .and. abs(t-(tstart+dt))<.1*dt ) lspec=.true.
!
!  Save spectrum snapshot.
!
    if (dspec/=impossible .or. itspec/=impossible_int) call powersnap(f)
    if (lroot.and.(dspec/=impossible .or. itspec/=impossible_int).and.lspec) print*, 'gen_output powersnap'
!
!  Save global variables.
!
    if (isaveglobal/=0) then
      if ((mod(it,isaveglobal)==0) .and. (mglobal/=0)) &
        call output_globals('global.dat',f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
    endif

  endsubroutine gen_output
!***********************************************************************
  subroutine timeloop(f,df,p)
!
!  Do loop in time.
!
! 26-feb-24/MR: carved out from main_sub
!
  use Density,         only: boussinesq
  use Diagnostics,     only: write_1daverages_prepare, save_name, &
                             diagnostics_clean_up, write_2daverages_prepare
  use Filter,          only: rmwig, rmwig_xyaverage
  use Fixed_point,     only: fixed_points_prepare
  use Forcing,         only: addforce
  use ImplicitPhysics, only: calc_heatcond_ADI
  use IO,              only: output_globals
  use Magnetic,        only: rescaling_magnetic
  use Messages,        only: timing, fatal_error_local_collect
  use Mpicomm,         only: mpibcast_logical, mpiwtime, MPI_COMM_PENCIL, mpibarrier
  use Particles_main,  only: particles_rprint_list, particles_initialize_modules, &
                             particles_load_balance, particles_stochastic
  use Signal_handling, only: emergency_stop
  use Sub,             only: control_file_exists, calc_scl_factor
  use Testscalar,      only: rescaling_testscalar
  use Testfield,       only: rescaling_testfield
  use TestPerturb,     only: testperturb_begin, testperturb_finalize
  use Timeavg,         only: update_timeavgs
  use Timestep,        only: time_step
  use Slices,          only: wvid_prepare
  use Solid_Cells,     only: time_step_ogrid
  use Streamlines,     only: tracers_prepare
  use Snapshot,        only: powersnap_prepare
  use GPU,             only: gpu_set_dt
!$ use OMP_lib
!$ use General, only: signal_send, signal_wait
!
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (mx,my,mz,mvar) :: df
  type (pencil_case) :: p
!
  logical :: lstop=.false., timeover=.false., resubmit=.false., lreload_file, lreload_always_file, &
             lonemorestep=.false.
  real :: wall_clock_time=0., time_per_step=0., timer_for_timestep
  real(KIND=rkind8) :: time_this_diagnostic
  integer :: it_this_diagnostic
!
!TP: Due to df being always limited to a kernel on the Astaroth side we have to know the timestep before we do the rhs
!    compared to the cpu where it is sufficient to know it after the rhs calculations
!    so we take the timestep calculated at the start of the last timestep
!    initially there is no previous timestep so we have a extra call here for there always to be a valid previous timestep
!    no advancement happens here
!
  if (lgpu .and. lcourant_dt .and. ldt) call gpu_set_dt()

  Time_loop: do while (it<=nt)
!
!  Possibility to turn off logspacing for time series output
!
    if (.not.(lit1_logspacing.and.real(t)<tmax_logspacing)) &
      lout = (mod(it-1,it1) == 0) .and. (it > it1start)
!
    if (lout .or. emergency_stop) then
!
!  Exit do loop if file `STOP' exists.
!
      lstop=control_file_exists('STOP',DELETE=.true.)
      if (lstop .or. emergency_stop) then
        if (lroot) then
          print*
          if (emergency_stop) print*, 'Emergency stop requested'
          if (lstop) print*, 'Found STOP file'
        endif
        resubmit=control_file_exists('RESUBMIT',DELETE=.true.)
        if (resubmit) print*, 'Cannot be resubmitted'
!$      if (lfarray_copied) call signal_send(lhelper_perf,.true.)
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
      if (lreloading) call reload(f, lreload_file, lreload_always_file)
    endif
!
!  Calculate scale factor of the universe.
!  The results of this are used in are used in src/special/disp_current.f90.
!  This routine is independent of the routine src/special/Lambda_CDM.f90,
!  which integrates ascale, Hubble, and tphys that are used in other routines.
!
    if (lread_scl_factor_file_new) call calc_scl_factor
!
    if (lwrite_sound) then
      if ( .not.lout_sound .and. abs( t-tsound - dsound )<= 1.1*dt ) then
        lout_sound = .true.
        tsound = t
      endif
    endif
!
!  1-D timestep control
!
    if (lwrite_1daverages) then
      if (d1davg==impossible) then
        l1davg = (mod(it-1,it1d) == 0)
      else
        call write_1daverages_prepare(t == 0.0 .and. lwrite_ic)
      endif
    endif
!
!  Average removal timestep control.
!
    if (it_rmv==0) then
      lrmv=.true.
    else
      lrmv = (mod(it-1,it_rmv) == 0)
    endif
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
!  If we want to write out video data, wvid_prepare sets lvideo=.true.
!  This allows pde to prepare some of the data.
!
    if (lwrite_slices) then
      call wvid_prepare
      if (t == 0.0 .and. lwrite_ic) lvideo = .true.
    endif
!
    if (lwrite_2daverages) call write_2daverages_prepare(t == 0.0 .and. lwrite_ic)
!
!  Exit do loop if maximum simulation time is reached; allow one extra
!  step if any diagnostic output needs to be produced.
!
    if (t >= tmax) then
      if (lonemorestep .or. .not. (lout .or. lvideo .or. l2davg)) then
        if (lroot) print *, 'Maximum simulation time exceeded'
!$      if (lfarray_copied) call signal_send(lhelper_perf,.true.)
        exit Time_loop
      endif
      lonemorestep = .true.
    endif
!
!   Prepare for the writing of the tracers and the fixed points.
!
    if (lwrite_tracers) call tracers_prepare
    if (lwrite_fixed_points) call fixed_points_prepare
!
!  Find out which pencils to calculate at current time-step.
!
    lpencil = lpenc_requested
!  MR: the following should only be done in the first substep, shouldn't it?
!  AB: yes, so this is the right place, right?
!  MR: I mean, it should be reversed to lpencil = lpenc_requested for 2nd and further substeps.
!  AB: the line above says lpencil = lpenc_requested, and now you want to revert to it. How?
    if (l1davg.or.lout) lpencil=lpencil .or. lpenc_diagnos
    if (l2davg) lpencil=lpencil .or. lpenc_diagnos2d
    if (lvideo) lpencil=lpencil .or. lpenc_video
!
!  Save state vector prior to update for the (implicit) ADI scheme.
!
    if (lADI) f(:,:,:,iTTold)=f(:,:,:,iTT)
!
    if (ltestperturb) call testperturb_begin(f,df)
!
!  Decide here whether or not we will need a power spectrum.
!  At least for the graviational wave spectra, this requires
!  advance warning so the relevant components of the f-array
!  can be filled.
!
    if (dspec==impossible) then
      lspec = (mod(it-1,itspec) == 0)
    else
      call powersnap_prepare
    endif
!
!  Possibility of logspacing for time series output
!
    if (lit1_logspacing.and.real(t)<tmax_logspacing) lout=lspec 
!
!  Time advance.
!
    timer_for_timestep = mpiwtime()
    call time_step(f,df,p)
    time_in_timestep = time_in_timestep + mpiwtime()-timer_for_timestep
!    tdiagnos=t
!
!  If overlapping grids are used to get body-confined grid around the solids
!  in the flow, call time step on these grids.
!
    if (lsolid_cells) call time_step_ogrid(f)
!
    lrmv=.false.
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
!
    if (ltestperturb) call testperturb_finalize(f)
!
    if (lboussinesq) call boussinesq(f)
!
    if (lroot) icount=icount+1  !  reliable loop count even for premature exit
!
!  Update time averages and time integrals.
!
    if (ltavg) call update_timeavgs(f,dt)
!
!  Add forcing and/or do rescaling (if applicable).
!
    if (lforcing.and..not.lgpu) call addforce(f)
    if (lparticles_lyapunov) call particles_stochastic
!    if (lspecial) call special_stochastic
    if (lrescaling_magnetic)  call rescaling_magnetic(f)
    if (lrescaling_testfield) call rescaling_testfield(f)
    if (lrescaling_testscalar) call rescaling_testscalar(f)
!
!  Check wall clock time, for diagnostics and for user supplied simulation time
!  limit.
!
    if (lroot.and.(idiag_walltime/=0.or.max_walltime/=0.0)) then
      wall_clock_time=mpiwtime()-time1
      if (lout) call save_name(wall_clock_time,idiag_walltime)
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

    call gen_output(f)
!
!  Do exit when timestep has become too short.
!  This may indicate an MPI communication problem, so the data are useless
!  and won't be saved!
!
    if ((it<nt) .and. (dt<dtmin)) then
      if (lroot) write(*,*) ' Time step has become too short: dt = ', dt
      save_lastsnap=.false.
!$    if (lfarray_copied) call signal_send(lhelper_perf,.true.)
      exit Time_loop
    endif
!
!  Exit do loop if wall_clock_time has exceeded max_walltime.
!
    if (max_walltime>0.0) then
      if (lroot.and.(wall_clock_time>max_walltime)) timeover=.true.
      call mpibcast_logical(timeover,comm=MPI_COMM_PENCIL)
      if (timeover) then
        if (lroot) then
          print*
          print*, 'Maximum walltime exceeded'
        endif
!$      if (lfarray_copied) call signal_send(lhelper_perf,.true.)
        exit Time_loop
      endif
    endif
!
!  Fatal errors sometimes occur only on a specific processor. In that case all
!  processors must be informed about the problem before the code can stop.
!
    call fatal_error_local_collect
    call timing('run','at the end of Time_loop',INSTRUCT='finalize')
!
    it=it+1
    headt=.false.

!if (any(lmasterflags)) then
!write(60+iproc,*) it, lmasterflags
!flush(60+iproc)
!endif
!$  if (lfarray_copied .and. any(lmasterflags)) then
!$    lhelperflags = lmasterflags
!$    lmasterflags = .false.
!$    call signal_send(lhelper_perf,.true.)
!$    lfarray_copied = .false.
!$  else if(lfarray_copied) then
!$    lfarray_copied = .false.
!$    call signal_send(lhelper_perf,.false.)
!$  endif

  enddo Time_loop

!$ call signal_wait(lhelper_perf, .false.)
!$ call signal_send(lhelper_run,.false.)
!Important that this comes after synchronization
!$ lmultithread = .false.

  endsubroutine timeloop
!***********************************************************************
  subroutine run_start() bind(C)
!
!  8-mar-13/MR: changed calls to wsnap and rsnap to grant reference to f by
!               address
! 31-oct-13/MR: replaced rparam by read_all_init_pars
! 10-feb-14/MR: initialize_mpicomm now called before read_all_run_pars
! 13-feb-13/MR: call of wsnap_down added
! 7-feb-24/TP: made main_sub for easier multithreading
!
  use Boundcond,       only: update_ghosts, initialize_boundcond
  use Chemistry,       only: chemistry_clean_up
  use Diagnostics,     only: phiavg_norm, report_undefined_diagnostics, trim_averages,diagnostics_clean_up
  use Equ,             only: initialize_pencils, debug_imn_arrays, rhs_sum_time, before_boundary_sum_time,&
                             radtransfer_sum_time,time_spent_copying_and_waiting
  use FArrayManager,   only: farray_clean_up
  use Farray_alloc
  use General,         only: random_seed_wrapper, touch_file, itoa
!$ use General,        only: signal_init, get_cpu
  use Grid,            only: construct_grid, box_vol, grid_bound_data, set_coorsys_dimmask, &
                             construct_serial_arrays, coarsegrid_interp
  use Gpu,             only: load_farray_to_GPU, initialize_gpu
  use HDF5_IO,         only: init_hdf5, initialize_hdf5, wdim
  use File_io,         only: file_exists,delete_file
  use IO,              only: rgrid, wgrid, directory_names, rproc_bounds, wproc_bounds, output_globals, input_globals
  use Messages
  use Mpicomm
  use NSCBC,           only: NSCBC_clean_up
  use Param_IO,        only: read_all_init_pars, read_all_run_pars, write_all_run_pars, write_pencil_info, get_downpars
  use Particles_main
  use Pencil_check,    only: pencil_consistency_check
  use PointMasses,     only: pointmasses_read_snapshot, pointmasses_write_snapshot
  use Python,          only: python_init, python_finalize
  use Register
  use SharedVariables, only: sharedvars_clean_up
  use Signal_handling, only: signal_prepare
  use Slices,          only: setup_slices
  use Snapshot,        only: powersnap, rsnap, wsnap
  use Solid_Cells,     only: wsnap_ogrid
  use Special,         only: initialize_mult_special
  use Sub,             only: control_file_exists, get_nseed
  use Syscalls,        only: memusage, sizeof_real
  use Timeavg,         only: wsnap_timeavgs
  use Timestep,        only: initialize_timestep,after_substep_sum_time
  use Training,        only: training_time, inference_time
!
!$ use OMP_lib
!$ use, intrinsic :: iso_c_binding
!
  implicit none

  type (pencil_case) :: p

  character(len=fnlen) :: fproc_bounds
  real(KIND=rkind8) :: time2, tvar1
  real :: wall_clock_time=0.
  integer :: mvar_in
  integer :: memuse, memory, memcpu
  logical :: suppress_pencil_check=.false.
  logical :: lnoreset_tzero=.false.
  logical :: lprocbounds_exist
  integer, parameter :: num_helper_masters=1
  integer :: i,j
  integer :: master_core_id
  integer :: helper_core_id
  integer, dimension(max_threads_possible) :: tmp_core_ids
!
  lrun = .true.
!
!  Get processor numbers and define whether we are root.
!
  call mpicomm_init
!
!  Initialize OpenMP use
!
!$ include 'omp_init.h'
!
!  Initialize Python use.
!
  call python_init
!
!  Identify version.
!
  if (lroot) call svn_id('$Id$')
!
!  Initialize use of multiple special modules if relevant.
!
  call initialize_mult_special
!
!  Read parameters from start.x (set in start.in/default values; possibly overwritten by 'read_all_run_pars').
!
  call read_all_init_pars
!
!  Read parameters and output parameter list.
!
  call read_all_run_pars
!
!  Initialize the message subsystem, eg. color setting etc.
!
  call initialize_messages
!
!  Initialise MPI communication.
!
  call initialize_mpicomm
!
!  Initialise HDF5 library.
!
  call init_hdf5
!
!  Initialize HDF_IO module.
!
  call initialize_hdf5
!
!  Check whether quad precision is supported
!
  if (rkind16<0) call warning('run','quad precision not supported, switch to double')
  if (rkind16==rkind8) call warning('run','quad precision suppressed')
!
  if (any(downsampl>1) .or. mvar_down>0 .or. maux_down>0) then
!
!  If downsampling, calculate local start indices and number of data in
!  output for each direction; inner ghost zones are here disregarded
!
    ldownsampl = .true.
    if (dsnap_down<=0.) dsnap_down=dsnap
!
    call get_downpars(1,nx,ipx)
    call get_downpars(2,ny,ipy)
    call get_downpars(3,nz,ipz)
!
    if (any(ndown==0)) &
      call fatal_error('run','zero points in processor '//trim(itoa(iproc))//' for downsampling')
      ! MR: better a warning and continue w. ldownsampl=.false.?

    call mpiallreduce_sum_int(ndown(1),ngrid_down(1),comm=MPI_COMM_XBEAM)
    call mpiallreduce_sum_int(ndown(2),ngrid_down(2),comm=MPI_COMM_YBEAM)
    call mpiallreduce_sum_int(ndown(3),ngrid_down(3),comm=MPI_COMM_ZBEAM)

  endif
!
!  Set up directory names.
!
  call directory_names
!
!  Read coordinates (if luse_oldgrid=T, otherwise regenerate grid).
!  luse_oldgrid=T can be useful if nghost_read_fewer > 0,
!  i.e. if one is changing the order of spatial derivatives.
!  Also write dim.dat (important when reading smaller meshes, for example)
!  luse_oldgrid=.true. by default, and the values written at the end of
!  each var.dat file are ignored anyway and only used for postprocessing.
!
  call set_coorsys_dimmask
!
  fproc_bounds = trim(datadir) // "/proc_bounds.dat"
  inquire (file=fproc_bounds, exist=lprocbounds_exist)
  call mpibarrier
!
  if (luse_oldgrid .and. lprocbounds_exist) then
    if (ip<=6.and.lroot) print*, 'reading grid coordinates'
    call rgrid('grid.dat')
    call rproc_bounds(fproc_bounds)
    call construct_serial_arrays
    call grid_bound_data
  else
    if (luse_oldgrid) call warning("run", "reconstructing the grid because proc_bounds.dat is missing")
    if (luse_xyz1) Lxyz = xyz1-xyz0
    call construct_grid(x,y,z,dx,dy,dz)
    lprocbounds_exist = .false.    ! triggers wproc_bounds later
  endif
!
!  Store metadata was the run in double or single precision
!  We do this instead of e.g. reading it from dim.dat since dim.dat won't exist
!  if we use HDF5-IO
!
    if(lroot) then
      if (file_exists('data/SINGLE_PRECISION_RUN')) call delete_file('data/SINGLE_PRECISION_RUN')
      if (file_exists('data/DOUBLE_PRECISION_RUN')) call delete_file('data/DOUBLE_PRECISION_RUN')
      if(sizeof_real() < 8) then
        call touch_file('data/SINGLE_PRECISION_RUN')
      else
        call touch_file('data/DOUBLE_PRECISION_RUN')
      endif
    endif

!
!  Shorthands (global).
!
  x0 = xyz0(1) ; y0 = xyz0(2) ; z0 = xyz0(3)
  Lx = Lxyz(1) ; Ly = Lxyz(2) ; Lz = Lxyz(3)
!
!  Register physics (incl. particles) modules.
!
  call register_modules
  if (lparticles) call particles_register_modules
  call initialize
!
!
!  Inform about verbose level.
!
  if (lroot) print*, 'The verbose level is ip=', ip, ' (ldebug=', ldebug, ')'
!
!  Position of equator (if any).
!
  if (lequatory) yequator=xyz0(2)+0.5*Lxyz(2)
  if (lequatorz) zequator=xyz0(3)+0.5*Lxyz(3)
!
!  Print resolution and dimension of the simulation.
!
  if (lroot) then
    write(*,'(a,i1,a)') ' This is a ', dimensionality, '-D run'
    print*, 'nxgrid, nygrid, nzgrid=', nxgrid, nygrid, nzgrid
    print*, 'Lx, Ly, Lz=', Lxyz
    call box_vol
    if (lyinyang) then
      print*, '      Vbox(Yin,Yang)=', box_volume
      print*, '      total volume  =', 4./3.*pi*(xyz1(1)**3-xyz0(1)**3)
    else
      print*, '      Vbox=', box_volume
    endif
  endif
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
!  set slots to be written in downsampled f, if not specified by input, use mvar_io
!
  if (ldownsampl) then
    if (mvar_down<0) then
      mvar_down=mvar
    else
      mvar_down=min(mvar,mvar_down)
    endif
!
!  maux_down -1 default uses all maux, maux_down 0 no aux and maux_down > 0 use min(maux,maux_down)
!
    if (maux_down<0) then
      maux_down=maux
    else
      maux_down=min(maux,maux_down)
    endif
    if (mvar_down+maux_down==0) ldownsampl=.false.
    if (mvar_down<mvar.and.maux_down>0) then
      mvar_down=mvar
      call warning('run','mvar_down<mvar and maux_down>0 -> to avoid gap in variable sequence, mvar_down is set to mvar')
    endif
  endif
!
! Shall we read also auxiliary variables or fewer variables (ex: turbulence
! field with 0 species as an input file for a chemistry problem)?
!
  if (lread_aux) then
    mvar_in=mvar+maux
  else if (lread_less) then
    mvar_in=4        !MR: strange - why 4?
  else
    mvar_in=mvar
  endif
!
!  Get state length of random number generator and put the default value.
!  With lreset_seed (which is not the default) we can reset the seed during
!  the run. This is necessary when lreinitialize_uu=T, inituu='gaussian-noise'.
!
  call get_nseed(nseed)
  if (lreset_seed) then
    seed(1)=-((seed0-1812+1)*10+iproc_world)
    call random_seed_wrapper(PUT=seed,CHANNEL=1)
    call random_seed_wrapper(PUT=seed,CHANNEL=2)
  else
    call random_seed_wrapper(GET=seed,CHANNEL=1)
    seed = seed0
    call random_seed_wrapper(PUT=seed,CHANNEL=1)
!
    call random_seed_wrapper(GET=seed2,CHANNEL=2)
    seed2 = seed0
    call random_seed_wrapper(PUT=seed2,CHANNEL=2)
  endif
!
!  Write particle block dimensions to file (may have been changed for better
!  load balancing).
!
  if (lroot) then
    if (lparticles) call particles_write_block(trim(datadir)//'/bdim.dat')
  endif
!
!  Read data.
!  Snapshot data are saved in the data subdirectory.
!  This directory must exist, but may be linked to another disk.
!
  f=0.
  if (lroot .and. ldebug) print*, 'memusage before rsnap=', memusage()/1024., 'MBytes'
  if (lroot) tvar1=mpiwtime()
  call rsnap('var.dat',f,mvar_in,lread_nogrid)
  if (lroot) print*,'rsnap: read snapshot var.dat in ',mpiwtime()-tvar1,' seconds'
!
!  If we decided to use a new grid, we need to overwrite the data
!  that we just read in from var.dat. (Note that grid information
!  was also used above, so we really need to do it twice then.)
!
  if (.not.luse_oldgrid) call construct_grid(x,y,z,dx,dy,dz) !MR: already called
!
!  Initialize diagnostics and write indices to file.
!
  call rprint_list(LRESET=.false.)
  if (lparticles) call particles_rprint_list(.false.)
  if (lreport_undefined_diagnostics) call report_undefined_diagnostics
!
!  Read particle snapshot.
!
  if (lparticles) call read_snapshot_particles
!
!  Read point masses.
!
  if (lpointmasses) call pointmasses_read_snapshot('qvar.dat')
!
!  Set initial time to zero if requested. This is dangerous, however!
!  One may forget removing this entry after having set this once.
!  It is therefore safer to say lini_t_eq_zero_once=.true.,
!  which does the reset once once, unless NORESET_TZERO is removed.
!
  if (lini_t_eq_zero) t=0.0
!
!  Set initial time to zero if requested, but blocks further resets.
!  See detailed comment above.
!
  lnoreset_tzero=control_file_exists('NORESET_TZERO')
  if (lini_t_eq_zero_once.and..not.lnoreset_tzero) then
    call touch_file('NORESET_TZERO')
    t=0.0
  endif
!
!  Set last tsound output time
!
  if (lwrite_sound) then
    if (tsound<0.0) then
      ! if sound output starts new
      tsound=t
      ! output initial values
      lout_sound=.true.
    endif
  endif
!
!  The following is here to avoid division in sub.f90 for diagnostic
!  outputs of integrated values in the non equidistant case.
!  Do this even for uniform meshes, in which case xprim=dx, etc.
!  Remember that dx_1=0 for runs without extent in that direction.
!
  if (nxgrid==1) then; xprim=1.0; else; xprim=1./dx_1; endif
  if (nygrid==1) then; yprim=1.0; else; yprim=1./dy_1; endif
  if (nzgrid==1) then; zprim=1.0; else; zprim=1./dz_1; endif
!
!  Determine slice positions and whether slices are to be written on this
!  processor. This can only be done after the grid has been established.
!
  call setup_slices
!
!  Initialize the list of neighboring processes.
!
  call update_neighbors     !MR: Isn't this only needed for particles?
!
!  Allow modules to do any physics modules do parameter dependent
!  initialization. And final pre-timestepping setup.
!  (must be done before need_XXXX can be used, for example)
!
  if (lroot .and. ldebug) print*, 'memusage before initialize modules=', memusage()/1024., 'MBytes'
  call initialize_timestep
  call initialize_modules(f)
  call initialize_boundcond(f)
!
!TP: done here that at runtime_compilation we know which pencils are requested
!
!  Find out which pencils are needed and write information about required,
!  requested and diagnostic pencils to disc.
!
  call choose_pencils
  call write_pencil_info

  lpencil = lpenc_requested
  if (nt>0) call initialize_gpu(f)
!
  if (it1d==impossible_int) then
    it1d=it1
  else
    if (it1d<it1) call stop_it_if_any(lroot,'run: it1d smaller than it1')
  endif
!
!  Read global variables (if any).
!
  if (mglobal/=0 .and. lread_global) call input_globals('global.dat', &
      f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
!
!  Prepare particles.
!
  if (lparticles) call particles_initialize_modules(f)
!
!  Only after register it is possible to write the correct dim.dat
!  file with the correct number of variables.
!  No IO-module-controlled reading operations allowed beyond this point!
!
  call wgrid("grid.dat", lwrite=.not.(lprocbounds_exist .and. luse_oldgrid))
  if (.not.lprocbounds_exist) call wproc_bounds(fproc_bounds)

  if (.not.luse_oldgrid .or. lwrite_dim_again) then
    call wdim('dim.dat')
    if (ip<11 .and. lroot) then
      print*,'Lz=',Lz
      print*,'z=',z
    endif
  endif
!
!  Write data to file for IDL (param2.nml).
!
  call write_all_run_pars('IDL')
!
!  Write parameters to log file (done after reading var.dat, since we
!  want to output time t.
!
  call write_all_run_pars
!
!  Possible debug output (can only be done after "directory" is set).
!  Check whether mn array is correct.
!
  if (ip<=3) call debug_imn_arrays
!
  if (mglobal/=0) call output_globals('global.dat',f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
!
!  Update ghost zones, so rprint works corrected for at the first
!  time step even if we didn't read ghost zones.
!
  call update_ghosts(f)
!
!  Allow here for the possibility to have spectral output
!  from the first time. This works for all variables, except
!  for the GW stress, but the GW signal itself is at the
!  correct time; see the comment in the GW module.
!
  if (lspec_start .and. t==tstart) lspec=.true.
!
!  Save spectrum snapshot.
!
!  tdiagnos=t
  if (dspec/=impossible) call powersnap(f)
!
!  Initialize pencils in the pencil_case.
!
  if (lpencil_init) call initialize_pencils(p,0.0)
!
!  Perform pencil_case consistency check if requested.
!
  it = 1          ! needed for pencil check
  suppress_pencil_check = control_file_exists("NO-PENCIL-CHECK")
  if ( ((lpencil_check .and. .not. suppress_pencil_check) .or. &
        (.not.lpencil_check.and.lpencil_check_small)) .and. nt>0 ) then
    if (lgpu) then
      call warning('run',"Pencil consistency check not supported on GPUs. You can consider running it with a CPU-only build")
    else
      call pencil_consistency_check(f,df,p)
    endif
  endif
!
!  Globally catch eventual 'stop_it_if_any' call from single MPI ranks
!
  call stop_it_if_any(.false.,'')
!
!  Prepare signal catching
!
  call signal_prepare
!
!  Trim 1D-averages for times past the current time.
!
  call trim_averages

  !$ call mpibarrier
  !$omp parallel
  !$    core_ids(omp_get_thread_num()+1) = get_cpu()
  !$omp end parallel
  !$ call mpibarrier

!$omp parallel num_threads(num_helper_masters+1) &
!$omp copyin(dxmax_pencil,fname,fnamex,fnamey,fnamez,fnamer,fnamexy,fnamexz,fnamerz,fname_keep,fname_sound,ncountsz,phiavg_norm)
!
!TP: remove master id from core ids since no one should run on master core and make sure new core ids indexing start from 1
!$ if (omp_get_thread_num() == 1) helper_core_id = get_cpu()
!$omp barrier
!$ if (omp_get_thread_num() == 0) then
!$   master_core_id = get_cpu()
!$   tmp_core_ids = 0
!$   tmp_core_ids(1) = helper_core_id
!$   j = 2
!$   do i = 1,20
!$     if (core_ids(i) /= master_core_id .and. core_ids(i) /= helper_core_id) then
!$       tmp_core_ids(j) = core_ids(i)
!$       j = j +1
!$     endif
!$   enddo
!$   core_ids = tmp_core_ids
!$ endif
!$omp barrier
!
! Ensures that all MPI processes call create_communicators with the same thread.
!
!$  do i=1,num_helper_masters
!$omp barrier
!$    if (omp_get_thread_num() == i) call create_communicators
!$  enddo
!$omp barrier
    call mpibarrier
!$  if (omp_get_thread_num() == 0) then
!
!  Start timing for final timing statistics.
!  Initialize timestep diagnostics during the run (whether used or not,
!  see idiag_timeperstep).
!
!!$print*,"Master core: ",get_cpu(), iproc
      if (lroot) then
        icount=0
        it_last_diagnostic=icount
        time1=mpiwtime()
        time_last_diagnostic=time1
      endif
      if (nt>0) call timeloop(f,df,p)
!$  else
!$    if (nt>0) call helper_loop(f,p)
!$  endif
!$omp barrier
!$omp end parallel
!
  if (lroot) then
!
    time2=mpiwtime()
!
    print*
    print*, 'Simulation finished after ', icount, ' time-steps'
!
  endif
!
  if (.not.lnowrite .or. nt==0) then
!
!  Write data at end of run for restart.
!
    if (lroot) then
      print*
      write(*,*) 'Writing final snapshot at time t =', t
    endif
!
    if (ncoarse>1) then
      call update_ghosts(f)
      if (lcoarse) then
        call coarsegrid_interp(f)
        lcoarse=.false.
      endif
      call update_ghosts(f)
    endif

    if (save_lastsnap.or.nt==0) then
      if (lparticles) call write_snapshot_particles(f,ENUM=.false.)
      if (lpointmasses) call pointmasses_write_snapshot('qvar.dat',ENUM=.false.)
      if (lsolid_cells) call wsnap_ogrid('ogvar.dat',ENUM=.false.)
!
      call timing(message='start writing ',instruct='initialize',lforce=ltiming_io)
      call wsnap('var.dat',f,mvar_io,ENUM=.false.)
      call timing(message='end writing ',instruct='finalize',lforce=ltiming_io)
      call wsnap_timeavgs('timeavg.dat',ENUM=.false.)
!
!  dvar is written for analysis and debugging purposes only.
!
      if (ip<=11 .or. lwrite_dvar) then
        call wsnap('dvar.dat',df,mvar,ENUM=.false.,noghost=.true.)
        call particles_write_dsnapshot('dpvar.dat',f)
      endif
!
!  Write crash files before exiting if we haven't written var.dat already
!
    else
      call wsnap('crash.dat',f,mvar_io,ENUM=.false.)
      if (ip<=11) call wsnap('dcrash.dat',df,mvar,ENUM=.false.)
    endif
  endif
!
!  Save spectrum snapshot. The spectrum of the last timestep is written by
!  default, but this can be turned off by setting lwrite_last_powersnap=F
!  in run_pars or init_pars.
!
  if (save_lastsnap) then
    if (dspec/=impossible) call powersnap(f,lwrite_last_powersnap)
  endif
!
!  Print wall clock time and time per step and processor for diagnostic
!  purposes.
!
  if (lroot) then
    wall_clock_time=time2-time1
    print*
    write(*,'(A,1pG10.3,A,1pG11.4,A)') ' Wall clock time [hours] = ', wall_clock_time/3600.0, &
                                       ' (+/- ', real(mpiwtick())/3600.0, ')'
    if (it>1) then
      if (lparticles) then
        write(*,'(A,1pG10.3)') ' Wall clock time/timestep/(meshpoint+particle) [microsec] =', &
                               wall_clock_time/icount/(nw+npar/ncpus)/ncpus/1.0e-6
      else
        write(*,'(A,1pG14.7)') ' Wall clock time/timestep/meshpoint [microsec] =', &
                               wall_clock_time/icount/nw/ncpus/1.0e-6
        write(*,'(A,1pG14.7)') ' Wall clock time/timestep/local meshpoint [microsec] =', &
                               wall_clock_time/icount/nw/1.0e-6
        write(*,'(A,1pG14.7)') ' Rhs wall clock time/timestep/local meshpoint [microsec] =', &
                               rhs_sum_time/icount/nw/1.0e-6
! 2025-Dec-08/Kishore: Touko, I would suggest using the existing variable `ip` (which
! controls some other diagnostic outputs as well), rather than introducing a new
! flag.
        if (lverbose_performance_log) then
          write(*,'(A,1pG14.7)') &
            ' After substep wall clock time/timestep/local meshpoint [microsec] =', &
            after_substep_sum_time/icount/nw/1.0e-6
          write(*,'(A,1pG14.7)') &
            ' Before boundary wall clock time/timestep/local meshpoint [microsec] =', &
            before_boundary_sum_time/icount/nw/1.0e-6
          if(lradiation_ray) write(*,'(A,1pG14.7)') &
            ' Radtransfer wall clock time/timestep/local meshpoint [microsec] =', &
            radtransfer_sum_time/icount/nw/1.0e-6
          write(*,'(A,1pG14.7)') &
            ' Diagnostics wall clock time/timestep/local meshpoint [microsec] =', &
            time_doing_diagnostics/icount/nw/1.0e-6
          write(*,'(A,1pG14.7)') &
            ' Timestep wall clock time/timestep/local meshpoint [microsec] =', &
            time_in_timestep/icount/nw/1.0e-6
          write(*,'(A,1pG14.7)') &
            ' Timestep (without waiting) wall clock time/timestep/local meshpoint [microsec] =', &
            (time_in_timestep-time_spent_copying_and_waiting)/icount/nw/1.0e-6
          write(*,'(A,1pG14.7)') &
            ' Timestep wall clock time/timestep/meshpoint [microsec] =', &
            time_in_timestep/icount/nw/ncpus/1.0e-6
          write(*,'(A,1pG14.7)') &
            ' Timestep (without waiting) wall clock time/timestep/meshpoint [microsec] =', &
            (time_in_timestep-time_spent_copying_and_waiting)/icount/nw/ncpus/1.0e-6
          write(*,'(A,1pG14.7)') &
            ' Diagnostics wall clock time/timestep/meshpoint [microsec] =', &
            time_doing_diagnostics/icount/nw/ncpus/1.0e-6
          write(*,'(A,1pG14.7)') &
            ' Rhs wall clock time/timestep/meshpoint [microsec] =', &
            rhs_sum_time/icount/nw/ncpus/1.0e-6
          write(*,'(A,1pG14.7)') &
            ' After substep wall clock time/timestep/meshpoint [microsec] =', &
            after_substep_sum_time/icount/nw/ncpus/1.0e-6
          write(*,'(A,1pG14.7)') &
            ' Before boundary wall clock time/timestep/meshpoint [microsec] =', &
            before_boundary_sum_time/icount/nw/ncpus/1.0e-6
          if(lradiation_ray) write(*,'(A,1pG14.7)') &
            ' Radtransfer wall clock time/timestep/meshpoint [microsec] =', &
            radtransfer_sum_time/icount/nw/ncpus/1.0e-6
          if(lradiation_ray) write(*,'(A,1pG14.7)') &
            ' Radtransfer+Rhs wall clock time/timestep/meshpoint [microsec] =', &
            (rhs_sum_time+radtransfer_sum_time)/icount/nw/ncpus/1.0e-6

          if(ltraining.and.training_time>0) write(*,'(A,1pG14.7)') &
            ' Training wall clock time/timestep/meshpoint [microsec] =', &
            (training_time)/icount/nw/1.0e-6

          if(ltraining.and.inference_time>0) write(*,'(A,1pG14.7)') &
            ' Inference wall clock time/timestep/meshpoint [microsec] =', &
            (inference_time)/icount/nw/1.0e-6


        endif
      endif
    endif
  endif

  memuse=memusage()
  call mpireduce_max_int(memuse,memcpu)
  call mpireduce_sum_int(memuse,memory)
  if (lroot) then
    print'(1x,a,f9.3)', 'Maximum used memory per cpu [MBytes] = ', memcpu/1024.
    if (memory>1e6) then
      print'(1x,a,f12.3)', 'Maximum used memory [GBytes] = ', memory/1024.**2
    else
      print'(1x,a,f12.3)', 'Maximum used memory [MBytes] = ', memory/1024.
    endif
    print*
  endif
!
!  Give all modules the possibility to exit properly.
!
  call finalize_modules(f)
!
!  Write success file, if the requested simulation is complete.
!
  if ((it > nt) .or. (t > tmax)) call touch_file('COMPLETED')
  if (t > tmax) call touch_file('ENDTIME')
!
!  Stop MPI.
!
  call mpifinalize
  call python_finalize
!
!  Free any allocated memory.
!  MR: Is this needed? the program terminates anyway
!
  call diagnostics_clean_up
  call farray_clean_up
  call sharedvars_clean_up
  call chemistry_clean_up
  call NSCBC_clean_up
  if (lparticles) call particles_cleanup
  call finalize
!
  endsubroutine run_start
!***********************************************************************
endmodule Run_module
!***********************************************************************
program run

  use Run_module
!
  call run_start
!
endprogram run
