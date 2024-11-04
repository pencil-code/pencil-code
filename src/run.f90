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

contains
!***********************************************************************
subroutine helper_loop(f,p)
!
  use Equ, only: perform_diagnostics
!$ use General, only: signal_wait, signal_send
  use Snapshot, only: perform_powersnap, perform_wsnap_ext, perform_wsnap_down
!
  real, dimension (mx,my,mz,mfarray) :: f
  type (pencil_case) :: p
!
! 7-feb-24/TP: coded
!
!$  do while(lhelper_run)

!$    call signal_wait(lhelper_perf,lhelper_run)
!$    if (lhelper_run .and. lhelperflags(PERF_DIAGS)) then 
        !print*,"doing diag"
        call perform_diagnostics(f,p)
!$    else 
!$      lhelperflags(PERF_DIAGS) = .false.
!$    endif
!$    if (lhelper_run .and. lhelperflags(PERF_WSNAP)) then 
        !print*,"doing wsnap"
        call perform_wsnap_ext(f)
!$    else 
!$      lhelperflags(PERF_WSNAP) = .false.
!$    endif
!$    if (lhelper_run .and. lhelperflags(PERF_WSNAP_DOWN)) then 
        !print*,"doing down"
        call perform_wsnap_down(f)
!$    else 
!$      lhelperflags(PERF_WSNAP_DOWN) = .false.
!$    endif
!$    if (lhelper_run .and. lhelperflags(PERF_POWERSNAP)) then 
        !print*,"doing power"
        call perform_powersnap(f)
!$    else 
!$      lhelperflags(PERF_POWERSNAP) = .false.
!$    endif
!$    call signal_send(lhelper_perf,.false.)

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
    if (ldt) then
      dt=0.0
    else
      dtmp=dt
    endif
    call read_all_run_pars
    if (.not.ldt) dt=dtmp
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
if (lgpu.and.lroot) print*, 'vor copy_farray_from_GPU'
if (lroot) flush(6)
    if (lgpu) call copy_farray_from_GPU(f)
    call initialize_modules(f)
    call initialize_boundcond
    if (lparticles) call particles_initialize_modules(f)
    if (lgpu) then
if (lroot) print*, 'vor reload_GPU_config'
if (lroot) flush(6)
      call reload_GPU_config
      call load_farray_to_GPU(f)
if (lroot) print*, 'nach load_to_GPU'
if (lroot) flush(6)
    endif
    call choose_pencils
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
    use Particles_main,   only: write_snapshot_particles
    use PointMasses,     only: pointmasses_write_snapshot
    use Mpicomm,         only: mpiwtime
    use Sub,             only: control_file_exists
    use Solid_Cells,     only: wsnap_ogrid
    use Timeavg,         only: wsnap_timeavgs
    use Fixed_point,     only: wfixed_points
    use Streamlines,     only: wtracers

    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
    integer :: isave_shift=0, i
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
  use Mpicomm,         only: mpibcast_logical, mpiwtime, MPI_COMM_WORLD, mpibarrier
  use Particles_main,  only: particles_rprint_list, particles_initialize_modules, & 
                             particles_load_balance, particles_stochastic
  use Signal_handling, only: emergency_stop
  use Sub,             only: control_file_exists, calc_scl_factor
  use Testscalar,      only: rescaling_testscalar
  use Testfield,       only: rescaling_testfield
  use TestPerturb,     only: testperturb_begin, testperturb_finalize
  use Timeavg,         only: ltavg, update_timeavgs
  use Timestep,        only: time_step
  use Slices,          only: wvid_prepare
  use Solid_Cells,     only: time_step_ogrid
  use Streamlines,     only: tracers_prepare
  use Snapshot,        only: powersnap_prepare
!$ use OMP_lib
!$ use General, only: signal_send, signal_wait
!!$ use, intrinsic :: iso_fortran_env
!
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (mx,my,mz,mvar) :: df
  type (pencil_case) :: p
!
  logical :: lstop=.false., timeover=.false., resubmit=.false., lreload_file, lreload_always_file, &
             lonemorestep=.false.
  real :: wall_clock_time=0., time_per_step=0.
  real(KIND=rkind8) :: time_this_diagnostic
  integer :: it_this_diagnostic

  Time_loop: do while (it<=nt)
!
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
!  calculate scale factor of the universe
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
!  Time advance.
!
    call time_step(f,df,p)
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
      call mpibcast_logical(timeover,comm=MPI_COMM_WORLD)
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
!$  if (lfarray_copied) then
!$    lhelperflags = lmasterflags
!$    lmasterflags = .false.
!$    call signal_send(lhelper_perf,.true.)
!$    lfarray_copied = .false.
!$  endif

  enddo Time_loop
!$ lmultithread = .false.
!$ call signal_wait(lhelper_perf, .false.)
!$ call signal_send(lhelper_run,.false.)

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
  use Equ,             only: initialize_pencils, debug_imn_arrays
  use FArrayManager,   only: farray_clean_up
  use Farray_alloc
  use General,         only: random_seed_wrapper, touch_file, itoa
!$ use General,        only: signal_init,get_cpu
  use Grid,            only: construct_grid, box_vol, grid_bound_data, set_coorsys_dimmask, &
                             construct_serial_arrays, coarsegrid_interp
  use Gpu,             only: gpu_init, register_gpu, load_farray_to_GPU, initialize_gpu, finalize_gpu
  use HDF5_IO,         only: init_hdf5, initialize_hdf5, wdim
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
  use Snapshot,        only: powersnap, rsnap, powersnap_prepare, wsnap
  use Solid_Cells,     only: wsnap_ogrid
  use Special,         only: initialize_mult_special
  use Sub,             only: control_file_exists, get_nseed
  use Syscalls,        only: memusage
  use Timeavg,         only: wsnap_timeavgs
  use Timestep,        only: initialize_timestep
!
!$ use OMP_lib
!$ use, intrinsic :: iso_c_binding
!!$ use, intrinsic :: iso_fortran_env
!!$ use mt, only: wait_all_thread_pool, push_task, free_thread_pool, depend_on_all, default_task_type
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
  logical :: lexist
  integer, parameter :: num_helpers=1
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
!  Initialize GPU use and make threadpool.
!
  call gpu_init
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
  inquire (file=fproc_bounds, exist=lexist)
  call mpibarrier
!
  cgrid: if (luse_oldgrid .and. lexist) then
    if (ip<=6.and.lroot) print*, 'reading grid coordinates'
    call rgrid('grid.dat')
    call rproc_bounds(fproc_bounds)
    call construct_serial_arrays
    call grid_bound_data
  else cgrid
    if (luse_oldgrid) call warning("run.x", "reconstructing the grid")
    if (luse_xyz1) Lxyz = xyz1-xyz0
    call construct_grid(x,y,z,dx,dy,dz)
    call wgrid("grid.dat", lwrite=.true.)
    call wproc_bounds(fproc_bounds)
  endif cgrid
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
  call register_gpu(f)
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
  if (ip<=12.and.lroot) tvar1=mpiwtime()
  call rsnap('var.dat',f,mvar_in,lread_nogrid)
  if (ip<=12.and.lroot) print*,'rsnap: read snapshot var.dat in ',mpiwtime()-tvar1,' seconds'
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
  if (lparticles) then
    if (ip <= 6 .and. lroot) print *, "reading particle snapshot"
    call read_snapshot_particles
  endif
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
  call initialize_timestep
  call initialize_modules(f)
  call initialize_boundcond
  call initialize_gpu
! Load farray to gpu
  call load_farray_to_GPU(f)
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
!  file with the correct number of variables
!
  call wgrid('grid.dat')
  if (.not.luse_oldgrid .or. lwrite_dim_again) then
    call wdim('dim.dat')
    if (ip<11) print*,'Lz=',Lz
    if (ip<11) print*,'z=',z
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
!  Find out which pencils are needed and write information about required,
!  requested and diagnostic pencils to disc.
!
  call choose_pencils
  call write_pencil_info
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
  if (lspec_start .and. t==tstart) then
    lspec=.true.
  endif
!
!  Save spectrum snapshot.
!
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
      call warning('run',"Pencil consistency check not supported on the GPU. You can consider running it on a CPU-only compilation")
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
  !$omp parallel copyin(fname,fnamex,fnamey,fnamez,fnamer,fnamexy,fnamexz,fnamerz,fname_keep,fname_sound,ncountsz,phiavg_norm)
  !$ core_ids(omp_get_thread_num()+1) = get_cpu()
  !$omp end parallel
  !$ call mpibarrier
!
!$omp parallel num_threads(num_helpers+1) &
!$omp copyin(fname,fnamex,fnamey,fnamez,fnamer,fnamexy,fnamexz,fnamerz,fname_keep,fname_sound,ncountsz,phiavg_norm)
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
!$  do i=1,num_helpers
!$omp barrier
!$    if (omp_get_thread_num() == i) call create_communicators
!$  enddo
!$omp barrier
!$  if (omp_get_thread_num() == 0) then
!
!  Start timing for final timing statistics.
!  Initialize timestep diagnostics during the run (whether used or not,
!  see idiag_timeperstep).
!
      !print*,"Master core: ",get_cpu(), iproc
      if (lroot) then
        icount=0
        it_last_diagnostic=icount
        time1=mpiwtime()
        time_last_diagnostic=time1
      endif
      call timeloop(f,df,p)
!print*, 'nach timeloop', iproc
!flush(6)
!stop
!$  else
!$    call helper_loop(f,p)
!$  endif
!$omp barrier
!$omp end parallel
!print*, 'nach parallel', iproc
!flush(6)
!
  if (lroot) then
!
    time2=mpiwtime()
!
    print*
    print*, 'Simulation finished after ', icount, ' time-steps'
!
!  Write data at end of run for restart.
!
    print*
    write(*,*) 'Writing final snapshot at time t =', t
  endif
!
  if (.not.lnowrite) then

    if (ncoarse>1) then
      call update_ghosts(f)
      if (lcoarse) then
        call coarsegrid_interp(f)
        lcoarse=.false.
      endif
      call update_ghosts(f)
    endif

    if (save_lastsnap) then
      if (lparticles) call write_snapshot_particles(f,ENUM=.false.)
      if (lpointmasses) call pointmasses_write_snapshot('qvar.dat',ENUM=.false.)
      if (lsolid_cells) call wsnap_ogrid('ogvar.dat',ENUM=.false.)
!
      call wsnap('var.dat',f,mvar_io,ENUM=.false.)
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
  call finalize_gpu
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
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=2000
    integer(KIND=ikind8), dimension(n_pars) :: p_par

call copy_addr(ncoarse,p_par(1)) ! int
call copy_addr(lcoarse,p_par(2)) ! int
call copy_addr(lcoarse_mn,p_par(3)) ! int
call copy_addr(tslice,p_par(6)) 
call copy_addr(eps_rkf,p_par(7)) 
call copy_addr(eps_stiff,p_par(8)) 
call copy_addr(eps_rkf0,p_par(9)) 
call copy_addr(dsound,p_par(10)) 
call copy_addr(tsound,p_par(11)) 
call copy_addr(soundeps,p_par(12)) 
call copy_addr(lfix_unit_std,p_par(14)) ! int
call copy_addr(m,p_par(44)) ! int
call copy_addr(n,p_par(45)) ! int
call copy_addr(lfirstpoint,p_par(46)) ! int
call copy_addr(dxyz_2,p_par(47)) ! (nx)
call copy_addr(dxyz_4,p_par(48)) ! (nx)
call copy_addr(dxyz_6,p_par(49)) ! (nx)
call copy_addr(dvol,p_par(50)) ! (nx)
call copy_addr(dvol_glob,p_par(53)) 
call copy_addr(x,p_par(54)) ! (mx)
call copy_addr(dx_1,p_par(55)) ! (mx)
call copy_addr(dx2,p_par(56)) ! (mx)
call copy_addr(dx_tilde,p_par(57)) ! (mx)
call copy_addr(xprim,p_par(58)) ! (mx)
call copy_addr(dvol_x,p_par(59)) ! (mx)
call copy_addr(dvol1_x,p_par(60)) ! (mx)
call copy_addr(y,p_par(61)) ! (my)
call copy_addr(dy_1,p_par(62)) ! (my)
call copy_addr(dy2,p_par(63)) ! (my)
call copy_addr(dy_tilde,p_par(64)) ! (my)
call copy_addr(yprim,p_par(65)) ! (my)
call copy_addr(dvol_y,p_par(66)) ! (my)
call copy_addr(dvol1_y,p_par(67)) ! (my)
call copy_addr(z,p_par(68)) ! (mz)
call copy_addr(dz_1,p_par(69)) ! (mz)
call copy_addr(dz2,p_par(70)) ! (mz)
call copy_addr(dz_tilde,p_par(71)) ! (mz)
call copy_addr(zprim,p_par(72)) ! (mz)
call copy_addr(dvol_z,p_par(73)) ! (mz)
call copy_addr(dvol1_z,p_par(74)) ! (mz)
call copy_addr(dx,p_par(75)) 
call copy_addr(dy,p_par(76)) 
call copy_addr(dz,p_par(77)) 
call copy_addr(dxmin,p_par(78)) 
call copy_addr(dxmax,p_par(79)) 
call copy_addr(xgrid,p_par(80)) ! (nxgrid)
call copy_addr(dx1grid,p_par(81)) ! (nxgrid)
call copy_addr(dxtgrid,p_par(82)) ! (nxgrid)
call copy_addr(ygrid,p_par(83)) ! (nygrid)
call copy_addr(dy1grid,p_par(84)) ! (nygrid)
call copy_addr(dytgrid,p_par(85)) ! (nygrid)
call copy_addr(zgrid,p_par(86)) ! (nzgrid)
call copy_addr(dz1grid,p_par(87)) ! (nzgrid)
call copy_addr(dztgrid,p_par(88)) ! (nzgrid)
call copy_addr(xglobal,p_par(89)) ! (mxgrid)
call copy_addr(yglobal,p_par(90)) ! (mygrid)
call copy_addr(zglobal,p_par(91)) ! (mzgrid)
call copy_addr(kx_nyq,p_par(92)) 
call copy_addr(ky_nyq,p_par(93)) 
call copy_addr(kz_nyq,p_par(94)) 
call copy_addr(lcartesian_coords,p_par(96)) ! int
call copy_addr(lspherical_coords,p_par(97)) ! int
call copy_addr(lcylindrical_coords,p_par(98)) ! int
call copy_addr(lpipe_coords,p_par(99)) ! int
call copy_addr(lsphere_in_a_box,p_par(100)) ! int
call copy_addr(lcylinder_in_a_box,p_par(101)) ! int
call copy_addr(luse_latitude,p_par(102)) ! int
call copy_addr(luse_oldgrid,p_par(103)) ! int
call copy_addr(luse_xyz1,p_par(104)) ! int
call copy_addr(lcylindrical_gravity,p_par(105)) ! int
call copy_addr(luniform_z_mesh_aspect_ratio,p_par(106)) ! int
call copy_addr(lconcurrent,p_par(107)) ! int
call copy_addr(tag_foreign,p_par(108)) ! int
call copy_addr(lforeign,p_par(109)) ! int
call copy_addr(lforeign_comm_nblckg,p_par(110)) ! int
call copy_addr(lyinyang,p_par(111)) ! int
call copy_addr(lyang,p_par(112)) ! int
call copy_addr(lcutoff_corners,p_par(113)) ! int
call copy_addr(iyinyang_intpol_type,p_par(115)) ! int
call copy_addr(nzgrid_eff,p_par(116)) ! int
call copy_addr(yy_biquad_weights,p_par(117)) ! (4)
call copy_addr(nycut,p_par(118)) ! int
call copy_addr(nzcut,p_par(119)) ! int
call copy_addr(rel_dang,p_par(120)) 
call copy_addr(lcubed_sphere,p_par(121)) ! int
call copy_addr(drcyl,p_par(122)) 
call copy_addr(dsurfxy,p_par(123)) 
call copy_addr(dsurfyz,p_par(124)) 
call copy_addr(dsurfzx,p_par(125)) 
call copy_addr(r_mn,p_par(126)) ! (nx)
call copy_addr(r1_mn,p_par(127)) ! (nx)
call copy_addr(r2_mn,p_par(128)) ! (nx)
call copy_addr(r2_weight,p_par(129)) ! (nx)
call copy_addr(sinth,p_par(130)) ! (my)
call copy_addr(sin1th,p_par(131)) ! (my)
call copy_addr(sin2th,p_par(132)) ! (my)
call copy_addr(costh,p_par(133)) ! (my)
call copy_addr(cotth,p_par(134)) ! (my)
call copy_addr(sinth_weight,p_par(135)) ! (my)
call copy_addr(sinph,p_par(136)) ! (mz)
call copy_addr(cosph,p_par(137)) ! (mz)
call copy_addr(cos1th,p_par(138)) ! (my)
call copy_addr(tanth,p_par(139)) ! (my)
call copy_addr(sinth_weight_across_proc,p_par(140)) ! (nygrid)
call copy_addr(rcyl_mn,p_par(141)) ! (nx)
call copy_addr(rcyl_mn1,p_par(142)) ! (nx)
call copy_addr(rcyl_mn2,p_par(143)) ! (nx)
call copy_addr(rcyl_weight,p_par(144)) ! (nx)
call copy_addr(glncrosssec,p_par(145)) ! (nx)
call copy_addr(rcyl,p_par(146)) ! (nrcyl)
call copy_addr(x12,p_par(147)) ! (mx)
call copy_addr(y12,p_par(148)) ! (my)
call copy_addr(z12,p_par(149)) ! (mz)
call copy_addr(coeff_grid,p_par(150)) ! (3)
call copy_addr(dxi_fact,p_par(151)) ! (3)
call copy_addr(trans_width,p_par(152)) ! (3)
call copy_addr(zeta_grid0,p_par(153)) 
call copy_addr(xbot_slice,p_par(154)) 
call copy_addr(xtop_slice,p_par(155)) 
call copy_addr(ybot_slice,p_par(156)) 
call copy_addr(ytop_slice,p_par(157)) 
call copy_addr(zbot_slice,p_par(158)) 
call copy_addr(ztop_slice,p_par(159)) 
call copy_addr(r_rslice,p_par(160)) 
call copy_addr(nth_rslice,p_par(161)) ! int
call copy_addr(nph_rslice,p_par(162)) ! int
call copy_addr(glncrosssec0,p_par(163)) 
call copy_addr(crosssec_x1,p_par(164)) 
call copy_addr(crosssec_x2,p_par(165)) 
call copy_addr(crosssec_w,p_par(166)) 
call copy_addr(lignore_nonequi,p_par(167)) ! int
call copy_addr(lcart_equi,p_par(168)) ! int
call copy_addr(nghost_read_fewer,p_par(171)) ! int
call copy_addr(lall_onesided,p_par(173)) ! int
call copy_addr(lxyz,p_par(175)) ! (3)
call copy_addr(xyz0,p_par(176)) ! (3)
call copy_addr(xyz1,p_par(177)) ! (3)
call copy_addr(xyz_star,p_par(178)) ! (3)
call copy_addr(lxyz_loc,p_par(179)) ! (3)
call copy_addr(xyz0_loc,p_par(180)) ! (3)
call copy_addr(xyz1_loc,p_par(181)) ! (3)
call copy_addr(r_int,p_par(182)) 
call copy_addr(r_ext,p_par(183)) 
call copy_addr(r_int_border,p_par(184)) 
call copy_addr(r_ext_border,p_par(185)) 
call copy_addr(r_ref,p_par(186)) 
call copy_addr(rsmooth,p_par(187)) 
call copy_addr(box_volume,p_par(188)) 
call copy_addr(area_xy,p_par(189)) 
call copy_addr(area_yz,p_par(190)) 
call copy_addr(area_xz,p_par(191)) 
call copy_addr(mu0,p_par(192)) 
call copy_addr(mu01,p_par(193)) 
call copy_addr(lfirst,p_par(194)) ! int
call copy_addr(llast,p_par(195)) ! int
call copy_addr(ldt_paronly,p_par(196)) ! int
call copy_addr(ldt,p_par(197)) ! int
call copy_addr(tmax,p_par(199)) 
call copy_addr(tstart,p_par(200)) 
call copy_addr(max_walltime,p_par(201)) 
call copy_addr(dt_incr,p_par(202)) 
call copy_addr(dt0,p_par(203)) 
call copy_addr(cdt,p_par(204)) 
call copy_addr(cdts,p_par(205)) 
call copy_addr(cdtr,p_par(206)) 
call copy_addr(cdtc,p_par(207)) 
call copy_addr(cdt_poly,p_par(208)) 
call copy_addr(cdtv,p_par(209)) 
call copy_addr(cdtv2,p_par(210)) 
call copy_addr(cdtv3,p_par(211)) 
call copy_addr(cdtsrc,p_par(212)) 
call copy_addr(cdtf,p_par(213)) 
call copy_addr(ddt,p_par(214)) 
call copy_addr(dtinc,p_par(215)) 
call copy_addr(dtdec,p_par(216)) 
call copy_addr(dtmin,p_par(217)) 
call copy_addr(dtmax,p_par(218)) 
call copy_addr(dt_epsi,p_par(219)) 
call copy_addr(dt_ratio,p_par(220)) 
call copy_addr(nu_sts,p_par(221)) 
call copy_addr(permute_sts,p_par(222)) ! int
call copy_addr(ireset_tstart,p_par(223)) ! int
call copy_addr(nt,p_par(224)) ! int
call copy_addr(it,p_par(225)) ! int
call copy_addr(itorder,p_par(226)) ! int
call copy_addr(itsub,p_par(227)) ! int
call copy_addr(it_timing,p_par(228)) ! int
call copy_addr(it_rmv,p_par(229)) ! int
call copy_addr(nproc_comm,p_par(230)) ! int
call copy_addr(ix,p_par(231)) ! int
call copy_addr(iy,p_par(232)) ! int
call copy_addr(iy2,p_par(233)) ! int
call copy_addr(iz,p_par(234)) ! int
call copy_addr(iz2,p_par(235)) ! int
call copy_addr(iz3,p_par(236)) ! int
call copy_addr(iz4,p_par(237)) ! int
call copy_addr(ix_loc,p_par(238)) ! int
call copy_addr(iy_loc,p_par(239)) ! int
call copy_addr(iy2_loc,p_par(240)) ! int
call copy_addr(iz_loc,p_par(241)) ! int
call copy_addr(iz2_loc,p_par(242)) ! int
call copy_addr(iz3_loc,p_par(243)) ! int
call copy_addr(iz4_loc,p_par(244)) ! int
call copy_addr(iproc,p_par(245)) ! int
call copy_addr(ipx,p_par(246)) ! int
call copy_addr(ipy,p_par(247)) ! int
call copy_addr(ipz,p_par(248)) ! int
call copy_addr(iproc_world,p_par(249)) ! int
call copy_addr(ipatch,p_par(250)) ! int
call copy_addr(lprocz_slowest,p_par(251)) ! int
call copy_addr(lzorder,p_par(252)) ! int
call copy_addr(lmorton_curve,p_par(253)) ! int
call copy_addr(xlneigh,p_par(254)) ! int
call copy_addr(ylneigh,p_par(255)) ! int
call copy_addr(zlneigh,p_par(256)) ! int
call copy_addr(xuneigh,p_par(257)) ! int
call copy_addr(yuneigh,p_par(258)) ! int
call copy_addr(zuneigh,p_par(259)) ! int
call copy_addr(poleneigh,p_par(260)) ! int
call copy_addr(nprocx_node,p_par(261)) ! int
call copy_addr(nprocy_node,p_par(262)) ! int
call copy_addr(nprocz_node,p_par(263)) ! int
call copy_addr(num_after_timestep,p_par(264)) ! int
call copy_addr(ighosts_updated,p_par(265)) ! int
call copy_addr(x0,p_par(266)) 
call copy_addr(y0,p_par(267)) 
call copy_addr(z0,p_par(268)) 
call copy_addr(lx,p_par(269)) 
call copy_addr(ly,p_par(270)) 
call copy_addr(lz,p_par(271)) 
call copy_addr(wav1,p_par(272)) 
call copy_addr(wav1z,p_par(273)) 
call copy_addr(lini_t_eq_zero,p_par(274)) ! int
call copy_addr(lini_t_eq_zero_once,p_par(275)) ! int
call copy_addr(reac_chem,p_par(285)) ! (nx)
call copy_addr(reac_dust,p_par(286)) ! (nx)
call copy_addr(trelax_poly,p_par(287)) 
call copy_addr(reac_pchem,p_par(288)) 
call copy_addr(alpha_ts,p_par(289)) ! (5)
call copy_addr(beta_ts,p_par(290)) ! (5)
call copy_addr(dt_beta_ts,p_par(291)) ! (5)
call copy_addr(lfractional_tstep_advance,p_par(292)) ! int
call copy_addr(lfractional_tstep_negative,p_par(293)) ! int
call copy_addr(lmaxadvec_sum,p_par(294)) ! int
call copy_addr(old_cdtv,p_par(295)) ! int
call copy_addr(leps_fixed,p_par(296)) ! int
call copy_addr(lmaximal_cdtv,p_par(297)) ! int
call copy_addr(lmaximal_cdt,p_par(298)) ! int
call copy_addr(llsode,p_par(300)) ! int
call copy_addr(lstep1,p_par(301)) ! int
call copy_addr(lchemonly,p_par(302)) ! int
call copy_addr(lsplit_second,p_par(303)) ! int
call copy_addr(lsnap,p_par(315)) ! int
call copy_addr(lsnap_down,p_par(316)) ! int
call copy_addr(lspec,p_par(317)) ! int
call copy_addr(lspec_start,p_par(318)) ! int
call copy_addr(lspec_at_tplusdt,p_par(319)) ! int
call copy_addr(dsnap,p_par(320)) 
call copy_addr(dsnap_down,p_par(321)) 
call copy_addr(d1davg,p_par(322)) 
call copy_addr(d2davg,p_par(323)) 
call copy_addr(dvid,p_par(324)) 
call copy_addr(dspec,p_par(325)) 
call copy_addr(dtracers,p_par(326)) 
call copy_addr(dfixed_points,p_par(327)) 
call copy_addr(crash_file_dtmin_factor,p_par(328)) 
call copy_addr(farray_smooth_width,p_par(329)) ! int
call copy_addr(isave,p_par(330)) ! int
call copy_addr(ialive,p_par(331)) ! int
call copy_addr(isaveglobal,p_par(332)) ! int
call copy_addr(nv1_capitalvar,p_par(333)) ! int
call copy_addr(lwrite_ts_hdf5,p_par(334)) ! int
call copy_addr(lsave,p_par(335)) ! int
call copy_addr(lread_aux,p_par(336)) ! int
call copy_addr(lwrite_aux,p_par(337)) ! int
call copy_addr(lwrite_dvar,p_par(338)) ! int
call copy_addr(lenforce_maux_check,p_par(339)) ! int
call copy_addr(lwrite_avg1d_binary,p_par(340)) ! int
call copy_addr(lread_oldsnap,p_par(341)) ! int
call copy_addr(lwrite_var_anyway,p_par(342)) ! int
call copy_addr(lwrite_last_powersnap,p_par(343)) ! int
call copy_addr(lwrite_fsum,p_par(344)) ! int
call copy_addr(lread_oldsnap_rho2lnrho,p_par(345)) ! int
call copy_addr(lread_oldsnap_nomag,p_par(346)) ! int
call copy_addr(lread_oldsnap_lnrho2rho,p_par(347)) ! int
call copy_addr(lread_oldsnap_noshear,p_par(348)) ! int
call copy_addr(lread_oldsnap_nohydro,p_par(349)) ! int
call copy_addr(lread_oldsnap_nohydro_nomu5,p_par(350)) ! int
call copy_addr(lread_oldsnap_onlya,p_par(351)) ! int
call copy_addr(lread_oldsnap_mskipvar,p_par(352)) ! int
call copy_addr(lread_oldsnap_nohydro_efield,p_par(353)) ! int
call copy_addr(lread_oldsnap_nohydro_ekfield,p_par(354)) ! int
call copy_addr(ldivu_perp,p_par(355)) ! int
call copy_addr(lread_oldsnap_nopscalar,p_par(356)) ! int
call copy_addr(lread_oldsnap_notestfield,p_par(357)) ! int
call copy_addr(lread_oldsnap_notestflow,p_par(358)) ! int
call copy_addr(lread_oldsnap_notestscalar,p_par(359)) ! int
call copy_addr(lread_oldsnap_noisothmhd,p_par(360)) ! int
call copy_addr(lread_oldsnap_nosink,p_par(361)) ! int
call copy_addr(lnamelist_error,p_par(362)) ! int
call copy_addr(ltolerate_namelist_errors,p_par(363)) ! int
call copy_addr(lparam_nml,p_par(364)) ! int
call copy_addr(lwrite_dim_again,p_par(365)) ! int
call copy_addr(allproc_print,p_par(366)) ! int
call copy_addr(lproc_print,p_par(367)) ! int
call copy_addr(lseparate_persist,p_par(368)) ! int
call copy_addr(ldistribute_persist,p_par(369)) ! int
call copy_addr(lpersist,p_par(370)) ! int
call copy_addr(lomit_add_data,p_par(371)) ! int
call copy_addr(save_lastsnap,p_par(372)) ! int
call copy_addr(noghost_for_isave,p_par(373)) ! int
call copy_addr(ltec,p_par(374)) ! int
call copy_addr(lformat,p_par(375)) ! int
call copy_addr(lread_less,p_par(376)) ! int
call copy_addr(lread_nogrid,p_par(377)) ! int
call copy_addr(lread_global,p_par(378)) ! int
call copy_addr(loutput_varn_at_exact_tsnap,p_par(379)) ! int
call copy_addr(ldirect_access,p_par(380)) ! int
call copy_addr(lread_from_other_prec,p_par(381)) ! int
call copy_addr(ldownsampl,p_par(382)) ! int
call copy_addr(ldownsampling,p_par(383)) ! int
call copy_addr(lrepair_snap,p_par(384)) ! int
call copy_addr(linterpol_on_repair,p_par(385)) ! int
call copy_addr(lastaroth_output,p_par(386)) ! int
call copy_addr(lzaver_on_input,p_par(388)) ! int
call copy_addr(lfatal_num_vector_369,p_par(389)) ! int
call copy_addr(lsmooth_farray,p_par(390)) ! int
call copy_addr(lread_scl_factor_file,p_par(391)) ! int
call copy_addr(lread_scl_factor_file_new,p_par(392)) ! int
call copy_addr(scl_factor_target,p_par(393)) 
call copy_addr(hp_target,p_par(394)) 
call copy_addr(appa_target,p_par(395)) 
call copy_addr(wweos_target,p_par(396)) 
call copy_addr(ip,p_par(397)) ! int
call copy_addr(omega,p_par(398)) 
call copy_addr(theta,p_par(399)) 
call copy_addr(phi,p_par(400)) 
call copy_addr(qshear,p_par(401)) 
call copy_addr(sshear,p_par(402)) 
call copy_addr(deltay,p_par(403)) 
call copy_addr(ldensity_nolog,p_par(404)) ! int
call copy_addr(lreference_state,p_par(405)) ! int
call copy_addr(lfullvar_in_slices,p_par(406)) ! int
call copy_addr(lsubstract_reference_state,p_par(407)) ! int
call copy_addr(ldensity_linearstart,p_par(408)) ! int
call copy_addr(lforcing_cont,p_par(409)) ! int
call copy_addr(lwrite_slices,p_par(410)) ! int
call copy_addr(lwrite_1daverages,p_par(411)) ! int
call copy_addr(lwrite_2daverages,p_par(412)) ! int
call copy_addr(lwrite_tracers,p_par(413)) ! int
call copy_addr(lwrite_fixed_points,p_par(414)) ! int
call copy_addr(lwrite_sound,p_par(415)) ! int
call copy_addr(lwrite_slice_xy2,p_par(416)) ! int
call copy_addr(lwrite_slice_xy,p_par(417)) ! int
call copy_addr(lwrite_slice_xz,p_par(418)) ! int
call copy_addr(lwrite_slice_yz,p_par(419)) ! int
call copy_addr(lwrite_slice_xy3,p_par(420)) ! int
call copy_addr(lwrite_slice_xy4,p_par(421)) ! int
call copy_addr(lwrite_slice_xz2,p_par(422)) ! int
call copy_addr(lwrite_slice_r,p_par(423)) ! int
call copy_addr(lgravx,p_par(424)) ! int
call copy_addr(lgravy,p_par(425)) ! int
call copy_addr(lgravz,p_par(426)) ! int
call copy_addr(lgravx_gas,p_par(427)) ! int
call copy_addr(lgravy_gas,p_par(428)) ! int
call copy_addr(lgravz_gas,p_par(429)) ! int
call copy_addr(lgravx_dust,p_par(430)) ! int
call copy_addr(lgravy_dust,p_par(431)) ! int
call copy_addr(lgravz_dust,p_par(432)) ! int
call copy_addr(lgravr,p_par(433)) ! int
call copy_addr(lwrite_ic,p_par(434)) ! int
call copy_addr(lnowrite,p_par(435)) ! int
call copy_addr(lserial_io,p_par(436)) ! int
call copy_addr(lmodify,p_par(437)) ! int
call copy_addr(lroot,p_par(438)) ! int
call copy_addr(lcaproot,p_par(439)) ! int
call copy_addr(ldebug,p_par(440)) ! int
call copy_addr(lfft,p_par(441)) ! int
call copy_addr(lproc_pt,p_par(442)) ! int
call copy_addr(lproc_p2,p_par(443)) ! int
call copy_addr(lfirst_proc_x,p_par(444)) ! int
call copy_addr(lfirst_proc_y,p_par(445)) ! int
call copy_addr(lfirst_proc_z,p_par(446)) ! int
call copy_addr(lfirst_proc_xy,p_par(447)) ! int
call copy_addr(lfirst_proc_yz,p_par(448)) ! int
call copy_addr(lfirst_proc_xz,p_par(449)) ! int
call copy_addr(lfirst_proc_xyz,p_par(450)) ! int
call copy_addr(llast_proc_x,p_par(451)) ! int
call copy_addr(llast_proc_y,p_par(452)) ! int
call copy_addr(llast_proc_z,p_par(453)) ! int
call copy_addr(llast_proc_xy,p_par(454)) ! int
call copy_addr(llast_proc_yz,p_par(455)) ! int
call copy_addr(llast_proc_xz,p_par(456)) ! int
call copy_addr(llast_proc_xyz,p_par(457)) ! int
call copy_addr(lnorth_pole,p_par(458)) ! int
call copy_addr(lsouth_pole,p_par(459)) ! int
call copy_addr(lpscalar_nolog,p_par(460)) ! int
call copy_addr(lalpm,p_par(461)) ! int
call copy_addr(lalpm_alternate,p_par(462)) ! int
call copy_addr(ldustdensity_log,p_par(463)) ! int
call copy_addr(lmdvar,p_par(464)) ! int
call copy_addr(ldcore,p_par(465)) ! int
call copy_addr(lneutraldensity_nolog,p_par(466)) ! int
call copy_addr(lvisc_smag,p_par(467)) ! int
call copy_addr(lslope_limit_diff,p_par(468)) ! int
call copy_addr(ltemperature_nolog,p_par(469)) ! int
call copy_addr(ltestperturb,p_par(470)) ! int
call copy_addr(lweno_transport,p_par(471)) ! int
call copy_addr(lstart,p_par(472)) ! int
call copy_addr(lrun,p_par(473)) ! int
call copy_addr(lreloading,p_par(474)) ! int
call copy_addr(ladv_der_as_aux,p_par(475)) ! int
call copy_addr(lghostfold_usebspline,p_par(476)) ! int
call copy_addr(lcooling_ss_mz,p_par(477)) ! int
call copy_addr(lshock_heat,p_par(478)) ! int
call copy_addr(density_scale_factor,p_par(479)) 
call copy_addr(pretend_lntt,p_par(480)) ! int
call copy_addr(nvar,p_par(481)) ! int
call copy_addr(naux,p_par(482)) ! int
call copy_addr(naux_com,p_par(483)) ! int
call copy_addr(nscratch,p_par(484)) ! int
call copy_addr(nglobal,p_par(485)) ! int
call copy_addr(n_odevars,p_par(486)) ! int
call copy_addr(lode,p_par(487)) ! int
call copy_addr(ilnrho,p_par(488)) ! int
call copy_addr(irho,p_par(489)) ! int
call copy_addr(irho_b,p_par(490)) ! int
call copy_addr(iss_b,p_par(491)) ! int
call copy_addr(ipp,p_par(492)) ! int
call copy_addr(irhs,p_par(493)) ! int
call copy_addr(ittold,p_par(494)) ! int
call copy_addr(ipoly,p_par(495)) ! int
call copy_addr(ip11,p_par(496)) ! int
call copy_addr(ip12,p_par(497)) ! int
call copy_addr(ip13,p_par(498)) ! int
call copy_addr(ip21,p_par(499)) ! int
call copy_addr(ip22,p_par(500)) ! int
call copy_addr(ip23,p_par(501)) ! int
call copy_addr(ip31,p_par(502)) ! int
call copy_addr(ip32,p_par(503)) ! int
call copy_addr(ip33,p_par(504)) ! int
call copy_addr(ipoly_fr,p_par(505)) ! int
call copy_addr(iuu,p_par(506)) ! int
call copy_addr(iux,p_par(507)) ! int
call copy_addr(iuy,p_par(508)) ! int
call copy_addr(iuz,p_par(509)) ! int
call copy_addr(iss,p_par(510)) ! int
call copy_addr(iphiuu,p_par(511)) ! int
call copy_addr(ilorentz,p_par(512)) ! int
call copy_addr(iuu0,p_par(513)) ! int
call copy_addr(iu0x,p_par(514)) ! int
call copy_addr(iu0y,p_par(515)) ! int
call copy_addr(iu0z,p_par(516)) ! int
call copy_addr(ioo,p_par(517)) ! int
call copy_addr(iox,p_par(518)) ! int
call copy_addr(ioy,p_par(519)) ! int
call copy_addr(ioz,p_par(520)) ! int
call copy_addr(ivv,p_par(521)) ! int
call copy_addr(ivx,p_par(522)) ! int
call copy_addr(ivy,p_par(523)) ! int
call copy_addr(ivz,p_par(524)) ! int
call copy_addr(igradu11,p_par(525)) ! int
call copy_addr(igradu12,p_par(526)) ! int
call copy_addr(igradu13,p_par(527)) ! int
call copy_addr(igradu21,p_par(528)) ! int
call copy_addr(igradu22,p_par(529)) ! int
call copy_addr(igradu23,p_par(530)) ! int
call copy_addr(igradu31,p_par(531)) ! int
call copy_addr(igradu32,p_par(532)) ! int
call copy_addr(igradu33,p_par(533)) ! int
call copy_addr(ispecialvar,p_par(534)) ! int
call copy_addr(ispecialvar2,p_par(535)) ! int
call copy_addr(iuut,p_par(536)) ! int
call copy_addr(iuxt,p_par(537)) ! int
call copy_addr(iuyt,p_par(538)) ! int
call copy_addr(iuzt,p_par(539)) ! int
call copy_addr(ioot,p_par(540)) ! int
call copy_addr(ioxt,p_par(541)) ! int
call copy_addr(ioyt,p_par(542)) ! int
call copy_addr(iozt,p_par(543)) ! int
call copy_addr(iuust,p_par(544)) ! int
call copy_addr(iuxst,p_par(545)) ! int
call copy_addr(iuyst,p_par(546)) ! int
call copy_addr(iuzst,p_par(547)) ! int
call copy_addr(ioost,p_par(548)) ! int
call copy_addr(ioxst,p_par(549)) ! int
call copy_addr(ioyst,p_par(550)) ! int
call copy_addr(iozst,p_par(551)) ! int
call copy_addr(ibbt,p_par(552)) ! int
call copy_addr(ibxt,p_par(553)) ! int
call copy_addr(ibyt,p_par(554)) ! int
call copy_addr(ibzt,p_par(555)) ! int
call copy_addr(ijjt,p_par(556)) ! int
call copy_addr(ijxt,p_par(557)) ! int
call copy_addr(ijyt,p_par(558)) ! int
call copy_addr(ijzt,p_par(559)) ! int
call copy_addr(ijxb,p_par(560)) ! int
call copy_addr(ijxbx,p_par(561)) ! int
call copy_addr(ijxby,p_par(562)) ! int
call copy_addr(ijxbz,p_par(563)) ! int
call copy_addr(iuxb,p_par(564)) ! int
call copy_addr(iuxbx,p_par(565)) ! int
call copy_addr(iuxby,p_par(566)) ! int
call copy_addr(iuxbz,p_par(567)) ! int
call copy_addr(iugb,p_par(568)) ! int
call copy_addr(iugbx,p_par(569)) ! int
call copy_addr(iugby,p_par(570)) ! int
call copy_addr(iugbz,p_par(571)) ! int
call copy_addr(ibgu,p_par(572)) ! int
call copy_addr(ibgux,p_par(573)) ! int
call copy_addr(ibguy,p_par(574)) ! int
call copy_addr(ibguz,p_par(575)) ! int
call copy_addr(ibdivu,p_par(576)) ! int
call copy_addr(ibdivux,p_par(577)) ! int
call copy_addr(ibdivuy,p_par(578)) ! int
call copy_addr(ibdivuz,p_par(579)) ! int
call copy_addr(ibxf,p_par(580)) ! int
call copy_addr(ibyf,p_par(581)) ! int
call copy_addr(ibzf,p_par(582)) ! int
call copy_addr(ibbf,p_par(583)) ! int
call copy_addr(ipotself,p_par(584)) ! int
call copy_addr(iaa,p_par(585)) ! int
call copy_addr(iax,p_par(586)) ! int
call copy_addr(iay,p_par(587)) ! int
call copy_addr(iaz,p_par(588)) ! int
call copy_addr(ispx,p_par(589)) ! int
call copy_addr(ispy,p_par(590)) ! int
call copy_addr(ispz,p_par(591)) ! int
call copy_addr(ifcr,p_par(592)) ! int
call copy_addr(ifcrx,p_par(593)) ! int
call copy_addr(ifcry,p_par(594)) ! int
call copy_addr(ifcrz,p_par(595)) ! int
call copy_addr(ihij,p_par(596)) ! int
call copy_addr(igij,p_par(597)) ! int
call copy_addr(ihht,p_par(598)) ! int
call copy_addr(ihhx,p_par(599)) ! int
call copy_addr(iggt,p_par(600)) ! int
call copy_addr(iggx,p_par(601)) ! int
call copy_addr(istresst,p_par(602)) ! int
call copy_addr(istressx,p_par(603)) ! int
call copy_addr(istress_ij,p_par(604)) ! int
call copy_addr(ihhtim,p_par(605)) ! int
call copy_addr(ihhxim,p_par(606)) ! int
call copy_addr(iggtim,p_par(607)) ! int
call copy_addr(iggxim,p_par(608)) ! int
call copy_addr(istresstim,p_par(609)) ! int
call copy_addr(istressxim,p_par(610)) ! int
call copy_addr(iaatest,p_par(611)) ! int
call copy_addr(iaztestpq,p_par(612)) ! int
call copy_addr(iaxtest,p_par(613)) ! int
call copy_addr(iaytest,p_par(614)) ! int
call copy_addr(iaztest,p_par(615)) ! int
call copy_addr(iuutest,p_par(616)) ! int
call copy_addr(iuztestpq,p_par(617)) ! int
call copy_addr(ihhtestpq,p_par(618)) ! int
call copy_addr(iqx,p_par(619)) ! int
call copy_addr(iqy,p_par(620)) ! int
call copy_addr(iqz,p_par(621)) ! int
call copy_addr(iqq,p_par(622)) ! int
call copy_addr(ntestscalar,p_par(623)) ! int
call copy_addr(ntestfield,p_par(624)) ! int
call copy_addr(ntestflow,p_par(625)) ! int
call copy_addr(ntestlnrho,p_par(626)) ! int
call copy_addr(icctest,p_par(627)) ! int
call copy_addr(icctestpq,p_par(628)) ! int
call copy_addr(iug,p_par(629)) ! int
call copy_addr(iam,p_par(630)) ! int
call copy_addr(iamx,p_par(631)) ! int
call copy_addr(iamy,p_par(632)) ! int
call copy_addr(iamz,p_par(633)) ! int
call copy_addr(ivisc_heat,p_par(634)) ! int
call copy_addr(ibb,p_par(635)) ! int
call copy_addr(ibx,p_par(636)) ! int
call copy_addr(iby,p_par(637)) ! int
call copy_addr(ibz,p_par(638)) ! int
call copy_addr(ijj,p_par(639)) ! int
call copy_addr(ijx,p_par(640)) ! int
call copy_addr(ijy,p_par(641)) ! int
call copy_addr(ijz,p_par(642)) ! int
call copy_addr(ibb_sph,p_par(643)) ! int
call copy_addr(ibb_sphr,p_par(644)) ! int
call copy_addr(ibb_spht,p_par(645)) ! int
call copy_addr(ibb_sphp,p_par(646)) ! int
call copy_addr(inusmag,p_par(647)) ! int
call copy_addr(ietasmag,p_par(648)) ! int
call copy_addr(iaak,p_par(649)) ! int
call copy_addr(iaakim,p_par(650)) ! int
call copy_addr(ieek,p_par(651)) ! int
call copy_addr(ieekim,p_par(652)) ! int
call copy_addr(iee,p_par(653)) ! int
call copy_addr(iex,p_par(654)) ! int
call copy_addr(iey,p_par(655)) ! int
call copy_addr(iez,p_par(656)) ! int
call copy_addr(ialfven,p_par(657)) ! int
call copy_addr(iff_diff,p_par(658)) ! int
call copy_addr(iff_diff1,p_par(659)) ! int
call copy_addr(iff_diff2,p_par(660)) ! int
call copy_addr(iff_div_uu,p_par(661)) ! int
call copy_addr(iff_div_aa,p_par(662)) ! int
call copy_addr(iff_div_ss,p_par(663)) ! int
call copy_addr(iff_div_rho,p_par(664)) ! int
call copy_addr(iff_char_c,p_par(665)) ! int
call copy_addr(iff_heat,p_par(666)) ! int
call copy_addr(isld_char,p_par(667)) ! int
call copy_addr(ivisc_forc,p_par(668)) ! int
call copy_addr(ivisc_forcx,p_par(669)) ! int
call copy_addr(ivisc_forcy,p_par(670)) ! int
call copy_addr(ivisc_forcz,p_par(671)) ! int
call copy_addr(i_adv_der,p_par(672)) ! int
call copy_addr(i_adv_derx,p_par(673)) ! int
call copy_addr(i_adv_dery,p_par(674)) ! int
call copy_addr(i_adv_derz,p_par(675)) ! int
call copy_addr(iuxbtest,p_par(676)) ! int
call copy_addr(ijxbtest,p_par(677)) ! int
call copy_addr(iugutest,p_par(678)) ! int
call copy_addr(iughtest,p_par(679)) ! int
call copy_addr(isghtest,p_par(680)) ! int
call copy_addr(ishock,p_par(681)) ! int
call copy_addr(ishock_perp,p_par(682)) ! int
call copy_addr(iyh,p_par(683)) ! int
call copy_addr(ihypvis,p_par(684)) ! int
call copy_addr(ihypres,p_par(685)) ! int
call copy_addr(iecr,p_par(686)) ! int
call copy_addr(ismagorinsky,p_par(687)) ! int
call copy_addr(iviscosity,p_par(688)) ! int
call copy_addr(iqrad,p_par(689)) ! int
call copy_addr(israd,p_par(690)) ! int
call copy_addr(ilntt,p_par(691)) ! int
call copy_addr(itt,p_par(692)) ! int
call copy_addr(ikapparho,p_par(693)) ! int
call copy_addr(ikr_frad,p_par(694)) ! int
call copy_addr(ikr_fradx,p_par(695)) ! int
call copy_addr(ikr_frady,p_par(696)) ! int
call copy_addr(ikr_fradz,p_par(697)) ! int
call copy_addr(igpotselfx,p_par(698)) ! int
call copy_addr(igpotselfy,p_par(699)) ! int
call copy_addr(igpotselfz,p_par(700)) ! int
call copy_addr(icc,p_par(701)) ! int
call copy_addr(ilncc,p_par(702)) ! int
call copy_addr(ialpm,p_par(703)) ! int
call copy_addr(ietat,p_par(704)) ! int
call copy_addr(iacc,p_par(705)) ! int
call copy_addr(issat,p_par(706)) ! int
call copy_addr(ittc,p_par(707)) ! int
call copy_addr(itauascalar,p_par(708)) ! int
call copy_addr(iaphi,p_par(709)) ! int
call copy_addr(ibphi,p_par(710)) ! int
call copy_addr(ieth,p_par(711)) ! int
call copy_addr(idet,p_par(712)) ! int
call copy_addr(iinvgrid,p_par(713)) ! int
call copy_addr(iguij,p_par(714)) ! int
call copy_addr(igu11,p_par(715)) ! int
call copy_addr(igu12,p_par(716)) ! int
call copy_addr(igu13,p_par(717)) ! int
call copy_addr(igu21,p_par(718)) ! int
call copy_addr(igu22,p_par(719)) ! int
call copy_addr(igu23,p_par(720)) ! int
call copy_addr(igu31,p_par(721)) ! int
call copy_addr(igu32,p_par(722)) ! int
call copy_addr(igu33,p_par(723)) ! int
call copy_addr(icooling,p_par(724)) ! int
call copy_addr(inetheat,p_par(725)) ! int
call copy_addr(ilnrhon,p_par(726)) ! int
call copy_addr(irhon,p_par(727)) ! int
call copy_addr(irhoe,p_par(728)) ! int
call copy_addr(iuun,p_par(729)) ! int
call copy_addr(iunx,p_par(730)) ! int
call copy_addr(iuny,p_par(731)) ! int
call copy_addr(iunz,p_par(732)) ! int
call copy_addr(iglobal_bx_ext,p_par(733)) ! int
call copy_addr(iglobal_by_ext,p_par(734)) ! int
call copy_addr(iglobal_bz_ext,p_par(735)) ! int
call copy_addr(iglobal_ax_ext,p_par(736)) ! int
call copy_addr(iglobal_ay_ext,p_par(737)) ! int
call copy_addr(iglobal_az_ext,p_par(738)) ! int
call copy_addr(iglobal_lnrho0,p_par(739)) ! int
call copy_addr(iglobal_ss0,p_par(740)) ! int
call copy_addr(icp,p_par(741)) ! int
call copy_addr(igpx,p_par(742)) ! int
call copy_addr(igpy,p_par(743)) ! int
call copy_addr(irr,p_par(744)) ! int
call copy_addr(iss_run_aver,p_par(745)) ! int
call copy_addr(ifenth,p_par(746)) ! int
call copy_addr(iss_flucz,p_par(747)) ! int
call copy_addr(itt_flucz,p_par(748)) ! int
call copy_addr(irho_flucz,p_par(749)) ! int
call copy_addr(iuu_fluc,p_par(750)) ! int
call copy_addr(iuu_flucx,p_par(751)) ! int
call copy_addr(iuu_flucy,p_par(752)) ! int
call copy_addr(iuu_flucz,p_par(753)) ! int
call copy_addr(iuu_sph,p_par(754)) ! int
call copy_addr(iuu_sphr,p_par(755)) ! int
call copy_addr(iuu_spht,p_par(756)) ! int
call copy_addr(iuu_sphp,p_par(757)) ! int
call copy_addr(ics,p_par(758)) ! int
call copy_addr(imn,p_par(759)) ! int
call copy_addr(lglob,p_par(760)) ! int
call copy_addr(necessary_imn,p_par(761)) ! int
call copy_addr(penc0,p_par(762)) 
call copy_addr(lpencil_check,p_par(763)) ! int
call copy_addr(lpencil_check_small,p_par(764)) ! int
call copy_addr(lpencil_check_no_zeros,p_par(765)) ! int
call copy_addr(lpencil_init,p_par(766)) ! int
call copy_addr(lpencil_requested_swap,p_par(767)) ! int
call copy_addr(lpencil_diagnos_swap,p_par(768)) ! int
call copy_addr(lpencil_check_diagnos_opti,p_par(769)) ! int
call copy_addr(lpencil_check_at_work,p_par(770)) ! int
call copy_addr(ipencil_swap,p_par(771)) ! int
call copy_addr(it1,p_par(773)) ! int
call copy_addr(it1start,p_par(774)) ! int
call copy_addr(it1d,p_par(775)) ! int
call copy_addr(itspec,p_par(776)) ! int
call copy_addr(nname,p_par(777)) ! int
call copy_addr(nnamev,p_par(778)) ! int
call copy_addr(nnamexy,p_par(779)) ! int
call copy_addr(nnamexz,p_par(780)) ! int
call copy_addr(nnamerz,p_par(781)) ! int
call copy_addr(nnamez,p_par(782)) ! int
call copy_addr(nnamey,p_par(783)) ! int
call copy_addr(nnamex,p_par(784)) ! int
call copy_addr(nnamer,p_par(785)) ! int
call copy_addr(nname_sound,p_par(786)) ! int
call copy_addr(ncoords_sound,p_par(787)) ! int
call copy_addr(nr_directions,p_par(788)) ! int
call copy_addr(itdiagnos,p_par(789)) ! int
call copy_addr(tdiagnos,p_par(790)) 
call copy_addr(dtdiagnos,p_par(791)) 
call copy_addr(t1ddiagnos,p_par(792)) 
call copy_addr(t2davgfirst,p_par(793)) 
call copy_addr(eps_rkf_diagnos,p_par(794)) 
call copy_addr(fweight,p_par(795)) ! (mname)
call copy_addr(lout,p_par(816)) ! int
call copy_addr(headt,p_par(817)) ! int
call copy_addr(headtt,p_par(818)) ! int
call copy_addr(lrmv,p_par(819)) ! int
call copy_addr(ldiagnos,p_par(820)) ! int
call copy_addr(lvideo,p_par(821)) ! int
call copy_addr(lwrite_prof,p_par(822)) ! int
call copy_addr(lout_sound,p_par(823)) ! int
call copy_addr(ltracers,p_par(824)) ! int
call copy_addr(lfixed_points,p_par(825)) ! int
call copy_addr(l2davg,p_par(826)) ! int
call copy_addr(l2davgfirst,p_par(827)) ! int
call copy_addr(l1davg,p_par(828)) ! int
call copy_addr(l1davgfirst,p_par(829)) ! int
call copy_addr(l1dphiavg,p_par(830)) ! int
call copy_addr(lwrite_xyaverages,p_par(831)) ! int
call copy_addr(lwrite_xzaverages,p_par(832)) ! int
call copy_addr(lwrite_yzaverages,p_par(833)) ! int
call copy_addr(lwrite_phizaverages,p_par(834)) ! int
call copy_addr(lwrite_yaverages,p_par(835)) ! int
call copy_addr(lwrite_zaverages,p_par(836)) ! int
call copy_addr(lwrite_phiaverages,p_par(837)) ! int
call copy_addr(ldiagnos_need_zaverages,p_par(838)) ! int
call copy_addr(ltime_integrals,p_par(839)) ! int
call copy_addr(lreset_seed,p_par(840)) ! int
call copy_addr(lproper_averages,p_par(841)) ! int
call copy_addr(lav_smallx,p_par(843)) ! int
call copy_addr(loutside_avg,p_par(844)) ! int
call copy_addr(xav_max,p_par(845)) 
call copy_addr(nvol,p_par(846)) 
call copy_addr(nvol1,p_par(847)) 
call copy_addr(nseed,p_par(848)) ! int
call copy_addr(seed0,p_par(849)) ! int
call copy_addr(ichannel1,p_par(850)) ! int
call copy_addr(ichannel2,p_par(851)) ! int
call copy_addr(fran1,p_par(852)) ! (2)
call copy_addr(fran2,p_par(853)) ! (2)
call copy_addr(lseed_global,p_par(854)) ! int
call copy_addr(lseed_procdependent,p_par(855)) ! int
call copy_addr(yequator,p_par(856)) 
call copy_addr(zequator,p_par(857)) 
call copy_addr(lequatory,p_par(858)) ! int
call copy_addr(lequatorz,p_par(859)) ! int
call copy_addr(name_half_max,p_par(860)) ! int
call copy_addr(radius_diag,p_par(862)) 
call copy_addr(lpoint,p_par(863)) ! int
call copy_addr(mpoint,p_par(864)) ! int
call copy_addr(npoint,p_par(865)) ! int
call copy_addr(lpoint2,p_par(866)) ! int
call copy_addr(mpoint2,p_par(867)) ! int
call copy_addr(npoint2,p_par(868)) ! int
call copy_addr(iproc_pt,p_par(869)) ! int
call copy_addr(iproc_p2,p_par(870)) ! int
call copy_addr(idiag_it,p_par(871)) ! int
call copy_addr(idiag_t,p_par(872)) ! int
call copy_addr(idiag_dt,p_par(873)) ! int
call copy_addr(idiag_walltime,p_par(874)) ! int
call copy_addr(idiag_timeperstep,p_par(875)) ! int
call copy_addr(idiag_rcylmphi,p_par(876)) ! int
call copy_addr(idiag_phimphi,p_par(877)) ! int
call copy_addr(idiag_zmphi,p_par(878)) ! int
call copy_addr(idiag_rmphi,p_par(879)) ! int
call copy_addr(idiag_dtv,p_par(880)) ! int
call copy_addr(idiag_dtdiffus,p_par(881)) ! int
call copy_addr(idiag_dtdiffus2,p_par(882)) ! int
call copy_addr(idiag_dtdiffus3,p_par(883)) ! int
call copy_addr(idiag_rmesh,p_par(884)) ! int
call copy_addr(idiag_rmesh3,p_par(885)) ! int
call copy_addr(idiag_maxadvec,p_par(886)) ! int
call copy_addr(idiag_eps_rkf,p_par(887)) ! int
call copy_addr(lemergency_brake,p_par(888)) ! int
call copy_addr(lcopysnapshots_exp,p_par(889)) ! int
call copy_addr(lwrite_2d,p_par(890)) ! int
call copy_addr(lbidiagonal_derij,p_par(891)) ! int
call copy_addr(vel_spec,p_par(892)) ! int
call copy_addr(mag_spec,p_par(893)) ! int
call copy_addr(uxj_spec,p_par(894)) ! int
call copy_addr(vec_spec,p_par(895)) ! int
call copy_addr(j_spec,p_par(896)) ! int
call copy_addr(jb_spec,p_par(897)) ! int
call copy_addr(ja_spec,p_par(898)) ! int
call copy_addr(oo_spec,p_par(899)) ! int
call copy_addr(relvel_spec,p_par(900)) ! int
call copy_addr(vel_phispec,p_par(901)) ! int
call copy_addr(mag_phispec,p_par(902)) ! int
call copy_addr(uxj_phispec,p_par(903)) ! int
call copy_addr(vec_phispec,p_par(904)) ! int
call copy_addr(uxy_spec,p_par(905)) ! int
call copy_addr(bxy_spec,p_par(906)) ! int
call copy_addr(jxbxy_spec,p_par(907)) ! int
call copy_addr(ep_spec,p_par(910)) ! int
call copy_addr(nd_spec,p_par(911)) ! int
call copy_addr(ud_spec,p_par(912)) ! int
call copy_addr(abs_u_spec,p_par(913)) ! int
call copy_addr(ro_spec,p_par(914)) ! int
call copy_addr(tt_spec,p_par(915)) ! int
call copy_addr(ss_spec,p_par(916)) ! int
call copy_addr(cc_spec,p_par(917)) ! int
call copy_addr(cr_spec,p_par(918)) ! int
call copy_addr(sp_spec,p_par(919)) ! int
call copy_addr(ssp_spec,p_par(920)) ! int
call copy_addr(sssp_spec,p_par(921)) ! int
call copy_addr(mu_spec,p_par(922)) ! int
call copy_addr(lr_spec,p_par(923)) ! int
call copy_addr(r2u_spec,p_par(924)) ! int
call copy_addr(r3u_spec,p_par(925)) ! int
call copy_addr(oun_spec,p_par(926)) ! int
call copy_addr(np_spec,p_par(927)) ! int
call copy_addr(np_ap_spec,p_par(928)) ! int
call copy_addr(rhop_spec,p_par(929)) ! int
call copy_addr(ele_spec,p_par(930)) ! int
call copy_addr(pot_spec,p_par(931)) ! int
call copy_addr(ux_spec,p_par(932)) ! int
call copy_addr(uy_spec,p_par(933)) ! int
call copy_addr(uz_spec,p_par(934)) ! int
call copy_addr(a0_spec,p_par(935)) ! int
call copy_addr(ucp_spec,p_par(936)) ! int
call copy_addr(ou_spec,p_par(937)) ! int
call copy_addr(ab_spec,p_par(938)) ! int
call copy_addr(azbz_spec,p_par(939)) ! int
call copy_addr(uzs_spec,p_par(940)) ! int
call copy_addr(ub_spec,p_par(941)) ! int
call copy_addr(lor_spec,p_par(942)) ! int
call copy_addr(emf_spec,p_par(943)) ! int
call copy_addr(tra_spec,p_par(944)) ! int
call copy_addr(gws_spec,p_par(945)) ! int
call copy_addr(gwh_spec,p_par(946)) ! int
call copy_addr(gwm_spec,p_par(947)) ! int
call copy_addr(str_spec,p_par(948)) ! int
call copy_addr(stg_spec,p_par(949)) ! int
call copy_addr(gws_spec_boost,p_par(950)) ! int
call copy_addr(gwh_spec_boost,p_par(951)) ! int
call copy_addr(stt_spec,p_par(952)) ! int
call copy_addr(stx_spec,p_par(953)) ! int
call copy_addr(gwd_spec,p_par(954)) ! int
call copy_addr(gwe_spec,p_par(955)) ! int
call copy_addr(gwf_spec,p_par(956)) ! int
call copy_addr(gwg_spec,p_par(957)) ! int
call copy_addr(scl_spec,p_par(958)) ! int
call copy_addr(vct_spec,p_par(959)) ! int
call copy_addr(tpq_spec,p_par(960)) ! int
call copy_addr(tgw_spec,p_par(961)) ! int
call copy_addr(scl_spec_boost,p_par(962)) ! int
call copy_addr(vct_spec_boost,p_par(963)) ! int
call copy_addr(har_spec,p_par(964)) ! int
call copy_addr(hav_spec,p_par(965)) ! int
call copy_addr(bb2_spec,p_par(966)) ! int
call copy_addr(jj2_spec,p_par(967)) ! int
call copy_addr(b2_spec,p_par(968)) ! int
call copy_addr(oned,p_par(969)) ! int
call copy_addr(twod,p_par(970)) ! int
call copy_addr(ab_phispec,p_par(971)) ! int
call copy_addr(ou_phispec,p_par(972)) ! int
call copy_addr(rhocc_pdf,p_par(973)) ! int
call copy_addr(cc_pdf,p_par(974)) ! int
call copy_addr(lncc_pdf,p_par(975)) ! int
call copy_addr(gcc_pdf,p_par(976)) ! int
call copy_addr(lngcc_pdf,p_par(977)) ! int
call copy_addr(lnspecial_pdf,p_par(978)) ! int
call copy_addr(special_pdf,p_par(979)) ! int
call copy_addr(ang_jb_pdf1d,p_par(980)) ! int
call copy_addr(ang_ub_pdf1d,p_par(981)) ! int
call copy_addr(ang_ou_pdf1d,p_par(982)) ! int
call copy_addr(test_nonblocking,p_par(983)) ! int
call copy_addr(onedall,p_par(984)) ! int
call copy_addr(lsfu,p_par(985)) ! int
call copy_addr(lsfb,p_par(986)) ! int
call copy_addr(lsfz1,p_par(987)) ! int
call copy_addr(lsfz2,p_par(988)) ! int
call copy_addr(lsfflux,p_par(989)) ! int
call copy_addr(lpdfu,p_par(990)) ! int
call copy_addr(lpdfb,p_par(991)) ! int
call copy_addr(lpdfz1,p_par(992)) ! int
call copy_addr(lpdfz2,p_par(993)) ! int
call copy_addr(ou_omega,p_par(994)) ! int
call copy_addr(cor_uu,p_par(995)) ! int
call copy_addr(ab_kzspec,p_par(996)) ! int
call copy_addr(ou_kzspec,p_par(997)) ! int
call copy_addr(ou_polar,p_par(998)) ! int
call copy_addr(ab_polar,p_par(999)) ! int
call copy_addr(jb_polar,p_par(1000)) ! int
call copy_addr(uut_spec,p_par(1001)) ! int
call copy_addr(uut_polar,p_par(1002)) ! int
call copy_addr(ouout_spec,p_par(1003)) ! int
call copy_addr(ouout2_spec,p_par(1004)) ! int
call copy_addr(ouout_polar,p_par(1005)) ! int
call copy_addr(out_spec,p_par(1006)) ! int
call copy_addr(uot_spec,p_par(1007)) ! int
call copy_addr(saffman_ub,p_par(1008)) ! int
call copy_addr(saffman_mag,p_par(1009)) ! int
call copy_addr(saffman_mag_c,p_par(1010)) ! int
call copy_addr(saffman_aa,p_par(1011)) ! int
call copy_addr(saffman_aa_c,p_par(1012)) ! int
call copy_addr(saffman_bb,p_par(1013)) ! int
call copy_addr(uu_fft3d,p_par(1014)) ! int
call copy_addr(oo_fft3d,p_par(1015)) ! int
call copy_addr(bb_fft3d,p_par(1016)) ! int
call copy_addr(jj_fft3d,p_par(1017)) ! int
call copy_addr(uu_xkyz,p_par(1018)) ! int
call copy_addr(oo_xkyz,p_par(1019)) ! int
call copy_addr(bb_xkyz,p_par(1020)) ! int
call copy_addr(jj_xkyz,p_par(1021)) ! int
call copy_addr(uu_kx0z,p_par(1022)) ! int
call copy_addr(oo_kx0z,p_par(1023)) ! int
call copy_addr(bb_kx0z,p_par(1024)) ! int
call copy_addr(jj_kx0z,p_par(1025)) ! int
call copy_addr(bb_k00z,p_par(1026)) ! int
call copy_addr(ee_k00z,p_par(1027)) ! int
call copy_addr(gwt_fft3d,p_par(1028)) ! int
call copy_addr(em_specflux,p_par(1029)) ! int
call copy_addr(hm_specflux,p_par(1030)) ! int
call copy_addr(hc_specflux,p_par(1031)) ! int
call copy_addr(fbcx_bot,p_par(1032)) ! (mcom)
call copy_addr(fbcx_top,p_par(1033)) ! (mcom)
call copy_addr(fbcy_bot,p_par(1034)) ! (mcom)
call copy_addr(fbcy_top,p_par(1035)) ! (mcom)
call copy_addr(fbcz_bot,p_par(1036)) ! (mcom)
call copy_addr(fbcz_top,p_par(1037)) ! (mcom)
call copy_addr(lreset_boundary_values,p_par(1038)) ! int
call copy_addr(udrift_bc,p_par(1039)) 
call copy_addr(xfreeze_square,p_par(1049)) 
call copy_addr(yfreeze_square,p_par(1050)) 
call copy_addr(rfreeze_int,p_par(1051)) 
call copy_addr(rfreeze_ext,p_par(1052)) 
call copy_addr(wfreeze,p_par(1053)) 
call copy_addr(wfreeze_int,p_par(1054)) 
call copy_addr(wfreeze_ext,p_par(1055)) 
call copy_addr(wborder,p_par(1056)) 
call copy_addr(wborder_int,p_par(1057)) 
call copy_addr(wborder_ext,p_par(1058)) 
call copy_addr(tborder,p_par(1059)) 
call copy_addr(fshift_int,p_par(1060)) 
call copy_addr(fshift_ext,p_par(1061)) 
call copy_addr(theta_lower_border,p_par(1062)) 
call copy_addr(wborder_theta_lower,p_par(1063)) 
call copy_addr(theta_upper_border,p_par(1064)) 
call copy_addr(wborder_theta_upper,p_par(1065)) 
call copy_addr(fraction_tborder,p_par(1066)) 
call copy_addr(lmeridional_border_drive,p_par(1067)) ! int
call copy_addr(border_frac_x,p_par(1068)) ! (2)
call copy_addr(border_frac_y,p_par(1069)) ! (2)
call copy_addr(border_frac_z,p_par(1070)) ! (2)
call copy_addr(border_frac_r,p_par(1071)) ! (2)
call copy_addr(lborder_hyper_diff,p_par(1072)) ! int
call copy_addr(lfrozen_bcs_x,p_par(1073)) ! int
call copy_addr(lfrozen_bcs_y,p_par(1074)) ! int
call copy_addr(lfrozen_bcs_z,p_par(1075)) ! int
call copy_addr(lstop_on_ioerror,p_par(1081)) ! int
call copy_addr(aux_count,p_par(1083)) ! int
call copy_addr(mvar_io,p_par(1084)) ! int
call copy_addr(mvar_down,p_par(1085)) ! int
call copy_addr(maux_down,p_par(1086)) ! int
call copy_addr(mskipvar,p_par(1087)) ! int
call copy_addr(iwig,p_par(1088)) ! int
call copy_addr(nfilter,p_par(1089)) ! int
call copy_addr(awig,p_par(1090)) 
call copy_addr(lrmwig_rho,p_par(1091)) ! int
call copy_addr(lrmwig_full,p_par(1092)) ! int
call copy_addr(lrmwig_xyaverage,p_par(1093)) ! int
call copy_addr(init_loops,p_par(1094)) ! int
call copy_addr(lfold_df,p_par(1095)) ! int
call copy_addr(lfold_df_3points,p_par(1096)) ! int
call copy_addr(lshift_datacube_x,p_par(1097)) ! int
call copy_addr(lkinflow_as_aux,p_par(1099)) ! int
call copy_addr(ampl_kinflow_x,p_par(1100)) 
call copy_addr(ampl_kinflow_y,p_par(1101)) 
call copy_addr(ampl_kinflow_z,p_par(1102)) 
call copy_addr(kx_kinflow,p_par(1103)) 
call copy_addr(ky_kinflow,p_par(1104)) 
call copy_addr(kz_kinflow,p_par(1105)) 
call copy_addr(dtphase_kinflow,p_par(1106)) 
call copy_addr(lfargo_advection,p_par(1107)) ! int
call copy_addr(lcorotational_frame,p_par(1108)) ! int
call copy_addr(rcorot,p_par(1109)) 
call copy_addr(omega_corot,p_par(1110)) 
call copy_addr(llocal_iso,p_par(1111)) ! int
call copy_addr(lisotropic_advection,p_par(1112)) ! int
call copy_addr(lreport_undefined_diagnostics,p_par(1113)) ! int
call copy_addr(ttransient,p_par(1114)) 
call copy_addr(b_ell,p_par(1115)) 
call copy_addr(rbound,p_par(1116)) 
call copy_addr(grads0,p_par(1117)) 
call copy_addr(lmonolithic_io,p_par(1118)) ! int
call copy_addr(lrescaling_magnetic,p_par(1119)) ! int
call copy_addr(lrescaling_testscalar,p_par(1120)) ! int
call copy_addr(lrescaling_testfield,p_par(1121)) ! int
call copy_addr(re_mesh,p_par(1122)) 
call copy_addr(ldynamical_diffusion,p_par(1123)) ! int
call copy_addr(ldyndiff_useumax,p_par(1124)) ! int
call copy_addr(lstratz,p_par(1125)) ! int
call copy_addr(lnoghost_strati,p_par(1126)) ! int
call copy_addr(tau_aver1,p_par(1127)) 
call copy_addr(lambda5,p_par(1128)) 
call copy_addr(lmultithread,p_par(1129)) ! int
call copy_addr(l1dphiavg_save,p_par(1130)) ! int
call copy_addr(l1davgfirst_save,p_par(1131)) ! int
call copy_addr(ldiagnos_save,p_par(1132)) ! int
call copy_addr(l2davgfirst_save,p_par(1133)) ! int
call copy_addr(lout_save,p_par(1134)) ! int
call copy_addr(l1davg_save,p_par(1135)) ! int
call copy_addr(l2davg_save,p_par(1136)) ! int
call copy_addr(lout_sound_save,p_par(1137)) ! int
call copy_addr(lvideo_save,p_par(1138)) ! int
call copy_addr(lchemistry_diag_save,p_par(1139)) ! int
call copy_addr(t1ddiagnos_save,p_par(1141)) 
call copy_addr(t2davgfirst_save,p_par(1142)) 
call copy_addr(tslice_save,p_par(1143)) 
call copy_addr(tsound_save,p_par(1144)) 
call copy_addr(num_helper_threads,p_par(1145)) ! int
call copy_addr(thread_id,p_par(1146)) ! int
call copy_addr(dt,p_par(1147))
call copy_addr(l2,p_par(1148))  ! int 
call copy_addr(m2,p_par(1149))  ! int 
call copy_addr(n2,p_par(1150))  ! int 
call copy_addr(l2i,p_par(1151)) ! int
call copy_addr(m2i,p_par(1152)) ! int
call copy_addr(n2i,p_par(1153)) ! int
call copy_addr(ltest_bcs,p_par(1154)) !bool
call copy_addr(fbcx,p_par(1155)) ! (mcom) (2)
call copy_addr(fbcy,p_par(1156)) ! (mcom) (2)
call copy_addr(fbcz,p_par(1157)) ! (mcom) (2)

call copy_addr(fbcx_1,p_par(1158)) ! (mcom) (2)
call copy_addr(fbcy_1,p_par(1159)) ! (mcom) (2)
call copy_addr(fbcz_1,p_par(1160)) ! (mcom) (2)

call copy_addr(fbcx_2,p_par(1161)) ! (mcom) (2)
call copy_addr(fbcy_2,p_par(1162)) ! (mcom) (2)
call copy_addr(fbcz_2,p_par(1163)) ! (mcom) (2)
    endsubroutine pushpars2c
!***********************************************************************
endmodule Run_module
!***********************************************************************
program run

  use Run_module
!
  call run_start
!
endprogram run
!***********************************************************************
