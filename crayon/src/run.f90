! $Id: run.f90 13590 2010-04-02 10:29:26Z Bourdin.KIS $
!
!***********************************************************************
!
!  The Pencil Code is a high-order finite-difference code for compressible
!  hydrodynamic flows with magnetic fields. It is highly modular and can 
!  easily be adapted to different types of problems.
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
  use Diagnostics
  use Equ,             only: debug_imn_arrays,initialize_pencils
  use EquationOfState, only: ioninit,ioncalc
  use FArrayManager,   only: farray_clean_up
  use IO
  use Messages
  use Mpicomm
  use Param_IO
  use Pencil_check,    only: pencil_consistency_check
  use Register
  use SharedVariables, only: sharedvars_clean_up
  use Signal_handling, only: signal_prepare, emergency_stop
  use Slices
  use Snapshot
  use Sub
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
  integer :: icount, mvar_in
  integer :: it_last_diagnostic, it_this_diagnostic
  logical :: lstop=.false., timeover=.false., resubmit=.false.
  logical :: suppress_pencil_check=.false.
  logical :: lreload_file=.false., lreload_always_file=.false.
!
  lrun=.true.
!
!  Initialize the message subsystem, eg. color setting etc.
!
  call initialize_messages()
!
!  Initialize MPI and register physics modules.
!  (must be done before lroot can be used, for example)
!
  call register_modules()
!
!  Identify version.
!
  if (lroot) call svn_id( &
      '$Id: run.f90 13590 2010-04-02 10:29:26Z Bourdin.KIS $')
!
!  Read parameters from start.x (default values; may be overwritten by
!  read_runpars).
!
  call rparam()
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
!  Read parameters and output parameter list.
!
  call read_runpars()
  call rprint_list(LRESET=.false.)
!
!  Position of equator (if any).
!
  if (lequatory) yequator=xyz0(2)+0.5*Lxyz(2)
  if (lequatorz) zequator=xyz0(3)+0.5*Lxyz(3)
!
! and limits to xaveraging.
!
  if (lav_smallx) call init_xaver
!
!  Will we write all slots of f?
!
  if (lwrite_aux) then
    mvar_io=mvar+maux
  else
    mvar_io=mvar
  endif
!
! Shall we read also auxiliary variables?
!
  if (lread_aux) then
    mvar_in=mvar+maux
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
!  Set up directory names `directory' and `directory_snap'.
!
  call directory_names()
!
!  Get state length of random number generator.
!
  call get_nseed(nseed)
!
!  Read data.
!  Snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!  NOTE: for io_dist, rtime doesn't read the time, only for io_mpio.
!
  call rsnap(trim(directory_snap)//'/var.dat',f,mvar_in)
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
!  Initialize pencils in the pencil_case.
!
  if (lpencil_init) call initialize_pencils(p,0.0)
!
!  Perform pencil_case consistency check if requested.
!
  suppress_pencil_check = control_file_exists("NO-PENCIL-CHECK")
  call mpibcast_logical(suppress_pencil_check, 1)
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
      if (lroot) call fatal_error('run','it1d smaller than it1')
    endif
  endif
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
      call mpibcast_logical(lstop,1)
      if (lstop .or. t>tmax .or. emergency_stop) then
        if (lroot) then
          print*
          if (emergency_stop) print*, 'Emergency stop requested'
          if (lstop) print*, 'Found STOP file'
          if (t>tmax) print*, 'Maximum simulation time exceeded'
          resubmit=control_file_exists('RESUBMIT',DELETE=.true.)
          if (resubmit) then
            print*, 'Cannot be resubmitted'
          else
          endif
        endif
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
      lreload_file       =control_file_exists('RELOAD', DELETE=.true.)
      lreload_always_file=control_file_exists('RELOAD_ALWAYS')
      lreloading         =lreload_file .or. lreload_always_file
!
!  In some compilers (particularly pathf90) the file reload is being give
!  unit = 1 hence there is conflict during re-reading of parameters.
!  In this temporary fix, the RELOAD file is being removed just after it
!  has been seen, not after RELOAD-ing has been completed. There must
!  be a better solution.
!
      call mpibcast_logical(lreloading, 1)
!
      if (lreloading) then
        if (lroot) write(0,*) 'Found RELOAD file -- reloading parameters'
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
        if (lwrite_yaverages)    call yaverages_clean_up()
        if (lwrite_zaverages)    call zaverages_clean_up()
        call rprint_list(LRESET=.true.) !(Re-read output list)
        call initialize_modules(f,LSTARTING=.false.)
        call choose_pencils()
        call wparam2()
        if (lroot .and. lreload_file) call remove_file('RELOAD')
        lreload_file        = .false.
        lreload_always_file = .false.
        lreloading          = .false.
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
!  Time advance.
!
    call rk_2n(f,df,p)
!
!  Print diagnostic averages to screen and file.
!
    if (lout)   call prints()
    if (l1davg) call write_1daverages()
    if (l2davg) call write_2daverages()
!
    if (lroot) icount=icount+1  !  reliable loop count even for premature exit
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
    call wsnap(trim(directory_snap)//'/VAR',f, &
        mvar_io,ENUM=.true.,FLIST='varN.list')
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
        call wsnap(trim(directory_snap)//'/var.dat', &
                   f,mvar_io,ENUM=.false.,noghost=noghost_for_isave)
        call wtime(trim(directory)//'/time.dat',t)
      endif
    endif
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
  call wtime(trim(directory)//'/time.dat',t)
  if (save_lastsnap.and..not.lnowrite) then
    call wsnap(trim(directory_snap)//'/var.dat',f,mvar_io,ENUM=.false.)
!
!  dvar is written for analysis and debugging purposes only.
!
    if (ip<=11 .or. lwrite_dvar) then
      call wsnap(trim(directory)//'/dvar.dat',df,mvar,enum=.false., &
          noghost=.true.)
    endif
!
!  Write crash files before exiting if we haven't written var.dat already
!
  else if (save_lastsnap) then
    call wsnap(trim(directory_snap)//'/crash.dat',f,mvar_io,ENUM=.false.)
    if (ip<=11) call wsnap(trim(directory)//'/dcrash.dat',df,mvar,ENUM=.false.)
  endif
!
!  Print wall clock time and time per step and processor for diagnostic
!  purposes.
!
  if (lroot) then
    wall_clock_time=time2-time1
    print*
    write(*,'(A,1pG10.3,A,1pG8.2,A)') &
        ' Wall clock time [hours] = ', wall_clock_time/3600.0, &
        ' (+/- ', real(mpiwtick())/3600.0, ')'
    if (it>1) then
      write(*,'(A,1pG10.3)') &
           ' Wall clock time/timestep/meshpoint [microsec] =', &
           wall_clock_time/icount/nw/ncpus/1.0e-6
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
  call xyaverages_clean_up()
  call xzaverages_clean_up()
  call yzaverages_clean_up()
  if (lwrite_yaverages)    call yaverages_clean_up()
  if (lwrite_zaverages)    call zaverages_clean_up()
!
endprogram run
