! $Id$
!
!***********************************************************************
      program run
!
!  solve the isothermal hydro equation in 3-D
!
!-----------------------------------------------------------------------
!   1-apr-01/axel+wolf: coded
!  17-aug-01/axel: ghost layers implemented
!  11-sep-01/axel: adapted from burgers_phi
!
        use Cdata
        use General
        use Messages
        use Mpicomm
        use Sub
        use IO
        use Register
        use Global
        use Param_IO
        use Equ
        use Gravity
        use Slices
        use Print
        use Timestep
        use Snapshot
        use Boundcond
        use Filter
        use Power_spectrum
        use Timeavg
        use Interstellar,    only: check_SN
        use Shear
        use Testfield,       only: rescaling_testfield
        use TestPerturb,     only: testperturb_begin, testperturb_finalize
        use Forcing
        use EquationOfState
        use Dustvelocity,    only: init_uud
        use Dustdensity,     only: init_nd
        use NeutralVelocity, only: init_uun
        use NeutralDensity,  only: init_lnrhon
        use Magnetic,        only: rescaling_magnetic
        use Particles_main
        use Particles_nbody, only: particles_nbody_read_snapshot,&
                                   particles_nbody_write_snapshot
        use FArrayManager,   only: farray_clean_up
        use SharedVariables, only: sharedvars_clean_up
        use Entropy,         only: calc_heatcond_ADI
!
        implicit none
!
        real, dimension (mx,my,mz,mfarray) :: f
        real, dimension (mx,my,mz,mvar) :: df
        type (pencil_case) :: p
        double precision :: time1,time2
        integer :: count, ierr
        logical :: stop=.false.,timeover=.false.,resubmit=.false.
        logical :: lreinit_file=.false., exist=.false.
        logical :: lreload_file=.false., lreload_always_file=.false.
        real :: wall_clock_time=0., time_per_step=0.
        double precision :: time_last_diagnostic, time_this_diagnostic
        integer :: it_last_diagnostic,it_this_diagnostic
        integer :: i,ivar
        real, allocatable, dimension (:,:,:,:) :: finit
!
        lrun = .true.
!
        call initialize_messages()
!
!  initialize MPI and register physics modules
!  (must be done before lroot can be used, for example)
!
        call register_modules()
        if (lparticles) call particles_register_modules()
!
!  identify version
!
        if (lroot) call cvs_id( &
             "$Id$")
!
!  read parameters from start.x (default values; may be overwritten by
!  read_runpars)
!
        call rparam()
!
!  derived parameters (that may still be overwritten)
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
        xyz0_loc(1)=xyz0(1)
        xyz0_loc(2)=xyz0(2)+ipy*Lxyz_loc(2)
        xyz0_loc(3)=xyz0(3)+ipz*Lxyz_loc(3)
        xyz1_loc(1)=xyz1(1)
        xyz1_loc(2)=xyz0(2)+(ipy+1)*Lxyz_loc(2)
        xyz1_loc(3)=xyz0(3)+(ipz+1)*Lxyz_loc(3)
!
!  populate wavenumber arrays for fft and calculate nyquist wavenumber
!
        if (nxgrid/=1) then
          kx_fft=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*2*pi/Lx
          kx_ny = nxgrid/2 * 2*pi/Lx
        else
          kx_fft=0.0
          kx_ny = 0.0
        endif
!
        if (nygrid/=1) then
          ky_fft=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*2*pi/Ly
          ky_ny = nygrid/2 * 2*pi/Ly
        else
          ky_fft=0.0
          ky_ny = 0.0
        endif
!
        if (nzgrid/=1) then
          kz_fft=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*2*pi/Lz
          ky_ny = nzgrid/2 * 2*pi/Lz
        else
          kz_fft=0.0
          kz_ny = 0.0
        endif
!
!  read parameters and output parameter list
!
        call read_runpars()
        call rprint_list(LRESET=.false.)
!
!  inner radius for freezing variables defaults to r_min
!  Note: currently (July 2005), hydro.f90 uses a different approach:
!  r_int will override rdampint, which doesn't seem to make much sense (if
!  you want rdampint to be overridden, then don't specify it in the first
!  place)
        if (rfreeze_int == -impossible .and. r_int > epsi) &
             rfreeze_int = r_int
        if (rfreeze_ext == -impossible) rfreeze_ext = r_ext
!
!  Will we write all slots of f?
!
        if (lwrite_aux) then
          mvar_io = mvar+maux
        else
          mvar_io = mvar
        endif
!
!  print resolution and dimension of the simulation
!
        dimensionality=min(nxgrid-1,1)+min(nygrid-1,1)+min(nzgrid-1,1)
        if (lroot) write(*,'(a,i1,a)') 'This is a ',dimensionality,'-D run'
        if (lroot) print*, 'nxgrid,nygrid,nzgrid=',nxgrid,nygrid,nzgrid
!
!  check if we want to divide cdtv by dimensionality
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
!  set up directory names `directory' and `directory_snap'
!
        call directory_names()
!
!  get state length of random number generator
!
        call get_nseed(nseed)
!
!  read data
!  snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!  NOTE: for io_dist, rtime doesn't read the time, only for io_mpio.
!
!--     call rsnap(trim(directory_snap)//'/var.dat',f,mvar_io)
        call rsnap(trim(directory_snap)//'/var.dat',f,mvar)
        if (lparticles) &
           call particles_read_snapshot(trim(directory_snap)//'/pvar.dat')
        if (lparticles_nbody) &
             call particles_nbody_read_snapshot(&
             trim(datadir)//'/proc0/spvar.dat')

!  read time and global variables (if any)
!
        call rtime(trim(directory)//'/time.dat',t)
        if (lglobal) call rglobal()
        if (mglobal/=0)  &
            call input_globals(trim(directory_snap)//'/global.dat', &
            f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
!
!  read coordinates
!
        if (ip<=6.and.lroot) print*,'reading grid coordinates'
        call rgrid(trim(directory)//'/grid.dat')
!
!  read processor boundaries
!
        if (lparticles) then
          if (ip<=6.and.lroot) print*,'reading processor boundaries'
          call rproc_bounds(trim(directory)//'/proc_bounds.dat')
        endif
!
!  The following is here to avoid division in sub.f90 for diagnostic
!  outputs of integrated values in the non equidistant case
!
       if (.not. lequidist(1)) xprim=1./dx_1
       if (.not. lequidist(2)) yprim=1./dy_1
       if (.not. lequidist(3)) zprim=1./dz_1
!
!  Determine slice positions and whether slices are to be written on this
!  processor. This can only be done after the grid has been established.
!
        call setup_slices()
!
!  Write parameters to log file (done after reading var.dat, since we
!  want to output time t
!
        call print_runpars(FILE=trim(datadir)//'/params.log', &
                           ANNOTATION='Running')
!
!ajwm For backward compatibility... Read seed.dat if it exists
!ajwm this will be overridden however by any setting in the
!ajwm persistent data block of var.dat... Remove by August 05
!
        inquire(file=trim(directory)//'/seed.dat',exist=exist)
        if (exist.and.(lforcing .or. linterstellar)) then
           if (lroot.and.ip<14) print*, 'run: reading seed file'
           call inpui(trim(directory)//'/seed.dat',seed,nseed)
           call random_seed_wrapper(put=seed(1:nseed))
        endif
!
!  Allow modules to do any physics modules do parameter dependent
!  initialization. And final pre-timestepping setup.
!  (must be done before need_XXXX can be used, for example)
!
        call initialize_modules(f,LSTARTING=.false.)
!
!  initialize ionization array
!
        if (leos_ionization) call ioninit(f)
        if (leos_temperature_ionization) call ioncalc(f)
!
!  Prepare particles.
!
        if (lparticles) then
          call particles_rprint_list(.false.)
          call particles_initialize_modules(LSTARTING=.false.)
        endif
!
!  Allocate the finit array if lADI=.true.
!  Do this also when ltestperturb=.true.
!
        if (lADI) allocate(finit(mx,my,mz,mfarray))
!
!  Write data to file for IDL
!
        call wparam2()
!
!  possible debug output (can only be done after "directory" is set)
!  check whether mn array is correct
!
        if (ip<=3) call debug_imn_arrays
!
!  Find out which pencils are needed and write information about required,
!  requested and diagnostic pencils to disc
!
        call choose_pencils()
        call write_pencil_info()
!
        if (lglobal) call wglobal()
        if (mglobal/=0)  &
            call output_globals(trim(directory_snap)//'/global.dat', &
            f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
!
!  update ghost zones, so rprint works corrected for at the first
!  time step even if we didn't read ghost zones
!
        call update_ghosts(f)
!
!  save spectrum snapshot
!
        if (dspec/=impossible) call powersnap(f)
!
!  Initialize pencils in the pencil_case...
!
        if (lpencil_init) call initialize_pencils(p,0.0)
!
!  Perform pencil_case consistency check if requested
!
        if (lpencil_check) call pencil_consistency_check(f,df,p)
!
!  Start timing for final timing statistics
!  Initialize timestep diagnostics during the run (whether used or not,
!  see idiag_timeperstep)
!
        if (lroot) then
          time1 = mpiwtime()
          time_last_diagnostic = time1
          count = 0
          it_last_diagnostic = count
        endif
!        
        if (it1d==impossible_int) then 
          it1d=it1 
        else
          if (it1d < it1) then
            if (lroot) call stop_it("it1d smaller than it1")
          endif
        endif
!
!  Do loop in time
!
        Time_loop: do while (it<=nt)
          lout=mod(it-1,it1).eq.0
          l1dout=mod(it-1,it1d).eq.0
          if (lout) then
!
!  Exit do loop if file `STOP' exists
!
            if (lroot) inquire(FILE="STOP", EXIST=stop)
            call mpibcast_logical(stop, 1)
            if (stop.or.t>tmax) then
              if (lroot) then
                if (stop) print*, "done: found STOP file"
                if (t>tmax) print*, "done: t > tmax"
                call remove_file("STOP")
                inquire(FILE="RESUBMIT", EXIST=resubmit)
                if (resubmit) then 
                  print*, "Cannot be resubmitted"
                  call remove_file("RESUBMIT")
                else
                endif
              endif
              exit Time_loop
            endif
!
!  Exit do loop if wall_clock_time has exeeded max_walltime
!
            if (lroot.and.max_walltime > 0.0) then
              if (wall_clock_time > max_walltime) timeover=.true.
            endif
            call mpibcast_logical(timeover, 1)
            if (timeover) then
              if (lroot) print*, "done: max_walltime exceeded"
              exit Time_loop
            endif
!
!  Re-read parameters if file `RELOAD' exists; then remove the file
!  (this allows us to change parameters on the fly).
!  Re-read parameters if file `RELOAD_ALWAYS' exists; don't remove file
!  (only useful for debugging RELOAD issues).
!
            if (lroot) inquire(FILE="RELOAD", EXIST=lreload_file)
            if (lroot) inquire(FILE="RELOAD_ALWAYS", EXIST=lreload_always_file)
            lreloading = lreload_file .or. lreload_always_file
! In some compilers (particularly pathf90) the file reload is being give
! unit = 1 hence there is conflict during re-reading of parameters. 
! In this temporary fix, the RELOAD file is being removed just after it
! has been seen, not after RELOAD-ing has been completed. There must
! be a better solution. 
            if (lroot .and. lreload_file) call remove_file("RELOAD")
            call mpibcast_logical(lreloading, 1)
            if (lreloading) then
              if (lroot) write(0,*) 'Found RELOAD file -- reloading parameters'
!  Re-read configuration
              dt=0.
              call read_runpars(PRINT=.true.,FILE=.true.,ANNOTATION='Reloading')
              call rprint_list(LRESET=.true.) !(Re-read output list)
              call initialize_modules(f,LSTARTING=.false.)
              if (lparticles) then
                call particles_rprint_list(.false.)
                call particles_initialize_modules(LSTARTING=.false.)
              endif
              call choose_pencils()
              call wparam2()
              if (lroot .and. lreload_file) call remove_file("RELOAD")
              lreload_file        = .false.
              lreload_always_file = .false.
              lreloading          = .false.
            endif
!
!  Reinit variables found in `REINIT' file; then remove the file.
!
            if (lroot) inquire(FILE="REINIT", EXIST=lreinit_file)
            if (lroot .and. lreinit_file) then
              if (lroot) print*, 'Found REINIT file'
              open(1,file='REINIT',action='read',form='formatted')
              nreinit=1
!
!  Read variable names from REINIT file
!
              ierr=0
              do while (ierr == 0)
                read(1,'(A5)',IOSTAT=ierr) reinit_vars(nreinit)
                if (reinit_vars(nreinit) /= '') nreinit=nreinit+1
              enddo
              close(1)
              nreinit=nreinit-1
              lreinit=.true.
              if (lroot) call remove_file("REINIT")
            endif
            call mpibcast_logical(lreinit, 1)
            if (lreinit) then
              if (lroot) print*, 'Reiniting variables: ', reinit_vars(1:nreinit)
              call mpibcast_int(nreinit, 1)
              call mpibcast_char(reinit_vars, 10)
!
!  Reinit all variables present in reinit_vars array
!
              do ivar=1,nreinit
                select case (reinit_vars(ivar))
                  case ('uud')
                    f(:,:,:,iudx(1):iudz(ndustspec))=0.
                    call init_uud(f)
                  case ('nd')
                    f(:,:,:,ind)=0.
                    call init_nd(f)
                  case ('particles')
                    call particles_init(f)
                  case default
                    if (lroot) print*, &
                        'reinit: Skipping unknown variable ', &
                        reinit_vars(ivar)
                endselect
              enddo
              call choose_pencils()
              lreinit=.false.
              lreinit_file=.false.
              reinit_vars=''
            endif
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
            if (mod(it,iwig).eq.0) then
              if (lrmwig_xyaverage) call rmwig_xyaverage(f,ilnrho)
              if (lrmwig_full) call rmwig(f,df,ilnrho,ilnrho,awig)
              if (lrmwig_rho) call rmwig(f,df,ilnrho,ilnrho,awig,explog=.true.)
            endif
          endif
!
!  If we want to write out video data, wvid sets lvid=.true.
!  This allows pde to prepare some of the data
!
          if (lwrite_slices) call wvid_prepare()
          if (lwrite_2daverages) call write_2daverages_prepare()
!
!  Find out which pencils to calculate at current time-step
!
          lpencil = lpenc_requested
          if (lout)   lpencil=lpencil .or. lpenc_diagnos
          if (l2davg) lpencil=lpencil .or. lpenc_diagnos2d
          if (lvid)   lpencil=lpencil .or. lpenc_video
!
!  save state vector prior to update
!
          if (lADI)   finit=f
          if (ltestperturb) call testperturb_begin(f,df)
!
!  Time advance
!
          call rk_2n(f,df,p)
!
! 07-Sep-07/dintrans+gastine: implicit advance of the radiative diffusion
! in the temperature equation (using temperature_idealgas)
!
          if (lADI) call calc_heatcond_ADI(finit,f)
          if (ltestperturb) call testperturb_finalize(f)
!
          if (lroot) then
            count = count + 1     !  reliable loop count even for premature exit
          endif
!
!  Update time averages and time antegrals
!
          call diagnostic
          if (ltavg) call update_timeavgs(f,dt)
!
!  Add forcing and/or do rescaling (if applicable)
!
          if (lforcing) call addforce(f)
          if (lrescaling_magnetic)  call rescaling_magnetic(f)
          if (lrescaling_testfield) call rescaling_testfield(f)
!
!  Check for SNe, and update f if necessary (see interstellar.f90)
!
          if (linterstellar) call check_SN(f,df)
!
          if (lout.and.lroot.and.(idiag_walltime/=0 .or. max_walltime/=0.)) then
            time2=mpiwtime()
            wall_clock_time = (time2-time1)
            if (idiag_walltime/=0) &
                call save_name(wall_clock_time,idiag_walltime)
          endif
!
          if (lout.and.lroot.and.idiag_timeperstep/=0) then
              it_this_diagnostic = it
            time_this_diagnostic = mpiwtime()
            time_per_step = (time_this_diagnostic - time_last_diagnostic) &
                           /(  it_this_diagnostic -   it_last_diagnostic)
              it_last_diagnostic =   it_this_diagnostic
            time_last_diagnostic = time_this_diagnostic
            call save_name(time_per_step,idiag_timeperstep)
          endif
!
          if (lout) call prints()
!
!  In regular intervals, calculate certain averages and do other output
!
          call write_1daverages()
          call write_2daverages()
!
!  Setting ialive=1 can be useful on flaky machines!
!  Each processor writes it's processor number (if it is alive!)
!  Set ialive=0 to fully switch this off
!
          if (ialive /= 0) then
            if (mod(it,ialive)==0) &
                 call outpui(trim(directory)//'/alive.info', &
                 spread(it,1,1) ,1) !(all procs alive?)
          endif
          call wsnap(trim(directory_snap)//'/VAR',f, &
              mvar_io,ENUM=.true.,FLIST='varN.list')
          if (lparticles) &
              call particles_write_snapshot(trim(directory_snap)//'/PVAR', &
              ENUM=.true.,FLIST='pvarN.list')
!  This is weird... if I write only to the root, the other processors complain...
          if (lparticles_nbody) &
              call particles_nbody_write_snapshot(&
              trim(directory_snap)//'/SPVAR',&
              ENUM=.true.,FLIST='spvarN.list')
          call wsnap_timeavgs(trim(directory_snap)//'/TAVG',ENUM=.true., &
               FLIST='tavgN.list')
!
!  Write slices (for animation purposes)
!
          if (lvid.and.lwrite_slices) call wvid(f,trim(directory)//'/slice_')
!
!  Save snapshot every isnap steps in case the run gets interrupted.
!  The time needs also to be written
!
          if (isave/=0.and..not.lnowrite) then
            if (mod(it,isave)==0) then
              call wsnap(trim(directory_snap)//'/var.dat', &
                         f,mvar_io,ENUM=.false.)
              if (lparticles) &
                  call particles_write_snapshot( &
                  trim(directory_snap)//'/pvar.dat',ENUM=.false.)
              if (lparticles_nbody.and.lroot) &
                  call particles_nbody_write_snapshot( &
                  trim(datadir)//'/proc0/spvar.dat',ENUM=.false.)
              call wsnap_timeavgs(trim(directory_snap)//'/timeavg.dat', &
                                  ENUM=.false.)
              call wtime(trim(directory)//'/time.dat',t)
            endif
          endif
!
!  Save spectrum snapshot
!
          if (dspec/=impossible) call powersnap(f)
!
!  Save global variables.
!
          if (isaveglobal/=0) then
            if (mod(it,isaveglobal)==0) then
              if (lglobal) call wglobal()
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
          if ((it < nt) .and. (dt < dtmin)) then
            if (lroot) &
                write(0,*) 'run: Time step has become too short: dt = ', dt
            call timestep_autopsy()
            save_lastsnap=.false.
            exit Time_loop
          endif
!
!  Fatal errors sometimes occur only on a specific processor. In that case all
!  processors must be informed about the problem before the code can stop.
!
          call fatal_error_local_collect()
!
          it=it+1
          headt=.false.
        enddo Time_loop
        if (lroot) time2=mpiwtime()
!
!  Write data at end of run for restart.
!  dvar is written for analysis purposes only
!
        if (lroot) print*, 'Writing final snapshot for t=', t
        call wtime(trim(directory)//'/time.dat',t)
        if (save_lastsnap.and..not.lnowrite) then
          call wsnap(trim(directory_snap)//'/var.dat',f,mvar_io,ENUM=.false.)
          if (lparticles) call particles_write_snapshot( &
              trim(directory_snap)//'/pvar.dat',ENUM=.false.)
          if (lparticles_nbody.and.lroot) &
               call particles_nbody_write_snapshot( &
               trim(datadir)//'/proc0/spvar.dat',ENUM=.false.)
          call wsnap_timeavgs(trim(directory_snap)//'/timeavg.dat', &
                                  ENUM=.false.)
          if (ip<=11 .or. lwrite_dvar) then
            call wsnap(trim(directory)//'/dvar.dat',df,mvar, &
                       enum=.false.,noghost=.true.)
            call particles_write_dsnapshot(trim(directory)//'/dpvar.dat')
          endif
        else
          call wsnap(trim(directory_snap)//'/crash.dat',f,mvar_io,ENUM=.false.)
          if (ip<=11) &
               call wsnap(trim(directory)//'/dcrash.dat',df,mvar,ENUM=.false.)
        endif
!
!  Save spectrum snapshot
!
        if (save_lastsnap) then
          if (dspec /= impossible) call powersnap(f,.true.)
        endif
!
!  Print wall clock time and time per step and processor  for diagnostic purposes
!
        if (lroot) then
          wall_clock_time = time2-time1
          print*
          write(*,'(A,1pG10.3,A,1pG8.2,A)') &
               ' Wall clock time [hours] = ', wall_clock_time/3600., &
               ' (+/- ', real(mpiwtick())/3600.,')'
          if (it>1) write(*,'(A,1pG10.3)') &
               ' Wall clock time/timestep/meshpoint [microsec] =', &
               wall_clock_time/count/nw/ncpus/1e-6
          print*
        endif
        call mpifinalize
!
!
!  Free any allocated memory
!
        call farray_clean_up()
        call sharedvars_clean_up()
        if (lADI) deallocate(finit)
!
      endprogram run

