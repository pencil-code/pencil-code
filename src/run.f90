! $Id: run.f90,v 1.91 2002-10-02 10:37:04 dobler Exp $
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
        use Wsnaps
        use Boundcond
        use Power_spectrum
!
        implicit none
!
        real, dimension (mx,my,mz,mvar) :: f,df
        integer :: time1,time2,count,count_rate
        logical :: stop=.false.,reload=.false.
        real :: Wall_clock_time
!
!  initialize MPI and register physics modules
!  (must be done before lroot can be used, for example)
!
        call initialize
!
!  call signal handler (for compaq machine only)
!  currently disabled; want to put as include file
!
!!     call siginit
!!     call signonbrutal
!
!  identify version
!
        if (lroot) call cvs_id( &
             "$Id: run.f90,v 1.91 2002-10-02 10:37:04 dobler Exp $")
!
!  ix,iy,iz are indices for checking variables at some selected point
!  set default values (should work also for 1-D and 2-D runs)
!
        ix=1+(mx-1)/2; iy=1+(my-1)/2; iz=1+(mz-1)/2
!        dtmin=1e-6  !!(AB: this should be an input parameter, better dimless)
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
!  [Currently none]
!
!  read parameters and output parameter list
!
        call read_runpars(PRINT=.true.)
        call rprint_list(.false.)
!
!  print resolution
!
        if (lroot) print*, 'nxgrid,nygrid,nzgrid=',nxgrid,nygrid,nzgrid
!
!  read data
!  snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!
        if (ip<=6.and.lroot) print*,'reading var files'
        call input(trim(directory_snap)//'/var.dat',f,mvar,1)
        call rtime(trim(directory)//'/time.dat',t)
        call rglobal()      ! Read global variables (if there are)
!
!  read coordinates
!
        if (ip<=6.and.lroot) print*,'reading grid coordinates'
        call rgrid(trim(directory)//'/grid.dat')
!
        call get_nseed(nseed)   ! get state length of random number generator
!
!  run initialization of individual modules
!
        call ss_run_hook()      ! calculate radiative conductivity, etc.
        call forcing_run_hook() ! get random seed from file, ..
!
!  Write data to file for IDL
!
      call wparam2()
!
!  setup gravity (obtain coefficients cpot(1:5); initialize global array gg)
!
        if (lgravr) call setup_grav()
!
        call wglobal()

!
!  advance equations
!  NOTE: headt=.true. in order to print header titles
!
        if(lroot) then
          call system_clock(count_rate=count_rate)
          call system_clock(count=time1)
          count = 0
        endif
!
!  update ghost zones, so rprint works corrected for at the first
!  time step even if we didn't read ghost zones
!
        call update_ghosts(f)
!
!  do loop in time
!
        Time_loop: do it=1,nt
          lout=mod(it-1,it1).eq.0
          if (lout) then
            !
            ! exit DO loop if file `STOP' exists
            !
            if (lroot) inquire(FILE="STOP", EXIST=stop)
            call mpibcast_logical(stop, 1)
            if (stop.or.t>tmax) then
              if (lroot) then
                if (stop) write(0,*) "done: found STOP file"
                if (t>tmax) write(0,*) "done: t > tmax"
                call remove_file("STOP")
              endif
              exit Time_loop
            endif
            !
            !  re-read parameters if file `RELOAD' exists; then remove
            !  the file
            !
            if (lroot) inquire(FILE="RELOAD", EXIST=reload)
            call mpibcast_logical(reload, 1)
            if (reload) then
              if (lroot) write(0,*) "Found RELOAD file -- reloading parameters"
              call read_runpars(PRINT=.true.) !(Re-read configuration)
              call rprint_list(.true.) !(Re-read output list)
              call mpibarrier() ! all CPUs must read RELOAD before it is removed
              if (lroot) call remove_file("RELOAD")
              reload = .false.
            endif
          endif
          !
          !  remove wiggles in lnrho in sporadic time intervals
          !  Necessary if the Reynolds number is large.
          !  iwig=500 is a typical value. For iwig=0 no action is taken.
          !  (The two queries below must come separately on compaq machines.)
          !
          if (iwig/=0) then
            if (mod(it,iwig).eq.0) then
              if (lrmwig_xyaverage) call rmwig_xyaverage(f,ilnrho)
              if (lrmwig_full) call rmwig(f,df,ilnrho,awig)
              !call rmwig(f,df,ilnrho,explog=.true.)
            endif
          endif
          !
          !  time advance
          !
          call rk_2n(f,df)
          count = count + 1     !  reliable loop count even for premature exit
          !
          !  advance shear parameter and add forcing (if applicable)
          !
          if (lshear) call advance_shear
          if (lforcing) call addforce(f)
          !
          !  in regular intervals, calculate certain averages
          !  and do other output.
          !
          if(lout) call write_xyaverages
          if(lout) call write_zaverages
          if(lout) call prints
          if (ialive /= 0) then ! set ialive=0 to fully switch this off
            if (mod(it,ialive)==0) &
                 call outpui(trim(directory)//'/alive.info', &
                 spread(it,1,1) ,1) !(all procs alive?)
          endif
          call wsnap(trim(directory_snap)//'/VAR',f,.true.)
          call wvid(trim(directory))
          !
          !  save snapshot every isnap steps in case the run gets interrupted
          !
          if (isave /= 0) then
            if (mod(it,isave)==0) &
                 call wsnap(trim(directory_snap)//'/var.dat',f,.false.)
          endif
          !
          !  save spectrum snapshot
          !
          if ((t>spect) .AND. (dspect .NE. impossible)) then
             spect=spect+dspect 
             if (vel_spec) call power(f,'u')
             if (mag_spec) call power(f,'b')
             if (vec_spec) call power(f,'a')
          endif
          !
          headt=.false.
          if ((it < nt) .and. (dt < dtmin)) then
            write(0,*) 'Time step has become too short: dt = ', dt
            exit Time_loop
          endif

        enddo Time_loop
        if(lroot) call system_clock(count=time2)
!
!  write data at end of run for restart
!  dvar is written for analysis purposes only
!
        call wsnap(trim(directory_snap)//'/var.dat',f,.false.)
        call wtime(trim(directory)//'/time.dat',t)
        if (ip<=10) call wsnap(trim(directory)//'/dvar.dat',df,.false.)
!
!  save spectrum snapshot
!
        if (vel_spec) call power(f,'u')
        if (mag_spec) call power(f,'b')
        if (vec_spec) call power(f,'a')
!
!  write seed parameters (only if forcing is turned on)
!
        if (lforcing) then
          call random_seed_wrapper(get=seed(1:nseed))
          call outpui(trim(directory)//'/seed.dat',seed,nseed)
        endif
!
!  print wall clock time and time per step and processor
!  for diagnostic purposes
!
        if(lroot) then
          Wall_clock_time=(time2-time1)/real(count_rate)
          print*
          print*,'Wall clock time [sec]=',Wall_clock_time,' (+/- ', 1./count_rate,')'
          if (it>1) print*, 'Wall clock time/timestep/meshpoint [microsec]=', &
               Wall_clock_time/count/nw/ncpus/1e-6
          print*
        endif
        call mpifinalize
!
      endprogram run
