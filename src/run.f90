! $Id: run.f90,v 1.20 2002-05-02 20:02:27 brandenb Exp $
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
        use Register
        use Global
        use Forcing
        use Equ
        use Slices
        use Timestep
!
        implicit none
!
        real, dimension (mx,my,mz,mvar) :: f,df
        integer :: time1,time2,count_rate
!        real :: time1,time2     ! cpu_time can measure longer times than
                                ! system clock, but takes 15 seconds to
                                ! calibrate at startup with Intel F95
        logical :: stop=.false.,reload=.false.
!     
!  initialize MPI
!
!        call siginit
!        call signonbrutal
!
!  identify version
!
        if (lroot) call cvs_id( &
             "$RCSfile: run.f90,v $", &
             "$Revision: 1.20 $", &
             "$Date: 2002-05-02 20:02:27 $")
!
        call initialize         ! register modules, etc.
!
!  ix,iy,iz are indices for checking variables at some selected point
!  set default values
!
        ix=mx/2; iy=my/2; iz=mz/2
        dtmin=1.e-6
!
!  read in parameters
!  nt is the number of timesteps
!  dsnap is the time interval of snapshot output
!  it1 is the frequency of output of rms and max values
!  nu is viscosity
!
        call rparam             ! Read parameters from start.x;
                                ! these may be overwritten by cread
        call cread(PRINT=.true.)
!
!  read data
!  snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!
        if (ip<=6) print*,'reading var files'
        call input(trim(directory)//'/var.dat',f,mvar,1)
        call rglobal()      ! Read global variables (if there are)
!
!  read coordinates
!
        if (ip<=6) print*,'reading grid coordinates'
        call rgrid(trim(directory)//'/grid.dat')
!
!  read seed field parameters (only if forcing is turned on)
!
        if (iforce/=0) then
          if (lroot) print*,'reading seed file'
          call inpui(trim(directory)//'/seed.dat',seed,2)
          if (iproc < 10) print*,'iproc,seed=',iproc,seed
          call random_seed(put=seed)
        endif
!
!  warn about the damping term
!
        if (lroot .and. (dampu /= 0.) .and. (t < tdamp)) then 
          print*, 'Damping velocities until time ', tdamp
        endif
!
!  advance equations
!  NOTE: headt=.true. in order to print header titles
!
        headt=.true.
        if(lroot) then
          call system_clock(count_rate=count_rate)
          call system_clock(count=time1)
!          call cpu_time(time1)
          print*,'start time loop'
        endif
!
!  do loop in time
!
        Time_loop: do it=1,nt
          if (ip.le.10) print*,'it=',it
          lout=mod(it-1,it1).eq.0
          if (lout) then
            inquire(FILE="STOP", EXIST=stop)! Exit DO loop if the file
                                            ! `STOP' exists
            if (stop) then
              if (lroot) write(0,*) "Found STOP file -- quitting"
              exit Time_loop
            endif
!
!  Re-read parameters if a file `RELOAD' has appeared
!
            inquire(FILE="RELOAD", EXIST=reload)
            if (reload) then
              if (lroot) write(0,*) "Found RELOAD file -- reloading parameters"
              call cread(PRINT=.true.) ! Re-read configuration
              call remove_file("RELOAD")
              reload = .false.
            endif
          endif

          if (iforce==1) call forcing1
          if (iforce==2) call forcing2
!
!  time advance
!
          call rk_2n(f,df)
          if(lout) call prints
          call wsnap(trim(directory)//'/VAR',f)
          call wvid(trim(directory))
!
!  save snapshot every isnap steps in case the run gets interrupted
!
          if (mod(it,isave).eq.0) then
!  update ghost zones for var.dat (cheap, since done infrequently)
            call initiate_isendrcv_bdry(f)
            call finalise_isendrcv_bdry(f)
!  write data
            call output(trim(directory)//'/var.dat',f,mvar)
          endif
!
          headt=.false.
          if (it>=nt) exit Time_loop
          if (dt < dtmin) then
            write(0,*) 'Time step has become too short: dt = ', dt
            exit Time_loop
          endif
        enddo Time_loop
        if(lroot) call system_clock(count=time2)
!        if(lroot) call cpu_time(time2)
!
!  write data at end of run for restart
!  dvar is written for analysis purposes only
!
!  update ghost zones for var.dat (cheap, since done once)
        call initiate_isendrcv_bdry(f)
        call finalise_isendrcv_bdry(f)
        call output(trim(directory)//'/var.dat',f,mvar)
        if (ip<=10) call output(trim(directory)//'/dvar.dat',df,mvar)
!
!  write seed parameters (only if forcing is turned on)
!
        if (iforce/=0) then
          call random_seed(get=seed)
          call outpui(trim(directory)//'/seed.dat',seed,2)
        endif
!
!  print wall clock time and time per step and processor
!  for diagnostic purposes
!
        if(lroot) &
             print*,'Wall clock time=',(time2-time1)/real(count_rate), &
                    ' (+/- ', 1./count_rate,')'
             print*,'time/step/pt [microsec]=',(time2-time1)/real(count_rate)/(it-1)/mw/1e-6
        call mpifinalize
!
      endprogram run
