! $Id: run.f90,v 1.29 2002-05-27 12:04:32 dobler Exp $
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
        use Print
        use Timestep
!
        implicit none
        real, dimension (mx,my,mz,mvar) :: f,df
        integer :: time1,time2,count_rate
        logical :: stop=.false.,reload=.false.
        real :: Wall_clock_time
!
!  initialize MPI
!
!       call siginit
!       call signonbrutal
!
!  initialize MPI and register physics modules
!  (must be done before lroot can be used, for example)
!
        call initialize
!
!  identify version
!
        if (lroot) call cvs_id( &
             "$RCSfile: run.f90,v $", &
             "$Revision: 1.29 $", &
             "$Date: 2002-05-27 12:04:32 $")
!
!  ix,iy,iz are indices for checking variables at some selected point
!  set default values
!
        ix=mx/2; iy=my/2; iz=mz/2
        dtmin=1e-6  !!(AB: this should be made an input parameter, better dimless)
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
!  read the print parameter list
!
        call rprint_list
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
          if (iproc < 10) print*,'iproc,seed(1:2)=',iproc,seed
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
!
!  time advance
!
          call rk_2n(f,df)
          if (iforce/=0) call addforce(f)
          if(lout) call prints
          call wsnap(trim(directory)//'/VAR',f,.true.)
          call wvid(trim(directory))
!
!  save snapshot every isnap steps in case the run gets interrupted
!
          if (mod(it,isave).eq.0) call wsnap(trim(directory)//'/var.dat',f,.false.)
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
        call wsnap(trim(directory)//'/var.dat',f,.false.)
        if (ip<=10) call wsnap(trim(directory)//'/dvar.dat',df,.false.)
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
        if(lroot) then
          Wall_clock_time=(time2-time1)/real(count_rate)
          print*,'Wall clock time=',Wall_clock_time,' (+/- ', 1./count_rate,')'
          if (it>1) print*, 'time/step/pt [microsec]=', &
               Wall_clock_time/(it-1)/nw/ncpus/1e-6
        endif
        call mpifinalize
!
      endprogram run
