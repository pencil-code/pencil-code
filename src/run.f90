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
        integer :: time1,time2,count_rate
        real, dimension (mx,my,mz,mvar) :: f,df
!     
!  initialize MPI
!
        call siginit
        call signonbrutal
!
        call initialize         ! register modules, etc.
!
!  ix,iy,iz are indices for checking variables at some selected point
!  set default values
!
        ix=mx/2; iy=my/2; iz=mz/2
!
!  read in parameters
!  nt is the number of timesteps
!  dsnap is the time interval of snapshot output
!  it1 is the frequency of output of rms and max values
!  nu is viscosity
!
        call cread
        call cprint
!
!  timestep
!
        ldt=dt.lt.0.
        if (ldt) cdt=abs(dt)
!
!  read data
!  snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!
        if (ip<=6) print*,'read var files'
        call input(trim(directory)//'/var.dat',f,mvar,1)
        call read_global()
!
!  read seed field parameters (only if forcing is turned on)
!
        if (iforce/=0) then
          if (lroot) print*,'read seed file'
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
          print*,'start time loop'
          print*,'$Id: run.f90,v 1.6 2002-01-10 18:04:21 dobler Exp $'
        endif
!
!  do loop in time
!
        Time_loop: do it=1,nt
          if (ip.le.10) print*,'it=',it
          lout=mod(it-1,it1).eq.0
          if (lout) call cread
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
        enddo Time_loop
        if(lroot) call system_clock(count=time2)
!
!  write data at end of run for restart
!  dvar is written for analysis purposes only
!
!  update ghost zones for var.dat (cheap, since done once
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
        if(lroot) &
             print*,'Wall clock time=',(time2-time1)/real(count_rate), &
                    ' (+/- ', 1./count_rate,')'
        call mpifinalize
!
      endprogram run
