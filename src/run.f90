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
        use Mpicomm
        use Cdata
        use Forcing
        use Slices
        use Sub
        use Timestep
        use Equ
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
        call mpicomm_init
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
!
!  read seed field parameters (only if forcing is turned on)
!
        if (iforce/=0) then
          if (iproc==root) print*,'read seed file'
          call inpui(trim(directory)//'/seed.dat',seed,2)
          if (iproc < 10) print*,'iproc,seed=',iproc,seed
          call random_seed(put=seed)
        endif
!
!  advance equations
!  NOTE: headt=.true. in order to print header titles
!
        headt=.true.
        if(iproc==root) then
          call system_clock(count_rate=count_rate)
          call system_clock(count=time1)
          print*,'start time loop'
          print*,'$Id: run.f90,v 1.1.1.1 2001-11-01 15:07:56 dobler Exp $'
        endif
!
!  do loop in time
!
        Time_loop: do it=1,nt
          if (ip.le.10) print*,'it=',it
          lout=mod(it,it1).eq.0
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
          if (mod(it,isave).eq.0) &
               call output(trim(directory)//'/var.dat',f,mvar)
!
          headt=.false.
          if (it>=nt) exit Time_loop
        enddo Time_loop
        if(iproc==root) call system_clock(count=time2)
!
!  write data at end of run for restart
!  dvar is written for analysis purposes only
!
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
        if(iproc==root) &
             print*,'Wall clock time=',(time2-time1)/real(count_rate),count_rate
        call mpifinalize
!
      endprogram run
