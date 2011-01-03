! $Id$
!
!***********************************************************************
!
!  This helper program extracts a velocity field with the given cadence
!  and writes it in a 'driver/vel_field.dat' file suitable for reading
!  with the 'special/solar_corona.f90' module.
!
!***********************************************************************
program write_driver
!
  use Boundcond,       only: update_ghosts
  use Cdata
  use Diagnostics
  use Equ,             only: debug_imn_arrays,initialize_pencils
  use FArrayManager,   only: farray_clean_up
  use IO
  use Messages
  use Mpicomm
  use Param_IO
  use Register
  use SharedVariables, only: sharedvars_clean_up
  use Signal_handling, only: signal_prepare, emergency_stop
  use Snapshot, only: rsnap
  use Sub
  use Special
!
  implicit none
!
!
  ! SET THE CADENCE FOR 'vel_field.dat' HERE:
  real, parameter :: cadence = 15.0 ! [s]
!
!
  real, dimension (mx,my,mz,mfarray) :: f
  type (pencil_case) :: p
  double precision :: time1, time2
  double precision :: time_last_diagnostic, time_this_diagnostic
  real :: wall_clock_time=0.0, time_per_step=0.0, next_snapshot=0
  integer :: icount, i, mvar_in, rn, isnap=0
  integer :: it_last_diagnostic, it_this_diagnostic
  logical :: lstop=.false., timeover=.false., resubmit=.false.
!
  inquire (IOLENGTH=rn) 1.
!
  lrun=.true.
!
!  Get processor numbers and define whether we are root.
!
  call mpicomm_init
!
!  Identify version.
!
  if (lroot) call svn_id ('$Id$')
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
!
!  Define the lenergy logical
!
  lenergy=lentropy.or.ltemperature.or.lthermal_energy
!
!  Inform about verbose level.
!
  if (lroot) print*, 'The verbose level is ip=', ip, ' (ldebug=', ldebug, ')'
!
!  Call rprint_list to initialize diagnostics and write indices to file.
!
  call rprint_list(LRESET=.false.)
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
!  Read data.
!  Snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!  NOTE: for io_dist, rtime doesn't read the time, only for io_mpio.
!
  call rsnap(trim(directory_snap)//'/var.dat',f,mvar_in)
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
!  The following is here to avoid division in sub.f90 for diagnostic
!  outputs of integrated values in the non equidistant case.
!  Do this even for uniform meshes, in which case xprim=dx, etc.
!  Remember that dx_1=0 for runs without extent in that direction.
!
  if (nxgrid==1) then; xprim=1.0; else; xprim=1/dx_1; endif
  if (nygrid==1) then; yprim=1.0; else; yprim=1/dy_1; endif
  if (nzgrid==1) then; zprim=1.0; else; zprim=1/dz_1; endif
!
!  Allow modules to do any physics modules do parameter dependent
!  initialization. And final pre-timestepping setup.
!  (must be done before need_XXXX can be used, for example)
!
  call initialize_modules(f,LSTARTING=.false.)
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
    endif
!
!  Find out which pencils to calculate at current time-step.
!
    lpencil = lpenc_requested
!
!  Compute velocity filed.
!
    call special_before_boundary(f)
!
!  Save velocity filed.
!
    if (t >= next_snapshot) then
!
      write (*,*) 'snapshot: ', isnap, t
      write (*,*) 'Ux: ', minval (f(l1:l2,m1:m2,n1,iux)), maxval (f(l1:l2,m1:m2,n1,iux))
      write (*,*) 'Uy: ', minval (f(l1:l2,m1:m2,n1,iuy)), maxval (f(l1:l2,m1:m2,n1,iuy))
!
      open (17, file='driver/vel_field.dat', status='unknown', form='unformatted', recl=nxgrid*nygrid*rn, access='direct')
      write (17, rec=isnap*2+1) f(l1:l2,m1:m2,n1,iux)
      write (17, rec=isnap*2+2) f(l1:l2,m1:m2,n1,iuy)
      close (17)
!
      open (2, file='driver/vel_times.dat', status='unknown', form='unformatted', recl=rn, access='direct')
      write (2, rec=isnap+1) t
      close (2)
!
      isnap = isnap + 1
      next_snapshot = isnap * cadence * unit_time
    endif
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
!  Exit do loop if wall_clock_time has exceeded max_walltime.
!
    if (max_walltime>0.0) then
      timeover=(wall_clock_time>max_walltime)
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
!  Print diagnostic averages to screen and file.
!
    tdiagnos = t
    if (lout)   call prints()
!
!  Time advance.
!
    t = t + unit_time
!
    it = it + 1
    headt = .false.
!
    if (lroot) icount=icount+1  !  reliable loop count even for premature exit
!
  enddo Time_loop
!
  if (lroot) then
    print*
    print*, 'Simulation finished after ', icount, ' time-steps'
  endif
!
  if (lroot) time2=mpiwtime()
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
  call vnames_clean_up()
!
endprogram write_driver
