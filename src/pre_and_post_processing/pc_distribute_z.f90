! This tool distributes a global data cube into the proc-directories.
!
! $Id$
!***********************************************************************
program pc_distribute_z
!
  use Cdata
  use Cparam, only: fnlen
  use IO
  use Messages
  use Param_IO
  use Register
  use Sub
  use General, only: delete_file
!
  implicit none
!
  character (len=fnlen) :: filename
  character (len=*), parameter :: directory_in = 'data/allprocs'
!
  real, dimension (:,:,:,:), allocatable :: f
  real, dimension (mxgrid) :: gx
  real, dimension (mygrid) :: gy
  real, dimension (mzgrid) :: gz
  logical :: ex
  integer :: mvar_in, pz, pa, io_len, alloc_err
  integer(kind=8) :: rec_len
  real :: t_sp   ! t in single precision for backwards compatibility
!
  lstart = .true.
  lmpicomm = .false.
  ipx = 0
  ipy = 0
  ipz = 0
  ylneigh = 0
  zlneigh = 0
  yuneigh = 0
  zuneigh = 0
!
  deltay = 0.0   ! Shearing not available due to missing fseek in Fortran
!
  if (IO_strategy /= "collect_xy") call fatal_error ('pc_distribute_z', &
      "This tool only makes sense together with the 'io_collect_xy' module.")
!
  write (*,*) 'Please enter the filename to convert (eg. var.dat, VAR1, ...):'
  read (*,*) filename
!
!  Identify version.
!
  call svn_id( &
      '$Id$')
!
!  Initialize the message subsystem, eg. color setting etc.
!
  call initialize_messages
!
!  Read parameters from start.x (default values; overwritten by 'read_all_run_pars').
!
  call read_all_init_pars
!
!  Read parameters and output parameter list.
!
  call read_all_run_pars
!
!  Derived parameters (that may still be overwritten).
!  [might better be put into another routine, possibly in 'read_all_run_pars']
!
  x0 = xyz0(1)
  y0 = xyz0(2)
  z0 = xyz0(3)
  Lx = Lxyz(1)
  Ly = Lxyz(2)
  Lz = Lxyz(3)
!
! Calculate dimensionality
!
  dimensionality = min(nxgrid-1,1) + min(nygrid-1,1) + min(nzgrid-1,1)
!
!  Register physics modules.
!
  call register_modules
!
!  Define the lenergy logical
!
  lenergy = lentropy .or. ltemperature .or. lthermal_energy
!
  if (lwrite_aux .and. .not. lread_aux) then
    print *, ''
    print *, 'lwrite_aux=T but lread_aux=F'
    print *, 'The code will write the auxiliary variables to allprocs/VARN'
    print *, ' without having read them from proc*/VARN'
    print *, ''
    call fatal_error("pc_distribute_z","Stop and check")
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
  allocate (f (mxgrid,mygrid,mz,mfarray), stat=alloc_err)
  if (alloc_err /= 0) call fatal_error ('pc_distribute_z', 'Failed to allocate memory for f.', .true.)
!
!  Print resolution and dimension of the simulation.
!
  write (*,'(a,i1,a)') ' This is a ', dimensionality, '-D run'
  print *, 'nxgrid, nygrid, nzgrid=', nxgrid, nygrid, nzgrid
  print *, 'Lx, Ly, Lz=', Lxyz
  print *, '      Vbox=', Lxyz(1)*Lxyz(2)*Lxyz(3)
!
  inquire (file=trim(directory_in)//'/'//filename, exist=ex)
  if (.not. ex) call fatal_error ('pc_distribute_z', 'File not found: '//trim(directory_in)//'/'//filename, .true.)
  inquire (file=trim(directory_in)//'/grid.dat', exist=ex)
  if (.not. ex) call fatal_error ('pc_distribute_z', 'File not found: '//trim(directory_in)//'/grid.dat', .true.)
!
  ! read time:
  open (lun_input, FILE=trim(directory_in)//'/'//filename, FORM='unformatted', status='old', position='append')
  backspace(lun_input)
  read (lun_input) t_sp, gx, gy, gz, dx, dy, dz
  close (lun_input)
  t = t_sp
!
  inquire (IOLENGTH=io_len) t_sp
!
!  Allow modules to do any physics modules do parameter dependent
!  initialization. And final pre-timestepping setup.
!  (must be done before need_XXXX can be used, for example)
!
  call initialize_modules(f)

  call read_and_distribute(filename,f(:,:,:,1:mvar_io),.false.)
  if (mglobal>0) &
    call read_and_distribute('global.dat',f(:,:,:,mvar+maux+1:mvar+maux+mglobal),.true.)
!
!  Give all modules the possibility to exit properly.
!
  call finalize_modules (f)
!
!  Free any allocated memory.
!
  deallocate (f)
!
contains
!*************************************************************************************************
subroutine read_and_distribute(filename,f,lonly_farray)

  character(LEN=*) :: filename
  real, dimension(:,:,:,:) :: f
  logical :: lonly_farray

  integer :: mvar

  mvar = size(f,4)
  open (lun_input, FILE=trim(directory_in)//'/'//filename, access='direct', recl=mxgrid*mygrid*io_len, status='old')
!
! Loop over processors
!
  write (*,*) 'writing "'//trim(filename)//'" IPZ-layer:'
!
  do ipz = 0, nprocz-1
!
    write (*,*) ipz+1, " of ", nprocz
!
    f = huge(1.0)
!
    ! read xy-layer:
    do pa = 1, mvar
      do pz = 1, mz
        read (lun_input, rec=pz+ipz*nz+(pa-1)*mzgrid) f(:,:,pz,pa)
      enddo
    enddo
!
    iproc = ipx + ipy * nprocx + ipz * nprocx*nprocy
    lroot = (iproc==root)
!
!  Set up flags for leading processors in each possible direction and plane
!
    lfirst_proc_x = (ipx == 0)
    lfirst_proc_y = (ipy == 0)
    lfirst_proc_z = (ipz == 0)
    lfirst_proc_xy = lfirst_proc_x .and. lfirst_proc_y
    lfirst_proc_yz = lfirst_proc_y .and. lfirst_proc_z
    lfirst_proc_xz = lfirst_proc_x .and. lfirst_proc_z
    lfirst_proc_xyz = lfirst_proc_x .and. lfirst_proc_y .and. lfirst_proc_z
!
!  Set up flags for trailing processors in each possible direction and plane
!
    llast_proc_x = (ipx == nprocx-1)
    llast_proc_y = (ipy == nprocy-1)
    llast_proc_z = (ipz == nprocz-1)
    llast_proc_xy = llast_proc_x .and. llast_proc_y
    llast_proc_yz = llast_proc_y .and. llast_proc_z
    llast_proc_xz = llast_proc_x .and. llast_proc_z
    llast_proc_xyz = llast_proc_x .and. llast_proc_y .and. llast_proc_z
!
!  Set up directory names.
!
    call directory_names
!
!  Write data.
!
    call delete_file(trim(directory_snap)//'/'//filename)
    if (ldirect_access) then
      rec_len = int (mxgrid, kind=8) * int (mygrid, kind=8) * int (mz, kind=8) * mvar * io_len
      open (lun_output, FILE=trim(directory_snap)//'/'//filename, status='new', access='direct', recl=rec_len )
      write(lun_output, rec=1) f
      close(lun_output)
      open (lun_output, FILE=trim(directory_snap)//'/'//filename, status='old', form='unformatted', position='append')
    else
      open (lun_output, FILE=trim(directory_snap)//'/'//filename, status='new',form='unformatted')
      write(lun_output) f
    endif
!
    if (.not.lonly_farray) then
!
!  Write additional data.
!
      write(lun_output) t_sp
      if (lroot) write (lun_output) gx, gy, gz, dx, dy, dz
    endif
    close(lun_output)
!
  enddo
!
  close (lun_input)
  print *, 'Written snapshot for time t =', t
!
endsubroutine read_and_distribute
!*********************************************************************************************
endprogram pc_distribute_z
