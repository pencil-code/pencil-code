! This tool distributes a global data cube in into the proc-directories.
!
! $Id: pc_distribute.f90 17157 2011-07-02 23:59:13Z Bourdin.KIS $
!***********************************************************************
program pc_distribute
!
  use Cdata
  use Cparam, only: fnlen
  use Diagnostics
  use Filter
  use Grid, only: initialize_grid
  use IO
  use Messages
  use Param_IO
  use Register
  use Snapshot
  use Sub
!
  implicit none
!
  character (len=fnlen) :: filename
  character (len=*), parameter :: directory_in = 'data/allprocs'
!
  real, dimension (mx,my,mz,mfarray) :: f
  integer, parameter :: ngx=nxgrid+2*nghost, ngy=nygrid+2*nghost, ngz=nzgrid+2*nghost
  real, dimension (:,:,:,:), allocatable :: gf
  real, dimension (ngx) :: gx
  real, dimension (ngy) :: gy
  real, dimension (ngz) :: gz
  real :: dummy_dx, dummy_dy, dummy_dz
  logical :: ex
  integer :: mvar_in, bytes, pz, pa, alloc_err
  real :: t_sp   ! t in single precision for backwards compatibility
!
  lrun=.true.
  lmpicomm = .false.
  root = 0
  lroot = .true.
  ipx = 0
  ipy = 0
  ipz = 0
  ylneigh = 0
  zlneigh = 0
  yuneigh = 0
  zuneigh = 0
  llcorn = 0
  lucorn = 0
  uucorn = 0
  ulcorn = 0
!
  deltay = 0.0   ! Shearing not available due to missing fseek in Fortran
!
  inquire (IOLENGTH=bytes) 1.0
!
  if (lcollective_IO) call fatal_error ('pc_distribute', &
      "Distributing snapshots currently requires the distributed IO-module.")
!
  write (*,*) 'Please enter the filename to convert (eg. var.dat, VAR1, ...):'
  read (*,*) filename
!
!  Identify version.
!
  if (lroot) call svn_id( &
      '$Id: pc_downscale.f90 17157 2011-07-02 23:59:13Z Bourdin.KIS $')
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
!  Derived parameters (that may still be overwritten).
!  [might better be put into another routine, possibly even in rparam or
!  read_runpars]
!
  x0 = xyz0(1) ; y0 = xyz0(2) ; z0 = xyz0(3)
  Lx = Lxyz(1) ; Ly = Lxyz(2) ; Lz = Lxyz(3)
!
! Calculate dimensionality
!
  dimensionality=min(nxgrid-1,1)+min(nygrid-1,1)+min(nzgrid-1,1)
!
!  Register physics modules.
!
  call register_modules()
!
!  Define the lenergy logical
!
  lenergy=lentropy.or.ltemperature.or.lthermal_energy
!
  if (lwrite_aux .and. .not. lread_aux) then
    if (lroot) then
      print *, ''
      print *, 'lwrite_aux=T but lread_aux=F'
      print *, 'The code will write the auxiliary variables to allprocs/VARN'
      print *, ' without having read them from proc*/VARN'
      print *, ''
      call fatal_error("pc_distribute","Stop and check")
    endif
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
  allocate (gf (ngx,ngy,mz,mvar_io), stat=alloc_err)
  if (alloc_err /= 0) call fatal_error ('pc_distribute', 'Failed to allocate memory for gf.', .true.)
!
!  Print resolution and dimension of the simulation.
!
  if (lroot) write (*,'(a,i1,a)') ' This is a ', dimensionality, '-D run'
  if (lroot) print *, 'nxgrid, nygrid, nzgrid=', nxgrid, nygrid, nzgrid
  if (lroot) print *, 'Lx, Ly, Lz=', Lxyz
  if (lroot) print *, '      Vbox=', Lxyz(1)*Lxyz(2)*Lxyz(3)
!
  inquire (file=trim(directory_in)//'/'//filename, exist=ex)
  if (.not. ex) call fatal_error ('pc_distribute', 'File not found: '//trim(directory_in)//'/'//filename, .true.)
  inquire (file=trim(directory_in)//'/grid.dat', exist=ex)
  if (.not. ex) call fatal_error ('pc_distribute', 'File not found: '//trim(directory_in)//'/grid.dat', .true.)
!
  ! read time:
  open (lun_input, FILE=trim(directory_in)//'/grid.dat', FORM='unformatted', status='old')
  read (lun_input) t_sp, gx, gy, gz, dummy_dx, dummy_dy, dummy_dz
  if (lshear) read (lun_input) deltay
  close (lun_input)
  t = t_sp
!
  open (lun_input, FILE=trim(directory_in)//'/'//filename, access='direct', recl=ngx*ngy*bytes, status='old')
!
!  Allow modules to do any physics modules do parameter dependent
!  initialization. And final pre-timestepping setup.
!  (must be done before need_XXXX can be used, for example)
!
  lpencil_check_at_work = .true.
  call initialize_modules (f, LSTARTING=.true.)
  lpencil_check_at_work = .false.
!
! Loop over processors
!
  write (*,*) "IPZ-layer:"
!
  do ipz = 0, nprocz-1
!
    write (*,*) ipz+1, " of ", nprocz
!
    f = huge(1.0)
    gf = huge(1.0)
!
    ! read xy-layer:
    do pa = 1, mvar_io
      do pz = 1, mz
        read (lun_input, rec=pz+ipz*nz+(pa-1)*ngz) gf(:,:,pz,pa)
      enddo
    enddo
!
    do ipy = 0, nprocy-1
      do ipx = 0, nprocx-1
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
        call directory_names()
!
!  Read coordinates.
!
        if (ip<=6.and.lroot) print*, 'reading grid coordinates'
        call rgrid ('grid.dat')
!
! Size of box at local processor. The if-statement is for 
! backward compatibility.
!
        if (all(lequidist)) then 
          Lxyz_loc(1)=Lxyz(1)/nprocx
          Lxyz_loc(2)=Lxyz(2)/nprocy
          Lxyz_loc(3)=Lxyz(3)/nprocz
          xyz0_loc(1)=xyz0(1)+ipx*Lxyz_loc(1) 
          xyz0_loc(2)=xyz0(2)+ipy*Lxyz_loc(2)
          xyz0_loc(3)=xyz0(3)+ipz*Lxyz_loc(3)
          xyz1_loc(1)=xyz0_loc(1)+Lxyz_loc(1)
          xyz1_loc(2)=xyz0_loc(2)+Lxyz_loc(2)
          xyz1_loc(3)=xyz0_loc(3)+Lxyz_loc(3)
        else
          xyz0_loc(1)=x(l1)
          xyz0_loc(2)=y(m1)
          xyz0_loc(3)=z(n1)
          xyz1_loc(1)=x(l2)
          xyz1_loc(2)=y(m2)
          xyz1_loc(3)=z(n2)
          Lxyz_loc(1)=xyz1_loc(1) - xyz0_loc(1)
          Lxyz_loc(2)=xyz1_loc(2) - xyz0_loc(3)
          Lxyz_loc(3)=xyz1_loc(3) - xyz0_loc(3)
        endif
!
!  Need to re-initialize the local grid for each processor.
!
        call initialize_grid()
!
        ! distribute gf to f:
        f(:,:,:,1:mvar_io) = gf(1+ipx*nx:mx+ipx*nx,1+ipy*ny:my+ipy*ny,:,:)
        x = gx(1+ipx*nx:mx+ipx*nx)
        y = gy(1+ipy*ny:my+ipy*ny)
        z = gz(1+ipz*nz:mz+ipz*nz)
!
        ! write data:
        call wsnap (filename, f(:,:,:,1:mvar_io), mvar_io, enum=.false., noghost=.true.)
!
      enddo
    enddo
  enddo
!
  close (lun_input)
  print *, 'Writing snapshot for time t =', t
!
!  Gvie all modules the possibility to exit properly.
!
  call finalize_modules (f, .true.)
!
!  Free any allocated memory.
!
  call fnames_clean_up()
  call vnames_clean_up()
!
endprogram pc_distribute
