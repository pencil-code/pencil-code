! This is a very feature-limited tool to downscale a data cube in the
! horzontal directions (x and y), if these are periodic.
!
! $Id$
!***********************************************************************
program pc_downscale
!
  use Cdata
  use Diagnostics
  use Filter
  use IO
  use Messages
  use Param_IO
  use Register
  use Snapshot
  use Sub
!
  implicit none
!
  integer, parameter :: reduce=2
  character (len=80) :: filename
  character (len=*), parameter :: directory_out = 'data/allprocs'
!
  real, dimension (mx,my,mz,mfarray) :: f
!!!  real, dimension (mx,my,mz,mvar) :: df
  integer, parameter :: nrx=nxgrid/reduce+2*nghost, nry=nygrid/reduce+2*nghost, nrz=nzgrid+2*nghost
  real, dimension (nrx,nry,mz,mfarray) :: rf
  real, dimension (nrx) :: rx
  real, dimension (nry) :: ry
!!!  type (pencil_case) :: p
  integer :: mvar_in, bytes, px, py, pz, pa, start_pos, end_pos
  real, parameter :: inv_reduce_2 = 1.0 / reduce**2.0
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
  inquire (IOLENGTH=bytes) 1.0
!
  write (*,*) 'Please enter the filename to convert (eg. var.dat, VAR1, ...):'
  read (*,*) filename
!
!  Identify version.
!
  if (lroot) call svn_id( &
      '$Id$')
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
  if (.not. lperi(1) .and. (reduce /= 1)) call fatal_error ('run', 'reduction impossible in X: not periodic')
  if (.not. lperi(2) .and. (reduce /= 1)) call fatal_error ('run', 'reduction impossible in Y: not periodic')
  if (mod (nx, reduce) /= 0) call fatal_error ('run', 'NX not dividable by reduce factor')
  if (mod (ny, reduce) /= 0) call fatal_error ('run', 'NY not dividable by reduce factor')
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
! Loop over processors
!
  do ipz = 0, nprocz-1
!
    rf = huge(1.0)
    rx = huge(1.0)
    ry = huge(1.0)
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
!  Set up directory names `directory' and `directory_snap'.
!
        call directory_names()
!
!  Read coordinates.
!
        if (ip<=6.and.lroot) print*, 'reading grid coordinates'
        call rgrid(trim(directory)//'/grid.dat')
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
!  Read data.
!  Snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!  NOTE: for io_dist, rtime doesn't read the time, only for io_mpio.
!
        call rsnap(trim(directory_snap)//'/'//trim(filename),f,mvar_in)
!
!  Read time and global variables (if any).
!
        call rtime(trim(directory)//'/'//trim(filename),t)
        if (mglobal/=0)  &
            call input_globals(trim(directory_snap)//'/global.dat', &
            f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
!
!  Allow modules to do any physics modules do parameter dependent
!  initialization. And final pre-timestepping setup.
!  (must be done before need_XXXX can be used, for example)
!
        call initialize_modules(f,LSTARTING=.true.)
!
!  Find out which pencils are needed and write information about required,
!  requested and diagnostic pencils to disc.
!
        call choose_pencils()
!
        ! downscale f:
        do pa = 1, mfarray
          do pz = 1, mz
            do py = 0, ny-1, reduce
              do px = 0, nx-1, reduce
                rf(nghost+1+(px+ipx*nx)/reduce,nghost+1+(py+ipy*ny)/reduce,pz,pa) = &
                    sum (f(nghost+1+px:nghost+px+reduce,nghost+1+py:nghost+py+reduce,pz,pa)) * inv_reduce_2
              enddo
            enddo
          enddo
        enddo
!
        ! downscale x coordinates:
        do px = 0, nx-1, reduce
          rx(nghost+1+(px+ipx*nx)/reduce) = sum (x(nghost+1+px:nghost+px+reduce)) / reduce
        enddo
!
        ! downscale y coordinates:
        do py = 0, ny-1, reduce
          ry(nghost+1+(py+ipy*ny)/reduce) = sum (y(nghost+1+py:nghost+py+reduce)) / reduce
        enddo
!
      enddo
    enddo
!
    ! communicate ghost cells along the y direction:
    rf(nghost+1:nrx-nghost,           1:nghost,  :,:) = rf(nghost+1:nrx-nghost,nry-2*nghost+1:nry-nghost,:,:)
    rf(nghost+1:nrx-nghost,nry-nghost+1:nry,     :,:) = rf(nghost+1:nrx-nghost,      nghost+1:2*nghost,  :,:)
    ry(           1:nghost) = ry(nry-2*nghost+1:nry-nghost)
    ry(nry-nghost+1:nry   ) = ry(      nghost+1:2*nghost  )
!
    ! communicate ghost cells along the x direction:
    rf(           1:nghost,:,:,:) = rf(nrx-2*nghost+1:nrx-nghost,:,:,:)
    rf(nrx-nghost+1:nrx,   :,:,:) = rf(      nghost+1:2*nghost,  :,:,:)
    rx(           1:nghost) = rx(nrx-2*nghost+1:nrx-nghost)
    rx(nrx-nghost+1:nrx   ) = rx(      nghost+1:2*nghost  )
!
    ! write xy-layer:
    open(lun_output,FILE=trim(directory_out)//'/'//trim(filename),FORM='unformatted',access='direct',recl=nrx*nry*bytes)
    do pa = 1, mfarray
      start_pos = nghost + 1
      end_pos = nz
      if (ipz == 0) start_pos = 1
      if (ipz == nprocz-1) end_pos = mz
      do pz = start_pos, end_pos
        write(lun_output,rec=pz+ipz*nz+(pa-1)*nrz) rf(:,:,pz,pa)
      enddo
    enddo
    close(lun_output)
!
  enddo
!
  ! write additional data:
  open(lun_output,FILE=trim(directory_out)//'/'//trim(filename),FORM='unformatted',position='append')
  t_sp = t
  if (lshear) then
    write(lun_output) t_sp,rx,ry,z,dx*reduce,dy*reduce,dz,deltay
  else
    write(lun_output) t_sp,rx,ry,z,dx*reduce,dy*reduce,dz
  endif
  close(lun_output)
!
  print*, 'Writing snapshot for time t =', t
!
!  Free any allocated memory.
!
  call fnames_clean_up()
  call vnames_clean_up()
!
endprogram pc_downscale
