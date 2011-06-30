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
  use Mpicomm, only: initialize_mpicomm
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
  integer, parameter :: nrx=nxgrid/reduce+2*nghost, nry=nygrid/reduce+2*nghost
  real, dimension (nrx,nry,mz,mfarray) :: rf
  real, dimension (nrx) :: rx
  real, dimension (nry) :: ry
!!!  type (pencil_case) :: p
  integer :: mvar_in, bytes, px, py, pz, pa
  integer :: xb, xe, xb_loc, xe_loc, yb, ye, yb_loc, ye_loc, zb, ze, zb_loc, ze_loc
  real, parameter :: inv_reduce_2 = 1.0 / reduce**2.0
  real :: t_sp   ! t in single precision for backwards compatibility
!
  lrun=.true.
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
! Loop over processors
!
  do ipz = 0, nprocz-1
!
    rf = 0.0
!
    do ipy = 0, nprocy-1
      do ipx = 0, nprocx-1
!
    iproc = ipx + ipy * nprocx + ipz * nprocx*nprocy
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
!  Inform about verbose level.
!
  if (lroot) print*, 'The verbose level is ip=', ip, ' (ldebug=', ldebug, ')'
!
!  Call rprint_list to initialize diagnostics and write indices to file.
!
  call rprint_list(LRESET=.false.)
!
!  Position of equator (if any).
!
  if (lequatory) yequator=xyz0(2)+0.5*Lxyz(2)
  if (lequatorz) zequator=xyz0(3)+0.5*Lxyz(3)
!
!  Limits to xaveraging.
!
  if (lav_smallx) call init_xaver
!
!  Inner radius for freezing variables defaults to r_min.
!  Note: currently (July 2005), hydro.f90 uses a different approach:
!  r_int will override rdampint, which doesn't seem to make much sense (if
!  you want rdampint to be overridden, then don't specify it in the first
!  place).
!
  if (rfreeze_int==-impossible .and. r_int>epsi) rfreeze_int=r_int
  if (rfreeze_ext==-impossible) rfreeze_ext=r_ext
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
!  Get state length of random number generator.
!
  call get_nseed(nseed)
!
!  Read data.
!  Snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!  NOTE: for io_dist, rtime doesn't read the time, only for io_mpio.
!
  call rsnap(trim(directory_snap)//'/'//trim(filename),f,mvar_in)
!
  call get_nseed(nseed)
!
!  Read time and global variables (if any).
!
  call rtime(trim(directory)//'/'//trim(filename),t)
  if (mglobal/=0)  &
      call input_globals(trim(directory_snap)//'/global.dat', &
      f(:,:,:,mvar+maux+1:mvar+maux+mglobal),mglobal)
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
!  Find out which pencils are needed and write information about required,
!  requested and diagnostic pencils to disc.
!
  call choose_pencils()
!
      t_sp = t
!
      if (lfirst_proc_x) then
        xb = 0
        xe = xb + mx - 1
        xb_loc = 0
        xe_loc = mx - 1
      else
        xb = ipx * nx + nghost
        xe = xb + mx - 1 - nghost
        xb_loc = nghost
        xe_loc = mx - 1
      endif
!
      if (lfirst_proc_y) then
        yb = 0
        ye = yb + my - 1
        yb_loc = 0
        ye_loc = my - 1
      else
        yb = ipy * ny + nghost
        ye = yb + my - 1 - nghost
        yb_loc = nghost
        ye_loc = my - 1
      endif
!
      if (lfirst_proc_z) then
        zb = 0
        ze = zb + mz - 1
        zb_loc = 0
        ze_loc = mz - 1
      else
        zb = ipz * nz + nghost
        ze = zb + mz - 1 - nghost
        zb_loc = nghost
        ze_loc = mz - 1
      endif
!
        ! downscale f:
        do pa = 1, mfarray
          do pz = 1, mz
            do py = 1, nry
              do px = 1, nrx
                rf(px,py,pz,pa) = sum (f(nghost+(px-1)*reduce+1:nghost+px*reduce,nghost+(py-1)*reduce+1:nghost+py*reduce,pz,pa)) &
                    * inv_reduce_2
              enddo
            enddo
          enddo
        enddo
!
        ! downscale x coordinates:
        do px = 1, nrx
          rx(px) = sum (x(nghost+(px-1)*reduce+1:nghost+px*reduce)) / reduce**2.0
        enddo
!
        ! downscale y coordinates:
        do py = 1, nry
          ry(py) = sum (y(nghost+(py-1)*reduce+1:nghost+py*reduce)) / reduce**2.0
        enddo
!
      enddo
    enddo
!
    ! communicate ghost cells along the y direction:
    rf(nghost+1:mx-nghost,1:nghost,:,:) = rf(nghost+1:mx-nghost,my-2*nghost+1:my-nghost,:,:)
    rf(nghost+1:mx-nghost,nghost+1:2*nghost,:,:) = rf(nghost+1:mx-nghost,my-nghost+1:my,:,:)
    ! communicate ghost cells along the x direction:
    rf(1:nghost,:,:,:) = rf(mx-2*nghost+1:mx-nghost,:,:,:)
    rf(nghost+1:2*nghost,:,:,:) = rf(mx-nghost+1:mx,:,:,:)
!
    ! write xy-layer:
    open(lun_output,FILE=trim(directory_out)//'/'//trim(filename),FORM='unformatted',access='direct',recl=nrx*nry*bytes)
    do pa = 1, mfarray
      write(lun_output,rec=ipz+pa*mz) rf(:,:,:,pa)
    enddo
    close(lun_output)
!
  enddo
!
  ! write additional data:
  open(lun_output,FILE=trim(directory_out)//'/'//trim(filename),FORM='unformatted',position='append')
  if (lshear) then
    write(lun_output) t_sp,rx,ry,z,dx,dy,dz,deltay
  else
    write(lun_output) t_sp,rx,ry,z,dx,dy,dz
  endif
  close(lun_output)
!
  print*, 'Writing final snapshot at time t =', t
!
!  Free any allocated memory.
!
  call fnames_clean_up()
  call vnames_clean_up()
!
endprogram pc_downscale
