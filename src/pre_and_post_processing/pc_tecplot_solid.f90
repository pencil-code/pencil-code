! This is a tool to collect a distributed data cube and record them in 
! tecplot 
! files. The files named "OutputN.dat"(N=0~iproc_world) means this tecplot format
! data file is converted from var.dat in directory procN.
! This file is based on the pc_collect.f90 writen by .
!
! $Id: pc_tecplot.f90 10000 2014-04-17 22:12:50Z Zhuang $
!
program pc_tecplot_solid
!
  use Cdata
  use Cparam, only: fnlen, intlen
  use Diagnostics
  use Filter
  use General, only: itoa
  use Grid, only: initialize_grid,construct_grid,set_coorsys_dimmask
  use IO
  use Messages
  use Param_IO
  use Register
  use Snapshot
  use Sub
  use Syscalls, only: sizeof_real
!
  implicit none
!
  character (len=fnlen) :: filename
  character (len=4) :: chproc
  character (len=*), parameter :: directory_out = 'allprocs'
!
  real, dimension(mx,my,mz,mfarray) :: f
  real, dimension(:,:,:,:), allocatable :: gf
  real, dimension(:,:,:), allocatable :: sum_Y
  real, dimension (mxgrid) :: gx, gy, gz
  integer   ::  i, j, k, kchem, iroot
  logical :: ex
  integer :: mvar_in, io_len, alloc_err, start_pos, end_pos
  integer(kind=8) :: rec_len
  real :: t_sp, t_test   ! t in single precision for backwards compatibility
!
  lstart=.true.
  lmpicomm = .false.
  iroot = 0
  lroot = .true.
  ipx = 0
  ipy = 0
  ipz = 0
  ylneigh = 0
  zlneigh = 0
  yuneigh = 0
  zuneigh = 0
!
  inquire (IOLENGTH=io_len) 1.0
!
  if (IO_strategy == "collect") call fatal_error ('pc_tecplot', &
      "Snapshots are already collected, when using the 'io_collect' module.")
  if (IO_strategy == "MPI-IO") call fatal_error ('pc_tecplot', &
      "Snapshots are already collected, when using the MPI-IO module.")
!
  write (*,*) 'Please enter the filename to convert (eg. var.dat, VAR1, ...):'
  read (*,*) filename
!
!  Identify version.
!
  if (lroot) call svn_id( &
      '$Id: pc_tecplot.f90 18690 2014-04-16 22:12:50Z Zhuang $')
!
!  Initialize the message subsystem, eg. color setting etc.
!
  call initialize_messages
!
!  Read parameters from start.x (default values; may be overwritten by
!  read_runpars).
!
  call read_all_init_pars
  call set_coorsys_dimmask
  lstart=.false.; lrun=.true.
!
!  Read parameters and output parameter list.
!
  call read_all_run_pars
!
!  Derived parameters (that may still be overwritten).
!  [might better be put into another routine, possibly even in rparam or
!  read_runpars]
!
  x0 = xyz0(1) ; y0 = xyz0(2) ; z0 = xyz0(3)
  Lx = Lxyz(1) ; Ly = Lxyz(2) ; Lz = Lxyz(3)
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
    if (lroot) then
      print *, ''
      print *, 'lwrite_aux=T but lread_aux=F'
      print *, 'The code will write the auxiliary variables to allprocs/VARN'
      print *, ' without having read them from proc*/VARN'
      print *, ''
      call fatal_error("pc_tecplot","Stop and check")
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
  allocate (gf (mxgrid,mygrid,mz,mvar_io), stat=alloc_err)
  allocate (sum_Y(mxgrid,mygrid,mz))
  if (alloc_err /= 0) call fatal_error ('pc_tecplot', 'Failed to allocate memory for gf.', .true.)
!
!  Print resolution and dimension of the simulation.
!
  if (lroot) write (*,'(a,i1,a)') ' This is a ', dimensionality, '-D run'
  if (lroot) print *, 'nxgrid, nygrid, nzgrid=', nxgrid, nygrid, nzgrid
  if (lroot) print *, 'Lx, Ly, Lz=', Lxyz
  if (lroot) print *, 'Vbox=', Lxyz(1)*Lxyz(2)*Lxyz(3)
  if (lroot) write (*,*) 'mvar = ', mvar_io
!
  iproc_world = 0
  call directory_names
  inquire (file=trim(directory_dist)//'/'//filename, exist=ex)
  if (.not. ex) call fatal_error ('pc_tecplot', 'File not found: '//trim(directory_dist)//'/'//filename, .true.)
!
!  Allow modules to do any physics modules do parameter dependent
!  initialization. And final pre-timestepping setup.
!  (must be done before need_XXXX can be used, for example)
!
  call construct_grid(x,y,z,dx,dy,dz)
  call initialize_modules(f)
!
  if (IO_strategy == "collect_xy") then
!
    gz = huge(1.0)
    t_test = huge(1.0)
!
    do ipz = 0, nprocz-1
      write (*,*) "IPZ-layer:"
      write (*,*) ipz+1, " of ", nprocz
!
      f = huge(1.0)
      gf = huge(1.0)
!
      iproc_world = ipz * nprocx*nprocy
      lroot = (iproc_world==iroot)
      lfirst_proc_z = (ipz == 0)
      llast_proc_z = (ipz == nprocz-1)
!
! Take the shortcut, files are well prepared for direct combination
! Set up directory names 'directory' and 'directory_snap'
!
      call directory_names
!
! Read the data
!
      rec_len = int (mxgrid, kind=8) * int (mygrid, kind=8) * mz
      rec_len = rec_len * mvar_in * io_len
      open (lun_input, FILE=trim (directory_snap)//'/'//filename, access='direct', recl=rec_len, status='old')
      read (lun_input, rec=1) gf
      close (lun_input)
!
! Write data in tecplot data format
!
      write(chproc,'(I4)') iproc_world
      open (lun_output, FILE=trim(datadir)//trim(directory_out)//'/'//'output'//chproc//'.dat', status='replace', FORM='FORMATTED')
      write(lun_output,*) 'TITLE="Output"'
      write(lun_output,*) 'variables="x","y","z","u","v","w","rho"'
      if (dimensionality==2) then
        write(lun_output,*) 'zone,','I=',mxgrid,',J=',mygrid,',K=',nz,',F=POINT'
      elseif (dimensionality==3) then
        write(lun_output,*) 'zone,','I=',mxgrid,',J=',mygrid,',K=',mzgrid,',F=POINT' 
      endif
!
      if (dimensionality==2) then
        start_pos = nghost + 1
        end_pos = nghost + nz
      elseif (dimensionality==3) then
        start_pos = 1
        end_pos = mz
      endif
      do k = start_pos, end_pos
      do j=1,mygrid
      do i=1,mxgrid
          write (lun_output, 101) x(i),y(j),z(k),gf(i,j,k,1),gf(i,j,k,2),gf(i,j,k,3),gf(i,j,k,4)
101 format(1X,7e16.8)
      enddo
      enddo
      enddo
      close (lun_output)
    enddo
  endif
!
  if (IO_strategy == "dist") then
!
    do ipz = 0, nprocz-1
    do ipy = 0, nprocy-1
    do ipx = 0, nprocx-1
!
      iproc_world = ipx + ipy * nprocx + ipz * nprocx*nprocy
      lroot = (iproc_world==iroot)
      gx = huge(1.0)
      gy = huge(1.0)
      f  = huge(1.0)
      gf = huge(1.0)
!
!  Set up directory names `directory' and `directory_snap'.
!
      call directory_names
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
      call construct_grid(x,y,z,dx,dy,dz)
!
!  Read data.
!  Snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!
      call rsnap (filename, f(:,:,:,1:mvar_in), mvar_in, lread_nogrid)
      t_sp = t
      write (*,*) 'Reading: '//trim(directory_snap)//'/'//trim(filename)//'.'
!
      if (lroot) t_test = t_sp
      if (t_test /= t_sp) then
        write (*,*) 'ERROR: '//trim(directory_snap)//'/'//trim(filename)//' IS INCONSISTENT: t=', t_sp
        stop 1
      endif
!
! collect f in gf:
!
      gf(1+ipx*nx:mx+ipx*nx,1+ipy*ny:my+ipy*ny,1+ipz*nz:mz+ipz*nz,:) = f(:,:,:,1:mvar_io)
      if (.not. ldensity_nolog) then
        gf(1+ipx*nx:mx+ipx*nx,1+ipy*ny:my+ipy*ny,1+ipz*nz:mz+ipz*nz,ilnrho)= &
          exp(gf(1+ipx*nx:mx+ipx*nx,1+ipy*ny:my+ipy*ny,1+ipz*nz:mz+ipz*nz,ilnrho))
      endif
      if (.not. ltemperature_nolog) then
        gf(1+ipx*nx:mx+ipx*nx,1+ipy*ny:my+ipy*ny,1+ipz*nz:mz+ipz*nz,ilnTT)= &
          exp(gf(1+ipx*nx:mx+ipx*nx,1+ipy*ny:my+ipy*ny,1+ipz*nz:mz+ipz*nz,ilnTT))
      endif
!
      if (lchemistry) then
        sum_Y(1+ipx*nx:mx+ipx*nx,1+ipy*ny:my+ipy*ny,1+ipz*nz:mz+ipz*nz)=0.0
        do kchem=1, nchemspec
          sum_Y(1+ipx*nx:mx+ipx*nx,1+ipy*ny:my+ipy*ny,1+ipz*nz:mz+ipz*nz)= &
            sum_Y(1+ipx*nx:mx+ipx*nx,1+ipy*ny:my+ipy*ny,1+ipz*nz:mz+ipz*nz)+f(:,:,:,ichemspec(kchem))
        enddo
      endif
      write(chproc,'(I4)') iproc_world
      open (lun_output, FILE=trim(datadir)//'/'//trim(directory_out)//'/'//'output'//chproc//'.dat',&
           status='replace', FORM='FORMATTED')
      write(lun_output,*) 'TITLE="Output"'
!
      if (ilnTT>0 .and. lchemistry) then
        write(lun_output,*) 'variables="x","y","z","u","v","w","rho","T","N2","O2","CO","CO2","Sum"'
      elseif (ilnTT>0) then
        write(lun_output,*) 'variables="x","y","z","u","v","w","rho","T"'
      else
        write(lun_output,*) 'variables="x","y","z","u","v","w","rho"'
      endif
!
      if (dimensionality==2) then
        write(lun_output,*) 'zone,','I=',mx,',J=',my,',K=',nz,',F=POINT'
      elseif (dimensionality==3) then
        write(lun_output,*) 'zone,','I=',mx,',J=',my,',K=',mz,',F=POINT'
      endif
!
      if (dimensionality==2) then
        start_pos=1+nghost+ipz*nz
        end_pos=mz+ipz*nz-nghost
      elseif (dimensionality==3) then
        start_pos=1+ipz*nz
        end_pos=mz+ipz*nz
      endif
!
      do k = start_pos, end_pos
      do j = 1+ipy*ny, my+ipy*ny
      do i = 1+ipx*nx, mx+ipx*nx
        write (lun_output,102) x(i-ipx*nx),y(j-ipy*ny),z(k-ipz*nz),gf(i,j,k,1:mvar),sum_Y(i,j,k)
      enddo
      enddo
      enddo
102 format(1X,13e16.8)
      close (lun_output)
    enddo
    enddo
    enddo
  endif
!
  print *, 'Writing snapshot for time t =', t
!
!  Give all modules the possibility to exit properly.
!
  call finalize_modules (f)
!
!  Free any allocated memory.
!
  deallocate (gf)
  call fnames_clean_up
  call vnames_clean_up
!
end program pc_tecplot_solid
