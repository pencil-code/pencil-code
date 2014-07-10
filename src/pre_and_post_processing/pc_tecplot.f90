! This is a tool to collect a distributed data cube and record them in tecplot
! files. The files named "OutputN.dat"(N=0~iproc) means this tecplot format
! data file is converted from var.dat in directory procN.
! This file is based on the pc_collect.f90 writen by Bourdin.KIS.
!
! $Id$
!***********************************************************************
program pc_collect
!
  use Cdata
  use Cparam, only: fnlen, intlen
  use Diagnostics
  use Filter
  use General, only: itoa
  use Grid, only: initialize_grid
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
  character (len=intlen) :: chproc
  character (len=*), parameter :: directory_out = 'data/allprocs'
!
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (:,:,:,:), allocatable :: gf
  real, dimension (mxgrid) :: gx, gy, gz
  integer   ::  i, j, k
  logical :: ex
  integer :: mvar_in, io_len, pz, alloc_err, start_pos, end_pos
  integer(kind=8) :: rec_len, num_rec
  real :: t_sp, t_test   ! t in single precision for backwards compatibility
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
      '$Id$')
!
!  Initialize the message subsystem, eg. color setting etc.
!
  call initialize_messages()
!
!  Read parameters from start.x (default values; may be overwritten by
!  read_runpars).
!
  call read_startpars()
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
  dimensionality = min(nxgrid-1,1) + min(nygrid-1,1) + min(nzgrid-1,1)
!
!  Register physics modules.
!
  call register_modules()
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
  iproc = 0
  call directory_names()
  inquire (file=trim(directory_dist)//'/'//filename, exist=ex)
  if (.not. ex) call fatal_error ('pc_tecplot', 'File not found: '//trim(directory_dist)//'/'//filename, .true.)
!
!  Allow modules to do any physics modules do parameter dependent
!  initialization. And final pre-timestepping setup.
!  (must be done before need_XXXX can be used, for example)
!
  lpencil_check_at_work = .true.
  call initialize_modules(f, LSTARTING=.true.)
  lpencil_check_at_work = .false.
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
    iproc = ipz * nprocx*nprocy
    lroot = (iproc==root)
    lfirst_proc_z = (ipz == 0)
    llast_proc_z = (ipz == nprocz-1)
!
! Take the shortcut, files are well prepared for direct combination
! Set up directory names 'directory' and 'directory_snap'
!
      call directory_names()
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
      chproc=itoa(ipz)
      open (lun_output, FILE=trim(directory_out)//'/'//'output'//chproc//'.dat', status='replace', FORM='FORMATTED')
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
101     format(1X,7e16.8)
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
    iproc = ipx + ipy * nprocx + ipz * nprocx*nprocy
    write (*,*) "now reading",iproc
    lroot = (iproc==root)
    gx = huge(1.0)
    gy = huge(1.0)
    f  = huge(1.0)
    gf = huge(1.0)
!
!  Set up directory names `directory' and `directory_snap'.
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
!  Read data.
!  Snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!
    call rsnap (filename, f(:,:,:,1:mvar_in), mvar_in)
    t_sp = t
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
!
    chproc=itoa(iproc)
    open (lun_output, FILE=trim(directory_out)//'/'//'output'//chproc//'.dat', status='replace', FORM='FORMATTED')
    write(lun_output,*) 'TITLE="Output"'
    write(lun_output,*) 'variables="x","y","z","u","v","w","rho"'
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
      write (lun_output,101) x(i-ipx*nx),y(j-ipy*ny),z(k-ipz*nz),&
        gf(i,j,k,1),gf(i,j,k,2),gf(i,j,k,3),gf(i,j,k,4)
    enddo
    enddo
    enddo
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
  call finalize_modules (f, .true.)
!
!  Free any allocated memory.
!
  deallocate (gf)
  call fnames_clean_up()
  call vnames_clean_up()
!
endprogram pc_collect
