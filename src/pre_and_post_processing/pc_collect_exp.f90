! This is a tool to collect a distributed data cube in one file.
!
! $Id: pc_collect.f90 22542 2014-11-18 22:51:03Z Bourdin.KIS $
!***********************************************************************
program pc_collect
!
  use Cdata
  use Cparam, only: fnlen
  !!use Diagnostics
  !!use Filter
  !!use Grid, only: initialize_grid
  use IO
  use Messages
  use Param_IO
  !!use Register
  use Snapshot
  !!use Sub
  !!use Syscalls, only: sizeof_real
!
  implicit none
!
  character (len=fnlen) :: filename
  character (len=*), parameter :: directory_out = 'data/allprocs'
!
  real, dimension (:,:,:,:), allocatable :: gf
  real, dimension (mxgrid) :: gx, gdx_1, gdx_tilde
  real, dimension (mygrid) :: gy, gdy_1, gdy_tilde
  real, dimension (mzgrid) :: gz, gdz_1, gdz_tilde
  logical :: ex
  integer :: mvar_in, io_len, pz, pa, start_pos, end_pos, alloc_err
  integer(kind=8) :: rec_len, num_rec
  real :: t_sp, t_test   ! t in single precision for backwards compatibility
!
  lstart = .true.
  lmpicomm = .false.
  root = 0
!
  inquire (IOLENGTH=io_len) 1.0
!
  if (IO_strategy == "collect") call fatal_error ('pc_collect', &
      "Snapshots are already collected, when using the 'io_collect' module.")
  if (IO_strategy == "MPI-IO") call fatal_error ('pc_collect', &
      "Snapshots are already collected, when using the MPI-IO module.")
!
  write(*,*) 'Please enter the filename to convert (eg. var.dat, VAR1, ...):'
  read (*,*) filename
!
!  Identify version.
!
   call svn_id('$Id: pc_collect.f90 22542 2014-11-18 22:51:03Z Bourdin.KIS $')
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
! Calculate dimensionality
!
  dimensionality = min(nxgrid-1,1) + min(nygrid-1,1) + min(nzgrid-1,1)
!
  if (lwrite_aux .and. .not. lread_aux) then
    print *, ''
    print *, 'lwrite_aux=T but lread_aux=F'
    print *, 'The code will write the auxiliary variables to allprocs/VARN'
    print *, ' without having read them from proc*/VARN'
    print *, ''
    call fatal_error("pc_collect","Stop and check")
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
  allocate(gf(mxgrid,mygrid,mz,mvar_io), stat=alloc_err)
  if (alloc_err /= 0) call fatal_error ('pc_collect', 'Failed to allocate memory for gf.', .true.)
!
!  Print resolution and dimension of the simulation.
!
  write (*,'(a,i1,a)') ' This is a ', dimensionality, '-D run'
  print *, 'nxgrid, nygrid, nzgrid=', nxgrid, nygrid, nzgrid
  print *, 'Lx, Ly, Lz=', Lxyz
  print *, '      Vbox=', Lxyz(1)*Lxyz(2)*Lxyz(3)
!
  call directory_names()
  inquire (file=trim(directory_dist)//'/'//filename, exist=ex)
  if (.not. ex) call fatal_error ('pc_collect', 'File not found: '//trim(directory_dist)//'/'//filename, .true.)
  open (lun_output, FILE=trim(directory_out)//'/'//filename, status='replace', access='direct', recl=mxgrid*mygrid*io_len)
!
!  Allow modules to do any physics modules do parameter dependent
!  initialization. And final pre-timestepping setup.
!  (must be done before need_XXXX can be used, for example)
!
  !!!call initialize_modules(f)
!
! Loop over processors
!
  write (*,'("IPZ-layer:")')
!
  gf = 0.

  do ipz = 0, nprocz-1
!
    write (*,*) ipz+1, " of ", nprocz
!
    lfirst_proc_z = (ipz == 0)
    llast_proc_z  = (ipz == nprocz-1)
!
    if (IO_strategy == "collect_xy") then
!
      iproc = ipz * nprocx*nprocy
      lroot = (iproc==root)
!
! Take the shortcut, files are well prepared for direct combination
!
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
! Read additional information and check consistency of timestamp
!      
      open (lun_input, FILE=trim (directory_snap)//'/'//filename, FORM='unformatted', status='old', position='append')
      backspace(lun_input); if (lroot) backspace(lun_input)
      read (lun_input) t_sp
!
      if (lroot) then
        t_test = t_sp
        read (lun_input) gx, gy, gz, dx, dy, dz
      elseif (t_test /= t_sp) then
        write (*,*) 'ERROR: '//trim(directory_snap)//'/'//trim(filename)//' IS INCONSISTENT: t=', t_sp
        stop 1
      endif
      close (lun_input)

      t = t_sp
!
! Write xy-layer
!
      do pa = 1, mvar_io
!
        if (lfirst_proc_z) then
          start_pos = 1
        else
          start_pos = nghost + 1
        endif
        if (llast_proc_z) then
          end_pos = mz
        else
          end_pos = nghost + nz
        endif
!
        do pz = start_pos, end_pos
          write (lun_output, rec=pz+ipz*nz+(pa-1)*mzgrid) gf(:,:,pz,pa)
        enddo
!
      enddo
!
! That is all we have to do
!
      cycle
    endif
!
    do ipy = 0, nprocy-1
      do ipx = 0, nprocx-1
!
        iproc = ipx + ipy * nprocx + ipz * nprocx*nprocy
        lroot = (iproc==root)

! Set up directory names 'directory' and 'directory_snap'

        call directory_names
!
!  Read coordinates.
!
        if (ip<=6.and.lroot) print*, 'reading grid coordinates'
        call rgrid('grid.dat')
!
!  Need to re-initialize the local grid for each processor.
!
!!!        call initialize_grid()
!
!  Read data.
!  Snapshot data are saved in the tmp subdirectory.
!  This directory must exist, but may be linked to another disk.
!
        call rsnap(filename, gf(1+ipx*nx:mx+ipx*nx,1+ipy*ny:my+ipy*ny,:,1:mvar_in), mvar_in)
        t_sp = t
!
        if (lroot) then
          t_test = t_sp
        elseif (t_test /= t_sp) then
          write (*,*) 'ERROR: '//trim(directory_snap)//'/'//trim(filename)//' IS INCONSISTENT: t=', t_sp
          stop 1
        endif
!
! collect x coordinates:
!
        gx       (1+ipx*nx:mx+ipx*nx) = x
        gdx_1    (1+ipx*nx:mx+ipx*nx) = dx_1
        gdx_tilde(1+ipx*nx:mx+ipx*nx) = dx_tilde
!
! collect y coordinates:
!
        gy       (1+ipy*ny:my+ipy*ny) = y
        gdy_1    (1+ipy*ny:my+ipy*ny) = dy_1
        gdy_tilde(1+ipy*ny:my+ipy*ny) = dy_tilde
!
      enddo
    enddo
!
! collect z coordinates:
!
    gz       (1+ipz*nz:mz+ipz*nz) = z
    gdz_1    (1+ipz*nz:mz+ipz*nz) = dz_1
    gdz_tilde(1+ipz*nz:mz+ipz*nz) = dz_tilde
!
! write xy-layer:
!
    do pa = 1, mvar_io

      if (lfirst_proc_z) then
        start_pos = 1
      else
        start_pos = nghost + 1
      endif

      if (llast_proc_z ) then
        end_pos = mz
      else
        end_pos = nghost + nz
      endif

      do pz = start_pos, end_pos
        write(lun_output, rec=pz+ipz*nz+(pa-1)*mzgrid) gf(:,:,pz,pa)
      enddo

    enddo
!
  enddo
!
! write additional data:
!
  close(lun_output)
  open (lun_output, FILE=trim(directory_out)//'/'//filename, FORM='unformatted', position='append', status='old')
  write(lun_output) t, gx, gy, gz, dx, dy, dz
  close(lun_output)
!
  if (IO_strategy == 'dist') then
!
! write global grid:
!
    open (lun_output, FILE=trim(directory_out)//'/grid.dat', FORM='unformatted', status='replace')
    write(lun_output) t_sp, gx, gy, gz, dx, dy, dz
    write(lun_output) dx, dy, dz
    write(lun_output) Lxyz
    write(lun_output) gdx_1, gdy_1, gdz_1
    write(lun_output) gdx_tilde, gdy_tilde, gdz_tilde
    close(lun_output)
  endif
!
  print *, 'Writing snapshot for time t =', t
!
endprogram pc_collect
