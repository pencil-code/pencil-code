! This is a tool to collect distributed zaverages into one file.
!
! $Id: pc_meanfield_collect.f90 23199 2015-03-23 02:11:37Z st.tuomisto@gmail.com $
!***********************************************************************
program pc_meanfield_collect
!
  use Cdata
  use Cparam, only: fnlen
  use Diagnostics
  use Filter
  use Grid, only: initialize_grid
  use IO
  use Messages
  use Mpicomm
  use Param_IO
  use Register
  use Snapshot
  use Sub
  use General, only: backskip_to_time
  use Syscalls, only: sizeof_real
!
  implicit none
!
  character (len=fnlen) :: arg_meanfield, arg_analyzer, arg_paramater, datafile, inputfile
  character (len=*), parameter :: directory_out = 'data/allprocs'
!
  logical :: ex
  integer :: mvar_in, io_len, pz, pa, start_pos, end_pos, alloc_err, len_in_rec, ierr, &
             nargs, i, ndatapoints, naverages, ndim1, ndim2
  integer(kind=8) :: rec_len
  
  real, dimension(:,:,:,:)  , allocatable   :: averages, results
  real, dimension(:)        , allocatable   :: t_values
  real :: t_sp, t_test , t_datapoint ! t in single precision for backwards compatibility
  logical :: process_constraint

  inquire (IOLENGTH=io_len) t_sp


!
!  Identify version.
!
  if (lroot) call svn_id( &
      '$Id: pc_meanfield_collect.f90 23199 2015-03-23 02:11:37Z st.tuomisto@gmail.com $')
!
!  Initialize the message subsystem, eg. color setting etc.a

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
! Calculate dimensionality
!
  dimensionality = min(nxgrid-1,1) + min(nygrid-1,1) + min(nzgrid-1,1)
!
!  Register phierrysics modules.
!
  call register_modules()

!  Print resolution and dimension of the simulation.
!
  if (lroot) write (*,'(a,i1,a)') ' This is a ', dimensionality, '-D run'
  if (lroot) print *, 'nxgrid, nygrid, nzgrid=', nxgrid, nygrid, nzgrid
  if (lroot) print *, 'Lx, Ly, Lz=', Lxyz
  if (lroot) print *, '      Vbox=', Lxyz(1)*Lxyz(2)*Lxyz(3)
!
  
  call directory_names() 
  
  ! Process arguments and set analysis variables

  nargs = command_argument_count()
  if (nargs == 0) call fatal_error ('pc_meanfield_collect', 'No command line arguments supplied', .true.)
  
  call get_command_argument(1, arg_meanfield)
  call get_command_argument(2, arg_analyzer)

  if (arg_meanfield == 'zaver') then 
    datafile=trim(directory_dist)//'/zaverages.dat'
    inputfile='zaver.in'
    naverages = parallel_count_lines(inputfile)
    if (naverages == 0) call fatal_error ('pc_meanfield_collect', 'No z-averages present in zaver.in', .true.)
    process_constraint = (lfirst_proc_z.and.naverages>0) 
    ndim1 = nx
    ndim2 = ny
  end if

  ! Check data with processes that match the process_constraint


  if (process_constraint) then
    allocate(averages(ndim1,ndim2,naverages,1),stat=ierr)
    if (ierr > 0) then
      call fatal_error ('pc_meanfield_collect', 'Error allocating averages '//datafile,.true.)
    end if
    open(1,file=datafile, form='unformatted', action='read', status='old', IOSTAT=ierr)
    t_datapoint = 0
    averages = 0.0
    ndatapoints = 0

    ! Calculate the amount of data within average file
    datacheck_loop : do
      
      ! Read time values

      read(1, IOSTAT=ierr) t_datapoint
      if (ierr < 0) then
        ! End of file reached
        exit datacheck_loop
      else if (ierr /= 0) then
        call fatal_error ('pc_meanfield_collect', 'Error reading the file '//datafile,.true.)
      end if

      ! Read averages
      
      read(1,  IOSTAT=ierr) averages(:,:,:,1)
      ndatapoints = ndatapoints + 1
      if (ierr /= 0) then
        call fatal_error ('pc_meanfield_collect', 'Error reading the file '//datafile,.true.)
      end if
    end do datacheck_loop
    close(1)

    ! Deallocate averages in order to allocate them fully
    deallocate(averages)

    ! ---------------------------------------------------------
  
    ! Get data with processes that match the process_constraint
  
    if (lroot) print *, 'ndatapoints=', ndatapoints
    allocate(averages(ndim1,ndim2,naverages,ndatapoints),stat=ierr)
    if (ierr > 0) then
      call fatal_error ('pc_meanfield_collect', 'Error allocating averages '//datafile,.true.)
    end if
    allocate(t_values(ndatapoints),stat=ierr)
    if (ierr > 0) then
      call fatal_error ('pc_meanfield_collect', 'Error allocating t_values '//datafile,.true.)
    end if
    open(1,file=datafile, form='unformatted', action='read', status='old', IOSTAT=ierr)
    t_values = 0
    averages = 0.0

    ! Read data in

    dataread_loop : do i=1,ndatapoints
      
      ! Read time values

      read(1, IOSTAT=ierr) t_values(i)

      ! Read averages
      
      read(1,  IOSTAT=ierr) averages(:,:,:,i)

    end do dataread_loop

    if (lroot .and. ip<=10) then
      write (*,'(a,1000f10.5)') ' Datapoint t_values='//new_line('A'),t_values
    end if
    close(1)

    ! Do analysis for data

    if (avg_analyzer) then

    end if


    ! Deallocate averages
    deallocate(averages)
  end if
!
!  Give all modules the possibility to exit properly.
!
!  Free any allocated memory.
!
  call fnames_clean_up()
  call vnames_clean_up()
!
endprogram pc_meanfield_collect

