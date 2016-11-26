! This is a tool to collect distributed data into one HDF5 file.
!
! $Id: pc_h5collect.f90 23201 2015-04-06 02:11:37Z st.tuomisto@gmail.com $
!***********************************************************************


program pc_h5collect
!
  use Cdata
  use Cparam
  !use Diagnostics
  !use Filter
  use Grid, only: initialize_grid
  !use IO
  !use Messages
  !use Mpicomm
  !use Param_IO
  use Register
  !use Snapshot
  !use Sub
  use General, only: safe_character_assign, itoa, directory_names_std
  use Syscalls, only: sizeof_real
  use Analyzers
  use HDF5
  use H5DS
  !use mpi
!
  implicit none
!
  include 'mpif.h'
  external sizeof_c

  ! Input and output variables & parameters

  character (len=fnlen) :: datafile, datafilename,infile, datagroup, logfile, &
                            hdffile, field, outfile
  character (len=*), parameter :: directory_out = 'data/allprocs'
  character (len=*), parameter :: mfcollect_in_file = 'mfcollect.in'

  character (len=fnlen), dimension(:)  , allocatable :: datagroups, fields, &
                                                     analyzernames
  integer :: hdferr
  integer(HID_T) :: hdfmemtype, hdffiletype, hdf_fileid, hdf_datagroup, hdf_filespace
  integer(HID_T), dimension(:), allocatable :: hdf_fieldgroups
  integer(HID_T), dimension(:,:), allocatable :: hdf_spaces, hdf_chunkspaces, hdf_datasets, &
                                                 hdf_plist_ids

  integer(HSIZE_T) , dimension(:,:), allocatable :: hdf_dims_full, hdf_dims, hdf_maxdims
  integer(HSIZE_T) :: hdf_chunkdims(4), chunksize
  
  integer(HID_T) :: hdf_file_plist
                    
  integer(HSIZE_T) , dimension(4) :: hdf_stride, hdf_count, hdf_block, hdf_offsets

  
  
  logical :: hdf_exists
  
  
  ! Internal parameters
  integer, parameter            :: len_datagroups = 30, len_fields = 100, &
                                   len_analyzers = 30, &
                                   rkind = selected_real_kind(10,40)

  ! Job communication parameters and variables
  integer, parameter            :: command_shutdown   = 0, &
                                   command_analyze    = 1, &
                                   command_hdfopen    = 2, &
                                   command_hdfclose   = 3 
             
  integer  :: nprocs_needed, analysisproc, command, fake_proc, nprocs, &
              nmpiprocs, mpierr
  integer :: mpistat(MPI_STATUS_SIZE)

  logical :: initerror, runerror, received, filecheck
  
  integer, dimension(:), allocatable   :: fake_procs
  
  character (len=intlen) :: chproc
  
  ! Timing variables
  real(kind=rkind)     :: t_taken_full, t_taken_analysis, analysisstart

  ! Data read variables
  integer :: idatagroup, ifield, ierr, ntimesteps, nfields, &
             i, j, nanalyzers, &
             ianalyzer, maxtimesteps
             
  integer :: dims(3), dims_full(3)
  
  integer, dimension(:), allocatable :: avgdims
  integer, dimension(:,:), allocatable :: offsets

  real, dimension(:,:,:,:,:), allocatable   :: dataarray_full
  real, dimension(:,:,:,:)  , allocatable   :: dataarray
  real, dimension(:)        , allocatable   :: t_values

  integer            :: datalen, tlen, data_stride, &
                            data_start, tsteplen, pos, &
                            averagelen, resultlen, &
                            filesize_check
  integer            :: filesize
  
  integer , parameter   :: t_start = 5

  ! Analyzer parameters
  procedure(AnalyzerTemplate), pointer :: analyzerfunction
  
  interface ljustify
    function ljustifyI(i)
        integer*4, intent(in) :: i
        character(len=30)  :: ljustifyI
    end function ljustifyI
    function ljustifyL(i)
        logical, intent(in) :: i
        character(len=30)  :: ljustifyL
    end function ljustifyL
  end interface

  
  namelist /collect_config/ datagroups, ldebug
  namelist /xaver_config/ analyzernames, maxtimesteps
  namelist /zaver_config/ analyzernames, maxtimesteps


! Initialize MPI

  !call mpicomm_init()
  !call initialize_messages()
  !call initialize_mpicomm()
  call MPI_INIT(mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nmpiprocs, mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, iproc , mpierr)
  analysisstart = MPI_WTIME() 

  nullify(analyzerfunction)

  datalen   = -1
  maxtimesteps = huge(maxtimesteps)
  
  initerror = .false.
  runerror  = .false.
  received  = .false.
  ldebug    = .false.

  lroot = (iproc==root)

!
!  Identify version.
!
!  Read parameters from start.x (default values; may be overwritten by
!  read_runpars).

!

  !call read_all_init_pars()
!
!  Read parameters and output parameter list.
!
  !call read_all_run_pars()
!
!  Register phierrysics modules.
!
  call register_modules()

  call H5open_F(hdferr)

! Root part of the analysis, read config and co-ordinate other programs

  if (lroot) then

    write (*,'(a,i1,a)') ' This is a ', dimensionality, '-D run'
    print *, 'nxgrid, nygrid, nzgrid=', nxgrid, nygrid, nzgrid
    print *, 'Lx, Ly, Lz=', Lxyz
    print *, '      Vbox=', Lxyz(1)*Lxyz(2)*Lxyz(3)
    
    ! Collect information on 

    allocate(datagroups(len_datagroups))
    allocate(fields(len_fields))
    allocate(analyzernames(len_analyzers))
    allocate(fake_procs(nprocx*nprocy*nprocz))
    allocate(offsets(4,nprocx*nprocy*nprocz))
    datagroups       = ''
    fields          = ''
    analyzernames   = ''
    fake_procs     = 0

    open(UNIT=1, FILE=mfcollect_in_file, IOSTAT=ierr)
    if (ierr /= 0) then
      write(*,*) 'Error opening '//mfcollect_in_file
      initerror = .true.
    else
      read(1, NML=collect_config, IOSTAT=ierr)
      if (ierr /= 0) then
        write(*,*) 'Error reading main configuration'
        initerror = .true.
      end if
    end if

    call MPI_BCAST(initerror, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
    if (.not. initerror) then
      analysisproc = 1
      do idatagroup=1,len_datagroups

        datagroup = trim(datagroups(idatagroup))
        ! Field dependent part
        ! ----------------------------------------
        nprocs_needed = 0
        if (datagroup == 'xaver') then
          ! This needs to be written based on zaver example 
        else if (datagroup == 'zaver') then
          dims        = [nx, ny, 1]
          dims_full   = [nxgrid, nygrid, 1]
          infile      ='zaver.in'
          outfile = trim(datadir)//'/zaverages.h5'
          datafilename='zaverages.dat'

          ! Read configuration
          read(1, NML=zaver_config, IOSTAT=ierr)
          if (ierr /= 0) then
            write(*,*) 'Error reading zaver configuration'
            initerror = .true.
          end if
          
          ! Determine needed processors
          nprocs_needed = nprocx*nprocy
          j=1
          do fake_proc=0,nprocx*nprocy*nprocz
            if (lprocz_slowest) then
              ipx = modulo(fake_proc, nprocx)
              ipy = modulo(fake_proc/nprocx, nprocy)
              ipz = fake_proc/nprocxy
            else
              ipx = modulo(fake_proc, nprocx)
              ipy = fake_proc/nprocxy
              ipz = modulo(fake_proc/nprocx, nprocy)
            endif
            if (ipz == 0) then
              fake_procs(j) = fake_proc
              offsets(:,j) = [dims(1)*ipx, dims(2)*ipy,0,0]
              j = j + 1
            end if
          end do
        else
          exit
        end if
         
        ! Universal part
        ! ----------------------------------------
       
        ! Determine picked fields
        nfields = 0
        open(UNIT=2, FILE=infile, STATUS='old', IOSTAT=ierr)
        if (ierr /= 0) then
          write(*,*) 'Could not open '//infile
          initerror = .true.
          exit
        end if
        do ifield=1,len_fields
          read(2,*,iostat=ierr) field
          if (ierr == 0) then
            fields(ifield) = trim(field)
            nfields = nfields + 1
          else
            exit
          end if
        end do
        close(2)
        if (nfields == 0) then
          write(*,*) 'No fields found.'
          initerror=.true.
        end if

        ! Determine analyzers
        nanalyzers = 0
        do i=1,len_analyzers
          if (analyzernames(i) /= '') then
            do j=1,npossibleanalyzers
              if (analyzernames(i) == possibleanalyzers(j)) then
                nanalyzers = nanalyzers + 1
                exit
              end if
            end do
            if (nanalyzers < i) then
              write(*,*) 'Invalid analyzer found: '//analyzernames(i)
              initerror = .true.
            end if
          else
            exit
          end if
        end do
 
        ! Check datafile properties

        chproc = itoa(fake_procs(1))
        call directory_names_std()
        !call safe_character_assign (directory, &
        !     trim(datadir)//'/proc'//chproc)
        call safe_character_assign (datafile, &
             trim(directory)//'/'//trim(datafilename))
        open(3,file=datafile, access='stream',form='unformatted', action='read', status='old', iostat=ierr)
        if (ierr /= 0) then
          write(*,*) 'Error opening datafile: '//datafile
          initerror = .true.
        end if

        inquire(3, size=filesize)
        pos=1
        read(3, pos=pos, IOSTAT=ierr) tlen
        if (ierr /= 0) then
          write(*,*) 'Error reading tlen'
          initerror = .true.
        end if
        pos=9+tlen
        read(3, pos=pos, IOSTAT=ierr) datalen
        if (ierr /= 0) then
          write(*,*) 'Error reading datalen'
          initerror = .true.
        end if
        close(3)

        data_stride = datalen/(nfields*dims(1)*dims(2)*dims(3))

        ntimesteps  = int(floor(real(filesize)/(16+tlen+datalen)))
        write(*,'(a30,a30)') &
        ' Number of timesteps:          ', ljustify(ntimesteps)
        if (ntimesteps > maxtimesteps) then
               ntimesteps  = maxtimesteps
        end if
        data_start  = t_start + tlen + 8
        tsteplen    = 16 + tlen + datalen
        

        !write(*,*) 'root part:', iproc,dims(2),dims_full(2), dims(2)read, dims(2)runs, nprocy, averagelen, data_stride

        write(*,'(a30)')     ' Analysis information          '
        write(*,'(a30)')     ' ----------------------------- '
        write(*,'(a30,a)') &
        ' Average type:                 ', trim(datagroup)
        write(*,'(a30,a30)') &
        ' Dimension 1 size:             ', ljustify(dims(1))
        write(*,'(a30,a30)') &
        ' Dimension 2 size:             ', ljustify(dims(2))
        write(*,'(a30,a30)') &
        ' Dimension 3 size:             ', ljustify(dims(3))
        write(*,'(a30,a)') &
        ' Datafile:                     ', trim(datafile)
        write(*,'(a30,a30)') &
        ' Number of fields:             ', ljustify(nfields)
        write(*,'(a30,a30)') &
        ' File size (bytes):            ', ljustify(filesize)
        write(*,'(a30,a30)') &
        ' Number of timesteps:          ', ljustify(ntimesteps)
        write(*,'(a30,a30)') &
        ' Datapoints per timestep:      ', ljustify(datalen/data_stride)
        write(*,'(a30,a30)') &
        ' Single precision time values: ', ljustify(tlen==4)
        write(*,'(a30,a30)') &
        ' Single precision data values: ', ljustify(data_stride==4)

        !c

        write(*,*) 'Starting analysis...'
        analysisproc=1
        if ((nanalyzers>0) .and. (nfields > 0) .and. (nprocs_needed > 0) &
          .and. (.not. initerror)) then
          ! Send datasets to create
          command = command_hdfopen
          do i=1,nmpiprocs-1
            call MPI_SEND(command, 1, &
                          MPI_INTEGER, analysisproc, 0, MPI_COMM_WORLD, mpierr)
            analysisproc = analysisproc + 1
            if (analysisproc>=nmpiprocs) then
              analysisproc = 1
            end if
          end do
          call BCastInfo()
          call OpenH5File()

          ! Commence the analysis
          command = command_analyze
          do i=1, nprocs_needed
            fake_proc = fake_procs(i)
            call MPI_SEND(command, 1, &
                          MPI_INTEGER, analysisproc, 0, MPI_COMM_WORLD, mpierr)
            call MPI_SEND(fake_proc, 1, &
                         MPI_INTEGER, analysisproc, 1, MPI_COMM_WORLD, mpierr)
            call MPI_SEND(offsets(:,i), 4, &
                         MPI_INTEGER, analysisproc, 2, MPI_COMM_WORLD, mpierr)
            call MPI_RECV(received, 1, &
                          MPI_LOGICAL, analysisproc, 3, MPI_COMM_WORLD, mpistat, mpierr)
            analysisproc = analysisproc+1
            if (analysisproc>=nmpiprocs) then
              analysisproc = 1
            end if
          end do
          command = command_hdfclose
          do i=1,nmpiprocs-1
            call MPI_SEND(command, 1, MPI_INTEGER, analysisproc, 0, MPI_COMM_WORLD, mpierr)
            analysisproc = analysisproc + 1
            if (analysisproc>=nmpiprocs) then
              analysisproc = 1
            end if
          end do
          call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
          call CloseH5File()
        end if

        if (allocated(avgdims)) then
            deallocate(avgdims)
        end if
      end do
      ! Send shutdown to other fields
      command = command_shutdown
      do i=1,nmpiprocs-1
        call MPI_SEND(command, 1, MPI_INTEGER, analysisproc, 0, MPI_COMM_WORLD, mpierr)
        analysisproc = analysisproc + 1
        if (analysisproc>=nmpiprocs) then
          analysisproc = 1
        end if
      end do
    end if
    close(1)
    if (allocated(datagroups)) then
      deallocate(datagroups)
    end if
    if (allocated(fields)) then
      deallocate(fields)
    end if
    if (allocated(analyzernames)) then
      deallocate(analyzernames)
    end if
    if (allocated(fake_procs)) then
      deallocate(fake_procs)
    end if
  else
    ! Non-root part of the analysis. Do the analysis.
    call MPI_BCAST(initerror, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
    if (.not. initerror) then
      command = -1
      do
        call MPI_RECV(command, 1, &
          MPI_INTEGER, root, 0, MPI_COMM_WORLD, mpistat, mpierr)
        if (command == command_hdfopen) then
          call BCastInfo()
          call OpenH5File()
        else if (command == command_hdfclose) then
          call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
          call CloseH5File()
        else if (command == command_analyze) then
          ! MPI communication part
          call MPI_RECV(fake_proc, 1, &
                        MPI_INTEGER, root, 1, MPI_COMM_WORLD, mpistat, mpierr)
          if (.not. allocated(offsets)) then
            allocate(offsets(4,1))
          end if
          call MPI_RECV(offsets(:,1), 4, &
                        MPI_INTEGER, root, 2, MPI_COMM_WORLD, mpistat, mpierr)
          call MPI_SEND(received, 1, &
                        MPI_LOGICAL, root, 3, MPI_COMM_WORLD, mpierr)
          chproc = itoa (fake_proc)
          call safe_character_assign(directory, trim(datadir)//'/proc'//chproc)
          call safe_character_assign (datafile, &
               trim(directory)//'/'//trim(datafilename))

          ! Data reading part
          ! -------------------------------------------------------------
      
          call Analyze()
        else
          exit
        end if
      end do
    end if
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
  call H5close_F(hdferr)
  call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
  write(*,*) 'Worker ', iproc, ' finished.'
  write(*,*) 'Time taken: ', MPI_WTIME() - analysisstart
  call MPI_FINALIZE(mpierr)

  contains

    subroutine Analyze

      runerror = .false.

      ! Check datafile properties

      open(3,file=datafile, access='stream',form='unformatted', action='read', status='old', iostat=ierr)
      if (ierr /= 0) then
        write(*,*) 'Error opening datafile: '//datafile
        runerror = .true.
      end if
      
      inquire(3, size=filesize_check)
      if (filesize_check /= filesize) then
        write(*,*) 'File sizes differ between procs: '//chproc
        runerror = .true.
      end if
      close(3)


      if (ldebug) then
        ! Open logfile
        logfile = trim(directory)//'/'//trim(field)//'_h5collect.log'
        open(2, file=logfile, action='write', status='replace', iostat=ierr)
        ! Log analysis parameters

        write(2,'(a30)')     ' Analysis information          '
        write(2,'(a30)')     ' ----------------------------- '
        write(2,'(a30,a30)') &
        ' Worker:                       ', ljustify(fake_proc)
        write(2,'(a30,a30)') &
        ' Dimension 1 size:             ', ljustify(dims(1))
        write(2,'(a30,a30)') &
        ' Dimension 2 size:             ', ljustify(dims(2))
        write(2,'(a30,a)') &
        ' Datafile:                     ', trim(datafile)
        write(2,'(a30,a30)') &
        ' Number of fields:             ', ljustify(nfields)
        write(2,'(a30,i30)') &
        ' File size (bytes):            ', filesize
        write(2,'(a30,a30)') &
        ' Number of timesteps:          ', ljustify(ntimesteps)
        write(2,'(a30,a30)') &
        ' Datapoints per timestep:      ', ljustify(datalen/data_stride)
        write(2,'(a30,a30)') &
        ' Single precision time values: ', ljustify(tlen==4)
        write(2,'(a30,a30)') &
        ' Single precision data values: ', ljustify(data_stride==4)
        close(2)
      end if

      if (.not. allocated(t_values)) then
        allocate(t_values(ntimesteps),stat=ierr) 
        if (ierr /= 0) then
          write(*,*) 'Proc '//chproc//': Error allocating t_values'
          runerror = .true.
        end if
      end if
      if (.not. allocated(dataarray_full)) then
        allocate(dataarray_full(dims(1),dims(2),dims(3),nfields, ntimesteps),stat=ierr)
        if (ierr /= 0) then
          write(*,*) 'Proc '//chproc//': Error allocating dataarray_full'
          runerror = .true.
        end if
      end if
      if (.not. allocated(dataarray)) then
        allocate(dataarray(dims(1),dims(2),dims(3),nfields),stat=ierr)
        if (ierr /= 0) then
          write(*,*) 'Proc '//chproc//': Error allocating dataarray'
          runerror = .true.
        end if
      end if
 
      ! Do data loading
      ! -------------------------

      t_values    = 0

      t_taken_full = MPI_WTIME()
      open(3,file=datafile, form='unformatted', &
          action='read', status='old', IOSTAT=ierr)
      do i=1,ntimesteps
        read(3, IOSTAT=ierr) t_values(i)
        read(3, IOSTAT=ierr) dataarray_full(1:dims(1), 1:dims(2), 1:dims(3), 1:nfields, i)
      end do
      close(3)

      do ianalyzer=1,nanalyzers
        nullify(analyzerfunction)
        resultlen = ntimesteps
        call getAnalyzer(analyzernames(ianalyzer),analyzerfunction,resultlen)
        if (allocated(dataarray) .and. ((resultlen /= size(dataarray,dim=4)))) then
          deallocate(dataarray)
        end if
        if (.not. allocated(dataarray)) then
          allocate(dataarray(dims(1),dims(2),dims(3),resultlen),stat=ierr)
          if (ierr /= 0) then
            write(*,*) 'Proc '//chproc//': Error allocating dataarray'
            runerror = .true.
            exit
          end if
        end if

        do ifield=1,nfields
          dataarray(1:dims(1),1:dims(2),1:dims(3),1:resultlen) = analyzerfunction(dataarray_full(:,:,:,ifield,:), &
                                 dims(1), dims(2), dims(3), ntimesteps, resultlen)
          
          ! Parallel data writing to HDF file
          call H5Dget_space_F(hdf_datasets(ianalyzer,ifield), hdf_filespace, hdferr)
          hdf_offsets  = [ offsets(1,1), offsets(2,1), offsets(3,1), offsets(4,1) ]
          hdf_stride   = [ 1, 1, 1, 1 ]
          hdf_block    = [ dims(1), dims(2), dims(3), 1 ]
          hdf_count    = [ 1, 1, 1, resultlen ]

          call H5Sselect_hyperslab_F(hdf_filespace, H5S_SELECT_SET_F, hdf_offsets, hdf_count, &
                                     hdferr, hdf_stride, hdf_block)
          !call H5Sselect_hyperslab_F(hdf_filespace, H5S_SELECT_SET_F, hdf_offsets, hdf_count, &
          !                          hdferr)
          call H5Dwrite_F(hdf_datasets(ianalyzer,ifield), hdfmemtype, dataarray, &
                                     hdf_dims(:,ianalyzer), hdferr, &
                                     file_space_id  = hdf_filespace, &
                                     mem_space_id   = hdf_chunkspaces(ianalyzer, ifield))

          call H5Sclose_F(hdf_filespace, hdferr)


        end do

      end do

      t_taken_full = MPI_WTIME() - t_taken_full
      write(*,*) 'Worker '//trim(chproc)//' finished analyzing '//trim(datagroup)//'. Time taken: ' , t_taken_full

      !call H5SCLOSE_F(hdfspace, hdferr)
      
      if (allocated(dataarray)) then
        deallocate(dataarray)
      end if
      if (allocated(t_values)) then
        deallocate(t_values)
      end if
      if (allocated(dataarray_full)) then
        deallocate(dataarray_full)
      end if
      nullify(analyzerfunction)

    end subroutine Analyze

    subroutine BCastInfo
      !implicit none
      call MPI_BCAST(datagroup, fnlen, &
                MPI_CHARACTER, root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(datafilename, fnlen, &
                MPI_CHARACTER, root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(outfile, fnlen, &
                MPI_CHARACTER, root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(nfields, 1, &
                MPI_INTEGER,  root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(nanalyzers, 1, &
                MPI_INTEGER,  root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(dims, 3, &
                MPI_INTEGER,  root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(dims_full, 3, &
                MPI_INTEGER,  root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(datalen, 1, &
                MPI_INTEGER, root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(filesize, 1, &
                MPI_INTEGER, root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(tlen, 1, &
                MPI_INTEGER, root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(ntimesteps, 1, &
                MPI_INTEGER, root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(data_stride, 1, &
                MPI_INTEGER, root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(data_start, 1, &
                MPI_INTEGER, root, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(tsteplen, 1, &
                MPI_INTEGER, root, MPI_COMM_WORLD, mpierr)
      if (.not. lroot) then
        if (allocated(fields)) then
          deallocate(fields)
        end if
        allocate(fields(nfields))
        if (allocated(analyzernames)) then
          deallocate(analyzernames)
        end if
        allocate(analyzernames(nanalyzers))
      end if
      do ifield=1,nfields
        call MPI_BCAST(fields(ifield), fnlen, &
                  MPI_CHARACTER, root, MPI_COMM_WORLD, mpierr)
      end do
      do ianalyzer=1,nanalyzers
        call MPI_BCAST(analyzernames(ianalyzer), fnlen, &
                  MPI_CHARACTER, root, MPI_COMM_WORLD, mpierr)
      end do
      filecheck = .true.
      call MPI_REDUCE(filecheck, received, 1, &
                      MPI_LOGICAL, MPI_LAND, root, MPI_COMM_WORLD, mpistat, mpierr)
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    end subroutine BCastInfo

    subroutine OpenH5File()
 
      call H5Pcreate_F(H5P_FILE_ACCESS_F, hdf_file_plist, hdferr)
      call H5Pset_fapl_mpio_F(hdf_file_plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
      ! Open hdf5-file
      inquire(file=outfile, exist=hdf_exists)
      if (hdf_exists) then
        open(3, file=outfile, action='write', status='replace', iostat=ierr)
        close(3, status='delete')
      end if
      call H5Fcreate_F(outfile, H5F_ACC_TRUNC_F, hdf_fileid, hdferr, access_prp = hdf_file_plist)
      call H5Pclose_F(hdf_file_plist,hdferr)
       
      if (data_stride == 4) then
        hdfmemtype  = H5T_NATIVE_REAL
        hdffiletype = H5T_IEEE_F32LE
      else if (data_stride == 8) then
        hdfmemtype  = H5T_NATIVE_DOUBLE
        hdffiletype = H5T_IEEE_F64LE
      else
        write(*,*) 'Could not determine the data type.'
        initerror = .true.
      end if

      if (data_stride /= sizeof_real()) then
        write(*,*) 'Data real precision has different precision than the collector.'
        initerror = .true.
      end if

      ! Create dataset
      if (.not. initerror) then
        if (allocated(hdf_fieldgroups)) then
          deallocate(hdf_fieldgroups)
        end if
        if (allocated(hdf_spaces)) then
          deallocate(hdf_spaces)
        end if
        if (allocated(hdf_chunkspaces)) then
          deallocate(hdf_chunkspaces)
        end if
        if (allocated(hdf_datasets)) then
          deallocate(hdf_datasets)
        end if
        if (allocated(hdf_plist_ids)) then
          deallocate(hdf_plist_ids)
        end if
        if (allocated(hdf_dims_full)) then
          deallocate(hdf_dims_full)
        end if
        if (allocated(hdf_dims)) then
          deallocate(hdf_dims)
        end if
        !if (allocated(hdf_chunkdims)) then
        !  deallocate(hdf_chunkdims)
        !end if
        allocate(hdf_fieldgroups(nfields))
        allocate(hdf_spaces(nanalyzers,nfields))
        allocate(hdf_datasets(nanalyzers,nfields))
        allocate(hdf_chunkspaces(nanalyzers,nfields))
        allocate(hdf_plist_ids(nanalyzers,nfields))
        allocate(hdf_dims_full(4,nfields))
        allocate(hdf_dims(4,nfields))
        allocate(hdf_maxdims(4,nfields))
        !allocate(hdf_chunkdims(4,nfields))
        do ianalyzer=1,nanalyzers
          resultlen = ntimesteps
          call getAnalyzer(analyzernames(ianalyzer),analyzerfunction,resultlen)
          hdf_dims_full(:,ianalyzer) = [dims_full(1), dims_full(2), dims_full(3), resultlen]
          hdf_dims(:,ianalyzer)      = [dims(1), dims(2), dims(3), resultlen]
          hdf_maxdims(:,ianalyzer) = [dims_full(1), dims_full(2), dims_full(3), -1]
        end do
        chunksize = data_stride
        do i=1,3
          if (chunksize*dims(i) < 10*1024*1024) then
            hdf_chunkdims(i) = dims(i)
            chunksize = chunksize*dims(i)
          end if
        end do
        if (chunksize < 10*1024*1024) then
          hdf_chunkdims(4) = 10
        else
          hdf_chunkdims(4) = 1
        end if

        call H5Gcreate_F(hdf_fileid, "/"//trim(datagroup), hdf_datagroup, hdferr)
        if (lroot .and. ldebug) then
          write(*,*) 'CreatedG: ', "/"//trim(datagroup)
        end if
        do ifield=1,nfields
          call H5Gcreate_F(hdf_datagroup, "/"//trim(datagroup)//"/"//fields(ifield), &
                         hdf_fieldgroups(ifield), hdferr)
          if (lroot .and. ldebug) then
            write(*,*) 'CreatedG: ', "/"//trim(datagroup)//"/"//fields(ifield)
          end if
          do ianalyzer=1,nanalyzers
            call H5Screate_simple_F(4, hdf_dims_full(:,ianalyzer), &
                                    hdf_spaces(ianalyzer,ifield), hdferr, hdf_maxdims)
            call H5Screate_simple_F(4, hdf_dims(:,ianalyzer), &
                                    hdf_chunkspaces(ianalyzer,ifield), hdferr, hdf_maxdims)
            if (lroot .and. ldebug) then
              write(*,*) 'CreatedS: ', trim(analyzernames(ianalyzer))
            end if
            call H5Pcreate_F(H5P_DATASET_CREATE_F, hdf_plist_ids(ianalyzer,ifield), hdferr)
            call H5Pset_chunk_F(hdf_plist_ids(ianalyzer, ifield), 4, hdf_chunkdims, hdferr)
            !call H5Pset_deflate_F(hdf_plist_ids(ianalyzer,ifield), 6, hdferr)
            !call H5Pset_chunk_F(hdf_plist_ids(ianalyzer, ifield), 4, hdf_dims(:,ianalyzer), hdferr)
            call H5Dcreate_F(hdf_fieldgroups(ifield), trim(analyzernames(ianalyzer)), &
                             hdffiletype, hdf_spaces(ianalyzer,ifield), &
                             hdf_datasets(ianalyzer,ifield), hdferr, hdf_plist_ids(ianalyzer,ifield))
            call H5DSset_label_F(hdf_datasets(ianalyzer,ifield),4,'x', hdferr)
            call H5DSset_label_F(hdf_datasets(ianalyzer,ifield),3,'y', hdferr)
            call H5DSset_label_F(hdf_datasets(ianalyzer,ifield),2,'z', hdferr)
            call H5DSset_label_F(hdf_datasets(ianalyzer,ifield),1,'t', hdferr)
            
            call H5Sclose_F(hdf_spaces(ianalyzer,ifield), hdferr)
            if (lroot .and. ldebug) then
              write(*,*) 'CreatedD: ', trim(analyzernames(ianalyzer))
            end if
          end do
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
     
      end if

    end subroutine OpenH5File

    subroutine CloseH5File

      do ifield=1,nfields
        do ianalyzer=1,nanalyzers
          call H5Dclose_F(hdf_datasets(ianalyzer,ifield),hdferr)
          call H5Pclose_F(hdf_plist_ids(ianalyzer,ifield), hdferr)
          if (lroot .and. ldebug) then
            write(*,*) 'ClosedD: ', trim(analyzernames(ianalyzer))
          end if
          call H5Sclose_F(hdf_chunkspaces(ianalyzer,ifield), hdferr)
          if (lroot .and. ldebug) then
            write(*,*) 'ClosedS: ', trim(analyzernames(ianalyzer))
          end if
        end do
        call H5Gclose_F(hdf_fieldgroups(ifield), hdferr)
        if (lroot .and. ldebug) then
          write(*,*) 'ClosedG: ', "/"//trim(datagroup)//"/"//fields(ifield)
        end if
      end do
      call H5Gclose_F(hdf_datagroup, hdferr)
      if (lroot .and. ldebug) then
        write(*,*) 'ClosedG: ', "/"//trim(datagroup)
      end if
      
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
      call H5Fclose_F(hdf_fileid, hdferr)

    end subroutine CloseH5File

end program pc_h5collect

function ljustifyI(i)
    integer, intent(in) :: i
    character(len=30)  :: ljustifyI
    write(ljustifyI, '(i30)') i
    ljustifyI = adjustl(ljustifyI)
    return
end function ljustifyI

function ljustifyL(l)
    logical, intent(in) :: l
    character(len=30)  :: ljustifyL
    write(ljustifyL, '(l30)') l
    ljustifyL = adjustl(ljustifyL)
    return
end function ljustifyL
