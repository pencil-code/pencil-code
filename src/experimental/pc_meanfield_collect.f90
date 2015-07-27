! This is a tool to collect distributed zaverages into one file.
!
! $Id: pc_meanfield_collect.f90 23201 2015-04-06 02:11:37Z st.tuomisto@gmail.com $
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
  use General, only: safe_character_assign, itoa
  use Syscalls, only: sizeof_real
  use Analyzers
  use HDF5
!
  implicit none
!
  include 'mpif.h'

  ! Input and output variables & parameters

  character (len=fnlen) :: datafile, datafilename,infile, avgname, logfile, &
                            outputfile, hdffile, avgfield, phdffile
  character (len=*), parameter :: directory_out = 'data/allprocs'
  character (len=*), parameter :: mfcollect_in_file = 'mfcollect.in'

  character (len=fnlen), dimension(:)  , allocatable :: avgnames, avgfields, &
                                                     analyzernames
  namelist /collect_config/ avgnames
  namelist /xaver_config/ avgfields, ndim2read, analyzernames
  namelist /zaver_config/ avgfields, ndim2read, analyzernames

  integer :: hdferr
  integer(HID_T) :: hdffileid, hdfavggroup, hdffieldgroup, hdfspace, hdfdataset, hdfmemtype
  integer(HID_T) :: phdf_fileid, phdf_avggroup, phdf_memtype, phdf_dataspace
  integer(HID_T), dimension(:), allocatable :: phdf_fieldgroups
  integer(HID_T), dimension(:,:), allocatable :: phdf_spaces, phdf_chunkspaces, phdf_datasets, &
                                                 phdf_plist_ids

  integer(HSIZE_T) , dimension(:,:), allocatable :: phdf_dims, phdf_chunkdims
  
  integer(HID_T) :: plist_id
                    
  integer(HSIZE_T) , dimension(3) :: hdfdims
  integer(HSIZE_T) , dimension(3) :: phdf_stride, phdf_count, phdf_block, phdf_offsets
  logical :: hdfexists, phdf_exists
  
  
  ! Internal parameters
  integer, parameter            :: len_avgnames = 30, len_avgfields = 100, &
                                   len_analyzers = 30, &
                                   rkind = selected_real_kind(10,40)

  ! Job communication parameters and variables
  integer, parameter            :: command_shutdown   = 0, &
                                   command_analyze    = 1, &
                                   command_hdfopen    = 2, &
                                   command_hdfclose   = 3 
             
  integer  :: nprocs_needed, analysisproc, command, fake_proc, nprocs, &
              mpierr, mpistat

  logical :: initerror, runerror, received, filecheck
  
  integer, dimension(:), allocatable   :: fake_procs
  
  character (len=intlen) :: chproc
  
  ! Timing variables
  real(kind=rkind)     :: t1,t2,analysisstart

  ! Data read variables
  integer :: iavg, ierr, ntimesteps, navgs, ndim1, ndim2, &
             ndim1_full, ndim2_full, &
             naverages, ndim2read, i, j, k, nanalyzers, &
             ianalyzer
  
  integer, dimension(:), allocatable :: avgdims
  integer, dimension(:,:), allocatable :: offsets

  real, dimension(:,:,:,:)  , allocatable   :: tmparray2
  real, dimension(:,:,:)    , allocatable   :: dataarray, tmparray, tmparray3
  real, dimension(:)        , allocatable   :: t_values

  integer               :: filesize, datalen, tlen, data_stride, &
                           data_start, tsteplen, dim2stride, pos, &
                           avgnumber, averagelen, dim2runs, resultlen
  
  integer , parameter   :: t_start = 5

  ! Analyzer parameters
  procedure(AnalyzerTemplate), pointer :: analyzerfunction
  character(len=fnlen) :: analyzername

  nullify(analyzerfunction)

  initerror = .false.
  runerror = .false.
  received = .false.


! Initialize MPI

  !call mpicomm_init()
  !call initialize_messages()
  !call initialize_mpicomm()
  call MPI_INIT(mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, iproc , mpierr)
  analysisstart = MPI_WTIME() 

  lroot = (iproc==root)

!
!  Identify version.
!
  if (lroot) call svn_id( &
      '$Id: pc_meanfield_collect.f90 23201 2015-04-06 02:11:37Z st.tuomisto@gmail.com $')
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
  
  call directory_names()

  if (nprocs < 2) then
    if (lroot) then
        write (*,*) 'Program needs more than one process'
        initerror = .true.
    end if
  end if

! Initialize HDF5

  call H5open_F(hdferr)
  call H5Pcreate_F(H5P_FILE_ACCESS_F, plist_id, hdferr)
  call H5Pset_fapl_mpio_F(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
  phdffile = trim(datadir_snap)//'/mfcollect.h5'
  ! Open hdf5-file
  inquire(file=hdffile, exist=phdf_exists)
  if (phdf_exists) then
    open(11, file=phdffile, action='write', status='replace', iostat=ierr)
    close(11, status='delete')
  end if
  call H5Fcreate_F(phdffile, H5F_ACC_TRUNC_F, phdf_fileid, hdferr, access_prp = plist_id)
! Root part of the analysis, read config and co-ordinate other programs

  if (lroot) then

    write (*,'(a,i1,a)') ' This is a ', dimensionality, '-D run'
    print *, 'nxgrid, nygrid, nzgrid=', nxgrid, nygrid, nzgrid
    print *, 'Lx, Ly, Lz=', Lxyz
    print *, '      Vbox=', Lxyz(1)*Lxyz(2)*Lxyz(3)
    
    ! Collect information on 

    allocate(avgnames(len_avgnames))
    allocate(avgfields(len_avgfields))
    allocate(analyzernames(len_analyzers))
    allocate(fake_procs(nprocx*nprocy*nprocz))
    allocate(offsets(2,nprocx*nprocy*nprocz))
    avgnames    = ''
    avgfields   = ''
    analyzernames   = ''
    fake_procs  = 0

    open(UNIT=1, FILE=mfcollect_in_file, IOSTAT=ierr)
    if (ierr /= 0) then
      write(*,*) 'Error opening '//mfcollect_in_file
      initerror = .true.
    else
      read(1, NML=collect_config, IOSTAT=ierr)
      if (ierr /= 0) then
        write(*,*) 'Error reading configuration'
        initerror = .true.
      end if
    end if

    call MPI_BCAST(initerror, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
    if (.not. initerror) then
      analysisproc = 1
      do iavg=1,len_avgnames

        avgname = avgnames(iavg)
        ! Field dependent part
        ! ----------------------------------------
        nprocs_needed = 0
        if (avgname == 'xaver') then
          ! This needs to be written based on zaver example 
        else if (avgname == 'zaver') then
          ndim1 = nx
          ndim2 = ny
          ndim1_full = nxgrid
          ndim2_full = nygrid
          ndim2read=ny
          infile='zaver.in'
          datafilename='zaverages.dat'

          ! Read configuration
          read(1, NML=zaver_config, IOSTAT=ierr)
          if (ierr /= 0) then
            write(*,*) 'Error reading zaver configuration'
            initerror = .true.
          end if
      
          if ((ndim2read < 0) .or. (modulo(ndim2,ndim2read) /= 0)) then
            write(*,*) 'Invalid ndim2read given, should divide second dimension of average.'
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
              offsets(:,j) = [ndim1*ipx, ndim2*ipy]
              j = j + 1
            end if
          end do
        else
          exit
        end if
         
        ! Universal part
        ! ----------------------------------------
        
        ! Determine picked fields
        navgs = 0
        do i=1,len_avgfields
          if (avgfields(i) /= '') then
            navgs = navgs + 1
          else
            exit
          end if
        end do
        i = 0
        if (navgs > 0) then
          allocate(avgdims(navgs))
          avgdims = 0
          naverages = 0
          open(UNIT=2, FILE=infile, STATUS='old', IOSTAT=ierr)
          if (ierr /= 0) then
            write(*,*) 'Could not open '//infile
            initerror = .true.
            exit
          else
            do i=1,len_avgfields
              read(2,*,iostat=ierr) avgfield
              if (ierr /= 0) then
                naverages = i-1
                exit
              else 
                avgfield = trim(avgfield)
                do j=1,navgs
                  if (avgfield == trim(avgfields(j))) then
                    avgdims(j)  = i
                  end if
                end do
              end if
            end do
          end if
          close(2)
        else 
          write(*,*) 'No avgfields found.'
          initerror=.true.
        end if

        if (any(avgdims == 0)) then
          write(*,*) 'Could not find some of required fields'
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

        ! Send datasets to create
        if ((nanalyzers>0) .and. (navgs > 0) .and. (nprocs_needed > 0) &
          .and. (.not. initerror)) then
          command = command_hdfopen
          do i=1,nprocs-1
            call MPI_SEND(command, 1, &
                          MPI_INT, i, 0, MPI_COMM_WORLD, mpierr)
          end do
          call BCastInfo()
          call OpenH5Groups()
        end if

        ! Commence the analysis
        if ((nanalyzers>0) .and. (navgs > 0) .and. (nprocs_needed > 0) &
          .and. (.not. initerror)) then
          command = command_analyze
          do i=1, nprocs_needed
            fake_proc = fake_procs(i)
            call MPI_SEND(command, 1, &
                          MPI_INT, analysisproc, 0, MPI_COMM_WORLD, mpierr)
            call MPI_SEND(fake_proc, 1, &
                          MPI_INT, analysisproc, 1, MPI_COMM_WORLD, mpierr)
            call MPI_SEND(offsets(:,i), 2, &
                          MPI_INT, analysisproc, 2, MPI_COMM_WORLD, mpierr)
            call MPI_RECV(received, 1, &
                          MPI_LOGICAL, analysisproc, 3, MPI_COMM_WORLD, mpistat, mpierr)
            analysisproc = analysisproc+1
            if (analysisproc>=nprocs) then
              analysisproc = 1
            end if
          end do
          command = command_hdfclose
          do i=1,nprocs-1
            call MPI_SEND(command, 1, MPI_INT, analysisproc, 0, MPI_COMM_WORLD, mpierr)
            analysisproc = analysisproc + 1
            if (analysisproc>=nprocs) then
              analysisproc = 1
            end if
            call CloseH5Groups()
          end do
        end if

        if (allocated(avgdims)) then
            deallocate(avgdims)
        end if
      end do
      ! Send shutdown to other fields
      command = command_shutdown
      do i=1,nprocs-1
        call MPI_SEND(command, 1, MPI_INT, analysisproc, 0, MPI_COMM_WORLD, mpierr)
        analysisproc = analysisproc + 1
        if (analysisproc>=nprocs) then
          analysisproc = 1
        end if
      end do
    end if
    close(1)
    if (allocated(avgnames)) then
      deallocate(avgnames)
    end if
    if (allocated(avgfields)) then
      deallocate(avgfields)
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
                      MPI_INT, root, 0, MPI_COMM_WORLD, mpistat, mpierr)
        if (command == command_hdfopen) then
          call BCastInfo()
          call OpenH5Groups()
        else if (command == command_hdfclose) then
          call CloseH5Groups()
        else if (command == command_analyze) then
          ! MPI communication part
          call MPI_RECV(fake_proc, 1, &
                        MPI_INT, root, 1, MPI_COMM_WORLD, mpistat, mpierr)
          if (.not. allocated(offsets)) then
            allocate(offsets(2,1))
          end if
          call MPI_RECV(offsets(:,1), 2, &
                        MPI_INT, root, 2, MPI_COMM_WORLD, mpistat, mpierr)
          call MPI_SEND(received, 1, &
                        MPI_LOGICAL, root, 3, MPI_COMM_WORLD, mpierr)
          chproc = itoa (fake_proc)
          call safe_character_assign (directory_dist, &
                                                  trim (datadir_snap)//'/proc'//chproc)
          call safe_character_assign (datafile, &
                                                  trim(directory_dist)//'/'//trim(datafilename))
          if (lprocz_slowest) then
            ipx = modulo(fake_proc, nprocx)
            ipy = modulo(fake_proc/nprocx, nprocy)
            ipz = fake_proc/nprocxy
          else
            ipx = modulo(fake_proc, nprocx)
            ipy = fake_proc/nprocxy
            ipz = modulo(fake_proc/nprocx, nprocy)
          endif

          ! Data reading part
          ! -------------------------------------------------------------

          ! Open logfile
          logfile = trim(directory_dist)//'/mfcollect.log'
          open(2, file=logfile, action='write', status='replace', iostat=ierr)
          runerror = .false.

          ! Open hdf5-file
          hdffile = trim(directory_dist)//'/mfcollect.h5'
          inquire(file=hdffile, exist=hdfexists)
          if (hdfexists) then
            open(11, file=hdffile, action='write', status='replace', iostat=ierr)
            close(11, status='delete')
          end if

          call H5Fcreate_F(hdffile, H5F_ACC_TRUNC_F, hdffileid, hdferr)
          if (hdferr /= 0) then
            write(*,*) 'Could not open file:'//hdffile
            runerror = .true.
          end if

          call H5Gcreate_F(hdffileid, "/"//trim(avgname), hdfavggroup, hdferr)
          do iavg=1,navgs
            call H5Gcreate_F(hdfavggroup, "/"//trim(avgname)//"/"//avgfields(iavg), hdffieldgroup, hdferr)
            call H5Gclose_F(hdffieldgroup, hdferr)
          end do
          ! Check datafile properties

          open(1,file=datafile, access='stream',form='unformatted', action='read', status='old', iostat=ierr)
          if (ierr /= 0) then
            write(2,*) 'Error opening datafile: '//datafile
            runerror = .true.
          end if

          filesize = -1
          tlen    = -1
          datalen     = -1
          
          inquire(1, size=filesize)
          pos=1
          read(1, pos=pos, IOSTAT=ierr) tlen
          if (ierr /= 0) then
            write(2,*) 'Error reading tlen'
            runerror = .true.
          end if
          pos=9+tlen
          read(1, pos=pos, IOSTAT=ierr) datalen
          if (ierr /= 0) then
            write(2,*) 'Error reading datalen'
            runerror = .true.
          end if
          data_stride = datalen/(naverages*ndim1*ndim2)
          if (data_stride == 4) then
            hdfmemtype = H5T_NATIVE_REAL
          else if (data_stride == 8) then
            hdfmemtype = H5T_NATIVE_DOUBLE
          else
            write(*,*) 'Could not determine the data type.'
            runerror = .true.
          end if
          dim2stride  = data_stride*ndim1
          ntimesteps  = int(floor(real(filesize)/(16+tlen+datalen)))
          data_start  = t_start + tlen + 8
          tsteplen    = 16 + tlen + datalen
          averagelen  = datalen/naverages 
          dim2runs   = ndim2/ndim2read

          ! Log analysis parameters

          write(2,'(a30)')     ' Analysis information          '
          write(2,'(a30)')     ' ----------------------------- '
          write(2,'(a30,a30)') &
          ' Worker:                       ', ljustifyI(fake_proc)
          write(2,'(a30,a30)') &
          ' Dimension 1 size:             ', ljustifyI(ndim1)
          write(2,'(a30,a30)') &
          ' Dimension 2 size:             ', ljustifyI(ndim2)
          write(2,'(a30,a)') &
          ' Datafile:                     ', trim(datafile)
          write(2,'(a30,a30)') &
          ' Number of averages:           ', ljustifyI(naverages)
          write(2,'(a30,a30)') &
          ' File size (bytes):            ', ljustifyI(filesize)
          write(2,'(a30,a30)') &
          ' Number of timesteps:          ', ljustifyI(ntimesteps)
          write(2,'(a30,a30)') &
          ' Datapoints per timestep:      ', ljustifyI(datalen/data_stride)
          write(2,'(a30,a30)') &
          ' Single precision time values: ', ljustifyL(tlen==4)
          write(2,'(a30,a30)') &
          ' Single precision data values: ', ljustifyL(data_stride==4)
          close(1)

          allocate(t_values(ntimesteps),stat=ierr) 
          if (ierr /= 0) then
            write(2,*) 'Error allocating t_values'
            runerror = .true.
          end if
          allocate(tmparray(ndim1,ndim2read,ntimesteps),stat=ierr)
          if (ierr /= 0) then
            write(2,*) 'Error allocating tmparray'
            runerror = .true.
          end if
          allocate(tmparray2(ndim1,ndim2,naverages,ntimesteps),stat=ierr)
          if (ierr /= 0) then
            write(2,*) 'Error allocating tmparray2'
            runerror = .true.
          end if
          allocate(tmparray3(ndim1,ndim2,ntimesteps),stat=ierr)
          if (ierr /= 0) then
            write(2,*) 'Error allocating tmparray3'
            runerror = .true.
          end if

          ! Test values
          ! -------------------------
          outputfile = trim(directory_dist)//'/testvalues_test.dat'
          open(3, file=outputfile, action='write', status='replace', IOSTAT=ierr)

          t_values    = 0
          
          !write(*,*) 'Test output:'
          write(2,*) 'Test output:'
          
          t1 = MPI_WTIME()
          open(1,file=datafile, form='unformatted', action='read', status='old', IOSTAT=ierr)
          do i=1,ntimesteps
            read(1, IOSTAT=ierr) t_values(i)
            read(1, IOSTAT=ierr) tmparray2(:,:,:,i)
          end do
          
          tmparray3 = tmparray2(:,:,avgdims(1),:)
          do ianalyzer=1,nanalyzers
            nullify(analyzerfunction)
            call getAnalyzer(analyzernames(ianalyzer),analyzerfunction,resultlen)
            if (.not. allocated(dataarray)) then
              allocate(dataarray(ndim1,ndim2,resultlen),stat=ierr)
              !write(*,*) 'Allocated dataarray, dim3 is' , size(dataarray,dim=3)
              if (ierr /= 0) then
                write(2,*) 'Error allocating dataarray'
                runerror = .true.
                exit
              end if
            else if (resultlen /= size(dataarray,dim=3)) then
              deallocate(dataarray)
              allocate(dataarray(ndim1,ndim2,resultlen),stat=ierr)
              !write(*,*) size(dataarray,dim=3)
              if (ierr /= 0) then
                write(2,*) 'Error allocating dataarray'
                runerror = .true.
                exit
              end if
            end if
            hdfdims = [ ndim1, ndim2, resultlen ]
            dataarray = analyzerfunction(tmparray3, ndim1, ndim2, ntimesteps, resultlen)
            
!            do iavg=1,navgs
!              do j=1,ndim2
!                write(2,*) tmparray2(10:15,j,avgdims(iavg),1)
!                write(3,*) tmparray2(10:15,j,avgdims(iavg),1)
!              end do
!            end do
            do iavg=1,navgs
              do j=1,ndim2
                write(3,*) dataarray(1:ndim1,j,1)
              end do
              call H5Gopen_F(hdfavggroup, "/"//trim(avgname)//"/"//avgfields(iavg), hdffieldgroup, hdferr)
              call H5Screate_simple_F(3, hdfdims, hdfspace, hdferr)
              call H5Dcreate_F(hdffieldgroup, trim(analyzernames(ianalyzer))//'_test',H5T_IEEE_F32LE,hdfspace,hdfdataset,hdferr)
              call H5Dwrite_F(hdfdataset, hdfmemtype, dataarray, hdfdims, hdferr)
              call H5Dclose_F(hdfdataset,hdferr)
              call H5Sclose_F(hdfspace, hdferr)
              call H5Gclose_F(hdffieldgroup, hdferr)

              ! Parallel version
              call H5Dget_space_F(phdf_datasets(ianalyzer,iavg), phdf_dataspace, hdferr)
              phdf_offsets  = [ offsets(1,1), offsets(2,1), 0 ]
              phdf_count    = [ 1, 1, 1 ]
              phdf_stride   = [ 1, 1, 1 ]
              phdf_block    = [ ndim1, ndim2, resultlen ]
              write(*,*) fake_proc, phdf_offsets
              call H5Sselect_hyperslab_F(phdf_dataspace, H5S_SELECT_SET_F, phdf_offsets, phdf_count, &
                                         hdferr,phdf_stride, phdf_block)
              call H5Dwrite_F(phdf_datasets(ianalyzer,iavg), hdfmemtype, dataarray, &
                                         phdf_dims(:,ianalyzer), hdferr, &
                                         file_space_id  = phdf_dataspace, &
                                         mem_space_id   = phdf_chunkspaces(ianalyzer, iavg))

              call H5Sclose_F(phdf_dataspace, hdferr)
            end do

          end do
          close(1)
          close(3)
          write(2,*) t_values(:ntimesteps)
          if (allocated(dataarray)) then
            deallocate(dataarray)
          end if

          t1 = Timer(t1)
          write(*,*) 'Time taken: ' , t1
          write(2,*) 'Time taken: ' , t1
          
          ! Stream I/O values
          ! -------------------------
          outputfile = trim(directory_dist)//'/testvalues_stream.dat'
          open(3, file=outputfile, action='write', status='replace', IOSTAT=ierr)

          t_values    = 0

          !write(*,*) 'Stream I/O output:'
          write(2,*) 'Stream I/O output:'
          
          t2 = MPI_WTIME()
          open(1,file=datafile, access='stream',form='unformatted', &
              action='read', status='old', IOSTAT=ierr)
          
          pos = t_start
          do i=1,ntimesteps
            read(1, pos=pos, IOSTAT=ierr) t_values(i)
            pos = pos + tsteplen
          end do
          
          do ianalyzer=1,nanalyzers
            nullify(analyzerfunction)
            call getAnalyzer(analyzernames(ianalyzer),analyzerfunction,resultlen)
            if (.not. allocated(dataarray)) then
              allocate(dataarray(ndim1,ndim2,resultlen),stat=ierr)
              !write(*,*) 'Allocated dataarray, dim3 is' , size(dataarray,dim=3)
              if (ierr /= 0) then
                write(2,*) 'Error allocating dataarray'
                runerror = .true.
                exit
              end if
            else if (resultlen /= size(dataarray,dim=3)) then
              deallocate(dataarray)
              allocate(dataarray(ndim1,ndim2,resultlen),stat=ierr)
              !write(*,*) size(dataarray,dim=3)
              if (ierr /= 0) then
                write(2,*) 'Error allocating dataarray'
                runerror = .true.
                exit
              end if
            end if
            do iavg=1,navgs
              do j=1,dim2runs
                pos = data_start+(avgdims(iavg)-1)*averagelen+(j-1)*dim2stride*ndim2read
                do i=1,ntimesteps
                  read(1, pos=pos, IOSTAT=ierr) tmparray(:,:,i)
                  pos = pos + tsteplen
                end do
                !write(*,*) j, (j-1)*ndim2read+1, j*ndim2read
                dataarray(:,(j-1)*ndim2read+1:j*ndim2read,:) = analyzerfunction(tmparray, ndim1, ndim2read, ntimesteps, resultlen)
!                do k=1,ndim2read
!                  write(2,*) tmparray(10:15,k,1)
!                  write(3,*) tmparray(10:15,k,1)
!                end do
              end do
            end do
            do iavg=1,navgs
              do j=1,ndim2
                write(3,*) dataarray(1:ndim1,j,1)
              end do
              call H5Gopen_F(hdfavggroup, "/"//trim(avgname)//"/"//avgfields(iavg), hdffieldgroup, hdferr)
              call H5Screate_simple_F(3, hdfdims, hdfspace, hdferr)
              call H5Dcreate_F(hdffieldgroup, trim(analyzernames(ianalyzer)),H5T_IEEE_F32LE,hdfspace,hdfdataset,hdferr)
              call H5Dwrite_F(hdfdataset, hdfmemtype, dataarray, hdfdims, hdferr)
              call H5Dclose_F(hdfdataset,hdferr)
              call H5Sclose_F(hdfspace, hdferr)
              call H5Gclose_F(hdffieldgroup, hdferr)
            end do
          end do
          write(2,*) t_values(:ntimesteps)
          
          if (allocated(dataarray)) then
            deallocate(dataarray)
          end if

          t2 = Timer(t2)
          write(*,*) 'Time taken: ' , t2
          write(2,*) 'Time taken: ' , t2
          
          write(*,*) 'Speedup: ', t1/t2
          write(2,*) 'Speedup: ', t1/t2

          close(1)
          close(3)

          !call H5SCLOSE_F(hdfspace, hdferr)
          call H5Gclose_F(hdfavggroup, hdferr)
          call H5Fclose_F(hdffileid, hdferr)
          
          if (allocated(dataarray)) then
            deallocate(dataarray)
          end if
          if (allocated(t_values)) then
            deallocate(t_values)
          end if
          if (allocated(tmparray)) then
            deallocate(tmparray)
          end if
          if (allocated(tmparray2)) then
            deallocate(tmparray2)
          end if
          if (allocated(tmparray3)) then
            deallocate(tmparray3)
          end if
          nullify(analyzerfunction)
          close(2)
        else
          write(*,*) 'Received shutdown command at ', iproc
          exit
        end if
      end do
    end if
  end if

  write(*,*) 'Time taken: ', Timer(analysisstart)

  call H5Pclose_F(plist_id,hdferr)
  call H5Fclose_F(phdf_fileid, hdferr)
  call fnames_clean_up()
  call vnames_clean_up()
  call H5close_f(hdferr)
  call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
  call MPI_FINALIZE(mpierr)

  contains

    subroutine BCastInfo()
      call MPI_BCAST(navgs, 1, &
                MPI_INT,  0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(naverages, 1, &
                MPI_INT,  0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(nanalyzers, 1, &
                MPI_INT,  0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(ndim1, 1, &
                MPI_INT,  0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(ndim2, 1, &
                MPI_INT,  0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(ndim2read, 1, &
                MPI_INT,  0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(ndim1_full, 1, &
                MPI_INT,  0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(ndim2_full, 1, &
                MPI_INT,  0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(avgname, fnlen, &
                MPI_CHAR, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(datafilename, fnlen, &
                MPI_CHAR, 0, MPI_COMM_WORLD, mpierr)
      write(*,*) iproc, navgs, naverages, nanalyzers, ndim1, ndim2, ndim2read, trim(avgname), ' ', trim(datafilename)
      if (.not. lroot) then
        if (allocated(avgdims)) then
          deallocate(avgdims)
        end if
        allocate(avgdims(navgs))
        if (allocated(avgfields)) then
          deallocate(avgfields)
        end if
        allocate(avgfields(navgs))
        if (allocated(analyzernames)) then
          deallocate(analyzernames)
        end if
        allocate(analyzernames(nanalyzers))
      end if
      call MPI_BCAST(avgdims, navgs, &
                MPI_INT,  0, MPI_COMM_WORLD, mpierr)
      do j=1,navgs
        call MPI_BCAST(avgfields(j), fnlen, &
                  MPI_CHAR, 0, MPI_COMM_WORLD, mpierr)
      end do
      do j=1,nanalyzers
        call MPI_BCAST(analyzernames(j), fnlen, &
                  MPI_CHAR, 0, MPI_COMM_WORLD, mpierr)
      end do
      filecheck = .true.
      call MPI_REDUCE(filecheck, received, 1, &
                      MPI_LOGICAL, MPI_LAND, 0, MPI_COMM_WORLD, mpistat, mpierr)

    end subroutine BCastInfo

    subroutine OpenH5Groups()
      ! Create dataset

      if (allocated(phdf_fieldgroups)) then
        deallocate(phdf_fieldgroups)
      end if
      if (allocated(phdf_spaces)) then
        deallocate(phdf_spaces)
      end if
      if (allocated(phdf_chunkspaces)) then
        deallocate(phdf_chunkspaces)
      end if
      if (allocated(phdf_datasets)) then
        deallocate(phdf_datasets)
      end if
      if (allocated(phdf_plist_ids)) then
        deallocate(phdf_plist_ids)
      end if
      if (allocated(phdf_dims)) then
        deallocate(phdf_dims)
      end if
      if (allocated(phdf_chunkdims)) then
        deallocate(phdf_chunkdims)
      end if
      allocate(phdf_fieldgroups(navgs))
      allocate(phdf_spaces(nanalyzers,navgs))
      allocate(phdf_datasets(nanalyzers,navgs))
      allocate(phdf_chunkspaces(nanalyzers,navgs))
      allocate(phdf_plist_ids(nanalyzers,navgs))
      allocate(phdf_dims(3,naverages))
      allocate(phdf_chunkdims(3,naverages))
      do ianalyzer=1,nanalyzers
        call getAnalyzer(analyzernames(ianalyzer),analyzerfunction,resultlen)
        phdf_dims(:,ianalyzer) = [ndim1_full, ndim2_full, resultlen]
        phdf_chunkdims(:,ianalyzer) = [ndim1, ndim2, resultlen]
      end do
      
      call H5Gcreate_F(phdf_fileid, "/"//trim(avgname), phdf_avggroup, hdferr)
      if (lroot) then
        write(*,*) 'CreatedG:', "/"//trim(avgname)
      end if
      do iavg=1,navgs
        call H5Gcreate_F(phdf_avggroup, "/"//trim(avgname)//"/"//avgfields(iavg), &
                       phdf_fieldgroups(iavg), hdferr)
        if (lroot) then
          write(*,*) 'CreatedG:', "/"//trim(avgname)//"/"//avgfields(iavg)
        end if
        do ianalyzer=1,nanalyzers
          write(*,*) phdf_dims(:,ianalyzer)
          write(*,*) phdf_chunkdims(:,ianalyzer)
          call H5Screate_simple_F(3, phdf_dims(:,ianalyzer), phdf_spaces(ianalyzer,iavg), hdferr)
          call H5Screate_simple_F(3, phdf_chunkdims(:,ianalyzer), phdf_chunkspaces(ianalyzer,iavg), hdferr)
          if (lroot) then
            write(*,*) 'CreatedS:', iavg, ianalyzer
          end if
          call H5Pcreate_F(H5P_DATASET_CREATE_F, phdf_plist_ids(ianalyzer,iavg), hdferr)
          call H5Pset_chunk_f(phdf_plist_ids(ianalyzer, iavg), 3, phdf_chunkdims(:,ianalyzer), hdferr)
          call H5Dcreate_F(phdf_fieldgroups(iavg), trim(analyzernames(ianalyzer)), &
                           H5T_IEEE_F32LE, phdf_spaces(ianalyzer,iavg), &
                           phdf_datasets(ianalyzer,iavg), hdferr, phdf_plist_ids(ianalyzer,iavg))
          
          call H5Sclose_F(phdf_spaces(ianalyzer,iavg), hdferr)
          if (lroot) then
            write(*,*) 'CreatedD:', iavg, ianalyzer
          end if
        end do
      end do

      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)

    end subroutine OpenH5Groups

    subroutine CloseH5Groups

      do iavg=1,navgs
        do ianalyzer=1,nanalyzers
          call H5Dclose_F(phdf_datasets(ianalyzer,iavg),hdferr)
          call H5Pclose_F(phdf_plist_ids(ianalyzer,iavg), hdferr)
          if (lroot) then
            write(*,*) 'ClosedD:', iavg, ianalyzer
          end if
          call H5Sclose_F(phdf_chunkspaces(ianalyzer,iavg), hdferr)
          if (lroot) then
            write(*,*) 'ClosedS:', iavg, ianalyzer
          end if
        end do
        call H5Gclose_F(phdf_fieldgroups(iavg), hdferr)
        if (lroot) then
          write(*,*) 'ClosedG:', "/"//trim(avgname)//"/"//avgfields(iavg)
        end if
      end do
      call H5Gclose_F(phdf_avggroup, hdferr)
      if (lroot) then
        write(*,*) 'ClosedG:', "/"//trim(avgname)
      end if

    end subroutine CloseH5Groups

    real(kind=rkind) function Timer(t)
        real(kind=rkind),intent(in) :: t
        Timer = MPI_WTIME()-t
        return
    end function Timer

    character(len=fnlen) function ljustifyI(i)
        integer, intent(in) :: i
        write(ljustifyI, '(i20)') i
        ljustifyI = adjustl(ljustifyI)
        return
    end function ljustifyI

    character(len=fnlen) function ljustifyL(l)
        logical, intent(in) :: l
        write(ljustifyL, '(l20)') l
        ljustifyL = adjustl(ljustifyL)
        return
    end function ljustifyL

end program pc_meanfield_collect
