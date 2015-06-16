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

  character (len=fnlen) :: datafile, inputfile, avgname, logfile, &
                            outputfile, hdffile, avgfield
  character (len=*), parameter :: directory_out = 'data/allprocs'
  character (len=*), parameter :: mfcollect_in_file = 'mfcollect.in'

  character (len=fnlen), dimension(:)  , allocatable :: avgnames, avgfields, &
                                                     analyzernames
  namelist /collect_config/ avgnames
  namelist /xaver_config/ avgfields, ndim2read, analyzernames
  namelist /zaver_config/ avgfields, ndim2read, analyzernames

  integer :: hdferr
  integer(HID_T) :: hdffileid, hdfavggroup, hdffieldgroup, hdfspace, hdfdataset, hdfmemtype
                    
  integer(HSIZE_T), dimension(3) :: hdfdims
  logical :: hdfexists
  
  
  ! Internal parameters
  integer, parameter            :: len_avgnames = 30, len_avgfields = 100, &
                                   len_analyzers = 30, &
                                   rkind = selected_real_kind(10,40)

  ! Job communication parameters and variables
  integer, parameter            :: command_shutdown = 0, &
                                   command_analyze = 1
             
  integer  :: nprocs_needed, analysisproc, command, fake_proc, nprocs, &
              mpierr, mpistat

  logical :: initerror, runerror, received, filecheck
  
  integer, dimension(:), allocatable :: fake_procs
  
  character (len=intlen) :: chproc
  
  ! Timing variables
  real(kind=rkind)     :: t1,t2,analysisstart

  ! Data read variables
  integer :: iavg, ierr, ntimesteps, navgs, ndim1, ndim2, &
             naverages, ndim2read, i, j, k, nanalyzers, &
             ianalyzer
  
  integer, dimension(:), allocatable :: avgdims

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

! Initialize HDF5

  call H5open_f(hdferr)

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
          ! This needs to be written based on zaver 
        else if (avgname == 'zaver') then
          ndim1 = nx
          ndim2 = ny
          ndim2read=ny
          inputfile='zaver.in'
          datafile='zaverages.dat'

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
          open(UNIT=2, FILE=inputfile, STATUS='old', IOSTAT=ierr)
          if (ierr /= 0) then
            write(*,*) 'Could not open '//inputfile
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

        ! Commence the analysis
        if ((nanalyzers>0) .and. (navgs > 0) .and. (nprocs_needed > 0) &
          .and. (.not. initerror)) then
          command = command_analyze
          do i=1, nprocs_needed
            fake_proc = fake_procs(i)
            write(*,'(a30,a5)') &
            ' Send 1: command to            ',ljustifyI(analysisproc)
            call MPI_SEND(command, 1, &
                          MPI_INT, analysisproc, 0, MPI_COMM_WORLD, mpierr)
            write(*,'(a30,a5)') &
            ' Send 2: fake_proc to          ',ljustifyI(analysisproc)
            call MPI_SEND(fake_proc, 1, &
                          MPI_INT, analysisproc, 1, MPI_COMM_WORLD, mpierr)
            write(*,'(a30,a5)') &
            ' Send 3: naverages to          ',ljustifyI(analysisproc)
            call MPI_SEND(naverages, 1, &
                          MPI_INT, analysisproc, 2, MPI_COMM_WORLD, mpierr)
            write(*,'(a30,a5)') &
            ' Send 4: ndim1 to              ',ljustifyI(analysisproc)
            call MPI_SEND(ndim1, 1, &
                          MPI_INT, analysisproc, 3, MPI_COMM_WORLD, mpierr)
            write(*,'(a30,a5)') &
            ' Send 5: ndim2 to              ',ljustifyI(analysisproc)
            call MPI_SEND(ndim2, 1, &
                          MPI_INT, analysisproc, 4, MPI_COMM_WORLD, mpierr)
            write(*,'(a30,a5)') &
            ' Send 5: ndim2read to          ',ljustifyI(analysisproc)
            call MPI_SEND(ndim2read, 1, &
                          MPI_INT, analysisproc, 5, MPI_COMM_WORLD, mpierr)
            write(*,'(a30,a5)') &
            ' Send 6: datafile to           ',ljustifyI(analysisproc)
            call MPI_SEND(datafile, fnlen, &
                          MPI_CHAR, analysisproc, 6, MPI_COMM_WORLD, mpierr)
            write(*,'(a30,a5)') &
            ' Send 7: navgs to              ',ljustifyI(analysisproc)
            call MPI_SEND(navgs, 1, &
                          MPI_INT, analysisproc, 7, MPI_COMM_WORLD, mpierr)
            do j=1,navgs
              write(*,'(a8,i1,a20,a5)') &
              ' Send 8-',j,': average to          ',ljustifyI(analysisproc)
              call MPI_SEND(avgfields(j), fnlen, &
                            MPI_CHAR, analysisproc, 8*j, MPI_COMM_WORLD, mpierr)
            end do
            write(*,'(a30,a5)') &
            ' Send 9: navgs to              ',ljustifyI(analysisproc)
            call MPI_SEND(avgdims, navgs, &
                          MPI_INT, analysisproc, 9, MPI_COMM_WORLD, mpierr)
            write(*,'(a30,a5)') &
            ' Send 10: nanalyzers to         ',ljustifyI(analysisproc)
            call MPI_SEND(nanalyzers, 1, &
                          MPI_INT, analysisproc, 10, MPI_COMM_WORLD, mpierr)
            do j=1,nanalyzers
              write(*,'(a9,i1,a20,a5)') &
              ' Send 11-',j,': analyzer to         ',ljustifyI(analysisproc)
              call MPI_SEND(analyzernames(j), fnlen, &
                            MPI_CHAR, analysisproc, 11*j, MPI_COMM_WORLD, mpierr)
            end do
            call MPI_SEND(avgname, fnlen, &
                          MPI_CHAR, analysisproc, 12, MPI_COMM_WORLD, mpierr)
            call MPI_RECV(received, 1, &
                          MPI_LOGICAL, analysisproc, 0, MPI_COMM_WORLD, mpistat, mpierr)
            write(*,'(a30,a5)') &
            ' Receive last: received from   ',ljustifyI(analysisproc)
            analysisproc = analysisproc+1
            if (analysisproc>=nprocs) then
              analysisproc = 1
            end if
          end do
        end if
        if (allocated(avgdims)) then
            deallocate(avgdims)
        end if
      end do
      ! Send shutdown to other fields
      command = command_shutdown
      do i=1,nprocs-1
        write(*,*) 'Sent ', command, ' to ', analysisproc
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
      command = 0
      do
        call MPI_RECV(command, 1, &
                      MPI_INT, root, 0, MPI_COMM_WORLD, mpistat, mpierr)
        if (command == command_analyze) then
          ! MPI communication part
          call MPI_RECV(fake_proc, 1, &
                        MPI_INT, root, 1, MPI_COMM_WORLD, mpistat, mpierr)
          call MPI_RECV(naverages, 1, &
                        MPI_INT, root, 2, MPI_COMM_WORLD, mpistat, mpierr)
          call MPI_RECV(ndim1, 1, &
                        MPI_INT, root, 3, MPI_COMM_WORLD, mpistat, mpierr)
          call MPI_RECV(ndim2, 1, &
                        MPI_INT, root, 4, MPI_COMM_WORLD, mpistat, mpierr)
          call MPI_RECV(ndim2read, 1, &
                        MPI_INT, root, 5, MPI_COMM_WORLD, mpistat, mpierr)
          call MPI_RECV(datafile, fnlen, &
                        MPI_CHAR, root, 6, MPI_COMM_WORLD, mpistat, mpierr)
          call MPI_RECV(navgs, 1, &
                        MPI_INT, root, 7, MPI_COMM_WORLD, mpistat, mpierr)
          allocate(avgdims(navgs))
          allocate(avgfields(navgs))
          do j=1,navgs
            call MPI_RECV(avgfields(j), fnlen, &
                          MPI_CHARACTER, root, 8*j, MPI_COMM_WORLD, mpistat, mpierr)
          end do
          call MPI_RECV(avgdims, navgs, &
                        MPI_INT, root, 9, MPI_COMM_WORLD, mpistat, mpierr)
          call MPI_RECV(nanalyzers, 1, &
                        MPI_INT, root, 10, MPI_COMM_WORLD, mpistat, mpierr)
          allocate(analyzernames(nanalyzers))
          do j=1,nanalyzers
            call MPI_RECV(analyzernames(j),fnlen, &
                          MPI_CHARACTER, root, 11*j, MPI_COMM_WORLD, mpistat, mpierr)
          end do
          call MPI_RECV(avgname, fnlen, &
                        MPI_CHAR, root, 12, MPI_COMM_WORLD, mpistat, mpierr)
          call MPI_SEND(received, 1, &
                        MPI_LOGICAL, root, 0, MPI_COMM_WORLD, mpierr)
          chproc = itoa (fake_proc)
          call safe_character_assign (directory_dist, &
                                                  trim (datadir_snap)//'/proc'//chproc)
          call safe_character_assign (datafile, &
                                                  trim(directory_dist)//'/'//trim(datafile))

          
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
          inquire(file=hdffile, exist=hdfexists)

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
          
          write(*,*) 'Test output:'
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
            call setAnalyzer(analyzernames(ianalyzer),analyzerfunction,resultlen)
            if (.not. allocated(dataarray)) then
              allocate(dataarray(ndim1,ndim2,resultlen),stat=ierr)
              write(*,*) 'Allocated dataarray, dim3 is' , size(dataarray,dim=3)
              if (ierr /= 0) then
                write(2,*) 'Error allocating dataarray'
                runerror = .true.
                exit
              end if
            else if (resultlen /= size(dataarray,dim=3)) then
              deallocate(dataarray)
              allocate(dataarray(ndim1,ndim2,resultlen),stat=ierr)
              write(*,*) size(dataarray,dim=3)
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

          write(*,*) 'Stream I/O output:'
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
            call setAnalyzer(analyzernames(ianalyzer),analyzerfunction,resultlen)
            if (.not. allocated(dataarray)) then
              allocate(dataarray(ndim1,ndim2,resultlen),stat=ierr)
              write(*,*) 'Allocated dataarray, dim3 is' , size(dataarray,dim=3)
              if (ierr /= 0) then
                write(2,*) 'Error allocating dataarray'
                runerror = .true.
                exit
              end if
            else if (resultlen /= size(dataarray,dim=3)) then
              deallocate(dataarray)
              allocate(dataarray(ndim1,ndim2,resultlen),stat=ierr)
              write(*,*) size(dataarray,dim=3)
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
                write(*,*) j, (j-1)*ndim2read+1, j*ndim2read
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
          if (allocated(analyzernames)) then
            deallocate(analyzernames)
          end if
          if (allocated(avgfields)) then
            deallocate(avgfields)
          end if
          if (allocated(avgdims)) then
            deallocate(avgdims)
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
  call fnames_clean_up()
  call vnames_clean_up()
  call H5close_f(hdferr)
  call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
  call MPI_FINALIZE(mpierr)

  contains

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
