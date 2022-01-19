! vi: syntax=fortran
!
! Code fragment to be used in both read_snap_single and read_snap_double - type agnostic
!
      real :: t_test   ! t in single precision for backwards compatibility
      logical :: ltest,lok,l0,lmail,lmail_
      integer :: len1d,ios,iosr,ii,iia,jj,kk,iip,ind,ind1,kka,kke,jja,jje,iie
      integer :: nprocz_src, nzgrid_src, mz_src, ipz_src
      character(LEN=fnlen) :: readdir, srcdir
      character(LEN=linelen) :: message
      character(LEN=1024) :: mailstr
      character(LEN=20) :: cjobid
      integer, parameter :: lun_input1=89

      if (lserial_io) call start_serialize

      readdir = directory_snap
      lok=.false.; l0=.true.; message=''; lmail=.false.
      if (ivar_omit(1)>0) then
        if (ivar_omit(2)<=0) ivar_omit(2)=nv
      endif
!
      if (lzaver_on_input) call prepare_zaver_on_input
!
      if (lrepair_snap.or.snaplink/='') then
        kka=max(0,ipz-1); kke=min(ipz+1,nprocz-1)
        jja=max(0,ipy-1); jje=min(ipy+1,nprocy-1)
        iie=1
      else
        kka=1; kke=1; jja=1; jje=1
      endif

kloop:do kk=kka,kke
jloop:  do jj=jja,jje
          iia=max(0,ipx-1)
          if (l0) iia=iia-1                           ! add one iteration for the original path
          if (lrepair_snap.or.snaplink/='') then
            iie=min(nprocx-1,ipx+1)
          else
            iie=iia
          endif
iloop:    do ii=iia,iie
            if (.not.l0) then
              iip=find_proc(ii,jj,kk)
              if (iip/=iproc) then
                ind=index(directory_snap,'proc')
                ind1=index(directory_snap(ind:),'/')
                if (ind1>0) then
                  readdir=directory_snap(:ind+3)//trim(itoa(iip))//directory_snap(ind1+ind-1:)
                else
                  readdir=directory_snap(:ind+3)//trim(itoa(iip))
                endif
              else
                cycle iloop
              endif
            endif
            open (lun_input, FILE=trim(readdir)//'/'//trim(file), FORM='unformatted', status='old', iostat=ios)
            if (ios==0) then

              if (islink(trim(readdir)//'/'//trim(file))) &
                 snaplink=trim(readdir)//'/'//trim(file)
!      if (ip<=8) print *, 'read_snap: open, mx,my,mz,nv=', mx, my, mz, nv

              if (lwrite_2d) then

                if (nx == 1) then
                  if (ivar_omit(1)>0) allocate(tmp_omit(1,my,mz,ivar_omit(1):ivar_omit(2)))
                  if (ivar_omit(1)==1) then
                    read (lun_input,iostat=iosr) tmp_omit,a(nghost+1,:,:,:)
                  elseif (ivar_omit(1)>1) then
                    read (lun_input,iostat=iosr) a(nghost+1,:,:,:ivar_omit(1)-1),tmp_omit,a(nghost+1,:,:,ivar_omit(1):)
                  else
                    read (lun_input,iostat=iosr) a(nghost+1,:,:,:)
                  endif
                elseif (ny == 1) then
                  if (ivar_omit(1)>0) allocate(tmp_omit(mx,1,mz,ivar_omit(1):ivar_omit(2)))
                  if (ivar_omit(1)==1) then
                    read (lun_input,iostat=iosr) tmp_omit,a(:,nghost+1,:,:)
                  elseif (ivar_omit(1)>1) then
                    read (lun_input,iostat=iosr) a(:,nghost+1,:,:ivar_omit(1)-1),tmp_omit,a(:,nghost+1,:,ivar_omit(1):)
                  else
                    read (lun_input,iostat=iosr) a(:,nghost+1,:,:)
                  endif
                elseif (nz == 1) then
                  if (ivar_omit(1)>0) allocate(tmp_omit(mx,my,1,ivar_omit(1):ivar_omit(2)))
                  if (ivar_omit(1)==1) then
                    read (lun_input,iostat=iosr) tmp_omit,a(:,:,nghost+1,:)
                  elseif (ivar_omit(1)>1) then
                    read (lun_input,iostat=iosr) a(:,:,nghost+1,:ivar_omit(1)-1),tmp_omit,a(:,:,nghost+1,ivar_omit(1):)
                  else
                    read (lun_input,iostat=iosr) a(:,:,nghost+1,:)
                  endif
                else
                  call fatal_error ('read_snap', 'lwrite_2d used for 3-D simulation!')
                endif

              else
!
                if (nghost_read_fewer==0) then

                  if (lzaver_on_input) then
!
!  Read data from all processors in z direction in source directory and average them.
!
                    allocate(tmp(mx,my,mz_src,nv))
                    if (ivar_omit(1)>0) allocate(tmp_omit(mx,my,mz_src,ivar_omit(1):ivar_omit(2)))
                    iosr=read_part(lun_input,tmp,tmp_omit)
                    a(:,:,nghost+1,:)=sum(tmp(:,:,nghost+1:mz_src-nghost,:),3)
                    do ipz_src=1,nprocz_src-1
                      open(lun_input1,FILE=trim(srcdir)//'/data/proc'//trim(itoa(iproc+ipz_src*nprocxy))// &
                                      '/'//trim(file), FORM='unformatted', status='old', iostat=ios)
                      iosr=read_part(lun_input1,tmp,tmp_omit)
                      a(:,:,nghost+1,:)=a(:,:,nghost+1,:)+sum(tmp(:,:,nghost+1:mz_src-nghost,:),3)
                      close(lun_input1)
                    enddo
                    a(:,:,nghost+1,:)=a(:,:,nghost+1,:)/nzgrid_src
                  else
                    if (ivar_omit(1)>0) allocate(tmp_omit(mx,my,mz,ivar_omit(1):ivar_omit(2)))
                    iosr=read_part(lun_input,a,tmp_omit)
                  endif
!
                elseif (nghost_read_fewer>0) then
!
!  Possibility of reading data with different numbers of ghost zones.
!  In that case, one must regenerate the mesh with luse_oldgrid=T.
!
                  if (ivar_omit(1)>0) allocate(tmp_omit(mx-2*nghost_read_fewer,my-2*nghost_read_fewer, &
                                                        mz-2*nghost_read_fewer,ivar_omit(1):ivar_omit(2)))
                  if (ivar_omit(1)==1) then
                    read (lun_input,iostat=iosr) tmp_omit, &
                         a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                           1+nghost_read_fewer:my-nghost_read_fewer, &
                           1+nghost_read_fewer:mz-nghost_read_fewer, :)
                  elseif (ivar_omit(1)>1) then
                    read (lun_input,iostat=iosr) &
                         a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                           1+nghost_read_fewer:my-nghost_read_fewer, &
                           1+nghost_read_fewer:mz-nghost_read_fewer,:ivar_omit(1)-1),tmp_omit, &
                         a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                           1+nghost_read_fewer:my-nghost_read_fewer, &
                           1+nghost_read_fewer:mz-nghost_read_fewer, ivar_omit(1):)
                  else
                    read (lun_input,iostat=iosr) &
                         a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                           1+nghost_read_fewer:my-nghost_read_fewer, &
                           1+nghost_read_fewer:mz-nghost_read_fewer,:)
                  endif
!
!  The following 3 possibilities allow us to replicate 1-D data input
!  in x (nghost_read_fewer=-1), y (-2), or z (-3) correspondingly.
!
                else
                  len1d=2*nghost+1
                  if (nghost_read_fewer==-1) then
                    if (ivar_omit(1)>0) allocate(tmp_omit(mx,len1d,len1d,ivar_omit(1):ivar_omit(2)))
                    if (ivar_omit(1)==1) then
                      read (lun_input,iostat=iosr) tmp_omit, &
                                       a(:,1:len1d,1:len1d,:)
                    elseif (ivar_omit(1)>1) then
                      read (lun_input,iostat=iosr) a(:,1:len1d,1:len1d,:ivar_omit(1)-1),tmp_omit, &
                                       a(:,1:len1d,1:len1d, ivar_omit(1):)
                    else
                      read (lun_input,iostat=iosr) a(:,1:len1d,1:len1d,:)
                    endif
                    a=spread(spread(a(:,m1,n1,:),2,my),3,mz)
                  elseif (nghost_read_fewer==-2) then
                    if (ivar_omit(1)>0) allocate(tmp_omit(len1d,my,len1d,ivar_omit(1):ivar_omit(2)))
                    if (ivar_omit(1)==1) then
                      read (lun_input,iostat=iosr) tmp_omit, &
                                       a(1:len1d,:,1:len1d,:)
                    elseif (ivar_omit(1)>1) then
                      read (lun_input,iostat=iosr) a(1:len1d,:,1:len1d,:ivar_omit(1)-1),tmp_omit, &
                                       a(1:len1d,:,1:len1d, ivar_omit(1):)
                    else
                      read (lun_input,iostat=iosr) a(1:len1d,:,1:len1d,:)
                    endif
                    a=spread(spread(a(l1,:,n1,:),1,mx),3,mz)
                  elseif (nghost_read_fewer==-3) then
                    if (ivar_omit(1)>0) allocate(tmp_omit(len1d,len1d,mz,ivar_omit(1):ivar_omit(2)))
                    if (ivar_omit(1)==1) then
                      read (lun_input,iostat=iosr) tmp_omit, &
                                       a(1:len1d,1:len1d,:,:)
                    elseif (ivar_omit(1)>1) then
                      read (lun_input,iostat=iosr) a(1:len1d,1:len1d,:,:ivar_omit(1)-1),tmp_omit, &
                                       a(1:len1d,1:len1d,:, ivar_omit(1):)
                    else
                      read (lun_input,iostat=iosr) a(1:len1d,1:len1d,:,:)
                    endif
                    a=spread(spread(a(l1,m1,:,:),1,mx),2,my)
                  else
                    call fatal_error('read_snap','nghost_read_fewer must be >=0')
                  endif
                endif
              endif
              if (iosr/=0) then
                close(lun_input)
                l0=.false.
                cycle iloop
              endif
              if (ip <= 8) print *, 'read_snap: read ', file
              if (mode == 1) then
!
!  Check whether we want to read deltay from snapshot.
!
                if (lshear.and..not.lread_oldsnap_noshear) then
                  read (lun_input,iostat=iosr) t_sp, x, y, z, dx, dy, dz, deltay
                else
                  if (nghost_read_fewer==0) then
                    read (lun_input,iostat=iosr) t_sp, x, y, z, dx, dy, dz
                  elseif (nghost_read_fewer>0) then
                    read (lun_input,iostat=iosr) t_sp
                  endif
                endif
!
!  Verify the read values for x, y and z.
!
                if (ip <= 3) print *, 'read_snap: x=', x
                if (ip <= 3) print *, 'read_snap: y=', y
                if (ip <= 3) print *, 'read_snap: z=', z
              else
                read (lun_input,iostat=iosr) t_sp
              endif
              if (iosr/=0) then
                close(lun_input)
                l0=.false.
                cycle iloop
              endif
              lok=.true.
              if (.not.l0) then
                message=' reading '//trim(readdir)//'/'//trim(file)// &
                        ' instead of '//trim(directory_snap)//'/'//trim(file)
                call warning('read_snap', trim(message), iproc)
                open(11,file=trim(datadir)//'/jobid.dat',position='append')
                write(11,*) trim(message) 
                close(11)
                if (mailaddress/='') lmail=.true.
              endif
              exit kloop
            endif
            l0=.false.
          enddo iloop
        enddo jloop
      enddo kloop

      if (.not.lok) &
        call fatal_error('read_snap','neither '//trim(directory_snap)//'/'//trim(file)// &
                         ' nor any rescue file could be opened ', .true.)
!
!  Verify consistency of the snapshots regarding their timestamp,
!  unless ireset_tstart=T, in which case we reset all times to tstart.
!
      if ((ireset_tstart == 0) .or. (tstart == impossible)) then
!
        t_test = t_sp
        call mpibcast_real(t_test,comm=MPI_COMM_PENCIL)
        call mpiallreduce_or((t_test /= t_sp) .and. .not. lread_from_other_prec &
                             .or. (abs(t_test-t_sp) > 1.e-6),ltest, MPI_COMM_PENCIL)
!
!  If timestamp deviates at any processor
!
        if (ltest) then
          if (ireset_tstart > 0) then
!
!  If reset of tstart enabled and tstart unspecified, use minimum of all t_sp
!
            if (ireset_tstart == MINT) then
              call mpiallreduce_min(t_sp,t_red,MPI_COMM_PENCIL)
              if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                     ' Using (min) t=', t_red,'with ireset_tstart=', MINT,'.'
            elseif (ireset_tstart >= MAXT) then
              call mpiallreduce_max(t_sp,t_red,MPI_COMM_PENCIL)
              if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                     ' Using (max) t=', t_red,'with ireset_tstart=', MAXT,'.'
            endif
            tstart = t_red
            if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT. Using t=', tstart, '.'
          else
            write (*,*) 'ERROR: '//trim(directory_snap)//'/'//trim(file)// &
                        ' IS INCONSISTENT: t=', t_sp
            call stop_it('read_snap')
          endif
        else
          tstart = t_sp
        endif
!
      endif
!
!  Set time or overwrite it by a given value.
!
      if (ireset_tstart > 0) then
        t = tstart
      else
        t = t_sp
      endif
!
!  Verify the read value for t.
!
      if (ip <= 3) print *, 'read_snap: t=', t
!
      if (lrepair_snap.or.snaplink/='') then
!
!  Send "repair" warning to user.
!
        call mpiallreduce_or(lmail,lmail_)
        if (lmail_) then
          if (lroot) then
            mailstr=trim(message)
            do ipp=1,ncpus-1
              call mpirecv_char(message,ipp,ipp)
              if (message/='') mailstr=trim(mailstr)//trim(message)
            enddo
          else
            call mpisend_char(message,root,iproc)
          endif
          if (lroot) then
            if (file_exists(trim(datadir)//'/jobid.dat')) then
              open(1,file=trim(datadir)//'/jobid.dat',position='append')
              backspace(1)
              read(1,*) cjobid
              close(1) 
              mailstr=', job '//trim(cjobid)//'. WARNING - read_snap:'//trim(mailstr)
            else
              mailstr='. WARNING - read_snap:'//trim(mailstr)
            endif
          
            call system_cmd( &
              'echo From `pwd`'//trim(mailstr)//"|"//trim(mailcmd)//" -s 'PencilCode Message' "//trim(mailaddress)//' >& /dev/null')
          endif
        endif
      endif

contains
!-------------------------------------------------------------------------------------------------------------------
      subroutine prepare_zaver_on_input
!
!  Prepares "z-averaging input" from a snapshot in a different (the source) run directory, 
!  which is inferred from the contents of the link that readdir//file needs to be.
!  Only usable for spherical coordinates by now.
!
!  31-mar-21/MR: coded
!
        character(LEN=20) :: cres
        character(LEN=2*fnlen) :: link

        if (.not.islink((trim(readdir)//'/'//trim(file)))) then
          call fatal_error('prepare_zaver_on_input',trim(readdir)//'/'//trim(file)// &
                                                    ' must be symbolic link for lzaver_on_input=T')
        else
          if (lactive_dimension(3)) &
            call fatal_error('prepare_zaver_on_input','z must be degenerate for lzaver_on_input=T')

          if (readlink(trim(readdir)//'/'//trim(file),link)) then
            srcdir=link(1:index(trim(link),'data')-1)
          else
            call fatal_error('prepare_zaver_on_input','file pointed to by symbolic link '//trim(readdir)//'/'//trim(file)// &
                             ' does not exist or string "link" too short')
          endif
!
!  Obtain nproc[xyz], n[xyz]grid from source run directory and check compatibility with target run.
!
          call extract_str("grep '^[^!]*nprocx *= *[1-9][0-9]*' "//trim(srcdir)// &
                           "/src/cparam.local | sed -e's/^.*nprocx *= *\([1-9][0-9]*\).*$/\1/'",cres)
          if (nprocx/=atoi(cres)) call fatal_error('prepare_zaver_on_input','non-matching nprocx in '//trim(srcdir))
          call extract_str("grep '^[^!]*nprocy *= *[1-9][0-9]*' "//trim(srcdir)// &
                           "/src/cparam.local | sed -e's/^.*nprocy *= *\([1-9][0-9]*\).*$/\1/'",cres)
          if (nprocy/=atoi(cres)) call fatal_error('prepare_zaver_on_input','non-matching nprocy in '//trim(srcdir))

          call extract_str("grep '^[^!]*nxgrid *= *[1-9][0-9]*' "//trim(srcdir)// &
                           "/src/cparam.local | sed -e's/^.*nxgrid *= *\([1-9][0-9]*\).*$/\1/'",cres)
          if (nxgrid/=atoi(cres)) call fatal_error('prepare_zaver_on_input','non-matching nxgrid in '//trim(srcdir))
          call extract_str("grep '^[^!]*nygrid *= *[1-9][0-9]*' "//trim(srcdir)// &
                           "/src/cparam.local | sed -e's/^.*nygrid *= *\([1-9][0-9]*\).*$/\1/'",cres)
          if (nygrid/=atoi(cres)) call fatal_error('prepare_zaver_on_input','non-matching nygrid in '//trim(srcdir))

          call extract_str("grep '^[^!]*nprocz *= *[1-9][0-9]*' "//trim(srcdir)// &
                           "/src/cparam.local | sed -e's/^.*nprocz *= *\([1-9][0-9]*\).*$/\1/'",cres)
          nprocz_src=atoi(cres)
          call extract_str("grep '^[^!]*nzgrid *= *[1-9][0-9]*' "//trim(srcdir)// &
                           "/src/cparam.local | sed -e's/^.*nzgrid *= *\([1-9][0-9]*\).*$/\1/'",cres)
          nzgrid_src=atoi(cres)
          mz_src=nzgrid_src/nprocz_src+2*nghost

        endif

      endsubroutine prepare_zaver_on_input
!-------------------------------------------------------------------------------------------------------------------
