!
! Code fragment to be used in both read_snap_single and read_snap_double - type agnostic
!
      real :: t_test   ! t in single precision for backwards compatibility
      logical :: ltest,lok,l0
      integer :: len1d,ios,iosr,ii,iia,jj,kk,iip,ind,ind1
      character(LEN=fnlen) :: readdir
!
      if (lserial_io) call start_serialize

      readdir = directory_snap
      lok=.false.; l0=.true.
!
kloop:do kk=max(0,ipz-1),min(ipz+1,nprocz-1)
jloop:  do jj=max(0,ipy-1),min(ipy+1,nprocy-1)
          iia=max(0,ipx-1)
          if (l0) iia=iia-1             ! add one iteration for the original path
iloop:    do ii=iia,min(nprocx-1,ipx+1)
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
                  if (ivar_omit1>0) allocate(tmp(1,my,mz,ivar_omit1:ivar_omit2))
                  if (ivar_omit1==1) then
                    read (lun_input,iostat=iosr) tmp,a(nghost+1,:,:,:)
                  elseif (ivar_omit1>1) then
                    read (lun_input,iostat=iosr) a(nghost+1,:,:,:ivar_omit1-1),tmp,a(nghost+1,:,:,ivar_omit1:)
                  else
                    read (lun_input,iostat=iosr) a(nghost+1,:,:,:)
                  endif
                elseif (ny == 1) then
                  if (ivar_omit1>0) allocate(tmp(mx,1,mz,ivar_omit1:ivar_omit2))
                  if (ivar_omit1==1) then
                    read (lun_input,iostat=iosr) tmp,a(:,nghost+1,:,:)
                  elseif (ivar_omit1>1) then
                    read (lun_input,iostat=iosr) a(:,nghost+1,:,:ivar_omit1-1),tmp,a(:,nghost+1,:,ivar_omit1:)
                  else
                    read (lun_input,iostat=iosr) a(:,nghost+1,:,:)
                  endif
                elseif (nz == 1) then
                  if (ivar_omit1>0) allocate(tmp(mx,my,1,ivar_omit1:ivar_omit2))
                  if (ivar_omit1==1) then
                    read (lun_input,iostat=iosr) tmp,a(:,:,nghost+1,:)
                  elseif (ivar_omit1>1) then
                    read (lun_input,iostat=iosr) a(:,:,nghost+1,:ivar_omit1-1),tmp,a(:,:,nghost+1,ivar_omit1:)
                  else
                    read (lun_input,iostat=iosr) a(:,:,nghost+1,:)
                  endif
                else
                  call fatal_error ('read_snap', 'lwrite_2d used for 3-D simulation!')
                endif
              else
!
                if (ivar_omit1>0) allocate(tmp(mx,my,mz,ivar_omit1:ivar_omit2))
!
!  Possibility of reading data with different numbers of ghost zones.
!  In that case, one must regenerate the mesh with luse_oldgrid=T.
!
                if (nghost_read_fewer==0) then
                  if (ivar_omit1==1) then
                    read (lun_input,iostat=iosr) tmp,a
                  elseif (ivar_omit1>1) then
                    read (lun_input,iostat=iosr) a(:,:,:,:ivar_omit1-1),tmp,a(:,:,:,ivar_omit1:)
                  else
                    read (lun_input,iostat=iosr) a
                  endif
                elseif (nghost_read_fewer>0) then
                  if (ivar_omit1==1) then
                    read (lun_input,iostat=iosr) tmp, &
                         a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                           1+nghost_read_fewer:my-nghost_read_fewer, &
                           1+nghost_read_fewer:mz-nghost_read_fewer, :)
                  elseif (ivar_omit1>1) then
                    read (lun_input,iostat=iosr) &
                         a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                           1+nghost_read_fewer:my-nghost_read_fewer, &
                           1+nghost_read_fewer:mz-nghost_read_fewer,:ivar_omit1-1),tmp, &
                         a(1+nghost_read_fewer:mx-nghost_read_fewer, &
                           1+nghost_read_fewer:my-nghost_read_fewer, &
                           1+nghost_read_fewer:mz-nghost_read_fewer, ivar_omit1:)
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
                    if (ivar_omit1==1) then
                      read (lun_input,iostat=iosr) tmp, &
                                       a(:,1:len1d,1:len1d,:)
                    elseif (ivar_omit1>1) then
                      read (lun_input,iostat=iosr) a(:,1:len1d,1:len1d,:ivar_omit1-1),tmp, &
                                       a(:,1:len1d,1:len1d, ivar_omit1:)
                    else
                      read (lun_input,iostat=iosr) a(:,1:len1d,1:len1d,:)
                    endif
                    a=spread(spread(a(:,m1,n1,:),2,my),3,mz)
                  elseif (nghost_read_fewer==-2) then
                    if (ivar_omit1==1) then
                      read (lun_input,iostat=iosr) tmp, &
                                       a(1:len1d,:,1:len1d,:)
                    elseif (ivar_omit1>1) then
                      read (lun_input,iostat=iosr) a(1:len1d,:,1:len1d,:ivar_omit1-1),tmp, &
                                       a(1:len1d,:,1:len1d, ivar_omit1:)
                    else
                      read (lun_input,iostat=iosr) a(1:len1d,:,1:len1d,:)
                    endif
                    a=spread(spread(a(l1,:,n1,:),1,mx),3,mz)
                  elseif (nghost_read_fewer==-3) then
                    if (ivar_omit1==1) then
                      read (lun_input,iostat=iosr) tmp, &
                                       a(1:len1d,1:len1d,:,:)
                    elseif (ivar_omit1>1) then
                      read (lun_input,iostat=iosr) a(1:len1d,1:len1d,:,:ivar_omit1-1),tmp, &
                                       a(1:len1d,1:len1d,:, ivar_omit1:)
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
              if (.not.l0) call warning('read_snap', 'reading '//trim(readdir)//'/'//trim(file)// &
                                        ' instead of '//trim(directory_snap)//'/'//trim(file),iproc)
              exit kloop
            endif
            l0=.false.
          enddo iloop
        enddo jloop
      enddo kloop

      if (.not.lok) &
        call fatal_error('read_snap','neither '//trim(directory_snap)//'/'//trim(file)// &
                         ' nor any rescue file could be opened ')
!
!  Verify consistency of the snapshots regarding their timestamp,
!  unless ireset_tstart=T, in which case we reset all times to tstart.
!
      if ((ireset_tstart == 0) .or. (tstart == impossible)) then
!
        t_test = t_sp
        call mpibcast_real(t_test,comm=MPI_COMM_WORLD)
        call mpiallreduce_or((t_test /= t_sp) .and. .not. lread_from_other_prec &
                             .or. (abs(t_test-t_sp) > 1.e-6),ltest, MPI_COMM_WORLD)
!
!  If timestamp deviates at any processor
!
        if (ltest) then
          if (ireset_tstart > 0) then
!
!  If reset of tstart enabled and tstart unspecified, use minimum of all t_sp
!
            if (ireset_tstart == MINT) then
              call mpiallreduce_min(t_sp,t_red,MPI_COMM_WORLD)
              if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                     ' Using (min) t=', t_red,'with ireset_tstart=', MINT,'.'
            elseif (ireset_tstart >= MAXT) then
              call mpiallreduce_max(t_sp,t_red,MPI_COMM_WORLD)
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
