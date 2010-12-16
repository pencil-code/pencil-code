! $Id$
!
!***********************************************************************
      program rvid_box
!
!  read and combine slices from individual processor directories data
!  /procN, write them to data/, where the can be used by used by
!  rvid_box.pro.Adapated from read_videosslices.f90 but reads
!  all files without wating for a user input.
!
!  14-dez-10/bing: coded 
!
      use Cparam
      use General
!
      implicit none
!
!
      real, dimension (nxgrid,nygrid) :: xy,xy2,xy3,xy4
      real, dimension (nxgrid,nzgrid) :: xz
      real, dimension (nygrid,nzgrid) :: yz
!
      real, dimension (nx,ny) :: xy_loc,xy2_loc,xy3_loc,xy4_loc
      real, dimension (nx,nz) :: xz_loc
      real, dimension (ny,nz) :: yz_loc
!
      integer :: ipx,ipy,ipz,iproc,it,stat
      integer :: ipy1,ipx1,ipz1,ipz2,ipz3,ipz4
      integer :: lun_r1=11,lun_r2=12,lun_r3=13,lun_r4=14,lun_r5=15,lun_r6=16
      integer :: lun_w1=21,lun_w2=22,lun_w3=23,lun_w4=24,lun_w5=25,lun_w6=26
      integer :: lun_pos=1,lun_video=2
      integer :: nevery=1,sindex
      integer :: isep1=0,idummy
      logical :: eof=.false.
      logical :: err=.false.
      logical :: lread_slice_xy,lread_slice_xy2,lread_slice_xy3
      logical :: lread_slice_xy4,lread_slice_xz,lread_slice_yz
      real :: t
      real :: slice_xpos=0., slice_ypos=0., slice_zpos=0., slice_z2pos=0.
      real :: slice_z3pos=0., slice_z4pos=0.
!
      character (len=120) :: file='',fullname='',wfile='',directory=''
      character (len=120) :: datadir='data',path=''
      character (len=5) :: chproc='',nindex=''
      character (len=20) :: field='lnrho',field2='',cfield=''
!
      logical :: exists,lwrite=.true.
!
      real :: min_xy_loc,min_xy2_loc,min_xz_loc,min_yz_loc
      real :: max_xy_loc,max_xy2_loc,max_xz_loc,max_yz_loc
      real :: min_xy3_loc,min_xy4_loc
      real :: max_xy3_loc,max_xy4_loc
!
!  initialize minimum and maximum values for each plane
!
      min_xy_loc=huge(min_xy_loc); max_xy_loc=-huge(max_xy_loc)
      min_xy2_loc=huge(min_xy2_loc); max_xy2_loc=-huge(max_xy2_loc)
      min_xz_loc=huge(min_xz_loc); max_xz_loc=-huge(max_xz_loc)
      min_yz_loc=huge(min_yz_loc); max_yz_loc=-huge(max_yz_loc)
      min_xy3_loc=huge(min_xy3_loc); max_xy3_loc=-huge(max_xy3_loc)
      min_xy4_loc=huge(min_xy4_loc); max_xy4_loc=-huge(max_xy4_loc)
!
!  read name of the field from video.in
!
      inquire(file='video.in',exist=exists)
      if (exists) then
        open(lun_video,file='video.in')
      else
        print*,'ERROR: video.in not found'
        STOP 1
      endif
!
!
!  Loop over aller processors to find the positions of the slices.
!  Therefore read all procN/slice_postions.dat
!
        ipz1=-1; ipz2=-1; ipz3=-1
        ipz4=-1; ipy1=-1; ipx1=-1
!
        do ipx=0,nprocx-1
          do ipy=0,nprocy-1
            do ipz=0,nprocz-1
              iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
              call chn(iproc,chproc,'directory_names')
              call safe_character_assign(directory, trim(datadir)//'/proc'//chproc)
!  Check for existence first 
              inquire(FILE=trim(directory)//'/slice_position.dat',EXIST=exists)
              if (.not.exists) then
                print*,'slice_position.dat for iproc=',iproc,'not found!'
                STOP 1
              endif
              open(lun_pos,file=trim(directory)//'/slice_position.dat',STATUS='unknown')
              read(lun_pos,'(l5,i5)') lread_slice_xy,idummy
              read(lun_pos,'(l5,i5)') lread_slice_xy2,idummy
              read(lun_pos,'(l5,i5)') lread_slice_xy3,idummy
              read(lun_pos,'(l5,i5)') lread_slice_xy4,idummy
              read(lun_pos,'(l5,i5)') lread_slice_xz,idummy
              read(lun_pos,'(l5,i5)') lread_slice_yz,idummy
              close(lun_pos)            
              if (lread_slice_xy) ipz1=ipz
              if (lread_slice_xy2) ipz2=ipz
              if (lread_slice_xy3) ipz3=ipz
              if (lread_slice_xy4) ipz4=ipz
              if (lread_slice_xz) ipy1=ipy
              if (lread_slice_yz) ipx1=ipx
            enddo
          enddo
        enddo
!
      stat=0
      sindex=1
!
      do while (stat==0)
        if (sindex==1) then
          read(lun_video,*,iostat=stat) cfield
!
!  read stride from internal reader!
!
          isep1=index(cfield,' ')
          field2=cfield(1:isep1)
          field=''
        endif
!
!  In case of a vector slice loop over the indices.
        if (field2=='uu'.or.field2=='oo'.or.field2=='uud'.or. &
            field2=='jj'.or.field2=='bb'.or.field2=='poynting'.or. &
            field2=='aa'.or.field2=='Frad'.or.field2=='bb1'.or. &
            field2=='bb11'.or.field2=='uu11') then
          call chn(sindex,nindex)
          call safe_character_assign(field,trim(field2)//trim(nindex))
!
          if (sindex<3) then
            sindex=sindex+1
          else
            sindex=1
          endif
        else
          field=field2
        endif
        write(*,*) "Reading next: ",field
!
!  loop through all times
!  reset error to false at each time step
!
        it = 1
        eof=.false.
        do while (.not.eof)
!
!  Check whether this is a time where we want to write data,
!  or whether we want to skip it.
!
          lwrite=(mod(it,nevery)==0)
!
!  First xy-plane:
!
          if (ipz1/=-1) then
            ipz=ipz1
            do ipy=0,nprocy-1
              do ipx=0,nprocx-1
                iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
                call chn(iproc,chproc,'rvid_box: xy-plane')
                call safe_character_assign(path,trim(datadir)//'/proc'//chproc)
                call safe_character_assign(file,'/slice_'//trim(field)//'.xy')
                call safe_character_assign(fullname,trim(path)//trim(file))
                inquire(FILE=trim(fullname),EXIST=exists)
                if (.not.exists) then
                  print*,"Slice not found: ", fullname
                  xy(:,1+ipy*ny:ny+ipy*ny)=0.
                  goto 998
                endif
                call read_slice(trim(fullname),xy_loc,slice_zpos,nx,ny,t,it,lun_r1,eof,err)
                min_xy_loc=min(min_xy_loc,minval(xy_loc))
                max_xy_loc=max(max_xy_loc,maxval(xy_loc))
                xy(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy_loc
              enddo
            enddo
            call safe_character_assign(wfile,trim(datadir)//trim(file))
            call append_slice(trim(wfile),xy,slice_zpos,nxgrid,nygrid,t,it,lun_w1,lwrite)
          endif
!
!  Continui only of not leof
!
          if (.not.eof) then
!  Second xy-plane:
!
            if (ipz2/=-1) then
              ipz=ipz2
              do ipy=0,nprocy-1
                do ipx=0,nprocx-1
                  iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
                  call chn(iproc,chproc,'rvid_box: xy2 plane')
                  call safe_character_assign(path,trim(datadir)//'/proc'//chproc)
                  call safe_character_assign(file,'/slice_'//trim(field)//'.xy2')
                  call safe_character_assign(fullname,trim(path)//trim(file))
                  inquire(FILE=trim(fullname),EXIST=exists)
                  if (.not.exists) then
                    print*,"Slice not found: ", fullname
                    xy2(:,1+ipy*ny:ny+ipy*ny)=0.
                    goto 998
                  endif
                  call read_slice(trim(fullname),xy2_loc,slice_z2pos,nx,ny,t,it,lun_r2,eof,err)                
                  min_xy2_loc=min(min_xy2_loc,minval(xy2_loc))
                  max_xy2_loc=max(max_xy2_loc,maxval(xy2_loc))
                  xy2(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy2_loc
                enddo
              enddo
              call safe_character_assign(wfile,trim(datadir)//trim(file))
              call append_slice(trim(wfile),xy2,slice_z2pos,nxgrid,nygrid,t,it,lun_w2,lwrite)
            endif
!
!  Third xy-plane:
!
            if (ipz3/=-1) then
              ipz=ipz3
              do ipy=0,nprocy-1
                do ipx=0,nprocx-1
                  iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
                  call chn(iproc,chproc,'rvid_box: xy3-plane')
                  call safe_character_assign(path,trim(datadir)//'/proc'//chproc)
                  call safe_character_assign(file,'/slice_'//trim(field)//'.xy3')
                  call safe_character_assign(fullname,trim(path)//trim(file))
                  inquire(FILE=trim(fullname),EXIST=exists)
                  if (.not.exists) then
                    print*,"Slice not found: ", fullname
                    xy3(:,1+ipy*ny:ny+ipy*ny)=0.
                    goto 998
                  endif
                  call read_slice(trim(fullname),xy3_loc,slice_z3pos,nx,ny,t,it,lun_r3,eof,err)
                  min_xy3_loc=min(min_xy3_loc,minval(xy3_loc))
                  max_xy3_loc=max(max_xy3_loc,maxval(xy3_loc))
                  xy3(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy3_loc
                enddo
              enddo
              call safe_character_assign(wfile,trim(datadir)//trim(file))
              call append_slice(trim(wfile),xy3,slice_z3pos,nxgrid,nygrid,t,it,lun_w3,lwrite)
            endif
!
!  Fourth xy-plane:
!
            if (ipz4/=-1) then
              ipz=ipz4
              do ipy=0,nprocy-1
                do ipx=0,nprocx-1
                  iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
                  call chn(iproc,chproc,'rvid_box: xy4-plane')
                  call safe_character_assign(path,trim(datadir)//'/proc'//chproc)
                  call safe_character_assign(file,'/slice_'//trim(field)//'.xy4')
                  call safe_character_assign(fullname,trim(path)//trim(file))
                  inquire(FILE=trim(fullname),EXIST=exists)
                  if (.not.exists) then
                    print*,"Slice not found: ", fullname
                    xy4(:,1+ipy*ny:ny+ipy*ny)=0.
                    goto 998
                  endif
                  call read_slice(trim(fullname),xy4_loc,slice_z4pos,nx,ny,t,it,lun_r4,eof,err)
                  min_xy4_loc=min(min_xy4_loc,minval(xy4_loc))
                  max_xy4_loc=max(max_xy4_loc,maxval(xy4_loc))
                  xy4(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy4_loc
                enddo
              enddo
              call safe_character_assign(wfile,trim(datadir)//trim(file))
              call append_slice(trim(wfile),xy4,slice_z4pos,nxgrid,nygrid,t,it,lun_w4,lwrite)
            endif
!
!  First xz-plane:
!
            if (ipy1/=-1) then
              ipy=ipy1
              do ipz=0,nprocz-1
                do ipx=0,nprocx-1
                  iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
                  call chn(iproc,chproc,'rvid_box: xz-plane')
                  call safe_character_assign(path,trim(datadir)//'/proc'//chproc)
                  call safe_character_assign(file,'/slice_'//trim(field)//'.xz')
                  call safe_character_assign(fullname,trim(path)//trim(file))
                  inquire(FILE=trim(fullname),EXIST=exists)
                  if (.not.exists) then
                    print*,"Slice not found: ", fullname
                    xz(:,1+ipz*nz:nz+ipz*nz)=0.
                    goto 998
                  endif
                  call read_slice(trim(fullname),xz_loc,slice_ypos,nx,nz,t,it,lun_r5,eof,err)
                  min_xz_loc=min(min_xz_loc,minval(xz_loc))
                  max_xz_loc=max(max_xz_loc,maxval(xz_loc))
                  xz(1+ipx*nx:nx+ipx*nx,1+ipz*nz:nz+ipz*nz)=xz_loc
                enddo
              enddo
              call safe_character_assign(wfile,trim(datadir)//trim(file))
              call append_slice(trim(wfile),xz,slice_ypos,nxgrid,nzgrid,t,it,lun_w5,lwrite)
            endif
!
!  First yz-plane:
!
            if (ipx1/=-1) then 
              ipx=ipx1
              do ipz=0,nprocz-1
                do ipy=0,nprocy-1
                  iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
                  call chn(iproc,chproc,'rvid_box: yz-plane')
                  call safe_character_assign(path,trim(datadir)//'/proc'//chproc)
                  call safe_character_assign(file,'/slice_'//trim(field)//'.yz')
                  call safe_character_assign(fullname,trim(path)//trim(file))
                  inquire(FILE=trim(fullname),EXIST=exists)
                  if (.not.exists) then
                    print*,"Slice not found: ", fullname
                    yz(1+ipy*ny:ny+ipy*ny,1+ipz*nz:nz+ipz*nz)=0.
                    goto 998
                  endif
                  call read_slice(trim(fullname),yz_loc,slice_xpos,ny,nz,t,it,lun_r6,eof,err)
                  min_yz_loc=min(min_yz_loc,minval(yz_loc))
                  max_yz_loc=max(max_yz_loc,maxval(yz_loc))
                  if (eof) goto 998
                  yz(1+ipy*ny:ny+ipy*ny,1+ipz*nz:nz+ipz*nz)=yz_loc
                enddo
              enddo
              call safe_character_assign(wfile,trim(datadir)//trim(file))
              call append_slice(trim(wfile),yz,slice_xpos,nygrid,nzgrid,t,it,lun_w6,lwrite)
            endif
!
!  We didnt reached eof 
            it = it+1
          else
! we reached and of time series. Close all open files
            close(lun_r1)
            close(lun_r2)
            close(lun_r3)
            close(lun_r4)
            close(lun_r5)
            close(lun_r6)            
            close(lun_w1)
            close(lun_w2)
            close(lun_w3)
            close(lun_w4)
            close(lun_w5)
            close(lun_w6)            
          endif

        enddo
!
998     continue
!
      enddo
!
      close(lun_video)
!
    endprogram rvid_box
!***********************************************************************
    subroutine read_slice(file,a,pos,ndim1,ndim2,t,it,lun,eof,err)
!
! read an existing slice file
!
!  12-nov-02/axel: coded
!
      integer :: ndim1,ndim2
      character (len=*) :: file
      real, dimension (ndim1,ndim2) :: a
      integer :: it,lun
      logical :: eof,err
      real :: t,pos
!
      if (it==1) open(lun,file=file,status='old',form='unformatted')
!
      pos=0.  ! By default (i.e. if missing from record)
      read(lun,end=999,err=998) a,t,pos
      goto 900
!
!  error: suspect wrong record length
!
998   read(lun,end=999,err=997) a,t
      goto 900
!
!  still an error, avoid this time
!
997   err=.true.
      goto 900
!
!  when end of file
!
999   eof=.true.
!
900   continue
!
    endsubroutine read_slice
!***********************************************************************
    subroutine append_slice(file,a,pos,ndim1,ndim2,t,it,lun,lwrite)
!
!  append to a slice file
!
!  12-nov-02/axel: coded
!
      integer :: ndim1,ndim2
      character (len=*) :: file
      real, dimension (ndim1,ndim2) :: a
      logical :: lwrite
      integer :: it,lun
      real :: t, pos
!
      if (it==1) open(lun,file=file,form='unformatted')
      if (lwrite) write(lun) a,t,pos
!
    endsubroutine append_slice
!***********************************************************************
