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
      integer :: ipx,ipy,ipz,iproc,it,nt=999999,stat
      integer :: ipy1,ipx1,ipz1,ipz2,ipz3,ipz4
      integer :: lun,lun1=1,lun2=2,lun3=3,lun4=4,lun5=5,lun6=6
      integer :: lun_pos=21,lun_video=19
      integer :: nevery=1,sindex
      integer :: isep1=0,isep2=0,idummy
      logical :: eof=.false.
      logical :: err=.false.,err_timestep=.false.
      logical :: lread_slice_xy,lread_slice_xy2,lread_slice_xy3
      logical :: lread_slice_xy4,lread_slice_xz,lread_slice_yz
      real :: t
      real :: slice_xpos=0., slice_ypos=0., slice_zpos=0., slice_z2pos=0.
      real :: slice_z3pos=0., slice_z4pos=0.
!
      character (len=120) :: file='',fullname='',wfile='',directory=''
      character (len=120) :: datadir='data',path='',cfield=''
      character (len=5) :: chproc='',nindex=''
      character (len=20) :: field='lnrho',field2=''
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
        goto 999
      endif

      stat=0
      sindex=1
      do while (stat==0)
        if (sindex==1) then
          read(lun_video,*,iostat=stat) cfield
!
!  read stride from internal reader!
!
          isep1=index(cfield,' ')
          isep2=len(cfield)
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
          if (sindex<3) then; sindex=sindex+1; else; sindex=1; endif
        else
          field=field2
        endif
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
                goto 999
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
!  loop through all times
!  reset error to false at each time step
!
        do it=1,nt
          err_timestep=.false.
          lun=10
!
!  Check whether this is a time where we want to write data,
!  or whether we want to skip it.
!
          lwrite=(mod(it,nevery)==0)
!
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
                call read_slice(trim(fullname),xy2_loc,slice_z2pos,nx,ny,t,it,lun,eof,err)
                
                min_xy2_loc=min(min_xy2_loc,minval(xy2_loc))
                max_xy2_loc=max(max_xy2_loc,maxval(xy2_loc))
                if (eof) goto 998
                xy2(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy2_loc
              enddo
            enddo
            call safe_character_assign(wfile,trim(datadir)//trim(file))
            err_timestep=err
            if (.not.err_timestep) then
              call append_slice(trim(wfile),xy2,slice_z2pos,nxgrid,nygrid,t,it,lun1,lwrite)
            else
              print*,'skip writing because of error; t=',t
            endif
          endif
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
                  call read_slice(trim(fullname),xy_loc,slice_zpos,nx,ny,t,it,lun,eof,err)
                  min_xy_loc=min(min_xy_loc,minval(xy_loc))
                  max_xy_loc=max(max_xy_loc,maxval(xy_loc))
                  if (eof) goto 998
                  xy(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy_loc
                enddo
              enddo
              call safe_character_assign(wfile,trim(datadir)//trim(file))
              err_timestep=err_timestep.or.err
              if (.not.err_timestep) then
                call append_slice(trim(wfile),xy,slice_zpos,nxgrid,nygrid,t,it,lun2,lwrite)
              else
                print*,'skip writing because of error; t=',t
              endif
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
                  call read_slice(trim(fullname),xz_loc,slice_ypos,nx,nz,t,it,lun,eof,err)
                  min_xz_loc=min(min_xz_loc,minval(xz_loc))
                  max_xz_loc=max(max_xz_loc,maxval(xz_loc))
                  if (eof) goto 998
                  xz(1+ipx*nx:nx+ipx*nx,1+ipz*nz:nz+ipz*nz)=xz_loc
                enddo
              enddo
              call safe_character_assign(wfile,trim(datadir)//trim(file))
              err_timestep=err_timestep.or.err
              if (.not.err_timestep) then
                call append_slice(trim(wfile),xz,slice_ypos,nxgrid,nzgrid,t,it,lun3,lwrite)
              else
                print*,'skip writing because of error; t=',t
              endif
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
                  call read_slice(trim(fullname),yz_loc,slice_xpos,ny,nz,t,it,lun,eof,err)
                  min_yz_loc=min(min_yz_loc,minval(yz_loc))
                  max_yz_loc=max(max_yz_loc,maxval(yz_loc))
                  if (eof) goto 998
                  yz(1+ipy*ny:ny+ipy*ny,1+ipz*nz:nz+ipz*nz)=yz_loc
                enddo
              enddo
              call safe_character_assign(wfile,trim(datadir)//trim(file))
              err_timestep=err_timestep.or.err
              if (.not.err_timestep) then
                call append_slice(trim(wfile),yz,slice_xpos,nygrid,nzgrid,t,it,lun4,lwrite)
              else
                print*,'skip writing because of error; t=',t
              endif
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
                  call read_slice(trim(fullname),xy3_loc,slice_z3pos,nx,ny,t,it,lun,eof,err)
                  min_xy3_loc=min(min_xy3_loc,minval(xy3_loc))
                  max_xy3_loc=max(max_xy3_loc,maxval(xy3_loc))
                  if (eof) goto 998
                  xy3(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy3_loc
                enddo
              enddo
              call safe_character_assign(wfile,trim(datadir)//trim(file))
              err_timestep=err_timestep.or.err
              if (.not.err_timestep) then
                call append_slice(trim(wfile),xy3,slice_z3pos,nxgrid,nygrid,t,it,lun5,lwrite)
              else
                print*,'skip writing because of error; t=',t
              endif
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
                  call read_slice(trim(fullname),xy4_loc,slice_z4pos,nx,ny,t,it,lun,eof,err)
                  min_xy4_loc=min(min_xy4_loc,minval(xy4_loc))
                  max_xy4_loc=max(max_xy4_loc,maxval(xy4_loc))
                  if (eof) goto 998
                  xy4(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy4_loc
                enddo
              enddo
              call safe_character_assign(wfile,trim(datadir)//trim(file))
              err_timestep=err_timestep.or.err
              if (.not.err_timestep) then
                call append_slice(trim(wfile),xy4,slice_z4pos,nxgrid,nygrid,t,it,lun6,lwrite)
              else
                print*,'skip writing because of error; t=',t
              endif
            endif
!
!  confirm writing only if lwrite=.true.
!
          enddo
998       continue
          eof=.false.
          err=.false.
          err_timestep=.false.
          close(lun)
          close(lun1)
          close(lun2)
          close(lun3)
          close(lun4)
          close(lun5)
          close(lun6)

        enddo
!
999     continue
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
      lun=lun+1
      goto 900
!
!  error: suspect wrong record length
!
998   read(lun,end=999,err=997) a,t
      lun=lun+1
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
