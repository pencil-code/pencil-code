! $Id$
!
!***********************************************************************
      program rvid_box
!
!  read and combine slices from individual processor directories data
!  /procN, write them to data/, where they can be used by rvid_box.pro
!
!  13-nov-02/axel: coded
!  22-sep-07/axel: changed Xy to xy2, to be compatible with Mac
!
      use Cparam
      use General
!
      implicit none
!
      real, dimension (nxgrid,nygrid) :: xy,xy2,xy3,xy4
      real, dimension (nxgrid,nzgrid) :: xz
      real, dimension (nygrid,nzgrid) :: yz
!
      real, dimension (nx,ny) :: xy_loc,xy2_loc,xy3_loc,xy4_loc
      real, dimension (nx,nz) :: xz_loc
      real, dimension (ny,nz) :: yz_loc
!
      integer :: ipx,ipy,ipz,iproc,it,nt=999999
      integer :: ipy1,ipx1,ipz1,ipz2,ipz3,ipz4
      integer :: lun,lun1=1,lun2=2,lun3=3,lun4=4,lun5=5,lun6=6
      integer :: itdebug=1,nevery=1,nproc_tot
      integer :: isep1=0,isep2=0,idummy
      logical :: eof=.false.
      logical :: err=.false.,err_timestep=.false.
      logical :: lread_slice_xy,lread_slice_xy2,lread_slice_xy3
      logical :: lread_slice_xy4,lread_slice_xz,lread_slice_yz
      real :: t
      real :: slice_xpos=0., slice_ypos=0., slice_zpos=0., slice_z2pos=0.
      real :: slice_z3pos=0., slice_z4pos=0.
!
      character (len=fnlen) :: file='',fullname='',wfile='',directory=''
      character (len=fnlen) :: datadir='data',path='',cfield=''
      character (len=20) :: field='lnrho'
!
      logical :: exists, lwritten_something=.false.,lwrite=.true.
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
!  read name of the field (must coincide with file extension)
!
      write(*,'(a)',ADVANCE='NO') 'enter variable (lnrho, uu1, ..., bb3) and stride (e.g. 10): '
      read(*,'(a)') cfield
!
!  read stride from internal reader
!
      isep1=index(cfield,' ')
      isep2=len(cfield)
      field=cfield(1:isep1)
      if (cfield(isep1+1:isep2)/=' ') read(cfield(isep1:isep2),*) nevery
!
! loop over all processors to find the positions of the slices.
! Therefore read all slice_positions.dat
!
      ipz1=-1; ipz2=-1; ipz3=-1
      ipz4=-1; ipy1=-1; ipx1=-1
!
      do ipx=0,nprocx-1
        do ipy=0,nprocy-1
          do ipz=0,nprocz-1
            iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
            call safe_character_assign(directory, trim(datadir)//'/proc'//itoa(iproc))
! check for existence first
            inquire(FILE=trim(directory)//'/slice_position.dat',EXIST=exists)
            if (.not.exists) then
              print*,'slice_position.dat for iproc=',iproc,'not found!'
              goto 999
            endif
            open(1,file=trim(directory)//'/slice_position.dat',STATUS='unknown')
            read(1,'(l5,i5)') lread_slice_xy,idummy
            read(1,'(l5,i5)') lread_slice_xy2,idummy
            read(1,'(l5,i5)') lread_slice_xy3,idummy
            read(1,'(l5,i5)') lread_slice_xy4,idummy
            read(1,'(l5,i5)') lread_slice_xz,idummy
            read(1,'(l5,i5)') lread_slice_yz,idummy
            close(1)
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
!  Need to reset
      if (ipz1/=-1) lread_slice_xy=.true.
      if (ipz2/=-1) lread_slice_xy2=.true.
      if (ipz3/=-1) lread_slice_xy3=.true.
      if (ipz4/=-1) lread_slice_xy4=.true.
      if (ipy1/=-1) lread_slice_xz=.true.
      if (ipx1/=-1) lread_slice_yz=.true.
!
      nproc_tot=0
      if (lread_slice_xy) nproc_tot=nproc_tot+nprocx*nprocy
      if (lread_slice_xy2) nproc_tot=nproc_tot+nprocx*nprocy
      if (lread_slice_xy3) nproc_tot=nproc_tot+nprocx*nprocy
      if (lread_slice_xy4) nproc_tot=nproc_tot+nprocx*nprocy
      if (lread_slice_xz)  nproc_tot=nproc_tot+nprocx*nprocz
      if (lread_slice_yz)  nproc_tot=nproc_tot+nprocy*nprocz
!
      write(*,*) 'Number of file units which will be used:',nproc_tot
      write(*,*) 'In case of restricted number of file units (system dependend)'
      write(*,*) 'use ./src/read_all_videofiles.x instead'
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
            call safe_character_assign(path,trim(datadir)//'/proc'//itoa(iproc))
            call safe_character_assign(file,'/slice_'//trim(field)//'.xy2')
            call safe_character_assign(fullname,trim(path)//trim(file))
            if (it<=itdebug) print*,trim(fullname)
            inquire(FILE=trim(fullname),EXIST=exists)
            if (.not.exists) then
              print*,"Slice not found: ", fullname
              xy2(:,1+ipy*ny:ny+ipy*ny)=0.
              goto 999
            endif
            call read_slice(trim(fullname),xy2_loc,slice_z2pos,nx,ny,t,it,lun,eof,err)
            min_xy2_loc=min(min_xy2_loc,minval(xy2_loc))
            max_xy2_loc=max(max_xy2_loc,maxval(xy2_loc))
            if (eof) goto 999
            xy2(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy2_loc
          enddo
        enddo
        call safe_character_assign(wfile,trim(datadir)//trim(file))
        err_timestep=err
        if (.not.err_timestep) then
          call append_slice(trim(wfile),xy2,slice_z2pos,nxgrid,nygrid,t,it,lun1,lwrite)
          lwritten_something=.true.
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
            call safe_character_assign(path,trim(datadir)//'/proc'//itoa(iproc))
            call safe_character_assign(file,'/slice_'//trim(field)//'.xy')
            call safe_character_assign(fullname,trim(path)//trim(file))
            if (it<=itdebug) print*,trim(fullname)
            inquire(FILE=trim(fullname),EXIST=exists)
            if (.not.exists) then
              print*,"Slice not found: ", fullname
              xy(:,1+ipy*ny:ny+ipy*ny)=0.
              goto 999
            endif
            call read_slice(trim(fullname),xy_loc,slice_zpos,nx,ny,t,it,lun,eof,err)
            min_xy_loc=min(min_xy_loc,minval(xy_loc))
            max_xy_loc=max(max_xy_loc,maxval(xy_loc))
            if (eof) goto 999
            xy(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy_loc
          enddo
        enddo
        call safe_character_assign(wfile,trim(datadir)//trim(file))
        err_timestep=err_timestep.or.err
        if (.not.err_timestep) then
          call append_slice(trim(wfile),xy,slice_zpos,nxgrid,nygrid,t,it,lun2,lwrite)
          lwritten_something=.true.
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
            call safe_character_assign(path,trim(datadir)//'/proc'//itoa(iproc))
            call safe_character_assign(file,'/slice_'//trim(field)//'.xz')
            call safe_character_assign(fullname,trim(path)//trim(file))
            if (it<=itdebug) print*,trim(fullname)
            inquire(FILE=trim(fullname),EXIST=exists)
            if (.not.exists) then
              print*,"Slice not found: ", fullname
              xz(:,1+ipz*nz:nz+ipz*nz)=0.
              goto 999
            endif
            call read_slice(trim(fullname),xz_loc,slice_ypos,nx,nz,t,it,lun,eof,err)
            min_xz_loc=min(min_xz_loc,minval(xz_loc))
            max_xz_loc=max(max_xz_loc,maxval(xz_loc))
            if (eof) goto 999
            xz(1+ipx*nx:nx+ipx*nx,1+ipz*nz:nz+ipz*nz)=xz_loc
          enddo
        enddo
        call safe_character_assign(wfile,trim(datadir)//trim(file))
        err_timestep=err_timestep.or.err
        if (.not.err_timestep) then
          call append_slice(trim(wfile),xz,slice_ypos,nxgrid,nzgrid,t,it,lun3,lwrite)
          lwritten_something=.true.
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
            call safe_character_assign(path,trim(datadir)//'/proc'//itoa(iproc))
            call safe_character_assign(file,'/slice_'//trim(field)//'.yz')
            call safe_character_assign(fullname,trim(path)//trim(file))
            if (it<=itdebug) print*,trim(fullname)
            inquire(FILE=trim(fullname),EXIST=exists)
            if (.not.exists) then
              print*,"Slice not found: ", fullname
              yz(1+ipy*ny:ny+ipy*ny,1+ipz*nz:nz+ipz*nz)=0.
              goto 999
            endif
            call read_slice(trim(fullname),yz_loc,slice_xpos,ny,nz,t,it,lun,eof,err)
            min_yz_loc=min(min_yz_loc,minval(yz_loc))
            max_yz_loc=max(max_yz_loc,maxval(yz_loc))
            if (eof) goto 999
            yz(1+ipy*ny:ny+ipy*ny,1+ipz*nz:nz+ipz*nz)=yz_loc
          enddo
        enddo
        call safe_character_assign(wfile,trim(datadir)//trim(file))
        err_timestep=err_timestep.or.err
        if (.not.err_timestep) then
          call append_slice(trim(wfile),yz,slice_xpos,nygrid,nzgrid,t,it,lun4,lwrite)
          lwritten_something=.true.
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
            call safe_character_assign(path,trim(datadir)//'/proc'//itoa(iproc))
            call safe_character_assign(file,'/slice_'//trim(field)//'.xy3')
            call safe_character_assign(fullname,trim(path)//trim(file))
            if (it<=itdebug) print*,trim(fullname)
            inquire(FILE=trim(fullname),EXIST=exists)
            if (.not.exists) then
              print*,"Slice not found: ", fullname
              xy3(:,1+ipy*ny:ny+ipy*ny)=0.
              goto 999
            endif
            call read_slice(trim(fullname),xy3_loc,slice_z3pos,nx,ny,t,it,lun,eof,err)
            min_xy3_loc=min(min_xy3_loc,minval(xy3_loc))
            max_xy3_loc=max(max_xy3_loc,maxval(xy3_loc))
            if (eof) goto 999
            xy3(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy3_loc
          enddo
        enddo
        call safe_character_assign(wfile,trim(datadir)//trim(file))
        err_timestep=err
        if (.not.err_timestep) then
          call append_slice(trim(wfile),xy3,slice_z3pos,nxgrid,nygrid,t,it,lun5,lwrite)
          lwritten_something=.true.
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
            call safe_character_assign(path,trim(datadir)//'/proc'//itoa(iproc))
            call safe_character_assign(file,'/slice_'//trim(field)//'.xy4')
            call safe_character_assign(fullname,trim(path)//trim(file))
            if (it<=itdebug) print*,trim(fullname)
            inquire(FILE=trim(fullname),EXIST=exists)
            if (.not.exists) then
              print*,"Slice not found: ", fullname
              xy4(:,1+ipy*ny:ny+ipy*ny)=0.
              goto 999
            endif
            call read_slice(trim(fullname),xy4_loc,slice_z4pos,nx,ny,t,it,lun,eof,err)
            min_xy4_loc=min(min_xy4_loc,minval(xy4_loc))
            max_xy4_loc=max(max_xy4_loc,maxval(xy4_loc))
            if (eof) goto 999
            xy4(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny)=xy4_loc
          enddo
        enddo
        call safe_character_assign(wfile,trim(datadir)//trim(file))
        err_timestep=err
        if (.not.err_timestep) then
          call append_slice(trim(wfile),xy4,slice_z4pos,nxgrid,nygrid,t,it,lun6,lwrite)
          lwritten_something=.true.
        else
          print*,'skip writing because of error; t=',t
        endif
      endif
!
!  confirm writing only if lwrite=.true.
!
      if (lwrite) then
        print*,'written full set of slices at t=',t,min_xy_loc,max_xy_loc
      endif
      enddo
!
999   continue
!
      if (lwritten_something) then
        print*,'last file read: ',trim(fullname)
        print*,'-------------------------------------------------'
        print*,'minimum and maximum values:'
        if (ipz1/=-1) print*,'xy-plane:',min_xy_loc,max_xy_loc
        if (ipz2/=-1) print*,'xy2-plane:',min_xy2_loc,max_xy2_loc
        if (ipy1/=-1) print*,'xz-plane:',min_xz_loc,max_xz_loc
        if (ipz3/=-1) print*,'xy3-plane:',min_xy3_loc,max_xy3_loc
        if (ipz4/=-1) print*,'xy4-plane:',min_xy4_loc,max_xy4_loc
        if (ipx1/=-1) print*,'yz-plane:',min_yz_loc,max_yz_loc
        print*,'-------------------------------------------------'
        print*,'finished OK'
      endif
!
      select case (trim(field))
        case ('ux','uy','uz','bx','by','bz','Fradx','Frady','Fradz','ax','ay','az','ox','oy','oz')
          print*,""
          print*,"*****************************************************************************"
          print*,"******                WARNING DEPRECATED SLICE NAME                    ******"
          print*,"*****************************************************************************"
          print*,"*** The slice name '"//trim(field)//"' is deprecated and soon will not be ***"
          print*,"*** supported any longer                                                  ***"
          print*,"*** New slice names are formed by taking the name specified in video.in   ***"
          print*,"*** eg. uu and in the case of vector or other multiple slices appending   ***"
          print*,"*** a number. For example the slice 'ux' is now 'uu1' and the slice 'uz'  ***"
          print*,"*** is now called 'uu3', 'ay'->'aa2' etc. Similarly for aa, bb, oo, uu    ***"
          print*,"*** and Frad slices.                                                      ***"
          print*,"*****************************************************************************"
          print*,""
      endselect
!
    endprogram rvid_box
!***********************************************************************
    subroutine read_ipz_position(file,ipz)
!
!  reads just one number from file
!
!  19-nov-06/axel: coded
!
      character (len=*) :: file
      integer :: ipz,lun=1
!
      open(lun,file=file,status='old',err=98)
      read(lun,*) ipz
      close(lun)
      goto 99
!
!  escape procedure of file doesn't exist
!
98    print*,";;;"
      print*,";;; data/z*_procnum.dat files don't exist."
      print*,";;; Type (e.g. by cut+paste):"
      print*,";;;    cp data/proc*/z*_procnum.dat data"
      print*,";;;"
      stop
!
99    end
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
