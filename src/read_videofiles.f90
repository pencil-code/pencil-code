! $Id: read_videofiles.f90,v 1.6 2003-09-02 13:43:12 dobler Exp $

!***********************************************************************
      program rvid_box
!
!  read slices that are necessary to construct the
!  four files used by rvid_box.pro
!
!  13-nov-02/axel: coded
!
      use Cparam
      use General
!
      implicit none
!
      real, dimension (nxgrid,nygrid) :: xy,xy2
      real, dimension (nxgrid,nzgrid) :: xz
      real, dimension (nygrid,nzgrid) :: yz
!
      real, dimension (nx,ny) :: xy_loc,xy2_loc
      real, dimension (nx,nz) :: xz_loc
      real, dimension (ny,nz) :: yz_loc
!
      integer :: ipy,ipz,iproc,it,nt=999999
      integer :: lun,lun1=1,lun2=2,lun3=3,lun4=4
      logical :: eof=.false.
      real :: t
      real :: slice_xpos=0., slice_ypos=0., slice_zpos=0., slice_z2pos=0.
!
      character (len=120) :: file='',fullname='',wfile=''
      character (len=120) :: datadir='data',path=''
      character (len=5) :: chproc=''
      character (len=20) :: field='lnrho'
!
!  read file name
!
      !call getarg (1,field)
      print*,'enter name of variable (lnrho, ux, ..., bz):'
      read*,field
!
      do it=1,nt
        lun=10
!
!  Top Xy-plane:
!  need data where ipz=nprocz-1
!
!tony NOT NECESSARILY SINCE SLICES CAN BE AT ARBITRARY z MESH POINT
      ipz=nprocz-1
      do ipy=0,nprocy-1
        iproc=ipy+nprocy*ipz
        call chn(iproc,chproc,'rvid_box: top xy')
        call safe_character_assign(path,trim(datadir)//'/proc'//chproc)
        call safe_character_assign(file,'/slice_'//trim(field)//'.Xy')
        call safe_character_assign(fullname,trim(path)//trim(file))
        if(it==1) print*,trim(fullname)
        call rslice(trim(fullname),xy2_loc,slice_z2pos,nx,ny,t,it,lun,eof)
        if(eof) goto 999
        xy2(:,1+ipy*ny:ny+ipy*ny)=xy2_loc
      enddo
      call safe_character_assign(wfile,trim(datadir)//trim(file))
      call wslice(trim(wfile),xy2,slice_z2pos,nxgrid,nygrid,t,it,lun1)
!
!  Bottom xy-plane:
!  need data where ipz=0
!
!tony NOT NECESSARILY SINCE SLICES CAN BE AT ARBITRARY z MESH POINT
      ipz=0
      do ipy=0,nprocy-1
        iproc=ipy+nprocy*ipz
        call chn(iproc,chproc,'rvid_box: bottom xy')
        call safe_character_assign(path,trim(datadir)//'/proc'//chproc)
        call safe_character_assign(file,'/slice_'//trim(field)//'.xy')
        call safe_character_assign(fullname,trim(path)//trim(file))
        if(it==1) print*,trim(fullname)
        call rslice(trim(fullname),xy_loc,slice_zpos,nx,ny,t,it,lun,eof)
        if(eof) goto 999
        xy(:,1+ipy*ny:ny+ipy*ny)=xy_loc
      enddo
      call safe_character_assign(wfile,trim(datadir)//trim(file))
      call wslice(trim(wfile),xy,slice_zpos,nxgrid,nygrid,t,it,lun2)
!
!  Front xz-plane:
!  need data where ipy=0
!
!tony NOT NECESSARILY SINCE SLICES CAN BE AT ARBITRARY y MESH POINT
      ipy=0
      do ipz=0,nprocz-1
        iproc=ipy+nprocy*ipz
        call chn(iproc,chproc,'rvid_box: front xz')
        call safe_character_assign(path,trim(datadir)//'/proc'//chproc)
        call safe_character_assign(file,'/slice_'//trim(field)//'.xz')
        call safe_character_assign(fullname,trim(path)//trim(file))
        if(it==1) print*,trim(fullname)
        call rslice(trim(fullname),xz_loc,slice_ypos,nx,nz,t,it,lun,eof)
        if(eof) goto 999
        xz(:,1+ipz*nz:nz+ipz*nz)=xz_loc
      enddo
      call safe_character_assign(wfile,trim(datadir)//trim(file))
      call wslice(trim(wfile),xz,slice_ypos,nxgrid,nzgrid,t,it,lun3)
!
!  Left side yz-plane:
!  need data where ipx=0
!
      do ipz=0,nprocz-1
      do ipy=0,nprocy-1
        iproc=ipy+nprocy*ipz
        call chn(iproc,chproc,'rvid_box: left yz')
        call safe_character_assign(path,trim(datadir)//'/proc'//chproc)
        call safe_character_assign(file,'/slice_'//trim(field)//'.yz')
        call safe_character_assign(fullname,trim(path)//trim(file))
        if(it==1) print*,trim(fullname)
        call rslice(trim(fullname),yz_loc,slice_xpos,ny,nz,t,it,lun,eof)
        if(eof) goto 999
        yz(1+ipy*ny:ny+ipy*ny,1+ipz*nz:nz+ipz*nz)=yz_loc
      enddo
      enddo
      call safe_character_assign(wfile,trim(datadir)//trim(file))
      call wslice(trim(wfile),yz,slice_xpos,nygrid,nzgrid,t,it,lun4)
!
      print*,'written full set of slices at t=',t
      enddo
999   continue
      print*,'finished OK'
      end
!***********************************************************************
    subroutine rslice(file,a,pos,ndim1,ndim2,t,it,lun,eof)
!
!  appending to an existing slice file
!
!  12-nov-02/axel: coded
!
      integer :: ndim1,ndim2
      character (len=*) :: file
      real, dimension (ndim1,ndim2) :: a
      integer :: it,lun
      logical :: eof
      real :: t,pos
!
      if(it==1) open(lun,file=file,form='unformatted')

      pos=0.  ! By default (i.e. if missing from record)
      read(lun,end=999,err=998) a,t,pos
      lun=lun+1
      goto 900
!
!  error: suspect wrong record length
!
998   read(lun,end=999) a,t
      lun=lun+1
      goto 900
!
!  when end of file
!
999   eof=.true.
!
900   continue
    endsubroutine rslice
!***********************************************************************
    subroutine wslice(file,a,pos,ndim1,ndim2,t,it,lun)
!
!  appending to an existing slice file
!
!  12-nov-02/axel: coded
!
      integer :: ndim1,ndim2
      character (len=*) :: file
      real, dimension (ndim1,ndim2) :: a
      integer :: it,lun
      real :: t, pos

!
      if(it==1) open(lun,file=file,form='unformatted')
      write(lun) a,t,pos
!
    endsubroutine wslice
!***********************************************************************
