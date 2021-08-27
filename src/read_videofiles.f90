! $Id$
!
!***********************************************************************
program read_videofiles
!
!  read and combine slices from individual processor directories data
!  /procN, write them to data/, where they can be used by rvid_box.pro
!
!  13-nov-02/axel: coded
!  22-sep-07/axel: changed Xy to xy2, to be compatible with Mac
!  01-aug-12/Bourdin.KIS: rewritten
!  10-apr-16/MR: modifications for Yin-Yang grid
!
  use Cparam
  use General
!
  implicit none
!
  integer :: ipx,ipy,ipz,iproc,it,num_slices,num_frames
  integer :: ipy1=-1,ipy2=-1,ipx1=-1,ipz1=-1,ipz2=-1,ipz3=-1,ipz4=-1
  integer, parameter :: lun=10
  integer :: itdebug=1,n_every=1
  integer :: isep1=0,isep2=0,nyy=1,iyy,ninds
  logical :: lyinyang=.false.
  real :: t
!
  character (len=fnlen) :: file='',fullname='',directory=''
  character (len=fnlen) :: datadir='data',path='',cfield=''
  character (len=labellen) :: field='lnrho', cn_every=''
!
  logical :: exists, lread_slice, lwritten_something=.false.
!
  real :: min_xy,min_xy2,min_xy3,min_xy4,min_xz,min_yz,min_xz2
  real :: max_xy,max_xy2,max_xy3,max_xy4,max_xz,max_yz,max_xz2
!
! Grid data
!
  real, dimension(nygrid) :: y, costh, sinth
  real, dimension(nzgrid) :: z, cosph, sinph
  real :: dy, dz

  integer, dimension(:), allocatable :: inds
  real, dimension(:,:), allocatable :: yz,yzyang
  character, dimension(-1:0) :: trufal=(/'F','T' /)
!
!  Read name of the field (must coincide with file extension)
!
  if (iargc()==0) then
    write(*,'(a)',ADVANCE='NO') 'enter variable (lnrho, uu1, ..., bb3) and stride (e.g. 10): '
    read(*,'(a)') cfield
!
!  Read stride from internal file.
!
    isep1=index(cfield,' ')
    isep2=len(cfield)
    field=cfield(1:isep1)
    if (cfield(isep1+1:isep2)/=' ') read(cfield(isep1:isep2),*) n_every
  else
!
!  Read field and stride command line.
!
    call getarg(1,field)
    call getarg(2,cn_every)
    if (cn_every/='') read(cn_every,*,err=11) n_every
    goto 12
11  print*, 'Invalid value for stride -- set to 1!!!'
    n_every=1
12  continue
  endif

  if (n_every <= 0) then
    print*, 'Invalid value for stride -- set to 1!!!'
    n_every=1
  endif
!
!  Find out the number of available slices and frames to read.
!
  open (lun,file=trim(datadir)//'/tvid.dat',form='formatted',STATUS='old')
  read (lun,*) t, num_slices
  close (lun)
  num_frames = 0
  do it = 1, num_slices
    if (mod(it, n_every) == 0) num_frames = num_frames + 1
  enddo
  if (num_frames == 0) then
    print *, 'No frames are to be written with these settings!'
    stop
  endif
!
!  Loop over all processors to find the positions of the slices.
!  iyy-loop: try to read from 2*ncpus procs as run could also be a Yin-Yang one.
!
  do iyy=0,0   !ncpus,ncpus
    do ipx=0,nprocx-1
      do ipy=0,nprocy-1
        do ipz=0,nprocz-1
          iproc=find_proc(ipx,ipy,ipz)+iyy
          call safe_character_assign(directory, trim(datadir)//'/proc'//itoa(iproc))
          ! check for existence first
          inquire(FILE=trim(directory)//'/slice_position.dat',EXIST=exists)
          if (exists) then
!
!  File found for iyy=ncpus, i.e. iproc=ncpus --> a Yin-Yang grid supposed.
!
            if (iyy>0.and..not.lyinyang) then
              nyy=2
              lyinyang=.true.
              print*, 'This run is treated as a Yin-Yang one!'
            endif
          else
            if (iyy==0 .or. lyinyang) then
              print *, 'slice_position.dat for iproc=', iproc, 'not found!'
              stop
            else
              exit
            endif
          endif
          open(lun,file=trim(directory)//'/slice_position.dat',form='formatted',STATUS='old')
          read(lun,*) lread_slice          ! xy
          if (lread_slice) ipz1=ipz
          read(lun,*) lread_slice          ! xy2
          if (lread_slice) ipz2=ipz
          read(lun,*) lread_slice          ! xy3
          if (lread_slice) ipz3=ipz
          read(lun,*) lread_slice          ! xy4
          if (lread_slice) ipz4=ipz
          read(lun,*) lread_slice          ! xz
          if (lread_slice) ipy1=ipy
          read(lun,*) lread_slice          ! xz2
          if (lread_slice) ipy2=ipy
          read(lun,*) lread_slice          ! yz
          if (lread_slice) ipx1=ipx
          close(lun)
        enddo
      enddo
    enddo
  enddo
!
!  XY-planes:
!
  call read_slice(ipz1,'xy', min_xy, max_xy)
  call read_slice(ipz2,'xy2',min_xy2,max_xy2)
  call read_slice(ipz3,'xy3',min_xy3,max_xy3)
  call read_slice(ipz4,'xy4',min_xy4,max_xy4)
!
!  XZ-planes:
!
  call read_slice(ipy1,'xz', min_xz, max_xz)
  call read_slice(ipy2,'xz2',min_xz2,max_xz2)
!
!  YZ-plane:
!
  if (lyinyang) then
    if (nygrid>1.and.nzgrid>1) then
!
!  If Yin-Yang grid generate (Yin or Yang) grid (assumed uniform) and merge both into yz.
!
      allocate(yzyang(2,nyzgrid),yz(2,2*nyzgrid),inds(nyzgrid))
      dy=pi/2./max(nygrid-1,1)
      dz=3.*pi/2./max(nzgrid-1,1)
      y=(indgen(nygrid)-1)*dy+pi/4.
      z=(indgen(nzgrid)-1)*dz+pi/4
      costh=cos(y); cosph=cos(z); sinth=sin(y); sinph=sin(z)
      call yin2yang_coors(costh,sinth,cosph,sinph,yzyang)
      ninds=merge_yin_yang(y,z,dy,dz,yzyang,yz,inds)
!
!  Hand over merged grid to allow for merging of read-in data.
!
      call read_slice(ipx1,'yz',min_yz,max_yz,yz(:,1:nyzgrid+ninds),inds(1:ninds))
    else
      stop 'Yin-Yang requires non-zero extent in y and z directions'
    endif
  else
    call read_slice(ipx1,'yz',min_yz,max_yz)
  endif
!
!  Print summary.
!
  if (lwritten_something) then

    open (lun,file=trim(datadir)//'/slice_position.dat',form='formatted',STATUS='replace')
    write(lun,*) trufal(min(ipz1,0))
    write(lun,*) trufal(min(ipz2,0))
    write(lun,*) trufal(min(ipz3,0))
    write(lun,*) trufal(min(ipz4,0))
    write(lun,*) trufal(min(ipy1,0))
    write(lun,*) trufal(min(ipy2,0))
    write(lun,*) trufal(min(ipx1,0))
    close(lun)

    print *,'last file read: ',trim(fullname)
    print *,'-------------------------------------------------'
    print *,'minimum and maximum values:'
    if (ipz1/=-1) print *,' xy-plane:',min_xy,max_xy
    if (ipz2/=-1) print *,'xy2-plane:',min_xy2,max_xy2
    if (ipz3/=-1) print *,'xy3-plane:',min_xy3,max_xy3
    if (ipz4/=-1) print *,'xy4-plane:',min_xy4,max_xy4
    if (ipy1/=-1) print *,' xz-plane:',min_xz,max_xz
    if (ipy2/=-1) print *,'xz2-plane:',min_xz2,max_xz2
    if (ipx1/=-1) print *,' yz-plane:',min_yz,max_yz
    print *,'-------------------------------------------------'
    print *,'finished OK'
  endif
!
  select case (trim(field))
    case ('ux','uy','uz','bx','by','bz','Fradx','Frady','Fradz','ax','ay','az','ox','oy','oz')
      print *,""
      print *,"*****************************************************************************"
      print *,"******                WARNING DEPRECATED SLICE NAME                    ******"
      print *,"*****************************************************************************"
      print *,"*** The slice name '"//trim(field)//"'"
      print *,"*** is deprecated and soon will not be supported any longer               ***"
      print *,"*** New slice names are formed by taking the name specified in video.in   ***"
      print *,"*** eg. uu and in the case of vector or other multiple slices appending   ***"
      print *,"*** a number. For example the slice 'ux' is now 'uu1' and the slice 'uz'  ***"
      print *,"*** is now called 'uu3', 'ay'->'aa2' etc. Similarly for aa, bb, oo, uu    ***"
      print *,"*** and Frad slices.                                                      ***"
      print *,"*****************************************************************************"
      print *,""
  endselect
!
  contains
!***********************************************************************
    subroutine read_slice(ipxyz,suffix,glob_min,glob_max,yz,inds)
!
!  Read the slices of a given variable.
!
!  01-Aug-12/Bourdin.KIS: rewritten
!  10-apr-16/MR: modifications for Yin-Yang grid
!
      integer,                           intent(in) :: ipxyz
      character (len=*),                 intent(in) :: suffix
      real,                              intent(out):: glob_min, glob_max
      real,    dimension(:,:), optional, intent(in) :: yz
      integer, dimension(:),   optional, intent(in) :: inds
!
      integer :: frame, slice, ndim1, ndim2, glob_ndim1, glob_ndim2
      integer :: ipx_start, ipx_end, ipy_start, ipy_end, ipz_start, ipz_end
      integer :: i,ind,ninds
      real, dimension (num_frames) :: times
      real, dimension (:,:), allocatable :: loc_slice
      real, dimension (:,:,:,:), allocatable :: glob_slice
      real, dimension (:,:), allocatable :: glob_slice_yy
      real :: slice_pos
      logical :: lexists
!
      if (ipxyz < 0) return
      print *, "read_slice: "//trim(suffix)
!
      ipx_start = 0
      ipx_end = nprocx-1
      ipy_start = 0
      ipy_end = nprocy-1
      ipz_start = 0
      ipz_end = nprocz-1
      if (suffix(1:2) == 'xz') then
        ipy_start = ipxyz
        ipy_end = ipxyz
        ndim1 = nx
        ndim2 = nz
        glob_ndim1 = nxgrid
        glob_ndim2 = nzgrid
      elseif (suffix(1:2) == 'yz') then
        ipx_start = ipxyz
        ipx_end = ipxyz
        ndim1 = ny
        ndim2 = nz
        glob_ndim1 = nygrid
        glob_ndim2 = nzgrid
        if (lyinyang) then
          ninds=size(inds)
          allocate(glob_slice_yy(2*nyzgrid,num_frames))
        endif
      else
        ipz_start = ipxyz
        ipz_end = ipxyz
        ndim1 = nx
        ndim2 = ny
        glob_ndim1 = nxgrid
        glob_ndim2 = nygrid
      endif
      allocate (loc_slice(ndim1,ndim2), glob_slice(glob_ndim1,glob_ndim2,num_frames,0:nyy-1))
!
!  Try to read 
!
      do iyy=0,nyy-1
      do ipz=ipz_start, ipz_end
        do ipy=ipy_start, ipy_end
          do ipx=ipx_start, ipx_end
            iproc=find_proc(ipx,ipy,ipz)+iyy*ncpus
            call safe_character_assign(path,trim(datadir)//'/proc'//itoa(iproc))
            call safe_character_assign(file,'/slice_'//trim(field)//'.'//trim(suffix))
            call safe_character_assign(fullname,trim(path)//trim(file))
!
            if (it <= itdebug) print *, trim(fullname)
            inquire (FILE=trim(fullname), EXIST=lexists)
            if (.not. lexists) then
              print *, "Slice not found: ", fullname
              deallocate (loc_slice, glob_slice)
              return
            endif
!
            open (lun, file=fullname, status='old', form='unformatted')
            ! Loop over frames.
            do frame = 1, num_frames
              ! Loop over available slices.
              do slice = 1, n_every
                it = slice + (frame-1) * n_every
                if (read_slice_file(lun,loc_slice,slice_pos)) then
                  deallocate (loc_slice, glob_slice)
                  return
                endif
              enddo
!
              times(frame) = t
              if (suffix(1:2) == 'xz') then
                glob_slice(1+ipx*nx:nx+ipx*nx,1+ipz*nz:nz+ipz*nz,frame,iyy) = loc_slice
              elseif (suffix(1:2) == 'yz') then
                glob_slice(1+ipy*ny:ny+ipy*ny,1+ipz*nz:nz+ipz*nz,frame,iyy) = loc_slice
              else
                glob_slice(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny,frame,iyy) = loc_slice
              endif
            enddo
            close (lun)
          enddo
        enddo
      enddo
      enddo
!
      call safe_character_assign(fullname,trim(datadir)//trim(file))
      if (lyinyang.and.suffix(1:2) == 'yz') then
        ind=1; iyy=0
        do i=1,2*nygrid
          glob_slice_yy(ind:ind+nzgrid-1,:) = glob_slice(i-iyy*nygrid,:,:,iyy)
          ind=ind+nzgrid
          if (i==nygrid) iyy=1
        enddo
        glob_slice_yy(nyzgrid+1:nyzgrid+ninds,:) = glob_slice_yy(nyzgrid+inds,:)

        glob_min = minval(glob_slice_yy)
        glob_max = maxval(glob_slice_yy)
        call write_slices(fullname,suffix,glob_slice,times,slice_pos,yz,glob_slice_yy(1:nyzgrid+ninds,:))
      else
        glob_min = minval(glob_slice)
        glob_max = maxval(glob_slice)
        call write_slices(fullname,suffix,glob_slice,times,slice_pos)
      endif
      lwritten_something = .true.
      deallocate (loc_slice, glob_slice)
!
    endsubroutine read_slice
!***********************************************************************
    function read_slice_file(lun,a,pos)
!
! Read an existing slice file
!
!  12-nov-02/axel: coded
!  01-Aug-2012/Bourdin.KIS: rewritten
!
      logical :: read_slice_file
      integer, intent(in) :: lun
      real, dimension (:,:), intent(out) :: a
      real, intent(out) :: pos
!
      integer :: ierr
!
      read(lun,iostat=ierr) a, t, pos
      if (ierr > 0) then
        ! IO-error: suspect old file format
        read(lun,iostat=ierr) a, t
        pos = 0. ! Fill missing position value
      endif
      read_slice_file = (ierr > 0) ! Unresolvable IO-error?
!
    endfunction read_slice_file
!***********************************************************************
    subroutine write_slices(filename,suffix,data,times,pos,yz,data_yy)
!
!  Write slices file
!
!  01-Aug-2012/Bourdin.KIS: rewritten
!
      character (len=*),               intent(in) :: filename, suffix
      real, dimension (:,:,:,:),       intent(in) :: data
      real, dimension (:),             intent(in) :: times
      real,                            intent(in) :: pos
      real, dimension (:,:), optional, intent(in) :: yz,data_yy
!
      integer :: frame
!
      open (lun,file=filename,form='unformatted',status='replace')
      if (present(yz)) then
        write (lun) size(yz,2)
        write (lun) yz
      endif
      do frame = 1, num_frames
        if (present(yz)) then
          write (lun) data_yy(:,frame), times(frame), pos
        else
          write (lun) data(:,:,frame,1), times(frame), pos
        endif
        print *, 'written full set of '//trim(suffix)//'-slices at t=', times(frame)
      enddo

      close (lun)
!
    endsubroutine write_slices
!***********************************************************************
endprogram read_videofiles
!***********************************************************************
