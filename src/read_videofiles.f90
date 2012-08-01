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
!  01-Aug-2012/Bourdin.KIS: rewritten
!
  use Cparam
  use General
!
  implicit none
!
  integer :: ipx,ipy,ipz,iproc,it,num_slices,num_frames
  integer :: ipy1=-1,ipx1=-1,ipz1=-1,ipz2=-1,ipz3=-1,ipz4=-1
  integer, parameter :: lun=10
  integer :: itdebug=1,n_every=1
  integer :: isep1=0,isep2=0,idummy
  real :: t
!
  character (len=fnlen) :: file='',fullname='',wfile='',directory=''
  character (len=fnlen) :: datadir='data',path='',cfield=''
  character (len=labellen) :: field='lnrho'
!
  logical :: lsuccess, exists, lread_slice, lwritten_something=.false., lwrite=.true.
!
  real :: min_xy,min_xy2,min_xy3,min_xy4,min_xz,min_yz
  real :: max_xy,max_xy2,max_xy3,max_xy4,max_xz,max_yz
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
  if (cfield(isep1+1:isep2)/=' ') read(cfield(isep1:isep2),*) n_every
!
!  Find out the number of available slices and frames to read.
!
  open (lun,file=trim(datadir)//'/tvid.dat',form='formatted',STATUS='old')
  read (lun,*) t, num_slices
  close (lun)
  num_slices = num_slices - 1
  num_frames = 0
  do it = 1, num_slices, n_every
    if (mod(it, n_every) == 0) num_frames = num_frames + 1
  enddo
!
!  loop over all processors to find the positions of the slices.
!  Therefore read all slice_positions.dat
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
        if (.not. exists) then
          print*,'slice_position.dat for iproc=',iproc,'not found!'
          stop
        endif
        open(lun,file=trim(directory)//'/slice_position.dat',form='formatted',STATUS='old')
        read(lun,'(l5,i5)') lread_slice,idummy
        if (lread_slice) ipz1=ipz
        read(lun,'(l5,i5)') lread_slice,idummy
        if (lread_slice) ipz2=ipz
        read(lun,'(l5,i5)') lread_slice,idummy
        if (lread_slice) ipz3=ipz
        read(lun,'(l5,i5)') lread_slice,idummy
        if (lread_slice) ipz4=ipz
        read(lun,'(l5,i5)') lread_slice,idummy
        if (lread_slice) ipy1=ipy
        read(lun,'(l5,i5)') lread_slice,idummy
        if (lread_slice) ipx1=ipx
        close(lun)
      enddo
    enddo
  enddo
!
    ! XY-planes:
  call read_slice(ipz1,'xy',min_xy,max_xy)
  call read_slice(ipz2,'xy2',min_xy2,max_xy2)
  call read_slice(ipz3,'xy3',min_xy3,max_xy3)
  call read_slice(ipz4,'xy4',min_xy4,max_xy4)
!
  ! XZ-plane:
  call read_slice(ipy1,'xz',min_xz,max_xz)
!
  ! YZ-plane:
  call read_slice(ipx1,'yz',min_yz,max_yz)
!
  if (lwritten_something) then
    print *,'last file read: ',trim(fullname)
    print *,'-------------------------------------------------'
    print *,'minimum and maximum values:'
    if (ipz1/=-1) print *,' xy-plane:',min_xy,max_xy
    if (ipz2/=-1) print *,'xy2-plane:',min_xy2,max_xy2
    if (ipz3/=-1) print *,'xy3-plane:',min_xy3,max_xy3
    if (ipz4/=-1) print *,'xy4-plane:',min_xy4,max_xy4
    if (ipy1/=-1) print *,' xz-plane:',min_xz,max_xz
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
    subroutine read_slice(ipxyz,suffix,glob_min,glob_max)
!
!  Read the slices of a given variable.
!
!  01-Aug-2012/Bourdin.KIS: rewritten
!
      integer, intent(in) :: ipxyz
      character (len=*), intent(in) :: suffix
      real, intent(out) :: glob_min, glob_max
!
      integer :: frame, slice, ndim1, ndim2, glob_ndim1, glob_ndim2, ipx, ipy, ipz
      integer :: ipx_start, ipx_end, ipy_start, ipy_end, ipz_start, ipz_end
      real, dimension (:), allocatable :: times
      real, dimension (:,:), allocatable :: loc_slice
      real, dimension (:,:,:), allocatable :: glob_slice
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
      else
        ipz_start = ipxyz
        ipz_end = ipxyz
        ndim1 = nx
        ndim2 = ny
        glob_ndim1 = nxgrid
        glob_ndim2 = nygrid
      endif
      allocate (times(num_frames), loc_slice(ndim1,ndim2), glob_slice(glob_ndim1,glob_ndim2,num_frames))
!
      do ipz=ipz_start, ipz_end
        do ipy=ipy_start, ipy_end
          do ipx=ipx_start, ipx_end
            iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
            call safe_character_assign(path,trim(datadir)//'/proc'//itoa(iproc))
            call safe_character_assign(file,'/slice_'//trim(field)//'.'//trim(suffix))
            call safe_character_assign(fullname,trim(path)//trim(file))
!
            if (it <= itdebug) print *, trim(fullname)
            inquire (FILE=trim(fullname), EXIST=lexists)
            if (.not. lexists) then
              print *, "Slice not found: ", fullname
              deallocate (times, loc_slice, glob_slice)
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
                  deallocate (times, loc_slice, glob_slice)
                  return
                endif
              enddo
!
              times(frame) = t
              if (suffix(1:2) == 'xz') then
                glob_slice(1+ipx*nx:nx+ipx*nx,1+ipz*nz:nz+ipz*nz,frame) = loc_slice
              elseif (suffix(1:2) == 'yz') then
                glob_slice(1+ipy*ny:ny+ipy*ny,1+ipz*nz:nz+ipz*nz,frame) = loc_slice
              else
                glob_slice(1+ipx*nx:nx+ipx*nx,1+ipy*ny:ny+ipy*ny,frame) = loc_slice
              endif
            enddo
            close (lun)
          enddo
        enddo
      enddo
!
      glob_min = minval(glob_slice)
      glob_max = maxval(glob_slice)
!
      call safe_character_assign(fullname,trim(datadir)//trim(file))
      call write_slices(fullname,glob_slice,times,slice_pos)
      lwritten_something = .true.
      do frame = 1, num_frames
        print *, 'written full set of '//trim(suffix)//'-slices at t=', times(frame)
      enddo
      deallocate (times, loc_slice, glob_slice)
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
        ! IO-error: suspect wrong record length
        read(lun,iostat=ierr) a, t
        pos = 0. ! Fill missing value
      endif
      read_slice_file = (ierr > 0) ! Unresolvable IO-error?
!
    endfunction read_slice_file
!***********************************************************************
    subroutine write_slices(filename,data,times,pos)
!
!  Write slices file
!
!  01-Aug-2012/Bourdin.KIS: rewritten
!
      character (len=*), intent(in) :: filename
      real, dimension (:,:,:), intent(in) :: data
      real, dimension (:), intent(in) :: times
      real, intent(in) :: pos
!
      integer :: frame
!
      open (lun,file=filename,form='unformatted',status='replace')
      do frame = 1, num_frames
        write (lun) data(:,:,frame), times(frame), pos
      enddo
      close (lun)
!
    endsubroutine write_slices
!***********************************************************************
endprogram read_videofiles
!***********************************************************************

