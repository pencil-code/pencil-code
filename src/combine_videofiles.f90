! $Id$
!
!***********************************************************************
    program combine_videofiles
!
!  combine slices for two different frequencies and combine
!  using a suitable color table layout.
!
!  24-feb-09/axel: adapted from read_videofiles.x
!
      use Cparam
      use General
!
      implicit none
!
      real, dimension (nxgrid,nygrid) :: xy,xy1,xy2
      real, dimension (nxgrid,nzgrid) :: xz,xz1,xz2
      real, dimension (nygrid,nzgrid) :: yz,yz1,yz2
!
      integer :: it
      integer :: lun0,lun1,lun2
      logical :: eof
      real :: t,fac1,fac2,maxval1,maxval2
!
      character (len=120) :: dir,wfile,rfile1,rfile2
      character (len=20) :: field
!
!  read name of the field (must coincide with file extension)
!
      field = 'Jrad'
      write(*,'(a)',ADVANCE='NO') 'enter name of variable (e.g. Jrad): '
      read*,field
!
      fac1 = 1.0
      fac2 = 1.0
!      write(*,'(a)',ADVANCE='NO') 'enter two factors for multiplication (e.g. 1.0, 1.0): '
!      read*,fac1,fac2 ! [PABourdin] What is this? Please document and fix above question.
!
      dir='data/slice_'
      dir='data/proc0/slice_'
!
!  loop through all times and convert xy, xz, and yz files
!  reset the lun to 10 every time. This guarantees unique luns every time round
!
      it=0
      do while (.true.)
        it=it + 1
        lun0=0
        lun1=10
        lun2=20
        call safe_character_assign(rfile1,trim(dir)//trim(field)//'1'//'.xy2')
        call safe_character_assign(rfile2,trim(dir)//trim(field)//'2'//'.xy2')
        call read_slice(trim(rfile1),xy1,nxgrid,nygrid,t,it,lun1,eof); if (eof) exit
        call read_slice(trim(rfile2),xy2,nxgrid,nygrid,t,it,lun2,eof); if (eof) exit
        maxval1=max(maxval1,maxval(fac1*xy1))
        maxval2=max(maxval2,maxval(fac2*xy2))
        xy=int(16*fac1*xy1)+16*int(16*fac2*xy2)
        xy=16*int(16*fac1*xy1)+int(16*fac2*xy2)
        call safe_character_assign(wfile,trim(dir)//trim(field)//'.xy2')
        call append_slice(trim(wfile),xy,nxgrid,nygrid,t,it,lun0)
!
        call safe_character_assign(rfile1,trim(dir)//trim(field)//'1'//'.xy')
        call safe_character_assign(rfile2,trim(dir)//trim(field)//'2'//'.xy')
        call read_slice(trim(rfile1),xy1,nxgrid,nygrid,t,it,lun1,eof); if (eof) exit
        call read_slice(trim(rfile2),xy2,nxgrid,nygrid,t,it,lun2,eof); if (eof) exit
        maxval1=max(maxval1,maxval(fac1*xy1))
        maxval2=max(maxval2,maxval(fac2*xy2))
        xy=int(16*fac1*xy1)+16*int(16*fac2*xy2)
        xy=16*int(16*fac1*xy1)+int(16*fac2*xy2)
        call safe_character_assign(wfile,trim(dir)//trim(field)//'.xy')
        call append_slice(trim(wfile),xy,nxgrid,nygrid,t,it,lun0)
!
        call safe_character_assign(rfile1,trim(dir)//trim(field)//'1'//'.xz')
        call safe_character_assign(rfile2,trim(dir)//trim(field)//'2'//'.xz')
        call read_slice(trim(rfile1),xz1,nxgrid,nzgrid,t,it,lun1,eof); if (eof) exit
        call read_slice(trim(rfile2),xz2,nxgrid,nzgrid,t,it,lun2,eof); if (eof) exit
        maxval1=max(maxval1,maxval(fac1*xz1))
        maxval2=max(maxval2,maxval(fac2*xz2))
        xz=int(16*fac1*xz1)+16*int(16*fac2*xz2)
        xz=16*int(16*fac1*xz1)+int(16*fac2*xz2)
        call safe_character_assign(wfile,trim(dir)//trim(field)//'.xz')
        call append_slice(trim(wfile),xz,nxgrid,nzgrid,t,it,lun0)
!
        call safe_character_assign(rfile1,trim(dir)//trim(field)//'1'//'.yz')
        call safe_character_assign(rfile2,trim(dir)//trim(field)//'2'//'.yz')
        call read_slice(trim(rfile1),yz1,nygrid,nzgrid,t,it,lun1,eof); if (eof) exit
        call read_slice(trim(rfile2),yz2,nygrid,nzgrid,t,it,lun2,eof); if (eof) exit
        maxval1=max(maxval1,maxval(fac1*yz1))
        maxval2=max(maxval2,maxval(fac2*yz2))
        yz=int(16*fac1*yz1)+16*int(16*fac2*yz2)
        yz=16*int(16*fac1*yz1)+int(16*fac2*yz2)
        call safe_character_assign(wfile,trim(dir)//trim(field)//'.yz')
        call append_slice(trim(wfile),yz,nygrid,nzgrid,t,it,lun0)
!
      enddo
!
      print*,'maxval1=',maxval1
      print*,'maxval2=',maxval2
    end
!***********************************************************************
    subroutine read_slice(file,a,ndim1,ndim2,t,it,lun,eof)
!
! read an existing slice file
!
!  12-nov-02/axel: coded
!
      integer :: ndim1,ndim2
      character (len=*) :: file
      real, dimension (ndim1,ndim2) :: a
      integer :: it,lun
      logical, intent(out) :: eof
      real :: t
      integer :: io_code
!
      eof=.false.
      if (it==1) open(lun,file=file,status='old',form='unformatted')
      read(lun,iostat=io_code) a,t
      if (io_code < 0) then
        ! end of file
        eof=.true.
      else
        lun=lun+1
      endif
!
    endsubroutine read_slice
!***********************************************************************
    subroutine append_slice(file,a,ndim1,ndim2,t,it,lun)
!
!  append to a slice file
!
!  12-nov-02/axel: coded
!
      integer :: ndim1,ndim2
      character (len=*) :: file
      real, dimension (ndim1,ndim2) :: a
      integer :: it,lun
      real :: t
      real, parameter :: pos=0.
!
      if (it==1) open(lun,file=file,form='unformatted')
      write(lun) a,t,pos
      lun=lun+1
!
    endsubroutine append_slice
!***********************************************************************
