      program slices
      implicit none
! Converts 2D slices from var.dat snapshot to the Dx readable 
! binary data files.
! The data are read in again by .dx routines on regular rectangular
! mesh. The dx_slices.f90 should be copied to the 
! directory of the run. Compile dx_slices.f90 and run it. 
! .net files should be linked to run directory (from pencil-code/dx).
! Then, Data Explorer is ready to be started from run directory.
! This slices routine is better to be integrated into the main code - to
! do in the future.

! Include grid sizes
      include 'src/cparam.local'
      include 'src/cparam.inc'

! This was taken from cparam.f90 file
  integer, parameter :: nghost=3
  integer, parameter :: nx=nxgrid,ny=nygrid/nprocy,nz=nzgrid/nprocz,nw=nx*ny*nz
  integer, parameter :: mx=nx+2*nghost,l1=1+nghost,l2=mx-nghost
  integer, parameter :: my=ny+2*nghost,m1=1+nghost,m2=my-nghost
  integer, parameter :: mz=nz+2*nghost,n1=1+nghost,n2=mz-nghost
  integer, parameter :: mw=mx*my*mz
  integer, parameter :: l1i=l1+nghost-1,l2i=l2-nghost+1
  integer, parameter :: m1i=m1+nghost-1,m2i=m2-nghost+1
  integer, parameter :: n1i=n1+nghost-1,n2i=n2-nghost+1

! These are definitions intrincis for todx_magn.f90
      real, dimension(mx,my,mz,mvar) :: f
      real, dimension(mx) :: x
      real, dimension(my) :: y
      real, dimension(mz) :: z    
      real, dimension(mx,my,mz,3) :: bbb
      real :: x00,y00,z00,t,dx,dy,dz
      integer :: i,j,k,n
      integer :: iux,iuy,iuz,ilnrho,iss,iax,iay,iaz
      character (len=2) :: sltype
      integer :: isltype,isl,nfile

! Read snapshot file
      open(unit=17, file='./tmp/proc0/var.dat', form='unformatted')
      read(17) f
      read(17) t,x,y,z,dx,dy,dz
      close(17)
      print*,'t=',t,' dx=',dx,' dy=',dy,' dz=',dz
      print*,'mx=',mx,' my=',my,' mz=',mz

! Define here the position of the corner of the box formed by only 
! true grid points
      x00=x(nghost+1)
      y00=y(nghost+1)
      z00=z(nghost+1)
      print*, 'x0=',x00,' y0=',y00,' z0=',z00
      print*,'mvar=',mvar

! Read in the type and position of the slice
      print*,'type of slice, xy, yz, or xz ?'
      read*, sltype
      print*,'number of the slice files (must be from 1 to 9) ?'
      read*, nfile
      select case(sltype)
         case('xy') 
            isltype=1
            print*,'position of the slice on z-axis of true grid, isl ?'
            print*,'must be 1 .le. isl .le.',nzgrid
            read*,isl
            if(isl.lt.1.or.isl.gt.nzgrid) then
               print*,'isl is out of grid limits'
               stop
            endif   

         case('yz')
            isltype=2
            print*,'position of the slice on x-axis of true grid, isl ?'
            print*,'must be 1 .le. isl .le.',nxgrid
            read*,isl
            if(isl.lt.1.or.isl.gt.nxgrid) then
               print*,'isl is out of grid limits'
               stop
            endif           

         case('xz')
            isltype=3
            print*,'position of the slice on y-axis of true grid, isl ?'
            print*,'must be 1 .le. isl .le.',nygrid
            read*,isl
            if(isl.lt.1.or.isl.gt.nygrid) then
               print*,'isl is out of grid limits'
               stop
            endif

         case default
        !
        !  Catch unknown values
        !
        print*,'No such value for sltype: ', trim(sltype)
        stop
      endselect     

!!
!!! Write out slices -- global select 
!!
      select case(sltype)

!!!!!!!
!!! Case of xy-slice
!!!!!!!
      case('xy')

! Write out parallel velocity field in the slice
      iux=1
      iuy=2
      iuz=3
      open(unit=18, file='./tmp/uuu_dx_sl'//char(nfile+48)//'.bin', &
        form='unformatted')
      write(18) (((real(f(i,j,isl,n)),n=iux,iuy),j=nghost+1,my-nghost), &
        i=nghost+1,mx-nghost)
      close(18)

! Write out ln density field for Dx
      if(mvar.ge.4) then
      ilnrho=4
      open(unit=19, file='./tmp/lnrho_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(19) ((real(f(i,j,isl,ilnrho)),j=nghost+1,my-nghost), & 
        i=nghost+1,mx-nghost)
      close(19)
      endif

! Write out entropy field for Dx
      if(mvar.ge.5) then
      iss=5
      open(unit=20, file='./tmp/ss_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(20) ((real(f(i,j,isl,iss)),j=nghost+1,my-nghost), & 
        i=nghost+1,mx-nghost)
      close(20)
      endif

! Create .dx scripts describing the data in Data Explorer format
      open(unit=21, file='./tmp/uuu_sl'//char(nfile+48)//'.dx')
      write(21,*) 'object 1 class gridpositions counts ', mx-2*nghost, &
        my-2*nghost
      write(21,*) 'origin  ', x00, y00 
      write(21,*) 'delta  ', dx,'  0   '
      write(21,*) 'delta  ', '  0  ',dy
      write(21,*) '  '
      write(21,*) 'object 2 class gridconnections counts ', mx-2*nghost, &
        my-2*nghost
      write(21,*) '  '
      write(21,*) 'object 3 class array type float rank 1 shape 2'
      write(21,*) 'items ', (mx-2*nghost)*(my-2*nghost), ' lsb binary'
      write(21,*) 'data file "uuu_dx_sl'//char(nfile+48)//'.bin", 4'
      write(21,*) '  '
      write(21,*) 'attribute "dep" string "positions" '
      write(21,*) '  '
      write(21,*) 'object "regular positions regular connections" class field'
      write(21,*) 'component "positions" value 1'
      write(21,*) 'component "connections" value 2'
      write(21,*) 'component "data" value 3'
      write(21,*) 'end'
      close(21)

      if(mvar.ge.4) then         
      open(unit=22, file='./tmp/lnrho_sl'//char(nfile+48)//'.dx')
      write(22,*) 'object 1 class gridpositions counts ', mx-2*nghost,  &
        my-2*nghost
      write(22,*) 'origin  ', x00, y00 
      write(22,*) 'delta  ', dx,'  0   '
      write(22,*) 'delta  ', '  0  ',dy
      write(22,*) '  '
      write(22,*) 'object 2 class gridconnections counts ', mx-2*nghost,  &
        my-2*nghost
      write(22,*) '  '
      write(22,*) 'object 3 class array type float rank 0'
      write(22,*) 'items ', (mx-2*nghost)*(my-2*nghost), ' lsb binary'
      write(22,*) 'data file "lnrho_dx_sl'//char(nfile+48)//'.bin", 4'
      write(22,*) '  '
      write(22,*) 'attribute "dep" string "positions" '
      write(22,*) '  '
      write(22,*) 'object "regular positions regular connections" class field'
      write(22,*) 'component "positions" value 1'
      write(22,*) 'component "connections" value 2'
      write(22,*) 'component "data" value 3'
      write(22,*) 'end'
      close(22)
      endif

      if(mvar.ge.5) then
      open(unit=23, file='./tmp/ss_sl'//char(nfile+48)//'.dx')
      write(23,*) 'object 1 class gridpositions counts ', mx-2*nghost, &
        my-2*nghost
      write(23,*) 'origin  ', x00, y00 
      write(23,*) 'delta  ', dx,'  0   '
      write(23,*) 'delta  ', '  0  ',dy
      write(23,*) '  '
      write(23,*) 'object 2 class gridconnections counts ', mx-2*nghost, &
        my-2*nghost
      write(23,*) '  '
      write(23,*) 'object 3 class array type float rank 0'
      write(23,*) 'items ', (mx-2*nghost)*(my-2*nghost), ' lsb binary'
      write(23,*) 'data file "ss_dx_sl'//char(nfile+48)//'.bin", 4'
      write(23,*) '  '
      write(23,*) 'attribute "dep" string "positions" '
      write(23,*) '  '
      write(23,*) 'object "regular positions regular connections" class field'
      write(23,*) 'component "positions" value 1'
      write(23,*) 'component "connections" value 2'
      write(23,*) 'component "data" value 3'
      write(23,*) 'end'
      close(23)
      endif



      if(mvar.ge.6) then
! Process magnetic part if mvar>5: no kinematic dynamo yet, needs to be added
! Write out 2D slice of parallel vector potential for Dx
      iax=6
      iay=7
      iaz=8
      open(unit=18, file='./tmp/aaa_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(18) (((real(f(i,j,isl,n)),n=iax,iay),j=nghost+1,my-nghost), & 
        i=nghost+1,mx-nghost)
      close(18)
      
! Calculate magnetic field, using central differences.
      do i=2,mx-1
         do j=2,my-1
            do k=2,mz-1
               bbb(i,j,k,1)=(f(i,j+1,k,iaz)-f(i,j-1,k,iaz))/2./dy-  &
       (f(i,j,k+1,iay)-f(i,j,k-1,iay))/2./dz
               bbb(i,j,k,2)=(f(i,j,k+1,iax)-f(i,j,k-1,iax))/2./dz-  &
       (f(i+1,j,k,iaz)-f(i-1,j,k,iaz))/2./dx
               bbb(i,j,k,3)=(f(i+1,j,k,iay)-f(i-1,j,k,iay))/2./dx-  &
       (f(i,j+1,k,iax)-f(i,j-1,k,iax))/2./dy               
            enddo
         enddo
      enddo   

! Write out parallel magnetic field for Dx
      open(unit=18, file='./tmp/bbb_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(18) (((real(bbb(i,j,isl,n)),n=1,2),j=nghost+1,my-nghost), & 
        i=nghost+1,mx-nghost)
      close(18)

! Create  .dx scripts describing magnetic data 
      open(unit=24, file='./tmp/aaa_sl'//char(nfile+48)//'.dx')
      write(24,*) 'object 1 class gridpositions counts ', mx-2*nghost,  &
        my-2*nghost
      write(24,*) 'origin  ', x00, y00 
      write(24,*) 'delta  ', dx,'  0   '
      write(24,*) 'delta  ', '  0  ',dy
      write(24,*) '  '
      write(24,*) 'object 2 class gridconnections counts ', mx-2*nghost,  &
        my-2*nghost
      write(24,*) '  '
      write(24,*) 'object 3 class array type float rank 1 shape 2'
      write(24,*) 'items ', (mx-2*nghost)*(my-2*nghost), ' lsb binary'
      write(24,*) 'data file "aaa_dx_sl'//char(nfile+48)//'.bin", 4'
      write(24,*) '  '
      write(24,*) 'attribute "dep" string "positions" '
      write(24,*) '  '
      write(24,*) 'object "regular positions regular connections" class field'
      write(24,*) 'component "positions" value 1'
      write(24,*) 'component "connections" value 2'
      write(24,*) 'component "data" value 3'
      write(24,*) 'end'
      close(24)

      open(unit=25, file='./tmp/bbb_sl'//char(nfile+48)//'.dx')
      write(25,*) 'object 1 class gridpositions counts ', mx-2*nghost, &
        my-2*nghost
      write(25,*) 'origin  ', x00, y00 
      write(25,*) 'delta  ', dx,'  0   '
      write(25,*) 'delta  ', '  0  ',dy
      write(25,*) '  '
      write(25,*) 'object 2 class gridconnections counts ', mx-2*nghost, &
        my-2*nghost
      write(25,*) '  '
      write(25,*) 'object 3 class array type float rank 1 shape 2'
      write(25,*) 'items ', (mx-2*nghost)*(my-2*nghost), ' lsb binary'
      write(25,*) 'data file "bbb_dx_sl'//char(nfile+48)//'.bin", 4'
      write(25,*) '  '
      write(25,*) 'attribute "dep" string "positions" '
      write(25,*) '  '
      write(25,*) 'object "regular positions regular connections" class field'
      write(25,*) 'component "positions" value 1'
      write(25,*) 'component "connections" value 2'
      write(25,*) 'component "data" value 3'
      write(25,*) 'end'
      close(25)


      endif

!!!!!!!!
!!!  Case of yz-slice
!!!!!!!!
      case('yz')

! Write out parallel velocity field in the slice
      iux=1
      iuy=2
      iuz=3
      open(unit=18, file='./tmp/uuu_dx_sl'//char(nfile+48)//'.bin', &
        form='unformatted')
      write(18) (((real(f(isl,j,k,n)),n=iuy,iuz),k=nghost+1,mz-nghost), &
        j=nghost+1,my-nghost)
      close(18)

! Write out ln density field for Dx
      if(mvar.ge.4) then
      ilnrho=4
      open(unit=19, file='./tmp/lnrho_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(19) ((real(f(isl,j,k,ilnrho)),k=nghost+1,mz-nghost), & 
        j=nghost+1,my-nghost)
      close(19)
      endif

! Write out entropy field for Dx
      if(mvar.ge.5) then
      iss=5
      open(unit=20, file='./tmp/ss_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(20) ((real(f(isl,j,k,iss)),k=nghost+1,mz-nghost), & 
        j=nghost+1,my-nghost)
      close(20)
      endif

! Create .dx scripts describing the data in Data Explorer format
      open(unit=21, file='./tmp/uuu_sl'//char(nfile+48)//'.dx')
      write(21,*) 'object 1 class gridpositions counts ', my-2*nghost, &
        mz-2*nghost
      write(21,*) 'origin  ', y00, z00 
      write(21,*) 'delta  ', dy,'  0   '
      write(21,*) 'delta  ', '  0  ',dz
      write(21,*) '  '
      write(21,*) 'object 2 class gridconnections counts ', my-2*nghost, &
        mz-2*nghost
      write(21,*) '  '
      write(21,*) 'object 3 class array type float rank 1 shape 2'
      write(21,*) 'items ', (my-2*nghost)*(mz-2*nghost), ' lsb binary'
      write(21,*) 'data file "uuu_dx_sl'//char(nfile+48)//'.bin", 4'
      write(21,*) '  '
      write(21,*) 'attribute "dep" string "positions" '
      write(21,*) '  '
      write(21,*) 'object "regular positions regular connections" class field'
      write(21,*) 'component "positions" value 1'
      write(21,*) 'component "connections" value 2'
      write(21,*) 'component "data" value 3'
      write(21,*) 'end'
      close(21)

      if(mvar.ge.4) then         
      open(unit=22, file='./tmp/lnrho_sl'//char(nfile+48)//'.dx')
      write(22,*) 'object 1 class gridpositions counts ', my-2*nghost,  &
        mz-2*nghost
      write(22,*) 'origin  ', y00, z00 
      write(22,*) 'delta  ', dy,'  0   '
      write(22,*) 'delta  ', '  0  ',dz
      write(22,*) '  '
      write(22,*) 'object 2 class gridconnections counts ', my-2*nghost,  &
        mz-2*nghost
      write(22,*) '  '
      write(22,*) 'object 3 class array type float rank 0'
      write(22,*) 'items ', (my-2*nghost)*(mz-2*nghost), ' lsb binary'
      write(22,*) 'data file "lnrho_dx_sl'//char(nfile+48)//'.bin", 4'
      write(22,*) '  '
      write(22,*) 'attribute "dep" string "positions" '
      write(22,*) '  '
      write(22,*) 'object "regular positions regular connections" class field'
      write(22,*) 'component "positions" value 1'
      write(22,*) 'component "connections" value 2'
      write(22,*) 'component "data" value 3'
      write(22,*) 'end'
      close(22)
      endif

      if(mvar.ge.5) then
      open(unit=23, file='./tmp/ss_sl'//char(nfile+48)//'.dx')
      write(23,*) 'object 1 class gridpositions counts ', my-2*nghost, &
        mz-2*nghost
      write(23,*) 'origin  ', y00, z00 
      write(23,*) 'delta  ', dy,'  0   '
      write(23,*) 'delta  ', '  0  ',dz
      write(23,*) '  '
      write(23,*) 'object 2 class gridconnections counts ', my-2*nghost, &
        mz-2*nghost
      write(23,*) '  '
      write(23,*) 'object 3 class array type float rank 0'
      write(23,*) 'items ', (my-2*nghost)*(mz-2*nghost), ' lsb binary'
      write(23,*) 'data file "ss_dx_sl'//char(nfile+48)//'.bin", 4'
      write(23,*) '  '
      write(23,*) 'attribute "dep" string "positions" '
      write(23,*) '  '
      write(23,*) 'object "regular positions regular connections" class field'
      write(23,*) 'component "positions" value 1'
      write(23,*) 'component "connections" value 2'
      write(23,*) 'component "data" value 3'
      write(23,*) 'end'
      close(23)
      endif



      if(mvar.ge.6) then
! Process magnetic part if mvar>5: no kinematic dynamo yet, needs to be added
! Write out 2D slice of parallel vector potential for Dx
      iax=6
      iay=7
      iaz=8
      open(unit=18, file='./tmp/aaa_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(18) (((real(f(isl,j,k,n)),n=iay,iaz),k=nghost+1,mz-nghost), & 
        j=nghost+1,my-nghost)
      close(18)

! Calculate magnetic field, using central differences.
      do i=2,mx-1
         do j=2,my-1
            do k=2,mz-1
               bbb(i,j,k,1)=(f(i,j+1,k,iaz)-f(i,j-1,k,iaz))/2./dy-  &
       (f(i,j,k+1,iay)-f(i,j,k-1,iay))/2./dz
               bbb(i,j,k,2)=(f(i,j,k+1,iax)-f(i,j,k-1,iax))/2./dz-  &
       (f(i+1,j,k,iaz)-f(i-1,j,k,iaz))/2./dx
               bbb(i,j,k,3)=(f(i+1,j,k,iay)-f(i-1,j,k,iay))/2./dx-  &
       (f(i,j+1,k,iax)-f(i,j-1,k,iax))/2./dy               
            enddo
         enddo
      enddo   

! Write out parallel magnetic field for Dx
      open(unit=18, file='./tmp/bbb_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(18) (((real(bbb(isl,j,k,n)),n=2,3),k=nghost+1,mz-nghost), & 
        j=nghost+1,my-nghost)
      close(18)

! Create  .dx scripts describing magnetic data 
      open(unit=24, file='./tmp/aaa_sl'//char(nfile+48)//'.dx')
      write(24,*) 'object 1 class gridpositions counts ', my-2*nghost,  &
        mz-2*nghost
      write(24,*) 'origin  ', y00, z00 
      write(24,*) 'delta  ', dy,'  0   '
      write(24,*) 'delta  ', '  0  ',dz
      write(24,*) '  '
      write(24,*) 'object 2 class gridconnections counts ', my-2*nghost,  &
        mz-2*nghost
      write(24,*) '  '
      write(24,*) 'object 3 class array type float rank 1 shape 2'
      write(24,*) 'items ', (my-2*nghost)*(mz-2*nghost), ' lsb binary'
      write(24,*) 'data file "aaa_dx_sl'//char(nfile+48)//'.bin", 4'
      write(24,*) '  '
      write(24,*) 'attribute "dep" string "positions" '
      write(24,*) '  '
      write(24,*) 'object "regular positions regular connections" class field'
      write(24,*) 'component "positions" value 1'
      write(24,*) 'component "connections" value 2'
      write(24,*) 'component "data" value 3'
      write(24,*) 'end'
      close(24)

      open(unit=25, file='./tmp/bbb_sl'//char(nfile+48)//'.dx')
      write(25,*) 'object 1 class gridpositions counts ', my-2*nghost, &
        mz-2*nghost
      write(25,*) 'origin  ', y00, z00 
      write(25,*) 'delta  ', dy,'  0   '
      write(25,*) 'delta  ', '  0  ',dz
      write(25,*) '  '
      write(25,*) 'object 2 class gridconnections counts ', my-2*nghost, &
        mz-2*nghost
      write(25,*) '  '
      write(25,*) 'object 3 class array type float rank 1 shape 2'
      write(25,*) 'items ', (my-2*nghost)*(mz-2*nghost), ' lsb binary'
      write(25,*) 'data file "bbb_dx_sl'//char(nfile+48)//'.bin", 4'
      write(25,*) '  '
      write(25,*) 'attribute "dep" string "positions" '
      write(25,*) '  '
      write(25,*) 'object "regular positions regular connections" class field'
      write(25,*) 'component "positions" value 1'
      write(25,*) 'component "connections" value 2'
      write(25,*) 'component "data" value 3'
      write(25,*) 'end'
      close(25)


      endif


!!!!!!!!
!!!  Case of xz-slice
!!!!!!!!
      case('xz')

! Write out parallel velocity field in the slice
      iux=1
      iuy=2
      iuz=3
      open(unit=18, file='./tmp/uuu_dx_sl'//char(nfile+48)//'.bin', &
        form='unformatted')
      write(18) (((real(f(i,isl,k,n)),n=iux,iuz,2),k=nghost+1,mz-nghost), &
        i=nghost+1,mx-nghost)
      close(18)

! Write out ln density field for Dx
      if(mvar.ge.4) then
      ilnrho=4
      open(unit=19, file='./tmp/lnrho_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(19) ((real(f(i,isl,k,ilnrho)),k=nghost+1,mz-nghost), & 
        i=nghost+1,mx-nghost)
      close(19)
      endif

! Write out entropy field for Dx
      if(mvar.ge.5) then
      iss=5
      open(unit=20, file='./tmp/ss_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(20) ((real(f(i,isl,k,iss)),k=nghost+1,mz-nghost), & 
        i=nghost+1,mx-nghost)
      close(20)
      endif

! Create .dx scripts describing the data in Data Explorer format
      open(unit=21, file='./tmp/uuu_sl'//char(nfile+48)//'.dx')
      write(21,*) 'object 1 class gridpositions counts ', mx-2*nghost, &
        mz-2*nghost
      write(21,*) 'origin  ', x00, z00 
      write(21,*) 'delta  ', dx,'  0   '
      write(21,*) 'delta  ', '  0  ',dz
      write(21,*) '  '
      write(21,*) 'object 2 class gridconnections counts ', mx-2*nghost, &
        mz-2*nghost
      write(21,*) '  '
      write(21,*) 'object 3 class array type float rank 1 shape 2'
      write(21,*) 'items ', (mx-2*nghost)*(mz-2*nghost), ' lsb binary'
      write(21,*) 'data file "uuu_dx_sl'//char(nfile+48)//'.bin", 4'
      write(21,*) '  '
      write(21,*) 'attribute "dep" string "positions" '
      write(21,*) '  '
      write(21,*) 'object "regular positions regular connections" class field'
      write(21,*) 'component "positions" value 1'
      write(21,*) 'component "connections" value 2'
      write(21,*) 'component "data" value 3'
      write(21,*) 'end'
      close(21)

      if(mvar.ge.4) then         
      open(unit=22, file='./tmp/lnrho_sl'//char(nfile+48)//'.dx')
      write(22,*) 'object 1 class gridpositions counts ', mx-2*nghost,  &
        mz-2*nghost
      write(22,*) 'origin  ', x00, z00 
      write(22,*) 'delta  ', dx,'  0   '
      write(22,*) 'delta  ', '  0  ',dz
      write(22,*) '  '
      write(22,*) 'object 2 class gridconnections counts ', mx-2*nghost,  &
        mz-2*nghost
      write(22,*) '  '
      write(22,*) 'object 3 class array type float rank 0'
      write(22,*) 'items ', (mx-2*nghost)*(mz-2*nghost), ' lsb binary'
      write(22,*) 'data file "lnrho_dx_sl'//char(nfile+48)//'.bin", 4'
      write(22,*) '  '
      write(22,*) 'attribute "dep" string "positions" '
      write(22,*) '  '
      write(22,*) 'object "regular positions regular connections" class field'
      write(22,*) 'component "positions" value 1'
      write(22,*) 'component "connections" value 2'
      write(22,*) 'component "data" value 3'
      write(22,*) 'end'
      close(22)
      endif

      if(mvar.ge.5) then
      open(unit=23, file='./tmp/ss_sl'//char(nfile+48)//'.dx')
      write(23,*) 'object 1 class gridpositions counts ', mx-2*nghost, &
        mz-2*nghost
      write(23,*) 'origin  ', x00, z00 
      write(23,*) 'delta  ', dx,'  0   '
      write(23,*) 'delta  ', '  0  ',dz
      write(23,*) '  '
      write(23,*) 'object 2 class gridconnections counts ', mx-2*nghost, &
        mz-2*nghost
      write(23,*) '  '
      write(23,*) 'object 3 class array type float rank 0'
      write(23,*) 'items ', (mx-2*nghost)*(mz-2*nghost), ' lsb binary'
      write(23,*) 'data file "ss_dx_sl'//char(nfile+48)//'.bin", 4'
      write(23,*) '  '
      write(23,*) 'attribute "dep" string "positions" '
      write(23,*) '  '
      write(23,*) 'object "regular positions regular connections" class field'
      write(23,*) 'component "positions" value 1'
      write(23,*) 'component "connections" value 2'
      write(23,*) 'component "data" value 3'
      write(23,*) 'end'
      close(23)
      endif


      if(mvar.ge.6) then
! Process magnetic part if mvar>5: no kinematic dynamo yet, needs to be added
! Write out 2D slice of parallel vector potential for Dx
      iax=6
      iay=7
      iaz=8
      open(unit=18, file='./tmp/aaa_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(18) (((real(f(i,isl,k,n)),n=iax,iaz,2),k=nghost+1,mz-nghost), & 
        i=nghost+1,mx-nghost)
      close(18)
      
! Calculate magnetic field, using central differences.
      do i=2,mx-1
         do j=2,my-1
            do k=2,mz-1
               bbb(i,j,k,1)=(f(i,j+1,k,iaz)-f(i,j-1,k,iaz))/2./dy-  &
       (f(i,j,k+1,iay)-f(i,j,k-1,iay))/2./dz
               bbb(i,j,k,2)=(f(i,j,k+1,iax)-f(i,j,k-1,iax))/2./dz-  &
       (f(i+1,j,k,iaz)-f(i-1,j,k,iaz))/2./dx
               bbb(i,j,k,3)=(f(i+1,j,k,iay)-f(i-1,j,k,iay))/2./dx-  &
       (f(i,j+1,k,iax)-f(i,j-1,k,iax))/2./dy               
            enddo
         enddo
      enddo   

! Write out parallel magnetic field for Dx
      open(unit=18, file='./tmp/bbb_dx_sl'//char(nfile+48)//'.bin', & 
        form='unformatted')
      write(18) (((real(bbb(i,isl,k,n)),n=1,3,2),k=nghost+1,mz-nghost), & 
        i=nghost+1,mx-nghost)
      close(18)

! Create  .dx scripts describing magnetic data 
      open(unit=24, file='./tmp/aaa_sl'//char(nfile+48)//'.dx')
      write(24,*) 'object 1 class gridpositions counts ', mx-2*nghost,  &
        mz-2*nghost
      write(24,*) 'origin  ', x00, z00 
      write(24,*) 'delta  ', dx,'  0   '
      write(24,*) 'delta  ', '  0  ',dz
      write(24,*) '  '
      write(24,*) 'object 2 class gridconnections counts ', mx-2*nghost,  &
        mz-2*nghost
      write(24,*) '  '
      write(24,*) 'object 3 class array type float rank 1 shape 2'
      write(24,*) 'items ', (mx-2*nghost)*(mz-2*nghost), ' lsb binary'
      write(24,*) 'data file "aaa_dx_sl'//char(nfile+48)//'.bin", 4'
      write(24,*) '  '
      write(24,*) 'attribute "dep" string "positions" '
      write(24,*) '  '
      write(24,*) 'object "regular positions regular connections" class field'
      write(24,*) 'component "positions" value 1'
      write(24,*) 'component "connections" value 2'
      write(24,*) 'component "data" value 3'
      write(24,*) 'end'
      close(24)

      open(unit=25, file='./tmp/bbb_sl'//char(nfile+48)//'.dx')
      write(25,*) 'object 1 class gridpositions counts ', mx-2*nghost, &
        mz-2*nghost
      write(25,*) 'origin  ', x00, z00 
      write(25,*) 'delta  ', dx,'  0   '
      write(25,*) 'delta  ', '  0  ',dz
      write(25,*) '  '
      write(25,*) 'object 2 class gridconnections counts ', mx-2*nghost, &
        mz-2*nghost
      write(25,*) '  '
      write(25,*) 'object 3 class array type float rank 1 shape 2'
      write(25,*) 'items ', (mx-2*nghost)*(mz-2*nghost), ' lsb binary'
      write(25,*) 'data file "bbb_dx_sl'//char(nfile+48)//'.bin", 4'
      write(25,*) '  '
      write(25,*) 'attribute "dep" string "positions" '
      write(25,*) '  '
      write(25,*) 'object "regular positions regular connections" class field'
      write(25,*) 'component "positions" value 1'
      write(25,*) 'component "connections" value 2'
      write(25,*) 'component "data" value 3'
      write(25,*) 'end'
      close(25)


      endif


!!!!! Finishing choosing the slices
      case default
        !
        !  Catch unknown values
        !
        print*,'No such value for sltype: ', trim(sltype)
        stop
      endselect     


      stop
      end

