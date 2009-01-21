! $Id$

!***********************************************************************
      program read_xyaver_alpha
!
!  read xy averages
!
!  27-aug-07/axel: coded
!
      implicit none
!
      integer :: i, it, nt=1000000, iz
      integer, parameter :: nzgrid=16, njtest=2
      real, parameter :: pi=3.14159265358979324D0
      real, dimension (nzgrid,3,njtest) :: Eijqz
      real, dimension (nzgrid,2,2) :: alpijz1, etaijkz1
      real, dimension (nzgrid) :: z, cz, sz
      real :: t, dz,alp11,alp21,eta11,eta21
!
!  prepare matrix data
!
      dz=2.*pi/nzgrid
      do iz=1,nzgrid; z(iz)=(iz-.5-nzgrid/2)*dz; enddo
      cz=cos(z)
      sz=sin(z)
!
!  open output file
!
      open(2,file='alphaeta.dat')
!
!  read until end of file
!
      open(1,file='data/xyaverages.dat')
      do it=1,nt
        read(1,'(1pe12.5)',end=999) t
        read(1,'(1p,8e13.5)',end=999) Eijqz
!
!  new
!
!     do i=1,3
!        alpijz1(:,i,1)=+cz*Eijqz(:,i,1)+sz*Eijqz(:,i,2)
!       etaijkz1(:,i,1)=-sz*Eijqz(:,i,1)+cz*Eijqz(:,i,2)
!     enddo
!
!  old
!
      do i=1,3
         alpijz1(:,i,1)=sz*Eijqz(:,i,1)+cz*Eijqz(:,i,2)
        etaijkz1(:,i,1)=cz*Eijqz(:,i,1)-sz*Eijqz(:,i,2)
      enddo
!
      alp11=sum(alpijz1(:,1,1))/nzgrid
      alp21=sum(alpijz1(:,2,1))/nzgrid
      eta11=sum(etaijkz1(:,1,1))/nzgrid
      eta21=sum(etaijkz1(:,2,1))/nzgrid
      write(2,*) t,alp11,alp21,eta11,eta21
!
      enddo
      close(1)
      close(2)
!
999   continue
      print*,'finished reading, t=',t
      end
