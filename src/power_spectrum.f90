! $Id: power_spectrum.f90,v 1.5 2002-10-02 20:11:14 dobler Exp $
!
!  reads in full snapshot and calculates power spetrum of u
!
!-----------------------------------------------------------------------
!   3-sep-02/axel+nils: coded
!   5-sep-02/axel: loop first over all points, then distribute to k-shells
!   23-sep-02/nils: adapted from postproc/src/power_spectrum.f90
!

module  power_spectrum
  !
  use Cdata
  use General
  use Mpicomm
  use Sub
  !
  implicit none
  !
  contains

!***********************************************************************
    subroutine power(f,sp)
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,im,in
  real, dimension (mx,my,mz,mvar) :: f
  real, dimension(nx,ny,nz) :: a1,b1,a2,b2,a3,b3
  real, dimension(nx,3) :: bb
  real, dimension(nk) :: spectrum=0.,spectrum_sum=0
  integer, dimension(nxgrid) :: kx
  integer, dimension(nygrid) :: ky
  integer, dimension(nzgrid) :: kz
  character (len=1) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call cvs_id( &
       "$Id: power_spectrum.f90,v 1.5 2002-10-02 20:11:14 dobler Exp $")
  !
  !  In fft, real and imaginary parts are handled separately.
  !  Initialize real part a1-a3; and put imaginary part, b1-b3, to zero
  !
  if (sp=='u') then
     a1=f(l1:l2,m1:m2,n1:n2,iux)
     a2=f(l1:l2,m1:m2,n1:n2,iuy)
     a3=f(l1:l2,m1:m2,n1:n2,iuz)
  elseif (sp=='b') then
     do n=n1,n2
        do m=m1,m2
           call curl(f,iaa,bb)
           im=m-nghost
           in=n-nghost
           a1(:,im,in)=bb(:,1)
           a2(:,im,in)=bb(:,2)
           a3(:,im,in)=bb(:,3)
        enddo
     enddo
  elseif (sp=='a') then
     a1=f(l1:l2,m1:m2,n1:n2,iax)
     a2=f(l1:l2,m1:m2,n1:n2,iay)
     a3=f(l1:l2,m1:m2,n1:n2,iaz)
  else
     print*,'There are no such sp=',sp
  endif
  b1=0
  b2=0
  b3=0
  !
  !  Doing the Fourier transform
  !
  call transform(a1,a2,a3,b1,b2,b3)
  !
  !  define wave vector
  !
  kx=cshift((/(i-(nxgrid-1)/2,i=0,nxgrid-1)/),+(nxgrid-1)/2)*2*pi/Lx
  ky=cshift((/(i-(nygrid-1)/2,i=0,nygrid-1)/),+(nygrid-1)/2)*2*pi/Ly
  kz=cshift((/(i-(nzgrid-1)/2,i=0,nzgrid-1)/),+(nzgrid-1)/2)*2*pi/Lz
  !
  !  integration over shells
  !
  spectrum=0
  if(lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
  do ikz=1,nz
     do iky=1,ny
        do ikx=1,nx
           k=nint(sqrt(float(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2)))
           if(k>=0 .and. k<=(nk-1)) spectrum(k+1)=spectrum(k+1) &
                +a1(ikx,iky,ikz)**2+b1(ikx,iky,ikz)**2 &
                +a2(ikx,iky,ikz)**2+b2(ikx,iky,ikz)**2 &
                +a3(ikx,iky,ikz)**2+b3(ikx,iky,ikz)**2
        enddo
     enddo
  enddo
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !
!
!  append to diagnostics file
!
  if (iproc==root) then
     if (ip<10) print*,'Writing power spectra of variable',sp &
          ,'to ',trim(datadir)//'/power'//trim(sp)//'.dat'
     open(1,file=trim(datadir)//'/power'//trim(sp)//'.dat',position='append')
     write(1,*) spectrum_sum
     close(1)
  endif
  !
end subroutine power
!***********************************************************************
end module power_spectrum


