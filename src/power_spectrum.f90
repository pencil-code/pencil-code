! $Id: power_spectrum.f90,v 1.12 2002-10-22 16:58:19 brandenb Exp $
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
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
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
       "$Id: power_spectrum.f90,v 1.12 2002-10-22 16:58:19 brandenb Exp $")
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
  !    Stopping the run if FFT=nofft
  !
  if (.NOT. lfft) call stop_it( 'You need FFT=fft in Makefile.local in order to get dynamical power spectrum')
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
     spectrum_sum=.5*spectrum_sum
     open(1,file=trim(datadir)//'/power'//trim(sp)//'.dat',position='append')
     write(1,*) t
     write(1,*) spectrum_sum
     close(1)
  endif
  !
  endsubroutine power
!***********************************************************************
  subroutine powerhel(f,sp)
!
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec
  real, dimension (mx,my,mz,mvar) :: f
  real, dimension(nx,ny,nz) :: a_re,a_im,b_re,b_im
  real, dimension(nx) :: bbi
  real, dimension(nk) :: spectrum=0.,spectrum_sum=0
  real, dimension(nk) :: spectrum_hel=0.,spectrum_hel_sum=0
  integer, dimension(nxgrid) :: kx
  integer, dimension(nygrid) :: ky
  integer, dimension(nzgrid) :: kz
  character (len=3) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call cvs_id( &
       "$Id: power_spectrum.f90,v 1.12 2002-10-22 16:58:19 brandenb Exp $")
  !
  !    Stopping the run if FFT=nofft
  !
  if (.NOT. lfft) call stop_it( 'You need FFT=fft in Makefile.local to get power spectra')
  !
  !  define wave vector
  !
  kx=cshift((/(i-(nxgrid-1)/2,i=0,nxgrid-1)/),+(nxgrid-1)/2)*2*pi/Lx
  ky=cshift((/(i-(nygrid-1)/2,i=0,nygrid-1)/),+(nygrid-1)/2)*2*pi/Ly
  kz=cshift((/(i-(nzgrid-1)/2,i=0,nzgrid-1)/),+(nzgrid-1)/2)*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  spectrum=0
  spectrum_hel=0
  !
  !  loop over all the components
  !
  do ivec=1,3
    !
    !  In fft, real and imaginary parts are handled separately.
    !  For "kin", calculate spectra of <uk^2> and <ok.uk>
    !  For "mag", calculate spectra of <bk^2> and <ak.bk>
    !
    if (sp=='kin') then
      do n=n1,n2
        do m=m1,m2
          call curli(f,iuu,bbi,ivec)
          im=m-nghost
          in=n-nghost
          a_re(:,im,in)=bbi  !(this corrsponds to vorticity)
        enddo
      enddo
      b_re=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
      a_im=0.
      b_im=0.
    elseif (sp=='mag') then
      do n=n1,n2
        do m=m1,m2
          call curli(f,iaa,bbi,ivec)
          im=m-nghost
          in=n-nghost
          b_re(:,im,in)=bbi  !(this corrsponds to magnetic field)
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)
      a_im=0.
      b_im=0.
    endif
    !
    !  Doing the Fourier transform
    !
    call transform_i(a_re,a_im)
    call transform_i(b_re,b_im)
    !
    !  integration over shells
    !
    if(lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k=nint(sqrt(float(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2)))
          if(k>=0 .and. k<=(nk-1)) then
            spectrum(k+1)=spectrum(k+1) &
               +b_re(ikx,iky,ikz)**2 &
               +b_im(ikx,iky,ikz)**2
            spectrum_hel(k+1)=spectrum_hel(k+1) &
               +a_re(ikx,iky,ikz)*b_re(ikx,iky,ikz) &
               +a_im(ikx,iky,ikz)*b_im(ikx,iky,ikz)
          endif
        enddo
      enddo
    enddo
    !
  enddo !(from loop over ivec)
  !
  !  Summing up the results from the different processors
  !  The result is available only on root
  !
  call mpireduce_sum(spectrum,spectrum_sum,nk)
  call mpireduce_sum(spectrum_hel,spectrum_hel_sum,nk)
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !
  !  append to diagnostics file
  !
  if (iproc==root) then
    if (ip<10) print*,'Writing power spectra of variable',sp &
         ,'to ',trim(datadir)//'/power'//trim(sp)//'.dat'
    spectrum_sum=.5*spectrum_sum
    !
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    write(1,*) t
    write(1,*) spectrum_sum
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    write(1,*) t
    write(1,*) spectrum_hel_sum
    close(1)
  endif
  !
  endsubroutine powerhel
!***********************************************************************
    subroutine powersnap(a)
!
!  Write a snapshot of power spectrum 
!
!  30-sep-97/axel: coded
!  07-oct-02/nils: adapted from wsnap
!  08-oct-02/tony: expanded file to handle 120 character datadir // '/tspec.dat'
!
      use Io
!
      real, dimension (mx,my,mz,mvar) :: a
      character (len=135) :: file
      character (len=4) :: ch
      logical lspec
      integer, save :: ifirst,nspec
      real, save :: tspec
!
!  Output snapshot with label in 'tpower' time intervals
!  file keeps the information about number and time of last snapshot
!
      file=trim(datadir)//'/tspec.dat'
!
!  at first call, need to initialize tsnap
!  tsnap calculated in out1, but only available to root processor
!
      if (ifirst==0) then
         call out1 (trim(file),tspec,nspec,dspec,t)
         ifirst=1
      endif
!
!  Check whether we want to output power snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently)
!
      call out2 (trim(file),tspec,nspec,dspec,t,lspec,ch,.false.)
      if (lspec) then
         if (vel_spec) call power(a,'u')
         if (mag_spec) call power(a,'b')
         if (vec_spec) call power(a,'a')
         if (ab_spec) call powerhel(a,'mag')
         if (ou_spec) call powerhel(a,'kin')
      endif
!
    endsubroutine powersnap
!***********************************************************************
end module power_spectrum
