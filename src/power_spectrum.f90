! $Id: power_spectrum.f90,v 1.21 2002-10-30 20:20:47 brandenb Exp $
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
!  Calculate power spectra (on shperical shells) of the variable
!  specified by `sp'.
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,im,in
  real, dimension (mx,my,mz,mvar) :: f
  real, dimension(nx,ny,nz) :: a1,b1,a2,b2,a3,b3
  real, dimension(nx,3) :: bb
  real, dimension(nk) :: spectrum=0.,spectrum_sum=0
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=1) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call cvs_id( &
       "$Id: power_spectrum.f90,v 1.21 2002-10-30 20:20:47 brandenb Exp $")
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
  if(.NOT.lfft) call stop_it('Need FFT=fft in Makefile.local to get spectra!')
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*2*pi/Lz
  !
  !  integration over shells
  !
  spectrum=0
  if(lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
  do ikz=1,nz
     do iky=1,ny
        do ikx=1,nx
           k=nint(sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2))
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
!  Calculate power and helicity spectra (on shperical shells) of the
!  variable specified by `sp', i.e. either the spectra of uu and kinetic
!  helicity, or those of bb and magnetic helicity..
!  Since this routine is only used at the end of a time step,
!  one could in principle reuse the df array for memory purposes.
!
  integer, parameter :: nk=nx/2
  integer :: i,k,ikx,iky,ikz,im,in,ivec
  real, dimension (mx,my,mz,mvar) :: f
  real, dimension(nx,ny,nz) :: a_re,a_im,b_re,b_im
  real, dimension(nx) :: bbi
  real, dimension(nk) :: spectrum=0.,spectrum_sum=0
  real, dimension(nk) :: spectrumhel=0.,spectrumhel_sum=0
  real, dimension(nxgrid) :: kx
  real, dimension(nygrid) :: ky
  real, dimension(nzgrid) :: kz
  character (len=3) :: sp
  !
  !  identify version
  !
  if (lroot .AND. ip<10) call cvs_id( &
       "$Id: power_spectrum.f90,v 1.21 2002-10-30 20:20:47 brandenb Exp $")
  !
  !    Stopping the run if FFT=nofft
  !
  if(.NOT.lfft) call stop_it('Need FFT=fft in Makefile.local to get spectra!')
  !
  !  Define wave vector, defined here for the *full* mesh.
  !  Each processor will see only part of it.
  !
  kx=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*2*pi/Lx
  ky=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*2*pi/Ly
  kz=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*2*pi/Lz
  !
  !  initialize power spectrum to zero
  !
  spectrum=0.
  spectrumhel=0.
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
          a_re(:,im,in)=bbi  !(this corresponds to vorticity)
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
          b_re(:,im,in)=bbi  !(this corresponds to magnetic field)
        enddo
      enddo
      a_re=f(l1:l2,m1:m2,n1:n2,iaa+ivec-1)
      a_im=0.
      b_im=0.
    endif
    !
    !  Doing the Fourier transform
    !
    select case (fft_switch)
    case ('fft_nr')
      call transform_nr(a_re,a_im)
      call transform_nr(b_re,b_im)
    case ('fftpack')
      call transform_fftpack(a_re,a_im)
      call transform_fftpack(b_re,b_im)
    case ('Singleton')
      call transform_i(a_re,a_im)
      call transform_i(b_re,b_im)
    case default
      call stop_it("powerhel: no fft_switch chosen")
    endselect
    !
    !  integration over shells
    !
    if(lroot .AND. ip<10) print*,'fft done; now integrate over shells...'
    do ikz=1,nz
      do iky=1,ny
        do ikx=1,nx
          k=nint(sqrt(kx(ikx)**2+ky(iky+ipy*ny)**2+kz(ikz+ipz*nz)**2))
          if(k>=0 .and. k<=(nk-1)) then
            spectrum(k+1)=spectrum(k+1) &
               +b_re(ikx,iky,ikz)**2 &
               +b_im(ikx,iky,ikz)**2
            spectrumhel(k+1)=spectrumhel(k+1) &
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
  call mpireduce_sum(spectrumhel,spectrumhel_sum,nk)
  !
  !  on root processor, write global result to file
  !  multiply by 1/2, so \int E(k) dk = (1/2) <u^2>
  !  ok for helicity, so \int F(k) dk = <o.u> = 1/2 <o*.u+o.u*>
  !
  !  append to diagnostics file
  !
  if (iproc==root) then
    if (ip<10) print*,'Writing power spectrum ',sp &
         ,' to ',trim(datadir)//'/power_'//trim(sp)//'.dat'
    !
    spectrum_sum=.5*spectrum_sum
    open(1,file=trim(datadir)//'/power_'//trim(sp)//'.dat',position='append')
    write(1,*) t
    write(1,*) spectrum_sum
    close(1)
    !
    open(1,file=trim(datadir)//'/powerhel_'//trim(sp)//'.dat',position='append')
    write(1,*) t
    write(1,*) spectrumhel_sum
    close(1)
  endif
  !
  endsubroutine powerhel
!***********************************************************************

end module power_spectrum
