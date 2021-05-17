!***********************************************************************
subroutine initialize_fourier()
!
!  Populate wavenumber arrays for fft and calculate Nyquist wavenumber.
!
!  17-May-2021/PABourdin: moved here from start.f90 and run.f90
!
      use Cdata
!
      integer :: i
!
      if (nxgrid/=1) then
        kx_fft=cshift((/(i-(nxgrid+1)/2,i=0,nxgrid-1)/),+(nxgrid+1)/2)*2*pi/Lx
        kx_fft2=kx_fft**2
        kx_nyq=nxgrid/2 * 2*pi/Lx
      else
        kx_fft=0.0
        kx_nyq=0.0
      endif
!
      if (nygrid/=1) then
        ky_fft=cshift((/(i-(nygrid+1)/2,i=0,nygrid-1)/),+(nygrid+1)/2)*2*pi/Ly
        ky_fft2=ky_fft**2
        ky_nyq=nygrid/2 * 2*pi/Ly
      else
        ky_fft=0.0
        ky_nyq=0.0
      endif
!
      if (nzgrid/=1) then
        kz_fft=cshift((/(i-(nzgrid+1)/2,i=0,nzgrid-1)/),+(nzgrid+1)/2)*2*pi/Lz
        kz_fft2=kz_fft**2
        kz_nyq=nzgrid/2 * 2*pi/Lz
      else
        kz_fft=0.0
        kz_nyq=0.0
      endif
!
    endsubroutine initialize_fourier
!***********************************************************************

