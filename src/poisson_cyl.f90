! $Id: poisson_cyl.f90,v 1.4 2008-04-23 18:53:11 wlyra Exp $

!
!  This module solves the Poisson equation in cylindrical coordinates
!    (d^2/dr^2 +1/r*d/dr + 1/r^2*d^2/dy^2 + d^2/dz^2) f = RHS(x,y,z)
!
!  Another file was coded for the module because an extra tridimensional
!  array bessel_grid is needed to store the bessel functions of the grid.
!  Calculating them at every time would be time-consuming and unecessary
!  Instead, they are calculated in start time and stored in the memory. 
!
!  For now, this file only solves the 2d cylindrical poisson equation, 
!  with Hankel transforms, following the theory outlined by Toomre (1962)
!  and in Binney and Tremaine's "Galactic Dynamics" p.74-76   
!
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpoisson=.true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Poisson

  use Cdata
  use Cparam
  use Fourier
  use Messages

  implicit none

  real :: kmax=0.0
  logical :: lrazor_thin=.true.,lsolve_bessel=.true.,lsolve_cyl2cart=.false.
  logical :: lsolve_direct=.false.

  include 'poisson.h'

  namelist /poisson_init_pars/ &
       kmax,lrazor_thin,lsolve_bessel,lsolve_cyl2cart,lsolve_direct
  namelist /poisson_run_pars/ &
       kmax,lrazor_thin,lsolve_bessel,lsolve_cyl2cart,lsolve_direct

  real, dimension(nx,ny,nx,nygrid) :: green_grid
  real, dimension(nx,nx,ny) :: bessel_grid
  real, dimension(nx) :: rad,kr_fft
  real, dimension(ny) :: tht
  integer, dimension(nygrid) :: m_fft
  integer :: nr,nth,nkr,nkt,nthgrid
  integer :: nktgrid,nkhgrid
  real :: dr,dkr,dr1,dth,dth1
  real :: r0,theta0,rn,theta1

  contains

!***********************************************************************
    subroutine initialize_poisson()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-oct-07/anders: adapted
!
      integer :: i
!
      if (coord_system/='cylindric') then
        if (lroot) print*, 'poisson_cyl: '//&
             'this module is only for cylindrical runs'
        call fatal_error('initialize_poisson','')
      endif
!
      if (lrazor_thin) then
        if (nzgrid/=1) then
          if (lroot) print*, 'initialize_poisson: '//&
               'razor-thin approximation only works with nzgrid==1'
          call fatal_error('initialize_poisson','')
        endif
      else
        if (lroot) print*, 'initialize_poisson: '//&
             'not yet implemented for 3D cylindrical runs'
        call fatal_error('initialize_poisson','')
      endif
!
      if (lsolve_cyl2cart) then 
        lsolve_bessel=.false.
        lsolve_direct=.false.
      endif
!
      if (lsolve_bessel) then 
        lsolve_cyl2cart=.false.
        lsolve_direct=.false.
      endif
!
      if (lsolve_direct) then 
        lsolve_bessel=.false.
        lsolve_cyl2cart=.false.
      endif
!
      if ((.not.lsolve_bessel).and.(.not.lsolve_cyl2cart)&
           .and.(.not.lsolve_direct)) then
        if (lroot) print*,'initialize_poisson: '//&
             'neither lsolve_bessel nor lsolve_cyl2cart nor lsolve_direct'//&
             'are switched on. Choose one of them'
        call fatal_error('initialize_poisson','')
      endif
!
! Keep the notation consistent
!
      rad=x(l1:l2)     ; tht=y(m1:m2)
      nr =nx           ; nth=ny
      r0=rad(1)        ; theta0=xyz0(2)+.5*dth 
      rn=rad(nr)       ; theta1=xyz1(2)-.5*dth
      dr=dx            ; dth=dy
      dr1=1./dr        ; dth1=1./dth
      nkr=nr           ; nkt=ny
      nktgrid=nygrid   ; nthgrid=nygrid
!
! Pre-calculate the radial wavenumbers
!           
      do i=1,nkr
        kr_fft(i)=.5*(i-1)/(nr*dr)
      enddo
      dkr=kr_fft(2)-kr_fft(1)
!
! Azimuthal wavenumbers (integers)
!
      m_fft=cshift((/(i-(nktgrid+1)/2,i=0,nktgrid-1)/),+(nktgrid+1)/2)
!
! Pre-calculate the bessel functions of the grid....
! This is VERY time consuming
!
      if (lsolve_bessel) &
           call calculate_cross_bessel_functions
!
      if (lsolve_direct) &
           call calculate_cross_green_functions
!
    endsubroutine initialize_poisson
!***********************************************************************
    subroutine inverse_laplacian(phi)
!
!  Dispatch solving the Poisson equation to inverse_laplacian_fft
!  or inverse_laplacian_semispectral, based on the boundary conditions
!
!  17-jul-2007/wolf: coded wrapper
!
      real, dimension (nx,ny,nz) :: phi
!
      intent(inout) :: phi
!
      if (lsolve_bessel) then
        call inverse_laplacian_bessel(phi)
      else if (lsolve_cyl2cart) then
        call inverse_laplacian_cyl2cart(phi)
      else if (lsolve_direct) then
        call inverse_laplacian_directsum(phi)
      else
        call fatal_error("inverse_laplacian","no solving method given")
      endif
!
    endsubroutine inverse_laplacian
!***********************************************************************
    subroutine inverse_laplacian_cyl2cart(phi)
!
!  Solve the 2D Poisson equation in cylindrical coordinates
!  by transforming to a periodic cartesian grid before 
!  Fourier transforming. 
!
!  This subroutine is faster than inverse_laplacian_bessel
!  for low resolutions and low number of processors.
!  But it gets slower when increasing any of them, 
!  due to the great amount of communication used. 
!
!  The frequent broadcast of a big array gave problems at the 
!  PIA cluster in Heidelberg after some thousands of time-steps. 
!  The problem was probably due to memory leaking. No problem occured 
!  at the UPPMAX cluster in Uppsala. So beware that using this broadcasting
!  extravaganza subroutine might not work in all clusters. 
!
!  01-12-07/wlad: coded
!  28-02-08/wlad: merged the serial and mpi versions
!
      use Mpicomm
!
      real, dimension (nx,ny,nz)  :: phi
      real, dimension (2*nx,2*ny) :: nphi,nb1
!
      real, dimension(nygrid) :: theta_serial
!
      real, dimension(2*nx)     :: xc,kkx_fft
      real, dimension(2*ny)     :: yc
      real, dimension(2*nygrid) :: kky_fft,yserial
!
! For communication
!
      real, dimension (nx,ny)         :: cross_proc
      real, dimension (2*nx,2*ny)     :: cross_proc_big
      real, dimension (nx,nygrid)     :: cross
      real, dimension (2*nx,2*nygrid) :: crossbig
!
! Cheap stuff
!
      real    :: x0,xn,y0,yn,dxc,dyc,dxc1,dyc1,Lxn,Lyn
      real    :: theta,radius,k2,xp,yp
      real    :: delr,delp,fr,fp,delx,dely,fx,fy
      real    :: p1,p2,p3,p4,interp_pot
!
      integer :: ix1,ix2,iy1,iy2,ir1,ir2,ip1,ip2
      integer :: i,j,ikx,iky,ir,im,ido,iup,ith
      integer :: nnx,nny,nnghost
      integer :: nnxgrid,nnygrid
!
      if (nx/=nygrid) &
           call fatal_error("","currently only works for nx=nygrid")
      if (nzgrid/=1)  &
           call fatal_error("","currently only works for 2D simulations")
!
! Expanded cartesian axes
!
      nnx=2*nx         ; nny=2*ny
      xn=2*rad(nr)     ; x0=-xn
      yn=xn            ; y0=-yn
      nnxgrid=2*nxgrid ; nnygrid=2*nygrid
!
      do i=1,nnx
        xc(i)=1.*(i-1)        /(nnxgrid-1)*(xn-x0)+x0
      enddo
      do m=1,nny
        yc(m)=1.*(m-1+ipy*nny)/(nnygrid-1)*(yn-y0)+y0
      enddo      
! 
      dxc=xc(2)-xc(1)  ; dyc=dxc
      dxc1=1/dxc       ; dyc1=1/dyc
!
! Now transform to Cartesian grid
!
      nnghost=npoint-nghost
!
      if (lmpicomm) then
!
! All processors send its density array to the root processor
!
        if (.not.lroot) then
          call mpisend_real(phi(:,:,nnghost),(/nx,ny/),root,111)
        else
          cross_proc=phi(:,:,nnghost)
!
! The root processor receives all arrays and
! stores them in a single big array of dimension
! nx*nygrid
!
          do j=0,ncpus-1
            if (j/=0) call mpirecv_real(cross_proc,(/nx,ny/),j,111)
            ido= j  * ny + 1
            iup=(j+1)*ny
            cross(:,ido:iup)=cross_proc
          enddo
        endif
!
! Broadcast the density field to all processors
!
        call mpibcast_real(cross,(/nx,nygrid/))
!
      else
!
! For serial runs, ny=nygrid, so just copy the density
!
        cross(:,1:ny)=phi(:,1:ny,nnghost)
!
      endif
!
! Need the serial theta later in order to compute the
! azimuthal displacement in parallel
!
      do i=1,nygrid
        theta_serial(i)=1.*(i-1)/(nthgrid-1)*(theta1-theta0)+theta0
      enddo
!
! Now transform the grid just like we would do in a serial run
!
      do m=1,nny
        do i=1,nnx
          radius=sqrt(xc(i)**2+yc(m)**2)
          if ((radius.ge.r0).and.(radius.le.rn)) then
            ir1=floor((radius-r0)*dr1 ) +1;ir2=ir1+1
            delr=radius-rad(ir1)
!
! this should never happen, but is here for warning
!
            if (ir1.lt.1 ) call fatal_error("","cyl2cart: ir1<1")
            if (ir2.gt.nr) call fatal_error("","cyl2cart: ir2>nr")
!
            theta=atan2(yc(m),xc(i))
            ip1=floor((theta - theta0)*dth1)+1;ip2=ip1+1
            if (ip1==0) then
              ip1=nthgrid
              delp=theta-theta_serial(ip1) + 2*pi
            else
              delp=theta-theta_serial(ip1)
            endif
            if (ip2==nthgrid+1) ip2=1
!
! Bilinear interpolation
!
            !p1=phi(ir1,ip1,nnghost);p2=phi(ir2,ip1,nnghost)
            !p3=phi(ir1,ip2,nnghost);p4=phi(ir2,ip2,nnghost)
            p1=cross(ir1,ip1);p2=cross(ir2,ip1)
            p3=cross(ir1,ip2);p4=cross(ir2,ip2)
!
            fr=delr*dr1
            fp=delp*dth1
!
            nphi(i,m)=fr*fp*(p1-p2-p3+p4) + fr*(p2-p1) + fp*(p3-p1) + p1
          else
            nphi(i,m)=0.
          endif
        enddo
      enddo
!
!  The right-hand-side of the Poisson equation is purely real.
!
      nb1=0.0
!
!  Forward transform (to k-space).
!
      call fourier_transform_xy_xy_other(nphi,nb1)
!
!  Solve Poisson equation
!
      Lxn=2*xc(nnx);Lyn=Lxn
!
      kkx_fft=cshift((/(i-(nnxgrid+1)/2,i=0,nnxgrid-1)/),+(nnxgrid+1)/2)*2*pi/Lxn
      kky_fft=cshift((/(i-(nnygrid+1)/2,i=0,nnygrid-1)/),+(nnygrid+1)/2)*2*pi/Lyn
!
      do iky=1,nny
        do ikx=1,nnx
          if ((kkx_fft(ikx)==0.0) .and. (kky_fft(iky+ipy*nny)==0.0)) then
            nphi(ikx,iky) = 0.0
            nb1(ikx,iky) = 0.0
          else
            if (.not.lrazor_thin) then
              call fatal_error("","3d case not implemented yet")
!
!  Razor-thin approximation. Here we solve the equation
!    del2Phi=4*pi*G*Sigma(x,y)*delta(z)
!  The solution at scale k=(kx,ky) is
!    Phi(x,y,z)=-(2*pi*G/|k|)*Sigma(x,y)*exp[i*(kx*x+ky*y)-|k|*|z|]
!
            else
              k2 = (kkx_fft(ikx)**2+kky_fft(iky+ipy*nny)**2)
              nphi(ikx,iky) = -.5*nphi(ikx,iky) / sqrt(k2)
              nb1(ikx,iky)  = -.5*nb1(ikx,iky)  / sqrt(k2)
            endif
          endif
!
!  Limit |k| < kmax
!
          if (kmax>0.0) then
            if (sqrt(k2)>=kmax) then
              nphi(ikx,iky) = 0.
              nb1(ikx,iky) = 0.
            endif
          endif
        enddo
      enddo
!
!  Inverse transform (to real space).
!
      call fourier_transform_xy_xy_other(nphi,nb1,linv=.true.)
!
      if (lmpicomm) then
        if (.not.lroot) then
          call mpisend_real(nphi,(/nnx,nny/),root,222)
        else
          cross_proc_big=nphi
          !the root processor receives all arrays and
          !stores them in a single big array of dimension
          !nx*nygrid
          do j=0,ncpus-1
            if (j/=0) call mpirecv_real(cross_proc_big,(/nnx,nny/),j,222)
            ido= j   *nny+1
            iup=(j+1)*nny
            crossbig(:,ido:iup)=cross_proc_big
          enddo
        endif
        call mpibcast_real(crossbig,(/nnx,nnygrid/))
      else
        crossbig(:,1:nny)=nphi(:,1:nny)
      endif
!
!  Convert back to cylindrical
!
      yserial(1:nygrid)=xc(1:nygrid)
      do ith=1,Nth
        do ir=1,Nr
!
          xp=rad(ir)*cos(tht(ith))
          yp=rad(ir)*sin(tht(ith))
!
          ix1 = floor((xp-x0)*dxc1)+1 ; ix2 = ix1+1
          iy1 = floor((yp-y0)*dyc1)+1 ; iy2 = iy1+1
!
          if (ix1 .lt.  1)      call fatal_error("","ix1 lt 1")
          if (iy1 .lt.  1)      call fatal_error("","iy1 lt 1")
          if (ix2 .gt. nnxgrid) call fatal_error("","ix2 gt nnxgrid")
          if (iy2 .gt. nnygrid) call fatal_error("","iy2 gt nnygrid")
!
          delx=xp-     xc(ix1);fx=delx*dxc1
          dely=yp-yserial(iy1);fy=dely*dyc1
!
! Bilinear interpolation
!
          p1=crossbig(ix1,iy1);p2=crossbig(ix2,iy1)
          p3=crossbig(ix1,iy2);p4=crossbig(ix2,iy2)
!
          interp_pot=fx*fy*(p1-p2-p3+p4) + fx*(p2-p1) + fy*(p3-p1) + p1
!
          do n=1,nz
            phi(ir,ith,n)=interp_pot
          enddo
!
        enddo
      enddo
!
    endsubroutine inverse_laplacian_cyl2cart
!***********************************************************************
    subroutine inverse_laplacian_bessel(phi)
!
!  Solve the 2D Poisson equation in cylindrical coordinates
!
!  This beautiful and elegant theory for calculating the 
!  potential of thin disks using bessel functions and hankel 
!  transforms is not so useful numerically because of the 
!  ammount of integrations to perform. 
!
!  In any case, we managed to optimize the computations so that 
!  only 30% of the computational time is spent in this routine. 
!  The number is only slightly dependant on parallelization and 
!  resolution, as opposed to inverse_laplacian_cyl2cart.
!  
!  For example, having this 
!
!   do ikr=1,nkr
!     SS(ikr)=sum(bessel_grid(2:nr-1,ikr,ikt)*sigma_tilde_rad(2:nr-1))+&
!          .5*(bessel_grid( 1,ikr,ikt)*sigma_tilde_rad(1)             +&
!              bessel_grid(nr,ikr,ikt)*sigma_tilde_rad(nr))
!   enddo
!
!  instead of 
!
!    do ikr=1,nkr
!      tmp=bessel_grid(:,ikr,ikt)*sigma_tilde*rad
!      SS(ikr)=sum(tmp(2:nr-1))+.5*(tmp(1)+tmp(nr))
!    enddo
!
!  lead to a factor 2 speed up. That's a lot for a subroutine
!  that usually takes most of the computing time. So, think (twice) 
!  before you modify it.
!
!  At every wavelength, the density-potential pair in fourier space is
!
!   Phi_tilde_k = exp(-k|z|) Jm(k*r) ; Sigma_tilde_k =-k/(2piG) Jm(k*r)  
!
!  So, the total potential and the total spectral density is 
!
!   Phi_tilde  =          Int [ S(k) Jm(k*r) exp(-k|z|) ]dk  
!
!   Sigma_tilde=-1/(2piG) Int [ S(k) Jm(k*r) k] dk
! 
!  The function above is the Hankel transform of S(k), so S(k) 
!  is the inverse Hankel transform of Sigma_tilde
!
!   S(k)=-2piG Int[ Jm(k*r) Sigma_tilde(r,m) r] dr 
!
!  06-03-08/wlad: coded
!
      use Mpicomm
!
      real, dimension (nx,ny,nz)  :: phi,b1
      complex, dimension(nx) :: SS,sigma_tilde_rad,tmp,tmp2
      integer :: i,ir,ikr,ikt,nnghost
      real :: fac
!
      if (nx/=nygrid) &
           call fatal_error("","currently only works for nx=nygrid")
      if (nzgrid/=1)  &
           call fatal_error("","currently only works for 2D simulations")
      nnghost=npoint-nghost
!
! Fourier transform in theta 
! 
      call transp(phi,'y'); b1=0.
      call fourier_transform_x(phi,b1)
      call transp(phi,'y');call transp(b1,'y')
!
! SS is the hankel transform of the density
!
      fac=-.5*dr*dkr !actually, -2*pi*G = -.5*rhs_poisson_const
      do ikt=1,nkt 
!       
! Hankel transform of the fourier-transformed density
!
        sigma_tilde_rad=rad*cmplx(phi(:,ikt,nnghost),b1(:,ikt,nnghost))
        do ikr=1,nkr
          SS(ikr)=sum(bessel_grid(2:nr-1,ikr,ikt)*sigma_tilde_rad(2:nr-1))+&
               .5*(bessel_grid( 1,ikr,ikt)*sigma_tilde_rad(1)             +&
                   bessel_grid(nr,ikr,ikt)*sigma_tilde_rad(nr))
        enddo
!
! Sum up all the contributions to Phi
!
        do ir=1,nr
          tmp(ir)=sum(bessel_grid(ir,2:nkr-1,ikt)*SS(2:nkr-1))+&
               .5*(bessel_grid(ir,  1,ikt)*SS(1)+             &
                   bessel_grid(ir,nkr,ikt)*SS(nkr))
        enddo
        do n=1,nz
          phi(:,ikt,n)=real(tmp)
           b1(:,ikt,n)=aimag(tmp)
        enddo
      enddo
!
! Transform back to real space
!
      call transp(phi,'y');call transp(b1,'y')
      call fourier_transform_x(phi,b1,linv=.true.)
      call transp(phi,'y')
!
      phi=phi*fac
!
    endsubroutine inverse_laplacian_bessel
!***********************************************************************
    subroutine inverse_laplacian_directsum(phi)
!
!  Solve the 2D Poisson equation in cylindrical coordinates
!
!  Direct summation phi=-G Int(rho/|r-r'| dV)
!
!  23-04-08/wlad: coded
!
      use Mpicomm
!
      real, dimension (nx,ny,nz)  :: phi
      real, dimension (nx,nygrid) :: cross,integrand
      real, dimension (nx,ny)     :: tmp,cross_proc
      real, dimension (nx)        :: intr
      integer :: ith,ir,ikr,ikt,imn
      integer :: nnghost,i,j,ido,iup
      real :: fac
!
      if (nzgrid/=1)  &
           call fatal_error("","currently only works for 2D simulations")
      nnghost=npoint-nghost

      !transfer fac to the green grid
      !fac=-.25*pi_1 !actually, -G = rhs_poisson_const/4*pi

      if (lmpicomm) then
!
! All processors send its density array to the root processor
!
        if (.not.lroot) then
          call mpisend_real(phi(:,:,nnghost),(/nx,ny/),root,111)
        else
          cross_proc=phi(:,:,nnghost)
!
! The root processor receives all arrays and
! stores them in a single big array of dimension
! nx*nygrid
!
          do j=0,ncpus-1
            if (j/=0) call mpirecv_real(cross_proc,(/nx,ny/),j,111)
            ido= j  * ny + 1
            iup=(j+1)*ny
            cross(:,ido:iup)=cross_proc
          enddo
        endif
!
! Broadcast the density field to all processors
!
        call mpibcast_real(cross,(/nx,nygrid/))
!
      else
!
! For serial runs, ny=nygrid, so just copy the density
!
        cross(:,1:ny)=phi(:,1:ny,nnghost)
!
      endif
!
! Now integrate through direct summation
!
      do ir=1,nr 
      do ith=1,nth
!
! the scaled green function already has the r*dr*dth term
!
        integrand=cross*green_grid(ir,ith,:,:)
!
        phi(ir,ith,1:nz)=sum(                    &
             sum(integrand(2:nr-1,:)) +          &
             .5*(integrand(1,:)+integrand(nr,:)) &
                       )
!
      enddo
      enddo
!
    endsubroutine inverse_laplacian_directsum
!***********************************************************************
    subroutine calculate_cross_bessel_functions
!
!  Calculate the Bessel functions related to the 
!  cylindrical grid
!
!     bessel_grid(ir,ikr,m)=J_m(kr(ikr)*rad(ir))
!
!  06-03-08/wlad: coded
!
      use General,only: besselj_nu_int
!
      real    :: tmp,arg
      integer :: ir,ikr,ikt
!
      if (lroot) &
        print*,'Pre-calculating the Bessel functions to '//&
        'solve the Poisson equation'
!
      do ikt=1,nkt
        do ir=1,nr;do ikr=1,nkr
          arg=kr_fft(ikr)*rad(ir)
          call besselj_nu_int(tmp,m_fft(ikt+ipy*nkt),arg)
          bessel_grid(ir,ikr,ikt)=tmp
        enddo;enddo
      enddo
!
      if (lroot) &
           print*,'calculated all the needed functions'
!
    endsubroutine calculate_cross_bessel_functions
!***********************************************************************
    subroutine calculate_cross_green_functions
!
!  Calculate the Green functions related to the 
!  cylindrical grid (modified by the jacobian, the 
!  gravitaional costant and the grid elements to 
!  ease the amount of calculations in runtime)
!
!     green_grid(ir,ip,ir',ip')=-G./|r-r'| * r*dr*dth
!
!  06-03-08/wlad: coded
!
      use Mpicomm
!
      real, dimension(nygrid) :: tht_serial
      real, dimension(ny) :: tht_proc
      real    :: jacobian,tmp,Delta,fac
      integer :: ir,ith,ikr,ikt,i,ido,iup,j
!
      if (lroot) &
        print*,'Pre-calculating the Green functions to '//&
        'solve the Poisson equation'
!
      fac=-.25*pi_1 ! -G = -rhs_poisson_const/4*pi
!
! Serial theta to compute the azimuthal displacement in parallel
!
      if (lmpicomm) then
!
! All processors send its density array to the root processor
!
        if (.not.lroot) then
          call mpisend_real(tht,ny,root,111)
        else
!
          tht_proc=tht
!
! The root processor receives all arrays and
! stores them in a single array of dimension nygrid
!
          do j=0,ncpus-1
            if (j/=0) call mpirecv_real(tht_proc,ny,j,111)
            ido= j  * ny + 1
            iup=(j+1)*ny
            tht_serial(ido:iup)=tht_proc
          enddo
        endif
!
! Broadcast the serial thetas to all processors
!
        call mpibcast_real(tht_serial,nygrid)
!
      else
!
! For serial runs, ny=nygrid, so just copy the array
!
        tht_serial(1:ny)=tht(1:ny)
!
      endif
!
! Define the smoothing length as the minimum resolution element present
!
      Delta=min(dr,dth)
!
      do ir =1,nr;do ith=1,nth
      do ikr=1,nr;do ikt=1,nthgrid
!
        jacobian=rad(ikr)*dr*dth
!
        tmp=sqrt(Delta**2 + rad(ir)**2 + rad(ikr)**2 - &
             2*rad(ir)*rad(ikr)*cos(tht(ith)-tht_serial(ikt)))
!
        green_grid(ir,ith,ikr,ikt)= fac*jacobian/tmp
!
      enddo;enddo
      enddo;enddo
!
      if (lroot) &
           print*,'calculated all the needed green functions'
!
    endsubroutine calculate_cross_green_functions
!***********************************************************************
    subroutine read_poisson_init_pars(unit,iostat)
!
!  Read Poisson init parameters.
!
!  17-oct-2007/anders: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=poisson_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=poisson_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_poisson_init_pars
!***********************************************************************
    subroutine write_poisson_init_pars(unit)
!
!  Write Poisson init parameters.
!
!  17-oct-2007/anders: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=poisson_init_pars)
!
    endsubroutine write_poisson_init_pars
!***********************************************************************
    subroutine read_poisson_run_pars(unit,iostat)
!
!  Read Poisson run parameters.
!
!  17-oct-2007/anders: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=poisson_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=poisson_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_Poisson_run_pars
!***********************************************************************
    subroutine write_poisson_run_pars(unit)
!
!  Write Poisson run parameters.
!
!  17-oct-2007/anders: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=poisson_run_pars)
!
    endsubroutine write_poisson_run_pars
!***********************************************************************
endmodule Poisson
