! $Id$

!
!  This module solves the Poisson equation in cylindrical coordinates
!    (d^2/dr^2 +1/r*d/dr + 1/r^2*d^2/dy^2 + d^2/dz^2) f = RHS(x,y,z)
!
!  Another method was coded for the module because an extra tridimensional
!  array bessel_grid is needed to store the bessel functions of the grid.
!  Calculating them at every time would be time-consuming and unecessary
!  Instead, they are calculated in start time and stored in the memory. 
!
!  For now, this method only solves the 2d cylindrical poisson equation, 
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

  real :: kmax=0.0,iteration_threshold=1e-3,rjac=1.
  logical :: lrazor_thin=.true.,lsolve_bessel=.false.,lsolve_cyl2cart=.false.
  logical :: lsolve_relax_sor=.false.
  character (len=labellen) :: ipoisson_method='nothing'

  integer, parameter :: mmax=8 !eight harmonics for the azimuthal direction

  include 'poisson.h'

  namelist /poisson_init_pars/ &
       kmax,lrazor_thin,ipoisson_method,iteration_threshold
  namelist /poisson_run_pars/ &
       kmax,lrazor_thin,ipoisson_method,iteration_threshold
!
! For the Bessel and Hankel transforms, the grid of Bessel functions
!
  real, dimension(nx,nx,ny) :: bessel_grid
!
! For the 3D mesh relaxation, the functions for the border
!
  real, dimension(nx,nz,nxgrid,nzgrid,0:mmax) :: Legendre_Qmod
  real, dimension(ny,nygrid,0:mmax) :: fourier_cosine_terms
!
! For correcting self-acceleration
!
  real, dimension(nx,ny,nz) :: phi_previous_step,rhs_previous_step
  real, dimension(nx) :: rad,kr_fft,sqrtrad_1,rad1
  real, dimension(ny) :: tht
  real, dimension(nz) :: zed
  integer, dimension(nygrid) :: m_fft
  integer :: nr,nth,nkr,nkt,nthgrid
  integer :: nktgrid,nkhgrid
  real :: dr,dkr,dr1,dth,dth1,dz1
  real :: r0,theta0,rn,theta1

  contains

!***********************************************************************
    subroutine initialize_poisson
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
      select case (ipoisson_method)

      case ('bessel')
        if (lroot) print*,'Selecting the cylindrical '//&
             'Poisson solver that employs Bessel functions'
        lsolve_bessel    =.true.

      case ('cyl2cart')
        if (lroot) print*,'Selecting the cylindrical '//&
             'Poisson solver that transforms to a periodic '//&
             'Cartesian grid and applies Fourier transforms there'
        lsolve_cyl2cart  =.true.

      case ('sor')
        if (lroot) print*,'Selecting the cylindrical '//&
             'Poisson solver that performs mesh-relaxation with SOR'
        lsolve_relax_sor=.true.

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'initialize_poisson: '//&
             'No such value for ipoisson_method: ',&
             trim(ipoisson_method)
        call fatal_error('initialize_poisson','')
!
      endselect
!
      if (lrazor_thin) then
        if (nzgrid/=1) then
          if (lroot) print*, 'initialize_poisson: '//&
               'razor-thin approximation only works with nzgrid==1'
          call fatal_error('initialize_poisson','')
        endif
      else
        if (.not.lsolve_relax_sor) then 
          if (lroot) print*, 'initialize_poisson: '//&
               'not yet implemented for 3D cylindrical runs'
          call fatal_error('initialize_poisson','')
        endif
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
! For the 3D with SOR
!
      zed=z(n1:n2)     ; dz1=1./dz
      sqrtrad_1=1./sqrt(rad)
      rad1 = 1./rad
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
! This ones below are VERY time consuming
! Pre-calculate the special functions of the grid....
!
      if (lsolve_bessel) &
           call calculate_cross_bessel_functions
!
      if (lsolve_relax_sor) &
        call calculate_cross_legendre_functions
!
!
! Spectral radius and omega. This spectral radius is 
! only for an equally spaced square grid. 
!
      rjac=(cos(pi/nr) + (dr/dz)**2*cos(pi/nz))/(1+(dr/dz)**2)
!
    endsubroutine initialize_poisson
!***********************************************************************
    subroutine inverse_laplacian(f,phi)
!
!  Dispatch solving the Poisson equation to inverse_laplacian_fft
!  or inverse_laplacian_semispectral, based on the boundary conditions
!
!  17-jul-2007/wolf: coded wrapper
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz) :: phi
!
      intent(inout) :: phi
!
      if (lsolve_bessel) then
        call inverse_laplacian_bessel(phi)
      else if (lsolve_cyl2cart) then
        call inverse_laplacian_cyl2cart(phi)
      else if (lsolve_relax_sor) then
        call inverse_laplacian_sor(f,phi)
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
           call fatal_error("inverse_laplacian_cyl2cart","currently only works for nx=nygrid")
      if (nzgrid/=1)  &
           call fatal_error("inverse_laplacian_cyl2cart","currently only works for 2D simulations")
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
          if ((radius>=r0).and.(radius<=rn)) then
            ir1=floor((radius-r0)*dr1 ) +1;ir2=ir1+1
            delr=radius-rad(ir1)
!
! this should never happen, but is here for warning
!
            if (ir1<1 ) call fatal_error("cyl2cart","ir1<1")
            if (ir2>nr) call fatal_error("cyl2cart","ir2>nr")
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
              call fatal_error("inverse_laplacian_cyl2cart","3d case not implemented yet")
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
          if (ix1 <  1)      call fatal_error("cyl2cart","ix1 lt 1")
          if (iy1 <  1)      call fatal_error("cyl2cart","iy1 lt 1")
          if (ix2 > nnxgrid) call fatal_error("cyl2cart","ix2 gt nnxgrid")
          if (iy2 > nnygrid) call fatal_error("cyl2cart","iy2 gt nnygrid")
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
      if (nzgrid/=1)  &
           call fatal_error("inverse_laplacian_bessel",&
           "currently only works for 2D simulations")
!
      nnghost=npoint-nghost
!
! Fast Fourier transform in theta
! 
      b1=0.
      call fourier_transform_y(phi,b1)
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
      call fourier_transform_y(phi,b1,linv=.true.)
!
      phi=phi*fac
!
    endsubroutine inverse_laplacian_bessel
!***********************************************************************
    subroutine inverse_laplacian_sor(f,phi)
!
!  Solve the 3D Poisson equation in cylindrical coordinates by 
!  using mesh-relaxation with the SOR algorithm. The borders have
!  to be pre-computed. For this, we use the analytical expression
!  of Cohl and Tohline (1999) using harmonic expansions in Legendre 
!  functions.
! 
!  This is very expensive and does not need to be done every time. Just 
!  when the potential has changed by some small threshold
!
!  23-04-08/wlad: coded
!
      use Mpicomm
      use Sub, only: del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny,nz)  :: phi,rhs,b1,b1_rhs
      real, dimension (nx) :: norm0,norm
      real, dimension (nx) :: a_band,b_band,c_band,d_band
      real, dimension (nx) :: e_band,b_band1,del2phi
      integer :: ir,ith,iz,im,i
      logical :: lfirst_timestep
      logical, dimension(nx,nz) :: lupdate
      logical, dimension(nx,ny,nz) :: lupdate_grid
      logical, save :: lfirstcall=.true.
      logical :: lverbose
!
      lverbose=lroot.and.(lfirstcall.or.(ip<=8))
!
      if (lverbose) print*,'Solving the Poisson equation with SOR'
!
      if (nzgrid==1)  &
           call fatal_error("inverse_laplacian_sor","This method uses the "//&
           "discretized Poisson equation. It cannot be used for 2d in the r-phi plane")
!
      if (nprocz/=1) &
           call fatal_error("inverse_laplacian_sor","Not yet implemented for z "//&
           "parallelization. Put all processors in the azimuthal direction")
!
! The other poisson solvers get the rhs as phi, and overwrite it.
! This poisson solver needs to store the rhs first.
!     
      rhs=phi
! 
! Upon re-starting, get the stored potential
!
      if ((it==1).and.(t/=0)) &
           phi_previous_step=f(l1:l2,m1:m2,n1:n2,ipotself)
!
! Get the logical for updating the grid. Check if the quantity 
! del2phi changed significantly compared to rhs. To avoid non-convergence
! due to discretization errors, compare the rhs from the previous time
! with the one from this timestep. 
!
      if (lfirstcall) then
        lupdate_grid=.true.
      else
        do m=m1,m2;do n=n1,n2
          ith=m-nghost;iz=n-nghost
          call del2(f,ipotself,del2phi)
          norm0=abs(del2phi - rhs_previous_step(:,ith,iz))
          norm=abs(del2phi - rhs(:,ith,iz))
          do i=l1,l2
            ir=i-nghost
            if (abs(norm(ir)-norm0(ir))/norm(ir) < 1e-3) then 
              lupdate_grid(ir,ith,iz)=.false.  
            else
              lupdate_grid(ir,ith,iz)=.true.  
            endif
          enddo
        enddo;enddo
      endif
!     
! Save the rhs from the previous time step
!   
      rhs_previous_step=rhs
!
! Start the solver. Get the borders
!
      if (lverbose) print*,'calculating the border'
      call get_border_values(rhs,phi,lupdate_grid,lverbose)
      if (lverbose) print*,'done calculating the border'
!
! For the other time-steps, the potential is known, so use discretization 
! of the 5-point formula (in fourier space) to determine the potential
!
      if (.not.lfirstcall) then
!
! Get the phi from the previous time step, apart from the newly 
! updated boundaries
!
        if (nprocz >= 2) then 
          if (lfirst_proc_z) then 
            phi(2:nr-1,:,2:nz)=phi_previous_step(2:nr-1,:,2:nz)
          elseif (llast_proc_z) then 
            phi(2:nr-1,:,1:nz-1)=phi_previous_step(2:nr-1,:,1:nz-1)
          else
            phi(2:nr-1,:,:)=phi_previous_step(2:nr-1,:,:)
          endif
        else
          phi(2:nr-1,:,2:nz-1)=phi_previous_step(2:nr-1,:,2:nz-1)
        endif
!
      endif
!
! Fourier transform the potential and the right-hand-side
!
      call fourier_transform_y(phi,b1)
      call fourier_transform_y(rhs,b1_rhs)
!
      if (lverbose) print*,'fourier transforms done'
!
! Solve the five point matrix in fourier space
!
      a_band= dr1**2 - .5*dr1*rad1
      c_band= dr1**2 + .5*dr1*rad1
      d_band= dz1**2
      e_band= dz1**2
!
      if (lverbose) &
           print*,'initializing the iterations '//&
           'to solve the Poisson equation in the grid'
!
      do im=1,nkt
!
        b_band=-2*(dr1**2 + dz1**2) - (m_fft(im+ipy*nkt)*rad1)**2
        b_band1=1/b_band
!          
! check the points that need updating, then call the five point solver, in both
! real and imaginary parts
!
        call five_point_solver(phi(:,im,:),rhs(:,im,:), &
             a_band,b_band,b_band1,c_band,d_band,e_band,&
             lupdate_grid(:,im,:),lverbose)
!
        call five_point_solver(b1(:,im,:),b1_rhs(:,im,:),&
             a_band,b_band,b_band1,c_band,d_band,e_band, &
             lupdate_grid(:,im,:),lverbose)
!
      enddo
!
      if (lverbose) print*,'done with the iterative solution in the grid'
!
! Fourier transform back to real space
!
      call fourier_transform_y(phi,b1,linv=.true.)
!
      if (lfirstcall) lfirstcall=.false.
!
! Save the phi from the previous step
!
      phi_previous_step=phi
!
    endsubroutine inverse_laplacian_sor
!*************************************************************************
    subroutine five_point_solver(lhs,rhs,a_band,b_band,b_band1,&
         c_band,d_band,e_band,lupdate,lverbose)
!
! Invert a five point matrix 
!
!  u(i,j) = 1/b*(a*u(i-1,j) + c*u(i+1,j) + d*u(i,j+1) + 
!                e*u(i,j-1) + b*u(i,j)   - rhs(i,j))
!
! using Chebychev (checkboard, red and black) acceleration 
!
      real, dimension (nx,nz) :: lhs,rhs,lhs_old
      real, dimension (nx) :: a_band,b_band,c_band,d_band,e_band,b_band1
      real :: omega,norm,sig,norm_old,anorm,resid
      integer :: n,i,iteration,skipped
      logical, dimension(nx,nz) :: lupdate
      logical :: lverbose
!
! Starters
!
      sig=1e5 ; norm=1e5
      lhs_old=lhs
      iteration=0
!
      if (lverbose) &
           print*,'five_point_solver: Jacobian of the spectral matrix=',rjac
!
      do while (sig > iteration_threshold)
        iteration=iteration+1
        if (mod(iteration,2)/=0) skipped=0
        do n=2,nz-1
          do i=2,nr-1
            !chebychev : odd-even ordering
            if (mod(n+i,2) /= (mod(iteration,2))) then
              if (lupdate(i,n)) then
                resid=  a_band(i)*lhs(i-1,n)+      &
                        c_band(i)*lhs(i+1,n)+      &
                        d_band(i)*lhs(i,n+1)+      &
                        e_band(i)*lhs(i,n-1)+      & 
                        b_band(i)*lhs(i,n) - rhs(i,n)
                anorm=anorm+abs(resid)
                lhs(i,n)=lhs(i,n)-omega*resid*b_band1(i)
              else
                skipped=skipped+1
              endif
            endif
          enddo
        enddo
        
        if (iteration==1) then 
          omega=1./(1-.5 *rjac**2)
        else
          omega=1./(1-.25*rjac**2*omega)
        endif
!
! error of the iteration
!
        norm_old=norm
        norm=sum((lhs-lhs_old)**2)
        sig=abs(norm-norm_old)/max(norm,epsi)
!
        if (iteration > 1000) then 
          print*,'maximum number of iterations exceeded'
          call fatal_error("five_point_solver","")
        endif
      enddo
!
      if (lverbose)  then
        print*,'fps: skipped ',skipped,' of ',(nr-2)*(nz-2)
        print*,'number of iterations',iteration
      endif
!
    endsubroutine five_point_solver
!***********************************************************************
    subroutine get_border_values(rhs,phi,lupdate_grid,lverbose)
!
! Calculate the value of the potential in the borders of the grid  
!
! 28-apr-08/wlad: coded
!
      use Mpicomm
!
      real, dimension (nx,ny,nz) :: rhs,phi,rhs_proc
      real, dimension (nx,nygrid,nzgrid)  :: rhs_serial
      integer :: j,jy,jz,iydo,iyup,izdo,izup
      integer :: ir,ith,iz,skipped
      real    :: potential
      logical, dimension(nx,ny,nz) :: lupdate_grid 
      logical :: lcall_from_self
      logical :: lverbose 
!
! Construct the serial density
!
      if (lmpicomm) then
        if (.not.lroot) then
          call mpisend_real(rhs,(/nx,ny,nz/),root,111)
        else
          do jy=0,nprocy-1
            do jz=0,nprocz-1
              j=jy+nprocy*jz
              if (j/=0) then 
                call mpirecv_real(rhs_proc,(/nx,ny,nz/),j,111)
              else
                rhs_proc=rhs
              endif
              iydo= jy  * ny + 1 ; izdo= jz  * nz + 1 
              iyup=(jy+1)*ny     ; izup=(jz+1)*nz     
              rhs_serial(:,iydo:iyup,izdo:izup)=rhs_proc
            enddo
          enddo
        endif
        call mpibcast_real(rhs_serial,(/nx,nygrid,nzgrid/))
      else
        rhs_serial(:,1:ny,1:nz)=rhs(:,1:ny,1:nz)
      endif
!
! Integrate the borders. Count the points skipped because 
! they do not need update, for monitoring purposes. 
!
      skipped=0
      if (lverbose) & 
           print*,'get_border_values: integrating the border values only'
! Recompute the border. Start with ir=1
      do iz=1,nz;do ith=1,nth
        if (lupdate_grid(1,ith,iz)) then 
          call integrate_border(rhs_serial,1,ith,iz,potential)
          phi(1,ith,iz)=potential
        else
          skipped=skipped+1
          phi(1,ith,iz)=phi_previous_step(1,ith,iz)
        endif
      enddo;enddo
      if (lverbose) print*,'done for ir=1'
! ir=nr
      do iz=1,nz;do ith=1,nth
        if (lupdate_grid(nr,ith,iz)) then 
          call integrate_border(rhs_serial,nr,ith,iz,potential)
          phi(nr,ith,iz)=potential
        else
          skipped=skipped+1
          phi(nr,ith,iz)=phi_previous_step(nr,ith,iz)
        endif
      enddo;enddo
      if (lverbose) print*,'done for ir=nr'
!
! Vertical. Take into account parallelization. Just do it
! for the first and last z-processor. Start with z=1 for the
! first processor.
!
      if (lfirst_proc_z) then 
        do ir=2,nr-1;do ith=1,nth
          if (lupdate_grid(ir,ith,1)) then 
            call integrate_border(rhs_serial,ir,ith,1,potential)
            phi(ir,ith,1)=potential
          else
            skipped=skipped+1
            phi(ir,ith,1)=phi_previous_step(ir,ith,1)
          endif
        enddo;enddo
      endif
      if (lverbose) print*,'done for iz=1'
! z=nz for the last processor
      if (llast_proc_z) then 
        do ir=2,nr-1;do ith=1,nth
          if (lupdate_grid(ir,ith,nz)) then 
            call integrate_border(rhs_serial,ir,ith,nz,potential)
            phi(ir,ith,nz)=potential
          else
            skipped=skipped+1
            phi(ir,ith,nz)=phi_previous_step(ir,ith,nz)
          endif
        enddo;enddo
      endif
      if (lverbose) print*,'done for iz=nz'
!
      if (lverbose) &
           print*,'border: skipped ',skipped,' of ',2*nth*(nz+nr-2)
!
    endsubroutine get_border_values
!***********************************************************************
    subroutine integrate_border(rhs_serial,ir,ith,iz,potential)
!
! Get the potential in the borders by direct integration (see routine
! calculate_cross_legendre_functions below) 
!
! Phi(r,tht,z) = -G/(pi*sqrt(r)) * &
!                 Int r'dr'dz'dth'/sqrt(r') * rho' * &
!                 Sum_m=0^infty eps(m)*Q_(m-1/2)*cos(m*(tht-tht'))
! 
!                        _               infty
!                G      /  d3x' rho(x')  ___
! Phi(x) = _  _______   |  ___________   \   e(m) Q (chi) cos[m(t-t')]
!                  __   |      __        /__       m-1/2
!             pi \/r   _/    \/r'         m=0
!
!
! This is time consuming and only done for the border, if they need updating. 
! In the grid, we iteratively solve the discrete Poisson equation via 
! Chebychev acceleration. 
!
! 28-apr-08/wlad: coded
!
      real, dimension(nx,nygrid,nzgrid) :: rhs_serial
      real, dimension (nx) :: intr
      real, dimension (nygrid) :: intp
      real, dimension (nzgrid) :: intz
!
      real :: summation_over_harmonics,potential,fac
      integer :: ir,ith,iz,ikr,ikt,ikz,im
!
! as rhs already the 4*pi*G factor built in, the factor -G/pi in front of
! the integral becomes 1/(4*pi^2)
!
      fac=-.25*pi_1**2*sqrtrad_1(ir)
!
      do ikz=1,nzgrid
        do ikt=1,nthgrid
          do ikr=1,nr
!
            summation_over_harmonics=0
! The Legendre function here is modified by including
! the jacobian, the division by sqrt(rad) and the neumann
! epsilon factor
            do im=0,mmax
                summation_over_harmonics=summation_over_harmonics+&
                     Legendre_Qmod(ir,iz,ikr,ikz,im)*&
                     fourier_cosine_terms(ith,ikt,im)
            enddo
            intr(ikr)=rhs_serial(ikr,ikt,ikz)*summation_over_harmonics
          enddo
          intp(ikt)=sum(intr(2:nr-1))+.5*(intr(1)+intr(nr))
        enddo
        intz(ikz)=sum(intp)  !this is phi-periodic
      enddo
!
      potential=fac*(sum(intz(2:nzgrid-1))+.5*(intz(1)+intz(nzgrid)))
!  
    endsubroutine integrate_border
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
      integer :: ir,ikr,ikt,nw,count
!
      if (lroot) &
        print*,'Pre-calculating the Bessel functions to '//&
        'solve the Poisson equation'
!
      nw=10*int(nr*nkr*nkt/10.) ; count=0
!
      do ikt=1,nkt
        do ir=1,nr;do ikr=1,nkr
          arg=kr_fft(ikr)*rad(ir)
          call besselj_nu_int(tmp,m_fft(ikt+ipy*nkt),arg,&
               loversample=.false.)
          bessel_grid(ir,ikr,ikt)=tmp
!
! Print progress
!
          if (lroot) then
            count=count+1
            if (mod(10*count,nw)==0) print*, 100*count/nw,'% done'
          endif
!
        enddo;enddo
      enddo
!
      if (lroot) &
           print*,'Calculated all the needed Bessel functions'
!
    endsubroutine calculate_cross_bessel_functions
!***********************************************************************
    subroutine calculate_cross_legendre_functions
! 
! Calculate the Legendre functions of second kind related to
! the cylindrical grid. They are needed to solve the potential
! by the analytical expression (Cohl & Tohline 1999)
!
! Phi(r,tht,z) = -G/(pi*sqrt(r)) * &
!                 Int r'dr'dz'dth'/sqrt(r') * rho' * &
!                 Sum_m=0^infty eps(m)*Q_(m-1/2)*cos(m*(tht-tht'))
! 
!                        _               infty
!                G      /  d3x' rho(x')  ___
! Phi(x) = _  _______   |  ___________   \   e(m) Q (chi) cos[m(t-t')]
!                  __   |      __        /__       m-1/2
!             pi \/r   _/    \/r'         m=0
!
! The Neumann factor e(m) is 1 for m=0 and 2 for other harmonics, 
! chi is a function of the grid
!          
!     chi= [r^2 + r'^2 + (z-z')^2]/(2rr')
!
! And the Legendre functions can be found by a recurrency relation
!
!  Q_(-1/2) = mu*K(mu)
!  Q_( 1/2) = chi*mu*K(mu) - (1+chi)*mu*E(mu)
!
!                4*(m-1)                   2m-3
!  Q_(m-1/2) =   _______  chi*Q_(m-3/2) - ______  Q_(m-5/2) 
!
!                 2m-1                     2m-1
!
! where mu=sqrt(2./(1+chi)) and K and E are the first and second
! kind complete elliptical integrals. The summation is truncated 
! at a maximum harmonic mmax
!
! 28-apr-08/wlad: coded
!
      use General, only: calc_complete_ellints
!
      real, dimension(0:mmax)   :: Legendre_Q
      real, dimension(nzgrid) :: zed_serial
      real, dimension(nygrid) :: tht_serial
      integer, dimension(mmax) :: neumann_factor_eps
      real    :: chi,mu,Kappa_mu,E_mu,jac
      integer :: ir,ith,iz,ikr,ikt,ikz
      integer :: j,im,ith_serial,nw,count
!
      if (lroot) &
        print*,'Pre-calculating the half-integer Legendre functions '//&
        'of second kind to solve the Poisson equation'     
!
! Get serial z's
!
      if (nprocz /= 1) then
        call get_serial_array(zed,zed_serial,'z')
      else
        zed_serial(1:nz)=zed(1:nz)
      endif
!
      nw=10*int(nr**2*nz*nzgrid/10) ; count=0
!
      do ir=1,nr;do iz=1,nz        !for all points in this processor 
      do ikr=1,nr;do ikz=1,nzgrid  !integrate over the whole r-z grid
!
! Jacobian
!
        jac=rad(ikr)*dr*dth*dz 
!
! Calculate the elliptic integrals
! 
        chi=(rad(ir)**2+rad(ikr)**2+(zed(iz)-zed_serial(ikz))**2)/&
             (2*rad(ir)*rad(ikr))
        mu=sqrt(2./(1+chi))
        call calc_complete_ellints(mu,Kappa_mu,E_mu,&
             loversample=.false.)
!
! Calculate the Legendre functions for each harmonic
!
        Legendre_Q(0)=    mu*Kappa_mu
        Legendre_Q(1)=chi*mu*Kappa_mu - (1+chi)*mu*E_mu
!
        do im=0,mmax
          if (im == 0) then
            neumann_factor_eps(im)=1
          else
            neumann_factor_eps(im)=2
            if (im >= 2) then 
              Legendre_Q(im)= &
                   4*(im-1)/(2*im-1)*chi*Legendre_Q(im-1) - &
                   (2*im-3)/(2*im-1)    *Legendre_Q(im-2)
            endif
          endif
!
! modify the Legendre function by multiplying the jacobian
! to speed up runtime. Also stick in the neumann factor and 
! the 1/sqrt(rad) factor
!
          Legendre_Qmod(ir,iz,ikr,ikz,im)=&
               Legendre_Q(im)*neumann_factor_eps(im)*jac/sqrt(rad(ikr))
!
        enddo
! 
! Finished grid integration of the Legendre functions
!
        if (lroot) then
          count=count+1
          if (mod(10*count,nw)==0) print*, 100*count/nw,'% done'
        endif
!
      enddo;enddo
      enddo;enddo
!
! Now the co-sines in the azimuthal direction
!
      if (nprocy /= 1) then
        call get_serial_array(tht,tht_serial,'y')
      else
        tht_serial(1:ny)=tht(1:ny)
      endif
!
      if (lroot) print*,'Start calculating the co-sines'
!
      do ith=1,nth;do ikt=1,nthgrid  
        do im=0,mmax
          fourier_cosine_terms(ith,ikt,im)=&
               cos(im*(tht(ith)-tht_serial(ikt)))
        enddo
      enddo;enddo
!
      if (lroot) &
           print*,'Calculated all the needed Legendre functions'
!
    endsubroutine calculate_cross_legendre_functions
!***********************************************************************
    subroutine get_serial_array(array,array_serial,var)
!
! Feed in the array present on the processor, and return 
! (broadcast) a serial array with the elements of that 
! array in all processors. Used for getting all y or z 
! elements. 
!
! 28-apr-08/wlad: coded
!
      use Mpicomm
!
      real, dimension(:) :: array
      real, dimension(:) :: array_serial
      real, dimension(:), allocatable :: array_proc
      
      integer :: ido,iup,j,nk,nkgrid,nprock,jy,jz
      character :: var
!
      if (var=='y') then 
        nk=ny ; nkgrid=nygrid ; nprock=nprocy
      else if (var=='z') then
        nk=nz ; nkgrid=nzgrid ; nprock=nprocz
      else
        print*,'var=',var
        call stop_it("you can only call it for var='y' or 'z'")
      endif
!
      allocate(array_proc(nk))
!
      if (lmpicomm) then
!
! All processors send its array to the root processor
!
        if (.not.lroot) then
          call mpisend_real(array,nk,root,111)
        else
!
          array_proc=array
!
! The root processor receives all arrays and
! stores them in a single array of dimension nkgrid
!
          do jy=0,nprocy-1
            do jz=0,nprocz-1
              !serial index of the processor
              j=jy+nprocy*jz
              if (j/=0) call mpirecv_real(array_proc,nk,j,111)
!
              if (var=='y') then 
                ido= jy  * nk + 1
                iup=(jy+1)*nk
                array_serial(ido:iup)=array_proc
              endif
!
              if (var=='z') then 
                ido= jz  * nk + 1
                iup=(jz+1)*nk
                array_serial(ido:iup)=array_proc
              endif
!
            enddo
          enddo
        endif
!
! Broadcast the serial array to all processors
!
        call mpibcast_real(array_serial,nkgrid)
!
      else
!
! For serial runs, nk=nkgrid, so just copy the array
!
        array_serial(1:nk)=array(1:nk)
!
      endif
!
    endsubroutine get_serial_array
!***********************************************************************
    subroutine read_poisson_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=poisson_init_pars, IOSTAT=iostat)
!
    endsubroutine read_poisson_init_pars
!***********************************************************************
    subroutine write_poisson_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=poisson_init_pars)
!
    endsubroutine write_poisson_init_pars
!***********************************************************************
    subroutine read_poisson_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=poisson_run_pars, IOSTAT=iostat)
!
    endsubroutine read_poisson_run_pars
!***********************************************************************
    subroutine write_poisson_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=poisson_run_pars)
!
    endsubroutine write_poisson_run_pars
!***********************************************************************
    subroutine get_acceleration(acceleration)
!
      use General, only: keep_compiler_quiet
!
      real, dimension(nx,ny,nz,3), intent(out) :: acceleration           !should I (CAN I?) make this allocatable?
!
      call keep_compiler_quiet(acceleration)
!
    endsubroutine get_acceleration
!***********************************************************************
endmodule Poisson
